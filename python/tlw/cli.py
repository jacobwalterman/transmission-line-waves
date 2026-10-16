import argparse
import tomllib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3d projection)

import _transmission_line_waves as engine


# ---------------------------------------------------------------- config

@dataclass
class Config:
    params: engine.Params
    l1: float
    l2: float
    l3: float
    l4: float
    times: np.ndarray


def make_times(spec):
    start, stop, num = spec["start"], spec["stop"], spec["num"]
    spacing = spec.get("spacing", "linear")
    if spacing == "linear":
        return np.linspace(start, stop, num)
    if spacing == "log":
        return np.geomspace(start, stop, num)
    raise ValueError(f"unknown time spacing {spacing!r}")


def load_config(path, overrides=None):
    with open(path, "rb") as f:
        v = tomllib.load(f)
    if overrides:
        v.update(overrides)

    times = make_times(v.pop("times"))
    l = [v.pop(f"l{i}") for i in range(1, 5)]
    h = [l[0] / (v["n1"] + 1)]
    h += [(l[i] - l[i - 1]) / (v[f"n{i + 1}"] + 1) for i in range(1, 4)]
    for i in range(4):
        v[f"h{i + 1}"] = h[i]

    return Config(params=engine.Params(**v), l1=l[0], l2=l[1], l3=l[2], l4=l[3], times=times)


# ------------------------------------------------------- grids and state

def build_grids(params):
    """Returns ((w1..w4 grids), (v1..v4 grids)).

    w1 uses i=1..n1+1; w2/w3/w4 start at i=0. v-grids all start at i=0.
    """
    p = params
    start2 = (p.n1 + 1) * p.h1
    start3 = start2 + (p.n2 + 1) * p.h2
    start4 = start3 + (p.n3 + 1) * p.h3

    def grid(start, h, first, stop):
        return start + h * np.arange(first, stop)

    w_grids = (
        grid(0.0, p.h1, 1, p.n1 + 2),
        grid(start2, p.h2, 0, p.n2 + 1),
        grid(start3, p.h3, 0, p.n3 + 2),
        grid(start4, p.h4, 0, p.n4 + 1),
    )
    v_grids = (
        grid(0.0, p.h1, 0, p.n1 + 2),
        grid(start2, p.h2, 0, p.n2 + 1),
        grid(start3, p.h3, 0, p.n3 + 2),
        grid(start4, p.h4, 0, p.n4 + 2),
    )
    return w_grids, v_grids


def split_state(state, params):
    """Splits the raw [w1..w4,v1..v4] state and patches in the boundary
    values that demonstrate[] defines outside the state vector proper:
    w1[0]=0, w2[n2+1]=w3[0], w4[n4+1]=0, v2[n2+1]=(c3^2/c2^2)*v3[0]."""
    p = params
    lengths = [p.n1 + 1, p.n2 + 1, p.n3 + 2, p.n4 + 1,
               p.n1 + 2, p.n2 + 1, p.n3 + 2, p.n4 + 2]
    w1, w2, w3, w4, v1, v2, v3, v4 = np.split(state, np.cumsum(lengths)[:-1])

    w1 = np.concatenate([[0.0], w1])
    w2 = np.concatenate([w2, [w3[0]]])
    w4 = np.concatenate([w4, [0.0]])
    v2 = np.concatenate([v2, [(p.c3 ** 2 / p.c2 ** 2) * v3[0]]])

    return w1, w2, w3, w4, v1, v2, v3, v4


# ----------------------------------------------------- initial conditions

def create_box(grid, l_first, l_last):
    box_length = l_last - l_first
    left = l_first + box_length / 4
    right = l_last - box_length / 4
    return np.where((grid >= left) & (grid <= right), 1.0, 0.0)


def create_sin(grid, l4):
    return np.sin(2 * np.pi * grid / l4)


def create_initial_conditions(init_cond, config):
    p = config.params
    w_grids, v_grids = build_grids(p)
    bounds = [(0, config.l1), (config.l1, config.l2), (config.l2, config.l3), (config.l3, config.l4)]
    signs = [-1.0, 1.0, -1.0, 1.0]

    if init_cond == "box":
        w = [s * create_box(g, *b) for s, g, b in zip(signs, w_grids, bounds)]
        v = [s * create_box(g, *b) for s, g, b in zip(signs, v_grids, bounds)]
    elif init_cond == "flat":
        flat_signs = [-1.0, 1.0, 1.0, -1.0]
        w = [np.full(len(g), s) for s, g in zip(flat_signs, w_grids)]
        v = [np.full(len(g), s) for s, g in zip(flat_signs, v_grids)]
    elif init_cond == "sin":
        w = [create_sin(g, config.l4) for g in w_grids]
        v = [create_sin(g, config.l4) for g in v_grids]
    else:
        raise ValueError(f"unknown init_cond {init_cond!r}")

    return np.concatenate(w + v)


# --------------------------------------------------------------- physics

def energy_of_section(h, c, w, v):
    wj, wj1 = w[:-1], w[1:]
    vj, vj1 = v[:-1], v[1:]
    return (h / 2) * np.sum(((wj1 + wj) / 2) ** 2 + c * ((vj1 + vj) / 2) ** 2)


def total_energy(params, w1, w2, w3, w4, v1, v2, v3, v4):
    p = params
    return (
        energy_of_section(p.h1, p.c1, w1, v1)
        + energy_of_section(p.h2, p.c2, w2, v2)
        + energy_of_section(p.h3, p.c3, w3, v3)
        + energy_of_section(p.h4, p.c4, w4, v4)
    )


def solve(config, init_cond, times):
    """Returns the raw state vector at each time in `times`."""
    initial_conditions = create_initial_conditions(init_cond, config)
    es = engine.decompose(config.params)
    coefficients = engine.coefficients(es, initial_conditions)
    # evolve() drops modes, so t=0 is returned exactly rather than reconstructed.
    return [initial_conditions if t == 0 else engine.evolve(es, coefficients, t) for t in times]


# -------------------------------------------------------------- plotting

def plot_surface_3d(ax, x, times, vals, label, cmap):
    """Matches generateALittlePlot: a ListPlot3D-style surface with mesh
    lines at every grid point and every recorded time step."""
    X, T = np.meshgrid(x, times)
    ax.plot_surface(
        X, T, vals, cmap=cmap, edgecolor="k", linewidth=0.15,
        antialiased=True, rstride=1, cstride=1,
    )
    ax.set_xlabel("x (m)", fontsize=10)
    ax.set_ylabel("t (sec)", fontsize=10)
    ax.set_title(label, fontsize=13, fontweight="bold")
    ax.view_init(elev=25, azim=-60)


def plot_eigenvalues(eigenvalues, drop=4, out_path="eigenvalues.png"):
    """Scatter eigenvalues in the complex plane.

    eigenvalues: complex 1D array, sorted by descending magnitude
    drop: number of trailing modes that evolve() discards
    """
    ev = np.asarray(eigenvalues)
    kept = ev[:len(ev) - drop] if drop > 0 else ev
    dropped = ev[len(ev) - drop:] if drop > 0 else ev[:0]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(kept.real, kept.imag, s=12, label="kept")
    if len(dropped):
        ax.scatter(dropped.real, dropped.imag, s=30, marker="x",
                   color="red", label=f"dropped ({drop})")

    ax.axvline(0, color="gray", linewidth=0.8)
    ax.set_xlabel("Re(λ)")
    ax.set_ylabel("Im(λ)")
    ax.set_title("Eigenvalues in the complex plane")
    ax.legend()
    ax.set_xlim(left=kept.real.min() * 0.05, right=0.05 * abs(kept.real.min()))
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# -------------------------------------------------------------- commands

def run_demo(config, init_cond="box", out_path="transmission_line.png"):
    params, times = config.params, config.times
    (g1w, g2w, g3w, g4w), (g1v, g2v, g3v, g4v) = build_grids(params)
    # x-coordinates matching the patched w/v arrays (w1 gets a leading
    # 0-node, w2/v2 get a trailing node shared with body 3's start).
    x_w = np.concatenate([[0.0], g1w, g2w, [g3w[0]], g3w, g4w, [config.l4]])
    x_v = np.concatenate([g1v, g2v, [g3v[0]], g3v, g4v])

    w_profiles, v_profiles, energies = [], [], []
    for state in solve(config, init_cond, times):
        sections = split_state(state, params)
        w_profiles.append(np.concatenate(sections[:4]))
        v_profiles.append(np.concatenate(sections[4:]))
        energies.append(total_energy(params, *sections))

    normalized_energy = np.array(energies) / energies[0]

    fig = plt.figure(figsize=(19, 6))
    ax1 = fig.add_subplot(1, 3, 1, projection="3d")
    ax2 = fig.add_subplot(1, 3, 2, projection="3d")
    ax3 = fig.add_subplot(1, 3, 3)

    plot_surface_3d(ax1, x_w, times, np.array(w_profiles), "w1 through w4", "coolwarm")
    plot_surface_3d(ax2, x_v, times, np.array(v_profiles), "v1 through v4", "coolwarm")

    ax3.plot(times, normalized_energy)
    ax3.set_xlabel("t (sec)")
    ax3.set_ylabel("E(t) / E(0)")
    ax3.set_title("Normalized energy")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"wrote {out_path}")


def run_spectrum(config, out_path="eigenvalues.png"):
    es = engine.decompose(config.params)
    plot_eigenvalues(es.sortedEigenvalues(), out_path=out_path)
    print(f"wrote {out_path}")


DEFAULT_CONFIG = Path(__file__).parent / "params.toml"


def main():
    parser = argparse.ArgumentParser(prog="tlw", description="Transmission-line wave solver")
    parser.add_argument("init", nargs="?", default="box", choices=["box", "flat", "sin"],
                        help="initial condition (default: box)")
    parser.add_argument("--out", default="transmission_line.png", help="output PNG path")
    parser.add_argument("--config", default=None, help="params TOML file (default: built-in)")
    args = parser.parse_args()
    run_demo(load_config(args.config or DEFAULT_CONFIG), args.init, args.out)


if __name__ == "__main__":
    main()
