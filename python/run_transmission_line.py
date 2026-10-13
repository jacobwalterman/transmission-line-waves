import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3d projection)

import _transmission_line_waves


def create_box(grids, l_first, l_last):
    box_length = l_last - l_first
    left = l_first + box_length / 4
    right = l_last - box_length / 4
    return np.where((grids >= left) & (grids <= right), 1.0, 0.0)


def create_sin(grids, l4):
    return np.sin(2 * np.pi * grids / l4)


def build_grids(params):
    # w-grids: w1 uses i=1..n1+1, w2/w3/w4 use i=0..(n+1 or n)
    first_grids_w = np.array([i * params.h1 for i in range(1, params.n1 + 2)])
    second_grids_w = np.array([(params.n1 + 1) * params.h1 + i * params.h2 for i in range(0, params.n2 + 1)])
    third_grids_w = np.array(
        [(params.n1 + 1) * params.h1 + (params.n2 + 1) * params.h2 + i * params.h3 for i in range(0, params.n3 + 2)]
    )
    fourth_grids_w = np.array(
        [(params.n1 + 1) * params.h1 + (params.n2 + 1) * params.h2 + (params.n3 + 1) * params.h3 + i * params.h4
         for i in range(0, params.n4 + 1)]
    )

    # v-grids: all i=0..n+1
    first_grids_v = np.array([i * params.h1 for i in range(0, params.n1 + 2)])
    second_grids_v = np.array([(params.n1 + 1) * params.h1 + i * params.h2 for i in range(0, params.n2 + 1)])
    third_grids_v = np.array(
        [(params.n1 + 1) * params.h1 + (params.n2 + 1) * params.h2 + i * params.h3 for i in range(0, params.n3 + 2)]
    )
    fourth_grids_v = np.array(
        [(params.n1 + 1) * params.h1 + (params.n2 + 1) * params.h2 + (params.n3 + 1) * params.h3 + i * params.h4
         for i in range(0, params.n4 + 2)]
    )
    return (
        (first_grids_w, second_grids_w, third_grids_w, fourth_grids_w),
        (first_grids_v, second_grids_v, third_grids_v, fourth_grids_v),
    )


def create_initial_conditions(init_cond, config):
    (first_grids_w, second_grids_w, third_grids_w, fourth_grids_w), \
        (first_grids_v, second_grids_v, third_grids_v, fourth_grids_v) = \
        build_grids(config.params)

    if init_cond == "box":
        w1 = -create_box(first_grids_w, 0, config.l1)
        w2 = create_box(second_grids_w, config.l1, config.l2)
        w3 = -create_box(third_grids_w, config.l2, config.l3)
        w4 = create_box(fourth_grids_w, config.l3, config.l4)
        v1 = -create_box(first_grids_v, 0, config.l1)
        v2 = create_box(second_grids_v, config.l1, config.l2)
        v3 = -create_box(third_grids_v, config.l2, config.l3)
        v4 = create_box(fourth_grids_v, config.l3, config.l4)
    elif init_cond == "flat":
        w1 = np.full(config.params.n1 + 1, -1.0)
        w2 = np.full(config.params.n2 + 1, 1.0)
        w3 = np.full(config.params.n3 + 2, 1.0)
        w4 = np.full(config.params.n4 + 1, -1.0)
        v1 = np.full(config.params.n1 + 2, -1.0)
        v2 = np.full(config.params.n2 + 1, 1.0)
        v3 = np.full(config.params.n3 + 2, 1.0)
        v4 = np.full(config.params.n4 + 2, -1.0)
    elif init_cond == "sin":
        w1 = create_sin(first_grids_w, config.l4)
        w2 = create_sin(second_grids_w, config.l4)
        w3 = create_sin(third_grids_w, config.l4)
        w4 = create_sin(fourth_grids_w, config.l4)
        v1 = create_sin(first_grids_v, config.l4)
        v2 = create_sin(second_grids_v, config.l4)
        v3 = create_sin(third_grids_v, config.l4)
        v4 = create_sin(fourth_grids_v, config.l4)
    else:
        raise ValueError(f"unknown init_cond {init_cond!r}")

    return np.concatenate([w1, w2, w3, w4, v1, v2, v3, v4])


def split_state(state, params):
    """Splits the raw [w1..w4,v1..v4] state and patches in the boundary
    values that demonstrate[] defines outside the state vector proper:
    w1[0]=0, w2[n2+1]=w3[0], w4[n4+1]=0, v2[n2+1]=(c3^2/c2^2)*v3[0]."""
    w1_len, w2_len, w3_len, w4_len = params.n1 + 1, params.n2 + 1, params.n3 + 2, params.n4 + 1
    v1_len, v2_len, v3_len, v4_len = params.n1 + 2, params.n2 + 1, params.n3 + 2, params.n4 + 2

    i = 0
    w1 = state[i:i + w1_len]; i += w1_len
    w2 = state[i:i + w2_len]; i += w2_len
    w3 = state[i:i + w3_len]; i += w3_len
    w4 = state[i:i + w4_len]; i += w4_len
    v1 = state[i:i + v1_len]; i += v1_len
    v2 = state[i:i + v2_len]; i += v2_len
    v3 = state[i:i + v3_len]; i += v3_len
    v4 = state[i:i + v4_len]; i += v4_len

    w1 = np.concatenate([[0.0], w1])
    w2 = np.concatenate([w2, [w3[0]]])
    w4 = np.concatenate([w4, [0.0]])
    v2 = np.concatenate([v2, [(params.c3 ** 2 / params.c2 ** 2) * v3[0]]])

    return w1, w2, w3, w4, v1, v2, v3, v4


def energy_of_section(h, c, w, v):
    wj, wj1 = w[:-1], w[1:]
    vj, vj1 = v[:-1], v[1:]
    return (h / 2) * np.sum(((wj1 + wj) / 2) ** 2 + c * ((vj1 + vj) / 2) ** 2)


def total_energy(params, w1, w2, w3, w4, v1, v2, v3, v4):
    return (
        energy_of_section(params.h1, params.c1, w1, v1)
        + energy_of_section(params.h2, params.c2, w2, v2)
        + energy_of_section(params.h3, params.c3, w3, v3)
        + energy_of_section(params.h4, params.c4, w4, v4)
    )


def little_plot_3d(ax, x, times, vals, label, cmap):
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

def run_demo(config, init_cond="box", out_path="transmission_line.png"):
    # Defaults from the Mathematica Manipulate[] front end.

    (g1w, g2w, g3w, g4w), (g1v, g2v, g3v, g4v) = build_grids(config.params)
    initial_conditions = create_initial_conditions(init_cond, config)
    # x-coordinates matching the patched w/v arrays (w1 gets a leading
    # 0-node, w2/v2 get a trailing node shared with body 3's start).
    x_w = np.concatenate([[0.0], g1w, g2w, [g3w[0]], g3w, g4w, [config.l4]])
    x_v = np.concatenate([g1v, g2v, [g3v[0]], g3v, g4v])

    end_time = 4.0
    time_iterations = 40
    times = np.linspace(0, end_time, time_iterations + 1)

    w_profiles, v_profiles, energies = [], [], []
    for t in times:
        state = run2(t)
        w1, w2, w3, w4, v1, v2, v3, v4 = split_state(state, config.params)
        w_profiles.append(np.concatenate([w1, w2, w3, w4]))
        v_profiles.append(np.concatenate([v1, v2, v3, v4]))
        energies.append(
            total_energy(params, w1, w2, w3, w4, v1, v2, v3, v4)
        )

    w_profiles = np.array(w_profiles)
    v_profiles = np.array(v_profiles)
    energies = np.array(energies)
    normalized_energy = energies / energies[0]

    fig = plt.figure(figsize=(19, 6))
    ax1 = fig.add_subplot(1, 3, 1, projection="3d")
    ax2 = fig.add_subplot(1, 3, 2, projection="3d")
    ax3 = fig.add_subplot(1, 3, 3)

    little_plot_3d(ax1, x_w, times, w_profiles, "w1 through w4", "coolwarm")
    little_plot_3d(ax2, x_v, times, v_profiles, "v1 through v4", "coolwarm")

    ax3.plot(times, normalized_energy)
    ax3.set_xlabel("t (sec)"); ax3.set_ylabel("E(t) / E(0)")
    ax3.set_title("Normalized energy")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(f"wrote {out_path}")


def run2(t):
    n1, n2, n3, n4 = 15, 20, 25, 15
    l1, l2, l3, l4 = 1.0, 3.0, 7.0, 7.5
    h1 = l1 / (n1 + 1)
    h2 = (l2 - l1) / (n2 + 1)
    h3 = (l3 - l2) / (n3 + 1)
    h4 = (l4 - l3) / (n4 + 1)
    c1, c2, c3, c4 = 7.0, 4.0, 3.0, 4.0
    k1, k2, k3, k4 = 1.0, 2.0, 3.0, 4.0
    
    initial_conditions = create_initial_conditions(
        "sin", n1, n2, n3, n4, h1, h2, h3, h4, l1, l2, l3, l4
    )
    if t == 0:
        return initial_conditions
    params = _transmission_line_waves.Params(n1, n2, n3, n4, h1, h2, h3, h4, c1, c2, c3, c4, k1, k2, k3, k4)
    es = _transmission_line_waves.decompose(params)
    fourier_coefficients = _transmission_line_waves.coefficients(es, initial_conditions)
    state = _transmission_line_waves.evolve(es, fourier_coefficients, t)
    return state

import argparse

import tomllib
import _transmission_line_waves as engine

from dataclasses import dataclass

@dataclass
class Config:
    params: engine.Params
    l1: float
    l2: float
    l3: float
    l4: float

def load_config(path, overrides=None):
    with open(path, "rb") as f:
        v = tomllib.load(f)
    if overrides:
        v.update(overrides)

    l = [v.pop(f"l{i}") for i in range(1, 5)]
    h = [l[0] / (v["n1"] + 1)]
    h += [(l[i] - l[i-1]) / (v[f"n{i+1}"] + 1) for i in range(1, 4)]
    for i in range(4):
        v[f"h{i+1}"] = h[i]

    return Config(params=engine.Params(**v), l1=l[0], l2=l[1], l3=l[2], l4=l[3])

import numpy as np
import matplotlib.pyplot as plt

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
    fig.tight_layout()
    ax.set_xlim(left=kept.real.min() * 0.05, right=0.05 * abs(kept.real.min()))
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(prog="tlw")
    sub = parser.add_subparsers(dest="command", required=True)

    demo = sub.add_parser("demo", help="run the plotting demo")
    demo.add_argument("--init", default="box", help="initial condition name")
    demo.add_argument("--out", default="transmission_line.png")

    s= sub.add_parser("spectrum", help="calculate the eigenvalues")
    args = parser.parse_args()

    config = load_config("params.toml")
    # if args.command == "demo":
    #     run_demo(params, args.init, args.out)
    es = _transmission_line_waves.decompose(config.params)
    if args.command == "spectrum":
        eigenvalues = es.sortedEigenvalues()
        plot_eigenvalues(eigenvalues)
        return
    # init_cond = "sin"
    # initial_conditions = create_initial_conditions(init_cond, config)
    # t = 4
    # fourier_coefficients = _transmission_line_waves.coefficients(es, initial_conditions)
    # state = _transmission_line_waves.evolve(es, fourier_coefficients, t)
if __name__ == "__main__":
    main()
