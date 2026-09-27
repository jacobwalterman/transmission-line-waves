import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3d projection)

BUILD_PATH = Path(__file__).resolve().parent.parent / "build"
sys.path.insert(0, str(BUILD_PATH))

import _transmission_line_waves


def create_box(grids, l_first, l_last):
    box_length = l_last - l_first
    left = l_first + box_length / 4
    right = l_last - box_length / 4
    return np.where((grids >= left) & (grids <= right), 1.0, 0.0)


def create_sin(grids, l4):
    return np.sin(2 * np.pi * grids / l4)


def build_grids(n1, n2, n3, n4, h1, h2, h3, h4):
    # w-grids: w1 uses i=1..n1+1, w2/w3/w4 use i=0..(n+1 or n)
    first_grids_w = np.array([i * h1 for i in range(1, n1 + 2)])
    second_grids_w = np.array([(n1 + 1) * h1 + i * h2 for i in range(0, n2 + 1)])
    third_grids_w = np.array(
        [(n1 + 1) * h1 + (n2 + 1) * h2 + i * h3 for i in range(0, n3 + 2)]
    )
    fourth_grids_w = np.array(
        [(n1 + 1) * h1 + (n2 + 1) * h2 + (n3 + 1) * h3 + i * h4
         for i in range(0, n4 + 1)]
    )

    # v-grids: all i=0..n+1
    first_grids_v = np.array([i * h1 for i in range(0, n1 + 2)])
    second_grids_v = np.array([(n1 + 1) * h1 + i * h2 for i in range(0, n2 + 1)])
    third_grids_v = np.array(
        [(n1 + 1) * h1 + (n2 + 1) * h2 + i * h3 for i in range(0, n3 + 2)]
    )
    fourth_grids_v = np.array(
        [(n1 + 1) * h1 + (n2 + 1) * h2 + (n3 + 1) * h3 + i * h4
         for i in range(0, n4 + 2)]
    )

    return (
        (first_grids_w, second_grids_w, third_grids_w, fourth_grids_w),
        (first_grids_v, second_grids_v, third_grids_v, fourth_grids_v),
    )


def create_initial_conditions(init_cond, n1, n2, n3, n4, h1, h2, h3, h4,
                               l1, l2, l3, l4):
    (first_grids_w, second_grids_w, third_grids_w, fourth_grids_w), \
        (first_grids_v, second_grids_v, third_grids_v, fourth_grids_v) = \
        build_grids(n1, n2, n3, n4, h1, h2, h3, h4)

    if init_cond == "box":
        w1 = -create_box(first_grids_w, 0, l1)
        w2 = create_box(second_grids_w, l1, l2)
        w3 = -create_box(third_grids_w, l2, l3)
        w4 = create_box(fourth_grids_w, l3, l4)
        v1 = -create_box(first_grids_v, 0, l1)
        v2 = create_box(second_grids_v, l1, l2)
        v3 = -create_box(third_grids_v, l2, l3)
        v4 = create_box(fourth_grids_v, l3, l4)
    elif init_cond == "flat":
        w1 = np.full(n1 + 1, -1.0)
        w2 = np.full(n2 + 1, 1.0)
        w3 = np.full(n3 + 2, 1.0)
        w4 = np.full(n4 + 1, -1.0)
        v1 = np.full(n1 + 2, -1.0)
        v2 = np.full(n2 + 1, 1.0)
        v3 = np.full(n3 + 2, 1.0)
        v4 = np.full(n4 + 2, -1.0)
    elif init_cond == "sin":
        w1 = create_sin(first_grids_w, l4)
        w2 = create_sin(second_grids_w, l4)
        w3 = create_sin(third_grids_w, l4)
        w4 = create_sin(fourth_grids_w, l4)
        v1 = create_sin(first_grids_v, l4)
        v2 = create_sin(second_grids_v, l4)
        v3 = create_sin(third_grids_v, l4)
        v4 = create_sin(fourth_grids_v, l4)
    else:
        raise ValueError(f"unknown init_cond {init_cond!r}")

    return np.concatenate([w1, w2, w3, w4, v1, v2, v3, v4])


def split_state(state, n1, n2, n3, n4, c2, c3):
    """Splits the raw [w1..w4,v1..v4] state and patches in the boundary
    values that demonstrate[] defines outside the state vector proper:
    w1[0]=0, w2[n2+1]=w3[0], w4[n4+1]=0, v2[n2+1]=(c3^2/c2^2)*v3[0]."""
    w1_len, w2_len, w3_len, w4_len = n1 + 1, n2 + 1, n3 + 2, n4 + 1
    v1_len, v2_len, v3_len, v4_len = n1 + 2, n2 + 1, n3 + 2, n4 + 2

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
    v2 = np.concatenate([v2, [(c3 ** 2 / c2 ** 2) * v3[0]]])

    return w1, w2, w3, w4, v1, v2, v3, v4


def energy_of_section(h, c, w, v):
    wj, wj1 = w[:-1], w[1:]
    vj, vj1 = v[:-1], v[1:]
    return (h / 2) * np.sum(((wj1 + wj) / 2) ** 2 + c * ((vj1 + vj) / 2) ** 2)


def total_energy(h1, h2, h3, h4, c1, c2, c3, c4, w1, w2, w3, w4, v1, v2, v3, v4):
    return (
        energy_of_section(h1, c1, w1, v1)
        + energy_of_section(h2, c2, w2, v2)
        + energy_of_section(h3, c3, w3, v3)
        + energy_of_section(h4, c4, w4, v4)
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


def run(init_cond="box", out_path="transmission_line.png"):
    # Defaults from the Mathematica Manipulate[] front end.
    n1, n2, n3, n4 = 15, 20, 25, 15
    l1, l2, l3, l4 = 1.0, 3.0, 7.0, 7.5
    h1 = l1 / (n1 + 1)
    h2 = (l2 - l1) / (n2 + 1)
    h3 = (l3 - l2) / (n3 + 1)
    h4 = (l4 - l3) / (n4 + 1)
    c1, c2, c3, c4 = 7.0, 4.0, 3.0, 4.0
    k1, k2, k3, k4 = 1.0, 2.0, 3.0, 4.0

    (g1w, g2w, g3w, g4w), (g1v, g2v, g3v, g4v) = build_grids(
        n1, n2, n3, n4, h1, h2, h3, h4
    )
    initial_conditions = create_initial_conditions(
        init_cond, n1, n2, n3, n4, h1, h2, h3, h4, l1, l2, l3, l4
    )
    # x-coordinates matching the patched w/v arrays (w1 gets a leading
    # 0-node, w2/v2 get a trailing node shared with body 3's start).
    x_w = np.concatenate([[0.0], g1w, g2w, [g3w[0]], g3w, g4w, [l4]])
    x_v = np.concatenate([g1v, g2v, [g3v[0]], g3v, g4v])

    end_time = 4.0
    time_iterations = 40
    times = np.linspace(0, end_time, time_iterations + 1)

    w_profiles, v_profiles, energies = [], [], []
    for t in times:
        state = np.array(
            _transmission_line_waves.evaluate_transmission_line_wave(
                n1, n2, n3, n4, h1, h2, h3, h4, c1, c2, c3, c4, k1, k2, k3, k4,
                initial_conditions, t,
            )
        )
        w1, w2, w3, w4, v1, v2, v3, v4 = split_state(state, n1, n2, n3, n4, c2, c3)
        w_profiles.append(np.concatenate([w1, w2, w3, w4]))
        v_profiles.append(np.concatenate([v1, v2, v3, v4]))
        energies.append(
            total_energy(h1, h2, h3, h4, c1, c2, c3, c4, w1, w2, w3, w4, v1, v2, v3, v4)
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


if __name__ == "__main__":
    init_cond = sys.argv[1] if len(sys.argv) > 1 else "box"
    out_path = sys.argv[2] if len(sys.argv) > 2 else "transmission_line.png"
    run(init_cond, out_path)
