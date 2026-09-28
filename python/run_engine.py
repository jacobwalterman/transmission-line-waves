import numpy as np
import matplotlib.pyplot as plt

import _transmission_line_waves

def run_engine(n, t):
    return _transmission_line_waves.evaluate_solution(n, t)

if __name__ == "__main__":
    n = 10
    t = 0;
    values = np.array(run_engine(n, t))

    values = values.real

    half = len(values) // 2
    positions = values[:half]
    velocities = values[half:]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))

    ax1.plot(positions, marker="o")
    ax1.set_title("Positions")
    ax1.set_xlabel("Grid index")
    ax1.set_ylabel("Position")

    ax2.plot(velocities, marker="o", color="tab:orange")
    ax2.set_title("Velocities")
    ax2.set_xlabel("Grid index")
    ax2.set_ylabel("Velocity")

    plt.tight_layout()
    plt.show()
