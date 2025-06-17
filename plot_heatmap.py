import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

df_implicit = pd.read_csv(
    "build/Debug/solutionImplicit.csv",
    header=0,
    delimiter=" ",
    names=["t", "x", "y", "T"],
)

df_explicit = pd.read_csv(
    "build/Debug/solutionExplicit.csv",
    header=0,
    delimiter=" ",
    names=["t", "x", "y", "T"],
)

times = df_implicit["t"].unique()
times.sort()

Tmin = 0  # df["T"].min()
Tmax = 5  # df["T"].max()


fig, ax = plt.subplots(figsize=(10, 8))

plt.colorbar(
    ax.pcolormesh(
        (df_implicit["T"],),
        shading="auto",
        vmax=Tmax,
        vmin=Tmin,
        cmap="plasma",
    )
)


def update(frame):
    ax.clear()

    current_time = times[frame]
    data_implicit = df_implicit[df_implicit["t"] == current_time]
    data_explicit = df_explicit[df_explicit["t"] == current_time]

    x = data_implicit["x"]
    y = data_implicit["y"]
    T_implicit = np.array(data_implicit["T"], dtype=np.float64)
    T_explicit = np.array(data_explicit["T"], dtype=np.float64)

    T_difference = np.abs(T_implicit - T_explicit)

    T_to_display = T_difference

    heatmap = ax.scatter(
        x, y, c=T_to_display, cmap="plasma", vmax=Tmax, vmin=Tmin, marker="s"
    )

    ax.set_title(f"t = {current_time:.2f}")

    return (heatmap,)


ani = FuncAnimation(
    fig, update, frames=len(times), interval=100, blit=False, repeat=True
)

plt.tight_layout()
plt.show()
