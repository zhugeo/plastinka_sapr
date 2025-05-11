import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

df = pd.read_csv(
    "build/Debug/solutionImplicit.csv",
    header=0,
    delimiter=" ",
    names=["t", "x", "y", "T"],
)

times = df["t"].unique()
times.sort()

Tmin = -10  # df["T"].min()
Tmax = 100  # df["T"].max()

fig, ax = plt.subplots(figsize=(10, 8))

plt.colorbar(
    ax.pcolormesh(
        (df["T"],),
        shading="auto",
        vmax=Tmax,
        vmin=Tmin,
        cmap="plasma",
    )
)


def update(frame):
    ax.clear()

    current_time = times[frame]
    frame_data = df[df["t"] == current_time]

    x = frame_data["x"]
    y = frame_data["y"]
    z = frame_data["T"]

    heatmap = ax.scatter(x, y, c=z, cmap="plasma", vmax=Tmax, vmin=Tmin)

    ax.set_title(f"t = {current_time:.2f}")

    return (heatmap,)


ani = FuncAnimation(
    fig, update, frames=len(times), interval=100, blit=False, repeat=True
)

plt.tight_layout()
plt.show()
