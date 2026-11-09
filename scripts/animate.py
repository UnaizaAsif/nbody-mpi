"""Animate the N-body trajectories saved by ./bin/traj.

Makes results/nbody.gif: a 3D view (slowly rotating) and the x-y projection
with short trails, one frame per saved step.

    ./bin/traj 100 1000 10          # writes results/trajectory.csv
    python3 scripts/animate.py      # writes results/nbody.gif
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

d = pd.read_csv("results/trajectory.csv").sort_values(["step", "body"])
steps = np.sort(d.step.unique())
N = d.body.nunique()
pos = d[["x", "y", "z"]].to_numpy().reshape(len(steps), N, 3)

LIM = 2.5          # view box; ejected bodies leave it, which is what happens physically
TRAIL = 8          # number of past frames drawn as a trail

fig = plt.figure(figsize=(11, 5.2), facecolor="#0B1D3A")
ax3 = fig.add_subplot(1, 2, 1, projection="3d")
ax2 = fig.add_subplot(1, 2, 2)
for ax in (ax3, ax2):
    ax.set_facecolor("#0B1D3A")
for a in (ax3.xaxis, ax3.yaxis, ax3.zaxis):
    a.set_pane_color((0.04, 0.11, 0.23, 1.0))
    a.label.set_color("white")
ax3.tick_params(colors="0.7", labelsize=7)
ax2.tick_params(colors="0.7", labelsize=8)
for s in ax2.spines.values():
    s.set_color("0.5")

ax3.set_xlim(-LIM, LIM); ax3.set_ylim(-LIM, LIM); ax3.set_zlim(-LIM, LIM)
ax3.set_xlabel("x"); ax3.set_ylabel("y"); ax3.set_zlabel("z")
ax2.set_xlim(-LIM, LIM); ax2.set_ylim(-LIM, LIM); ax2.set_aspect("equal")
ax2.set_xlabel("x", color="white"); ax2.set_ylabel("y", color="white")
ax2.set_title("x-y projection with trails", color="white", fontsize=11)

sc3 = ax3.scatter([], [], [], s=10, c="#F2C14E", depthshade=False)
trails = ax2.plot(np.full((2, N), np.nan), np.full((2, N), np.nan), color="#5B8BD0", lw=0.7, alpha=0.6)
sc2 = ax2.scatter([], [], s=12, c="#F2C14E", zorder=3)
title = fig.suptitle("", color="white", fontsize=13)


def update(k):
    p = pos[k]
    sc3._offsets3d = (p[:, 0], p[:, 1], p[:, 2])
    ax3.view_init(elev=22, azim=30 + k * 1.2)
    sc2.set_offsets(p[:, :2])
    lo = max(0, k - TRAIL)
    for i, line in enumerate(trails):
        line.set_data(pos[lo:k + 1, i, 0], pos[lo:k + 1, i, 1])
    title.set_text(f"N = {N} bodies, step {steps[k]} of {steps[-1]}  (direct O(N²), velocity-Verlet, dt = 0.01)")
    return [sc3, sc2, title, *trails]


anim = FuncAnimation(fig, update, frames=len(steps), interval=80, blit=False)
anim.save("results/nbody.gif", writer=PillowWriter(fps=12), dpi=80)
print("Saved results/nbody.gif")
