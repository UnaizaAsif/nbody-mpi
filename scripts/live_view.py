"""Real-time viewer for bin/mpi_live.

Reads frames from a file or named pipe while the MPI simulation is running and
redraws the bodies, coloured by the MPI rank that owns them.

    ./scripts/live.sh 500 3000 5 4        # easiest: N=500, 3000 steps, draw every 5, P=4
    python3 scripts/live_view.py PIPE     # what live.sh runs

Close the window to stop watching. Set LIVE_SAVE=file.png to save the last
frame instead of opening a window (for machines without a display).
"""
import os
import struct
import sys
import numpy as np
import matplotlib

SAVE = os.environ.get("LIVE_SAVE")
if SAVE:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    sys.exit("Usage: python3 scripts/live_view.py PIPE_OR_FILE")
f = open(sys.argv[1], "rb")


def read(fmt):
    size = struct.calcsize(fmt)
    data = f.read(size)
    if len(data) < size:
        return None
    return struct.unpack(fmt, data)


N, P = read("2i")
per_rank = read(f"{P}i")
owner = np.repeat(np.arange(P), per_rank)
COLORS = ["#4C9BE8", "#F2A541", "#6BCB77", "#E8615A", "#B28DFF", "#4DD0E1", "#F48FB1", "#C0CA33"]
colors = [COLORS[r % len(COLORS)] for r in owner]
LIM = 2.0

plt.style.use("dark_background")
fig = plt.figure(figsize=(11, 5.5))
ax3 = fig.add_subplot(1, 2, 1, projection="3d")
ax2 = fig.add_subplot(1, 2, 2)
ax3.set_xlim(-LIM, LIM); ax3.set_ylim(-LIM, LIM); ax3.set_zlim(-LIM, LIM)
for a in (ax3.xaxis, ax3.yaxis, ax3.zaxis):
    a.set_pane_color((0.12, 0.12, 0.15, 1.0))
ax3.set_xlabel("x"); ax3.set_ylabel("y"); ax3.set_zlabel("z")
ax2.set_xlim(-LIM, LIM); ax2.set_ylim(-LIM, LIM); ax2.set_aspect("equal")
ax2.set_xlabel("x"); ax2.set_ylabel("y"); ax2.set_title("x-y projection")
sc3 = ax3.scatter([], [], [], s=6, c=[], depthshade=False)
sc2 = ax2.scatter(np.zeros(N), np.zeros(N), s=6, c=colors)
for r in range(P):
    ax2.scatter([], [], c=COLORS[r % len(COLORS)], label=f"rank {r}: {per_rank[r]} bodies")
ax2.legend(loc="upper right", fontsize=8, markerscale=2)
title = fig.suptitle("", fontsize=12)
info = fig.text(0.5, 0.02, "", ha="center", fontsize=10, color="0.8")

frames = 0
while True:
    hdr = read("i")
    if hdr is None or hdr[0] == -1:
        break
    step = hdr[0]
    sim_t, comm_t = read("2d")
    raw = f.read(8 * 3 * N)
    if len(raw) < 8 * 3 * N:
        break
    pos = np.frombuffer(raw, dtype=np.float64).reshape(N, 3)
    frames += 1
    sc3._offsets3d = (pos[:, 0], pos[:, 1], pos[:, 2])
    sc3.set_color(colors)
    ax3.view_init(elev=22, azim=30 + frames * 0.5)
    sc2.set_offsets(pos[:, :2])
    title.set_text(f"Live MPI N-body: N = {N}, P = {P}, step {step}")
    share = 100 * comm_t / sim_t if sim_t > 0 else 0
    info.set_text(f"simulated time {step * 0.01:.2f}   |   wall time {sim_t:.1f} s   |   "
                  f"communication {share:.0f}% of rank 0 time")
    if SAVE:
        continue
    if not plt.fignum_exists(fig.number):
        break
    plt.pause(0.001)

if SAVE:
    fig.savefig(SAVE, dpi=90)
    print(f"Saved {SAVE} after {frames} frames")
else:
    print(f"Stream ended after {frames} frames")
    if plt.fignum_exists(fig.number):
        plt.show()
