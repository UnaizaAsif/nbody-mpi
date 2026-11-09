"""Extra figures for the report.

1. results/time_vs_N.png       total time vs N (log-log), sequential and MPI P=1,2,4,8
2. results/decomposition.png   block decomposition for N=1003, P=4 and the Allgatherv step

Run after plot.py (needs results/seq.csv and results/mpi.csv):
    python3 scripts/plot_extra.py
"""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

main_N = [100, 500, 1000, 5000, 10000]

# ---------- 1. Total time vs N ----------
seq = pd.read_csv("results/seq.csv").groupby("N")["total"].median()
mpi = pd.read_csv("results/mpi.csv").groupby(["N", "P"])["total"].median().reset_index()

plt.figure()
plt.plot(main_N, [seq[n] for n in main_N], "k^-", label="sequential")
for P in [1, 2, 4, 8]:
    d = mpi[(mpi.P == P) & (mpi.N.isin(main_N))].sort_values("N")
    plt.plot(d.N, d.total, "o-", label=f"MPI, P={P}")
plt.xscale("log"); plt.yscale("log")
plt.xlabel("Number of bodies N"); plt.ylabel("Total time for 100 steps (s)")
plt.title("Total time vs N (direct method)")
plt.grid(True, which="both", alpha=0.3); plt.legend()
plt.savefig("results/time_vs_N.png", dpi=150)
print("Saved results/time_vs_N.png")

# ---------- 2. Decomposition diagram ----------
N, P = 1003, 4
base, rem = N // P, N % P
counts = [base + (1 if r < rem else 0) for r in range(P)]
displs = [sum(counts[:r]) for r in range(P)]
colors = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]

fig, ax = plt.subplots(figsize=(10, 5.2))
ax.set_xlim(0, 10); ax.set_ylim(0, 5.2); ax.axis("off")

# Global position array split into blocks
ax.text(5, 4.95, f"Position array of all N = {N} bodies (x, y, z interleaved, 24N bytes)",
        ha="center", fontsize=11, weight="bold")
x0, width = 0.5, 9.0
for r in range(P):
    w = width * counts[r] / N
    x = x0 + width * displs[r] / N
    ax.add_patch(plt.Rectangle((x, 4.15), w, 0.55, color=colors[r], alpha=0.85))
    ax.text(x + w / 2, 4.42, f"[{displs[r]}, {displs[r] + counts[r]})\n{counts[r]} bodies",
            ha="center", va="center", fontsize=9, color="white", weight="bold")

# Ranks
for r in range(P):
    cx = 1.0 + r * 2.35
    ax.add_patch(FancyBboxPatch((cx, 1.75), 1.75, 1.45, boxstyle="round,pad=0.05",
                                fc="white", ec=colors[r], lw=2))
    ax.text(cx + 0.875, 2.95, f"Rank {r}", ha="center", fontsize=11, weight="bold", color=colors[r])
    ax.text(cx + 0.875, 2.42, f"updates {counts[r]} bodies\n(own v, a)\nholds all N positions",
            ha="center", va="center", fontsize=8.5)
    xb = x0 + width * (displs[r] + counts[r] / 2) / N
    ax.add_patch(FancyArrowPatch((xb, 4.12), (cx + 0.875, 3.27), arrowstyle="-|>",
                                 mutation_scale=12, color=colors[r], lw=1.5))
    ax.add_patch(FancyArrowPatch((cx + 0.875, 1.72), (cx + 0.875, 1.05), arrowstyle="-|>",
                                 mutation_scale=12, color="0.3", lw=1.5))

# Allgatherv bar
ax.add_patch(FancyBboxPatch((0.9, 0.45), 9.0 - 0.6, 0.55, boxstyle="round,pad=0.05",
                            fc="#EEEEEE", ec="0.3", lw=1.5))
ax.text(5.1, 0.72, "MPI_Allgatherv (MPI_IN_PLACE): every rank receives all updated positions",
        ha="center", va="center", fontsize=10, weight="bold")
ax.text(5.1, 0.12, f"base = N / P = {base},  rem = N mod P = {rem}:  ranks r < rem get base + 1  "
        f"->  counts = {counts}", ha="center", fontsize=9, style="italic")
plt.savefig("results/decomposition.png", dpi=150, bbox_inches="tight")
print("Saved results/decomposition.png")
