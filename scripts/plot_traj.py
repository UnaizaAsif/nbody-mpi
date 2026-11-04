import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = pd.read_csv("results/trajectory.csv")
last = d.step.max()
fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

for ax, (u, v) in zip(axes, [("x", "y"), ("x", "z")]):
    for _, b in d.groupby("body"):
        ax.plot(b[u], b[v], lw=0.6, alpha=0.5)
    s = d[d.step == 0]
    e = d[d.step == last]
    ax.scatter(s[u], s[v], c="green", s=14, label="start (step 0)", zorder=3)
    ax.scatter(e[u], e[v], c="red",   s=14, label=f"end (step {last})", zorder=3)
    ax.set_xlabel(u); ax.set_ylabel(v)
    ax.set_title(f"{u}-{v} projection, N={d.body.nunique()} bodies")
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, alpha=0.3); ax.legend()

plt.tight_layout()
plt.savefig("results/trajectories.png", dpi=150)
print("Saved results/trajectories.png")
