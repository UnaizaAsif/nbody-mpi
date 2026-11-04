import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = pd.read_csv("results/compare.csv")
med = d.groupby(["version", "N", "P"])[["total", "comp", "comm"]].median().reset_index()
a = med[med.version == "inter"].drop(columns="version")
b = med[med.version == "3calls"].drop(columns="version")
t = a.merge(b, on=["N", "P"], suffixes=("_inter", "_3calls"))
t["comm_ratio"] = t["comm_3calls"] / t["comm_inter"]      # >1 means 3 calls is slower
t["total_ratio"] = t["total_3calls"] / t["total_inter"]
t.to_csv("results/compare_summary.csv", index=False)
print(t[["N", "P", "comm_inter", "comm_3calls", "comm_ratio",
         "total_inter", "total_3calls", "total_ratio"]].round(4).to_string())

fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=False)
for ax, P in zip(axes, [2, 4, 8]):
    s = t[t.P == P]
    ax.plot(s.N, s.comm_3calls, "o-", label="3 calls")
    ax.plot(s.N, s.comm_inter, "s-", label="1 interleaved call")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("N"); ax.set_ylabel("Communication time (s)")
    ax.set_title(f"P = {P}"); ax.grid(True, which="both", alpha=0.3); ax.legend()
plt.tight_layout()
plt.savefig("results/comm_compare.png", dpi=150)
print("Saved results/comm_compare.png")
