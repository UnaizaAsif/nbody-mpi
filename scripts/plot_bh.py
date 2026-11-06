import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = pd.read_csv("results/bh.csv")
med = d.groupby(["version", "N", "P"])[["total", "comp", "comm", "imbalance"]].median().reset_index()

ds = med[med.version == "direct_seq"].set_index("N")["total"]
bs = med[med.version == "bh_seq"].set_index("N")["total"]

# Table 1: sequential direct vs sequential Barnes-Hut
t1 = pd.DataFrame({"direct_seq": ds, "bh_seq": bs})
t1["bh_gain"] = t1.direct_seq / t1.bh_seq
t1.to_csv("results/bh_seq_summary.csv")
print("Sequential: direct vs Barnes-Hut (seconds, median)")
print(t1.round(3).to_string())

# Table 2: MPI Barnes-Hut scaling
m = med[med.version == "bh_mpi"].copy()
m["speedup_vs_bhseq"] = m.N.map(bs) / m.total
m["efficiency"] = m.speedup_vs_bhseq / m.P
m["speedup_vs_direct_seq"] = m.N.map(ds) / m.total
m["comm_frac"] = m.comm / m.total
m.to_csv("results/bh_mpi_summary.csv", index=False)
print("\nMPI Barnes-Hut")
print(m[["N", "P", "total", "comp", "comm", "imbalance", "speedup_vs_bhseq",
         "efficiency", "speedup_vs_direct_seq", "comm_frac"]].round(3).to_string())

# Figure 1: time vs N, direct vs BH (log-log)
plt.figure()
plt.plot(t1.index, t1.direct_seq, "o-", label="direct O(N^2), sequential")
plt.plot(t1.index, t1.bh_seq, "s-", label="Barnes-Hut, sequential")
plt.xscale("log"); plt.yscale("log")
plt.xlabel("N"); plt.ylabel("Time for 100 steps (s)")
plt.title("Direct vs Barnes-Hut (sequential)")
plt.grid(True, which="both", alpha=0.3); plt.legend()
plt.savefig("results/bh_vs_direct.png", dpi=150)

# Figure 2: MPI Barnes-Hut speedup vs P
plt.figure()
for N in sorted(m.N.unique()):
    s = m[m.N == N]
    plt.plot(s.P, s.speedup_vs_bhseq, "o-", label=f"N={N}")
plt.plot([1, 2, 4, 8], [1, 2, 4, 8], "k--", label="ideal")
plt.xlabel("Processes (P)"); plt.ylabel("Speedup vs sequential Barnes-Hut")
plt.title("MPI Barnes-Hut speedup")
plt.grid(True); plt.legend()
plt.savefig("results/bh_mpi_speedup.png", dpi=150)

# Figure 3: compute vs communication for the largest N
N = m.N.max()
s = m[m.N == N].sort_values("P")
plt.figure()
plt.bar(s.P.astype(str), s.comp, label="compute")
plt.bar(s.P.astype(str), s.comm, bottom=s.comp, label="communication")
plt.xlabel("Processes (P)"); plt.ylabel("Time (s)")
plt.title(f"MPI Barnes-Hut: compute vs communication, N={N}")
plt.legend()
plt.savefig("results/bh_breakdown.png", dpi=150)
print("\nSaved bh_vs_direct.png, bh_mpi_speedup.png, bh_breakdown.png")
