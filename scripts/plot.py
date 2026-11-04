import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

seq = pd.read_csv("results/seq.csv").groupby("N")["total"].median().reset_index()
seq = seq.rename(columns={"total": "tseq"})
mpi = pd.read_csv("results/mpi.csv").groupby(["N", "P"]).median(numeric_only=True).reset_index()
mpi = mpi.merge(seq, on="N")
mpi["speedup"] = mpi["tseq"] / mpi["total"]
mpi["efficiency"] = mpi["speedup"] / mpi["P"]
mpi["comm_frac"] = mpi["comm"] / mpi["total"]
mpi.to_csv("results/summary.csv", index=False)
print(mpi[["N","P","total","comp","comm","imbalance","speedup","efficiency","comm_frac"]].round(4).to_string())

main_N = [100, 500, 1000, 5000, 10000]
Ps = [1, 2, 4, 8]

# 1. Speedup
plt.figure()
for N in main_N:
    d = mpi[mpi.N == N]; plt.plot(d.P, d.speedup, "o-", label=f"N={N}")
plt.plot(Ps, Ps, "k--", label="ideal")
plt.xlabel("Processes (P)"); plt.ylabel("Speedup"); plt.legend(); plt.grid(True)
plt.title("Speedup vs processes"); plt.savefig("results/speedup.png", dpi=150)

# 2. Efficiency
plt.figure()
for N in main_N:
    d = mpi[mpi.N == N]; plt.plot(d.P, d.efficiency, "o-", label=f"N={N}")
plt.axhline(1.0, color="k", ls="--")
plt.xlabel("Processes (P)"); plt.ylabel("Efficiency"); plt.legend(); plt.grid(True)
plt.title("Parallel efficiency"); plt.savefig("results/efficiency.png", dpi=150)

# 3. Compute vs communication (stacked), one figure per N
for N in main_N:
    d = mpi[mpi.N == N].sort_values("P")
    plt.figure()
    plt.bar(d.P.astype(str), d.comp, label="compute")
    plt.bar(d.P.astype(str), d.comm, bottom=d.comp, label="communication")
    plt.xlabel("Processes (P)"); plt.ylabel("Time (s)"); plt.legend()
    plt.title(f"Compute vs communication, N={N}")
    plt.savefig(f"results/breakdown_N{N}.png", dpi=150)

# 4. Load imbalance (even vs uneven N)
plt.figure()
for N in [1000, 1003, 5000, 5001]:
    d = mpi[mpi.N == N]; plt.plot(d.P, d.imbalance, "o-", label=f"N={N}")
plt.xlabel("Processes (P)"); plt.ylabel("max/avg compute time"); plt.legend(); plt.grid(True)
plt.title("Load imbalance"); plt.savefig("results/imbalance.png", dpi=150)
print("Graphs saved in results/")
