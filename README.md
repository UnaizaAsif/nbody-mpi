# Distributed N-Body Simulation (MPI)

CT-488 Distributed Computing, Complex Computing Activity.
Group members: Laiba Minhaj (CT-23001), Wania Masood (CT-23002), Unaiza Asif (CT-23008)

## 1. What this project does

Simulates N bodies attracting each other by gravity using the direct O(N^2)
method. There is a sequential baseline and an MPI version that splits the
bodies across processes. Both use the same physics code (`src/common.c`),
and the MPI results are verified against the sequential results.

As an advanced extension, the project also contains a Barnes-Hut octree
version (O(N log N), approximate) in both sequential and MPI form
(`src/bh.c`, `src/seq_bh.c`, `src/mpi_bh.c`). The direct method stays the
baseline for all correctness and speedup measurements.

## 2. Environment used for the reported results

- CPU: Intel Core i5-8350U @ 1.70 GHz, 4 vCPUs (VMware Workstation)
- RAM: 3.8 GiB
- OS: Ubuntu 24.04 LTS
- Compiler: gcc 13.3.0, flags `-O3 -march=native -Wall`
- MPI: Open MPI 4.1.6
- Python: pandas, matplotlib (plots only)

## 3. Install

```bash
sudo apt update
sudo apt install -y build-essential openmpi-bin libopenmpi-dev \
    python3-pip python3-pandas python3-matplotlib
```

## 4. Build

```bash
./scripts/build.sh
```

This builds seven programs in `bin/`:

| Program | Purpose |
|---|---|
| `seq` | sequential direct baseline |
| `mpi_nbody` | MPI direct version (positions sent in one interleaved Allgatherv) |
| `mpi_nbody_3calls` | MPI variant (x, y, z in separate arrays, 3 Allgatherv calls) |
| `traj` | saves trajectories for plotting |
| `seq_bh` | sequential Barnes-Hut |
| `mpi_bh` | MPI Barnes-Hut (same decomposition and Allgatherv as `mpi_nbody`) |
| `bh_accuracy` | force error of Barnes-Hut against the direct method for several theta |

## 5. Run

```bash
./bin/seq N STEPS [outfile]
mpirun --oversubscribe -np P ./bin/mpi_nbody N STEPS [outfile]
mpirun --oversubscribe -np P ./bin/mpi_nbody_3calls N STEPS [outfile]

./bin/seq_bh N STEPS THETA [outfile]
mpirun --oversubscribe -np P ./bin/mpi_bh N STEPS THETA [outfile]
./bin/bh_accuracy N
./bin/traj N STEPS EVERY [outfile]
```

`outfile` (optional) receives the final positions, used by the verification scripts.
`THETA` is the Barnes-Hut opening angle (0.5 in all experiments).
The MPI programs print the block owned by each rank to stderr at startup.

Example:

```bash
./bin/seq 1000 100
mpirun --oversubscribe -np 4 ./bin/mpi_nbody 1000 100
```

Output (CSV, one line):
- seq: `N,steps,time,energy_start,energy_end`
- mpi (`mpi_nbody`, `mpi_nbody_3calls`, `mpi_bh`): `N,P,steps,total,compute_max,comm_max,imbalance,energy_start,energy_end`
- seq_bh: same as seq

`--oversubscribe` is needed when P is larger than the number of cores (P=8 on 4 vCPUs).

## 6. Simulation parameters (fixed in `src/common.h`)

| Parameter | Value |
|---|---|
| G | 1.0 |
| time step dt | 0.01 |
| softening epsilon | 0.05 |
| random seed | 42 |
| body mass | 1/N |
| initial positions | uniform in [-1, 1]^3 |
| initial velocities | uniform in [-0.1, 0.1] per component |
| integrator | velocity-Verlet (leapfrog) |
| steps used in experiments | 100 |
| Barnes-Hut opening angle theta | 0.5 |

## 7. Data decomposition

Bodies are split into contiguous blocks. With base = N / P and rem = N % P,
the first `rem` ranks get one extra body. Example N=1003, P=4: 251, 251, 251, 250.
Every rank holds all positions; each rank updates only its own block, then
`MPI_Allgatherv` shares the new positions (the "v" variant handles unequal blocks).

## 8. Verify correctness

```bash
./scripts/verify.sh
```

Compares final positions of the MPI runs against the sequential run for
N = 100, 1000, 1003, 5001 and P = 1, 2, 4, 8. Tolerance: 1e-9.
Expected: all lines PASS (observed max difference: 0).

```bash
./scripts/verify_bh.sh
```

Checks MPI Barnes-Hut against sequential Barnes-Hut (tolerance 1e-9, observed 0)
for N = 1000, 1003, 5001 and P = 1, 2, 4, 8, and sequential Barnes-Hut against
the direct method (tolerance 1e-2, because Barnes-Hut is approximate;
observed 1.5e-3 to 2.2e-3).

## 9. Reproduce the experiments

```bash
./scripts/run_seq.sh        # sequential timings, 5 repetitions
./scripts/run_mpi.sh        # MPI timings, N x P grid, 5 repetitions
python3 scripts/plot.py     # summary.csv + speedup/efficiency/breakdown/imbalance graphs

./scripts/run_compare.sh    # 1 interleaved call vs 3 calls
python3 scripts/compare.py  # compare_summary.csv + comm_compare.png

./bin/traj 100 1000 10
python3 scripts/plot_traj.py   # trajectories.png
python3 scripts/plot_extra.py  # time_vs_N.png + decomposition.png

# Barnes-Hut extension
./bin/bh_accuracy 1000         # force error table for several theta
./scripts/run_bh_seq.sh        # direct vs Barnes-Hut, sequential, back to back
./scripts/run_bh.sh            # direct seq, BH seq, BH MPI (P=1,2,4,8), N up to 20000
python3 scripts/plot_bh.py     # bh_seq_summary.csv, bh_mpi_summary.csv + 3 graphs
```

Timings on a shared VM drift between sessions, so every speedup is computed
from runs measured in the same session (back to back). For Barnes-Hut MPI use
`speedup_vs_bhseq`; `speedup_vs_direct_seq` mixes in a slow direct baseline and
should not be quoted.

Run each script once, with no other heavy programs running. Medians over
repetitions are reported. Raw logs are in `results/` (`seq.csv`, `mpi.csv`,
`compare.csv`, `bh.csv`, `bh_seq_check.csv`).

## 10. Folder structure

```
nbody/
  src/        common.h common.c seq.c mpi_nbody.c mpi_nbody_3calls.c traj.c
              bh.h bh.c seq_bh.c mpi_bh.c bh_accuracy.c
  scripts/    build.sh verify.sh verify.py run_seq.sh run_mpi.sh
              run_compare.sh plot.py compare.py plot_traj.py plot_extra.py
              verify_bh.sh run_bh_seq.sh run_bh.sh run_direct_again.sh plot_bh.py
  bin/        compiled programs
  results/    CSV logs and graphs
  README.md
```
