# Distributed N-Body Simulation (MPI)

CT-488 Distributed Computing, Complex Computing Activity.
Group members: Laiba Minhaj (CT-23001), Wania Masood (CT-23002), Unaiza Asif (CT-23008)

## 1. What this project does

Simulates N bodies attracting each other by gravity using the direct O(N^2)
method. There is a sequential baseline and an MPI version that splits the
bodies across processes. Both use the same physics code (`src/common.c`),
and the MPI results are verified against the sequential results.

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

This builds four programs in `bin/`:

| Program | Purpose |
|---|---|
| `seq` | sequential baseline |
| `mpi_nbody` | MPI version (positions sent in one interleaved Allgatherv) |
| `mpi_nbody_3calls` | MPI variant (x, y, z in separate arrays, 3 Allgatherv calls) |
| `traj` | saves trajectories for plotting |

## 5. Run

```bash
./bin/seq N STEPS [outfile]
mpirun --oversubscribe -np P ./bin/mpi_nbody N STEPS [outfile]
```

Example:

```bash
./bin/seq 1000 100
mpirun --oversubscribe -np 4 ./bin/mpi_nbody 1000 100
```

Output (CSV, one line):
- seq: `N,steps,time,energy_start,energy_end`
- mpi: `N,P,steps,total,compute_max,comm_max,imbalance,energy_start,energy_end`

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

## 9. Reproduce the experiments

```bash
./scripts/run_seq.sh        # sequential timings, 5 repetitions
./scripts/run_mpi.sh        # MPI timings, N x P grid, 5 repetitions
python3 scripts/plot.py     # summary.csv + speedup/efficiency/breakdown/imbalance graphs

./scripts/run_compare.sh    # 1 interleaved call vs 3 calls
python3 scripts/compare.py  # compare_summary.csv + comm_compare.png

./bin/traj 100 1000 10
python3 scripts/plot_traj.py   # trajectories.png
```

Run each script once, with no other heavy programs running. Medians over
repetitions are reported. Raw logs are in `results/` (`seq.csv`, `mpi.csv`,
`compare.csv`).

## 10. Folder structure

```
nbody/
  src/        common.h common.c seq.c mpi_nbody.c mpi_nbody_3calls.c traj.c
  scripts/    build.sh verify.sh verify.py run_seq.sh run_mpi.sh
              run_compare.sh plot.py compare.py plot_traj.py
  bin/        compiled programs
  results/    CSV logs and graphs
  README.md
```
