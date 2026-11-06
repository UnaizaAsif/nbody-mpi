#!/bin/bash
set -e
FLAGS="-O3 -march=native -Wall"
gcc   $FLAGS -o bin/seq        src/seq.c        src/common.c -lm
mpicc $FLAGS -o bin/mpi_nbody  src/mpi_nbody.c  src/common.c -lm
mpicc $FLAGS -o bin/mpi_nbody_3calls src/mpi_nbody_3calls.c src/common.c -lm
gcc   $FLAGS -o bin/traj       src/traj.c       src/common.c -lm
gcc   $FLAGS -o bin/seq_bh      src/seq_bh.c      src/bh.c src/common.c -lm
gcc   $FLAGS -o bin/bh_accuracy src/bh_accuracy.c src/bh.c src/common.c -lm
mpicc $FLAGS -o bin/mpi_bh      src/mpi_bh.c      src/bh.c src/common.c -lm
echo "Build OK"
