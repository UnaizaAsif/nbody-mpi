#!/bin/bash
STEPS=100
for N in 100 1000 1003 5001; do
  ./bin/seq $N $STEPS results/seq_$N.txt > /dev/null
  for P in 1 2 4 8; do
    mpirun --oversubscribe -np $P ./bin/mpi_nbody $N $STEPS results/mpi_${N}_$P.txt > /dev/null 2>&1
    echo -n "N=$N P=$P : "
    python3 scripts/verify.py results/seq_$N.txt results/mpi_${N}_$P.txt
  done
done
