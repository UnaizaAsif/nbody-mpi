#!/bin/bash
STEPS=100
THETA=0.5
for N in 1000 1003 5001; do
  ./bin/seq_bh $N $STEPS $THETA results/bhseq_$N.txt > /dev/null
  for P in 1 2 4 8; do
    mpirun --oversubscribe -np $P ./bin/mpi_bh $N $STEPS $THETA results/bhmpi_${N}_$P.txt > /dev/null 2>&1
    echo -n "BH-MPI vs BH-seq  N=$N P=$P : "
    python3 scripts/verify.py results/bhseq_$N.txt results/bhmpi_${N}_$P.txt
  done
  echo -n "BH-seq vs DIRECT  N=$N (tol 1e-2) : "
  python3 scripts/verify.py results/seq_$N.txt results/bhseq_$N.txt 1e-2
done
