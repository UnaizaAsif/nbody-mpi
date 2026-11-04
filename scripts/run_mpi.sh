#!/bin/bash
STEPS=100
echo "N,P,steps,total,comp,comm,imbalance,e0,e1,rep" > results/mpi.csv
for rep in 1 2 3 4 5; do
  for N in 100 500 1000 5000 10000 1003 5001; do
    for P in 1 2 4 8; do
      echo "mpi N=$N P=$P rep=$rep"
      echo "$(mpirun --oversubscribe -np $P ./bin/mpi_nbody $N $STEPS 2>/dev/null),$rep" >> results/mpi.csv
    done
  done
done
echo "DONE"
