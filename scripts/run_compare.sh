#!/bin/bash
STEPS=100
echo "version,N,P,steps,total,comp,comm,imbalance,e0,e1,rep" > results/compare.csv
for rep in 1 2 3 4 5; do
  for N in 100 500 1000 5000 10000; do
    for P in 2 4 8; do
      echo "N=$N P=$P rep=$rep"
      echo "inter,$(mpirun --oversubscribe -np $P ./bin/mpi_nbody $N $STEPS 2>/dev/null),$rep" >> results/compare.csv
      echo "3calls,$(mpirun --oversubscribe -np $P ./bin/mpi_nbody_3calls $N $STEPS 2>/dev/null),$rep" >> results/compare.csv
    done
  done
done
echo "DONE"
