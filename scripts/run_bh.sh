#!/bin/bash
STEPS=100
THETA=0.5
echo "version,N,P,steps,total,comp,comm,imbalance,e0,e1,rep" > results/bh.csv
for rep in 1 2 3; do
  for N in 1000 5000 10000 20000; do
    echo "N=$N rep=$rep"
    s=$(./bin/seq $N $STEPS)
    echo "direct_seq,$N,1,$STEPS,$(echo $s | cut -d, -f3),,,,$(echo $s | cut -d, -f4),$(echo $s | cut -d, -f5),$rep" >> results/bh.csv
    s=$(./bin/seq_bh $N $STEPS $THETA)
    echo "bh_seq,$N,1,$STEPS,$(echo $s | cut -d, -f3),,,,$(echo $s | cut -d, -f4),$(echo $s | cut -d, -f5),$rep" >> results/bh.csv
    for P in 1 2 4 8; do
      echo "bh_mpi,$(mpirun --oversubscribe -np $P ./bin/mpi_bh $N $STEPS $THETA 2>/dev/null),$rep" >> results/bh.csv
    done
  done
done
echo "DONE"
