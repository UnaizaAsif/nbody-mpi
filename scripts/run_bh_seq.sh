#!/bin/bash
STEPS=100
THETA=0.5
echo "version,N,total,rep" > results/bh_seq_check.csv
for rep in 1 2 3; do
  for N in 1000 5000 10000; do
    echo "N=$N rep=$rep"
    s=$(./bin/seq $N $STEPS);                echo "direct,$N,$(echo $s | cut -d, -f3),$rep" >> results/bh_seq_check.csv
    s=$(./bin/seq_bh $N $STEPS $THETA);      echo "bh,$N,$(echo $s | cut -d, -f3),$rep" >> results/bh_seq_check.csv
  done
done
echo "DONE"
