#!/bin/bash
STEPS=100
echo "N,steps,total,e0,e1,rep" > results/direct_again.csv
for rep in 1 2 3; do
  for N in 1000 5000 10000; do
    echo "seq N=$N rep=$rep"
    echo "$(./bin/seq $N $STEPS),$rep" >> results/direct_again.csv
  done
done
echo "DONE"

