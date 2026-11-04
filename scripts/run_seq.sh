#!/bin/bash
STEPS=100
echo "N,steps,total,e0,e1,rep" > results/seq.csv
for rep in 1 2 3 4 5; do
  for N in 100 500 1000 5000 10000 1003 5001; do
    echo "seq N=$N rep=$rep"
    echo "$(./bin/seq $N $STEPS),$rep" >> results/seq.csv
  done
done
echo "DONE"
