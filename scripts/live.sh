#!/bin/bash
# Real-time visualization of the MPI simulation.
# Usage: ./scripts/live.sh [N] [STEPS] [EVERY] [P]
# Example: ./scripts/live.sh 500 3000 5 4
N=${1:-500}
STEPS=${2:-3000}
EVERY=${3:-5}
P=${4:-4}
PIPE=/tmp/nbody_live_$$
mkfifo "$PIPE" || exit 1
trap 'rm -f "$PIPE"' EXIT

python3 scripts/live_view.py "$PIPE" &
VIEWER=$!
mpirun --oversubscribe -np "$P" ./bin/mpi_live "$N" "$STEPS" "$EVERY" "$PIPE"
wait $VIEWER
