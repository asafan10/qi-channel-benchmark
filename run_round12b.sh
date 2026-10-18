#!/bin/sh
# PT-QUBO grid extended below the best upper temperature (PTTUNE4, run alone; tuning is scored in sweeps, not time),
# then the pools of the new best setting re-run under the usual two-process load
cd "$(dirname "$0")"
python3 bench_round12.py PTPOOL 0-15 > results/r12b_p1.log 2>&1 &
python3 bench_round12.py PTPOOL 15-30 > results/r12b_p2.log 2>&1 &
wait
