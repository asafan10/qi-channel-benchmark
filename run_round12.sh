#!/bin/sh
# two processes under the usual two-process load; results/r12_*
cd "$(dirname "$0")"
( python3 bench_round12.py PTTUNE; python3 bench_round12.py POLSTART 0-15; python3 bench_round12.py PTPOOL 0-15; python3 bench_round12.py LONGPOL 0-15 ) > results/r12_p1.log 2>&1 &
( python3 bench_round12.py POLSTART 15-30; python3 bench_round12.py PTPOOL 15-30; python3 bench_round12.py LONGPOL 15-30 ) > results/r12_p2.log 2>&1 &
wait
