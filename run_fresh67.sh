#!/bin/sh
# two processes, one per core, each running all solvers interleaved on its own instances; results/fresh67_*.jsonl
cd "$(dirname "$0")"
( python3 bench_fresh67.py 3 200 0-5 results/fresh67_K3_200.jsonl; python3 bench_fresh67.py 4 200 0-3 results/fresh67_K4_200a.jsonl ) > results/fresh67_p1.log 2>&1 &
( python3 bench_fresh67.py 3 500 0-5 results/fresh67_K3_500.jsonl; python3 bench_fresh67.py 4 200 3-5 results/fresh67_K4_200b.jsonl ) > results/fresh67_p2.log 2>&1 &
wait
