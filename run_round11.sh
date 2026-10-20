#!/bin/sh
# two processes under the usual two-process load; results/r11_*.jsonl
cd "$(dirname "$0")"
( python3 bench_round11.py BSB256 3 200 0-5 results/r11_bsb256.jsonl; python3 bench_round11.py BSB256 4 200 0-5 results/r11_bsb256.jsonl ) > results/r11_p1.log 2>&1 &
( python3 bench_round11.py SAREP 100 0-5 10 results/r11_sarep.jsonl; python3 bench_round11.py BSB256 3 500 0-5 results/r11_bsb256_500.jsonl ) > results/r11_p2.log 2>&1 &
wait
