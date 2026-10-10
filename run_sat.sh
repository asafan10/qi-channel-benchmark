#!/bin/bash
cd /home/claude/qi
python3 bench_sa_tuned.py 3 200 results/satuned_K3.jsonl 0-30 >> results/satuned.log 2>&1
python3 bench_sa_tuned.py 4 200 results/satuned_K4.jsonl 0-30 >> results/satuned.log 2>&1
python3 bench_sa_tuned.py 3 500 results/satuned_K3.jsonl 0-10 >> results/satuned.log 2>&1
echo DONE >> results/satuned.log
