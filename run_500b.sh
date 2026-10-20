#!/bin/bash
cd /home/claude/qi
python3 bench_sa_tuned.py 3 500 results/satuned_K3.jsonl 10-30 >> results/satuned_500.log 2>&1
echo DONE >> results/satuned_500.log
