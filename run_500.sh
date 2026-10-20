#!/bin/bash
cd /home/claude/qi
python3 bench_hybrid.py 3 500 results/hybrid_K3.jsonl 10-30 >> results/hybrid_500.log 2>&1
python3 bench_hybrid_ext.py 3 500 results/hybext_K3.jsonl 10-30 >> results/hybrid_500.log 2>&1
echo DONE >> results/hybrid_500.log
