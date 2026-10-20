#!/bin/bash
cd /home/claude/qi
python3 bench_hybrid.py 3 200 results/hybrid_K3.jsonl 10-30 >> results/hybrid_more2.log 2>&1
python3 bench_hybrid.py 4 200 results/hybrid_K4.jsonl 10-30 >> results/hybrid_more2.log 2>&1
python3 bench_hybrid_ext.py 3 200 results/hybext_K3.jsonl 10-30 >> results/hybext_more2.log 2>&1
python3 bench_hybrid_ext.py 4 200 results/hybext_K4.jsonl 10-30 >> results/hybext_more2.log 2>&1
echo DONE >> results/hybrid_more2.log
