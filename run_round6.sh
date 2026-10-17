#!/bin/bash
# tuned SA-QUBO (TUNE, then POOLs) and breakout-local-search pools, always two processes
cd "$(dirname "$0")"
printf '%s\n' "saq_tuned.py TUNE" "bls.py P12 0-15 3" "saq_tuned.py POOL 0-5 3" "saq_tuned.py POOL 0-5 4" "bls.py P12 15-30 3" "bls.py P12 0-15 4" "bls.py P12 15-30 4" |
  xargs -P 2 -I{} sh -c 'python3 bench_{} >> "results/round6_$(echo {} | tr " " _ ).log" 2>&1'
