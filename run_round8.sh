#!/bin/bash
# pools on the 30 irregular instances for the solvers whose irregular-tuned setting differs from the regular-tuned one
cd "$(dirname "$0")"
printf '%s\n' "POOL 0-15 sqa" "POOL 0-15 saq" "POOL 15-30 sqa" "POOL 15-30 saq" |
  xargs -P 2 -I{} sh -c 'python3 bench_irrtune.py {} >> "results/r8_$(echo {} | tr " " _ ).log" 2>&1'
