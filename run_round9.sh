#!/bin/bash
# pools on the 30 irregular instances for bSB with its full P5 grid tuned on the irregular held-out instances
cd "$(dirname "$0")"
printf '%s\n' "POOL 0-15 bsbp5" "POOL 15-30 bsbp5" |
  xargs -P 2 -I{} sh -c 'python3 bench_irrtune.py {} >> "results/r9_$(echo {} | tr " " _ ).log" 2>&1'
