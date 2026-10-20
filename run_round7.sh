#!/bin/bash
# bSB at penalty 0.06 and at its tuned 0.12 on the instances of Table 7; tuned SQA and bit-flip SA on the 30 irregular
# instances; always two processes
cd "$(dirname "$0")"
printf '%s\n' "BSBP 0-5 3 0.06" "BSBP 0-5 4 0.06" "BSBP 0-5 3 0.12" "BSBP 0-5 4 0.12" "IRR 0-15 sqa" "IRR 15-30 sqa" "IRR 0-15 saq" "IRR 15-30 saq" |
  xargs -P 2 -I{} sh -c 'python3 bench_lowpen.py {} >> "results/r7_$(echo {} | tr " " _ ).log" 2>&1'
