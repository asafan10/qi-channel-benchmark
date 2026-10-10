# Experiments B-D (reference includes paper's final_reference.json: False)


## Table D. Fixed-budget restart protocol P1 on one machine (K=3): gap to best-known (%), mean ± std over 5 instances [mean restarts]

| n | budget | SQA | SA-native | SA-QUBO | bSB |
|---|---|---|---|---|---|
| 20 | 3 s | 0.0 ± 0.0 [16] | 0.0 ± 0.0 [969] | 1.0 ± 1.9 [472] | 0.2 ± 0.5 [49] |
| 50 | 3 s | 8.3 ± 1.3 [5] | 0.0 ± 0.0 [350] | 22.9 ± 3.5 [147] | 1.6 ± 2.0 [29] |
| 100 | 6 s | 12.2 ± 3.8 [5] | 5.4 ± 2.1 [319] | 33.5 ± 9.4 [147] | 5.2 ± 3.1 [34] |
| 200 | 12 s | 14.2 ± 2.0 [5] | 8.8 ± 1.2 [289] | 42.8 ± 1.5 [140] | 9.0 ± 2.7 [37] |
| 500 | 12 s | 16.4 ± 1.2 [2] | 12.2 ± 1.0 [99] | 44.7 ± 3.1 [67] | 12.2 ± 0.9 [15] |

## Table C. Independent SB implementation (simulated-bifurcation 2.0.0, 128 agents, held-out-tuned): single-run gap (%), median over instances [range], raw feasibility, time per run

| K | n | variant | median gap % [range] | raw feasibility | time/run (s) | our bSB x1 single-run median gap % |
|---|---|---|---|---|---|---|
| 3 | 200 | pkg-bSB | 27.0 [11.7, 30.1] | 0.01 | 3.4 | 13.1 |
| 3 | 200 | pkg-dSB-heated | 41.3 [41.1, 42.4] | 1.00 | 4.5 | 13.1 |
| 3 | 500 | pkg-bSB | 17.5 [11.3, 29.3] | 0.00 | 16.3 | 14.0 |
| 3 | 500 | pkg-dSB-heated | 43.5 [41.6, 45.4] | 1.00 | 21.4 | 14.0 |
| 4 | 200 | pkg-bSB | 1.6 [0.1, 3.4] | 0.00 | 5.5 | 11.6 |
| 4 | 200 | pkg-dSB-heated | 89.3 [81.0, 104.5] | 1.00 | 7.5 | 11.6 |

## Table B. Certified lower bounds and decoder-only baselines: mean over instances

| K | n | best-known f_ref | SDP lower bound | certified gap of f_ref to optimum, mean [max] (%) | repair of empty vector (scan-order greedy), gap % | degree-ordered greedy, gap % |
|---|---|---|---|---|---|---|
| 3 | 20 | 5.17 | 5.12 | 1.0 [2.9] | 18.8 | 26.2 |
| 3 | 50 | 16.51 | 16.33 | 1.0 [2.0] | 27.8 | 32.3 |
| 3 | 100 | 38.10 | 37.38 | 1.9 [4.4] | 18.7 | 33.3 |
| 3 | 200 | 82.75 | 80.90 | 2.3 [4.2] | 23.6 | 38.7 |
| 3 | 500 | 224.95 | - | - | 20.0 | 35.3 |
| 4 | 50 | 8.68 | 8.16 | 6.3 [7.2] | 9.7 | 25.6 |
| 4 | 100 | 21.03 | 19.22 | 9.4 [10.7] | 3.0 | 22.7 |
| 4 | 200 | 47.45 | 42.68 | 11.2 [11.5] | 1.8 | 26.5 |
