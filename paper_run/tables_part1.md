# Final tables (all gaps scored against the single final best-known reference per instance)


## Table 1. Reference drift: how much lower the final best-known is than the original benchmark's best-known

| K | n | mean drop (%) | max drop (%) | instances improved |
|---|---|---|---|---|
| 3 | 20 | 0.00 | 0.00 | 0/5 |
| 3 | 50 | 0.00 | 0.00 | 0/5 |
| 3 | 100 | 1.96 | 3.33 | 4/5 |
| 3 | 200 | 7.50 | 10.41 | 5/5 |
| 3 | 500 | 10.57 | 11.95 | 5/5 |
| 4 | 50 | 0.21 | 1.06 | 1/5 |
| 4 | 100 | 0.15 | 0.64 | 3/5 |
| 4 | 200 | 3.38 | 3.97 | 5/5 |

## Table 2. Fixed-budget (3/6/12 s) restart protocol: gap to final best-known (%), mean ± std over 5 instances


**K=3**

| n | Greedy | Tabu | SA-native | CP-SAT | SA-QUBO | PT-QUBO | bSB | dSB |
|---|---|---|---|---|---|---|---|---|
| 20 | 26.2 ± 2.7 | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.0 ± 0.0 | 0.7 ± 1.5 | 10.1 ± 13.3 | 0.2 ± 0.5 | 0.8 ± 1.0 |
| 50 | 32.3 ± 4.2 | 0.0 ± 0.0 | 1.2 ± 1.7 | 32.3 ± 4.2 | 20.2 ± 3.8 | 32.6 ± 6.9 | 1.5 ± 1.8 | 17.5 ± 4.9 |
| 100 | 33.7 ± 2.4 | 2.2 ± 1.5 | 6.5 ± 1.9 | 33.7 ± 2.4 | 36.5 ± 8.4 | 52.8 ± 10.2 | 7.0 ± 2.5 | 31.6 ± 7.6 |
| 200 | 39.2 ± 3.9 | 12.3 ± 1.5 | 9.6 ± 2.3 | 39.2 ± 3.9 | 43.0 ± 2.8 | 58.9 ± 3.9 | 8.6 ± 1.9 | 33.5 ± 2.4 |
| 500 | 35.3 ± 3.6 | 21.8 ± 2.4 | 12.1 ± 1.6 | 35.3 ± 3.6 | 43.6 ± 2.8 | 61.7 ± 4.0 | 12.6 ± 0.7 | 34.9 ± 2.2 |

**K=4**

| n | Greedy | Tabu | SA-native | CP-SAT | SA-QUBO | PT-QUBO | bSB | dSB |
|---|---|---|---|---|---|---|---|---|
| 50 | 25.6 ± 6.6 | 0.2 ± 0.4 | 3.3 ± 0.9 | 25.6 ± 6.6 | 36.3 ± 5.9 | 52.5 ± 9.6 | 10.0 ± 1.1 | 45.5 ± 7.9 |
| 100 | 22.7 ± 3.0 | 0.2 ± 0.2 | 2.7 ± 1.3 | 22.7 ± 3.0 | 55.2 ± 6.0 | 71.8 ± 10.1 | 8.5 ± 0.9 | 56.7 ± 12.7 |
| 200 | 27.1 ± 1.9 | 3.5 ± 0.5 | 4.1 ± 0.5 | 27.1 ± 1.9 | 59.8 ± 5.1 | 83.7 ± 3.7 | 9.7 ± 0.4 | 67.3 ± 6.1 |

## Table 3. One-hot penalty sweep (K=3, n=100, 6 s per cell): gap to final best-known (%), mean over 5 instances

Native-space SA (6 s restarts): 3.2%

| scale rho | SA-QUBO gap % | bSB gap % | SA-QUBO raw feasibility | bSB raw feasibility |
|---|---|---|---|---|
| 1.0 | 95.2 | 21.7 | 1.00 | 1.00 |
| 0.5 | 66.3 | 10.7 | 1.00 | 1.00 |
| 0.25 | 36.3 | 5.3 | 1.00 | 1.00 |
| 0.12 | 10.4 | 6.0 | 1.00 | 1.00 |
| 0.06 | 8.7 | 5.5 | 0.92 | 0.95 |

## Table 4. Original (as-published, corrected) solvers vs re-implementations (K=3): gap to final best-known (%), 5 instances

| n | solver | gap % | raw one-hot feasibility | time (s) |
|---|---|---|---|---|
| 20 | orig-SBA | 18.8 | 0.00 | 0.0 |
| 20 | orig-DA | 126.2 | 0.81 | 1.3 |
| 20 | orig-PA | 97.4 | 0.50 | 2.9 |
| 20 | new-SA-native | 0.0 | 1.00 | 3.0 |
| 20 | new-bSB | 0.2 | 1.00 | 3.0 |
| 20 | new-SA-QUBO | 0.0 | 1.00 | 3.0 |
| 50 | orig-SBA | 27.8 | 0.00 | 0.0 |
| 50 | orig-DA | 116.6 | 0.85 | 3.4 |
| 50 | orig-PA | 111.4 | 0.47 | 7.7 |
| 50 | new-SA-native | 1.3 | 1.00 | 3.0 |
| 50 | new-bSB | 2.8 | 1.00 | 3.0 |
| 50 | new-SA-QUBO | 19.4 | 1.00 | 3.0 |
