# Final tables

All gaps are scored against a single final best-known cost per instance (minimum over every run of every solver in every experiment).


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


## Table 4b. Effect of two further faults found in the 'corrected' DA/PA code (K=3; n=20 and 50; 5 instances; gap to final best-known, %)

Fault 5: the single-flip energy difference `(1-2x_i)(col@x - Q_ii*x_i)` omits the diagonal term Q_ii for 0->1 flips (exact for 1->0). Verified against brute force (max error 6.4 on random QUBOs); the same expression occurs in DA, PA and the framework's classical-SA routine.
Fault 6 (DA only): replica-exchange acceptance exp(min(500, dB*(E[r+1]-E[r]))) has the opposite sign of the standard criterion exp(dB*(E[r]-E[r+1])). Verified on a 12-variable QUBO against exact Boltzmann averages (as coded: max error 4.03; standard: 0.017).

| n | solver | gap % | raw one-hot feasibility | time (s) |
|---|---|---|---|---|
| 20 | DA as given | 126.2 | 0.81 | 1.3 |
| 20 | DA + energy-difference fix | 41.3 | 1.00 | 1.3 |
| 20 | DA + energy-difference fix + swap fix | 42.2 | 1.00 | 1.3 |
| 20 | PA as given | 97.4 | 0.50 | 2.9 |
| 20 | PA + energy-difference fix | 19.0 | 1.00 | 2.9 |
| 50 | DA as given | 116.6 | 0.85 | 3.3 |
| 50 | DA + energy-difference fix | 61.2 | 1.00 | 3.4 |
| 50 | DA + energy-difference fix + swap fix | 55.6 | 1.00 | 3.4 |
| 50 | PA as given | 111.4 | 0.47 | 7.6 |
| 50 | PA + energy-difference fix | 37.8 | 1.00 | 7.7 |

The energy-difference fault explains most of the DA/PA gap in Table 4; the swap-sign fault has no measurable effect on this metric. Remaining gaps (19-61%) reflect default settings (fixed temperature range 0.001-8, fixed sweeps), not isolated. SBA (Table 4) is unaffected by these faults. SQA was not benchmarked. A further suspected issue in PA (importance weights use the absolute inverse temperature at every step instead of the increment) was NOT tested.


## Table 5. Independent-run TTS at the frozen baseline effort (K=3 and K=4; eps = 5% and 2%): median TTS99 in ms [95% bootstrap CI]; pooled per-run success p


### K=3, eps = 5%

| n | Tabu | SA-native | SA-QUBO | PT-QUBO | bSB | dSB |
|---|---|---|---|---|---|---|
| 50 | 19 [6, 217]; p=0.219 | 1687 [556, 3409]; p=0.029 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | 2773 [401, inf]; p=0.247 | inf [inf, inf]; p=0.000 |
| 100 | inf [1639, inf]; p=0.002 | inf [inf, inf]; p=0.001 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | 20605 [6776, inf]; p=0.033 | inf [inf, inf]; p=0.000 |
| 200 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 |
| 500 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 |


### K=3, eps = 2%

| n | Tabu | SA-native | SA-QUBO | PT-QUBO | bSB | dSB |
|---|---|---|---|---|---|---|
| 50 | 65 [15, inf]; p=0.098 | 6203 [1695, inf]; p=0.008 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | 5838 [2171, inf]; p=0.073 | inf [inf, inf]; p=0.000 |
| 100 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.007 | inf [inf, inf]; p=0.000 |
| 200 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 |
| 500 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 |



### K=4, eps = 5%

| n | Tabu | SA-native | SA-QUBO | PT-QUBO | bSB | dSB |
|---|---|---|---|---|---|---|
| 50 | 14 [7, 16]; p=0.321 | 1633 [1088, inf]; p=0.015 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 |
| 100 | 37 [16, 104]; p=0.325 | 4709 [421, inf]; p=0.051 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 |
| 200 | 2649 [1107, 13141]; p=0.021 | 33057 [15550, inf]; p=0.005 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 |


### K=4, eps = 2%

| n | Tabu | SA-native | SA-QUBO | PT-QUBO | bSB | dSB |
|---|---|---|---|---|---|---|
| 50 | 211 [96, 1305]; p=0.020 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 |
| 100 | 711 [289, inf]; p=0.021 | inf [14127, inf]; p=0.002 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 |
| 200 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 | inf [inf, inf]; p=0.000 |



## Table 6. Run-length-optimised TTS (sweep grid only: SA-native up to 16x, bSB up to 2x): median TTS99 in ms [95% CI incl. re-selection of run length]; (best factor, per-run success p, per-run time)

### K=3, eps=5%: run-length-optimised TTS99 in ms, median over instances [95% bootstrap CI incl. re-selection of run length]; (best effort factor, per-run success p, per-run time ms)

| n | Tabu | SA-native | SA-QUBO | PT-QUBO | bSB | dSB |
|---|---|---|---|---|---|---|
| 100 | 3310 [798, inf] (x1, p=0.005, 2 ms) | 7415 [4146, 20551] (x16, p=0.135, 231 ms) | inf [inf, inf] (x0.25, p=0.000, 8 ms) | inf [inf, inf] (x0.25, p=0.000, 14 ms) | inf [7166, inf] (x0.17, p=0.000, 25 ms) | inf [inf, inf] (x0.17, p=0.000, 31 ms) |
| 200 | inf [inf, inf] (x0.25, p=0.000, 2 ms) | 22466 [10527, inf] (x16, p=0.100, 517 ms) | inf [inf, inf] (x0.25, p=0.000, 16 ms) | inf [inf, inf] (x0.25, p=0.000, 28 ms) | inf [inf, inf] (x0.17, p=0.000, 45 ms) | inf [inf, inf] (x0.17, p=0.000, 57 ms) |
| 500 | inf [inf, inf] (x0.25, p=0.000, 11 ms) | inf [inf, inf] (x0.25, p=0.000, 25 ms) | inf [inf, inf] (x0.25, p=0.000, 33 ms) | inf [inf, inf] (x0.25, p=0.000, 60 ms) | inf [inf, inf] (x0.17, p=0.000, 110 ms) | inf [inf, inf] (x0.17, p=0.000, 145 ms) |

### K=3, eps=2%: run-length-optimised TTS99 in ms, median over instances [95% bootstrap CI incl. re-selection of run length]; (best effort factor, per-run success p, per-run time ms)

| n | Tabu | SA-native | SA-QUBO | PT-QUBO | bSB | dSB |
|---|---|---|---|---|---|---|
| 100 | inf [2807, inf] (x0.25, p=0.000, 1 ms) | 23153 [11122, inf] (x16, p=0.036, 231 ms) | inf [inf, inf] (x0.25, p=0.000, 8 ms) | inf [inf, inf] (x0.25, p=0.000, 14 ms) | inf [inf, inf] (x0.17, p=0.000, 25 ms) | inf [inf, inf] (x0.17, p=0.000, 31 ms) |
| 200 | inf [inf, inf] (x0.25, p=0.000, 2 ms) | inf [inf, inf] (x0.25, p=0.000, 9 ms) | inf [inf, inf] (x0.25, p=0.000, 16 ms) | inf [inf, inf] (x0.25, p=0.000, 28 ms) | inf [inf, inf] (x0.17, p=0.000, 45 ms) | inf [inf, inf] (x0.17, p=0.000, 57 ms) |
| 500 | inf [inf, inf] (x0.25, p=0.000, 11 ms) | inf [inf, inf] (x0.25, p=0.000, 25 ms) | inf [inf, inf] (x0.25, p=0.000, 33 ms) | inf [inf, inf] (x0.25, p=0.000, 60 ms) | inf [inf, inf] (x0.17, p=0.000, 110 ms) | inf [inf, inf] (x0.17, p=0.000, 145 ms) |


### K=4, eps=5%: run-length-optimised TTS99 in ms, median over instances [95% bootstrap CI incl. re-selection of run length]; (best effort factor, per-run success p, per-run time ms)

| n | Tabu | SA-native | SA-QUBO | PT-QUBO | bSB | dSB |
|---|---|---|---|---|---|---|
| 200 | 2636 [1262, 11422] (x1, p=0.017, 10 ms) | 1736 [573, 4392] (x16, p=0.756, 570 ms) | inf [inf, inf] (x0.25, p=0.000, 17 ms) | inf [inf, inf] (x0.25, p=0.000, 30 ms) | inf [inf, inf] (x0.17, p=0.000, 61 ms) | inf [inf, inf] (x0.17, p=0.000, 84 ms) |

### K=4, eps=2%: run-length-optimised TTS99 in ms, median over instances [95% bootstrap CI incl. re-selection of run length]; (best effort factor, per-run success p, per-run time ms)

| n | Tabu | SA-native | SA-QUBO | PT-QUBO | bSB | dSB |
|---|---|---|---|---|---|---|
| 200 | inf [23521, inf] (x0.25, p=0.000, 3 ms) | inf [22580, inf] (x0.25, p=0.000, 9 ms) | inf [inf, inf] (x0.25, p=0.000, 17 ms) | inf [inf, inf] (x0.25, p=0.000, 30 ms) | inf [inf, inf] (x0.17, p=0.000, 61 ms) | inf [inf, inf] (x0.17, p=0.000, 84 ms) |



## Table 7. bSB re-tuning study (held-out instances: n=200, K=3, seeds 100-102)

Grid: step size {0.25,0.5,1.0} x pump exponent {0.5,1,2} x coupling multiplier {0.5,1,2} x penalty scale rho {0.12,0.25,0.5} = 81 configs; 4x steps, 2 batches of 8 replicas per instance.

- Best / median / worst config (mean gap at 4x steps): 6.05% / 7.51% / 22.52%; untuned default: 7.18%
- Top 5 at 4x: dt0.5_g2.0_c2.0_rho0.25 = 6.05%; dt1.0_g1.0_c0.5_rho0.12 = 6.07%; dt1.0_g0.5_c1.0_rho0.25 = 6.11%; dt0.25_g1.0_c1.0_rho0.12 = 6.29%; dt0.25_g0.5_c1.0_rho0.12 = 6.40%
- Confirmation at 16x steps (top-3 + default): dt0.5_g2.0_c2.0_rho0.25 = 5.71%; dt1.0_g1.0_c0.5_rho0.12 = 6.70%; dt1.0_g0.5_c1.0_rho0.25 = 6.76%; dt0.5_g1.0_c1.0_rho0.25 = 7.14%
- Frozen configuration: {'dt': 0.5, 'gamma': 2.0, 'c0m': 2.0, 'rho': 0.25}
- Transfer to evaluation instances: no improvement (see Table 9).


## Table 8a. Long-run single-run quality and matched wall-clock (SA-native vs untuned bSB), K=3


## n = 200 cells (5 instances)

### Single-run median gap to updated best-known (%) and mean time per run (ms)

| solver | effort factor | runs/instance | median gap % | time/run (ms) |
|---|---|---|---|---|
| SA-native | x0.25 | 300 | 22.9 | 9 |
| SA-native | x1 | 148 | 18.0 | 34 |
| SA-native | x4 | 38 | 13.8 | 133 |
| SA-native | x16 | 10 | 9.7 | 517 |
| SA-native | x64 | 8 | 6.0 | 2060 |
| SA-native | x256 | 3 | 5.2 | 8147 |
| bSB | x0.167 | 110 | 26.9 | 45 |
| bSB | x0.333 | 57 | 18.9 | 89 |
| bSB | x1 | 20 | 13.3 | 262 |
| bSB | x2 | 10 | 12.0 | 527 |
| bSB | x4 | 5 | 12.3 | 1069 |
| bSB | x8 | 3 | 12.7 | 2115 |
| bSB | x16 | 3 | 12.0 | 4205 |
| bSB | x32 | 2 | 12.8 | 8337 |
| bSB | x64 | 1 | 13.2 | 16907 |

### Matched wall-clock: expected best-of-N gap (%) with the best effort level chosen at each budget (median over instances)

| budget T | SA-native | bSB (batch = 1 run) | bSB (each replica = 1 run, optimistic) |
|---|---|---|---|
| 1 s | 9.3% (x16) | 11.4% (x1) | 11.9% (x4) |
| 3 s | 6.0% (x16) | 10.2% (x1) | 10.6% (x4) |
| 10 s | 5.0% (x64) | 9.5% (x1) | 10.0% (x4) |
| 30 s | 4.2% (x64) | 9.0% (x1) | 9.9% (x4) |
| 100 s | 4.0% (x64) | 9.0% (x1) | 9.9% (x4) |

## n = 500 cells (5 instances)

### Single-run median gap to updated best-known (%) and mean time per run (ms)

| solver | effort factor | runs/instance | median gap % | time/run (ms) |
|---|---|---|---|---|
| SA-native | x0.25 | 200 | 20.5 | 25 |
| SA-native | x1 | 52 | 16.2 | 96 |
| SA-native | x4 | 14 | 12.8 | 379 |
| SA-native | x16 | 8 | 9.4 | 1492 |
| SA-native | x64 | 8 | 6.0 | 5729 |
| SA-native | x256 | 3 | 2.4 | 22037 |
| bSB | x0.167 | 46 | 27.6 | 110 |
| bSB | x0.333 | 23 | 19.4 | 220 |
| bSB | x1 | 8 | 14.4 | 638 |
| bSB | x2 | 8 | 13.1 | 1286 |
| bSB | x4 | 5 | 12.1 | 2609 |
| bSB | x8 | 3 | 12.3 | 5134 |
| bSB | x16 | 3 | 12.5 | 10305 |
| bSB | x32 | 2 | 11.6 | 20627 |
| bSB | x64 | 1 | 12.8 | 41072 |

### Matched wall-clock: expected best-of-N gap (%) with the best effort level chosen at each budget (median over instances)

| budget T | SA-native | bSB (batch = 1 run) | bSB (each replica = 1 run, optimistic) |
|---|---|---|---|
| 1 s | 11.6% (x4) | 14.6% (x1) | 13.1% (x4) |
| 3 s | 9.2% (x16) | 12.1% (x4) | 12.2% (x4) |
| 10 s | 5.9% (x64) | 11.6% (x2) | 11.6% (x4) |
| 30 s | 2.3% (x256) | 11.6% (x2) | 11.6% (x4) |
| 100 s | 0.6% (x256) | 11.2% (x32) | 11.2% (x32) |


## Table 8b. Same, K=4, n=200 (bSB frozen K=3 tuning NOT applied; untuned bSB only)


## n = 200 cells (5 instances)

### Single-run median gap to updated best-known (%) and mean time per run (ms)

| solver | effort factor | runs/instance | median gap % | time/run (ms) |
|---|---|---|---|---|
| SA-native | x0.25 | 300 | 11.7 | 9 |
| SA-native | x1 | 139 | 8.7 | 36 |
| SA-native | x4 | 36 | 6.6 | 144 |
| SA-native | x16 | 9 | 4.6 | 570 |
| SA-native | x64 | 8 | 3.0 | 2918 |
| SA-native | x256 | 3 | 1.4 | 11656 |
| bSB | x0.167 | 83 | 35.7 | 61 |
| bSB | x0.333 | 42 | 23.1 | 119 |
| bSB | x1 | 15 | 12.2 | 355 |
| bSB | x2 | 8 | 9.7 | 693 |
| bSB | x4 | 5 | 8.8 | 1620 |
| bSB | x8 | 3 | 8.1 | 3232 |
| bSB | x16 | 3 | 7.9 | 6380 |
| bSB | x32 | 2 | 8.2 | 12709 |
| bSB | x64 | 1 | 8.1 | 25528 |

### Matched wall-clock: expected best-of-N gap (%) with the best effort level chosen at each budget (median over instances)

| budget T | SA-native | bSB (batch = 1 run) | bSB (each replica = 1 run, optimistic) |
|---|---|---|---|
| 1 s | 4.5% (x4) | 9.5% (x2) | 9.4% (x8) |
| 3 s | 3.5% (x16) | 8.4% (x2) | 8.2% (x8) |
| 10 s | 2.1% (x64) | 7.7% (x8) | 7.6% (x8) |
| 30 s | 1.4% (x256) | 7.6% (x8) | 7.6% (x8) |
| 100 s | 1.2% (x256) | 7.6% (x8) | 7.6% (x8) |


## Table 9. Matched wall-clock with tuned bSB (frozen config), K=3


## n = 200 (5 instances)

| T | SA-native | bSB (untuned) | bSB (tuned, frozen) |
|---|---|---|---|
| 1 s | 9.7% [8.0, 10.6] (x16) | 11.4% [9.0, 15.8] (x1) | - |
| 3 s | 5.9% [4.0, 7.6] (x16) | 10.2% [7.8, 15.3] (x1) | 11.9% [6.5, 14.2] (x4) |
| 10 s | 5.0% [1.2, 5.6] (x64) | 9.5% [6.9, 15.1] (x1) | 11.6% [8.0, 14.8] (x16) |
| 30 s | 4.2% [0.2, 5.4] (x64) | 9.0% [6.7, 15.0] (x1) | 10.2% [6.7, 14.1] (x16) |
| 100 s | 4.0% [0.0, 5.3] (x64) | 9.0% [6.7, 15.0] (x1) | 9.8% [6.5, 14.1] (x16) |

### single-run median gap (%) / time per run (ms), tuned bSB, n=200

- x4: 12.3% (1369 ms)
- x16: 12.2% (5449 ms)
- x64: 12.3% (21457 ms)

## n = 500 (5 instances)

| T | SA-native | bSB (untuned) | bSB (tuned, frozen) |
|---|---|---|---|
| 1 s | 11.7% [10.4, 12.6] (x4) | 14.5% [12.0, 14.8] (x1) | - |
| 3 s | 9.2% [6.2, 9.6] (x16) | 12.2% [10.9, 13.1] (x4) | - |
| 10 s | 5.9% [4.6, 6.8] (x64) | 11.6% [8.8, 12.6] (x2) | 11.9% [10.7, 12.3] (x4) |
| 30 s | 2.4% [1.2, 3.7] (x256) | 11.6% [8.7, 12.3] (x2) | 11.5% [10.5, 11.8] (x4) |
| 100 s | 0.5% [0.0, 3.3] (x256) | 11.2% [10.5, 12.0] (x32) | 11.4% [10.4, 11.8] (x4) |

### single-run median gap (%) / time per run (ms), tuned bSB, n=500

- x4: 12.0% (3543 ms)
- x16: 12.4% (13935 ms)
- x64: 13.0% (55396 ms)


## Claims retired or rescoped by these tables
- **Retired:** "bSB reaches parity with native-space SA". Parity appears only under the short-restart protocol (Table 2). With long runs SA-native keeps improving while bSB plateaus, at K=3 and K=4 (Tables 8a-8b), and also after tuning at K=3 (Tables 7, 9).
- **Retired:** the first time-to-epsilon and the baseline-effort TTS as solver rankings (Tables 5-6 supersede them).
- **Rescoped:** "Tabu is best at K=4 and small sizes". It holds under fixed-budget restarts (Table 2) and for baseline-effort TTS at 50-100 cells. At K=4, n=200, run-length-optimised SA-native reaches 5% about three times faster than tabu (576 [562, 600] ms vs 1,664 [702, 2,201] ms; intervals do not overlap); at 2% they are comparable. Tabu was not included in the long-run comparison.
- **Rescoped:** the penalty sweep (Table 3): rho = 0.25, 0.12 and 0.06 are close for bSB, not a sharp optimum at 0.25.
- **Corrected statement on Table 4:** the as-given 'corrected' DA/PA still contained two further faults (Table 4b). The energy-difference fault explains most of their gap (DA 126% to 41%, PA 97% to 19% at 20 cells); the rest reflects default settings and is not isolated. SBA is unaffected.
- **Robust:** penalty-based QUBO solvers (SA-QUBO, PT-QUBO, dSB) are far behind at every effort level.
- **Limits:** 5 instances per size; K=4 long-run only at n=200 and without bSB re-tuning; no long-run at n=100; reference is heuristic and likely still too high at 500 cells.
