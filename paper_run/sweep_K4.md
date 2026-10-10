
### K=4, eps=5%: run-length-optimised TTS99 in ms, median over instances [95% bootstrap CI incl. re-selection of run length]; (best effort factor, per-run success p, per-run time ms)

| n | Tabu | SA-native | SA-QUBO | PT-QUBO | bSB | dSB |
|---|---|---|---|---|---|---|
| 200 | 2636 [1262, 11422] (x1, p=0.017, 10 ms) | 1736 [573, 4392] (x16, p=0.756, 570 ms) | inf [inf, inf] (x0.25, p=0.000, 17 ms) | inf [inf, inf] (x0.25, p=0.000, 30 ms) | inf [inf, inf] (x0.17, p=0.000, 61 ms) | inf [inf, inf] (x0.17, p=0.000, 84 ms) |

### K=4, eps=2%: run-length-optimised TTS99 in ms, median over instances [95% bootstrap CI incl. re-selection of run length]; (best effort factor, per-run success p, per-run time ms)

| n | Tabu | SA-native | SA-QUBO | PT-QUBO | bSB | dSB |
|---|---|---|---|---|---|---|
| 200 | inf [23521, inf] (x0.25, p=0.000, 3 ms) | inf [22580, inf] (x0.25, p=0.000, 9 ms) | inf [inf, inf] (x0.25, p=0.000, 17 ms) | inf [inf, inf] (x0.25, p=0.000, 30 ms) | inf [inf, inf] (x0.17, p=0.000, 61 ms) | inf [inf, inf] (x0.17, p=0.000, 84 ms) |

### K=4: median single-run gap to best-known (%) by effort factor (per-run time in ms)


**n=200**

| solver | x0.2500 | x1.0000 | x4.0000 | x16.0000 |
|---|---|---|---|---|
| Tabu (factors 0.25,1,4,16) | 11.4% (3 ms) | 10.1% (10 ms) | 8.5% (38 ms) | 6.8% (153 ms) |
| SA-native (factors 0.25,1,4,16) | 11.7% (9 ms) | 8.7% (36 ms) | 6.6% (144 ms) | 4.6% (570 ms) |
| SA-QUBO (factors 0.25,1,4,16) | 88.9% (17 ms) | 74.8% (67 ms) | 63.7% (268 ms) | 58.1% (1060 ms) |
| PT-QUBO (factors 0.25,1,4,16) | 114.4% (30 ms) | 92.9% (116 ms) | 74.5% (445 ms) | 63.6% (1785 ms) |
| bSB (factors 0.17,0.33,1,2) | 35.7% (61 ms) | 23.1% (119 ms) | 12.2% (355 ms) | 9.7% (693 ms) |
| dSB (factors 0.17,0.33,1,2) | 106.4% (84 ms) | 95.7% (171 ms) | 77.1% (502 ms) | 64.7% (1013 ms) |
