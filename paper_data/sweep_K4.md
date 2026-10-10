
### K=4, eps=5%: run-length-optimised TTS99 in ms, median over instances [95% bootstrap CI incl. re-selection of run length]; (best effort factor, per-run success p, per-run time ms)

| n | Tabu | SA-native | SA-QUBO | PT-QUBO | bSB | dSB |
|---|---|---|---|---|---|---|
| 200 | 1664 [702, 2201] (x1, p=0.035, 10 ms) | 576 [562, 600] (x16, p=0.978, 570 ms) | inf [inf, inf] (x0.25, p=0.000, 17 ms) | inf [inf, inf] (x0.25, p=0.000, 30 ms) | inf [inf, inf] (x0.17, p=0.000, 61 ms) | inf [inf, inf] (x0.17, p=0.000, 84 ms) |

### K=4, eps=2%: run-length-optimised TTS99 in ms, median over instances [95% bootstrap CI incl. re-selection of run length]; (best effort factor, per-run success p, per-run time ms)

| n | Tabu | SA-native | SA-QUBO | PT-QUBO | bSB | dSB |
|---|---|---|---|---|---|---|
| 200 | 23308 [11393, inf] (x16, p=0.024, 153 ms) | 22519 [10508, inf] (x16, p=0.089, 570 ms) | inf [inf, inf] (x0.25, p=0.000, 17 ms) | inf [inf, inf] (x0.25, p=0.000, 30 ms) | inf [inf, inf] (x0.17, p=0.000, 61 ms) | inf [inf, inf] (x0.17, p=0.000, 84 ms) |

### K=4: median single-run gap to best-known (%) by effort factor (per-run time in ms)


**n=200**

| solver | x0.2500 | x1.0000 | x4.0000 | x16.0000 |
|---|---|---|---|---|
| Tabu (factors 0.25,1,4,16) | 10.7% (3 ms) | 9.3% (10 ms) | 7.9% (38 ms) | 6.3% (153 ms) |
| SA-native (factors 0.25,1,4,16) | 10.5% (9 ms) | 7.3% (36 ms) | 5.3% (144 ms) | 3.6% (570 ms) |
| SA-QUBO (factors 0.25,1,4,16) | 88.9% (17 ms) | 74.8% (67 ms) | 63.7% (268 ms) | 56.6% (1060 ms) |
| PT-QUBO (factors 0.25,1,4,16) | 113.6% (30 ms) | 91.5% (116 ms) | 74.2% (445 ms) | 63.6% (1785 ms) |
| bSB (factors 0.17,0.33,1,2) | 35.0% (61 ms) | 21.5% (119 ms) | 11.2% (355 ms) | 9.0% (693 ms) |
| dSB (factors 0.17,0.33,1,2) | 102.9% (84 ms) | 93.4% (171 ms) | 74.8% (502 ms) | 63.9% (1013 ms) |
