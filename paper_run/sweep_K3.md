
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

### K=3: median single-run gap to best-known (%) by effort factor (per-run time in ms)


**n=100**

| solver | x0.2500 | x1.0000 | x4.0000 | x16.0000 |
|---|---|---|---|---|
| Tabu (factors 0.25,1,4,16) | 19.0% (1 ms) | 16.7% (2 ms) | 14.8% (9 ms) | 13.3% (37 ms) |
| SA-native (factors 0.25,1,4,16) | 21.6% (4 ms) | 17.2% (15 ms) | 12.6% (58 ms) | 9.8% (231 ms) |
| SA-QUBO (factors 0.25,1,4,16) | 64.4% (8 ms) | 52.6% (32 ms) | 45.0% (125 ms) | 38.8% (495 ms) |
| PT-QUBO (factors 0.25,1,4,16) | 78.4% (14 ms) | 62.2% (52 ms) | 49.9% (205 ms) | 41.2% (820 ms) |
| bSB (factors 0.17,0.33,1,2) | 25.0% (25 ms) | 17.2% (48 ms) | 11.9% (142 ms) | 9.5% (280 ms) |
| dSB (factors 0.17,0.33,1,2) | 52.0% (31 ms) | 44.5% (62 ms) | 36.8% (177 ms) | 33.8% (353 ms) |

**n=200**

| solver | x0.2500 | x1.0000 | x4.0000 | x16.0000 |
|---|---|---|---|---|
| Tabu (factors 0.25,1,4,16) | 24.1% (2 ms) | 23.1% (8 ms) | 21.1% (31 ms) | 20.2% (124 ms) |
| SA-native (factors 0.25,1,4,16) | 22.9% (9 ms) | 18.0% (34 ms) | 13.8% (133 ms) | 9.7% (517 ms) |
| SA-QUBO (factors 0.25,1,4,16) | 63.5% (16 ms) | 53.1% (63 ms) | 46.2% (253 ms) | 40.0% (1011 ms) |
| PT-QUBO (factors 0.25,1,4,16) | 79.0% (28 ms) | 64.9% (106 ms) | 52.5% (414 ms) | 45.9% (1658 ms) |
| bSB (factors 0.17,0.33,1,2) | 26.9% (45 ms) | 18.9% (89 ms) | 13.3% (262 ms) | 12.0% (527 ms) |
| dSB (factors 0.17,0.33,1,2) | 53.7% (57 ms) | 46.6% (114 ms) | 39.9% (334 ms) | 35.6% (676 ms) |

**n=500**

| solver | x0.2500 | x1.0000 | x4.0000 | x16.0000 |
|---|---|---|---|---|
| Tabu (factors 0.25,1,4,16) | 28.8% (11 ms) | 28.4% (44 ms) | 28.8% (176 ms) | 28.1% (695 ms) |
| SA-native (factors 0.25,1,4,16) | 20.5% (25 ms) | 16.2% (96 ms) | 12.8% (379 ms) | 9.4% (1492 ms) |
| SA-QUBO (factors 0.25,1,4,16) | 59.8% (33 ms) | 51.3% (132 ms) | 44.5% (526 ms) | 38.6% (2117 ms) |
| PT-QUBO (factors 0.25,1,4,16) | 79.5% (60 ms) | 66.3% (224 ms) | 55.0% (880 ms) | 46.8% (3558 ms) |
| bSB (factors 0.17,0.33,1,2) | 27.6% (110 ms) | 19.4% (220 ms) | 14.4% (638 ms) | 13.1% (1286 ms) |
| dSB (factors 0.17,0.33,1,2) | 54.1% (146 ms) | 45.9% (297 ms) | 38.7% (871 ms) | 34.7% (1689 ms) |
