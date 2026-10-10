# Experiment A: native-space refinement after bSB decoding

Reference: minimum over every run of every experiment available (refs.py). Same machine for all timings.


## K=3, n=200 (30 instances)

### Single-run median gap (%) and mean time per run (ms)

| solver | effort | median gap % | time/run (ms) | runs/instance |
|---|---|---|---|---|
| SA-native | x1 | 17.7 | 41 | 40 |
| SA-native | x4 | 13.7 | 158 | 20 |
| SA-native | x16 | 10.0 | 620 | 10 |
| SA-native | x64 | 6.2 | 2440 | 5 |
| SA-native | x256 | 3.1 | 9649 | 3 |
| SA-tuned | x1 | 16.5 | 35 | 40 |
| SA-tuned | x4 | 12.8 | 136 | 20 |
| SA-tuned | x16 | 8.8 | 520 | 10 |
| SA-tuned | x64 | 5.8 | 2020 | 5 |
| SA-tuned | x256 | 2.1 | 7877 | 3 |
| bSB | x1 | 13.2 | 317 | 6 |
| bSB | x4 | 12.1 | 1267 | 4 |
| bSB | x16 | 11.4 | 4913 | 2 |
| bSB+descent | x1 | 13.1 | 317 | 6 |
| bSB+descent | x4 | 12.1 | 1267 | 4 |
| bSB+descent | x16 | 11.4 | 4913 | 2 |
| bSB+polish | x1+1k | 12.3 | 348 | 6 |
| bSB+polish | x1+4k | 11.0 | 439 | 6 |
| bSB+polish | x1+16k | 8.0 | 785 | 6 |
| bSB+polish | x4+1k | 11.5 | 1298 | 4 |
| bSB+polish | x4+4k | 10.4 | 1390 | 4 |
| bSB+polish | x4+16k | 8.8 | 1755 | 4 |
| bSB+polish | x16+1k | 10.1 | 4944 | 2 |
| bSB+polish | x16+4k | 8.8 | 5032 | 2 |
| bSB+polish | x16+16k | 7.5 | 5385 | 2 |
| bSB+polish | x1+16k (ext) | 8.9 | 805 | 6 |
| bSB+polish | x1+64k (ext) | 5.8 | 2181 | 6 |

### Matched wall-clock: expected best-of-N gap (%), median over instances [95% bootstrap CI], best effort per budget

| T | SA-native | SA-tuned | bSB | bSB+descent | bSB+polish |
|---|---|---|---|---|---|
| 1 s | 9.8 [8.9, 10.6] (x16) | 8.7 [7.8, 9.4] (x16) | 12.4 [12.0, 12.8] (x1) | 12.4 [12.0, 12.6] (x1) | 7.9 [7.0, 8.6] (x1+16k) |
| 3 s | 5.9 [5.4, 6.9] (x64) | 5.9 [5.1, 6.4] (x64) | 11.6 [10.9, 12.2] (x1) | 11.6 [10.8, 12.1] (x1) | 5.8 [5.4, 6.6] (x1+64k (ext)) |
| 10 s | 4.0 [3.2, 5.0] (x64) | 2.6 [1.5, 3.5] (x256) | 10.9 [10.4, 11.8] (x16) | 10.9 [10.4, 11.8] (x16) | 4.2 [3.1, 4.7] (x1+64k (ext)) |
| 30 s | 1.8 [1.0, 2.3] (x256) | 1.3 [0.8, 2.6] (x256) | 10.5 [9.9, 11.4] (x16) | 10.5 [9.9, 11.4] (x16) | 3.1 [2.3, 4.1] (x1+64k (ext)) |
| 100 s | 0.3 [0.0, 1.1] (x256) | 0.7 [0.3, 2.0] (x256) | 10.5 [9.9, 11.3] (x16) | 10.5 [9.8, 11.3] (x16) | 3.0 [2.1, 4.0] (x1+64k (ext)) |

## K=3, n=500 (30 instances)

### Single-run median gap (%) and mean time per run (ms)

| solver | effort | median gap % | time/run (ms) | runs/instance |
|---|---|---|---|---|
| SA-native | x1 | 15.5 | 116 | 40 |
| SA-native | x4 | 11.5 | 449 | 20 |
| SA-native | x16 | 8.1 | 1743 | 10 |
| SA-native | x64 | 5.0 | 6733 | 5 |
| SA-native | x256 | 1.7 | 26045 | 3 |
| SA-tuned | x1 | 14.6 | 104 | 40 |
| SA-tuned | x4 | 11.0 | 397 | 20 |
| SA-tuned | x16 | 7.3 | 1510 | 10 |
| SA-tuned | x64 | 4.1 | 5743 | 5 |
| SA-tuned | x256 | 1.7 | 21904 | 3 |
| bSB | x1 | 13.2 | 752 | 6 |
| bSB | x4 | 11.4 | 2974 | 4 |
| bSB | x16 | 10.9 | 11757 | 2 |
| bSB+descent | x1 | 13.0 | 752 | 6 |
| bSB+descent | x4 | 11.4 | 2974 | 4 |
| bSB+descent | x16 | 10.9 | 11757 | 2 |
| bSB+polish | x1+1k | 11.6 | 848 | 6 |
| bSB+polish | x1+4k | 9.6 | 1124 | 6 |
| bSB+polish | x1+16k | 7.5 | 2182 | 6 |
| bSB+polish | x4+1k | 10.2 | 3069 | 4 |
| bSB+polish | x4+4k | 9.2 | 3344 | 4 |
| bSB+polish | x4+16k | 7.0 | 4396 | 4 |
| bSB+polish | x16+1k | 10.1 | 11851 | 2 |
| bSB+polish | x16+4k | 8.9 | 12125 | 2 |
| bSB+polish | x16+16k | 6.6 | 13172 | 2 |
| bSB+polish | x1+16k (ext) | 7.4 | 2184 | 6 |
| bSB+polish | x1+64k (ext) | 4.6 | 6157 | 6 |

### Matched wall-clock: expected best-of-N gap (%), median over instances [95% bootstrap CI], best effort per budget

| T | SA-native | SA-tuned | bSB | bSB+descent | bSB+polish |
|---|---|---|---|---|---|
| 1 s | 10.4 [9.8, 11.4] (x4) | 10.1 [9.3, 10.8] (x4) | 13.0 [12.2, 13.7] (x1) | 12.9 [12.2, 13.7] (x1) | 11.4 [10.8, 12.0] (x1+1k) |
| 3 s | 7.8 [7.5, 8.7] (x16) | 7.0 [6.5, 7.7] (x16) | 12.1 [11.5, 13.1] (x1) | 12.1 [11.4, 13.1] (x1) | 7.2 [6.6, 8.2] (x1+16k (ext)) |
| 10 s | 4.8 [4.0, 5.5] (x64) | 4.4 [3.7, 5.0] (x64) | 10.8 [10.0, 11.8] (x4) | 10.8 [10.0, 11.8] (x4) | 4.5 [4.1, 5.4] (x1+64k (ext)) |
| 30 s | 1.8 [1.5, 2.5] (x256) | 1.7 [1.2, 1.9] (x256) | 10.6 [9.9, 11.5] (x4) | 10.6 [9.9, 11.5] (x4) | 3.6 [2.7, 4.2] (x1+64k (ext)) |
| 100 s | 0.8 [0.5, 1.8] (x256) | 0.6 [0.3, 0.9] (x256) | 10.5 [9.5, 11.1] (x16) | 10.5 [9.5, 11.1] (x16) | 3.3 [2.4, 3.8] (x1+64k (ext)) |

## K=4, n=200 (30 instances)

### Single-run median gap (%) and mean time per run (ms)

| solver | effort | median gap % | time/run (ms) | runs/instance |
|---|---|---|---|---|
| SA-native | x1 | 8.1 | 43 | 40 |
| SA-native | x4 | 5.9 | 168 | 20 |
| SA-native | x16 | 4.1 | 668 | 10 |
| SA-native | x64 | 2.8 | 2651 | 5 |
| SA-native | x256 | 1.2 | 10598 | 3 |
| SA-tuned | x1 | 8.3 | 37 | 40 |
| SA-tuned | x4 | 5.3 | 143 | 20 |
| SA-tuned | x16 | 3.4 | 562 | 10 |
| SA-tuned | x64 | 1.8 | 2256 | 5 |
| SA-tuned | x256 | 0.8 | 8897 | 3 |
| bSB | x1 | 11.7 | 405 | 6 |
| bSB | x4 | 8.4 | 1575 | 4 |
| bSB | x16 | 7.9 | 6227 | 2 |
| bSB+descent | x1 | 10.3 | 405 | 6 |
| bSB+descent | x4 | 8.2 | 1575 | 4 |
| bSB+descent | x16 | 7.8 | 6227 | 2 |
| bSB+polish | x1+1k | 6.9 | 444 | 6 |
| bSB+polish | x1+4k | 5.2 | 557 | 6 |
| bSB+polish | x1+16k | 3.6 | 1013 | 6 |
| bSB+polish | x4+1k | 7.1 | 1614 | 4 |
| bSB+polish | x4+4k | 5.1 | 1727 | 4 |
| bSB+polish | x4+16k | 3.7 | 2177 | 4 |
| bSB+polish | x16+1k | 7.2 | 6265 | 2 |
| bSB+polish | x16+4k | 5.1 | 6379 | 2 |
| bSB+polish | x16+16k | 3.7 | 6833 | 2 |
| bSB+polish | x1+16k (ext) | 3.6 | 1032 | 6 |
| bSB+polish | x1+64k (ext) | 2.3 | 2868 | 6 |

### Matched wall-clock: expected best-of-N gap (%), median over instances [95% bootstrap CI], best effort per budget

| T | SA-native | SA-tuned | bSB | bSB+descent | bSB+polish |
|---|---|---|---|---|---|
| 1 s | 4.2 [3.9, 4.3] (x16) | 3.4 [3.2, 3.6] (x16) | 11.3 [10.9, 11.5] (x1) | 9.8 [9.6, 10.1] (x1) | 5.2 [4.9, 5.6] (x1+4k) |
| 3 s | 2.7 [2.4, 3.1] (x64) | 1.9 [1.7, 2.1] (x64) | 8.4 [7.9, 8.8] (x4) | 8.3 [7.8, 8.7] (x4) | 2.8 [2.4, 2.9] (x1+16k (ext)) |
| 10 s | 2.0 [1.8, 2.4] (x64) | 0.8 [0.6, 1.0] (x256) | 7.5 [7.2, 8.1] (x4) | 7.4 [6.9, 7.9] (x4) | 1.7 [1.4, 1.9] (x1+64k (ext)) |
| 30 s | 1.2 [0.8, 1.3] (x256) | 0.4 [0.2, 0.6] (x256) | 7.3 [7.1, 8.0] (x4) | 7.2 [6.8, 7.8] (x4) | 1.2 [1.0, 1.5] (x1+64k (ext)) |
| 100 s | 0.6 [0.3, 1.0] (x256) | 0.2 [0.0, 0.3] (x256) | 7.3 [7.1, 8.0] (x4) | 7.2 [6.8, 7.8] (x4) | 1.1 [0.9, 1.4] (x1+64k (ext)) |
