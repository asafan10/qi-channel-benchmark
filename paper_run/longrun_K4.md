# Symmetric long-run comparison, K=4


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
