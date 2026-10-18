# Tables 6 and 7: before and after the fresh P4 run

The earlier tables are those of 12 October 2026, 13:32 +03:00, the last version before the run was restarted; their code and data are those of commit d7b722c (12 October 2026, 08:30 +03:00), since nothing but the text changed in between.
The current tables come from commit dc54af3 (13 October 2026, 04:26 +03:00), which added the fresh pools (results/fresh67_*.jsonl) and their summary.
Every later commit carries the same fresh tables. The fresh run (12 October 23:56 to 13 October 03:30, machine B, two processes, one per core) used per-instance seeds, but pool sizes and the order of draws depend on measured run times, so a rerun reproduces the tables statistically, not bit for bit. An earlier attempt that stopped after one instance is kept separately (results/fresh67_aborted) and was not pooled.

## Before the fresh run (data of commit d7b722c)

```
Table 6. Median single-run gap to best-known (%) and mean time per run, by effort.
  -----------------------------------------------------------------------------------------------------------
  Setting          Solver       × 1             × 4             × 16            × 64            × 256
  ---------------- ----------- --------------- --------------- --------------- --------------- --------------
  K = 3, n = 200   SA-native   18.0 (34 ms)    13.8 (133 ms)   9.7 (0.5 s)     6.0 (2.1 s)     5.2 (8.1 s)
                   bSB         13.3 (262 ms)   12.3 (1.1 s)    12.0 (4.2 s)    13.2 (16.9 s)   not run
  K = 3, n = 500   SA-native   16.2 (96 ms)    12.8 (379 ms)   9.4 (1.5 s)     6.0 (5.7 s)     2.4 (22.0 s)
                   bSB         14.4 (638 ms)   12.1 (2.6 s)    12.5 (10.3 s)   12.8 (41.1 s)   not run
  K = 4, n = 200   SA-native   8.7 (36 ms)     6.6 (144 ms)    4.6 (0.6 s)     3.0 (2.9 s)     1.4 (11.7 s)
                   bSB         12.2 (355 ms)   8.8 (1.6 s)     7.9 (6.4 s)     8.1 (25.5 s)    not run
  -----------------------------------------------------------------------------------------------------------
Multiples of the baseline: 1,000 SA sweeps, 3,000 bSB steps. Five instances.

Table 7. Matched wall-clock comparison according to Eq. (22).
  -----------------------------------------------------------------------
  Budget            K = 3, n = 200    K = 3, n = 500    K = 4, n = 200
  ----------------- ----------------- ----------------- -----------------
  1 s               9.7 / 11.4        11.7 / 14.5       4.6 / 9.5
  3 s               5.9 / 10.2        9.2 / 12.2        3.6 / 8.4
  10 s              5.0 / 9.5         5.9 / 11.6        2.1 / 7.6
  30 s              4.2 / 9.0         2.4 / 11.6        1.3 / 7.6
  100 s             4.0 / 9.0         0.5 / 11.2        1.2 / 7.6
  -----------------------------------------------------------------------
Expected best-of-N gap to best-known (%), native SA / untuned bSB, median over five instances, with the best effort chosen at each budget.
```

## After the fresh run (commit dc54af3)

```
Table 6. Median single-run gap to best-known (%) and mean time per run, by effort.
  -------------------------------------------------------------------------------------------------------
  K, n     Solver           × 1             × 4             × 16            × 64            × 256
  -------- --------------- --------------- --------------- --------------- --------------- --------------
  3, 200   SA (tuned T₀)   –               12.7 (125 ms)   9.4 (487 ms)    6.0 (1.9 s)     3.2 (7.6 s)
  3, 200   SA (T₀ = 2)     –               13.5 (145 ms)   9.9 (568 ms)    6.5 (2.2 s)     3.7 (9.0 s)
  3, 200   bSB             13.7 (297 ms)   12.1 (1.1 s)    12.4 (4.6 s)    12.6 (18.7 s)   –
  3, 200   bSB (tuned)     –               12.0 (1.2 s)    12.3 (4.9 s)    –               –
  3, 500   SA (tuned T₀)   –               11.4 (374 ms)   8.6 (1.4 s)     5.6 (5.5 s)     2.2 (21.6 s)
  3, 500   SA (T₀ = 2)     –               12.0 (417 ms)   8.5 (1.6 s)     5.4 (6.2 s)     3.6 (24.9 s)
  3, 500   bSB             14.3 (722 ms)   12.0 (2.8 s)    12.4 (11.3 s)   12.2 (45.7 s)   –
  3, 500   bSB (tuned)     –               12.7 (3.2 s)    12.8 (12.6 s)   –               –
  4, 200   SA (tuned T₀)   –               5.8 (131 ms)    3.6 (516 ms)    2.3 (2.1 s)     1.3 (8.3 s)
  4, 200   SA (T₀ = 2)     –               6.3 (150 ms)    4.4 (599 ms)    3.0 (2.4 s)     1.6 (9.7 s)
  4, 200   bSB             12.3 (363 ms)   8.5 (1.4 s)     8.1 (5.6 s)     7.6 (22.4 s)    –
  -------------------------------------------------------------------------------------------------------
Protocol P4, fresh run pools on one machine.
Multiples of the baseline: 1,000 SA sweeps, 3,000 bSB steps; five instances; –: not run.

Table 7. Matched wall-clock comparison according to Eq. (22).
  -------------------------------------------------------------------------
  Budget     K = 3, n = 200       K = 3, n = 500       K = 4, n = 200
  ---------- -------------------- -------------------- --------------------
  1 s        8.3 / 10.0 / 12.2    10.5 / 11.2 / 14.2   3.6 / 4.5 / 11.6
  3 s        5.8 / 6.6 / 10.9     7.5 / 8.7 / 12.2     2.2 / 3.1 / 8.0
  10 s       3.0 / 4.2 / 10.1     5.4 / 5.5 / 11.4     1.3 / 1.7 / 7.3
  30 s       0.7 / 2.5 / 11.1     2.6 / 3.3 / 10.6     0.6 / 1.1 / 7.0
  100 s      −0.3 / 0.8 / 11.1    0.9 / 2.1 / 11.2     −0.2 / 0.7 / 6.7
  -------------------------------------------------------------------------
Exact expected best-of-N gap to best-known (%), native SA with the tuned and with the frozen starting temperature / untuned bSB, median over five instances, best eligible effort at each budget (intervals and tuned bSB in Table S7).
[]
Figure 1. Matched wall-clock comparison at K = 3 (200 and 500 cells): exact expected best-of-N gap to best-known against wall-clock budget for native SA with the tuned and the frozen starting temperature, untuned bSB and tuned bSB, with 95% bootstrap bands over five instances; only efforts with at least 2N stored runs are used.
[]
Figure 2. Matched wall-clock comparison at K = 4 (200 cells). Same conventions as Figure 1; bSB is untuned.
Two caveats apply. First, the best-known cost includes native SA’s own long runs, so SA values close to 0% at large budgets are partly limited by the reference and should be read as a trend; Section 5.9 repeats the tests against a reference without those runs. Second, the timings are single-core CPU measurements; bSB is dominated by sparse matrix products and could be faster on a GPU, which we did not measure.
5.3 Re-tuning bSB
To test whether the plateau is an artefact of untuned parameters, we tuned bSB on held-out instances (protocol P5; details and Tables S6 and S7 in Additional file 1). Over an 81-configuration grid of step size, pump exponent, coupling multiplier and penalty scale, the best configuration improved on the default by about one point on the tuning instances, but the improvement did not transfer: on the evaluation instances the tuned and untuned variants overlap within their intervals at every budget, and the tuned single-run gap stays at 12.0–12.8% at 4 and 16 times the baseline steps. A cross-validation on ten further held-out instances (seeds 110–119, 200 cells, 12,000 steps, two batches each, scored against the better of two 64,000-sweep native SA runs) shows that the tuning did not generalise: the five best grid configurations reached 4.8–6.1% against 4.8% for the default, the best tuning-set configuration beat the default on 5 of 10 instances, and the order of the five changed (Spearman correlation 0.20 between tuning and validation scores, which with five configurations is only indicative). With two or three tuning instances per method, tuning is underpowered here, so this shows that we found no setting that lifts the plateau, not that none exists; Section 5.11 shows that on unconstrained Max-Cut the same code plateaus only near the best cut. Other SB variants, such as thermal-fluctuation versions, were not tested.
5.4 Time to solution
Run-length-optimised time to solution (protocols P2 and P3, Eq. (21); Table S8 in Additional file 1) agrees with the matched comparison but is unstable at our scale: from 200 cells a single run rarely comes within 5% of best-known, so Sections 5.2 and 5.6 compare solvers by the expected best of N runs, Eq. (22).
5.5 Penalty design
Table S13 in Additional file 1 sweeps the one-hot penalty scale ρ at 100 and 200 cells (6 s and 12 s of restarts per solver and instance, common seeds, one machine). The penalty has a first-order effect on penalty-based solvers: SA on the QUBO falls from a gap of 95.8% at ρ = 1 (the conventional choice satisfying Proposition 1) to 39.3% at ρ = 0.25, the value used in the main experiments, and to 7.4% at ρ = 0.06 at 100 cells (109.9%, 42.0% and 10.1% at 200 cells). bSB falls from 22.6% to 6.1% at ρ = 0.25 and stays at 5.3–5.8% below it, while native SA, which has no penalty, reaches 3.9% in the same time. Raw feasibility stays 1.00 down to ρ = 0.12 and is 0.91–0.95 at ρ = 0.06. Because the SA-QUBO schedule starts at 0.5λ, the sweep measures penalty and schedule together; tuned for each penalty, bit-flip SA at 200 cells spans 9.5–35.2% (Section 5.11), so part of the factor of ten is the schedule. Under this short-restart protocol bSB at its best penalty is level with native SA within the seed noise (8.0% against 9.7% at 200 cells), as in Section 5.1, and the parity disappears with longer runs (Section 5.2). The main experiments keep ρ = 0.25, fixed before this sweep, so their results for the penalty-based solvers are a lower bound on what these solvers reach with a smaller penalty; for bSB the choice matters little, since below ρ = 0.25 its gap changed by less than one point. At ρ = 0.25 and 500 cells the penalty-based solvers remain far behind under P1 (43.6% for SA-QUBO, 61.7% for PT-QUBO and 34.9% for dSB). Comparisons between penalty-based solvers are therefore only meaningful when the penalty and the schedule are reported and tuned.
5.6 Refinement after decoding
Table 8 and Figure 3 report protocol P6 on 30 instances at 200 cells and 30 at 500 cells, with every solver re-run on one machine. The 1-opt descent leaves bSB essentially unchanged: the decoded state of the selected replica is already a 1-opt local minimum, and the descent makes on average 0.5–1.6 moves at K = 3 and 0.5–5.8 at K = 4. bSB therefore stops in 1-opt local minima that are still about 10% above the best plan, which suggests that what it lacks is not local improvement but moves between distant plans, which a finite-temperature native search supplies (Section 6.2 tests one form of this). A native-space anneal started from the decoded state removes most of the plateau. At 100 s the gap of bSB alone is 10.5% (200 cells) and 10.5% (500 cells) at K = 3 and 7.3% at K = 4, whereas bSB followed by the anneal reaches 3.0%, 3.3% and 1.1%.
Against native SA the picture depends on the budget. Up to 10 s the combination and native SA are level within their intervals at K = 3 (for example 4.2 [3.1, 4.7] against 4.0 [3.2, 5.0] at 200 cells and 10 s). At K = 4 the combination is ahead of native SA with its frozen schedule at 10 s (1.7 [1.4, 1.9] against 2.0 [1.8, 2.4]), a descriptive difference, because the frozen schedule's effort there fails the pool-size rule of Section 5.9; tuned native SA is significantly better at 10 s, so the advantage disappears when native SA’s starting temperature is tuned as carefully as the polish. At K = 3 and 30 s or more, native SA alone has the lower median (0.8 [0.5, 1.8] against 3.3 [2.4, 3.8] at 500 cells and 100 s), and Section 5.9 tests these differences up to 10 s; at 30 s and more the estimates are close to the best stored run and are descriptive. A per-run comparison shows how much the bSB stage contributes: bSB at baseline steps plus 64,000 polish sweeps reaches a median single-run gap of 5.8% at 200 cells and 4.6% at 500 cells in 2.2 s and 6.2 s, against 6.2% and 5.0% for native SA with 64,000 sweeps in 2.4 s and 6.7 s; at K = 4 the figures are 2.3% and 2.8%. The bSB stage therefore provides a modest head start for the anneal, worth about half a point to a point per run at comparable time.
With fresh pools, Table 7 and Table 8 can be compared at long budgets: native SA with the frozen schedule at K = 3, 200 cells and 100 s reaches 0.8% on the five instances of Table 7, against 0.3% over the 30 instances of Table 8, whose long-budget values rest on few stored runs (Section 5.9).
Table 8. Refinement after decoding (protocol P6).
  -------------------------------------------------------------------------------------------------
  K    n     T       SA-native          bSB                 bSB + descent       bSB + polish
  ---- ----- ------- ------------------ ------------------- ------------------- -------------------
  3    200   1 s     9.8 [8.9, 10.6]    12.4 [12.0, 12.8]   12.4 [12.0, 12.6]   7.9 [7.0, 8.6]
  3    200   3 s     5.9 [5.4, 6.9]     11.6 [10.9, 12.2]   11.6 [10.8, 12.1]   5.8 [5.4, 6.6]
  3    200   10 s    4.0 [3.2, 5.0]     10.9 [10.4, 11.8]   10.9 [10.4, 11.8]   4.2 [3.1, 4.7]
  3    200   30 s    1.8 [1.0, 2.3]     10.5 [9.9, 11.4]    10.5 [9.9, 11.4]    3.1 [2.3, 4.1]
  3    200   100 s   0.3 [0.0, 1.1]     10.5 [9.9, 11.3]    10.5 [9.8, 11.3]    3.0 [2.1, 4.0]
  3    500   1 s     10.4 [9.8, 11.4]   13.0 [12.2, 13.7]   12.9 [12.2, 13.7]   11.4 [10.8, 12.0]
  3    500   3 s     7.8 [7.5, 8.7]     12.1 [11.5, 13.1]   12.1 [11.4, 13.1]   7.2 [6.6, 8.2]
  3    500   10 s    4.8 [4.0, 5.5]     10.8 [10.0, 11.8]   10.8 [10.0, 11.8]   4.5 [4.1, 5.4]
  3    500   30 s    1.8 [1.5, 2.5]     10.6 [9.9, 11.5]    10.6 [9.9, 11.5]    3.6 [2.7, 4.2]
  3    500   100 s   0.8 [0.5, 1.8]     10.5 [9.5, 11.1]    10.5 [9.5, 11.1]    3.3 [2.4, 3.8]
  4    200   1 s     4.2 [3.9, 4.3]     11.3 [10.9, 11.5]   9.8 [9.6, 10.1]     5.2 [4.9, 5.6]
  4    200   3 s     2.7 [2.4, 3.1]     8.4 [7.9, 8.8]      8.3 [7.8, 8.7]      2.8 [2.4, 2.9]
  4    200   10 s    2.0 [1.8, 2.4]     7.5 [7.2, 8.1]      7.4 [6.9, 7.9]      1.7 [1.4, 1.9]
  4    200   30 s    1.2 [0.8, 1.3]     7.3 [7.1, 8.0]      7.2 [6.8, 7.8]      1.2 [1.0, 1.5]
  4    200   100 s   0.6 [0.3, 1.0]     7.3 [7.1, 8.0]      7.2 [6.8, 7.8]      1.1 [0.9, 1.4]
  -------------------------------------------------------------------------------------------------
Expected best-of-N gap to best-known (%) at matched wall-clock budget T, median over 30 instances at 200 cells and 30 at 500 cells [95% bootstrap interval], with the best effort chosen at each budget. All solvers on one machine. Every 100 s cell and most 30 s cells rest on efforts with fewer than 2N stored runs on some instance and approximate the best stored run rather than a best of N (Section 5.9).
```

