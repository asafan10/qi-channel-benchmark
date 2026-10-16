# Pre-specified analysis plan: protocol P6 on the irregular layout at four channels

Written 2026-10-16 about 13:00 (Cairo), before bench_sigma_k4.py produced any result. Committed to the repository before the analysis is run.

Purpose: the irregular-layout comparison of Section 5.10 rests on one setting (jitter 0.3, K = 3, 200 cells). This run repeats it at
K = 4 with every setting transferred unchanged, to test whether its verdicts hold at a second channel count.

Data: the runs of bench_sigma_k4.py only (30 instances, seeds 0-29, make_layout(200, seed, 4, 0.3)), timed on one machine under a
two-process load (seeds 0-14 and 15-29). Pools are fixed in advance (see the script) and are not extended after any result is seen.

Reference per instance: the best plan of any run of bench_sigma_k4.py on that instance, together with the reuse plan; breakout local search
(bench_bls.py, run on the same instances afterwards) is reported against this reference but not added to it, as in Table 12.

Estimator and eligibility: as in PREREG_sigma03_ext.md (exact expected best of N without replacement; an effort is eligible at budget T
only if 1 <= N and m >= 2N on every instance; per solver and budget the eligible effort with the lowest median gap over instances).

Primary test: tuned native SA (T0 = 0.15) against bSB with refinement at T in {1, 3, 10, 30, 100} s. Equivalence: two one-sided Wilcoxon
signed-rank tests on the paired gap differences with margin 0.5 points, Holm over the tested budgets; equivalence at T holds if and only
if the Holm-adjusted p value is below 0.05. Difference tests: two-sided Wilcoxon, Holm over the tested budgets. "Equivalent up to X s"
stops at the first budget that is untested or not equivalent. Either outcome is reported.

Secondary (separate families): tuned native SA against bSB alone (two-sided Wilcoxon, Holm over budgets); frozen native SA against bSB
with refinement (same). Descriptive only: the smallest margin at which equivalence holds at each budget, and breakout local search.
