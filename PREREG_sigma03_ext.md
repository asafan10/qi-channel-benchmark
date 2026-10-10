# Pre-specified analysis plan for the re-analysis of protocol P6 on the irregular layout

This is a pre-specified analysis plan, not a preregistration: it was written after the referee's critique (round 2) and while the runs were in progress.

Version 1 written 2026-10-10 about 13:50 (Cairo), before bench_sigma_ext.py produced any result.
Amendment 1, 2026-10-10 about 14:30 (Cairo), after two of the 30 instances had finished and after the analysis code had been test-run on those two
instances to check that it executes (no decision depends on two instances): the rules for budgets without an eligible effort, below, were added.
Both versions are committed to the repository with timestamps before the full analysis is run.

Data: the fresh runs of bench_sigma_ext.py only (30 instances, jitter 0.3, 200 cells, K = 3), timed on one machine under a
two-process load. The runs of bench_sigma.py are not mixed in, except through the reference.

Reference per instance: the best plan of any run of bench_sigma.py and bench_sigma_ext.py on it, including the reuse plans
(primary); secondary: the same with the semidefinite-rounding plans added.

Estimator: expected best of N runs, N = floor(T / mean run time of the effort on that instance), computed exactly without replacement
from the m stored runs: E = sum_k f_(k) C(m-k, N-1) / C(m, N), f_(1) <= ... <= f_(m). An effort is eligible at budget T only if
1 <= N and m >= 2N for every instance. Per solver and budget, the eligible effort with the lowest median gap over instances is used.

Primary test: tuned native SA against bSB with refinement at T in {1, 3, 10, 30, 100} s. Two one-sided Wilcoxon signed-rank tests on
the paired gap differences (points of gap) with margin 0.5 points; Holm correction over the five budgets. Equivalence at T holds if and
only if the Holm-adjusted TOST p value is below 0.05. Difference tests: two-sided Wilcoxon, Holm over the five budgets.
The claim "equivalent up to X s" uses the largest X such that equivalence holds at every budget up to X. Either outcome is reported.

Budgets without an eligible effort (amendment 1): if no effort of a solver satisfies m >= 2N on every instance at a budget, that cell is
reported as "–" (not estimable) and the comparisons involving it are not run. Such a budget counts as untested, not as a failure:
"equivalent up to X s" is the largest X such that equivalence holds at every tested budget up to X, and the text states which budgets were
untested. The Holm family is the set of tested budgets of that comparison (at most five).

Separate families: the comparison with the semidefinite-rounding plans in the reference, the cross-fitted variants and the equal-time comparison
with semidefinite rounding are each reported as their own families (the last is exploratory and uncorrected).

Amendment 2, 2026-10-10 about 14:35 (Cairo), before the full analysis, after a referee's comments:
- "Equivalent up to X s" stops at the first budget, in increasing order, at which the comparison is either untested (no eligible effort for
  either solver) or not equivalent. A verdict is reported only at budgets where both solvers have an eligible effort.
- The per-instance eligibility counts (instances with N >= 1 and m >= 2N, per solver, effort and budget) are reported.
- The wording "still improving slowly" for bSB alone is used when its best-of-N gap falls by more than 0.1 point from 30 s to 100 s. This
  threshold was chosen after the first two instances had been seen; it is descriptive and is not tested.
- The per-run times of the two processes (seeds 0-14 and 15-29) are compared per effort, and any drift in load is reported.
