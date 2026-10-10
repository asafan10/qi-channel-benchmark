# Irregular-layout rerun: results for checking (before the rebuild)

Plan: `PREREG_sigma03_ext.md` (v1 13:50, amendment 1 14:24, amendment 2 14:27; commits 46209b4, 7dea688). Data: 30 instances, jitter 0.3, 200 cells, K = 3, fresh pools only. Estimator: exact expected best of N without replacement; an effort is eligible at a budget only if m >= 2N on every instance.

## 1. Equivalence table (primary reference)

Rule outcome: **equivalent up to 30 s**; 100 s untested (no eligible effort for either solver). With the semidefinite-rounding plans in the reference (separate family): up to 30 s.

| Comparison | T (s) | median diff (points; negative favours SA) | 90% CI | SA better on | Holm p (difference) | Holm p (TOST, margin 0.5) |
|---|---|---|---|---|---|---|
| SA-tuned − bSB+polish | 1 | -0.24 | [-0.28, -0.15] | 22/30 | 0.00061 | 1.5e-07 |
| SA-tuned − bSB+polish | 3 | -0.15 | [-0.24, -0.02] | 21/30 | 0.31 | 1.5e-07 |
| SA-tuned − bSB+polish | 10 | -0.07 | [-0.19, +0.06] | 17/30 | 1 | 9.3e-08 |
| SA-tuned − bSB+polish | 30 | -0.01 | [-0.04, +0.09] | 16/30 | 1 | 1.5e-07 |
| SA-frozen − bSB+polish | 1 | +0.04 | [-0.06, +0.18] | 13/30 | 0.37 | 2.8e-06 |
| SA-frozen − bSB+polish | 3 | +0.14 | [+0.08, +0.23] | 6/30 | 0.0013 | 2.8e-06 |
| SA-frozen − bSB+polish | 10 | +0.20 | [+0.11, +0.29] | 6/30 | 0.0015 | 1.1e-05 |
| SA-frozen − bSB+polish | 30 | +0.20 | [+0.05, +0.32] | 9/30 | 0.0064 | 6.2e-05 |
| SA-tuned − bSB | 1 | -2.31 | [-2.78, -1.97] | 30/30 | 7.5e-09 | 1 |
| SA-tuned − bSB | 3 | -2.49 | [-2.91, -2.17] | 30/30 | 7.5e-09 | 1 |
| SA-tuned − bSB | 10 | -2.42 | [-3.01, -2.16] | 30/30 | 7.5e-09 | 1 |
| SA-tuned − bSB | 30 | -3.33 | [-3.52, -2.88] | 30/30 | 7.5e-09 | 1 |

Expected best-of-N gap (%, median over instances [95% bootstrap CI], effort chosen):

| Solver | 1 s | 3 s | 10 s | 30 s | 100 s |
|---|---|---|---|---|---|
| SA-tuned | 2.76 [2.67, 3.01] x4 | 2.05 [1.92, 2.28] x16 | 1.46 [1.32, 1.69] x64 | 1.03 [0.84, 1.32] x256 | – |
| SA-frozen | 3.15 [2.94, 3.40] x4 | 2.34 [2.19, 2.56] x16 | 1.82 [1.46, 1.97] x64 | 1.19 [1.01, 1.42] x256 | – |
| bSB | 5.07 [4.78, 5.65] x1 | 4.42 [4.17, 5.02] x1 | 3.90 [3.69, 4.46] x4 | 4.18 [3.76, 4.48] x16 | – |
| bSB+polish | 3.03 [2.89, 3.19] x1+16k | 2.13 [2.01, 2.25] x1+16k | 1.52 [1.31, 1.63] x1+64k | 0.96 [0.84, 1.18] x1+64k | – |

## 2. Pools and per-instance eligibility (m runs stored; N range over instances; instances with N >= 1 and m >= 2N)

| Solver | effort | m | mean run time (s) | 1 s | 3 s | 10 s | 30 s | 100 s |
|---|---|---|---|---|---|---|---|---|
| SA-tuned | x1 | 50 | 0.033 | N 22–36: 5/30 | N 68–110: 0/30 | N 227–367: 0/30 | N 683–1103: 0/30 | N 2278–3678: 0/30 |
| SA-tuned | x4 | 40 | 0.132 | N 6–9: 30/30 | N 18–27: 4/30 | N 61–92: 0/30 | N 185–277: 0/30 | N 617–926: 0/30 |
| SA-tuned | x16 | 34 | 0.520 | N 1–2: 30/30 | N 4–6: 30/30 | N 15–20: 4/30 | N 47–61: 0/30 | N 157–205: 0/30 |
| SA-tuned | x64 | 26 | 2.069 | N 0–0: 0/30 | N 1–1: 30/30 | N 4–5: 30/30 | N 12–15: 4/30 | N 40–52: 0/30 |
| SA-tuned | x256 | 20 | 8.250 | N 0–0: 0/30 | N 0–0: 0/30 | N 1–1: 30/30 | N 3–3: 30/30 | N 10–13: 4/30 |
| SA-frozen | x1 | 50 | 0.037 | N 21–33: 6/30 | N 65–101: 0/30 | N 217–337: 0/30 | N 653–1013: 0/30 | N 2177–3378: 0/30 |
| SA-frozen | x4 | 40 | 0.146 | N 5–8: 30/30 | N 17–24: 16/30 | N 56–80: 0/30 | N 170–240: 0/30 | N 567–801: 0/30 |
| SA-frozen | x16 | 34 | 0.577 | N 1–1: 30/30 | N 4–5: 30/30 | N 14–18: 22/30 | N 43–55: 0/30 | N 146–186: 0/30 |
| SA-frozen | x64 | 26 | 2.317 | N 0–0: 0/30 | N 1–1: 30/30 | N 3–4: 30/30 | N 10–13: 30/30 | N 35–45: 0/30 |
| SA-frozen | x256 | 20 | 9.202 | N 0–0: 0/30 | N 0–0: 0/30 | N 0–1: 26/30 | N 2–3: 30/30 | N 9–11: 5/30 |
| bSB | x1 | 78 | 0.270 | N 2–4: 30/30 | N 8–12: 30/30 | N 29–40: 29/30 | N 88–120: 0/30 | N 296–403: 0/30 |
| bSB | x4 | 48 | 1.044 | N 0–1: 10/30 | N 2–3: 30/30 | N 7–10: 30/30 | N 23–30: 4/30 | N 77–102: 0/30 |
| bSB | x16 | 40 | 4.109 | N 0–0: 0/30 | N 0–0: 0/30 | N 1–2: 30/30 | N 5–7: 30/30 | N 19–26: 4/30 |
| bSB+polish | x1+4k | 78 | 0.392 | N 2–2: 30/30 | N 6–8: 30/30 | N 20–27: 30/30 | N 62–82: 0/30 | N 206–274: 0/30 |
| bSB+polish | x1+16k | 78 | 0.756 | N 1–1: 30/30 | N 3–4: 30/30 | N 10–14: 30/30 | N 32–42: 4/30 | N 108–142: 0/30 |
| bSB+polish | x1+64k | 78 | 2.204 | N 0–0: 0/30 | N 1–1: 30/30 | N 3–4: 30/30 | N 11–14: 30/30 | N 37–48: 4/30 |

## 3. Timing check

- Mean per-run times of the two processes (seeds 0–14 vs 15–29) agree within 2.2% for every effort.
- **But the load drifted over time.** The four instances run first (seeds 0, 1, 15, 16, before the workspace paused at about 14:30) ran 21–30% slower per run than the other 26 (SA x256: 9.70 vs 8.03 s; polish x1+64k: 2.62 vs 2.14 s; bSB x16: 5.14 vs 3.95 s). They overlapped with short diagnostic and test-analysis jobs I ran on the same two cores; the slowdown was not the same for every solver (bSB 30%, SA 21%).
- Consequences: (a) the pools were sized from those slow first instances, so on the other 26 the runs were faster, N at 100 s exceeded m/2, and 100 s became untested; (b) as a post hoc sensitivity check (not in the plan), dropping the four instances leaves every verdict unchanged:
  equivalent up to 30 s; tuned − polish at 30 s -0.01 points (TOST Holm p 1e-06).
- The runs were interrupted once (workspace paused about 14:30 to 15:58) and restarted; every run records its own time, so the pause itself does not enter any measurement.

## 4. Other results that feed the text

- bSB alone: single-run median gap 5.8%, 5.0%, 5.1% at 1, 4, 16 times its steps; best of N 3.9% at 10 s (x4) and 4.2% at 30 s (x16); not estimable at 100 s. The 0.1-point 'still improving' rule cannot be applied (no 100 s value), so the text reports the 10–30 s values without either word.
- Semidefinite rounding at equal time vs tuned SA (exploratory, uncorrected): rounding + descent alone 13/30, +0.15 points, p = 0.42 (no difference); with the 20 SA polish runs 26/30, -0.47 points, p = 6e-05. The polish does the work, so the method is renamed 'SDP-initialised native SA'. Times include the 10^-4 solve on all 30 instances (median 35 s, max 55 s).
- Penalty/interference force ratio (median over spins, our integrator, rho 0.25, 3,000 steps), at steps 0 / 750 / 1500 / 2250 / 3000: regular 1.34 / 0.67 / 0.55 / 2.20 / 2.21; jitter 0.3 1.37 / 0.70 / 0.71 / 3.31 / 3.31.
- Certified gap of the new reference: 6.5% on average.
