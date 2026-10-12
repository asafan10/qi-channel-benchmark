# Benchmarking quantum-inspired Ising solvers for wireless channel planning

Code and raw results for A. Saafan, "Benchmarking Ballistic Simulated Bifurcation on a Penalty-QUBO Formulation of Wireless Channel Planning: Domain Baselines, Fair Protocols and Partial Certification" (submitted to *EPJ Quantum Technology*). Table and figure numbers below are those of the submitted manuscript (Tables B1–B9 are in Appendix B).

The G-set files `gset/G1`, `gset/G22` and `gset/G32` are the standard instances in the usual edge-list format, converted from the public copies in the eco-dqn repository (github.com/tomdbar/eco-dqn, `_graphs/benchmarks`); G1 was checked edge by edge against an independent copy (github.com/Open-Finance-Lab/RLSolver). Their best-known cuts (11,624, 13,359, 1,410) are the published values.

All instances are regenerated from their seeds (`bench_ran_qubo.make_instance`), so no instance files are needed. Every raw run is in `results/` and `paper_data/`; the analysis scripts read only those files.

## Environment

- Protocols P1–P5: Intel Core i7-5500U at 2.40 GHz, 16 GB RAM, 64-bit Windows, Python 3.14.
- Protocols P6–P11 and all checks of Sections 5.5–5.10: Intel Xeon at 2.10 GHz, 2 cores, 8 GB RAM, Python 3.13, NumPy 2.5, SciPy 1.18, Numba 0.68, OR-Tools 9.15, PyTorch 2.14 (CPU), simulated-bifurcation 2.0.0, CVXPY 1.9 with SCS. One single-threaded process per core.
- `pip install -r requirements.txt`. Import OR-Tools before CVXPY (a HiGHS library clash otherwise; the scripts already do this).

## Where each result comes from

| Paper item | Experiment script(s) | Raw output | Analysis |
|---|---|---|---|
| Table 2 (instances; Table 1 summarises the literature) | `bench_ran_qubo.py` (`make_instance`) | — | `paper_run/compile_tables.py` |
| Table 3 (solver settings) | `bench_ran_qubo.py` (`CFG`) | `paper_data/frozen_cfg.json` | — |
| Tables S1 and S2 (Supplementary Section S1; faults of the earlier implementation) | `make_fixed_variants.py` (writes `qrc_fix1.py`, `qrc_fix2.py`), `asp_fix.py` (run inside `paper_run/`), `de_check.py`, `da_swap_check.py`; audited code: `quantum_ran_complete.py`, `quantum_inspired_ran.py`, `quantum_optimisation.py` | `paper_run/asp_fix.json` | `paper_run/compile_tables.py` |
| Tables 5–7, B1–B4, Figures 1–2 (protocols P1–P5) | `bench_ran_qubo.py` (P1), `bench_tts.py` (P2), `bench_sweep.py` (P3), `bench_longrun.py` (P4), tuning grid (P5) | `paper_data/*.jsonl`, `paper_data/tuning_grid.json` | `update_paper_ref.py`, then in `paper_run/`: `final_tables.py`, `analyse_longrun.py`, `analyse_tuned.py`, `analyse_tuned_ieee.py`, `compile_tables.py` |
| Table 8 (penalty sweeps, 100 and 200 cells) | `bench_round4.py` (100 cells), `bench_round5.py` part A (200 cells) | `results/p1_100.jsonl`, `results/pen200.jsonl` | `docx_work/build_revised.py` |
| Section 5.3 (cross-validated bSB tuning) | `bench_round5.py` part B | `results/bsb_cv.jsonl` | `docx_work/r5_edits.py` |
| Table 9, Figure 3 (refinement after decoding, P6) | `pilot_polish.py`, `bench_hybrid.py`, `bench_hybrid_ext.py`, `bench_sa_tuned.py` | `results/hybrid_K*.jsonl`, `hybext_K*.jsonl`, `satuned_K*.jsonl` | `analyse_hybrid.py` |
| Table 10 (SQA, P8) | `sqa.py`, `test_sqa.py` (brute-force check), `sqa_pilot.py`, `bench_sqa.py` | `results/sqa_K3.jsonl` | `docx_work/build_revised.py` (`sqa_row`) |
| Section 5.7 (independent SB package, P8) | `sbpkg_pilot.py`, `bench_sbpkg.py`, `bench_core.py` (diagnosis), `bench_round4.py` (bisection) | `results/sbpkg_*.jsonl`, `core_sbdiag.jsonl` | `analyse_extra.py` |
| Table 11 (dual-corrected semidefinite bounds, P9) | `bound_probe.py` (CP-SAT bound), `bench_sdp.py` (SCS), `bench_sdp_rig.py` (dual-corrected bound, SCS tolerance 1e-6 up to 200 cells), `run_sdp500b.py` (500 cells, tolerance 1e-4) | `results/sdp_*.jsonl`, `results/sdp_rig.jsonl` | `refs.py` (`sdp_bounds`) |
| Table 12 (decoder-only and domain baselines, P7, P10) | `decoder_floor.py`, `reuse_baseline.py`, `bench_reuse_polish.py` | `results/decoder_floor.json`, `reuse.json`, `reuse_polish.jsonl` | `refs.py` (`reuse_data`) |
| Table 13, Section 5.9 (tuned SA, paired Wilcoxon tests with Holm correction, weaker reference, work and depth) | `bench_core.py` (tuning of the starting temperature), `bench_sa_tuned.py` | `results/core_*.jsonl`, `satuned_K*.jsonl` | `analyse_hybrid.py` → `results/robust.json`; `analyse_ops.py` → `results/ops.json` |
| Table 14 (layout sweep, P11) | `bench_jitter.py` | `results/jitter.jsonl` | `docx_work/build_revised.py` |
| Section 5.10, Table 15 (full protocol on the irregular layout, jitter 0.3, 30 instances; fresh pools sized so that m >= 2N, decision rule in `PREREG_sigma03_ext.md`) | `bench_sigma_ext.py 0-15` and `15-30` (two processes), then `bench_sigma_ext2.py 0-15` and `15-30` (pool extension for 100 s, plan amendment 3), `bench_sigma.py` (reference and bounds) | `results/sigma03_ext.jsonl`, `results/sigma03_ext2.jsonl`, `results/sigma03.jsonl` | run `sdp_round_analysis.py` first, then `analyse_sigma_ext.py` → `results/sigma03.json`, `results/sigma03_sdpref.json`, `results/sigma03_ext.json`, and the semidefinite-rounding block of `results/sdp_round.json`; `EXCL_EARLY=1 analyse_sigma_ext.py` → `results/sigma03_ext_excl_early.json` (post hoc timing check) (the first-round analysis of `bench_sigma.py` alone is `analyse_sigma.py` → `results/sigma03_v1.json`) |
| Section 5.11, Table 16 (move-set and encoding controls) | `bench_formulation.py`, `bench_formulation2.py` (one-hot-preserving move, rho = 0.06) | `results/formulation.jsonl`, `results/formulation2.jsonl` | `docx_work/r6_edits.py` |
| Section 5.11, Table B6 (bSB on unconstrained Max-Cut: G1, G22, G32 and generated G-set families) | `gset_check.py` (reads `gset/G<k>`; `--gset-only` for the files alone) | `results/gset_check.jsonl` (generated-instance rows alone: `gset_check_generated.jsonl`) | `docx_work/r6_edits.py` |
| Sections 5.8 and 5.10, Table S3 (semidefinite rounding with native local search) | `sdp_round.py 20,50,100`, `sdp_round.py 200 0,1,2,3,4 1e-4`, `SIGMA=0.3 sdp_round.py 200 0,...,29 1e-4` | `results/sdp_round.jsonl`, `results/sdp_round_sigma0.3.jsonl` | `sdp_round_analysis.py` → `results/sdp_round.json`; `SDPREF=1 analyse_sigma.py` → `results/sigma03_sdpref.json` (tests with the rounding added to the reference) |
| Section 5.7, Supplementary Section S2 (cause of the SB package collapse) | `pkg_diag.py` (energy check, coupling scales, scale swap), `pkg_ablate.py` (16-way ablation of the package's integrator) | `results/pkg_diag.json`, `results/pkg_ablate.json` | `docx_work/r6_edits.py`, `docx_work/r8_edits.py` |
| Section 5.6 (Table 9 on seeds 0-4; pool sizes of the estimator) | — | — | `table9_subset.py` → `results/table9_subset.json` |
| Section 5.11, Table 16 (tuned single-bit-flip SA, paired with bSB at matched time) | `bench_round8.py T` | `results/qubo_tuned.jsonl` | `docx_work/r6_edits.py` |
| Section 5.11 (penalty versus interference forces inside our integrator) | `force_ratio.py` | `results/force_ratio.json` | `docx_work/r8_edits.py` |
| Section 5.11, Table B8 (planted colourings: our bSB, the SB package, native SA) | `bench_round8.py P`, `bench_round8q.py` | `results/planted.jsonl`, `results/planted_small.jsonl` | `docx_work/r6_edits.py` |
| Section 5.11 (sign of the field: 15 small instances solved by enumeration, field as in Eq. (19) and flipped) | `field_sign_check.py` | `results/field_sign_check.json` | `docx_work/r9_edits.py` |
| Sections 5.8, 5.11 and 6.2 (round-3 probes: A bSB started from the semidefinite rounding; B heated bSB, field scale, replica count; C order of the repair; D replica diversity on G22) | `bench_round3.py A`, `B`, `C`, `D` | `results/round3_A.json` … `round3_D.json` | `summarise_round3.py` → `results/round3_summary.json`; `docx_work/r9_edits.py` |
| Figures 1 and 2 (matched wall-clock, P4; same computation as Table 7) | — | `paper_run/*.jsonl` | in `paper_run/`: `plot_matched_ieee.py 3`, `plot_matched_ieee.py 4` → `figures/fig_K3.png`, `figures/fig_K4.png` |
| Figure 4 (full protocol on the irregular layout) | — | `results/sigma03.json` | `plot_sigma03.py` → `figures/sigma03_matched.png` |
| Section 5.11 (initial amplitude of the oscillators: replica collapse on G22, diversity on the irregular layout) | `bench_amp.py G`, `bench_amp.py L` | `results/amp_G.json`, `results/amp_L.json` | `docx_work/r9_edits.py` |
| Section 5.6 (why Tables 7 and 9 differ at 100 s on seeds 0-4) | — | `paper_run/longrun_K3.jsonl`, `paper_run/sweep_K3.jsonl`, `results/hybrid_K3.jsonl` | `table7_vs_table9.py` → `results/table7_vs_table9.json` |
| Section 5.9, daggers in Tables 13 and B7 (pool-size check of the regular-layout P6 tests: exact estimator, efforts with m >= 2N only) | — | `results/hybrid_K*.jsonl`, `hybext_K*.jsonl`, `satuned_K*.jsonl` | run `analyse_hybrid.py` first, then `p6_eligibility.py` → `results/p6_eligibility.json` |
| Section 5.9 (robustness of the paired tests: cross-fitted effort, raw-cost and bound-normalised differences) | `run_sdp200_more.py 3`, `run_sdp200_more.py 4` (bounds for seeds 5-29) | `results/sdp_rig_p6.jsonl` | `analyse_hybrid.py`, `analyse_sigma.py` |
| Table B7 (raw costs, paired differences with bootstrap intervals) | — | — | `analyse_hybrid.py` → `results/robust.json` |
| Table B5 (path-loss exponent and interference radius) | `bench_round5.py` part C | `results/eta_dmax.jsonl` | `docx_work/r5_edits.py` |

The best-known reference used for every gap is built by `refs.build_ref()` from all runs in `results/` and `paper_data/`; `build_ref(no_native_sa=True)` gives the weaker reference of Section 5.9.

## Reproducing the main statistics

From the raw results, without re-running any solver:

```
python3 analyse_hybrid.py     # Table 9, Figure 3, paired tests and Holm-adjusted p values (results/robust.json)
python3 analyse_ops.py        # work-matched comparison, CPU cost per coupling evaluation, depth (results/ops.json)
python3 decoder_floor.py      # decoder-only baselines
cd paper_run && python3 ../asp_fix.py   # Table A1 (about 3 min)
```

These were checked in a clean copy of this repository: `robust.json` and `ops.json` are reproduced exactly (bootstrap resampling is seeded), and `asp_fix.py` reproduces `asp_fix.json` to within 4e-12.

Re-running the experiments: `bench_hybrid.py 3 200,500 results/hybrid_K3.jsonl 0-30` (and `4 200` for K = 4), `bench_hybrid_ext.py` likewise, `bench_sa_tuned.py`, `bench_sqa.py 20,50,100,200,500 results/sqa_K3.jsonl`, `bench_sdp_rig.py`, `bench_jitter.py`, `bench_round4.py`, `bench_round5.py`. Protocol P6 takes about 3 h per setting on one core; one rigorous bound takes 12–20 min at 200 cells (tolerance 1e-6) and 6–8 min at 500 cells (tolerance 1e-4; the dual-corrected bound is valid at any tolerance).

## Manuscript build

`manuscript/` holds the submitted manuscript and its Supplementary Information as Word files and as Markdown, so the quoted numbers can be checked without Word.

`docx_work/` contains the scripts that compute every number quoted in the text from the raw results and assemble the manuscript. They edit the original draft, `docx_work/orig.docx`, which is included so that they run end to end:

```
cd docx_work && python3 build_revised.py && python3 build_epj.py   # needs pandoc; writes docx_work/epj/Saafan_EPJQT_manuscript.docx and Saafan_EPJQT_supplementary_information.docx
```

Run in a clean copy of this repository, this reproduces the text of `manuscript/Saafan_EPJQT_manuscript.docx` paragraph for paragraph.

`gpu_bsb.py` is a GPU variant of bSB that was not used for any reported result.
