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
| Section 5.9 (reference bias at every budget: tests against the reference without native SA and on raw costs, regular and irregular layouts) | — | as Tables 13 and 15 | run `analyse_hybrid.py` first, then `nosa_ref_check.py` → `results/nosa_ref_check.json` |
| Section 5.7 (basis of the 2-point reading rule: bootstrap half-widths of Table 7) | — | `paper_run/*.jsonl` | in `paper_run/`: `two_point_rule.py` → `results/two_point_rule.json` |
| Figures 1 and 2 (matched wall-clock, P4; same computation as Table 7) | — | `paper_run/*.jsonl` | in `paper_run/`: `plot_matched_ieee.py 3`, `plot_matched_ieee.py 4` → `figures/fig_K3.png`, `figures/fig_K4.png` |
| Figure 4 (full protocol on the irregular layout) | — | `results/sigma03.json` | `plot_sigma03.py` → `figures/sigma03_matched.png` |
| Section 5.11 (initial amplitude of the oscillators: replica collapse on G22, diversity on the irregular layout) | `bench_amp.py G`, `bench_amp.py L` | `results/amp_G.json`, `results/amp_L.json` | `docx_work/r9_edits.py` |
| Section 5.6 (why Tables 7 and 9 differ at 100 s on seeds 0-4) | — | `paper_run/longrun_K3.jsonl`, `paper_run/sweep_K3.jsonl`, `results/hybrid_K3.jsonl` | `table7_vs_table9.py` → `results/table7_vs_table9.json` |
| Section 5.9, daggers in Tables 13 and B7 (pool-size check of the regular-layout P6 tests: exact estimator, efforts with m >= 2N only) | — | `results/hybrid_K*.jsonl`, `hybext_K*.jsonl`, `satuned_K*.jsonl` | run `analyse_hybrid.py` first, then `p6_eligibility.py` → `results/p6_eligibility.json` |
| Section 5.9 (robustness of the paired tests: cross-fitted effort, raw-cost and bound-normalised differences) | `run_sdp200_more.py 3`, `run_sdp200_more.py 4` (bounds for seeds 5-29) | `results/sdp_rig_p6.jsonl` | `analyse_hybrid.py`, `analyse_sigma.py` |
| Table B7 (raw costs, paired differences with bootstrap intervals) | — | — | `analyse_hybrid.py` → `results/robust.json` |
| Table B5 (path-loss exponent and interference radius) | `bench_round5.py` part C | `results/eta_dmax.jsonl` | `docx_work/r5_edits.py` |

## Experiments added in revision (numbering of the final manuscript)

The table above uses the section and table numbers of an earlier draft. The experiments below were added in revision and use the final numbering.

| Result in the manuscript | Run | Raw records | Analysis |
|---|---|---|---|
| Tables 6, 7, S7 and Figures 1-2 (protocol P4, fresh pools on machine B) | `bench_fresh67.py` (two processes), `topup_fresh67.py` | `results/fresh67_*.jsonl` | `analyse_fresh67.py` |
| Section 5.11, encoding controls (K = 2 Max-Cut, domain wall, per-cell penalty), Tables S17-S18 | `bench_round6.py` K2, DW, PC; `bench_round6b.py` | `results/round6_*.json`, `results/round6b_*.json` | `analyse_round6.py` |
| Section 5.3, heated bSB and K = 4 re-tuning; Section 6.2, spectral coupling scale; polish from random plans | `bench_round7.py` C0, POL, K4T, HEAT, HEATP | `results/round7_*.json` | `analyse_round7.py` |
| Section 6.2, decoding all replicas and best-visited states; Section 5.11, one-hot bSB at K = 2; Section 5.6, polish sweep cost; constraint-preserving bSB tuning | `bench_round9.py` REP, REP4, OH2, POLT, CPT, CPTB, CPP | `results/round9_*` | `analyse_round9.py` |
| Section 5.3, compiled constraint-preserving bSB at matched time and on planted colourings | `bench_round10.py` NB, PL, CPPN | `results/round10_*` | `analyse_round10.py` |
| Sections 5.2 and 5.10, breakout local search at matched time | `bench_bls.py` T7, T12 | `results/bls_*.jsonl` | `analyse_bls.py` |
| Section 5.10 and Table S21, protocol P6 on the irregular layout at K = 4 (analysis plan `PREREG_irregular_K4.md`, committed before the runs); breakout local search on the same instances | `bench_sigma_k4.py 0-15` and `15-30` (two processes); `bench_bls.py T12K4 0-15` and `15-30` | `results/sigma03_K4.jsonl`, `results/bls_T12_K4.jsonl` | `analyse_sigma_k4.py` → `results/sigma03_K4_summary.json` |
| Section 5.7 and Table S22, simulated quantum annealing tuned on held-out instances, at matched wall-clock time | `bench_sqa_tuned.py TUNE`, then `TUNE2` (grid extension), then `POOL 0-5 3` and `POOL 0-5 4` (two processes) | `results/sqa_tune.json` (first grid alone: `results/sqa_tune_grid1.json`), `results/sqa_pool_K3.jsonl`, `results/sqa_pool_K4.jsonl` | `analyse_sqa_tuned.py` → `results/sqa_tuned_summary.json` |
| Sections 5.11 and 6.1 and Table S23, single-bit-flip SA on the QUBO tuned on held-out instances (penalty scale and schedule), at matched wall-clock time | `bench_saq_tuned.py TUNE`, then `POOL 0-5 3` and `POOL 0-5 4`; all round-6 jobs are queued two at a time by `run_round6.sh` | `results/saq_tune.json`, `results/saq_pool_K3.jsonl`, `results/saq_pool_K4.jsonl` | `analyse_saq_tuned.py` → `results/saq_tuned_summary.json` |
| Section 5.10 and Table S24, breakout local search in the paired protocol of the irregular layout (K = 3 and 4, 30 instances each) | `bench_bls.py P12 a-b K` (18 runs of 30 s per instance; seeds 20-24 at K = 4 were rerun with `redo` because their first runs had one process to themselves; the analysis keeps the last record) | `results/bls_P12_K3.jsonl`, `results/bls_P12_K4.jsonl` and the 100 s runs of `results/bls_T12*.jsonl` | `analyse_bls_p12.py` → `results/bls_p12_summary.json` |
| Section 5.11 and Tables S22, S25, penalty tuning extended for bit-flip SA (0.03, 0.015) and done for bSB (0.015-0.25) on the held-out instances; bSB at 0.06 and 0.12 at matched time on the instances of Table 7 | `bench_lowpen.py SAQ2`, `BSBT`, `BSBP 0-5 K rho` (queued two at a time by `run_round7.sh`) | `results/saq_tune.json` (first grid: `results/saq_tune_grid1.json`), `results/bsb_rho_tune.json`, `results/bsb_rho*_K*.jsonl` | `analyse_lowpen.py` (a) → `results/lowpen_summary.json` |
| Sections 5.7, 5.11, 6.1 and Table S26, tuned SQA and tuned bit-flip SA on the 30 irregular instances, paired tests against tuned native SA | `bench_lowpen.py IRR a-b sqa` and `IRR a-b saq` | `results/irr_sqa_K3.jsonl`, `results/irr_saq_K3.jsonl` | `analyse_lowpen.py` (b) |
| Sections 5.11 and 6.1 and Table S26, every solver tuned on two held-out instances of the irregular layout (seeds 100-101); pools on the 30 test instances for SQA and bit-flip SA, whose irregular-tuned settings differ | `bench_irrtune.py SAQ`, `SQA`, `BSB`, `SA`, then `POOL a-b sqa` and `POOL a-b saq` (`run_round8.sh`; bit-flip SA seeds 28-29 rerun with `redo` under two-process load, the analysis keeps the last record) | `results/irrtune_*.json`, `results/irrpool_*_K3.jsonl` | `analyse_lowpen.py` (b) |
| Section 5.11 and Tables S26-S27, bSB's full P5 grid (81 configurations) tuned on the irregular held-out instances, pools of the re-tuned bSB on the 30 test instances; equivalence tests (margin 0.5) for every row of Table S26 | `bench_irrtune.py BSBP5`, then `POOL a-b bsbp5` (`run_round9.sh`) | `results/irrtune_bsbp5.json`, `results/irrpool_bsbp5_K3.jsonl` | `analyse_lowpen.py` (b, c) |
| Section 5.2, Table 7 against the reference without native SA (values in parentheses in Table 7) | - | fresh P4 pools | `nosa_table7.py` → `results/nosa_table7.json` |
| Section 5.6, Table 8 against the reference without native SA (second line of each cell; seventh review, second report) | - | P6 pools | `table8_nosa.py` → `results/table8_nosa.json` |
| Section 5.2 and Table 6, bSB at 256 times its baseline steps on the five P4 instances; Section 5.7, ten repeats of protocol P1 for native SA at 100 cells (2-point threshold) | `bench_round11.py BSB256` and `SAREP` (two at a time by `run_round11.sh`) | `results/r11_bsb256*.jsonl`, `results/r11_sarep.jsonl` | `analyse_round11.py` → `results/round11_summary.json` |
| Sections 5.10 and 5.11, Table 14 (ninth to eleventh reviews): tuned PT-QUBO on the irregular layout, the polish from six starts, bSB with a polish of 256,000 sweeps | `bench_round12.py PTTUNE`, `PTTUNE2`, `PTTUNE3` (grid extension to 122 settings; the grid before it is kept in `results/r12_pttune_grid2.json`), `PTTUNE4` (extension below the best upper temperature, 158 settings; the grid before it in `results/r12_pttune_grid3.json`), `PTPOOL` (re-run for the new best setting by `run_round12b.sh`; the pools of the previous best are kept in `results/r12_ptpool_grid3.jsonl`), `POLSTART`, `LONGPOL` (two at a time by `run_round12.sh`) | `results/r12_*.json(l)` | `analyse_round12.py` → `results/round12_summary.json` |
| Section 5.10, fixed reuse pattern against the per-instance best pattern | `fixed_reuse.py` | - | `results/fixed_reuse.json` |
| Section 5.7, minimal reproduction of the SB package collapse | `pkg_collapse_repro.py` | - | - |

The text of these results is written into the manuscript by `docx_work/build_epj.py` from `docx_work/r39_results.json`, which `analyse_round9.py`, `analyse_round10.py` and `analyse_bls.py` produce in that order; `combine_r43.py` then adds the sentences generated by `analyse_sigma_k4.py` and `analyse_sqa_tuned.py`, and `combine_r44.py` those generated by `analyse_saq_tuned.py`, `analyse_bls_p12.py` and `analyse_lowpen.py` (run in that order). `docx_work/r46_edits.py`, run by `build_epj.py`, adds the values of `results/nosa_table7.json` and `results/table8_nosa.json` to Tables 7 and 8 and moves the full text of Sections 5.3, 5.5, the package part of 5.7, the move-set paragraph of 6.1 and the probes of 6.2 to Section S24, leaving summaries in the main text; `docx_work/r47_edits.py` makes the minor fixes of the third report, `docx_work/r50_edits.py` the changes for the eighth review, `docx_work/r51_edits.py` those for the ninth to eleventh (Table 14, conventions, new run results), `docx_work/r52_edits.py` those of their second reports (confirmatory families named, 100 s dropped from summary sentences), `docx_work/r53_edits.py` those of their third reports (bSB checks first in the future work, conclusions by solver), and `docx_work/r54_edits.py` the five-instance wording (no "level within the 2-point threshold") and a clickable DOI link (non-breaking hyphens) for every reference that has one; the PT-QUBO stopping rule is stated in `r51_edits.py`. Supplementary tables are numbered in order of appearance by the build.

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

`manuscript/` holds the submitted manuscript and its Supplementary Information (Sections S1–S13, including the former Appendix B) as Word files and as Markdown, so the quoted numbers can be checked without Word.

`docx_work/` contains the scripts that compute every number quoted in the text from the raw results and assemble the manuscript. They edit the original draft, `docx_work/orig.docx`, which is included so that they run end to end:

```
cd docx_work && python3 build_revised.py && python3 build_epj.py   # needs pandoc; writes docx_work/epj/Saafan_EPJQT_manuscript.docx and Saafan_EPJQT_supplementary_information.docx
```

Run in a clean copy of this repository, this reproduces the text of `manuscript/Saafan_EPJQT_manuscript.docx` paragraph for paragraph.

`gpu_bsb.py` is a GPU variant of bSB that was not used for any reported result.
