# Benchmarking quantum-inspired Ising solvers for wireless channel planning

Code and raw results for A. Saafan, "Benchmarking Ballistic Simulated Bifurcation on a Penalty-QUBO Formulation of Wireless Channel Planning: Domain Baselines and Matched Protocols" (submitted to *EPJ Quantum Technology*). Table, figure and section numbers below are those of the submitted manuscript; tables and sections labelled S are in its Additional file 1.

The G-set files `gset/G1`, `gset/G22` and `gset/G32` are the standard instances in the usual edge-list format, converted from the public copies in the eco-dqn repository (github.com/tomdbar/eco-dqn, `_graphs/benchmarks`); G1 was checked edge by edge against an independent copy (github.com/Open-Finance-Lab/RLSolver). Their best-known cuts (11,624, 13,359, 1,410) are the published values.

All instances are regenerated from their seeds (`bench_ran_qubo.make_instance`, `bench_jitter.make_layout`), so no instance files are needed. Every raw run is in `results/` and `paper_data/`; the analysis scripts read only those files. Script names that contain `round` refer to the order in which the experiments were run, not to sections of the paper.

## Environment

- Machine A (protocols P1–P3 and P5): Intel Core i7-5500U at 2.40 GHz, 16 GB RAM, 64-bit Windows, Python 3.14.
- Machine B (protocols P4 and P6–P11 and the checks of Sections 5.5–5.11): Intel Xeon at 2.10 GHz, two cores, 8 GB RAM, 64-bit Linux, Python 3.13, NumPy 2.5, SciPy 1.18, Numba 0.68, OR-Tools 9.15, PyTorch 2.14 (CPU), simulated-bifurcation 2.0.0, CVXPY 1.9 with SCS. One single-threaded process per core; timed runs were made two at a time so that both cores carried the same load.
- `pip install -r requirements.txt`. Import OR-Tools before CVXPY (a HiGHS library clash otherwise; the scripts already do this).

## Where each result comes from

| Result in the manuscript | Run | Raw records | Analysis |
|---|---|---|---|
| Table 2 (instance statistics) | `bench_ran_qubo.py` (`make_instance`) | — | `analyse_round12.py` (instance statistics) |
| Table 3 (solver settings) | `bench_ran_qubo.py` (`CFG`) | `paper_data/frozen_cfg.json` | — |
| Tables S1 and S2, Section S1 (faults of the earlier implementation) | `make_fixed_variants.py` (writes `qrc_fix1.py`, `qrc_fix2.py`), `asp_fix.py` (run inside `paper_run/`), `de_check.py`, `da_swap_check.py`; audited code: `quantum_ran_complete.py`, `quantum_inspired_ran.py`, `quantum_optimisation.py` | `paper_run/asp_fix.json` | `paper_run/compile_tables.py` |
| Table 5 and Tables S6, S7, S9 (protocols P1–P3 and P5) | `bench_ran_qubo.py` (P1), `bench_tts.py` (P2), `bench_sweep.py` (P3), tuning grid (P5) | `paper_data/*.jsonl`, `paper_data/tuning_grid.json` | `update_paper_ref.py`, then in `paper_run/`: `final_tables.py`, `analyse_tuned.py`, `analyse_tuned_ieee.py`, `compile_tables.py` |
| Tables 6, 7 and S8, Figures 1 and 2 (protocol P4, fresh pools on machine B) | `bench_fresh67.py` (two processes), `topup_fresh67.py` | `results/fresh67_*.jsonl` | `analyse_fresh67.py` → `paper_run/figures/fig_K3.*`, `fig_K4.*` |
| Table 7, values in parentheses (reference without native SA) | — | fresh P4 pools | `nosa_table7.py` → `results/nosa_table7.json` |
| Table S14 and Section 5.5 (one-hot penalty sweep, 100 and 200 cells) | `bench_round4.py` (100 cells), `bench_round5.py` part A (200 cells) | `results/p1_100.jsonl`, `results/pen200.jsonl` | `docx_work/build_revised.py` |
| Section 5.3 (cross-validated bSB tuning) | `bench_round5.py` part B | `results/bsb_cv.jsonl` | `docx_work/r5_edits.py` |
| Table 8 and Figure 3 (refinement after decoding, P6) | `pilot_polish.py`, `bench_hybrid.py`, `bench_hybrid_ext.py`, `bench_sa_tuned.py` | `results/hybrid_K*.jsonl`, `hybext_K*.jsonl`, `satuned_K*.jsonl` | `analyse_hybrid.py` → `figures/hybrid_matched_wallclock.*` |
| Table 8, second line of each entry (reference without native SA) | — | P6 pools | `table8_nosa.py` → `results/table8_nosa.json` |
| Section 5.6 (Table 8 restricted to seeds 0–4) | — | — | `table9_subset.py` → `results/table9_subset.json` |
| Table S15 (simulated quantum annealing under P1, protocol P8) | `sqa.py`, `test_sqa.py` (brute-force check), `sqa_pilot.py`, `bench_sqa.py` | `results/sqa_K3.jsonl` | `docx_work/build_revised.py` (`sqa_row`) |
| Section 5.7 and Section S2 (independent SB package, P8) | `sbpkg_pilot.py`, `bench_sbpkg.py`, `bench_core.py` (diagnosis), `bench_round4.py` (bisection), `pkg_diag.py`, `pkg_ablate.py` (16-way ablation of the package's integrator) | `results/sbpkg_*.jsonl`, `core_sbdiag.jsonl`, `pkg_diag.json`, `pkg_ablate.json` | `analyse_extra.py`, `docx_work/r6_edits.py`, `docx_work/r8_edits.py` |
| Section 5.7 (minimal reproduction of the package collapse) | `pkg_collapse_repro.py` | — | — |
| Table 9 (dual-corrected semidefinite bounds, P9) | `bound_probe.py` (CP-SAT bound), `bench_sdp.py` (SCS), `bench_sdp_rig.py` (dual-corrected bound, tolerance 1e-6 up to 200 cells), `run_sdp500b.py` (500 cells, tolerance 1e-4) | `results/sdp_*.jsonl`, `results/sdp_rig.jsonl` | `refs.py` (`sdp_bounds`) |
| Table 10 (decoder-only and domain baselines, P7, P10) | `decoder_floor.py`, `reuse_baseline.py`, `bench_reuse_polish.py` | `results/decoder_floor.json`, `reuse.json`, `reuse_polish.jsonl` | `refs.py` (`reuse_data`) |
| Table S20 (fixed reuse pattern against the per-instance best pattern) | `fixed_reuse.py` | — | `results/fixed_reuse.json` |
| Table S3, Sections 5.8, 5.10 and S4 (semidefinite rounding with native local search) | `sdp_round.py 20,50,100`, `sdp_round.py 200 0,1,2,3,4 1e-4`, `SIGMA=0.3 sdp_round.py 200 0,...,29 1e-4` | `results/sdp_round.jsonl`, `results/sdp_round_sigma0.3.jsonl` | `sdp_round_analysis.py` → `results/sdp_round.json`; `SDPREF=1 analyse_sigma.py` → `results/sigma03_sdpref.json` |
| Table 11, Table S12, Sections 5.9 and S3 (tuned native SA, paired Wilcoxon tests with Holm correction, work and depth) | `bench_core.py` (tuning of the starting temperature), `bench_sa_tuned.py` | `results/core_*.jsonl`, `satuned_K*.jsonl` | `analyse_hybrid.py` → `results/robust.json`; `analyse_ops.py` → `results/ops.json` |
| Section S19 (daggers of Tables 11 and S12: pool-size check of the regular-layout tests) | — | P6 pools | after `analyse_hybrid.py`: `p6_eligibility.py` → `results/p6_eligibility.json` |
| Section S19 (cross-fitted effort, raw-cost and bound-normalised differences) | `run_sdp200_more.py 3`, `run_sdp200_more.py 4` (bounds for seeds 5–29) | `results/sdp_rig_p6.jsonl` | `analyse_hybrid.py`, `analyse_sigma.py` |
| Section S19 (tests against the reference without native SA and on raw costs, both layouts) | — | as Tables 11 and 12 | after `analyse_hybrid.py`: `nosa_ref_check.py` → `results/nosa_ref_check.json` |
| Table S16 and Section 5.10 (layout sweep, P11) | `bench_jitter.py` | `results/jitter.jsonl` | `docx_work/build_revised.py` |
| Table S10 (path-loss exponent and interference radius) | `bench_round5.py` part C | `results/eta_dmax.jsonl` | `docx_work/r5_edits.py` |
| Table 12, Figure 4, Section S5 and Table S4 (full protocol on the irregular layout, jitter 0.3, 30 instances; analysis plan `PREREG_sigma03_ext.md`) | `bench_sigma_ext.py 0-15` and `15-30` (two processes), then `bench_sigma_ext2.py 0-15` and `15-30` (pool extension for 100 s, plan amendment 3), `bench_sigma.py` (reference and bounds) | `results/sigma03_ext.jsonl`, `results/sigma03_ext2.jsonl`, `results/sigma03.jsonl` | `sdp_round_analysis.py` first, then `analyse_sigma_ext.py` → `results/sigma03.json`, `results/sigma03_ext.json`; `plot_sigma03.py` → `figures/sigma03_matched.*`; first analysis (Table S4): `analyse_sigma.py` → `results/sigma03_v1.json`; `EXCL_EARLY=1 analyse_sigma_ext.py` → timing check of Section S5 |
| Table S21 (the same protocol at four channels; analysis plan `PREREG_irregular_K4.md`, committed before the runs) | `bench_sigma_k4.py 0-15` and `15-30` (two processes) | `results/sigma03_K4.jsonl` | `analyse_sigma_k4.py` → `results/sigma03_K4_summary.json` |
| Sections 5.2 and 5.10, Table S24 (breakout local search at matched time and in the paired protocol) | `bench_bls.py T7`, `T12`, `T12K4`, `P12 a-b K` (18 runs of 30 s per instance; seeds 20–24 at K = 4 were rerun with `redo`, the analysis keeps the last record) | `results/bls_*.jsonl` | `analyse_bls.py`, `analyse_bls_p12.py` → `results/bls_p12_summary.json` |
| Sections 5.10 and 5.11, Table 14 (tuned PT-QUBO on the irregular layout, the polish from six starts, bSB with a polish of 256,000 sweeps) | `bench_round12.py PTTUNE`, `PTTUNE2`, `PTTUNE3`, `PTTUNE4` (grid extensions; the earlier grids are kept in `results/r12_pttune_grid2.json` and `results/r12_pttune_grid3.json`), `PTPOOL` (pools of the final best setting, `run_round12b.sh`; the pools of the previous best in `results/r12_ptpool_grid3.jsonl`), `POLSTART`, `LONGPOL` (`run_round12.sh`) | `results/r12_*.json(l)` | `analyse_round12.py` → `results/round12_summary.json` |
| Table S17 and Section S18 (move-set controls; tuned single-bit-flip SA paired with bSB) | `bench_formulation.py`, `bench_formulation2.py`, `bench_round8.py T` | `results/formulation.jsonl`, `formulation2.jsonl`, `qubo_tuned.jsonl` | `docx_work/r6_edits.py` |
| Tables S18 and S19, Sections 5.11 and S21 (encoding controls: two channels, domain wall, per-cell penalty; phase domains) | `bench_round6.py K2`, `DW`, `PC`; `bench_round6b.py` | `results/round6_*.json`, `round6b_*.json` | `analyse_round6.py` |
| Table S11, Sections 5.11 and S11 (bSB on unconstrained Max-Cut: G1, G22, G32 and generated G-set families; replica diversity) | `gset_check.py` (reads `gset/G<k>`), `bench_amp.py G`, `bench_amp.py L` | `results/gset_check.jsonl`, `amp_G.json`, `amp_L.json` | `docx_work/r6_edits.py`, `docx_work/r9_edits.py` |
| Table S13 and Section S13 (planted colourings; penalty and interference forces; sign of the field) | `bench_round8.py P`, `bench_round8q.py`, `force_ratio.py`, `field_sign_check.py` | `results/planted.jsonl`, `planted_small.jsonl`, `force_ratio.json`, `field_sign_check.json` | `docx_work/r6_edits.py`, `r8_edits.py`, `r9_edits.py` |
| Sections 5.8, 5.11, 6.2, S17 and S24 (probes: A bSB started from the semidefinite rounding; B heated bSB, field scale, replica count; C order of the repair; D replica diversity on G22) | `bench_round3.py A`, `B`, `C`, `D` | `results/round3_A.json` … `round3_D.json` | `summarise_round3.py` → `results/round3_summary.json` |
| Sections 5.3, 6.2 and S23 (heated bSB, re-tuning at four channels, spectral coupling scale, polish from random plans) | `bench_round7.py C0`, `POL`, `K4T`, `HEAT`, `HEATP` | `results/round7_*.json` | `analyse_round7.py` |
| Sections 5.3, 5.6, 5.11, 6.2 and S23 (all replicas and best-visited states; one-hot bSB at two channels; cost of the polish; constraint-preserving bSB) | `bench_round9.py REP`, `REP4`, `OH2`, `POLT`, `CPT`, `CPTB`, `CPP`; `bench_round10.py NB`, `PL`, `CPPN` | `results/round9_*`, `round10_*` | `analyse_round9.py`, `analyse_round10.py` |
| Table S22 and Section 5.7 (simulated quantum annealing tuned on held-out instances) | `bench_sqa_tuned.py TUNE`, `TUNE2` (grid extension), `POOL 0-5 3`, `POOL 0-5 4` | `results/sqa_tune.json` (first grid: `sqa_tune_grid1.json`), `sqa_pool_K*.jsonl` | `analyse_sqa_tuned.py` → `results/sqa_tuned_summary.json` |
| Tables S23 and S25, Sections 5.11 and 6.1 (single-bit-flip SA on the QUBO and bSB with tuned penalty scales) | `bench_saq_tuned.py TUNE`, `POOL 0-5 3`, `POOL 0-5 4` (`run_round6.sh`); `bench_lowpen.py SAQ2`, `BSBT`, `BSBP 0-5 K rho` (`run_round7.sh`) | `results/saq_tune.json` (first grid: `saq_tune_grid1.json`), `saq_pool_K*.jsonl`, `bsb_rho_tune.json`, `bsb_rho*_K*.jsonl` | `analyse_saq_tuned.py`, `analyse_lowpen.py` (a) |
| Tables S26 and S27, Section 5.11 (single-bit annealers and re-tuned bSB on the irregular layout, settings transferred or tuned on seeds 100–101) | `bench_lowpen.py IRR a-b sqa`, `IRR a-b saq`; `bench_irrtune.py SAQ`, `SQA`, `BSB`, `SA`, `BSBP5`, then `POOL a-b sqa`, `saq`, `bsbp5` (`run_round8.sh`, `run_round9.sh`; bit-flip SA seeds 28–29 rerun with `redo`, the analysis keeps the last record) | `results/irr_*_K3.jsonl`, `irrtune_*.json`, `irrpool_*_K3.jsonl` | `analyse_lowpen.py` (b, c) |
| Table 6 (bSB at 256 times its baseline steps) and Section 5.7 (ten repeats of protocol P1 for native SA; basis of the 2-point reading threshold) | `bench_round11.py BSB256`, `SAREP` (`run_round11.sh`) | `results/r11_bsb256*.jsonl`, `r11_sarep.jsonl` | `analyse_round11.py` → `results/round11_summary.json`; in `paper_run/`: `two_point_rule.py` |

The best-known reference used for every gap is built by `refs.build_ref()` from all runs in `results/` and `paper_data/`; `build_ref(no_native_sa=True)` gives the weaker reference of Section 5.9.

Figures 1–4 share one style file, `figstyle.py` (one colour, marker, line style and name per method; axes in the symbols of the paper).

## Reproducing the main statistics

From the raw results, without re-running any solver:

```
python3 analyse_hybrid.py     # Table 8, Figure 3, paired tests and Holm-adjusted p values (results/robust.json)
python3 analyse_fresh67.py    # Tables 6 and 7, Figures 1 and 2
python3 plot_sigma03.py       # Figure 4
python3 analyse_ops.py        # work-matched comparison, CPU cost per coupling evaluation, depth (results/ops.json)
python3 decoder_floor.py      # decoder-only baselines
cd paper_run && python3 ../asp_fix.py   # Tables S1-S2 (about 3 min)
```

These were checked in a clean copy of this repository: `robust.json` and `ops.json` are reproduced exactly (bootstrap resampling is seeded), the four figures are reproduced pixel for pixel, and `asp_fix.py` reproduces `asp_fix.json` to within 4e-12.

Re-running the experiments: `bench_hybrid.py 3 200,500 results/hybrid_K3.jsonl 0-30` (and `4 200` for K = 4), `bench_hybrid_ext.py` likewise, `bench_sa_tuned.py`, `bench_sqa.py 20,50,100,200,500 results/sqa_K3.jsonl`, `bench_sdp_rig.py`, `bench_jitter.py`, `bench_round4.py`, `bench_round5.py`. Protocol P6 takes about 3 h per setting on one core; one dual-corrected bound takes 12–20 min at 200 cells (tolerance 1e-6) and 6–8 min at 500 cells (tolerance 1e-4; the bound is valid at any tolerance). In protocol P4 the number of runs per pool depends on measured run times, so a rerun reproduces Tables 6, 7 and S8 statistically, not bit for bit.

## Manuscript build

`manuscript/` holds the submitted manuscript and its Additional file 1 (Sections S1–S24) as Word files and as plain text, so that every quoted number can be checked without Word.

`docx_work/` contains the scripts that compute every number quoted in the text from the raw results and assemble the manuscript. They edit the original draft, `docx_work/orig.docx`, which is included so that they run end to end:

```
python3 analyse_round11.py && python3 analyse_round12.py && python3 table8_nosa.py && python3 combine_r43.py && python3 combine_r44.py
cd docx_work && python3 build_revised.py && python3 build_epj.py   # writes docx_work/epj/Saafan_EPJQT_manuscript.docx and Saafan_EPJQT_supplementary_information.docx
EPJ_MARK=1 EPJ_OUT=marked.docx EPJ_OUT_S=marked_S.docx python3 build_epj.py   # the same with the last edit pass marked (bold text, highlighted symbols)
```

The text of the results is written by `build_epj.py` from `docx_work/r39_results.json` (produced by `analyse_round9.py`, `analyse_round10.py` and `analyse_bls.py`, then `combine_r43.py` and `combine_r44.py`) and by the edit passes `docx_work/r46_edits.py` to `r55_edits.py`, run in that order; `r55_edits.py` is the final pass on notation (one symbol per quantity, bold vectors and matrices, upright function names), tables and wording. Run in a clean copy of this repository, the build reproduces `manuscript/Saafan_EPJQT_manuscript.docx` and the supplement exactly.

`gpu_bsb.py` is a GPU variant of bSB that was not used for any reported result.
