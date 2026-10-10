from pathlib import Path
parts=[Path("tables_part1.md").read_text().split("\n",1)[1]]

import json, collections
import numpy as np
_rows = json.load(open("asp_fix.json")); _g = collections.defaultdict(list)
for r in _rows: _g[(r["n"], r["solver"])].append(r)
_names = {"DA as-given": "DA as given", "DA + dE fix": "DA + energy-difference fix", "DA + dE fix + swap fix": "DA + energy-difference fix + swap fix", "PA as-given": "PA as given", "PA + dE fix": "PA + energy-difference fix"}
t4b = ["\n## Table 4b. Effect of two further faults found in the 'corrected' DA/PA code (K=3; n=20 and 50; 5 instances; gap to final best-known, %)\n",
"Fault 5: the single-flip energy difference `(1-2x_i)(col@x - Q_ii*x_i)` omits the diagonal term Q_ii for 0->1 flips (exact for 1->0). Verified against brute force (max error 6.4 on random QUBOs); the same expression occurs in DA, PA and the framework's classical-SA routine.",
"Fault 6 (DA only): replica-exchange acceptance exp(min(500, dB*(E[r+1]-E[r]))) has the opposite sign of the standard criterion exp(dB*(E[r]-E[r+1])). Verified on a 12-variable QUBO against exact Boltzmann averages (as coded: max error 4.03; standard: 0.017).\n",
"| n | solver | gap % | raw one-hot feasibility | time (s) |", "|---|---|---|---|---|"]
for n in (20, 50):
    for k in _names:
        v = _g[(n, k)]
        t4b.append(f"| {n} | {_names[k]} | {np.mean([r['gap'] for r in v]):.1f} | {np.mean([r['feas'] for r in v]):.2f} | {np.mean([r['t'] for r in v]):.1f} |")
t4b.append("\nThe energy-difference fault explains most of the DA/PA gap in Table 4; the swap-sign fault has no measurable effect on this metric. Remaining gaps (19-61%) reflect default settings (fixed temperature range 0.001-8, fixed sweeps), not isolated. SBA (Table 4) is unaffected by these faults. SQA was not benchmarked. A further suspected issue in PA (importance weights use the absolute inverse temperature at every step instead of the increment) was NOT tested.\n")

parts.append("\n".join(t4b))
parts.append("\n## Table 5. Independent-run TTS at the frozen baseline effort (K=3 and K=4; eps = 5% and 2%): median TTS99 in ms [95% bootstrap CI]; pooled per-run success p\n")
for K in (3,4):
    t=Path(f"tts_K{K}.md").read_text()
    for e in ("5%","2%"):
        i=t.index(f"epsilon = {e}"); j=t.find("###",i+5)
        seg=t[i:j if j>0 else None]; parts.append(f"\n### K={K}, eps = {e}\n"+"\n".join(seg.split("\n")[1:]))
parts.append("\n## Table 6. Run-length-optimised TTS (sweep grid only: SA-native up to 16x, bSB up to 2x): median TTS99 in ms [95% CI incl. re-selection of run length]; (best factor, per-run success p, per-run time)\n")
for K in (3,4):
    t=Path(f"sweep_K{K}.md").read_text(); i=t.index("### K=%d, eps=5%%"%K); j=t.index("### K=%d: median single-run"%K); parts.append(t[i:j])
parts.append(Path("tables_tuning.md").read_text())
parts.append("\n## Table 8a. Long-run single-run quality and matched wall-clock (SA-native vs untuned bSB), K=3\n"+Path("longrun_K3.md").read_text().split("\n",1)[1])
parts.append("\n## Table 8b. Same, K=4, n=200 (bSB frozen K=3 tuning NOT applied; untuned bSB only)\n"+Path("longrun_K4.md").read_text().split("\n",1)[1])
parts.append("\n## Table 9. Matched wall-clock with tuned bSB (frozen config), K=3\n"+Path("tuned_K3.md").read_text().split("\n",1)[1])
parts.append("""
## Claims retired or rescoped by these tables
- **Retired:** "bSB reaches parity with native-space SA". Parity appears only under the short-restart protocol (Table 2). With long runs SA-native keeps improving while bSB plateaus, at K=3 and K=4 (Tables 8a-8b), and also after tuning at K=3 (Tables 7, 9).
- **Retired:** the first time-to-epsilon and the baseline-effort TTS as solver rankings (Tables 5-6 supersede them).
- **Rescoped:** "Tabu is best at K=4 and small sizes". It holds under fixed-budget restarts (Table 2) and for baseline-effort TTS at 50-100 cells. At K=4, n=200, run-length-optimised SA-native reaches 5% about three times faster than tabu (576 [562, 600] ms vs 1,664 [702, 2,201] ms; intervals do not overlap); at 2% they are comparable. Tabu was not included in the long-run comparison.
- **Rescoped:** the penalty sweep (Table 3): rho = 0.25, 0.12 and 0.06 are close for bSB, not a sharp optimum at 0.25.
- **Corrected statement on Table 4:** the as-given 'corrected' DA/PA still contained two further faults (Table 4b). The energy-difference fault explains most of their gap (DA 126% to 41%, PA 97% to 19% at 20 cells); the rest reflects default settings and is not isolated. SBA is unaffected.
- **Robust:** penalty-based QUBO solvers (SA-QUBO, PT-QUBO, dSB) are far behind at every effort level.
- **Limits:** 5 instances per size; K=4 long-run only at n=200 and without bSB re-tuning; no long-run at n=100; reference is heuristic and likely still too high at 500 cells.
""")
Path("FINAL_TABLES.md").write_text("# Final tables\n\nAll gaps are scored against a single final best-known cost per instance (minimum over every run of every solver in every experiment).\n"+"\n".join(parts))
print(len(Path("FINAL_TABLES.md").read_text().splitlines()),"lines")
