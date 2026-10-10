"""Tables for Experiments B (SDP bounds + decoder-only baselines), C (independent SB package) and D (SQA)."""
import json, glob, numpy as np
from pathlib import Path
from refs import build_ref, sdp_bounds, HAS_PAPER_REF
from exp_common import setup, decode_eval, interference, greedy

ref = build_ref(); lb = sdp_bounds(); out = [f"# Experiments B-D (reference includes paper's final_reference.json: {HAS_PAPER_REF})\n"]

# ---- D: SQA under P1, same machine
rows = [json.loads(l) for l in open("results/sqa_K3.jsonl")]
S = ["SQA", "SA-native", "SA-QUBO", "bSB"]
out += ["\n## Table D. Fixed-budget restart protocol P1 on one machine (K=3): gap to best-known (%), mean ± std over 5 instances [mean restarts]\n",
        "| n | budget | " + " | ".join(S) + " |", "|---|---|" + "---|" * len(S)]
for n in sorted({r["n"] for r in rows}):
    sel = [r for r in rows if r["n"] == n]; cells = []
    for s in S:
        g = [(r["res"][s]["E"] / ref[(3, n, r["seed"])] - 1) * 100 for r in sel]
        cells.append(f"{np.mean(g):.1f} ± {np.std(g):.1f} [{np.mean([r['res'][s]['restarts'] for r in sel]):.0f}]")
    out.append(f"| {n} | {sel[0]['budget']:.0f} s | " + " | ".join(cells) + " |")

# ---- C: independent SB package
out += ["\n## Table C. Independent SB implementation (simulated-bifurcation 2.0.0, 128 agents, held-out-tuned): single-run gap (%), median over instances [range], raw feasibility, time per run\n",
        "| K | n | variant | median gap % [range] | raw feasibility | time/run (s) | our bSB x1 single-run median gap % |", "|---|---|---|---|---|---|---|"]
hyb = {}
for f in glob.glob("results/hybrid_K*.jsonl"):
    for l in open(f):
        r = json.loads(l); hyb[(r["K"], r["n"], r["seed"])] = r
for f in sorted(glob.glob("results/sbpkg_K*.jsonl")):
    rr = [json.loads(l) for l in open(f)]
    for n in sorted({r["n"] for r in rr}):
        sel = [r for r in rr if r["n"] == n]; K = sel[0]["K"]
        ours = np.median([np.median([(x["f_sel"] / ref[(K, n, r["seed"])] - 1) * 100 for x in hyb[(K, n, r["seed"])]["runs"]["bSB"]["1"]]) for r in sel])
        for v in sel[0]["runs"]:
            g = [np.median([(x["f"] / ref[(K, n, r["seed"])] - 1) * 100 for x in r["runs"][v]]) for r in sel]
            phi = np.mean([x["phi"] for r in sel for x in r["runs"][v]]); t = np.mean([x["t"] for r in sel for x in r["runs"][v]])
            out.append(f"| {K} | {n} | {v} | {np.median(g):.1f} [{min(g):.1f}, {max(g):.1f}] | {phi:.2f} | {t:.1f} | {ours:.1f} |")

# ---- B: SDP bounds and decoder-only baselines
out += ["\n## Table B. Certified lower bounds and decoder-only baselines: mean over instances\n",
        "| K | n | best-known f_ref | SDP lower bound | certified gap of f_ref to optimum, mean [max] (%) | repair of empty vector (scan-order greedy), gap % | degree-ordered greedy, gap % |",
        "|---|---|---|---|---|---|---|"]
for K, ns in ((3, (20, 50, 100, 200, 500)), (4, (50, 100, 200))):
    for n in ns:
        keys = [(K, n, s) for s in range(5) if (K, n, s) in ref]
        if not keys: continue
        fr = [ref[k] for k in keys]; L = [lb.get(k) for k in keys]
        g0 = []; gg = []
        for (k_, n_, s) in keys:
            W, csr, a, J, lam, jj = setup(n, s, K)
            g0.append((decode_eval(np.zeros(n * K, dtype=np.int64), K, csr)[0] / ref[(K, n, s)] - 1) * 100)
            gg.append((interference(greedy(csr, n, K), *csr) / ref[(K, n, s)] - 1) * 100)
        if all(x is not None for x in L):
            cg = [(f / l - 1) * 100 for f, l in zip(fr, L)]
            bl = f"{np.mean(L):.2f}"; cgs = f"{np.mean(cg):.1f} [{max(cg):.1f}]"
        else:
            bl = cgs = "-"
        out.append(f"| {K} | {n} | {np.mean(fr):.2f} | {bl} | {cgs} | {np.mean(g0):.1f} | {np.mean(gg):.1f} |")
txt = "\n".join(out) + "\n"; Path("results/extra.md").write_text(txt); print(txt)
