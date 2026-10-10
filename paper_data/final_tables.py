"""final_tables.py -- recompute all tables against ONE final reference per (K, n, seed):
   the minimum cost found by ANY run of ANY solver in ANY experiment of this project."""
import json, sys, collections
import numpy as np
from pathlib import Path

def jl(f): return [json.loads(l) for l in open(f)] if Path(f).exists() else []

# ---------------------------------------------------------------- final reference
def build_ref(K):
    ref = {}; first = {}
    def upd(k, v): ref[k] = min(ref.get(k, np.inf), v)
    res_files = {3: ["results_K3.jsonl", "results_K3_fix.jsonl"], 4: ["results_K4.jsonl"]}[K]
    merged = {}
    for f in res_files:
        for r in jl(f):
            k = (r["n"], r["seed"])
            if k in merged: merged[k]["res"].update(r["res"])
            else: merged[k] = r
    for k, r in merged.items():
        first[k] = min(v["E"] for v in r["res"].values()); upd(k, first[k])
    for f in (f"tts_K{K}.jsonl", f"sweep_K{K}.jsonl", f"refrun_K{K}.jsonl", f"longrun_K{K}.jsonl", f"tuned_K{K}.jsonl"):
        for r in jl(f):
            upd((r["n"], r["seed"]), r["new_best"])
    if K == 3:
        for r in json.load(open("lam_sweep.json")): pass  # lam sweep: n=100 seeds 0-4
        for r in json.load(open("asp.json")): upd((r["n"], r["seed"]), r["E"])
    return ref, first, merged

FINAL = {K: build_ref(K) for K in (3, 4)}
# lam_sweep costs (n=100, K=3)
for r in json.load(open("lam_sweep.json")):
    k = (100, r["seed"]); FINAL[3][0][k] = min(FINAL[3][0][k], r["E"])
json.dump({f"K{K}": {f"{n}_{s}": v for (n, s), v in FINAL[K][0].items()} for K in (3, 4)}, open("final_reference.json", "w"))

ORDER = ["Greedy", "Tabu", "SA-native", "CP-SAT", "SA-QUBO", "PT-QUBO", "bSB", "dSB"]
out = ["# Final tables (all gaps scored against the single final best-known reference per instance)\n"]

# ---------------------------------------------------------------- T1: reference drift
out.append("\n## Table 1. Reference drift: how much lower the final best-known is than the original benchmark's best-known\n")
out.append("| K | n | mean drop (%) | max drop (%) | instances improved |"); out.append("|---|---|---|---|---|")
for K in (3, 4):
    ref, first, _ = FINAL[K]
    for n in sorted({k[0] for k in first}):
        d = [100 * (first[(n, s)] - ref[(n, s)]) / first[(n, s)] for s in range(5)]
        out.append(f"| {K} | {n} | {np.mean(d):.2f} | {np.max(d):.2f} | {sum(x > 1e-9 for x in d)}/5 |")

# ---------------------------------------------------------------- T2: fixed-budget restarts
out.append("\n## Table 2. Fixed-budget (3/6/12 s) restart protocol: gap to final best-known (%), mean ± std over 5 instances\n")
for K in (3, 4):
    ref, first, merged = FINAL[K]; sizes = sorted({k[0] for k in merged})
    out.append(f"\n**K={K}**\n"); out.append("| n | " + " | ".join(ORDER) + " |"); out.append("|---|" + "---|" * len(ORDER))
    for n in sizes:
        cells = []
        for s in ORDER:
            g = [100 * (merged[(n, sd)]["res"][s]["E"] - ref[(n, sd)]) / ref[(n, sd)] for sd in range(5)]
            cells.append(f"{np.mean(g):.1f} ± {np.std(g):.1f}")
        out.append(f"| {n} | " + " | ".join(cells) + " |")

# ---------------------------------------------------------------- T3: penalty sweep
out.append("\n## Table 3. One-hot penalty sweep (K=3, n=100, 6 s per cell): gap to final best-known (%), mean over 5 instances\n")
d = json.load(open("lam_sweep.json")); ref3 = FINAL[3][0]
g = collections.defaultdict(list); fe = collections.defaultdict(list)
for r in d:
    k = (r["solver"], r["lam_scale"]); g[k].append(100 * (r["E"] - ref3[(100, r["seed"])]) / ref3[(100, r["seed"])]); fe[k].append(r["feas"])
out.append(f"Native-space SA (6 s restarts): {np.mean(g[('SA-native', None)]):.1f}%\n")
out.append("| scale rho | SA-QUBO gap % | bSB gap % | SA-QUBO raw feasibility | bSB raw feasibility |"); out.append("|---|---|---|---|---|")
for ls in [1.0, 0.5, 0.25, 0.12, 0.06]:
    out.append(f"| {ls} | {np.mean(g[('SA-QUBO', ls)]):.1f} | {np.mean(g[('bSB', ls)]):.1f} | {np.mean(fe[('SA-QUBO', ls)]):.2f} | {np.mean(fe[('bSB', ls)]):.2f} |")

# ---------------------------------------------------------------- T4: as-published solvers
out.append("\n## Table 4. Original (as-published, corrected) solvers vs re-implementations (K=3): gap to final best-known (%), 5 instances\n")
a = json.load(open("asp.json")); gg = collections.defaultdict(list)
for r in a: gg[(r["n"], r["solver"])].append((100 * (r["E"] - ref3[(r["n"], r["seed"])]) / ref3[(r["n"], r["seed"])], r["feas"], r["t"]))
out.append("| n | solver | gap % | raw one-hot feasibility | time (s) |"); out.append("|---|---|---|---|---|")
for n in (20, 50):
    for s in ["orig-SBA", "orig-DA", "orig-PA", "new-SA-native", "new-bSB", "new-SA-QUBO"]:
        v = np.array(gg[(n, s)]); out.append(f"| {n} | {s} | {v[:,0].mean():.1f} | {v[:,1].mean():.2f} | {v[:,2].mean():.1f} |")

Path("tables_part1.md").write_text("\n".join(out) + "\n"); print("\n".join(out))

# ------------------------------------------- hand over the final reference to the other analysis scripts
import analyse_tts, analyse_sweep
for K in (3, 4):
    ref = FINAL[K][0]
    rows = analyse_tts.load(f"tts_K{K}.jsonl")
    for r in rows: r["ref"] = ref[(r["n"], r["seed"])]
    analyse_tts.analyse(rows, K)                      # writes tts_K{K}.md
analyse_sweep.load_refs = lambda K: FINAL[K][0]
for K in (3, 4): analyse_sweep.analyse(K)             # writes sweep_K{K}.md
