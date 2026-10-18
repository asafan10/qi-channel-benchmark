"""Tables S10, S11, S13, S14, S15, S16 and S17 of Additional file 1, and the numbers of Sections 5.3, 5.5 and 5.7 that come from the
same records, computed from the raw results. Run from the repository root:

    python3 paper_tables.py            -> paper_tables/Table_S*.csv and results/section_numbers.json

Column headings and symbols are those of the paper (K channels, n cells, rho penalty scale, phi raw feasibility, eta path-loss
exponent, dₘₐₓ interference radius, zeta layout jitter, lambda one-hot penalty weight). Gaps are in per cent of the reference
(best-known) cost of refs.py unless stated otherwise. The other tables of the paper are written by the analysis scripts listed in
README.md; paper_tables/ also holds every numerical table of the paper as published, for comparison."""
import csv, glob, json, os
import numpy as np
from scipy.stats import spearmanr
from refs import build_ref

REF = build_ref()
OUT = "paper_tables"; os.makedirs(OUT, exist_ok=True)
MINUS = lambda s: s.replace("-", "−")
rows_of = lambda f: [json.loads(l) for l in open(f)]


def write(name, header, rows):
    with open(os.path.join(OUT, f"Table_{name}.csv"), "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerows([header] + rows)
    print(f"Table {name}: {len(rows)} rows")


# ---------------------------------------------------------------- Table S10: path-loss exponent and interference radius
ED = {}
for r in rows_of("results/eta_dmax.jsonl"): ED.setdefault((r["eta"], r["dmax"]), []).append(r)
rows = []
for (eta, dm), rs in sorted(ED.items()):
    best = lambda r: min(r["reuse"], r["reuse_sa"], *r["sa"])
    g = lambda v: float(np.mean([(v(r) / best(r) - 1) * 100 for r in rs]))
    rows.append([f"{eta:g}", f"{dm:g}", f"{g(lambda r: r['reuse']):.1f}", f"{g(lambda r: min(r['sa'])):.1f}", f"{g(lambda r: r['reuse_sa']):.1f}",
                 f"{sum(r['reuse'] < min(r['sa']) for r in rs)}/{len(rs)}"])
write("S10", ["η", "dₘₐₓ", "Reuse", "SA-native", "Reuse + SA", "Reuse better"], rows)

# ---------------------------------------------------------------- Table S11: bSB on unconstrained Max-Cut
G = rows_of("results/gset_check.jsonl")
rows = []
for g in sorted(G, key=lambda g: (not g["name"].startswith("G"), g["n"], g["name"])):
    b = g["bsb"]
    rows.append([g["name"], f"{g['n']:,}", f"{g['ref']:,.0f}"] + [f"{b[s]['gap_best']:.2f} / {b[s]['gap_median']:.2f}" for s in ("1000", "3000", "10000", "30000")])
write("S11", ["Instance", "n", "Best cut", "1,000", "3,000", "10,000", "30,000"], rows)

# ---------------------------------------------------------------- Table S13: planted 3-colourings (optimum 0; cost in % of the instance scale)
PQ = rows_of("results/planted_small.jsonl")
rows = []
for n in (50, 100, 200):
    s = [x for x in PQ if x["n"] == n]
    bm = lambda k, f: float(np.mean([f(x["res"][k]) / x["scale"] * 100 for x in s]))
    rows.append([str(n), f"{bm('bSB|3000', min):.2f} / {bm('bSB|3000', np.median):.2f}", f"{sum(min(x['res']['bSB|3000']) < 1e-9 for x in s)}/{len(s)}",
                 f"{bm('pkg-bSB', min):.2f} / {bm('pkg-bSB', np.median):.2f}", f"{np.mean([np.mean(x['res']['pkg-phi']) for x in s]):.2f}",
                 f"{abs(bm('SA-native|64000', min)):.2f} / {bm('SA-native|64000', np.median):.2f}", f"{sum(min(x['res']['SA-native|64000']) < 1e-9 for x in s)}/{len(s)}"])
write("S13", ["n", "bSB (ours)", "Optimum", "SB package", "Package ϕ", "Native SA", "Optimum"], rows)

# ---------------------------------------------------------------- Table S14: one-hot penalty sweep (protocol P1, K = 3)
PEN = {}
for f, n in (("results/p1_100.jsonl", 100), ("results/pen200.jsonl", 200)):
    for r in rows_of(f):
        for k, v in r["res"].items(): PEN.setdefault((n, k), []).append(((v["E"] / REF[(3, n, r["seed"])] - 1) * 100, v["restarts"], v["feas"]))
pm = lambda n, k, i: float(np.mean([x[i] for x in PEN[(n, k)]]))
rows = []
for n in (100, 200):
    for rho in ("1.0", "0.5", "0.25", "0.12", "0.06"):
        rows.append([str(n) if rho == "1.0" else "", f"{float(rho):g}", f"{pm(n, 'SA-QUBO|' + rho, 0):.1f} [{pm(n, 'SA-QUBO|' + rho, 1):.0f}]",
                     f"{pm(n, 'bSB|' + rho, 0):.1f} [{pm(n, 'bSB|' + rho, 1):.0f}]", f"{pm(n, 'SA-QUBO|' + rho, 2):.2f}", f"{pm(n, 'bSB|' + rho, 2):.2f}"])
write("S14", ["n", "ρ", "SA-QUBO gap", "bSB gap", "SA-QUBO ϕ", "bSB ϕ"], rows)

# ---------------------------------------------------------------- Table S15: protocol P1 with simulated quantum annealing (K = 3)
SQA = rows_of("results/sqa_K3.jsonl")
rows = []
for n in (20, 50, 100, 200, 500):
    sel = [r for r in SQA if r["n"] == n]; cells = []
    for s in ("SQA", "SA-native", "SA-QUBO", "bSB"):
        g = [(r["res"][s]["E"] / REF[(3, n, r["seed"])] - 1) * 100 for r in sel]
        cells.append(f"{np.mean(g):.1f} ± {np.std(g):.1f} [{np.mean([r['res'][s]['restarts'] for r in sel]):.0f}]")
    rows.append([str(n), f"{sel[0]['budget']:.0f} s"] + cells)
write("S15", ["n", "Budget", "SQA", "SA-native", "SA-QUBO", "bSB"], rows)

# ---------------------------------------------------------------- Table S16: layout sweep (protocol P11, 200 cells)
JIT = {}
for r in rows_of("results/jitter.jsonl"):
    best = min(min(v) for v in r["res"].values())
    for k, v in r["res"].items(): JIT.setdefault((r["layout"], r["K"], k), []).append((min(v) / best - 1) * 100)
rows = []
for lay, lab in (("0.15", "0.15"), ("0.3", "0.3"), ("0.5", "0.5"), ("poisson", "Poisson")):
    for K in (3, 4):
        rows.append([lab, str(K)] + [f"{np.mean(JIT[(lay, K, k)]):.1f}" for k in ("reuse", "SA x256", "reuse+SA", "bSB+polish", "bSB from reuse+polish")])
write("S16", ["ζ", "K", "Reuse", "SA, 256,000 sweeps", "Reuse + SA", "bSB + polish", "bSB from reuse + polish"], rows)

# ---------------------------------------------------------------- Table S17: move-set and encoding controls (K = 3, 200 cells)
F1 = rows_of("results/formulation.jsonl"); F2 = rows_of("results/formulation2.jsonl")
fg = lambda F, nm: float(np.median([np.mean([(x["f"] / r["ref"] - 1) * 100 for x in r["res"][nm]]) for r in F]))
QT = rows_of("results/qubo_tuned.jsonl")
SCH = {"own": "0.5λ to 0.15", "0.5lam-0.02": "0.5λ to 0.02", "0.5-0.02": "0.5 to 0.02", "1.0-0.02": "1.0 to 0.02", "2.0-0.02": "2.0 to 0.02"}
rows = []
for lab, key, F in (("Native SA", "native", F1), ("Native SA, empty cells allowed at penalty λ (ρ=0.25)", "native+empty", F1),
                    ("One-hot QUBO, one-hot-preserving move, ρ=0.25", "qubo-swap-rho0.25", F2), ("One-hot QUBO, one-hot-preserving move, ρ=0.06", "qubo-swap-rho0.06", F2),
                    ("One-hot QUBO, single-bit flips, ρ=0.25", "qubo-sameT-rho0.25", F1), ("One-hot QUBO, single-bit flips, ρ=0.12", "qubo-sameT-rho0.12", F1),
                    ("One-hot QUBO, single-bit flips, ρ=0.06", "qubo-sameT-rho0.06", F2), ("One-hot QUBO, single-bit flips, own schedule, ρ=0.25", "qubo-orig", F1)):
    rows.append([lab] + [f"{fg(F, f'{key}|{S}'):.1f}" for S in (4000, 16000, 64000)])
for rho in (0.25, 0.12, 0.06):
    R = [x for x in QT if x["rho"] == rho]
    sa_mean = float(np.median([(np.mean(x["sa"]) / REF[(3, 200, x["seed"])] - 1) * 100 for x in R]))
    rows.append([f"One-hot QUBO, single-bit flips, tuned schedule ({SCH[R[0]['sched']]}), ρ={rho}", "–", "–", f"{sa_mean:.1f}"])
write("S17", ["Annealer", "4,000", "16,000", "64,000"], rows)

# ---------------------------------------------------------------- numbers of Sections 5.3, 5.5 and 5.7
SN = {}
CV = rows_of("results/bsb_cv.jsonl"); cfgs = list(CV[0]["res"]); tops = [c for c in cfgs if c != "default"]
cvg = {c: float(np.mean([(np.mean(r["res"][c]) / r["sa_best"] - 1) * 100 for r in CV])) for c in cfgs}
rk = sorted(tops, key=lambda c: CV[0]["tuning_score"][c])
SN["5.3_cross_validation"] = dict(instances=len(CV), top5_gap_range=[min(cvg[c] for c in tops), max(cvg[c] for c in tops)], default_gap=cvg["default"],
                                  best_tuned_beats_default=int(sum(np.mean(r["res"][rk[0]]) < np.mean(r["res"]["default"]) for r in CV)),
                                  spearman_tuning_vs_validation=float(spearmanr([CV[0]["tuning_score"][c] for c in tops], [cvg[c] for c in tops]).correlation))
pw = [r["res"]["bSB|" + rho]["E"] < r["res"]["SA-QUBO|" + rho]["E"] for f in ("results/p1_100.jsonl", "results/pen200.jsonl") for r in rows_of(f)
      for rho in ("1.0", "0.5", "0.25", "0.12", "0.06")]
SN["5.5_penalty_sweep"] = dict(bSB_below_SA_QUBO=int(sum(pw)), comparisons=len(pw),
                               best_rho_SA_QUBO_200=min(("0.25", "0.12", "0.06"), key=lambda r: pm(200, "SA-QUBO|" + r, 0)),
                               best_rho_bSB_200=min(("0.5", "0.25", "0.12", "0.06"), key=lambda r: pm(200, "bSB|" + r, 0)))
PKG = {}
for f in glob.glob("results/sbpkg_K*.jsonl"):
    for r in rows_of(f): PKG[(r["K"], r["n"], r["seed"])] = r
def pk(K, n, var):
    sel = [r for k, r in PKG.items() if k[:2] == (K, n)]
    g = [float(np.median([(x["f"] / REF[(K, n, r["seed"])] - 1) * 100 for x in r["runs"][var]])) for r in sel]
    return dict(median_gap=float(np.median(g)), min_gap=min(g), max_gap=max(g), phi=float(np.mean([x["phi"] for r in sel for x in r["runs"][var]])))
SN["5.7_sb_package"] = {f"K{K}_n{n}_{v}": pk(K, n, v) for K, n in ((3, 200), (3, 500), (4, 200)) for v in ("pkg-bSB", "pkg-dSB-heated")}
BIS = rows_of("results/sbpkg_bisect.jsonl")
ok = [r for r in BIS if r["setting"] == "held-out" and r["n"] <= 100]
SN["5.7_sb_package"]["max_energy_excess_pct_up_to_100_cells"] = max((r["E_pkg"] - r["E_ours"]) / abs(r["E_ours"]) * 100 for r in ok)
json.dump(SN, open("results/section_numbers.json", "w"), indent=1, ensure_ascii=False)
print("results/section_numbers.json")
