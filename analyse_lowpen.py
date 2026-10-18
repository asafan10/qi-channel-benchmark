"""Penalty tuning for bSB and bit-flip SA; single-bit annealers on the irregular layout.
 (a) bSB at penalty scale 0.06 (the bit-flip SA optimum) and 0.12 (its own optimum on the held-out instances) at matched wall-clock time on
     the instances of Table 7, next to tuned bit-flip SA, tuned native SA and default bSB (rho 0.25); same estimator, eligibility rule,
     reference and repair (raw feasibility reported).
 (b) tuned SQA and tuned bit-flip SA on the 30 irregular instances of Table 12 (K = 3): exact expected best of N, eligible efforts, the
     reference of Table 12 (which they do not enter), paired two-sided Wilcoxon tests against tuned native SA, Holm over the ten tests (post hoc).
 -> results/lowpen_summary.json"""
import json, glob, io, contextlib, numpy as np
from math import comb
from scipy.stats import wilcoxon
from refs import build_ref

BUD = [1, 3, 10, 30, 100]
f1 = lambda x: "–" if x is None else f"{x:.1f}".replace("-", "−")


def exp_min_exact(f, N):
    f = np.sort(np.asarray(f)); m = len(f)
    if N < 1 or 2 * N > m: return np.nan
    w = np.array([comb(m - k, N - 1) for k in range(1, m + 1)], dtype=float) / comb(m, N)
    return float(np.dot(w, f))


def best_eff(P_, seeds, ref, T):
    best = None
    for e in P_[seeds[0]]:
        v = np.array([(exp_min_exact([f for f, _ in P_[s][e]], int(T // np.mean([t for _, t in P_[s][e]]))) / ref[s] - 1) * 100 for s in seeds])
        if np.isnan(v).any(): continue
        if best is None or np.median(v) < np.median(best[1]): best = (e, v)
    return best


def holm(ps):
    o = sorted(range(len(ps)), key=lambda i: ps[i]); out = [0.0] * len(ps); mx = 0.0
    for k, i in enumerate(o): mx = max(mx, min(1.0, (len(ps) - k) * ps[i])); out[i] = mx
    return out


S = {"a": {}, "b": {}}
# ---------------------------------------------------------------- (a)
REF = build_ref(); fr = {}
for f in glob.glob("results/fresh67_K*.jsonl"):
    for l in open(f): r = json.loads(l); fr[(r["K"], r["n"], r["seed"])] = r
for l in open("results/fresh67_topup.jsonl"):
    t = json.loads(l)
    for e, v in t["runs"].items(): fr[(t["K"], t["n"], t["seed"])]["runs"][t["kind"]][e] += v
rows_a = []
for K in (3, 4):
    seeds = list(range(5)); ref = {s: REF[(K, 200, s)] for s in seeds}; pools = {}; phi = {}
    for rho in ("0.06", "0.12"):
        rr = {json.loads(l)["seed"]: json.loads(l) for l in open(f"results/bsb_rho{rho}_K{K}.jsonl")}
        pools[f"bSB rho {rho}"] = {s: {e: [(x["f"], x["t"]) for x in v] for e, v in rr[s]["runs"].items()} for s in seeds}
        phi[rho] = dict(mean=float(np.mean([x["phi"] for s in seeds for v in rr[s]["runs"].values() for x in v])),
                        min=float(min(x["phi"] for s in seeds for v in rr[s]["runs"].values() for x in v)))
    sq = {json.loads(l)["seed"]: json.loads(l) for l in open(f"results/saq_pool_K{K}.jsonl")}
    pools["SA-QUBO"] = {s: {e: [(x["f"], x["t"]) for x in v] for e, v in sq[s]["runs"].items()} for s in seeds}
    for sv in ("SA-tuned", "bSB"):
        pools[sv] = {s: {e: [(x["f"], x["t"]) for x in v] for e, v in fr[(K, 200, s)]["runs"][sv].items()} for s in seeds}
    med = {sv: {} for sv in pools}; wins = {}
    for sv, P_ in pools.items():
        for T in BUD:
            b = best_eff(P_, seeds, ref, T)
            if b: med[sv][T] = b[1]
    for T in BUD:
        rows_a.append([str(K), f"{T} s"] + [f1(float(np.median(med[sv][T]))) if T in med[sv] else "–"
                                            for sv in ("bSB", "bSB rho 0.12", "bSB rho 0.06", "SA-QUBO", "SA-tuned")])
    S["a"][str(K)] = dict(phi=phi, med={sv: {str(T): float(np.median(v)) for T, v in d.items()} for sv, d in med.items()},
                          wins_vs_saq={str(T): int(sum(a < b for a, b in zip(med["SA-QUBO"][T], med["bSB rho 0.06"][T])))
                                       for T in BUD if T in med["SA-QUBO"] and T in med["bSB rho 0.06"]})
S["a"]["header"] = ["$K$", "Budget", "bSB, ρ = 0.25", "bSB, ρ = 0.12", "bSB, ρ = 0.06", "SA-QUBO, ρ = 0.06", "SA-native (tuned)"]
S["a"]["rows"] = rows_a; S["a"]["widths"] = [500, 800, 1300, 1300, 1300, 1500, 1500]

# ---------------------------------------------------------------- (b)
src = open("analyse_sigma_ext.py").read(); ns = {"__name__": "defs"}
with contextlib.redirect_stdout(io.StringIO()): exec(src[:src.index('OUT = {"seeds"')], ns)
REF12 = ns["REF"]; SEEDS = list(ns["SEEDS"]); _, per = ns["table"](REF12)
tests = []; rows_b = []; res = {}
import os as _os
SRC_B = [("results/irr_sqa_K3.jsonl", "SQA", "regular"), ("results/irr_saq_K3.jsonl", "SA-QUBO", "regular"),
         ("results/irrpool_sqa_K3.jsonl", "SQA", "irregular"), ("results/irrpool_saq_K3.jsonl", "SA-QUBO", "irregular"),
         ("results/irrpool_bsbp5_K3.jsonl", "bSB", "irregular")]
for path, name, tuned in SRC_B:
    if not _os.path.exists(path): continue
    rr = {}
    for l in open(path): r = json.loads(l); rr[r["seed"]] = r
    seeds = [s for s in SEEDS if s in rr]; idx = [SEEDS.index(s) for s in seeds]
    if len(seeds) < len(SEEDS): print("incomplete:", path, len(seeds)); continue
    P_ = {s: {e: [(x["f"], x["t"]) for x in v] for e, v in rr[s]["runs"].items()} for s in seeds}
    res[f"{name}|{tuned}"] = dict(n=len(seeds), cfg=rr[seeds[0]]["cfg"], phi_mean=float(np.mean([x["phi"] for s in seeds for v in rr[s]["runs"].values() for x in v])))
    for T in BUD:
        b = best_eff(P_, seeds, REF12, T)
        if b is None or ("SA-tuned", T) not in per: continue
        sa = np.asarray(per[("SA-tuned", T)])[idx]; dl = b[1] - sa
        bs = np.asarray(per[("bSB", T)])[idx] if ("bSB", T) in per else None
        tests.append(dict(solver=name, tuned=tuned, T=T, med=float(np.median(b[1])), sa=float(np.median(sa)), effort=b[0], med_diff=float(np.median(dl)),
                          bsb=float(np.median(bs)) if bs is not None else None, wins_vs_bsb=int((b[1] < bs).sum()) if bs is not None else None,
                          wins=int((dl < 0).sum()), losses=int((dl > 0).sum()), n=len(dl), p=float(wilcoxon(dl).pvalue) if np.any(dl != 0) else 1.0,
                          p_tost=float(max(wilcoxon(dl + 0.5, alternative="greater").pvalue, wilcoxon(dl - 0.5, alternative="less").pvalue))))
for t, h in zip(tests, holm([t["p"] for t in tests])): t["p_holm"] = h
for t, h in zip(tests, holm([t["p_tost"] for t in tests])): t["p_tost_holm"] = h
S["b_bsb_alone"] = {str(T): dict(med=float(np.median(per[("bSB", T)])), sa=float(np.median(per[("SA-tuned", T)])),
                                 med_diff=float(np.median(np.asarray(per[("bSB", T)]) - np.asarray(per[("SA-tuned", T)]))))
                    for T in BUD if ("bSB", T) in per and ("SA-tuned", T) in per}
pf = lambda p: "< 0.001" if p < 0.001 else f"{p:.3f}"
for t in tests:
    rows_b.append([t["solver"], t["tuned"], f"{t['T']} s", f1(t["med"]), f1(t["sa"]), f1(t["bsb"]), f"{t['med_diff']:+.2f}".replace("-", "−"),
                   f"{t['wins']}/{t['losses']}", pf(t["p_holm"]), pf(t["p_tost_holm"])])
S["b"] = dict(res=res, tests=tests, rows=rows_b, header=["Solver", "Tuned on", "Budget", "Gap", "SA-native (tuned)", "bSB", "Difference", "Lower/higher", "p, difference", "p, equivalence"],
              widths=[1000, 950, 700, 700, 1100, 650, 950, 1000, 900, 950])
_it = {w: json.load(open(f"results/irrtune_{w}.json")) for w in ("saq", "sqa", "bsb", "sa")} if _os.path.exists("results/irrtune_sa.json") else None
# (c) the held-out tuning on the irregular layout: grid size, the setting used before (default or transferred) and the best setting
if _it and _os.path.exists("results/irrtune_bsbp5.json"):
    _it["bsbp5"] = json.load(open("results/irrtune_bsbp5.json"))
    def _g(w, cfg):
        k = json.dumps(cfg, sort_keys=True); return _it[w]["grid"][k]["mean_gap"] if k in _it[w]["grid"] else None
    _desc = lambda c: ", ".join(f"{k} {v}" for k, v in c.items())
    used = {"sa": dict(T0=0.6), "bsb": dict(rho=0.25), "bsbp5": dict(dt=0.5, gamma=1.0, c0m=1.0, rho=0.25),
            "saq": json.load(open("results/saq_tune.json"))["3"]["best"], "sqa": json.load(open("results/sqa_tune.json"))["3"]["best"]}
    names = {"sa": ("Native SA", "starting temperature"), "bsb": ("bSB", "penalty scale"), "bsbp5": ("bSB", "step size, pump exponent, coupling multiplier, penalty scale"),
             "saq": ("Bit-flip SA", "penalty scale, starting and final temperature"), "sqa": ("SQA", "Trotter slices, temperature, initial field, penalty scale")}
    rows_c = []
    for w in ("sa", "bsb", "bsbp5", "saq", "sqa"):
        gu = _g(w, used[w]); gb = _it[w]["best_gap"]
        rows_c.append([names[w][0], f"{names[w][1]} ({len(_it[w]['grid'])})", f1(gu) if gu is not None else "–", _desc(_it[w]["best"]).replace("lam/2", "λ/2"), f1(gb)])
    S["c"] = dict(rows=rows_c, header=["Solver", "Grid (configurations)", "Setting used before", "Best setting", "Best"], widths=[1100, 2900, 1200, 2700, 700])
json.dump(S, open("results/lowpen_summary.json", "w"), indent=1, ensure_ascii=False)
for r in rows_a: print(r)
for r in rows_b: print(r)
print(res)
