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
for which, name in (("sqa", "SQA"), ("saq", "SA-QUBO")):
    rr = {}
    for l in open(f"results/irr_{which}_K3.jsonl"): r = json.loads(l); rr[r["seed"]] = r
    seeds = [s for s in SEEDS if s in rr]; idx = [SEEDS.index(s) for s in seeds]
    P_ = {s: {e: [(x["f"], x["t"]) for x in v] for e, v in rr[s]["runs"].items()} for s in seeds}
    res[name] = dict(n=len(seeds), phi_mean=float(np.mean([x["phi"] for s in seeds for v in rr[s]["runs"].values() for x in v])), T={})
    for T in BUD:
        b = best_eff(P_, seeds, REF12, T)
        if b is None or ("SA-tuned", T) not in per: continue
        sa = np.asarray(per[("SA-tuned", T)])[idx]; dl = b[1] - sa
        tests.append(dict(solver=name, T=T, med=float(np.median(b[1])), sa=float(np.median(sa)), effort=b[0], med_diff=float(np.median(dl)),
                          wins=int((dl < 0).sum()), losses=int((dl > 0).sum()), n=len(dl), p=float(wilcoxon(dl).pvalue) if np.any(dl != 0) else 1.0))
for t, h in zip(tests, holm([t["p"] for t in tests])): t["p_holm"] = h
pf = lambda p: "< 0.001" if p < 0.001 else f"{p:.3f}"
for t in tests:
    rows_b.append([t["solver"], f"{t['T']} s", f1(t["med"]), f1(t["sa"]), f"{t['med_diff']:+.2f}".replace("-", "−"), f"{t['wins']}/{t['losses']}", pf(t["p_holm"])])
S["b"] = dict(res=res, tests=tests, rows=rows_b, header=["Solver", "Budget", "Gap", "SA-native (tuned)", "Difference", "Lower/higher", "p"],
              widths=[1300, 800, 1000, 1500, 1200, 1300, 1000])
A3, A4 = S["a"]["3"], S["a"]["4"]
RB = S["b"]["res"]
json.dump(S, open("results/lowpen_summary.json", "w"), indent=1, ensure_ascii=False)
for r in rows_a: print(r)
for r in rows_b: print(r)
print({k: v for k, v in res.items()}, {K: S["a"][K]["phi"] for K in ("3", "4")})
