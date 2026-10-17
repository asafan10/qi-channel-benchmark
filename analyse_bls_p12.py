"""Breakout local search (BLS) inside the paired protocol of Table 12 (K = 3) and Table S21 (K = 4) on the irregular
layout. Pools: the 18 runs of 30 s of bench_bls.py P12 (anytime records at 1, 3, 10, 30 s) and the two 100 s runs of T12/T12K4 (records at
1-100 s); a run of length e is one sample of effort e. Estimator, eligibility rule (m >= 2N on every instance), choice of the best eligible
effort, reference and tuned native SA values exactly as in analyse_sigma_ext.py (K = 3) and analyse_sigma_k4.py (K = 4); BLS does not enter
the reference. Paired two-sided Wilcoxon tests of BLS minus tuned native SA and of BLS minus bSB with refinement, Holm over all 20 tests
(post hoc family).  -> results/bls_p12_summary.json"""
import json, io, contextlib, numpy as np
from math import comb
from scipy.stats import wilcoxon

BUD = [1, 3, 10, 30, 100]


def load_ns(path, stop):
    src = open(path).read(); ns = {"__name__": "defs"}
    with contextlib.redirect_stdout(io.StringIO()): exec(src[:src.index(stop)], ns)
    return ns


NS = {3: load_ns("analyse_sigma_ext.py", 'OUT = {"seeds"'), 4: load_ns("analyse_sigma_k4.py", "res, per = table(REF)")}


def exp_min_exact(f, N):
    f = np.sort(np.asarray(f)); m = len(f)
    if N < 1 or N > m: return np.nan
    w = np.array([comb(m - k, N - 1) for k in range(1, m + 1)], dtype=float) / comb(m, N)
    return float(np.dot(w, f))


def holm(ps):
    o = sorted(range(len(ps)), key=lambda i: ps[i]); out = [0.0] * len(ps); mx = 0.0
    for k, i in enumerate(o): mx = max(mx, min(1.0, (len(ps) - k) * ps[i])); out[i] = mx
    return out


S = {"K": {}}; rows = []; allt = {"SA-tuned": [], "bSB+polish": []}
for K in (3, 4):
    ns = NS[K]; REF = ns["REF"]
    _, per = ns["table"](REF) if K == 4 else ns["table"](REF)
    pools = {}
    for l in open(f"results/bls_P12_K{K}.jsonl"):
        r = json.loads(l); pools[r["seed"]] = list(r["runs"])        # last record per instance (seeds 20-24 at K = 4 were rerun under two-process load)
    for l in open("results/bls_T12.jsonl" if K == 3 else "results/bls_T12_K4.jsonl"):
        r = json.loads(l)
        if r["seed"] in pools and len(pools[r["seed"]]) == 18: pools[r["seed"]].extend(r["runs"])
    seeds = [s for s in (ns["SEEDS"] if K == 3 else ns["SEEDS"]) if s in pools and len(pools[s]) == 20]
    full = len(seeds) == len(ns["SEEDS"])
    idx = [list(ns["SEEDS"]).index(s) for s in seeds]
    res = {}; perb = {}
    for T in BUD:
        cand = {}
        for e in BUD:
            vals = []; ok = True
            for s in seeds:
                sm = [r[str(e)] for r in pools[s] if str(e) in r]; N = int(T // e)
                if N < 1 or len(sm) < 2 * N: ok = False; break
                vals.append((exp_min_exact(sm, N) / REF[s] - 1) * 100)
            if ok: cand[e] = np.array(vals)
        e_ = min(cand, key=lambda e: np.median(cand[e])); perb[T] = cand[e_]
        rng = np.random.default_rng(T); bs = [np.median(rng.choice(cand[e_], len(seeds))) for _ in range(2000)]
        res[str(T)] = dict(med=float(np.median(cand[e_])), ci=np.percentile(bs, [2.5, 97.5]).tolist(), effort=e_,
                           below_ref=int((cand[e_] < 0).sum()))
    tst = {}
    for other in ("SA-tuned", "bSB+polish"):
        tst[other] = []
        for T in BUD:
            if (other, T) not in per: continue
            dl = perb[T] - np.asarray(per[(other, T)])[idx]
            tst[other].append(dict(T=T, med_diff=float(np.median(dl)), wins_bls=int((dl < 0).sum()), n=len(dl),
                                   p=float(wilcoxon(dl).pvalue) if np.any(dl != 0) else 1.0, other_med=float(np.median(np.asarray(per[(other, T)])[idx]))))
            allt[other].append(tst[other][-1] | {"K": K})
    S["K"][str(K)] = dict(n=len(seeds), complete=full, res=res, tests=tst)
_all = allt["SA-tuned"] + allt["bSB+polish"]
for t, h in zip(_all, holm([t["p"] for t in _all])): t["p_holm"] = h
S["holm"] = allt
f1 = lambda x: f"{x:.1f}".replace("-", "−"); f2 = lambda x: f"{x:+.2f}".replace("-", "−")
pf = lambda p: "< 0.001" if p < 0.001 else f"{p:.3f}"
for K in (3, 4):
    for T in BUD:
        r = S["K"][str(K)]["res"][str(T)]
        ts = next((t for t in allt["SA-tuned"] if t["K"] == K and t["T"] == T), None)
        tp = next((t for t in allt["bSB+polish"] if t["K"] == K and t["T"] == T), None)
        rows.append([str(K), f"{T} s", f"{f1(r['med'])} [{f1(r['ci'][0])}, {f1(r['ci'][1])}]", f"{r['effort']} s",
                     f"{f1(ts['other_med'])}" if ts else "–", f"{f2(ts['med_diff'])}; {pf(ts['p_holm'])}" if ts else "–",
                     f"{f2(tp['med_diff'])}; {pf(tp['p_holm'])}" if tp else "–"])
S["rows"] = rows
S["header"] = ["$K$", "Budget", "BLS", "BLS run length", "SA-native (tuned)", "BLS minus tuned SA; p", "BLS minus bSB + polish; p"]
S["widths"] = [500, 800, 1600, 1100, 1300, 1700, 1700]
json.dump(S, open("results/bls_p12_summary.json", "w"), indent=1, ensure_ascii=False)
for r in rows: print(r)
print({K: (S["K"][K]["n"], S["K"][K]["complete"]) for K in S["K"]})
