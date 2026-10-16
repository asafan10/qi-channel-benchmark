"""Protocol P6 on the irregular layout at K = 4 (bench_sigma_k4.py), analysed as fixed in PREREG_irregular_K4.md: exact expected best of N
without replacement, efforts eligible only where m >= 2N on every instance, paired Wilcoxon and TOST (margin 0.5 points) with Holm over
the tested budgets. Reference per instance: best plan of any run of bench_sigma_k4.py and the reuse plan. Breakout local search
(results/bls_T12_K4.jsonl, if present) is reported against the same reference, not added to it.
-> results/sigma03_K4_summary.json (table rows)"""
import json, os, itertools, math
import numpy as np
from math import comb
from scipy.stats import wilcoxon
from exp_common import interference, csr_from_W
from bench_jitter import make_layout

MARGIN = 0.5; BUDGETS = [1, 3, 10, 30, 100]; BOOT = 2000; K = 4
EXT = {r["seed"]: r for r in map(json.loads, open("results/sigma03_K4.jsonl"))}
SEEDS = sorted(EXT)


def reuse_cost(seed):
    n = 200; W = make_layout(n, seed, K, 0.3); csr = csr_from_W(W)
    cols = int(math.ceil(math.sqrt(n * 1.15))); idx = np.arange(n); c = idx % cols; r = idx // cols; q = c - (r // 2); best = np.inf
    for a, b, a2, b2 in itertools.product(range(2), repeat=4):
        p = 2 * ((a * q + b * r) % 2) + ((a2 * q + b2 * r) % 2)
        if len(set(p.tolist())) == 4: best = min(best, float(interference(p.astype(np.int64), *csr)))
    return best


def samples(r, which):
    R = r["runs"]; d = {}
    if which in ("SA-tuned", "SA-frozen"):
        for f, recs in R[which].items(): d[f"x{f}"] = [(x["f"], x["t"]) for x in recs]
    elif which == "bSB":
        for f, recs in R["bSB"].items(): d[f"x{f}"] = [(x["f_sel"], x["t"]) for x in recs]
    else:
        for s in R["bSB"]["1"][0]["pol"]:
            d[f"x1+{int(s) // 1000}k"] = [(x["pol"][s]["f"], x["t"] + x["pol"][s]["t"]) for x in R["bSB"]["1"]]
    return d


SOL = ["SA-frozen", "SA-tuned", "bSB", "bSB+polish"]
RU = {s: reuse_cost(s) for s in SEEDS}
REF = {s: min([RU[s]] + [f for w_ in SOL for runs in samples(EXT[s], w_).values() for f, _ in runs]) for s in SEEDS}


def exp_min_exact(f, N):
    f = np.sort(np.asarray(f)); m = len(f)
    if N < 1 or N > m: return np.nan
    w = np.array([comb(m - k, N - 1) for k in range(1, m + 1)], dtype=float) / comb(m, N)
    return float(np.dot(w, f))


def table(ref):
    res, per = {}, {}
    for s_ in SOL:
        for T in BUDGETS:
            cand = {}
            for e in samples(EXT[SEEDS[0]], s_):
                vals = []; ok = True
                for sd in SEEDS:
                    runs = samples(EXT[sd], s_)[e]; N = int(T // np.mean([x[1] for x in runs]))
                    if N < 1 or len(runs) < 2 * N: ok = False; break
                    vals.append((exp_min_exact([x[0] for x in runs], N) / ref[sd] - 1) * 100)
                if ok: cand[e] = np.array(vals)
            if not cand: res[f"{s_}|{T}"] = None; continue
            e_ = min(cand, key=lambda e: np.median(cand[e])); v = cand[e_]
            rng = np.random.default_rng(T); bs = [np.median(rng.choice(v, len(v))) for _ in range(BOOT)]
            res[f"{s_}|{T}"] = [float(np.median(v)), *np.percentile(bs, [2.5, 97.5]).tolist(), e_]; per[(s_, T)] = v
    return res, per


def holm(ps):
    o = sorted(range(len(ps)), key=lambda i: ps[i]); out = [0.0] * len(ps); mx = 0.0
    for k, i in enumerate(o): mx = max(mx, min(1.0, (len(ps) - k) * ps[i])); out[i] = mx
    return out


def tests(per, a_, b_):
    rows = []
    for T in BUDGETS:
        if (a_, T) not in per or (b_, T) not in per: continue
        dl = per[(a_, T)] - per[(b_, T)]
        rng = np.random.default_rng(100 + T); bs = [np.median(rng.choice(dl, len(dl))) for _ in range(BOOT)]
        p = float(wilcoxon(dl).pvalue) if np.any(dl != 0) else 1.0
        pt = float(max(wilcoxon(dl + MARGIN, alternative="greater").pvalue, wilcoxon(dl - MARGIN, alternative="less").pvalue))
        lo = next((m_ / 100 for m_ in range(1, 301) if max(wilcoxon(dl + m_ / 100, alternative="greater").pvalue,
                                                           wilcoxon(dl - m_ / 100, alternative="less").pvalue) < 0.05), None)
        rows.append(dict(T=T, med_diff=float(np.median(dl)), ci95=np.percentile(bs, [2.5, 97.5]).tolist(), wins_a=int((dl < 0).sum()), N=len(dl),
                         p=p, p_tost=pt, smallest_margin=lo))
    for k, h in zip(rows, holm([r["p"] for r in rows])): k["p_holm"] = h
    for k, h in zip(rows, holm([r["p_tost"] for r in rows])): k["p_tost_holm"] = h
    return rows


res, per = table(REF)
T_pol = tests(per, "SA-tuned", "bSB+polish"); T_fro = tests(per, "SA-frozen", "bSB+polish"); T_bsb = tests(per, "SA-tuned", "bSB")
upto = 0
for T in BUDGETS:
    t_ = next((t for t in T_pol if t["T"] == T), None)
    if t_ is not None and t_["p_tost_holm"] < 0.05: upto = T
    else: break
phi = [x["phi"] for sd in SEEDS for f in EXT[sd]["runs"]["bSB"] for x in EXT[sd]["runs"]["bSB"][f]]
bls = None
if os.path.exists("results/bls_T12_K4.jsonl"):
    seen = {}
    for l in open("results/bls_T12_K4.jsonl"):
        r = json.loads(l); seen.setdefault(r["seed"], r)
    bls = {str(T): float(np.median([(np.mean([x[str(T)] for x in r["runs"]]) / REF[r["seed"]] - 1) * 100 for r in seen.values()])) for T in BUDGETS}
S = dict(n=len(SEEDS), res=res, tuned_vs_polish=T_pol, frozen_vs_polish=T_fro, tuned_vs_bsb=T_bsb, equivalent_up_to=upto,
         reuse_gap=float(np.median([(RU[s] / REF[s] - 1) * 100 for s in SEEDS])), phi_min=float(min(phi)), bls=bls)

# ---- table rows (Table S21)
g = lambda key: res.get(key)
cell = lambda key: "–" if g(key) is None else f"{g(key)[0]:.1f} [{g(key)[1]:.1f}, {g(key)[2]:.1f}]".replace("-", "−")
def dcell(rows, T):
    t_ = next((t for t in rows if t["T"] == T), None)
    if t_ is None: return "–"
    pt = "< 0.001" if t_["p_tost_holm"] < 0.001 else f"{t_['p_tost_holm']:.3f}"
    return f"{t_['med_diff']:+.2f}; {pt}".replace("-", "−")
rows = []
for T in BUDGETS:
    rows.append([f"{T} s"] + [cell(f"{s_}|{T}") for s_ in SOL] + [dcell(T_pol, T)] + ([f"{bls[str(T)]:.1f}".replace("-", "−")] if bls else []))
header = ["Budget", "SA-native (frozen)", "SA-native (tuned)", "bSB", "bSB + polish", "Tuned SA minus bSB + polish; p, equivalence"] + (["BLS"] if bls else [])
widths = [900, 1350, 1350, 1500, 1350, 1650] + ([800] if bls else [])
tp = {t["T"]: t for t in T_pol}
eq_set = [t["T"] for t in T_pol if t["p_tost_holm"] < 0.05]
S["equivalent_set"] = eq_set
S["header"], S["rows"], S["widths"] = header, rows, widths
json.dump(S, open("results/sigma03_K4_summary.json", "w"), indent=1, ensure_ascii=False)
print(json.dumps({k: S[k] for k in ("n", "equivalent_up_to", "reuse_gap", "phi_min", "bls")}, indent=1, ensure_ascii=False))
for r_ in rows: print(r_)
