"""Analysis of bench_round12.py on the irregular layout of Table 12 (K = 3, 200 cells, 30 instances), with the
reference, estimator and eligibility rule of Table 12 (the new runs do not enter the reference):
 (a) tuned PT-QUBO: exact expected best of N at 1-30 s, best eligible effort, paired two-sided Wilcoxon tests against tuned native SA,
     Holm over the four budgets (post hoc); compared descriptively with SQA tuned on the same instances
 (b) the polish from six starts: median over instances of the mean single-run gap after 16,000 and 64,000 polish sweeps; paired Wilcoxon
     tests of the bSB start against each other start at both lengths, Holm over the ten tests (post hoc)
 (c) bSB x1 followed by 256,000 polish sweeps: expected best of N at 30 and 100 s against tuned native SA and against the capped polish of
     Table 12 (paired Wilcoxon, Holm over the four tests, post hoc)
 (d) instance statistics at 200 cells on the near-regular layout (Table 2) and the intra-cell to first-neighbour coupling ratio
-> results/round12_summary.json"""
import io, json, os, contextlib, numpy as np
from math import comb
from scipy.stats import wilcoxon

src = open("analyse_sigma_ext.py").read(); ns = {"__name__": "defs"}
with contextlib.redirect_stdout(io.StringIO()): exec(src[:src.index('OUT = {"seeds"')], ns)
REF = ns["REF"]; SEEDS = list(ns["SEEDS"]); _, per = ns["table"](REF)
LP = json.load(open("results/lowpen_summary.json"))


def emin(f, N):
    f = np.sort(np.asarray(f)); m = len(f)
    if N < 1 or 2 * N > m: return np.nan
    w = np.array([comb(m - k, N - 1) for k in range(1, m + 1)], dtype=float) / comb(m, N)
    return float(np.dot(w, f))


def holm(ps):
    o = sorted(range(len(ps)), key=lambda i: ps[i]); out = [0.0] * len(ps); mx = 0.0
    for k, i in enumerate(o): mx = max(mx, min(1.0, (len(ps) - k) * ps[i])); out[i] = mx
    return out


def best_eff(P_, seeds, T):
    best = None
    for e in P_[seeds[0]]:
        v = np.array([(emin([f for f, _ in P_[s][e]], int(T // np.mean([t for _, t in P_[s][e]]))) / REF[s] - 1) * 100 for s in seeds])
        if np.isnan(v).any(): continue
        if best is None or np.median(v) < np.median(best[1]): best = (e, v)
    return best


def load(path):
    rr = {}
    for l in open(path): r = json.loads(l); rr[r["seed"]] = r
    return rr


out = {}
# (a) PT-QUBO
tune = json.load(open("results/r12_pttune.json")); pt = load("results/r12_ptpool.jsonl")
assert all(s in pt for s in SEEDS), "PT pools incomplete"
P_ = {s: {e: [(x["f"], x["t"]) for x in v] for e, v in pt[s]["runs"].items()} for s in SEEDS}
tests = []
for T in (1, 3, 10, 30):
    b = best_eff(P_, SEEDS, T)
    if b is None: continue
    sa = np.asarray(per[("SA-tuned", T)]); d = b[1] - sa
    sq = next((t for t in LP["b"]["tests"] if t["solver"] == "SQA" and t["tuned"] == "irregular" and t["T"] == T), None)
    tests.append(dict(T=T, effort=b[0], med=float(np.median(b[1])), sa=float(np.median(sa)), med_diff=float(np.median(d)),
                      wins=int((d < 0).sum()), p=float(wilcoxon(d).pvalue), sqa=sq["med"] if sq else None))
for t, h in zip(tests, holm([t["p"] for t in tests])): t["p_holm"] = h
_cfgs = [g["cfg"] for g in tune["grid"].values()]
def _side(par):   # tested values below and above the best setting's value, other parameters free (lam/2 upper temperatures excluded)
    v = [c[par] for c in _cfgs if not isinstance(c[par], str)]; b = tune["best"][par]
    return dict(below=any(x < b for x in v), above=any(x > b for x in v))
_edge = {par: _side(par) for par in ("R", "rho", "Tmin", "Tmax")}
out["pt"] = dict(cfg=tune["best"], grid=len(tune["grid"]), edge=_edge, extensions=sum(bool(tune.get(k)) for k in ("extended", "extended2", "extended3")), prev=dict(cfg=tune.get("best_before"), gap=next((g["mean_gap"] for g in tune["grid"].values() if g["cfg"] == tune.get("best_before")), None)), heldout_gap=tune["best_gap"],
                 phi=float(np.mean([x["phi"] for s in SEEDS for v in pt[s]["runs"].values() for x in v])), tests=tests)
# (b) polish from other starts
ps = load("results/r12_polstart.jsonl"); assert all(s in ps for s in SEEDS), "polish-start runs incomplete"
names = ["random", "row-order", "greedy", "DSATUR", "SDP rounding", "bSB"]
gap = lambda s, f: (f / REF[s] - 1) * 100
pol = {}
for nm in names:
    pol[nm] = {}
    for sw in ("16000", "64000"):
        v = np.array([np.mean([gap(s, x["f"]) for x in ps[s]["runs"][nm][sw]]) for s in SEEDS])
        pol[nm][sw] = dict(med=float(np.median(v)), per=v.tolist(), t=float(np.mean([x["t"] for s in SEEDS for x in ps[s]["runs"][nm][sw]])))
    if nm in ps[SEEDS[0]]["start_cost"]:
        pol[nm]["start"] = float(np.median([gap(s, ps[s]["start_cost"][nm]) for s in SEEDS]))
    elif nm == "bSB":
        pol[nm]["start"] = float(np.median([np.mean([gap(s, x["f_start"]) for x in ps[s]["runs"]["bSB"]["64000"]]) for s in SEEDS]))
ptests = []
for sw in ("16000", "64000"):
    for nm in names[:-1]:
        d = np.array(pol["bSB"][sw]["per"]) - np.array(pol[nm][sw]["per"])
        ptests.append(dict(sweeps=int(sw), other=nm, med_diff=float(np.median(d)), wins_bsb=int((d < 0).sum()), p=float(wilcoxon(d).pvalue)))
for t, h in zip(ptests, holm([t["p"] for t in ptests])): t["p_holm"] = h
out["polstart"] = dict(starts=pol, tests=ptests, t_sdp=float(np.median([ps[s]["t_sdp"] for s in SEEDS])))
# (c) the long polish
lp = load("results/r12_longpol.jsonl"); assert all(s in lp for s in SEEDS), "long-polish runs incomplete"
ltests = []
for T in (30, 100):
    v = np.array([gap(s, emin([x["f"] for x in lp[s]["runs"]], int(T // np.mean([x["t"] for x in lp[s]["runs"]])))) for s in SEEDS])
    if np.isnan(v).any(): continue
    for other in ("SA-tuned", "bSB+polish"):
        o = np.asarray(per[(other, T)]); d = v - o
        ltests.append(dict(T=T, other=other, med=float(np.median(v)), other_med=float(np.median(o)), med_diff=float(np.median(d)),
                           wins=int((d < 0).sum()), p=float(wilcoxon(d).pvalue)))
for t, h in zip(ltests, holm([t["p"] for t in ltests])): t["p_holm"] = h
out["longpol"] = dict(t=float(np.mean([x["t"] for s in SEEDS for x in lp[s]["runs"]])),
                      single=float(np.median([np.median([gap(s, x["f"]) for x in lp[s]["runs"]]) for s in SEEDS])), tests=ltests)
# (d) instance statistics (near-regular layout, 200 cells, seeds 0-4)
from bench_ran_qubo import make_instance, build_qubo
nn, hl, deg, dw = [], [], [], []
for s in range(5):
    W = make_instance(200, s, 3); a, J, lam = build_qubo(W, 3, 0.25); hl.append(lam / 2)
    nn.append(W[np.triu(W, 1) > 1.37 ** -3.5].mean())          # first shell: distance below the midpoint of 1 and sqrt(3)
    deg.append(float((W > 0).sum(1).mean())); dw.append(float(W.sum(1).max()))
_bysize = {}
for n_ in (20, 50, 100, 200, 500):   # Table 2: mean degree and largest weighted degree per size (near-regular layout, seeds 0-4)
    dg_, dw_ = [], []
    for s in range(5):
        W_ = make_instance(n_, s, 3); dg_.append(float((W_ > 0).sum(1).mean())); dw_.append(float(W_.sum(1).max()))
    _bysize[n_] = dict(deg=float(np.mean(dg_)), deg_lo=min(dg_), deg_hi=max(dg_), dw_lo=min(dw_), dw_hi=max(dw_))
out["instance_stats"] = dict(by_size=_bysize, half_lam=float(np.median(hl)), nn_w=float(np.mean(nn)), ratio=float(np.median(np.array(hl) / (np.array(nn) / 4))),
                             deg=float(np.mean(deg)), deg_lo=min(deg), deg_hi=max(deg), dw_lo=min(dw), dw_hi=max(dw))
json.dump(out, open("results/round12_summary.json", "w"), indent=1)
print(json.dumps({k: (v if k != "polstart" else {n: {kk: (vv if kk != "16000" and kk != "64000" else round(vv["med"], 2)) for kk, vv in d.items()} for n, d in v["starts"].items()}) for k, v in out.items()}, indent=1, default=str)[:4000])
