"""Hardware-neutral comparison (CPU wall-clock may favour annealing).
Work of a run = number of coupling evaluations (multiply-adds over non-zero couplings):
  native SA      : sweeps x n proposals, each scanning the neighbours of one cell -> sweeps x 2|E|
  bSB            : steps x replicas x nnz(J)        (one sparse matrix-vector product per replica and step)
  bSB + polish   : bSB work + polish sweeps x 2|E|
Outputs results/ops.json: CPU time per coupling evaluation, work-matched expected best-of-N gaps, paired tests of
native SA against bSB at matched wall-clock time (own Holm family), and hardware-independent plateau figures."""
import json, glob, numpy as np
from scipy.stats import wilcoxon
from refs import build_ref
from exp_common import setup

rng = np.random.default_rng(17); ref = build_ref()
rows = {}
for f in glob.glob("results/hybrid_K*.jsonl"):
    for l in open(f):
        r = json.loads(l); rows[(r["K"], r["n"], r["seed"])] = r
ext = {}; sat = {}
for f in glob.glob("results/hybext_K*.jsonl"):
    for l in open(f):
        r = json.loads(l); ext[(r["K"], r["n"], r["seed"])] = r["runs"]
for f in glob.glob("results/satuned_K*.jsonl"):
    for l in open(f):
        r = json.loads(l); sat[(r["K"], r["n"], r["seed"])] = r["runs"]
SIZE = {}
for k in rows:
    K, n, s = k; W, csr, a, J, lam, jj = setup(n, s, K)
    SIZE[k] = (len(csr[1]), J.nnz)                     # 2|E| directed neighbour entries, nnz of the QUBO coupling matrix


def runs(k, which):
    """-> {effort: (costs, mean time, work per run)}"""
    r = rows[k]; e2, nnz = SIZE[k]; out = {}
    if which in ("SA-native", "SA-tuned"):
        src = r["runs"]["SA"] if which == "SA-native" else sat.get(k, {})
        for f, recs in src.items():
            out[f"x{f}"] = (np.array([x["f"] for x in recs]), np.mean([x["t"] for x in recs]), int(f) * 1000 * e2)
    for f, recs in r["runs"]["bSB"].items():
        wb = int(f) * 3000 * 8 * nnz
        if which == "bSB":
            out[f"x{f}"] = (np.array([x["f_sel"] for x in recs]), np.mean([x["t"] for x in recs]), wb)
        if which == "bSB+polish":
            for s in recs[0]["pol"]:
                out[f"x{f}+{s}"] = (np.array([x["pol"][s]["f"] for x in recs]), np.mean([x["t"] + x["pol"][s]["t"] for x in recs]), wb + int(s) * e2)
    if which == "bSB+polish" and k in ext:
        for s in ext[k][0]["pol"]:
            out[f"x1+{s}e"] = (np.array([x["pol"][s]["f"] for x in ext[k]]), np.mean([x["t"] + x["pol"][s]["t"] for x in ext[k]]), 3000 * 8 * nnz + int(s) * e2)
    return out


def exp_min(f, N, B=300):
    N = min(N, 1000); return float(np.mean([rng.choice(f, size=N, replace=True).min() for _ in range(B)]))


def matched(keys, which, budget, by):          # by: 1 = time, 2 = work
    per = {}
    for k in keys:
        for e, v in runs(k, which).items():
            N = int(budget // v[by]); per.setdefault(e, {})[k] = (exp_min(v[0], N) / ref[k] - 1) * 100 if N >= 1 else np.nan
    best = None
    for e, d in per.items():
        if len(d) < len(keys): continue
        v = np.array([d[k] for k in keys])
        if np.isnan(v).any(): continue
        if best is None or np.median(v) < np.median(best[1]): best = (e, v)
    return best


groups = sorted({(k[0], k[1]) for k in rows}); out = {"ns_per_eval": {}, "work": {}, "tests_time_sa_vs_bsb": [], "plateau": {}}
for K, n in groups:
    keys = sorted(k for k in rows if k[0] == K and k[1] == n)
    # CPU cost per coupling evaluation
    for which in ("SA-native", "bSB"):
        rr = [v for k in keys for v in runs(k, which).values()]
        out["ns_per_eval"][f"{K}|{n}|{which}"] = float(np.median([v[1] / v[2] * 1e9 for v in rr]))
    # work-matched expected best-of-N
    for Wb in (1e9, 1e10, 1e11):
        for which in ("SA-native", "SA-tuned", "bSB", "bSB+polish"):
            b = matched(keys, which, Wb, 2)
            out["work"][f"{K}|{n}|{which}|{Wb:.0e}"] = None if b is None else [float(np.median(b[1])), b[0]]
    # paired tests SA-native vs bSB at matched wall-clock time
    for T in (1, 3, 10, 30, 100):
        a_ = matched(keys, "SA-native", T, 1); b_ = matched(keys, "bSB", T, 1)
        if a_ and b_:
            dlt = a_[1] - b_[1]
            out["tests_time_sa_vs_bsb"].append(dict(K=K, n=n, T=T, p=float(wilcoxon(dlt).pvalue), med=float(np.median(dlt)), wins=int((dlt < 0).sum()), N=len(dlt)))
    # hardware-independent plateau: best bSB plan over all batches of all efforts, vs one native SA run of 64k sweeps
    best_b = [(min(x["f_sel"] for v in rows[k]["runs"]["bSB"].values() for x in v) / ref[k] - 1) * 100 for k in keys]
    nb = [sum(len(v) for v in rows[k]["runs"]["bSB"].values()) for k in keys]
    sa64 = [np.median([(x["f"] / ref[k] - 1) * 100 for x in rows[k]["runs"]["SA"]["64"]]) for k in keys]
    out["plateau"][f"{K}|{n}"] = dict(best_bsb=float(np.median(best_b)), n_batches=int(np.median(nb)), sa64_single=float(np.median(sa64)),
                                      depth_sa64=64000 * n, depth_bsb=3000)
t = out["tests_time_sa_vs_bsb"]; order = sorted(range(len(t)), key=lambda i: t[i]["p"]); mx = 0.0
for rank, i in enumerate(order):
    mx = max(mx, min(1.0, (len(t) - rank) * t[i]["p"])); t[i]["p_holm"] = mx
json.dump(out, open("results/ops.json", "w"), indent=1)
print(json.dumps(out["ns_per_eval"], indent=0)); print(json.dumps(out["work"], indent=0)); print(json.dumps(out["plateau"], indent=0))
for x in t: print(x["K"], x["n"], x["T"], "holm %.2g med %.2f wins %d/%d" % (x["p_holm"], x["med"], x["wins"], x["N"]))
