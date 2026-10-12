"""Reference-bias check at every budget. The paired tests of native SA against bSB with refinement are repeated
(1) with every gap measured against the reference built without native SA's runs (refs.build_ref(no_native_sa=True)), and
(2) on raw costs, which use no reference at all, for the regular layouts (protocol P6, Table 13, same efforts) and for the irregular
layout (Table 15, same eligible efforts; there the SA-free reference is the best bSB or bSB-alone run of the instance).
Verdicts (Holm within each family) are compared with the primary analysis. -> results/nosa_ref_check.json"""
import io, json, contextlib, numpy as np
from scipy.stats import wilcoxon
with contextlib.redirect_stdout(io.StringIO()):
    exec(open("analyse_hybrid.py").read())


def holm(ps):
    o = sorted(range(len(ps)), key=lambda i: ps[i]); out = [0.0] * len(ps); mx = 0.0
    for r, i in enumerate(o): mx = max(mx, min(1.0, (len(ps) - r) * ps[i])); out[i] = mx
    return out


def verdict(p, med): return "ns" if p >= 0.05 else ("SA" if med < 0 else "polish")


out = {"regular": [], "irregular": []}
ROB = {(t["K"], t["n"], t["T"], t["a"]): t for t in tests}
rows_ = []
for K, n in groups:
    for T in BUDGETS:
        for a_ in ("SA-native", "SA-tuned"):
            A = INST.get((K, n, a_, T)); P = INST.get((K, n, "bSB+polish", T))
            if not A or not P: continue
            seeds = A[1]; ga, gp = np.array(A[2]), np.array(P[2])
            ca = np.array([ref[(K, n, s)] * (1 + g / 100) for s, g in zip(seeds, ga)]); cp = np.array([ref[(K, n, s)] * (1 + g / 100) for s, g in zip(seeds, gp)])
            rn = np.array([ref_nosa[(K, n, s)] for s in seeds])
            dn = (ca - cp) / rn * 100; dr = ca - cp
            rows_.append(dict(K=K, n=n, T=T, a=a_, med_nosa=float(np.median(dn)), p_nosa=float(wilcoxon(dn).pvalue), med_raw=float(np.median(dr)),
                              p_raw=float(wilcoxon(dr).pvalue), p_primary=ROB[(K, n, T, a_)]["p_holm"], med_primary=ROB[(K, n, T, a_)]["med_diff"],
                              gap_nosa_a=float(np.median((ca / rn - 1) * 100)), gap_nosa_pol=float(np.median((cp / rn - 1) * 100))))
for key in ("p_nosa", "p_raw"):
    for r, h in zip(rows_, holm([r[key] for r in rows_])): r[key + "_holm"] = h
for r in rows_:
    r["v_primary"] = verdict(r["p_primary"], r["med_primary"]); r["v_nosa"] = verdict(r["p_nosa_holm"], r["med_nosa"]); r["v_raw"] = verdict(r["p_raw_holm"], r["med_raw"])
out["regular"] = rows_
# irregular layout
src = open("analyse_sigma_ext.py").read(); exec(src[:src.index("LB = {")])
SG = json.load(open("results/sigma03.json")); eff = {k: v[3] for k, v in SG["res"].items()}
REF_NOSA = {}
for s in SEEDS:
    v = [x[0] for e in samples(EXT[s], "bSB").values() for x in e]
    REF_NOSA[s] = min(v)
from math import comb
def exp_min_exact(f, N):
    f = np.sort(np.asarray(f)); m = len(f)
    if N < 1 or N > m: return np.nan
    w_ = np.array([comb(m - k, N - 1) for k in range(1, m + 1)], dtype=float) / comb(m, N); return float(np.dot(w_, f))
def cost(s, which, T):
    e = eff[f"{which}|{T}"]; d = samples(EXT[s], {"SA-native": "SA-frozen"}.get(which, which))[e]
    return exp_min_exact([x[0] for x in d], int(T // np.mean([x[1] for x in d])))
irr = []
for T in BUDGETS:
    for a_ in ("SA-tuned", "SA-native"):
        ca = np.array([cost(s, a_, T) for s in SEEDS]); cp = np.array([cost(s, "bSB+polish", T) for s in SEEDS]); rn = np.array([REF_NOSA[s] for s in SEEDS])
        dn = (ca - cp) / rn * 100; dr = ca - cp
        irr.append(dict(T=T, a=a_, med_nosa=float(np.median(dn)), p_nosa=float(wilcoxon(dn).pvalue), med_raw=float(np.median(dr)), p_raw=float(wilcoxon(dr).pvalue),
                        gap_nosa_a=float(np.median((ca / rn - 1) * 100)), gap_nosa_pol=float(np.median((cp / rn - 1) * 100))))
for key in ("p_nosa", "p_raw"):
    for a_ in ("SA-tuned", "SA-native"):
        sub = [r for r in irr if r["a"] == a_]
        for r, h in zip(sub, holm([r[key] for r in sub])): r[key + "_holm"] = h
tp = {(t["a"], t["T"]): t for t in SG["tests_vs_polish"]}
for r in irr:
    t = tp[({"SA-native": "SA-native"}.get(r["a"], r["a"]), r["T"])]
    r["v_primary"] = verdict(t["p_holm"], t["med_diff"]); r["v_nosa"] = verdict(r["p_nosa_holm"], r["med_nosa"]); r["v_raw"] = verdict(r["p_raw_holm"], r["med_raw"])
out["irregular"] = irr
out["summary"] = {k: dict(n=len(out[k]), same_nosa=sum(r["v_nosa"] == r["v_primary"] for r in out[k]), same_raw=sum(r["v_raw"] == r["v_primary"] for r in out[k]))
                  for k in ("regular", "irregular")}
json.dump(out, open("results/nosa_ref_check.json", "w"), indent=1)
print(out["summary"])
for k in ("regular", "irregular"):
    for r in out[k]:
        if r["v_nosa"] != r["v_primary"] or r["v_raw"] != r["v_primary"]: print(k, {x: r[x] for x in r if x.startswith(("K", "n", "T", "a", "v_", "med", "p_"))})
