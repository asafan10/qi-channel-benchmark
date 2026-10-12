"""Pool-size check of the regular-layout P6 comparisons.
Applies the rule of the irregular-layout analysis (PREREG_sigma03_ext.md): the exact expected best of N without replacement, and an effort
is used at a budget only if its pool holds m >= 2N runs on every instance. For each setting, solver and budget: the effort used in the paper
(resampling with replacement, any pool), its pool and N, whether it is eligible, the best eligible effort, and the paired tests of native SA
(frozen, tuned) against bSB with refinement restricted to eligible efforts (Wilcoxon, Holm over the tests that remain).
-> results/p6_eligibility.json"""
import io, json, contextlib, numpy as np
from math import comb
from scipy.stats import wilcoxon
with contextlib.redirect_stdout(io.StringIO()):
    exec(open("analyse_hybrid.py").read())


def exp_min_exact(f, N):
    f = np.sort(np.asarray(f)); m = len(f)
    if N < 1 or N > m: return np.nan
    w = np.array([comb(m - k, N - 1) for k in range(1, m + 1)], dtype=float) / comb(m, N)
    return float(np.dot(w, f))


ROB = {(t["K"], t["n"], t["T"], t["a"]): t for t in json.load(open("results/robust.json"))["tests"]}
out = {"rows": {}, "tests": []}; EL = {}
for K, n in groups:
    sel = sorted([r for r in rows if r["K"] == K and r["n"] == n], key=lambda r: r["seed"])
    for s in ("SA-native", "SA-tuned", "bSB", "bSB+polish"):
        if s == "SA-tuned" and not all(r.get("sat") for r in sel): continue
        for T in BUDGETS:
            used = res.get((K, n, s, T)); info = {}
            for e in samples(sel[0], s):
                ms = [len(samples(r, s)[e][0]) for r in sel]; Ns = [int(T // samples(r, s)[e][1]) for r in sel]
                info[e] = dict(m_min=min(ms), N_max=max(Ns), N_min=min(Ns), eligible=all(N >= 1 and m >= 2 * N for m, N in zip(ms, Ns)))
            best = None
            for e, v in info.items():
                if not v["eligible"]: continue
                g = np.array([(exp_min_exact(samples(r, s)[e][0], int(T // samples(r, s)[e][1])) / ref[(K, n, r["seed"])] - 1) * 100 for r in sel])
                if best is None or np.median(g) < best[1]: best = (e, float(np.median(g)), g)
            ue = used[3] if used else None
            out["rows"][f"{K}|{n}|{s}|{T}"] = dict(effort_used=ue, used_eligible=bool(ue and info[ue]["eligible"]),
                                                    m_used=info[ue]["m_min"] if ue else None, N_used=info[ue]["N_max"] if ue else None,
                                                    best_eligible=best[0] if best else None, gap_eligible=best[1] if best else None,
                                                    gap_used=float(used[0]) if used else None)
            if best: EL[(K, n, s, T)] = best[2]
for K, n in groups:
    for T in BUDGETS:
        for a_ in ("SA-native", "SA-tuned"):
            if (K, n, a_, T) not in EL or (K, n, "bSB+polish", T) not in EL: continue
            d = EL[(K, n, a_, T)] - EL[(K, n, "bSB+polish", T)]
            out["tests"].append(dict(K=K, n=n, T=T, a=a_, med_diff=float(np.median(d)), wins_a=int((d < 0).sum()),
                                     p=float(wilcoxon(d).pvalue) if np.any(d != 0) else 1.0, p_holm_paper=ROB[(K, n, T, a_)]["p_holm"],
                                     med_diff_paper=ROB[(K, n, T, a_)]["med_diff"]))
tt = out["tests"]; order = sorted(range(len(tt)), key=lambda i: tt[i]["p"]); mx = 0.0
for rank, i in enumerate(order): mx = max(mx, min(1.0, (len(tt) - rank) * tt[i]["p"])); tt[i]["p_holm"] = mx
nrow = len(out["rows"]); nel = sum(v["used_eligible"] for v in out["rows"].values())
out["summary"] = dict(rows=nrow, rows_used_effort_eligible=nel, tests_paper=len(ROB), tests_eligible=len(tt))
json.dump(out, open("results/p6_eligibility.json", "w"), indent=1)
print(out["summary"])
for k, v in out["rows"].items(): print(k, v)
for t in tt: print(t["K"], t["n"], t["T"], t["a"], f"{t['med_diff']:+.2f} wins {t['wins_a']} p_holm {t['p_holm']:.3f} | paper {t['med_diff_paper']:+.2f} {t['p_holm_paper']:.3f}")
