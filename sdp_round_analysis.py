"""Summary of sdp_round.py. -> results/sdp_round.json
Regular layouts (K = 3, seeds 0-4): rounded plan vs the solver-only reference, the reuse-inclusive reference and the dual-corrected bound.
Irregular layout (sigma = 0.3, seeds 0-29): rounded plan vs the reference of Section 4.11 (best plan of every P6/P10 run), the bound,
and paired Wilcoxon tests against tuned native SA and bSB with refinement at 30 s and 100 s (expected best-of-N of analyse_sigma.py)."""
import json, io, contextlib, numpy as np
from scipy.stats import wilcoxon
import refs

OUTD = {"regular": {}, "sigma": {}}
ref = refs.build_ref(); RU = refs.reuse_data()
lb = {(d["n"], d["seed"]): d["lb_f_rig"] for d in map(json.loads, open("results/sdp_rig.jsonl")) if d["K"] == 3}
REG = [json.loads(l) for l in open("results/sdp_round.jsonl")]
for n in sorted({r["n"] for r in REG}):
    rs = [r for r in REG if r["n"] == n]
    f = {r["seed"]: r["f_best"] for r in rs}
    g = [(f[s] / ref[(3, n, s)] - 1) * 100 for s in f]
    ru = [min(RU[(3, n, s)]) for s in f]
    OUTD["regular"][n] = dict(N=len(rs), gap_ref_median=float(np.median(g)), gap_ref_min=float(min(g)), n_improved=int(sum(x < -1e-7 for x in g)),
                             gap_lb_mean=float(np.mean([(f[s] / lb[(n, s)] - 1) * 100 for s in f])), gap_lb_max=float(max((f[s] / lb[(n, s)] - 1) * 100 for s in f)),
                             vs_reuse_mean=float(np.mean([(f[s] / x - 1) * 100 for s, x in zip(f, ru)])), n_equal_reuse=int(sum(abs(f[s] - x) < 1e-6 for s, x in zip(f, ru))),
                             t_total_mean=float(np.mean([r["t_sdp"] + r["t_round"] + r["t_polish"] for r in rs])))
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    exec(open("analyse_sigma.py").read())
SR = {d["seed"]: d for d in map(json.loads, open("results/sdp_round_sigma0.3.jsonl"))}
seeds_ = [s for s in sorted(SR) if s in ref]
gs = np.array([(SR[s]["f_best"] / ref[s] - 1) * 100 for s in seeds_]); tt = np.array([SR[s]["t_sdp"] + SR[s]["t_round"] + SR[s]["t_polish"] for s in seeds_])
o = dict(N=len(seeds_), gap_ref_median=float(np.median(gs)), gap_ref_min=float(gs.min()), gap_ref_max=float(gs.max()), n_improved=int((gs < -1e-7).sum()),
         gap_lb_mean=float(np.mean([(SR[s]["f_best"] / lb[s] - 1) * 100 for s in seeds_])), ref_lb_mean=float(np.mean([(ref[s] / lb[s] - 1) * 100 for s in seeds_])),
         t_median=float(np.median(tt)), t_max=float(tt.max()), tests=[])
for other in ("SA-tuned", "SA-native", "bSB+polish"):
    for T in (30, 100):
        v = np.asarray(INST[(other, T)]); m = dict(zip([r["seed"] for r in R], v))
        dl = np.array([gs[i] - m[s] for i, s in enumerate(seeds_)])
        o["tests"].append(dict(vs=other, T=T, med_diff=float(np.median(dl)), wins_sdp=int((dl < 0).sum()), p=float(wilcoxon(dl).pvalue), other_median=float(np.median(v))))
OUTD["sigma"] = o
json.dump(OUTD, open("results/sdp_round.json", "w"), indent=1)
print(json.dumps(OUTD, indent=1))
