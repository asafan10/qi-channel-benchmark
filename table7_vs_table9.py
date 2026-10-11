"""Why Table 7 (protocol P4, machine A) and Table 9 (protocol P6, machine B) differ for native SA at K = 3, 200 cells, 100 s on the
same five instances (seeds 0-4). Both use the final best-known reference. For each effort: pool size, time per run,
single-run gaps, the best stored run per instance, and a comparison of the run distributions centred per instance.
-> results/table7_vs_table9.json"""
import json, numpy as np, refs
from scipy.stats import mannwhitneyu
REF = refs.build_ref(); ref = lambda s: REF[(3, 200, s)]
L = {r["seed"]: r for r in map(json.loads, open("paper_run/longrun_K3.jsonl")) if r["n"] == 200}
SW = {r["seed"]: r for r in map(json.loads, open("paper_run/sweep_K3.jsonl")) if r["n"] == 200}
H = {r["seed"]: r for r in map(json.loads, open("results/hybrid_K3.jsonl")) if r["n"] == 200 and r["seed"] < 5}


def t7(s):
    d = {float(fk): (np.array(c["f"]), float(np.mean(c["t"]))) for fk, c in SW[s]["runs"]["SA-native"].items()}
    d.update({float(fk): (np.array([x["f_sel"] for x in recs]), float(np.mean([x["t"] for x in recs]))) for fk, recs in L[s]["runs"]["SA-native"].items()})
    return d


def t9(s): return {float(f): (np.array([x["f"] for x in recs]), float(np.mean([x["t"] for x in recs]))) for f, recs in H[s]["runs"]["SA"].items()}


out = {}
for name, src in (("table7", t7), ("table9", t9)):
    for fk in sorted(src(0)):
        g = [(src(s)[fk][0] / ref(s) - 1) * 100 for s in range(5)]
        out[f"{name}|x{fk:g}"] = dict(pool=[len(x) for x in g], t_run=float(np.mean([src(s)[fk][1] for s in range(5)])),
                                      N_at_100s=int(100 // np.mean([src(s)[fk][1] for s in range(5)])),
                                      run_mean=float(np.mean(np.concatenate(g))), single_median=float(np.median([np.median(x) for x in g])),
                                      best_stored_median=float(np.median([x.min() for x in g])), runs=[x.round(3).tolist() for x in g])
for fk in (64.0, 256.0):
    a = [(t7(s)[fk][0] / ref(s) - 1) * 100 for s in range(5)]; b = [(t9(s)[fk][0] / ref(s) - 1) * 100 for s in range(5)]
    c = [np.median(np.concatenate([x, y])) for x, y in zip(a, b)]
    out[f"compare|x{fk:g}"] = dict(mean_t7=float(np.mean(np.concatenate(a))), mean_t9=float(np.mean(np.concatenate(b))),
                                   p_mwu_centred=float(mannwhitneyu(np.concatenate([x - m for x, m in zip(a, c)]), np.concatenate([y - m for y, m in zip(b, c)])).pvalue))
json.dump(out, open("results/table7_vs_table9.json", "w"), indent=1)
for k, v in out.items(): print(k, {kk: vv for kk, vv in v.items() if kk != "runs"})
