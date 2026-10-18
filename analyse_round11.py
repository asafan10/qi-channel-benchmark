"""(a) bSB at 256 times its baseline steps (Table 6), (b) ten repeats of protocol P1 for native SA at 100 cells
(the 2-point reading threshold). -> results/round11_summary.json"""
import json, os, numpy as np
from refs import build_ref
REF = build_ref(); out = {"bsb256": {}, "sarep": {}}
rows = []
for f in ("results/r11_bsb256.jsonl", "results/r11_bsb256_500.jsonl"):
    if os.path.exists(f): rows += [json.loads(l) for l in open(f)]
F67 = json.load(open("results/fresh67_summary.json"))
for K, n in ((3, 200), (3, 500), (4, 200)):
    rr = {r["seed"]: r for r in rows if (r["K"], r["n"]) == (K, n)}
    if not rr: continue
    per = {s: float(np.median([(x["f"] / REF[(K, n, s)] - 1) * 100 for x in r["runs"]])) for s, r in rr.items()}
    out["bsb256"][f"{K}_{n}"] = dict(n_inst=len(rr), gap=float(np.median(list(per.values()))), per=per,
                                     t=float(np.mean([np.mean([x["t"] for x in r["runs"]]) for r in rr.values()])),
                                     phi=float(np.mean([x["phi"] for r in rr.values() for x in r["runs"]])),
                                     gap64=F67["table6"][f"{K}_{n}"]["bSB"]["64"]["gap"] if f"{K}_{n}" in F67.get("table6", {}) else None)
sa = [json.loads(l) for l in open("results/r11_sarep.jsonl")]
reps = sorted({r["rep"] for r in sa}); seeds = sorted({r["seed"] for r in sa})
g = {(r["rep"], r["seed"]): (r["f"] / REF[(3, 100, r["seed"])] - 1) * 100 for r in sa}
mean_r = [float(np.mean([g[(k, s)] for s in seeds])) for k in reps]
med_r = [float(np.median([g[(k, s)] for s in seeds])) for k in reps]
out["sarep"] = dict(reps=len(reps), seeds=seeds, mean_by_rep=mean_r, median_by_rep=med_r,
                    mean_sd=float(np.std(mean_r, ddof=1)), mean_range=[min(mean_r), max(mean_r)],
                    median_sd=float(np.std(med_r, ddof=1)), median_range=[min(med_r), max(med_r)],
                    per_instance_sd=[float(np.std([g[(k, s)] for k in reps], ddof=1)) for s in seeds],
                    restarts=float(np.mean([r["restarts"] for r in sa])))
json.dump(out, open("results/round11_summary.json", "w"), indent=1)
print(json.dumps(out, indent=1)[:3000])
