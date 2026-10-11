"""Summary of the round-3 probes (bench_round3.py A-D) -> results/round3_summary.json.
Gaps are relative to the best plan found for each irregular-layout instance (the reference of Table 15, analyse_sigma_ext.py)."""
import json, re
import numpy as np

src = open("analyse_sigma_ext.py").read()
exec(src[:src.index("LB = {")])                      # REF: same per-instance reference as Table 15
out = {}
g = lambda f, s: (f / REF[s] - 1) * 100

A = json.load(open("results/round3_A.json"))["rows"]
out["A"] = dict(
    n=len(A),
    start_gap_median=float(np.median([g(r["f_start"], r["seed"]) for r in A])),
    by_a0={a0: dict(median_gap=float(np.median([g(min(r["runs"][a0]), r["seed"]) for r in A])),
                    median_change_vs_start=float(np.median([(min(r["runs"][a0]) / r["f_start"] - 1) * 100 for r in A])),
                    improved=sum(min(r["runs"][a0]) < r["f_start"] - 1e-9 for r in A),
                    kept=sum(abs(min(r["runs"][a0]) - r["f_start"]) <= 1e-9 for r in A),
                    worse=sum(min(r["runs"][a0]) > r["f_start"] + 1e-9 for r in A),
                    batches_identical=all(max(r["runs"][a0]) - min(r["runs"][a0]) < 1e-9 for r in A))
               for a0 in ("0.0", "0.5", "0.9")})

B = json.load(open("results/round3_B.json"))["rows"]; names = list(B[0]["res"])
out["B"] = {k: dict(median_gap=float(np.median([g(min(r["res"][k]["f"]), r["seed"]) for r in B])),
                    change_vs_default=[float((min(r["res"][k]["f"]) / min(r["res"]["default"]["f"]) - 1) * 100) for r in B],
                    phi=float(np.mean([np.mean(r["res"][k]["phi"]) for r in B])), t=float(np.mean([r["res"][k]["t"] for r in B])))
            for k in names}

C = json.load(open("results/round3_C.json"))["rows"]
dev = {k: [abs(np.mean(r["res"][k]["f_perm"]) / r["res"][k]["f_index"] - 1) * 100 for r in C] for k in ("SA-QUBO", "bSB")}
spread = {k: [(max(r["res"][k]["f_perm"] + [r["res"][k]["f_index"]]) / min(r["res"][k]["f_perm"] + [r["res"][k]["f_index"]]) - 1) * 100 for r in C]
          for k in ("SA-QUBO", "bSB")}
out["C"] = dict(n=len(C), phi_range={k: [min(r["res"][k]["phi"] for r in C), max(r["res"][k]["phi"] for r in C)] for k in dev},
                max_mean_dev_pct={k: float(max(v)) for k, v in dev.items()}, max_spread_pct={k: float(max(v)) for k, v in spread.items()},
                bsb_better_than_saqubo=sum(r["res"]["bSB"]["f_index"] < r["res"]["SA-QUBO"]["f_index"] for r in C))
out["D"] = json.load(open("results/round3_D.json"))
json.dump(out, open("results/round3_summary.json", "w"), indent=1)
print(json.dumps(out, indent=1))
