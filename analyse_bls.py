"""Breakout local search at matched wall-clock time (bench_bls.py) against tuned native SA and bSB (Table 7 pools,
exact expected best of N) and against tuned native SA and bSB with refinement on the irregular layout (Table 12). The BLS value at budget T
is the mean over runs of the best cost found within T s of one run, i.e. the expected result of one run of length T.
-> results/bls_summary.json"""
import json, glob, numpy as np
from refs import build_ref

BUD = ["1", "3", "10", "30", "100"]; REF = build_ref(); S = {}
med = lambda v: float(np.median(v))
T7 = json.load(open("results/fresh67_summary.json"))["table7"]
for K, n in ((3, 200), (3, 500), (4, 200)):
    rows = [json.loads(l) for l in open(f"results/bls_T7_K{K}_{n}.jsonl")]
    g = {T: [(np.mean([r_[T] for r_ in r["runs"]]) / REF[(K, n, r["seed"])] - 1) * 100 for r in rows] for T in BUD}
    S[f"{K}_{n}"] = dict(n=len(rows), bls={T: med(v) for T, v in g.items()},
                         sa_tuned={T: T7[f"{K}_{n}"]["SA-tuned"][T]["gap"] for T in BUD if T in T7[f"{K}_{n}"]["SA-tuned"]},
                         bsb={T: T7[f"{K}_{n}"]["bSB"][T]["gap"] for T in BUD if T in T7[f"{K}_{n}"]["bSB"]})
src = open("analyse_sigma_ext.py").read(); gl = {}; exec(src[:src.index("LB = {")], gl); REFI = gl["REF"]
_seen = {}
for l in open("results/bls_T12.jsonl"):
    r = json.loads(l); _seen.setdefault(r["seed"], r)          # first record per seed (two processes shared the last seeds)
rows = list(_seen.values())
g12 = {T: [(np.mean([r_[T] for r_ in r["runs"]]) / REFI[r["seed"]] - 1) * 100 for r in rows] for T in BUD}
P = json.load(open("results/sigma03_ext.json"))["primary"]["res"]
S["T12"] = dict(n=len(rows), bls={T: med(v) for T, v in g12.items()}, below_ref={T: int(sum(x < 0 for x in v)) for T, v in g12.items()},
                sa_tuned={T: P[f"SA-tuned|{T}"][0] for T in BUD}, polish={T: P[f"bSB+polish|{T}"][0] for T in BUD if f"bSB+polish|{T}" in P})
json.dump(S, open("results/bls_summary.json", "w"), indent=1)
print(json.dumps(S, indent=1))
