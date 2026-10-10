"""Dual-corrected semidefinite bounds for the remaining P6 instances at 200 cells (seeds 5-29, K = 3 and 4), SCS tolerance 1e-4,
used to normalise paired differences by a quantity no solver influences (Section 4.10). -> results/sdp_rig_p6.jsonl"""
import json, sys
from ortools.sat.python import cp_model  # noqa: F401
from bench_sdp_rig import rig
from bench_ran_qubo import make_instance
K = int(sys.argv[1]); done = set()
try:
    for l in open("results/sdp_rig_p6.jsonl"): r = json.loads(l); done.add((r["K"], r["n"], r["seed"]))
except FileNotFoundError: pass
for s in range(5, 30):
    if (K, 200, s) in done: continue
    r = rig(make_instance(200, s, K), K, eps=1e-4)
    with open("results/sdp_rig_p6.jsonl", "a") as fh: fh.write(json.dumps(dict(K=K, n=200, seed=s, eps=1e-4, **r)) + "\n")
    print(K, s, r["lb_f_rig"], r["t"], flush=True)
