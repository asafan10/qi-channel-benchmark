"""Rigorous semidefinite bounds at 500 cells (K = 3) with SCS tolerance 1e-4. The dual-corrected bound is valid at any
tolerance; at 200 cells 1e-4 was 42 times faster than 1e-6 and 0.4% looser."""
import json, sys
from ortools.sat.python import cp_model  # noqa: F401  (import before cvxpy)
from bench_sdp_rig import rig
from bench_ran_qubo import make_instance
for s in map(int, sys.argv[1].split(",")):
    r = rig(make_instance(500, s, 3), 3, eps=1e-4)
    with open("results/sdp_rig.jsonl", "a") as fh: fh.write(json.dumps(dict(K=3, n=500, seed=s, eps=1e-4, **r)) + "\n")
    print(500, s, r["lb_f_rig"], r["slack"], r["t"], flush=True)
