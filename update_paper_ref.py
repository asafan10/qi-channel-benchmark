"""Refresh the 'revision' records in paper_run/refrun_K*.jsonl with the current combined reference, then the paper's own
pipeline (run_all.sh) rescored every table against it."""
import json
from refs import build_ref
ref = build_ref(); changed = 0
for K in (3, 4):
    p = f"paper_run/refrun_K{K}.jsonl"
    recs = [json.loads(l) for l in open(p)]
    orig = [r for r in recs if r.get("source") != "revision"]; old = {(r["n"], r["seed"]): r["new_best"] for r in recs if r.get("source") == "revision"}
    with open(p, "w") as fh:
        for r in orig: fh.write(json.dumps(r) + "\n")
        for r in orig:
            v = ref[(K, r["n"], r["seed"])]
            if abs(v - old.get((r["n"], r["seed"]), v)) > 1e-9: changed += 1; print("changed", K, r["n"], r["seed"], old.get((r["n"], r["seed"])), v)
            fh.write(json.dumps(dict(n=r["n"], seed=r["seed"], K=K, prev_bks=v, new_best=v, source="revision")) + "\n")
print("records changed:", changed)
