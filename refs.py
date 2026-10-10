"""Single best-known reference per (K, n, seed): minimum cost found by any run of any experiment available here
(additional experiments + the paper's tts_K*.jsonl + final_reference.json if supplied). Also loads SDP lower bounds."""
import json, glob, os
import numpy as np


def build_ref(no_native_sa=False):
    """no_native_sa=True: reference from runs of every solver EXCEPT native-space SA (no SA-native runs, no SA polish,
    no reference runs or stored best-knowns, whose provenance includes SA). Used to show the conclusions do not depend
    on scoring native SA against its own runs."""
    ref = {}
    def upd(k, v): ref[k] = min(ref.get(k, np.inf), float(v))
    if no_native_sa:
        return _ref_no_sa(upd, ref)
    for f in glob.glob("results/hybrid_K*.jsonl") + glob.glob("results/hybext_K*.jsonl") + glob.glob("results/satuned_K*.jsonl"):
        for l in open(f):
            r = json.loads(l); upd((r["K"], r["n"], r["seed"]), r["best"])
    for f in glob.glob("results/core_pol_t0.jsonl"):
        for l in open(f):
            r = json.loads(l)
            for v in r["f_sel"]: upd((r["K"], r["n"], r["seed"]), v)
            for v in r["runs"].values():
                for x in v: upd((r["K"], r["n"], r["seed"]), x["f"])
    for f in glob.glob("results/core_sbdiag.jsonl"):
        for l in open(f):
            r = json.loads(l)
            for v in r["runs"].values():
                for x in v: upd((r["K"], r["n"], r["seed"]), x["f"])
    for f in glob.glob("results/sqa_K*.jsonl"):
        for l in open(f):
            r = json.loads(l)
            for v in r["res"].values(): upd((r["K"], r["n"], r["seed"]), v["E"])
    for f in glob.glob("results/sbpkg_K*.jsonl"):
        for l in open(f):
            r = json.loads(l)
            for v in r["runs"].values():
                for x in v: upd((r["K"], r["n"], r["seed"]), x["f"])
    for K in (3, 4):
        p = f"paper_data/tts_K{K}.jsonl"
        if os.path.exists(p):
            for l in open(p):
                r = json.loads(l); upd((K, r["n"], r["seed"]), r["prev_bks"]); upd((K, r["n"], r["seed"]), r["new_best"])
    for p in glob.glob("paper_data/longrun_K*.jsonl") + glob.glob("paper_data/tuned_K*.jsonl") + glob.glob("paper_data/sweep_K*.jsonl"):
        for l in open(p):
            r = json.loads(l); upd((r["K"], r["n"], r["seed"]), r["new_best"])
    for p in glob.glob("paper_data/refrun_K*.jsonl"):
        for l in open(p):
            r = json.loads(l); upd((r["K"], r["n"], r["seed"]), r["prev_bks"]); upd((r["K"], r["n"], r["seed"]), r["new_best"])
    for p in glob.glob("paper_data/results_K*.jsonl"):
        for l in open(p):
            r = json.loads(l); K = r.get("K", int(p.split("_K")[1][0]))
            for v in r["res"].values(): upd((K, r["n"], r["seed"]), v["E"])
    p = "paper_data/final_reference.json"
    if os.path.exists(p):
        for Kk, d in json.load(open(p)).items():
            for ns, v in d.items():
                n, s = map(int, ns.split("_")); upd((int(Kk[1:]), n, s), v)
    return ref


def _ref_no_sa(upd, ref):
    for f in glob.glob("results/hybrid_K*.jsonl"):
        for l in open(f):
            r = json.loads(l); k = (r["K"], r["n"], r["seed"])
            for v in r["runs"]["bSB"].values():
                for x in v: upd(k, x["f_sel"]); upd(k, x["f_desc"])
    for f in glob.glob("results/hybext_K*.jsonl"):
        for l in open(f):
            r = json.loads(l); k = (r["K"], r["n"], r["seed"])
            for x in r["runs"]: upd(k, x["f_sel"])
    for f in glob.glob("results/core_pol_t0.jsonl"):
        for l in open(f):
            r = json.loads(l)
            for v in r["f_sel"]: upd((r["K"], r["n"], r["seed"]), v)
    for f in glob.glob("results/core_sbdiag.jsonl") + glob.glob("results/sbpkg_K*.jsonl"):
        for l in open(f):
            r = json.loads(l)
            for v in r["runs"].values():
                for x in v: upd((r["K"], r["n"], r["seed"]), x["f"])
    for f in glob.glob("results/sqa_K*.jsonl"):
        for l in open(f):
            r = json.loads(l)
            for s, v in r["res"].items():
                if s != "SA-native": upd((r["K"], r["n"], r["seed"]), v["E"])
    for K in (3, 4):
        p = f"paper_data/tts_K{K}.jsonl"
        if os.path.exists(p):
            for l in open(p):
                r = json.loads(l)
                for s, v in r["runs"].items():
                    if s != "SA-native":
                        for x in v["f"]: upd((K, r["n"], r["seed"]), x)
    for p in glob.glob("paper_data/sweep_K*.jsonl"):
        for l in open(p):
            r = json.loads(l)
            for s, v in r["runs"].items():
                if s == "SA-native": continue
                for c in v.values():
                    for x in c["f"]: upd((r["K"], r["n"], r["seed"]), x)
    for p in glob.glob("paper_data/longrun_K*.jsonl"):
        for l in open(p):
            r = json.loads(l)
            for c in r["runs"]["bSB"].values():
                for x in c: upd((r["K"], r["n"], r["seed"]), x["f_sel"])
    for p in glob.glob("paper_data/tuned_K*.jsonl"):
        for l in open(p):
            r = json.loads(l)
            for c in r["runs"].values():
                for x in c: upd((r["K"], r["n"], r["seed"]), x["f_sel"])
    for p in glob.glob("paper_data/results_K*.jsonl"):
        for l in open(p):
            r = json.loads(l); K = r.get("K", int(p.split("_K")[1][0]))
            for s, v in r["res"].items():
                if s != "SA-native": upd((K, r["n"], r["seed"]), v["E"])
    return ref


def sdp_bounds(rigorous=True):
    """SDP lower bounds on f*; with rigorous=True the dual-corrected bounds of bench_sdp_rig.py replace the SCS objectives."""
    lb = {}
    for f in glob.glob("results/sdp_K*.jsonl"):
        for l in open(f):
            r = json.loads(l); lb[(r["K"], r["n"], r["seed"])] = r["lb_f"]
    if rigorous and os.path.exists("results/sdp_rig.jsonl"):
        for l in open("results/sdp_rig.jsonl"):
            r = json.loads(l); lb[(r["K"], r["n"], r["seed"])] = r["lb_f_rig"]
    return lb


def reuse_data():
    """(K, n, seed) -> (reuse-pattern cost, best of reuse + native-SA polish)"""
    out = {}
    if os.path.exists("results/reuse_polish.jsonl"):
        for l in open("results/reuse_polish.jsonl"):
            r = json.loads(l); out[(r["K"], r["n"], r["seed"])] = (r["reuse"], r["best"])
    return out


HAS_PAPER_REF = os.path.exists("paper_data/final_reference.json")


def provenance(K, n, seed, tol=1e-7):
    """Which solver families attain the combined best-known cost of an instance."""
    target = build_ref()[(K, n, seed)]; src = set()
    hit = lambda v: abs(v - target) <= tol * max(1.0, abs(target))
    def scan(path, fn):
        for p in glob.glob(path):
            for l in open(p):
                r = json.loads(l)
                if (r.get("K", K), r["n"], r["seed"]) != (K, n, seed): continue
                fn(r)
    def hyb(r):
        if any(hit(x["f"]) for v in r["runs"]["SA"].values() for x in v): src.add("native SA")
        if any(hit(p["f"]) for v in r["runs"]["bSB"].values() for x in v for p in x["pol"].values()): src.add("bSB with refinement")
        if any(hit(x["f_sel"]) or hit(x["f_desc"]) for v in r["runs"]["bSB"].values() for x in v): src.add("bSB")
    scan(f"results/hybrid_K{K}.jsonl", hyb)
    scan(f"results/hybext_K{K}.jsonl", lambda r: src.add("bSB with refinement") if any(hit(p["f"]) for x in r["runs"] for p in x["pol"].values()) else None)
    scan(f"results/satuned_K{K}.jsonl", lambda r: src.add("native SA") if any(hit(x["f"]) for v in r["runs"].values() for x in v) else None)
    scan("results/core_pol_t0.jsonl", lambda r: src.add("bSB with refinement") if any(hit(x["f"]) for v in r["runs"].values() for x in v) else None)
    scan(f"paper_data/refrun_K{K}.jsonl", lambda r: src.add("native SA") if hit(r["new_best"]) else None)
    def lr(r):
        if any(hit(x["f_sel"]) for v in r["runs"]["SA-native"].values() for x in v): src.add("native SA")
        if any(hit(x["f_sel"]) for v in r["runs"]["bSB"].values() for x in v): src.add("bSB")
    scan(f"paper_data/longrun_K{K}.jsonl", lr)
    def sw(r):
        for s, v in r["runs"].items():
            for c in v.values():
                if any(hit(x) for x in c["f"]): src.add({"SA-native": "native SA", "Tabu": "tabu search"}.get(s, s))
    scan(f"paper_data/sweep_K{K}.jsonl", sw)
    def tt(r):
        for s, v in r["runs"].items():
            if any(hit(x) for x in v["f"]): src.add({"SA-native": "native SA", "Tabu": "tabu search"}.get(s, s))
    scan(f"paper_data/tts_K{K}.jsonl", tt)
    def rk(r):
        for s, v in r["res"].items():
            if hit(v["E"]): src.add({"SA-native": "native SA", "Tabu": "tabu search"}.get(s, s))
    scan(f"paper_data/results_K{K}.jsonl", rk)
    return src
