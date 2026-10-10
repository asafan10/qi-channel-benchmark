"""Extension of Experiment A: bSB x1 (6 batches) + longer native polish (16k and 64k sweeps), so the polish grid is not the limit."""
import sys, time, numpy as np
from exp_common import *
POL=[16000,64000]
def run(n,seed,K):
    W, csr, a, J, lam, jj = setup(n, seed, K); rng=np.random.default_rng(52000+seed); recs=[]
    for b in range(6):
        t0=time.perf_counter(); Xb=sb_batch(a,J,CFG["sb_steps"],CFG["sb_R"],rng,False)
        es=[qubo_energy(x,a,*jj) for x in Xb]; ch,_=repair(Xb[int(np.argmin(es))].astype(np.int64),K,*csr); fsel=interference(ch,*csr); tsb=time.perf_counter()-t0
        r=dict(f_sel=float(fsel),t=tsb,pol={})
        for s in POL:
            t0=time.perf_counter(); v=sa_native(ch,K,*csr,0.3,0.02,s,810000+1000*b+s%997)[1]; r["pol"][str(s)]=dict(f=float(v),t=time.perf_counter()-t0)
        recs.append(r)
    return recs
if __name__=="__main__":
    SEEDS = range(*map(int, sys.argv[4].split("-"))) if len(sys.argv) > 4 else range(5)
    K=int(sys.argv[1]); sizes=[int(x) for x in sys.argv[2].split(",")]; outf=sys.argv[3]; warm(K); done=done_keys(outf)
    for n in sizes:
        for seed in SEEDS:
            if (n,seed,K) in done: continue
            t0=time.perf_counter(); recs=run(n,seed,K)
            append_jsonl(outf,dict(n=n,seed=seed,K=K,best=min(min(p["f"] for p in r["pol"].values()) for r in recs),runs=recs))
            print(K,n,seed,"%.0fs"%(time.perf_counter()-t0),flush=True)
