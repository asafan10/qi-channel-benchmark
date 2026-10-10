"""Experiment D: path-integral SQA (sqa.py; P=16, T=0.15, Gamma 1.0->0.01, 4000 sweeps; chosen on held-out seeds 99,100, n=100)
under the paper's fixed-budget restart protocol P1 (3 s n<=50, 6 s n<=150, 12 s otherwise), K=3, seeds 0-4,
with SA-native, SA-QUBO and bSB re-run on the same machine as anchors (frozen CFG)."""
import sys, time, numpy as np
from exp_common import *
from sqa import sqa_solve
def anytime(step, budget):
    t0=time.perf_counter(); best=np.inf; feas=None; k=0
    while time.perf_counter()-t0<budget:
        e,f=step(k); k+=1
        if e<best: best,feas=e,f
    return dict(E=float(best),feas=float(feas),restarts=k,t=time.perf_counter()-t0)
if __name__=="__main__":
    K=3; sizes=[int(x) for x in sys.argv[1].split(",")]; outf=sys.argv[2]; warm(K); done=done_keys(outf)
    for n in sizes:
        budget=3.0 if n<=50 else (6.0 if n<=150 else 12.0)
        for seed in range(5):
            if (n,seed,K) in done: continue
            W,csr,a,J,lam,jj=setup(n,seed,K); rng=np.random.default_rng(61000+seed); N=n*K
            T0,T1=CFG["T0_scale"]*lam,CFG["T1"]; res={}
            res["SQA"]=anytime(lambda k: decode_eval(sqa_solve(a,J,16,0.15,1.0,0.01,4000,rng,k)[0],K,csr),budget)
            res["SA-native"]=anytime(lambda k:(sa_native(rng.integers(0,K,n),K,*csr,2.0,0.02,CFG["native_sweeps"],k)[1],1.0),budget)
            res["SA-QUBO"]=anytime(lambda k: decode_eval(sa_qubo(rng.integers(0,2,N),a,*jj,T0,T1,CFG["sweeps"],k)[0],K,csr),budget)
            def sbstep(k):
                Xb=sb_batch(a,J,CFG["sb_steps"],CFG["sb_R"],rng,False); es=[qubo_energy(x,a,*jj) for x in Xb]
                return decode_eval(Xb[int(np.argmin(es))],K,csr)
            res["bSB"]=anytime(sbstep,budget)
            append_jsonl(outf,dict(n=n,seed=seed,K=K,budget=budget,res=res))
            print(n,seed,{k:(round(v["E"],2),v["restarts"]) for k,v in res.items()},flush=True)
