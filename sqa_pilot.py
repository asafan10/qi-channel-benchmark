# SQA hyper-parameters chosen on the paper's pilot instance (n=100, K=3, seed 99) and on seed 100, then frozen.
import itertools, time, json, numpy as np
from exp_common import *
from sqa import sqa_solve
warm(3); K=3; rows=[]
inst={s:setup(100,s,K) for s in (99,100)}
for P,T,G0,sw in itertools.product((8,16),(0.05,0.15,0.4),(1.0,3.0,8.0),(200,1000)):
    fs=[];ts=[]
    for s,(W,csr,a,J,lam,jj) in inst.items():
        rng=np.random.default_rng(1)
        for k in range(2):
            t0=time.perf_counter(); x,e=sqa_solve(a,J,P,T,G0,0.01,sw,rng,k); ts.append(time.perf_counter()-t0); fs.append(decode_eval(x,K,csr)[0])
    rows.append(dict(P=P,T=T,G0=G0,sweeps=sw,f=float(np.mean(fs)),t=float(np.mean(ts))))
    print(P,T,G0,sw,"mean f %.2f  t %.2fs"%(np.mean(fs),np.mean(ts)),flush=True)
json.dump(rows,open("results/sqa_pilot.json","w"))
b=min(rows,key=lambda r:r["f"]); print("BEST",b)
