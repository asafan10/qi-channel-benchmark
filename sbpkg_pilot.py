# Held-out tuning of the simulated-bifurcation package (seeds 100,101; K=3, n=200); 128 agents; Ising form.
import time, numpy as np, torch, simulated_bifurcation as sb, warnings, json, itertools
warnings.filterwarnings("ignore")
from exp_common import *
torch.set_num_threads(1)
def prob(n,seed,K):
    W, csr, a, J, lam, jj = setup(n, seed, K)
    Jm=(J*0.25).toarray(); h=a/2+np.asarray(J.sum(1)).ravel()/4
    return csr, torch.tensor(Jm/2,dtype=torch.float32), torch.tensor(h,dtype=torch.float32)
P={s:prob(200,s,3) for s in (100,101)}
rows=[]
for mode,heated,ts,ps in itertools.product(("ballistic","discrete"),(False,True),(0.1,0.25,0.5),(0.01,0.003,0.001)):
    steps=int(round(1.0/(ts*ps)))
    if steps>20000: continue
    sb.set_env(time_step=ts,pressure_slope=ps)
    fs=[];tt=[]
    for s,(csr,Q,l) in P.items():
        torch.manual_seed(1); t0=time.perf_counter()
        S,E=sb.minimize(Q,l,domain="spin",agents=128,max_steps=steps,mode=mode,heated=heated,verbose=False,best_only=True,early_stopping=False)
        tt.append(time.perf_counter()-t0); x=((S.cpu().numpy()+1)//2).astype(np.int64); fs.append(decode_eval(x,3,csr)[0])
    rows.append(dict(mode=mode,heated=heated,ts=ts,ps=ps,steps=steps,f=fs,t=tt))
    print(mode,heated,ts,ps,steps,[round(f,2) for f in fs],"%.1fs"%np.mean(tt),flush=True)
sb.reset_env()
json.dump(rows,open("results/sbpkg_pilot.json","w"))
