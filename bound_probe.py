from ortools.sat.python import cp_model
import time, numpy as np, cvxpy as cp
from exp_common import *
K=3
for n in (50,100):
    W,csr,a,J,lam,jj=setup(n,0,K); Wtot=np.triu(W,1).sum()
    # CP-SAT, 1 worker, 60 s, hint from a long SA run
    ch,best=sa_native(np.random.default_rng(0).integers(0,K,n),K,*csr,2.0,0.02,16000,1)
    m=cp_model.CpModel(); SC=10000; x=[[m.NewBoolVar("") for c in range(K)] for i in range(n)]
    for i in range(n): m.AddExactlyOne(x[i])
    obj=[]; ii,jx=np.nonzero(np.triu(W,1))
    for i,j in zip(ii,jx):
        w=int(round(W[i,j]*SC)); z=m.NewBoolVar("")
        for c in range(K): m.AddBoolOr([x[i][c].Not(),x[j][c].Not(),z])
        obj.append(w*z)
    for i in range(n):
        for c in range(K): m.AddHint(x[i][c],int(ch[i]==c))
    m.Minimize(sum(obj)); s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=60; s.parameters.num_workers=1
    t0=time.perf_counter(); st=s.Solve(m)
    print(f"n={n}: SA best {best:.2f}; CP-SAT obj {s.ObjectiveValue()/SC:.2f} bound {s.BestObjectiveBound()/SC:.2f} status {s.StatusName(st)} {time.perf_counter()-t0:.0f}s",flush=True)
    # SDP (Frieze-Jerrum): max-K-cut <= (K-1)/K * sum w (1 - X_ij); min f >= W - that
    X=cp.Variable((n,n),PSD=True); Wu=np.triu(W,1)
    cons=[cp.diag(X)==1, X>=-1/(K-1)]
    prob=cp.Problem(cp.Maximize((K-1)/K*cp.sum(cp.multiply(Wu,1-X))),cons); t0=time.perf_counter(); prob.solve(solver="SCS",eps=1e-5)
    print(f"   SDP lower bound on f*: {Wtot-prob.value:.2f}  (W={Wtot:.2f})  {time.perf_counter()-t0:.0f}s",flush=True)
