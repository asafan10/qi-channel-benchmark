"""
gpu_bsb.py -- to be run by the authors on a CUDA GPU (not available in this environment).

Ports the paper's bSB (bench_ran_qubo.sb_batch, identical update rule, dt = 0.5, linear pump, inelastic walls,
c0 = 0.5 sqrt(N-1)/||Jm||_F) to PyTorch so that many replicas run in parallel on the GPU, then compares at matched wall-clock
with the CPU native SA of the paper.

Usage:   python gpu_bsb.py 3 200,500 gpu_K3.jsonl         (K, sizes, output)
Output per instance: for R in {8, 128, 1024} replicas and steps in {3000, 12000}: GPU wall-clock per batch (synchronised)
and the decoded, repaired cost of the best-energy replica, plus the best cost over all replicas.
Analyse with the matched-wall-clock estimator of the paper (Eq. 25), treating one batch as one run.
"""
import sys, time, json, math
import numpy as np, torch
from bench_ran_qubo import make_instance, csr_from_W, build_qubo, decode_eval, qubo_energy, CFG

dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def sb_batch_torch(a, J, steps, R, seed, dt=0.5, dtype=torch.float32):
    g = torch.Generator(device=dev).manual_seed(seed)
    N = len(a)
    Jm = (J * 0.25).tocoo()
    Jt = torch.sparse_coo_tensor(np.vstack([Jm.row, Jm.col]), Jm.data, (N, N), dtype=dtype, device=dev).coalesce().to_sparse_csr()
    h = torch.tensor(a / 2 + np.asarray(J.sum(1)).ravel() / 4, dtype=dtype, device=dev)[:, None]
    scale = math.sqrt(float((Jm.data ** 2).sum())) + 1e-12; c0 = 0.5 * math.sqrt(N - 1) / scale
    X = (torch.rand((N, R), generator=g, device=dev, dtype=dtype) - 0.5) * 0.2
    Y = (torch.rand((N, R), generator=g, device=dev, dtype=dtype) - 0.5) * 0.2
    for t in range(steps):
        at = t / steps
        force = -(Jt @ X + h)
        Y += dt * (-(1.0 - at) * X + c0 * force)
        X += dt * Y
        m = X.abs() > 1
        X = torch.where(m, X.sign(), X); Y = torch.where(m, torch.zeros_like(Y), Y)
    return (X > 0).to(torch.int64).T.cpu().numpy()


if __name__ == "__main__":
    K = int(sys.argv[1]); sizes = [int(x) for x in sys.argv[2].split(",")]; outf = sys.argv[3]
    print("device:", dev, torch.cuda.get_device_name(0) if dev.type == "cuda" else "")
    for n in sizes:
        for seed in range(5):
            W = make_instance(n, seed, K); csr = csr_from_W(W); a, J, lam = build_qubo(W, K, CFG["lam_scale"])
            jj = (J.indptr.astype(np.int64), J.indices.astype(np.int64), J.data.astype(np.float64))
            sb_batch_torch(a, J, 50, 8, 0)  # warm-up
            rec = []
            for R in (8, 128, 1024):
                for steps in (3000, 12000):
                    for k in range(3):
                        if dev.type == "cuda": torch.cuda.synchronize()
                        t0 = time.perf_counter(); Xb = sb_batch_torch(a, J, steps, R, 1000 * seed + k)
                        if dev.type == "cuda": torch.cuda.synchronize()
                        t = time.perf_counter() - t0
                        es = [qubo_energy(x, a, *jj) for x in Xb]
                        f_sel = float(decode_eval(Xb[int(np.argmin(es))], K, csr)[0])
                        f_all = min(float(decode_eval(x, K, csr)[0]) for x in Xb[: min(R, 256)])
                        rec.append(dict(R=R, steps=steps, t=t, f_sel=f_sel, f_best_any=f_all))
                        print(n, seed, R, steps, "%.3fs" % t, round(f_sel, 2), flush=True)
            with open(outf, "a") as fh: fh.write(json.dumps(dict(n=n, seed=seed, K=K, runs=rec)) + "\n")
