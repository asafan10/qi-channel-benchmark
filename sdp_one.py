import sys, json
from bench_sdp import sdp_bound, make_instance
K, n, s, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
r = sdp_bound(make_instance(n, s, K), K)
open(out, "a").write(json.dumps(dict(n=n, seed=s, K=K, **r)) + "\n"); print("done", r["lb_f"])
