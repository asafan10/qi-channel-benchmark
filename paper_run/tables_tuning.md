
## Table 7. bSB re-tuning study (held-out instances: n=200, K=3, seeds 100-102)

Grid: step size {0.25,0.5,1.0} x pump exponent {0.5,1,2} x coupling multiplier {0.5,1,2} x penalty scale rho {0.12,0.25,0.5} = 81 configs; 4x steps, 2 batches of 8 replicas per instance.

- Best / median / worst config (mean gap at 4x steps): 6.05% / 7.51% / 22.52%; untuned default: 7.18%
- Top 5 at 4x: dt0.5_g2.0_c2.0_rho0.25 = 6.05%; dt1.0_g1.0_c0.5_rho0.12 = 6.07%; dt1.0_g0.5_c1.0_rho0.25 = 6.11%; dt0.25_g1.0_c1.0_rho0.12 = 6.29%; dt0.25_g0.5_c1.0_rho0.12 = 6.40%
- Confirmation at 16x steps (top-3 + default): dt0.5_g2.0_c2.0_rho0.25 = 5.71%; dt1.0_g1.0_c0.5_rho0.12 = 6.70%; dt1.0_g0.5_c1.0_rho0.25 = 6.76%; dt0.5_g1.0_c1.0_rho0.25 = 7.14%
- Frozen configuration: {'dt': 0.5, 'gamma': 2.0, 'c0m': 2.0, 'rho': 0.25}
- Transfer to evaluation instances: no improvement (see Table 9).
