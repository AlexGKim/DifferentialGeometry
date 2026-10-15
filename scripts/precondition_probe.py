"""Does the Fisher rotation buy what it predicts?  Measured in minutes.

    python scripts/precondition_probe.py --from runs/stage3-dense

A fit started from the SALT2 template spends its warmup travelling to the
posterior, and travelling is slow whatever metric it carries: the rotation is
built where the posterior is and says nothing about the road there.  So the
question "does the rotation work" is a question about the sampling phase, and
it is answered by starting at the posterior and watching a few hundred draws
rather than by waiting out a full fit.

Each chain starts at the mean of a previous run, in the current
parameterization, with the rotation built at that same point.  What to compare
against, from runs of this model at these orders: a dense metric gave a median
of 2047 leapfrogs per draw, and a diagonal metric without the rotation sat
60 percent of the time at the cap of 4095.
"""

import argparse
import time
from pathlib import Path

import numpy as np

from dgsn.data.sample import analysis_sample, stage3_sample
from dgsn.stan.hier_data import build
from dgsn.stan.precondition import fisher, point_from_run, rotation
from dgsn.stan.run import sample_four

K_PSI, K_Q = 12, 8


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="src", default="runs/stage3-dense")
    ap.add_argument("--out", default="runs/stage3-probe")
    ap.add_argument("--iter", type=int, default=200)
    ap.add_argument("--max-treedepth", type=int, default=12)
    ap.add_argument("--no-rotation", action="store_true", help="the control: same start, R = identity")
    args = ap.parse_args()

    sample = analysis_sample()
    names = stage3_sample()
    data = build(names, K_PSI, K_Q, sample=sample)
    raw, hyper, per_sn = point_from_run(args.src, K_PSI, K_Q, len(names))

    R = None if args.no_rotation else rotation(fisher(data, raw, hyper, per_sn, K_PSI, K_Q))
    data = build(names, K_PSI, K_Q, sample=sample, R=R)

    # the same point, in the coordinates this model samples
    kbq = np.exp(np.log(np.arange(1, K_Q)).mean())
    kbp = np.exp(np.log(np.arange(1, K_PSI)).mean())
    ini = dict(
        w_coef=raw if R is None else np.linalg.solve(R, raw),
        # the source run predates the pivot: its tau is the scale at mode one
        tau_q=hyper["tau_q"] / kbq ** (hyper["alpha_q_m"] + 0.5),
        tau_psi=hyper["tau_psi"] / kbp ** (hyper["alpha_psi_m"] + 0.5),
        alpha_q_m=hyper["alpha_q_m"], alpha_psi_m=hyper["alpha_psi_m"],
        sigma_g=hyper["sigma_g"], sigma_r=hyper["sigma_r"],
        mu_p=float(np.mean(per_sn["x_p"])), mu_psi=float(np.mean(per_sn["x_psi"])),
        mu_delta=float(np.mean(per_sn["delta_t"])),
        s_delta=max(float(np.std(per_sn["delta_t"])), 1e-3),
        m_peak=per_sn["m_peak"], colr=per_sn["colr"], x_s=per_sn["x_s"],
    )
    ini["u_p"] = per_sn["x_p"] - ini["mu_p"]
    ini["u_psi"] = per_sn["x_psi"] - ini["mu_psi"]
    ini["u_del"] = (per_sn["delta_t"] - ini["mu_delta"]) / ini["s_delta"]
    print(f"{'no rotation (control)' if R is None else 'Fisher rotation'}, warm start from {args.src}")

    t = time.time()
    fit = sample_four("hier", data, ini, Path(args.out), seed=11, threads_per_chain=1,
                      iter_warmup=args.iter, iter_sampling=args.iter,
                      max_treedepth=args.max_treedepth, refresh=50)
    d = fit.draws_pd()
    print(f"\n{time.time()-t:.0f} s wall for {args.iter} warmup + {args.iter} sampling, four chains")
    print(f"  leapfrogs per draw: median {d.n_leapfrog__.median():.0f}, mean {d.n_leapfrog__.mean():.0f}")
    print(f"  treedepth median {d.treedepth__.median():.0f}, at the cap {100*(d.treedepth__ >= args.max_treedepth).mean():.0f}%")
    print(f"  divergent {d.divergent__.sum():.0f} of {len(d)}, stepsize {d.stepsize__.iloc[0]:.2e}")
    print(f"  baselines at these orders: dense metric 2047 leapfrogs; diagonal without rotation 60% at 4095")


if __name__ == "__main__":
    main()
