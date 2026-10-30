"""Which prior on Eq. (coefprior)'s scales gives the better geometry.

The stage-4 fit at N = 442 never sampled: 1200 of 1200 transitions at the
treedepth cap, within-chain split-Rhat 3.5 to 4.8 on lp__, and traces that a
straight line in draw index accounts for to 95-99%.  The funnel is the
candidate -- on runs/stage3-common-split, this ladder's reference,
corr(log tau_psi, log||d||) = +0.826 with ||d|| at Rhat 1.371 and ESS 9.6 --
so the question is whether removing it moves the stepsize.  Three treatments,
one program, differing in nothing but the prior:

    shape 0   tau and alpha sampled; hier_common.stan exactly.  Verified
              bit-identical on lp__ and all 162 gradient coordinates, so the
              reference run below stands in for this arm and is not re-run.
    shape 1   tau and alpha held, normal.  runs/stage3-common-fixtau measured
              ||d|| ESS 18 -> 71 but Rhat 1.004 -> 1.273; this re-measures it
              at the current parameterization as a harness control.
    shape 2   the same held scales, Student-t of nu degrees of freedom.  The
              heavy tail is the point: it removes the four funnel coordinates
              without the quadratic cost of departing from a held scale that
              reverted shape 1.

Why N = 30 and K = (24, 16).  Those are the reference run's, so its posterior
is the answer to compare against -- Section "Traps in measuring" requires the
floor be measured with the exact answer substituted before any of this is
read as a result.  K matches the stage-4 run too, so a verdict transfers.

Cost.  The reference took 49-60 ks a chain at mean treedepth 10-11.9, so each
arm here is of order fifteen hours and the three run together.  Wall time is
therefore contended and not comparable between arms; leapfrogs per effective
draw is, and is what the report prints.

    python scripts/scale_ladder.py --check-only
    python scripts/scale_ladder.py
    python scripts/scale_ladder.py --report
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from dgsn.data.sample import analysis_sample, median_x1, stage3_sample
from dgsn.salt.globals import fit_globals
from dgsn.stan import data as D
from dgsn.stan.hier_data import build, inits
from dgsn.stan.run import chain_csvs, log_prob, sample_four

sys.path.insert(0, str(Path(__file__).parent))
from stage3_hier import gradient_check  # noqa: E402

PROG = "hier_common_scale"
K_PSI, K_Q = 24, 16
REFERENCE = Path("runs/stage3-common-split")      # shape 0, already measured

# name, prior_shape, nu
ARMS = [
    ("fixnorm", 1, 1.0),      # held scales, normal: the fixtau control
    ("cauchy", 2, 1.0),       # held scales, Cauchy
    ("t4", 2, 4.0),           # held scales, Student-t nu = 4
]


def arm_data(names, sample, shape, nu):
    return build(names, K_PSI, K_Q, sample=sample, mode=False, trav=False,
                 prior_shape=shape, nu_coef=nu)


def arm_inits(names, sample, G, x1_ref, shape):
    """``inits`` writes a scalar for each scale; shape 0 declares them as
    length-one arrays and shapes 1 and 2 declare nothing, so the point is
    adapted to the declaration rather than the generator being branched."""
    ini = inits(names, K_PSI, K_Q, x1_ref, sample=sample, globals_=G,
                mode=False, trav=False)
    for k in ("tau_q", "tau_psi", "alpha_q_m", "alpha_psi_m"):
        ini[k] = [float(ini[k])] if shape == 0 else []
    return ini


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-only", action="store_true",
                    help="coordinate counts and the gradient check, no sampling")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--iter", type=int, default=500, help="the reference run's 500")
    ap.add_argument("--max-treedepth", type=int, default=12)
    ap.add_argument("--arms", default=None,
                    help="comma-separated subset of " + ",".join(a[0] for a in ARMS))
    args = ap.parse_args()

    arms = ARMS if args.arms is None else [a for a in ARMS if a[0] in args.arms.split(",")]
    if not arms:
        sys.exit(f"no arm matched {args.arms!r}; have {[a[0] for a in ARMS]}")

    if args.report:
        from scale_report import report
        report(arms, REFERENCE, K_PSI, K_Q)
        return

    names = stage3_sample()
    sample = analysis_sample()
    x1_ref = median_x1()
    G = fit_globals(K_PSI, K_Q, x1_ref)
    print(f"{PROG}: {len(names)} objects, K = ({K_PSI}, {K_Q}), "
          f"{args.iter}/{args.iter}, max treedepth {args.max_treedepth}")
    print(f"held scales: sigma_q {D.TAU_Q} sigma_psi {D.TAU_PSI} "
          f"alpha_q {D.ALPHA_Q} alpha_psi {D.ALPHA_PSI} "
          f"(runs/stage3-common-split; see stan/hier_common_scale.stan on why "
          f"holding them at measured values is legitimate here)")
    print(f"shape 0 is not re-run: {REFERENCE} is that arm, "
          f"verified identical on lp__ and all gradient coordinates\n")

    pids = {}
    for tag, shape, nu in arms:
        data = arm_data(names, sample, shape, nu)
        ini = arm_inits(names, sample, G, x1_ref, shape)
        lp0, worst = gradient_check(data, ini, name=PROG)
        n_coord = len(log_prob(PROG, data, ini)[1])
        out = Path(f"runs/ladder-{tag}")
        print(f"{tag:9s} shape {shape} nu {nu:g}: {n_coord} coordinates, "
              f"lp {lp0:.2f}, gradient worst {worst:.1e}")
        if args.check_only:
            continue
        res = sample_four(PROG, data, ini, out, seed=1, threads_per_chain=1,
                          iter_warmup=args.iter, iter_sampling=args.iter,
                          max_treedepth=args.max_treedepth, refresh=10, detach=True)
        pids[tag] = res
        print(f"{'':9s} launched {out}, pids {res}")

    if args.check_only:
        return
    print("\nprogress:  grep -h Iteration runs/ladder-*/chain-*.log | tail")
    print("when done: python scripts/scale_ladder.py --report")


if __name__ == "__main__":
    main()
