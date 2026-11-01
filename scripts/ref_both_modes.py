"""The (24, 16) reference after the softening change, started in both modes.

    python scripts/ref_both_modes.py --check-only   transfer both starts, print lp
    python scripts/ref_both_modes.py                four chains, four jobs, detached
    python scripts/ref_both_modes.py --report       diagnostics of the finished run

Why both modes.  runs/stage3-common-halfnorm straddles two solutions 63 nats
apart -- chains 2 and 4 at tau_psi 2.5 with kappa0 near 1600, chains 1 and 3 at
tau_psi 22 with kappa0 near 9700 -- and runs/stage3-common-split, which supplied
every number in src/dgsn/stan/data.py, never leaves the lower one.  Four chains
from the Salt2 start would most likely land in the lower mode again and the run
would not have looked at the other, so chains 1 and 2 start at the best draw of
the lower mode and chains 3 and 4 at the best draw of the higher.  Whether the
upper chains stay or migrate is the measurement.

What changed since the reference, and nothing else: the phase taper is a weight
on the epoch's term in Eq. (like) rather than a factor on the predicted flux,
and the phase clamp is at eps = 1e-3 rather than 1e-2, so the prediction is
modified for no epoch that exists.  Same program, same orders, same objects,
same seed, same iterations, same treedepth, R the identity throughout.

The reference itself was run at eps = 1e-2 and so is not a controlled
comparison for the clamp; it is still the comparison for everything else, and
the arclength reach is reported here because the two modes differ in it by more
than a magnitude.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from dgsn.data.sample import analysis_sample, median_x1, stage3_sample
from dgsn.salt.globals import fit_globals
from dgsn.stan import data as D
from dgsn.stan.hier_data import build, inits_from_run
from dgsn.stan.run import sample_four

PROG = "hier_common"
K_PSI, K_Q = 24, 16
LOWER = Path("runs/stage3-common-split")       # mode A: tau_psi 2.2-2.7
UPPER = Path("runs/stage3-common-halfnorm")    # straddles; its best draw is mode B
N_EPOCHS = 1568                                # the reference's epoch count, asserted


def starts(names, sample, G, x1_ref, data, verbose=True):
    """Two transferred points, lower mode first.

    Carried through inits_from_run, which reads the relabeled blocks and never
    w_coef -- a rotation mixes w_coef's slots, and UPPER was run rotated, so
    assembling from it would invert nothing visibly and cost tens of thousands
    of nats.  It prints lp at the transferred point, which is the only check
    that catches a bad transfer; a count of blocks carried cannot.
    """
    out = []
    for src in (LOWER, UPPER):
        ini = inits_from_run(src, names, K_PSI, K_Q, x1_ref, mode=False, trav=False,
                             sample=sample, globals_=G, R=None, program=PROG,
                             data=data, verbose=verbose)
        out.append(ini)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-only", action="store_true")
    ap.add_argument("--report", nargs="?", const="", default=None)
    ap.add_argument("--out", default="runs/ref-k2416-eps3")
    ap.add_argument("--iter", type=int, default=500)
    ap.add_argument("--max-treedepth", type=int, default=12)
    args = ap.parse_args()

    sample = analysis_sample()
    names = stage3_sample()
    x1_ref = median_x1()
    G = fit_globals(K_PSI, K_Q, x1_ref)
    data = build(names, K_PSI, K_Q, sample=sample, mode=False, trav=False)

    assert data["eps"] == D.EPS == 1e-3, f"phase clamp is {data['eps']}, not 1e-3"
    assert data["N"] == N_EPOCHS, f"{data['N']} epochs, not the reference's {N_EPOCHS}"
    print(f"{PROG}: {len(names)} objects, {data['N']} epochs, K = ({K_PSI}, {K_Q})")
    print(f"phase clamp eps {data['eps']:g}, taper eps_t {data['eps_t']:g}, "
          f"arclength clamp eps_s {data['eps_s']:g}")

    if args.report is not None:
        from cmdstanpy import from_csv
        from dgsn.stan.run import chain_csvs
        sys.path.insert(0, "scripts")
        from stage3_hier import report
        fit = from_csv([str(c) for c in chain_csvs(args.report or args.out)])
        report(fit, names, sample, data)
        return

    lower, upper = starts(names, sample, G, x1_ref, data)
    # R is the identity here, so raw = w_coef and the turning block is its tail;
    # ||d|| is what the two modes differ in by an order of magnitude.
    for lab, ini in (("lower", lower), ("upper", upper)):
        w = np.asarray(ini["w_coef"], dtype=float)
        print(f"  {lab}: tau_psi {float(ini['tau_psi']):7.3f}  "
              f"|b| {np.linalg.norm(w[:K_Q - 1]):7.2f}  "
              f"|d| {np.linalg.norm(w[K_Q - 1:]):8.2f}  "
              f"b_0 {float(ini['b_0']):+.3f}")

    if args.check_only:
        return

    out = Path(args.out)
    # save_warmup because the question is whether the upper chains stay, and a
    # chain that is in the lower mode by the end of warmup has to be told from
    # one the adaptation drove there.  It changes no transition.
    pids = sample_four(PROG, data, [lower, lower, upper, upper], out, seed=1,
                       threads_per_chain=1, iter_warmup=args.iter,
                       iter_sampling=args.iter, max_treedepth=args.max_treedepth,
                       refresh=10, save_warmup=True, detach=True)
    print(f"launched {out}, chains 1-2 from the lower mode and 3-4 from the upper, "
          f"pids {pids}")


if __name__ == "__main__":
    main()
