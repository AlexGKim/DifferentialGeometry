"""How many modes of the global functions the data can actually see.

    python scripts/visible_modes.py --k-psi 12 24 48 --k-q 8 16

The truncation orders set the size of the only block the sampler finds hard:
every stiff direction of the hierarchy lives in the coefficients.  A mode the
data cannot distinguish from zero costs sampling time and contributes a
near-null direction, without improving any fit --- so the useful order is
bounded by how many singular values of the forward map's Jacobian stand above
the noise, and that is measurable in seconds, without running a fit.

The count below is the number of singular values of ``d log f_i / d coef_k``,
weighted by each epoch's own error, that exceed a fraction of the largest.
Section 5.1 sets the orders by a held-out log score over a grid of fits; this
is the cheap bound that says which points of that grid can possibly differ.
"""

import argparse

import numpy as np

from dgsn.data.sample import analysis_sample, median_x1, stage3_sample
from dgsn.salt.globals import fit_globals
from dgsn.stan.hier_data import build, inits
from dgsn.stan.precondition import jacobian


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k-psi", type=int, nargs="+", default=[12, 24])
    ap.add_argument("--k-q", type=int, nargs="+", default=[8, 16])
    ap.add_argument("--n-sn", type=int, default=10)
    ap.add_argument("--cut", type=float, default=0.01, help="fraction of the largest singular value")
    args = ap.parse_args()

    sample = analysis_sample()
    names = stage3_sample()[:args.n_sn]
    x1_ref = median_x1()
    print(f"{len(names)} objects; a mode counts as visible when its singular value exceeds "
          f"{args.cut:g} of the largest\n")
    print(" K_psi  K_q    turning visible   traversal visible   cond(turning)")
    for kp in args.k_psi:
        for kq in args.k_q:
            data = build(names, kp, kq, sample=sample)
            ini = inits(names, kp, kq, x1_ref, sample=sample, globals_=fit_globals(kp, kq, x1_ref))
            raw = np.asarray(ini["w_coef"], dtype=float)
            hyper = {k: ini[k] for k in ("alpha_q_m", "alpha_psi_m", "tau_q", "tau_psi")}
            # pieces() assembles x_s from s_bar and x_s_dev, so this follows the
            # parameterization rather than restating it
            from dgsn.stan.precondition import pieces
            per = pieces(ini, kp, kq)[2]
            J, f = jacobian(data, raw, hyper, per, kp, kq)
            # each epoch in the error the likelihood gives it, Eq. (like),
            # which is now the archive's alone: the fitted scatter that used
            # to enter the weight here is gone, so a bright epoch does count
            # for as much as its small sigma_f says it does
            sig_f = np.asarray(data["sig_f"], dtype=float)
            w = f / sig_f
            Jw = w[:, None] * J
            n1, n3 = kq - 1, kp - 1
            # the mean functions: no gauge normalization, so no null direction.
            # (z_a and z_c carry one each by construction, their radial
            # directions being invisible to a = at/||at||.)
            trav = np.linalg.svd(Jw[:, n1:2 * n1], compute_uv=False)
            turn = np.linalg.svd(Jw[:, 2 * n1 + n3:], compute_uv=False)
            print(f" {kp:5d} {kq:4d}      {int((turn > args.cut*turn[0]).sum()):3d} of {n3:<3d}"
                  f"           {int((trav > args.cut*trav[0]).sum()):3d} of {n1:<3d}"
                  f"          {turn[0]/turn[-1]:.1e}")


if __name__ == "__main__":
    main()
