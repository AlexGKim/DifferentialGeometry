"""Floor check on <Delta> = 0, and on the two flat placement priors.

The claim under test is that the constraint costs no likelihood: a common
shift of every Delta is a shift of the phase origin, and q_0 is a free
function of v that reshapes to absorb it.  If that is right, the mode of the
constrained program reaches the same total log-likelihood as the free one on
the same data.  If it is not, the basis is too stiff to absorb the shift and
the constraint is buying identifiability with fit.

Likelihood and not log density: the two programs carry different priors by
construction, so only sum(lp_sn), which is prior-free and reported by both
in generated quantities, is comparable between them.

Run from the repository root.  Takes a few minutes an arm at polish=False.
"""
import sys
import time

import numpy as np

sys.path.insert(0, "src")
from dgsn.data.sample import analysis_sample, median_x1, stage3_sample
from dgsn.salt.globals import fit_globals
from dgsn.stan.hier_data import build, inits
from dgsn.stan.run import optimize

K_PSI, K_Q = 24, 16


def total_loglike(mode: dict) -> float:
    """sum of lp_sn, the per-supernova log-likelihood of generated quantities."""
    return sum(v for k, v in mode.items() if k.startswith("lp_sn["))


def to_free(ini: dict, n: int) -> dict:
    """``inits`` speaks the constrained program's coordinates, which are the
    ladder's; this rewrites one of its starts in the free program's.  The
    point itself is unchanged --- the same colr and the same delta_t reach
    the likelihood either way --- so the two arms start together and the
    comparison is of optima and not of starts."""
    out = {k: v for k, v in ini.items() if k not in ("c_bar", "c_dev", "delta_t")}
    out["colr"] = float(ini["c_bar"]) + np.asarray(ini["c_dev"], dtype=float)
    out["mu_delta"] = 0.0
    out["s_delta"] = 1.0
    out["u_del"] = np.asarray(ini["delta_t"], dtype=float)
    return out


def main() -> None:
    sample = analysis_sample()
    names = stage3_sample()
    G = fit_globals(K_PSI, K_Q, median_x1())
    data = build(names, K_PSI, K_Q, sample=sample, mode=False, trav=False)
    ini = inits(names, K_PSI, K_Q, median_x1(), sample=sample, globals_=G,
                mode=False, trav=False)
    n = len(names)
    print(f"{n} supernovae, {data['N']} epochs, K = ({K_PSI}, {K_Q})\n")

    out = {}
    for tag, prog, start in (("Delta free      ", "hier_common_free", to_free(ini, n)),
                             ("<Delta> = 0     ", "hier_common", ini)):
        t0 = time.time()
        m = optimize(prog, data, start, jacobian=False, polish=False)
        ll = total_loglike(m)
        dbar = (float(m["mu_delta"]) if "mu_delta" in m
                else float(np.mean([m[f"delta_t[{i}]"] for i in range(1, n + 1)])))
        spread = float(np.std([m[f"delta_t[{i}]"] for i in range(1, n + 1)]))
        out[tag] = ll
        # Eq. (like) is now the archive's errors alone, so lp_sn is exactly
        # -chi^2/2 and the fit quality reads straight off it.
        print(f"{tag}  sum lp_sn = {ll:14.4f}   chi2/epoch = {-2*ll/data['N']:7.3f}   "
              f"<Delta> = {dbar:+.4f} d   sd(Delta) = {spread:.4f} d   [{time.time()-t0:.0f} s]")

    d = out["<Delta> = 0     "] - out["Delta free      "]
    print(f"\nconstraint costs {-d:+.4f} nats of log-likelihood over {data['N']} epochs "
          f"({-d / data['N']:+.2e} an epoch)")
    print("A cost of order a nat or less is the constraint being free: q_0 absorbs\n"
          "the shift.  A cost of many nats is the basis failing to, and would mean\n"
          "<Delta> = 0 is not a gauge at this K_q.")


if __name__ == "__main__":
    main()
