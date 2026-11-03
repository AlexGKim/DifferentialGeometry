"""Floor check on <Delta> = 0, and on the two flat placement priors.

<Delta> = 0 is not a gauge fixing, and this file measures what it costs
rather than testing whether it is free.  It was written the other way round,
on the premise that a common shift of every Delta is a shift of the phase
origin which q_0 reshapes to absorb.  That premise is false by construction.
Eq. (travlog) integrates the rate from zero, so s(0) = 0 for every q_0 and
every (x_s, x_p); Eq. (quadrature) integrates the tangent from zero, so
gamma(0) = (0, 0); and Eq. (curvlin) pins Psi_0(0) = Psi_1(0) = 0, so
psi(0) = -pi/2 exactly, which is the statement that p = 0 is the critical
point of M_g.  Absorbing a common shift a would need s(0) = s(a), and the
magnitude translation that comes with it is absorbable -- (1,1) and (0,-1)
span the plane, so m_peak and colr take it -- but psi(0) = -pi/2 is not:
model-definition.tex's gauge enumeration says so in words, the turning
admitting no additive constant, Psi_0(0) = Psi_1(0) = 0 having already pinned
it.  So the phase origin is pinned to the model's own g-band extremum, and
<Delta> is the physical offset from t_0 there: model-definition.tex's prior
section, "t_0 is B maximum and tmax is g maximum, and the offset between them
is measured by the fit, centre and sign included, rather than assumed."

Nor is there anywhere else for the common mode to go.  d_dev is a
sum_to_zero_vector, so with free_mean = 0 the sum of delta_t is exactly zero
and no coordinate in the program carries a common phase offset.  q_0 cannot:
its degree-zero term enters x_s_tot, the log rate, and a rate shift is a
dilation of the phase axis, not a translation of it.  Six hier_common runs at
thirty objects put <Delta> at -0.52 to -0.63 d with a posterior sd of 0.029 to
0.041, and runs/stage4-precon-k2416 at 442 puts it at -0.44 with a sd of
0.002 to 0.006.  The constraint is therefore a displacement of the phase
origin by about half a day, held against a posterior that locates it to a few
hundredths, and the number below is its price in fit.

Likelihood and not log density: the arms carry different numbers of prior
terms, so only sum(lp_sn), which is prior-free and reported by both in
generated quantities, is comparable between them.

Both arms are hier_common and differ only in free_mean, so the comparison is
of <Delta> and of nothing else.  They were two programs until 2026-11-02 --
hier_common_free against hier_common -- which stopped testing the constraint
when hier_common gained a mu_delta of its own: both arms then sampled a free
mean under a label saying one did not, and the 1.96 nats this file is credited
with ceased to be reproducible by it while six source files and the note went
on quoting it.  The two programs also disagreed about how the deviations are
parameterized, so even before that the arms differed in more than the mean.

Run from the repository root.  Takes a few minutes an arm at polish=False.
"""
import sys
import time

import numpy as np

sys.path.insert(0, "src")
from dgsn.data.sample import analysis_sample, median_x1, stage3_sample
from dgsn.salt.globals import fit_globals
from dgsn.stan.hier_data import build, for_free_mean, inits
from dgsn.stan.run import optimize

K_PSI, K_Q = 24, 16


def total_loglike(mode: dict) -> float:
    """sum of lp_sn, the per-supernova log-likelihood of generated quantities."""
    return sum(v for k, v in mode.items() if k.startswith("lp_sn["))



def main() -> None:
    sample = analysis_sample()
    names = stage3_sample()
    G = fit_globals(K_PSI, K_Q, median_x1())
    data = build(names, K_PSI, K_Q, sample=sample, mode=False, trav=False)
    ini = inits(names, K_PSI, K_Q, median_x1(), sample=sample, globals_=G,
                mode=False, trav=False)
    n = len(names)
    print(f"{n} supernovae, {data['N']} epochs, K = ({K_PSI}, {K_Q})\n")

    # The arms differ in <Delta> and in nothing else.  They used to be two
    # programs, hier_common_free against hier_common, and that stopped testing
    # the constraint the day hier_common gained its own mu_delta: both arms
    # then sampled a free mean, the label stayed, and the 1.96 nats this file
    # is credited with became unreproducible by it -- while six source files
    # and the note went on quoting the number.  One program and one data flag
    # cannot drift apart that way, and it isolates <Delta> exactly, where the
    # two programs also disagreed about how the deviations are parameterized.
    d_free = build(names, K_PSI, K_Q, sample=sample, mode=False, trav=False,
                   free_mean=True)
    d_zero = build(names, K_PSI, K_Q, sample=sample, mode=False, trav=False,
                   free_mean=False)
    out = {}
    for tag, prog, start, data in (
            ("Delta free      ", "hier_common", for_free_mean(ini, d_free), d_free),
            ("<Delta> = 0     ", "hier_common", for_free_mean(ini, d_zero), d_zero)):
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
    print("This is a price, not a test of a gauge: s(0) = 0, gamma(0) = 0 and\n"
          "psi(0) = -pi/2 are identities, so no q_0 absorbs a common phase shift,\n"
          "and d_dev being sum-to-zero leaves the common mode nowhere else to go.\n"
          "Read it against the posterior, which locates <Delta> near -0.55 d to\n"
          "within a few hundredths over thirty objects.")


if __name__ == "__main__":
    main()
