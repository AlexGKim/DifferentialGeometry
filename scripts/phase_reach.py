"""How far outside the modelling window the Delta box can carry an epoch.

    python scripts/phase_reach.py

Section 9.3 of docs/tex/implementation.tex rests on these numbers.  The
softening exists for draws that put an epoch outside the nominal phase window,
and whether any draw can is a question about the support rather than about the
posterior: delta_t = mu_delta + d_dev, with mu_delta declared in [-3, 3] and
the deviations walled at the same, so Delta reaches [-6, 6] for one object
while the modelling window exceeds the data window by only five days at each
end.

Two earlier statements of the argument -- the comment at src/dgsn/stan/data.py
and Sections 4 and 9.3 of the notes -- used mu_delta's bound alone, concluded a
reach of [-17.98, 42.99] and so claimed the clamp could never act on an epoch
that exists.  It can.  The realized draws stay well inside, which is a
different and weaker statement.

The epoch set comes from hier_data.build, so it is the fit's own selection and
not a second one assembled here.
"""

import numpy as np

from dgsn.data.sample import analysis_sample, stage3_sample
from dgsn.stan import data as D
from dgsn.stan.hier_data import build

BOX = 3.0        # the declared bound on mu_delta, and the wall on d_dev


def main():
    sample = analysis_sample()
    names = stage3_sample()
    data = build(names, D.K_PSI, D.K_Q, sample=sample, mode=False, trav=False)
    dt, z = data["dt"], data["z"]
    start, stop = data["start"], data["stop"]
    pmin, pmax = D.MODEL_WINDOW
    reach = 2 * BOX

    print(f"modelling window ({pmin:g}, {pmax:g}) d, data window "
          f"{D.DATA_WINDOW[0]:g} to {D.DATA_WINDOW[1]:g}, "
          f"{data['N']} epochs over {data['N_sn']} objects")
    print(f"Delta = mu_delta + d_dev, each bounded by {BOX:g}, so "
          f"|Delta| <= {reach:g}\n")
    print(f"{'sn':>18} {'z':>7} {'phase at Delta=0':>18} "
          f"{'at |Delta|=' + format(reach, 'g'):>19} "
          f"{'exits at |Delta|':>17}")
    worst, who = np.inf, None
    out = []
    for n, name in enumerate(names):
        s = dt[start[n] - 1:stop[n]]
        zz = float(z[n])
        p0 = (s.min() / (1 + zz), s.max() / (1 + zz))
        pr = ((s.min() - reach) / (1 + zz), (s.max() + reach) / (1 + zz))
        early = s.min() - pmin * (1 + zz)     # Delta the earliest epoch tolerates
        late = pmax * (1 + zz) - s.max()      # and the latest
        lim = min(early, late)
        out.append((pr, lim))
        if lim < worst:
            worst, who = lim, name
        flag = "  exits" if (pr[0] < pmin or pr[1] > pmax) else ""
        print(f"{name:>18} {zz:7.4f} {p0[0]:8.2f} {p0[1]:8.2f} "
              f"{pr[0]:9.2f} {pr[1]:9.2f} {lim:17.3f}{flag}")

    lo = min(p[0] for p, _ in out)
    hi = max(p[1] for p, _ in out)
    nlo = sum(1 for p, _ in out if p[0] < pmin)
    nhi = sum(1 for p, _ in out if p[1] > pmax)
    print(f"\nreach on the support             [{lo:.2f}, {hi:.2f}] d")
    print(f"objects that can be carried out  {nlo} early, {nhi} late, "
          f"of {len(names)}")
    print(f"binding object                   {who}, at |Delta| = {worst:.3f} d")
    lims = np.array([lim for _, lim in out])
    print(f"window crossed before the bound binds    "
          f"{(lims < reach).sum()} of {len(names)} objects, "
          f"by {reach - lims.min():.3f} d at the earliest")
    print(f"keeping every epoch inside needs a wall at "
          f"+-{worst - BOX:.2f}; it is at +-{BOX:g}")


if __name__ == "__main__":
    main()
