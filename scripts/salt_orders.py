"""Regenerate what Section 2.3 quotes from a script "not in the repository":
the truncation orders a SALT2 track asks for, now with the traversal
expanded in v, and the quadrature accuracy of Section 5.1 remeasured in v.

    python scripts/salt_orders.py [--quick]

Each track is fit by least squares at each pair of orders and the truncated
model evaluated by forward.stan; the error is the maximum over the data
window, in magnitudes, at two grid resolutions.
"""

import sys

import numpy as np

from dgsn.salt.globals import fit_track, refine_track, track_error
from dgsn.salt.tracks import salt_track
from dgsn.stan import data as D

X1 = (-2.0, -1.0, 0.0, 1.0, 1.7)
K_PSI = (12, 24, 48)
K_Q = (8, 16, 32, 64)
quick = "--quick" in sys.argv


def main():
    n = 2000 if quick else 10000
    tracks = {x1: salt_track(x1, n_per_side=n) for x1 in X1}
    for x1, tr in tracks.items():
        print(f"x1={x1:+.1f}: covered [{tr['p_range'][0]:.2f}, {tr['p_range'][1]:.2f}] d, "
              f"s(-15)={np.interp(-15, tr['p'], tr['s']):.2f}, s(+40)={np.interp(40, tr['p'], tr['s']):.2f}, "
              f"mean log rate in v {np.trapezoid(tr['log_rate'], tr['v']) / (tr['v'][-1] - tr['v'][0]):.2f}")

    print("\nTurning alone: max |gamma(s) - track| over the data window (mag), K_q = 64")
    print("   x1  " + "".join(f"  K_psi={k:<3d}" for k in K_PSI))
    fits = {}
    for x1 in X1:
        row = []
        for kp in K_PSI:
            fits[x1, kp] = fit_track(tracks[x1], kp, 64)
            row.append(track_error(tracks[x1], fits[x1, kp]["psi"], fits[x1, kp]["q"],
                                   fits[x1, kp]["x_s"], kp, 64)["curve"])
        print(f"  {x1:+.1f} " + "".join(f"  {e:9.4f} " for e in row))

    print("\nLight curve, fit in magnitude space: max |m_model - m_track| over the data window (mag)")
    print("  the least-squares start refined by track_fit.stan; m = 10 | 20 panels per side")
    print("   x1   K_psi " + "".join(f"        K_q={k:<3d}    " for k in K_Q))
    refined = {}
    for x1 in X1:
        for kp in K_PSI:
            row = []
            for kq in K_Q:
                f = fit_track(tracks[x1], kp, kq)
                r = refine_track(tracks[x1], f, kp, kq)
                refined[x1, kp, kq] = r
                e2 = track_error(tracks[x1], r["psi"], r["q"], r["x_s"], kp, kq, 2 * D.M_PER_SIDE,
                                 placement=(r["m_peak"], r["colr"], r["delta_t"]))
                row.append((r["lightcurve"], e2["lightcurve"]))
            print(f"  {x1:+.1f}   {kp:3d}  " + "".join(f"  {a:8.4f} {b:8.4f}" for a, b in row))

    print("\nTraversal of the refined fit: max |s_model(p) - s_track(p)| over the data window (mag), K_psi = 48")
    print("   x1  " + "".join(f"   K_q={k:<3d}" for k in K_Q))
    for x1 in X1:
        row = []
        for kq in K_Q:
            r = refined[x1, 48, kq]
            row.append(track_error(tracks[x1], r["psi"], r["q"], r["x_s"], 48, kq,
                                   placement=(r["m_peak"], r["colr"], r["delta_t"]))["traversal"])
        print(f"  {x1:+.1f} " + "".join(f"  {e:8.4f} " for e in row))

    print("\nQuadrature: |m - m=40| at K_psi=24, K_q=16, over the data window, x1 = 0")
    f = refined[0.0, 24, 16]
    pl = (f["m_peak"], f["colr"], f["delta_t"])
    ref = track_error(tracks[0.0], f["psi"], f["q"], f["x_s"], 24, 16, 40, placement=pl)
    for m in (5, 10, 20):
        e = track_error(tracks[0.0], f["psi"], f["q"], f["x_s"], 24, 16, m, placement=pl)
        print(f"  m={m:2d}: traversal {np.abs(e['s_of_p'] - ref['s_of_p']).max():.1e} mag, "
              f"curve {np.abs(e['gamma'] - ref['gamma']).max():.1e} mag, "
              f"light curve {np.abs(e['mag'] - ref['mag']).max():.1e} mag")


if __name__ == "__main__":
    main()
