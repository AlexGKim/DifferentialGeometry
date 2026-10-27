"""Psi_0(u) for the two tau_psi solutions of a split run.

    python scripts/plot_turning_modes.py --out runs/stage3-common-halfnorm

The turning function is expanded in Legendre polynomials of
u = asinh(s/delta), so u is the coordinate the coefficients live in and the
one that shows what the basis is doing.  forward.stan reports psi at
requested arclengths and psi = -pi/2 + Psi_0 by Eq. (cturn), so the curves
here are its psi shifted by pi/2 --- the same convention as
scripts/plot_traversal.py.

Three panels: each mode against the SALT2 x1 = 0 turning function on its own
y-scale, then both on one scale, because the two differ by about six times in
range and a shared axis flattens the smaller one to a line.  The SALT2 track
is an illustration of what the basis has to represent and is never fitted
here, as Section 5 requires.

The rug along the bottom is where the fit's 1568 epochs actually land in u,
got by asking forward.stan for s at each epoch's rest-frame phase.  It is the
only honest way to see which part of Psi_0 the data speak to.
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))

from dgsn.data.epochs import likelihood_epochs  # noqa: E402
from dgsn.data.sample import analysis_sample, stage3_sample  # noqa: E402
from plot_fit import GRID_KEYS, evaluate, read_draws  # noqa: E402
from plot_modes import split_chains  # noqa: E402

LOW, HIGH = "#1f77b4", "#d62728"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/stage3-common-halfnorm")
    ap.add_argument("--n-draws", type=int, default=10)
    ap.add_argument("--nu", type=int, default=1600)
    args = ap.parse_args()

    out = Path(args.out)
    ga = json.load(open(out / "data.json"))
    grid = {k: ga[k] for k in GRID_KEYS}
    grid["xg"] = ga.get("xg", list(np.polynomial.legendre.leggauss(int(grid["J"]))[0]))
    grid["wg"] = ga.get("wg", list(np.polynomial.legendre.leggauss(int(grid["J"]))[1]))
    kq, kpsi = int(grid["Kq"]), int(grid["Kpsi"])
    dl = float(grid["delta"])

    d = read_draws(out)
    lo_ch, hi_ch, med = split_chains(d)
    print("per-chain median tau_psi: " + ", ".join(f"ch{c} {med[c]:.2f}" for c in sorted(med)))
    print(f"low {lo_ch}  high {hi_ch}")
    groups = [("low", lo_ch, LOW), ("high", hi_ch, HIGH)]

    # the u grid: uniform in u, which is the grid the basis is smooth on.  A
    # grid uniform in s has du ~ 1 at its first step and misses the structure
    # near the anchor entirely.
    ulo, uhi = np.arcsinh(grid["gsmin"] / dl), np.arcsinh(grid["gsmax"] / dl)
    U = np.linspace(ulo, uhi, args.nu)
    s_eval = dl * np.sinh(U)
    p_eval = np.linspace(grid["pmin"] + 0.5, grid["pmax"] - 0.5, 8)

    sample = analysis_sample()
    names = stage3_sample()
    rng = np.random.default_rng(0)

    curves = {}
    for label, chs, _ in groups:
        sub = d[d.chain.isin(chs)].reset_index(drop=True)
        pick = rng.choice(len(sub), size=min(args.n_draws, len(sub)), replace=False)
        curves[label] = [np.asarray(evaluate(grid, sub.iloc[i], 1, 0.03, kpsi, kq,
                                             p_eval, s_eval)["psi"]) + np.pi / 2
                         for i in pick]
        tv = [np.abs(np.diff(c)).sum() for c in curves[label]]
        rg = [c.max() - c.min() for c in curves[label]]
        print(f"  {label:4s}: Psi_0 range {np.mean(rg):6.2f} rad, "
              f"total variation {np.mean(tv):7.2f} rad ({np.mean(tv)/(2*np.pi):.2f} turns)")

    # where the epochs land in u, at one draw of the low mode
    sub = d[d.chain.isin(lo_ch)].reset_index(drop=True)
    draw = sub.iloc[0]
    u_ep = []
    for i, nm in enumerate(names, start=1):
        ep = likelihood_epochs(nm, sample.loc[nm].t0, sample.loc[nm].redshift)
        ph = (ep["dt"] - float(draw[f"delta_t.{i}"])) / (1 + ep["z"])
        sp = np.asarray(evaluate(grid, draw, i, ep["z"], kpsi, kq, ph, s_eval[:8])["s_of_p"])
        u_ep.append(np.arcsinh(sp / dl))
    u_ep = np.concatenate(u_ep)
    print(f"\n{len(u_ep)} epochs land in u in [{u_ep.min():+.2f}, {u_ep.max():+.2f}]; "
          f"middle 90% within [{np.percentile(u_ep,5):+.2f}, {np.percentile(u_ep,95):+.2f}]")

    # the SALT2 x1 = 0 turning function, for scale only
    try:
        from plot_traversal import turning_on
        psi_t = turning_on(U, kpsi)
        print(f"SALT2 x1=0 Psi_0 range {psi_t.max()-psi_t.min():.2f} rad, "
              f"total variation {np.abs(np.diff(psi_t)).sum():.2f} rad")
    except Exception as e:            # the reference is illustration, not a requirement
        psi_t = None
        print(f"(SALT2 reference unavailable: {type(e).__name__}: {e})")

    fig, ax = plt.subplots(1, 3, figsize=(18, 4.6))
    panels = [("low", ax[0], LOW), ("high", ax[1], HIGH)]
    for label, a, colour in panels:
        for j, c in enumerate(curves[label]):
            a.plot(U, c, "-", lw=0.9, alpha=0.6, color=colour,
                   label=f"$\\tau_\\psi$ {label}" if j == 0 else None)
        a.set_title(f"$\\tau_\\psi$ {label}: "
                    f"{np.mean([x.max()-x.min() for x in curves[label]]):.1f} rad range, "
                    f"{np.mean([np.abs(np.diff(x)).sum() for x in curves[label]])/(2*np.pi):.1f} turns",
                    fontsize=10)
    for j, c in enumerate(curves["low"]):
        ax[2].plot(U, c, "-", lw=0.9, alpha=0.6, color=LOW, label="low" if j == 0 else None)
    for j, c in enumerate(curves["high"]):
        ax[2].plot(U, c, "-", lw=0.9, alpha=0.6, color=HIGH, label="high" if j == 0 else None)
    ax[2].set_title("both, one scale", fontsize=10)

    for a in ax:
        if psi_t is not None:
            a.plot(U, psi_t, lw=2.0, color="0.3", zorder=1, label="SALT2 $x_1=0$")
        a.axvline(0, color="0.8", lw=0.6)
        a.axhline(0, color="0.8", lw=0.6)
        lo_, hi_ = a.get_ylim()
        a.plot(u_ep, np.full_like(u_ep, lo_ + 0.02 * (hi_ - lo_)), "|", ms=4,
               color="0.45", alpha=0.25, zorder=0)
        a.set_ylim(lo_, hi_)
        a.set_xlabel(r"$u=\mathrm{asinh}(s/\delta)$,  $\delta=%.3f$" % dl)
        a.set_ylabel(r"$\Psi_0$ (rad)")
        a.legend(fontsize=8, frameon=False)
        sec = a.secondary_xaxis("top", functions=(lambda u: u, lambda u: u))
        ticks = [-1, -0.1, -0.01, 0, 0.01, 0.1, 1]
        sec.set_xticks([np.arcsinh(t / dl) for t in ticks])
        sec.set_xticklabels([f"{t:g}" for t in ticks], fontsize=7)
        sec.set_xlabel("arclength $s$ at those $u$", fontsize=8)
    fig.suptitle(r"$\Psi_0(u)$ for the two $\tau_\psi$ solutions"
                 "   (rug: where the fit's 1568 epochs land in $u$)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path = out / "modes_turning.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
