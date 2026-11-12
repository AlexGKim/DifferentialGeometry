"""What the two departure modes do to a light curve, from runs/stage3-flat-k88-n60.

    python scripts/plot_full_curves.py

Every curve here is evaluated by stan/full_curves.stan, a generated-quantities
program that reuses hier_full.stan's own data, gauge and transformed
parameters and calls log_flux from dgsn_functions.stan.  Python reads the
draws and draws the axes; it evaluates no basis, no quadrature and no clamp.

The object shown is the population mean one: x_s = 0, which is what
sum_to_zero_vector makes the mean, colour at c_bar, m_peak at the sample mean,
delta_t = 0 and z = 0, so the abscissa is rest-frame phase directly.  The two
modes enter log_flux exactly as they do in the fit, through b + x_p a for the
traversal and d + x_psi c for the turning.

Amplitudes are that draw's own x_p and x_psi against that draw's own a and c,
so the sign pair (a, x_p) -> (-a, -x_p) cancels and no relabelling is needed.
The band is the object-to-object spread, mean +- 1 sd of the amplitude across
the sixty supernovae; the posterior uncertainty on the mean curve is far
narrower and is quoted in the caption instead of drawn.
"""
import glob
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import cmdstanpy

RUN = "runs/stage3-flat-k88-n60"
OUT = "docs/tex/flat-departure-modes.png"
# Epochs run -15.0 to +40.0 d rest-frame; the shaded strips carry no data.
G, R = "#117733", "#CC3311"          # g and r, colourblind-safe
NG = 161


def draws():
    """Curves on the phase grid, one per post-warmup draw, evaluated by Stan."""
    m = cmdstanpy.CmdStanModel(exe_file="stan/full_curves",
                               stan_file="stan/full_curves.stan")
    gq = m.generate_quantities(data=f"{RUN}/data.json",
                               previous_fit=sorted(glob.glob(f"{RUN}/hier_flat-*.csv")))
    d = gq.draws(concat_chains=True)
    idx = {n: j for j, n in enumerate(gq.column_names)}
    get = lambda k: d[:, [idx[f"{k}[{i}]"] if f"{k}[{i}]" in idx else idx[f"{k}.{i}"]
                          for i in range(1, 2 * NG + 1)]]
    return {k: get(k) for k in ("mag_ref", "mag_qlo", "mag_qhi",
                               "mag_plo", "mag_phi", "mag_zero")}, d, idx


def band(x):
    """Posterior median and 68 percent interval down the draw axis."""
    return np.median(x, 0), np.quantile(x, 0.16, 0), np.quantile(x, 0.84, 0)


def main():
    C, d, idx = draws()
    p = np.linspace(-20.0, 45.0, NG)
    g, r = slice(0, NG), slice(NG, 2 * NG)

    xp = d[:, idx["x_p_bar"]]
    xpsi = d[:, idx["x_psi_bar"]]
    xps = d[:, idx["x_p_sd"]]
    sg, sr = 0.0352, 0.0500          # the fitted per-band model error of this run
    xpss = d[:, idx["x_psi_sd"]]

    fig, ax = plt.subplots(2, 2, figsize=(11.0, 8.2))

    # --- the light curves themselves -------------------------------------
    a0 = ax[0, 0]
    for sl, col, lab in ((g, G, "$g$"), (r, R, "$r$")):
        b0, lo, hi = band(C["mag_ref"][:, sl])
        a0.plot(p, b0, color=col, lw=2.0, label=f"{lab}, population mean object")
        a0.fill_between(p, lo, hi, color=col, alpha=0.18, lw=0)
        for klo, khi in (("mag_qlo", "mag_qhi"), ("mag_plo", "mag_phi")):
            a0.plot(p, np.median(C[klo][:, sl], 0), color=col, lw=0.8, alpha=0.75)
            a0.plot(p, np.median(C[khi][:, sl], 0), color=col, lw=0.8, alpha=0.75)
    a0.axvspan(-20, -15, color="0.85", lw=0, zorder=0)
    a0.axvspan(+40, +45, color="0.85", lw=0, zorder=0)
    a0.invert_yaxis()
    a0.set_xlabel("rest-frame phase (d)")
    a0.set_ylabel("magnitude")
    a0.set_title("Mean object (heavy), each mode at $\\pm1$ object sd (light)")
    a0.legend(fontsize=8, frameon=False)

    # --- colour ----------------------------------------------------------
    a1 = ax[0, 1]
    cb, lo, hi = band(C["mag_ref"][:, g] - C["mag_ref"][:, r])
    a1.plot(p, cb, color="k", lw=2.0, label="population mean object")
    a1.fill_between(p, lo, hi, color="k", alpha=0.15, lw=0)
    for k, ls, lab in (("mag_qlo", "--", "$x_p \\mp$ 1 sd"), ("mag_qhi", "--", None),
                       ("mag_plo", ":", "$x_\\psi \\mp$ 1 sd"), ("mag_phi", ":", None)):
        a1.plot(p, np.median(C[k][:, g] - C[k][:, r], 0), color="k", lw=1.0,
                ls=ls, alpha=0.8, label=lab)
    a1.axvspan(-20, -15, color="0.85", lw=0, zorder=0)
    a1.axvspan(+40, +45, color="0.85", lw=0, zorder=0)
    a1.set_xlabel("rest-frame phase (d)")
    a1.set_ylabel("$g - r$ (mag)")
    a1.set_title("Colour")
    a1.legend(fontsize=8, frameon=False)

    # --- the two perturbations -------------------------------------------
    panels = ((ax[1, 0], "mag_qlo", "mag_qhi",
               r"traversal departure  $x_p\,q_1$", xp, xps),
              (ax[1, 1], "mag_plo", "mag_phi",
               r"turning departure  $x_\psi\,\Psi_1$", xpsi, xpss))
    for a, klo, khi, title, amp, asd in panels:
        for sl, col, lab in ((g, G, "$g$"), (r, R, "$r$")):
            lo = C[klo][:, sl] - C["mag_ref"][:, sl]
            hi = C[khi][:, sl] - C["mag_ref"][:, sl]
            a.plot(p, np.median(hi, 0), color=col, lw=2.0, label=f"{lab}, $+1$ sd")
            a.plot(p, np.median(lo, 0), color=col, lw=1.4, ls="--",
                   label=f"{lab}, $-1$ sd")
            q = np.quantile(hi, [0.16, 0.84], axis=0)
            a.fill_between(p, q[0], q[1], color=col, alpha=0.18, lw=0)
        a.axhline(0, color="0.6", lw=0.8)
        a.axvspan(-20, -15, color="0.85", lw=0, zorder=0)
        a.axvspan(+40, +45, color="0.85", lw=0, zorder=0)
        a.set_xlabel("rest-frame phase (d)")
        a.set_ylabel(r"$\Delta$ magnitude from the mean object")
        a.set_title(f"{title}\nmean {amp.mean():+.2f}, object sd "
                    f"{asd.mean():.2f}; shaded = 68% posterior on $+1$ sd")
        a.legend(fontsize=8, frameon=False)

    fig.tight_layout()
    fig.savefig(OUT, dpi=160)
    print("wrote", OUT)

    # numbers for the caption and the report
    for key, title in (("mag_qhi", "q_1 at +1 object sd"),
                       ("mag_phi", "Psi_1 at +1 object sd"),
                       ("mag_zero", "turning both modes off entirely")):
        for sl, lab in ((g, "g"), (r, "r")):
            dm = C[key][:, sl] - C["mag_ref"][:, sl]
            md = np.median(dm, 0)
            sig = sg if lab == "g" else sr
            for a, b in ((-15.0, 40.0), (-10.0, 35.0)):
                w = (p >= a) & (p <= b)
                j = int(np.argmax(np.abs(np.where(w, md, 0.0))))
                rms = np.sqrt((md[w] ** 2).mean())
                print(f"  {title:<34} {lab}  [{a:+5.0f},{b:+5.0f}] d: "
                      f"peak |dm| {abs(md[j]):6.3f} at {p[j]:+6.1f} d, "
                      f"rms {rms:6.3f} mag = {rms/sig:6.1f} x sigma_{lab}")
    base = band(C["mag_ref"][:, g])
    print(f"  posterior 68% width of the base g curve at peak: "
          f"{(base[2]-base[1])[np.argmin(base[0])]:.4f} mag")


if __name__ == "__main__":
    main()
