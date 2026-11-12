"""The two shared functions of runs/stage3-flat-k88-n60, each with its mode.

    python scripts/plot_full_geometry.py

Both decompositions are exact and additive, which is why the components can be
drawn on the same axes as the total:

    ln ds/dp = x_s_tot + q_0(p) + x_p q_1(p)
    psi(s)   = -pi/2   + Psi_0(s) + x_psi Psi_1(s)

stan/full_geometry.stan evaluates every one of them with hier_full.stan's own
basis, gauge and quadrature.  Python reads the draws and draws the axes.

The gauge normalises the modes, ||a|| = 1 and c' Gs c = 1, so q_1 and Psi_1 are
shapes and the whole amplitude sits in x_p and x_psi.  Curves are drawn at the
population mean amplitude; the light envelope is +-1 object sd.
"""
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import cmdstanpy

RUN = "runs/stage3-flat-k88-n60"
OUT = "docs/tex/flat-geometry.png"
C0, C1 = "#004488", "#BB5566"        # mean function, departure mode
PLO, PHI = -15.0, 40.0               # the rest-frame phases epochs actually reach


def draws():
    m = cmdstanpy.CmdStanModel(exe_file="stan/full_geometry",
                               stan_file="stan/full_geometry.stan")
    gq = m.generate_quantities(data=f"{RUN}/data.json",
                               previous_fit=sorted(glob.glob(f"{RUN}/hier_flat-*.csv")))
    d = gq.draws(concat_chains=True)
    idx = {n: j for j, n in enumerate(gq.column_names)}

    def vec(k, n):
        return d[:, [idx[f"{k}[{i}]"] if f"{k}[{i}]" in idx else idx[f"{k}.{i}"]
                     for i in range(1, n + 1)]]
    return d, idx, vec


def med(x):
    return np.median(x, 0)


def main():
    d, idx, vec = draws()
    ngp, ngu = 161, 241
    p = med(vec("p_out", ngp))
    u = med(vec("u_out", ngu))
    s = med(vec("s_out", ngu))
    q0, q1 = vec("q0_p", ngp), vec("q1_p", ngp)
    P0f, P1f = vec("Psi0_u", ngu), vec("Psi1_u", ngu)
    s_ref, s_lo, s_hi = vec("s_ref", ngp), vec("s_qlo", ngp), vec("s_qhi", ngp)
    u_traj = vec("u_traj", ngp)
    xp, xpsi = d[:, idx["x_p_bar"]], d[:, idx["x_psi_bar"]]
    xps, xpss = d[:, idx["x_p_sd"]], d[:, idx["x_psi_sd"]]
    xst = d[:, idx["xst"]]

    inwin = (p >= PLO) & (p <= PHI)
    ur = (med(u_traj)[inwin].min(), med(u_traj)[inwin].max())
    seen = (u >= ur[0]) & (u <= ur[1])

    fig, ax = plt.subplots(2, 2, figsize=(11.4, 8.4))

    # ---- psi(s) and its two pieces, against s --------------------------
    a = ax[0, 0]
    tot = -np.pi / 2 + P0f + xpsi[:, None] * P1f
    a.plot(s, med(tot), color="k", lw=2.2, label=r"$\psi = -\pi/2 + \Psi_0 + x_\psi\Psi_1$")
    a.plot(s, med(-np.pi / 2 + P0f), color=C0, lw=1.6,
           label=r"$-\pi/2 + \Psi_0$  (mean turning)")
    a.plot(s, med(xpsi[:, None] * P1f), color=C1, lw=1.6,
           label=r"$x_\psi\Psi_1$  (departure, at the mean $x_\psi$)")
    lo = med(-np.pi / 2 + P0f + (xpsi - xpss)[:, None] * P1f)
    hi = med(-np.pi / 2 + P0f + (xpsi + xpss)[:, None] * P1f)
    a.fill_between(s, lo, hi, color="k", alpha=0.12, lw=0)
    a.axhline(-np.pi / 2, color="0.7", lw=0.8, ls=":")
    a.axvspan(s[seen].min(), s[seen].max(), color="#FFDD88", alpha=0.35, lw=0, zorder=0)
    a.set_xlim(s[seen].min() - 1.5, s[seen].max() + 1.5)
    a.set_xlabel("arclength $s$ (mag)")
    a.set_ylabel(r"$\psi$ (rad)")
    a.set_title("Turning.  Shaded: the arclengths the epochs reach")
    a.legend(fontsize=7.5, frameon=False, loc="best")

    # ---- the same, against u, the coordinate the basis lives in --------
    a = ax[0, 1]
    a.plot(u, med(P0f), color=C0, lw=1.8, label=r"$\Psi_0(u)$")
    a.plot(u, med(xpsi[:, None] * P1f), color=C1, lw=1.8, label=r"$x_\psi\Psi_1(u)$")
    a.plot(u, med(P1f), color=C1, lw=1.0, ls="--",
           label=r"$\Psi_1(u)$, the gauge-normalised shape")
    a.axvspan(ur[0], ur[1], color="#FFDD88", alpha=0.35, lw=0, zorder=0)
    a.axhline(0, color="0.7", lw=0.8)
    a.set_xlabel(r"$u = \mathrm{asinh}(s/\delta)$")
    a.set_ylabel("rad")
    a.set_title("Turning, in its expansion coordinate")
    a.legend(fontsize=7.5, frameon=False)

    # ---- the log rate and its two pieces -------------------------------
    a = ax[1, 0]
    tot = xst[:, None] + q0 + xp[:, None] * q1
    a.plot(p, med(tot), color="k", lw=2.2, label=r"$\ln ds/dp = x_s^{tot} + q_0 + x_p q_1$")
    a.plot(p, med(xst[:, None] + q0), color=C0, lw=1.6, label=r"$x_s^{tot} + q_0$")
    a.plot(p, med(xp[:, None] * q1), color=C1, lw=1.6,
           label=r"$x_p q_1$  (departure, at the mean $x_p$)")
    lo = med(xst[:, None] + q0 + (xp - xps)[:, None] * q1)
    hi = med(xst[:, None] + q0 + (xp + xps)[:, None] * q1)
    a.fill_between(p, lo, hi, color="k", alpha=0.12, lw=0)
    a.axhline(0, color="0.7", lw=0.8)
    a.axvspan(PLO, PHI, color="#FFDD88", alpha=0.35, lw=0, zorder=0)
    a.set_ylim(np.floor(med(tot)[inwin].min()) - 2.0,
               np.ceil(med(xp[:, None] * q1)[inwin].max()) + 0.5)
    a.set_xlabel("rest-frame phase $p$ (d)")
    a.set_ylabel(r"$\ln\,ds/dp$")
    a.set_title("Traversal rate, where the components add")
    a.legend(fontsize=7.5, frameon=False)

    # ---- s(p) itself ----------------------------------------------------
    a = ax[1, 1]
    a.plot(p, med(s_ref), color="k", lw=2.2, label=r"$s(p)$, mean object")
    a.plot(p, med(s_hi), color=C1, lw=1.3, ls="-", label=r"$x_p + 1$ object sd")
    a.plot(p, med(s_lo), color=C1, lw=1.3, ls="--", label=r"$x_p - 1$ object sd")
    a.fill_between(p, np.minimum(med(s_lo), med(s_hi)),
                   np.maximum(med(s_lo), med(s_hi)), color=C1, alpha=0.15, lw=0)
    a.axhline(0, color="0.7", lw=0.8)
    a.axvspan(PLO, PHI, color="#FFDD88", alpha=0.35, lw=0, zorder=0)
    # Outside the sampled window s plunges past -30 mag, which is the clamp
    # region and carries no epoch; hold the scale to the arclengths reached.
    lim = max(abs(med(s_ref)[inwin]).max(), abs(med(s_hi)[inwin]).max(),
              abs(med(s_lo)[inwin]).max())
    a.set_ylim(-1.35 * lim, 1.35 * lim)
    a.set_xlabel("rest-frame phase $p$ (d)")
    a.set_ylabel("arclength $s$ (mag)")
    a.set_title(r"Traversal.  $s(0)=0$ pins the anchor")
    a.legend(fontsize=7.5, frameon=False)

    fig.tight_layout()
    fig.savefig(OUT, dpi=160)
    print("wrote", OUT)

    w = inwin
    print(f"  arclengths reached over [{PLO:+.0f},{PHI:+.0f}] d: "
          f"s in [{med(s_ref)[w].min():+.3f}, {med(s_ref)[w].max():+.3f}] mag, "
          f"u in [{ur[0]:+.2f}, {ur[1]:+.2f}]")
    print(f"  x_s_tot = {xst.mean():+.4f};  mean |q_0| over the window "
          f"{np.abs(med(q0)[w]).mean():.4f}")
    print(f"  departure vs mean, in the log rate: "
          f"rms(x_p q_1) {np.sqrt((med(xp[:,None]*q1)[w]**2).mean()):.4f} against "
          f"rms(q_0) {np.sqrt((med(q0)[w]**2).mean()):.4f}")
    sv = seen
    print(f"  departure vs mean, in the turning: "
          f"rms(x_psi Psi_1) {np.sqrt((med(xpsi[:,None]*P1f)[sv]**2).mean()):.4f} rad against "
          f"rms(Psi_0) {np.sqrt((med(P0f)[sv]**2).mean()):.4f} rad")
    print(f"  total turning swing over the reached arclengths: "
          f"{med(-np.pi/2 + P0f + xpsi[:,None]*P1f)[sv].max() - med(-np.pi/2 + P0f + xpsi[:,None]*P1f)[sv].min():.3f} rad")


if __name__ == "__main__":
    main()
