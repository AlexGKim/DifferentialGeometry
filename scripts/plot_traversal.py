"""The two figures of Section 9 of docs/tex/implementation.tex.

    python scripts/plot_traversal.py

Why they exist.  The traversal is expanded in Legendre polynomials of
v = asinh(p/delta_p) rather than of phase, and integrated by a rule whose
nodes are laid in v.  Both choices are hard to motivate in words and
immediate in a picture: the rate is a cusp in phase and a smooth broad curve
in v, and the rule puts a fifth of its nodes in the two panels either side of
zero, which is where the cusp is.

The SALT2 track at x1 = 0 is used as an illustration of what the basis has
to represent, not as anything fitted: the standing commitment is that
templates quantify systematics and never fit, and nothing here reaches a
likelihood.  Stan evaluates the model throughout, through forward.stan.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dgsn.salt.globals import fit_track
from dgsn.salt.tracks import salt_track
from dgsn.stan import data as D
from dgsn.stan.run import fixed_param

# PNG, not PDF: docs/tex/*.pdf is gitignored as a build product, where the
# note's figures are tracked as PNG.

KP = 24
OUT = "docs/tex"


def rate_on(V, kq):
    """ln(ds/dp) and ds/dp of the K_q fit, at phases p = delta_p sinh(V)."""
    tr = salt_track(0.0, n_per_side=4000)
    fit = fit_track(tr, KP, kq)
    P = D.DELTA_P * np.sinh(V)
    d = D.grid_data(k_psi=KP, k_q=kq)
    d.update(a=np.zeros(kq - 1), b=np.asarray(fit["q"]), c=np.zeros(KP - 1),
             d=np.asarray(fit["psi"]))
    d.update(D.params_data(x_s=float(fit["x_s"]), x_p=0.0, x_psi=0.0))
    d.update(z=0.0)
    d.update(D.eval_points(p_eval=P, s_eval=np.array([0.0, 1.0]),
                           dt=np.array([0.0, 0.0]), band=np.array([1, 2])))
    s = np.asarray(fixed_param("forward", d)["s_of_p"], dtype=float)
    r = np.diff(s) / np.diff(P)
    return r, np.log(r), float(fit["x_s"]), s


def turning_on(U, kpsi):
    """Psi_0 at arclengths s = delta sinh(U), from the K_psi fit.

    forward.stan reports psi at the requested arclengths, and
    psi = -pi/2 + Psi_0 by Eq. (cturn), so Psi_0 = psi + pi/2."""
    tr = salt_track(0.0, n_per_side=4000)
    fit = fit_track(tr, kpsi, 16)
    S = D.DELTA * np.sinh(U)
    d = D.grid_data(k_psi=kpsi, k_q=16)
    d.update(a=np.zeros(15), b=np.asarray(fit["q"]), c=np.zeros(kpsi - 1),
             d=np.asarray(fit["psi"]))
    d.update(D.params_data(x_s=float(fit["x_s"]), x_p=0.0, x_psi=0.0))
    d.update(z=0.0)
    d.update(D.eval_points(p_eval=np.array([0.0, 1.0]), s_eval=S,
                           dt=np.array([0.0, 0.0]), band=np.array([1, 2])))
    return np.asarray(fixed_param("forward", d)["psi"], dtype=float) + np.pi / 2


def track_rate(V):
    tr = salt_track(0.0, n_per_side=4000)
    P = D.DELTA_P * np.sinh(V)
    inw = (tr["p"] >= D.MODEL_WINDOW[0]) & (tr["p"] <= D.MODEL_WINDOW[1])
    s = np.interp(P, tr["p"][inw], tr["s"][inw])
    r = np.diff(s) / np.diff(P)
    return r, np.log(r), s


def main():
    vlo, vhi = np.arcsinh(np.array(D.MODEL_WINDOW) / D.DELTA_P)
    wlo, whi = np.arcsinh(np.array(D.DATA_WINDOW) / D.DELTA_P)
    V = np.linspace(vlo + 1e-6, vhi - 1e-6, 6001)
    vm = 0.5 * (V[1:] + V[:-1])
    pm = D.DELTA_P * np.sinh(vm)
    r_t, y_t, s_t = track_rate(V)
    fits = {kq: rate_on(V, kq) for kq in (8, 16)}
    k = (vm >= wlo) & (vm <= whi)

    # ---- Figure 1: the two shared functions, each in its own coordinate ----
    ulo, uhi = np.arcsinh(np.array(D.ARC_RANGE) / D.DELTA)
    tr = salt_track(0.0, n_per_side=4000)
    u_tr = np.arcsinh(tr["s"] / D.DELTA)
    inu = (u_tr >= ulo) & (u_tr <= uhi)
    U = np.linspace(ulo + 1e-6, uhi - 1e-6, 4001)
    psi_t = np.interp(U, u_tr[inu], tr["psi"][inu] + np.pi / 2)
    turns = {kp: turning_on(U, kp) for kp in (12, 24)}
    # where the data window lands in u
    s_lo = np.interp(D.DATA_WINDOW[0], tr["p"], tr["s"])
    s_hi = np.interp(D.DATA_WINDOW[1], tr["p"], tr["s"])
    u_lo, u_hi = np.arcsinh(s_lo / D.DELTA), np.arcsinh(s_hi / D.DELTA)
    print(f"   the data window maps to s in [{s_lo:+.2f}, {s_hi:+.2f}] mag, "
          f"u in [{u_lo:+.2f}, {u_hi:+.2f}]")
    ku = (U >= u_lo) & (U <= u_hi)
    for kp in turns:
        print(f"   K_psi={kp}: max |residual| in Psi_0 over the data window's arclengths "
              f"{np.abs(turns[kp][ku] - psi_t[ku]).max():.4f} rad")

    fig, ax = plt.subplots(1, 2, figsize=(10.6, 4.3))
    a = ax[0]
    a.plot(U, psi_t, lw=2.2, color="0.25", label="SALT2 $x_{1}=0$")
    for kp, col in ((12, "#ff7f0e"), (24, "#d62728")):
        a.plot(U, turns[kp], lw=1.3, ls="--", color=col, label=f"$K_\\psi={kp}$")
    for e in (u_lo, u_hi):
        a.axvline(e, color="0.4", lw=1.0, ls="-.")
    a.axvline(0, color="0.8", lw=0.6)
    a.axhline(0, color="0.8", lw=0.6)
    a.set_xlabel(r"$u=\mathrm{asinh}(s/\delta)$")
    a.set_ylabel(r"$\Psi_0(u) = \psi + \pi/2$   (rad)")
    # the K_psi = 12 fit runs away past the data window's right edge; let it
    # leave the frame rather than flatten everything the panel is about
    a.set_ylim(psi_t.min() - 0.7, psi_t.max() + 1.0)
    a.set_title("the turning, in $u$", fontsize=10)
    a.legend(fontsize=8, frameon=False, loc="upper left")
    a.grid(alpha=0.15)

    a = ax[1]
    a.plot(vm, y_t, lw=2.2, color="0.25", label="SALT2 $x_{1}=0$")
    for kq, col in ((8, "#ff7f0e"), (16, "#d62728")):
        a.plot(vm, fits[kq][1], lw=1.3, ls="--", color=col, label=f"$K_q={kq}$")
    a.axhline(fits[16][2], color="#2ca02c", lw=1.0, ls=":",
              label=f"$x_s={fits[16][2]:.2f}$")
    for e in (wlo, whi):
        a.axvline(e, color="0.4", lw=1.0, ls="-.")
    a.axvline(0, color="0.8", lw=0.6)
    a.set_xlabel(r"$v=\mathrm{asinh}(p/\delta_p)$")
    a.set_ylabel(r"$\ln(\mathrm{d}s/\mathrm{d}p)=x_s+q_0(v)$")
    a.set_title("the traversal, in $v$", fontsize=10)
    a.legend(fontsize=8, frameon=False, loc="upper center", ncol=2)
    a.grid(alpha=0.15)
    fig.tight_layout()
    fig.savefig(f"{OUT}/traversal-coordinates.png", dpi=200)
    print(f"wrote {OUT}/traversal-coordinates.png")
    for kq in fits:
        print(f"   K_q={kq}: max |residual| in ln rate over the data window "
              f"{np.abs(fits[kq][1][k] - y_t[k]).max():.3f}, "
              f"rms {np.sqrt(np.mean((fits[kq][1][k] - y_t[k])**2)):.4f}")

    # ---- Figure 2: the integrand and the rule ------------------------------
    m, J = D.M_PER_SIDE, D.J_POINTS
    edges = np.concatenate([vlo * (m - np.arange(m + 1)) / m,
                            vhi * np.arange(1, m + 1) / m])
    xg, _ = D.gl_table(J)
    nodes = np.concatenate([0.5 * (edges[q + 1] + edges[q])
                            + 0.5 * (edges[q + 1] - edges[q]) * xg for q in range(2 * m)])
    jac = D.DELTA_P * np.cosh(vm)
    fig, a = plt.subplots(figsize=(8.4, 4.0))
    a.semilogy(vm, r_t * jac, lw=2.2, color="0.25",
               label=r"$\mathrm{d}s/\mathrm{d}v=(\mathrm{d}s/\mathrm{d}p)\,\delta_p\cosh v$")
    a.semilogy(vm, fits[16][0] * jac, lw=1.3, ls="--", color="#d62728", label="$K_q=16$")
    for e in edges:
        a.axvline(e, color="#1f77b4", lw=0.7, alpha=0.5)
    a.plot(nodes, np.full_like(nodes, a.get_ylim()[0] * 1.4), "|", ms=8, color="#1f77b4")
    for e in (wlo, whi):
        a.axvline(e, color="0.4", lw=1.0, ls="-.")
    a.plot([], [], color="#1f77b4", lw=0.7, label=f"$2m={2*m}$ panels, $J={J}$ nodes each")
    a.plot([], [], color="0.4", lw=1.0, ls="-.", label="data window")
    a.set_xlabel(r"$v=\mathrm{asinh}(p/\delta_p)$")
    a.set_ylabel(r"$\mathrm{d}s/\mathrm{d}v$  (mag)")
    a.set_title("the integrand the quadrature sees, and where its nodes fall", fontsize=10)
    a.legend(fontsize=8, frameon=False, loc="upper center")
    a.grid(alpha=0.15, which="both")
    fig.tight_layout()
    fig.savefig(f"{OUT}/traversal-quadrature.png", dpi=200)
    print(f"wrote {OUT}/traversal-quadrature.png")
    inner = (nodes > edges[m - 1]) & (nodes < edges[m + 1])
    print(f"   {int(inner.sum())} of {2*m*J} nodes in the two panels either side of zero")


if __name__ == "__main__":
    main()
