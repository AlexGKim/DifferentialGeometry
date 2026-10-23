"""The two figures of Section 9 of docs/tex/implementation.tex.

    python scripts/plot_traversal.py

Why they exist.  The traversal is expanded in Legendre polynomials of
v = asinh(p/delta_p) rather than of phase, and integrated by a rule whose
nodes are laid in v.  Both choices are hard to motivate in words and
immediate in a picture: the rate is a cusp in phase and a smooth broad curve
in v, and the rule puts a fifth of its nodes in the two panels either side of
zero, which is where the cusp is.

The \salt\ track at x1 = 0 is used as an illustration of what the basis has
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

    # ---- Figure 1: the same function in the two coordinates ----------------
    fig, ax = plt.subplots(1, 2, figsize=(10.6, 4.3))
    a = ax[0]
    a.semilogy(pm, r_t, lw=2.2, color="0.25")
    a.semilogy(pm, fits[16][0], lw=1.3, ls="--", color="#d62728")
    j = np.argmin(r_t[k])
    a.plot(pm[k][j], r_t[k][j], "o", ms=6, mfc="none", color="#2ca02c")
    a.annotate(f"{r_t[k][j]:.3f} mag/d\nat {pm[k][j]:+.1f} d", (pm[k][j], r_t[k][j]),
               textcoords="offset points", xytext=(12, 4), fontsize=8, color="#2ca02c")
    for e in D.DATA_WINDOW:
        a.axvline(e, color="0.4", lw=1.0, ls="-.")
    a.set_xlabel(r"rest-frame phase $p$  (d)")
    a.set_ylabel(r"$\mathrm{d}s/\mathrm{d}p$  (mag/d)")
    a.set_title("in phase: a cusp a day wide", fontsize=10)
    a.grid(alpha=0.15, which="both")

    a = ax[1]
    a.plot(vm, y_t, lw=2.2, color="0.25", label=r"SALT2 $x_{1}=0$")
    for kq, col in ((8, "#ff7f0e"), (16, "#d62728")):
        a.plot(vm, fits[kq][1], lw=1.3, ls="--", color=col, label=f"$K_q={kq}$")
    a.axhline(fits[16][2], color="#2ca02c", lw=1.0, ls=":",
              label=f"$x_s={fits[16][2]:.2f}$")
    for e in (wlo, whi):
        a.axvline(e, color="0.4", lw=1.0, ls="-.")
    a.set_xlabel(r"$v=\mathrm{asinh}(p/\delta_p)$")
    a.set_ylabel(r"$\ln(\mathrm{d}s/\mathrm{d}p)=x_s+q_0(v)$")
    a.set_title("in $v$: smooth, and what the basis fits", fontsize=10)
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
