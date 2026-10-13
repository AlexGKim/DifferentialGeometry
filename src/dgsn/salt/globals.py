"""The global functions from SALT2 tracks: the start of the fit (Section 5)
and the data-supplied globals of stage 2.

The basis is evaluated by ``forward.stan`` at the track's points and the
coefficients follow by linear least squares on the tangent angle in ``u`` and
the log rate in ``v``, then are refined by fitting the track's light curve in
magnitude space (``track_fit.stan``), which is the objective that matters:
the log rate at the anchor is deep but contributes little arclength, so a fit
that resolves it there at the expense of the late phases reproduces the track
poorly.
``b, d`` come from the track at the sample-median ``x1``; ``a, c`` from the
difference between the tracks at ``x1 = +1`` and ``-1``, which the gauge then
normalizes.  Orthogonalizing the mean against the departure mode moves the
mean's component along it into the per-supernova amplitude, so the reference
track is reproduced at the ``x_p, x_psi`` returned here and not at zero.
"""

from __future__ import annotations

import numpy as np

from dgsn.stan import data as D
from dgsn.stan.run import fixed_param


def design_at(track: dict, k_psi: int, k_q: int, m_per_side: int = D.M_PER_SIDE) -> dict:
    """``Xu`` at the track's arclengths and ``Xv`` at its phases, plus the
    fixed grid's basis and weights for the sign conditions."""
    data = D.grid_data(k_psi=k_psi, k_q=k_q, m_per_side=m_per_side)
    data.update(D.zero_globals(k_psi, k_q))
    data.update(D.params_data())
    data.update(z=0.0)
    inside = (track["s"] > D.ARC_RANGE[0]) & (track["s"] < D.ARC_RANGE[1])
    data.update(D.eval_points(p_eval=track["p"],
                              s_eval=np.concatenate([track["s"][inside], [D.GAUGE_SPAN[1]]])))
    out = fixed_param("forward", data)
    return dict(inside=inside, Xu=out["Xu"][:-1], Xu_gsmax=out["Xu"][-1], Xv=out["Xv"],
                Xv_nodes=out["Xv_nodes_out"], v_nodes=out["v_nodes_out"], w_v=out["w_v_out"],
                Gs=out["Gs_out"])


def _weights(x: np.ndarray) -> np.ndarray:
    """Square roots of the trapezoid widths, so least squares in these weights
    approximates the L2 fit under the flat measure of ``x``."""
    dx = np.gradient(x)
    return np.sqrt(np.abs(dx))


# The ridge of the prior's form applied to the template fits, lambda k^{2 alpha}
# c_k^2.  A template reaches s ~ +3.3 at +40 d while the interval runs to +8,
# so a quarter of the turning interval holds no data and the unregularized
# fit leaves arbitrary coefficients there, which the difference of two such
# fits then inherits as the departure mode.  The scales are set against the
# data terms: the weighted normal matrix of the least squares is close to
# the identity on constrained modes, and the magnitude-space fit at 0.01 mag
# weights a constrained mode at ~1e6 per unit coefficient; both ridges shrink
# the highest mode by a few percent and an unconstrained one to zero.
RIDGE_ALPHA = 2.0
RIDGE_LSTSQ = 1e-7
RIDGE_MAG = 1e-4


def _ridge_rows(n: int, scale: float) -> np.ndarray:
    k = np.arange(1, n + 1)
    return np.diag(np.sqrt(scale * k ** (2 * RIDGE_ALPHA)))


def fit_track(track: dict, k_psi: int, k_q: int, m_per_side: int = D.M_PER_SIDE) -> dict:
    """Raw coefficients of one track: ``psi + pi/2`` on the turning basis in
    ``u``, ``log ds/dp`` less its mean in ``v`` on the traversal basis, with
    the ridge above."""
    des = design_at(track, k_psi, k_q, m_per_side)
    inside = des["inside"]
    wu = _weights(track["u"][inside])
    wv = _weights(track["v"])
    y_psi = track["psi"][inside] + np.pi / 2
    vlo, vhi = np.arcsinh(np.array(D.MODEL_WINDOW) / D.DELTA_P)
    x_s = float(np.trapezoid(track["log_rate"], track["v"]) / (vhi - vlo))
    y_rate = track["log_rate"] - x_s
    Au = np.vstack([wu[:, None] * des["Xu"], _ridge_rows(k_psi - 1, RIDGE_LSTSQ)])
    Av = np.vstack([wv[:, None] * des["Xv"], _ridge_rows(k_q - 1, RIDGE_LSTSQ)])
    coef_psi = np.linalg.lstsq(Au, np.concatenate([wu * y_psi, np.zeros(k_psi - 1)]), rcond=None)[0]
    coef_q = np.linalg.lstsq(Av, np.concatenate([wv * y_rate, np.zeros(k_q - 1)]), rcond=None)[0]
    return dict(psi=coef_psi, q=coef_q, x_s=x_s,
                rms_psi=float(np.sqrt(np.mean((des["Xu"] @ coef_psi - y_psi) ** 2))),
                rms_rate=float(np.sqrt(np.mean((des["Xv"] @ coef_q - y_rate) ** 2))),
                design=des)


def fix_gauge(a_raw, b_raw, c_raw, d_raw, design: dict) -> dict:
    """Eq. (gauge) and the two signs of Section 2.3, on the coefficients.

    The traversal pair under the flat measure in ``v``, where the basis is
    orthonormal and the conditions are the Euclidean ones; the turning pair
    under ``ds``, the measure the light curve sees, through the Gram matrix
    ``G_s`` of the basis in that measure."""
    a = a_raw / np.linalg.norm(a_raw)
    x_p = float(b_raw @ a)
    b = b_raw - x_p * a
    Gs = design["Gs"]
    c = c_raw / np.sqrt(c_raw @ Gs @ c_raw)
    x_psi = float(d_raw @ Gs @ c)
    d = d_raw - x_psi * c
    # int_0^{v_max} q1 dv > 0, on the fixed nodes right of zero
    right = design["v_nodes"] > 0
    if float(design["w_v"][right] @ (design["Xv_nodes"][right] @ a)) < 0:
        a, x_p = -a, -x_p
    # Psi_1 > 0 at the gauge span's upper edge
    if float(design["Xu_gsmax"] @ c) < 0:
        c, x_psi = -c, -x_psi
    return dict(a=a, b=b, c=c, d=d, x_p=x_p, x_psi=x_psi)


def fit_globals(k_psi: int = D.K_PSI, k_q: int = D.K_Q, x1_ref: float = 0.0,
                m_per_side: int = D.M_PER_SIDE) -> dict:
    """The four gauge-fixed global functions from SALT2, and the point
    ``(x_s, x_p, x_psi)`` at which they reproduce the reference track."""
    from dgsn.salt.tracks import salt_track

    fits = {}
    for x1 in (x1_ref, x1_ref + 1.0, x1_ref - 1.0):
        tr = salt_track(x1)
        fits[x1] = refine_track(tr, fit_track(tr, k_psi, k_q, m_per_side), k_psi, k_q, m_per_side)
    ref, plus, minus = fits[x1_ref], fits[x1_ref + 1.0], fits[x1_ref - 1.0]
    out = fix_gauge(plus["q"] - minus["q"], ref["q"], plus["psi"] - minus["psi"], ref["psi"],
                    ref["design"])
    out.update(x_s=ref["x_s"], x1_ref=x1_ref, k_psi=k_psi, k_q=k_q,
               lightcurve=ref["lightcurve"], lightcurve_rms=ref["lightcurve_rms"],
               plus_lightcurve=plus["lightcurve"], minus_lightcurve=minus["lightcurve"])
    return out


def track_error(track: dict, coef_psi, coef_q, x_s: float, k_psi: int, k_q: int,
                m_per_side: int = D.M_PER_SIDE, window=D.DATA_WINDOW, placement=None) -> dict:
    """How well truncated coefficients reproduce a track, in magnitudes, over
    ``window``: the traversal error ``s_model(p) - s_track(p)``, the curve
    error ``gamma(s_track) - (g, r)`` at the track's own arclengths, and the
    light-curve error with both compounded, the quantity the note's orders
    are set by.  The coefficients are raw (ungauged): a track is one
    supernova and the gauge is a relabeling."""
    p, s = track["p"], track["s"]
    inw = (p >= window[0]) & (p <= window[1])
    ins = inw & (s > D.ARC_RANGE[0]) & (s < D.ARC_RANGE[1])
    i0 = int(np.argmin(np.abs(p)))
    g0, r0 = track["g"][i0], track["r"][i0]
    m_peak, colr, delta_t = placement if placement is not None else (g0, g0 - r0, 0.0)
    data = D.grid_data(k_psi=k_psi, k_q=k_q, m_per_side=m_per_side)
    data.update(a=np.zeros(k_q - 1), b=np.asarray(coef_q), c=np.zeros(k_psi - 1),
                d=np.asarray(coef_psi))
    data.update(D.params_data(m_peak=m_peak, colr=colr, x_s=x_s, delta_t=delta_t))
    data.update(z=0.0)
    dt = np.concatenate([p[inw], p[inw]])
    band = np.concatenate([np.ones(inw.sum(), int), 2 * np.ones(inw.sum(), int)])
    data.update(D.eval_points(p_eval=p[inw], s_eval=s[ins], dt=dt, band=band))
    out = fixed_param("forward", data)
    mag = D.ZP - out["lf"][:, 0] / (0.4 * np.log(10))
    err_lc = np.abs(mag - np.concatenate([track["g"][inw], track["r"][inw]]))
    err_curve = np.abs(out["gamma"] - np.stack([track["g"][ins] - g0, track["r"][ins] - r0], 1))
    return dict(traversal=float(np.abs(out["s_of_p"] - s[inw]).max()),
                curve=float(err_curve.max()),
                lightcurve=float(err_lc.max()),
                lightcurve_rms=float(np.sqrt(np.mean(err_lc ** 2))),
                s_of_p=out["s_of_p"], gamma=out["gamma"], mag=mag)


def track_epochs(track: dict, n_per_band: int = 200, window=D.DATA_WINDOW,
                 sig_mag: float = 0.01) -> dict:
    """A template track as noise-free epochs over ``window``, both bands,
    sampled uniformly in ``v``, with errors proportional to flux so that the
    fit weights every epoch equally in magnitude."""
    p, g, r = track["p"], track["g"], track["r"]
    inw = (p >= window[0]) & (p <= window[1])
    idx = np.unique(np.linspace(0, inw.sum() - 1, n_per_band).round().astype(int))
    pw = p[inw][idx]
    f = 10 ** (-0.4 * (np.concatenate([g[inw][idx], r[inw][idx]]) - D.ZP))
    return dict(dt=np.concatenate([pw, pw]), f=f, sig_f=0.4 * np.log(10) * sig_mag * f,
                band=np.concatenate([np.ones(len(pw), int), 2 * np.ones(len(pw), int)]), z=0.0)


def refine_track(track: dict, fit: dict, k_psi: int, k_q: int,
                 m_per_side: int = D.M_PER_SIDE, window=D.DATA_WINDOW) -> dict:
    """Refine a least-squares start by fitting the track's light curve in
    magnitude space with ``track_fit.stan``: ``b, d`` and the placement, the
    departure modes off.  Returns the raw coefficients, the placement and
    the maximum and rms light-curve error over ``window`` in magnitudes."""
    from dgsn.stan.run import optimize

    ep = track_epochs(track, window=window)
    data = D.grid_data(k_psi=k_psi, k_q=k_q, m_per_side=m_per_side)
    data.update(D.epochs_data(ep["dt"], ep["f"], ep["sig_f"], ep["band"], 0.0, 0.0, 0.0))
    data.update(ridge=RIDGE_MAG, ridge_alpha=RIDGE_ALPHA)
    i0 = int(np.argmin(np.abs(track["p"])))
    start = dict(m_peak=float(track["g"][i0]), colr=float(track["g"][i0] - track["r"][i0]),
                 x_s=fit["x_s"], delta_t=0.0, b=fit["q"], d=fit["psi"])
    # CmdStan's default tolerances: the tables want 1e-4 mag, not the 1e-12
    # the exactness tests ask of optimize().
    # A stalled line search on a hard track (fast decliners at low order)
    # still leaves the best point found, which is the number the tables want.
    # Where the least-squares seed rings --- exp of a polynomial fit to a log
    # rate that diverges at the template's edge --- the start is not finite;
    # a constant rate at the track's mean log rate is then the seed.
    opts = dict(polish=False, tol_rel_grad=1e7, tol_grad=1e-8, tol_param=1e-8,
                tol_obj=1e-12, tol_rel_obj=1e4, iter=2000)
    try:
        out = optimize("track_fit", data, start, **opts)
    except RuntimeError as err:
        if "Initialization failed" not in str(err):
            raise
        out = optimize("track_fit", data, dict(start, b=np.zeros(k_q - 1)), **opts)
    b = np.array([out[f"b[{k}]"] for k in range(1, k_q)])
    d = np.array([out[f"d[{k}]"] for k in range(1, k_psi)])
    err = track_error(track, d, b, out["x_s"], k_psi, k_q, m_per_side, window,
                      placement=(out["m_peak"], out["colr"], out["delta_t"]))
    return dict(q=b, psi=d, x_s=out["x_s"], m_peak=out["m_peak"], colr=out["colr"],
                delta_t=out["delta_t"], lightcurve=err["lightcurve"],
                lightcurve_rms=err["lightcurve_rms"], design=fit["design"])
