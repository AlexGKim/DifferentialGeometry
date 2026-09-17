"""Template light curves, from ``sncosmo``: the SALT2 tracks that place the
start of the fit and supply the stage-2 global functions, the per-supernova
starting point, and the Hsiao light curves the floors of Section 3 are
measured on.

Everything here is synthesis: phases, magnitudes and their derivatives.  The
bases and the forward map are Stan's (``dgsn.stan``); nothing here fits.
Initialization is exempt from the rule that SALT2 stays downstream
(Section 5), and the floors are pushed through the model before any fit.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid

from dgsn.stan import data as D

BANDS = ("ztfg", "ztfr")
R_V = 3.1


def _model(source: str, ebv: float = 0.0):
    import sncosmo

    if ebv:
        return sncosmo.Model(source=source, effects=[sncosmo.F99Dust(r_v=R_V)],
                             effect_names=["mw"], effect_frames=["obs"])
    return sncosmo.Model(source=source)


def phase_grid(n_per_side: int = 10000) -> np.ndarray:
    """Phases uniform in ``v`` on each side of zero, zero included, so that
    the sampling follows the resolution the traversal basis has."""
    vlo, vhi = np.arcsinh(np.array(D.MODEL_WINDOW) / D.DELTA_P)
    v = np.concatenate([np.linspace(vlo, 0.0, n_per_side + 1)[:-1],
                        np.linspace(0.0, vhi, n_per_side + 1)])
    return D.DELTA_P * np.sinh(v)


def g_maximum(model, t_guess: float = 0.0, half_width: float = 5.0) -> float:
    """The time of the ``g`` maximum, refined by a parabola through the
    minimum magnitude on a fine grid."""
    t = np.linspace(t_guess - half_width, t_guess + half_width, 20001)
    g = model.bandmag("ztfg", "ab", t)
    i = int(np.argmin(g))
    i = min(max(i, 1), len(t) - 2)
    y0, y1, y2 = g[i - 1], g[i], g[i + 1]
    return float(t[i] + 0.5 * (t[1] - t[0]) * (y0 - y2) / (y0 - 2 * y1 + y2))


def track(model, t_max: float, n_per_side: int = 10000) -> dict:
    """A template's track through magnitude space, as the model sees it.

    Phase is measured from the ``g`` maximum.  Returns the phases, the two
    magnitudes, the arclength from the anchor, the tangent angle with
    ``psi(0) = -pi/2``, the log rate, and the two expansion coordinates.
    """
    p = phase_grid(n_per_side)
    with np.errstate(divide="ignore", invalid="ignore"):
        g = model.bandmag("ztfg", "ab", t_max + p)
        r = model.bandmag("ztfr", "ab", t_max + p)
    # A template returns zero flux, hence an infinite magnitude, where it has
    # no coverage: SALT2 begins 20 d before B maximum, which is before the
    # modelling window's start when measured from g maximum.  The track is
    # what the template does cover; the basis fit does not reach past it.
    ok = np.isfinite(g) & np.isfinite(r)
    p, g, r = p[ok], g[ok], r[ok]
    dg, dr = np.gradient(g, p), np.gradient(r, p)
    speed = np.hypot(dg, dr)
    s = cumulative_trapezoid(speed, p, initial=0.0)
    i0 = int(np.argmin(np.abs(p)))
    s = s - s[i0]
    psi = np.unwrap(np.arctan2(dr, dg))
    psi = psi - psi[i0] - np.pi / 2               # the convention of Section 2.1
    r_slope_sign = float(np.sign(dr[i0]))         # must be negative: r still brightening
    return dict(p=p, g=g, r=r, s=s, psi=psi, log_rate=np.log(speed),
                u=np.arcsinh(s / D.DELTA), v=np.arcsinh(p / D.DELTA_P),
                r_slope_sign=r_slope_sign, psi_at_anchor=float(psi[i0]),
                p_range=(float(p[0]), float(p[-1])))


def salt_track(x1: float, n_per_side: int = 10000) -> dict:
    """The SALT2 track at ``x1``, ``c = 0``, ``z = 0``."""
    m = _model("salt2")
    m.set(z=0.0, t0=0.0, x0=1.0, x1=x1, c=0.0)
    return track(m, g_maximum(m), n_per_side)


def per_sn_start(row) -> dict:
    """Where a supernova's chain starts (Section 5): ``delta_t = 0``,
    ``m_peak`` and ``colr`` at what the archive's SALT2 fit implies on the
    observed scale, Galactic dust included since the data carry it,
    ``x_s`` at the sample value, and the departure amplitudes at zero."""
    m = _model("salt2", ebv=float(row.mwebv) if np.isfinite(row.mwebv) else 0.0)
    m.set(z=float(row.redshift), t0=float(row.t0), x0=float(row.x0),
          x1=float(row.x1), c=float(row.c))
    if m.param_names.count("mwebv"):
        m.set(mwebv=float(row.mwebv))
    t_g = g_maximum(m, float(row.t0), 5.0 * (1 + float(row.redshift)))
    g = float(m.bandmag("ztfg", "ab", t_g))
    r = float(m.bandmag("ztfr", "ab", t_g))
    return D.params_data(m_peak=g, colr=g - r, x_s=D.XS_BAR, x_p=0.0, x_psi=0.0, delta_t=0.0)


def hsiao_epochs(z: float = 0.0, ebv: float = 0.0, cadence: float = 2.0,
                 snr_peak: float = 50.0, m_peak: float = 18.0) -> dict:
    """Noise-free epochs from the Hsiao template on a ZTF-like cadence, in the
    observer frame, for the floors of Section 3.  Errors are set at a fixed
    signal-to-noise at peak, so that the fit weights epochs as ZTF would.
    Returns the epoch set as ``single_sn.stan`` takes it, with ``t0`` at the
    ``g`` maximum of the ``z = 0``, dust-free template."""
    m = _model("hsiao", ebv=ebv)
    m.set(z=z, t0=0.0)
    if ebv:
        m.set(mwebv=ebv)
    m.set_source_peakmag(m_peak, "ztfg", "ab")
    ref = _model("hsiao")
    ref.set(z=0.0, t0=0.0)
    t_g = g_maximum(ref)
    p_obs = np.arange(D.DATA_WINDOW[0], D.DATA_WINDOW[1] + 1e-9, cadence) * (1 + z)
    dt = np.concatenate([p_obs, p_obs + 0.5 * cadence])
    band = np.concatenate([np.ones(len(p_obs), int), 2 * np.ones(len(p_obs), int)])
    f = np.empty(len(dt))
    for b in (1, 2):
        sel = band == b
        f[sel] = m.bandflux(BANDS[b - 1], t_g + dt[sel], zp=D.ZP, zpsys="ab")
    sig_f = np.full(len(dt), f.max() / snr_peak)
    return dict(dt=dt, f=f, sig_f=sig_f, band=band, z=z)
