"""The declared constants of the note, and the data dictionaries Stan takes.

Every number here is one the note states in Sections 2.3, 4 and 5.1; nothing
is derived.  The bases, the grid and the forward map live in ``stan/`` and are
never evaluated in Python.
"""

from __future__ import annotations

import numpy as np

# Section 2.3
# The arclength scale below which turning cannot be seen, and so need not be
# resolved: a feature confined to |s| < L moves a unit-speed curve by at most
# L for the epochs inside it, and for those outside contributes a constant
# translation per branch, which m_peak and colr absorb exactly.  Measured by
# pushing a localized bump through the forward map and letting the placement
# absorb what it can, the detection threshold is 5e-3 mag; at K_psi = 12 the
# basis cannot form a feature narrower than ~1e-2 anyway.  At the old 3e-4 a
# third of the turning basis covered |s| < 0.006 mag, where 12 of 1568 epochs
# lie.
DELTA = 5e-3                 # asinh scale in arclength, mag, Eq. (umap)
DELTA_P = 1.0                # asinh scale in phase, days, Eq. (vmap)
MODEL_WINDOW = (-20.0, 45.0)  # the modelling window, days
DATA_WINDOW = (-15.0, 40.0)   # the selection window of Section 3, days
# Where <psi_1, psi_1> is taken, mag, Section 2.3.  The span must be CONTAINED
# WITHIN the arclengths the sample reaches, never contain them: amplitude Psi_1
# puts inside the span but beyond the data counts toward the norm and costs the
# likelihood nothing, so the normalization is spent on invisible structure and
# x_psi grows to compensate.  Measured both ways -- a span inside the reach put
# 95% of the norm where supernovae are and gave x_psi a posterior sd of 0.67; a
# span 1.2 mag wider than the reach put 28% there, Psi_1 spiking at the unseen
# edge, and x_psi's sd rose to 1.42, wider than its own prior.
GAUGE_SPAN = (-6.0, 5.0)
ARC_RANGE = (-8.0, 8.0)      # [s_min, s_max], mag, Section 3: the 442-object fit reaches -6.99 at -15 d and +7.35 at +40 d

# Section 4
EPS = 0.01                   # phase clamp width, fraction of the window
EPS_T = 0.02                 # taper turn-off, fraction of the window
EPS_S = 1e-3                 # arclength clamp width, fraction of the interval
ZP = 30.0                    # the archive zero point

# Section 5.1
# Panels on each side of zero, M = 10 in all.  The rule was accurate to
# 4e-7 mag at M = 20 against a scatter of 0.038; halving it gives 2e-4, still
# 190 times below the noise, and the per-object node sweeps are about
# two thirds of a gradient.  Raising delta shortens the u range, which
# narrows every panel further at fixed M.
M_PER_SIDE = 5
J_POINTS = 8                 # Gauss-Legendre points per panel

# Section 5
# A \salt\-derived mean log rate in v, at the sample-median x_1.  No longer a
# prior centre for the hierarchy, which names the mean as s_bar and constrains
# the per-supernova deviations to sum to zero; kept for single_sn.stan, where
# one object has no population to sum over, and as a starting value.
XS_BAR = -3.0
K_PSI = 24                   # turning order, Section 2.3
K_Q = 16                     # traversal order at stage 2, remeasured in v by scripts/salt_orders.py


def gl_table(n: int = J_POINTS) -> tuple[np.ndarray, np.ndarray]:
    """The ``n``-point Gauss-Legendre rule on ``[-1, 1]``."""
    return np.polynomial.legendre.leggauss(n)


def grid_data(k_psi: int = K_PSI, k_q: int = K_Q, m_per_side: int = M_PER_SIDE,
              j_points: int = J_POINTS) -> dict:
    """The block ``stan/dgsn_grid_data.stan`` reads."""
    xg, wg = gl_table(j_points)
    return dict(
        Kpsi=int(k_psi), Kq=int(k_q), m=int(m_per_side), J=int(j_points),
        xg=xg, wg=wg,
        delta=DELTA, delta_p=DELTA_P,
        pmin=MODEL_WINDOW[0], pmax=MODEL_WINDOW[1],
        smin=ARC_RANGE[0], smax=ARC_RANGE[1],
        gsmin=GAUGE_SPAN[0], gsmax=GAUGE_SPAN[1],
        eps=EPS, eps_t=EPS_T, eps_s=EPS_S, zp=ZP,
    )


def zero_globals(k_psi: int = K_PSI, k_q: int = K_Q) -> dict:
    """``a = b = c = d = 0``: the reference supernova is a straight line
    walked at unit rate.  The starting point of the closed-form tests."""
    return dict(a=np.zeros(k_q - 1), b=np.zeros(k_q - 1),
                c=np.zeros(k_psi - 1), d=np.zeros(k_psi - 1))


def params_data(m_peak: float = 0.0, colr: float = 0.0, x_s: float = 0.0,
                x_p: float = 0.0, x_psi: float = 0.0, delta_t: float = 0.0) -> dict:
    """The six of Table 1 as ``forward.stan`` data, and as ``single_sn`` inits."""
    return dict(m_peak=float(m_peak), colr=float(colr), x_s=float(x_s),
                x_p=float(x_p), x_psi=float(x_psi), delta_t=float(delta_t))


def eval_points(p_eval=(), s_eval=(), dt=(), band=()) -> dict:
    """Where ``forward.stan`` is asked for ``s(p)``, ``gamma(s)`` and fluxes."""
    p_eval = np.atleast_1d(np.asarray(p_eval, dtype=float))
    s_eval = np.atleast_1d(np.asarray(s_eval, dtype=float))
    dt = np.atleast_1d(np.asarray(dt, dtype=float))
    band = np.atleast_1d(np.asarray(band, dtype=int))
    return dict(Np=len(p_eval), p_eval=p_eval, Ns=len(s_eval), s_eval=s_eval,
                N=len(dt), dt=dt, band=band)


def epochs_data(dt, f, sig_f, band, z: float) -> dict:
    """The epoch set of Section 4 as ``single_sn.stan`` data."""
    dt = np.asarray(dt, dtype=float)
    return dict(N=len(dt), dt=dt, f=np.asarray(f, dtype=float),
                sig_f=np.asarray(sig_f, dtype=float), band=np.asarray(band, dtype=int),
                z=float(z))


def prior_data(prior_on: bool = True, mu_p: float = 0.0, mu_psi: float = 0.0,
               mu_delta: float = 0.0, s_delta: float = 1.0, xs_bar: float = XS_BAR,
               m_bright: float = 0.0) -> dict:
    """The priors of Section 5 with the population hyperparameters as data.

    One supernova identifies none of the hyperparameters; stage 3 samples them.
    ``prior_on = False`` is for the stage-1 fits that must hold to numerical
    precision.
    """
    return dict(prior_on=int(bool(prior_on)), mu_p=float(mu_p), mu_psi=float(mu_psi),
                mu_delta=float(mu_delta), s_delta=float(s_delta), xs_bar=float(xs_bar),
                m_bright=float(m_bright))
