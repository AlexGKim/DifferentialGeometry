"""Data and initial values for ``stan/hier.stan`` (stages 3 and 4).

Epochs are ragged by supernova and passed flat with 1-based inclusive
``start``/``stop``.  Initialization follows Section 5: per supernova as
stage 2, the global functions from the ridge-regularized SALT2 fits, with the
non-centred ``z`` recovered from the gauge-fixed coefficients at the initial
``tau`` and ``alpha``.
"""

from __future__ import annotations

import numpy as np

from dgsn.data.epochs import likelihood_epochs
from dgsn.data.sample import analysis_sample
from dgsn.salt.globals import fit_globals
from dgsn.salt.tracks import per_sn_start
from dgsn.stan import data as D

INIT_ALPHA = 2.0        # alpha at the start; the ridge of the SALT2 fits uses the same
INIT_SIGMA = 0.05       # mag; a scatter the first iterations can live with


def build(names, k_psi: int, k_q: int, grainsize: int = 1, m_per_side: int = D.M_PER_SIDE,
          sample=None) -> dict:
    """The data block of ``hier.stan`` for the objects ``names``."""
    sample = analysis_sample() if sample is None else sample
    dt, f, sig_f, band, z, m_bright, start, stop = [], [], [], [], [], [], [], []
    pos = 1
    for name in names:
        row = sample.loc[name]
        ep = likelihood_epochs(name, row.t0, row.redshift)
        n = len(ep["dt"])
        dt.append(ep["dt"]); f.append(ep["f"]); sig_f.append(ep["sig_f"]); band.append(ep["band"])
        z.append(ep["z"]); m_bright.append(ep["m_bright"])
        start.append(pos); stop.append(pos + n - 1); pos += n
    data = D.grid_data(k_psi=k_psi, k_q=k_q, m_per_side=m_per_side)
    data.update(
        N_sn=len(names), N=pos - 1,
        dt=np.concatenate(dt), f=np.concatenate(f), sig_f=np.concatenate(sig_f),
        band=np.concatenate(band).astype(int), start=np.array(start), stop=np.array(stop),
        z=np.array(z), m_bright=np.array(m_bright), xs_bar=D.XS_BAR, grainsize=int(grainsize),
    )
    return data


def inits(names, k_psi: int, k_q: int, x1_ref: float, sample=None, globals_=None,
          tau: dict | None = None) -> dict:
    """Initial values: per-supernova from SALT2 and the sample, global
    functions from the SALT2 fits.  ``tau`` defaults to the norm the SALT2
    coefficients imply at ``INIT_ALPHA``, so that ``z`` starts near unit scale
    --- the centred pair starts at the coefficients themselves, and ``tau``
    at the width their prior wants."""
    sample = analysis_sample() if sample is None else sample
    G = globals_ if globals_ is not None else fit_globals(k_psi, k_q, x1_ref)
    kq, kp = np.arange(1, k_q), np.arange(1, k_psi)
    wq, wp = kq ** -INIT_ALPHA, kp ** -INIT_ALPHA
    if tau is None:
        # the width the centred prior wants at the start: tau^2 = mean((b_k / k^-alpha)^2);
        # a and c do not enter, the gauge having removed their scale
        tau = dict(q=float(np.sqrt(np.mean((G["b"] / wq) ** 2))),
                   psi=float(np.sqrt(np.mean((G["d"] / wp) ** 2))))
    per = [per_sn_start(sample.loc[n]) for n in names]
    return dict(
        z_a=G["a"] / wq, b_til=np.asarray(G["b"], dtype=float),
        z_c=G["c"] / wp, d_til=np.asarray(G["d"], dtype=float),
        tau_q=tau["q"], tau_psi=tau["psi"],
        alpha_q_m=INIT_ALPHA - 0.5, alpha_psi_m=INIT_ALPHA - 0.5,
        sigma_g=INIT_SIGMA, sigma_r=INIT_SIGMA,
        mu_p=G["x_p"], mu_psi=G["x_psi"], mu_delta=0.0, s_delta=1.0,
        m_peak=np.array([p["m_peak"] for p in per]),
        colr=np.array([p["colr"] for p in per]),
        x_s=np.full(len(names), G["x_s"]),
        u_p=np.zeros(len(names)),      # x_p starts at mu_p = the SALT2 reference
        u_psi=np.zeros(len(names)), u_del=np.zeros(len(names)),
    )


PARAM_NAMES = ["m_peak", "colr", "x_s", "x_p", "x_psi", "delta_t"]
