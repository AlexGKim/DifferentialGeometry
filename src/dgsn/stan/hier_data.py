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
          sample=None, R=None, mode: bool = True) -> dict:
    """The data block of ``hier.stan`` for the objects ``names``.

    ``mode=False`` is for ``hier_nomode.stan``, where Eq. (curvlin) carries no
    ``x_psi Psi_1`` term, so the coefficient block is shorter by ``Kpsi - 1``
    and ``R`` with it."""
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
        z=np.array(z), m_bright=np.array(m_bright), grainsize=int(grainsize),
        R=np.eye(n_coef(k_psi, k_q, mode)) if R is None else np.asarray(R, dtype=float),
    )
    return data


def n_coef(k_psi: int, k_q: int, mode: bool = True) -> int:
    """The rotated block: z_a, b_til, z_c, d_til laid end to end, or without
    z_c when the turning mode is absent."""
    return 2 * (k_q - 1) + (2 if mode else 1) * (k_psi - 1)


def pack(*vectors) -> np.ndarray:
    """The coefficient coordinates in the order the program slices them:
    ``z_a, b_til, z_c, d_til``, or ``z_a, b_til, d_til`` without the mode."""
    return np.concatenate([np.asarray(v, dtype=float) for v in vectors])


def unpack(raw, k_psi: int, k_q: int, mode: bool = True) -> dict:
    n1, n3 = k_q - 1, k_psi - 1
    raw = np.asarray(raw, dtype=float)
    if not mode:
        return dict(z_a=raw[:n1], b_til=raw[n1:2*n1], d_til=raw[2*n1:2*n1+n3])
    return dict(z_a=raw[:n1], b_til=raw[n1:2*n1], z_c=raw[2*n1:2*n1+n3], d_til=raw[2*n1+n3:])


def inits(names, k_psi: int, k_q: int, x1_ref: float, sample=None, globals_=None,
          tau: dict | None = None, R=None, mode: bool = True) -> dict:
    """Initial values: per-supernova from SALT2 and the sample, global
    functions from the SALT2 fits.  ``tau`` defaults to the norm the SALT2
    coefficients imply at ``INIT_ALPHA``, so that ``z`` starts near unit scale
    --- the centred pair starts at the coefficients themselves, and ``tau``
    at the width their prior wants."""
    sample = analysis_sample() if sample is None else sample
    G = globals_ if globals_ is not None else fit_globals(k_psi, k_q, x1_ref)
    kq, kp = np.arange(1, k_q), np.arange(1, k_psi)
    # the pivot of hier.stan: tau is the scale at the geometric-mean mode
    wq = (kq / np.exp(np.log(kq).mean())) ** -INIT_ALPHA
    wp = (kp / np.exp(np.log(kp).mean())) ** -INIT_ALPHA
    if tau is None:
        # the width the centred prior wants at the start: tau^2 = mean((b_k / k^-alpha)^2);
        # a and c do not enter, the gauge having removed their scale
        tau = dict(q=float(np.sqrt(np.mean((G["b"] / wq) ** 2))),
                   psi=float(np.sqrt(np.mean((G["d"] / wp) ** 2))))
    per = [per_sn_start(sample.loc[n]) for n in names]
    raw = (pack(G["a"] / wq, G["b"], G["c"] / wp, G["d"]) if mode
           else pack(G["a"] / wq, G["b"], G["d"]))
    out = dict(
        # w = R^-1 raw, so that the start is the same point whatever rotation
        # is in force
        w_coef=raw if R is None else np.linalg.solve(np.asarray(R, dtype=float), raw),
        tau_q=tau["q"], tau_psi=tau["psi"],
        alpha_q_m=INIT_ALPHA - 0.5, alpha_psi_m=INIT_ALPHA - 0.5,
        sigma_g=INIT_SIGMA, sigma_r=INIT_SIGMA,
        mu_p=G["x_p"], mu_delta=0.0, s_delta=1.0,
        m_peak=np.array([p["m_peak"] for p in per]),
        colr=np.array([p["colr"] for p in per]),
        s_bar=G["x_s"],                # the mean log rate; the deviations carry the rest
        x_s_dev=np.zeros(len(names)),  # sums to zero, as sum_to_zero_vector requires
        u_p=np.zeros(len(names)),      # x_p starts at mu_p = the SALT2 reference
        u_del=np.zeros(len(names)),
    )
    if mode:
        out.update(mu_psi=G["x_psi"], u_psi=np.zeros(len(names)))
    return out


def _from_mode(mode: dict, name: str):
    """One parameter out of an ``optimize`` result, scalar or vector.

    CmdStan reports ``tau_q`` but ``u_del[1]``, so the shape has to be
    discovered rather than assumed."""
    if name in mode:
        return float(mode[name])
    vals = []
    i = 1
    while f"{name}[{i}]" in mode:
        vals.append(mode[f"{name}[{i}]"])
        i += 1
    if not vals:
        raise KeyError(f"{name} not in the optimize result")
    return np.array(vals, dtype=float)


def inits_from_optimize(program: str, data: dict, start: dict, k_psi: int, k_q: int,
                        mode: bool = True, verbose: bool = True,
                        polish: bool = True) -> dict:
    """Initial values from Stan's own optimizer on the data, re-placed on the
    ``s_delta`` ridge.

    Why this exists.  Section 5 starts the chains at a \\salt\\ fit, which the
    third standing commitment permits as its one exemption, licensed by
    dispersed-start agreement.  A start found by optimizing the model's own
    log density on the observed epochs needs no exemption: the only template
    left is the one that starts the optimizer, a dependence on a basin rather
    than on a point.

    Why the re-placement.  The constrained-scale mode is *not finite* in
    ``s_delta``.  Scaling ``s_delta`` by any factor while dividing ``u_del``
    by it leaves ``delta_t``, and so the likelihood, exactly unchanged, while
    the prior ``-1/2 sum u_del^2`` rises toward zero; measured along that
    ridge the log density gains 7.3 nats monotonically and is flat to machine
    precision by a factor of thirty, which is where the optimizer stops --- it
    returned ``s_delta`` 5.48 against a posterior 0.186.  Starting a sampler
    there mis-scales the per-supernova timing block by thirty.  Taking
    ``u_del`` to unit rms and ``s_delta`` to match moves along the ridge to
    where the non-centred pair is well scaled, and changes no fitted
    magnitude: ``delta_t`` is preserved exactly.

    ``polish`` is ``optimize``'s restart loop, which exists to resolve the
    direction the light curve is nearly flat along.  On this posterior it
    cannot terminate in ``s_delta``: that direction is flat in likelihood and
    rising in prior, so each restart walks further up the ridge, by fractions
    of a nat, and the re-placement then discards the distance gained.  Four of
    the five passes are spent that way, at about three minutes each on thirty
    objects, so ``polish=False`` is the better default when the mode is only
    wanted as a starting point.

    It is the only pair needing this.  The gauge removes the scale of ``a``
    and ``c``, so ``tau_q`` and ``tau_psi`` are informed by the centred
    ``b_til`` and ``d_til`` and land within a factor of 1.7 of their
    posteriors, as do both scatters.
    """
    from dgsn.stan.run import optimize

    m = optimize(program, data, start, jacobian=False, polish=polish)
    # s_bar and x_s_dev, not x_s: the optimizer reports x_s too, as a
    # transformed parameter, but an inits dict has to name what is sampled.
    keys = ["w_coef", "tau_q", "tau_psi", "alpha_q_m", "alpha_psi_m", "sigma_g", "sigma_r",
            "mu_p", "mu_delta", "s_delta", "s_bar", "x_s_dev", "u_p", "u_del",
            "m_peak", "colr"]
    if mode:
        keys += ["mu_psi", "u_psi"]
    out = {k: _from_mode(m, k) for k in keys}
    # back onto the constraint surface after the roundtrip through text
    out["x_s_dev"] = out["x_s_dev"] - np.mean(out["x_s_dev"])

    u = np.asarray(out["u_del"], dtype=float)
    rms = float(np.sqrt(np.mean(u ** 2)))
    if rms > 0:
        s_old = float(out["s_delta"])
        out["s_delta"] = s_old * rms          # delta_t = mu_delta + s_delta u_del is fixed
        out["u_del"] = u / rms
        if verbose:
            print(f"optimize start: lp {float(m['lp__']):.1f}; "
                  f"s_delta re-placed {s_old:.4g} -> {out['s_delta']:.4g} "
                  f"(rms u_del {rms:.4g} -> 1), delta_t unchanged")
    elif verbose:
        print(f"optimize start: lp {float(m['lp__']):.1f}; u_del is identically zero, "
              "so s_delta is left as found")
    return out


PARAM_NAMES = ["m_peak", "colr", "x_s", "x_p", "x_psi", "delta_t"]
