"""A rotation of the coefficient coordinates, from the Fisher information.

All of the hierarchy's stiffness is in the coefficients: of the 226 sampled
coordinates of a thirty-object fit, the seven stiffest directions put 96 to
100 percent of their weight on the 46 global ones, and of those the 36
coefficients carry the whole of it --- whitening them takes the condition
number from 2.4e5 to 7.0e3, where whitening the per-object blocks changes
nothing at all.  The stiffest directions are sign-alternating combinations of
the turning coefficients: the near-null directions of a truncated polynomial
basis seen through an integration, which smooths them away.

The metric is the Fisher information of the forward map, ``J' Sigma^-1 J``,
with the prior precision added.  Not the Hessian: measured at the \\salt\\
start, at the posterior mean and at both optimize modes, the Hessian of this
posterior is indefinite --- it is nowhere locally Gaussian --- whereas
``J' Sigma^-1 J`` is positive semi-definite by construction and needs no mode.

``J`` comes from ``stan/jacobian.stan`` by central differences, so the gauge
applied is the one the fit uses; a preconditioner built from a second copy of
that algebra would be silently wrong.  The rotation is a reparameterization:
``R`` is constant, its Jacobian is constant, and the posterior is untouched.
"""

from __future__ import annotations

import numpy as np

from dgsn.stan.hier_data import n_coef, pack, unpack
from dgsn.stan.run import fixed_param


def _flux(data: dict, raw: np.ndarray, hyper: dict, per_sn: dict, k_psi: int, k_q: int) -> np.ndarray:
    d = {k: v for k, v in data.items() if k not in ("grainsize", "R", "f", "sig_f", "m_bright", "xs_bar")}
    d.update(unpack(raw, k_psi, k_q))
    d.update(alpha_q_m=float(hyper["alpha_q_m"]), alpha_psi_m=float(hyper["alpha_psi_m"]))
    d.update({k: np.asarray(per_sn[k], dtype=float) for k in
              ("m_peak", "colr", "x_s", "x_p", "x_psi", "delta_t")})
    return fixed_param("jacobian", d)["log_flux_out"]


def jacobian(data: dict, raw: np.ndarray, hyper: dict, per_sn: dict, k_psi: int, k_q: int,
             h: float = 1e-3) -> tuple[np.ndarray, np.ndarray]:
    """``d log f_i / d raw_k`` at the given point, and the fluxes there.

    Central differences: ``2 n_coef + 1`` evaluations of the forward map over
    the whole epoch set, which is seconds at the sizes these fits run at."""
    n = n_coef(k_psi, k_q)
    base = _flux(data, raw, hyper, per_sn, k_psi, k_q)
    J = np.empty((len(base), n))
    for k in range(n):
        up, dn = raw.astype(float).copy(), raw.astype(float).copy()
        up[k] += h
        dn[k] -= h
        J[:, k] = (_flux(data, up, hyper, per_sn, k_psi, k_q)
                   - _flux(data, dn, hyper, per_sn, k_psi, k_q)) / (2 * h)
    return J, np.exp(base)


def fisher(data: dict, raw: np.ndarray, hyper: dict, per_sn: dict, k_psi: int, k_q: int,
           h: float = 1e-3) -> np.ndarray:
    """``J' Sigma^-1 J`` plus the prior precision, over the coefficients.

    The likelihood is Gaussian in flux, so with ``J`` the derivative of the
    *log* flux the weight of epoch ``i`` is ``(f_i / sigma_i)^2``.  The
    variance carries the intrinsic scatter of Eq. (like), as the fit does."""
    J, f = jacobian(data, raw, hyper, per_sn, k_psi, k_q, h)
    sig_f = np.asarray(data["sig_f"], dtype=float)
    band = np.asarray(data["band"], dtype=int)
    sb = np.where(band == 1, float(hyper["sigma_g"]), float(hyper["sigma_r"]))
    var = sig_f ** 2 + (0.4 * np.log(10) * sb * f) ** 2
    W = f ** 2 / var
    F = J.T @ (W[:, None] * J)

    # prior precision: std_normal on z_a, z_c; Eq. (coefprior) at the pivot on
    # b_til, d_til, which is the only part of the prior that is not unit
    kq, kp = np.arange(1, k_q), np.arange(1, k_psi)
    wq = float(hyper["tau_q"]) * (kq / np.exp(np.log(kq).mean())) ** -(float(hyper["alpha_q_m"]) + 0.5)
    wp = float(hyper["tau_psi"]) * (kp / np.exp(np.log(kp).mean())) ** -(float(hyper["alpha_psi_m"]) + 0.5)
    prior = pack(np.ones(k_q - 1), wq ** -2, np.ones(k_psi - 1), wp ** -2)
    return F + np.diag(prior)


def rotation(F: np.ndarray) -> np.ndarray:
    """``R`` such that coordinates ``w`` with ``raw = R w`` have unit metric.

    With ``F = L L'``, taking ``R = L'^-1`` gives ``R' F R = I``."""
    L = np.linalg.cholesky(F)
    return np.linalg.inv(L.T)


def build(data: dict, inits: dict, k_psi: int, k_q: int, h: float = 1e-3) -> np.ndarray:
    """The rotation at a start point, from that point's own ``inits`` dict."""
    raw = np.asarray(inits["w_coef"], dtype=float)      # identity rotation at the start
    hyper = {k: inits[k] for k in ("alpha_q_m", "alpha_psi_m", "tau_q", "tau_psi",
                                   "sigma_g", "sigma_r")}
    per_sn = {k: inits[k] for k in ("m_peak", "colr", "x_s")}
    per_sn["x_p"] = inits["mu_p"] + np.asarray(inits["u_p"], dtype=float)
    per_sn["x_psi"] = inits["mu_psi"] + np.asarray(inits["u_psi"], dtype=float)
    per_sn["delta_t"] = inits["mu_delta"] + inits["s_delta"] * np.asarray(inits["u_del"], dtype=float)
    return rotation(fisher(data, raw, hyper, per_sn, k_psi, k_q, h))


def point_from_run(out_dir, k_psi: int, k_q: int, n_sn: int) -> tuple[np.ndarray, dict, dict]:
    """The posterior mean of a finished run, as the pieces this module wants.

    The rotation must be built where the posterior is, not where the chain
    starts: built at the \\salt\\ start, whose log density is ten thousand
    below the bulk, it made the conditioning worse by seven orders of
    magnitude; built at the posterior mean of the same model it takes the
    whole posterior from 485 steps to 88."""
    import pandas as pd
    from pathlib import Path

    frames = []
    for csv in sorted(Path(out_dir).glob("hier-?.csv")):
        hdr, rows = None, []
        for line in open(csv):
            if line.startswith("#"):
                continue
            if hdr is None:
                hdr = line.strip().split(",")
                continue
            p = line.strip().split(",")
            if len(p) == len(hdr):
                rows.append(p)
        frames.append(pd.DataFrame(np.asarray(rows[:-1], dtype=float), columns=hdr))
    mean = pd.concat(frames, ignore_index=True).mean()
    raw = pack([mean[f"z_a.{i}"] for i in range(1, k_q)],
               [mean[f"b_til.{i}"] for i in range(1, k_q)],
               [mean[f"z_c.{i}"] for i in range(1, k_psi)],
               [mean[f"d_til.{i}"] for i in range(1, k_psi)])
    hyper = {k: float(mean[k]) for k in
             ("alpha_q_m", "alpha_psi_m", "tau_q", "tau_psi", "sigma_g", "sigma_r")}
    per_sn = {v: np.array([mean[f"{v}.{n}"] for n in range(1, n_sn + 1)])
              for v in ("m_peak", "colr", "x_s", "x_p", "x_psi", "delta_t")}
    return raw, hyper, per_sn


def source_size(out_dir) -> int:
    """How many supernovae the run in ``out_dir`` held."""
    from pathlib import Path

    csv = sorted(Path(out_dir).glob("hier-?.csv"))[0]
    for line in open(csv):
        if line.startswith("#"):
            continue
        names = line.strip().split(",")
        return sum(1 for c in names if c.startswith("m_peak."))
    raise ValueError(f"no header in {out_dir}")


def from_run(out_dir, data: dict, k_psi: int, k_q: int, per_sn: dict | None = None) -> np.ndarray:
    """``R`` built at the posterior mean of the run in ``out_dir``.

    A rotation may be reused across sample sizes: the coefficients it turns
    are the same 36 however many supernovae there are, and the Fisher
    information over them grows close to proportionally with the sample, a
    uniform scaling that Stan's diagonal metric absorbs by itself.  What
    might not carry is the *shape*, more objects constraining different
    combinations, and that is a question for measurement.

    When the source is smaller than the fit being built, its globals are
    used and the per-supernova point must be supplied by the caller --- the
    source knows nothing of the objects it never saw."""
    n_src = source_size(out_dir)
    raw, hyper, src_per = point_from_run(out_dir, k_psi, k_q, n_src)
    n = int(data["N_sn"])
    if n != n_src:
        if per_sn is None:
            raise ValueError(
                f"{out_dir} holds {n_src} supernovae and this fit has {n}; "
                "pass per_sn for the objects the source never saw")
        src_per = per_sn
    return rotation(fisher(data, raw, hyper, src_per, k_psi, k_q))
