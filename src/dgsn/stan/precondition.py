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

Where it is evaluated decides whether it helps.  Built at the \\salt\\ start,
ten thousand nats below the bulk, it made the conditioning worse by seven
orders of magnitude; built at the joint mode, where the variances collapse, it
measured worse than no rotation at all.  Built at a previous run's posterior
mean it took the posterior from 485 steps to 88.

``J`` comes from ``stan/jacobian.stan`` by central differences, so the gauge
applied is the one the fit uses; a preconditioner built from a second copy of
that algebra would be silently wrong.  The rotation is a reparameterization:
``R`` is constant, its Jacobian is constant, and the posterior is untouched.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from dgsn.stan.hier_data import n_coef, pack, unpack
from dgsn.stan.run import fixed_param, log_prob

HYPER = ("alpha_q_m", "alpha_psi_m", "tau_q", "tau_psi", "sigma_g", "sigma_r")


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


def pieces(point: dict, k_psi: int, k_q: int, R_src=None) -> tuple[np.ndarray, dict, dict]:
    """``raw``, the hyperparameters and the per-supernova point of a full
    parameter point, in the form ``fisher`` wants.

    ``R_src`` is the rotation the point was sampled under, since ``w_coef``
    are rotated coordinates and ``raw = R w_coef``; ``None`` means identity."""
    w = np.atleast_1d(np.asarray(point["w_coef"], dtype=float))
    raw = w if R_src is None else np.asarray(R_src, dtype=float) @ w
    hyper = {k: float(point[k]) for k in HYPER}
    per = {k: np.asarray(point[k], dtype=float) for k in ("m_peak", "colr", "x_s")}
    per["x_p"] = float(point["mu_p"]) + np.asarray(point["u_p"], dtype=float)
    per["x_psi"] = float(point["mu_psi"]) + np.asarray(point["u_psi"], dtype=float)
    per["delta_t"] = (float(point["mu_delta"])
                      + float(point["s_delta"]) * np.asarray(point["u_del"], dtype=float))
    return raw, hyper, per


def posterior_mean(out_dir) -> dict:
    """The posterior mean of every sampled parameter of a finished run.

    Pooled over chains, which is sound even when they sit in opposite sign
    branches: the gauge normalizes what the coefficients carry, so a partial
    cancellation of ``z_c`` is a change of scale and the point stays in the
    bulk --- measured at 17.17 against the run's own 17.16 in ``-2 lp`` per
    epoch, with one chain of four in the opposite branch of ``psi``."""
    from dgsn.stan.run import model

    rows, hdr = [], None
    for csv in sorted(Path(out_dir).glob("hier-?.csv")):
        h, r = None, []
        for line in open(csv):
            if line.startswith("#"):
                continue
            if h is None:
                h = line.strip().split(",")
                continue
            p = line.strip().split(",")
            if len(p) == len(h):
                r.append(p)
        hdr = h
        rows.append(np.asarray(r[:-1], dtype=float))        # the last line may be partial
    if not rows:
        raise ValueError(f"no chain CSVs in {out_dir}")
    mean = np.vstack(rows).mean(0)
    col = {c: i for i, c in enumerate(hdr)}
    out = {}
    for name in model("hier").src_info()["parameters"]:
        if name in col:
            out[name] = float(mean[col[name]])
        else:
            v = [mean[col[k]] for k in (f"{name}.{i}" for i in range(1, len(col) + 1)) if k in col]
            if not v:
                raise KeyError(f"{name} not in {out_dir}")
            out[name] = np.array(v)
    return out


def source_rotation(out_dir, k_psi: int, k_q: int) -> np.ndarray | None:
    """The rotation the run in ``out_dir`` was sampled under, or ``None``."""
    R = json.load(open(Path(out_dir) / "data.json")).get("R")
    if R is None:
        return None
    R = np.asarray(R, dtype=float)
    return None if np.allclose(R, np.eye(len(R))) else R


def carry(out_dir, data: dict, k_psi: int, k_q: int) -> tuple[dict, float, float]:
    """A finished run's point moved into the gauge span of ``data``.

    The span enters the model only through ``G_s``, and there only twice: it
    rescales ``c``, which is normalized through ``G_s``, and it changes the
    multiple of ``c`` that the orthogonalization takes out of ``d``.  Both
    displacements lie along ``c`` --- along ``Psi_1`` --- and the turning is
    ``psi = Psibar(d) + x_psi Psi_1(c)``, so the affine change
    ``x_psi -> lambda x_psi + mu`` absorbs both exactly.  Fitting those two
    numbers against Stan's own log density recovers the whole displacement:
    moving the span from (-7, 7.5) to (-6, 5) cost 5353 nats, of which an
    affine carry at ``lambda`` 1.040, ``mu`` +1.494 recovered 101 percent.

    Python searches two numbers here and evaluates nothing: every log density
    is Stan's.  When the spans agree the carry is the identity and no search
    is run."""
    from scipy.optimize import minimize

    point = posterior_mean(out_dir)
    src = json.load(open(Path(out_dir) / "data.json"))
    if (float(src["gsmin"]), float(src["gsmax"])) == (float(data["gsmin"]), float(data["gsmax"])):
        return point, 1.0, 0.0

    mu_psi, u_psi = float(point["mu_psi"]), np.asarray(point["u_psi"], dtype=float)

    def shifted(lam, mu):
        return dict(point, mu_psi=lam * mu_psi + mu, u_psi=lam * u_psi)

    res = minimize(lambda t: -log_prob("hier", data, shifted(*t))[0], [1.0, 0.0],
                   method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-2, maxiter=200))
    lam, mu = float(res.x[0]), float(res.x[1])
    return shifted(lam, mu), lam, mu


def build(data: dict, inits: dict, k_psi: int, k_q: int, h: float = 1e-3) -> np.ndarray:
    """The rotation at a start point, from that point's own ``inits`` dict."""
    raw, hyper, per_sn = pieces(inits, k_psi, k_q)      # identity rotation at the start
    return rotation(fisher(data, raw, hyper, per_sn, k_psi, k_q, h))


def from_optimize(data: dict, inits: dict, k_psi: int, k_q: int, name: str = "hier") -> np.ndarray:
    """``R`` built at the mode of the unadjusted density, found by ``optimize``.

    Kept as a measurement, not as a route.  The argument for it was that a
    rotation needs only a point where the model fits the data, which
    ``optimize`` reaches in minutes where a pilot fit takes hours.  Measured,
    it is worse than no rotation at all: 4095 leapfrogs a draw against 2047
    unpreconditioned, and 48 percent of transitions divergent.

    The reason is not that the variances collapse there, which is what this
    docstring said until it was measured: ``sigma_g`` and ``sigma_r`` land
    within a fifth of their posteriors and the taus within a factor of 1.7.
    What runs away is ``s_delta``, at 5.48 against a posterior 0.186.  It is
    the one non-centred pair left --- the gauge removes ``a``'s and ``c``'s
    scale, so the taus are informed by the centred ``b_til`` and ``d_til`` ---
    and scaling ``s_delta`` by any factor while dividing ``u_del`` by it
    leaves ``delta_t``, and so the likelihood, untouched while the prior
    ``-1/2 sum u_del^2`` rises toward zero.  Measured along that ridge the log
    density gains 7.3 nats monotonically and plateaus by a factor of 30, which
    is where the optimizer stopped: the constrained-scale mode is not finite
    in ``s_delta``.  A Fisher built there mis-scales the per-supernova timing
    block by thirty.

    Use ``from_run``, which carries a previous posterior across a change of
    basis.  To start a *sampler* from this mode instead, re-place it on the
    ridge first --- ``s_delta`` to sd(``delta_t`` - ``mu_delta``) and ``u_del``
    to match --- which keeps the fit and fixes the scaling.
    """
    from dgsn.stan.run import optimize

    mode = optimize(name, data, inits, jacobian=False)
    n = int(data["N_sn"])
    point = {k: float(mode[k]) for k in HYPER + ("mu_p", "mu_psi", "mu_delta", "s_delta")}
    point["w_coef"] = np.array([mode[f"w_coef[{i}]"] for i in range(1, n_coef(k_psi, k_q) + 1)])
    for v in ("m_peak", "colr", "x_s", "u_p", "u_psi", "u_del"):
        point[v] = np.array([mode[f"{v}[{i}]"] for i in range(1, n + 1)])
    raw, hyper, per_sn = pieces(point, k_psi, k_q)
    return rotation(fisher(data, raw, hyper, per_sn, k_psi, k_q))


def source_size(out_dir) -> int:
    """How many supernovae the run in ``out_dir`` held."""
    csv = sorted(Path(out_dir).glob("hier-?.csv"))[0]
    for line in open(csv):
        if line.startswith("#"):
            continue
        return sum(1 for c in line.strip().split(",") if c.startswith("m_peak."))
    raise ValueError(f"no header in {out_dir}")


def from_run(out_dir, data: dict, k_psi: int, k_q: int, per_sn: dict | None = None,
             verbose: bool = True) -> np.ndarray:
    """``R`` built at the posterior mean of the run in ``out_dir``.

    A rotation may be reused across sample sizes: the coefficients it turns
    are the same 36 however many supernovae there are, and the Fisher
    information over them grows close to proportionally with the sample, a
    uniform scaling that Stan's diagonal metric absorbs by itself.  What
    might not carry is the *shape*, more objects constraining different
    combinations, and that is a question for measurement.

    A change of gauge span is carried by ``carry``.  A change of sample is
    not: the source knows nothing of the objects it never saw, so ``per_sn``
    must then be supplied, and because the carry touches only ``x_psi`` ---
    which those objects take from their own start, already in the new span ---
    a supplied ``per_sn`` needs no correction."""
    point, lam, mu = carry(out_dir, data, k_psi, k_q)
    raw, hyper, src_per = pieces(point, k_psi, k_q, source_rotation(out_dir, k_psi, k_q))
    n, n_src = int(data["N_sn"]), source_size(out_dir)
    if per_sn is not None:
        src_per = {k: np.asarray(v, dtype=float) for k, v in per_sn.items()}
    elif n != n_src:
        raise ValueError(
            f"{out_dir} holds {n_src} supernovae and this fit has {n}; "
            "pass per_sn for the objects the source never saw")
    if verbose and (lam, mu) != (1.0, 0.0):
        lp = log_prob("hier", data, point)[0]
        print(f"carried across the gauge span: lambda {lam:.4f}, mu {mu:+.4f}; "
              f"lp at the carried point {lp:.1f} ({-2 * lp / int(data['N']):.2f} per epoch)")
    return rotation(fisher(data, raw, hyper, src_per, k_psi, k_q))
