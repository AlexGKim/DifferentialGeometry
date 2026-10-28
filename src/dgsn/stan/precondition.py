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

from dgsn.stan.hier_data import _v_span, n_coef, pack, unpack
from dgsn.stan.run import chain_csvs, fixed_param, log_prob

HYPER = ("alpha_q_m", "alpha_psi_m", "tau_q", "tau_psi")


def _flux(data: dict, raw: np.ndarray, hyper: dict, per_sn: dict, k_psi: int, k_q: int,
          mode: bool = True, trav: bool = True) -> np.ndarray:
    d = {k: v for k, v in data.items() if k not in ("grainsize", "R", "f", "sig_f", "m_bright", "xs_bar")}
    parts = unpack(raw, k_psi, k_q, mode, trav)
    # jacobian.stan declares all four blocks; the absent ones are passed as
    # zeros and ignored, has_mode telling it the gauge does not apply
    d.update(z_a=parts.get("z_a", np.zeros(k_q - 1)), b_til=parts["b_til"],
             z_c=parts.get("z_c", np.zeros(k_psi - 1)), d_til=parts["d_til"],
             has_mode=int(bool(mode and trav)))
    d.update(alpha_q_m=float(hyper["alpha_q_m"]), alpha_psi_m=float(hyper["alpha_psi_m"]))
    d.update({k: np.asarray(per_sn[k], dtype=float) for k in
              ("m_peak", "colr", "x_s", "x_p", "x_psi", "delta_t")})
    return fixed_param("jacobian", d)["log_flux_out"]


def jacobian(data: dict, raw: np.ndarray, hyper: dict, per_sn: dict, k_psi: int, k_q: int,
             h: float = 1e-3, mode: bool = True, trav: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """``d log f_i / d raw_k`` at the given point, and the fluxes there.

    Central differences: ``2 n_coef + 1`` evaluations of the forward map over
    the whole epoch set, which is seconds at the sizes these fits run at."""
    n = n_coef(k_psi, k_q, mode, trav)
    base = _flux(data, raw, hyper, per_sn, k_psi, k_q, mode, trav)
    J = np.empty((len(base), n))
    for k in range(n):
        up, dn = raw.astype(float).copy(), raw.astype(float).copy()
        up[k] += h
        dn[k] -= h
        J[:, k] = (_flux(data, up, hyper, per_sn, k_psi, k_q, mode, trav)
                   - _flux(data, dn, hyper, per_sn, k_psi, k_q, mode, trav)) / (2 * h)
    return J, np.exp(base)


def fisher(data: dict, raw: np.ndarray, hyper: dict, per_sn: dict, k_psi: int, k_q: int,
           h: float = 1e-3, mode: bool = True, trav: bool = True) -> np.ndarray:
    """``J' Sigma^-1 J`` plus the prior precision, over the coefficients.

    The likelihood is Gaussian in flux, so with ``J`` the derivative of the
    *log* flux the weight of epoch ``i`` is ``(f_i / sigma_{f,i})^2``.  The
    variance is the archive's error alone, as the fit's now is: Eq. (like)
    carried a second, fitted term and no longer does."""
    J, f = jacobian(data, raw, hyper, per_sn, k_psi, k_q, h, mode, trav)
    sig_f = np.asarray(data["sig_f"], dtype=float)
    W = (f / sig_f) ** 2
    F = J.T @ (W[:, None] * J)

    # prior precision: std_normal on z_a, z_c; Eq. (coefprior) at the pivot on
    # b_til, d_til, which is the only part of the prior that is not unit
    kq, kp = np.arange(1, k_q), np.arange(1, k_psi)
    # sampled in the hierarchies, but a program that holds them as data puts
    # them in the data dict instead; take whichever is there
    tq = float(hyper["tau_q"]) if "tau_q" in hyper else float(data["tau_q"])
    tp = float(hyper["tau_psi"]) if "tau_psi" in hyper else float(data["tau_psi"])
    wq = tq * (kq / np.exp(np.log(kq).mean())) ** -(float(hyper["alpha_q_m"]) + 0.5)
    wp = tp * (kp / np.exp(np.log(kp).mean())) ** -(float(hyper["alpha_psi_m"]) + 0.5)
    prior = (pack(np.ones(k_q - 1), wq ** -2, np.ones(k_psi - 1), wp ** -2)
             if mode and trav else pack(wq ** -2, wp ** -2))
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
    # Prefer the relabeled blocks.  w_coef is the sampled vector and is *not*
    # relabeled --- it cannot be, since R mixes its slots --- so a point built
    # from w_coef together with relabeled mu_p and u_p has a's sign against
    # x_p's, inverting the traversal.  Measured, that cost 38109 nats of
    # likelihood at a point that should have transferred exactly.
    blocks = [k for k in ("z_a", "b_til", "z_c", "d_til") if k in point]
    if blocks:
        raw = pack(*[point[k] for k in blocks])
    else:
        w = np.atleast_1d(np.asarray(point["w_coef"], dtype=float))
        raw = w if R_src is None else np.asarray(R_src, dtype=float) @ w
    hyper = {k: float(point[k]) for k in HYPER if k in point}
    per = {"m_peak": np.asarray(point["m_peak"], dtype=float),
           "delta_t": (np.asarray(point["delta_t"], dtype=float)
                       if "delta_t" in point
                       else float(point["mu_delta"])
                            + np.asarray(point["d_dev"], dtype=float))}
    # x_s_tot and colr are transformed parameters in a finished run's output and
    # have to be assembled from their pieces in an inits dict.  For x_s it is
    # the total that reaches the forward map, x_s + b_0 P_0, not the sampled
    # deviation: taking the deviation would evaluate the Fisher information at
    # a point the model never visits, off by the whole mean log rate.
    per["x_s"] = (np.asarray(point["x_s_tot"], dtype=float) if "x_s_tot" in point
                  else np.asarray(point["x_s"], dtype=float)
                       + float(point["b_0"]) / np.sqrt(_v_span()))
    per["colr"] = (np.asarray(point["colr"], dtype=float) if "colr" in point
                   else float(point["c_bar"]) + np.asarray(point["c_dev"], dtype=float))
    # The departure amplitudes exist only where the modes do.  hier_common has
    # neither, and passes literal zeros to the forward map in their place.
    nsn = len(per["m_peak"])
    per["x_p"] = (float(point["mu_p"]) + np.asarray(point["u_p"], dtype=float)
                  if "mu_p" in point else np.zeros(nsn))
    per["x_psi"] = (float(point["mu_psi"]) + np.asarray(point["u_psi"], dtype=float)
                    if "mu_psi" in point else np.zeros(nsn))
    return raw, hyper, per


def _draws(out_dir) -> tuple[np.ndarray, list]:
    """Every post-warmup draw of a run, relabeled, with the column names.

    Section 2.3 leaves each pair a sign and ``hier.stan`` reports which branch
    a draw is in.  Relabeling here, once, is what makes anything downstream
    able to pool: without it a run with two chains in each branch cancels its
    signed blocks, ``||mean(a)||`` falling to 0.071 where every draw has
    ``||a|| = 1``, and a rotation or a transfer built on that mean lands
    fifteen thousand nats below where the run itself sat."""
    rows, hdr = [], None
    for csv in chain_csvs(out_dir):
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
        rows.append(np.asarray(r[:-1], dtype=float))     # the last line may be partial
    if not rows:
        raise ValueError(f"no chain CSVs in {out_dir}")
    X = np.vstack(rows)
    col = {c: i for i, c in enumerate(hdr)}
    for flag, signed in (("flip_q", ("z_a", "mu_p", "u_p")),
                         ("flip_psi", ("z_c", "mu_psi", "u_psi"))):
        if flag not in col:
            continue
        flip = X[:, col[flag]] > 0.5
        if not flip.any():
            continue
        idx = [i for nm in signed for c, i in col.items()
               if c == nm or c.startswith(nm + ".")]
        X[np.ix_(flip, idx)] *= -1
    return X, hdr


def _named(row, col) -> dict:
    """One row of the draw table as a parameter dict."""
    from dgsn.stan.run import model

    out = {}
    # z_a, b_til, z_c, d_til are harvested alongside w_coef because the sign
    # relabeling acts on *them*: negating a slot of w_coef would be wrong
    # whenever R is not the identity, since the rotation mixes the slots.
    # Reading the relabeled blocks and re-packing is correct either way.
    for name in list(model("hier").src_info()["parameters"]) + \
            ["x_s", "x_p", "x_psi", "delta_t", "z_a", "b_til", "z_c", "d_til"]:
        if name in col:
            out[name] = float(row[col[name]])
            continue
        v = [row[col[k]] for k in (f"{name}.{i}" for i in range(1, len(col) + 1)) if k in col]
        if v:
            out[name] = np.array(v)
    return out


def best_draw(out_dir) -> dict:
    """The highest-density draw of a finished run, relabeled.

    For a starting point this is preferable to the mean of the draws, and the
    reason is not subtle: the model is nonlinear, so a parameter-wise average
    need not be a point the model fits at all.  Measured on an unconverged run
    --- R-hat 1.44, effective size 8 --- the mean sat at ``-2 lp`` per epoch
    32.3 where every draw in it sat near 17.0.  A draw cannot fail that way.
    """
    X, hdr = _draws(out_dir)
    col = {c: i for i, c in enumerate(hdr)}
    return _named(X[int(np.argmax(X[:, col["lp__"]]))], col)


def posterior_mean(out_dir) -> dict:
    """The posterior mean of every parameter of a finished run, relabeled.

    Sound only because ``_draws`` relabels first.  Prefer ``best_draw`` for a
    starting point or a rotation's reference: the model is nonlinear, so a
    parameter-wise mean need not be a point the model fits, and on an
    unconverged run it is not --- ``-2 lp`` per epoch 32.3 at the mean where
    every draw in it sat near 17.0."""
    X, hdr = _draws(out_dir)
    col = {c: i for i, c in enumerate(hdr)}
    return _named(X.mean(0), col)


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
    was where the optimizer stopped: the constrained-scale mode was not
    finite in ``s_delta``, and a Fisher built there mis-scaled the
    per-supernova timing block by thirty.  That ridge went with the
    population of ``delta_t``, which is now a ``sum_to_zero_vector`` with no
    scale of its own, so the warning no longer applies and this has not been
    re-measured since.  ``from_run``, which carries a previous posterior
    across a change of basis, remains the better-tested route.
    """
    from dgsn.stan.run import optimize

    mode = optimize(name, data, inits, jacobian=False)
    n = int(data["N_sn"])
    point = {k: float(mode[k]) for k in HYPER + ("mu_p", "mu_psi", "b_0", "c_bar", "mu_delta")}
    point["w_coef"] = np.array([mode[f"w_coef[{i}]"] for i in range(1, n_coef(k_psi, k_q) + 1)])
    for v in ("m_peak", "x_s", "c_dev", "d_dev", "u_p", "u_psi"):
        point[v] = np.array([mode[f"{v}[{i}]"] for i in range(1, n + 1)])
    raw, hyper, per_sn = pieces(point, k_psi, k_q)
    return rotation(fisher(data, raw, hyper, per_sn, k_psi, k_q))


def source_size(out_dir) -> int:
    """How many supernovae the run in ``out_dir`` held."""
    csv = chain_csvs(out_dir)[0]
    for line in open(csv):
        if line.startswith("#"):
            continue
        return sum(1 for c in line.strip().split(",") if c.startswith("m_peak."))
    raise ValueError(f"no header in {out_dir}")


def from_run(out_dir, data: dict, k_psi: int, k_q: int, per_sn: dict | None = None,
             verbose: bool = True) -> np.ndarray:
    """``R`` built at the posterior mean of the run in ``out_dir``.

    The coefficients it turns are the same however many supernovae there are,
    and the Fisher information over them grows close to proportionally with
    the sample, a uniform scaling that Stan's diagonal metric absorbs by
    itself.  That is the argument for reuse across sample sizes, not a
    measurement of it: ``runs/TRANSFER-VERDICT.txt`` has the rotation built
    where the model fits the data, per-supernova parameters included, so it
    is not a portable artifact, and whether its *shape* survives a change of
    sample size --- more objects constraining different combinations --- was
    still open when this was written.  Across a change of ``K_psi`` or
    ``K_q`` there is no reuse at all; the guard below says why.

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
    # The layout follows the data: hier.stan carries four coefficient blocks,
    # hier_nomode three, hier_common two.  Which one is read off the length,
    # and a length belonging to none of them is an error rather than a guess.
    # unpack slices a prefix, so a source at a larger K used to be accepted
    # silently: a (24, 16) source read at (12, 8) put b_til[7:15] where d_til
    # belongs, filling the turning coefficients from the rate ones, and the
    # run then printed a plausible condition number and passed the gradient
    # check, which tests Stan's gradient against its own log density and not
    # the rotation's provenance.
    layouts = {n_coef(k_psi, k_q, m, t): (m, t)
               for m, t in ((True, True), (False, True), (False, False))}
    if len(raw) not in layouts:
        raise ValueError(
            f"{out_dir} carries {len(raw)} coefficients and this fit wants "
            f"{sorted(layouts)} at K = ({k_psi}, {k_q}); a rotation cannot be "
            "carried across a change of K_psi or K_q, so build the source at "
            "this fit's K or run the fit at the source's")
    mode, trav = layouts[len(raw)]
    return rotation(fisher(data, raw, hyper, src_per, k_psi, k_q, mode=mode, trav=trav))
