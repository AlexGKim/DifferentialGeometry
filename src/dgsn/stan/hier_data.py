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
DELTA_BOUND = 3.0       # d; must match <lower=-3, upper=3> on mu_delta and the
                        # uniform wall on d_dev in the hierarchies


def _v_span() -> float:
    """``vhi - vlo`` of Eq. (vmap), which sets P_0 = 1/sqrt(vhi - vlo).

    Stan builds this in ``transformed data``; the same two ``asinh`` calls are
    repeated here only to place ``b_0``, never to evaluate the basis."""
    return float(np.arcsinh(D.MODEL_WINDOW[1] / D.DELTA_P)
                 - np.arcsinh(D.MODEL_WINDOW[0] / D.DELTA_P))


def build(names, k_psi: int, k_q: int, grainsize: int = 1, m_per_side: int = D.M_PER_SIDE,
          sample=None, R=None, mode: bool = True, trav: bool = True) -> dict:
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
        R=np.eye(n_coef(k_psi, k_q, mode, trav)) if R is None else np.asarray(R, dtype=float),
    )
    return data


def n_coef(k_psi: int, k_q: int, mode: bool = True, trav: bool = True) -> int:
    """The rotated block: z_a, b_til, z_c, d_til laid end to end, dropping
    z_c when the turning mode is absent and z_a when the traversal mode is.

    ``mode`` is Psi_1 of Eq. (curvlin), ``trav`` is q_1 of Eq. (travlog)."""
    return (2 if trav else 1) * (k_q - 1) + (2 if mode else 1) * (k_psi - 1)


def pack(*vectors) -> np.ndarray:
    """The coefficient coordinates in the order the program slices them:
    ``z_a, b_til, z_c, d_til``, or ``z_a, b_til, d_til`` without the mode."""
    return np.concatenate([np.asarray(v, dtype=float) for v in vectors])


def unpack(raw, k_psi: int, k_q: int, mode: bool = True, trav: bool = True) -> dict:
    n1, n3 = k_q - 1, k_psi - 1
    raw = np.asarray(raw, dtype=float)
    if not trav:                       # b_til, d_til alone
        return dict(b_til=raw[:n1], d_til=raw[n1:n1+n3])
    if not mode:
        return dict(z_a=raw[:n1], b_til=raw[n1:2*n1], d_til=raw[2*n1:2*n1+n3])
    return dict(z_a=raw[:n1], b_til=raw[n1:2*n1], z_c=raw[2*n1:2*n1+n3], d_til=raw[2*n1+n3:])


def inits(names, k_psi: int, k_q: int, x1_ref: float, sample=None, globals_=None,
          tau: dict | None = None, R=None, mode: bool = True,
          trav: bool = True) -> dict:
    """Initial values: per-supernova from SALT2 and the sample, global
    functions from the SALT2 fits.  The centred pair starts at the coefficients
    themselves; ``z`` is put on unit scale by the same power law.  ``tau`` is
    no longer an initial value --- Eq. (ccoef)'s scales are ``data``, fixed at
    ``D.TAU_Q`` and ``D.TAU_PSI`` --- and the argument is kept only so that a
    caller can still override the weights used for ``a`` and ``c``."""
    sample = analysis_sample() if sample is None else sample
    G = globals_ if globals_ is not None else fit_globals(k_psi, k_q, x1_ref)
    kq, kp = np.arange(1, k_q), np.arange(1, k_psi)
    # the pivot of hier.stan: tau is the scale at the geometric-mean mode
    wq = (kq / np.exp(np.log(kq).mean())) ** -INIT_ALPHA
    wp = (kp / np.exp(np.log(kp).mean())) ** -INIT_ALPHA
    if tau is None:
        # the measured scales of D.TAU_Q, D.TAU_PSI; starting there costs the
        # chain nothing and is where a fit at this K has been seen to land
        tau = dict(q=float(D.TAU_Q), psi=float(D.TAU_PSI))
    per = [per_sn_start(sample.loc[n]) for n in names]
    colr0 = np.array([q["colr"] for q in per])
    if not trav:
        raw = pack(G["b"], G["d"])     # no z_a, no z_c
    elif mode:
        raw = pack(G["a"] / wq, G["b"], G["c"] / wp, G["d"])
    else:
        raw = pack(G["a"] / wq, G["b"], G["d"])
    out = dict(
        # w = R^-1 raw, so that the start is the same point whatever rotation
        # is in force
        w_coef=raw if R is None else np.linalg.solve(np.asarray(R, dtype=float), raw),
        tau_q=tau["q"], tau_psi=tau["psi"],
        alpha_q_m=INIT_ALPHA - 0.5, alpha_psi_m=INIT_ALPHA - 0.5,
        m_peak=np.array([p["m_peak"] for p in per]),
        # the anchor colour as a named mean and deviations about it; the
        # start is the same point in colr either way
        c_bar=float(colr0.mean()),
        c_dev=colr0 - colr0.mean(),
        # the mean log rate is q_0's k = 0 coefficient now, so the SALT2 value
        # enters as b_0 = sqrt(vhi - vlo) * x_s rather than as a named mean
        b_0=float(G["x_s"]) * np.sqrt(_v_span()),
        x_s=np.zeros(len(names)),      # sums to zero, as sum_to_zero_vector requires
        mu_delta=0.0,                  # the near gauge, named so the metric can
        d_dev=np.zeros(len(names)),    # scale it; deviations sum to zero
    )
    if mode:
        out.update(mu_psi=G["x_psi"], u_psi=np.zeros(len(names)))
    if trav:
        out.update(mu_p=G["x_p"], u_p=np.zeros(len(names)))
    return out


def _from_mode(mode: dict, name: str):
    """One parameter out of an ``optimize`` result, scalar or vector.

    CmdStan reports ``alpha_q_m`` but ``d_dev[1]``, so the shape has to be
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
                        mode: bool = True, trav: bool = True, verbose: bool = True,
                        polish: bool = True) -> dict:
    """Initial values from Stan's own optimizer on the data.

    Why this exists.  Section 5 starts the chains at a \\salt\\ fit, which the
    third standing commitment permits as its one exemption, licensed by
    dispersed-start agreement.  A start found by optimizing the model's own
    log density on the observed epochs needs no exemption: the only template
    left is the one that starts the optimizer, a dependence on a basin rather
    than on a point.

    What used to be here.  The mode of the previous parameterization was not
    finite in ``s_delta``: scaling it while dividing ``u_del`` by the same
    factor left ``delta_t``, and so the likelihood, exactly unchanged while
    the prior ``-1/2 sum u_del^2`` rose toward zero, so the optimizer walked
    up that ridge and returned ``s_delta`` 5.48 against a posterior 0.186.
    Every start had to be re-placed along it.  The ridge was a ridge in the
    scale of a population of ``delta_t``, and that population is gone ---
    ``delta_t`` is now a ``sum_to_zero_vector`` with no scale of its own ---
    so the re-placement has gone with it.

    This voids the reason ``polish`` defaults to False at the call sites.
    That default was not about cost: ``optimize``'s restart loop could not
    terminate in ``s_delta``, each restart walking further up the ridge by
    fractions of a nat for about three minutes, and the re-placement then
    discarded the distance gained.  With the ridge removed the restarts may
    well converge, and whether they earn their time is now an open question
    to be measured rather than an answered one.

    ``tau_q`` and ``tau_psi`` need no start at all now, being ``data``.
    """
    from dgsn.stan.run import optimize

    m = optimize(program, data, start, jacobian=False, polish=polish)
    # s_bar and x_s_dev, not x_s: the optimizer reports x_s too, as a
    # transformed parameter, but an inits dict has to name what is sampled.
    keys = ["w_coef", "tau_q", "tau_psi", "alpha_q_m", "alpha_psi_m",
            "b_0", "x_s", "c_bar", "c_dev", "mu_delta", "d_dev", "m_peak"]
    if mode:
        keys += ["mu_psi", "u_psi"]
    if trav:
        keys += ["mu_p", "u_p"]
    out = {k: _from_mode(m, k) for k in keys}
    # Back onto the constraint surfaces after the roundtrip through text: the
    # mode is written at 18 significant figures, whose sum is no longer
    # exactly zero, and sum_to_zero_vector requires that it is.
    for k in ("x_s", "c_dev", "d_dev"):
        v = np.asarray(out[k], dtype=float)
        out[k] = v - v.mean()
    out["mu_delta"] = float(np.clip(out["mu_delta"],
                                    -DELTA_BOUND + 1e-6, DELTA_BOUND - 1e-6))
    if verbose:
        print(f"optimize start: lp {float(m['lp__']):.1f}")
    return out


def inits_from_run(out_dir, names, k_psi: int, k_q: int, x1_ref: float,
                   mode: bool = True, trav: bool = True, sample=None, globals_=None,
                   R=None, program: str | None = None, data: dict | None = None,
                   use_mean: bool = False, verbose: bool = True) -> dict:
    """Start a richer model at a simpler one's posterior mean.

    The three hierarchy programs are nested, and the embedding is exact: a
    departure mode contributes to the likelihood only through
    ``x_p q_1`` or ``x_psi Psi_1``, so setting the amplitudes to zero ---
    ``mu_p = u_p = 0``, ``mu_psi = u_psi = 0`` --- reproduces the simpler
    model's prediction whatever the mode vector is.  The mode vector itself
    cannot be set to zero: the gauge normalizes it, and ``at / ||at||`` is
    undefined there, so ``z_a`` and ``z_c`` keep their \\salt\\ starting
    values and are simply not used at the start.

    Why bother.  A converged posterior of a simpler model is a better start
    than the \\salt\\ fit or an optimizer's mode: it is a point a sampler has
    actually occupied, so the scatters, the coefficient scales and every
    per-supernova placement are where the data put them rather than where a
    template or a maximum put them.
    """
    from dgsn.stan.precondition import best_draw, posterior_mean

    src = posterior_mean(out_dir) if use_mean else best_draw(out_dir)
    out = inits(names, k_psi, k_q, x1_ref, sample=sample, globals_=globals_, R=R,
                mode=mode, trav=trav)
    # carried across: everything the simpler model also had
    carried = []
    for k in list(out):
        if k in src and k != "w_coef":
            v = np.asarray(src[k], dtype=float)
            if np.shape(v) == np.shape(np.asarray(out[k], dtype=float)):
                out[k] = v if np.ndim(out[k]) else float(v)
                carried.append(k)
    # A source that predates the constrained placements carries colr whole.
    # Split it, or c_bar and c_dev keep their \salt\ values while everything
    # around them comes from the posterior.
    if "colr" in src and "c_bar" in out:
        c = np.asarray(src["colr"], dtype=float)
        if c.shape == np.shape(out["c_dev"]):
            out["c_bar"], out["c_dev"] = float(c.mean()), c - c.mean()
            carried.append("colr -> c_bar, c_dev")
    for k in ("x_s", "c_dev", "d_dev"):
        if k in out:
            out[k] = np.asarray(out[k], dtype=float) - np.mean(out[k])
    # delta_t is declared on a box now, not on the sum-to-zero hyperplane, so
    # it is clipped rather than projected: a source drawn under a wider prior
    # can sit outside the bound, and a start outside its own support is
    # rejected.  The largest |Delta| ever drawn is 2.43 d, so this normally
    # does nothing.
    # a source that predates the split carries delta_t whole: separate it
    if "delta_t" in src and "mu_delta" in out:
        dt = np.asarray(src["delta_t"], dtype=float)
        if dt.shape == np.shape(out["d_dev"]):
            out["mu_delta"] = float(np.clip(dt.mean(), -DELTA_BOUND + 1e-6,
                                            DELTA_BOUND - 1e-6))
            out["d_dev"] = dt - dt.mean()
            carried.append("delta_t as mu_delta + d_dev")

    # the coefficient block, by name, so a shorter source lands in the right slots
    # From the relabeled blocks, not from w_coef: w_coef is the sampled vector
    # and the relabeling cannot act on it, R mixing its slots, so a point
    # assembled from w_coef beside relabeled mu_p and u_p has a's sign against
    # x_p's.  That inverted the traversal and cost 38109 nats.
    s_un = {k: np.asarray(src[k], dtype=float)
            for k in ("z_a", "b_til", "z_c", "d_til") if k in src}
    if not s_un:
        s_un = unpack(np.asarray(src["w_coef"], dtype=float), k_psi, k_q,
                      mode="z_c" in _slots(out_dir), trav="z_a" in _slots(out_dir))
    d_un = unpack(np.asarray(out["w_coef"], dtype=float), k_psi, k_q, mode, trav)
    for key in d_un:
        if key in s_un:
            d_un[key] = s_un[key]
            carried.append(key)
    raw = pack(*[d_un[k] for k in ("z_a", "b_til", "z_c", "d_til") if k in d_un])
    out["w_coef"] = raw if R is None else np.linalg.solve(np.asarray(R, dtype=float), raw)
    # A source from before Eq. (xsmean) carries x_s itself; split it the way
    # the gauge now does, which is exactly how the constraint defines it.
    if "x_s" in src and "b_0" in out and "b_0" not in src:
        xs = np.asarray(src["x_s"], dtype=float)
        out["b_0"] = float(xs.mean()) * np.sqrt(_v_span())
        out["x_s"] = xs - xs.mean()
        carried.append("x_s as b_0 + x_s")
    # the new modes start unused
    zeroed = [k for k in ("mu_p", "mu_psi") if k in out and k not in src]
    for k in zeroed:
        out[k] = 0.0
        out["u" + ("_p" if k == "mu_p" else "_psi")] = np.zeros(len(names))
    if verbose:
        msg = (f"started from {out_dir}: carried {len(set(carried))} blocks"
               + (f", new modes zeroed: {', '.join(zeroed)}" if zeroed else ""))
        # A count of blocks says nothing about fidelity: a transfer that
        # silently drops x_s reports nineteen blocks carried and lands fifteen
        # thousand nats below the Salt2 start.  Report the log density, which
        # cannot be fooled that way.
        if program is not None and data is not None:
            from dgsn.stan.run import log_prob
            msg += f"; lp at the transferred point {log_prob(program, data, out)[0]:.1f}"
        print(msg)
    return out


def _slots(out_dir) -> set:
    """Which coefficient blocks the run in ``out_dir`` sampled."""
    from dgsn.stan.run import chain_csvs

    for line in open(chain_csvs(out_dir)[0]):
        if line.startswith("#"):
            continue
        cols = line.strip().split(",")
        return {n for n in ("z_a", "b_til", "z_c", "d_til")
                if any(c.startswith(n + ".") for c in cols)}
    raise ValueError(f"no header yet in {out_dir}: Stan writes its configuration "
                     "comments before the column names, so a run still in early warmup "
                     "has no header to read")


PARAM_NAMES = ["m_peak", "colr", "x_s", "x_p", "x_psi", "delta_t"]
