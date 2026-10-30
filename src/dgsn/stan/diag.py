"""Split-Rhat and ESS over draws read straight from CmdStan's CSVs.

Why this exists rather than ``fit.summary()``.  A ladder of arms has to be
read while it is still running, and ``from_csv`` wants complete files; these
functions take whatever rows are on disk.  They are also the only way to get
either number for a derived quantity such as ``||d||``, which is not a column.

Vehtari et al. (2021) as CmdStan implements it, rank normalization
included: ``ess_bulk`` normalizes the pooled draws before measuring the
autocorrelation, and ``rhat`` takes the larger of the rank-normalized and
the folded rank-normalized split statistics.  Both are needed to match
``stansummary``.  The plain versions are kept as ``ess`` and ``rhat_plain``
because they are the ones that answer "is this chain stationary in itself",
rank normalization being a device for the pooled estimate.

``scripts/check_diag.py`` is the validation: it runs ``stansummary`` over a
run and compares column by column.
"""

from __future__ import annotations

import numpy as np


def _autocov(x: np.ndarray) -> np.ndarray:
    n = len(x)
    m = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(x - x.mean(), m)
    return np.fft.irfft(f * np.conjugate(f), m)[:n].real / n


def _split(chains: np.ndarray) -> np.ndarray:
    """(M, N) -> (2M, N/2).  Stan measures both statistics on split chains, and
    for ESS it matters: a chain that drifts monotonically has most of its
    autocorrelation removed once each half carries its own mean, which is why
    an unsplit ESS reads about half of CmdStan's on a drifting column."""
    c = np.asarray(chains, float)
    h = c.shape[1] // 2
    return np.concatenate([c[:, :h], c[:, c.shape[1] - h:]], axis=0)


def ess(chains) -> float:
    """``chains`` of shape (M, N).  Geyer's initial-positive sequence on the
    combined autocorrelation of the split chains, with the between-chain
    variance folded in."""
    chains = _split(chains)
    m_chain, n = chains.shape
    if m_chain < 2:
        raise ValueError("ESS needs at least two chains")
    acov = np.stack([_autocov(c) for c in chains])
    chain_var = acov[:, 0] * n / (n - 1.0)
    mean_var = chain_var.mean()
    var_plus = mean_var * (n - 1.0) / n + chains.mean(1).var(ddof=1)
    if var_plus <= 0:                      # a held constant: no variance at all
        return np.nan
    rho = np.zeros(n)
    rho[0] = 1.0
    rho_even, rho_odd = 1.0, 1.0 - (mean_var - acov[:, 1].mean()) / var_plus
    if n > 1:
        rho[1] = rho_odd
    t = 1
    while t + 2 < n - 1 and (rho_even + rho_odd) > 0:
        r_e = 1.0 - (mean_var - acov[:, t + 1].mean()) / var_plus
        r_o = 1.0 - (mean_var - acov[:, t + 2].mean()) / var_plus
        if r_e + r_o <= 0:
            break
        rho_even, rho_odd = r_e, r_o
        rho[t + 1], rho[t + 2] = r_e, r_o
        t += 2
    max_t = t
    t = 1
    while t <= max_t - 3:                  # initial monotone sequence
        if rho[t + 1] + rho[t + 2] > rho[t - 1] + rho[t]:
            rho[t + 1] = rho[t + 2] = (rho[t - 1] + rho[t]) / 2.0
        t += 2
    tau = -1.0 + 2.0 * rho[: max_t + 1].sum()
    if max_t + 1 < n:
        tau += rho[max_t + 1]
    return m_chain * n / max(tau, 1.0 / np.log10(m_chain * n))


def rhat_plain(chains) -> float:
    """Split-Rhat on the draws as they are: each chain halved, so a chain that
    drifts within itself fails on its own."""
    c = _split(chains)
    n = c.shape[1]
    w = c.var(1, ddof=1).mean()
    if w <= 0:
        return np.nan
    return float(np.sqrt(((n - 1) / n * w + c.mean(1).var(ddof=1)) / w))


def _rank_normalize(chains: np.ndarray) -> np.ndarray:
    """Pooled average ranks through the normal quantile, Eq. (14) of Vehtari
    et al. (2021).  This is what makes the statistics below invariant to a
    monotone transform and finite on a distribution without a mean."""
    from scipy.special import ndtri
    from scipy.stats import rankdata
    flat = np.asarray(chains, float).ravel()
    r = rankdata(flat).reshape(np.shape(chains))
    return ndtri((r - 3.0 / 8.0) / (flat.size - 0.25))


def rhat(chains) -> float:
    """CmdStan's published R_hat: the larger of the rank-normalized split-Rhat
    and the same on the folded draws |x - median(x)|, which catches a
    disagreement in scale that the location statistic misses."""
    c = np.asarray(chains, float)
    if c.var() <= 0:
        return np.nan
    bulk = rhat_plain(_rank_normalize(c))
    folded = rhat_plain(_rank_normalize(np.abs(c - np.median(c))))
    return float(np.nanmax([bulk, folded]))


def ess_bulk(chains) -> float:
    """CmdStan's ESS_bulk: ``ess`` of the rank-normalized draws."""
    c = np.asarray(chains, float)
    if c.var() <= 0:
        return np.nan
    return ess(_rank_normalize(c))


def read_csvs(paths, columns=None):
    """``{name: (M, N) array}`` from CmdStan CSVs, truncated to the shortest
    chain so the arrays are rectangular.  Comment lines are skipped, so a file
    a chain is still writing reads as the rows it has."""
    paths = [str(p) for p in paths]
    names = None
    for line in open(paths[0]):
        if line.startswith("lp__"):
            names = line.strip().split(",")
            break
    if names is None:
        raise ValueError(f"no header in {paths[0]}")
    keep = names if columns is None else [c for c in columns if c in names]
    idx = [names.index(c) for c in keep]
    chains = []
    for p in paths:
        rows = []
        for line in open(p):
            if line.startswith("#") or line.startswith("lp__"):
                continue
            parts = line.strip().split(",")
            if len(parts) < len(names):     # a partial final row
                continue
            rows.append([float(parts[j]) for j in idx])
        if rows:
            chains.append(np.array(rows))
    if len(chains) < 2:
        raise ValueError(f"only {len(chains)} chain(s) with draws under {paths[0]}")
    n = min(len(c) for c in chains)
    stacked = np.stack([c[:n] for c in chains])          # (M, N, len(keep))
    return {c: stacked[:, :, j] for j, c in enumerate(keep)}, n
