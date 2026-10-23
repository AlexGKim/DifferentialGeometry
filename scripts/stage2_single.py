"""Stage 2 of Section 5.2: one supernova, six parameters, the global
functions supplied as data from the SALT2 fits.

    python scripts/stage2_single.py [--check-only] [--name ZTF..]

--check-only runs the build check (single_sn's fluxes against forward.stan
for the same inputs) and the gradient against finite differences at the
start, and stops.  Otherwise four chains are run as four jobs and the
posterior reported against the archive.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

from dgsn.data.epochs import likelihood_epochs
from dgsn.data.sample import analysis_sample, median_x1
from dgsn.salt.globals import fit_globals
from dgsn.salt.tracks import per_sn_start, salt_track
from dgsn.stan import data as D
from dgsn.stan.run import fixed_param, log_prob, sample_four

NAMES = ["m_peak", "colr", "x_s", "x_p", "x_psi", "delta_t"]


def pick_object(sample, min_g_prepeak=3):
    """Near the median epoch count and the median x1, among objects with at
    least ``min_g_prepeak`` good g epochs before t0, so that m_peak is tested
    where it is measured rather than extrapolated from the decline."""
    n = sample.ztfg + sample.ztfr
    score = ((n - n.median()) / n.std()) ** 2 + ((sample.x1 - sample.x1.median()) / sample.x1.std()) ** 2
    for name in score.sort_values().index:
        row = sample.loc[name]
        ep = likelihood_epochs(name, row.t0, row.redshift)
        if ep["n_g_prepeak"] >= min_g_prepeak:
            return row
    raise RuntimeError("no object with enough pre-peak g epochs")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-only", action="store_true")
    ap.add_argument("--name", default=None)
    ap.add_argument("--out", default="runs/stage2")
    args = ap.parse_args()

    sample = analysis_sample()
    row = sample.loc[args.name] if args.name else pick_object(sample)
    name = row.name
    ep = likelihood_epochs(name, row.t0, row.redshift)
    print(f"{name}: z = {row.redshift:.4f}, x1 = {row.x1:+.2f}, c = {row.c:+.3f}, epochs g {ep['n_g']} r {ep['n_r']} "
          f"(g before t0: {ep['n_g_prepeak']}), m_bright {ep['m_bright']:.2f}")

    x1_ref = median_x1()
    G = fit_globals(D.K_PSI, D.K_Q, x1_ref)
    print(f"globals from SALT2 at x1 = {x1_ref:.3f}; reference at x_s {G['x_s']:.3f}, "
          f"x_p {G['x_p']:.3f}, x_psi {G['x_psi']:.3f}; likelihood on the archive's errors alone")

    start = per_sn_start(row)
    start.update(x_s=G["x_s"], x_p=G["x_p"], x_psi=G["x_psi"])
    print("start:", {k: round(v, 3) for k, v in start.items()})

    data = D.grid_data()
    data.update({k: G[k] for k in ("a", "b", "c", "d")})
    data.update(D.epochs_data(ep["dt"], ep["f"], ep["sig_f"], ep["band"], ep["z"]))
    data.update(D.prior_data(prior_on=True, mu_p=G["x_p"], mu_psi=G["x_psi"], mu_delta=0.0,
                             s_delta=1.0, xs_bar=G["x_s"], m_bright=ep["m_bright"]))

    # Build check: the same forward map through the other program.
    fdata = D.grid_data()
    fdata.update({k: G[k] for k in ("a", "b", "c", "d")})
    fdata.update(start)
    fdata.update(z=ep["z"])
    fdata.update(D.eval_points(dt=ep["dt"], band=ep["band"]))
    f_forward = fixed_param("forward", fdata)["f_hat"]
    fit0 = sample_one_fixed(data, start)
    print(f"build check: max |f_hat(single_sn) / f_hat(forward) - 1| = {np.abs(fit0 / f_forward - 1).max():.1e}")

    lp, grad = log_prob("single_sn", data, start)
    fd = np.empty(6)
    for i, k in enumerate(NAMES):
        up, dn = dict(start), dict(start)
        up[k] += 1e-5
        dn[k] -= 1e-5
        fd[i] = (log_prob("single_sn", data, up)[0] - log_prob("single_sn", data, dn)[0]) / 2e-5
    print(f"gradient check at the start: lp = {lp:.3f}, max relative disagreement with finite differences "
          f"{(np.abs(grad - fd) / np.maximum(np.abs(fd), 1e-3)).max():.1e}")
    if args.check_only:
        return

    out = Path(args.out) / name
    fit = sample_four("single_sn", data, start, out, seed=1)
    print(fit.diagnose())
    summ = fit.summary()
    print(summ.loc[NAMES + ["s_pmin", "s_pmax", "kappa0", "frac_tapered", "frac_clamped"]].to_string())
    print(f"\ndelta_t = t_max - t0: {summ.loc['delta_t', 'Mean']:+.3f} +- {summ.loc['delta_t', 'StdDev']:.3f} d")
    print(f"m_peak {summ.loc['m_peak', 'Mean']:.3f} vs SALT2-implied {start['m_peak']:.3f}; "
          f"colr {summ.loc['colr', 'Mean']:.3f} vs {start['colr']:.3f}")
    draws = fit.draws_pd()[NAMES].to_numpy()
    corr = np.corrcoef(draws, rowvar=False)
    print("\nposterior correlations:")
    print("          " + "".join(f"{k:>8}" for k in NAMES))
    for i, k in enumerate(NAMES):
        print(f"  {k:>8}" + "".join(f"{corr[i, j]:8.2f}" for j in range(len(NAMES))))
    f_hat = fit.stan_variable("f_hat").mean(0)
    r = (ep["f"] - f_hat) / ep["sig_f"]
    phase = (ep["dt"] - summ.loc["delta_t", "Mean"]) / (1 + ep["z"])
    print("\nresiduals / sigma_f by phase bin (g | r):")
    for lo, hi in ((-15, -5), (-5, 5), (5, 15), (15, 25), (25, 40)):
        inb = (phase >= lo) & (phase < hi)
        g, rr = r[inb & (ep["band"] == 1)], r[inb & (ep["band"] == 2)]
        print(f"  [{lo:+3d},{hi:+3d}) d: g n={len(g):2d} mean {g.mean() if len(g) else np.nan:+.2f} rms {np.sqrt(np.mean(g**2)) if len(g) else np.nan:.2f} | "
              f"r n={len(rr):2d} mean {rr.mean() if len(rr) else np.nan:+.2f} rms {np.sqrt(np.mean(rr**2)) if len(rr) else np.nan:.2f}")


def sample_one_fixed(data, start):
    """single_sn's generated quantities at the start: fixed_param on the sampler program."""
    from dgsn.stan.run import model
    fit = model("single_sn").sample(data=data, inits=start, fixed_param=True, chains=1, iter_warmup=0,
                                    iter_sampling=1, adapt_engaged=False, sig_figs=18, show_progress=False)
    return np.asarray(fit.stan_variable("f_hat"))[0]


if __name__ == "__main__":
    main()
