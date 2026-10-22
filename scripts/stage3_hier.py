"""Stage 3 of Section 5.2: the full model at the smallest orders on a few
tens of supernovae.  First contact with the hierarchy.

    python scripts/stage3_hier.py --check-only     gradient vs finite differences; cost
    python scripts/stage3_hier.py                  four chains, four jobs
    python scripts/stage3_hier.py --dispersed      four chains from jittered starts
    python scripts/stage3_hier.py --report [dir]   diagnostics of an existing run

Stage 4 is scripts/stage4_full.py, the same program on the training sample.
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
from cmdstanpy import from_csv, write_stan_json

from dgsn.data.sample import analysis_sample, median_x1, stage3_sample, x1_fifth
from dgsn.salt.globals import fit_globals
from dgsn.stan import data as D
from dgsn.stan.hier_data import build, inits, inits_from_optimize
from dgsn.stan.run import log_prob, model, sample_four, seconds_per_gradient

K_PSI, K_Q = 12, 8      # the smallest orders; override with --k-psi/--k-q
BOUNDED = {"tau_q", "tau_psi", "alpha_q_m", "alpha_psi_m", "sigma_g", "sigma_r", "s_delta"}
GLOBAL_SCALARS = ["tau_q", "tau_psi", "alpha_q", "alpha_psi", "sigma_g", "sigma_r",
                  "mu_p_rel", "mu_psi_rel", "mu_delta", "s_delta"]


def gradient_check(data, point, n_coords=20, h=1e-5, seed=0, name="hier"):
    """Central differences of the log density in a random subset of coordinates."""
    lp0, grad = log_prob(name, data, point)
    # flatten the parameter dict in Stan's order to index the gradient
    names = list(model(name).src_info()["parameters"].keys())
    # Stan's gradient is on the unconstrained scale, where a sum_to_zero_vector
    # of length N carries N - 1 coordinates though it is reported as N.  Index
    # by the constrained length and every block after it is off by one, which
    # showed up here as a finite-difference mismatch of 87.
    short = {"x_s_dev": 1}
    flat = []
    for n in names:
        v = np.atleast_1d(np.asarray(point[n], dtype=float))
        flat += [(n, i) for i in range(len(v) - short.get(n, 0))]
    rng = np.random.default_rng(seed)
    # x_s_dev is a sum_to_zero_vector: moving one element alone leaves the
    # constraint violated and Stan rejects the point.  x_s reaches the forward
    # map as s_bar + x_s_dev and s_bar is checked, so the path is covered.
    movable = [j for j, (nm, _) in enumerate(flat) if nm != "x_s_dev"]
    idx = rng.choice(movable, min(n_coords, len(movable)), replace=False)
    worst = 0.0
    for j in idx:
        n, i = flat[j]
        up, dn = {k: np.array(v, dtype=float, copy=True) for k, v in point.items()}, \
                 {k: np.array(v, dtype=float, copy=True) for k, v in point.items()}
        if np.ndim(up[n]) == 0:
            up[n] = up[n] + h; dn[n] = dn[n] - h
        else:
            up[n][i] += h; dn[n][i] -= h
        fd = (log_prob(name, data, up)[0] - log_prob(name, data, dn)[0]) / (2 * h)
        # log_prob's gradient is on the unconstrained scale: for a parameter
        # bounded below at zero that is theta d/dtheta.
        if n in BOUNDED:
            fd *= float(np.atleast_1d(point[n])[i])
        rel = abs(grad[j] - fd) / max(abs(fd), 1e-3)
        worst = max(worst, rel)
    return lp0, worst


def cost(data, point, name="hier"):
    """Seconds per log_prob call at the full epoch set and at one epoch per
    object: the difference is the epoch loops, the remainder the node sweeps
    plus process launch.  Launch is measured separately from a run with
    N_sn = 1 and one epoch, and removed."""
    def timeit(d, reps=5):
        t = time.time()
        for _ in range(reps):
            log_prob(name, d, point)
        return (time.time() - t) / reps
    full = timeit(data)
    first = data["start"]
    one = dict(data, N=data["N_sn"], dt=data["dt"][first - 1], f=data["f"][first - 1],
               sig_f=data["sig_f"][first - 1], band=data["band"][first - 1],
               start=np.arange(1, data["N_sn"] + 1), stop=np.arange(1, data["N_sn"] + 1))
    nodes = timeit(one)
    return full, nodes


def report(fit, names, sample, data):
    print(fit.diagnose())
    summ = fit.summary()
    print("\nglobal scalars:")
    # hier_nomode.stan has no mu_psi, c or x_psi, so take what is present
    have = [g for g in GLOBAL_SCALARS if g in summ.index]
    print(summ.loc[have, ["Mean", "StdDev", "5%", "95%", "ESS_bulk", "R_hat"]].to_string())
    for v in ("a_rel", "c_rel", "b", "d"):
        rows = summ.loc[[r for r in summ.index if r.startswith(v + "[")]]
        if len(rows):
            print(f"{v}: R_hat max {rows.R_hat.max():.3f}, ESS_bulk min {rows.ESS_bulk.min():.0f}")
    d = fit.draws_pd()
    chains = sorted(d.chain__.unique())
    branch = f"flip_q {[int(d[d.chain__ == c].flip_q.mean().round()) for c in chains]}"
    switch = int((d.groupby('chain__').flip_q.nunique() > 1).sum())
    if "flip_psi" in d.columns:
        branch += f", flip_psi {[int(d[d.chain__ == c].flip_psi.mean().round()) for c in chains]}"
        switch += int((d.groupby('chain__').flip_psi.nunique() > 1).sum())
    print(f"\nsign branches per chain: {branch}; draws switching within a chain: {switch}")
    print(f"frac_tapered mean {d.frac_tapered.mean():.2e} max {d.frac_tapered.max():.2e}; frac_clamped mean {d.frac_clamped.mean():.2e} max {d.frac_clamped.max():.2e}")
    print(f"seconds per gradient (upper bound): {seconds_per_gradient(fit):.2e}; treedepth median {d.treedepth__.median():.0f}, leapfrogs/draw median {d.n_leapfrog__.median():.0f}")
    # per-SN: x_psi width vs prior, timing-block correlations, interval check
    xs = fit.stan_variable("x_s"); xp = fit.stan_variable("x_p_rel"); dl = fit.stan_variable("delta_t")
    has_mode = "x_psi_rel[1]" in summ.index
    xpsi = fit.stan_variable("x_psi_rel") if has_mode else np.zeros_like(xp)
    corr = lambda u, v: np.array([np.corrcoef(u[:, i], v[:, i])[0, 1] for i in range(u.shape[1])])
    print(f"\nper-SN posterior sd of"
          + (f" x_psi: median {np.median(xpsi.std(0)):.2f} (prior 1);" if has_mode else "")
          + f" x_p {np.median(xp.std(0)):.2f}; delta_t {np.median(dl.std(0)):.2f} d")
    print(f"per-SN correlations, median over objects: x_s-x_p {np.median(corr(xs, xp)):+.2f}, "
          f"x_p-delta {np.median(corr(xp, dl)):+.2f}, x_s-delta {np.median(corr(xs, dl)):+.2f}"
          + (f", x_psi-delta {np.median(corr(xpsi, dl)):+.2f}" if has_mode else ""))
    smin, smax = fit.stan_variable("s_pmin").mean(0), fit.stan_variable("s_pmax").mean(0)
    print(f"s(-15 d) range over objects [{smin.min():.2f}, {smin.max():.2f}], s(+40 d) [{smax.min():.2f}, {smax.max():.2f}] against {D.ARC_RANGE}")
    lp_sn = fit.stan_variable("lp_sn").mean(0)
    n_ep = data["stop"] - data["start"] + 1
    per_ep = -2 * lp_sn / n_ep
    fifth = x1_fifth(sample.loc[names])
    print("\n-2 lp per epoch by x1 fifth (mean over objects):", " ".join(f"{q}: {per_ep[fifth.to_numpy() == q].mean():.2f}" for q in range(5)))
    worst = np.argsort(per_ep)[-3:][::-1]
    print("worst objects:", ", ".join(f"{names[i]} ({per_ep[i]:.1f}, x1 {sample.loc[names[i]].x1:+.2f})" for i in worst))
    print(f"delta_t: mean over objects {dl.mean():+.3f} d, sd of object means {dl.mean(0).std():.3f}; mu_delta {summ.loc['mu_delta', 'Mean']:+.3f}, s_delta {summ.loc['s_delta', 'Mean']:.3f}")


# The parameters block of hier.stan in declaration order, and which are
# log-transformed to the unconstrained space.
PARAMS = ["z_a", "b_til", "z_c", "d_til", "tau_q", "tau_psi", "alpha_q_m", "alpha_psi_m",
          "sigma_g", "sigma_r", "mu_p", "mu_psi", "mu_delta", "s_delta",
          "m_peak", "colr", "x_s", "u_p", "u_psi", "u_del"]


def warm_start(run_dir):
    """Last draw of each chain as its init, and the mean within-chain
    covariance of the unconstrained draws as an inverse metric.  Within,
    not pooled: chains on opposite sign branches differ in raw z_c."""
    import pandas as pd
    chains = [pd.read_csv(f, comment="#") for f in sorted(Path(run_dir).glob("hier-?.csv"))]
    cols = [c for p in PARAMS for c in chains[0].columns if c == p or c.startswith(p + ".")]
    assert len(cols) == len(set(cols))
    unc = []
    for d in chains:
        u = d[cols].to_numpy(dtype=float).copy()
        for j, c in enumerate(cols):
            if c.split(".")[0] in BOUNDED:
                u[:, j] = np.log(u[:, j])
        unc.append(u - u.mean(0))
    cov = sum(u.T @ u for u in unc) / (sum(len(u) for u in unc) - len(unc))
    inits = []
    for d in chains:
        last = d.iloc[-1]
        inits.append({p: (float(last[p]) if p in d.columns
                          else last[[c for c in cols if c.startswith(p + ".")]].to_numpy(dtype=float))
                      for p in PARAMS})
    print(f"warm start from {run_dir}: {len(cols)} unconstrained, {sum(len(u) for u in unc)} draws; "
          f"inverse-metric condition number {np.linalg.cond(cov):.1e}")
    return inits, cov


def dispersed_starts(ini, R, k_psi: int, k_q: int, mode: bool, n: int = 4, seed: int = 7):
    """``n`` starts jittered about ``ini``, for the dispersed-start check.

    Agreement between chains from dispersed starts is what checks the \\salt\\
    initialization (Section 5), so the dispersion has to reach the parameters
    that decide where a chain goes.
    """
    # The coefficients are jittered where their priors are written --- in
    # z_a, b_til, z_c, d_til --- and not in the rotated coordinates the
    # sampler uses, so the dispersion means the same thing whatever
    # rotation is in force; the result is turned back through R.  This
    # block addressed ini["z_a"] directly until 2026-10-22, which has not
    # existed since the coefficients were packed into w_coef, so every
    # dispersed run raised KeyError before reaching the sampler.
    from dgsn.stan.hier_data import n_coef, pack, unpack

    rng = np.random.default_rng(seed)
    Rm = np.eye(n_coef(k_psi, k_q, mode)) if R is None else np.asarray(R, dtype=float)
    raw0 = Rm @ np.asarray(ini["w_coef"], dtype=float)
    starts = []
    for _ in range(n):
        s = {kk: np.array(v, dtype=float, copy=True) if np.ndim(v) else float(v)
             for kk, v in ini.items()}
        co = unpack(raw0, k_psi, k_q, mode)
        # the non-centred pair is unit scale under its own prior
        for key in ("z_a", "z_c"):
            if key in co:
                co[key] = co[key] + 0.5 * rng.normal(size=len(co[key]))
        # the centred pair in coefficient units: half the prior width at the start
        for key, tau, alpha_m, K in (("b_til", ini["tau_q"], ini["alpha_q_m"], k_q),
                                     ("d_til", ini["tau_psi"], ini["alpha_psi_m"], k_psi)):
            kk = np.arange(1, K)
            w = tau * (kk / np.exp(np.log(kk).mean())) ** -(alpha_m + 0.5)
            co[key] = co[key] + 0.5 * w * rng.normal(size=K - 1)
        raw = (pack(co["z_a"], co["b_til"], co["z_c"], co["d_til"]) if mode
               else pack(co["z_a"], co["b_til"], co["d_til"]))
        s["w_coef"] = raw if R is None else np.linalg.solve(Rm, raw)
        # The common locations are dispersed too, not only the deviations
        # about them: mu_p and mu_delta are the slowest-mixing directions
        # measured, so four chains sharing them would leave the check
        # blind to exactly what it is meant to test.
        if "x_s_dev" in s:
            # stay on the sum-to-zero hyperplane, or the start is invalid
            j = 0.3 * rng.normal(size=len(s["x_s_dev"]))
            s["x_s_dev"] = s["x_s_dev"] + j - j.mean()
        for key, w in (("mu_p", 1.0), ("mu_psi", 1.0), ("mu_delta", 0.5),
                       ("s_bar", 0.3)):
            if key in s:
                s[key] = s[key] + w * float(rng.normal())
        for key, w in (("m_peak", 0.3), ("colr", 0.3), ("x_s", 0.3),
                       ("u_p", 1.0), ("u_psi", 1.0), ("u_del", 1.0)):
            if key in s:
                s[key] = s[key] + w * rng.normal(size=len(s[key]))
        starts.append(s)
    return starts


def main():
    global K_PSI, K_Q
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-only", action="store_true")
    ap.add_argument("--dispersed", action="store_true")
    ap.add_argument("--report", nargs="?", const="", default=None, help="read CSVs from this directory (default --out)")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--out", default="runs/stage3")
    ap.add_argument("--iter", type=int, default=None, help="warmup and sampling iterations each (default Stan's)")
    ap.add_argument("--max-treedepth", type=int, default=12)
    ap.add_argument("--k-psi", type=int, default=K_PSI)
    ap.add_argument("--k-q", type=int, default=K_Q)
    ap.add_argument("--dense", action="store_true", help="dense_e metric, adapted by Stan")
    ap.add_argument("--adapt-delta", type=float, default=None)
    ap.add_argument("--warm-from", default=None,
                    help="probe: start each chain at the last draw of a previous run's chain, with that "
                         "run's within-chain covariance as a fixed dense metric (step size alone adapted)")
    ap.add_argument("--precondition-optimize", action="store_true",
                    help="build the rotation at the optimize mode, which needs no previous run; "
                         "the only option after a change of basis")
    ap.add_argument("--init-optimize", action="store_true",
                    help="start the chains at Stan's own optimum on the data, re-placed on "
                         "the s_delta ridge, instead of at the SALT2 fit; needs no template "
                         "except to start the optimizer")
    ap.add_argument("--no-mode", action="store_true",
                    help="hier_nomode.stan: Eq. (curvlin) without x_psi Psi_1, so one "
                         "shared curve and no turning mode")
    ap.add_argument("--precondition-from", default=None, metavar="DIR",
                    help="rotate the coefficients by the Fisher information at the posterior "
                         "mean of the run in DIR; built at a start instead it makes things worse")
    args = ap.parse_args()

    K_PSI, K_Q = args.k_psi, args.k_q
    prog = "hier_nomode" if args.no_mode else "hier"
    mode = not args.no_mode
    sample = analysis_sample()
    names = stage3_sample()
    x1_ref = median_x1()
    G = fit_globals(K_PSI, K_Q, x1_ref)
    k_psi, k_q = args.k_psi, args.k_q
    data = build(names, k_psi, k_q, sample=sample, mode=mode)
    ini = inits(names, k_psi, k_q, x1_ref, sample=sample, globals_=G, mode=mode)
    R = None
    if args.precondition_optimize:
        from dgsn.stan.precondition import from_optimize
        import time as _t
        t0 = _t.time()
        R = from_optimize(data, ini, k_psi, k_q)
        data = build(names, k_psi, k_q, sample=sample, R=R)
        print(f"preconditioned at the optimize mode in {_t.time()-t0:.0f} s, "
              f"rotation condition {np.linalg.cond(R):.1e}")
    elif args.precondition_from:
        from dgsn.stan.precondition import from_run
        R = from_run(args.precondition_from, data, k_psi, k_q)
        data = build(names, k_psi, k_q, sample=sample, R=R)
        print(f"preconditioned from {args.precondition_from}: "
              f"rotation condition {np.linalg.cond(R):.1e}")
    ini = inits(names, k_psi, k_q, x1_ref, sample=sample, globals_=G, R=R, mode=mode)
    if args.init_optimize:
        # the SALT2 fit starts the optimizer; the optimum starts the chains.
        # polish=False: its restarts cannot terminate in s_delta, flat in
        # likelihood and rising in prior, so they walk the ridge for four
        # passes and the re-placement discards the distance -- 182 s against
        # fifteen minutes, at the same point.
        ini = inits_from_optimize(prog, data, ini, k_psi, k_q, mode=mode, polish=False)
    print(f"{prog}: {len(names)} objects, {data['N']} epochs, K = ({K_PSI}, {K_Q}); x1 range [{sample.loc[names].x1.min():+.2f}, {sample.loc[names].x1.max():+.2f}]")
    print(f"start: tau_q {ini['tau_q']:.2f} tau_psi {ini['tau_psi']:.2f}, "
          f"|w_coef| {np.linalg.norm(ini['w_coef']):.2f}, "
          f"mu_p {ini['mu_p']:.2f}"
          + (f" mu_psi {ini['mu_psi']:.2f}" if 'mu_psi' in ini else ''))

    if args.report is not None:
        fit = from_csv([str(c) for c in sorted(Path(args.report or args.out).glob("hier-?.csv"))])
        report(fit, names, sample, data)
        return

    lp0, worst = gradient_check(data, ini, name=prog)
    print(f"lp at start {lp0:.2f}; gradient vs finite differences, worst of 20 coordinates: {worst:.1e}")
    full, nodes = cost(data, ini, name=prog)
    per_ep = (full - nodes) / (data["N"] - data["N_sn"])
    print(f"cost per log_prob call: {full:.3f} s at {data['N']} epochs, {nodes:.3f} s at one epoch per object (nodes + launch) "
          f"-> {per_ep * 1e6:.0f} us per epoch; projected per gradient for 442 objects at ~31400 epochs: "
          f"{per_ep * 31400 + (nodes - 0.0) * 442 / data['N_sn']:.2f} s (upper bound: node term includes launch)")
    if args.check_only:
        return

    starts = ini
    if args.dispersed:
        starts = dispersed_starts(ini, R, K_PSI, K_Q, mode)
    out = Path(args.out + ("-dispersed" if args.dispersed else ""))
    metric, metric_file, warmup, refresh = ("dense_e" if args.dense else None), None, args.iter, 100
    if args.warm_from:
        starts, inv_metric = warm_start(args.warm_from)
        out.mkdir(parents=True, exist_ok=True)
        metric_file = out / "inv_metric.json"
        write_stan_json(str(metric_file), {"inv_metric": inv_metric})
        # below 20 warmup iterations Stan adapts the step size and not the metric
        metric, warmup, refresh = "dense_e", 19, 10
    t = time.time()
    fit = sample_four(prog, data, starts, out, seed=1, threads_per_chain=args.threads,
                      iter_warmup=warmup, iter_sampling=args.iter, max_treedepth=args.max_treedepth,
                      metric=metric, metric_file=metric_file, adapt_delta=args.adapt_delta,
                      refresh=refresh)
    print(f"\nfour chains in {time.time() - t:.0f} s wall")
    report(fit, names, sample, data)


if __name__ == "__main__":
    main()
