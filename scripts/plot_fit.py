"""Plots for examining a hierarchical fit: the global parameters, and what the
model says a supernova does.

    python scripts/plot_fit.py --out runs/stage4-smoke
    python scripts/plot_fit.py --out runs/stage3-dense --n-draws 60

Two figures land in the run directory.

``corner.png`` is the ten global scalars, with the four chains drawn in
different colours rather than pooled. Pooling hides what matters in a run that
has not converged: chains sitting in different places look like a wide
posterior. Drawn this way, disagreement is visible as separated clouds.

``lightcurves.png`` is, for each of a few objects, the track through magnitude
space on the left and the light curve against the data on the right, one line
per posterior draw. The track is the object the model is actually about --- a
curve in the plane of the two magnitudes, which the supernova walks along ---
and the turn near the anchor, where both bands peak within about a day, is
what the turning basis exists to resolve.

The model is never evaluated in Python: every curve here comes from
``forward.stan`` at one draw's parameters, through ``fixed_param``.  The grid
constants come from the run's own ``data.json``, not from today's defaults, so
that a run made before a constant changed still plots correctly.  Data are
shown in native flux, as the archive gives them and as the likelihood sees
them; nothing is converted to magnitudes.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from dgsn.data.epochs import likelihood_epochs  # noqa: E402
from dgsn.data.sample import analysis_sample, training_split  # noqa: E402
from dgsn.stan.run import chain_csvs, fixed_param  # noqa: E402

GLOBALS = ["tau_q", "tau_psi", "alpha_q", "alpha_psi",
           "mu_p_rel", "mu_psi_rel", "c_bar", "b_0"]
LOGGED = {"tau_q", "tau_psi"}                 # plotted as log10, their range spanning decades
CHAIN_COLOURS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
GRID_KEYS = ["Kpsi", "Kq", "m", "J", "delta", "delta_p", "pmin", "pmax",
             "smin", "smax", "gsmin", "gsmax", "eps", "eps_t", "eps_s", "zp"]


def read_draws(out: Path) -> pd.DataFrame:
    """Post-warmup draws of every chain, with a chain label."""
    frames = []
    for k in chain_csvs(out):
        hdr, rows = None, []
        for line in open(k):
            if line.startswith("#"):
                continue
            if hdr is None:
                hdr = line.strip().split(",")
                continue
            p = line.strip().split(",")
            if len(p) == len(hdr):
                rows.append(p)
        f = pd.DataFrame(np.asarray(rows[:-1], dtype=float), columns=hdr)
        f["chain"] = int(k.stem.split("-")[1])
        frames.append(f)
    return pd.concat(frames, ignore_index=True)


def corner(d: pd.DataFrame, path: Path) -> None:
    """Lower-triangle scatter with per-chain marginals on the diagonal."""
    cols = [c for c in GLOBALS if c in d.columns]
    n = len(cols)
    lab = {c: (f"log10 {c}" if c in LOGGED else c) for c in cols}
    val = {c: (np.log10(d[c]) if c in LOGGED else d[c]) for c in cols}
    fig, ax = plt.subplots(n, n, figsize=(2.0 * n, 2.0 * n))
    for i in range(n):
        for j in range(n):
            a = ax[i, j]
            if j > i:
                a.axis("off")
                continue
            for ch, colour in zip(sorted(d.chain.unique()), CHAIN_COLOURS):
                m = d.chain == ch
                if i == j:
                    a.hist(val[cols[i]][m], bins=30, histtype="step", color=colour,
                           density=True, lw=1.1)
                else:
                    a.plot(val[cols[j]][m], val[cols[i]][m], ".", ms=1.2, alpha=0.35,
                           color=colour, rasterized=True)
            if i == n - 1:
                a.set_xlabel(lab[cols[j]], fontsize=8)
            else:
                a.set_xticklabels([])
            if j == 0 and i != 0:
                a.set_ylabel(lab[cols[i]], fontsize=8)
            else:
                a.set_yticklabels([])
            a.tick_params(labelsize=6)
    handles = [plt.Line2D([], [], color=c, lw=2) for c in CHAIN_COLOURS[:d.chain.nunique()]]
    fig.legend(handles, [f"chain {c}" for c in sorted(d.chain.unique())],
               loc="upper right", fontsize=11, frameon=False)
    fig.suptitle("Global parameters, by chain", fontsize=13, x=0.35, y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"wrote {path}")


def evaluate(grid: dict, draw: pd.Series, idx: int, z: float, k_psi: int, k_q: int,
             p_eval: np.ndarray, s_eval: np.ndarray,
             dt_obs: np.ndarray | None = None,
             band_obs: np.ndarray | None = None) -> dict:
    """One draw's curve and light curve, from forward.stan.

    ``dt_obs``/``band_obs`` add the observed epochs to the same call, so that
    ``f_obs`` is the model at the data in one forward evaluation rather than
    an interpolation off the plotting grid.

    ``idx`` is the supernova's 1-based position in the fit.  The sign-relabeled
    amplitudes are used where the fit provides them, so that chains in opposite
    branches plot on top of one another rather than mirrored."""
    def vec(name, n):
        return np.array([draw[f"{name}.{i}"] for i in range(1, n)])

    # The mode-free ladder (--no-trav-mode) has no a, c, x_p or x_psi at all:
    # forward.stan still wants them, and zeros reproduce what hier_common
    # evaluates, which passes a_zero, c_zero and 0, 0 into log_flux.
    has_mode = "a.1" in draw or "a_rel.1" in draw
    a = vec("a_rel" if "a_rel.1" in draw else "a", k_q) if has_mode else np.zeros(k_q - 1)
    c = vec("c_rel" if "c_rel.1" in draw else "c", k_psi) if has_mode else np.zeros(k_psi - 1)
    x_p = draw.get(f"x_p_rel.{idx}", draw.get(f"x_p.{idx}", 0.0)) if has_mode else 0.0
    x_psi = draw.get(f"x_psi_rel.{idx}", draw.get(f"x_psi.{idx}", 0.0)) if has_mode else 0.0
    # magnitudes at rest-frame phase p are asked for as dt = p (1 + z) + delta_t,
    # the model carrying dt back to p by Eq. (modeltime)
    dt = np.concatenate([p_eval * (1 + z) + draw[f"delta_t.{idx}"]] * 2)
    band = np.concatenate([np.ones(len(p_eval), int), 2 * np.ones(len(p_eval), int)])
    if dt_obs is not None:
        dt = np.concatenate([dt, np.asarray(dt_obs, float)])
        band = np.concatenate([band, np.asarray(band_obs, int)])
    data = dict(grid)
    data.update(a=a, b=vec("b", k_q), c=c, d=vec("d", k_psi),
                m_peak=float(draw[f"m_peak.{idx}"]), colr=float(draw[f"colr.{idx}"]),
                x_s=float(draw.get(f"x_s_tot.{idx}", draw.get(f"x_s.{idx}"))),
                x_p=float(x_p), x_psi=float(x_psi),
                delta_t=float(draw[f"delta_t.{idx}"]), z=float(z),
                Np=len(p_eval), p_eval=p_eval, Ns=len(s_eval), s_eval=s_eval,
                N=len(dt), band=band, dt=dt)
    out = fixed_param("forward", data)
    f = out["f_hat"]
    return dict(gamma=out["gamma"], s_of_p=out["s_of_p"], psi=out["psi"], f_g=f[:len(p_eval)],
                f_r=f[len(p_eval):2 * len(p_eval)],
                f_obs=f[2 * len(p_eval):] if dt_obs is not None else None,
                m_peak=float(draw[f"m_peak.{idx}"]),
                colr=float(draw[f"colr.{idx}"]))


def lightcurves(d: pd.DataFrame, names: list[str], index: dict,
                sample, grid: dict, n_draws: int, path: Path) -> None:
    k_psi, k_q = int(grid["Kpsi"]), int(grid["Kq"])
    p_eval = np.linspace(grid["pmin"] + 0.5, grid["pmax"] - 0.5, 160)
    s_eval = np.linspace(grid["gsmin"], grid["gsmax"], 220)
    rng = np.random.default_rng(0)
    pick = rng.choice(len(d), size=min(n_draws, len(d)), replace=False)

    fig, ax = plt.subplots(len(names), 2, figsize=(13, 3.6 * len(names)))
    ax = np.atleast_2d(ax)
    for row, name in enumerate(names):
        r = sample.loc[name]
        ep = likelihood_epochs(name, r.t0, r.redshift)
        idx = index[name]
        curves = [evaluate(grid, d.iloc[i], idx, ep["z"], k_psi, k_q, p_eval, s_eval)
                  for i in pick]

        at = ax[row, 0]
        for cv in curves:
            at.plot(cv["gamma"][:, 0] + cv["m_peak"],
                    cv["gamma"][:, 1] + cv["m_peak"] - cv["colr"],
                    "-", lw=0.6, alpha=0.25, color="#1f77b4")
        at.plot([], [], "-", color="#1f77b4", label=f"{len(curves)} draws")
        at.invert_xaxis(); at.invert_yaxis()
        at.set_xlabel("$m_g$"); at.set_ylabel("$m_r$")
        at.set_title(f"{name}: the track through magnitude space", fontsize=10)
        at.legend(fontsize=8, frameon=False)

        al = ax[row, 1]
        for cv in curves:
            al.plot(p_eval, cv["f_g"], "-", lw=0.6, alpha=0.25, color="#1f77b4")
            al.plot(p_eval, cv["f_r"], "-", lw=0.6, alpha=0.25, color="#d62728")
        phase = (ep["dt"] - np.median(d[f"delta_t.{idx}"])) / (1 + ep["z"])
        for bnd, colour, mk in ((1, "#1f77b4", "o"), (2, "#d62728", "s")):
            m = ep["band"] == bnd
            al.errorbar(phase[m], ep["f"][m], yerr=ep["sig_f"][m], fmt=mk, ms=3.5,
                        color=colour, lw=0.8, capsize=0, zorder=3,
                        label=f"{'g' if bnd == 1 else 'r'} ({int(m.sum())} epochs)")
        al.axhline(0, color="0.7", lw=0.6)
        al.set_xlabel("rest-frame phase (d)"); al.set_ylabel("flux (zp 30)")
        al.set_title(f"$x_1$ {r.x1:+.2f}, $z$ {r.redshift:.4f}, "
                     f"$-2\\,$lp/epoch {-2 * np.median(d[f'lp_sn.{idx}']) / len(ep['dt']):.1f}",
                     fontsize=10)
        al.legend(fontsize=8, frameon=False)
    fig.suptitle("Posterior draws against the data, in the flux the likelihood sees", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/stage4-smoke")
    ap.add_argument("--n-draws", type=int, default=40)
    ap.add_argument("--names", nargs="*", default=None,
                    help="objects to plot (default: best, median and worst by -2 lp per epoch)")
    args = ap.parse_args()

    out = Path(args.out)
    grid_all = json.load(open(out / "data.json"))
    grid = {k: grid_all[k] for k in GRID_KEYS}
    grid["xg"] = grid_all.get("xg", list(np.polynomial.legendre.leggauss(int(grid["J"]))[0]))
    grid["wg"] = grid_all.get("wg", list(np.polynomial.legendre.leggauss(int(grid["J"]))[1]))
    d = read_draws(out)
    print(f"{len(d)} draws over {d.chain.nunique()} chains from {out}")

    corner(d, out / "corner.png")

    # the fit's object list, in the order the data were built
    sample = analysis_sample()
    train, _ = training_split()
    n_sn = int(grid_all["N_sn"])
    if n_sn == len(train):
        names_fit = train
    else:
        import sys
        sys.path.insert(0, str(Path(__file__).parent))
        from stage4_full import stratified
        names_fit = stratified(train, n_sn, sample)
    index = {nm: i + 1 for i, nm in enumerate(names_fit)}

    if args.names:
        names = args.names
    else:
        n_ep = np.array(grid_all["stop"]) - np.array(grid_all["start"]) + 1
        per_ep = np.array([-2 * np.median(d[f"lp_sn.{i}"]) / n_ep[i - 1]
                           for i in range(1, n_sn + 1)])
        order = np.argsort(per_ep)
        names = [names_fit[order[0]], names_fit[order[len(order) // 2]], names_fit[order[-1]]]
        print("objects: best, median and worst by -2 lp per epoch -> " + ", ".join(names))

    lightcurves(d, names, index, sample, grid, args.n_draws, out / "lightcurves.png")


if __name__ == "__main__":
    main()
