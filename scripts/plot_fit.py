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
from dgsn.data.sample import (analysis_sample, stage3_sample,  # noqa: E402
                              training_split)
from dgsn.stan.run import chain_csvs, fixed_param  # noqa: E402

# Every global scalar any of these programs samples; the list is filtered
# against the run's own columns, so a program that has none of a given one
# simply drops that panel.  mu_delta, kappa0 and the two sigma_b are
# hier_base's and hier_full's and were missing here, which left a hier_base
# corner plot showing two panels out of six.
GLOBALS = ["tau_q", "tau_psi", "alpha_q", "alpha_psi",
           "mu_p_rel", "mu_psi_rel", "mu_x_p", "mu_x_psi",
           "c_bar", "b_0", "mu_delta", "kappa0", "sigma_b.1", "sigma_b.2"]
LOGGED = {"tau_q", "tau_psi"}                 # plotted as log10, their range spanning decades
CHAIN_COLOURS = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"]
GRID_KEYS = ["Kpsi", "Kq", "m", "J", "delta", "delta_p", "pmin", "pmax",
             "smin", "smax", "gsmin", "gsmax", "eps", "eps_t", "eps_s", "zp"]


def read_draws(out: Path, keep: set[int] | None = None) -> pd.DataFrame:
    """Post-warmup draws of every chain, with a chain label.

    ``keep``, when given, is the set of chain numbers to read; the rest are
    skipped.  Excluding a chain is a claim about the run, so it is a command
    line argument and never a default: see ``--chains``.

    A run made with ``save_warmup=1`` writes its warmup iterations into the
    same CSV, so keeping every row silently pools adaptation with the
    posterior.  That is not a cosmetic difference: on ladder-base-k812-n60
    the warmup transient puts 80 percent of iterations 1-100 at b_0 < -10.5
    and NONE of the sampling draws there, so the pooled corner plot shows a
    second cluster and a third, and the fit has neither.  Both counts are
    read from the CSV's own header rather than assumed, so a run written
    without save_warmup keeps every row as before."""
    frames = []
    for k in chain_csvs(out):
        if keep is not None and int(k.stem.split("-")[1]) not in keep:
            continue
        hdr, rows, warm, saved = None, [], 0, False
        for line in open(k):
            if line.startswith("#"):
                t = line.lstrip("# ").replace(" ", "")
                if t.startswith("num_warmup="):
                    warm = int(t.split("=")[1].split("(")[0])
                elif t.startswith("save_warmup="):
                    saved = t.split("=")[1].split("(")[0] not in ("0", "false")
                continue
            if hdr is None:
                hdr = line.strip().split(",")
                continue
            p = line.strip().split(",")
            if len(p) == len(hdr):
                rows.append(p)
        rows = rows[:-1]
        if saved and warm:
            if len(rows) <= warm:
                raise SystemExit(
                    f"{k.name} has {len(rows)} rows and {warm} of them are warmup; "
                    f"the run has not left adaptation and there is nothing to plot")
            rows = rows[warm:]
        f = pd.DataFrame(np.asarray(rows, dtype=float), columns=hdr)
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
            for ch in sorted(d.chain.unique()):
                colour = CHAIN_COLOURS[(int(ch) - 1) % len(CHAIN_COLOURS)]
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
    used = sorted(int(c) for c in d.chain.unique())
    handles = [plt.Line2D([], [], color=CHAIN_COLOURS[(c - 1) % len(CHAIN_COLOURS)], lw=2)
               for c in used]
    fig.legend(handles, [f"chain {c}" for c in used],
               loc="upper right", fontsize=11, frameon=False)
    fig.suptitle(f"Global parameters, by chain  (chains {', '.join(map(str, used))})",
                 fontsize=13, x=0.35, y=0.98)
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
                sample, grid: dict, n_draws: int, path: Path,
                span: tuple[np.ndarray, np.ndarray] | None = None) -> None:
    k_psi, k_q = int(grid["Kpsi"]), int(grid["Kq"])
    p_eval = np.linspace(grid["pmin"] + 0.5, grid["pmax"] - 0.5, 160)
    rng = np.random.default_rng(0)
    pick = rng.choice(len(d), size=min(n_draws, len(d)), replace=False)

    fig, ax = plt.subplots(len(names), 3, figsize=(18.5, 3.6 * len(names)))
    ax = np.atleast_2d(ax)
    for row, name in enumerate(names):
        r = sample.loc[name]
        ep = likelihood_epochs(name, r.t0, r.redshift)
        idx = index[name]
        # Nothing is drawn outside the phase the object was observed over.
        # The model is defined well beyond it -- the phase clamp runs to
        # [pmin, pmax] and the arclength clamp to [smin, smax] -- but what it
        # says there is not constrained by this object's data, and drawing it
        # invents features: ZTF19aadoocn carried a flux peak of 1.5e5 at
        # p = -12.7, three days before its first epoch, larger than anything
        # the data show.  Those bounds are not the binding ones; the data are.
        # The tracks take each draw's own epoch reach rather than a pooled
        # median.  s(p) is monotonic, so each mask is one
        # contiguous stretch and a line through it joins nothing across a gap.
        phase = (ep["dt"] - np.median(d[f"delta_t.{idx}"])) / (1 + ep["z"])
        p_lo, p_hi = float(phase.min()), float(phase.max())
        p_ok = (p_eval >= p_lo) & (p_eval <= p_hi)
        # s_pmin and s_pmax are NOT the epochs' reach, whatever this file used
        # to say: hier_full.stan:518 sets them to arclength([pmin + 5,
        # pmax - 5]), the modelling window inset by five days, which for
        # ZTF19aadoocn is s in [-4.1, 3.5] against an epoch reach of about a
        # third of that.  The arclengths the epochs actually sit at are
        # s_epoch, sliced by this object's start and stop.
        if span is not None:
            lo_e, hi_e = int(span[0][idx - 1]), int(span[1][idx - 1])
            cols = [f"s_epoch.{j}" for j in range(lo_e, hi_e + 1)]
            se = d.iloc[pick][cols].to_numpy()
            s_bounds = list(zip(se.min(1), se.max(1)))
        else:
            s_bounds = [(float(d.iloc[i][f"s_pmin.{idx}"]),
                         float(d.iloc[i][f"s_pmax.{idx}"])) for i in pick]
        # The arclength grid is this object's, not the gauge span's.  A
        # well-sampled object reaches past [gsmin, gsmax] -- ZTF20aaodkvl
        # reaches s = -7.6 against a gauge span starting at -6 -- and a grid
        # fixed to the gauge span would cut its track at -6 and show a limit
        # of the plot as though it were a limit of the data.  Clipped to
        # [smin, smax], beyond which psi is not defined at all.
        s_eval = np.linspace(max(float(grid["smin"]), min(a for a, _ in s_bounds)),
                             min(float(grid["smax"]), max(b for _, b in s_bounds)), 220)
        curves = [evaluate(grid, d.iloc[i], idx, ep["z"], k_psi, k_q, p_eval, s_eval)
                  for i in pick]
        for (a_s, b_s), cv in zip(s_bounds, curves):
            cv["p_ok"] = p_ok
            cv["s_ok"] = (s_eval >= a_s) & (s_eval <= b_s)
        s_lo = float(np.median([a for a, _ in s_bounds]))
        s_hi = float(np.median([b for _, b in s_bounds]))

        # Eq. (modelmag): m_g = gamma_0 + m_peak, m_r = gamma_1 + m_peak - colr,
        # so the colour is gamma_0 - gamma_1 + colr and carries m_peak not at
        # all -- the colour-magnitude track is the curve's own shape plus the
        # two placement parameters, and nothing else.
        ac = ax[row, 0]
        for cv in curves:
            gr = cv["gamma"][:, 0] - cv["gamma"][:, 1] + cv["colr"]
            mr = cv["gamma"][:, 1] + cv["m_peak"] - cv["colr"]
            ok = cv["s_ok"]
            ac.plot(gr[ok], mr[ok], "-", lw=0.7, alpha=0.3, color="#6a3d9a")
        # the anchor, s = 0: gamma = (0,0), so the colour is colr exactly
        ac.plot([cv["colr"] for cv in curves],
                [cv["m_peak"] - cv["colr"] for cv in curves],
                ".", ms=3, color="#e31a1c", zorder=4, label="anchor $s=0$")
        ac.plot([], [], "-", color="#6a3d9a", lw=1.4,
                label=f"$s\\in[{s_lo:.1f},{s_hi:.1f}]$, the data's reach")

        ac.invert_yaxis()
        ac.set_xlabel("$g-r$ (mag)"); ac.set_ylabel("$r$ (mag)")
        ac.set_title(f"{name}: colour-magnitude track", fontsize=10)
        ac.legend(fontsize=7, frameon=False)

        at = ax[row, 1]
        for cv in curves:
            ok = cv["s_ok"]
            at.plot((cv["gamma"][:, 0] + cv["m_peak"])[ok],
                    (cv["gamma"][:, 1] + cv["m_peak"] - cv["colr"])[ok],
                    "-", lw=0.6, alpha=0.25, color="#1f77b4")
        at.plot([], [], "-", color="#1f77b4", label=f"{len(curves)} draws")
        at.invert_xaxis(); at.invert_yaxis()
        at.set_xlabel("$m_g$"); at.set_ylabel("$m_r$")
        at.set_title(f"{name}: the track through magnitude space", fontsize=10)
        at.legend(fontsize=8, frameon=False)

        al = ax[row, 2]
        for cv in curves:
            ok = cv["p_ok"]
            al.plot(p_eval[ok], cv["f_g"][ok], "-", lw=0.6, alpha=0.25, color="#1f77b4")
            al.plot(p_eval[ok], cv["f_r"][ok], "-", lw=0.6, alpha=0.25, color="#d62728")
        al.set_xlim(p_lo - 1.5, p_hi + 1.5)
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
    fig.suptitle("Posterior draws: the track in magnitude space, and the data in "
                 "the flux the likelihood sees", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"wrote {path}")


def resolve_names(n_sn: int, n_ep_run: np.ndarray, train: list[str],
                  sample) -> list[str]:
    """The fit's object list, confirmed against ``data.json``'s epoch counts.

    Several selections are the same length -- ``stage3_sample()`` and a
    thirty-object ``stratified()`` share one object of thirty -- so the list
    is chosen by reproducing ``stop - start + 1`` exactly, not by which
    script is presumed to have made the run.  A run whose selection matches
    nothing on offer raises instead of plotting the wrong photometry."""
    import sys
    sys.path.insert(0, str(Path(__file__).parent))
    from stage4_full import stratified

    cands = [("training_split", train)]
    if n_sn != len(train):
        cands.append(("stage3_sample", list(stage3_sample())))
        # stage3_hier.py builds its sample as stage3_sample(n_per_fifth=...),
        # so a run of 5k objects is the k-per-fifth selection and NOT the
        # default one.  Recovered the same way stage3_hier.py's own
        # _source_names recovers it, n_per_fifth = n_sn // 5.  Still a
        # candidate and not a guess: the epoch-count test below is what
        # accepts it, exactly as for every other entry here.
        if n_sn % 5 == 0 and n_sn // 5 != 6:
            cands.append((f"stage3_sample(n_per_fifth={n_sn // 5})",
                          list(stage3_sample(n_per_fifth=n_sn // 5))))
        cands.append(("stratified", list(stratified(train, n_sn, sample))))
    tried = []
    for label, nm in cands:
        if len(nm) != n_sn:
            tried.append(f"{label} (length {len(nm)})")
            continue
        cnt = np.array([len(likelihood_epochs(x, sample.loc[x].t0,
                                              sample.loc[x].redshift)["dt"])
                        for x in nm])
        if np.array_equal(cnt, n_ep_run):
            print(f"object list: {label} ({n_sn} objects, "
                  f"{int(n_ep_run.sum())} epochs, counts match data.json)")
            return nm
        tried.append(f"{label} (epochs {int(cnt.sum())} vs {int(n_ep_run.sum())})")
    raise SystemExit("cannot identify the fit's object list from data.json; tried "
                     + "; ".join(tried))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/stage4-smoke")
    ap.add_argument("--n-draws", type=int, default=40)
    ap.add_argument("--names", nargs="*", default=None,
                    help="objects to plot (default: best, median and worst by -2 lp per epoch)")
    ap.add_argument("--chains", type=int, nargs="*", default=None,
                    help="chain numbers to plot (default: all).  Dropping a chain asserts "
                         "it is not sampling the same distribution; the figure titles record "
                         "which chains were used so the claim travels with the plot.")
    ap.add_argument("--suffix", default="",
                    help="appended to the output filenames, to keep a restricted-chain "
                         "figure beside the full one rather than overwriting it")
    args = ap.parse_args()

    out = Path(args.out)
    grid_all = json.load(open(out / "data.json"))
    grid = {k: grid_all[k] for k in GRID_KEYS}
    grid["xg"] = grid_all.get("xg", list(np.polynomial.legendre.leggauss(int(grid["J"]))[0]))
    grid["wg"] = grid_all.get("wg", list(np.polynomial.legendre.leggauss(int(grid["J"]))[1]))
    keep = set(args.chains) if args.chains else None
    d = read_draws(out, keep)
    if keep is not None:
        missing = keep - set(d.chain.unique())
        if missing:
            raise SystemExit(f"--chains asked for {sorted(missing)}, which this run does not have")
    used = sorted(d.chain.unique())
    print(f"{len(d)} draws over {d.chain.nunique()} chains from {out}  (chains {used})")

    sfx = args.suffix
    corner(d, out / f"corner{sfx}.png")

    # the fit's object list, in the order the data were built, identified by
    # the one thing data.json records about it: how many epochs each object
    # contributed.  Guessing instead is how this went wrong -- a stage-3 run
    # takes stage3_sample() and a stage-4 subsample takes stratified(), both
    # are thirty objects long, and they overlap in one.  Every per-object
    # panel was then the wrong supernova's photometry against the right
    # index's parameters, and nothing in the figure said so.
    sample = analysis_sample()
    train, _ = training_split()
    n_sn = int(grid_all["N_sn"])
    n_ep_run = np.asarray(grid_all["stop"]) - np.asarray(grid_all["start"]) + 1
    names_fit = resolve_names(n_sn, n_ep_run, train, sample)
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

    lightcurves(d, names, index, sample, grid, args.n_draws,
                out / f"lightcurves{sfx}.png",
                span=(np.asarray(grid_all["start"]), np.asarray(grid_all["stop"])))


if __name__ == "__main__":
    main()
