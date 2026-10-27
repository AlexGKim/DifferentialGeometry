"""Light-curve fits from each of the two posterior modes, side by side.

    python scripts/plot_modes.py --out runs/stage3-common-halfnorm

``runs/stage3-common-halfnorm`` split two-and-two: chains 2 and 4 sit at
tau_psi ~ 2.6 with ||d|| ~ 14, chains 1 and 3 at tau_psi ~ 22 with ||d|| ~ 130
and a higher lp__.  This asks the only question that settles which is a
supernova: what does each solution say the light curve does, against the data
the likelihood sees.

Chains are grouped by their own tau_psi --- the largest gap between the sorted
per-chain medians --- rather than by a hard-coded chain number, so the script
survives a rerun that lands the modes on different chains.

The model is never evaluated in Python.  Every curve comes from
``forward.stan`` at one draw's parameters through ``fixed_param``, reusing
``plot_fit.evaluate``.  Flux is native, as the archive gives it and the
likelihood sees it, and the phase axis is observed days from the catalogued
t0, which is the one axis both modes share: each carries its own delta_t, so
a rest-frame phase axis would move the data under the curves.
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))

from dgsn.data.epochs import likelihood_epochs  # noqa: E402
from dgsn.data.sample import analysis_sample, stage3_sample, training_split  # noqa: E402
from plot_fit import GRID_KEYS, evaluate, read_draws  # noqa: E402

LOW, HIGH = "#1f77b4", "#d62728"        # the two modes
BAND = {1: "g", 2: "r"}


def split_chains(d, key="tau_psi"):
    """The two groups, by the largest gap in the sorted per-chain medians."""
    med = {int(c): float(np.median(d[d.chain == c][key])) for c in sorted(d.chain.unique())}
    order = sorted(med, key=lambda c: med[c])
    vals = [med[c] for c in order]
    gaps = np.diff(vals)
    if len(gaps) == 0:
        raise SystemExit("one chain only; nothing to compare")
    cut = int(np.argmax(gaps)) + 1
    return order[:cut], order[cut:], med


def objects_of(grid_all, sample):
    """The fit's object list, in the order the data were built.

    Stage 3 and stage 4 pick their objects differently --- ``stage3_sample``
    takes six per fifth of x1, ``stratified`` subsets the training split ---
    and a run records neither choice, only the resulting epoch counts.  Taking
    the wrong list silently pairs each supernova's parameters with another's
    photometry, so every candidate is checked against ``start``/``stop`` from
    the run's own data before it is used."""
    train, _ = training_split()
    n_sn = int(grid_all["N_sn"])
    want = np.asarray(grid_all["stop"]) - np.asarray(grid_all["start"]) + 1

    def candidates():
        yield "stage3_sample", stage3_sample()
        yield "training_split", train
        from stage4_full import stratified
        yield "stratified", stratified(train, n_sn, sample)

    tried = []
    for label, names in candidates():
        if len(names) != n_sn:
            tried.append(f"{label} ({len(names)} objects)")
            continue
        got = np.array([len(likelihood_epochs(nm, sample.loc[nm].t0,
                                              sample.loc[nm].redshift)["dt"])
                        for nm in names])
        if np.array_equal(got, want):
            print(f"object list: {label}, epoch counts match data.json")
            return list(names)
        tried.append(f"{label} ({int((got != want).sum())} of {n_sn} epoch counts differ)")
    raise SystemExit("no candidate object list matches the run's epoch counts: "
                     + "; ".join(tried))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/stage3-common-halfnorm")
    ap.add_argument("--n-draws", type=int, default=8, help="draws per mode per object")
    ap.add_argument("--n-objects", type=int, default=4)
    ap.add_argument("--names", nargs="*", default=None)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    out = Path(args.out)
    grid_all = json.load(open(out / "data.json"))
    grid = {k: grid_all[k] for k in GRID_KEYS}
    grid["xg"] = grid_all.get("xg", list(np.polynomial.legendre.leggauss(int(grid["J"]))[0]))
    grid["wg"] = grid_all.get("wg", list(np.polynomial.legendre.leggauss(int(grid["J"]))[1]))
    k_psi, k_q = int(grid["Kpsi"]), int(grid["Kq"])

    d = read_draws(out)
    lo_ch, hi_ch, med = split_chains(d)
    print(f"{len(d)} draws, {d.chain.nunique()} chains from {out}")
    print("per-chain median tau_psi: " + ", ".join(f"ch{c} {med[c]:.2f}" for c in sorted(med)))
    print(f"low-tau_psi group: chains {lo_ch}   high-tau_psi group: chains {hi_ch}")
    groups = [("low", lo_ch, LOW), ("high", hi_ch, HIGH)]

    sample = analysis_sample()
    names_fit = objects_of(grid_all, sample)
    n_sn = int(grid_all["N_sn"])
    n_ep = np.asarray(grid_all["stop"]) - np.asarray(grid_all["start"]) + 1

    # -2 lp per epoch in each mode, and the objects where they disagree most
    per = {}
    for label, chs, _ in groups:
        sub = d[d.chain.isin(chs)]
        per[label] = np.array([-2 * np.median(sub[f"lp_sn.{i}"]) / n_ep[i - 1]
                               for i in range(1, n_sn + 1)])
    gap = per["high"] - per["low"]
    print(f"\n-2 lp/epoch, summed over the fit: low {(per['low']*n_ep).sum():.0f}, "
          f"high {(per['high']*n_ep).sum():.0f}  (lower is a better fit)")
    rank = np.argsort(-np.abs(gap))
    if args.names:
        names = args.names
    else:
        names = [names_fit[i] for i in rank[:args.n_objects]]
        print("objects: largest |high - low| in -2 lp/epoch -> " + ", ".join(names))
    print("\n  object            x1      z   -2lp/ep low   high    high-low")
    for nm in names:
        i = names_fit.index(nm)
        print(f"  {nm:16s} {sample.loc[nm].x1:+5.2f} {sample.loc[nm].redshift:.4f}"
              f"   {per['low'][i]:9.2f} {per['high'][i]:6.2f} {gap[i]:+10.2f}")

    p_eval = np.linspace(grid["pmin"] + 0.5, grid["pmax"] - 0.5, 160)
    s_eval = np.linspace(grid["gsmin"], grid["gsmax"], 220)
    rng = np.random.default_rng(args.seed)

    fig, ax = plt.subplots(len(names), 3, figsize=(18.5, 3.7 * len(names)))
    ax = np.atleast_2d(ax)
    for row, name in enumerate(names):
        r = sample.loc[name]
        ep = likelihood_epochs(name, r.t0, r.redshift)
        idx = names_fit.index(name) + 1
        at, al, ar = ax[row, 0], ax[row, 1], ax[row, 2]
        for label, chs, colour in groups:
            sub = d[d.chain.isin(chs)].reset_index(drop=True)
            pick = rng.choice(len(sub), size=min(args.n_draws, len(sub)), replace=False)
            fobs = []
            for j, i in enumerate(pick):
                draw = sub.iloc[i]
                cv = evaluate(grid, draw, idx, ep["z"], k_psi, k_q, p_eval, s_eval,
                              dt_obs=ep["dt"], band_obs=ep["band"])
                fobs.append(cv["f_obs"])
                # observed days from the catalogued t0: the axis both modes share
                dt = p_eval * (1 + ep["z"]) + float(draw[f"delta_t.{idx}"])
                at.plot(cv["gamma"][:, 0] + cv["m_peak"],
                        cv["gamma"][:, 1] + cv["m_peak"] - cv["colr"],
                        "-", lw=0.8, alpha=0.45, color=colour,
                        label=f"$\\tau_\\psi$ {label}" if j == 0 else None)
                al.plot(dt, cv["f_g"], "-", lw=0.9, alpha=0.5, color=colour)
                al.plot(dt, cv["f_r"], "--", lw=0.9, alpha=0.5, color=colour)
            # residual at the posterior-median model, in units of the ZTF error
            fhat = np.median(np.vstack(fobs), axis=0)
            res = (ep["f"] - fhat) / ep["sig_f"]
            for bnd, mk in ((1, "o"), (2, "s")):
                m = ep["band"] == bnd
                ar.plot(ep["dt"][m], res[m], mk, ms=3.5, alpha=0.75, color=colour,
                        label=f"{label} {BAND[bnd]}  $\\chi^2$/ep {np.mean(res[m]**2):.1f}")
        at.invert_xaxis(); at.invert_yaxis()
        at.set_xlabel("$m_g$"); at.set_ylabel("$m_r$")
        at.set_title(f"{name}: track through magnitude space", fontsize=10)
        at.legend(fontsize=8, frameon=False)

        for bnd, mk in ((1, "o"), (2, "s")):
            m = ep["band"] == bnd
            al.errorbar(ep["dt"][m], ep["f"][m], yerr=ep["sig_f"][m], fmt=mk, ms=3.5,
                        color="0.25" if bnd == 1 else "0.55", lw=0.8, capsize=0, zorder=3,
                        label=f"{BAND[bnd]} data ({int(m.sum())} ep)")
        al.axhline(0, color="0.7", lw=0.6)
        al.set_xlim(ep["dt"].min() - 5, ep["dt"].max() + 5)
        al.set_xlabel("observed days from catalogued $t_0$")
        al.set_ylabel("flux (zp 30)")
        i0 = idx - 1
        al.set_title(f"$x_1$ {r.x1:+.2f}, $z$ {r.redshift:.4f};  "
                     f"$-2\\,$lp/epoch  low {per['low'][i0]:.1f}  high {per['high'][i0]:.1f}"
                     f"   (solid $g$, dashed $r$)", fontsize=10)
        al.legend(fontsize=8, frameon=False)

        for y in (-1, 1):
            ar.axhline(y, color="0.8", lw=0.6, ls=":")
        ar.axhline(0, color="0.6", lw=0.6)
        ar.set_xlim(*al.get_xlim())
        ar.set_xlabel("observed days from catalogued $t_0$")
        ar.set_ylabel(r"$(f - \hat f)\,/\,\sigma_f$")
        ar.set_title("residual at the median model (circles $g$, squares $r$)", fontsize=10)
        ar.legend(fontsize=7, frameon=False, ncol=2)
    fig.suptitle("The two $\\tau_\\psi$ solutions against the same data, in native flux",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    path = out / "modes_lightcurves.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
