"""Read the arms of scripts/scale_ladder.py side by side.

Reads partial CSVs, so it works on a running ladder; the draw count is
printed with every arm because a comparison across unequal draws is not one.

The figure of merit is leapfrogs per effective draw, not wall time and not
ESS: the arms share a machine and contend for it, so seconds are not
comparable, while a leapfrog is the same unit of work in every arm.  The
reference is runs/stage3-common-split, which is shape 0 at these settings.
"""

from __future__ import annotations

import glob
from pathlib import Path

import numpy as np

from dgsn.stan.diag import ess_bulk, read_csvs, rhat, rhat_plain

# the name a quantity has in each program: hier_common samples the scales, so
# calls them tau_q; hier_common_scale may hold them, so names the resolved
# value tau_q_eff.  alpha_q and alpha_psi carry the same name in both.
ALIAS = {"tau_q": ("tau_q", "tau_q_eff"), "tau_psi": ("tau_psi", "tau_psi_eff")}
SAMPLER = ["lp__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__"]
SCALARS = ["tau_q", "tau_psi", "alpha_q", "alpha_psi", "kappa0"]


def _read(run: Path, k_psi: int, k_q: int):
    csvs = sorted(glob.glob(str(run / "*-[0-9].csv")))
    if not csvs:
        return None
    header = None
    for line in open(csvs[0]):
        if line.startswith("lp__"):
            header = set(line.strip().split(","))
            break
    if header is None:
        return None
    want = list(SAMPLER)
    resolved = {}
    for key in SCALARS:
        for cand in ALIAS.get(key, (key,)):
            if cand in header:
                resolved[key] = cand
                want.append(cand)
                break
    dk = [f"d_til.{k}" for k in range(1, k_psi - 1 + 1) if f"d_til.{k}" in header]
    bk = [f"b_til.{k}" for k in range(1, k_q - 1 + 1) if f"b_til.{k}" in header]
    try:
        draws, n = read_csvs(csvs, want + dk + bk)
    except ValueError:
        return None
    out = {key: draws[col] for key, col in resolved.items() if col in draws}
    out.update({k: draws[k] for k in SAMPLER if k in draws})
    out["||d||"] = np.sqrt(sum(draws[k] ** 2 for k in dk))
    out["||b||"] = np.sqrt(sum(draws[k] ** 2 for k in bk))
    return out, n, len(csvs)


def _held(v: np.ndarray) -> bool:
    return bool(np.var(v) <= 0)


def report(arms, reference: Path, k_psi: int, k_q: int):
    runs = [("shape0 (reference)", Path(reference))]
    runs += [(tag, Path(f"runs/ladder-{tag}")) for tag, _, _ in arms]
    got = []
    for tag, run in runs:
        r = _read(run, k_psi, k_q)
        if r is None:
            print(f"{tag:20s} {run}: no draws yet")
            continue
        d, n, nch = r
        got.append((tag, run, d, n, nch))
    if not got:
        print("nothing to report")
        return

    print(f"\n{'arm':20s} {'draws':>12s} {'stepsize':>22s} {'depth':>6s} {'leapfrog':>9s} {'div':>5s}")
    for tag, run, d, n, nch in got:
        ss = d["stepsize__"][:, 0]
        print(f"{tag:20s} {nch} x {n:<8d} "
              f"{np.array2string(ss, precision=4, max_line_width=200):>22s} "
              f"{d['treedepth__'].mean():6.2f} {d['n_leapfrog__'].mean():9.0f} "
              f"{int(d['divergent__'].sum()):5d}")
    for tag, run, d, n, nch in got:
        ss = d["stepsize__"][:, 0]
        print(f"{'':20s} {tag}: stepsize spread {ss.max()/ss.min():.2f}x, "
              f"{100*np.mean(d['treedepth__'] >= d['treedepth__'].max()):.0f}% at the "
              f"deepest tree it reached ({int(d['treedepth__'].max())})")

    print(f"\n{'quantity':12s}" + "".join(f"{t:>26s}" for t, *_ in got))
    print(f"{'':12s}" + "".join(f"{'mean  R_hat  ESS  lf/eff':>26s}" for _ in got))
    for key in SCALARS + ["||b||", "||d||", "lp__"]:
        row = f"{key:12s}"
        for tag, run, d, n, nch in got:
            if key not in d:
                row += f"{'--':>26s}"
                continue
            v = d[key]
            lf = d["n_leapfrog__"].sum()
            if _held(v):
                row += f"{f'{v.ravel()[0]:.4g} held':>26s}"
            else:
                e = ess_bulk(v)
                row += (f"{v.mean():9.4g}{rhat(v):7.3f}{e:7.1f}"
                        f"{lf/e:9.0f}" if np.isfinite(e) and e > 0
                        else f"{v.mean():9.4g}{rhat(v):7.3f}{'--':>7s}{'--':>9s}")
        print(row)

    print("\nthe funnel, and whether each arm removed it:")
    for tag, run, d, n, nch in got:
        if "tau_psi" not in d:
            continue
        tau, nd = d["tau_psi"].ravel(), d["||d||"].ravel()
        if _held(tau):
            print(f"  {tag:20s} tau_psi held at {tau[0]:.4g}: the direction is gone, "
                  f"so there is no correlation to measure")
        else:
            c = np.corrcoef(np.log(tau), np.log(nd))[0, 1]
            print(f"  {tag:20s} corr(log tau_psi, log||d||) = {c:+.4f}")

    print("\nwithin-chain stationarity (plain split-Rhat per chain; drift fails on its own):")
    for tag, run, d, n, nch in got:
        parts = []
        for key in ("lp__", "||d||"):
            if key in d and not _held(d[key]):
                parts.append(f"{key} " + " ".join(
                    f"{rhat_plain(d[key][i:i+1]):.2f}" for i in range(d[key].shape[0])))
        print(f"  {tag:20s} " + ";  ".join(parts))
