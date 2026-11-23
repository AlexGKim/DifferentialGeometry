"""Split R-hat, reported separately for global and per-supernova coordinates.

The convergence gate for this pair of runs is on the GLOBAL parameters only,
at R-hat < 1.05.  Per-object coordinates are reported but do not gate.
"""
import glob, sys
import numpy as np

# Every block that is one-per-supernova.  Anything else is global.
PER_SN = ("m_peak", "colr", "c_dev", "x_s", "x_s_tot", "d_dev", "delta_t",
          "x_p", "x_psi", "x_p_dev", "x_psi_dev", "x_s_dev",
          "x_p_rel", "x_psi_rel", "u_p", "u_psi", "u_del",
          "lp_sn", "ll_sn", "chi2_sn", "s_pmin", "s_pmax", "s_epoch", "resid_z")


def load(path, warm=500):
    names, rows = None, []
    for ln in open(path):
        if ln.startswith("#"):
            continue
        p = ln.rstrip("\n").split(",")
        if names is None:
            names = p; continue
        if len(p) != len(names):
            continue
        rows.append([float(x) for x in p])
    return names, np.asarray(rows)[warm:]


def split_rhat(ch):
    m, n = ch.shape
    h = n // 2
    s = np.concatenate([ch[:, :h], ch[:, h:2 * h]], 0)
    w = s.var(1, ddof=1).mean()
    b = h * s.mean(1).var(ddof=1)
    return np.sqrt(((h - 1) / h * w + b / h) / w) if w > 0 else np.nan


def main(run, prog, bar=1.05):
    fs = sorted(glob.glob(f"{run}/{prog}-*.csv"))
    if not fs:
        print(f"no chains in {run}"); return 2
    dat = [load(f) for f in fs]
    names = dat[0][0]
    n = min(d[1].shape[0] for d in dat)
    if n < 50:
        print(f"only {n} post-warmup draws; too few to judge"); return 2
    X = np.stack([d[1][:n] for d in dat])
    glo, per = [], []
    for j, nm in enumerate(names):
        if nm.endswith("__") or X[:, :, j].std() == 0:
            continue
        r = split_rhat(X[:, :, j])
        if not np.isfinite(r):
            continue
        (per if nm.split(".")[0] in PER_SN else glo).append((r, nm))
    glo.sort(reverse=True); per.sort(reverse=True)
    g = np.array([r for r, _ in glo]); p = np.array([r for r, _ in per])
    print(f"{run}: {[d[1].shape[0] for d in dat]} post-warmup, common window {n}")
    print(f"  GLOBAL   {len(g):5d} coords  max {g.max():.4f}  "
          f">1.05 {(g > 1.05).sum()}  >1.01 {(g > 1.01).sum()}")
    print("    worst:", ", ".join(f"{nm} {r:.3f}" for r, nm in glo[:8]))
    if len(p):
        print(f"  PER-SN   {len(p):5d} coords  max {p.max():.4f}  "
              f">1.05 {(p > 1.05).sum()}  >1.01 {(p > 1.01).sum()}  (does not gate)")
        print("    worst:", ", ".join(f"{nm} {r:.3f}" for r, nm in per[:6]))
    ok = g.max() < bar
    print(f"\n  GATE: global max R-hat {g.max():.4f} {'<' if ok else '>='} {bar}  "
          f"-> {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2],
                  float(sys.argv[3]) if len(sys.argv) > 3 else 1.05))
