"""The floors of Section 3, pushed through the truncated model before any
fit: the K-correction at z = 0.05 and the Galactic extinction drift at
E(B-V) = 0.3, each as what it does to x_psi and to the per-band residual
that the fit would have to absorb, Eq. (like) carrying no scatter of its own.

    python scripts/floors.py

Each is a difference between two fits of noise-free Hsiao light curves, so
that the template's own mismatch with the SALT2-fitted globals cancels.
"""

import numpy as np

from dgsn.data.sample import median_x1
from dgsn.salt.globals import fit_globals
from dgsn.salt.tracks import _model, g_maximum, hsiao_epochs
from dgsn.stan import data as D
from dgsn.stan.run import fixed_param, optimize

NAMES = ["m_peak", "colr", "x_s", "x_p", "x_psi", "delta_t"]


def fit_epochs(G, epochs, start):
    data = D.grid_data()
    data.update({k: G[k] for k in ("a", "b", "c", "d")})
    data.update(D.epochs_data(epochs["dt"], epochs["f"], epochs["sig_f"], epochs["band"],
                              epochs["z"]))
    data.update(D.prior_data(prior_on=False))
    fit = optimize("single_sn", data, start)
    fp = {k: fit[k] for k in NAMES}
    fdata = D.grid_data()
    fdata.update({k: G[k] for k in ("a", "b", "c", "d")})
    fdata.update(D.params_data(**fp))
    fdata.update(z=epochs["z"])
    fdata.update(D.eval_points(dt=epochs["dt"], band=epochs["band"]))
    f_hat = fixed_param("forward", fdata)["f_hat"]
    resid = -2.5 * np.log10(f_hat / epochs["f"])
    rms = {b: float(np.sqrt(np.mean(resid[epochs["band"] == i] ** 2))) for i, b in ((1, "g"), (2, "r"))}
    return fp, rms


def report(label, base, test):
    (p0, r0), (p1, r1) = base, test
    print(f"\n{label}")
    print("  parameter shifts:  " + "  ".join(f"{k} {p1[k] - p0[k]:+.4f}" for k in NAMES))
    print(f"  residual rms (mag): g {r0['g']:.4f} -> {r1['g']:.4f}   r {r0['r']:.4f} -> {r1['r']:.4f}")


def main():
    x1_ref = median_x1()
    G = fit_globals(D.K_PSI, D.K_Q, x1_ref)
    print(f"globals from SALT2 at x1 = {x1_ref:.3f}: x_s {G['x_s']:.3f}, reference at x_p {G['x_p']:.3f}, x_psi {G['x_psi']:.3f}")
    start = D.params_data(m_peak=18.0, colr=0.0, x_s=G["x_s"], x_p=G["x_p"], x_psi=G["x_psi"])

    base = fit_epochs(G, hsiao_epochs(z=0.0), start)
    print("\nHsiao at z = 0 against the SALT2-fitted globals (the template mismatch, which the floors are measured against):")
    print("  " + "  ".join(f"{k} {base[0][k]:.4f}" for k in NAMES) + f"   rms g {base[1]['g']:.4f} r {base[1]['r']:.4f}")

    kcorr = fit_epochs(G, hsiao_epochs(z=0.05), dict(start, m_peak=base[0]["m_peak"]))
    report("K-correction floor: Hsiao observed at z = 0.05 (phases de-redshifted by the model) vs z = 0", base, kcorr)

    # Extinction: the drift is what remains after the constant translation.
    ebv = 0.3
    dusty = hsiao_epochs(z=0.0, ebv=ebv)
    m_d, m_0 = _model("hsiao", ebv=ebv), _model("hsiao")
    m_d.set(z=0.0, t0=0.0, mwebv=ebv)
    m_0.set(z=0.0, t0=0.0)
    m_d.set_source_peakmag(18.0, "ztfg", "ab")
    m_0.set_source_peakmag(18.0, "ztfg", "ab")
    t_g = g_maximum(m_0)
    A = {b: -2.5 * np.log10(m_d.bandflux(n, t_g) / m_0.bandflux(n, t_g)) for b, n in ((1, "ztfg"), (2, "ztfr"))}
    print(f"\nextinction at peak for E(B-V) = {ebv}: A_g {A[1]:.4f}, A_r {A[2]:.4f}")
    clean = hsiao_epochs(z=0.0)
    translated = dict(clean, f=clean["f"] * 10 ** (-0.4 * np.where(clean["band"] == 1, A[1], A[2])))
    fit_t = fit_epochs(G, translated, dict(start, m_peak=18.0 + A[1]))
    fit_d = fit_epochs(G, dusty, dict(start, m_peak=18.0 + A[1]))
    report(f"Extinction drift floor: Hsiao with F99 dust at E(B-V) = {ebv} vs the same translated by the peak (A_g, A_r)", fit_t, fit_d)


if __name__ == "__main__":
    main()
