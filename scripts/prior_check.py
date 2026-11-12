"""Check the prior terms Eq. (coefprior) and Eq. (hypernorm) introduce.

Those terms never reach the likelihood, so the transfer gate of
``model-definition.tex`` Section 5 --- ``sum ll_sn`` unchanged at the carried
point --- cannot see an error in any of them.  A wrong hyperplane
normalisation, a dropped Jacobian or a flipped sign would pass that gate and
silently bias ``a``, ``c``, ``tau_q`` and ``tau_psi``.

The checks exploit the fact that the new terms live on directions along which
the likelihood is *exactly* invariant, so the log density can be isolated from
it without a second implementation of the likelihood:

A. Adding ``lam a`` to ``b_til`` changes nothing the forward map reads --- that
   is what makes the direction redundant in the first place --- so the whole
   change in ``lp`` must be the unit prior on the discarded component.
   Likewise ``lam c`` added to ``d_til`` under ``Gs``.
B. Scaling ``z_c`` leaves ``c = ct / sqrt(c' Gs c)`` unchanged, so again the
   whole change must be the standard normal on ``z_c``.
C. ``tau_q`` and ``tau_psi`` enter no coefficient the forward map reads, only
   the widths, so scaling one gives a change in ``lp`` that is pure prior ---
   and its ``-(K - 2) log lambda`` term is the claim that tau now normalises
   ``K - 2`` coefficients rather than ``K - 1``.
D. The declared density in ``b_til`` and in ``d_til`` integrates to one,
   independent of ``a`` and ``c``.  This is what pins the hyperplane
   normalisation and the Jacobian together: with either wrong, the integral
   picks up a dependence on the hyperplane's orientation.

Run as ``python scripts/prior_check.py``.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dgsn.data import sample as S                     # noqa: E402
from dgsn.stan import hier_data                       # noqa: E402
from dgsn.stan.run import fixed_param, log_prob       # noqa: E402

PROGRAM = "hier_modes"
TOL = 1e-7


def build(n_per_fifth: int = 3, k: int = 8):
    names = S.stage3_sample(n_per_fifth=n_per_fifth)
    data = hier_data.build(names, k_psi=k, k_q=k)
    inits = hier_data.inits(names, k_psi=k, k_q=k, x1_ref=0.0)
    return data, inits, names


def unpack_w(params, k_q, k_psi):
    """z_a, b_til, z_c, d_til out of w_coef, with R the identity."""
    w = np.asarray(params["w_coef"], dtype=float)
    n_a = n_b = k_q - 1
    n_c = k_psi - 1
    return (w[:n_a], w[n_a:n_a + n_b], w[n_a + n_b:n_a + n_b + n_c],
            w[n_a + n_b + n_c:])


def repack_w(z_a, b_til, z_c, d_til):
    return np.concatenate([z_a, b_til, z_c, d_til])


def gauge(z_a, b_til, z_c, d_til, Gs, k_q, k_psi, a_m, p_m):
    kq = np.arange(1.0, k_q)
    kp = np.arange(1.0, k_psi)
    at = kq ** -(a_m + 0.5) * z_a
    ct = kp ** -(p_m + 0.5) * z_c
    a = at / np.sqrt(at @ at)
    b = b_til - (b_til @ a) * a
    c = ct / np.sqrt(ct @ Gs @ ct)
    d = d_til - (d_til @ Gs @ c) * c
    return a, b, c, d


def main() -> int:
    data, inits, names = build()
    k_q, k_psi = int(data["Kq"]), int(data["Kpsi"])
    Gs = np.asarray(fixed_param("gauge_gram", data)["gram_s"], dtype=float)
    p0 = {k: (np.asarray(v, dtype=float).tolist() if isinstance(v, (list, np.ndarray))
              else float(v)) for k, v in inits.items()}
    lp0, _ = log_prob(PROGRAM, data, p0)
    z_a, b_til, z_c, d_til = unpack_w(p0, k_q, k_psi)
    a_m, p_m = float(p0["alpha_q_m"]), float(p0["alpha_psi_m"])
    a, b, c, d = gauge(z_a, b_til, z_c, d_til, Gs, k_q, k_psi, a_m, p_m)
    print(f"{len(names)} objects, {data['N']} epochs, K = ({k_q}, {k_psi}); "
          f"lp at the start {lp0:.10f}")
    print(f"gauge at the start: |a| = {np.sqrt(a @ a):.12f}, b.a = {b @ a:+.3e}, "
          f"c'Gs c = {c @ Gs @ c:.12f}, d'Gs c = {d @ Gs @ c:+.3e}")

    def std_lpdf(x):
        return -0.5 * x * x - 0.5 * np.log(2 * np.pi)

    fails = []

    # ---- A: the two redundant coefficient directions ----
    print("\nA. the redundant directions: likelihood invariant, unit prior")
    for lam in (0.3, -1.7, 4.0):
        p = dict(p0)
        p["w_coef"] = repack_w(z_a, b_til + lam * a, z_c, d_til).tolist()
        lp, _ = log_prob(PROGRAM, data, p)
        s0, s1 = b_til @ a, (b_til + lam * a) @ a
        want = std_lpdf(s1) - std_lpdf(s0)
        ok = abs((lp - lp0) - want) < TOL
        fails.append(not ok)
        print(f"   b_til += {lam:+5.2f} a : dlp {lp - lp0:+14.9f}  want {want:+14.9f}  "
              f"{'ok' if ok else 'FAIL'}")
    for lam in (0.3, -1.7, 4.0):
        p = dict(p0)
        p["w_coef"] = repack_w(z_a, b_til, z_c, d_til + lam * c).tolist()
        lp, _ = log_prob(PROGRAM, data, p)
        s0, s1 = d_til @ Gs @ c, (d_til + lam * c) @ Gs @ c
        want = std_lpdf(s1) - std_lpdf(s0)
        ok = abs((lp - lp0) - want) < TOL
        fails.append(not ok)
        print(f"   d_til += {lam:+5.2f} c : dlp {lp - lp0:+14.9f}  want {want:+14.9f}  "
              f"{'ok' if ok else 'FAIL'}")

    # ---- B: the two scale directions the gauge removes ----
    print("\nB. the gauge-removed scales: likelihood invariant, z standard normal")
    for lam in (1.3, 0.6):
        p = dict(p0)
        p["w_coef"] = repack_w(lam * z_a, b_til, z_c, d_til).tolist()
        lp, _ = log_prob(PROGRAM, data, p)
        want = -0.5 * (lam ** 2 - 1.0) * (z_a @ z_a)
        ok = abs((lp - lp0) - want) < TOL
        fails.append(not ok)
        print(f"   z_a *= {lam:4.2f} : dlp {lp - lp0:+14.9f}  want {want:+14.9f}  "
              f"{'ok' if ok else 'FAIL'}")
        p = dict(p0)
        p["w_coef"] = repack_w(z_a, b_til, lam * z_c, d_til).tolist()
        lp, _ = log_prob(PROGRAM, data, p)
        want = -0.5 * (lam ** 2 - 1.0) * (z_c @ z_c)
        ok = abs((lp - lp0) - want) < TOL
        fails.append(not ok)
        print(f"   z_c *= {lam:4.2f} : dlp {lp - lp0:+14.9f}  want {want:+14.9f}  "
              f"{'ok' if ok else 'FAIL'}")

    # ---- C: tau scales only the widths, so -(K-2) log lambda ----
    print("\nC. tau: pure prior, and it normalises K-2 coefficients not K-1")
    kqv = np.arange(1.0, k_q)
    kpv = np.arange(1.0, k_psi)
    wq = float(p0["tau_q"]) * (kqv / np.exp(np.log(kqv).mean())) ** -(a_m + 0.5)
    wp = float(p0["tau_psi"]) * (kpv / np.exp(np.log(kpv).mean())) ** -(p_m + 0.5)
    for lam in (1.4, 0.7):
        for nm, w, coef, kk in (("tau_q", wq, b, k_q), ("tau_psi", wp, d, k_psi)):
            p = dict(p0)
            t0 = float(p0[nm])
            p[nm] = t0 * lam
            lp, _ = log_prob(PROGRAM, data, p)
            quad = -0.5 * (coef / w) @ (coef / w)
            half0 = -0.5 * (t0 / 10.0) ** 2
            half1 = -0.5 * (t0 * lam / 10.0) ** 2
            want = (quad * (1.0 / lam ** 2 - 1.0) - (kk - 2) * np.log(lam)
                    + half1 - half0)
            ok = abs((lp - lp0) - want) < TOL
            fails.append(not ok)
            print(f"   {nm} *= {lam:4.2f} : dlp {lp - lp0:+14.9f}  want {want:+14.9f}  "
                  f"{'ok' if ok else 'FAIL'}   [-(K-2)log lam = "
                  f"{-(kk - 2) * np.log(lam):+.6f}]")

    # ---- D: the declared density integrates to one, whatever the hyperplane ----
    print("\nD. normalisation: the declared density in b_til and d_til integrates")
    print("   to one independent of the hyperplane's orientation (Monte Carlo,")
    print("   importance sampling from the unconstrained Eq. (coefprior) normal)")
    rng = np.random.default_rng(11)
    M = 400000
    for trial in range(3):
        if trial == 0:
            av, cv = a, c
        else:
            za = rng.standard_normal(k_q - 1)
            zc = rng.standard_normal(k_psi - 1)
            av, _, cv, _ = gauge(za, b_til, zc, d_til, Gs, k_q, k_psi, a_m, p_m)
        # b_til: proposal N(0, diag(wq^2)); target = Eq.(coefprior) on the
        # projection, normalised on the hyperplane, times a unit normal on the
        # discarded component.  Jacobian is 1 here.
        X = rng.standard_normal((M, k_q - 1)) * wq
        s = X @ av
        B = X - np.outer(s, av)
        logt = (-0.5 * ((B / wq) ** 2).sum(1) - np.log(wq).sum()
                + 0.5 * np.log(((av * wq) ** 2).sum())
                - 0.5 * (k_q - 2) * np.log(2 * np.pi)
                + std_lpdf(s))
        logp = (-0.5 * ((X / wq) ** 2).sum(1) - np.log(wq).sum()
                - 0.5 * (k_q - 1) * np.log(2 * np.pi))
        Ib = np.exp(logt - logp).mean()
        eb = np.exp(logt - logp).std(ddof=1) / np.sqrt(M)
        # d_til: the same with n = Gs c, plus the Jacobian log||Gs c||
        n = Gs @ cv
        nn = np.sqrt(n @ n)
        Y = rng.standard_normal((M, k_psi - 1)) * wp
        sd = Y @ n
        D = Y - np.outer(sd, cv)
        logt = (-0.5 * ((D / wp) ** 2).sum(1) - np.log(wp).sum()
                + 0.5 * np.log((((n / nn) * wp) ** 2).sum()) + np.log(nn)
                - 0.5 * (k_psi - 2) * np.log(2 * np.pi)
                + std_lpdf(sd))
        logp = (-0.5 * ((Y / wp) ** 2).sum(1) - np.log(wp).sum()
                - 0.5 * (k_psi - 1) * np.log(2 * np.pi))
        Id = np.exp(logt - logp).mean()
        ed = np.exp(logt - logp).std(ddof=1) / np.sqrt(M)
        okb, okd = abs(Ib - 1) < 5 * eb + 2e-3, abs(Id - 1) < 5 * ed + 2e-3
        fails.append(not (okb and okd))
        print(f"   draw {trial}: |Gs c| = {nn:7.4f}   "
              f"int(b_til) = {Ib:.5f} +- {eb:.5f} {'ok' if okb else 'FAIL'}   "
              f"int(d_til) = {Id:.5f} +- {ed:.5f} {'ok' if okd else 'FAIL'}")

    bad = sum(fails)
    print(f"\n{len(fails) - bad} of {len(fails)} checks passed")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
