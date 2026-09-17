"""The forward map of ``forward.stan`` at fixed parameters: the pins that must
hold exactly, the scalings that have closed forms, and the bases."""

import numpy as np
import pytest

from dgsn.stan import data as D
from dgsn.stan.run import fixed_param


def run_forward(m_per_side=D.M_PER_SIDE, k_psi=D.K_PSI, k_q=D.K_Q, globals_=None,
                params=None, **points):
    data = D.grid_data(k_psi=k_psi, k_q=k_q, m_per_side=m_per_side)
    data.update(globals_ or D.zero_globals(k_psi, k_q))
    data.update(params or D.params_data())
    data.update(z=0.0)
    data.update(D.eval_points(**points))
    return fixed_param("forward", data)


def random_globals(k_psi=D.K_PSI, k_q=D.K_Q, seed=0):
    """Smooth-looking coefficients: decaying with mode index, gauge-fixed."""
    rng = np.random.default_rng(seed)
    kq, kp = np.arange(1, k_q), np.arange(1, k_psi)
    a = rng.normal(size=k_q - 1) * kq ** -1.5
    b = rng.normal(size=k_q - 1) * kq ** -1.5
    c = rng.normal(size=k_psi - 1) * kp ** -1.5
    d = rng.normal(size=k_psi - 1) * kp ** -1.5
    a /= np.linalg.norm(a); b -= (b @ a) * a
    c /= np.linalg.norm(c); d -= (d @ c) * c
    return dict(a=a, b=b, c=c, d=d)


def test_anchor_exact():
    """s(0) = 0, gamma(0) = 0, psi(0) = -pi/2, exactly, at arbitrary globals."""
    out = run_forward(globals_=random_globals(), params=D.params_data(x_s=-1.8, x_p=0.7, x_psi=-0.4),
                      p_eval=[0.0], s_eval=[0.0])
    assert out["s_of_p"][0] == 0.0
    assert out["gamma"][0, 0] == 0.0 and out["gamma"][0, 1] == 0.0
    assert out["psi"][0] == -np.pi / 2


@pytest.mark.parametrize("m_per_side", [D.M_PER_SIDE, 2 * D.M_PER_SIDE])
def test_xs_scales_arclength(m_per_side):
    """x_s scales s(p) by e^{x_s}, at every phase, at arbitrary globals."""
    g = random_globals()
    p = np.linspace(-20, 45, 27)
    s0 = run_forward(m_per_side, globals_=g, params=D.params_data(x_s=0.0, x_p=0.3), p_eval=p)["s_of_p"]
    s1 = run_forward(m_per_side, globals_=g, params=D.params_data(x_s=-1.8, x_p=0.3), p_eval=p)["s_of_p"]
    rel = np.abs(s1 - np.exp(-1.8) * s0) / np.abs(np.exp(-1.8) * s0).max()
    print(f"x_s scaling m={m_per_side}: max relative error {rel.max():.2e}")
    assert rel.max() < 1e-14


@pytest.mark.parametrize("m_per_side", [D.M_PER_SIDE, 2 * D.M_PER_SIDE])
def test_unit_rate_is_identity(m_per_side):
    """a = b = 0, x_s = 0: s = p on the unclamped map."""
    p = np.linspace(-20, 45, 27)
    out = run_forward(m_per_side, p_eval=p)
    assert np.abs(out["s_of_p"] - p).max() < 1e-12


def test_bases_orthonormal():
    """Both bases orthonormal under their own flat measure; every turning mode
    zero at s = 0; the transform triangular; the shifted Gram conditioned as
    the note says.

    The inner products are taken with a 64-point Gauss-Legendre rule over the
    whole interval, exact for products of the modes.  The model's own grid is
    not: 8 points per panel integrate a product of two degree-23 polynomials
    only to 7e-6, which is a fact about that grid and not about the basis."""
    x, w = np.polynomial.legendre.leggauss(64)
    ulo, uhi = np.arcsinh(np.array(D.ARC_RANGE) / D.DELTA)
    vlo, vhi = np.arcsinh(np.array(D.MODEL_WINDOW) / D.DELTA_P)
    u = 0.5 * (uhi + ulo) + 0.5 * (uhi - ulo) * x
    v = 0.5 * (vhi + vlo) + 0.5 * (vhi - vlo) * x
    out = run_forward(s_eval=np.concatenate([[0.0], D.DELTA * np.sinh(u)]),
                      p_eval=D.DELTA_P * np.sinh(v))
    Xu, Xv = out["Xu"][1:], out["Xv"]
    Gu = Xu.T @ (0.5 * (uhi - ulo) * w[:, None] * Xu)
    Gv = Xv.T @ (0.5 * (vhi - vlo) * w[:, None] * Xv)
    eu = np.abs(Gu - np.eye(D.K_PSI - 1)).max()
    ev = np.abs(Gv - np.eye(D.K_Q - 1)).max()
    print(f"orthonormality: turning {eu:.2e}, traversal {ev:.2e}")
    assert eu < 1e-12 and ev < 1e-12
    assert np.abs(out["Xu"][0]).max() == 0.0          # vanishes at s = 0, exactly
    T = out["T_out"]
    assert np.abs(np.tril(T, -1)).max() == 0.0        # upper triangular: ascending degree
    # T = R^{-T} with R R^T the Gram matrix: cond(Gram) = cond(T)^2
    cond_gram = np.linalg.cond(T) ** 2
    print(f"shifted Gram condition number at K_psi = {D.K_PSI}: {cond_gram:.1f}")
    assert 10 < cond_gram < 20


def test_turning_gram_under_ds():
    """The Gram matrix of the turning basis under ds, exported by
    forward.stan, against the same integral taken with a 400-point
    Gauss-Legendre rule in u on the exported basis."""
    x, w = np.polynomial.legendre.leggauss(400)
    ulo, uhi = np.arcsinh(np.array(D.ARC_RANGE) / D.DELTA)
    u = 0.5 * (uhi + ulo) + 0.5 * (uhi - ulo) * x
    out = run_forward(s_eval=D.DELTA * np.sinh(u))
    Xu = out["Xu"]
    G = Xu.T @ ((0.5 * (uhi - ulo) * w * D.DELTA * np.cosh(u))[:, None] * Xu)
    err = np.abs(out["Gs_out"] - G).max() / np.abs(G).max()
    print(f"G_s: max relative error {err:.2e}; diagonal range [{np.diag(G).min():.2e}, {np.diag(G).max():.2e}]")
    assert err < 1e-10
    assert np.allclose(out["Gs_out"], out["Gs_out"].T)


def test_kappa0_matches_finite_difference():
    """psi'(0) from the Legendre derivative against a difference of psi."""
    g = random_globals()
    h = 1e-7
    out = run_forward(globals_=g, params=D.params_data(x_psi=0.5), s_eval=[-h, h])
    fd = (out["psi"][1] - out["psi"][0]) / (2 * h)
    print(f"kappa0 {out['kappa0']:.6f} vs finite difference {fd:.6f}")
    assert abs(out["kappa0"] - fd) < 1e-5 * max(1.0, abs(fd))


def test_flux_matches_magnitudes():
    """Eq. (predflux): the flux at an epoch is 10^{-0.4(m - zp)} of the
    magnitude the curve gives there, with the placement added, inside the
    window where the taper is unity."""
    g = random_globals()
    par = D.params_data(m_peak=18.0, colr=0.1, x_s=-1.8, x_p=0.2, x_psi=-0.3, delta_t=0.5)
    dt = np.array([0.5, 3.5, -4.5, 20.5])
    band = np.array([1, 2, 1, 2])
    out = run_forward(globals_=g, params=par, dt=dt, band=band, p_eval=dt - 0.5)
    # The same s by both routes, up to the phase clamp's exponential tail:
    # eps e^{-w/eps} of the window, 7e-11 d at p = -5, 3e-4 d at the data
    # window edge, which is what "exponential accuracy inside" costs.
    assert np.abs(out["lf"][:, 2] - out["s_of_p"]).max() < 1e-10
    s_hat = out["lf"][:, 3]                                   # after the arclength clamp
    print(f"clamp displacement: {np.abs(s_hat - out['s_of_p']).max():.2e} mag")
    out2 = run_forward(globals_=g, params=par, s_eval=s_hat)
    mag = np.where(band == 1, out2["gamma"][:, 0] + 18.0, out2["gamma"][:, 1] + 18.0 - 0.1)
    f_expected = 10 ** (-0.4 * (mag - D.ZP))
    assert np.abs(out["f_hat"] / f_expected - 1).max() < 1e-13
    assert np.abs(out["lf"][:, 1]).max() < 1e-12       # log taper is zero inside
