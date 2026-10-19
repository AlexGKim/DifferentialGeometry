"""The quadrature, clamps and taper on closed forms, through ``checks.stan``.

Each quadrature check runs at two resolutions, ``m`` and ``2m`` panels per
side: one number is a claim, two that agree are a measurement.
"""

import numpy as np
import pytest

from dgsn.stan import data as D
from dgsn.stan.run import fixed_param

K_CIRCLE = 0.3   # rad/mag: the turning rate of a SALT2 track away from the anchor


def run_checks(m_per_side, k=K_CIRCLE, s_eval=(), p_eval=(), x_eval=(), w_eval=()):
    data = D.grid_data(m_per_side=m_per_side)
    s_eval, p_eval = np.atleast_1d(s_eval), np.atleast_1d(p_eval)
    x_eval, w_eval = np.atleast_1d(x_eval), np.atleast_1d(w_eval)
    data.update(k=k, Ns=len(s_eval), s_eval=s_eval, Np=len(p_eval), p_eval=p_eval,
                Nx=len(x_eval), x_eval=x_eval, Nw=len(w_eval), w_eval=w_eval)
    return fixed_param("checks", data)


@pytest.mark.parametrize("m_per_side", [D.M_PER_SIDE, 2 * D.M_PER_SIDE])
def test_circle(m_per_side):
    """Constant psi' = k gives gamma(s) = ((1 - cos ks)/k, -sin(ks)/k)."""
    s = np.linspace(D.ARC_RANGE[0] + 0.1, D.ARC_RANGE[1] - 0.1, 41)
    out = run_checks(m_per_side, s_eval=s)
    exact = np.stack([(1 - np.cos(K_CIRCLE * s)) / K_CIRCLE,
                      -np.sin(K_CIRCLE * s) / K_CIRCLE], axis=1)
    err = np.abs(out["circle"] - exact).max()
    print(f"circle m={m_per_side}: max error {err:.2e}")
    # Measured on the re-scaled basis ([-8, 8], delta 5e-3): 4.5e-9 at the
    # production m = 5 and 5.0e-13 at m = 10.  Both are far under the 0.038 mag
    # scatter the fit works against, which is what the panel count is chosen
    # by; the 1e-7 mag of Section 5.1 was a tighter requirement than the data
    # can use.  (1.1e-9 and 8e-14 at m = 10 and 20 on the earlier [-12, 6].)
    assert err < 1e-7


def test_circle_anchor_exact():
    """gamma(0) = 0 exactly: zero is a panel edge, not a partial panel."""
    out = run_checks(D.M_PER_SIDE, s_eval=[0.0])
    assert out["circle"][0, 0] == 0.0 and out["circle"][0, 1] == 0.0


@pytest.mark.parametrize("m_per_side", [D.M_PER_SIDE, 2 * D.M_PER_SIDE])
def test_unit_rate(m_per_side):
    """ds/dp = 1 gives s = p, through the v grid and its Jacobian."""
    p = np.concatenate([np.linspace(-20, 45, 27), [0.0, -15.0, 40.0]])
    out = run_checks(m_per_side, p_eval=p)
    err = np.abs(out["s_unit"] - p).max()
    print(f"unit rate m={m_per_side}: max error {err:.2e}")
    assert err < 1e-12
    assert out["s_unit"][np.where(p == 0.0)[0][0]] == 0.0


@pytest.mark.parametrize("m_per_side", [D.M_PER_SIDE, 2 * D.M_PER_SIDE])
def test_smooth_integrand(m_per_side):
    """int_0^x exp(0.3t) cos(3t) dt against its antiderivative on the v grid."""
    vlo, vhi = np.arcsinh(np.array(D.MODEL_WINDOW) / D.DELTA_P)
    x = np.linspace(vlo, vhi, 33)
    out = run_checks(m_per_side, x_eval=x)
    F = lambda t: np.exp(0.3 * t) * (0.3 * np.cos(3 * t) + 3 * np.sin(3 * t)) / (0.09 + 9)
    err = np.abs(out["I_test"] - (F(x) - F(0.0))).max()
    print(f"smooth m={m_per_side}: max error {err:.2e}")
    assert err < 1e-13


def test_clamp_and_taper():
    """Section 4: the clamp is inert inside and saturates outside; the taper is
    unity at the data window edge, half a percent down at the modelling window
    edge, e^-1 at eps_t outside, gone at 2 eps_t; the softplus never overflows."""
    margin = 5.0 / (D.MODEL_WINDOW[1] - D.MODEL_WINDOW[0])   # the data window edge
    w = np.array([0.5, 0.2, margin, 1 - margin, 0.0, 1.0, -D.EPS_T, 1 + D.EPS_T,
                  -2 * D.EPS_T, -1.0, 2.0, -100.0, 100.0])
    out = run_checks(D.M_PER_SIDE, w_eval=w)
    c, t, sp = out["clamped"], out["taper"], out["sp"]
    # inert inside: agreement to exponential accuracy, eps e^{-w/eps}
    assert abs(c[0] - 0.5) < 1e-15 and abs(c[1] - 0.2) < D.EPS * np.exp(-0.2 / D.EPS) * 1.01
    # the edge itself is moved inward by eps ln 2 of the range
    assert abs(c[4] - D.EPS * np.log(2)) < 1e-12 and abs(c[5] - (1 - D.EPS * np.log(2))) < 1e-12
    # saturation outside
    assert abs(c[9]) < 1e-12 and abs(c[10] - 1) < 1e-12
    assert abs(c[11]) < 1e-15 and abs(c[12] - 1) < 1e-15
    # taper
    assert t[2] > 0.99 and t[3] > 0.99                      # data window edge: not touched
    assert abs(t[4] - np.exp(-(0.1 * np.log(2)) ** 2)) < 1e-12  # modelling edge: half a percent
    assert abs(t[6] - np.exp(-1)) < 1e-3 and abs(t[7] - np.exp(-1)) < 1e-3
    assert t[8] < np.exp(-4) + 1e-6                        # gone at 2 eps_t
    assert t[11] == 0.0 and t[12] == 0.0                   # far outside: underflow, not NaN
    # softplus: finite and right at x/eps = +-1e4
    assert np.isfinite(sp).all()
    assert abs(sp[12] - 100.0) < 1e-12 and sp[11] == 0.0
