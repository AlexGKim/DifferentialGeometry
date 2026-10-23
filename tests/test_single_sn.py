"""The fitter of stage 1, through ``single_sn.stan``: the dilation and
translation tests of Section 5.2 as gradients that must vanish at the exact
parameters and fits that must return there, and the gradient against finite
differences of the log density itself.

The dilation test is exact only for a reference that walks at a constant
rate, ``q_0 = 0``: with ``q_0`` a function of phase, dilating time rescales
its argument and the copy leaves the family the model spans at ``x_p = 0``.
So it is run with ``b = 0``.  The translation test uses the full SALT2 fit.
"""

import numpy as np
import pytest

from dgsn.salt.globals import fit_globals
from dgsn.stan import data as D
from dgsn.stan.run import fixed_param, log_prob, optimize

NAMES = ["m_peak", "colr", "x_s", "x_p", "x_psi", "delta_t"]


@pytest.fixture(scope="module")
def salt_globals():
    return fit_globals(k_psi=D.K_PSI, k_q=D.K_Q, x1_ref=0.0)


def synth_epochs(G, params, dilation=1.0, shift=(0.0, 0.0), cadence=1.5, snr=50.0,
                 p_range=(-14.5, 39.5), eps=D.EPS):
    """Noise-free epochs of the model at ``params``, observed at times dilated
    by ``dilation`` and translated by ``shift`` magnitudes per band."""
    p = np.arange(p_range[0], p_range[1], cadence)
    dt_rest = np.concatenate([p, p + 0.5 * cadence])
    band = np.concatenate([np.ones(len(p), int), 2 * np.ones(len(p), int)])
    data = D.grid_data()
    data["eps"] = eps
    data.update({k: G[k] for k in ("a", "b", "c", "d")})
    data.update(params)
    data.update(z=0.0)
    data.update(D.eval_points(dt=dt_rest, band=band))
    f = fixed_param("forward", data)["f_hat"]
    f = f * 10 ** (-0.4 * np.where(band == 1, shift[0], shift[1]))
    return dict(dt=dt_rest * dilation, f=f, sig_f=np.full(len(f), f.max() / snr), band=band, z=0.0)


def fit_data(G, epochs, prior_on=False, eps=D.EPS):
    data = D.grid_data()
    data["eps"] = eps
    data.update({k: G[k] for k in ("a", "b", "c", "d")})
    data.update(D.epochs_data(epochs["dt"], epochs["f"], epochs["sig_f"], epochs["band"],
                              epochs["z"]))
    data.update(D.prior_data(prior_on=prior_on, mu_p=G["x_p"], mu_psi=G["x_psi"],
                             xs_bar=G["x_s"], m_bright=18.0))
    return data


def perturbed(truth):
    return {k: truth[k] + dv for k, dv in
            zip(NAMES, [0.2, 0.1, 0.3, 0.5, 0.5, 0.4])}


def check_exact_fit(data, truth, tol=1e-6):
    """The sharp check is the gradient at the exact parameters; the fit from
    a perturbed start must then return there.  x_psi is held to 1e-3 rather
    than 1e-6: the light curve is nearly flat along it (the departure mode is
    concentrated at the turn, where little arclength is), so the optimizer
    resolves it only to where -2 ln L stops changing --- which is what the
    last assertion measures directly."""
    lp_truth, grad = log_prob("single_sn", data, truth)
    print(f"gradient at the exact parameters: {np.abs(grad).max():.2e}")
    assert np.abs(grad).max() < 1e-6
    fit = optimize("single_sn", data, perturbed(truth))
    err = {k: fit[k] - truth[k] for k in NAMES}
    print("fit - truth:", {k: f"{v:.1e}" for k, v in err.items()})
    # 1e-3 is the optimizer's reach, not the model's: the gradient above and
    # the -2 ln L below are the precision statements.
    assert max(abs(v) for v in err.values()) < 1e-3
    lp_fit, _ = log_prob("single_sn", data, {k: fit[k] for k in NAMES})
    print(f"-2 ln L at the fit above its value at the truth: {2 * (lp_truth - lp_fit):.1e}")
    assert 2 * (lp_truth - lp_fit) < 1e-6


def test_dilation(salt_globals):
    """A time-dilated copy of a constant-rate reference fits with
    x_p = x_psi = 0 and x_s = -log(dilation), the rest unchanged.

    The copy's epochs must stay inside the modelling window, and the phase
    clamp's exponential tail --- eps e^{-w/eps} of the window, 1e-4 d at the
    data-window edge with eps = 0.01 --- displaces the original and the copy
    by different amounts, so the clamp width is set to 1e-3 here, where the
    tail is below rounding.  The clamp has its own tests."""
    G = dict(salt_globals, b=np.zeros_like(salt_globals["b"]))
    lam = 1.1
    truth0 = D.params_data(m_peak=18.0, colr=0.05, x_s=-2.0, x_p=0.0, x_psi=0.0, delta_t=0.0)
    epochs = synth_epochs(G, truth0, dilation=lam, p_range=(-14.0, 39.0), eps=1e-3)
    truth = dict(truth0, x_s=-2.0 - np.log(lam))
    check_exact_fit(fit_data(G, epochs, eps=1e-3), truth)


def test_translation(salt_globals):
    """A constant vector added to both bands moves m_peak by its g component
    and colr by its g - r component, and nothing else."""
    G = salt_globals
    truth0 = D.params_data(m_peak=18.0, colr=0.05, x_s=G["x_s"], x_p=G["x_p"], x_psi=G["x_psi"])
    dg, dr = 0.3, -0.2
    epochs = synth_epochs(G, truth0, shift=(dg, dr))
    truth = dict(truth0, m_peak=18.0 + dg, colr=0.05 + (dg - dr))
    check_exact_fit(fit_data(G, epochs), truth)


def test_gradient_vs_finite_differences(salt_globals):
    """Stan's gradient against central differences of its own log density, at
    a point away from the optimum, priors and scatter on."""
    G = salt_globals
    truth0 = D.params_data(m_peak=18.0, colr=0.05, x_s=G["x_s"], x_p=G["x_p"], x_psi=G["x_psi"])
    epochs = synth_epochs(G, truth0)
    data = fit_data(G, epochs, prior_on=True)
    point = perturbed(truth0)
    _, grad = log_prob("single_sn", data, point)
    for h in (1e-4, 1e-5):
        fd = np.empty(len(NAMES))
        for i, k in enumerate(NAMES):
            up, dn = dict(point), dict(point)
            up[k] += h
            dn[k] -= h
            fd[i] = (log_prob("single_sn", data, up)[0] - log_prob("single_sn", data, dn)[0]) / (2 * h)
        rel = np.abs(grad - fd) / np.maximum(np.abs(fd), 1e-3)
        print(f"h={h:.0e}: max relative disagreement {rel.max():.2e}")
        assert rel.max() < 1e-5
