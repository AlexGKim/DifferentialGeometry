"""The hierarchy program against stage 2, its gauge, its sign relabeling,
its reduce_sum, and the split."""

import numpy as np
import pytest

from dgsn.data.sample import (N_VALID_PER_FIFTH, analysis_sample, median_x1, stage3_sample,
                              training_split, x1_fifth)
from dgsn.salt.globals import fit_globals
from dgsn.stan import data as D
from dgsn.stan.hier_data import build, inits
from dgsn.stan.run import log_prob, model

K_PSI, K_Q = 12, 8


@pytest.fixture(scope="module")
def sample():
    return analysis_sample()


@pytest.fixture(scope="module")
def names(sample):
    return stage3_sample()[:3]


@pytest.fixture(scope="module")
def G():
    return fit_globals(K_PSI, K_Q, median_x1())


def test_split_is_the_note_split(sample):
    train, valid = training_split()
    assert len(train) == 442 and len(valid) == 110
    assert not set(train) & set(valid) and set(train) | set(valid) == set(sample.index)
    per_fifth = x1_fifth(sample.loc[valid]).value_counts().sort_index().tolist()
    assert per_fifth == [N_VALID_PER_FIFTH] * 5
    assert training_split() == (train, valid)          # deterministic


def test_hier_matches_single_sn_on_one_object(sample, names, G):
    """With the globals pinned at the SALT2 values and the hyperparameters at
    stage 2's, the per-object term of hier.stan is single_sn.stan's: the two
    log densities differ by the priors on z, alpha, and the pinned constants,
    which are evaluated here and removed."""
    from dgsn.data.epochs import likelihood_epochs

    name = names[0]
    row = sample.loc[name]
    ep = likelihood_epochs(name, row.t0, row.redshift)
    sigma = 0.05
    ini = inits([name], K_PSI, K_Q, median_x1(), sample=sample, globals_=G)
    ini.update(sigma_g=sigma, sigma_r=sigma, mu_delta=0.0, s_delta=1.0)
    point = {k: ini[k] for k in ini}
    # single_sn at the same per-SN point and priors
    d2 = D.grid_data(k_psi=K_PSI, k_q=K_Q)
    d2.update({k: G[k] for k in "abcd"})
    d2.update(D.epochs_data(ep["dt"], ep["f"], ep["sig_f"], ep["band"], ep["z"], sigma, sigma))
    d2.update(D.prior_data(prior_on=True, mu_p=ini["mu_p"], mu_psi=ini["mu_psi"], mu_delta=0.0,
                           s_delta=1.0, xs_bar=D.XS_BAR, m_bright=ep["m_bright"]))
    p2 = dict(m_peak=float(ini["m_peak"][0]), colr=float(ini["colr"][0]),
              x_s=float(ini["x_s"][0]), x_p=ini["mu_p"] + float(ini["u_p"][0]),
              x_psi=ini["mu_psi"] + float(ini["u_psi"][0]),
              delta_t=ini["mu_delta"] + ini["s_delta"] * float(ini["u_del"][0]))
    lp2, _ = log_prob("single_sn", d2, p2)
    dh = build([name], K_PSI, K_Q, sample=sample)
    lph, _ = log_prob("hier", dh, point)
    # hier's extras: std_normal on z_a, z_c (up to constants Stan drops: -z^2/2),
    # normal(0, tau k^-alpha) on the centred b_til, d_til (kept terms
    # -log w - b^2/2w^2), exponential(1) on alpha - 1/2; the per-SN priors are
    # the same in both.  The gauge-fixed a..d in hier equal G's, so the
    # likelihoods agree.
    extra = -0.5 * sum(float(np.sum(np.square(ini[k]))) for k in ("z_a", "z_c"))
    for k, tau, K in (("b_til", ini["tau_q"], K_Q), ("d_til", ini["tau_psi"], K_PSI)):
        kk = np.arange(1, K)
        w = tau * (kk / np.exp(np.log(kk).mean())) ** -2.0
        extra += float(np.sum(-np.log(w) - 0.5 * (ini[k] / w) ** 2))
    extra += -(ini["alpha_q_m"] + ini["alpha_psi_m"])
    print(f"single_sn lp {lp2:.6f}, hier lp - extras {lph - extra:.6f}")
    assert abs((lph - extra) - lp2) < 1e-6


def test_gauge_and_signs(sample, names, G):
    """From a fixed_param draw at random z: the gauge identities hold, and the
    sign flags relabel a, c and every x_p, x_psi jointly, leaving lp unchanged."""
    dh = build(names, K_PSI, K_Q, sample=sample)
    ini = inits(names, K_PSI, K_Q, median_x1(), sample=sample, globals_=G)
    rng = np.random.default_rng(3)
    for k in ("z_a", "b_til", "z_c", "d_til", "u_p", "u_psi"):
        ini[k] = rng.normal(size=len(ini[k]))
    fit = model("hier").sample(data=dh, inits=ini, fixed_param=True, chains=1, iter_warmup=0,
                               iter_sampling=1, adapt_engaged=False, sig_figs=18,
                               show_progress=False)
    v = {k: np.asarray(fit.stan_variable(k))[0] for k in ("a", "b", "c", "d", "flip_q", "flip_psi",
                                                          "a_rel", "c_rel", "x_p_rel", "x_psi_rel")}
    # G_s from forward.stan at the same grid
    fd = D.grid_data(k_psi=K_PSI, k_q=K_Q); fd.update(D.zero_globals(K_PSI, K_Q)); fd.update(D.params_data())
    fd.update(z=0.0); fd.update(D.eval_points())
    from dgsn.stan.run import fixed_param
    Gs = fixed_param("forward", fd)["Gs_out"]
    assert abs(v["a"] @ v["a"] - 1) < 1e-14 and abs(v["b"] @ v["a"]) < 1e-14
    assert abs(v["c"] @ Gs @ v["c"] - 1) < 1e-13 and abs(v["d"] @ Gs @ v["c"]) < 1e-13
    sq, sp = (-1 if v["flip_q"] else 1), (-1 if v["flip_psi"] else 1)
    x_p = ini["mu_p"] + ini["u_p"]
    x_psi = ini["mu_psi"] + ini["u_psi"]
    assert np.array_equal(v["a_rel"], sq * v["a"]) and np.array_equal(v["c_rel"], sp * v["c"])
    assert np.allclose(v["x_p_rel"], sq * x_p, atol=0) and np.allclose(v["x_psi_rel"], sp * x_psi, atol=0)
    # lp invariant under the joint relabeling on the raw scale: negate z_a with
    # (mu_p, u_p), z_c with (mu_psi, u_psi)
    lp0, _ = log_prob("hier", dh, ini)
    flipped = dict(ini, z_a=-ini["z_a"], u_p=-ini["u_p"], mu_p=-ini["mu_p"],
                   z_c=-ini["z_c"], u_psi=-ini["u_psi"], mu_psi=-ini["mu_psi"])
    # b, d are invariant under the flips (b = b_til - (b_til.a)a is even in a)
    lp1, _ = log_prob("hier", dh, flipped)
    print(f"flip flags {v['flip_q']:.0f} {v['flip_psi']:.0f}; lp {lp0:.6f} vs relabeled {lp1:.6f}")
    assert abs(lp0 - lp1) < 1e-8


def test_reduce_sum_grainsize(sample, names, G):
    ini = inits(names, K_PSI, K_Q, median_x1(), sample=sample, globals_=G)
    lps = []
    for gs in (1, len(names)):
        dh = build(names, K_PSI, K_Q, grainsize=gs, sample=sample)
        lps.append(log_prob("hier", dh, ini)[0])
    assert abs(lps[0] - lps[1]) < 1e-10 * max(1.0, abs(lps[0]))
