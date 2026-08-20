# Design: Differential-Geometric SN Ia Light-Curve Model

**Date:** 2026-08-20
**Status:** Approved design, implementation not started

The scientific content of this model — the geometry, the likelihood, the
degeneracy analysis, the validation argument — lives in
[`docs/tex/model-definition.tex`](../../tex/model-definition.tex). This
document is the **engineering** design: module boundaries, interfaces, data
flow, testing strategy, and risks. It does not restate the science.

---

## 1. Scope

Build a trainable implementation of the geometric light-curve model and the
five-rung model ladder, evaluated on held-out Hubble residuals against a SALT2
benchmark, using ZTF SN Ia DR2 photometry at `z < 0.05`.

**Out of scope for this iteration:** cosmological parameter inference, the
third band (proprietary), spectroscopic data, host-galaxy correlations.

---

## 2. Module boundaries

Each module has one purpose, a stated interface, and is testable alone.

### `dgsn.data`

Converts ZTF SN Ia DR2 into padded, masked arrays ready for `vmap`.

```
load_sample(config)      -> SampleMeta   # per-SN: id, z, salt2 x1/c/t0
load_photometry(config)  -> PhotBatch
```

`PhotBatch` is a frozen dataclass of JAX arrays, all shape `(N, K)` where `K`
is the padded max epochs per SN:

| field | meaning |
|---|---|
| `time` | observer-frame MJD |
| `flux`, `flux_err` | native ZTF difference-imaging flux; negatives retained |
| `zp` | per-epoch zeropoint |
| `band_idx` | integer index into the band list |
| `mask` | `True` for real epochs, `False` for padding |

**Depends on:** `ztfidr`, `$ZTFIDRPATH`. No other module imports `ztfidr`.

The quality cuts and the `z < 0.05` selection live here and nowhere else, so
that changing selection is a one-file change.

### `dgsn.geometry`

Pure geometry. Knows nothing about supernovae, flux, or photometry.

```
frenet.integrate(kappa_fn, s_grid, n_dim) -> gamma            # (len(s), n)
direct.curve_and_invariants(gamma_fn, s)  -> (gamma, kappa)
reparam.s_of_t(t, t_max, tscale, a1, a2)  -> s
reparam.ds_dt(...)                        -> ds/dt            # for the constraint
invariants.gauge_anchor(gamma, s)         -> s_offset
```

`frenet.integrate` solves the generalized Frenet system in `R^n` via diffrax
with an adjoint. `direct` is the independent oracle: parameterise `gamma(s)`
with a network and recover `kappa` by `jax.jacfwd`. The two must agree.

### `dgsn.model`

```
kappa_net.KappaNet(n_dim, width, depth)   # equinox Module
  __call__(s, theta) -> kappa             # (n_dim - 1,)

forward.predict_flux(params, latents, batch) -> flux_pred
likelihood.chi2(params, latents, batch)      -> scalar
```

`forward` composes `reparam` → `frenet` → magnitude → flux → mask. The
magnitude-to-flux conversion is the *only* place the two spaces meet, and the
one place `mu` is applied.

### `dgsn.train`

Auto-decoder. Network weights and the `(N, n_latent)` latent array are a single
optax-optimised pytree.

```
train(config) -> TrainState
```

Implements the staged schedule: converge rung L0, then unfreeze latents one at
a time. Rung membership is a config field, not a code branch.

### `dgsn.infer`

Per-SN fits with the network frozen; numpyro posteriors over latents. Used for
the "does this parameter earn its place" question, which is about posterior
widths, not point estimates.

### `dgsn.eval`

```
ladder.run(configs)          -> LadderResult
crossval.leave_one_out(...)  -> residuals
hubble.residual_scatter(...) -> scatter, alpha, beta
```

Leave-one-out follows SALT2's jackknife: retrain excluding one SN, refit only
that SN's own latents, examine its residuals.

---

## 3. Data flow

```
ztfidr ──► dgsn.data ──► PhotBatch (padded, masked)
                              │
                              ▼
   latents (N,7) ──► reparam.s_of_t ──► s
                                        │
   KappaNet(s, theta) ──► frenet.integrate ──► gamma(s)  [magnitudes]
                                        │
                              + mu, then 10^(-0.4(m - zp))
                                        │
                                        ▼
                              chi2 against native flux
                                        │
                              optax ──► weights + latents
```

Two properties this enforces: data are never converted to magnitudes, and SALT2
never appears anywhere in this path.

---

## 4. Configuration

One YAML per ladder rung in `configs/`, differing only in which latents are
free and which regularisers are active. Shared keys: phase window, band list,
redshift cut, network width/depth, solver tolerances, optimiser schedule.

The phase window `[-15, +40]` d is a config value precisely so sensitivity to
it can be tested later.

---

## 5. Testing strategy

**Geometry (exact, no data).** These are the tests that matter most, because
they check the mathematics rather than the fit quality.

- A circle of radius `R` integrates to constant `kappa = 1/R`; recovered
  numerically to solver tolerance.
- A straight line gives `kappa = 0`.
- **Isometry invariance:** translating `gamma` leaves `kappa` bitwise-close.
  This is the central claim of the magnitude-space choice and must be a test.
- **Reparameterisation invariance:** changing `tscale, a1, a2` leaves the set
  of points `{gamma(s)}` unchanged. This is the shape/timing separation claim.
- **Oracle agreement:** `frenet.integrate` and `direct.curve_and_invariants`
  agree on `kappa` to tolerance. Regression test against silent solver bugs.
- `ds/dt > 0` is enforced across the fit window for sampled `(a1, a2)`.

**Model.** Round-trip recovery: generate synthetic photometry from known
latents with realistic ZTF cadence and noise, confirm training recovers them
within uncertainty. This is the only use of simulation, and it tests the
machinery, not the science.

**Data.** Masking correctness — padded epochs contribute exactly zero to
`chi2`. Cuts are applied idempotently.

**Evaluation.** Leave-one-out is genuinely leaving one out: assert the held-out
SN's data never enters the training loss.

---

## 6. Risks

| Risk | Mitigation |
|---|---|
| Timing warp degenerate with shape | Shrink `a1,a2` toward zero; report correlation with `theta` as a headline diagnostic. A large correlation is a negative result to publish, not to hide. |
| `kappa -> 0` at inflections breaks the frame | Use an inflection-robust frame construction; report where along `s` curvature approaches zero. |
| Dust forced into `theta` | Rung L2c exists precisely to measure this. |
| 7 latents overfit | All scoring is out-of-sample. No in-sample number is ever reported as evidence. |
| `$ZTFIDRPATH` access unavailable | Blocks all real-data work. Confirm access before implementation begins. |
| JAX debugging cost | `direct.py` oracle runs eagerly and small; debug there first. |

---

## 7. Open item

**ZTF collaboration access.** `ztfidr` is an interface to a password-protected
repository. Everything downstream of `dgsn.data` can be built and tested
against synthetic photometry, but no scientific result is possible without it.
This should be resolved before implementation starts.
