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

**Depends on:** `ztfcosmo` (or the local archive). No other module imports it.

The quality cuts and the `z < 0.05` selection live here and nowhere else, so
that changing selection is a one-file change. Three things this module must get
right, all verified against the archive:

- **Flag cut is `flag & 31 == 0`**, not `flag == 0`. The flag is a bitmask;
  only bits `[1,2,4,8,16]` indicate bad photometry. Cutting on `flag == 0`
  preferentially keeps faint baseline epochs and destroys the sample.
- **Milky Way extinction** (`mwebv`) is known and removed deterministically
  before fitting. It is not a free parameter.
- **Missing data**: 16 SNe passing metadata cuts have no light-curve file, and
  36 catalogue rows have NaN `t0`. Both are dropped here, with counts logged.

Band selection is a config field, since the same loader serves both the 599-SN
`g,r` primary sample and the 177-SN `g,r,i` torsion subsample.

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
  This is what magnitude space buys for the *rigid* nuisances — distance,
  absolute luminosity, peculiar velocity, per-band calibration. It does **not**
  cover dust; see the next test for why that matters.
- **First variation of curvature (the complement, with an exact expected value).**
  `model-definition.tex` now states the perturbation formula rather than describing
  the effect in words: for `delta_gamma = phi*T + psi*N`,

  ```
  Delta_kappa = psi'' + kappa^2 * psi + kappa' * phi
  ```

  so applying `gamma(s) -> gamma(s) + c*u(s)` must change `kappa` by the value this
  predicts, not merely change it. Assert agreement with the closed form to solver
  tolerance — this is a sharper test than "must differ", and it is the numerical
  check on a formula the note now leans on in three places (the rigid/deformation
  split, the `ubar`-drops-out claim, and the `Delta kappa` diagnostic).

  Two corollaries to assert alongside it, both cheap:
  - **Constant `u` gives exactly zero.** The formula annihilates translations
    identically, which is the whole content of the magnitude-space argument. This
    replaces the separate isometry test's role for dust.
  - **Only `du` contributes.** For `u = ubar + du`, feeding `ubar` alone leaves
    `kappa` unchanged and the entire change comes from `du`. A sign or normalisation
    error in that split would silently inflate every downstream dust diagnostic.

  Use a **deliberately exaggerated** `u(s)` for the amplitude, and assert
  non-constancy of `u` in the fixture. A *realistic* `u` puts `Delta kappa` near
  solver tolerance — the real phase variation is second order — so the test would be
  measuring the ODE integrator rather than the geometry. Decouple the mechanism from
  the physical amplitude; the amplitude question belongs to the `Delta kappa`
  diagnostic below.
- **Reparameterisation invariance:** changing `tscale, a1, a2` leaves the set
  of points `{gamma(s)}` unchanged. This is the shape/timing separation claim.
- **Oracle agreement:** `frenet.integrate` and `direct.curve_and_invariants`
  agree on `kappa` to tolerance. Regression test against silent solver bugs.
- `ds/dt > 0` is enforced across the fit window for sampled `(a1, a2)`.
- **Gauge anchor:** at `s = 0` the reconstructed curve satisfies `dm_g/ds = 0`
  and `T(0) = ±(0,1)`. Two things this test must **not** assert. The sign is
  empirical, not gauge — it records which band peaks first, is expected to be
  `(0,-1)`, and is carried in config; assert against the configured value, never
  a hard-coded `+1`. And "one condition fixes both gauges" holds only for `n=2`:
  in `R^3` the frame is an element of `SO(3)` and two conditions remain unfixed,
  so the `n=3` case must be skipped rather than asserted.
- **Regularity — coincident maxima break arclength:** construct a synthetic curve
  whose bands peak at the same epoch and assert `ds/dt -> 0` there, i.e. that
  unit-speed parameterisation fails rather than silently returning garbage.
  Non-coincident band maxima is a condition on the data, and the code should
  detect its violation.
- **Total turning:** for a synthetic hairpin curve with known asymptotes along
  `(1,1)`, `integrate(kappa) ds` recovers `pi` to tolerance. Validates the sum
  rule machinery on a case where the answer is known.

**Model.** Round-trip recovery: generate synthetic photometry from known
latents with realistic ZTF cadence and noise, confirm training recovers them
within uncertainty. This is the only use of simulation, and it tests the
machinery, not the science.

**Data.** Masking correctness — padded epochs contribute exactly zero to
`chi2`. Cuts are applied idempotently.

**Evaluation.** Leave-one-out is genuinely leaving one out: assert the held-out
SN's data never enters the training loss.

**Identification (rung L2c).** Two diagnostics that are part of the deliverable,
not optional extras:

- Population `Corr(c, theta)` under the identifying assumption that line-of-sight
  dust is independent of explosion physics, computed **both with and without** the
  condition imposed. These are *reported numbers*, not pass/fail assertions — the
  test asserts only that both are computed and recorded. The unimposed value is the
  interesting one: small means the data determine the split and the assumption is
  merely a check; large means it is load-bearing.
- Error-weighted overlap `rho^2 = (d' W P d) / (d' W d)` with
  `P = J (J' W J)^-1 J' W`, `J` the Jacobian columns `df/dtheta_1, df/dtheta_2`,
  `W = diag(1/sigma_i^2)`, over the epochs actually observed. Two details the
  implementation must get right, both now stated in `model-definition.tex`: build it
  in **flux**, where the errors are Gaussian, via
  `df/dtheta = -0.4 ln10 * f * dm/dtheta`; and form `d` from `du`, the varying part,
  **not** the full `u` — the mean part is a translation that cannot be confused with
  shape, so including it inflates `rho^2`. Unit-test the overlap function on
  constructed cases where the answer is known (orthogonal ⇒ 0, parallel ⇒ 1), since
  the astrophysical value has no ground truth to check against.
- `Delta kappa` from `du(s)` versus the `kappa` range spanned over fitted `theta`,
  **stratified by fitted `c`** rather than sample-averaged. This is the direct test
  of the second-order claim, and the stratification is the point: the hierarchy is
  expected to degrade at the reddened end, and an average would hide exactly the
  regime where it fails. Assert the stratification exists — a single scalar here is
  a bug, not a result.

---

## 6. Risks

| Risk | Mitigation |
|---|---|
| Timing warp degenerate with shape | Shrink `a1,a2` toward zero; report correlation with `theta` as a headline diagnostic. A large correlation is a negative result to publish, not to hide. |
| `kappa -> 0` at inflections breaks the frame | Use an inflection-robust frame construction; report where along `s` curvature approaches zero. |
| Dust forced into `theta` | Not preventable by geometry — dust and intrinsic diversity are both deformations. But only the *varying* part `c*du(s)` competes with shape, and that is second order where the signal is first, so the expectation is that this risk is small. Rung L2c plus `Corr(c, theta)` measures it, stratified by `c` since the hierarchy weakens for red objects. Harmless for the Hubble-diagram result either way. |
| Good news on shape obscures the colour degeneracy | Separate risk, opposite sign. `c*ubar` is indistinguishable from a per-SN intrinsic colour offset at *first* order, and no order counting helps. Do not let a clean shape-channel result be reported as though dust were solved. Quote the two channels separately. |
| 7 latents overfit | All scoring is out-of-sample. No in-sample number is ever reported as evidence. |
| Wrong flag cut silently guts the sample | Regression test asserting the `z<0.05` primary sample has 599 SNe at >=5 good `g` and `r` epochs. |
| Torsion subsample too small (177) | Report it as a separate analysis with its own uncertainties; do not pool with the primary sample. |
| `kappa` extrapolates nonsensically outside the window | Report `integral kappa ds` against the sum-rule value `pi`. A trained `kappa` that extrapolates to something far from `pi` is unphysical even if it fits the window well. Diagnostic only — never a hard constraint, since the asymptotes are unobserved. |
| JAX debugging cost | `direct.py` oracle runs eagerly and small; debug there first. |

---

## 7. Sample sizes

Measured from the archive on 2026-08-20 (`z<0.05`, both quality flags, epochs
with `flag & 31 == 0` and SNR > 5 in the phase window):

| min epochs/band | g | r | i | g&r | g&r&i |
|---|---|---|---|---|---|
| >=3 | 636 | 639 | 216 | 636 | 215 |
| >=5 | 606 | 631 | 179 | **599** | **177** |
| >=10 | 490 | 537 | 82 | 465 | 80 |

The `g&r&i` column is the torsion subsample: the only configuration in which
`theta_2` can be tested as a genuine second invariant rather than a second
direction of curvature variation.
