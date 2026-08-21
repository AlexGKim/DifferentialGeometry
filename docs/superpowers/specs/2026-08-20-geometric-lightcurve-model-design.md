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

~~**Pre-explosion (null) epochs are a second epoch set with its own rules**~~ — added
2026-08-21, **removed 2026-08-22** with `CLAUDE.md` decision 8. There is no `is_null`
mask, no null epoch set, and no separate null cut: the likelihood epoch set is the
in-window detections, full stop. What survives from that work is narrower and still
required — `PhotBatch` carries per-band `flux_offset` and `offset_unc`, because
`flux_offset` is subtracted deterministically (the same rule as `1+z` and `mwebv`) and
`offset_unc` is rank-one, which matters for *faint detections late in the window* whether
or not nulls are used. Those column readings remain **not yet verified against
`ztfcosmo`** and that verification is still a prerequisite of writing this module. Null
handling returns only if the parked turn-on rung does.

### `dgsn.geometry`

Pure geometry. Knows nothing about supernovae, flux, or photometry.

```
frenet.integrate(kappa_fn, s_grid, n_dim) -> gamma            # (len(s), n)
direct.curve_and_invariants(gamma_fn, s)  -> (gamma, kappa)
reparam.s_of_t(t, t_max, tscale, a1, a2, z) -> s              # cubic, no singularity
reparam.ds_dt(...)                        -> ds/dt            # for the constraint
```

**Removed 2026-08-22** with `CLAUDE.md` decision 8: the `A, u_expl` arguments to
`s_of_t` and its logarithmic singularity; `reparam.t_expl`; and
`invariants.g_max_arclength`. The last is gone because `s = 0` *is* the `g`-maximum by the
`eq:frameanchor` gauge, so `t_max` maps to `s = 0` and there is no root-find anywhere in
the model.

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
one place `mu` is applied. It is also what absorbs the `m -> +inf` divergence: the
conversion extends continuously to `f = 0`, so a faint prediction is an ordinary value of
the forward model rather than a branch. The divergence itself is never reached — it sits at
`|s| = infinity`, outside the bounded arclength domain.

`likelihood.chi2` is **not diagonal**. Per band,
`C = diag(sigma^2) + offset_unc^2 · 1·1^T`, inverted in closed form by
Sherman–Morrison; `flux_offset` is subtracted from the data, not fitted. Adds no
parameter. Negligible beside bright detections, it matters for the faint epochs late in
the window, where a diagonal treatment overstates their joint information. A diagonal
implementation is a bug, not an approximation — see `CLAUDE.md`.

### `dgsn.train`

Auto-decoder. Network weights and the `(N, n_latent)` latent array are a single
optax-optimised pytree.

```
train(config) -> TrainState
```

Implements the staged schedule: converge rung L0 (`mu, c, w, t_max` — SALT2's
parameter count), then unfreeze latents one at
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
   latents (N,8) ──► reparam.s_of_t ──► s
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
it can be tested later. **Two** windows, not three (revised 2026-08-22, when the
three-window split went with decision 8): the *arclength domain* `[s_min, s_max]`, and the
*selection window* `[-15,+40]` d, which also delimits the likelihood epoch set.

The **arclength domain** is a config value with a checkable condition attached: it must be
wide enough that `reparam.s_of_t` stays inside it across the selection window for every
sampled parameter value. That is a runtime assertion, not an assumption — a fit that walks
off the end of the domain has no forward model, and must fail loudly rather than
extrapolate. There are no `A` or `u_expl` hyperparameters and no `t_expl`; both went with
decision 8, along with the frame-sign entry and the global orientation `Phi`. For `n = 2`
the orientation is **gauge**, fixed by `eq:frameanchor`, so there is nothing to configure
and nothing to fit; for `n >= 3` there are `n(n-1)/2 - 1` global orientation parameters
(two for `n = 3`), fitted, not configured.

The **reddening direction `e_c`** is a unit vector in the `(m_g, m_r)` plane whose status
depends on the rung (revised twice on 2026-08-22). Below L2c it is **exactly
unidentifiable** — free per-SN `mu` along `(1,1)` and `c` along `e_c` span `R^2`, so
changing `e_c` is undone by an invertible relabelling of `(mu, c)` — so at L0–L2 it is a
**config value fixed at `(1,-1)/sqrt2`**, and any attempt to fit it there is a bug. At
**L2c** it becomes a genuine **single sample-wide fit parameter**, the two-band analogue of
SALT2's fitted `beta`, given content by the phase-varying `du(s)`; `(1,-1)/sqrt2` is then
the nested fixed-`R_V` idealisation. Per-SN `c` is the amplitude along `e_c`; only the phase
variation `du(s)` is precomputed from the extinction law and template SED, with `e_c`
setting its mean direction. A per-SN reddening direction is **not** fitted at any rung (it
would reintroduce a forbidden per-SN rotation). So the two per-SN phase-independent shifts
are `mu` along the physics-fixed `(1,1)` and `c` along `e_c`.

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
- **The frame anchor (revised 2026-08-22 — orientation is gauge for `n = 2`).** Assert
  the `eq:frameanchor` gauge holds on the integrated template: `dm_g/ds = 0` at `s = 0`
  with `T(0) = (0,-1)` and `kappa(0) > 0`, and that this is a **minimum of magnitude**
  (`m_g'' = kappa(0) > 0`), i.e. peak brightness — a sign error here silently anchors on the
  faintest point. Assert it costs no generality: for a random turning `kappa`, a vertical
  tangent exists and the anchor is reachable by shifting the network origin. Two things
  this test must **not** do. It must not assert or record a frame sign / peak ordering:
  that prediction is **dropped**, there is no `sign(kappa)` ordering diagnostic, and no
  config entry carries a sign. And it must not try to recover a global orientation `Phi`
  for `n = 2` — there is none.
- **The orientation/origin redundancy (new 2026-08-22).** The sharpest test of the
  counting: integrating `kappa(.;theta)` from initial tangent angle `psi` must give the
  **same set of points** as integrating `kappa(.+a;theta)` from `psi + int_0^a kappa`, to
  solver tolerance, for random `a`. Assert it, then assert the consequence — a
  synthetic round-trip must **fail** to recover an injected global rotation for `n = 2`
  (it is gauge), and must recover exactly `n(n-1)/2 - 1` of them for `n >= 3` (two for
  `n = 3`). A test that recovers a two-band `Phi` is detecting a bug in the gauge fixing.
- **Reddening direction `e_c` (revised twice 2026-08-22 — gauge below L2c).** Below L2c
  assert `e_c` is **exactly unidentifiable**: the Fisher matrix over
  `(mu, c, e_c-angle)` must have a null direction, and two fits with different fixed `e_c`
  (not parallel to `(1,1)`) must reach the *same* likelihood with `(mu, c)` related by the
  predicted invertible relabelling. A round-trip that appears to recover `e_c` at L0–L2 is a
  bug. At **L2c** assert the opposite: inject a synthetic `e_c != (1,-1)/sqrt2`, generate a
  population with a spread of per-SN `c`, and confirm a round-trip recovers it — the
  phase-varying `du(s)` is what breaks the flat direction, so also assert the null direction
  *reappears* if `du(s)` is set to zero. Assert `e_c` is **not** a per-SN latent at any rung
  (fitting it per SN would recover a forbidden per-SN rotation), and that at L2c the imported
  `du(s)` uses `e_c` for its mean direction.
- **Regularity — coincident maxima break arclength:** construct a synthetic curve
  whose bands peak at the same epoch and assert `ds/dt -> 0` there, i.e. that
  unit-speed parameterisation fails rather than silently returning garbage.
  Non-coincident band maxima is a condition on the data, and the code should
  detect its violation.
- **Total turning (integrator machinery only).** For a synthetic hairpin curve with
  known asymptotes along `(1,1)`, `integrate(kappa) ds` recovers `pi` to tolerance.
  This validates the integrator on a case where the answer is known — it is **not** a
  test of a model sum rule, which is now qualitative (`CLAUDE.md` decision 5). The
  branch-resolved turning-split test is **deleted** along with the branch table and the
  peak-ordering prediction (revised 2026-08-22).
- **Vanishing curvature, incidental only (revised 2026-08-22).** There is no by-design
  zero-curvature set any more, so this tests the *incidental* case. For `n = 2`, assert
  `frenet.py` integrates cleanly through an isolated `kappa = 0` — the rotation ODE is
  regular there — and that **`direct.py` returns exactly `0`**, not "small", from what is a
  `0/0` form in its `gamma''` recovery; then assert the two implementations still agree
  across the crossing. For `n >= 3`, assert the code **reports** `tau` as unidentifiable
  wherever `kappa_1 = 0`, not merely noisy: there is no osculating plane, so whatever the
  network emits is meaningless. That is a reporting assertion, not a numerical one, and it
  requires the frame be carried by parallel transport (**Bishop**), which is itself worth an
  assertion.
- **Arclength domain coverage (new 2026-08-22).** Assert that `reparam.s_of_t` maps the
  whole selection window inside the configured `[s_min, s_max]` for every sampled
  parameter value, and that a deliberately narrow domain raises rather than silently
  extrapolating. This replaces the deleted continuity-at-`t_expl` test as the guard on the
  traversal.
- **`s`-origin / `mu` null direction.** Assert that with `s = 0` anchored at the
  `g`-maximum (`eq:frameanchor`) the Fisher matrix has no null direction in
  `(s_offset, mu, c)`, and that removing the anchor produces one. This is an exactly flat
  direction, so it is cheaply detectable and catastrophic if reintroduced. Note the anchor
  must be an **intrinsic** feature of the curve for this to work — the old
  "end of the zero-curvature segment" anchor did not remove the flat direction, because
  nothing forbade the curvature staying zero past `s = 0`.

**Deleted 2026-08-22 with decision 8**, and not to be reinstated while the turn-on rung
stays parked: the straight-segment oracle test, the power-law rise test
(`alpha_X = 0.4 ln10 e_X A`, the fireball case `A = 3.071`), the continuity-at-`t_expl`
test, and the `s_g(theta)` root-find tests (uniqueness and differentiability through the
implicit function theorem). Six tests in, six tests out.

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
| `kappa -> 0` breaks the frame | **Incidental only** since decision 8 was struck (2026-08-22) — there is no by-design zero-curvature set. Harmless for `n=2`: `eq:frenet2d` is a rotation ODE regular at `kappa = 0` and `N` is defined as the fixed 90-degree rotation of `T`. For `n >= 3` the osculating plane does not exist, so a **Bishop / parallel-transport** frame is required — not merely an "inflection-robust" one — and `tau` must be reported unidentifiable there. Either way: report where along `s` curvature approaches zero, near the secondary maximum in particular. |
| Dust forced into `theta` | Not preventable by geometry — dust and intrinsic diversity are both deformations. But only the *varying* part `c*du(s)` competes with shape, and that is second order where the signal is first, so the expectation is that this risk is small. Rung L2c plus `Corr(c, theta)` measures it, stratified by `c` since the hierarchy weakens for red objects. Harmless for the Hubble-diagram result either way. |
| Good news on shape obscures the colour degeneracy | Separate risk, opposite sign. `c*ubar` is indistinguishable from a per-SN intrinsic colour offset at *first* order, and no order counting helps. Do not let a clean shape-channel result be reported as though dust were solved. Quote the two channels separately. |
| 8 per-SN parameters overfit | All scoring is out-of-sample. No in-sample number is ever reported as evidence. Note L0 carries only 4, matching SALT2's count, so the base rung is not over-parameterised at all. |
| Wrong flag cut silently guts the sample | Regression test asserting the `z<0.05` primary sample has 599 SNe at >=5 good `g` and `r` epochs. |
| Torsion subsample too small (177) | Report it as a separate analysis with its own uncertainties; do not pool with the primary sample. |
| `kappa` extrapolates nonsensically outside the window | Report `integral kappa ds` against the qualitative value `≈ pi`. A trained `kappa` extrapolating to something far from `pi` is suspect, but the total is not pinned exactly — **both** asymptotes lie outside the bounded arclength domain and neither is observed (`CLAUDE.md` decision 5). Diagnostic only, never a hard constraint. |
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
