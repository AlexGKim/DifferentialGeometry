# DifferentialGeometry

Empirical modelling of Type Ia supernova multi-band light curves using
differential geometry.

## The idea in one paragraph

A supernova observed in `n` bands supplies, at each epoch, a point in `R^n`.
Over the event those points trace a **curve**. Parameterise it by arclength `s`
and the curve is determined, up to a rigid motion, by `n-1` curvature
invariants `kappa_i(s)`. The model splits a light curve into two independent
pieces: `kappa(s)` carries **shape**, and a time-to-arclength map `s(t)` carries
**timing**. A reparameterisation of time changes how fast a supernova moves
along its curve but does not change the curve. That exact separation is the
scientific claim being tested.

The classical model to beat is SALT2 (Guy et al. 2007), whose `x1` entangles
shape and timing in a single linear coefficient.

## Current configuration

| | |
|---|---|
| Bands | ZTF `g`, `r` (public survey). `i` is proprietary — out of scope. |
| Ambient space | `R^2`, so **one** invariant `kappa(s)`. No torsion. |
| Code generality | Written for general `R^n` so a third band promotes the model without restructuring. |
| Data | ZTF SN Ia DR2 via `ztfidr`, cut to `z < 0.05`. |
| Curve space | **Magnitude**. Likelihood evaluated in **native flux**. |
| Phase window | Rest-frame `[-15, +40]` d (config parameter). |

## Per-SN parameters (baseline, 7)

`mu` (normalization) · `theta_1`, `theta_2` (shape, condition the network) ·
`sigma`, `t_max` (timing) · `a_1`, `a_2` (timing warp)

SALT2 uses four. The excess is deliberate and is the object of study, not a
free choice.

## Non-negotiable design decisions

These were settled deliberately. Do not silently revise them.

1. **Curve in magnitude space, likelihood in flux.** In magnitude space both
   distance (`(1,1)/sqrt2`) and colour (`(1,-1)/sqrt2`) are *translations*,
   hence isometries, hence `kappa` is exactly invariant. In flux space distance
   is an isotropic dilation and extinction an *anisotropic* one — neither is an
   isometry, and the claim that `kappa` is intrinsic collapses. But ZTF errors
   are Gaussian in flux and negative fluxes are meaningful, so the likelihood
   must use native flux. **Never convert data to magnitudes.**

2. **SALT2 is strictly downstream.** Used only as an independent benchmark on
   the same objects. Never to preprocess, K-correct, interpolate, or
   initialise. Using it upstream would make the independence claim false.

3. **`z < 0.05`.** A two-band model has no SED and cannot self-consistently
   K-correct. The cut buys independence at the cost of sample size. Do not
   relax it to gain statistics without revisiting the whole K-correction story.

4. **No per-SN rotation.** A rotation of `R^n` mixes photometric bands, which
   have fixed physical identities. The initial point and frame are global.

5. **Peculiar velocity is not in the likelihood.** `mu` is free per SN and
   absorbs it exactly. It enters only at the Hubble-diagram stage.

6. **Evaluation is out-of-sample, always.** With 7 latents and a conditioned
   network, good in-sample fits are guaranteed and carry no evidential weight.

## The deliverable is a ladder, not a fit

| Rung | Free per SN | Question |
|---|---|---|
| L0 | `mu, sigma, t_max` | is a rigid template enough? |
| L1 | `+ theta_1` | does one shape parameter earn its place? |
| L2 | `+ theta_2` | does a second? |
| L2c | L2 `+ c` translation | does dust need its own direction? |
| L3 | `+ a_1, a_2` | does nonlinear timing earn its place? |

Scored on **held-out Hubble residual scatter** under cross-validation
(leave-one-out jackknife, following SALT2's own validation method).

The ladder is also the **training schedule**: fit L0 to convergence, then
introduce latents one at a time. SALT2 stages its components the same way.

## Stack

JAX · diffrax (Frenet ODE, adjoint) · equinox (networks) · optax ·
numpyro (posteriors, later).

Training is an **auto-decoder**: network weights and the `N x 7` array of
per-SN latents live in one pytree and are optimised jointly.

Two implementations of the geometry, which must agree on `kappa`:
- `geometry/frenet.py` — diffrax ODE integration. **Primary.**
- `geometry/direct.py` — parameterise `gamma(s)` directly, recover `kappa` by
  autodiff. Fast correctness oracle and regression test.

## Layout

```
src/dgsn/
  data/      ztfidr ingest, quality cuts, z<0.05 selection
  geometry/  frenet.py (R^n), direct.py, reparam.py, invariants.py
  model/     kappa_net.py, forward.py, likelihood.py
  train/     auto-decoder loop, staged schedule
  infer/     per-SN fits, numpyro posteriors
  eval/      ladder, cross-validation, Hubble residuals
docs/tex/    LaTeX science docs (model-definition, salt2-distillation)
docs/        Markdown engineering docs (plan.md)
configs/     one YAML per ladder rung
references/  source papers
```

## Documentation convention

- **Scientific content → LaTeX** in `docs/tex/`. Always `\input{macros}` so
  notation stays consistent; cite from `refs.bib` (entries sourced from NASA
  ADS via the `nasa-ads` MCP server).
- **Engineering content → Markdown**. This file, `README.md`, `docs/plan.md`.

## Known limitations to keep in view

- **Two bands admit no torsion.** `theta_2` is a second direction of variation
  in `kappa`, not the torsion originally envisaged. Do not describe it as
  torsion in writing.
- **Timing warp vs shape degeneracy.** Over a finite, noisy window, a warp in
  `s(t)` and a change in `kappa` can mimic each other. Report the fitted
  correlation between `(a_1,a_2)` and `theta` as a diagnostic. A large
  correlation means the shape/timing separation is not being realised — that is
  a negative result to report, not to absorb.
- **Curvature degenerates at inflections.** Where `kappa -> 0` the classical
  Frenet frame is ill-defined; the `(m_g, m_r)` path may approach this near the
  secondary maximum. Use a frame construction robust to inflections and report
  where `kappa` approaches zero.
- **Gauge must be fixed** before latents are interpretable: `s=0` anchored to
  the `g`-band maximum; latents normalised to zero mean, unit variance; sign
  fixed so `theta_1` correlates positively with light-curve width.

## Data access prerequisite

`ztfidr` is only an interface. It reads from `$ZTFIDRPATH` pointing at the
`ztfcosmoidr/dr2` repository, which is **password-protected** ZTF collaboration
data. The package alone is not sufficient — collaboration access is required.
`sample.data` carries SALT2 fits (`x1`, `c`, `t0`) for the downstream benchmark
and supports `get_data(redshift_range=...)` for the `z<0.05` cut.

## Working agreements

- Scientific claims in docs must be traceable to a citation or to code that
  produced the number. No asserted values.
- Report negative results plainly. A ladder rung that fails to earn its
  parameter is a result.
- When a design decision above is genuinely wrong, say so and change it
  explicitly — do not work around it.
