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
| Data | ZTF SN Ia DR2 (`ztfcosmo`), `z < 0.05`, both quality flags. |
| Primary sample | **599 SNe** with >=5 good `g` and `r` epochs in the phase window. Ambient space `R^2`, **one** invariant `kappa(s)`, no torsion. |
| Torsion subsample | **177 SNe** that additionally have >=5 good `i` epochs. Ambient space `R^3`, so `kappa(s)` **and** `tau(s)`. |
| Code generality | Written for general `R^n`. This is not speculative — the `R^3` subsample exists and is the reason. |
| Curve space | **Magnitude**. Likelihood evaluated in **native flux**. |
| Phase window | Rest-frame `[-15, +40]` d (config parameter). |

Coverage measured directly from the archive (2026-08-20), at `z<0.05` with
`lccoverage_flag` and `fitquality_flag` set, counting epochs with
`flag & 31 == 0` and SNR > 5 inside the phase window:

| min epochs/band | g | r | i | g&r | g&r&i |
|---|---|---|---|---|---|
| >=3 | 636 | 639 | 216 | 636 | 215 |
| >=5 | 606 | 631 | 179 | 599 | 177 |
| >=10 | 490 | 537 | 82 | 465 | 80 |

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

- **Two bands admit no torsion.** In the primary `g,r` analysis `theta_2` is a
  second direction of variation in `kappa`, *not* torsion. Do not call it
  torsion in writing. Torsion is only meaningful on the 177-SN `g,r,i`
  subsample, where it is a genuine second invariant.
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

## Data details

A local copy of the archive is at `data/ztfsniadr2_lite.zip` (1.4 GB, gitignored
along with everything else under `data/`). The maintained interface is
[`ztfcosmo`](https://github.com/ZwickyTransientFacility/ztfcosmo), which needs
no credentials — it reads remotely from `ztfcosmo.in2p3.fr`, or from
`$ZTFCOSMODIR` if pointed at a local copy.

**`tables/snia_data.csv`** — 3628 rows. Columns include `ztfname`, `redshift`,
SALT2 fits `t0, x0, x1, c` with full covariances, `mwebv`, `sn_type`,
`lccoverage_flag`, `fitquality_flag`.

**`lightcurves/<name>_lc.csv`** — whitespace-delimited, `#` comment header.
Columns: `mjd filter flux flux_err ZP flag mag mag_err field_id rcid
flux_offset offset_unc err_scale in_baseline`. Filters are `ztfg/ztfr/ztfi`;
`ZP = 30`; flux is difference-imaging flux with negatives retained; `mag = 99`
marks a non-detection.

**Flag handling — do not get this wrong.** `flag` is a bitmask, and `flag == 0`
is *not* the right cut. The official `ztfcosmo` default excludes only bits
`[1,2,4,8,16]`, i.e. **`flag & 31 == 0`**:

| bit | meaning |
|---|---|
| 1 | `flux_err == 0`, unphysical error |
| 2 | `chi2dof > 3`, extreme outlier |
| 4 | `cloudy > 1` |
| 8 | `infobits > 0` |
| 16 | `mag_lim < 19.3` |

Bits 32–1024 encode seeing, field, moon illumination, airmass and detection
significance, and are informational — they are **not** excluded. Cutting on
`flag == 0` throws away most genuine detections (it preferentially keeps faint
baseline epochs) and silently reduces the usable sample to near zero.

**Milky Way extinction.** `mwebv` is known per SN, so MW reddening is a
*deterministic* translation in magnitude space and must be removed up front,
not fitted. Only host-galaxy reddening is a free parameter.

**Caveat:** 16 of the 669 SNe passing the metadata cuts have no light-curve
file in the archive, and 36 rows sample-wide have NaN `t0`. Handle both in
`dgsn.data` rather than downstream.

## Working agreements

- Scientific claims in docs must be traceable to a citation or to code that
  produced the number. No asserted values.
- Report negative results plainly. A ladder rung that fails to earn its
  parameter is a result.
- When a design decision above is genuinely wrong, say so and change it
  explicitly — do not work around it.
