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

**Known exception to "SALT2 strictly downstream":** the phase window used to
*count* these epochs was placed with the archive's SALT2 `t0`, because the window
needs a date of maximum and `t_max` is what the model fits. Selection only — no
photometry is transformed or initialised by it. Recount against an SED-free
anchor (window on the brightest good `g` epoch) to confirm the counts barely
move, and handle the 36 NaN-`t0` rows with that estimator or exclude them.

## Per-SN parameters (baseline, 7)

`mu` (normalization) · `theta_1`, `theta_2` (shape, condition the network) ·
`w`, `t_max` (timing) · `a_1`, `a_2` (timing warp)

`w` is the timing scale in **rest-frame days**; `t_max` is in **observer-frame
MJD** and is always free per SN, never taken from the archive's SALT2 `t0`.
`sigma` is reserved for flux measurement errors and must not be used for the
timing scale. Time dilation is deterministic: rest-frame phase is
`p = (t - t_max)/(1+z)`, divided out up front like `mwebv`, so `w` is a genuine
rest-frame stretch and not a repackaging of `(1+z)`.

SALT2 uses four. The excess is deliberate and is the object of study, not a
free choice.

## Non-negotiable design decisions

These were settled deliberately. Do not silently revise them.

1. **Curve in magnitude space, likelihood in flux.** In magnitude space distance
   (`(1,1)/sqrt2`) and colour (`(1,-1)/sqrt2`) are *translations*, hence
   isometries, hence `kappa` is invariant. In flux space distance is an isotropic
   dilation and extinction an *anisotropic* one — neither is an isometry, and the
   claim that `kappa` is intrinsic collapses. Extinction, not distance, is what
   forces the choice: an isotropic dilation is repairable by rescaling, an
   anisotropic one takes a circle to an ellipse and is not. But ZTF errors are
   Gaussian in flux and negative fluxes are meaningful, so the likelihood must
   use native flux. **Never convert data to magnitudes.**

   **Two logical types, not one axis. Know which side an effect is on.**
   - **Rigid displacements — removed exactly.** Distance, absolute luminosity,
     peculiar velocity, per-band calibration error. Each is a fixed property of
     the object or the telescope, so it applies the *same* displacement at every
     phase: a map of the plane, hence a translation, hence invisible to `kappa`
     with no residual. Note this is *not* about being achromatic — a zeropoint
     error is strongly chromatic and still exactly removed.
   - **Deformations — not removed by any geometry.** Dust *and intrinsic SN
     diversity*, i.e. the signal itself. Both displace different parts of the
     curve by different amounts, so both bend it. **Do not use
     "phase-independent vs phase-dependent" as the nuisance/signal criterion**:
     it sorts dust onto the same side as the thing `kappa(s;theta)` exists to
     measure. Separation can only come from population-level assumptions, which
     are substantive and falsifiable — never from the choice of ambient space.

   **Same type ≠ same size. The shape contamination is second order.** Split
   `c*u(s) = c*ubar + c*du(s)` with `ubar = <u>` and `<du> = 0`. The first term is
   a constant vector, hence a translation, hence removed **exactly**. Only `c*du(s)`
   bends the curve, and it is small twice over: in the reddening amplitude `c`, and
   in the fractional phase variation `|du|/|ubar|`, which is small because a
   band-integrated extinction responds only weakly to SED evolution inside the
   filter. So there is a **hierarchy**: intrinsic diversity enters `kappa` at first
   order (it *is* the signal), dust at second. **Expect the phase-dependent variance
   the geometry sees to be dominated by intrinsic SN dispersion.** This is what
   makes `theta` interpretable despite the shared category — but it is a
   *quantitative* claim to be measured, not asserted, and the diagnostics below
   exist to test it.

   Two caveats. (i) Amplitude-dependent: `c` is not small for the reddest objects,
   so **stratify every diagnostic by fitted `c`** rather than sample-averaging; the
   hierarchy degrades at the reddened end. (ii) **Shape channel only.** In the
   *colour* channel dust and intrinsic colour mix at **first** order — both shift
   `ubar` — and that degeneracy is not second order, is not resolved here, and is
   the same one SALT2 confronts through `beta`.

   **So what magnitude space actually buys:** (i) distance, the quantity being
   inferred, is made *exactly* orthogonal to shape — this is why `mu` can be free
   per SN and absorb peculiar velocity exactly; (ii) dust is *reduced* from a
   shape-destroying anisotropic dilation to a small deformation about a rigid
   part. Reduced, not removed. Do not overclaim `kappa` as nuisance-free.

   **Dust in detail.** Band extinction `A_X(p)` is a flux-weighted average of the
   monochromatic law over the filter, weighted by the SN's own SED — which
   evolves with phase. So dust is a *phase-dependent displacement*
   `gamma -> gamma + c*u(s)`, not a translation, and it does change `kappa` and
   arclength. Three consequences:
   - **`z < 0.05` does not help.** This is a *rest-frame* effect, identical at
     `z=0`. Low `z` suppresses the K-correction branch of the missing-SED
     problem and does nothing to the dust branch.
   - **`mwebv` removal is approximate too.** `E(B-V)` is known; `A_g(p)`,
     `A_r(p)` are not without an SED. Milder (low Galactic columns), not exact.
   - **L2c fits the amplitude `c` of a precomputed `u(s)`**, not a constant
     vector — a constant absorbs only `c*ubar` and leaves `c*du(s)` to contaminate
     `theta`, which is the failure L2c exists to prevent. Same parameter count.
     But by the order counting this is **cheap insurance, not a necessity**: a
     constant would capture most of the displacement, and the refinement matters
     only for the reddest objects. Adopted because it costs nothing. Cost:
     `kappa` is not intrinsic under dust; reddening becomes a known deformation
     in the forward model, not a quotiented symmetry.

   **Mixing with the signal is tolerable for the distance indicator — do not
   over-escalate it.** Since dust and intrinsic diversity are the same type,
   `theta` will carry some reddening. That does *not* break the deliverable: for
   held-out Hubble residual scatter, isolating physical causes is not required,
   and the fitted Tripp `alpha, beta` absorb the mixture. SALT2 is in the same
   position — its `c` mixes dust with intrinsic colour, which is why `beta` is
   fitted rather than set to `R_V + 1`. The requirement is only that the mixture
   be *stable across the sample*; `z < 0.05` gives too little redshift span for a
   trend, so there is no route to cosmological bias **here**. Revisit before any
   wide-redshift application. Contamination threatens *interpretation* of
   `theta`, not the distance indicator — report the two conclusions separately.

   **Identifying assumption, if `theta` is to be interpreted.** Line-of-sight dust
   is a property of the intervening ISM with no causal link to the explosion, so
   impose **zero sample correlation between fitted `c` and fitted `theta`**. This
   is an *identifying assumption* — substantive, falsifiable, possibly wrong —
   and is a different kind of thing from the latent zero-mean/unit-covariance
   normalisation, which is a costless gauge. Known to be only approximate: the
   host mass step couples SN properties to host environment, and dust column is
   environmental.

   Whether the assumption is load-bearing depends on the order counting. If the
   hierarchy holds, the data fix the shape channel at first order and the condition
   is a **check** — `Corr(c, theta)` should come out small whether or not it is
   imposed. If the hierarchy fails, the condition holds the analysis up and
   conclusions about `theta` rest on it. So **impose it and also report the
   unimposed value**.

   **Required diagnostics before interpreting `theta`:**
   - `Corr(c, theta)` over the population, computed **both with and without** the
     identifying condition imposed. A large unconstrained value falsifies the order
     counting.
   - **Error-weighted overlap** of `du(s)` — the *varying* part, **not** the full
     `u(s)`; the mean part is a translation and cannot be confused with shape, so
     using full `u` overstates the overlap — with
     `span{dgamma/dtheta_1, dgamma/dtheta_2}`, weighted by flux errors at the
     epochs *actually observed*. This is the honest degeneracy statement: it asks
     whether this cadence and these error bars can tell them apart, not whether
     they differ in principle.
   - `Delta kappa(s)` from propagating the `du(s)` deformation through Frenet,
     compared to the `kappa` range the network spans over fitted `theta`. **This is
     the direct test of the order counting.** Report it **stratified by `c`**. Small
     everywhere ⇒ picture confirmed, `theta` interpretable. Order unity at any `c`
     ⇒ the shape latents are substantially measuring dust there: a negative result
     for interpretation, to be reported as one.

   An external SED template (e.g. Hsiao) is allowed for *quantifying these
   systematics only* — never to fit, initialise, or K-correct.

   **SALT2 is structurally better placed on dust, and that is an honest cost.**
   SALT2's colour *law* is phase-independent by construction, but SALT2 carries an
   SED, so `A_X(p)` comes out phase-dependent automatically. This model has no SED
   and must import `u(s)` from an external template.

   **Metric is a convention.** `kappa` presupposes a metric on magnitude space
   and none is physically preferred. `diag(1,1)` is a choice; it does not affect
   invariance (translations are isometries of any translation-invariant metric)
   but does set what "shape" means. Declared in the gauge-fixing section.

2. **SALT2 is strictly downstream.** Used only as an independent benchmark on
   the same objects. Never to preprocess, K-correct, interpolate, or
   initialise. Using it upstream would make the independence claim false.

3. **`z < 0.05`.** A two-band model has no SED and cannot self-consistently
   K-correct. The cut buys independence at the cost of sample size. Do not
   relax it to gain statistics without revisiting the whole K-correction story.

4. **No per-SN rotation, and the frame is fixed by the g-maximum.** Two
   separate claims, both needed.
   *(a)* Rotation is not a physical degree of freedom: distance and zeropoint
   errors are **translations**, stretch is a **reparameterisation**, and dust is a
   **phase-dependent displacement** (see decision 1 — *not* a translation).
   Nothing acts as a rotation, so none is fitted.
   *(b)* The template frame is fixed by the same anchor as the arclength
   origin — at the `g`-maximum `dm_g/ds = 0`, so `T(0)` is `+/-(0,1)`, with the
   sign empirical per *(c)* below (expected `(0,-1)`). One condition fixes both
   gauges, costs no generality, and is observable inside
   the phase window. **This one-condition claim holds only for `n=2`**: fixing
   the frame means fixing an element of `SO(n)`, which has `n(n-1)/2` parameters,
   so the `g`-max condition is enough in `R^2` but leaves two conditions unfixed
   in the `R^3` torsion subsample. Those are an **open question**, not settled.
   *(c)* The residual sign `T(0) = (0,-1)` vs `(0,+1)` is **empirical, not
   gauge**. Arclength increases with time and brightening means *decreasing*
   magnitude, so `(0,-1)` says `r` is still brightening at `g`-max and `(0,+1)`
   says it is already fading — i.e. the sign records **which band peaks first**.
   Determine it by direct measurement (per-band polynomial fit near each peak,
   record the order); carry it in config, do not hard-code. Do **not** infer it
   from effective wavelength: peak epoch is *not* monotonic in wavelength across
   UV to NIR. The sum rule gives no independent handle — the *sense* of turning
   (`kappa>0` vs `kappa<0`) is exactly equivalent to the peak ordering, so it
   restates the question. Since the frame is global and no rotation is fitted, a
   mixed-ordering sample would be evidence *against* the no-rotation design;
   report the minority fraction and the distribution of max-to-max separation.
   *(d)* **Coincident maxima break arclength.** If all bands peak at the same
   epoch then `||dgamma/dp|| = 0`, so `ds/dt = 0` and unit-speed parameterisation
   fails. Non-coincident band maxima is a regularity condition on the data.
   *(e)* The anchor is `g`-maximum, **not** `B`-maximum. Locating `B` max needs
   an SED. So `t_max` is not comparable to SALT2 `t0`; measure the offset
   distribution, assume neither its centre nor its sign.

5. **Curvature sum rule (`∫ kappa ds = pi mod 2pi`).** Flux vanishes before
   explosion and after, so `m -> +inf` in every band at both ends: the curve is
   a **hairpin** whose two asymptotes are parallel rays along `(1,1)`. The
   tangent turns from `-(1,1)/sqrt2` to `+(1,1)/sqrt2`, a net turning of `pi`.
   Corollaries: `kappa -> 0` at both extremes, and the asymptotic ray
   separation is the terminal colour.

   The asymptotes lie **outside** the `[-15,+40]` d window, so this constrains
   the *extrapolation* of a trained `kappa`, not the fit. Use it as a
   validation diagnostic, or at most a weak prior. **Never impose it as a hard
   constraint** on a fit that cannot see the asymptotes.

6. **Peculiar velocity is not in the likelihood.** `mu` is free per SN and
   absorbs it exactly. It enters only at the Hubble-diagram stage.

7. **Evaluation is out-of-sample, always.** With 7 latents and a conditioned
   network, good in-sample fits are guaranteed and carry no evidential weight.

## The deliverable is a ladder, not a fit

| Rung | Free per SN | Question |
|---|---|---|
| L0 | `mu, w, t_max` | is a rigid template enough? |
| L1 | `+ theta_1` | does one shape parameter earn its place? |
| L2 | `+ theta_2` | does a second? |
| L2c | L2 `+ c` amplitude of `u(s)` | does dust need its own direction? |
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
docs/superpowers/specs/
             Markdown engineering design docs
references/  source papers
```

Planned but **not yet created**: `configs/` (one YAML per ladder rung). Everything
under `src/dgsn/` is currently empty scaffolding — the project is at the
design-document stage.

## Documentation convention

- **Scientific content → LaTeX** in `docs/tex/`. Always `\input{macros}` so
  notation stays consistent; cite from `refs.bib` (entries sourced from NASA
  ADS via the `nasa-ads` MCP server).
- **Engineering content → Markdown**. This file, `README.md`, and the design docs
  under `docs/superpowers/specs/`.
- **`docs/model-note-backlog.md` is a deliberate exception** to "scientific content →
  LaTeX". See the pair of purposes below.

### The two documents and what each is for

These purposes are settled and are the reason the note looks the way it does. Read
them before editing either file.

**`docs/tex/model-definition.tex` — a working note, for the author, to understand the
model.** It is not a paper and has no external audience yet. Four consequences, all
load-bearing:

- **Prose for an astronomer, substance in mathematics.** Definitions, equations, short
  propositions. When a point needs a paragraph of argument to land, that is the signal
  either that it belongs in the backlog or that the claim should be restated as a
  formula. Converting argument to mathematics is usually both shorter *and* stronger —
  e.g. the first variation `Δκ = ψ'' + κ²ψ + κ'φ` replaced four paragraphs on rigid
  displacements versus deformations, and does more than they did, since translations
  drop out of it identically.
- **No abstract, deliberately.** The note is evolving; an abstract is maintenance
  overhead with no reader. **Do not add one back** until asked.
- **Length is a constraint, not an accident.** It was cut from 18 pages to 14 on
  2026-08-20 precisely because it had grown too long to serve its purpose. Do not let
  it grow back. A new subsection is a cost that has to be justified.
- **Say each thing once.** The bloat that forced the cut was the same argument restated
  in four or five sections. Cross-reference instead; if a section needs a point already
  made, it references the equation and adds only what is new there.

**`docs/model-note-backlog.md` — a staging area, not a document.** Argumentative
material raised while the note evolved, parked so the note can stay short. Each entry
records the claim, why it was raised, and whether it is settled or open. Traffic runs
both ways and neither direction is automatic: material moves *out* of the note when it
turns out to be argument rather than model, and *back in* only when it is both settled
and expressible as mathematics. Do not let it become a second note.

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
- **Dust vs shape degeneracy.** The same structure but expected milder: dust and
  intrinsic diversity are both deformations, so no geometric argument separates
  them (decision 1), yet only `c*du(s)` competes with shape and that is second
  order where intrinsic diversity is first. `c` is degenerate against `theta` to
  the extent `du(s)` lies in `span{dgamma/dtheta}` at the observed epochs.
  Expected weak; the `Corr(c, theta) = 0` assumption is a check on that, not the
  thing holding it up. Consequences are **asymmetric** — harmless for the distance
  indicator, a limitation for interpreting `theta`. Report both, stratified by `c`.
- **Dust vs intrinsic colour is the first-order degeneracy, and is unresolved.**
  `c*ubar` and a per-SN intrinsic colour offset are the same displacement, so `c`
  conflates the two exactly as SALT2's `c` does. The order counting does *not*
  help here. Do not let the good news about the shape channel obscure this.
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
