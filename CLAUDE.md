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

## Per-SN parameters (baseline, 8)

`mu` (normalization) · `c` (colour) · `theta_1`, `theta_2` (shape, condition the
network) · `w`, `t_max` (timing) · `a_1`, `a_2` (timing warp)

**Why exactly two thetas (revised 2026-08-21).** Not because two was requested, and
not fixed by the geometry — `n` fixes the number of invariant *functions* (`n-1`), and
says nothing about how many latents condition them. The count comes from **SALT2
parity plus one increment**. SALT2 describes a SN with one time-dependent shape
coefficient `x1` and one time-independent colour `c`. Here colour is *not* a latent: it
is a translation, and translations leave `kappa` untouched, so it conditions nothing.
That leaves one latent for the SALT2 shape freedom and one more — the extra
time-dependent DOF this model exists to test. `K = 2` is therefore a working value with
a derivation, and whether the second earns its place is what the ladder answers.

**`c` is a baseline parameter, not an L2c addition (revised 2026-08-21).** Previously
the baseline had no free colour and `c` appeared only at L2c. That was wrong on its own
terms: `eq:tripp` was standardising on `beta*theta_2`, so `theta_2` was silently doing
SALT2's `c` job while also being billed as the new degree of freedom. Now `c` is free at
every rung. What changes along the ladder is the *function it multiplies*: split
`c*u(s) = c*ubar + c*du(s)`, the baseline uses the constant `ubar` (a pure translation,
exactly SALT2's `c`), and **L2c restores the full `u(s)` at no extra parameter**. L2c
tests whether the second-order `c*du(s)` term matters.

**The reddening direction `e_c`: gauge below L2c, a global parameter at L2c (revised
twice on 2026-08-22).** Previously the phase-independent reddening translation was fixed
along `(1,-1)` (i.e. `perp (1,1)`) and `ubar`, hence its direction, was imported from the
extinction law. The first revision made it a **single sample-wide fit parameter `e_c`**, a
unit vector in the `(m_g, m_r)` plane — the two-band analogue of SALT2 fitting `beta`
rather than adopting `R_V`, with `(1,-1)/sqrt2` the nested fixed-`R_V` idealisation. The
second revision sharpened this: **below L2c `e_c` is exactly unidentifiable**. With `mu`
along `(1,1)` and `c` along `e_c` both free per SN, the two span `R^2` for any `e_c` not
parallel to `(1,1)`, so changing the direction is undone exactly by an invertible
relabelling of `(mu, c)` — a **flat likelihood direction, not a slow one**. So `e_c` is
*gauge* at L0–L2 (fix it at `(1,-1)/sqrt2`; fitting it there is meaningless) and becomes a
genuine global *parameter* only at **L2c**, where the phase-varying `du(s)` gives it
content. Per-SN `c` remains the *amplitude* along `e_c`; only
the *phase variation* `du(s)` is still imported from the template (`Fitzpatrick1999` +
`Hsiao2007`) at L2c, with `e_c` setting its mean direction. A **per-SN** reddening
direction is *not* fitted — it would be a per-SN rotation of the colour axis with no
physical referent, the same reason no per-SN rotation is fitted (decision 4a). So the
two per-SN phase-independent shifts are `mu` along `(1,1)` (distance, all bands equal,
direction fixed by physics) and `c` along `e_c` (reddening, direction fitted globally).

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
   - **L2c fits the amplitude `c` of `u(s)`**, not a constant
     vector — a constant absorbs only `c*ubar` and leaves `c*du(s)` to contaminate
     `theta`, which is the failure L2c exists to prevent. Same parameter count. (Since
     the 2026-08-22 revision the mean *direction* of `u` is the fitted global `e_c`;
     only the phase-variation *shape* `du(s)` is precomputed.)
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

   An external SED template (`Hsiao2007`) and extinction law (`Fitzpatrick1999`)
   are allowed for *quantifying these systematics only* — constructing `u(s)`,
   evaluating `Delta kappa` — never to fit, initialise, or K-correct.

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

4. **No per-SN rotation; the orientation is gauge for `n = 2` and `n(n-1)/2 - 1` global
   parameters for `n >= 3`.** *(b)* and *(c)* were **revised twice on 2026-08-22**. The
   first revision made the orientation a fitted global parameter `Phi`; the second, made
   when the zero-curvature segment was dropped (decision 8, struck below), showed that
   **for `n = 2` no such parameter exists** — orientation is pure gauge. `Phi` is gone
   from the model. *(a)*, *(d)* and *(e)* survive unchanged — no *per-SN* rotation is
   fitted, and none of this reopens that.
   *(a)* No **per-SN** rotation is a physical degree of freedom: distance and zeropoint
   errors are **translations**, stretch is a **reparameterisation**, and dust is a
   **phase-dependent displacement** (see decision 1 — *not* a translation). None acts
   as a per-SN rotation, so none is fitted per SN.
   *(b)* **The orientation/origin redundancy, and what it leaves.** `kappa(.;theta)`
   integrated with initial tangent angle `psi` generates the *same curve* as
   `kappa(.+a;theta)` integrated with `psi + int_0^a kappa`. So the network's own
   arclength origin and the frame orientation carry **one joint redundancy**: of the
   `n(n-1)/2` rotational numbers in the Frenet initial data, one is absorbed by sliding
   the network origin, leaving **`n(n-1)/2 - 1` global orientation parameters** — **zero
   for `n = 2`**, **two for `n = 3`**. This is why `Phi` was not added capacity but a
   coordinate on a redundancy, and it **closes** the former `R^3` open frame question as
   a definite count rather than a missing condition.
   The arclength origin is fixed instead at the **principal `g`-maximum**,
   `dm_g/ds|_{s=0} = 0` with `T(0) = (0,-1)` and `kappa(0) > 0` (`eq:frameanchor` in the
   note). This is **gauge**: every turning planar curve has a vertical tangent, so it
   costs no generality, and it does *not* constrain which band peaks first. The sign
   check: `m_g'' = -sin(phi) kappa = kappa(0) > 0`, a **minimum of magnitude**, i.e. peak
   brightness. The mirror convention `T(0) = (0,+1)`, `kappa(0) < 0` is equivalent.
   Because `s = 0` is the `g`-maximum, `t_max` maps to `s = 0` and **there is no
   `s_g(theta)` root-find** anywhere in the model.
   *(c)* **The "which band peaks first" / frame-sign prediction is dropped.** It is not a
   quantity of interest. Do not reintroduce a `sign(kappa)` peak-ordering diagnostic or a
   max-to-max ordering report; the branch table of decision 5 is deleted for the same
   reason.
   *(d)* **Coincident maxima break arclength.** If all bands peak at the same
   epoch then `||dgamma/dp|| = 0`, so `ds/dt = 0` and unit-speed parameterisation
   fails. Non-coincident band maxima is a regularity condition on the data.
   *(e)* The anchor is `g`-maximum, **not** `B`-maximum. Locating `B` max needs
   an SED. So `t_max` is not comparable to SALT2 `t0`; measure the offset
   distribution, assume neither its centre nor its sign.

5. **Curvature sum rule — qualitative only (`∫ kappa ds ≈ pi`).** Flux vanishes before
   explosion and after, so `m -> +inf` in every band at both ends: the curve is a
   **hairpin**, both ends escaping to infinity, and the tangent reverses between them,
   turning through **approximately `pi`**. Corollaries that do *not* depend on the
   escape direction: `kappa -> 0` at both extremes, and the asymptotic ray separation is
   the terminal colour — predicted only *relative* to the early colour, since the free
   `c` shifts both together.

   **Demoted to a qualitative diagnostic (revised twice on 2026-08-22).** The old exact
   `∫ kappa ds = pi (mod 2pi)`, the exact `(1,1)` asymptote directions, and the branch
   table splitting the integral at the `g`-maximum are all **deleted** — the branch table
   was peak-ordering machinery, and 4(c) drops that. With decision 8 struck, **both**
   asymptotes now lie outside the model domain (`eq:domain`, a bounded interval of
   arclength) and neither is ever observed, so **nothing here is imposed at all**:
   hairpin, `kappa -> 0` at the extremes, `∫ kappa ds ≈ pi` and the terminal colour are
   qualitative validation statements about a trained `kappa`'s *extrapolation*.
   **Never impose any of them as a constraint.** By unit speed, `f = 0` requires
   `|s| = infinity` — the dark phase is a single ideal point in the closure, never a point
   of the domain.

6. **Peculiar velocity is not in the likelihood.** `mu` is free per SN and
   absorbs it exactly. It enters only at the Hubble-diagram stage.

7. **Evaluation is out-of-sample, always.** With 8 per-SN parameters and a conditioned
   network, good in-sample fits are guaranteed and carry no evidential weight.

8. ~~**A built-in early zero-curvature segment, a singular traversal, and
   pre-explosion epochs in the likelihood**~~ — **STRUCK IN FULL on 2026-08-22, one day
   after being added.** All six parts (a)–(f) are out of the model. It was listed
   non-negotiable, so it is struck explicitly rather than worked around.

   **Why.** The segment's own defining equation was inconsistent:
   `s in (-inf, s_end]` with `kappa_i = 0` for `s <= 0` truncates the *late* end while
   the unit-speed argument that motivated the segment forces `(-inf, +inf)`. Worse, it
   did not do the job it was carried for: nothing forbade `kappa == 0` on `(0, s_1]`
   too, so the "end of the segment" was not a well-defined arclength origin and the
   origin-versus-translation flat direction it was supposed to remove survived. And the
   dark phase does not need a semi-infinite *ray*: by unit speed `f = 0` is the single
   ideal point `s = -infinity` in the closure. The semi-infinite ray was occupied by the
   *early rise*, which is finite in time and infinite in arclength — a fact about the
   traversal, not about the curve.

   **What replaces it.** A **bounded** arclength domain `s in [s_min, s_max]`
   (`eq:domain`), held in config and checked wide enough to cover the selection window
   for every sampled parameter value. Nothing is imposed on `kappa_i` inside it. The
   arclength origin is fixed by decision 4(b)'s `g`-maximum condition, which is an
   *intrinsic* feature of the curve and therefore actually fixes the gauge. The
   traversal is the plain cubic `s(t) = u[1 + a_1 u + a_2 u^2]`, `u = (t-t_max)/(w(1+z))`,
   with no logarithmic term. There is no `t_expl` in the model and no pre-explosion epoch
   in the likelihood; two windows, not three (model support = `eq:domain`; the selection
   window `[-15,+40]` d also delimits the likelihood epoch set).

   **The physics is parked, not deleted.** The turn-on — semi-infinite domain, singular
   traversal, derived or free `t_expl`, pre-explosion nulls, and the early power-law rise
   with `alpha_g/alpha_r` as a check on the orientation — is a coherent **candidate later
   rung**, recorded as open in `docs/model-note-backlog.md` and scored on held-out
   residuals like any other rung.

   **What survives independently of the segment**, and is now recorded where it belongs:
   the `n = 2` regularity of the Frenet rotation ODE at `kappa = 0` (`N` is *defined* as
   the fixed 90-degree rotation of `T`, so only `N ∝ gamma''/|gamma''|` fails); and for
   `n >= 3`, `tau` **undefined** — not merely ill-conditioned — wherever `kappa_1`
   vanishes, requiring a parallel-transport (**Bishop**) frame and an explicit
   unidentifiability report. Both are in Known limitations below.

## The deliverable is a ladder, not a fit

| Rung | Free per SN | # | Question |
|---|---|---|---|
| L0 | `mu, c, w, t_max` | 4 | **same count as SALT2** — does the mechanism alone win? |
| L1 | `+ theta_1` | 5 | does one shape parameter earn its place? |
| L2 | `+ theta_2` | 6 | does a second? |
| L2c | L2, `c` upgraded to full `u(s)` | 6 | does phase-dependent dust matter? |
| L3 | `+ a_1, a_2` | 8 | does nonlinear timing earn its place? |

**L0 is a controlled experiment, not a starting point.** It carries four parameters in
one-to-one correspondence with SALT2's `(x0, x1, c, t0)` and differs in *mechanism*
alone: SALT2 produces stretch-like variation from an additive component, L0 from an
exact reparameterisation. L0 versus SALT2 is therefore a like-for-like test at fixed
parameter count, before any new freedom enters. Everything above L0 asks whether added
freedom pays.

**The mechanism difference, stated precisely.** SALT2 has *no* way to reparameterise
time. A stretch `p -> p/s` leaves the curve pointwise unchanged and alters only its
traversal; SALT2 instead uses the additive `x1*M1(p,lambda)`. These agree only to first
order — `m(p/s) = m(p) - (s-1)*p*m'(p) + O((s-1)^2)`, so `M1 ∝ -p dM0/dp`. This is not
a reconstruction after the fact: **SALT2 initialises `M1` as exactly that finite
difference**, the sequence at `s=1.1` minus that at `s=1`. The published `x1 -> s`
conversion needs a cubic, and that cubic *is* the residual nonlinearity the
linearisation leaves behind.

Two consequences, both load-bearing:
- The dominant mode of SN Ia diversity is the one this model represents **exactly** and
  SALT2 **linearises**. That is the strongest structural argument for the approach.
- With stretch removed into `w`, what remains for `theta_1` is only shape variation no
  reparameterisation can produce, so **`theta_1` should be small**. Concrete check:
  fitted `w` should track `s = 0.98 + 0.091*x1 + 0.003*x1^2 - 0.00075*x1^3`, and the
  scatter about that relation measures what the linearisation misses. Run this before
  building the full pipeline — it is cheap and it is the sharpest early test of whether
  separating shape from timing buys anything.

Standardisation is `mu_corr = mu - alpha*theta_1 - beta*c - gamma*theta_2 + M`. The
first three terms are SALT2's own, one-to-one. **`gamma` is the headline number** — the
coefficient of the one DOF SALT2 lacks.

Scored on **held-out Hubble residual scatter** under cross-validation
(leave-one-out jackknife, following SALT2's own validation method).

The ladder is also the **training schedule**: fit L0 to convergence, then
introduce latents one at a time. SALT2 stages its components the same way.

## Stack

JAX · diffrax (Frenet ODE, adjoint) · equinox (networks) · optax ·
numpyro (posteriors, later).

Training is an **auto-decoder**: network weights and the `N x 8` array of
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
  it grow back. A new subsection is a cost that has to be justified. The 2026-08-21
  restructure into six sections deleted the ladder and limitations sections and came
  out **net-zero at 14 pages** (13 of body plus one of bibliography), missing its own
  target of *shorter*: the new §2 material — the singular traversal, the three windows,
  the rank-one covariance, the gauge table — spent the whole saving. Next time the
  budget must be checked against a page count, not against deleted sections. The
  2026-08-22 removal of the zero-curvature segment took it from 15 pages back to **14**,
  which is the standing figure to check against.
- **The note is now structured as a Model section, not an article** (2026-08-21): §1
  motivation, §2 the model (curve → placement → traversal → observables → gauge), then
  §3–§6 as supporting material. Argument goes behind the model, never interleaved with
  it, and the constitutive chain must read in dependency order with no forward
  references into the supporting sections.
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
- **Curvature degenerates at inflections.** With decision 8 struck there is no
  by-design zero-curvature set, so this is once again purely about the **incidental**
  near-zero curvature the `(m_g, m_r)` path may approach near the secondary maximum, plus
  whatever the trained `kappa` does near the domain ends. "Report where along `s` the
  curvature approaches zero" is still required. What it costs: for `n = 2`, **nothing** —
  the Frenet system is a rotation ODE, regular at `kappa = 0`, because `N` is *defined* as
  the fixed 90-degree rotation of `T` independently of `gamma''`; only
  `N ∝ gamma''/|gamma''|` fails, so integrating the ODE is the definition that survives.
  For `n >= 3` it is worse than ill-conditioning: where `kappa_1 = 0` there is no
  osculating plane, so `tau` is **undefined**, whatever the network emits there is
  **unidentifiable and must be reported as such**, and the frame must be carried by
  parallel transport (a **Bishop frame**), not by the classical construction and not
  merely an "inflection-robust" one.
- **Gauge, assumption and prediction must be distinguished** before latents are
  interpretable. Four kinds; the classification was revised again on **2026-08-22** when
  decision 8 was struck, and the net movement is *toward* gauge.
  - *Gauge* (costless, not testable): `s = 0` at the **principal `g`-maximum**,
    `dm_g/ds|_{s=0} = 0` with `T(0) = (0,-1)`, `kappa(0) > 0` (decision 4b) — this
    replaces "end of the zero-curvature segment" and, unlike it, actually fixes the
    origin; the **frame orientation** for `n = 2`, by the same condition; the metric
    `diag(1,1)`; latents normalised to zero mean and unit variance, with the sign fixed so
    `theta_1` correlates positively with light-curve width; position fixed by normalising
    the template to `gamma(0) = 0`, which makes `mu` *exactly* the peak `g` magnitude and
    `c` the colour at that epoch; and the **reddening direction `e_c` below L2c**, where it
    is exactly unidentifiable.
  - *Assumption* (substantive, falsifiable): `Corr(c, theta) = 0`. **`t_expl` scaling with
    `w` has left this list** — there is no `t_expl` in the model since decision 8 was
    struck. **Frame orientation is not here either**, in either direction: it is gauge for
    `n = 2` and a parameter for `n >= 3`, never an assumption.
  - *Prediction* (checked, not chosen): colour **evolution**; the hairpin, `∫ kappa ds ≈
    pi` and the terminal colour, all as *qualitative* statements about extrapolation
    outside `eq:domain` (decision 5). **Not** the early power-law index — that left with
    decision 8 and is now part of the parked turn-on rung. **Not** "which band peaks
    first" — dropped by 4(c).
  - *Parameter*: `c`, the amplitude along `e_c`, free per SN (**revised 2026-08-21** from a
    template function of `theta`); `e_c` itself **at L2c only**, one global unit vector for
    the whole sample, not per-SN; and for `n >= 3` the **`n(n-1)/2 - 1` global orientation
    numbers** (two for `n = 3`), one fewer than `SO(n)` because of the orientation/origin
    redundancy of 4(b). With `c` free, neither early nor terminal colour is predicted
    absolutely; only their *difference* is, since `c` displaces the whole curve and
    cancels. The model predicts colour evolution and leaves the colour zero point free —
    structurally the same split SALT2 makes between a fixed colour law and a fitted `c`.
- **The network's arclength origin and the frame orientation are one joint redundancy.**
  `kappa(.;theta)` with initial tangent angle `psi` generates the same curve as
  `kappa(.+a;theta)` with `psi + int_0^a kappa`. Two consequences, and the second replaces
  the old straight-segment degeneracy bullet: (i) for `n = 2` there is **no** free
  orientation parameter at all, and for `n >= 3` there are `n(n-1)/2 - 1`; (ii) the
  arclength origin must be fixed at an **intrinsic feature** of the curve — the
  `g`-maximum — not at a location defined by a construction, or a genuine null direction of
  the likelihood survives. It must be fixed, not regularised. `mu` stays an **all-band**
  offset — a `g`-only offset would mix the distance and colour channels and would revise
  decision 1.

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

**Pre-explosion (null) epochs — kept on file, but downstream of a parked rung.** With
decision 8 struck there are **no pre-explosion epochs in the likelihood**, so none of the
following is needed by the current model. It is retained because the turn-on rung is parked
rather than abandoned, and because the second item (rank-one `offset_unc`) applies to
**faint detections late in the window** too, where it is still live. Everything below rests
on readings of the archive columns that are **inferred from the column list and NOT YET
VERIFIED against `ztfcosmo`**. Verify before writing any of them into a document or relying
on them in code; "no asserted values" applies to the data model as much as to results.

- **"Expected flux is exactly zero" is false.** The expected pre-explosion flux is
  `flux_offset`, and zero only after it is subtracted. Subtract it
  deterministically — the same rule as `1+z` and `mwebv`.
- **`offset_unc` is rank-one**, one offset per light curve common to every epoch of a
  band, so **`chi2` is not diagonal**: `C = diag(sigma^2) + offset_unc^2 · 1·1^T`.
  Negligible beside a bright detection, **dominant** across a run of nulls. A diagonal
  treatment makes their joint constraint appear to tighten as `1/sqrt(N_pre)` when the
  true floor is `offset_unc` and does not shrink — overstating the information by a
  factor growing with `N_pre` and biasing `t_expl` in the sign of `flux_offset`.
  Marginalise analytically by Sherman–Morrison; adds no parameter.
- **Exclude `in_baseline` epochs.** They are the data the zero level was estimated
  from, so their residuals are shrunk toward zero by construction; using them as
  constraints on `f = 0` double-counts, giving an over-tight `t_expl` and an
  understated chi-squared.
- **Drop the `SNR > 5` cut, and `flag & 31 == 0` is the wrong cut here — a change in
  kind, not a tightening.** Expected SNR is zero, so an SNR threshold keeps only
  upward noise excursions: a one-sided selection on the noise realisation of the very
  quantity being fitted, biasing `t_expl` early. Select nulls on **provenance and
  error only, never on measured flux**. And the recorded flag warning runs the other
  way — it protects *detections*. For a null the failure mode inverts: *depth* is
  benign (a shallow epoch has large `sigma` and self-weights down under a correct
  Gaussian flux likelihood, so bit 16 costs information, not correctness), while
  *bias* is fatal, and bits 32–1024 (seeing, field, moon, airmass) are precisely the
  bits correlated with systematic offsets. So the null epoch set needs its own cut,
  derived for that purpose. Pin down whether `sigma` is pre- or post-`err_scale`,
  which matters far more at zero signal than for detections.

The **selection window and its `SNR > 5` criterion are untouched**, so the 599/177
counts and the regression test asserting 599 stand.

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
