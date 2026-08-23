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

## This file does not define the model

**Pruned 2026-08-23.** `CLAUDE.md` used to carry the per-SN parameter list, the ladder,
the gauge/assumption/prediction classification, and a long chain of modelling arguments.
All of it has moved to `docs/model-note-backlog.md` ("Moved out of `CLAUDE.md`"), and the
model itself is defined in `docs/tex/model-definition.tex`.

The reason is concrete rather than tidiness. The model is a work in progress, and while
two files both asserted it they drifted: this file simultaneously fixed the reddening
direction at `(1,-1)/sqrt2` **and** claimed that `gamma(0) = 0` makes `mu` exactly the peak
`g` magnitude, which are incompatible — the second holds only at `e_c = (0,-1)`. That
inconsistency was live for a day and was propagated into the note before it was caught.

So: **do not record modelling choices here.** Parameter counts, gauge fixings, ladder
rungs, and what a symbol means all belong in the note, beside the mathematics that makes
them true. What stays below is what does not change when the model does — the data, the
samples, the tooling, and how the documents are kept.

## Standing methodological commitments

These constrain how the work is done, not what the model is. They are stable; revising
one is a deliberate act, not a side effect.

1. **Curve in magnitude space, likelihood in native flux.** In magnitude space distance
   and reddening act as *translations*, hence isometries, hence `kappa` is invariant. In
   flux space distance is an isotropic dilation and extinction an *anisotropic* one —
   neither is an isometry, and the claim that `kappa` is intrinsic collapses. Extinction,
   not distance, forces the choice: an isotropic dilation is repairable by rescaling, an
   anisotropic one takes a circle to an ellipse and is not. But ZTF errors are Gaussian in
   flux and negative fluxes are meaningful, so the likelihood must use native flux.
   **Never convert data to magnitudes.** The full argument — including why dust is only
   *reduced* rather than removed, and why "phase-independent versus phase-dependent" is
   the wrong nuisance/signal criterion — is in the backlog and must not be re-derived here.

2. **SALT2 is strictly downstream.** An independent benchmark on the same objects, never
   used to preprocess, K-correct, interpolate, or initialise. Using it upstream would make
   the independence claim false. One recorded exception, below: the epoch-counting window.

3. **`z < 0.05`.** A two-band model has no SED and cannot self-consistently K-correct. The
   cut buys independence at the cost of sample size. Do not relax it to gain statistics
   without revisiting the whole K-correction story.

4. **Evaluation is out-of-sample, always.** With per-SN parameters in the high single
   digits and a conditioned network, good in-sample fits are guaranteed and carry no
   evidential weight.

5. **External SED templates and extinction laws are for quantifying systematics only.**
   `Hsiao2007` and `Fitzpatrick1999` may be used to construct the reddening displacement
   and to evaluate its effect on `kappa` — never to fit, initialise, or K-correct.

6. **Code is written for general `R^n`.** Not speculative: the `R^3` subsample below
   exists and is the reason.

## Samples

| | |
|---|---|
| Data | ZTF SN Ia DR2 (`ztfcosmo`), `z < 0.05`, both quality flags. |
| Primary sample | **599 SNe** with >=5 good `g` and `r` epochs in the phase window. Ambient space `R^2`, one invariant. |
| Torsion subsample | **177 SNe** that additionally have >=5 good `i` epochs. Ambient space `R^3`, so a genuine second invariant. |
| Phase window | Rest-frame `[-15, +40]` d (config parameter). |

Coverage measured directly from the archive (2026-08-20), at `z<0.05` with
`lccoverage_flag` and `fitquality_flag` set, counting epochs with
`flag & 31 == 0` and SNR > 5 inside the phase window:

| min epochs/band | g | r | i | g&r | g&r&i |
|---|---|---|---|---|---|
| >=3 | 636 | 639 | 216 | 636 | 215 |
| >=5 | 606 | 631 | 179 | 599 | 177 |
| >=10 | 490 | 537 | 82 | 465 | 80 |

Two bands admit **no torsion**: in the primary `g,r` analysis a second shape latent is a
second direction of variation in `kappa`, not torsion, and must not be called torsion in
writing.

**Known exception to "SALT2 strictly downstream":** the phase window used to
*count* these epochs was placed with the archive's SALT2 `t0`, because the window
needs a date of maximum and the date of maximum is what the model fits. Selection only —
no photometry is transformed or initialised by it. Recount against an SED-free
anchor (window on the brightest good `g` epoch) to confirm the counts barely
move, and handle the 36 NaN-`t0` rows with that estimator or exclude them.

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
not fitted. Only host-galaxy reddening is a free quantity. Note that what is known is
`E(B-V)`, not the band extinctions themselves; converting needs an SED, so this removal
is an approximation too — milder, given the low Galactic columns at these latitudes, but
not exact.

**Caveat:** 16 of the 669 SNe passing the metadata cuts have no light-curve
file in the archive, and 36 rows sample-wide have NaN `t0`. Handle both in
`dgsn.data` rather than downstream.

### Null epochs — on file, not currently used

No pre-explosion epochs enter the likelihood as the model now stands; this is retained
because the turn-on rung is parked rather than abandoned (see the backlog), and because
the rank-one `offset_unc` item applies to **faint detections late in the window** too,
where it is live. Everything here rests on readings of the archive columns that are
**inferred from the column list and NOT YET VERIFIED against `ztfcosmo`**. Verify before
writing any of them into a document or relying on them in code; "no asserted values"
applies to the data model as much as to results.

- **"Expected flux is exactly zero" is false.** The expected pre-explosion flux is
  `flux_offset`, and zero only after it is subtracted. Subtract it deterministically —
  the same rule as `1+z` and `mwebv`.
- **`offset_unc` is rank-one**, one offset per light curve common to every epoch of a
  band, so **`chi2` is not diagonal**: `C = diag(sigma^2) + offset_unc^2 · 1·1^T`.
  Negligible beside a bright detection, **dominant** across a run of nulls. A diagonal
  treatment makes their joint constraint appear to tighten as `1/sqrt(N)` when the true
  floor is `offset_unc` and does not shrink. Marginalise analytically by
  Sherman–Morrison; adds no parameter.
- **Exclude `in_baseline` epochs.** They are the data the zero level was estimated from,
  so their residuals are shrunk toward zero by construction; using them as constraints on
  `f = 0` double-counts.
- **A null epoch set needs its own cut.** Drop `SNR > 5` — expected SNR is zero, so a
  threshold keeps only upward noise excursions, a one-sided selection on the noise
  realisation of the very quantity being fitted. And `flag & 31 == 0` is the wrong cut
  here, a change in kind rather than a tightening: it protects *detections*, whereas for
  a null *depth* is benign (a shallow epoch has large `sigma` and self-weights down under
  a correct Gaussian flux likelihood) while *bias* is fatal, and bits 32–1024 (seeing,
  field, moon, airmass) are precisely the bits correlated with systematic offsets. Select
  nulls on **provenance and error only, never on measured flux**. Pin down whether
  `sigma` is pre- or post-`err_scale`, which matters far more at zero signal.

The **selection window and its `SNR > 5` criterion are untouched**, so the 599/177
counts and the regression test asserting 599 stand.

## Stack

JAX · diffrax (Frenet ODE, adjoint) · equinox (networks) · optax ·
numpyro (posteriors, later).

Training is an **auto-decoder**: network weights and the per-SN latent array live in one
pytree and are optimised jointly.

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
  target of *shorter*: the new §2 material spent the whole saving. Next time the
  budget must be checked against a page count, not against deleted sections. The
  2026-08-22 removal of the zero-curvature segment took it from 15 pages back to **14**,
  which is the standing figure to check against.
- **The note is structured as a Model section, not an article** (2026-08-21): §1
  motivation, §2 the model (curve → placement → traversal → observables → gauge), then
  the supporting sections. Argument goes behind the model, never interleaved with
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

## Working agreements

- Scientific claims in docs must be traceable to a citation or to code that
  produced the number. No asserted values.
- Report negative results plainly. A ladder rung that fails to earn its
  parameter is a result.
- When a commitment above is genuinely wrong, say so and change it explicitly — do not
  work around it.
- Modelling choices go in the note, not here. See "This file does not define the model".
