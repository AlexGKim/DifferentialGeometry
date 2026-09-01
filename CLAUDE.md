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

The project is exploratory. This file records what does not change when the model
does — the data, the sample, the tooling — and deliberately does not impose much
else. **Modelling choices belong in `docs/tex/model-definition.tex`**, beside the
mathematics that makes them true, not here; two files asserting the model have
drifted before.

## Standing commitments

Only these three. Others were retired on 2026-08-28 and are recorded in the
backlog under "Retired from `CLAUDE.md` as standing rules", with the reasoning
intact should any need restoring — notably **SALT2 strictly downstream** and
**evaluation is out-of-sample**, both of which must be re-established before any
result is claimed.

1. **Curve in magnitude space, likelihood in native flux.** In magnitude space
   distance and reddening act as translations, hence isometries, hence `kappa` is
   invariant. In flux space extinction is an *anisotropic* dilation, which takes a
   circle to an ellipse and is not repairable — so the claim that `kappa` is
   intrinsic collapses. But ZTF errors are Gaussian in flux and negative fluxes are
   meaningful, so the likelihood must use native flux. **Never convert data to
   magnitudes.** Full argument in the backlog.

2. **`z < 0.05`.** A two-band model has no SED and cannot self-consistently
   K-correct; the cut keeps the correction small enough to neglect. Relaxing it
   means revisiting the K-correction story.

3. **External SED templates and extinction laws quantify systematics only.**
   `Hsiao2007` and `Fitzpatrick1999` may be used to construct the reddening
   displacement and evaluate its effect on `kappa` — never to fit or initialise.

## Sample

`g` and `r` only; the `i` band and the `R^3` torsion subsample were dropped on
2026-08-28. Regenerate every number below with `python -m dgsn.data.sample`.

| | |
|---|---|
| Data | ZTF SN Ia DR2 (`ztfcosmo`), `z < 0.05`, both quality flags `== 1`. |
| Selected | **652 SNe**; 14 have no light-curve file, leaving **638**. |
| Analysis sample | **598 SNe** with >=5 good epochs in *each* of `g` and `r`. |
| Phase window | Rest-frame `[-15, +40]` d. |

| min epochs/band | g | r | i | g&r |
|---|---|---|---|---|
| >=3 | 635 | 638 | 215 | 635 |
| >=5 | 605 | 630 | 178 | 598 |
| >=10 | 489 | 536 | 82 | 464 |

The 14 without light-curve files are **exactly** the 14 with NaN SALT2 `t0` — the
absence is of photometry, not of shape.

Two bands admit **no torsion**: a second shape latent is a second direction of
variation in `kappa`, not torsion, and must not be called torsion in writing.

## Data details

A local copy of the archive is at `data/ztfsniadr2_lite.zip` (1.4 GB, gitignored
along with everything else under `data/`), unpacked at
`data/raw/ztfsniadr2_lite/`. The maintained interface is
[`ztfcosmo`](https://github.com/ZwickyTransientFacility/ztfcosmo), which needs no
credentials — it reads remotely from `ztfcosmo.in2p3.fr`, or from `$ZTFCOSMODIR`.

**`tables/snia_data.csv`** — 3628 rows. Columns include `ztfname`, `redshift`,
SALT2 fits `t0, x0, x1, c` with full covariances, `mwebv`, `sn_type`,
`lccoverage_flag`, `fitquality_flag`.

**`lightcurves/<name>_lc.csv`** — whitespace-delimited, `#` comment header.
Columns: `mjd filter flux flux_err ZP flag mag mag_err field_id rcid
flux_offset offset_unc err_scale in_baseline`. Filters are `ztfg/ztfr/ztfi`;
`ZP = 30`; flux is difference-imaging flux with negatives retained; `mag = 99`
marks a non-detection.

**Quality flags are NaN for 18 objects, and `!= 0` admits them.** `NaN != 0` is
True in pandas, so `lccoverage_flag != 0` silently keeps objects whose quality is
unknown. This produced a wrong selection of 669 that stood for a week. Require
`== 1`.

**Carry the selection through to the counting stage.** The epoch counts were once
produced by iterating the light-curve files on disk rather than the selected
objects, so the `z < 0.05` cut was applied when the file list was built but not
again when epochs were counted. One supernova at `z = 0.0925` (ZTF18abkhdxe, 38
good `g` and 122 good `r` epochs) passed every threshold and inflated all six
counts by one. Selection is an object list, not a directory listing.

**In a light curve, `filter` is a DataFrame method.** `lc.filter == "ztfg"`
compares a bound method to a string and is False for every row, silently yielding
zero epochs. Use `lc["filter"]`.

**Epoch `flag` is a bitmask, and `flag == 0` is *not* the right cut.** The
official `ztfcosmo` default excludes only bits `[1,2,4,8,16]`, i.e.
**`flag & 31 == 0`**:

| bit | meaning |
|---|---|
| 1 | `flux_err == 0`, unphysical error |
| 2 | `chi2dof > 3`, extreme outlier |
| 4 | `cloudy > 1` |
| 8 | `infobits > 0` |
| 16 | `mag_lim < 19.3` |

Bits 32–1024 encode seeing, field, moon illumination, airmass and detection
significance, and are informational — they are **not** excluded. Cutting on
`flag == 0` throws away most genuine detections and silently reduces the usable
sample to near zero.

**Milky Way extinction is not removed from the photometry.** The data reach the
fit as observed; the Galactic term is absorbed by the placement parameters and
taken out downstream. Why that is legitimate is in the note — do not restate it
here. `mwebv` is the SFD `E(B-V)` at the SN coordinates. Not a small-extinction
sample: median 0.043, 24% above 0.1, 5.5% above 0.3, max 1.06, a fifth at
`|b| < 20°`. No latitude or reddening cut is imposed.

**Do not apply the 0.86 rescaling twice.** Schlafly & Finkbeiner (2011)
coefficients (`A_g = 3.303 E(B-V)`, `A_r = 2.285 E(B-V)` for SDSS-like bands;
compute for ZTF's own filters rather than borrowing) are defined against the
**raw** SFD map value and already include the ~14% recalibration. Feed them
uncorrected `mwebv`. Rescaling the map by 0.86 *and* using SF11 coefficients
double-counts; the older SFD98 coefficients (`R_g, R_r = 3.793, 2.751`)
systematically over-correct. Bites downstream, where the correction is applied.

## Stack

JAX · diffrax (Frenet ODE, adjoint) · optax · numpyro (posteriors, later).

No neural networks: the global functions are basis expansions, so there is
nothing for `equinox` to do. Training is an **auto-decoder** — basis
coefficients and the per-SN latent array live in one pytree and are optimised
jointly — with a linear decoder.

Two implementations of the geometry, which must agree on `kappa`:
- `geometry/frenet.py` — diffrax ODE integration. **Primary.**
- `geometry/direct.py` — parameterise `gamma(s)` directly, recover `kappa` by
  autodiff. Fast correctness oracle and regression test.

## Layout

```
src/dgsn/
  data/      sample.py — selection, quality cuts, epoch counting
  geometry/  frenet.py, direct.py, reparam.py, invariants.py
  model/     kappa_net.py, forward.py, likelihood.py
  train/     auto-decoder loop, staged schedule
  infer/     per-SN fits, numpyro posteriors
  eval/      cross-validation, Hubble residuals
docs/tex/    LaTeX science docs (model-definition, salt2-distillation)
references/  source papers
```

Apart from `data/sample.py`, everything under `src/dgsn/` is empty scaffolding.

## Documentation

- **Scientific content → LaTeX** in `docs/tex/`. Always `\input{macros}`; cite
  from `refs.bib` (entries sourced from NASA ADS, never typed from memory).
- **Engineering content → Markdown**, including this file.
- **`docs/tex/model-definition.tex`** is a working note for the author, not a
  paper. Substance in mathematics, prose for an astronomer; no abstract; say each
  thing once. When a point needs a paragraph of argument, that is the signal it
  belongs in the backlog or should be restated as a formula.
- **`docs/model-note-backlog.md`** is a staging area, not a second note.
  Argumentative material parked so the note stays short. Material moves out of the
  note when it is argument rather than model, and back only when it is both settled
  and expressible as mathematics.

## Working agreements

- Claims in docs must be traceable to a citation or to code that produced the
  number. No asserted values.
- Report negative results plainly. A parameter that fails to earn its place is a
  result.
- When a commitment above is genuinely wrong, say so and change it explicitly.
