# DifferentialGeometry

Empirical modelling of Type Ia supernova multi-band light curves using
differential geometry. `g` and `r` only.

**The model is defined in `docs/tex/model-definition.tex`, which overrides this
file.** Here: the data, the traps, the tooling.

## Standing commitments

1. **Curve in magnitude space, likelihood in native flux.** Magnitude space
   because distance and reddening act there as translations, so `kappa` is
   invariant; flux because ZTF errors are Gaussian there and negative fluxes are
   meaningful. **Never convert data to magnitudes.**
2. **`z < 0.05`.** A two-band model has no SED and cannot self-consistently
   K-correct; the cut keeps the correction small enough to neglect.
3. **External SED templates and extinction laws quantify systematics only** —
   never to fit. Initialisation is exempt: it sets where chains start, not what
   they converge to, so SALT2 places the starting point (`model-definition.tex`
   §5) and dispersed-start agreement is what checks it.

Others were retired 2026-08-28 and are in the backlog. **SALT2 strictly
downstream** and **evaluation is out-of-sample** must be re-established before
any result is claimed.

## Data

ZTF SN Ia DR2 through [`ztfcosmo`](https://github.com/ZwickyTransientFacility/ztfcosmo),
which needs no credentials — it reads from `ztfcosmo.in2p3.fr` or `$ZTFCOSMODIR`.
Local copy unpacked at `data/raw/ztfsniadr2_lite/`; everything under `data/` is
gitignored.

**`tables/snia_data.csv`** — 3628 rows: `ztfname`, `redshift`, SALT2
`t0, x0, x1, c` with covariances, `mwebv`, `sn_type`, `lccoverage_flag`,
`fitquality_flag`.

**`lightcurves/<name>_lc.csv`** — whitespace-delimited under a `#` header:
`mjd filter flux flux_err ZP flag mag mag_err field_id rcid flux_offset
offset_unc err_scale in_baseline`. Filters `ztfg/ztfr/ztfi`; `ZP = 30`;
difference-imaging flux with negatives retained; `mag = 99` is a non-detection.

Selection and every count come from `python -m dgsn.data.sample` — regenerate
rather than quote from memory.

Milky Way extinction is **not** removed: the data reach the fit as observed and
the Galactic term is taken out downstream. `mwebv` is the SFD `E(B-V)`. Not a
small-extinction sample — median 0.042, 23% above 0.1, max 1.06.

## Traps

Every one of these has silently produced a wrong number.

**Quality flags are NaN for 18 objects, and `!= 0` admits them.** `NaN != 0` is
True in pandas, so `lccoverage_flag != 0` keeps objects whose quality is
unknown. Gave a wrong selection of 669 that stood for a week. Require `== 1`.

**Carry the selection through to the counting stage.** Counting epochs by
iterating light-curve files on disk rather than the selected objects applied the
`z < 0.05` cut when the file list was built but not again when epochs were
counted. One supernova at `z = 0.0925` passed every threshold and inflated all
six counts by one. Selection is an object list, not a directory listing.

**`sn_type != "snia-pec"` does not select normal SNe Ia.** One 91bg/86G object
is typed `snia`, with the peculiarity recorded only in `sub_type`. Select
positively on `sn_type == "snia-cosmo"`, which is DR2's own cosmology
classification and admits `norm`, `91t` and `99aa` and no peculiar subtype.

**`filter` is a DataFrame method.** `lc.filter == "ztfg"` compares a bound
method to a string, is False for every row, and yields zero epochs. Use
`lc["filter"]`.

**Epoch `flag` is a bitmask and `flag == 0` is the wrong cut.** The `ztfcosmo`
default excludes only bits 1, 2, 4, 8, 16 — **`flag & 31 == 0`** — being
`flux_err == 0`, `chi2dof > 3`, `cloudy > 1`, `infobits > 0`, `mag_lim < 19.3`.
Bits 32–1024 encode seeing, field, moon, airmass and are informational. Cutting
on `flag == 0` discards most genuine detections.

**Do not apply the 0.86 rescaling twice.** Schlafly & Finkbeiner (2011)
coefficients are defined against the **raw** SFD map and already carry the ~14%
recalibration, so feed them uncorrected `mwebv`. Rescaling the map *and* using
SF11 coefficients double-counts; the older SFD98 values over-correct. Compute
the coefficients for ZTF's own filters rather than borrowing SDSS-like ones.

## Stack

JAX · optax · numpyro (posteriors, later). No neural networks: the global
functions are basis expansions. Training is an auto-decoder — basis coefficients
and the per-SN latent array in one pytree, optimised jointly.

Two independent implementations of the geometry that must agree on `kappa`:
quadrature of the Frenet relations, and autodiff of a directly parameterised
`gamma(s)`.

`src/dgsn/` is empty scaffolding apart from `data/sample.py`.

## Conventions

- Science in LaTeX under `docs/tex/`; `\input{macros}`; cite from `refs.bib`,
  entries taken from ADS and never typed from memory.
- `model-definition.tex` is a working note, not a paper: substance in
  mathematics, no abstract, say each thing once. An argument needing a paragraph
  belongs in `docs/model-note-backlog.md`, which is staging and not a second note.
- No asserted values — claims trace to a citation or to code that produced them.
- Report negative results plainly. A parameter that fails to earn its place is a
  result.
