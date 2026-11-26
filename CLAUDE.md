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

## Traps in measuring, not in the data

These are about the harness rather than the archive. Each has produced a wrong
number here, and each was reported before being caught.

**Measure the floor first.** Run the pipeline with the exact answer substituted
for the thing being estimated — the true `kappa`, the true log rate — and see
what error remains. If it does not vanish, the harness is the limit and every
number downstream is noise about the harness. Four separate conclusions about
basis size were withdrawn for want of this one check, including "no `delta`
works at low `x1`", which was a phase grid too coarse to resolve the turn.

**Two resolutions, always.** One number is a claim; two that agree are a
measurement. A near-cusp defeats any uniform grid, and the grid that binds is
rarely the one last thought about — refining the phase grid says nothing if the
arclength grid is the limit. Report the resolution alongside the result so
under-resolution is visible rather than inferred.

**Verify the edit landed.** `str.replace` and `sed` fail silently on a mismatch,
and a patch applied to a wrapper's source does nothing when the constant lives
in the module it `exec`s. Twice, "refinements" reran the identical calculation.
Results identical to many digits across a supposed change are the tell: treat
them as evidence that nothing changed, not as evidence of convergence.

**`stepsize__` is in the adapted metric's units and is not comparable across
programs.** The metric differs between runs, so a larger `eps` can mean a better
geometry or a differently scaled metric and nothing else. Two claims were
withdrawn for this on 2026-11-25: "the comoving frame's step is 3.4 to 6.4x
larger", which compared a mid-warmup value to a finished one, and "7.5e-3 is 2.6
to 4.2x above every hier_rot chain", where `eps * n_leapfrog__` — the trajectory
length in each chain's own units, the quantity that is comparable — read 13.29
against hier_rot's 12.97. Compare `eps * n_leapfrog__`, or compare within one
program.

**A falling `n_leapfrog__` beside a rising `divergent__` is truncation, not
speed.** A divergence ends tree growth, so a failing chain does fewer gradients
per draw and reports more draws per hour. `hier_rot` chain 1 fell from 4095 to
980 leapfrogs while its divergence rate went 0.14 to 0.99, and `rot_white` chain
4 ran 1.5x faster per draw than its siblings at 72 percent divergent. Count
gradients per draw you can keep, not per draw.

**At `max_depth` saturation the leapfrog count carries no information.** When
`cap(12)` is near 1 every clean draw costs 4095 and the divergent ones cost
almost nothing, so gradients per clean draw is just 4095/(1 - div) and
chain-to-chain agreement in it is arithmetic, not a virtue: `hier_rot`'s 4097,
4097, 4096, 4130 was read as uniformity when it was saturation. It also means
the trajectory length the posterior actually wants is unmeasured and larger than
4095 — raising `max_depth` is the only way to learn it.

**`accept_stat__` is worthless once the metric is frozen.** After a divergence
the retained subtree is short and local and accepts well, so the statistic
improves as the chain fails. `hier_rot` chain 1's sampling phase reported mean
0.829 to 0.927 with tenth percentile 0.478 to 0.814 while 93 percent of its
draws were divergent. Read `divergent__`.

**`~` drops the normalising constant; `target +=` does not.** A predicted `lp`
difference was out by exactly 17.4598 = 3.5 log(2 pi) + 6 log(2 pi) because
`z_a ~ std_normal()` omits the constant while the hand calculation included it.
A constant offset equal across all chains is the signature.

**The warmup schedule decides which draws are comparable.** `init_buffer` 75,
`window` 25 doubling, `term_buffer` 50, so for `num_warmup = 500` the metric
updates at the end of 100, 150, 250 and 450 and is frozen from 451. Draws either
side of one of those boundaries are different samplers. `hier_rot` chain 1 broke
at exactly 451 — divergence 0.14 to 0.78, leapfrogs 4095 to 1320 — when the
metric froze and dual averaging restarted with 50 iterations to settle a step
size on a bimodal acceptance signal.

**Two runs of one model differ by 2x, so warmup-window comparisons cannot
resolve 1.5x.** `runs/white-k812-n60` and `runs/whiteR-k812-n60` are the same
model, same seed, same inits, differing only in the order of two floating-point
operations. They part company at draws 8 to 12 with `lp` differing in the
seventh figure and by draw 42 differ by 1.99x in gradients per clean draw. Four
chains and tens of draws is not a measurement of a parameterisation; the finished
run's effective sample size per gradient is. The Hessian of `lp` at stored draws
is the cheap and reproducible alternative — two finite-difference step sizes
agreeing to five figures, 2.7 percent coefficient of variation over sixteen
draws — and it predicted a reparameterisation's realised spectrum to four
figures where the sampling comparison could not see it at all.

## Stack

Stan, and only Stan, for the model: the forward map, the likelihood and the
sampler are one `functions` block under `stan/`, `#include`d by each program
(`model-definition.tex` §5). Python (`src/dgsn/`) builds data, drives
`cmdstanpy` and reads output; it evaluates no basis, quadrature or clamp.
`sncosmo` synthesises the SALT2 and Hsiao templates that place the start of
the fit and measure the floors. No neural networks: the global functions are
basis expansions.

`README.md` has the environment and the stage-by-stage run instructions.

## Conventions

- Science in LaTeX under `docs/tex/`; `\input{macros}`; cite from `refs.bib`,
  entries taken from ADS and never typed from memory.
- `model-definition.tex` is a working note, not a paper: substance in
  mathematics, no abstract, say each thing once. An argument needing a paragraph
  belongs in `docs/model-note-backlog.md`, which is staging and not a second note.
- No asserted values — claims trace to a citation or to code that produced them.
- Report negative results plainly. A parameter that fails to earn its place is a
  result.
