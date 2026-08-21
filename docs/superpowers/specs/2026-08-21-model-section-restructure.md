# Restructure `model-definition.tex` as a Model section, with the zero-flux boundary built in

**Date:** 2026-08-21
**Status:** Implemented. The note, `CLAUDE.md`, this directory's `2026-08-20` design
doc, and `docs/model-note-backlog.md` were all updated in the same pass.

This is the design record for that pass, kept because it holds material the outcome does
not: the alternative architecture that was weighed and rejected, why the choice still
stands, the honest cost of it, and what was ruled out of scope. The *conclusions* live in
`CLAUDE.md`; this file is the reasoning behind them.

**One target was missed.** The note was expected to come out shorter than 14 pages. It
came out at 14 (13 of body plus one of bibliography), unchanged from before: deleting the
ladder and limitations sections funded the new §2 exactly and no more.

## Context

The note reads as a whole article: nine sections in which the model is interleaved with the
arguments defending it. The author wants it to read as the **Model section** of an article —
the model stated as a hierarchy, everything supporting it demoted behind.

The hierarchy:

1. A light curve is a **curve** plus the **traversal** of that curve.
2. Curves are described by invariant functions such as **curvature**, represented by
   parameterised ML models.
3. The **orientation and position** of the curve are set by the physical fact that flux is
   zero before explosion and tends to zero after.
4. The **traversal** is a parameterised function, again an ML model.
5. That predicts **magnitude versus time per band**; with per-supernova **date of
   `g`-maximum** and **peak `g` magnitude** it predicts the observables, giving the
   likelihood.

Two further requirements make this more than a reshuffle:

- Light curves **include pre-explosion epochs whose expected flux is zero**, entering the
  likelihood.
- `kappa(s)` has a **built-in early zero-curvature segment**, by construction rather than
  learned.

**A correction to the tempting justification.** Pre-explosion epochs do *not* let the fit
"see" the early asymptote: their predicted flux is zero for every value of every parameter,
so they carry no Fisher information about the geometry. What they constrain is the
**turn-on epoch**. What legitimises item 3 is different and narrower — a built-in straight
segment puts the early ray *inside the model domain*, so the boundary condition is an
initial condition again rather than data at `s = -infinity`. The late
(radioactive-tail) asymptote stays outside and stays a diagnostic.

## Decisions taken (author, 2026-08-21)

**Orientation is anchored at the early ray.** Integrate the Frenet system from the straight
segment with `T = -(1,1)/sqrt(2)`. Consequences:

- The arclength origin moves to the **end of the zero-curvature segment**: domain
  `s in (-inf, s_end]` with `kappa(s;theta) = 0` for `s <= 0`.
- The `g`-maximum is located implicitly as the root `s_g(theta)` of `dm_g/ds = 0` — a
  `theta`-dependent root-find, differentiable via the implicit function theorem.
- **The frame sign becomes a prediction.** Previously `T(0) = (0,±1)` was an empirical
  config value; now it follows from `sign(kappa)`, so "which band peaks first" is
  *predicted* and checked rather than supplied. A strict improvement; decision 4(c) is
  rewritten.
- `t_max` remains the parameter and the `g`-max epoch, so the traversal must place
  `s(t_max) = s_g(theta)`. This is the accepted cost.

**The explosion epoch is derived, not fitted:**
`t_expl = t_max + w u_expl(theta) (1+z)`, so the per-SN count stays at seven. The built-in
segment is what supplies the fixed template location that makes this inversion possible —
requirement 2 is why requirement 1 costs zero parameters. Assumption: dark-phase duration
scales with `w`, i.e. rise time is a stretch property. Falsifiable; a free per-SN `t_expl`
is a candidate later rung, scored on held-out residuals like everything else.

### The strongest argument for this choice, which was not on the table when it was made

**It makes real progress on the `R^3` open question.** `dm_g/ds = 0` supplies **1** of the
`n(n-1)/2 = 3` parameters of `SO(3)`. Fixing `T` on the early segment to `-(1,1,1)/sqrt(3)`
supplies **2** of 3, leaving only the rotation about that axis. Decision 4(b) currently
records the `R^3` frame as "an **open question**, not settled"; this halves it. Record this
in the rewritten 4(b) — it is the reason the change is worth its cost, and it is a
*favourable* rewrite of the `n=3` remark, not merely a forced one.

### The dissent, and why the choice still stands

The review recommends a third architecture: keep the `g`-max anchor, keep `kappa = 0` early
as a template feature, and let the early tangent **come out** of the trained model, checked
against `-(1,1)/sqrt(2)` as a validation number. Its case is that this delivers the same
testable early-ray statement at zero cost — no revision to decisions 4(b), 4(c) or 5, no
reclassification of gauge as assumption, no new degeneracy, no root-find. Two of its three
objections are answered:

- *"The built-in segment is a finite straight piece at finite magnitude, not the
  asymptote."* Answered by result (i): the segment is **semi-infinite**, so `s -> -infinity`
  *is* in the closure of the domain. This is why (i) is load-bearing rather than a technical
  footnote.
- *"Pre-explosion data are direction-blind."* Conceded, and already recorded — but the
  direction is a boundary condition, not something the nulls were ever going to supply.
- *"It converts a gauge into an identifying assumption."* **Not** answered. This is the real
  cost and is recorded below.

Against it: Option 3 makes no progress on `R^3`, and leaves the frame sign a config input
rather than a prediction. The author's choice stands; the cost is paid explicitly.

### The honest cost of anchoring at the ray

Fixing the ray direction to `-(1,1)/sqrt(2)` is **not** a costless gauge choice, and the
note must not present it as one. A straight segment along the diagonal means
`d(m_g - m_r)/ds = 0`: **early colour is constant**, equivalently the two bands share an
early power-law index. That is a falsifiable astrophysical claim. Two mitigating facts to
state alongside it: the segment covers only the earliest phases, and ZTF photometry there is
sparse enough that the assumption is weakly tested by this data — which is both a defence
and a caution.

Contrast with the anchor it replaces: `T(0) = (0,±1)` at the `g`-maximum was genuinely
costless, since it only labelled `s = 0`. Swapping a costless condition for a substantive
one is the real content of this change, and decision 4(b) must say so.

### The gauge divides between the two anchors

| Gauge | Fixed by | Note |
|---|---|---|
| Orientation (`SO(2)`) | early ray, `T = -(1,1)/sqrt(2)` | inside the domain; substantive, not costless |
| Arclength origin | end of the straight segment | replaces the `g`-max anchor |
| Position, parallel `(1,1)` | `mu`, free per SN | absorbs distance and peculiar velocity exactly |
| Position, perpendicular `(1,1)` | early-phase colour | constant along the ray; a template quantity |

**A new exact degeneracy, and the reason the origin must be fixed where it is.** On a
straight segment, translating the arclength origin is a displacement *along the ray*, i.e.
along `(1,1)` — which is exactly what `mu` does. So the `s`-origin and `mu` span an exactly
flat likelihood direction. Anchoring `s = 0` at the end of the built-in segment is what
removes it, and the note must say so; left unfixed, the fit has a genuine null direction.
This is new content that the record does not cover.

**`mu` as peak `g` magnitude needs one extra step.** With `s = 0` no longer at the
`g`-maximum, `gamma_g(s_g(theta))` is `theta`-dependent, so `mu_peak = mu + gamma_g(0)` no
longer differs from `mu` by a global constant. Fix by normalising the template per `theta`:
translate the integrated curve so `m_g = 0` at its own `g`-maximum. That is a deterministic
function of `theta`, not a new parameter, and makes `mu` *exactly* the peak `g` magnitude.
Keep `mu` an all-band offset — a `g`-only offset would mix the distance and colour channels
and revise decision 1. Recommend the perpendicular offset be a template function of
`theta` rather than one global constant; global would assert every supernova shares the
same early colour, and `theta`-dependence costs nothing since it is template, not per-SN.

## The mathematics to state

All three verified. (i) is unit speed alone; (ii) checked numerically —
`0.4 ln(10)/sqrt(2) = 0.6512`, and fireball `alpha = 2` gives `A = 3.071`; (iii) derived
independently from `integral kappa ds = Delta phi`, which also corrected it to be
branch-dependent. Each is stated in the note with its own check attached.

**(i) `f = 0` is unreachable at finite arclength, so the traversal must be singular.**
`f = 0` means `m = +infinity`; the curve is unit-speed so
`||gamma(s) - gamma(0)|| <= |s|`, hence `|s| = infinity`. The cubic `eq:reparam` reaches
`s = -infinity` only at `t = -infinity`, so **exactly zero predicted flux is unreachable at
any finite epoch**, and a naive gate (`f = 0` before `t_expl`, geometry after) makes
predicted flux *discontinuous* at `t_expl`. The resolution consistent with the record is a
semi-infinite zero-curvature segment traversed by a map with a logarithmic singularity, so
`s -> -infinity` as `t -> t_expl^+`. This **replaces `eq:reparam`** and requires restating
its monotonicity clause.

**(ii) A straight ray plus a logarithmic map is a power-law rise.** On a ray of direction
`-(1,1)/sqrt(2)`, `m_X = c_X - s/sqrt(2)`, so `s = A ln(t - t_expl) + B` gives

```
f_X  ~  (t - t_expl)^alpha ,    alpha = 0.4 ln(10) A / sqrt(2)  ~=  0.651 A
```

— the standard early power law, index set by the coefficient of the log; fireball
`alpha = 2` gives `A ~= 3.07`. This is the payoff: the built-in segment is not an
architectural hack but reproduces the observed early rise, with `alpha` as a template
quantity.

**(iii) The sum rule splits at the `g`-maximum, and the split is branch-dependent.** With
`T = (cos phi, sin phi)`, `integral kappa ds = Delta phi` exactly. The early asymptote is
`phi = 225 deg`; `T(s_g) = (0,-1)` is `270 deg`, `T(s_g) = (0,+1)` is `90 deg`. So

| branch | pre-max | post-max | total |
|---|---|---|---|
| `kappa > 0` | `+pi/4` | `+3pi/4` | `+pi` |
| `kappa < 0` | `-3pi/4` | `-pi/4` | `-pi = pi (mod 2pi)` |

**The pre-maximum turning integral is not a branch-independent constant.** Two consequences.
First, `pi/4` may only be quoted together with the branch it belongs to; an unqualified
`pi/4` would be wrong for half the possible orderings. Second, this is exactly why the
integral must **not** be imposed: imposing `pi/4` for every `theta` would *assert* the peak
ordering, so no supernova could come out with the minority ordering and decision 4(c)'s
mixed-ordering diagnostic would return nothing by construction. Report it as a validation
number, stratified by the predicted `sign(kappa)`.

## Target structure

**Six** sections replacing nine; all argument behind the model.

| § | Content | Source |
|---|---|---|
| 1 | Motivation | trim current §1 |
| **2** | **The model** — the deliverable of this pass | below |
| 3 | Why magnitude space | current §2.3, argument only |
| 4 | Phase-dependent dust and its diagnostics | current §2.4, argument + diagnostics |
| 5 | Samples | current §4 |
| 6 | Relation to SALT2, and distance | current §5 incl. `eq:tripp` |

**The ladder and the limitations sections come out** (current §8, §9). Neither is model
definition, and the note is a working note about the model. Nothing is lost from the durable
record: the L0–L3 table is already in `CLAUDE.md`'s "The deliverable is a ladder" and the
validation method in the spec, and every limitations entry is already in `CLAUDE.md`'s "Known
limitations to keep in view". Two things this requires:

- **Resolve the seven inbound cross-references, do not leave them dangling** — `sec:ladder`
  (3) and `sec:limitations` (4). Each becomes either an inline clause or a pointer to
  `CLAUDE.md`; a `\ref` to a deleted label is a build failure, which the verification step
  catches.
- **Park the two entries that were new content, in the backlog**, since the note is no longer
  their home: the `n >= 3` `tau`-unidentifiability finding, and the re-scoped near-zero
  curvature item. Both also go into `CLAUDE.md` per the edit list below, which is where they
  are actually load-bearing.

This is also most of the length funding — the new §2 material is paid for outright, and the
consolidations below become headroom rather than a constraint. Expect the note to come out
**shorter** than 14 pages, not net-zero.

§2, one subsection per level:

- **2.1 A light curve is a curve and a traversal.** `gamma(p) = (m_1,...,m_n)`; the
  two-object split; the declaration that the curve lives in magnitude space with metric
  `diag(1,1)` — declaration only, argument cross-referenced to §3.
- **2.2 The curve.** Arclength, `eq:frenet`, `eq:frenet2d`, torsion for `n=3`,
  `eq:kappanet`, and the **built-in zero-curvature segment** with its frame remark: on
  `kappa = 0` the classical normal is undefined on an *interval*, so there is no continuous
  extension by limits; the frame is carried by the ODE, where `kappa = 0` gives `T, N`
  constant.
- **2.3 Where the curve sits: the zero-flux boundary.** Hairpin; early ray carrying
  orientation, with the constant-early-colour assumption stated as such; position parallel
  to `(1,1)` deferred to `mu`; perpendicular separation as terminal colour, an *output*;
  `eq:sumrule` and split (iii); the frame sign now predicted.
- **2.4 The traversal.** The singular map replacing `eq:reparam`, result (ii), `ds/dt > 0`,
  the deterministic `1+z` division, and `t_expl` derived.
- **2.5 From curve to observables.** `eq:forward`; the flux conversion carrying the
  divergence, with `f = 0` as its image; `eq:chi2` **with the rank-one term** below;
  the parameter table.
- **2.6 Gauge.** Short table: each condition, whether gauge / substantive / predicted, and
  where fixed.

### Parameter table, restructured

| Level | Parameters |
|---|---|
| 2.2 the curve | `theta_1, theta_2` |
| 2.3 placement | **none** — worth stating outright |
| 2.4 traversal | `w, a_1, a_2` |
| 2.5 observables | `t_max, mu` |
| L2c only | `c` |

Still seven, each attached to its layer. `t_max` mechanically enters the traversal even
though it is *interpreted* at 2.5; `c` is a deformation of the integrated curve, belongs to
no level, and surfaces only in `eq:forward`.

## Pre-explosion epochs: what the archive actually says

The phrase "expected flux is exactly zero" is **false at the archive's own data model**, and
this materially changes `eq:chi2`. Three findings:

- **`flux_offset ± offset_unc`** means the difference-imaging zero level is neither exactly
  zero nor exactly known. Expected pre-explosion flux is `flux_offset`, and zero only after
  it is subtracted.
- **`offset_unc` is rank-one** — one offset per light curve, common to every epoch of a
  band. Negligible for bright detections; **dominant** across a long run of pre-explosion
  epochs. In a diagonal `eq:chi2` their joint constraint appears to tighten as
  `1/sqrt(N_pre)`, whereas the true floor is `offset_unc` and does not shrink. So a diagonal
  treatment **overstates their information by a factor growing with `N_pre` and biases
  `t_expl` in the sign of `flux_offset`**. Fix: subtract `flux_offset` deterministically and
  marginalise `offset_unc` analytically as a known rank-one covariance (Sherman–Morrison) —
  the same rule already applied to `1+z` and `mwebv`: known effects are removed, not fitted.
  Keeps the count at seven.
- **`in_baseline` epochs must be excluded.** They are the data from which the zero level was
  estimated, so their residuals are shrunk toward zero by construction; using them as
  constraints on `f = 0` double-counts, giving an over-tight `t_expl` and an understated
  chi-squared.

**All three readings of these columns are inferred from the column list and must be verified
against `ztfcosmo` before being written anywhere** — "no asserted values" applies to the data
model as much as to results. This verification is a prerequisite of the data-side edits, not
of the note restructuring, so it does not gate §2.

Also: **drop the `SNR > 5` cut for pre-explosion epochs.** Expected SNR there is zero, so
requiring `SNR > 5` keeps only upward noise excursions — a one-sided selection on the noise
realisation of the very quantity being fitted, biasing `t_expl` early. Select these epochs on
**provenance and error only, never on measured flux**.

**And `flag & 31 == 0` is the wrong cut for nulls — a change in kind, not a tightening.**
The recorded warning runs the other way: it protects *detections*, since `flag == 0`
"preferentially keeps faint baseline epochs". For a null the failure mode inverts.
*Depth* is the benign problem — a shallow epoch has large `sigma` and self-weights down under
a correct Gaussian flux likelihood, so bit 16 costs information, not correctness. *Bias* is
the fatal one: the whole signal is "is this flux consistent with zero", and bits 32–1024
(seeing, field, moon, airmass), recorded as "informational — **not** excluded", are precisely
the bits correlated with systematic offsets. A small offset negligible beside a bright
detection is fatal beside a null. So the null epoch set needs its own cut, derived for that
purpose. Pin down whether `sigma` is pre- or post-`err_scale`, which matters far more at zero
signal than for detections.

**Three windows, named separately.** The note currently conflates one window serving as the
curve's domain, the likelihood's epoch set, and the selection criterion. Split: *model
support* extends before `-15` d; *likelihood epoch set* is in-window detections **plus**
pre-explosion epochs; *selection window* stays `[-15,+40]` d. **Pre-explosion epochs must
not count toward ">= 5 good epochs per band"** — they constrain no geometry, so counting
them would admit supernovae with no geometric constraint, and would invalidate the 599/177
counts and the regression test asserting 599.

Most of this is data-ingest engineering and belongs in `CLAUDE.md`'s "Data details" and the
spec, **not** in the LaTeX note. Only the three-window split and the non-diagonal covariance
are model statements.

## Length: the funding source

New material is paid for by consolidating duplications already located:

- `g`-max anchor stated four times (428–442, 651–663, 787–792, 835–841); "one condition
  fixes both gauges" three times. §7's recaps duplicate §6 almost entirely.
- Nuisance-effects table twice (110–125, 630–642).
- Flux-vs-magnitude rule four times (196–200, 340–345, 507–512, 767–769).
- `eq:phase` never referenced and contained in `eq:reparam` — delete.
- Parameter inventory in four places (459–473, 474–476, 580–594, 850–862).
- Dust-vs-shape asymmetry three times (308–313 + 346–354, 817–833, 901–909).
- `n=3` needs two more conditions; frame sign empirical; no-effect-acts-as-rotation — each
  stated twice or three times.

With §8 and §9 gone these consolidations are headroom rather than the budget, so aim for a
note **shorter than 14 pages**. No abstract, and **no opening list of the five levels** — that
is an abstract by another name and breaks "say each thing once".

## Text that becomes false and must be amended, not merely extended

- **Lines 762–769**: "The divergence of `m` at both ends is also why the phase window exists
  at all: flux space would keep the curve bounded ... but would forfeit `eq:transinv`." Once
  the likelihood includes epochs *at* the divergence, the window is no longer the answer to
  it; the answer is that the forward map `m -> f` absorbs the divergence with `f = 0` as its
  image.
- **Line 663**: "unlike the asymptotic constraint of Section~\ref{sec:sumrule}, it is
  observable inside the phase window" — the stated ground for preferring the `g`-max anchor,
  now reversed. Also the subsection title at **line 651**, "The frame is fixed by the
  $g$-maximum".
- **Lines 777–785**: the gauge/assumption inventory, whose line 782–783 asserts "the
  arclength origin, the metric, the latent normalisation and **the frame orientation are
  gauge**". Frame orientation must move to the assumption side. Likewise the "Frame
  orientation" paragraph at **lines 835–841**.
- **Lines 414–426**: `eq:reparam` and its "positivity is definitional, not regulatory"
  justification, superseded by the singular map.
- **Lines 428–442**: the chain-rule derivation that the fitted `t_max` *is* the `g`-max
  epoch, which assumed `s(t_max) = 0`.
- **Lines 911–916** (the vanishing-curvature limitation): the section it sits in is deleted,
  but the content splits **three ways** and must be routed, not dropped — and the current
  wording is wrong for `n=2` and understated for `n=3`. It is currently contingent — "may pass
  near such points" — whereas a built-in segment makes `kappa = 0` on a set of positive
  measure, for every `theta`, always.
  - **`n = 2`: harmless, and the recorded limitation does not bite at all.** `eq:frenet2d` is
    a rotation ODE, **nonsingular at `kappa = 0`**, because in the plane `N` is *defined* as
    the fixed 90-degree rotation of `T`, globally and independently of `gamma''`. Only the
    classical definition `N ∝ gamma''/|gamma''|` fails. The note already writes the frame in
    rotation form, so the primary sample needs one sentence saying so — into §2.2.
  - **`n >= 3`: a genuinely new limitation, worse than "ill-defined frame".** Where
    `kappa = 0` the osculating plane does not exist, so `tau` is not merely ill-conditioned
    but **undefined** — whatever `tau(s)` the network emits on the straight segment is
    **unidentifiable**. This compounds the open `R^3` frame question and requires a
    rotation-minimising (Bishop / parallel-transport) frame, not merely an "inflection-robust"
    one. New content: **`CLAUDE.md` does not record this.**
  - The *incidental* near-zero curvature at the secondary maximum remains a limitation, since
    only that case makes "report where `kappa` approaches zero" meaningful — but it now lands
    in `CLAUDE.md` and the backlog rather than in the note.

  Keep line 715's "distinct from the vanishing-`kappa` degeneracy" clause verbatim — it
  becomes more load-bearing, not less. Decision 4(d) is untouched: regularity is about
  *speed*, and a straight segment is perfectly regular.

## Cross-references that must survive

Referenced three or more times: `sec:phasedust` (9), `eq:transinv` (8), `eq:firstvar` (4),
`eq:orders` (4), `sec:samples` (4), `sec:rotation` (4), `eq:frameanchor` (4), `sec:gauge`
(4), `eq:frenet` (3), `eq:transcol` (3), `eq:bandext` (3), `eq:dustcurve` (3), `eq:dkappa`
(3), `eq:reparam` (3).

`sec:limitations` (4) and `sec:ladder` (3) are **deleted**, so their inbound references must
be rewritten inline — see the target-structure note above. `eq:frameanchor` (4) survives as a
label but changes meaning: it is no longer the frame condition, only the `g`-maximum
condition that locates `s_g(theta)`. Check each of its four call sites reads correctly under
the new meaning rather than assuming the reference still works.

Two couplings to fix rather than preserve:

- `eq:dustcurve` is **constitutive** (used by `eq:forward` and L2c) but sits inside the
  argumentative §2.4. Move the definition to §2.5, leave the argument in §4.
- All four references to `sec:rotation` come from MODEL blocks though §6 is mostly argument.
  Split: boundary condition and frame to §2.3, the no-rotation argument to §3.

Safe to delete (never referenced): `sec:motivation`, `eq:salt2`, `eq:pertdecomp`,
`sec:magspacemath`, `eq:phase`, `eq:chi2`, `sec:framesign`. The last two of those starred
labels would render wrongly under any `\ref` anyway.

## Consequences to record elsewhere

- **`CLAUDE.md`** — six specific edits, each an explicit revision, not a workaround:
  1. **Decision 4 heading**, "the frame is fixed by the g-maximum" → fixed on the early
     zero-curvature segment. State explicitly that **4(a) is untouched** — no rotation is
     still fitted — so the revision is not over-read.
  2. **Decision 4(b), essentially in full.** "One condition fixes both gauges, costs no
     generality, and is observable inside the phase window" — every clause becomes false.
     Rewrite the `n=3` remark *favourably*: 2 of 3 conditions, not 1 of 3.
  3. **Decision 4(c).** "carry it in config, do not hard-code" → the sign becomes a model
     output; the config entry disappears and the external measurement becomes a validation
     target.
  4. **Decision 5.** State that the asymptotic *tangent* is imposed at finite `s` (and why
     that is acceptable: result (i) puts `s = -infinity` in the closure) while the *integral*
     remains diagnostic-only and branch-dependent. Without this the extension is silent.
  5. **"Gauge must be fixed" bullet** under Known limitations: frame orientation moves out of
     the gauge list into assumptions; add the new `s`-origin/`mu` degeneracy and its fixing
     condition.
  6. Add the zero-curvature segment, pre-explosion epochs, derived `t_expl`, the
     constant-early-colour assumption, the `n>=3` `tau` unidentifiability, and the
     `flux_offset` / `offset_unc` / `in_baseline` / null-cut handling under "Data details".

  Note also: requiring pre-explosion coverage would define a **new subsample whose count must
  be measured from the archive, not asserted**. Out of scope here; do not quote a number.
- **Spec**: `frenet.py` handles a zero segment for free in `R^2` (the rotation ODE is regular
  there), but **`direct.py` must return exactly `0`** from what is a `0/0` form on a straight
  segment — a new requirement on the correctness oracle and a new regression test. Also:
  upgrade the "inflection-robust frame" risk row to a **Bishop / parallel-transport** frame
  requirement for `n >= 3`, and add an assertion that `tau` is reported as **unidentifiable**
  (not merely noisy) on the built-in segment. Add tests for the power-law rise (ii) and the
  branch-resolved turning split (iii). The line-270 regression test asserting 599 stands
  unchanged, because the selection window and its `SNR > 5` criterion are untouched.
- **Backlog**: nothing qualifies for *return* — every settled entry the new §2 needs is
  already back in mathematical form. Traffic runs the other way this time: park the two
  limitations entries that were new content (`tau` unidentifiable for `n >= 3`; incidental
  near-zero `kappa`) plus anything from the deleted §8/§9 that is not already in `CLAUDE.md`.
  Record for each whether it is settled or open, per the backlog's own convention. Resist
  re-narrating levels 2, 3 and 5 in prose; that is the failure mode that forced the 18→14 page
  cut.

## Verification

- `pdflatex` twice plus `bibtex`; zero errors, zero undefined references, zero undefined
  citations. Grep every deleted `\label` for surviving references first.
- Check the constitutive chain reads in dependency order with no forward references into
  supporting sections: `gamma` def → `eq:frenet`/`eq:frenet2d` → `eq:kappanet` → boundary
  conditions → traversal → `eq:dustcurve` → `eq:forward` → flux → `eq:chi2`.
- Verify (ii) by checking the `alpha`–`A` constant by hand in the note itself, and (iii) by
  confirming the halves sum to `pi`, so the note carries its own checks.
- Confirm no surviving `\ref` to `sec:ladder` or `sec:limitations`, and that each of the seven
  former call sites still reads as a complete sentence.
- Report before/after line and page count; confirm the note is shorter than 14 pages.
- No code to run — `src/` is empty scaffolding.

## Out of scope

Implementing any of it in `src/`; adding Hsiao/Fitzpatrick to `refs.bib` (ADS auth expired —
make no citation rather than fabricating a bibkey); revisiting `z < 0.05` or
SALT2-strictly-downstream.
