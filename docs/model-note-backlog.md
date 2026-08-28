# Model note backlog

A holding file, not a document. Material raised during the evolution of
[`docs/tex/model-definition.tex`](tex/model-definition.tex) that is worth keeping but
does not belong inline while the note is still moving. Each item records the claim, why
it was raised, and whether it is settled or open.

The note is a working document for the author, whose substance is mathematics and whose
shortness is a constraint — so material moves **here** when it turns out to be argument
rather than model, and moves **back** only when it is both settled and expressible as a
formula. Not everything here is destined to return; some of it exists so the note does
not have to carry it. See the documentation convention in `CLAUDE.md`, of which this
file is a deliberate exception: it is a staging area, not a scientific document.

---

## Dust and the signal are the same logical type

**Settled.** This is the correction that reorganised the whole magnitude-space argument,
and it is worth stating once at length because the wrong version is seductive.

There are two logical types of effect, not one axis:

- **Rigid displacements.** Luminosity distance, absolute luminosity, peculiar velocity,
  per-band calibration error. Each is a fixed property of the object, its distance, or
  the telescope, so it applies the *same* displacement to every point of the curve at
  every phase. Such an effect is a map of the plane the curve is drawn on — the plane
  slides, the drawing on it does not change — and `κ` is blind to it exactly, with no
  approximation and no residual.
- **Deformations.** Host and Galactic dust, *and intrinsic SN diversity, i.e. the signal
  itself.* Both displace different parts of the curve by different amounts, and that
  bends it. No choice of ambient space and no choice of invariant separates them: a
  deformation is a deformation.

**The criterion that must not be used** is "phase-independent versus phase-dependent".
It sorts dust onto the same side as the thing `κ(s;θ)` exists to measure, so it does not
distinguish nuisance from signal at all. Intrinsic diversity *is* a phase-dependent
physical effect.

**The criterion that must also not be used** is chromaticity. A zeropoint error differs
in every band, is strongly chromatic, and still leaves `κ` untouched. This is why the
per-band calibration row is worth keeping in the note's table — it kills the tempting
"κ survives distance because distance is achromatic" story.

Dust fails to be a map of the plane not because it is chromatic but because two
supernovae passing through the same point `(m_g, m_r)` at different phases suffer
*different* displacements: the absorption a broad filter delivers depends on the spectrum
inside that filter, and that spectrum evolves.

**Consequence for the programme.** Whatever separation between dust and shape is achieved
must come from modelling restrictions and population-level assumptions, which are
substantive and falsifiable — never from the choice of ambient space.

---

## What magnitude space actually buys, and what it does not

**Settled.** Not a nuisance-free shape description. Two concrete things:

1. The quantity actually being inferred — distance — is made *exactly* orthogonal to the
   shape description. This is what allows `μ` to be free per supernova and to absorb
   peculiar velocity exactly, and it is the reason the construction can be a distance
   indicator at all.
2. Dust is *reduced* from a transformation that destroys shape outright — in flux space
   an anisotropic dilation, which takes a circle to an ellipse — to a deformation about a
   rigid part, of which the geometry removes the rigid part exactly.

Reduced, not removed. Do not overclaim `κ` as nuisance-free.

**A corollary about what is not delivered.** Since `κ` ignores *all* additive offsets, it
carries no information about absolute magnitude or absolute colour, only about how they
evolve. Both must be supplied separately: absolute magnitude by `μ`, absolute colour by
the globally fixed initial point `γ(0)`, or by `c` in rung L2c. A limitation on what the
geometry can deliver, not a defect in it.

---

## The order counting, and its two caveats

**Settled as an expectation; the numbers are open.** Same *type* is not the same *size*,
and the sizes decide whether `θ` is usable.

Split `c·u(s) = c·ū + c·δu(s)` with `ū = ⟨u⟩` and `⟨δu⟩ = 0`. The first term is a constant
vector, hence a translation, hence invisible to `κ` exactly. Only the second bends the
curve, and it is small twice over: in the reddening amplitude `c`, and in the fractional
phase variation `‖δu‖/‖ū‖`, which is small because a band-integrated extinction responds
only weakly to SED evolution inside the filter — the effective wavelength shifts by a
small fraction of the filter width, and `A_X` shifts by a correspondingly small fraction
of itself.

So there is a **hierarchy**: intrinsic diversity enters `κ` at first order (it *is* the
signal), dust at second. Expect the phase-dependent variance the geometry sees to be
dominated by intrinsic SN dispersion. This is a quantitative claim to be measured, not a
structural one to be asserted.

**Caveat 1 — amplitude dependence.** `c` is not small for the reddest objects, so the
hierarchy degrades where extinction is largest. Every diagnostic must be stratified *as a
function of* fitted `c`, never sample-averaged: an average hides exactly the regime where
the argument fails.

**Caveat 2 — shape channel only.** In the *colour* channel dust and intrinsic colour mix
at **first** order, since both shift `ū`. That degeneracy is not second order, is not
resolved here, and is the same one SALT2 confronts through `β`. Do not let a clean
shape-channel result be reported as though dust were solved.

---

## Why mixing with the signal is tolerable for the distance indicator

**Settled.** `θ₁, θ₂` will carry some dust — a second-order amount by the order counting,
but nonzero. This does *not* threaten the deliverable, and would not even if the order
counting failed.

The deliverable is a distance indicator scored on held-out Hubble residual scatter, and
for that purpose isolating physical causes is not required: if `θ₁` partly encodes
reddening, the fitted Tripp coefficients `α, β` absorb the mixture, and the residual
scatter — which is what is being measured — is unaffected by how the two contributions are
labelled. SALT2 is in exactly this position and always has been: its `c` mixes dust with
intrinsic colour, and `β` is fitted from the data rather than set to the `R_V + 1` of a
pure extinction law.

**The condition required** is only that the mixture be *stable across the sample*, since
a mixture varying systematically with redshift would propagate into cosmology. Here the
`z < 0.05` cut earns its keep for a second, independent reason: the sample spans too
little redshift for such a trend to exist, so there is no route from this contamination to
a cosmological bias. **That protection is specific to this analysis** and must be
revisited before the model is applied over a wide redshift range.

Contamination threatens *interpretation* of `θ`, not the distance indicator. Report the
two conclusions separately. Do not over-escalate this into "the interpretation programme
fails".

---

## The identifying assumption `Corr(c, θ) = 0`

**Substantive and falsifiable; status open until measured.**

Interpreting `θ` as explosion physics is a stronger goal than the deliverable. The
supporting assumption is physical: line-of-sight dust is a property of the intervening
ISM and has no causal connection to the explosion, so the two should be uncorrelated over
the population. Imposing zero sample correlation between fitted `c` and fitted `θ` fixes
the split.

This is an **identifying assumption** — substantive, falsifiable, possibly wrong — and is
of a different kind from the latent zero-mean/unit-covariance normalisation, which is a
costless gauge. It is known to be only approximate: the host-galaxy mass step couples SN
properties to host environment, and dust column is an environmental property, so a
genuine correlation is expected at some level.

**Whether it is load-bearing depends on the order counting.** If the hierarchy holds, the
shape channel is determined by the data at first order and the condition functions as a
*check* — `Corr(c, θ)` should come out small whether or not it is imposed, and finding
that it does is evidence for the whole picture. If the hierarchy fails, the condition
becomes load-bearing and conclusions about `θ` rest on an assumption rather than on the
photometry. Which obtains is measured, not chosen — so **impose it and also report the
unimposed value.**

---

## Why L2c fits an amplitude of `u(s)` rather than a constant vector

**Settled, with an honest downgrade.** A constant translation removes `c·ū` and leaves
`c·δu(s)` free to be absorbed by `θ₁, θ₂` — exactly the contamination L2c exists to
prevent, reached by a remedy that does not extend to it. So L2c is defined with the
phase-dependent displacement: `u` is precomputed once from an extinction law and a
template SED and held fixed, while its amplitude `c` is free per supernova. Same
parameter count.

**But by the order counting this is cheap insurance, not a necessity.** A constant vector
would capture most of the reddening displacement, and the refinement matters only for the
reddest objects, at the level the `Δκ` diagnostic will quantify. It is adopted because it
costs nothing, not because the fit would fail without it.

**What it costs.** The clean statement that `κ` is intrinsic under dust. Reddening becomes
a known deformation fitted in the forward model rather than a symmetry quotiented out by
the geometry.

**A second-order approximation inside the second-order term.** `u` depends on the SED,
which depends on the shape being measured, so the dust deformation and intrinsic colour
evolution are physically *coupled*, not merely similar in direction. Holding `u` fixed is
therefore itself an approximation — of the same second-order size, since it is an error in
an already-second-order term.

---

## SALT2 is structurally better placed on dust

**Settled; an honest cost of the two-band design.** SALT2's colour *law* `CL(λ)` is
phase-independent by construction, but SALT2 carries an SED, so band extinction `A_X(p)`
comes out phase-dependent automatically. This model has no SED and must import `u(s)` from
an external template. The comparison runs against the geometric model here and should be
conceded plainly rather than buried.

The mitigating point: SALT2's `c` mixes dust with intrinsic colour, so the first-order
colour degeneracy is a shared problem rather than a new defect of this model.

---

## The redshift cut does not help with dust

**Settled.** Phase-dependent extinction is a *rest-frame* effect: exactly as large at
`z = 0` as at `z = 0.1`. Phase-dependent extinction and K-corrections are the same
underlying problem — band-integrated magnitudes depend on an SED the model does not have —
but only the K-correction branch is suppressed by low redshift. The dust branch survives
the cut untouched.

The same caveat applies to the Milky Way correction: `mwebv` is known, but converting it
to `A_g(p)` and `A_r(p)` is not possible without an SED, so subtracting a constant offset
is an approximation there too — smaller only because Galactic columns at high latitude are
low.

---

## Parked from the deleted ladder and limitations sections (2026-08-21)

**Housekeeping, not argument.** The note was restructured as a Model section and its
ladder and limitations sections were removed — neither is model definition. Nothing was
lost from the durable record: the L0–L3 table lives in `CLAUDE.md`'s "The deliverable is
a ladder", the cross-validation method in the spec, and every limitations entry in
`CLAUDE.md`'s "Known limitations to keep in view". Two entries below were *new* content
with no other home at the time, and are recorded here as well as in `CLAUDE.md`.

*Superseded 2026-08-23.* `CLAUDE.md` no longer holds either section — it was pruned of
modelling choices, and both moved **here**, to "Moved out of `CLAUDE.md`" at the end of
this file. The cross-references in the paragraph above are stale and are kept only to
show where the material used to sit.

**Resist re-narrating the model in prose here.** Levels 2, 3 and 5 of the hierarchy are
now stated mathematically in the note; restating them argumentatively is the failure mode
that forced the 18→14 page cut.

### `tau` is undefined, not merely ill-conditioned, wherever `kappa_1 = 0`

**Settled.** *Rewritten 2026-08-22*: the built-in segment is gone (see the turn-on rung
entry below), so `kappa = 0` is no longer a set of positive measure by construction and
the claim is once again about **incidental** zeros. The mathematics is unchanged and does
not depend on how the zero arises.

For `n = 2` this costs **nothing**, and the recorded "curvature degenerates at
inflections" limitation does not bite at all: the Frenet system in the plane is a
rotation ODE, nonsingular at `kappa = 0`, because `N` is *defined* as the fixed
90-degree rotation of `T`, globally and independently of `gamma''`. Only the classical
definition `N ∝ gamma''/|gamma''|` fails. Integrating the ODE rather than
differentiating the curve is therefore the definition that survives, not a numerical
convenience.

For `n >= 3` it is a genuine limitation. Where `kappa_1 = 0` the osculating plane does not
exist, so `tau` is not ill-conditioned but **undefined**: whatever `tau(s)` the network
emits there is **unidentifiable** and must be reported as such. This requires a
rotation-minimising (Bishop / parallel-transport) frame, not merely an
"inflection-robust" one. It no longer compounds the `R^3` frame question, which is
**closed** — see the orientation/origin redundancy entry below.

### Incidental near-zero curvature remains a limitation

**Open.** Since 2026-08-22 this is the *only* case, the by-design segment having been
dropped: the `(m_g, m_r)` path may approach
`kappa = 0` *incidentally* near the secondary maximum of redder bands. Only this case
makes the instruction "report where along `s` curvature approaches zero" meaningful, so
the instruction survives — but it now belongs to `CLAUDE.md` and the code, not to the
note. Whether it actually occurs in the fitted window is an empirical question that no
number yet answers.

### Regularity is about speed, and is untouched

**Settled.** Decision 4(d) — coincident band maxima give `||dgamma/dp|| = 0`, hence
`ds/dt = 0`, and unit-speed parameterisation fails — is unaffected by any of the above.
A straight segment is perfectly regular: it has no turning, but it has unit speed. The
note's clause distinguishing the regularity condition from the vanishing-`kappa`
degeneracy became *more* load-bearing under the restructure, not less, and was kept.

---

## Why the latent count is `K = 2`

**Superseded 2026-08-23 — there is one shape latent, not two.** The note carries a single
`x_kappa`, and that is correct. The argument below is kept because its *structure* survives
and was reused: SALT2 parity plus a single increment. What changed is where the second
increment went. It is not a second shape latent but `x_p` in the traversal, so the model
still carries two degrees of freedom beyond SALT2 parity while `kappa` is conditioned by
one parameter alone. Read what follows as the reasoning that led there, not as the count.

*Original entry follows.* **Settled as a working value; the empirical question is what the
ladder answers.** Moved out of the note on 2026-08-21 — it is a justification of a
modelling choice, i.e. argument, not model definition.

The latent count is not fixed by the geometry: `n` fixes the number of curvature functions,
`n-1`, and says nothing about how many parameters condition them. So `K` is a modelling
choice, and ultimately an empirical one.

The working value `K = 2` is set by parity with SALT2 plus a single increment. SALT2
describes a supernova with one time-dependent shape coefficient `x_1` and one
time-independent colour `c`. In the present framework the colour is not a latent at all: it
is a translation, and translations leave `κ` untouched (`eq:transinv`), so it conditions
nothing. That leaves one latent to carry the SALT2 shape freedom and one more — the extra
time-dependent degree of freedom this model exists to test. Hence `θ = (θ₁, θ₂)`, and
whether the second earns its place is answered by the ladder, not asserted.

---

## Use of an external SED template

**Settled.** An external SED template (Hsiao et al. 2007) and a standard extinction law
(Fitzpatrick 1999) are allowed for *quantifying these systematics only* — constructing
`u(s)`, evaluating `Δκ`. Never to fit, initialise, or K-correct the photometry. Used this
way it does not compromise the independence from SALT2.

**Closed 2026-08-21.** `refs.bib` now carries `Hsiao2007` (`2007ApJ...663.1187H`) and
`Fitzpatrick1999` (`1999PASP..111...63F`), both exported from NASA ADS once the token was
renewed, and both cited in the note at the point where `u(s)` is defined. Volume and page
were filled in from the bibcodes, which encode them — ADS's short export omits both, as
the header of `refs.bib` already warns.

---

## The turn-on rung: singular traversal, `t_expl`, pre-explosion nulls, power-law rise

**Open — a candidate later rung.** Added to the model as `CLAUDE.md` decision 8 on
2026-08-21 and **struck in full on 2026-08-22**, one day later. Parked here because the
physics is right and the machinery is coherent; what failed was its place in the *current*
model.

**Why it was struck.** Three things, of which the first is fatal on its own terms.

1. **The defining equation was inconsistent.** `s in (-inf, s_end]` with `kappa_i = 0` for
   `s <= 0` truncates the *late* end, while the unit-speed argument that motivated the
   segment ("both ends escape to infinity") forces `(-inf, +inf)`. The same paragraph both
   called the one-sided statement the two-sided one and claimed the terminal colour was an
   *output* while placing the late asymptote outside the domain.
2. **It did not fix the gauge it was carried for.** Nothing forbade `kappa == 0` on
   `(0, s_1]` as well, so "the end of the segment" was not a well-defined arclength origin,
   and the origin-versus-translation flat direction it existed to remove survived it. The
   replacement — anchoring `s = 0` at the `g`-maximum, an **intrinsic** feature — does fix it.
3. **The dark phase never needed a semi-infinite ray.** By unit speed,
   `f = 0 <=> m = +inf <=> |s| = infinity`, so the dark phase is the **single ideal point**
   `s = -infinity` in the closure of the domain, not an interval. The semi-infinite ray was
   occupied by the *early rise*, which is finite in time and infinite in arclength — a fact
   about the traversal, not about the curve. Checked: on a ray with band indices
   `alpha_g, alpha_r`, `ds/d(Delta t) = (2.5/ln10) sqrt(alpha_g^2 + alpha_r^2)/Delta t`, so
   `s = A ln(Delta t) + const` with exactly the `A` of the log term the segment introduced —
   about `A ln10 ≈ 7` arclength units per decade of `Delta t`. The identification is right;
   it just does not require prescribing `kappa = 0` anywhere.

**What the rung would be.** A semi-infinite arclength domain with `kappa_i = 0` imposed
early; a traversal `s(t) = ... + A ln[(u - u_expl)/(-u_expl)]` singular as `t -> t_expl+`,
so predicted flux is continuous at turn-on rather than gated; `t_expl` either derived as
`t_max + w u_expl (1+z)` (the assumption being that dark-phase duration scales with `w`,
i.e. rise time is a stretch property) or free per SN; pre-explosion epochs in the
likelihood, which carry **zero** Fisher information about the geometry and constrain
`t_expl` only, so they must never count toward the ">=5 good epochs per band" selection;
and the payoff, an early power-law rise `f_X ~ (t - t_expl)^alpha_X` with
`alpha_X = 0.4 ln10 e_X A` — `0.651 A` on the diagonal ray, so fireball `alpha = 2` needs
`A ≈ 3.07`, and the ratio `alpha_g/alpha_r` is the observable form of the ray direction.

**What returning it would cost.** The data-side work is nontrivial and is recorded in
`CLAUDE.md`'s "Pre-explosion (null) epochs" block, still marked **not verified against
`ztfcosmo`**: `flux_offset` subtraction, rank-one `offset_unc`, excluding `in_baseline`
epochs, and a null-specific quality cut (the detection cut `flag & 31 == 0` is the *wrong*
cut for a null, and an `SNR > 5` cut on a null is a one-sided selection on the noise). It
would also define a further subsample, whose size must be **measured**, not asserted.
Scored on held-out residuals like any other rung.

---

## The orientation/origin redundancy, and why `n = 2` has no free orientation

**Settled, 2026-08-22.** This is the result that closed the long-standing "open `R^3` frame
question", and it did so by turning it into a count rather than answering it.

Integrating `kappa(.;theta)` from initial tangent angle `psi` gives the **same curve** as
integrating `kappa(.+a;theta)` from `psi + int_0^a kappa`. So the network's own arclength
origin and the frame orientation are not independent: they carry **one joint redundancy**.
Of the `n(n-1)/2` rotational numbers in the Frenet initial data, one is absorbed by sliding
the network origin, leaving

```
n(n-1)/2 - 1   global orientation parameters
```

— **zero** for `n = 2`, **two** for `n = 3`.

Three consequences. (i) The global orientation `Phi`, introduced earlier the same day as a
sample-wide fit parameter, is not added capacity but a *coordinate* on this redundancy; for
two bands it does not exist, and orientation is **pure gauge**. (ii) The old `R^3` question
("fixing a direction supplies only 2 of 3 — what fixes the third?") is malformed: only two
are ever free. (iii) The arclength origin must therefore be anchored at an **intrinsic**
feature of the curve. The choice is the principal `g`-maximum, `dm_g/ds|_{s=0} = 0` with
`T(0) = (0,-1)` and `kappa(0) > 0`. It costs no generality — every turning planar curve has
a vertical tangent — and it does **not** constrain which band peaks first. Sign check, worth
keeping because it is easy to get backwards: `m_g'' = -sin(phi) kappa = kappa(0) > 0`, a
**minimum of magnitude**, i.e. peak brightness. The mirror convention `T(0) = (0,+1)`,
`kappa(0) < 0` is equivalent. A useful side effect: `t_max` maps to `s = 0`, so the
`s_g(theta)` root-find disappears from the model entirely.

---

## Retired from `CLAUDE.md` as standing rules (2026-08-28)

`CLAUDE.md` was pared down because it had accumulated more rules than an exploratory
phase can carry. The three below were removed **as enforced commitments**, not because
they are wrong. Each is recorded here so the reasoning is recoverable and so re-adopting
one is a deliberate act.

### SALT2 strictly downstream

Was: *an independent benchmark on the same objects, never used to preprocess,
K-correct, interpolate, or initialise. Using it upstream would make the independence
claim false.*

This is what backed the project's independence claim, and it is the one whose removal has
teeth: nothing now prevents SALT2 entering upstream by accident. The known breach is
already on file — the epoch-counting window is placed with the archive's SALT2 `t0`,
because the window needs a date of maximum and the date of maximum is what the model
fits. Selection only; no photometry is transformed or initialised by it.

**Consequence of retiring it:** independence is no longer protected by default and has to
be re-established deliberately before any result is claimed. The check that settles the
window breach is a recount against an SED-free anchor, comparing the *symmetric
difference* of the two samples rather than the counts, since equal counts can hide
different membership.

### Evaluation is out-of-sample, always

Was: *with per-SN parameters in the high single digits and a conditioned network, good
in-sample fits are guaranteed and carry no evidential weight.*

The argument is unchanged and still correct; it is retired only as a standing rule, on
the grounds that during exploration in-sample fits are useful diagnostics even though
they are not evidence. Nothing about the reasoning expires — restore it before any
result is reported.

### Code is written for general `R^n`

Was justified entirely by the `R^3` torsion subsample, which was dropped on 2026-08-28
when the analysis was restricted to `g` and `r`. With no `R^3` sample the justification
is gone, so the commitment went with it. The `n = 2` case wants its own code path
regardless: the signed-curvature convention (normal as the rotated tangent) is specific
to the plane and is not a specialisation of the general-`n` Frenet apparatus, where
`kappa >= 0` is forced. See [the orientation/origin
redundancy](#the-orientationorigin-redundancy-and-why-n--2-has-no-free-orientation).

---

## Null epochs — parked, not currently used

Moved out of `CLAUDE.md` on 2026-08-28. No pre-explosion epochs enter the likelihood as
the model stands, and everything below rests on readings of the archive columns that are
**inferred from the column list and NOT YET VERIFIED against `ztfcosmo`**.


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

---

## The basis and anchoring redundancies, and the conventions that remove them

**Settled, 2026-08-27.** Moved out of the note, where it had been the subsection
`Degeneracies`. The note now carries only the closing paragraph of the parameterisation
section: the counting table and the statement that no exactly flat direction survives.
Everything below is the derivation behind that table — argument rather than model, and by
2026-08-27 a second telling of conventions the model section already adopts inline
(`eq:travgauge`, `eq:curvgauge`). It is kept because the count is the kind of claim that
gets challenged, and because two of its entries record a choice that could have gone the
other way.

**The general count.** Each expansion is linear in global functions with per-supernova
coefficients, and the model sees only the span, never the basis. A change of basis undone
by a compensating change of coefficients leaves every light curve unaltered. An expansion
of rank `n` with `k` coefficients pinned carries

```
n^2 - kn   flat directions
```

— with nothing pinned the whole of `GL(n,R)` is absorbed; pinning a coefficient forbids the
`n` transformations that would move it.

**Two routes, and which is available is forced rather than chosen.** The `n^2 - kn`
conditions come either from *declaring the basis outright*, possible only when its elements
are known a priori, or from *normalising the population of coefficients*, which is the only
route when the basis functions are learned. Placement takes the first and cannot take the
second, because `<mu>` is not free to normalise — it carries the cosmology. Stretch and
curvature deviation have no absolute scale, so normalising them costs nothing.

**Placement.** Since `1` and `e_c` span `R^2`, replacing `e_c` by `alpha*1 + beta*e_c` is
undone exactly by `(mu, c) -> (mu + alpha*c, beta*c)`. Any direction not parallel to `1`
gives the same model, so `e_c = (0,-1)` is declared rather than fitted. Both basis vectors
declared is four conditions, the whole redundancy. See also
[`e_c` is exactly unidentifiable below L2c](#e_c-is-exactly-unidentifiable-below-l2c).

**Curvature.** For any `eta` and `lambda != 0`,

```
(kappa_0, kappa_1, x_kappa) -> (kappa_0 + eta*kappa_1, lambda*kappa_1,
                                (x_kappa - eta)/lambda)
```

leaves `kappa` unchanged for every supernova. Only these two directions arise, not four,
because the coefficient of `kappa_0` is pinned at unity: rescaling `kappa_0`, or adding
`kappa_0` to `kappa_1`, would move it off unity and so is not available. Removed by
`eq:curvgauge` — `<x_kappa> = 0` fixes `eta`, `Var(x_kappa) = 1` fixes `|lambda|`, and the
sign condition `Cov(x, x_kappa) > 0` removes `lambda < 0`, there being no other way to
distinguish `kappa_1` from `-kappa_1`.

**Traversal.** Here the redundancy is the full four, because *both* coefficients of
`p(s) = x*p_0(s) + x_p*p_1(s)` are free per supernova. The four directions act linearly:

```
p_0 -> nu*p_0,          x   -> x/nu
p_1 -> lambda*p_1,      x_p -> x_p/lambda
p_0 -> p_0 + eta*p_1,   x_p -> x_p - eta*x
p_1 -> p_1 + rho*p_0,   x   -> x - rho*x_p
```

`eq:travgauge` fixes exactly one each of `nu, eta, lambda, rho`, the last by
`rho = Cov(x, x_p)/Var(x_p)`. In the singular limit `p_1 ∝ p_0` the family collapses to one
dimension and no convention helps.

**The one live choice: `Cov(x, x_p) = 0` is not Fisher orthogonality.** `Cov(x, x_p) = 0`
decorrelates the fitted coefficients *across the sample*, which is what makes `x_p`
interpretable as timing variation not already carried by stretch. Requiring `p_1` to be
orthogonal to `p_0` under the Fisher weight would instead decorrelate the *estimates* of
`x` and `x_p` within one supernova. The two fix `rho` to different values and cannot both
hold. The population condition is the gauge; the Fisher condition is a matter of
conditioning, to be handled by the basis used for fitting.

**Anchoring — not basis redundancies.** These concern where `s` and `p` are pinned, and are
recorded because a convention mistaken for a measurement is the error these entries exist to
prevent.

- *Arclength origin against frame orientation.* Integrating `kappa(.+a)` from
  `psi_0 + int_0^a kappa` reproduces the same curve, so origin and orientation carry one
  joint redundancy. Fixing `psi_0 = -pi/2` removes it: a shift now requires
  `int_0^a kappa = 0`, so generically `a = 0`. The anchor is self-enforcing, since
  `psi(0) = -pi/2` forces `m_g'(0) = 0`. The count behind this — `n(n-1)/2 - 1` free
  orientation parameters, zero for `n = 2` — is in
  [The orientation/origin redundancy](#the-orientationorigin-redundancy-and-why-n--2-has-no-free-orientation).
- *Phase origin against the traversal.* A per-supernova shift in `p` would be absorbed by
  the traversal were `p(0)` free, leaving any fitted epoch of maximum undetermined.
  `p_0(0) = p_1(0) = 0` removes it. Equivalently, it keeps the constant function out of
  `span{p_0, p_1}`: were the constant admitted, its coefficient would be exactly such a
  shift. This is a statement about the *dimension* of the span, not about a basis within it.

---

## `e_c` is exactly unidentifiable below L2c

**Settled, 2026-08-22.** The reddening direction was made a sample-wide fit parameter
earlier the same day; the same afternoon showed that below L2c there is nothing to fit.

With `mu` free per SN along `(1,1)` and `c` free per SN along `e_c`, the two span `R^2` for
any `e_c` not parallel to `(1,1)`. So changing `e_c` is undone **exactly** by an invertible
relabelling of `(mu, c)`: a flat likelihood direction, not a slow one. At L0–L2, `e_c` is
therefore **gauge** — fix it by fiat and do not fit it; a round-trip that appears
to recover it is a bug.

*Amended 2026-08-23.* Which value it is fixed at is **not** recorded here or in `CLAUDE.md`,
and the former instruction to fix it at `(1,-1)/sqrt2` is withdrawn. The choice is exactly
gauge for the fit but **not** for the meaning of `mu` and `c`: with `gamma(0) = 0` the
predicted peak magnitudes are `m_g(0) = mu + c e_g`, `m_r(0) = mu + c e_r`, so "`mu` is the
peak `g` magnitude" holds only when `e_g = 0` and "`c` is the peak colour" only when
`e_g - e_r = 1`, i.e. both only at `e_c = (0,-1)`; while `(1,-1)/sqrt2` instead makes `mu`
the peak *mean* magnitude and the colour shift `c*sqrt2`. Because the value and the prose
interpreting it must agree, both now live together in the note. Pinning the value in two
files is how the two drifted apart.

*Amended again 2026-08-23.* "Gains content only at L2c" is **too narrow** — it is a
statement about the per-SN light-curve likelihood alone. There are two channels, and only
the first is flat.

- **Per-SN light-curve fit.** Exactly flat in the direction, as above. Nothing to fit, at
  any rung.
- **Population standardisation.** Dust moves every reddened object along a *common* line in
  the `(mu, c)` plane: `(A_g, A_r) = A_g*1 + E(g-r)*e_c`, giving slope
  `dmu/dc = R_g = A_g/E(g-r)`. One supernova cannot determine that slope; the sample can.
  This is the two-band counterpart of SALT2 fitting `beta` rather than adopting `R_V`, and
  it works at **every rung, L0 included**, with no phase dependence required.

So the reddening direction is a fitted quantity throughout — fitted at the Hubble-diagram
stage, not per SN. The `du(s)` route of L2c is the harder channel and is no longer the only
one; what L2c adds is sensitivity to the *phase variation* of the displacement, not the
existence of a determinable direction.

Two things this does not change. The **basis** used to parameterise per-SN placement stays
gauge — `(1,1)` with `(0,-1)` is a naming convention, and the withdrawal of the instruction
to fix `e_c` at `(1,-1)/sqrt2` stands. And what the fitted slope measures is the
**effective** mixture of dust and intrinsic colour, exactly as SALT2's `beta` does, not
`R_g` for dust alone — the first-order dust/intrinsic-colour degeneracy is untouched by any
of this.

It gains content only at **L2c**, through the phase-varying `du(s)`, which is not a
translation and so cannot be absorbed by any relabelling of two amplitudes. There it is a
genuine global parameter, the two-band analogue of SALT2's fitted `beta`, with
`(1,-1)/sqrt2` the nested fixed-`R_V` idealisation.

The general lesson is worth stating once, because it recurs: **a direction is only
identifiable if something in the model varies along it in a way the free amplitudes cannot
mimic.** Two free per-SN amplitudes spanning the plane make every phase-*independent*
direction question vacuous.

---

## Moved out of `CLAUDE.md` (2026-08-23)

**Housekeeping.** `CLAUDE.md` was pruned of modelling choices. It now carries only what is
durable and operational — data handling, sample definitions, stack, layout, documentation
convention, working agreements, and a short list of standing methodological commitments —
because the model itself is a work in progress and the note is where it is being worked
out. Two files asserting the model is how they drift apart; `e_c` above is the worked
example.

**Everything in this section is OPEN, not settled**, and one conflict has since been
resolved against it.

*Resolved 2026-08-23: there is one shape coefficient, not two.* The note carries a single
`x_kappa` in a linear expansion `kappa = kappa_0 + x_kappa*kappa_1`, and that is correct;
the two conditioning latents `theta_1, theta_2` are withdrawn, along with the `K = 2`
derivation recorded earlier in this file. What the model gained instead is a second *new*
per-SN parameter of a different kind — `x_p`, in the traversal
`p(s; x_s, x_p) = x_s(p_0 + x_p p_1)` — so the count of new degrees of freedom beyond
SALT2 parity is unchanged at two, but one of them acts on shape and the other on timing,
where previously both acted on shape.

The ladder below has not been rewritten to match, and **L1/L2 as stated no longer have
referents**. The obvious reading is L0 = `(mu, c, x_s, t_max)`, still SALT2 parity, then
one rung per new parameter, `+x_kappa` and `+x_p`, in whichever order is to be argued.
That reading is *not* adopted here — it is a modelling decision for the note. Treat the
table as history until it is made.

### The ladder

| Rung | Free per SN | # | Question |
|---|---|---|---|
| L0 | `mu, c, w, t_max` | 4 | **same count as SALT2** — does the mechanism alone win? |
| L1 | `+ theta_1` | 5 | does one shape parameter earn its place? |
| L2 | `+ theta_2` | 6 | does a second? |
| L2c | L2, `c` upgraded to full `u(s)` | 6 | does phase-dependent dust matter? |
| L3 | `+ a_1, a_2` | 8 | does nonlinear timing earn its place? |

**L0 is a controlled experiment, not a starting point.** It carries four parameters in
one-to-one correspondence with SALT2's `(x0, x1, c, t0)` and differs in *mechanism* alone:
SALT2 produces stretch-like variation from an additive component, L0 from an exact
reparameterisation. L0 versus SALT2 is a like-for-like test at fixed parameter count,
before any new freedom enters. Everything above L0 asks whether added freedom pays.

Standardisation is `mu_corr = mu - alpha*theta_1 - beta*c - gamma*theta_2 + M`; the first
three terms are SALT2's own, one-to-one, and **`gamma` is the headline number** — the
coefficient of the one DOF SALT2 lacks. Scored on **held-out Hubble residual scatter**
under leave-one-out cross-validation, following SALT2's own validation method. The ladder
is also the **training schedule**: fit L0 to convergence, then introduce latents one at a
time, as SALT2 stages its components.

### The mechanism difference, and the `w` versus `x1` check

SALT2 has *no* way to reparameterise time. A stretch `p -> p/s` leaves the curve pointwise
unchanged and alters only its traversal; SALT2 instead uses the additive `x1*M1(p,lambda)`.
These agree only to first order — `m(p/s) = m(p) - (s-1)*p*m'(p) + O((s-1)^2)`, so
`M1 ∝ -p dM0/dp`. This is not a reconstruction after the fact: **SALT2 initialises `M1` as
exactly that finite difference**, the sequence at `s=1.1` minus that at `s=1`. The published
`x1 -> s` conversion needs a cubic, and that cubic *is* the residual nonlinearity the
linearisation leaves behind.

Two consequences:

- The dominant mode of SN Ia diversity is the one this model represents **exactly** and
  SALT2 **linearises**. That is the strongest structural argument for the approach.
- With stretch removed into `w`, what remains for the first shape latent is only variation
  no reparameterisation can produce, so **it should be small**. Concrete check: fitted `w`
  should track `s = 0.98 + 0.091*x1 + 0.003*x1^2 - 0.00075*x1^3`, and the scatter about that
  relation measures what the linearisation misses. Cheap, and the sharpest early test of
  whether separating shape from timing buys anything — worth running before the full
  pipeline exists.

### Per-SN parameters and timing conventions

Baseline eight: `mu` (normalization) · `c` (colour) · `theta_1, theta_2` (shape, condition
the network) · `w, t_max` (timing) · `a_1, a_2` (timing warp). SALT2 uses four; the excess
is the object of study.

Conventions that should survive whatever the parameter count becomes: `w` is a timing scale
in **rest-frame days**; `t_max` is **observer-frame MJD**, always free per SN and never
taken from the archive's SALT2 `t0`; `sigma` is reserved for flux measurement errors and
must not be used for the timing scale; time dilation is deterministic, `p = (t - t_max)/(1+z)`
divided out up front like `mwebv`, so `w` is a genuine rest-frame stretch and not a
repackaging of `(1+z)`. The traversal is the cubic `s(t) = u[1 + a_1 u + a_2 u^2]` with
`u = (t-t_max)/(w(1+z))`, on a bounded arclength domain held in config.

### Curvature sum rule — qualitative only, never imposed

Flux vanishes before explosion and after, so `m -> +inf` in every band at both ends: the
curve is a **hairpin**, both ends escaping to infinity, the tangent reversing between them
and turning through **approximately `pi`**. Corollaries independent of the escape direction:
`kappa -> 0` at both extremes, and the asymptotic ray separation is the terminal colour —
predicted only *relative* to the early colour, since a free `c` shifts both together.

Both asymptotes lie outside the bounded model domain and neither is ever observed, so
**nothing here is imposed at all**: the hairpin, `kappa -> 0` at the extremes,
`∫ kappa ds ≈ pi` and the terminal colour are qualitative validation statements about a
trained `kappa`'s *extrapolation*. **Never impose any of them as a constraint.** By unit
speed, `f = 0` requires `|s| = infinity` — the dark phase is a single ideal point in the
closure, never a point of the domain. The earlier exact forms — `∫ kappa ds = pi (mod 2pi)`,
exact `(1,1)` asymptote directions, and a branch table splitting the integral at the
`g`-maximum — are deleted, the branch table because peak ordering is not a quantity of
interest.

### No per-SN rotation

No **per-SN** rotation is a physical degree of freedom: distance and zeropoint errors are
**translations**, stretch is a **reparameterisation**, and dust is a **phase-dependent
displacement** — not a translation. None acts as a per-SN rotation, so none is fitted per
SN. This is also why a per-SN reddening *direction* is not fitted: it would be a per-SN
rotation of the colour axis with no physical referent.

The anchor is `g`-maximum, **not** `B`-maximum, since locating `B` max needs an SED. So
`t_max` is not comparable to SALT2's `t0`; measure the offset distribution and assume
neither its centre nor its sign.

### Gauge, assumption, prediction, parameter

The four-way classification, as it stood when it left `CLAUDE.md`. Note that the note's
current `psi_0` paragraph has already moved band ordering out of *gauge* and into stated
physical knowledge, which this list does not yet reflect.

- *Gauge* (costless, not testable): the arclength origin at the principal `g`-maximum; the
  frame orientation for `n = 2`, by the same condition; the metric `diag(1,1)`; latents
  normalised to zero mean and unit variance with the sign fixed so the first correlates
  positively with light-curve width; position fixed by `gamma(0) = 0`; and `e_c` below L2c.
- *Assumption* (substantive, falsifiable): `Corr(c, theta) = 0`.
- *Prediction* (checked, not chosen): colour **evolution**; the hairpin, `∫ kappa ds ≈ pi`
  and the terminal colour, as qualitative statements about extrapolation outside the domain.
  **Not** "which band peaks first" while that is an input rather than an output.
- *Parameter*: `c`, the amplitude along `e_c`, free per SN; `e_c` itself **at L2c only**, one
  global unit vector for the sample, not per-SN; and for `n >= 3` the `n(n-1)/2 - 1` global
  orientation numbers, one fewer than `SO(n)` by the orientation/origin redundancy.

With `c` free, neither early nor terminal colour is predicted absolutely; only their
*difference* is, since `c` displaces the whole curve and cancels. The model predicts colour
evolution and leaves the colour zero point free — structurally the split SALT2 makes between
a fixed colour law and a fitted `c`.

### Degeneracies to report

- **Two bands admit no torsion.** In the primary `g,r` analysis a second shape latent is a
  second direction of variation in `kappa`, *not* torsion. Do not call it torsion in
  writing. Torsion is meaningful only on the 177-SN `g,r,i` subsample.
- **Timing warp versus shape.** Over a finite, noisy window a warp in `s(t)` and a change in
  `kappa` can mimic each other. Report the fitted correlation between `(a_1, a_2)` and the
  shape latents. A large correlation means the shape/timing separation is not being
  realised — a negative result to report, not to absorb.
- **Dust versus shape.** Same structure, expected milder, for the reasons in the order-counting
  entry above. Report stratified by fitted `c`; consequences are asymmetric — harmless for the
  distance indicator, a limitation for interpreting the shape latents.
- **Dust versus intrinsic colour** is the first-order degeneracy and is unresolved; `c`
  conflates the two exactly as SALT2's `c` does. The order counting does not help here.
- **Peculiar velocity** is not in the likelihood: `mu` is free per SN and absorbs it exactly,
  entering only at the Hubble-diagram stage.
