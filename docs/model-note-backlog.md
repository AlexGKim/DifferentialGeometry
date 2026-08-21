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

**Settled as a working value; the empirical question is what the ladder answers.** Moved
out of the note on 2026-08-21 — it is a justification of a modelling choice, i.e. argument,
not model definition.

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

## `e_c` is exactly unidentifiable below L2c

**Settled, 2026-08-22.** The reddening direction was made a sample-wide fit parameter
earlier the same day; the same afternoon showed that below L2c there is nothing to fit.

With `mu` free per SN along `(1,1)` and `c` free per SN along `e_c`, the two span `R^2` for
any `e_c` not parallel to `(1,1)`. So changing `e_c` is undone **exactly** by an invertible
relabelling of `(mu, c)`: a flat likelihood direction, not a slow one. At L0–L2, `e_c` is
therefore **gauge** — fix it at `(1,-1)/sqrt2` and do not fit it; a round-trip that appears
to recover it is a bug.

It gains content only at **L2c**, through the phase-varying `du(s)`, which is not a
translation and so cannot be absorbed by any relabelling of two amplitudes. There it is a
genuine global parameter, the two-band analogue of SALT2's fitted `beta`, with
`(1,-1)/sqrt2` the nested fixed-`R_V` idealisation.

The general lesson is worth stating once, because it recurs: **a direction is only
identifiable if something in the model varies along it in a way the free amplitudes cannot
mimic.** Two free per-SN amplitudes spanning the plane make every phase-*independent*
direction question vacuous.
