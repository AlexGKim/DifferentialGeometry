# Runs ledger

`runs/` is gitignored, so a run is 11 MB of CSVs on one disk with nothing in the
repository saying what produced it. This file is that record: for each run worth
deriving from, the inputs by hash, the command, and the diagnostics claimed, so a
later re-run can be *checked* rather than trusted. Entries are append-only; a
superseded entry is marked, not deleted.

What belongs here is a run later work starts from or measures against. Ordinary
exploratory runs do not — `docs/model-note-backlog.md` is where their lessons go.

---

## stage3-k88-free — the first converged fit

`runs/stage3-k88-free`, finished 2026-11-03.

The first fit here that CmdStan's `diagnose()` clears without a complaint. It is
a reference point and a starting condition for later work. It is **not** a
result, and nothing in it establishes that `K = (8, 8)` is the right order —
see Known gaps.

### Inputs

| | |
|---|---|
| commit | `024fd46` — the code state the run was launched from; `plot_fit.py` changed after it and does not affect the fit |
| program | `stan/hier_common.stan` — Eq. (curvlin) and Eq. (travlog) with no departure modes |
| `free_mean` | 1, so ⟨Δ⟩ is sampled |
| orders | `Kq = 8`, `Kpsi = 8` |
| selection | `dgsn.data.sample.stage3_sample()` — 30 objects, 1568 epochs |
| `data.json` | sha256 `c3c4b953b3216e89…` |
| `inits-1..4.json` | sha256 `f45ea0f42fb1d98f…` — **all four identical**, one SALT2-placed start |
| CmdStan | 2.39.0 |

The object list is identified by reproducing `data.json`'s per-object epoch
counts, not by which script is presumed to have made the run:
`scripts/plot_fit.py:resolve_names` does this and raises otherwise. It matters —
`stage3_sample()` and a thirty-object `stratified()` are both thirty long and
share one object, and the plots were wrong until this was checked.

### Sampler

```
seed 1, 4 chains, 500 warmup + 500 sampling, max_treedepth 11,
diag_e, adapt_delta at CmdStan's default
```

```bash
python scripts/stage3_hier.py --no-trav-mode --k-psi 8 --k-q 8 \
  --iter 500 --max-treedepth 11 --out runs/stage3-k88-free
```

### Diagnostics

| | |
|---|---|
| `diagnose()` | no problems detected |
| R-hat, max | 1.01425 at `lp_sn[18]`, a generated quantity |
| R-hat, worst sampled parameter | 1.00991 at `x_s[5]` |
| treedepth | mean 8.02, max 9, none at the cap |
| divergences | 0 of 2000 |
| realised `sqrt(cond)` | 119 over 138 sampled coordinates |
| wall time | 6615–7293 s warmup per chain |

`sqrt(cond)` is taken over the sampled coordinates with one slot dropped from
each of the three `sum_to_zero_vector` blocks (`x_s`, `c_dev`, `d_dev`).
Including the redundant slots, or the transformed parameters, leaves exact zero
eigenvalues and reports 1.2e17.

### Values

| | |
|---|---|
| `sum lp_sn` | −5769.297 mean, −5745.255 best draw |
| chi²/epoch | 7.359 |
| ⟨Δ⟩ | −0.3951 ± 0.0959 d, 4.1 sd from zero |
| population spread of Δ | 0.3029 d |
| `c_bar` | +0.0728 ± 0.0024 |
| `b_0` | −7.99 ± 0.05 |
| `tau_psi` | 2.046 ± 0.712 — single mode, no sign of the two-mode split at (24,16) |
| `alpha_q` | 1.219 mean; `alpha_psi` 0.923 mean, neither against its ½ floor |
| `kappa0` | positive throughout, so Eq. (peakcond) holds |

### Known gaps

1. **No dispersed-start check.** All four `inits` are byte-identical, so R-hat
   here measures momentum-seed dispersion and not start dispersion. The SALT2
   start is exempt from "templates never to fit" only because dispersed-start
   agreement checks it (CLAUDE.md), and that check has not been run. This
   model was bimodal at (24, 16), where agreement proved to be within-mode.
2. **One resolution.** The `K_psi` ladder at 12 and 16 was not run, so "two
   resolutions, always" is unmet on the parameter that matters most here.
3. **Fit quality.** 1680 nats behind `runs/ref-k2416-eps3` at (24, 16) on the
   same 30 objects and 1568 epochs, and the misfit is visible rather than
   statistical: ZTF18abkhcwl's r-band shoulder at phase 18–27 is systematically
   under-predicted. The reference never converged and is down to two readable
   chains, so its likelihood is from a truncated sampler, but the gap is far too
   large to be an artifact of that.
4. **No out-of-sample score.** `src/dgsn/eval/` is empty, so nothing here says
   this predicts better than anything else. "SALT2 strictly downstream" and
   "evaluation is out-of-sample" both remain to be re-established.

### Companion

`runs/stage3-k88-dzero` is the same configuration with `free_mean = 0`, run to
price ⟨Δ⟩ = 0 and nothing else. It also converged; the constraint costs 8.9
nats of `sum lp_sn` at the best draw and buys no mixing. That figure replaces
the 1.96 nats quoted around the source tree, which came from a comparison in
which both arms sampled a free mean.

---

## stage3-k88-n60 — the same orders at twice the objects

`runs/stage3-k88-n60`, finished 2026-11-04.

The second fit here that `diagnose()` clears without a complaint, and the first
at sixty objects. It is a reference point and a starting condition. It is
**not** a result: its fit quality is materially worse than
`stage3-k88-free`'s and the misfit is systematic rather than statistical —
see Known gaps.

### Inputs

| | |
|---|---|
| commit | `05dbb73` plus `scripts/stage3_hier.py`'s `--save-warmup` flag, committed immediately after as part of this work |
| program | `stan/hier_common.stan` — unchanged from `stage3-k88-free` |
| `free_mean` | 1, so ⟨Δ⟩ is sampled |
| orders | `Kq = 8`, `Kpsi = 8` — identical to `stage3-k88-free` |
| selection | `stage3_sample(n_per_fifth=12)` — 60 objects, 3682 epochs |
| `data.json` | sha256 `00e0cd556ecfd05f…` |
| `inits-1..4.json` | sha256 `98b7a2d8c4859351…` — **all four identical** |
| initialised from | `runs/stage3-k88-free` — 30 objects carried by name, 30 placed from \salt; all 14 coefficients carried directly, no padding, since the orders match |
| CmdStan | 2.39.0 |

The object list is identified by reproducing `data.json`'s per-object epoch
counts, as `scripts/plot_fit.py:resolve_names` does. Note that `resolve_names`
calls `stage3_sample()` at its default `n_per_fifth = 6` and so cannot identify
this run; the check here passed `n_per_fifth = 12` explicitly.

### Sampler

```
seed 1, 4 chains, 500 warmup + 500 sampling, max_treedepth 12,
save_warmup on, diag_e, adapt_delta at CmdStan's default
```

```bash
python scripts/stage3_hier.py --no-trav-mode --k-psi 8 --k-q 8 \
  --n-per-fifth 12 --init-from runs/stage3-k88-free \
  --iter 500 --max-treedepth 12 --save-warmup --out runs/stage3-k88-n60
```

### Diagnostics

| | |
|---|---|
| `diagnose()` | no problems detected |
| R-hat | 0 of 653 coordinates above 1.01; median 1.0000, worst 1.0084 at `lp_sn[8]` |
| treedepth | mean 9.48, max 10, none at the cap of 12 |
| divergences | 0 of 2000 |
| ESS | minimum 787 at `mu_delta`, median 1904 |
| realised `sqrt(cond)` | 1396 over 258 sampled coordinates |
| wall time | 23502–26347 s warmup, 14281–16371 s sampling, 37783–42553 s total |

`sqrt(cond)` is taken over exactly the coordinates the `parameters` block
declares, with one slot dropped from each of the three `sum_to_zero_vector`
blocks (`x_s`, `c_dev`, `d_dev`): $14 + 4 + 60 + 1 + 59 + 1 + 59 + 59 = 258$,
the same accounting that gives 138 at thirty objects. Including the
transformed parameters or the redundant slots leaves exact zero eigenvalues
and reports a NaN.

**Warmup divergences are not a defect here, and nearly caused a wrong call.**
47 to 62 divergences per chain arrived during warmup and persisted with the
step size apparently settled, which was read as a funnel. All of them were
gone post-warmup — zero in 2000 draws — once the term buffer finished tuning
the step size against the final metric. Divergence counts, like tree depth,
are not to be read from warmup; see the backlog entry of the same date.

### Values

| | |
|---|---|
| `sum lp_sn` | −32615.589 mean, −32580.192 best draw |
| chi²/epoch | 17.716 over 3682 epochs |
| ⟨Δ⟩ | **+0.4118 ± 0.0428 d** — positive, where `stage3-k88-free` gave −0.3951 ± 0.0959 |
| population spread of Δ | 0.2873 d, per-object means in $[-0.321, +1.029]$ |
| posterior sd of Δ | 0.0560 d, median over objects |
| `c_bar` | +0.1106 ± 0.0017, against +0.0728 ± 0.0024 at thirty |
| `b_0` | −8.1440 ± 0.0222 |
| `tau_psi` | 2.5277 ± 0.8831 — single mode, agreeing across all four chains |
| `tau_q` | 1.2299 ± 0.4511 |
| `alpha_q` | 1.1994 ± 0.4572; `alpha_psi` 0.8585 ± 0.2690, neither against its ½ floor |
| `kappa0` | +6.4374 ± 1.0133, positive throughout, so Eq. (peakcond) holds |

### What doubling the objects cost, at fixed order

| | 30 objects (`k88-free`) | 60 objects (this run) |
|---|---|---|
| sampled coordinates | 138 | 258 |
| `sqrt(cond)` | 119 | 1396 |
| treedepth | mean 8.00, max 8 | mean 9.48, max 10 |
| step size | 1.64e−2 | 5.5–5.8e−3 |
| divergences | 0 of 2000 | 0 of 2000 |
| chi²/epoch | 7.359 | 17.716 |

`sqrt(cond)` rose by a factor of 11.7 while the realised trajectory grew by
about 1.5 doublings, i.e. 2.8× the leapfrogs — so $\sqrt{\mathrm{cond}}$
is a loose upper guide here, over-predicting the leapfrog count by 2× (1396
against the 714 that depth 9.48 implies). Worth remembering before any
$\sqrt{\mathrm{cond}}$ estimate is used as a wall-clock budget.

### Known gaps

1. **Fit quality, and it is systematic.** chi²/epoch 17.716 against 7.359 at
   thirty. The 30 carried objects are the *identical* objects and epochs and
   now fit at 9.205, a 25% degradation, so the global functions are being
   pulled; the 30 new objects fit at 24.029. The misfit is a tail, not a
   uniform loss — the median object is at 7.98, close to `k88-free`'s 7.359,
   while one object holds 29% of the total chi² and eleven hold 69%:

   | name | epochs | chi²/epoch | z | x1 | c |
   |---|---|---|---|---|---|
   | ZTF18aasdted | 222 | 85.79 | 0.0181 | +0.71 | +0.152 |
   | ZTF18aaxcntm | 247 | 39.32 | 0.0268 | −1.44 | +0.041 |
   | ZTF19aarlefl | 22 | 33.76 | 0.0435 | +0.18 | +0.003 |
   | ZTF20abeywdn | 80 | 32.48 | 0.0254 | +0.78 | −0.130 |

   `ZTF18aasdted` is the lowest-redshift object in the sample with 222 epochs,
   so it carries the smallest error bars anywhere in the fit and dominates
   chi² because the model's systematic error there exceeds its statistical
   error. `(8, 8)` is converging on a sample it cannot describe.
2. **⟨Δ⟩ changed sign.** +0.4118 ± 0.0428 here against −0.3951 ± 0.0959 in
   `stage3-k88-free`, on the same program at the same orders with only the
   object count changed, and all four chains agreeing to ±0.006 so it is not a
   mode artifact. Zero remains excluded in both, so nothing here disturbs
   ⟨Δ⟩ ≠ 0; what it disturbs is the *value*, which does not transfer across
   sample size and is plausibly absorbing systematic error from the badly-fit
   objects above. Any quoted figure for ⟨Δ⟩ must carry its sample and orders.
3. **No dispersed-start check.** All four `inits` are byte-identical, so R-hat
   here measures momentum-seed dispersion and not start dispersion — the same
   gap as `stage3-k88-free`. Momentum seeds alone split
   `runs/stage3-k1216-n60` into two basins, so this is not a formality.
4. **One resolution in K.** `(8, 8)` only. The ladder measured +1653 nats from
   `(8, 8)` to `(16, 24)` at thirty objects, so the fit quality above is
   expected and the order is the suspect; but `(12, 16)` saturated the depth
   cap at these sixty objects and `(12, 12)` is unresolved rather than
   refuted. "Two resolutions, always" is unmet.
5. **No out-of-sample score.** `src/dgsn/eval/` is still empty. "SALT2
   strictly downstream" and "evaluation is out-of-sample" both remain to be
   re-established.

---

## stage3-k88-n60-floor — a fitted per-band model error, and what it bought

`runs/stage3-k88-n60-floor`, finished 2026-11-05.

The same orders, the same sixty objects and the same `data.json` as
`stage3-k88-n60` — sha256 `00e0cd556ecfd05f…` in both — with one change to the
likelihood: a fitted per-band Gaussian model error added in quadrature,
Eq. (like). It is a controlled comparison of that one term.

It fixes what it was built to fix and it is substantially cheaper to sample.
It is **not** converged to the bar `stage3-k88-n60` met, and the geometry it
finds is a long way from that run's — see Known gaps.

### Inputs

| | |
|---|---|
| commit | the floor added to `stan/dgsn_functions.stan`, `stan/hier_common.stan`, `src/dgsn/stan/hier_data.py` on top of `bb2bc93` |
| program | `stan/hier_common.stan` with `sigma_b` sampled and `partial_sum_floor` in `reduce_sum` |
| `free_mean` | 1, so ⟨Δ⟩ is sampled |
| orders | `Kq = 8`, `Kpsi = 8` |
| selection | `stage3_sample(n_per_fifth=12)` — 60 objects, 3682 epochs |
| `data.json` | sha256 `00e0cd556ecfd05f…` — **identical to `stage3-k88-n60`** |
| `inits-1..4.json` | sha256 `54f24c7744d6c0a8…` — all four identical, `sigma_b` entering at 0.05 mag in each |
| initialised from | `runs/stage3-k88-n60` |
| CmdStan | 2.39.0 |

### Sampler

```
seed 1, 4 chains, 500 warmup + 500 sampling, max_treedepth 12,
save_warmup on, diag_e, adapt_delta at CmdStan's default
```

```bash
python scripts/stage3_hier.py --no-trav-mode --k-psi 8 --k-q 8 \
  --n-per-fifth 12 --init-from runs/stage3-k88-n60 \
  --iter 500 --max-treedepth 12 --save-warmup --out runs/stage3-k88-n60-floor
```

### Diagnostics

| | |
|---|---|
| `diagnose()` | no problems detected — **but see the R-hat row; `diagnose()` complains at 1.05, not at this project's 1.01** |
| R-hat | **62 of 778 coordinates above 1.01**; median 1.0024, worst 1.0141 at `w_coef[6]` |
| treedepth | 8.00 in every one of 2000 draws, max 8, cap 12 |
| divergences | 0 of 2000 |
| ESS | minimum 364 at `b_0`, median 1422 |
| realised `sqrt(cond)` | 248 over 260 sampled coordinates |
| wall time | 17635 s warmup, 6936 s sampling, 23739–24571 s total |

The 260 coordinates are the 258 of `stage3-k88-n60` plus `sigma_b[1..2]`, taken
from the `parameters` block with one slot dropped from each
`sum_to_zero_vector`. `sqrt(cond)` was read at three draw counts — 310 at 596,
301 at 646, 248 at 2000 — so it is still falling with draws and 248 is an upper
bound, as 1396 is for the unfloored run.

**Tree depth of exactly 8 in every draw is not a cap.** Every draw also spends
exactly $2^8-1 = 255$ leapfrogs, which is the signature of a depth cap, and was
read as one here before `stage3-k88-n60` settled it: that run occupied depths 9
and 10 and *nothing else* in 2000 draws. Once the step size is frozen after
warmup this model's depth is concentrated in one or two bins, and
`n_leapfrog__` is exactly $2^d - 1$ in 100% of this run's draws and 82–84% of
the unfloored run's, so it carries nothing beyond the depth.

### Values

| | |
|---|---|
| `sigma_b` (g) | **0.0469 ± 0.0014 mag** |
| `sigma_b` (r) | **0.0679 ± 0.0014 mag** — 45% larger than g, so the per-band form earns its place |
| chi²/epoch on $V$ | 1.1089 over 3682 epochs; per object, median 1.035, max 2.733 |
| chi²/epoch on archive errors | 40.538, against 17.716 unfloored |
| `sum lp_sn` | −74630.811 mean, −63616.994 best draw (archive errors, no floor, so comparable with the unfloored −32615.589) |
| ⟨Δ⟩ | +0.4345 ± 0.1979 d — the mean agrees with the unfloored +0.4118, the error bar is 4.6× wider |
| population spread of Δ | 0.3460 d, per-object means in $[-0.551, +0.947]$ |
| posterior sd of Δ | 0.2251 d median over objects, against 0.0560 unfloored |
| `c_bar` | +0.0922 ± 0.0065, against +0.1106 ± 0.0017 |
| `b_0` | −7.7681 ± 0.1416, against −8.1440 ± 0.0222 |
| `tau_psi` | 2.2028 ± 0.7846 — single mode, all four chains |
| `tau_q` | 0.9036 ± 0.3429 |
| `alpha_q` | 1.6645 ± 0.6675; `alpha_psi` 0.9099 ± 0.3058 |
| `kappa0` | **+52.7097 ± 5.3852**, positive in 100% of draws, against +6.4374 ± 1.0133 unfloored |

### What the floor fixed

The term was added because the misfit in `stage3-k88-n60` tracked how well each
supernova was *observed* rather than anything about the supernova — the
signature of model error being absorbed by whichever objects had the smallest
error bars. Measured on all 2000 draws:

| Spearman correlation of chi²/epoch with | unfloored | floored |
|---|---|---|
| peak S/N | **+0.862** | **−0.086** |
| median `sig_f`/`f` | −0.797 | +0.134 |
| number of epochs | +0.316 | +0.004 |

The dependence is not reduced, it is gone. One object held 29% of the total
chi² before; it now takes five to hold 33.1%.

### What the floor cost, and what it bought

| | `stage3-k88-n60` | this run | factor |
|---|---|---|---|
| `sqrt(cond)` at 2000 draws | 1396 | 248 | 5.6× easier |
| leapfrogs per draw | ~767 | 255 | 3.0× fewer |
| step size | 5.5e−3 | 1.41–1.60e−2 | 2.7× larger |
| treedepth | 9.50, bins 9–10 | 8.00, bin 8 only | |
| wall time | 37783–42553 s | 23739–24571 s | 1.6–1.8× faster |
| divergences | 0 of 2000 | 0 of 2000 | |
| coordinates above R-hat 1.01 | 0 of 653 | 62 of 778 | **worse** |
| minimum ESS | 787 | 364 | **worse** |

Better conditioning and worse mixing at once. The misses concentrate in the
per-object phase quantities — 30 of the 62 are `s_pmin`, 19 are `delta_t` — and
`mu_delta` itself sits at R-hat 1.0113 with ESS 535. That is coherent with the
mechanism: the floor caps how tightly any single epoch can pin a phase, which
is exactly how it flattens the stiff directions, and the same cap leaves the
per-object phases less determined. The 4.6× inflation of ⟨Δ⟩'s error bar and
the 4× inflation of Δ's per-object posterior sd are the same effect priced in
the parameter it matters for.

Also note that `sqrt(cond)` now predicts the leapfrog count to 18% (248 against
255) where it over-predicted by 1.8× unfloored. The earlier warning that
`sqrt(cond)` is a loose guide was measured on a posterior whose stiffness came
from a handful of over-precise objects; with those capped, the diagonal-metric
argument applies.

### Known gaps

1. **Not converged by this project's bar.** `diagnose()` passing means only
   that nothing exceeds 1.05. Sixty-two coordinates exceed 1.01 and the
   minimum ESS is 364. This is a reference point, not a result, and a longer
   run or a reparameterisation of the per-object phase block is wanted before
   anything is claimed from it.

2. **`kappa0` moved by a factor of eight**, from +6.44 ± 1.01 to +52.71 ± 5.38,
   and `b_0` and `c_bar` moved by 17 and 11 of the unfloored run's standard
   deviations. The floor is not a cosmetic reweighting; it relocates the
   geometry. Which location is right is not settled here, and `lp__` cannot
   settle it because the two runs have different likelihoods.

3. **A per-band constant absorbs phase-structured error uniformly.** A
   well-fitted `sigma_b` says the error budget is honest, not that the curve is
   right. Per-epoch residuals as a generated quantity are what would show
   *where* the model fails; they are not implemented.

4. **No out-of-sample score**, so nothing here says the floored model predicts
   better. With a fitted scatter this matters more than before: held-out
   chi²/σ_f² no longer ranks orders, and selection must use the held-out log
   score with its $\sum \ln V$ term. `src/dgsn/eval/` is still empty.

## stage3-k1216-n60-floor — restoring the orders, and why it failed

`runs/stage3-k1216-n60-floor`, launched 2026-11-06, **stopped by hand 2026-11-07**
at 529/928/545/564 draws of 1000. A negative result, kept because it localises
the obstruction.

The hypothesis under test was that the Eq. (like) model error removes the
convergence failure of `stage3-k1216-n60`, which was abandoned at 22 h a chain
with 99.6% of transitions at a `max_treedepth` of 11 and a 2–2 split in
`mu_delta` against `c_bar`. Three things differed from that run, not one: the
likelihood, `max_treedepth` 12 rather than 11, and a start taken from a
converged sixty-object fit rather than a thirty-object one.

**The hypothesis is false as stated.** Both failure modes survived. One thing
did improve, and it is attributable to the model error.

### Inputs

| | |
|---|---|
| commit | `199978a` plus the uncommitted `k_a`/`k_c` additions; the program itself is `199978a`'s |
| program | `stan/hier_common.stan` — no departure modes, `sigma_b` sampled |
| `free_mean` | 1, so ⟨Δ⟩ is sampled |
| orders | `Kq = 12`, `Kpsi = 16` |
| selection | `stage3_sample(n_per_fifth=12)` — 60 objects, 3682 epochs |
| `data.json` | sha256 `36751d7858ef9a40…` |
| `inits-1..4.json` | sha256 `7ae4f684f983112e…` — all four identical, so no dispersed-start check |
| initialised from | `runs/stage3-k88-n60-floor` |
| sampled coordinates | 272, against the reference's 260 |

### Sampler end state, post-warmup

| | ch 1 | ch 2 | ch 3 | ch 4 |
|---|---|---|---|---|
| draws past warmup | 29 | 428 | 45 | 64 |
| fraction at depth ≥ 11 | **1.00** | 0.39 | **1.00** | **1.00** |
| fraction at the cap of 12 | 0.70 | 0.00 | 0.10 | 0.00 |
| mean `n_leapfrog__` | 3564 | 1598 | 2535 | 2081 |
| frozen step size | 1.48e-3 | 2.95e-3 | 1.56e-3 | 1.66e-3 |
| divergences in warmup | 75 | 67 | 79 | 76 |
| divergences post-warmup | 0 | 0 | 0 | 0 |
| `mu_delta` | −1.76 | **+0.23** | −1.81 | −1.77 |
| `c_bar` | +0.032 | **+0.130** | +0.050 | +0.028 |

Against the reference `stage3-k88-n60-floor`: depth 8.00 in all 2000 draws with
nothing at the cap, step size 1.4–1.6e-2, `mu_delta` +0.4345. **The step size
froze an order of magnitude smaller and three chains sit at depth ≥ 11 in every
post-warmup draw.** The extra treedepth granted over `stage3-k1216-n60` was
spent, not saved.

The warmup divergence counts in the table above are **not** a failure
signature: `stage3-k88-n60-floor` converged with 65/68/71/63 over its own first
500 rows and `stage3-k1212-n60-floor` failed with 61/69/64/68, so all three are
the same to within their scatter. See the correction in
`stage3-k1212-n60-floor`. The discriminating columns are step size, treedepth
and leapfrog count.

`treedepth__` *means* are a trap here and were briefly misread as healthy: a
divergent transition terminates at once and is recorded at depth 0–6, so a
window with half its draws at the cap of 12 and five divergences reports a mean
of 9.92 while using 2667 leapfrogs. Read the histogram, or the leapfrog count,
never the mean depth alone.

### The split returned, as 1–3

Chains 1, 3, 4 fell to `mu_delta` ≈ −1.8 by draw 150 and stayed; chain 2
random-walks over ±1.5 across the whole of warmup and never settles. Chain 2 is
also the only chain not pinned at depth, and ran 1.7× faster.

What improved: the `lp__` gap between the groups is 6–14 nats, against 58 nats
in `stage3-k1216-n60`, and the majority now holds the higher-density region
(`lp__` −32134/−32126/−32129 against chain 2's −32140). Narrowing that gap is a
real effect of the model error. It is not sufficient.

R̂ over the 275 CSV columns of the sampled blocks, dropping one slot per
`sum_to_zero_vector`, on the 27 common post-warmup draws: **9.59 over four
chains, and 3.01 over chains 1, 3, 4 alone**, with 120 of 275 above 1.01 inside
that agreeing subset and every worst offender a `w_coef` slot. The three chains
that agree on the anchor disagree on the global shape. At n = 27 the R̂ noise
floor is near 1.1, so 3.0 is not a small-sample artefact.

### Where the obstruction is — the Fisher spectrum over the coefficients

Computed at draws of this run, with the Eq. (like) weight `f²/V` rather than
the archive weight (see Known gaps), prior precision included:

| block | coefficients | sqrt(cond) | λ_min |
|---|---|---|---|
| `b_til` (q₀, Kq = 12) | 11 | 55.2 | 1.07e+2 |
| `d_til` (Ψ₀, Kpsi = 16) | 15 | 720.1 | 7.31e−1 |
| both | 26 | 966.3 | 4.96e−1 |

**The cost is entirely the turning order.** The four lowest eigendirections put
99.0–99.9% of their weight on `d_til`; the lowest, at λ = 0.496, loads on
alternating odd-index turning coefficients (`d7` −0.52, `d5` +0.41, `d3` −0.40,
`d9` +0.34, `d11` −0.29) — an oscillatory combination the data cannot resolve.
`Kq = 12` is nearly free: on its own it conditions at 55, against 43.9 for the
whole coefficient block at `(8,8)`.

Truncating Ψ₀ while holding `Kq = 12`, at three points of this run:

| | Kpsi=8 | Kpsi=10 | Kpsi=12 | Kpsi=14 | Kpsi=16 |
|---|---|---|---|---|---|
| best `lp__` | 106 | 172 | 271 | **805** | 966 |
| median `lp__` | 76 | 136 | 195 | **331** | 605 |
| best with `mu_delta` > 0 | 94 | 176 | 240 | **558** | 683 |
| λ_min, best `lp__` | 37.7 | 14.5 | 6.15 | **0.71** | 0.50 |

**The cliff is between Kpsi = 12 and 14**, at all three points, and it survives
forcing `tau_psi` to the reference run's 1.267 in place of this run's inflated
value (260 → 590 instead of 271 → 805). So it is a property of the data and the
basis, not of the hyperparameter. `Kpsi = 12` is the highest turning order that
stays conditioned — and it is the order at which `stage3-fisher` and
`stage3-dense` did detect both departure modes.

### The hyperparameters press the boundary

Over all draws:

| | mean | sd | min | max | fraction below 0.05 |
|---|---|---|---|---|---|
| `tau_q` | 0.6311 | 0.1845 | 0.2442 | 1.7393 | 0.000 |
| `tau_psi` | 5.6500 | 2.4313 | 0.6399 | 14.6377 | 0.000 |
| `alpha_q_m` | 0.5222 | 0.3902 | 0.0007 | 2.6480 | 0.067 |
| `alpha_psi_m` | 0.2663 | 0.2234 | 0.0000 | 1.7949 | **0.141** |

`alpha_psi_m` is declared `real<lower=0>`, so a posterior that wants zero runs
to −∞ on the unconstrained log scale. It reaches 0.0000 and spends 14% of draws
below 0.05, against 0.212 at `(8,8)`; `tau_psi` meanwhile inflated 4.5× to 5.65
from 1.267. Together these flatten the coefficient prior exactly where the
near-null direction lives, and the flat direction in turn pulls the exponent
down — the two findings above are one feedback, not two. Whether `alpha_psi_m`
should be bounded away from zero is a question for `model-definition.tex`, not
for the driver.

### A fixed Fisher rotation is not the remedy

`precondition.rotation` gives `R'FR = I` at its build point by construction.
Built at the best draw and evaluated elsewhere in the same posterior:

| evaluated at | sqrt(cond) unrotated | rotated by `R` |
|---|---|---|
| the build draw | 966.3 | 1.000 |
| median `lp__` draw | 605.0 | 240.7 |
| lowest `lp__` draw | 518.7 | **907.4** |

At the lowest-density draw the rotation leaves the geometry **worse than it
found it**. The curvature varies too much across this posterior for one linear
rotation, so `--precondition-from` is measured not to be the first remedy to
reach for at these orders. This removes the leading item from the fallback list.

### Known gaps

- **`precondition.fisher` weights by the archive error alone.** Its docstring
  states this as current — "the variance is the archive's error alone, as the
  fit's now is" — which `199978a` made false when it added the fitted per-band
  term. `sigma_b` is absent from `HYPER`, so `pieces` never returns it and
  `fisher` cannot see it. The Fisher spectra above therefore had to be recomputed
  by hand with `V_i = σ_f,i² + (0.4 ln10 · σ_b(i) · f_i)²`; under the stale
  archive weight the same point reads 1227.1 rather than 966.3. Both
  `--precondition-from` and `--precondition-optimize` build their rotation from
  the superseded likelihood.
- **`precondition._draws` returned 2558 of the 2566 rows** of this interrupted
  run, i.e. warmup included, despite documenting itself as post-warmup only. The
  points used above are legitimate points of the space and the highest-`lp__`
  one is not a warmup artefact, but `best_draw` on an unfinished run is not
  selecting what its docstring says.
- Starts were byte-identical across chains, so nothing here is a dispersed-start
  check.
- Nothing is out of sample; `src/dgsn/eval/` is still empty.
- The run never reached 1000 draws, so none of its parameter values are quoted
  here as measurements of anything but its own failure.

## gauge-nullity — the gauge construction costs sqrt(cond) proportional to tau_psi

A measurement across archived runs, **not a run**. Made 2026-11-07 while
diagnosing `stage3-k1216-n60-floor`, to answer whether fixing Eq. (gauge) by
projection is itself a cause of the convergence failures.

### Only the mode-bearing programs fix the gauge at all

`stan/hier_common.stan` sets `b = b_til; d = d_til;` — no normalisation, no
projection. With one function per pair there is no gauge freedom to fix, q₀
having no partner to be orthogonal to. So **no failure of `hier_common` can be
attributed to the gauge**, which disposes of the question for
`stage3-k1216-n60-floor` and for the `(12,12)` run that followed it.

`stan/hier.stan` and `stan/hier_modes.stan` carry the construction in full:
`a = at/‖at‖`, `b = b_til − (b_til·a)a`, `c = ct/√(ct'G_s ct)`,
`d = d_til − (d_til'G_s c)c`.

### The projection leaves exactly four redundant coordinates

Those four maps are invariant under: the scale of `at` (hence of `z_a`),
`b_til → b_til + s·a`, the scale of `ct`, and `d_til → d_til + s·c`. So the
likelihood has four directions of exactly zero curvature. Measured as the
nullity of `J'WJ` with no prior, at the highest- and median-`lp__` draws of
three archived runs:

| run | (Kq, Kpsi) | coefficients | nullity, best `lp__` | nullity, median `lp__` |
|---|---|---|---|---|
| `stage3-fisher` | (8,12) | 36 | 4 | 4 |
| `stage4-100` | (8,12) | 36 | 4 | 4 |
| `stage3-k2416` | (16,24) | 76 | 4 | 4 |

Four at every order and every point, independent of the coefficient count, and
the four predicted generators carry Rayleigh quotients of 1e-16 to 4e-10
relative to λ_max — they *are* that null space.

### λ_min of the coefficient Fisher is the prior on a coordinate the data cannot see

Adding the prior precision back, the smallest eigenvalue of the full
coefficient Fisher is carried by **one** of the four — `d_til`'s `G_s`-component
along `c` — and equals the prior precision along it:

| run | point | τ_ψ | λ_min | prior precision along `d_til along c` | \|cos\| with the λ_min eigenvector |
|---|---|---|---|---|---|
| `stage3-fisher` | best | 1.182 | 9.9713e−01 | 1.0021e+00 | 0.9956 |
| `stage3-fisher` | median | 6.187 | 3.2002e−02 | 3.2034e−02 | 1.0000 |
| `stage4-100` | best | 3.213 | 1.7173e−01 | 1.7175e−01 | 1.0000 |
| `stage4-100` | median | 40.483 | 1.3781e−03 | 1.3781e−03 | 1.0000 |

The likelihood contributes exactly nothing there, so the prior sets it alone, and
the prior's scale on `d_til` is `tau_psi`. Hence **λ_min ∝ τ_ψ⁻² and
sqrt(cond) ∝ τ_ψ** — a redundant coordinate, invisible to the data and
discarded by the projection, controls the step size, with its width set by a
sampled hyperparameter that drifts large. At the `stage4-100` median draw τ_ψ
had reached 40.5 and sqrt(cond) 1.8e6.

What removing the four directions buys, as the ratio of sqrt(cond) with them to
without:

| run | (Kq, Kpsi) | best `lp__` | median `lp__` |
|---|---|---|---|
| `stage3-fisher` | (8,12) | 1.00 | **5.62** |
| `stage4-100` | (8,12) | 1.38 | **9.84** |
| `stage3-k2416` | (16,24) | 1.00 | 1.00 |

**Two points were essential and one would have misled.** At the best-`lp__`
draw the ratio is 1.00 and the honest reading would have been "the gauge costs
nothing"; that conclusion was in fact drawn and withdrawn here. τ_ψ happens to
sit near 1 at the best draw of both `(8,12)` runs, which is exactly where the
tau-scaled prior precision is of order the unit scale and the effect vanishes.

### This is not the funnel the non-centring already fixed

`corr(log τ, log‖z‖)`, the recorded diagnostic for the stage-3 funnel, over all
draws:

| run | corr(log τ_q, log‖z_a‖) | corr(log τ_ψ, log‖z_c‖) |
|---|---|---|
| `stage3-fisher` | +0.032 | +0.064 |
| `stage3-dense` | +0.007 | −0.011 |
| `stage4-100` | −0.198 | +0.057 |
| `stage3-k2416` | −0.018 | +0.028 |

All but one below 0.07, so centring the data-pinned `b̃, d̃` and non-centring the
gauge-scaled `ã, c̃` is holding. The coupling found here is a different one: it
is between τ_ψ and the *redundant* component of `d_til`, which non-centring does
not touch because that component is not scaled by `z`.

### The fix is exact and cheap

The component is discarded — `d = d_til − (d_til'G_s c)c` throws it away — so it
need not carry the coefficient prior at all. Decomposing `d_til` into its
`G_s`-orthogonal part plus `s·c` and giving `s` a unit-scale prior pins λ_min at
1 regardless of τ_ψ, at no cost in expressiveness, since nothing downstream reads
`s`. The same applies to `b_til along a` with τ_q, whose prior precision was
measured at 1.26, 3.22, 0.32 and 0.13 across the four points above — smaller
swings, same mechanism. Removing the coordinate outright would be cleaner still
but needs a basis of `c`'s complement that varies with `c`.

This is a change to how Eq. (gauge) is realised, so it belongs to
`model-definition.tex` before it belongs to a program.

### A trap, recorded because it produced a wrong number here

`precondition.pieces` assembles `x_p` from `mu_p` and `u_p`, not from a column
named `x_p`. A point built by reading `x_p` and `x_psi` directly leaves both
amplitudes at zero, so `b + x_p a = b` and `d + x_psi c = d`, the mode vectors
never reach the forward map, and **the whole `z_a` and `z_c` blocks go null**.
That reported a nullity of 18 of 36 and 38 of 76 — in each case exactly
`dim(z_a) + dim(z_c)`, which is the tell. Assert `rms(x_p) > 0` before trusting
any spectrum from this path. Relatedly, `precondition._named` harvests only a
curated name list and silently omits `colr`, `x_s_tot` and `b_0`, so a point
built with it from an older run raises `KeyError` rather than returning a
partial point — the safer failure of the two.

## stage3-k1212-n60-floor — the conditioning cliff is Kpsi, and what is left is multimodality

> **Corrected 2026-11-08, after the entry was first written the same day.** This
> entry originally concluded that what remained after the conditioning cliff was
> a funnel in `(log tau_psi, d_til)`, and proposed non-centring as the fix. That
> was wrong, and the sections below have been rewritten. The posterior is
> multimodal, the barriers are measured, and `tau_psi` labels which mode a chain
> landed in rather than a ridge it slides along. The observations in the tables
> stand; the mechanism they were read as supporting does not. What produced the
> error is recorded under "How the funnel reading went wrong", because the
> pattern is reusable and cost most of a day.

`runs/stage3-k1212-n60-floor`, launched 2026-11-07 09:20, **chain 4 stopped by
hand 2026-11-08 10:47** at 817 rows of 1000; chains 1–3 ran to completion. A
negative result on convergence and a positive one on where the obstruction is
not.

Launched to split the two failures of `stage3-k1216-n60-floor` apart by changing
one thing, `Kpsi` 16 → 12, with `Kq` left at 12. The truncation scan in that
entry put the cliff between `Kpsi` 12 and 14 from the Fisher spectrum alone;
this run tests the prediction in a sampler.

**The prediction holds.** The conditioning failure is `Kpsi`, and at `Kpsi = 12`
it is gone: step sizes and treedepths return to the converged reference's range.
The run still does not converge, for a different and now identifiable reason.

### Inputs

| | |
|---|---|
| commit | `199978a` plus the uncommitted `k_a`/`k_c` additions; the program itself is `199978a`'s |
| program | `stan/hier_common.stan` — no departure modes, `sigma_b` sampled |
| `free_mean` | 1, so ⟨Δ⟩ is sampled |
| orders | `Kq = 12`, `Kpsi = 12` |
| selection | `stage3_sample(n_per_fifth=12)` — 60 objects, 3682 epochs |
| `data.json` | sha256 `8faa8d006936e746…` |
| `inits-1..4.json` | sha256 `401588fdd2b40595…` — all four identical, so again no dispersed-start check |
| initialised from | `runs/stage3-k88-n60-floor`; 15 blocks carried, lp at the transferred point −32426.30 |
| gradient check | worst of 20 coordinates 1.9e-07 |
| sampled coordinates | 268 |
| wall time | 37395 s, 58395 s, 41075 s for chains 1–3; chain 4 stopped at ~26 h |

### Sampler end state, each chain's own post-warmup draws

| | ch 1 | ch 3 | ch 2 | ch 4 |
|---|---|---|---|---|
| draws past warmup | 500 | 500 | 500 | 316 |
| mean `treedepth__` | 9.00 | 9.00 | 10.00 | **12.00** |
| fraction at the cap of 12 | 0.00 | 0.00 | 0.00 | **1.00** |
| mean `n_leapfrog__` | 511 | 520 | 1023 | 4095 |
| frozen step size | 8.56e-3 | 6.60e-3 | 4.64e-3 | 1.05e-3 |
| divergences post-warmup | 0 | 0 | 0 | 0 |
| `tau_q` | 0.566 | 0.515 | 0.665 | 0.632 |
| `tau_psi` | 1.853 | 1.979 | 4.319 | 5.646 |
| `alpha_q_m` | 0.675 | 0.753 | 0.370 | 0.426 |
| `alpha_psi_m` | 0.236 | 0.212 | 0.226 | 0.245 |
| `mu_delta` | −2.2246 | −2.7507 | −2.2232 | −2.4858 |
| `c_bar` | +0.0508 | +0.0681 | +0.0204 | +0.0583 |
| `sigma_g` | 0.0474 | 0.0473 | 0.0492 | 0.0492 |
| `sigma_r` | 0.0615 | 0.0612 | 0.0589 | 0.0582 |
| χ²/epoch over V | 1.1059 | 1.1051 | 1.1051 | 1.1057 |
| `lp__` | −32223.0 | −32213.7 | −32207.1 | **−32195.5** |

Columns are ordered by `tau_psi`, not by chain id, because that is the order
everything else follows.

Against `stage3-k88-n60-floor`: depth 8.00, step size 1.4–1.6e-2. Against
`stage3-k1216-n60-floor`: depth ≥ 11 in every post-warmup draw on three chains
of four, step size 1.5–1.7e-3. **Three of the four chains here are at depth 9–10
with step sizes within a factor of two of the converged reference.** `tau_psi`
is 1.9–5.6 against that run's 4.6–8.1, and `alpha_psi_m` no longer presses zero:
0.212–0.245, with the reference at 0.389–0.428 and nothing near the boundary.

### It does not converge, and it is not converging

| | over the common 316-draw window |
|---|---|
| sampled coordinates | 616 |
| max split R-hat | **4.368** |
| above 1.01 | 521 |
| above 1.05 | 467 |

The project's bar is 1.01 and `diagnose()`'s is 1.05; this misses both by a
factor, not a margin. **More draws would not have helped, and that is measured
rather than asserted.** Split R-hat's limit is set by the between-chain term `B`
against the within-chain term `W`; on the worst coordinate `d.5`:

| draws used | W | B | **B/W** | split R-hat |
|---|---|---|---|---|
| 50 | 0.167 | 232 | **1392** | 7.53 |
| 100 | 0.169 | 520 | **3071** | 7.90 |
| 150 | 0.313 | 704 | **2250** | 5.57 |
| 200 | 0.464 | 792 | **1707** | 4.25 |
| 250 | 0.393 | 920 | **2344** | 4.44 |
| 296 | 0.316 | 1042 | **3301** | 4.83 |

`B/W` ends higher than it began. The fall in R-hat from 7.90 to 4.25 between 100
and 200 draws is `W` growing, not `B` shrinking, and `W` has since fallen back
while `B` kept climbing. With the offsets real and fixed the asymptote is
√(1 + B/W) ≈ 57. The count above 1.01 sat at 521–528 of 616 in every window from
50 draws to 316 and never moved. **A falling R-hat is not evidence of
convergence unless B/W is what is falling.**

### At least three modes, with the barriers measured

**The four chains are in separate modes, and `tau_psi` labels which one.** Walk
the full parameter vector in a straight line from one chain's median draw to
another's and evaluate Eq. (18)'s log-likelihood along it. The test calibrates
itself: for two independent draws from the *same* mode in `D` dimensions the
midpoint sits **above** both endpoints by `D/4`, and this run samples `D = 268`
coordinates, so the same-mode prediction is +67 nats.

| path | midpoint vs the better endpoint | same-mode prediction | implied barrier |
|---|---|---|---|
| ch1 ↔ ch3 | **+72.7** | +67 | none — one mode |
| ch3 ↔ ch2 | −27.4 | +67 | ≈ 94 |
| ch1 ↔ ch2 | −107.5 | +67 | ≈ 174 |
| ch2 ↔ ch4 | −526.8 | +67 | ≈ 594 |
| ch3 ↔ ch4 | −1230.4 | +67 | ≈ 1297 |
| ch1 ↔ ch4 | −1366.6 | +67 | ≈ 1434 |

`ch1 ↔ ch3` reproduces the same-mode prediction to 5.7 nats of 67, which is what
makes the contrast trustworthy; the lower four rows are nothing like it. A
straight-line path bounds the barrier from above, so a lower curved route between
modes is not excluded; the contrast against the calibrated same-mode row is the
result, not the absolute heights.

**The control: the same walk on `stage3-k88-n60-floor`, which is converged.**
All six paths between its four chains give a midpoint gain of +54.7 to +65.7
nats against its own `D/4 = 65` (260 sampled coordinates), with `tau_psi`
spanning 2.174–2.246 and `‖d_til‖` 6.72–6.79. That run is unimodal and the walk
says so, reproducing what R-hat 1.014 already established. **So the test is
calibrated against a known answer, and the multimodality is a `Kpsi = 12`
phenomenon rather than a property of sixty objects.**

**The barriers are a property of the population, not of this sample**, which is
what decides whether more objects would remove them. Decomposing each midpoint
dip object by object — a noise barrier being a sum of mixed signs and a real one
a vote most objects join:

| path | objects voting against the midpoint | dip / its own sampling error |
|---|---|---|
| ch1 → ch3 | 12 of 60 | **−5.7** (midpoint better: one mode) |
| ch3 → ch2 | 36 of 60 | 1.2 |
| ch1 → ch2 | 36 of 60 | 3.5 |
| ch2 → ch4 | 47 of 60 | 5.6 |
| ch3 → ch4 | 53 of 60 | 6.0 |
| ch1 → ch4 | 54 of 60 | 5.8 |

So the separations that survive are `{ch1, ch3}` against `{ch4}` and `{ch2}`
against `{ch4}` — **two modes, not three**, since ch2 against ch3 is 1.2 standard
errors. `‖d_til‖` is 5.83, 6.49, 9.00 and 14.72 across the four chains and chain
4 holds the highest likelihood (−32294.86 against −32324.41). Because the
log-likelihood is a sum over objects, a barrier the objects agree on scales with
the sample: more supernovae sharpen these modes and deepen the barriers rather
than dissolving them.

This is the two-modes-in-`tau_psi` phenomenon already on record for the
mode-bearing programs, now at `hier_common` and with barrier heights attached.

The observations below stand as observations, and they are what the mode labels
look like from inside each chain. The turning coefficients are an affine function
of `tau_psi` across the four chains pooled — which is a between-mode regression,
not a ridge:

| coefficient | slope against `tau_psi` | corr |
|---|---|---|
| `d.1` | −0.814 | −0.845 |
| `d.3` | +0.980 | +0.823 |
| `d.5` | −1.110 | −0.834 |
| `d.8` | +0.656 | +0.844 |
| `d.11` | −0.189 | −0.738 |

| chain | `tau_psi` mean | `tau_psi` range | step size | depth |
|---|---|---|---|---|
| 1 | 1.60 | [0.68, 4.31] | 8.6e-3 | 9.00 |
| 3 | 2.06 | [0.98, 5.27] | 6.6e-3 | 9.00 |
| 2 | 3.19 | [1.05, 8.04] | 4.6e-3 | 10.00 |
| 4 | 5.69 | [2.82, 14.54] | 1.1e-3 | 12.00 |

The signs alternate with index, the same alternating-odd-index pattern the Fisher
λ_min eigenvector loaded on at `Kpsi = 16` (d7 −0.520, d5 +0.405, d3 −0.402,
d9 +0.343, d11 −0.287). The modes therefore differ from one another along that
pattern, and the chain in the highest-`tau_psi` mode is also the stiffest: chain
4 was pinned at the cap for every one of its 316 sampling draws. What is *not*
established by this table is any motion along a ridge — pooling draws from four
modes produces a between-mode regression whatever the within-mode geometry is.

**Chain 1 is not stationary within itself**, so within its own mode it had also
not finished:

| within-chain split on `d.5` | first half | second half | split R-hat |
|---|---|---|---|
| ch 1 | −0.164 | −2.494 | **4.28** |
| ch 2 | −3.801 | −4.291 | 1.12 |
| ch 3 | −2.780 | −2.916 | 1.01 |
| ch 4 | −8.092 | −7.576 | 1.21 |

Chains 2–4 are locally stationary, each inside its own mode, which is the worst
case for R-hat: well-behaved locally and wrong globally, so no additional
sampling repairs it.

### The prior is not the mechanism, and that is measured

The obvious suspect is `stan/hier_common.stan:275`,

```stan
d_til ~ normal(0, tau_psi * pow(kpsi / kbar_psi, -(alpha_psi_m + 0.5)));
```

whose precision goes as `tau_psi⁻²`, and the program's own comment at line 95
records `corr(log tau_psi, log‖d‖) = +0.82` as this posterior's worst direction.
**It is not the cause here.** Solve the generalised eigenproblem
`F_like u = λ P u` on the `d_til` block, with `P` the diagonal prior precision
and `F_like` the Fisher with the per-object block profiled out:

| run | draw | `tau_psi` | generalised spectrum, ascending | directions with λ < 1 |
|---|---|---|---|---|
| `stage3-k1212-n60-floor` | best lp | 3.169 | 1.6e5 … 2.8e8 | **0 of 11** |
| `stage3-k1212-n60-floor` | median lp | 2.563 | 1.1e5 … 1.1e8 | **0 of 11** |
| `stage3-k88-n60-floor` | best lp | 2.058 | 1.3e5 … 2.7e7 | **0 of 7** |
| `stage3-k88-n60-floor` | median lp | 1.692 | 5.6e4 … 6.4e7 | **0 of 7** |

The likelihood dominates the prior by five to eight orders of magnitude in every
direction of the block, at both orders and both draws. **Non-centring is the
remedy for prior-dominated coordinates and there are none**, so it is not the
fix, and neither is any reparameterisation of Eq. (coefprior)'s scale. The
`λ_min ∝ τ_ψ⁻²` identity of `gauge-nullity` applies to a gauge-redundant
direction the data cannot see; nothing here is that direction.

Nor is the width a per-object or traversal degeneracy. Widths along the leading
principal direction of the chain's own `d_til` motion, at the median-lp draw:

| | sd | × conditional | share of the log gap |
|---|---|---|---|
| conditioning on everything | 0.0107 | 1.0 | 0% |
| profiling the 240 per-object parameters | 0.0225 | 2.1 | 14% |
| profiling `b_til` | 0.0131 | 1.2 | 4% |
| profiling `b_til` **and** per-object | 0.0624 | 5.9 | 33% |
| sampled spread within chain 1 | 2.4554 | 230.4 | 100% |

Everything available to a linear analysis accounts for a factor of 6 out of 230.
The residue is the between-mode separation, which no local matrix can see.

### How the funnel reading went wrong

Recorded because it cost most of a day and the pattern recurs.

1. Four chains at four `tau_psi` values with overlapping ranges and a clean
   affine `d_til`–`tau_psi` relation (|corr| 0.74–0.85) read as one ridge. Pooled
   across modes, that relation is what multimodality always looks like.
2. A conditional curvature was then compared against a marginal spread: the
   Fisher width along the direction is 0.021 and the chain's spread is 2.46, and
   the gap was attributed first to the hyperparameters, then to the per-object
   block. Measured, those account for 1.5× and 2.1×.
3. The step that looked decisive and was not: chain 1 travels 2.46 along the
   direction and loses 16.3 nats where a Gaussian of width 0.021 would lose 6707,
   which was read as the Fisher overstating the curvature 400-fold. A direct
   1-D profile of the log density settles it the other way — half-width 0.043 at
   a 0.5-nat drop, in agreement with the Fisher to a factor of two. **A chain can
   span many conditional widths without any curvature being wrong, if it is not
   sampling one mode.**

The general lesson is the cheap test that ends it. **Before explaining a wide
marginal, walk between the chains.** The midpoint of two same-mode draws is
`D/4` nats *above* both endpoints, so the walk is self-calibrating and costs one
forward map per step. It should have been the first measurement, not the sixth.

### What differs between the modes, and what does not

| quantity | corr with `tau_psi` |
|---|---|
| `sigma_r` | **−0.619** |
| `sigma_g` | +0.419 |
| `mu_delta` | −0.115 |
| `c_bar` | +0.026 |
| χ²/epoch over V | +0.031 |

**`sigma_r` varies between the modes, so the reduction it shows cannot be quoted
as a result.** `sigma_r` fell to 0.0582–0.0615 against the reference's
0.0678–0.0679, and `sigma_g` moved only 0.0470 → 0.0473–0.0492, at χ²/epoch over
V essentially unchanged (1.1051–1.1059 against 1.1084–1.1101) — so the restored
orders do appear to absorb structure the per-band error was carrying in r, with
the residuals shrinking alongside V rather than V alone. But it tracks the mode
label at −0.62, so what is measured is that the modes disagree about `sigma_r`,
not that `Kpsi = 12` lowers it. It must be re-read once one mode is known to be
the right one.

**`mu_delta` and `c_bar` do not distinguish the modes**, at −0.115 and +0.026, so
their behaviour is a separate matter. All four chains put `mu_delta` at −2.22 to −2.75
against the reference's +0.408 to +0.484 — a genuine order effect and not a ridge
artefact. ⟨Δ⟩ was already recorded as disagreeing in sign between the thirty- and
sixty-object samples; the orders now disagree on it too, and by far more than any
quoted precision. This is unexplained.

### Two corrections to earlier entries

**Warmup divergence counts do not separate these runs and must not be read as a
failure signature.** The `stage3-k1216-n60-floor` entry quotes 75/67/79/76 as
part of its end state. Over the same first 500 rows:

| run | warmup divergences |
|---|---|
| `stage3-k88-n60-floor`, converged | 65 68 71 63 |
| `stage3-k1216-n60-floor`, failed | 75 67 79 76 |
| `stage3-k1212-n60-floor`, failed | 61 69 64 68 |

All three are the same to within their own scatter. The discriminating columns
are step size, treedepth and leapfrog count.

**`tau_psi = 1.267` and `mu_delta = +0.4345`, quoted from
`stage3-k88-n60-floor` throughout this ledger, are best-draw values, not
posterior means.** The post-warmup means are `tau_psi` 2.17–2.25 and `mu_delta`
+0.408–0.484. Comparing a mean to a best draw made `(12,12)`'s `tau_psi` look
inflated when it is sitting on the reference. Means go with means.

### Known gaps

- **Nothing here is converged, so no value above is a measurement of the
  posterior.** The tables describe where four chains sat, which is what localises
  the obstruction and nothing more.
- **No remedy has been tested, and reparameterisation is ruled out rather than
  untried.** With 0 of 11 directions prior-dominated there is nothing for
  non-centring to act on. Multimodality is addressed by where chains start and
  how they are tempered, not by coordinates: dispersed starts,
  `scripts/ref_both_modes.py`, or starting every chain in the highest-likelihood
  mode and checking it is stable there.
- **Which mode is right is not known, and the margin is not significant.**
  Chain 4's mode leads chain 1's by 29.6 nats of Eq. (18) log-likelihood, but
  only 33 of 60 objects favour it and the per-object scatter is 3.17 against a
  mean of 0.493, so the margin is **1.2** times its own sampling error over the
  object population. Against chain 3 it is 1.0, against chain 2 it is 0.2. The
  lead comes from a handful of objects, not from the sample.

  An earlier version of this entry said the full 30150-epoch training sample
  would resolve this "decisively — about 240 nats if it holds". That scaled the
  margin and not its error. Both go as the sample, the margin as `N` and the
  error as `sqrt(N)`, so 442 objects buys `sqrt(442/60) = 2.7` in significance:
  1.2 becomes about **3.3 σ**. Worth having, not decisive.
- **The barriers, unlike the margin, are a population property**, so more data
  will not remove the modes — it sharpens them. 47 to 54 objects of 60 vote
  against the midpoint on the three paths involving chain 4, at 5.6 to 6.0
  standard errors, and the barrier scales as the sample: ch1↔ch4 would be near
  9300 nats at 442 objects against 1263 here. The modes cannot merge either,
  sitting at `‖d_til‖` 5.83 and 14.72.
- **The defensible mode count is two, not three.** The ch2↔ch3 separation is
  1.2 standard errors on 36 of 60 objects, so those two may be one mode with a
  shallow dip between the draws; ch1↔ch2 is 3.5. What survives is {ch1, ch3}
  against {ch4}, and {ch2} against {ch4}. And it remains a lower bound: four
  chains from byte-identical inits is the worst possible design for counting
  modes.
- **Chain 4 contributes 316 draws against the others' 500**, so the common
  window used for R-hat is 316 and the per-chain table mixes two lengths.
- **All four inits are identical**, so the dispersed-start check that
  `model-definition.tex` §5 requires of an initialisation from \salt\ is still
  not done, in this run or any other on record.
- **`Kpsi` between 12 and 14 is bracketed but not resolved.** The cliff was
  located by a Fisher scan at three points and confirmed at 12 by this run; 13
  and 14 are untested in a sampler.
- Evaluation remains in-sample; `src/dgsn/eval/` is empty.

## stage3-abox-k88-n60 — the exponential prior on the decay exponents, removed

`runs/stage3-abox-k88-n60`, finished 2026-11-12. **The reference fit of
`model-v1.2`.**

A controlled comparison of one prior against `stage3-k88-n60-floor`, the
`model-v1.1` reference. The whole of the difference, with comments stripped:

```
<   real<lower=0> alpha_q_m;          >   real<lower=0, upper=8> alpha_q_m;
<   real<lower=0> alpha_psi_m;        >   real<lower=0, upper=8> alpha_psi_m;
<   alpha_q_m ~ exponential(1);       (removed)
<   alpha_psi_m ~ exponential(1);     (removed)
```

Everything else — the forward map, the gauge, the likelihood, the quadrature,
the data, the orders, the seed — is `hier_common.stan`'s.

**Convergence is not lost.** On every diagnostic it matches or beats v1.1, and
no physical quantity moves. What moves is α itself, which is the point: the
`exponential(1)` was pulling both decay exponents down by most of a posterior
standard deviation and doing nothing else.

### Inputs

| | |
|---|---|
| commit | `7fa2216`, where the change entered as `stan/hier_abox.stan`; folded into `stan/hier_common.stan` immediately after, so that `--no-trav-mode` gives model-v1.2 and the old program is reached by the `model-v1.1` tag |
| program | `stan/hier_abox.stan` as run; now `stan/hier_common.stan`, whose code differs from the `model-v1.1` tag's by the four lines above and nothing else |
| `free_mean` | 1, so ⟨Δ⟩ is sampled |
| orders | `Kq = 8`, `Kpsi = 8` |
| selection | `stage3_sample(n_per_fifth=12)` — 60 objects, 3682 epochs |
| `data.json` | sha256 `3a017877acbd860d…` against v1.1's `00e0cd556ecfd05f…`. **Every key the program reads is bit-identical**; the dict carries two extra keys, `Ka` and `Kc`, which the program does not declare and Stan therefore ignores. Checked key by key, not by sha. |
| `inits-1..4.json` | sha256 `312a7d110d5b9e3d…` — all four identical |
| initialised from | `runs/stage3-k88-n60-floor`, its best-`lp__` draw |
| metric | **none transferred, deliberately** — v1.1 adapted its own from a unit start, and a warm metric would pre-solve the conditioning this run exists to test |
| CmdStan | 2.39.0 |

### Sampler

```
seed 1, 4 chains, 500 warmup + 500 sampling, max_treedepth 12,
save_warmup on, diag_e, adapt_delta at CmdStan's default
```

```bash
# as run, when the change lived in its own program:
python scripts/stage3_hier.py --no-trav-mode --alpha-box --k-psi 8 --k-q 8 \
  --n-per-fifth 12 --init-from runs/stage3-k88-n60-floor \
  --iter 500 --max-treedepth 12 --save-warmup --out runs/stage3-abox-k88-n60

# the same run at HEAD, --alpha-box having been retired with the duplicate:
python scripts/stage3_hier.py --no-trav-mode --k-psi 8 --k-q 8 \
  --n-per-fifth 12 --init-from runs/stage3-k88-n60-floor \
  --iter 500 --max-treedepth 12 --save-warmup --out runs/stage3-abox-k88-n60
```

### Diagnostics, against v1.1

| | v1.1 (`Exp(1)`) | this run (flat box) |
|---|---|---|
| max split R̂, 776 coordinates | 1.014 | **1.011** |
| coordinates > 1.05 | 0 | 0 |
| coordinates > 1.01 | 60 | **2** |
| post-warmup divergences | 0, 0, 0, 0 | **0, 0, 0, 0** |
| treedepth | 8 on every draw | 8; two chains touch 9 |
| median post-warmup leapfrog | 255 | 255 |
| step size | 1.41–1.60e-2 | 1.27–1.37e-2 |
| min ESS over 24 globals | 375 (`b_0`) | 313 (`b_0`) |
| median ESS over those 24 | 985 | 781 |
| gradients, 4 chains incl. warmup | 2,368,026 | 2,399,043 |

Wall was 53,699–55,090 s per chain against v1.1's 24,571 s, but **that is
contention and not the model**: three runs shared twelve cores, load average 23.
The gradient counts above are the cost that belongs to the model, and they agree
to 1.3%.

### What it measured

α rises, and nothing else moves. Shifts in units of the pooled posterior sd:

| | v1.1 | this run | shift |
|---|---|---|---|
| α_q | 1.665 ± 0.668 | **2.198 ± 0.725** | **+0.77σ** |
| α_ψ | 0.910 ± 0.306 | **1.047 ± 0.366** | **+0.41σ** |
| κ₀ | 52.710 ± 5.384 | 52.676 ± 5.540 | −0.01σ |
| b₀ | −7.768 ± 0.142 | −7.767 ± 0.133 | +0.01σ |
| σ_g | 0.0469 ± 0.0014 | 0.0469 ± 0.0014 | −0.01σ |
| σ_r | 0.0679 ± 0.0014 | 0.0679 ± 0.0015 | −0.01σ |
| b₁…b₇ | — | — | all within 0.16σ |
| d₁…d₇ | — | — | all within 0.06σ |

The box never binds: the largest draw is α_q = 5.00 and α_ψ = 3.12 against a
ceiling of 8.5, about ten posterior sd above either.

Read the sign carefully. `exponential(1)` on `alpha_m` decays with `alpha_m`,
so it pulls α toward its ½ floor — toward *rougher* q₀ and Ψ₀. The data prefer
smoother ones. The prior was not a tail bound; it was a preference, and it acted.

### What it does not establish

- **It starts from v1.1's solution.** `--init-from` takes the best-`lp__` draw,
  so this shows the solution *holds* without the prior, not that it is *found*
  from a cold, SALT2-placed start. A cold run is the stronger test and has not
  been done.
- **ESS is ~20% worse**, consistently across every b coefficient, at equal
  gradient count. The direction is uniform enough that it is unlikely to be
  noise: α's posterior is wider, and α sets the b prior widths. Reported as a
  cost, not explained away.
- Everything in-sample, as ever. `src/dgsn/eval/` is still empty.

## stage3-k88-n255-warm — the sample widened to 255, with a transferred metric

`runs/stage3-k88-n255-warm`, finished 2026-11-12.

Answers whether the sample can be widened: **it can.** `hier_common.stan`, 255
objects, 19027 epochs, `data.json` sha256 `61c7d163d6b5d7eb…`, initialised and
preconditioned from `runs/stage3-k88-n60-floor` with `--init-buffer 15`.

```bash
python scripts/stage3_hier.py --no-trav-mode --k-psi 8 --k-q 8 --n-per-fifth 51 \
  --init-from runs/stage3-k88-n60-floor --metric-from runs/stage3-k88-n60-floor \
  --init-buffer 15 --iter 500 --max-treedepth 12 --save-warmup \
  --out runs/stage3-k88-n255-warm
```

| | this run (255) | v1.1 (60) |
|---|---|---|
| max split R̂ | **1.017** (3117 coords) | 1.014 (776) |
| coordinates > 1.05 | 0 | 0 |
| coordinates > 1.01 | 8 | 60 |
| post-warmup divergences | **0, 0, 0, 0** | 0, 0, 0, 0 |
| treedepth | 9 on every draw | 8 on every draw |
| median post-warmup leapfrog | 511 | 255 |
| step size | 8.49–9.55e-3 | 1.41–1.60e-2 |

One doubling more than v1.1 for 4.25× the objects and 5.17× the epochs. The
metric transfer is what made it affordable: the cold attempt sat pinned at the
treedepth cap, and at matched iteration 11 the warm start ran at a step size
179× larger. Wall 230,195–236,415 s per chain, under contention.

## stage3-flat-k88-n60 — the hierarchy removed and the amplitudes reparameterised

`runs/stage3-flat-k88-n60`, finished 2026-11-12. **Does not clear the bar.**

`stan/hier_flat.stan`, renamed `stan/hier_full.stan` after this run (the run's CSVs
keep the old basename): τ and α dropped entirely, `at = z_a` and `ct = z_c` with
only the gauge-fixing normals left on b and d, *and* `x_p`, `x_psi` moved from
`mu + u` to flat boxed ±20. 390 sampled dimensions.

| | this run | v1.1 | `hier_modes` (stopped) |
|---|---|---|---|
| max split R̂ | **1.077** (8433 coords) | 1.014 | 1.234 on the traversal globals, run stopped before completion |
| coordinates > 1.05 | **49** | 0 | — |
| post-warmup divergences | 6, 7, 4, 2 | 0 | 17–51 |
| treedepth | 7–11, median 11 | 8 | 12, 83–90% at cap |
| median post-warmup leapfrog | 2047 | 255 | 2748–3913 |
| step size | 2.62–3.24e-3 | 1.41–1.60e-2 | — |

Much better than `hier_modes` and still not converged. The residue is `d_dev`
(19 of 60 objects above 1.05, worst 1.077) and `delta_t` (7 of 60); the
traversal globals `a` and `b` sit at 1.069 and 1.070.

**It changed two things at once** — the hierarchy removal and a 120-coordinate
amplitude reparameterisation — so it attributes neither. `stage3-abox-k88-n60`
now partly resolves this from the other side: removing the α prior *alone* costs
nothing geometrically (depth 8, unchanged), so this run's depth-11 geometry is
**not** attributable to the hierarchy removal.

One correction worth keeping. At 67 post-warmup draws `x_p.10` read R̂ 1.400 and
looked like a per-object failure; at 500 draws it reads 1.029. The 67-draw floor,
measured on runs known to converge — max 1.107 on n255, 1.056 on v1.1 — is what
said not to interpret it. Measure the floor first.
