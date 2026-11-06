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
