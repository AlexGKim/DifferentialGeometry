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
