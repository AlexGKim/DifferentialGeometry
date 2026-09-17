# DifferentialGeometry

Empirical modelling of Type Ia supernova `g`, `r` light curves as curves in
magnitude space. Each supernova traces a curve; the curve's shape lives in its
tangent angle as a function of arclength, and the timing in how fast the
supernova walks along it, so that shape and timing are separated by
construction rather than entangled in one parameter as in SALT2. The model is
defined in [`docs/tex/model-definition.tex`](docs/tex/model-definition.tex),
which is the authority for everything here.

## Documents

| | |
|---|---|
| [`docs/tex/model-definition.tex`](docs/tex/model-definition.tex) | The model, the sample, the likelihood, the inference and its staged bring-up |
| [`docs/model-note-backlog.md`](docs/model-note-backlog.md) | Reasoning that would not fit inline, and retired decisions |
| [`CLAUDE.md`](CLAUDE.md) | The data, the traps, the tooling |

Build the note:

```bash
cd docs/tex && latexmk -pdf model-definition
```

## Environment

The forward map, the likelihood and the sampler are Stan (`stan/`); Python
builds data and reads output, and evaluates nothing of the model itself.

```bash
conda create -n dgsn -c conda-forge python=3.12 numpy scipy pandas matplotlib pytest cmdstan cmdstanpy sncosmo
conda activate dgsn
pip install -e .
```

CmdStan is found through the environment's activation (`$CMDSTAN`). The first
use of each Stan program compiles it, about ten seconds each. `sncosmo` fetches
the SALT2 and Hsiao templates on first use.

Data: the ZTF SN Ia DR2 `lite` release unpacked at `data/raw/ztfsniadr2_lite/`,
or wherever `$ZTFCOSMODIR` points. Everything under `data/` is gitignored.

## Layout

```
stan/dgsn_functions.stan   the model, once; #included by the programs
stan/forward.stan          the forward map at fixed parameters (fixed_param)
stan/checks.stan           closed-form checks of the quadrature and clamps (tests only)
stan/single_sn.stan        one supernova, six parameters, globals as data
src/dgsn/stan/             constants of the note (data.py); compile-once drivers (run.py)
src/dgsn/salt/             SALT2 and Hsiao synthesis (tracks.py); the global functions from SALT2 (globals.py)
src/dgsn/data/             selection and counts (sample.py); the likelihood epoch set (epochs.py)
scripts/                   the stage runs below
tests/                     stage 1, as pytest
```

## Running the stages

The stages are those of the note's Section 5.2. Each tests what the one before
it could not.

### Stage 1 — the forward map against closed forms

```bash
pytest -q -s
```

About three minutes; `-s` shows the measured errors. What runs:

- `tests/test_checks.py` — the quadrature on a circle (constant `psi'`), a unit
  rate (`s = p`), and a smooth integrand with a known antiderivative; the
  clamp, taper and softplus at their edges. Each at two grid resolutions.
- `tests/test_forward.py` — `s(0) = 0`, `gamma(0) = 0`, `psi(0) = -pi/2`
  exactly; `x_s` scales arclength by `e^{x_s}`; both bases orthonormal, the
  turning basis zero at `s = 0` and triangular; `psi'(0)`; flux against
  magnitude.
- `tests/test_single_sn.py` — the dilation and translation tests as a gradient
  that vanishes at the exact parameters and a fit that returns there; Stan's
  gradient against finite differences of its own log density.

Pass: every closed form at the tolerance recorded in the test (machine
precision, or the measured quadrature floor where the test says so), and the
two resolutions agreeing.

The floors of Section 3, which the note assigns to this stage:

```bash
python scripts/floors.py
```

About one minute. Pushes noise-free Hsiao light curves through the truncated
model with SALT2-fitted globals — at `z = 0.05` against `z = 0` for the
K-correction, with F99 dust at `E(B-V) = 0.3` against a constant translation
for the extinction drift — and reports what each does to the six parameters
and to the per-band residual.

### Regenerating what Section 2.3 and 5.1 quote

```bash
python scripts/salt_orders.py            # add --quick for a coarser track grid
```

About ten minutes. The truncation orders a SALT2 track asks for at each `x1`,
with the traversal in `v`, and the quadrature accuracy against panels doubled.
The track fits carry a ridge of the form of the note's eq. (coefprior)
(`RIDGE_*` in `dgsn/salt/globals.py`); the two grid resolutions must agree for
a cell to count. `K_q = 64` at `m = 10` and `x1 = -2` at any order do not: the
first is `exp` of a degree-63 polynomial exploiting quadrature error, the
second a track whose early rise the turning fit does not reproduce (0.19 mag
at every `K_psi`).

### Stage 2 — one supernova, six parameters

```bash
python scripts/stage2_single.py --check-only
python scripts/stage2_single.py [--name ZTF18...] [--sigma 0.05]
```

`--check-only` (seconds): the sampler program's fluxes against `forward.stan`
for the same inputs, and its gradient against finite differences at the
SALT2-placed start. The fit (a minute or two) runs four chains as four jobs,
with the global functions from the SALT2 fits at the sample-median `x1` and
the population hyperparameters as data, and reports `t_max - t0`, `m_peak`
and `colr` against the archive, the residuals by phase, and the interval and
turning checks of the note. Output under `runs/stage2/<name>/`.

Pass: the chains agree with no divergences, and `t_max - t0` is a fraction of
a day.

### Stages 3–5

Not yet implemented: the hierarchy on a few tens of supernovae; the full
training sample; the grid over truncation orders.

## Regenerating the sample counts

```bash
python -m dgsn.data.sample --tex
```

Writes `docs/tex/sample-counts.tex`, which the note `\input`s. Regenerate
rather than quote from memory.

## Status

Stages 1 and 2 implemented and passing (2026-09-17). The note's numbers that
came from a script outside the repository — the orders in `v`, the quadrature
accuracy in `v`, the floors — are now produced by `scripts/`.
