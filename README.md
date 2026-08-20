# DifferentialGeometry

Empirical modelling of Type Ia supernova multi-band light curves using
differential geometry.

A supernova observed in *n* bands traces a **curve** through band-magnitude
space. Parameterised by arclength, that curve is determined up to a rigid
motion by its curvature invariants. This project models `kappa(s)` with a
neural network and the time-to-arclength map `s(t)` separately, so that
**shape** and **timing** are cleanly separated — unlike SALT2, whose `x1`
entangles them.

The question under test: does this separation describe Type Ia diversity with
fewer parameters than SALT2?

## Documents

| | |
|---|---|
| [`docs/tex/model-definition.tex`](docs/tex/model-definition.tex) | The model: geometry, likelihood, degeneracies, validation plan |
| [`docs/tex/salt2-distillation.tex`](docs/tex/salt2-distillation.tex) | SALT2 (Guy et al. 2007) distilled to what this analysis needs |
| [`CLAUDE.md`](CLAUDE.md) | Design decisions and working context |

Build the science documents:

```bash
cd docs/tex
pdflatex model-definition && bibtex model-definition && pdflatex model-definition && pdflatex model-definition
```

## Install

```bash
pip install -e ".[dev,data,benchmark]"
```

## Data

Light curves come from ZTF SN Ia DR2 via
[`ztfcosmo`](https://github.com/ZwickyTransientFacility/ztfcosmo), which needs
no credentials — it reads remotely from `ztfcosmo.in2p3.fr`, or from a local
copy:

```bash
export ZTFCOSMODIR=/path/to/ztfsniadr2
```

The sample is cut to `z < 0.05`, where K-corrections are small enough that a
two-band model with no SED remains self-consistent. Two samples result:

| | SNe | ambient space | invariants |
|---|---|---|---|
| Primary | 599 | `R^2` (`g`,`r`) | `kappa(s)` |
| Torsion subsample | 177 | `R^3` (`g`,`r`,`i`) | `kappa(s)`, `tau(s)` |

Requiring at least five epochs per band with SNR > 5 and `flag & 31 == 0`
inside the rest-frame phase window `[-15, +40]` d.

## Status

Design complete; implementation not yet started.
