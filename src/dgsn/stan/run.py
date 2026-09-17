"""Compile once, then drive the programs under ``stan/``.

Four entry points, one per way the note uses Stan: ``fixed_param`` for the
forward map at fixed parameters, ``optimize`` for the stage-1 fits,
``log_prob`` for gradients, ``sample_four`` for the four chains of Section 5,
run as four separate jobs.
"""

from __future__ import annotations

from pathlib import Path

import logging

import numpy as np
from cmdstanpy import CmdStanMCMC, CmdStanModel, from_csv

# A restarted L-BFGS that cannot move at its first iteration is at the
# optimum, and cmdstanpy logs that as an error; real failures still raise.
logging.getLogger("cmdstanpy").setLevel(logging.CRITICAL)

STAN_DIR = Path(__file__).resolve().parents[3] / "stan"
_models: dict[str, CmdStanModel] = {}


def model(name: str) -> CmdStanModel:
    """The compiled program ``stan/<name>.stan``; compiled on first use."""
    if name not in _models:
        _models[name] = CmdStanModel(
            stan_file=STAN_DIR / f"{name}.stan",
            stanc_options={"include-paths": [str(STAN_DIR)]},
        )
    return _models[name]


def fixed_param(name: str, data: dict) -> dict[str, np.ndarray]:
    """One draw of ``generated quantities`` at the data's parameters.

    Written at 18 significant figures, the most CmdStan allows, so that the
    closed-form comparisons are limited by arithmetic and not by the CSV.
    """
    fit = model(name).sample(
        data=data, fixed_param=True, chains=1, iter_warmup=0, iter_sampling=1, adapt_engaged=False,
        sig_figs=18, show_progress=False, show_console=False,
    )
    return {k: np.asarray(fit.stan_variable(k))[0] for k in fit.metadata.stan_vars}


def optimize(name: str, data: dict, inits: dict, polish: bool = True, **kw) -> dict[str, float]:
    """Mode of the log density; the fitter of the stage-1 tests and floors.

    L-BFGS, restarted from its own result while the restart still moves it:
    a restart discards the Hessian approximation, which is what lets it
    resolve a direction the light curve is nearly flat along (x_psi).
    Newton is not used --- its step along that direction is enormous and it
    diverges to NaN.  Written back at 18 significant figures: CmdStan's
    default of 6 quantises the mode."""
    opts = dict(tol_rel_grad=1e3, tol_grad=1e-12, tol_param=1e-12, tol_obj=1e-14,
                tol_rel_obj=1e2, iter=5000, history_size=20, require_converged=False)
    opts.update(kw)

    def run(start):
        # A line search that cannot decrease further is convergence at the
        # tolerances asked for, not an error.
        fit = model(name).optimize(
            data=data, inits=start, algorithm="lbfgs", jacobian=False, sig_figs=18,
            show_console=False, **opts,
        )
        return {k: float(v) for k, v in fit.optimized_params_dict.items()}

    out = run(inits)
    if polish:
        params = list(model(name).src_info()["parameters"].keys())
        for _ in range(4):
            nxt = run({k: out[k] for k in params})
            moved = max(abs(nxt[k] - out[k]) for k in params)
            if not np.isfinite(nxt["lp__"]) or nxt["lp__"] < out["lp__"]:
                break
            out = nxt
            if moved < 1e-9:
                break
    return out


def log_prob(name: str, data: dict, params: dict) -> tuple[float, np.ndarray]:
    """Log density and its gradient at ``params``, on the unconstrained scale
    with the Jacobian off; all six parameters here are unbounded so the two
    scales coincide."""
    df = model(name).log_prob(params=params, data=data, jacobian=False, sig_figs=18)
    lp = float(df["lp__"].iloc[0])
    grad = df.drop(columns="lp__").iloc[0].to_numpy(dtype=float)
    return lp, grad


def sample_four(name: str, data: dict, inits: dict, out_dir: Path, seed: int = 1,
                **kw) -> CmdStanMCMC:
    """Four chains as four separate jobs (Section 5), then combined."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    csvs = []
    for k in range(1, 5):
        fit = model(name).sample(
            data=data, inits=inits, chains=1, chain_ids=[k], seed=seed,
            sig_figs=10, output_dir=str(out_dir), show_progress=False, **kw,
        )
        csvs += list(fit.runset.csv_files)
    return from_csv(csvs)
