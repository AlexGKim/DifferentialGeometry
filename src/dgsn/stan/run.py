"""Compile once, then drive the programs under ``stan/``.

Four entry points, one per way the note uses Stan: ``fixed_param`` for the
forward map at fixed parameters, ``optimize`` for the stage-1 fits,
``log_prob`` for gradients, ``sample_four`` for the four chains of Section 5,
run as four separate jobs.
"""

from __future__ import annotations

from pathlib import Path

import logging
import sys

import numpy as np
from cmdstanpy import CmdStanMCMC, CmdStanModel, from_csv

# A restarted L-BFGS that cannot move at its first iteration is at the
# optimum, and cmdstanpy logs that as an error; real failures still raise.
logging.getLogger("cmdstanpy").setLevel(logging.CRITICAL)

STAN_DIR = Path(__file__).resolve().parents[3] / "stan"
_models: dict[str, CmdStanModel] = {}


THREADED = {"hier", "hier_nomode"}   # programs compiled with STAN_THREADS, for reduce_sum


def _cpp_options(name: str) -> dict | None:
    """Make variables for the build.

    On macOS the compiler is the system one rather than the environment's.
    The conda toolchain ships a linker that cannot read the libSystem stub in
    recent SDKs --- after an OS upgrade to Darwin 27 it failed to resolve
    ``printf`` in a three-line program, so no Stan model would build, while
    ``/usr/bin/clang++`` built them unchanged.  Models already compiled keep
    working, which is why this surfaces only on the first edit after an
    upgrade."""
    opts: dict = {}
    if name in THREADED:
        opts["STAN_THREADS"] = True
    if sys.platform == "darwin" and Path("/usr/bin/clang++").exists():
        opts["CXX"] = "/usr/bin/clang++"
    return opts or None


def model(name: str) -> CmdStanModel:
    """The compiled program ``stan/<name>.stan``; compiled on first use."""
    if name not in _models:
        _models[name] = CmdStanModel(
            stan_file=STAN_DIR / f"{name}.stan",
            stanc_options={"include-paths": [str(STAN_DIR)]},
            cpp_options=_cpp_options(name),
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


def optimize(name: str, data: dict, inits: dict, polish: bool = True,
             jacobian: bool = False, **kw) -> dict[str, float]:
    """Mode of the log density; the fitter of the stage-1 tests and floors.

    L-BFGS, restarted from its own result while the restart still moves it:
    a restart discards the Hessian approximation, which is what lets it
    resolve a direction the light curve is nearly flat along (x_psi).
    Newton is not used --- its step along that direction is enormous and it
    diverges to NaN.  Written back at 18 significant figures: CmdStan's
    default of 6 quantises the mode.

    ``jacobian`` selects which density is maximized: ``False`` gives the mode
    on the constrained scale, which is what a point estimate means, and
    ``True`` the mode of the density the sampler moves on, which is the point
    whose curvature can precondition it."""
    opts = dict(tol_rel_grad=1e3, tol_grad=1e-12, tol_param=1e-12, tol_obj=1e-14,
                tol_rel_obj=1e2, iter=5000, history_size=20, require_converged=False)
    opts.update(kw)

    def run(start):
        # A line search that cannot decrease further is convergence at the
        # tolerances asked for, not an error.
        fit = model(name).optimize(
            data=data, inits=start, algorithm="lbfgs", jacobian=jacobian, sig_figs=18,
            show_console=False, **opts,
        )
        return {k: float(v) for k, v in fit.optimized_params_dict.items()}

    out = run(inits)
    if polish:
        params = list(model(name).src_info()["parameters"].keys())

        def as_inits(flat: dict) -> dict:
            """CmdStan reports a vector parameter flat, as ``x[1]``, ``x[2]``,
            ...; a restart needs it whole again."""
            ini = {}
            for p in params:
                if p in flat:
                    ini[p] = flat[p]
                else:
                    keys = sorted((k for k in flat if k.startswith(p + "[")),
                                  key=lambda k: int(k[len(p) + 1:-1]))
                    ini[p] = np.array([flat[k] for k in keys])
            return ini

        for _ in range(4):
            nxt = run(as_inits(out))
            moved = max(float(np.max(np.abs(np.asarray(nxt[k]) - np.asarray(out[k]))))
                        for k in out if k != "lp__")
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


def sample_four(name: str, data: dict, inits, out_dir: Path, seed: int = 1,
                threads_per_chain: int = 1, parallel: bool = True, iter_warmup: int | None = None,
                iter_sampling: int | None = None, max_treedepth: int | None = None,
                metric: str | None = None, metric_file: Path | None = None,
                adapt_delta: float | None = None, refresh: int = 100,
                detach: bool = False, **kw) -> CmdStanMCMC | list[int]:
    """Four chains as four separate jobs (Section 5), then combined.

    Each chain is its own CmdStan process with its own output file
    ``<out_dir>/<name>-<k>.csv`` and log ``chain-<k>.log``; four cmdstanpy
    ``sample`` calls from threads share a timestamped filename and interleave.
    ``inits`` may be one dict or a list of four.  ``metric`` is CmdStan's
    ``diag_e`` (its default) or ``dense_e``; ``metric_file`` a JSON
    ``inv_metric`` over the unconstrained parameters to start from.

    ``detach`` puts each chain in its own session so that a run of days
    outlives the process that launched it, and returns the four pids rather
    than a fit; read such a run with ``from_csv`` when it finishes."""
    import subprocess

    from cmdstanpy import write_stan_json

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    inits4 = inits if isinstance(inits, list) else [inits] * 4
    write_stan_json(str(out_dir / "data.json"), data)
    exe = model(name).exe_file
    procs = []
    for k in range(1, 5):
        write_stan_json(str(out_dir / f"inits-{k}.json"), inits4[k - 1])
        cmd = [str(exe), "sample"]
        if iter_warmup is not None:
            cmd.append(f"num_warmup={iter_warmup}")
        if iter_sampling is not None:
            cmd.append(f"num_samples={iter_sampling}")
        if adapt_delta is not None:
            cmd += ["adapt", f"delta={adapt_delta}"]
        if max_treedepth is not None or metric is not None:
            cmd.append("algorithm=hmc")
            if metric is not None:
                cmd.append(f"metric={metric}")
            if metric_file is not None:
                cmd.append(f"metric_file={metric_file}")
            if max_treedepth is not None:
                cmd += ["engine=nuts", f"max_depth={max_treedepth}"]
        cmd += ["data", f"file={out_dir / 'data.json'}", f"init={out_dir / f'inits-{k}.json'}",
                "random", f"seed={seed}", f"id={k}",
                "output", f"file={out_dir / f'{name}-{k}.csv'}", "sig_figs=10", f"refresh={refresh}"]
        if name in THREADED:
            cmd.append(f"num_threads={threads_per_chain}")
        log = open(out_dir / f"chain-{k}.log", "w")
        p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT,
                             start_new_session=detach)
        procs.append(p)
        if not parallel and not detach:
            p.wait()
    if detach:
        return [p.pid for p in procs]
    for p in procs:
        p.wait()
    failed = [k for k, p in zip(range(1, 5), procs) if p.returncode != 0]
    if failed:
        raise RuntimeError(f"chains {failed} failed; see {out_dir}/chain-*.log")
    return from_csv([str(out_dir / f"{name}-{k}.csv") for k in range(1, 5)])


def seconds_per_gradient(fit: CmdStanMCMC) -> float:
    """Total sampling time over total leapfrogs, from the CSVs' timing lines
    and ``n_leapfrog__``; an upper bound, since warmup leapfrogs are not
    counted against warmup time separately."""
    total = 0.0
    for c in fit.runset.csv_files:
        for line in open(c):
            if "seconds (Total)" in line:
                total += float(line.split()[1])
    return total / float(fit.draws_pd()["n_leapfrog__"].sum())
