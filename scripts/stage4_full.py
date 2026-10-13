"""Stage 4 of Section 5.2: the full training sample at one pair of orders.
The same program as stage 3 on the 442 training objects, run detached.

    python scripts/stage4_full.py --check-only     projected wall time from stage 3
    python scripts/stage4_full.py                  launch four chains, three threads each
    python scripts/stage4_full.py --report         diagnostics of whatever CSVs exist
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from cmdstanpy import from_csv

from dgsn.data.sample import analysis_sample, median_x1, training_split
from dgsn.salt.globals import fit_globals
from dgsn.stan.hier_data import build, inits
from dgsn.stan.run import log_prob, model

sys.path.insert(0, str(Path(__file__).parent))
from stage3_hier import cost, report  # noqa: E402

K_PSI, K_Q = 24, 16
THREADS = 3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-only", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--out", default="runs/stage4")
    ap.add_argument("--leapfrogs", type=float, default=None,
                    help="leapfrogs per iteration to project with (default: read from runs/stage3)")
    ap.add_argument("--iter", type=int, default=1000)
    args = ap.parse_args()

    sample = analysis_sample()
    names, _ = training_split()
    x1_ref = median_x1()
    G = fit_globals(K_PSI, K_Q, x1_ref)
    data = build(names, K_PSI, K_Q, grainsize=max(1, len(names) // (4 * THREADS)), sample=sample)
    ini = inits(names, K_PSI, K_Q, x1_ref, sample=sample, globals_=G)
    print(f"{len(names)} training objects, {data['N']} epochs, K = ({K_PSI}, {K_Q}), grainsize {data['grainsize']}")

    out = Path(args.out)
    if args.report:
        csvs = sorted(out.glob("hier-?.csv"))
        print(f"{len(csvs)} chain files")
        fit = from_csv([str(c) for c in csvs])
        report(fit, names, sample, data)
        return

    full, nodes = cost(data, ini)
    per_grad = full - 0.03            # process launch, measured ~30 ms at stage 3
    leap = args.leapfrogs
    if leap is None:
        s3 = sorted(Path("runs/stage3").glob("hier-?.csv"))
        if s3:
            d = from_csv([str(c) for c in s3]).draws_pd()
            leap = float(d["n_leapfrog__"].mean())
        else:
            leap = 31.0
    hours = per_grad * leap * 2 * args.iter / 3600
    print(f"log_prob {full:.2f} s per call single-threaded -> ~{per_grad:.2f} s per gradient; "
          f"{leap:.0f} leapfrogs per iteration -> {hours:.0f} h per chain single-threaded, "
          f"~{hours / THREADS * 1.3:.0f} h at {THREADS} threads (reduce_sum efficiency ~75%)")
    if args.check_only:
        return

    out.mkdir(parents=True, exist_ok=True)
    import json
    from cmdstanpy import write_stan_json
    write_stan_json(str(out / "data.json"), data)
    write_stan_json(str(out / "inits.json"), ini)
    exe = model("hier").exe_file
    procs = []
    for k in range(1, 5):
        cmd = [str(exe), "sample", f"num_warmup={args.iter}", f"num_samples={args.iter}",
               "data", f"file={out / 'data.json'}", f"init={out / 'inits.json'}",
               "random", "seed=1", f"id={k}", "output", f"file={out / f'hier-{k}.csv'}",
               "sig_figs=10", "refresh=10", f"num_threads={THREADS}"]
        log = open(out / f"chain-{k}.log", "w")
        procs.append(subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, start_new_session=True))
        print("launched chain", k, "pid", procs[-1].pid)
    print(f"logs under {out}/chain-*.log; python scripts/stage4_full.py --report while it runs")


if __name__ == "__main__":
    main()
