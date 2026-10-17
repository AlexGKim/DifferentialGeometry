"""Stage 4 of Section 5.2: the hierarchy on the training sample.

The same program as stage 3 on more objects, launched detached so that a run
of days outlives the shell that started it.

    python scripts/stage4_full.py --n-sn 100 --check-only   cost and gradient check
    python scripts/stage4_full.py --n-sn 100 --iter 300 --max-treedepth 10
    python scripts/stage4_full.py --iter 300 --max-treedepth 8
    python scripts/stage4_full.py --report --out runs/stage4-smoke

On capping and short runs.  A dense metric, which is what made stage 3
converge, does not scale here: the full sample carries 442 * 6 + 36 + 10
parameters, so a dense metric is a matrix of order 2700 squared, costs that
many flops in a matrix-vector product every leapfrog, and cannot be estimated
from fewer warmup draws than it has rows.  Stage 4 therefore runs on the
diagonal metric, under which thirty objects already spent most transitions at
maximum treedepth, and the cost of a fit is leapfrogs times the cost of a
gradient --- of order a second at the full sample.  A converged stage 4 is not
reachable that way, so the runs below are capped in treedepth and short, and
their purpose is to show the program works at full size and to measure the
geometry there.  Read R-hat and the effective sample sizes they print as
descriptions of a truncated sampler, not as convergence.
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
from cmdstanpy import from_csv

from dgsn.data.sample import analysis_sample, median_x1, training_split, x1_fifth
from dgsn.salt.globals import fit_globals
from dgsn.stan.hier_data import build, inits
from dgsn.stan.run import sample_four

sys.path.insert(0, str(Path(__file__).parent))
from stage3_hier import cost, gradient_check, report  # noqa: E402

K_PSI, K_Q = 12, 8          # the smallest orders, as stage 3; raise only once a fit is healthy


def stratified(names: list[str], n: int, sample) -> list[str]:
    """``n`` of ``names`` taken round-robin across the fifths of ``x1``, so a
    subsample spans the range along which the difficulty of this model varies
    (Section 5.1) rather than whatever the name order gives."""
    if n >= len(names):
        return names
    fifth = x1_fifth(sample.loc[names])
    buckets = [[nm for nm in names if fifth[nm] == q] for q in range(5)]
    out: list[str] = []
    i = 0
    while len(out) < n:
        b = buckets[i % 5]
        j = i // 5
        if j < len(b):
            out.append(b[j])
        i += 1
        if i > 5 * max(len(b) for b in buckets):
            break
    return out[:n]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-only", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--out", default=None, help="default runs/stage4 or runs/stage4-<n_sn>")
    ap.add_argument("--n-sn", type=int, default=None, help="objects, stratified in x1 (default: all 442)")
    ap.add_argument("--k-psi", type=int, default=K_PSI)
    ap.add_argument("--k-q", type=int, default=K_Q)
    ap.add_argument("--threads", type=int, default=3)
    ap.add_argument("--iter", type=int, default=300, help="warmup and sampling iterations each")
    ap.add_argument("--max-treedepth", type=int, default=8)
    ap.add_argument("--adapt-delta", type=float, default=None)
    ap.add_argument("--dense", action="store_true", help="dense_e; see the note at the top of this file")
    ap.add_argument("--wait", action="store_true", help="block until the chains finish, then report")
    ap.add_argument("--precondition-from", default=None, metavar="DIR",
                    help="rotate the coefficients by the Fisher information at the posterior "
                         "mean of the run in DIR; built at a start instead it makes things worse")
    args = ap.parse_args()

    sample = analysis_sample()
    train, _ = training_split()
    names = stratified(train, args.n_sn, sample) if args.n_sn else train
    out = Path(args.out or (f"runs/stage4-{len(names)}" if args.n_sn else "runs/stage4"))

    # grainsize: one slice per thread per chain, so reduce_sum splits the
    # sample rather than handing one worker everything
    grain = max(1, len(names) // max(1, args.threads))
    data = build(names, args.k_psi, args.k_q, grainsize=grain, sample=sample)
    print(f"{len(names)} objects, {data['N']} epochs, K = ({args.k_psi}, {args.k_q}), "
          f"grainsize {grain}, x1 range [{sample.loc[names].x1.min():+.2f}, {sample.loc[names].x1.max():+.2f}]")

    if args.report:
        csvs = sorted(out.glob("hier-?.csv"))
        print(f"{len(csvs)} chain files in {out}")
        report(from_csv([str(c) for c in csvs]), names, sample, data)
        return

    x1_ref = median_x1()
    G = fit_globals(args.k_psi, args.k_q, x1_ref)
    ini = inits(names, args.k_psi, args.k_q, x1_ref, sample=sample, globals_=G)

    R = None
    if args.precondition_from:
        from dgsn.stan.precondition import from_run, source_size
        # The source run knows nothing of objects it never saw, so the point
        # at which the sensitivities are taken uses its globals and this
        # fit's own start for the per-supernova parameters.
        per = {k: np.asarray(ini[k], dtype=float) for k in ("m_peak", "colr", "x_s")}
        per["x_p"] = ini["mu_p"] + np.asarray(ini["u_p"], dtype=float)
        per["x_psi"] = ini["mu_psi"] + np.asarray(ini["u_psi"], dtype=float)
        per["delta_t"] = ini["mu_delta"] + ini["s_delta"] * np.asarray(ini["u_del"], dtype=float)
        R = from_run(args.precondition_from, data, args.k_psi, args.k_q, per_sn=per)
        data = build(names, args.k_psi, args.k_q, grainsize=grain, sample=sample, R=R)
        ini = inits(names, args.k_psi, args.k_q, x1_ref, sample=sample, globals_=G, R=R)
        print(f"preconditioned from {args.precondition_from} "
              f"({source_size(args.precondition_from)} objects -> {len(names)}), "
              f"rotation condition {np.linalg.cond(R):.1e}")

    lp0, worst = gradient_check(data, ini)
    print(f"lp at start {lp0:.2f}; gradient vs finite differences, worst of 20 coordinates: {worst:.1e}")
    full, nodes = cost(data, ini)
    per_grad = full - 0.03                      # process launch, ~30 ms
    leap = 2 ** args.max_treedepth - 1
    # Threading buys far less than its thread count once four chains contend for
    # six performance cores: measured at 100 objects, a lone single-threaded
    # chain ran at 5.7 gradients/s and four concurrent three-thread chains at
    # 6.9 each, so the wall-clock factor over this single-call timing is ~1.2,
    # not threads * 0.8.  reduce_sum waits on its slowest slice.
    contend = 1.2 if args.threads > 1 else 1.0
    print(f"log_prob {full:.3f} s per call ({nodes:.3f} s at one epoch per object) -> ~{per_grad:.2f} s "
          f"per gradient single-threaded; at most {leap} leapfrogs per draw and {2 * args.iter} iterations, "
          f"that is up to {per_grad * leap * 2 * args.iter / 3600:.0f} h per chain single-threaded, "
          f"~{per_grad * leap * 2 * args.iter / 3600 / contend:.0f} h with four chains at {args.threads} threads each")
    if args.check_only:
        return

    t = time.time()
    res = sample_four("hier", data, ini, out, seed=1, threads_per_chain=args.threads,
                      iter_warmup=args.iter, iter_sampling=args.iter,
                      max_treedepth=args.max_treedepth,
                      metric="dense_e" if args.dense else None,
                      adapt_delta=args.adapt_delta, refresh=10, detach=not args.wait)
    if not args.wait:
        print(f"launched four chains, pids {res}; logs {out}/chain-*.log")
        print(f"progress: grep -h Iteration {out}/chain-*.log | tail -4")
        print(f"when done: python scripts/stage4_full.py --report --out {out}"
              + (f" --n-sn {args.n_sn}" if args.n_sn else "")
              + f" --k-psi {args.k_psi} --k-q {args.k_q}")
        return
    print(f"\nfour chains in {time.time() - t:.0f} s wall")
    report(res, names, sample, data)


if __name__ == "__main__":
    main()
