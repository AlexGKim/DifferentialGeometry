"""Validate src/dgsn/stan/diag.py against CmdStan's own stansummary.

The docstring of that module claims its ESS_bulk and R_hat are CmdStan's.
This is the claim, run: stansummary over a finished run, column by column
against diag.py on the same CSVs.  stansummary prints two significant
figures by default, so the comparison is to the digits it gives.

    python scripts/check_diag.py runs/stage3-common-halfnorm
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
from cmdstanpy import cmdstan_path

from dgsn.stan.diag import ess_bulk, read_csvs, rhat


def main():
    run = Path(sys.argv[1] if len(sys.argv) > 1 else "runs/stage3-common-halfnorm")
    csvs = sorted(run.glob("*-[0-9].csv"))
    if not csvs:
        sys.exit(f"no chain CSVs under {run}")
    exe = Path(cmdstan_path()) / "bin" / "stansummary"
    out = subprocess.run([str(exe), "--sig_figs=6", *map(str, csvs)],
                         capture_output=True, text=True)
    if out.returncode:
        sys.exit(f"stansummary failed: {out.stderr.strip() or out.stdout.strip()}\n"
                 "a chain killed mid-row leaves a short line it will not read; "
                 "diag.py skips such a line, which is why it can read a live run")
    head, ref = None, {}
    for line in out.stdout.splitlines():
        p = line.split()
        if p[:1] == ["Mean"] or (head is None and "ESS_bulk" in p):
            head = p
            continue
        if head and len(p) == len(head) + 1:
            try:
                ref[p[0]] = dict(zip(head, [float(x) for x in p[1:]]))
            except ValueError:
                pass
    if not ref:
        sys.exit("could not parse stansummary output")
    scalars = [k for k, v in ref.items()
               if "[" not in k and np.isfinite(v.get("R_hat", np.nan))]
    draws, n = read_csvs(csvs, scalars)
    print(f"{run}: {len(draws)} scalar columns, {n} draws x {len(csvs)} chains\n")
    print(f"{'column':16s} {'ESS_bulk stansummary':>21s} {'diag.py':>10s} "
          f"{'R_hat stansummary':>19s} {'diag.py':>10s}")
    worst_e = worst_r = 0.0
    for k in sorted(draws):
        e_ref, r_ref = ref[k]["ESS_bulk"], ref[k]["R_hat"]
        e_got, r_got = ess_bulk(draws[k]), rhat(draws[k])
        if not (np.isfinite(e_ref) and np.isfinite(e_got)):
            continue
        de = abs(e_got - e_ref) / max(e_ref, 1e-12)
        dr = abs(r_got - r_ref)
        worst_e, worst_r = max(worst_e, de), max(worst_r, dr)
        print(f"{k:16s} {e_ref:21.2f} {e_got:10.2f} {r_ref:19.4f} {r_got:10.4f}")
    print(f"\nworst relative ESS_bulk difference {worst_e:.2%}; "
          f"worst absolute R_hat difference {worst_r:.4f}")
    ok = worst_e < 0.02 and worst_r < 0.01
    print("PASS" if ok else "FAIL -- diag.py is not reproducing stansummary")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
