"""Selection and epoch counting for the ZTF SN Ia DR2 sample.

Produces the numbers quoted in the Sample section of
``docs/tex/model-definition.tex``.  Run as a module to regenerate them::

    python -m dgsn.data.sample

The archive root is taken from ``$ZTFCOSMODIR`` if set, otherwise from the
local unpacked copy under ``data/raw/``.

A note on the quality flags, because it has already caused one wrong number.
``lccoverage_flag`` and ``fitquality_flag`` are floats and are NaN for 18
objects.  Testing them with ``!= 0`` admits those NaNs, since ``NaN != 0`` is
True; that is where the previously quoted count of 669 came from.  Require
equality with 1.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Epochs are counted in this rest-frame phase window, at this signal-to-noise.
PHASE_WINDOW = (-15.0, 40.0)
MIN_SNR = 5.0
MAX_REDSHIFT = 0.05

# The ztfcosmo default: exclude bits 1, 2, 4, 8, 16 and no others.  Bits 32
# and above encode seeing, field, moon and airmass and are informational.
QUALITY_MASK = 31

BANDS = ("ztfg", "ztfr", "ztfi")

# DR2's own cosmology classification.  It contains norm, 91t and 99aa and no
# peculiar subtype; snia-pec holds 91bg, Iax, 03fg, 02es, 00cx and 18byg.  Note
# that ``sn_type != "snia-pec"`` is NOT equivalent: one 91bg/86G object is typed
# ``snia``, with the peculiarity recorded only in ``sub_type``.
SN_TYPE = "snia-cosmo"


def archive_root() -> Path:
    """Locate the DR2 archive."""
    env = os.environ.get("ZTFCOSMODIR")
    if env:
        return Path(env)
    return Path(__file__).resolve().parents[3] / "data/raw/ztfsniadr2_lite"


def metadata(root: Path | None = None) -> pd.DataFrame:
    """The full ``snia_data`` table, one row per catalogued object."""
    root = root or archive_root()
    return pd.read_csv(root / "tables/snia_data.csv")


def select(meta: pd.DataFrame) -> pd.DataFrame:
    """Normal SNe Ia passing the redshift cut with both quality flags set."""
    keep = (
        (meta.redshift < MAX_REDSHIFT)
        & (meta.lccoverage_flag == 1)
        & (meta.fitquality_flag == 1)
        & (meta.sn_type == SN_TYPE)
    )
    return meta.loc[keep]


def light_curve(name: str, root: Path | None = None) -> pd.DataFrame:
    """One light curve.  Whitespace delimited under a ``#`` comment header."""
    root = root or archive_root()
    return pd.read_csv(
        root / f"lightcurves/{name}_lc.csv", sep=r"\s+", comment="#"
    )


def good_epochs(name: str, t0: float, z: float, root: Path | None = None) -> dict[str, int]:
    """Count good epochs per band inside the phase window.

    ``filter`` is indexed rather than accessed as an attribute: ``lc.filter``
    resolves to the DataFrame method and silently compares unequal to every
    band name.
    """
    lc = light_curve(name, root)
    phase = (lc.mjd - t0) / (1 + z)
    ok = (
        ((lc.flag.astype(int) & QUALITY_MASK) == 0)
        & (lc.flux / lc.flux_err > MIN_SNR)
        & (phase >= PHASE_WINDOW[0])
        & (phase <= PHASE_WINDOW[1])
    )
    return {b: int((ok & (lc["filter"] == b)).sum()) for b in BANDS}


def epoch_counts(sub: pd.DataFrame, root: Path | None = None) -> pd.DataFrame:
    """Per-band good-epoch counts for every selected object with photometry."""
    root = root or archive_root()
    present = {p.name[:-7] for p in (root / "lightcurves").glob("*_lc.csv")}
    rows: list[dict[str, int]] = []
    names: list[str] = []
    for _, r in sub.iterrows():
        if r.ztfname not in present or np.isnan(r.t0):
            continue
        rows.append(good_epochs(r.ztfname, r.t0, r.redshift, root))
        names.append(r.ztfname)
    return pd.DataFrame(rows, index=pd.Index(names))


def analysis_sample(min_epochs: int = 5, root: Path | None = None) -> pd.DataFrame:
    """The objects the analysis uses: selected, with photometry, and at least
    ``min_epochs`` good epochs in each of ``g`` and ``r``.  Metadata rows,
    with the per-band counts joined on."""
    root = root or archive_root()
    sub = select(metadata(root)).set_index("ztfname")
    n = epoch_counts(sub.reset_index(), root)
    keep = n.index[(n.ztfg >= min_epochs) & (n.ztfr >= min_epochs)]
    return sub.loc[keep].join(n.loc[keep])


def median_x1(root: Path | None = None) -> float:
    """The sample-median SALT2 ``x1``: where the reference track of Section 5
    is taken."""
    return float(analysis_sample(root=root).x1.median())


# The training/validation split of Section 5.1: 22 per fifth of x1 withheld,
# 110 in all, from the 552.  One seed, so that stages 3, 4 and 5 see the same
# split; the hash printed by ``main`` is what a log is checked against.
SPLIT_SEED = 1
N_VALID_PER_FIFTH = 22


def training_split(seed: int = SPLIT_SEED, root: Path | None = None) -> tuple[list[str], list[str]]:
    """``(train, valid)`` names.  The analysis sample ordered by x1 and cut
    into fifths; ``N_VALID_PER_FIFTH`` drawn from each without replacement."""
    s = analysis_sample(root=root).sort_values("x1")
    fifths = np.array_split(s.index.to_numpy(), 5)
    rng = np.random.default_rng(seed)
    valid: list[str] = []
    for names in fifths:
        valid += sorted(rng.choice(names, N_VALID_PER_FIFTH, replace=False).tolist())
    train = sorted(set(s.index) - set(valid))
    return train, sorted(valid)


def split_hash(names: list[str]) -> str:
    """A short hash of a sorted name list, for the log."""
    import hashlib

    return hashlib.sha1("\n".join(sorted(names)).encode()).hexdigest()[:10]


def x1_fifth(sample: pd.DataFrame) -> pd.Series:
    """Which fifth of the analysis sample's x1 distribution each object is in, 0..4."""
    order = analysis_sample().x1.sort_values()
    edges = np.quantile(order, [0.2, 0.4, 0.6, 0.8])
    return pd.Series(np.searchsorted(edges, sample.x1.to_numpy()), index=sample.index)


def stage3_sample(n_per_fifth: int = 6, seed: int = SPLIT_SEED, min_g_prepeak: int = 3,
                  root: Path | None = None) -> list[str]:
    """``n_per_fifth`` training objects per fifth of x1, preferring those with
    pre-peak g coverage, as the stage-2 pick does; deterministic."""
    from dgsn.data.epochs import likelihood_epochs

    train, _ = training_split(seed, root)
    s = analysis_sample(root=root).loc[train]
    fifth = x1_fifth(s)
    rng = np.random.default_rng(seed)
    out: list[str] = []
    for q in range(5):
        names = s.index[fifth == q].to_numpy()
        rng.shuffle(names)
        good = [n for n in names
                if likelihood_epochs(n, s.loc[n].t0, s.loc[n].redshift, root)["n_g_prepeak"] >= min_g_prepeak]
        pick = (good + [n for n in names if n not in good])[:n_per_fifth]
        out += pick
    return out


def write_tex(counts: dict[str, int], path: Path) -> None:
    """Emit the counts as LaTeX macros for ``\\input`` by the note.

    The note quotes these numbers in running text; writing them from here is
    what keeps the two from drifting when the selection changes.
    """
    lines = ["% Generated by dgsn.data.sample --tex.  Do not edit.\n"]
    lines += [f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in counts.items()]
    path.write_text("".join(lines))


def main() -> None:
    root = archive_root()
    meta = metadata(root)
    sub = select(meta)
    present = {p.name[:-7] for p in (root / "lightcurves").glob("*_lc.csv")}
    missing = set(sub.ztfname) - present
    nan_t0 = set(sub[sub.t0.isna()].ztfname)

    print(f"archive                        {root}")
    print(f"catalogued objects             {len(meta)}")
    print(f"z < {MAX_REDSHIFT}, both flags set     {len(sub)}")
    print(f"  without a light-curve file   {len(missing)}")
    print(f"  with NaN SALT2 t0            {len(nan_t0)}")
    print(f"  those two sets identical     {missing == nan_t0}")

    n = epoch_counts(sub, root)
    print(f"objects with photometry        {len(n)}")
    print()
    print(f"{'min epochs/band':>16}  {'g':>5} {'r':>5} {'i':>5} {'g&r':>5}")
    for m in (3, 5, 10):
        print(
            f"{'>= ' + str(m):>16}  {int((n.ztfg >= m).sum()):5d}"
            f" {int((n.ztfr >= m).sum()):5d} {int((n.ztfi >= m).sum()):5d}"
            f" {int(((n.ztfg >= m) & (n.ztfr >= m)).sum()):5d}"
        )

    mw = sub.mwebv.dropna()
    print()
    print(f"mwebv median {mw.median():.3f}, 90th pct {mw.quantile(0.9):.3f}, max {mw.max():.3f}")
    for t in (0.1, 0.2, 0.3):
        print(f"  fraction above {t:.1f}          {(mw > t).mean():.3f}  ({int((mw > t).sum())})")

    train, valid = training_split(root=root)
    print()
    print(f"split (seed {SPLIT_SEED}): {len(train)} training, {len(valid)} validation; "
          f"hash train {split_hash(train)}, valid {split_hash(valid)}")

    if "--tex" in sys.argv:

        def both(m: int) -> int:
            return int(((n.ztfg >= m) & (n.ztfr >= m)).sum())

        out = Path(__file__).resolve().parents[3] / "docs/tex/sample-counts.tex"
        write_tex(
            {
                "Ncatalogued": len(meta),
                "Nselected": len(sub),
                "Nnolc": len(missing),
                "Nphot": len(n),
                "Nanalysis": both(5),
                "Nthree": both(3),
                "Nten": both(10),
            },
            out,
        )
        print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
