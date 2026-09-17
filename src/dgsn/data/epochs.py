"""The epoch set the likelihood sees (Section 4), per supernova.

Distinct from the counting in ``sample.py``: no signal-to-noise cut, every
``g`` and ``r`` epoch with the quality mask clear, inside the rest-frame
selection window placed on the archive's ``t0``, and not among the epochs
the baseline was estimated from.  The ``flux`` column already carries the
baseline correction --- every light-curve header reads
``baseline_corr (applied)`` --- so it is taken as published.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from dgsn.data.sample import QUALITY_MASK, light_curve
from dgsn.stan import data as D

BAND_INDEX = {"ztfg": 1, "ztfr": 2}


def likelihood_epochs(name: str, t0: float, z: float, root: Path | None = None) -> dict:
    """``dt = mjd - t0``, flux, its error and the band index (1 = g, 2 = r),
    plus ``m_bright``, the magnitude of the brightest good ``g`` epoch at
    signal-to-noise above five: the SED-free anchor the ``m_peak`` prior is
    centred on."""
    lc = light_curve(name, root)
    phase = (lc.mjd - t0) / (1 + z)
    ok = (
        lc["filter"].isin(BAND_INDEX)
        & ((lc.flag.astype(int) & QUALITY_MASK) == 0)
        & (lc.in_baseline.astype(int) == 0)
        & (phase >= D.DATA_WINDOW[0])
        & (phase <= D.DATA_WINDOW[1])
    )
    sel = lc.loc[ok]
    g_det = sel[(sel["filter"] == "ztfg") & (sel.flux / sel.flux_err > 5.0)]
    m_bright = float(D.ZP - 2.5 * np.log10(g_det.flux.max())) if len(g_det) else np.nan
    return dict(
        dt=(sel.mjd - t0).to_numpy(dtype=float),
        f=sel.flux.to_numpy(dtype=float),
        sig_f=sel.flux_err.to_numpy(dtype=float),
        band=sel["filter"].map(BAND_INDEX).to_numpy(dtype=int),
        z=float(z),
        m_bright=m_bright,
        n_g=int((sel["filter"] == "ztfg").sum()),
        n_r=int((sel["filter"] == "ztfr").sum()),
        n_g_prepeak=int(((sel["filter"] == "ztfg") & (sel.mjd < t0)).sum()),
    )
