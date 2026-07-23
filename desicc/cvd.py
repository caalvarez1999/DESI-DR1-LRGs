"""Velocity-dispersion (CVD) correction for Lick/D4000 indices.

Thin wrapper around aux.corrections.C_VD: it looks up, for each galaxy's
own velocity dispersion, a precomputed correction (and its uncertainty) from
a grid built on MILES stellar templates broadened to a range of velocity
dispersions, and applies it index by index.
"""

import numpy as np
from astropy.table import Table

from aux.corrections import C_VD

from .config import CVD_C_FUNCS, CVD_DC_FUNCS, CVD_UVES_MAX, CVD_UVES_STEP, LICK_INDICES


def apply_cvd(pyl: Table, vd) -> dict:
    """Return {index: values, dindex: values}, corrected for velocity dispersion.

    `pyl` must have every column in LICK_INDICES (and its d-prefixed
    uncertainty); `vd` is the per-row velocity dispersion (km/s), aligned
    with `pyl`. Rows with vd or index values that are NaN come back as NaN,
    same as C_VD's own behavior.
    """
    uves = np.arange(0.0, CVD_UVES_MAX, CVD_UVES_STEP)
    C = np.load(CVD_C_FUNCS, allow_pickle=True)
    dC = np.load(CVD_DC_FUNCS, allow_pickle=True)

    nobj = len(pyl)
    nidx = len(LICK_INDICES)
    raw = np.zeros((nobj, nidx, 2), dtype=float)
    for j, name in enumerate(LICK_INDICES):
        raw[:, j, 0] = np.asarray(pyl[name], dtype=float)
        raw[:, j, 1] = np.asarray(pyl["d" + name], dtype=float)

    corrected = C_VD(raw, np.asarray(vd, dtype=float), uves, C, dC, LICK_INDICES)

    out = {}
    for j, name in enumerate(LICK_INDICES):
        out[name] = corrected[:, j, 0]
        out["d" + name] = corrected[:, j, 1]
    return out
