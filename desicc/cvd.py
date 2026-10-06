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


def _load_cvd_grids():
    uves = np.arange(0.0, CVD_UVES_MAX, CVD_UVES_STEP)
    C = np.load(CVD_C_FUNCS, allow_pickle=True)
    dC = np.load(CVD_DC_FUNCS, allow_pickle=True)
    return uves, C, dC


def apply_cvd(pyl: Table, vd) -> dict:
    """Return {index: values, dindex: values}, corrected for velocity dispersion.

    `pyl` must have every column in LICK_INDICES (and its d-prefixed
    uncertainty); `vd` is the per-row velocity dispersion (km/s), aligned
    with `pyl`. Rows with vd or index values that are NaN come back as NaN,
    same as C_VD's own behavior.
    """
    uves, C, dC = _load_cvd_grids()

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


def apply_cvd_named(values: dict, vd) -> dict:
    """Correct an arbitrary subset of LICK_INDICES-named quantities for velocity
    dispersion, reusing each named index's own MILES correction curve.

    `values` is {name: (value_array, error_array)} for names that must all be
    in LICK_INDICES (e.g. a fastspecfit emission-line EW cross-matched onto a
    same-named Lick index's bandpass, see pipeline/05_fastspeclines.py). Every
    other slot in the underlying 33-index array C_VD expects is padded with
    NaN and dropped again on the way out, so this applies exactly the same
    per-index correction factor as apply_cvd would for that index, without
    requiring values for the other 32.

    Returns {name: values, dname: values} for exactly the names passed in.
    """
    uves, C, dC = _load_cvd_grids()

    vd = np.asarray(vd, dtype=float)
    nobj = len(vd)
    nidx = len(LICK_INDICES)
    raw = np.full((nobj, nidx, 2), np.nan, dtype=float)
    for name, (value, error) in values.items():
        j = LICK_INDICES.index(name)
        raw[:, j, 0] = np.asarray(value, dtype=float)
        raw[:, j, 1] = np.asarray(error, dtype=float)

    corrected = C_VD(raw, vd, uves, C, dC, LICK_INDICES)

    out = {}
    for name in values:
        j = LICK_INDICES.index(name)
        out[name] = corrected[:, j, 0]
        out["d" + name] = corrected[:, j, 1]
    return out
