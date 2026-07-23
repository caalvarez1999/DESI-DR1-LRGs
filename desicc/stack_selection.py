"""Select which galaxies from PARENT4.fits (CVD-corrected) go into a given
vd/archaeology/cosmology stack, before stage 10 builds it.

The redshift limit per vd group comes from a precomputed lookup: for a given
scaling relation between velocity dispersion and stellar age (the
"archaeology") and a cosmology (age -> redshift), it's the highest redshift
at which a galaxy of that vd's minimum plausible age could still be observed
-- i.e. how far back this vd group lets you look. Those lookups live in
archaeology/limitesarchaeologyTMJ_<COSMOLOGY>_<ARCHAEOLOGY>.npy (one file per
cosmology/archaeology combination, each mapping every vd group to a zcut).

Beyond the redshift cut, a galaxy also has to have a finished, good-quality
Lick measurement (done=1, qual>=0.5), show no significant nebular emission in
[OII], Hb, [OIII]5007 (Ha is excused if unmeasured -- emission dilutes/biases
absorption-line indices), and a CaII H/K ratio under 1.2 (removes stars and
AGN-contaminated spectra that slipped through the galaxy-type cut in stage 00).
"""

import numpy as np
from astropy.table import Table

from .catalog import find_column, parse_vd_bin
from .config import ARCHAEOLOGY_DIR

EMISSION_INDICES = ("OII", "Hb", "OIII5007")
HA_INDEX = "Ha"
CA_H, CA_K = "CaIIH", "CaIIK"
CA_RATIO_MAX = 1.2
QUAL_MIN = 0.5


def load_zcut(cosmology: str, archaeology: str, vd_bin: str) -> float:
    path = ARCHAEOLOGY_DIR / f"limitesarchaeologyTMJ_{cosmology}_{archaeology}.npy"
    if not path.exists():
        raise FileNotFoundError(f"No archaeology/cosmology lookup at {path}")

    d = np.load(path, allow_pickle=True).item()
    groups = np.asarray(d["grupo"]).astype(str)
    zcuts = np.asarray(d["redshift"]).astype(float)

    match = np.where(groups == vd_bin)[0]
    if len(match) != 1:
        raise ValueError(f"vd bin '{vd_bin}' not found (or ambiguous) in {path}")
    return float(zcuts[match[0]])


def _passes_ew_cut(ew: np.ndarray, err: np.ndarray) -> np.ndarray:
    """EW < 0 means emission. Accept EW > -2, or -5 < EW <= -2 if it's not a
    significant detection (|EW|/|err| < 2)."""
    return (ew > -2.0) | ((ew > -5.0) & (np.abs(ew) / np.abs(err) < 2.0))


def emission_mask(pyl: Table) -> np.ndarray:
    n = len(pyl)
    ok = np.ones(n, dtype=bool)

    for name in EMISSION_INDICES:
        ew = np.asarray(pyl[name], dtype=float)
        dew = np.asarray(pyl["d" + name], dtype=float)
        valid = np.isfinite(ew) & np.isfinite(dew) & (dew != 0.0)
        good = np.zeros(n, dtype=bool)
        good[valid] = _passes_ew_cut(ew[valid], dew[valid])
        ok &= valid & good

    ew_ha = np.asarray(pyl[HA_INDEX], dtype=float)
    dew_ha = np.asarray(pyl["d" + HA_INDEX], dtype=float)
    valid_ha = np.isfinite(ew_ha) & np.isfinite(dew_ha) & (dew_ha != 0.0)
    good_ha = np.zeros(n, dtype=bool)
    good_ha[valid_ha] = _passes_ew_cut(ew_ha[valid_ha], dew_ha[valid_ha])

    return ok & (~valid_ha | good_ha)


def ca_ratio_mask(pyl: Table) -> np.ndarray:
    ca_h = np.asarray(pyl[CA_H], dtype=float)
    ca_k = np.asarray(pyl[CA_K], dtype=float)
    valid = np.isfinite(ca_h) & np.isfinite(ca_k) & (ca_k != 0.0)

    ratio = np.full(len(pyl), np.inf, dtype=float)
    ratio[valid] = ca_h[valid] / ca_k[valid]
    return valid & (ratio < CA_RATIO_MAX)


def select_for_stacking(meta: Table, pyl: Table, vd_bin: str, zcut: float) -> np.ndarray:
    """Boolean mask into (meta, pyl): rows eligible for this vd/archaeology/
    cosmology stack. meta and pyl must be row-aligned (e.g. PARENT4.fits)."""
    vd_min, vd_max = parse_vd_bin(vd_bin)
    vd_col = find_column(meta, "vd")
    z_col = find_column(meta, "z1")
    done_col = find_column(pyl, "done")
    qual_col = find_column(pyl, "qual")

    vd = np.asarray(meta[vd_col], dtype=float)
    z1 = np.asarray(meta[z_col], dtype=float)
    done = np.asarray(pyl[done_col], dtype=float)
    qual = np.asarray(pyl[qual_col], dtype=float)

    in_vd = (vd >= vd_min) & (vd < vd_max)
    in_z = np.isfinite(z1) & (z1 <= zcut)
    done_ok = np.isfinite(done) & (done.astype(int) == 1)
    qual_ok = np.isfinite(qual) & (qual >= QUAL_MIN)

    return in_vd & in_z & done_ok & qual_ok & emission_mask(pyl) & ca_ratio_mask(pyl)
