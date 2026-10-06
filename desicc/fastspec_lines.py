"""Cross-match PARENT tables against fastspecfit's FASTSPEC HDU by targetid,
and pull out emission-line EW values for the lines that double up with one
of pylick's Lick indices. See pipeline/05_fastspeclines.py.

The fastspec-iron FASTSPEC HDU is a huge, wide (~900-column), row-major FITS
table (tens of GB) -- fancy-indexing the whole astropy FITS_rec (`data[idx]`)
forces it to materialize every column for the selected rows and reliably
runs a normal workstation out of memory. read_fastspec_ew_columns reads the
handful of columns this stage actually needs, one bounded row-chunk at a
time, so peak memory stays proportional to the columns we keep, not the
whole table.
"""

import numpy as np
from astropy.io import fits

from .config import FASTSPEC_EW_SOURCE_COLUMNS

DEFAULT_CHUNKSIZE = 200_000


def raw_fastspec_columns(names=None) -> list:
    """The underlying FASTSPEC column names (plus their _IVAR companions)
    needed to build every name in `names` (default: every FASTSPEC_EW_LINES
    name) via FASTSPEC_EW_SOURCE_COLUMNS, TARGETID included."""
    if names is None:
        names = FASTSPEC_EW_SOURCE_COLUMNS.keys()
    cols = {"TARGETID"}
    for name in names:
        for col in FASTSPEC_EW_SOURCE_COLUMNS[name]:
            cols.add(col)
            cols.add(f"{col}_IVAR")
    return sorted(cols)


def read_fastspec_columns(path, colnames, chunksize=DEFAULT_CHUNKSIZE) -> dict:
    """Read exactly `colnames` off the FASTSPEC HDU at `path`, chunked by row
    range so no more than `chunksize` rows' worth of the table is ever
    resident at once. Returns {colname: full-length 1-D array}.
    """
    with fits.open(path, memmap=True) as hdul:
        rec = hdul["FASTSPEC"].data
        nrows = len(rec)

        out = {}
        for name in colnames:
            dtype = "i8" if name == "TARGETID" else "f8"
            out[name] = np.empty(nrows, dtype=dtype)

        for start in range(0, nrows, chunksize):
            end = min(start + chunksize, nrows)
            chunk = rec[start:end]
            for name in colnames:
                out[name][start:end] = chunk[name]

    return out


def crossmatch_index(target_ids, source_ids) -> np.ndarray:
    """Return, for each id in `target_ids`, its row index in `source_ids`.

    Both are assumed unique (true of targetid in both PARENT tables and the
    fastspec-iron catalog). Raises if any `target_ids` value has no match.
    """
    target_ids = np.asarray(target_ids)
    source_ids = np.asarray(source_ids)

    order = np.argsort(source_ids)
    sorted_source = source_ids[order]
    pos = np.searchsorted(sorted_source, target_ids)
    pos = np.clip(pos, 0, len(sorted_source) - 1)
    idx = order[pos]

    missing = source_ids[idx] != target_ids
    if missing.any():
        raise ValueError(
            f"{missing.sum()} targetid value(s) have no match in the fastspec catalog"
        )
    return idx


def extract_fastspec_ew(columns: dict, idx, names) -> dict:
    """Return {name: (value, error)} for each of `names` (must be keys of
    FASTSPEC_EW_SOURCE_COLUMNS), reindexing `columns` (as read by
    read_fastspec_columns) with `idx` to match the target row order.

    A name backed by more than one fastspecfit column (the [OII] 3726/3729
    doublet) is summed; its error is the quadrature sum of the components'
    1/sqrt(ivar), NaN if any component's ivar isn't positive -- same
    ivar-invalid convention used throughout desicc (e.g. stage 00's vd/Dn4000
    errors).
    """
    idx = np.asarray(idx)
    out = {}
    for name in names:
        cols = FASTSPEC_EW_SOURCE_COLUMNS[name]

        value = np.zeros(len(idx), dtype=float)
        variance = np.zeros(len(idx), dtype=float)
        valid = np.ones(len(idx), dtype=bool)
        for col in cols:
            value += columns[col][idx]
            ivar = columns[f"{col}_IVAR"][idx]
            ok = ivar > 0
            variance += np.where(ok, 1.0 / ivar, 0.0)
            valid &= ok

        error = np.where(valid, np.sqrt(variance), np.nan)
        out[name] = (value, error)
    return out
