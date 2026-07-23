"""Helpers for reading PARENT tables and parsing pipeline inputs."""

import uuid as uuid_mod

import numpy as np
from astropy.table import Table

# Values that show up where a real identifier is missing -- either because
# a left join found no match, or because we pre-filled the column with a
# placeholder in stage 00.
BAD_STRING_PLACEHOLDERS = {"", "--", "0.0", "nan", "none", "null"}


def decode(x) -> str:
    if x is None:
        return ""
    if isinstance(x, bytes):
        return x.decode("utf-8", errors="ignore").strip()
    return str(x).strip()


def is_placeholder(s) -> bool:
    return decode(s).lower() in BAD_STRING_PLACEHOLDERS


def is_valid_uuid(s) -> bool:
    s = decode(s)
    if is_placeholder(s):
        return False
    try:
        uuid_mod.UUID(s)
        return True
    except ValueError:
        return False


def find_column(table: Table, name: str) -> str:
    """Case-insensitive column lookup -- FITS readers don't all agree on case."""
    lname = name.lower()
    for c in table.colnames:
        if c.lower() == lname:
            return c
    raise KeyError(f"Column '{name}' not found. Available: {table.colnames}")


def parse_vd_bin(vd_bin: str) -> tuple[int, int]:
    """Parse a 'vd bin' string like '160180' into (160, 180): 160 <= vd < 180."""
    s = vd_bin.strip()
    if len(s) != 6 or not s.isdigit():
        raise ValueError('Expected a 6-digit bin, e.g. "160180" for 160 <= vd < 180')
    return int(s[:3]), int(s[3:])


def normalize_specid(x) -> int:
    """specid sometimes arrives as a string with stray characters; keep only the digits."""
    digits = "".join(ch for ch in decode(x) if ch.isdigit())
    return int(digits) if digits else 0


def to_str_and_valid(col, fill: str = "") -> tuple[np.ndarray, np.ndarray]:
    """Convert a table column to (string array, boolean valid mask).

    Handles masked columns (e.g. rows a left join didn't match) and the
    placeholder strings above; `valid` is True only where the value is a
    real, usable result.
    """
    arr = np.array(col)
    if isinstance(arr, np.ma.MaskedArray):
        valid = ~arr.mask
        arr_str = arr.astype(str).filled(fill)
    else:
        valid = np.ones(len(arr), dtype=bool)
        arr_str = arr.astype(str)

    arr_str = np.char.strip(arr_str)
    valid &= ~np.array([is_placeholder(s) for s in arr_str])
    return arr_str, valid
