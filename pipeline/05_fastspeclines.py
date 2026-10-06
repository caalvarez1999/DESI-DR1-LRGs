#!/usr/bin/env python3
"""Stage 05 - cross-match PARENT4 against fastspecfit's FASTSPEC HDU and pull
the emission-line EW values that double up with one of pylick's Lick indices,
CVD-corrected the same way PARENT4's own Lick indices are.

Of the 33 Lick/D4000 indices pylick measures (desicc/config.py:LICK_INDICES),
only 9 have an actual emission line fastspecfit fits an EW for -- the rest
(CaIIK, CaIIH, D4000, Dn4000, and every Fe/Mg/TiO/CN/Ca molecular-or-metal
index) are stellar-continuum absorption features with no fastspecfit
counterpart. See desicc/config.py:FASTSPEC_EW_LINES/FASTSPEC_EW_SOURCE_COLUMNS
for the exact mapping (OII sums the [OII] 3726/3729 doublet; HdA/HdF, HgA/HgF
and Hb0/Hb are two bandpass definitions of the same fastspecfit line each).

Cross-match is by targetid (unique in both PARENT4 and the fastspec-iron
catalog, so a plain lookup, no ambiguity). Each matched EW value is then
corrected with the same per-index MILES velocity-dispersion correction curve
PARENT4 used for that Lick index (desicc.cvd.apply_cvd_named) -- i.e. the
HALPHA_EW value picked up for "Ha" is multiplied by the same C(vd) factor
that corrected PARENT4's own "Ha" Lick index for that galaxy.

Output: PARENT4.fits's METADATA and PYLICK HDUs, unchanged, plus a new
'FASTSPEC9LINESCVD' HDU (targetid + 9 CVD-corrected EW values and their
uncertainties).

Example:
    python3 pipeline/05_fastspeclines.py
"""

import argparse
import sys
from pathlib import Path

from astropy.io import fits
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from desicc.catalog import find_column
from desicc.config import (
    FASTSPEC_EW_HDU_NAME,
    FASTSPEC_EW_LINES,
    PARENT_FASTSPEC_CVD,
    PARENT_FINAL_CVD,
    RAW_CATALOG,
)
from desicc.cvd import apply_cvd_named
from desicc.fastspec_lines import crossmatch_index, extract_fastspec_ew, raw_fastspec_columns, read_fastspec_columns


def build_fastspec_hdu(targetid, corrected: dict) -> fits.BinTableHDU:
    cols = [fits.Column(name="targetid", format="K", array=targetid)]
    for name in FASTSPEC_EW_LINES:
        cols.append(fits.Column(name=name, format="E", array=corrected[name].astype("f4")))
        cols.append(fits.Column(name=f"d{name}", format="E", array=corrected[f"d{name}"].astype("f4")))
    return fits.BinTableHDU.from_columns(cols, name=FASTSPEC_EW_HDU_NAME)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--parent", type=Path, default=PARENT_FINAL_CVD, help="PARENT4.fits (CVD-corrected)")
    parser.add_argument("--fastspec", type=Path, default=RAW_CATALOG, help="fastspec-iron FITS file")
    parser.add_argument("--output", type=Path, default=PARENT_FASTSPEC_CVD)
    args = parser.parse_args()

    with fits.open(args.parent, memmap=False) as hdul:
        meta = Table(hdul["METADATA"].data)
        pyl = Table(hdul["PYLICK"].data)
        meta_hdu = hdul["METADATA"].copy()
        pyl_hdu = hdul["PYLICK"].copy()

    colnames = raw_fastspec_columns(FASTSPEC_EW_LINES)
    columns = read_fastspec_columns(args.fastspec, colnames)
    idx = crossmatch_index(meta["targetid"], columns["TARGETID"])
    raw = extract_fastspec_ew(columns, idx, FASTSPEC_EW_LINES)

    vd_col = find_column(meta, "vd")
    corrected = apply_cvd_named(raw, meta[vd_col])

    fastspec_hdu = build_fastspec_hdu(meta["targetid"], corrected)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out = fits.HDUList([fits.PrimaryHDU(), meta_hdu, pyl_hdu, fastspec_hdu])
    out.writeto(args.output, overwrite=True)
    print(f"Cross-matched {len(meta)} targets against {args.fastspec}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
