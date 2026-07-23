#!/usr/bin/env python3
"""Stage 00 - build the working PARENT table from the DESI fastspec-iron catalog.

The fastspec-iron release ships METADATA and FASTSPEC as two HDUs of the
same file, one row per target, in the same order. An earlier version of
this script treated them as two independent catalogs and joined them on
TARGETID -- which silently dropped a handful of galaxies whenever the join
wasn't perfect. Since the two HDUs are already row-aligned, the fix is to
not join at all: just apply the same boolean mask to both.

Selection applied:
  - LRG target, not also flagged ELG/QSO
  - SPECTYPE == GALAXY
  - Redshift agrees across METADATA (Z vs Z_RR) and FASTSPEC (Z)
  - Velocity dispersion and its uncertainty are both finite
  - VD_MIN < velocity dispersion < VD_MAX (desicc/config.py)

Output: METADATA (galaxy properties, SPARCL ids left blank) + an empty
PYLICK HDU (index columns pre-allocated as NaN, filled in stage 03).
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from desicc.config import LICK_INDICES, PARENT_MERGED, RAW_CATALOG, VD_MAX, VD_MIN


def build_mask(meta, fast) -> np.ndarray:
    desi_target = meta["DESI_TARGET"]
    is_lrg = (desi_target & (1 << 0)) != 0
    is_elg = (desi_target & (1 << 1)) != 0
    is_qso = (desi_target & (1 << 2)) != 0
    is_galaxy = meta["SPECTYPE"] == "GALAXY"

    z_consistent = (
        (meta["ZWARN"] == 0)
        & ~np.isnan(meta["Z"]) & ~np.isnan(meta["Z_RR"])
        & (np.abs(meta["Z"] - meta["Z_RR"]) / (1.0 + meta["Z_RR"]) < 0.001)
        & (np.abs(fast["Z"] - meta["Z_RR"]) / (1.0 + meta["Z_RR"]) < 0.001)
        & ~np.isnan(fast["Z"])
    )

    dvd = np.where(fast["VDISP_IVAR"] > 0, 1.0 / np.sqrt(fast["VDISP_IVAR"]), np.nan)
    has_vdisp = ~np.isnan(fast["VDISP"]) & ~np.isnan(dvd)
    in_vd_range = (fast["VDISP"] > VD_MIN) & (fast["VDISP"] < VD_MAX)

    return is_lrg & ~is_elg & ~is_qso & is_galaxy & z_consistent & has_vdisp & in_vd_range


def build_parent_hdulist(meta, fast, mask: np.ndarray) -> fits.HDUList:
    n = int(mask.sum())
    targetid = meta["TARGETID"][mask]

    vd_ivar = fast["VDISP_IVAR"][mask]
    dvd = np.where(vd_ivar > 0, 1.0 / np.sqrt(vd_ivar), np.nan)
    dn4000_ivar = fast["DN4000_IVAR"][mask]
    ddn4000 = np.where(dn4000_ivar > 0, 1.0 / np.sqrt(dn4000_ivar), np.nan)

    meta_cols = [
        fits.Column(name="targetid", format="K", array=targetid),
        fits.Column(name="specid", format="17A", array=np.full(n, "", dtype="U17")),
        fits.Column(name="sparclid", format="36A", array=np.full(n, "", dtype="U36")),
        fits.Column(name="DR", format="8A", array=np.full(n, "", dtype="U8")),
        fits.Column(name="survey", format="A4", array=fast["SURVEY"][mask]),
        fits.Column(name="program", format="A4", array=fast["PROGRAM"][mask]),
        fits.Column(name="healpix", format="J", array=fast["HEALPIX"][mask]),
        fits.Column(name="RA", format="D", array=meta["RA"][mask]),
        fits.Column(name="DEC", format="D", array=meta["DEC"][mask]),
        fits.Column(name="desitarget", format="K", array=meta["DESI_TARGET"][mask]),
        fits.Column(name="spectype", format="A10", array=meta["SPECTYPE"][mask]),
        fits.Column(name="subtype", format="A10", array=meta["SUBTYPE"][mask]),
        fits.Column(name="z1", format="D", array=fast["Z"][mask]),
        fits.Column(name="vd", format="E", array=fast["VDISP"][mask]),
        fits.Column(name="dvd", format="E", array=dvd),
        fits.Column(name="FSF_SFR", format="E", array=fast["SFR"][mask]),
        fits.Column(name="FSF_Dn4000", format="E", array=fast["DN4000"][mask]),
        fits.Column(name="FSF_dDn4000", format="E", array=ddn4000),
        fits.Column(name="snmedian", format="E", array=np.full(n, np.nan, dtype="f4")),
    ]
    meta_hdu = fits.BinTableHDU.from_columns(meta_cols, name="METADATA")

    pylick_cols = [
        fits.Column(name="targetid", format="K", array=targetid),
        fits.Column(name="done", format="J", array=np.zeros(n, dtype="i4")),
        fits.Column(name="qual", format="E", array=np.full(n, np.nan, dtype="f4")),
    ]
    for name in LICK_INDICES:
        pylick_cols.append(fits.Column(name=name, format="E", array=np.full(n, np.nan, dtype="f4")))
        pylick_cols.append(fits.Column(name=f"d{name}", format="E", array=np.full(n, np.nan, dtype="f4")))
    pylick_hdu = fits.BinTableHDU.from_columns(pylick_cols, name="PYLICK")

    return fits.HDUList([fits.PrimaryHDU(), meta_hdu, pylick_hdu])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, default=RAW_CATALOG, help="fastspec-iron FITS file")
    parser.add_argument("--output", type=Path, default=PARENT_MERGED)
    args = parser.parse_args()

    with fits.open(args.input, memmap=True) as hdul:
        meta = hdul["METADATA"].data
        fast = hdul["FASTSPEC"].data

        mask = build_mask(meta, fast)
        print(f"{mask.sum()} / {len(mask)} targets pass the selection")

        hdul_out = build_parent_hdulist(meta, fast, mask)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    hdul_out.writeto(args.output, overwrite=True)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
