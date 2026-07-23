#!/usr/bin/env python3
"""Stage 04 - apply the velocity-dispersion (CVD) correction to PARENT3.

Reads the current PARENT3.fits (written and updated by stage 03 across vd
bins) and writes PARENT4.fits: the same table with every Lick/D4000
index corrected for each galaxy's own velocity dispersion
(aux.corrections.C_VD). This is a full, stateless recompute from whatever
is currently in PARENT3.fits -- safe to re-run any time, e.g. after stage 03
has filled in more bins.

Example:
    python3 pipeline/04_cvd_correction.py
"""

import argparse
import sys
from pathlib import Path

from astropy.io import fits
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from desicc.catalog import find_column
from desicc.config import PARENT_FINAL, PARENT_FINAL_CVD
from desicc.cvd import apply_cvd


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, default=PARENT_FINAL)
    parser.add_argument("--output", type=Path, default=PARENT_FINAL_CVD)
    args = parser.parse_args()

    with fits.open(args.input, memmap=False) as hdul:
        meta = Table(hdul["METADATA"].data)
        pyl = Table(hdul["PYLICK"].data)

    vd_col = find_column(meta, "vd")
    corrected = apply_cvd(pyl, meta[vd_col])
    for name, values in corrected.items():
        pyl[name] = values

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out = fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU(meta.as_array(), name="METADATA"),
        fits.BinTableHDU(pyl.as_array(), name="PYLICK"),
    ])
    out.writeto(args.output, overwrite=True)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
