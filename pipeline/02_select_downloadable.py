#!/usr/bin/env python3
"""Stage 02 - keep only the rows that have a usable sparclid.

A left join can leave 'sparclid' empty (or a placeholder like "--") for
targets SPARCL didn't recognize; those can't be retrieved in stage 03, so
they're dropped here. METADATA and PYLICK are filtered together to stay
row-aligned.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from desicc.catalog import decode, find_column, is_placeholder
from desicc.config import PARENT_DOWNLOADABLE, PARENT_WITH_IDS


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, default=PARENT_WITH_IDS)
    parser.add_argument("--output", type=Path, default=PARENT_DOWNLOADABLE)
    args = parser.parse_args()

    with fits.open(args.input, memmap=False) as hdul:
        meta = Table(hdul["METADATA"].data)
        pylick = Table(hdul["PYLICK"].data)

        sparclid_col = find_column(meta, "sparclid")
        sparclid = np.array([decode(s) for s in meta[sparclid_col]])
        good = ~np.array([is_placeholder(s) for s in sparclid])

        n_total = len(meta)
        print(f"Downloadable: {int(good.sum())} / {n_total}  (discarded: {n_total - int(good.sum())})")

        out = fits.HDUList([
            fits.PrimaryHDU(),
            fits.BinTableHDU(data=meta[good].as_array(), name="METADATA"),
            fits.BinTableHDU(data=pylick[good].as_array(), name="PYLICK"),
        ])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.writeto(args.output, overwrite=True)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
