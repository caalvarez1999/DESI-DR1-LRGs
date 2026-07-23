#!/usr/bin/env python3
"""Stage 12 - apply the velocity-dispersion (CVD) correction to a stack
summary table (stage 11's output).

Same correction as stage 04 (aux.corrections.C_VD), applied to the stack's
own mean vd instead of a single galaxy's. A full, stateless recompute --
safe to re-run any time, e.g. after stage 10 has produced more stacks and
stage 11 has been re-run to pick them up.

Example:
    python3 pipeline/12_stackcvd.py --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 100
"""

import argparse
import sys
from pathlib import Path

from astropy.io import fits
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from desicc.catalog import find_column
from desicc.config import ARCHAEOLOGY_CHOICES, COSMOLOGY_CHOICES, DEFAULT_TARGET_SN, stacked_table_paths
from desicc.cvd import apply_cvd


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vd-bin", help='e.g. "200225" for 200 <= vd < 225 -- used to build default paths')
    parser.add_argument("--archaeology", choices=ARCHAEOLOGY_CHOICES)
    parser.add_argument("--cosmology", choices=COSMOLOGY_CHOICES)
    parser.add_argument("--target-sn", type=float, default=DEFAULT_TARGET_SN)
    parser.add_argument("--input", type=Path, default=None, help="overrides the path built from --vd-bin/--archaeology/--cosmology")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    if args.input is None:
        if not (args.vd_bin and args.archaeology and args.cosmology):
            parser.error("either --input, or all of --vd-bin/--archaeology/--cosmology, are required")
        default_input, default_output = stacked_table_paths(args.vd_bin, args.archaeology, args.cosmology, args.target_sn)
        input_path = default_input
        output_path = args.output or default_output
    else:
        input_path = args.input
        output_path = args.output or input_path.with_name(input_path.stem + "_CVD.fits")

    with fits.open(input_path, memmap=False) as hdul:
        meta = Table(hdul["METADATA"].data)
        pyl = Table(hdul["PYLICK"].data)

    vd_col = find_column(meta, "vd")
    corrected = apply_cvd(pyl, meta[vd_col])
    for name, values in corrected.items():
        pyl[name] = values

    output_path.parent.mkdir(parents=True, exist_ok=True)
    out = fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU(meta.as_array(), name="METADATA"),
        fits.BinTableHDU(pyl.as_array(), name="PYLICK"),
    ])
    out.writeto(output_path, overwrite=True)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
