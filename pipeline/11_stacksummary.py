#!/usr/bin/env python3
"""Stage 11 - collect one vd/archaeology/cosmology group's stack FITS files
(stage 10's output) into a single summary table, the stacked-spectra
equivalent of PARENT3.fits.

Reads every stack<NNNN>.fits under
spectra/stacks/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/ for the given
--target-sn, and writes one row per stack into METADATA (stack filename,
ngals, vd/dvd, z1/dz1, snmedian -- the properties of the stack itself, not
present in any individual-galaxy file) + PYLICK (the stack filename,
replicated from METADATA as the row identifier -- same convention as
targetid in the per-galaxy PYLICK tables -- plus the Lick/D4000 indices
measured on that stack's composite spectrum).

Example:
    python3 pipeline/11_stacksummary.py --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 100
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from desicc.config import (
    ARCHAEOLOGY_CHOICES,
    COSMOLOGY_CHOICES,
    DEFAULT_TARGET_SN,
    LICK_INDICES,
    STACKS_ROOT,
    stacked_table_paths,
)


def read_stack(path: Path) -> dict:
    with fits.open(path, memmap=False) as hdul:
        meta = hdul["METADATA"].data
        pyl = hdul["PYLICK"].data
        row = {
            "stack": path.name,
            "ngals": int(meta["ngals"][0]),
            "vd": float(meta["vd"][0]),
            "dvd": float(meta["dvd"][0]),
            "z1": float(meta["z1"][0]),
            "dz1": float(meta["dz1"][0]),
            "snmedian": float(meta["snmedian"][0]),
        }
        for name in LICK_INDICES:
            row[name] = float(pyl[name][0])
            row["d" + name] = float(pyl["d" + name][0])
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vd-bin", required=True, help='e.g. "200225" for 200 <= vd < 225')
    parser.add_argument("--archaeology", required=True, choices=ARCHAEOLOGY_CHOICES)
    parser.add_argument("--cosmology", required=True, choices=COSMOLOGY_CHOICES)
    parser.add_argument("--target-sn", type=float, default=DEFAULT_TARGET_SN)
    parser.add_argument("--stacks-root", type=Path, default=STACKS_ROOT)
    parser.add_argument("--output", type=Path, default=None,
                         help="default: TABLAS/STACKed/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/"
                              "stackedDESI_vd<bin>_SN<sn>_<ARCHAEOLOGY>_<COSMOLOGY>.fits")
    args = parser.parse_args()

    stack_dir = args.stacks_root / f"vd{args.vd_bin}" / args.archaeology / args.cosmology
    file_prefix = f"SN{args.target_sn:g}"
    stack_files = sorted(stack_dir.glob(f"{file_prefix}_stack*.fits"))

    print(f"{len(stack_files)} stack FITS found in {stack_dir} matching '{file_prefix}_stack*.fits'")
    if not stack_files:
        print("Nothing to summarize.")
        return

    rows = [read_stack(p) for p in stack_files]

    meta = Table()
    meta["stack"] = [r["stack"] for r in rows]
    for col in ("ngals", "vd", "dvd", "z1", "dz1", "snmedian"):
        meta[col] = np.array([r[col] for r in rows])

    pyl = Table()
    pyl["stack"] = meta["stack"]
    for name in LICK_INDICES:
        pyl[name] = np.array([r[name] for r in rows], dtype="f8")
        pyl["d" + name] = np.array([r["d" + name] for r in rows], dtype="f8")

    output = args.output
    if output is None:
        output, _ = stacked_table_paths(args.vd_bin, args.archaeology, args.cosmology, args.target_sn)

    output.parent.mkdir(parents=True, exist_ok=True)
    hdul = fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU(meta.as_array(), name="METADATA"),
        fits.BinTableHDU(pyl.as_array(), name="PYLICK"),
    ])
    hdul.writeto(output, overwrite=True)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
