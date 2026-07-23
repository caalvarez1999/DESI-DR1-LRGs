#!/usr/bin/env python3
"""Stage 10 - build S/N-targeted stacks for one vd/archaeology/cosmology group.

Selects the single galaxies (from PARENT4.fits, i.e. the CVD-corrected
catalog) eligible for this group -- in the vd bin, below the archaeology
+cosmology redshift cut, done/qual OK, no significant emission, CaII H/K
ratio OK (desicc/stack_selection.py) -- sorts them by redshift, and consumes
them into consecutive stacks that each reach --target-sn (desicc/stacking.py).

Output goes to spectra/stacks/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/, one
FITS per stack (METADATA + SPECTRUM + PYLICK). This is a long-running
process -- it checkpoints after every stack to a resume file in the same
directory, and picks up where it left off if interrupted and re-run.

Example:
    python3 pipeline/10_stackmaker.py --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 100
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
    PARENT_FINAL_CVD,
    SINGLEGALAXIES_ROOT,
    STACKS_ROOT,
)
from desicc.stack_selection import load_zcut, select_for_stacking
from desicc.stacking import run as run_stacking


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vd-bin", required=True, help='e.g. "200225" for 200 <= vd < 225')
    parser.add_argument("--archaeology", required=True, choices=ARCHAEOLOGY_CHOICES)
    parser.add_argument("--cosmology", required=True, choices=COSMOLOGY_CHOICES)
    parser.add_argument("--target-sn", type=float, default=DEFAULT_TARGET_SN)
    parser.add_argument("--parent", type=Path, default=PARENT_FINAL_CVD,
                         help="CVD-corrected catalog to select single galaxies from (default: %(default)s)")
    parser.add_argument("--singlegalaxies-root", type=Path, default=SINGLEGALAXIES_ROOT)
    parser.add_argument("--stacks-root", type=Path, default=STACKS_ROOT)
    args = parser.parse_args()

    zcut = load_zcut(args.cosmology, args.archaeology, args.vd_bin)
    print(f"zcut for vd{args.vd_bin} / {args.archaeology} / {args.cosmology}: z <= {zcut:.4f}")

    with fits.open(args.parent, memmap=False) as hdul:
        meta = Table(hdul["METADATA"].data)
        pyl = Table(hdul["PYLICK"].data)

    sel = select_for_stacking(meta, pyl, args.vd_bin, zcut)
    print(f"{int(sel.sum())} / {len(meta)} galaxies eligible for this stack")
    if not np.any(sel):
        print("Nothing eligible, nothing to do.")
        return

    rows = meta[sel]
    rows = rows[np.argsort(rows["z1"])]

    spectra_dir = args.singlegalaxies_root / f"vd{args.vd_bin}"
    out_dir = args.stacks_root / f"vd{args.vd_bin}" / args.archaeology / args.cosmology
    file_prefix = f"SN{args.target_sn:g}"
    resume_path = out_dir / f"{file_prefix}_resume.json"

    stats = run_stacking(
        rows, spectra_dir, out_dir,
        target_sn=args.target_sn, resume_path=resume_path, file_prefix=file_prefix,
    )

    print(f"\nDone: {stats['n_stacks']} stacks written to {out_dir}")


if __name__ == "__main__":
    main()
