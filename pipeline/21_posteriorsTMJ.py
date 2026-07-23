#!/usr/bin/env python3
"""Stage 21 - collect one vd/archaeology/cosmology group's TMJ MCMC chains
(stage 20's output) into a single posteriors table, adding a TMJ HDU on top
of the METADATA+PYLICK structure carried over from the CVD-corrected stack
summary table (stage 12's output).

For every stack, if it was fittable and its burned-in chain exists
(chains/TMJ/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/chainDESIFULLTMJ<grid>_stack
<NNNN>_..._burnin.txt), reduces it to a mode + 16th/84th-percentile spread
per parameter (desicc/tmj.summarize_one_stack, ported from
MYLIBS/posteriors.py's allmodalgkde). Otherwise the row is filled with NaN
and TMJ.valid=0, same convention as the original 9_FULLTMJ.py.

Example:
    python3 pipeline/21_posteriorsTMJ.py --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 100
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
    CHAINS_ROOT,
    COSMOLOGY_CHOICES,
    DEFAULT_TARGET_SN,
    DEFAULT_TMJ_GRID,
    modelled_table_path,
    stacked_table_paths,
    tmj_tabla_path,
)
from desicc.tmj import chain_basename, fittable_mask, stack_number, summarize_one_stack

TMJ_PARAMS = ("t", "Z", "afe")
TMJ_COLS = []
for _p in TMJ_PARAMS:
    TMJ_COLS += [_p, f"d{_p}_low", f"d{_p}_up", f"d{_p}", f"{_p}_npeaks"]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vd-bin", required=True, help='e.g. "200225" for 200 <= vd < 225')
    parser.add_argument("--archaeology", required=True, choices=ARCHAEOLOGY_CHOICES)
    parser.add_argument("--cosmology", required=True, choices=COSMOLOGY_CHOICES)
    parser.add_argument("--target-sn", type=float, default=DEFAULT_TARGET_SN)
    parser.add_argument("--grid", default=DEFAULT_TMJ_GRID, help="TMJ model grid tag, e.g. custom250924b")
    parser.add_argument("--input", type=Path, default=None,
                         help="CVD-corrected stack summary table (default: stage 12's default output path)")
    parser.add_argument("--chains-root", type=Path, default=CHAINS_ROOT)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    sn_tag = f"{round(args.target_sn):03d}"

    input_path = args.input
    if input_path is None:
        _, input_path = stacked_table_paths(args.vd_bin, args.archaeology, args.cosmology, args.target_sn)

    with fits.open(input_path, memmap=False) as hdul:
        meta = Table(hdul["METADATA"].data)
        pyl = Table(hdul["PYLICK"].data)

    fittable = fittable_mask(meta, pyl)
    chains_dir = args.chains_root / f"vd{args.vd_bin}" / args.archaeology / args.cosmology
    tabla = None

    rows = []
    n_valid = 0
    for i in range(len(meta)):
        stack_no = stack_number(meta["stack"][i])
        base = chains_dir / chain_basename(
            grid=args.grid, stack_no=stack_no, vd_bin=args.vd_bin, sn_tag=sn_tag,
            archaeology=args.archaeology, cosmology=args.cosmology,
        )
        has_chain = Path(str(base) + "_burnin.txt").exists()
        row = {"stack": meta["stack"][i], "valid": int(bool(fittable[i]) and has_chain)}
        if row["valid"]:
            if tabla is None:
                tabla = np.load(tmj_tabla_path(args.grid))
            row.update(summarize_one_stack(base, tabla))
            n_valid += 1
        else:
            for col in TMJ_COLS:
                row[col] = np.nan
        rows.append(row)

    print(f"{n_valid} / {len(meta)} stacks have a valid TMJ posterior")

    tmj = Table()
    tmj["stack"] = [r["stack"] for r in rows]
    for col in TMJ_COLS:
        # *_npeaks stays float: it's NaN (not 0) for invalid/unfit rows,
        # same convention as the original 9_FULLTMJ.py.
        tmj[col] = np.array([r[col] for r in rows], dtype="f8")
    tmj["valid"] = np.array([r["valid"] for r in rows], dtype="i8")

    output = args.output
    if output is None:
        output = modelled_table_path(args.vd_bin, args.archaeology, args.cosmology, args.target_sn, args.grid)

    output.parent.mkdir(parents=True, exist_ok=True)
    hdul = fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU(meta.as_array(), name="METADATA"),
        fits.BinTableHDU(pyl.as_array(), name="PYLICK"),
        fits.BinTableHDU(tmj.as_array(), name="TMJ"),
    ])
    hdul.writeto(output, overwrite=True)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
