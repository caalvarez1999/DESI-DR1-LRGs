#!/usr/bin/env python3
"""Stage 20 - MCMC-fit the TMJ (Thomas, Maraston & Johansson 2011) stellar
population model to every fittable stack in one vd/archaeology/cosmology
group's CVD-corrected summary table (stage 12's output).

For each stack with >=1 galaxy and finite values/errors on every index in
TMJ_FIT_INDICES (desicc/config.py), runs the MCMC fit against the --grid
model (desicc/tmj.fit_one_stack, ported from STACKING/MCMC2_MCMC.py +
STACKING/MCMC3_burnin_maker.py), writing a burned-in chain to
chains/TMJ/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/chainDESIFULLTMJ<grid>_stack
<NNNN>_vd<bin>_SN<sn>_<ARCHAEOLOGY>_<COSMOLOGY>_burnin.txt (the raw .h5 chain
is a large intermediate, deleted once burned in). A stack whose _burnin.txt
already exists is skipped -- safe to re-run any time to pick up new stacks
(stage 10 output) or resume after being interrupted.

Example:
    python3 pipeline/20_SPSfitTMJ.py --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 100
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
    DEFAULT_TMJ_NITER,
    DEFAULT_TMJ_NWALKERS,
    stacked_table_paths,
    tmj_model_path,
    tmj_tabla_path,
)
from desicc.tmj import chain_basename, fit_one_stack, fittable_mask, measured_indices, stack_number


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
    parser.add_argument("--nwalkers", type=int, default=DEFAULT_TMJ_NWALKERS)
    parser.add_argument("--niter", type=int, default=DEFAULT_TMJ_NITER)
    parser.add_argument("--overwrite", action="store_true", help="re-fit stacks that already have a burnin chain")
    args = parser.parse_args()

    sn_tag = f"{round(args.target_sn):03d}"

    input_path = args.input
    if input_path is None:
        _, input_path = stacked_table_paths(args.vd_bin, args.archaeology, args.cosmology, args.target_sn)

    with fits.open(input_path, memmap=False) as hdul:
        meta = Table(hdul["METADATA"].data)
        pyl = Table(hdul["PYLICK"].data)

    ok = fittable_mask(meta, pyl)
    print(f"{int(ok.sum())} / {len(meta)} stacks fittable in {input_path}")
    if not np.any(ok):
        print("Nothing to fit.")
        return

    model = np.load(tmj_model_path(args.grid))
    tabla = np.load(tmj_tabla_path(args.grid))
    T = np.unique(tabla[:, 0])
    Z = np.unique(tabla[:, 1])
    A = np.unique(tabla[:, 2])
    print(f"model grid {args.grid}: {model.shape}, T: {len(T)} values, Z: {len(Z)} values, afe: {len(A)} values")

    chains_dir = args.chains_root / f"vd{args.vd_bin}" / args.archaeology / args.cosmology
    chains_dir.mkdir(parents=True, exist_ok=True)

    idx = np.where(ok)[0]
    for n, i in enumerate(idx, start=1):
        stack_no = stack_number(meta["stack"][i])
        base = chains_dir / chain_basename(
            grid=args.grid, stack_no=stack_no, vd_bin=args.vd_bin, sn_tag=sn_tag,
            archaeology=args.archaeology, cosmology=args.cosmology,
        )
        if not args.overwrite and Path(str(base) + "_burnin.txt").exists():
            print(f"[{n}/{len(idx)}] stack{stack_no:04d}: already fit, skipping")
            continue

        print(f"[{n}/{len(idx)}] stack{stack_no:04d}: fitting...")
        medidas = measured_indices(pyl[i])
        P = fit_one_stack(base, medidas, model, T, Z, A, nwalkers=args.nwalkers, niter=args.niter)
        print(f"[{n}/{len(idx)}] stack{stack_no:04d}: done, chi2 p-value={P:.4g}")

    print(f"\nDone: chains written to {chains_dir}")


if __name__ == "__main__":
    main()
