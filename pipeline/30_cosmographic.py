#!/usr/bin/env python3
"""Stage 30 - MCMC-fit the cosmographic (Taylor-expansion) t(z) model,
combining several vd groups' TMJ posteriors tables (stage 21's output) into
a joint H(z0)/q(z0)/j(z0) + per-group age-at-z0 fit.

Reads TABLAS/modelled/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/posteriorsDESI_TMJ_
<grid>_..._SN<sn>.fits for every --vd-bin given, stacks them, keeps only
reliable rows (valid stack, z1/dz1>50, dz1<0.01, t/dt>5, dt<1 --
desicc.cosmographic.select_fittable, same cuts as 13_tz_cosmographic.py),
and MCMC-fits a --order Taylor expansion of t(z) around --z0 (default: the
sample's median z1) for all selected groups jointly
(desicc/cosmographic.py, aux/mcmc_cosmographic.py, aux/taylor_tz[_wide].py
-- ported from ARTICLE3/13_tz_cosmographic.py + aux13_MCMC.py +
MYLIBS/TAYLOR_emcee_tz[_widepriors].py). The burned-in chain is written to
chains/cosmographic/<name>_burnin.txt; a chain whose file already exists is
skipped unless --overwrite is given.

Example:
    python3 pipeline/30_cosmographic.py --vd-bins 200225 225250 250280 280320 320355 \\
        --archaeology ALVAREZ --cosmology PLANCK --target-sn 150 --priors wide
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from desicc.config import (
    ARCHAEOLOGY_CHOICES,
    COSMOGRAPHIC_CHAINS_ROOT,
    COSMOGRAPHIC_PRIORS_CHOICES,
    COSMOLOGY_CHOICES,
    DEFAULT_COSMOGRAPHIC_NITER,
    DEFAULT_COSMOGRAPHIC_NWALKERS,
    DEFAULT_COSMOGRAPHIC_ORDER,
    DEFAULT_COSMOGRAPHIC_PRIORS,
    DEFAULT_TARGET_SN,
    DEFAULT_TMJ_GRID,
)
from desicc.cosmographic import fit, load_group_tables, select_fittable


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vd-bins", nargs="+", required=True,
                         help='vd groups to fit jointly, e.g. "200225 225250 250280 280320 320355"')
    parser.add_argument("--archaeology", required=True, choices=ARCHAEOLOGY_CHOICES)
    parser.add_argument("--cosmology", required=True, choices=COSMOLOGY_CHOICES)
    parser.add_argument("--target-sn", type=float, default=DEFAULT_TARGET_SN)
    parser.add_argument("--grid", default=DEFAULT_TMJ_GRID, help="TMJ model grid tag, e.g. custom250924b")
    parser.add_argument("--priors", choices=COSMOGRAPHIC_PRIORS_CHOICES, default=DEFAULT_COSMOGRAPHIC_PRIORS)
    parser.add_argument("--order", type=int, default=DEFAULT_COSMOGRAPHIC_ORDER)
    parser.add_argument("--z0", type=float, default=None, help="default: median z1 of the fittable sample")
    parser.add_argument("--nwalkers", type=int, default=DEFAULT_COSMOGRAPHIC_NWALKERS)
    parser.add_argument("--niter", type=int, default=DEFAULT_COSMOGRAPHIC_NITER)
    parser.add_argument("--chains-root", type=Path, default=COSMOGRAPHIC_CHAINS_ROOT)
    parser.add_argument("--name", default=None, help="chain basename (default: built from grid/vd-bins/SN/archaeology/cosmology/priors)")
    parser.add_argument("--overwrite", action="store_true", help="re-fit even if the burnin chain already exists")
    args = parser.parse_args()

    sn_tag = f"{round(args.target_sn):03d}"

    res = load_group_tables(args.vd_bins, args.archaeology, args.cosmology, args.target_sn, args.grid)
    print(f"{len(res)} stacks across {len(args.vd_bins)} vd groups")

    sel = select_fittable(res)
    print(f"{len(sel)} / {len(res)} stacks pass the reliability cuts")
    if len(sel) <= 1:
        print("Not enough data to fit.")
        return

    z0 = args.z0 if args.z0 is not None else float(np.median(sel["z1"]))
    print(f"z0 = {z0:.6f}")

    name = args.name
    if name is None:
        vd_label = "-".join(args.vd_bins)
        name = f"cosmographicTMJ{args.grid}_vd{vd_label}_SN{sn_tag}_{args.archaeology}_{args.cosmology}_{args.priors}"

    burnin_path = args.chains_root / f"{name}_burnin.txt"
    if not args.overwrite and burnin_path.exists():
        print(f"{burnin_path} already exists, skipping (use --overwrite to re-fit).")
        return

    out_path = fit(
        sel, z0=z0, priors=args.priors, order=args.order, name=name,
        nwalkers=args.nwalkers, niter=args.niter, chains_dir=args.chains_root,
    )
    print(f"\nDone: burned-in chain written to {out_path}")


if __name__ == "__main__":
    main()
