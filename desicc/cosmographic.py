"""Stage 30 - cosmographic (Taylor-expansion) fit of t(z) combining several
vd groups' TMJ ages (stage 21's output) to get H(z0), q(z0), j(z0) plus a
per-group age-at-z0 nuisance parameter.

Ported from ARTICLE3/13_tz_cosmographic.py's single (non z-binned) fit path:
reads each selected vd group's TMJ posteriors table, stacks them, applies
the same reliability cuts, and MCMC-fits the Taylor expansion
(aux/mcmc_cosmographic.py, aux/taylor_tz[_wide].py -- ported from
ARTICLE3/aux13_MCMC.py + MYLIBS/TAYLOR_emcee_tz[_widepriors].py), then cuts
the burn-in (aux/burnin_cosmographic.py, ported from
ARTICLE3/aux13_burnin_maker.py).
"""

import re
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table, hstack, vstack

from .config import modelled_table_path

GROUP_RE = re.compile(r"_vd(\d+)_")


def group_of(stack_name: str) -> str:
    m = GROUP_RE.search(str(stack_name))
    if not m:
        raise ValueError(f"can't find a vd<bin> group in {stack_name!r}")
    return m.group(1)


def var_names(groups) -> list:
    return [f"age{g}" for g in groups] + ["Hz0", "qz0", "jz0"]


def load_group_tables(vd_bins, archaeology: str, cosmology: str, target_sn: float, grid: str) -> Table:
    """Reads and vstacks each vd group's stage 21 TMJ posteriors table
    (METADATA + TMJ, row-aligned) into one flat table."""
    tables = []
    for vd_bin in vd_bins:
        path = modelled_table_path(vd_bin, archaeology, cosmology, target_sn, grid)
        with fits.open(path, memmap=False) as hdul:
            meta = Table(hdul["METADATA"].data)
            tmj = Table(hdul["TMJ"].data)
        tmj_no_stack = tmj[[c for c in tmj.colnames if c != "stack"]]
        tables.append(hstack([meta, tmj_no_stack]))
    return vstack(tables)


def select_fittable(res: Table) -> Table:
    """Same reliability cuts as 13_tz_cosmographic.py's RESRES."""
    ok = (
        (res["valid"] > 0.0)
        & (res["z1"] / res["dz1"] > 50.0)
        & (res["dz1"] < 0.01)
        & (res["t"] / res["dt"] > 5.0)
        & (res["dt"] < 1.0)
    )
    return res[ok]


def fit(res_sel: Table, *, z0: float, priors: str, order: int, name: str,
        nwalkers: int, niter: int, chains_dir: Path) -> Path:
    """Runs the MCMC (writes {chains_dir}/{name}.h5), burns it in (writes
    {name}_burnin.txt), removes the raw .h5 and un-burned .txt, and returns
    the burned-in chain's path."""
    from aux.burnin_cosmographic import burnin
    from aux.mcmc_cosmographic import MCMC_one

    groups = [group_of(s) for s in res_sel["stack"]]
    dynamic_var_nms = var_names(sorted(set(groups)))

    chains_dir = Path(chains_dir)
    chains_dir.mkdir(parents=True, exist_ok=True)

    MCMC_one(
        np.asarray(res_sel["z1"], dtype=float), z0,
        np.asarray(res_sel["t"], dtype=float), np.asarray(res_sel["dt"], dtype=float),
        groups, priors=priors, ORDER=order, name=name, nwalkers=nwalkers, Niter=niter,
        var_nms=dynamic_var_nms, chains_dir=str(chains_dir) + "/",
    )

    base = chains_dir / name
    burnin(str(base), txt="yes", ndim=len(dynamic_var_nms), nwalkers=nwalkers)
    Path(str(base) + ".h5").unlink()
    Path(str(base) + ".txt").unlink()
    return Path(str(base) + "_burnin.txt")
