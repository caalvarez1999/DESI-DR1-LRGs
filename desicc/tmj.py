"""Stage 20/21 - TMJ (Thomas, Maraston & Johansson 2011) stellar-population-
synthesis fit to a stack's Lick indices: an MCMC fit of age/metallicity/
[alpha/Fe] against a precomputed model grid (aux/mcmc_tmj.py, ported from
STACKING/MCMC2_MCMC.py), followed by a burn-in cut (aux/burnin.py, ported
from STACKING/MCMC3_burnin_maker.py) and a mode + 16th/84th-percentile
summary of the burned-in chain (aux/posteriors.py, ported from
MYLIBS/posteriors.py).

TMJ_FIT_INDICES (the 13 indices the fit itself uses) is unrelated to the
ARCHAEOLOGY scaling-relation choice used elsewhere in this pipeline, despite
both ultimately coming from an author's variable named "johansson" in two
different, unrelated contexts -- see desicc/config.py.
"""

import re
from pathlib import Path

import numpy as np
from astropy.table import Table

from aux.burnin import burnin
from aux.mcmc_tmj import MCMC_one
from aux.posteriors import allmodalgkde
from .config import TMJ_ERROR_INFLATION, TMJ_FIT_INDICES, TMJ_MODEL_INDICES

STACK_RE = re.compile(r"stack(\d+)")


def stack_number(stack_name: str) -> int:
    m = STACK_RE.search(str(stack_name))
    if not m:
        raise ValueError(f"can't find a stackNNNN number in {stack_name!r}")
    return int(m.group(1))


def chain_basename(*, grid: str, stack_no: int, vd_bin: str, sn_tag: str,
                    archaeology: str, cosmology: str) -> str:
    return (f"chainDESIFULLTMJ{grid}_stack{stack_no:04d}_vd{vd_bin}_SN{sn_tag}"
            f"_{archaeology}_{cosmology}")


def fittable_mask(meta: Table, pyl: Table) -> np.ndarray:
    """True where a stack has >=1 galaxy and every TMJ_FIT_INDICES value/error
    is finite -- same condition 9_FULLTMJ.py used to build `selgood`."""
    ok = np.asarray(meta["ngals"]) >= 1
    for name in TMJ_FIT_INDICES:
        ok &= np.isfinite(np.asarray(pyl[name], dtype=float))
        ok &= np.isfinite(np.asarray(pyl["d" + name], dtype=float))
    return ok


def measured_indices(pyl_row) -> np.ndarray:
    """(25, 2) array of [value, error*TMJ_ERROR_INFLATION], one row per
    TMJ_MODEL_INDICES entry -- the `medidas` MCMC_one expects."""
    out = np.zeros((len(TMJ_MODEL_INDICES), 2))
    for i, name in enumerate(TMJ_MODEL_INDICES):
        out[i, 0] = pyl_row[name]
        out[i, 1] = pyl_row["d" + name] * TMJ_ERROR_INFLATION
    return out


def fit_one_stack(chain_base: Path, medidas: np.ndarray, model: np.ndarray,
                   T: np.ndarray, Z: np.ndarray, A: np.ndarray, *,
                   nwalkers: int, niter: int) -> float:
    """Runs the MCMC (writes chain_base.h5), burns it in (writes
    chain_base_burnin.txt), removes the raw .h5, and returns the fit's chi2
    p-value (a diagnostic, not otherwise used downstream)."""
    P = MCMC_one(
        namefile=str(chain_base), medidas=medidas, elmodelo=model,
        T_ARRAY=T, Z_ARRAY=Z, AFE_ARRAY=A, lista=TMJ_FIT_INDICES,
        nwalkers=nwalkers, Niter=niter, span_pylick=25,
    )
    burnin(str(chain_base), nwalkers=nwalkers)
    Path(str(chain_base) + ".h5").unlink()
    return float(P)


def summarize_one_stack(chain_base: Path, tabla: np.ndarray) -> dict:
    """Reads chain_base_burnin.txt, returns mode/16-84th-percentile spread/
    peak-count for t, Z, afe -- same reduction as 9_FULLTMJ.py's OUT table."""
    out, extra = allmodalgkde(str(chain_base), tabla)
    row = {}
    for ipar, par in enumerate(("t", "Z", "afe")):
        mode, p16, med, p84 = out[ipar, :]
        dlow = abs(mode - p16)
        dup = abs(mode - p84)
        row[par] = float(mode)
        row["d" + par + "_low"] = float(dlow)
        row["d" + par + "_up"] = float(dup)
        row["d" + par] = float(dlow + dup)
        row[par + "_npeaks"] = int(extra["n_peaks"][ipar])
    return row
