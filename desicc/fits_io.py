"""Build and validate the per-galaxy FITS files (METADATA + SPECTRUM + PYLICK)."""

from datetime import datetime, UTC

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

from .catalog import decode, normalize_specid
from .config import LICK_INDICES


def compute_snmedian(flux, ivar) -> float:
    flux = np.asarray(flux, dtype=float)
    ivar = np.asarray(ivar, dtype=float)

    dflux = np.full_like(ivar, np.nan)
    good_ivar = ivar > 0
    dflux[good_ivar] = ivar[good_ivar] ** -0.5

    valid = np.isfinite(flux) & np.isfinite(dflux) & (dflux > 0)
    if not np.any(valid):
        return float("nan")
    return float(np.nanmedian(flux[valid] / dflux[valid]))


def spectrum_hdu(sparcl_row: pd.Series) -> fits.BinTableHDU:
    wavelength = np.asarray(sparcl_row["wavelength"], dtype=float)
    flux = np.asarray(sparcl_row["flux"], dtype=float)
    ivar = np.asarray(sparcl_row["ivar"], dtype=float)
    mask = np.asarray(sparcl_row["mask"])
    wave_sigma = np.asarray(sparcl_row["wave_sigma"], dtype=float)

    dflux = np.full_like(ivar, np.nan)
    good_ivar = ivar > 0
    dflux[good_ivar] = ivar[good_ivar] ** -0.5

    table = Table()
    table["W"] = wavelength
    table["F"] = flux
    table["dF"] = dflux
    table["mask"] = mask
    table["S"] = wave_sigma
    return fits.BinTableHDU(table, name="SPECTRUM")


def metadata_hdu(parent_row, cols: dict, snmedian: float) -> fits.BinTableHDU:
    table = Table()
    table["targetid"] = [int(parent_row[cols["targetid"]])]
    table["specid"] = [normalize_specid(parent_row[cols["specid"]])]
    table["sparclid"] = [decode(parent_row[cols["sparclid"]])]
    table["healpix"] = [int(parent_row[cols["healpix"]])]
    table["program"] = [decode(parent_row[cols["program"]])]
    table["survey"] = [decode(parent_row[cols["survey"]])]
    table["DR"] = [decode(parent_row[cols["DR"]])]
    table["FSF_SFR"] = [float(parent_row[cols["FSF_SFR"]])]
    table["FSF_Dn4000"] = [float(parent_row[cols["FSF_Dn4000"]])]
    table["FSF_dDn4000"] = [float(parent_row[cols["FSF_dDn4000"]])]
    table["z1"] = [float(parent_row[cols["z1"]])]
    table["vd"] = [float(parent_row[cols["vd"]])]
    table["dvd"] = [float(parent_row[cols["dvd"]])]
    table["snmedian"] = [float(snmedian)]
    return fits.BinTableHDU(table, name="METADATA")


def pylick_hdu(targetid: int, qual: int, values: dict, *, nsampling: int,
               apply_resolution: bool, sigma_resolution: float) -> fits.BinTableHDU:
    table = Table()
    table["targetid"] = [int(targetid)]
    table["qual"] = [int(qual)]
    for name in LICK_INDICES:
        table[name] = [float(values.get(name, np.nan))]
        table["d" + name] = [float(values.get("d" + name, np.nan))]

    hdu = fits.BinTableHDU(table, name="PYLICK")
    hdu.header["PYL_TS"] = (datetime.now(UTC).isoformat(timespec="seconds"), "UTC timestamp")
    hdu.header["NSAMPL"] = (int(nsampling), "measure_indices nsampling")
    hdu.header["RESOL"] = (bool(apply_resolution), "Resolution correction applied")
    hdu.header["SIGRES"] = (float(sigma_resolution), "Sigma resolution (Angstrom), if RESOL")
    return hdu


def is_valid_galaxy_fits(path) -> bool:
    """True if `path` has the HDUs and columns the rest of the pipeline relies on.

    Used by the "review" existing-file policy in stage 03 to tell a
    genuinely finished file apart from one left over from a crashed or
    half-written run.
    """
    try:
        with fits.open(path, memmap=False) as hdul:
            names = {h.name.upper() for h in hdul if getattr(h, "name", None)}
            if not {"METADATA", "SPECTRUM", "PYLICK"} <= names:
                return False

            meta_cols = set(hdul["METADATA"].data.columns.names)
            if not {"targetid", "z1", "snmedian"} <= meta_cols:
                return False

            spec_cols = set(hdul["SPECTRUM"].data.columns.names)
            if not {"W", "F", "dF", "mask", "S"} <= spec_cols:
                return False

            pyl_cols = set(hdul["PYLICK"].data.columns.names)
            if not {"targetid", "qual"} <= pyl_cols:
                return False
            for name in LICK_INDICES:
                if name not in pyl_cols or f"d{name}" not in pyl_cols:
                    return False
        return True
    except Exception:
        return False
