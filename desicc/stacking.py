"""Core logic for stage 10: build S/N-targeted stacks of single-galaxy
spectra within one vd/archaeology/cosmology group.

Galaxies (already selected by desicc.stack_selection and sorted by redshift)
are consumed in order: starting from the bluest remaining galaxy, nfinder()
binary-searches for the smallest number of consecutive galaxies whose
composite spectrum reaches the target S/N, then composite_final_miles()
rebuilds that same group as the final stack -- normalizing each spectrum in
the rest-frame window, degrading each one to MILES resolution (C_resol) only
at this final stage (skipped during the S/N search, for speed), and
median-combining. Lick indices are then measured on the resulting spectrum.

This is a long-running, resumable process: run() checkpoints (how many
galaxies have been consumed, how many stacks written) to a JSON file after
every stack, and picks the checkpoint back up on the next call, so it's safe
to interrupt and re-run.
"""

import json
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table

from aux.corrections import C_resol, C_toair
from aux.emcee_lick import I_numberopen
from aux.indices_singlegalaxies import measure_indices

from .config import DEFAULT_SIGMA_RESOLUTION, LICK_INDICES, RESTFRAME_WINDOW, STACK_NORM_WINDOW_RF


def _read_single_spectrum(path: Path):
    with fits.open(path, memmap=False) as hdul:
        data = hdul["SPECTRUM"].data
        w = np.asarray(data["W"], dtype=float)
        f = np.asarray(data["F"], dtype=float)
        df = np.asarray(data["dF"], dtype=float)
        s = np.asarray(data["S"], dtype=float)
        mask = np.asarray(data["mask"])

    valid = (mask == 0) & np.isfinite(w) & np.isfinite(f) & np.isfinite(df)
    return w, f, df, s, valid


def _normalized_restframe(w, f, df, s, valid, z: float, wave_grid: np.ndarray):
    w, f, df, s = w[valid], f[valid], df[valid], s[valid]
    w_rf = C_toair(w) / (1.0 + z)

    lo, hi = STACK_NORM_WINDOW_RF
    norm_window = (w_rf > lo) & (w_rf < hi)
    if np.sum(norm_window) < 3:
        return None

    med = np.nanmedian(f[norm_window])
    if not np.isfinite(med) or med <= 0:
        return None

    f_i = np.interp(wave_grid, w_rf, f / med, left=np.nan, right=np.nan)
    df_i = np.interp(wave_grid, w_rf, df / med, left=np.nan, right=np.nan)
    s_i = np.interp(wave_grid, w_rf, s, left=np.nan, right=np.nan)
    return f_i, df_i, s_i


def quick_sn_estimate(rows: Table, spectra_dir: Path, *, dl: float = 1.0) -> float:
    """Median S/N of the quick (undegraded) composite of `rows` -- used to
    drive nfinder's search, not the final stack."""
    lo, hi = RESTFRAME_WINDOW
    wave_grid = np.arange(lo, hi + dl, dl)

    fluxes = []
    for row in rows:
        path = spectra_dir / f"{int(row['targetid'])}.fits"
        if not path.exists():
            continue
        try:
            w, f, df, s, valid = _read_single_spectrum(path)
        except Exception:
            continue
        if np.sum(valid) < 10:
            continue

        result = _normalized_restframe(w, f, df, s, valid, float(row["z1"]), wave_grid)
        if result is None:
            continue
        fluxes.append(result[0])

    if len(fluxes) < 2:
        return float(rows["snmedian"][0]) if len(rows) else float("nan")

    fluxes = np.asarray(fluxes)
    f_stack = np.nanmedian(fluxes, axis=0)
    n = len(fluxes)
    df_stack = 1.4826 * np.nanmedian(np.abs(fluxes - f_stack), axis=0) / np.sqrt(n)

    ok = np.isfinite(f_stack) & np.isfinite(df_stack) & (df_stack > 0)
    if not np.any(ok):
        return 0.001
    return float(np.nanmedian(f_stack[ok] / df_stack[ok]))


def nfinder(target_sn: float, rows: Table, start: int, spectra_dir: Path, *,
            max_iters: int = 200, tol: float = 1.0) -> np.ndarray:
    """Smallest run of consecutive rows (from `start`) whose quick composite
    reaches target_sn, via doubling + bisection. Returns row indices into
    `rows` (not `start`-relative)."""
    max_extend = len(rows) - start

    def sn_of(n: int) -> float:
        return quick_sn_estimate(rows[start:min(start + n, len(rows))], spectra_dir)

    n_lo = 1
    sn = sn_of(1)
    if not np.isfinite(sn):
        return np.arange(start, start + 1)

    n_hi = n_lo
    while sn < target_sn and n_hi < max_extend:
        n_lo = n_hi
        n_hi = min(max_extend, max(2, n_hi * 2))
        sn = sn_of(n_hi)
        if not np.isfinite(sn):
            break

    n_mid = n_hi
    for _ in range(max_iters):
        n_mid = (n_lo + n_hi) // 2
        sn_mid = sn_of(n_mid)
        if not np.isfinite(sn_mid):
            n_hi = max(1, n_mid - 1)
            continue
        if abs(sn_mid - target_sn) < tol or n_hi - n_lo <= 1:
            break
        if sn_mid < target_sn:
            n_lo = n_mid
        else:
            n_hi = n_mid

    return np.arange(start, min(start + max(1, n_mid), len(rows)))


def composite_final_miles(rows: Table, spectra_dir: Path, *,
                           sigma_resolution: float = DEFAULT_SIGMA_RESOLUTION):
    """Build the final stack: normalize + resolution-degrade + median-combine
    every galaxy in `rows`. Returns (W, F, dF, S, mask, used_targetids)."""
    lo, hi = RESTFRAME_WINDOW
    wave_grid = np.arange(lo, hi + 0.1, 0.1)

    fluxes, sigmas, used_ids = [], [], []
    for row in rows:
        tid = int(row["targetid"])
        path = spectra_dir / f"{tid}.fits"
        if not path.exists():
            continue
        try:
            w, f, df, s, valid = _read_single_spectrum(path)
        except Exception:
            continue
        if np.sum(valid) < 10:
            continue

        result = _normalized_restframe(w, f, df, s, valid, float(row["z1"]), wave_grid)
        if result is None:
            continue
        f_i, _df_i, s_i = result
        if np.all(np.isnan(f_i)):
            continue

        fluxes.append(f_i)
        sigmas.append(s_i)
        used_ids.append(tid)

    n_used = len(fluxes)
    if n_used == 0:
        n = len(wave_grid)
        empty = np.full(n, np.nan)
        return wave_grid, empty, empty, empty, np.zeros(n, dtype=bool), np.array([], dtype=np.int64)

    fluxes = np.asarray(fluxes)
    sigmas = np.asarray(sigmas)

    f_stack = np.nanmedian(fluxes, axis=0)
    s_stack = np.nanmedian(sigmas, axis=0)
    df_stack = 1.4826 * np.nanmedian(np.abs(fluxes - f_stack), axis=0) / np.sqrt(n_used)

    n_contrib = np.sum(np.isfinite(fluxes), axis=0)
    mask = (n_contrib > 0) & np.isfinite(f_stack) & np.isfinite(df_stack) & np.isfinite(s_stack) & (df_stack > 0)

    s_filled = np.where(np.isfinite(s_stack), s_stack, sigma_resolution)
    _, f_deg, df_deg = C_resol(wave_grid, 0.0, f_stack, df_stack, s_filled, sigma_resolution)

    return wave_grid, f_deg, df_deg, s_filled, mask, np.asarray(used_ids, dtype=np.int64)


def measure_stack_indices(w, f, df, mask, *, nsampling: int = 1000):
    ok = mask & np.isfinite(w) & np.isfinite(f) & np.isfinite(df) & (df > 0)
    w_ok, f_ok, df_ok = w[ok], f[ok], df[ok]
    order = np.argsort(w_ok)
    w_ok, f_ok, df_ok = w_ok[order], f_ok[order], df_ok[order]

    index_ids = I_numberopen(LICK_INDICES)
    values, errors, _, _ = measure_indices(
        w_ok, f_ok, df_ok, SN=False, truesampling=True, nsampling=nsampling, index_list=index_ids,
    )
    return np.asarray(values, dtype=float), np.asarray(errors, dtype=float)


def build_stack_hdulist(*, wave, flux, dflux, sigma, mask, ngals, vd, dvd, z1, dz1,
                         snmedian, targetids, index_values, index_errors) -> fits.HDUList:
    meta = Table()
    meta["ngals"] = [int(ngals)]
    meta["vd"] = [float(vd)]
    meta["dvd"] = [float(dvd)]
    meta["z1"] = [float(z1)]
    meta["dz1"] = [float(dz1)]
    meta["snmedian"] = [float(snmedian)]
    meta["targetids"] = [json.dumps([int(t) for t in targetids])]
    meta_hdu = fits.BinTableHDU(meta, name="METADATA")

    spec = Table()
    spec["W"] = np.asarray(wave, dtype=float)
    spec["F"] = np.asarray(flux, dtype=float)
    spec["dF"] = np.asarray(dflux, dtype=float)
    spec["S"] = np.asarray(sigma, dtype=float)
    spec["mask"] = np.asarray(mask, dtype=bool)
    spec_hdu = fits.BinTableHDU(spec, name="SPECTRUM")

    pyl = Table()
    for name, val, err in zip(LICK_INDICES, index_values, index_errors):
        pyl[name] = [float(val)]
        pyl["d" + name] = [float(err)]
    pyl_hdu = fits.BinTableHDU(pyl, name="PYLICK")

    return fits.HDUList([fits.PrimaryHDU(), meta_hdu, spec_hdu, pyl_hdu])


def run(rows: Table, spectra_dir: Path, out_dir: Path, *, target_sn: float,
        resume_path: Path, file_prefix: str = "stack",
        sigma_resolution: float = DEFAULT_SIGMA_RESOLUTION,
        nsampling: int = 1000) -> dict:
    """Consume `rows` (must already be sorted by z1) into stacks until
    exhausted, writing one FITS per stack to out_dir and checkpointing to
    resume_path after each one. Safe to interrupt and re-call."""
    out_dir.mkdir(parents=True, exist_ok=True)

    if resume_path.exists():
        state = json.loads(resume_path.read_text())
        n = int(state.get("n", 0))
        n_stacks = int(state.get("n_stacks", 0))
        print(f"Resuming: n={n}, n_stacks={n_stacks}")
    else:
        n, n_stacks = 0, 0

    t0 = time.time()
    while n < len(rows):
        idxs = nfinder(target_sn, rows, n, spectra_dir)
        if len(idxs) < 2:
            n += 1
            continue

        subset = rows[idxs]
        wave, flux, dflux, sigma, mask, used_ids = composite_final_miles(
            subset, spectra_dir, sigma_resolution=sigma_resolution,
        )
        if len(used_ids) < 2:
            # NMAD uncertainty (and hence the valid-pixel mask) is undefined
            # with fewer than 2 contributing spectra -- skip and retry.
            n += 1
            continue

        used_rows = subset[np.isin(np.asarray(subset["targetid"]), used_ids)]

        sn_pix = np.full_like(flux, np.nan)
        ok = mask & np.isfinite(flux) & np.isfinite(dflux) & (dflux > 0)
        sn_pix[ok] = flux[ok] / dflux[ok]
        snmedian = float(np.nanmedian(sn_pix[ok])) if np.any(ok) else float("nan")

        index_values, index_errors = measure_stack_indices(wave, flux, dflux, mask, nsampling=nsampling)

        n_stacks += 1
        stack_path = out_dir / f"{file_prefix}_stack{n_stacks:04d}.fits"
        hdul = build_stack_hdulist(
            wave=wave, flux=flux, dflux=dflux, sigma=sigma, mask=mask,
            ngals=len(used_ids),
            vd=float(np.nanmean(used_rows["vd"])), dvd=float(np.nanstd(used_rows["vd"])),
            z1=float(np.nanmean(used_rows["z1"])), dz1=float(np.nanstd(used_rows["z1"])),
            snmedian=snmedian, targetids=used_ids,
            index_values=index_values, index_errors=index_errors,
        )
        hdul.writeto(stack_path, overwrite=True)

        n = int(idxs[-1] + 1)
        resume_path.write_text(json.dumps({"n": n, "n_stacks": n_stacks}))
        print(f"stack {n_stacks} -> {stack_path.name} | ngals={len(used_ids)} | "
              f"S/N~{snmedian:.1f} | {n}/{len(rows)} galaxies consumed | "
              f"{time.time() - t0:.0f}s elapsed")

    resume_path.unlink(missing_ok=True)
    return {"n_stacks": n_stacks, "n_galaxies": len(rows)}
