"""Core logic for stage 03: download spectra for one velocity-dispersion bin,
compute snmedian, measure Lick indices, and fold the result into PARENT3.fits.

Shared between the interactive and argparse entry points in pipeline/ so
they can't drift apart -- the only difference between them is how they
collect `vd_bin` and `policy` from the user.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

from . import fits_io, lick_measurement
from .catalog import decode, find_column, is_valid_uuid, parse_vd_bin
from .config import (
    DEFAULT_APPLY_RESOLUTION,
    DEFAULT_BATCH_SIZE_SPARCL_RETRIEVE,
    DEFAULT_NSAMPLING,
    DEFAULT_SIGMA_RESOLUTION,
    LICK_INDICES,
    PARENT_DOWNLOADABLE,
    PARENT_FINAL,
    SPARCL_SPECTRUM_FIELDS,
    SINGLEGALAXIES_ROOT,
)
from .sparcl_client import login as sparcl_login

PARENT_COLUMNS = [
    "targetid", "specid", "sparclid", "DR", "survey", "program", "healpix",
    "z1", "vd", "dvd", "FSF_SFR", "FSF_Dn4000", "FSF_dDn4000",
]

# review:  redo only what's missing or fails is_valid_galaxy_fits()
# replace: (re)download everything in the bin, regardless of what's on disk
# skip:    keep any file that already exists, no content check at all
EXISTING_POLICIES = ("review", "replace", "skip")


def load_parent_columns(parent_path: Path) -> tuple[Table, dict]:
    parent = Table.read(parent_path, hdu="METADATA")
    cols = {name: find_column(parent, name) for name in PARENT_COLUMNS}
    return parent, cols


def rows_in_vd_bin(parent: Table, cols: dict, vd_bin: str) -> np.ndarray:
    vd_min, vd_max = parse_vd_bin(vd_bin)
    vd = np.asarray(parent[cols["vd"]], dtype=float)
    return np.where((vd >= vd_min) & (vd < vd_max))[0]


def select_rows_to_process(parent: Table, cols: dict, idxs: np.ndarray,
                            spectra_dir: Path, policy: str) -> list:
    if policy == "replace":
        return [parent[i] for i in idxs]

    needed, kept = [], 0
    for i in idxs:
        tid = int(parent[i][cols["targetid"]])
        path = spectra_dir / f"{tid}.fits"
        if policy == "review":
            ok = path.exists() and fits_io.is_valid_galaxy_fits(path)
        else:  # skip
            ok = path.exists()

        if ok:
            kept += 1
        else:
            needed.append(parent[i])

    print(f"Already on disk (kept as-is): {kept} | to download: {len(needed)}")
    return needed


def download_and_measure(rows: list, cols: dict, out_dir: Path, client, *,
                          batch_size=DEFAULT_BATCH_SIZE_SPARCL_RETRIEVE,
                          nsampling=DEFAULT_NSAMPLING,
                          apply_resolution=DEFAULT_APPLY_RESOLUTION,
                          sigma_resolution=DEFAULT_SIGMA_RESOLUTION):
    if not rows:
        print("Nothing to download.")
        return

    out_dir.mkdir(parents=True, exist_ok=True)

    uuids, uuid_rows = [], []
    for row in rows:
        u = decode(row[cols["sparclid"]])
        if is_valid_uuid(u):
            uuids.append(u)
            uuid_rows.append(row)
        else:
            print(f"Skipping targetid={int(row[cols['targetid']])}: invalid sparclid '{u}'")

    print(f"{len(uuid_rows)}/{len(rows)} rows have a usable sparclid")

    n_written = 0
    for start in range(0, len(uuids), batch_size):
        batch_uuids = uuids[start:start + batch_size]
        batch_rows = uuid_rows[start:start + batch_size]
        print(f"SPARCL batch {start}-{start + len(batch_uuids) - 1} ({len(batch_uuids)} objects)")

        try:
            records = client.retrieve(uuid_list=batch_uuids, include=SPARCL_SPECTRUM_FIELDS)
            df = pd.DataFrame(records.records)
        except Exception as exc:
            print(f"  retrieve() failed: {exc}")
            continue

        if df.empty:
            print("  empty batch returned by SPARCL")
            continue

        by_uuid = {decode(r[cols["sparclid"]]): r for r in batch_rows}
        returned = set()

        for _, spec_row in df.iterrows():
            u = decode(spec_row.get("sparcl_id"))
            returned.add(u)
            parent_row = by_uuid.get(u)
            if parent_row is None:
                continue

            tid = int(parent_row[cols["targetid"]])
            snmedian = fits_io.compute_snmedian(spec_row["flux"], spec_row["ivar"])
            meta_hdu = fits_io.metadata_hdu(parent_row, cols, snmedian)
            spec_hdu = fits_io.spectrum_hdu(spec_row)
            spec_table = spec_hdu.data

            qual, values = lick_measurement.measure(
                np.asarray(spec_table["W"], dtype=float),
                np.asarray(spec_table["F"], dtype=float),
                np.asarray(spec_table["dF"], dtype=float),
                np.asarray(spec_table["S"], dtype=float),
                np.asarray(spec_table["mask"]),
                zgal=float(meta_hdu.data["z1"][0]),
                sn=float(meta_hdu.data["snmedian"][0]),
                nsampling=nsampling,
                apply_resolution=apply_resolution,
                sigma_resolution=sigma_resolution,
            )
            pyl_hdu = fits_io.pylick_hdu(
                tid, qual, values,
                nsampling=nsampling, apply_resolution=apply_resolution, sigma_resolution=sigma_resolution,
            )

            hdul = fits.HDUList([fits.PrimaryHDU(), meta_hdu, spec_hdu, pyl_hdu])
            hdul.writeto(out_dir / f"{tid}.fits", overwrite=True)
            n_written += 1

        missing = [decode(r[cols["sparclid"]]) for r in batch_rows if decode(r[cols["sparclid"]]) not in returned]
        if missing:
            print(f"  SPARCL didn't return {len(missing)} of this batch (e.g. {missing[:5]})")

        print(f"  wrote {n_written} FITS so far")

    print(f"Done: {n_written} files written to {out_dir}")


def load_or_seed_parent_final(path: Path, seed: Path) -> tuple[Table, Table]:
    """Read PARENT3.fits if it already exists, otherwise seed it from PARENT2.fits."""
    src = path if path.exists() else seed
    with fits.open(src, memmap=False) as hdul:
        meta = Table(hdul["METADATA"].data)
        pyl = Table(hdul["PYLICK"].data)
    return meta, pyl


def fold_bin_into_parent(meta: Table, pyl: Table, vd_bin: str, spectra_dir: Path) -> dict:
    """Copy every <targetid>.fits under spectra_dir into the matching row of
    (meta, pyl), in place. Rows already marked done=1 are left untouched, so
    this is safe to re-run across bins and re-runs of the same bin.
    """
    tid_col = find_column(meta, "targetid")
    vd_col = find_column(meta, "vd")
    sn_col = find_column(meta, "snmedian")
    done_col = find_column(pyl, "done")
    qual_col = find_column(pyl, "qual")

    vd_min, vd_max = parse_vd_bin(vd_bin)
    vd = np.asarray(meta[vd_col], dtype=float)
    in_bin = np.where((vd >= vd_min) & (vd < vd_max))[0]

    stats = {"updated": 0, "skipped_done": 0, "skipped_missing": 0, "skipped_bad": 0}

    for row in in_bin:
        if pyl[done_col][row] > 0.99:
            stats["skipped_done"] += 1
            continue

        tid = int(meta[tid_col][row])
        path = spectra_dir / f"{tid}.fits"
        if not path.exists():
            stats["skipped_missing"] += 1
            continue
        if not fits_io.is_valid_galaxy_fits(path):
            stats["skipped_bad"] += 1
            continue

        with fits.open(path, memmap=False) as hdul:
            g_meta = hdul["METADATA"].data
            g_pyl = hdul["PYLICK"].data

            meta[sn_col][row] = g_meta["snmedian"][0]
            pyl[done_col][row] = 1
            pyl[qual_col][row] = g_pyl["qual"][0]
            for name in LICK_INDICES:
                pyl[name][row] = g_pyl[name][0]
                pyl["d" + name][row] = g_pyl["d" + name][0]

        stats["updated"] += 1

    return stats


def write_parent_table(meta: Table, pyl: Table, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    hdul = fits.HDUList([
        fits.PrimaryHDU(),
        fits.BinTableHDU(meta.as_array(), name="METADATA"),
        fits.BinTableHDU(pyl.as_array(), name="PYLICK"),
    ])
    hdul.writeto(path, overwrite=True)


def run(vd_bin: str, policy: str, *, parent_path: Path = PARENT_DOWNLOADABLE,
        spectra_root: Path = SINGLEGALAXIES_ROOT, parent_final_path: Path = PARENT_FINAL,
        **download_kwargs) -> dict:
    if policy not in EXISTING_POLICIES:
        raise ValueError(f"policy must be one of {EXISTING_POLICIES}, got {policy!r}")

    spectra_dir = spectra_root / f"vd{vd_bin}"
    parent, cols = load_parent_columns(parent_path)
    idxs = rows_in_vd_bin(parent, cols, vd_bin)
    print(f"{len(idxs)} targets in vd bin {vd_bin}")

    if len(idxs) == 0:
        print("No targets in this bin, nothing to do.")
    else:
        rows = select_rows_to_process(parent, cols, idxs, spectra_dir, policy)
        if rows:
            client = sparcl_login()
            download_and_measure(rows, cols, spectra_dir, client, **download_kwargs)
        else:
            print("Everything already in place, nothing to do.")

    meta, pyl = load_or_seed_parent_final(parent_final_path, parent_path)
    stats = fold_bin_into_parent(meta, pyl, vd_bin, spectra_dir)
    write_parent_table(meta, pyl, parent_final_path)
    print(f"Wrote {parent_final_path}")

    return stats
