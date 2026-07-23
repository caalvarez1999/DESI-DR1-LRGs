#!/usr/bin/env python3
"""Stage 01 - look up specid/sparclid/DR for every target via SPARCL.

fastspec-iron doesn't carry SPARCL's own identifiers, so we query SPARCL by
TARGETID to get them -- stage 03 needs sparclid to retrieve each spectrum.
Results are cached to disk so re-running this stage after a partial failure
doesn't re-query objects we already have; pass --force-requery to ignore it.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Column, Table, join

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from desicc.catalog import find_column, to_str_and_valid
from desicc.config import (
    DEFAULT_BATCH_SIZE_SPARCL_FIND,
    PARENT_MERGED,
    PARENT_WITH_IDS,
    SPARCL_IDS_CACHE,
)
from desicc.sparcl_client import batched_find, login


def query_sparcl_ids(targetids, cache_path: Path, force_requery: bool) -> Table:
    if cache_path.exists() and not force_requery:
        print(f"Reading cached SPARCL ids from {cache_path}")
        return Table.read(cache_path)

    client = login()
    print("Querying SPARCL for targetid -> specid/sparcl_id/_dr ...")
    df = batched_find(
        client,
        [str(t) for t in targetids],
        outfields=["targetid", "specid", "sparcl_id", "_dr"],
        batch_size=DEFAULT_BATCH_SIZE_SPARCL_FIND,
    )
    df = df.drop_duplicates(subset="targetid", keep="first")
    df = df.rename(columns={"specid": "specid_sparcl", "sparcl_id": "sparclid_sparcl", "_dr": "DR_sparcl"})

    table = Table.from_pandas(df[["targetid", "specid_sparcl", "sparclid_sparcl", "DR_sparcl"]])
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    table.write(cache_path, overwrite=True)
    print(f"Cached {len(table)} SPARCL ids to {cache_path}")
    return table


def merge_ids_into_parent(parent: Table, sparcl_ids: Table) -> Table:
    joined = join(parent, sparcl_ids, keys="targetid", join_type="left")

    for src, dst, width in (
        ("specid_sparcl", "specid", 17),
        ("sparclid_sparcl", "sparclid", 36),
        ("DR_sparcl", "DR", 8),
    ):
        if src not in joined.colnames:
            continue
        values, valid = to_str_and_valid(joined[src])

        new_col = Column(np.full(len(joined), "", dtype=f"S{width}"), name=dst)
        new_col[valid] = np.array(values[valid], dtype=f"S{width}")
        joined.replace_column(dst, new_col)
        joined.remove_column(src)

    return joined


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, default=PARENT_MERGED)
    parser.add_argument("--output", type=Path, default=PARENT_WITH_IDS)
    parser.add_argument("--cache", type=Path, default=SPARCL_IDS_CACHE)
    parser.add_argument("--force-requery", action="store_true", help="Ignore the cache and query SPARCL again")
    args = parser.parse_args()

    with fits.open(args.input, memmap=False) as hdul:
        parent = Table(hdul["METADATA"].data)
        pylick = hdul["PYLICK"].data

        targetid_col = find_column(parent, "targetid")
        sparcl_ids = query_sparcl_ids(parent[targetid_col], args.cache, args.force_requery)
        merged = merge_ids_into_parent(parent, sparcl_ids)

        out = fits.HDUList([
            fits.PrimaryHDU(),
            fits.BinTableHDU(data=merged.as_array(), name="METADATA"),
            fits.BinTableHDU(data=pylick, name="PYLICK"),
        ])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.writeto(args.output, overwrite=True)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
