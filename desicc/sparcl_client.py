"""Thin wrapper around the SPARCL client: login, plus batched/retried queries.

Kept separate from the pipeline logic so both `find` (stage 01, id lookup)
and `retrieve` (stage 03, spectrum download) share the same login and retry
behaviour instead of each script reimplementing it slightly differently.
"""

import os
import time

import pandas as pd
from sparcl.client import SparclClient

from .config import DEFAULT_MAX_RETRIES, DEFAULT_RETRY_SLEEP


def login() -> SparclClient:
    user = os.getenv("SPARCL_USER")
    pwd = os.getenv("SPARCL_PASS")
    if not user or not pwd:
        raise RuntimeError("Set SPARCL_USER and SPARCL_PASS in the environment before running this.")

    client = SparclClient()
    client.login(user, pwd)
    return client


def batched_find(client, targetids, outfields, batch_size, max_retries=DEFAULT_MAX_RETRIES,
                  retry_sleep=DEFAULT_RETRY_SLEEP) -> pd.DataFrame:
    """Query client.find() in batches, retrying a failed batch before giving up on it.

    A handful of flaky batches shouldn't abort a query over the whole
    catalog, so failures are logged and skipped rather than raised -- the
    caller sees whatever came back from the batches that did succeed.
    """
    results = []
    for start in range(0, len(targetids), batch_size):
        subset = targetids[start:start + batch_size]
        for attempt in range(1, max_retries + 1):
            try:
                res = client.find(constraints={"targetid": subset}, outfields=outfields)
                results.append(pd.DataFrame(res.records))
                break
            except Exception as exc:
                print(f"find() batch {start}, attempt {attempt}/{max_retries}: {exc}")
                if attempt < max_retries:
                    time.sleep(retry_sleep)
                else:
                    print(f"  giving up on batch {start} after {max_retries} attempts")

    if not results:
        raise RuntimeError("SPARCL returned no results for any batch.")
    return pd.concat(results, ignore_index=True)
