#!/usr/bin/env python3
"""Stage 03 (interactive) - same as 03_download_spectra.py, but asks for the
vd bin and the existing-file policy via input() instead of CLI flags.

Downloads and measures spectra for the bin, then folds the result into
PARENT3.fits (seeded from PARENT2.fits the first time it doesn't exist yet).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from desicc.download import run


def ask_existing_policy() -> str:
    review = input("Check FITS files already on disk before (re)downloading? [Y/n] ").strip().lower()
    if review in ("", "y", "yes"):
        return "review"

    replace = input("Redownload everything in this bin, replacing what's already there? [y/N] ").strip().lower()
    return "replace" if replace in ("y", "yes") else "skip"


def main():
    vd_bin = input('Velocity dispersion bin, e.g. "200225" for 200 <= vd < 225: ').strip()
    policy = ask_existing_policy()
    stats = run(vd_bin, policy)

    print("\n===== FOLD INTO PARENT3 =====")
    print(f"Updated             : {stats['updated']}")
    print(f"Skipped (done==1)   : {stats['skipped_done']}")
    print(f"Skipped (missing)   : {stats['skipped_missing']}")
    print(f"Skipped (bad FITS)  : {stats['skipped_bad']}")
    print("==============================")


if __name__ == "__main__":
    main()
