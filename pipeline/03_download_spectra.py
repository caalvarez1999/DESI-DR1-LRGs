#!/usr/bin/env python3
"""Stage 03 (CLI) - download spectra for one velocity-dispersion bin,
compute snmedian, measure Lick indices for each galaxy, and fold the result
into PARENT3.fits.

Per-galaxy FITS files are written under --spectra-root as before; PARENT3.fits
is seeded from PARENT2.fits the first time it doesn't exist yet, then updated
in place across runs -- rows already marked done=1 are left untouched, so
re-running this on the same bin (or a bin you've already processed) only
fills in what's new.

Example:
    python3 pipeline/03_download_spectra.py --vd-bin 200225 --existing review

See pipeline/03_download_spectra_input.py for the same thing driven by
input() prompts instead of flags.
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from desicc.config import PARENT_DOWNLOADABLE, PARENT_FINAL, SINGLEGALAXIES_ROOT
from desicc.download import EXISTING_POLICIES, run


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vd-bin", required=True, help='e.g. "200225" for 200 <= vd < 225')
    parser.add_argument(
        "--existing", choices=EXISTING_POLICIES, default="review",
        help="review: redo only what's missing or invalid on disk (default) | "
             "replace: redownload everything in the bin | "
             "skip: only fetch files that don't exist yet, no content check",
    )
    parser.add_argument("--parent", type=Path, default=PARENT_DOWNLOADABLE)
    parser.add_argument("--spectra-root", type=Path, default=SINGLEGALAXIES_ROOT)
    parser.add_argument("--parent-final", type=Path, default=PARENT_FINAL,
                         help="PARENT3 table to fold this bin's measurements into (default: %(default)s)")
    args = parser.parse_args()

    stats = run(
        args.vd_bin, args.existing,
        parent_path=args.parent, spectra_root=args.spectra_root,
        parent_final_path=args.parent_final,
    )

    print("\n===== FOLD INTO PARENT3 =====")
    print(f"Updated             : {stats['updated']}")
    print(f"Skipped (done==1)   : {stats['skipped_done']}")
    print(f"Skipped (missing)   : {stats['skipped_missing']}")
    print(f"Skipped (bad FITS)  : {stats['skipped_bad']}")
    print("==============================")


if __name__ == "__main__":
    main()
