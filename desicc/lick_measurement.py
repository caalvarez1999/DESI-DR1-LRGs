"""Lick-index measurement for a single galaxy spectrum.

This wraps the author's aux/pylick toolchain (air-wavelength conversion,
resolution matching, Monte Carlo index sampling), vendored in this repo
under aux/ and pylick/ -- see README.md for what's there and why.
"""

import numpy as np

from aux.corrections import C_toair, C_resol
from aux.emcee_lick import I_numberopen
from aux.indices_singlegalaxies import measure_indices

from .config import LICK_INDICES, RESTFRAME_WINDOW


def measure(w, f, df, s, mask, *, zgal: float, sn: float, nsampling: int,
            apply_resolution: bool, sigma_resolution: float) -> tuple[int, dict]:
    """Return (qual, {index: value, dindex: uncertainty}) for one spectrum.

    qual=0 means either the spectrum failed the S/N-driven quality cut, or
    the resolution-matching step raised on it; every index is left as NaN
    in that case.
    """
    values = {name: np.nan for name in LICK_INDICES}
    values.update({f"d{name}": np.nan for name in LICK_INDICES})

    good = mask == 0
    w, f, df, s = w[good], f[good], df[good], s[good]

    w_rest = C_toair(w) / (1.0 + zgal)
    lo, hi = RESTFRAME_WINDOW
    cut = (w_rest >= lo) & (w_rest <= hi)
    w_cut, f_cut, df_cut, s_cut = w_rest[cut], f[cut], df[cut], s[cut]

    qual = 1
    if np.sum(np.abs(df_cut)) < 0.1 or len(w_cut) < 0.05 * len(good) or len(w) < 0.1 * len(good):
        qual = 0

    if apply_resolution and qual:
        try:
            w_cut, f_cut, df_cut = C_resol(w_cut, 0.0, f_cut, df_cut, s_cut, sigma_resolution)
        except Exception:
            qual = 0
            w_cut = f_cut = df_cut = np.array([])

    if qual and len(w_cut) > 0:
        index_ids = I_numberopen(LICK_INDICES)
        try:
            vals, dvals, _, _ = measure_indices(
                w_cut, f_cut, df_cut,
                SN=float(sn), truesampling=True, nsampling=nsampling, index_list=index_ids,
            )
        except Exception:
            vals = np.full(len(LICK_INDICES), np.nan)
            dvals = np.full(len(LICK_INDICES), np.nan)

        for j, name in enumerate(LICK_INDICES):
            values[name] = float(vals[j]) if np.isfinite(vals[j]) else np.nan
            values[f"d{name}"] = float(dvals[j]) if np.isfinite(dvals[j]) else np.nan

    return qual, values
