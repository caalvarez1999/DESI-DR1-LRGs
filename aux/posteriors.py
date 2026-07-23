# -*- coding: utf-8 -*-
"""
CARLOS ALONSO ÁLVAREZ
Extract data from posteriors
"""

# Vendored from the author's personal library (MYLIBS/posteriors.py),
# trimmed to just allmodalgkde() -- what 9_FULLTMJ.py/
# pipeline/21_posteriorsTMJ.py use to turn a burned-in chain into a
# mode+16/84th-percentile summary per parameter. Left out: all(), the older
# beam-smoothed histogram version of the same idea; allmedianps(),
# all_COSMOGRAPHIC(), all_dARCHAEOLOGIC(), marginal(), peak_and_HMx() --
# none used by the TMJ pipeline stage. The module-level
# `from getdist import loadMCSamples, plots` in the original file is also
# dropped: unused by allmodalgkde().

import numpy as np

def allmodalgkde(fullfile, bmodel, prom_rel=0.05, min_sep_rel=0.02, return_peaks=True):
    from scipy.stats import gaussian_kde
    from scipy.signal import find_peaks, peak_prominences

    """
    Calcula modo (KDE), percentiles y diagnóstico de multimodalidad.

    Parámetros
    ----------
    fullfile : str
        Prefijo del archivo de cadena (se asume fullfile+'_burnin.txt').
    bmodel : array
        Archivo .npy con la malla / modelo (antes tenías bmodel; aquí se usa este argumento).
    prom_rel : float
        Prominencia mínima relativa respecto al máximo de densidad (0.05 = 5%).
    min_sep_rel : float
        Separación mínima relativa en x para considerar picos distintos, como fracción de (xmax-xmin).
    return_peaks : bool
        Si True, devuelve también detalles de picos (posiciones, prominencias, alturas).

    Returns
    -------
    out : (npar, 4) array
        Columnas: [mode, p16, med, p84]
    extra : dict (opcional)
        Por parámetro: número de picos, pico principal/secundario, razón de alturas, etc.
    """
    chain = np.loadtxt(fullfile + '_burnin.txt')
    fullstring = chain[:, 2:]  # igual que tú
    npar = fullstring.shape[1]

    out = np.zeros((npar, 4))
    extra = {
        "n_peaks": np.zeros(npar, dtype=int),
        "is_multimodal": np.zeros(npar, dtype=bool),
        "main_peak_x": np.full(npar, np.nan),
        "second_peak_x": np.full(npar, np.nan),
        "height_ratio_2nd_to_1st": np.full(npar, np.nan),
        "prom_ratio_2nd_to_1st": np.full(npar, np.nan),
    }

    # opcional: almacenar arrays variables (listas)
    peak_details = []  # cada elemento: dict con 'x', 'height', 'prom'

    for i in range(npar):
        sample = fullstring[:, i]

        kde = gaussian_kde(sample, bw_method="scott")

        # grilla: usa la de tu modelo pero fuerza que sea uniforme para find_peaks
        grid_raw = np.unique(bmodel[:, i])
        xmin, xmax = grid_raw.min(), grid_raw.max()

        # densidad evaluada en grilla uniforme (más estable para picos)
        # tamaño: al menos 512 puntos o igual al tamaño de tu grid_raw si es grande
        ngrid = max(512, grid_raw.size)
        grid = np.linspace(xmin, xmax, ngrid)

        dens = kde(grid)
        imax = np.argmax(dens)
        mode = grid[imax]

        p16 = np.percentile(sample, 16)
        med = np.median(sample)
        p84 = np.percentile(sample, 84)

        out[i, :] = [mode, p16, med, p84]

        # --- detección de picos ---
        prom_abs = prom_rel * dens.max()

        # separación mínima en puntos: fracción del rango convertida a índice
        min_sep_x = min_sep_rel * (xmax - xmin)
        dx = grid[1] - grid[0]
        min_sep_pts = max(1, int(np.round(min_sep_x / dx)))

        peaks, props = find_peaks(dens, prominence=prom_abs, distance=min_sep_pts)

        n_peaks = len(peaks)
        extra["n_peaks"][i] = n_peaks
        extra["is_multimodal"][i] = (n_peaks >= 2)

        if n_peaks >= 1:
            # ordena picos por altura
            heights = dens[peaks]
            order = np.argsort(heights)[::-1]
            peaks_sorted = peaks[order]
            heights_sorted = heights[order]

            # prominencias (ya las da find_peaks en props si usas prominence)
            proms_sorted = props["prominences"][order] if "prominences" in props else peak_prominences(dens, peaks_sorted)[0]

            main_idx = peaks_sorted[0]
            extra["main_peak_x"][i] = grid[main_idx]

            if n_peaks >= 2:
                print('error multimodal')
                second_idx = peaks_sorted[1]
                extra["second_peak_x"][i] = grid[second_idx]
                extra["height_ratio_2nd_to_1st"][i] = heights_sorted[1] / heights_sorted[0]
                extra["prom_ratio_2nd_to_1st"][i] = proms_sorted[1] / proms_sorted[0]

            if return_peaks:
                peak_details.append({
                    "x": grid[peaks_sorted],
                    "height": heights_sorted,
                    "prom": proms_sorted
                })
        else:
            if return_peaks:
                peak_details.append({"x": np.array([]), "height": np.array([]), "prom": np.array([])})

    if return_peaks:
        extra["peak_details"] = peak_details

    return out, extra
