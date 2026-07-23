# spectra/

Every spectrum in this repo comes from DESI, so there's no separate `DESI/`
level -- two subfolders, both checked into the repo as an empty skeleton
(`.gitkeep` files) so the layout is visible without having to run anything --
the actual FITS files are **not** included (see the "Data & storage" section
of the top-level README).

## singlegalaxies/

Where stage 03 (`pipeline/03_download_spectra.py`) writes one FITS file per
galaxy: `spectra/singlegalaxies/vd<bin>/<targetid>.fits`.

The `vd<bin>/` folders cover the working sample's velocity-dispersion range,
`VD_MIN < vd < VD_MAX` (200-400, `desicc/config.py`, applied in stage 00) in
the same bins used in the ARTICLE3 analysis this pipeline reproduces:
`vd200225`, `vd225250`, `vd250280`, `vd280320`, `vd320355`, `vd355400`.

To populate a bin, run:

```bash
python3 pipeline/03_download_spectra.py --vd-bin 200225 --existing review
```

## stacks/

Where stage 10 (`pipeline/10_stackmaker.py`) writes one FITS file per stack:
`spectra/stacks/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/SN<target>_stack<NNNN>.fits`,
plus a `SN<target>_resume.json` checkpoint while a run is in progress.

The `vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/` folders cover every combination of
the same 6 vd bins above, the 3 archaeologies (`ALVAREZ`, `JOHANSSON`,
`THOMAS`), and the 2 cosmologies (`PLANCK`, `RIESS`) -- see
`ARCHAEOLOGY_CHOICES`/`COSMOLOGY_CHOICES` in `desicc/config.py`.

To populate one, run:

```bash
python3 pipeline/10_stackmaker.py --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 100
```

Stage 11 then collects a group's stack FITS files into one summary table,
and stage 12 CVD-corrects it -- both written to
`TABLAS/STACKed/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/` (see the top-level
README).
