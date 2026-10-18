# DESI-DR1 Cosmic Chronometers pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21535100.svg)](https://doi.org/10.5281/zenodo.21535100)

This repository contains the pipeline used in *Cosmography with DESI-DR1 Cosmic
Chronometers: Direct H(z) measurements from Luminous Red Galaxy ages*
(Álvarez et al. 2026). Starting from the DESI DR1 FastSpecFit value-added
catalogue (`fastspec-iron`), it:

1. selects Luminous Red Galaxies (LRGs) with a reliable stellar velocity
   dispersion in the range 200 < σ < 400 km/s;
2. retrieves their spectra from [SPARCL](https://astrosparcl.datalab.noirlab.edu/)
   and measures Lick/D4000 spectral indices on each of them;
3. builds stacked spectra at a target S/N for each velocity-dispersion group;
4. fits the Thomas, Maraston & Johansson (2011; TMJ) simple stellar population
   model to the indices of each stack to obtain its age, metallicity and
   [α/Fe];
5. fits a pivotal-redshift cosmographic expansion of the age–redshift relation,
   t(z), jointly to all the velocity-dispersion groups to constrain
   H(z₀), q(z₀) and j(z₀).

The large catalogue-level tables (`PARENT0`–`PARENT4.fits`) are available on
Zenodo: [doi:10.5281/zenodo.21535100](https://doi.org/10.5281/zenodo.21535100).

## Citation

If you use this code or the associated data products, please cite:

> Álvarez, C. A., Cueli, M. M., Haridasu, B. S., et al. 2026, *Cosmography
> with DESI-DR1 Cosmic Chronometers: Direct H(z) measurements from Luminous
> Red Galaxy ages*, A&A (submitted)

and, for the `PARENT` tables, the Zenodo record
[doi:10.5281/zenodo.21535100](https://doi.org/10.5281/zenodo.21535100).

This work uses DESI DR1 data; please also follow the DESI data release
acknowledgement and citation guidelines if you use them.

## Pipeline

Each catalogue-level stage (00–04) is named after the `PARENTX.fits` table it
produces.

| stage | script | description |
|---|---|---|
| 00 | `pipeline/00_merge_catalogs.py` | Apply the sample selection to the FastSpecFit catalogue and write `PARENT0.fits` |
| 01 | `pipeline/01_query_sparcl_ids.py` | Retrieve the SPARCL identifiers (`specid`, `sparclid`, `DR`) of each target and write `PARENT1.fits` |
| 02 | `pipeline/02_select_downloadable.py` | Remove the targets not found in SPARCL and write `PARENT2.fits` |
| 03 | `pipeline/03_download_spectra.py` | Download the spectra of one velocity-dispersion bin, measure their Lick indices, write one FITS file per galaxy and add the results to `PARENT3.fits` |
| 04 | `pipeline/04_cvd_correction.py` | Apply the velocity-dispersion correction to `PARENT3.fits` and write `PARENT4.fits` |
| 10 | `pipeline/10_stackmaker.py` | Select the eligible galaxies of one group from `PARENT4.fits` and build stacks at the target S/N, one FITS file per stack |
| 11 | `pipeline/11_stacksummary.py` | Collect the stacks of one group into a summary table |
| 12 | `pipeline/12_stackcvd.py` | Apply the velocity-dispersion correction to that summary table |
| 20 | `pipeline/20_SPSfitTMJ.py` | Fit the TMJ model to every stack with an MCMC sampler and write one chain per stack |
| 21 | `pipeline/21_posteriorsTMJ.py` | Collect the chains of one group into a table of posterior summaries |
| 30 | `pipeline/30_cosmographic.py` | Fit the cosmographic t(z) model jointly to several groups and write the MCMC chain |

Run the stages in order from the repository root. For the baseline
configuration of the paper (S/N = 150, age–σ relation of Álvarez et al. 2025,
Planck cosmology):

```bash
python3 pipeline/00_merge_catalogs.py
python3 pipeline/01_query_sparcl_ids.py
python3 pipeline/02_select_downloadable.py
python3 pipeline/03_download_spectra.py --vd-bin 200225 --existing review
python3 pipeline/04_cvd_correction.py
python3 pipeline/10_stackmaker.py   --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 150
python3 pipeline/11_stacksummary.py --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 150
python3 pipeline/12_stackcvd.py     --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 150
python3 pipeline/20_SPSfitTMJ.py    --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 150
python3 pipeline/21_posteriorsTMJ.py --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 150
python3 pipeline/30_cosmographic.py --vd-bins 200225 225250 250280 280320 320355 --archaeology ALVAREZ --cosmology PLANCK --target-sn 150 --priors wide
```

Stages 03 and 10–21 are run once per velocity-dispersion bin
(`200225`, `225250`, `250280`, `280320`, `320355`, `355400`).

### Catalogue and spectra (00–04)

Stage 03 has an interactive version, `03_download_spectra_input.py`, which
asks for the velocity-dispersion bin and the policy for existing files instead
of reading them from flags; both use the same code (`desicc/download.py`).

The `--existing` option controls how spectra already on disk are treated:

- `review` (default): re-download only the files that are missing or fail a
  structural check (`desicc.fits_io.is_valid_galaxy_fits`);
- `replace`: re-download all the files of the bin;
- `skip`: download only the files that do not exist, without checking the
  existing ones.

After downloading and measuring the spectra of a bin, stage 03 adds the
results to `PARENT3.fits` (created from `PARENT2.fits` on the first run and
updated in place afterwards). Rows already marked as `done=1` in the `PYLICK`
table are left untouched, so the stage can be resumed or re-run bin by bin.

Stage 04 always regenerates `PARENT4.fits` from the current content of
`PARENT3.fits`, so it can be re-run at any time.

### Stacking (10–12)

A *group* is defined by a velocity-dispersion bin and three further choices:

- `--archaeology`: the age–velocity dispersion relation used for the
  redshift cut (`ALVAREZ`, `JOHANSSON` or `THOMAS`, after Álvarez et al. 2025,
  Johansson et al. 2012 and Thomas et al. 2010);
- `--cosmology`: the cosmology used to convert ages into redshifts (`PLANCK`
  or `RIESS`);
- `--target-sn`: the spectral S/N that each stack must reach.

Stages 10–12 take the same four flags and process one group at a time.

Stage 10 first reads the maximum redshift allowed for the bin, given the
chosen relation and cosmology, from
`archaeology/limitesarchaeologyTMJ_<COSMOLOGY>_<ARCHAEOLOGY>.npy`
(the *archaeological coherence cut* of the paper; `desicc/stack_selection.load_zcut`).
It then selects the galaxies of `PARENT4.fits` within that redshift and
velocity-dispersion range that have a valid measurement (`done=1`,
`qual>=0.5`), no significant nebular emission and a CaII H/K ratio below 1.2
(`desicc.stack_selection.select_for_stacking`). The selected galaxies are
sorted by redshift and combined into consecutive stacks that reach the target
S/N (`desicc/stacking.py`). This stage is long; it saves a checkpoint
(`SN<target>_resume.json`) after each stack and resumes from it if
interrupted.

Stage 11 collects the stack files of the group
(`spectra/stacks/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/SN<target>_stack<NNNN>.fits`)
into a single summary table,
`TABLAS/STACKed/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/stackedDESI_vd<bin>_SN<sn>_<ARCHAEOLOGY>_<COSMOLOGY>.fits`,
where `<sn>` is the target S/N zero-padded to three digits (e.g. `050`, `150`).
Its `METADATA` table contains one row per stack (`stack`, `ngals`, `vd`/`dvd`,
`z1`/`dz1`, `snmedian`), and its `PYLICK` table contains the indices measured
on each stacked spectrum, with `stack` as row identifier. Stage 12 applies the
velocity-dispersion correction to this table, using the mean velocity
dispersion of each stack, and writes the corresponding `_CVD.fits` file.

### TMJ fit (20–21)

Stages 20–21 take the same group flags as stages 10–12 and read the
corrected summary table produced by stage 12. The `--grid` option selects the
TMJ model grid (default and only grid included: `custom250924b`).

Stage 20 selects the stacks with at least one galaxy and finite values and
uncertainties for all the indices used in the fit (`TMJ_FIT_INDICES` in
`desicc/config.py`, the 13 Lick indices of the paper), and fits the age (`t`),
metallicity (`Z`) and [α/Fe] (`afe`) of each one with an MCMC sampler
(`desicc/tmj.py`, `aux/mcmc_tmj.py`). The burn-in is then removed
(`aux/burnin.py`), and the chain of each stack is written to
`chains/TMJ/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/chainDESIFULLTMJ<grid>_stack<NNNN>_vd<bin>_SN<sn>_<ARCHAEOLOGY>_<COSMOLOGY>_burnin.txt`.
Stacks whose chain already exists are skipped, so the stage can be resumed.

Stage 21 collects the chains of the group into a single table,
`TABLAS/modelled/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/posteriorsDESI_TMJ_<grid>_vd<bin>_SN<sn>_<ARCHAEOLOGY>_<COSMOLOGY>.fits`,
with the same `METADATA` and `PYLICK` tables as stages 11–12 plus a `TMJ`
table with one row per stack. For each of `t`, `Z` and `afe`, it gives the
mode of the posterior and its 16th–84th percentile range (`<par>`,
`d<par>_low`, `d<par>_up`, `d<par>`, `<par>_npeaks`; `desicc/tmj.py`,
`aux/posteriors.py`), and a `valid` flag (1 if the stack was fitted).

### Cosmographic fit (30)

Stage 30 takes `--vd-bins` (one or more groups, e.g.
`200225 225250 250280 280320 320355`), the `--archaeology`, `--cosmology`,
`--target-sn` and `--grid` options of stages 20–21 (used to locate the
corresponding tables), and `--priors` (`narrow` or `wide`), which sets the
prior ranges on `Hz0`, `qz0` and `jz0` (`aux/taylor_tz.py`,
`aux/taylor_tz_wide.py`).

It combines the posterior tables of the selected groups and keeps only the
stacks that pass the quality flags of the paper (`valid`, `z1/dz1 > 50`,
`dz1 < 0.01`, `t/dt > 5`, `dt < 1` Gyr; `desicc.cosmographic.select_fittable`).
It then fits an expansion of t(z) of order `--order` around `--z0` (by
default, the median redshift of the selected stacks) jointly to all the
groups, with common `Hz0`, `qz0` and `jz0` and one nuisance age parameter per
group (`desicc/cosmographic.py`, `aux/mcmc_cosmographic.py`). After removing
the burn-in (`aux/burnin_cosmographic.py`), the chain is written to
`chains/cosmographic/<name>_burnin.txt`, where `<name>` defaults to
`cosmographicTMJ<grid>_vd<bins>_SN<sn>_<ARCHAEOLOGY>_<COSMOLOGY>_<priors>` and
can be changed with `--name`. Existing chains are not overwritten unless
`--overwrite` is given.

All stages also accept options to read from and write to custom locations
(`--input`/`--output` for most stages; `--parent`, `--spectra-root` and
`--parent-final` for stage 03; `--singlegalaxies-root`, `--stacks-root`,
`--chains-root` and `--output` for stages 10–21; `--chains-root` and `--name`
for stage 30).

## Repository layout

```
desicc/                  library code (configuration, SPARCL client, FITS I/O, Lick-index measurement, velocity-dispersion correction, stacking, fits)
pipeline/                numbered scripts of each stage
aux/                     auxiliary modules and data used by the pipeline (see Included third-party code and data)
pylick/                  index-measurement package (see Included third-party code and data)
archaeology/             redshift limits of the archaeological coherence cut, read by stage 10
tablasDESI/              location of the FastSpecFit catalogue (not included)
TABLAS/PARENT/           location of the PARENT0-4 tables (available on Zenodo)
TABLAS/STACKed/          stack summary tables (stages 11-12), one folder per group
TABLAS/modelled/         TMJ posterior tables (stage 21), one folder per group
spectra/singlegalaxies/  per-galaxy FITS files (stage 03), one folder per velocity-dispersion bin (not included)
spectra/stacks/          per-stack FITS files (stage 10), one folder per group (not included)
chains/TMJ/              MCMC chains of stage 20, one folder per group (not included)
chains/cosmographic/     MCMC chains of stage 30 (not included)
stackstosinglegalaxies.fits  list of the DESI targets included in each baseline stack
```

## Included third-party code and data

The repository is self-contained: apart from the FastSpecFit catalogue and the
spectra downloaded from SPARCL, no external code or data are needed.

- `pylick/`: a copy of the pyLick package (Borghi et al. 2022) used to measure
  the spectral indices. The only modification is that its auxiliary files
  (`tableall.dat`, `plotbelli.style`) are read from `aux/`, relative to the
  working directory, so the scripts must be run from the repository root.
- `aux/corrections.py`, `aux/emcee_lick.py`, `aux/indices_singlegalaxies.py`:
  air-wavelength conversion, resolution matching, velocity-dispersion
  correction and index-measurement functions, used by
  `desicc/lick_measurement.py`, `desicc/cvd.py` and `desicc/stacking.py`.
- `aux/CVD/MILES/{C,dC}funs_MILESstars.npy`: velocity-dispersion correction
  functions and their uncertainties, computed from MILES stellar spectra
  (Appendix B of the paper).
- `aux/TMJ/TMJ_MILES_custom250924b.npy` and
  `aux/TMJ/tabla/TMJ_MILES_custom250924b_tabla.npy`: interpolated grid of TMJ
  model predictions (age × metallicity × [α/Fe] × Lick index) and its lookup
  table.
- `aux/mcmc_tmj.py`, `aux/burnin.py`, `aux/posteriors.py`: MCMC fit, burn-in
  removal and posterior summary of the TMJ fit.
- `aux/mcmc_cosmographic.py`, `aux/burnin_cosmographic.py`,
  `aux/taylor_tz.py`, `aux/taylor_tz_wide.py`: MCMC fit, burn-in removal and
  log-posteriors (narrow and wide priors) of the cosmographic fit.

## Setup

Python 3.11 or later and the packages in `requirements.txt`:

```bash
pip install -r requirements.txt
```

Stages 01 and 03 require SPARCL credentials:

```bash
export SPARCL_USER=...
export SPARCL_PASS=...
```

Run all the scripts from the repository root, since `aux/` is read relative to
the working directory.

## Data

### Not included in this repository

- **FastSpecFit catalogue** (`tablasDESI/v2.1/fastspec-iron-main-dark.fits`,
  the input of stage 00, set by `RAW_CATALOG` in `desicc/config.py`). It is a
  public DESI DR1 data product and can be downloaded from
  <https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v2.1/catalogs/>.
- **`PARENT0`–`PARENT4.fits`**, the outputs of stages 00–04 (0.7–1.2 GB each).
  They exceed the GitHub file-size limit and are available on Zenodo:
  [doi:10.5281/zenodo.21535100](https://doi.org/10.5281/zenodo.21535100).
  To use them with the pipeline, place them in `TABLAS/PARENT/`.
- **Per-galaxy spectra** (`spectra/singlegalaxies/`), per-stack spectra
  (`spectra/stacks/`) and **MCMC chains** (`chains/`). They are too large or
  too numerous to be tracked, and they can be regenerated with the pipeline;
  only the folder structure is kept (see `spectra/README.md`). The individual
  spectra are publicly available from DESI DR1 through SPARCL.
- **SPARCL query cache** (`sparcl_ids_cache.fits`), which only speeds up
  re-runs of stage 01.

### Included in this repository

- **`TABLAS/STACKed/`**: stack summary tables (stages 11–12) for all the
  configurations explored in the paper: 6 velocity-dispersion bins × 3
  age–σ relations × 2 cosmologies × 11 target S/N levels (50–150 in steps of
  10), before and after the velocity-dispersion correction.
- **`TABLAS/modelled/`**: TMJ posterior tables (stage 21) for the same
  configurations, obtained with the `custom250924b` grid.
- **`archaeology/limitesarchaeologyTMJ_<COSMOLOGY>_<ARCHAEOLOGY>.npy`**:
  redshift limits of the archaeological coherence cut for each combination
  of age–σ relation and cosmology, read by stage 10.
- **`stackstosinglegalaxies.fits`**: a `STACKS` table with the name of each
  stack (`vd<bin>_stack<NNNN>`), its number of galaxies (`ngals`) and the
  DESI `targetid`s of the galaxies that compose it, for the 394 stacks of the
  baseline configuration (S/N = 150, `ALVAREZ`, `PLANCK`) in the five groups
  used in the cosmographic fit (200–355 km/s). Together with the `PARENT`
  tables, it allows each stack to be traced back to its individual DESI
  spectra.

All the tables can be read directly, e.g. with `astropy.table.Table.read` or
TOPCAT.
