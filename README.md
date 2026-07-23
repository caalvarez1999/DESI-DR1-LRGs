# DESI Cosmic Chronometers pipeline

Builds the working galaxy sample used in the Cosmic Chronometers analysis
(ARTICLE3): starting from the DESI fastspecfit ("fastspec-iron") value-added
catalog, select LRGs with a well-measured velocity dispersion in
`VD_MIN < vd < VD_MAX` (200-400 km/s, `desicc/config.py`), fetch their
spectra from [SPARCL](https://astrosparcl.datalab.noirlab.edu/), and measure
Lick/D4000 spectral indices on each one.

This covers the catalog-to-spectra stages (00-04) -- building the sample,
fetching spectra, measuring Lick indices per galaxy and folding those into a
single table, and applying the velocity-dispersion correction -- the
stacking stages (10-12): building S/N-targeted composite spectra per
vd/archaeology/cosmology group, summarizing them, and CVD-correcting the
summary -- the TMJ stellar-population fit stages (20-21): MCMC-fitting the
Thomas, Maraston & Johansson (2011) SSP model to each stack's indices to get
an age/metallicity/[alpha/Fe] posterior -- and the cosmographic fit (30):
MCMC-fitting a Taylor expansion of t(z) across several vd groups' ages to
get H(z0)/q(z0)/j(z0), the actual cosmochronometric measurement this whole
pipeline is built towards.

## Pipeline

Each stage is named after the `PARENTX.fits` it produces.

| stage | script | does |
|---|---|---|
| 00 | `pipeline/00_merge_catalogs.py` | Apply the sample selection to the raw catalog, write `PARENT0.fits` |
| 01 | `pipeline/01_query_sparcl_ids.py` | Look up each target's `specid`/`sparclid`/`DR` on SPARCL, write `PARENT1.fits` |
| 02 | `pipeline/02_select_downloadable.py` | Drop targets SPARCL couldn't identify, write `PARENT2.fits` |
| 03 | `pipeline/03_download_spectra.py` | Download spectra for one velocity-dispersion bin, measure Lick indices, write one FITS per galaxy, and fold the result into `PARENT3.fits` |
| 04 | `pipeline/04_cvd_correction.py` | Apply the velocity-dispersion correction to `PARENT3.fits`, write `PARENT4.fits` |
| 10 | `pipeline/10_stackmaker.py` | For one vd/archaeology/cosmology group, select eligible galaxies from `PARENT4.fits` and build S/N-targeted stacks, one FITS per stack |
| 11 | `pipeline/11_stacksummary.py` | Collect that group's stack FITS files into one summary table (`TABLAS/STACKed/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/stackedDESI_....fits`) |
| 12 | `pipeline/12_stackcvd.py` | Apply the velocity-dispersion correction to that summary table |
| 20 | `pipeline/20_SPSfitTMJ.py` | MCMC-fit the TMJ SSP model to every fittable stack, write a burned-in chain per stack |
| 21 | `pipeline/21_posteriorsTMJ.py` | Collect that group's chains into one posteriors table (`TABLAS/modelled/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/posteriorsDESI_TMJ_....fits`) |
| 30 | `pipeline/30_cosmographic.py` | MCMC-fit the Taylor-expansion t(z) model jointly across several vd groups' posteriors tables, write a burned-in chain |

Run them in order from the repo root:

```bash
python3 pipeline/00_merge_catalogs.py
python3 pipeline/01_query_sparcl_ids.py
python3 pipeline/02_select_downloadable.py
python3 pipeline/03_download_spectra.py --vd-bin 200225 --existing review
python3 pipeline/04_cvd_correction.py
python3 pipeline/10_stackmaker.py --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 100
python3 pipeline/11_stacksummary.py --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 100
python3 pipeline/12_stackcvd.py --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 100
python3 pipeline/20_SPSfitTMJ.py --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 100
python3 pipeline/21_posteriorsTMJ.py --vd-bin 200225 --archaeology ALVAREZ --cosmology PLANCK --target-sn 100
python3 pipeline/30_cosmographic.py --vd-bins 200225 225250 250280 280320 320355 --archaeology ALVAREZ --cosmology PLANCK --target-sn 150 --priors wide
```

Stage 03 also has an interactive twin, `03_download_spectra_input.py`, which
asks for the vd bin and the existing-file policy via prompts instead of
flags -- same underlying logic (`desicc/download.py`), just a different
front door.

`--existing` (or the equivalent prompt) controls what happens to spectra
already on disk for the requested bin:

- `review` (default): re-download only files that are missing or fail a
  structural check (`desicc.fits_io.is_valid_galaxy_fits`) -- the normal
  "top up the sample" mode.
- `replace`: re-download everything in the bin regardless of what's there.
- `skip`: only fetch files that don't exist at all, no content check.

Stage 03 is meant to be re-run once per vd bin: after downloading and
measuring, it folds that bin's per-galaxy FITS files into `PARENT3.fits`
(seeded from `PARENT2.fits` the first time it doesn't exist yet, then
updated in place across runs). Rows already marked `done=1` in `PYLICK` are
left untouched, so re-running a bin only fills in what's new.

Stage 04 is a full, stateless recompute: it always regenerates
`PARENT4.fits` from whatever is currently in `PARENT3.fits`, so it's
safe to re-run any time, e.g. after stage 03 has filled in more bins.

### Stacking (10-12)

A "group" is a choice of vd bin, `--archaeology` (the scaling relation
between velocity dispersion and stellar age: `ALVAREZ`/`JOHANSSON`/`THOMAS`),
`--cosmology` (age -> redshift: `PLANCK`/`RIESS`), and `--target-sn` (the
S/N each stack should reach). All three stages take the same four flags, and
operate on one group at a time.

Stage 10 first works out, from `archaeology/limitesarchaeologyTMJ_<COSMOLOGY>
_<ARCHAEOLOGY>.npy`, the maximum redshift a galaxy in this vd bin can have
and still plausibly be old enough to matter for this archaeology+cosmology
(`desicc/stack_selection.load_zcut`) -- i.e. how far back this vd group lets
you look. It then selects galaxies from `PARENT4.fits` within that redshift
and vd range, with a finished good-quality measurement (`done=1`,
`qual>=0.5`), no significant nebular emission, and a CaII H/K ratio under 1.2
(`desicc.stack_selection.select_for_stacking` -- ported from
`6_FINALSELECTION_of_singlegalaxies_by_groups.py`), sorts them by redshift,
and consumes them into consecutive stacks that each reach `--target-sn`
(`desicc/stacking.py` -- ported from `7_STACKINGfull.py`). This is a
long-running process -- it checkpoints to a `SN<target>_resume.json` file in
the output directory after every stack, and picks the checkpoint back up if
interrupted and re-run.

Stage 11 then collects that group's stack FITS files (`spectra/stacks/
vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/SN<target>_stack<NNNN>.fits`) into a single
summary table, `TABLAS/STACKed/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/
stackedDESI_vd<bin>_SN<sn>_<ARCHAEOLOGY>_<COSMOLOGY>.fits` (`<sn>` is
`--target-sn` zero-padded to 3 digits, e.g. `100` -> `100`, `50` -> `050`) --
the stacked-spectra equivalent of `PARENT3.fits`: one row per stack in
METADATA (`stack`, `ngals`, `vd`/`dvd`, `z1`/`dz1`, `snmedian` -- properties
of the stack itself, not present in any individual-galaxy file) and PYLICK
(`stack`, replicated from METADATA as the row identifier -- same convention
as `targetid` in the per-galaxy PYLICK tables -- plus the indices measured
on that stack's composite spectrum). Stage 12 applies the velocity-dispersion
correction to that summary table, the same way stage 04 does for
`PARENT3.fits`, using each stack's own mean vd, and writes the
`..._CVD.fits` sibling in the same directory -- a full, stateless recompute,
safe to re-run any time stage 11 has picked up more stacks.

### TMJ fit (20-21)

Stages 20-21 take the same `--vd-bin`/`--archaeology`/`--cosmology`/
`--target-sn` group flags as stages 10-12, reading stage 12's CVD-corrected
summary table, plus `--grid` (default `custom250924b`, the only TMJ model
grid vendored in `aux/TMJ/` -- see Vendored code below).

Stage 20 selects the stacks with >=1 galaxy and finite values/errors on
every index in `TMJ_FIT_INDICES` (`desicc/config.py` -- a fixed set of 13
Lick indices from the original analysis, unrelated to the `--archaeology`
scaling-relation choice despite a similar name upstream), and for each one
MCMC-fits age (`t`)/metallicity (`Z`)/`[alpha/Fe]` (`afe`) against the model
grid (`desicc/tmj.py`, `aux/mcmc_tmj.py` -- ported from
`STACKING/MCMC2_MCMC.py`), then cuts the burn-in
(`aux/burnin.py` -- ported from `STACKING/MCMC3_burnin_maker.py`). Each
stack's burned-in chain is written to `chains/TMJ/vd<bin>/<ARCHAEOLOGY>/
<COSMOLOGY>/chainDESIFULLTMJ<grid>_stack<NNNN>_vd<bin>_SN<sn>_<ARCHAEOLOGY>
_<COSMOLOGY>_burnin.txt` (`<NNNN>` is the stack's own number, read back from
its `stack` column, e.g. `..._stack0001.fits` -> `0001`); the much larger
raw `.h5` chain is deleted right after. This is a long-running process (one
MCMC run per stack) -- a stack whose `_burnin.txt` already exists is skipped,
so it's safe to re-run any time or after being interrupted.

Stage 21 then collects that group's chains into a single posteriors table,
`TABLAS/modelled/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/posteriorsDESI_TMJ_<grid>
_vd<bin>_SN<sn>_<ARCHAEOLOGY>_<COSMOLOGY>.fits` -- the same METADATA+PYLICK
structure as stage 11/12's table, plus a TMJ HDU with one row per stack:
`stack`, then for each of `t`/`Z`/`afe` the KDE mode plus 16th/84th-percentile
spread (`<par>`, `d<par>_low`, `d<par>_up`, `d<par>`, `<par>_npeaks` --
`desicc/tmj.py`, `aux/posteriors.py` -- ported from `MYLIBS/posteriors.py`),
and `valid` (1 if the stack was fittable and its chain exists, 0 -- with the
rest of the row NaN -- otherwise). A full, stateless recompute, safe to
re-run any time stage 20 has produced more chains.

### Cosmographic fit (30)

Stage 30 takes `--vd-bins` (one or more groups, e.g. `200225 225250 250280
280320 320355`), `--archaeology`/`--cosmology`/`--target-sn`/`--grid` (same
meaning as stages 20-21, used to find each group's stage 21 table), and
`--priors` (`narrow` or `wide` -- the prior ranges on `Hz0`/`qz0`/`jz0`,
`aux/taylor_tz.py`/`aux/taylor_tz_wide.py`).

It reads and stacks the selected groups' TMJ posteriors tables, keeps only
reliable rows (`valid`, `z1/dz1 > 50`, `dz1 < 0.01`, `t/dt > 5`, `dt < 1` --
`desicc.cosmographic.select_fittable`, ported from
`13_tz_cosmographic.py`'s `RESRES`), and MCMC-fits a `--order`-th order
Taylor expansion of t(z) around `--z0` (default: the sample's median z1)
jointly across every selected group -- one shared `Hz0`/`qz0`/`jz0` plus one
`age<vd-bin>` nuisance parameter per group (`desicc/cosmographic.py`,
`aux/mcmc_cosmographic.py` -- ported from `13_tz_cosmographic.py` +
`aux13_MCMC.py`), then cuts the burn-in (`aux/burnin_cosmographic.py` --
ported from `aux13_burnin_maker.py`). The burned-in chain is written to
`chains/cosmographic/<name>_burnin.txt` (`<name>` defaults to
`cosmographicTMJ<grid>_vd<bins->joined->with->dashes>_SN<sn>_<ARCHAEOLOGY>
_<COSMOLOGY>_<priors>`, overridable with `--name`); the raw `.h5` chain and
the un-burned `.txt` are deleted right after. A chain whose `_burnin.txt`
already exists is skipped unless `--overwrite` is given.

## Layout

```
desicc/                    shared library code (config, SPARCL client, FITS I/O, Lick measurement, CVD correction, stacking)
pipeline/                  the numbered entry-point scripts above
aux/                       vendored Lick-index/CVD primitives and data this pipeline actually uses (see Vendored code below)
pylick/                    vendored index-measurement package aux/ builds on
archaeology/               archaeology+cosmology -> redshift-cut lookups stage 10 reads -- included (see Data & storage)
tablasDESI/                where the raw fastspec-iron VAC would go -- not included (see Data & storage)
TABLAS/PARENT/             PARENTX.fits summary tables -- included (see Data & storage)
TABLAS/STACKed/             stage 11/12 stack summary tables, one subfolder per vd/archaeology/cosmology group -- included (see Data & storage)
TABLAS/modelled/            stage 21 TMJ posteriors tables, one subfolder per vd/archaeology/cosmology group -- included (see Data & storage)
spectra/singlegalaxies/    per-galaxy FITS files (stage 03), one subfolder per vd bin -- folder skeleton only, files not included
spectra/stacks/             per-stack FITS files (stage 10), one subfolder per vd/archaeology/cosmology group -- folder skeleton only, files not included
chains/TMJ/                 stage 20 MCMC chains, one subfolder per vd/archaeology/cosmology group -- folder skeleton only, files not included
chains/cosmographic/        stage 30 MCMC chains -- folder skeleton only, files not included
stackstosinglegalaxies.fits  stack -> targetids lookup for the baseline stacks -- included (see Data & storage)
```

Every spectrum here comes from DESI, so `spectra/` has no separate `DESI/`
level (see `spectra/README.md`).

## Vendored code

`aux/` and `pylick/` are real files in this repo, not symlinks to the
author's personal astronomy library -- the pipeline runs standalone, no
other repo needed (besides the raw catalog and whatever spectra/galaxies it
downloads itself).

Only the pieces this pipeline actually calls are vendored, copied verbatim
from the author's personal library (a 37 MB library covering many unrelated
analyses, referred to below as "the original library" -- not part of this
repo):

- `aux/corrections.py`, `aux/emcee_lick.py`, `aux/indices_singlegalaxies.py`
  -- `C_toair`, `C_resol`, `C_VD`, `I_numberopen`/`I_typeopen`,
  `measure_indices`, used by `desicc/lick_measurement.py`, `desicc/cvd.py`,
  and `desicc/stacking.py`. In the original library these live under a package
  called `MYLIBS`; only its internal absolute import
  (`import MYLIBS.emcee_lick`, now `import aux.emcee_lick`) was changed to
  match, everything else is untouched.
- `aux/CVD/MILES/{C,dC}funs_MILESstars.npy` -- the velocity-dispersion
  correction grids `desicc/cvd.py` looks up, precomputed on MILES stellar
  templates.
- `pylick/` (all of it) -- the index-measurement package
  `indices_singlegalaxies.py` calls into.
- `aux/tableall.dat`, `aux/plotbelli.style` -- the only two files pylick
  actually opens (from `pylick/_config.py`'s `dir_lib = './aux/'`, resolved
  relative to the process's current working directory, so `aux/` has to
  live wherever you run `python3 pipeline/...` from -- normally the repo
  root, same as here). In the original library these come from a separate
  `libs/` data folder; `pylick/_config.py` and the hardcoded path in
  `pylick/indices.py` were repointed at `aux/` to match.
- `aux/mcmc_tmj.py` (`MCMC_one`, from `STACKING/MCMC2_MCMC.py`),
  `aux/burnin.py` (`burnin`, from `STACKING/MCMC3_burnin_maker.py`),
  `aux/posteriors.py` (`allmodalgkde`, from `MYLIBS/posteriors.py`) -- the
  MCMC fit, burn-in cut, and posterior summary `desicc/tmj.py` and
  `pipeline/20_SPSfitTMJ.py`/`21_posteriorsTMJ.py` use. Each file only keeps
  the one function this pipeline calls (plus its direct dependencies) out of
  several fit variants the original files defined; see the header comment in
  each for exactly what was left out and why (mainly a `Z`/`[alpha/Fe]`-only
  fit variant unused here, and, in `MCMC2_MCMC.py`'s case, ~1.4 GB of
  module-level `np.load(...)` calls for a *different*, unused fit variant in
  the same file).
- `aux/TMJ/TMJ_MILES_custom250924b.npy`, `aux/TMJ/tabla/
  TMJ_MILES_custom250924b_tabla.npy` -- the TMJ SSP model grid (age x
  metallicity x `[alpha/Fe]` x Lick index) and its flattened lookup table,
  precomputed on MILES stellar templates. Only this one grid is vendored --
  it is the only one with full, consistent coverage in the historical
  `TABLAS/modelled/` results (see Data & storage); other grid tags exist in
  the original analysis but only cover one-off test combinations.
- `aux/mcmc_cosmographic.py` (`MCMC_one`, from `ARTICLE3/aux13_MCMC.py`),
  `aux/burnin_cosmographic.py` (`burnin`, from
  `ARTICLE3/aux13_burnin_maker.py`), `aux/taylor_tz.py`/`aux/taylor_tz_wide.py`
  (the narrow/wide-priors log-posteriors, from
  `MYLIBS/TAYLOR_emcee_tz.py`/`TAYLOR_emcee_tz_widepriors.py`) -- the
  cosmographic MCMC fit, burn-in cut, and Taylor-expansion model
  `desicc/cosmographic.py`/`pipeline/30_cosmographic.py` use. The wide
  variant isn't just wider priors -- its `C3` Taylor coefficient also
  differs from the narrow variant's, an inconsistency present in the
  original library and kept as-is (see the header comment in
  `aux/taylor_tz_wide.py`).

If you update from the author's original library, re-copy these same files
into `aux/`/`pylick/` rather than re-symlinking the whole thing, to keep the
repo self-contained.

## Setup

Python 3.11+, plus:

```bash
pip install -r requirements.txt
```

SPARCL credentials, for stages 01 and 03:

```bash
export SPARCL_USER=...
export SPARCL_PASS=...
```

`aux/` and `pylick/` (see Vendored code above) ship with the repo, so no
further setup is needed for those -- just run everything from the repo
root, since `aux/` is read relative to the current working directory.

## Data & storage

This repo is meant to be a blueprint of the pipeline, not a data dump --
what's tracked here is either code, or tables small and useful enough to
read on their own. Deliberately *not* included:

- **`tablasDESI/`** -- `RAW_CATALOG` in `desicc/config.py` points at the raw
  fastspec-iron VAC (`tablasDESI/v2.1/fastspec-iron-main-dark.fits`), stage
  00's input. It isn't distributed here: it's DESI's own data product,
  already published in the official DESI data release, so duplicating it in
  this repo would just be redundant. Only the empty `v2.1/`/`v3.0/`
  directory skeleton is kept, to show where it goes.
- **`TABLAS/PARENT/PARENT{0,1,2,3,4}.fits`** -- the actual output of every
  catalog-to-spectra stage, real and complete (see below) -- but at
  700MB-1.2GB each they exceed GitHub's hard 100MB-per-file push limit, so
  they stay local-only rather than pulling in Git LFS for just these five
  files. Regenerate them by running stages 00-04 (needs `tablasDESI/` and
  SPARCL/spectra access, see Setup below).
- **`spectra/singlegalaxies/vd<bin>/<targetid>.fits`** -- the
  per-galaxy spectra and Lick measurements stage 03 writes. There's one file
  per galaxy (order of 10^5-10^6 across all bins), which is both far too
  much data and far too many files to track in git.
- **`spectra/stacks/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/*.fits`** -- the
  per-stack composite spectra stage 10 writes. Same reasoning: regenerable,
  and there can be many per group.
- **`chains/TMJ/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/*_burnin.txt`** -- stage
  20's burned-in MCMC chains, one per stack. Same reasoning: regenerable,
  and potentially large in aggregate (one full MCMC run each).
- **`chains/cosmographic/*_burnin.txt`** -- stage 30's burned-in MCMC
  chains. Same reasoning.
- **`sparcl_ids_cache.fits`** -- stage 01's SPARCL query cache. Purely a
  speed-up for re-runs, not informative on its own.

For `spectra/`, `chains/TMJ/` and `chains/cosmographic/`, only the directory
skeleton is kept (see `spectra/README.md`), so it's clear where each stage
puts its output; the files themselves are regenerated by running the
pipeline.

What *is* included:

- Everything else in `TABLAS/`, i.e. `TABLAS/STACKed/` and
  `TABLAS/modelled/` below -- the useful checkpoints small enough to check
  in. Each one is a real, complete summary of the sample (or the stacks, or
  the fits) at that point in the pipeline, readable directly (with
  `astropy.table.Table.read`, TOPCAT, etc.) without re-running anything, and
  together with the (locally-kept) PARENT tables they're the actual
  scientific content this repo produces -- worth keeping even though the
  repo otherwise avoids checking in data.
- **`TABLAS/STACKed/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/stackedDESI_....fits`**
  -- stage 11/12 stack summary tables, same idea as the PARENT tables but
  for the stacks. This includes the full historical set of stacks from the
  original ARTICLE3 analysis (all 6 in-range vd bins x 3 archaeologies x 2
  cosmologies x 11 target-S/N values, raw and CVD-corrected), copied and
  converted from ARTICLE3's `.ecsv` tables into this repo's METADATA+PYLICK
  `.fits` convention -- not (re-)produced by running stages 10-12 in this
  repo, since that would mean re-downloading and re-stacking the full sample.
- **`TABLAS/modelled/vd<bin>/<ARCHAEOLOGY>/<COSMOLOGY>/posteriorsDESI_TMJ_....fits`**
  -- stage 21 TMJ posteriors tables, same idea again but for the TMJ fit.
  Includes the full historical `custom250924b` set (all 6 in-range vd bins x
  3 archaeologies x 2 cosmologies x 11 target-S/N values), copied and
  converted from ARTICLE3's `.ecsv` tables the same way as `TABLAS/STACKed/`
  -- not (re-)produced by running stages 20-21 in this repo, since that
  would mean re-running ~400 individual MCMC fits.
- **`archaeology/limitesarchaeologyTMJ_<COSMOLOGY>_<ARCHAEOLOGY>.npy`** --
  tiny (~1 KB each) precomputed lookups stage 10 needs to even select a
  sample, not something regenerated by any stage in this repo.
- **`stackstosinglegalaxies.fits`** -- a single `STACKS` HDU with `stack`
  (`vd<bin>_stack<NNNN>`), `ngals`, and `targetids` (the variable-length
  array of single-galaxy `targetid`s that went into that stack), for the
  394 "baseline" stacks (SN150, PLANCK, ALVAREZ, the same configuration
  used in the README's cosmographic cross-check) across the 5 vd groups
  200225-320355. Built once from the real per-stack FITS files this
  pipeline's stage 10 originally produced for that run
  (`WRK/spectra/DESI/STACKS/onlyfinal/vd<bin>/`, outside this repo -- each
  one's own METADATA already has a `targetids` column, a JSON-encoded list,
  see `desicc/stacking.py`). Individual stack FITS files aren't tracked in
  git (see `spectra/stacks/` above), and stage 11's summary table doesn't
  carry `targetids` over from them, so this is the only place in the repo a
  stack name can be traced back to the actual DESI targets it was built from.

Every stage also takes `--input`/`--output` (or, for stage 03,
`--parent`/`--spectra-root`/`--parent-final`; for stage 04, `--input`/
`--output`; for stages 10-12 and 20-21, `--singlegalaxies-root`/
`--stacks-root`/`--chains-root`/`--output`; for stage 30, `--chains-root`/
`--name`, on top of `--vd-bin`/`--archaeology`/`--cosmology`/`--target-sn`,
or `--vd-bins` for stage 30) flags if you'd rather keep the data somewhere
else. Actually running the pipeline (downloading the raw catalog, querying
SPARCL, fetching spectra) needs external services -- that part can't be
made to run standalone -- but no other repo or library is required.
