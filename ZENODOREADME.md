# ZENODO README FOR PARENTX.fits [X == 0, 1, 2, 3, 4]

Tables are built using the following stages of the pipeline freely available in https://github.com/caalvarez1999/DESI-DR1-LRGs . Note that these tables were built using the https://data.desi.lbl.gov/public/dr1/vac/dr1/fastspecfit/iron/v2.1/catalogs/fastspec-iron-main-dark.fits 

| 00 | `pipeline/00_merge_catalogs.py` | Apply the sample selection to the raw catalog, write `PARENT0.fits` |
| 01 | `pipeline/01_query_sparcl_ids.py` | Look up each target's `specid`/`sparclid`/`DR` on SPARCL, write `PARENT1.fits` |
| 02 | `pipeline/02_select_downloadable.py` | Drop targets SPARCL couldn't identify, write `PARENT2.fits` |
| 03 | `pipeline/03_download_spectra.py` | Download spectra for one velocity-dispersion bin, measure Lick indices, write one FITS per galaxy, and fold the result into `PARENT3.fits` |
| 04 | `pipeline/04_cvd_correction.py` | Apply the velocity-dispersion correction to `PARENT3.fits`, write `PARENT4.fits` |

## Internal structure

Every `PARENTX.fits` is a 3-HDU file: an empty `PRIMARY` HDU plus two binary
tables, `METADATA` and `PYLICK`, always row-aligned by `targetid` (row *i* of
`METADATA` and row *i* of `PYLICK` describe the same galaxy). All five stages
share the exact same set of columns end to end -- nothing is ever added or
removed from the schema after stage 00. What changes from one `PARENTX.fits`
to the next is (a) which rows are present, (b) which cells are still the
NaN/empty placeholder stage 00 allocated versus a real measured value, and,
for the `PYLICK` indices only, (c) the storage precision (`E`/float32 up to
`PARENT3.fits`, `D`/float64 in `PARENT4.fits`).

Row count only drops once, at stage 02, and is stable after that:

| file | rows |
|---|---|
| `PARENT0.fits` | 1,855,626 |
| `PARENT1.fits` | 1,855,626 |
| `PARENT2.fits` | 1,728,366 |
| `PARENT3.fits` | 1,728,366 |
| `PARENT4.fits` | 1,728,366 |

### METADATA (19 columns)

| column | format | unit | description |
|---|---|---|---|
| `targetid` | K (int64) | -- | DESI `TARGETID`, primary key |
| `specid` | 17A | -- | SPARCL `specid` |
| `sparclid` | 36A | -- | SPARCL UUID for this spectrum |
| `DR` | 8A | -- | SPARCL data release tag (e.g. `DESI-DR1`) |
| `survey` | 4A | -- | DESI survey (`main`) |
| `program` | 4A | -- | DESI program (`dark`/`bright`) |
| `healpix` | J (int32) | -- | DESI healpix pixel |
| `RA` | D | deg | Right ascension |
| `DEC` | D | deg | Declination |
| `desitarget` | K (int64) | -- | `DESI_TARGET` bitmask (bit 0 = LRG, used in the sample selection) |
| `spectype` | 10A | -- | Redrock `SPECTYPE` (`GALAXY` for every row, by selection) |
| `subtype` | 10A | -- | Redrock `SUBTYPE` |
| `z1` | D | -- | Redshift (fastspecfit `Z`) |
| `vd` | E | km/s | Stellar velocity dispersion (fastspecfit `VDISP`) |
| `dvd` | E | km/s | Uncertainty on `vd`, `1/sqrt(VDISP_IVAR)` |
| `FSF_SFR` | E | Msun/yr | fastspecfit star-formation rate |
| `FSF_Dn4000` | E | -- | fastspecfit narrow 4000A break |
| `FSF_dDn4000` | E | -- | Uncertainty on `FSF_Dn4000`, `1/sqrt(DN4000_IVAR)` |
| `snmedian` | E | -- | Median per-pixel S/N of this galaxy's own downloaded spectrum, in the rest-frame index window |

### PYLICK (69 columns)

| column | format | unit | description |
|---|---|---|---|
| `targetid` | K (int64) | -- | Row identifier, mirrors `METADATA.targetid` |
| `done` | J (int32) | -- | 0/1: has this row's Lick measurement been folded in from a per-galaxy FITS yet |
| `qual` | E/D | -- | 0/1 quality flag from the index measurement (0 = spectrum failed the S/N cut or resolution matching raised; indices left NaN) |
| 32 Lick/D4000 indices + their `d<name>` uncertainty (64 columns) | E up to `PARENT3.fits`, D in `PARENT4.fits` | see below | see below |

Lick/D4000 indices, in `desicc/config.py:LICK_INDICES` order, with the units
from `aux/tableall.dat` (`A` = angstrom equivalent width, `mag` =
magnitude, `break` = dimensionless flux-ratio break):

`OII`(A), `CaIIK`(A), `CaIIH`(A), `D4000`(break), `Dn4000`(break), `HdA`(A), `HdF`(A), `CN1`(mag), `CN2`(mag), `Ca4227`(A), `G4300`(A), `HgA`(A), `HgF`(A), `Fe4383`(A), `Ca4455`(A), `Fe4531`(A), `C24668`(A), `Hb0`(A), `Hb`(A), `OIII5007`(A), `Fe5015`(A), `Mg1`(mag), `Mg2`(mag), `Mgb`(A), `Fe5270`(A), `Fe5335`(A), `Fe5406`(A), `Fe5709`(A), `Fe5782`(A), `NaD`(A), `TiO1`(mag), `TiO2`(mag), `Ha`(A)

Each has a companion `d<name>` column with the same unit, carrying its
Monte-Carlo uncertainty.

### What each stage adds or fills in

- **00** (`00_merge_catalogs.py:build_parent_hdulist`) -- creates both HDUs
  from the raw fastspec-iron catalog, for every row passing `build_mask`
  (LRG & not ELG/QSO, `SPECTYPE == GALAXY`, redshift consistent between
  `METADATA`/`FASTSPEC`, finite `vd`/`dvd` with `VD_MIN < vd < VD_MAX`).
  `METADATA` is fully populated except `specid`/`sparclid`/`DR` (written as
  `""`) and `snmedian` (`NaN`). `PYLICK` is fully allocated with `done=0`,
  `qual=NaN`, and every index/`d`-index set to `NaN` -- no measurement has
  happened yet, this stage only reserves the columns.

- **01** (`01_query_sparcl_ids.py:query_sparcl_ids` + `merge_ids_into_parent`)
  -- looks up `targetid -> specid/sparcl_id/_dr` on SPARCL (cached to
  `PARENT/sparcl_ids_cache.fits`) and left-joins the result onto `METADATA`.
  Only `specid`/`sparclid`/`DR` change; row count and `PYLICK` are untouched.

- **02** (`02_select_downloadable.py:main`) -- boolean-masks `METADATA` and
  `PYLICK` together on "`sparclid` isn't a placeholder"
  (`desicc/catalog.py:is_placeholder`), dropping the targets SPARCL didn't
  recognize (1,855,626 -> 1,728,366 rows in the tables shipped here). No
  column is touched, only rows are removed, so both HDUs stay row-aligned.

- **03** (`desicc/download.py:run` -> `download_and_measure` ->
  `fold_bin_into_parent`) -- run once per velocity-dispersion bin. For every
  row in that bin it downloads the spectrum from SPARCL, computes
  `snmedian` (`desicc/fits_io.py:compute_snmedian`) and measures the 32
  Lick/D4000 indices (`desicc/lick_measurement.py:measure`, wrapping the
  vendored `pylick`/`aux` toolchain), writes one per-galaxy FITS
  (`spectra/singlegalaxies/vd<bin>/<targetid>.fits`, with its own
  `METADATA`+`SPECTRUM`+`PYLICK`), then copies that galaxy's `snmedian`
  into `PARENT3`'s `METADATA` and its `qual`/32 indices/32 uncertainties
  into `PARENT3`'s `PYLICK`, setting `done=1`. Rows already at `done=1` are
  skipped, so this stage is resumable and can be re-run bin by bin;
  `PARENT3.fits` is seeded from `PARENT2.fits` the first time it doesn't
  exist and then updated in place across runs.

- **04** (`04_cvd_correction.py:main` -> `desicc/cvd.py:apply_cvd` ->
  `aux/corrections.py:C_VD`) -- a stateless, full recompute over whatever is
  currently in `PARENT3.fits`: for every row it looks up, from a MILES
  stellar-template grid, the correction (and its uncertainty) appropriate
  for that galaxy's own `vd`, and applies it index by index to all 32
  Lick/D4000 indices and their `d`-index uncertainties, replacing the raw
  values with velocity-dispersion-corrected ones. No columns are added;
  `METADATA` is copied through unchanged and the `PYLICK` index columns are
  upgraded from `E` (float32) to `D` (float64) as they're overwritten.

