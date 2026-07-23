"""Constants and default paths shared by every pipeline stage.

Every path here is relative to the repo root and can be overridden with
CLI flags in the pipeline scripts -- these are just the defaults for a
fresh checkout (see README).
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Stage 00 input: the DESI fastspecfit ("fastspec-iron") VAC, not distributed
# with this repo. METADATA and FASTSPEC in this file are row-aligned (one
# row per target, same order in both HDUs) -- there is no join to do here.
RAW_CATALOG = ROOT / "tablasDESI" / "v2.1" / "fastspec-iron-main-dark.fits"

# Velocity-dispersion range the sample is restricted to, applied in stage 00.
# Bins outside this range (e.g. the historical vd160180/vd180200) aren't part
# of the working sample -- see the stacking analysis this pipeline feeds.
VD_MIN = 200.0  # km/s, exclusive
VD_MAX = 400.0  # km/s, exclusive

PARENT_DIR = ROOT / "TABLAS" / "PARENT"
PARENT_MERGED = PARENT_DIR / "PARENT0.fits"        # stage 00 output
PARENT_WITH_IDS = PARENT_DIR / "PARENT1.fits"      # stage 01 output
PARENT_DOWNLOADABLE = PARENT_DIR / "PARENT2.fits"  # stage 02 output
PARENT_FINAL = PARENT_DIR / "PARENT3.fits"         # stage 03 output, raw indices
PARENT_FINAL_CVD = PARENT_DIR / "PARENT4.fits"     # stage 04 output, CVD-corrected

SPARCL_IDS_CACHE = PARENT_DIR / "sparcl_ids_cache.fits"

SPECTRA_ROOT = ROOT / "spectra"                        # every spectrum in this repo is DESI's
SINGLEGALAXIES_ROOT = SPECTRA_ROOT / "singlegalaxies"  # stage 03 output, one FITS per galaxy
STACKS_ROOT = SPECTRA_ROOT / "stacks"                  # stage 10 output, one FITS per stack

STACKS_TABLE_DIR = ROOT / "TABLAS" / "STACKed"         # stage 11/12 output (summary tables)

def stacked_table_paths(vd_bin: str, archaeology: str, cosmology: str, target_sn: float) -> tuple[Path, Path]:
    """(raw, CVD-corrected) paths for a group's stage 11/12 summary table."""
    sn_tag = f"{round(target_sn):03d}"
    out_dir = STACKS_TABLE_DIR / f"vd{vd_bin}" / archaeology / cosmology
    base = f"stackedDESI_vd{vd_bin}_SN{sn_tag}_{archaeology}_{cosmology}"
    return out_dir / f"{base}.fits", out_dir / f"{base}_CVD.fits"

# Velocity-dispersion (CVD) correction grids: index correction as a function
# of imposed velocity dispersion, precomputed on MILES stellar templates.
CVD_DIR = ROOT / "aux" / "CVD" / "MILES"
CVD_C_FUNCS = CVD_DIR / "Cfuns_MILESstars.npy"
CVD_DC_FUNCS = CVD_DIR / "dCfuns_MILESstars.npy"
CVD_UVES_MAX = 500.0   # km/s, exclusive
CVD_UVES_STEP = 0.1    # km/s, grid spacing the C/dC funcs were computed at

# Archaeological scaling relation (age(vd)) + cosmology (age -> redshift),
# used to work out the maximum redshift a galaxy can have and still be used
# in a given vd-group's stack. See desicc/stack_selection.py.
ARCHAEOLOGY_DIR = ROOT / "archaeology"
COSMOLOGY_CHOICES = ("PLANCK", "RIESS")
ARCHAEOLOGY_CHOICES = ("ALVAREZ", "JOHANSSON", "THOMAS")

# Stacking: rest-frame window used to normalize each spectrum before
# combining, and default target S/N for the stacks (see desicc/stacking.py).
STACK_NORM_WINDOW_RF = (4200.0, 4400.0)
DEFAULT_TARGET_SN = 100.0

# TMJ (Thomas, Maraston & Johansson 2011) stellar-population-synthesis fit
# to a stack's Lick indices, see desicc/tmj.py and
# pipeline/20_SPSfitTMJ.py / pipeline/21_posteriorsTMJ.py.
MODELLED_TABLE_DIR = ROOT / "TABLAS" / "modelled"  # stage 21 output (summary tables)
CHAINS_ROOT = ROOT / "chains" / "TMJ"              # stage 20 output (MCMC chains)

TMJ_MODEL_DIR = ROOT / "aux" / "TMJ"
DEFAULT_TMJ_GRID = "custom250924b"  # only grid vendored in aux/TMJ/ -- see README

def tmj_model_path(grid: str = DEFAULT_TMJ_GRID) -> Path:
    return TMJ_MODEL_DIR / f"TMJ_MILES_{grid}.npy"

def tmj_tabla_path(grid: str = DEFAULT_TMJ_GRID) -> Path:
    return TMJ_MODEL_DIR / "tabla" / f"TMJ_MILES_{grid}_tabla.npy"

def modelled_table_path(vd_bin: str, archaeology: str, cosmology: str, target_sn: float,
                         grid: str = DEFAULT_TMJ_GRID) -> Path:
    """Path for a group's stage 21 (TMJ posteriors) summary table."""
    sn_tag = f"{round(target_sn):03d}"
    out_dir = MODELLED_TABLE_DIR / f"vd{vd_bin}" / archaeology / cosmology
    return out_dir / f"posteriorsDESI_TMJ_{grid}_vd{vd_bin}_SN{sn_tag}_{archaeology}_{cosmology}.fits"

# Lick indices the TMJ fit itself is run on -- a fixed choice from the
# original analysis (unrelated to the ARCHAEOLOGY scaling-relation choice
# used elsewhere in this pipeline, despite the similar name of the original
# variable it came from: `list_johansson` in 9_FULLTMJ.py).
TMJ_FIT_INDICES = [
    "HdA", "HdF", "G4300", "HgA", "HgF", "Fe4383", "Fe4531", "Hb", "Mg2",
    "Mgb", "Fe5270", "Fe5335", "Fe5406",
]

# The 25-index axis TMJ_MODEL_DIR's model grids are tabulated on, in
# aux.emcee_lick.Lick_indices order -- span_pylick=25 in aux/mcmc_tmj.py.
TMJ_MODEL_INDICES = [
    "HdA", "HdF", "CN1", "CN2", "Ca4227", "G4300", "HgA", "HgF", "Fe4383",
    "Ca4455", "Fe4531", "C24668", "Hb", "Fe5015", "Mg1", "Mg2", "Mgb",
    "Fe5270", "Fe5335", "Fe5406", "Fe5709", "Fe5782", "NaD", "TiO1", "TiO2",
]

# Systematic-error inflation applied to Lick index errors before the TMJ
# fit, carried over unchanged from 9_FULLTMJ.py.
TMJ_ERROR_INFLATION = 10.0

DEFAULT_TMJ_NWALKERS = 50
DEFAULT_TMJ_NITER = 10000

# Cosmographic (Taylor-expansion) fit of t(z) across vd groups, see
# desicc/cosmographic.py and pipeline/30_cosmographic.py.
COSMOGRAPHIC_CHAINS_ROOT = ROOT / "chains" / "cosmographic"
COSMOGRAPHIC_PRIORS_CHOICES = ("narrow", "wide")
DEFAULT_COSMOGRAPHIC_PRIORS = "narrow"
DEFAULT_COSMOGRAPHIC_ORDER = 3
DEFAULT_COSMOGRAPHIC_NWALKERS = 40
DEFAULT_COSMOGRAPHIC_NITER = 100000

# Lick/D4000 indices measured on every spectrum, in pylick's naming.
LICK_INDICES = [
    "OII", "CaIIK", "CaIIH", "D4000", "Dn4000", "HdA", "HdF", "CN1", "CN2",
    "Ca4227", "G4300", "HgA", "HgF", "Fe4383", "Ca4455", "Fe4531", "C24668",
    "Hb0", "Hb", "OIII5007", "Fe5015", "Mg1", "Mg2", "Mgb", "Fe5270",
    "Fe5335", "Fe5406", "Fe5709", "Fe5782", "NaD", "TiO1", "TiO2", "Ha",
]

SPARCL_SPECTRUM_FIELDS = [
    "sparcl_id", "specid", "targetid",
    "flux", "wavelength", "ivar", "mask", "wave_sigma",
]

# Rest-frame wavelength window used for index measurement (Angstrom).
RESTFRAME_WINDOW = (3600.0, 6700.0)

DEFAULT_NSAMPLING = 100
DEFAULT_APPLY_RESOLUTION = True
DEFAULT_SIGMA_RESOLUTION = 2.5 / 2.355  # target resolution (sigma), Angstrom

DEFAULT_BATCH_SIZE_SPARCL_FIND = 500
DEFAULT_BATCH_SIZE_SPARCL_RETRIEVE = 50
DEFAULT_MAX_RETRIES = 20
DEFAULT_RETRY_SLEEP = 1.0  # seconds
