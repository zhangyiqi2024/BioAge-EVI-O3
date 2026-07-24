# DATA/analysis_ready — restricted data placeholder

**The processed analytical dataset (`analysis_dataset.dta`) is NOT shipped in
this package.** It is derived from individual-level records of the China
Health and Retirement Longitudinal Study (CHARLS), which are third-party
restricted data. Under the CHARLS data use agreement, individual-level data
may not be redistributed to third parties; every user must obtain access
directly from CHARLS. This restriction applies to the data itself, not merely
to personal identifiers, so de-identified versions are also excluded. The
analytical dataset additionally contains spatially resolved environmental
exposures (500-m / 1-km residential buffers) whose spatial detail prevents
full anonymisation.

## How to obtain the data

1. Apply for CHARLS data access at http://charls.pku.edu.cn (free, subject to
   approval by the CHARLS team).
2. Obtain the environmental grids from their public repositories (see the
   manuscript Data Availability statement): MODIS MOD13A1 EVI (NASA LAADS),
   ChinaHighO₃ MDA8 O₃ (National Tibetan Plateau Data Center), nighttime
   lights, DEM, MERRA-2 inversion frequency.
3. Rebuild `analysis_dataset.dta` with the shipped pipeline:
   `CODE/00_keep_used_variables.do` then `CODE/01_imputation_ipw.do`
   (Stata 18). The R master run then reproduces every downstream output from
   the dataset. All scripts run from relative paths and skip gracefully when
   restricted inputs are absent.

`codebook.csv` and `variable_rename_mapping.csv` in this folder describe the
analytical variables and are safe to share.
