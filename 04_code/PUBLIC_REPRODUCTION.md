# BioAge, EVI within residential buffers and ozone 🔬

Code and aggregate results for **Repeated-measure associations of EVI within residential buffers and ozone with biological age acceleration in China**. The revision release branch is
`release-2026-10-13-fig3` in `zhangyiqi2024/BioAge-EVI-O3`; use the Git commit and
`RELEASE_MANIFEST.json` to identify a particular release. The historical Zenodo
archive is not changed by running this pipeline.

## ▶️ One-command reviewer reproduction

Use Python 3.11 and the packages in `04_code/requirements_public.txt`:

```shell
python -m pip install -r 04_code/requirements_public.txt
python 04_code/run_pipeline.py --mode public
```

The public mode is also the default. The entry point discovers its own location,
so it works when launched from any working directory. Only `03_source_data/`
and `04_code/` are required. MODIS/ChinaHighO3 rasters, CHARLS participant records,
R, Stata, Word, manuscript files and archived figure images are not required.
The two directories must remain siblings. No data download or network access
occurs during reproduction.

## 📊 Outputs and scope

- `04_code/output/main/`: main Figures 1–4, PNG/TIFF/PDF; Figures 2–4 also have vector SVG.
- `04_code/output/supplementary/`: Supplementary Figures S1–S9; S9 also has editable PDF/SVG.
- `04_code/output/tables/`: Table 1 and 35 SI table panels, as HTML/CSV/JSON;
  open `index.html` or print `all_tables.html` to PDF.
- `04_code/output/logs/`: per-step logs, software versions, run status and
  `public_reproduction_status.json` with source/output hashes and cell checks.

An alternative output directory can be set with `--output-dir PATH`.
Excel is optional: install openpyxl and run
`python 04_code/PYTHON/render_public_tables.py --xlsx`.

This path reconstructs the displayed results from disclosed aggregate source
data. It does not refit participant-level models or independently validate the
archived estimates. Curves use evaluated prediction grids; distribution panels
use density/box summaries and binned exposures. Individual observations and
outliers are not distributed. Table CSV/JSON preserve displayed precision and
footnotes; each has provenance distinguishing an archived result file from a
value documented in the manuscript. This distinction also applies to the eight
tables lacking a separate machine-readable historical result file.

S5 now uses the retained 5,914-observation summaries in Tables S13–S14. S6
explicitly distinguishes its 2,785-person baseline age/sex-adjusted GAM curves
from its historical 6,002-observation quantile g-computation sample (ages 45–94).
Neither is the primary 5,914-observation model.

## 📂 Data and code guide

`03_source_data/README.md` describes the English source dictionaries and
provenance records. `03_source_data/tables/manifest.json` lists all physical
tables and checksums. `05_documentation/public_table_crosswalk.csv` and
`04_code/PYTHON/audit_public_reproduction.py` map paper items to sources.
`RELEASE_MANIFEST.json` records every release file. Files in `02_figures/` and
`06_tables/` are comparison archives, not public-run dependencies.

Optional restricted-data analysis is separate:

```shell
python 04_code/run_pipeline.py --mode analysis --stata StataMP-64.exe --rscript Rscript
```

It requires independently authorised inputs specified in
`04_code/DATA/BUILD_SCOPE.md`. The analysis mode requires independently authorised restricted inputs described
in `04_code/DATA/BUILD_SCOPE.md`. The public mode uses only the supplied aggregate
files and does not refit participant-level models or redistribute restricted records.
R 4.4.3 on the audit workstation has a documented UCRT exit failure; the separate
public Python path is unaffected. See `04_code/LOCAL_VERIFICATION.md`.

The LightGBM split is participant-level, but initial imputation and climate PCA
precede splitting. The historical `nested` filenames do not establish fully
nested validation. MDE is a power quantity, not an exclusion bound.

## 🔐 Licence and data access

Code/documentation: MIT (see LICENSE). Public tables/figure coordinates are
author-provided aggregate study results; cite this study when reusing them.
Third-party data retain their own licences. No individual-level DTA/CSV,
identifiers, row-level predictions or SHAP values are included. Obtain CHARLS
microdata under its own access agreement. Use this curated release directory
for upload, not the local parent workspace containing restricted data.

The structure follows the single-entry-point and figure/table crosswalk
conventions in Yusaku Horiuchi's replication-package-guide and FAIR principles.

The package reconstructs 13 figures and 36 physical table panels.

## 📊 Current figure composition

Figure 2 combines the 1-km joint exposure distribution and marginal histograms,
four unweighted dual-exposure estimates, and observed BAA distributions at both
buffers. Only the 1-km estimates are primary; 500-m estimates assess buffer
sensitivity. Descriptive mean points have no naive mean confidence intervals;
the source interval columns remain unchanged. Suppressed counts of 1–4 remain
distinct from zero, so displayed cell counts need not sum to the model-frame
total. Exposure coverage and continuous models use 5,914 observations; the
descriptive BAA panels use 6,472 observations per buffer.

Figure 3 retains the adjusted joint-pattern means (6,471 observations) and all
18 continuous-model sensitivity estimates: 16 non-MI specifications and both
1-km MI rows. All nonsignificant IPW EVI results remain displayed. The non-MI
models use robust variance and 5,914 observations from 3,154 participants;
MI uses model-based variance and Rubin pooling across 20 imputations, with its
sample counts unspecified in the aggregate table. Figure S9 contains all 12
joint contrast rows and both DiD precision rows, including nominal P, MDE and
TOST. These displays do not establish effect modification. No significance
threshold, estimate or nominal P value is changed.

The new renderers are `04_code/PYTHON/render_primary_figure.py` for Figure 2
and `04_code/PYTHON/render_joint_sensitivity.py` for Figure 3 and Figure S9.
They produce editable PDF/SVG and 600-dpi PNG/TIFF from public aggregate inputs.

