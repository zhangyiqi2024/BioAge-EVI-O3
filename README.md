# BioAge, residential greenness and ozone 🔬

Code and aggregate results for **Joint associations of residential greenness
and ozone with biological ageing**. This is the public GitHub release branch `release-2026-10-13-fig3`.
The historical Zenodo archive is not updated by this upload.

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

- `04_code/output/main/`: main Figures 1–4, PNG/TIFF/PDF; Figure 4 also has vector SVG.
- Figure 3 combines adjusted four-group means, formal contrasts and DiD precision diagnostics at both residential buffers.
- `04_code/output/supplementary/`: Supplementary Figures S1–S8.
- `04_code/output/tables/`: Table 1 and 34 SI table panels, as HTML/CSV/JSON;
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
`04_code/DATA/BUILD_SCOPE.md`. The original primary BAA construction/calibration
chain and exact local ozone release remain incompletely documented. Downstream
plot reconstruction does not resolve those scientific provenance limitations.
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
