# Analysis reproduction 🔬

This code accompanies **Joint associations of residential greenness and ozone with biological ageing**.

## Scope and inputs

The primary analysis is the mutually adjusted continuous EVI/O3 model in
CODE/43_primary_continuous_models.R. The supplied analysis-level outcome is
BAA_PhenoAge, a stored study-specific biological-age acceleration measure.
The original calibration is unavailable; the variable name does not verify a
standard PhenoAge implementation or distinguish residualisation from subtraction.
The dataset has 6,568 observations before model-specific filters. The
primary complete-case frame has 5,914 observations / 3,154 participants.
Metadata translation does not recalculate outcomes, weights or exposure values.

The GitHub release excludes individual-level CHARLS data and row-level
prediction/SHAP outputs. Each user must obtain CHARLS permission independently.
See DATA/BUILD_SCOPE.md for required files and reconstruction limits.
An authorised local analysis input belongs at
DATA/analysis_ready/analysis_dataset.dta.

## Run

From the folder containing this README:

    python run_pipeline.py --mode analysis --stata StataMP-64.exe --rscript Rscript

Executable paths may be supplied explicitly. Paths for inputs, outputs and
scripts are relative to this folder, and the entry point discovers its own
location. The command executes the Stata master, verifies corrected joint
weights, and executes all R analyses through 54. A failed process or Stata
error yields a non-zero exit. Optional restricted inputs are reported as SKIP.
The run logs, exact versions, status table and elapsed time are under logs/.

For an environment preflight only:

    python run_pipeline.py --mode check --rscript Rscript

For reviewer reproduction of every displayed figure/table from public aggregates:

    python run_pipeline.py --mode public

Only this code folder and sibling `03_source_data/` are needed. Install
`requirements_public.txt` in Python 3.11. Outputs are `output/main`,
`output/supplementary`, `output/tables` and `output/logs`. The public path
reconstructs 12 figures and 35 physical tables without fitting models or
rebuilding the primary BAA. See PUBLIC_REPRODUCTION.md.

## Environment

Audited on Windows with Python 3.13, R 4.4.3 and Stata 18 MP.
See LOCAL_VERIFICATION.md: a workstation R/UCRT failure prevents certification
of a clean complete R run. The public Python figure path is tested separately.
Use an existing conda environment; no new virtual environment is required.
Restricted-data utility dependencies are in requirements_python.txt; the public path uses requirements_public.txt. Supplementary plotting
uses requirements_plotting.txt in an existing conda environment with Matplotlib.
Direct R dependencies
and exact tested versions are listed in requirements_r_versions.csv; requirements.R
installs missing packages from CRAN. The version manifest is a record of the
tested environment, not an automatic historical-version installer.
Stata uses estout/esttab, coefplot, ivreg2, ranktest and standard commands;
see STATA_REQUIREMENTS.md. Stochastic steps set seed 20240101 unless explicitly
documented by a script. Version/thread-dependent LightGBM drift remains possible.

## Current results and filenames

The original script numbers and output filenames are retained so older source
CSVs remain resolvable. A filename containing v4/v5/nested is historical metadata:
it is not a claim of a current version, nested cross-validation or causal validity.
The LightGBM validation uses participant-level splitting, but initial imputation
and climate PCA precede the split. MDE is a power quantity, not a confidence bound.

The primary exposure inputs represent four preceding years. The study records
describe CHAP ozone inputs as 1-km grids, and the official CHAP V2 archive covers
2000 onward, but the exact locally used release/checksums have not been recovered.
Do not infer raw-grid provenance from the analysis buffer labels (500 m / 1 km).
The historical GitHub/Zenodo release mentioned in the manuscript predates this
local revision; no upload or DOI version update is performed here.

## Privacy and licence

The code is MIT-licensed (see the release LICENSE). Restricted CHARLS-derived
records are not covered by that licence and cannot be redistributed by this
package. Do not publish DATA/*.dta, DATA/restricted/, logs, per-person predictions,
per-person SHAP values, or local reports containing participant identifiers.
