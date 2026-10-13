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
reconstructs 13 figures (4 main and 9 supplementary) and 35 physical tables without fitting models or
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
The revision release branch is `release-2026-10-13-fig3` in
`zhangyiqi2024/BioAge-EVI-O3`. The default GitHub branch and historical Zenodo
archive represent earlier releases. Running this code does not publish a release
or update a DOI version.

## Privacy and licence

The code is MIT-licensed (see the release LICENSE). Restricted CHARLS-derived
records are not covered by that licence and cannot be redistributed by this
package. Do not publish DATA/*.dta, DATA/restricted/, logs, per-person predictions,
per-person SHAP values, or local reports containing participant identifiers.

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

