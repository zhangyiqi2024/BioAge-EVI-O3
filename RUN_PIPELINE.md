# Pipeline Runner

Working directory: `04_code/` (project root). Structure: `CODE/  DATA/  output/  logs/`

The analysis pipeline is a single sequence numbered **00–42** in analysis order
(see the script inventory in `README.md`). **All paths are relative** — the
package can be placed anywhere on disk.

- Stata steps: 00–01 (data preparation, pre-built), 02, 05–07, 20–26, 28–29, 32–33, 37
- R steps: 03–04, 08–19, 27, 30–31, 34–36, 38–42

## Full pipeline

```bash
cd path/to/final_Submission/04_code     # the folder containing CODE/, DATA/, output/

# Step 1: Stata steps (Stata 17+; adjust the executable name to your install)
StataMP-64.exe /e do CODE/00_master_run.do      # Windows batch
# or, interactively:  do CODE/00_master_run.do

# Step 2: R steps (run after Stata — several steps use its outputs)
Rscript CODE/00_master_run.R
```

Both masters fail fast with a clear message if the working directory is not
the project root. No script changes the working directory.

## Restricted-data steps (automatically guarded)

CHARLS raw data require an approved application (http://charls.pku.edu.cn) and
are **not** part of this package. Steps that need them check for the files and
are **skipped with a message** when absent; every other step runs end-to-end
from the shipped `DATA/analysis_ready/analysis_dataset.dta` alone:

| Step | Needs | When skipped |
|---|---|---|
| 00–01 | raw CHARLS waves (data preparation) | always — `analysis_dataset.dta` is pre-built |
| 37 | raw Harmonized CHARLS D + 2020 Exit/Sample_Info | raw files absent |
| 03 | raw merged `data_生物年龄.dta` | raw file absent |
| 38 | `analysis_dataset_mortality_baseline.dta` (from 37) or raw mortality files | both absent |
| 41 | raw merged/PhenoAge/LAST files | any absent |

The raw files are expected three levels above `04_code/` (i.e.
`<repo>/DATA/...`, the BioAge-EVI-O3 repository layout); with them in place
the guarded steps run automatically.

Notes:
- Steps 00–01 (`00_keep_used_variables.do`, `01_imputation_ipw.do`) are retained
  for documentation; the shipped dataset already contains all imputed values,
  IPW weights and joint exposure groups.
- Step 04 writes Figure S1 one level up into `02_figures/supplementary/`
  (created if missing).
- Step 11 (Figure 3) annotates the DiD rows with the MDE/TOST table produced by
  `45_mde_tost.R`; if that table is absent (e.g. a fresh master-run sequence),
  11 now invokes `CODE/45_mde_tost.R` in a subprocess first, so the figure
  renders identically in any script order (corrective pass 2026-07-24).
- v5 scripts 43–53 run after the masters in the order listed in `README.md`;
  44/46 skip cleanly when the restricted raw files are absent.
- Reproduction caveat: LightGBM outputs (steps 34/36/49-era exploratory ML)
  drift at the ~0.1% level across runs even within lightgbm 4.6.0, and PDF/DOCX
  renders embed creation timestamps — byte-identity is expected for CSV
  outputs and PNG/SVG/TIFF figures, not for PDF/DOCX binaries.

## Key single-analysis runs

```bash
# Primary joint association (four-group contrasts, DiD; MS Table 2 / Fig. 3 source)
Rscript CODE/08_joint_association.R

# Main figures
Rscript CODE/09_figure1_boxplots.R          # Figure 1
Rscript CODE/10_figure2_meta_analysis.R     # Figure 2
Rscript CODE/11_figure3_joint_exposure.R    # Figure 3 (requires 08 outputs)

# Instrumental variables (SI Tables S2, S11 + IV forest/full tables)
Rscript CODE/30_iv_diagnostics.R
Rscript CODE/31_iv_full_tables.R
```

See `RUN_FIGURE3.md` for the Figure 3 end-to-end documentation.

## Packaged outputs

`output/tables/` ships a **curated 29-file subset** of the key result tables
(CSV/DOCX/STER) cited in the manuscript and SI; the copies under
`03_source_data/` and `06_tables/` are byte-identical mirrors of these files
(enforced by `../.docx_fix/final_verify_all.py`). A full pipeline run produces
**70+ files** in `output/tables/` (plus figures and logs) — intermediates,
`.ster` estimate stores and log files that are not individually cited are
re-generated on demand rather than shipped.
