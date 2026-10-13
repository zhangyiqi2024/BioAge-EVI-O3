# Restricted-input and preparation contract

## Reviewer public-reproduction contract

For peer review, the reproducible starting point is the distributed
`03_source_data/` aggregate CSV set and the relative-path entry point
`04_code/run_pipeline.py --mode public`. Reviewers do **not** need the original
MODIS or ChinaHighO₃ raster files. This public path regenerates all four main
figures, eight supplementary figures and 35 physical tables from aggregate data;
it does not refit restricted-data models or rebuild the stored biological-age
score. Figure 2 uses exported density/box geometry. S5 uses the retained
S13/S14 summaries; S6/S7 use source-verified curve/calibration coordinates.
Individual dots in Figure 2, Figure 4c and individual S4 outliers are omitted.
The revised official figures use the same public renderers. S6 explicitly
retains its separately documented historical exploratory model frames.
The machine-readable status is written to
`04_code/output/logs/public_reproduction_status.json`. See
`04_code/PUBLIC_REPRODUCTION.md` for the current reviewer command and outputs.

## What can be reproduced

With the authorised analysis-level input, the package reproduces downstream
association models, precision analyses and sensitivity analyses. English metadata
conversion is executable in PYTHON/prepare_analysis_data.py and preserves every
numeric value, row order and missing value exactly (round-trip assertions).
The codebook contains English labels and category codes, never participant examples.
The public mapping file is intentionally empty; raw-language conversion belongs
to a private access-controlled adapter and is not distributed.

## What is not supplied

This package does not contain the complete original biomarker calibration,
residence linkage, environmental raster extraction or 15-year exposure-pattern
construction workflow. The biological-age outcome and exposure summaries are
supplied upstream inputs; numerical downstream reproduction is not validation of
these constructs or of causal assumptions. Do not replace the outcome or claim
raw-CHARLS-to-paper reconstruction from this package.

## Local conversion

    python PYTHON/prepare_analysis_data.py --input authorised_input.dta --output DATA/analysis_ready/analysis_dataset.dta --audit logs/data_conversion.json

Input and output must differ. The public converter accepts English-named inputs;
unmapped variables or untranslated category labels cause a clear error. This is
deliberate: no category meaning is inferred. Raw-language conversion must occur
in the private author adapter before this public step. All inputs and derived
microdata remain restricted.

To export paper-facing names (for example `stored_baa`, `chronological_age`,
`residential_evi_1km_mean4y` and `ozone_1km_mean4y`), run:

    python PYTHON/export_canonical_data.py

The restricted canonical DTA/CSV is written to DATA/canonical/. The public
canonical_schema_crosswalk.csv maps all 205 canonical names to their internal
analysis names. To use an authorised canonical DTA with the historical scripts:

    python PYTHON/export_canonical_data.py --restore-internal DATA/canonical/analysis_canonical.dta --restored-output DATA/analysis_ready/analysis_dataset_restored.dta

The DTA adapter verifies a one-to-one name map; export asserts exact value and
missingness preservation on round trip. Do not upload either microdata format.

## Optional preparation

For an analysis-level pre-imputation input that already contains the supplied
outcome, lagged exposures and longitudinal exposure groups:
1. Convert names/metadata to DATA/restricted/analysis_input_english.dta.
2. Run CODE/00_keep_used_variables.do to produce selected_input.dta.
3. Run CODE/01_imputation_ipw.do to produce analysis_dataset_prepared_candidate.dta.
4. Compare candidate and archived data before any promotion. Preparation cannot
recreate unavailable upstream outcomes or missing longitudinal exposure fields.
The corrected nnet algorithm is separately verified by CODE/01b_verify_joint_ipw.R;
all three archived joint weights reproduce exactly with numeric-coded covariates,
P(observed category), and R type-7 1st/99th-percentile winsorisation.
Stata candidate calculations use Stata's quantile convention and optimizer and
are not asserted byte-identical to archived R-derived weights.

## Optional restricted derivatives

Place these authorised, English-named files under DATA/restricted/:
- biomarker_merged.dta: merged sample-selection input.
- phenoage_qc.dta: supplied post-quality-control outcome reference.
- master_merged.dta: covariates and exposure history for attrition analysis.
- biomarker_selected.dta: selected biomarker file for split-sample validation.
- merged_all_waves.dta: full merged data for sequential exposure selection.
- charls/Harmonized_CHARLS_D/H_CHARLS_D_Data.dta.
- charls/2020charls/Exit_Module.dta and Sample_Infor.dta.

These are distinct upstream derivatives, not interchangeable copies of the
analysis input. In particular, script 41 expects master_merged.dta greenness
values on the historical 0-10 scale and explicitly divides by 10; analysis-level
EVI is on the 0-1 scale. Script 44 expects the documented CHARLS biomarker names
and units and fits an auxiliary score retaining a cystatin-C albumin proxy.
Its results cannot establish that the stored primary BAA was generated that way.
The exact column selections and transformations are explicit at each loader.
Do not substitute files or infer missing unit conversions from filenames.

Missing optional files cause their documented analyses to be skipped. Existing
archived tables from those analyses are descriptive release artifacts, not
newly reproduced outputs. The exact original ozone input release and checksums
remain unavailable; official CHAP V2 coverage alone does not identify local inputs.
