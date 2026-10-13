# Public aggregate source data 📊

This directory and sibling `04_code/` are sufficient for reviewer reconstruction.
Run `python 04_code/run_pipeline.py --mode public` from the package root, or
launch that entry point by path from any working directory.

- Figure prediction grids, coefficients, density polygons, quartiles, calibration
  points and counts are stored as CSV. They are aggregate drawing inputs, not
  participant-level observations.
- `FIGURE_SOURCE_CODEBOOK.md`, `figure_source_schema.json` and
  `figure_source_provenance.json` document the recovered distribution/S6/S7
  sources and their exact model/sample definitions. Their `current_composition`
  mappings describe current Figure 2, Figure 3 and Figure S9 use while retaining
  historical export provenance.
- `tables/README.md` explains the 35 displayed table sources. Each CSV has
  JSON layout/notes/provenance. `tables/manifest.json` includes checksums.
- Other CSVs retain their historical analysis-result filenames (`v4`, `v5`,
  `nested`) for traceability. Those suffixes are not claims of current versions
  or nested cross-validation. EVI is unitless; O3 uses micrograms per cubic
  metre; BAA and BAA contrasts use study-defined years. Model-specific scales
  such as log exposure are stated in table notes.

Distribution layers omit participant points and individual outliers. Suppressed
histogram counts of 1–4 are missing rather than zero. Missing curve limits or
contour cells preserve the original evaluated range/mask.

Reconstruction preserves the displayed archived values; it does not independently
recalculate model estimates. Eight tables lack a separate historical result CSV
and are explicitly marked `documented_displayed_values` in their provenance.
Original calibration of the stored primary BAA and exact local ozone release
remain incomplete provenance issues, outside this public drawing path.

The source results are author-provided aggregates; cite the study on reuse.
No standard open-data licence is asserted for third-party data. The MIT licence
in `04_code/LICENSE` covers original code and documentation only. Individual
CHARLS records and original environmental rasters are not distributed.

## 📊 Current panel sources

The public pipeline renders 13 figures (4 main and 9 supplementary). Figure 2a
reuses the 1-km exposure-bin layers historically named for Figure 4. Figure 2b
uses four unweighted dual-exposure results: the 1-km primary model and 500-m
buffer sensitivity. Figure 2c–d reuse archived density/box geometry and mean
points for all four groups (6,472 observations per buffer). Source mean CI
columns are preserved but not plotted because they ignore repeated observations
within participants. Suppressed counts stay missing and are not reconstructed.

Figure 3a–b use the 6,471-observation joint-model means. Figure 3c–d use all
16 rows of `primary_continuous_models_v5.csv` and both rows of
`mi_sensitivity_v5.csv`. MI native-unit confidence limits are scaled by the
archived 1-km IQR without changing nominal P values. The non-MI frame is 5,914
observations from 3,154 participants; MI counts are unspecified in the source.
Figure S9 uses all rows of `Figure3b_contrasts.csv` and
`Figure3b_DiD_precision.csv`, retaining the interaction results and precision
limitations. Filenames remain unchanged for provenance traceability.

