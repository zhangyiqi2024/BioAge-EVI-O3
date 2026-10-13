# Public aggregate source data 📊

This directory and sibling `04_code/` are sufficient for reviewer reconstruction.
Run `python 04_code/run_pipeline.py --mode public` from the package root, or
launch that entry point by path from any working directory.

- Figure prediction grids, coefficients, density polygons, quartiles, calibration
  points and counts are stored as CSV. They are aggregate drawing inputs, not
  participant-level observations.
- `FIGURE_SOURCE_CODEBOOK.md`, `figure_source_schema.json` and
  `figure_source_provenance.json` document the recovered distribution/S6/S7
  sources and their exact model/sample definitions.
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
