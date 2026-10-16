# GitHub reviewer-materials package

Release: `2026.10.17-reviewer-materials`  
Branch: `release-2026-10-17-reviewer-materials`  
Tag: `v2026.10.17-reviewer-materials`  
Repository: https://github.com/zhangyiqi2024/BioAge-EVI-O3

This release contains the public materials needed to inspect and reproduce the
reported aggregate figures and tables. The manuscript, supplementary
information, cover letter, author declarations, and other journal-upload
documents are intentionally kept outside GitHub in the local submission
package. They are excluded from this release by design.

## Start here

From the package root, install the public Python dependencies and run:

```shell
python -m pip install -r 04_code/requirements_public.txt
python 04_code/run_pipeline.py --mode public
```

The public path uses only `03_source_data/` and `04_code/`. It does not
download data, access the network, refit participant-level models, or require
R, Stata, Word, CHARLS microdata, or original environmental rasters.

## Package map

| Directory | Reviewer use |
|---|---|
| `02_figures/` | Main and supplementary publication figures in raster and vector formats |
| `03_source_data/` | Aggregate-only drawing inputs, table sources, schemas, and provenance |
| `04_code/` | End-to-end relative-path analysis code, public reproduction entry point, and requirements |
| `05_documentation/` | Figure/table crosswalks and public-source documentation |
| `06_tables/` | Curated aggregate result tables and presentation-ready table files |
| `RELEASE_MANIFEST.json` | Byte counts and SHA-256 hashes for every published file |
| `SHA256SUMS.txt` | Plain-text checksum list for quick integrity checks |

## Evidence map

`05_documentation/reviewer_analysis_map.csv` classifies the published outputs as
primary, supportive, or exploratory and maps each figure/table to its public
source and rendering script. The aggregate source data preserve displayed
precision and provenance; they do not expose individual observations.

The primary evidence is the mutually adjusted continuous-exposure analysis.
Four-group joint patterns, buffer comparisons, and precision analyses are
descriptive or supportive. IV, mediation, machine-learning, spline, and other
extended modules remain exploratory or supportive and should not be read as
additional confirmatory tests.

## Data access and limits

CHARLS participant records and the original environmental rasters remain under
their source access agreements and are not redistributed. The public tables
and figure coordinates are author-provided aggregates intended to reproduce the
displayed outputs. They do not support reconstruction of individual-level
models or independent validation of upstream exposure processing.

EVI is described as EVI within residential buffers throughout this package;
the release does not use “residential EVI” or “residential greenness” as a
variable name. The study uses repeated measurements from 2011 and 2015, so
results are repeated-measure or longitudinal associations with the stated
temporal and causal limits.

## Licence

Original code and documentation are released under the MIT licence in
`LICENSE`. Third-party data and software retain their own licences. No rights
to redistribute restricted participant records are granted.
