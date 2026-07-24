# Figure 3 — End-to-end reproducible pipeline

## Pipeline overview

```
CODE/08_joint_association.R           (R ≥ 4.0; primary joint-association models)
        │
        ├─► output/tables/joint_association_adjusted_means_v4.csv
        └─► output/tables/joint_association_contrasts_v4.csv
                │
                ▼
CODE/11_figure3_joint_exposure.R      (R ≥ 4.0)
                │
                ├─► output/figures/Figure3_joint_exposure_enhanced.svg    (primary, editable text)
                ├─► output/figures/Figure3_joint_exposure_enhanced.pdf    (cairo_pdf)
                ├─► output/figures/Figure3_joint_exposure_enhanced.tiff   (600 dpi LZW)
                └─► output/figures/Figure3_joint_exposure_enhanced.png    (300 dpi preview)
```

## Step 1 — Run the joint-association models (R)

```bash
cd "E:/0_Research/0_Paper/BioAge_EVI_O3/Submission/final_Submission/04_code"
Rscript CODE/08_joint_association.R
```

IPW-weighted linear mixed models (WeMix; participant random intercept; robust
variance) with marginal standardization. This generates:
- `output/tables/joint_association_adjusted_means_v4.csv`
- `output/tables/joint_association_contrasts_v4.csv`

## Step 2 — Render the enhanced two-panel Figure 3 (R)

```bash
cd "E:/0_Research/0_Paper/BioAge_EVI_O3/Submission/final_Submission/04_code"
Rscript CODE/11_figure3_joint_exposure.R
```

This reads the CSVs from Step 1 and generates all four formats above.

## R package dependencies

```r
install.packages(c("ggplot2", "dplyr", "patchwork", "scales", "svglite", "ragg"))
```

## Design notes (nature-figure v2.0 compliant)

- **Archetype**: `quantitative grid` — two panels of aligned numerical evidence
- **Dimensions**: 183 × 210 mm (Nature double-column)
- **Font**: Arial, 7 pt base (axis text 6.5 pt, legend 6.2 pt, tags 9 pt bold)
- **Palette**: Wong 2011 colourblind-safe (Nature Methods)
- **Caption**: **Not embedded on the figure** — printed to console for manual insertion into manuscript text. Follows Nature legend conventions (Fig. N | title → panel descriptions → statistics → Source Data boilerplate).

## Figure 3 panel description

| Panel | Content | Data source |
|---|---|---|
| **a** | Adjusted marginal mean BAA (years) by joint EVI × O₃ group. Dots = point estimates; error bars = 95% CI. Two facets (1 km / 500 m buffers). n = 6,471. | `joint_association_adjusted_means_v4.csv` |
| **b** | Forest plot of contrasts vs. reference group (Low EVI / low O₃) plus DiD interaction. Horizontal bars = 95% CI. Numeric annotations = estimate [95% CI]. Orange = p < 0.05, grey = p ≥ 0.05. | `joint_association_contrasts_v4.csv` |

## Manuscript figure legend

The R script prints a Nature-style legend to the console. Copy it into the manuscript:

> **Fig. 3 | Adjusted BAA across joint EVI × O₃ exposure-pattern groups.**
> **a,** Model-adjusted mean BAA across the four joint exposure-pattern groups at 1-km and 500-m residential buffers. Points are covariate-adjusted marginal means from linear mixed-effects models; error bars denote 95% confidence intervals (n = 6,471 observations). Marker shape and colour identify the exposure-pattern group.
> **b,** Forest plot of adjusted contrasts versus the reference group (Low EVI / low O₃) plus the difference-in-differences interaction, which estimates departure from additivity of mean BAA differences. Error bars denote 95% CI. \*p < 0.05. These are adjusted associations and should not be interpreted as evidence that EVI causally offsets O₃-related BAA. Source data are provided as a Source Data file.
