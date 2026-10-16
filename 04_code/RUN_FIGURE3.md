# Figure 3 and Supplementary Figure S9 📊

The current Figure 3 has four panels. Panels a–b show adjusted marginal mean
BAA for four joint exposure patterns at 1-km and 500-m residential buffers
(6,471 observations). Panels c–d show all 18 continuous-model estimates:
16 single/dual-exposure, unweighted/IPW, two-buffer specifications and both
1-km multiple-imputation (MI) estimates. The primary model is the unweighted
1-km dual-exposure model; 500-m models assess buffer sensitivity.

All nonsignificant IPW EVI estimates and the positive 500-m IPW point estimates
remain displayed. Non-MI continuous models use 5,914 observations from 3,154
participants and robust variance. MI uses 20 imputations, model-based variance
and Rubin pooling; its sample counts are not specified in the aggregate table.
Do not assign the primary sample count to MI. MI confidence limits are
converted from native exposure units to the same archived 1-km IQR scale.

Figure S9 preserves all 12 joint contrast rows and both DiD precision rows.
The low-EVI/high-O3 contrast against reference and the O3 contrast within low
EVI are algebraically identical; displaying both does not create independent
evidence. The ±1-year band is an equivalence margin. MDE is an approximate
80%-power quantity, not an exclusion bound; TOST does not establish
equivalence at alpha = 0.05.

## ▶️ Public rendering

Run the complete figure/table pipeline from the package root:

```shell
python 04_code/run_pipeline.py --mode public
```

Or render only these figures to chosen output folders:

```shell
python 04_code/PYTHON/render_joint_sensitivity.py --root . --out output/main --si-out output/supplementary
```

`--out` and `--si-out` may be supplied separately. Outputs are
`Figure3_joint_patterns` and `FigS9_joint_contrasts_precision`, each as editable
PDF/SVG and 600-dpi PNG/TIFF. The public renderer uses only aggregate CSVs and
does not refit models, change P values or select results by significance.

## 📂 Statistical sources

- `CODE/08_joint_association.R`: adjusted means and contrasts, disclosed as
  `03_source_data/Figure3a_adjusted_means.csv` and `Figure3b_contrasts.csv`.
- `CODE/43_primary_continuous_models.R`: all 16 non-MI specifications in
  `03_source_data/primary_continuous_models_v5.csv`.
- `CODE/47_mi_sensitivity.R`: both MI estimates in
  `03_source_data/mi_sensitivity_v5.csv`.
- `CODE/45_mde_tost.R`: DiD precision in
  `03_source_data/Figure3b_DiD_precision.csv`.
- `PYTHON/render_joint_sensitivity.py`: public aggregate-to-display rendering.

The older R-rendered contrast plot remains an analysis intermediate. Figure 2
shows the four core unweighted estimates as association anchors; the complete
sensitivity family is shown in Figure 3. Different panels use different
estimands and sample frames and do not establish effect modification.
