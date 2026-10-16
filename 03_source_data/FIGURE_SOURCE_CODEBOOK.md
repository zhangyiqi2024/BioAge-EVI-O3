# Aggregate figure-source codebook

These additional source files enable public figure rendering using only `03_source_data` and `04_code`.
They contain aggregated drawing geometry and statistics, with no individual records. Missing histogram counts flagged as suppressed represent counts of 1-4. Other missing values preserve original plot boundaries or masks.

The existing Figure 4 curve and masked-contour CSVs remain authoritative. This codebook describes the newly exported supporting layers, S6/S7 geometry and the formal raincloud geometry.

## Provenance and interpretation

- S6 GAM curves use 2785 baseline observations, age/sex adjustment and predictions for a male aged 58.8 years. They are exploratory, not the primary multivariable models.
- S6 qgcomp weights use 6002 observations aged 45-94; the original exploratory model includes 88 observations outside the primary age range. Full-precision coefficients independently reproduce the displayed -49.6%/+50.4% weights.
- S7 observed calibration proportions are crude death-indicator means, not censoring-adjusted Kaplan-Meier risk estimates.
- Individual raincloud dots and Figure 4 scatter points are omitted; Figure 4 exposure distributions use aggregate binned counts. Counts 1-4 are suppressed. Statistical curves, boxplots, means and confidence intervals preserve formal source geometry.
- Author-side R processes produced EXPORT_COMPLETE and validated exports, then returned Windows 3221225477 (0xc0000005). This is not a clean R process exit. Public redraw reads CSVs and does not invoke these R fits.

## FigS6_archived_display_weights.csv

Rounded signed percentages shown in the archived S6 left panel. Values are independently verified by the full-precision qgcomp output in FigS6_signed_weights.csv.

Rows: 2. SHA-256: `c413be82dea964f42b4dc38dec1ddc54c2baa6c1d9c4eff6306afd7b108621a0`.

| Variable | Type | Missing | Description |
|---|---|---:|---|
| exposure | str | 0 | Exposure name: EVI or O3. |
| weight | float64 | 0 | Rounded signed contribution as a percentage of the sum of absolute exposure coefficients. |

## FigS6_gam_curve_coordinates.csv

S6 right panels: single-exposure GAM curves at age 58.8 and female=0; baseline 2011, ages 46-89, N=2785 per model. x is EVI (unitless) or ozone (micrograms per cubic metre); fit and intervals are BAA years.

Rows: 200. SHA-256: `6fd588eeca984c7be0810bc1a28b9e8b09776623f9cebf46f6a5bdf8d184eb90`.

| Variable | Type | Missing | Description |
|---|---|---:|---|
| exposure | str | 0 | Exposure name: EVI or O3. |
| x | float64 | 0 | Exposure value or categorical plotting coordinate, as specified for the file. |
| fit | float64 | 0 | Model-fitted mean BAA, in years. |
| lci | float64 | 0 | Lower pointwise 95% confidence limit on the response scale. |
| uci | float64 | 0 | Upper pointwise 95% confidence limit on the response scale. |
| se | float64 | 0 | Standard error of the corresponding estimate on its response scale. |
| model | str | 0 | Model specification in English. |

## FigS6_gam_metadata.csv

S6 GAM sample and smoothing specification. The 3197 baseline observations include 412 with missing exposure, leaving 2785 in each model. GAMs adjust for age and sex, use thin plate splines with k=10 and GCV selection, and are distinct from the primary models.

Rows: 2. SHA-256: `31a2b6c0ce4e506635f3d51dbb2de7f6bfc8f63d76959742ebe164a8b42d010c`.

| Variable | Type | Missing | Description |
|---|---|---:|---|
| exposure | str | 0 | Exposure name: EVI or O3. |
| n | int64 | 0 | Number of observations or participants in the documented frame. |
| baseline_year | int64 | 0 | Calendar year used for baseline selection. |
| prediction_age | float64 | 0 | Age fixed for model predictions, in years. |
| prediction_female | int64 | 0 | Female indicator fixed for predictions: 0=male, 1=female. |
| smooth_basis | str | 0 | GAM smoothing basis. |
| smooth_k | int64 | 0 | Basis dimension supplied to the smoother. |
| smooth_edf | float64 | 0 | Estimated effective degrees of freedom of the exposure smoother. |
| estimation_method | str | 0 | GAM smoothing-parameter selection method. |

## FigS6_qgcomp_metadata.csv

Sample metadata for the legacy exploratory qgcomp panel; its sample differs from the primary N=5914 frame. Ages are in years and counts are observations.

Rows: 1. SHA-256: `7df4ba078e9f3c56d33e316ed7f1f0230ac20a15a89df4a7a9ccf62abce1ae2e`.

| Variable | Type | Missing | Description |
|---|---|---:|---|
| n | int64 | 0 | Number of observations or participants in the documented frame. |
| age_min | int64 | 0 | Minimum observed age, in years. |
| age_max | int64 | 0 | Maximum observed age, in years. |
| n_outside_primary_age_range | int64 | 0 | Number of observations with age outside 46-89 years. |
| n_2011 | int64 | 0 | Number of observations in the 2011 wave. |
| n_2015 | int64 | 0 | Number of observations in the 2015 wave. |
| exposure_quantiles | int64 | 0 | Number of exposure quantile categories. |
| age_filter_applied | bool | 0 | Whether the primary age-range filter was applied to this exploratory model. |

## FigS6_signed_weights.csv

S6 left panel: fully adjusted dual-exposure quantile g-computation, q=4; signed coefficient divided by sum of absolute coefficients, times 100. N=6002 ages 45-94; the legacy exploratory script does not apply the primary age filter, so 88 observations lie outside ages 46-89.

Rows: 2. SHA-256: `cfa050f645227b96cda23d9c90781910730e2e0337657ae49197485a8f653852`.

| Variable | Type | Missing | Description |
|---|---|---:|---|
| exposure | str | 0 | Exposure name: EVI or O3. |
| coefficient | float64 | 0 | Change in BAA years per one exposure-quartile score increase. |
| se | float64 | 0 | Standard error of the corresponding estimate on its response scale. |
| ci_lower | float64 | 0 | Lower normal-approximation 95% confidence limit, in BAA years. |
| ci_upper | float64 | 0 | Upper normal-approximation 95% confidence limit, in BAA years. |
| p_value | float64 | 0 | Two-sided model p-value. |
| signed_weight_percent | float64 | 0 | Full-precision signed coefficient divided by sum of absolute exposure coefficients, times 100. |
| n | int64 | 0 | Number of observations or participants in the documented frame. |
| model | str | 0 | Model specification in English. |

## FigS7_calibration_annotations.csv

S7 model metrics for N=3284 baseline participants with 202 deaths; horizon is 7 years. C-statistics, Brier scores and calibration coefficients are dimensionless.

Rows: 10. SHA-256: `33f22e82dabe31d15abf3703f16737132f4eb6a623171761c4e8b646e95f96ac`.

| Variable | Type | Missing | Description |
|---|---|---:|---|
| metric | str | 0 | Name of the model-performance or sample-size metric. |
| value | float64 | 0 | Numeric value of the named metric; units follow the metric name and file notes. |

## FigS7_calibration_deciles.csv

S7 calibration decile means; predicted 7-year risk and observed crude death-indicator proportion. Observed values are not censoring-adjusted Kaplan-Meier risks. n is the number of participants in each decile.

Rows: 10. SHA-256: `4a16a037ffb1eb0b1e428698c085a64e631fc1d778ae046c9c29c98e24ce8f05`.

| Variable | Type | Missing | Description |
|---|---|---:|---|
| bin | int64 | 0 | Predicted-risk decile index, from 1 to 10. |
| pred_risk | float64 | 0 | Mean predicted 7-year death risk in the decile, as a proportion. |
| obs_risk | float64 | 0 | Crude mean observed death indicator in the decile, as a proportion. |
| n | int64 | 0 | Number of observations or participants in the documented frame. |

## FigS7_calibration_smoother.csv

S7 smoothing geometry extracted directly from ggplot_build. Rates and standard errors are on a proportion scale. Missing boundary confidence limits preserve the original plot's scale handling; do not impute them.

Rows: 80. SHA-256: `0c54d0f0892deee166e3379590c9d2af6592c47bb0b00afc3e5d054d6b5f9102`.

| Variable | Type | Missing | Description |
|---|---|---:|---|
| predicted_risk | float64 | 0 | Predicted-risk coordinate for the calibration smoother, as a proportion. |
| fitted_observed_rate | float64 | 2 | Smoothed crude observed death proportion. |
| lci | float64 | 8 | Lower pointwise 95% confidence limit on the response scale. |
| uci | float64 | 0 | Upper pointwise 95% confidence limit on the response scale. |
| se | float64 | 0 | Standard error of the corresponding estimate on its response scale. |

## Figure2a_raincloud_annotations.csv

Formal group means, standard deviations, standard errors, normal-approximation 95% CIs and text positions. Outcome statistics and text_y are BAA years; x is categorical.

Rows: 8. SHA-256: `1a0eb5d713622f76c037bdb3475402eb3016748b0e8658fdc870d3bb0f4bb0d0`.

| Variable | Type | Missing | Description |
|---|---|---:|---|
| group | str | 0 | Four-category joint EVI/ozone trajectory group. |
| n | int64 | 0 | Number of observations or participants in the documented frame. |
| mean | float64 | 0 | Arithmetic mean BAA, in years. |
| sd | float64 | 0 | Sample standard deviation of BAA, in years. |
| scale | str | 0 | Residential buffer: 1KM or 500m. |
| se | float64 | 0 | Standard error of the corresponding estimate on its response scale. |
| ci_lower | float64 | 0 | Lower normal-approximation 95% confidence limit, in BAA years. |
| ci_upper | float64 | 0 | Upper normal-approximation 95% confidence limit, in BAA years. |
| x | int64 | 0 | Exposure value or categorical plotting coordinate, as specified for the file. |
| lab | str | 0 | Preformatted mean and 95% CI annotation from the formal figure. |
| text_y | float64 | 0 | Vertical BAA-year coordinate of the annotation; shifted away from zero where specified by the original source. |

## Figure2a_raincloud_boxplot_statistics.csv

Formal boxplot summaries, in BAA years. x is the plotted category centre minus 0.055. Group order at x=1,2,3,4 is low EVI/high ozone, low/low, high/high, high EVI/low ozone.

Rows: 8. SHA-256: `97fc25ce9f6dde986dfe78b84a224e79dbf289220cda472088dd619a3d7cfca2`.

| Variable | Type | Missing | Description |
|---|---|---:|---|
| x | float64 | 0 | Exposure value or categorical plotting coordinate, as specified for the file. |
| whisker_low | float64 | 0 | Lower boxplot whisker, in BAA years. |
| q1 | float64 | 0 | First quartile of BAA, in years. |
| median | float64 | 0 | Median BAA, in years. |
| q3 | float64 | 0 | Third quartile of BAA, in years. |
| whisker_high | float64 | 0 | Upper boxplot whisker, in BAA years. |
| group | str | 0 | Four-category joint EVI/ozone trajectory group. |
| scale | str | 0 | Residential buffer: 1KM or 500m. |

## Figure2a_raincloud_polygons.csv

Exact formal half-density polygons: 1024 vertices per group and buffer (512 density points plus 512 closing vertices). x is a categorical plot coordinate; baa_years is the outcome in years. Preserve row order within each polygon.

Rows: 8192. SHA-256: `411fd16189d377819ca8fe07068f4b14d0cc2bf85a43dec2c98905704eee1382`.

| Variable | Type | Missing | Description |
|---|---|---:|---|
| x | float64 | 0 | Exposure value or categorical plotting coordinate, as specified for the file. |
| baa_years | float64 | 0 | Biological-age acceleration coordinate, in years. |
| group | str | 0 | Four-category joint EVI/ozone trajectory group. |
| scale | str | 0 | Residential buffer: 1KM or 500m. |

## Figure2a_rain_histograms.csv

Public replacement for individual rain marks: full-group 2-year BAA bins from -40 to 70. Original picture used up to 600 sampled individual points per group. Counts 1-4 are suppressed, so visible bin counts need not sum to group n. All input observations fall within the displayed range.

Rows: 440. SHA-256: `7ca6a0760ade3da214228728e448e1f3a4af54a1c08eb0f486027dbe31ff4236`.

| Variable | Type | Missing | Description |
|---|---|---:|---|
| scale | str | 0 | Residential buffer: 1KM or 500m. |
| group | str | 0 | Four-category joint EVI/ozone trajectory group. |
| bin_index | int64 | 0 | One-based histogram bin index. |
| left | float64 | 0 | Left bin edge in the relevant exposure or BAA units. |
| right | float64 | 0 | Right bin edge; the last histogram bin includes its right edge. |
| center | float64 | 0 | Bin midpoint. |
| n_observations | float64 | 66 | Number of observations; missing when a bin count of 1-4 is suppressed. |
| suppressed_small_count | bool | 0 | True indicates a count between 1 and 4; its numeric value is intentionally missing. |

## Figure4_plot_metadata.csv

Primary complete-case frame N=5914; exposure quantiles and ranges used by Figure 4. EVI is unitless and ozone is in micrograms per cubic metre. Quantiles use R type 7 / NumPy linear interpolation.

Rows: 2. SHA-256: `7332c4084e339974545087fc6fdf906baacdf28b3370be28a91151e8f4f21292`.

| Variable | Type | Missing | Description |
|---|---|---:|---|
| exposure | str | 0 | Exposure name: EVI or O3. |
| n_observations | int64 | 0 | Number of observations; missing when a bin count of 1-4 is suppressed. |
| q01 | float64 | 0 | 1st exposure percentile. |
| q25 | float64 | 0 | 25th exposure percentile. |
| q50 | float64 | 0 | 50th exposure percentile. |
| q75 | float64 | 0 | 75th exposure percentile. |
| q99 | float64 | 0 | 99th exposure percentile. |
| minimum | float64 | 0 | Minimum exposure value. |
| maximum | float64 | 0 | Maximum exposure value. |
| quantile_method | str | 0 | Algorithm used to compute the exposure quantiles. |

## Figure4_exposure_histograms.csv

40 exposure bins for each Figure 4 exposure; aggregate replacement for individual rugs. Counts 1-4 are suppressed. Exposure units follow Figure4_plot_metadata.csv.

Rows: 80. SHA-256: `e134c8504cd3709357e756f8c4a6c1ebdc4ad364a4965a50e74ba2d2b4d07413`.

| Variable | Type | Missing | Description |
|---|---|---:|---|
| exposure | str | 0 | Exposure name: EVI or O3. |
| bin_index | int64 | 0 | One-based histogram bin index. |
| left | float64 | 0 | Left bin edge in the relevant exposure or BAA units. |
| right | float64 | 0 | Right bin edge; the last histogram bin includes its right edge. |
| center | float64 | 0 | Bin midpoint. |
| n_observations | float64 | 1 | Number of observations; missing when a bin count of 1-4 is suppressed. |
| suppressed_small_count | bool | 0 | True indicates a count between 1 and 4; its numeric value is intentionally missing. |

## Figure4_joint_density_bins.csv

20 by 20 exposure-bin counts, replacing sampled individual scatter. EVI is unitless; ozone is in micrograms per cubic metre. Counts 1-4 are suppressed. Existing Figure4_joint_contour_v5.csv retains the original 70 by 70 prediction grid with 490 masked cells.

Rows: 400. SHA-256: `9185cc11654b24aca9b3ebe7b2d0389fd64d0bbd8be055b331a4a96dcf4b77f0`.

| Variable | Type | Missing | Description |
|---|---|---:|---|
| evi_left | float64 | 0 | Left EVI bin edge, unitless. |
| evi_right | float64 | 0 | Right EVI bin edge, unitless. |
| evi_center | float64 | 0 | EVI bin midpoint, unitless. |
| ozone_left | float64 | 0 | Left ozone bin edge, micrograms per cubic metre. |
| ozone_right | float64 | 0 | Right ozone bin edge, micrograms per cubic metre. |
| ozone_center | float64 | 0 | Ozone bin midpoint, micrograms per cubic metre. |
| n_observations | float64 | 7 | Number of observations; missing when a bin count of 1-4 is suppressed. |
| suppressed_small_count | bool | 0 | True indicates a count between 1 and 4; its numeric value is intentionally missing. |

## Current composition: Figures 2–3 and S9

Historical source names and export provenance above are preserved. The JSON
`current_composition` mappings identify current panel use; no numeric source
table was edited for this composition.

| Figure/panel | Public source files | Current interpretation |
|---|---|---|
| 2a | Figure4_joint_density_bins.csv; Figure4_exposure_histograms.csv; Figure4_plot_metadata.csv | Joint 1-km exposure counts with marginal histograms; N=5,914 observations. Counts 1–4 remain suppressed, not zero. |
| 2b | Figure2_continuous_forest.csv | Four unweighted dual-model rows; 1 km primary, 500 m sensitivity; per-buffer IQR, robust 95% CI; N=5,914 observations/3,154 participants. |
| 2c–d | Figure2a_raincloud_polygons.csv; Figure2a_raincloud_boxplot_statistics.csv; Figure2a_raincloud_annotations.csv | Unadjusted 1-km/500-m distributions, N=6,472 per buffer. Geometry is translated by group name to the common group order. Means are points; naive mean CIs remain in the source but are not plotted. |
| 3a–b | Figure3a_adjusted_means.csv | Adjusted joint-pattern means and 95% CI, N=6,471 observations. |
| 3c–d | primary_continuous_models_v5.csv; mi_sensitivity_v5.csv | Every non-MI specification (16 rows) and both MI rows (2); no P screening. Non-MI N=5,914/3,154 participants, robust variance; MI counts unspecified, model-based variance and Rubin pooling. |
| S9a | Figure3b_contrasts.csv | All 12 contrasts; duplicated parameterisations are disclosed and are not independent evidence. |
| S9b | Figure3b_contrasts.csv; Figure3b_DiD_precision.csv | DiD and 95% CI, nominal P, approximate 80% MDE, and TOST at a ±1-year equivalence margin. |

For MI, the source CI columns are in native exposure units. Multiplication by
`mi_beta_iqr / mi_beta` converts them to the archived 1-km IQR scale, without
altering P values or refitting. The source-derived IQR is checked against the
continuous-model IQR. MI addresses remaining missingness only and does not
recover uncertainty from prior deterministic imputation. MDE is not an exclusion
bound; interaction and equivalence remain unestablished at alpha = 0.05.

