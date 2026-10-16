"""Export aggregate drawing coordinates from authorised local study inputs.

This author-side utility prepares source CSVs for public figure rendering.
Reviewers need only the exported CSVs, not the restricted input or R refitting.
The original R scripts define the model/geometry; no image is digitised.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
import pandas as pd
import pyreadstat

ROOT = Path(__file__).resolve().parents[2]
GROUPS = ["Low EVI / low O3", "Low EVI / high O3",
          "High EVI / low O3", "High EVI / high O3"]

R_APPEND = {
    "S6Formal": r'''
coordinates <- data.frame(exposure=dr_df$exposure,x=dr_df$x,fit=dr_df$fit,
 lci=dr_df$fit-1.96*dr_df$se,uci=dr_df$fit+1.96*dr_df$se,se=dr_df$se,
 model="Baseline single-exposure GAM; adjusted for age and sex; GCV.Cp")
write.csv(coordinates,"exports/FigS6_gam_curve_coordinates.csv",row.names=FALSE)
metadata <- data.frame(exposure=c("EVI","O3"),n=c(nobs(gam_evi),nobs(gam_o3)),
 baseline_year=2011,prediction_age=58.8,prediction_female=0,
 smooth_basis="Thin plate regression spline",smooth_k=10,
 smooth_edf=c(summary(gam_evi)$s.table[1,"edf"],summary(gam_o3)$s.table[1,"edf"]),
 estimation_method=c(gam_evi$method,gam_o3$method))
write.csv(metadata,"exports/FigS6_gam_metadata.csv",row.names=FALSE)
write.csv(weights_df[,c("exposure","weight")],"exports/FigS6_archived_display_weights.csv",row.names=FALSE)
''',
    "S6Weights": r'''
weights <- data.frame(exposure=exposure_labels,coefficient=as.numeric(beta_vec),
 se=as.numeric(se_vec),ci_lower=as.numeric(ci_low),ci_upper=as.numeric(ci_high),
 p_value=as.numeric(p_vec),signed_weight_percent=as.numeric(w_signed),
 n=nrow(dat_qg),model="Fully adjusted dual-exposure quantile g-computation; q=4")
write.csv(weights,"exports/FigS6_signed_weights.csv",row.names=FALSE)
metadata <- data.frame(n=nrow(dat_qg),age_min=min(dat_qg$age),age_max=max(dat_qg$age),
 n_outside_primary_age_range=sum(dat_qg$age<46|dat_qg$age>89),
 n_2011=sum(dat_qg$year==2011),n_2015=sum(dat_qg$year==2015),
 exposure_quantiles=4,age_filter_applied=FALSE)
write.csv(metadata,"exports/FigS6_qgcomp_metadata.csv",row.names=FALSE)
''',
    "S6": r'''
coordinates <- rbind(
 data.frame(exposure="EVI", x=dose_response_gpsm$EVI_exposure,
 fit=dose_response_gpsm$predicted_BAA, lci=dose_response_gpsm$lower,
 uci=dose_response_gpsm$upper, model="GPS weighted natural spline; df=2"),
 data.frame(exposure="O3", x=dose_response_o3$O3_exposure,
 fit=dose_response_o3$predicted_BAA, lci=dose_response_o3$lower,
 uci=dose_response_o3$upper, model="Generalized additive model; REML"))
write.csv(coordinates,"exports/FigS6_dose_response_coordinates.csv",row.names=FALSE)
write.csv(threshold_results,"exports/FigS6_threshold_verification.csv",row.names=FALSE)
write.csv(data.frame(exposure=c("EVI","O3"),n=c(nrow(df_model_evi),nrow(df_model_o3))),
 "exports/FigS6_model_frames.csv",row.names=FALSE)
''',
    "S7": r'''
write.csv(plot_df,"exports/FigS7_calibration_deciles.csv",row.names=FALSE)
built <- ggplot_build(p_cal)
smooth <- built$data[[3]][,c("x","y","ymin","ymax","se")]
names(smooth) <- c("predicted_risk","fitted_observed_rate","lci","uci","se")
write.csv(smooth,"exports/FigS7_calibration_smoother.csv",row.names=FALSE)
details <- data.frame(metric=c("prediction_horizon_years","n_observations","n_deaths",
 "c_baa_alone","c_baa_plus_age","c_age_alone","calibration_slope",
 "calibration_intercept","brier_score","scaled_brier_score"),
 value=c(predict_horizon,nrow(df_cal),sum(df_cal[[mortality_var]]),
 c_pheno_alone,c_pheno_age,c_age_alone,cal_slope,cal_intercept,brier_score,brier_scaled))
write.csv(details,"exports/FigS7_calibration_annotations.csv",row.names=FALSE)
''',
    "Figure2a": r'''
density_rows <- list(); box_rows <- list()
for(sc in c("1KM","500m")) {
 panel <- if(sc=="1KM") p_a else p_b
 built <- ggplot_build(panel)
 den <- built$data[[1]][,c("y","density","scaled","count","n","violinwidth","group")]
 den$scale <- sc; den$group <- group_levels[den$group]
 names(den)[1] <- "baa_years"
 density_rows[[sc]] <- den
 boxes <- built$data[[2]][,c("ymin","lower","middle","upper","ymax","group")]
 boxes$scale <- sc; boxes$group <- group_levels[boxes$group]
 names(boxes)[1:5] <- c("whisker_low","q1","median","q3","whisker_high")
 box_rows[[sc]] <- boxes
}
write.csv(do.call(rbind,density_rows),"exports/Figure2a_violin_coordinates.csv",row.names=FALSE)
write.csv(do.call(rbind,box_rows),"exports/Figure2a_boxplot_statistics.csv",row.names=FALSE)
write.csv(as.data.frame(stats_all),"exports/Figure2a_group_summaries.csv",row.names=FALSE)
''',
    "Raincloud": r'''
polygon_rows <- list(); box_rows <- list(); annotation_rows <- list()
for(sc in c("1KM","500m")) {
 panel <- if(sc=="1KM") p_a else p_b
 clouds <- panel$layers[[1]]$data[,c("x","y","group")]
 names(clouds)[2] <- "baa_years"
 clouds$scale <- sc
 polygon_rows[[sc]] <- clouds
 built <- ggplot_build(panel)
 boxes <- built$data[[3]][,c("x","ymin","lower","middle","upper","ymax","group")]
 boxes$scale <- sc; boxes$group <- group_levels[boxes$group]
 names(boxes)[2:6] <- c("whisker_low","q1","median","q3","whisker_high")
 box_rows[[sc]] <- boxes
 anno <- panel$layers[[5]]$data
 annotation_rows[[sc]] <- as.data.frame(anno)
}
write.csv(do.call(rbind,polygon_rows),"exports/Figure2a_raincloud_polygons.csv",row.names=FALSE)
write.csv(do.call(rbind,box_rows),"exports/Figure2a_raincloud_boxplot_statistics.csv",row.names=FALSE)
write.csv(do.call(rbind,annotation_rows),"exports/Figure2a_raincloud_annotations.csv",row.names=FALSE)
''',
}
SOURCE_SCRIPT = {"S6": "19_dose_response.R", "S7": "38_phenoage_validation.R",
                 "Figure2a": "50_figure1_real_distribution.R", "S6Weights": "15_qgcomp_final.R"}

PUBLIC_FIGURES = ["S6Formal", "S6Weights", "S7", "Raincloud", "Figure4"]
ALTERNATIVE_FILES = ["Figure2a_violin_coordinates.csv", "Figure2a_boxplot_statistics.csv",
                     "Figure2a_group_summaries.csv", "FigS6_dose_response_coordinates.csv",
                     "FigS6_model_frames.csv", "FigS6_threshold_verification.csv"]

FILE_NOTES = {
    "FigS6_gam_curve_coordinates.csv": "S6 right panels: single-exposure GAM curves at age 58.8 and female=0; baseline 2011, ages 46-89, N=2785 per model. x is EVI (unitless) or ozone (micrograms per cubic metre); fit and intervals are BAA years.",
    "FigS6_gam_metadata.csv": "S6 GAM sample and smoothing specification. The 3197 baseline observations include 412 with missing exposure, leaving 2785 in each model. GAMs adjust for age and sex, use thin plate splines with k=10 and GCV selection, and are distinct from the primary models.",
    "FigS6_archived_display_weights.csv": "Rounded signed percentages shown in the archived S6 left panel. Values are independently verified by the full-precision qgcomp output in FigS6_signed_weights.csv.",
    "FigS6_signed_weights.csv": "S6 left panel: fully adjusted dual-exposure quantile g-computation, q=4; signed coefficient divided by sum of absolute coefficients, times 100. N=6002 ages 45-94; the legacy exploratory script does not apply the primary age filter, so 88 observations lie outside ages 46-89.",
    "FigS6_qgcomp_metadata.csv": "Sample metadata for the legacy exploratory qgcomp panel; its sample differs from the primary N=5914 frame. Ages are in years and counts are observations.",
    "FigS7_calibration_deciles.csv": "S7 calibration decile means; predicted 7-year risk and observed crude death-indicator proportion. Observed values are not censoring-adjusted Kaplan-Meier risks. n is the number of participants in each decile.",
    "FigS7_calibration_smoother.csv": "S7 smoothing geometry extracted directly from ggplot_build. Rates and standard errors are on a proportion scale. Missing boundary confidence limits preserve the original plot's scale handling; do not impute them.",
    "FigS7_calibration_annotations.csv": "S7 model metrics for N=3284 baseline participants with 202 deaths; horizon is 7 years. C-statistics, Brier scores and calibration coefficients are dimensionless.",
    "Figure2a_raincloud_polygons.csv": "Exact formal half-density polygons: 1024 vertices per group and buffer (512 density points plus 512 closing vertices). x is a categorical plot coordinate; baa_years is the outcome in years. Preserve row order within each polygon.",
    "Figure2a_raincloud_boxplot_statistics.csv": "Formal boxplot summaries, in BAA years. x is the plotted category centre minus 0.055. Group order at x=1,2,3,4 is low EVI/high ozone, low/low, high/high, high EVI/low ozone.",
    "Figure2a_raincloud_annotations.csv": "Formal group means, standard deviations, standard errors, normal-approximation 95% CIs and text positions. Outcome statistics and text_y are BAA years; x is categorical.",
    "Figure2a_rain_histograms.csv": "Public replacement for individual rain marks: full-group 2-year BAA bins from -40 to 70. Original picture used up to 600 sampled individual points per group. Counts 1-4 are suppressed, so visible bin counts need not sum to group n. All input observations fall within the displayed range.",
    "Figure4_plot_metadata.csv": "Primary complete-case frame N=5914; exposure quantiles and ranges used by Figure 4. EVI is unitless and ozone is in micrograms per cubic metre. Quantiles use R type 7 / NumPy linear interpolation.",
    "Figure4_exposure_histograms.csv": "40 exposure bins for each Figure 4 exposure; aggregate replacement for individual rugs. Counts 1-4 are suppressed. Exposure units follow Figure4_plot_metadata.csv.",
    "Figure4_joint_density_bins.csv": "20 by 20 exposure-bin counts, replacing sampled individual scatter. EVI is unitless; ozone is in micrograms per cubic metre. Counts 1-4 are suppressed. Existing Figure4_joint_contour_v5.csv retains the original 70 by 70 prediction grid with 490 masked cells.",
}

COLUMN_NOTES = {
    "exposure": "Exposure name: EVI or O3.", "group": "Four-category joint EVI/ozone trajectory group.",
    "scale": "Residential buffer: 1KM or 500m.", "x": "Exposure value or categorical plotting coordinate, as specified for the file.",
    "baa_years": "Biological-age acceleration coordinate, in years.", "fit": "Model-fitted mean BAA, in years.",
    "lci": "Lower pointwise 95% confidence limit on the response scale.", "uci": "Upper pointwise 95% confidence limit on the response scale.",
    "se": "Standard error of the corresponding estimate on its response scale.",
    "model": "Model specification in English.", "n": "Number of observations or participants in the documented frame.",
    "n_observations": "Number of observations; missing when a bin count of 1-4 is suppressed.",
    "baseline_year": "Calendar year used for baseline selection.", "prediction_age": "Age fixed for model predictions, in years.",
    "prediction_female": "Female indicator fixed for predictions: 0=male, 1=female.",
    "smooth_basis": "GAM smoothing basis.", "smooth_k": "Basis dimension supplied to the smoother.",
    "smooth_edf": "Estimated effective degrees of freedom of the exposure smoother.",
    "estimation_method": "GAM smoothing-parameter selection method.",
    "weight": "Rounded signed contribution as a percentage of the sum of absolute exposure coefficients.",
    "signed_weight_percent": "Full-precision signed coefficient divided by sum of absolute exposure coefficients, times 100.",
    "coefficient": "Change in BAA years per one exposure-quartile score increase.",
    "ci_lower": "Lower normal-approximation 95% confidence limit, in BAA years.",
    "ci_upper": "Upper normal-approximation 95% confidence limit, in BAA years.",
    "p_value": "Two-sided model p-value.", "age_min": "Minimum observed age, in years.",
    "age_max": "Maximum observed age, in years.", "n_outside_primary_age_range": "Number of observations with age outside 46-89 years.",
    "n_2011": "Number of observations in the 2011 wave.", "n_2015": "Number of observations in the 2015 wave.",
    "exposure_quantiles": "Number of exposure quantile categories.", "age_filter_applied": "Whether the primary age-range filter was applied to this exploratory model.",
    "bin": "Predicted-risk decile index, from 1 to 10.", "pred_risk": "Mean predicted 7-year death risk in the decile, as a proportion.",
    "obs_risk": "Crude mean observed death indicator in the decile, as a proportion.",
    "predicted_risk": "Predicted-risk coordinate for the calibration smoother, as a proportion.",
    "fitted_observed_rate": "Smoothed crude observed death proportion.",
    "metric": "Name of the model-performance or sample-size metric.", "value": "Numeric value of the named metric; units follow the metric name and file notes.",
    "whisker_low": "Lower boxplot whisker, in BAA years.", "whisker_high": "Upper boxplot whisker, in BAA years.",
    "q1": "First quartile of BAA, in years.", "q3": "Third quartile of BAA, in years.",
    "median": "Median BAA, in years.", "mean": "Arithmetic mean BAA, in years.",
    "sd": "Sample standard deviation of BAA, in years.", "lab": "Preformatted mean and 95% CI annotation from the formal figure.",
    "text_y": "Vertical BAA-year coordinate of the annotation; shifted away from zero where specified by the original source.",
    "bin_index": "One-based histogram bin index.", "left": "Left bin edge in the relevant exposure or BAA units.",
    "right": "Right bin edge; the last histogram bin includes its right edge.", "center": "Bin midpoint.",
    "suppressed_small_count": "True indicates a count between 1 and 4; its numeric value is intentionally missing.",
    "q01": "1st exposure percentile.", "q25": "25th exposure percentile.", "q50": "50th exposure percentile.",
    "q75": "75th exposure percentile.", "q99": "99th exposure percentile.",
    "minimum": "Minimum exposure value.", "maximum": "Maximum exposure value.",
    "quantile_method": "Algorithm used to compute the exposure quantiles.",
    "evi_left": "Left EVI bin edge, unitless.", "evi_right": "Right EVI bin edge, unitless.",
    "evi_center": "EVI bin midpoint, unitless.", "ozone_left": "Left ozone bin edge, micrograms per cubic metre.",
    "ozone_right": "Right ozone bin edge, micrograms per cubic metre.", "ozone_center": "Ozone bin midpoint, micrograms per cubic metre.",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def export_r(figure: str, root: Path, scratch: Path, rscript: str,
             analysis: Path, mortality: Path | None,
             raincloud_source: Path | None = None,
             s6_source: Path | None = None, s6_input: Path | None = None) -> dict:
    work = scratch / figure
    for folder in ["CODE", "DATA/analysis_ready", "output/tables", "output/figures", "logs", "exports"]:
        (work / folder).mkdir(parents=True, exist_ok=True)
    marker = work / "exports/export_complete.txt"
    marker.unlink(missing_ok=True)
    if figure == "S6Formal":
        if not s6_source or not s6_source.is_file() or not s6_input or not s6_input.is_file():
            raise ValueError("Formal S6 export requires --s6-source and --s6-input")
        script = s6_source
        content = script.read_text(encoding="utf-8-sig")
        header = content[:content.index("# =====================================================================\n# FIGURE 1:")]
        section = content[content.index("# FIGURE 4: Dose-response"):content.index("# FIGURE 5:")]
        adapted = header + section
        adapted, n_read = re.subn(r'df <- read.csv\("[^\"]+"\)',
                                 'df <- read.csv("DATA/analysis_ready/analysis_dataset.csv")', adapted)
        adapted, n_save = re.subn(r'save_tiff\(fig4, "[^\"]+", 183, 100\)',
                                 'save_tiff(fig4, "output/figures/S6_formal.tiff", 183, 100)', adapted)
        if (n_read, n_save) != (1, 1):
            raise ValueError("Unexpected S6 source IO layout")
        shutil.copy2(s6_input, work / "DATA/analysis_ready/analysis_dataset.csv")
        (work / "CODE" / script.name).write_text(adapted, encoding="utf-8")
    elif figure == "Raincloud":
        if raincloud_source is None or not raincloud_source.is_file():
            raise ValueError("Raincloud export requires --raincloud-source pointing to source script 71")
        script = raincloud_source
        adapted = script.read_text(encoding="utf-8-sig")
        for variable, value in [("B01", "."), ("DTA", "DATA/analysis_ready/analysis_dataset.dta")]:
            adapted, matches = re.subn(r"^" + variable + r" <- .*$",
                                       f'{variable} <- "{value}"', adapted, flags=re.MULTILINE)
            if matches != 1:
                raise ValueError(f"Expected one {variable} IO assignment, found {matches}")
        for old, new in [("08_\u5206\u6790\u811a\u672c", "logs"),
                         ("09_\u5206\u6790\u7ed3\u679c", "output/tables"),
                         ("10_\u56fe\u8868\u4ea7\u51fa", "output/figures")]:
            adapted = adapted.replace('"' + old + '"', '"' + new + '"')
        (work / "CODE" / script.name).write_text(adapted, encoding="utf-8")
    else:
        script = root / "04_code/CODE" / SOURCE_SCRIPT[figure]
        shutil.copy2(script, work / "CODE" / script.name)
    shutil.copy2(analysis, work / "DATA/analysis_ready/analysis_dataset.dta")
    if figure == "S7":
        if mortality is None or not mortality.is_file():
            return {"status": "missing_mortality_input"}
        shutil.copy2(mortality, work / "DATA/analysis_ready/analysis_dataset_mortality_baseline.dta")
    wrapper = 'source("CODE/' + script.name + '", encoding="UTF-8")\n'
    wrapper += R_APPEND[figure]
    wrapper += '\nwriteLines(capture.output(sessionInfo()), "exports/r_session_info.txt")\n'
    wrapper += '\nwriteLines("EXPORT_COMPLETE", "exports/export_complete.txt")\n'
    (work / "export_coordinates.R").write_text(wrapper, encoding="utf-8")
    env = os.environ.copy()
    env.update(LANG="English_United States.utf8", LC_ALL="English_United States.utf8")
    with (work / "execution.log").open("w", encoding="utf-8") as log:
        run = subprocess.run([rscript, "--vanilla", "export_coordinates.R"], cwd=work,
                             env=env, stdout=log, stderr=subprocess.STDOUT)
    result = {"script": script.name, "script_sha256": sha(script),
              "executed_script_sha256": sha(work / "CODE" / script.name),
              "process_returncode": run.returncode,
              "export_complete": (work / "exports/export_complete.txt").is_file(),
              "status": "process_ok" if run.returncode == 0 else "process_failed"}
    if not result["export_complete"]:
        return result
    if figure == "S6Formal":
        from PIL import Image
        actual_image = np.asarray(Image.open(work / "output/figures/S6_formal.tiff").convert("RGB"))
        formal_image = np.asarray(Image.open(root / "02_figures/supplementary/FigS6_dose_response_qgcomp.png").convert("RGB"))
        result["render_pixels_identical"] = actual_image.shape == formal_image.shape and np.array_equal(actual_image, formal_image)
        result["source_csv_sha256"] = sha(s6_input)
        result["io_adaptation"] = "Executed original header and independent Figure 4 section only; redirected CSV input and TIFF output to scratch."
        result["model_specification"] = "Baseline 2011, age 46-89, single-exposure mgcv GAM adjusted for age and female indicator; prediction at age 58.8 and female=0."
        if not result["render_pixels_identical"]:
            result["render_mean_absolute_pixel_difference"] = float(np.abs(actual_image.astype(float) - formal_image).mean()) if actual_image.shape == formal_image.shape else None
    elif figure == "S6Weights":
        weights = pd.read_csv(work / "exports/FigS6_signed_weights.csv")
        np.testing.assert_allclose(np.round(weights.signed_weight_percent, 1), [-49.6, 50.4], rtol=0, atol=0.000001)
        np.testing.assert_allclose(np.round(weights.coefficient, 3), [-0.531, 0.539], rtol=0, atol=0.000001)
        result["formal_weight_percentages_match"] = True
        result["archived_qgcomp_coefficients_match"] = True
    elif figure == "S6":
        actual = pd.read_csv(work / "exports/FigS6_threshold_verification.csv")
        expected = pd.read_csv(root / "06_tables/threshold_results_EVI_O3.csv")
        pd.testing.assert_frame_equal(actual, expected, check_dtype=False, rtol=1e-8, atol=1e-10)
        result["archived_thresholds_match"] = True
    elif figure == "S7":
        actual = pd.read_csv(work / "exports/FigS7_calibration_annotations.csv").set_index("metric").value
        expected = pd.read_csv(root / "03_source_data/PhenoAge_validation_metrics.csv").set_index("Metric").Value
        metric_map = {"prediction_horizon_years": "Prediction horizon (years)",
                      "n_observations": "N (baseline)", "n_deaths": "N deaths",
                      "c_baa_alone": "C-statistic (BAA_PhenoAge alone)",
                      "c_baa_plus_age": "C-statistic (BAA_PhenoAge + age)",
                      "c_age_alone": "C-statistic (age alone)",
                      "calibration_slope": "Calibration slope",
                      "calibration_intercept": "Calibration intercept",
                      "brier_score": "Brier score", "scaled_brier_score": "Scaled Brier score"}
        for new, old in metric_map.items():
            np.testing.assert_allclose(actual[new], expected[old], rtol=0, atol=0.000051)
        result["metrics_match_archived_rounding"] = True
        result["render_sha256"] = sha(work / "output/figures/phenoage_calibration_plot.png")
        result["formal_render_sha256"] = sha(root / "02_figures/supplementary/FigS7_phenoage_calibration.png")
        result["render_byte_identical"] = result["render_sha256"] == result["formal_render_sha256"]
    elif figure == "Figure2a":
        actual = pd.read_csv(work / "exports/Figure2a_group_summaries.csv")
        expected = pd.read_csv(root / "03_source_data/Figure1_observed_distribution.csv")
        expected = expected.set_index(["scale", "group"])
        actual = actual.set_index(["scale", "group"]).loc[expected.index]
        np.testing.assert_allclose(actual.n, expected.n, rtol=0, atol=0)
        for column in expected.columns.drop("n"):
            np.testing.assert_allclose(actual[column], expected[column], rtol=0, atol=0.000501)
        result["group_summaries_match_archived_rounding"] = True
        image = work / "output/figures/Figure1_observed_distribution.png"
        formal = root / "02_figures/supplementary/panel_sources/Figure1_raincloud.png"
        result["render_byte_identical"] = sha(image) == sha(formal)
        result["use_note"] = "Full-violin geometry is superseded by Raincloud source 71 for the formal panel."
    elif figure == "Raincloud":
        result["io_only_adaptation"] = "Reassigned B01 and DTA; redirected three output directories to scratch."
        actual = pd.read_csv(work / "exports/Figure2a_raincloud_annotations.csv")
        expected = pd.read_csv(root / "03_source_data/Figure1_observed_distribution.csv")
        expected = expected.set_index(["scale", "group"])
        actual = actual.set_index(["scale", "group"]).loc[expected.index]
        np.testing.assert_allclose(actual.n, expected.n, rtol=0, atol=0)
        for column in ["mean", "sd", "se", "ci_lower", "ci_upper"]:
            np.testing.assert_allclose(actual[column], expected[column], rtol=0, atol=0.000501)
        result["group_summaries_match_archived_rounding"] = True
        image = work / "output/figures/Figure1_raincloud_final_v3.png"
        formal = root / "02_figures/supplementary/panel_sources/Figure1_raincloud.png"
        result["render_sha256"] = sha(image)
        result["formal_render_sha256"] = sha(formal)
        result["render_byte_identical"] = sha(image) == sha(formal)
        if not result["render_byte_identical"]:
            raise ValueError("Raincloud source 71 did not reproduce the formal panel byte-identically")
        result["rain_histogram"] = export_rain_histogram(analysis, work / "exports")
    result["exports"] = []
    for path in (work / "exports").glob("*.csv"):
        frame = pd.read_csv(path)
        if {"id", "id_num", "participant_id"}.intersection(frame.columns):
            raise ValueError("Unexpected identifier in aggregate output")
        target = root / "03_source_data" / path.name
        shutil.copy2(path, target)
        result["exports"].append({"file": path.name, "rows": len(frame), "sha256": sha(path)})
    if run.returncode:
        result["status"] = "coordinates_exported_and_checked_but_process_exit_failed"
    return result


def export_rain_histogram(analysis: Path, destination: Path) -> dict:
    """Replace original individual rain marks by full-group aggregate bins."""
    frame, _ = pyreadstat.read_dta(str(analysis), apply_value_formats=False)
    frame = frame.loc[frame.age.between(46, 89)].copy()
    rows, checks = [], []
    edges = np.linspace(-40, 70, 56)
    expected = {"1KM": [2273, 1584, 1555, 1060], "500m": [2720, 2178, 1108, 466]}
    for scale in ["1KM", "500m"]:
        for group_id, group in enumerate(GROUPS, start=1):
            values = frame.loc[frame["joint_exp_" + scale] == group_id, "BAA_PhenoAge"].dropna()
            if len(values) != expected[scale][group_id - 1]:
                raise ValueError(f"Unexpected raincloud group size: {scale}, {group}, {len(values)}")
            counts, _ = np.histogram(values, bins=edges)
            checks.append(dict(scale=scale, group=group, n=len(values),
                               within_plot_range=int(counts.sum()),
                               below_plot_range=int((values < edges[0]).sum()),
                               above_plot_range=int((values > edges[-1]).sum())))
            for i, count in enumerate(counts):
                rows.append(dict(scale=scale, group=group, bin_index=i + 1,
                                 left=edges[i], right=edges[i + 1], center=(edges[i] + edges[i + 1]) / 2,
                                 n_observations=int(count) if count == 0 or count >= 5 else np.nan,
                                 suppressed_small_count=bool(0 < count < 5)))
    pd.DataFrame(rows).to_csv(destination / "Figure2a_rain_histograms.csv", index=False)
    return {"display_change": "Full-group aggregate counts replace the original random sample of up to 600 individual rain marks per group; 2-year bins; counts 1-4 suppressed.",
            "checks": checks}


def export_figure4(root: Path, analysis: Path) -> dict:
    frame, _ = pyreadstat.read_dta(str(analysis), apply_value_formats=False)
    frame = frame.loc[frame.age.between(46, 89)].copy()
    frame["ln_hh_consumption_total"] = np.log(frame.hh_consumption_total.where(frame.hh_consumption_total > 0))
    columns = ["BAA_PhenoAge", "EVI_1KM_mean4y", "O3_1KM_mean4y", "id_num", "age", "sex",
               "education", "ln_hh_consumption_total", "retirement", "pension_insurance",
               "smoking_current", "alcohol_current", "TMP", "PER", "SNOW", "DEM", "lat",
               "bmi", "hypertension", "dyslipidemia", "year"]
    model = frame.dropna(subset=columns)
    if len(model) != 5914:
        raise ValueError(f"Figure 4 primary frame changed: {len(model)}")
    curves = pd.read_csv(root / "03_source_data/Figure4_exposure_response_v5.csv")
    contour = pd.read_csv(root / "03_source_data/Figure4_joint_contour_v5.csv")
    if len(contour) != 4900 or contour.fit.isna().sum() != 490:
        raise ValueError("Archived masked Figure 4 contour differs from 70x70/490-mask contract")
    metadata, bins = [], []
    arrays = []
    for exposure, variable in [("EVI", "EVI_1KM_mean4y"), ("O3", "O3_1KM_mean4y")]:
        values = model[variable].to_numpy()
        arrays.append(values)
        q = np.quantile(values, [.01, .25, .5, .75, .99], method="linear")
        plotted = curves.loc[curves.exposure == exposure, "x"]
        np.testing.assert_allclose([plotted.min(), plotted.max()], q[[0, 4]], rtol=1e-12)
        metadata.append(dict(exposure=exposure, n_observations=len(values),
                             q01=q[0], q25=q[1], q50=q[2], q75=q[3], q99=q[4],
                             minimum=values.min(), maximum=values.max(),
                             quantile_method="R type 7 / NumPy linear"))
        count, edges = np.histogram(values, bins=40)
        for i, n in enumerate(count):
            bins.append(dict(exposure=exposure, bin_index=i+1, left=edges[i], right=edges[i+1],
                             center=(edges[i]+edges[i+1])/2,
                             n_observations=int(n) if n == 0 or n >= 5 else np.nan,
                             suppressed_small_count=bool(0 < n < 5)))
    count2d, xedges, yedges = np.histogram2d(arrays[0], arrays[1], bins=20)
    joint = []
    for i in range(len(xedges)-1):
        for j in range(len(yedges)-1):
            n = int(count2d[i, j])
            joint.append(dict(evi_left=xedges[i], evi_right=xedges[i+1],
                              evi_center=(xedges[i]+xedges[i+1])/2,
                              ozone_left=yedges[j], ozone_right=yedges[j+1],
                              ozone_center=(yedges[j]+yedges[j+1])/2,
                              n_observations=n if n == 0 or n >= 5 else np.nan,
                              suppressed_small_count=bool(0 < n < 5)))
    exports = []
    for name, rows in [("Figure4_plot_metadata.csv", metadata),
                       ("Figure4_exposure_histograms.csv", bins),
                       ("Figure4_joint_density_bins.csv", joint)]:
        path = root / "03_source_data" / name
        pd.DataFrame(rows).to_csv(path, index=False, float_format="%.17g")
        exports.append({"file": name, "rows": len(rows), "sha256": sha(path)})
    return {"status": "verified_aggregate_layers_exported", "n_observations": len(model),
            "curve_endpoint_quantiles_match": True, "contour_grid_rows": len(contour),
            "masked_contour_cells": int(contour.fit.isna().sum()),
            "contour_sha256": sha(root / "03_source_data/Figure4_joint_contour_v5.csv"),
            "display_change": "Individual rugs and sampled points replaced by binned counts; counts 1-4 suppressed.",
            "exports": exports}


def finalize_public_metadata(root: Path, scratch: Path, report: dict) -> None:
    """Publish aggregate-only metadata, and archive alternative drawing inputs."""
    schemas, provenance_figures = [], {}
    for figure in PUBLIC_FIGURES:
        result = report["figures"][figure]
        if "exports" not in result:
            raise ValueError(f"No completed coordinate export for {figure}")
        provenance_figures[figure] = {key: result[key] for key in [
            "script", "script_sha256", "executed_script_sha256", "process_returncode",
            "export_complete", "status", "render_pixels_identical", "render_byte_identical",
            "render_sha256", "formal_render_sha256", "source_csv_sha256",
            "formal_weight_percentages_match", "archived_qgcomp_coefficients_match",
            "group_summaries_match_archived_rounding", "metrics_match_archived_rounding",
            "curve_endpoint_quantiles_match", "n_observations", "contour_grid_rows",
            "masked_contour_cells", "contour_sha256", "display_change", "model_specification",
        ] if key in result}
        provenance_figures[figure]["exports"] = result["exports"]
        for item in result["exports"]:
            path = root / "03_source_data" / item["file"]
            frame = pd.read_csv(path)
            if sha(path) != item["sha256"] or len(frame) != item["rows"]:
                raise ValueError(f"Export changed since verification: {path.name}")
            if any("\u4e00" <= char <= "\u9fff" for char in path.read_text(encoding="utf-8-sig")):
                raise ValueError(f"Non-English source-data content: {path.name}")
            if {"id", "id_num", "participant_id", "lon", "lat", "address"}.intersection(frame.columns):
                raise ValueError(f"Unexpected individual variable: {path.name}")
            if "suppressed_small_count" in frame:
                suppressed = frame.suppressed_small_count
                if not frame.loc[suppressed, "n_observations"].isna().all():
                    raise ValueError(f"Small-count suppression changed: {path.name}")
                counts = frame.loc[~suppressed, "n_observations"]
                if not ((counts == 0) | (counts >= 5)).all():
                    raise ValueError(f"Unmasked small count: {path.name}")
            schemas.append({"figure": figure, "file": path.name, "rows": len(frame),
                            "description": FILE_NOTES[path.name], "sha256": sha(path),
                            "record_type": "aggregate_drawing_coordinates_or_statistics",
                            "columns": [{"name": column, "dtype": str(frame[column].dtype),
                                         "missing": int(frame[column].isna().sum()),
                                         "description": COLUMN_NOTES[column]} for column in frame]})
    source = (root / "03_source_data").resolve()
    archive = (scratch / "alternative_sources_not_for_public_rendering").resolve()
    if not archive.is_relative_to(root.resolve()) or not source.is_relative_to(root.resolve()):
        raise ValueError("Archive and source paths must stay inside the project workspace")
    archive.mkdir(parents=True, exist_ok=True)
    for filename in ALTERNATIVE_FILES:
        old = source / filename
        if old.is_file():
            target = archive / filename
            if target.exists():
                if sha(old) != sha(target):
                    raise ValueError(f"Archive collision for {filename}")
                old.unlink()
            else:
                shutil.move(str(old), str(target))
    provenance = {
        "schema_version": "1.0", "purpose": "Public redraw of formal statistical geometry from aggregate CSVs.",
        "export_utility": "04_code/PYTHON/export_public_figure_sources.py",
        "data_scope": "No participant identifiers, individual observations, individual predictions or individual SHAP values are included in these exports.",
        "analysis_input_sha256": report["analysis_input_sha256"],
        "mortality_input_sha256": report.get("mortality_input_sha256"),
        "figures": provenance_figures,
        "important_notes": [
            "S6 GAM curves use 2785 baseline observations, age/sex adjustment and predictions for a male aged 58.8 years. They are exploratory, not the primary multivariable models.",
            "S6 qgcomp weights use 6002 observations aged 45-94; the original exploratory model includes 88 observations outside the primary age range. Full-precision coefficients independently reproduce the displayed -49.6%/+50.4% weights.",
            "S7 observed calibration proportions are crude death-indicator means, not censoring-adjusted Kaplan-Meier risk estimates.",
            "Raincloud individual dots and Figure 4 rugs/scatter are replaced with aggregate binned counts. Counts 1-4 are suppressed. Statistical curves, boxplots, means and confidence intervals preserve formal source geometry.",
            "Author-side R processes produced EXPORT_COMPLETE and validated exports, then returned Windows 3221225477 (0xc0000005). This is not a clean R process exit. Public redraw reads CSVs and does not invoke these R fits.",
        ],
    }
    (source / "figure_source_schema.json").write_text(json.dumps(schemas, indent=2), encoding="utf-8")
    (source / "figure_source_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    lines = ["# Aggregate figure-source codebook", "",
             "These additional source files enable public figure rendering using only `03_source_data` and `04_code`.",
             "They contain aggregated drawing geometry and statistics, with no individual records. Missing histogram counts flagged as suppressed represent counts of 1-4. Other missing values preserve original plot boundaries or masks.", "",
             "The existing Figure 4 curve and masked-contour CSVs remain authoritative. This codebook describes the newly exported supporting layers, S6/S7 geometry and the formal raincloud geometry.", "",
             "## Provenance and interpretation", ""]
    lines.extend("- " + note for note in provenance["important_notes"])
    for schema in schemas:
        lines.extend(["", "## " + schema["file"], "", schema["description"], "",
                      f"Rows: {schema['rows']}. SHA-256: `{schema['sha256']}`.", "",
                      "| Variable | Type | Missing | Description |", "|---|---|---:|---|"])
        lines.extend(f"| {col['name']} | {col['dtype']} | {col['missing']} | {col['description']} |" for col in schema["columns"])
    (source / "FIGURE_SOURCE_CODEBOOK.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--analysis-input", type=Path)
    parser.add_argument("--mortality-input", type=Path)
    parser.add_argument("--raincloud-source", type=Path)
    parser.add_argument("--s6-source", type=Path)
    parser.add_argument("--s6-input", type=Path)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--finalize-existing", action="store_true",
                        help="Validate existing exports, create public metadata and archive alternative sources; no R refitting.")
    parser.add_argument("--figures", nargs="+", choices=["S6", "S6Formal", "S6Weights", "S7", "Figure2a", "Figure4", "Raincloud"],
                        default=["S6Formal", "S6Weights", "S7", "Raincloud", "Figure4"],
                        help="Figure2a and S6 refer to alternative source scripts, not the formal archived panels.")
    args = parser.parse_args()
    root, scratch = args.root.resolve(), args.work_dir.resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    report_path = scratch / "figure_source_export_report.json"
    if args.finalize_existing:
        report = json.loads(report_path.read_text(encoding="utf-8"))
        finalize_public_metadata(root, scratch, report)
        print("Public source metadata validated and written; alternative sources archived.")
        return
    analysis = (args.analysis_input or root / "04_code/DATA/analysis_ready/analysis_dataset.dta").resolve()
    mortality = args.mortality_input.resolve() if args.mortality_input else None
    report = json.loads(report_path.read_text(encoding="utf-8")) if report_path.exists() else {"figures": {}}
    report["analysis_input_sha256"] = sha(analysis)
    if mortality and mortality.is_file():
        report["mortality_input_sha256"] = sha(mortality)
    for figure in args.figures:
        print(f"Exporting {figure}", flush=True)
        try:
            result = export_figure4(root, analysis) if figure == "Figure4" else export_r(
                figure, root, scratch, args.rscript, analysis, mortality,
                args.raincloud_source.resolve() if args.raincloud_source else None,
                args.s6_source.resolve() if args.s6_source else None,
                args.s6_input.resolve() if args.s6_input else None)
        except Exception as error:
            result = {"status": "failed", "error": str(error)}
        report["figures"][figure] = result
        report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps({"figure": figure, "status": result["status"]}), flush=True)
    if any(x.get("status") in {"failed", "process_failed", "missing_mortality_input"}
           for x in report["figures"].values()):
        raise SystemExit(1)
    if all(figure in report["figures"] for figure in PUBLIC_FIGURES):
        finalize_public_metadata(root, scratch, report)


if __name__ == "__main__":
    main()
