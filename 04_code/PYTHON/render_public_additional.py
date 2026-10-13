"""Render publication figures from public aggregate coordinates only."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from PIL import Image


MM = 1 / 25.4


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_figure(fig, output: Path, stem: str, dpi: int = 600) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("png", "tiff", "svg", "pdf"):
        path = output / f"{stem}.{ext}"
        kwargs = {"dpi": dpi, "bbox_inches": "tight", "facecolor": "white"}
        if ext == "tiff":
            kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
        if ext == "pdf":
            kwargs["metadata"] = {"CreationDate": None, "ModDate": None}
        if ext == "svg":
            kwargs["metadata"] = {"Date": None}
        fig.savefig(path, **kwargs)
        paths.append({"path": path.name, "sha256": sha256(path), "bytes": path.stat().st_size})
    plt.close(fig)
    return {"files": paths, "dpi": dpi}


def render_figure4(root: Path, output: Path) -> dict:
    curves_path = root / "03_source_data/Figure4_exposure_response_v5.csv"
    contour_path = root / "03_source_data/Figure4_joint_contour_v5.csv"
    curves = pd.read_csv(curves_path)
    contour = pd.read_csv(contour_path)
    required_curve = {"exposure", "x", "fit", "lin", "lci", "uci"}
    required_contour = {"EVI", "O3", "fit"}
    if not required_curve.issubset(curves.columns):
        raise ValueError("Figure 4 curve source columns are incomplete")
    if not required_contour.issubset(contour.columns):
        raise ValueError("Figure 4 contour source columns are incomplete")

    metadata = pd.read_csv(root / "03_source_data/Figure4_plot_metadata.csv").set_index("exposure")
    if set(metadata.index) != {"EVI", "O3"} or not metadata.n_observations.eq(5914).all():
        raise ValueError("Figure 4 must use the primary 5,914-observation frame")
    histograms = pd.read_csv(root / "03_source_data/Figure4_exposure_histograms.csv")
    fig, axes = plt.subplots(1, 3, figsize=(9.0, 3.55), gridspec_kw={"width_ratios": [1, 1, 1.12]})
    for axis, exposure, colour, label in zip(
        axes[:2], ["EVI", "O3"], ["#009E73", "#D55E00"], ["4-year mean EVI (1 km)", "4-year mean O₃ (μg m⁻³)"]
    ):
        d = curves[curves.exposure == exposure].sort_values("x")
        if len(d) == 0 or not (d.lci <= d.uci).all():
            raise ValueError("Invalid curve data")
        axis.axvspan(metadata.loc[exposure, "q25"], metadata.loc[exposure, "q75"], color="#D9D9D9", alpha=0.35)
        axis.fill_between(d.x, d.lci, d.uci, color=colour, alpha=0.18, linewidth=0)
        axis.plot(d.x, d.fit, color=colour, linewidth=1.6, label="Spline estimate")
        axis.plot(d.x, d.lin, color="#333333", linestyle="--", linewidth=0.9, label="Primary linear estimate")
        axis.set_xlabel(label)
        axis.set_ylabel("Predicted BAA (years)")
        axis.grid(axis="y", color="#E5E5E5", linewidth=0.5)
        axis.spines[["top", "right"]].set_visible(False)
        axis.legend(frameon=False, fontsize=6, loc="best")
        bins = histograms[(histograms.exposure == exposure) & (histograms.n_observations > 0)]
        # These are bin counts, never individual exposure coordinates.
        axis.vlines(bins.center, 0, 0.07*bins.n_observations/bins.n_observations.max(),
                    transform=axis.get_xaxis_transform(), color=colour, linewidth=0.7)

    pivot = contour.pivot(index="O3", columns="EVI", values="fit").sort_index()
    x = pivot.columns.to_numpy(dtype=float)
    y = pivot.index.to_numpy(dtype=float)
    z = pivot.to_numpy(dtype=float)
    masked = np.ma.masked_invalid(z)
    image = axes[2].pcolormesh(x, y, masked, shading="auto", cmap="RdYlGn_r")
    valid = np.isfinite(z)
    if valid.any():
        levels = np.linspace(float(np.nanmin(z)), float(np.nanmax(z)), 6)
        axes[2].contour(x, y, masked, levels=levels, colors="white", linewidths=0.45)
    axes[2].set_xlabel("4-year mean EVI (1 km)")
    axes[2].set_ylabel("4-year mean O₃ (μg m⁻³)")
    axes[2].set_title("Joint additive prediction surface", fontsize=8)
    fig.colorbar(ScalarMappable(norm=Normalize(vmin=np.nanmin(z), vmax=np.nanmax(z)), cmap="RdYlGn_r"), ax=axes[2], pad=0.02, label="Predicted BAA (years)")
    for letter, axis in zip("abc", axes):
        axis.text(-0.12, 1.04, letter, transform=axis.transAxes, fontweight="bold", fontsize=11)
    fig.suptitle("Exposure–response associations of EVI and O₃ with biological ageing", fontsize=10, fontweight="bold")
    fig.text(0.5, 0.01, "N = 5,914 observations. Shading: IQR; lower ticks: binned exposure counts. Grey cells: outside the observed data cloud.", ha="center", fontsize=6.5)
    axes[2].set_facecolor("#EEEEEE")
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    return {"sources": [curves_path.relative_to(root).as_posix(), contour_path.relative_to(root).as_posix()],
            "presentation": "Aggregate distribution ticks; no individual point layer; archived grid mask preserved.",
            "geometry_checks_pass": True,
            **save_figure(fig, output, "Figure4_exposure_response")}


def render_s1(root: Path, output: Path) -> dict:
    sys.path.insert(0, str(root / "04_code/PYTHON"))
    import compose_main_figures as composer

    canvas = Image.new("RGB", (1700, 2750), "white")
    draw = composer.ImageDraw.Draw(canvas)
    composer.draw_flow_panel(draw, root, 35, 35, 1630, 2650)
    path = output / "FigS1_participant_flow.png"
    output.mkdir(parents=True, exist_ok=True)
    composer.save(canvas, output, "FigS1_participant_flow")
    return {"source": "03_source_data/FigureS1_participant_flow.csv", "path": path.name,
            "sha256": sha256(path), "bytes": path.stat().st_size}


def render_s8(root: Path, output: Path) -> dict:
    shap_path = root / "03_source_data/FigS8_lightgbm_nested_shap_groups.csv"
    imp_path = root / "03_source_data/lightgbm_nested_importance_v5.csv"
    shap = pd.read_csv(shap_path).sort_values("share_pct", ascending=True)
    imp = pd.read_csv(imp_path).sort_values("r2_drop_mean", ascending=True).tail(10)
    fig, axes = plt.subplots(1, 2, figsize=(180 * MM, 88 * MM), gridspec_kw={"width_ratios": [1, 1.05]})
    axes[0].barh(shap.group, shap.share_pct, color="#2C86B8")
    axes[0].set_xlabel("Share of mean |SHAP| (%)")
    axes[0].set_title("a   Grouped SHAP attribution", loc="left", fontweight="bold")
    axes[0].grid(axis="x", color="#E5E5E5", linewidth=0.5)
    labels = {"climate_pc5": "Climate PC5", "climate_pc1": "Climate PC1",
              "O3_1KM_mean4y": "O₃ (1 km, 4-year mean)", "bmi": "BMI",
              "pension_insurance": "Pension insurance", "year": "Survey year", "sex": "Sex",
              "dyslipidemia": "Dyslipidaemia", "hypertension": "Hypertension",
              "alcohol_current": "Current alcohol use"}
    axes[1].barh(imp.feature.map(lambda value: labels.get(value, value)), imp.r2_drop_mean, color="#E9A923")
    axes[1].set_xlabel("Permutation importance (test-set R² drop)")
    axes[1].set_title("b   Test-set permutation importance", loc="left", fontweight="bold")
    axes[1].grid(axis="x", color="#E5E5E5", linewidth=0.5)
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(labelsize=7)
    fig.suptitle("Participant-level LightGBM evaluation", fontsize=10, fontweight="bold")
    fig.text(0.5, 0.01, "Aggregate reconstruction from held-out summary outputs; no row-level predictions or SHAP values are distributed.", ha="center", fontsize=6.5)
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    return {"sources": [shap_path.relative_to(root).as_posix(), imp_path.relative_to(root).as_posix()],
            "geometry_checks_pass": True,
            **save_figure(fig, output, "FigS8_lightgbm_nested")}


def render_s6(root: Path, output: Path) -> dict:
    path = root / "03_source_data/FigS6_gam_curve_coordinates.csv"
    data = pd.read_csv(path)
    metadata = pd.read_csv(root / "03_source_data/FigS6_gam_metadata.csv").set_index("exposure")
    weights = pd.read_csv(root / "03_source_data/FigS6_signed_weights.csv")
    legacy = pd.read_csv(root / "03_source_data/FigS6_qgcomp_metadata.csv")
    if len(legacy) != 1 or not weights.n.eq(6002).all() or int(legacy.iloc[0]['n']) != 6002:
        raise ValueError("Unexpected exploratory quantile g-computation frame")
    if (int(legacy.iloc[0].age_min), int(legacy.iloc[0].age_max)) != (45, 94):
        raise ValueError("Historical exploratory age range changed")
    if not metadata.n.eq(2785).all() or not metadata.prediction_age.eq(58.8).all() or not metadata.prediction_female.eq(0).all():
        raise ValueError("Unexpected baseline GAM prediction frame")
    if len(data) != 200 or not (data.lci <= data.uci).all():
        raise ValueError("Expected 200 dose-response coordinates with valid intervals")
    fig, axes = plt.subplots(1, 3, figsize=(9.4, 3.8), gridspec_kw={"width_ratios": [1, 1, 0.85]})
    for axis, exposure, colour in zip(axes[:2], ["EVI", "O3"], ["#009E73", "#D55E00"]):
        d = data[data.exposure == exposure].sort_values("x")
        if len(d) != 100: raise ValueError("Expected 100 coordinates per exposure")
        axis.fill_between(d.x, d.lci, d.uci, color=colour, alpha=0.18)
        line, = axis.plot(d.x, d.fit, color=colour, linewidth=1.5)
        np.testing.assert_allclose(line.get_ydata(), d.fit)
        axis.set_xlabel("EVI" if exposure == "EVI" else "O₃ (μg m⁻³)")
        axis.set_ylabel("Predicted BAA (years)")
        axis.set_title(("a  EVI" if exposure == "EVI" else "b  O₃"), loc="left", fontweight="bold")
        axis.grid(axis="y", color="#E8E8E8", linewidth=0.5)
        axis.spines[["top", "right"]].set_visible(False)
    np.testing.assert_allclose(weights.signed_weight_percent,
                               100 * weights.coefficient / weights.coefficient.abs().sum())
    axes[2].barh(weights.exposure.str.replace("O3", "O₃"), weights.signed_weight_percent,
                 color=["#009E73", "#D55E00"])
    axes[2].axvline(0, color="#777777", linewidth=0.8)
    axes[2].set(xlim=(-65, 65), xlabel="Signed relative coefficient (%)")
    axes[2].set_title("c  Quantile g-computation", loc="left", fontweight="bold")
    axes[2].spines[["top", "right"]].set_visible(False)
    for index, value in enumerate(weights.signed_weight_percent):
        axes[2].text(value / 2, index, f"{value:+.1f}%", color="white", ha="center", va="center", fontsize=8, fontweight="bold")
    fig.suptitle("Exploratory exposure–response curves and signed coefficients", fontsize=10)
    fig.text(0.5, 0.08, f"a–b: baseline GAMs, n = {int(metadata.loc['EVI', 'n']):,}; age- and sex-adjusted; evaluated at age 58.8 years, male.", ha="center", fontsize=6.5)
    fig.text(0.5, 0.035, "c: historical exploratory sample, n = 6,002, ages 45–94; distinct from the primary analysis. Shading: pointwise 95% CI.", ha="center", fontsize=6.5)
    fig.tight_layout(rect=(0, 0.14, 1, 0.94))
    return {"source": path.relative_to(root).as_posix(), "coordinates": len(data), "geometry_checks_pass": True,
            **save_figure(fig, output, "FigS6_dose_response_qgcomp")}


def render_s7(root: Path, output: Path) -> dict:
    source = root / "03_source_data"
    deciles = pd.read_csv(source / "FigS7_calibration_deciles.csv")
    smooth = pd.read_csv(source / "FigS7_calibration_smoother.csv")
    annotations = pd.read_csv(source / "FigS7_calibration_annotations.csv").set_index("metric").value
    bounded = smooth.dropna(subset=["lci", "uci"])
    if len(deciles) != 10 or len(bounded) < 2 or not (bounded.lci <= bounded.uci).all():
        raise ValueError("Invalid calibration source")
    if not deciles[["pred_risk", "obs_risk"]].ge(0).all().all():
        raise ValueError("Negative calibration probability")
    fig, axis = plt.subplots(figsize=(5.4, 4.4))
    axis.plot([0, 1], [0, 1], color="#777777", linestyle="--", linewidth=0.9, label="Perfect calibration")
    points = axis.scatter(deciles.pred_risk, deciles.obs_risk, s=30, color="#0072B2", zorder=3, label="Risk deciles")
    np.testing.assert_allclose(points.get_offsets(), deciles[["pred_risk", "obs_risk"]].to_numpy())
    axis.fill_between(smooth.predicted_risk, smooth.lci, smooth.uci, color="#0072B2", alpha=0.15)
    axis.plot(smooth.predicted_risk, smooth.fitted_observed_rate, color="#0072B2", linewidth=1.3, label="Linear calibration fit")
    upper = min(1.0, 1.1*max(deciles.pred_risk.max(), deciles.obs_risk.max()))
    axis.set(xlim=(0, upper), ylim=(0, upper), xlabel="Predicted mortality risk", ylabel="Observed mortality proportion")
    axis.set_title(f"Apparent mortality calibration\nC = {annotations['c_baa_plus_age']:.3f}; slope = {annotations['calibration_slope']:.3f}; scaled Brier = {annotations['scaled_brier_score']:.3f}", fontsize=9)
    axis.legend(frameon=False, fontsize=7, loc="upper left")
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(color="#EEEEEE", linewidth=0.5)
    from matplotlib.ticker import PercentFormatter
    axis.xaxis.set_major_formatter(PercentFormatter(1))
    axis.yaxis.set_major_formatter(PercentFormatter(1))
    fig.text(0.5, 0.035, f"N = {int(annotations['n_observations']):,}; {int(annotations['n_deaths'])} deaths; {int(annotations['prediction_horizon_years'])}-year horizon; 10 risk groups.", ha="center", fontsize=6.5)
    fig.text(0.5, 0.01, "Analytical-cohort evaluation; not external validation of the stored age score.", ha="center", fontsize=6.5)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    return {"sources": ["03_source_data/FigS7_calibration_deciles.csv", "03_source_data/FigS7_calibration_smoother.csv"],
            "deciles": 10, "geometry_checks_pass": True, **save_figure(fig, output, "FigS7_phenoage_calibration")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--main-out", type=Path, required=True)
    parser.add_argument("--si-out", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.titlesize": 8,
                         "svg.hashsalt": "bioage-public", "pdf.fonttype": 42, "ps.fonttype": 42})
    report = {
        "Figure4": render_figure4(root, args.main_out.resolve()),
        "S1": render_s1(root, args.si_out.resolve()),
        "S8": render_s8(root, args.si_out.resolve()),
        "S6": render_s6(root, args.si_out.resolve()),
        "S7": render_s7(root, args.si_out.resolve()),
    }
    (args.si_out.resolve() / "additional_public_render_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
