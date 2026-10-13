"""Render source-verified Supplementary Figures S2-S4 as native graphics.

Run from any working directory. The package root is inferred from this script.
An explicit output directory prevents accidental replacement of formal figures.
S4 requires public aggregate box statistics. Such reconstructions omit individual
outlier points and explicitly state this difference from the full-data figure.
S5 is rendered separately by render_public_s5.py from the retained S13/S14
summaries. No model is fitted and no individual value is fabricated.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.text import Text
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
from PIL import Image


MM = 1 / 25.4
COLORS = ["#367A9A", "#CA9054", "#769B70", "#9B81B0"]
AGE_GROUPS = ["<=60 years", "61-70 years", "71-80 years", ">80 years", "All observations"]
PATTERNS = ["both_improved", "only_o3_improved"]
PATTERN_LABELS = ["EVI and ozone shifted", "Ozone only shifted"]
ADJUSTMENTS = ["Unadjusted", "Age and sex", "Multivariable"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_info(path: Path, root: Path) -> dict:
    return {"path": path.relative_to(root).as_posix(), "sha256": sha256(path)}


def style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.65, "legend.frameon": False,
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "savefig.facecolor": "white", "figure.facecolor": "white",
        "text.usetex": False,
    })


def load_coefficients(path: Path, expected_rows: int) -> pd.DataFrame:
    data = pd.read_csv(path, dtype=str, keep_default_na=False)
    if len(data) != expected_rows:
        raise ValueError(f"Expected {expected_rows} source rows, found {len(data)}")
    coefficient = data["Coefficient (SE)"].str.extract(r"^(-?\d+\.\d+) \((\d+\.\d+)\)$")
    intervals = data["95% CI"].str.extract(r"^\[(-?\d+\.\d+), (-?\d+\.\d+)\]$")
    if coefficient.isna().any().any() or intervals.isna().any().any():
        raise ValueError("Unrecognised coefficient or CI format")
    data["estimate"] = coefficient[0].astype(float)
    data["standard_error"] = coefficient[1].astype(float)
    data["ci_lower"] = intervals[0].astype(float)
    data["ci_upper"] = intervals[1].astype(float)
    if not ((data.ci_lower <= data.estimate) & (data.estimate <= data.ci_upper)).all():
        raise ValueError("A source estimate lies outside its confidence interval")
    if set(data.Buffer) != {"1KM"}:
        raise ValueError("Unexpected exposure buffer")
    return data


def style_axis(axis) -> None:
    axis.set_axisbelow(True)
    axis.grid(axis="y", color="#E3E3E3", linewidth=0.55)
    axis.axhline(0, color="#333333", linestyle=(0, (4, 3)), linewidth=0.8)
    axis.tick_params(length=3, width=0.65)


def export(fig, stem: str, output: Path, dpi: int) -> dict:
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    canvas = fig.bbox
    clipped = []
    for item in fig.findobj(match=Text):
        if not item.get_visible() or not item.get_text():
            continue
        bounds = item.get_window_extent(renderer)
        if bounds.width and bounds.height and (
            bounds.x0 < -1 or bounds.y0 < -1
            or bounds.x1 > canvas.x1 + 1 or bounds.y1 > canvas.y1 + 1
        ):
            clipped.append(item.get_text())
    if clipped:
        raise ValueError(f"Text outside figure canvas: {clipped}")
    exports = []
    for extension in ["png", "tiff", "svg", "pdf"]:
        path = output / f"{stem}.{extension}"
        options = {"dpi": dpi}
        if extension == "tiff":
            options["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(path, **options)
        record = {"filename": path.name, "sha256": sha256(path), "bytes": path.stat().st_size}
        if extension in ["png", "tiff"]:
            with Image.open(path) as image:
                record.update({"pixels": list(image.size), "dpi": [float(value) for value in image.info.get("dpi", [])]})
                expected = np.floor(fig.get_size_inches() * dpi).astype(int)
                if np.max(np.abs(np.array(image.size) - expected)) > 1:
                    raise AssertionError("Unexpected raster dimensions")
        if extension == "svg":
            svg = path.read_text(encoding="utf-8")
            if "<image" in svg or "<text" not in svg:
                raise AssertionError("SVG must contain editable text and no raster image")
            record.update({"editable_text": True, "embedded_rasters": 0})
        exports.append(record)
    plt.close(fig)
    return {"width_mm": 180, "height_mm": round(fig.get_size_inches()[1] / MM),
            "text_within_canvas": True, "exports": exports}


def plot_bar(axis, x: float, row: pd.Series, width: float, color: str) -> None:
    rectangle = axis.bar(x, row.estimate, width=width, color=color, edgecolor="white", linewidth=0.5)[0]
    error = axis.errorbar(x, row.estimate,
                         yerr=[[row.estimate - row.ci_lower], [row.ci_upper - row.estimate]],
                         fmt="none", ecolor="#222222", elinewidth=0.8, capsize=2.5, capthick=0.8)
    assert rectangle.get_height() == row.estimate
    vertices = error.lines[2][0].get_segments()[0]
    np.testing.assert_allclose(vertices[:, 1], [row.ci_lower, row.ci_upper], rtol=0, atol=1e-15)


def render_s2(root: Path, output: Path, dpi: int) -> dict:
    path = root / "03_source_data/quartile_scores_detailed_results.csv"
    data = load_coefficients(path, 15)
    data.to_csv(output / "FigS2_verified_plot_values.csv", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(180 * MM, 108 * MM), sharey=True)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.74, bottom=0.27, wspace=0.17)
    fig.suptitle("Quartile-score associations with biological age acceleration", fontsize=10, y=0.97)
    panel_specs = [("EVI_score", "EVI score"), ("mixed_score", "Mixed score"), ("O3_score", "Ozone score")]
    plotted = 0
    for label, axis, (exposure, title) in zip("abc", axes, panel_specs):
        for index, adjustment in enumerate(ADJUSTMENTS):
            rows = data[(data.Exposure == exposure) & (data.Adjustment == adjustment)]
            for _, row in rows.iterrows():
                dual = row["Model Type"] == "Dual-scores"
                x = index if exposure == "mixed_score" else index + (-0.19 if dual else 0.19)
                plot_bar(axis, x, row, 0.72 if exposure == "mixed_score" else 0.36, COLORS[1 if dual else 0])
                plotted += 1
        axis.set_title(f"{label}   {title}", loc="left", pad=9)
        axis.set_xticks(range(3), ["Unadjusted", "Age and sex", "Multivariable"], rotation=36, ha="right")
        axis.set_xlim(-0.55, 2.55)
        axis.set_ylim(-1.08, 0.07)
        axis.yaxis.set_major_locator(MultipleLocator(0.25))
        style_axis(axis)
    axes[0].set_ylabel("BAA difference per score point (years)")
    fig.legend(handles=[Patch(facecolor=COLORS[0], label="Single-score model"),
                        Patch(facecolor=COLORS[1], label="Dual-score model")],
               loc="upper center", bbox_to_anchor=(0.5, 0.90), ncol=2)
    fig.text(0.10, 0.115, "Higher scores indicate higher greenness or lower ozone; mixed score = EVI score + ozone score.", fontsize=7)
    fig.text(0.10, 0.064, "Bars: coefficients; error bars: 95% confidence intervals. Exposure buffer: 1 km.", fontsize=7)
    assert plotted == 15
    return {"source": source_info(path, root), "plot_rows": plotted, "precision_decimals": 3,
            "geometry_equals_source": True, **export(fig, "FigS2_quartile_scores", output, dpi)}


def render_s3(root: Path, output: Path, dpi: int) -> dict:
    path = root / "03_source_data/iqr_models_detailed_results.csv"
    data = load_coefficients(path, 8)
    data.to_csv(output / "FigS3_verified_plot_values.csv", index=False)
    fig, axis = plt.subplots(figsize=(180 * MM, 125 * MM))
    fig.subplots_adjust(left=0.13, right=0.98, top=0.85, bottom=0.28)
    fig.suptitle("Continuous-exposure associations with biological age acceleration", fontsize=10, y=0.97)
    plotted = 0
    for x, exposure in enumerate(["EVI_iqr", "O3_iqr"]):
        for j in range(4):
            adjustment = ADJUSTMENTS[min(j, 2)]
            model = "Dual-exposure" if j == 3 else "Single-exposure"
            rows = data[(data.Exposure == exposure) & (data.Adjustment == adjustment) & (data["Model Type"] == model)]
            if len(rows) != 1:
                raise ValueError("Ambiguous S3 model/adjustment combination")
            plot_bar(axis, x + (j - 1.5) * 0.17, rows.iloc[0], 0.16, COLORS[j])
            plotted += 1
    axis.set_xticks([0, 1], ["EVI", "Ozone"])
    axis.set_ylim(-1.2, 1.25)
    axis.set_xlim(-0.52, 1.52)
    axis.set_ylabel("BAA difference per IQR increase (years)")
    axis.yaxis.set_major_locator(MultipleLocator(0.5))
    style_axis(axis)
    labels = ["Unadjusted", "Age and sex adjusted", "Multivariable adjusted", "Dual exposure, multivariable adjusted"]
    fig.legend(handles=[Patch(facecolor=color, label=label) for color, label in zip(COLORS, labels)],
               loc="upper center", bbox_to_anchor=(0.55, 0.21), ncol=2, fontsize=7.5)
    fig.text(0.13, 0.077, "Four-year mean exposure within 1 km; error bars: model-based 95% confidence intervals.", fontsize=7)
    fig.text(0.13, 0.038, "These supplementary REML models are distinct from the primary robust-variance models.", fontsize=7)
    assert plotted == 8
    return {"source": source_info(path, root), "plot_rows": plotted, "precision_decimals": 3,
            "geometry_equals_source": True, **export(fig, "FigS3_iqr_models", output, dpi)}


def render_s4(root: Path, output: Path, dpi: int) -> dict:
    aggregate_path = root / "03_source_data/FigS4_counterfactual_age_aggregates.csv"
    figures, audit_rows = [], []
    aggregate_only = aggregate_path.is_file()
    identity_error = None
    if aggregate_only:
        path = aggregate_path
        table = pd.read_csv(path)
        required = {"age_group", "pattern", "n_observations", "med", "q1", "q3", "whislo", "whishi"}
        if not required.issubset(table.columns) or len(table) != 10:
            raise ValueError("S4 aggregate source requires ten complete age-pattern groups")
        if table.duplicated(["age_group", "pattern"]).any():
            raise ValueError("Duplicate S4 aggregate age-pattern groups")
        table = table.set_index(["age_group", "pattern"]).loc[
            [(age, pattern) for age in AGE_GROUPS for pattern in PATTERNS]
        ].reset_index()
        for _, row in table.iterrows():
            stats = {key: float(row[key]) for key in ["med", "q1", "q3", "whislo", "whishi"]}
            if not stats["whislo"] <= stats["q1"] <= stats["med"] <= stats["q3"] <= stats["whishi"]:
                raise ValueError("Invalid S4 box-statistic ordering")
            figures.append({**stats, "fliers": []})
    else:
        raise FileNotFoundError("Required public S4 aggregate source is missing: " + str(aggregate_path))
    table.to_csv(output / "FigS4_verified_box_statistics.csv", index=False)
    expected_counts = [2631, 525, 1544, 344, 631, 124, 92, 23, 4898, 1016]
    if table.n_observations.tolist() != expected_counts:
        raise AssertionError("S4 age-pattern counts do not match the archived frame")
    fig, axis = plt.subplots(figsize=(180 * MM, 125 * MM))
    fig.subplots_adjust(left=0.12, right=0.99, top=0.77, bottom=0.25)
    fig.suptitle("Linear mixed-model BAA scenario contrasts by age and exposure pattern", fontsize=9.5, y=0.97)
    positions = [i + offset for i in range(5) for offset in [-0.18, 0.18]]
    artists = axis.bxp(figures, positions=positions, widths=0.28, patch_artist=True,
                      manage_ticks=False, showfliers=not aggregate_only,
                      medianprops={"color": "#1A1A1A", "linewidth": 1.0},
                      flierprops={"marker": "o", "markersize": 2, "markerfacecolor": "none",
                                  "markeredgecolor": "#666666", "markeredgewidth": 0.5},
                      whiskerprops={"linewidth": 0.75}, capprops={"linewidth": 0.75})
    for j, box in enumerate(artists["boxes"]):
        box.set_facecolor(COLORS[j % 2])
        box.set_edgecolor("#333333")
        box.set_linewidth(0.75)
        np.testing.assert_allclose(np.unique(box.get_path().vertices[:, 1]),
                                   [figures[j]["q1"], figures[j]["q3"]], rtol=0, atol=1e-15)
        np.testing.assert_allclose(artists["medians"][j].get_ydata(), figures[j]["med"], rtol=0, atol=1e-15)
        axis.text(positions[j], 6.02, f"n={expected_counts[j]:,}", ha="center", va="center", fontsize=6.8)
    axis.set_xticks(range(5), ["≤60 years", "61–70 years", "71–80 years", ">80 years", "All observations"])
    axis.set_xlim(-0.55, 4.55)
    axis.set_ylim(0.65, 6.25)
    axis.set_ylabel("Current minus scenario-predicted BAA (years)")
    axis.yaxis.set_major_locator(MultipleLocator(1))
    axis.grid(axis="y", color="#E3E3E3", linewidth=0.55)
    axis.set_axisbelow(True)
    fig.legend(handles=[Patch(facecolor=COLORS[i], label=label) for i, label in enumerate(PATTERN_LABELS)],
               loc="upper center", bbox_to_anchor=(0.55, 0.90), ncol=2, fontsize=8)
    fig.text(0.12, 0.16, r"Scenario: EVI below 0.40 set to 0.40; ozone above 60 μg m$^{-3}$ set to 60 μg m$^{-3}$.", fontsize=7)
    box_note = ("Aggregate reconstruction: median, IQR and 1.5-IQR whiskers; individual outliers not displayed."
                if aggregate_only else
                "Boxes: median and IQR; whiskers: observed values within 1.5 IQR; circles: outliers.")
    fig.text(0.12, 0.11, box_note, fontsize=7)
    fig.text(0.12, 0.06, "N counts participant-wave observations aged 46–89 years. Descriptive predictions, not causal effects.", fontsize=7)
    return {"source": source_info(path, root), "n_observations": 5914, "box_groups": len(figures),
            "counts": expected_counts, "quartile_rule": "Stata default; averaged_inverted_cdf",
            "aggregate_reconstruction": aggregate_only, "individual_outliers_displayed": not aggregate_only,
            "prediction_identity_max_abs_difference": identity_error,
            "geometry_equals_source": True, **export(fig, "FigS4_counterfactual_age", output, dpi)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--dpi", type=int, choices=[300, 600], default=600)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output_dir.resolve()
    formal = root / "02_figures"
    if output == formal or formal in output.parents:
        raise ValueError("Render candidates separately; do not overwrite formal figures")
    output.mkdir(parents=True, exist_ok=True)
    style()
    report = {"backend": "Python/matplotlib", "matplotlib": matplotlib.__version__,
              "numpy": np.__version__, "pandas": pd.__version__, "dpi": args.dpi,
              "S2": render_s2(root, output, args.dpi),
              "S3": render_s3(root, output, args.dpi),
              "S4": render_s4(root, output, args.dpi)}
    (output / "figure_source_comparison.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("Rendered S2-S4 directly from verified aggregate sources.")


if __name__ == "__main__":
    main()
