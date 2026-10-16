"""Draw the documented LightGBM scenario summaries in Tables S13 and S14.

This replaces an irreproducible historical image with the retained 5,914-row
analysis summaries. It does not refit a model or recreate individual points.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

AGES = ["<=60 years", "61-70 years", "71-80 years", ">80 years", "All participants"]
AGE_LABELS = ["≤60 years", "61–70 years", "71–80 years", ">80 years", "All observations"]
PATTERNS = ["Both exposures adverse", "Only O3 adverse"]


def render(root: Path, output: Path) -> dict:
    source = root / "03_source_data"
    mean_path = source / "s13_scenario_reductions_by_age_pattern_v4.csv"
    distribution_path = source / "s14_scenario_reduction_distribution_v4.csv"
    means = pd.read_csv(mean_path).set_index(["AgeGroup", "Pattern"])
    distribution = pd.read_csv(distribution_path).set_index(["AgeGroup", "Pattern"])
    expected_index = [(age, pattern) for age in AGES for pattern in PATTERNS]
    if set(means.index) != set(expected_index) or set(distribution.index) != set(expected_index):
        raise ValueError("Expected ten age-by-exposure-pattern aggregate rows")
    if int(means.loc["All participants", "N"].sum()) != 5914:
        raise ValueError("Retained LightGBM scenario sample must contain 5,914 observations")
    np.testing.assert_array_equal(means.loc[expected_index, "N"], distribution.loc[expected_index, "N"])
    np.testing.assert_allclose(means.CI_lower, means.Mean_reduction - 1.96 * means.SE)
    np.testing.assert_allclose(means.CI_upper, means.Mean_reduction + 1.96 * means.SE)
    if not ((distribution.Q1 <= distribution.Median) & (distribution.Median <= distribution.Q3)).all():
        raise ValueError("Invalid interquartile summary")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9,
                         "svg.hashsalt": "bioage-public", "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 4.5), sharey=True)
    colours = ["#0072B2", "#D55E00"]
    for j, pattern in enumerate(PATTERNS):
        order = [(age, pattern) for age in AGES]
        m = means.loc[order]
        d = distribution.loc[order]
        y = np.arange(5) + (-0.15 if j == 0 else 0.15)
        label = pattern.replace("O3", "O₃")
        for axis, values, lower, upper in [
            (axes[0], m.Mean_reduction, m.CI_lower, m.CI_upper),
            (axes[1], d.Median, d.Q1, d.Q3),
        ]:
            artists = axis.errorbar(values, y, xerr=[values-lower, upper-values],
                                    fmt="o" if j == 0 else "s", markersize=4.5, capsize=3,
                                    color=colours[j], linewidth=1, label=label)
            np.testing.assert_allclose(np.asarray(artists.lines[0].get_xdata(), dtype=float), values.to_numpy(dtype=float))
    axes[0].set_yticks(range(5), AGE_LABELS)
    axes[0].invert_yaxis()
    axes[0].set_title("a  Mean and 95% confidence interval", loc="left", fontweight="bold")
    axes[1].set_title("b  Median and interquartile range", loc="left", fontweight="bold")
    for axis in axes:
        axis.axvline(0, color="#999999", linewidth=0.7, linestyle="--")
        axis.grid(axis="y", color="#E6E6E6", linewidth=0.5)
        axis.spines[["top", "right"]].set_visible(False)
        axis.set_xlabel("Current minus scenario-predicted BAA (years)")
    fig.suptitle("LightGBM scenario summaries by age and exposure pattern", fontsize=10, y=0.98)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.56, 0.91),
               ncol=2, frameon=False)
    fig.subplots_adjust(left=0.17, right=0.97, top=0.77, bottom=0.30, wspace=0.14)
    fig.text(0.17, 0.15, "5,914 observations: 4,898 with both exposures adverse; 1,016 with only O₃ adverse.", fontsize=7)
    fig.text(0.17, 0.10, "Scenario targets: EVI ≥0.40 and O₃ ≤60 μg m⁻³. Sources: Supplementary Tables S13–S14.", fontsize=7)
    fig.text(0.17, 0.05, "Descriptive model outputs; negligible predictive skill; no causal interpretation.", fontsize=7)
    output.mkdir(parents=True, exist_ok=True)
    files = []
    for ext in ["png", "tiff", "pdf", "svg"]:
        path = output / f"FigS5_lightgbm_scenarios.{ext}"
        opts = {"dpi": 600, "facecolor": "white"}
        if ext == "tiff": opts["pil_kwargs"] = {"compression": "tiff_lzw"}
        if ext == "pdf": opts["metadata"] = {"CreationDate": None, "ModDate": None}
        if ext == "svg": opts["metadata"] = {"Date": None}
        fig.savefig(path, **opts)
        files.append({"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    plt.close(fig)
    report = {"sources": [mean_path.relative_to(root).as_posix(), distribution_path.relative_to(root).as_posix()],
              "n_observations": 5914, "aggregate_rows": 10, "geometry_checks_pass": True,
              "individual_records_distributed": False, "files": files}
    (output / "FigS5_render_verification.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(render(args.root.resolve(), args.output_dir.resolve()), indent=2))
