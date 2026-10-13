"""Reconstruct the observed-distribution panel from aggregate geometry."""
from __future__ import annotations
from io import BytesIO
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

GROUPS = ["Low EVI / high O3", "Low EVI / low O3", "High EVI / high O3", "High EVI / low O3"]


def raincloud_panel(root: Path) -> Image.Image:
    source = root / "03_source_data"
    polygons = pd.read_csv(source / "Figure2a_raincloud_polygons.csv")
    boxes = pd.read_csv(source / "Figure2a_raincloud_boxplot_statistics.csv").set_index(["scale", "group"])
    summaries = pd.read_csv(source / "Figure2a_raincloud_annotations.csv").set_index(["scale", "group"])
    if len(boxes) != 8 or len(summaries) != 8:
        raise ValueError("Expected eight public distribution summaries")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
    fig, axes = plt.subplots(2, 1, figsize=(10.8, 7.8), sharex=True, sharey=True)
    colours = ["#E69F00", "#0072B2", "#CC79A7", "#009E73"]
    for axis, scale in zip(axes, ["1KM", "500m"]):
        for j, (group, colour) in enumerate(zip(GROUPS, colours), 1):
            d = polygons[(polygons.scale == scale) & (polygons.group == group)]
            if len(d) == 0: raise ValueError("Missing raincloud density polygon")
            axis.fill(d.x, d.baa_years, facecolor=colour, alpha=0.58, edgecolor="none")
            b = boxes.loc[(scale, group)]
            m = summaries.loc[(scale, group)]
            if not b.whisker_low <= b.q1 <= b['median'] <= b.q3 <= b.whisker_high:
                raise ValueError("Invalid boxplot statistics")
            axis.bxp([{"med": b['median'], "q1": b.q1, "q3": b.q3,
                       "whislo": b.whisker_low, "whishi": b.whisker_high, "fliers": []}],
                     positions=[j-0.11], widths=0.11, showfliers=False, manage_ticks=False,
                     patch_artist=True, boxprops={"facecolor": "white", "linewidth": 1},
                     medianprops={"color": "black"}, whiskerprops={"color": "black"})
            axis.errorbar(j-0.11, m['mean'], yerr=[[m['mean']-m.ci_lower], [m.ci_upper-m['mean']]],
                          fmt="o", markersize=4, markerfacecolor="white", markeredgecolor="black",
                          ecolor="black", capsize=3, zorder=4)
            axis.text(j, 67, f"n = {int(m['n']):,}", ha="center", fontsize=9)
            axis.text(j+0.04, 49, f"Mean {m['mean']:+.2f}\n[{m.ci_lower:+.2f}, {m.ci_upper:+.2f}]",
                      ha="left", va="top", fontsize=7)
        if int(summaries.loc[scale, "n"].sum()) != 6472:
            raise ValueError("Observed distribution sample differs from 6,472 observations")
        axis.axhline(0, color="#888888", linestyle="--", linewidth=0.8)
        axis.set(ylim=(-40, 74), xlim=(0.55, 4.55), ylabel="Observed BAA (years)")
        axis.set_title("1 km residential buffer" if scale == "1KM" else "500 m residential buffer", loc="left")
        axis.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xticks(range(1, 5), [g.replace(" / ", "\n").replace("O3", "O₃") for g in GROUPS])
    fig.subplots_adjust(left=0.08, right=0.99, top=0.95, bottom=0.14, hspace=0.20)
    fig.text(0.5, 0.035, "Densities and box summaries; white points show means with 95% CI. Individual observations are not displayed.",
             ha="center", fontsize=8)
    stream = BytesIO()
    fig.savefig(stream, format="png", dpi=300, facecolor="white")
    plt.close(fig)
    stream.seek(0)
    return Image.open(stream).convert("RGB")
