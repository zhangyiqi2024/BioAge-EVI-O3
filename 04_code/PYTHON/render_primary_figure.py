"""Render Figure 2 from public exposure bins and archived model summaries.

No participant records, new fits, reconstructed small counts or simulated
observations are used. PDF/SVG text is editable; PNG/TIFF exports use 600 dpi.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd


STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 8,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "svg.hashsalt": "primary-associations-public",
    "savefig.facecolor": "white",
}
GROUPS = ["Low EVI / low O3", "Low EVI / high O3",
          "High EVI / low O3", "High EVI / high O3"]
GROUP_COLORS = ["#32779d", "#c39438", "#31847c", "#a36789"]
OLD_CENTRES = {"Low EVI / high O3": 1, "Low EVI / low O3": 2,
               "High EVI / high O3": 3, "High EVI / low O3": 4}
GROUP_LABELS = ["Low EVI\nlow O$_3$", "Low EVI\nhigh O$_3$",
                "High EVI\nlow O$_3$", "High EVI\nhigh O$_3$"]
SOURCE_FILES = [
    "Figure4_joint_density_bins.csv", "Figure4_exposure_histograms.csv",
    "Figure4_plot_metadata.csv", "Figure2_continuous_forest.csv",
    "Figure2a_raincloud_polygons.csv", "Figure2a_raincloud_boxplot_statistics.csv",
    "Figure2a_raincloud_annotations.csv",
]


def _suppression(frame: pd.DataFrame) -> np.ndarray:
    flags = frame.suppressed_small_count.astype(str).str.lower()
    if not flags.isin(["true", "false"]).all():
        raise ValueError("Small-count suppression flags must be explicit booleans")
    suppressed = flags.eq("true").to_numpy()
    counts = frame.n_observations.to_numpy(dtype=float)
    if not np.array_equal(np.isnan(counts), suppressed):
        raise ValueError("Only suppressed cells may have missing counts")
    released = counts[~suppressed]
    if not (np.isfinite(released).all() and (released >= 0).all()
            and np.equal(released, np.floor(released)).all()
            and ((released == 0) | (released >= 5)).all()):
        raise ValueError("Released counts must be zero or integers of at least five")
    return suppressed


def load_sources(root: Path) -> dict[str, pd.DataFrame]:
    """Load and validate every displayed source without changing its values."""
    frames = [pd.read_csv(Path(root) / "03_source_data" / name) for name in SOURCE_FILES]
    joint, hist, meta, forest, polygons, boxes, summaries = frames
    if len(joint) != 400 or joint.duplicated(["evi_left", "ozone_left"]).any():
        raise ValueError("Expected one released row for each of 20 by 20 joint bins")
    if joint.evi_left.nunique() != 20 or joint.ozone_left.nunique() != 20:
        raise ValueError("Joint exposure grid must have 20 bins along each axis")
    _suppression(joint)
    _suppression(hist)
    if len(hist) != 80 or set(hist.exposure) != {"EVI", "O3"}:
        raise ValueError("Expected both forty-bin marginal exposure histograms")
    if any(len(hist[hist.exposure.eq(e)]) != 40 for e in ("EVI", "O3")):
        raise ValueError("Each marginal histogram must contain forty bins")
    if len(meta) != 2 or not meta.n_observations.eq(5914).all():
        raise ValueError("Primary exposure frame must be 5,914 observations")
    core = forest[forest.weight_spec.eq("unweighted") & forest.model.eq("dual-exposure")].copy()
    expected_models = {(e, s) for e in ("EVI", "O3") for s in ("1KM", "500m")}
    if len(core) != 4 or set(zip(core.exposure, core.scale)) != expected_models:
        raise ValueError("Expected all four unweighted dual-exposure coefficients")
    if not (core.n.eq(5914).all() and core.n_participants.eq(3154).all()):
        raise ValueError("Primary models must retain their observation/participant counts")
    values = core[["beta_iqr", "lci_iqr", "uci_iqr", "p"]].to_numpy()
    if not (np.isfinite(values).all() and core.p.between(0, 1).all()
            and core.lci_iqr.le(core.beta_iqr).all() and core.uci_iqr.ge(core.beta_iqr).all()):
        raise ValueError("Invalid model estimate, interval or P value")
    expected_groups = {(s, g) for s in ("1KM", "500m") for g in GROUPS}
    for data in (boxes, summaries):
        if len(data) != 8 or set(zip(data.scale, data.group)) != expected_groups:
            raise ValueError("Expected eight complete group summaries")
    for scale, group in expected_groups:
        shape = polygons[polygons.scale.eq(scale) & polygons.group.eq(group)]
        old = OLD_CENTRES[group]
        if len(shape) != 1024 or not shape.x.between(old + .029, old + .271).all():
            raise ValueError("Archived half-density geometry has changed")
        if not np.isclose(shape.baa_years.min(), -40) or not np.isclose(shape.baa_years.max(), 70):
            raise ValueError("Archive outcome range must remain -40 to 70 years")
        box = boxes[boxes.scale.eq(scale) & boxes.group.eq(group)].iloc[0]
        if not box.whisker_low <= box.q1 <= box["median"] <= box.q3 <= box.whisker_high:
            raise ValueError("Box-summary ordering is invalid")
    if not summaries.groupby("scale").n.sum().eq(6472).all():
        raise ValueError("Each descriptive buffer must contain 6,472 observations")
    return dict(joint=joint, hist=hist, metadata=meta, models=core,
                polygons=polygons, boxes=boxes, summaries=summaries)


def _heading(fig, x: float, y: float, letter: str, title: str, subtitle: str) -> None:
    fig.text(x, y, letter, fontsize=10, weight="bold", va="top")
    fig.text(x + .026, y, title, fontsize=8, weight="bold", va="top")
    fig.text(x + .026, y - .027, subtitle, fontsize=6.5, color="#52616b", va="top")


def _draw_exposure(fig, sources: dict) -> None:
    joint, hist = sources["joint"], sources["hist"]
    ax = fig.add_axes([.103, .605, .302, .222])
    top = fig.add_axes([.103, .835, .302, .054], sharex=ax)
    side = fig.add_axes([.413, .605, .059, .222], sharey=ax)
    cmap = LinearSegmentedColormap.from_list("exposure_counts", ["#edf3f6", "#83b2bf", "#2e718b", "#193d59"])
    norm = Normalize(0, 240)
    for r in joint.itertuples(index=False):
        suppressed = str(r.suppressed_small_count).lower() == "true"
        # Missing cells are explicitly hatched; true zero is genuinely white.
        color = "#eeeeee" if suppressed else ("white" if r.n_observations == 0 else cmap(norm(r.n_observations)))
        rect = Rectangle((r.evi_left, r.ozone_left), r.evi_right-r.evi_left,
                         r.ozone_right-r.ozone_left, facecolor=color,
                         edgecolor="#949ba1" if suppressed else "none",
                         hatch="////" if suppressed else None, linewidth=.15)
        ax.add_patch(rect)
    ax.set(xlim=(joint.evi_left.min(), joint.evi_right.max()),
           ylim=(joint.ozone_left.min(), joint.ozone_right.max()),
           xlabel="EVI", ylabel="O$_3$ (µg m$^{-3}$)",
           xticks=[.2, .3, .4, .5], yticks=[70, 80, 90])
    ax.tick_params(length=2.5, pad=2)
    for e, marginal in (("EVI", top), ("O3", side)):
        h = hist[hist.exposure.eq(e)]
        released = h[~h.suppressed_small_count.astype(str).str.lower().eq("true")]
        suppressed = h[h.suppressed_small_count.astype(str).str.lower().eq("true")]
        if e == "EVI":
            marginal.bar(released.left, released.n_observations, width=released.right-released.left,
                         align="edge", color="#7398a7", linewidth=0)
            for r in suppressed.itertuples(index=False):
                marginal.plot(r.center, 0, marker="x", color="#7b858d", ms=4.2,
                              mew=.8, linestyle="none", clip_on=False, zorder=5)
            marginal.set_ylim(0, float(released.n_observations.max()) * 1.05)
            marginal.set_ylabel("Count", fontsize=5.8, labelpad=1)
            marginal.set_yticks([0, 500])
            marginal.tick_params(axis="x", bottom=False, labelbottom=False)
        else:
            marginal.barh(released.left, released.n_observations, height=released.right-released.left,
                          align="edge", color="#a5b6be", linewidth=0)
            for r in suppressed.itertuples(index=False):
                # A nonquantitative marker avoids implying an imputed bar length.
                marginal.plot(0, r.center, marker="x", color="#7b858d", ms=4.2,
                              mew=.8, linestyle="none", clip_on=False, zorder=5)
            marginal.set_xlim(0, float(released.n_observations.max()) * 1.05)
            marginal.set_xticks([0, 500])
            marginal.set_xlabel("Count", fontsize=5.8, labelpad=1)
            marginal.tick_params(axis="y", left=False, labelleft=False)
        marginal.tick_params(labelsize=5.5, length=2, pad=1)
        marginal.spines[["top", "right"]].set_visible(False)
    cax = fig.add_axes([.105, .540, .282, .010])
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal",
                      ticks=[0, 60, 120, 180, 240])
    bar.outline.set_visible(False)
    bar.ax.tick_params(labelsize=5.5, length=2, pad=1)
    bar.set_label("Observations per joint exposure bin", fontsize=6, labelpad=1)
    handles = [Patch(facecolor="white", edgecolor="#c0c7cc", label="Zero"),
               Patch(facecolor="#eeeeee", edgecolor="#949ba1", hatch="////", label="1–4 (joint)"),
               Line2D([], [], marker="x", color="#7b858d", linestyle="none", markersize=4,
                      label="1–4 (margin)")]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(.080, .484),
               ncol=3, fontsize=5.8, frameon=False, handlelength=1.3,
               handleheight=1.2, columnspacing=1.0)


def _draw_models(fig, models: pd.DataFrame) -> None:
    ax = fig.add_axes([.562, .578, .409, .306])
    ax.set(xlim=(-1.86, 1.75), ylim=(-.78, 5.8), yticks=[])
    ax.axvspan(-1.15, 0, color="#eef5f2", ymin=.065, ymax=.91, zorder=0)
    ax.plot([0, 0], [-.3, 5.23], color="#a8b1b6", lw=.7, zorder=1)
    ax.text(-1.81, 5.55, "1 km · primary", weight="bold", fontsize=7, color="#263f4e")
    ax.text(-1.81, 2.40, "500 m · buffer sensitivity", weight="bold", fontsize=7, color="#263f4e")
    ax.text(1.5, 5.55, "P", fontsize=7, weight="bold", ha="center")
    for scale, ys in (("1KM", [4.75, 3.52]), ("500m", [1.57, .34])):
        for exposure, y in zip(("EVI", "O3"), ys):
            r = models[models.scale.eq(scale) & models.exposure.eq(exposure)].iloc[0]
            color = "#31847c" if exposure == "EVI" else "#b68437"
            marker = "o" if scale == "1KM" else "s"
            ax.text(-1.81, y, "EVI" if exposure == "EVI" else "O$_3$", va="center", fontsize=7, color=color, weight="bold")
            ax.errorbar(r.beta_iqr, y, xerr=[[r.beta_iqr-r.lci_iqr], [r.uci_iqr-r.beta_iqr]],
                        fmt=marker, markersize=5, color=color, capsize=3,
                        elinewidth=1.1, markerfacecolor=color if scale == "1KM" else "white",
                        markeredgewidth=1.1, zorder=3)
            ax.text(0, y-.45, f"{r.beta_iqr:+.2f} [{r.lci_iqr:+.2f}, {r.uci_iqr:+.2f}]",
                    fontsize=6.1, ha="center", color="#43535e")
            ax.text(1.5, y, f"{r.p:.3f}", fontsize=6.5, ha="center", va="center")
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_bounds(-1.2, 1.2)
    ax.set_xticks([-1, -.5, 0, .5, 1], ["−1", "−0.5", "0", "0.5", "1"])
    ax.tick_params(axis="x", length=2.5, pad=3)
    ax.set_xlabel("BAA difference per IQR (years)", labelpad=4)
    fig.text(.574, .505, "Mutually adjusted; unweighted models", fontsize=6.5, color="#52616b")
    fig.text(.574, .484, "Points and whiskers: estimates and 95% CI", fontsize=6.5, color="#52616b")


def _draw_observed(ax, sources: dict, scale: str) -> None:
    polygons, boxes, summaries = (sources[key] for key in ("polygons", "boxes", "summaries"))
    for centre, (group, color) in enumerate(zip(GROUPS, GROUP_COLORS), 1):
        p = polygons[polygons.scale.eq(scale) & polygons.group.eq(group)]
        # Translate only the category coordinate; preserve vertex order and BAA.
        new_x = p.x.to_numpy() - OLD_CENTRES[group] + centre
        ax.fill(new_x, p.baa_years.to_numpy(), color=color, alpha=.64, lw=0)
        b = boxes[boxes.scale.eq(scale) & boxes.group.eq(group)].iloc[0]
        s = summaries[summaries.scale.eq(scale) & summaries.group.eq(group)].iloc[0]
        box_x = b.x - OLD_CENTRES[group] + centre
        ax.bxp([{"med": b["median"], "q1": b.q1, "q3": b.q3,
                 "whislo": b.whisker_low, "whishi": b.whisker_high, "fliers": []}],
               positions=[box_x], widths=.13, manage_ticks=False, showfliers=False,
               patch_artist=True, boxprops={"facecolor": "white", "edgecolor": color, "linewidth": .85},
               medianprops={"color": "#34424c", "linewidth": .8},
               whiskerprops={"color": "#53606a", "linewidth": .7},
               capprops={"color": "#53606a", "linewidth": .7})
        ax.plot(box_x, s["mean"], "o", ms=3.7, mfc="white", mec="#253641", mew=.8, zorder=4)
        ax.text(centre, 67, f"n = {int(s.n):,}", fontsize=6.2, ha="center", color="#44545e")
    ax.set(xlim=(.50, 4.48), ylim=(-40, 73), xticks=range(1, 5),
           xticklabels=GROUP_LABELS, yticks=[-40, -20, 0, 20, 40, 60])
    ax.set_ylabel("Observed BAA (years)")
    ax.axhline(0, color="#adb7bd", ls=(0, (3, 3)), lw=.7, zorder=0)
    ax.grid(axis="y", color="#edf0f2", lw=.55)
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", length=0, pad=6, labelsize=6.2)
    ax.tick_params(axis="y", length=2.5, pad=2)


def render_figure2(root: Path, out: Path) -> list[Path]:
    """Render the approved four-panel Figure 2 from public aggregates only."""
    sources = load_sources(Path(root))
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(183 / 25.4, 194 / 25.4))
        fig.text(.065, .977, "Exposure distributions and primary associations", fontsize=10.5, weight="bold", va="top")
        _heading(fig, .065, .932, "a", "Joint exposure distribution, 1 km", "Primary model frame: 5,914 observations")
        _heading(fig, .550, .932, "b", "Continuous-exposure associations", "5,914 observations from 3,154 participants")
        _draw_exposure(fig, sources)
        _draw_models(fig, sources["models"])
        _heading(fig, .065, .442, "c", "Observed BAA, 1 km", "Descriptive frame: 6,472 observations")
        _heading(fig, .550, .442, "d", "Observed BAA, 500 m", "Descriptive frame: 6,472 observations")
        _draw_observed(fig.add_axes([.102, .141, .370, .244]), sources, "1KM")
        _draw_observed(fig.add_axes([.587, .141, .370, .244]), sources, "500m")
        fig.text(.065, .065, "a: hatched cells and × marginal markers indicate counts of 1–4; these counts are suppressed.", fontsize=6.5, color="#52616b")
        fig.text(.065, .043, "c,d: archived half-densities and boxes; white points show means. Individual observations are not displayed.", fontsize=6.5, color="#52616b")
        paths = []
        for extension in ("pdf", "svg", "png", "tiff"):
            settings = {"pil_kwargs": {"compression": "tiff_lzw"}} if extension == "tiff" else {}
            if extension == "pdf":
                settings["metadata"] = {"CreationDate": None, "ModDate": None}
            elif extension == "svg":
                settings["metadata"] = {"Date": None}
            path = out / f"Figure2_primary_associations.{extension}"
            fig.savefig(path, dpi=600, **settings)
            paths.append(path)
        plt.close(fig)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    render_figure2(args.root.resolve(), args.out.resolve())
    print("Rendered Figure 2 from public aggregates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
