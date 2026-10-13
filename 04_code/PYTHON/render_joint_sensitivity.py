"""Render joint-pattern associations and complete sensitivity summaries.

Public aggregate inputs are read without model refitting or significance-based
selection. PDF and SVG retain editable text; PNG and TIFF use 600 dpi.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "svg.hashsalt": "joint-patterns-public",
    "savefig.facecolor": "white",
}
GROUP_COLORS = ["#32779d", "#c39438", "#31847c", "#a36789"]
GROUP_MARKERS = ["o", "D", "^", "s"]
GROUP_LABELS = [
    "Low EVI\nlow O$_3$", "Low EVI\nhigh O$_3$",
    "High EVI\nlow O$_3$", "High EVI\nhigh O$_3$",
]
SPECIFICATIONS = [
    ("unweighted_dual", "Primary, dual", 9.5),
    ("unweighted_single", "Single exposure", 8.5),
    ("ipw_dual", "IPW, dual", 7.5),
    ("ipw_single", "IPW, single", 6.5),
    ("mi", "MI, dual (20 sets)", 5.5),
    ("unweighted_dual", "Unweighted, dual", 3.5),
    ("unweighted_single", "Single exposure", 2.5),
    ("ipw_dual", "IPW, dual", 1.5),
    ("ipw_single", "IPW, single", 0.5),
]
CONTRASTS = [
    ("Low EVI/high O3 vs reference", "Low EVI / high O$_3$ vs ref."),
    ("High EVI/low O3 vs reference", "High EVI / low O$_3$ vs ref."),
    ("High EVI/high O3 vs reference", "High EVI / high O$_3$ vs ref."),
    ("O3 contrast within low EVI", "O$_3$ contrast within low EVI"),
    ("O3 contrast within high EVI", "O$_3$ contrast within high EVI"),
    ("Difference-in-differences interaction", "DiD interaction"),
]


def _check_intervals(frame: pd.DataFrame, estimate: str = "estimate") -> None:
    values = frame[[estimate, "lci", "uci"]].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("All displayed estimates and confidence bounds must be finite")
    if not ((frame.lci <= frame[estimate]) & (frame[estimate] <= frame.uci)).all():
        raise ValueError("Confidence bounds do not contain their estimates")


def load_continuous_models(root: Path) -> pd.DataFrame:
    """Return every archived continuous-model row, with MI on the IQR scale."""
    source = Path(root) / "03_source_data"
    primary = pd.read_csv(source / "primary_continuous_models_v5.csv")
    mi = pd.read_csv(source / "mi_sensitivity_v5.csv")
    keys = ["exposure", "scale", "weight_spec", "model"]
    expected = {
        (e, s, w, m)
        for e in ("EVI", "O3") for s in ("1KM", "500m")
        for w in ("unweighted", "ipw")
        for m in ("single-exposure", "dual-exposure")
    }
    actual = set(primary[keys].itertuples(index=False, name=None))
    if len(primary) != 16 or actual != expected or primary.duplicated(keys).any():
        raise ValueError("Expected all 16 unique archived continuous-model specifications")
    if len(mi) != 2 or set(mi.exposure) != {"EVI", "O3"}:
        raise ValueError("Expected both archived 1-km multiple-imputation estimates")
    if not (primary.n.eq(5914).all() and primary.n_participants.eq(3154).all()):
        raise ValueError("Continuous-model frame differs from the figure sample annotation")

    rows = []
    for r in primary.itertuples(index=False):
        specification = r.weight_spec + ("_dual" if r.model == "dual-exposure" else "_single")
        rows.append({
            "exposure": r.exposure, "scale": r.scale, "specification": specification,
            "estimate": r.beta_iqr, "lci": r.lci_iqr, "uci": r.uci_iqr, "p": r.p,
            "n_observations": r.n, "n_participants": r.n_participants,
            "iqr": r.exposure_iqr, "variance": "robust sandwich",
            "source": "primary_continuous_models_v5.csv",
        })
    for r in mi.itertuples(index=False):
        # Positive unit conversion only: the retained MI table stores native-unit
        # confidence bounds and both native-unit and per-IQR point estimates.
        iqr = r.mi_beta_iqr / r.mi_beta
        reference = primary[(primary.exposure == r.exposure) & primary.scale.eq("1KM")].iloc[0]
        if not (iqr > 0 and np.isclose(iqr, reference.exposure_iqr, atol=1e-10, rtol=0)):
            raise ValueError("MI IQR scaling differs from the archived primary IQR")
        rows.append({
            "exposure": r.exposure, "scale": "1KM", "specification": "mi",
            "estimate": r.mi_beta_iqr, "lci": r.mi_lci * iqr,
            "uci": r.mi_uci * iqr, "p": r.mi_p,
            "n_observations": np.nan, "n_participants": np.nan,
            "iqr": iqr, "variance": "Rubin pooled; model-based; archived normal CI",
            "source": "mi_sensitivity_v5.csv",
        })
    data = pd.DataFrame(rows)
    _check_intervals(data)
    if not data.p.between(0, 1).all():
        raise ValueError("P values must lie between zero and one")
    if not np.array_equal(data.iloc[:16].p.to_numpy(), primary.p.to_numpy()):
        raise ValueError("Primary P values changed during preparation")
    if not np.array_equal(data.iloc[16:].p.to_numpy(), mi.mi_p.to_numpy()):
        raise ValueError("MI P values changed during preparation")
    return data


def _load_means(root: Path) -> pd.DataFrame:
    means = pd.read_csv(Path(root) / "03_source_data" / "Figure3a_adjusted_means.csv")
    if len(means) != 8 or not means.n.eq(6471).all():
        raise ValueError("Expected eight joint-pattern means from 6,471 observations")
    for scale in ("1KM", "500m"):
        if sorted(means.loc[means.scale.eq(scale), "group_order"]) != [1, 2, 3, 4]:
            raise ValueError("Joint-pattern group ordering is incomplete")
    _check_intervals(means, "adjusted_mean")
    return means


def _draw_means(ax, means: pd.DataFrame, scale: str, letter: str) -> None:
    frame = means[means.scale.eq(scale)].sort_values("group_order")
    for i, r in enumerate(frame.itertuples(index=False)):
        ax.errorbar(i, r.adjusted_mean,
                    yerr=[[r.adjusted_mean - r.lci], [r.uci - r.adjusted_mean]],
                    fmt=GROUP_MARKERS[i], color=GROUP_COLORS[i], ecolor="#70787e",
                    elinewidth=1, capsize=3, markersize=5, zorder=3)
    ax.axhline(0, color="#9ba4ac", lw=0.8)
    ax.set(xlim=(-0.45, 3.45), ylim=(-3.15, 2.25), xticks=range(4),
           xticklabels=GROUP_LABELS, yticks=[-3, -2, -1, 0, 1, 2])
    ax.grid(axis="y", color="#e9edef", lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(("1 km" if scale == "1KM" else "500 m") + " joint exposure patterns",
                 loc="left", fontweight="bold", pad=12)
    ax.text(-0.13, 1.08, letter, transform=ax.transAxes, fontweight="bold", fontsize=10)
    ax.text(0, 1.005, "Adjusted means; joint-model frame n = 6,471",
            transform=ax.transAxes, fontsize=6.3, color="#51616d")
    ax.set_ylabel("Adjusted mean BAA (years)" if scale == "1KM" else "")


def _draw_sensitivity(ax, data: pd.DataFrame, exposure: str, letter: str) -> None:
    ax.set_title(("EVI" if exposure == "EVI" else "O$_3$") + ": continuous-model sensitivity",
                 loc="left", fontweight="bold", pad=10)
    ax.text(-0.05, 1.055, letter, transform=ax.transAxes, fontweight="bold", fontsize=10)
    ax.set(xlim=(-3.5, 2.5), ylim=(-0.1, 11.3))
    ax.axvline(0, ymin=0.03, ymax=0.86, color="#bdc5ca", lw=0.8)
    ax.axhspan(4.8, 10.45, color="#f5f7f8", zorder=0)
    ax.text(-3.4, 10.65, "1 km", weight="bold", color="#243f50", fontsize=7)
    ax.text(-3.4, 4.35, "500 m", weight="bold", color="#243f50", fontsize=7)
    ax.text(2.1, 10.65, "P", weight="bold", fontsize=7)
    for i, (spec, label, y) in enumerate(SPECIFICATIONS):
        scale = "1KM" if i < 5 else "500m"
        r = data[data.exposure.eq(exposure) & data.scale.eq(scale) & data.specification.eq(spec)].iloc[0]
        primary = spec == "unweighted_dual"
        color = "#243f50" if primary else ("#846b92" if spec.startswith("ipw") else "#347f83")
        marker = "s" if primary else ("D" if spec == "mi" else "o")
        ax.errorbar(r.estimate, y, xerr=[[r.estimate - r.lci], [r.uci - r.estimate]],
                    fmt=marker, color=color, markersize=4.4 if primary else 3.8,
                    markerfacecolor=color if primary else "white", markeredgewidth=1,
                    capsize=2.4, elinewidth=0.9, zorder=3)
        ax.text(-3.4, y, label, ha="left", va="center", fontsize=6.35,
                fontweight="bold" if primary else "normal")
        ax.text(2.1, y, f"{r.p:.3f}", va="center", ha="left", fontsize=6.5)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_bounds(-1.2, 1.8)
    ax.set_xticks([-1, 0, 1])
    ax.set_xlabel("BAA difference per IQR (years)", labelpad=5)
    ax.tick_params(axis="y", length=0)


def _export(fig, out: Path, stem: str) -> list[Path]:
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("pdf", "svg", "png", "tiff"):
        extra = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {}
        if ext == "pdf":
            extra["metadata"] = {"CreationDate": None, "ModDate": None}
        elif ext == "svg":
            extra["metadata"] = {"Date": None}
        path = out / f"{stem}.{ext}"
        fig.savefig(path, dpi=600, **extra)
        paths.append(path)
    plt.close(fig)
    return paths


def render_figure3(root: Path, out: Path) -> list[Path]:
    """Render the four-panel main figure from the complete public result set."""
    data = load_continuous_models(root)
    means = _load_means(root)
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(183 / 25.4, 199 / 25.4))
        grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.38], left=0.085, right=0.975,
                               top=0.895, bottom=0.195, wspace=0.19, hspace=0.46)
        fig.text(0.085, 0.975, "Joint exposure patterns and model sensitivity",
                 fontsize=11, weight="bold", va="top")
        _draw_means(fig.add_subplot(grid[0, 0]), means, "1KM", "a")
        _draw_means(fig.add_subplot(grid[0, 1]), means, "500m", "b")
        _draw_sensitivity(fig.add_subplot(grid[1, 0]), data, "EVI", "c")
        _draw_sensitivity(fig.add_subplot(grid[1, 1]), data, "O3", "d")
        fig.text(0.085, 0.100, "Whiskers: 95% CI. Primary/IPW: robust variance; 5,914 observations from 3,154 participants.", fontsize=6.6)
        fig.text(0.085, 0.079, "MI: 20 imputations, model-based variance with Rubin pooling; MI sample counts are unspecified.", fontsize=6.6)
        fig.text(0.085, 0.058, "DiD P = 0.605 / 0.779 (1 km / 500 m); interaction and precision results are shown in Fig. S9.", fontsize=6.6)
        fig.text(0.085, 0.037, "Panels use different estimands and sample frames; these associations do not establish effect modification.", fontsize=6.6)
        return _export(fig, out, "Figure3_joint_patterns")


def load_joint_inference(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    source = Path(root) / "03_source_data"
    contrasts = pd.read_csv(source / "Figure3b_contrasts.csv")
    precision = pd.read_csv(source / "Figure3b_DiD_precision.csv")
    expected = {(scale, key) for scale in ("1KM", "500m") for key, _ in CONTRASTS}
    if len(contrasts) != 12 or set(zip(contrasts.scale, contrasts.contrast)) != expected:
        raise ValueError("Expected all 12 joint-model contrast rows")
    if len(precision) != 2 or set(precision.scale) != {"1KM", "500m"}:
        raise ValueError("Expected precision estimates for both spatial buffers")
    if not contrasts.n.eq(6471).all() or not contrasts.p.between(0, 1).all():
        raise ValueError("Joint-model sample or nominal P values are invalid")
    _check_intervals(contrasts)
    _check_intervals(precision)
    for r in precision.itertuples(index=False):
        d = contrasts[contrasts.scale.eq(r.scale) & contrasts.contrast.eq(CONTRASTS[-1][0])].iloc[0]
        if not np.allclose([r.estimate, r.lci, r.uci], [d.estimate, d.lci, d.uci], atol=1e-12, rtol=0):
            raise ValueError("DiD contrast and precision source files disagree")
    return contrasts, precision


def _draw_contrasts(ax, contrasts: pd.DataFrame) -> None:
    ax.set(xlim=(-8.8, 8.9), ylim=(-0.2, 15.9))
    ax.text(-8.8, 15.75, "a", weight="bold", fontsize=10)
    ax.text(-8.05, 15.75, "Adjusted contrasts and interaction", weight="bold", fontsize=8.5)
    ax.text(4.2, 14.75, "Estimate [95% CI]", fontsize=6.7, weight="bold")
    ax.text(8.0, 14.75, "P", fontsize=6.7, weight="bold")
    ax.axvline(0, ymin=0.02, ymax=0.88, lw=0.7, color="#b7c1c8")
    for scale, first_y in (("1KM", 13.25), ("500m", 5.75)):
        ax.text(-8.6, first_y + 0.85, "1 km" if scale == "1KM" else "500 m",
                weight="bold", fontsize=7, color="#243f50")
        for index, (key, label) in enumerate(CONTRASTS):
            r = contrasts[contrasts.scale.eq(scale) & contrasts.contrast.eq(key)].iloc[0]
            y = first_y - index
            did = key.startswith("Difference")
            color = "#c48741" if did else ("#347f83" if "within" in key else "#365f7d")
            if index % 2 == 0:
                ax.axhspan(y - 0.43, y + 0.43, color="#f5f7f8", zorder=0)
            ax.errorbar(r.estimate, y, xerr=[[r.estimate - r.lci], [r.uci - r.estimate]],
                        fmt="^" if did else "o", color=color, markersize=4.3,
                        capsize=2.7, elinewidth=0.9, zorder=3)
            ax.text(-8.6, y, label, fontsize=6.5, va="center")
            ax.text(4.2, y, f"{r.estimate:+.2f} [{r.lci:+.2f}, {r.uci:+.2f}]",
                    fontsize=6.5, va="center")
            ax.text(8.0, y, f"{r.p:.3f}", fontsize=6.5, va="center")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_bounds(-3.5, 3.5)
    ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])
    ax.set_xlabel("Adjusted contrast in BAA (years)", labelpad=4)


def _draw_precision(ax, precision: pd.DataFrame, contrasts: pd.DataFrame) -> None:
    ax.set(xlim=(-4.7, 8.6), ylim=(-0.8, 3.4))
    ax.text(-4.7, 3.2, "b", fontsize=10, weight="bold")
    ax.text(-4.1, 3.2, "DiD precision", fontsize=8.5, weight="bold")
    ax.axvspan(-1, 1, ymin=0.17, ymax=0.73, color="#f3f0df", zorder=0)
    ax.text(0, 2.45, "±1-year equivalence margin", fontsize=6.3, color="#796f3e", ha="center")
    ax.plot([0, 0], [-0.05, 2.2], lw=0.7, color="#8b969e", zorder=1)
    for bound in (-1, 1):
        ax.plot([bound, bound], [-0.05, 2.2], lw=0.6, color="#b8ae74", ls="--")
    for xpos, label in ((3.5, "DiD P"), (5.0, "80% MDE (y)"), (7.3, "TOST P")):
        ax.text(xpos, 2.35, label, fontsize=6.4, weight="bold", ha="center")
    for index, scale in enumerate(("1KM", "500m")):
        r = precision[precision.scale.eq(scale)].iloc[0]
        p = contrasts[contrasts.scale.eq(scale) & contrasts.contrast.eq(CONTRASTS[-1][0])].iloc[0].p
        y = 1.7 - index * 1.45
        ax.text(-4.6, y, "1 km" if scale == "1KM" else "500 m", fontsize=6.6, weight="bold", va="center")
        ax.errorbar(r.estimate, y, xerr=[[r.estimate - r.lci], [r.uci - r.estimate]],
                    fmt="^", color="#c48741", markersize=5, capsize=3, elinewidth=1, zorder=3)
        ax.text(0, y - 0.47, f"{r.estimate:+.2f} [{r.lci:+.2f}, {r.uci:+.2f}]",
                fontsize=6.5, ha="center")
        for xpos, value in ((3.5, f"{p:.3f}"), (5.0, f"{r.mde_80pct:.2f}"), (7.3, f"{r.tost_p:.3f}")):
            ax.text(xpos, y, value, fontsize=6.6, ha="center", va="center")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_bounds(-3.5, 3.5)
    ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])
    ax.set_xlabel("DiD estimate (years)", labelpad=4)


def render_joint_inference(root: Path, out: Path) -> list[Path]:
    """Preserve all joint contrasts and interaction precision in Figure S9."""
    contrasts, precision = load_joint_inference(root)
    with plt.rc_context(STYLE):
        fig = plt.figure(figsize=(183 / 25.4, 199 / 25.4))
        fig.text(0.055, 0.975, "Joint exposure contrasts and interaction precision",
                 fontsize=10.5, weight="bold", va="top")
        _draw_contrasts(fig.add_axes([0.055, 0.37, 0.91, 0.56]), contrasts)
        _draw_precision(fig.add_axes([0.055, 0.12, 0.91, 0.20]), precision, contrasts)
        fig.text(0.055, 0.056, "Whiskers: 95% CI. Joint-model frame: 6,471 observations. Reference: low EVI / low O$_3$.", fontsize=6.5)
        fig.text(0.055, 0.035, "MDE: minimum detectable interaction at 80% power, not an exclusion bound. Equivalence was not established.", fontsize=6.5)
        return _export(fig, out, "FigS9_joint_contrasts_precision")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--out", type=Path, help="Output directory for the main Figure 3")
    parser.add_argument("--si-out", type=Path, help="Output directory for Supplementary Figure S9")
    args = parser.parse_args()
    if args.out is None and args.si_out is None:
        parser.error("Specify --out and/or --si-out")
    if args.out is not None:
        render_figure3(args.root.resolve(), args.out.resolve())
        print("Rendered Figure 3 from public aggregates")
    if args.si_out is not None:
        render_joint_inference(args.root.resolve(), args.si_out.resolve())
        print("Rendered Figure S9 from public aggregates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
