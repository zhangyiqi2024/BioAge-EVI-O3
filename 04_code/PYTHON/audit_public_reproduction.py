"""Verify public figures and displayed tables using only the release inputs.

Reconstruction starts at aggregate results, not restricted participant records.
Missing inputs, missing outputs or altered table cells produce a nonzero exit.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path

from PIL import Image
from render_public_tables import CellReader, display_cells, load_sources

STATES = ("REGENERATED_FROM_PUBLIC_SOURCE", "VERIFIED_AGAINST_ARCHIVED_OUTPUT",
          "ARCHIVED_ONLY", "NOT_REPRODUCIBLE_FROM_RELEASE")

# Every distributed source used by each figure is listed explicitly.
FIGURES = [
    ("Figure 1", "main", "Figure1_study_design", ["FigureS1_participant_flow.csv"],
     "Aggregate flow counts; timing and conceptual DAG encoded in the renderer."),
    ("Figure 2", "main", "Figure2_primary_associations",
     ["Figure2_continuous_forest.csv", "Figure2a_raincloud_polygons.csv",
      "Figure2a_raincloud_boxplot_statistics.csv", "Figure2a_raincloud_annotations.csv",
      "Figure4_joint_density_bins.csv", "Figure4_exposure_histograms.csv", "Figure4_plot_metadata.csv"],
     "Joint exposure counts with suppressed small cells, four core unweighted continuous estimates, and observed BAA density/box/mean summaries; no individual points."),
    ("Figure 3", "main", "Figure3_joint_patterns",
     ["Figure3a_adjusted_means.csv", "primary_continuous_models_v5.csv", "mi_sensitivity_v5.csv"],
     "Four-group adjusted means and all 18 archived continuous/MI estimates; different sample frames and variance estimators are disclosed."),
    ("Figure 4", "main", "Figure4_exposure_response",
     ["Figure4_exposure_response_v5.csv", "Figure4_joint_contour_v5.csv",
      "Figure4_plot_metadata.csv", "Figure4_exposure_histograms.csv"],
     "Evaluated spline and linear grids, IQRs and binned exposures; restricted point layer omitted."),
    ("Supplementary Figure S1", "supplementary", "FigS1_participant_flow",
     ["FigureS1_participant_flow.csv"], "Aggregate participant-flow counts."),
    ("Supplementary Figure S2", "supplementary", "FigS2_quartile_scores",
     ["quartile_scores_detailed_results.csv"], "Quartile-score coefficients and confidence intervals."),
    ("Supplementary Figure S3", "supplementary", "FigS3_iqr_models",
     ["iqr_models_detailed_results.csv"], "Per-IQR coefficients and confidence intervals."),
    ("Supplementary Figure S4", "supplementary", "FigS4_counterfactual_age",
     ["FigS4_counterfactual_age_aggregates.csv"], "Age-group box/whisker summaries; individual outliers omitted."),
    ("Supplementary Figure S5", "supplementary", "FigS5_lightgbm_scenarios",
     ["s13_scenario_reductions_by_age_pattern_v4.csv", "s14_scenario_reduction_distribution_v4.csv"],
     "Retained 5,914-observation scenario means, confidence intervals and quartiles from Tables S13-S14."),
    ("Supplementary Figure S6", "supplementary", "FigS6_dose_response_qgcomp",
     ["FigS6_gam_curve_coordinates.csv", "FigS6_gam_metadata.csv", "FigS6_signed_weights.csv", "FigS6_qgcomp_metadata.csv"],
     "Baseline age/sex-adjusted GAM curves and separately estimated signed exposure weights."),
    ("Supplementary Figure S7", "supplementary", "FigS7_phenoage_calibration",
     ["FigS7_calibration_deciles.csv", "FigS7_calibration_smoother.csv", "FigS7_calibration_annotations.csv"],
     "Risk-decile calibration and archived smoother; apparent validation, not external validation."),
    ("Supplementary Figure S8", "supplementary", "FigS8_lightgbm_nested",
     ["FigS8_lightgbm_nested_shap_groups.csv", "lightgbm_nested_importance_v5.csv"],
     "Aggregate grouped SHAP and permutation importance; no row-level SHAP values."),
    ("Supplementary Figure S9", "supplementary", "FigS9_joint_contrasts_precision",
     ["Figure3b_contrasts.csv", "Figure3b_DiD_precision.csv"],
     "Complete joint-group and within-stratum contrasts, DiD interaction, MDE and TOST; formerly main Figure 3c,d."),
]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def file_record(path: Path, root: Path) -> dict:
    return {"path": Path(os.path.relpath(path, root)).as_posix(),
            "sha256": digest(path), "bytes": path.stat().st_size}


def audit(root: Path, main_output: Path, si_output: Path, table_output: Path) -> dict:
    figures, tables, errors = [], [], []
    for item, folder, stem, names, note in FIGURES:
        inputs = [root / "03_source_data" / name for name in names]
        generated = (main_output if folder == "main" else si_output) / (stem + ".png")
        archived = root / "02_figures" / folder / (stem + ".png")
        missing = [Path(os.path.relpath(p, root)).as_posix() for p in inputs + [generated] if not p.is_file()]
        record = {"item": item, "notes": note, "model_refitting": False,
                  "sources": [file_record(p, root) for p in inputs if p.is_file()]}
        if missing:
            record.update(state="NOT_REPRODUCIBLE_FROM_RELEASE", missing=missing)
            errors.append(f"{item}: missing source or output")
        else:
            with Image.open(generated) as picture:
                picture.verify()
            record.update(state="REGENERATED_FROM_PUBLIC_SOURCE", generated=file_record(generated, root))
        if archived.is_file():
            record["archived"] = file_record(archived, root)
            if generated.is_file():
                record["byte_identical_to_archived"] = digest(generated) == digest(archived)
        figures.append(record)
    source = root / "03_source_data/tables"
    for metadata, rows in load_sources(source):
        item = metadata["table_id"]
        output = table_output / (item + ".html")
        record = {"item": item, "model_refitting": False,
                  "source": file_record(source / (item + ".csv"), root),
                  "notes": "Displayed archived table values reconstructed from public CSV/JSON."}
        if not output.is_file():
            record["state"] = "NOT_REPRODUCIBLE_FROM_RELEASE"
            errors.append(f"{item}: missing HTML")
        else:
            reader = CellReader()
            reader.feed(output.read_text(encoding="utf-8"))
            expected = {(item, ri, ci): value for ri, ci, value in display_cells(metadata, rows)}
            if reader.cells != expected:
                errors.append(f"{item}: generated cells differ from public source")
            record.update(state="REGENERATED_FROM_PUBLIC_SOURCE" if reader.cells == expected else "NOT_REPRODUCIBLE_FROM_RELEASE",
                          cells_match=reader.cells == expected, cells_verified=len(expected),
                          generated=file_record(output, root))
        tables.append(record)
    if len(figures) != 13 or len(tables) != 36:
        errors.append("Expected 13 figures and 36 physical tables")
    return {
        "status": "PASS" if not errors else "FAIL",
        "contract": "Reviewer reproduction uses only 03_source_data and 04_code. Raw environmental rasters and restricted participant data are not required.",
        "scope": "Aggregate-to-display reconstruction, not upstream data construction or statistical model refitting.",
        "states": STATES, "figures": figures, "tables": tables, "errors": errors,
        "summary": {"figure_counts_by_state": dict(Counter(r["state"] for r in figures)),
                    "table_counts_by_state": dict(Counter(r["state"] for r in tables)),
                    "table_cells_verified": sum(r.get("cells_verified", 0) for r in tables)},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("root", "main-output", "si-output", "table-output", "output"):
        parser.add_argument("--" + option, type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.root.resolve(), args.main_output.resolve(), args.si_output.resolve(), args.table_output.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], **report["summary"], "errors": report["errors"]}, indent=2))
    if report["errors"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
