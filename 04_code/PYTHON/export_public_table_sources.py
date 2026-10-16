"""Author-side export of manuscript tables to public, aggregate-only CSV.

The exported CSV files preserve the manuscript's displayed strings, including
rounding and uncertainty notation. Linked aggregate result files retain their
original precision. This is a display/provenance export, not model refitting.
Reviewers use render_public_tables.py; they do not need the author DOCX files.
Only the Python standard library is required.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[2]
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
W = "{" + NS["w"] + "}"
TITLE_RE = re.compile(r"^(?:Supplementary )?Table\s+(S?\d+)\b", re.I)
CJK_RE = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff]")
NUMBER_RE = re.compile(r"[−+\-]?\d+(?:,\d{3})*(?:\.\d+)?(?:[eE][+\-]?\d+)?")
PRIVATE_HEADERS = {
    "id", "id_num", "participant_id", "participant_id_numeric", "charls_id",
    "harm_id", "w5_id", "merge_id", "lat", "lon", "latitude", "longitude",
}

# These are related aggregate analyses, not a claim that all displayed cells
# are generated directly from the linked files. Missing CSVs fail explicitly.
SOURCES = {
    "Table_1": ("02_descriptive_stats.do", []),
    "Table_S01A": ("12_quartile_scores.R", ["quartile_scores_detailed_results.csv"]),
    "Table_S01B": ("15_qgcomp_final.R", ["gcomp_core_exposures_summary.csv"]),
    "Table_S02": ("30_iv_diagnostics.R;31_iv_full_tables.R", ["iv_estimates_v4.csv"]),
    "Table_S03": ("32_mediation_analysis.do", ["mediation_analysis_results.csv"]),
    "Table_S04": ("03_cohort_flow_replication.R;41_attrition_analysis.R", ["s4_cohort_flow_v4.csv"]),
    "Table_S05": ("22_nonlinearity_tests.do", []),
    "Table_S06": ("19_dose_response.R", ["threshold_results_EVI_O3.csv"]),
    "Table_S07": ("23_interaction_analysis.do", []),
    "Table_S08": ("24_logistic_sensitivity.do", []),
    "Table_S09": ("16_iqr_models.R", ["iqr_models_detailed_results.csv"]),
    "Table_S10": ("25_delta_baa_sensitivity.do", []),
    "Table_S11": ("30_iv_diagnostics.R;31_iv_full_tables.R", ["iv_validity_tests_v4.csv"]),
    "Table_S12": ("18_qgcomp_heterogeneity.R", []),
    "Table_S13": ("42_predictive_scenario_si_tables.R", ["s13_scenario_reductions_by_age_pattern_v4.csv"]),
    "Table_S14": ("42_predictive_scenario_si_tables.R", ["s14_scenario_reduction_distribution_v4.csv"]),
    "Table_S15": ("05_adjusted_baa_ba.do", ["detailed_results_summary.csv"]),
    "Table_S16": ("08_joint_association.R", ["joint_association_contrasts_v4.csv"]),
    "Table_S17": ("36_predictive_scenarios.R", ["predictive_model_performance_v4.csv", "predictive_scenario_contrasts_v4.csv"]),
    "Table_S18": ("41_attrition_analysis.R", ["s18_excluded_vs_included_v4.csv"]),
    "Table_S19": ("", []),
    "Table_S20": ("26_kmeans_robustness.do;27_kmeans_unweighted_verification.R", []),
    "Table_S21": ("40_clinical_benchmarks.R;43_primary_continuous_models.R", ["s21_clinical_benchmarks_v4.csv", "primary_continuous_models_v5.csv"]),
    "Table_S22": ("38_phenoage_validation.R", ["PhenoAge_validation_metrics.csv"]),
    "Table_S23A": ("41_attrition_analysis.R", ["attrition_baseline_comparison_v4.csv"]),
    "Table_S23B": ("41_attrition_analysis.R", ["attrition_ipaw_joint_v4.csv"]),
    "Table_S23C": ("41_attrition_analysis.R", ["attrition_pattern_mixture_v4.csv"]),
    "Table_S24": ("43_primary_continuous_models.R", ["primary_continuous_models_v5.csv"]),
    "Table_S25": ("44_outcome_validation_splitsample.R", ["outcome_validation_splitsample_v5.csv"]),
    "Table_S26": ("45_mde_tost.R", ["did_mde_tost_v5.csv"]),
    "Table_S27": ("46_selection_exposure_table.R", ["selection_exposure_by_step_v5.csv"]),
    "Table_S28": ("47_mi_sensitivity.R", ["mi_sensitivity_v5.csv"]),
    "Table_S29": ("48_weight_diagnostics.R", ["weight_distribution_v5.csv", "weight_balance_smd_v5.csv"]),
    "Table_S30A": ("49_lightgbm_nested.R", ["lightgbm_nested_performance_v5.csv"]),
    "Table_S30B": ("49_lightgbm_nested.R", ["lightgbm_nested_shap_groups_v5.csv"]),
    "Table_S31": ("55_copollutant_sensitivity.R", ["copollutant_sensitivity_v1.csv"]),
}
AUTHOR_EVIDENCE = {
    "Table_1": ["06_tables/descriptive_stats_baseline.xlsx"],
    "Table_S05": ["06_tables/Nonlinearity_Test_Results_Formatted.docx"],
    "Table_S07": ["06_tables/EVI_O3_interaction_summary_professional.docx"],
    "Table_S10": ["06_tables/results_summary_DeltaBAA.docx"],
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def plain_text(element: ET.Element) -> str:
    """Read visible runs, retaining explicit breaks without reading field code."""
    pieces = []
    for node in element.iter():
        if node.tag == W + "t":
            pieces.append(node.text or "")
        elif node.tag == W + "tab":
            pieces.append("\t")
        elif node.tag in {W + "br", W + "cr"} and node.get(W + "type") != "page":
            pieces.append("\n")
    return "".join(pieces).strip()


def paragraph_texts(elements: list[ET.Element]) -> list[str]:
    return [plain_text(p) for p in elements if p.tag == W + "p" and plain_text(p)]


def extract_grid(table: ET.Element) -> tuple[list[list[str]], list[list[dict]], int]:
    """Keep raw Word cells and explicit merged-cell structure in a rectangular CSV."""
    ncols = len(table.findall("w:tblGrid/w:gridCol", NS))
    rows, layout = [], []
    header_count = 0
    for row_no, row in enumerate(table.findall("w:tr", NS)):
        before = row.find("w:trPr/w:gridBefore", NS)
        col = int(before.get(W + "val", "0")) if before is not None else 0
        data = [""] * ncols
        cells = []
        for cell in row.findall("w:tc", NS):
            if cell.find(".//w:tbl", NS) is not None:
                raise ValueError("Nested table needs an explicit export policy")
            span_element = cell.find("w:tcPr/w:gridSpan", NS)
            span = int(span_element.get(W + "val", "1")) if span_element is not None else 1
            value = "\n".join(plain_text(p) for p in cell.findall("w:p", NS)).strip()
            if col + span > ncols:
                raise ValueError("Cell extends beyond the declared Word table grid")
            merge = cell.find("w:tcPr/w:vMerge", NS)
            merge_value = merge.get(W + "val", "continue") if merge is not None else None
            data[col] = value
            cells.append({"column": col, "colspan": span, "vertical_merge": merge_value})
            col += span
        repeat = row.find("w:trPr/w:tblHeader", NS)
        if repeat is not None and repeat.get(W + "val", "1") not in {"0", "false"}:
            if header_count != row_no:
                raise ValueError("Non-contiguous repeated table headers")
            header_count += 1
        rows.append(data)
        layout.append(cells)
    if not rows or not ncols:
        raise ValueError("Empty table")
    return rows, layout, header_count or 1


def extract_tables(path: Path, kind: str) -> list[dict]:
    with ZipFile(path) as archive:
        if archive.testzip() is not None:
            raise ValueError(f"Invalid DOCX archive: {path.name}")
        document = ET.fromstring(archive.read("word/document.xml"))
    if document.findall(".//w:ins", NS) or document.findall(".//w:del", NS):
        raise ValueError("Accept or resolve tracked changes before exporting tables")
    body = document.find("w:body", NS)
    if body is None:
        raise ValueError("Missing document body")
    blocks = list(body)
    table_positions = [i for i, b in enumerate(blocks) if b.tag == W + "tbl"]
    extracted = []
    current_title = ""
    current_number = ""
    last_title_position = -1
    for physical_index, position in enumerate(table_positions, 1):
        for i in range(position):
            if blocks[i].tag == W + "p":
                value = plain_text(blocks[i])
                match = TITLE_RE.match(value)
                if match:
                    current_title, current_number = value, match.group(1).upper()
                    last_title_position = i
        if not current_title or not current_number:
            raise ValueError(f"Table {physical_index} in {path.name} has no caption")
        previous_table = table_positions[physical_index - 2] if physical_index > 1 else -1
        prefix_blocks = blocks[max(last_title_position + 1, previous_table + 1):position]
        prefix = paragraph_texts(prefix_blocks)
        # The table's caption may include the panel label, or a separate paragraph may do so.
        panel_match = re.search(r"\bPanel\s+([A-Za-z])\b", current_title, re.I)
        panel = panel_match.group(1).upper() if panel_match else ""
        panel_title = ""
        for i, value in enumerate(prefix):
            match = re.match(r"^Panel\s+([A-Za-z])\b", value, re.I)
            if match:
                panel, panel_title = match.group(1).upper(), value
                prefix = prefix[i + 1:]
                break
        if kind == "main":
            table_id = "Table_" + current_number
        else:
            table_id = f"Table_S{int(current_number[1:]):02d}{panel}"
        following = []
        for element in blocks[position + 1:]:
            if element.tag != W + "p":
                if element.tag == W + "tbl":
                    break
                continue
            value = plain_text(element)
            if not value:
                # Images delimit a table's following notes in the main manuscript.
                if element.find(".//w:drawing", NS) is not None or element.find(".//w:pict", NS) is not None:
                    break
                continue
            if TITLE_RE.match(value) or re.match(r"^(?:Panel\s+[A-Za-z]\b|Supplementary Note|Figure\s+\d)", value):
                break
            if kind == "main" and following and not value.startswith("Abbreviations:"):
                break
            following.append(value)
        data, layout, headers = extract_grid(blocks[position])
        extracted.append({
            "table_id": table_id, "display_number": current_number,
            "title": current_title, "panel": panel, "panel_title": panel_title,
            "notes_before": prefix, "notes_after": following,
            "document_kind": kind, "document_name": path.name,
            "document_sha256": sha256(path), "physical_table_index": physical_index,
            "body_element_index": position, "rows": data, "cell_layout": layout,
            "header_row_count": headers, "column_count": len(data[0]),
        })
    # S1 and S30 place their shared explanation after the last physical panel.
    for number in {t["display_number"] for t in extracted}:
        group = [t for t in extracted if t["display_number"] == number]
        if len(group) > 1:
            for table in group:
                table["shared_notes"] = list(group[-1]["notes_after"]) if not table["notes_after"] else []
    return extracted


def check_english(value: object, location: str) -> None:
    if CJK_RE.search(json.dumps(value, ensure_ascii=False)):
        raise ValueError(f"Non-English CJK content found in {location}")


def locate_aggregate(root: Path, basename: str) -> Path:
    for subdir in ("03_source_data", "06_tables", "04_code/output/tables"):
        path = root / subdir / basename
        if path.is_file():
            return path
    raise FileNotFoundError(f"Missing allowlisted aggregate result: {basename}")


def link_aggregate(root: Path, output: Path, basename: str) -> dict:
    source = locate_aggregate(root, basename)
    with source.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.reader(stream))
    if not rows or len(rows) > 250:
        raise ValueError(f"Unexpected aggregate size: {basename}")
    if any(h.strip().lower() in PRIVATE_HEADERS for h in rows[0]):
        raise ValueError(f"Identifier or individual-coordinate header in {basename}")
    check_english(rows, basename)
    # Keep source bytes and full precision; no automatic numeric conversion.
    destination = output / "archived_results" / basename
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return {
        "public_path": destination.relative_to(root).as_posix(),
        "author_archive_path": source.relative_to(root).as_posix(),
        "sha256": sha256(destination), "row_count": len(rows) - 1,
        "relationship": "related_aggregate_result; not asserted to derive every displayed cell",
    }


def read_display_csv(path: Path) -> list[list[str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.reader(stream))
    return rows[1:]


def export(root: Path, output: Path, report_dir: Path, crosswalk: Path) -> dict:
    inputs = [
        (root / "01_manuscript/Zhang_main_manuscript.docx", "main", 1),
        (root / "01_manuscript/Zhang_supplementary_information.docx", "supplementary", 35),
    ]
    tables = []
    for path, kind, expected in inputs:
        current = extract_tables(path, kind)
        if len(current) != expected:
            raise ValueError(f"Expected {expected} tables in {path.name}; found {len(current)}")
        tables.extend(current)
    if set(t["table_id"] for t in tables) != set(SOURCES):
        raise ValueError("Extracted table IDs do not match the explicit 36-table inventory")
    if len({t["table_id"] for t in tables}) != len(tables):
        raise ValueError("Duplicate table IDs")
    output.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    crosswalk.parent.mkdir(parents=True, exist_ok=True)
    source_cache = {}
    manifest_tables, verification, crosswalk_rows = [], [], []
    for source_table in tables:
        table = dict(source_table)
        table_id = table["table_id"]
        check_english(table, table_id)
        data = table.pop("rows")
        csv_path = output / (table_id + ".csv")
        with csv_path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow([f"column_{i + 1:02d}" for i in range(table["column_count"])])
            writer.writerows(data)
        restored = read_display_csv(csv_path)
        if data != restored:
            raise AssertionError(f"CSV round-trip mismatch: {table_id}")
        # Numeric comparison is cell-specific, not a permissive document-wide bag.
        numeric_cells = sum(bool(NUMBER_RE.search(v)) for row in data for v in row)
        numeric_tokens = sum(len(NUMBER_RE.findall(v)) for row in data for v in row)
        for ri, row in enumerate(data):
            for ci, value in enumerate(row):
                if NUMBER_RE.findall(value) != NUMBER_RE.findall(restored[ri][ci]):
                    raise AssertionError(f"Numeric cell mismatch: {table_id} row {ri + 1} col {ci + 1}")
        script_names, aggregate_names = SOURCES[table_id]
        linked = []
        for basename in aggregate_names:
            if basename not in source_cache:
                source_cache[basename] = link_aggregate(root, output, basename)
            linked.append(source_cache[basename])
        table.update({
            "schema_version": 1,
            "csv_file": csv_path.name, "csv_sha256": sha256(csv_path),
            "row_count_including_display_headers": len(data),
            "csv_encoding": "UTF-8", "values_are_strings": True,
            "csv_header_role": "storage column keys; displayed headers are the first header_row_count data rows",
            "export_method": "documented_displayed_values",
            "provenance_status": "archived_aggregate_sources_linked" if linked else "documented_displayed_values",
            "inference_scope": "Reconstruction of the authors' archived displayed results; no model refitting or independent verification of estimates.",
            "aggregate_only": True,
            "related_analysis_scripts": ["04_code/CODE/" + s for s in script_names.split(";") if s],
            "linked_aggregate_results": linked,
            "other_author_evidence": AUTHOR_EVIDENCE.get(table_id, []),
        })
        meta_path = output / (table_id + ".json")
        write_json(meta_path, table)
        manifest_tables.append({
            "table_id": table_id, "csv_file": csv_path.name, "metadata_file": meta_path.name,
            "csv_sha256": sha256(csv_path), "metadata_sha256": sha256(meta_path),
            "title": table["title"], "panel_title": table["panel_title"],
        })
        verification.append({
            "table_id": table_id, "document": table["document_name"],
            "physical_table_index": table["physical_table_index"],
            "rows": len(data), "columns": table["column_count"],
            "rectangular_cells_checked": sum(len(row) for row in data),
            "physical_cells_checked": sum(len(row) for row in table["cell_layout"]),
            "numeric_cells_checked": numeric_cells, "numeric_tokens_checked": numeric_tokens,
            "text_cell_match": True, "numeric_cell_match": True,
            "notes_and_title_preserved": True,
            "source_document_sha256": table["document_sha256"],
        })
        crosswalk_rows.append({
            "table_id": table_id, "title": table["title"], "panel_title": table["panel_title"],
            "csv_path": csv_path.relative_to(root).as_posix(),
            "metadata_path": meta_path.relative_to(root).as_posix(),
            "renderer": "04_code/PYTHON/render_public_tables.py",
            "html_output": f"04_code/output/public_tables/{table_id}.html",
            "related_analysis_scripts": ";".join(table["related_analysis_scripts"]),
            "linked_public_aggregate_results": ";".join(x["public_path"] for x in linked),
            "provenance_status": table["provenance_status"],
            "export_method": "documented_displayed_values",
            "source_document": table["document_name"],
            "physical_table_index": table["physical_table_index"],
            "verification_scope": "all displayed labels and numerical strings; not model re-estimation",
        })
    # Detect concurrent manuscript edits rather than silently mixing versions.
    for path, _, _ in inputs:
        expected_hash = next(t["document_sha256"] for t in tables if t["document_name"] == path.name)
        if sha256(path) != expected_hash:
            raise RuntimeError(f"Source document changed during export: {path.name}")
    manifest = {
        "schema_version": 1, "table_count": len(tables), "main_table_count": 1,
        "supplementary_physical_table_count": 35,
        "data_scope": "published table cells and aggregate model summaries only; no participant records",
        "reconstruction_scope": "archived_displayed_results_not_model_refitting",
        "tables": manifest_tables,
    }
    write_json(output / "manifest.json", manifest)
    with crosswalk.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(crosswalk_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(crosswalk_rows)
    totals = {
        "status": "PASS", "tables_verified": len(verification),
        "physical_cells_checked": sum(v["physical_cells_checked"] for v in verification),
        "numeric_cells_checked": sum(v["numeric_cells_checked"] for v in verification),
        "numeric_tokens_checked": sum(v["numeric_tokens_checked"] for v in verification),
        "linked_aggregate_csv_files": len(source_cache),
        "provenance_status_counts": dict(Counter(x["provenance_status"] for x in crosswalk_rows)),
        "tables": verification,
        "boundary": "Exact table-display export does not validate model identification, source-outcome derivation, or original parameter estimation.",
    }
    write_json(report_dir / "table_export_verification.json", totals)
    (output / "README.md").write_text(
        "# Public table sources 📊\n\n"
        "This directory contains Table 1 and all 35 supplementary table panels "
        "as aggregate-only CSV with per-table JSON metadata. No participant records are included.\n\n"
        "## Rebuild locally 🛠️\n\n"
        "From the repository root, run:\n\n"
        "```text\npython 04_code/PYTHON/render_public_tables.py\n```\n\n"
        "The renderer uses only the Python standard library and this directory. "
        "It writes individual HTML files, a combined HTML document, a searchable index, "
        "and CSV copies to `04_code/output/public_tables/`. The author DOCX files are not required. "
        "Use `--output` to choose another output directory or `--xlsx` to add a workbook "
        "when openpyxl is installed. Paths are resolved from the script location or `--root`.\n\n"
        "## Meaning of the files 🔎\n\n"
        "CSV cells are strings so trailing zeroes, confidence intervals and significance "
        "markers retain the manuscript's displayed precision. Generic `column_01` keys "
        "are storage headers; displayed table headers occupy the first data row(s). "
        "JSON records titles, pre-table and post-table notes, shared panel notes, merged "
        "cells, header rows, source-document hashes and related analysis scripts.\n\n"
        "The author-side exporter is `export_public_table_sources.py`; it reads DOCX "
        "only when the authors refresh the public sources. Each export is verified "
        "cell by cell against the source document. `archived_results/` contains "
        "allowlisted related aggregate CSVs at their archived precision. Linked results "
        "do not imply that every displayed cell was derived anew from those files. "
        "Tables with no archived aggregate CSV link are explicitly marked "
        "`documented_displayed_values`.\n\n"
        "## Interpretation 🧭\n\n"
        "This is reconstruction of the authors' archived result tables. It does not "
        "refit models or independently validate estimates. All methodological limitations "
        "and denominator qualifications in the manuscript table notes are retained.\n",
        encoding="utf-8",
    )
    return totals


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path, help="Default: ROOT/03_source_data/tables")
    parser.add_argument("--report-dir", type=Path, help="Default: ROOT/04_code/output/table_source_export")
    parser.add_argument("--crosswalk", type=Path, help="Default: ROOT/05_documentation/public_table_crosswalk.csv")
    args = parser.parse_args()
    root = args.root.resolve()
    def destination(value: Path | None, default: str) -> Path:
        return (value if value is not None and value.is_absolute() else root / (value or default)).resolve()
    result = export(
        root, destination(args.output, "03_source_data/tables"),
        destination(args.report_dir, "04_code/output/table_source_export"),
        destination(args.crosswalk, "05_documentation/public_table_crosswalk.csv"),
    )
    print(json.dumps({k: v for k, v in result.items() if k != "tables"}, indent=2))


if __name__ == "__main__":
    main()
