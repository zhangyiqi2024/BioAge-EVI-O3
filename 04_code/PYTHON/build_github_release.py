"""Build an allowlisted public package without restricted microdata.

Run from any working directory. Existing generated releases require --replace.
The builder does not publish, upload or change the analysis input.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[2]
RESTRICTED_HEADERS = {
    "id", "id_num", "participant_id", "participant_id_numeric", "charls_id",
    "harm_id", "w5_id", "merge_id", "latitude", "longitude", "lat", "lon",
}
AGGREGATE_CSV = {
    "03_source_data": """
        detailed_results_summary.csv FigS8_lightgbm_nested_performance.csv
        FigS8_lightgbm_nested_shap_groups.csv Figure1_observed_distribution.csv
        Figure2_continuous_forest.csv Figure3a_adjusted_means.csv
        Figure3b_contrasts.csv Figure3b_DiD_precision.csv
        Figure4a_scenario_contrasts.csv Figure4_exposure_response_v5.csv
        Figure4_joint_contour_v5.csv Figure4_model_performance.csv
        FigureS1_participant_flow.csv gcomp_core_exposures_summary.csv
        iqr_models_detailed_results.csv mediation_analysis_results.csv
        ml_counterfactual_diagnostics.csv ml_extrapolation_analysis.csv
        PhenoAge_validation_metrics.csv FigS4_counterfactual_age_aggregates.csv
    """.split(),
    "06_tables": """
        attrition_baseline_comparison_v4.csv attrition_ipaw_joint_v4.csv
        attrition_pattern_mixture_v4.csv bmi_robustness_v5.csv
        Complete_Results_with_Meta_Final.csv Detailed_Results_Final.csv
        detailed_results_summary.csv did_mde_tost_v5.csv
        dual_exposure_iv_results_pca.csv Exposure_Effects_Summary.csv
        gcomp_core_exposures_summary.csv gcomp_core_exposures_summary_1KM_500m.csv
        iqr_models_detailed_results.csv iqr_two_scales_detailed_results.csv
        iv_estimates_v4.csv iv_validity_tests_v4.csv
        joint_association_adjusted_means_v4.csv joint_association_contrasts_v4.csv
        lightgbm_nested_counterfactual_v5.csv lightgbm_nested_importance_v5.csv
        lightgbm_nested_performance_v5.csv lightgbm_nested_shap_features_v5.csv
        lightgbm_nested_shap_groups_v5.csv mediation_analysis_results.csv
        mi_sensitivity_v5.csv outcome_validation_splitsample_v5.csv
        predictive_model_performance_v4.csv predictive_scenario_contrasts_v4.csv
        primary_continuous_models_v5.csv quartile_scores_detailed_results.csv
        s13_scenario_reductions_by_age_pattern_v4.csv
        s14_scenario_reduction_distribution_v4.csv s18_excluded_vs_included_v4.csv
        s21_clinical_benchmarks_v4.csv s4_cohort_flow_v4.csv
        selection_exposure_by_step_v5.csv single_exposure_iv_results_pca.csv
        threshold_results_EVI_O3.csv weight_balance_smd_v5.csv weight_distribution_v5.csv
    """.split(),
}
MIT = """MIT License

Copyright (c) 2026 YiQi Zhang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Scope: this licence applies to original analysis code and documentation.
CHARLS microdata, environmental datasets, and third-party software remain
subject to their own licences and access agreements. No rights to redistribute
restricted participant records are granted by this licence.
"""
CITATION = """cff-version: 1.2.0
message: "If you use the code or aggregate source data, please cite the manuscript and this release."
title: "Code and aggregate source data for Repeated-measure associations of EVI within residential buffers and ozone with biological age acceleration in China"
type: software
authors:
  - family-names: Zhang
    given-names: YiQi
  - family-names: Ying
    given-names: Jun
  - family-names: Zhang
    given-names: Ruiwu
  - family-names: Zhang
    given-names: Youyou
repository-code: "https://github.com/zhangyiqi2024/BioAge-EVI-O3"
license: MIT
version: "2026.10.17-reviewer-materials"
date-released: "2026-10-17"
abstract: "Aggregate-only reproduction package for the manuscript Repeated-measure associations of EVI within residential buffers and ozone with biological age acceleration in China. The public Python path rebuilds the displayed main and supplementary figures and tables without distributing individual-level CHARLS records or environmental rasters."
"""
README = "# BioAge, EVI within residential buffers and ozone 🔬\n\nCode and aggregate results for **Repeated-measure associations of EVI within residential buffers and ozone with biological age acceleration in China**. The revision release branch is\n`release-2026-10-17-reviewer-materials` in `zhangyiqi2024/BioAge-EVI-O3`; use the Git commit and\n`RELEASE_MANIFEST.json` to identify a particular release. The historical Zenodo\narchive is not changed by running this pipeline.\n\n## ▶️ One-command reviewer reproduction\n\nUse Python 3.11 and the packages in `04_code/requirements_public.txt`:\n\n```shell\npython -m pip install -r 04_code/requirements_public.txt\npython 04_code/run_pipeline.py --mode public\n```\n\nThe public mode is also the default. The entry point discovers its own location,\nso it works when launched from any working directory. Only `03_source_data/`\nand `04_code/` are required. MODIS/ChinaHighO3 rasters, CHARLS participant records,\nR, Stata, Word, manuscript files and archived figure images are not required.\nThe two directories must remain siblings. No data download or network access\noccurs during reproduction.\n\n## 📊 Outputs and scope\n\n- `04_code/output/main/`: main Figures 1–4, PNG/TIFF/PDF; Figures 2–4 also have vector SVG.\n- `04_code/output/supplementary/`: Supplementary Figures S1–S9; S9 also has editable PDF/SVG.\n- `04_code/output/tables/`: Table 1 and 35 SI table panels, as HTML/CSV/JSON;\n  open `index.html` or print `all_tables.html` to PDF.\n- `04_code/output/logs/`: per-step logs, software versions, run status and\n  `public_reproduction_status.json` with source/output hashes and cell checks.\n\nAn alternative output directory can be set with `--output-dir PATH`.\nExcel is optional: install openpyxl and run\n`python 04_code/PYTHON/render_public_tables.py --xlsx`.\n\nThis path reconstructs the displayed results from disclosed aggregate source\ndata. It does not refit participant-level models or independently validate the\narchived estimates. Curves use evaluated prediction grids; distribution panels\nuse density/box summaries and binned exposures. Individual observations and\noutliers are not distributed. Table CSV/JSON preserve displayed precision and\nfootnotes; each has provenance distinguishing an archived result file from a\nvalue documented in the manuscript. This distinction also applies to the eight\ntables lacking a separate machine-readable historical result file.\n\nS5 now uses the retained 5,914-observation summaries in Tables S13–S14. S6\nexplicitly distinguishes its 2,785-person baseline age/sex-adjusted GAM curves\nfrom its historical 6,002-observation quantile g-computation sample (ages 45–94).\nNeither is the primary 5,914-observation model.\n\n## 📂 Data and code guide\n\n`03_source_data/README.md` describes the English source dictionaries and\nprovenance records. `03_source_data/tables/manifest.json` lists all physical\ntables and checksums. `05_documentation/public_table_crosswalk.csv` and\n`04_code/PYTHON/audit_public_reproduction.py` map paper items to sources.\n`RELEASE_MANIFEST.json` records every release file. Files in `02_figures/` and\n`06_tables/` are comparison archives, not public-run dependencies.\n\nOptional restricted-data analysis is separate:\n\n```shell\npython 04_code/run_pipeline.py --mode analysis --stata StataMP-64.exe --rscript Rscript\n```\n\nIt requires independently authorised inputs specified in\n`04_code/DATA/BUILD_SCOPE.md`. The public mode uses only the supplied aggregate\nfiles and does not refit participant-level models or redistribute restricted records.\nR 4.4.3 on the audit workstation has a documented UCRT exit failure; the separate\npublic Python path is unaffected. See `04_code/LOCAL_VERIFICATION.md`.\n\nThe LightGBM split is participant-level, but initial imputation and climate PCA\nprecede splitting. The historical `nested` filenames do not establish fully\nnested validation. MDE is a power quantity, not an exclusion bound.\n\n## 🔐 Licence and data access\n\nCode/documentation: MIT (see LICENSE). Public tables/figure coordinates are\nauthor-provided aggregate study results; cite this study when reusing them.\nThird-party data retain their own licences. No individual-level DTA/CSV,\nidentifiers, row-level predictions or SHAP values are included. Obtain CHARLS\nmicrodata under its own access agreement. Use this curated release directory\nfor upload, not the local parent workspace containing restricted data.\n\nThe structure follows the single-entry-point and figure/table crosswalk\nconventions in Yusaku Horiuchi's replication-package-guide and FAIR principles.\n\nThe package reconstructs 13 figures and 36 physical table panels.\n\n## 📊 Current figure composition\n\nFigure 2 combines the 1-km joint exposure distribution and marginal histograms,\nfour unweighted dual-exposure estimates, and observed BAA distributions at both\nbuffers. Only the 1-km estimates are primary; 500-m estimates assess buffer\nsensitivity. Descriptive mean points have no naive mean confidence intervals;\nthe source interval columns remain unchanged. Suppressed counts of 1–4 remain\ndistinct from zero, so displayed cell counts need not sum to the model-frame\ntotal. Exposure coverage and continuous models use 5,914 observations; the\ndescriptive BAA panels use 6,472 observations per buffer.\n\nFigure 3 retains the adjusted joint-pattern means (6,471 observations) and all\n18 continuous-model sensitivity estimates: 16 non-MI specifications and both\n1-km MI rows. All nonsignificant IPW EVI results remain displayed. The non-MI\nmodels use robust variance and 5,914 observations from 3,154 participants;\nMI uses model-based variance and Rubin pooling across 20 imputations, with its\nsample counts unspecified in the aggregate table. Figure S9 contains all 12\njoint contrast rows and both DiD precision rows, including nominal P, MDE and\nTOST. These displays do not establish effect modification. No significance\nthreshold, estimate or nominal P value is changed.\n\nThe new renderers are `04_code/PYTHON/render_primary_figure.py` for Figure 2\nand `04_code/PYTHON/render_joint_sensitivity.py` for Figure 3 and Figure S9.\nThey produce editable PDF/SVG and 600-dpi PNG/TIFF from public aggregate inputs.\n\n"


README += ('\n## 📌 Fixed version\n\n'
           'This figure revision is release `2026.10.17-reviewer-materials`, preserved at '
           '[tag v2026.10.17-reviewer-materials](https://github.com/zhangyiqi2024/BioAge-EVI-O3/tree/v2026.10.17-reviewer-materials). '
           'Use the tag to retrieve the exact revision cited in the manuscript.\n')


def audit_aggregate(path: Path) -> None:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, [])
    forbidden = RESTRICTED_HEADERS.intersection(x.strip().lower() for x in header)
    if forbidden:
        raise ValueError(f"Restricted columns in {path.name}: {sorted(forbidden)}")
    if any(re.search(r"(?i)(participant_prediction|detailed_lightgbm|shap_values)", x)
           for x in [path.name]):
        raise ValueError(f"Row-level artifact is not allowed: {path.name}")


def build(root: Path, target: Path, replace: bool) -> dict:
    root, target = root.resolve(), target.resolve()
    if target == root or root not in target.parents:
        raise ValueError("Release output must be a strict child of the package root.")
    if target.exists():
        if not replace or not (target / "RELEASE_MANIFEST.json").is_file():
            raise ValueError("Existing output is not a replaceable generated release.")
        # Absolute descendant and generated-manifest checks precede deletion.
        shutil.rmtree(target)
    target.mkdir(parents=True)

    def copy(relative: str) -> None:
        src, dst = root / relative, target / relative
        if not src.is_file():
            raise FileNotFoundError(src)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)

    core = [
        "04_code/README.md", "04_code/LICENSE", "04_code/RUN_PIPELINE.md", "04_code/RUN_FIGURE3.md",
        "04_code/requirements.R", "04_code/requirements_r_versions.csv",
        "04_code/requirements_python.txt", "04_code/requirements_plotting.txt",
        "04_code/STATA_REQUIREMENTS.md", "04_code/run_pipeline.py",
        "04_code/LOCAL_VERIFICATION.md", "04_code/PUBLIC_REPRODUCTION.md", "04_code/DATA/BUILD_SCOPE.md",
        "04_code/DATA/source_name_mapping.json", "04_code/DATA/canonical_schema_crosswalk.csv",
        "04_code/DATA/analysis_ready/README.md", "04_code/DATA/analysis_ready/codebook.csv",
        "05_documentation/figure_crosswalk.csv", "05_documentation/si_crosswalk.csv",
        "05_documentation/public_table_crosswalk.csv", "05_documentation/public_figure_crosswalk.csv", "04_code/requirements_public.txt",
    ]
    for relative in core:
        copy(relative)
    for path in (root / "04_code/CODE").iterdir():
        if path.is_file() and path.suffix.lower() in {".r", ".do"}:
            copy(path.relative_to(root).as_posix())
    public_python = [
        "README.md", "prepare_analysis_data.py", "export_canonical_data.py",
        "compose_main_figures.py", "render_primary_figure.py",
        "render_joint_sensitivity.py", "render_supplementary_figures.py",
        "render_public_additional.py", "audit_public_reproduction.py",
        "render_public_distributions.py", "render_public_s5.py", "render_public_tables.py",
        "export_public_figure_sources.py", "export_public_table_sources.py",
        "build_github_release.py", "verify_submission_references.py",
        "verify_manuscript_si.py",
    ]
    for name in public_python:
        copy("04_code/PYTHON/" + name)
    for folder, names in AGGREGATE_CSV.items():
        for name in names:
            path = root / folder / name
            if not path.exists() and name == "FigS4_counterfactual_age_aggregates.csv":
                continue
            audit_aggregate(path)
            copy(path.relative_to(root).as_posix())
    # Every CSV in 03_source_data is an allowlisted aggregate source. Keeping
    # this directory self-contained means the public renderer never needs the
    # restricted analysis tree or a presentation-table duplicate.
    for path in sorted((root / "03_source_data").glob("*.csv")):
        audit_aggregate(path)
        copy(path.relative_to(root).as_posix())
    for path in sorted((root / "03_source_data").rglob("*")):
        if not path.is_file() or path.suffix.lower() not in {".csv", ".json", ".md"}:
            continue
        if path.suffix.lower() == ".csv":
            audit_aggregate(path)
        copy(path.relative_to(root).as_posix())
    # Reviewed aggregate-only presentation tables; never copy arbitrary binaries.
    for name in [
        "descriptive_stats_baseline.xlsx", "Detailed_Ordered_Results_Table.xlsx",
        "EVI_O3_interaction_summary_professional.docx", "IV_Detailed_Results.xlsx",
        "IV_Summary_Table_and_Figure.docx", "IV_Validity_Tests_Concise.docx",
        "IV_Validity_Tests_Concise_Data.xlsx", "Nonlinearity_Test_Results_Formatted.docx",
        "Ordered_Exposure_Effects_Summary_Table.xlsx", "results_summary_DeltaBAA.docx",
    ]:
        copy("06_tables/" + name)
    for folder, pattern in [("02_figures/main", "Figure*.*"),
                            ("02_figures/supplementary", "FigS*.*")]:
        for path in (root / folder).glob(pattern):
            if path.suffix.lower() in {".png", ".tiff", ".pdf", ".svg"}:
                copy(path.relative_to(root).as_posix())
    (target / "README.md").write_text(README, encoding="utf-8")
    (target / "LICENSE").write_text(MIT, encoding="utf-8")
    (target / "CITATION.cff").write_text(CITATION, encoding="utf-8")
    # Preserve audited bytes across Git checkouts on all operating systems.
    (target / ".gitattributes").write_text("* -text\n", encoding="utf-8")
    (target / ".gitignore").write_text(
        "*.dta\n*.rds\n*.RData\n*.Rhistory\n**/DATA/restricted/\n"
        "**/DATA/analysis_ready/*\n!**/DATA/analysis_ready/README.md\n"
        "!**/DATA/analysis_ready/codebook.csv\n"
        "**/DATA/canonical/\n**/logs/\n**/output/\n__pycache__/\n*.pyc\n.env\n",
        encoding="utf-8")
    files = []
    for path in sorted(target.rglob("*")):
        if path.is_file():
            files.append({"path": path.relative_to(target).as_posix(),
                          "bytes": path.stat().st_size,
                          "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    result = {
        "status": "GITHUB_RELEASE_BRANCH",
        "repository": "https://github.com/zhangyiqi2024/BioAge-EVI-O3",
        "branch": "release-2026-10-17-reviewer-materials",
        "version": "2026.10.17-reviewer-materials",
        "tag": "v2026.10.17-reviewer-materials",
        "files": files,
        "excluded": ["individual-level data", "DTA/RDS", "logs", "caches",
                     "row-level predictions and SHAP values", "local audit records"],
    }
    (target / "RELEASE_MANIFEST.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return {"output": str(target), "files": len(files), "bytes": sum(f["bytes"] for f in files)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()
    print(json.dumps(build(args.root, args.output or args.root / "07_github_release", args.replace), indent=2))


if __name__ == "__main__":
    main()
