# Python tools 🐍

For reviewer reconstruction of all 13 figures and 35 displayed table panels:

```shell
python 04_code/run_pipeline.py --mode public
```

Use Python 3.11 and `04_code/requirements_public.txt`. This path needs only
`03_source_data/` and `04_code/`, with no author-side source construction or
model refitting. See `04_code/PUBLIC_REPRODUCTION.md` for output paths and limits.

Use an existing conda environment. The data conversion, schema adapter, figure
composition and release-builder scripts resolve their package root from their
own path. The R/Stata model implementations remain in CODE to preserve the
reported estimation methods and weight/variance definitions.

From the package root:

```shell
python 04_code/PYTHON/prepare_analysis_data.py --input authorised_input.dta --output 04_code/DATA/analysis_ready/analysis_dataset.dta --audit audit/data_conversion.json
python 04_code/PYTHON/export_canonical_data.py
python 04_code/PYTHON/compose_main_figures.py --root . --out audit/main_figures
python 04_code/PYTHON/render_primary_figure.py --root . --out audit/main_figures
python 04_code/PYTHON/render_joint_sensitivity.py --root . --out audit/main_figures --si-out audit/supplementary_figures
python 04_code/PYTHON/build_github_release.py
```

The conversion preserves values and missingness; it does not reconstruct the
missing upstream biological-age calibration or environmental linkage. Both
internal and canonical microdata are restricted and excluded from public builds.
The builder requires `--replace` to replace a previous generated release.

Supplementary plotting needs the separate optional plotting dependencies:

```shell
python 04_code/PYTHON/render_supplementary_figures.py --root . --output-dir audit/supplementary_figures
```

For authorised local copies of the manuscript and SI (not distributed in this
code-only public package), save validation reports explicitly outside source files:

```shell
python 04_code/PYTHON/verify_submission_references.py path/to/manuscript.docx --output audit/references.json
python 04_code/PYTHON/verify_manuscript_si.py --root . --manuscript path/to/manuscript.docx --si path/to/si.docx --output audit/manuscript_si.json
```

Historical document-editing utilities in the local author workspace are not
part of the public analysis allowlist. On the audit workstation all temporary
and verification outputs were directed to the task-specific temporary folder.

## 📊 Current figure modules

`render_primary_figure.py` creates the four-panel Figure 2 from exposure-bin
counts, four unweighted dual-model estimates, and descriptive density/box/mean
summaries without mean confidence intervals. `render_joint_sensitivity.py`
creates Figure 3 from adjusted means and all 18 continuous-model estimates,
and Figure S9 from the complete contrasts and DiD precision. Their outputs
include editable PDF/SVG and 600-dpi PNG/TIFF. No model is refitted or selected
by P value; the inference threshold remains alpha = 0.05.

