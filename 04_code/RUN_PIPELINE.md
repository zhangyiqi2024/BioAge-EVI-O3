# Pipeline runner ▶️

Reviewer entry point (default mode):

```shell
python path/to/04_code/run_pipeline.py --mode public
```

Install `requirements_public.txt` in a Python 3.11 environment first. Inputs
are only the sibling `03_source_data/` folder and this `04_code/` folder.
All 4 main figures, 9 SI figures, 1 main table and 35 SI table panels are
rebuilt from public aggregates. No original raster or participant record is
required. This is displayed-result reconstruction, not statistical refitting.

Outputs: `04_code/output/main/`, `supplementary/`, `tables/` and `logs/`.
Use `--output-dir PATH` to redirect all public outputs. Every failure exits
nonzero; inspect `logs/pipeline_result.json`, `session_info.json` and
`public_reproduction_status.json`. The runner overwrites stale success status
at startup and records failure if a step fails.

See `PUBLIC_REPRODUCTION.md` for source provenance, output reading and scope.

## Optional restricted-data analysis

```shell
python path/to/04_code/run_pipeline.py --mode analysis --stata StataMP-64.exe --rscript Rscript
python path/to/04_code/run_pipeline.py --mode check --rscript Rscript
```

Analysis mode requires `DATA/analysis_ready/analysis_dataset.dta`. It checks R
requirements, runs the Stata master, verifies joint weights and runs R steps.
Optional restricted inputs produce explicit SKIP entries. See DATA/BUILD_SCOPE.md
and LOCAL_VERIFICATION.md for the unresolved upstream BAA construction and local
R/UCRT exit error. Analysis-mode logs are under `04_code/logs/`.

The 00/01 helpers are not the complete upstream outcome/geospatial workflow.
They write candidate inputs and never overwrite the canonical dataset.

## 📊 Figure renderers

`PYTHON/render_primary_figure.py` renders the four-panel Figure 2.
`PYTHON/render_joint_sensitivity.py` renders Figure 3 and Figure S9. Figures
2–3 and S9 include editable PDF/SVG and 600-dpi PNG/TIFF; Figure 4 also has
vector SVG. The complete pipeline has 13 figures and 36 physical table panels.
See `RUN_FIGURE3.md` for model coverage and statistical scope.

