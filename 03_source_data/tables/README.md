# Public table sources 📊

This directory contains Table 1 and all 34 supplementary table panels as aggregate-only CSV with per-table JSON metadata. No participant records are included.

## Rebuild locally 🛠️

From the repository root, run:

```text
python 04_code/PYTHON/render_public_tables.py
```

The renderer uses only the Python standard library and this directory. It writes individual HTML files, a combined HTML document, a searchable index, and CSV copies to `04_code/output/public_tables/`. The author DOCX files are not required. Use `--output` to choose another output directory or `--xlsx` to add a workbook when openpyxl is installed. Paths are resolved from the script location or `--root`.

## Meaning of the files 🔎

CSV cells are strings so trailing zeroes, confidence intervals and significance markers retain the manuscript's displayed precision. Generic `column_01` keys are storage headers; displayed table headers occupy the first data row(s). JSON records titles, pre-table and post-table notes, shared panel notes, merged cells, header rows, source-document hashes and related analysis scripts.

The author-side exporter is `export_public_table_sources.py`; it reads DOCX only when the authors refresh the public sources. Each export is verified cell by cell against the source document. `archived_results/` contains allowlisted related aggregate CSVs at their archived precision. Linked results do not imply that every displayed cell was derived anew from those files. Tables with no archived aggregate CSV link are explicitly marked `documented_displayed_values`.

## Interpretation 🧭

This is reconstruction of the authors' archived result tables. It does not refit models or independently validate estimates. All methodological limitations and denominator qualifications in the manuscript table notes are retained.
