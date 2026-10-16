# Local verification record

## Current public acceptance: 2026-10-13

Only `03_source_data` and `04_code` were copied into an empty test directory.
The entry point was launched from a different working directory using Python
3.11.15, Matplotlib 3.11.1, pandas 2.3.3, NumPy 2.4.6 and Pillow 12.3.0.
All six pipeline steps exited zero in about 33 seconds. All 12 generated PNGs
were byte-identical to the revised formal figure files. All 36 table panels
and 1,916 displayed cells passed source-to-HTML verification. The copy contained
no DOCX, original rasters, individual-level data, RDS or pre-existing outputs.
Removing the required S4 aggregate source produced the expected failure;
the renderer cannot fall back to restricted data.

This certifies aggregate-to-display reconstruction. It does not certify
statistical refitting, primary BAA construction or the exact original ozone
input version. The archived R/UCRT issue below is outside the public path.

## Earlier restricted-data and partial-render checks

Audit date: 2026-10-12, Windows, R 4.4.3, Stata 18 MP, existing conda Python.

- English metadata conversion preserves all 6,568 rows, 205 variables and
  numerical values/missingness exactly. Original DTA is retained privately.
- Canonical English aliases are one-to-one and reversible, with exact DTA
  numeric round-trip assertions. Codebook contains no participant examples.
- The three archived joint propensity weights have zero maximum absolute
  difference from the specified nnet reconstruction, including missingness.
- The complete Stata analysis master and repaired wrapper completed in an
  independent scratch copy with recorded Stata status 0. Restricted mortality
  inputs were absent and that optional step was skipped explicitly.
- The 00/01 preparation helpers ran successfully on a copied authorised
  analysis-level input; they produced a separate candidate and Stata status 0.
  This smoke test does not establish raw biomarker/exposure reconstruction.
- The public Python path ran from a copied release with no restricted data.
  All three rebuilt main-figure PNGs were byte-identical to their archives.
  The release scan found no DTA/RDS/cache/log files and no CJK in code.
- Supplementary S2-S4 rendered successfully using only three aggregate CSVs
  in an isolated no-DTA copy. S4's ten group summaries preserve all box/whisker
  statistics but omit individual outlier points, as labelled on the public plot.
- A full clean R pipeline is **not certified**. Loading `cli` or `rlang` alone
  (without data or model code) causes process exit 0xc0000005 in Windows
  `ucrtbase.dll` on this workstation. Rscript/Rterm 4.4.1/4.4.2/4.4.3, --vanilla,
  independent input copies and clean locale/PATH probes reproduce the issue.
  Plain R, Rcpp, nnet and glue exit normally. The runner deliberately fails on
  nonzero exits and does not reinterpret a written output as overall success.
- Independently completed joint-model and primary continuous-model computations
  reproduce the archived means, contrasts and all 16 primary/sensitivity rows.
  Their R processes nevertheless fail on exit. MDE/TOST and the corrected S21
  outputs also match their source specifications. These numerical checks do
  not amount to an end-to-end successful R run.
- Optional raw-data-dependent analyses are not newly reproduced when their
  documented restricted inputs are absent.

The local audit folder (excluded from GitHub) records commands, process exit
codes, comparisons and runtime probes. A public rendering smoke test is separate
from statistical refitting; it cannot validate unavailable outcome calibration.
