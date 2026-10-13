# Current main Figure 3

The current public composite has four panels. Panels a-b show adjusted marginal mean BAA for the four joint exposure patterns at 1-km and 500-m residential buffers. Panel c shows adjusted pairwise and within-stratum contrasts plus the difference-in-differences (DiD) interaction. Panel d shows the DiD estimate with 95% CI, the 80% minimum detectable interaction and the prespecified ±1-year margin. The MDE is a power quantity, and the TOST results do not establish equivalence.

Statistical sources:

- `CODE/08_joint_association.R` creates the adjusted means and contrast tables.
- `CODE/45_mde_tost.R` creates the DiD precision table.
- `PYTHON/compose_main_figures.py` composes the public display from `03_source_data/Figure3a_adjusted_means.csv`, `Figure3b_contrasts.csv` and `Figure3b_DiD_precision.csv`.

The public composite is an aggregate-to-display reconstruction. It does not refit the underlying statistical models. The older R-rendered contrast output remains an analysis intermediate and is not the canonical submitted Figure 3.
