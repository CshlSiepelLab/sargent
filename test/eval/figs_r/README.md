# Manuscript figures — R / ggplot2

These R scripts are the manuscript figure deliverables.  They render to
`~/newwriting/arg-mapping/figs2/`.  Each script reads from a small TSV
cache built by a sibling Python preprocessor (or, for Figs 2 / 3, directly
from the relevant `est.*` outputs and `*_v2*_results.csv` data files in
the parent `test/eval/`).

## Files

- `theme.R` — shared theme: Tableau-10 palette, Helvetica via cairo_pdf,
  base 8 pt, sargent_theme() + arm/flavour helpers.  Sourced by every
  fig*.R script.
- `fig2.R` — within-population h² recovery + r²_calib vs typed coverage
  (5 cov × 10 seeds × 4 arms).  Reads `est.predict.json` and
  `est.summary.json` under `~/test-sargent/fig2_n600/`.
- `fig3.R` — cross-ancestry portability, 2 panels: target-injection axis
  at F_ST = 0.16 (Panel A) and F_ST ladder at mix_B = 0.20 (Panel B).
  5 arms including SARGENT-W.  Reads `block1_1_v2_results.csv` and
  `block1_1_v2_variants_results.csv` under `test/eval/`.
- `fig3_supp.R` — companion supplementary figure: F_ST ladder at
  mix_B = 0 (no target injection).  Single panel; 5 arms.
- `fig4.R` — localization headline, 4 panels (Manhattan / variance
  capture / AUC across architectures / fused-lasso polygenic-vs-
  oligogenic).  Reads cache TSVs built by `preprocess_fig4.py`.
- `fig5.R` — architecture characterization (placeholder Panel A
  pending Block 2.2; Panels B = WHEN axis from Block 2.1 v2 data; C =
  HOW axis variance-capture curves from Block 2.5.2).  Reads cache TSVs
  built by `preprocess_fig5.py`.
- `preprocess_fig4.py` — emits `_cache/fig4_panel{A,B,C,D}.tsv` from
  case + est files (uses tskit for tree midpoints in Panel A).
- `preprocess_fig5.py` — emits `_cache/fig5_panel{B,C}.tsv` from the
  Block 2.1 v2 and Block 2.5.2 outputs.
- `fig4_overnight.sh` — driver that re-fits the 40 cells used by Fig 4
  Panels A / B / C with the recommended `--ws2-l2 0.5` default.

## Output convention

- Cache: `_cache/` (gitignored; regenerable from the preprocessors).
- Output: `~/newwriting/arg-mapping/figs2/` (PDFs + PNGs at 6.5 in wide;
  letter paper with 1 in margins).

## Conventions

- Palette: vine-style Tableau-10 (defined once in `theme.R`).
  Arms: BayesR red, GBLUP blue, SARGENT-U green, SARGENT-W teal,
  SARGENT purple.
- Error bars: light, ±1 SE, alpha 0.55, small per-arm jitter on x so
  coincident pairs (SARGENT/BayesR, SARGENT-U/SARGENT-W) are
  distinguishable.
- Theme: `theme_minimal(base_size = 8)`; axis text 7 pt; titles 9 pt
  bold; panel tags ("A", "B", ...) 11 pt bold in the top-left of each
  panel.
- Save helper `save_fig()` in `theme.R` writes PDF (cairo) + PNG (300
  dpi, white background) at the given width × height.

## Reproduce

```bash
# Each figure can be regenerated independently after its data lands
Rscript test/eval/figs_r/fig2.R
Rscript test/eval/figs_r/fig3.R
Rscript test/eval/figs_r/fig3_supp.R
mamba run -n sargent python test/eval/figs_r/preprocess_fig4.py
Rscript test/eval/figs_r/fig4.R
mamba run -n sargent python test/eval/figs_r/preprocess_fig5.py
Rscript test/eval/figs_r/fig5.R
```
