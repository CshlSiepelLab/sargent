# Manuscript Figure 3 — cross-ancestry portability

The manuscript deliverable is **`test/eval/figs_r/fig3.R`** (R/ggplot2),
which renders to `~/newwriting/arg-mapping/figs2/fig3.{pdf,png}`. This
`plot.py` is the original Python prototype, preserved for reference; it
implements an earlier panel design (typed-coverage axis at mix_B=0) that
was superseded on 2026-06-02 (pm) by the target-injection layout below.

## Current R version (test/eval/figs_r/fig3.R)

Two panels at 6.5 × 3.3 in, using the shared theme at
`test/eval/figs_r/theme.R` (Tableau-10 palette, Helvetica via
`cairo_pdf`, base 8 pt).

- **Panel A — Target-injection axis at F_ST = 0.16, typed = 0.70.**
  Grouped bars over mix_B ∈ {0.00, 0.05, 0.20} for 5 arms (BayesR,
  GBLUP, SARGENT, SARGENT-U, SARGENT-W) with light error bars.  Story:
  without target injection cross-ancestry prediction is dismal; with
  mix_B = 0.20 SARGENT-U lifts to r²_calib ≈ 0.22 while SARGENT-with-β
  stays at ≈ 0.07 — confirming β is the cross-ancestry liability.
  F_ST = 0.16 is annotated in-panel.
- **Panel B — F_ST ladder at typed = 0.70, mix_B = 0.20 (the predictive
  regime).**  Lines for the same 5 arms with light per-arm error bars
  (small interleaved jitter so coincident pairs are distinguishable).
  Story: SARGENT-U dominates across F_ST ∈ {0.07, 0.12, 0.16, 0.30}.
  SARGENT-W tracks SARGENT-U closely (rescue is essentially complete
  with target injection).  SARGENT-with-β collapses with BayesR at
  F_ST ≥ 0.16.

## Companion supplementary figure

`test/eval/figs_r/fig3_supp.R` renders the same F_ST ladder but at
**mix_B = 0** (no target injection).  Same 5 arms.  Shows the absolute
collapse without target individuals and the partial rescue verdict
(SARGENT-W < SARGENT-U by 0.01–0.03).  Output at
`~/newwriting/arg-mapping/figs2/fig3_supp.{pdf,png}`.

## Naming

SARGENT-U = pure-ARG-kinship BLUP, no β, uniform per-tree weights
(equivalent to ARG-BLUP-uniform).  SARGENT-W = `sargent --arg-only` =
per-tree-weighted ARG without β.  See `test/eval/block1_1_design.md`
§ "SARGENT cross-ancestry rescue test" for the full verdict.

## Data sources

- `test/eval/block1_1_v2_results.csv` (server-side Block 1.1 v2,
  commit `ad348d5`; 4 F_ST × 3 mix × 2 typed × 5 seeds = 120 cells).
- `test/eval/block1_1_v2_variants_results.csv` (SARGENT-W rescue
  test, commit `b7f9d73`; same 120 cells; schema adds
  `sargentW_r2c`).

## Reproduce

```bash
Rscript test/eval/figs_r/fig3.R          # 2-panel main figure
Rscript test/eval/figs_r/fig3_supp.R     # mix_B = 0 supplementary
```

## Python legacy

`test/eval/fig3/plot.py` implements the pre-2026-06-02 layout
(typed-coverage axis at mix_B = 0).  Kept for reference but not the
manuscript deliverable.
