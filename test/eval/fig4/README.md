# Manuscript Figure 4 — localization (the headline)

SARGENT's central new contribution: **region-scale localization of
phenotypic variance at recombination-aware resolution**.

The manuscript deliverable is **`test/eval/figs_r/fig4.R`** (R/ggplot2),
which renders to `~/newwriting/arg-mapping/figs2/fig4.{pdf,png}`. This
`plot.py` is the original Python prototype; it implements an earlier
3-panel layout superseded on 2026-06-03 by the 4-panel R version below.

## Current R version (test/eval/figs_r/fig4.R)

Four panels (A row, BCD row) at 7.0 × 5.4 in.  Shared theme at
`test/eval/figs_r/theme.R`.

- **Panel A — Manhattan-style ŵ²_s.** Estimated per-tree variance
  weight along the chromosome on the **best-case Block 1.2 cell**
  `c1000_dpanmictic_s2` (sparse architecture, panmictic, contrast
  1000), fitted with the **recommended default `--ws2-l2 0.5`**
  (fused-lasso TV smoothing).  Truth shown as an orange step line;
  estimator as a thin blue line + dots.  Gold shading marks true
  active blocks.  Linear y-axis (the inference-noise floor falls
  below visibility).
- **Panel B — Variance-capture curves.** Cumulative-vs-top-k for
  sparse vs diffuse architectures (Phase 10 admixed cells, fitted
  with `--ws2-l2 0.5`).  Truth + estimator overlaid per architecture.
- **Panel C — Localization AUC across (architecture × demography).**
  30-cell Block 1.2 v1 grid, fitted with `--ws2-l2 0.5`.  Grouped
  bars (admixed purple vs panmictic green) at contrast ∈ {1000, 100,
  10} ("Active-region variance enrichment" axis).  AUC values:
  contrast=1000 admixed 0.859, panmictic 0.765; contrast=10 admixed
  0.659, panmictic 0.540.  Idealized ceiling overlaid as diamonds.
- **Panel D — Fused-lasso opt-in AUC contrast.** Polygenic
  (bf = 0.5, λ₂ = 0.5) vs oligogenic (bf = 0.005, λ₁ = 0.5) from 10
  Block 2.5.2 cells re-fit on the new binary.  Δ AUC = +0.036
  polygenic, +0.053 oligogenic — targeted improvement aligned with
  oligogenic architecture.

## Fused-lasso framing

The default estimator for Fig 4 production is **`--ws2-l2 0.5`** —
the fused-lasso TV smoothing prior on log ŵ²_s.  Smoke test on
`c1000_dpanmictic_s2` (the Panel A cell) confirmed L2 = 0.5 is the
optimal default for contiguous-block architectures: AUC 0.837 (no
penalty) → 0.858 (L2 = 0.5), with no improvement at L2 = 1.0.  L1
penalty was tested and rejected at this contrast (over-shrinks active
blocks).  Panel D additionally shows L1 = 0.5 as the **opt-in
oligogenic specialization** (helps bf ≤ 0.005).

## Data sources

- Panel A: `~/test-sargent/block1_2_v1/c1000_dpanmictic_s2/case.*`
  (truth) + `~/test-sargent/fig4_l2/block12_c1000_dpanmictic_s2/est.*`
  (estimator).
- Panel B: `~/test-sargent/phase10_admix/{sparse,diffuse}_s{1..5}/case.*`
  + `~/test-sargent/fig4_l2/phase10_{sparse,diffuse}_s{1..5}/est.*`.
- Panel C: `~/test-sargent/block1_2_v1/c*_d*_s*/case.*` + the
  matching `~/test-sargent/fig4_l2/block12_c*_d*_s*/est.*` cells (30).
- Panel D: `~/test-sargent/block2_5_2/bf{0.5,0.005}_dadmixed_s{1..5}/`
  default outputs + `~/test-sargent/fig4_paneld_fl/bf{0.5,0.005}_s*_l{1,2}_0.5/`
  fused-lasso re-fits.

## Driver

The `--ws2-l2 0.5` overnight re-fit of Panels A/B/C cells is at
`test/eval/figs_r/fig4_overnight.sh`.  Panel D's fused-lasso fits use
the same binary with `--ws2-l1 0.5` (oligogenic) or `--ws2-l2 0.5`
(polygenic).

## Reproduce

```bash
# After fig4_overnight.sh + Panel D L1/L2 fits have landed
mamba run -n sargent python test/eval/figs_r/preprocess_fig4.py
Rscript test/eval/figs_r/fig4.R
```

## Python legacy

`test/eval/fig4/plot.py` implements the pre-2026-06-03 3-panel layout
(Panel A on an admixed cell with the legacy estimator; no Panel D).
Kept for reference but not the manuscript deliverable.
