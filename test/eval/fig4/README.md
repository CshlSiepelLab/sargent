# Manuscript Figure 4 — localization (the headline)

The reframed manuscript's central new contribution: **SARGENT identifies
which regions of the genome carry trait variance, at recombination-aware
resolution**.

Three panels:

- **Panel A — Manhattan-style per-tree variance.** Estimated ŵ²_s along
  a representative simulated chromosome, with the true active blocks
  shaded. Single Block 1.2 cell (sparse + admixed + contrast 1000),
  seed showing a clean recovery.
- **Panel B — Variance-capture curves.** Cumulative ARG-variance vs
  top-k trees, for two architectures (sparse vs diffuse) and both
  truth and estimator. Demonstrates that the framework picks up
  architecture differences qualitatively, with the estimator
  systematically more concentrated than truth (a methodological caveat
  acknowledged in the caption).
- **Panel C — Localization AUC across architectures.** Estimator AUC
  (ranking active trees above null trees) as a function of contrast
  and demography. From the 30-cell Block 1.2 v1 sweep. The headline
  number: estimator above chance in 5 of 6 cells, with AUC 0.51–0.81.

## Files

- `plot.py` — generates `fig4.{pdf,png}`. Reads:
  - Panel A: `~/test-sargent/block1_2_v1/c1000_dadmixed_s1/case.*` (default
    representative cell; override via `--cell`).
  - Panel B: `~/test-sargent/phase10_admix/{sparse,diffuse}_s1/...` for
    the cross-architecture comparison.
  - Panel C: `~/test-sargent/block1_2_v1/*/est.oracle.json` for the AUC grid.

## Reproduce

```bash
# Data is already on the laptop (Phase 10 + Block 1.2 v1).  No sweep
# needed — this just builds the figure from existing outputs.
python test/eval/fig4/plot.py     # produces fig4.{pdf,png}
```

Override defaults via env / args: `--cell <dir>` for Panel A
representative cell; `--base-b12 <dir>` for Block 1.2 root;
`--base-phase10 <dir>` for Phase 10 root.
