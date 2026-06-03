# Manuscript Figure 2 — within-population: h² and prediction vs typed coverage

Revised 2026-06-02 (mixture-axis story dropped; non-identifiability of the
SNP/ARG split moves to supplement). Both panels share the same 50 cells.

## Panels

- **Panel A — ĥ² recovery vs typed coverage.** Points + ±1 SE error bars
  for BayesR (`h²_snp` from `est.summary.json`), GBLUP (`gblup_h2`),
  SARGENT-U (`arg_blup_uniform_h2`), and SARGENT (`h²` from
  `est.summary.json`). Truth dashed at 0.60.
- **Panel B — Held-out r²_calib vs typed coverage.** Same four arms,
  same cells. Points + ±1 SE error bars.

Coverage grid: `f_typed ∈ {1.00, 0.70, 0.30, 0.10, 0.05}` (extends below
the rare_variant baseline to amplify the typed-SNP collapse). Sparse
architecture, panmictic, single population, `h²_SNP = h²_ARG = 0.3`,
`n_train = 600`, `L = 300 kb`, 10 seeds per coverage level.

## Naming

- **SARGENT** = full pipeline (BayesR fixed effect + per-tree-weighted
  ARG random effect).
- **SARGENT-U** = the **ARG-BLUP-uniform** special case (no fixed
  effect, uniform per-tree weights — pure ARG-kinship BLUP).
- **GBLUP** = typed-SNP GRM BLUP.
- **BayesR** = typed-SNP sparse fixed-effect PRS.

## Files

- `run_sweep.sh` — parallel driver (5 cov × 10 seeds = 50 cells,
  4-way parallel via `xargs -P`). Per cell: `sim → sargent → predict.py`.
- `plot.py` — generates `fig2.{pdf,png}` from
  `est.predict.json` (R²s + GBLUP/SARGENT-U h²) and `est.summary.json`
  (BayesR h², SARGENT total h²).

## Reproduce

```bash
bash test/eval/fig2/run_sweep.sh
python test/eval/fig2/plot.py --base-rva ~/test-sargent/fig2_n600
```

Override defaults via env: `N_TRAIN=800 N_SEEDS=10 N_PAR=4`
`COVERAGES="1.0 0.7 0.3 0.1 0.05"` `ARCH=sparse` `BASE=...`.
