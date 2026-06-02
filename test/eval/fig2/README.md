# Manuscript Figure 2 — variance decomposition

Reframed Fig 2 (per 2026-06-01 evening design discussion). One figure,
two panels, panmictic single-population sims throughout. Story:

> SARGENT recovers total `h²` accurately while sequentially prioritising
> the typed-SNP fixed-effect component. As true random share grows or
> typed coverage drops, the estimated SNP-vs-ARG split shifts in the
> right direction even though the partition is not strictly identifiable.

## Panels

- **Panel A — mixture sweep at full coverage.** 5 true mixtures of
  `h²_SNP` and `h²_ARG` (total = 0.6): all-SNP, heavy-SNP, balanced,
  heavy-ARG, all-ARG. 10 seeds each. Stacked bars showing SARGENT's
  estimated `h²_SNP` + `h²_ARG`, with the simulated `h²_ARG` share
  marked as a reference line. Plus overlay of BayesR (typed-SNP only)
  and ARG-BLUP-uniform (ARG only) total h².
- **Panel B — coverage sweep at balanced mixture.** Reuses the
  `~/test-sargent/rare_variant/` data (4 coverages × 2 architectures × 5
  seeds; `h²_SNP = h²_ARG = 0.3`). Lines: SARGENT `h²_SNP`,
  SARGENT `h²_ARG`, SARGENT total, BayesR total, ARG-BLUP-uniform total.
  Shaded ribbon = ±1 SE over seeds.

## Files

- `run_sweep.sh` — Panel A mixture sweep driver (50 cells; parallelism via
  `xargs -P`).
- `plot.py` — generates both panels into `fig2.{pdf,png}`. Reads:
  - Panel A: `~/test-sargent/fig2/mix*_s*/est.summary.json`
  - Panel B: `~/test-sargent/rare_variant/a*_cv*_d*_s*/est.summary.json`
- `aggregate.py` — extracts per-cell `h²` values into a combined JSON for
  plotting (called by `plot.py`).

## Reproduce

```bash
bash test/eval/fig2/run_sweep.sh        # ~2h with 6-way parallelism
python test/eval/fig2/plot.py           # produces fig2.{pdf,png}
```

Override defaults via env: `N_PAR=4 N_SEEDS=10 H2_TOTAL=0.6`
`BASE_FIG2=...` `BASE_RVA=...`.
