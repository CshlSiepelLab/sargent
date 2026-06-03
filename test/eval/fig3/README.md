# Manuscript Figure 3 — cross-ancestry portability

Per 2026-06-02 plan (revised Option 3 in the Fig 2 + Fig 3 reconsideration).

Two panels:

- **Panel A — R² vs typed coverage at F_ST = 0.16, mix_B = 0.**
  Grouped bars over typed = {0.10, 0.70} for the four arms (SARGENT,
  SARGENT-U, GBLUP, BayesR). SARGENT-U = the SARGENT-U special
  case of SARGENT (no fixed effect, uniform per-tree weights — pure
  ARG-kinship BLUP). Shows SARGENT-U is coverage-invariant while BayesR
  / GBLUP / SARGENT all collapse at low coverage. This is the
  load-bearing comparison: the ARG-vs-typed gap widens ~3-4× as typed
  coverage drops (Block 1.1 v2 verdict).
- **Panel B — R² vs F_ST at typed = 0.70, mix_B = 0.**
  Line plot over F_ST ∈ {0.07, 0.12, 0.16, 0.30} with ±1 SE ribbons.
  Shows the ARG framework retains its advantage across the entire
  ladder; BayesR and SARGENT (which inherits the BayesR fixed effect)
  collapse hardest. Frame: advantage *preserved* across F_ST, not
  *widening* with F_ST (the original headline did not survive v2).

Note on the framework representative: **SARGENT-U**, not
SARGENT, is the figure's framework representative. SARGENT inherits
the BayesR fixed-effect liability cross-ancestry; the within-method
SARGENT-vs-SARGENT-U null moves to the supplement.

## Files

- `plot.py` — generates `fig3.{pdf,png}` from
  `test/eval/block1_1_v2_results.csv` (committed in `ad348d5`).
- Data: `test/eval/block1_1_v2_results.csv` (120 rows, schema
  `split, fst, mix, typed, seed, <arm>_r2, <arm>_r2c, h2_snp, h2_arg`
  for arms ∈ {bayesR, gblup, argunif, sargent}).

## Reproduce

```bash
python test/eval/fig3/plot.py     # produces fig3.{pdf,png}
```

Override defaults via flags:

```bash
python test/eval/fig3/plot.py \
    --panel-a-fst 0.30 --panel-b-typed 0.10 \
    --out /tmp/fig3_alt
```
