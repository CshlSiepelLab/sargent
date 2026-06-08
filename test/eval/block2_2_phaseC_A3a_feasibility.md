# Block 2.2 Phase C — A3a (allele-age) feasibility

**Date:** 2026-06-08. **Branch:** `phaseC-a3a-feasibility` (not pushed).
**Verdict: FEASIBLE — PASS.** A multi-component AI-REML on per-age-bin ARG
kinship recovers the per-bin h² **diagonal** that the Phase-B per-edge posterior
could not, and recovers the **bulk of total h²** (no collapse). One scaling
caveat on deep bins is documented below.

## Stage 1 — multi-component REML engine, calibrated against GCTA

`test/eval/stratified_reml.py`: AI-REML (average information) + EM warm-start
for `y = Xβ + Σ_k u_k + ε`, σ²≥0 (unconstrained reported alongside), SEs from the
AI matrix. GCTA v1.94.1 is the calibration partner (`write_gcta_grm` +
`--reml`/`--reml --mgrm` wrapper).

| test | setup | result | verdict |
|---|---|---|---|
| **T1** | single-component vs GCTA `--reml`, n=300, 5 seeds | max \|σ̂²_g(mc) − σ̂²_g(GCTA)\| = **1.5e-4**; h² matches to 3 dp | PASS (<1e-3) |
| **T2** | two-component recovery + GCTA `--mgrm`, truth (0.3, 0.2) | mean σ̂²₁=0.302±0.118, σ̂²₂=0.236±0.203; **exact GCTA agreement** every seed | PASS |
| **T3** | boundary, true σ²₂=0 | constrained σ̂²₂=0.000, unconstrained −0.071, no crash | PASS |
| **T4** | **real** ARG kinship (fig5 cells, cond(K)~9e5), h²_true=0.6 | recovers 0.616/0.407/0.734; **mc == GCTA to 1e-5** | PASS |

T2's per-component SE is ~0.11–0.14 at n=300 (the under-power the design flagged);
means land on truth, point estimates are noisy → A3a reports mean over 3 seeds.

## Stage 2 — A3a: per-age-bin ARG-REML on the Fig-5 Panel B cells

`test/eval/stratified_reml_a3a.py`. Per cell (9 cells = 3 archs × 3 seeds,
n=300, L=200 kb, h²_arg=0.85): per-cell age tertiles → per-bin diploid kinship
via the canonical assembler `kinship_assembly.build_per_bin_kinship`
(K_b = Σ_{e∈bin}(ℓ_e/2)·m_e·m_eᵀ) → standardise y, normalise each K_b to mean
diagonal 1 (GREML-MS convention) → `reml_mc` for σ²_ε + σ²_r K_r + σ²_m K_m +
σ²_a K_a. Truth-side per-bin σ² from `case.true_ws2_per_edge.tsv`
(Σ ws2_e·ℓ_e per bin), same tertiles.

**Guardrail (design risk):** Σ_b K_b vs K_ARG — **max rel err 3.9e-15** across
all 9 cells. The assembler is exact; mapping is correct.

### Per-cell σ̂²_b

| cell | σ̂²_recent | σ̂²_middle | σ̂²_ancient | σ̂²_ε | total ĥ² | conv |
|---|---|---|---|---|---|---|
| age_ancient_s1 | 0.026 | 0.000 | 2.799 | 0.131 | 0.956 | ✓ |
| age_ancient_s2 | 0.016 | 0.000 | 3.089 | 0.162 | 0.950 | ✓ |
| age_ancient_s3 | 0.005 | 0.029 | 2.169 | 0.136 | 0.942 | ✓ |
| age_middle_s1 | 0.015 | 0.954 | 0.000 | 0.083 | 0.921 | ✓ |
| age_middle_s2 | 0.050 | 0.786 | 0.007 | 0.099 | 0.895 | ✓ |
| age_middle_s3 | 0.028 | 0.816 | 0.080 | 0.105 | 0.898 | ✓ |
| age_recent_s1 | 0.909 | 0.028 | 0.018 | 0.071 | 0.931 | ✓ |
| age_recent_s2 | 0.788 | 0.013 | 0.000 | 0.145 | 0.847 | ✓ |
| age_recent_s3 | 0.953 | 0.100 | 0.000 | 0.049 | 0.956 | ✓ |

### Diagonal matrix (mean over 3 seeds)

| true arch | recent | middle | ancient | diag−2nd | PASS |
|---|---|---|---|---|---|
| age_recent | **0.884** | 0.047 | 0.006 | +0.836 | ✓ |
| age_middle | 0.031 | **0.852** | 0.029 | +0.821 | ✓ |
| age_ancient | 0.016 | 0.010 | **2.686** | +2.670 | ✓ |

## Pass/fail

- **Diagonal recovery — PASS.** Every architecture's matching bin dominates by
  ≥ 0.10 (actually ≥ 0.82). Clean recovery of the per-edge architecture that the
  Phase-B per-edge posterior could not localize.
- **Total ĥ² ≥ 0.5 — PASS** (mean 0.922; simulated 0.85). The Phase-B
  collapsed-h² failure mode is resolved.
- **Sum-rule — PASS** (3.9e-15).
- **Engine vs GCTA — PASS** (T1–T4).

## Surprises / caveats (honest report)

1. **Deep-bin σ̂² scale inflation.** The ancient architecture's σ̂²_ancient reads
   ~2.7 (> 1), which is not a literal h². Cause: ancient edges are deep, their
   descendant-indicator vectors cover most individuals, so **K_ancient is
   near-rank-1**; mean-diagonal normalization then does not map σ̂²_b to h² units
   for deep/rank-concentrated bins. The **diagonal pattern is robust** (ancient
   clearly dominates), but absolute per-bin h² for deep bins needs a better
   normalization (e.g. variance-explained `σ̂²_b·tr(K_b)/n / Vp`, or report
   per-bin fractions). This is a reporting refinement, not a method failure —
   the σ²_b matrix is firmly diagonal and total h² did not collapse.
2. **n=300 two-component SE** ~0.1 (Stage-1 T2). A3a's three-seed means are
   clean because the architectures are near-degenerate (one bin carries almost
   all variance), so the signal ≫ SE. Mixed architectures would need more seeds
   or larger n.

## Next steps (require Adam's go — do not start)

Acts I–II (calibration vs SNP-stratified REML on functional annotations), then
A3b (frequency) / A3c (population). Per the brief, A3a feasibility reports here
first; nothing downstream starts without authorization.

## Reproduce

```bash
PYBIN=/DATA/home/asiepel/miniforge3/envs/sargent/bin/python
# Stage 1 calibration (needs gcta64 on the env PATH):
$PYBIN test/eval/stratified_reml.py
# Stage 2 A3a (needs the 9 fig5_panelB_edges cells):
$PYBIN test/eval/stratified_reml_a3a.py
# -> test/eval/figs_r/_cache/fig5_panelB_A3a.tsv
```
