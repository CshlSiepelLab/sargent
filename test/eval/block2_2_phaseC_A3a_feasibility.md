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
(K_b = Σ_{e∈bin}(ℓ_e/2)·m_e·m_eᵀ) → standardise y, **center each K_b (H K H,
H = I − 11ᵀ/n)** then normalise to mean diagonal 1 → `reml_mc` for σ²_ε +
σ²_r K_r + σ²_m K_m + σ²_a K_a → report h²_b = σ̂²_b/Vp. Truth-side per-bin σ²
from `case.true_ws2_per_edge.tsv` (Σ ws2_e·ℓ_e per bin), same tertiles.

**Why centering (the deep-bin fix):** ancient edges are deep, so their
descendant indicators align with the all-ones vector, which REML projects out
via the intercept. Without centering, that mass inflated σ̂²_ancient (~2.7,
Vp ≫ Var(y)) — an off-scale but still diagonal estimate. Centering each K_b
removes the intercept-aligned component before the fit; σ̂²_ancient drops to
0.835 and all bins read as clean h² ∈ [0,1]. (This also matches the centered
SNP-GRMs used in Acts I–II, so ARG-REML and SNP-REML stay comparable.)

**Guardrail (design risk):** Σ_b K_b vs K_ARG — **max rel err 3.9e-15** across
all 9 cells. The assembler is exact; mapping is correct.

### Per-cell ĥ²_b (centered)

| cell | ĥ²_recent | ĥ²_middle | ĥ²_ancient | σ̂²_ε | total ĥ² | conv |
|---|---|---|---|---|---|---|
| age_ancient_s1 | 0.022 | 0.000 | 0.866 | 0.131 | 0.888 | ✓ |
| age_ancient_s2 | 0.014 | 0.000 | 0.845 | 0.162 | 0.860 | ✓ |
| age_ancient_s3 | 0.005 | 0.035 | 0.795 | 0.136 | 0.835 | ✓ |
| age_middle_s1 | 0.014 | 0.906 | 0.000 | 0.083 | 0.921 | ✓ |
| age_middle_s2 | 0.053 | 0.838 | 0.002 | 0.099 | 0.894 | ✓ |
| age_middle_s3 | 0.029 | 0.837 | 0.025 | 0.105 | 0.891 | ✓ |
| age_recent_s1 | 0.896 | 0.028 | 0.006 | 0.071 | 0.930 | ✓ |
| age_recent_s2 | 0.832 | 0.014 | 0.000 | 0.145 | 0.846 | ✓ |
| age_recent_s3 | 0.865 | 0.091 | 0.000 | 0.049 | 0.956 | ✓ |

### Diagonal matrix (mean ĥ²_b over 3 seeds)

| true arch | recent | middle | ancient | diag−2nd | PASS |
|---|---|---|---|---|---|
| age_recent | **0.864** | 0.044 | 0.002 | +0.820 | ✓ |
| age_middle | 0.032 | **0.861** | 0.009 | +0.828 | ✓ |
| age_ancient | 0.014 | 0.012 | **0.835** | +0.821 | ✓ |

Every architecture recovers ĥ² ≈ 0.85 in its matching bin (≈ the simulated
0.85), off-diagonals ≤ 0.04. Clean, symmetric diagonal.

## Pass/fail

- **Diagonal recovery — PASS.** Every architecture's matching bin dominates by
  ≥ 0.10 (actually ≥ 0.82). Clean recovery of the per-edge architecture that the
  Phase-B per-edge posterior could not localize.
- **Total ĥ² ≥ 0.5 — PASS** (mean 0.891; simulated 0.85). The Phase-B
  collapsed-h² failure mode is resolved.
- **Sum-rule — PASS** (3.9e-15).
- **Engine vs GCTA — PASS** (T1–T4).

## Surprises / caveats (honest report)

1. **Deep-bin σ̂² scale inflation — found and FIXED.** Before centering, the
   ancient architecture's σ̂²_ancient read ~2.7 (> 1, not a literal h²): ancient
   edges are deep, their descendant indicators align with the all-ones vector,
   and mean-diagonal normalization did not map σ̂²_b to h² units for such
   rank-concentrated bins. **Fix: center each K_b (H K H) before normalizing**
   (standard GREML GRM processing) — this removes the intercept-aligned
   component the REML projects out anyway, and σ̂²_ancient drops to 0.835, all
   bins clean h² ∈ [0,1]. The diagonal was always robust; centering makes the
   absolute scale correct too. (Same centering will apply to the SNP-GRMs in
   Acts I–II, keeping ARG-REML and SNP-REML comparable.)
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
