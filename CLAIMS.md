# SARGENT — Claims table

Working triage of claims the manuscript might make, against current evidence.
Each claim is bucketed:

- **HELD** — supported by existing data; ready for the manuscript.
- **PLAUSIBLE / UNTESTED** — mechanistically motivated; needs an experiment.
- **FALSIFIED / WEAK** — undermined by current evidence; drop or substantially reframe.
- **OUT OF SCOPE** — defensible but not for this paper.

A claim's bucket determines what Block-1 experiment, if any, it generates.
*Draft for Adam's review (2026-05-29); please edit buckets, add/remove claims,
and flag anything I've miscategorized.*

A second annotation, **[MINOR]**, marks claims that are not central to the
story — either not currently a strong claim in the manuscript, or easily
removed without weakening the narrative. These do not need to drive
experiments and can be reframed or dropped quietly during Block 3 rework.
Per Adam (2026-05-29): the obviation of the finite-mixture `wₛ²` prior is
already well established and does not need to be dwelt on.

Tally (draft): HELD 13 · PLAUSIBLE/UNTESTED 12 · FALSIFIED/WEAK 11 · OUT OF SCOPE 2.
Of the FALSIFIED/WEAK group, 5 are flagged **[MINOR]**: A2, B3, C2, C4, F2.

---

## A. Variance components (total `h²`, decomposition)

| # | Claim | Bucket | Evidence held | Evidence missing | Source |
|---|---|---|---|---|---|
| A1 | SARGENT recovers total `h²` accurately. | **HELD** | Phase 6 n=1000: h²_total est 0.81–0.85 vs truth 0.76–0.90 (within 3–15%). SARGENT `--arg-only` on Phase 10: bias +0.06–0.10. | Direct comparison vs GCTA / eGRM / tslmm on identical input. | ROADMAP Phase 6; `phase10_oracle.md`. |
| A2 | SARGENT separates `h²` into SNP (fixed) vs ARG (random) components. | **FALSIFIED · [MINOR]** | Phase 1 concession: not identifiable in principle (X is an observed subset of the genealogy that K_ARG encodes in expectation). | — (it is structurally non-identifiable; no experiment will resolve it). Not a strong manuscript claim; reframe Stage-1/Stage-2 as an "operational prioritization" rather than a decomposition. | ROADMAP Phase 1. |
| A3 | Per-segment `wₛ²` heterogeneity improves `h²` recovery vs uniform-weight ARG-LMM. | **FALSIFIED** | Phase 10: uniform ≈ oracle ≈ estimated `h²` to ~2nd decimal across both diffuse and sparse architectures. | — (the BLUP absorbs kinship into σ_u², so this is a structural negative). | `phase10_oracle.md` (h² addendum). |
| A4 | SARGENT recovers per-region (partitioned) `h²` respecting recombination boundaries. | **PLAUSIBLE / UNTESTED** | Phase 10 localization ceiling 0.84 + estimator 0.74 in one architecture imply the per-region track has signal. | Quantitative per-region `h²` recovery vs truth; comparison vs S-LDSC / MESC. | Block 1.2 sweep + new "partitioned-h² recovery" metric. |

## B. Prediction

| # | Claim | Bucket | Evidence held | Evidence missing | Source |
|---|---|---|---|---|---|
| B1 | SARGENT (sequential) beats BayesR at low typing coverage. | **HELD** | Phase 3 f=0.03: SARGENT/ARG-BLUP-uni R²=0.58 vs BayesR 0.07 (5 seeds). | — (held empirically, but caveat: the gain is the ARG kinship's, not per-segment weights'). | ROADMAP Phase 3 overnight; `summary.md`. |
| B2 | SARGENT beats ARG-BLUP-uniform within-population. | **FALSIFIED** | Phase 10 sparse + admixed: oracle ceiling +0.023 R² over uniform (paired SE 0.006); estimator at or below uniform (−0.008). Phase 3 within-population: SARGENT ≈ ARG-BLUP-uniform. | — (the ceiling itself is small; no experiment will lift it materially). | `phase10_oracle.md`. |
| B3 | Joint-EM polish (`--em N`) gains 3–7% R² over sequential within-population. | **WEAK · [MINOR]** | Phase 3 numerical gain 3–7%; but seed SDs ±0.15–0.27 dwarf the effect. | Not a load-bearing claim; drop or fold into a brief "Stage-2 polish converges further" Methods note. | ROADMAP Phase 3. |
| B4 | SARGENT improves cross-ancestry PRS portability. | **PLAUSIBLE / UNTESTED** | Phase 3 low-coverage win is a *proxy* for portability (typed SNPs miss causal); mechanism for ARG kinship to help is sound. | Cross-population transfer experiment (train pop A, predict pop B). | Block 1.1 (the prediction headline). |
| B5 | SARGENT extends polygenic prediction to deeply diverged populations (Neandertal-like). | **PLAUSIBLE / UNTESTED** | Mechanism: typed shared variants port their effects; ARG random component shrinks under Brownian-motion prior. | Simulation-based proof of concept at archaic divergence depths. | Beyond Block 1.1; tentative future. |
| B6 | SARGENT's ARG kinship recovers prediction value from untyped (SV-like) variants. | **MIXED / REFRAME** | LD-leaky regime (Phase 3 low-coverage): picks up missing variance via tag-LD. Clean rare-MAF regime (Phase 6b/c Class S sweep): no recovery (n=200 underpowered for anyone). | Adequately-powered (large-n) test of rare-variant recovery in the clean regime. | `class_ps_sweep.txt`; ROADMAP Phase 6b. |
| B7 | SARGENT's hybrid (typed-SNP fixed effects + ARG random) is the right object for transfer. | **PLAUSIBLE / UNTESTED** | Mechanistic: typed shared variants port their fixed effects directly; ARG carries ancestry-varying remainder. No within-pop result distinguishes the hybrid from ARG-BLUP-uniform. | Empirical separation of SARGENT from ARG-BLUP-uniform in a transfer setting. | Block 1.1 (the linchpin). |

## C. Localization / partitioned-`h²`

| # | Claim | Bucket | Evidence held | Evidence missing | Source |
|---|---|---|---|---|---|
| C1 | SARGENT identifies active genealogical regions carrying missing heritability, at *region* scale. | **HELD (one architecture)** | Phase 10 sparse + admixed: estimator AUC 0.74 (range 0.59–0.90), ceiling 0.84. | Robustness across architectures (sparsity, concentration, panmictic vs admixed). | `phase10_oracle.md`. |
| C2 | SARGENT localizes at *single-tree* resolution. | **FALSIFIED · [MINOR]** | Earlier per-tree localization debugging: AUC ≈ 0.5 at n=100–300; adjacent-tree correlation + information-starvation cap fine-scale resolution. | — (structural barrier). Not a strong manuscript claim; reframe localization as *region-scale* (C1/C3/A4) without ever promising single-tree. | `project_next_steps_predict_localize.md`. |
| C3 | SARGENT's local-`h²` track is sharper than fixed-window partitioning (S-LDSC / MESC) because it respects recombination boundaries. | **PLAUSIBLE / UNTESTED** | Conceptually compelling: ARG segments are recombination-defined; fixed windows are not. | Direct comparison vs S-LDSC / MESC / partitioned-LDSC on the same data. | Block 1.2 stretch goal. |
| C4 | The current per-tree REML estimator localizes reliably across signal strengths. | **WEAK / UNTESTED · [MINOR]** | One architecture tested. Within that architecture, est AUC has large seed variance (0.59–0.90). | Multi-architecture sweep (Block 1.2). Estimator quality is a Block-1.2/Block-2 *internal* benchmark, not a manuscript claim per se — the paper-facing claim is C1 (region-scale localization works), not "this specific estimator is the right one." | Block 1.2. |
| C5 | A sparse-fused-lasso prior on log-`wₛ²` closes the localization gap (estimator → ceiling AUC). | **PLAUSIBLE / UNTESTED** | Mechanism: regularization concentrates the limited information budget; the prior matches the favorable architecture. | Build and validate against the Phase 10 ceiling. | Block 2 (Phase 11; conditional). |

## D. Genetic architecture / polygenicity

| # | Claim | Bucket | Evidence held | Evidence missing | Source |
|---|---|---|---|---|---|
| D1 | SARGENT estimates polygenicity (π_k from BayesR mixture) accurately. | **WEAK** | π_k is biased by the same fixed/random cannibalization that biases `h²_snp` / `h²_arg`. | Targeted test on simulations with known polygenicity. | Phase 1 / `project_stage2_algorithm.md`. |
| D2 | SARGENT distinguishes polygenicity regimes (sparse vs diffuse). | **PLAUSIBLE / UNTESTED** | The mixture prior is structurally capable of this. | A 2–3 regime comparison. | Old Phase 3 design (untested at depth); low priority. |

## E. Structured / admixed populations

| # | Claim | Bucket | Evidence held | Evidence missing | Source |
|---|---|---|---|---|---|
| E1 | SARGENT is particularly well suited for structured / admixed populations. | **PLAUSIBLE / UNTESTED** | Phase 10 used admixed sims and confirms the mechanism, but only as a setup for the oracle diagnostic — no head-to-head methods comparison under structure. | Comparison of methods specifically in structured/admixed cohorts. | Block 1.1 (the same experiment that does B4). |
| E2 | SARGENT prediction works in admixed individuals (African Americans, Latin Americans) via local-ancestry-aware ARG kinship. | **PLAUSIBLE / UNTESTED** | Admixture pulse model already in `simulate.py`; mechanism aligned with funded scope (All of Us / HARC). | Real-data extension or admixed-sample simulation with explicit local-ancestry annotation. | Block 1.1 extension. |

## F. Methodological / mathematical

| # | Claim | Bucket | Evidence held | Evidence missing | Source |
|---|---|---|---|---|---|
| F1 | Mean-field VB with the collapsed (E[V]) bound gives sensible posteriors. | **HELD** | Standard derivation; algorithm works as intended in tests. | — | `main.tex` Methods; production code. |
| F2 | The per-segment mixture-prior update for `wₛ²` (as derived in `main.tex`) drives accurate weight recovery. | **FALSIFIED in code · [MINOR]** | The derivation in `main.tex` is correct, but the *implementation* abandons it (signal too small in `S ≫ n`) for point-estimate REML + post-hoc ρ. Obviation of the finite-mixture `wₛ²` prior is already well established and does not need to be dwelt on. | — (manuscript ↔ code mismatch). Methods fix: replace the mixture derivation with the point-estimate REML + post-hoc ρ that the code actually runs; no need to motivate at length. | `main.tex` eq. for `rho_{sm}`; `variational.c:575`. |
| F3 | The sequential inference strategy mitigates fixed/random competition. | **PARTIAL / REFRAME** | Empirically: Stage 1 does prioritize SNP fixed effects (h²_snp inflates, h²_arg shrinks in sequential mode). | The "decomposition" is not actually identifiable (A2). Should be reframed as "an operational choice for prioritization," not "recovers the true split." | ROADMAP Phase 1. |
| F4 | SARGENT reduces to BayesR, GREML, ARG-LMM, phylogenetic LMM in special cases. | **HELD (cosmetic)** | Structural property of the model. | — | `main.tex` Results overview. |
| F5 | The post-hoc logistic ρ map gives interpretable per-segment activation probabilities. | **HELD** | Implemented and stable; matches the prior-mean anchor. | Calibration check (does ρ ≈ posterior probability of activation?). | `variational.c:621–653`. |

## G. Scalability / computational

| # | Claim | Bucket | Evidence held | Evidence missing | Source |
|---|---|---|---|---|---|
| G1 | SARGENT inference is linear-time / scalable to biobank cohorts. | **FALSIFIED (current implementation)** | Phase 12: per-iteration cost is `O(n³)` (CG iters × matvec × S, all ~linear in n). Phase 6 n=1000: ~10h wall. | Preconditioner + warm-start CG (Block 4) recovers ~O(n²); fast tree-sequence matvec (Lee/Zhu) goes lower. Biobank-scale test never run. | Phase 12 (now Block 4). |
| G2 | Hutchinson trace estimation enables efficient variance-component fitting. | **HELD** | Works as designed; default 10 probes adequate at tested scales. | — | `hutchinson.c`. |
| G3 | CG-based V-solve avoids dense factorization. | **HELD** | No `dpotrf` / `dgesv` in production path. | — (but iterations grow with n without a preconditioner — see G1). | `conjgrad.c`. |
| G4 | Robust to ARG inference error (close relatives dominate). | **PLAUSIBLE / UNTESTED** | Manuscript heuristic. Phase 3+ ran on simplified true ARGs, not inferred. | A `tsinfer`/`tsdate`-inferred-ARG sensitivity test. | `main.tex` "Effect of ARG inference" (placeholder). |

## H. Comparisons / unification

| # | Claim | Bucket | Evidence held | Evidence missing | Source |
|---|---|---|---|---|---|
| H1 | SARGENT improves over BayesR (SNP-only sparse regression). | **HELD (low coverage)** | Phase 3 f=0.03. | Stronger comparison under realistic missing-variant scenarios. | Phase 3. |
| H2 | SARGENT improves over GBLUP (standard genomic GRM). | **HELD (low coverage)** | Phase 3 f=0.03. | — | Phase 3. |
| H3 | SARGENT improves over ARG-BLUP-uniform / eGRM / tslmm within-population. | **FALSIFIED** | Phase 10 + Phase 3 — no within-population gap. | — | Phase 10. |
| H4 | SARGENT improves over multi-kernel LMM (mkLMM) / windowed local GRM / hapla. | **UNTESTED** | None implemented as baselines. | Implement minimal versions of the 2–3 most relevant baselines, run side-by-side. | `main.tex` comparison table (lists 9; need triage). |
| H5 | SARGENT (C) is faster/more accurate than Nurdan's R implementation. | **UNTESTED** | Phase 3.5 smoke ran; no side-by-side. | Block 1 or parking-lot; gated on availability of Nurdan's code. | Old Phase 6d. |
| H6 | SARGENT is the right object for transfer / portability *because* it integrates fixed + random (vs eGRM which uses ARG kinship alone). | **PLAUSIBLE / UNTESTED** | Mechanism is sound; no empirical demonstration. | The Block 1.1 transfer experiment, specifically including a per-arm decomposition (BayesR / ARG-BLUP-uniform / SARGENT). | Block 1.1. |

## I. Out of scope (for this paper)

- I1. Real-data analysis on a major biobank (UK Biobank, AoU). *(Future work.)*
- I2. Fine-grained polygenic-score deployment for clinical reporting. *(Future work.)*

---

## What this drives in Block 1

Re-grouping the **PLAUSIBLE / UNTESTED** claims by experiment:

- **Block 1.1 — PRS portability / cross-ancestry transfer**: covers B4, B7, E1, E2, H6, partly B5.
- **Block 1.2 — Localization robustness sweep**: covers A4, C1 (breadth), C4, partly C3.
- **Block 1.3 — `predict.py:fit_blup` audit**: prerequisite for any honest `h²` claim (A1 absolute number).
- **Possibly Block 2 — Phase 11 (sparse fused lasso)**: gated on Block 1.2 outcome → C5, and likely needed to fully establish C3 / A4.
- **Deferred (parking lot)**: B5 (Neandertal), G4 (ARG-inference sensitivity), H4 (extra-baselines triage), H5 (Nurdan comparison).

## Manuscript-side consequences (for next week's Block 3)

**Load-bearing drops/reframes** (these *matter* for the new story):
- **A3** (per-segment `wₛ²` improves `h²`) — drop; reframe `wₛ²` as a *localization* lever, not a `h²`-recovery lever.
- **B2** (SARGENT > ARG-BLUP-uniform within-pop) — drop; reposition the prediction headline around portability (B4/B7/E1/E2/H6).
- **G1** (linear-time scalability) — drop the claim; flag computational efficiency as future work (Block 4).

**Light reframes** ([MINOR] — these were never strong manuscript claims; trim without dwelling):
- **A2** (fixed/random split) — replace with "operational prioritization" framing in a single sentence.
- **B3** (`--em` polish gain) — drop the headline number; keep as one-line Methods note.
- **C2** (single-tree resolution) — never promise tree-level; commit to *region-scale* throughout (C1/C3/A4).
- **C4** (per-tree REML estimator quality) — internal to Block-1.2/Block-2; do not surface as a paper claim.
- **F2** (finite-mixture `wₛ²` prior) — replace the derivation in Methods with the point-estimate REML the code actually runs; no need to motivate the switch at length (obviation is well established).

**Reframe as primary contributions**: C1/C3/A4 (region-scale localization, partitioned `h²`); B4/E1/E2/H6 (portability via ARG kinship + sparse fixed effects).

**Comparison table**: triage H4 from 9 baselines to 2–3 essential ones.
