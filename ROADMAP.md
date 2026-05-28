# SARGENT — Development Roadmap

Last updated: 2026-05-29.

This document captures the agreed plan for finalising SARGENT. The
priority order is **correctness first, scale later**. We finish each
algorithmic mode at small scale, validate it, then tackle the engineering
optimisations needed for biobank scale.

> **2026-05-27 strategic re-evaluation.** A full review of code, results,
> and manuscript concluded that the distinctive per-segment `ws²` machinery
> has not yet beaten ARG-BLUP-uniform on *any* experiment, and that the
> fixed/random decomposition and fine-scale localization claims are not
> supported. The two concrete follow-ups agreed are **Phase 10 (oracle
> diagnostic)** and **Phase 11 (hierarchical-centering + sparse-fused-lasso
> `ws²` redesign)** — see the new sections at the end of this file. These
> now take priority over the open items in Phases 6b/6c/7.
>
> **2026-05-28 — Phase 10 verdict.** Under genuinely optimistic conditions
> (admixed demography F_ST≈0.26, sparse 5%-active concentrated blocks,
> contrast 1000, 5 seeds at n=300+100), the oracle diagnostic gives a sharp
> answer: prediction ceiling oracle−uniform = **+0.023 ± 0.013 R²** (small
> but reliably positive; estimator gap −0.008, i.e. current per-tree REML
> captures *none* of it), while the localization ceiling = **0.844 AUC**
> with the current estimator at **0.74** (clear ~0.10 AUC headroom). The
> reframe this forces: **per-segment `ws²` is a localization tool, not a
> prediction tool.** Phase 11 (sparse fused lasso) is worth building, but
> as a *localization / local-`h²`* method, with prediction held at
> uniform-equivalent as a safety floor. The manuscript's prediction story
> belongs to the ARG kinship + structure-aware demography (transfer /
> admixture / Neandertal), not to weight heterogeneity. See
> `test/eval/phase10_oracle.md` and the updated Phase 11 below.
>
> **2026-05-29 — reorganisation.** Replaced the old chronological phase
> list with a goal-driven **Active plan: Blocks 0–5** below. Priority is
> now (1) home in on a defensible "what the method provides vs. does
> not" — see `CLAIMS.md`; (2) run the minimum experiments to back that
> story; (3) rework the manuscript (next week, Adam-led). Lower-priority
> items from the old phases have been moved to **Parking lot — deferred
> items** at the end of the active plan, not deleted. The completed and
> retired phases remain in place below for historical reference.

If this session is interrupted, restart from the next un-checked item
in **Block 0 / Block 1** in the Active plan immediately below.

---

## Active plan (2026-05-29 onward): Blocks 0–5

### Block 0 — Claims table (immediate, lightweight)

**Goal**: a triage of every claim the manuscript might make, bucketed into
*held / plausible-untested / falsified-weak / out-of-scope*, with the
evidence we hold and the gap each one creates. The exercise sets the
experimental agenda for Block 1 by exposing exactly which untested claims
we want to support.

- [x] Draft `CLAIMS.md` (claims A1–H6 + I; tallies: HELD 13, plausible 12,
      falsified 11). 2026-05-29.
- [ ] Adam reviews and edits buckets / claims; we converge on the final
      set of "plausible / untested" claims that drive Block 1.

### Block 1 — Critical experiments to back the tight story

#### Block 1.1 — PRS portability / cross-ancestry transfer

**Goal**: the prediction headline. Train in pop A, predict in pop B (and
mixed-train, predict in A and B separately) across a range of split times
/ genetic distances. Backs claims B4, B7, E1, E2, H6.

- [ ] Extend `scripts/simulate.py` to allow sampling *separately* from
      the two source populations (`--demography admixed` currently
      samples only from ADM). Add `--sample-from {ADM,A,B,mixed}` or
      equivalent.
- [ ] Driver `test/eval/portability_sweep.sh`: simulate (sweep split time
      ∈ {1k, 2k, 4k, 8k} gen) × (train: pop A; test: pop A, B, mixed)
      × seeds.
- [ ] Score BayesR, GBLUP, ARG-BLUP-uniform, sequential SARGENT, joint-EM
      SARGENT. Per-arm OOS R² with paired CIs across seeds.
- [ ] Write-up `test/eval/portability.md` with the verdict on whether the
      hybrid distinguishes itself in transfer (the linchpin claim).

#### Block 1.2 — Localization robustness sweep

**Goal**: extend Phase 10's single architecture to a grid so the
localization claim is robust, not anecdotal. Backs claims A4, C1, C4
(and gates Block 2).

- [ ] Sweep grid: demography ∈ {panmictic, admixed (split 6k)} × sparsity
      ∈ {1%, 5%, 20% of trees active} × contrast ∈ {10, 100, 1000} × ≥5
      seeds. (Reuse `oracle_ws2.{py,sh}`.)
- [ ] Report estimator AUC vs Phase 10's localization ceiling across the
      grid (one heatmap).
- [ ] Add per-region `h²` recovery metric: regress the truth-per-block
      `h²` on the estimator's local-`h²` track and report `R²` + slope.
      (Backs A4.)
- [ ] Write-up `test/eval/localization_sweep.md` with the verdict on
      whether the current per-tree REML is good enough across
      architectures, or whether Phase 11 is needed.

#### Block 1.3 — `predict.py:fit_blup` h² bias audit

**Goal**: explain and ideally fix the +0.25 h² over-estimation surfaced by
Phase 10. Prerequisite for any absolute h² number in the manuscript. Backs
honest reporting of A1, and (re-)validates Phase 3 / Phase 6 / Phase 6b
tables. Details already drafted under **Cross-cutting concerns →
`predict.py fit_blup` upward h² bias audit**.

- [ ] Reproduce the bias on a clean GRM from random genotypes (no admixture).
- [ ] Implement and test the proper REML projection
      `P = V⁻¹ − V⁻¹X(X'V⁻¹X)⁻¹X'V⁻¹`.
- [ ] Sweep `h2_init` / `max_iter`; cross-check against `gcta --reml`.
- [ ] Annotate or correct historical eval write-ups once root cause is in
      hand.

### Block 2 — Phase 11 (sparse fused lasso) — **gated on Block 1.2**

**Goal**: only build if Block 1.2 shows the current estimator's
localization is inadequate under realistic / weaker / diffuse architectures.
If the current method is robust enough across the grid, defer Phase 11.

Detailed design lives at `## Phase 11 — ws² redesign…` below; if Block 1.2
triggers it, that section becomes the implementation plan. Success
criterion: closes ≥ 50 % of the estimator→ceiling AUC gap, no prediction
regression, diffuse-control selects zero blocks.

### Block 3 — Manuscript rework (next week, Adam-led)

**Goal**: triggered once Blocks 0 + 1 are far enough along to know what's
defensible. Major restructure, not incremental edit.

- [ ] Skeleton + figure plan + drop-list / rebuild-list (from `CLAIMS.md`).
- [ ] **Methods**: fix the manuscript ↔ code mismatch on the `wₛ²` update
      (paper still derives a mixture update; code uses point-estimate
      REML + post-hoc ρ).  Document the sequential strategy as an
      operational choice for prioritisation, not a true decomposition.
- [ ] Comparison-table triage: of the 9 baselines currently listed,
      decide which 2–3 are essential to implement, drop the rest.
- [ ] New figures: region-scale local-`h²` Manhattan; cross-ancestry
      portability; possibly a reframed conceptual figure that emphasizes
      the localization story.

### Block 4 — Computational efficiency (deferred until story is locked)

**Goal**: only needed if the manuscript makes biobank-scale claims; if the
contribution is framed methodologically with scale-up as follow-up work,
this stays deferred. Was Phase 12.

- [ ] Verify (instrument `cg_state_t.last_iters`) that CG iterations
      grow ~linearly in n.
- [ ] Jacobi (then low-rank/Nyström if needed) preconditioner.
- [ ] Warm-start CG from previous EM iteration's `r̃`.
- [ ] Sub-linear `max_iter` cap; smarter Hutchinson; sub-quadratic ARG matvec
      (Lee 2025 / Zhu 2026 direction); revisit Phase 5.5 OpenMP after
      preconditioner reshapes the profile.

### Block 5 — Publication readiness (last)

Was old Phase 9; unchanged.

- [ ] Documentation: README + `sargent_help.txt` consolidation.
- [ ] Automated tests / CI for `test/lib/check_tree_matvec`, `smoke.sh`.
- [ ] Dead-code / consistency pass.
- [ ] Final Nurdan-pipeline end-to-end check.

### Parking lot — deferred items (still relevant, not now)

Items from retired phases preserved for future revisit. Each links to
where the original detail lives (below in this doc).

| Item | What it is | Why deferred | Revisit when |
|---|---|---|---|
| **Phase 6b open follow-ups** (unobs-MAF-conditional already done; remaining: `predict.py` h²_unobs reporting, `sv_sweep.sh` aggregator, scaled-t SV, per-effect-scale variation) | Refinements of the SV-flavoured simulator | Phase 10 verdict made the SV/missing-h² story a *localization* story, not a *recovery* story; the original SV-sweep design is partly superseded. | If localization claim broadens to specifically "ungenotyped-SV regions" and we need richer SV simulations to back it. |
| **Phase 6c — methods comparison on SV sims** (`sv_sweep.sh`, aggregator) | A within-pop methods comparison on the SV grid | Superseded; merged into the broader claim-driven comparison in Block 1.1 / 3. | If reviewers ask for an SV-specific comparison panel. |
| **Phase 6d — Nurdan R comparison** | Side-by-side benchmark of SARGENT (C) vs Nurdan's R implementation | Lower priority than the methodological story; depends on availability of Nurdan's code. | If accuracy parity needs to be demonstrated (e.g. for the manuscript appendix). |
| **Phase 7 — per-window localization (MoM-based)** | Was the localization plan; partly merged into Phase 10 / Block 1.2; the MoM module was deleted (2026-05-20). | Phase 10 supersedes; per-tree resolution is structurally limited (claim C2 falsified). | Unlikely to revisit; absorbed into Block 1.2 + Phase 11 if triggered. |
| **D2 (CLAIMS) — distinguishing polygenicity regimes** | A 2–3-regime polygenicity inference test | Low priority; π_k is biased (D1) so this would be a downstream test of a biased estimator. | After Block 1.3 fixes the absolute calibration. |
| **G4 (CLAIMS) — robustness to inferred ARGs** | Sensitivity test under `tsinfer`/`tsdate`-inferred ARGs | Manuscript heuristic claim that we'd want to back empirically eventually | Before any real-data run that uses inferred ARGs. |
| **H4 (CLAIMS) — extra comparison baselines** (mkLMM, LMM-winGRM, hapla) | Manuscript currently lists 9 comparators; only 3–4 are essential | Implementation effort high; story can be made with 2–3 baselines (BayesR, GBLUP, ARG-BLUP-uniform plus one external). | If reviewers ask for a specific comparator. |
| **Phase 5.5 leftovers** — smarter Hutchinson trace; OpenMP follow-up | Performance refinements after the preconditioner | Lower priority than the preconditioner pair (in Block 4). | Once Block 4 lands and the profile is re-shot. |
| **B5 (CLAIMS) — Neandertal/archaic prediction** | Striking but un-validatable application | Validation is fundamentally hard (no Neandertal phenotypes); paper's rigor lives in Blocks 1.1 + 1.2. | If the rest of the paper holds and we want a capstone application. |

---

## Current state of the code (as of 2026-05-19, commit `49a8109`)

**Algorithm modes already implemented:**

- **Sequential** (default): `var_fixed` (Stage 1 BayesR with V = σ²I) followed
  by `var_random` (Stage 2 per-tree REML EM with full V, plus post-hoc rho).
  Polyak-Ruppert averaging on σ²/ws² in `var_random`.  Matches the
  manuscript's "Sequential Inference Strategy" with the Stage-2 update form
  modified (Methods revision pending in `main.tex`).
- **Joint EM** (`--em <N>`): `var_joint` — interleaved coordinate descent.
  Constant V-aware precision `v_jj ≈ tr(V⁻¹)` (BOLT-LMM-style). No PR
  averaging (so the iter-final β, ws², σ², and var(r) are internally
  consistent). Used for prediction; biased for variance components.
- **Uniform-random** (`--uniform-arg`): orthogonal toggle on Stage 2 (in
  both sequential and joint modes) — replaces per-tree ws² with a single
  shared σ_u² on the uniformly-weighted kinship. ARG-BLUP equivalent.

**Auxiliary outputs:**
- `.beta.tsv` — per-SNP γ_jk and E[β_j]
- `.w2.tsv` — per-tree rho_m and E[w_s²]
- `.summary.json` — coherent record of h²_fixed, h²_random, h², σ² from
  iter-final state. **Authoritative variance-component output. Read this
  instead of parsing stderr.**

**Removed flags / dead code (all cleaned up):**
- `--horseshoe` — only changed post-hoc rho's φ_ref; cosmetic.
- `--pca` — was reserved and errored at runtime.
- `var_fixed_with_offset` — was for block-wise EM driver that has been
  replaced by interleaved var_joint.
- `extract_h2.py` — superseded by `.summary.json`.

**Test infrastructure**: `test/eval/` has `predict.py` (OOS prediction for
BayesR / GBLUP / ARG-BLUP-uniform / SARGENT), `run_suite.py` (reads
.summary.json), `extract_ws2_localization.py`, plus drivers
`overnight_sweep.sh`, `medium_tests.sh`. Test data lives in
`~/test-sargent/`.

**Known limitations** (will be addressed by this roadmap):
1. Variance partition is biased in both production modes: sequential
   over-attributes to fixed, joint EM over-attributes to random.
2. Polygenicity (π_k from BayesR) is biased by the same cannibalisation.
3. No MoM estimator: nothing in the code currently gives unbiased h².
4. No per-window / per-region variance breakdown.
5. Tree-iteration cost is O(n·S) per V-matvec (not the
   O(|ARG edges|) achievable with incremental tree-diff iteration).

---

## Phase 1 — Heritability reporting via sequential SARGENT (rename + reframe)

**Decision (2026-05-20)**: We accept that fixed-vs-random variance
decomposition is fundamentally not identifiable in this model (the K_X
and K_g matrices share structure by construction). Instead of trying to
estimate an inaccessible truth, we **reframe** the algorithm's outputs:

- $h^2_{\text{snp}}$ = variance explained by the BayesR fit on the typed
  SNP panel. Operationally defined as $1 - \mathrm{var}(r)/\mathrm{var}(y)$
  after the Stage 1 BayesR fit.
- $h^2_{\text{arg}}$ = **additional** variance captured by the ARG-based
  random effect after the SNP fit. "Missing heritability recovered."
- $h^2 = h^2_{\text{snp}} + h^2_{\text{arg}}$.

Under this framing, sequential SARGENT (Stage 1 then Stage 2) is the
*right* algorithm, not a buggy decomposer. The output of
`.summary.json` is exactly the answer the algorithm gives for this
sequential-prioritisation question.

**Naming convention applied across the codebase (2026-05-20)**:
- `sg->h2_fixed` → `sg->h2_snp`
- `sg->h2_random` → `sg->h2_arg`
- `sg->uniform_random` → `sg->uniform_arg`
- `--fixed-only` (`-F`) → `--snp-only` (`-S`)
- `--random-only` (`-R`) → `--arg-only` (`-A`)
- `--uniform-random` (`-U`) → `--uniform-arg` (`-U`)
- `simulate.py --h2-fixed` / `--h2-random` → `--h2-snp` / `--h2-arg`
- params.json fields `h2_fixed_target/realised` → `h2_snp_target/realised` etc.
- Test case labels: `fixed_only` → `snp_only`, `rand_only` → `arg_only`, etc.

### MoM comparison result (2026-05-20): DISCARDED

Built single-GRM MoM (closed-form 2-VC system on $K_g + \sigma^2 I$),
compared against sequential SARGENT on 3 seeds at full BED coverage
(f=1.0):

| Seed | Realised h²_snp + h²_arg | Sequential h² | Single-GRM MoM h² |
|---|---|---|---|
| s42 | 0.31 + 0.31 (≈0.6) | 0.87 | **0.05** |
| s43 | 0.36 + 0.36 (≈0.7) | 0.85 | **0.01** |
| s44 | 0.53 + 0.53 (≈1.0) | 0.72 | **0.04** |

MoM is consistently 10–20× too low. Root cause: $K_g$ is effectively
rank ≈ 2 (the deep coalescent branches dominate, so $K_g \approx
c\cdot\mathbf{1}_n\mathbf{1}_n^T$ with small structured residual). On
centred $y$, $K_g y \approx 0$, so $y^T K_g y$ provides almost no
signal. Single-GRM MoM with raw ARG kinship is not useful.

**Decision**: discard MoM. The `var_mom` C function, the `--mode
varcomp` CLI option, and all related test infrastructure have been
removed (commit pending). `--mode` now accepts only `predict` (default)
and `localize` (stub for Phase 7). The reframed sequential SARGENT
output via `.summary.json` is the canonical answer for heritability
reporting.

**If we ever need a real total-h² estimator**, the K_g matrix needs to
be **de-meaned** first (project out the constant component) so it
captures only the structured kinship deviations. Standard GREML
implementations do this. Worth considering for a future revisit.

---

## Phase 2 — `--mode` CLI structure (status: implemented in skeletal form)

**Goal**: future-proof the CLI for the localize extension while keeping
the current behaviour simple.

### CLI

```
sargent ...                       # default: --mode predict
sargent --mode predict   ...       # current main pipeline
sargent --mode localize  ...       # window-based per-region h² (Phase 7; stub)
```

### Default behaviour

`--mode predict` (default):
- Runs the sequential Stage 1 (BayesR on typed SNPs) → Stage 2 (per-tree
  REML EM on the ARG) pipeline. **Re-defaulted 2026-05-21** after the
  Phase 3.5 Nurdan smoke test exposed joint EM's σ²-collapse at scale
  (see Phase 4.5 below); the brief Phase 3 default flip was reverted.
- `--em N` (N ≥ 1) opts into joint variational coordinate descent.
  EXPERIMENTAL until Phase 4.5; runs fine on debug-scale data but
  degenerates on cases with ≳1000 local trees.
- `--sequential` is now a no-op kept for explicit-documentation and
  forward-compat with scripts written during the brief joint-EM default.
- `--snp-only` / `--arg-only` restrict to one stage and conflict with
  `--em N>0`.
- Outputs `.beta.tsv`, `.w2.tsv`, `.summary.json` with `h2_snp`,
  `h2_arg`, `h2`, `sigma2`.

`--mode localize` (Phase 7, not yet implemented):
- Per-window decomposition of $h^2_{\text{arg}}$ across the genome.
- Outputs `.local_h2.tsv` per-window summaries.

### Status

- [x] `--mode` flag in CLI, default `predict` (commit 2026-05-20)
- [x] `--mode varcomp` removed after the MoM evaluation showed it didn't help
- [x] Help text and `sargent_help.txt` updated
- [x] `--snp-only` / `-S` and `--arg-only` / `-A` replace the older
      `--fixed-only` / `--random-only` flags
- [x] `--uniform-arg` / `-U` replaces `--uniform-random`
- [ ] Update `test/eval/run_suite.py` and predict pipeline to use new
      `h2_snp` / `h2_arg` field names throughout (mostly done by sed
      pass on 2026-05-20; verify by running a small case)

---

## Phase 3 — Smaller-scale validation — DONE 2026-05-20

**Goal**: confirm sequential and joint-EM modes give expected results
before scaling up.  MoM dropped from this phase (module deleted after
empirical failure; see "MoM comparison result" above).

### Test sets

- **Debug scale**: n=300, L=100 kb, regime=sparse, 5 h² splits × 5 seeds.
  Regenerated post-RNG-fix.
- **h² total sweep**: at one BED coverage (f=0.5), vary h² split across
  snp_only … arg_only at total h²=0.8.

### Metrics

`predict`:
- OOS R² across {BayesR, GBLUP, ARG-BLUP-uniform, sequential SARGENT,
  joint-EM SARGENT}.
- Mean ± SD across seeds.

`localize` — deferred to Phase 7.

### Deliverables

- [x] Re-run overnight suite with cleaned-up code, regenerate `summary.md`
- [x] One-page comparison: sequential vs joint EM (`phase3_summary.md`)
- [x] Diagnostic figures (bed-coverage R² curve; h² decomposition with
  total bar) under `test-sargent/overnight/figs/`

### Findings (2026-05-20)

- Low-coverage regression is **gone**.  At f=0.03 both SARGENT variants
  recover R²≈0.58, matching ARG-BLUP-uniform; BayesR/GBLUP collapse to ~0.07.
- Joint EM > sequential by 3–7% R² across all sweep cells, and is robust
  to single-seed blow-ups that hurt sequential (e.g. seed 43 f=0.6).
- ~~Joint EM is now the default for `--mode predict`~~ — flipped on
  2026-05-20 based on this debug-scale result, then **reverted on
  2026-05-21** after the Phase 3.5 Nurdan smoke test exposed joint EM's
  σ²-collapse at scale.  See Phase 4.5 (joint EM redesign).
- Decompositions are unidentifiable in opposite directions (sequential
  over-attributes to h²_snp, joint EM over-attributes to h²_arg), but
  inferred h²_total agrees with realised truth within 5–15% for both
  on debug-scale data.

---

## Phase 3.5 — Postdoc-pipeline compatibility

**Goal**: make sargent + the eval scripts ingest input files produced
by Nurdan's simulation pipeline (under `~/test-sargent/nurdan/`) so
she can test SARGENT against her own data in parallel.

### File-format inventory (from `~/test-sargent/nurdan/...`)

- **Genotypes**: PLINK `.bed/.bim/.fam` (already supported) *and* a CSV
  variant (`all.genotype.csv`).  Use the BED.
- **Tree sequence**: `all.ts` (binary tskit format, our standard).
  Newick `.nwk` also provided but we don't need it.
- **Phenotype file**: `all.vc_*.pheno`, format
  `0 tsk_<i> y1 y2 y3 y4 y5` — five independent replicate phenotypes
  per individual; pick one column.  Current `read_phenotypes` takes
  the *last* whitespace-separated field, which would pick the 5th
  replicate.  Needs an explicit column selector.
- **Truth file**: `all.vc_*.truth.csv` with columns
  `phenotype,h2_fix,h2_rand,sigma2,h2_total,n_active_kernels,orthogonalized`
  — still uses old `h2_fix`/`h2_rand` naming; ingest layer translates.
- **SNP effects**: `all.vc_*.snp_effects.csv` (long form,
  `phenotype,snp_id,effect_size`).
- **Kernel weights**: `all.vc_*.kernel_weights.csv` (long form,
  `phenotype,kernel,gamma,pi`).
- **Train/test split convention**: indv 1–150 train, 151–200 test
  (different from our 2/3 split).

### Work items

- [ ] Add `--pheno-col <N>` (1-indexed) option to `sargent` for the
  multi-column phenotype file format.
- [ ] Add a small adapter script `scripts/import_nurdan.py` that
  converts a Nurdan case directory into our `case.*` layout so it can
  be fed to `predict.py` without modification.  Output: `case.phen`
  (one column), `case.trees`, `case.bed/bim/fam`, `case.params.json`
  (translated from `truth.csv`), `case.true_beta.tsv`,
  `case.true_ws2.tsv`.
- [ ] Confirm sargent runs on one Nurdan case end-to-end and the
  recovered h²_snp / h²_arg are within tolerance of `truth.csv`.

### Deliverables

- [ ] `scripts/import_nurdan.py` committed
- [ ] One example case under `test/eval/nurdan_smoke/` (paths only;
  data lives outside repo)
- [ ] `test/eval/phase3_5_summary.md` with R² numbers on at least
  one Nurdan case vs the manuscript's expected result

### Notes

- Highest external-dependency phase; do early to unblock the postdoc.
- Implementation should not touch the model code; pure ingest layer.

---

## Phase 4 — Profiling (desktop track)

**Goal**: identify the top three hot functions in sargent so we know
where to optimise.  Owned by the desktop Claude (macOS, `samply`).

### Tools

- **Primary**: `samply record ./build/src/progs/sargent ...` then open
  the resulting `.json.gz` profile in the Firefox Profiler UI.
- Backup: Xcode `Instruments.app` Time Profiler template if samply
  has issues with `dlopen`'d tskit symbols.
- Build with debug info: re-configure cmake with
  `-DCMAKE_BUILD_TYPE=RelWithDebInfo` so samply can resolve symbols.

### Targets to profile

Representative case: balanced h², bed_coverage=1.0, n=300, L=300 kb
(medium scale; see `test/sim/make_test_suite.py`). Both sequential
and joint EM modes (default invocation + `--sequential`).

Expected hot spots (to confirm):
1. `tree_cov_matvec` — called once per tree per V-matvec.
2. `compute_quad_all` — called once per outer iter, walks all trees.
3. `hutch_trace_all` — called once per outer iter, nprobe walks.
4. `cg_solve` — does ~30–50 V-matvecs per call.
5. `geno_read_variant` — called per SNP per outer iter in joint EM.

### Deliverables

- [ ] `test/profiling/profile_setup.sh` — simulates a fixed
  medium-scale case (seeded) and runs sargent under samply for each
  of {sequential, joint EM}
- [ ] `test/profiling/PROFILE.md` with measured breakdown of wall
  time across the candidate hot spots, screenshot of each flamegraph
- [ ] Recommendations for what to optimise, fed into Phase 5

---

## Phase 4.5 — Joint EM redesign at scale — DONE 2026-05-21

**Goal**: make `--em N` mode usable on cases with ~thousands of local
trees.  Originally joint EM σ² collapsed to ~0, h²_arg → 1, h²_snp → 0
on the Nurdan smoke case (n=200, S=2739, p=2488), independent of the
iteration count.

**What landed** (warm-start + damping + h²_arg fix + divergence revert):
- `--em N` now means "run sequential SARGENT first, then N polish
  iterations of var_joint starting from sequential's posterior".
  Joint EM is a *polish* step, not a cold-start algorithm.
- σ² update damped (ω=0.25, ratio clipped to [0.5, 2.0], floor raised
  to 0.001·var_y) -- prevents the Newton step from overshooting when
  V is near-singular.
- β update damped (ω=0.5) -- halves the cumulative drift from stale r̃
  across the β-sweep.
- h²_arg now reports `var(E[u | y, β])/var(y)` using the BLUP closed
  form `E[u] = r − σ²·r̃`, replacing the bogus
  `(var_r − σ²)/var_y` formula that saturated at 1 when β fit y.
- Divergence detection: if any polish iter shows `rel_dsig > 0.05`
  (BayesR mixture-snap signature) the iter is reverted to the start-
  of-iter snapshot and the polish stops.  Captures the genuine
  improvement up to the snap; bails before degeneration.

**End-to-end results** (vs sequential-only):

| Case            | Sequential                       | Phase 4.5 polish                              |
|-----------------|----------------------------------|-----------------------------------------------|
| s42_balanced    | h²=0.857  h²_snp=0.857  h²_arg=0  σ²=0.188 | h²=0.941  h²_snp=0.734  h²_arg=0.207  σ²=0.189 (revert at iter 9) |
| Nurdan          | h²=0.826  h²_snp=0.768  h²_arg=0.057  σ²=0.174 | h²=0.800  h²_snp=0.483  h²_arg=0.317  σ²=0.162 (clean 20 iters) |
| Truth (both)    | h²=0.80   h²_snp=0.40   h²_arg=0.40   σ²=0.20  | |

Phase 4.5 produces a more balanced SNP/ARG decomposition than sequential
in both cases.  On Nurdan it's now usable; on s42 the revert keeps
the genuine improvement (h²_arg=0.21 vs sequential's 0) and discards
the polish-induced explosion at iter 9.

**Diagnosis of the underlying issue** (for posterity / Phase 4.5 v2):

The σ² update `σ² *= A_0/B_0` is the textbook REML EM Newton step.
The dynamics are unstable when BayesR's mixture posterior snaps
discontinuously (a SNP suddenly favouring a different effect-size
component because the residual structure shifted).  When that
happens within an outer iter, the β-sweep produces a large
discontinuous β change, r changes a lot, the A_0/B_0 ratio computed
from the *pre-sweep* state is no longer a valid Newton step for σ²,
and σ² overshoots.  Once V's condition number is wrong the next
iter's CG produces a wild r_tilde, and the cascade is fast.

Sequential SARGENT avoids this because Stage 2 holds β fixed; only
ws² and σ² update, both with continuous likelihoods.

A real fix would smooth the BayesR mixture posterior (annealing
temperature, soft mixture-component weights) -- bigger surgery,
deferred.  The warm-start + revert combination is the pragmatic
Phase 4.5 fix that makes the existing `--em N` mode usable.

---

### Original problem statement (pre-Phase-4.5)

Originally joint EM σ² collapsed to ~0, h²_arg → 1, h²_snp → 0 on the
Nurdan smoke case (n=200, S=2739, p=2488), independent of the
iteration count.  Sequential SARGENT on the same case gives
h²_total=0.826 (within 3% of truth), so the failure was specific to
the joint-EM dynamics, not the underlying model.

### What we know from the diagnostic (2026-05-21)

Three independent issues, in stacked order of cause:

1. **β-sweep instability** — `r̃ = V^{-1} r` is computed at the start
   of each outer iter and held stale through the entire β-sweep.  With
   p ≫ n_train (2488 SNPs in LD), all SNPs see the same residual and
   cumulatively overshoot.  Confirmed: iter 1 alone produces residual
   variance 23× the phenotype variance.  A naïve fix
   (`β = (1-ω)β_old + ω·β_raw` with ω=0.5) stabilises the per-SNP
   trajectory but does NOT fix issue 2.
2. **σ² collapse** — the multiplicative update `σ² *= A_0/B_0` with
   `A_0 = ‖V^{-1} r‖²` and `B_0 ≈ tr(V^{-1})` (Hutchinson) can swing 2
   orders of magnitude per iter when V is near-singular, and finds a
   fixed point at σ² ≈ 0 once β is back under control.  Damping σ²
   would slow the descent but not change the fixed point.
3. **h²_arg reporting bug** — `h²_arg = (var_r - σ²)/var_y` reports as
   the leftover after subtracting σ², regardless of whether the random
   effect `u` is actually fitting anything.  When β=0 and σ²→0,
   h²_arg→1 by definition even with ws²_mean ≈ 1e-8 (i.e. `u` is
   essentially zero).  Independent of (1) and (2).

### Proposed redesign

Stacking patches on the current dynamics risks finding a fourth issue
underneath.  Better to rewrite the σ² and h²_arg paths from first
principles, mirroring what var_random does:

- **σ² update**: replace multiplicative `A_0/B_0` with a direct
  residual-variance update that conditions on the current u prediction:
  `σ² = ‖y - X β - E[u]‖² / n`.  This is what var_random does
  internally and is empirically stable on Nurdan-scale data.
- **h²_arg**: report `var(E[u])/var(y)` where `E[u] = Σ_s w_s G_s r̃`,
  not the residual subtraction.
- **β-sweep**: either (a) multiple inner β-sweeps with `r̃` recomputed
  between them, or (b) use the locally-updated `r` instead of the stale
  `r̃` (forfeits the V-correction but is stable; equivalent to
  classical BayesR with σ² that already accounts for u).
- **Warm-start option**: as a fallback, accept that joint EM is a
  *polish* step initialised from sequential's posterior, not a
  cold-start algorithm.  Sequential runs first, joint EM runs N≤10 iters
  on top.

### Acceptance criteria

- [ ] Joint EM converges (σ² stays positive, h²_arg < 1) on the
  Nurdan smoke case
- [ ] At debug scale, joint EM still wins by ≥3% R² over sequential
  (so we don't regress on the Phase 3 finding)
- [ ] At Nurdan scale, joint EM matches or beats sequential on h²_total
  recovery
- [ ] Re-default `--mode predict` to joint EM if both prior bullets hold

### Dependencies

- Phase 4 profiling (where is the per-iter cost?) — already done.
- Independent of Phase 5 (incremental matvec).
- Should be scheduled before Phase 6 (large-scale validation), since
  scaling up only matters once both modes work.

---

## Phase 5 — Incremental tree matvec

**Goal**: replace tree-by-tree O(n·S) traversals with O(|ARG edges|)
incremental updates, à la Palamara group / Lee et al. Single biggest
scalability win.

### Mechanism

Maintain $S[u]$ (subtree sum) and $\text{acc}[u]$ (accumulated path
sum) as global state across the tree iteration. At each tree-to-next-tree
step, tskit's edge-diff machinery tells us which edges were removed and
which were added. For each removed edge $(p, c)$, propagate $-S[c]$ up
the path from $p$ to the root, updating $\text{acc}$ along the way.
For each added edge, propagate $+S[c]$ similarly.

This keeps the matvec output $\text{out}[i]$ valid at every step.

### References to consult

- **Palamara group**: eGRM implementation (Fan, Bhatia et al.). Probably
  has a tree-incremental matvec. GPL-3 licensed; we may need to extract
  the algorithm not the code if our license differs.
- **Lee et al.**: tree-based heritability estimation, may have a similar
  implementation.
- **tskit's own examples** of incremental tree statistics (e.g.
  `tsk_treeseq_genealogical_nearest_neighbours`) for code patterns.

### Implementation steps

1. Add a paranoia mode: keep current `tree_cov_matvec` available and add
   a flag to compare incremental output to the naive output on each call.
2. Refactor `v_matvec` to use a stateful object that:
   - Initialises with full first-tree traversal.
   - Steps forward via edge-diff updates.
   - Returns updated $\text{out}[i]$ values.
3. Verify incremental matches naive on small ARGs (paranoia mode).
4. Once correct, remove paranoia mode (or keep behind a debug flag).
5. Same treatment for `compute_quad_all` and `hutch_trace_all`, both of
   which iterate over trees similarly.

### Estimated effort

1–2 weeks of careful work for the core, plus testing. Could be
significantly shortened by reusing an existing implementation.

### Deliverables

- [ ] Refactored `tree_cov_matvec` (incremental)
- [ ] Verified against naive on at least 100 ARG configurations
- [ ] Benchmark: speedup vs current code at debug / small / medium scale
- [ ] Updated `PROFILE.md` showing new bottleneck distribution

---

## Phase 5.5 — Algorithmic / parallelism follow-ups (deferred)

Captured 2026-05-21 after Phase 5 chunk 4 landed only 1.17–1.29× end-to-end.
The matvec hot path is near its arithmetic floor; the bigger gains live
above the matvec layer.  Not actioning these now — pick them up after
the chunk-4 re-profile says which one matters.

### Preconditioned CG

CG currently runs 50–100 inner iterations per V-solve.  V = Σ ws² G_s +
σ² I; a Jacobi preconditioner (M = diag(V)) is cheap to maintain
incrementally (we already have S[u] which gives diag(G_s); diag(V) is a
weighted sum across trees).  Hope: halve the iteration count.

- **Where it pays**: every V-solve in CG, ~50 matvecs per outer iter.
  Multiplies into hutch_trace_all (nprobe × CG_solve).
- **Risk**: preconditioner construction must stay O(n) per V-solve; if
  it becomes O(n × trees) it costs as much as it saves.

### Smarter Hutchinson trace estimation

Today: vanilla Hutchinson, nprobe Rademacher probes, variance O(1/√K).
Alternatives:

- **Hutch++** (Meyer-Musco-Musco-Woodruff 2021): combines a low-rank
  approximation with residual Hutchinson.  Order-of-magnitude variance
  reduction for matrices with rapid spectral decay (which V is, when
  ws² is sparse — i.e. most of the time for sargent).
- **XTrace / XNysTrace** (Epperly-Tropp-Webber 2024): adaptive variants
  with better tail bounds.

Both could plausibly drop nprobe from 30 to ~10 with equivalent
posterior variance.

### Better EM convergence

Joint EM (`var_joint`) and per-tree REML EM (`var_random`) currently
run a fixed max-iter.  Plausible wins:

- Track `‖θ_k - θ_{k-1}‖` and bail when stable.  Empirically the joint
  EM h² values are pinned by iter 15–20 even though we run 50.
- Polyak-Ruppert averaging is already in `var_random`; could add it to
  `var_joint`.

These cost only `if` statements and one extra norm computation per iter.

### OpenMP over per-sample walks

**Current state (Phase 5.5, shipped 2026-05-21):** A `#pragma omp parallel for`
is already in `tree_cov_matvec_inc_readout`, gated on `n_hap > 1000`.
Measured speedup on n=1000 (n_hap=2000): ~1.1–1.3× end-to-end with default
threads (observed via `top`: processes barely exceed 100% CPU even with 4+
threads set).

**Why so little gain — two structural problems:**

1. **Per-tree fork-join overhead.**  The parallel region is entered once per
   tree per matvec call (~2946 trees × ~50–200 CG iterations = millions of
   spawns).  Each spawn costs ~1–10 µs; the work per thread is comparable
   (~250 haps × depth-10 = 2500 multiply-adds).  Overhead dominates.

2. **Only readout is parallelised; advance is serial.**  `_advance` (edge-diff
   S[u] updates) has chained pointer dependencies through `our_parent[]` and
   cannot be trivially parallelised.  It accounts for a non-trivial fraction
   of walker time at low per-tree readout work.

**Correct fix:** hoist the `omp parallel` to `tree_walker_run_inc` so the
thread team is created **once per matvec call**, then reuse it across all
trees with `#pragma omp for` inside readout:

```c
#pragma omp parallel  /* once per matvec — keep team alive across trees */
{
    while (advance tree) {
        #pragma omp single  /* serial advance; team waits */
        tree_cov_matvec_inc_advance(inc, &tree);
        #pragma omp for     /* parallel readout; reuses team */
        for (hi = 0; hi < n_hap; hi++) { ... }
        #pragma omp single  /* serial visit callback */
        visit(gv, n, ...);
    }
}
```

Expected gain from eliminating per-tree fork overhead: approach the 1.83×
(8-thread) speedup on readout that was measured in Phase 5.5 benchmarks.
End-to-end gain still limited by the serial advance fraction, but
likely 1.4–1.7× realistically.

**Work plan (explicit deliverable)**:

- [ ] Restructure `tree_walker_run_inc` (src/lib/trees.c) so the
  `#pragma omp parallel` wraps the tree-iteration loop body.  Trees
  still advance sequentially (the edge-diff updates S[] / our_parent[]
  in order); use `#pragma omp single` to bracket `_advance` and the
  visit callback, and `#pragma omp for` to parallelise the per-sample
  readout.
- [ ] Drop the `if(n_hap > 1000)` threshold inside `_readout` (no
  longer needed once fork-join overhead is one-time per matvec).
- [ ] Validate via the harness (`test/lib/check_tree_matvec`) — output
  must be bit-identical to single-thread.
- [ ] Bench wall-time on s42_balanced (n_hap=400), medium (n_hap=2000),
  and Nurdan (n_hap=400, ~2700 trees).  Compare against the pre-hoist
  numbers in `PROFILE_PHASE5.md` / `PHASE6_LARGER_SCALE.md`.
- [ ] If realised speedup is below 1.3× on medium with 8 threads,
  re-profile to find the next bottleneck before declaring it done.

**Sequencing**: schedule after Phase 6 closes out (we have the baseline
numbers now to gauge against).

### σ² damping convention (var_fixed vs var_joint)

After Phase 4.5 (var_joint) and PR #15 (var_fixed) shipped, the two
functions use different σ² damping schemes.  This is a design choice,
not an inconsistency to "fix":

- **var_fixed**: σ² update is the MLE direct form `σ² ← ||r||²/n`.
  Bounded above by var(y) by construction.  Damping rule: full step
  on descent, half step on ascent.  Cheap and adequate because the
  proposal is already in a sane range.

- **var_joint**: σ² update is the Newton-style form `σ² *= A_0/B_0`,
  where `A_0 = ||V⁻¹r||²` and `B_0 = tr(V⁻¹)` (Hutchinson estimate).
  The ratio `A_0/B_0` is UNBOUNDED — when V is near-singular CG can
  produce a wild r̃ and the ratio spikes to 100+.  Damping rule:
  symmetric ω=0.25 multiplicative + per-iter ratio clipped to
  [0.5, 2.0] + floor at 0.001·var(y).  Asymmetric-only would allow
  arbitrarily large ascent steps when V is ill-conditioned.

Both schemes target the same root cause (BayesR mixture-snap cascade)
but the input dynamics differ.  Code comments in both functions
cross-reference this design note.

**Goal**: confirm SARGENT works at "moderately large" scales after the
incremental matvec + Phase 4.5 + OpenMP have all landed.

### Test sets

- **Medium**: n=1000, L=1 Mb, ~3000 trees, bed_coverage=1.0, balanced
  h².  Already simulated; desktop running the multi-thread × mode matrix.
- **Larger**: n=2000, L=2 Mb, ~6000 trees, bed_coverage=1.0, balanced h².
  Tentative; desktop will fork this if Medium results look clean.
- **Original target** (deferred): n=1000-2000, L=5-10 Mb (10k+ trees).
  Only feasible after additional perf work.

### Metrics

For each (case × mode × threads):
- Wall time end-to-end.
- h² values (snp/arg/total) vs realised truth in `params.json`.
- OpenMP speedup ratio (`OMP_NUM_THREADS=1` baseline vs default).
- Phase 4.5 polish behaviour (clean vs divergence-revert).

### Deliverables

- [ ] `test/profiling/PHASE6_LARGER_SCALE.md` from the desktop run.
- [ ] OpenMP investigation if the multi-thread speedup is below
  expectation (microbench predicted ~1.5-1.8× on n_haps=2000).
- [ ] Update this ROADMAP entry with measured numbers.

---

## Phase 6b — SV-flavoured unobserved-variants simulation (simulator done 2026-05-21)

> **Status 2026-05-29: PARTIALLY ABSORBED.** Core simulator done (PRs #12,
> #13, plus MAF-conditional carve-out in #17). Remaining open follow-ups
> (h²_unobs reporting in predict.py, sv_sweep aggregator, scaled-t SV,
> per-effect-scale variation) are moved to the **Parking lot** above —
> Phase 10's verdict (`ws²` as localization, not SV recovery) made these
> lower priority. Revisit if Block 1.2 motivates richer SV simulations.

**Goal**: missing-heritability simulation in which typed SNPs are
present at reasonable density *and* a separate class of unobserved
variants is generated with larger / more variable effect sizes
(stand-in for structural variants).  The expectation is that the ARG
random component carries the unobserved-variant contribution, giving
SARGENT a clearer win over BayesR than the bed-coverage-subsampling
scheme.

Causal classes:
- **Class S (SNPs)**: typed, BayesR-prior effect sizes
  (spike-and-slab, mostly small).
- **Class U (unobserved)**: not in the BED, larger effects from
  Gaussian or Laplace at `--unobs-effect-scale × √max(σ²_mix)`
  (default 5×).  Stand-in for SVs and rare large-effect alleles.

### What's done (PRs #12 + #13)

- `scripts/simulate.py` extended with `--h2-unobs`, `--n-unobs-causal`,
  `--unobs-effect-dist {gaussian,laplace}`, `--unobs-effect-scale`.
- Class U sites carved out of the polymorphic pool BEFORE causal_idx
  and force-excluded from the BED panel regardless of `--bed-coverage`.
- `params.json` carries `h2_unobs_target / h2_unobs_realised`,
  `n_unobs_causal`, `unobs_effect_dist`, `unobs_effect_scale`.
- Companion `--causal-fraction` flag (PR #12) defaults to 1% (SBayesR
  canonical), so SNP-causal scales with the BED-eligible variant count.

### Open follow-ups (deferred per user review)

- [ ] `predict.py` / `run_suite.py` recognising `h2_unobs_target` and
  reporting SV recovery alongside h²_snp / h²_arg.
- [ ] `test/eval/sv_sweep.sh` — sweep over h²_S / h²_U grid at fixed
  coverage and density.  Aggregates as the Phase 3 overnight sweep does.
- [ ] MAF-conditional SV sampling (real SVs are MAF-skewed; currently
  uniform draw from the polymorphic pool).
- [ ] Scaled-t SV distribution as a third option besides Gaussian
  and Laplace.
- [ ] Per-effect-scale variation (each SV draws its own scale rather
  than sharing one).

### Success criterion (validated in Phase 6c)

At "small h²_S, large h²_U, dense SNPs": SARGENT should beat both
BayesR-only and ARG-BLUP-uniform.  If it doesn't, that signals the
per-segment ws² structure isn't buying us anything beyond a single
ARG-wide variance component.

---

## Phase 6c — Method comparison on Phase 6b SV simulations

> **Status 2026-05-29: RETIRED (merged into Block 1.1 / Block 3).** The
> within-population SV methods comparison is superseded by the broader
> claim-driven comparison structure: Block 1.1 (PRS portability) covers
> the prediction comparison, Block 3 / `CLAIMS.md` H4 covers the
> manuscript-table baseline triage. Section retained below for context.

**Goal**: systematic comparison of SARGENT vs baseline methods on
the SV-flavoured simulator, to demonstrate (or refute) the polygenic
recovery story.

### Methods to compare

For each Phase 6b case (h²_S / h²_U grid):

- BayesR-only (SNP fit; current `predict.py` baseline).
- GBLUP (single variance component on standard kinship).
- ARG-BLUP-uniform (single ARG random component, no per-tree weights).
- **Sequential SARGENT** (current default).
- **Joint EM polish** (Phase 4.5 `--em N`).

Each compared on:

- Prediction R² (held-out test set).
- h²_snp / h²_arg recovery vs realised truth.
- h²_unobs recovery (specifically: how much of the SV-class variance
  does the ARG random component pick up?).
- Wall time.

### Deliverables

- [ ] Method-comparison harness (`test/eval/sv_sweep.sh`) — Phase 6b
  follow-up to extend.
- [ ] Aggregator script (`overnight_aggregate.py` equivalent for
  Phase 6c outputs).
- [ ] One-page summary (`test/eval/phase6c_summary.md`) with the
  comparison plot.

### Success criterion

SARGENT (sequential or polish) Pareto-dominates the simpler baselines
on the high-h²_U cases, and matches them on the low-h²_U cases.  If
it doesn't, this is the empirical falsification of the
per-segment-weights approach.

### Dependencies

- Phase 6b simulator (done) + the predict.py / sv_sweep.sh follow-ups
  listed above.

---

## Phase 6d — Comparison with Nurdan's R implementation

> **Status 2026-05-29: DEFERRED to parking lot.** Lower priority than the
> methodological story (Block 0 / 1) and depends on availability of
> Nurdan's code. Revisit if reviewer parity demonstration is needed.

**Goal**: head-to-head benchmark of SARGENT (C) against Nurdan's R
implementation on the same input cases.  Two axes:

1. **Accuracy**: do the two implementations produce comparable h²
   estimates and prediction R² on the same input?
2. **Speed**: SARGENT's wall time is currently slower than expected
   for a compiled C implementation.  Find out whether it actually
   beats the R code; if not, why.

### Work items

- [ ] Obtain a runnable copy of Nurdan's R code and its required
  R packages (likely `phasta` / `tskit-r` / custom).
- [ ] Identify a small set of canonical test cases that both
  implementations can run unchanged (likely the existing Nurdan
  smoke case is the natural starting point).
- [ ] Run both implementations on each case; record h² values,
  prediction R² (if applicable), and wall time.
- [ ] Tabulate side-by-side; document any algorithmic differences
  that might explain divergence.

### Hypotheses to test

- C should be faster than R for the matvec-dominated cost (chunk 4
  delivered 2× over Phase 4 baseline; R interpreter overhead per
  tree-walk should be much larger).
- If C is NOT faster, suspect: (a) we're doing more iterations than
  R is, (b) we're doing more probes / matvecs per iteration, (c)
  the R code has a numerical trick we're missing.

### Deliverables

- [ ] Side-by-side benchmark table in `test/eval/phase6d_r_vs_c.md`.
- [ ] If C is slower, identify the root cause and either fix it or
  document it.

---

## Phase 7 — Localization of missing heritability

> **Status 2026-05-29: RETIRED — superseded by Phase 10 verdict + Block
> 1.2.** The per-tree resolution claim is structurally limited (CLAIMS
> C2 falsified); region-scale localization is held (C1) and now lives in
> Block 1.2 (robustness sweep) ± Block 2 (Phase 11). The MoM module this
> phase originally relied on was already deleted 2026-05-20.

**Goal**: per-window estimates of h²_random vs h²_fixed across the
genome, highlighting regions where the ARG captures variance the BED
misses.

### Algorithm

Window-based MoM with one $\sigma^2_w$ per window:

$$V = \sigma^2_\beta K_X + \sum_w \sigma^2_w K^{(w)} + \sigma^2 I$$

Solve the (W+2)-dimensional moment system. For typical genome partitions
W ≈ 1000–3000.

Two implementations:
- **Full**: tr(K^(w) K^(w')) for all w, w' pairs — O(W²) traces, costly.
- **Diagonal-block approximation**: assume off-diagonal blocks ~0
  (windows are approximately independent at Mb scale). Decouples to
  O(W) traces.

The diagonal-block version is likely sufficient and dramatically
cheaper.

### Outputs

- `.local_h2.tsv`: tab-separated per-window report with columns
  `window_start`, `window_end`, `n_trees`, `h2_snp`, `h2_arg`,
  `h2_missing` (where missing = random − fixed contribution at this window).

### Estimated effort

1–2 weeks after Phase 5 (incremental matvec) is in place.

### Deliverables

- [ ] `src/lib/mom.c` extended to per-window MoM
- [ ] CLI: `--mode localize`
- [ ] Test: simulator generates phenotype with known per-window
  contribution; localize recovers it
- [ ] Visualisation script (Python, optional) producing a Manhattan-like
  plot of local missing heritability

---

## Phase 8 — Manuscript update

> **Status 2026-05-29: REPLACED by Block 3.** The original "bring the
> manuscript in line with Phases 1–7" framing presumed incremental
> updates; the Phase 10 verdict + CLAIMS table audit require major
> restructure (drop claims A2/A3/B2/C2/G1; reframe around localization +
> portability). Block 3 is the active version.

**Goal**: bring the manuscript in line with all the features and
findings from Phases 1-7 once the validation rounds are complete.

### Sections needing attention

- Joint EM mode: document the Phase 4.5 warm-start design, the BayesR
  mixture-snap diagnosis, and the divergence-revert safety net.  Note
  the known limitation (cold-start joint EM is unreliable at scale).
- Identifiability framing: capture the "h²_snp = SNP-prioritised fit,
  h²_arg = additional ARG-recovered variance" interpretation that
  Phase 1 settled on.
- SV-flavoured simulator (Phase 6b) as a test bed for missing-
  heritability recovery; results from Phase 6c if they support the
  per-segment-weights design.
- Method-comparison story from Phase 6c (SARGENT vs BayesR / GBLUP /
  ARG-BLUP-uniform).
- Drop the MoM section (module deleted on 2026-05-20; never shipped).
- Localization section (only after Phase 7).

### Dependencies

- Should follow Phase 6c (method comparison) and Phase 6d (R-vs-C
  benchmark) so the empirical claims are grounded.
- Phase 7 (localization) can be a follow-up paper if too much to
  fit in this revision.

---

## Phase 9 — Cleanup before publication

> **Status 2026-05-29: still applies; renumbered as Block 5.** Content
> below remains the working checklist.

**Goal**: code quality and reproducibility pass before the manuscript
goes out / the code goes public.

### Items

- [ ] Documentation: README.md + sargent_help.txt currently fragmented;
  accumulated flags (`--em`, `--sequential`, `--pheno-col`,
  `--causal-fraction`, `--h2-unobs`, `--unobs-effect-{dist,scale}`,
  `--n-unobs-causal`) need a coherent top-level summary.
- [ ] Automated test / CI: `test/lib/check_tree_matvec` and `smoke.sh`
  exist but don't run on push.  GitHub Actions workflow that runs at
  least the smoke harness + a tiny `--sequential` smoke test.
- [ ] Dead-code / consistency pass: months of fast development have
  left some dead branches and inconsistent naming.  One careful read
  to catch them.
- [ ] Verify Nurdan's pipeline still works end-to-end with
  `import_nurdan.py` after all the post-Phase-3.5 changes.

---

## Cross-cutting concerns

### Manuscript revisions

In progress (Adam editing manually):
- [x] Revision 1: Stage 2 algorithm description updated to reflect
  per-segment REML EM (in `main.tex`, uncommitted).
- [x] Revision 3: PR averaging documented (in `main.tex`, uncommitted).
- [ ] Joint EM mode (Phase 4.5) — see Phase 8.
- [ ] MoM section: drop (module deleted 2026-05-20).
- [ ] Localization section — Phase 7.
- [ ] SV simulator + method comparison — Phase 8.

### predict.py `fit_blup` upward h² bias — audit (queued 2026-05-28)

Phase 10 surfaced that `predict.py`'s `fit_blup` (single-component multiplicative-
EM REML used for GBLUP / ARG-BLUP-uniform / SARGENT evaluation) systematically
**over-estimates `h²` by ~+0.25** on the n=300 admixed cases (truth ≈ 0.50, fits
0.75–0.81). The bias is independent of kinship choice (uniform ≈ oracle ≈
estimated) — so it is a `fit_blup` problem, not a per-segment-weights problem.

This contaminates every prediction-eval h² number we have ever reported
(Phase 3 overnight, Phase 6 medium-scale, Phase 6b/6c). SARGENT's own
`--arg-only` `h²_arg` is much less biased (~+0.06 to +0.10), suggesting the
C code's per-tree EM + PR averaging is doing something better — leading
candidates for the cause are (a) ML- vs REML-style score (no projection P
to remove the fixed-effects subspace beyond y-centering), (b) the
dominant rank-1 component of the aggregate ARG kinship interacting poorly
with the single-K parameterisation, (c) too few EM iterations (`max_iter=50`)
under bad conditioning.

- [ ] Reproduce on a clean known-truth case (simple GRM from random
      genotypes, no admixture) — does the bias persist?
- [ ] Implement the proper REML projection `P = V⁻¹ − V⁻¹X(X'V⁻¹X)⁻¹X'V⁻¹`
      (with X = intercept), compare against the current centred-`y` ML-style.
- [ ] Sweep `h2_init` and `max_iter` — fixed-point vs convergence?
- [ ] Cross-check against a reference REML (e.g., `gcta --reml` on the same
      kinship + phenotype).
- [ ] Annotate the historical eval write-ups with the bias direction once
      understood, and fix `fit_blup` (or replace) before any new prediction
      tables go into the manuscript.

### Test data conventions

- `~/test-sargent/` is the data directory (gitignored, separate from
  the repo).
- Simulator: `scripts/simulate.py` with `--regime sparse` (default,
  3-component spike-and-slab) and `--bed-coverage` (RNG-independent
  from phenotype generation since commit 0b6b4b3).
- Test suite drivers: `test/eval/overnight_sweep.sh` and
  `test/eval/medium_tests.sh`.

### Open questions to revisit

- Whether to ship MoM alongside an external dependency (e.g.,
  pcgc-s, RHE-mc) — currently planning to ship our own implementation
  for self-containedness.
- Whether to expose the Hutchinson `nprobe` setting per mode — MoM may
  want a different default than the iterative algorithms.
- How to handle missing genotypes (currently the BED reader expects
  complete data).

---

## Phase 10 — Oracle diagnostic: is per-segment `ws²` worth rescuing? — DONE 2026-05-28

**Outcome.** Two ceilings measured in the optimistic regime (admixed,
F_ST≈0.26, sparse 5%-active concentrated blocks, 5 seeds at n=300+100):

- **Prediction ceiling** (oracle − uniform R²) = **+0.023 ± 0.013**
  (paired SE 0.006; all 5 sparse seeds positive). Estimator gap = −0.008.
  → Per-segment weights are *not* a meaningful prediction lever even with
  a perfect estimator. Bounded headroom ≈ +0.02 R² total.
- **Localization ceiling** (idealized AUC at uniform-null V) = **0.844**
  (range 0.78–0.90). Current per-tree REML estimator = **0.74** mean
  (range 0.59–0.90). → Active blocks *are* identifiable under structure,
  and there is ~0.10 AUC of headroom for a better estimator.

Diffuse control: oracle = uniform to the digit on all 5 seeds. ✓

**Verdict (forces the Phase 11 reframe below):**

> **Per-segment `ws²` is a localization tool, not a prediction tool.**

Detail and discussion: `test/eval/phase10_oracle.md`.

**Goal (as completed)**: settle, cheaply and decisively, *before* any new
modeling, whether heterogeneous per-segment weights can ever beat ARG-BLUP-
uniform — and if so, whether the gap is a *ceiling* problem (the truth is
flat) or an *estimation* problem (the truth is heterogeneous but we can't
recover it).

### Why (diagnosis from the 2026-05-27 re-evaluation)

In every experiment to date, SARGENT ≈ ARG-BLUP-uniform. The corrected
explanation is **information starvation, not collinearity**: each local tree
spans a tiny slice of the genome and carries ~`h²/S` of trait variance, so
the per-segment Fisher information is `~1/S²` (the manuscript's `O(n²/S²)`
score), while the phenotype supplies only ~`n` independent constraints to
split across `S` weights → `~n/S` of information per weight, hopeless when
`S ≫ n`. This holds *whether or not* the `Gₛ` are collinear — individual
`Gₛ` are in fact rich and distinct (only adjacent trees are strongly
correlated, decaying over a few recombinations); it is the *aggregate*
kinship that is near-rank-1, not each `Gₛ`.

Two independent reasons SARGENT ≈ uniform, which this phase separates:
1. **Estimation starvation** — per-tree `ws²` is unrecoverable → fitted
   weights collapse to uniform-plus-noise.
2. **Architecture flatness** — we have only ever simulated *diffuse*
   architectures, for which the *true* `ws²` is ~constant and even an oracle
   ties uniform. The favorable regime (sparse large-effect regions in a
   large-Nₑ, weakly-structured ARG) has never been tested.

### Experiment

Simulate in a **large-Nₑ, weakly-structured** coalescent (rich, distinct
local trees). Two architectures:
- **(i) diffuse** — effects spread ~uniformly across segments (current default).
- **(ii) sparse large-effect** — a handful of contiguous regions carry most
  of `h²_arg`; rest ~null. This is the regime the method *should* favour.

For each, compare prediction (OOS R²) **and** local-`h²` recovery under:
- **oracle** `ws²` (the true per-segment weights from the simulation),
- **uniform** `ws²` (ARG-BLUP-uniform),
- **estimated** `ws²` (current per-tree REML EM).

### Decision rule

- **oracle ≈ uniform even on (ii)** → per-segment idea is *fundamentally*
  limited for prediction; retire it for prediction, keep only as coarse
  local-`h²` inference. Pivot the paper to the transfer/kinship story.
- **oracle ≫ uniform on (ii) but estimated ≈ uniform** → it is an
  *estimation* problem with a real ceiling to chase → proceed to Phase 11.

### Deliverables

- [x] `scripts/simulate.py`: `--ws2-architecture {diffuse,sparse}` knob
      (+ block params); also `--demography {panmictic,admixed}` for the
      optimistic regime. Emits true per-tree `ws²` to `.full.true_ws2.tsv`.
- [x] `test/eval/oracle_ws2.{sh,py}` — runs uniform / oracle / estimated
      paired ARG-BLUP via a single-pass multi-weighting kinship; reports
      OOS R², local-`h²` recovery, **and the localization ceiling**
      (idealized per-tree score `Aₛ/tr(Gₛ)` at the homogeneous-null V,
      `idealized_loc_auc`).
- [x] `test/eval/phase10_oracle.md` — ceiling table + verdict.

> Shares scaffolding with the cross-distance "linchpin" simulation
> (two-population coalescent, tunable split time) used to test whether the
> SNP+ARG hybrid beats plain ARG-BLUP *in transfer* — build them together.

---

## Phase 11 — `ws²` redesign as a **localization** tool: hierarchical centering + sparse fused lasso

> **Status 2026-05-29: gated on Block 1.2 outcome (now Block 2 in the
> Active plan).** Build only if Block 1.2's robustness sweep shows the
> current per-tree REML estimator's localization is inadequate under
> realistic / weak / diffuse architectures. Detailed design below remains
> the implementation plan if triggered.

**Goal (reframed by Phase 10).** Replace the per-tree REML EM with a
parameterization + prior whose purpose is **local-`h²` localization, not
prediction enhancement**. Concretely: close the ~0.10 AUC gap between the
current estimator (0.74) and the localization ceiling (0.84) on the sparse
admixed regime, while degrading *exactly* to ARG-BLUP-uniform when no
heterogeneity is supported — so prediction stays at the uniform baseline
(which Phase 10 showed is already near the +0.023 R² ceiling).

The prediction ceiling itself (≤ ~+0.02 R²) is not the target: chasing it
gains little. The product of Phase 11 is a data-driven local-`h²` track
(piecewise-constant, recombination-respecting), not a prediction PRS.

### Model

Multiplicative hierarchical centering (well-identified global scale + local
deviation), in log space for positivity:

```
log wₛ² = μ + θₛ ,        wₛ² = w̄² · exp(θₛ) ,   w̄² = exp(μ)
```

Prior: a **sparse fused lasso** on the deviation track θ along the genome:

```
-log p(θ) = λ₁ Σ_s |θₛ|  +  λ₂ Σ_s c_s |θₛ − θ_{s−1}|
```

- `λ₁` (sparsity on levels) pins most `θₛ` to exactly 0 → those segments sit
  at the global scale `w̄²`; also anchors θ at 0 so μ is identified.
- `λ₂` (sparsity on differences) makes nonzero deviations form a few
  contiguous blocks with sharp change-points.
- **ARG-aware fusion weights** `c_s ∝ ⟨Gₛ, G_{s−1}⟩` (or shared-edge count /
  1 − SPR-distance): fuse hard where adjacent trees are nearly identical,
  let θ jump cheaply *exactly* where the genealogy changes (recombination
  boundaries). The segmentation is emergent, not chosen.

Net prior shape: "genome-wide ARG variance `w̄²`, plus a handful of
contiguous regions of elevated/depressed local heritability" — exactly the
sparse-large-effect regime Phase 10 tests.

### Self-certifying property

`λ₁,λ₂ → ∞` gives `θ ≡ 0` = ARG-BLUP-uniform exactly. So the model is a
strict generalization that cannot lose to uniform under CV-tuned penalties,
and **the fitted number of nonzero blocks is itself the answer** to "does
per-segment weighting buy anything?" (0 blocks ⇒ uniform was right.)

### Inference (reuses existing Stage 2 machinery)

Per outer iteration, the expensive quantities are unchanged — one CG solve
for `r̃ = V⁻¹r` and Hutchinson for `Bₛ = tr(V⁻¹Gₛ)`:
1. **Global scale** μ: the existing `--uniform-arg` update
   `w̄² ← w̄²·(Σ Aₛ)/(Σ Bₛ)`.
2. **Deviations** θ: a proximal-gradient step
   `θ ← prox_{η Ω}(θ + η g)`, where `gₛ = (wₛ²/2)(Aₛ − Bₛ)` (the old REML
   score rescaled to log space) and `Ω` is the sparse-fused-lasso penalty.
   The prox has an exact `O(S)` 1-D solution (Condat / taut-string DP) —
   pure scalar work, no linear algebra. (MM alternative: reweighted-ℓ₂ on
   differences → a tridiagonal `O(S)` solve, approximate zeros.)
3. Recompute `wₛ² = w̄²·exp(θₛ)`; next CG/Hutchinson pass.

Per-iteration cost ≈ current Stage 2; only the cheap score→weights map
changes. Generalized MM (stochastic REML score, non-convex likelihood) —
same character as the current algorithm.

### Hyperparameters & output

- `λ₁, λ₂`: CV on OOS R², or — cleaner — choose so the **df = number of
  nonzero blocks** stays `≪ n`, directly enforcing the "`≲ n` estimable
  targets" requirement.
- Output: piecewise-constant `log wₛ²` track = a data-driven local-`h²`
  Manhattan plot whose blocks respect recombination boundaries (the
  manuscript `localh` figure in an *estimable* form), plus
  `K = Σ wₛ² Gₛ` for BLUP prediction that reduces to ARG-BLUP when no blocks
  are found.
- **Bonus (complexity):** operating on a few effective blocks instead of `S`
  segments also trims the per-iteration trace/score work — ties into the
  Phase 5.5 efficiency thread.

### Deliverables

- [ ] `--ws2-prior {reml,sparse-fused}` in `sargent`; the sparse-fused path
      with `--lambda1 --lambda2` (or a single strength + ratio).
- [ ] Sparse-fused-lasso prox (1-D, `O(S)`) in `src/lib/`.
- [ ] ARG-aware fusion weights `c_s` from the tree sequence.
- [ ] **Validate against Phase 10 ceilings on the sparse + admixed regime**:
      does sparse-fused close the localization AUC gap 0.74 → 0.84+ (the
      primary target), while staying at or above ARG-BLUP-uniform R² on
      both sparse and diffuse architectures (the safety floor)?
- [ ] `test/eval/phase11_sparse_fused.md` reporting:
      (i) localization AUC vs the Phase 10 ceiling, per seed;
      (ii) prediction R² vs uniform (must not regress);
      (iii) the fitted block track on a representative seed.

### Success criterion

- Primary: localization AUC closes ≥ 50 % of the estimator–ceiling gap
  (0.74 → ≥ 0.79 mean over 5 sparse + admixed seeds) **without** dropping
  prediction R² below ARG-BLUP-uniform.
- Secondary: on the diffuse control, the fused-lasso fit selects zero
  blocks (`λ`-tuned) and reproduces ARG-BLUP-uniform exactly.

### Future extension

Fully-Bayesian sibling: spike-and-slab on `θₛ` with a Markov (Potts) prior on
the indicators — the fused lasso is its convex, fast relaxation. Build only
if posterior uncertainty on the blocks is wanted.

---

## Phase 12 — Complexity recovery: the CG cost is O(n³)

> **Status 2026-05-29: renumbered as Block 4 in the Active plan; deferred
> until story is locked.** Required only if the manuscript makes
> biobank-scale claims. Diagnosis and implementation plan below remain
> current.

**Goal**: restore near-quadratic (ideally near-linear) scaling of the full-`V`
inference path (Stage 2 / joint EM).

### Diagnosis (2026-05-27)

The cubic is **not** a dense solve and **not** a solver regression — there is
no Cholesky/inverse/eigendecomposition anywhere, and `conjgrad.c` has had
`max_iter = n` with no preconditioner since the first commit. The cube is the
product of three separately-linear factors that all scale with sample size.
Per Stage-2 / joint-EM iteration:

```
(nprobe + 1) CG solves  ×  (CG iterations)  ×  (one V·x matvec)
```

| Factor | Scaling | Why |
|---|---|---|
| CG iterations | ~O(n) | `max_iter = n`, **no preconditioner**; unpreconditioned CG on an ill-conditioned `V` (small-σ² regime) drifts toward the cap. |
| matvec walks all samples | O(n) | `tree_walker_run_inc` per-sample walk-up (83–91% of wall time in Phase 5 profiling). |
| number of local trees S | ~O(n) | S grows ~linearly with n (772 @ n=300, ~3000 @ n=1000); matvec sums over all S. |

`CG_iters(n) × matvec(n·S) = n·n·n` → **O(n³)**. (Phase 6's observed ~40×
slowdown 300→1000 decomposed as `n_hap × n_trees × (CG/n) ≈ 43×` — that *is*
the cube.) Stage 1 (BayesR, `V = σ²I`) is genuinely linear with no CG; the
cube only appears once the full-`V` CG path is exercised. The profiler
attributing 83–91% to the matvec is misleading: the matvec is cheap *per
call* but is called `(nprobe+1) × CG_iters ≈ n·nprobe` times — cut the call
count and you cut the cube.

### Near-term (do these — recovers ~O(n²))

- [ ] **Verify first**: `cg_state_t.last_iters` already exists; log CG
      iteration counts at n=200 vs n=1000 to confirm they grow ~linearly
      (expected) vs stay flat. Tells us whether the preconditioner or the
      matvec is the bigger win.
- [ ] **Preconditioner** (highest leverage; the hook already exists, just
      `NULL`). Jacobi (diagonal of `V`) first; consider a low-rank /
      Nyström / deflation preconditioner if Jacobi is weak. Bounds CG
      iterations to O(1)–O(√κ) *independent of n* → removes a factor of n.
      (Supersedes the "Preconditioned CG" bullet under Phase 5.5.)
- [ ] **Warm-start CG** from the previous EM iteration's `r̃`. The solution
      moves slowly across EM iterations → far fewer iterations. Trivial;
      stacks with the preconditioner.

### Later (captured; recovers further toward near-linear)

- [ ] **Cap `max_iter` sub-linearly** (e.g. O(√n)) with the already-relaxed
      Hutchinson tolerance — cheap safety net once the preconditioner lands.
- [ ] **Sub-quadratic matvec**: `n·S ≈ n²` is intrinsic to an *exact*
      V-matvec; the Lee 2025 / Zhu 2026 fast tree-sequence matvec algorithms
      (the Phase 5 direction) get below `O(n·S)` via edge-diff structure →
      push the remaining factor toward near-linear.
- [ ] **Reduce CG-solve count / smarter trace estimation**: `nprobe` is a
      direct multiplier on Hutchinson; variance-reduced or low-rank trace
      estimators cut it (ties to the Phase 5.5 "Smarter Hutchinson" item).
- [ ] **Coarsening payoff**: the Phase 11 sparse-fused-lasso reduces `S` to a
      few effective blocks, trimming the per-iteration trace/score work — an
      incidental complexity win from the modeling redesign.
- [ ] **OpenMP over the per-sample walk** (Phase 5.5 item; was not delivering
      a real speedup — revisit after the preconditioner changes the profile).
