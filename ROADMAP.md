# SARGENT — Development Roadmap

Last updated: 2026-05-19.

This document captures the agreed plan for finalising SARGENT. The
priority order is **correctness first, scale later**. We finish each
algorithmic mode at small scale, validate it, then tackle the engineering
optimisations needed for biobank scale.

If this session is interrupted, restart from the next un-checked item.

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

Worth doing after Phase 6 validation is complete.

---

## Phase 6 — Larger-scale validation

**Goal**: confirm SARGENT works at "moderately large" scales after the
incremental matvec.

### Test sets

- **Larger scale**: n=1000, L=5 Mb (~10k trees), bed_coverage ∈
  {1.0, 0.3, 0.1}, balanced h², 3–5 seeds.
- **X-larger** (stretch goal): n=2000, L=10 Mb. Probably only feasible
  with the incremental matvec optimisation.

### Metrics

Same as Phase 3 for each mode, plus:
- Wall time per case for each mode (pre- and post-incremental optimisation).
- Memory footprint.

### Deliverables

- [ ] Larger-scale benchmark report in `summary.md`
- [ ] Demonstrated that prediction R² and h² values don't shift
  with scale (within multi-seed noise)

---

## Phase 6b — SV-flavoured unobserved-variants simulation

**Goal**: a more realistic missing-heritability simulation in which
typed SNPs are present at reasonable density *and* a separate class of
unobserved variants is generated with larger / more variable effect
sizes (a stand-in for structural variants).  The expectation is that
the ARG random component carries the unobserved-variant contribution,
giving SARGENT a clearer win over BayesR than the current
bed-coverage-subsampling scheme.

### Why this is interesting

Current simulator: all causal variants are common SNPs; missing
heritability is induced by dropping a fraction $f$ from the BED.  This
penalises BayesR but the "unseen" variants have the same effect-size
distribution as the seen ones.

Proposed: two distinct causal classes:
- **Class S (SNPs)**: typed, BayesR-prior effect sizes (small-to-moderate,
  spike-and-slab).
- **Class U (unobserved)**: not in the BED, larger and more dispersed
  effect-size distribution (e.g. Laplace or scaled-t with heavier tail);
  represents structural variants and rare large-effect alleles.

Heritability split: choose $h^2_S$, $h^2_U$, $\sigma^2$.  Class U is
exactly the heritability the ARG component must recover.

### Work items

- [ ] Extend `scripts/simulate.py`: new flags `--h2-unobs <h2>` and
  `--unobs-effect-dist {laplace,t,custom}` and `--unobs-effect-scale`.
- [ ] Class U sites are sampled from the tree-sequence variants but
  excluded from the BED before writing.
- [ ] `params.json` carries `h2_unobs_target` / `h2_unobs_realised`
  alongside the existing fields.
- [ ] Update `predict.py` to recognise the new params field (no
  algorithm change — SARGENT's job is to recover this without knowing
  the partition).
- [ ] Add `test/eval/sv_sweep.sh` — a sweep over $h^2_S / h^2_U$ at
  fixed coverage and density.  Aggregate as before.

### Success criterion

At the parameter combination "small h²_S, large h²_U, dense SNPs",
SARGENT should beat both BayesR and ARG-BLUP-uniform.  If it doesn't,
that's a signal that the per-segment $w_s^2$ structure isn't actually
buying us anything beyond a single ARG-wide variance component, and
we'd need to rethink the random model.

### Dependencies

- Sits after Phase 6 (large-scale validation) in the *committed*
  ordering, but can be developed earlier as desktop work.

---

## Phase 7 — Localization of missing heritability

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

## Cross-cutting concerns

### Manuscript revisions

In progress (Adam editing manually):
- [x] Revision 1: Stage 2 algorithm description updated to reflect
  per-segment REML EM (in `main.tex`, uncommitted).
- [x] Revision 3: PR averaging documented (in `main.tex`, uncommitted).
- [ ] Add description of the joint-EM mode (if we want to keep it
  documented).
- [ ] Add MoM section once Phase 1 is done.
- [ ] Add localization section once Phase 7 is done.

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
