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

The per-sample walk-up `for hi in 0..n_haps: ...` in
`tree_cov_matvec_inc_readout` is embarrassingly parallel — each sample's
accumulation is independent.  At n_haps=2000 (medium) or biobank scale,
this is 4–8× on a multicore laptop, more on a workstation.

- Requires no new state; OpenMP `#pragma omp parallel for` over the
  outer loop, `reduction(+:out[:n_indiv])` for the accumulator.
- Build already links libomp.dylib on macOS.

This is the biggest theoretical multiplier currently sitting unused.
Worth doing *after* the re-profile if matvec stays in the top three.

---

## Phase 6 — Larger-scale validation (in progress)

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
