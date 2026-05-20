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
- Runs sequential SARGENT (Stage 1 + Stage 2) by default.
- Or joint variational coordinate descent if `--em N` is given.
- `--snp-only` / `--arg-only` restrict to one stage.
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

## Phase 3 — Smaller-scale validation

**Goal**: confirm that all three modes give expected results before
scaling up.

### Test sets

- **Debug scale**: n=300, L=100 kb, regime=sparse, 5 h² splits × 5 seeds.
  Already generated for prior work; regenerate after `simulate.py` RNG fix.
- **Small scale**: n=200, L=500 kb, 5 h² splits × 5 seeds.
- **h² total sweep**: at one BED coverage (f=0.5) and one h² split
  (balanced), vary total h² ∈ {0.3, 0.5, 0.8}. Confirms the picture
  doesn't shift at lower heritability.

### Metrics per mode

`predict`:
- OOS R² across {BayesR, GBLUP, ARG-BLUP-uniform, SARGENT} (existing
  `predict.py`).
- Mean ± SD across seeds.

`varcomp`:
- MoM h²_fixed, h²_random vs simulator realised, mean ± SD across seeds.
- Polygenicity: estimated π_k vs simulator pi_mix.
- Comparison row: MoM vs sequential vs joint EM, all three reported.

`localize` — deferred to Phase 7.

### Deliverables

- [ ] Re-run overnight suite with new code post-Phase 1, regenerate
  `summary.md` with MoM rows added
- [ ] One-page comparison: sequential vs joint EM vs MoM
- [ ] If MoM is competitive (≤10% bias across all simulated cases),
  declare Phase 3 done

---

## Phase 4 — Profiling

**Goal**: identify the top three hot functions in sargent so we know
where to optimise.

### Tools

- macOS: `Instruments.app` Time Profiler template (sample-based).
- Or `gprof` if we add `-pg` to CFLAGS in `build/`.
- Or just `clock_gettime()` instrumentation inside the candidate
  hot functions.

### Targets to profile

Representative case: balanced h², bed_coverage=1.0, n=300, L=300 kb
(medium scale). Both sequential and joint EM modes.

Expected hot spots (to confirm):
1. `tree_cov_matvec` — called once per tree per V-matvec.
2. `compute_quad_all` — called once per outer iter, walks all trees.
3. `hutch_trace_all` — called once per outer iter, nprobe walks.
4. `cg_solve` — does ~30–50 V-matvecs per call.
5. `geno_read_variant` — called per SNP per outer iter in joint EM.

### Deliverables

- [ ] `PROFILE.md` with measured breakdown of wall time across the
  candidate hot spots
- [ ] Recommendations for what to optimise (informed by phase 5 plan
  below)

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
- [ ] Demonstrated that prediction R² and MoM h² values don't shift
  with scale (within multi-seed noise)

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
