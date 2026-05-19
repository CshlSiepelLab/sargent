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
- **Uniform-random** (`--uniform-random`): orthogonal toggle on Stage 2 (in
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

## Phase 1 — MoM variance component estimation

**Goal**: unbiased global h²_fixed, h²_random, σ² estimates that don't
suffer from the fixed/random cannibalisation. Also unbiased polygenicity
(BayesR π_k) by residualising y against the MoM random component before
the BayesR fit.

### Algorithm

Two random-effect model: $y = X\beta + u + \varepsilon$ with $u \sim
N(0, \sigma^2_g K_g)$, $\varepsilon \sim N(0, \sigma^2 I)$, and
$\beta \sim N(0, \sigma^2_\beta I)$ in the moment-equation approximation
($K_X = XX^T/p$).

Three moment equations from $E[y^T M y] = \mathrm{tr}(MV)$ for
$M \in \{I, K_X, K_g\}$:

$$
\begin{pmatrix}
\mathrm{tr}(K_X) & \mathrm{tr}(K_g) & n \\
\mathrm{tr}(K_X^2) & \mathrm{tr}(K_X K_g) & \mathrm{tr}(K_X) \\
\mathrm{tr}(K_X K_g) & \mathrm{tr}(K_g^2) & \mathrm{tr}(K_g)
\end{pmatrix}
\begin{pmatrix} \sigma^2_\beta \\ \sigma^2_g \\ \sigma^2 \end{pmatrix}
=
\begin{pmatrix} y^T y \\ y^T K_X y \\ y^T K_g y \end{pmatrix}
$$

Solve by least squares with non-negativity constraint on the variance
components (NNLS).

### Trace estimates

- `tr(K_X) = sum_j ||x_j||² / p ≈ n` (exact; no estimation needed).
- `tr(K_g) ≈ S · E[tr(G_s)/n] · n / S = E[tr(G_s)/n] · n` — actually we want
  the average per-tree trace, which the existing `estimate_g_scale` already
  computes during `sargent_new`.
- `tr(K_X^2)`, `tr(K_g^2)`, `tr(K_X K_g)`: stochastic Hutchinson estimates
  via z ~ Rademacher, computing z^T M z for the relevant M's.

Reuse existing CG + matvec infrastructure. Need `K_X v` matvec via the
geno reader, which we already have. Need `K_g v` matvec, which is
`v_matvec` with uniform ws²[s] = 1/S and σ² = 0.

### Polygenicity (the post-MoM BayesR refit)

After MoM gives $\hat\sigma^2_g$ and $\hat\sigma^2$:

1. Construct $V_{\text{mom}} = \hat\sigma^2_g K_g + \hat\sigma^2 I$.
2. Compute $\hat u = \hat\sigma^2_g K_g V_{\text{mom}}^{-1} y$ (one CG solve).
3. Form residualised response $\tilde y = y - \hat u$.
4. Run var_fixed on $\tilde y$ with $\sigma^2 = \hat\sigma^2$ held fixed.
5. Read off π_k from the converged Stage 1 state.

Output goes to `.pi.tsv` (per-component mixture proportions) and the
post-MoM γ_jk and β̂ to `.beta.tsv`. Total work: 1 CG solve + 1 var_fixed
run beyond MoM itself.

### Code structure

New file `src/lib/mom.c`:

```c
typedef struct {
    double h2_fixed, h2_random, sigma2;
    double sigma2_beta, sigma2_g;        // raw variance components
    int n_active_pi;                     // polygenicity
    double pi[MAX_K];                    // mixture proportions
} mom_result_t;

mom_result_t *var_mom(sargent_t *sg, int silent);
void mom_result_free(mom_result_t *r);
```

Driver in sargent.c reads `--mode varcomp` and dispatches to `var_mom`,
then `var_fixed` on the residualised y for polygenicity, then writes
`.summary.json` + `.pi.tsv` + `.beta.tsv`.

### Testing strategy

- Validate on the `overnight/multiseed` suite (5 seeds × 5 BED levels).
  MoM h² should be approximately constant across BED coverage levels (since
  it doesn't use β to estimate variance) — unlike sequential SARGENT's
  cannibalisation pattern.
- Compare MoM h² to simulator realised h². Expect within ~10% accuracy
  even at debug scale.
- Compare polygenicity (π_k) to simulator's `pi_mix`. Should be close.

### Estimated effort

3–5 days for global MoM + polygenicity, including a parallel
implementation in Python for validation against the C version.

### Deliverables

- [ ] `src/lib/mom.c` and matching header
- [ ] CLI: `--mode varcomp` dispatches to MoM + polygenicity
- [ ] Output: `.summary.json` augmented with MoM h² fields and
  `mode: varcomp`; new `.pi.tsv`
- [ ] Validation: MoM h² matches realised h² ± 10% across seeds in
  multiseed bed-coverage sweep
- [ ] Update `ROADMAP.md` to reflect completion

---

## Phase 2 — `--mode` CLI structure

**Goal**: separate three operational use cases with one binary and one
flag. Default is `predict` (most common).

### CLI

```
sargent --mode predict   ...   # joint EM, prediction-optimised  (default)
sargent --mode varcomp   ...   # MoM + polygenicity refit, unbiased h²
sargent --mode localize  ...   # window-based MoM (Phase 7)
```

### Default behaviour

`--mode predict` (default):
- Runs joint EM by default, unless `--fixed-only` or `--random-only` are
  specified (then sequential).
- Outputs `.beta.tsv`, `.w2.tsv`, `.summary.json` (h² values from joint
  EM, with caveat noted).

`--mode varcomp`:
- Runs `var_mom` for variance components.
- Then runs `var_fixed` on MoM-residualised y for polygenicity.
- Outputs `.summary.json` (MoM-based h²) and `.pi.tsv`.
- `.beta.tsv` is from the post-MoM BayesR (interpretable as fixed-effect
  posterior conditional on the MoM random component).

`--mode localize` (Phase 7):
- Runs window-based MoM.
- Outputs `.local_h2.tsv` per-window decomposition.

### Implementation notes

- The `--em` and `--fixed-only` / `--random-only` flags become invalid
  outside `--mode predict`. Sargent should error early with a clear message
  if a user combines them.
- `--mode predict --em 0` (no EM, sequential) remains available for users
  who specifically want the original two-stage algorithm.

### Deliverables

- [ ] `--mode` flag in CLI, default `predict`
- [ ] Validate that legacy invocations (no `--mode`) still work
- [ ] Update help text and `sargent_help.txt`
- [ ] Update test/eval/run_suite.py to pass the right `--mode`

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
  `window_start`, `window_end`, `n_trees`, `h2_fixed`, `h2_random`,
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
