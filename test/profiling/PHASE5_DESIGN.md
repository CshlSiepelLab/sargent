# Phase 5 — Incremental tree-matvec: design

**Goal**: replace per-tree O(n_nodes) two-pass traversal with O(edge_diff)
incremental updates across the tree sequence.  Profile (Phase 4) shows 96.6 %
of wall time is in `tsk_tree_postorder_from` + `tsk_tree_preorder_from` +
sargent's own loop in `tree_cov_matvec`; edge-insert/remove operations are
only ~2 %.  Replacing the traversal cost with the edge-diff cost is a
projected ≥30× speedup on the dominant operation.

## The algorithm we're optimising

`tree_cov_matvec(out, in, tree)` computes `out = G_s * in` where `G_s` is
the n_indiv × n_indiv kinship matrix induced by local tree s
(`G_s[i,j] = ∑ branch lengths shared above i and j`).

Two-pass tree traversal (current implementation in `src/lib/trees.c`):

```
1. Post-order: S[u] = subtree sum of in over samples below u
   leaf sample h:  S[u] = 0.5 · in[indiv(h)]
   internal u:     S[u] = Σ_children c  S[c]

2. Pre-order: acc[u] = Σ_ancestors v of u  branch_length(v) · S[v]
   root r:     acc[r] = 0
   non-root u: acc[u] = acc[parent(u)] + branch_length(u) · S[u]

3. Read off: out[i] = acc[hap_a(i)] + acc[hap_b(i)]
```

## Three callers

All three iterate trees externally and call `tree_cov_matvec` per tree.

| Caller            | What it accumulates           | State                     |
|-------------------|--------------------------------|---------------------------|
| `v_matvec`        | `out += ws2[s] · G_s in`      | n-vector `out`            |
| `compute_quad_all`| `A[s] = in' G_s in`           | per-tree scalar (binned)  |
| `hutch_trace_all` | `B[s], B_0` via Rademacher    | per-tree scalar (binned)  |

Phase-5 refactor introduces a **stateful tree-walker API** so each caller
plugs in a per-tree callback while the walker handles iteration and
incremental state.

```c
typedef void (*tree_visit_fn)(double *gv,      /* G_s · in (output of matvec) */
                              int      n,       /* n_indiv */
                              int      s_orig,  /* original tree index */
                              int      sp,      /* pruned segment index */
                              void    *ctx);

void tree_walker_run(const tsk_treeseq_t *ts,
                     const int *hap_to_indiv,
                     int n_indiv,
                     int S_orig,
                     const int *prune_map,
                     const double *in,    /* matvec input vector */
                     double *gv,          /* scratch n-vector */
                     tree_visit_fn visit, void *ctx);
```

`v_matvec`, `compute_quad_all`, and `hutch_trace_all` all become thin
shells around `tree_walker_run` with custom `visit` callbacks.

## Incremental state — what we maintain across trees

We maintain `S[u]` and **directly the output contribution** `out_haploid[h]`
for each sample h (rather than materialising `acc[u]` for all nodes).

Define `out_h(h) = Σ_ancestors v of h  branch_length(v) · S[v]` for sample h.
Then `gv[i] = out_h(hap_a(i)) + out_h(hap_b(i))`.

When an edge (p, c) changes:

- `S[c]` itself is unchanged (subtree of c is the same; only its attachment
  point moves).  This is the key insight that makes incremental viable.
- `S[anc]` for ancestors `anc` of `p` changes by `±S[c]` (subtree of c is
  no longer/now under them).
- `out_h(h)` for samples h **below c** changes by `± Δacc(c)`, where
  `Δacc(c)` is the difference in `Σ branch_length(v) · S[v]` along the
  path from c up to (and not including) the LCA of old-parent and
  new-parent.

So per edge change we update:
- `S[u]` along the ancestor path of `p`        →  O(depth)
- `out_h(h)` for samples below c              →  O(samples_below_c)

Total work across the tree sequence:
`O(num_edge_changes · (depth + avg_samples_below_changed_edge))`,
which is empirically << `O(num_trees · n_nodes)`.

## Edge-diff iteration

tskit's `tsk_tree_first` / `tsk_tree_next` machinery already exposes the
edge diffs internally via the `edge_in`/`edge_out` lists per tree boundary.
We'll use those to drive incremental state, then visit the per-tree
callback once `S[]`, `out_h[]` are in a consistent state for tree s.

## Validation strategy — paranoia harness

This is the highest-risk piece of work in the project.  Incremental update
correctness is notoriously hard to get right, and silent bugs would show
up as subtle inference drift, not crashes.

The plan:

1. Add `tree_cov_matvec_naive` (the current implementation) and a new
   `tree_cov_matvec_incremental` (the new one); both are exposed by the
   library, side by side, throughout development.
2. Add a debug-flag `--paranoid-matvec` to sargent that calls both
   implementations for every matvec and asserts agreement at `L_inf < 1e-8`.
3. Run Phase 3 sweep before and after; numerical agreement to machine
   precision is the success criterion.
4. Once stable across at least the full Phase 3 sweep + one Nurdan run,
   drop `tree_cov_matvec_naive` from the hot path (keep it as a
   `--paranoid-matvec` cross-check or move to test code).

## Scope plan, in landable chunks

### Chunk 1 — Walker API refactor, no behaviour change

- Add `tree_walker_run` calling the existing `tree_cov_matvec`
  internally per tree.
- Rewrite `v_matvec`, `compute_quad_all`, and `estimate_g_scale` to
  call `tree_walker_run` with the appropriate visit callback.
- Smoke + regression: sargent output bit-identical to pre-refactor.

`hutch_trace_all` is **deliberately left untouched** in this chunk:
its inner loop interleaves per-tree matvecs across multiple probe
vectors, so refactoring it through the single-`in` walker would change
floating-point accumulation order and break the bit-identical
acceptance criterion.  Chunk 2 either generalises the walker to a
batch-`in` API or accepts the FP drift.

**Acceptance**: sargent output bit-identical to pre-refactor on the
s42_balanced debug case for both sequential and joint-EM-50 modes.

### Chunk 2 — Paranoia harness + incremental `S[]`

Two parts; ship together.

- **Reference implementation**: add `tree_cov_matvec_pairwise` that
  builds G_s explicitly (O(n² · n_nodes) per tree, fine for the
  small cases we test against) and multiplies it against `in`.
- **Harness**: a test driver (under `test/lib/` with its own tiny
  build target) that loads a small tree sequence and asserts
  `tree_cov_matvec ≡ tree_cov_matvec_pairwise` at machine precision,
  for 100 random `in` vectors across all trees.  This is the trust
  baseline before we add incremental.
- **Incremental S[u]**: maintain `S[u]` across `tsk_tree_next` using
  the edge-diff machinery; recompute `out_h[]` via the existing
  pre-order pass per tree.  Modest speedup (saves the post-order
  pass).

**Acceptance**: harness passes for both the naive (now reference)
implementation and the new incremental-`S[]` version on a 100-tree
case.

### Chunk 3 — Incremental `out_h[]`

- Update `out_h[]` per edge diff using the Δacc rule.
- Drop the per-tree pre-order pass entirely.
- This is the >30× speedup.

**Acceptance**: paranoia harness passes on the Phase 3 multi-seed
sweep (40 cases × multiple matvecs each), and end-to-end sargent
output matches naive to `L_inf < 1e-10`.

### Chunk 4 — Drop naive from hot path, profile

- Keep naive only behind a `--paranoid-matvec` debug flag.
- Re-profile with `samply`; expect tskit traversal cost to fall from
  70 % to <5 % of wall time.

## Open questions for later

- **`compute_quad_all` and the prune_map**: when an edge change updates
  `S[anc]` along an ancestor path, the affected pruned-segment scalar
  `A[sp]` for an old tree that has already been visited is, by
  construction, already final.  But during the visit of tree s we must
  finish processing all incremental updates before reading off
  `A[sp] += in' gv`.  The walker API enforces this ordering naturally.
- **`v_matvec` warm-start within CG**: if CG calls `v_matvec` with
  nearby `in` vectors, can we maintain `S[u]` across calls?  No — `S[u]`
  depends on `in`, which changes between matvecs.  We re-initialise at
  every matvec call.  The cost we save is purely the per-tree traversal
  amortised across the tree sequence; not across matvec calls.
- **Pruned tree handling**: with `VAR_N_SEGMENTS = 0` we have S_pruned ==
  S_orig (no pruning).  Incremental works the same way regardless; the
  prune_map just changes which scalar bin we accumulate into.

## Files changed (chunk 1)

- `src/lib/trees.c` — add `tree_walker_run` next to existing functions.
- `src/lib/variational.c` — `v_matvec` and `compute_quad_all` rewritten
  as walker callbacks.
- `include/trees.h` (new or expanded) — walker API declaration.
- `test/profiling/PHASE5_DESIGN.md` — this file.
