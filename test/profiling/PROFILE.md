# Phase 4 Profiling Report

**Date:** 2026-05-20  
**Branch:** desk/profiling  
**Build:** RelWithDebInfo (`/opt/homebrew/opt/llvm/bin/clang`)  
**Profiler:** samply 0.13.1 @ 1000 Hz  
**Profiles:** `~/test-sargent/profiling/profile_joint_em.json.gz` and `profile_sequential.json.gz`  

Load in Firefox Profiler: https://profiler.firefox.com/ → File → Load from file.

---

## Test case

```
n=300 individuals, L=300 kb, 772 trees, 871 SNPs
h²_snp = h²_arg = 0.4 (balanced), bed_coverage = 1.0, seed = 42
```

Generated via `profile_setup.sh` using `scripts/simulate.py`.

---

## Wall-time breakdown

Symbol resolution: `nm -n` on both `sargent` (base `0x100000000`) and
`libtskit.dylib`. Percentages are **self-time** (leaf-frame samples).

### Joint EM (`--em 50`), wall time ≈ 317 s

| % wall | Self samples | Function | Library |
|-------:|-------------:|----------|---------|
| 40.3% | 129,118 | `tsk_tree_postorder_from` | libtskit |
| 30.4% | 97,310 | `tsk_tree_preorder_from` | libtskit |
| 25.9% | 83,159 | `tree_cov_matvec` | sargent |
| 0.5% | 1,601 | `tsk_tree_insert_edge` | libtskit |
| 0.5% | 1,554 | `tsk_tree_remove_branch` | libtskit |
| 0.4% | 1,430 | `tsk_tree_remove_edge` | libtskit |
| 0.4% | 1,183 | `v_matvec` | sargent |
| 0.4% | 1,172 | `tsk_tree_insert_branch` | libtskit |
| 0.2% | 642 | `tsk_tree_position_next` | libtskit |
| < 0.1% | — | everything else | — |

**Group summary (joint EM):**
- tskit tree traversal (`postorder` + `preorder`): **70.7%**
- sargent `tree_cov_matvec` (own loop code + readout): **25.9%**
- tskit tree stepping (insert/remove edge/branch): **1.9%**
- `v_matvec` loop overhead: **0.4%**
- `cg_solve`, `hutch_trace_all`, `compute_quad_all`, `geno_read_variant`: **< 0.1% combined**

### Sequential (`--sequential`), wall time ≈ 220 s

| % wall | Self samples | Function | Library |
|-------:|-------------:|----------|---------|
| 40.1% | 88,349 | `tsk_tree_postorder_from` | libtskit |
| 30.7% | 67,598 | `tsk_tree_preorder_from` | libtskit |
| 25.6% | 56,452 | `tree_cov_matvec` | sargent |
| 0.5% | 1,125 | `tsk_tree_remove_branch` | libtskit |
| 0.5% | 1,083 | `tsk_tree_insert_edge` | libtskit |
| 0.5% | 1,071 | `tsk_tree_remove_edge` | libtskit |
| 0.4% | 796 | `tsk_tree_insert_branch` | libtskit |
| 0.3% | 746 | `v_matvec` | sargent |
| 0.2% | 518 | `tsk_tree_position_next` | libtskit |
| < 0.1% | — | everything else | — |

The profile is essentially identical between modes. Sequential is faster because it does fewer
total V-matvec calls (stage 1 BayesR + one CG pass in stage 2, no joint EM loop).

---

## What's surprising

**1. tskit's traversal API dominates, not sargent code.**
The expected hot spots (`cg_solve`, `hutch_trace_all`, `compute_quad_all`,
`geno_read_variant`) each account for < 0.1% of self-time. The profiler attributes 96.6%
of wall time to tree traversal: `tsk_tree_postorder_from` (40%) and
`tsk_tree_preorder_from` (30%) from libtskit, plus sargent's own loop overhead in
`tree_cov_matvec` (26%). The bottleneck is the **tskit API call overhead per tree per
matvec**, not the algorithmic structure of the EM iterations.

**2. Tree stepping is already cheap.**
The edge-diff updates (`tsk_tree_insert_edge`, `tsk_tree_remove_branch`, etc.) used when
advancing from tree to tree account for only ~2% of wall time. The expensive part is
rebuilding the post/pre-order traversal for every matvec on every tree, not the tree
adjacency bookkeeping.

**3. Hutch and CG are free wrappers.**
`hutch_trace_all` itself (excluding callees) registers at < 0.1%, and `cg_solve` is
similarly invisible. All their cost flows through `v_matvec` → `tree_cov_matvec` →
tskit traversals. The computation is entirely memory-bandwidth-bound tree traversal.

---

## Implications for Phase 5

The incremental tree matvec (Phase 5) is even more important than the roadmap
anticipated. The current implementation calls `tsk_tree_postorder_from` and
`tsk_tree_preorder_from` once **per tree per V-matvec call**. With S=772 trees,
nprobe=30 Hutchinson probes, and ~40 CG iterations per probe solve:

```
Total traversal calls ≈ 50 outer iters × (30 probes × 40 CG × 772) ≈ 46M traversals
```

The incremental approach maintains `S[u]` (subtree sum) and `acc[u]` (accumulated path
sum) as running state across the tree sequence. At each tree boundary tskit's edge-diff
machinery fires, updating only the changed edges — visible in the profile as the ~2%
edge-insert/remove cost. Phase 5 replaces the 70.7% traversal cost with that 2%
edge-diff cost, a potential **35× speedup** on the dominant operation.

**Phase 5 attack order:**
1. `tree_cov_matvec` in `trees.c` — replace the two tskit traversal calls with an
   incremental update driven by `tsk_tree_edge_diffs_init` / `tsk_tree_seek`.
2. Once `v_matvec` is incremental, `hutch_trace_all` and `compute_quad_all` benefit
   automatically since they call `v_matvec` (or `tree_cov_matvec` directly).
3. No changes needed to `cg_solve`, `var_joint`, `var_random`, or any EM logic.
