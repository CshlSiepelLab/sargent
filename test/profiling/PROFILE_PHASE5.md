# Phase 5 Profile — Post-Chunk-4 Re-profile

**Date**: 2026-05-21  
**Branch**: main (commit f834678, post-chunk-4 incremental matvec)  
**Tool**: samply 1000 Hz sampling + `analyze_profile.py`  
**Case**: same Phase 4 case (n=300, L=300 kb, 772 trees, 871 SNPs, h²_snp=h²_arg=0.4, seed 42)  
**Binary**: `build/src/progs/sargent` (-DCMAKE_BUILD_TYPE=RelWithDebInfo)

---

## Wall time

| Mode          | Phase 4   | Phase 5   | Speedup |
|---------------|----------:|----------:|--------:|
| sequential    |   220.2 s |   105.1 s |  2.10×  |
| joint EM-50   |   320.5 s |   160.9 s |  1.99×  |

(Phase 4 times are from `test/profiling/PROFILE.md`; same case, same flags.)

---

## Top functions by self-time — sequential

| %wall | function |
|------:|---------|
| 83.4% | `sargent::tree_walker_run_inc` |
|  4.6% | `libtskit::tsk_tree_postorder_from` |
|  3.7% | `libtskit::tsk_tree_preorder_from` |
|  3.3% | `sargent::tree_cov_matvec` |
|  0.9% | `libtskit::tsk_tree_remove_branch` |
|  0.9% | `libtskit::tsk_tree_insert_edge` |
|  0.8% | `libtskit::tsk_tree_remove_edge` |
|  0.7% | `libtskit::tsk_tree_insert_branch` |
|  0.4% | `sargent::v_matvec_visit` |
|  0.3% | `libtskit::tsk_tree_position_next` |

## Top functions by self-time — joint EM-50

| %wall | function |
|------:|---------|
| 90.8% | `sargent::tree_walker_run_inc` |
|  1.6% | `libtskit::tsk_tree_postorder_from` |
|  1.2% | `libtskit::tsk_tree_preorder_from` |
|  1.0% | `sargent::tree_cov_matvec` |
|  1.0% | `libtskit::tsk_tree_insert_edge` |
|  1.0% | `libtskit::tsk_tree_remove_branch` |
|  0.9% | `libtskit::tsk_tree_remove_edge` |
|  0.8% | `libtskit::tsk_tree_insert_branch` |
|  0.4% | `sargent::v_matvec_visit` |
|  0.4% | `libtskit::tsk_tree_position_next` |

---

## Diagnosis

Phase 5 chunks 1–4 successfully displaced the Phase 4 bottleneck.  In Phase 4,
70.8% of wall time was inside libtskit's `tsk_tree_postorder_from` +
`tsk_tree_preorder_from`; in Phase 5 those functions together account for only
8.3% (sequential) and 2.8% (joint EM), giving ~2× end-to-end speedup on this
772-tree case.

The new dominant cost is **`tree_walker_run_inc`** (83–91%), which is the
incremental walker's per-sample walk-up: for every tree in every matvec call,
each haplotype walks from leaf to root through `our_parent[]` accumulating
`bl[u] * S[u]`.  This is O(n_haps × avg_depth) per tree — comparable arithmetic
to the old tskit traversal, but now inside sargent's own code.  The remaining
~8% libtskit cost is the edge-diff bookkeeping (`insert/remove_edge/branch`,
`tsk_tree_position_next`) that drives the incremental state forward.

The 3.3%/1.0% still showing for `sargent::tree_cov_matvec` reflects
`hutch_trace_all`, which (by design in chunk 1) was not ported to the
incremental path and still calls the naive per-tree traversal once per probe
vector.  Batching hutch_trace_all across probes (chunk 3b) would eliminate
this cost and amortise the edge-diff overhead across K probes.

**Next targets in priority order:**

1. **Chunk 3b — batch hutch_trace_all**: walk the tree sequence once for all K
   probes instead of K times.  Eliminates the remaining 8.3%/2.8% libtskit
   traversal and amortises edge-diff updates.  Expected ~1.1–1.2× additional
   wall-time gain on this case.

2. **Per-sample walk-up in `tree_walker_run_inc`**: the 83–91% remaining cost.
   Options: (a) ancestor-path cache (see PHASE5_DESIGN.md — doesn't scale to
   biobank), (b) SIMD vectorisation of the walk-up inner loop, (c) thinning the
   `our_parent[]` pointer chain by caching a fixed-depth skip pointer.

---

## Phase 4 baseline (same case, for reference)

| %wall | function |
|------:|---------|
| 40.3% | `libtskit::tsk_tree_postorder_from` |
| 30.4% | `libtskit::tsk_tree_preorder_from` |
| 25.9% | `sargent::tree_cov_matvec` |
|  0.5% | `libtskit::tsk_tree_insert_edge` |
|  0.5% | `libtskit::tsk_tree_remove_branch` |

---

## How to reproduce

```bash
# Profile files (not committed; kept under ~/test-sargent/profiling/):
#   ~/test-sargent/profiling/phase5/profile_sequential.json.gz
#   ~/test-sargent/profiling/phase5/profile_joint_em.json.gz

# Analyze
python3 test/profiling/analyze_profile.py \
    ~/test-sargent/profiling/phase5/profile_sequential.json.gz \
    ~/test-sargent/profiling/phase5/profile_joint_em.json.gz

# Re-run profiles from scratch
SARGENT=build/src/progs/sargent
CASE=~/test-sargent/profiling/case
DATA_DIR=~/test-sargent/profiling/phase5
mkdir -p $DATA_DIR

samply record --save-only -o $DATA_DIR/profile_sequential.json.gz \
    $SARGENT --silent --sequential --nprobe 10 \
    --outroot $DATA_DIR/run_sequential \
    $CASE.phen $CASE.trees $CASE

samply record --save-only -o $DATA_DIR/profile_joint_em.json.gz \
    $SARGENT --silent --em 50 --nprobe 10 \
    --outroot $DATA_DIR/run_joint_em \
    $CASE.phen $CASE.trees $CASE
```

Symbol resolution is automated by `analyze_profile.py` (same directory).
