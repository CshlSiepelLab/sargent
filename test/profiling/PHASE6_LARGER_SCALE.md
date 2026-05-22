# Phase 6 — Larger-Scale Validation

**Date**: 2026-05-21/22  
**Branch**: desk/phase6-larger-scale  
**Binary**: `build/src/progs/sargent` (commit b7329c4, RelWithDebInfo, OpenMP enabled)  
**Cases**: medium (n=1000, L=1 Mb), seeds 1, 2, 3, 42  

---

## Cases

| param | value |
|-------|-------|
| n_indiv | 1000 |
| L | 1 Mb |
| n_causal | 20 |
| trees | 2946–3037 (per seed) |
| SNPs | 3249–3340 (per seed) |
| h²_snp (truth) | 0.379–0.450 |
| h²_arg (truth) | 0.379–0.450 |
| σ² (truth) | 0.200 |

---

## Wall time — sequential mode (`--sequential`)

Runs were launched in parallel (4 simultaneous processes, same machine).
OpenMP parallelism is partially effective (per-tree fork-join overhead —
see ROADMAP Phase 5.5 for the structural fix deferred to later).

| seed | wall time | user time |
|-----:|----------:|----------:|
|   s1 |  683 min  |   77 min  |
|   s2 |  648 min  |   72 min  |
|   s3 |  665 min  |   72 min  |
|  s42 |  546 min  |   61 min  |

The gap between wall and user time reflects 4-way resource contention;
effective per-run parallelism was minimal.  User time of ~60–77 min per run
represents the actual CPU work; this scales as O(n × n_trees × CG_max_iter)
relative to the Phase 5 case (n=300, ~105 s wall/user).

Scaling ratio: (2000/600) × (3000/772) × (1000/300) ≈ 43× — consistent with
observed ~40× increase in user time (77 min vs ~2 min estimated Stage 2 cost
in the Phase 5 case).

---

## Estimates vs truth — sequential mode

Sequential mode runs Stage 1 (BayesR) then Stage 2 (ARG REML) on residuals.
When h²_snp ≈ h²_arg, Stage 1 captures most total genetic variance, leaving
little ARG signal for Stage 2.  The bias below is expected for sequential mode
and motivates joint EM.

| seed | h²_snp est | h²_snp truth | h²_arg est | h²_arg truth | σ² est | σ² truth |
|-----:|-----------:|-------------:|-----------:|-------------:|-------:|---------:|
|   s1 |   0.831    |   0.379      |   0.017    |   0.379      |  0.161 |  0.200   |
|   s2 |   0.792    |   0.401      |   0.021    |   0.401      |  0.187 |  0.200   |
|   s3 |   0.799    |   0.450      |   0.020    |   0.450      |  0.161 |  0.200   |
|  s42 |   0.805    |   0.410      |   0.014    |   0.410      |  0.176 |  0.200   |

**h²_total** (h²_snp + h²_arg) estimated 0.812–0.848 vs truth 0.758–0.900 — the
total genetic variance estimate is in the right ballpark, but the split between
SNP and ARG components is badly wrong.  This is the known sequential-mode
pathology when both components are large and overlapping.

---

## Stage 1 (var_fixed) stability — divergence fix verified

Prior to this branch, n=1000 seed 42 diverged in Stage 1: σ² blew up from
0.12 to 4.72 over 31 iterations, driving h²_snp → −3.84.

Fix (commit b7329c4): beta damping (½-step per sweep) + asymmetric σ² damping
(downward: full step; upward: half step).  Analogous to the Phase 4.5 fix in
`var_joint` (commit bdce0e6).

Post-fix Stage 1 behaviour on all seeds:

| seed | Stage 1 iters | final σ² | final h²_snp | rel_dsig at stop |
|-----:|--------------:|---------:|-------------:|-----------------:|
|   s1 |   72          |  0.1788  |    0.831     |   9.6e-05        |
|   s2 |   72          |  0.2079  |    0.792     |   5.9e-05        |
|   s3 |   98          |  0.1790  |    0.799     |   9.8e-05        |
|  s42 |   96          |  0.1898  |    0.806     |   4.2e-05        |

Smooth monotone convergence; no blow-up on any seed.

---

## Stage 2 (var_random) stability

Stage 2 ran 100 EM iterations (VAR_RANDOM_MAX_ITER) with Polyak-Ruppert
averaging starting at burn-in=40.  h²_rand converged to ~0.014–0.022 on all
seeds — low, as expected given the residuals after Stage 1 have minimal ARG
signal.  No numerical instability observed.

---

## OpenMP inefficiency (observation)

All runs used `OMP_NUM_THREADS` default (system threads).  `top` showed most
sargent processes at or just above 100% CPU, indicating near-serial execution
despite OpenMP.  Root cause: `#pragma omp parallel for` is entered once per
tree inside `tree_cov_matvec_inc_readout`, spawning/joining a thread team
~3000× per matvec.  At n_hap=2000, work per team launch (~2000 ops) is
comparable to spawn overhead.

The correct fix — hoisting the parallel region to `tree_walker_run_inc` — is
documented in ROADMAP Phase 5.5 and deferred.

---

## How to reproduce

```bash
SARGENT=build/src/progs/sargent
BASE=~/test-sargent/phase6/medium_multi_seed

# Seeds 1, 2, 3 were simulated with:
#   python3 test/simulate/simulate.py --seed <S> --n 1000 --length 1e6 \
#       --h2-snp 0.4 --h2-arg 0.4 --out $BASE/s<S>

# Seed 42 reuses ~/test-sargent/phase5_scale/medium (same params)

for s in 1 2 3; do
  ( $SARGENT --silent --sequential --nprobe 10 \
      --outroot $BASE/s${s} \
      $BASE/s${s}.phen $BASE/s${s}.trees $BASE/s${s} \
      > $BASE/s${s}.log 2>&1
    echo "done $s"
  ) &
done

# s42
( $SARGENT --silent --sequential --nprobe 10 \
    --outroot $BASE/s42 \
    ~/test-sargent/phase5_scale/medium.phen \
    ~/test-sargent/phase5_scale/medium.trees \
    ~/test-sargent/phase5_scale/medium \
    > $BASE/s42.log 2>&1
  echo "done s42"
) &

wait
```
