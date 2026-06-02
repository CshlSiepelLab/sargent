# Block 2.1 — Allele-age architecture in `simulate.py`

Date: 2026-05-31.
Status: v1 implemented + qualitatively validated.

Adds the **WHEN** axis of architecture characterization (per ROADMAP
Block 2.5 framing). Lets the simulator produce traits whose per-variant
effect sizes systematically vary with allele age, so downstream
SARGENT estimators can be tested on age-stratified architectures.

## CLI additions to `scripts/simulate.py`

```
--age-tau-coef TAU      per-variant β² scales as (age / age_ref)^TAU
                        before h²_snp normalisation
                        TAU = 0     no age dependence (default; legacy)
                        TAU < 0     young variants larger (purifying-
                                    selection-driven; O'Connor / Zeng)
                        TAU > 0     ancient variants larger (drift /
                                    neutral background)
--age-tau-ref REF       reference age in generations (default 10000;
                        ≈ 4·Ne for Ne=1e4)
```

The total h²_snp is preserved by the downstream `scale_f` rescaling;
TAU only changes the *shape* of the per-variant effect-size distribution
over allele age.

## Per-variant age outputs

`<prefix>.true_beta.tsv` now has an `age` column (generations).
`<prefix>.params.json` includes `age_tau_coef` + `age_tau_ref`.

The age of a variant is the time of its mutation (from
`var.site.mutations[0].time`); falls back to the parent-node time on the
local tree if the mutation time is unknown.

## Per-age-bin h² helper: `test/eval/age_partition.py`

For analysis on the simulated traits. Reads `true_beta.tsv` + `params.json`
and produces per-age-bin h² partition (Σ_{j∈bin} β²_j · Var(x_j)).

```
python test/eval/age_partition.py --case-prefix CASE [--sargent-prefix EST]
    [--bins 3 | --quantiles 0.33,0.67 | --age-edges 1,1000,10000,100000]
    [--out PATH.json]
```

Default: equal-quantile bins of the observed causal-variant ages.

When `--sargent-prefix` is supplied, also computes an *estimated* partition
from sargent's `.beta.tsv` — Block 2.5.2 will use this to compare
recovered vs true age-stratified h².

## Validation (Block 2.1 smoke tests)

Same simulation (n=50, L=100kb, h²_snp=0.5, causal-fraction=0.05,
20 causal variants), three values of TAU, equal-quantile 3-bin partition:

| TAU | recent bin (n=7) | middle bin (n=6) | ancient bin (n=7) |
|---|---|---|---|
| −1.0 | **0.887** | 0.099 | 0.015 |
|  0.0 | 0.057 | 0.634 | 0.309 |
| +1.0 | 0.001 | 0.244 | **0.755** |

(Fractions of total β² · Var(x) per bin. With TAU < 0, recent variants
dominate; with TAU > 0, ancient variants dominate. TAU = 0 gives the
unaltered BayesR draw.)

The bin edges are common across the three (computed from the TAU=0 sim's
ages); the variance moves between bins as TAU changes, with the total
h²_snp held to ~0.5 by the downstream rescaling.

## Next applications

- **Block 2.5.3** (joint distribution age × polygenicity): variance-capture
  curves stratified by age bin. With TAU set, the WHEN-stratified
  variance-capture is the manuscript joint-distribution figure.
- **Block 2.2 (B), (C)** (ARG-only / hybrid LDSC analog): allele age becomes
  a non-positional annotation, the first example LDSC cannot do. The
  per-age-bin partition here is the validation of the ARG-LDSC framework.
- **Selection-effect demonstration**: compare TAU = −1 (selection-favoured
  young/rare) vs TAU = 0 (neutral) on a single architecture; can SARGENT
  detect the difference at the per-tree variance level?

## Reproduce

```bash
cd /tmp && mkdir b21 && cd b21
for tau in -1 0 1; do
  python ~/git-src/sargent/scripts/simulate.py \
      --n-indiv 50 --n-test 20 --length 100000 \
      --h2-snp 0.5 --h2-arg 0.0 \
      --causal-fraction 0.05 --bed-coverage 1.0 \
      --age-tau-coef $tau --age-tau-ref 1000 \
      --seed 1 --out smk_t$tau
  python ~/git-src/sargent/test/eval/age_partition.py \
      --case-prefix smk_t$tau --bins 3
done
```

## v2 — ARG-side age decomposition, n-scaling test (2026-06-01)

Per Adam's 2026-06-01 directive: "Should this test use SNPs at all or
would it be better to do it based on the ARG alone? Maybe that would
require some work to justify but it might still be needed — as you say,
the SNP part of the model doesn't really know about age." → The SNP-side
test (v1, above) put the age signal on BayesR effect sizes, but BayesR
has no age awareness. v2 moves the test entirely to the ARG side: the
per-tree variance `ŵ²_s` is binned by the *mean branch age of the
local tree* (`per_tree_mean_branch_age` in `arch_metrics.py`), and the
per-bin sum of `ŵ²_s · ℓ_s` is compared against truth.

### Setup

- Trait: `h²_arg = 0.5`, `h²_snp = 0`, uniform `ws²` (no artificial age
  stratification — we're testing whether the *natural* age distribution
  of variance is recovered).
- Demography: panmictic, single pop.
- Sims: `L = 100 kb`, `Ne = 1e4`. 3 sample sizes × 3 seeds = 9 cells.
- Bins: recent (≤5k gen), middle (5k–20k gen), ancient (>20k gen).
- Driver: `test/eval/age_decomp_nscan.sh` (BASE=`~/test-sargent/block2_1_v2_nscan`).

### Results — per-bin h² fractions (mean ± sd over 3 seeds)

| n | bin | truth frac | est frac | bias (est − truth) | mean trees |
|---|---|---|---|---|---|
| 300  | recent  | 0.022 ± 0.027 | 0.181 ± 0.314 | +0.160 | 29.7  |
| 300  | middle  | 0.719 ± 0.200 | 0.618 ± 0.340 | -0.101 | 194.7 |
| 300  | ancient | 0.259 ± 0.227 | 0.200 ± 0.255 | -0.059 | 32.7  |
| 600  | recent  | 0.009 ± 0.008 | 0.005 ± 0.007 | -0.004 | 22.3  |
| 600  | middle  | 0.639 ± 0.329 | 0.817 ± 0.192 | +0.178 | 249.0 |
| 600  | ancient | 0.352 ± 0.322 | 0.177 ± 0.187 | -0.174 | 33.7  |
| 1000 | recent  | 0.052 ± 0.045 | 0.090 ± 0.154 | +0.038 | 54.7  |
| 1000 | middle  | 0.856 ± 0.129 | 0.774 ± 0.215 | -0.083 | 223.0 |
| 1000 | ancient | 0.092 ± 0.088 | 0.136 ± 0.067 | +0.044 | 40.7  |

### Total L1 bias (Σ|est − truth| across bins, mean over seeds) vs n

| n | recent | middle | ancient | total L1 bias |
|---|---|---|---|---|
| 300  | +0.160 | -0.101 | -0.059 | 0.319 |
| 600  | -0.004 | +0.178 | -0.174 | 0.356 |
| 1000 | +0.038 | -0.083 | +0.044 | 0.165 |

Wall: 156 min. n=600 is paradoxically slightly worse than n=300 in
total L1 bias, but the spread across seeds at both n=300 and n=600 is
huge — this is small-sample noise across 3 seeds, not a real trend.
The clearer pattern: **n=1000 ≈ halves the L1 bias relative to n=300**.

### Interpretation

- The middle (5k–20k gen) bin carries the bulk of the recoverable
  variance under uniform `ws²` (since most segment-weighted ARG
  variance lives in this depth band at Ne=1e4); at n=1000 the
  estimator captures 0.774 vs truth 0.856 — close.
- The recent (≤5k gen) bin is intrinsically low-variance (truth ≤ 0.05
  across all n) and noisy on the estimator side; not interpretable
  per cell, but the bias trend ≤ 0.04 at n=1000 is acceptable.
- The ancient (>20k gen) bin is also small-mass but the bias drops
  from -0.174 (n=600) to +0.044 (n=1000).

**Conclusion**: per-tree-mean-age binning *is* a viable
ARG-only WHEN-axis primitive — but it needs n≥1000 to escape
seed-to-seed noise at L=100kb. For the manuscript Fig 5 (Block 2.5
WHEN axis), we should run at n≥1000 and likely 5+ seeds for the
ribbons to be tight. At n=1000 the estimator total `ĥ²` was close to
the simulated 0.5 (visible in `est.summary.json` if needed for the
draft figure).

### Files

- `test/eval/age_decomp_nscan.sh` — driver.
- `test/eval/arch_metrics.py::per_tree_age_binned_h2()` — the
  estimator-side aggregation; matches the truth-side computation by
  using the **same** bin edges for both truth and estimator.
- Output: `~/test-sargent/block2_1_v2_nscan/n{300,600,1000}_s{1,2,3}/`,
  `~/test-sargent/block2_1_v2_nscan/run.log` (aggregated tables).

### Reproduce

```bash
BASE=$HOME/test-sargent/block2_1_v2_nscan \
  bash test/eval/age_decomp_nscan.sh
```

(Resumable: skips cells with `est.summary.json` already present.)

## Caveats

- For variants NOT in the BED panel, `Var(x_j)` is approximated as 1
  (standardised). This is correct for causal variants that are in the BED;
  for off-BED causal variants the approximation can underweight (in
  practice the realised h² in the partition already accounts for this
  because `simulate.py`'s `scale_f` uses standardised dosages).
- The age field assumes msprime's mutation model emits `time` (current
  default); if unavailable, falls back to the parent-node time in the
  local tree at the variant position.
- TAU's effect interacts with the BayesR mixture: the variance prior
  classes (π_k, σ²_k) are sampled FIRST, then scaled by age. To recreate
  the "varies-with-age" effect-size *distribution* (vs *multiplicative
  scaling*), a deeper refactor would draw σ²_k itself age-dependently.
  v1 uses the multiplicative scaling for simplicity; rebound for v2 if
  needed.
