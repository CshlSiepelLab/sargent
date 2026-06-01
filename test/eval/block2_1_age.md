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
