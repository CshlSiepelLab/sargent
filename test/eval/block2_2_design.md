# Block 2.2 — Partitioned heritability across functional annotations

Goal: the manuscript's **Fig 5 Panel A** (the "where does the heritability
live?" head-to-head). Per annotation, four numbers that (approximately) sum to
100%: **A (typed S-LDSC) | B (ARG-only) | C-typed | C-untyped**. The headline
story-cell is **B > A** / **C-untyped > 0** — a real rare-variant fraction the
typed-SNP method misses.

Plan: **Stage 1** = a v1 Fig 5 Panel A draft (this doc), **Stage 2** = full
validation across regimes + diffuse control + recovery-bias (later;
`partitioned_h2.md`).

## The three estimators (all yield per-annotation h²_c)

| flavour | data it sees | what it captures |
|---|---|---|
| **A** typed-SNP S-LDSC | common typed SNPs (GWAS sumstats + LD ref) | annotation h² tagged by common typed variants — **misses rare** |
| **B** ARG-only | the ARG + phenotype (SARGENT ŵ²_s) | full annotation h² incl. rare (trees encode every variant) |
| **C** hybrid | typed β̂ + ARG ŵ²_s jointly | a **(typed, untyped)** pair per annotation |

**Per-annotation untyped fraction = B − A** (subtraction). This is the *primary*
number: A and B each use one signal, so each is internally clean — no fragile
β-vs-ws² separation. SARGENT's *internal* joint typed/untyped split (flavour C)
is reported only as a **cross-check**; if it disagrees with B−A that is itself
the honest finding (the separation is unstable — same difficulty seen in the
Block 1.1 h²_snp divergence), and we fall back to the clean A-vs-B comparison.

## Two load-bearing design decisions

1. **Single causal architecture drives all three flavours.** One set of causal
   variants, positioned in annotations, generates the phenotype. Typed causal
   (on the array) → β̂ → A / C-typed. Untyped/rare causal (off-array) → only the
   ARG → B / C-untyped. So A, B, C all measure the *same* underlying truth.
   (Rejected: enriching the ARG ws² independently of the causal variants — that
   would make B measure something unrelated to A, and the head-to-head
   meaningless.)

2. **Enrichment via causal DENSITY, not effect size.** A variant's causal
   probability ∝ τ_{annotation}; effect sizes are constant. So enriched
   annotations carry τ-fold more causal variants per base — keeping the signal
   polygenic (which S-LDSC and the ws²-regression both need) while giving
   per-annotation h² ∝ base_frac × τ. (Rejected: scaling β² by τ — produces a
   few huge-effect spikes, near-monogenic, bad for both estimators.)

## Phase A — `simulate.py --annotation-genes` (DONE, verified)

A simple, interpretable gene model (not abstract tiles):

- **Exons**: single-exon "genes", length ~lognormal(mean `--gene-mean-exon-bp`,
  default 150 bp), placed with exponential intergenic gaps so exonic bases ≈
  `--gene-exon-frac` (default 0.015).
- **Promoters**: `--gene-promoter-bp` (default 100) immediately upstream of each
  exon → promoter bases ≈ exon_frac·promoter_bp/mean_exon_bp ≈ 1% automatically.
- **Intergenic**: the complement (~97.5%), baseline.
- **Enrichment**: exon `--gene-tau-exon` (default 10), promoter
  `= tau_exon · --gene-tau-promoter-ratio` (default 1/3 → 3.3), intergenic 1.

Outputs (only when `--annotation-genes` is set; fully backward-compatible
otherwise):

- `<out>.annotations.tsv` — per-variant `snp_idx, bp_pos, annotation, in_bed,
  is_causal`.
- `<out>.annot_map.tsv` — the contiguous full-genome partition `start_bp,
  end_bp, annotation` (for flavour A per-SNP annotation + flavour B per-tree
  mapping).
- `params.json` adds `annotation_genes`, the `gene_*` params, `annotation_ids`,
  and **`h2_by_annotation`**: per category the realized `base_frac, tau_c,
  n_causal, h2_typed, h2_untyped, h2_total` — the **ground truth** every
  estimator is scored against (LD-free additive Σβ² partition, the quantity
  S-LDSC targets).

### Phase A verification (smoke tests)

- Base fractions exact: exon 1.51% / promoter 0.99% / intergenic 97.51%.
- `annot_map.tsv` contiguous, no gaps/overlaps, covers full genome.
- Density enrichment ∝ τ: promoter 3.2–4.1× (target 3.3, exact); exon 6.6–8.6×
  (target 10) — undershoot at high causal-fraction is **saturation** (10× on
  1.5% of bases demands a large share of exon variants be causal); approaches 10
  at causal-fraction ≤ 0.015. Realized enrichment is recorded as truth either
  way.
- Typed/untyped split populates with `--bed-coverage < 1` and/or `--h2-unobs`.
- Backward-compatible: no flag → no annotation files, `annotation_genes=False`.

## Phase B (next) — reference sim + three estimators

- **B1**: one ~5 Mb gene-annotated reference sim; `sargent --ws2-l2 0.5` →
  `est.beta.tsv` + `est.w2.tsv`.
- **B2**: `test/eval/run_ldsc.py` (A; partitioned-GREML proxy + note if the
  Finucane S-LDSC binary isn't accessible), `test/eval/arg_ldsc.py` (B;
  regress ŵ²_s on per-tree annotation indicators, report the matrix condition
  number — degeneracy is the structural risk to watch), `test/eval/hybrid_ldsc.py`
  (C; joint split as cross-check vs B−A).
- **B3**: tidy `test/eval/block2_2_v1_results.csv` —
  `annotation, true_tau, est_h2c, est_se, estimator_flavour` (A|B|C-typed|C-untyped).
  Commit + push for Adam to render Fig 5 Panel A v1.

## Reproduce (Phase A smoke test)

```bash
PYBIN=/DATA/home/asiepel/miniforge3/envs/sargent/bin/python
$PYBIN scripts/simulate.py --annotation-genes --length 5e6 --n-indiv 400 \
  --h2-snp 0.3 --h2-arg 0.2 --h2-unobs 0.1 --n-unobs-causal 60 --unobs-maf-max 0.05 \
  --bed-coverage 0.6 --causal-fraction 0.02 --seed 1 --out /tmp/g
# inspect params.json -> h2_by_annotation
```
