# Block 2.2 — Partitioned heritability across functional annotations

Goal: the manuscript's **Fig 5 Panel A** (the "where does the heritability
live?" head-to-head). Per annotation, four numbers that (approximately) sum to
100%: **A (typed S-LDSC) | B (ARG-only) | C-typed | C-untyped**. The headline
story-cell is **B > A** / **C-untyped > 0** — a real rare-variant fraction the
typed-SNP method misses.

Plan: **Stage 1** = a v1 Fig 5 Panel A draft (this doc), **Stage 2** = full
validation across regimes + diffuse control + recovery-bias (later;
`partitioned_h2.md`).

> **STATUS (2026-06-04): Stage-1 premise FALSIFIED (honest negative, partial
> rescue left open).** The full three-flavour pipeline was built and runs
> end-to-end, but the science does not hold: SARGENT's per-tree ŵ² does **not
> localize heritability to the correct trees** (corr ≈ 0 with true per-tree
> variance, due to adjacent-tree collinearity), so per-annotation partitioning
> recovers at most a weak, heavily-compressed signal. Flavours A and C also
> break. See **"Verdict — Phase B"** at the bottom. One untested lever (high
> n / very few trees + low ŵ²-ℓ2) could give a partial rescue; not yet run.
> Manuscript implications deferred to a separate session (Adam).

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

## Phase B pilot (n=400, 1.5 Mb, big-signal) — findings

End-to-end pilot of all three flavours (exon 5% τ=10, promoter τ=3.3, equal
effects, common-only typed panel via --bed-maf-min 0.05, h2 split 0.35 typed /
0.25 rare-untyped). Truth clean: exon h2_total 0.148 / promoter 0.043 /
intergenic 0.354; ~75% untyped.

Result table (truth | B | A | C-typed | C-untyped=B−A):

| anno | tot | B(ARG) | typ | A(typed) | C-typed | C-unt(B−A) |
|---|---|---|---|---|---|---|
| exon | .148 | .063 | .034 | .140 | .000 | −.077 |
| promoter | .043 | .042 | .008 | −.030 | .000 | .072 |
| intergenic | .354 | .429 | .098 | .150 | .004 | .279 |

Findings:
1. **Flavour B attribution method matters decisively.** Span-FRACTION weighting
   of ws² returns ~base_frac·h2_arg (NO enrichment) — a 150 bp exon gets ~3% of
   its kb tree's ws² even when the exon variant elevated that tree. The
   per-base RATE regression (E[ws²_s]=Σ_c ρ_c·bp_{s,c}) detects enrichment:
   promoter recovered exactly (.042 vs .043), exon compressed (.063 vs .148).
   Diagnostic: exon-containing trees carry 1.77× the ws² of others (signal is
   real); the per-tree resolution smears a sub-tree annotation, so recovered
   enrichment ≈ 2.7× when truth is 7.8× (documented compression). **Use the
   regression.**
2. **Flavour A (typed HE proxy) over-captures clustered annotations.** Exon
   A-estimate (.140) ≈ exon TOTAL truth (.148), not typed-only (.034): typed
   common SNPs in a tight 150 bp exon are in strong local LD with the exon's
   rare variants and TAG them. So the array does NOT miss the rare fraction in
   clustered annotations → B−A goes negative for exon. (Promoter A is pure noise
   at n=400.) This is partly a real effect, partly n=400 power.
3. **C-typed ≈ 0**: SARGENT shrinks the typed β̂ to ~0 at n=400, so the hybrid
   typed slice is empty and the internal typed/untyped split is unusable — the
   Block-1.1-style typed-vs-ARG separation difficulty, as predicted.

Implication: the ARG primitive (B) is the viable leg; the typed leg (A) and the
B−A subtraction break on (a) n=400 noise and (b) the real fact that clustered
annotations are tagged by common SNPs. Open question for the manuscript framing:
partitioned h² may work cleanly only for tree-scale-or-coarser annotations
(gene bodies / regulatory domains), not individual 150 bp exons.

## Verdict — Phase B: honest negative (2026-06-04)

**Bottom line: fine-scale ARG-based partitioned heritability (the Fig 5 Panel A
headline as designed) is not supported with SARGENT's per-tree ŵ² primitive at
any scale we can afford. The root cause is that ŵ² does not localize
heritability to the correct trees — corr(ŵ², true per-tree variance) ≈ 0 —
because adjacent ARG trees are near-identical genealogies (collinear), making
per-tree variances ill-conditioned; the ŵ²-ℓ2 penalty then shrinks them toward
uniform. This is the Phase 10 conclusion ("ŵ² is a localization tool, not a
prediction tool") sharpened: for partitioning it barely localizes at all.**

### What was built (works; reusable; the negative is about the science not the plumbing)
- `simulate.py --annotation-genes` + `--bed-maf-min`: gene-structured
  annotations (single-exon genes + upstream promoters + intergenic) with
  per-base causal-DENSITY enrichment τ, and a common-only (array-like) typed
  panel. Ground truth `h2_by_annotation` (LD-free Σβ² partition, typed/untyped
  split) recorded in params.json. Verified correct (Phase A).
- `partition_lib.py` + `run_ldsc.py` (A, typed HE/GREML proxy with an S-LDSC
  stub for Stage 2) + `arg_ldsc.py` (B, per-base ŵ²-rate regression) +
  `hybrid_ldsc.py` (C) + `aggregate_b22.py` (→ Fig 5 CSV). Full pipeline runs
  end-to-end and is correct.

### Evidence chain (all configurations consistent)

| config | n | n_trees | n/n_tr | annotation | corr(ŵ²,true var) | genic/exon B recovery |
|---|---|---|---|---|---|---|
| pilot  | 400 | 4097 | 0.10 | 150 bp exon  | exon-trees 1.77× ŵ² | truth 7.8× → 2.7× (compressed) |
| coarse | 400 | 4097 | 0.10 | 5 kb genic   | 0.001  | truth 2.3× → 0.85× (depleted!) |
| short  | 400 | 273  | 1.47 | 1.5 kb genic | −0.015 | truth 4.1× → 1.47× (compressed) |

Raising n/n_trees 15× (0.10 → 1.47) did **not** restore the correlation ⇒
identifiability is not the binding constraint; collinearity is.

### Why each flavour fails
- **B (ARG-only).** Span-FRACTION weighting of ŵ² returns ≈ base_frac·ĥ²_arg
  (no enrichment) — a 150 bp exon gets ~3% of its kb tree's ŵ² even when the
  exon variant elevated that tree. The per-base RATE regression
  (E[ŵ²_s]=Σ_c ρ_c·bp_{s,c}) is the correct method and recovers some signal,
  but ŵ² itself does not track the true per-tree variance (corr ≈ 0), so the
  aggregate enrichment is weak and heavily compressed. **Root failure.**
- **A (typed HE/GREML proxy).** Unstable at n=400 (promoter h² estimated
  negative). And in tight clustered annotations (150 bp exons) the common typed
  SNPs are in strong local LD with the exon's rare variants and TAG them, so the
  array does NOT miss the rare fraction there — the B−A subtraction goes
  negative. A real effect (clustered annotations are well-tagged), not only n=400
  noise.
- **C (hybrid).** SARGENT shrinks the typed β̂ to ≈0 at this scale, so C-typed
  is empty and the internal typed/untyped split is unusable — the same
  typed-vs-ARG separation difficulty as the Block 1.1 h²_snp divergence.

### Mechanism
Per-tree variance estimation is over-determined in the count sense (~n²/2
phenotype pairs ≫ n_trees parameters) but **ill-conditioned**: neighbouring
trees share almost all of their genealogy, so their per-tree kinship matrices
are nearly collinear and the individual per-tree variances cannot be separated.
The ŵ²-ℓ2 regularisation resolves the ill-conditioning by shrinking toward a
flat profile, which destroys exactly the localization that partitioning needs.
More individuals raise SNR but leave the collinearity untouched — hence
corr(ŵ², true) ≈ 0 across a 15× change in n/n_trees.

### Open levers (possible PARTIAL rescue — untested, left for a future session)
1. **High n / very few trees + low ŵ²-ℓ2.** n≈1500, L≈100 kb (n/n_trees ≈ 6)
   with `--ws2-l2` → 0 (minimise shrinkage toward uniform). The single most
   favourable regime; one ~30–60 min fit would settle whether a narrow operating
   point exists. Read: collinearity makes a clean rescue unlikely, but this is
   the honest remaining test.
2. **Coarser / non-positional partition.** A 2-way genic/intergenic split at
   large scale, or the non-positional allele-age partition (Block 2.1's
   per-branch primitive), where the per-tree resolution demand is far lower than
   150 bp–5 kb positional annotations.

### Scope note
This bears only on the partitioning APPLICATION (Block 2.2 / Fig 5 Panel A WHERE
axis). The localization result (Block 1.2 / Phase 10) and the heritability /
prediction results are unaffected. Manuscript implications to be considered
separately (Adam, another session).

### Reproduce (coarse premise check + the localization diagnostic)
```bash
PYBIN=/DATA/home/asiepel/miniforge3/envs/sargent/bin/python
SG=build/src/progs/sargent; BASE=$HOME/test-sargent/block2_2_v1/coarse
$PYBIN scripts/simulate.py --annotation-genes --gene-mean-exon-bp 5000 \
  --gene-exon-frac 0.20 --gene-promoter-bp 2000 --gene-tau-exon 8 \
  --length 1.5e6 --n-indiv 400 --h2-snp 0.5 --h2-arg 0.0 \
  --sigma2-mix 0 1.0 --pi-mix 0 1 --causal-fraction 0.012 --seed 1 --out $BASE/case
$SG --silent --nprobe 10 --arg-only --ws2-l2 0.5 --outroot $BASE/est_B \
  $BASE/case.phen $BASE/case.trees $BASE/case
$PYBIN test/eval/arg_ldsc.py --case-prefix $BASE/case --trees $BASE/case.trees \
  --est-prefix $BASE/est_B --out $BASE/est
# localization diagnostic (corr ŵ² vs true per-tree variance): /tmp/b22_short.sh
```
