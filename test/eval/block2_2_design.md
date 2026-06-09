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

---

## Phase C — Direct stratified GREML on per-bin ARG kinship (2026-06-08 pivot)

> **2026-06-09 status:** The engine works (A3a feasibility PASS, commit
> `1afe903`; A3a writeup at `block2_2_phaseC_A3a_feasibility.md`).  The
> **identifiability audit** that followed reframes Phase C: **Acts I-II
> (functional-annotation calibration against SNP-REML) are DEAD** by the
> same positional-collinearity wall Phase B hit (per-annotation GRM corr
> 0.73–0.86; multi-component REML collapses toward equal allocation).
> Only the **evolutionary axes** are identifiable — A3a (age) and A3b
> (frequency) are the spine; A3c (population) is at-risk pending a
> balanced-cohort test.  ROADMAP entry `6f8cf9d` and REFRAME_NOTES v6.7
> have the full audit.  Acts I-II section below is **preserved for the
> record** but not on the implementation path.

**Method change.** Phase B failed because it relied on **ŵ²-then-regress**:
SARGENT estimated a per-segment variance ŵ²_s, and per-annotation h² was
recovered via per-base regression of ŵ²_s on annotation coverage. The
mechanism that broke it (`corr(ŵ², true per-tree var) ≈ 0`) was not a coding
bug — it is a structural property of per-tree variance estimation under
adjacent-tree collinearity, and the ŵ²-ℓ2 regulariser that fixes the
ill-conditioning destroys the localization that the regression needs.

**Phase C bypasses ŵ² entirely.** Instead of estimating per-segment
variances and aggregating them, we partition **edges** by annotation
directly, build one ARG-kinship matrix per bin, and fit a **multi-component
REML** for the per-bin variance components. Standard machinery (Yang et al.
2011 GREML-MS, Speed et al. 2017 LDAK-Cat, GCTA `--reml --mgrm`), exact ML
under the model, no per-tree intermediate.

### Why this avoids the Phase B failure

Phase B asked the data to do something it cannot: localize each tree's
contribution to phenotype variance. With adjacent trees sharing >99% of
genealogy, the per-tree variances are jointly ill-conditioned — the
estimator can identify *sums* across nearby trees but not the trees
individually.

Phase C asks a strictly easier question: identify variance components
across **partitions of edges** that are aggregated *before* estimation.
K_b = Σ_{e ∈ bin_b} (ℓ_e/2)·m_e·m_eᵀ is exact, deterministic, and well-
conditioned as long as the bins are chunky enough that K_b's eigenstructure
is well-separated from the other K_b'. The estimator is the standard
multi-component REML — the same one used by every published partitioned-
heritability method. It does not need ŵ² to localize; it does not need
per-tree resolution at all.

| | Phase B (regression on ŵ²) | Phase C (direct stratified REML) |
|---|---|---|
| Per-segment variance estimate needed? | **Yes** — ŵ²_s must localize | **No** — bypass entirely |
| Per-bin kinship K_b | Implicit (via regression) | **Explicit** (sum over edges) |
| Estimator | OLS regression of ŵ² on annotation coverage | Multi-component REML on (K_1, …, K_K) |
| Identifiability | Collinear (Phase B failure) | Well-posed for chunky bins |
| Literature precedent | None | Yang/Speed multi-component GREML |

### Validation strategy: simulation first, real data deferred

The new method needs careful in-silico validation before any real-data
application. We validate against **SNP-stratified REML** (not S-LDSC):
same statistical machinery (REML), same data requirements (individual-
level), same sample-size regime (n ≥ 200 works). The only thing differing
between SNP-REML and ARG-REML is how K_b is constructed — from typed
genotypes vs from ARG edges. Any gap is attributable to feature space, not
to estimator differences.

S-LDSC is parked as a future comparison once we deploy at biobank n; in our
n ≤ 1000 sim regime its asymptotic regression is too noisy to be a fair
calibration partner.

### Two-act experimental arc

> **2026-06-09 update:** Acts I-II (functional-annotation calibration
> against SNP-REML) are **DEAD** — the per-annotation GRMs are mutually
> collinear (corr 0.73–0.86; multi-component REML collapses toward equal
> allocation).  Same positional wall Phase B hit.  Sections below preserved
> for the record; **only Act III is on the path.**  Act III itself
> re-prioritized: A3a (age) + A3b (frequency) are the spine, A3c
> (population) is at-risk pending balanced-cohort test.

**Act I (Acts I-II) — Validation against SNP-stratified REML** [DEAD,
2026-06-09] (calibration + gap-opening + bias-correction story; site-level
annotations).

**Act III — ARG-distinctive evolutionary annotations** (age, frequency,
population; only-ARG-can-do).  **This is now the entire path.**

Acts I and II share data (one annotated sim + one estimator pipeline, three
typing regimes). Act III is three separate sub-experiments.

### Acts I-II: three regimes against SNP-stratified REML [DEAD — preserved for record]

All three regimes use one simulation framework (`simulate.py
--annotation-genes`, already in place from Phase A) with three different
**typed-SNP panel ascertainment** policies. Causal variants and per-
annotation h² truth are fixed across regimes; only the typed/untyped split
moves.

#### Regime 1 — Full typing (concordance)

- **Typed panel**: every causal variant present.
- **Prediction**: SNP-REML and ARG-REML give the same σ²_b per annotation.
- **Pass criterion**: per-bin |σ²_b,SNP − σ²_b,ARG| < 0.05 (averaged over
  seeds), for every annotation.
- **Message**: under the conditions where SNP-REML is well-specified,
  ARG-REML recovers the same answer. *No statistical loss from changing
  the feature space.*

#### Regime 2 — Partial typing (gap-opening)

- **Typed panel**: random 50% of variants typed (no AF correlation —
  isolates the feature-coverage axis from the AF-bias axis).
- **Prediction**: SNP-REML under-estimates total h²; per-annotation
  estimates compressed toward zero. ARG-REML recovers full truth (uses
  all variation).
- **Pass criterion**: ARG-REML closer to truth than SNP-REML on all
  annotations; mean gap on the rare-enriched annotation ≥ 0.10.
- **Message**: as soon as typed coverage is incomplete, ARG-REML captures
  heritability SNP-REML cannot see.

#### Regime 3 — MAF-stratified typing (bias-correction)

- **Typed panel**: variants typed iff MAF ≥ 5% (array-like ascertainment).
- **Prediction**: SNP-REML *mis-allocates* — over-estimates annotations
  with low rare-variant load (intergenic) by absorbing common-only
  signal, under-estimates annotations enriched for rare variants
  (exon/LoF-like). ARG-REML unbiased on both.
- **Pass criterion**: sign of bias error opposite between enriched and
  unenriched annotations in SNP-REML; ARG-REML unbiased within sim
  noise (|bias| < 0.05).
- **Message**: SNP-REML partitioned-h² estimates are not just under-
  powered for rare-enriched annotations — they're systematically *mis-
  allocated*. ARG-REML fixes the bias, not just the variance.

#### Act I/II simulation parameters (provisional)

- n = 400, L = 5 Mb (matches Phase A reference sim; tractable REML).
- Annotation model: gene-structured (`--annotation-genes`); exon
  τ_e = 10, promoter τ_p ≈ 3.3, intergenic baseline (existing Phase A
  defaults).
- h² split: total h² = 0.5; per-annotation h² truth from
  `params.json:h2_by_annotation`.
- Seeds per regime: 5 (enough for SE on h² estimates with σ ≈ 0.03).
- 3 regimes × 5 seeds = **15 cells per Acts I-II**.

#### Site-to-edge annotation mapping

For site-level (genomic) annotations like exon/promoter, each ARG edge
carries 0 or more mutations, each at a specific genomic position with a
specific annotation. Mapping choice: **per-mutation rank-1
decomposition** — each mutation contributes one rank-1 kinship piece
(ℓ_e/2)·m_e·m_eᵀ to the bin of its annotation; edges with mutations in
multiple annotations contribute proportionally to each bin's K_b.

This is the cleanest mapping (each *mutation* sits in exactly one
annotation, just as S-LDSC and SNP-REML are per-mutation). It also
reduces to SNP-REML in the limit where every mutation is typed: the same
K_b assembled from the same mutations gives the same kinship matrix.

(Rejected: edge-majority assignment, which is ambiguous for edges with
mutations in multiple annotations; edge-weight = annotation-fraction,
which double-counts via shared edges.)

### Act III: ARG-distinctive evolutionary annotations

Three sub-experiments, each on its own simulation framework. Each tests
whether ARG-REML recovers a per-bin h² diagonal under a per-bin
architecture — the same recovery test as Fig 5 Panel B Phase B, but with
the stratified-REML aggregation instead of the per-edge posterior.

#### A3a — Allele age (edge midpoint age)

- **Simulation**: per-edge architectures with per-cell tertile age bins —
  **reuse the existing Fig 5 Panel B sims** (9 cells, n=300, L=200kb,
  3 archs × 3 seeds). No new simulations needed.
- **Estimator**: ARG-REML with K_recent / K_middle / K_ancient built from
  edges in each age tertile.
- **Pass criterion**: per-bin σ̂²_b matches simulated h² diagonal; per-cell
  total h² ≥ 0.5 (i.e., method recovers the bulk of h² even when
  architecture is per-edge — the failure mode of the per-edge posterior
  in Phase B was missing total h², not just mis-allocation).
- **Headline plot**: Fig 5 Panel B redo with a recovered diagonal.

#### A3b — Allele frequency (edge DAF)

- **Simulation**: new sweep, per-AF-bin architectures. Causal effects
  concentrated in edges with DAF in [0, 5%], [5%, 25%], [25%, 50%]
  respectively, three architectures.
- Suggested parameters: n=400, L=300 kb, 3 archs × 3 seeds = 9 cells.
- **Estimator**: ARG-REML with K_rare / K_low / K_common built from
  edges by descendant count.
- **Pass criterion**: per-bin σ̂²_b matches diagonal; **bonus**: head-
  to-head against SNP-REML on AF-binned typed genotypes — ARG-REML
  should outperform on the rare bin (rare variants partly untyped).
- **Headline plot**: per-AF-bin σ̂²_b recovery, ARG vs SNP.

#### A3c — Population identity (per-edge population label)

- **Simulation**: multi-population sim, e.g., 2-pop split or OOA-style.
  Reuse Block 1.1 v3 cells, or new stdpopsim cells if Block 1.1's are not
  per-edge-annotated.
- Suggested parameters: n=400 admixed, OOA or 2-pop with admixture pulse;
  edges labeled by descendant population. Causal effects in three bins:
  pre-split shared, post-split pop-A, post-split pop-B.
- **Estimator**: ARG-REML with K_shared / K_popA / K_popB built from
  edges by population assignment.
- **Pass criterion**: per-bin σ̂²_b matches diagonal; per-population h²
  recovery aligns with cross-ancestry portability findings from Block
  1.1 (shared variance portable; pop-specific variance not).
- **Headline plot**: per-population σ̂²_b recovery; ties back to Fig 6
  (cross-ancestry portability) with a mechanistic decomposition.

### Engineering checklist

Code that exists and we can reuse:
- `simulate.py --annotation-genes` (Phase A, verified).
- Per-edge streaming (used by `edge_h2.py`).
- Per-tree REML EM (`predict.py`).
- Phase B's `partition_lib.py` — site→edge mapping primitives (the
  mapping itself is reusable; the regression-on-ŵ² wrapper is not).

New code to write:
- **Multi-component REML** (in-house, or wrap GCTA). In-house first for
  control; GCTA wrapper added for review-defensibility.
  - Inputs: list of K_b matrices (as dense or as edge-stream specs), y, X.
  - Outputs: σ²_ε, σ²_1...σ²_K with REML log-lik and per-component SE.
  - Use AI-REML (average information) for efficiency at small K.
- **K_b assembly from edge stream + annotation lookup**. For site
  annotations: per-mutation iteration with site→annotation lookup. For
  evolutionary annotations: per-edge property (age, descendant count,
  population label) computed during the same edge pass.
- **SNP-REML head-to-head wrapper**. Build K_b,SNP from typed
  genotypes restricted to annotation b; same multi-component REML
  engine.
- **Aggregator** (`aggregate_b22c.py`) — tidy CSV across regimes/seeds
  with columns: `regime, annotation, true_h2, est_h2_ARG, est_h2_SNP,
  estimator_se, seed`.

Estimated implementation work: ~3–5 days of coding + validation.

### Pass/fail criteria summary

Acts I-II all pass → Phase C validated against SNP-REML on site
annotations. Greenlight Act III.

Act III all pass → manuscript-ready partitioned-h² framework with
ARG-distinctive evolutionary annotations.

Any fail mode triggers a write-up + diagnostic before the next regime.

### Open risks

- **Multi-component REML at n=400, K=3 may be under-powered.** σ²_b
  estimates can have SE ~0.05–0.10 per bin. Mitigation: 5+ seeds per
  cell, and report mean ± SE not point estimates. If still noisy, scale
  to n=1000.
- **Site-to-edge mapping has subtleties.** Mutations on the same edge in
  different annotations split the edge's kinship contribution across
  bins — fine in principle, but the resulting K_b are sums of *fractional*
  rank-1 pieces, which is numerically unusual. Validate by sanity-
  checking that Σ_b K_b = K_ARG (total ARG kinship) within numerical
  tolerance.
- **Population-edge assignment is tree-by-tree.** An edge may have
  descendants in multiple populations (admixture, incomplete sorting).
  Mapping choice: weight edge contribution to K_b by fraction of
  descendants in population b. Validate that this aggregates sensibly.
- **REML convergence with shrinkage near zero.** When a bin's true σ²_b
  is small, REML can hit the σ² ≥ 0 boundary and report 0. Standard
  fix: report unconstrained estimates alongside constrained for
  diagnostic clarity.

### Suggested implementation order

1. **Prototype A3a (age) first** — reuses existing Fig 5 Panel B sims, no
   new simulation needed. Smallest engineering ask. Tests the
   multi-component REML on a known-architecture per-edge sim. Sanity
   check: does ARG-REML recover the diagonal that the per-edge posterior
   couldn't? If yes, framework validated; if no, the new method shares
   Phase B's failure mode and we need to diagnose.
2. **Acts I-II validation against SNP-REML** — builds the calibration
   arc with the multi-component REML now known to work.
3. **A3b (frequency)** — new simulations, modest engineering.
4. **A3c (population)** — adapts multi-population sims; bridges to
   Block 1.1 / Fig 6.

Total Phase C effort: ~2–3 weeks of implementation + sweeps. Lighter than
Phase B's regression apparatus, because the estimator is standard
multi-component REML rather than a custom regression-on-ŵ² pipeline.

### Real data: explicitly deferred

Once Phase C passes in-silico validation across Acts I-II and III, a
real-data deployment (UK Biobank, AllOfUs) becomes the natural next paper.
Defer until simulation validation is complete. Engineering for biobank
scale (stochastic trace estimators, CG solves) is well-trodden territory
but not on the methods-paper critical path.

### Downstream follow-on: revisit Fig 5 Panel C

> **2026-06-09 update:** A3a landed PASS (commit `1afe903`).  The
> identifiability audit eliminates one of the two candidates below:
> **spatial-bin variance-capture (a) is dead** for the same positional-
> collinearity reason that killed Acts I-II — coarse spatial bins are
> still positional, and a multi-component REML on coarsely-binned
> contiguous-tree GRMs will hit the same near-identical-kinship wall.
> Only candidate **(b) frequency-tertile** survives, and it largely
> *collapses into Panel B (A3b)* under the v6.7 figure restructure
> (REFRAME_NOTES v6.7).  Panel C is therefore likely to be either
> retired or repurposed as a per-stratum *graded-architecture* recovery
> demonstration, not a standalone variance-capture panel.  Defer the
> final call to post-A3b graded-recovery report.

The current Fig 5 Panel C (variance-capture curves from Block 2.5.2)
has an honest but weak story: truth-side separates the three bf regimes
(0.005 / 0.05 / 0.5) by ~100× in M_e, but estimator-side ŵ² is flat at
M_e ≈ 12–20 across all regimes — the panel as drawn shows the gap
without recovering it.  Two candidate v2 directions were proposed
post-A3a; the identifiability audit settles them:

- **(a) [DEAD 2026-06-09] Spatial-bin variance-capture under GREML.**
  Partition the genome into ~K coarse spatial bins, build K_b per bin,
  fit multi-component REML, sort σ̂²_b descending, cumulative-sum →
  variance-capture curve.  The per-bin K_b are still *positional* and
  will be mutually collinear by the same argument that killed Acts I-II
  (relatedness is global; every genomic subset estimates the same
  genome-average relationship).  Coarser bins won't help; collinearity
  worsens at scale.

- **(b) Frequency-tertile variance-capture using A3b.**  Reframe the
  HOW axis from "spatial polygenicity" to "rare-vs-common variance
  concentration" — variance-capture curves binned by edge DAF instead
  of by spatial position.  **Surviving candidate**, but under the v6.7
  figure restructure (Panel B now = frequency-stratified GREML), this
  candidate largely *is* Panel B.  A separate Panel C only makes sense
  if it adds a distinct view (e.g., cumulative σ²_b within each
  frequency stratum, showing concentration vs spread).

**Likely outcome:** Panel C is retired or repurposed.  Current panel
stays as fallback only if no clean Panel C v2 emerges.  Decision:
post-A3b graded-recovery report.

### Reproduce (provisional — implementation pending)

```bash
# A3a (age) prototype — re-uses existing Fig 5 Panel B cells
PYBIN=/Users/acs/mambaforge/envs/sargent/bin/python
BASE=$HOME/test-sargent/fig5_panelB_edges
$PYBIN test/eval/stratified_reml.py \
  --base $BASE --annotation age-tertile \
  --out test/eval/figs_r/_cache/block2_2_phaseC_A3a.tsv

# Acts I-II — concordance + gap-opening + bias-correction
$PYBIN test/eval/stratified_reml_sweep.py \
  --base $HOME/test-sargent/block2_2_phaseC \
  --regimes full,partial,maf-strat \
  --n-seeds 5 \
  --out test/eval/figs_r/_cache/block2_2_phaseC_actsI-II.tsv
```

(These scripts do not exist yet — they're the engineering deliverables
of Phase C.)
