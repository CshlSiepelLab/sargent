# Block 1.1 — Cross-ancestry PRS portability (design v1)

Date: 2026-05-29.
Status: design locked (pending Adam's sign-off on simulate.py scaffolding).
Owners: Claude (sim prep + execution), Adam (decision review).
Backs claims: B4, B7, E1, E2, H6 (see `CLAIMS.md`).

## Goal

Test whether SARGENT (BayesR fixed + per-segment-weighted ARG random) outperforms
ARG-BLUP-uniform, GBLUP, and BayesR when training and test individuals come from
genetically diverged populations, as a function of how much target-population
data the training set contains.

This is the **new prediction headline**. If SARGENT does not win here, the
prediction half of the manuscript collapses to localization-only.

## Design

### Demography (msprime)

Two-population split, A and B, F_ST ≈ 0.15 (mid EUR-vs-EAS regime):

- Both populations: `N_e = 10⁴`.
- Split time: 3500 generations ago.
- Recombination rate: 1 cM/Mb.
- Mutation rate: 1.25 × 10⁻⁸.
- Sequence length: 1 Mb in v1.

Implemented as new `--demography transfer` in `scripts/simulate.py`. No
recent admixture pulse (that is the separate `admixed` setting).

### Training-set composition ladder

Train n = 500 individuals, mixture of A and B:

| level | f_B (train) | n_A | n_B (train) |
|---|---|---|---|
| pure-transfer | 0.00 | 500 | 0   |
| trace        | 0.05 | 475 | 25  |
| moderate     | 0.20 | 400 | 100 |

Test set: 200 individuals, **always 100 % B**.

(*This is the central experimental knob — it answers both "does genealogy alone
transfer?" (f_B=0) and "does SARGENT extract more value from a small in-target
sample than baselines?" (f_B>0) in one swept figure.*)

### Genetic architecture

- h²_total = 0.6.
- 1 % causal variants from the typed panel (≈ 100–200 causal SNPs at 1 Mb).
- BayesR-style mixture priors with σ² classes {0, 10⁻⁴, 10⁻³, 10⁻²}.
- **β shared perfectly across populations** (ρ_β = 1) in v1. Any portability
  gap is therefore attributable to LD / allele-freq / ARG divergence, not
  to true effect heterogeneity. ρ_β < 1 is a v2 expansion.
- No per-segment-`ws²` heterogeneity in v1 (the architecture for the
  *prediction* experiment is `diffuse`; the *localization* sweep — Block 1.2 —
  is where heterogeneity is varied).

### Typed coverage (two levels)

| level | f_typed | rationale |
|---|---|---|
| high | 0.7 | "good genotyping" — closer to imputed-array regime. |
| low  | 0.1 | LD-leaky regime, matches Phase 3's headline. |

### Arms (the comparison that wins or loses the prediction story)

| arm | tool | implementation |
|---|---|---|
| BayesR  | typed-SNP sparse-regression PRS | `sargent --no-arg` (or equivalent) |
| GBLUP   | typed-SNP GRM BLUP              | `predict.py:fit_blup` on G_typed |
| ARG-BLUP-uni | uniform ARG kinship BLUP    | `predict.py:fit_blup` on K_ARG (no weights) |
| SARGENT | BayesR fixed + per-segment-weighted ARG random | `sargent` full pipeline |

For SARGENT specifically, the post-fit decomposition records the fixed-effect
contribution and the random-effect contribution to test-set R² separately —
this is what backs H6 directly.

### Sample sizes

- n_train = 500 (total; ladder controls A/B split).
- n_test  = 200 (all B).
- 5 seeds per (ladder × coverage) cell.

### Metrics

- Within-A R²: train A→A on a held-out 200 A, as a sanity baseline.
- **Cross R²: train→B-test** — the headline. Reported per arm per cell.
- **Portability ratio**: R²(train→B-test) / R²(within-A train→A-test).
- For SARGENT only: fixed-only R², random-only R², full R² (additive
  decomposition).

## Sequencing

### v1 (this run, ~1-2 days wall, no prior decisions to back up)

3 mix levels × 2 coverages × 4 arms × 5 seeds = **120 fits**. Yes-or-no answer
on B4/B7/H6.

- **If positive** (SARGENT > ARG-BLUP-uni *and* the SARGENT gap widens with
  f_B): expand to v2.
- **If null** (no separation at any ladder level): prediction story collapses,
  manuscript becomes localization-only, freeing cycles for Block 2 redirect.

### v2 (expand if v1 positive)

- F_ST ladder (split time = 1500 / 3500 / 7000 → F_ST ≈ 0.05 / 0.15 / 0.25).
- Admixed scenario (recent pulse, α = 0.5, train on admixed cohort → predict on admixed).
- ρ_β < 1 (true effect heterogeneity).
- Add eGRM as an explicit baseline.

## Decision log (locked 2026-05-29)

| # | Question | Choice | Rationale |
|---|---|---|---|
| 1 | F_ST headline | 0.15 (split 3500 gen) | "Real" cross-ancestry; not the easy 0.25 case. |
| 2 | ρ_β | 1.0 in v1 | Cleaner mechanism; ρ<1 is a v2 follow-up. |
| 3 | n_train | 500 | 3× Phase 10; trades off vs seeds (5) and length (1 Mb). |
| 4 | Admixed | v2 only | Keeps v1 small; admixed adds design tax. |
| 5 | eGRM baseline | v2 only | Mechanistically ≈ ARG-BLUP-uniform; B2 says it won't change v1 verdict. |
| 6 | **Train composition** | **Ladder 0/5/20 % B** | Per Adam (2026-05-29): a single point conflates two claims; a ladder gives both the pure-transfer verdict and the deployment-regime verdict in one figure. |

## Reproduce (planned interface)

```bash
DEMOG=transfer SPLIT=3500 \
  N_TRAIN=500 N_TEST=200 LENGTH=1000000 H2=0.6 P_CAUSAL=0.01 \
  MIX_FRAC_B="0.00 0.05 0.20" \
  TYPED_FRACS="0.7 0.1" \
  SEEDS="1 2 3 4 5" \
  BASE=/path/to/block1_1 \
  bash test/eval/portability.sh
```

Outputs land per cell at `<BASE>/mix<f_B>_typed<f>_s<seed>/` with one
`<arm>.json` per arm and aggregate at `<BASE>/portability_table.md`.

## Open implementation tasks

- [x] `scripts/simulate.py`: `--demography transfer` + `--transfer-split-time` +
      `--mix-fraction-B`. Emits `.train.anc.tsv` / `.test.anc.tsv`. Smoke-tested.
- [x] `test/eval/portability.sh`: driver + aggregator. Resumable. Smoke-tested
      single-cell (42s at n=80, L=150kb) and 2-cell ladder.
- [x] No separate `portability.py` needed — `predict.py` already runs all 4 arms
      (BayesR / GBLUP / ARG-BLUP-uniform / SARGENT). Driver invokes it directly.
- [ ] Hold v1 sims behind explicit "go" from Adam — no execution before sign-off.

## Launching v1 on Linux server / workstation (SSH-friendly)

The full v1 sweep is **120 fits**, estimated 5–15 h wall on a 32-core box, or
several days on 8 cores. Designed to run detached so the SSH session can drop.

### One-time server setup (do once per machine)

```bash
ssh me@server
cd ~/path/to/sargent          # or `git clone https://github.com/asiepel/sargent.git`
git fetch origin && git checkout main && git pull --ff-only

# Build sargent (requires the `sargent` mamba env with meson/ninja)
mamba activate sargent
cmake -B build -S .
cmake --build build -j
ls -l build/src/progs/sargent  # expect an executable

# Verify python deps
python -c "import msprime, tskit, numpy; print('ok')"
```

### Launch v1 in a detached tmux session

```bash
# Pick where outputs go. Plan ~5-10 GB of disk for the full v1.
BASE=$HOME/test-sargent/block1_1_v1

tmux new -s b11 -d "
  cd ~/path/to/sargent && \
  mamba activate sargent && \
  BASE=$BASE bash test/eval/portability.sh 2>&1 | tee $BASE/run.log
"
```

(`tmux` keeps the run alive when SSH drops; reattach later with `tmux attach -t b11`.)

### Check progress at any time

```bash
ssh me@server
tail -f $HOME/test-sargent/block1_1_v1/run.log   # live progress
# or, count completed cells:
ls $HOME/test-sargent/block1_1_v1/mix*_typed*_s*/est.predict.json | wc -l
```

Progress lines look like `[37/120] mix=0.05 typed=0.70 seed=3 -> .../mix0.05_typed0.70_s3` followed by `done in 290s (avg 280s/cell, ETA 392m)`.

### Resume after interruption

The driver checks for `est.predict.json` per cell and skips any that exist —
just rerun the same `bash test/eval/portability.sh` command (re-attach to tmux
or run a fresh launch). Completed cells are not redone.

### Output you care about

- `$BASE/portability_table.md` — written at end. Three sections: per-seed,
  per-cell summary (mean ± sd), and the headline portability curve (SARGENT −
  ARG-BLUP-uniform gap vs mix fraction).
- `$BASE/mix*_typed*_s*/est.predict.json` — raw per-cell results.
- `$BASE/mix*_typed*_s*/{sim,sargent,predict}.log` — per-cell logs for
  debugging failures.

### Pilot before committing the full run

Recommended: run a 6-cell pilot first (~30-45 min) to confirm nothing's broken
on the server, then launch the full sweep:

```bash
BASE=$HOME/test-sargent/block1_1_pilot \
  N_TRAIN=200 N_TEST=80 LENGTH=300000 \
  MIX_FRAC_B="0.00 0.05 0.20" TYPED_FRACS="0.7" SEEDS="1 2" \
  bash test/eval/portability.sh
```

If the pilot produces a sensible `portability_table.md`, fire the full v1
with the default env vars (3 mix × 2 typed × 5 seeds = 120 cells).

## Interpretation (v1, 2026-05-31) — FINAL (30/30 cells)

Status: all 30 cells complete. (4 cells originally hung in Stage 2 — see the
CG-cap fix, PR #21; they were then lost once to an external OOM on the shared
box, and finally completed cleanly. The CG cap is validated accuracy-neutral:
cap=64 vs uncapped on a hard cell agree to Δh2_arg=1.2e-3, and a well-conditioned
cell reproduces uncapped to 5 sig figs.) **Caveat:** those same 4 cells expose a
*separate* fixed-effect fitting divergence (see "Fixed-effect divergence" below);
the headline verdict is reported on the 26 well-behaved cells and is unchanged
whether or not the 4 are included.

**Metric note.** Read these on `r2_calib` (squared Pearson correlation), not
raw R². In the cross-ancestry regime raw R² is destroyed by *calibration*,
not discrimination: a SNP rare in train (pop A) but common in test (pop B),
standardised by train mean/std, blows up the fixed-effect prediction variance
and sends raw R² to −3 even when the prediction is positively correlated with
the phenotype. `r2_calib` (= R² after the optimal linear recalibration)
isolates discrimination and is the fair cross-method metric. See the
`predict.py` / `portability.sh` patch (PR #21).

### Result

Mean `r2_calib` (squared correlation) by target-mix, on the **26 well-behaved
cells** (the 4 fixed-effect-divergent cells excluded; see below):

| mix_B | GBLUP | ARG-BLUP-uni | SARGENT | **SARGENT − ARG-BLUP-uni** |
|---|---|---|---|---|
| 0.00 (n=10) | 0.048 | **0.069** | 0.053 | **−0.016** |
| 0.05 (n=9)  | 0.038 | **0.096** | 0.047 | **−0.048** |
| 0.20 (n=7)  | 0.171 | **0.221** | 0.141 | **−0.080** |

The verdict is unchanged if the 4 divergent cells are kept (all-30 gap:
−0.016 / −0.057 / −0.117 — more negative at mix=0.20, same direction). So the
conclusion is robust to that choice.

Two facts: (1) prediction transfers and scales sensibly — `r2_calib` rises
with target-population training fraction (sanity ✓); (2) **SARGENT does not
beat ARG-BLUP-uniform — it is consistently slightly below, and the gap widens
(more negative) with f_B.** ARG-BLUP-uniform is the best arm at every level,
beating both the SNP-GRM (GBLUP) and SARGENT's per-segment-weighted variant.

The verdict rule required SARGENT **>** ARG-BLUP-uniform with the gap widening
**in SARGENT's favour**. We observe the opposite ⇒ **null (mildly negative)
prediction verdict.**

### Fixed-effect divergence (separate finding, flagged for follow-up)

Completing the 4 formerly-hung cells (`mix0.20_typed0.70` s1/s4/s5,
`mix0.05_typed0.10` s4 — all extreme cross-ancestry) surfaced a **fixed-vs-random
`h²` decomposition divergence in sargent's Stage-1 fit**: the *simulated* truth
is `h²_snp ≈ h²_arg ≈ 0.30`, but sargent fits **`h²_snp ≈ −3`, `h²_arg ≈ +3.8`**
(total stays ~sane). The garbage fixed effect collapses BayesR/SARGENT point
prediction to ~0 `r2_calib` on these cells, while ARG-BLUP-uniform (no fixed
effect) is unaffected — which is why the divergent cells make the mix=0.20 gap
look worse.

This is **pre-existing and independent of the CG-cap fix** (the uncapped run
reproduces `h²_snp = −3.13` to the digit; Stage 1 is untouched by PR #21) — the
cap merely let these cells finish and reveal it. It is the fixed-vs-random
non-identifiability (sargent's own note: "the classical fixed-vs-random
decomposition is not identifiable") blowing up in the extreme-divergence regime.
Tracked as a separate issue (#22); diagnosis + fix in
`test/eval/h2_decomp_divergence.md`. **Update (2026-06-01):** a trust-region
cap on the Stage-1 fixed-effect variance (`VAR_FIXED_VFCAP=0.6`) recovers these
cells' h²_snp from ≈−3 to sane near-truth values (harmlessly to good cells), but
end-to-end validation shows it does **not** recover their prediction (SARGENT
r²_calib stays ≈0; the bounded betas remain ancestry-confounded). So the cells
stay excluded from the prediction verdict; robust cross-ancestry prediction on
them awaits the ancestry-covariate fix.

### Why this is a robust null, not a wrong-regime artifact

The architecture actually simulated in v1 (driver passes neither `--regime`
nor `--ws2-architecture`, so both default): the ARG random effect is
`regime=sparse` — a spike-and-slab draw (95% null / 4% medium / 1% large ⇒
~5% of trees active), **scattered i.i.d.**, with `ws2_architecture=null` (no
contiguous large-effect blocks). NB this is a discrepancy with the design
text above, which described the prediction architecture as "diffuse": the
runs were sparse-*scattered*, not diffuse. Either way, v1 did **not** simulate
the *contiguous large-effect-block* landscape (`--ws2-architecture sparse`,
the case `build_ws2_architecture` calls "the regime that should favour
per-segment weights if anything does").

That favourable regime, however, was already tested in **Phase 10** with the
**oracle** (true) `ws²` (see `phase10_oracle.md`): 3 contiguous blocks, 5%
active, deep admixture (F_ST≈0.26). Even with *perfect knowledge* of the
weights, the prediction ceiling over uniform was only **+0.023 R²**
(95% CI +0.01…+0.04, n=5), and the *estimated*-weight arm fell **below**
uniform (−0.008). So:

- v1 here: estimated weights, scattered-sparse, **cross-ancestry** ⇒ SARGENT
  ≤ uniform, gap widening with f_B.
- Phase 10: oracle weights, contiguous-block, **admixed** ⇒ ceiling +0.023,
  estimated < uniform.

Both regimes agree. The prediction lever is **structurally small even with a
perfect estimator**, and the per-tree REML does not reach it. This confirms
the Phase 10 verdict cross-ancestry:

> **Per-segment `ws²` is a localization tool, not a prediction tool.**

### Consequence

The prediction half of the manuscript does not gain from per-segment weights;
ARG-BLUP-uniform is the prediction workhorse (with the right demography). The
per-segment-weights contribution earns its keep in **localization / partitioned
`h²`** (Phase 10: idealized AUC 0.84 vs estimator 0.74 — real headroom), which
is where the Block 2 fused-lasso should aim. Per the design's own decision rule,
this **redirects cycles from the prediction sweep to Block 2.**

### Caveats / what would still be genuinely new

- 26/30 cells; the capped re-run cell `mix0.20_typed0.10_s5` carries a Stage-1
  `h2_snp=−0.12` anomaly (fixed/random identifiability for `typed=0.10`),
  unrelated to the CG-cap fix; it does not drive the verdict (excluding it,
  the headline is unchanged: −0.016 / −0.048 / −0.071).
- Not excluded by either v1 or Phase 10: a *higher-contrast* contiguous-block
  architecture, or **ancestry-divergent effects** (ρ_β < 1, where some segments
  transfer better across ancestry than others). That is the one mechanism that
  could in principle give per-segment weighting a portability edge, and it is a
  v2 / Block 1.2 question — but Phase 10's tight CI makes a large swing unlikely.

## v2 results (F_ST ladder, 2026-06-02) — PROVISIONAL (framing Adam-led)

Run: 4 F_ST × 3 mix_B × 2 typed × 5 seeds = **120/120 cells, 0 failures**, at
n_train=500, n_test=200, L=1Mb (tmux `b11v2`, memory-aware parallel dispatcher).
Binary carries the CG-iteration cap and the vf-cap (`VAR_FIXED_VFCAP=0.6`).
Per-cell data: `test/eval/block1_1_v2_results.csv`; full table:
`<BASE>/portability_table.md` § "Portability surface by F_ST". Metric is
`r2_calib` (squared correlation) — raw R² is destroyed cross-ancestry (the
aggregator's raw-R² section shows a −220 blow-up at mix=0.05 from one
mis-calibrated cell).

### vf-cap held at scale: 0 / 120 cells with h²_snp < 0

The Stage-1 divergence is fully contained across the whole ladder, **including
F_ST=0.30**. So Block 1.1a (ancestry covariates) is **not** required for
h²-reporting robustness. (It remains relevant for SARGENT's *prediction* — see
below — but the −3 catastrophe is gone.)

### Portability surface — mean r²_calib by F_ST (30 cells/split)

| F_ST | BayesR | GBLUP | **ARG-uni** | SARGENT | ARG-uni−GBLUP | ARG-uni−BayesR |
|---|---|---|---|---|---|---|
| 0.07 | 0.112 | 0.112 | **0.155** | 0.140 | +0.043 | +0.043 |
| 0.12 | 0.098 | 0.119 | **0.172** | 0.124 | +0.053 | +0.074 |
| 0.16 | 0.054 | 0.091 | **0.130** | 0.066 | +0.039 | +0.076 |
| 0.30 | 0.037 | 0.062 | **0.089** | 0.048 | +0.027 | +0.052 |

- All arms decline with F_ST (expected). **ARG-BLUP-uniform is the best arm at
  every F_ST.**
- Retention (F_ST 0.07→0.30): BayesR 33%, SARGENT 34%, GBLUP 55%, ARG-uni 57%.
  **BayesR collapses hardest — and SARGENT collapses with it**, because it
  inherits the BayesR fixed effect (the cross-ancestry confounding the vf-cap
  bounds but does not fix). GBLUP ports about as well as ARG-uni.
- SARGENT < ARG-uni at every F_ST → **within-method null holds** (consistent
  with v1; SARGENT−ARG-uni = −0.024 / −0.026 / −0.076 by mix).

### The headline claim "framework−typed gap widens with F_ST" is NOT supported

Framework−typed (framework=max(ARG-uni,SARGENT), typed=max(BayesR,GBLUP),
per-cell paired): +0.036 / +0.043 / +0.031 / +0.025 at F_ST 0.07/0.12/0.16/0.30.
The ARG framework **beats typed-SNP methods at every F_ST** (~+0.03–0.05), but
the gap is **flat-to-slightly-narrowing** (peaks at mid-divergence F_ST≈0.12).
A figure built on the F_ST axis shows a roughly flat gap, not a widening one.

### Where the framework advantage DOES widen: typed coverage, not F_ST

ARG-uni is **coverage-invariant** (it uses the ARG kinship, not typed SNPs);
GBLUP degrades at low coverage:

| F_ST | ARG-uni − GBLUP @ typed=0.70 | ARG-uni − GBLUP @ typed=0.10 |
|---|---|---|
| 0.07 | +0.018 | +0.068 |
| 0.30 | +0.016 | +0.039 |

The gap is ~3–4× larger at low coverage, at every F_ST. Target injection helps
ARG methods most: at F_ST=0.30 mix=0.20, ARG-uni recovers to 0.182 while BayesR
stays at 0.066.

### Verdict (provisional)

- **Holds:** ARG framework > typed-SNP across the ladder; BayesR (sparse PRS)
  collapses worst; **ARG-BLUP-uniform is the robust prediction workhorse**; the
  advantage is largest in the **low-coverage** and **with-target-injection**
  regimes.
- **Does not hold:** the "framework advantage widens with F_ST" headline.
- **Manuscript implication:** lead the cross-ancestry figure with the
  **typed-coverage axis** (and BayesR's collapse), not F_ST; and use
  **ARG-BLUP-uniform**, not SARGENT, as the framework's prediction representative
  (SARGENT's BayesR fixed effect is the cross-ancestry liability). Framing is
  Adam-led — these numbers are captured here; the figure narrative is Block 3.

## SARGENT cross-ancestry rescue test (SARGENT-W, 2026-06-03)

**Question:** the v2 surface showed SARGENT (full) collapses to BayesR levels
cross-ancestry (it inherits the β fixed-effect liability). Can it be rescued by
dropping β? **SARGENT-W** = `sargent --arg-only` = per-tree ŵ² with **no β
fixed effect** (predict path = K_sg·α only). Re-fit on the existing 120 v2 cells
(no re-sim); **120/120, 0 failures**; ŵ²-only so the h²_snp divergence cannot
occur (W-arm h²_arg ∈ [0.45, 0.70], 0 negative). The EM variant (`--em 5`) was
scoped out (≥5× cost, anticipated negative).

NB: the SARGENT-W cells were fit with the **new EMMA Stage-2 binary** (perf PR),
the v2 baselines with the old Hutchinson binary. A 3-cell pilot confirmed they
agree within ~0.01 r²_calib (Δh²_arg ≤ 0.006, Δr²_calib ≤ 0.009; ws² correlates
0.92–0.98 on signal cells), so the comparison is valid.

### Result — r²_calib by F_ST (the three SARGENT flavours)

| F_ST | **SARGENT-U**(uniform ŵ², no β) | SARGENT-W (per-tree ŵ², no β) | SARGENT-full (β + per-tree ŵ²) | W−full | W−U |
|---|---|---|---|---|---|
| 0.07 | **0.155** | 0.143 | 0.140 | +0.003 | −0.012 |
| 0.12 | **0.172** | 0.145 | 0.124 | +0.021 | −0.027 |
| 0.16 | **0.130** | 0.102 | 0.066 | +0.036 | −0.028 |
| 0.30 | **0.089** | 0.070 | 0.048 | +0.022 | −0.019 |

(SARGENT-U = ARG-BLUP-uniform from the v2 surface.)

### Verdict: the rescue is PARTIAL — SARGENT-U remains the cross-pop representative

- **Dropping β helps (or is neutral) at every F_ST** — SARGENT-W ≥ SARGENT-full
  by +0.00 to +0.04, biggest at mid-divergence (F_ST=0.16). So β is confirmed
  the cross-ancestry liability the collapse came from.
- **But SARGENT-W never reaches SARGENT-U** — it stays −0.01 to −0.03 below
  uniform at every F_ST. The **per-tree ŵ² weighting is a residual drag below
  uniform**, independent of β (consistent with Phase 10 / Block 1.2 / v1:
  ŵ² is a localization primitive, not a prediction primitive).
- **SARGENT-U (uniform ŵ², no β) is the best SARGENT flavour at every F_ST.**

**Both ingredients are cross-ancestry liabilities for prediction**: β (sharply,
at mid-F_ST) and per-tree ŵ² (mildly, throughout). The clean manuscript framing
the data supports: **SARGENT-U is the cross-population framework representative;
SARGENT-with-β is the within-population representative.** Use ARG-BLUP-uniform
for the cross-ancestry figure, not SARGENT-full or SARGENT-W.

Per-cell data: `test/eval/block1_1_v2_variants_results.csv`
(`sargentW_r2`, `sargentW_r2c`, `sargentW_h2arg`). Driver:
`test/eval/sargent_variants_v2.sh` (VARIANTS=W). Provisional framing — Adam-led.

## Block 1.1 v3 — Realistic demographies (infrastructure done 2026-06-05)

**Status:** infrastructure complete and smoke-tested. Sweeps pending.

### Goal

Re-run the v1/v2 verdict under standard pop-gen-aware demographic models
(stdpopsim HomSap catalog: Gravel OOA, Browning AmericanAdmixture,
Ragsdale archaic).  Required for journal submission where reviewers
expect realistic simulations (AJHG, Genome Research, Nature Genetics).

### Infrastructure

`scripts/simulate.py` gains `--demography stdpopsim` with three
companion flags:

| flag | purpose |
|---|---|
| `--stdpopsim-model NAME` | stdpopsim HomSap model ID (e.g. `OutOfAfrica_3G09`) |
| `--train-pop NAME` | population sampled for training (e.g. `CEU`) |
| `--test-pop NAME` | population sampled for testing (e.g. `YRI`, `ADMIX`) |

`--mix-fraction-B` still controls the target-population injection fraction
(same as v1/v2 `transfer` demography).  Individual layout in the tree
sequence is identical to `transfer`:
`[0..n_train_A-1]` from `train_pop`, `[n_train_A..n_total-1]` from `test_pop`.
All downstream code (sargent, predict.py, portability.sh aggregator) is
unchanged.

`params.json` additions:
- `stdpopsim_model`, `train_pop`, `test_pop` — model / population provenance
- `fst_train_test` — Hudson F_ST between all train-pop and test-pop samples
  (computed via `ts.Fst()` from the full post-mutation tree sequence)

stdpopsim 0.3.0 installed in the `sargent` mamba env.  Models available:

| model ID | pops relevant to v3 |
|---|---|
| `OutOfAfrica_3G09` | YRI, CEU, CHB |
| `AmericanAdmixture_4B18` | AFR, EUR, ASIA, ADMIX |
| `OutOfAfricaArchaicAdmixture_5R19` | YRI, CEU, CHB (modern only) |

### Smoke-test results (2026-06-05, n_train=100, n_test=50, L=300 kb, seed=1)

4 cells, PASS=4, FAIL=0.  Driver: `test/eval/block1_1_v3_smoke.sh`.

| cell | model | train→test | F_ST | ARG-BLUP-uni r²_calib |
|---|---|---|---|---|
| OOA_3G09_CEU_YRI | OutOfAfrica_3G09 | CEU→YRI mix=0 | 0.0842 | 0.0007 |
| OOA_3G09_CEU_YRI_mix05 | OutOfAfrica_3G09 | CEU→YRI mix=0.05 | 0.0986 | 0.1832 |
| ADM_4B18_EUR_ADMIX | AmericanAdmixture_4B18 | EUR→ADMIX mix=0 | 0.0110 | 0.2285 |
| ARC_5R19_CEU_YRI | OutOfAfricaArchaicAdmixture_5R19 | CEU→YRI mix=0 | 0.1052 | 0.0175 |
| ARC_5R19_CEU_Nean1000 | OutOfAfricaArchaicAdmixture_5R19 | CEU→Nean@t=1000gen | 0.3919 | 0.0004 |

Notes:
- Pure CEU→YRI transfer (mix=0): F_ST≈0.08, prediction near zero — consistent
  with v2 pure-transfer behaviour.
- 5% YRI injection into training: r²_calib 0.001 → 0.183 at this tiny n,
  showing target-population injection mechanism works under realistic OOA demography.
- AmericanAdmixture EUR→ADMIX: low F_ST (0.011) as expected; r²_calib 0.23
  at n=100 (European ancestry dominates admixed test cohort).
- Ragsdale archaic model simulates identically to OOA for modern-pop samples.
- **Neanderthal probe**: F_ST≈0.39–0.54 (varies with L and seed), well above
  the modern-human ceiling (~0.15 in Gravel OOA). ARG-BLUP-uni r²_calib≈0
  as expected — prediction at extreme divergence. Used via
  `--test-pop Neanderthal --test-pop-time 1000`. Scenario 2 (archaic-ancestry
  mapping with non-zero `--train-pop-time`) uses the same flags.

### Planned sweeps (pending authorization)

See ROADMAP.md § "Block 1.1 v3" for the full grid:

| sub | demography | cells |
|---|---|---|
| v3.OOA | `OutOfAfrica_3G09`, 3 pop-pairs × 2 typed × 2 mix_B × 5 seeds | 60 |
| v3.ADM-tr | `AmericanAdmixture_4B18`, EUR train → ADMIX test | 20 |
| v3.ADM-te | `AmericanAdmixture_4B18`, test from ADMIX | 30 |
| v3.ARC | `OutOfAfricaArchaicAdmixture_5R19` (stretch) | 5–10 |

The existing `portability.sh` driver is compatible without modification.
For the stdpopsim cells, call:

```bash
BASE=$HOME/test-sargent/block1_1_v3_ooa
tmux new -s b11v3 -d "
  BASE=$BASE \
  DEMOG_ARGS='--demography stdpopsim --stdpopsim-model OutOfAfrica_3G09 --train-pop CEU --test-pop YRI' \
  bash test/eval/portability_stdpopsim.sh 2>&1 | tee $BASE/run.log
"
```

(A thin `portability_stdpopsim.sh` wrapper around `portability.sh` that
substitutes `--demography stdpopsim ...` for `--demography transfer ...`
is needed; write it at launch time.)
