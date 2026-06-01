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
Tracked as a separate issue; diagnosis in `test/eval/h2_decomp_divergence.md`.

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
