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
