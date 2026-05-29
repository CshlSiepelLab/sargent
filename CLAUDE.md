# SARGENT — Claude Code project guide

This file is auto-loaded by Claude Code in this repo. It is the entry point
for any session, including server-side sessions launched via SSH.

## What this project is

**SARGENT** (Sparse ARG-based ENcoding of Traits): an ARG-based polygenic
heritability inference tool. C core (production estimator), Python eval
pipeline (sim + scoring). PI: Adam Siepel (CSHL).

The current methodological focus and remaining open questions live in:

- **`ROADMAP.md`** — active plan (Blocks 0–5) at the top; preserved
  historical phases below. **Read this first** to understand current priorities.
- **`CLAIMS.md`** — manuscript-claim triage: HELD / PLAUSIBLE-UNTESTED /
  FALSIFIED-WEAK / OUT-OF-SCOPE, with `[MINOR]` markers on the
  trim-quietly group.
- **`test/eval/phase10_oracle.md`** — Phase 10 oracle diagnostic
  (per-segment ws² ceiling). The result that defined the current story
  (verdict: ws² is a localization tool, not a prediction tool).
- **`test/eval/block1_1_design.md`** — Block 1.1 (cross-ancestry PRS
  portability) design + SSH-launch instructions. **This is the active
  experiment.**
- **`test/eval/fit_blup_audit.md`** *(may not exist yet)* — Block 1.3
  audit of the +0.25 h² bias in `predict.py:fit_blup`. Was being run by
  an agent in the macOS-side session — check whether the file exists
  before assuming the audit is pending.

## Current state (as of 2026-05-29, late)

Most recent merge: `f0ee0b9` — `phase10_oracle.md` h² addendum refreshed
with patched (centered) numbers.

Active work blocks:

- **Block 1.1** (prediction headline): cross-ancestry PRS portability.
  Pipeline complete (`portability.sh`); v1 launched on Adam's Linux
  server (evolgen). 120 cells in flight. Plotting helper at
  `test/eval/plot_portability.py` writes the headline curve + per-arm bars
  when results land.
- **Block 1.2** (localization robustness sweep): scaffolding complete
  (`localization_sweep.sh`); v1 may be running on the laptop in parallel
  (30 cells, ~2 h wall) — check `~/test-sargent/block1_2_v1/run.log`.
  Design doc: `test/eval/block1_2_design.md`.
- **Block 1.3** (`fit_blup` h² bias): **audit done + patch applied**.
  `test/eval/predict.py` now double-centers K (see `center_kinship`);
  fit_blup h² bias dropped from +0.26 → +0.05 (audit's predicted floor)
  on all 10 Phase 10 sims. R² unchanged. Audit: `test/eval/fit_blup_audit.md`.
  Patch commit: `ebf757c`.
- **Block 0** (claims triage): Adam signed off (2026-05-29). `CLAIMS.md`
  is current; bucket reassignments only if further evidence comes in.

Deferred: Block 2 (sparse fused lasso, gated on Block 1.2 outcome),
Block 3 (manuscript, Adam-led, next week), Block 4 (efficiency),
Block 5 (publication readiness).

## Build / environment

```bash
# Build
mamba activate sargent           # provides cmake, meson, ninja, msprime, tskit
cmake -B build -S .
cmake --build build -j
ls build/src/progs/sargent       # the binary

# Python deps check
python -c "import msprime, tskit, numpy; print('ok')"
```

Python in scripts/eval assumes the `sargent` mamba env. The driver scripts
auto-detect `mamba run -n sargent python` if available, else fall back to
plain `python`. To force a specific interpreter, export `PYBIN=/path/to/python`.

## Common tasks

### Launch Block 1.1 v1 sweep

See `test/eval/block1_1_design.md` § "Launching v1 on Linux server /
workstation". Short form:

```bash
# Pilot first (30-45 min on a workstation; sanity-check the setup)
BASE=$HOME/test-sargent/block1_1_pilot \
  N_TRAIN=200 N_TEST=80 LENGTH=300000 \
  MIX_FRAC_B="0.00 0.05 0.20" TYPED_FRACS="0.7" SEEDS="1 2" \
  bash test/eval/portability.sh

# Then v1 (120 cells, 5-15h wall) in tmux
BASE=$HOME/test-sargent/block1_1_v1
mkdir -p $BASE
tmux new -s b11 -d "BASE=$BASE bash test/eval/portability.sh 2>&1 | tee $BASE/run.log"
```

Driver is resumable: skips any cell whose `est.predict.json` already exists.

### Monitor progress

```bash
tail -f $BASE/run.log
ls $BASE/mix*_typed*_s*/est.predict.json | wc -l   # completed cells / 120
```

### Interpret results (when v1 finishes)

`$BASE/portability_table.md` has three sections:

1. Per-seed table (every cell).
2. Summary by (mix × typed): mean ± sd over seeds, plus the paired
   `SARGENT − ARG-BLUP-uni` gap.
3. Headline curve: SARGENT − ARG-BLUP-uni gap as a function of mix_B.

Verdict rule (per `block1_1_design.md`):

- **Positive verdict**: SARGENT > ARG-BLUP-uniform *and* the gap widens with
  f_B (target-pop training fraction). Triggers v2 expansion (F_ST ladder,
  admixed scenario, ρ_β<1, eGRM baseline).
- **Null verdict**: no SARGENT gap at any ladder level. Prediction story
  collapses to localization-only; freeing cycles for Block 2 redirect.

When v1 reports, write findings into a new section at the bottom of
`block1_1_design.md` (mirror `phase10_oracle.md` § "Interpretation") —
do not start manuscript rewrites yet (that's Adam-led, Block 3, next week).

### Launching Block 1.3 from scratch (if the audit was not done previously)

If `test/eval/fit_blup_audit.md` does not exist, the audit is still pending.
The brief: locate `test/eval/predict.py:fit_blup`, diagnose the +0.25
systematic upward bias in the fitted h² (observed in Phase 10 across all
three BLUP arms; SARGENT's own `--arg-only` h² has bias only +0.06–0.10).
Candidates: (a) ML vs proper REML projection, (b) rank-1 dominance of ARG
kinship, (c) kinship normalization, (d) EM convergence. **Do not modify
production code** without Adam's approval — produce a diagnosis and recommended
fix, audit only.

## Guardrails

Adam is the PI; you are doing his bidding. Hold to these:

- **Never push to `main` or open a PR without Adam's explicit "go".** Commits
  are okay if asked; PR + merge requires "PR + merge" or similar.
- **Long simulation runs (Block 1.1 v1, Block 1.2, etc.) require explicit
  authorization.** Pilots and smoke tests are fine to launch.
- **Don't modify production C code (`src/`) for diagnostic purposes.** Write
  diagnosis scripts in `test/eval/` or `/tmp/` instead.
- **Don't claim "tests pass" or "results look good" without showing the
  numbers.** Adam values quantitative evidence over narrative.
- When in doubt, propose and ask before acting.

## Where to write things

- **Findings / analysis from a finished experiment** → append to that
  experiment's `*.md` (e.g., new section in `block1_1_design.md` or
  `phase10_oracle.md`).
- **New claims to triage** → propose an edit to `CLAIMS.md`.
- **Strategic decisions / new blocks** → propose an edit to `ROADMAP.md`.
- **Scratch / one-off diagnostics** → `/tmp/` (not in the repo).
- **Reusable diagnostic infrastructure** → `test/eval/`.

## Conventions

- Per-cell output convention: `<BASE>/<cell-name>/{case.*, est.*, *.log,
  est.predict.json}` (one directory per parameter cell).
- Per-cell parameter encoding in directory name, e.g.
  `mix0.05_typed0.70_s3` — the aggregator regex in `portability.sh`
  expects this.
- `simulate.py` outputs follow the prefix convention `<prefix>.{phen,
  trees, full.trees, bed, bim, fam, params.json, true_beta.tsv,
  true_ws2.tsv, train.anc.tsv, test.anc.tsv}`.
- Eval `*.md` docs use `Reproduce` blocks with exact env-var invocations
  at the end (see `phase10_oracle.md` for the template).

## If something feels off

If results are surprising, parameters seem wrong, or you suspect a bug
masking a result: stop, write up what you see in a short message, and
ask before propagating downstream. Adam has been burned by overconfident
intermediate steps before; an extra clarifying question costs less than
a wrong-direction sweep.
