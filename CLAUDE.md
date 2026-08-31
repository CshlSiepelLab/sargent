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
- **`test/eval/fit_blup_audit.md`** — Block 1.3 audit of the +0.25 h²
  bias in `predict.py:fit_blup`. Done (commits `53af237`, `ebf757c`);
  root cause: rank-1 constant offset in ARG kinship; fix: double-center K.

## Current state (as of 2026-06-03)

Branch context: `main`. All active PRs (#21, #23, #24) merged.

Active work blocks:

- **Block 1.1** (cross-ancestry PRS portability): v1 + v2 + SARGENT-W
  rescue test all done. **Verdict: ARG framework (SARGENT-U) beats
  typed-SNP methods across the F_ST ladder; per-segment ŵ² and β are
  both cross-ancestry liabilities vs ARG-BLUP-uniform.** The manuscript
  contribution is "framework vs typed-SNP methods", not "SARGENT vs
  ARG-BLUP-uniform". See `test/eval/block1_1_design.md` §§ "v1 results",
  "v2 results", "SARGENT cross-ancestry rescue test".
- **Block 1.2** (localization robustness sweep): v1 done on laptop
  (commit `2c19ba3`). **Verdict: estimator above chance in 5/6 cells**;
  estimator-ceiling gap small (0.04-0.11) across the grid. C1 HELD;
  C4 [MINOR] → HELD; C5 weakens. See
  `test/eval/block1_2_design.md` § "v1 results".
- **Block 1.3** (`fit_blup` h² bias): **audit + patch done** (commits
  `53af237`, `ebf757c`). `predict.py` double-centers K
  (`center_kinship`); h² bias +0.26 → +0.05; R² unchanged.
- **Block 0** (claims triage): Adam signed off 2026-05-29. Bucket
  moves C1 → HELD; C4 → HELD [MINOR]; C5 → WEAK (Block 1.2) and
  D2 → FALSIFIED at current estimator (Block 2.5.2) still awaiting
  Adam's sign-off on `CLAIMS.md` patch.

**Block 2 is now the Applications arc** (replaces the deferred sparse
fused lasso, demoted to parking lot). Localization confirmed as the
SARGENT-distinctive strength; the applications direction builds on
the per-tree variance primitive. Sub-blocks:

- **2.1** Allele-age decomposition (highest priority new experiment)
- **2.2** Partitioned heritability across functional annotations,
  three-estimator framework: (A) typed-only S-LDSC, (B) ARG-only,
  (C) hybrid with typed × untyped split per annotation
- **2.3** ARG-aware fine-mapping (likely separate paper)
- **2.4** "Linkage without pedigrees" framing for the manuscript
- **2.5** Polygenicity / effect-size architecture (M_e,ARG,
  variance-capture curve, per-branch effect-size — the HOW axis)

Block 3 (manuscript, Adam-led, next week), Block 4 (efficiency),
Block 5 (publication readiness) follow.

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

### NEVER kill processes by command-line pattern

`pkill -f <pat>`, `pgrep -f <pat>` and `awk '/<pat>/ {print $1}'` all match **the
shell running the command**, because that shell's own command line contains the
pattern. The result is that the session kills itself (exit 143/144) and, worse,
can kill unrelated processes that happen to match. This has happened repeatedly
in this project, including once to a *different session's* jobs.

**The rule:** never put the target's name in the same command that kills.

```bash
# WRONG -- matches this very shell, kills the session
pkill -f myjob.sh
for p in $(pgrep -f myjob.py); do kill $p; done

# RIGHT -- three steps, and every PID verified before it is signalled
pgrep -u "$USER" -x sargent                 # exact comm, never -f
cat /proc/<PID>/cmdline | tr '\0' ' '       # confirm it is the intended job
kill -TERM <PID>                            # explicit numeric PID, separate call
```

Use `pgrep -x <comm>` (exact process name) rather than `-f`. When a pattern is
genuinely needed, read `/proc/<pid>/cmdline` for each candidate and skip `$$`.
Prefer PID files written at launch. Verify identity (uid, `comm`, cmdline)
immediately before signalling, and never kill anything you did not start.

**`comm` is truncated to 15 characters.** `pgrep -x` matches against `comm`, so
an executable with a longer name silently matches *nothing* -- and "no output"
reads exactly like "not running":

```bash
pgrep -u "$USER" -x check_tree_matvec   # -> nothing; the name is 17 chars
pgrep -u "$USER" -x check_tree_matv     # -> the actual PIDs
```

Check `cat /proc/<pid>/comm` for the real (truncated) name before concluding a
job has finished or died. A wrong "nothing is running" has twice led to
relaunching work that was already in flight, and would be far worse if it led
to a `-f` fallback.

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
