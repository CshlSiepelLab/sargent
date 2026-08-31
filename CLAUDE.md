# SARGENT — Claude Code project guide

This file is auto-loaded by Claude Code in this repo. It is the entry point
for any session, including server-side sessions launched via SSH.

## What this project is

**SARGENT** (Sparse ARG-based ENcoding of Traits): an ARG-based polygenic
heritability and prediction tool. C core (production estimator) plus its tests.
PI: Adam Siepel (CSHL).

**This repo is the software only.** Commit `3739ec2` (2026-06-29, "Make sargent
pure software development") moved every research document out. Earlier versions
of this file pointed at `ROADMAP.md`, `CLAIMS.md`,
`test/eval/phase10_oracle.md`, `test/eval/block1_1_design.md` and
`test/eval/fit_blup_audit.md` — **all five are gone from this repo.** Do not go
looking for them here.

The live methodological documents are in the sibling repo
`~/git-src/sargent-paper`:

- **`notes/PHI_PROGRAM.md`** — the φ axis (φ_src = h²_untyped/h²_total as the
  design knob, φ_tag as the measured quantity). The current experimental frame.
- **`notes/SPEED_PLAN.md`** — kinship-build performance, including what was
  measured and which two confident diagnoses were wrong.
- **`analysis/RUNTIME.md`** — per-method runtimes and the Table S1 numbers.
- **`notes/CLI_REDESIGN.md`** — the `--model` flavour CLI.

The manuscript lives in `~/newwriting/arg-mapping` (`main.tex`,
`supplement.tex`, `figs/`).

## Current state (as of 2026-08-31)

Branch `main`, **sargent v0.9.0**, everything pushed.

What the estimator does now, and what changed recently:

- **`--model <flavour>`** selects one arm and runs only its prerequisites:
  `snp-lmm`, `sargent-f`, `sargent-ur`, `sargent-u`, `sargent-w`,
  `sargent-wf`, `all`. Before this, every invocation computed every arm, so no
  per-method cost could be attributed. `all` is the FIGURE set and deliberately
  excludes `sargent-w`/`sargent-wf`.
- **h² is a realized-variance quantity**, adapted from BSLMM's PVE (eq. 13–15):
  the variance of the *sum* Xβ̃+u over a model-implied denominator, not the
  plug-in σ²_u/(σ²_u+σ²_e). The plug-in inflates under kernel misspecification;
  mean |bias| went from +0.104/+0.170 to 0.009/0.012.
- **The ARG kinship build is ~8x faster** (v0.8.0). The readout re-derived each
  sample's root path once per descendant; a time-ordered prefix sweep plus
  column-level parallelism took the benchmark cell 185.05 s → 22.75 s at the
  production 2-thread config. The tsLMM gap went 34x → 4.0x.
- **SARGENT-R is retired** (v0.9.0). It was the single best arm in 13.2% of
  6452 cells with a median gain over SARGENT-U of exactly 0.0000. Localization
  moved to the per-branch score test (Fig 5), which never used it. The ws²
  TABLE is still written — Blocks 2.1/2.2 read `est.w2.tsv` via the legacy
  `--arg-only` path.

**Figure status.** The within-population and cross-ancestry φ figures are
current. `main.tex:1153` still carries the caveat that *the density supplement
and the reconstructed-ARG figures show the plug-in h²*; the reconstructed-ARG
half is being cleared. **The density figure plots `sargent` and `sargent_r`
series that the retired arm produced** — re-running those cells at v0.9.0 would
delete two series from a supplementary figure, so it needs a decision first,
not a re-run.

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

### Build and check a version

```bash
cmake -B build -S . && cmake --build build -j
build/src/progs/sargent --version
```

**The version string is read at CMake CONFIGURE time**, so a plain rebuild after
a version bump keeps reporting the old number. If it looks wrong, delete
`build/CMakeCache.txt` and reconfigure. This has caused a 420-cell run to be
launched against a binary reporting 0.6.0.

### Time a method

`~/git-src/sargent-paper/analysis/timing_*.sh`. Measure every row in ONE pass on
an idle machine — an unrelated job on other cores moved BSLMM by 42% between
passes, which is how a published table came to have rows that were not
comparable to one another.

### Fit one cell

```bash
sargent --silent --nprobe 10 --threads 2 --predict --model all \
        --outroot <dir>/est <dir>/case.phen <dir>/case.trees <dir>/case
```

Aggregators prefer `est_pve2` > `est_pve` > `est_joint` (see
`analysis/branch_select/aggregate_fig2_phi.py`), so a re-fit that should
supersede an older pass writes `est_pve2`.

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

This repo holds the C estimator and its tests, nothing else. Analysis prose and
eval infrastructure belong in `~/git-src/sargent-paper`.

- **Findings from a finished experiment** → append a section to the relevant
  `sargent-paper/notes/*.md` (e.g. `PHI_PROGRAM.md`, `SPEED_PLAN.md`), or
  `analysis/RUNTIME.md` for timings.
- **Manuscript text and figures** → `~/newwriting/arg-mapping`.
- **Reusable analysis or benchmark scripts** → `sargent-paper/analysis/`,
  committed. A `Reproduce` block naming a path in a session scratchpad is worse
  than none: `RUNTIME.md` cited `run_ext_cell.sh` for weeks and that script has
  never existed.
- **Scratch / one-off diagnostics** → the session scratchpad (not any repo).
- **Changes to the C core** → here, with a regression against the previous
  release before merging.

## Conventions

- Per-cell output convention: `<BASE>/<cell-name>/{case.*, est.*, *.log,
  est.predict.json}` (one directory per parameter cell).
- Per-cell parameter encoding in directory name, e.g.
  `mix0.05_typed0.70_s3`, or `phig064_s1` on the φ grid — the aggregator regexes
  expect these.
- `simulate.py` outputs follow the prefix convention `<prefix>.{phen,
  trees, full.trees, bed, bim, fam, params.json, true_beta.tsv,
  true_ws2.tsv, train.anc.tsv, test.anc.tsv}`.
- Analysis `*.md` docs end with a `Reproduce` block giving exact invocations,
  naming COMMITTED scripts only.
- Aggregators prefer `est_pve2` > `est_pve` > `est_joint`, so a re-fit meant to
  supersede an earlier pass writes `est_pve2` and leaves the old one intact.

## If something feels off

If results are surprising, parameters seem wrong, or you suspect a bug
masking a result: stop, write up what you see in a short message, and
ask before propagating downstream. Adam has been burned by overconfident
intermediate steps before; an extra clarifying question costs less than
a wrong-direction sweep.
