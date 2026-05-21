# Phase 3.5 — Postdoc-pipeline compatibility

**Status**: ingest pipeline working end-to-end as of 2026-05-20.
Sequential SARGENT smoke test passed on a Nurdan Fig1 case;
joint-EM smoke test pending.

## Goal

Let SARGENT ingest the Nurdan simulation-pipeline output directly so the
postdoc can run experiments in parallel without an in-house simulator
detour.

## Design

The translation is intentionally thin — sargent reads PLINK BED and
tskit `.ts` natively, and we now read multi-column phenotype files
natively too. A one-time symlink + metadata pass is all that's required.

### Native support added to sargent

- `--pheno-col N` (`-C N`) on `sargent`: pick a 1-indexed column from the
  phenotype file. Default (`0`) reads the *last* whitespace-separated
  field on each line (legacy single-column behaviour).
- Implementation: `read_phenotypes` in `src/progs/sargent.c` was extended
  to take a column index and selects either the N-th or last field.
- Smoke test for the flag: `sargent --pheno-col 3 ...` on a Nurdan
  `all.vc_balanced.pheno` (replicate 1 of 5) returns the same value as
  the column-1 truth would expect.

### Adapter: `scripts/import_nurdan.py`

Creates a `case.*` layout for `predict.py` / `run_suite.py`:

| Output | Contents |
|---|---|
| `case.bed/bim/fam` | symlink to `data/all.{bed,bim,fam}` |
| `case.trees` | symlink to `data/all.ts` |
| `case.phen` | symlink to `pheno/all.<trait>.pheno` (multi-column; use `--pheno-col`) |
| `case.params.json` | translated from `pheno/all.<trait>.truth.csv` |

No data is copied. The adapter prints the exact sargent command to use
(including the right `--pheno-col`, which is `replicate + 2` because of
the two leading ID columns in the `.pheno` format).

### What we deliberately did NOT do

- **Train/test split**: the postdoc's convention is "indv 1–150 train,
  151–200 test", but we don't yet have BED/TS subsetting in the
  pipeline. Deferred until OOS comparison is needed; first prove sargent
  produces sensible in-sample fits on her data.
- **Per-tree truth from `kernel_weights.csv`**: Nurdan's pipeline
  deduplicates trees and stores weights at the unique-kernel level (163
  unique vs 2739 raw trees in the example case). Mapping back to
  per-tree weights via `Cmatrix_unique_idx.csv` is straightforward but
  unnecessary for first-pass validation.

## Smoke test (Fig1, 200 ind, 1 Mb, balanced h²=0.4/0.4)

Source: `~/test-sargent/nurdan/Simulation_Fig1_200ind_CausalRatio_0.1_Seed6903_1Mb/`,
trait `vc_balanced`, replicate 1.

```bash
python scripts/import_nurdan.py \
    --src ~/test-sargent/nurdan/Simulation_Fig1_200ind_CausalRatio_0.1_Seed6903_1Mb \
    --trait vc_balanced --pheno-col 1 \
    --out /tmp/nurdan_smoke/case
sargent --pheno-col 3 \
    /tmp/nurdan_smoke/case.phen \
    /tmp/nurdan_smoke/case.trees \
    /tmp/nurdan_smoke/case
```

### Results

Truth: h²_snp=0.4, h²_arg=0.4, σ²=0.2, h²_total=0.8 (200 indv × 2739
trees × 2488 SNPs).  Sequential SARGENT recovers total h² to within 3%;
SNP-vs-ARG split is biased toward SNP, as expected from Phase 3
(sequential always parks ~0.75–0.85 on h²_snp, irrespective of truth).

| Mode | h²_snp | h²_arg | h²_total | σ² | wall time |
|---|---|---|---|---|---|
| sequential   | 0.768 | 0.057 | **0.826** | 0.174 | 21m23s |
| joint EM 20  | TODO  | TODO  | TODO      | TODO  | TODO (running) |

A `--em 50` run on this case was killed after 3h21m without finishing —
the default 50 outer iterations is too aggressive for cases of this
size.  Smoke-test convention going forward: use `--em 20` on Nurdan-scale
inputs.  Phase 4 profiling on the desktop will tell us where the per-iter
cost lives and inform whether the default needs to scale with n_trees.

## How the postdoc should use this

```bash
# One-time per case directory:
python scripts/import_nurdan.py \
    --src <her case dir> \
    --trait vc_balanced \
    --pheno-col <1..5> \
    --out <some prefix>/case

# Then for each fit:
sargent --pheno-col <printed in adapter output> \
    <prefix>/case.phen <prefix>/case.trees <prefix>/case
```

Run the same `import_nurdan.py` invocation with `--pheno-col 2`, `3`, …
to point at each replicate.

## Open items

- [ ] Train/test split (when she wants OOS prediction): add `--n-test N`
  back to the adapter and use PLINK `--keep` + tskit `simplify()` to
  produce subset filesets.
- [ ] Accuracy comparison harness: read `snp_effects.csv` and the
  expanded per-tree weights from `kernel_weights.csv` ×
  `Cmatrix_unique_idx.csv`, compare to sargent's posterior β and ws².
- [ ] Add a case in `test/eval/run_suite.py` for Nurdan-imported cases
  so the same aggregator picks them up.
