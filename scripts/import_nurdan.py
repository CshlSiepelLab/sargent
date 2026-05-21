#!/usr/bin/env python3
"""Import a Nurdan-pipeline simulation case into the sargent case.* layout.

The translation is intentionally thin: we DO NOT rewrite genotypes or
phenotypes.  Sargent reads PLINK BED, tskit .ts, and multi-column
phenotype files natively (use sargent --pheno-col N to pick a column).
All this script does is:

  * symlink the BED triple and the tree sequence into the case.* layout
    that predict.py / run_suite.py expect
  * read truth.csv for the chosen phenotype replicate and write a
    case.params.json with the field names sargent / predict.py use

Usage:
    python import_nurdan.py
        --src   <Nurdan case dir>
        --trait <vc_balanced | vc_highfix_lowrand | ...>
        --pheno-col <1..K>      # which replicate column (within truth.csv);
                                #   passed verbatim to sargent's --pheno-col
                                #   after accounting for the 2-column ID prefix
        --out   <case prefix>

Notes
-----
- Nurdan .pheno format is "family_id  indv_id  y1 y2 ... yK".  The
  --pheno-col semantics here are 1-indexed *over the phenotype replicates*
  (not over the file columns); the truth.csv "phenotype" column uses the
  same indexing.  The corresponding sargent flag is --pheno-col (k+2).
- We do NOT split train/test.  Run sargent in-sample first; once the
  postdoc has a use case for OOS prediction we can add BED/TS subsetting.
"""
import argparse
import csv
import json
from pathlib import Path
import sys


def parse_truth(path: Path, phenotype_1indexed: int):
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if int(row["phenotype"]) == phenotype_1indexed:
                return row
    sys.exit(f"ERROR: phenotype {phenotype_1indexed} not found in {path}")


def n_indiv_from_fam(fam_path: Path) -> int:
    with open(fam_path) as f:
        return sum(1 for line in f if line.strip())


def n_variants_from_bim(bim_path: Path) -> int:
    with open(bim_path) as f:
        return sum(1 for line in f if line.strip())


def symlink_or_replace(src: Path, dst: Path):
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src.resolve())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True,
                    help="Nurdan case directory (containing data/, pheno/, kernels/)")
    ap.add_argument("--trait", required=True,
                    help="Trait variant name (e.g. vc_balanced, vc_highfix_lowrand)")
    ap.add_argument("--pheno-col", type=int, required=True,
                    help="Phenotype replicate to use, 1-indexed (matches truth.csv 'phenotype')")
    ap.add_argument("--out", required=True,
                    help="Output prefix (e.g. /path/to/case)")
    args = ap.parse_args()

    src = Path(args.src).resolve()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    bed   = src / "data"  / "all.bed"
    bim   = src / "data"  / "all.bim"
    fam   = src / "data"  / "all.fam"
    ts    = src / "data"  / "all.ts"
    pheno = src / "pheno" / f"all.{args.trait}.pheno"
    truth = src / "pheno" / f"all.{args.trait}.truth.csv"

    for p in (bed, bim, fam, ts, pheno, truth):
        if not p.exists():
            sys.exit(f"ERROR: missing {p}")

    n_indiv    = n_indiv_from_fam(fam)
    n_variants = n_variants_from_bim(bim)

    truth_row = parse_truth(truth, args.pheno_col)
    h2_snp_target = float(truth_row["h2_fix"])
    h2_arg_target = float(truth_row["h2_rand"])
    sigma2        = float(truth_row["sigma2"])
    h2_total      = float(truth_row["h2_total"])

    # Symlink genotypes, tree sequence, AND the raw multi-column .pheno file.
    # Users pass case.phen + --pheno-col N to sargent.
    symlink_or_replace(bed,   Path(f"{args.out}.bed"))
    symlink_or_replace(bim,   Path(f"{args.out}.bim"))
    symlink_or_replace(fam,   Path(f"{args.out}.fam"))
    symlink_or_replace(ts,    Path(f"{args.out}.trees"))
    symlink_or_replace(pheno, Path(f"{args.out}.phen"))

    # The 2 leading ID columns in the .pheno file mean replicate k maps to
    # sargent --pheno-col (k + 2).
    sargent_pheno_col = args.pheno_col + 2

    params = {
        "n_indiv":         n_indiv,
        "n_test":          0,
        "n_total":         n_indiv,
        "n_variants":      n_variants,
        "h2_snp_target":   h2_snp_target,
        "h2_arg_target":   h2_arg_target,
        "h2_total_target": h2_total,
        "sigma2_target":   sigma2,
        "source":          "nurdan_import",
        "source_case":     str(src),
        "trait":           args.trait,
        "pheno_replicate": args.pheno_col,
        "sargent_pheno_col": sargent_pheno_col,
        "n_active_kernels_truth": int(truth_row["n_active_kernels"]),
    }
    with open(f"{args.out}.params.json", "w") as f:
        json.dump(params, f, indent=2)

    print(f"wrote case at {args.out}.*")
    print(f"  n_indiv    = {n_indiv}")
    print(f"  n_variants = {n_variants}")
    print(f"  truth      = h2_snp={h2_snp_target}  h2_arg={h2_arg_target}  "
          f"sigma2={sigma2}  h2_total={h2_total}")
    print(f"")
    print(f"Run sargent with:")
    print(f"  sargent --pheno-col {sargent_pheno_col} \\")
    print(f"          {args.out}.phen {args.out}.trees {args.out}")


if __name__ == "__main__":
    main()
