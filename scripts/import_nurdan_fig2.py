#!/usr/bin/env python3
"""Import a Nurdan Figure 2 (missing-data) simulation into sargent case.* layout.

Figure 2 sims differ from Figure 1:
  - phenotype files: all.pheno{1,2,3}.pheno  (no trait name, no truth.csv)
  - pheno_col selects a column within the phenotype file (3 = first replicate)
  - optional SNP downsampling: use --ds-factor N to pick all.N.down.bed/bim/fam
  - tree sequence is ALWAYS all.ts (full ARG regardless of SNP downsampling)

Usage:
    python import_nurdan_fig2.py
        --src      <Nurdan Fig2 case dir>
        --pheno    <1|2|3>        # which pheno file (all.pheno1.pheno etc.)
        --pheno-col <int>         # 1-indexed column in that file (3=rep1, 4=rep2, 5=rep3)
        --ds-factor <0|1|2|4|10> # SNP downsampling factor (0=all SNPs, default 0)
        --h2-truth  <float>       # known simulation h2_snp (for params.json)
        --out       <case prefix>
"""
import argparse
import json
from pathlib import Path
import sys


def symlink_or_replace(src: Path, dst: Path):
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    dst.symlink_to(src.resolve())


def n_indiv_from_fam(fam: Path) -> int:
    return sum(1 for line in fam.read_text().splitlines() if line.strip())


def n_variants_from_bim(bim: Path) -> int:
    return sum(1 for line in bim.read_text().splitlines() if line.strip())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src",       required=True)
    ap.add_argument("--pheno",     type=int, default=1,
                    help="Pheno file index: 1 -> all.pheno1.pheno (default 1)")
    ap.add_argument("--pheno-col", type=int, default=3,
                    help="1-indexed column within the pheno file (default 3 = first replicate)")
    ap.add_argument("--ds-factor", type=int, default=0,
                    help="SNP downsampling factor (0=all SNPs; 1/2/4/10=fractions). "
                         "Tree sequence is always the full all.ts.")
    ap.add_argument("--h2-truth",  type=float, default=None,
                    help="Known h2_snp from simulation (e.g. 0.8 for pheno1)")
    ap.add_argument("--out",       required=True)
    args = ap.parse_args()

    src = Path(args.src).resolve()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Tree sequence is always the full one
    ts  = src / "data" / "all.ts"

    # Genotype BED — either full or downsampled
    if args.ds_factor == 0:
        bed = src / "data" / "all.bed"
        bim = src / "data" / "all.bim"
        fam = src / "data" / "all.fam"
    else:
        bed = src / "data" / f"all.{args.ds_factor}.down.bed"
        bim = src / "data" / f"all.{args.ds_factor}.down.bim"
        fam = src / "data" / f"all.{args.ds_factor}.down.fam"

    pheno = src / "pheno" / f"all.pheno{args.pheno}.pheno"

    for p in (ts, bed, bim, fam, pheno):
        if not p.exists():
            sys.exit(f"ERROR: missing {p}")

    n_indiv    = n_indiv_from_fam(fam)
    n_variants = n_variants_from_bim(bim)

    symlink_or_replace(bed,   Path(f"{args.out}.bed"))
    symlink_or_replace(bim,   Path(f"{args.out}.bim"))
    symlink_or_replace(fam,   Path(f"{args.out}.fam"))
    symlink_or_replace(ts,    Path(f"{args.out}.trees"))
    symlink_or_replace(pheno, Path(f"{args.out}.phen"))

    params = {
        "n_indiv":         n_indiv,
        "n_test":          0,
        "n_total":         n_indiv,
        "n_variants":      n_variants,
        "h2_snp_target":   args.h2_truth if args.h2_truth is not None else None,
        "h2_arg_target":   0.0,       # Fig2 phenotypes are pure-SNP
        "h2_total_target": args.h2_truth if args.h2_truth is not None else None,
        "sigma2_target":   None,
        "source":          "nurdan_fig2_import",
        "source_case":     str(src),
        "pheno_file":      args.pheno,
        "pheno_col_file":  args.pheno_col,   # column index in the .phen file
        "ds_factor":       args.ds_factor,
        "sargent_pheno_col": args.pheno_col, # same: col in file = col for sargent
    }
    with open(f"{args.out}.params.json", "w") as f:
        json.dump(params, f, indent=2)

    print(f"wrote case at {args.out}.*")
    print(f"  n_indiv    = {n_indiv}")
    print(f"  n_variants = {n_variants} (ds_factor={args.ds_factor})")
    print(f"  tree seq   = all.ts (full, {n_indiv} individuals)")
    print(f"  h2_truth   = {args.h2_truth} (pure SNP, h2_arg=0)")
    print()
    print(f"Run sargent with:")
    print(f"  sargent --pheno-col {args.pheno_col} \\")
    print(f"          {args.out}.phen {args.out}.trees {args.out} \\")
    print(f"          -o {args.out}")


if __name__ == "__main__":
    main()
