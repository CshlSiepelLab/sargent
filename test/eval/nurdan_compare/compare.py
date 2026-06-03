#!/usr/bin/env python3
"""compare.py — tabulate SARGENT vs MKLMM results from JSON outputs.

Usage:
    python compare.py <sargent_summary.json> <nurdan_result.json>
"""
import json, sys, pathlib

def load(path):
    return json.loads(pathlib.Path(path).read_text())

def fmt(v):
    if v is None: return "  n/a  "
    return f"{v:7.4f}"

def main():
    if len(sys.argv) < 3:
        sys.exit("Usage: compare.py <sargent_summary.json> <nurdan_result.json>")
    s = load(sys.argv[1])
    n = load(sys.argv[2])

    truth_snp   = n.get("h2_snp_truth",   "?")
    truth_arg   = n.get("h2_arg_truth",   "?")
    truth_total = n.get("h2_total_truth", "?")
    truth_s2    = n.get("sigma2_truth",   "?")

    print(f"\n{'='*60}")
    print(f"  SARGENT vs MKLMM (ARG-LMM) — {n.get('trait','')} rep{n.get('pheno_col','')}")
    print(f"  n_indiv={n.get('n_indiv','?')}  "
          f"n_snps_sargent={s.get('n_variants','?')}  "
          f"n_snps_mklmm={n.get('n_snps','?')}  "
          f"n_trees_sargent={s.get('n_trees','?')}  "
          f"n_kernels_mklmm={n.get('n_kernels','?')}")
    print(f"{'='*60}")
    print(f"{'Quantity':<18}  {'Truth':>7}  {'SARGENT':>7}  {'MKLMM':>7}  {'diff(S-M)':>9}")
    print(f"{'-'*60}")

    rows = [
        ("h2_snp",   truth_snp,   s.get("h2_snp"),       n.get("h2_snp_est")),
        ("h2_arg",   truth_arg,   s.get("h2_arg"),        n.get("h2_arg_est")),
        ("h2_total", truth_total, s.get("h2"),            n.get("h2_total_est")),
        ("sigma2",   truth_s2,    s.get("sigma2"),        n.get("sigma2_est")),
    ]
    for label, truth, sv, mv in rows:
        diff = f"{sv-mv:+8.4f}" if sv is not None and mv is not None else "    n/a "
        t = fmt(truth) if isinstance(truth, float) else f"  {truth}  "
        print(f"  {label:<16}  {t}  {fmt(sv)}  {fmt(mv)}  {diff}")

    print(f"{'-'*60}")
    wt_s = s.get("wall_s") or s.get("wall_total_s")  # sargent summary doesn't include wall time
    wt_m = n.get("wall_total_s")
    wt_fit_m = n.get("wall_fit_s")
    print(f"  {'wall_total (s)':<16}  {'':>7}  {'(see below)':>7}  {fmt(wt_m)}  ")
    print(f"  {'wall_fit (s)':<16}  {'':>7}  {'':>7}  {fmt(wt_fit_m)}  ")
    print(f"\n  MKLMM backfit iters: {n.get('n_backfit_iters','?')}")
    print(f"  MKLMM active kernels (est/truth): "
          f"{n.get('n_active_kernels_est','?')} / {n.get('n_active_kernels_truth','?')}")
    print(f"{'='*60}\n")

if __name__ == "__main__":
    main()
