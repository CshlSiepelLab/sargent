#!/usr/bin/env python3
"""
run_suite.py — run SARGENT on every dataset in the test suite and collect metrics.

For each case in <suitedir>/<scale>/<case_name>/ it:
  1. Runs sargent with --ts-geno (or optionally with BED files)
  2. Parses stderr for sigma2, h2_fixed, h2_random
  3. Computes ws2_corr, beta_corr, and variance accuracy
  4. Writes a JSON result per case and a summary TSV

Usage:
    python run_suite.py [--suitedir DIR] [--sargent PATH] [--scales small medium]
                        [--cases ...] [--nprobe N] [--prior {bb,horseshoe}]
                        [--pi0 VAL] [--ts-geno] [--outdir DIR]
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
from metrics import (load_params, load_true_ws2, load_true_beta,
                     load_inferred_ws2, load_inferred_beta, load_inferred_rho,
                     ws2_metrics, beta_metrics, variance_metrics)

CASE_NAMES = ["fixed_only", "fixed_heavy", "balanced", "rand_heavy", "rand_only"]


def find_sargent(given):
    if given:
        return given
    repo_root = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
    candidates = [
        os.path.join(repo_root, "build", "src", "progs", "sargent"),
        os.path.join(repo_root, "build", "sargent"),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    # Try PATH
    import shutil
    found = shutil.which("sargent")
    if found:
        return found
    return None


def parse_stderr(stderr_text):
    """Extract sigma2, h2_fixed, h2_random from sargent stderr.

    Matches both per-iteration format (sigma2=X.Y, h2_fixed=X.Y) and the
    final summary block format (sigma2 : X.Y, h2.fixed : X.Y, h2.random : X.Y).
    Latest match wins (so summary-block values override per-iteration ones).
    """
    result = {}
    for line in stderr_text.splitlines():
        m = re.search(r'h2[\._]fixed\s*[=:]\s*([0-9.eE+-]+)', line)
        if m:
            result["h2_fixed"] = float(m.group(1))
        m = re.search(r'h2[\._]random\s*[=:]\s*([0-9.eE+-]+)', line)
        if m:
            result["h2_random"] = float(m.group(1))
        m = re.search(r'sigma2\s*[=:]\s*([0-9.eE+-]+)', line)
        if m:
            result["sigma2"] = float(m.group(1))
        # Plain "h2 : X.Y" (must not match h2_fixed/h2_random/h2.fixed/h2.random)
        m = re.search(r'(?<![._a-zA-Z])h2\s*:\s*([0-9.eE+-]+)', line)
        if m and "fixed" not in line and "random" not in line:
            result["h2"] = float(m.group(1))
    return result


def run_one(case_dir, case_name, sargent_bin, args, outdir):
    prefix = os.path.join(case_dir, case_name)
    out_prefix = os.path.join(outdir, case_name)

    if not os.path.exists(prefix + ".trees"):
        print(f"  [skip] {case_name}: no .trees file found", flush=True)
        return None

    params = load_params(prefix)
    has_fixed  = params["h2_fixed_target"]  > 0
    has_random = params["h2_random_target"] > 0

    # Note: do NOT pass --silent — we need the summary block in stderr.
    cmd = [sargent_bin]
    if args.ts_geno:
        cmd += ["--ts-geno"]
    if args.nprobe > 0:
        cmd += ["--nprobe", str(args.nprobe)]
    if args.prior == "horseshoe":
        cmd += ["--horseshoe"]
    elif args.pi0 is not None:
        cmd += ["--pi0", str(args.pi0)]
    if not has_random:
        cmd += ["--fixed-only"]
    elif not has_fixed:
        cmd += ["--random-only"]
    cmd += ["--outroot", out_prefix]
    cmd += [prefix + ".phen", prefix + ".trees"]
    if not args.ts_geno:
        cmd += [prefix]

    t0 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - t0

    if result.returncode != 0:
        print(f"  ERROR: sargent failed for {case_name} (exit {result.returncode})",
              file=sys.stderr)
        print(result.stderr[-2000:], file=sys.stderr)
        return {"case": case_name, "error": result.stderr[-500:], "elapsed_s": elapsed}

    parsed = parse_stderr(result.stderr + result.stdout)

    record = {
        "case":        case_name,
        "scale":       args._scale,
        "regime":      params.get("regime", "unknown"),
        "prior":       args.prior,
        "elapsed_s":   round(elapsed, 2),
        "n_indiv":     params["n_indiv"],
        "length_bp":   int(params["length"]),
        "n_trees":     params["n_trees"],
        "n_variants":  params["n_variants"],
        "n_causal":    params["n_causal"],
        "h2_fixed_target":   params["h2_fixed_target"],
        "h2_random_target":  params["h2_random_target"],
        "sigma2_true":       params["sigma2"],
    }

    # Variance metrics (from sargent summary output)
    record["sigma2_inferred"]   = parsed.get("sigma2", float("nan"))
    record["h2_fixed_inferred"] = parsed.get("h2_fixed", float("nan"))
    record["h2_rand_inferred"]  = parsed.get("h2_random", float("nan"))

    vm = variance_metrics(
        sigma2_inferred=record["sigma2_inferred"],
        sigma2_true=record["sigma2_true"],
        h2_fixed_inferred=record["h2_fixed_inferred"],
        h2_fixed_true=params["h2_fixed_target"],
        h2_rand_inferred=record["h2_rand_inferred"],
        h2_rand_true=params["h2_random_target"],
    )
    record.update(vm)

    # ws2 metrics
    if has_random and os.path.exists(out_prefix + ".w2.tsv"):
        true_ws2     = load_true_ws2(prefix)
        inferred_ws2 = load_inferred_ws2(out_prefix)
        inferred_rho = load_inferred_rho(out_prefix)
        if len(inferred_ws2) == len(true_ws2):
            record.update(ws2_metrics(true_ws2, inferred_ws2, inferred_rho))
        else:
            record["ws2_error"] = (f"length mismatch: {len(inferred_ws2)} "
                                   f"vs {len(true_ws2)}")

    # beta metrics
    if has_fixed and os.path.exists(out_prefix + ".beta.tsv"):
        true_beta     = load_true_beta(prefix, params["n_variants"])
        inferred_beta = load_inferred_beta(out_prefix)
        if len(inferred_beta) == params["n_variants"]:
            record.update(beta_metrics(true_beta, inferred_beta))
        else:
            record["beta_error"] = (f"length mismatch: {len(inferred_beta)} "
                                    f"vs {params['n_variants']}")

    return record


def print_summary(records):
    fields = ["case", "sigma2_inferred", "sigma2_true", "sigma2_relerr",
              "h2_fixed_inferred", "h2_rand_inferred",
              "ws2_corr_active", "beta_corr_causal", "elapsed_s"]

    header = "\t".join(fields)
    print(header, flush=True)
    for r in records:
        if r is None:
            continue
        row = []
        for f in fields:
            v = r.get(f, "")
            if isinstance(v, float):
                row.append(f"{v:.4g}")
            else:
                row.append(str(v))
        print("\t".join(row), flush=True)


def write_markdown(records, path, args):
    """Write a human-readable markdown report."""
    if not records:
        return
    with open(path, "w") as f:
        f.write(f"# SARGENT test-suite report\n\n")
        f.write(f"- Prior: **{args.prior}**\n")
        if args.pi0 is not None:
            f.write(f"- pi0: **{args.pi0}**\n")
        f.write(f"- nprobe: **{args.nprobe}**\n")
        f.write(f"- ts-geno: **{args.ts_geno}**\n\n")

        # Group by scale
        scales = []
        for r in records:
            if r and r.get("scale") not in scales:
                scales.append(r["scale"])

        for scale in scales:
            scale_recs = [r for r in records if r and r.get("scale") == scale]
            if not scale_recs:
                continue
            ex = scale_recs[0]
            f.write(f"## Scale `{scale}`  (n_indiv={ex['n_indiv']}, "
                    f"length={ex['length_bp']:,} bp, regime={ex['regime']})\n\n")

            # Variance accuracy
            f.write("### Variance components\n\n")
            f.write("| case | trees | vars | h²f true | h²f inf | h²r true | h²r inf "
                    "| σ² true | σ² inf | elapsed |\n")
            f.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
            for r in scale_recs:
                f.write(f"| {r['case']} | {r['n_trees']} | {r['n_variants']} | "
                        f"{r['h2_fixed_target']:.2f} | {fmt(r.get('h2_fixed_inferred'))} | "
                        f"{r['h2_random_target']:.2f} | {fmt(r.get('h2_rand_inferred'))} | "
                        f"{r['sigma2_true']:.2f} | {fmt(r.get('sigma2_inferred'))} | "
                        f"{r['elapsed_s']:.1f}s |\n")
            f.write("\n")

            # Random-effects (ws²) recovery
            f.write("### Random-effects (ws²) recovery\n\n")
            f.write("| case | n_active true | n_called (rho≥0.5) | "
                    "AUC_rho | AUC_ws² | precision | recall | "
                    "ws²_corr_active |\n")
            f.write("|---|---:|---:|---:|---:|---:|---:|---:|\n")
            for r in scale_recs:
                if "n_active_true" not in r:
                    continue
                f.write(f"| {r['case']} | {r['n_active_true']} | "
                        f"{r.get('n_called_active','-')} | "
                        f"{fmt(r.get('auc_rho'))} | "
                        f"{fmt(r.get('auc_ws2'))} | "
                        f"{fmt(r.get('precision_rho'))} | "
                        f"{fmt(r.get('recall_rho'))} | "
                        f"{fmt(r.get('ws2_corr_active'))} |\n")
            f.write("\n")

            # Fixed-effects (beta) recovery
            f.write("### Fixed-effects (β) recovery\n\n")
            f.write("| case | n_causal | β_corr_causal | β_corr_all | β_MSE_causal |\n")
            f.write("|---|---:|---:|---:|---:|\n")
            for r in scale_recs:
                if "beta_corr_causal" not in r and "beta_corr_all" not in r:
                    continue
                f.write(f"| {r['case']} | {r.get('n_causal','-')} | "
                        f"{fmt(r.get('beta_corr_causal'))} | "
                        f"{fmt(r.get('beta_corr_all'))} | "
                        f"{fmt(r.get('beta_mse_causal'))} |\n")
            f.write("\n")


def fmt(x):
    if x is None:
        return "-"
    try:
        if isinstance(x, float):
            if x != x:   # NaN
                return "nan"
            return f"{x:.3g}"
    except Exception:
        pass
    return str(x)


def write_tsv(records, path):
    if not records:
        return
    all_keys = []
    for r in records:
        if r:
            for k in r:
                if k not in all_keys:
                    all_keys.append(k)
    with open(path, "w") as f:
        f.write("\t".join(all_keys) + "\n")
        for r in records:
            if r is None:
                continue
            row = []
            for k in all_keys:
                v = r.get(k, "")
                if isinstance(v, float):
                    row.append(f"{v:.6g}")
                else:
                    row.append(str(v))
            f.write("\t".join(row) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--suitedir", default=None,
                    help="base directory of test suite (default: ../../../test-sargent/suite)")
    ap.add_argument("--sargent",  default=None, help="path to sargent binary")
    ap.add_argument("--scales",   nargs="+", default=["small"],
                    choices=["debug", "small", "medium"],
                    help="scales to evaluate")
    ap.add_argument("--cases",    nargs="+", default=CASE_NAMES,
                    choices=CASE_NAMES, help="which cases to run")
    ap.add_argument("--nprobe",   type=int, default=10,
                    help="Hutchinson probe count (default 10)")
    ap.add_argument("--prior",    default="bb", choices=["bb", "horseshoe"],
                    help="Stage 2 prior (default: bb)")
    ap.add_argument("--pi0",      type=float, default=None,
                    help="pi0 value for BB prior (default: none)")
    ap.add_argument("--ts-geno",  action="store_true",
                    help="use tree-sequence genotypes instead of BED files")
    ap.add_argument("--outdir",   default=None,
                    help="directory for sargent output files (default: <suitedir>/results)")
    args = ap.parse_args()

    if args.suitedir is None:
        args.suitedir = os.path.join(os.path.expanduser("~"), "test-sargent", "suite")

    sargent_bin = find_sargent(args.sargent)
    if not sargent_bin or not os.path.isfile(sargent_bin):
        sys.exit(f"ERROR: cannot find sargent binary (tried: {sargent_bin})")

    if args.outdir is None:
        args.outdir = os.path.join(args.suitedir, "results")

    print(f"Suite dir : {args.suitedir}", flush=True)
    print(f"Sargent   : {sargent_bin}", flush=True)
    print(f"Output    : {args.outdir}", flush=True)
    print(flush=True)

    all_records = []
    for scale in args.scales:
        args._scale = scale
        results_dir = os.path.join(args.outdir, scale)
        os.makedirs(results_dir, exist_ok=True)

        print(f"=== Scale: {scale} ===", flush=True)
        scale_records = []

        for case_name in args.cases:
            case_dir = os.path.join(args.suitedir, scale, case_name)
            if not os.path.isdir(case_dir):
                print(f"  [skip] {case_name}: directory not found", flush=True)
                continue

            print(f"  Running {case_name} ...", flush=True)
            rec = run_one(case_dir, case_name, sargent_bin, args, results_dir)
            if rec:
                scale_records.append(rec)
                # Save individual JSON result
                json_path = os.path.join(results_dir, case_name + ".result.json")
                with open(json_path, "w") as f:
                    json.dump(rec, f, indent=2)

        all_records.extend(scale_records)

        # Write per-scale TSV summary
        tsv_path = os.path.join(results_dir, "summary.tsv")
        write_tsv(scale_records, tsv_path)
        print(f"  Written: {tsv_path}", flush=True)
        print(flush=True)

    print("=== Summary ===", flush=True)
    print_summary(all_records)

    # Write combined TSV + markdown report
    combined_path = os.path.join(args.outdir, "all_results.tsv")
    write_tsv(all_records, combined_path)
    md_path = os.path.join(args.outdir, "report.md")
    write_markdown(all_records, md_path, args)
    print(f"\nAll results: {combined_path}", flush=True)
    print(f"Markdown report: {md_path}", flush=True)


if __name__ == "__main__":
    main()
