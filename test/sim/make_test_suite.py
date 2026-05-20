#!/usr/bin/env python3
"""
make_test_suite.py — generate a structured suite of simulation datasets for SARGENT testing.

Produces five fixed:random h² splits at two scales (small / medium), all with total h²=0.8.
Each dataset is a subdirectory of <outdir>/<scale>/<case_name>/.

Usage:
    python make_test_suite.py [--outdir DIR] [--scales small medium] [--seeds-offset N]

Output structure:
    <outdir>/small/snp_only/     n=100, h2_snp=0.80, h2_arg=0.00
    <outdir>/small/snp_heavy/    n=100, h2_snp=0.60, h2_arg=0.20
    <outdir>/small/balanced/       n=100, h2_snp=0.40, h2_arg=0.40
    <outdir>/small/arg_heavy/     n=100, h2_snp=0.20, h2_arg=0.60
    <outdir>/small/arg_only/      n=100, h2_snp=0.00, h2_arg=0.80
    <outdir>/medium/*              n=500, same splits
"""

import argparse
import json
import os
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SIMULATE_PY = os.path.join(SCRIPT_DIR, "..", "..", "scripts", "simulate.py")

# Five h² splits; total h²=0.8, sigma2=0.2
CASES = [
    {"name": "snp_only",  "h2_snp": 0.80, "h2_arg": 0.00},
    {"name": "snp_heavy", "h2_snp": 0.60, "h2_arg": 0.20},
    {"name": "balanced",    "h2_snp": 0.40, "h2_arg": 0.40},
    {"name": "arg_heavy",  "h2_snp": 0.20, "h2_arg": 0.60},
    {"name": "arg_only",   "h2_snp": 0.00, "h2_arg": 0.80},
]

SCALES = {
    "debug": {
        "n_indiv":     300,
        "length":      100_000,
        "recomb_rate": 1e-8,
        "mut_rate":    1e-8,
        "Ne":          10_000,
        "n_causal":    5,
        "seed_base":   1,
    },
    "small": {
        "n_indiv":     100,
        "length":      500_000,
        "recomb_rate": 1e-8,
        "mut_rate":    1e-8,
        "Ne":          10_000,
        "n_causal":    10,
        "seed_base":   101,
    },
    "medium": {
        "n_indiv":     500,
        "length":      2_000_000,
        "recomb_rate": 1e-8,
        "mut_rate":    1e-8,
        "Ne":          10_000,
        "n_causal":    20,
        "seed_base":   201,
    },
}


def generate_one(case, scale_cfg, outdir, seed_offset, regime, n_test, dry_run, python_bin):
    name  = case["name"]
    scale_name = scale_cfg["_scale_name"]
    dest  = os.path.join(outdir, scale_name, name)
    os.makedirs(dest, exist_ok=True)
    prefix = os.path.join(dest, name)

    # Skip if already generated (params file exists)
    if os.path.exists(prefix + ".params.json") and not dry_run:
        print(f"  [skip] {scale_name}/{name} — already exists", flush=True)
        return

    seed = scale_cfg["seed_base"] + seed_offset + list(CASES).index(case)

    cmd = [
        python_bin, SIMULATE_PY,
        "--n-indiv",     str(scale_cfg["n_indiv"]),
        "--n-test",      str(n_test),
        "--length",      str(scale_cfg["length"]),
        "--recomb-rate", str(scale_cfg["recomb_rate"]),
        "--mut-rate",    str(scale_cfg["mut_rate"]),
        "--Ne",          str(scale_cfg["Ne"]),
        "--n-causal",    str(scale_cfg["n_causal"]),
        "--h2-snp",    str(case["h2_snp"]),
        "--h2-arg",   str(case["h2_arg"]),
        "--regime",      regime,
        "--seed",        str(seed),
        "--out",         prefix,
    ]

    if dry_run:
        print("  [dry] " + " ".join(cmd), flush=True)
        return

    print(f"  Generating {scale_name}/{name} (seed={seed}) ...", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: simulate.py failed for {scale_name}/{name}:", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        return

    # Load and echo realised h²
    with open(prefix + ".params.json") as f:
        p = json.load(f)
    print(f"    n_trees={p['n_trees']}, n_variants={p['n_variants']}, "
          f"h2_snp_realised={p['h2_snp_realised']:.3f}, "
          f"h2_arg_realised={p['h2_arg_realised']:.3f}",
          flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir",       default=None,
                    help="output base directory (default: ../../../test-sargent/suite)")
    ap.add_argument("--scales",       nargs="+", default=["small", "medium"],
                    choices=["debug", "small", "medium"],
                    help="which scales to generate (default: small medium)")
    ap.add_argument("--regime",       default="sparse",
                    choices=["sparse", "diffuse"],
                    help="ws² mixture regime passed to simulate.py (default: sparse)")
    ap.add_argument("--n-test",       type=int, default=0,
                    help="held-out test individuals per case (default: 0)")
    ap.add_argument("--seeds-offset", type=int, default=0,
                    help="add this to all seeds (for replicate runs)")
    ap.add_argument("--dry-run",      action="store_true",
                    help="print commands without executing")
    ap.add_argument("--python",       default=sys.executable,
                    help="Python interpreter to use for simulate.py")
    args = ap.parse_args()

    if args.outdir is None:
        # default: ~/test-sargent/suite
        args.outdir = os.path.join(os.path.expanduser("~"), "test-sargent", "suite")

    print(f"Output directory: {args.outdir}", flush=True)
    print(f"Scales: {args.scales}", flush=True)
    print(f"simulate.py: {SIMULATE_PY}", flush=True)
    print(flush=True)

    for scale_name in args.scales:
        cfg = dict(SCALES[scale_name])
        cfg["_scale_name"] = scale_name
        print(f"=== Scale: {scale_name} (n={cfg['n_indiv']}, L={cfg['length']:,}bp, "
              f"regime={args.regime}) ===", flush=True)
        for case in CASES:
            generate_one(case, cfg, args.outdir, args.seeds_offset, args.regime,
                         args.n_test, args.dry_run, args.python)
        print(flush=True)

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
