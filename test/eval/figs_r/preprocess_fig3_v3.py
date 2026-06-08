#!/usr/bin/env python3
"""
preprocess_fig3_v3.py — emit block1_1_v3_results.csv from the v3 cell sweep.

Walks every cell under $BASE (default ~/test-sargent/block1_1_v3), parses the
cell name for (model, train_pop, test_pop, test_time, mix, typed, seed),
reads est.predict.json for r2_calib per arm, and case.params.json for the
stdpopsim-computed F_ST.  Emits a flat CSV in the same column shape as the
v2 results file so fig3.R can drop in.

CLI:
  python preprocess_fig3_v3.py [--base BASE] [--out OUT]
"""
import argparse
import csv
import glob
import json
import math
import os
import re
import sys


CELL_PAT = re.compile(
    r"([^_]+)_([^x]+)x([^_@]+)(?:@t(\d+))?_mix([\d.]+)_typed([\d.]+)_s(\d+)"
)


def r2c(arm_d):
    """Extract r²_calib from one arm's est.predict.json sub-dict; fall back to r²."""
    if not isinstance(arm_d, dict):
        return None
    v = arm_d.get("r2_calib")
    if isinstance(v, (int, float)) and not math.isnan(float(v)):
        return float(v)
    rr = arm_d.get("r")
    if isinstance(rr, (int, float)) and not math.isnan(float(rr)):
        return float(rr) * float(rr)
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=os.path.expanduser(
        "~/test-sargent/block1_1_v3"))
    ap.add_argument("--out", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "block1_1_v3_results.csv"))
    args = ap.parse_args()

    pred_files = sorted(glob.glob(os.path.join(args.base, "*_s*", "est.predict.json")))
    if not pred_files:
        print(f"No cells found under {args.base}", file=sys.stderr)
        sys.exit(1)

    rows = []
    skipped = 0
    for jf in pred_files:
        cell = os.path.basename(os.path.dirname(jf))
        m = CELL_PAT.match(cell)
        if not m:
            skipped += 1
            continue
        mtag, train_pop, test_pop = m.group(1), m.group(2), m.group(3)
        test_time = int(m.group(4)) if m.group(4) else 0
        mix, typed, seed = float(m.group(5)), float(m.group(6)), int(m.group(7))

        with open(jf) as fh:
            d = json.load(fh)

        pjf = os.path.join(os.path.dirname(jf), "case.params.json")
        params = json.load(open(pjf)) if os.path.exists(pjf) else {}

        rows.append({
            "cell":         cell,
            "model":        mtag,
            "train_pop":    train_pop,
            "test_pop":     test_pop,
            "test_time":    test_time,
            "fst":          params.get("fst_train_test"),
            "mix":          mix,
            "typed":        typed,
            "seed":         seed,
            "bayesR_r2c":   r2c(d.get("bayesR")),
            "gblup_r2c":    r2c(d.get("gblup")),
            "argunif_r2c":  r2c(d.get("arg_blup_uniform")),
            "sargent_r2c":  r2c(d.get("sargent")),
        })

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fields = ["cell", "model", "train_pop", "test_pop", "test_time",
              "fst", "mix", "typed", "seed",
              "bayesR_r2c", "gblup_r2c", "argunif_r2c", "sargent_r2c"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"Wrote {args.out}  ({len(rows)} cells; {skipped} skipped)")

    # Tiny sanity summary
    by_model = {}
    for r in rows:
        by_model.setdefault(r["model"], []).append(r)
    for mt, rs in sorted(by_model.items()):
        fsts = [r["fst"] for r in rs if isinstance(r["fst"], (int, float))]
        print(f"  {mt}: {len(rs)} cells; F_ST range "
              f"[{min(fsts):.3f}, {max(fsts):.3f}]" if fsts else f"  {mt}: {len(rs)} cells; F_ST n/a")


if __name__ == "__main__":
    main()
