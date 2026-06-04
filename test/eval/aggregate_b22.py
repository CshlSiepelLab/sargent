#!/usr/bin/env python
"""Block 2.2 Phase B3 — aggregate the three flavours into the Fig 5 Panel A CSV.

One row per (annotation x estimator_flavour) with columns:
  annotation, base_frac, true_tau, true_h2, est_h2c, est_se, estimator_flavour
where estimator_flavour in {A, B, C-typed, C-untyped} and true_h2 is the truth
that flavour targets (typed for A / C-typed, total for B, untyped for C-untyped).
"""
import argparse
import csv
import json
import math
import numpy as np
import partition_lib as P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case-prefix", required=True)
    ap.add_argument("--a-json", required=True)
    ap.add_argument("--b-json", required=True)
    ap.add_argument("--c-json", required=True)
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args()

    ann_ids = P.annotation_ids(args.case_prefix)
    truth = P.read_params(args.case_prefix)["h2_by_annotation"]
    A = json.load(open(args.a_json))["annotations"]
    B = json.load(open(args.b_json))["annotations"]
    C = json.load(open(args.c_json))["annotations"]

    rows = []
    for a in ann_ids:
        t = truth[a]
        base, tau = t["base_frac"], t["tau_c"]
        # A — targets typed truth
        rows.append([a, base, tau, t["h2_typed"], A[a]["h2c"], A[a]["se"], "A"])
        # B — targets total truth
        rows.append([a, base, tau, t["h2_total"], B[a]["h2c"], B[a]["se"], "B"])
        # C-typed — targets typed truth
        rows.append([a, base, tau, t["h2_typed"], C[a]["c_typed"], float("nan"), "C-typed"])
        # C-untyped — targets untyped truth; primary = subtraction, se propagated
        se_unt = math.sqrt(A[a]["se"]**2 + B[a]["se"]**2)
        rows.append([a, base, tau, t["h2_untyped"], C[a]["c_untyped_subtraction"],
                     se_unt, "C-untyped"])

    with open(args.out_csv, "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["annotation", "base_frac", "true_tau", "true_h2",
                     "est_h2c", "est_se", "estimator_flavour"])
        for r in rows:
            wr.writerow([r[0]] + [f"{x:.6g}" if isinstance(x, float) else x for x in r[1:]])
    print(f"wrote {len(rows)} rows -> {args.out_csv}")
    # quick console view
    print(f"\n{'anno':<11}{'flavour':<11}{'true_h2':>9}{'est_h2c':>9}{'est_se':>8}")
    for r in rows:
        se = "" if (isinstance(r[5], float) and np.isnan(r[5])) else f"{r[5]:.4f}"
        print(f"{r[0]:<11}{r[6]:<11}{r[3]:>9.4f}{r[4]:>9.4f}{se:>8}")


if __name__ == "__main__":
    main()
