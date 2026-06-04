#!/usr/bin/env python
"""Block 2.2 flavour C — hybrid (typed + ARG) partitioned heritability.

From the joint SARGENT fit (est_C: beta-hat for typed SNPs + per-tree ws^2),
produce a (typed, untyped) pair per annotation.

  C-typed_c     = sum_{j in c, typed} beta_hat_j^2          (typed-SNP variance)
  C-untyped_c   : TWO estimates, reported side by side
     (primary, robust)   = B_c - A_c        (the subtraction; see design doc)
     (cross-check)        = h2_arg(est_C) * frac_c(est_C ws^2) - C-typed_c
                            (the joint fit's own residual-ARG split)

If the two C-untyped estimates agree, the joint typed/untyped decomposition is
trustworthy; if they diverge, that is the honest finding (the internal split is
unstable -- same difficulty as the Block 1.1 h2_snp divergence), and the
subtraction is the number to trust.

Output: <out>.C.json
"""
import argparse
import json
import numpy as np
import partition_lib as P


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-prefix", required=True)
    ap.add_argument("--trees", required=True)
    ap.add_argument("--est-prefix", required=True, help="SARGENT full-fit prefix (est_C)")
    ap.add_argument("--a-json", required=True, help="flavour A output (.A.json)")
    ap.add_argument("--b-json", required=True, help="flavour B output (.B.json)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ann_ids = P.annotation_ids(args.case_prefix)
    params = P.read_params(args.case_prefix)
    n_bed = params["n_variants"]

    # C-typed: partition sum beta_hat^2 by annotation
    beta = P.read_betas(args.est_prefix, n_bed)
    col_ann = P.bed_column_annotations(args.case_prefix)
    if len(beta) != len(col_ann):
        # tolerate length differences by truncating to the common prefix
        m = min(len(beta), len(col_ann))
        beta, col_ann = beta[:m], col_ann[:m]
    c_typed = {a: float((beta[col_ann == a] ** 2).sum()) for a in ann_ids}

    # internal residual-ARG split from est_C ws^2 (same per-base regression as B)
    w = P.read_w2(args.est_prefix)
    L = P.tree_annotation_bp(args.trees, args.case_prefix, ann_ids)
    summ = P.read_summary(args.est_prefix)
    h2_arg_C = float(summ.get("h2_arg", 0.0))
    arg_h2, _, _ = P.regress_attribute(w, L, h2_arg_C)
    arg_by_ann = {a: float(arg_h2[i]) for i, a in enumerate(ann_ids)}

    A = json.load(open(args.a_json))["annotations"]
    B = json.load(open(args.b_json))["annotations"]

    out = {"flavour": "C", "estimator": "hybrid typed beta + ARG ws2",
           "h2_arg_C": h2_arg_C, "annotations": {}}
    for a in ann_ids:
        ct = c_typed[a]
        unt_sub = float(B[a]["h2c"] - A[a]["h2c"])              # primary
        unt_int = float(arg_by_ann[a])                          # cross-check (residual ARG)
        out["annotations"][a] = {
            "c_typed": ct,
            "c_untyped_subtraction": unt_sub,
            "c_untyped_internal": unt_int,
            "c_total_subtraction": ct + unt_sub,
        }
    json.dump(out, open(args.out + ".C.json", "w"), indent=2)
    print("[C] " + "  ".join(
        f"{a}: typ={out['annotations'][a]['c_typed']:.3f} "
        f"unt(B-A)={out['annotations'][a]['c_untyped_subtraction']:.3f} "
        f"unt(int)={out['annotations'][a]['c_untyped_internal']:.3f}"
        for a in ann_ids))


if __name__ == "__main__":
    main()
