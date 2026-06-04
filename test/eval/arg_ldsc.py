#!/usr/bin/env python
"""Block 2.2 flavour B — ARG-only partitioned heritability.

Reads SARGENT's per-tree posterior ws^2 (E_ws2) and the per-tree annotation
matrix M (fraction of each tree's span in each annotation), and attributes the
total ARG heritability across annotations.  This is the SARGENT-native primitive
that no other ARG-LMM paper has implemented.

Estimate (per-base variance-rate regression):
    fit  E[ws2_s] = sum_c rho_c L[s,c]   (L[s,c] = bp of annotation c in tree s)
    then h2_c^B = rho_c * (total bp of c), renormalised to SARGENT's h2_arg.

This lets an annotation be ENRICHED per base.  Span-FRACTION weighting (the
naive alternative) cannot: it returns ~ base_frac * h2_arg regardless of
enrichment, because a 150 bp exon gets ~3% of its kb-scale tree's ws2 even when
the exon variant is what elevated that tree.  NOTE: per-tree ws2 still COMPRESSES
fine-scale enrichment (a sub-tree annotation's signal smears across its tree),
so recovered enrichment underestimates the truth -- documented in the writeup.

Diagnostic: cond(L^T L) -- the structural-degeneracy check flagged in the design.
SE: delete-a-block jackknife over (contiguous) trees.
Output: <out>.B.json
"""
import argparse
import json
import numpy as np
import partition_lib as P


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-prefix", required=True)
    ap.add_argument("--trees", required=True, help="<case>.trees the est was fit on")
    ap.add_argument("--est-prefix", required=True, help="SARGENT --arg-only output prefix (est_B)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--jk-blocks", type=int, default=50)
    args = ap.parse_args()

    ann_ids = P.annotation_ids(args.case_prefix)
    w = P.read_w2(args.est_prefix)
    L = P.tree_annotation_bp(args.trees, args.case_prefix, ann_ids)
    if len(w) != L.shape[0]:
        raise SystemExit(f"tree count mismatch: w2 has {len(w)} trees, L has {L.shape[0]}")
    summ = P.read_summary(args.est_prefix)
    h2_arg = float(summ.get("h2_arg", summ.get("h2", 0.0)))

    h2c, rho, cond = P.regress_attribute(w, L, h2_arg)

    # jackknife over contiguous tree blocks
    n_tr = len(w)
    blocks = np.array_split(np.arange(n_tr), args.jk_blocks)
    jk = []
    for b in blocks:
        keep = np.ones(n_tr, bool); keep[b] = False
        hb, _, _ = P.regress_attribute(w[keep], L[keep], h2_arg)
        jk.append(hb)
    jk = np.array(jk); J = len(blocks)
    se = np.sqrt((J - 1) / J * np.sum((jk - jk.mean(0))**2, axis=0))

    bp_tot = L.sum(0)
    out = {"flavour": "B", "estimator": "ARG-only per-base ws2 rate regression",
           "h2_arg_total": h2_arg, "design_condition_number": cond,
           "n_trees": n_tr, "annotations": {}}
    for i, a in enumerate(ann_ids):
        out["annotations"][a] = {"h2c": float(h2c[i]), "se": float(se[i]),
                                 "rho_perbp": float(rho[i]),
                                 "ann_bp_frac": float(bp_tot[i] / bp_tot.sum())}
    json.dump(out, open(args.out + ".B.json", "w"), indent=2)
    print(f"[B] h2_arg={h2_arg:.3f} cond(MtM)={cond:.2e}  " +
          "  ".join(f"{a}:{out['annotations'][a]['h2c']:.4f}±{out['annotations'][a]['se']:.4f}"
                    for a in ann_ids))


if __name__ == "__main__":
    main()
