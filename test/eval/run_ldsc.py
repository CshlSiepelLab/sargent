#!/usr/bin/env python
"""Block 2.2 flavour A — typed-SNP-only partitioned heritability.

PROXY for stratified LD-score regression (Finucane et al. 2015).  The real
S-LDSC tool is not installed in this environment, and at the n~800 in-sample
scale of the v1 reference simulation S-LDSC would be badly underpowered (it is
designed for N>>10k GWAS summary stats + an external LD reference), so a direct
run would penalise it for sample size rather than for its actual rare-variant
blind spot.  Instead we use the heritability-partition method that is principled
at this scale: multi-component Haseman-Elston (method-of-moments GREML) on
per-annotation GRMs built from TYPED SNPs only.  This is exactly "what a
typed-SNP method can recover", which is what flavour A is meant to represent.

Stage 2 will install real Finucane S-LDSC and run it on a scaled-up,
properly-powered simulation; the function `sldsc_stub()` below documents that
interface.

Method:
  y standardised; for each annotation c build K_c = Z_c Z_c^T / m_c from typed
  SNPs in c (Z standardised).  Regress off-diagonal products y_i y_j on
  {K_c[i,j]} jointly (no intercept) -> per-annotation h2_c.  SE by delete-a-block
  jackknife over individuals.

Output: JSON {annotation: {h2c, se}}  -> <out>.A.json
"""
import argparse
import json
import sys
import numpy as np
import partition_lib as P


def sldsc_stub(*_a, **_k):
    """Stage-2 placeholder for real Finucane S-LDSC.

    Interface (to implement against `ldsc.py` once installed):
      inputs : per-typed-SNP GWAS summary stats (.sumstats), per-annotation
               LD scores (.l2.ldscore) from a matched LD reference panel,
               regression weights.
      call   : ldsc.py --h2 sumstats --ref-ld annot --w-ld weights --overlap-annot
      output : per-annotation h2_c + jackknife SE from the .results file.
    """
    raise NotImplementedError("real S-LDSC deferred to Stage 2 (see module docstring)")


def he_partition(y, K_list):
    """Multi-component Haseman-Elston: solve E[y_i y_j] = sum_c h2_c K_c[i,j].

    Returns h2 (len C) using off-diagonal pairs, plus the design condition no.
    """
    n = len(y)
    iu = np.triu_indices(n, k=1)
    resp = (np.outer(y, y))[iu]                       # y_i y_j, i<j
    X = np.column_stack([K[iu] for K in K_list])      # (npairs, C)
    XtX = X.T @ X
    beta = np.linalg.solve(XtX, X.T @ resp)
    cond = np.linalg.cond(XtX)
    return beta, cond


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case-prefix", required=True)
    ap.add_argument("--out", required=True, help="output prefix (writes <out>.A.json)")
    ap.add_argument("--jk-blocks", type=int, default=20, help="jackknife blocks over individuals")
    args = ap.parse_args()

    ann_ids = P.annotation_ids(args.case_prefix)
    y = P.read_phen(args.case_prefix)
    Z, col_ann = P.typed_genotypes_std(args.case_prefix)
    n = len(y)

    # Per-annotation GRM from typed SNPs in that annotation.
    present = [a for a in ann_ids if (col_ann == a).sum() > 0]
    K_list, m_c = [], {}
    for a in present:
        Za = Z[col_ann == a]
        m_c[a] = Za.shape[0]
        K_list.append((Za.T @ Za) / Za.shape[0])      # (n, n)

    h2, cond = he_partition(y, K_list)

    # delete-a-block jackknife over individuals
    rng = np.random.default_rng(0)
    perm = rng.permutation(n)
    blocks = np.array_split(perm, args.jk_blocks)
    jk = []
    for b in blocks:
        keep = np.ones(n, bool); keep[b] = False
        yk = y[keep]
        Kk = [K[np.ix_(keep, keep)] for K in K_list]
        hb, _ = he_partition(yk, Kk)
        jk.append(hb)
    jk = np.array(jk)
    J = len(blocks)
    se = np.sqrt((J - 1) / J * np.sum((jk - jk.mean(0))**2, axis=0))

    out = {"flavour": "A", "estimator": "typed-only HE/GREML proxy (S-LDSC stand-in)",
           "design_condition_number": float(cond), "n": n,
           "annotations": {}}
    for i, a in enumerate(present):
        out["annotations"][a] = {"h2c": float(h2[i]), "se": float(se[i]),
                                 "m_typed": int(m_c[a])}
    for a in ann_ids:
        out["annotations"].setdefault(a, {"h2c": 0.0, "se": 0.0, "m_typed": 0})
    json.dump(out, open(args.out + ".A.json", "w"), indent=2)
    print(f"[A] cond(XtX)={cond:.2e}  " +
          "  ".join(f"{a}:{out['annotations'][a]['h2c']:.4f}±{out['annotations'][a]['se']:.4f}"
                    for a in ann_ids))


if __name__ == "__main__":
    main()
