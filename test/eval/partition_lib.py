#!/usr/bin/env python
"""Shared helpers for Block 2.2 partitioned-heritability estimators.

Reads the simulate.py --annotation-genes outputs (case.annotations.tsv,
case.annot_map.tsv, case.params.json) and SARGENT outputs (est.w2.tsv,
est.summary.json, est.beta.tsv), and provides the per-tree annotation-fraction
matrix that flavours B and C regress against.
"""
import json
import numpy as np


# ---- simulate.py annotation outputs ------------------------------------------

def read_params(case_prefix):
    return json.load(open(case_prefix + ".params.json"))


def annotation_ids(case_prefix):
    p = read_params(case_prefix)
    ids = p.get("annotation_ids")
    if ids is None:
        raise ValueError(f"{case_prefix}.params.json has no annotation_ids "
                         "(was it simulated with --annotation-genes?)")
    return list(ids)


def read_annotations(case_prefix):
    """Per full-panel variant table.

    Returns dict of np arrays: snp_idx, bp_pos, annotation (str), in_bed (bool),
    is_causal (bool), all length p (full panel, snp_idx order).
    """
    snp_idx, bp, ann, inbed, caus = [], [], [], [], []
    with open(case_prefix + ".annotations.tsv") as f:
        next(f)
        for line in f:
            a, b, c, d, e = line.rstrip("\n").split("\t")
            snp_idx.append(int(a)); bp.append(int(b)); ann.append(c)
            inbed.append(d == "1"); caus.append(e == "1")
    return {"snp_idx": np.array(snp_idx), "bp_pos": np.array(bp),
            "annotation": np.array(ann), "in_bed": np.array(inbed),
            "is_causal": np.array(caus)}


def read_annot_map(case_prefix):
    """The contiguous genomic partition: arrays start, end, ann (str), sorted."""
    s, e, a = [], [], []
    with open(case_prefix + ".annot_map.tsv") as f:
        next(f)
        for line in f:
            x, y, z = line.rstrip("\n").split("\t")
            s.append(int(x)); e.append(int(y)); a.append(z)
    return np.array(s, float), np.array(e, float), np.array(a)


def bed_column_annotations(case_prefix):
    """Annotation label of each BED (typed) column, in BED column order.

    BED columns are the in_bed variants in snp_idx order (simulate.py
    write_plink convention), so this is just the annotation of in_bed rows.
    """
    t = read_annotations(case_prefix)
    return t["annotation"][t["in_bed"]]


# ---- per-tree annotation-fraction matrix (flavours B / C) ---------------------

def tree_annotation_bp(trees_file, case_prefix, ann_ids):
    """L (n_trees x n_ann): base pairs of each annotation in each tree's span.

    Row s sums to tree s's span.  Tree order matches est.w2.tsv (tskit order).
    This bp matrix (not fractions) is what flavours B/C regress ws^2 against:
    a per-base variance-rate model E[ws2_s] = sum_c rho_c L[s,c] lets an
    annotation be ENRICHED per base, which fraction-weighting cannot detect.
    """
    import tskit
    ts = tskit.load(trees_file)
    ms, me, ma = read_annot_map(case_prefix)
    ann_index = {a: i for i, a in enumerate(ann_ids)}
    ma_idx = np.array([ann_index[a] for a in ma])
    n_ann = len(ann_ids)
    L = np.zeros((ts.num_trees, n_ann))
    for s, tr in enumerate(ts.trees()):
        l, r = tr.interval.left, tr.interval.right
        if r <= l:
            continue
        lo = np.searchsorted(me, l, side="right")
        hi = np.searchsorted(ms, r, side="left")
        for k in range(lo, hi):
            ov = min(r, me[k]) - max(l, ms[k])
            if ov > 0:
                L[s, ma_idx[k]] += ov
    return L


def regress_attribute(w, L, total):
    """Per-base variance-rate attribution of ws^2 across annotations.

    Fit E[w_s] = sum_c rho_c L[s,c] (non-negative LS), then annotation c's
    heritability = rho_c * (total bp of c), renormalised to `total` (SARGENT's
    h2_arg).  Returns (h2 per annotation, rho per-bp rates, cond(L^T L)).
    """
    rho, *_ = np.linalg.lstsq(L, w, rcond=None)
    rho = np.clip(rho, 0.0, None)
    h2 = rho * L.sum(0)
    s = h2.sum()
    h2 = (h2 / s * total) if s > 0 else h2
    cond = float(np.linalg.cond(L.T @ L))
    return h2, rho, cond


# ---- SARGENT outputs ---------------------------------------------------------

def read_w2(est_prefix):
    """Per-tree posterior E_ws2 (np array, tree order)."""
    w = []
    with open(est_prefix + ".w2.tsv") as f:
        next(f)
        for line in f:
            w.append(float(line.rstrip("\n").split("\t")[-1]))   # E_ws2 col
    return np.array(w)


def read_summary(est_prefix):
    return json.load(open(est_prefix + ".summary.json"))


def read_betas(est_prefix, n_bed):
    """Per-typed-SNP fitted beta (np array length n_bed, BED column order).

    SARGENT writes est.beta.tsv; if absent (e.g. --arg-only) returns zeros.
    """
    import os
    path = est_prefix + ".beta.tsv"
    if not os.path.exists(path):
        return np.zeros(n_bed)
    vals = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        # last numeric column is the (posterior) beta; tolerate idx/beta layouts
        bcol = len(header) - 1
        for line in f:
            parts = line.rstrip("\n").split("\t")
            vals.append(float(parts[bcol]))
    v = np.array(vals)
    if len(v) != n_bed:
        # if an index column made it 2-wide etc., fall back to last col already used
        pass
    return v


# ---- phenotype + typed genotypes (flavour A) ---------------------------------

def read_phen(case_prefix):
    y = np.loadtxt(case_prefix + ".phen")
    return (y - y.mean()) / (y.std() if y.std() > 0 else 1.0)


def typed_genotypes_std(case_prefix):
    """Standardised typed-SNP dosage (n_bed x n) and per-column annotation.

    Reuses predict.py:read_bed.  Returns (Z, col_ann) with Z standardised to
    zero mean / unit variance per SNP (monomorphic columns dropped).
    """
    import importlib.util, os
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location("predict", os.path.join(here, "predict.py"))
    # predict.py imports heavy deps at module load; only need read_bed, so exec guarded
    pred = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pred)
    p = read_params(case_prefix)
    n = p["n_indiv"]; n_bed = p["n_variants"]
    dos = pred.read_bed(case_prefix, n, n_bed)          # (n_bed, n) raw 0/1/2
    col_ann = bed_column_annotations(case_prefix)
    mu = np.nanmean(dos, axis=1, keepdims=True)
    sd = np.nanstd(dos, axis=1, keepdims=True)
    keep = (sd[:, 0] > 0)
    Z = (dos[keep] - mu[keep]) / sd[keep]
    Z = np.nan_to_num(Z)
    return Z, col_ann[keep]
