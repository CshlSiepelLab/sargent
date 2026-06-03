#!/usr/bin/env python3
"""
preprocess_fig4.py — emit tidy TSVs for fig4.R to consume.

Run once after sargent/predict outputs land; re-run only if the
underlying cells change. Outputs go to test/eval/figs_r/_cache/fig4_*.tsv.

Generated files:
  fig4_panelA.tsv    tree_idx, pos_mid_kb, ws2_est, ws2_truth, active
  fig4_panelB.tsv    arch, cell, kind (truth|est), k_frac, captured_frac
  fig4_panelC.tsv    contrast, demog, seed, auc, ceiling
  fig4_panelD.tsv    bf, mode (default|fused_lasso), seed, auc, fl_tag
"""

import argparse
import glob
import json
import os
import sys

import numpy as np
import tskit


# ---------- Panel A ----------

def panel_a(cell_dir, out_tsv):
    case = os.path.join(cell_dir, "case")
    est  = os.path.join(cell_dir, "est")
    if not os.path.exists(est + ".w2.tsv"):
        print(f"  [skip A] no {est}.w2.tsv", file=sys.stderr); return
    w2_est   = np.loadtxt(est + ".w2.tsv",          skiprows=1)[:, -1]
    truth_p  = case + ".true_ws2.tsv"
    if not os.path.exists(truth_p):
        truth_p = case + ".full.true_ws2.tsv"
    w2_truth = np.loadtxt(truth_p,                  skiprows=1)[:, -1]
    ts       = tskit.load(case + ".trees")
    breaks   = list(ts.breakpoints())
    mids_kb  = np.array([0.5 * (breaks[i] + breaks[i + 1]) for i in range(ts.num_trees)]) / 1e3

    n = min(len(w2_truth), len(w2_est), len(mids_kb))
    w2_truth = w2_truth[:n]; w2_est = w2_est[:n]; mids_kb = mids_kb[:n]
    active = (w2_truth > w2_truth.min() + 1e-12).astype(int)

    with open(out_tsv, "w") as f:
        f.write("tree_idx\tpos_mid_kb\tws2_est\tws2_truth\tactive\n")
        for i in range(n):
            f.write(f"{i}\t{mids_kb[i]:.6g}\t{w2_est[i]:.6g}\t{w2_truth[i]:.6g}\t{active[i]}\n")
    print(f"  wrote {out_tsv}  ({n} trees)")


# ---------- Panel B ----------

def _capture(w):
    w = np.sort(np.asarray(w, dtype=float))[::-1]
    cum = np.cumsum(w)
    tot = cum[-1] if cum[-1] > 0 else 1.0
    k = np.arange(1, len(w) + 1) / len(w)
    return k, cum / tot


def panel_b(base_phase10, out_tsv):
    rows = []
    for arch in ("sparse", "diffuse"):
        for cd in sorted(glob.glob(os.path.join(base_phase10, f"{arch}_s*"))):
            case = os.path.join(cd, "case")
            est  = os.path.join(cd, "est")
            if not os.path.exists(est + ".w2.tsv"): continue
            tp = case + ".full.true_ws2.tsv"
            if not os.path.exists(tp): tp = case + ".true_ws2.tsv"
            w_t = np.loadtxt(tp,                skiprows=1)[:, -1]
            w_e = np.loadtxt(est + ".w2.tsv",   skiprows=1)[:, -1]
            cell = os.path.basename(cd)
            for kind, w in [("truth", w_t), ("est", w_e)]:
                k, f = _capture(w)
                for ki, fi in zip(k, f):
                    rows.append((arch, cell, kind, float(ki), float(fi)))
    with open(out_tsv, "w") as f:
        f.write("arch\tcell\tkind\tk_frac\tcaptured_frac\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    print(f"  wrote {out_tsv}  ({len(rows)} rows)")


# ---------- Panel C ----------

def panel_c(base_b12, out_tsv):
    rows = []
    for jf in sorted(glob.glob(os.path.join(base_b12, "c*_d*_s*", "est.oracle.json"))):
        cell = os.path.basename(os.path.dirname(jf))
        # cell name like c1000_dadmixed_s1
        try:
            parts = cell.split("_")
            contrast = float(parts[0][1:])
            demog    = parts[1][1:]  # "admixed" or "panmictic"
            seed     = int(parts[2][1:])
        except (IndexError, ValueError):
            continue
        d = json.load(open(jf))
        auc  = d.get("ws2_recovery", {}).get("auc_active")
        ceil = d.get("loc_ceiling_auc")
        if not isinstance(auc, (int, float)): continue
        rows.append((contrast, demog, seed, float(auc),
                     float(ceil) if isinstance(ceil, (int, float)) else ""))
    with open(out_tsv, "w") as f:
        f.write("contrast\tdemog\tseed\tauc\tceiling\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    print(f"  wrote {out_tsv}  ({len(rows)} rows)")


# ---------- Panel D ----------

def _auc(scores, labels):
    pos = scores[labels == 1]; neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0: return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), float); ranks[order] = np.arange(1, len(scores) + 1)
    _, inv, counts = np.unique(scores, return_inverse=True, return_counts=True)
    tm = np.zeros(len(counts)); np.add.at(tm, inv, ranks); tm /= counts
    ranks = tm[inv]
    r_pos = ranks[labels == 1].sum()
    u = r_pos - len(pos) * (len(pos) + 1) / 2.0
    return float(u / (len(pos) * len(neg)))


def _truth_active(case_root, est_len):
    for suf in (".true_ws2.tsv", ".full.true_ws2.tsv"):
        p = os.path.join(case_root, "case" + suf)
        if not os.path.exists(p): continue
        wt = np.loadtxt(p, skiprows=1)[:, -1]
        if len(wt) == est_len:
            return wt, (wt > wt.min() + 1e-12).astype(int)
    return None, None


def panel_d(base_b252, fl_runs_dir, b26_real_dir, out_tsv):
    conditions = [("0.5", "lambda2=0.5"), ("0.005", "lambda1=0.5")]
    rows = []
    for bf, fl_tag in conditions:
        for seed in (1, 2, 3, 4, 5):
            case_root = os.path.join(base_b252, f"bf{bf}_dadmixed_s{seed}")
            unreg_p = os.path.join(case_root, "est.w2.tsv")
            if os.path.exists(unreg_p):
                w_u = np.loadtxt(unreg_p, skiprows=1)[:, -1]
                wt, active = _truth_active(case_root, len(w_u))
                if wt is not None and active.sum() > 0:
                    rows.append((bf, "default", seed, _auc(w_u, active), fl_tag))
            if seed == 1:
                fl_p = os.path.join(b26_real_dir,
                                    "bf05.l1_0.0_l2_0.5.w2.tsv" if bf == "0.5"
                                    else "bf005.l1_0.5_l2_0.0.w2.tsv")
            else:
                fl_p = os.path.join(fl_runs_dir, f"bf{bf}_fl_s{seed}.w2.tsv")
            if os.path.exists(fl_p):
                w_f = np.loadtxt(fl_p, skiprows=1)[:, -1]
                wt, active = _truth_active(case_root, len(w_f))
                if wt is not None and active.sum() > 0:
                    rows.append((bf, "fused_lasso", seed, _auc(w_f, active), fl_tag))
    with open(out_tsv, "w") as f:
        f.write("bf\tmode\tseed\tauc\tfl_tag\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    print(f"  wrote {out_tsv}  ({len(rows)} rows)")


# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache-dir", default=os.path.join(
        os.path.dirname(__file__), "_cache"))
    ap.add_argument("--cell-a", default=os.path.expanduser(
        "~/test-sargent/block1_2_v1/c1000_dadmixed_s1"))
    ap.add_argument("--base-phase10", default=os.path.expanduser(
        "~/test-sargent/phase10_admix"))
    ap.add_argument("--base-b12", default=os.path.expanduser(
        "~/test-sargent/block1_2_v1"))
    ap.add_argument("--base-b252", default=os.path.expanduser(
        "~/test-sargent/block2_5_2"))
    ap.add_argument("--fl-runs-dir", default="/tmp/fig4_fl_runs")
    ap.add_argument("--b26-real-dir", default="/tmp/b26_real")
    ap.add_argument("--panels", default="abcd",
                    help="Which panels to compute (subset of 'abcd').")
    args = ap.parse_args()

    os.makedirs(args.cache_dir, exist_ok=True)
    if "a" in args.panels: panel_a(args.cell_a,        os.path.join(args.cache_dir, "fig4_panelA.tsv"))
    if "b" in args.panels: panel_b(args.base_phase10,  os.path.join(args.cache_dir, "fig4_panelB.tsv"))
    if "c" in args.panels: panel_c(args.base_b12,      os.path.join(args.cache_dir, "fig4_panelC.tsv"))
    if "d" in args.panels: panel_d(args.base_b252, args.fl_runs_dir, args.b26_real_dir,
                                    os.path.join(args.cache_dir, "fig4_panelD.tsv"))


if __name__ == "__main__":
    main()
