#!/usr/bin/env python
"""Block 2.2 Phase C — A3a (allele-age) feasibility prototype (Stage 2).

For each Fig-5 Panel B per-edge-architecture cell, build per-age-tertile ARG
kinship matrices with the canonical assembler (test/eval/kinship_assembly.py),
fit a multi-component REML for the per-bin variance components (Stage-1 engine
in stratified_reml.py), and test whether ARG-REML recovers the per-bin h²
diagonal that the Phase-B per-edge posterior could not.

Pass: for each true architecture, the matching observed-bin sigma2_hat is the
largest of the three by >= 0.10 (mean over 3 seeds).
Bonus: per-cell total h_hat^2 = Vg/Vp >= 0.5 (vs simulated 0.85).

Usage: python stratified_reml_a3a.py [--base DIR] [--out TSV]
"""
import argparse
import csv
import glob
import os

import numpy as np
import tskit

import kinship_assembly as ka
import stratified_reml as srm

BINS = ["recent", "middle", "ancient"]


def truth_per_bin(cell, q33, q67):
    """Truth-side per-bin sigma2 = Σ_{e∈bin} ws2_e·ℓ_e from
    case.true_ws2_per_edge.tsv (cols: tree_idx child_node parent_age child_age
    length ws2_e), binned on the same age tertiles as the estimator."""
    path = os.path.join(cell, "case.true_ws2_per_edge.tsv")
    out = {b: 0.0 for b in BINS}
    if not os.path.exists(path):
        return out
    d = np.loadtxt(path, skiprows=1, ndmin=2)
    mids = 0.5 * (d[:, 2] + d[:, 3])
    ell, ws2 = d[:, 4], d[:, 5]
    for m, l, w in zip(mids, ell, ws2):
        b = "recent" if m <= q33 else ("middle" if m <= q67 else "ancient")
        out[b] += w * l
    return out


def run_cell(cell):
    ts = tskit.load(os.path.join(cell, "case.trees"))
    n = ts.num_individuals
    y = np.loadtxt(os.path.join(cell, "case.phen"))
    y = (y - y.mean()) / (y.std() or 1.0)          # unit variance -> σ̂²_b read as h²_b
    q33, q67 = ka._per_cell_tertiles(ts)

    Kb = ka.build_per_bin_kinship(ts, ka.bin_by_age_tertile(ts), n)
    Ks = [Kb[b] for b in BINS]

    # sum-rule guardrail: Σ_b K_b must equal K_ARG
    Ksum = sum(Ks)
    Karg = ka._reference_k_arg(ts, n)
    rel_err = float(np.linalg.norm(Ksum - Karg) / (np.linalg.norm(Karg) + 1e-12))

    # GREML-MS convention: normalise EACH K_b to mean diagonal 1, so σ̂²_b is the
    # variance explained by bin b and (y standardised) σ̂²_b ≈ h²_b, Σ_b σ̂²_b ≈ total h².
    Ksn = [K / (float(np.mean(np.diag(K))) or 1.0) for K in Ks]

    res = srm.reml_mc(y, Ksn, em_warmup=3, max_iter=200)
    truth = truth_per_bin(cell, q33, q67)
    return res, truth, rel_err, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=os.path.expanduser("~/test-sargent/fig5_panelB_edges"))
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                  "figs_r/_cache/fig5_panelB_A3a.tsv"))
    args = ap.parse_args()
    cells = sorted(glob.glob(os.path.join(args.base, "*_s*")))
    if not cells:
        raise SystemExit(f"no cells under {args.base}")

    rows = []
    print(f"A3a — per-age-bin ARG-REML on {len(cells)} cells\n")
    print(f"{'cell':>24} {'sumrule':>8} {'sig_rec':>8} {'sig_mid':>8} {'sig_anc':>8} "
          f"{'sig_e':>7} {'tot_h2':>7} {'conv':>5}")
    for cell in cells:
        res, truth, rel_err, n = run_cell(cell)
        name = os.path.basename(cell)
        arch = name.rsplit("_s", 1)[0]
        seed = int(name.rsplit("_s", 1)[1])
        s = res["sigma2"]
        print(f"{name:>24} {rel_err:>8.1e} {s[0]:>8.3f} {s[1]:>8.3f} {s[2]:>8.3f} "
              f"{res['sigma2_e']:>7.3f} {res['total_h2']:>7.3f} {str(res['converged']):>5}")
        for i, b in enumerate(BINS):
            rows.append({"true_arch": arch, "observed_bin": b, "seed": seed,
                         "sigma2_hat": s[i], "se": res["se"][i],
                         "sigma2_e": res["sigma2_e"], "total_h2": res["total_h2"],
                         "truth_per_bin": truth[b], "sumrule_relerr": rel_err})

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["true_arch", "observed_bin", "seed",
                                          "sigma2_hat", "se", "sigma2_e", "total_h2",
                                          "truth_per_bin", "sumrule_relerr"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {args.out} ({len(rows)} rows)")

    # 3x3 mean sigma2_hat matrix (arch x observed bin) + pass check
    archs = ["age_recent_edges", "age_middle_edges", "age_ancient_edges"]
    print("\nMean sigma2_hat by (true_arch, observed_bin)  [diagonal should dominate]:")
    print(f"  {'true_arch':>20} {'recent':>9} {'middle':>9} {'ancient':>9}  {'diag-2nd':>9} {'PASS':>5}")
    overall = True
    for ai, a in enumerate(archs):
        means = {b: np.mean([r["sigma2_hat"] for r in rows
                             if r["true_arch"] == a and r["observed_bin"] == b]) for b in BINS}
        diagb = BINS[ai]
        diag = means[diagb]
        second = max(v for b, v in means.items() if b != diagb)
        ok = (diag - second) >= 0.10
        overall = overall and ok
        print(f"  {a:>20} {means['recent']:>9.3f} {means['middle']:>9.3f} "
              f"{means['ancient']:>9.3f}  {diag-second:>+9.3f} {str(ok):>5}")
    tot = np.mean([r["total_h2"] for r in rows]) if rows else 0.0
    print(f"\n  diagonal PASS (each diag - 2nd >= 0.10): {overall}")
    print(f"  mean total_h2 = {tot:.3f}  (bonus PASS if >= 0.5; simulated 0.85)")
    maxrel = max(r["sumrule_relerr"] for r in rows)
    print(f"  max sum-rule rel err = {maxrel:.1e}  (Σ_b K_b vs K_ARG; should be ~0)")


if __name__ == "__main__":
    main()
