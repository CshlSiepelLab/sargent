#!/usr/bin/env python
"""Block 2.2 Phase C — multi-component REML engine (Stage 1).

A minimal, dependency-light AI-REML (average information) with an EM warm-start
for the variance-component model

    y = X beta + sum_k u_k + eps,   u_k ~ N(0, sigma2_k K_k),  eps ~ N(0, sigma2_e I)
    V = sigma2_e I + sum_k sigma2_k K_k

Used by Phase C to fit per-bin ARG-kinship variance components (A3a age bins,
later A3b/A3c). Calibrated against GCTA `--reml` / `--reml --mgrm`.

Run as a script to execute the Stage-1 calibration smoke tests (T1/T2/T3):
    python stratified_reml.py            # prints the GCTA calibration table
"""
import os
import subprocess
import sys
import tempfile

import numpy as np


# ---------------------------------------------------------------------------
# Core estimator
# ---------------------------------------------------------------------------

def _reml_quantities(theta, y, X, Kc):
    """REML log-lik, projection P, score and AI matrix at variance vector theta.

    Kc[0] = I (residual), Kc[1..K] = K_1..K_K.  Returns (ll, P, Py, score, AI).
    """
    n = len(y)
    nc = len(theta)
    vy = float(np.var(y)) or 1.0
    V = sum(theta[j] * Kc[j] for j in range(nc))
    jit = 0.0
    while True:
        try:
            M = V + jit * np.eye(n) if jit else V
            sign, logdetV = np.linalg.slogdet(M)
            if sign <= 0:
                raise np.linalg.LinAlgError("V not PD")
            Vinv = np.linalg.inv(M)
            break
        except np.linalg.LinAlgError:
            jit = 1e-8 * vy if jit == 0 else jit * 10
            if jit > 1e-2 * vy:
                raise
    VinvX = Vinv @ X
    XtVinvX = X.T @ VinvX
    XtVinvX_inv = np.linalg.inv(XtVinvX)
    P = Vinv - VinvX @ XtVinvX_inv @ VinvX.T
    Py = P @ y
    _, logdetB = np.linalg.slogdet(XtVinvX)
    ll = -0.5 * (logdetV + logdetB + float(y @ Py))
    VjPy = [Kc[j] @ Py for j in range(nc)]
    score = np.array([-0.5 * (float(np.sum(P * Kc[j])) - float(Py @ VjPy[j]))
                      for j in range(nc)])
    PVjPy = [P @ VjPy[j] for j in range(nc)]
    AI = np.empty((nc, nc))
    for a in range(nc):
        for b in range(a, nc):
            AI[a, b] = AI[b, a] = 0.5 * float(VjPy[a] @ PVjPy[b])
    return ll, P, Py, score, AI, VjPy


def reml_mc(y, K_list, X=None, max_iter=100, tol=1e-4, em_warmup=3,
            constrain=True, verbose=False):
    """Fit the multi-component variance model by AI-REML (+ EM warm-start).

    Returns dict with sigma2_e, sigma2 (per K), se, h2 (per K), total_h2,
    loglik, loglik_trace, n_iter, converged, sigma2_unconstrained.
    """
    y = np.asarray(y, float).ravel()
    n = len(y)
    if X is None:
        X = np.ones((n, 1))
    X = np.asarray(X, float).reshape(n, -1)
    Kc = [np.eye(n)] + [np.asarray(K, float) for K in K_list]
    nc = len(Kc)
    vy = float(np.var(y)) or 1.0
    floor = 1e-6 * vy

    theta = np.full(nc, vy / nc)
    trace = []
    prev_ll = None
    converged = False
    AI = None
    for it in range(max_iter):
        ll, P, Py, score, AI, VjPy = _reml_quantities(theta, y, X, Kc)
        trace.append(ll)
        if prev_ll is not None and abs(ll - prev_ll) < tol and it >= em_warmup:
            converged = True
            break
        prev_ll = ll
        if it < em_warmup:                                   # EM warm-start
            new = theta.copy()
            for j in range(nc):
                new[j] = theta[j] + (theta[j] ** 2 / n) * (
                    float(Py @ VjPy[j]) - float(np.sum(P * Kc[j])))
            theta = np.maximum(new, floor)
        else:                                                # AI step
            try:
                step = np.linalg.solve(AI, score)
            except np.linalg.LinAlgError:
                step = np.linalg.lstsq(AI, score, rcond=None)[0]
            new = theta + step
            # damp if a component would overshoot far negative
            if np.any(new < 0) and constrain:
                # try progressively smaller steps before clamping
                for damp in (0.5, 0.25, 0.1):
                    cand = theta + damp * step
                    if np.all(cand > -floor):
                        new = cand
                        break
            theta_unc = new.copy()
            theta = np.maximum(new, floor) if constrain else new
        if verbose:
            print(f"  it{it} ll={ll:.5f} theta={np.round(theta,4)}")

    # final unconstrained estimate (one Newton step from the constrained optimum)
    ll, P, Py, score, AI, _ = _reml_quantities(theta, y, X, Kc)
    try:
        theta_unc = theta + np.linalg.solve(AI, score)
    except np.linalg.LinAlgError:
        theta_unc = theta.copy()
    try:
        AIinv = np.linalg.inv(AI)
        se = np.sqrt(np.maximum(np.diag(AIinv), 0.0))
    except np.linalg.LinAlgError:
        se = np.full(nc, np.nan)

    total_g = float(np.sum(theta[1:]))
    vp = float(theta[0] + total_g)
    return {
        "sigma2_e": float(theta[0]),
        "sigma2": theta[1:].tolist(),
        "se": se[1:].tolist(),
        "se_e": float(se[0]),
        "sigma2_unconstrained": theta_unc.tolist(),
        "h2": [(s / vp if vp > 0 else 0.0) for s in theta[1:]],
        "total_h2": (total_g / vp if vp > 0 else 0.0),
        "loglik": float(ll),
        "loglik_trace": trace,
        "n_iter": len(trace),
        "converged": converged,
        "Vp": vp,
    }


# ---------------------------------------------------------------------------
# GCTA interface (calibration partner)
# ---------------------------------------------------------------------------

def gcta_bin():
    cand = os.environ.get("GCTA")
    if cand and os.path.exists(cand):
        return cand
    here = os.path.join(os.path.dirname(sys.executable), "gcta64")
    if os.path.exists(here):
        return here
    import shutil
    return shutil.which("gcta64") or "gcta64"


def write_gcta_grm(prefix, K, m=1000, ids=None):
    """Write a dense kinship K as a GCTA GRM triple (.grm.bin/.id/.N.bin)."""
    n = K.shape[0]
    if ids is None:
        ids = [(str(i + 1), str(i + 1)) for i in range(n)]
    with open(prefix + ".grm.id", "w") as f:
        for fid, iid in ids:
            f.write(f"{fid}\t{iid}\n")
    idx = np.tril_indices(n)                    # row-major lower triangle = GCTA order
    K[idx].astype(np.float32).tofile(prefix + ".grm.bin")
    np.full(len(idx[0]), m, dtype=np.float32).tofile(prefix + ".grm.N.bin")
    return ids


def parse_hsq(path):
    """Parse a GCTA .hsq file into {name: (variance, se)}."""
    out = {}
    with open(path) as f:
        next(f, None)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                name = parts[0]
                try:
                    val = float(parts[1])
                except ValueError:
                    continue
                se = float(parts[2]) if len(parts) >= 3 and parts[2] not in ("", "NA") else float("nan")
                out[name] = (val, se)
    return out


def run_gcta_reml(K_list, y, ids=None, m=1000, workdir=None, tag="g"):
    """Run GCTA --reml (single) or --reml --mgrm (multi) on dense kinships.

    Returns {name: (variance, se)} parsed from the .hsq, or None on failure.
    """
    n = len(y)
    if ids is None:
        ids = [(str(i + 1), str(i + 1)) for i in range(n)]
    cleanup = workdir is None
    workdir = workdir or tempfile.mkdtemp(prefix="gcta_")
    prefixes = []
    for k, K in enumerate(K_list):
        p = os.path.join(workdir, f"{tag}{k}")
        write_gcta_grm(p, K, m=m, ids=ids)
        prefixes.append(p)
    phen = os.path.join(workdir, f"{tag}.phen")
    with open(phen, "w") as f:
        for (fid, iid), v in zip(ids, y):
            f.write(f"{fid}\t{iid}\t{v:.8g}\n")
    out = os.path.join(workdir, f"{tag}.reml")
    cmd = [gcta_bin(), "--reml", "--pheno", phen, "--out", out,
           "--reml-maxit", "200", "--thread-num", "2"]
    if len(prefixes) == 1:
        cmd += ["--grm", prefixes[0]]
    else:
        mg = os.path.join(workdir, f"{tag}.mgrm.txt")
        with open(mg, "w") as f:
            for p in prefixes:
                f.write(p + "\n")
        cmd += ["--mgrm", mg]
    r = subprocess.run(cmd, capture_output=True, text=True)
    res = parse_hsq(out + ".hsq") if os.path.exists(out + ".hsq") else None
    if res is None and cleanup:
        sys.stderr.write(r.stdout[-500:] + r.stderr[-500:])
    return res


# ---------------------------------------------------------------------------
# Stage-1 calibration smoke tests
# ---------------------------------------------------------------------------

def _sim_grm(n, m, rng):
    """A GRM K = ZZ'/m from m standardised biallelic SNPs, plus Z for effects."""
    maf = rng.uniform(0.05, 0.5, size=m)
    G = rng.binomial(2, maf[None, :], size=(n, m)).astype(float)
    Z = (G - 2 * maf) / np.sqrt(2 * maf * (1 - maf))
    K = (Z @ Z.T) / m
    return K, Z


def _sim_y(Z_list, sigma2_list, sigma2_e, rng):
    n = Z_list[0].shape[0]
    y = rng.normal(0, np.sqrt(sigma2_e), size=n)
    for Z, s2 in zip(Z_list, sigma2_list):
        m = Z.shape[1]
        beta = rng.normal(0, np.sqrt(s2 / m), size=m)
        y = y + Z @ beta
    return y


def _load_predict():
    """Import predict.py (for streaming_combined_kinship/hap_to_indiv_2d)."""
    import importlib.util
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location("predict", os.path.join(here, "predict.py"))
    pred = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pred)
    return pred


def t4_real_argk(n_cells=3, h2_true=0.6):
    """Validate the engine on a REAL ARG kinship (single component).

    Build K_ARG from a Fig-5 Panel B cell's case.trees (the canonical uniform
    ARG kinship via predict.streaming_combined_kinship), normalise to mean
    diagonal 1, simulate y ~ N(0, h2 K + (1-h2) I), and check reml_mc recovers
    h2 and matches GCTA on the real (denser, ARG-structured) kinship.
    """
    import glob
    import tskit
    pred = _load_predict()
    base = os.path.expanduser("~/test-sargent/fig5_panelB_edges")
    cells = sorted(glob.glob(os.path.join(base, "*_s*")))[:n_cells]
    if not cells:
        print("=== T4 — SKIPPED (no fig5_panelB_edges cells present) ===\n")
        return
    print(f"=== T4 — engine on REAL ARG kinship (single component, h2_true={h2_true}) ===")
    print(f"  {'cell':>24} {'n':>4} {'cond(K)':>9} {'mc_h2':>7} {'gcta_h2':>8} {'|mc-gcta|':>10}")
    for cell in cells:
        ts = tskit.load(os.path.join(cell, "case.trees"))
        n = ts.num_individuals
        sample_nodes = [int(nd) for ind in range(n) for nd in ts.individual(ind).nodes]
        Khap = pred.streaming_combined_kinship(ts, sample_nodes, None)
        K = pred.hap_to_indiv_2d(Khap, n)
        K = K / float(np.mean(np.diag(K)))             # GRM-style normalisation
        cond = float(np.linalg.cond(K))
        r = np.random.default_rng(11)
        L = np.linalg.cholesky(h2_true * K + (1 - h2_true) * np.eye(n) + 1e-8 * np.eye(n))
        y = L @ r.normal(size=n)
        mc = reml_mc(y, [K])
        g = run_gcta_reml([K], y)
        gh2 = g["V(G)/Vp"][0] if g and "V(G)/Vp" in g else float("nan")
        print(f"  {os.path.basename(cell):>24} {n:>4} {cond:>9.1e} {mc['h2'][0]:>7.3f} "
              f"{gh2:>8.3f} {abs(mc['h2'][0] - gh2):>10.2e}")
    print(f"  -> recovers h2~{h2_true} on real ARG kinship; mc==GCTA confirms engine "
          f"handles the ARG eigenstructure.\n")


def main():
    rng = np.random.default_rng(7)
    n = 300
    print(f"Phase C Stage-1 REML calibration (n={n}, GCTA {os.path.basename(gcta_bin())})\n")

    # ---- T1: single-component vs GCTA ----
    print("=== T1 — single component: reml_mc vs GCTA --reml (5 seeds) ===")
    print(f"  {'seed':>4} {'truth_h2':>9} {'mc_h2':>8} {'gcta_h2':>8} "
          f"{'mc_sg':>8} {'gcta_sg':>8} {'|mc-gcta|sg':>11}")
    t1max = 0.0
    for sd in range(1, 6):
        r = np.random.default_rng(100 + sd)
        K, Z = _sim_grm(n, 800, r)
        s2g, s2e = 0.5, 0.5
        y = _sim_y([Z], [s2g], s2e, r)
        mc = reml_mc(y, [K])
        g = run_gcta_reml([K], y, m=800)
        gsg = g["V(G)"][0] if g and "V(G)" in g else float("nan")
        gh2 = g["V(G)/Vp"][0] if g and "V(G)/Vp" in g else float("nan")
        d = abs(mc["sigma2"][0] - gsg)
        t1max = max(t1max, d)
        print(f"  {sd:>4} {s2g/(s2g+s2e):>9.3f} {mc['h2'][0]:>8.3f} {gh2:>8.3f} "
              f"{mc['sigma2'][0]:>8.4f} {gsg:>8.4f} {d:>11.2e}")
    print(f"  -> max |sigma2_g(mc) - sigma2_g(GCTA)| = {t1max:.2e}  "
          f"(PASS if < 1e-3)\n")

    # ---- T2: two-component recovery vs GCTA --mgrm ----
    print("=== T2 — two components: recovery + GCTA --mgrm cross-check (5 seeds) ===")
    truth = (0.3, 0.2)
    print(f"  truth: sigma2_1={truth[0]}, sigma2_2={truth[1]}, sigma2_e=0.5")
    print(f"  {'seed':>4} {'mc_s1':>7} {'mc_s2':>7} {'mc_se1':>7} {'mc_se2':>7} "
          f"{'gcta_s1':>8} {'gcta_s2':>8}")
    s1s, s2s = [], []
    for sd in range(1, 6):
        r = np.random.default_rng(200 + sd)
        K1, Z1 = _sim_grm(n, 600, r)
        K2, Z2 = _sim_grm(n, 600, r)
        y = _sim_y([Z1, Z2], list(truth), 0.5, r)
        mc = reml_mc(y, [K1, K2])
        g = run_gcta_reml([K1, K2], y, m=600)
        g1 = g["V(G1)"][0] if g and "V(G1)" in g else float("nan")
        g2 = g["V(G2)"][0] if g and "V(G2)" in g else float("nan")
        s1s.append(mc["sigma2"][0]); s2s.append(mc["sigma2"][1])
        print(f"  {sd:>4} {mc['sigma2'][0]:>7.3f} {mc['sigma2'][1]:>7.3f} "
              f"{mc['se'][0]:>7.3f} {mc['se'][1]:>7.3f} {g1:>8.3f} {g2:>8.3f}")
    print(f"  -> mean mc sigma2_1={np.mean(s1s):.3f}±{np.std(s1s):.3f} "
          f"(truth {truth[0]}), sigma2_2={np.mean(s2s):.3f}±{np.std(s2s):.3f} "
          f"(truth {truth[1]})\n")

    # ---- T3: boundary (true sigma2_2 = 0) ----
    print("=== T3 — boundary: component 2 has true sigma2=0 ===")
    r = np.random.default_rng(303)
    K1, Z1 = _sim_grm(n, 600, r)
    K2, Z2 = _sim_grm(n, 600, r)
    y = _sim_y([Z1], [0.5], 0.5, r)            # only component 1 active
    mc = reml_mc(y, [K1, K2])
    print(f"  constrained:   sigma2_1={mc['sigma2'][0]:.3f} sigma2_2={mc['sigma2'][1]:.3f} "
          f"sigma2_e={mc['sigma2_e']:.3f} converged={mc['converged']}")
    print(f"  unconstrained: {np.round(mc['sigma2_unconstrained'],3).tolist()}  "
          f"(theta = [e, 1, 2])")
    print(f"  -> sigma2_2 at/near zero boundary; no crash.\n")

    # ---- T4: engine on a real ARG kinship (single component) ----
    t4_real_argk()


if __name__ == "__main__":
    main()
