#!/bin/bash
# medium_tests.sh — run a small medium-scale follow-up once the overnight
# sweep at /Users/acs/test-sargent/overnight/summary.md completes.
#
# 3 BED-coverage levels (f = 1.0, 0.3, 0.1) on the balanced h² case,
# single seed (42), both sequential and joint-EM SARGENT.
#
# Medium scale: n_train=300, n_test=150, L=300kb (~ 500 trees, 15 causal SNPs)
# — about 4-5x more compute than the overnight debug scale.

set -e

SARGENT=/Users/acs/git-src/sargent/build/src/progs/sargent
SIM=/Users/acs/git-src/sargent/scripts/simulate.py
PRED=/Users/acs/git-src/sargent/test/eval/predict.py
SUMMARY=/Users/acs/test-sargent/overnight/summary.md
BASE=/Users/acs/test-sargent/overnight/medium

# Wait for the overnight script's summary.md to exist and contain its Status
# section (which is the last thing the aggregator writes).
until [ -f "$SUMMARY" ] && grep -q "## Status" "$SUMMARY"; do
    sleep 60
done
# Pad slightly to ensure the overnight script has fully exited
sleep 10

start_time=$(date +%s)
log() {
    elapsed=$(($(date +%s) - start_time))
    h=$((elapsed / 3600))
    m=$(((elapsed % 3600) / 60))
    s=$((elapsed % 60))
    printf "[%dh%02dm%02ds] %s\n" "$h" "$m" "$s" "$*"
}

mkdir -p "$BASE"
log "=== Medium-scale follow-up (n_train=300, L=300kb, single seed=42) ==="

for f in 1.0 0.3 0.1; do
    dir="$BASE/f${f}"
    mkdir -p "$dir"
    log "  medium f=$f: simulating"
    mamba run -n sargent python "$SIM" \
        --n-indiv 300 --n-test 150 --length 300000 \
        --recomb-rate 1e-8 --mut-rate 1e-8 --Ne 10000 \
        --h2-snp 0.40 --h2-arg 0.40 --n-causal 15 \
        --bed-coverage "$f" --seed 42 \
        --out "$dir/case" > "$dir/sim.log" 2>&1
    log "  medium f=$f: sequential SARGENT"
    $SARGENT --silent --nprobe 10 --sequential --outroot "$dir/seq" \
        "$dir/case.phen" "$dir/case.trees" "$dir/case" \
        > "$dir/seq_sargent.log" 2>&1
    mamba run -n sargent python "$PRED" \
        --case-prefix "$dir/case" --sargent-prefix "$dir/seq" \
        > "$dir/seq_predict.log" 2>&1
    log "  medium f=$f: joint-EM SARGENT"
    $SARGENT --silent --nprobe 10 --em 50 --outroot "$dir/em" \
        "$dir/case.phen" "$dir/case.trees" "$dir/case" \
        > "$dir/em_sargent.log" 2>&1
    mamba run -n sargent python "$PRED" \
        --case-prefix "$dir/case" --sargent-prefix "$dir/em" \
        > "$dir/em_predict.log" 2>&1
done

log "=== Aggregating medium results and appending to summary.md ==="
mamba run -n sargent python - <<'PY' >> "$SUMMARY"
import json, os
print()
print("## Medium-scale follow-up (single seed=42)")
print()
print("n_train=300, n_test=150, length=300kb, regime=sparse, balanced "
      "(h²_fixed=h²_random=0.40), 15 causal SNPs.")
print()
print("| f | BayesR | GBLUP | ARG-BLUP-uni | Seq SARGENT | **Joint EM SARGENT** |")
print("|---:|:---:|:---:|:---:|:---:|:---:|")
for f in ("1.0", "0.3", "0.1"):
    base = f"/Users/acs/test-sargent/overnight/medium/f{f}"
    seq = json.load(open(f"{base}/seq.predict.json")) if os.path.exists(f"{base}/seq.predict.json") else None
    em  = json.load(open(f"{base}/em.predict.json"))  if os.path.exists(f"{base}/em.predict.json")  else None
    def r2(d, k):
        if d is None or k not in d: return "-"
        return f"{d[k]['r2']:+.3f}"
    if seq:
        print(f"| {f} | {r2(seq,'bayesR')} | {r2(seq,'gblup')} | {r2(seq,'arg_blup_uniform')} | "
              f"{r2(seq,'sargent')} | **{r2(em,'sargent') if em else '-'}** |")
    else:
        print(f"| {f} | - | - | - | - | - |")
PY

log "=== Medium tests done ==="
