#!/bin/bash
# overnight_sweep.sh — multi-seed bed-coverage sweep + h² sweep at f=0.5.
#
# Runs:
#   1. Multi-seed bed-coverage sweep on the balanced h² case across 5 seeds,
#      5 BED-coverage levels, both sequential and joint-EM SARGENT.
#   2. h² sweep at f=0.5 across 5 h² splits, 3 seeds, both methods.
#
# Output: /Users/acs/test-sargent/overnight/{multiseed,h2sweep}/<case>/{seq,em}.predict.json
#         /Users/acs/test-sargent/overnight/summary.md (final aggregated table)

set -e

SARGENT=/Users/acs/git-src/sargent/build/src/progs/sargent
SIM=/Users/acs/git-src/sargent/scripts/simulate.py
PRED=/Users/acs/git-src/sargent/test/eval/predict.py
AGG=/Users/acs/git-src/sargent/test/eval/overnight_aggregate.py
BASE=/Users/acs/test-sargent/overnight
rm -rf "$BASE" && mkdir -p "$BASE"

start_time=$(date +%s)
log() {
    elapsed=$(($(date +%s) - start_time))
    h=$((elapsed / 3600))
    m=$(((elapsed % 3600) / 60))
    s=$((elapsed % 60))
    printf "[%dh%02dm%02ds] %s\n" "$h" "$m" "$s" "$*"
}

# Run sequential + predict
run_seq() {
    local dir=$1
    $SARGENT --silent --nprobe 10 --sequential --outroot "$dir/seq" \
        "$dir/case.phen" "$dir/case.trees" "$dir/case" \
        > "$dir/seq_sargent.log" 2>&1
    mamba run -n sargent python "$PRED" \
        --case-prefix "$dir/case" --sargent-prefix "$dir/seq" \
        > "$dir/seq_predict.log" 2>&1
}

# Run joint EM + predict
run_em() {
    local dir=$1
    $SARGENT --silent --nprobe 10 --em 50 --outroot "$dir/em" \
        "$dir/case.phen" "$dir/case.trees" "$dir/case" \
        > "$dir/em_sargent.log" 2>&1
    mamba run -n sargent python "$PRED" \
        --case-prefix "$dir/case" --sargent-prefix "$dir/em" \
        > "$dir/em_predict.log" 2>&1
}

# ===== Multi-seed bed-coverage sweep (balanced h2 case) =====
log "=== Multi-seed bed-coverage sweep (balanced h² split) ==="
for seed in 42 43 44 45 46; do
    for f in 1.0 0.6 0.3 0.1 0.03; do
        dir="$BASE/multiseed/s${seed}_f${f}"
        mkdir -p "$dir"
        log "  multiseed seed=$seed f=$f: simulating"
        mamba run -n sargent python "$SIM" \
            --n-indiv 200 --n-test 100 --length 100000 \
            --recomb-rate 1e-8 --mut-rate 1e-8 --Ne 10000 \
            --h2-snp 0.40 --h2-arg 0.40 --n-causal 10 \
            --bed-coverage "$f" --seed "$seed" \
            --out "$dir/case" > "$dir/sim.log" 2>&1
        log "  multiseed seed=$seed f=$f: sequential SARGENT"
        run_seq "$dir"
        log "  multiseed seed=$seed f=$f: joint-EM SARGENT"
        run_em "$dir"
    done
done

# ===== h² sweep at f=0.5 =====
log "=== h² sweep at f=0.5 ==="
for seed in 42 43 44; do
    for case_spec in \
        "snp_only:0.80:0.00" \
        "snp_heavy:0.60:0.20" \
        "balanced:0.40:0.40" \
        "arg_heavy:0.20:0.60" \
        "arg_only:0.00:0.80"; do
        IFS=':' read -r name h2f h2r <<< "$case_spec"
        dir="$BASE/h2sweep/s${seed}_${name}"
        mkdir -p "$dir"
        log "  h2sweep seed=$seed $name (h2f=$h2f h2r=$h2r): simulating"
        mamba run -n sargent python "$SIM" \
            --n-indiv 200 --n-test 100 --length 100000 \
            --recomb-rate 1e-8 --mut-rate 1e-8 --Ne 10000 \
            --h2-snp "$h2f" --h2-arg "$h2r" --n-causal 10 \
            --bed-coverage 0.5 --seed "$seed" \
            --out "$dir/case" > "$dir/sim.log" 2>&1
        log "  h2sweep seed=$seed $name: sequential SARGENT"
        run_seq "$dir"
        log "  h2sweep seed=$seed $name: joint-EM SARGENT"
        run_em "$dir"
    done
done

log "=== Aggregating ==="
mamba run -n sargent python "$AGG" --base "$BASE" > "$BASE/summary.md" 2> "$BASE/summary.err"

log "=== Done ==="
log "Summary: $BASE/summary.md"
