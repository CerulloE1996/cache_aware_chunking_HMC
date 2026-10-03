#!/bin/bash
##
## ---- Paper 1 (3 Oct 2026): BayesMVP's AVX-512 maths functions vs AVX2 vs stan::math (same installed build, same
##      chunked case: N = 50,000), 1 chain (100 chunks, CPU 20) and 80 chains (500 chunks, CPUs 16-95, one thread
##      per physical core, no SMT), two repeats each; two-run difference timing as in E4.
##
D="$(cd "$(dirname "$0")" && pwd)"; cd "$(dirname "$D")"
LOG="$D/log.txt"
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" >> "$LOG"; }
for repeat in 1 2; do
    for vect_type in AVX512 AVX2 Stan; do
        for spec in "1 100 30 20" "80 500 50 16-95"; do
            set -- $spec
            label="${vect_type}_chains$1_chunks$2_rep$repeat"
            log "START $label"
            MECH_VECT_TYPE=$vect_type MECH_ALGORITHM=MD_BayesMVP MECH_N=50000 MECH_CHUNKS=$2 MECH_CHAINS=$1 \
            MECH_THREADS_PER_CHAIN=1 MECH_N_ITER=$3 MECH_N_ITER_SHORT=2 MECH_SEED=1000 MECH_LABEL="$label" \
            MECH_RESULTS_FILE="$D/times.csv" OMP_NUM_THREADS=1 \
                nice -n 19 taskset -c $4 timeout 1800 Rscript "$D/mechanism_case_vect_type.R" >> "$D/R_output.txt" 2>&1 < /dev/null
            log "END   $label exit=$?"
        done
    done
done
log "ALL DONE"
