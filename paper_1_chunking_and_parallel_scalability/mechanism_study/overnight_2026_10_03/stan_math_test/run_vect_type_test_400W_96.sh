#!/bin/bash
##
## ---- Paper 1 (3 Oct 2026): the maths-function test (AVX-512 vs AVX2 vs stan::math, chunked BayesMVP, N = 50,000)
##      repeated after the CPU package power limit was raised from 360 W to 400 W (HSMP, by Enzo): 1 chain
##      (100 chunks, CPU 20), 96 chains (500 chunks, CPUs 0-95: one chain per physical core) and 180 chains (500 chunks, all CPUs), two repeats
##      each. CPU temperature and package power are logged every second (power_400W.csv; phase = run label).
##
D="$(cd "$(dirname "$0")" && pwd)"; M="$(dirname "$(dirname "$D")")"; cd "$(dirname "$D")"
LOG="$D/log_400W.txt"
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" >> "$LOG"; }
PHASE="$D/phase_400W.txt"; echo "start" > "$PHASE"
taskset -c 191 "$M/temp_power_log" "$D/power_400W.csv" "$PHASE" /sys/class/hwmon/hwmon5 &
trap 'rm -f "$PHASE"' EXIT
for repeat in 1 2; do
    for vect_type in AVX512 AVX2 Stan; do
        for spec in "1 100 30 20" "96 500 50 0-95" "180 500 50 0-191"; do
            set -- $spec
            label="${vect_type}_chains$1_chunks$2_rep$repeat"
            echo "$label" > "$PHASE"; log "START $label"
            MECH_VECT_TYPE=$vect_type MECH_ALGORITHM=MD_BayesMVP MECH_N=50000 MECH_CHUNKS=$2 MECH_CHAINS=$1 \
            MECH_THREADS_PER_CHAIN=1 MECH_N_ITER=$3 MECH_N_ITER_SHORT=2 MECH_SEED=1000 MECH_LABEL="$label" \
            MECH_RESULTS_FILE="$D/times_400W.csv" OMP_NUM_THREADS=1 \
                nice -n 19 taskset -c $4 timeout 1800 Rscript "$D/mechanism_case_vect_type.R" \
                >> "$D/R_output_400W.txt" 2>&1 < /dev/null
            log "END   $label exit=$?"
        done
    done
done
echo "idle_after" > "$PHASE"; sleep 5
log "ALL DONE"
