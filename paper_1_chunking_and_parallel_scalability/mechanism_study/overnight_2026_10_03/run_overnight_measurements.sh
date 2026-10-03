#!/bin/bash
##
## ---- Paper 1: overnight measurements on the local-HPC (set up on 3 Oct 2026, 03:40).
##      Waits until the machine is free (1-min load average below 12 for 20 consecutive minutes),
##      starting no later than LATEST_START (HHMM, default 05:30) and not before EARLIEST_START (default
##      00:00); on 3 Oct only after the other runs started that night have been seen (or with SEEN_BUSY=1).
##      Then it runs:
##        1. long CPU temperature runs: N = 50,000, 180 chains with one thread each; each workload runs for
##           ~600 s of sampling, after the CPU has cooled to within 2 C (Tctl) of the idle baseline;
##        2. the rest of the top-down floating-point pass (stopped by hand at 03:25 on 3 Oct), and the
##           floating-point width pass (scalar / 128 / 256 / 512-bit uops) of all five workloads;
##        3. AVX-512 vs AVX2-only builds of BayesMVP (both from the same source snapshot), with the same
##           chunked case timed with each build (E4 timing: short and long runs).
##
O="$(cd "$(dirname "$0")" && pwd)"; M="$(dirname "$O")"; cd "$M"
LOG="$O/log.txt"
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" >> "$LOG"; }
load_now() { cut -d' ' -f1 /proc/loadavg; }
load_below() { awk -v current_load="$(load_now)" -v limit="$1" 'BEGIN { exit !(current_load < limit) }'; }
load_above() { awk -v current_load="$(load_now)" -v limit="$1" 'BEGIN { exit !(current_load > limit) }'; }
log "waiting for a free machine"
##
## ---- 0. Wait for a free machine at night ------------------------------------------------------------------
##
seen_busy=${SEEN_BUSY:-0}; free_minutes=0   ## SEEN_BUSY=1: the other runs have already been seen
earliest_start=$((10#${EARLIEST_START:-0})); latest_start=$((10#${LATEST_START:-530}))
while true; do
    load_above 20 && seen_busy=1
    [ "$(date +%Y-%m-%d)" != "2026-10-03" ] && seen_busy=1
    time_now=$((10#$(date +%H%M)))
    if [ $seen_busy = 1 ] && [ $time_now -le $latest_start ] && [ $time_now -ge $earliest_start ] \
       && load_below 12; then
        free_minutes=$((free_minutes + 1))
    else
        free_minutes=0
    fi
    [ $free_minutes -ge 20 ] && break
    sleep 60
done
log "machine free for 20 min; starting"
##
## ---- 1. Long CPU temperature runs --------------------------------------------------------------------------
##
T="$O/temperature_long"
echo "idle_baseline" > "$T/phase.txt"
trap 'rm -f "$T/phase.txt"' EXIT
taskset -c 191 "$M/temp_power_log" "$T/temperatures.csv" "$T/phase.txt" /sys/class/hwmon/hwmon5 &
sleep 300
baseline_Tctl=$(awk -F, 'NR > 1 && $2 == "idle_baseline" { s += $3; n++ } END { printf "%.2f", s / n }' \
                    "$T/temperatures.csv")
log "idle baseline Tctl = $baseline_Tctl"
##
## Cool-down: wait until the mean Tctl of the last 30 s is within 2 C of the baseline (120-900 s):
fn_cool_down() {   ## label
    echo "idle_before_$1" > "$T/phase.txt"
    local waited=0 recent
    while true; do
        sleep 10; waited=$((waited + 10))
        recent=$(tail -30 "$T/temperatures.csv" | awk -F, '{ s += $3; n++ } END { printf "%.2f", s / n }')
        if [ $waited -ge 120 ] && awk -v r="$recent" -v b="$baseline_Tctl" 'BEGIN { exit !(r <= b + 2) }'; then
            break
        fi
        if [ $waited -ge 900 ]; then
            log "cool-down before $1 ended at 900 s (Tctl $recent, baseline $baseline_Tctl)"; break
        fi
    done
    until load_below 12; do log "waiting before $1: load $(load_now)"; sleep 60; done
}
##
fn_temperature_R() {   ## label, algorithm, chunks, iterations per sampling call
    fn_cool_down "$1"
    echo "$1" > "$T/phase.txt"; log "START $1"
    MECH_ALGORITHM=$2 MECH_N=50000 MECH_CHUNKS=$3 MECH_CHAINS=180 MECH_THREADS_PER_CHAIN=1 MECH_N_ITER=$4 \
    MECH_SECONDS=600 MECH_SEED=1000 MECH_LABEL=$1 MECH_RESULTS_FILE="$T/calls.csv" OMP_NUM_THREADS=1 \
        timeout 2400 ./umc_stat "$T/umc.csv" "$1" -- ./pmc_stat "$T/counts.csv" "$1" -- \
        Rscript mechanism_case_timed.R >> "$T/R_output.txt" 2>&1 < /dev/null
    log "END   $1 exit=$?"
}
##
fn_cool_down Mplus_standard
echo "Mplus_standard" > "$T/phase.txt"; log "START Mplus_standard"
( cd "$T/mplus_run" && timeout 2400 "$M/umc_stat" "$T/umc.csv" Mplus_standard -- \
      "$M/pmc_stat" "$T/counts.csv" Mplus_standard -- /snap/bin/mpdemo model.inp \
      >> "$T/R_output.txt" 2>&1 < /dev/null )
log "END   Mplus_standard exit=$?"
fn_temperature_R MD_BayesMVP_chunking_500   MD_BayesMVP          500 100
fn_temperature_R AD_Stan_tape_chunked_250   AD_Stan_tape_chunked 250 30
fn_temperature_R MD_BayesMVP_chunks1        MD_BayesMVP          1   8
fn_temperature_R AD_Stan_chunks1            AD_Stan              1   12
echo "idle_after" > "$T/phase.txt"; sleep 300
rm -f "$T/phase.txt"; sleep 3
log "temperature runs done"
##
## ---- 2. Floating-point passes (pmc_stat_topdown) ----------------------------------------------------------
##
S="$M/topdown_study_2026_10_03"
fn_counter_pass() {   ## event set, label, algorithm, chunks, iterations
    until load_below 12; do log "waiting before $2 $1: load $(load_now)"; sleep 60; done
    log "START $2 $1"
    PMC_EVENT_SET=$1 MECH_ALGORITHM=$3 MECH_N=50000 MECH_CHUNKS=$4 MECH_CHAINS=180 MECH_THREADS_PER_CHAIN=1 \
    MECH_N_ITER=$5 MECH_N_ITER_SHORT=2 MECH_SEED=1000 MECH_LABEL=$2 MECH_RESULTS_FILE="$S/times_$1.csv" \
    OMP_NUM_THREADS=1 \
        timeout 1800 ./pmc_stat_topdown "$S/counts_$1.csv" "$2" -- Rscript mechanism_case.R \
        >> "$S/R_output.txt" 2>&1 < /dev/null
    log "END   $2 $1 exit=$?"
}
fn_counter_pass fp      MD_BayesMVP_chunking_500 MD_BayesMVP          500 150
fn_counter_pass fp      AD_Stan_chunks1          AD_Stan              1   18
fn_counter_pass fp      AD_Stan_tape_chunked_250 AD_Stan_tape_chunked 250 45
log "START Mplus_standard fpwidth"
( cd "$S/mplus_run" && PMC_EVENT_SET=fpwidth timeout 1800 "$M/pmc_stat_topdown" "$S/counts_fpwidth.csv" \
      Mplus_standard -- /snap/bin/mpdemo model.inp >> "$S/R_output.txt" 2>&1 < /dev/null )
log "END   Mplus_standard fpwidth exit=$?"
fn_counter_pass fpwidth MD_BayesMVP_chunks1      MD_BayesMVP          1   12
fn_counter_pass fpwidth MD_BayesMVP_chunking_500 MD_BayesMVP          500 150
fn_counter_pass fpwidth AD_Stan_chunks1          AD_Stan              1   18
fn_counter_pass fpwidth AD_Stan_tape_chunked_250 AD_Stan_tape_chunked 250 45
##
## ---- 3. AVX-512 vs AVX2-only builds of BayesMVP -------------------------------------------------------------
##
A="$O/avx_test"
avx512_pattern='^  override AVX_FLAGS += -mavx -mavx2 -mavx512f -mavx512cd -mavx512bw -mavx512dq -mavx512vl$'
for variant in avx512 avx2; do
    rm -rf "$A/src_$variant" "$A/lib_$variant"
    cp -r "$A/BayesMVP_source_snapshot_0340" "$A/src_$variant"; mkdir -p "$A/lib_$variant"
    rm -f "$A/src_$variant"/src/*.o "$A/src_$variant"/src/*.so
    if [ $variant = avx2 ]; then
        sed -i 's/^override CPU_BASE_FLAGS = -O3  -march=native  -mtune=native$/& -mno-avx512f/' \
            "$A/src_$variant/src/Makevars"
        sed -i "s/$avx512_pattern/  override AVX_FLAGS += -mavx -mavx2/" "$A/src_$variant/src/Makevars"
    fi
    log "flags $variant: $(grep -h '^override CPU_BASE_FLAGS\|^  override AVX_FLAGS += ' \
                               "$A/src_$variant/src/Makevars" | tr '\n' '|')"
    ( cd "$A" && MAKEFLAGS=-j32 timeout 3600 R CMD INSTALL --preclean --library="$A/lib_$variant" "src_$variant" \
          > "$A/install_$variant.txt" 2>&1 )
    log "build $variant exit=$? zmm_instructions=$(objdump -d "$A/lib_$variant/BayesMVP/libs/BayesMVP.so" \
                                                    2>/dev/null | grep -c zmm)"
    log "R_LIBS=$A/lib_$variant loads BayesMVP from: \
$(R_LIBS="$A/lib_$variant" Rscript -e 'cat(find.package("BayesMVP"))' 2>&1)"
done
for repeat in 1 2; do
    for variant in avx512 avx2; do
        for spec in "180 500 50" "1 100 30"; do
            set -- $spec
            label="${variant}_chains$1_chunks$2_rep$repeat"
            until load_below 12; do log "waiting before $label: load $(load_now)"; sleep 60; done
            log "START $label"
            R_LIBS="$A/lib_$variant" MECH_ALGORITHM=MD_BayesMVP MECH_N=50000 MECH_CHUNKS=$2 MECH_CHAINS=$1 \
            MECH_THREADS_PER_CHAIN=1 MECH_N_ITER=$3 MECH_N_ITER_SHORT=2 MECH_SEED=1000 MECH_LABEL="$label" \
            MECH_RESULTS_FILE="$A/times.csv" OMP_NUM_THREADS=1 \
                timeout 1800 Rscript mechanism_case.R >> "$A/R_output.txt" 2>&1 < /dev/null
            log "END   $label exit=$?"
        done
    done
done
log "ALL DONE"























