#!/bin/bash
##
## ---- Paper 1 (3 Oct 2026): continuous full-load CPU temperature runs on the local-HPC (follow-up to part 1 of
##      run_overnight_measurements.sh). There, each workload ran as a loop of short sampling calls (~50 s), with
##      ~5-8 s of single-threaded R work between calls, so the CPU cooled briefly every minute. Here, chunked BayesMVP
##      runs as ONE sampling call of ~600 s (N = 50,000, 180 chains with one thread each), after the CPU has cooled to
##      within 2 C (Tctl) of idle. (The Stan arms keep the full nuisance trace, ~0.43 GB per iteration for 180 chains at
##      N = 50,000, so one ~600 s Stan call would need ~140 GB of RAM; they are not repeated here.)
##      Starts after run_overnight_measurements.sh has finished, and only while the desktop has been idle
##      for at least 30 minutes.
##
O="$(cd "$(dirname "$0")" && pwd)"; M="$(dirname "$O")"; cd "$M"
C="$O/temperature_continuous"; mkdir -p "$C"
LOG="$O/log_continuous.txt"
log() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*" >> "$LOG"; }
load_now() { cut -d' ' -f1 /proc/loadavg; }
load_below() { awk -v current_load="$(load_now)" -v limit="$1" 'BEGIN { exit !(current_load < limit) }'; }
x_idle_seconds() { DISPLAY=:1 python3 "$O/x_idle_seconds.py" 2>/dev/null || echo 0; }
##
log "waiting for run_overnight_measurements.sh to finish"
until grep -q "ALL DONE" "$O/log.txt"; do sleep 60; done
until [ "$(x_idle_seconds)" -ge 1800 ] && load_below 12; do
    log "waiting: desktop idle $(x_idle_seconds) s, load $(load_now)"; sleep 120
done
log "starting"
##
echo "idle_baseline" > "$C/phase.txt"
trap 'rm -f "$C/phase.txt"' EXIT
taskset -c 191 "$M/temp_power_log" "$C/temperatures.csv" "$C/phase.txt" /sys/class/hwmon/hwmon5 &
sleep 300
baseline_Tctl=$(awk -F, 'NR > 1 && $2 == "idle_baseline" { s += $3; n++ } END { printf "%.2f", s / n }' \
                    "$C/temperatures.csv")
log "idle baseline Tctl = $baseline_Tctl"
##
fn_cool_down() {   ## label
    echo "idle_before_$1" > "$C/phase.txt"
    local waited=0 recent
    while true; do
        sleep 10; waited=$((waited + 10))
        recent=$(tail -30 "$C/temperatures.csv" | awk -F, '{ s += $3; n++ } END { printf "%.2f", s / n }')
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
fn_continuous_R() {   ## label, algorithm, chunks, iterations of the single sampling call
    fn_cool_down "$1"
    echo "$1" > "$C/phase.txt"; log "START $1"
    MECH_ALGORITHM=$2 MECH_N=50000 MECH_CHUNKS=$3 MECH_CHAINS=180 MECH_THREADS_PER_CHAIN=1 MECH_N_ITER=$4 \
    MECH_SECONDS=1 MECH_SEED=1000 MECH_LABEL=$1 MECH_RESULTS_FILE="$C/calls.csv" OMP_NUM_THREADS=1 \
        timeout 2400 ./umc_stat "$C/umc.csv" "$1" -- ./pmc_stat "$C/counts.csv" "$1" -- \
        Rscript mechanism_case_timed.R >> "$C/R_output.txt" 2>&1 < /dev/null
    log "END   $1 exit=$?"
}
##
## ~600 s in one call: 0.50 s per iteration (part 1 of the overnight runs):
fn_continuous_R MD_BayesMVP_chunking_500_continuous MD_BayesMVP 500 1200
echo "idle_after" > "$C/phase.txt"; sleep 120
rm -f "$C/phase.txt"; sleep 3
log "ALL DONE"
























