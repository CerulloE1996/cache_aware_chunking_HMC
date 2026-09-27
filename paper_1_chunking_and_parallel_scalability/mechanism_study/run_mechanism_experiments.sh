#!/bin/bash
##
## ---- Runs every case in cases_<device>.csv under pmc_stat, one fresh R process per case, skipping cases already done.
##      Hardware counts: mechanism_counts_<device>.csv; sampling times: mechanism_times_<device>.csv; log: mechanism_log_<device>.txt
##      Usage: bash run_mechanism_experiments.sh <HPC|Laptop> [cases file] [umc]
##      With "umc" (HPC only), each case also runs under umc_stat, which records DRAM traffic at the memory controllers
##      (mechanism_umc_<device>.csv).
##
cd "$(dirname "$0")"
device="$1"
cases_file="${2:-cases_${device}.csv}"
use_umc="${3:-}"
umc_file="$PWD/mechanism_umc_${device}.csv"
counts_file="$PWD/mechanism_counts_${device}.csv"
times_file="$PWD/mechanism_times_${device}.csv"
log_file="$PWD/mechanism_log_${device}.txt"
[ -f "$counts_file" ] || echo "label,wall_seconds,exit_status,max_rss_kb,cycles,instructions,fills_local_L2,fills_local_L3,fills_other_CCX,fills_DRAM" > "$counts_file"
tail -n +2 "$cases_file" | tr -d '\r"' | while IFS=, read -r label experiment algorithm N chunks chains tpc n_iter cpus; do
    if [ -f "$times_file" ] && grep -q "^\"$label\"," "$times_file"; then continue; fi
    echo "$(date '+%H:%M:%S') START $label (load $(cut -d' ' -f1 /proc/loadavg))" >> "$log_file"
    pin=""; [ "$cpus" != "all" ] && pin="taskset -c $cpus"
    umc_prefix=""; [ "$use_umc" = "umc" ] && umc_prefix="./umc_stat $umc_file $label --"
    MECH_ALGORITHM=$algorithm MECH_N=$N MECH_CHUNKS=$chunks MECH_CHAINS=$chains MECH_THREADS_PER_CHAIN=$tpc \
    MECH_N_ITER=$n_iter MECH_N_ITER_SHORT=2 MECH_SEED=1000 MECH_LABEL=$label MECH_RESULTS_FILE="$times_file" OMP_NUM_THREADS=1 \
        timeout 3600 $pin $umc_prefix ./pmc_stat "$counts_file" "$label" -- Rscript mechanism_case.R >> "$log_file.R_output" 2>&1 < /dev/null
    status=$?
    echo "$(date '+%H:%M:%S') END   $label exit=$status" >> "$log_file"
done
echo "$(date '+%H:%M:%S') ALL DONE" >> "$log_file"
