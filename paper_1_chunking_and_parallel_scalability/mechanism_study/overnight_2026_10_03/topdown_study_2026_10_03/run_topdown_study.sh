#!/bin/bash
##
## ---- Paper 1: Zen 4 top-down (memory-bound vs core-bound) and floating-point op rate of each workload at full SMT load,
##      local-HPC, N = 50,000, 180 chains with one thread each (E4 settings). Two passes per workload (PMC_EVENT_SET=topdown, fp),
##      each with six counters (no multiplexing). R workloads are windowed to the sampling call by mechanism_case.R's snapshots.
##      Waits for the temperature study to finish first, so the two do not disturb each other.
##
S="$(cd "$(dirname "$0")" && pwd)"; M="$(dirname "$S")"; cd "$M"
until grep -q "ALL DONE" "$M/temperature_study_2026_10_03/log.txt" 2>/dev/null; do sleep 20; done
for set in topdown fp; do
    echo "$(date '+%H:%M:%S') START Mplus_standard $set" >> "$S/log.txt"
    ( cd "$S/mplus_run" && PMC_EVENT_SET=$set timeout 1800 "$M/pmc_stat_topdown" "$S/counts_$set.csv" Mplus_standard -- \
          /snap/bin/mpdemo model.inp >> "$S/R_output.txt" 2>&1 < /dev/null )
    for spec in "MD_BayesMVP_chunks1 MD_BayesMVP 1 12" "MD_BayesMVP_chunking_500 MD_BayesMVP 500 150" \
                "AD_Stan_chunks1 AD_Stan 1 18" "AD_Stan_tape_chunked_250 AD_Stan_tape_chunked 250 45"; do
        set -- $spec
        echo "$(date '+%H:%M:%S') START $1 $set" >> "$S/log.txt"
        PMC_EVENT_SET=$set MECH_ALGORITHM=$2 MECH_N=50000 MECH_CHUNKS=$3 MECH_CHAINS=180 MECH_THREADS_PER_CHAIN=1 MECH_N_ITER=$4 \
        MECH_N_ITER_SHORT=2 MECH_SEED=1000 MECH_LABEL=$1 MECH_RESULTS_FILE="$S/times_$set.csv" OMP_NUM_THREADS=1 \
            timeout 1800 ./pmc_stat_topdown "$S/counts_$set.csv" "$1" -- Rscript mechanism_case.R >> "$S/R_output.txt" 2>&1 < /dev/null
        echo "$(date '+%H:%M:%S') END   $1 $set exit=$?" >> "$S/log.txt"
    done
done
echo "$(date '+%H:%M:%S') ALL DONE" >> "$S/log.txt"
