#!/bin/bash
## Stan working set per individual: build and run the probe harness (run from working_set/stan).
## Requires BridgeStan 2.6.2 (with its Stan Math and TBB) and R with cmdstanr; BS is the BridgeStan directory.
set -e
BS=${BS:-$HOME/.bridgestan/bridgestan-2.6.2}
MODEL=../../stan_models/LC_MVP_bin_PartialLog_v5_reduce_sum_static.stan
##
## ---- 1. Model C++ from the paper's .stan file (stanc 2.36.0; --allow-undefined, no optimisation flag, as in the benchmarks),
##         then the three probes (model_probes.diff) applied to a copy
cp $MODEL model_rs.stan
$BS/bin/stanc --allow-undefined --o=model_rs.hpp model_rs.stan
cp model_rs.hpp model_rs_probed.hpp
patch model_rs_probed.hpp model_probes.diff
##
## ---- 2. Data and parameter JSON files (uses fn_paper1_stan_data() and the seed-123 datasets of the benchmarks)
mkdir -p json
Rscript make_json.R
##
## ---- 3. Build: malloc interposition (C) and the harness (C++), with the BridgeStan make/local flags of the local-HPC
gcc -O2 -c heap_hook.c -o heap_hook.o
g++ -O3 -march=native -mtune=native -mfma -mavx -mavx2 -mavx512f -mavx512vl -mavx512dq \
    -fno-math-errno -fno-signed-zeros -fno-trapping-math -DNDEBUG -DBOOST_DISABLE_ASSERTS \
    -Wno-deprecated-declarations -Wno-sign-compare -Wno-ignored-attributes -std=c++17 -pthread -D_REENTRANT \
    -Wno-class-memaccess -DSTAN_THREADS \
    -I $BS/stan/lib/stan_math/lib/tbb_2020.3/include -I $BS/stan/src -I $BS/stan/lib/rapidjson_1.1.0/ \
    -I $BS/stan/lib/stan_math/ -I $BS/stan/lib/stan_math/lib/eigen_3.4.0 -I $BS/stan/lib/stan_math/lib/boost_1.84.0 \
    -I $BS/stan/lib/stan_math/lib/sundials_6.1.1/include -I $BS/stan/lib/stan_math/lib/sundials_6.1.1/src/sundials -I . \
    ws_harness.cpp heap_hook.o -o ws_harness \
    -Wl,-L,$BS/stan/lib/stan_math/lib/tbb -Wl,-rpath,$BS/stan/lib/stan_math/lib/tbb -ltbb
##
## ---- 4. Run: one chunk of n rows (chunk_size = N) for n = 500 to 20,000, plus N = 20,000 in 4 and 32 partial sums;
##         parameter points A (u ~ U(0, 1), true Se/Sp/prevalence and correlations) and B (all unconstrained parameters 0);
##         three repetitions each (the third is used)
rm -f results.txt time.txt
for pt in A B; do
    for c in "500 500" "1000 1000" "2500 2500" "5000 5000" "10000 10000" "20000 20000" "20000 5000" "20000 1000"; do
        set -- $c
        /usr/bin/time -f "TIME label=${pt}_N$1_cs$2 wall=%e maxrss_kb=%M" \
            ./ws_harness json/data_N$1_cs$2.json json/init_${pt}_N$1.json ${pt}_N$1_cs$2 3 >> results.txt 2>> time.txt
    done
done
##
## ---- 5. Per-individual slopes (fit_summary.txt)
python3 fit_slopes.py results.txt > fit_summary.txt
cat fit_summary.txt
