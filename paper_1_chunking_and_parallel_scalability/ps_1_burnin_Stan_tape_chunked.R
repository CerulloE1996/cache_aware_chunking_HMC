
## ps_1_burnin_Stan_tape_chunked.R
##
## PS1 (burn-in) for the Stan TAPE-CHUNKED model: the same burn-in timing as the BayesMVP burn-in study
## (ps_1_burnin_optimizing_N_chunks_and_WCP_threads.R), applied to LC_MVP_bin_PartialLog_v5_reduce_sum_static.stan, the model
## of the AD_Stan_tape_chunked and AD_Stan_WCP arms of the Paper 1 sampling study (alg_paper_1_chunking_WCP_par_scaling.R).
##
## Settings are READ from the two existing runners when this file runs (neither runner is run or edited):
##   - burn-in chain counts, WCP threads per chain (keyed by N and chain count), L, eps, timed iterations per N, n_runs,
##     thread limit, shared tau_ii: from ps_1_burnin_optimizing_N_chunks_and_WCP_threads.R (burnin_benchmark_settings);
##   - chunk counts per N: stan_chunk_candidates of alg_paper_1_chunking_WCP_par_scaling.R.
## num_chunks = the number of tape chunks (reduce_sum_static partial sums, chunk_size = ceiling(N / num_chunks)).
## WCP = reduce_sum threads within a chain, from a TBB pool of n_chains x n_WCP threads (see the header of
## R_fn_run_ps_1_burnin_benchmark_Stan_tape_chunked.R). WCP = 1 is the AD_Stan_tape_chunked arm; WCP > 1 is the AD_Stan_WCP arm.
##
## Output (one file per N, saved after every configuration/run; a restarted run skips what is saved):
##   burnin_outputs/<HPC_ or Laptop_>ps1_burnin_benchmark_Stan_tape_chunked_N<N>_L20_n_runs3.rds
## with the BayesMVP burn-in columns plus implementation = "Stan_tape_chunked", configuration_set, stan_chunk_size,
## TBB_pool_threads, time_worker_setup_seconds and the per-iteration times (sec_per_iter_each).
##
## Run it on an otherwise IDLE machine - anything else running makes the timings meaningless.
##
## Headless controls (environment variables; unset = the defaults below):
##   PS1_BURNIN_N                  N's to run, comma-separated (default: all four)
##   PS1_BURNIN_CONFIGURATION_SET  "chunk_grid" (default), "WCP_only" or "all"
##   PS1_BURNIN_OUTPUT_DIR         output folder (default: burnin_outputs)
##   PS1_BURNIN_DRY_RUN            "TRUE" prints the grid and the configuration counts only (no model, no MCMC)
##   PS1_BURNIN_SMOKE_TEST         "TRUE" runs ONE tiny configuration: N = 500, 4 chains, WCP 1, 10 chunks, 1 run, 2 timed iterations
##
# rm(list = ls())
# .rs.restartR()
##
##
## ------- Set options / paths:  ------------------------------------------------------------------------------------------------------
##
{
      options(scipen = 999)
      n_total_threads <- parallel::detectCores()
      computer <- ifelse(n_total_threads > 16, "Local_HPC", "Laptop")
      ##
      algorithm_study_dir <- path.expand("~/Documents/Work/PhD_work/Alg_paper_analysis")
      ps1_dir             <- file.path(algorithm_study_dir, "paper_1_chunking_and_parallel_scalability")
      ##
      burnin_runner_file_path   <- file.path(ps1_dir, "ps_1_burnin_optimizing_N_chunks_and_WCP_threads.R")
      sampling_runner_file_path <- file.path(ps1_dir, "alg_paper_1_chunking_WCP_par_scaling.R")
      ##
      burnin_output_dir <- Sys.getenv("PS1_BURNIN_OUTPUT_DIR", unset = file.path(ps1_dir, "burnin_outputs"))
      json_dir          <- file.path(tempdir(), "ps1_burnin_Stan_tape_chunked_json")
      ##
      run_as_dry_run    <- identical(toupper(Sys.getenv("PS1_BURNIN_DRY_RUN",    unset = "FALSE")), "TRUE")
      run_as_smoke_test <- identical(toupper(Sys.getenv("PS1_BURNIN_SMOKE_TEST", unset = "FALSE")), "TRUE")
      configuration_set_to_run <- Sys.getenv("PS1_BURNIN_CONFIGURATION_SET", unset = "chunk_grid")
      if (!configuration_set_to_run %in% c("chunk_grid", "WCP_only", "all")) {
            stop("PS1_BURNIN_CONFIGURATION_SET must be chunk_grid, WCP_only or all.")
      }
}
##
## ------- Functions (defines only):  --------------------------------------------------------------------------------------------------
##
{
      `%>%` <- magrittr::`%>%`   ## used by R_fn_ps1_burnin_make_configuration_grid()
      ##
      source(file.path(ps1_dir, "R_fn_run_ps_1_burnin_benchmark.R"))
      source(file.path(ps1_dir, "R_fn_run_ps_1_burnin_benchmark_Stan_tape_chunked.R"))
      ##
      ## ---- the sampling study's helpers (compiled-model cache, Stan data list), in their own environment:
      ##
      paper1_sampling_helpers <- new.env(parent = globalenv())
      sys.source(file = file.path(ps1_dir, "R_fns_alg_paper_1_chunking_WCP_par_scaling.R"), envir = paper1_sampling_helpers)
}
##
## ------- Benchmark settings, read from the two runners. EDIT the runners, not these lines:  -----------------------------------------
##
{
      burnin_benchmark_settings <- R_fn_ps1_burnin_read_BayesMVP_burnin_settings( burnin_runner_file_path = burnin_runner_file_path,
                                                                                  computer                = computer)
      ##
      BayesMVP_burnin_chunks_given_N <- burnin_benchmark_settings$num_chunks_burnin_vec_given_N
      stan_chunk_candidates <- R_fn_ps1_burnin_read_assignment_from_runner( runner_file_path = sampling_runner_file_path,
                                                                            object_name      = "stan_chunk_candidates")
      burnin_benchmark_settings$num_chunks_burnin_vec_given_N <- stan_chunk_candidates
      ##
      ## ---- Stan inits as in the sampling runner (stan_via_NicoStan$init_radius = 0.01):
      ##
      burnin_benchmark_settings$stan_init_radius <- 0.01
      ##
      ## ---- NicoStan build: the Stan arm needs the load-once burn-in worker (one Stan model + data per chain, constructed once per
      ##      worker). The installed NicoStan headers show which worker the installed library was built from; the library md5 is
      ##      saved in every row. PS1_BURNIN_ALLOW_PER_ITERATION_MODEL_RELOAD = "TRUE" allows the old build (tests only).
      ##
      NicoStan_burnin_header_path <- file.path(find.package("NicoStan"), "include", "NicoStan", "runtime", "MCMC",
                                               "EHMC_burn_multi_thread_samp_fns_RCPP.hpp")
      burnin_benchmark_settings$Stan_model_loading <- if (file.exists(NicoStan_burnin_header_path) &&
                                                          any(grepl(pattern = "Stan_models_per_chain", x = readLines(NicoStan_burnin_header_path),
                                                                    fixed = TRUE))) "once_per_worker" else "every_iteration"
      burnin_benchmark_settings$NicoStan_library_md5 <- unname(tools::md5sum(file.path(find.package("NicoStan"), "libs", "NicoStan.so")))
      if (burnin_benchmark_settings$Stan_model_loading != "once_per_worker" &&
          !identical(Sys.getenv("PS1_BURNIN_ALLOW_PER_ITERATION_MODEL_RELOAD"), "TRUE")) {
            stop("The installed NicoStan (", find.package("NicoStan"), ") reloads the Stan model in every burn-in iteration; ",
                 "install the load-once build before running the Stan tape-chunked burn-in.")
      }
      ##
      ## ---- Timed burn-in iterations per N: NULL = the BayesMVP burn-in study's n_timed_iters_given_N (same trajectory lengths,
      ##      because tau_ii is drawn from the same seeds). A list here, e.g. list("500" = 200, "2500" = 40, "10000" = 10, "50000" = 5),
      ##      replaces it for the Stan model only.
      ##
      # stan_n_timed_iters_given_N <- NULL
      ## the saved Stan runs at N = 10,000 and 50,000 used 5 and 3 timed iterations (3 runs each, 30 Sep 2026):
      stan_n_timed_iters_given_N <- list("10000" = 5, "50000" = 3)
      ##
      ## ---- headless override, e.g. PS1_BURNIN_STAN_N_TIMED_ITERS="500=200,2500=40,10000=5,50000=3":
      ##
      if (nzchar(Sys.getenv("PS1_BURNIN_STAN_N_TIMED_ITERS"))) {
            timed_iters_pairs <- strsplit(x = strsplit(x = Sys.getenv("PS1_BURNIN_STAN_N_TIMED_ITERS"), split = ",", fixed = TRUE)[[1]],
                                          split = "=", fixed = TRUE)
            stan_n_timed_iters_given_N <- stats::setNames(object = lapply(X = timed_iters_pairs, FUN = function(pair) as.numeric(pair[2])),
                                                          nm     = vapply(X = timed_iters_pairs, FUN = function(pair) pair[1], FUN.VALUE = character(1)))
      }
      if (!is.null(stan_n_timed_iters_given_N)) {
            burnin_benchmark_settings$n_timed_iters_given_N <- utils::modifyList(x = burnin_benchmark_settings$n_timed_iters_given_N,
                                                                                val = stan_n_timed_iters_given_N)
      }
      ##
      N_vec_to_benchmark <- c(500, 2500, 10000, 50000)
      if (nzchar(Sys.getenv("PS1_BURNIN_N"))) {
            N_vec_to_benchmark <- as.numeric(strsplit(x = Sys.getenv("PS1_BURNIN_N"), split = ",", fixed = TRUE)[[1]])
      }
      ##
      if (run_as_smoke_test) {
            burnin_benchmark_settings$n_chains_burnin_vec <- 4
            burnin_benchmark_settings$n_threads_WCP_burnin_vec_given_N <- list("500" = list("4" = 1))
            burnin_benchmark_settings$num_chunks_burnin_vec_given_N    <- list("500" = 10)
            burnin_benchmark_settings$n_timed_iters_given_N            <- list("500" = 2)
            burnin_benchmark_settings$n_runs <- 1
            N_vec_to_benchmark <- 500
            configuration_set_to_run <- "chunk_grid"
      }
      ##
      message(NicoStan::colourise(text = paste0("PS1 burn-in (Stan tape-chunked) on ", computer, ": chains ",
                                                paste(burnin_benchmark_settings$n_chains_burnin_vec, collapse = "/"),
                                                ", L = ", burnin_benchmark_settings$L_main, ", eps = ", burnin_benchmark_settings$eps_main,
                                                ", n_runs = ", burnin_benchmark_settings$n_runs, ", thread limit = ",
                                                burnin_benchmark_settings$n_total_threads, ", configuration set = ", configuration_set_to_run,
                                                ", Stan model loading = ", burnin_benchmark_settings$Stan_model_loading,
                                                ", NicoStan.so md5 = ", burnin_benchmark_settings$NicoStan_library_md5,
                                                if (run_as_smoke_test) " (SMOKE TEST)" else ""),
                                  fg = "cyan"))
      for (N_key in as.character(N_vec_to_benchmark)) {
            message(NicoStan::colourise(text = paste0("  N = ", N_key, ": Stan chunks ",
                                                      paste(burnin_benchmark_settings$num_chunks_burnin_vec_given_N[[N_key]], collapse = ", "),
                                                      " | BayesMVP burn-in chunks ", paste(BayesMVP_burnin_chunks_given_N[[N_key]], collapse = ", "),
                                                      " | timed iterations ", burnin_benchmark_settings$n_timed_iters_given_N[[N_key]]),
                                        fg = "cyan"))
      }
}
##
## ------- Configuration grid (settings only; no model initialisation or MCMC):  --------------------------------------------------------
##
{
      grid_outputs <- R_fn_ps1_burnin_make_configuration_grid( burnin_benchmark_settings = burnin_benchmark_settings,
                                                               N_vec                     = N_vec_to_benchmark)
      configuration_grid <- R_fn_ps1_burnin_add_configuration_sets( configuration_grid        = grid_outputs$configuration_grid,
                                                                    burnin_benchmark_settings = burnin_benchmark_settings)
      if (configuration_set_to_run != "all") {
            configuration_grid <- configuration_grid[configuration_grid$configuration_set == configuration_set_to_run, , drop = FALSE]
      }
      ##
      distinct_configurations <- unique(configuration_grid[, c("N", "configuration_set", "n_chains_burnin",
                                                               "num_chunks_burnin", "n_threads_WCP_burnin")])
      ##
      ## ---- a selection can be empty (e.g. no WCP-only configuration at N = 500 on the laptop, where every WCP count is a chunk count):
      ##
      if (nrow(distinct_configurations) > 0) {
            configuration_counts <- stats::aggregate( x   = list(n_configurations = rep(1, nrow(distinct_configurations))),
                                                      by  = distinct_configurations[, c("N", "configuration_set", "n_chains_burnin")],
                                                      FUN = sum)
            print(x = configuration_counts[order(configuration_counts$N, configuration_counts$configuration_set,
                                                 configuration_counts$n_chains_burnin), ], row.names = FALSE)
      }
      message(NicoStan::colourise(text = paste0("Selected: ", nrow(distinct_configurations), " configurations, ",
                                                nrow(configuration_grid), " configuration/runs."),
                                  fg = "cyan"))
}
##
## ------- Data, compiled model and run:  -----------------------------------------------------------------------------------------------
##
if (!run_as_dry_run && nrow(configuration_grid) > 0) {
      ##
      ## ---- the same datasets as the BayesMVP burn-in study (COVID-19 DGM, parallel-scaling N's, seed 123):
      ##
      simulator_environment <- new.env(parent = globalenv())
      sys.source(file = file.path(algorithm_study_dir, "0_utilities", "shared_functions", "R_fn_sim_bin_COVID_19_LC_MVP_data.R"),
                 envir = simulator_environment)
      true_vals_list <- simulator_environment$R_fn_simulate_binary_LC_MVP_data_COVID_19( study_type              = "algorithm_parallel_scaling_tests",
                                                                                        seed                    = 123,
                                                                                        corr_force_positive_DGM = FALSE)
      N_vec   <- true_vals_list$N_vec
      y_list  <- true_vals_list$y_list
      ##
      stan_data_by_N <- list()
      for (N in N_vec_to_benchmark) {
            stan_data_by_N[[as.character(N)]] <- paper1_sampling_helpers$fn_paper1_stan_data(y = y_list[[match(N, N_vec)]])
      }
      ##
      ## ---- the compiled Stan tape-chunked model of the sampling study (machine cache; recompiled only if its fingerprint changed):
      ##
      stan_models <- paper1_sampling_helpers$fn_paper1_compile_stan_via_NicoStan(
            settings   = list( algorithm_study_dir  = algorithm_study_dir,
                               device               = if (computer == "Local_HPC") "HPC" else "Laptop",
                               stan_model_cache_dir = file.path(ps1_dir, "stan_model_cache", Sys.info()[["nodename"]]),
                               force_recompile      = FALSE,
                               stan_via_NicoStan    = list( step_size         = 0.00001,
                                                            fixed_L           = 1,
                                                            randomize_tau     = FALSE,
                                                            metric_shape_main = "diag",
                                                            init_radius       = 0.01)),
            algorithms = "AD_Stan_tape_chunked")
      ##
      burnin_benchmark_results <- R_fn_run_ps_1_burnin_benchmark_Stan_tape_chunked(
            y_list                    = y_list[match(N_vec_to_benchmark, N_vec)],
            N_vec                     = N_vec_to_benchmark,
            burnin_benchmark_settings = burnin_benchmark_settings,
            configuration_grid        = configuration_grid,
            stan_data_by_N            = stan_data_by_N,
            model_specification       = stan_models[["AD_Stan_tape_chunked"]],
            output_dir                = burnin_output_dir,
            json_dir                  = json_dir)
}
##
message(NicoStan::colourise(text = "PS1_BURNIN_STAN_TAPE_CHUNKED_DONE", fg = "green"))
























