##
## R_fn_run_ps_1_burnin_benchmark_Stan_tape_chunked.R
##
## -| ------------------------------    Pilot study 1 (burn-in) - the Stan tape-chunked model, timed like BayesMVP ---------------------
##
## Extension of R_fn_run_ps_1_burnin_benchmark.R (which is NOT edited): the same burn-in timing, applied to the Stan
## tape-chunked model of the Paper 1 sampling study (LC_MVP_bin_PartialLog_v5_reduce_sum_static.stan, the model of the
## AD_Stan_tape_chunked and AD_Stan_WCP arms), run by NicoStan's persistent burn-in worker through BridgeStan.
##
## What is timed is exactly what R_fn_ps1_time_burnin_iters() times for BayesMVP: the persistent burn-in worker
## (fn_create_persistent_burnin_worker + fn_persistent_burnin_run_one_iter), all burn-in chains advanced one joint
## diffusion-HMC iteration in parallel, with
##   - eps = 1e-5 and L = 20 (tau = L * eps), tau_ii ~ U(0, 2 tau), shared across chains, same seeds per run;
##   - the same burn-in chain counts, WCP thread counts, timed iterations and repeats as the BayesMVP burn-in study;
##   - num_chunks = the number of tape chunks (reduce_sum_static partial sums): chunk_size = ceiling(N / num_chunks),
##     exactly as fn_paper1_stan_partition_columns() sets it for the sampling study.
##
## WCP for the Stan model: the within-chain threads are reduce_sum threads. They come from the shared TBB pool,
## which is set to n_chains_burnin x n_threads_WCP_burnin threads - as the sampling study does for AD_Stan_WCP /
## AD_Stan_tape_chunked (RcppParallel::setThreadOptions(numThreads = n_chains x threads_per_chain)), and as NicoStan's
## own burn-in does for Stan models (burnin_TBB_pool_equals_n_chains is forced to FALSE there). n_threads_WCP_burnin = 1
## is the AD_Stan_tape_chunked arm (tape chunking, one thread per chain); n_threads_WCP_burnin > 1 is the AD_Stan_WCP arm
## (the same model and chunks, plus reduce_sum threads). For Stan the n_threads_WCP argument of the worker only sizes the
## (unused) built-in-model workspaces; the thread count is set by the TBB pool.
##
## Each burn-in iteration of NicoStan's worker loads the Stan model and its JSON data once per chain (the same C++ path
## the real NicoStan burn-in uses for Stan models), so that cost is inside the timed iteration, as it is in a real burn-in.
##
## The compiled model, its cache and the Stan data list are those of the sampling study: this file uses
## fn_paper1_compile_stan_via_NicoStan(), fn_paper1_stan_data() and fn_paper1_COVID_data() from
## R_fns_alg_paper_1_chunking_WCP_par_scaling.R (sourced into its own environment; not edited).
##



##
## ---- Read one list assignment (e.g. stan_chunk_candidates) from a runner file WITHOUT running the runner:
##      every "name <- value" assignment is found in the parsed file (inside braces and if/else blocks too) and the LAST one
##      is evaluated in an empty environment, so the grids always follow the runner's current (uncommented) settings.
##
R_fn_ps1_burnin_read_assignment_from_runner <- function( runner_file_path,
                                                          object_name
) {

        if (!file.exists(runner_file_path)) stop("Runner file not found: ", runner_file_path)
        parsed_runner <- parse(file = runner_file_path, keep.source = FALSE)
        found_values <- list()
        ##
        collect_assignments <- function(expression_to_search) {
              if (is.call(expression_to_search)) {
                    call_name <- as.character(expression_to_search[[1]])[1]
                    if (call_name %in% c("<-", "=") && length(expression_to_search) == 3 &&
                        is.name(expression_to_search[[2]]) && identical(as.character(expression_to_search[[2]]), object_name)) {
                          found_values[[length(found_values) + 1]] <<- expression_to_search[[3]]
                    }
                    for (argument_index in seq_along(expression_to_search)[-1]) {
                          argument_expression <- expression_to_search[[argument_index]]
                          if (!missing(argument_expression)) collect_assignments(argument_expression)
                    }
              }
              invisible(NULL)
        }
        for (top_level_expression in as.list(parsed_runner)) collect_assignments(top_level_expression)
        ##
        if (length(found_values) == 0) stop("No assignment to '", object_name, "' found in ", runner_file_path)
        return(eval(expr = found_values[[length(found_values)]], envir = new.env(parent = baseenv())))

}
##
## ---- Read the BayesMVP burn-in settings block (burnin_benchmark_settings) of ps_1_burnin_optimizing_N_chunks_and_WCP_threads.R
##      without running the runner: the top-level brace block that creates burnin_benchmark_settings is evaluated with the
##      runner's own device name ("Local_HPC" or "Laptop"). It only assigns lists and vectors.
##
R_fn_ps1_burnin_read_BayesMVP_burnin_settings <- function( burnin_runner_file_path,
                                                            computer
) {

        if (!file.exists(burnin_runner_file_path)) stop("Burn-in runner not found: ", burnin_runner_file_path)
        parsed_runner <- parse(file = burnin_runner_file_path, keep.source = FALSE)
        settings_block_index <- which(vapply(X = as.list(parsed_runner), FUN = function(top_level_expression) {
              is.call(top_level_expression) && identical(as.character(top_level_expression[[1]]), "{") &&
                  any(grepl(pattern = "burnin_benchmark_settings <- list()", x = deparse(top_level_expression), fixed = TRUE))
        }, FUN.VALUE = logical(1)))
        if (length(settings_block_index) != 1) {
              stop("Expected exactly one settings block creating burnin_benchmark_settings in ", burnin_runner_file_path,
                   "; found ", length(settings_block_index), ".")
        }
        settings_environment <- new.env(parent = baseenv())
        assign(x = "computer", value = computer, envir = settings_environment)
        eval(expr = parsed_runner[[settings_block_index]], envir = settings_environment)
        ##
        return(get(x = "burnin_benchmark_settings", envir = settings_environment))

}
##
## ---- Add the WCP-only configurations of the sampling runner (include_NicoStan_WCP_only / include_BayesMVP_WCP_only):
##      one chunk per WCP thread (num_chunks = n_threads_WCP), for every WCP count > 1, where that chunk count is not in the
##      chunk grid already. Returned rows carry configuration_set = "WCP_only"; chunk-grid rows carry "chunk_grid".
##
R_fn_ps1_burnin_add_configuration_sets <- function( configuration_grid,
                                                    burnin_benchmark_settings
) {

        configuration_grid$configuration_set <- "chunk_grid"
        WCP_only_rows <- list()
        ##
        for (N in unique(configuration_grid$N)) {
              N_key <- as.character(N)
              chunk_grid_for_N <- burnin_benchmark_settings[["num_chunks_burnin_vec_given_N"]][[N_key]]
              for (n_chains_burnin in burnin_benchmark_settings[["n_chains_burnin_vec"]]) {
                    WCP_vec <- burnin_benchmark_settings[["n_threads_WCP_burnin_vec_given_N"]][[N_key]][[as.character(n_chains_burnin)]]
                    WCP_vec <- WCP_vec[WCP_vec > 1 & !(WCP_vec %in% chunk_grid_for_N)]
                    WCP_vec <- WCP_vec[n_chains_burnin * WCP_vec <= burnin_benchmark_settings[["n_total_threads"]]]
                    if (length(WCP_vec) == 0) next
                    WCP_only_rows[[length(WCP_only_rows) + 1]] <- dplyr::mutate(
                                      .data = tidyr::expand_grid( run_number           = seq_len(burnin_benchmark_settings$n_runs),
                                                                  n_threads_WCP_burnin = WCP_vec),
                                      n_chains_burnin   = n_chains_burnin,
                                      num_chunks_burnin = .data$n_threads_WCP_burnin,
                                      N                 = N,
                                      n_threads_total   = n_chains_burnin * .data$n_threads_WCP_burnin,
                                      configuration_set = "WCP_only")
              }
        }
        ##
        all_rows <- dplyr::bind_rows(configuration_grid, dplyr::bind_rows(WCP_only_rows))
        all_rows <- dplyr::select(.data = all_rows, "N", "n_chains_burnin", "num_chunks_burnin", "n_threads_WCP_burnin",
                                  "n_threads_total", "run_number", "configuration_set")
        ##
        return(dplyr::arrange(.data = all_rows, .data$N, .data$configuration_set, .data$run_number))

}
##
## ---- Rows of a saved result file that already complete a configuration/run with the same timing settings (resume):
##
R_fn_ps1_burnin_completed_keys <- function( results_file_path,
                                            L_main,
                                            n_timed_iters
) {

        key_columns <- c("N", "n_chains_burnin", "num_chunks_burnin", "n_threads_WCP_burnin", "run_number")
        if (!file.exists(results_file_path)) return(character(0))
        saved_results <- readRDS(file = results_file_path)
        if (!is.data.frame(saved_results) || nrow(saved_results) == 0) return(character(0))
        completed_rows <- saved_results[is.finite(saved_results$total_timed_seconds) & saved_results$total_timed_seconds > 0 &
                                        saved_results$L_main == L_main & saved_results$n_timed_iters == n_timed_iters, , drop = FALSE]
        ##
        return(do.call(what = paste, args = c(unname(as.list(completed_rows[key_columns])), sep = "_")))

}
##
R_fn_ps1_burnin_configuration_keys <- function( configuration_grid ) {

        return(paste(configuration_grid$N, configuration_grid$n_chains_burnin, configuration_grid$num_chunks_burnin,
                     configuration_grid$n_threads_WCP_burnin, configuration_grid$run_number, sep = "_"))

}
##
## ---- Save by writing a temporary file and renaming it, so an interrupted save never leaves a truncated .rds:
##
R_fn_ps1_burnin_save_atomically <- function( object_to_save,
                                             results_file_path
) {

        temporary_file_path <- paste0(results_file_path, ".tmp_", Sys.getpid())
        saveRDS(object = object_to_save, file = temporary_file_path)
        if (!file.rename(from = temporary_file_path, to = results_file_path)) {
              stop("Could not rename ", temporary_file_path, " to ", results_file_path)
        }
        invisible(results_file_path)

}




##
## ---- The Stan tape-chunked model for one N and one chunk count (chunk_size is Stan DATA, so every chunk count needs its own
##      initialisation, as in fn_paper1_run_stan_via_NicoStan()). The JSON is copied to a private file per (N, chunk count),
##      because the model-cache JSON of a later initialisation may overwrite the earlier one.
##
R_fn_ps1_burnin_Stan_make_model_inputs <- function( y,
                                                    num_chunks_burnin,
                                                    stan_data_without_chunk_size,
                                                    model_specification,
                                                    json_dir,
                                                    MCMC_seed = 1000
) {

        N <- nrow(y)
        stan_data <- stan_data_without_chunk_size
        stan_data$chunk_size <- as.integer(ceiling(N / num_chunks_burnin))   ## = fn_paper1_stan_partition_columns()
        ##
        initialised_model <- NicoStan:::initialise_model( Model_type           = "Stan",
                                                          stream               = MCMC_seed,
                                                          sample_nuisance      = TRUE,
                                                          n_nuisance_override  = nrow(y) * ncol(y),
                                                          model_args_list      = list(y = y),
                                                          Stan_data_list       = stan_data,
                                                          Stan_model_file_path = model_specification$stan_source_path,
                                                          stanc_args           = model_specification$compile_stanc_args,
                                                          make_args            = if (is.null(model_specification$compile_make_args))
                                                                                     model_specification$make_args else
                                                                                     model_specification$compile_make_args)
        ##
        model_info <- initialised_model$bs_model$model_info()
        if (!any(grepl(pattern = "STAN_THREADS[[:space:]]*=[[:space:]]*true", x = model_info, ignore.case = TRUE))) {
              stop("The Stan tape-chunked model must be compiled with STAN_THREADS=true.")
        }
        if (initialised_model$n_nuisance != nrow(y) * ncol(y) ||
            initialised_model$n_params_main != 2 * choose(n = ncol(y), k = 2) + 2 * ncol(y) + 1) {
              stop("The Stan evaluator dimensions do not match the Paper 1 model.")
        }
        if (!all(startsWith(x = initialised_model$bs_main_param_names[seq_len(initialised_model$n_nuisance)], prefix = "u_raw."))) {
              stop("The Paper 1 nuisance block must be the first Stan parameter declaration.")
        }
        ordered_beta_indices <- match(x = c("beta_vec.1", paste0("beta_vec.", ncol(y) + 1)),
                                      table = initialised_model$bs_main_param_names) - initialised_model$n_nuisance
        if (anyNA(ordered_beta_indices)) stop("Could not locate the Stan class-order coordinates.")
        ##
        dir.create(path = json_dir, recursive = TRUE, showWarnings = FALSE)
        json_file_path <- file.path(json_dir, paste0("Stan_tape_chunked_N_", N, "_num_chunks_", num_chunks_burnin,
                                                     "_chunk_size_", stan_data$chunk_size, ".json"))
        if (!file.copy(from = initialised_model$json_file_path, to = json_file_path, overwrite = TRUE)) {
              stop("Could not copy the Stan data JSON to ", json_file_path)
        }
        ##
        Model_args_as_Rcpp_List <- initialised_model$Model_args_as_Rcpp_List
        Model_args_as_Rcpp_List$json_file_path <- normalizePath(path = json_file_path, mustWork = TRUE)
        Model_args_as_Rcpp_List$N <- as.integer(N)
        Model_args_as_Rcpp_List$n_tests <- as.integer(ncol(y))
        ##
        return(list( y                        = y,
                     N                        = N,
                     num_chunks_burnin        = num_chunks_burnin,
                     stan_chunk_size          = stan_data$chunk_size,
                     n_params_main            = initialised_model$n_params_main,
                     n_nuisance               = initialised_model$n_nuisance,
                     ordered_beta_indices     = ordered_beta_indices,
                     Model_args_as_Rcpp_List  = Model_args_as_Rcpp_List))

}
##
## ---- Time the burn-in iteration for ONE (n_chains_burnin, num_chunks, n_threads_WCP) configuration of the Stan model.
##      Mirrors R_fn_ps1_time_burnin_iters(); the differences are the Stan model, its inits and the TBB pool (see the header).
##
R_fn_ps1_time_burnin_iters_Stan <- function( model_inputs,
                                             n_chains_burnin,
                                             n_threads_WCP_burnin,
                                             L_main = 20,
                                             eps_main = 1e-5,
                                             n_untimed_burnin_iter_before_timing = 0,
                                             n_timed_iters = 10,
                                             seed = 1,
                                             n_total_threads = NULL,
                                             share_tau_ii_across_chains_in_burnin = TRUE,
                                             init_radius = 0.01
) {

        thread_limit <- R_fn_ps1_burnin_thread_limit(n_total_threads)
        if (n_chains_burnin * n_threads_WCP_burnin > thread_limit) {
              stop("Requested ", n_chains_burnin, " burn-in chains x ", n_threads_WCP_burnin, " WCP threads = ",
                   n_chains_burnin * n_threads_WCP_burnin, " threads, exceeding the thread limit of ", thread_limit, ".")
        }
        time_worker_setup_start <- proc.time()[["elapsed"]]
        n_params_main <- model_inputs$n_params_main
        n_nuisance    <- model_inputs$n_nuisance
        ##
        ## ---- eps tiny, trajectory length L leapfrog steps on average (tau_ii ~ U(0, 2 tau) as in burn-in), as for BayesMVP:
        ##
        EHMC_args_as_Rcpp_List <- NicoStan:::init_EHMC_args_as_Rcpp_List( diffusion_HMC            = TRUE,
                                                                          diffusion_HMC_integrator = "kick_flow_kick")
        EHMC_args_as_Rcpp_List$share_tau_ii_across_chains <- share_tau_ii_across_chains_in_burnin
        EHMC_args_as_Rcpp_List$eps_main <- eps_main
        EHMC_args_as_Rcpp_List$tau_main <- L_main * eps_main
        EHMC_args_as_Rcpp_List$eps_us   <- eps_main
        EHMC_args_as_Rcpp_List$tau_us   <- L_main * eps_main
        ##
        ## ---- Metric: the same default dense-main metric initialisation as the BayesMVP burn-in benchmark:
        ##
        EHMC_Metric_as_Rcpp_List <- NicoStan:::init_EHMC_Metric_as_Rcpp_List( n_params_main     = n_params_main,
                                                                              n_nuisance        = n_nuisance,
                                                                              metric_shape_main = "dense")
        ##
        ## ---- Inits: the Stan arms' inits of the sampling study (uniform on (-init_radius, init_radius) on the unconstrained
        ##      scale, seeded; the two class-order beta coordinates sorted), because the BayesMVP inits are in BayesMVP's own
        ##      parameterisation:
        ##
        initial_values <- withr::with_seed(seed = seed, code = {
              list( main     = matrix(data = stats::runif(n = n_params_main * n_chains_burnin, min = -init_radius, max = init_radius),
                                      nrow = n_params_main, ncol = n_chains_burnin),
                    nuisance = matrix(data = stats::runif(n = n_nuisance * n_chains_burnin, min = -init_radius, max = init_radius),
                                      nrow = n_nuisance, ncol = n_chains_burnin))
        })
        beta_first  <- initial_values$main[model_inputs$ordered_beta_indices[1], ]
        beta_second <- initial_values$main[model_inputs$ordered_beta_indices[2], ]
        initial_values$main[model_inputs$ordered_beta_indices[1], ] <- pmin(beta_first, beta_second)
        initial_values$main[model_inputs$ordered_beta_indices[2], ] <- pmax(beta_first, beta_second)
        ##
        ## ---- TBB pool of n_chains_burnin x n_threads_WCP_burnin threads: chains and reduce_sum share it (see the header):
        ##
        TBB_pool_threads <- n_chains_burnin * n_threads_WCP_burnin
        RcppParallel::setThreadOptions(numThreads = TBB_pool_threads)
        ##
        burnin_worker_pointer <- NicoStan:::fn_create_persistent_burnin_worker( n_threads_R              = n_chains_burnin,
                                                                                partitioned_HMC_R        = FALSE,
                                                                                diffusion_HMC_R          = TRUE,
                                                                                Model_type_R             = "Stan",
                                                                                sample_nuisance_R        = TRUE,
                                                                                force_autodiff_R         = TRUE,
                                                                                force_PartialLog_R       = FALSE,
                                                                                multi_attempts_R         = FALSE,
                                                                                y_Eigen_R                = model_inputs$y,
                                                                                Model_args_as_Rcpp_List  = model_inputs$Model_args_as_Rcpp_List,
                                                                                EHMC_args_as_Rcpp_List   = EHMC_args_as_Rcpp_List,
                                                                                EHMC_Metric_as_Rcpp_List = EHMC_Metric_as_Rcpp_List,
                                                                                n_threads_WCP            = n_threads_WCP_burnin)
        NicoStan:::fn_persistent_burnin_update_adaptation( burnin_worker_pointer,
                                                           EHMC_args_as_Rcpp_List,
                                                           EHMC_Metric_as_Rcpp_List)
        NicoStan:::fn_persistent_burnin_set_theta( burnin_worker_pointer,
                                                   initial_values$main,
                                                   initial_values$nuisance)
        time_worker_setup_seconds <- proc.time()[["elapsed"]] - time_worker_setup_start
        ##
        ## ---- Optional untimed iterations, then the timed iterations (worker construction above is outside the timer):
        ##
        sec_per_iter <- rep(NA_real_, n_timed_iters)
        n_divs <- 0
        ##
        for (iter_index in seq_len(n_untimed_burnin_iter_before_timing + n_timed_iters)) {

                  iter_start_time <- proc.time()[["elapsed"]]
                  burnin_iter_outputs <- NicoStan:::fn_persistent_burnin_run_one_iter( burnin_worker_pointer,
                                                                                       seed + iter_index,
                                                                                       iter_index)
                  iter_wall_time <- proc.time()[["elapsed"]] - iter_start_time
                  ##
                  if (iter_index > n_untimed_burnin_iter_before_timing) {

                          sec_per_iter[iter_index - n_untimed_burnin_iter_before_timing] <- iter_wall_time
                          n_divs_this_iter <- sum(burnin_iter_outputs$other_main_out_vector_all_chains_output_to_R[2, ])
                          n_divs <- n_divs + n_divs_this_iter
                          ##
                          if (n_divs_this_iter > 0) {
                                stop("Burn-in PS1 (Stan) encountered ", n_divs_this_iter, " divergent transition(s) in timed iter ",
                                     iter_index - n_untimed_burnin_iter_before_timing,
                                     ". The configuration is invalid for a timing comparison.")
                          }

                  }

        }
        rm(burnin_worker_pointer) ; invisible(gc(verbose = FALSE))
        ##
        return(list( sec_per_iter              = sec_per_iter,
                     total_timed_seconds       = sum(sec_per_iter),
                     n_divs                    = n_divs,
                     TBB_pool_threads          = TBB_pool_threads,
                     time_worker_setup_seconds = time_worker_setup_seconds))

}
##
## ---- Run the Stan tape-chunked burn-in grid for the selected N's; one .rds per N, saved after EVERY configuration/run, and
##      configurations/runs already in that file (same L and timed iterations) are skipped, so a restarted run resumes.
##
R_fn_run_ps_1_burnin_benchmark_Stan_tape_chunked <- function( y_list,
                                                              N_vec,
                                                              burnin_benchmark_settings,
                                                              configuration_grid,
                                                              stan_data_by_N,
                                                              model_specification,
                                                              output_dir,
                                                              json_dir
) {

        if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
        device_prefix <- if (burnin_benchmark_settings$device == "Laptop") "Laptop_" else "HPC_"
        all_N_results <- list()
        ##
        for (dataset_position in seq_along(N_vec)) {

              N <- N_vec[dataset_position]
              y_for_this_N <- y_list[[dataset_position]]
              if (!is.matrix(y_for_this_N) || nrow(y_for_this_N) != N) stop("Dataset ", dataset_position, " does not have N = ", N, " rows.")
              n_timed_iters <- burnin_benchmark_settings[["n_timed_iters_given_N"]][[as.character(N)]]
              results_file_path <- file.path(output_dir, paste0(device_prefix, "ps1_burnin_benchmark_Stan_tape_chunked_N", N,
                                                                "_L", burnin_benchmark_settings$L_main,
                                                                "_n_runs", burnin_benchmark_settings$n_runs, ".rds"))
              ##
              configuration_grid_for_N <- dplyr::filter(.data = configuration_grid, .data$N == .env$N)
              completed_keys <- R_fn_ps1_burnin_completed_keys( results_file_path = results_file_path,
                                                                L_main            = burnin_benchmark_settings$L_main,
                                                                n_timed_iters     = n_timed_iters)
              is_completed <- R_fn_ps1_burnin_configuration_keys(configuration_grid_for_N) %in% completed_keys
              message(NicoStan::colourise(text = paste0("Stan tape-chunked burn-in, N = ", N, ": ", nrow(configuration_grid_for_N),
                                                        " configuration/run(s) selected, ", sum(is_completed),
                                                        " already saved (skipped), ", sum(!is_completed), " to run. File: ",
                                                        results_file_path),
                                          fg = "cyan"))
              configuration_grid_for_N <- configuration_grid_for_N[!is_completed, , drop = FALSE]
              if (nrow(configuration_grid_for_N) == 0) next
              ##
              ## ---- One initialisation per chunk count (outside every timer), reused by all chains/WCP/runs:
              ##
              model_inputs_by_chunks <- list()
              for (num_chunks_burnin in sort(unique(configuration_grid_for_N$num_chunks_burnin))) {
                    time_model_initialisation_start <- proc.time()[["elapsed"]]
                    model_inputs_by_chunks[[as.character(num_chunks_burnin)]] <- R_fn_ps1_burnin_Stan_make_model_inputs(
                          y                            = y_for_this_N,
                          num_chunks_burnin            = num_chunks_burnin,
                          stan_data_without_chunk_size = stan_data_by_N[[as.character(N)]],
                          model_specification          = model_specification,
                          json_dir                     = json_dir)
                    message(NicoStan::colourise(text = paste0("Initialised the Stan tape-chunked model: N = ", N, ", num_chunks = ",
                                                              num_chunks_burnin, ", chunk_size = ",
                                                              model_inputs_by_chunks[[as.character(num_chunks_burnin)]]$stan_chunk_size,
                                                              " (", formatC(x = proc.time()[["elapsed"]] - time_model_initialisation_start,
                                                                            format = "f", digits = 1), " s)"),
                                                fg = "cyan"))
              }
              ##
              for (configuration_index in seq_len(nrow(configuration_grid_for_N))) {

                    configuration <- dplyr::slice(.data = configuration_grid_for_N, configuration_index)
                    model_inputs <- model_inputs_by_chunks[[as.character(configuration$num_chunks_burnin)]]
                    timing_outputs <- R_fn_ps1_time_burnin_iters_Stan( model_inputs          = model_inputs,
                                                                       n_chains_burnin       = configuration$n_chains_burnin,
                                                                       n_threads_WCP_burnin  = configuration$n_threads_WCP_burnin,
                                                                       L_main                = burnin_benchmark_settings$L_main,
                                                                       eps_main              = burnin_benchmark_settings$eps_main,
                                                                       n_untimed_burnin_iter_before_timing = burnin_benchmark_settings$n_untimed_burnin_iter_before_timing,
                                                                       n_timed_iters         = n_timed_iters,
                                                                       seed                  = 1000 * configuration$run_number,
                                                                       n_total_threads       = burnin_benchmark_settings$n_total_threads,
                                                                       share_tau_ii_across_chains_in_burnin = burnin_benchmark_settings$share_tau_ii_across_chains_in_burnin,
                                                                       init_radius           = burnin_benchmark_settings$stan_init_radius)
                    ##
                    new_row <- tibble::tibble(
                          device                      = burnin_benchmark_settings$device,
                          N                           = N,
                          n_chains_burnin             = configuration$n_chains_burnin,
                          num_chunks_burnin           = configuration$num_chunks_burnin,
                          n_threads_WCP_burnin        = configuration$n_threads_WCP_burnin,
                          n_threads_total             = configuration$n_chains_burnin * configuration$n_threads_WCP_burnin,
                          run_number                  = configuration$run_number,
                          L_main                      = burnin_benchmark_settings$L_main,
                          n_timed_iters               = n_timed_iters,
                          total_timed_seconds         = timing_outputs$total_timed_seconds,
                          median_sec_per_iter         = stats::median(timing_outputs$sec_per_iter),
                          mean_sec_per_iter           = mean(timing_outputs$sec_per_iter),
                          n_divs                      = timing_outputs$n_divs,
                          ##
                          implementation              = "Stan_tape_chunked",
                          configuration_set           = configuration$configuration_set,
                          stan_chunk_size             = model_inputs$stan_chunk_size,
                          eps_main                    = burnin_benchmark_settings$eps_main,
                          TBB_pool_threads            = timing_outputs$TBB_pool_threads,
                          time_worker_setup_seconds   = timing_outputs$time_worker_setup_seconds,
                          sec_per_iter_each           = paste(formatC(x = timing_outputs$sec_per_iter, format = "f", digits = 6), collapse = ","),
                          Stan_model_loading          = burnin_benchmark_settings$Stan_model_loading,
                          NicoStan_library_md5        = burnin_benchmark_settings$NicoStan_library_md5,
                          time_completed              = format(x = Sys.time(), format = "%Y-%m-%d %H:%M:%S"))
                    ##
                    saved_results <- if (file.exists(results_file_path)) readRDS(file = results_file_path) else NULL
                    R_fn_ps1_burnin_save_atomically( object_to_save    = dplyr::bind_rows(saved_results, new_row),
                                                     results_file_path = results_file_path)
                    ##
                    progress_message <- paste0( "Stan tape-chunked | N = ", formatC(x = N, format = "d", width = 6),
                                                " | chains ", formatC(x = configuration$n_chains_burnin, format = "d", width = 2),
                                                " | chunks ", formatC(x = configuration$num_chunks_burnin, format = "d", width = 4),
                                                " | n_WCP (threads/chain) ", formatC(x = configuration$n_threads_WCP_burnin, format = "d", width = 3),
                                                " | total threads ", formatC(x = new_row$n_threads_total, format = "d", width = 3),
                                                " | run ", configuration$run_number,
                                                " | ", formatC(x = timing_outputs$total_timed_seconds, format = "f", digits = 4),
                                                " s total | ", formatC(x = new_row$mean_sec_per_iter, format = "f", digits = 4),
                                                " s / iter (mean) | setup ", formatC(x = timing_outputs$time_worker_setup_seconds,
                                                                                     format = "f", digits = 2), " s",
                                                " | ", configuration_index, "/", nrow(configuration_grid_for_N))
                    message(NicoStan::colourise(text = progress_message, fg = "cyan"))

              }
              ##
              all_N_results[[as.character(N)]] <- readRDS(file = results_file_path)
              message(NicoStan::colourise(text = paste0("saved: ", results_file_path), fg = "green"))

        }
        ##
        return(dplyr::bind_rows(all_N_results))

}
##
## ---- BayesMVP top-up file of a canonical result file (the report reads only the canonical file names, not this one):
##
R_fn_ps1_burnin_top_up_file_path <- function( results_file_path ) {

        return(sub(pattern = "\\.rds$", replacement = "_top_up_missing_configurations.rds", x = results_file_path))

}
##
## ---- Merge the top-up rows that the canonical file does not hold yet into it. The merge is written only after checking that
##      every existing row is preserved unchanged (same count, identical values) and that no configuration/run is duplicated.
##
R_fn_ps1_burnin_merge_top_up_into_canonical <- function( results_file_path ) {

        top_up_file_path <- R_fn_ps1_burnin_top_up_file_path(results_file_path)
        if (!file.exists(top_up_file_path)) return(invisible(0))
        top_up_rows <- readRDS(file = top_up_file_path)
        old_rows <- if (file.exists(results_file_path)) readRDS(file = results_file_path) else top_up_rows[0, , drop = FALSE]
        key_columns <- c("N", "n_chains_burnin", "num_chunks_burnin", "n_threads_WCP_burnin", "run_number")
        old_keys    <- do.call(what = paste, args = c(unname(as.list(old_rows[key_columns])), sep = "_"))
        top_up_keys <- do.call(what = paste, args = c(unname(as.list(top_up_rows[key_columns])), sep = "_"))
        rows_to_add <- top_up_rows[!(top_up_keys %in% old_keys), intersect(names(old_rows), names(top_up_rows)), drop = FALSE]
        if (nrow(rows_to_add) == 0) return(invisible(0))
        ##
        merged_rows <- dplyr::bind_rows(old_rows, rows_to_add)
        merged_keys <- do.call(what = paste, args = c(unname(as.list(merged_rows[key_columns])), sep = "_"))
        old_rows_preserved <- nrow(merged_rows) == nrow(old_rows) + nrow(rows_to_add) &&
                              identical(names(merged_rows), names(old_rows)) &&
                              identical(as.data.frame(merged_rows[seq_len(nrow(old_rows)), , drop = FALSE]), as.data.frame(old_rows))
        if (!old_rows_preserved || anyDuplicated(merged_keys) > 0) {
              stop("Merge check failed for ", results_file_path, ": the existing ", nrow(old_rows),
                   " row(s) would not be preserved unchanged, or a configuration/run would be duplicated. Nothing was written.")
        }
        R_fn_ps1_burnin_save_atomically( object_to_save    = merged_rows,
                                         results_file_path = results_file_path)
        message(NicoStan::colourise(text = paste0("merged ", nrow(rows_to_add), " top-up row(s) into ", basename(results_file_path),
                                                  " (", nrow(old_rows), " existing rows preserved; now ", nrow(merged_rows), ")"),
                                    fg = "green"))
        return(invisible(nrow(rows_to_add)))

}
##
## ---- BayesMVP: run ONLY the configurations missing from an existing burn-in result file. Every new row is first appended to a
##      separate top-up file (<canonical name>_top_up_missing_configurations.rds), then merged into the canonical file by
##      R_fn_ps1_burnin_merge_top_up_into_canonical(), which checks that the existing rows are preserved. The canonical file is copied
##      to backup_dir once before its first change (an existing backup is never overwritten). Uses the unchanged
##      R_fn_ps1_burnin_make_model_inputs() and R_fn_ps1_time_burnin_iters(); rows carry exactly the existing file's columns.
##
R_fn_run_ps_1_burnin_benchmark_BayesMVP_missing <- function( y,
                                                             N,
                                                             n_tests,
                                                             burnin_benchmark_settings,
                                                             missing_configuration_grid,
                                                             results_file_path,
                                                             backup_dir
) {

        if (nrow(missing_configuration_grid) == 0) {
              message(NicoStan::colourise(text = paste0("BayesMVP burn-in, N = ", N, ": nothing missing in ", results_file_path), fg = "green"))
              return(invisible(NULL))
        }
        if (file.exists(results_file_path)) {
              dir.create(path = backup_dir, recursive = TRUE, showWarnings = FALSE)
              backup_file_path <- file.path(backup_dir, basename(results_file_path))
              if (!file.exists(backup_file_path)) {
                    if (!file.copy(from = results_file_path, to = backup_file_path, overwrite = FALSE, copy.date = TRUE)) {
                          stop("Could not back up ", results_file_path, " to ", backup_file_path)
                    }
                    message(NicoStan::colourise(text = paste0("Backed up ", results_file_path, " -> ", backup_file_path), fg = "green"))
              }
        }
        top_up_file_path <- R_fn_ps1_burnin_top_up_file_path(results_file_path)
        existing_columns <- if (file.exists(results_file_path)) names(readRDS(file = results_file_path)) else NULL
        ##
        model_inputs <- R_fn_ps1_burnin_make_model_inputs( y              = y,
                                                           n_tests        = n_tests,
                                                           SIMD_vect_type = burnin_benchmark_settings$SIMD_vect_type)
        n_timed_iters <- burnin_benchmark_settings[["n_timed_iters_given_N"]][[as.character(N)]]
        missing_configuration_grid <- dplyr::arrange(.data = missing_configuration_grid, .data$run_number)
        ##
        for (configuration_index in seq_len(nrow(missing_configuration_grid))) {

              configuration <- dplyr::slice(.data = missing_configuration_grid, configuration_index)
              timing_outputs <- R_fn_ps1_time_burnin_iters( model_inputs          = model_inputs,
                                                            n_chains_burnin       = configuration$n_chains_burnin,
                                                            num_chunks_burnin     = configuration$num_chunks_burnin,
                                                            n_threads_WCP_burnin  = configuration$n_threads_WCP_burnin,
                                                            L_main                = burnin_benchmark_settings$L_main,
                                                            eps_main              = burnin_benchmark_settings$eps_main,
                                                            n_untimed_burnin_iter_before_timing = burnin_benchmark_settings$n_untimed_burnin_iter_before_timing,
                                                            n_timed_iters         = n_timed_iters,
                                                            seed                  = 1000 * configuration$run_number,
                                                            n_total_threads       = burnin_benchmark_settings$n_total_threads,
                                                            share_tau_ii_across_chains_in_burnin = burnin_benchmark_settings$share_tau_ii_across_chains_in_burnin,
                                                            burnin_TBB_pool_equals_n_chains      = burnin_benchmark_settings$burnin_TBB_pool_equals_n_chains)
              new_row <- tibble::tibble(
                    device                          = burnin_benchmark_settings$device,
                    N                               = N,
                    n_chains_burnin                 = configuration$n_chains_burnin,
                    num_chunks_burnin               = configuration$num_chunks_burnin,
                    n_threads_WCP_burnin            = configuration$n_threads_WCP_burnin,
                    n_threads_total                 = configuration$n_chains_burnin * configuration$n_threads_WCP_burnin,
                    run_number                      = configuration$run_number,
                    L_main                          = burnin_benchmark_settings$L_main,
                    n_timed_iters                   = n_timed_iters,
                    total_timed_seconds             = timing_outputs$total_timed_seconds,
                    median_sec_per_iter             = stats::median(timing_outputs$sec_per_iter),
                    mean_sec_per_iter               = mean(timing_outputs$sec_per_iter),
                    n_divs                          = timing_outputs$n_divs)
              if (!is.null(existing_columns)) new_row <- new_row[, intersect(existing_columns, names(new_row)), drop = FALSE]
              ##
              ## ---- 1. the top-up file; 2. the checked merge into the canonical file:
              ##
              saved_top_up_rows <- if (file.exists(top_up_file_path)) readRDS(file = top_up_file_path) else NULL
              R_fn_ps1_burnin_save_atomically( object_to_save    = dplyr::bind_rows(saved_top_up_rows, new_row),
                                               results_file_path = top_up_file_path)
              R_fn_ps1_burnin_merge_top_up_into_canonical(results_file_path = results_file_path)
              ##
              progress_message <- paste0( "BayesMVP (missing) | N = ", formatC(x = N, format = "d", width = 6),
                                          " | chains ", formatC(x = configuration$n_chains_burnin, format = "d", width = 2),
                                          " | chunks ", formatC(x = configuration$num_chunks_burnin, format = "d", width = 4),
                                          " | n_WCP (threads/chain) ", formatC(x = configuration$n_threads_WCP_burnin, format = "d", width = 3),
                                          " | total threads ", formatC(x = new_row$n_threads_total, format = "d", width = 3),
                                          " | run ", configuration$run_number,
                                          " | ", formatC(x = timing_outputs$total_timed_seconds, format = "f", digits = 4),
                                          " s total | ", formatC(x = new_row$mean_sec_per_iter, format = "f", digits = 4),
                                          " s / iter (mean) | ", configuration_index, "/", nrow(missing_configuration_grid))
              message(NicoStan::colourise(text = progress_message, fg = "cyan"))

        }
        message(NicoStan::colourise(text = paste0("merged into: ", results_file_path), fg = "green"))
        invisible(readRDS(file = results_file_path))

}
























