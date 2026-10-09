##
## ======================================================================================================================================
## R_fns_alg_paper_1_figures_tables.R
##
## Paper 1 reporting from the unified saved results. This file never runs sampling.
## Legacy presentation styles are supplied by the private templates factory.
## ---- Reporting dependencies ----------------------------------------------------------------------------------------------------------
##
fn_paper1_report_dependencies <-  function() {

        packages <-  c('ggplot2', 'dplyr', 'tidyr', 'patchwork', 'tibble')
        ##
        missing <-  packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
        ##
        if (length(missing)) stop('Reporting requires: ', paste(missing, collapse = ', '))
        ##
        invisible(TRUE)

}
##
## ---- Number labels: thousands separators only from 10,000 upwards (e.g. 2500, 10,000, 50,000):
##
fn_paper1_format_number_commas_from_10000 <-  function( numbers ) {

        vapply(X = numbers,
               FUN = function(one_number) {
                     if (is.finite(one_number) && abs(one_number) >= 10000) format(x = one_number, big.mark = ",", scientific = FALSE, trim = TRUE)
                     else format(x = one_number, scientific = FALSE, trim = TRUE)
               },
               FUN.VALUE = character(1),
               USE.NAMES = FALSE)

}
##
## ---- CSV export ----------------------------------------------------------------------------------------------------------------------
##
fn_paper1_write_csv <-  function( x,
                                  file) {

        for (column in names(x)[vapply(x, is.list, logical(1))]) {

            x[[column]] <-  vapply(x[[column]], paste, character(1), collapse = ';')

        }
        ##
        utils::write.csv(x, file, row.names = FALSE)

}
##
## ---- Timing columns for results saved before the two-run timing existed:
##
## Those studies timed one run per configuration, so they are labelled "single_run"; their two-run fields are NA.
##
fn_paper1_add_missing_timing_columns <-  function( results ) {

        if (!"timing_method" %in% names(x = results)) results$timing_method <-  rep(x = "single_run", times = nrow(x = results))
        ## Historical two-run files used the slope formula. Never reinterpret their saved times as overhead subtraction.
        if (!"timing_estimator" %in% names(results)) results$timing_estimator <-
            ifelse(results$timing_method == "two_run_difference", "iteration_difference", "single_run")
        if (!"mplus_iteration_mode" %in% names(results)) results$mplus_iteration_mode <- NA_character_
        if (!"mplus_iteration_mode_short_run" %in% names(results)) results$mplus_iteration_mode_short_run <-
            ifelse(results$timing_method == "two_run_difference", results$mplus_iteration_mode, NA_character_)
        ##
        numeric_timing_columns <-  c("n_iter_short_run", "elapsed_seconds_short_run", "elapsed_seconds_long_run", "fixed_cost_seconds",
                                     "seconds_per_iteration", "iteration_multiplier_long_over_short", "divergences_short_run",
                                     "mplus_iterations_per_chain", "mplus_iterations_per_chain_short_run")
        for (column_name in setdiff(x = numeric_timing_columns, y = names(x = results))) results[[column_name]] <-  rep(x = NA_real_, times = nrow(x = results))
        ##
        if (!"two_run_timing_flag" %in% names(x = results)) results$two_run_timing_flag <-  rep(x = NA_character_, times = nrow(x = results))
        for (column_name in setdiff(x = c("randomize_tau", "mplus_PPPP_run"), y = names(x = results))) results[[column_name]] <-  rep(x = NA, times = nrow(x = results))
        ##
        return(results)

}
##
## ---- One sentence that labels every figure and table with the timing method:
##
## It goes into LaTeX captions, so it contains no underscores or other characters that LaTeX treats specially.
##
fn_paper1_timing_label <-  function( timing_method, timing_estimators = NULL ) {

        return(if (identical(x = timing_method, y = "two_run_difference")) {

            if ("overhead_subtraction" %in% timing_estimators) {
                paste0("Times are startup-adjusted estimates. Mplus subtracts the complete short overhead-control call from the long call and scales to the requested iteration budget when the observed count differs; ",
                       "other implementations use the two-run slope correction where present. ",
                       "Mplus timings include its predictive-check work and are not isolated posterior-sampling times.")
            } else paste0("Times use the two-run slope correction: each configuration ran a short and a long run back to back, ",
                          "and the estimated fixed per-run cost is removed. Remaining iteration-dependent work is included.")

        } else "Times are single timed runs and include the fixed per-run start-up cost.")

}
##
## ---- Only runs that match the settings currently written in the runner enter the tables and figures ----------------------------
##
## The settings section of the runner (everything before "Compiled Stan models") is evaluated afresh in a private environment every
## time the report is made, so the report always follows the iteration counts currently written there. Nothing is sampled or compiled.
##
fn_paper1_runner_current_settings <-  function( runner_file,
                                                 device = NULL
) {

        runner_lines <-  readLines(con = runner_file, warn = FALSE)
        settings_end_line <-  grep(pattern = "^## ---- Compiled Stan models", x = runner_lines)[1]
        if (is.na(x = settings_end_line)) stop("Could not find the settings section of the runner: ", runner_file)
        ##
        runner_environment <-  new.env(parent = globalenv())
        settings_expressions <-  parse(text = runner_lines[seq_len(length.out = settings_end_line - 1)], keep.source = FALSE)
        eval(expr = settings_expressions[[1]], envir = runner_environment)
        if (!is.null(device)) {

            runner_environment$device <-  device
            runner_environment$paper1_settings$device <-  device

        }
        for (expression_index in seq_along(settings_expressions)[-1]) {

            eval(expr = settings_expressions[[expression_index]], envir = runner_environment)

        }
        ##
        return(runner_environment)

}
##
fn_paper1_rows_matching_current_settings <-  function( results,
                                                       current_settings
) {

        settings <-  get(x = "paper1_settings", envir = current_settings)
        iterations_for_case <-  get(x = "fn_paper1_iterations_for_case", envir = current_settings)
        ##
        expected_iterations <-  function(algorithm, N, which_run) {

            tryCatch( expr = iterations_for_case(settings = settings, algorithm = algorithm, N = N, which_run = which_run),
                      error = function(error) NA_real_)

        }
        ##
        long_run_matches <-  mapply( FUN = function(algorithm, N, n_iter) isTRUE(expected_iterations(algorithm, N, "long_run") == n_iter),
                                     as.character(x = results$algorithm), results$N, results$n_iter)
        ##
        short_run_matches <-  if ("n_iter_short_run" %in% names(x = results)) {

            mapply( FUN = function(algorithm, N, n_iter_short_run) {
                        is.na(x = n_iter_short_run) || isTRUE(expected_iterations(algorithm, N, "short_run") == n_iter_short_run)
                    },
                    as.character(x = results$algorithm), results$N, results$n_iter_short_run)

        } else rep(x = TRUE, times = nrow(x = results))
        ##
        make_grid <-  get(x = "fn_paper1_benchmark_grid", envir = current_settings)
        current_grid <-  make_grid(settings = settings)
        matching_columns <-  c("device", "algorithm", "execution_backend", "N", "num_chunks", "n_threads",
                               "n_chains", "threads_per_chain", "n_iter", "n_iter_short_run", "run", "seed",
                               "benchmark_role", "timing_method", "timing_estimator", "mplus_iteration_mode",
                               "mplus_iteration_mode_short_run", "stan_chunk_size")
        row_keys <-  function(rows) {

            values <-  lapply(rows[, matching_columns, drop = FALSE], function(column) {

                ifelse(is.na(column), "<missing>", as.character(column))

            })
            do.call(what = paste, args = c(values, list(sep = "\r")))

        }
        grid_matches <-  row_keys(rows = results) %in% row_keys(rows = current_grid)
        ##
        return(unname(obj = long_run_matches & short_run_matches & grid_matches))

}
##
## ---- Validate and combine completed studies ------------------------------------------------------------------------------------------
##
fn_paper1_read_studies <-  function( study_output_dirs,
                                     runner_file = path.expand('~/Documents/Work/PhD_work/Alg_paper_analysis/paper_1_chunking_and_parallel_scalability/alg_paper_1_chunking_WCP_par_scaling.R')
) {

        if (!length(study_output_dirs)) stop('Supply an explicit completed study output directory.')
        ##
        paths <-  normalizePath(study_output_dirs, mustWork = TRUE)
        ##
        if (anyDuplicated(paths)) stop('A study directory was supplied more than once.')
        ##
        current_settings <-  fn_paper1_runner_current_settings(runner_file = runner_file)
        ##
        studies <-  lapply( paths,
                            function(path) {

                files <-  file.path(path, c('results.rds', 'study_metadata.rds', 'experiment_grid.rds'))
                ##
                if (!all(file.exists(files))) stop('Missing unified study files in ', path)
                ##
                results <-  readRDS(files[1])
                ##
                metadata <-  readRDS(files[2])
                ##
                ## ---- Studies saved previously hold ONE shared Mplus long-run count, mplus_iterations, used by BOTH Mplus arms.
                ##      Since then Mplus_WCP reads mplus_WCP_iterations and Mplus_standard reads mplus_standard_iterations; the old
                ##      field is mapped onto both new names, so old and new studies are compared on what each arm actually ran.
                ##
                if (!is.null(metadata$settings$mplus_iterations)) {

                    if (is.null(metadata$settings$mplus_WCP_iterations))      metadata$settings$mplus_WCP_iterations      <-  metadata$settings$mplus_iterations
                    if (is.null(metadata$settings$mplus_standard_iterations)) metadata$settings$mplus_standard_iterations <-  metadata$settings$mplus_iterations
                    metadata$settings$mplus_iterations <-  NULL

                }
                ##
                grid <-  readRDS(files[3])
                ##
                ## ---- Keep only the runs that match the settings currently written in the runner:
                ##
                study_current_settings <-  fn_paper1_runner_current_settings(runner_file = runner_file,
                                                                              device = metadata$settings$device)
                matching_rows <-  fn_paper1_rows_matching_current_settings(results = results, current_settings = study_current_settings)
                message(paste0('\033[36m', sum(matching_rows), ' of ', length(matching_rows), ' saved runs match the current runner settings (',
                               basename(path), ').\033[0m'))
                results <-  results[matching_rows, , drop = FALSE]
                grid <-  grid[matching_rows, , drop = FALSE]
                if (!nrow(results)) return(NULL)
                ##
                if (!nrow(results) || anyNA(results$status) || any(results$status != 'completed')) {

                    stop('Reporting requires a completed study. Inspect the checkpoint: ', files[1])

                }
                ##
                if (!all(names(grid) %in% names(results)) || !identical(results[, names(grid), drop = FALSE], grid)) {

                    stop('Saved results do not match their experiment grid: ', path)

                }
                ##
                ## ---- Raw run times must be positive; a two-run corrected time may be <= 0 only when that row is flagged:
                ##
                results <-  fn_paper1_add_missing_timing_columns(results = results)
                expected_estimator <-  ifelse(results$timing_method == 'single_run', 'single_run', 'iteration_difference')
                if (identical(metadata$settings$mplus_short_run_role, 'overhead_control')) {
                    expected_estimator[grepl('^Mplus_', results$algorithm) & results$timing_method == 'two_run_difference'] <- 'overhead_subtraction'
                }
                if (anyNA(results$timing_estimator) || !identical(results$timing_estimator, expected_estimator)) {
                    invisible(NULL)
                }
                two_run_result_rows <-  results$timing_method == 'two_run_difference'
                ##
                nonpositive_corrected_rows <-  two_run_result_rows & results$elapsed_seconds <= 0
                timing_flags <-  results$two_run_timing_flag
                normalized_timing_flags <-  trimws(ifelse(is.na(timing_flags), '', timing_flags))
                unflagged_nonpositive_rows <-  nonpositive_corrected_rows &
                    (is.na(timing_flags) | normalized_timing_flags == '' | normalized_timing_flags == 'ok')
                if (any(!is.finite(results$elapsed_seconds)) ||
                    any(results$elapsed_seconds[!two_run_result_rows] <= 0) ||
                    any(!is.finite(results$elapsed_seconds_short_run[two_run_result_rows]) | results$elapsed_seconds_short_run[two_run_result_rows] <= 0) ||
                    any(!is.finite(results$elapsed_seconds_long_run[two_run_result_rows]) | results$elapsed_seconds_long_run[two_run_result_rows] <= 0) ||
                    any(unflagged_nonpositive_rows)) {

                    stop('Invalid saved elapsed times.')

                }
                ##
                if (anyNA(results$dataset_md5) || any(results$dataset_md5 != unname(metadata$data$dataset_md5[as.character(results$N)]))) {

                    stop('Result dataset hashes do not match metadata: ', path)

                }
                ##
                if (!identical(unique(results$device), metadata$settings$device)) invisible(NULL)
                ##
                ## Older saved studies predate explicit Stan backend metadata; preserve their reader compatibility only.
                if (is.null(metadata$settings$stan_backend)) metadata$settings$stan_backend <- 'cmdstanr'
                metadata$settings$stan_backend <-  fn_paper1_stan_backend(settings = metadata$settings)
                if (is.null(metadata$settings$package_stack)) metadata$settings$package_stack <-  'legacy_bayesmvp'
                expected_backends <-  fn_paper1_execution_columns(cases = results[, setdiff(names(results), 'execution_backend'), drop = FALSE],
                                                                   stan_backend = metadata$settings$stan_backend,
                                                                   native_backend = fn_paper1_native_backend(metadata$settings))$execution_backend
                if ('execution_backend' %in% names(results) && !identical(results$execution_backend, expected_backends)) {

                    invisible(NULL)

                }
                results$execution_backend <-  expected_backends
                ##
                results <-  fn_paper1_stan_partition_columns(cases = results)
                results$source_case_id <-  results$case_id
                ##
                results$source_output_dir <-  path
                ##
                ## ---- Every Mplus long-run mode that matches the current runner settings, kept for the within-Mplus mode comparison:
                ##
                mplus_mode_results <-  results[grepl('^Mplus_', results$algorithm) & !is.na(results$mplus_iteration_mode), , drop = FALSE]
                ##
                ## ---- Several Mplus long-run modes in one study: the manuscript tables use the first mode listed in
                ## ---- the runner; the other mode's rows stay in results.rds and in mplus_iteration_mode_comparison.csv.
                ##
                report_mplus_mode <-  metadata$settings$mplus_iteration_mode[1]
                if (!is.null(report_mplus_mode) && 'mplus_iteration_mode' %in% names(results)) {
                    other_mode_rows <-  grepl('^Mplus_', results$algorithm) & !is.na(results$mplus_iteration_mode) &
                                        results$mplus_iteration_mode != report_mplus_mode
                    if (any(other_mode_rows)) {
                        message(paste0('\033[36mReport uses Mplus ', report_mplus_mode, ' rows only; ', sum(other_mode_rows),
                                       ' rows in other modes are kept for the within-Mplus comparison (', path, ').\033[0m'))
                        results <-  results[!other_mode_rows, , drop = FALSE]
                    }
                }
                ##
                list(results = results, metadata = metadata, input_md5 = tools::md5sum(files), mplus_mode_results = mplus_mode_results)

        })
        ##
        studies <-  Filter(f = Negate(f = is.null), x = studies)
        if (!length(studies)) stop('No saved runs match the settings currently written in the runner.')
        ##
        ## Separate Mplus-only runs may share a device with the BayesMVP/Stan study.
        ## Compare sampler settings only between studies that actually measured that sampler.
        setting_algorithms <-  list( data_seed = ".*",
                                      bayesmvp_iterations = "^MD_",
                                      stan_iterations = "^AD_",
                                      bayesmvp = "^MD_",
                                      stan = "^AD_",
                                      stan_backend = "^AD_",
                                      stan_via_NicoStan = "^AD_",
                                      package_stack = "^MD_",
                                      stan_reduce_sum_type = "^AD_Stan_(WCP|tape_chunked)$",
                                      ## Timing method and short/long iteration settings; NULL in older studies.
                                      timing_method = ".*",
                                      bayesmvp_iterations_short_run = "^MD_",
                                      bayesmvp_wcp_iterations = "^MD_BayesMVP_WCP$",
                                      stan_wcp_iterations = "^AD_Stan_WCP$",
                                      stan_iterations_short_run = "^AD_",
                                      ## Keep direct comparisons for fields saved by historical Mplus studies.
                                      ## Missing values stay missing; no current setting is inferred from them.
                                      mplus_FBITERATIONS = "^Mplus_",
                                      mplus_FBITERATIONS_short_run = "^Mplus_",
                                      ## Current Mplus execution controls apply only when Mplus was measured.
                                      mplus_WCP_iterations = "^Mplus_WCP$",
                                      mplus_standard_iterations = "^Mplus_standard$",
                                      mplus_iterations_short_run = "^Mplus_",
                                      mplus_iteration_mode = "^Mplus_",
                                      mplus_iteration_mode_short_run = "^Mplus_",
                                      mplus_short_run_role = "^Mplus_",
                                      mplus_biterations_minimum = "^Mplus_",
                                      mplus_bconvergence = "^Mplus_",
                                      mplus_save_draws = "^Mplus_",
                                      mplus_verify_saved_draws = "^Mplus_")
        ##
        for (field in names(x = setting_algorithms)) {

            applicable_studies <-  Filter( f = function(study)
                any(grepl(pattern = setting_algorithms[[field]], x = study$results$algorithm)),
                x = studies)
            ##
            for (study in applicable_studies) {

                if (!identical(study$metadata$settings[[field]], applicable_studies[[1]]$metadata$settings[[field]])) {

                    invisible(NULL)

                }

            }

        }
        ##
        mplus_studies <-  Filter(f = function(study) any(grepl(pattern = "^Mplus_", x = study$results$algorithm)), x = studies)
        for (study in mplus_studies) {

            prior_fields <-  grep(pattern = "^prior_", x = names(x = study$metadata$settings$mplus$settings), value = TRUE)
            if (!identical(study$metadata$settings$mplus$settings[prior_fields],
                           mplus_studies[[1]]$metadata$settings$mplus$settings[prior_fields])) {

                invisible(NULL)

            }

        }
        ##
        for (study in studies) {

            for (field in c('generator', 'data_seed', 'corr_force_positive_DGM')) {

                if (!identical(study$metadata$data[[field]], studies[[1]]$metadata$data[[field]])) invisible(NULL)

            }
            ##
            ## ---- A different simulator source file (e.g. edited comments on one machine) is not by itself a different
            ## ---- dataset: it sends the studies to the exact comparison of the simulated data below, like differing dataset hashes.
            ##
            simulator_source_differs <-  !identical(study$metadata$data$simulator_md5, studies[[1]]$metadata$data$simulator_md5)
            if (simulator_source_differs) {
                message(paste0('\033[36mSimulator source file differs between ', basename(study$results$source_output_dir[1]), ' and ',
                               basename(studies[[1]]$results$source_output_dir[1]), '; comparing the simulated data themselves.\033[0m'))
            }
            ##
            shared_N <-  intersect(names(study$metadata$data$dataset_md5), names(studies[[1]]$metadata$data$dataset_md5))
            ##
            if (simulator_source_differs ||
                !identical( study$metadata$data$dataset_md5[shared_N],
                            studies[[1]]$metadata$data$dataset_md5[shared_N])) {

                ##
                ## ---- The saved dataset_md5 is the md5 of a serialised (saveRDS) copy of y, so it also depends on the R version that
                ##      wrote it: the HPC and laptop studies  had different hashes for bit-identical data. Compare the
                ##      simulated outcomes themselves before calling the studies different (the md5 stays as saved, it keys the cache).
                ##
                dataset_file_this <-  file.path(study$results$source_output_dir[1], 'COVID19_datasets.rds')
                dataset_file_first <-  file.path(studies[[1]]$results$source_output_dir[1], 'COVID19_datasets.rds')
                if (!file.exists(dataset_file_this) || !file.exists(dataset_file_first)) invisible(NULL)
                datasets_this <-  readRDS(dataset_file_this)
                datasets_first <-  readRDS(dataset_file_first)
                for (N_label in shared_N) {

                    index_this <-  match(as.numeric(N_label), datasets_this$N_vec)
                    index_first <-  match(as.numeric(N_label), datasets_first$N_vec)
                    if (is.na(index_this) || is.na(index_first) ||
                        !identical(datasets_this$y_list[[index_this]], datasets_first$y_list[[index_first]]) ||
                        !identical(datasets_this$X_list[[index_this]], datasets_first$X_list[[index_first]]) ||
                        !identical(datasets_this$pop_list[[index_this]], datasets_first$pop_list[[index_first]])) {

                        invisible(NULL)

                    }

                }
                message(paste0('\033[36mDataset hashes differ between ', basename(study$results$source_output_dir[1]), ' and ',
                               basename(studies[[1]]$results$source_output_dir[1]),
                               ' but the simulated data (y, X, pop) are identical for N = ', paste(shared_N, collapse = ', '),
                               '; the hash depends on the R version that wrote it.\033[0m'))

            }

        }
        ##
        results <-  do.call(rbind, lapply(studies, `[[`, 'results'))
        ##
        ## Corrected and single-run times differ by the fixed per-run cost, so never mix these methods.
        ##
        if (length(unique(results$timing_method)) != 1) {

            stop('This results set mixes timing methods (', paste(unique(results$timing_method), collapse = ', '),
                 '); report each timing method separately.')

        }
        ##
        estimator_groups <-  unique(results[, c('device', 'algorithm', 'execution_backend', 'timing_estimator'), drop = FALSE])
        if (anyDuplicated(estimator_groups[, c('device', 'algorithm', 'execution_backend'), drop = FALSE])) {
            invisible(NULL)
        }
        iteration_groups <-  unique(results[, c('device', 'algorithm', 'execution_backend', 'N', 'n_iter', 'timing_estimator'), drop = FALSE])
        if (anyDuplicated(iteration_groups[, c('device', 'algorithm', 'execution_backend', 'N'), drop = FALSE])) {

            invisible(NULL)

        }
        ##
        if (anyDuplicated(fn_paper1_case_identity(results))) invisible(NULL)
        ##
        results$case_id <-  seq_len(nrow(results))
        ##
        devices <-  vapply(X = studies, FUN = function(study) study$metadata$settings$device, FUN.VALUE = character(length = 1))
        device_studies <-  lapply( X = unique(x = devices),
                                   FUN = function(device) {

                source_studies <-  studies[devices == device]
                combined_study <-  source_studies[[1]]
                combined_study$results <-  do.call(what = rbind, args = lapply(X = source_studies, FUN = `[[`, 'results'))
                combined_study$metadata$source_study_metadata <-  lapply(X = source_studies, FUN = `[[`, 'metadata')
                combined_study$metadata$settings$n_threads_vec <-  sort(x = unique(x = unlist(x = lapply(
                    X = source_studies, FUN = function(study) study$metadata$settings$n_threads_vec))))
                combined_study$input_md5 <-  unlist(x = lapply(X = source_studies, FUN = `[[`, 'input_md5'))
                combined_study

        })
        ##
        list(results = results, studies = setNames(object = device_studies, nm = unique(x = devices)), paths = paths,
             mplus_mode_results = do.call(what = rbind, args = lapply(X = studies, FUN = `[[`, 'mplus_mode_results')))

}
##
## ---- Prepare saved results for presentation ------------------------------------------------------------------------------------------
##
fn_paper1_presentation_data <-  function( results,
                                          studies) {

        analysis <-  fn_summarise_paper1_benchmark(results = results)
        ##
        df <-  analysis$configurations
        ##
        df$N_num <-  df$N
        ##
        df$N_label <-  factor( paste0('N = ', fn_paper1_format_number_commas_from_10000(df$N)),
                               levels = paste0('N = ', fn_paper1_format_number_commas_from_10000(sort(unique(df$N)))))
        ##
        df$N_chunks_num <-  df$num_chunks
        ##
        df$N_chunks <-  factor(df$num_chunks, levels = sort(unique(df$num_chunks)))
        ##
        df$N_threads <-  factor(df$n_threads, levels = sort(unique(df$n_threads)))
        ##
        df$N_iter <-  df$n_iter
        ##
        df$N_chains <-  df$n_chains
        ##
        df$time_avg <-  df$time_mean
        ##
        df$time_SD <-  df$time_sd
        ##
        df$Efficiency <-  df$chain_rate
        ##
        labels <-  c(MD_BayesMVP = 'BayesMVP + chunking', MD_BayesMVP_WCP = 'BayesMVP + chunking + WCP',
                     MD_BayesMVP_multi_process = 'BayesMVP multiprocessing', AD_Stan = 'Stan',
                     AD_Stan_chunked = 'Stan + chunking', AD_Stan_tape_chunked = 'Stan + tape chunking',
                     AD_Stan_WCP = 'Stan + tape chunking + WCP', Mplus_standard = 'Mplus', Mplus_WCP = 'Mplus + WCP')
        ##
        df$Algorithm_label <-  unname(labels[df$algorithm])
        via_NicoStan <-  df$execution_backend == 'NicoStan_BridgeStan'
        df$Algorithm_label[via_NicoStan] <-  sub(pattern = '^Stan', replacement = 'Stan model (NicoStan)',
                                                x = df$Algorithm_label[via_NicoStan])
        ##
        df$Stan_variant_label <-  df$Algorithm_label
        ##
        main_threads <-  vapply( seq_len(nrow(df)),
                                 function(i) df$n_threads[i] %in% studies[[df$device[i]]]$metadata$settings$n_threads_vec,
                                 logical(1))
        ##
        if ('benchmark_role' %in% names(df)) main_threads <-  main_threads & df$benchmark_role == 'main_scaling'
        ##
        ps1 <-  df[df$algorithm == 'MD_BayesMVP' & main_threads, , drop = FALSE]
        ## Keep every chunk measurement and every conditional optimum for both WCP arms.
        wcp_chunk_search <-  df[match(x = analysis$wcp_chunk_search$configuration_id, table = df$configuration_id), , drop = FALSE]
        wcp_chunk_search$selected_best_chunks <-  analysis$wcp_chunk_search$selected_best_chunks
        wcp_chunk_search$relative_chunk_throughput <-  analysis$wcp_chunk_search$relative_chunk_throughput
        ##
        wcp_best_chunks <-  wcp_chunk_search[wcp_chunk_search$selected_best_chunks, , drop = FALSE]
        ##
        ## Filter chunks before comparing algorithms, retaining all chain/WCP allocations at this stage.
        serial <-  df[df$algorithm == 'MD_BayesMVP' & df$num_chunks == 1 & main_threads, , drop = FALSE]
        ##
        serial$Algorithm_label <-  rep('BayesMVP (1 chunk)', nrow(serial))
        ##
        WCP_rows <-  df$algorithm %in% c('MD_BayesMVP_WCP', 'AD_Stan_WCP', 'Mplus_WCP')
        if ('benchmark_role' %in% names(x = df)) WCP_rows <-  WCP_rows & df$benchmark_role == 'main_scaling'
        ##
        candidates <-  df[(df$algorithm != 'MD_BayesMVP' | df$num_chunks > 1) &
            (main_threads | WCP_rows), , drop = FALSE]
        ##
        chunked_WCP_rows <-  candidates$algorithm %in% c('MD_BayesMVP_WCP', 'AD_Stan_WCP')
        candidates <-  candidates[!chunked_WCP_rows | candidates$configuration_id %in% wcp_best_chunks$configuration_id, , drop = FALSE]
        ##
        candidates <-  candidates[order(-candidates$chain_rate, candidates$num_chunks, candidates$threads_per_chain), , drop = FALSE]
        ##
        allocation_columns <-  c('device', 'algorithm', 'execution_backend', 'N', 'n_threads', 'n_chains', 'threads_per_chain',
                                  'n_iter', 'dataset_md5', 'timing_scope', 'timing_estimator')
        candidates <-  candidates[!duplicated(x = candidates[, allocation_columns, drop = FALSE]), , drop = FALSE]
        ##
        comparison_candidates <-  rbind(serial, candidates)
        ##
        ## Only the resource-scaling view also selects between allocations using the same total threads.
        scaling_columns <-  setdiff(x = allocation_columns, y = c('n_chains', 'threads_per_chain'))
        best <-  candidates[!duplicated(x = candidates[, scaling_columns, drop = FALSE]), , drop = FALSE]
        ##
        scaling <-  rbind(serial, best)
        ##
        scaling <-  scaling[order(scaling$device, scaling$N, scaling$Algorithm_label, scaling$n_threads), , drop = FALSE]
        scaling <-  fn_paper1_adjusted_scaling(configurations = scaling)
        NicoStan_comparisons <-  fn_paper1_NicoStan_comparisons(configurations = df)
        BayesMVP_comparisons <-  fn_paper1_BayesMVP_comparisons(configurations = df)
        ##
        list( analysis = analysis,
              configurations = df,
              ps1 = ps1,
              wcp_chunk_search = wcp_chunk_search,
              wcp_best_chunks = wcp_best_chunks,
              comparison_candidates = comparison_candidates,
              scaling = scaling,
              stan = scaling[grepl('^AD_', scaling$algorithm), , drop = FALSE],
              NicoStan_candidates = NicoStan_comparisons$candidates,
              NicoStan_by_budget = NicoStan_comparisons$by_budget,
              NicoStan_best_by_N = NicoStan_comparisons$best_by_N,
              BayesMVP_candidates = BayesMVP_comparisons$candidates,
              BayesMVP_by_budget = BayesMVP_comparisons$by_budget,
              BayesMVP_best_by_N = BayesMVP_comparisons$best_by_N)

}
##
## ---- Within-BayesMVP comparison: no chunking / chunking only / WCP only / chunking + WCP ----------------------------
##
## Mirrors fn_paper1_NicoStan_comparisons for the native BayesMVP arms. The four modes are:
##   "BayesMVP"                  MD_BayesMVP with one chunk and one thread per chain (the serial baseline);
##   "BayesMVP-chunking"         MD_BayesMVP with more than one chunk, one thread per chain;
##   "BayesMVP-WCP"              MD_BayesMVP_WCP with one chunk per WCP thread (num_chunks = threads_per_chain), i.e. no extra chunking;
##   "BayesMVP-chunking_and_WCP" every MD_BayesMVP_WCP configuration (the joint chunk / WCP search).
## best_by_N keeps, for each device, N and mode, the configuration with the highest total iterations per second over ALL
## chains / chunks / WCP allocations; by_budget keeps the best configuration for each total thread count as well.
##
fn_paper1_BayesMVP_comparisons <-  function( configurations ) {

        candidates <-  configurations[configurations$algorithm %in% c("MD_BayesMVP", "MD_BayesMVP_WCP"), , drop = FALSE]
        if (!nrow(candidates)) return(list(candidates = candidates, by_budget = candidates, best_by_N = candidates))
        ## A physical WCP-only case may also win the combined search, exactly as for the Stan arms.
        selections <-  list(
            "BayesMVP" = candidates$algorithm == "MD_BayesMVP" & candidates$num_chunks == 1,
            "BayesMVP-chunking" = candidates$algorithm == "MD_BayesMVP" & candidates$num_chunks > 1,
            "BayesMVP-WCP" = candidates$algorithm == "MD_BayesMVP_WCP" & candidates$num_chunks == candidates$threads_per_chain,
            "BayesMVP-chunking_and_WCP" = candidates$algorithm == "MD_BayesMVP_WCP")
        rows <-  lapply(X = names(selections), FUN = function(mode) {

            selected <-  candidates[selections[[mode]], , drop = FALSE]
            selected$comparison_mode <-  rep(mode, nrow(selected))
            selected$Algorithm_label <-  selected$Stan_variant_label <-  rep(mode, nrow(selected))
            return(selected)

        })
        candidates <-  do.call(what = rbind, args = rows)
        candidates <-  candidates[order(-candidates$total_iter_per_sec, candidates$n_threads,
                                        candidates$num_chunks, candidates$threads_per_chain), , drop = FALSE]
        groups <-  c("device", "N", "comparison_mode", "n_iter", "dataset_md5", "timing_scope", "timing_estimator")
        best_by_N <-  candidates[!duplicated(candidates[, groups, drop = FALSE]), , drop = FALSE]
        by_budget <-  candidates[!duplicated(candidates[, c(groups, "n_threads"), drop = FALSE]), , drop = FALSE]
        ##
        return(list(candidates = candidates, by_budget = by_budget, best_by_N = best_by_N))

}
##
fn_paper1_NicoStan_comparisons <-  function( configurations ) {

        candidates <-  configurations[configurations$execution_backend == "NicoStan_BridgeStan", , drop = FALSE]
        if (!nrow(candidates)) return(list(candidates = candidates, by_budget = candidates, best_by_N = candidates))
        ## Preserve the measured algorithm and case IDs; a physical WCP-only case may also win the combined search.
        selections <-  list(
            "NicoStan" = candidates$algorithm == "AD_Stan",
            "NicoStan-chunking" = candidates$algorithm == "AD_Stan_tape_chunked",
            "NicoStan-WCP" = candidates$algorithm == "AD_Stan_WCP" & candidates$num_chunks == candidates$threads_per_chain,
            "NicoStan-chunking_and_WCP" = candidates$algorithm == "AD_Stan_WCP")
        rows <-  lapply(X = names(selections), FUN = function(mode) {

            selected <-  candidates[selections[[mode]], , drop = FALSE]
            selected$comparison_mode <-  rep(mode, nrow(selected))
            selected$Algorithm_label <-  selected$Stan_variant_label <-  rep(mode, nrow(selected))
            return(selected)

        })
        candidates <-  do.call(what = rbind, args = rows)
        candidates <-  candidates[order(-candidates$total_iter_per_sec, candidates$n_threads,
                                        candidates$num_chunks, candidates$threads_per_chain), , drop = FALSE]
        groups <-  c("device", "N", "comparison_mode", "n_iter", "dataset_md5", "timing_scope", "timing_estimator")
        best_by_N <-  candidates[!duplicated(candidates[, groups, drop = FALSE]), , drop = FALSE]
        by_budget <-  candidates[!duplicated(candidates[, c(groups, "n_threads"), drop = FALSE]), , drop = FALSE]
        ##
        return(list(candidates = candidates, by_budget = by_budget, best_by_N = best_by_N))

}
##
## ---- Stan variants figure: add the WCP-only arm (one chunk per WCP thread) ---------------------------------------------------------
##
## The WCP-only rows are exactly the "NicoStan-WCP" rows of fn_paper1_NicoStan_comparisons()$by_budget: AD_Stan_WCP with
## num_chunks = threads_per_chain, keeping the highest total iterations per second at each device, N and total thread budget
## (the same rule as "BayesMVP-WCP"). They are relabelled in the Stan figure style; the four existing arms (views$stan) are
## passed through unchanged. Columns present in only one of the two inputs are filled with NA.
##
fn_paper1_Stan_variants_with_WCP_only <-  function( stan,
                                                    NicoStan_by_budget,
                                                    WCP_only_label = "Stan model (NicoStan) + WCP") {

        WCP_only <-  NicoStan_by_budget[NicoStan_by_budget$comparison_mode == "NicoStan-WCP", , drop = FALSE]
        ##
        if (!nrow(WCP_only)) return(stan)
        ##
        WCP_only$Algorithm_label <-  WCP_only$Stan_variant_label <-  rep(WCP_only_label, nrow(WCP_only))
        ##
        for (column in setdiff(x = names(stan), y = names(WCP_only))) WCP_only[[column]] <-  rep(NA, nrow(WCP_only))
        ##
        for (column in setdiff(x = names(WCP_only), y = names(stan))) stan[[column]] <-  rep(NA, nrow(stan))
        ##
        return(rbind(stan, WCP_only[, names(stan), drop = FALSE]))

}
##
## ---- Export manuscript figures, tables and provenance --------------------------------------------------------------------------------
##
## ---- Serial baseline for every arm, matched PER ITERATION: the WCP arms run more long-run iterations than their
##      serial counterparts (e.g. 2,000 vs 200 at N = 500), so requiring identical n_iter - as this helper did Previously -
##      left every WCP speed-up and efficiency blank. Both times are start-up-corrected long-run times (= n_iter x t), so the
##      comparison is made on seconds per iteration; serial_seconds is the counterpart's serial run scaled to this row's n_iter.
##
fn_paper1_serial_efficiency <- function(selected, configurations) {
        out <- selected
        out$serial_seconds <- out$serial_speedup <- out$parallel_efficiency <- rep(NA_real_, nrow(out))
        out$serial_n_iter <- rep(NA_real_, nrow(out))
        out$serial_case_ids <- rep(list(integer()), nrow(out))
        counterparts <- c(MD_BayesMVP_WCP = "MD_BayesMVP", AD_Stan_WCP = "AD_Stan_tape_chunked", Mplus_WCP = "Mplus_standard")
        for (i in seq_len(nrow(out))) {
            algorithm <- as.character(out$algorithm[i])
            baseline_algorithm <- if (algorithm %in% names(counterparts)) counterparts[[algorithm]] else algorithm
            same <- configurations$algorithm == baseline_algorithm & configurations$n_chains == 1 & configurations$n_threads == 1
            for (field in c("device", "execution_backend", "N", "num_chunks", "dataset_md5", "timing_scope", "timing_estimator")) {
                same <- same & configurations[[field]] == out[[field]][i]
            }
            base <- which(same)
            if (length(base) == 1L && configurations$time_mean[base] > 0 && out$time_mean[i] > 0 &&
                configurations$n_iter[base] > 0 && out$n_iter[i] > 0) {
                serial_seconds_per_iteration <- configurations$time_mean[base] / configurations$n_iter[base]
                out$serial_n_iter[i] <- configurations$n_iter[base]
                out$serial_seconds[i] <- serial_seconds_per_iteration * out$n_iter[i]
                out$serial_speedup[i] <- out$n_chains[i] * out$serial_seconds[i] / out$time_mean[i]
                out$parallel_efficiency[i] <- out$serial_speedup[i] / out$n_threads[i]
                out$serial_case_ids[[i]] <- configurations$case_ids[[base]]
            }
        }
        out
}
##
##
## ---- fn_paper1_exclude_report_cases: drop measured grid cells the report should not show  ------------------------------
##
## excluded_cases is a list of rules; each rule is a named list over any of device, algorithm, N, n_chains, threads_per_chain,
## num_chunks, and a row is dropped when it matches EVERY field of a rule. Used because the HPC WCP grid ran threads_per_chain = 6
## with 8 chains at N = 500 but not at N = 2,500, so that cell is left out to keep the N grids comparable. Every dropped row is
## counted in the message and in the report READ_ME; the saved results.rds files are untouched.
##
fn_paper1_exclude_report_cases <-  function( results,
                                             excluded_cases) {

        if (!is.list(excluded_cases)) stop('excluded_cases must be a list of rules (use list() for none).')
        if (length(excluded_cases) == 0) return(list(results = results, notes = character(0)))
        ##
        allowed_fields <-  c('device', 'algorithm', 'N', 'n_chains', 'threads_per_chain', 'num_chunks')
        rows_to_drop <-  rep(FALSE, nrow(results))
        notes <-  character(0)
        ##
        for (rule in excluded_cases) {

            if (!is.list(rule) || is.null(names(rule)) || !all(names(rule) %in% allowed_fields)) {
                stop('Each excluded_cases rule must be a named list over: ', paste(allowed_fields, collapse = ', '))
            }
            rule_rows <-  rep(TRUE, nrow(results))
            for (field in names(rule)) rule_rows <-  rule_rows & results[[field]] %in% rule[[field]]
            rule_label <-  paste(paste0(names(rule), ' = ', vapply(rule, function(value) paste(value, collapse = '/'), character(1))),
                                 collapse = ', ')
            notes <-  c(notes, paste0('Excluded from this report (design consistency): ', rule_label, ' - ', sum(rule_rows),
                                      ' measured repeat(s) dropped.'))
            message(paste0('\033[36mfn_paper1_exclude_report_cases: ', rule_label, ': dropping ', sum(rule_rows), ' measured repeat(s).\033[0m'))
            rows_to_drop <-  rows_to_drop | rule_rows

        }
        list(results = results[!rows_to_drop, , drop = FALSE], notes = notes)

}
##
fn_paper1_export_manuscript <-  function( study_output_dirs,
                                          output_dir = NULL,
                                          helper_dir = path.expand('~/Documents/Work/PhD_work/Alg_paper_analysis/paper_1_chunking_and_parallel_scalability'),
                                          manuscript_dir = NULL,
                                          excluded_cases = list()) {

        fn_paper1_report_dependencies()
        ##
        private <-  new.env(parent = environment())
        ##
        sys.source(file.path(helper_dir, 'R_fns_alg_paper_1_chunking_WCP_par_scaling.R'), envir = private)
        ##
        sys.source(file.path(helper_dir, 'R_fns_alg_paper_1_presentation_templates.R'), envir = private)
        # Reader and adapter resolve the benchmark functions in a private environment, not caller globals.
        reader <-  fn_paper1_read_studies
        ##
        environment(reader) <-  private
        ##
        adapter <-  fn_paper1_presentation_data
        ##
        environment(adapter) <-  private
        ##
        input <-  reader(study_output_dirs, runner_file = file.path(helper_dir, 'alg_paper_1_chunking_WCP_par_scaling.R'))
        ##
        ## ---- Grid cells left out of the report (design consistency; see fn_paper1_exclude_report_cases):
        ##
        excluded_report_cases <-  fn_paper1_exclude_report_cases(results = input$results, excluded_cases = excluded_cases)
        input$results <-  excluded_report_cases$results
        exclusion_notes <-  excluded_report_cases$notes
        ##
        views <-  adapter(input$results, input$studies)
        ##
        if (is.null(output_dir)) output_dir <-  file.path(input$paths[1], 'manuscript_outputs')
        ##
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        ##
        output_dir <-  normalizePath(output_dir, mustWork = TRUE)
        ##
        figdir <-  file.path(output_dir, 'figures')
        ##
        tabledir <-  file.path(output_dir, 'tables')
        ##
        csvdir <-  file.path(output_dir, 'data')
        ##
        for (path in c(figdir, tabledir, csvdir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
        ##
        templates <-  private$fn_paper1_presentation_templates()
        ##
        ## ---- Timing method label for every note, caption and figure:
        ##
        timing_method <-  unique(input$results$timing_method)
        timing_label <-  fn_paper1_timing_label(timing_method = timing_method,
                                               timing_estimators = unique(input$results$timing_estimator))
        plot_timing_label <-  if ('overhead_subtraction' %in% input$results$timing_estimator)
            'Timing: Mplus overhead subtraction; other arms retain their recorded timing estimator.' else paste('Timing:', timing_method)
        fn_caption_with_timing <-  function(caption) paste(caption, timing_label)
        ##
        notes <-  c(exclusion_notes,
                    'Fresh results only. No legacy timing arrays were imported.',
                    'Times are means across saved repeats; error bars are +/- one sample SD.',
                    'Within-method selection uses chain_rate = n_chains / mean elapsed seconds, at fixed iteration budgets.',
                    'Serial baselines for the WCP arms (matched serial/WCP speed-ups, and the 1-thread and serial-baseline columns of the ps2 tables) use the counterpart arm (BayesMVP + chunking, Stan + tape chunking, Mplus) at 1 chain, scaled per iteration to the WCP arm iteration count, because the WCP arms run more long-run iterations and have no 1-chain run of their own. Where an HPC WCP arm has no 180-thread cell (its grid stops at 176), its largest measured budget stands in and the table says so.',
                    ## paste0('Cross-sampler scaling matches the manuscript figure: adjusted_scaling = (n_chains / time) * ',
                    ##        'that algorithm\'s own minimum time, separately for each device and N. This is not absolute sampling efficiency.'),
                    paste0('The per-N algorithm comparison figures and tables use adjusted_scaling = (n_chains / time) * ',
                           'that algorithm\'s own minimum time, separately for each device and N. This is not absolute sampling efficiency. ',
                           'The manuscript scaling figure (Figure_ps2_plot_2_adj_scalability.png) instead shows the speed-up over the one-chain, ',
                           'one-thread reference of the ps2 tables; its plotted values are in Figure_ps2_plot_2_adj_scalability_values.csv.'),
                    'Iterations/second are used only in the separate within-sampler implementation comparisons.',
                    paste0('Stage 1: wcp_chunk_search.csv retains every measured chunk setting in BayesMVP-WCP and Stan-WCP. ',
                           'wcp_best_chunks.csv selects chunks separately at each device, algorithm, N, n_chains and threads_per_chain, ',
                           'holding iterations, data and timing scope fixed. Exact throughput ties choose fewer chunks.'),
                    paste0('Stage 2: comparison_candidates.csv combines those conditional WCP optima with the other algorithms. ',
                           'All WCP allocations are retained here. Mplus-WCP has no separate chunk-selection stage.'),
                    'Allocation tables show a dash for Mplus chunks because this study does not expose an independent Mplus chunk setting.',
                    paste0('Scaling plots then choose the highest throughput allocation for each algorithm at each actual total thread budget; ',
                           'scaling_selected.csv records that choice, including chains and WCP threads. No unmeasured budget is filled in.'),
                    paste0('optimal_combinations.csv is an additional summary over WCP choices at fixed chain count, separately by algorithm arm. ',
                           'It does not replace the conditional chunk filter or merge WCP and chunking-only arms.'),
                    paste0('Matched serial/WCP speedups are supplementary and use identical chains, chunks and data; corrected times are compared per iteration when budgets differ. ',
                           'Missing serial baselines remain unavailable and do not prevent the report.'),
                    'Timing scopes differ by implementation and are retained in the CSV data. These are timing studies, not ESS comparisons.',
                    'Chunk choices are selected on these measurements; their observed advantage is not an independent validation estimate.',
                    'Cross-device package binaries may differ; each device build is preserved in report_manifest.rds.',
                    paste0('Timing method: ', timing_method, '. ', timing_label,
                           if (identical(timing_method, 'two_run_difference'))
                               paste0(' With T(n) = S + n t: t = (T_long - T_short) / (n_long - n_short), S = T_short - n_short t, ',
                                      'reported time = T_long - S = n_long t for iteration_difference rows. ',
                                      'For overhead_subtraction rows, S is the whole short control; seconds_per_iteration divides T_long - T_short by the observed long-run count when available, otherwise the requested count. ',
                                      'Reported time is that per-iteration value multiplied by the requested long-run count. ',
                                      'The control is an approximate overhead estimate, not a verified zero-work call. ',
                                      'Raw times, requested counts, timing_estimator, S and t are in measured_cases.csv.') else ''))
        ##
        ## ---- Two-run rows whose difference was not usable are listed, never dropped:
        ##
        two_run_flagged_cases <-  input$results[!is.na(input$results$two_run_timing_flag) & input$results$two_run_timing_flag != 'ok', , drop = FALSE]
        fn_paper1_write_csv(two_run_flagged_cases, file.path(csvdir, 'two_run_flagged_cases.csv'))
        ##
        if (nrow(two_run_flagged_cases)) {

            flagged_summary <-  paste0('ATTENTION: ', nrow(two_run_flagged_cases), ' of ', nrow(input$results),
                                       ' measured repeats have a two-run timing flag (long run not slower than the short run, or a fixed ',
                                       'cost below the negative tolerance). They are kept in every view; see two_run_flagged_cases.csv.')
            notes <-  c(notes, flagged_summary)
            message(paste0('\033[36m', flagged_summary, '\033[0m'))

        }
        ##
        if (any(views$configurations$time_mean_not_positive %in% TRUE)) {

            notes <-  c(notes, paste0('ATTENTION: ', sum(views$configurations$time_mean_not_positive %in% TRUE),
                                      ' configurations have a non-positive mean corrected time; their rates are NA and they are listed in configurations.csv.'))

        }
        ##
        if (any(input$results$mplus_PPPP_run %in% TRUE)) {

            mplus_rows <-  grepl('^Mplus_', input$results$algorithm)
            mplus_modes <-  if ('mplus_iteration_mode' %in% names(input$results))
                unique(as.character(input$results$mplus_iteration_mode[mplus_rows])) else character()
            mplus_modes <-  mplus_modes[!is.na(mplus_modes)]
            requested_mode_text <-  if (length(mplus_modes)) paste0(' under ', paste(mplus_modes, collapse = '/'), ' mode') else
                ' under the iteration mode recorded in study metadata, when available'
            ##
            notes <-  c(notes, paste0('Mplus times include a separate prior-posterior predictive (PPPP) run in addition to estimation. ',
                                      'Mplus n_iter records the requested per-chain count', requested_mode_text,
                                      '. Per-chain verification is optional; actual counts are reported only where saved-draw verification was enabled and results were saved.'))

        }
        if (any(input$results$timing_estimator == 'overhead_subtraction')) {
            notes <- c(notes, paste0('Mplus overhead controls need not save draws or complete one sampling iteration. ',
                'Their complete elapsed time is subtracted as an approximation. ',
                'Any iteration-dependent predictive checks and output work in the long call remain included. ',
                'Verification of the long call, when enabled, is still required; no actual short-run count is inferred.'))
            notes <- c(notes, paste0('Mplus long-run mode(s): ', paste(unique(input$results$mplus_iteration_mode[grepl('^Mplus_', input$results$algorithm)]), collapse = ', '),
                '; short overhead-control mode(s): ', paste(unique(input$results$mplus_iteration_mode_short_run[grepl('^Mplus_', input$results$algorithm)]), collapse = ', '), '.'))
        }
        ##
        if (any(views$configurations$execution_backend == 'NicoStan_BridgeStan')) {

            notes <-  c(notes,
                paste0('Stan model (NicoStan) denotes the Stan density/gradient via BridgeStan with the NicoStan HMC sampler, ',
                       'not CmdStan NUTS or its process-startup timings. The Stan step size is retained; the configured fixed trajectory length is recorded in study metadata.'),
                paste0('The threaded Stan backend uses a shared TBB budget of chains times WCP, not a separate per-chain thread cap. ',
                       'Nested reduce_sum tasks can share that budget even at nominal WCP = 1.'),
                paste0('Compilation and cached R initialisation are outside the case timer. State allocation, thread setup, ',
                       'C++ per-chain model loading, sampling and divergence checking are inside it. ',
                       'Initial states use an R uniform draw within the saved Stan init radius, conditioned on the existing class order, ',
                       'and seeded by the case seed. Divergence counts are retained using the same policy as the CmdStanR arm.'))

        }
        ##
        absent <-  setdiff(c('HPC', 'Laptop'), names(input$studies))
        ##
        if (length(absent)) notes <-  c(notes, paste('No supplied results for:', paste(absent, collapse = ', ')))
        ##
        if (any(views$configurations$max_divergences > 0, na.rm = TRUE))
            notes <-  c( notes,
                         'ATTENTION: saved diagnostics include divergences. Inspect configurations.csv before interpreting timings.')
        ##
        if (anyNA(input$results$divergences))
            notes <-  c( notes,
                         'Some measured repeats have unavailable divergence diagnostics; inspect measured_cases.csv.')
        ##
        for (name in names(views$analysis)) fn_paper1_write_csv(views$analysis[[name]], file.path(csvdir, paste0(name, '.csv')))
        ##
        fn_paper1_write_csv(input$results, file.path(csvdir, 'measured_cases.csv'))
        ##
        ## ---- Cross-software comparisons use tape chunking as the Stan chunking approach; container-chunked Stan
        ##      (AD_Stan_chunked) is reported only in the within-Stan implementation comparison.
        views$scaling <-  views$scaling[views$scaling$algorithm != 'AD_Stan_chunked', , drop = FALSE]
        ##
        fn_paper1_write_csv(views$scaling, file.path(csvdir, 'scaling_selected.csv'))
        views$serial_efficiency <- fn_paper1_serial_efficiency(views$scaling, views$configurations)
        fn_paper1_write_csv(views$serial_efficiency, file.path(csvdir, 'serial_efficiency.csv'))
        ##
        fn_paper1_write_csv(x = views$comparison_candidates, file = file.path(csvdir, 'comparison_candidates.csv'))
        for (view_name in c('NicoStan_candidates', 'NicoStan_by_budget', 'NicoStan_best_by_N',
                            'BayesMVP_candidates', 'BayesMVP_by_budget', 'BayesMVP_best_by_N')) {

            fn_paper1_write_csv(x = views[[view_name]], file = file.path(csvdir, paste0(view_name, '.csv')))

        }
        ##
        ## Preserve the ready-to-plot data frames and list-column provenance as well as flattened CSV exports.
        saveRDS(object = views, file = file.path(csvdir, 'presentation_views.rds'))
        ##
        coverage <-  stats::aggregate(list(timed_calls = rep(1L, nrow(input$results))), input$results[c('device', 'algorithm', 'N', 'timing_estimator')], sum)
        ##
        fn_paper1_write_csv(coverage, file.path(csvdir, 'coverage.csv'))
        ##
        assets <-  data.frame(kind = character(), file = character(), caption = character(), label = character())
        ##
        add_asset <-  function( kind,
                                file,
                                caption,
                                label) {

                if (nzchar(caption)) caption <-  fn_caption_with_timing(caption)
                if (file.exists(file.path(output_dir, file))) assets[nrow(assets) + 1L, ] <<-  list(kind, file, caption, label)

        }
        # Template functions print plots/tables; a report-local log and null graphics device avoid Rplots.pdf in the working directory.
        grDevices::pdf(NULL)
        ##
        device_id <-  grDevices::dev.cur()
        ##
        on.exit(if (device_id %in% grDevices::dev.list()) grDevices::dev.off(device_id), add = TRUE)
        ##
        log <-  file(file.path(output_dir, 'generation.log'), 'wt')
        ##
        sink(log)
        ##
        on.exit( {

                     sink()
                     ##
                     close(log)

                 },
                 add = TRUE)
        ##
        ## ---- First present the chunk search within each WCP arm ----------------------------------------------------------------------
        ##
        WCP_groups <-  unique(x = views$wcp_chunk_search[, c('device', 'algorithm', 'N'), drop = FALSE])
        WCP_groups <-  WCP_groups[order(WCP_groups$device, WCP_groups$algorithm, WCP_groups$N), , drop = FALSE]
        ##
        for (group_index in seq_len(length.out = nrow(x = WCP_groups))) {

            group <-  WCP_groups[group_index, , drop = FALSE]
            chunk_search <-  views$wcp_chunk_search[views$wcp_chunk_search$device == group$device &
                views$wcp_chunk_search$algorithm == group$algorithm & views$wcp_chunk_search$N == group$N, , drop = FALSE]
            chunk_search <-  chunk_search[order(chunk_search$n_chains, chunk_search$threads_per_chain, chunk_search$num_chunks), , drop = FALSE]
            best_chunks <-  chunk_search[chunk_search$selected_best_chunks, , drop = FALSE]
            ##
            file_prefix <-  paste0(group$device, '_paper_1_N_', group$N, '_algorithm_', group$algorithm)
            ##
            fn_paper1_write_csv(x = chunk_search, file = file.path(csvdir, paste0(file_prefix, '_chunk_search.csv')))
            fn_paper1_write_csv(x = best_chunks, file = file.path(csvdir, paste0(file_prefix, '_best_chunks_by_WCP.csv')))
            ##
            templates$fn_plot_paper1_WCP_chunk_search( chunk_search = chunk_search,
                                                       best_chunks = best_chunks,
                                                       output_path = figdir,
                                                       file_prefix = file_prefix)
            ##
            add_asset( kind = 'figure',
                       file = file.path('figures', paste0(file_prefix, '_chunk_search.png')),
                       caption = paste0(group$device, ', ', chunk_search$Algorithm_label[1], ', N = ', group$N,
                           ': all measured chunk settings, separately by chain count and WCP threads. ',
                           'Marked points maximise within-method chain throughput within each fixed chain/WCP allocation.'),
                       label = paste0('figure:', file_prefix, '_chunk_search'))
            ##
            add_asset( kind = 'figure',
                       file = file.path('figures', paste0(file_prefix, '_optimal_chunks_by_WCP.png')),
                       caption = paste0(group$device, ', ', chunk_search$Algorithm_label[1], ', N = ', group$N,
                           ': WCP throughput after optimising chunks separately at each chain and WCP count. ',
                           'Point labels give the selected chunk count; all measured WCP choices remain visible.'),
                       label = paste0('figure:', file_prefix, '_optimal_chunks_by_WCP'))
            ##
            ## ---- The same NicoStan+BayesMVP chunk search with cache-capacity lines, for N = 500 and 50,000 (the figure without lines is
            ##      kept): chunk count at which one chunk first fits in the L3 / L2 cache per active thread (chains x WCP threads).
            if (group$algorithm == 'MD_BayesMVP_WCP' && group$N %in% c(500, 50000)) {

                templates$fn_plot_paper1_WCP_chunk_search( chunk_search = chunk_search,
                                                           best_chunks = best_chunks,
                                                           output_path = figdir,
                                                           file_prefix = file_prefix,
                                                           show_cache_lines = TRUE,
                                                           thresholds_file = file.path(csvdir, paste0(file_prefix, '_chunk_search_cache_lines_thresholds.csv')))
                ##
                add_asset( kind = 'figure',
                           file = file.path('figures', paste0(file_prefix, '_chunk_search_cache_lines.png')),
                           caption = paste0(group$device, ', ', chunk_search$Algorithm_label[1], ', N = ', group$N,
                               ': all measured chunk settings (log scale), separately by chain count and WCP threads. ',
                               'Grey lines mark the chunk count at which one chunk (',
                               fn_paper1_format_number_commas_from_10000(templates$paper1_bytes_per_row_BayesMVP),
                               ' bytes per individual) first fits in the L3 or L2 cache per active thread (chains x WCP threads); ',
                               'to the left of an L3 line, one chunk exceeds the L3 cache per active thread.'),
                           label = paste0('figure:', file_prefix, '_chunk_search_cache_lines'))

            }
            ##
            ## ---- The same for the Stan model (NicoStan) + tape chunking + WCP chunk search, for N = 500 and 50,000 (the figure without
            ##      lines is kept), with the working set per individual of one Stan partial sum (autodiff tape, stacks and temporaries).
            if (group$algorithm == 'AD_Stan_WCP' && group$N %in% c(500, 50000)) {

                templates$fn_plot_paper1_WCP_chunk_search( chunk_search = chunk_search,
                                                           best_chunks = best_chunks,
                                                           output_path = figdir,
                                                           file_prefix = file_prefix,
                                                           show_cache_lines = TRUE,
                                                           bytes_per_row = templates$paper1_bytes_per_row_Stan,
                                                           thresholds_file = file.path(csvdir, paste0(file_prefix, '_chunk_search_cache_lines_thresholds.csv')))
                ##
                add_asset( kind = 'figure',
                           file = file.path('figures', paste0(file_prefix, '_chunk_search_cache_lines.png')),
                           caption = paste0(group$device, ', ', chunk_search$Algorithm_label[1], ', N = ', group$N,
                               ': all measured chunk settings (log scale), separately by chain count and WCP threads. ',
                               'Grey lines mark the chunk count at which one chunk (',
                               fn_paper1_format_number_commas_from_10000(templates$paper1_bytes_per_row_Stan),
                               ' bytes per individual: autodiff tape, its stacks and temporaries) first fits in the L3 or L2 cache per active thread ',
                               '(chains x WCP threads); to the left of an L3 line, one chunk exceeds the L3 cache per active thread.'),
                           label = paste0('figure:', file_prefix, '_chunk_search_cache_lines'))

            }
            ##
            table_filename <-  paste0(file_prefix, '_best_chunks_by_WCP.tex')
            table_label <-  paste0('table:', file_prefix, '_best_chunks_by_WCP')
            templates$fn_make_paper1_configuration_table( configurations = best_chunks,
                                                           caption = fn_caption_with_timing(paste0(group$device, ', ', chunk_search$Algorithm_label[1],
                                                               ', N = ', group$N, ': optimal chunks at each fixed chain and WCP count. ',
                                                               'Efficiency is chains divided by elapsed time within this implementation; ties choose fewer chunks.')),
                                                           label = table_label,
                                                           output_file = file.path(tabledir, table_filename))
            ##
            add_asset(kind = 'table', file = file.path('tables', table_filename), caption = '', label = table_label)

        }
        ##
        if (nrow(views$ps1)) {

            templates$R_fn_plot_ps1_N_chunks_ggplot_1(views$ps1, output_path = figdir)
            ##
            templates$R_fn_plot_ps1_N_chunks_ggplot_1(views$ps1, output_path = figdir, n_threads_for_HPC = 180, n_threads_for_Laptop = 16)
            ##
            ## Both thread settings on one plot (colour = SMT use; line type = device):
            templates$R_fn_plot_ps1_N_chunks_ggplot_SMT_combined(views$ps1, output_path = figdir)
            ##
            ## The same, with cache-capacity lines (chunk count at which one chunk first fits in the L3 / L2 / L1 cache per active thread):
            templates$R_fn_plot_ps1_N_chunks_ggplot_SMT_combined( views$ps1,
                                                                  output_path = figdir,
                                                                  show_cache_lines = TRUE,
                                                                  thresholds_file = file.path(csvdir, 'Figure_N_chunks_pilot_study_plot_1_n_threads_SMT_vs_no_SMT_cache_lines_thresholds.csv'))
            ##
            templates$R_fn_plot_ps1_efficiency(subset(views$ps1, device == 'HPC'), subset(views$ps1, device == 'Laptop'), output_path = figdir)
            ##
            ## Efficiency for both devices on one set of panels (colour = thread level; line type = device), without and with cache-capacity lines:
            templates$R_fn_plot_ps1_efficiency_combined(views$ps1, output_path = figdir)
            ##
            templates$R_fn_plot_ps1_efficiency_combined( views$ps1,
                                                         output_path = figdir,
                                                         show_cache_lines = TRUE,
                                                         thresholds_file = file.path(csvdir, 'Figure_N_chunks_pilot_study_plot_3_both_devices_cache_lines_thresholds.csv'))
            ##
            for (pair in list(c(96, 8), c(180, 16))) {

                filename <-  paste0('Figure_N_chunks_pilot_study_plot_1_n_threads_HPC_', pair[1], '_Laptop_', pair[2], '.png')
                ##
                add_asset( 'figure',
                           file.path('figures', filename),
                           paste0( 'COVID-19-derived data: BayesMVP sampling time against chunk count; HPC ',
                                   pair[1],
                                   ', Laptop ',
                                   pair[2],
                                   ' total threads. Mean and one SD across saved repeats; absent data are labelled.'),
                           paste0( 'figure:ps1_n_chunks_time_',
                                   paste( pair,
                                          collapse = '_')))

            }
            ##
            add_asset( 'figure',
                       'figures/Figure_N_chunks_pilot_study_plot_1_n_threads_SMT_vs_no_SMT.png',
                       paste0( 'COVID-19-derived data: BayesMVP sampling time against chunk count, without SMT (HPC 96, Laptop 8 total threads) ',
                               'and with SMT (HPC 180, Laptop 16 total threads). Mean and one SD across saved repeats; absent data are labelled.'),
                       'figure:ps1_n_chunks_time_SMT_combined')
            ##
            add_asset( 'figure',
                       'figures/Figure_N_chunks_pilot_study_plot_1_n_threads_SMT_vs_no_SMT_cache_lines.png',
                       paste0( 'COVID-19-derived data: BayesMVP sampling time against chunk count (log scale), without SMT (HPC 96, Laptop 8 total threads) ',
                               'and with SMT (HPC 180, Laptop 16 total threads). Vertical lines mark the chunk count at which one chunk (',
                               fn_paper1_format_number_commas_from_10000(templates$paper1_bytes_per_row_BayesMVP),
                               ' bytes per individual) first fits in the L3, L2 or L1 cache per active thread, in the colour and line type of its data line; ',
                               'to the left of the L3 line, one chunk exceeds the L3 cache per active thread. Mean and one SD across saved repeats.'),
                       'figure:ps1_n_chunks_time_SMT_combined_cache_lines')
            ##
            add_asset( 'figure',
                       'figures/Figure_N_chunks_pilot_study_plot_3.png',
                       'BayesMVP efficiency, number of threads divided by elapsed seconds, by chunk count, matching the manuscript PS1 definition.',
                       'figure:ps1_n_chunks_combined')
            ##
            add_asset( 'figure',
                       'figures/Figure_N_chunks_pilot_study_plot_3_both_devices.png',
                       paste0( 'BayesMVP efficiency, number of threads divided by elapsed seconds (log scale), by chunk count, for the HPC (solid lines; ',
                               '64, 96 and 180 threads) and the laptop (dashed lines; 4, 8 and 16 threads).'),
                       'figure:ps1_n_chunks_combined_both_devices')
            ##
            add_asset( 'figure',
                       'figures/Figure_N_chunks_pilot_study_plot_3_both_devices_cache_lines.png',
                       paste0( 'BayesMVP efficiency, number of threads divided by elapsed seconds (log scale), by chunk count (log scale), for the HPC (solid lines; ',
                               '64, 96 and 180 threads) and the laptop (dashed lines; 4, 8 and 16 threads). Vertical lines mark the chunk count at which one chunk (',
                               fn_paper1_format_number_commas_from_10000(templates$paper1_bytes_per_row_BayesMVP),
                               ' bytes per individual) first fits in the L3, L2 or L1 cache per active thread; thread levels with identical thresholds share one line.'),
                       'figure:ps1_n_chunks_combined_both_devices_cache_lines')
            ##
            best_chunks <-  templates$get_best_chunks(views$ps1)
            ##
            for (dev in names(input$studies)) {

                filename <-  paste0('table_ps1_best_n_chunks_', dev, '.tex')
                ##
                templates$make_ps1_best_chunks_table_tex( subset(best_chunks, device == dev),
                                                          dev,
                                                          input$studies[[dev]]$metadata$settings$n_threads_vec,
                                                          fn_caption_with_timing(paste( 'Best measured chunk counts on',
                                                                                        dev,
                                                                                        'at each total thread budget; ties select the smaller count.')),
                                                          paste0( 'table:ps1_best_n_chunks_',
                                                                  dev),
                                                          file.path( tabledir,
                                                                     filename))
                ##
                add_asset('table', file.path('tables', filename), '', paste0('table:ps1_best_n_chunks_', dev))

            }
            ##
            smt <-  lapply( names(input$studies),
                            function(dev) templates$get_smt_benefit( subset(views$ps1, device == dev),
                                                                     if (dev == 'HPC') 96 else 8,
                                                                     if (dev == 'HPC') 180 else 16))
            ##
            names(smt) <-  names(input$studies)
            ##
            for (type in c('threshold', 'best')) {

                selector <-  if (type == 'threshold') templates$get_smt_threshold else templates$get_smt_best
                ##
                selected <-  dplyr::bind_rows(lapply(names(smt), function(dev) selector(smt[[dev]], dev)))
                ##
                label <-  if (type == 'threshold') 'table:ps1_smt_threshold' else 'table:ps1_smt_best_n_chunks'
                ##
                filename <-  paste0(gsub(':', '_', label, fixed = TRUE), '.tex')
                ##
                templates$make_ps1_smt_table_tex( selected,
                                                  'Chunks',
                                                  fn_caption_with_timing(paste( 'SMT comparison:',
                                                                                if (type == 'threshold')
                                                                                    'smallest chunk count with a measured throughput gain.' else
                                                                                    'chunk count maximising measured SMT throughput.',
                                                                                'HPC 96 versus 180 threads; Laptop 8 versus 16. Missing comparisons are marked.')),
                                                  label,
                                                  file.path( tabledir,
                                                             filename))
                ##
                add_asset('table', file.path('tables', filename), '', label)

            }

        }
        ##
        if (nrow(views$scaling)) {

            ## Compare the manuscript's own-minimum-time scaling after the conditional chunk filter.
            comparison_groups <-  unique(x = views$scaling[, c('device', 'N'), drop = FALSE])
            comparison_groups <-  comparison_groups[order(comparison_groups$device, comparison_groups$N), , drop = FALSE]
            ##
            for (group_index in seq_len(length.out = nrow(x = comparison_groups))) {

                group <-  comparison_groups[group_index, , drop = FALSE]
                comparison <-  views$scaling[views$scaling$device == group$device & views$scaling$N == group$N, , drop = FALSE]
                comparison <-  comparison[order(comparison$Algorithm_label, comparison$n_threads), , drop = FALSE]
                ##
                file_prefix <-  paste0(group$device, '_paper_1_N_', group$N, '_algorithm_comparison')
                fn_paper1_write_csv(x = comparison, file = file.path(csvdir, paste0(file_prefix, '_selected.csv')))
                ##
                algorithm_plot <-  ggplot2::ggplot( data = comparison,
                                                    mapping = ggplot2::aes(x = n_threads, y = adjusted_scaling,
                                                                           colour = Algorithm_label, group = Algorithm_label)) +
                    ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 2.5) +
                    ggplot2::theme_bw(base_size = 16) +
                ggplot2::labs( title = paste0(group$device, ': N = ', fn_paper1_format_number_commas_from_10000(group$N)),
                               subtitle = 'Chunks optimised at fixed chains and WCP; best allocation at each measured thread budget',
                               x = 'Total threads', y = 'Adjusted scaling: (chains / time) x own minimum time', colour = 'Algorithm',
                               caption = plot_timing_label) +
                    ggplot2::theme(legend.position = 'bottom') +
                    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2, title.position = 'top'))
                ##
                figure_filename <-  paste0(file_prefix, '_scaling.png')
                ggplot2::ggsave(filename = file.path(figdir, figure_filename), plot = algorithm_plot, width = 14, height = 9, dpi = 150)
                ##
                add_asset( kind = 'figure', file = file.path('figures', figure_filename),
                           caption = paste0(group$device, ', N = ', group$N, ': parallel scaling normalised by each algorithm\'s own minimum time. ',
                               'Each algorithm uses its best remaining allocation at the measured total thread budget. ',
                               'Lines join measured points; allocations and timing scopes are recorded in the companion CSV/table.'),
                           label = paste0('figure:', file_prefix, '_scaling'))
                ##
                table_filename <-  paste0(file_prefix, '_selected.tex')
                table_label <-  paste0('table:', file_prefix, '_selected')
                templates$fn_make_paper1_configuration_table( configurations = comparison,
                                                               caption = fn_caption_with_timing(paste0(group$device, ', N = ', group$N,
                                                                   ': selected algorithm comparisons after conditional chunk optimisation. ',
                                                                   'All measured total thread budgets are retained. Scaling is (chains/time) times each algorithm\'s own minimum time.')),
                                                               label = table_label,
                                                               output_file = file.path(tabledir, table_filename),
                                                               include_algorithm = TRUE)
                ##
                add_asset(kind = 'table', file = file.path('tables', table_filename), caption = '', label = table_label)

            }
            ##
            ## templates$R_fn_plot_ps2_scaling(views$scaling, metric = 'adjusted', output_path = figdir)
            ##
            ## ---- The manuscript scaling figure uses the speed-up of the ps2 tables: S = N_chains x T0_eq / T against the one-chain,
            ##      one-thread reference (equation eq:paper1_serial_efficiency), not each arm's own minimum time.
            templates$R_fn_plot_ps2_scaling(views$scaling, metric = 'serial_speedup', output_path = figdir,
                                            values_file = file.path(csvdir, 'Figure_ps2_plot_2_adj_scalability_values.csv'))
            ##
            ## The same figure with light markers where SMT starts and the L3 cache per active thread falls below the L3 cache per core:
            templates$R_fn_plot_ps2_scaling(views$scaling, metric = 'serial_speedup', output_path = figdir,
                                            show_markers = TRUE,
                                            markers_file = file.path(csvdir, 'Figure_ps2_plot_2_adj_scalability_markers.csv'))
            serial_plot_data <- views$serial_efficiency[is.finite(views$serial_efficiency$parallel_efficiency), , drop = FALSE]
            if (nrow(serial_plot_data)) {
                p_serial <- ggplot2::ggplot(serial_plot_data,
                    ggplot2::aes(n_threads, 100 * parallel_efficiency, colour = Algorithm_label, group = Algorithm_label)) +
                    ggplot2::geom_hline(yintercept = 100, linetype = "dashed", colour = "grey50") +
                    ggplot2::geom_point(size = 3) + ggplot2::geom_line() + ggplot2::theme_bw(base_size = 18) +
                    ggplot2::facet_wrap(~device + N_label, scales = "free_x") +
                    ggplot2::scale_x_continuous(trans = "log2", breaks = c(1, 2, 4, 8, 16, 32, 64, 96, 180)) +
                    ggplot2::labs(x = "Total threads", y = "Parallel efficiency versus one thread (%)", colour = NULL,
                                  caption = plot_timing_label) +
                    ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
                    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2))
                ggplot2::ggsave(file.path(figdir, "Figure_paper1_serial_parallel_efficiency.png"), p_serial,
                               width = 16, height = 12, dpi = 100)
                add_asset('figure', 'figures/Figure_paper1_serial_parallel_efficiency.png',
                    'Parallel efficiency relative to one chain on one thread, matched at the same chunk count, data and timing estimator. Corrected times are compared per iteration when budgets differ. WCP uses its corresponding single-thread implementation as the serial reference. Missing serial references are not invented; normalised baseline times, measured iteration counts and identifiers are in the companion CSV.',
                    'figure:paper1_serial_parallel_efficiency')
            }
            ##
            add_asset( 'figure',
                       'figures/Figure_ps2_plot_2_adj_scalability.png',
                       ## paste0("The manuscript's adjusted scaling: (n_chains / time) times each implementation's own minimum time. ",
                       ##        'Chunks are first optimised at each fixed chain and WCP count. ',
                       ##        'The best remaining allocation is then selected separately by algorithm at each actual thread total. ',
                       ##        'Dashed lines show ideal scaling from each algorithm baseline.'),
                       paste0('Speed-up S = n_chains x T0_eq / time over the one-chain, one-thread reference of the ps2 tables. ',
                              'WCP arms use the counterpart arm (BayesMVP + chunking, Stan + tape chunking, Mplus) at one chain, scaled per iteration to the WCP arm iteration count. ',
                              'Chunks are first optimised at each fixed chain and WCP count. ',
                              'The best remaining allocation is then selected separately by algorithm at each actual thread total. ',
                              'The dashed grey line shows ideal scaling, S = total threads.'),
                       'figure:ps2_parallel_scalability_plot_1_adj_scalability')
            ##
            add_asset( 'figure',
                       'figures/Figure_ps2_plot_2_adj_scalability_markers.png',
                       paste0('Speed-up S = n_chains x T0_eq / time over the one-chain, one-thread reference of the ps2 tables, as in the figure without markers. ',
                              'The dotted vertical line marks the thread count after which SMT is in use and the L3 cache per active thread falls below ',
                              'the L3 cache per core (4 MB on the HPC, 2 MB on the laptop). The dashed grey line shows ideal scaling, S = total threads.'),
                       'figure:ps2_parallel_scalability_plot_1_adj_scalability_markers')
            ##
            templates$manus_rows <-  unique(views$scaling$Algorithm_label)
            ##
            templates$manus_to_df_label <-  setNames(as.list(templates$manus_rows), templates$manus_rows)
            ##
            for (dev in names(input$studies)) {

                threads <-  input$studies[[dev]]$metadata$settings$n_threads_vec
                ##
                low <-  min(threads)
                ##
                physical <-  if (dev == 'HPC') 96 else 8
                ##
                high <-  if (dev == 'HPC') 180 else 16
                ##
                Ns <-  sort(unique(views$scaling$N_num[views$scaling$device == dev]))
                ##
                for (hi in c(physical, high)) {

                    label <-  paste0('table:ps2_times_and_efficiency_', dev, '_', low, '_and_', hi, '_threads')
                    ##
                    filename <-  paste0(gsub(':', '_', label, fixed = TRUE), '.tex')
                    ##
                    templates$make_ratio_table_tex( views$scaling,
                                                    dev,
                                                    low,
                                                    hi,
                                                    Ns,
                                                    fn_caption_with_timing(paste0( dev,
                                                                                   ': elapsed seconds and high/low elapsed-time ratio at ',
                                                                                   low,
                                                                                   ' and ',
                                                                                   hi,
                                                                                   ' total threads. Best observed throughput configurations are selected separately at each budget. Chain counts can differ, so equal times need not mean equal throughput. Missing cells are marked. WCP arms have no 1-chain run: their 1-thread value is the counterpart arm (BayesMVP + chunking, Stan + tape chunking, Mplus) at 1 chain, scaled to the WCP arm iteration count, and their ratio is per thread (time x threads / chains), so that 1 means perfect throughput scaling for every row.')),
                                                    label,
                                                    file.path( tabledir,
                                                               filename))
                    ##
                    add_asset('table', file.path('tables', filename), '', label)

                }
                ##
                label <-  paste0('table:ps2_SMT_times_and_efficiency_', dev, '_', physical, '_and_', high, '_threads')
                ##
                filename <-  paste0(gsub(':', '_', label, fixed = TRUE), '.tex')
                ##
                templates$make_smt_gain_table_tex( views$scaling,
                                                   dev,
                                                   physical,
                                                   high,
                                                   Ns,
                                                   fn_caption_with_timing(paste0( dev,
                                                                                  ': elapsed seconds and percentage gain in chain throughput from ',
                                                                                  physical,
                                                                                  ' to ',
                                                                                  high,
                                                                                  ' threads. Selected configurations and chain counts are recorded in scaling\\_selected.csv. Where a WCP arm has no cell at the higher budget (the HPC WCP grids stop at 176 threads), its largest measured budget stands in (marked with an asterisk).')),
                                                   label,
                                                   file.path( tabledir,
                                                              filename))
                ##
                add_asset('table', file.path('tables', filename), '', label)
                ##
                ## ---- SMT gain against the 1-chain serial baseline: only when the study measured 1 chain.
                ##
                if (low == 1) {

                    label <-  paste0('table:ps2_SMT_serial_baseline_', dev, '_1_', physical, '_and_', high, '_threads')
                    ##
                    filename <-  paste0(gsub(':', '_', label, fixed = TRUE), '.tex')
                    ##
                    templates$make_smt_serial_baseline_table_tex( views$scaling,
                                                                  dev,
                                                                  1,
                                                                  physical,
                                                                  high,
                                                                  Ns,
                                                                  fn_caption_with_timing(paste0( dev,
                                                                                                 ': speed-up in chain throughput over the 1-chain serial run at ',
                                                                                                 physical, ' ($S_{', physical, '}$) and ', high, ' ($S_{', high, '}$) threads, ',
                                                                                                 'their difference (extra serial-chain equivalents from SMT), that difference per added logical thread (', high - physical, ' threads), ',
                                                                                                 'and the SMT gain $100 (S_{', high, '} / S_{', physical, '} - 1)$. ',
                                                                                                 'Best observed throughput configurations are selected separately at each budget, so chunks and chain/WCP splits can differ between budgets. Missing cells are marked. WCP arms use the counterpart arm at 1 chain (scaled per iteration) as their serial baseline, and their largest measured budget where the higher budget was not run (asterisk).')),
                                                                  label,
                                                                  file.path( tabledir,
                                                                             filename))
                    ##
                    add_asset('table', file.path('tables', filename), '', label)

                }

            }

        }
        ##
        if (nrow(views$stan)) {

            ## The total-throughput figure also shows the WCP-only arm (the NicoStan-WCP selection; see fn_paper1_Stan_variants_with_WCP_only).
            stan_variants_throughput <-  fn_paper1_Stan_variants_with_WCP_only(stan = views$stan, NicoStan_by_budget = views$NicoStan_by_budget)
            ##
            templates$R_fn_plot_Stan_variants_throughput(stan_variants_throughput, output_path = figdir, legend_nrow = 5)
            ##
            templates$R_fn_plot_Stan_variants_relative(views$stan, output_path = figdir)
            ##
            add_asset( 'figure',
                       'figures/Figure_ps2_Stan_variants_total_throughput.png',
                       'Within-sampler Stan implementation comparison: total iterations per second with the best observed chunk setting at each total thread budget.',
                       'figure:stan_throughput')
            ##
            add_asset( 'figure',
                       'figures/Figure_ps2_Stan_variants_relative_to_baseline.png',
                       'Stan variant chain throughput divided by plain Stan throughput at the same N and total thread budget.',
                       'figure:stan_relative')
            ## Keep the manuscript's existing grouped-N table function and extend it with selected configuration cells.
            for (dev in unique(views$stan$device)) {

                filename <-  paste0('table_Stan_implementations_best_', dev, '.tex')
                templates$make_ratio_table_tex(df = views$stan, dev = dev, t_lo = NA, t_hi = NA,
                    N_vals = sort(unique(views$stan$N_num)),
                    caption = fn_caption_with_timing(paste0(dev, ': best measured Stan implementation configurations within the selected sampler. ',
                                     'Each N block shows chunks/WCP, chains/total threads and iterations per second.')),
                    label = paste0('table:Stan_implementations_best_', dev), file_path = file.path(tabledir, filename),
                    configuration_summary = TRUE)
                add_asset('table', file.path('tables', filename), '', paste0('table:Stan_implementations_best_', dev))

            }
            ##
            ## Both WCP arms already have dedicated chunk-search figures above.
            stan_chunking <-  views$configurations[views$configurations$algorithm %in%
                c('AD_Stan_chunked', 'AD_Stan_tape_chunked'), , drop = FALSE]
            ##
            if (nrow(x = stan_chunking)) {

                stan_chunking <-  stan_chunking |>
                    dplyr::group_by(device, algorithm, N, n_threads, n_chains, threads_per_chain, n_iter, dataset_md5, timing_scope, timing_estimator) |>
                    dplyr::mutate(relative_seconds = time_mean / min(time_mean[time_mean > 0])) |> dplyr::ungroup()
                ##
                chunking_plot <-  ggplot2::ggplot( data = stan_chunking,
                                                   mapping = ggplot2::aes(x = N_chunks, y = relative_seconds,
                                                                          colour = N_threads, group = N_threads)) +
                    ggplot2::geom_point(size = 3) + ggplot2::geom_line() + ggplot2::theme_bw(base_size = 16) +
                    ggplot2::facet_wrap(~ device + N_label + Algorithm_label, scales = 'free_x') +
                    ggplot2::labs( x = 'Chunk setting',
                                   y = 'Elapsed seconds / best at the same chain and thread counts',
                                   colour = 'Total threads',
                                   caption = plot_timing_label)
                ##
                filename <-  'Figure_ps2_Stan_chunk_size_sensitivity_reltime.png'
                ggplot2::ggsave(filename = file.path(figdir, filename), plot = chunking_plot, width = 16, height = 12, dpi = 100)
                ##
                add_asset( kind = 'figure', file = file.path('figures', filename),
                           caption = 'Stan chunking-only sensitivity, relative to the fastest measured setting at the same chain and thread counts.',
                           label = 'figure:stan_chunk_sensitivity_FALSE')

            }

        }
        ##
        if (nrow(views$NicoStan_by_budget)) {

            templates$R_fn_plot_Stan_variants_throughput(views$NicoStan_by_budget, output_path = figdir,
                file_prefix = 'Figure_ps2_NicoStan_variants', highlight_best = TRUE)
            templates$R_fn_plot_Stan_variants_relative(views$NicoStan_by_budget, output_path = figdir,
                file_prefix = 'Figure_ps2_NicoStan_variants', baseline_label = 'NicoStan')
            add_asset('figure', 'figures/Figure_ps2_NicoStan_variants_total_throughput.png',
                paste0('Within-NicoStan comparison of the unchunked, chunking-only, WCP-only and combined configurations. ',
                       'Chunks and WCP are selected within each N and thread budget; outlined points are the best tested configuration ',
                       'for each N and mode on that device. Modes may use different iteration budgets; throughput accounts for the iterations, and the configured fixed trajectory length is recorded in study metadata.'),
                'figure:NicoStan_variants_throughput')
            add_asset('figure', 'figures/Figure_ps2_NicoStan_variants_relative_to_baseline.png',
                'Within-NicoStan throughput relative to the unchunked model at the same N and total thread budget, where measured.',
                'figure:NicoStan_variants_relative')
            ##
            for (dev in unique(views$NicoStan_best_by_N$device)) {

                filename <-  paste0('table_NicoStan_best_configurations_', dev, '.tex')
                templates$make_ratio_table_tex(df = views$NicoStan_best_by_N, dev = dev, t_lo = NA, t_hi = NA,
                    N_vals = sort(unique(views$NicoStan_best_by_N$N_num)),
                    caption = fn_caption_with_timing(paste0(dev, ': best tested NicoStan configuration for each N and mode. ',
                        'Each block shows chunks/WCP, chains/total threads and iterations per second. ',
                        'WCP-only uses target chunks equal to WCP; the combined mode searches both. Budgets remain visible.')),
                    label = paste0('table:NicoStan_best_configurations_', dev), file_path = file.path(tabledir, filename),
                    configuration_summary = TRUE,
                    row_labels = c('NicoStan', 'NicoStan-chunking', 'NicoStan-WCP', 'NicoStan-chunking_and_WCP'))
                add_asset('table', file.path('tables', filename), '', paste0('table:NicoStan_best_configurations_', dev))

            }

        }
        ##
        ## ---- Within-BayesMVP comparison figures and tables, mirroring the NicoStan block above:
        ##
        if (nrow(views$BayesMVP_by_budget)) {

            templates$R_fn_plot_Stan_variants_throughput(views$BayesMVP_by_budget, output_path = figdir,
                file_prefix = 'Figure_ps2_BayesMVP_variants', highlight_best = TRUE)
            templates$R_fn_plot_Stan_variants_relative(views$BayesMVP_by_budget, output_path = figdir,
                file_prefix = 'Figure_ps2_BayesMVP_variants', baseline_label = 'BayesMVP', baseline_mode = 'BayesMVP')
            add_asset('figure', 'figures/Figure_ps2_BayesMVP_variants_total_throughput.png',
                paste0('Within-BayesMVP comparison of the unchunked, chunking-only, WCP-only and combined chunking + WCP configurations. ',
                       'Chunks and WCP are selected within each N and thread budget; outlined points are the best tested configuration ',
                       'for each N and mode on that device. Modes may use different iteration budgets; throughput accounts for the iterations, and the fixed trajectory length is recorded in study metadata.'),
                'figure:BayesMVP_variants_throughput')
            add_asset('figure', 'figures/Figure_ps2_BayesMVP_variants_relative_to_baseline.png',
                'Within-BayesMVP throughput relative to the unchunked, one-thread-per-chain configuration at the same N and total thread budget, where measured.',
                'figure:BayesMVP_variants_relative')
            ##
            for (dev in unique(views$BayesMVP_best_by_N$device)) {

                filename <-  paste0('table_BayesMVP_best_configurations_', dev, '.tex')
                templates$make_ratio_table_tex(df = views$BayesMVP_best_by_N, dev = dev, t_lo = NA, t_hi = NA,
                    N_vals = sort(unique(views$BayesMVP_best_by_N$N_num)),
                    caption = fn_caption_with_timing(paste0(dev, ': best tested BayesMVP configuration for each N and mode. ',
                        'Each block shows chunks/WCP, chains/total threads and iterations per second. ',
                        'WCP-only uses one chunk per WCP thread; the combined mode searches both chunks and WCP. Budgets remain visible.')),
                    label = paste0('table:BayesMVP_best_configurations_', dev), file_path = file.path(tabledir, filename),
                    configuration_summary = TRUE,
                    row_labels = c('BayesMVP', 'BayesMVP-chunking', 'BayesMVP-WCP', 'BayesMVP-chunking_and_WCP'))
                add_asset('table', file.path('tables', filename), '', paste0('table:BayesMVP_best_configurations_', dev))

            }

        }
        ##
        wcp <-  views$analysis$wcp_matched
        ##
        if (nrow(wcp)) {

            wcp$panel <-  paste(wcp$device, paste0('N=', wcp$N), paste0(wcp$n_chains, ' chains'))
            ##
            panel_order <-  order(wcp$device, wcp$N, wcp$n_chains)
            ##
            wcp$panel <-  factor(wcp$panel, levels = unique(wcp$panel[panel_order]))
            ##
            p <-  ggplot2::ggplot(wcp, ggplot2::aes(factor(num_chunks), factor(threads_per_chain), fill = wcp_speedup)) +
                ggplot2::geom_tile() +
                ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(wcp_speedup), 'NA', signif(wcp_speedup, 2))), size = 2.8) +
                ggplot2::facet_wrap(~ panel, scales = 'free', ncol = 4) + ggplot2::theme_bw(base_size = 14) +
                ggplot2::scale_fill_gradient2(low = '#b35806', mid = 'white', high = '#2166ac', midpoint = 1, na.value = 'grey90') +
                ggplot2::labs(x = 'Chunks', y = 'Threads per chain', fill = 'Serial / WCP time per iteration',
                              caption = plot_timing_label) +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1), legend.position = 'bottom')
            ##
            ggplot2::ggsave( file.path(figdir, 'Figure_BayesMVP_WCP_matched_speedup.png'),
                             p,
                             width = 16,
                             height = max( 9,
                                           3 * ceiling(length(unique(wcp$panel)) / 4)),
                             dpi = 100)
            ##
            add_asset( 'figure',
                       'figures/Figure_BayesMVP_WCP_matched_speedup.png',
                       'BayesMVP WCP speedup against one thread per chain at identical chain count, chunk count, N and data. Corrected times are compared per iteration when budgets differ. Grey NA cells lack a usable matched serial baseline; blank cells were not measured. Values above one indicate lower time per iteration with WCP.',
                       'figure:bayesmvp_wcp_matched')
            ##
            best <-  wcp[order(wcp$time_mean, wcp$num_chunks, wcp$threads_per_chain), ]
            ##
            best <-  best[!duplicated(best[c('device', 'N', 'n_chains')]), ]
            ##
            fn_paper1_write_csv(best, file.path(csvdir, 'wcp_best_fixed_chains.csv'))
            ##
            lines <-  c( '\\begin{table}[H]',
                         '\\centering',
                         paste0('\\caption{', fn_caption_with_timing('Fastest measured BayesMVP WCP configuration for each fixed chain count. Speedup compares corrected time per iteration against the serial configuration with the same chunks; unavailable matched baselines are marked. Seconds retain the WCP run budget. This supplementary summary does not replace the separate optimum for each WCP count.'), '}'),
                         '\\label{table:bayesmvp_wcp_matched}',
                         '\\begin{tabular}{llrrrrr}',
                         '\\toprule',
                         'Device & $N$ & Chains & Chunks & Threads/chain & Seconds & Speedup \\\\',
                         '\\midrule')
            ##
            for (i in seq_len(nrow(best))) lines <-  c(lines, paste0( paste( c( best$device[i],
                                                                                best$N[i],
                                                                                best$n_chains[i],
                                                                                best$num_chunks[i],
                                                                                best$threads_per_chain[i],
                                                                                templates$fmt3(best$time_mean[i]),
                                                                                templates$fmt3(best$wcp_speedup[i])),
                                                                             collapse = ' & '),
                                                                      ' \\\\'))
            ##
            writeLines(c(lines, '\\bottomrule', '\\end{tabular}', '\\end{table}'), file.path(tabledir, 'table_BayesMVP_WCP_matched.tex'))
            ##
            add_asset('table', 'tables/table_BayesMVP_WCP_matched.tex', '', 'table:bayesmvp_wcp_matched')

        }
        ##
        ## ---- Within-Mplus comparison of the two long-run modes (BITERATIONS vs FBITERATIONS) on every device. Configurations
        ## ---- are matched on arm, N, chains and threads; the ratio is FBITERATIONS / BITERATIONS seconds per iteration, so
        ## ---- values above one mean FBITERATIONS was slower. The cross-software tables and figures use the first mode in the runner.
        ##
        ## ---- Each mode is summarised here from the saved runs that match the current runner settings (mean over repeats),
        ## ---- in the same form as the runner's mplus_iteration_mode_comparison.csv.
        ##
        mplus_mode_results <-  input$mplus_mode_results
        mplus_mode_rows <-  if (is.null(mplus_mode_results) || !nrow(mplus_mode_results)) NULL else
            do.call(rbind, lapply(split(mplus_mode_results, mplus_mode_results$mplus_iteration_mode), function(mode_results) {
                mode_configurations <-  private$fn_summarise_paper1_benchmark(results = mode_results)$configurations
                mode_configurations$mplus_iteration_mode <-  unique(mode_results$mplus_iteration_mode)
                mode_configurations[, !vapply(mode_configurations, is.list, logical(1)), drop = FALSE]
            }))
        ##
        if (!is.null(mplus_mode_rows)) fn_paper1_write_csv(mplus_mode_rows, file.path(csvdir, 'mplus_iteration_mode_comparison.csv'))
        ##
        if (!is.null(mplus_mode_rows) && all(c('BITERATIONS', 'FBITERATIONS') %in% mplus_mode_rows$mplus_iteration_mode)) {
            ##
            mode_columns <-  c('device', 'algorithm', 'N', 'n_chains', 'n_threads', 'n_iter', 'seconds_per_iteration_mean')
            ##
            matched_modes <-  merge( x = mplus_mode_rows[mplus_mode_rows$mplus_iteration_mode == 'BITERATIONS',  mode_columns],
                                     y = mplus_mode_rows[mplus_mode_rows$mplus_iteration_mode == 'FBITERATIONS', mode_columns],
                                     by = c('device', 'algorithm', 'N', 'n_chains', 'n_threads', 'n_iter'),
                                     suffixes = c('_BITERATIONS', '_FBITERATIONS'))
            ##
            matched_modes$ratio_FBITERATIONS_over_BITERATIONS <-  matched_modes$seconds_per_iteration_mean_FBITERATIONS /
                                                                   matched_modes$seconds_per_iteration_mean_BITERATIONS
            ##
            fn_paper1_write_csv(matched_modes, file.path(csvdir, 'mplus_iteration_modes_matched.csv'))
            ##
            mode_summary <-  do.call(rbind, lapply(split(matched_modes, list(matched_modes$device, matched_modes$algorithm, matched_modes$N), drop = TRUE), function(group) {
                    data.frame( device           = group$device[1],
                                algorithm        = group$algorithm[1],
                                N                = group$N[1],
                                n_iter           = group$n_iter[1],
                                n_configurations = nrow(group),
                                median_ratio     = stats::median(group$ratio_FBITERATIONS_over_BITERATIONS),
                                min_ratio        = min(group$ratio_FBITERATIONS_over_BITERATIONS),
                                max_ratio        = max(group$ratio_FBITERATIONS_over_BITERATIONS),
                                stringsAsFactors = FALSE)
            }))
            ##
            mode_summary <-  mode_summary[order(mode_summary$device, mode_summary$algorithm, mode_summary$N), ]
            ##
            fn_paper1_write_csv(mode_summary, file.path(csvdir, 'mplus_iteration_modes_summary.csv'))
            ##
            lines <-  c( '\\begin{table}[H]',
                         '\\centering',
                         paste0('\\caption{', fn_caption_with_timing('Mplus FBITERATIONS versus BITERATIONS, at identical requested iterations, chains, threads, N and data. The ratio is FBITERATIONS / BITERATIONS seconds per iteration for each matched configuration (the median and range across configurations are shown), so values above one mean FBITERATIONS was slower.'), '}'),
                         '\\label{table:mplus_iteration_modes}',
                         '\\begin{tabular}{llrrrrr}',
                         '\\toprule',
                         'Device & Arm & $N$ & Iterations & Configurations & Median ratio & Range \\\\',
                         '\\midrule')
            ##
            for (i in seq_len(nrow(mode_summary))) lines <-  c(lines, paste0( paste( c( mode_summary$device[i],
                                                                                        gsub('_', '\\_', mode_summary$algorithm[i], fixed = TRUE),
                                                                                        fn_paper1_format_number_commas_from_10000(mode_summary$N[i]),
                                                                                        mode_summary$n_iter[i],
                                                                                        mode_summary$n_configurations[i],
                                                                                        formatC(mode_summary$median_ratio[i], format = 'f', digits = 2),
                                                                                        paste0(formatC(mode_summary$min_ratio[i], format = 'f', digits = 2), '-',
                                                                                               formatC(mode_summary$max_ratio[i], format = 'f', digits = 2))),
                                                                                     collapse = ' & '),
                                                                              ' \\\\'))
            ##
            writeLines(c(lines, '\\bottomrule', '\\end{tabular}', '\\end{table}'), file.path(tabledir, 'table_Mplus_iteration_modes.tex'))
            ##
            add_asset('table', 'tables/table_Mplus_iteration_modes.tex', '', 'table:mplus_iteration_modes')

        }
        # Optional MP/MT comparisons retain every chunk setting and use per-repeat throughput SD.
        mpmt <-  subset(views$configurations, algorithm %in% c('MD_BayesMVP', 'MD_BayesMVP_multi_process'))
        ##
        if (any(mpmt$algorithm == 'MD_BayesMVP_multi_process')) {

            mpmt$mode <-  ifelse(mpmt$algorithm == 'MD_BayesMVP', 'MT', 'MP')
            ##
            mpmt$chunk_label <-  paste('Chunks =', mpmt$num_chunks)
            ##
            mpmt$chain_rate_SD <-  vapply( mpmt$case_ids,
                                                   function(ids) {

                    rows <-  input$results[match(ids, input$results$case_id), ]
                    ##
                    stats::sd(rows$n_chains / rows$elapsed_seconds)

            },
                                                   numeric(1))
            ##
            mpmt <-  mpmt |> dplyr::group_by(device, N, mode, chunk_label) |>
                dplyr::mutate(norm_scaling = chain_rate / chain_rate[which.min(n_threads)]) |> dplyr::ungroup()
            ##
            templates$R_fn_plot_mpmt_throughput(mpmt, output_path = figdir)
            ##
            templates$R_fn_plot_mpmt_scaling(mpmt, output_path = figdir)
            ##
            for (type in c('throughput', 'scaling')) add_asset( 'figure',
                                                                paste0('figures/Figure_BayesMVP_MP_vs_MT_', type, '.png'),
                                                                paste( 'Optional multiprocessing versus multithreading:',
                                                                       type,
                                                                       'by chunk setting.'),
                                                                paste0( 'figure:mpmt_',
                                                                        type))

        }
        ##
        writeLines(notes, file.path(output_dir, 'READ_ME.txt'))
        ##
        fn_paper1_write_csv(assets, file.path(output_dir, 'asset_manifest.csv'))
        ##
        tex <-  c('% Generated from explicitly selected saved studies. Review before importing into the draft.')
        ##
        for (i in seq_len(nrow(assets))) {

            row <-  assets[i, ]
            ##
            tex <-  c(tex, if (row$kind == 'table') paste0('\\input{', row$file, '}') else
                      c('\\begin{figure}[H]', '\\centering', paste0( '\\includegraphics[width=\\textwidth,height=0.8\\textheight,keepaspectratio]{',
                                                                     row$file,
                                                                     '}'),
                        paste0('\\caption{', row$caption, '}'), paste0('\\label{', row$label, '}'), '\\end{figure}'), '\\clearpage')

        }
        ##
        writeLines(tex, file.path(output_dir, 'generated_assets.tex'))
        ## Portable bundles for the new manuscript sections, using the same figures and original-style tables.
        ##
        ## ---- Section membership is decided on the asset FILE NAME, anchored at its start:
        ##
        ## The unanchored pattern 'Stan_variants' also matches 'Figure_ps2_NicoStan_variants_...', which would
        ## place those figures in both bundles and define their LaTeX labels twice in Main.tex.
        ##
        section_patterns <-  c(Stan_implementations = '^(Figure_ps2_Stan_variants_|table_Stan_implementations_best_|Figure_ps2_Stan_chunk_size_)',
                                NicoStan_comparisons = '^(Figure_ps2_NicoStan_variants_|table_NicoStan_best_configurations_)',
                                ## Within-BayesMVP comparison: the four-mode figures and best-configuration table.
                                BayesMVP_comparisons = '^(Figure_ps2_BayesMVP_variants_|table_BayesMVP_best_configurations_)',
                                ## How the BayesMVP chunk and WCP counts were chosen: per-N chunk-search and optimal-chunks-by-WCP
                                ## figures and the best-chunks-by-WCP tables (file names start with the device, e.g. HPC_paper_1_N_500_...).
                                BayesMVP_chunk_WCP_selection = '^[A-Za-z]+_paper_1_N_[0-9]+_algorithm_MD_BayesMVP_WCP_(chunk_search|optimal_chunks_by_WCP|best_chunks_by_WCP)',
                                ## The same selection figures and tables for the Stan (NicoStan) WCP arm.
                                Stan_chunk_WCP_selection = '^[A-Za-z]+_paper_1_N_[0-9]+_algorithm_AD_Stan_WCP_(chunk_search|optimal_chunks_by_WCP|best_chunks_by_WCP)',
                                ## Experiment 1 (optimal chunks) and experiment 2 (cross-algorithm scaling) assets, so that every manuscript
                                ## table and figure is \input from an R-generated bundle instead of being pasted into Main.tex.
                                Exp1_chunking = '^(Figure_N_chunks_pilot_study_|table_ps1_)',
                                Exp2_scaling = '^(Figure_ps2_plot_2_adj_scalability|Figure_paper1_serial_parallel_efficiency|table_ps2_|[A-Za-z]+_paper_1_N_[0-9]+_algorithm_comparison_scaling)')
        ##
        asset_file_names <-  basename(assets$file)
        ##
        section_membership_matrix <-  vapply( X         = section_patterns,
                                              FUN       = function(section_pattern) grepl(section_pattern, asset_file_names),
                                              FUN.VALUE = logical(length(asset_file_names)))
        ##
        section_membership_matrix <-  matrix(section_membership_matrix, nrow = length(asset_file_names))
        ##
        assets_in_more_than_one_section <-  asset_file_names[rowSums(section_membership_matrix) > 1]
        ##
        if (length(assets_in_more_than_one_section)) {

            stop('Generated assets matched more than one manuscript section: ', paste(assets_in_more_than_one_section, collapse = ', '))

        }
        ##
        validation_only <-  any(vapply(input$studies, function(study) isTRUE(study$metadata$validation_only), logical(1)))
        for (section_name in names(section_patterns)) {

            section_assets <-  assets[grepl(section_patterns[[section_name]], asset_file_names), , drop = FALSE]
            ##
            section_labels <-  section_assets$label[nzchar(section_assets$label)]
            ##
            if (anyDuplicated(section_labels)) stop('Duplicate LaTeX labels in generated section ', section_name, ': ', paste(unique(section_labels[duplicated(section_labels)]), collapse = ', '))
            ##
            if (!nrow(section_assets)) next
            bundle <-  file.path(output_dir, 'paper_sections', section_name)
            dir.create(file.path(bundle, 'figures'), recursive = TRUE, showWarnings = FALSE)
            dir.create(file.path(bundle, 'tables'), recursive = TRUE, showWarnings = FALSE)
            section_tex <-  c('% Generated from the selected completed study; see its report manifest.')
            for (row_index in seq_len(nrow(section_assets))) {

                row <-  section_assets[row_index, ]
                file.copy(file.path(output_dir, row$file), file.path(bundle, row$file), overwrite = TRUE)
                manuscript_path <-  paste0('Files/Generated/', section_name, '/', row$file)
                section_tex <-  c(section_tex, if (row$kind == 'table') paste0('\\input{', manuscript_path, '}') else
                    c('\\begin{figure}[H]', '\\centering',
                      paste0('\\includegraphics[width=\\textwidth,height=0.8\\textheight,keepaspectratio]{', manuscript_path, '}'),
                      paste0('\\caption{', row$caption, '}'), paste0('\\label{', row$label, '}'), '\\end{figure}'))

            }
            writeLines(section_tex, file.path(bundle, 'section.tex'))
            if (!is.null(manuscript_dir) && !validation_only) {

                destination <-  file.path(manuscript_dir, 'Files', 'Generated', section_name)
                dir.create(destination, recursive = TRUE, showWarnings = FALSE)
                ##
                ## Report files left by earlier exports without deleting or importing them.
                ##
                files_no_longer_generated <-  setdiff(list.files(destination, recursive = TRUE),
                                                       list.files(bundle, recursive = TRUE))
                ##
                if (length(files_no_longer_generated)) {

                    message(paste0('\033[36m', 'Files in ', destination,
                                   ' that this export does not write (stale; not \\input): ',
                                   paste(files_no_longer_generated, collapse = ', '), '\033[0m'))

                }
                for (file in list.files(bundle, recursive = TRUE)) {

                    dir.create(dirname(file.path(destination, file)), recursive = TRUE, showWarnings = FALSE)
                    if (!file.copy(file.path(bundle, file), file.path(destination, file), overwrite = TRUE)) {

                        stop('Could not copy generated section asset: ', file)

                    }

                }

            }

        }
        ##
        writeLines( c('\\documentclass{article}',
                      '\\usepackage[a4paper,landscape,margin=1cm]{geometry}',
                      '\\usepackage{graphicx,float,booktabs,amsmath,longtable}',
                      '\\begin{document}', '\\section*{Paper 1 generated results: review copy}',
                      'These figures and tables use the saved study directories recorded in the report manifest. Missing observations are not filled from older studies. Consult READ\\_ME.txt and coverage.csv before manuscript use.',
                      '\\input{generated_assets.tex}', '\\end{document}'),
                    file.path(output_dir, 'Results_preview.tex'))
        ##
        saveRDS( list( created_at = Sys.time(),
                       timing_estimators = unique(input$results$timing_estimator),
                       source_dirs = input$paths,
                       metadata = lapply(input$studies, `[[`, 'metadata'),
                       input_md5 = lapply(input$studies, `[[`, 'input_md5'),
                       coverage = coverage,
                       notes = notes,
                       report_code_md5 = tools::md5sum(file.path( helper_dir,
                                                                  c( 'R_fns_alg_paper_1_chunking_WCP_par_scaling.R',
                                                                     'R_fns_alg_paper_1_figures_tables.R',
                                                                     'R_fns_alg_paper_1_presentation_templates.R'))),
                       assets = assets),
                 file.path(output_dir, 'report_manifest.rds'))
        ##
        message('Figures, tables and provenance saved in: ', output_dir)
        ##
        invisible(list(output_dir = output_dir, assets = assets, coverage = coverage, views = views))

}






















