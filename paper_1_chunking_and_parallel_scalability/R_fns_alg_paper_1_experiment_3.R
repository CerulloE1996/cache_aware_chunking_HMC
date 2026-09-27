##
## ======================================================================================================================================
## R_fns_alg_paper_1_experiment_3.R
##
## Paper 1, Experiment 3 (absolute efficiency): time to reach a target minimum ESS, for BayesMVP vs Stan (via NicoStan) vs
## Mplus, read from Paper 2's pilot studies (ps5 = Mplus, ps6 = Stan, ps7 = BayesMVP). This file never runs sampling; it
## only reads saved CSV/RDS summaries and writes a LaTeX table bundle, following the same paper_sections/<bundle>/
## {tables/*.tex, section.tex} convention as R_fns_alg_paper_1_figures_tables.R, so it \input's into Main.tex the same way.
##
## Per-run definition:
##   time_to_target_seconds = time_burnin_seconds + (time_sampling_seconds + time_summary_seconds) * (target_min_ESS / min_ESS_over_parameters_of_interest)
## then averaged (mean) over runs, separately for each (device, N, software); time_summary_seconds is the posterior-summary
## time where it was recorded (ps7) and zero otherwise. min_ESS_over_parameters_of_interest is the minimum
## ESS over Se / Sp / prevalence (never the raw parameter block), matching ps7's own "min_ESS" definition.
##
## Assumptions baked into the per-source readers below (none of these were spelled out in the source CSV column names,
## so they are stated here explicitly and repeated in the runner's printed self-check):
##   - ps5 (Mplus): time_sampling_seconds = time_sampling_mins * 60 when that column exists, else min_ESS / ESS_per_sec_sampling.
##                  time_burnin_seconds = time_total - time_sampling_seconds.
##   - ps6 (Stan):  time_burnin_seconds = the saved time_burnin column; time_sampling_seconds = the saved time_pb column
##                  ("post burn-in" sampling time). time_total in these files also includes model compilation and other
##                  fixed overhead (it is NOT time_burnin + time_pb + time_sampling for small N, where compilation
##                  dominates), so time_total is deliberately never used for the time-to-target formula. The saved per-run
##                  CSVs carry no time_summaries column, so time_pb alone is the sampling time.
##   - ps7 (BayesMVP): read via fn_ps7_summarise_directory(); the configuration ranked first (fastest to target) within
##                  each N is selected, and its saved RUNS (ps7_summary$runs: time_burnin, time_sampling, time_summaries,
##                  min_ESS per run) become one row each, so the mean over runs in fn_paper1_exp3_time_to_target()
##                  reproduces ps7's own mean_time_to_target_min_ESS for that configuration.
##
## ---- Software row labels used throughout (fixed, so every table and every join uses the identical strings) --------------------------
##
fn_paper1_exp3_software_labels <- function() {

        ## return(c( bayesmvp = "BayesMVP",
        return(c( bayesmvp = "NicoStan+BayesMVP",
                  stan     = "Stan (NUTS)",
                  mplus    = "Mplus (PX-Gibbs)"))

}


##
## ---- ps5 (Mplus): one CSV file -> one row per run ------------------------------------------------------------------------------------
##
#' @param ps5_csv_file_path   path to one ps5 "full_results" CSV (one or more N values in one file).
#' @param device_label        "HPC" or "Laptop" (or another explicit label); the caller decides this, it is never guessed
#'                             from the file path.
#' @return tibble, one row per run: device, software_label, N, run_index, time_burnin_seconds, time_sampling_seconds,
#'         min_ESS_over_parameters_of_interest, max_Rhat, source_file.
fn_paper1_exp3_rows_from_ps5_csv <- function( ps5_csv_file_path,
                                              device_label) {

        if (!isTRUE(is.character(device_label)) || length(device_label) != 1 || !nzchar(device_label)) {
              stop("fn_paper1_exp3_rows_from_ps5_csv: device_label must be one non-empty string.")
        }
        if (!file.exists(ps5_csv_file_path)) stop("fn_paper1_exp3_rows_from_ps5_csv: file does not exist: ", ps5_csv_file_path)
        ##
        ps5_raw_table <- utils::read.csv(file = ps5_csv_file_path, stringsAsFactors = FALSE)
        ##
        ##
        ## ---- required columns; stop() loudly rather than silently producing NA time-to-target rows later:
        ##
        required_column_names <- c("N", "run", "time_total", "max_Rhat", "min_ESS")
        missing_column_names   <- setdiff(required_column_names, names(ps5_raw_table))
        if (length(missing_column_names) > 0) {
              stop("fn_paper1_exp3_rows_from_ps5_csv: ", ps5_csv_file_path, " is missing column(s) ",
                   paste(missing_column_names, collapse = ", "))
        }
        ##
        ## ---- sampling seconds: time_sampling_mins * 60 when present, else back it out from ESS_per_sec_sampling:
        ##
        {
              if ("time_sampling_mins" %in% names(ps5_raw_table)) {
                    time_sampling_seconds_vector <- ps5_raw_table$time_sampling_mins * 60
                    sampling_seconds_source       <- "time_sampling_mins * 60"
              } else if ("ESS_per_sec_sampling" %in% names(ps5_raw_table)) {
                    time_sampling_seconds_vector <- ps5_raw_table$min_ESS / ps5_raw_table$ESS_per_sec_sampling
                    sampling_seconds_source       <- "min_ESS / ESS_per_sec_sampling"
              } else {
                    stop("fn_paper1_exp3_rows_from_ps5_csv: ", ps5_csv_file_path,
                         " has neither time_sampling_mins nor ESS_per_sec_sampling; cannot recover sampling seconds.")
              }
              ##
              time_burnin_seconds_vector <- ps5_raw_table$time_total - time_sampling_seconds_vector
              ##
              if (any(!is.finite(time_burnin_seconds_vector))) {
                    stop("fn_paper1_exp3_rows_from_ps5_csv: non-finite burn-in seconds recovered from ", ps5_csv_file_path,
                         " (sampling seconds source: ", sampling_seconds_source, ").")
              }
              if (any(time_burnin_seconds_vector < -1e-6)) {
                    stop("fn_paper1_exp3_rows_from_ps5_csv: negative burn-in seconds recovered from ", ps5_csv_file_path,
                         " (time_total < recovered sampling seconds); inspect the source CSV before trusting this file.")
              }
              time_burnin_seconds_vector <- pmax(time_burnin_seconds_vector, 0)
        }
        ##
        exp3_rows_from_ps5 <- tibble::tibble( device                               = device_label,
                                              software_label                       = unname(fn_paper1_exp3_software_labels()["mplus"]),
                                              N                                     = ps5_raw_table$N,
                                              run_index                            = ps5_raw_table$run,
                                              time_burnin_seconds                  = time_burnin_seconds_vector,
                                              time_sampling_seconds                = time_sampling_seconds_vector,
                                              min_ESS_over_parameters_of_interest  = ps5_raw_table$min_ESS,
                                              max_Rhat                             = ps5_raw_table$max_Rhat,
                                              source_file                          = basename(ps5_csv_file_path))
        ##
        return(exp3_rows_from_ps5)

}


##
## ---- ps6 (Stan via NicoStan): one or more CSV files (one per N) -> one row per run --------------------------------------------------
##
#' @param ps6_csv_file_paths  character vector of ps6 "full_results" CSV paths (typically one file per N).
#' @param device_label        explicit device label, as for fn_paper1_exp3_rows_from_ps5_csv().
#' @param select_fastest_burnin  if TRUE (and the CSVs carry an n_burnin column), keep only the runs of the burn-in length
#'                             with the lowest mean per-run time to the target ESS at each N - the same selection rule used
#'                             for BayesMVP (fastest configuration by mean time to target over its runs).
#' @param target_min_ESS_by_N  named numeric vector (names = N as character); required when select_fastest_burnin = TRUE.
#'                             When the CSVs also carry an arm column (all-arms production CSVs), the selection is over
#'                             (arm, burn-in length) cells at each N instead of burn-in length alone.
#' @return tibble, same columns as fn_paper1_exp3_rows_from_ps5_csv(), software_label = "Stan (NUTS)".
#'         When the CSVs carry a time_summaries column (the production PS6 runs), it is returned as time_summary_seconds and
#'         enters the time-to-target formula; older CSVs without it contribute zero, as before.
fn_paper1_exp3_rows_from_ps6_csv <- function( ps6_csv_file_paths,
                                              device_label,
                                              select_fastest_burnin = FALSE,
                                              target_min_ESS_by_N   = NULL) {

        if (!isTRUE(is.character(device_label)) || length(device_label) != 1 || !nzchar(device_label)) {
              stop("fn_paper1_exp3_rows_from_ps6_csv: device_label must be one non-empty string.")
        }
        if (length(ps6_csv_file_paths) == 0) stop("fn_paper1_exp3_rows_from_ps6_csv: ps6_csv_file_paths is empty.")
        missing_files <- ps6_csv_file_paths[!file.exists(ps6_csv_file_paths)]
        if (length(missing_files) > 0) stop("fn_paper1_exp3_rows_from_ps6_csv: file(s) do not exist: ", paste(missing_files, collapse = ", "))
        ##
        required_column_names <- c("N", "run", "time_burnin", "time_pb", "max_Rhat", "min_ESS")
        ##
        exp3_rows_from_ps6_list <- lapply( ps6_csv_file_paths,
                                           function(one_ps6_csv_file_path) {

              ps6_raw_table <- utils::read.csv(file = one_ps6_csv_file_path, stringsAsFactors = FALSE)
              ##
              missing_column_names <- setdiff(required_column_names, names(ps6_raw_table))
              if (length(missing_column_names) > 0) {
                    stop("fn_paper1_exp3_rows_from_ps6_csv: ", one_ps6_csv_file_path, " is missing column(s) ",
                         paste(missing_column_names, collapse = ", "))
              }
              ##
              ## time_pb ("post burn-in") is treated as the sampling-phase time; time_total is NOT used here because it
              ## also includes model compilation and other fixed overhead not part of the sampling process itself - see
              ## the file header comment above.
              ##
              if (any(!is.finite(ps6_raw_table$time_burnin)) || any(ps6_raw_table$time_burnin < 0) ||
                  any(!is.finite(ps6_raw_table$time_pb))     || any(ps6_raw_table$time_pb     < 0)) {
                    stop("fn_paper1_exp3_rows_from_ps6_csv: non-finite or negative time_burnin/time_pb in ", one_ps6_csv_file_path)
              }
              ##
              ps6_rows_one_file <- tibble::tibble( device                               = device_label,
                              software_label                       = unname(fn_paper1_exp3_software_labels()["stan"]),
                              N                                     = ps6_raw_table$N,
                              run_index                            = ps6_raw_table$run,
                              time_burnin_seconds                  = ps6_raw_table$time_burnin,
                              time_sampling_seconds                = ps6_raw_table$time_pb,
                              min_ESS_over_parameters_of_interest  = ps6_raw_table$min_ESS,
                              max_Rhat                             = ps6_raw_table$max_Rhat,
                              source_file                          = basename(one_ps6_csv_file_path))
              ##
              ## ---- production PS6 CSVs (ps6_BIN_Stan_production_full_results_N*.csv) carry the measured cmdstanr
              ##      post-processing time (time_summaries) and the burn-in length (n_burnin):
              ##
              if ("time_summaries" %in% names(ps6_raw_table)) {
                    if (any(!is.finite(ps6_raw_table$time_summaries)) || any(ps6_raw_table$time_summaries < 0)) {
                          stop("fn_paper1_exp3_rows_from_ps6_csv: non-finite or negative time_summaries in ", one_ps6_csv_file_path)
                    }
                    ps6_rows_one_file$time_summary_seconds <- ps6_raw_table$time_summaries
              }
              if ("n_burnin" %in% names(ps6_raw_table)) ps6_rows_one_file$n_burnin <- ps6_raw_table$n_burnin
              ##
              ## ---- all-arms production PS6 CSVs (ps6_BIN_Stan_production_ALL_ARMS_full_results_N*.csv) also carry the arm
              ##      (config_id | burn-in model | warm-up chains | sampling chains), so the selection runs over (arm, burn-in length):
              ##
              if ("arm" %in% names(ps6_raw_table)) ps6_rows_one_file$arm <- ps6_raw_table$arm
              ##
              ps6_rows_one_file

        })
        ##
        exp3_rows_from_ps6 <- dplyr::bind_rows(exp3_rows_from_ps6_list)
        ##
        ## ---- optional: keep only the burn-in length with the lowest mean per-run time to the target ESS at each N
        ##      (same per-run formula as fn_paper1_exp3_time_to_target()):
        ##
        if (isTRUE(select_fastest_burnin)) {

              if (!("n_burnin" %in% names(exp3_rows_from_ps6))) {
                    stop("fn_paper1_exp3_rows_from_ps6_csv: select_fastest_burnin = TRUE needs an n_burnin column in the ps6 CSV(s).")
              }
              if (is.null(target_min_ESS_by_N) || is.null(names(target_min_ESS_by_N))) {
                    stop("fn_paper1_exp3_rows_from_ps6_csv: select_fastest_burnin = TRUE needs target_min_ESS_by_N (named by N).")
              }
              N_without_target <- setdiff(unique(exp3_rows_from_ps6$N), as.numeric(names(target_min_ESS_by_N)))
              if (length(N_without_target) > 0) {
                    stop("fn_paper1_exp3_rows_from_ps6_csv: no target_min_ESS for N = ", paste(N_without_target, collapse = ", "))
              }
              ##
              summary_seconds <- if ("time_summary_seconds" %in% names(exp3_rows_from_ps6)) exp3_rows_from_ps6$time_summary_seconds else 0
              target_min_ESS  <- unname(target_min_ESS_by_N[as.character(exp3_rows_from_ps6$N)])
              exp3_rows_from_ps6$time_to_target_for_selection <- exp3_rows_from_ps6$time_burnin_seconds +
                  (exp3_rows_from_ps6$time_sampling_seconds + summary_seconds) * (target_min_ESS / exp3_rows_from_ps6$min_ESS_over_parameters_of_interest)
              ##
              ## burnin_means <- exp3_rows_from_ps6 |>
              ##     dplyr::group_by(.data$N, .data$n_burnin) |>
              ##     dplyr::summarise( n_runs                      = dplyr::n(),
              ##                       mean_time_to_target_seconds = mean(.data$time_to_target_for_selection),
              ##                       .groups                     = "drop") |>
              ##     dplyr::group_by(.data$N) |>
              ##     dplyr::mutate(selected = .data$mean_time_to_target_seconds == min(.data$mean_time_to_target_seconds)) |>
              ##     dplyr::ungroup()
              ##
              ## ---- cells = (N, arm, n_burnin) when the CSVs carry an arm column, (N, n_burnin) otherwise; the selected cell at
              ##      each N is the one with the lowest mean per-run time to the target ESS:
              ##
              selection_cell_columns <- if ("arm" %in% names(exp3_rows_from_ps6)) c("N", "arm", "n_burnin") else c("N", "n_burnin")
              ##
              burnin_means <- exp3_rows_from_ps6 |>
                  dplyr::group_by(dplyr::across(dplyr::all_of(selection_cell_columns))) |>
                  dplyr::summarise( n_runs                      = dplyr::n(),
                                    mean_time_to_target_seconds = mean(.data$time_to_target_for_selection),
                                    .groups                     = "drop") |>
                  dplyr::group_by(.data$N) |>
                  dplyr::mutate(selected = .data$mean_time_to_target_seconds == min(.data$mean_time_to_target_seconds)) |>
                  dplyr::ungroup()
              ##
              for (row_index in seq_len(nrow(burnin_means))) {
                    ## message(paste0("\033[36m", "ps6 (Stan) N = ", burnin_means$N[row_index], ", n_burnin = ", burnin_means$n_burnin[row_index],
                    ##                ": mean time to target = ", formatC(burnin_means$mean_time_to_target_seconds[row_index], format = "f", digits = 2),
                    ##                " s over ", burnin_means$n_runs[row_index], " runs",
                    ##                if (burnin_means$selected[row_index]) " (selected)" else "", "\033[0m"))
                    message(paste0("\033[36m", "ps6 (Stan) N = ", burnin_means$N[row_index],
                                   if ("arm" %in% names(burnin_means)) paste0(", arm = ", burnin_means$arm[row_index]) else "",
                                   ", n_burnin = ", burnin_means$n_burnin[row_index],
                                   ": mean time to target = ", formatC(burnin_means$mean_time_to_target_seconds[row_index], format = "f", digits = 2),
                                   " s over ", burnin_means$n_runs[row_index], " runs",
                                   if (burnin_means$selected[row_index]) " (selected)" else "", "\033[0m"))
              }
              ##
              ## selected_pairs <- burnin_means[burnin_means$selected, c("N", "n_burnin")]
              ## exp3_rows_from_ps6 <- dplyr::semi_join(x = exp3_rows_from_ps6, y = selected_pairs, by = c("N", "n_burnin"))
              selected_pairs <- burnin_means[burnin_means$selected, selection_cell_columns]
              exp3_rows_from_ps6 <- dplyr::semi_join(x = exp3_rows_from_ps6, y = selected_pairs, by = selection_cell_columns)
              exp3_rows_from_ps6$time_to_target_for_selection <- NULL

        }
        ##
        return(exp3_rows_from_ps6)

}


##
## ---- ps7 (BayesMVP): a directory of saved runs -> one row per N (the fastest-to-target configuration) -------------------------------
##
#' @param ps7_runs_directory        directory of saved ps7_run_* RDS fits for one device (fn_ps7_summarise_directory reads
#'                                   it non-recursively, so pass the folder that directly contains the run files).
#' @param device_label               explicit device label, as above.
#' @param target_min_ESS_by_N        named numeric vector, names = N as character (e.g. c("500" = 7000, "2500" = 2500)).
#' @param ps7_functions_directory    the ps_7_basic_MCMC_settings_BayesMVP/functions directory; every non-test .R file in
#'                                    it is sourced before extraction (fn_ps7_summarise_runs() needs R_fn_ps7_get_target_min_ESS
#'                                    and fn_ps7_parse_run_name() from these files; fn_ps7_extract_tau_sweep() itself pulls
#'                                    the nested-R-hat helpers from the installed BayesMVP package).
#' @param ps7_expected_run_indices  the repeat-run indices a configuration must have ALL of (default runs 1-3) to be eligible
#'                                    for selection; configurations with any of them missing are not ranked.
#' @return tibble, same columns as fn_paper1_exp3_rows_from_ps5_csv() (run_index is always NA: this is a per-configuration
#'         mean over that configuration's repeat runs, not a single run), software_label = "BayesMVP".
fn_paper1_exp3_rows_from_ps7_directory <- function( ps7_runs_directory,
                                                    device_label,
                                                    target_min_ESS_by_N,
                                                    ps7_functions_directory,
                                                    ps7_file_pattern = "^ps7_run_",
                                                    ps7_expected_run_indices = c(1, 2, 3)) {

        if (!isTRUE(is.character(device_label)) || length(device_label) != 1 || !nzchar(device_label)) {
              stop("fn_paper1_exp3_rows_from_ps7_directory: device_label must be one non-empty string.")
        }
        if (!dir.exists(ps7_runs_directory)) stop("fn_paper1_exp3_rows_from_ps7_directory: directory does not exist: ", ps7_runs_directory)
        if (is.null(names(target_min_ESS_by_N)) || any(!nzchar(names(target_min_ESS_by_N)))) {
              stop("fn_paper1_exp3_rows_from_ps7_directory: target_min_ESS_by_N must be a numeric vector named by N (e.g. c(\"2500\" = 2500)).")
        }
        if (!dir.exists(ps7_functions_directory)) stop("fn_paper1_exp3_rows_from_ps7_directory: functions directory does not exist: ", ps7_functions_directory)
        ##
        ##
        ## ---- fn_ps7_extract_tau_sweep.R uses the bare %>% pipe, so dplyr must be ATTACHED, not just namespaced:
        ##
        if (!requireNamespace("dplyr", quietly = TRUE)) stop("fn_paper1_exp3_rows_from_ps7_directory: package dplyr is required.")
        library(dplyr)
        ##
        ## ---- source every ps7 helper file (skip test_*.R, which call testthat and are not meant to be sourced):
        ##
        {
              ps7_function_file_paths <- list.files( path       = ps7_functions_directory,
                                                     pattern    = "\\.R$",
                                                     full.names = TRUE)
              ps7_function_file_paths <- ps7_function_file_paths[!grepl("^test_", basename(ps7_function_file_paths))]
              if (length(ps7_function_file_paths) == 0) stop("fn_paper1_exp3_rows_from_ps7_directory: no .R files found in ", ps7_functions_directory)
              for (one_ps7_function_file_path in ps7_function_file_paths) source(one_ps7_function_file_path, local = TRUE)
              if (!exists("fn_ps7_summarise_directory", mode = "function")) {
                    stop("fn_paper1_exp3_rows_from_ps7_directory: fn_ps7_summarise_directory() was not defined after sourcing ", ps7_functions_directory)
              }
        }
        ##
        message(paste0("\033[36m", "ps7: extracting saved runs from ", ps7_runs_directory,
                       " (this reads every ps7_run_* RDS file in that directory once; it can take a minute)", "\033[0m"))
        ##
        ## ---- ps7_file_pattern restricts the saved runs read (e.g. to one N when a directory also holds other N values):
        ps7_summary <- fn_ps7_summarise_directory( dir              = ps7_runs_directory,
                                                   pattern          = ps7_file_pattern,
                                                   target_min_ESS   = target_min_ESS_by_N,
                                                   sort_by          = "time_to_target",
                                                   sort_within_N    = TRUE,
                                                   verbose          = TRUE)
        ##
        ps7_configurations <- ps7_summary$configurations
        if (nrow(ps7_configurations) == 0) stop("fn_paper1_exp3_rows_from_ps7_directory: fn_ps7_summarise_directory() returned no configurations for ", ps7_runs_directory)
        ##
        ## ---- N is only a COLUMN of $configurations when it varies across the extracted runs (fn_ps7_summarise_runs()
        ##      labels rows only by settings that vary; a directory holding a single N puts N in $constant_settings instead):
        ##
        if (!("N" %in% names(ps7_configurations))) {
              constant_N <- ps7_summary$constant_settings$N
              if (is.null(constant_N) || length(constant_N) != 1 || !is.finite(constant_N)) {
                    stop("fn_paper1_exp3_rows_from_ps7_directory: N is neither a varying nor a constant setting in ", ps7_runs_directory)
              }
              ps7_configurations$N <- constant_N
        }
        ##
        N_found        <- sort(unique(ps7_configurations$N))
        N_with_target  <- as.numeric(names(target_min_ESS_by_N))
        N_without_target <- setdiff(N_found, N_with_target)
        if (length(N_without_target) > 0) {
              stop("fn_paper1_exp3_rows_from_ps7_directory: no target_min_ESS entry for N = ", paste(N_without_target, collapse = ", "),
                   " found in ", ps7_runs_directory, "; add it to target_min_ESS_by_N.")
        }
        message(paste0("\033[36m", "ps7: N values found in ", ps7_runs_directory, ": ", paste(N_found, collapse = ", "), "\033[0m"))
        ##
        ## best_configuration_per_N <- dplyr::filter(.data = ps7_configurations, .data$rank == 1)
        ##
        ## ---- Only COMPLETE configurations are eligible: every expected run index (ps7_expected_run_indices, runs 1-3) must be
        ##      present among that configuration's saved runs. A configuration with fewer runs (e.g. one still running) is left out,
        ##      and the fastest-to-target configuration is then re-ranked within N among the complete ones only.
        ##
        {
              ps7_runs_for_completeness <- ps7_summary$runs
              if (is.null(ps7_runs_for_completeness) || !all(c("configuration", "run") %in% names(ps7_runs_for_completeness))) {
                    stop("fn_paper1_exp3_rows_from_ps7_directory: fn_ps7_summarise_directory()$runs lacks the configuration / run columns.")
              }
              run_indices_by_configuration <- split(as.numeric(ps7_runs_for_completeness$run), ps7_runs_for_completeness$configuration)
              complete_configuration_keys  <- names(run_indices_by_configuration)[vapply(X         = run_indices_by_configuration,
                                                                                         FUN       = function(run_indices) all(ps7_expected_run_indices %in% run_indices),
                                                                                         FUN.VALUE = logical(1))]
              ##
              complete_configurations <- ps7_configurations[ps7_configurations$group_key %in% complete_configuration_keys &
                                                            is.finite(ps7_configurations$mean_time_to_target_min_ESS), , drop = FALSE]
              message(paste0("\033[36m", "ps7: ", nrow(complete_configurations), " of ", nrow(ps7_configurations),
                             " configuration(s) have all expected runs (", paste(ps7_expected_run_indices, collapse = ", "),
                             ") and a finite time to target; only these are ranked.", "\033[0m"))
              N_without_complete_configuration <- setdiff(N_found, unique(complete_configurations$N))
              if (length(N_without_complete_configuration) > 0) {
                    message(paste0("\033[36m", "ps7: no complete configuration for N = ", paste(N_without_complete_configuration, collapse = ", "),
                                   "; that N is left out of the BayesMVP rows.", "\033[0m"))
              }
              if (nrow(complete_configurations) == 0) stop("fn_paper1_exp3_rows_from_ps7_directory: no complete configuration in ", ps7_runs_directory)
              ##
              complete_configurations$rank_among_complete <- stats::ave(complete_configurations$mean_time_to_target_min_ESS,
                                                                         as.character(complete_configurations$N),
                                                                         FUN = function(values) rank(values, ties.method = "min"))
              best_configuration_per_N <- dplyr::filter(.data = complete_configurations, .data$rank_among_complete == 1)
              best_configuration_per_N <- best_configuration_per_N[!duplicated(as.character(best_configuration_per_N$N)), , drop = FALSE]
              ##
              for (row_index in seq_len(nrow(best_configuration_per_N))) {
                    message(paste0("\033[32m", "ps7: selected for N = ", best_configuration_per_N$N[row_index], ": ",
                                   best_configuration_per_N$group_key[row_index], " (", best_configuration_per_N$n_runs[row_index],
                                   " runs, mean time to target = ", formatC(best_configuration_per_N$mean_time_to_target_min_ESS[row_index],
                                                                              format = "f", digits = 2), " s)", "\033[0m"))
              }
        }
        ##
        ## ---- One row per SAVED RUN of the fastest configuration at each N (ps7_summary$runs; runs$configuration is the
        ##      configurations$group_key), so that the time to the target ESS is the mean of the per-run estimates and includes
        ##      the posterior-summary time, exactly as ps7's own time_to_target_min_ESS does.
        ps7_runs <- ps7_summary$runs
        if (is.null(ps7_runs) || !all(c("configuration", "time_burnin", "time_sampling", "time_summaries", "min_ESS") %in% names(ps7_runs))) {
              stop("fn_paper1_exp3_rows_from_ps7_directory: fn_ps7_summarise_directory()$runs lacks the per-run timing columns.")
        }
        if (!("N" %in% names(ps7_runs))) ps7_runs$N <- ps7_summary$constant_settings$N
        best_runs <- ps7_runs[ps7_runs$configuration %in% best_configuration_per_N$group_key, , drop = FALSE]
        if (nrow(best_runs) == 0) stop("fn_paper1_exp3_rows_from_ps7_directory: no saved runs match the best configuration(s).")
        ##
        exp3_rows_from_ps7 <- tibble::tibble( device                               = device_label,
                                              software_label                       = unname(fn_paper1_exp3_software_labels()["bayesmvp"]),
                                              N                                     = best_runs$N,
                                              run_index                            = if ("run" %in% names(best_runs)) as.numeric(best_runs$run) else NA_real_,
                                              time_burnin_seconds                  = best_runs$time_burnin,
                                              time_sampling_seconds                = best_runs$time_sampling,
                                              time_summary_seconds                 = best_runs$time_summaries,
                                              min_ESS_over_parameters_of_interest  = best_runs$min_ESS,
                                              max_Rhat                             = if ("max_Rhat" %in% names(best_runs)) best_runs$max_Rhat else NA_real_,
                                              source_file                          = paste0("ps7 best configuration run: ", best_runs$configuration))
        ##
        return(exp3_rows_from_ps7)

}


##
## ---- combine per-run rows from any/all of the three sources into the time-to-target summary ------------------------------------------
##
#' @param exp3_rows            rbind of fn_paper1_exp3_rows_from_ps5_csv() / _ps6_csv() / _ps7_directory() outputs.
#' @param target_min_ESS_by_N  named numeric vector, names = N as character; same target used by every software at a
#'                              given N (ps7's own target_min_ESS is not reused here - target_min_ESS_by_N is the single
#'                              source of truth for all three softwares, so the comparison is like-for-like).
#' @return tibble, one row per (device, software_label, N): target_min_ESS, n_runs_pooled, mean_time_to_target_seconds, median_time_to_target_seconds,
#'         min_time_to_target_seconds, max_time_to_target_seconds.
fn_paper1_exp3_time_to_target <- function( exp3_rows,
                                           target_min_ESS_by_N) {

        if (nrow(exp3_rows) == 0) stop("fn_paper1_exp3_time_to_target: exp3_rows has no rows.")
        if (is.null(names(target_min_ESS_by_N)) || any(!nzchar(names(target_min_ESS_by_N)))) {
              stop("fn_paper1_exp3_time_to_target: target_min_ESS_by_N must be a numeric vector named by N.")
        }
        ##
        N_without_target <- setdiff(unique(exp3_rows$N), as.numeric(names(target_min_ESS_by_N)))
        if (length(N_without_target) > 0) {
              stop("fn_paper1_exp3_time_to_target: no target_min_ESS for N = ", paste(N_without_target, collapse = ", "))
        }
        ##
        exp3_rows$target_min_ESS <- unname(target_min_ESS_by_N[as.character(exp3_rows$N)])
        ##
        if (any(exp3_rows$target_min_ESS <= 0) || any(exp3_rows$min_ESS_over_parameters_of_interest <= 0)) {
              stop("fn_paper1_exp3_time_to_target: target_min_ESS and min_ESS_over_parameters_of_interest must both be positive.")
        }
        ##
        ## ---- Posterior-summary time is scaled together with the sampling time (equation eq:paper1_time_to_target_ess); sources
        ##      without a recorded summary time (ps5 Mplus, ps6 Stan) contribute zero here, so their totals are lower bounds.
        if (!("time_summary_seconds" %in% names(exp3_rows))) exp3_rows$time_summary_seconds <- 0
        exp3_rows$time_summary_seconds[is.na(exp3_rows$time_summary_seconds)] <- 0
        ##
        exp3_rows$time_to_target_seconds <- exp3_rows$time_burnin_seconds +
            (exp3_rows$time_sampling_seconds + exp3_rows$time_summary_seconds) *
            (exp3_rows$target_min_ESS / exp3_rows$min_ESS_over_parameters_of_interest)
        ##
        ## ---- Mplus discards the first half of each chain as burn-in, so its burn-in grows with the run length: scale the whole run.
        is_mplus_row <- exp3_rows$software_label == unname(fn_paper1_exp3_software_labels()["mplus"])
        exp3_rows$time_to_target_seconds[is_mplus_row] <- (exp3_rows$time_burnin_seconds[is_mplus_row] + exp3_rows$time_sampling_seconds[is_mplus_row] +
                                                            exp3_rows$time_summary_seconds[is_mplus_row]) *
            (exp3_rows$target_min_ESS[is_mplus_row] / exp3_rows$min_ESS_over_parameters_of_interest[is_mplus_row])
        ##
        if (any(!is.finite(exp3_rows$time_to_target_seconds))) stop("fn_paper1_exp3_time_to_target: non-finite time_to_target_seconds produced.")
        ##
        exp3_summary <- exp3_rows |>
            dplyr::group_by(.data$device, .data$software_label, .data$N) |>
            dplyr::summarise( target_min_ESS                = dplyr::first(.data$target_min_ESS),
                              n_runs_pooled                  = dplyr::n(),
                              mean_time_to_target_seconds    = mean(.data$time_to_target_seconds),
                              median_time_to_target_seconds  = stats::median(.data$time_to_target_seconds),
                              min_time_to_target_seconds     = min(.data$time_to_target_seconds),
                              max_time_to_target_seconds     = max(.data$time_to_target_seconds),
                              .groups                        = "drop")
        ##
        return(exp3_summary)

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
## ---- format seconds: 3 significant figures, minutes in brackets when > 120 s ----------------------------------------------------------
##
fn_paper1_exp3_format_seconds <- function(seconds_value) {

        if (!is.finite(seconds_value)) return("---")
        ##
        seconds_text <- paste0(fn_paper1_format_number_commas_from_10000(signif(seconds_value, 3)), " s")
        ##
        if (seconds_value > 120) {
              minutes_text <- format(round(seconds_value / 60, 1), nsmall = 1, trim = TRUE)
              seconds_text <- paste0(seconds_text, " (", minutes_text, " min)")
        }
        ##
        return(seconds_text)

}


##
## ---- Experiment 3 LaTeX table for one device: rows BayesMVP / Stan (NUTS) / Mplus (PX-Gibbs), columns per N -------------------------
##
#' @param exp3_summary   the tibble returned by fn_paper1_exp3_time_to_target().
#' @param device_label   which device's rows of exp3_summary to use (one table per device).
#' @param N_values       N values to show, left to right, in this exact order.
#' @param table_caption  LaTeX caption text (no leading/trailing \\caption{}).
#' @param table_label    LaTeX \\label{} text (no leading/trailing braces).
#' @param output_file_path   where to write the .tex file.
#' @return invisibly, list(file_path, device, caption, label): fed straight into fn_paper1_exp3_write_bundle().
fn_paper1_exp3_table_tex <- function( exp3_summary,
                                      device_label,
                                      N_values,
                                      table_caption,
                                      table_label,
                                      output_file_path) {

        if (!isTRUE(is.character(device_label)) || length(device_label) != 1 || !nzchar(device_label)) {
              stop("fn_paper1_exp3_table_tex: device_label must be one non-empty string.")
        }
        if (length(N_values) == 0) stop("fn_paper1_exp3_table_tex: N_values is empty.")
        ##
        device_summary <- dplyr::filter(.data = exp3_summary, .data$device == device_label)
        ##
        software_labels_in_row_order <- unname(fn_paper1_exp3_software_labels()[c("bayesmvp", "stan", "mplus")])
        ##
        ## ---- one lookup: time (seconds) for every software x N cell, NA where nothing was found: --------------------------------------
        ##
        time_seconds_by_software_and_N <- sapply( X   = N_values,
                                                  FUN = function(one_N_value) {

              sapply( X   = software_labels_in_row_order,
                     FUN = function(one_software_label) {

                    matching_rows <- dplyr::filter( .data      = device_summary,
                                                    .data$N               == one_N_value,
                                                    .data$software_label  == one_software_label)
                    ##
                    if (nrow(matching_rows) == 0) return(NA_real_)
                    if (nrow(matching_rows) > 1)  stop("fn_paper1_exp3_table_tex: more than one summary row for device = ", device_label,
                                                       ", N = ", one_N_value, ", software = ", one_software_label)
                    return(matching_rows$mean_time_to_target_seconds)

              })

        })
        ##
        rownames(time_seconds_by_software_and_N) <- software_labels_in_row_order
        colnames(time_seconds_by_software_and_N) <- as.character(N_values)
        ##
        bayesmvp_row_name <- software_labels_in_row_order[1]
        ##
        ## ---- build the table body, bolding the fastest software in each N column: ---------------------------------------------------
        ##
        table_body_lines <- character(0)
        ##
        for (one_software_label in software_labels_in_row_order) {

              time_cells_text    <- character(length(N_values))
              speedup_cells_text <- character(length(N_values))
              ##
              for (N_index in seq_along(N_values)) {

                    one_N_value             <- N_values[N_index]
                    this_cell_time_seconds  <- time_seconds_by_software_and_N[one_software_label, N_index]
                    fastest_time_in_column  <- suppressWarnings(min(time_seconds_by_software_and_N[, N_index], na.rm = TRUE))
                    ##
                    formatted_time <- fn_paper1_exp3_format_seconds(this_cell_time_seconds)
                    if (is.finite(this_cell_time_seconds) && is.finite(fastest_time_in_column) &&
                        isTRUE(all.equal(this_cell_time_seconds, fastest_time_in_column))) {
                          formatted_time <- paste0("\\textbf{", formatted_time, "}")
                    }
                    time_cells_text[N_index] <- formatted_time
                    ##
                    bayesmvp_time_seconds <- time_seconds_by_software_and_N[bayesmvp_row_name, N_index]
                    ##
                    speedup_cells_text[N_index] <- if (identical(one_software_label, bayesmvp_row_name)) {
                          if (is.finite(this_cell_time_seconds)) "1$\\times$" else "---"
                    } else if (is.finite(this_cell_time_seconds) && is.finite(bayesmvp_time_seconds)) {
                          paste0(format(signif(this_cell_time_seconds / bayesmvp_time_seconds, 3), trim = TRUE), "$\\times$")
                    } else "---"

              }
              ##
              table_body_lines <- c(table_body_lines,
                                    paste0(paste(c(one_software_label, as.vector(rbind(time_cells_text, speedup_cells_text))), collapse = " & "), " \\\\"))

        }
        ##
        column_spec  <- paste0("l", strrep("rr", length(N_values)))
        header_N_row <- paste0(" & ", paste(paste0("\\multicolumn{2}{c}{$N = ", fn_paper1_format_number_commas_from_10000(N_values), "$}"), collapse = " & "), " \\\\")
        header_2_row <- paste0(" & ", paste(rep("Time & Speed-up", length(N_values)), collapse = " & "), " \\\\")
        ##
        table_tex_lines <- c( "\\begin{table}[H]",
                              "\\centering",
                              paste0("\\caption{", table_caption, "}"),
                              paste0("\\label{", table_label, "}"),
                              paste0("\\begin{tabular}{", column_spec, "}"),
                              "\\toprule",
                              header_N_row,
                              header_2_row,
                              "\\midrule",
                              table_body_lines,
                              "\\bottomrule",
                              "\\end{tabular}",
                              "\\end{table}")
        ##
        dir.create(dirname(output_file_path), recursive = TRUE, showWarnings = FALSE)
        writeLines(table_tex_lines, output_file_path)
        ##
        message(paste0("\033[36m", "Experiment 3 table for device = ", device_label, " written to: ", output_file_path, "\033[0m"))
        ##
        return(invisible(list( file_path = output_file_path,
                               device    = device_label,
                               caption   = table_caption,
                               label     = table_label)))

}


##
## ---- assemble one or more device tables into the Exp3_absolute_efficiency bundle (paper_sections/<bundle>/...) ----------------------
##
#' @param table_results   list of return values from fn_paper1_exp3_table_tex(), one per device, in the order they should
#'                         appear in section.tex.
#' @param output_dir      the R-side output directory (bundle written to <output_dir>/paper_sections/Exp3_absolute_efficiency).
#' @param manuscript_dir   if not NULL, also copied to <manuscript_dir>/Files/Generated/Exp3_absolute_efficiency, exactly as
#'                          fn_paper1_export_manuscript() does for the other bundles in R_fns_alg_paper_1_figures_tables.R.
#' @return invisibly, the bundle directory path.
fn_paper1_exp3_write_bundle <- function( table_results,
                                         output_dir,
                                         manuscript_dir) {

        if (length(table_results) == 0) stop("fn_paper1_exp3_write_bundle: table_results is empty.")
        ##
        bundle_dir <- file.path(output_dir, "paper_sections", "Exp3_absolute_efficiency")
        dir.create(file.path(bundle_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
        ##
        section_tex_lines <- c("% Generated by fn_paper1_exp3_write_bundle(); see R_fns_alg_paper_1_experiment_3.R.")
        ##
        for (one_table_result in table_results) {

              destination_file_name <- basename(one_table_result$file_path)
              destination_file_path <- file.path(bundle_dir, "tables", destination_file_name)
              ##
              ## fn_paper1_exp3_table_tex() may already have written straight into the bundle's own tables/ directory
              ## (the runner script does this); only copy when the source is genuinely a different file.
              ##
              source_is_destination <- file.exists(one_table_result$file_path) && file.exists(destination_file_path) &&
                  identical(normalizePath(one_table_result$file_path), normalizePath(destination_file_path))
              if (!source_is_destination && !file.copy(one_table_result$file_path, destination_file_path, overwrite = TRUE)) {
                    stop("fn_paper1_exp3_write_bundle: could not copy ", one_table_result$file_path, " into the bundle.")
              }
              ##
              manuscript_path <- paste0("Files/Generated/Exp3_absolute_efficiency/tables/", destination_file_name)
              section_tex_lines <- c(section_tex_lines, paste0("\\input{", manuscript_path, "}"))

        }
        ##
        writeLines(section_tex_lines, file.path(bundle_dir, "section.tex"))
        ##
        message(paste0("\033[36m", "Experiment 3 bundle written to: ", bundle_dir, "\033[0m"))
        ##
        if (!is.null(manuscript_dir)) {

              destination_dir <- file.path(manuscript_dir, "Files", "Generated", "Exp3_absolute_efficiency")
              dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
              ##
              files_no_longer_generated <- setdiff(list.files(destination_dir, recursive = TRUE),
                                                    list.files(bundle_dir, recursive = TRUE))
              ##
              if (length(files_no_longer_generated) > 0) {
                    message(paste0("\033[36m", "Files in ", destination_dir,
                                   " that this export does not write (stale; not \\input): ",
                                   paste(files_no_longer_generated, collapse = ", "), "\033[0m"))
              }
              ##
              for (one_bundle_file in list.files(bundle_dir, recursive = TRUE)) {

                    dir.create(dirname(file.path(destination_dir, one_bundle_file)), recursive = TRUE, showWarnings = FALSE)
                    if (!file.copy(file.path(bundle_dir, one_bundle_file), file.path(destination_dir, one_bundle_file), overwrite = TRUE)) {
                          stop("fn_paper1_exp3_write_bundle: could not copy ", one_bundle_file, " into ", destination_dir)
                    }

              }
              ##
              message(paste0("\033[36m", "Experiment 3 bundle copied to manuscript directory: ", destination_dir, "\033[0m"))

        }
        ##
        return(invisible(bundle_dir))

}






















