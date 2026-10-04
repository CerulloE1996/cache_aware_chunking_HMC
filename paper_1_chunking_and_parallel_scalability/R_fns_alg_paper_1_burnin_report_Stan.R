##
## =======================================================================================================================================
## R_fns_alg_paper_1_burnin_report_Stan.R
##
## Paper 1 burn-in report, extended with (i) the tape-chunked Stan burn-in study (local-HPC) and (ii) any additional NicoStan+BayesMVP
## burn-in files (e.g. the N = 50,000, N_chunks = 1000 runs). The helpers in R_fns_alg_paper_1_burnin_report.R are sourced by the
## runner (alg_paper_1_burnin_report_Stan.R) and reused unchanged; this file adds:
##
##   - fn_paper1_burnin_read_runner_settings():        the CURRENT settings of a burn-in runner (chains, chunks, WCP, L, timed
##                                                      iterations, runs, thread limit), read from the runner file at run time;
##   - fn_paper1_burnin_read_Stan_runner_settings():   the same for ps_1_burnin_Stan_tape_chunked.R (NicoStan+BayesMVP burn-in settings
##                                                      with the sampling runner's stan_chunk_candidates, as that runner reads them);
##   - fn_paper1_burnin_filter_time_completed():       optional: only Stan rows completed after a given time (e.g. after a rebuild);
##   - fn_paper1_burnin_expected_configurations():      every (N, chains, chunks, WCP) configuration those settings define;
##   - fn_paper1_burnin_list_files() / _read_files():   an explicit list of saved .rds files per implementation and device;
##   - fn_paper1_burnin_stop_on_duplicate_runs():       stops if two files hold the same configuration and run number;
##   - fn_paper1_burnin_filter_to_runner_settings():    keeps only runs matching the runner settings, and lists missing configurations;
##   - fn_paper1_burnin_best_by_cell():                 per (N, burn-in chains): best configuration and both speed-ups;
##   - fn_paper1_burnin_figure_Stan():                  the tape-chunked Stan figure, one for all devices (same layout as the
##                                                      NicoStan+BayesMVP figure, plus the best NicoStan+BayesMVP level of each panel);
##   - fn_paper1_burnin_placeholder_values():           the values of every \BurninTBD{ID} placeholder in the draft section;
##   - fn_paper1_burnin_checks():                       checks of the existing burn-in sentences against the current results;
##   - fn_paper1_burnin_fill_placeholders():            a filled copy of the draft section (the draft itself is not changed);
##   - fn_paper1_burnin_write_Stan_bundle():            writes the figures, CSV files, the research-output table and the filled draft.
##
## This file only defines functions. It never runs sampling, compiles anything, or edits the benchmark scripts or their saved .rds files.
##
## Column schema of the saved .rds files: as documented at the top of R_fns_alg_paper_1_burnin_report.R.
##


##
## ---- Small local helpers: ----------------------------------------------------------------------------------------------------------
##
fn_paper1_burnin_message <-  function( text ) {

        ##
        ## ---- cyan informational progress via NicoStan::colourise(), with a plain-text fallback:
        ##
        if (requireNamespace("NicoStan", quietly = TRUE)) {
              message(NicoStan::colourise(text = text, fg = "cyan"))
        } else {
              message(text)
        }

}
##
fn_paper1_burnin_chains_word <-  function( n_chains_burnin ) {

        chains_words <-  c("4" = "four", "8" = "eight", "16" = "sixteen", "2" = "two", "32" = "thirty-two")
        chains_key   <-  as.character(x = n_chains_burnin)
        if (!chains_key %in% names(x = chains_words)) return(chains_key)
        return(unname(obj = chains_words[chains_key]))

}
##
fn_paper1_burnin_N_key <-  function( N ) {

        return(format(x = N, scientific = FALSE, trim = TRUE))

}
##
fn_paper1_burnin_replace_fixed <-  function( text,
                                             pattern,
                                             replacement
) {

        ##
        ## ---- literal (not regular-expression) replacement of every occurrence, so that "$" and "\" in LaTeX values are kept as written:
        ##
        output_pieces  <-  character(0)
        remaining_text <-  text
        repeat {

              match_position <-  regexpr(pattern = pattern, text = remaining_text, fixed = TRUE)
              if (match_position == -1) break
              output_pieces  <-  c(output_pieces, substr(x = remaining_text, start = 1, stop = match_position - 1), replacement)
              remaining_text <-  substr(x = remaining_text, start = match_position + nchar(x = pattern), stop = nchar(x = remaining_text))

        }
        return(paste0(paste(output_pieces, collapse = ""), remaining_text))

}


##
## ---- fn_paper1_burnin_read_runner_settings: the CURRENT burnin_benchmark_settings of a burn-in runner, for one device: -----------------
##
## The runner's settings block (the top-level { ... } expression containing "burnin_benchmark_settings <- list()") is evaluated on its
## own, in a separate environment with `computer` set to "Local_HPC" or "Laptop"; nothing else in the runner is run (no setwd(), no
## packages, no data, no benchmark).
##
fn_paper1_burnin_read_runner_settings <-  function( runner_file_path,
                                                    device
) {

        if (!is.character(x = runner_file_path) || length(x = runner_file_path) != 1 || is.na(x = runner_file_path) ||
            !file.exists(runner_file_path)) {
              stop("fn_paper1_burnin_read_runner_settings: runner file not found: ", paste(runner_file_path, collapse = ", "))
        }
        ##
        runner_expressions <-  parse(file = runner_file_path, keep.source = FALSE)
        is_settings_block  <-  vapply( X         = seq_along(along.with = runner_expressions),
                                       FUN       = function(expression_index) {
                                             expression_text <-  deparse(expr = runner_expressions[[expression_index]], width.cutoff = 500)
                                             any(grepl(pattern = "burnin_benchmark_settings <- list()", x = expression_text, fixed = TRUE))
                                       },
                                       FUN.VALUE = logical(1))
        if (sum(is_settings_block) != 1) {
              stop("fn_paper1_burnin_read_runner_settings: expected exactly one top-level settings block containing ",
                   "'burnin_benchmark_settings <- list()' in ", runner_file_path, ", found ", sum(is_settings_block), ".")
        }
        ##
        settings_environment <-  new.env(parent = globalenv())
        assign(x = "computer", value = fn_paper1_burnin_device_field(device_argument = device), envir = settings_environment)
        tryCatch( expr  = eval(expr = runner_expressions[[which(is_settings_block)]], envir = settings_environment),
                  error = function(error) {
                        stop("fn_paper1_burnin_read_runner_settings: could not evaluate the settings block of ", runner_file_path,
                             " for device = '", device, "': ", conditionMessage(error))
                  })
        runner_settings <-  get(x = "burnin_benchmark_settings", envir = settings_environment)
        ##
        required_setting_names <-  c("n_chains_burnin_vec", "n_threads_WCP_burnin_vec_given_N", "num_chunks_burnin_vec_given_N",
                                     "n_timed_iters_given_N", "L_main", "n_runs", "n_total_threads")
        missing_setting_names  <-  setdiff(x = required_setting_names, y = names(x = runner_settings))
        if (length(x = missing_setting_names) > 0) {
              stop("fn_paper1_burnin_read_runner_settings: ", runner_file_path, " (device = '", device, "') does not set: ",
                   paste(missing_setting_names, collapse = ", "))
        }
        runner_settings$runner_file_path <-  runner_file_path
        ##
        fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_read_runner_settings: ", basename(path = runner_file_path), " (", device,
                                               "): chains = {", paste(runner_settings$n_chains_burnin_vec, collapse = ", "),
                                               "}, L = ", runner_settings$L_main, ", runs = ", runner_settings$n_runs,
                                               ", thread limit = ", runner_settings$n_total_threads))
        return(runner_settings)

}


##
## ---- fn_paper1_burnin_read_last_literal_assignment: the last "name <- value" of a runner whose value is NULL, TRUE/FALSE or a list/c: -
##
## Searches the whole parsed runner (inside braces and if/else too) for assignments to object_name (a name, or e.g.
## "paper1_settings$include_BayesMVP_WCP_only"), keeps only those whose right-hand side is NULL, TRUE/FALSE or a list(...) / c(...)
## literal, and evaluates the last of these in an environment that only sees base R. Returns NULL if none exists.
##
fn_paper1_burnin_read_last_literal_assignment <-  function( runner_file_path,
                                                            object_name
) {

        if (!file.exists(runner_file_path)) stop("fn_paper1_burnin_read_last_literal_assignment: runner file not found: ", runner_file_path)
        parsed_runner <-  parse(file = runner_file_path, keep.source = FALSE)
        found_values  <-  list()
        ##
        collect_assignments <-  function(expression_to_search) {

              if (is.call(x = expression_to_search)) {
                    call_name <-  as.character(x = expression_to_search[[1]])[1]
                    if (call_name %in% c("<-", "=") && length(x = expression_to_search) == 3 &&
                        identical(x = paste(deparse(expr = expression_to_search[[2]]), collapse = ""), y = object_name)) {
                          value_expression <-  expression_to_search[[3]]
                          is_literal <-  is.null(x = value_expression) || (is.logical(x = value_expression) && length(x = value_expression) == 1) ||
                                         (is.call(x = value_expression) && as.character(x = value_expression[[1]])[1] %in% c("list", "c"))
                          if (is_literal) found_values[[length(x = found_values) + 1]] <<-  list(value_expression = value_expression)
                    }
                    for (argument_index in seq_along(along.with = expression_to_search)[-1]) {
                          argument_expression <-  expression_to_search[[argument_index]]
                          if (!missing(argument_expression)) collect_assignments(argument_expression)
                    }
              }
              invisible(NULL)

        }
        for (top_level_expression in as.list(x = parsed_runner)) collect_assignments(top_level_expression)
        ##
        if (length(x = found_values) == 0) return(NULL)
        return(eval(expr = found_values[[length(x = found_values)]]$value_expression, envir = new.env(parent = baseenv())))

}


##
## ---- fn_paper1_burnin_read_Stan_runner_settings: the settings of ps_1_burnin_Stan_tape_chunked.R, for one device: ---------------------
##
## That runner takes its settings from two other runners when it runs, so the same is done here (nothing is run):
##   - chains, WCP grid, L, eps, timed iterations, runs, thread limit: the NicoStan+BayesMVP burn-in runner's settings block;
##   - chunk grid: the last stan_chunk_candidates assignment of the sampling runner (alg_paper_1_chunking_WCP_par_scaling.R);
##   - timed iterations: replaced where the Stan runner sets stan_n_timed_iters_given_N to a list, or where the environment variable
##     PS1_BURNIN_STAN_N_TIMED_ITERS (e.g. "500=200,2500=40,10000=5,50000=3") is set - as in the Stan runner.
##
fn_paper1_burnin_read_Stan_runner_settings <-  function( BayesMVP_runner_file_path,
                                                         sampling_runner_file_path,
                                                         Stan_runner_file_path,
                                                         device
) {

        if (!file.exists(Stan_runner_file_path)) stop("fn_paper1_burnin_read_Stan_runner_settings: Stan runner not found: ", Stan_runner_file_path)
        ##
        runner_settings <-  fn_paper1_burnin_read_runner_settings( runner_file_path = BayesMVP_runner_file_path,
                                                                   device           = device)
        ##
        stan_chunk_candidates <-  fn_paper1_burnin_read_last_literal_assignment( runner_file_path = sampling_runner_file_path,
                                                                                 object_name      = "stan_chunk_candidates")
        if (is.null(x = stan_chunk_candidates)) stop("fn_paper1_burnin_read_Stan_runner_settings: no stan_chunk_candidates in ", sampling_runner_file_path)
        runner_settings$num_chunks_burnin_vec_given_N <-  stan_chunk_candidates
        ##
        stan_n_timed_iters_given_N <-  fn_paper1_burnin_read_last_literal_assignment( runner_file_path = Stan_runner_file_path,
                                                                                      object_name      = "stan_n_timed_iters_given_N")
        timed_iters_environment_variable <-  Sys.getenv("PS1_BURNIN_STAN_N_TIMED_ITERS")
        if (nzchar(x = timed_iters_environment_variable)) {
              timed_iters_pairs <-  strsplit(x = strsplit(x = timed_iters_environment_variable, split = ",", fixed = TRUE)[[1]], split = "=", fixed = TRUE)
              stan_n_timed_iters_given_N <-  stats::setNames( object = lapply(X = timed_iters_pairs, FUN = function(pair) as.numeric(x = pair[2])),
                                                              nm     = vapply(X = timed_iters_pairs, FUN = function(pair) pair[1], FUN.VALUE = character(1)))
        }
        if (!is.null(x = stan_n_timed_iters_given_N)) {
              runner_settings$n_timed_iters_given_N <-  utils::modifyList(x = runner_settings$n_timed_iters_given_N, val = stan_n_timed_iters_given_N)
        }
        ##
        runner_settings$runner_file_path <-  Stan_runner_file_path
        fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_read_Stan_runner_settings (", device, "): Stan chunks ",
                                               paste0(names(x = stan_chunk_candidates), ": {", vapply(X = stan_chunk_candidates, FUN = paste, FUN.VALUE = character(1), collapse = ", "), "}", collapse = "; "),
                                               "; timed iterations ", paste0(names(x = runner_settings$n_timed_iters_given_N), " = ",
                                                                             unlist(x = runner_settings$n_timed_iters_given_N), collapse = ", ")))
        return(runner_settings)

}


##
## ---- fn_paper1_burnin_filter_time_completed: keep only rows completed at or after a given time (e.g. after a rebuild): -----------------
##
## Rows without a time_completed column (all NicoStan+BayesMVP files) are returned unchanged. minimum_time_completed = NULL keeps all rows.
##
fn_paper1_burnin_filter_time_completed <-  function( rows,
                                                     minimum_time_completed,
                                                     implementation,
                                                     device
) {

        if (is.null(x = minimum_time_completed) || nrow(x = rows) == 0) return(rows)
        if (!"time_completed" %in% names(x = rows)) {
              stop("fn_paper1_burnin_filter_time_completed: ", implementation, " (", device, ") rows have no time_completed column, ",
                   "so the minimum completion time cannot be applied.")
        }
        completion_times <-  as.POSIXct(x = rows$time_completed, format = "%Y-%m-%d %H:%M:%S")
        is_recent        <-  !is.na(x = completion_times) & completion_times >= as.POSIXct(x = minimum_time_completed, format = "%Y-%m-%d %H:%M:%S")
        fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_filter_time_completed: ", implementation, " (", device, "): kept ", sum(is_recent),
                                               " of ", nrow(x = rows), " row(s) completed at or after ", minimum_time_completed))
        return(rows[is_recent, , drop = FALSE])

}


##
## ---- fn_paper1_burnin_expected_configurations: every configuration the runner settings define (as the benchmark builds its grid): ----
##
## Invalid configurations are skipped exactly as in R_fn_ps1_burnin_make_configuration_grid(): N_WCP > N_chunks (a thread with no
## chunk), or N_chains x N_WCP above the runner's thread limit. With include_WCP_only = TRUE, the WCP-only configurations of the sampling
## runner are added as in R_fn_ps1_burnin_add_configuration_sets(): N_chunks = N_WCP for every N_WCP > 1 not already in the chunk grid.
##
fn_paper1_burnin_expected_configurations <-  function( runner_settings,
                                                       N_values,
                                                       include_WCP_only
) {

        expected_rows <-  list()
        ##
        for (N in N_values) {

              N_key         <-  fn_paper1_burnin_N_key(N = N)
              chunk_values  <-  runner_settings$num_chunks_burnin_vec_given_N[[N_key]]
              WCP_by_chains <-  runner_settings$n_threads_WCP_burnin_vec_given_N[[N_key]]
              n_timed_iters <-  runner_settings$n_timed_iters_given_N[[N_key]]
              ##
              if (is.null(x = chunk_values) || is.null(x = WCP_by_chains) || is.null(x = n_timed_iters)) {
                    fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_expected_configurations: N = ", N, " is not in the grid of ",
                                                           basename(path = runner_settings$runner_file_path), "; no configurations expected."))
                    next
              }
              ##
              for (n_chains_burnin in runner_settings$n_chains_burnin_vec) {

                    WCP_values <-  WCP_by_chains[[as.character(x = n_chains_burnin)]]
                    if (is.null(x = WCP_values)) {
                          stop("fn_paper1_burnin_expected_configurations: no WCP entry for N = ", N, ", ", n_chains_burnin, " burn-in chains in ",
                               runner_settings$runner_file_path)
                    }
                    one_grid <-  tidyr::expand_grid( num_chunks_burnin    = chunk_values,
                                                     n_threads_WCP_burnin = WCP_values)
                    is_valid <-  (one_grid$n_threads_WCP_burnin <= one_grid$num_chunks_burnin) &
                                 (n_chains_burnin * one_grid$n_threads_WCP_burnin <= runner_settings$n_total_threads)
                    one_grid <-  one_grid[is_valid, , drop = FALSE]
                    one_grid$configuration_set <-  rep(x = "chunk_grid", times = nrow(x = one_grid))
                    ##
                    if (isTRUE(x = include_WCP_only)) {
                          WCP_only_values <-  WCP_values[WCP_values > 1 & !(WCP_values %in% chunk_values) &
                                                         n_chains_burnin * WCP_values <= runner_settings$n_total_threads]
                          one_grid <-  dplyr::bind_rows(one_grid, tibble::tibble( num_chunks_burnin    = WCP_only_values,
                                                                                  n_threads_WCP_burnin = WCP_only_values,
                                                                                  configuration_set    = rep(x = "WCP_only", times = length(x = WCP_only_values))))
                    }
                    ##
                    expected_rows[[length(x = expected_rows) + 1]] <-  tibble::tibble( N                    = rep(x = N, times = nrow(x = one_grid)),
                                                                                       n_chains_burnin      = rep(x = n_chains_burnin, times = nrow(x = one_grid)),
                                                                                       num_chunks_burnin    = one_grid$num_chunks_burnin,
                                                                                       n_threads_WCP_burnin = one_grid$n_threads_WCP_burnin,
                                                                                       n_timed_iters        = rep(x = n_timed_iters, times = nrow(x = one_grid)),
                                                                                       L_main               = rep(x = runner_settings$L_main, times = nrow(x = one_grid)),
                                                                                       configuration_set    = one_grid$configuration_set)

              }

        }
        ##
        return(dplyr::bind_rows(expected_rows))

}


##
## ---- fn_paper1_burnin_list_files: the saved .rds files matching a regular expression, plus explicitly listed extra files: -------------
##
fn_paper1_burnin_list_files <-  function( output_dir,
                                          file_regex,
                                          extra_file_paths
) {

        if (!dir.exists(paths = output_dir)) {
              stop("fn_paper1_burnin_list_files: output_dir does not exist: ", output_dir)
        }
        matched_file_paths <-  list.files(path = output_dir, pattern = file_regex, full.names = TRUE)
        ##
        if (length(x = extra_file_paths) > 0) {
              missing_extra_files <-  extra_file_paths[!file.exists(extra_file_paths)]
              if (length(x = missing_extra_files) > 0) {
                    stop("fn_paper1_burnin_list_files: extra file(s) not found: ", paste(missing_extra_files, collapse = ", "))
              }
        }
        all_file_paths <-  unique(x = normalizePath(path = c(matched_file_paths, extra_file_paths), mustWork = TRUE))
        fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_list_files: ", length(x = all_file_paths), " file(s) for '", file_regex, "'",
                                               if (length(x = all_file_paths) > 0) paste0(": ", paste(basename(path = all_file_paths), collapse = ", ")) else ""))
        return(all_file_paths)

}


##
## ---- fn_paper1_burnin_read_files: read an explicit list of saved burn-in .rds files for ONE device and ONE implementation: ------------
##
## Same checks as fn_paper1_burnin_read_outputs() (older column names translated in memory, device field checked, incomplete rows
## skipped with a message, divergent rows kept and flagged), with an "implementation" column added.
##
fn_paper1_burnin_read_files <-  function( file_paths,
                                          device,
                                          implementation
) {

        expected_device_field <-  fn_paper1_burnin_device_field(device_argument = device)
        required_columns      <-  c("device", "N", "n_chains_burnin", "num_chunks_burnin", "n_threads_WCP_burnin",
                                    "n_threads_total", "run_number", "L_main", "n_timed_iters", "total_timed_seconds",
                                    "mean_sec_per_iter", "n_divs")
        previous_to_current_column_names <-  c( n_timed_iters       = "n_timed_iterations",
                                                mean_sec_per_iter   = "mean_seconds_per_iteration",
                                                median_sec_per_iter = "median_seconds_per_iteration",
                                                n_divs              = "n_divergent_transitions")
        ##
        per_file_rows <-  vector(mode = "list", length = length(x = file_paths))
        ##
        for (file_index in seq_along(along.with = file_paths)) {

              this_file_path <-  file_paths[file_index]
              this_result    <-  readRDS(file = this_file_path)
              if (!is.data.frame(x = this_result)) {
                    stop("fn_paper1_burnin_read_files: ", this_file_path, " does not contain a data frame.")
              }
              ##
              for (current_column_name in names(x = previous_to_current_column_names)) {

                    previous_column_name <-  previous_to_current_column_names[[current_column_name]]
                    if (!current_column_name %in% names(x = this_result) && previous_column_name %in% names(x = this_result)) {
                          this_result[[current_column_name]] <-  this_result[[previous_column_name]]
                    }

              }
              if (!all(required_columns %in% names(x = this_result))) {
                    stop("fn_paper1_burnin_read_files: ", this_file_path, " is missing required column(s): ",
                         paste(setdiff(x = required_columns, y = names(x = this_result)), collapse = ", "))
              }
              if (!"median_sec_per_iter" %in% names(x = this_result)) {
                    this_result$median_sec_per_iter <-  NA_real_
              }
              ##
              observed_device_values <-  unique(x = this_result$device)
              if (!all(observed_device_values == expected_device_field)) {
                    stop("fn_paper1_burnin_read_files: ", this_file_path, " was listed for device = '", device, "', but its 'device' column contains ",
                         paste(observed_device_values, collapse = ", "), " (expected only '", expected_device_field, "').")
              }
              ##
              usable_rows <-  complete.cases(this_result[required_columns]) &
                              rowSums(!is.finite(as.matrix(this_result[setdiff(x = required_columns, y = "device")]))) == 0 &
                              this_result$n_timed_iters > 0 & this_result$total_timed_seconds > 0 & this_result$mean_sec_per_iter > 0
              usable_rows[is.na(x = usable_rows)] <-  FALSE
              if (any(!usable_rows)) {
                    fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_read_files: skipping ", sum(!usable_rows),
                                                           " incomplete/invalid row(s) in ", basename(path = this_file_path)))
              }
              this_result <-  this_result[usable_rows, , drop = FALSE]
              ##
              this_result$converged <-  this_result$n_divs == 0
              if (any(!this_result$converged)) {
                    fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_read_files: ", sum(!this_result$converged),
                                                           " row(s) with divergent transitions in ", basename(path = this_file_path), " (kept, converged = FALSE)"))
              }
              this_result$source_file    <-  this_file_path
              this_result$implementation <-  implementation
              ##
              ## ---- optional columns kept when present (Stan files: completion time and configuration set):
              ##
              optional_columns <-  intersect(x = c("time_completed", "configuration_set", "Stan_model_loading", "sec_per_iter_each"), y = names(x = this_result))
              per_file_rows[[file_index]] <-  tibble::as_tibble(x = this_result[, c(required_columns, "median_sec_per_iter", "converged",
                                                                                     "source_file", "implementation", optional_columns)])

        }
        ##
        tidy_rows <-  dplyr::bind_rows(per_file_rows)
        fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_read_files: ", implementation, " (", device, "): read ", nrow(x = tidy_rows),
                                               " usable row(s) from ", length(x = file_paths), " file(s)"))
        return(tidy_rows)

}


##
## ---- fn_paper1_burnin_stop_on_duplicate_runs: two files must never hold the same configuration and run number: -------------------------
##
fn_paper1_burnin_stop_on_duplicate_runs <-  function( rows,
                                                      implementation
) {

        if (nrow(x = rows) == 0) return(invisible(TRUE))
        key_columns    <-  c("device", "N", "n_chains_burnin", "num_chunks_burnin", "n_threads_WCP_burnin", "run_number")
        duplicate_keys <-  rows %>%
              dplyr::group_by(dplyr::across(.cols = dplyr::all_of(x = key_columns))) %>%
              dplyr::summarise( n_rows       = dplyr::n(),
                                source_files = paste(unique(x = basename(path = .data$source_file)), collapse = " + "),
                                .groups      = "drop") %>%
              dplyr::filter(.data$n_rows > 1)
        if (nrow(x = duplicate_keys) > 0) {
              print(x = duplicate_keys, n = 50, width = Inf)
              stop("fn_paper1_burnin_stop_on_duplicate_runs: ", implementation, ": ", nrow(x = duplicate_keys),
                   " (configuration, run) key(s) appear more than once across the saved files (printed above). ",
                   "Decide which measurement to keep (e.g. move the older file out of the file list) before reporting.")
        }
        return(invisible(TRUE))

}


##
## ---- fn_paper1_burnin_filter_to_runner_settings: keep only runs matching the runner's current settings; list missing configurations: --
##
fn_paper1_burnin_filter_to_runner_settings <-  function( rows,
                                                         expected_configurations,
                                                         n_runs,
                                                         device,
                                                         implementation
) {

        key_columns  <-  c("N", "n_chains_burnin", "num_chunks_burnin", "n_threads_WCP_burnin", "n_timed_iters", "L_main")
        ##
        matched_rows <-  dplyr::semi_join(x = rows, y = expected_configurations, by = key_columns)
        matched_rows <-  matched_rows[matched_rows$run_number <= n_runs, , drop = FALSE]
        ##
        excluded_rows <-  dplyr::anti_join(x = rows, y = matched_rows, by = c(key_columns, "run_number"))
        if (nrow(x = excluded_rows) > 0) {
              excluded_configurations <-  unique(x = excluded_rows[, c(key_columns, "run_number"), drop = FALSE])
              fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_filter_to_runner_settings: ", implementation, " (", device, "): excluded ",
                                                     nrow(x = excluded_rows), " row(s) whose settings are not in the current runner grid (kept on disk, not reported):"))
              print(x = excluded_configurations, n = 50, width = Inf)
        }
        ##
        observed_runs <-  matched_rows %>%
              dplyr::group_by(dplyr::across(.cols = dplyr::all_of(x = key_columns))) %>%
              dplyr::summarise(n_runs_observed = dplyr::n_distinct(.data$run_number), .groups = "drop")
        coverage <-  dplyr::left_join(x = expected_configurations, y = observed_runs, by = key_columns)
        coverage$n_runs_observed[is.na(x = coverage$n_runs_observed)] <-  0
        ##
        missing_configurations <-  coverage[coverage$n_runs_observed < n_runs, , drop = FALSE]
        missing_configurations$device         <-  rep(x = device, times = nrow(x = missing_configurations))
        missing_configurations$implementation <-  rep(x = implementation, times = nrow(x = missing_configurations))
        missing_configurations$n_runs_expected <-  rep(x = n_runs, times = nrow(x = missing_configurations))
        if (nrow(x = missing_configurations) > 0) {
              fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_filter_to_runner_settings: ", implementation, " (", device, "): ",
                                                     nrow(x = missing_configurations), " configuration(s) in the runner grid have fewer than ", n_runs,
                                                     " saved run(s):"))
              print(x = missing_configurations, n = 50, width = Inf)
        }
        fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_filter_to_runner_settings: ", implementation, " (", device, "): kept ",
                                               nrow(x = matched_rows), " row(s)"))
        ##
        return(list(rows = matched_rows, missing_configurations = missing_configurations))

}


##
## ---- fn_paper1_burnin_keep_Stan_model_loading: keep only Stan rows of the given NicoStan build (column Stan_model_loading): --------------
##
fn_paper1_burnin_keep_Stan_model_loading <-  function( rows,
                                                       required_model_loading,
                                                       device
) {

        if (nrow(x = rows) == 0) return(rows)
        if (!"Stan_model_loading" %in% names(x = rows)) {
              stop("fn_paper1_burnin_keep_Stan_model_loading: Stan rows (", device, ") have no Stan_model_loading column.")
        }
        is_required_build <-  !is.na(x = rows$Stan_model_loading) & rows$Stan_model_loading == required_model_loading
        fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_keep_Stan_model_loading: Stan (", device, "): kept ", sum(is_required_build), " of ",
                                               nrow(x = rows), " row(s) with Stan_model_loading = '", required_model_loading, "'"))
        return(rows[is_required_build, , drop = FALSE])

}


##
## ---- fn_paper1_burnin_read_leapfrog_steps: leapfrog steps of every timed burn-in iteration (run_number, iteration, n_leapfrog_steps): ---
##
## The shared trajectory length depends only on the seed (1000 x run_number) and the iteration, so these counts hold for every N, chain
## count, chunk count, WCP count and implementation.
##
fn_paper1_burnin_read_leapfrog_steps <-  function( csv_path ) {

        if (!file.exists(csv_path)) stop("fn_paper1_burnin_read_leapfrog_steps: not found: ", csv_path)
        leapfrog_steps <-  utils::read.csv(file = csv_path, stringsAsFactors = FALSE)
        if (!all(c("run_number", "iteration", "n_leapfrog_steps") %in% names(x = leapfrog_steps))) {
              stop("fn_paper1_burnin_read_leapfrog_steps: ", csv_path, " needs columns run_number, iteration, n_leapfrog_steps.")
        }
        leapfrog_steps$n_leapfrog_steps <-  as.numeric(x = leapfrog_steps$n_leapfrog_steps)
        return(leapfrog_steps)

}


##
## ---- fn_paper1_burnin_add_sec_per_step: seconds per leapfrog step and mean leapfrog steps per timed iteration, per row: ------------------
##
## Timed seconds: sum of sec_per_iter_each where present (Stan files), otherwise total_timed_seconds (NicoStan+BayesMVP files); divided by
## the leapfrog steps of iterations 1 to n_timed_iters of that run. Needed because the Stan runs timed fewer iterations at N >= 10,000.
##
fn_paper1_burnin_add_sec_per_step <-  function( rows,
                                                leapfrog_steps
) {

        if (nrow(x = rows) == 0) return(rows)
        steps_sum <-  vapply(X = seq_len(length.out = nrow(x = rows)), FUN.VALUE = numeric(1), FUN = function(row_index) {
              run_steps <-  leapfrog_steps[leapfrog_steps$run_number == rows$run_number[row_index] &
                                           leapfrog_steps$iteration  <= rows$n_timed_iters[row_index], , drop = FALSE]
              if (nrow(x = run_steps) != rows$n_timed_iters[row_index]) {
                    stop("fn_paper1_burnin_add_sec_per_step: leapfrog steps missing for run ", rows$run_number[row_index],
                         ", iterations 1-", rows$n_timed_iters[row_index])
              }
              sum(run_steps$n_leapfrog_steps)
        })
        timed_seconds <-  rows$total_timed_seconds
        if ("sec_per_iter_each" %in% names(x = rows)) {
              has_each <-  !is.na(x = rows$sec_per_iter_each) & nzchar(x = rows$sec_per_iter_each)
              timed_seconds[has_each] <-  vapply(X = strsplit(x = rows$sec_per_iter_each[has_each], split = ",", fixed = TRUE),
                                                 FUN = function(each_seconds) sum(as.numeric(x = each_seconds)), FUN.VALUE = numeric(1))
              if (any(has_each)) {
                    relative_difference <-  max(abs(timed_seconds[has_each] - rows$total_timed_seconds[has_each]) / rows$total_timed_seconds[has_each])
                    fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_add_sec_per_step: largest relative difference between the sum of ",
                                                           "sec_per_iter_each and total_timed_seconds: ", formatC(x = relative_difference, digits = 3, format = "g")))
              }
        }
        rows$sec_per_step   <-  timed_seconds / steps_sum
        rows$steps_per_iter <-  steps_sum / rows$n_timed_iters
        return(rows)

}


##
## ---- fn_paper1_burnin_best_by_cell: per (N, burn-in chains), the best configuration and both speed-ups (as in the table helper): --------
##
fn_paper1_burnin_best_by_cell <-  function( rows,
                                            device,
                                            implementation
) {

        configuration_summary <-  fn_paper1_burnin_summarise_rows(rows = rows, device = device)
        cell_columns          <-  c("N", "n_chains_burnin")
        ##
        best_overall <-  configuration_summary %>%
              dplyr::group_by(dplyr::across(.cols = dplyr::all_of(x = cell_columns))) %>%
              dplyr::slice(which.min(.data$sec_per_iter)) %>%
              dplyr::ungroup() %>%
              dplyr::transmute( N                         = .data$N,
                                n_chains_burnin           = .data$n_chains_burnin,
                                best_num_chunks_burnin    = .data$num_chunks_burnin,
                                best_n_threads_WCP_burnin = .data$n_threads_WCP_burnin,
                                best_n_threads_total      = .data$n_threads_total,
                                best_sec_per_iter         = .data$sec_per_iter,
                                best_n_runs_available     = .data$n_runs_available)
        ##
        ## ---- seconds per leapfrog step and leapfrog steps per timed iteration of the best configuration (mean over runs), where available:
        ##
        if ("sec_per_step" %in% names(x = rows)) {
              step_summary <-  rows %>%
                    dplyr::group_by(dplyr::across(.cols = dplyr::all_of(x = c("N", "n_chains_burnin", "num_chunks_burnin", "n_threads_WCP_burnin")))) %>%
                    dplyr::summarise( best_sec_per_step        = mean(x = .data$sec_per_step),
                                      best_mean_steps_per_iter = mean(x = .data$steps_per_iter),
                                      best_n_timed_iters       = max(.data$n_timed_iters),
                                      .groups                  = "drop") %>%
                    dplyr::rename( best_num_chunks_burnin    = "num_chunks_burnin",
                                   best_n_threads_WCP_burnin = "n_threads_WCP_burnin")
              best_overall <-  dplyr::left_join(x = best_overall, y = step_summary,
                                                by = c("N", "n_chains_burnin", "best_num_chunks_burnin", "best_n_threads_WCP_burnin"))
        }
        ##
        ## ---- reference (i): N_chunks = 1, N_WCP = 1:
        ##
        no_chunking <-  configuration_summary %>%
              dplyr::filter(.data$num_chunks_burnin == 1, .data$n_threads_WCP_burnin == 1) %>%
              dplyr::transmute( N                        = .data$N,
                                n_chains_burnin          = .data$n_chains_burnin,
                                no_chunking_sec_per_iter = .data$sec_per_iter)
        ##
        ## ---- reference (ii): best N_chunks with N_WCP = 1:
        ##
        best_chunking_only <-  configuration_summary %>%
              dplyr::filter(.data$n_threads_WCP_burnin == 1) %>%
              dplyr::group_by(dplyr::across(.cols = dplyr::all_of(x = cell_columns))) %>%
              dplyr::slice(which.min(.data$sec_per_iter)) %>%
              dplyr::ungroup() %>%
              dplyr::transmute( N                               = .data$N,
                                n_chains_burnin                 = .data$n_chains_burnin,
                                best_chunking_only_num_chunks   = .data$num_chunks_burnin,
                                best_chunking_only_sec_per_iter = .data$sec_per_iter)
        ##
        best_by_cell <-  best_overall %>%
              dplyr::left_join(y = no_chunking,        by = cell_columns) %>%
              dplyr::left_join(y = best_chunking_only, by = cell_columns) %>%
              dplyr::mutate( speed_up_vs_1_chunk_WCP_1     = .data$no_chunking_sec_per_iter        / .data$best_sec_per_iter,
                             speed_up_vs_best_chunks_WCP_1 = .data$best_chunking_only_sec_per_iter / .data$best_sec_per_iter,
                             device                        = device,
                             implementation                = implementation) %>%
              dplyr::arrange(.data$N, .data$n_chains_burnin)
        ##
        if (any(is.na(x = best_by_cell$no_chunking_sec_per_iter))) {
              fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_best_by_cell: ", implementation, " (", device, "): no N_chunks = 1, N_WCP = 1 ",
                                                     "measurement for ", sum(is.na(x = best_by_cell$no_chunking_sec_per_iter)),
                                                     " (N, chains) cell(s); their total speed-up is n/a."))
        }
        ##
        return(list(configuration_summary = configuration_summary, best_by_cell = best_by_cell))

}


##
## ---- fn_paper1_burnin_figure_Stan_panel (private helper): one (N, burn-in chains) panel of fn_paper1_burnin_figure_Stan: ---------------
##
## As fn_paper1_burnin_figure_best_chunks_panel() in R_fns_alg_paper_1_burnin_report.R, plus a dotted grey line at the best measured
## NicoStan+BayesMVP seconds per burn-in iteration for the same (N, burn-in chains), and a three-line speed-up label.
##
fn_paper1_burnin_figure_Stan_panel <-  function( panel_best_by_WCP,
                                                 panel_no_chunking,
                                                 panel_best_overall,
                                                 panel_BayesMVP_best,
                                                 panel_speed_up_label,
                                                 show_N_strip,
                                                 show_chains_strip,
                                                 series_levels,
                                                 series_labels,
                                                 series_colours,
                                                 series_linetypes,
                                                 series_shapes,
                                                 best_overall_label,
                                                 panel_WCP_only = NULL
) {

        ##
        ## ---- x breaks at every N_WCP value measured in this panel; the tick labels are thinned (first and last always kept) so that
        ##      no two labels are closer than 11% of the panel's log2 range - the unlabelled ticks stay on the axis:
        ##
        ## measured_WCP_values <-  sort(x = unique(x = c(panel_best_by_WCP$n_threads_WCP_burnin, panel_no_chunking$n_threads_WCP_burnin)))
        measured_WCP_values <-  sort(x = unique(x = c(panel_best_by_WCP$n_threads_WCP_burnin, panel_no_chunking$n_threads_WCP_burnin,
                                                      if (is.null(x = panel_WCP_only)) NULL else panel_WCP_only$n_threads_WCP_burnin)))
        measured_WCP_log2   <-  log2(x = measured_WCP_values)
        minimum_log2_gap    <-  if (length(x = measured_WCP_values) > 1) 0.11 * diff(x = range(measured_WCP_log2)) else 0
        is_labelled_tick    <-  rep(x = FALSE, times = length(x = measured_WCP_values))
        is_labelled_tick[c(1, length(x = measured_WCP_values))] <-  TRUE
        last_labelled_log2  <-  measured_WCP_log2[1]
        for (tick_index in seq_along(along.with = measured_WCP_values)[-c(1, length(x = measured_WCP_values))]) {
              if (measured_WCP_log2[tick_index] - last_labelled_log2 >= minimum_log2_gap &&
                  measured_WCP_log2[length(x = measured_WCP_values)] - measured_WCP_log2[tick_index] >= minimum_log2_gap) {
                    is_labelled_tick[tick_index] <-  TRUE
                    last_labelled_log2 <-  measured_WCP_log2[tick_index]
              }
        }
        tick_labels <-  ifelse(test = is_labelled_tick, yes = as.character(x = measured_WCP_values), no = "")
        ##
        ## ---- log10 y breaks, thinned from the top down so that no two labels are closer than 9% of the panel's log10 range:
        ##
        y_breaks_function <-  function(y_limits) {

              candidate_breaks <-  scales::breaks_log(n = 6)(y_limits)
              candidate_breaks <-  sort(x = candidate_breaks[candidate_breaks >= min(y_limits) & candidate_breaks <= max(y_limits)], decreasing = TRUE)
              if (length(x = candidate_breaks) < 2) return(candidate_breaks)
              minimum_log10_gap <-  0.09 * diff(x = log10(x = range(y_limits)))
              kept_breaks <-  candidate_breaks[1]
              for (candidate_break in candidate_breaks[-1]) {
                    if (log10(x = kept_breaks[length(x = kept_breaks)]) - log10(x = candidate_break) >= minimum_log10_gap) {
                          kept_breaks <-  c(kept_breaks, candidate_break)
                    }
              }
              return(sort(x = kept_breaks))

        }
        ##
        ## ---- N_chunks labels below each point, above a point that is slower than its left-hand neighbour, further below the outlined best point:
        ##
        panel_best_by_WCP   <-  panel_best_by_WCP[order(panel_best_by_WCP$n_threads_WCP_burnin), , drop = FALSE]
        rises_from_left     <-  c(FALSE, diff(x = panel_best_by_WCP$sec_per_iter) > 0)
        is_best_overall     <-  panel_best_by_WCP$n_threads_WCP_burnin == panel_best_overall$n_threads_WCP_burnin[1]
        panel_best_by_WCP$label_vjust <-  ifelse(test = rises_from_left, yes = -1.0, no = ifelse(test = is_best_overall, yes = 3.0, no = 2.0))
        ##
        ## ---- N_chunks point labels thinned where they would overlap (e.g. on a flat tail): the best point and the first point are always
        ##      labelled; any other point is labelled only if, against every labelled point, it is at least 18% of the panel's log2 x range
        ##      away, or at least 15% of the panel's log10 y range (of the plotted values) above or below:
        ##
        point_log2        <-  log2(x = panel_best_by_WCP$n_threads_WCP_burnin)
        point_log10       <-  log10(x = panel_best_by_WCP$sec_per_iter)
        panel_log10_range <-  diff(x = range(log10(x = c(panel_best_by_WCP$sec_per_iter, panel_no_chunking$sec_per_iter,
                                                         if (is.null(x = panel_BayesMVP_best)) NULL else panel_BayesMVP_best$sec_per_iter))))
        minimum_label_log2_gap  <-  if (length(x = point_log2) > 1) 0.18 * diff(x = range(point_log2)) else 0
        minimum_label_log10_gap <-  0.15 * panel_log10_range
        is_labelled       <-  is_best_overall | seq_along(along.with = point_log2) == 1
        for (point_index in seq_along(along.with = point_log2)) {
              if (!is_labelled[point_index]) {
                    far_enough <-  abs(point_log2[point_index]  - point_log2[is_labelled])  >= minimum_label_log2_gap |
                                   abs(point_log10[point_index] - point_log10[is_labelled]) >= minimum_label_log10_gap
                    if (all(far_enough)) is_labelled[point_index] <-  TRUE
              }
        }
        panel_point_labels <-  panel_best_by_WCP[is_labelled, , drop = FALSE]
        ##
        best_overall_fill <-  c(NA)
        names(x = best_overall_fill) <-  best_overall_label
        ##
        panel_plot <-  ggplot2::ggplot()
        ##
        ## ---- dotted level of the best measured NicoStan+BayesMVP configuration (same N and burn-in chains):
        ##
        if (!is.null(x = panel_BayesMVP_best) && nrow(x = panel_BayesMVP_best) > 0) {
              panel_plot <-  panel_plot +
                    ggplot2::geom_hline( data      = panel_BayesMVP_best,
                                         mapping   = ggplot2::aes(yintercept = .data$sec_per_iter,
                                                                  colour     = .data$series,
                                                                  linetype   = .data$series),
                                         linewidth = 0.7)
        }
        ##
        ## ---- dashed reference level of N_chunks = 1, N_WCP = 1 (the standard Stan model), where measured:
        ##
        if (nrow(x = panel_no_chunking) > 0) {
              panel_plot <-  panel_plot +
                    ggplot2::geom_hline( data      = panel_no_chunking,
                                         mapping   = ggplot2::aes(yintercept = .data$sec_per_iter,
                                                                  colour     = .data$series,
                                                                  linetype   = .data$series),
                                         linewidth = 0.6)
        }
        ##
        ## ---- best N_chunks at each N_WCP:
        ##
        panel_plot <-  panel_plot +
              ggplot2::geom_line( data      = panel_best_by_WCP,
                                  mapping   = ggplot2::aes(x        = .data$n_threads_WCP_burnin,
                                                           y        = .data$sec_per_iter,
                                                           colour   = .data$series,
                                                           linetype = .data$series),
                                  linewidth = 0.8) +
              ggplot2::geom_point( data     = panel_best_by_WCP,
                                   mapping  = ggplot2::aes(x      = .data$n_threads_WCP_burnin,
                                                           y      = .data$sec_per_iter,
                                                           colour = .data$series,
                                                           shape  = .data$series),
                                   size     = 2.2) +
              ggplot2::geom_text( data        = panel_point_labels,
                                  mapping     = ggplot2::aes(x     = .data$n_threads_WCP_burnin,
                                                             y     = .data$sec_per_iter,
                                                             label = .data$num_chunks_burnin,
                                                             vjust = .data$label_vjust),
                                  colour      = series_colours[["best_chunks_at_each_WCP"]],
                                  size        = 3.3,
                                  show.legend = FALSE)
        ##
        ## ---- WCP-only (N_chunks = N_threads/chain), where measured:
        ##
        if (!is.null(x = panel_WCP_only) && nrow(x = panel_WCP_only) > 0) {
              panel_plot <-  panel_plot +
                    ggplot2::geom_line( data      = panel_WCP_only,
                                        mapping   = ggplot2::aes(x        = .data$n_threads_WCP_burnin,
                                                                 y        = .data$sec_per_iter,
                                                                 colour   = .data$series,
                                                                 linetype = .data$series),
                                        linewidth = 0.8) +
                    ggplot2::geom_point( data     = panel_WCP_only,
                                         mapping  = ggplot2::aes(x      = .data$n_threads_WCP_burnin,
                                                                 y      = .data$sec_per_iter,
                                                                 colour = .data$series,
                                                                 shape  = .data$series),
                                         size     = 2.6)
        }
        ##
        ## ---- N_chunks = 1, N_WCP = 1 measurement (drawn after the best line, so it stays visible where the two coincide):
        ##
        if (nrow(x = panel_no_chunking) > 0) {
              panel_plot <-  panel_plot +
                    ggplot2::geom_point( data    = panel_no_chunking,
                                         mapping = ggplot2::aes(x      = .data$n_threads_WCP_burnin,
                                                                y      = .data$sec_per_iter,
                                                                colour = .data$series,
                                                                shape  = .data$series),
                                         size    = 3.4,
                                         stroke  = 1.1)
        }
        ##
        ## ---- best measured (N_chunks, N_WCP) configuration, outlined, and the speed-up label:
        ##
        panel_plot <-  panel_plot +
              ggplot2::geom_point( data    = panel_best_overall,
                                   mapping = ggplot2::aes(x    = .data$n_threads_WCP_burnin,
                                                          y    = .data$sec_per_iter,
                                                          fill = .data$best_overall_series),
                                   shape   = 21,
                                   size    = 5.0,
                                   stroke  = 1.0,
                                   colour  = "black") +
              ggplot2::annotate( geom  = "text",
                                 x     = Inf,
                                 y     = Inf,
                                 label = panel_speed_up_label,
                                 hjust = 1.06,
                                 vjust = 1.20,
                                 size  = 3.4) +
              ggplot2::scale_x_continuous( transform    = "log2",
                                           breaks       = measured_WCP_values,
                                           labels       = tick_labels,
                                           minor_breaks = NULL,
                                           expand       = ggplot2::expansion(mult = c(0.08, 0.16))) +
              ggplot2::scale_y_log10( breaks = y_breaks_function,
                                      labels = function(y_breaks) trimws(x = formatC(x = y_breaks, digits = 3, format = "fg")),
                                      expand = ggplot2::expansion(mult = c(0.20, 0.75)),
                                      guide  = ggplot2::guide_axis(check.overlap = TRUE)) +
              ggplot2::scale_colour_manual(   name = NULL, values = series_colours,   breaks = series_levels, limits = series_levels, labels = series_labels) +
              ggplot2::scale_linetype_manual( name = NULL, values = series_linetypes, breaks = series_levels, limits = series_levels, labels = series_labels) +
              ggplot2::scale_shape_manual(    name = NULL, values = series_shapes,    breaks = series_levels, limits = series_levels, labels = series_labels,
                                              na.translate = FALSE) +
              ggplot2::scale_fill_manual(     name = NULL, values = best_overall_fill, na.value = NA) +
              ggplot2::guides( colour   = ggplot2::guide_legend(order = 1),
                               linetype = ggplot2::guide_legend(order = 1),
                               shape    = ggplot2::guide_legend(order = 1),
                               fill     = ggplot2::guide_legend(order = 2)) +
              ## ggplot2::labs( x = expression(N[WCP] ~ "(within-chain-parallelism threads per burn-in chain)"),
              ggplot2::labs( x = expression(N["threads/chain"] ~ "(threads of each burn-in chain)"),
                             y = "Seconds per burn-in iteration (log scale)") +
              ggplot2::theme_bw(base_size = 14) +
              ggplot2::theme( legend.position  = "bottom",
                              legend.direction = "vertical",
                              legend.key.width = ggplot2::unit(x = 1.6, units = "lines"))
        ##
        ## ---- strips on the outer edges of the grid only (N along the top row, burn-in chains down the right-hand column);
        ##      the device / burn-in chain strips are plotmath, so they are parsed:
        ##
        if (show_N_strip && show_chains_strip) {
              panel_plot <-  panel_plot + ggplot2::facet_grid(rows = ggplot2::vars(.data$chains_label), cols = ggplot2::vars(.data$N_label),
                                                              labeller = ggplot2::labeller(.rows = ggplot2::label_parsed))
        } else if (show_N_strip) {
              panel_plot <-  panel_plot + ggplot2::facet_grid(cols = ggplot2::vars(.data$N_label))
        } else if (show_chains_strip) {
              panel_plot <-  panel_plot + ggplot2::facet_grid(rows = ggplot2::vars(.data$chains_label),
                                                              labeller = ggplot2::labeller(.rows = ggplot2::label_parsed))
        }
        ##
        return(panel_plot)

}


##
## ---- fn_paper1_burnin_figure_Stan: tape-chunked Stan seconds per burn-in iteration vs N_WCP, per (device, chains) x N, with BayesMVP: ---
##
## One figure for all devices (rows = device x burn-in chains, columns = N, one ggplot per panel), with the layout and references of
## fn_paper1_burnin_figure_best_chunks():
##   - "Best N_chunks at each N_WCP" = solid line, points labelled with N_chunks; its N_WCP = 1 point is the best chunking-only setting;
##   - "Best measured NicoStan+BayesMVP configuration" = dotted grey line at BayesMVP's best seconds per leapfrog step times the mean
##     leapfrog steps per timed Stan iteration (i.e., BayesMVP's time for the same trajectory lengths; same N, device and burn-in chains);
##   - "No chunking, no WCP" (N_chunks = 1, N_WCP = 1) and the "total" label line only if such a run exists (the Stan chunk grid has none);
##   - label: "WCP" = best chunking-only / best, "BayesMVP" = Stan best / BayesMVP best seconds per leapfrog step.
##
fn_paper1_burnin_figure_Stan <-  function( Stan_rows_by_device,
                                           Stan_best_by_cell_by_device,
                                           BayesMVP_best_by_cell_by_device,
                                           file_path
) {

        device_full_name <-  c(HPC = "local-HPC", Laptop = "laptop")
        devices          <-  intersect(x = c("HPC", "Laptop"), y = names(x = Stan_rows_by_device))
        devices          <-  devices[vapply(X = devices, FUN = function(device) !is.null(x = Stan_rows_by_device[[device]]) && nrow(x = Stan_rows_by_device[[device]]) > 0,
                                            FUN.VALUE = logical(1))]
        if (length(x = devices) == 0) stop("fn_paper1_burnin_figure_Stan: no tape-chunked Stan rows for any device.")
        ##
        ## ---- per-device configuration summaries, with the row key (device x burn-in chains):
        ##
        configuration_summary <-  dplyr::bind_rows(lapply(X = devices, FUN = function(device) {
              device_summary <-  fn_paper1_burnin_summarise_rows(rows = Stan_rows_by_device[[device]], device = device)
              device_summary$device_selector <-  device
              device_summary
        }))
        ##
        N_values  <-  sort(x = unique(x = configuration_summary$N))
        row_keys  <-  unique(x = configuration_summary[order(match(x = configuration_summary$device_selector, table = devices), configuration_summary$n_chains_burnin),
                                                       c("device_selector", "n_chains_burnin")])
        ## row_labels <-  paste0(device_full_name[row_keys$device_selector], "\n", row_keys$n_chains_burnin, " burn-in chains")
        ## (plotmath strip labels: the device above N_chains = ..., parsed by the panel's labeller)
        ## row_labels <-  paste0('atop("', device_full_name[row_keys$device_selector], '", N[chains]==', row_keys$n_chains_burnin, '~"(burn-in)")')
        row_labels <-  paste0( 'atop("', device_full_name[row_keys$device_selector], '", N["burn_chains"]==',
                               row_keys$n_chains_burnin, ')')
        ##
        configuration_summary$N_label      <-  factor(x = paste0("N = ", fn_paper1_format_number_commas_from_10000(configuration_summary$N)),
                                                      levels = paste0("N = ", fn_paper1_format_number_commas_from_10000(N_values)))
        ## configuration_summary$chains_label <-  factor(x = paste0(device_full_name[configuration_summary$device_selector], "\n", configuration_summary$n_chains_burnin, " burn-in chains"),
        ##                                               levels = row_labels)
        ## configuration_summary$chains_label <-  factor(x = paste0('atop("', device_full_name[configuration_summary$device_selector], '", N[chains]==',
        ##                                                          configuration_summary$n_chains_burnin, '~"(burn-in)")'),
        ##                                               levels = row_labels)
        configuration_summary$chains_label <-  factor( x = paste0( 'atop("', device_full_name[configuration_summary$device_selector],
                                                                   '", N["burn_chains"]==',
                                                                   configuration_summary$n_chains_burnin, ')'),
                                                       levels = row_labels)
        ##
        no_chunking <-  configuration_summary[configuration_summary$num_chunks_burnin == 1 & configuration_summary$n_threads_WCP_burnin == 1, , drop = FALSE]
        ##
        include_no_chunking_series <-  nrow(x = no_chunking) > 0
        include_BayesMVP_series    <-  any(vapply(X = devices, FUN = function(device) !is.null(x = BayesMVP_best_by_cell_by_device[[device]]), FUN.VALUE = logical(1)))
        ##
        ## all_series_levels    <-  c("no_chunking", "best_chunks_at_each_WCP", "BayesMVP_best")
        ## all_series_labels    <-  list( no_chunking             = expression("No chunking, no WCP (" * N[chunks] * " = 1, " * N[WCP] * " = 1)"),
        ##                                best_chunks_at_each_WCP = expression("Best " * N[chunks] * " at each " * N[WCP] * " (point labels give " * N[chunks] * ")"),
        ##                                BayesMVP_best           = expression("Best measured NicoStan+BayesMVP configuration"))
        ## all_series_colours   <-  c(no_chunking = "#D55E00", best_chunks_at_each_WCP = "#0072B2", BayesMVP_best = "grey45")
        ## all_series_linetypes <-  c(no_chunking = "dashed",  best_chunks_at_each_WCP = "solid",   BayesMVP_best = "dotted")
        ## all_series_shapes    <-  c(no_chunking = 0,         best_chunks_at_each_WCP = 16,        BayesMVP_best = NA)
        ##
        ## series_levels      <-  all_series_levels[c(include_no_chunking_series, TRUE, include_BayesMVP_series)]
        ## (WCP-only, i.e. N_chunks = N_threads/chain, is its own series, drawn in green with triangles)
        WCP_only <-  configuration_summary[configuration_summary$num_chunks_burnin == configuration_summary$n_threads_WCP_burnin &
                                           configuration_summary$n_threads_WCP_burnin > 1, , drop = FALSE]
        include_WCP_only_series <-  nrow(x = WCP_only) > 0
        ##
        all_series_levels    <-  c("no_chunking", "best_chunks_at_each_WCP", "WCP_only", "BayesMVP_best")
        all_series_labels    <-  list( no_chunking             = expression("AD_Stan (" * N[chunks] * " = 1, " * N["threads/chain"] * " = 1)"),
                                       best_chunks_at_each_WCP = expression("Best " * N[chunks] * " at each " * N["threads/chain"] * " (point labels give " * N[chunks] * ")"),
                                       WCP_only                = expression("AD_Stan_WCP (WCP-only: " * N[chunks] * " = " * N["threads/chain"] * ")"),
                                       BayesMVP_best           = expression("Best measured NicoStan+BayesMVP configuration"))
        all_series_colours   <-  c(no_chunking = "#D55E00", best_chunks_at_each_WCP = "#0072B2", WCP_only = "#009E73", BayesMVP_best = "grey45")
        all_series_linetypes <-  c(no_chunking = "dashed",  best_chunks_at_each_WCP = "solid",   WCP_only = "22",      BayesMVP_best = "dotted")
        all_series_shapes    <-  c(no_chunking = 0,         best_chunks_at_each_WCP = 16,        WCP_only = 17,        BayesMVP_best = NA)
        ##
        series_levels      <-  all_series_levels[c(include_no_chunking_series, TRUE, include_WCP_only_series, include_BayesMVP_series)]
        series_labels      <-  do.call(what = c, args = unname(obj = all_series_labels[series_levels]))
        series_colours     <-  all_series_colours[series_levels]
        series_linetypes   <-  all_series_linetypes[series_levels]
        series_shapes      <-  all_series_shapes[series_levels]
        best_overall_label <-  "Best measured configuration"
        ##
        no_chunking$series <-  factor(x = rep(x = "no_chunking", times = nrow(x = no_chunking)), levels = all_series_levels)
        WCP_only$series    <-  factor(x = rep(x = "WCP_only",    times = nrow(x = WCP_only)),    levels = all_series_levels)
        ##
        best_by_WCP <-  configuration_summary %>%
              dplyr::group_by(dplyr::across(.cols = dplyr::all_of(x = c("device_selector", "N", "n_chains_burnin", "n_threads_WCP_burnin")))) %>%
              dplyr::slice(which.min(.data$sec_per_iter)) %>%
              dplyr::ungroup()
        best_by_WCP$series <-  factor(x = rep(x = "best_chunks_at_each_WCP", times = nrow(x = best_by_WCP)), levels = all_series_levels)
        ##
        best_overall <-  configuration_summary %>%
              dplyr::group_by(dplyr::across(.cols = dplyr::all_of(x = c("device_selector", "N", "n_chains_burnin")))) %>%
              dplyr::slice(which.min(.data$sec_per_iter)) %>%
              dplyr::ungroup()
        best_overall$best_overall_series <-  rep(x = best_overall_label, times = nrow(x = best_overall))
        ##
        panel_list    <-  list()
        panel_summary <-  list()
        ##
        for (row_index in seq_len(length.out = nrow(x = row_keys))) {
              for (N_index in seq_along(along.with = N_values)) {

                    this_device <-  row_keys$device_selector[row_index]
                    this_chains <-  row_keys$n_chains_burnin[row_index]
                    this_N      <-  N_values[N_index]
                    ##
                    in_panel <-  function(summary_rows) summary_rows$device_selector == this_device & summary_rows$N == this_N & summary_rows$n_chains_burnin == this_chains
                    panel_best_by_WCP  <-  best_by_WCP[in_panel(best_by_WCP), , drop = FALSE]
                    panel_no_chunking  <-  no_chunking[in_panel(no_chunking), , drop = FALSE]
                    panel_best_overall <-  best_overall[in_panel(best_overall), , drop = FALSE]
                    panel_WCP_only     <-  WCP_only[in_panel(WCP_only), , drop = FALSE]
                    ##
                    if (nrow(x = panel_best_by_WCP) == 0) {
                          fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_figure_Stan: no rows for ", this_device, ", N = ", this_N, ", ",
                                                                 this_chains, " burn-in chains; leaving that panel empty."))
                          panel_list[[length(x = panel_list) + 1]] <-  patchwork::plot_spacer()
                          next
                    }
                    ##
                    panel_BayesMVP_best   <-  NULL
                    BayesMVP_sec_per_iter <-  NA_real_
                    Stan_over_BayesMVP    <-  NA_real_
                    BayesMVP_best_by_cell <-  BayesMVP_best_by_cell_by_device[[this_device]]
                    Stan_best_by_cell     <-  Stan_best_by_cell_by_device[[this_device]]
                    if (!is.null(x = BayesMVP_best_by_cell) && !is.null(x = Stan_best_by_cell)) {
                          BayesMVP_cell <-  BayesMVP_best_by_cell[BayesMVP_best_by_cell$N == this_N & BayesMVP_best_by_cell$n_chains_burnin == this_chains, , drop = FALSE]
                          Stan_cell     <-  Stan_best_by_cell[Stan_best_by_cell$N == this_N & Stan_best_by_cell$n_chains_burnin == this_chains, , drop = FALSE]
                          if (nrow(x = BayesMVP_cell) == 1 && nrow(x = Stan_cell) == 1) {
                                ##
                                ## ---- BayesMVP's seconds for the trajectory lengths of the timed Stan iterations, and the per-step ratio:
                                ##
                                BayesMVP_sec_per_iter <-  BayesMVP_cell$best_sec_per_step * Stan_cell$best_mean_steps_per_iter
                                Stan_over_BayesMVP    <-  Stan_cell$best_sec_per_step / BayesMVP_cell$best_sec_per_step
                                panel_BayesMVP_best   <-  data.frame( sec_per_iter = BayesMVP_sec_per_iter,
                                                                      series       = factor(x = "BayesMVP_best", levels = all_series_levels),
                                                                      N_label      = panel_best_overall$N_label[1],
                                                                      chains_label = panel_best_overall$chains_label[1])
                          }
                    }
                    ##
                    ## ---- speed-ups, computed as in fn_paper1_burnin_best_table_tex(), and the ratio to BayesMVP:
                    ##
                    best_sec_per_iter         <-  panel_best_overall$sec_per_iter
                    best_chunking_only_sec    <-  panel_best_by_WCP$sec_per_iter[panel_best_by_WCP$n_threads_WCP_burnin == 1]
                    no_chunking_sec           <-  panel_no_chunking$sec_per_iter
                    speed_up_vs_1_chunk_WCP_1 <-  if (length(x = no_chunking_sec) == 1)        no_chunking_sec        / best_sec_per_iter else NA_real_
                    speed_up_vs_best_chunks   <-  if (length(x = best_chunking_only_sec) == 1) best_chunking_only_sec / best_sec_per_iter else NA_real_
                    ##
                    label_lines <-  character(0)
                    if (include_no_chunking_series) label_lines <-  c(label_lines, paste0("total ", fn_paper1_burnin_format_fold(fold_vector = speed_up_vs_1_chunk_WCP_1)))
                    label_lines <-  c(label_lines, paste0("WCP ", fn_paper1_burnin_format_fold(fold_vector = speed_up_vs_best_chunks)))
                    if (include_BayesMVP_series)    label_lines <-  c(label_lines, paste0("BayesMVP ", fn_paper1_burnin_format_fold(fold_vector = Stan_over_BayesMVP)))
                    panel_speed_up_label <-  paste(label_lines, collapse = "\n")
                    ##
                    panel_summary[[length(x = panel_summary) + 1]] <-  data.frame( device                     = this_device,
                                                                                   N                          = this_N,
                                                                                   n_chains_burnin            = this_chains,
                                                                                   best_num_chunks_burnin     = panel_best_overall$num_chunks_burnin,
                                                                                   best_n_threads_WCP_burnin  = panel_best_overall$n_threads_WCP_burnin,
                                                                                   best_sec_per_iter          = best_sec_per_iter,
                                                                                   speed_up_vs_1_chunk_WCP_1  = speed_up_vs_1_chunk_WCP_1,
                                                                                   speed_up_vs_best_chunks    = speed_up_vs_best_chunks,
                                                                                   BayesMVP_best_sec_per_iter = BayesMVP_sec_per_iter,
                                                                                   Stan_over_BayesMVP         = Stan_over_BayesMVP,
                                                                                   panel_speed_up_label       = panel_speed_up_label)
                    ##
                    panel_list[[length(x = panel_list) + 1]] <-  fn_paper1_burnin_figure_Stan_panel( panel_best_by_WCP    = panel_best_by_WCP,
                                                                                                     panel_no_chunking    = panel_no_chunking,
                                                                                                     panel_best_overall   = panel_best_overall,
                                                                                                     panel_BayesMVP_best  = panel_BayesMVP_best,
                                                                                                     panel_speed_up_label = panel_speed_up_label,
                                                                                                     show_N_strip         = (row_index == 1),
                                                                                                     show_chains_strip    = (N_index == length(x = N_values)),
                                                                                                     series_levels        = series_levels,
                                                                                                     series_labels        = series_labels,
                                                                                                     series_colours       = series_colours,
                                                                                                     series_linetypes     = series_linetypes,
                                                                                                     series_shapes        = series_shapes,
                                                                                                     best_overall_label   = best_overall_label,
                                                                                                     panel_WCP_only       = panel_WCP_only)

              }
        }
        ##
        ## ---- one collected legend at the bottom, one x-axis title and one y-axis title (as for the NicoStan+BayesMVP figure):
        ##
        burnin_figure <-  patchwork::wrap_plots(panel_list, ncol = length(x = N_values), byrow = TRUE) +
              patchwork::plot_layout(guides = "collect", axis_titles = "collect") &
              ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")
        ##
        dir.create(path = dirname(path = file_path), recursive = TRUE, showWarnings = FALSE)
        ggplot2::ggsave( filename = file_path,
                         plot     = burnin_figure,
                         width    = 10,
                         height   = 2.2 * nrow(x = row_keys) + 1.6,
                         dpi      = 150)
        fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_figure_Stan: wrote figure (", length(x = panel_summary), " panel(s), devices: ",
                                               paste(devices, collapse = ", "), ") to ", file_path))
        ##
        return(invisible(list(figure = burnin_figure, panel_summary = do.call(what = rbind, args = panel_summary))))

}


##
## ---- fn_paper1_burnin_placeholder_values: the value of every \BurninTBD{ID} placeholder (draft section and Part VI sentence): -----------
##
## Formats: speed-ups and ratios with 2 decimal places (as the figure labels and the existing text), seconds with 3 significant digits
## (fn_paper1_burnin_format_seconds_per_iter), chain counts as words (as the existing text: "with four burn-in chains").
##
fn_paper1_burnin_placeholder_values <-  function( Stan_best_by_cell_by_device,
                                                  BayesMVP_best_by_cell_by_device
) {

        fold_text    <-  function(value) if (length(x = value) == 1 && is.finite(value)) formatC(x = value, digits = 2, format = "f") else NA_character_
        seconds_text <-  function(value) if (length(x = value) == 1 && is.finite(value)) fn_paper1_burnin_format_seconds_per_iter(seconds_per_iter_vector = value) else NA_character_
        integer_text <-  function(value) if (length(x = value) == 1 && is.finite(value)) format(x = value, scientific = FALSE, trim = TRUE) else NA_character_
        range_text   <-  function(values, which_end) {
              finite_values <-  values[is.finite(values)]
              if (length(x = finite_values) == 0) return(NA_character_)
              fold_text(value = if (which_end == "min") min(finite_values) else max(finite_values))
        }
        cell_value   <-  function(best_by_cell, N, n_chains_burnin, column_name) {
              if (is.null(x = best_by_cell) || nrow(x = best_by_cell) == 0) return(NA_real_)
              cell_rows <-  best_by_cell[best_by_cell$N == N & best_by_cell$n_chains_burnin == n_chains_burnin, , drop = FALSE]
              if (nrow(x = cell_rows) != 1) return(NA_real_)
              return(as.numeric(x = cell_rows[[column_name]]))
        }
        Stan_over_BayesMVP <-  function(device) {
              Stan_best     <-  Stan_best_by_cell_by_device[[device]]
              BayesMVP_best <-  BayesMVP_best_by_cell_by_device[[device]]
              if (is.null(x = Stan_best) || is.null(x = BayesMVP_best)) return(numeric(0))
              both_implementations <-  dplyr::inner_join( x      = dplyr::select(.data = Stan_best,     dplyr::all_of(x = c("N", "n_chains_burnin", "best_sec_per_step"))),
                                                          y      = dplyr::select(.data = BayesMVP_best, dplyr::all_of(x = c("N", "n_chains_burnin", "best_sec_per_step"))),
                                                          by     = c("N", "n_chains_burnin"),
                                                          suffix = c("_Stan", "_BayesMVP"))
              return(both_implementations$best_sec_per_step_Stan / both_implementations$best_sec_per_step_BayesMVP)
        }
        ##
        values <-  list()
        add_value <-  function(id, value, source_description) {
              values[[length(x = values) + 1]] <<-  data.frame(id = id, value = value, source = source_description, stringsAsFactors = FALSE)
        }
        ##
        ## ---- NicoStan+BayesMVP (Part VI; the HPC values change only through the added N = 50,000 runs, the laptop values are checks):
        ##
        for (device in names(x = BayesMVP_best_by_cell_by_device)) {

              BayesMVP_best <-  BayesMVP_best_by_cell_by_device[[device]]
              if (is.null(x = BayesMVP_best) || nrow(x = BayesMVP_best) == 0) next
              add_value(paste0("BayesMVP_", device, "_total_min"), range_text(BayesMVP_best$speed_up_vs_1_chunk_WCP_1,     "min"), paste0(device, ": min over cells, best vs N_chunks = 1, N_WCP = 1"))
              add_value(paste0("BayesMVP_", device, "_total_max"), range_text(BayesMVP_best$speed_up_vs_1_chunk_WCP_1,     "max"), paste0(device, ": max over cells, best vs N_chunks = 1, N_WCP = 1"))
              add_value(paste0("BayesMVP_", device, "_WCP_min"),   range_text(BayesMVP_best$speed_up_vs_best_chunks_WCP_1, "min"), paste0(device, ": min over cells, best vs best N_chunks at N_WCP = 1"))
              add_value(paste0("BayesMVP_", device, "_WCP_max"),   range_text(BayesMVP_best$speed_up_vs_best_chunks_WCP_1, "max"), paste0(device, ": max over cells, best vs best N_chunks at N_WCP = 1"))
              if (device == "HPC") {
                    for (n_chains_burnin in c(4, 8, 16)) {
                          add_value(paste0("BayesMVP_HPC_N50000_WCP_c", n_chains_burnin),
                                    integer_text(cell_value(BayesMVP_best, 50000, n_chains_burnin, "best_n_threads_WCP_burnin")),
                                    paste0("HPC: best N_WCP, N = 50,000, ", n_chains_burnin, " chains"))
                    }
                    add_value("BayesMVP_HPC_N50000_c4_sec", seconds_text(cell_value(BayesMVP_best, 50000, 4, "best_sec_per_iter")),
                              "HPC: best seconds per burn-in iteration, N = 50,000, 4 chains")
              }

        }
        ##
        ## ---- tape-chunked Stan (Part VII):
        ##
        for (device in names(x = Stan_best_by_cell_by_device)) {

              Stan_best <-  Stan_best_by_cell_by_device[[device]]
              if (is.null(x = Stan_best) || nrow(x = Stan_best) == 0) next
              add_value(paste0("Stan_", device, "_WCP_min"), range_text(Stan_best$speed_up_vs_best_chunks_WCP_1, "min"), paste0(device, ": min over cells, best vs best N_chunks at N_WCP = 1"))
              add_value(paste0("Stan_", device, "_WCP_max"), range_text(Stan_best$speed_up_vs_best_chunks_WCP_1, "max"), paste0(device, ": max over cells, best vs best N_chunks at N_WCP = 1"))
              add_value(paste0("Stan_", device, "_WCP_selected_count"), integer_text(sum(Stan_best$best_n_threads_WCP_burnin > 1)),
                        paste0(device, ": cells whose best configuration has N_WCP > 1"))
              add_value(paste0("Stan_", device, "_n_cells"), integer_text(nrow(x = Stan_best)), paste0(device, ": number of (N, chains) cells"))
              ##
              ratios <-  Stan_over_BayesMVP(device = device)
              add_value(paste0("BayesMVP_vs_Stan_", device, "_min"), range_text(ratios, "min"), paste0(device, ": min over cells, Stan best / BayesMVP best"))
              add_value(paste0("BayesMVP_vs_Stan_", device, "_max"), range_text(ratios, "max"), paste0(device, ": max over cells, Stan best / BayesMVP best"))
              ##
              if (device == "HPC") {
                    if (any(is.finite(Stan_best$speed_up_vs_best_chunks_WCP_1))) {
                          largest_row <-  Stan_best[which.max(Stan_best$speed_up_vs_best_chunks_WCP_1), , drop = FALSE]
                          add_value("Stan_HPC_WCP_max_N",      paste0("$N = ", fn_paper1_format_number_commas_from_10000(largest_row$N), "$"), "HPC: cell of the largest WCP gain")
                          add_value("Stan_HPC_WCP_max_chains", fn_paper1_burnin_chains_word(n_chains_burnin = largest_row$n_chains_burnin), "HPC: cell of the largest WCP gain")
                    }
                    for (n_chains_burnin in c(4, 8, 16)) {
                          add_value(paste0("Stan_HPC_N50000_chunks_c", n_chains_burnin),
                                    integer_text(cell_value(Stan_best, 50000, n_chains_burnin, "best_num_chunks_burnin")),
                                    paste0("HPC: best N_chunks, N = 50,000, ", n_chains_burnin, " chains"))
                          add_value(paste0("Stan_HPC_N50000_WCP_c", n_chains_burnin),
                                    integer_text(cell_value(Stan_best, 50000, n_chains_burnin, "best_n_threads_WCP_burnin")),
                                    paste0("HPC: best N_WCP, N = 50,000, ", n_chains_burnin, " chains"))
                    }
                    add_value("Stan_HPC_N50000_c4_sec", seconds_text(cell_value(Stan_best, 50000, 4, "best_sec_per_iter")),
                              "HPC: best seconds per burn-in iteration, N = 50,000, 4 chains")
              }

        }
        ##
        if (length(x = values) == 0) return(data.frame(id = character(0), value = character(0), source = character(0)))
        return(do.call(what = rbind, args = values))

}


##
## ---- fn_paper1_burnin_checks: the existing burn-in sentences (Part VI, Discussion) and the direction words of Part VII: ----------------
##
fn_paper1_burnin_checks <-  function( Stan_best_by_cell_by_device,
                                      BayesMVP_best_by_cell_by_device,
                                      expected_configurations_by_implementation_and_device
) {

        checks <-  list()
        add_check <-  function(check, passed, detail) {
              checks[[length(x = checks) + 1]] <<-  data.frame(check = check, status = if (isTRUE(passed)) "OK" else "CHECK TEXT", detail = detail, stringsAsFactors = FALSE)
        }
        BayesMVP_best_by_cell_HPC <-  BayesMVP_best_by_cell_by_device[["HPC"]]
        ##
        if (!is.null(x = BayesMVP_best_by_cell_HPC) && nrow(x = BayesMVP_best_by_cell_HPC) > 0) {
              ##
              ## ---- Part VI: "at N = 10,000 and N = 50,000 on the HPC, the best setting falls from N_threads/chain = 16 with four burn-in chains,
              ##      to 8 with eight chains and 4 with sixteen chains":
              ##
              for (N in c(10000, 50000)) {
                    observed_WCP <-  vapply(X = c(4, 8, 16), FUN.VALUE = numeric(1), FUN = function(n_chains_burnin) {
                          cell_rows <-  BayesMVP_best_by_cell_HPC[BayesMVP_best_by_cell_HPC$N == N & BayesMVP_best_by_cell_HPC$n_chains_burnin == n_chains_burnin, , drop = FALSE]
                          if (nrow(x = cell_rows) == 1) as.numeric(x = cell_rows$best_n_threads_WCP_burnin) else NA_real_
                    })
                    add_check(paste0("Part VI: HPC best N_WCP at N = ", N, " is 16 / 8 / 4 (4 / 8 / 16 chains)"),
                              identical(x = as.numeric(observed_WCP), y = c(16, 8, 4)),
                              paste0("observed: ", paste(observed_WCP, collapse = " / ")))
              }
              ##
              ## ---- Part VI: "on both machines, the gains grow with N, and are largest at N = 50,000" (HPC part):
              ##
              for (speed_up_column in c("speed_up_vs_1_chunk_WCP_1", "speed_up_vs_best_chunks_WCP_1")) {
                    for (n_chains_burnin in sort(x = unique(x = BayesMVP_best_by_cell_HPC$n_chains_burnin))) {
                          chain_rows <-  BayesMVP_best_by_cell_HPC[BayesMVP_best_by_cell_HPC$n_chains_burnin == n_chains_burnin, , drop = FALSE]
                          chain_rows <-  chain_rows[order(chain_rows$N), , drop = FALSE]
                          speed_ups  <-  chain_rows[[speed_up_column]]
                          add_check(paste0("Part VI: HPC ", speed_up_column, " grows with N and is largest at N = 50,000 (", n_chains_burnin, " chains)"),
                                    all(diff(x = speed_ups) > 0) && chain_rows$N[which.max(speed_ups)] == 50000,
                                    paste0("by N: ", paste(formatC(x = speed_ups, digits = 2, format = "f"), collapse = ", ")))
                    }
              }
              ##
              ## ---- Discussion: "WCP was most useful when only a small number of burn-in chains could run, particularly at larger N":
              ##
              for (N in c(10000, 50000)) {
                    N_rows    <-  BayesMVP_best_by_cell_HPC[BayesMVP_best_by_cell_HPC$N == N, , drop = FALSE]
                    N_rows    <-  N_rows[order(N_rows$n_chains_burnin), , drop = FALSE]
                    WCP_gains <-  N_rows$speed_up_vs_best_chunks_WCP_1
                    add_check(paste0("Discussion: HPC WCP gain falls as burn-in chains rise (N = ", N, ")"),
                              length(x = WCP_gains) > 1 && all(diff(x = WCP_gains) < 0),
                              paste0("by chains (", paste(N_rows$n_chains_burnin, collapse = "/"), "): ",
                                     paste(formatC(x = WCP_gains, digits = 2, format = "f"), collapse = ", ")))
              }
        }
        ##
        ## ---- Part VI: "on the laptop, WCP is only selected (with N_threads/chain = 2) for four chains at N >= 2500":
        ##
        BayesMVP_best_by_cell_Laptop <-  BayesMVP_best_by_cell_by_device[["Laptop"]]
        if (!is.null(x = BayesMVP_best_by_cell_Laptop) && nrow(x = BayesMVP_best_by_cell_Laptop) > 0) {
              WCP_cells <-  BayesMVP_best_by_cell_Laptop[BayesMVP_best_by_cell_Laptop$best_n_threads_WCP_burnin > 1, , drop = FALSE]
              expected_cells <-  paste(c(2500, 10000, 50000), 4)
              add_check("Part VI: laptop WCP selected (N_WCP = 2) only for four chains at N >= 2500",
                        setequal(x = paste(WCP_cells$N, WCP_cells$n_chains_burnin), y = expected_cells) && all(WCP_cells$best_n_threads_WCP_burnin == 2),
                        paste0("WCP cells: ", paste0("N = ", WCP_cells$N, ", ", WCP_cells$n_chains_burnin, " chains, N_WCP = ", WCP_cells$best_n_threads_WCP_burnin, collapse = "; ")))
        }
        ##
        ## ---- Part VII: "BayesMVP ... x-fold to y-fold faster per burn-in iteration" needs Stan / BayesMVP > 1 in every cell:
        ##
        for (device in names(x = Stan_best_by_cell_by_device)) {

              Stan_best     <-  Stan_best_by_cell_by_device[[device]]
              BayesMVP_best <-  BayesMVP_best_by_cell_by_device[[device]]
              if (is.null(x = Stan_best) || is.null(x = BayesMVP_best)) next
              both_implementations <-  dplyr::inner_join( x      = dplyr::select(.data = Stan_best,     dplyr::all_of(x = c("N", "n_chains_burnin", "best_sec_per_step"))),
                                                          y      = dplyr::select(.data = BayesMVP_best, dplyr::all_of(x = c("N", "n_chains_burnin", "best_sec_per_step"))),
                                                          by     = c("N", "n_chains_burnin"),
                                                          suffix = c("_Stan", "_BayesMVP"))
              ratios <-  both_implementations$best_sec_per_step_Stan / both_implementations$best_sec_per_step_BayesMVP
              add_check(paste0("Part VII (", device, "): BayesMVP faster than Stan per burn-in iteration in every cell"),
                        length(x = ratios) > 0 && all(ratios > 1),
                        paste0("Stan / BayesMVP: ", paste(formatC(x = ratios, digits = 2, format = "f"), collapse = ", ")))

        }
        ##
        ## ---- selected N_chunks at the largest N_chunks in the runner grid (the optimum may lie beyond the tested range):
        ##
        best_by_cell_lists <-  list(BayesMVP = BayesMVP_best_by_cell_by_device, Stan = Stan_best_by_cell_by_device)
        for (implementation in names(x = best_by_cell_lists)) {
              for (device in names(x = best_by_cell_lists[[implementation]])) {

                    best_by_cell            <-  best_by_cell_lists[[implementation]][[device]]
                    expected_configurations <-  expected_configurations_by_implementation_and_device[[paste0(implementation, "_", device)]]
                    if (is.null(x = best_by_cell) || nrow(x = best_by_cell) == 0 || is.null(x = expected_configurations) || nrow(x = expected_configurations) == 0) next
                    largest_chunks_by_N <-  expected_configurations %>%
                          dplyr::group_by(dplyr::across(.cols = dplyr::all_of(x = "N"))) %>%
                          dplyr::summarise(largest_num_chunks = max(.data$num_chunks_burnin), .groups = "drop")
                    edge_rows <-  dplyr::inner_join(x = best_by_cell, y = largest_chunks_by_N, by = "N")
                    edge_rows <-  edge_rows[edge_rows$best_num_chunks_burnin == edge_rows$largest_num_chunks, , drop = FALSE]
                    add_check(paste0(implementation, " (", device, "): selected N_chunks is not the largest tested N_chunks"),
                              nrow(x = edge_rows) == 0,
                              if (nrow(x = edge_rows) == 0) "none at the edge" else paste0("at the edge: ", paste0("N = ", edge_rows$N, ", ", edge_rows$n_chains_burnin,
                                                                                                                 " chains (", edge_rows$best_num_chunks_burnin, " chunks)", collapse = "; ")))

              }
        }
        ##
        if (length(x = checks) == 0) return(data.frame(check = character(0), status = character(0), detail = character(0)))
        return(do.call(what = rbind, args = checks))

}



##
## ---- fn_paper1_burnin_fill_placeholders: a filled copy of a .tex file (the input file is not changed): ---------------------------------
##
fn_paper1_burnin_fill_placeholders <-  function( tex_path_in,
                                                 tex_path_out,
                                                 placeholder_values
) {

        if (!file.exists(tex_path_in)) stop("fn_paper1_burnin_fill_placeholders: not found: ", tex_path_in)
        if (normalizePath(path = tex_path_in) == normalizePath(path = tex_path_out, mustWork = FALSE)) {
              stop("fn_paper1_burnin_fill_placeholders: tex_path_out must differ from tex_path_in (the draft is never overwritten).")
        }
        ##
        filled_text <-  paste(readLines(con = tex_path_in, warn = FALSE), collapse = "\n")
        for (row_index in seq_len(length.out = nrow(x = placeholder_values))) {

              this_value <-  placeholder_values$value[row_index]
              if (is.na(x = this_value)) next
              filled_text <-  fn_paper1_burnin_replace_fixed( text        = filled_text,
                                                              pattern     = paste0("\\BurninTBD{", placeholder_values$id[row_index], "}"),
                                                              replacement = this_value)

        }
        ##
        remaining_placeholders <-  unique(x = regmatches(x = filled_text, m = gregexpr(pattern = "\\\\BurninTBD\\{[^}]*\\}", text = filled_text))[[1]])
        remaining_placeholders <-  setdiff(x = remaining_placeholders, y = "\\BurninTBD{#1}")
        writeLines(text = filled_text, con = tex_path_out)
        fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_fill_placeholders: wrote ", tex_path_out,
                                               if (length(x = remaining_placeholders) > 0) paste0("; still unfilled: ", paste(remaining_placeholders, collapse = ", ")) else "; all placeholders filled"))
        return(invisible(remaining_placeholders))

}


##
## ---- fn_paper1_burnin_write_Stan_table_tex (research output; not part of the paper): Stan and BayesMVP best configurations per cell: ---
##
fn_paper1_burnin_write_Stan_table_tex <-  function( Stan_best_by_cell,
                                                    BayesMVP_best_by_cell,
                                                    device,
                                                    file_path
) {

        device_full_name <-  c(HPC = "local-HPC", Laptop = "laptop")
        BayesMVP_columns <-  if (is.null(x = BayesMVP_best_by_cell)) NULL else dplyr::select(.data = BayesMVP_best_by_cell,
                                                                                              dplyr::all_of(x = c("N", "n_chains_burnin", "best_num_chunks_burnin",
                                                                                                                  "best_n_threads_WCP_burnin", "best_sec_per_iter")))
        both_implementations <-  if (is.null(x = BayesMVP_columns)) Stan_best_by_cell else dplyr::left_join( x      = Stan_best_by_cell,
                                                                                                             y      = BayesMVP_columns,
                                                                                                             by     = c("N", "n_chains_burnin"),
                                                                                                             suffix = c("", "_BayesMVP"))
        column_or_NA <-  function(row, column_name) if (column_name %in% names(x = row)) row[[column_name]] else NA
        ##
        table_lines <-  c( "\\begin{table}[H]",
                           "\\centering",
                           paste0("\\caption{Burn-in on the ", device_full_name[[device]], ": best measured configuration of tape-chunked Stan and of NicoStan+BayesMVP at each $N$ and burn-in chain count.}"),
                           "\\begin{tabular}{rrrrrrrrr}",
                           "\\hline",
                           " & & \\multicolumn{4}{c}{Tape-chunked Stan} & \\multicolumn{3}{c}{NicoStan+BayesMVP} \\\\",
                           "$N$ & Chains & Chunks & $N_{\\text{WCP}}$ & Sec./iter & WCP & Chunks & $N_{\\text{WCP}}$ & Sec./iter \\\\ \\hline")
        for (row_index in seq_len(length.out = nrow(x = both_implementations))) {

              this_row    <-  both_implementations[row_index, ]
              table_lines <-  c( table_lines,
                                 paste0( fn_paper1_format_number_commas_from_10000(this_row$N),
                                         " & ", this_row$n_chains_burnin,
                                         " & ", this_row$best_num_chunks_burnin,
                                         " & ", this_row$best_n_threads_WCP_burnin,
                                         " & ", fn_paper1_burnin_format_seconds_per_iter(seconds_per_iter_vector = this_row$best_sec_per_iter),
                                         " & ", fn_paper1_burnin_format_speed_up(speed_up_vector = this_row$speed_up_vs_best_chunks_WCP_1),
                                         " & ", column_or_NA(this_row, "best_num_chunks_burnin_BayesMVP"),
                                         " & ", column_or_NA(this_row, "best_n_threads_WCP_burnin_BayesMVP"),
                                         " & ", fn_paper1_burnin_format_seconds_per_iter(seconds_per_iter_vector = as.numeric(x = column_or_NA(this_row, "best_sec_per_iter_BayesMVP"))),
                                         " \\\\"))

        }
        table_lines <-  c(table_lines, "\\hline", "\\end{tabular}", "\\end{table}")
        dir.create(path = dirname(path = file_path), recursive = TRUE, showWarnings = FALSE)
        writeLines(text = table_lines, con = file_path)
        fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_write_Stan_table_tex: wrote ", file_path))
        return(invisible(both_implementations))

}


##
## ---- fn_paper1_burnin_write_Stan_bundle: figures, CSV files, research-output tables, filled draft; optional copy into a manuscript: -----
##
## Figures are written to <output_dir>/Burnin_chunks_WCP/figures with the file names the manuscript uses. If manuscript_dir is not NULL,
## they are also copied to <manuscript_dir>/Files/Supplement/assets/Burnin_chunks_WCP/figures - the folder that Main.tex \includegraphics
## reads (the older bundle writer copies to Files/Generated/Burnin_chunks_WCP, which Main.tex does not read).
##
fn_paper1_burnin_write_Stan_bundle <-  function( BayesMVP_rows_by_device,
                                                 Stan_rows_by_device,
                                                 expected_configurations_by_implementation_and_device,
                                                 missing_configurations,
                                                 output_dir,
                                                 manuscript_dir,
                                                 draft_tex_path
) {

        figure_dir <-  file.path(output_dir, "Burnin_chunks_WCP", "figures")
        values_dir <-  file.path(output_dir, "Burnin_chunks_WCP", "values")
        dir.create(path = figure_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(path = values_dir, recursive = TRUE, showWarnings = FALSE)
        figure_file_paths <-  character(0)
        ##
        ## ---- NicoStan+BayesMVP figures (existing function, unchanged), and best configurations per cell:
        ##
        BayesMVP_best_by_cell_by_device <-  list()
        for (device_selector in names(x = BayesMVP_rows_by_device)) {

              device_rows <-  BayesMVP_rows_by_device[[device_selector]]
              if (is.null(x = device_rows) || nrow(x = device_rows) == 0) next
              figure_file_path <-  file.path(figure_dir, paste0("figure_burnin_", device_selector, ".png"))
              fn_paper1_burnin_figure_best_chunks(rows = device_rows, device = device_selector, file_path = figure_file_path)
              figure_file_paths <-  c(figure_file_paths, figure_file_path)
              BayesMVP_best_by_cell_by_device[[device_selector]] <-  fn_paper1_burnin_best_by_cell( rows           = device_rows,
                                                                                                   device         = device_selector,
                                                                                                   implementation = "BayesMVP")$best_by_cell

        }
        ##
        ## ---- tape-chunked Stan (one figure for all devices with results):
        ##
        Stan_best_by_cell_by_device <-  list()
        for (device_selector in names(x = Stan_rows_by_device)) {

              device_rows <-  Stan_rows_by_device[[device_selector]]
              if (is.null(x = device_rows) || nrow(x = device_rows) == 0) next
              Stan_best_by_cell_by_device[[device_selector]] <-  fn_paper1_burnin_best_by_cell( rows           = device_rows,
                                                                                               device         = device_selector,
                                                                                               implementation = "Stan_tape_chunked")$best_by_cell
              fn_paper1_burnin_write_Stan_table_tex( Stan_best_by_cell     = Stan_best_by_cell_by_device[[device_selector]],
                                                     BayesMVP_best_by_cell = BayesMVP_best_by_cell_by_device[[device_selector]],
                                                     device                = device_selector,
                                                     file_path             = file.path(values_dir, paste0("table_burnin_best_Stan_and_BayesMVP_", device_selector,
                                                                                                          "_research_output.tex")))

        }
        if (length(x = Stan_best_by_cell_by_device) > 0) {
              figure_file_path <-  file.path(figure_dir, "figure_burnin_Stan_tape_chunked.png")
              fn_paper1_burnin_figure_Stan( Stan_rows_by_device             = Stan_rows_by_device,
                                            Stan_best_by_cell_by_device     = Stan_best_by_cell_by_device,
                                            BayesMVP_best_by_cell_by_device = BayesMVP_best_by_cell_by_device,
                                            file_path                       = figure_file_path)
              figure_file_paths <-  c(figure_file_paths, figure_file_path)
              ##
              ## ---- one tape-chunked Stan figure per device as well (the manuscript shows the local-HPC and the laptop separately,
              ##      as for the NicoStan+BayesMVP burn-in figures):
              ##
              for (device_selector in names(x = Stan_best_by_cell_by_device)) {
                    device_figure_file_path <-  file.path(figure_dir, paste0("figure_burnin_Stan_tape_chunked_", device_selector, ".png"))
                    fn_paper1_burnin_figure_Stan( Stan_rows_by_device             = Stan_rows_by_device[device_selector],
                                                  Stan_best_by_cell_by_device     = Stan_best_by_cell_by_device[device_selector],
                                                  BayesMVP_best_by_cell_by_device = BayesMVP_best_by_cell_by_device[device_selector],
                                                  file_path                       = device_figure_file_path)
                    figure_file_paths <-  c(figure_file_paths, device_figure_file_path)
              }
        } else {
              fn_paper1_burnin_message(text = "fn_paper1_burnin_write_Stan_bundle: no tape-chunked Stan rows; the Stan figure, table and values are skipped.")
        }
        ##
        ## ---- per-cell results, placeholder values and checks (CSV):
        ##
        all_best_by_cell <-  dplyr::bind_rows(c(BayesMVP_best_by_cell_by_device, Stan_best_by_cell_by_device))
        utils::write.csv(x = all_best_by_cell, file = file.path(values_dir, "burnin_best_by_cell_BayesMVP_and_Stan.csv"), row.names = FALSE)
        ##
        placeholder_values <-  fn_paper1_burnin_placeholder_values( Stan_best_by_cell_by_device     = Stan_best_by_cell_by_device,
                                                                    BayesMVP_best_by_cell_by_device = BayesMVP_best_by_cell_by_device)
        utils::write.csv(x = placeholder_values, file = file.path(values_dir, "burnin_placeholder_values.csv"), row.names = FALSE)
        print(x = placeholder_values)
        ##
        checks <-  fn_paper1_burnin_checks( Stan_best_by_cell_by_device                          = Stan_best_by_cell_by_device,
                                            BayesMVP_best_by_cell_by_device                      = BayesMVP_best_by_cell_by_device,
                                            expected_configurations_by_implementation_and_device = expected_configurations_by_implementation_and_device)
        utils::write.csv(x = checks, file = file.path(values_dir, "burnin_text_checks.csv"), row.names = FALSE)
        print(x = checks)
        ##
        ##
        ## ---- always written (an empty file = no configuration of the runner grids is missing):
        ##
        if (is.null(x = missing_configurations)) missing_configurations <-  data.frame()
        utils::write.csv(x = missing_configurations, file = file.path(values_dir, "burnin_missing_configurations.csv"), row.names = FALSE)
        ##
        ## ---- filled copy of the draft section:
        ##
        if (!is.null(x = draft_tex_path)) {
              fn_paper1_burnin_fill_placeholders( tex_path_in        = draft_tex_path,
                                                  tex_path_out       = file.path(values_dir, "draft_section_filled.tex"),
                                                  placeholder_values = placeholder_values)
        }
        ##
        ## ---- optional copy of the figures into a manuscript folder (the path Main.tex reads):
        ##
        if (!is.null(x = manuscript_dir)) {

              destination <-  file.path(manuscript_dir, "Files", "Supplement", "assets", "Burnin_chunks_WCP", "figures")
              dir.create(path = destination, recursive = TRUE, showWarnings = FALSE)
              for (figure_file_path in figure_file_paths) {
                    copy_succeeded <-  file.copy(from = figure_file_path, to = file.path(destination, basename(path = figure_file_path)), overwrite = TRUE)
                    if (!copy_succeeded) stop("fn_paper1_burnin_write_Stan_bundle: could not copy ", figure_file_path, " to ", destination)
              }
              fn_paper1_burnin_message(text = paste0("fn_paper1_burnin_write_Stan_bundle: copied ", length(x = figure_file_paths), " figure(s) to ", destination))

        }
        ##
        return(invisible(list( figure_file_paths  = figure_file_paths,
                               best_by_cell       = all_best_by_cell,
                               placeholder_values = placeholder_values,
                               checks             = checks)))

}























