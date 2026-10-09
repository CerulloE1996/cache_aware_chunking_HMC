##
## ==============================================================================================================
## R_fns_alg_paper_1_burnin_report_cmdstanr.R
##
## Paper 1, E1 burn-in, Stan via cmdstanr: functions for the burn-in figures of the Stan via cmdstanr shootout
## (ps_1_burnin_Stan_cmdstanr_shootout_4_chains.R), one figure per device, in the visual design of the
## NicoStan+BayesMVP burn-in figures (fn_paper1_burnin_figure_best_chunks() in R_fns_alg_paper_1_burnin_report.R)
## and of the Stan via NicoStan burn-in figures.
##
## This file only defines functions; alg_paper_1_burnin_report_cmdstanr.R is the runner. Nothing here runs
## sampling or compiles anything: the shootout runner is only parsed, for its settings, and its saved .rds files
## are only read.
##
## R_fns_alg_paper_1_burnin_report.R must be sourced first, for its number and fold-label helpers
## (fn_paper1_format_number_commas_from_10000() and fn_paper1_burnin_format_fold()).
##
## Column schema of the saved .rds files (one per device and N;
## see ps_1_burnin_Stan_cmdstanr_shootout_4_chains.R):
##   device ("HPC" or "Laptop"), N, algorithm, n_chains_burnin, n_threads_per_chain, n_threads_total,
##   num_chunks (NA for AD_Stan), chunk_size, run_number, n_iter, max_treedepth, mean_n_leapfrog,
##   time_wall_seconds, time_sampling_slowest_chain_seconds, time_sampling_mean_chain_seconds,
##   sec_per_step_slowest_chain, n_divergences.
##
## The figures and values use seconds per burn-in iteration of the slowest chain, i.e., per run,
## time_sampling_slowest_chain_seconds / n_iter (every configuration ran 2^max_treedepth - 1 = 15 leapfrog steps
## per iteration; mean_n_leapfrog 14.94-15), averaged over the runs of each configuration.
## ==============================================================================================================
##
##
##
## ---- fn_paper1_burnin_cmdstanr_runner_settings: the shootout runner's current settings: -----------------------
##
## Evaluates only the top-level settings assignments of the shootout runner (N grid, burn-in chain count, runs,
## maximum tree depth, timed iterations, N_chunks grid, N_threads/chain grid, output folder, and its file name
## and configuration grid functions) for one device, so that the report always follows the runner's current
## settings; the data simulation, model compilation and sampling steps of the runner are never evaluated.
##
fn_paper1_burnin_cmdstanr_runner_settings <-  function( runner_file_path,
                                                        device
) {

        if (!is.character(runner_file_path) || length(runner_file_path) != 1 || !file.exists(runner_file_path)) {
              stop(paste0( "fn_paper1_burnin_cmdstanr_runner_settings: runner_file_path must be an existing ",
                           "file, got '", paste0(runner_file_path, collapse = ", "), "'."))
        }
        if (!is.character(device) || length(device) != 1 || !(device %in% c("HPC", "Laptop"))) {
              stop("fn_paper1_burnin_cmdstanr_runner_settings: device must be 'HPC' or 'Laptop'.")
        }
        ##
        ## ---- the runner's settings assignments (everything else in the runner is skipped):
        ##
        settings_names <-  c( "algorithm_study_dir", "N_vec", "n_chains_burnin", "n_runs", "step_size",
                              "max_treedepth", "n_iter_given_N", "chunk_grid_given_N", "threads_per_chain_vec",
                              "output_dir", "fn_file_name", "fn_configuration_grid")
        ##
        settings_environment <-  new.env(parent = baseenv())
        assign(x = "device", value = device, envir = settings_environment)
        ##
        for (runner_expression in as.list(parse(file = runner_file_path, keep.source = FALSE))) {

              is_settings_assignment <-  is.call(runner_expression) &&
                                         identical(runner_expression[[1]], as.name("<-")) &&
                                         is.name(runner_expression[[2]]) &&
                                         as.character(runner_expression[[2]]) %in% settings_names
              ##
              if (is_settings_assignment) eval(expr = runner_expression, envir = settings_environment)

        }
        ##
        settings_found <-  vapply( settings_names,
                                   function(setting_name) exists( setting_name,
                                                                  envir    = settings_environment,
                                                                  inherits = FALSE),
                                   logical(1))
        if (!all(settings_found)) {
              stop(paste0( "fn_paper1_burnin_cmdstanr_runner_settings: setting(s) not found in ",
                           basename(runner_file_path), ": ",
                           paste0(settings_names[!settings_found], collapse = ", ")))
        }
        ##
        runner_settings <-  mget(x = settings_names, envir = settings_environment)
        runner_settings$device           <-  device
        runner_settings$runner_file_path <-  runner_file_path
        ##
        ## ---- the runner's own file name for each N:
        ##
        runner_settings$file_path_given_N <-  vapply( runner_settings$N_vec,
                                                      function(N) file.path( runner_settings$output_dir,
                                                                             runner_settings$fn_file_name(N)),
                                                      character(1))
        names(runner_settings$file_path_given_N) <-  as.character(runner_settings$N_vec)
        ##
        ## ---- the runner's own configuration grid (algorithm, N_threads/chain, N_chunks) for each N:
        ##
        runner_settings$configuration_grid <-  do.call(rbind, lapply(runner_settings$N_vec, function(N) {

              configurations <-  runner_settings$fn_configuration_grid(N)
              ##
              data.frame( N                   = N,
                          algorithm           = vapply( configurations,
                                                        function(configuration) configuration$algorithm,
                                                        character(1)),
                          n_threads_per_chain = vapply( configurations,
                                                        function(configuration) configuration$threads_per_chain,
                                                        numeric(1)),
                          num_chunks          = vapply( configurations,
                                                        function(configuration) {
                                                              as.numeric(configuration$num_chunks)
                                                        },
                                                        numeric(1)))

        }))
        ##
        return(runner_settings)

}
##
##
##
## ---- fn_paper1_burnin_cmdstanr_read_runs: saved runs of one device matching the runner's current settings: ----
##
## Keeps only the runs whose burn-in N_chains, maximum tree depth, timed iterations, run number and configuration
## (algorithm, N_threads/chain, N_chunks) match the runner's current settings and grid; every left-out run and
## every configuration with missing runs is reported, and the saved .rds files are never changed.
##
fn_paper1_burnin_cmdstanr_read_runs <-  function( runner_settings ) {

        ## was: required_columns <-  c( "device", "N", "algorithm", "n_chains_burnin", "n_threads_per_chain",
        ## was:                         "num_chunks", "run_number", "n_iter", "max_treedepth", "mean_n_leapfrog",
        ## was:                         "sec_per_step_slowest_chain", "n_divergences")
        required_columns <-  c( "device", "N", "algorithm", "n_chains_burnin", "n_threads_per_chain",
                                "num_chunks", "run_number", "n_iter", "max_treedepth", "mean_n_leapfrog",
                                "time_sampling_slowest_chain_seconds", "n_divergences")
        ##
        message_start <-  paste0("Stan via cmdstanr burn-in, ", runner_settings$device)
        ##
        runs_list <-  list()
        ##
        for (N in runner_settings$N_vec) {

              file_path <-  runner_settings$file_path_given_N[[as.character(N)]]
              ##
              if (!file.exists(file_path)) {
                    message(NicoStan::colourise( text = paste0( message_start, ", N = ", N, ": no saved runs (",
                                                                basename(file_path),
                                                                "); this panel is left out."),
                                                 fg   = "red"))
                    next
              }
              ##
              runs <-  readRDS(file = file_path)
              ##
              missing_columns <-  setdiff(required_columns, names(runs))
              if (length(missing_columns) > 0) {
                    stop(paste0( "fn_paper1_burnin_cmdstanr_read_runs: ", basename(file_path),
                                 " has no column(s) ", paste0(missing_columns, collapse = ", "), "."))
              }
              ##
              ## ---- runs that match the runner's current settings and configuration grid:
              ##
              grid_for_N <-  runner_settings$configuration_grid[runner_settings$configuration_grid$N == N, ]
              grid_keys  <-  paste( grid_for_N$algorithm, grid_for_N$n_threads_per_chain, grid_for_N$num_chunks,
                                    sep = "_")
              run_keys   <-  paste(runs$algorithm, runs$n_threads_per_chain, runs$num_chunks, sep = "_")
              ##
              matches_runner_settings <-  runs$device == runner_settings$device &
                                          runs$N == N &
                                          runs$n_chains_burnin == runner_settings$n_chains_burnin &
                                          runs$max_treedepth == runner_settings$max_treedepth &
                                          runs$n_iter == runner_settings$n_iter_given_N[[as.character(N)]] &
                                          runs$run_number %in% seq_len(runner_settings$n_runs) &
                                          run_keys %in% grid_keys
              matches_runner_settings[is.na(matches_runner_settings)] <-  FALSE
              ##
              if (any(!matches_runner_settings)) {
                    message(NicoStan::colourise( text = paste0( message_start, ", N = ", N, ": left out ",
                                                                sum(!matches_runner_settings), " run(s) whose ",
                                                                "settings differ from the runner's."),
                                                 fg   = "cyan"))
              }
              ##
              ## ---- every configuration of the runner's grid should have all of its runs:
              ##
              runs_per_configuration    <-  table(factor(run_keys[matches_runner_settings], levels = grid_keys))
              incomplete_configurations <-  names(runs_per_configuration)[ runs_per_configuration <
                                                                               runner_settings$n_runs]
              if (length(incomplete_configurations) > 0) {
                    message(NicoStan::colourise( text = paste0( message_start, ", N = ", N, ": fewer than ",
                                                                runner_settings$n_runs, " runs for ",
                                                                paste0( incomplete_configurations,
                                                                        collapse = ", ")),
                                                 fg   = "red"))
              }
              ##
              ## was: runs_list[[length(runs_list) + 1]] <-  runs[matches_runner_settings, required_columns,
              ## was:                                             drop = FALSE]
              runs_matched <-  runs[matches_runner_settings, required_columns, drop = FALSE]
              ##
              ## ---- seconds per burn-in iteration of the slowest chain, per run:
              ##
              runs_matched$sec_per_burnin_iteration_slowest_chain <-
                    runs_matched$time_sampling_slowest_chain_seconds / runs_matched$n_iter
              ##
              runs_list[[length(runs_list) + 1]] <-  runs_matched

        }
        ##
        if (length(runs_list) == 0) {
              stop(paste0( "fn_paper1_burnin_cmdstanr_read_runs: no saved runs for device '",
                           runner_settings$device, "' in ", runner_settings$output_dir, "."))
        }
        ##
        runs_all <-  do.call(rbind, runs_list)
        ##
        message(NicoStan::colourise( text = paste0( message_start, ": read ", nrow(runs_all), " run(s) for N = ",
                                                    paste0(sort(unique(runs_all$N)), collapse = ", "), "."),
                                     fg   = "cyan"))
        ##
        return(runs_all)

}
##
##
##
## ---- fn_paper1_burnin_cmdstanr_summarise_runs: mean over runs per configuration: ------------------------------
##
## was: ## Mean seconds per leapfrog step of the slowest chain over the runs of each (device, N, algorithm,
## was: ## N_threads/chain, N_chunks), as in ps_1_burnin_Stan_cmdstanr_shootout_summary.R.
## Mean seconds per burn-in iteration of the slowest chain over the runs of each (device, N, algorithm,
## N_threads/chain, N_chunks), as in ps_1_burnin_Stan_cmdstanr_shootout_summary.R.
##
fn_paper1_burnin_cmdstanr_summarise_runs <-  function( runs ) {

        grouping_columns <-  c("device", "N", "algorithm", "n_threads_per_chain", "num_chunks")
        ##
        configuration_summary <-  dplyr::group_by( runs,
                                                   dplyr::across(.cols = dplyr::all_of(grouping_columns)))
        ## was: configuration_summary <-  dplyr::summarise( configuration_summary,
        ## was:                                             sec_per_step     =
        ## was:                                                 mean(.data$sec_per_step_slowest_chain),
        ## was:                                             n_runs_available = dplyr::n(),
        ## was:                                             .groups          = "drop")
        configuration_summary <-  dplyr::summarise(
                                      configuration_summary,
                                      sec_per_burnin_iteration_slowest_chain =
                                          mean(.data$sec_per_burnin_iteration_slowest_chain),
                                      n_runs_available                       = dplyr::n(),
                                      .groups                                = "drop")
        ##
        return(as.data.frame(configuration_summary))

}
##
##
##
## ---- fn_paper1_burnin_cmdstanr_figure_values: per device and N, the values shown in each figure panel: --------
##
## was: ##   - the best measured configuration (lowest mean seconds per leapfrog step of the slowest chain, over
## was: ##     every configuration including AD_Stan) and every runner-up within 1% of it;
## was: ##   - the best tape chunking only configuration (AD_Stan_tape_chunked, N_threads/chain = 1) and the "WCP"
## was: ##     fold, i.e., its seconds per leapfrog step divided by that of the best measured configuration;
## was: ##   - AD_Stan (N_threads/chain = 1) and the "total" fold, i.e., its seconds per leapfrog step divided by
## was: ##     that of the best measured configuration, and the speed-up of tape chunking only relative to
## was: ##     AD_Stan.
##   - the best measured configuration (lowest mean seconds per burn-in iteration of the slowest chain, over
##     every configuration including AD_Stan) and every runner-up within 1% of it;
##   - the best tape chunking only configuration (AD_Stan_tape_chunked, N_threads/chain = 1) and the "WCP" fold,
##     i.e., its seconds per burn-in iteration divided by that of the best measured configuration;
##   - AD_Stan (N_threads/chain = 1) and the "total" fold, i.e., its seconds per burn-in iteration divided by
##     that of the best measured configuration, and the speed-up of tape chunking only relative to AD_Stan.
##
fn_paper1_burnin_cmdstanr_figure_values <-  function( configuration_summary ) {

        fn_first_or_NA <-  function(values) if (length(values) == 0) NA else values[1]
        ##
        values_rows <-  list()
        ##
        for (device in unique(configuration_summary$device)) {
              for (N in sort(unique(configuration_summary$N[configuration_summary$device == device]))) {

                    cell <-  configuration_summary[ configuration_summary$device == device &
                                                    configuration_summary$N == N, , drop = FALSE]
                    ## was: cell <-  cell[order(cell$sec_per_step), , drop = FALSE]
                    cell <-  cell[order(cell$sec_per_burnin_iteration_slowest_chain), , drop = FALSE]
                    ##
                    best_configuration <-  cell[1, , drop = FALSE]
                    ## was: best_sec_per_step  <-  best_configuration$sec_per_step
                    best_sec_per_burnin_iteration_slowest_chain <-
                          best_configuration$sec_per_burnin_iteration_slowest_chain
                    ##
                    AD_Stan_rows      <-  cell[cell$algorithm == "AD_Stan", , drop = FALSE]
                    tape_chunked_rows <-  cell[ cell$algorithm == "AD_Stan_tape_chunked" &
                                                cell$n_threads_per_chain == 1, , drop = FALSE]
                    ##
                    ## was: AD_Stan_sec_per_step           <-  fn_first_or_NA(AD_Stan_rows$sec_per_step)
                    ## was: best_tape_chunked_sec_per_step <-  fn_first_or_NA(tape_chunked_rows$sec_per_step)
                    AD_Stan_sec_per_burnin_iteration_slowest_chain           <-
                          fn_first_or_NA(AD_Stan_rows$sec_per_burnin_iteration_slowest_chain)
                    best_tape_chunked_sec_per_burnin_iteration_slowest_chain <-
                          fn_first_or_NA(tape_chunked_rows$sec_per_burnin_iteration_slowest_chain)
                    ##
                    ## ---- runners-up within 1% of the best measured configuration (essentially joint-best):
                    ##
                    runners_up <-  cell[-1, , drop = FALSE]
                    ## was: runners_up <-  runners_up[runners_up$sec_per_step <= 1.01 * best_sec_per_step, ,
                    ## was:                             drop = FALSE]
                    runners_up <-  runners_up[ runners_up$sec_per_burnin_iteration_slowest_chain <=
                                                   1.01 * best_sec_per_burnin_iteration_slowest_chain, ,
                                               drop = FALSE]
                    ##
                    within_1_percent_of_best <-  "none"
                    if (nrow(runners_up) > 0) {
                          runners_up_percent_slower_than_best <-
                                100 * ( runners_up$sec_per_burnin_iteration_slowest_chain /
                                            best_sec_per_burnin_iteration_slowest_chain - 1)
                          within_1_percent_of_best <-  paste0( runners_up$algorithm, " ",
                                                               runners_up$n_threads_per_chain, "/",
                                                               runners_up$num_chunks, " (+",
                                                               ## was: formatC( 100 * (runners_up$sec_per_step /
                                                               ## was:                     best_sec_per_step - 1),
                                                               formatC( runners_up_percent_slower_than_best,
                                                                        format = "f", digits = 2),
                                                               "%)",
                                                               collapse = "; ")
                    }
                    ##
                    ## ---- was (seconds per leapfrog step columns):
                    ##        best_sec_per_step                = best_sec_per_step,
                    ##        best_tape_chunked_sec_per_step   = best_tape_chunked_sec_per_step,
                    ##        WCP_fold                         = best_tape_chunked_sec_per_step /
                    ##                                               best_sec_per_step,
                    ##        AD_Stan_sec_per_step             = AD_Stan_sec_per_step,
                    ##        total_fold                       = AD_Stan_sec_per_step / best_sec_per_step,
                    ##        tape_chunked_speed_up_vs_AD_Stan = AD_Stan_sec_per_step /
                    ##                                               best_tape_chunked_sec_per_step,
                    ##
                    WCP_fold   <-  best_tape_chunked_sec_per_burnin_iteration_slowest_chain /
                                       best_sec_per_burnin_iteration_slowest_chain
                    total_fold <-  AD_Stan_sec_per_burnin_iteration_slowest_chain /
                                       best_sec_per_burnin_iteration_slowest_chain
                    tape_chunked_speed_up_vs_AD_Stan <-
                          AD_Stan_sec_per_burnin_iteration_slowest_chain /
                              best_tape_chunked_sec_per_burnin_iteration_slowest_chain
                    ##
                    values_rows[[length(values_rows) + 1]] <-  data.frame(
                          device                           = device,
                          N                                = N,
                          best_algorithm                   = best_configuration$algorithm,
                          best_n_threads_per_chain         = best_configuration$n_threads_per_chain,
                          best_num_chunks                  = best_configuration$num_chunks,
                          best_sec_per_burnin_iteration_slowest_chain              =
                              best_sec_per_burnin_iteration_slowest_chain,
                          best_tape_chunked_num_chunks     = fn_first_or_NA(tape_chunked_rows$num_chunks),
                          best_tape_chunked_sec_per_burnin_iteration_slowest_chain =
                              best_tape_chunked_sec_per_burnin_iteration_slowest_chain,
                          WCP_fold                         = WCP_fold,
                          AD_Stan_sec_per_burnin_iteration_slowest_chain           =
                              AD_Stan_sec_per_burnin_iteration_slowest_chain,
                          total_fold                       = total_fold,
                          tape_chunked_speed_up_vs_AD_Stan = tape_chunked_speed_up_vs_AD_Stan,
                          within_1_percent_of_best         = within_1_percent_of_best)

              }
        }
        ##
        return(do.call(rbind, values_rows))

}
##
##
##
## ---- fn_paper1_burnin_cmdstanr_figure_panel (private helper): one N panel of the figure: ---------------------
##
## Same look as fn_paper1_burnin_figure_best_chunks_panel() in R_fns_alg_paper_1_burnin_report.R:
## x = N_threads/chain on a log2 scale with a tick at every measured value (two rows of tick labels where adjacent
## was: ## values are close), y = seconds per leapfrog step of the slowest chain (log10), with this panel's own
## was: ## range.
## values are close), y = seconds per burn-in iteration of the slowest chain (log10), with this panel's own range.
## AD_Stan is an open square with a dashed line at its level, the best N_chunks at each N_threads/chain is a solid
## line (points labelled with N_chunks, leaving out any label that would overlap a label drawn before it),
## WCP-only (AD_Stan_WCP) is a green dashed line with triangles, and the best measured configuration is circled.
##
fn_paper1_burnin_cmdstanr_figure_panel <-  function( panel_best_by_threads,
                                                     panel_no_chunking,
                                                     panel_WCP_only,
                                                     panel_best_overall,
                                                     panel_fold_label,
                                                     show_right_strips,
                                                     series_levels,
                                                     series_labels,
                                                     best_overall_label,
                                                     y_axis_title
) {

        ##
        ## ---- x breaks at every N_threads/chain measured in this panel:
        ##
        measured_threads_values <-  sort(unique(c( panel_best_by_threads$n_threads_per_chain,
                                                   panel_no_chunking$n_threads_per_chain,
                                                   panel_WCP_only$n_threads_per_chain)))
        ##
        ## ---- two rows of x labels where adjacent measured values are close on the log2 axis (e.g. 16 and 22):
        ##
        measured_threads_log2 <-  log2(measured_threads_values)
        smallest_relative_gap <-  if (length(measured_threads_values) > 1) {
                                      min(diff(measured_threads_log2)) / diff(range(measured_threads_log2))
                                  } else 1
        x_axis_label_rows     <-  if (smallest_relative_gap < 0.10) 2 else 1
        ##
        ## ---- log10 y breaks, thinned from the top down so that no two labels are closer than 9% of the panel's
        ##      log10 range (as in fn_paper1_burnin_figure_best_chunks_panel()):
        ##
        y_breaks_function <-  function(y_limits) {

              candidate_breaks <-  scales::breaks_log(n = 6)(y_limits)
              candidate_breaks <-  sort( candidate_breaks[ candidate_breaks >= min(y_limits) &
                                                           candidate_breaks <= max(y_limits)],
                                         decreasing = TRUE)
              if (length(candidate_breaks) < 2) return(candidate_breaks)
              minimum_log10_gap <-  0.09 * diff(log10(range(y_limits)))
              kept_breaks <-  candidate_breaks[1]
              for (candidate_break in candidate_breaks[-1]) {
                    if (log10(kept_breaks[length(kept_breaks)]) - log10(candidate_break) >= minimum_log10_gap) {
                          kept_breaks <-  c(kept_breaks, candidate_break)
                    }
              }
              return(sort(kept_breaks))

        }
        ##
        ## ---- y labels in plain decimals (e.g. 0.0006 rather than 6e-04), without trailing zeros:
        ##
        y_labels_function <-  function(y_breaks) {

              y_labels <-  format(y_breaks, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
              y_labels[is.na(y_breaks)] <-  ""
              return(y_labels)

        }
        ##
        ## ---- N_chunks labels sit below each point, except above a point that is slower than its left-hand
        ##      neighbour (so that the label does not sit on the rising line segment) or slower than AD_Stan at
        ##      the same N_threads/chain (so that the label does not sit on the AD_Stan square below it); the
        ##      label of the circled best point sits further below, clear of its circle:
        ##
        panel_best_by_threads <-  panel_best_by_threads[ order(panel_best_by_threads$n_threads_per_chain), ,
                                                         drop = FALSE]
        ## was: rises_from_left <-  c(FALSE, diff(panel_best_by_threads$sec_per_step) > 0)
        rises_from_left <-  c(FALSE, diff(panel_best_by_threads$sec_per_burnin_iteration_slowest_chain) > 0)
        ##
        no_chunking_index    <-  match( panel_best_by_threads$n_threads_per_chain,
                                        panel_no_chunking$n_threads_per_chain)
        ## was: is_above_no_chunking <-  panel_best_by_threads$sec_per_step >
        ## was:                              panel_no_chunking$sec_per_step[no_chunking_index]
        is_above_no_chunking <-  panel_best_by_threads$sec_per_burnin_iteration_slowest_chain >
                                     panel_no_chunking$sec_per_burnin_iteration_slowest_chain[no_chunking_index]
        is_above_no_chunking[is.na(is_above_no_chunking)] <-  FALSE
        ##
        is_best_overall <-  panel_best_by_threads$algorithm  == panel_best_overall$algorithm[1] &
                            panel_best_by_threads$num_chunks == panel_best_overall$num_chunks[1] &
                            panel_best_by_threads$n_threads_per_chain == panel_best_overall$n_threads_per_chain[1]
        is_best_overall[is.na(is_best_overall)] <-  FALSE
        ##
        panel_best_by_threads$label_vjust <-  ifelse( rises_from_left | is_above_no_chunking, -1.0,
                                                      ifelse(is_best_overall, 3.0, 2.0))
        ##
        ## ---- labels in drawing order: the circled best point first, then from the smallest N_threads/chain
        ##      upwards; check_overlap then leaves out any label that would overlap a label drawn before it:
        ##
        label_order  <-  order(!is_best_overall, panel_best_by_threads$n_threads_per_chain)
        panel_labels <-  panel_best_by_threads[label_order, , drop = FALSE]
        ##
        ## ---- with options(paper1.Stan_N_chunks_as_partial_sums = TRUE), the labels give the number of partial sums which
        ##      reduce_sum_static() actually ran (R_fn_number_of_partial_sums_run_by_reduce_sum_static.R), not the requested N_chunks:
        if (isTRUE(x = getOption(x = "paper1.Stan_N_chunks_as_partial_sums")) && nrow(x = panel_labels) > 0) {
              panel_labels$num_chunks <-  fn_number_of_partial_sums_run_by_reduce_sum_static( N_units            = panel_labels$N,
                                                                              N_chunks_requested = panel_labels$num_chunks)
        }
        ##
        ## ---- series styles (no chunking = AD_Stan, best N_chunks at each N_threads/chain,
        ##      WCP-only = AD_Stan_WCP), as in fn_paper1_burnin_figure_best_chunks_panel():
        ##
        series_colours   <-  c("#D55E00", "#0072B2", "#009E73")
        series_linetypes <-  c("dashed",  "solid",   "22")
        series_shapes    <-  c(0,         16,        17)
        names(series_colours)   <-  series_levels
        names(series_linetypes) <-  series_levels
        names(series_shapes)    <-  series_levels
        ##
        best_overall_fill <-  c(NA)
        names(best_overall_fill) <-  best_overall_label
        ##
        ## ---- every y mapping below was .data$sec_per_step (seconds per leapfrog step of the slowest chain)
        ##      and is now .data$sec_per_burnin_iteration_slowest_chain:
        ##
        panel_plot <-  ggplot2::ggplot() +
              ##
              ## ---- dashed reference level of AD_Stan (N_chunks = 1, N_threads/chain = 1) across the panel:
              ##
              ggplot2::geom_hline( data      = panel_no_chunking,
                                   mapping   = ggplot2::aes( yintercept =
                                                                 .data$sec_per_burnin_iteration_slowest_chain,
                                                             colour     = .data$series,
                                                             linetype   = .data$series),
                                   linewidth = 0.6) +
              ##
              ## ---- best N_chunks at each N_threads/chain:
              ##
              ggplot2::geom_line( data      = panel_best_by_threads,
                                  mapping   = ggplot2::aes( x        = .data$n_threads_per_chain,
                                                            y        =
                                                                .data$sec_per_burnin_iteration_slowest_chain,
                                                            colour   = .data$series,
                                                            linetype = .data$series),
                                  linewidth = 0.8) +
              ggplot2::geom_point( data    = panel_best_by_threads,
                                   mapping = ggplot2::aes( x      = .data$n_threads_per_chain,
                                                           y      = .data$sec_per_burnin_iteration_slowest_chain,
                                                           colour = .data$series,
                                                           shape  = .data$series),
                                   size    = 2.2) +
              ggplot2::geom_text( data          = panel_labels,
                                  mapping       = ggplot2::aes( x     = .data$n_threads_per_chain,
                                                                y     =
                                                                    .data$sec_per_burnin_iteration_slowest_chain,
                                                                label = .data$num_chunks,
                                                                vjust = .data$label_vjust),
                                  colour        = series_colours[[2]],
                                  size          = 3.3,
                                  check_overlap = TRUE,
                                  show.legend   = FALSE) +
              ##
              ## ---- WCP-only (AD_Stan_WCP, i.e., N_chunks = N_threads/chain):
              ##
              ggplot2::geom_line( data      = panel_WCP_only,
                                  mapping   = ggplot2::aes( x        = .data$n_threads_per_chain,
                                                            y        =
                                                                .data$sec_per_burnin_iteration_slowest_chain,
                                                            colour   = .data$series,
                                                            linetype = .data$series),
                                  linewidth = 0.8) +
              ggplot2::geom_point( data    = panel_WCP_only,
                                   mapping = ggplot2::aes( x      = .data$n_threads_per_chain,
                                                           y      = .data$sec_per_burnin_iteration_slowest_chain,
                                                           colour = .data$series,
                                                           shape  = .data$series),
                                   size    = 2.6) +
              ##
              ## ---- AD_Stan measurement (drawn after the best line, so it stays visible where the two coincide):
              ##
              ggplot2::geom_point( data    = panel_no_chunking,
                                   mapping = ggplot2::aes( x      = .data$n_threads_per_chain,
                                                           y      = .data$sec_per_burnin_iteration_slowest_chain,
                                                           colour = .data$series,
                                                           shape  = .data$series),
                                   size    = 3.4,
                                   stroke  = 1.1) +
              ##
              ## ---- best measured configuration, circled:
              ##
              ggplot2::geom_point( data    = panel_best_overall,
                                   mapping = ggplot2::aes( x    = .data$n_threads_per_chain,
                                                           y    = .data$sec_per_burnin_iteration_slowest_chain,
                                                           fill = .data$best_overall_series),
                                   shape   = 21,
                                   size    = 5.0,
                                   stroke  = 1.0,
                                   colour  = "black") +
              ##
              ## ---- fold labels of the best measured configuration ("total" relative to AD_Stan, "WCP" relative
              ##      to the best tape chunking only configuration):
              ##
              ggplot2::annotate( geom  = "text",
                                 x     = Inf,
                                 y     = Inf,
                                 label = panel_fold_label,
                                 hjust = 1.06,
                                 vjust = 1.25,
                                 size  = 3.4) +
              ggplot2::scale_x_continuous( transform    = "log2",
                                           breaks       = measured_threads_values,
                                           labels       = measured_threads_values,
                                           minor_breaks = NULL,
                                           ## ---- 0.14 (0.08 in fn_paper1_burnin_figure_best_chunks_panel()), so
                                           ##      that 3- and 4-digit N_chunks labels at the smallest and
                                           ##      largest N_threads/chain stay inside the panel:
                                           expand       = ggplot2::expansion(mult = 0.14),
                                           guide        = ggplot2::guide_axis(n.dodge = x_axis_label_rows)) +
              ##
              ## ---- y expansion: 0.30 below (0.20 in fn_paper1_burnin_figure_best_chunks_panel()), so that the
              ##      N_chunks label under the circled best point stays inside the panel; 0.45 above, for the
              ##      fold labels:
              ##
              ggplot2::scale_y_log10( breaks = y_breaks_function,
                                      labels = y_labels_function,
                                      expand = ggplot2::expansion(mult = c(0.30, 0.45)),
                                      guide  = ggplot2::guide_axis(check.overlap = TRUE)) +
              ggplot2::scale_colour_manual(   name   = NULL,
                                              values = series_colours,
                                              breaks = series_levels,
                                              labels = series_labels) +
              ggplot2::scale_linetype_manual( name   = NULL,
                                              values = series_linetypes,
                                              breaks = series_levels,
                                              labels = series_labels) +
              ggplot2::scale_shape_manual(    name   = NULL,
                                              values = series_shapes,
                                              breaks = series_levels,
                                              labels = series_labels) +
              ggplot2::scale_fill_manual(     name     = NULL,
                                              values   = best_overall_fill,
                                              na.value = NA) +
              ggplot2::guides( colour   = ggplot2::guide_legend(order = 1),
                               linetype = ggplot2::guide_legend(order = 1),
                               shape    = ggplot2::guide_legend(order = 1),
                               fill     = ggplot2::guide_legend(order = 2)) +
              ggplot2::labs( x = expression(N["threads/chain"] ~ "(threads of each burn-in chain)"),
                             y = y_axis_title) +
              ggplot2::theme_bw(base_size = 14) +
              ggplot2::theme( legend.position  = "bottom",
                              legend.direction = "vertical",
                              legend.key.width = ggplot2::unit(x = 1.6, units = "lines"))
        ##
        ## ---- N strip on top of every panel; the burn-in N_chains strip (inner) and the device strip (outer),
        ##      both plotmath and so parsed, on the right of the last panel only:
        ##
        if (show_right_strips) {
              panel_plot <-  panel_plot +
                    ggplot2::facet_grid( rows     = ggplot2::vars(.data$device_label, .data$chains_label),
                                         cols     = ggplot2::vars(.data$N_label),
                                         labeller = ggplot2::labeller(.rows = ggplot2::label_parsed))
        } else {
              panel_plot <-  panel_plot +
                    ggplot2::facet_grid(cols = ggplot2::vars(.data$N_label))
        }
        ##
        return(panel_plot)

}
##
##
##
## was: ## ---- fn_paper1_burnin_cmdstanr_figure: seconds per leapfrog step against N_threads/chain, one panel per
## was: ##      N: ----
## ---- fn_paper1_burnin_cmdstanr_figure: seconds per burn-in iteration against N_threads/chain, per N: ---------
##
## One row of panels (the shootout fixes the burn-in N_chains), one panel per N, each with its own x and y range;
## the fold labels come from fn_paper1_burnin_cmdstanr_figure_values(), so the figure and the values file agree.
## Each panel is its own ggplot and patchwork assembles the row, as in fn_paper1_burnin_figure_best_chunks().
##
fn_paper1_burnin_cmdstanr_figure <-  function( configuration_summary,
                                               figure_values,
                                               device,
                                               n_chains_burnin,
                                               file_path,
                                               width_inches,
                                               height_inches,
                                               dpi,
                                               show_total_fold_label
) {

        for (helper_name in c("fn_paper1_format_number_commas_from_10000", "fn_paper1_burnin_format_fold")) {
              if (!exists(helper_name, mode = "function")) {
                    stop(paste0( "fn_paper1_burnin_cmdstanr_figure: ", helper_name, "() not found; source ",
                                 "R_fns_alg_paper_1_burnin_report.R first."))
              }
        }
        ##
        device_rows <-  configuration_summary[configuration_summary$device == device, , drop = FALSE]
        if (nrow(device_rows) == 0) {
              stop(paste0("fn_paper1_burnin_cmdstanr_figure: no configurations for device '", device, "'."))
        }
        ##
        N_values <-  sort(unique(device_rows$N))
        ##
        ## ---- strip labels in the paper's notation (the burn-in N_chains and device strips are plotmath):
        ##
        N_label_levels <-  paste0("N = ", fn_paper1_format_number_commas_from_10000(N_values))
        ##
        device_rows$N_label      <-  factor( paste0( "N = ",
                                                     fn_paper1_format_number_commas_from_10000(device_rows$N)),
                                             levels = N_label_levels)
        device_rows$chains_label <-  paste0('N["burn_chains"]==', n_chains_burnin)
        device_rows$device_label <-  if (device == "HPC") '"local-HPC"' else '"laptop"'
        ##
        series_levels <-  c("no_chunking", "best_chunks_at_each_threads_per_chain", "WCP_only")
        series_labels <-  c( expression("AD_Stan (" * N[chunks] * " = 1, " * N["threads/chain"] * " = 1)"),
                             expression("Best " * N[chunks] * " at each " * N["threads/chain"] *
                                        " (point labels give " * N[chunks] * ")"),
                             # expression("AD_Stan_WCP (WCP-only: " * N[chunks] * " = " * N["threads/chain"] * ")"))
                             if (isTRUE(x = getOption(x = "paper1.Stan_N_chunks_as_partial_sums"))) {
                                 expression("AD_Stan_WCP (WCP-only: requested " * N[chunks] * " = " * N["threads/chain"] * ")")
                             } else expression("AD_Stan_WCP (WCP-only: " * N[chunks] * " = " * N["threads/chain"] * ")"))
        best_overall_label <-  "Best measured configuration"
        ##
        ## was: y_axis_title <-  "Seconds per leapfrog step\n(slowest chain, log scale)"
        y_axis_title <-  "Seconds per burn-in iteration\n(slowest chain, log scale)"
        ##
        ## ---- AD_Stan (no chunking, N_threads/chain = 1):
        ##
        no_chunking <-  device_rows[device_rows$algorithm == "AD_Stan", , drop = FALSE]
        no_chunking$series <-  factor("no_chunking", levels = series_levels)
        ##
        ## ---- best N_chunks at each (N, N_threads/chain) over AD_Stan_tape_chunked (N_threads/chain = 1),
        ##      AD_Stan_WCP and AD_Stan_WCP_chunking; its N_threads/chain = 1 point is the best tape chunking
        ##      only configuration:
        ##
        best_by_threads <-  dplyr::group_by( device_rows[device_rows$algorithm != "AD_Stan", , drop = FALSE],
                                             .data$N, .data$n_threads_per_chain)
        ## was: best_by_threads <-  dplyr::slice(best_by_threads, which.min(.data$sec_per_step))
        best_by_threads <-  dplyr::slice( best_by_threads,
                                          which.min(.data$sec_per_burnin_iteration_slowest_chain))
        best_by_threads <-  as.data.frame(dplyr::ungroup(best_by_threads))
        best_by_threads$series <-  factor("best_chunks_at_each_threads_per_chain", levels = series_levels)
        ##
        ## ---- WCP-only (AD_Stan_WCP, N_chunks = N_threads/chain > 1):
        ##
        WCP_only <-  device_rows[device_rows$algorithm == "AD_Stan_WCP", , drop = FALSE]
        WCP_only$series <-  factor("WCP_only", levels = series_levels)
        ##
        ## ---- best measured configuration per N (over every configuration, as in the shootout summary):
        ##
        best_overall <-  dplyr::group_by(device_rows, .data$N)
        ## was: best_overall <-  dplyr::slice(best_overall, which.min(.data$sec_per_step))
        best_overall <-  dplyr::slice(best_overall, which.min(.data$sec_per_burnin_iteration_slowest_chain))
        best_overall <-  as.data.frame(dplyr::ungroup(best_overall))
        best_overall$best_overall_series <-  best_overall_label
        ##
        ## ---- one panel per N:
        ##
        panel_list <-  list()
        ##
        for (N_index in seq_along(N_values)) {

              this_N <-  N_values[N_index]
              ##
              panel_values <-  figure_values[ figure_values$device == device & figure_values$N == this_N, ,
                                              drop = FALSE]
              if (nrow(panel_values) != 1) {
                    stop(paste0( "fn_paper1_burnin_cmdstanr_figure: expected one row of figure values for ",
                                 device, ", N = ", this_N, ", found ", nrow(panel_values), "."))
              }
              ##
              WCP_fold_text   <-  fn_paper1_burnin_format_fold(fold_vector = panel_values$WCP_fold)
              total_fold_text <-  fn_paper1_burnin_format_fold(fold_vector = panel_values$total_fold)
              ##
              panel_fold_label <-  paste0("WCP ", WCP_fold_text)
              if (show_total_fold_label) {
                    panel_fold_label <-  paste0("total ", total_fold_text, "\n", panel_fold_label)
              }
              ##
              panel_list[[length(panel_list) + 1]] <-  fn_paper1_burnin_cmdstanr_figure_panel(
                    panel_best_by_threads = best_by_threads[best_by_threads$N == this_N, , drop = FALSE],
                    panel_no_chunking     = no_chunking[no_chunking$N == this_N, , drop = FALSE],
                    panel_WCP_only        = WCP_only[WCP_only$N == this_N, , drop = FALSE],
                    panel_best_overall    = best_overall[best_overall$N == this_N, , drop = FALSE],
                    panel_fold_label      = panel_fold_label,
                    show_right_strips     = (N_index == length(N_values)),
                    series_levels         = series_levels,
                    series_labels         = series_labels,
                    best_overall_label    = best_overall_label,
                    y_axis_title          = y_axis_title)

        }
        ##
        ## ---- assemble the row: one collected legend at the bottom, one x-axis title and one y-axis title:
        ##
        burnin_figure <-  patchwork::wrap_plots(panel_list, ncol = length(N_values), byrow = TRUE) +
              patchwork::plot_layout(guides = "collect", axis_titles = "collect") &
              ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")
        ##
        dir.create(path = dirname(file_path), recursive = TRUE, showWarnings = FALSE)
        ggplot2::ggsave( filename = file_path,
                         plot     = burnin_figure,
                         width    = width_inches,
                         height   = height_inches,
                         dpi      = dpi)
        ##
        message(NicoStan::colourise( text = paste0( "Stan via cmdstanr burn-in, ", device, ": wrote ",
                                                    length(panel_list), " panel(s) to ", file_path),
                                     fg   = "green"))
        ##
        return(invisible(list(figure = burnin_figure, file_path = file_path)))

}























