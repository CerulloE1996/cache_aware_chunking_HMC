##
## =======================================================================================================================================
## R_fns_alg_paper_1_burnin_report.R
##
## Paper 1 burn-in report: turns the saved PS1 burn-in benchmark outputs (chunks x WCP x burn-in chain count, one seconds-per-iteration
## measurement per configuration and run) into the manuscript bundle Files/Generated/Burnin_chunks_WCP.
##
## This file only defines functions; alg_paper_1_burnin_report.R in this directory is the runner that supplies every setting explicitly
## and calls them. This file never runs sampling, compiles anything, or touches ps_1_burnin_optimizing_N_chunks_and_WCP_threads.R /
## R_fn_run_ps_1_burnin_benchmark.R - it only reads the .rds files those scripts already saved.
##
## Column schema of the saved .rds files (see R_fn_run_ps_1_burnin_benchmark.R, R_fn_run_ps_1_burnin_benchmark() and
## R_fn_ps1_burnin_read_results()):
##   device, N, n_chains_burnin, num_chunks_burnin, n_threads_WCP_burnin, n_threads_total, run_number, L_main,
##   n_timed_iters, total_timed_seconds, median_sec_per_iter, mean_sec_per_iter, n_divs
## with device stored as "Local_HPC" (HPC) or "Laptop" (Laptop) - NOT "HPC" - and one .rds file per N, named
## "<HPC_ or Laptop_><ps1_burnin_benchmark_N><N>_L<L_main>_n_runs<n_runs>.rds" inside the supplied output_dir.
##

require(dplyr)
require(tidyr)
require(tibble)
require(ggplot2)

##
## ---- Small local helpers (file-writing and formatting only; no defaults, no global lookups): ------------------------------------------
##
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
fn_paper1_burnin_device_field <-  function( device_argument ) {

        ##
        ## ---- Translate the runner's "HPC"/"Laptop" selector into the canonical "device" value the .rds rows actually carry:
        ##
        if (!is.character(x = device_argument) || length(x = device_argument) != 1 || is.na(x = device_argument)) {
              stop("fn_paper1_burnin_device_field: device_argument must be a single non-missing character string.")
        }
        if (device_argument == "HPC") {
              return("Local_HPC")
        } else if (device_argument == "Laptop") {
              return("Laptop")
        } else {
              stop("fn_paper1_burnin_device_field: device_argument must be 'HPC' or 'Laptop', got '", device_argument, "'.")
        }

}
##
fn_paper1_burnin_file_prefix <-  function( device_argument ) {

        if (device_argument == "HPC") {
              return("HPC_")
        } else if (device_argument == "Laptop") {
              return("Laptop_")
        } else {
              stop("fn_paper1_burnin_file_prefix: device_argument must be 'HPC' or 'Laptop', got '", device_argument, "'.")
        }

}
##
fn_paper1_burnin_format_seconds_per_iter <-  function( seconds_per_iter_vector ) {

        formatted_strings <-  formatC(x = seconds_per_iter_vector, digits = 3, format = "g")
        formatted_strings[is.na(x = seconds_per_iter_vector)] <-  "---"
        return(formatted_strings)

}
##
fn_paper1_burnin_format_speed_up <-  function( speed_up_vector ) {

        formatted_strings <-  paste0(formatC(x = speed_up_vector, digits = 2, format = "f"), "$\\times$")
        formatted_strings[is.na(x = speed_up_vector)] <-  "---"
        return(formatted_strings)

}
##
## ---- Plain-text speed-up labels for figures (e.g. "4.23-fold"), same 2 decimal places as fn_paper1_burnin_format_speed_up:
##
fn_paper1_burnin_format_fold <-  function( fold_vector ) {

        formatted_strings <-  paste0(formatC(x = fold_vector, digits = 2, format = "f"), "-fold")
        formatted_strings[is.na(x = fold_vector)] <-  "n/a"
        return(formatted_strings)

}

##
## ---- fn_paper1_burnin_read_outputs: read every saved PS1 burn-in .rds for ONE device into tidy rows: ----------------------------------
##
fn_paper1_burnin_read_outputs <-  function( output_dir,
                                            device
) {

        ##
        ## ---- Argument checks (helpers here take no defaults; every setting is explicit and self-checking): -----------------------------
        ##
        if (missing(x = output_dir) || !is.character(x = output_dir) || length(x = output_dir) != 1 || is.na(x = output_dir)) {
              stop("fn_paper1_burnin_read_outputs: output_dir must be supplied as a single non-missing character string.")
        }
        if (missing(x = device)) {
              stop("fn_paper1_burnin_read_outputs: device must be supplied explicitly ('HPC' or 'Laptop').")
        }
        if (!dir.exists(paths = output_dir)) {
              stop("fn_paper1_burnin_read_outputs: output_dir does not exist: ", output_dir)
        }
        expected_device_field <-  fn_paper1_burnin_device_field(device_argument = device)
        file_prefix           <-  fn_paper1_burnin_file_prefix(device_argument = device)
        ##
        required_columns <-  c("device", "N", "n_chains_burnin", "num_chunks_burnin", "n_threads_WCP_burnin",
                               "n_threads_total", "run_number", "L_main", "n_timed_iters", "total_timed_seconds",
                               "mean_sec_per_iter", "n_divs")
        ##
        ## ---- Discover every saved N for this device (no hard-coded N grid; whatever is on disk is read): -------------------------------
        ##
        candidate_file_paths <-  list.files( path       = output_dir,
                                             pattern    = paste0("^", file_prefix, "ps1_burnin_benchmark_N[0-9]+_L[0-9]+_n_runs[0-9]+\\.rds$"),
                                             full.names = TRUE)
        ##
        if (length(x = candidate_file_paths) == 0) {
              message(paste0("\033[36m", "fn_paper1_burnin_read_outputs: no saved '", device,
                             "' burn-in .rds files found under ", output_dir, "; returning zero-row tidy rows.", "\033[0m"))
              empty_rows <-  tibble::tibble( device               = character(0),
                                             N                    = numeric(0),
                                             n_chains_burnin      = numeric(0),
                                             num_chunks_burnin    = numeric(0),
                                             n_threads_WCP_burnin = numeric(0),
                                             n_threads_total      = numeric(0),
                                             run_number           = numeric(0),
                                             L_main               = numeric(0),
                                             n_timed_iters        = numeric(0),
                                             total_timed_seconds  = numeric(0),
                                             mean_sec_per_iter    = numeric(0),
                                             median_sec_per_iter  = numeric(0),
                                             n_divs               = numeric(0),
                                             source_file           = character(0))
              return(empty_rows)
        }
        ##
        per_file_rows <-  vector(mode = "list", length = length(x = candidate_file_paths))
        ##
        for (file_index in seq_along(along.with = candidate_file_paths)) {

              this_file_path <-  candidate_file_paths[file_index]
              this_result    <-  readRDS(file = this_file_path)
              ##
              if (!is.data.frame(x = this_result)) {
                    stop("fn_paper1_burnin_read_outputs: ", this_file_path, " does not contain a data frame - inspect it before using it.")
              }
              ##
              ## ---- Older saved.rds files (e.g. the N500/N50000 HPC files on disk as of) used earlier column names;
              ## translate those in memory only, exactly as R_fn_ps1_burnin_read_results() does, without touching the file on disk:
              ##
              previous_to_current_column_names <-  c( n_timed_iters       = "n_timed_iterations",
                                                       mean_sec_per_iter   = "mean_seconds_per_iteration",
                                                       median_sec_per_iter = "median_seconds_per_iteration",
                                                       n_divs              = "n_divergent_transitions")
              for (current_column_name in names(x = previous_to_current_column_names)) {

                    previous_column_name <-  previous_to_current_column_names[[current_column_name]]
                    if (!current_column_name %in% names(x = this_result) && previous_column_name %in% names(x = this_result)) {
                          this_result[[current_column_name]] <-  this_result[[previous_column_name]]
                          message(paste0("\033[36m", "fn_paper1_burnin_read_outputs: ", this_file_path, " used the older column name '",
                                         previous_column_name, "'; read as '", current_column_name, "'.", "\033[0m"))
                    }

              }
              if (!all(required_columns %in% names(x = this_result))) {
                    stop("fn_paper1_burnin_read_outputs: ", this_file_path, " is missing required column(s): ",
                        paste(setdiff(x = required_columns, y = names(x = this_result)), collapse = ", "))
              }
              if (!"median_sec_per_iter" %in% names(x = this_result)) {
                    this_result$median_sec_per_iter <-  NA_real_
              }
              ##
              ## ---- Self-check: every row in a "<device>_..." file must carry the expected canonical device value: ------------------------
              ##
              observed_device_values <-  unique(x = this_result$device)
              if (!all(observed_device_values == expected_device_field)) {
                    stop("fn_paper1_burnin_read_outputs: ", this_file_path, " was matched to device = '", device,
                        "' by its filename, but its saved 'device' column contains ", paste(observed_device_values, collapse = ", "),
                        " (expected only '", expected_device_field, "'). Fix the mismatch before trusting this report.")
              }
              ##
              usable_rows <-  complete.cases(this_result[required_columns]) &
                              rowSums(!is.finite(as.matrix(this_result[setdiff(x = required_columns, y = "device")]))) == 0 &
                              this_result$n_timed_iters > 0 & this_result$total_timed_seconds > 0 &
                              this_result$mean_sec_per_iter > 0
              usable_rows[is.na(x = usable_rows)] <-  FALSE
              if (any(!usable_rows)) {
                    message(paste0("\033[36m", "fn_paper1_burnin_read_outputs: skipping ", sum(!usable_rows),
                                   " incomplete/invalid row(s) in ", this_file_path, "\033[0m"))
              }
              this_result <-  dplyr::filter(.data = this_result, usable_rows)
              ##
              ## ---- Divergent rows are KEPT and flagged, never dropped: convergence is reported as its own column. -----------------
              ##
              this_result$converged <-  this_result$n_divs == 0
              if (any(!this_result$converged)) {
                    message(paste0("\033[36m", "fn_paper1_burnin_read_outputs: ", sum(!this_result$converged),
                                   " row(s) with divergent transitions in ", this_file_path, " (kept, converged = FALSE)", "\033[0m"))
              }
              this_result$source_file <-  this_file_path
              ##
              per_file_rows[[file_index]] <-  tibble::as_tibble(x = this_result[, c(required_columns, "median_sec_per_iter", "converged", "source_file")])

        }
        ##
        tidy_rows <-  dplyr::bind_rows(per_file_rows)
        message(paste0("\033[36m", "fn_paper1_burnin_read_outputs: read ", nrow(x = tidy_rows), " usable row(s) for device = '", device,
                       "' from ", length(x = candidate_file_paths), " file(s) under ", output_dir, "\033[0m"))
        return(tidy_rows)

}

##
## ---- fn_paper1_burnin_exclude_configurations: drop grid cells the study ran but the report should not show: ------------
##
## excluded_WCP_given_N is a list keyed by N (as character), each element the N_WCP values to drop at that N, e.g. list("500" = 6):
## the HPC ps1 grid ran N_WCP = 6 at N = 500 for every burn-in chain count but skipped it for 8 chains at N = 2,500, so it is left
## out of the report at N = 500 to keep the two N grids comparable. Every dropped row is counted in the message; nothing is hidden
## silently, and the saved .rds outputs are untouched.
##
fn_paper1_burnin_exclude_configurations <-  function( rows,
                                                      excluded_WCP_given_N) {

        if (!is.list(x = excluded_WCP_given_N)) {
              stop("fn_paper1_burnin_exclude_configurations: excluded_WCP_given_N must be a list keyed by N (use list() for none).")
        }
        if (length(x = excluded_WCP_given_N) == 0) return(rows)
        ##
        rows_to_drop <-  rep(x = FALSE, times = nrow(x = rows))
        for (N_label in names(x = excluded_WCP_given_N)) {

              excluded_WCP_values <-  excluded_WCP_given_N[[N_label]]
              rows_to_drop <-  rows_to_drop | (rows$N == as.numeric(x = N_label) & rows$n_threads_WCP_burnin %in% excluded_WCP_values)
              message(paste0("\033[36m", "fn_paper1_burnin_exclude_configurations: N = ", N_label, ", N_WCP in {",
                             paste(excluded_WCP_values, collapse = ", "), "}: dropping ",
                             sum(rows$N == as.numeric(x = N_label) & rows$n_threads_WCP_burnin %in% excluded_WCP_values),
                             " row(s) from the report.", "\033[0m"))

        }
        return(rows[!rows_to_drop, , drop = FALSE])

}
##
## ---- fn_paper1_burnin_summarise_rows (private helper): median over repeats within (N, chains, chunks, WCP): ---------------------------
##
fn_paper1_burnin_summarise_rows <-  function( rows,
                                              device
) {

        expected_device_field <-  fn_paper1_burnin_device_field(device_argument = device)
        ##
        device_rows <-  dplyr::filter(.data = rows, .data$device == expected_device_field)
        if (nrow(x = device_rows) == 0) {
              stop("fn_paper1_burnin_summarise_rows: no rows for device = '", device, "' (expected device field '",
                  expected_device_field, "') - check that fn_paper1_burnin_read_outputs was called for this device.")
        }
        ##
        configuration_columns <-  c("device", "N", "n_chains_burnin", "num_chunks_burnin", "n_threads_WCP_burnin", "n_threads_total")
        ##
        configuration_summary <-  device_rows %>%
              dplyr::group_by(dplyr::across(.cols = dplyr::all_of(x = configuration_columns))) %>%
              dplyr::summarise( sec_per_iter      = mean(x = .data$mean_sec_per_iter, na.rm = TRUE),
                                n_runs_available  = dplyr::n_distinct(.data$run_number, na.rm = TRUE),
                                .groups           = "drop")
        ##
        return(configuration_summary)

}

##
## ---- fn_paper1_burnin_best_table_tex: per (N, burn-in chain count), the best configuration and both speed-ups: ------------------------
##
fn_paper1_burnin_best_table_tex <-  function( rows,
                                              device,
                                              caption,
                                              label,
                                              file_path
) {

        for (required_string_argument in list(device = device, caption = caption, label = label, file_path = file_path)) {
              if (!is.character(x = required_string_argument) || length(x = required_string_argument) != 1 ||
                  is.na(x = required_string_argument)) {
                    stop("fn_paper1_burnin_best_table_tex: device, caption, label and file_path must all be single non-missing strings.")
              }
        }
        ##
        configuration_summary <-  fn_paper1_burnin_summarise_rows(rows = rows, device = device)
        ##
        ## ---- Best overall configuration per (N, chains):
        ##
        best_overall <-  configuration_summary %>%
              dplyr::group_by(.data$N, .data$n_chains_burnin) %>%
              dplyr::slice(which.min(.data$sec_per_iter)) %>%
              dplyr::ungroup() %>%
              dplyr::rename( best_num_chunks_burnin    = num_chunks_burnin,
                             best_n_threads_WCP_burnin = n_threads_WCP_burnin,
                             best_n_threads_total      = n_threads_total,
                             best_sec_per_iter          = sec_per_iter,
                             best_n_runs_available      = n_runs_available) %>%
              dplyr::select(N, n_chains_burnin, best_num_chunks_burnin, best_n_threads_WCP_burnin,
                            best_n_threads_total, best_sec_per_iter, best_n_runs_available)
        ##
        ## ---- Baseline (1 chunk, WCP = 1) per (N, chains) - what chunking and WCP are both compared against:
        ##
        baseline_serial <-  configuration_summary %>%
              dplyr::filter(.data$num_chunks_burnin == 1, .data$n_threads_WCP_burnin == 1) %>%
              dplyr::select(N, n_chains_burnin, baseline_sec_per_iter = sec_per_iter)
        ##
        ## ---- Best chunking-only configuration (WCP = 1, any chunk count) per (N, chains) - isolates chunking's own gain:
        ##
        best_chunking_only <-  configuration_summary %>%
              dplyr::filter(.data$n_threads_WCP_burnin == 1) %>%
              dplyr::group_by(.data$N, .data$n_chains_burnin) %>%
              dplyr::slice(which.min(.data$sec_per_iter)) %>%
              dplyr::ungroup() %>%
              dplyr::select(N, n_chains_burnin, best_chunking_only_sec_per_iter = sec_per_iter)
        ##
        table_rows <-  best_overall %>%
              dplyr::left_join(y = baseline_serial,     by = c("N", "n_chains_burnin")) %>%
              dplyr::left_join(y = best_chunking_only,   by = c("N", "n_chains_burnin")) %>%
              dplyr::mutate( speed_up_vs_1_chunk_WCP_1        = .data$baseline_sec_per_iter        / .data$best_sec_per_iter,
                             speed_up_vs_best_chunks_WCP_1     = .data$best_chunking_only_sec_per_iter / .data$best_sec_per_iter) %>%
              dplyr::arrange(.data$N, .data$n_chains_burnin)
        ##
        missing_baseline_combinations <-  table_rows %>% dplyr::filter(is.na(x = .data$baseline_sec_per_iter))
        if (nrow(x = missing_baseline_combinations) > 0) {
              message(paste0("\033[36m", "fn_paper1_burnin_best_table_tex: no (1 chunk, WCP = 1) baseline measurement for ",
                             nrow(x = missing_baseline_combinations),
                             " (N, chains) combination(s) on device = '", device, "'; their speed-up vs 1 chunk/WCP 1 is left blank.",
                             "\033[0m"))
        }
        ##
        ## ---- Emit the LaTeX table (booktabs-free, matching the other Paper 1 tables in this directory): -------------------------------
        ##
        table_lines <-  c( "\\begin{table}[H]",
                           "\\centering",
                           paste0("\\caption{", caption, "}"),
                           paste0("\\label{", label, "}"),
                           "\\begin{tabular}{rrrrrrrr}",
                           "\\hline",
                           paste0("$N$ & Burn-in & Best & Best & Total & Sec./ & Speed-up vs & Speed-up vs \\\\"),
                           paste0(" & chains & chunks & $N_{\\text{WCP}}$ & threads & iter & 1 chunk/WCP 1 & best chunks (WCP 1) \\\\ \\hline"))
        ##
        for (row_index in seq_len(length.out = nrow(x = table_rows))) {

              this_row <-  table_rows[row_index, ]
              table_lines <-  c( table_lines,
                                 paste0( fn_paper1_format_number_commas_from_10000(this_row$N),
                                        " & ", this_row$n_chains_burnin,
                                        " & ", this_row$best_num_chunks_burnin,
                                        " & ", this_row$best_n_threads_WCP_burnin,
                                        " & ", this_row$best_n_threads_total,
                                        " & ", fn_paper1_burnin_format_seconds_per_iter(seconds_per_iter_vector = this_row$best_sec_per_iter),
                                        " & ", fn_paper1_burnin_format_speed_up(speed_up_vector = this_row$speed_up_vs_1_chunk_WCP_1),
                                        " & ", fn_paper1_burnin_format_speed_up(speed_up_vector = this_row$speed_up_vs_best_chunks_WCP_1),
                                        " \\\\"))

        }
        table_lines <-  c(table_lines, "\\hline", "\\end{tabular}", "\\end{table}")
        ##
        dir.create(path = dirname(path = file_path), recursive = TRUE, showWarnings = FALSE)
        writeLines(text = table_lines, con = file_path)
        message(paste0("\033[36m", "fn_paper1_burnin_best_table_tex: wrote ", nrow(x = table_rows), " row(s) to ", file_path, "\033[0m"))
        ##
        return(invisible(table_rows))

}

##
## ---- fn_paper1_burnin_figure: seconds per burn-in iteration vs N_WCP, one line per chunk count, facets N x burn-in chains, log-y: ------
##
fn_paper1_burnin_figure <-  function( rows,
                                      device,
                                      file_path
) {

        for (required_string_argument in list(device = device, file_path = file_path)) {
              if (!is.character(x = required_string_argument) || length(x = required_string_argument) != 1 ||
                  is.na(x = required_string_argument)) {
                    stop("fn_paper1_burnin_figure: device and file_path must both be single non-missing strings.")
              }
        }
        ##
        configuration_summary <-  fn_paper1_burnin_summarise_rows(rows = rows, device = device)
        ##
        plot_data <-  configuration_summary %>%
              dplyr::mutate( num_chunks_burnin_label = factor(x = .data$num_chunks_burnin,
                                                              levels = sort(x = unique(x = .data$num_chunks_burnin))),
                             N_label                  = factor(x = paste0("N = ", fn_paper1_format_number_commas_from_10000(.data$N)),
                                                                levels = paste0("N = ", fn_paper1_format_number_commas_from_10000(sort(x = unique(x = .data$N))))),
                             chains_label              = factor(x = paste0(.data$n_chains_burnin, " burn-in chains"),
                                                                levels = paste0(sort(x = unique(x = .data$n_chains_burnin)),
                                                                                " burn-in chains")))
        ##
        burnin_figure <-  ggplot2::ggplot( data    = plot_data,
                                           mapping = ggplot2::aes(x = .data$n_threads_WCP_burnin, y = .data$sec_per_iter,
                                                                  colour = .data$num_chunks_burnin_label,
                                                                  group  = .data$num_chunks_burnin_label)) +
              ggplot2::geom_line(linewidth = 0.7) +
              ggplot2::geom_point(size = 2) +
              ggplot2::scale_y_log10() +
              ggplot2::facet_grid(rows = ggplot2::vars(.data$chains_label), cols = ggplot2::vars(.data$N_label), scales = "free_y") +
              ggplot2::labs( x      = expression(N[WCP] ~ "(within-chain-parallelism threads per burn-in chain)"),
                             y      = "Seconds per burn-in iteration (log scale)",
                             colour = expression(N[chunks])) +
              ggplot2::theme_bw(base_size = 14) +
              ggplot2::theme(legend.position = "bottom")
        ##
        dir.create(path = dirname(path = file_path), recursive = TRUE, showWarnings = FALSE)
        ggplot2::ggsave( filename = file_path,
                         plot     = burnin_figure,
                         width    = 10,
                         height   = 3 * dplyr::n_distinct(plot_data$chains_label) + 1,
                         dpi      = 150)
        message(paste0("\033[36m", "fn_paper1_burnin_figure: wrote figure (", nrow(x = plot_data), " point(s)) to ", file_path, "\033[0m"))
        ##
        return(invisible(burnin_figure))

}

##
## ---- fn_paper1_burnin_figure_best_chunks_panel (private helper): one (N, burn-in chains) panel of fn_paper1_burnin_figure_best_chunks: ---
##
## x = N_WCP on a log2 scale with breaks only at the N_WCP values measured in this panel; y = seconds per burn-in iteration (log10),
## with this panel's own range. The N_chunks = 1, N_WCP = 1 measurement is an open square with a dashed line at its level, the best
## N_chunks at each N_WCP is a solid line (each point labelled with its N_chunks), and the best measured configuration is outlined.
##
fn_paper1_burnin_figure_best_chunks_panel <-  function( panel_best_by_WCP,
                                                        panel_no_chunking,
                                                        panel_best_overall,
                                                        panel_speed_up_label,
                                                        show_N_strip,
                                                        show_chains_strip,
                                                        series_levels,
                                                        series_labels,
                                                        best_overall_label,
                                                        panel_WCP_only = NULL
) {

        ##
        ## ---- x breaks only at the N_WCP values measured in this panel:
        ##
        ## measured_WCP_values <-  sort(x = unique(x = c(panel_best_by_WCP$n_threads_WCP_burnin, panel_no_chunking$n_threads_WCP_burnin)))
        measured_WCP_values <-  sort(x = unique(x = c(panel_best_by_WCP$n_threads_WCP_burnin, panel_no_chunking$n_threads_WCP_burnin,
                                                      if (is.null(x = panel_WCP_only)) NULL else panel_WCP_only$n_threads_WCP_burnin)))
        ##
        ## ---- two rows of x labels where adjacent measured N_WCP values are close on the log2 axis (e.g. 16, 24, 32, 44):
        ##
        measured_WCP_log2     <-  log2(x = measured_WCP_values)
        smallest_relative_gap <-  if (length(x = measured_WCP_values) > 1) min(diff(x = measured_WCP_log2)) / diff(x = range(measured_WCP_log2)) else 1
        x_axis_label_rows     <-  if (smallest_relative_gap < 0.10) 2 else 1
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
        ## ---- N_chunks labels sit below each point, except above a point that is slower than its left-hand neighbour (so that the label
        ##      does not sit on the rising line segment); the label of the outlined best point sits further below, clear of its circle:
        ##
        panel_best_by_WCP   <-  panel_best_by_WCP[order(panel_best_by_WCP$n_threads_WCP_burnin), , drop = FALSE]
        rises_from_left     <-  c(FALSE, diff(x = panel_best_by_WCP$sec_per_iter) > 0)
        is_best_overall     <-  panel_best_by_WCP$n_threads_WCP_burnin == panel_best_overall$n_threads_WCP_burnin[1]
        panel_best_by_WCP$label_vjust <-  ifelse(test = rises_from_left, yes = -1.0, no = ifelse(test = is_best_overall, yes = 3.0, no = 2.0))
        ##
        ## ---- Series styles (first level = no chunking reference, second level = best N_chunks at each N_WCP):
        ##
        ## series_colours   <-  c("#D55E00", "#0072B2")
        ## series_linetypes <-  c("dashed",  "solid")
        ## series_shapes    <-  c(0,         16)
        ## (third level = WCP-only, i.e. N_chunks = N_threads/chain, drawn in green with triangles)
        series_colours   <-  c("#D55E00", "#0072B2", "#009E73")[seq_along(along.with = series_levels)]
        series_linetypes <-  c("dashed",  "solid",   "22")[seq_along(along.with = series_levels)]
        series_shapes    <-  c(0,         16,        17)[seq_along(along.with = series_levels)]
        names(x = series_colours)   <-  series_levels
        names(x = series_linetypes) <-  series_levels
        names(x = series_shapes)    <-  series_levels
        ##
        best_overall_fill <-  c(NA)
        names(x = best_overall_fill) <-  best_overall_label
        ##
        panel_plot <-  ggplot2::ggplot() +
              ##
              ## ---- dashed reference level of N_chunks = 1, N_WCP = 1 across the whole panel:
              ##
              ggplot2::geom_hline( data        = panel_no_chunking,
                                   mapping     = ggplot2::aes(yintercept = .data$sec_per_iter,
                                                              colour     = .data$series,
                                                              linetype   = .data$series),
                                   linewidth   = 0.6) +
              ##
              ## ---- best N_chunks at each N_WCP:
              ##
              ggplot2::geom_line( data        = panel_best_by_WCP,
                                  mapping     = ggplot2::aes(x        = .data$n_threads_WCP_burnin,
                                                             y        = .data$sec_per_iter,
                                                             colour   = .data$series,
                                                             linetype = .data$series),
                                  linewidth   = 0.8) +
              ggplot2::geom_point( data        = panel_best_by_WCP,
                                   mapping     = ggplot2::aes(x      = .data$n_threads_WCP_burnin,
                                                              y      = .data$sec_per_iter,
                                                              colour = .data$series,
                                                              shape  = .data$series),
                                   size        = 2.2) +
              ggplot2::geom_text( data        = panel_best_by_WCP,
                                  mapping     = ggplot2::aes(x     = .data$n_threads_WCP_burnin,
                                                             y     = .data$sec_per_iter,
                                                             label = .data$num_chunks_burnin,
                                                             vjust = .data$label_vjust),
                                  colour      = series_colours[[2]],
                                  size        = 3.3,
                                  show.legend = FALSE) +
              ##
              ## ---- WCP-only (N_chunks = N_threads/chain), where measured:
              ##
              (if (is.null(x = panel_WCP_only) || nrow(x = panel_WCP_only) == 0) NULL else
                    list( ggplot2::geom_line( data        = panel_WCP_only,
                                              mapping     = ggplot2::aes(x        = .data$n_threads_WCP_burnin,
                                                                         y        = .data$sec_per_iter,
                                                                         colour   = .data$series,
                                                                         linetype = .data$series),
                                              linewidth   = 0.8),
                          ggplot2::geom_point( data        = panel_WCP_only,
                                               mapping     = ggplot2::aes(x      = .data$n_threads_WCP_burnin,
                                                                          y      = .data$sec_per_iter,
                                                                          colour = .data$series,
                                                                          shape  = .data$series),
                                               size        = 2.6))) +
              ##
              ## ---- N_chunks = 1, N_WCP = 1 measurement (drawn after the best line, so it stays visible where the two coincide):
              ##
              ggplot2::geom_point( data        = panel_no_chunking,
                                   mapping     = ggplot2::aes(x      = .data$n_threads_WCP_burnin,
                                                              y      = .data$sec_per_iter,
                                                              colour = .data$series,
                                                              shape  = .data$series),
                                   size        = 3.4,
                                   stroke      = 1.1) +
              ##
              ## ---- best measured (N_chunks, N_WCP) configuration, outlined:
              ##
              ggplot2::geom_point( data        = panel_best_overall,
                                   mapping     = ggplot2::aes(x    = .data$n_threads_WCP_burnin,
                                                              y    = .data$sec_per_iter,
                                                              fill = .data$best_overall_series),
                                   shape       = 21,
                                   size        = 5.0,
                                   stroke      = 1.0,
                                   colour      = "black") +
              ##
              ## ---- speed-ups of the best configuration relative to both references (as in fn_paper1_burnin_best_table_tex):
              ##
              ggplot2::annotate( geom  = "text",
                                 x     = Inf,
                                 y     = Inf,
                                 label = panel_speed_up_label,
                                 hjust = 1.06,
                                 vjust = 1.25,
                                 size  = 3.4) +
              ggplot2::scale_x_continuous( transform    = "log2",
                                           breaks       = measured_WCP_values,
                                           labels       = measured_WCP_values,
                                           minor_breaks = NULL,
                                           expand       = ggplot2::expansion(mult = 0.08),
                                           guide        = ggplot2::guide_axis(n.dodge = x_axis_label_rows)) +
              ggplot2::scale_y_log10( breaks = y_breaks_function,
                                      expand = ggplot2::expansion(mult = c(0.20, 0.45)),
                                      guide  = ggplot2::guide_axis(check.overlap = TRUE)) +
              ggplot2::scale_colour_manual(   name = NULL, values = series_colours,   breaks = series_levels, labels = series_labels) +
              ggplot2::scale_linetype_manual( name = NULL, values = series_linetypes, breaks = series_levels, labels = series_labels) +
              ggplot2::scale_shape_manual(    name = NULL, values = series_shapes,    breaks = series_levels, labels = series_labels) +
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
        ##      the burn-in chain strips are plotmath (N_chains = ...), so they are parsed:
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
## ---- fn_paper1_burnin_figure_best_chunks: seconds per burn-in iteration vs N_WCP against both burn-in references, per (N, chains): ------
##
## Manuscript version of the burn-in figure (fn_paper1_burnin_figure above, one line per chunk count on shared axes, is kept unchanged).
## Rows = burn-in chains, columns = N; each panel has its own x and y range. Per panel:
##   - "No chunking, no WCP" = the N_chunks = 1, N_WCP = 1 measurement (open square, dashed line at its level) = reference (i);
##   - "Best N_chunks at each N_WCP" = the fastest measured chunk count at each N_WCP (solid line, points labelled with N_chunks);
##     its N_WCP = 1 point is the best chunking-only setting = reference (ii);
##   - the best measured (N_chunks, N_WCP) configuration (outlined point), and its speed-up relative to (i) ("total") and to (ii) ("WCP"),
##     computed exactly as in fn_paper1_burnin_best_table_tex().
## facet_grid() shares y along a row and x along a column, so each panel is its own ggplot and patchwork assembles the grid.
##
fn_paper1_burnin_figure_best_chunks <-  function( rows,
                                                  device,
                                                  file_path
) {

        for (required_string_argument in list(device = device, file_path = file_path)) {
              if (!is.character(x = required_string_argument) || length(x = required_string_argument) != 1 ||
                  is.na(x = required_string_argument)) {
                    stop("fn_paper1_burnin_figure_best_chunks: device and file_path must both be single non-missing strings.")
              }
        }
        ##
        configuration_summary <-  fn_paper1_burnin_summarise_rows(rows = rows, device = device)
        ##
        N_values      <-  sort(x = unique(x = configuration_summary$N))
        chains_values <-  sort(x = unique(x = configuration_summary$n_chains_burnin))
        ##
        configuration_summary <-  dplyr::mutate( .data        = configuration_summary,
                                                 N_label      = factor(x      = paste0("N = ", fn_paper1_format_number_commas_from_10000(.data$N)),
                                                                       levels = paste0("N = ", fn_paper1_format_number_commas_from_10000(N_values))),
                                                 ## chains_label = factor(x      = paste0(.data$n_chains_burnin, " burn-in chains"),
                                                 ##                       levels = paste0(chains_values, " burn-in chains")))
                                                 ## chains_label = factor(x      = paste0("N[chains]==", .data$n_chains_burnin, '~"(burn-in)"'),
                                                 ##                       levels = paste0("N[chains]==", chains_values, '~"(burn-in)"')))
                                                 chains_label = factor(x      = paste0('N["burn_chains"]==', .data$n_chains_burnin),
                                                                       levels = paste0('N["burn_chains"]==', chains_values)))
        ##
        ## series_levels      <-  c("no_chunking", "best_chunks_at_each_WCP")
        ## series_labels      <-  c(expression("No chunking, no WCP (" * N[chunks] * " = 1, " * N[WCP] * " = 1)"),
        ##                          expression("Best " * N[chunks] * " at each " * N[WCP] * " (point labels give " * N[chunks] * ")"))
        series_levels      <-  c("no_chunking", "best_chunks_at_each_WCP", "WCP_only")
        series_labels      <-  c(expression("MD_BayesMVP (" * N[chunks] * " = 1, " * N["threads/chain"] * " = 1)"),
                                 expression("Best " * N[chunks] * " at each " * N["threads/chain"] * " (point labels give " * N[chunks] * ")"),
                                 expression("MD_BayesMVP_WCP (WCP-only: " * N[chunks] * " = " * N["threads/chain"] * ")"))
        best_overall_label <-  "Best measured configuration"
        ##
        ## ---- Reference (i): N_chunks = 1, N_WCP = 1 per (N, chains):
        ##
        no_chunking <-  dplyr::filter(.data = configuration_summary, .data$num_chunks_burnin == 1, .data$n_threads_WCP_burnin == 1)
        no_chunking <-  dplyr::mutate(.data = no_chunking, series = factor(x = "no_chunking", levels = series_levels))
        ##
        ## ---- Best N_chunks at each (N, chains, N_WCP); its N_WCP = 1 row is reference (ii), the best chunking-only setting:
        ##
        best_by_WCP <-  dplyr::group_by(.data = configuration_summary, .data$N, .data$n_chains_burnin, .data$n_threads_WCP_burnin)
        best_by_WCP <-  dplyr::slice(.data = best_by_WCP, which.min(.data$sec_per_iter))
        best_by_WCP <-  dplyr::ungroup(x = best_by_WCP)
        best_by_WCP <-  dplyr::mutate(.data = best_by_WCP, series = factor(x = "best_chunks_at_each_WCP", levels = series_levels))
        ##
        ## ---- WCP-only (N_chunks = N_threads/chain, with N_threads/chain > 1) per (N, chains, N_WCP), where measured:
        ##
        WCP_only <-  dplyr::filter(.data = configuration_summary, .data$num_chunks_burnin == .data$n_threads_WCP_burnin,
                                   .data$n_threads_WCP_burnin > 1)
        WCP_only <-  dplyr::mutate(.data = WCP_only, series = factor(x = "WCP_only", levels = series_levels))
        ##
        ## ---- Best measured (N_chunks, N_WCP) per (N, chains):
        ##
        best_overall <-  dplyr::group_by(.data = configuration_summary, .data$N, .data$n_chains_burnin)
        best_overall <-  dplyr::slice(.data = best_overall, which.min(.data$sec_per_iter))
        best_overall <-  dplyr::ungroup(x = best_overall)
        best_overall <-  dplyr::mutate(.data = best_overall, best_overall_series = best_overall_label)
        ##
        ## ---- One panel per (chains, N), row-major (rows = burn-in chains, columns = N):
        ##
        panel_list    <-  list()
        panel_summary <-  list()
        ##
        for (chains_index in seq_along(along.with = chains_values)) {
              for (N_index in seq_along(along.with = N_values)) {

                    this_chains <-  chains_values[chains_index]
                    this_N      <-  N_values[N_index]
                    ##
                    panel_best_by_WCP  <-  best_by_WCP[best_by_WCP$N   == this_N & best_by_WCP$n_chains_burnin  == this_chains, , drop = FALSE]
                    panel_no_chunking  <-  no_chunking[no_chunking$N   == this_N & no_chunking$n_chains_burnin  == this_chains, , drop = FALSE]
                    panel_best_overall <-  best_overall[best_overall$N == this_N & best_overall$n_chains_burnin == this_chains, , drop = FALSE]
                    panel_WCP_only     <-  WCP_only[WCP_only$N         == this_N & WCP_only$n_chains_burnin     == this_chains, , drop = FALSE]
                    ##
                    if (nrow(x = panel_best_by_WCP) == 0) {
                          message(paste0("\033[36m", "fn_paper1_burnin_figure_best_chunks: no rows for N = ", this_N, ", ", this_chains,
                                         " burn-in chains on device = '", device, "'; leaving that panel empty.", "\033[0m"))
                          panel_list[[length(x = panel_list) + 1]] <-  patchwork::plot_spacer()
                          next
                    }
                    ##
                    ## ---- speed-ups, computed exactly as in fn_paper1_burnin_best_table_tex():
                    ##
                    best_sec_per_iter          <-  panel_best_overall$sec_per_iter
                    best_chunking_only_sec     <-  panel_best_by_WCP$sec_per_iter[panel_best_by_WCP$n_threads_WCP_burnin == 1]
                    no_chunking_sec            <-  panel_no_chunking$sec_per_iter
                    speed_up_vs_1_chunk_WCP_1  <-  if (length(x = no_chunking_sec) == 1)        no_chunking_sec        / best_sec_per_iter else NA_real_
                    speed_up_vs_best_chunks    <-  if (length(x = best_chunking_only_sec) == 1) best_chunking_only_sec / best_sec_per_iter else NA_real_
                    ##
                    panel_speed_up_label <-  paste0("total ", fn_paper1_burnin_format_fold(fold_vector = speed_up_vs_1_chunk_WCP_1), "\n",
                                                    "WCP ",   fn_paper1_burnin_format_fold(fold_vector = speed_up_vs_best_chunks))
                    ##
                    panel_summary[[length(x = panel_summary) + 1]] <-  data.frame( device                    = device,
                                                                                   N                         = this_N,
                                                                                   n_chains_burnin           = this_chains,
                                                                                   best_num_chunks_burnin    = panel_best_overall$num_chunks_burnin,
                                                                                   best_n_threads_WCP_burnin = panel_best_overall$n_threads_WCP_burnin,
                                                                                   best_sec_per_iter         = best_sec_per_iter,
                                                                                   no_chunking_sec_per_iter  = if (length(x = no_chunking_sec) == 1) no_chunking_sec else NA_real_,
                                                                                   best_chunking_only_sec    = if (length(x = best_chunking_only_sec) == 1) best_chunking_only_sec else NA_real_,
                                                                                   speed_up_vs_1_chunk_WCP_1 = speed_up_vs_1_chunk_WCP_1,
                                                                                   speed_up_vs_best_chunks   = speed_up_vs_best_chunks,
                                                                                   panel_speed_up_label      = panel_speed_up_label)
                    ##
                    panel_list[[length(x = panel_list) + 1]] <-  fn_paper1_burnin_figure_best_chunks_panel( panel_best_by_WCP    = panel_best_by_WCP,
                                                                                                            panel_no_chunking    = panel_no_chunking,
                                                                                                            panel_best_overall   = panel_best_overall,
                                                                                                            panel_speed_up_label = panel_speed_up_label,
                                                                                                            show_N_strip         = (chains_index == 1),
                                                                                                            show_chains_strip    = (N_index == length(x = N_values)),
                                                                                                            series_levels        = series_levels,
                                                                                                            series_labels        = series_labels,
                                                                                                            best_overall_label   = best_overall_label,
                                                                                                            panel_WCP_only       = panel_WCP_only)

              }
        }
        ##
        ## ---- Assemble the grid: one collected legend at the bottom, one x-axis title and one y-axis title:
        ##
        burnin_figure <-  patchwork::wrap_plots(panel_list, ncol = length(x = N_values), byrow = TRUE) +
              patchwork::plot_layout(guides = "collect", axis_titles = "collect") &
              ggplot2::theme(legend.position = "bottom", legend.box = "horizontal")
        ##
        dir.create(path = dirname(path = file_path), recursive = TRUE, showWarnings = FALSE)
        ggplot2::ggsave( filename = file_path,
                         plot     = burnin_figure,
                         width    = 10,
                         height   = 3 * length(x = chains_values) + 1.6,
                         dpi      = 150)
        message(paste0("\033[36m", "fn_paper1_burnin_figure_best_chunks: wrote figure (", length(x = panel_summary), " panel(s)) to ", file_path, "\033[0m"))
        ##
        return(invisible(list(figure = burnin_figure, panel_summary = do.call(what = rbind, args = panel_summary))))

}

##
## ---- fn_paper1_burnin_write_bundle: writes paper_sections/Burnin_chunks_WCP/{figures,tables,section.tex}, optionally to manuscript_dir: -
##
fn_paper1_burnin_write_bundle <-  function( rows,
                                            devices,
                                            output_dir,
                                            manuscript_dir
) {

        if (!is.character(x = devices) || length(x = devices) == 0 || anyNA(x = devices) || !all(devices %in% c("HPC", "Laptop"))) {
              stop("fn_paper1_burnin_write_bundle: devices must be a non-empty character vector of 'HPC' and/or 'Laptop'.")
        }
        if (!is.character(x = output_dir) || length(x = output_dir) != 1 || is.na(x = output_dir)) {
              stop("fn_paper1_burnin_write_bundle: output_dir must be a single non-missing character string.")
        }
        if (!is.null(x = manuscript_dir)) {
              if (!is.character(x = manuscript_dir) || length(x = manuscript_dir) != 1 || is.na(x = manuscript_dir)) {
                    stop("fn_paper1_burnin_write_bundle: manuscript_dir must be NULL or a single non-missing character string.")
              }
        }
        ##
        bundle_name <-  "Burnin_chunks_WCP"
        bundle_dir  <-  file.path(output_dir, "paper_sections", bundle_name)
        figure_dir  <-  file.path(bundle_dir, "figures")
        table_dir   <-  file.path(bundle_dir, "tables")
        dir.create(path = figure_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(path = table_dir,  recursive = TRUE, showWarnings = FALSE)
        ##
        device_full_name <-  c(HPC = "local HPC", Laptop = "laptop")
        ##
        section_tex_lines <-  c("% Generated by fn_paper1_burnin_write_bundle() from the saved PS1 burn-in benchmark outputs.")
        ##
        for (device_selector in devices) {

              expected_device_field <-  fn_paper1_burnin_device_field(device_argument = device_selector)
              device_rows            <-  dplyr::filter(.data = rows, .data$device == expected_device_field)
              if (nrow(x = device_rows) == 0) {
                    message(paste0("\033[36m", "fn_paper1_burnin_write_bundle: no rows for device = '", device_selector,
                                   "'; skipping its figure and table.", "\033[0m"))
                    next
              }
              ##
              table_file_name <-  paste0("table_burnin_best_", device_selector, ".tex")
              table_label      <-  paste0("table:burnin_best_", device_selector)
              fn_paper1_burnin_best_table_tex( rows      = rows,
                                               device    = device_selector,
                                               caption   = paste0("BayesMVP burn-in, ", device_full_name[[device_selector]],
                                                                  ": best measured (chunks, $N_{\\text{WCP}}$) at each $N$ and burn-in",
                                                                  " chain count, by mean seconds per burn-in iteration, with speed-up",
                                                                  " relative to a single chunk at $N_{\\text{WCP}} = 1$, and relative to",
                                                                  " the best chunk count at $N_{\\text{WCP}} = 1$ (chunking alone, before WCP)."),
                                               label     = table_label,
                                               file_path = file.path(table_dir, table_file_name))
              ##
              figure_file_name <-  paste0("figure_burnin_", device_selector, ".png")
              figure_label       <-  paste0("figure:burnin_", device_selector)
              ## fn_paper1_burnin_figure( rows      = rows,
              ##                          device    = device_selector,
              ##                          file_path = file.path(figure_dir, figure_file_name))
              fn_paper1_burnin_figure_best_chunks( rows      = rows,
                                                   device    = device_selector,
                                                   file_path = file.path(figure_dir, figure_file_name))
              ##
              section_tex_lines <-  c( section_tex_lines,
                                       "\\begin{figure}[H]",
                                       "\\centering",
                                       paste0("\\includegraphics[width=\\textwidth,height=0.8\\textheight,keepaspectratio]{",
                                             "Files/Generated/", bundle_name, "/figures/", figure_file_name, "}"),
                                       ## paste0("\\caption{BayesMVP burn-in, ", device_full_name[[device_selector]],
                                       ##       ": seconds per burn-in iteration against $N_{\\text{WCP}}$, one line per chunk count,",
                                       ##       " facetted by $N$ and burn-in chain count (log-scale $y$-axis).}"),
                                       paste0("\\caption{BayesMVP burn-in, ", device_full_name[[device_selector]],
                                             ": seconds per burn-in iteration against $N_{\\text{WCP}}$ for each $N$ and burn-in chain count,",
                                             " for $N_{\\text{chunks}} = 1$ at $N_{\\text{WCP}} = 1$ and for the best measured $N_{\\text{chunks}}$",
                                             " at each $N_{\\text{WCP}}$ (log-scale axes).}"),
                                       paste0("\\label{", figure_label, "}"),
                                       "\\end{figure}",
                                       paste0("\\input{Files/Generated/", bundle_name, "/tables/", table_file_name, "}"))

        }
        ##
        writeLines(text = section_tex_lines, con = file.path(bundle_dir, "section.tex"))
        message(paste0("\033[36m", "fn_paper1_burnin_write_bundle: wrote ", file.path(bundle_dir, "section.tex"), "\033[0m"))
        ##
        if (!is.null(x = manuscript_dir)) {

              destination <-  file.path(manuscript_dir, "Files", "Generated", bundle_name)
              dir.create(path = destination, recursive = TRUE, showWarnings = FALSE)
              for (relative_file in list.files(path = bundle_dir, recursive = TRUE)) {

                    dir.create(path = dirname(path = file.path(destination, relative_file)), recursive = TRUE, showWarnings = FALSE)
                    copy_succeeded <-  file.copy( from      = file.path(bundle_dir, relative_file),
                                                  to        = file.path(destination, relative_file),
                                                  overwrite = TRUE)
                    if (!copy_succeeded) {
                          stop("fn_paper1_burnin_write_bundle: could not copy ", relative_file, " to ", destination)
                    }

              }
              message(paste0("\033[36m", "fn_paper1_burnin_write_bundle: copied bundle to ", destination, "\033[0m"))

        }
        ##
        return(invisible(list(bundle_dir = bundle_dir, manuscript_dir = manuscript_dir)))

}






















