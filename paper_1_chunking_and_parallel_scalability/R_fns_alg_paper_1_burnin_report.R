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
              fn_paper1_burnin_figure( rows      = rows,
                                       device    = device_selector,
                                       file_path = file.path(figure_dir, figure_file_name))
              ##
              section_tex_lines <-  c( section_tex_lines,
                                       "\\begin{figure}[H]",
                                       "\\centering",
                                       paste0("\\includegraphics[width=\\textwidth,height=0.8\\textheight,keepaspectratio]{",
                                             "Files/Generated/", bundle_name, "/figures/", figure_file_name, "}"),
                                       paste0("\\caption{BayesMVP burn-in, ", device_full_name[[device_selector]],
                                             ": seconds per burn-in iteration against $N_{\\text{WCP}}$, one line per chunk count,",
                                             " facetted by $N$ and burn-in chain count (log-scale $y$-axis).}"),
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






















