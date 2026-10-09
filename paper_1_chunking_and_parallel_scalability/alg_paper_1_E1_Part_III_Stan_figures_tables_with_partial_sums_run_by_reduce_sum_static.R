##
## ======================================================================================================================================
## make_E1_Part_I_Stan_twin_outputs.R
## (this copy: alg_paper_1_E1_Part_III_Stan_figures_tables_with_partial_sums_run_by_reduce_sum_static.R, which gives the
##  N_chunks of AD_Stan_tape_chunked as the number of partial sums which reduce_sum_static() actually ran; run with the
##  arguments <output folder> [requested] - "requested" keeps the requested N_chunks, to check that the earlier outputs are
##  reproduced)
##
## Stan (tape-chunked, chunking-only) twins of the Paper 1 E1 Part I outputs (NicoStan+BayesMVP, chunking only, one thread per chain):
##   - table_ps1_best_n_chunks_Stan_HPC.tex, table_ps1_best_n_chunks_Stan_Laptop.tex
##   - table_ps1_smt_threshold_Stan.tex, table_ps1_smt_best_n_chunks_Stan.tex
##   - Figure_N_chunks_pilot_study_plot_1_n_threads_SMT_vs_no_SMT_cache_lines_Stan.png
##   - Figure_N_chunks_pilot_study_plot_3_both_devices_cache_lines_Stan.png
##   - numbers.csv
##
## Arm: AD_Stan_tape_chunked (one thread per chain; chains = PS2 chain counts), with the standard Stan model AD_Stan as its
## N_chunks = 1 point. Reads the saved results only (no sampling); runs are filtered to the settings currently written in the runner.
##
## ---- Paths -----------------------------------------------------------------------------------------------------------------------------
##
## Input: the presentation views saved by the latest manuscript export (manuscript_outputs_final_both_devices_one_thread_reference,
## 2026-09-27), built from paper_1_computational_outputs/HPC and /Laptop results.rds (md5 5a2b1886... and a79d77fe..., as recorded in
## its report_manifest.rds). HPC/results.rds itself is being rewritten by the running grid-extension study, so it is not read here.
## A copy of the views and manifest is kept in code/input_snapshot.
##
## Re-run (part3_stan_v2): input switched to manuscript_outputs_final_both_devices_extended_chunk_grid/data/presentation_views.rds
## (2026-09-29 export after the extended Stan and BayesMVP chunk grids were run on both machines; laptop part from Laptop_merged_2026_09_29).
## code/input_snapshot still holds the pre-extension input for reference.
##
{
      paper1_analysis_dir <-  path.expand("~/Documents/Work/PhD_work/Alg_paper_analysis/paper_1_chunking_and_parallel_scalability")
      ##
      paper1_runner_file <-  file.path(paper1_analysis_dir, "alg_paper_1_chunking_WCP_par_scaling.R")
      ##
      paper1_export_dir <-  file.path(paper1_analysis_dir, "paper_1_computational_outputs", "manuscript_outputs_final_both_devices_narrow_WCP_2026_10_03")
      ##
      ##
      # Stan_twin_code_dir <-  "/home/enzocerullo/Documents/Work/PhD_work/Alg_papers_LaTeX/paper_1_v46_2026_10_03/verify3/P3a_generator/Stan_code"
      Stan_twin_code_dir <-  file.path(paper1_analysis_dir, "E1_Part_III_Stan_code")   ## the same three files, copied here
      ##
      # Stan_twin_output_dir <-  "/home/enzocerullo/Documents/Work/PhD_work/Alg_papers_LaTeX/paper_1_v46_2026_10_03/verify3/P3a_generator/output"
      Stan_twin_output_dir <-  commandArgs(trailingOnly = TRUE)[1]
      ##
      show_requested_N_chunks <-  identical(x = commandArgs(trailingOnly = TRUE)[2], y = "requested")
      ##
      dir.create(Stan_twin_output_dir, recursive = TRUE, showWarnings = FALSE)
}
##
## ---- Saved presentation views, report functions (copied) and the templates (adapted copy) -----------------------------------------------
##
{
      source(file.path(Stan_twin_code_dir, "R_fns_alg_paper_1_figures_tables.R"), local = TRUE)
      ##
      fn_paper1_report_dependencies()
      ##
      private_report_environment <-  new.env(parent = environment())
      ##
      sys.source(file.path(Stan_twin_code_dir, "R_fns_alg_paper_1_chunking_WCP_par_scaling.R"), envir = private_report_environment)
      ##
      sys.source(file.path(Stan_twin_code_dir, "R_fns_alg_paper_1_presentation_templates_Stan_twin.R"), envir = private_report_environment)
      ##
      templates <-  private_report_environment$fn_paper1_presentation_templates()
      ##
      presentation_views_file <-  file.path(paper1_export_dir, "data", "presentation_views.rds")
      ##
      message(paste0("\033[36m", "Reading ", presentation_views_file, " (md5 ", unname(tools::md5sum(presentation_views_file)), ")", "\033[0m"))
      ##
      presentation_views <-  readRDS(presentation_views_file)
      ##
      report_manifest <-  readRDS(file.path(paper1_export_dir, "report_manifest.rds"))
}
##
## ---- Keep only configurations that match the settings currently written in the runner (read at run time) -----------------------------
##
{
      all_configurations <-  presentation_views$configurations
      ##
      configuration_key_columns <-  c("device", "algorithm", "N", "num_chunks", "n_threads", "n_chains", "threads_per_chain", "n_iter")
      ##
      configuration_keys <-  function(rows) do.call(paste, c(lapply(rows[, configuration_key_columns, drop = FALSE], as.character), list(sep = "|")))
      ##
      current_runner_grid_keys <-  character()
      current_runner_n_threads_vec <-  list()
      ##
      for (device_name in unique(all_configurations$device)) {

            current_runner_environment <-  fn_paper1_runner_current_settings(runner_file = paper1_runner_file, device = device_name)
            ##
            current_runner_settings <-  get(x = "paper1_settings", envir = current_runner_environment)
            ##
            current_runner_grid <-  get(x = "fn_paper1_benchmark_grid", envir = current_runner_environment)(settings = current_runner_settings)
            current_runner_grid$device <-  device_name
            ##
            current_runner_grid_keys <-  c(current_runner_grid_keys, configuration_keys(current_runner_grid))
            ##
            ## Stan chunking-only configurations in the current runner grid (for the list of configurations not yet measured):
            current_runner_Stan_chunking_grid <-  unique(current_runner_grid[current_runner_grid$algorithm %in% c("AD_Stan", "AD_Stan_tape_chunked"),
                                                                             configuration_key_columns, drop = FALSE])
            current_runner_Stan_chunking_grid_all_devices <-  if (exists("current_runner_Stan_chunking_grid_all_devices"))
                rbind(current_runner_Stan_chunking_grid_all_devices, current_runner_Stan_chunking_grid) else current_runner_Stan_chunking_grid
            current_runner_n_threads_vec[[device_name]] <-  current_runner_settings$n_threads_vec

      }
      ##
      configuration_matches_current_runner <-  configuration_keys(all_configurations) %in% current_runner_grid_keys
      ##
      message(paste0("\033[36m", sum(configuration_matches_current_runner), " of ", length(configuration_matches_current_runner),
                     " saved configurations match the current runner settings.", "\033[0m"))
      ##
      all_configurations <-  all_configurations[configuration_matches_current_runner, , drop = FALSE]
}
##
## ---- Stan chunking-only view: AD_Stan_tape_chunked plus AD_Stan as the N_chunks = 1 point ---------------------------------------------
##
{
      ## (same definition of the main thread budgets as the ps1 view in fn_paper1_presentation_data)
      configuration_on_main_thread_budget <-  vapply( seq_len(nrow(all_configurations)),
                                                      function(i) all_configurations$n_threads[i] %in%
                                                          current_runner_n_threads_vec[[all_configurations$device[i]]],
                                                      logical(1))
      ##
      configuration_on_main_thread_budget <-  configuration_on_main_thread_budget & all_configurations$benchmark_role == "main_scaling"
      ##
      ps1_Stan <-  all_configurations[ all_configurations$algorithm %in% c("AD_Stan", "AD_Stan_tape_chunked") &
                                       configuration_on_main_thread_budget, , drop = FALSE]
      ##
      if (any(ps1_Stan$algorithm == "AD_Stan" & ps1_Stan$num_chunks != 1)) stop("AD_Stan rows with num_chunks other than 1.")
      if (any(ps1_Stan$threads_per_chain != 1) || any(ps1_Stan$n_chains != ps1_Stan$n_threads)) stop("Stan chunking-only rows must use one thread per chain.")
      if (anyDuplicated(ps1_Stan[, c("device", "N", "n_threads", "num_chunks")])) stop("Duplicate Stan chunking-only configurations.")
      ##
      ## N_chunks of AD_Stan_tape_chunked as the number of partial sums which reduce_sum_static() actually ran:
      if (!show_requested_N_chunks) {
            source(file.path(paper1_analysis_dir, "R_fn_number_of_partial_sums_run_by_reduce_sum_static.R"), local = TRUE)
            ps1_Stan$num_chunks_requested <-  ps1_Stan$num_chunks
            ps1_Stan$num_chunks <-  fn_number_of_partial_sums_run_by_reduce_sum_static( N_units            = ps1_Stan$N,
                                                                                        N_chunks_requested = ps1_Stan$num_chunks)
            # if (anyDuplicated(ps1_Stan[, c("device", "N", "n_threads", "num_chunks")])) stop("Two requested N_chunks ran as the same partial sums.")
            ##
            ## Two requested N_chunks which ran as the same partial sums are the same configuration timed twice (on the
            ## local-HPC with one thread, N_chunks = 6 and 8 at N = 500 and 2500, both 8 partial sums): they are shown as one
            ## configuration, with the mean of their timings.
            same_partial_sums_key <-  paste(ps1_Stan$device, ps1_Stan$N, ps1_Stan$n_threads, ps1_Stan$num_chunks, sep = "|")
            for (key in unique(x = same_partial_sums_key[duplicated(x = same_partial_sums_key)])) {
                  rows_with_key <-  which(same_partial_sums_key == key)
                  for (timing_column in intersect(x = c("time_mean", "time_median", "chain_rate", "total_iter_per_sec",
                                                        "seconds_per_iteration_mean", "seconds_per_iteration_median"),
                                                  y = names(x = ps1_Stan))) {
                        ps1_Stan[[timing_column]][rows_with_key[1]] <-  mean(x = ps1_Stan[[timing_column]][rows_with_key])
                  }
                  message(paste0("\033[36m", "Same partial sums (", key, "): requested N_chunks ",
                                 paste(ps1_Stan$num_chunks_requested[rows_with_key], collapse = " and "),
                                 " shown as one configuration (mean timing).", "\033[0m"))
            }
            ps1_Stan <-  ps1_Stan[!duplicated(x = same_partial_sums_key), , drop = FALSE]
      }
      ##
      ps1_Stan$N_chunks_num <-  ps1_Stan$num_chunks
      ps1_Stan$N_chunks <-  factor(ps1_Stan$num_chunks, levels = sort(unique(ps1_Stan$num_chunks)))
      ps1_Stan$N_threads <-  factor(ps1_Stan$n_threads, levels = sort(unique(ps1_Stan$n_threads)))
      ps1_Stan$N_label <-  factor(as.character(ps1_Stan$N_label),
                                  levels = paste0("N = ", fn_paper1_format_number_commas_from_10000(sort(unique(ps1_Stan$N)))))
      ps1_Stan <-  ps1_Stan[order(ps1_Stan$device, ps1_Stan$N, ps1_Stan$n_threads, ps1_Stan$num_chunks), , drop = FALSE]
      ##
      fn_paper1_write_csv(ps1_Stan, file.path(Stan_twin_output_dir, "ps1_Stan_configurations_used.csv"))
      ##
      ## Configurations in the current runner grid with no saved measurement in the input (e.g. the extended Stan chunk grid):
      Stan_chunking_configurations_not_measured <-  current_runner_Stan_chunking_grid_all_devices[
          !configuration_keys(current_runner_Stan_chunking_grid_all_devices) %in% configuration_keys(ps1_Stan), , drop = FALSE]
      Stan_chunking_configurations_not_measured <-  Stan_chunking_configurations_not_measured[
          order(Stan_chunking_configurations_not_measured$device, Stan_chunking_configurations_not_measured$algorithm,
                Stan_chunking_configurations_not_measured$N, Stan_chunking_configurations_not_measured$n_threads,
                Stan_chunking_configurations_not_measured$num_chunks), , drop = FALSE]
      utils::write.csv(Stan_chunking_configurations_not_measured,
                       file.path(Stan_twin_output_dir, "Stan_chunking_configurations_in_current_runner_grid_not_measured.csv"), row.names = FALSE)
      ##
      message(paste0("\033[36m", nrow(Stan_chunking_configurations_not_measured),
                     " Stan chunking-only configurations in the current runner grid have no saved measurement in the input.", "\033[0m"))
      ##
      message(paste0("\033[36m", "Stan chunking-only configurations: ", nrow(ps1_Stan),
                     " (AD_Stan: ", sum(ps1_Stan$algorithm == "AD_Stan"),
                     ", AD_Stan_tape_chunked: ", sum(ps1_Stan$algorithm == "AD_Stan_tape_chunked"), ")", "\033[0m"))
}
##
## ---- Figures (same styling as the BayesMVP figures; 19,152 bytes per individual for the cache lines) -----------------------------------
##
{
      grDevices::pdf(NULL)
      ##
      templates$R_fn_plot_ps1_N_chunks_ggplot_SMT_combined( ps1_Stan,
                                                            output_path = Stan_twin_output_dir,
                                                            show_cache_lines = TRUE,
                                                            bytes_per_row = templates$paper1_bytes_per_row_Stan,
                                                            thresholds_file = file.path(Stan_twin_output_dir,
                                                                'Figure_N_chunks_pilot_study_plot_1_n_threads_SMT_vs_no_SMT_cache_lines_Stan_thresholds.csv'),
                                                            output_file_suffix = "_Stan")
      ##
      templates$R_fn_plot_ps1_efficiency_combined( ps1_Stan,
                                                   output_path = Stan_twin_output_dir,
                                                   show_cache_lines = TRUE,
                                                   bytes_per_row = templates$paper1_bytes_per_row_Stan,
                                                   thresholds_file = file.path(Stan_twin_output_dir,
                                                       'Figure_N_chunks_pilot_study_plot_3_both_devices_cache_lines_Stan_thresholds.csv'),
                                                   output_file_suffix = "_Stan")
      ##
      grDevices::dev.off()
}
##
## ---- Best N_chunks by device, N and N_threads (efficiency = chains / time; ties choose the smaller count) ------------------------------
##
{
      best_chunks_Stan <-  as.data.frame(templates$get_best_chunks(ps1_Stan))
      ##
      tape_chunked_rows <-  ps1_Stan[ps1_Stan$algorithm == "AD_Stan_tape_chunked", , drop = FALSE]
      ##
      cell_keys <-  function(rows) paste(rows$device, rows$N, rows$n_threads, sep = "|")
      ##
      largest_tested_chunks_by_cell <-  tapply(tape_chunked_rows$num_chunks, cell_keys(tape_chunked_rows), max)
      tested_chunks_by_cell <-  tapply(ps1_Stan$num_chunks, cell_keys(ps1_Stan), function(chunks) paste(sort(chunks), collapse = ";"))
      ##
      best_chunks_Stan$largest_tested_N_chunks <-  unname(largest_tested_chunks_by_cell[cell_keys(best_chunks_Stan)])
      best_chunks_Stan$best_at_largest_tested_chunks <-  best_chunks_Stan$best_N_chunks == best_chunks_Stan$largest_tested_N_chunks
      ##
      thread_vec_by_device <-  current_runner_n_threads_vec
      ##
      device_caption_names <-  c(HPC = "the HPC", Laptop = "the laptop")
      ##
      for (device_name in names(thread_vec_by_device)) {

            templates$make_ps1_best_chunks_table_tex_manuscript_layout(
                best_chunks_df = best_chunks_Stan[best_chunks_Stan$device == device_name, , drop = FALSE],
                thread_vec = thread_vec_by_device[[device_name]],
                caption = paste0( "Optimal $N_{\\text{chunks}}$ for tape-chunked Stan (\\texttt{AD\\_Stan\\_tape\\_chunked}) for ",
                                  device_caption_names[[device_name]],
                                  " at each $N_{\\text{threads}}$ and dataset size ",
                                  "($N_{\\text{chunks}} = 1$ is the standard Stan model, \\texttt{AD\\_Stan}; ",
                                  "$^{\\dagger}$ = the largest $N_{\\text{chunks}}$ tested for that $N$)."),
                label = paste0("table:ps1_best_n_chunks_Stan_", device_name),
                file_path = file.path(Stan_twin_output_dir, paste0("table_ps1_best_n_chunks_Stan_", device_name, ".tex")))

      }
}
##
## ---- SMT tables (HPC 96 vs 180 threads; laptop 8 vs 16) ----------------------------------------------------------------------------------
##
{
      SMT_thread_pairs <-  list(HPC = c(96, 180), Laptop = c(8, 16))
      ##
      SMT_benefit_Stan <-  lapply( names(SMT_thread_pairs),
                                   function(device_name) templates$get_smt_benefit( ps1_Stan[ps1_Stan$device == device_name, , drop = FALSE],
                                                                                    SMT_thread_pairs[[device_name]][1],
                                                                                    SMT_thread_pairs[[device_name]][2]))
      names(SMT_benefit_Stan) <-  names(SMT_thread_pairs)
      ##
      SMT_threshold_Stan <-  as.data.frame(dplyr::bind_rows(lapply(names(SMT_benefit_Stan),
                                                                   function(device_name) templates$get_smt_threshold(SMT_benefit_Stan[[device_name]], device_name))))
      ##
      SMT_best_Stan <-  as.data.frame(dplyr::bind_rows(lapply(names(SMT_benefit_Stan),
                                                              function(device_name) templates$get_smt_best(SMT_benefit_Stan[[device_name]], device_name))))
      ##
      Stan_caption_prefix <-  "Tape-chunked Stan (\\texttt{AD\\_Stan\\_tape\\_chunked}; $N_{\\text{chunks}} = 1$ is the standard Stan model, \\texttt{AD\\_Stan}): "
      ##
      templates$make_ps1_smt_table_tex_manuscript_layout(
          SMT_threshold_Stan,
          "Min.\\ $N_{\\text{chunks}}$",
          paste0( Stan_caption_prefix,
                  "the minimum $N_{\\text{chunks}}$ required for SMT to improve efficiency, \n",
                  "compared to physical cores alone (96 vs 180 threads on local-HPC; 8 vs 16 on laptop),\n",
                  "with the corresponding efficiency gain from using SMT."),
          "table:ps1_smt_threshold_Stan",
          file.path(Stan_twin_output_dir, "table_ps1_smt_threshold_Stan.tex"))
      ##
      templates$make_ps1_smt_table_tex_manuscript_layout(
          SMT_best_Stan,
          "Best\\ $N_{\\text{chunks}}$",
          paste0( Stan_caption_prefix,
                  "the $N_{\\text{chunks}}$ values which maximise efficiency with SMT for each $N$,\n",
                  "comparing efficiency with and without SMT at the same $N_{\\text{chunks}}$\n",
                  "(96 vs 180 threads on HPC; 8 vs 16 on laptop),\n",
                  "with the corresponding efficiency gain from using SMT."),
          "table:ps1_smt_best_n_chunks_Stan",
          file.path(Stan_twin_output_dir, "table_ps1_smt_best_n_chunks_Stan.tex"))
}
##
## ---- numbers.csv: per device, N and N_threads, plus the SMT gains -------------------------------------------------------------------------
##
{
      best_chunk_records <-  do.call(rbind, lapply(split(ps1_Stan, cell_keys(ps1_Stan)), function(cell) {

            cell <-  cell[order(-cell$Efficiency, cell$num_chunks), , drop = FALSE]     ## ties: smaller count first
            ##
            no_chunking_row <-  cell[cell$algorithm == "AD_Stan", , drop = FALSE]
            tape_chunked_cell <-  cell[cell$algorithm == "AD_Stan_tape_chunked", , drop = FALSE]
            ##
            time_no_chunking <-  if (nrow(no_chunking_row)) no_chunking_row$time_mean[1] else NA_real_
            ##
            data.frame( record_type                                        = "best_N_chunks",
                        device                                             = cell$device[1],
                        N                                                  = cell$N[1],
                        n_threads                                          = cell$n_threads[1],
                        n_chains                                           = cell$n_chains[1],
                        n_iter_long_run                                    = cell$n_iter[1],
                        tested_N_chunks                                    = paste(sort(cell$num_chunks), collapse = ";"),
                        time_no_chunking_AD_Stan_seconds                   = time_no_chunking,
                        time_no_chunking_AD_Stan_SD_seconds                = if (nrow(no_chunking_row)) no_chunking_row$time_sd[1] else NA_real_,
                        best_N_chunks                                      = cell$num_chunks[1],
                        best_is_AD_Stan_no_chunking                        = cell$algorithm[1] == "AD_Stan",
                        time_best_seconds                                  = cell$time_mean[1],
                        time_best_SD_seconds                               = cell$time_sd[1],
                        n_repeats_best                                     = cell$n_repeats[1],
                        speed_up_best_vs_no_chunking                       = time_no_chunking / cell$time_mean[1],
                        second_best_N_chunks                               = if (nrow(cell) > 1) cell$num_chunks[2] else NA_real_,
                        time_second_best_seconds                           = if (nrow(cell) > 1) cell$time_mean[2] else NA_real_,
                        best_tape_chunked_N_chunks                         = if (nrow(tape_chunked_cell)) tape_chunked_cell$num_chunks[1] else NA_real_,
                        time_best_tape_chunked_seconds                     = if (nrow(tape_chunked_cell)) tape_chunked_cell$time_mean[1] else NA_real_,
                        speed_up_best_tape_chunked_vs_no_chunking          = if (nrow(tape_chunked_cell)) time_no_chunking / tape_chunked_cell$time_mean[1] else NA_real_,
                        largest_tested_N_chunks                            = if (nrow(tape_chunked_cell)) max(tape_chunked_cell$num_chunks) else NA_real_,
                        best_at_largest_tested_N_chunks                    = if (nrow(tape_chunked_cell)) cell$num_chunks[1] == max(tape_chunked_cell$num_chunks) else NA,
                        best_tape_chunked_at_largest_tested_N_chunks       = if (nrow(tape_chunked_cell)) tape_chunked_cell$num_chunks[1] == max(tape_chunked_cell$num_chunks) else NA,
                        efficiency_best_chains_per_second                  = cell$Efficiency[1],
                        max_divergences_best                               = cell$max_divergences[1],
                        stringsAsFactors = FALSE)

      }))
      ##
      best_chunk_records <-  best_chunk_records[order(best_chunk_records$device, best_chunk_records$N, best_chunk_records$n_threads), , drop = FALSE]
      ##
      SMT_records <-  dplyr::bind_rows( dplyr::mutate(SMT_threshold_Stan, record_type = "SMT_threshold_min_N_chunks_with_gain"),
                                        dplyr::mutate(SMT_best_Stan,      record_type = "SMT_best_N_chunks_with_SMT"))
      SMT_records <-  data.frame( record_type                              = SMT_records$record_type,
                                  device                                   = SMT_records$device,
                                  N                                        = SMT_records$N,
                                  n_threads_without_SMT                    = unname(vapply(SMT_records$device, function(d) SMT_thread_pairs[[d]][1], numeric(1))),
                                  n_threads_with_SMT                       = unname(vapply(SMT_records$device, function(d) SMT_thread_pairs[[d]][2], numeric(1))),
                                  SMT_N_chunks                             = SMT_records$N_chunks,
                                  efficiency_without_SMT_chains_per_second = SMT_records$eff_physical,
                                  efficiency_with_SMT_chains_per_second    = SMT_records$eff_smt,
                                  SMT_gain_percent                         = SMT_records$gain_pct,
                                  stringsAsFactors = FALSE)
      ##
      numbers_Stan <-  dplyr::bind_rows(best_chunk_records, SMT_records)
      ##
      utils::write.csv(numbers_Stan, file.path(Stan_twin_output_dir, "numbers.csv"), row.names = FALSE)
      ##
      ## Full SMT comparison at every chunk count (for checking the two SMT tables):
      utils::write.csv( dplyr::bind_rows(lapply(names(SMT_benefit_Stan), function(device_name) dplyr::mutate(SMT_benefit_Stan[[device_name]], device = device_name))),
                        file.path(Stan_twin_output_dir, "SMT_comparison_every_N_chunks_Stan.csv"),
                        row.names = FALSE)
      ##
      message(paste0("\033[32m", "Stan twin outputs written to: ", Stan_twin_output_dir, "\033[0m"))
}
























