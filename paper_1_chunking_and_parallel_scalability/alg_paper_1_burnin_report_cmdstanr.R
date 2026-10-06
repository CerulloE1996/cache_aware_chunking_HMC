##
## ==============================================================================================================
## alg_paper_1_burnin_report_cmdstanr.R
##
## Runner for the burn-in figures of the Stan via cmdstanr burn-in study (E1, Stan via cmdstanr), one per device:
## was: ## seconds per leapfrog step of the slowest chain against N_threads/chain, with the burn-in N_chains of
## was: ## the shootout (4), for AD_Stan, AD_Stan_tape_chunked, AD_Stan_WCP and AD_Stan_WCP_chunking, in the
## was: ## visual design of the NicoStan+BayesMVP and Stan via NicoStan burn-in figures.
## seconds per burn-in iteration of the slowest chain against N_threads/chain, with the burn-in N_chains of the
## shootout (4), for AD_Stan, AD_Stan_tape_chunked, AD_Stan_WCP and AD_Stan_WCP_chunking, in the visual design of
## the NicoStan+BayesMVP and Stan via NicoStan burn-in figures.
##
## Reads the saved runs of ps_1_burnin_Stan_cmdstanr_shootout_4_chains.R (never re-run here), keeping only the
## runs that match that runner's current settings (read from the runner itself), and writes:
##   - report_outputs/Burnin_chunks_WCP/figures/figure_burnin_Stan_cmdstanr_HPC.png and
##     report_outputs/Burnin_chunks_WCP/figures/figure_burnin_Stan_cmdstanr_Laptop.png;
##   - report_outputs/Burnin_chunks_WCP/values/burnin_Stan_cmdstanr_figure_values.csv;
## then checks the best configuration per device and N against burnin_Stan_cmdstanr_shootout_summary.csv
## (written by ps_1_burnin_Stan_cmdstanr_shootout_summary.R).
## ==============================================================================================================
##
##
##
## ---- Settings: -----------------------------------------------------------------------------------------------
##
paper_1_dir <-  file.path( Sys.getenv("HOME"), "Documents", "Work", "PhD_work", "Alg_paper_analysis",
                           "paper_1_chunking_and_parallel_scalability")
##
shootout_runner_file_path <-  file.path(paper_1_dir, "ps_1_burnin_Stan_cmdstanr_shootout_4_chains.R")
##
devices <-  c("HPC", "Laptop")
##
figure_dir <-  file.path(paper_1_dir, "report_outputs", "Burnin_chunks_WCP", "figures")
values_dir <-  file.path(paper_1_dir, "report_outputs", "Burnin_chunks_WCP", "values")
##
values_file_path  <-  file.path(values_dir, "burnin_Stan_cmdstanr_figure_values.csv")
summary_file_path <-  file.path(values_dir, "burnin_Stan_cmdstanr_shootout_summary.csv")
##
## ---- the same width and dpi as the other burn-in figures (figure_burnin_HPC.png and figure_burnin_Laptop.png:
##      10 in wide at 150 dpi, with 3 in per row of panels plus 1.6 in for the strips, axis titles and legend);
##      the shootout has one row of panels (burn-in N_chains = 4):
##
figure_width_inches  <-  10
figure_height_inches <-  3 * 1 + 1.6
figure_dpi           <-  150
##
## ---- "total" fold label (AD_Stan relative to the best measured configuration) above the "WCP" fold label,
##      as in the NicoStan+BayesMVP burn-in figures (FALSE = "WCP" fold label only):
##
show_total_fold_label <-  TRUE
##
##
##
## ---- Functions (definitions only; no sampling, no compilation): -----------------------------------------------
##
source(file.path(paper_1_dir, "R_fns_alg_paper_1_burnin_report.R"))
source(file.path(paper_1_dir, "R_fns_alg_paper_1_burnin_report_cmdstanr.R"))
##
##
##
## ---- One figure per device: ---------------------------------------------------------------------------------
##
figure_values_list <-  list()
##
for (device in devices) {

      runner_settings <-  fn_paper1_burnin_cmdstanr_runner_settings( runner_file_path = shootout_runner_file_path,
                                                                     device           = device)
      ##
      runs <-  fn_paper1_burnin_cmdstanr_read_runs(runner_settings = runner_settings)
      ##
      configuration_summary <-  fn_paper1_burnin_cmdstanr_summarise_runs(runs = runs)
      ##
      figure_values <-  fn_paper1_burnin_cmdstanr_figure_values(configuration_summary = configuration_summary)
      ##
      fn_paper1_burnin_cmdstanr_figure( configuration_summary = configuration_summary,
                                        figure_values         = figure_values,
                                        device                = device,
                                        n_chains_burnin       = runner_settings$n_chains_burnin,
                                        file_path             = file.path( figure_dir,
                                                                           paste0( "figure_burnin_Stan_cmdstanr_",
                                                                                   device, ".png")),
                                        width_inches          = figure_width_inches,
                                        height_inches         = figure_height_inches,
                                        dpi                   = figure_dpi,
                                        show_total_fold_label = show_total_fold_label)
      ##
      figure_values_list[[device]] <-  figure_values

}
##
##
##
## ---- Values file (per device and N): -------------------------------------------------------------------------
##
figure_values_all <-  do.call(rbind, figure_values_list)
rownames(figure_values_all) <-  NULL
##
dir.create(path = values_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(x = figure_values_all, file = values_file_path, row.names = FALSE)
##
message(NicoStan::colourise(text = paste0("Wrote ", values_file_path), fg = "green"))
##
##
##
## was: ## ---- Check against the shootout summary (best configuration and its seconds per leapfrog step): ------
## ---- Check against the shootout summary (best configuration and its seconds per burn-in iteration): ----------
##
if (!file.exists(summary_file_path)) {

      message(NicoStan::colourise( text = paste0("No shootout summary to check against: ", summary_file_path),
                                   fg   = "red"))

} else {

      shootout_summary <-  utils::read.csv(summary_file_path)
      ##
      for (row_index in seq_len(nrow(figure_values_all))) {

            figure_row  <-  figure_values_all[row_index, ]
            summary_row <-  shootout_summary[ shootout_summary$device == figure_row$device &
                                              shootout_summary$N == figure_row$N, , drop = FALSE]
            ##
            check_start <-  paste0("Check, ", figure_row$device, ", N = ", figure_row$N, ": ")
            ##
            if (nrow(summary_row) != 1) {
                  message(NicoStan::colourise( text = paste0(check_start, "not in the shootout summary."),
                                               fg   = "red"))
                  next
            }
            ##
            ## ---- the summary stores N_chunks = 0 for AD_Stan (NA here):
            ##
            figure_best_num_chunks <-  if (is.na(figure_row$best_num_chunks)) 0 else figure_row$best_num_chunks
            ##
            same_best_configuration <-  figure_row$best_algorithm == summary_row$best_algorithm &&
                                        figure_row$best_n_threads_per_chain ==
                                            summary_row$best_n_threads_per_chain &&
                                        figure_best_num_chunks == summary_row$best_num_chunks
            ##
            ## was: relative_difference_sec_per_step <-  abs( figure_row$best_sec_per_step /
            ## was:                                               summary_row$best_sec_per_step - 1)
            relative_difference_sec_per_burnin_iteration_slowest_chain <-
                  abs( figure_row$best_sec_per_burnin_iteration_slowest_chain /
                           summary_row$best_sec_per_burnin_iteration_slowest_chain - 1)
            relative_difference_total_fold   <-  abs( figure_row$total_fold /
                                                          summary_row$speed_up_vs_AD_Stan_4x1 - 1)
            ##
            ## was: all_match <-  same_best_configuration &&
            ## was:               relative_difference_sec_per_step < 1e-9 &&
            all_match <-  same_best_configuration &&
                          relative_difference_sec_per_burnin_iteration_slowest_chain < 1e-9 &&
                          relative_difference_total_fold < 1e-9
            ##
            message(NicoStan::colourise( text = paste0( check_start,
                                                        if (all_match) "matches" else "MISMATCH", " (figure ",
                                                        figure_row$best_algorithm, " ",
                                                        figure_row$best_n_threads_per_chain, "/",
                                                        figure_row$best_num_chunks, ", summary ",
                                                        summary_row$best_algorithm, " ",
                                                        summary_row$best_n_threads_per_chain, "/",
                                                        summary_row$best_num_chunks, ")"),
                                         fg   = if (all_match) "green" else "red"))

      }

}























