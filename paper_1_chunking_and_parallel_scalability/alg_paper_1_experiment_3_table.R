##
## =====================================================================================================================================
## alg_paper_1_experiment_3_table.R
##
## Paper 1, Experiment 3 (absolute efficiency): builds the BayesMVP vs Stan (NUTS) vs Mplus (PX-Gibbs) time-to-target-ESS
## table from Paper 2's pilot studies (ps5 = Mplus, ps6 = Stan, ps7 = BayesMVP) and writes the Exp3_absolute_efficiency
## bundle. Every setting below is explicit - no defaults live in R_fns_alg_paper_1_experiment_3.R - so re-running this
## after the ps1/ps7 rerun means editing the paths below and nothing else.
##
## ---- Explicit locations -----------------------------------------------------------------------------------------------------------
##
algorithm_study_dir <-  path.expand("~/Documents/Work/PhD_work/Alg_paper_analysis")
##
helper_dir <-  file.path(algorithm_study_dir, "paper_1_chunking_and_parallel_scalability")
##
pilot_study_dir <-  file.path(algorithm_study_dir, "1_appendix_pilot_studies")
##
ps7_functions_dir <-  file.path(pilot_study_dir, "ps_7_basic_MCMC_settings_BayesMVP", "functions")
##
##
## ---- Explicit per-device source files/directories -----------------------------------------------------------------------------------
##
## Specify every device with available data; leave a device's entries as character(0) / NULL to skip it. The Laptop ps5
## file below is the current N = 500/2500/10000 export; point it at the final file after the ps1/ps7 rerun.
##
ps5_csv_file_path_by_device <-  list(
      Laptop = file.path(pilot_study_dir, "ps_5_basic_MCMC_settings_Mplus", "outputs",
                          paste0("Laptop_ps5_full_results_N500-2500-10000_WCPFALSE_popprevTRUE_n_threads16_n_chains16_",
                                 "n_fb_iter1000_n_thin10pi_Wish14_14pi_p-1.48_0.65_0.65pi_lab-1.5_0.5pi_lft-1_0.5pi_nd2_0.5.csv")),
      HPC    = file.path(pilot_study_dir, "ps_5_basic_MCMC_settings_Mplus", "outputs",
                          "ps5_full_results_N500-2500-10000_WCPFALSE_popprevTRUE_n_threads180_n_chains180_n_fb_iter1000_n_thin10pi_Wish14_14pi_p-1.48_0.65_0.65pi_lab-1.5_0.5pi_lft-1_0.5pi_nd2_0.5.csv"))
##
## ps6_csv_file_paths_by_device <-  list(
##       Laptop = character(0),
##       HPC    = file.path(pilot_study_dir, "ps_6_basic_MCMC_settings_Stan", "outputs", "DGM_3",
##                           paste0("ps6_without_partial_log_full_results_N", c(500, 2500, 10000), ".csv")))
##
## ---- Production PS6 Stan runs (cmdstanr NUTS, tape-chunked reduce_sum model with one thread per chain, 96 chains at N = 500
##      and 64 chains at N = 2500 and 10,000, 250 sampling iterations, burn-in lengths 250 / 500 / 1000, three runs each),
##      one CSV per N built from the 27 ps6_BIN_Stan_N*_run*.RDS files; these carry the measured posterior-summary time.
##
## ps6_csv_file_paths_by_device <-  list(
##       Laptop = character(0),
##       HPC    = file.path(pilot_study_dir, "ps_6_basic_MCMC_settings_Stan", "outputs",
##                           paste0("ps6_BIN_Stan_production_full_results_N", c(500, 2500, 10000), ".csv")))
##
## ---- All binary HPC PS6 production arms, one CSV per N with an arm column (config_id | burn-in model | warm-up chains |
##      sampling chains): the "standard" arm above (27 runs) plus the "asymburn4" arm (27 runs: 4 warm-up chains with the
##      unpartitioned AD_Stan model, then 96 / 64 sampling chains with the tape-chunked model, one thread per chain).
##      With select_fastest_burnin = TRUE the Stan row is the (arm, burn-in length) cell with the lowest mean time to target.
##
ps6_csv_file_paths_by_device <-  list(
      Laptop = character(0),
      HPC    = file.path(pilot_study_dir, "ps_6_basic_MCMC_settings_Stan", "outputs",
                          paste0("ps6_BIN_Stan_production_ALL_ARMS_full_results_N", c(500, 2500, 10000), ".csv")))
##
## ---- Stan burn-in length: as for BayesMVP, keep the burn-in length with the lowest mean time to the target ESS (three runs)
##      at each N; FALSE pools every run in the CSV(s).
##
ps6_select_fastest_burnin_by_device <-  list( Laptop = FALSE,
                                              HPC    = TRUE)
##
ps7_runs_directory_by_device <-  list(
      Laptop = NA_character_,
      HPC    = file.path(pilot_study_dir, "ps_7_basic_MCMC_settings_BayesMVP", "outputs", "DGM_3"))
##
## ---- N values with no ps7 runs under the current settings are read from the old-configuration runs (250 burn-in,
##      100 sampling iterations) kept in the backup folder; only the listed N values are taken from there.
##
ps7_old_configuration_runs_directory_by_device <-  list(
      Laptop = NA_character_,
      HPC    = file.path(pilot_study_dir, "ps_7_basic_MCMC_settings_BayesMVP", "outputs", "DGM_3", "2nd_Sept_2026_backup"))
# ps7_old_configuration_N_values_by_device <-  list( Laptop = numeric(0),
#                                                    HPC    = c(500))
## Paper results use only runs matching the CURRENT runner settings, so no N is taken from the backup folder:
ps7_old_configuration_N_values_by_device <-  list( Laptop = numeric(0),
                                                   HPC    = numeric(0))
##
## ---- Only ps7 runs with the runner's CURRENT number of sampling iterations are used; n_iter is read from the active line of
##      the ps7 runner at run time (never hard-coded here), and matched through the "_it<n_iter>_" file-name token:
##
ps7_runner_file <-  file.path(pilot_study_dir, "ps_7_basic_MCMC_settings_BayesMVP", "ps_7_MCMC_settings_BayesMVP.R")
ps7_runner_n_iter_lines <-  grep( pattern = "^[[:space:]]*n_iter[[:space:]]*=[[:space:]]*[0-9]+[[:space:]]*,",
                                  x = readLines(con = ps7_runner_file),
                                  value = TRUE)
if (length(ps7_runner_n_iter_lines) != 1) {
      stop(paste0("Expected exactly one active 'n_iter = <number>,' line in ", ps7_runner_file, "; found ", length(ps7_runner_n_iter_lines), "."))
}
ps7_current_n_iter <-  as.numeric(sub(pattern = "^[[:space:]]*n_iter[[:space:]]*=[[:space:]]*([0-9]+).*$", replacement = "\\1",
                                      x = ps7_runner_n_iter_lines))
message(paste0("\033[36m", "ps7: using only runs with n_iter = ", ps7_current_n_iter, " (the current setting in ", basename(ps7_runner_file), ")", "\033[0m"))
##
## ---- Target minimum ESS per N (over Se / Sp / prevalence), same target used for every software at that N ---------------------------
##
## Suggested from the ps7 pilot study; update if the rerun changes the target policy.
##
target_min_ESS_by_N <-  c( "500"   = 7000,
                           "2500"  = 2500,
                           "10000" = 1000,
                           "50000" = 1000)
##
## ---- N values to show, left to right, per device (only N values present in that device's sources are usable) -----------------------
##
N_values_by_device <-  list( Laptop = c(500, 2500, 10000),
                             HPC    = c(500, 2500, 10000, 50000))
##
## ---- Text written in the table for any value that is not available (the manuscript's [TBD] placeholder macro);
##      every such cell is also listed in the console, never filled in:
##
exp3_placeholder_text <-  "\\EthreeTBD{}"
##
## ---- Output locations -----------------------------------------------------------------------------------------------------------
##
output_dir <-  file.path(helper_dir, "paper_1_computational_outputs", "Exp3_absolute_efficiency_generated")
##
manuscript_dir <-  NULL   ## set explicitly to the paper_1 LaTeX directory to also copy the bundle there.
##
## ---- Load the helper functions -----------------------------------------------------------------------------------------------------
##
source(file.path(helper_dir, "R_fns_alg_paper_1_experiment_3.R"), local = TRUE)
##
##
## ---- Build one table per device that has any source, then bundle them together ------------------------------------------------------
##
table_results <-  list()
##
for (device_label in names(ps5_csv_file_path_by_device)) {

      ps5_csv_file_path   <-  ps5_csv_file_path_by_device[[device_label]]
      ps6_csv_file_paths  <-  ps6_csv_file_paths_by_device[[device_label]]
      ps7_runs_directory  <-  ps7_runs_directory_by_device[[device_label]]
      ##
      exp3_rows_for_this_device <-  list()
      ##
      if (length(ps5_csv_file_path) == 1 && !is.na(ps5_csv_file_path) && nzchar(ps5_csv_file_path)) {

            message(paste0("\033[36m", "Reading ps5 (Mplus) for device = ", device_label, ": ", ps5_csv_file_path, "\033[0m"))
            exp3_rows_for_this_device$mplus <-  fn_paper1_exp3_rows_from_ps5_csv( ps5_csv_file_path = ps5_csv_file_path,
                                                                                  device_label      = device_label)

      }
      ##
      if (length(ps6_csv_file_paths) > 0 && !anyNA(ps6_csv_file_paths)) {

            existing_ps6_csv_file_paths <-  ps6_csv_file_paths[file.exists(ps6_csv_file_paths)]
            missing_ps6_csv_file_paths  <-  setdiff(ps6_csv_file_paths, existing_ps6_csv_file_paths)
            if (length(missing_ps6_csv_file_paths) > 0) {
                  message(paste0("\033[36m", "Skipping missing ps6 (Stan) file(s) for device = ", device_label, ": ",
                                 paste(missing_ps6_csv_file_paths, collapse = ", "), "\033[0m"))
            }
            if (length(existing_ps6_csv_file_paths) > 0) {
                  message(paste0("\033[36m", "Reading ps6 (Stan) for device = ", device_label, ": ",
                                 paste(existing_ps6_csv_file_paths, collapse = ", "), "\033[0m"))
                  ## exp3_rows_for_this_device$stan <-  fn_paper1_exp3_rows_from_ps6_csv( ps6_csv_file_paths = existing_ps6_csv_file_paths,
                  ##                                                                      device_label       = device_label)
                  exp3_rows_for_this_device$stan <-  fn_paper1_exp3_rows_from_ps6_csv( ps6_csv_file_paths    = existing_ps6_csv_file_paths,
                                                                                       device_label          = device_label,
                                                                                       select_fastest_burnin = isTRUE(ps6_select_fastest_burnin_by_device[[device_label]]),
                                                                                       target_min_ESS_by_N   = target_min_ESS_by_N)
            }

      }
      ##
      if (length(ps7_runs_directory) == 1 && !is.na(ps7_runs_directory) && nzchar(ps7_runs_directory) && dir.exists(ps7_runs_directory)) {

            message(paste0("\033[36m", "Reading ps7 (BayesMVP) for device = ", device_label, ": ", ps7_runs_directory, "\033[0m"))
            exp3_rows_for_this_device$bayesmvp <-  fn_paper1_exp3_rows_from_ps7_directory( ps7_runs_directory       = ps7_runs_directory,
                                                                                            device_label             = device_label,
                                                                                            target_min_ESS_by_N      = target_min_ESS_by_N,
                                                                                            ps7_functions_directory  = ps7_functions_dir,
                                                                                            ps7_file_pattern         = paste0("^ps7_run_.*_it", ps7_current_n_iter, "_"))

      }
      ##
      ps7_old_runs_directory <-  ps7_old_configuration_runs_directory_by_device[[device_label]]
      ps7_old_N_values       <-  ps7_old_configuration_N_values_by_device[[device_label]]
      ##
      if (length(ps7_old_N_values) > 0 && length(ps7_old_runs_directory) == 1 && !is.na(ps7_old_runs_directory) &&
          nzchar(ps7_old_runs_directory) && dir.exists(ps7_old_runs_directory)) {

            message(paste0("\033[36m", "Reading ps7 (BayesMVP, old configuration) for device = ", device_label, ", N = ",
                           paste(ps7_old_N_values, collapse = ", "), ": ", ps7_old_runs_directory, "\033[0m"))
            ps7_old_rows <-  fn_paper1_exp3_rows_from_ps7_directory( ps7_runs_directory       = ps7_old_runs_directory,
                                                                     device_label             = device_label,
                                                                     target_min_ESS_by_N      = target_min_ESS_by_N,
                                                                     ps7_functions_directory  = ps7_functions_dir,
                                                                     ps7_file_pattern         = paste0("^ps7_run_.*_N(", paste(ps7_old_N_values, collapse = "|"), ")_"))
            ps7_old_rows <-  ps7_old_rows[ps7_old_rows$N %in% ps7_old_N_values, , drop = FALSE]
            ps7_old_rows <-  ps7_old_rows[!(ps7_old_rows$N %in% exp3_rows_for_this_device$bayesmvp$N), , drop = FALSE]
            if (nrow(ps7_old_rows) > 0) exp3_rows_for_this_device$bayesmvp_old_configuration <-  ps7_old_rows

      }
      ##
      if (length(exp3_rows_for_this_device) == 0) {
            message(paste0("\033[36m", "No Experiment 3 sources for device = ", device_label, "; skipping.", "\033[0m"))
            next
      }
      ##
      exp3_rows_combined <-  dplyr::bind_rows(exp3_rows_for_this_device)
      ##
      exp3_summary <-  fn_paper1_exp3_time_to_target( exp3_rows           = exp3_rows_combined,
                                                       target_min_ESS_by_N = target_min_ESS_by_N)
      ##
      N_values_present_for_this_device <-  intersect(N_values_by_device[[device_label]], unique(exp3_summary$N))
      if (length(N_values_present_for_this_device) == 0) {
            message(paste0("\033[36m", "No usable N values for device = ", device_label, "; skipping.", "\033[0m"))
            next
      }
      ##
      table_results[[device_label]] <-  fn_paper1_exp3_table_tex(
            exp3_summary      = exp3_summary,
            device_label      = device_label,
            N_values          = N_values_present_for_this_device,
            table_caption     = paste0("Experiment 3: mean time (seconds) to reach the target minimum ESS over Se, Sp and ",
                                       "prevalence, and speed-up relative to BayesMVP, device = ", device_label, "."),
            ## table_label       = paste0("table:exp3_absolute_efficiency_", tolower(device_label)),
            table_label       = paste0("S:table:exp3_absolute_efficiency_", tolower(device_label)),
            placeholder_text  = exp3_placeholder_text,
            output_file_path  = file.path(output_dir, "paper_sections", "Exp3_absolute_efficiency", "tables",
                                          paste0("table_Exp3_absolute_efficiency_", device_label, ".tex")))

}
##
if (length(table_results) == 0) stop("alg_paper_1_experiment_3_table.R: no device produced a table; check the source paths above.")
##
bundle_dir <-  fn_paper1_exp3_write_bundle( table_results  = table_results,
                                            output_dir     = output_dir,
                                            manuscript_dir = manuscript_dir)
##
message(paste0("\033[36m", "Experiment 3 bundle written to: ", bundle_dir, "\033[0m"))






















