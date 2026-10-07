##
## ===============================================================================================================
## alg_paper_1_burnin_table_selected_configurations.R
##
## Table table:paper1_burnin_selected_configurations of Paper 1 (E1 Part VI): for each model (NicoStan+BayesMVP,
## the Stan model via NicoStan, the Stan model via cmdstanr), device and N, the configuration with the LOWEST mean
## seconds per (joint) burn-in iteration across repetitions, over all burn-in chain counts, with its near-ties
## (within 1%, or within run-to-run noise) and its seconds per burn-in iteration (per burn-in iteration of the
## slowest chain for cmdstanr).
##
## Reads the saved burn-in .rds files (never edits them), keeps only the runs that match the CURRENT settings of
## each NicoStan burn-in runner (read from the runner files at run time, as alg_paper_1_burnin_report_Stan.R
## does), and writes into a staging folder (never over report_outputs/Burnin_chunks_WCP/figures or values):
##   - table_burnin_selected_configurations_rows.tex   (the body rows of the table, with the seconds rows),
##   - burnin_selected_configurations_lowest_time.csv  (every selection, its type, grid edges and near-ties).
## ===============================================================================================================
##
##
##
## ------- Set options / paths. EDIT HERE:  ---------------------------------------------------------------------
##
{
      ps1_dir <-  file.path( Sys.getenv("HOME"), "Documents", "Work", "PhD_work", "Alg_paper_analysis",
                             "paper_1_chunking_and_parallel_scalability")
      ##
      ## ---- runners whose CURRENT settings decide which saved NicoStan burn-in runs are used:
      ##
      BayesMVP_runner_file <-  file.path(ps1_dir, "ps_1_burnin_optimizing_N_chunks_and_WCP_threads.R")
      Stan_runner_file     <-  file.path(ps1_dir, "ps_1_burnin_Stan_tape_chunked.R")
      sampling_runner_file <-  file.path(ps1_dir, "alg_paper_1_chunking_WCP_par_scaling.R")
      ##
      burnin_output_dir <-  file.path(ps1_dir, "burnin_outputs")
      N_vec             <-  c(500, 2500, 10000, 50000)
      devices           <-  c("HPC", "Laptop")
      ##
      ## ---- near-tie threshold (relative gap in mean seconds per burn-in iteration):
      ##
      near_tie_relative_gap <-  0.01
      ##
      ## ---- Stan via NicoStan rows of the build that loads the Stan model once per burn-in worker only:
      ##
      Stan_required_model_loading <-  "once_per_worker"
      ##
      ## ---- staging folder (the existing report outputs are not overwritten):
      ##
      staging_dir <-  file.path( ps1_dir, "report_outputs", "Burnin_chunks_WCP",
                                 "staging_selected_configurations_lowest_time")
}
##
##
##
## ------- Functions (defines only; no sampling, no compilation):  ----------------------------------------------
##
{
      `%>%` <-  magrittr::`%>%`   ## used inside the report helpers
      source(file.path(ps1_dir, "R_fns_alg_paper_1_burnin_report.R"))
      source(file.path(ps1_dir, "R_fns_alg_paper_1_burnin_report_Stan.R"))
}
##
##
##
## ------- Read the NicoStan burn-in runs (NicoStan+BayesMVP and the Stan model via NicoStan):  -----------------
##
{
      include_WCP_only <-  list(
            "BayesMVP" = isTRUE(x = fn_paper1_burnin_read_last_literal_assignment(
                                         runner_file_path = sampling_runner_file,
                                         object_name      = "paper1_settings$include_BayesMVP_WCP_only")),
            "Stan"     = isTRUE(x = fn_paper1_burnin_read_last_literal_assignment(
                                         runner_file_path = sampling_runner_file,
                                         object_name      = "paper1_settings$include_NicoStan_WCP_only")))
      BayesMVP_sampling_chunks <-  fn_paper1_burnin_read_last_literal_assignment(
                                         runner_file_path = sampling_runner_file,
                                         object_name      = "bayesmvp_chunk_candidates")
      file_regex <-  list(
            "BayesMVP" = "ps1_burnin_benchmark_N[0-9]+_L[0-9]+_n_runs[0-9]+\\.rds$",
            "Stan"     = "ps1_burnin_benchmark_Stan_tape_chunked_N[0-9]+_L[0-9]+_n_runs[0-9]+\\.rds$")
      ##
      NicoStan_rows <-  list("BayesMVP" = list(), "Stan" = list())
      for (implementation in c("BayesMVP", "Stan")) {
            for (device in devices) {

                  if (implementation == "BayesMVP") {
                        runner_settings <-  fn_paper1_burnin_read_runner_settings(
                                                    runner_file_path = BayesMVP_runner_file,
                                                    device           = device)
                        for (N_key in names(x = runner_settings$num_chunks_burnin_vec_given_N)) {
                              runner_settings$num_chunks_burnin_vec_given_N[[N_key]] <-  sort(x = unique(x = c(
                                    runner_settings$num_chunks_burnin_vec_given_N[[N_key]],
                                    BayesMVP_sampling_chunks[[N_key]])))
                        }
                  } else {
                        runner_settings <-  fn_paper1_burnin_read_Stan_runner_settings(
                                                    BayesMVP_runner_file_path = BayesMVP_runner_file,
                                                    sampling_runner_file_path = sampling_runner_file,
                                                    Stan_runner_file_path     = Stan_runner_file,
                                                    device                    = device)
                  }
                  expected_configurations <-  fn_paper1_burnin_expected_configurations(
                                                    runner_settings  = runner_settings,
                                                    N_values         = N_vec,
                                                    include_WCP_only = include_WCP_only[[implementation]])
                  file_prefix <-  fn_paper1_burnin_file_prefix(device_argument = device)
                  file_paths  <-  fn_paper1_burnin_list_files(
                                        output_dir       = burnin_output_dir,
                                        file_regex       = paste0("^", file_prefix, file_regex[[implementation]]),
                                        extra_file_paths = character(0))
                  device_rows <-  fn_paper1_burnin_read_files( file_paths     = file_paths,
                                                               device         = device,
                                                               implementation = implementation)
                  device_rows <-  device_rows[device_rows$N %in% N_vec, , drop = FALSE]
                  if (implementation == "Stan") {
                        device_rows <-  fn_paper1_burnin_keep_Stan_model_loading(
                                              rows                   = device_rows,
                                              required_model_loading = Stan_required_model_loading,
                                              device                 = device)
                  }
                  fn_paper1_burnin_stop_on_duplicate_runs(rows = device_rows, implementation = implementation)
                  device_rows <-  fn_paper1_burnin_filter_to_runner_settings(
                                        rows                    = device_rows,
                                        expected_configurations = expected_configurations,
                                        n_runs                  = runner_settings$n_runs,
                                        device                  = device,
                                        implementation          = implementation)$rows
                  ##
                  NicoStan_rows[[implementation]][[device]] <-  data.frame(
                        N                   = device_rows$N,
                        n_chains_burnin     = device_rows$n_chains_burnin,
                        n_threads_per_chain = device_rows$n_threads_WCP_burnin,
                        num_chunks          = device_rows$num_chunks_burnin,
                        run_number          = device_rows$run_number,
                        sec_per_iter        = device_rows$mean_sec_per_iter)

            }
      }
}
##
##
##
## ------- Read the Stan via cmdstanr burn-in runs (N_burn_chains = 4; per iteration of the slowest chain):  ---
##
{
      cmdstanr_rows <-  list()
      for (device in devices) {

            device_files <-  file.path( burnin_output_dir,
                                        paste0( device,
                                                "_ps1_burnin_Stan_cmdstanr_shootout_nchains4_treedepth4_N",
                                                N_vec, "_n_runs3.rds"))
            device_runs  <-  do.call(what = rbind, args = lapply(X = device_files, FUN = readRDS))
            cmdstanr_rows[[device]] <-  data.frame(
                  N                   = device_runs$N,
                  n_chains_burnin     = device_runs$n_chains_burnin,
                  n_threads_per_chain = device_runs$n_threads_per_chain,
                  num_chunks          = ifelse( test = is.na(x = device_runs$num_chunks),  ## AD_Stan: 1 chunk
                                                yes  = 1,
                                                no   = device_runs$num_chunks),
                  run_number          = device_runs$run_number,
                  sec_per_iter        = device_runs$time_sampling_slowest_chain_seconds / device_runs$n_iter)

      }
}
##
##
##
## ------- Select the fastest configuration per model, device and N, and write the table rows into staging:  -----
##
{
      rows_by_model <-  list( "NicoStan+BayesMVP"            = NicoStan_rows[["BayesMVP"]],
                              "Stan model (via NicoStan)"    = NicoStan_rows[["Stan"]],
                              "Stan model (via cmdstanr)"    = cmdstanr_rows)
      ##
      selections_by_model <-  list()
      all_selections      <-  list()
      for (model_title in names(x = rows_by_model)) {
            for (device in devices) {

                  selections <-  fn_paper1_burnin_select_lowest_time_over_all_burn_chains(
                                       rows                  = rows_by_model[[model_title]][[device]],
                                       near_tie_relative_gap = near_tie_relative_gap)
                  selections_by_model[[model_title]][[device]] <-  selections
                  all_selections[[length(x = all_selections) + 1]] <-  cbind( model  = model_title,
                                                                             device = device,
                                                                             selections)

            }
      }
      all_selections <-  do.call(what = rbind, args = all_selections)
      ##
      dir.create(path = staging_dir, recursive = TRUE, showWarnings = FALSE)
      utils::write.csv( x         = all_selections,
                        file      = file.path(staging_dir, "burnin_selected_configurations_lowest_time.csv"),
                        row.names = FALSE)
      fn_paper1_burnin_selected_configurations_table_rows_tex(
            selections_by_model = selections_by_model,
            N_values            = N_vec,
            file_path           = file.path(staging_dir, "table_burnin_selected_configurations_rows.tex"))
      ##
      print(x = all_selections[, c("model", "device", "N", "configuration_label", "configuration_type",
                                   "mean_sec_per_iter", "near_ties_within_1_percent",
                                   "near_ties_within_run_to_run_noise")])
}
























