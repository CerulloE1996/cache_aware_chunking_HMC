##
## alg_paper_1_burnin_report_Stan.R
## (this copy: alg_paper_1_burnin_report_Stan_with_partial_sums_run_by_reduce_sum_static.R - the same report, written to the
##  folder given as its first argument, with the Stan N_chunks labels given as the number of partial sums which
##  reduce_sum_static() actually ran; with "requested" as its second argument, the requested N_chunks as before)
##
## Runner for the Paper 1 burn-in report with the tape-chunked Stan burn-in study added to the NicoStan+BayesMVP study.
## Reads the saved PS1 burn-in .rds outputs (never edits them), keeps only the runs that match the CURRENT settings of each burn-in
## runner (read from the runner files at run time), and writes:
##   - figure_burnin_HPC.png / figure_burnin_Laptop.png           (NicoStan+BayesMVP; existing figure function, unchanged),
##   - figure_burnin_Stan_tape_chunked.png                         (tape-chunked Stan, all devices, with the best NicoStan+BayesMVP level),
##   - burnin_best_by_cell_BayesMVP_and_Stan.csv                   (best configuration and speed-ups per device, N and burn-in chains),
##   - burnin_placeholder_values.csv / draft_section_filled.tex    (values of every \BurninTBD{ID}; the draft itself is not changed),
##   - burnin_text_checks.csv                                      (the existing burn-in sentences checked against the results),
##   - burnin_missing_configurations.csv                           (runner-grid configurations with fewer saved runs than n_runs).
## If burnin_report_manuscript_dir is set, the figures are also copied to <dir>/Files/Supplement/assets/Burnin_chunks_WCP/figures,
## the folder that Main.tex reads.
##
## Every setting is explicit below; the helper functions take no defaults. The environment variable BURNIN_REPORT_IMPLEMENTATIONS
## (e.g. "BayesMVP") restricts the implementations reported.
##

rm(list = ls())
##
options(paper1.Stan_N_chunks_as_partial_sums = !identical(x = commandArgs(trailingOnly = TRUE)[2], y = "requested"))
##
## ------- Set options / paths. EDIT HERE:  --------------------------------------------------------------------------------------------
##
{
      ps1_dir <- "/home/enzocerullo/Documents/Work/PhD_work/Alg_paper_analysis/paper_1_chunking_and_parallel_scalability"
      ##
      ## ---- folder holding this runner, R_fns_alg_paper_1_burnin_report_Stan.R and draft_section.tex:
      ##
      burnin_report_Stan_dir <- ps1_dir   ## the Stan burn-in report functions now sit next to this script
      ##
      ## ---- runners whose CURRENT settings decide which saved runs are reported:
      ##
      burnin_report_BayesMVP_runner_file <- file.path(ps1_dir, "ps_1_burnin_optimizing_N_chunks_and_WCP_threads.R")
      burnin_report_Stan_runner_file     <- file.path(ps1_dir, "ps_1_burnin_Stan_tape_chunked.R")
      burnin_report_sampling_runner_file <- file.path(ps1_dir, "alg_paper_1_chunking_WCP_par_scaling.R")
      ##
      ## ---- implementations and devices reported:
      ##
      burnin_report_devices_by_implementation <- list(
            "BayesMVP" = c("HPC", "Laptop"),
            "Stan"     = c("HPC", "Laptop"))
      burnin_report_implementations <- names(x = burnin_report_devices_by_implementation)
      if (nzchar(Sys.getenv("BURNIN_REPORT_IMPLEMENTATIONS"))) {
            burnin_report_implementations <- strsplit(x = Sys.getenv("BURNIN_REPORT_IMPLEMENTATIONS"), split = ",", fixed = TRUE)[[1]]
      }
      if (nzchar(Sys.getenv("BURNIN_REPORT_STAN_DEVICES"))) {
            burnin_report_devices_by_implementation[["Stan"]] <- strsplit(x = Sys.getenv("BURNIN_REPORT_STAN_DEVICES"), split = ",", fixed = TRUE)[[1]]
      }
      ##
      ## ---- Stan rows of the NicoStan build that loads the Stan model once per burn-in worker only (column Stan_model_loading):
      ##
      burnin_report_Stan_required_model_loading <- "once_per_worker"
      ##
      ## ---- leapfrog steps of every timed burn-in iteration (for per-leapfrog-step comparisons of Stan and NicoStan+BayesMVP):
      ##
      burnin_report_leapfrog_steps_csv <- file.path(ps1_dir, "burnin_outputs", "leapfrog_steps_per_timed_iteration_BayesMVP.csv")
      ##
      ## ---- where the saved .rds files live, and the file-name pattern of each implementation (the "HPC_" / "Laptop_" prefix is added):
      ##
      burnin_report_output_dir <- file.path(ps1_dir, "burnin_outputs")
      ##
      burnin_report_file_regex_by_implementation <- list(
            "BayesMVP" = "ps1_burnin_benchmark_N[0-9]+_L[0-9]+_n_runs[0-9]+\\.rds$",
            "Stan"     = "ps1_burnin_benchmark_Stan_tape_chunked_N[0-9]+_L[0-9]+_n_runs[0-9]+\\.rds$")
      ##
      ## ---- additional saved files read on top of the pattern matches. character(0) = none:
      ##
      burnin_report_extra_files <- list(
            "BayesMVP" = list("HPC" = character(0), "Laptop" = character(0)),
            "Stan"     = list("HPC" = character(0), "Laptop" = character(0)))
      ##
      ## ---- Stan rows completed before this time are left out (e.g. runs of an earlier build), "YYYY-mm-dd HH:MM:SS"; NULL = keep all:
      ##
      burnin_report_Stan_minimum_time_completed <- list("HPC" = NULL, "Laptop" = NULL)
      ##
      ## ---- N values to report, and grid cells left out of the report (none; alg_paper_1_burnin_report.R still leaves out N_WCP = 6 at N = 500):
      ##
      burnin_report_N_vec                <- c(500, 2500, 10000, 50000)
      ## burnin_report_excluded_WCP_given_N <- list("500" = c(6))
      ## ---- none: N_WCP = 6 at N = 500 (16 chains) is in the current runner grid and was run 3 times, so it is reported:
      burnin_report_excluded_WCP_given_N <- list()
      ##
      ## ---- outputs: figures / values folder, manuscript copy (the v32 package in the scratch folder; NULL = no copy), draft section:
      ##
      burnin_report_results_dir    <- file.path(burnin_report_Stan_dir, "report_outputs")
      burnin_report_results_dir    <- commandArgs(trailingOnly = TRUE)[1]   ## this copy: its own output folder
      burnin_report_manuscript_dir <- file.path(burnin_report_Stan_dir, "v32_package")
      burnin_report_manuscript_dir <- NULL   ## this copy: no manuscript copy
      burnin_report_draft_tex_path <- file.path(burnin_report_Stan_dir, "draft_section.tex")
}
##
## ------- Functions (defines only; no sampling, no compilation):  ---------------------------------------------------------------------
##
{
      source(file.path(ps1_dir, "R_fns_alg_paper_1_burnin_report.R"))
      source(file.path(burnin_report_Stan_dir, "R_fns_alg_paper_1_burnin_report_Stan.R"))
      source(file.path(ps1_dir, "R_fn_number_of_partial_sums_run_by_reduce_sum_static.R"))
}
##
## ------- Read, check and filter every implementation's saved burn-in outputs:  -------------------------------------------------------
##
{
      burnin_report_leapfrog_steps <- fn_paper1_burnin_read_leapfrog_steps(csv_path = burnin_report_leapfrog_steps_csv)
      ##
      ## ---- WCP-only configurations (N_chunks = N_WCP) and the NicoStan+BayesMVP chunk counts, as in the sampling runner:
      ##
      burnin_report_include_WCP_only <- list(
            "BayesMVP" = isTRUE(x = fn_paper1_burnin_read_last_literal_assignment( runner_file_path = burnin_report_sampling_runner_file,
                                                                                    object_name      = "paper1_settings$include_BayesMVP_WCP_only")),
            "Stan"     = isTRUE(x = fn_paper1_burnin_read_last_literal_assignment( runner_file_path = burnin_report_sampling_runner_file,
                                                                                    object_name      = "paper1_settings$include_NicoStan_WCP_only")))
      burnin_report_BayesMVP_sampling_chunks <- fn_paper1_burnin_read_last_literal_assignment( runner_file_path = burnin_report_sampling_runner_file,
                                                                                              object_name      = "bayesmvp_chunk_candidates")
      fn_paper1_burnin_message(text = paste0("alg_paper_1_burnin_report_Stan.R: WCP-only configurations included: BayesMVP = ",
                                             burnin_report_include_WCP_only[["BayesMVP"]], ", Stan = ", burnin_report_include_WCP_only[["Stan"]]))
      ##
      burnin_report_rows                    <- list("BayesMVP" = list(), "Stan" = list())
      burnin_report_expected_configurations <- list()
      burnin_report_missing_configurations  <- list()
      ##
      for (implementation in burnin_report_implementations) {
            for (device_to_read in burnin_report_devices_by_implementation[[implementation]]) {

                  if (implementation == "BayesMVP") {
                        runner_settings <- fn_paper1_burnin_read_runner_settings( runner_file_path = burnin_report_BayesMVP_runner_file,
                                                                                  device           = device_to_read)
                        ##
                        ## ---- chunk counts: the burn-in runner's list together with any count of the sampling runner's BayesMVP grid it lacks:
                        ##
                        for (N_key in names(x = runner_settings$num_chunks_burnin_vec_given_N)) {
                              runner_settings$num_chunks_burnin_vec_given_N[[N_key]] <- sort(x = unique(x = c(runner_settings$num_chunks_burnin_vec_given_N[[N_key]],
                                                                                                              burnin_report_BayesMVP_sampling_chunks[[N_key]])))
                        }
                  } else {
                        runner_settings <- fn_paper1_burnin_read_Stan_runner_settings( BayesMVP_runner_file_path = burnin_report_BayesMVP_runner_file,
                                                                                       sampling_runner_file_path = burnin_report_sampling_runner_file,
                                                                                       Stan_runner_file_path     = burnin_report_Stan_runner_file,
                                                                                       device                    = device_to_read)
                  }
                  expected_configurations <- fn_paper1_burnin_expected_configurations( runner_settings  = runner_settings,
                                                                                        N_values         = burnin_report_N_vec,
                                                                                        include_WCP_only = burnin_report_include_WCP_only[[implementation]])
                  burnin_report_expected_configurations[[paste0(implementation, "_", device_to_read)]] <- expected_configurations
                  ##
                  file_paths <- fn_paper1_burnin_list_files( output_dir       = burnin_report_output_dir,
                                                             file_regex       = paste0("^", fn_paper1_burnin_file_prefix(device_argument = device_to_read),
                                                                                       burnin_report_file_regex_by_implementation[[implementation]]),
                                                             extra_file_paths = burnin_report_extra_files[[implementation]][[device_to_read]])
                  if (length(x = file_paths) == 0) next
                  ##
                  device_rows <- fn_paper1_burnin_read_files( file_paths     = file_paths,
                                                              device         = device_to_read,
                                                              implementation = implementation)
                  device_rows <- device_rows[device_rows$N %in% burnin_report_N_vec, , drop = FALSE]
                  if (implementation == "Stan") {
                        device_rows <- fn_paper1_burnin_keep_Stan_model_loading( rows                   = device_rows,
                                                                                 required_model_loading = burnin_report_Stan_required_model_loading,
                                                                                 device                 = device_to_read)
                        device_rows <- fn_paper1_burnin_filter_time_completed( rows                   = device_rows,
                                                                               minimum_time_completed = burnin_report_Stan_minimum_time_completed[[device_to_read]],
                                                                               implementation         = implementation,
                                                                               device                 = device_to_read)
                  }
                  fn_paper1_burnin_stop_on_duplicate_runs(rows = device_rows, implementation = implementation)
                  ##
                  filtered_outputs <- fn_paper1_burnin_filter_to_runner_settings( rows                    = device_rows,
                                                                                  expected_configurations = expected_configurations,
                                                                                  n_runs                  = runner_settings$n_runs,
                                                                                  device                  = device_to_read,
                                                                                  implementation          = implementation)
                  burnin_report_missing_configurations[[paste0(implementation, "_", device_to_read)]] <- filtered_outputs$missing_configurations
                  ##
                  device_rows <- fn_paper1_burnin_exclude_configurations( rows                 = filtered_outputs$rows,
                                                                          excluded_WCP_given_N = burnin_report_excluded_WCP_given_N)
                  burnin_report_rows[[implementation]][[device_to_read]] <- fn_paper1_burnin_add_sec_per_step( rows           = device_rows,
                                                                                                               leapfrog_steps = burnin_report_leapfrog_steps)

            }
      }
      ##
      burnin_report_missing_configurations <- dplyr::bind_rows(burnin_report_missing_configurations)
}
##
## ------- Write the figures, values, checks and the filled draft:  --------------------------------------------------------------------
##
{
      burnin_report_Stan_bundle <- fn_paper1_burnin_write_Stan_bundle( BayesMVP_rows_by_device                              = burnin_report_rows[["BayesMVP"]],
                                                                       Stan_rows_by_device                                  = burnin_report_rows[["Stan"]],
                                                                       expected_configurations_by_implementation_and_device = burnin_report_expected_configurations,
                                                                       missing_configurations                               = burnin_report_missing_configurations,
                                                                       output_dir                                           = burnin_report_results_dir,
                                                                       manuscript_dir                                       = burnin_report_manuscript_dir,
                                                                       draft_tex_path                                       = burnin_report_draft_tex_path)
      ##
      fn_paper1_burnin_message(text = paste0("alg_paper_1_burnin_report_Stan.R: outputs written to ", burnin_report_results_dir))
}























