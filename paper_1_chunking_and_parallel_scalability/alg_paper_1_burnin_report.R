##
## alg_paper_1_burnin_report.R
##
## Runner for the Paper 1 burn-in study report (chunks x WCP x burn-in chain count). Reads the saved PS1 burn-in .rds
## outputs from ps_1_burnin_optimizing_N_chunks_and_WCP_threads.R / R_fn_run_ps_1_burnin_benchmark.R (NOT edited here or by
## R_fns_alg_paper_1_burnin_report.R) and writes the Files/Generated/Burnin_chunks_WCP manuscript bundle.
##
## Every setting below is explicit - the helper functions in R_fns_alg_paper_1_burnin_report.R take no defaults and stop()
## if something required is missing, so edit this file (not the helpers) to change what gets read or where it is written.
##

rm(list = ls())
##
## ------- Set options / paths. EDIT HERE:  --------------------------------------------------------------------------------------------
##
{
      ps1_dir <- "/home/enzocerullo/Documents/Work/PhD_work/Alg_paper_analysis/paper_1_chunking_and_parallel_scalability"
      ##
      ## ---- devices to include in this report, and where each one's saved burn-in .rds files live:
      ##
      burnin_report_devices <- c("HPC", "Laptop")
      ##
      burnin_report_output_dir_by_device <- list(
            "HPC"    = file.path(ps1_dir, "burnin_outputs"),
            "Laptop" = file.path(ps1_dir, "burnin_outputs"))
      ##
      ## ---- where the bundle (paper_sections/Burnin_chunks_WCP/{figures,tables,section.tex}) is written:
      ##
      burnin_report_output_dir <- file.path(ps1_dir, "manuscript_outputs_burnin")
      ##
      ## ---- manuscript Files/Generated destination; NULL skips copying there (e.g. for a scratch/test run):
      ##
      burnin_report_manuscript_dir <- "/home/enzocerullo/Documents/Work/PhD_work/Alg_papers_LaTeX/paper_1_chunking_and_parallel_scalability"
      ##
      ## ---- grid cells to leave OUT of the report, keyed by N: N_WCP = 6 was run at N = 500 for every chain count but not for
      ##      8 chains at N = 2,500, so it is dropped at N = 500 to keep the N grids comparable. list() = none.
      ##
      burnin_report_excluded_WCP_given_N <- list("500" = c(6))
      ##
      ## ---- N values to report: the  burn-in rerun measured N = 500 / 2,500 / 10,000 on both devices; the HPC
      ##      N = 50,000 file on disk is from the build and is left out rather than mixed in.
      ##
      burnin_report_N_vec <- c(500, 2500, 10000)
}
##
## ------- Functions (defines only; no sampling, no compilation):  ---------------------------------------------------------------------
##
{
      source(file.path(ps1_dir, "R_fns_alg_paper_1_burnin_report.R"))
}
##
## ------- Read every requested device's saved burn-in outputs and bind them into one tidy table:  -------------------------------------
##
{
      burnin_report_rows_by_device <- vector(mode = "list", length = length(x = burnin_report_devices))
      names(x = burnin_report_rows_by_device) <- burnin_report_devices
      ##
      for (device_to_read in burnin_report_devices) {

            device_output_dir <- burnin_report_output_dir_by_device[[device_to_read]]
            if (is.null(x = device_output_dir)) {
                  stop("alg_paper_1_burnin_report.R: no burnin_report_output_dir_by_device entry for device = '", device_to_read, "'.")
            }
            burnin_report_rows_by_device[[device_to_read]] <- fn_paper1_burnin_read_outputs( output_dir = device_output_dir,
                                                                                              device     = device_to_read)

      }
      ##
      burnin_report_rows <- dplyr::bind_rows(burnin_report_rows_by_device)
      burnin_report_rows <- burnin_report_rows[burnin_report_rows$N %in% burnin_report_N_vec, , drop = FALSE]
      message(paste0("\033[36m", "alg_paper_1_burnin_report.R: kept N = ", paste(burnin_report_N_vec, collapse = ", "),
                     " (", nrow(burnin_report_rows), " rows)", "\033[0m"))
      ##
      burnin_report_rows <- fn_paper1_burnin_exclude_configurations( rows                 = burnin_report_rows,
                                                                     excluded_WCP_given_N = burnin_report_excluded_WCP_given_N)
      message(paste0("\033[36m", "alg_paper_1_burnin_report.R: read ", nrow(x = burnin_report_rows), " total row(s) across device(s): ",
                     paste(burnin_report_devices, collapse = ", "), "\033[0m"))
}
##
## ------- Write the manuscript bundle (figures, tables, section.tex; copied to manuscript_dir if set):  -------------------------------
##
{
      burnin_report_bundle <- fn_paper1_burnin_write_bundle( rows           = burnin_report_rows,
                                                              devices        = burnin_report_devices,
                                                              output_dir     = burnin_report_output_dir,
                                                              manuscript_dir = burnin_report_manuscript_dir)
      ##
      message(paste0("\033[36m", "alg_paper_1_burnin_report.R: bundle written to ", burnin_report_bundle$bundle_dir, "\033[0m"))
}






















