##
## ======================================================================================================================================
## alg_paper_1_Stan_chunk_search_figures_with_partial_sums_run_by_reduce_sum_static.R
##
## Paper 1: the four Stan chunk-search figures (AD_Stan_WCP vs. AD_Stan_WCP_chunking; local-HPC and laptop, N = 500 and 50,000), with
## N_chunks given as the number of partial sums which reduce_sum_static() actually ran
## (R_fn_number_of_partial_sums_run_by_reduce_sum_static.R). Reads the saved presentation views only (no sampling).
##
## With show_requested_N_chunks = TRUE, the figures are drawn with the requested N_chunks instead (as in the paper before the
## relabelling), to check that this script reproduces the earlier figures.
##
{
      paper1_analysis_dir <-  path.expand("~/Documents/Work/PhD_work/Alg_paper_analysis/paper_1_chunking_and_parallel_scalability")
      ##
      presentation_views_file <-  file.path(paper1_analysis_dir, "paper_1_computational_outputs",
                                             "manuscript_outputs_final_both_devices_narrow_WCP_2026_10_03", "data", "presentation_views.rds")
      ##
      args <-  commandArgs(trailingOnly = TRUE)
      output_dir <-  if (length(x = args) >= 1) args[1] else file.path(paper1_analysis_dir, "report_outputs", "Stan_partial_sums_figures")
      show_requested_N_chunks <-  length(x = args) >= 2 && identical(x = args[2], y = "requested")
      ##
      dir.create(path = output_dir, recursive = TRUE, showWarnings = FALSE)
}
##
## ---- Report functions, presentation templates and the partial-sum count ---------------------------------------------------------------
##
{
      source(file.path(paper1_analysis_dir, "R_fns_alg_paper_1_figures_tables.R"), local = TRUE)
      fn_paper1_report_dependencies()
      ##
      private <-  new.env(parent = environment())
      sys.source(file.path(paper1_analysis_dir, "R_fns_alg_paper_1_presentation_templates.R"), envir = private)
      templates <-  private$fn_paper1_presentation_templates()
      ##
      source(file.path(paper1_analysis_dir, "R_fn_number_of_partial_sums_run_by_reduce_sum_static.R"), local = TRUE)
      ##
      views <-  readRDS(file = presentation_views_file)
}
##
## ---- One figure for each device and N --------------------------------------------------------------------------------------------------
##
for (device in c("HPC", "Laptop")) {

      for (N in c(500, 50000)) {

            chunk_search <-  views$wcp_chunk_search[views$wcp_chunk_search$device == device &
                                                    views$wcp_chunk_search$algorithm == "AD_Stan_WCP" &
                                                    views$wcp_chunk_search$N == N, , drop = FALSE]
            chunk_search <-  chunk_search[order(chunk_search$n_chains, chunk_search$threads_per_chain, chunk_search$num_chunks), ,
                                          drop = FALSE]
            ##
            if (!show_requested_N_chunks) {

                  chunk_search$num_chunks_requested <-  chunk_search$num_chunks
                  chunk_search$num_chunks <-  fn_number_of_partial_sums_run_by_reduce_sum_static( N_units            = chunk_search$N,
                                                                                                  N_chunks_requested = chunk_search$num_chunks)

            }
            ##
            best_chunks <-  chunk_search[chunk_search$selected_best_chunks, , drop = FALSE]
            ##
            file_prefix <-  paste0(device, "_paper_1_N_", N, "_algorithm_AD_Stan_WCP")
            ##
            templates$fn_plot_paper1_WCP_chunk_search( chunk_search    = chunk_search,
                                                       best_chunks     = best_chunks,
                                                       output_path     = output_dir,
                                                       file_prefix     = file_prefix,
                                                       show_cache_lines = TRUE,
                                                       bytes_per_row   = templates$paper1_bytes_per_row_Stan,
                                                       thresholds_file = file.path(output_dir, paste0(file_prefix,
                                                                                                      "_chunk_search_cache_lines_thresholds.csv")))
            ##
            message(paste0("\033[36m", "Saved: ", file.path(output_dir, paste0(file_prefix, "_chunk_search_cache_lines.png")), "\033[0m"))

      }

}
























