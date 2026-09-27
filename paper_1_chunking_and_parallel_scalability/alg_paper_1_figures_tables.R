##
## =====================================================================================================================================
## alg_paper_1_figures_tables.R
##
## Paper 1: regenerate figures and tables from saved results without sampling.
## Defaults to the same stable device folder as the main runner. Supply both device paths to combine studies.
## With character(0), RStudio opens a chooser: select results.rds inside the study directory.
## To combine HPC and Laptop, put both explicit directory paths in this vector.
## Reports first show chunk optimisation within both WCP arms at fixed chains/WCP,
## then compare the filtered configurations with the other algorithms.
##
## ---- Select completed studies and output locations ------------------------------------------------------------------------------------
##
study_output_dirs <-  file.path(path.expand('~/Documents/Work/PhD_work/Alg_paper_analysis'),
                                'paper_1_chunking_and_parallel_scalability', 'paper_1_computational_outputs',
                                if (parallel::detectCores() > 16) 'HPC' else 'Laptop')
##
## Example: study_output_dirs <- c('/path/to/paper_1_computational_outputs/HPC', '/path/to/paper_1_computational_outputs/Laptop')
##
report_output_dir <-  NULL  ## NULL reuses manuscript_outputs in the selected study directory.
##
algorithm_study_dir <-  path.expand('~/Documents/Work/PhD_work/Alg_paper_analysis')
##
helper_dir <-  file.path(algorithm_study_dir, 'paper_1_chunking_and_parallel_scalability')
##
## ---- Load the reporting functions and select saved results -----------------------------------------------------------------------------
##
source(file.path(helper_dir, 'R_fns_alg_paper_1_figures_tables.R'), local = TRUE)
##
if (!length(study_output_dirs)) {

    if (!interactive()) stop('Set study_output_dirs to the completed study directory/directories.')
    ##
    message('Select results.rds from the completed Paper 1 study.')
    ##
    selected_results_file <-  file.choose()
    ##
    if (basename(selected_results_file) != 'results.rds') stop('Select the unified study results.rds file.')
    ##
    study_output_dirs <-  dirname(selected_results_file)

}
##
## ---- Export figures and tables from the selected studies -------------------------------------------------------------------------------
##
##
## ---- Grid cells to leave OUT of the report (design consistency; the saved results are untouched, dropped rows are counted in the
##      READ_ME). The HPC WCP grid ran threads_per_chain = 6 with 8 chains at N = 500 but not at N = 2,500:
##
report_excluded_cases <-  list( list(device = "HPC", N = 500, n_chains = 8, threads_per_chain = 6))
##
paper1_report <-  fn_paper1_export_manuscript( study_output_dirs,
                                               output_dir = report_output_dir,
                                               helper_dir = helper_dir,
                                               manuscript_dir = file.path(dirname(dirname(helper_dir)), "Alg_papers_LaTeX",
                                                                          "paper_1_chunking_and_parallel_scalability"),
                                               excluded_cases = report_excluded_cases)
##
message(paste0('\033[36mReport saved in: ', paper1_report$output_dir, '\033[0m'))






















