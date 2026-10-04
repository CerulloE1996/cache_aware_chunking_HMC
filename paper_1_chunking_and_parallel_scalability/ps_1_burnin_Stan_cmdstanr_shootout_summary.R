##
## ===============================================================================================================
## ps_1_burnin_Stan_cmdstanr_shootout_summary.R
##
## Summarises the Stan via cmdstanr burn-in shootout (ps_1_burnin_Stan_cmdstanr_shootout_4_chains.R)
## for E1 Part VII:
##   - the fastest configuration per device and N (seconds per leapfrog step of the slowest chain),
##     its speed-up relative to AD_Stan with N_chains x N_threads/chain = 4 x 1,
##     and every configuration within 1% of it;
##   - the best configuration of each implementation (AD_Stan, tape chunking, WCP-only, WCP + chunking);
##   - the ratio of its seconds per leapfrog step to that of the best Stan via NicoStan configuration
##     with burn-in N_chains = 4
##     (report_outputs/Burnin_chunks_WCP/values/burnin_best_by_cell_BayesMVP_and_Stan.csv);
##   - the two rows of table table:paper1_burnin_selected_configurations, given as N_threads/chain / N_chunks.
## ===============================================================================================================
##
##
##
## ---- Settings: ------------------------------------------------------------------------------------------------
##
paper_1_dir <-  file.path( Sys.getenv("HOME"), "Documents", "Work", "PhD_work", "Alg_paper_analysis",
                           "paper_1_chunking_and_parallel_scalability")
##
N_vec <-  c(500, 2500, 10000, 50000)
##
NicoStan_values_file <-  file.path( paper_1_dir, "report_outputs", "Burnin_chunks_WCP", "values",
                                    "burnin_best_by_cell_BayesMVP_and_Stan.csv")
##
fn_shootout_file <-  function(device, N) {
        file.path( paper_1_dir, "burnin_outputs",
                   paste0( device, "_ps1_burnin_Stan_cmdstanr_shootout_nchains4_treedepth4_N", N,
                           "_n_runs3.rds"))
}
##
NicoStan_values <-  utils::read.csv(NicoStan_values_file)
NicoStan_values <-  NicoStan_values[ NicoStan_values$implementation == "Stan_tape_chunked" &
                                     NicoStan_values$n_chains_burnin == 4, ]
##
##
##
## ---- Per device and N: ----------------------------------------------------------------------------------------
##
summary_rows <-  list()
##
for (device in c("HPC", "Laptop")) {
      for (N in N_vec) {

            file_path <-  fn_shootout_file(device = device, N = N)
            if (!file.exists(file_path)) { message(paste0(device, ", N = ", N, ": not run yet")) ; next }
            ##
            runs <-  readRDS(file_path)
            runs$num_chunks[is.na(runs$num_chunks)] <-  0
            ##
            configurations <-  stats::aggregate( sec_per_step_slowest_chain ~
                                                     algorithm + n_threads_per_chain + num_chunks,
                                                 data = runs, FUN = mean)
            configurations <-  configurations[order(configurations$sec_per_step_slowest_chain), ]
            ##
            sec_per_step_AD_Stan <-  configurations$sec_per_step_slowest_chain[
                                          configurations$algorithm == "AD_Stan"]
            best <-  configurations[1, ]
            ##
            within_1_percent <-  configurations[ configurations$sec_per_step_slowest_chain <=
                                                     1.01 * best$sec_per_step_slowest_chain, ]
            ##
            NicoStan_best <-  NicoStan_values[NicoStan_values$device == device & NicoStan_values$N == N, ]
            ##
            best_of_each <-  do.call(rbind, lapply( split(configurations, configurations$algorithm),
                                                    function(rows) rows[1, ]))
            ##
            message(paste0( "\n", device, ", N = ", N, ": fastest ", best$algorithm,
                            " 4 x ", best$n_threads_per_chain,
                            ", N_chunks ", best$num_chunks, " (",
                            formatC( sec_per_step_AD_Stan / best$sec_per_step_slowest_chain,
                                     format = "f", digits = 2),
                            "x vs AD_Stan 4 x 1); NicoStan best ",
                            NicoStan_best$best_n_threads_WCP_burnin, "/",
                            NicoStan_best$best_num_chunks_burnin, ", NicoStan sec/step / cmdstanr sec/step = ",
                            formatC(NicoStan_best$best_sec_per_step / best$sec_per_step_slowest_chain,
                                    format = "f", digits = 2)))
            message(paste0( "  within 1%: ",
                            paste0( within_1_percent$algorithm, " 4 x ", within_1_percent$n_threads_per_chain,
                                    "/", within_1_percent$num_chunks, collapse = "; ")))
            for (row_index in seq_len(nrow(best_of_each))) {
                  row <-  best_of_each[row_index, ]
                  message(paste0( "  best ", row$algorithm, ": 4 x ", row$n_threads_per_chain, "/",
                                  row$num_chunks, ", ",
                                  formatC( sec_per_step_AD_Stan / row$sec_per_step_slowest_chain,
                                           format = "f", digits = 2),
                                  "x"))
            }
            ##
            summary_rows[[length(summary_rows) + 1]] <-  data.frame(
                  device = device,
                  N = N,
                  best_algorithm = best$algorithm,
                  best_n_threads_per_chain = best$n_threads_per_chain,
                  best_num_chunks = best$num_chunks,
                  best_sec_per_step = best$sec_per_step_slowest_chain,
                  speed_up_vs_AD_Stan_4x1 = sec_per_step_AD_Stan / best$sec_per_step_slowest_chain,
                  n_within_1_percent = nrow(within_1_percent),
                  NicoStan_best_n_threads_per_chain = NicoStan_best$best_n_threads_WCP_burnin,
                  NicoStan_best_num_chunks = NicoStan_best$best_num_chunks_burnin,
                  NicoStan_over_cmdstanr_sec_per_step = NicoStan_best$best_sec_per_step /
                                                            best$sec_per_step_slowest_chain)

      }
}
##
summary_table <-  do.call(rbind, summary_rows)
##
utils::write.csv( x = summary_table,
                  file = file.path(paper_1_dir, "report_outputs", "Burnin_chunks_WCP", "values",
                                   "burnin_Stan_cmdstanr_shootout_summary.csv"),
                  row.names = FALSE)
##
##
##
## ---- Rows of table table:paper1_burnin_selected_configurations (dagger = WCP-only): --------------------------
##
for (device in c("HPC", "Laptop")) {
      rows <-  summary_table[summary_table$device == device, ]
      entries <-  paste0( rows$best_n_threads_per_chain, "/", rows$best_num_chunks,
                          ifelse(rows$best_algorithm == "AD_Stan_WCP", "$^{\\dagger}$", ""))
      row_start <-  if (device == "HPC") "Stan (cmdstanr)   & local-HPC & 4  & " else
                                         "                  & laptop    & 4  & "
      message(paste0( row_start,
                      paste0(entries, collapse = " & "), " \\\\"))
}
























