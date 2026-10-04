##
## ===============================================================================================================
## ps_1_burnin_Stan_cmdstanr_shootout_4_chains.R
##
## Burn-in configuration shootout for Stan (cmdstanr) at a FIXED burn-in chain count, on the COVID-19 data.
## This is an updated re-run of the PS2 supplement "BURNIN-CONFIG SHOOTOUT AT FIXED n_chains = 4"
## (legacy/ps_2_parallel_scaling_vs_Mplus_Stan/ps_2_parallel_scaling_vs_Mplus_Stan_v2.R),
## whose results were never saved; it informs the Stan burn-in configuration of E3 / Paper 2.
##
## Changes from the legacy shootout:
##   - the current Paper 1 Stan models: AD_Stan, and the reduce_sum model for tape chunking (N_threads/chain = 1),
##     WCP-only (N_chunks = N_threads/chain) and WCP + chunking (N_chunks > N_threads/chain);
##   - max_treedepth = 4, i.e., 15 leapfrog steps per iteration with the fixed step size, the closest NUTS match
##     to the mean of 20 leapfrog steps of the Stan-via-NicoStan burn-in study (E1 Part VII),
##     whose timed iterations (200, 40, 5 and 3 at N = 500, 2500, 10,000 and 50,000) are also used here;
##   - every candidate N_chunks of the Stan grid (table ps1_pilot_study_params_Stan)
##     with N_chunks >= N_threads/chain, plus WCP-only,
##     at N_threads/chain in {1, 2, 4, 8, 16, 22, 44} (local-HPC) or {1, 2, 4} (laptop);
##   - the main endpoint is seconds per leapfrog step of the slowest chain (cmdstan's per-chain sampling time),
##     since the burn-in chains advance together; the wall-clock time of each $sample() call is also kept.
##
## One .rds file per device and N, named by its settings; an existing file is reused, never re-run.
## ===============================================================================================================
##
##
##
## ---- Settings: ------------------------------------------------------------------------------------------------
##
algorithm_study_dir <-  file.path(Sys.getenv("HOME"), "Documents", "Work", "PhD_work", "Alg_paper_analysis")
##
device <-  if (parallel::detectCores() < 17) "Laptop" else "HPC"
##
N_vec <-  c(500, 2500, 10000, 50000)
##
n_chains_burnin <-  4
n_runs          <-  3
step_size       <-  0.00001
max_treedepth   <-  4
init_value      <-  0.01
##
n_iter_given_N <-  list("500" = 200, "2500" = 40, "10000" = 5, "50000" = 3)
##
## ---- Stan N_chunks grid (table ps1_pilot_study_params_Stan):
##
chunk_grid_given_N <-  list( "500"   = c(2, 4, 10),
                             "2500"  = c(4, 10, 25, 50, 100, 250),
                             "10000" = c(4, 10, 25, 50, 100, 250, 500, 1000),
                             "50000" = c(10, 25, 50, 100, 250, 500, 1000, 2000, 5000))
##
threads_per_chain_vec <-  if (device == "HPC") c(1, 2, 4, 8, 16, 22, 44) else c(1, 2, 4)
##
output_dir <-  file.path(algorithm_study_dir, "paper_1_chunking_and_parallel_scalability", "burnin_outputs")
##
fn_file_name <-  function(N) {
        paste0( device, "_ps1_burnin_Stan_cmdstanr_shootout_nchains", n_chains_burnin,
                "_treedepth", max_treedepth, "_N", N, "_n_runs", n_runs, ".rds")
}
##
##
##
## ---- Data (COVID-19 simulator, seed 123, as in the legacy shootout) and the compiled Stan models: -------------
##
source(file.path( algorithm_study_dir,
                  "paper_1_chunking_and_parallel_scalability",
                  "R_fns_alg_paper_1_chunking_WCP_par_scaling.R"))
##
true_vals_list <-  fn_paper1_COVID_data( algorithm_study_dir = algorithm_study_dir,
                                        N_vec = N_vec,
                                        seed = 123)
##
Stan_model_obj <-  fn_paper1_compile_stan( settings = list( algorithm_study_dir = algorithm_study_dir,
                                                            force_recompile = FALSE,
                                                            device = device),
                                           algorithms = c("AD_Stan", "AD_Stan_tape_chunked"))
##
##
##
## ---- Configuration grid for one N: ----------------------------------------------------------------------------
##
fn_configuration_grid <-  function(N) {
        configurations <-  list(list( algorithm = "AD_Stan", model = "AD_Stan",
                                      threads_per_chain = 1, num_chunks = NA))
        ##
        for (threads_per_chain in threads_per_chain_vec) {
              num_chunks_vec <-  chunk_grid_given_N[[as.character(N)]]
              num_chunks_vec <-  num_chunks_vec[num_chunks_vec >= threads_per_chain]
              if (threads_per_chain > 1) num_chunks_vec <-  sort(unique(c(threads_per_chain, num_chunks_vec)))
              ##
              for (num_chunks in num_chunks_vec) {
                    algorithm <-  if (threads_per_chain == 1) {
                        "AD_Stan_tape_chunked"
                    } else if (num_chunks == threads_per_chain) {
                        "AD_Stan_WCP"
                    } else "AD_Stan_WCP_chunking"
                    ##
                    configurations[[length(configurations) + 1]] <-  list( algorithm = algorithm,
                                                                          model = "AD_Stan_tape_chunked",
                                                                          threads_per_chain = threads_per_chain,
                                                                          num_chunks = num_chunks)
              }
        }
        ##
        return(configurations)
}
##
##
##
## ---- Run every N (smallest first), saving one file per N: -----------------------------------------------------
##
for (N in N_vec) {

      file_path <-  file.path(output_dir, fn_file_name(N))
      ##
      if (file.exists(file_path)) {
            message(paste0("N = ", N, ": reusing saved runs in ", basename(file_path)))
            next
      }
      ##
      n_iter <-  n_iter_given_N[[as.character(N)]]
      ##
      y <-  true_vals_list$y_list[[which(N_vec == N)]]
      ##
      Stan_data_base <-  fn_paper1_stan_data(y = y)
      ##
      configurations <-  fn_configuration_grid(N = N)
      ##
      message(paste0("N = ", N, ": ", length(configurations), " configurations x ", n_runs, " runs, ", n_iter,
                     " iterations, ", device, ", ", format(Sys.time(), "%H:%M:%S")))
      ##
      results_list <-  list()
      ##
      for (configuration in configurations) {

            Stan_data <-  Stan_data_base
            ##
            chunk_size <-  if (is.na(configuration$num_chunks)) NA else
                               max(1, ceiling(N / configuration$num_chunks))
            ##
            if (!is.na(chunk_size)) Stan_data$chunk_size <-  as.integer(chunk_size)
            ##
            n_threads_total <-  n_chains_burnin * configuration$threads_per_chain
            ##
            for (run_number in seq_len(n_runs)) {

                  time_start <-  proc.time()[["elapsed"]]
                  ##
                  fit <-  Stan_model_obj[[configuration$model]]$sample(
                                data              = Stan_data,
                                seed              = run_number * 1000 + n_threads_total,
                                chains            = n_chains_burnin,
                                parallel_chains   = n_chains_burnin,
                                threads_per_chain = if (configuration$model == "AD_Stan") NULL else
                                                        configuration$threads_per_chain,
                                iter_warmup       = 0,
                                iter_sampling     = n_iter,
                                adapt_engaged     = FALSE,
                                step_size         = step_size,
                                metric            = "diag_e",
                                init              = init_value,
                                refresh           = 0,
                                show_messages     = FALSE,
                                save_warmup       = FALSE,
                                max_treedepth     = max_treedepth)
                  ##
                  time_wall_seconds <-  proc.time()[["elapsed"]] - time_start
                  ##
                  ## ---- after the clock: per-chain sampling times and leapfrog steps
                  ##
                  time_sampling_each_chain <-  fit$time()$chains$sampling
                  diagnostics <-  fit$sampler_diagnostics(format = "draws_df")
                  mean_n_leapfrog <-  mean(diagnostics$n_leapfrog__)
                  ##
                  results_list[[length(results_list) + 1]] <-  data.frame(
                        device = device,
                        N = N,
                        algorithm = configuration$algorithm,
                        n_chains_burnin = n_chains_burnin,
                        n_threads_per_chain = configuration$threads_per_chain,
                        n_threads_total = n_threads_total,
                        num_chunks = configuration$num_chunks,
                        chunk_size = chunk_size,
                        run_number = run_number,
                        n_iter = n_iter,
                        max_treedepth = max_treedepth,
                        mean_n_leapfrog = mean_n_leapfrog,
                        time_wall_seconds = time_wall_seconds,
                        time_sampling_slowest_chain_seconds = max(time_sampling_each_chain),
                        time_sampling_mean_chain_seconds = mean(time_sampling_each_chain),
                        sec_per_step_slowest_chain = max(time_sampling_each_chain) / (n_iter * mean_n_leapfrog),
                        n_divergences = sum(diagnostics$divergent__))
                  ##
                  rm(fit) ; invisible(gc())

            }
            ##
            last_rows <-  utils::tail(results_list, n_runs)
            message(paste0( "  ", configuration$algorithm, " ",
                            n_chains_burnin, "x", configuration$threads_per_chain,
                            ", N_chunks ", configuration$num_chunks, ": sec/step (slowest chain) ",
                            formatC( mean(vapply(last_rows, function(r) r$sec_per_step_slowest_chain,
                                                 numeric(1))),
                                     format = "e", digits = 3),
                            ", wall ",
                            formatC( mean(vapply(last_rows, function(r) r$time_wall_seconds, numeric(1))),
                                     format = "f", digits = 2), " s"))

      }
      ##
      results <-  do.call(rbind, results_list)
      saveRDS(object = results, file = file_path)
      ##
      message(paste0("N = ", N, ": saved ", basename(file_path), ", ", format(Sys.time(), "%H:%M:%S")))

}
##
##
##
## ---- Summary: the fastest configuration per N, and the best of each implementation: ---------------------------
##
for (N in N_vec) {

      file_path <-  file.path(output_dir, fn_file_name(N))
      if (!file.exists(file_path)) next
      ##
      results <-  readRDS(file_path)
      ##
      summary_table <-  stats::aggregate( cbind(sec_per_step_slowest_chain, time_wall_seconds) ~
                                              algorithm + n_threads_per_chain + num_chunks,
                                          data = transform( results,
                                                            num_chunks = ifelse( is.na(num_chunks),
                                                                                 0, num_chunks)),
                                          FUN = mean)
      summary_table <-  summary_table[order(summary_table$sec_per_step_slowest_chain), ]
      ##
      reference <-  summary_table$sec_per_step_slowest_chain[summary_table$algorithm == "AD_Stan"]
      ##
      message(paste0( "\n", device, ", N = ", N,
                      ": best of each implementation (speed-up vs AD_Stan 4x1):"))
      for (algorithm in unique(summary_table$algorithm)) {
            best <-  summary_table[summary_table$algorithm == algorithm, ][1, ]
            message(paste0( "  ", algorithm, " ", n_chains_burnin, "x", best$n_threads_per_chain,
                            ", N_chunks ", best$num_chunks, ": ",
                            formatC(best$sec_per_step_slowest_chain, format = "e", digits = 3), " s/step (",
                            formatC(reference / best$sec_per_step_slowest_chain, format = "f", digits = 2), "x)"))
      }

}
























