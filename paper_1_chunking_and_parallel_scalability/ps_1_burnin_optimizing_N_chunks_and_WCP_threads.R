
## ps_1_burnin_optimizing_N_chunks_and_WCP_threads.R
##
## PS1 (burn-in): the empirically optimal (num_chunks_burnin, n_threads_WCP_burnin) for BayesMVP's BURN-IN.
##
## ps_1_optimizing_N_chunks.R finds the best settings for SAMPLING (many chains in parallel -> WCP = 1, many chains).
## Burn-in is different: n_chains_burnin (selected per device) chains, advanced ONE iteration at a time, so its wall time is the
## latency of a single iteration - which splitting each chain's lp/grad over n_threads_WCP OpenMP threads can cut.
##
## What is timed: the persistent burn-in worker the real burn-in calls every iteration (all burn-in chains, one joint
## diffusion-HMC iteration, in parallel), with
##   - eps = 1e-5 (tiny -> no divergences; ESS is irrelevant here),
##   - L = 20 leapfrog steps (tau = L * eps; the value the sampling ps1 uses), tau_ii ~ U(0, 2 tau) as in burn-in,
##   - the same seeds for every configuration, so every configuration draws the same trajectory lengths.
## Output: seconds per burn-in iteration for every (N, n_chains_burnin, num_chunks_burnin, n_threads_WCP_burnin).
##
## Run it on an otherwise IDLE machine - anything else running (e.g. a ps7 sweep) makes the timings meaningless.
##
# rm(list = ls())
# .rs.restartR()
## FALSE loads and plots saved results only. TRUE also runs the benchmark first.
run_burnin_benchmark <- TRUE
##
##
## ------- Set options / paths:  ------------------------------------------------------------------------------------------------------
##
{
      options(scipen = 999)
      n_total_threads <- parallel::detectCores()
      computer <- ifelse(n_total_threads > 16, "Local_HPC", "Laptop")
      ##
      if (computer == "Local_HPC") {
        algorithm_study_dir <- "/home/enzocerullo/Documents/Work/PhD_work/Alg_paper_analysis"
        pkg_R_dir           <- "/home/enzocerullo/Documents/Work/PhD_work/R_packages/NicoStan/inst/NicoStan/R"   ## NicoStan: the sampler R code moved here
        bayesmvp_ext_R_dir  <- "/home/enzocerullo/Documents/Work/PhD_work/R_packages/BayesMVP/inst/BayesMVP/R"   ## BayesMVP extension: hard-coded LC-MVP model inits
      } else {
        algorithm_study_dir <- "/home/enzo/Documents/Work/PhD_work/Alg_paper_analysis"
        pkg_R_dir           <- "/home/enzo/Documents/Work/PhD_work/R_packages/NicoStan/inst/NicoStan/R"   ## NicoStan: the sampler R code moved here
        bayesmvp_ext_R_dir  <- "/home/enzo/Documents/Work/PhD_work/R_packages/BayesMVP/inst/BayesMVP/R"   ## BayesMVP extension: hard-coded LC-MVP model inits
      }
      setwd(algorithm_study_dir)
      ##
      ps1_dir            <- file.path(algorithm_study_dir, "paper_1_chunking_and_parallel_scalability")
      burnin_output_dir  <- file.path(ps1_dir, "burnin_outputs")
}
##
## ------- Packages and functions (the live BayesMVP R code, as ps7 uses; the C++ comes from the installed BayesMVP):  ---------------
##
if (run_burnin_benchmark) {
      source(file.path(algorithm_study_dir, 
                       "0_utilities", 
                       "shared_configs", 
                       "load_R_packages.R"))
      {
          require(Rcpp)
          require(RcppParallel)
          require(NicoStan)
          require(BayesMVP)
          require(ggplot2)
          require(dplyr)
          require(tibble)
          require(tidyr)
      }
      ##
      ##
      ## ---- The package R files are NO LONGER sourced here: since the NicoStan split, initialise_model() and
      ##      friends look up objects that live inside the package namespace (e.g. ".nicostan_model_types", set up by the
      ##      model backend), so a copy sourced into the global environment shadows the installed function and then fails
      ##      with "object '.nicostan_model_types' not found". The three functions this runner uses - initialise_model(),
      ##      init_EHMC_args_as_Rcpp_List() and init_EHMC_Metric_as_Rcpp_List() - are all exported by the installed BayesMVP
      ##      (which loads NicoStan), exactly as ps7 uses them. The old sourcing loop is kept below for reference.
      ##
      # for (package_R_file in c("R_fns_misc_BayesMVP.R",
      #                          "R_fn_initialise_model.R",
      #                          "R_fn_burnin_helper_fns.R")) {
      #       source(file.path(pkg_R_dir, package_R_file))
      # }
      # source(file.path(bayesmvp_ext_R_dir, "R_fns_init_hard_coded_models.R"))
      ##
      source(file.path(ps1_dir, "R_fn_run_ps_1_burnin_benchmark.R"))
      source(file.path(algorithm_study_dir, "0_utilities", "shared_functions", "R_fn_sim_bin_COVID_19_LC_MVP_data.R"))
}
##
## ------- Data: the same datasets as the sampling ps1 (COVID-19 DGM, parallel-scaling N's, seed 123):  ------------------------------
##
if (run_burnin_benchmark) {
      true_vals_list <- R_fn_simulate_binary_LC_MVP_data_COVID_19( study_type              = "algorithm_parallel_scaling_tests",
                                                                   seed                    = 123,
                                                                   corr_force_positive_DGM = FALSE)
      N_vec   <- true_vals_list$N_vec
      y_list  <- true_vals_list$y_list
      n_tests <- true_vals_list$n_tests
      print(N_vec)
}
##
## ------- Benchmark settings. EDIT HERE:  -------------------------------------------------------------------------------------------
##
{
      burnin_benchmark_settings <- list()
      ##
      burnin_benchmark_settings$device <- computer
      ##
      ## ---- fixed trajectory: L leapfrog steps (tau = L * eps), tiny eps. L = 20 matches the sampling ps1:
      ##
      burnin_benchmark_settings$L_main   <- 20
      burnin_benchmark_settings$eps_main <- 0.00001
      ##
      ## ---- device-specific burn-in chain counts; select chunks/WCP separately for each count:
      ##
      if (computer == "Local_HPC") {
            
                burnin_benchmark_settings$n_chains_burnin_vec <- c(4, 8, 16)
                burnin_benchmark_settings$n_total_threads <- 180               ## 96 cores + SMT (as in the sampling ps1)
                ##
                burnin_benchmark_settings$SIMD_vect_type  <- "AVX512"
                ##
                ## Lists are keyed by N, then burn-in chain count. Total threads = chains x WCP.
                ## Four chains: WCP 24/45 gives 96/180 threads. Eight chains: WCP 12/22 gives 96/176 threads.
                ##
                burnin_benchmark_settings$n_threads_WCP_burnin_vec_given_N <- list(
                      "500"   = list("4" = c(1, 2, 4, 8),     ## 4*c(1, 2, 4, 8)    =  4 8 16 32
                                     "8" = c(1, 2, 4, 8),     ## 8*c(1, 2, 4, 8)    =  8 16 32 64
                                    "16" = c(1, 2, 4, 6, 8)), ## 6*c(1, 2, 4, 6, 8) = 16  32  64  96 128
                      ##
                      "2500"  = list("4" = c(1, 2, 4, 8),
                                     "8" = c(1, 2, 4, 8),
                                    "16" = c(1, 2, 4, 6, 8)),
                      ##
                      "10000" = list("4" = c(1, 2, 4, 8, 16, 24, 32, 44),
                                     "8" = c(1, 2, 4, 8, 12, 16, 22),
                                    "16" = c(1, 2, 4, 6, 8, 11)),
                      ##
                      "50000" = list("4" = c(1, 2, 4, 8, 16, 24, 32, 44),
                                     "8" = c(1, 2, 4, 8, 12, 16, 22),
                                    "16" = c(1, 2, 4, 6, 8, 11)))
                ##
                # 4*c(1, 2, 4, 8)
                # c(4, 8*c(1, 2, 4, 8, 12, 16)) ## Note:  at N = 2500, the max. n_chunks = 10, so the helper will skip (12, 16)
                # c(4, 8, 16*c(1, 2, 4, 6, 8, 11)) ## Note:  at N = 2500, the max. n_chunks = 10, so the helper will skip (11)
                # ##
                # 4*c(1, 2, 4, 8, 16, 24, 32, 44)
                # c(4, 8*c(1, 2, 4, 8, 12, 16, 22))
                # c(4, 8, 16*c(1, 2, 4, 6, 8, 11))
            
      } else { ## Laptop
        
                burnin_benchmark_settings$n_chains_burnin_vec <- c(4, 8)
                burnin_benchmark_settings$n_total_threads <- 16                ## 8 cores / 16 threads
                ##
                burnin_benchmark_settings$SIMD_vect_type  <- "AVX2"
                ##
                burnin_benchmark_settings$n_threads_WCP_burnin_vec_given_N <- list(
                      "500"   = list("4" = c(1, 2, 4),
                                     "8" = c(1, 2)),
                      ##
                      "2500"  = list("4" = c(1, 2, 4),
                                     "8" = c(1, 2)),
                      ##
                      "10000" = list("4" = c(1, 2, 4),
                                     "8" = c(1, 2)),
                      ##
                      "50000" = list("4" = c(1, 2, 4),
                                     "8" = c(1, 2)))
            
      }
      ##
      ## ---- chunk counts per N (configurations with n_threads_WCP > num_chunks are skipped - a thread needs >= 1 chunk):
      ##
      burnin_benchmark_settings$num_chunks_burnin_vec_given_N <- list( "500"   = c(1, 2, 4, 10),
                                                                       "2500"  = c(1, 2, 4, 10, 25),
                                                                       "10000" = c(1, 4, 10, 25, 50, 100),
                                                                       "50000" = c(1, 4, 10, 25, 50, 100, 250, 500))
      ##
      ## ---- burn-in iterations timed per configuration (ONLY burn-in is run - no sampling), and repeats.
      ##      n_untimed_burnin_iter_before_timing = extra burn-in iterations run BEFORE the timer starts (0 = time every one):
      ##
      burnin_benchmark_settings$n_timed_iters_given_N <- list(     "500"   = 200, 
                                                                   "2500"  = 40,
                                                                   "10000" = 10,
                                                                   "50000" = 10)
      burnin_benchmark_settings$n_untimed_burnin_iter_before_timing <- 0
      burnin_benchmark_settings$n_runs <- 3
      ##
      burnin_benchmark_settings$share_tau_ii_across_chains_in_burnin <- TRUE  ## FALSE restores independent trajectory draws.
      burnin_benchmark_settings$burnin_TBB_pool_equals_n_chains      <- TRUE  ## FALSE = TBB pool of n_chains x n_WCP threads (the old setting).
      ##
      ## ---- which N's to run:
      ##
      # N_vec_to_benchmark <- c(500, 2500, 10000, 50000)
      # # N_vec_to_benchmark <- N_vec
      N_vec_to_benchmark <- c(500, 2500, 10000)
}
##
## ------- Preview the thread grid (settings only; no model initialization or MCMC):  -------------------------------------------------
##
if (run_burnin_benchmark) {
      burnin_thread_plot_settings <- list(axis_text_size      = 11,
                                          axis_title_size     = 13,
                                          strip_text_size     = 12,
                                          axis_tick_length_mm = 2,
                                          axis_tick_linewidth = 0.4,
                                          facet_ncol          = 2)
      burnin_thread_grid_plot <- R_fn_ps1_burnin_plot_thread_grid(
            burnin_benchmark_settings = burnin_benchmark_settings,
            N_vec                     = N_vec_to_benchmark,
            plot_settings             = burnin_thread_plot_settings)
      print(x = burnin_thread_grid_plot)
}
##
## ------- Run (one .rds per N in outputs/burnin_benchmarks):  ------------------------------------------------------------------------
##
if (run_burnin_benchmark) {
      burnin_benchmark_results <- R_fn_run_ps_1_burnin_benchmark( y_list                    = y_list[match(N_vec_to_benchmark, N_vec)],
                                                                  N_vec                     = N_vec_to_benchmark,
                                                                  n_tests                   = n_tests,
                                                                  burnin_benchmark_settings = burnin_benchmark_settings,
                                                                  output_dir                = burnin_output_dir)
}
##
## ------- LOAD SAVED RESULTS AND PLOT: run this section independently; no benchmark/data setup needed -------------------------------
##
{
      ## These settings select saved files only. Missing devices/Ns and incomplete repetitions are allowed.
      burnin_analysis_study_dir <- if (parallel::detectCores() > 16) {
            "/home/enzocerullo/Documents/Work/PhD_work/Alg_paper_analysis"
      } else {
            "/home/enzo/Documents/Work/PhD_work/Alg_paper_analysis"
      }
      burnin_analysis_ps1_dir <- file.path(burnin_analysis_study_dir, "paper_1_chunking_and_parallel_scalability")
      ##
      burnin_analysis_output_dir <- file.path(burnin_analysis_ps1_dir, "burnin_outputs")
      burnin_analysis_N_vec <- c(500, 2500, 10000, 50000)
      burnin_analysis_device_prefixes <- c("HPC_", "Laptop_")
      burnin_analysis_L_main <- 20
      burnin_analysis_n_runs <- 3                    ## Filename suffix; fewer completed repetitions are OK.
      ##
      require(ggplot2)
      require(dplyr)
      require(tibble)
      require(tidyr)
      ##
      source(file = file.path(burnin_analysis_ps1_dir, 
                              "R_fn_run_ps_1_burnin_benchmark.R"))
      ##
      ## ---- Plot appearance: text sizes in points, tick length in mm, image dimensions in inches. ----------------------------------
      ##
      burnin_plot_settings <- list( base_size            = 20,
                                    axis_text_size       = 14,
                                    axis_title_size      = 16,
                                    ##
                                    axis_text_x_angle    = 45,
                                    axis_text_x_hjust    = 1,
                                    ##
                                    axis_tick_length_mm  = 2,
                                    axis_tick_linewidth  = 0.4,
                                    strip_text_size      = 16,
                                    ##
                                    legend_text_size     = 16,
                                    legend_title_size    = 16,
                                    legend_position      = "bottom",
                                    ##
                                    plot_title_size      = 16,
                                    plot_subtitle_size   = 12,
                                    plot_caption_size    = 12,
                                    ##
                                    point_size           = 5,
                                    line_width           = 1.00,
                                    ##
                                    facet_ncol           = 3,          ## Numeric chain order: 4, 8, 16.
                                    facet_scales         = "free",     ## "fixed", "free_x", "free_y", or "free".
                                    ##
                                    y_scale              = "linear",   ## Zero requires "linear"; "log10" is also supported.
                                    y_limits             = c(0, NA),   ## Automatic upper end; use c(0, 0.1) for a chosen range.
                                    y_n_breaks           = 6,          ## Approximate tick count PER PANEL when both overrides below are NULL.
                                    y_break_increment    = NULL,       ## Optional fixed spacing in seconds; overrides y_n_breaks.
                                    y_breaks             = NULL,       ## Optional explicit tick vector; overrides the count and increment.
                                    y_label_digits       = NULL,       ## Automatic decimal precision for each panel's tick spacing.
                                    ##
                                    width_inches         = 14,
                                    height_inches        = NULL,       ## NULL selects height from the number of panels.
                                    dpi                  = 150)
      ##
      ## Both plots start with the same settings; these copies can be edited separately.
      ##
      
      ## Recorded console output, retained as comments.
      ##       e: 41 × 14
      ##    device      N n_chains_burnin num_chunks_burnin n_threads_WCP_burnin n_threads_total total_timed_seconds mean_sec_per_iter
      ##    <chr>   <dbl>           <dbl>             <dbl>                <dbl>           <dbl>               <dbl>             <dbl>
      ##  1 Local_…  2500               8                10                    8              64               0.458            0.0115
      ##  2 Local_…  2500               4                10                    6              24               0.469            0.0117
      ##  3 Local_…  2500               4                 4                    4              16               0.471            0.0118
      ##  4 Local_…  2500               8                 4                    4              32               0.485            0.0121
      ##  5 Local_…  2500               4                10                    8              32               0.496            0.0124
      ##  6 Local_…  2500               4                10                    4              16               0.503            0.0126
      ##  7 Local_…  2500               8                10                    4              32               0.509            0.0127
      ##  8 Local_…  2500              16                10                    4              64               0.519            0.0130
      ##  9 Local_…  2500              16                 4                    4              64               0.549            0.0137
      ## 10 Local_…  2500               4                 5                    4              16               0.551            0.0138
      ## 11 Local_…  2500               8                 5                    4              32               0.552            0.0138
      ## 12 Local_…  2500              16                 5                    4              64               0.577            0.0144
      ## 13 Local_…  2500              16                10                    6              96               0.595            0.0149
      ## 14 Local_…  2500               8                 4                    2              16               0.619            0.0155
      ## 15 Local_…  2500               4                 4                    2               8               0.622            0.0156
      ## 16 Local_…  2500               4                10                    2               8               0.635            0.0159
      ## 17 Local_…  2500               8                 2                    2              16               0.638            0.0159
      ## 18 Local_…  2500               4                 2                    2               8               0.638            0.0159
      ## 19 Local_…  2500               8                10                    2              16               0.638            0.0160
      ## 20 Local_…  2500              16                 4                    2              32               0.656            0.0164
      ## 21 Local_…  2500              16                10                    2              32               0.662            0.0165
      ## 22 Local_…  2500              16                10                    8             128               0.670            0.0167
      ## 23 Local_…  2500              16                 2                    2              32               0.671            0.0168

      burnin_WCP_plot_settings    <- burnin_plot_settings
      burnin_chunks_plot_settings <- burnin_plot_settings
      ##
      burnin_benchmark_results <- R_fn_ps1_burnin_read_results(
            output_dir      = burnin_analysis_output_dir,
            N_vec           = burnin_analysis_N_vec,
            device_prefixes = burnin_analysis_device_prefixes,
            L_main          = burnin_analysis_L_main,
            n_runs          = burnin_analysis_n_runs)
      ##
      # if (nrow(burnin_benchmark_results) > 0) {
        
          burnin_benchmark_summary <- R_fn_ps1_burnin_summarise( burnin_benchmark_results = burnin_benchmark_results)
          ##
          str(burnin_benchmark_summary)
          config_summary <- burnin_benchmark_summary$configuration_summary ## %>% print(width = 9000, n = 10)
          config_summary %>% print(width = 9000, n = 10)
          config_summary
          ##
          config_summary %>% 
            dplyr::filter(N == 500) %>%
            dplyr::arrange(total_timed_seconds) %>% 
            print(width = 125, n = 25)
          ##
          config_summary %>% 
            dplyr::filter(N == 2500) %>%
            dplyr::arrange(total_timed_seconds) %>% 
            print(width = 125, n = 25)
          ##
          config_summary %>% 
            dplyr::filter(N == 10000) %>%
            dplyr::arrange(total_timed_seconds) %>% 
            print(width = 125, n = 25)
          ##
          config_summary %>% 
            dplyr::filter(N == 50000) %>%
            dplyr::arrange(total_timed_seconds) %>% 
            print(width = 125, n = 25)
          
          ##
          cat("\n==== best observed configuration per device / N / chain count (available runs only) ====\n")
          burnin_benchmark_summary$best_configuration_per_N %>%
            dplyr::select(device, N, n_chains_burnin, num_chunks_burnin, n_threads_WCP_burnin,
                          n_threads_total, n_runs_available, sec_per_iter, speed_up_vs_WCP_1, projected_sec_for_250_iter) %>%
            print(n = Inf, width = Inf)
          ##
          cat("\n==== all configurations ====\n")
          config_summary %>%
            dplyr::mutate(sec_per_iter = signif(x = sec_per_iter, digits = 3),
                          speed_up_vs_WCP_1 = round(x = speed_up_vs_WCP_1, digits = 2)) %>%
            print(n = Inf, width = Inf)
          ##
          ## ---- Burn-in Plot 1:
          ##
          burnin_benchmark_plot <- R_fn_ps1_burnin_plot( burnin_benchmark_summary = burnin_benchmark_summary,
                                                         plot_settings            = burnin_WCP_plot_settings)
          print(burnin_benchmark_plot)
          ##
          ## ---- Burn-in Plot 2:
          ##
          burnin_benchmark_chunk_plot <- R_fn_ps1_burnin_plot( burnin_benchmark_summary = burnin_benchmark_summary,
                                                               x_variable               = "num_chunks_burnin",
                                                               plot_settings            = burnin_chunks_plot_settings)
          print(x = burnin_benchmark_chunk_plot)
          ##
          burnin_analysis_panel_count <-  config_summary %>%
            dplyr::distinct(device, N, n_chains_burnin) %>%
            nrow()
          burnin_WCP_plot_height <- burnin_WCP_plot_settings$height_inches
          if (is.null(burnin_WCP_plot_height)) {
                burnin_WCP_plot_height <- max(4, 3 * ceiling(burnin_analysis_panel_count / burnin_WCP_plot_settings$facet_ncol))
          }
          burnin_chunks_plot_height <- burnin_chunks_plot_settings$height_inches
          if (is.null(burnin_chunks_plot_height)) {
                burnin_chunks_plot_height <- max(4, 3 * ceiling(burnin_analysis_panel_count / burnin_chunks_plot_settings$facet_ncol))
          }
          burnin_analysis_file_stem <- paste0("ps1_burnin_saved_results_", paste(burnin_analysis_device_prefixes, collapse = ""),
                                              "L", burnin_analysis_L_main, "_n_runs", burnin_analysis_n_runs)
          ##
          ggsave(filename = file.path(burnin_analysis_output_dir, paste0(burnin_analysis_file_stem, "_WCP.png")),
                 plot     = burnin_benchmark_plot,
                 width    = burnin_WCP_plot_settings$width_inches,
                 height   = burnin_WCP_plot_height,
                 dpi      = burnin_WCP_plot_settings$dpi)
          ##
          ggsave(filename = file.path(burnin_analysis_output_dir, paste0(burnin_analysis_file_stem, "_chunks.png")),
                 plot     = burnin_benchmark_chunk_plot,
                 width    = burnin_chunks_plot_settings$width_inches,
                 height   = burnin_chunks_plot_height,
                 dpi      = burnin_chunks_plot_settings$dpi)
          
      # } else {
      #       message("No usable saved burn-in results found for these file settings; nothing to plot yet.")
      # }
}






















