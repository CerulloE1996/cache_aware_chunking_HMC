##
## R_fn_run_ps_1_burnin_benchmark.R
##
## -| ------------------------------    Pilot study 1 (burn-in) - optimal (num_chunks_burnin, n_threads_WCP_burnin) for BayesMVP ----------
##
## The sampling benchmark (ps_1_optimizing_N_chunks.R) times the post-burn-in sampler, where many chains run in
## parallel and WCP = 1 is best. Burn-in is different: only n_chains_burnin (e.g. 4) chains, advanced ONE iter at
## a time, so its wall time is the LATENCY of one iter - which within-chain parallelism (WCP > 1) can cut.
##
## These functions time exactly what the burn-in runs every iter: the persistent burn-in worker
## (fn_create_persistent_burnin_worker + fn_persistent_burnin_run_one_iter), i.e. all burn-in chains advanced one
## joint diffusion-HMC iter in parallel, each chain's lp/grad split over n_threads_WCP OpenMP threads.
## eps is tiny and L fixed, so nothing diverges and every configuration does the same work:
##   - eps_main = 1e-5, tau_main = L * eps (L = 20, the value the sampling benchmark uses);
##   - tau_ii ~ U(0, 2 * tau) per chain and iter, exactly as in burn-in; the seeds are the same for every
##     configuration, so every configuration draws the SAME trajectory lengths.
## Only the C++ iter is timed (the R-side adaptation between iters does not depend on chunks/threads).
##


##
## ---- Model inputs for one N (initialised ONCE per N; num_chunks is then set per configuration, as in the sampling ps1):
##
R_fn_ps1_burnin_make_model_inputs <- function( y,
                                               n_tests,
                                               SIMD_vect_type,
                                               MCMC_seed = 1000
) {

        N <- nrow(y)
        ##
        init_object <- BayesMVP:::initialise_model( Model_type            = "LC_MVP",
                                         stream                = MCMC_seed,
                                         sample_nuisance       = TRUE,
                                         n_nuisance_override   = NULL,
                                         model_args_list       = list( y         = y,
                                                                       n_tests   = n_tests,
                                                                       n_class   = 2,
                                                                       N         = N,
                                                                       vect_type = SIMD_vect_type),
                                         compile               = TRUE,
                                         force_recompile       = FALSE,
                                         cmdstanr_model_fit_obj = NULL,
                                         Stan_data_list        = NULL,
                                         Stan_model_file_path  = NULL,
                                         Stan_cpp_user_header  = NULL,
                                         Stan_cpp_flags        = NULL,
                                         stanc_args            = NULL)
        ##
        Model_args_as_Rcpp_List <- init_object$Model_args_as_Rcpp_List
        Model_args_as_Rcpp_List$model_so_file  <- "none"     ## as R_fn_sample_model does for built-in models
        Model_args_as_Rcpp_List$json_file_path <- "none"
        ##
        n_params_main <- Model_args_as_Rcpp_List$n_params_main
        n_nuisance    <- Model_args_as_Rcpp_List$n_nuisance
        ##
        return(list( y                        = y,
                     N                        = N,
                     n_params_main            = n_params_main,
                     n_nuisance               = n_nuisance,
                     Model_args_as_Rcpp_List  = Model_args_as_Rcpp_List))

}
##
## ---- Time the burn-in iter for ONE (n_chains_burnin, num_chunks, n_threads_WCP) configuration:
##
R_fn_ps1_burnin_thread_limit <- function(n_total_threads = NULL) {

        machine_threads <- parallel::detectCores(logical = TRUE)
        if (length(machine_threads) != 1 || !is.finite(machine_threads) || machine_threads < 1) {
              stop("Cannot determine the machine's hardware thread count; burn-in PS1 cannot validate its thread budget.")
        }
        if (is.null(n_total_threads)) return(machine_threads)
        if (!is.numeric(n_total_threads) || length(n_total_threads) != 1 ||
            !is.finite(n_total_threads) || n_total_threads < 1 || n_total_threads != floor(n_total_threads)) {
              stop("n_total_threads must be one positive integer.")
        }
        return(min(n_total_threads, machine_threads))

}



R_fn_ps1_time_burnin_iters <- function( model_inputs,
                                             n_chains_burnin,
                                             num_chunks_burnin,
                                             n_threads_WCP_burnin,
                                             L_main = 20,
                                             eps_main = 1e-5,
                                             n_untimed_burnin_iter_before_timing = 0,
                                             n_timed_iters = 10,
                                             seed = 1,
                                             n_total_threads = NULL,
                                             share_tau_ii_across_chains_in_burnin = TRUE,
                                             burnin_TBB_pool_equals_n_chains = FALSE
) {

        if (!is.logical(share_tau_ii_across_chains_in_burnin) || length(share_tau_ii_across_chains_in_burnin) != 1 ||
            is.na(share_tau_ii_across_chains_in_burnin)) {
              stop("share_tau_ii_across_chains_in_burnin must be TRUE or FALSE.")
        }
        if (!is.logical(burnin_TBB_pool_equals_n_chains) || length(burnin_TBB_pool_equals_n_chains) != 1 ||
            is.na(burnin_TBB_pool_equals_n_chains)) {
              stop("burnin_TBB_pool_equals_n_chains must be TRUE or FALSE.")
        }
        thread_limit <- R_fn_ps1_burnin_thread_limit(n_total_threads)
        for (thread_setting in list(n_chains_burnin = n_chains_burnin,
                                     n_threads_WCP_burnin = n_threads_WCP_burnin)) {
              if (!is.numeric(thread_setting) || length(thread_setting) != 1 ||
                  !is.finite(thread_setting) || thread_setting < 1 || thread_setting != floor(thread_setting)) {
                    stop("n_chains_burnin and n_threads_WCP_burnin must each be one positive integer.")
              }
        }
        if (n_chains_burnin * n_threads_WCP_burnin > thread_limit) {
              stop("Requested ", n_chains_burnin, " burn-in chains x ", n_threads_WCP_burnin,
                   " WCP threads = ", n_chains_burnin * n_threads_WCP_burnin,
                   " threads, exceeding the effective machine/configured limit of ", thread_limit, ".")
        }

        n_params_main <- model_inputs$n_params_main
        n_nuisance    <- model_inputs$n_nuisance
        ##
        ## ---- Model args with this chunk count:
        ##
        Model_args_as_Rcpp_List <- model_inputs$Model_args_as_Rcpp_List
        Model_args_as_Rcpp_List$Model_args_ints[4] <- num_chunks_burnin
        ##
        ## ---- eps tiny, trajectory length fixed at L leapfrog steps (on average; tau_ii ~ U(0, 2 tau) as in burn-in):
        ##
        EHMC_args_as_Rcpp_List <- init_EHMC_args_as_Rcpp_List( diffusion_HMC            = TRUE,
                                                               diffusion_HMC_integrator = "kick_flow_kick")
        EHMC_args_as_Rcpp_List$share_tau_ii_across_chains <- share_tau_ii_across_chains_in_burnin
        EHMC_args_as_Rcpp_List$eps_main <- eps_main
        EHMC_args_as_Rcpp_List$tau_main <- L_main * eps_main
        EHMC_args_as_Rcpp_List$eps_us   <- eps_main
        EHMC_args_as_Rcpp_List$tau_us   <- L_main * eps_main
        ##
        ## ---- Metric: use the same default dense-main metric initialisation as the original sampling PS1.
        ##      Do not invent an artificial nuisance mass or centre for the benchmark.
        ##
        EHMC_Metric_as_Rcpp_List <- init_EHMC_Metric_as_Rcpp_List( n_params_main     = n_params_main,
                                                                    n_nuisance        = n_nuisance,
                                                                    metric_shape_main = "dense")
        ##
        ## ---- Inits (as in the sampling ps1): corrs 0.01, betas -1 (class 1) / +1 (class 2), prevalence 0.20; nuisance 0.01:
        ##
        theta_main_vectors_all_chains_input_from_R <- matrix(0.01, nrow = n_params_main, ncol = n_chains_burnin)
        n_corrs            <- 2 * choose(ncol(model_inputs$y), 2)
        n_covariates_total <- n_params_main - n_corrs - 1
        theta_main_vectors_all_chains_input_from_R[(n_corrs + 1):(n_corrs + n_covariates_total / 2), ] <- -1
        theta_main_vectors_all_chains_input_from_R[(n_corrs + 1 + n_covariates_total / 2):(n_corrs + n_covariates_total), ] <- +1
        theta_main_vectors_all_chains_input_from_R[n_params_main, ] <- atanh(2 * 0.20 - 1)
        theta_us_vectors_all_chains_input_from_R <- matrix(0.01, nrow = n_nuisance, ncol = n_chains_burnin)
        ##
        ## ---- Same thread pool as R_fn_sample_model sets before burn-in: n_threads_WCP_burnin x n_chains_burnin threads,
        ##      or exactly n_chains_burnin with burnin_TBB_pool_equals_n_chains = TRUE (each chain's OpenMP team still
        ##      gets its n_threads_WCP_burnin threads):
        ##
        RcppParallel::setThreadOptions(numThreads = if (burnin_TBB_pool_equals_n_chains) n_chains_burnin else n_threads_WCP_burnin * n_chains_burnin)
        ##
        burnin_worker_pointer <- BayesMVP:::fn_create_persistent_burnin_worker( n_threads_R              = n_chains_burnin,
                                                                                 partitioned_HMC_R        = FALSE,
                                                                                 diffusion_HMC_R          = TRUE,
                                                                                 Model_type_R             = "LC_MVP",
                                                                                 sample_nuisance_R        = TRUE,
                                                                                 force_autodiff_R         = FALSE,
                                                                                 force_PartialLog_R       = FALSE,
                                                                                 multi_attempts_R         = TRUE,
                                                                                 y_Eigen_R                = model_inputs$y,
                                                                                 Model_args_as_Rcpp_List  = Model_args_as_Rcpp_List,
                                                                                 EHMC_args_as_Rcpp_List   = EHMC_args_as_Rcpp_List,
                                                                                 EHMC_Metric_as_Rcpp_List = EHMC_Metric_as_Rcpp_List,
                                                                                 n_threads_WCP            = n_threads_WCP_burnin)
        BayesMVP:::fn_persistent_burnin_update_adaptation( burnin_worker_pointer,
                                                           EHMC_args_as_Rcpp_List,
                                                           EHMC_Metric_as_Rcpp_List)
        BayesMVP:::fn_persistent_burnin_set_theta( burnin_worker_pointer,
                                                   theta_main_vectors_all_chains_input_from_R,
                                                   theta_us_vectors_all_chains_input_from_R)
        ##
        ## ---- Optional untimed HMC iters, then the timed iters.
        ##      The benchmark driver sets the untimed count to zero; worker/thread-pool construction above is already outside the timer.
        ##
        sec_per_iter <- rep(NA_real_, n_timed_iters)
        n_divs <- 0
        ##
        for (iter_index in seq_len(n_untimed_burnin_iter_before_timing + n_timed_iters)) {
          
                  iter_start_time <- proc.time()[["elapsed"]]
                  burnin_iter_outputs <- BayesMVP:::fn_persistent_burnin_run_one_iter(  burnin_worker_pointer,
                                                                                        seed + iter_index,
                                                                                        iter_index)
                  iter_wall_time <- proc.time()[["elapsed"]] - iter_start_time
                  ##
                  if (iter_index > n_untimed_burnin_iter_before_timing) {
                      
                          sec_per_iter[iter_index - n_untimed_burnin_iter_before_timing] <- iter_wall_time
                          n_divs_this_iter <- sum(burnin_iter_outputs$other_main_out_vector_all_chains_output_to_R[2, ])
                          n_divs <- n_divs + n_divs_this_iter
                          ##
                          if (n_divs_this_iter > 0) {
                                stop("Burn-in PS1 encountered ", n_divs_this_iter,
                                     " divergent transition(s) in timed iter ",
                                     iter_index - n_untimed_burnin_iter_before_timing,
                                     ". The configuration is invalid for a timing comparison.")
                          }
                          
                  }
              
        }
        rm(burnin_worker_pointer) ; invisible(gc(verbose = FALSE))
        ##
        return(list( sec_per_iter   = sec_per_iter,
                     total_timed_seconds     = sum(sec_per_iter),
                     n_divs = n_divs))

}


##
## ---- Build and validate the thread/chunk grid before any model initialization: ----------------------------------------------------
##
R_fn_ps1_burnin_make_configuration_grid <-  function( burnin_benchmark_settings,
                                                      N_vec
) {

        thread_limit <-  R_fn_ps1_burnin_thread_limit(n_total_threads = burnin_benchmark_settings[["n_total_threads"]])
        n_chains_burnin_vec <-  burnin_benchmark_settings[["n_chains_burnin_vec"]]
        if (!is.numeric(n_chains_burnin_vec) || length(n_chains_burnin_vec) == 0 ||
            any(!is.finite(n_chains_burnin_vec)) || any(n_chains_burnin_vec < 1) ||
            any(n_chains_burnin_vec != floor(n_chains_burnin_vec)))
        {
            stop("n_chains_burnin_vec must contain positive integers.")
        }
        if (!is.numeric(N_vec) || length(N_vec) == 0 || any(!is.finite(N_vec)) ||
            any(N_vec < 1) || any(N_vec != floor(N_vec)))
        {
            stop("N_vec must contain positive integers.")
        }
        configuration_grids_by_N <-  list()
        ##
        for (N in unique(x = N_vec)) {
              N_key <- as.character(N)
              if (is.null(burnin_benchmark_settings[["num_chunks_burnin_vec_given_N"]][[N_key]])) {
                    stop("No num_chunks_burnin candidates were supplied for N = ", N, ".")
              }
              if (is.null(burnin_benchmark_settings[["n_threads_WCP_burnin_vec_given_N"]][[N_key]])) {
                    stop("No n_threads_WCP_burnin candidates were supplied for N = ", N, ".")
              }
              if (is.null(burnin_benchmark_settings[["n_timed_iters_given_N"]][[N_key]])) {
                    stop("No n_timed_iters value was supplied for N = ", N, ".")
              }
              num_chunks_burnin_vec    <- burnin_benchmark_settings[["num_chunks_burnin_vec_given_N"]][[N_key]]
              n_threads_WCP_burnin_by_chains <- burnin_benchmark_settings[["n_threads_WCP_burnin_vec_given_N"]][[N_key]]
              ##
              ## ---- every valid configuration; invalid ones are skipped, not silently dropped (listed below):
              ##
              if (!is.list(n_threads_WCP_burnin_by_chains)) {
                    stop("WCP candidates for N = ", N, " must be a list keyed by burn-in chain count.")
              }
              configuration_grids_by_chains <- list()
              for (n_chains_burnin in n_chains_burnin_vec) {
                    chain_count_key <- as.character(n_chains_burnin)
                    n_threads_WCP_burnin_vec <- n_threads_WCP_burnin_by_chains[[chain_count_key]]
                    if (is.null(n_threads_WCP_burnin_vec)) {
                          stop("Missing WCP entry [[\"", N_key, "\"]][[\"", chain_count_key,
                               "\"]] in n_threads_WCP_burnin_vec_given_N.")
                    }
                    if (!is.numeric(n_threads_WCP_burnin_vec) || length(n_threads_WCP_burnin_vec) == 0 ||
                        any(!is.finite(n_threads_WCP_burnin_vec)) || any(n_threads_WCP_burnin_vec < 1) ||
                        any(n_threads_WCP_burnin_vec != floor(n_threads_WCP_burnin_vec))) {
                          stop("Supply positive integer WCP candidates for N = ", N, ", n_chains_burnin = ", n_chains_burnin, ".")
                    }
                    ## Reverse the expansion axes, then restore the columns, preserving the original first-axis-fastest run order.
                    configuration_grids_by_chains[[chain_count_key]] <- tidyr::expand_grid(
                          run_number           = seq_len(burnin_benchmark_settings$n_runs),
                          n_threads_WCP_burnin = n_threads_WCP_burnin_vec,
                          num_chunks_burnin    = num_chunks_burnin_vec,
                          n_chains_burnin      = n_chains_burnin) %>%
                        dplyr::select(n_chains_burnin, num_chunks_burnin, n_threads_WCP_burnin, run_number)
              }
              configuration_grid <- dplyr::bind_rows(configuration_grids_by_chains)
              is_valid_configuration <- (configuration_grid$n_threads_WCP_burnin <= configuration_grid$num_chunks_burnin) &
                                        (configuration_grid$n_chains_burnin * configuration_grid$n_threads_WCP_burnin <= thread_limit)
              skipped_configurations <- configuration_grid %>%
                  dplyr::filter(!is_valid_configuration) %>%
                  dplyr::distinct(n_chains_burnin, num_chunks_burnin, n_threads_WCP_burnin)
              if (nrow(skipped_configurations) > 0) {
                    message("N = ", N, ": skipping ", nrow(skipped_configurations), " configuration(s) with n_threads_WCP > num_chunks ",
                            "(threads with no chunk) or n_chains x n_threads_WCP > ", thread_limit, " threads.")
              }
              configuration_grid <- dplyr::filter(.data = configuration_grid, is_valid_configuration)
              if (nrow(configuration_grid) == 0) {
                    stop("No valid burn-in PS1 configurations for N = ", N, " within the thread limit of ", thread_limit, ".")
              }
              configuration_grid <- dplyr::mutate(.data = configuration_grid, N = .env$N,
                                                    n_threads_total = .data$n_chains_burnin * .data$n_threads_WCP_burnin)
              configuration_grids_by_N[[N_key]] <- dplyr::arrange(.data = configuration_grid, .data$run_number)
        }
        ##
        return(list(configuration_grid = dplyr::bind_rows(configuration_grids_by_N),
                    thread_limit       = thread_limit))
}

##
## ---- Plot appearance: explicit settings shared by the preview and saved-results plots: --------------------------------------------
##
R_fn_ps1_burnin_plot_settings <-  function(plot_settings = list()) {
  
        defaults <-  list(base_size              = 12,
                          axis_text_size         = 11,
                          axis_title_size        = 13,
                          axis_text_x_angle      = 45,
                          axis_text_x_hjust      = 1,
                          axis_tick_length_mm    = 2,
                          axis_tick_linewidth    = 0.4,
                          strip_text_size        = 12,
                          legend_text_size       = 11,
                          legend_title_size      = 12,
                          plot_title_size        = 14,
                          plot_subtitle_size     = 11,
                          plot_caption_size      = 10,
                          legend_position        = "bottom",
                          point_size             = 2,
                          line_width             = 0.6,
                          facet_ncol             = 3,
                          facet_scales           = "free_y",
                          preview_facet_scales   = "free_x",
                          y_scale                = "log10",
                          y_limits               = NULL,
                          y_n_breaks             = 6,
                          y_break_increment      = NULL,
                          y_breaks               = NULL,
                          y_label_digits         = NULL,
                          width_inches           = 14,
                          height_inches          = NULL,
                          dpi                    = 150)
        ##
        unknown_settings <-  setdiff(x = names(plot_settings), y = names(defaults))
        ##
        if (length(unknown_settings) > 0) { 
             stop("Unknown plot setting(s): ", paste(unknown_settings, collapse = ", "))
        }
        ##
        return(utils::modifyList(x = defaults, val = plot_settings, keep.null = TRUE))
        
}

##
R_fn_ps1_burnin_plot_theme <-  function(plot_settings = list()) {
        require(ggplot2)
        ##
        plot_settings <-  R_fn_ps1_burnin_plot_settings( plot_settings = plot_settings)
        ##
        theme_bw(base_size =  plot_settings$base_size) +
        theme(axis.text         = element_text(size = plot_settings$axis_text_size),
              axis.text.x       = element_text(angle = plot_settings$axis_text_x_angle,
                                               hjust = plot_settings$axis_text_x_hjust),
              axis.title        = element_text(size = plot_settings$axis_title_size),
              axis.ticks.length = grid::unit(x = plot_settings$axis_tick_length_mm, units = "mm"),
              axis.ticks        = element_line(linewidth = plot_settings$axis_tick_linewidth),
              strip.text        = element_text(size = plot_settings$strip_text_size),
              legend.text       = element_text(size = plot_settings$legend_text_size),
              legend.title      = element_text(size = plot_settings$legend_title_size),
              legend.position   = plot_settings$legend_position,
              plot.title        = element_text(size = plot_settings$plot_title_size),
              plot.subtitle     = element_text(size = plot_settings$plot_subtitle_size),
              plot.caption      = element_text(size = plot_settings$plot_caption_size))
}

##
## ---- Preview WCP x chains without constructing a model or running MCMC: ------------------------------------------------------------
##
R_fn_ps1_burnin_plot_thread_grid <-  function( burnin_benchmark_settings,
                                               N_vec,
                                               plot_settings = list()
) {
  
        require(ggplot2)
        ##
        plot_settings <- R_fn_ps1_burnin_plot_settings( plot_settings = plot_settings)
        ##
        grid_outputs <-   R_fn_ps1_burnin_make_configuration_grid(
                          burnin_benchmark_settings = burnin_benchmark_settings,
                          N_vec                     = N_vec)
        ## Collapse chunk counts/repetitions: each point is a runnable (N, chains, WCP) combination.
        thread_grid <-  unique(x = grid_outputs$configuration_grid[, c("N",
                                                                       "n_chains_burnin",
                                                                       "n_threads_WCP_burnin",
                                                                       "n_threads_total")])
        ##
        for (variable_name in c("N", "n_chains_burnin", "n_threads_WCP_burnin")) {
              numeric_values <- as.numeric(as.character(thread_grid[[variable_name]]))
              thread_grid[[variable_name]] <- factor(x = numeric_values, levels = sort(unique(numeric_values)))
        }
        ##
        ggplot(data = thread_grid,
               mapping = aes(x      = n_threads_WCP_burnin,
                             y      = n_threads_total,
                             colour = n_chains_burnin,
                             group  = n_chains_burnin)) +
            geom_hline(yintercept = grid_outputs$thread_limit, linetype = "dashed") +
            geom_line(linewidth = plot_settings$line_width) +
            geom_point(size = plot_settings$point_size) +
            scale_x_discrete(drop = TRUE) +
            facet_wrap(facets = ~ N, scales = plot_settings$preview_facet_scales,
                       ncol = plot_settings$facet_ncol, labeller = label_both) +
            scale_y_continuous(limits = c(0, grid_outputs$thread_limit),
                               breaks = sort(x = unique(x = c(pretty(x = c(0, grid_outputs$thread_limit)),
                                                               grid_outputs$thread_limit)))) +
            labs(title    = paste0(burnin_benchmark_settings$device, ": burn-in thread grid"),
                 subtitle = paste0("Total threads = chains x WCP; dashed line = effective limit (",
                                   grid_outputs$thread_limit, ")."),
                 x        = "WCP threads per chain",
                 y        = "Total threads",
                 colour   = "Burn-in chains",
                 caption  = "Only combinations supported by at least one selected chunk count are shown.") +
            R_fn_ps1_burnin_plot_theme(plot_settings = plot_settings)
}

##
## ---- Run the whole grid (N x n_chains_burnin x num_chunks_burnin x n_threads_WCP_burnin x runs) and save per N:
##
R_fn_run_ps_1_burnin_benchmark <- function( y_list,
                                            N_vec,
                                            n_tests,
                                            burnin_benchmark_settings,
                                            output_dir
) {

        ## Shared trajectory lengths are the PS1 burn-in default; an explicit FALSE retains independent draws.
        share_tau_ii_across_chains_in_burnin <- burnin_benchmark_settings[["share_tau_ii_across_chains_in_burnin"]]
        if (is.null(share_tau_ii_across_chains_in_burnin)) share_tau_ii_across_chains_in_burnin <- TRUE
        if (!is.logical(share_tau_ii_across_chains_in_burnin) || length(share_tau_ii_across_chains_in_burnin) != 1 ||
            is.na(share_tau_ii_across_chains_in_burnin)) {
              stop("share_tau_ii_across_chains_in_burnin must be TRUE or FALSE.")
        }
        ## TBB pool of exactly n_chains_burnin threads (as R_fn_sample_model with burnin_TBB_pool_equals_n_chains = TRUE);
        ## absent/FALSE = n_threads_WCP_burnin x n_chains_burnin threads:
        burnin_TBB_pool_equals_n_chains <- burnin_benchmark_settings[["burnin_TBB_pool_equals_n_chains"]]
        if (is.null(burnin_TBB_pool_equals_n_chains)) burnin_TBB_pool_equals_n_chains <- FALSE
        if (!is.logical(burnin_TBB_pool_equals_n_chains) || length(burnin_TBB_pool_equals_n_chains) != 1 ||
            is.na(burnin_TBB_pool_equals_n_chains)) {
              stop("burnin_TBB_pool_equals_n_chains must be TRUE or FALSE.")
        }
        ##
        grid_outputs <- R_fn_ps1_burnin_make_configuration_grid(
              burnin_benchmark_settings = burnin_benchmark_settings,
              N_vec                     = N_vec)
        thread_limit <- grid_outputs$thread_limit
        if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
        all_N_results <- list()
        ##
        ## ---- loop over POSITIONS, so the i-th N always goes with the i-th dataset:
        ##
        for (dataset_position in seq_along(N_vec)) {

              N <- N_vec[dataset_position]
              y_for_this_N <- y_list[[dataset_position]]
              if (!is.matrix(y_for_this_N) || nrow(y_for_this_N) != N) {
                    stop("dataset ", dataset_position, " does not have N = ", N, " rows - check N_vec_to_benchmark against N_vec.")
              }
              model_inputs <- R_fn_ps1_burnin_make_model_inputs( y              = y_for_this_N,
                                                                 n_tests        = n_tests,
                                                                 SIMD_vect_type = burnin_benchmark_settings$SIMD_vect_type)
              n_timed_iters <- burnin_benchmark_settings[["n_timed_iters_given_N"]][[as.character(N)]]
              configuration_grid <- dplyr::filter(.data = grid_outputs$configuration_grid, .data$N == .env$N)
              ##
              ## ---- runs in the OUTER position so a partial run still covers the whole grid once:
              ##
              configuration_grid <- dplyr::arrange(.data = configuration_grid, .data$run_number)
              per_configuration_rows <- vector("list", nrow(configuration_grid))
              ##
              for (configuration_index in seq_len(nrow(configuration_grid))) {
                    configuration <- dplyr::slice(.data = configuration_grid, configuration_index)
                    timing_outputs <- R_fn_ps1_time_burnin_iters( model_inputs          = model_inputs,
                                                                       n_chains_burnin       = configuration$n_chains_burnin,
                                                                       num_chunks_burnin     = configuration$num_chunks_burnin,
                                                                       n_threads_WCP_burnin  = configuration$n_threads_WCP_burnin,
                                                                       L_main                = burnin_benchmark_settings$L_main,
                                                                       eps_main              = burnin_benchmark_settings$eps_main,
                                                                       n_untimed_burnin_iter_before_timing   = burnin_benchmark_settings$n_untimed_burnin_iter_before_timing,
                                                                       n_timed_iters    = n_timed_iters,
                                                                       seed                  = 1000 * configuration$run_number,
                                                                       n_total_threads       = thread_limit,
                                                                       share_tau_ii_across_chains_in_burnin = share_tau_ii_across_chains_in_burnin,
                                                                       burnin_TBB_pool_equals_n_chains      = burnin_TBB_pool_equals_n_chains)
                    per_configuration_rows[[configuration_index]] <- tibble::tibble(
                          device                          = burnin_benchmark_settings$device,
                          N                               = N,
                          n_chains_burnin                 = configuration$n_chains_burnin,
                          num_chunks_burnin               = configuration$num_chunks_burnin,
                          n_threads_WCP_burnin            = configuration$n_threads_WCP_burnin,
                          n_threads_total                 = configuration$n_chains_burnin * configuration$n_threads_WCP_burnin,
                          run_number                      = configuration$run_number,
                          L_main                          = burnin_benchmark_settings$L_main,
                          n_timed_iters              = n_timed_iters,
                          total_timed_seconds             = timing_outputs$total_timed_seconds,
                          median_sec_per_iter    = median(timing_outputs$sec_per_iter),
                          mean_sec_per_iter      = mean(timing_outputs$sec_per_iter),
                          n_divs         = timing_outputs$n_divs)
                    progress_message <-  paste0( "N = ", formatC(x = N, format = "d", width = 6),
                                                  " | chains ", formatC(x = configuration$n_chains_burnin, format = "d", width = 2),
                                                  " | chunks ", formatC(x = configuration$num_chunks_burnin, format = "d", width = 4),
                                                  " | n_WCP (threads/chain) ",
                                                  formatC(x = configuration$n_threads_WCP_burnin, format = "d", width = 3),
                                                  " | total threads ",
                                                  formatC( x = configuration$n_chains_burnin * configuration$n_threads_WCP_burnin,
                                                           format = "d", width = 3),
                                                  " | run ", configuration$run_number,
                                                  " | ", formatC( x = timing_outputs$total_timed_seconds,
                                                                  format = "f", digits = 4, decimal.mark = "."),
                                                  " s total | ", formatC( x = mean(x = timing_outputs$sec_per_iter),
                                                                          format = "f", digits = 4, decimal.mark = "."),
                                                  " s / iter (mean)",
                                                  if (timing_outputs$n_divs > 0) "   <-- DIVERGENCES (eps too big?)" else "")
                    ##
                    cat( NicoStan::colourise( text = progress_message,
                                              fg = if (timing_outputs$n_divs > 0) "red" else "cyan"),
                         "\n", sep = "")
              }
              ##
              N_results <- dplyr::bind_rows(per_configuration_rows)
              results_file_path <- file.path(output_dir, paste0(if (burnin_benchmark_settings$device == "Laptop") "Laptop_" else "HPC_",
                                                                "ps1_burnin_benchmark_N", N, "_L", burnin_benchmark_settings$L_main,
                                                                "_n_runs", burnin_benchmark_settings$n_runs, ".rds"))
              saveRDS(N_results, results_file_path)
              message(NicoStan::colourise(text = paste0("saved: ", results_file_path), fg = "green"))
              all_N_results[[as.character(N)]] <- N_results

        }
        ##
        return(dplyr::bind_rows(all_N_results))

}


##
## ---- Read available saved rows, allowing missing files, configurations and repetitions: ------------------------------------------
##
R_fn_ps1_burnin_read_results <-  function(  output_dir,
                                            N_vec,
                                            device_prefixes,
                                            L_main,
                                            n_runs
) {
        saved_results <-  list()
        required_columns <-  c("device", "N", "n_chains_burnin", "num_chunks_burnin", "n_threads_WCP_burnin",
                               "n_threads_total", "run_number", "L_main", "n_timed_iters", "total_timed_seconds",
                               "mean_sec_per_iter", "n_divs")
        ##
        for (device_prefix in device_prefixes)
        {
            for (N in N_vec)
            {
                results_file_path <-  file.path(output_dir, paste0(device_prefix, "ps1_burnin_benchmark_N", N,
                                                                    "_L", L_main, "_n_runs", n_runs, ".rds"))
                if (!file.exists(results_file_path))
                {
                    message("Missing (skipped): ", results_file_path)
                    next
                }
                results <-  tryCatch(expr = readRDS(file = results_file_path),
                                     error = function(error)
                                     {
                                         warning("Could not read ", results_file_path, ": ", conditionMessage(error), call. = FALSE)
                                         NULL
                                     })
                if (is.null(results)) next
                ##
                ## Older RDS files keep their original column names. Translate only in memory;
                ## files already using the current names pass through unchanged.
                ##
                if (is.data.frame(results))
                {
                    previous_column_names <-  c(n_timed_iters       = "n_timed_iterations",
                                                 mean_sec_per_iter   = "mean_seconds_per_iteration",
                                                 median_sec_per_iter = "median_seconds_per_iteration",
                                                 n_divs              = "n_divergent_transitions")
                    for (current_column_name in names(previous_column_names))
                    {
                        previous_column_name <-  previous_column_names[[current_column_name]]
                        if (!current_column_name %in% names(results) && previous_column_name %in% names(results))
                        {
                            results[[current_column_name]] <-  results[[previous_column_name]]
                        }
                    }
                }
                if (!is.data.frame(results) || !all(required_columns %in% names(results)))
                {
                    warning("Unexpected burn-in results format (skipped): ", results_file_path, call. = FALSE)
                    next
                }
                numeric_columns <-  setdiff(x = required_columns, y = "device")
                if (!all(vapply(X = results[numeric_columns], FUN = is.numeric, FUN.VALUE = logical(1))))
                {
                    warning("Non-numeric burn-in result fields (skipped): ", results_file_path, call. = FALSE)
                    next
                }
                usable_rows <-  complete.cases(results[required_columns]) &
                                rowSums(!is.finite(as.matrix(results[numeric_columns]))) == 0 &
                                results$n_timed_iters > 0 & results$total_timed_seconds > 0 &
                                results$mean_sec_per_iter > 0 & results$n_divs == 0 &
                                results$N == N & results$L_main == L_main
                usable_rows[is.na(usable_rows)] <-  FALSE
                if (any(!usable_rows)) message("Skipped ", sum(!usable_rows), " incomplete/invalid timing row(s): ", results_file_path)
                results <-  dplyr::filter(.data = tibble::as_tibble(x = results), usable_rows)
                if (nrow(results) == 0) next
                if (!"median_sec_per_iter" %in% names(results)) results$median_sec_per_iter <-  NA_real_
                saved_results[[results_file_path]] <-  dplyr::select(.data = results, dplyr::all_of(c(required_columns, "median_sec_per_iter")))
                message("Loaded ", nrow(results), " completed configuration/run row(s): ", results_file_path)
            }
        }
        ##
        if (length(saved_results) == 0) return(tibble::tibble())
        return(dplyr::bind_rows(saved_results))
}

##
## ---- Summarise: per (N, n_chains_burnin), median over available runs; best configuration; speed-up over WCP = 1:
##
R_fn_ps1_burnin_summarise <- function(burnin_benchmark_results) {

        configuration_columns <-  c("device", "N", "n_chains_burnin", "num_chunks_burnin", "n_threads_WCP_burnin", "n_threads_total")
        configuration_summary <-  tibble::as_tibble(x = burnin_benchmark_results) %>%
            dplyr::filter(dplyr::if_all(dplyr::all_of(configuration_columns), ~ !is.na(.x))) %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(configuration_columns))) %>%
            dplyr::summarise(dplyr::across(c(total_timed_seconds, mean_sec_per_iter),
                                           ~ mean(x = .x, na.rm = TRUE)),
                             median_sec_per_iter = stats::median(x = .data$median_sec_per_iter, na.rm = TRUE),
                             n_runs_available = dplyr::n_distinct(.data$run_number, na.rm = TRUE),
                             .groups = "drop") %>%
            dplyr::filter(.data$n_runs_available > 0) %>%
            dplyr::mutate(sec_per_iter = .data$mean_sec_per_iter)
        ##
        ## ---- speed-up relative to the best single-threaded (WCP = 1) configuration of the same N and n_chains:
        ##
        if (any(configuration_summary$n_threads_WCP_burnin == 1)) {
              best_single_thread_seconds <- configuration_summary %>%
                  dplyr::filter(.data$n_threads_WCP_burnin == 1, !is.na(.data$total_timed_seconds)) %>%
                  dplyr::group_by(.data$device, .data$N, .data$n_chains_burnin) %>%
                  dplyr::summarise(best_WCP_1_total_timed_seconds = min(.data$total_timed_seconds), .groups = "drop")
              configuration_summary <- dplyr::left_join(x = configuration_summary, y = best_single_thread_seconds,
                                                          by = c("device", "N", "n_chains_burnin"))
        } else {
              configuration_summary <- dplyr::mutate(.data = configuration_summary, best_WCP_1_total_timed_seconds = NA_real_)
        }
        configuration_summary <- dplyr::mutate(.data = configuration_summary,
                                                speed_up_vs_WCP_1 = .data$best_WCP_1_total_timed_seconds / .data$total_timed_seconds)
        ## projected C++ kernel time of a 250-iter burn-in at this L (real burn-in L varies with the adaptation):
        configuration_summary <- dplyr::mutate(.data = configuration_summary, projected_sec_for_250_iter = 250 * .data$sec_per_iter)
        ##
        best_configuration_per_N <- configuration_summary %>%
            dplyr::arrange(.data$device, .data$N, .data$n_chains_burnin, .data$num_chunks_burnin, .data$n_threads_WCP_burnin) %>%
            dplyr::group_by(.data$device, .data$N, .data$n_chains_burnin) %>%
            dplyr::slice(which.min(.data$total_timed_seconds)) %>%
            dplyr::ungroup() %>%
            dplyr::arrange(.data$n_chains_burnin, .data$N, .data$device)
        ##
        return(list( configuration_summary     = dplyr::arrange(.data = configuration_summary, .data$N, .data$n_chains_burnin,
                                                                 .data$num_chunks_burnin, .data$n_threads_WCP_burnin),
                     best_configuration_per_N  = best_configuration_per_N))

}


##
## ---- Plot: seconds per burn-in iter vs n_threads_WCP_burnin, one line per num_chunks_burnin, one panel per N:
##
R_fn_ps1_burnin_plot <- function( burnin_benchmark_summary,
                                  x_variable = "n_threads_WCP_burnin",
                                  plot_settings = list()
) {

        require(ggplot2)
        ##
        plot_settings <- R_fn_ps1_burnin_plot_settings( plot_settings = plot_settings)
        ##
        configuration_summary <- tibble::as_tibble(x = burnin_benchmark_summary$configuration_summary)
        ##
        x_variable <- match.arg(arg = x_variable, choices = c("n_threads_WCP_burnin", "num_chunks_burnin"))
        line_variable <- if (x_variable == "n_threads_WCP_burnin") "num_chunks_burnin" else "n_threads_WCP_burnin"
        ##
        ## Numeric factor levels control axis labels, legend entries and panel order independently of row order.
        ##
        configuration_summary <- configuration_summary %>%
            dplyr::mutate(dplyr::across(c(N, n_chains_burnin, n_threads_WCP_burnin, num_chunks_burnin), function(values) {
                numeric_values <- as.numeric(x = as.character(x = values))
                factor(x = numeric_values, levels = sort(x = unique(x = numeric_values)))
            }))
        ##
        panel_order <- configuration_summary %>%
            dplyr::mutate(panel_row_number = dplyr::row_number()) %>%
            dplyr::arrange(device, N, n_chains_burnin) %>%
            dplyr::pull(panel_row_number)
        ##
        panel_labels <- paste0(configuration_summary$device, ": N = ", configuration_summary$N,
                               " (", configuration_summary$n_chains_burnin, " burn-in chains)")
        ##
        configuration_summary <- configuration_summary %>%
            dplyr::mutate(panel_label = factor(x = panel_labels, levels = unique(x = panel_labels[panel_order])),
                          plot_x = .data[[x_variable]],
                          plot_group = .data[[line_variable]])
        ##
        line_plot_data <- configuration_summary %>%
            dplyr::group_by(panel_label, plot_group) %>%
            dplyr::filter(dplyr::n() > 1) %>%
            dplyr::ungroup()
        y_scale <- match.arg(arg = plot_settings$y_scale, choices = c("log10", "linear"))
        ##
        ## Explicit tick values take precedence over a regular increment; NULL keeps automatic ticks.
        ## Limits may contain NA to keep that end automatic, e.g. c(0, NA).
        ##
        y_limits <- plot_settings$y_limits
        if (!is.null(y_limits)) {
              if (!is.numeric(y_limits) || length(y_limits) != 2 || any(is.infinite(y_limits)) ||
                  (!anyNA(y_limits) && y_limits[1] >= y_limits[2])) {
                    stop("y_limits must be NULL or two increasing numbers (NA allows an automatic endpoint).")
              }
              if (y_scale == "log10" && any(y_limits <= 0, na.rm = TRUE)) {
                    stop("Use y_scale = 'linear' to include zero in y_limits.")
              }
        }
        y_breaks <- plot_settings$y_breaks
        if (!is.null(y_breaks)) {
              if (!is.numeric(y_breaks) || any(!is.finite(y_breaks)) ||
                  (y_scale == "log10" && any(y_breaks <= 0))) {
                    stop("y_breaks must contain finite numbers, all positive when y_scale = 'log10'.")
              }
              y_breaks <- sort(unique(y_breaks))
        } else if (!is.null(plot_settings$y_break_increment)) {
              y_break_increment <- plot_settings$y_break_increment
              if (!is.numeric(y_break_increment) || length(y_break_increment) != 1 ||
                  !is.finite(y_break_increment) || y_break_increment <= 0) {
                    stop("y_break_increment must be a single positive number or NULL.")
              }
              y_breaks <- function(limits) {
                    first_tick <- ceiling(limits[1] / y_break_increment - 1e-9)
                    last_tick <- floor(limits[2] / y_break_increment + 1e-9)
                    if (first_tick > last_tick) return(numeric(0))
                    ticks <- seq(from = first_tick, to = last_tick) * y_break_increment
                    if (y_scale == "log10") ticks <- ticks[ticks > 0]
                    return(ticks)
              }
        } else {
              if (!is.numeric(plot_settings$y_n_breaks) || length(plot_settings$y_n_breaks) != 1 ||
                  !is.finite(plot_settings$y_n_breaks) || plot_settings$y_n_breaks < 2 ||
                  plot_settings$y_n_breaks != floor(plot_settings$y_n_breaks)) {
                    stop("y_n_breaks must be a single integer of at least 2.")
              }
              y_breaks <- ggplot2::waiver()
        }
        y_labels <- ggplot2::waiver()
        if (!is.null(plot_settings$y_label_digits)) {
              y_label_digits <- plot_settings$y_label_digits
              if (!is.numeric(y_label_digits) || length(y_label_digits) != 1 || !is.finite(y_label_digits) ||
                  y_label_digits < 0 || y_label_digits != floor(y_label_digits)) {
                    stop("y_label_digits must be a single nonnegative integer or NULL.")
              }
              y_labels <- function(values) formatC(x = values, format = "f", digits = y_label_digits)
        }
        ##
        ggplot( configuration_summary,
                aes( x      = plot_x,
                     y      = sec_per_iter,
                     colour = plot_group,
                     group  = plot_group)) +
          geom_line(data = line_plot_data, linewidth = plot_settings$line_width) +
          geom_point(size = plot_settings$point_size) +
          scale_x_discrete(drop = TRUE) +
          scale_colour_discrete(drop = FALSE) +
          scale_y_continuous(trans = if (y_scale == "log10") "log10" else "identity",
                             limits = y_limits, breaks = y_breaks, labels = y_labels,
                             n.breaks = if (inherits(y_breaks, "waiver")) plot_settings$y_n_breaks else NULL,
                             oob = scales::oob_keep, expand = expansion(mult = c(0, 0.05)),
                             guide = guide_axis(check.overlap = FALSE)) +
          facet_wrap(facets = ~ panel_label, scales = plot_settings$facet_scales, ncol = plot_settings$facet_ncol,
                     axes = "all_y", axis.labels = "all_y") +
          labs( x      = if (x_variable == "n_threads_WCP_burnin") "n_WCP (threads/chain)" else "Number of chunks",
                y      = "seconds per burn-in iter (mean over runs)",
                colour = if (line_variable == "num_chunks_burnin") "Number of chunks" else "n_WCP (threads/chain)",
                caption = "Available completed runs only; repetition counts are reported in n_runs_available.") +
          R_fn_ps1_burnin_plot_theme( plot_settings = plot_settings)

}






















