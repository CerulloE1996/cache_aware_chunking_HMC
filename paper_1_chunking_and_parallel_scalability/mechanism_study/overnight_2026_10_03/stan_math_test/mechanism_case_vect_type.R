##
## =====================================================================================================================================
## mechanism_case_vect_type.R (copy of mechanism_case.R, 3 Oct 2026)
##
## As mechanism_case.R, plus MECH_VECT_TYPE (AVX512, AVX2 or Stan): the vectorisation of BayesMVP's own maths
## functions (exp, log, log-sum-exp, tanh, Phi, log Phi, inv Phi), set in Model_args_strings before the case runs.
## "Stan" uses the stan::math scalar functions in the same hand-written log-posterior and gradient.
##
## Runs ONE Paper 1 sampling case with the main runner's own settings and case function (fn_paper1_run_case), so the sampler,
## model, data, step size and trajectory length are exactly those of the study. It is launched by run_mechanism_experiments.sh
## under pmc_stat, which counts the CPU hardware events of the whole process. The case is given by environment variables:
##
##     MECH_ALGORITHM          MD_BayesMVP, MD_BayesMVP_WCP, AD_Stan, AD_Stan_tape_chunked or AD_Stan_WCP
##     MECH_N                  500, 2500, 10000 or 50000
##     MECH_CHUNKS             number of chunks
##     MECH_CHAINS             number of chains
##     MECH_THREADS_PER_CHAIN  WCP threads per chain (1 for the non-WCP arms)
##     MECH_N_ITER             number of sampling iterations
##     MECH_SEED               seed
##     MECH_LABEL              label written to the results file
##     MECH_RESULTS_FILE       CSV file the elapsed sampling time is appended to
##
## Nothing is written to the study's output, cache or model folders: compiled Stan models go to mechanism_study/stan_model_cache.
##
{
      mechanism_algorithm <-  Sys.getenv("MECH_ALGORITHM")
      mechanism_N <-  as.numeric(Sys.getenv("MECH_N"))
      mechanism_chunks <-  as.numeric(Sys.getenv("MECH_CHUNKS"))
      mechanism_chains <-  as.numeric(Sys.getenv("MECH_CHAINS"))
      mechanism_threads_per_chain <-  as.numeric(Sys.getenv("MECH_THREADS_PER_CHAIN", unset = "1"))
      mechanism_n_iter <-  as.numeric(Sys.getenv("MECH_N_ITER"))
      mechanism_seed <-  as.numeric(Sys.getenv("MECH_SEED", unset = "1000"))
      mechanism_label <-  Sys.getenv("MECH_LABEL")
      mechanism_results_file <-  Sys.getenv("MECH_RESULTS_FILE")
}
##
## ---- Settings, helpers and simulated data from the main runner (evaluated up to its dry-run line; nothing is run there): -------------
##
{
      algorithm_study_dir <-  path.expand("~/Documents/Work/PhD_work/Alg_paper_analysis")
      paper1_dir <-  file.path(algorithm_study_dir, "paper_1_chunking_and_parallel_scalability")
      mechanism_dir <-  file.path(paper1_dir, "mechanism_study")
      setwd(paper1_dir)
      ##
      runner_expressions <-  parse(file = "alg_paper_1_chunking_WCP_par_scaling.R", keep.source = TRUE)
      for (expression_index in seq_along(runner_expressions)) {
            expression_text <-  paste(deparse(runner_expressions[[expression_index]]), collapse = " ")
            if (grepl("dry_run = TRUE", expression_text, fixed = TRUE)) break
            eval(runner_expressions[[expression_index]], envir = globalenv())
            if (expression_index == 1) run_benchmark <-  FALSE
      }
      ##
      paper1_plan <-  suppressMessages(fn_run_paper1_benchmark( settings       = paper1_settings,
                                                                dry_run        = TRUE,
                                                                Stan_model_obj = NULL))
      mechanism_settings <-  paper1_plan$settings
      ##
      ## Compiled Stan models for these experiments live in their own folder (the study's model cache is left untouched):
      mechanism_settings$stan_model_cache_dir <-  file.path(mechanism_dir, "stan_model_cache", Sys.info()[["nodename"]])
      dir.create(path = mechanism_settings$stan_model_cache_dir, recursive = TRUE, showWarnings = FALSE)
      ##
      study_output_dir <-  file.path(paper1_dir, "paper_1_computational_outputs", mechanism_settings$device)
      simulated_data <-  readRDS(file = file.path(study_output_dir, "COVID19_datasets.rds"))
      y <-  simulated_data$y_list[[match(mechanism_N, simulated_data$N_vec)]]
}
##
## ---- The case (same columns as the main runner's grid): ------------------------------------------------------------------------------
##
{
      mechanism_case <-  data.frame( device            = mechanism_settings$device,
                                     algorithm         = mechanism_algorithm,
                                     N                 = mechanism_N,
                                     num_chunks        = mechanism_chunks,
                                     n_threads         = mechanism_chains * mechanism_threads_per_chain,
                                     n_chains          = mechanism_chains,
                                     threads_per_chain = mechanism_threads_per_chain,
                                     n_iter            = mechanism_n_iter,
                                     n_iter_short_run  = 1,
                                     timing_method     = mechanism_settings$timing_method,
                                     run               = 1,
                                     base_seed         = mechanism_seed,
                                     seed              = mechanism_seed,
                                     benchmark_role    = "mechanism_study",
                                     execution_backend = if (grepl("^MD_", mechanism_algorithm)) "BayesMVP_NicoStan" else "NicoStan_BridgeStan",
                                     case_id           = 1,
                                     timing_estimator  = "single_run",
                                     timing_run_label  = mechanism_label,
                                     stringsAsFactors  = FALSE)
      ##
      if (grepl("^AD_Stan", mechanism_algorithm) && mechanism_algorithm != "AD_Stan") {
            mechanism_case$stan_chunk_size <-  ceiling(mechanism_N / mechanism_chunks)
      }
}
##
## ---- Runtime objects the case function needs (as built by the main runner), then the case itself: -----------------------------------
##
{
      if (grepl("^MD_", mechanism_algorithm)) {

            mechanism_runtime <-  list( bayesmvp_model_args = new.env(parent = emptyenv()),
                                        bayesmvp_fixed_L    = mechanism_settings$bayesmvp$L_main)
            ##
            ## ---- Vectorisation of BayesMVP's maths functions (MECH_VECT_TYPE): the model arguments are built
            ##      here, as in fn_paper1_run_case, with the vect_type strings set, and put in the runtime cache
            ##      that fn_paper1_run_case reads (it then only sets the chunk count and J_grad_option):
            mechanism_vect_type <-  Sys.getenv("MECH_VECT_TYPE", unset = "AVX512")
            stopifnot(mechanism_vect_type %in% c("AVX512", "AVX2", "Stan"))
            vect_type_model <-  BayesMVP:::initialise_model( Model_type = "LC_MVP", stream = mechanism_seed,
                                                             sample_nuisance = TRUE, n_nuisance_override = NULL,
                                                             model_args_list = list(y = y), compile = TRUE,
                                                             force_recompile = FALSE, cmdstanr_model_fit_obj = NULL,
                                                             Stan_data_list = NULL, Stan_model_file_path = NULL,
                                                             Stan_cpp_user_header = NULL, Stan_cpp_flags = NULL,
                                                             stanc_args = NULL)
            vect_type_args <-  vect_type_model$Model_args_as_Rcpp_List
            vect_type_args$Model_args_strings[c(1, 4:11)] <-  mechanism_vect_type
            message(paste0("\033[36m", "vect_type strings: ",
                           paste(vect_type_args$Model_args_strings[c(1, 4:11)], collapse = ", "), "\033[0m"))
            mechanism_runtime$bayesmvp_model_args[[as.character(mechanism_N)]] <-  vect_type_args

      } else {

            mechanism_models <-  fn_paper1_compile_stan_via_NicoStan( settings   = mechanism_settings,
                                                                      algorithms = mechanism_algorithm)
            mechanism_runtime <-  list( stan_models             = mechanism_models,
                                        stan_data               = stats::setNames(object = list(fn_paper1_stan_data(y = y)),
                                                                                  nm = as.character(mechanism_N)),
                                        stan_via_NicoStan_cache = new.env(parent = emptyenv()))

      }
      ##
      ## ---- Counter snapshots: pmc_stat appends the running totals when it receives SIGUSR1, so each sampling call below is the
      ##      difference of two snapshots (R start-up, package loading and model initialisation are excluded):
      fn_mechanism_snapshot <-  function() {
            pmc_stat_pid <-  as.integer(Sys.getenv("PMC_STAT_PID", unset = "0"))
            if (pmc_stat_pid > 0) tools::pskill(pid = pmc_stat_pid, signal = tools::SIGUSR1)
            umc_stat_pid <-  as.integer(Sys.getenv("UMC_STAT_PID", unset = "0"))   ## memory-controller counters, when umc_stat is used
            if (umc_stat_pid > 0) tools::pskill(pid = umc_stat_pid, signal = tools::SIGUSR1)
            Sys.sleep(0.25)
      }
      ##
      fn_mechanism_run <-  function(n_iter) {
            case_now <-  mechanism_case
            case_now$n_iter <-  n_iter
            fn_paper1_run_case( case     = case_now,
                                y        = y,
                                settings = mechanism_settings,
                                runtime  = mechanism_runtime)
      }
      ##
      mechanism_n_iter_short <-  as.numeric(Sys.getenv("MECH_N_ITER_SHORT", unset = "2"))
      ##
      invisible(fn_mechanism_run(n_iter = 1))                              ## untimed warm-up: model initialisation, first-touch memory
      fn_mechanism_snapshot()                                              ## snapshot 1
      mechanism_output_short <-  fn_mechanism_run(n_iter = mechanism_n_iter_short)
      fn_mechanism_snapshot()                                              ## snapshot 2
      mechanism_output_long <-  fn_mechanism_run(n_iter = mechanism_n_iter)
      fn_mechanism_snapshot()                                              ## snapshot 3
      ##
      message(paste0("\033[36mMechanism case ", mechanism_label, ": elapsed sampling time = ",
                     formatC(mechanism_output_short$elapsed_seconds, format = "f", digits = 3), " s (", mechanism_n_iter_short, " iterations), ",
                     formatC(mechanism_output_long$elapsed_seconds, format = "f", digits = 3), " s (", mechanism_n_iter, " iterations)\033[0m"))
      ##
      mechanism_result <-  data.frame( label                     = mechanism_label,
                                       device                    = mechanism_settings$device,
                                       algorithm                 = mechanism_algorithm,
                                       N                         = mechanism_N,
                                       num_chunks                = mechanism_chunks,
                                       n_chains                  = mechanism_chains,
                                       threads_per_chain         = mechanism_threads_per_chain,
                                       n_iter_short              = mechanism_n_iter_short,
                                       n_iter_long               = mechanism_n_iter,
                                       elapsed_seconds_short_run = mechanism_output_short$elapsed_seconds,
                                       elapsed_seconds_long_run  = mechanism_output_long$elapsed_seconds,
                                       stringsAsFactors          = FALSE)
      utils::write.table( x = mechanism_result, file = mechanism_results_file, sep = ",", row.names = FALSE,
                          col.names = !file.exists(mechanism_results_file), append = file.exists(mechanism_results_file))
}
























