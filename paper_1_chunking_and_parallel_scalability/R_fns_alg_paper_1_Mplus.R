##
## ========================================================================================================================================
## R_fns_alg_paper_1_Mplus.R
##
## Original PS2 intercept-only six-test LC-MVP specification, using the unified Paper 1 data and case settings.
##
## ---- Prepare inspectable input/data files without starting Mplus -----------------------------------------------------------------------
##
fn_paper1_prepare_Mplus <-  function( y,
                                      N,
                                      MCMC_seed,
                                      settings
) {

        if (!requireNamespace(package = "MplusAutomation", quietly = TRUE)) stop("MplusAutomation is required for Mplus.")
        ##
        if (!is.matrix(x = y) || nrow(x = y) != N || ncol(x = y) != 6L || anyNA(x = y) || !all(y %in% c(0, 1))) {

            stop("Paper 1 Mplus requires the same N-by-six binary data as the other implementations.")

        }
        ##
        for (field in c("n_chains", "n_WCP", "n_threads", "n_fb_iter", "n_thin", "total_thread_limit")) {

            value <-  settings[[field]]
            if (length(x = value) != 1L || !is.finite(x = value) || value < 1 || value != floor(x = value)) {

                stop("Supply a positive integer Mplus setting: ", field)

            }

        }
        ##
        if (!is.character(settings$iteration_mode) || length(settings$iteration_mode) != 1L ||
            !settings$iteration_mode %in% c("FBITERATIONS", "BITERATIONS")) stop("Choose Mplus iteration_mode explicitly.")
        if (length(settings$biterations_minimum) != 1L || !is.finite(settings$biterations_minimum) ||
            settings$biterations_minimum < 0 || settings$biterations_minimum != floor(settings$biterations_minimum) ||
            (settings$iteration_mode == "BITERATIONS" && settings$biterations_minimum >= settings$n_fb_iter)) {
            stop("Mplus biterations_minimum must be a non-negative integer smaller than the requested maximum.")
        }
        if (length(settings$bconvergence) != 1L || !is.numeric(settings$bconvergence) ||
            !is.finite(settings$bconvergence) || settings$bconvergence < 0) {
            stop("Set Mplus bconvergence explicitly to a non-negative number.")
        }
        for (field in c("save_draws", "verify_saved_draws")) {
            if (!is.logical(settings[[field]]) || length(settings[[field]]) != 1L || is.na(settings[[field]])) stop("Set Mplus ", field, " explicitly.")
        }
        if (settings$verify_saved_draws && !settings$save_draws) stop("Saved-draw verification requires save_draws = TRUE.")
        ##
        ## ---- FBITERATIONS block-size validation:
        ##
        ## FBITERATIONS: Mplus 8.10 runs 100 x floor(FBITERATIONS / 100) iterations per chain, and at least 100: FBITERATIONS
        ## = 2, 4, 40, 80 and 150 all ran 100 (bparam.dat counts, identical draws for 2 / 4 / 40).
        ## BITERATIONS uses requested caps; CPU time alone does not establish actual iteration counts.
        ## Permit short BITERATIONS requests without imposing the FBITERATIONS restriction.
        ##
        if (identical(settings$iteration_mode, "FBITERATIONS") && settings$n_fb_iter %% 100 != 0) {

            stop("Mplus FBITERATIONS = ", settings$n_fb_iter, " is not a multiple of 100. ",
                 "The tested Mplus build rounds FBITERATIONS requests to blocks of 100. Use a multiple of 100 for FBITERATIONS.")

        }
        ##
        if (settings$n_threads != settings$n_chains * settings$n_WCP ||
            (isTRUE(x = settings$WCP) && settings$n_WCP <= 1) ||
            (!isTRUE(x = settings$WCP) && settings$n_WCP != 1)) {

            stop("Mplus requires n_threads = n_chains * n_WCP, with n_WCP = 1 for standard and greater than one for WCP.")

        }
        ##
        if (settings$n_threads > min(settings$total_thread_limit, parallel::detectCores())) {

            stop("Mplus processor request exceeds the configured or available thread limit.")

        }
        ##
        required_priors <-  c("prior_IW_nd", "prior_IW_d", "prior_prev_alpha", "prior_prev_beta",
                              "prior_beta_mean_test1_nd", "prior_beta_mean_test1_d",
                              "prior_beta_sd_test1_nd", "prior_beta_sd_test1_d")
        if (any(!required_priors %in% names(x = settings))) stop("Missing original PS2 Mplus prior settings.")
        ##
        if (is.null(x = settings$output_dir)) stop("Supply the Mplus case output directory.")
        dir.create(path = settings$output_dir, recursive = TRUE, showWarnings = FALSE)
        ##
        data <-  as.data.frame(x = y)
        names(x = data) <-  paste0("u", seq_len(length.out = 6L))
        ##
        iteration_control <-  if (settings$iteration_mode == "BITERATIONS") {
            paste0("BITERATIONS = ", settings$n_fb_iter, " (", settings$biterations_minimum, ");\n")
        } else paste0("FBITERATIONS = ", settings$n_fb_iter, ";\n")
        ##
        analysis_text <-  paste0("ESTIMATOR = BAYES;\n",
                                 "CHAINS = ", settings$n_chains, ";\n",
                                 "PROCESSORS = ", settings$n_threads, ";\n",
                                 "TYPE = MIXTURE;\n",
                                 iteration_control,
                                 "BCONVERGENCE = ", settings$bconvergence, ";\n",
                                 "THIN = ", settings$n_thin, ";\n",
                                 "STSEED = ", MCMC_seed, ";\n",
                                 "OPTSEED = ", MCMC_seed, ";\n",
                                 "MCSEED = ", MCMC_seed, ";\n",
                                 "BSEED = ", MCMC_seed, ";")
        ##
        model_text <-  paste0("%OVERALL%\n",
                              "[C#1*-1] (p43);\n",
                              "u1-u6 WITH u1-u6*0 (p1-p15);\n",
                              "[u1$1-u6$1*-1] (p16-p21);\n",
                              "%C#2%\n",
                              "u1-u6 WITH u1-u6*0 (p22-p36);\n",
                              "[u1$1-u6$1*+1] (p37-p42);")
        ##
        priors_text <-  paste0("p1-p15 ~ IW(0.0001, ", settings$prior_IW_d, ");\n",
                               "p22-p36 ~ IW(0.0001, ", settings$prior_IW_nd, ");\n",
                               "p43 ~ D(", settings$prior_prev_alpha, ", ", settings$prior_prev_beta, ");\n",
                               "p16 ~ N(", settings$prior_beta_mean_test1_d, ", ", settings$prior_beta_sd_test1_d^2, ");\n",
                               "p17-p21 ~ N(0, 1);\n",
                               "p37 ~ N(", settings$prior_beta_mean_test1_nd, ", ", settings$prior_beta_sd_test1_nd^2, ");\n",
                               "p38-p42 ~ N(0, 1);")
        ##
        model <-  MplusAutomation::mplusObject( TITLE = "Paper 1 PS2 LC-MVP parallel scaling;",
                                                VARIABLE = "NAMES = u1 u2 u3 u4 u5 u6; CATEGORICAL = u1-u6; CLASSES = C(2);",
                                                ANALYSIS = analysis_text,
                                                MODEL = model_text,
                                                MODELPRIORS = priors_text,
                                                SAVEDATA = if (settings$save_draws) "bparameters = bparam.dat;" else NULL,
                                                rdata = data,
                                                autov = FALSE,
                                                quiet = TRUE)
        ##
        case_directory <-  normalizePath(path = settings$output_dir, mustWork = TRUE)
        input_file <-  file.path(case_directory, "model.inp")
        data_file <-  file.path(case_directory, "data.dat")
        ## Keep paths inside the Mplus input short and relative to its case directory.
        previous_directory <-  getwd()
        on.exit(expr = setwd(dir = previous_directory), add = TRUE)
        setwd(dir = case_directory)
        ## run = 0 writes the input/data only. The command is not executed here.
        written_model <-  MplusAutomation::mplusModeler( object = model,
                                                        dataout = "data.dat",
                                                        modelout = "model.inp",
                                                        run = 0,
                                                        check = FALSE,
                                                        Mplus_command = "mpdemo",
                                                        writeData = "always",
                                                        hashfilename = FALSE,
                                                        quiet = TRUE)
        ##
        return(list(input_file = input_file, data_file = data_file, model = written_model))

}
##
## ---- Execute one prepared case; the shared Paper 1 executor measures the elapsed call -------------------------------------------------
##
fn_paper1_run_Mplus <-  function( y,
                                  N,
                                  MCMC_seed,
                                  settings
) {

        prepared <-  fn_paper1_prepare_Mplus(y = y, N = N, MCMC_seed = MCMC_seed, settings = settings)
        ##
        command <-  settings$Mplus_command
        if (is.null(x = command)) command <-  MplusAutomation::detectMplus()
        ##
        MplusAutomation::runModels( target = prepared$input_file,
                                    recursive = FALSE,
                                    showOutput = FALSE,
                                    replaceOutfile = "always",
                                    logFile = file.path(settings$output_dir, "Mplus_run.log"),
                                    Mplus_command = command,
                                    quiet = TRUE)
        ##
        output_file <-  sub(pattern = "\\.inp$", replacement = ".out", x = prepared$input_file)
        if (!file.exists(file = output_file)) stop("Mplus produced no output: ", output_file)
        output_text <-  readLines(con = output_file, warn = FALSE)
        ##
        ## This is a timing study: posterior convergence is not an acceptance condition.
        if (any(grepl(pattern = "*** ERROR", x = output_text, fixed = TRUE)) ||
            !any(grepl(pattern = "Elapsed Time", x = output_text, fixed = TRUE))) {
            stop("Mplus reported an input/runtime error or no completed timing record. Inspect: ", output_file)
        }
        ##
        return(list(input_file = prepared$input_file, output_file = output_file))

}
##
## ---- Verify the requested iteration control outside the case clock -------------------------------------------------------------------
##
## The output maximum is a requested cap, not an observed iteration counter. Ordinary timing runs do not require
## convergence or saved draws. Optional bparam.dat checking is available for diagnostics at THIN = 1.
## Informative normal threshold priors also trigger a PPPP computation; its diagnostic need not finish successfully
## in these short timing runs. Engine/input errors remain failures.
##
fn_paper1_verify_Mplus_run <-  function( output_file,
                                         n_fb_iter,
                                         n_chains,
                                         iteration_mode,
                                         verify_saved_draws
) {

        output_text <-  readLines(con = output_file, warn = FALSE)
        requested_label <-  if (iteration_mode == "FBITERATIONS") "Fixed number of iterations" else "Maximum number of iterations"
        reported_line <-  grep(requested_label, output_text, value = TRUE, fixed = TRUE)
        if (length(reported_line) != 1L || as.integer(sub(".*?([0-9]+)[[:space:]]*$", "\\1", reported_line)) != n_fb_iter) {
            stop("Mplus did not report the requested iteration control. Inspect: ", output_file)
        }
        ##
        ## ---- Mplus runs at least 200 FBITERATIONS iterations for a SINGLE chain:
        ##
        ## Measured on the HPC and the laptop: CHAINS = 1 with FBITERATIONS = 100 saves 200 draws, while 200 / 300 / 500 save
        ## 200 / 300 / 500 (Mplus splits the one chain into two halves for its convergence statistic and needs 100 per half);
        ## 2 and 4 chains save 100 per chain for FBITERATIONS = 100; BITERATIONS = 100 with 1 chain saves 100. So only the
        ## N = 50,000 long run (100 iterations) is affected. The expected count below is what the timing must use.
        ##
        n_iterations_expected <-  if (identical(iteration_mode, "FBITERATIONS") && n_chains == 1L) max(n_fb_iter, 200L) else n_fb_iter
        if (n_iterations_expected != n_fb_iter) {
            message(paste0("\033[36m", "    Mplus single-chain FBITERATIONS run: ", n_fb_iter, " requested, ",
                           n_iterations_expected, " iterations actually run (Mplus runs at least 200 for one chain); ",
                           "the timing uses the actual count.", "\033[0m"))
        }
        if (!verify_saved_draws) {
            return(list(iterations_per_chain = n_iterations_expected,
                        iterations_source    = "requested_count_and_single_chain_rule",
                        PPPP_run = any(grepl("Prior Posterior Predictive P-Value", output_text, fixed = TRUE))))
        }
        ## Optional diagnostic only; ordinary Paper 1 timing does not require saved draws.
        bparam_file <-  file.path(dirname(path = output_file), "bparam.dat")
        if (!file.exists(file = bparam_file)) stop("Mplus saved no bparam.dat to verify the iteration count: ", bparam_file)
        if (file.info(bparam_file)$size == 0) {
            stop("Mplus saved an empty bparam.dat; this sampling run's iteration count cannot be verified: ", bparam_file)
        }
        ##
        first_bparam_lines <-  readLines(con = bparam_file, n = 64, warn = FALSE)
        first_bparam_line <-  first_bparam_lines[nzchar(trimws(first_bparam_lines))][1]
        if (is.na(first_bparam_line)) stop("Mplus saved no draw records to verify: ", bparam_file)
        n_bparam_columns <-  length(x = scan(text = first_bparam_line, what = numeric(), quiet = TRUE))
        ## Single-chain BPARAMETERS files have only an iteration index; multiple-chain files prefix chain and iteration.
        n_index_columns <-  if (n_chains == 1L) 1L else 2L
        if (n_bparam_columns <= n_index_columns) stop("Mplus saved a malformed draw record: ", bparam_file)
        chain_and_iteration <-  utils::read.table( file       = bparam_file,
                                                   colClasses = c(rep("integer", n_index_columns), rep("NULL", n_bparam_columns - n_index_columns)))
        if (n_chains == 1L) chain_and_iteration <-  data.frame(chain = 1L, iteration = chain_and_iteration[[1]])
        names(x = chain_and_iteration) <-  c("chain", "iteration")
        iterations_per_chain <-  table(chain_and_iteration$chain)
        ##
        correct_indices <-  !anyNA(chain_and_iteration) &&
            identical(sort(unique(chain_and_iteration$chain)), seq_len(n_chains)) &&
            all(vapply(split(chain_and_iteration$iteration, chain_and_iteration$chain),
                       function(iterations) identical(sort(iterations), seq_len(n_iterations_expected)), logical(1)))
        if (!correct_indices) {

            stop("Mplus did not run the expected iterations: requested ", iteration_mode, " = ", n_fb_iter, " for ", n_chains,
                 " chains (expected ", n_iterations_expected, " per chain), bparam.dat has ", length(x = iterations_per_chain),
                 " chains with ", paste(unique(x = as.integer(x = iterations_per_chain)), collapse = "/"), " iterations. Inspect: ", bparam_file)

        }
        ##
        output_text <-  readLines(con = output_file, warn = FALSE)
        ##
        return(list( iterations_per_chain = n_iterations_expected,
                     iterations_source    = "saved_draws",
                     PPPP_run             = any(grepl(pattern = "Prior Posterior Predictive P-Value", x = output_text, fixed = TRUE))))

}






















