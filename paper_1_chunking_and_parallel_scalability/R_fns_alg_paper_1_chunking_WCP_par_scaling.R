##
## =======================================================================================================================================
## R_fns_alg_paper_1_chunking_WCP_par_scaling.R
##
## One computational study for Paper 1: chunk optimisation and parallel scaling use the SAME measurements.
## Functions used by alg_paper_1_chunking_WCP_par_scaling.R in this directory.
## In RStudio, open that runner, edit its settings, and click Source to run the study.
## This helper only defines functions; the runner controls execution and contains the editable study settings.
## The separate PS1 burn-in benchmark is outside this study.
##
## ---- Stan execution backend and explicit result provenance ----------------------------------------------------------------------------
##
fn_paper1_stan_backend <-  function( settings ) {

        backend <-  settings$stan_backend
        ##
        if (!is.character(x = backend) || length(x = backend) != 1L || is.na(x = backend) ||
            !backend %in% c("cmdstanr", "NicoStan")) stop("stan_backend must be 'cmdstanr' or 'NicoStan'.")
        ##
        return(backend)

}
##
## ---- Stan runtime receipt: one real TBB mapping is required outside every timed call -----------------------------------------------
##
fn_paper1_distinct_tbb_paths_from_maps <-  function( map_lines ) {

        path_start <-  regexpr(pattern = "[[:space:]]/", text = map_lines)
        mapped_paths <-  ifelse(test = path_start > 0L,
                                yes = substring(text = map_lines, first = path_start + 1L), no = "")
        deleted_paths <-  unique(x = sub(pattern = " \\(deleted\\)$", replacement = "", x =
                                         mapped_paths[grepl(pattern = " \\(deleted\\)$", x = mapped_paths)]))
        mapped_paths <-  sub(pattern = " \\(deleted\\)$", replacement = "", x = mapped_paths)
        mapped_paths <-  mapped_paths[nzchar(x = mapped_paths)]
        mapped_paths <-  mapped_paths[grepl(pattern = "^libtbb(\\.so|$)", x = basename(path = mapped_paths))]
        unavailable_paths <-  unique(x = mapped_paths[!file.exists(mapped_paths)])
        mapped_paths <-  unique(x = normalizePath(path = mapped_paths[file.exists(mapped_paths)], mustWork = TRUE))
        list(real_paths = mapped_paths,
             unavailable_paths = unavailable_paths,
             deleted_paths = deleted_paths)

}
##
fn_paper1_stan_runtime_receipt <-  function( stage,
                                             require_single_tbb = TRUE
) {

        system_name <-  unname(obj = Sys.info()[["sysname"]])
        receipt <-  list( receipt_version = 1L,
                          stage = as.character(x = stage),
                          system_name = system_name,
                          mapped_tbb_count = NA_integer_,
                          mapped_tbb_path = NA_character_,
                          mapped_tbb_md5 = NA_character_,
                          unavailable_tbb_path = NA_character_,
                          deleted_tbb_path = NA_character_,
                          verification = "not_applicable_non_Linux")
        ## /proc/self/maps records the libraries actually mapped by this R process, including dependencies loaded by a Stan evaluator.
        if (!identical(x = system_name, y = "Linux")) return(receipt)
        ##
        maps_file <-  "/proc/self/maps"
        if (!file.exists(maps_file)) stop("Linux Stan runtime receipt cannot read ", maps_file, ".")
        map_lines <-  readLines(con = maps_file, warn = FALSE)
        mapped_tbb_paths <-  fn_paper1_distinct_tbb_paths_from_maps(map_lines = map_lines)
        mapped_paths <-  mapped_tbb_paths$real_paths
        unavailable_paths <-  mapped_tbb_paths$unavailable_paths
        deleted_paths <-  mapped_tbb_paths$deleted_paths
        ##
        receipt$mapped_tbb_count <-  as.integer(x = length(x = mapped_paths) + length(x = unavailable_paths))
        receipt$unavailable_tbb_path <-  if (length(x = unavailable_paths)) paste(unavailable_paths, collapse = "; ") else NA_character_
        receipt$deleted_tbb_path <-  if (length(x = deleted_paths)) paste(deleted_paths, collapse = "; ") else NA_character_
        if (length(x = deleted_paths)) {

            receipt$verification <-  "deleted_tbb_mapping"
            if (isTRUE(x = require_single_tbb)) {

                stop("Stan/NicoStan timing found a deleted mapped libtbb path after ", stage,
                     ": ", paste(deleted_paths, collapse = "; "), ".")

            }

        } else if (length(x = unavailable_paths)) {

            receipt$verification <-  "invalid_tbb_mapping_path"
            if (isTRUE(x = require_single_tbb)) {

                stop("Stan/NicoStan timing found a mapped libtbb path that is no longer file-backed after ", stage,
                     ": ", paste(unavailable_paths, collapse = "; "), ".")

            }

        } else if (length(x = mapped_paths) == 1L) {

            receipt$mapped_tbb_path <-  mapped_paths
            receipt$mapped_tbb_md5 <-  unname(obj = tools::md5sum(files = mapped_paths))
            receipt$verification <-  "verified_single_tbb"

        } else {

            receipt$verification <-  "invalid_tbb_mapping_count"
            if (isTRUE(x = require_single_tbb)) {

                stop("Stan/NicoStan timing requires exactly one distinct real libtbb mapping after ", stage,
                     "; observed ", length(x = mapped_paths), ": ",
                     if (length(x = mapped_paths)) paste(mapped_paths, collapse = "; ") else "none", ".")

            }

        }
        ##
        return(receipt)

}
##
fn_paper1_stan_load_policy_fingerprint <-  function() {

        package_name <-  "NicoStan"
        package_namespace <-  asNamespace(ns = package_name)
        package_path <-  getNamespaceInfo(ns = package_namespace, which = "path")
        namespace_file <-  file.path(package_path, "NAMESPACE")
        function_md5 <-  function(function_name) {

                function_value <-  get0(x = function_name, envir = package_namespace,
                                         mode = "function", inherits = FALSE)
                if (is.null(x = function_value)) return(NA_character_)
                fn_paper1_object_md5(object = deparse(expr = function_value, control = "all"))

        }
        ## The installed NAMESPACE and the evaluated load hooks identify the R policy that selected the process runtime.
        return(list( package = package_name,
                     package_path = normalizePath(path = package_path, mustWork = TRUE),
                     NAMESPACE_md5 = if (file.exists(file = namespace_file)) unname(tools::md5sum(files = namespace_file)) else NA_character_,
                     onLoad_definition_md5 = function_md5(function_name = ".onLoad"),
                     setup_env_onload_definition_md5 = function_md5(function_name = ".setup_env_onload")))

}
##
fn_paper1_execution_columns <-  function( cases,
                                          stan_backend,
                                          native_backend
) {

        if (!"execution_backend" %in% names(x = cases)) {

            if (is.null(stan_backend) || is.null(native_backend)) stop("Missing execution_backend: supply explicit backend provenance.")

            cases$execution_backend <-  ifelse(test = grepl(pattern = "^AD_", x = cases$algorithm),
                                               yes = if (stan_backend == "NicoStan") "NicoStan_BridgeStan" else "cmdstanr",
                                               no = ifelse(test = grepl(pattern = "^MD_", x = cases$algorithm),
                                                           yes = native_backend, no = "mplus"))

        }
        ##
        return(cases)

}
##
fn_paper1_native_backend <-  function( settings ) {

        return (if (identical(x = settings$package_stack, y = "NicoStan_BayesMVP")) "BayesMVP_NicoStan" else "bayesmvp_native")

}
##
fn_paper1_load_split_packages <-  function( settings ) {

        if (!is.null(x = settings$package_library)) {

            library_path <-  normalizePath(path = settings$package_library, mustWork = TRUE)
            for (package in intersect(x = c("NicoStan", "BayesMVP"), y = loadedNamespaces())) {

                loaded_path <-  getNamespaceInfo(ns = asNamespace(ns = package), which = "path")
                if (normalizePath(path = dirname(path = loaded_path)) != library_path) {

                    stop(package, " is already loaded from another library. Restart R and Source the runner before loading it.")

                }

            }
            .libPaths(new = c(library_path, .libPaths()))

        }
        ## Load TBB first, then the shared sampler, then its specialised model provider.
        for (package in c("RcppParallel", "NicoStan", "BayesMVP")) {

            if (!requireNamespace(package = package, quietly = TRUE)) stop("Install the split package dependency: ", package)

        }
        for (package in c("NicoStan", "BayesMVP")) {

            if (utils::packageVersion(pkg = package) < numeric_version(x = "0.1.9000")) {

                stop("Paper 1 now requires the new NicoStan + BayesMVP packages; the loaded ", package, " is the older package.")

            }
            message(paste0("Using ", package, " ", utils::packageVersion(pkg = package), " from ", find.package(package = package)))

        }
        ##
        return(invisible(x = TRUE))

}
##
fn_paper1_stan_via_NicoStan_settings <-  function( settings ) {

        stan <-  settings$stan_via_NicoStan
        if (length(x = stan$step_size) != 1L || !is.finite(x = stan$step_size) || stan$step_size <= 0) {

            stop("Supply a positive Stan step_size for the NicoStan backend.")

        }
        ##
        if (length(stan$fixed_L) != 1L || !is.finite(stan$fixed_L) ||
            stan$fixed_L < 1 || stan$fixed_L != floor(stan$fixed_L)) {

            stop("Set stan_via_NicoStan$fixed_L explicitly to a positive integer for Stan-via-NicoStan.")

        }
        if (!identical(stan$randomize_tau, FALSE)) stop("Paper 1 requires stan_via_NicoStan$randomize_tau = FALSE for fixed L!")
        ##
        if (length(x = stan$init_radius) != 1L || !is.numeric(x = stan$init_radius) || !is.finite(x = stan$init_radius) || stan$init_radius < 0) {

            stop("Stan-via-NicoStan currently requires a non-negative numeric init radius.")

        }
        ##
        if (length(x = stan$metric_shape_main) != 1L || !stan$metric_shape_main %in% c("diag", "dense")) {

            stop("Unsupported Stan metric for the NicoStan backend.")

        }
        ##
        return(list( step_size = stan$step_size,
                     fixed_L = as.integer(x = stan$fixed_L),
                     randomize_tau = stan$randomize_tau,
                     metric_shape_main = stan$metric_shape_main,
                     init_radius = stan$init_radius,
                     initialisation = "R_uniform_unconstrained_conditioned_on_Stan_class_order_seeded_by_case_seed",
                     thread_policy = "shared_TBB_total_budget_chains_times_WCP",
                     make_args = list("STAN_THREADS=true", "PRECOMPILED_HEADERS=false")))

}
##
## ---- Check the trajectory length with the SAME compiled converter that the sampler call will use:
##
## The Stan arms sample through NicoStan's library and the BayesMVP arms through BayesMVP's own library (which
## compiles its own copy of the R-list converter), so package_name selects the library whose converter is checked.
## Returns the checked values so every result row can record the L and randomize_tau that reached C++.
##
fn_paper1_check_fixed_trajectory <-  function(  EHMC_args,
                                                fixed_L,
                                                package_name
) {

        check_function <-  get0(x = "Rcpp_fn_check_joint_trajectory", envir = asNamespace(ns = package_name),
                                 mode = "function", inherits = FALSE)
        if (is.null(x = check_function)) {

            stop("The installed ", package_name, " predates fixed-L support. Rebuild/install the updated ", package_name,
                 " package, then restart R before running this benchmark arm.")

        }
        ##
        checked <-  check_function(EHMC_args_as_Rcpp_List = EHMC_args)
        if (!identical(x = checked$randomize_tau, y = FALSE) ||
            !identical(as.numeric(x = checked$L_main), as.numeric(x = fixed_L))) {

            stop("The installed ", package_name, " C++ sampler did not accept the requested fixed trajectory length: requested L = ",
                 fixed_L, " with randomize_tau = FALSE, C++ received L = ", checked$L_main,
                 " with randomize_tau = ", checked$randomize_tau, ".")

        }
        ##
        return(invisible(x = checked))

}
##
fn_paper1_stan_HMC_args <-  function( configuration ) {

        EHMC_args <-  NicoStan:::init_EHMC_args_as_Rcpp_List(diffusion_HMC = FALSE)
        EHMC_args$eps_main <-  configuration$step_size
        EHMC_args$tau_main <-  configuration$fixed_L * configuration$step_size
        EHMC_args$tau_main_ii <-  EHMC_args$tau_main
        EHMC_args$randomize_tau <-  configuration$randomize_tau
        fn_paper1_check_fixed_trajectory(EHMC_args = EHMC_args, fixed_L = configuration$fixed_L, package_name = "NicoStan")
        ##
        return(EHMC_args)

}
##
## ---- Timing method: one timed run, or a short and a long run back to back:
##
## Every timed region also contains a fixed per-run cost S (state allocation, thread setup, per-chain Stan model
## loading, Mplus start-up, input writing and output parsing). With T(n) = S + n t, a SHORT run (n_low iterations)
## and a LONG run (n_high) with the same seed and configuration give t = (T_high - T_low) / (n_high - n_low),
## S = T_low - n_low t, and the sampling-only time of the long run, T_high - S = n_high t.
## "single_run" keeps the previous behaviour: one timed run at the long-run iteration count, S included.
## Provenance: introduced (a two-point linear fit of time on iterations, not a
## published benchmark method).
##
fn_paper1_algorithm_family <-  function( algorithm ) {

        algorithm_family <-  ifelse( test = grepl(pattern = "^AD_", x = algorithm),
                                     yes  = "stan",
                                     no   = ifelse( test = grepl(pattern = "^MD_", x = algorithm),
                                                    yes  = "bayesmvp",
                                                    no   = "mplus"))
        ##
        return(algorithm_family)

}
##
fn_paper1_iteration_setting_names <-  function() {

        return(list( bayesmvp = c( long_run  = "bayesmvp_iterations",
                                   short_run = "bayesmvp_iterations_short_run"),
                     stan     = c( long_run  = "stan_iterations",
                                   short_run = "stan_iterations_short_run"),
                     mplus    = c( long_run  = "mplus_WCP_iterations",
                                   short_run = "mplus_iterations_short_run")))

}
##
# fn_paper1_iterations_for_case <-  function( settings,
#                                             algorithm,
#                                             N,
#                                             which_run
# ) {
fn_paper1_iterations_for_case <-  function( settings,
                                            algorithm,
                                            N,
                                            which_run,
                                            n_chains = NA,
                                            num_chunks = NA,
                                            threads_per_chain = NA
) {

        algorithm_family <-  fn_paper1_algorithm_family(algorithm = algorithm)
        setting_name <-  fn_paper1_iteration_setting_names()[[algorithm_family]][[which_run]]
        ##
        ## ---- The BayesMVP WCP arm has its own LONG-run count: its runs use 4-16 chains x many threads, so at the
        ## standard arm's counts each long run lasted 0.2-0.4 s and repeats disagreed by 25-50% (OS scheduling jitter). The runner
        ## must set bayesmvp_wcp_iterations explicitly (about 10x the standard counts at N <= 10,000, 4x at N = 50,000, so that
        ## every WCP long run lasts at least ~2 s). The short run keeps bayesmvp_iterations_short_run.
        ##
        ## ---- The Stan WCP arm gets the SAME treatment: the runner sets the per-N stan_wcp_iterations separately
        ## from bayesmvp_wcp_iterations. Previously the Stan WCP arm wrongly used the standard stan_iterations.
        ##
        if (identical(algorithm, "AD_Stan_WCP") && identical(which_run, "long_run")) {

            setting_name <-  "stan_wcp_iterations"
            if (is.null(x = settings$stan_wcp_iterations)) stop("Set paper1_settings$stan_wcp_iterations in the runner (one long-run count per N for the Stan WCP arm).")

        }
        if (identical(algorithm, "MD_BayesMVP_WCP") && identical(which_run, "long_run")) {

            setting_name <-  "bayesmvp_wcp_iterations"
            if (is.null(x = settings$bayesmvp_wcp_iterations)) stop("Set paper1_settings$bayesmvp_wcp_iterations in the runner (one long-run count per N for the BayesMVP WCP arm).")

        }
        ##
        ## ---- The Mplus_standard arm has its own, much SMALLER long-run count, as in the original PS1 script
        ## (Mplus_standard FBITERATIONS = 5 x N_iter, Mplus_WCP = 30 x N_iter). Mplus_WCP reads mplus_WCP_iterations (the shared mplus_iterations Previously).
        ##
        if (identical(algorithm, "Mplus_standard") && identical(which_run, "long_run")) {

            setting_name <-  "mplus_standard_iterations"
            if (is.null(x = settings$mplus_standard_iterations)) stop("Set paper1_settings$mplus_standard_iterations in the runner (one long-run count per N for the Mplus_standard arm).")

        }
        n_iterations <-  unname(obj = settings[[setting_name]][as.character(x = N)])
        ##
        ## ---- 2026-10-03: the narrow-WCP allocations with many chains (n_chains >= wcp_many_chains_minimum;
        ##      up to 90 chains x 2 threads on the HPC) do 2-6x the work per iteration of the 4-16-chain WCP runs
        ##      that the WCP counts were sized for, so they read their own long-run counts; WCP-only cases
        ##      (N_chunks = N_threads/chain, the slowest, memory-bound ones) read
        ##      wcp_many_chains_iterations_WCP_only where set.
        ##      Callers that do not pass n_chains (e.g. the settings check) get the per-arm counts, as before.
        ##
        many_chain_counts <-  settings$wcp_many_chains_iterations[[algorithm]]
        if (identical(which_run, "long_run") && !is.null(x = many_chain_counts) &&
            length(x = n_chains) == 1 && !is.na(x = n_chains) && n_chains >= settings$wcp_many_chains_minimum) {

            setting_name <-  "wcp_many_chains_iterations"
            n_iterations <-  unname(obj = many_chain_counts[as.character(x = N)])
            WCP_only_counts <-  settings$wcp_many_chains_iterations_WCP_only[[algorithm]]
            if (!is.null(x = WCP_only_counts) && length(x = num_chunks) == 1 && !is.na(x = num_chunks) &&
                length(x = threads_per_chain) == 1 && !is.na(x = threads_per_chain) &&
                num_chunks == threads_per_chain) {

                setting_name <-  "wcp_many_chains_iterations_WCP_only"
                n_iterations <-  unname(obj = WCP_only_counts[as.character(x = N)])

            }

        }
        mplus_mode <-  if (identical(which_run, "short_run")) settings$mplus_iteration_mode_short_run else settings$mplus_iteration_mode
        ##
        if (length(x = n_iterations) != 1 || is.na(x = n_iterations) || !is.finite(x = n_iterations) ||
            n_iterations < 1 || n_iterations != floor(x = n_iterations)) {

            stop("Missing/invalid ", setting_name, " for ", algorithm, " at N = ", N, ": supply one positive integer per N.")

        }
        ##
        ## Mplus 8.10 runs 100 x floor(FBITERATIONS / 100) iterations per chain, and at least 100: FBITERATIONS = 2, 4, 40,
        ## 80 or 150 all ran 100 (bparam.dat counts and identical draws, measured), so other values never take effect.
        ## BITERATIONS uses requested caps; CPU time alone does not establish actual iteration counts.
        ## Permit short BITERATIONS requests without imposing the FBITERATIONS restriction.
        ##
        if (algorithm_family == "mplus" && "FBITERATIONS" %in% mplus_mode && n_iterations %% 100 != 0) {

            stop(setting_name, " = ", n_iterations, " for ", algorithm, " at N = ", N, " is not a multiple of 100. ",
                 "The tested Mplus build rounds FBITERATIONS requests to blocks of 100. Use a multiple of 100 for FBITERATIONS.")

        }
        if (algorithm_family == "mplus" && "BITERATIONS" %in% mplus_mode &&
            settings$mplus_biterations_minimum >= n_iterations) {
            stop("mplus_biterations_minimum must be smaller than ", setting_name, " for ", algorithm, " at N = ", N, ".")
        }
        ##
        return(n_iterations)

}
##
fn_paper1_timing_settings <-  function( settings,
                                        algorithms
) {

        timing_method <-  settings$timing_method
        ##
        if (is.null(x = timing_method) || !is.character(x = timing_method) || length(x = timing_method) != 1 ||
            is.na(x = timing_method) || !timing_method %in% c("two_run_difference", "single_run")) {

            stop("Set settings$timing_method to 'two_run_difference' or 'single_run'.")

        }
        ##
        tolerance_fraction <-  settings$two_run_negative_fixed_cost_tolerance_fraction
        ##
        if (timing_method == "two_run_difference" &&
            (is.null(x = tolerance_fraction) || length(x = tolerance_fraction) != 1 || !is.finite(x = tolerance_fraction) ||
             tolerance_fraction <= 0 || tolerance_fraction >= 1)) {

            stop("Set settings$two_run_negative_fixed_cost_tolerance_fraction to one number between 0 and 1.")

        }
        ##
        if (any(grepl("^Mplus_", algorithms))) {
            if (timing_method == "two_run_difference" &&
                (length(settings$mplus_short_run_role) != 1L || is.na(settings$mplus_short_run_role) ||
                 !settings$mplus_short_run_role %in% c("overhead_control", "iteration_pair"))) {
                stop("Set mplus_short_run_role to 'overhead_control' or 'iteration_pair' in the runner.")
            }
            ## One long-run mode, or both: c("FBITERATIONS", "BITERATIONS") runs every Mplus case in each mode;
            ## the first mode listed is the one used for the main analysis and manuscript tables.
            if (length(settings$mplus_iteration_mode) < 1L || anyNA(settings$mplus_iteration_mode) || anyDuplicated(settings$mplus_iteration_mode) ||
                !all(settings$mplus_iteration_mode %in% c("BITERATIONS", "FBITERATIONS"))) stop("Set mplus_iteration_mode explicitly in the runner.")
            if (timing_method == "two_run_difference") {
                if (!identical(settings$mplus_iteration_mode_short_run, "BITERATIONS")) {
                    stop("The Mplus short overhead run must use BITERATIONS; only the long-run mode is selectable.")
                }
                if (any(settings$mplus_iteration_mode != settings$mplus_iteration_mode_short_run) &&
                    !identical(settings$mplus_short_run_role, "overhead_control")) {
                    stop("Different Mplus short/long modes require mplus_short_run_role = 'overhead_control'.")
                }
            }
            if (length(settings$mplus_bconvergence) != 1L || !is.numeric(settings$mplus_bconvergence) ||
                !is.finite(settings$mplus_bconvergence) || settings$mplus_bconvergence < 0) {
                stop("Set mplus_bconvergence explicitly in the runner to a non-negative number.")
            }
            if (length(settings$mplus_biterations_minimum) != 1L || !is.finite(settings$mplus_biterations_minimum) ||
                settings$mplus_biterations_minimum < 0 || settings$mplus_biterations_minimum != floor(settings$mplus_biterations_minimum)) {
                stop("Set mplus_biterations_minimum explicitly in the runner to a non-negative integer.")
            }
            for (field in c("mplus_save_draws", "mplus_verify_saved_draws")) {
                if (length(settings[[field]]) != 1L || !is.logical(settings[[field]]) || is.na(settings[[field]])) stop("Set ", field, " explicitly in the runner.")
            }
            if (settings$mplus_verify_saved_draws && !settings$mplus_save_draws) stop("Saved-draw verification requires mplus_save_draws = TRUE.")
        }
        ##
        ## Check every selected arm and N before any sampling, so a mistyped or missing knob cannot pass silently.
        ##
        for (algorithm in unique(x = algorithms)) {

            selected_N <- settings$N_by_algorithm[[algorithm]]
            if (!length(selected_N) || any(!is.finite(selected_N)) || any(selected_N < 1) || any(selected_N != floor(selected_N))) {
                stop("Supply positive integer N_by_algorithm values for ", algorithm, " in the runner.")
            }
            for (N in unique(x = selected_N)) {

                n_iterations_long_run <-  fn_paper1_iterations_for_case( settings  = settings,
                                                                          algorithm = algorithm,
                                                                          N         = N,
                                                                          which_run = "long_run")
                ##
                if (timing_method == "two_run_difference") {

                    n_iterations_short_run <-  fn_paper1_iterations_for_case( settings  = settings,
                                                                               algorithm = algorithm,
                                                                               N         = N,
                                                                               which_run = "short_run")
                    ##
                    if (n_iterations_short_run >= n_iterations_long_run) {

                        stop("The short run must use fewer iterations than the long run for ", algorithm, " at N = ", N,
                             " (short = ", n_iterations_short_run, ", long = ", n_iterations_long_run, ").")

                    }

                }

            }

        }
        ##
        ##
        ## ---- Untimed warm-up before each timed pair: the runner must set it explicitly.
        ##
        warm_up <-  settings$untimed_warm_up_run_before_timing
        if (is.null(x = warm_up) || !is.logical(x = warm_up) || length(x = warm_up) != 1 || is.na(x = warm_up)) {

            stop("Set settings$untimed_warm_up_run_before_timing to TRUE or FALSE in the runner.")

        }
        message(paste0("\033[36m", "Untimed 1-iteration warm-up call before every timed pair: ", warm_up, "\033[0m"))
        ##
        return(list( timing_method      = timing_method,
                     tolerance_fraction = if (timing_method == "two_run_difference") tolerance_fraction else NA_real_,
                     untimed_warm_up_run_before_timing = warm_up))

}
##
## ---- Fields every case runner may return; custom test runners may omit the optional ones:
##
fn_paper1_runner_value_or_missing <-  function( runner_output,
                                                field_name,
                                                missing_value
) {

        field_value <-  runner_output[[field_name]]
        ##
        return(if (is.null(x = field_value)) missing_value else field_value)

}
##
## ---- Time one configuration with the selected timing method:
##
fn_paper1_time_case_with_timing_method <-  function( case,
                                                     y,
                                                     settings,
                                                     runtime,
                                                     case_runner,
                                                     timing_settings
) {

        fn_check_raw_elapsed_time <-  function( runner_output,
                                                run_label
        ) {

                if (length(x = runner_output$elapsed_seconds) != 1L || !is.finite(x = runner_output$elapsed_seconds) ||
                    runner_output$elapsed_seconds <= 0) stop("The case runner returned an invalid elapsed time for the ", run_label, ".")

        }
        ##
        if (timing_settings$timing_method == "single_run") {

            single_run_case <-  case
            single_run_case$timing_run_label <-  "single_run"
            ##
            single_run_output <-  case_runner(case = single_run_case, y = y, settings = settings, runtime = runtime)
            fn_check_raw_elapsed_time(runner_output = single_run_output, run_label = "single run")
            ##
            single_run_output$elapsed_seconds_short_run <-  NA_real_
            single_run_output$elapsed_seconds_long_run <-  single_run_output$elapsed_seconds
            single_run_output$fixed_cost_seconds <-  NA_real_
            single_run_output$seconds_per_iteration <-  NA_real_
            single_run_output$two_run_timing_flag <-  NA_character_
            single_run_output$timing_estimator <-  "single_run"
            single_run_output$divergences_short_run <-  NA_real_
            single_run_output$mplus_iterations_per_chain_short_run <-  NA_real_
            ##
            return(single_run_output)

        }
        ##
        n_iterations_short_run <-  case$n_iter_short_run
        n_iterations_long_run <-  case$n_iter
        ##
        short_run_case <-  case
        short_run_case$n_iter <-  n_iterations_short_run
        short_run_case$timing_run_label <-  paste0("short_run_", n_iterations_short_run, "_iterations")
        ##
        long_run_case <-  case
        long_run_case$timing_run_label <-  paste0("long_run_", n_iterations_long_run, "_iterations")
        ##
        ## ---- Untimed warm-up call: same configuration, 1 iteration, result discarded, so that the short
        ## run is not the cold first call of the pair (start-up at N = 50,000 with 64+ chains is 4-11 s and colder the
        ## first time; without this T_long < T_short happened). Mplus runs it in its short-run mode (see is_short_call).
        ##
        if (isTRUE(x = timing_settings$untimed_warm_up_run_before_timing)) {

            warm_up_case <-  case
            warm_up_case$n_iter <-  1
            warm_up_case$timing_run_label <-  "warm_up_1_iterations_untimed"
            warm_up_output <-  case_runner(case = warm_up_case, y = y, settings = settings, runtime = runtime)
            message(paste0("\033[36m", "    warm-up call (1 iteration, untimed): ",
                           signif(x = warm_up_output$elapsed_seconds, digits = 4), " s", "\033[0m"))

        }
        ##
        ## Short run first, then the long run, back to back with the same seed and configuration.
        ##
        short_run_output <-  case_runner(case = short_run_case, y = y, settings = settings, runtime = runtime)
        fn_check_raw_elapsed_time(runner_output = short_run_output, run_label = "short run")
        ##
        ## ---- Print each run's time as soon as it finishes, so progress is visible while the study runs:
        ##
        message(paste0("\033[36m", "    short run (", n_iterations_short_run, " iterations): ",
                       signif(x = short_run_output$elapsed_seconds, digits = 4), " s", "\033[0m"))
        ##
        long_run_output <-  case_runner(case = long_run_case, y = y, settings = settings, runtime = runtime)
        fn_check_raw_elapsed_time(runner_output = long_run_output, run_label = "long run")
        ##
        message(paste0("\033[36m", "    long run  (", n_iterations_long_run, " iterations): ",
                       signif(x = long_run_output$elapsed_seconds, digits = 4), " s", "\033[0m"))
        ##
        ## Both runs must have used the same trajectory settings, or their difference is not n t.
        ##
        if (!identical(x = fn_paper1_runner_value_or_missing(short_run_output, "mean_L", NA_real_),
                       y = fn_paper1_runner_value_or_missing(long_run_output, "mean_L", NA_real_)) ||
            !identical(x = fn_paper1_runner_value_or_missing(short_run_output, "randomize_tau", NA),
                       y = fn_paper1_runner_value_or_missing(long_run_output, "randomize_tau", NA))) {

            stop("The short and long runs of case ", case$case_id, " used different trajectory settings.")

        }
        ##
        elapsed_seconds_short_run <-  short_run_output$elapsed_seconds
        elapsed_seconds_long_run <-  long_run_output$elapsed_seconds
        ##
        overhead_control <-  grepl("^Mplus_", case$algorithm) && identical(settings$mplus_short_run_role, "overhead_control")
        if (overhead_control) {
            ## The complete short call estimates overhead. Its actual sampling count is not assumed to equal its cap.
            fixed_cost_seconds <-  elapsed_seconds_short_run
            ##
            ## ---- Use the iterations Mplus actually ran: a single chain under FBITERATIONS runs 2 x the request
            ## (see fn_paper1_verify_Mplus_run), so the per-iteration time divides by the actual count and the reported
            ## time is that per-iteration time for the REQUESTED count, keeping every chain count on the same iteration budget.
            ##
            n_iterations_long_run_actual <-  fn_paper1_runner_value_or_missing(long_run_output, "mplus_iterations_per_chain", NA_real_)
            if (!is.finite(n_iterations_long_run_actual) || n_iterations_long_run_actual < 1) n_iterations_long_run_actual <-  n_iterations_long_run
            if (n_iterations_long_run_actual != n_iterations_long_run) {
                message(paste0("\033[36m", "    long run actually ran ", n_iterations_long_run_actual, " iterations per chain (requested ",
                               n_iterations_long_run, "); reported time is scaled to the requested count.", "\033[0m"))
            }
            seconds_per_iteration <-  (elapsed_seconds_long_run - elapsed_seconds_short_run) / n_iterations_long_run_actual
            reported_seconds <-  seconds_per_iteration * n_iterations_long_run
            timing_estimator <-  "overhead_subtraction"
        } else {
            seconds_per_iteration <-  (elapsed_seconds_long_run - elapsed_seconds_short_run) / (n_iterations_long_run - n_iterations_short_run)
            fixed_cost_seconds <-  elapsed_seconds_short_run - n_iterations_short_run * seconds_per_iteration
            reported_seconds <-  elapsed_seconds_long_run - fixed_cost_seconds
            timing_estimator <-  "iteration_difference"
        }
        ##
        ## Keep every number; flagged rows are reported, never dropped.
        ##
        two_run_timing_flag <-  if (elapsed_seconds_long_run <= elapsed_seconds_short_run) {

            "long_run_not_slower_than_short_run"

        } else if (fixed_cost_seconds < -timing_settings$tolerance_fraction * elapsed_seconds_long_run) {
            ##
            ## Judged against the long-run time: a negative S adds |S| to the reported time T_long - S, so this bounds
            ## the size of that error. A short run of one iteration is too short a yardstick (milliseconds of jitter).
            ##

            "fixed_cost_negative_beyond_tolerance"

        } else "ok"
        ##
        message(paste0("\033[36m", if (overhead_control) "    per requested long-run iteration = " else "    per iteration = ",
                       signif(x = seconds_per_iteration, digits = 4), " s",
                       if (overhead_control) "; overhead control = " else "; fixed cost = ", signif(x = fixed_cost_seconds, digits = 4), " s",
                       "; reported time (long run - fixed cost) = ", signif(x = reported_seconds, digits = 4), " s",
                       "\033[0m"))
        ##
        if (two_run_timing_flag != "ok") {

            message(paste0("\033[36m", "Two-run timing flag for case ", case$case_id, ": ", two_run_timing_flag,
                           " (T_short = ", signif(x = elapsed_seconds_short_run, digits = 4), " s at ", n_iterations_short_run,
                           " iterations, T_long = ", signif(x = elapsed_seconds_long_run, digits = 4), " s at ", n_iterations_long_run,
                           " iterations). The row is kept and flagged.", "\033[0m"))

        }
        ##
        timed_output <-  long_run_output
        timed_output$elapsed_seconds <-  reported_seconds
        timed_output$elapsed_seconds_short_run <-  elapsed_seconds_short_run
        timed_output$elapsed_seconds_long_run <-  elapsed_seconds_long_run
        timed_output$fixed_cost_seconds <-  fixed_cost_seconds
        timed_output$seconds_per_iteration <-  seconds_per_iteration
        timed_output$two_run_timing_flag <-  two_run_timing_flag
        timed_output$timing_estimator <-  timing_estimator
        timed_output$divergences_short_run <-  short_run_output$divergences
        timed_output$mplus_iterations_per_chain_short_run <-  fn_paper1_runner_value_or_missing( runner_output = short_run_output,
                                                                                                 field_name    = "mplus_iterations_per_chain",
                                                                                                 missing_value = NA_real_)
        ##
        return(timed_output)

}
##
## ---- Run one complete Stan timing case in a fresh R process ---------------------------------------------------------------------------
##
## The untimed warm-up, short run and long run stay together in one child so the paired timing contract is unchanged.
## Only plain prepared-model cache entries return to the parent; BridgeStan pointers and native traces never cross the process boundary.
## This releases retained native memory between cases; it does not reduce the live-memory peak inside the largest case.
##
fn_paper1_should_isolate_Stan_case <-  function( case,
                                                settings,
                                                case_runner
) {

        return(isTRUE(x = settings$run_each_Stan_case_in_fresh_R_process) &&
               identical(x = fn_paper1_stan_backend(settings = settings), y = "NicoStan") &&
               grepl(pattern = "^AD_", x = case$algorithm) &&
               identical(x = case_runner, y = fn_paper1_run_case))

}
##
fn_paper1_time_Stan_case_in_fresh_R_process <-  function( case,
                                                          y,
                                                          settings,
                                                          runtime,
                                                          timing_settings
) {

        if (!requireNamespace(package = "callr", quietly = TRUE)) {

            stop("callr is required when run_each_Stan_case_in_fresh_R_process = TRUE.")

        }
        if (!is.environment(x = runtime$stan_via_NicoStan_cache)) {

            stop("Fresh Stan case execution requires the execution-local prepared-model cache.")

        }
        ## Pass only the selected plain model/data objects and settings used by the NicoStan Stan case runner.
        worker_settings <-  list( algorithm_study_dir = settings$algorithm_study_dir,
                                  package_stack = settings$package_stack,
                                  package_library = settings$package_library,
                                  total_thread_limit = settings$total_thread_limit,
                                  stan_backend = settings$stan_backend,
                                  stan_via_NicoStan = settings$stan_via_NicoStan,
                                  output_dir = settings$output_dir)
        model_specification <-  runtime$stan_models[[as.character(x = case$algorithm)]]
        stan_data <-  runtime$stan_data[[as.character(x = case$N)]]
        if (!is.list(x = model_specification) || !is.list(x = stan_data)) {

            stop("Fresh Stan case execution could not resolve the selected model specification or Stan data.")

        }
        prepared_model_cache <-  as.list(x = runtime$stan_via_NicoStan_cache, all.names = TRUE)
        helper_file <-  file.path(settings$algorithm_study_dir,
                                  "paper_1_chunking_and_parallel_scalability",
                                  "R_fns_alg_paper_1_chunking_WCP_par_scaling.R")
        if (!file.exists(file = helper_file)) stop("Missing Paper 1 helper for the fresh Stan case process: ", helper_file)
        parent_pid <-  Sys.getpid()
        parent_libpath <-  .libPaths()
        parent_working_dir <-  getwd()
        ## Self-contained worker: source the same helper, rebuild the small cache environment, then run the complete timing pair.
        worker <-  function(case,
                            y,
                            settings,
                            model_specification,
                            stan_data,
                            prepared_model_cache,
                            timing_settings,
                            helper_file,
                            libpath,
                            working_dir,
                            parent_pid) {

                setwd(dir = working_dir)
                .libPaths(new = libpath)
                worker_environment <-  new.env(parent = globalenv())
                sys.source(file = helper_file, envir = worker_environment)
                worker_environment$fn_paper1_load_split_packages(settings = settings)
                worker_runtime <-  list( stan_models = stats::setNames(object = list(model_specification), nm = as.character(x = case$algorithm)),
                                         stan_data = stats::setNames(object = list(stan_data), nm = as.character(x = case$N)),
                                         stan_via_NicoStan_cache = list2env(x = prepared_model_cache,
                                                                           envir = new.env(parent = emptyenv())))
                ## A prepared cache entry skips BridgeStan's R6 initialisation. In that path only, pin the generated model
                ## library through R's loader on the controlling thread so native worker dlclose calls cannot unload it early.
                stan_partition <-  worker_environment$fn_paper1_stan_partition_columns(cases = case)
                chunk_size <-  if (case$algorithm == "AD_Stan") "none" else as.character(x = stan_partition$stan_chunk_size)
                model_key <-  paste0(model_specification$stan_source_md5, "_N_", case$N, "_chunk_size_", chunk_size)
                prepared_model <-  prepared_model_cache[[model_key]]
                if (!is.null(x = prepared_model)) {

                    model_library_path <-  prepared_model$model_args$model_so_file
                    if (is.null(x = model_library_path)) model_library_path <-  model_specification$library_path
                    if (!is.character(x = model_library_path) || length(x = model_library_path) != 1L ||
                        is.na(x = model_library_path) || !file.exists(file = model_library_path)) {

                        stop("Fresh Stan case process could not resolve the compiled model library to pin.")

                    }
                    assign(x = ".paper1_pinned_Stan_model_DLL",
                           value = dyn.load(x = normalizePath(path = model_library_path, mustWork = TRUE),
                                            local = TRUE, now = TRUE),
                           envir = .GlobalEnv)

                }
                output <-  worker_environment$fn_paper1_time_case_with_timing_method(
                                 case = case,
                                 y = y,
                                 settings = settings,
                                 runtime = worker_runtime,
                                 case_runner = worker_environment$fn_paper1_run_case,
                                 timing_settings = timing_settings)
                output$execution_process_mode <-  "fresh_R_process_per_complete_Stan_case"
                output$execution_process_pid <-  Sys.getpid()
                output$execution_process_parent_pid <-  parent_pid
                list( output = output,
                      prepared_model_cache = as.list(x = worker_runtime$stan_via_NicoStan_cache, all.names = TRUE))

        }
        environment(worker) <-  baseenv()
        message(paste0("\033[36mStarting fresh R process for complete Stan case ", case$case_id,
                       " (warm-up, short run and long run).\033[0m"))
        child_process <-  NULL
        on.exit({

            if (!is.null(x = child_process) && isTRUE(x = child_process$is_alive())) {

                try(expr = child_process$kill_tree(), silent = TRUE)

            }

        }, add = TRUE)
        child_process <-  tryCatch({

            callr::r_bg( func = worker,
                         args = list(case = case,
                                     y = y,
                                     settings = worker_settings,
                                     model_specification = model_specification,
                                     stan_data = stan_data,
                                     prepared_model_cache = prepared_model_cache,
                                     timing_settings = timing_settings,
                                     helper_file = helper_file,
                                     libpath = parent_libpath,
                                     working_dir = parent_working_dir,
                                     parent_pid = parent_pid),
                         libpath = parent_libpath,
                         wd = parent_working_dir,
                         stdout = "",
                         stderr = "",
                         system_profile = FALSE,
                         user_profile = FALSE,
                         supervise = TRUE,
                         error = "error")

        }, error = function(condition) {

            stop(paste0("Could not start fresh Stan case process for case ", case$case_id, ": ", conditionMessage(condition)),
                 call. = FALSE)

        })
        tryCatch({

            child_process$wait()

        }, interrupt = function(condition) {

            if (isTRUE(x = child_process$is_alive())) try(expr = child_process$kill_tree(), silent = TRUE)
            stop(paste0("Fresh Stan case process was interrupted for case ", case$case_id, "."), call. = FALSE)

        })
        child_exit_status <-  child_process$get_exit_status()
        if (length(x = child_exit_status) != 1L || is.na(x = child_exit_status) || child_exit_status != 0L) {

            stop(paste0("Fresh Stan case process for case ", case$case_id,
                        " did not exit cleanly (exit status ",
                        if (length(x = child_exit_status) == 1L && !is.na(x = child_exit_status)) child_exit_status else "unavailable",
                        "). Its returned timing result and prepared-model cache were rejected."),
                 call. = FALSE)

        }
        child_result <-  tryCatch({

            child_process$get_result()

        }, error = function(condition) {

            stop(paste0("Fresh Stan case process failed for case ", case$case_id, ": ", conditionMessage(condition)),
                 call. = FALSE)

        })
        if (!is.list(x = child_result) || !is.list(x = child_result$output) ||
            !is.list(x = child_result$prepared_model_cache)) {

            stop("Fresh Stan case process returned an incomplete result.")

        }
        ## Preserve the expensive prepared-model cache across children. Its entries contain only model arguments, dimensions and paths.
        for (cache_name in names(x = child_result$prepared_model_cache)) {

            assign(x = cache_name, value = child_result$prepared_model_cache[[cache_name]],
                   envir = runtime$stan_via_NicoStan_cache)

        }
        message(paste0("\033[36mFresh R process for Stan case ", case$case_id,
                       " exited; its native worker memory has been released.\033[0m"))
        return(child_result$output)

}
##
## ---- One COVID-19 simulation shared by every implementation ---------------------------------------------------------------------------
##
fn_paper1_object_md5 <-  function( object ) {

        temporary_file <-  tempfile(pattern = "paper1_signature_")
        ##
        on.exit(expr = unlink(x = temporary_file), add = TRUE)
        ##
        saveRDS(object = object, file = temporary_file, compress = FALSE, version = 2)
        ##
        return(unname(obj = tools::md5sum(files = temporary_file)))

}
##
fn_paper1_COVID_data <-  function( algorithm_study_dir,
                                   N_vec,
                                   seed
) {

        simulator_file <-  file.path( algorithm_study_dir,
                                      "0_utilities",
                                      "shared_functions",
                                      "R_fn_sim_bin_COVID_19_LC_MVP_data.R")
        ##
        simulator_environment <-  new.env(parent = globalenv())
        ##
        sys.source(file = simulator_file, envir = simulator_environment)
        ## Generate the full original N grid BEFORE selecting subsets, preserving the seed-123 datasets.
        simulated_data <-  simulator_environment$R_fn_simulate_binary_LC_MVP_data_COVID_19( study_type = "algorithm_parallel_scaling_tests",
                                                                                            seed = seed,
                                                                                            corr_force_positive_DGM = FALSE)
        ##
        dataset_positions <-  match(x = N_vec, table = simulated_data$N_vec)
        ##
        if (anyNA(x = dataset_positions) || anyDuplicated(x = N_vec)) stop("Requested N values must be unique COVID simulator sizes.")
        ## Keep all per-dataset fields aligned, including populations and latent truth used by downstream readers.
        dataset_fields <-  c("y_list", "X_list", "pop_list", "Sigma_nd_true_observed_list", "Sigma_d_true_observed_list",
                             "prev_true_observed_list", "Se_true_observed_list", "Sp_true_observed_list")
        ##
        for (field in dataset_fields) simulated_data[[field]] <-  simulated_data[[field]][dataset_positions]
        ##
        simulated_data$N_vec <-  N_vec
        ##
        simulated_data$y_binary_list <-  simulated_data$y_list
        ##
        for (dataset_index in seq_along(along.with = N_vec)) {

            y <-  simulated_data$y_list[[dataset_index]]
            ##
            if (nrow(x = y) != N_vec[dataset_index] || ncol(x = y) != 6L || anyNA(x = y) ||
                !all(y %in% c(0, 1))) stop("COVID benchmark data must have N rows and six binary outcomes.")

        }
        ##
        simulated_data$benchmark_data <-  list( generator = "R_fn_simulate_binary_LC_MVP_data_COVID_19",
                                                data_seed = seed,
                                                corr_force_positive_DGM = FALSE,
                                                simulator_md5 = unname(obj = tools::md5sum(files = simulator_file)),
                                                dataset_md5 = setNames( object = vapply( X = simulated_data$y_list,
                                                                                         FUN = fn_paper1_object_md5,
                                                                                         FUN.VALUE = character(length = 1)),
                                                                        nm = as.character(x = N_vec)))
        ##
        return(simulated_data)

}

## ---- A single experiment grid: no PS1/PS2 identity and no repeated matching cells -------------------------------------------------------
##
fn_paper1_case_identity <-  function( grid ) {

        grid <-  fn_paper1_execution_columns(cases = grid, stan_backend = NULL, native_backend = NULL)
        identity_columns <-  c("device", "algorithm", "execution_backend", "N", "num_chunks", "n_chains", "threads_per_chain",
                               "n_iter", "run", "seed")
        if ("mplus_iteration_mode" %in% names(x = grid)) identity_columns <-  c(identity_columns, "mplus_iteration_mode")
        return(grid[, identity_columns, drop = FALSE])

}
##
fn_paper1_stan_partition_columns <-  function( cases ) {

        ## Drop the old exported column when reading historical results; all Stan arms now use one name.
        cases$stan_grainsize <-  NULL
        cases$stan_chunk_size <-  NA_integer_
        ##
        stan_chunked_cases <-  cases$algorithm %in% c("AD_Stan_chunked", "AD_Stan_tape_chunked", "AD_Stan_WCP")
        cases$stan_chunk_size[stan_chunked_cases] <-  as.integer(x = ceiling(x = cases$N[stan_chunked_cases] /
                                                                                cases$num_chunks[stan_chunked_cases]))
        ##
        return(cases)

}
##
fn_paper1_benchmark_grid <-  function( settings ) {

        stan_backend <-  fn_paper1_stan_backend(settings = settings)
        ##
        supported <-  c("MD_BayesMVP", "MD_BayesMVP_multi_process", "MD_BayesMVP_WCP",
                        "AD_Stan", "AD_Stan_chunked", "AD_Stan_tape_chunked", "AD_Stan_WCP",
                        "Mplus_standard", "Mplus_WCP")
        ##
        if (!length(x = settings$algorithms) || any(!settings$algorithms %in% supported)) stop("Unsupported or empty algorithm selection.")
        ##
        timing_settings <-  fn_paper1_timing_settings(settings = settings, algorithms = settings$algorithms)
        ##
        if (length(x = settings$n_runs) != 1L || !is.finite(x = settings$n_runs) ||
            settings$n_runs < 1 || settings$n_runs != floor(x = settings$n_runs)) stop("n_runs must be a positive integer.")
        ##
        if (length(x = settings$total_thread_limit) != 1L || !is.finite(x = settings$total_thread_limit) ||
            settings$total_thread_limit < 1 || settings$total_thread_limit != floor(x = settings$total_thread_limit)) {

            stop("total_thread_limit must be a positive integer.")

        }
        ##
        rows <-  list()
        skipped_cases <-  0
        ##
        for (algorithm in unique(x = settings$algorithms)) {

            algorithm_n_runs <-  if (grepl(pattern = "^Mplus_", x = algorithm) ) settings$mplus_n_runs else settings$n_runs
            if (length(algorithm_n_runs) != 1L || !is.finite(algorithm_n_runs) || algorithm_n_runs < 1 ||
                algorithm_n_runs != floor(algorithm_n_runs)) stop("Invalid repeat count for ", algorithm, ".")
            ##
            chain_counts <-  settings$n_chains_by_algorithm[[algorithm]]
            sample_sizes <-  settings$N_by_algorithm[[algorithm]]
            wcp_algorithm <-  algorithm %in% c("MD_BayesMVP_WCP", "AD_Stan_WCP", "Mplus_WCP")
            chunked_wcp_algorithm <-  algorithm %in% c("MD_BayesMVP_WCP", "AD_Stan_WCP")
            ##
            if (!is.numeric(x = chain_counts) || !length(x = chain_counts) || anyNA(x = chain_counts) ||
                any(!is.finite(x = chain_counts)) || any(chain_counts < 1 | chain_counts != floor(x = chain_counts)) ||
                anyDuplicated(x = chain_counts)) stop("Supply unique positive integer chain counts for ", algorithm)
            ##
            if (!length(x = sample_sizes) || any(!sample_sizes %in% settings$N_vec)) stop("Missing/unsupported N selection for ", algorithm)
            ##
            for (N in unique(x = sample_sizes)) {

                chunks <-  settings$chunks_by_algorithm[[algorithm]][[as.character(x = N)]]
                if (algorithm %in% c("AD_Stan", "Mplus_standard", "Mplus_WCP")) chunks <-  1
                ##
                if (!is.numeric(x = chunks) || !length(x = chunks) || anyNA(x = chunks) ||
                    any(!is.finite(x = chunks)) || any(chunks < 1 | chunks != floor(x = chunks))) {

                    stop("Missing/invalid chunk grid for ", algorithm, " at N = ", N)

                }
                ##
                if (wcp_algorithm) {

                    wcp_candidates_by_chain <-  settings$threads_per_chain_by_algorithm[[algorithm]][[as.character(x = N)]]
                    ##
                    if (!is.list(x = wcp_candidates_by_chain) || is.null(x = names(x = wcp_candidates_by_chain)) ||
                        anyNA(x = names(x = wcp_candidates_by_chain)) || any(names(x = wcp_candidates_by_chain) == "") ||
                        anyDuplicated(x = names(x = wcp_candidates_by_chain)) ||
                        !all(as.character(x = chain_counts) %in% names(x = wcp_candidates_by_chain))) {

                        stop("Supply WCP candidate lists keyed by each selected chain count for ", algorithm, " at N = ", N)

                    }

                }
                ##
                run_order <-  if (grepl(pattern = "^MD_", x = algorithm)) {

                    expand.grid(chain_index = seq_along(along.with = chain_counts), run = seq_len(length.out = algorithm_n_runs))

                } else {

                    expand.grid(run = seq_len(length.out = algorithm_n_runs), chain_index = seq_along(along.with = chain_counts))

                }
                ##
                rows_before_N <-  length(x = rows)
                parallel_case_found <-  !wcp_algorithm
                ##
                ## ---- Extra chunk counts run with ONE chain only (settings$serial_reference_extra_chunks, keyed by N), so that every
                ##      WCP-only configuration (N_chunks = N_threads_per_chain) has a one-chain run at the same chunk count for the
                ##      matched-serial efficiency. Counts already in the chunk grid are not duplicated.
                serial_reference_extra_chunks <-  if (algorithm %in% c("MD_BayesMVP", "AD_Stan_tape_chunked")) {

                    setdiff(x = settings$serial_reference_extra_chunks[[as.character(x = N)]], y = chunks)

                } else NULL
                ##
                if (length(x = serial_reference_extra_chunks) && (!is.numeric(x = serial_reference_extra_chunks) ||
                    anyNA(x = serial_reference_extra_chunks) || any(serial_reference_extra_chunks < 1) ||
                    any(serial_reference_extra_chunks != floor(x = serial_reference_extra_chunks)))) {

                    stop("Invalid serial_reference_extra_chunks for ", algorithm, " at N = ", N)

                }
                ##
                for (num_chunks in unique(x = c(chunks, serial_reference_extra_chunks))) {

                    for (run_index in seq_len(length.out = nrow(x = run_order))) {

                        chain_index <-  run_order$chain_index[run_index]
                        repeat_index <-  run_order$run[run_index]
                        n_chains <-  chain_counts[chain_index]
                        ##
                        if (num_chunks %in% serial_reference_extra_chunks && n_chains != 1) next
                        ##
                        wcp_candidates <-  if (wcp_algorithm) {

                            wcp_candidates_by_chain[[as.character(x = n_chains)]]

                        } else 1
                        ##
                        if (!is.numeric(x = wcp_candidates) || !length(x = wcp_candidates) || anyNA(x = wcp_candidates) ||
                            any(!is.finite(x = wcp_candidates)) ||
                            any(wcp_candidates < 1 | wcp_candidates != floor(x = wcp_candidates))) {

                            stop("Supply positive integer WCP candidates for ", algorithm, " at N = ", N, ", n_chains = ", n_chains)

                        }
                        ## One thread per chain belongs to the existing non-WCP arms.
                        if (wcp_algorithm && any(wcp_candidates == 1)) {

                            stop("WCP candidates must be greater than one; use the non-WCP arm for WCP=1: ",
                                 algorithm, " at N = ", N, ", n_chains = ", n_chains)

                        }
                        ##
                        feasible <-  n_chains * wcp_candidates <= settings$total_thread_limit
                        ##
                        if (chunked_wcp_algorithm) feasible <-  feasible & wcp_candidates <= num_chunks
                        ##
                        skipped_cases <-  skipped_cases + sum(!feasible)
                        ##
                        for (threads_per_chain in unique(x = wcp_candidates[feasible])) {

                            n_threads <-  n_chains * threads_per_chain
                            measured_algorithm <-  algorithm
                            benchmark_role <-  "main_scaling"
                            seed_index <-  chain_index
                            ##
                            ## Retain PS2's Stan seed offset wherever the actual thread total was an original PS2 budget.
                            if (algorithm == "AD_Stan_WCP") {

                                ps2_thread_index <-  match(x = n_threads, table = settings$n_threads_vec)
                                if (!is.na(x = ps2_thread_index)) seed_index <-  ps2_thread_index

                            }
                            ##
                            if (wcp_algorithm) parallel_case_found <-  TRUE
                            ##
                            ##
                            ## ---- Iterations of the (long) timed run and, for two-run timing, of the short run:
                            ##
                            ## Mplus reads its explicit iteration setting and selected mode. The old rule FBITERATIONS = 2 x bayesmvp_iterations
                            ## gave 80 / 20 / 4, which Mplus silently ran as 100 iterations per chain.
                            ##
                            # n_iter <-  fn_paper1_iterations_for_case( settings  = settings,
                            #                                           algorithm = algorithm,
                            #                                           N         = N,
                            #                                           which_run = "long_run")
                            n_iter <-  fn_paper1_iterations_for_case( settings          = settings,
                                                                      algorithm         = algorithm,
                                                                      N                 = N,
                                                                      which_run         = "long_run",
                                                                      n_chains          = n_chains,
                                                                      num_chunks        = num_chunks,
                                                                      threads_per_chain = threads_per_chain)
                            ##
                            n_iter_short_run <-  if (timing_settings$timing_method == "two_run_difference") {

                                fn_paper1_iterations_for_case( settings  = settings,
                                                               algorithm = algorithm,
                                                               N         = N,
                                                               which_run = "short_run")

                            } else NA_real_
                            ##
                            base_seed <-  repeat_index * 1000
                            actual_seed <-  base_seed + if (grepl(pattern = "^AD_", x = algorithm)) seed_index * 1000 else 0
                            ##
                            rows[[length(x = rows) + 1]] <-  data.frame( device = settings$device,
                                                                       algorithm = measured_algorithm,
                                                                       N = N,
                                                                       num_chunks = num_chunks,
                                                                       n_threads = n_threads,
                                                                       n_chains = n_chains,
                                                                       threads_per_chain = threads_per_chain,
                                                                       n_iter = n_iter,
                                                                       n_iter_short_run = n_iter_short_run,
                                                                       timing_method = timing_settings$timing_method,
                                                                       run = repeat_index,
                                                                       thread_index = seed_index,
                                                                       base_seed = base_seed,
                                                                       seed = actual_seed,
                                                                       benchmark_role = benchmark_role,
                                                                       stringsAsFactors = FALSE)

                        }

                    }

                }
                ##
                if (length(x = rows) == rows_before_N || !parallel_case_found) {

                    stop("No feasible chunk/WCP combinations for ", algorithm, " at N = ", N)

                }

            }

        }
        ##
        grid <-  do.call(what = rbind, args = rows)
        ## WCP-only uses one target partition per WCP thread, without an additional chunk search.
        ## Reuse overlapping measurements; add only the missing sparse WCP-only configurations.
        if (stan_backend == "NicoStan" && isTRUE(x = settings$include_NicoStan_WCP_only)) {

            wcp_only <-  grid[grid$algorithm == "AD_Stan_WCP", , drop = FALSE]
            wcp_only$num_chunks <-  wcp_only$threads_per_chain
            ## 2026-10-03: the long-run count can depend on the allocation (many-chain WCP-only cases; see
            ## fn_paper1_iterations_for_case), so it is recomputed for these rows (unchanged for 4-16 chains):
            wcp_only$n_iter <-  mapply( FUN = function(N, n_chains, num_chunks, threads_per_chain) {
                                          fn_paper1_iterations_for_case( settings          = settings,
                                                                         algorithm         = "AD_Stan_WCP",
                                                                         N                 = N,
                                                                         which_run         = "long_run",
                                                                         n_chains          = n_chains,
                                                                         num_chunks        = num_chunks,
                                                                         threads_per_chain = threads_per_chain)
                                      },
                                      wcp_only$N, wcp_only$n_chains,
                                      wcp_only$num_chunks, wcp_only$threads_per_chain)
            grid <-  rbind(grid, wcp_only)
            grid <-  grid[order(match(x = grid$algorithm, table = settings$algorithms)), , drop = FALSE]

        }
        ##
        ## ---- BayesMVP WCP-only rows:
        ##
        ## Same rule as the Stan arm: one chunk per WCP thread (num_chunks = n_WCP), so that EVERY WCP count has a WCP-only
        ## measurement. Without this, WCP-only exists only where the WCP count happens to be in the chunk candidate list
        ## (2 and 4 at small N, 4 at large N), and the within-BayesMVP comparison (no chunking / chunking only / WCP only /
        ## chunking + WCP, best per N) cannot select a WCP-only optimum. Overlapping measurements are deduplicated below.
        ##
        if (isTRUE(x = settings$include_BayesMVP_WCP_only)) {

            bayesmvp_wcp_only <-  grid[grid$algorithm == "MD_BayesMVP_WCP", , drop = FALSE]
            bayesmvp_wcp_only$num_chunks <-  bayesmvp_wcp_only$threads_per_chain
            ## 2026-10-03: the long-run count can depend on the allocation (many-chain WCP-only cases; see
            ## fn_paper1_iterations_for_case), so it is recomputed for these rows (unchanged for 4-16 chains):
            bayesmvp_wcp_only$n_iter <-  mapply( FUN = function(N, n_chains, num_chunks, threads_per_chain) {
                                          fn_paper1_iterations_for_case( settings          = settings,
                                                                         algorithm         = "MD_BayesMVP_WCP",
                                                                         N                 = N,
                                                                         which_run         = "long_run",
                                                                         n_chains          = n_chains,
                                                                         num_chunks        = num_chunks,
                                                                         threads_per_chain = threads_per_chain)
                                      },
                                      bayesmvp_wcp_only$N, bayesmvp_wcp_only$n_chains,
                                      bayesmvp_wcp_only$num_chunks, bayesmvp_wcp_only$threads_per_chain)
            grid <-  rbind(grid, bayesmvp_wcp_only)
            grid <-  grid[order(match(x = grid$algorithm, table = settings$algorithms)), , drop = FALSE]

        }
        grid <-  fn_paper1_execution_columns(cases = grid, stan_backend = stan_backend,
                                             native_backend = fn_paper1_native_backend(settings = settings))
        ## Each physical configuration/repeat is measured once.
        grid <-  grid[!duplicated(x = fn_paper1_case_identity(grid = grid)), , drop = FALSE]
        rownames(x = grid) <-  NULL
        grid$case_id <-  seq_len(length.out = nrow(x = grid))
        grid$timing_estimator <-  if (timing_settings$timing_method == "single_run") "single_run" else
            ifelse(grepl("^Mplus_", grid$algorithm) & identical(settings$mplus_short_run_role, "overhead_control"),
                   "overhead_subtraction", "iteration_difference")
        grid$mplus_iteration_mode <-  ifelse(grepl("^Mplus_", grid$algorithm), settings$mplus_iteration_mode[1], NA_character_)
        ##
        ## ---- Several Mplus long-run modes: every Mplus case is repeated once per extra mode, so FBITERATIONS
        ## ---- and BITERATIONS runs sit side by side in results.rds, each with its own cache key.
        ##
        for (extra_mplus_mode in settings$mplus_iteration_mode[-1]) {

            extra_mplus_rows <-  grid[grepl("^Mplus_", grid$algorithm) & grid$mplus_iteration_mode == settings$mplus_iteration_mode[1] &
                                      !is.na(grid$mplus_iteration_mode), , drop = FALSE]
            if (nrow(x = extra_mplus_rows)) {

                extra_mplus_rows$mplus_iteration_mode <-  extra_mplus_mode
                grid <-  rbind(grid, extra_mplus_rows)

            }

        }
        grid <-  grid[order(match(x = grid$algorithm, table = settings$algorithms)), , drop = FALSE]
        rownames(x = grid) <-  NULL
        grid$case_id <-  seq_len(length.out = nrow(x = grid))
        grid$mplus_iteration_mode_short_run <-  if (timing_settings$timing_method == "two_run_difference")
            ifelse(grepl("^Mplus_", grid$algorithm), settings$mplus_iteration_mode_short_run, NA_character_) else NA_character_
        grid <-  fn_paper1_stan_partition_columns(cases = grid)
        ##
        if (skipped_cases > 0) message("Excluded ", skipped_cases, " candidate repeats exceeding the thread limit or chunk count.")
        ##
        return(grid)

}
##
## ---- Stan preparation: retain the existing PS2 model, priors and compilation flags ----------------------------------------------------
##
fn_paper1_stan_data <-  function( y ) {

        N <-  nrow(x = y)
        ##
        n_tests <-  ncol(x = y)
        ##
        prior_mean <-  prior_sd <-  design <-  list()
        ##
        for (class_index in 1:2) {

            prior_mean[[class_index]] <-  matrix(data = 0, nrow = 1, ncol = n_tests)
            ##
            prior_sd[[class_index]] <-  matrix(data = 1, nrow = 1, ncol = n_tests)

        }
        ##
        prior_mean[[1]][1, 1] <-  -2.10
        ##
        prior_mean[[2]][1, 1] <-  0.40
        ##
        prior_sd[[1]][1, 1] <-  0.45
        ##
        prior_sd[[2]][1, 1] <-  0.375
        ##
        for (test_index in seq_len(length.out = n_tests)) design[[test_index]] <-  matrix(data = 1, nrow = N, ncol = 1)
        ##
        return(list( N = N,
                     n_tests = n_tests,
                     y = y,
                     n_class = 2L,
                     n_pops = 1L,
                     pop = rep(x = 1L, times = N),
                     n_covariates_max_nd = 1L,
                     n_covariates_max_d = 1L,
                     n_covariates_max = 1L,
                     X_nd = design,
                     X_d = design,
                     n_covs_per_outcome = array(data = 1L, dim = c(2, n_tests)),
                     corr_force_positive = 0L,
                     known_num = 0L,
                     overflow_threshold = 5,
                     underflow_threshold = -5,
                     prior_only = 0L,
                     prior_beta_mean = prior_mean,
                     prior_beta_sd = prior_sd,
                     prior_LKJ = matrix(data = c(12, 3), ncol = 1),
                     prior_p_alpha = matrix(data = 5, ncol = 1),
                     prior_p_beta = matrix(data = 10, ncol = 1),
                     Phi_type = 1L,
                     handle_numerical_issues = 1L,
                     fully_vectorised = 1L))

}
##
fn_paper1_compile_stan <-  function( settings,
                                     algorithms
) {

        compile_environment <-  new.env(parent = globalenv())
        ##
        sys.source( file = file.path(settings$algorithm_study_dir, "0_utilities", "shared_functions", "R_fn_compile_Stan_model.R"),
                    envir = compile_environment)
        ##
        model_names <-  c(AD_Stan = "LC_MVP_bin_PartialLog_v5.stan",
                          AD_Stan_chunked = "LC_MVP_bin_PartialLog_v5_chunked.stan",
                          AD_Stan_tape_chunked = "LC_MVP_bin_PartialLog_v5_reduce_sum_static.stan",
                          AD_Stan_WCP = "LC_MVP_bin_PartialLog_v5_reduce_sum_static.stan")
        ##
        compiled <-  models <-  list()
        ##
        for (algorithm in algorithms) {

            file_name <-  model_names[[algorithm]]
            ##
            if (is.null(x = compiled[[file_name]])) {

                threaded <-  grepl(pattern = "reduce_sum", x = file_name)
                ##
                model_file <-  if (threaded) {

                    file.path(settings$algorithm_study_dir, "paper_1_chunking_and_parallel_scalability", "stan_models", file_name)

                } else system.file("stan_models", file_name, package = "BayesMVP", mustWork = TRUE)
                ##
                if (!file.exists(file = model_file)) stop("Missing Paper 1 Stan model: ", model_file)
                ##
                output <-  compile_environment$R_fn_compile_Stan_model( force_recompile = settings$force_recompile,
                                                                        stan_threads = threaded,
                                                                        Stan_model_file_path = model_file,
                                                                        custom_cpp_user_header_file_path = NULL,
                                                                        set_custom_optimised_CXX_CPP_flags = TRUE,
                                                                        CXX_COMPILER_PATH = if (.Platform$OS.type == "unix")
                                                                            "/opt/AMD/aocc-compiler-5.0.0/bin/clang++" else "g++",
                                                                        CPP_COMPILER_PATH = if (.Platform$OS.type == "unix")
                                                                            "/opt/AMD/aocc-compiler-5.0.0/bin/clang" else "gcc",
                                                                        MATH_FLAGS = "-fno-math-errno -fno-signed-zeros -fno-trapping-math",
                                                                        FMA_FLAGS = "-mfma",
                                                                        AVX_FLAGS = if (settings$device == "HPC")
                                                                            "-mavx -mavx2 -mavx512f -mavx512vl -mavx512dq" else "-mavx -mavx2",
                                                                        THREAD_FLAGS = if (threaded)
                                                                            "-pthread -D_REENTRANT -DSTAN_THREADS" else "-D_REENTRANT")
                ##
                if (threaded && !isTRUE(x = cmdstanr:::model_compile_info(output$mod$exe_file())[["STAN_THREADS"]])) {

                    stop("The reduce_sum model must be compiled with STAN_THREADS=true.")

                }
                ##
                compiled[[file_name]] <-  output$mod

            }
            ##
            models[[algorithm]] <-  compiled[[file_name]]

        }
        ##
        return(models)

}
##
## ---- Compile Stan evaluators once and reuse data-specific initialisation in the threaded backend -------------------------------------
##
fn_paper1_compile_stan_via_NicoStan <-  function( settings,
                                                  algorithms
) {

        if (!requireNamespace(package = "bridgestan", quietly = TRUE)) stop("The NicoStan backend requires bridgestan.")
        configuration <-  fn_paper1_stan_via_NicoStan_settings(settings = settings)
        page_release_linux <- identical(x = Sys.info()[["sysname"]], y = "Linux")
        if (page_release_linux) {
            original_bridgestan_path <- getFromNamespace(
                "get_bridgestan_path", "bridgestan")(download = FALSE)
            page_release_runtime <- NicoStan:::R_fn_bridge_memory_prepare_runtime(
                original_bridgestan_path = original_bridgestan_path)
            bridgestan::set_bridgestan_path(page_release_runtime$bridgestan_path)
            on.exit(bridgestan::set_bridgestan_path(original_bridgestan_path), add = TRUE)
        }
        ##
        model_names <-  c(AD_Stan = "LC_MVP_bin_PartialLog_v5.stan",
                          AD_Stan_chunked = "LC_MVP_bin_PartialLog_v5_chunked.stan",
                          AD_Stan_tape_chunked = "LC_MVP_bin_PartialLog_v5_reduce_sum_static.stan",
                          AD_Stan_WCP = "LC_MVP_bin_PartialLog_v5_reduce_sum_static.stan")
        ##
        ## ---- Compiled-model cache:
        ## The libraries live in settings$stan_model_cache_dir (set in the runner, one folder per machine) instead of the
        ## timestamped study folder, so they are compiled once and reused by later studies. A model is recompiled only when
        ## its .stan source, the make arguments or the BridgeStan version change (checked by a fingerprint file per model).
        ##
        if (!is.character(x = settings$stan_model_cache_dir) || length(x = settings$stan_model_cache_dir) != 1 || !nzchar(x = settings$stan_model_cache_dir)) {

            stop("Set paper1_settings$stan_model_cache_dir in the Paper 1 runner.")

        }
        model_directory <-  settings$stan_model_cache_dir
        dir.create(path = model_directory, recursive = TRUE, showWarnings = FALSE)
        compiled <-  models <-  list()
        ##
        for (algorithm in algorithms) {

            file_name <-  model_names[[algorithm]]
            if (is.null(x = file_name)) stop("Unknown Stan arm: ", algorithm)
            ##
            if (is.null(x = compiled[[file_name]])) {

                original_file <-  if (grepl(pattern = "reduce_sum", x = file_name)) {

                    file.path(settings$algorithm_study_dir, "paper_1_chunking_and_parallel_scalability", "stan_models", file_name)

                } else system.file("stan_models", file_name, package = "BayesMVP", mustWork = TRUE)
                ##
                model_file <-  file.path(model_directory, file_name)
                model_library_cached <-  paste0(tools::file_path_sans_ext(x = model_file), "_model.so")
                compile_fingerprint_file <-  paste0(tools::file_path_sans_ext(x = model_file), "_compile_fingerprint.txt")
                compile_make_args <- configuration$make_args
                compile_stanc_args <- NULL
                memory_compile_options <- NULL
                if (page_release_linux) {
                    memory_compile_options <- NicoStan:::R_fn_bridge_memory_compile_options(
                        model_file = model_file, make_args = configuration$make_args)
                    compile_make_args <- memory_compile_options$make_args
                    compile_stanc_args <- memory_compile_options$stanc_args
                }
                compile_fingerprint <-  paste0( "stan_source_md5 = ", unname(obj = tools::md5sum(files = original_file)),
                                                "; make_args = ", paste(unlist(x = compile_make_args), collapse = " "),
                                                "; bridgestan_version = ", as.character(x = utils::packageVersion(pkg = "bridgestan")))
                if (page_release_linux) {
                    compile_fingerprint <- paste0(compile_fingerprint,
                        "; memory_release_abi = v1; memory_release_header_md5 = ",
                        memory_compile_options$release_header_md5,
                        "; page_release_runtime_key = ", page_release_runtime$runtime_key)
                }
                ##
                model_has_memory_release <- !page_release_linux ||
                    NicoStan:::R_fn_bridge_memory_has_export(model_library_cached)
                cached_library_is_valid <-  model_has_memory_release && file.exists(model_library_cached) &&
                                            file.exists(compile_fingerprint_file) &&
                                            identical(x = readLines(con = compile_fingerprint_file, warn = FALSE), y = compile_fingerprint)
                ##
                if (cached_library_is_valid) {

                    message(paste0("Reusing compiled Stan evaluator for NicoStan: ", file_name, " (", model_directory, ")"))
                    model_library <-  model_library_cached

                } else {

                    ## Copy the source only when recompiling: rewriting it would change its timestamp and make BridgeStan rebuild.
                    if (!file.copy(from = original_file, to = model_file, overwrite = TRUE)) stop("Could not copy Stan source: ", original_file)
                    ## STAN_THREADS is required for concurrent autodiff in every arm, including models without reduce_sum.
                    message(paste0("Compiling Stan evaluator for NicoStan: ", file_name))
                    unlink(x = c(model_library_cached, paste0(tools::file_path_sans_ext(x = model_file), ".o")))
                    model_library <-  bridgestan::compile_model(stan_file = model_file,
                                                                stanc_args = compile_stanc_args,
                                                                make_args = compile_make_args)
                    if (page_release_linux && !NicoStan:::R_fn_bridge_memory_has_export(model_library)) {
                        stop("Compiled model is missing nicostan_bs_release_thread_memory_v1: ", model_library)
                    }
                    writeLines(text = compile_fingerprint, con = compile_fingerprint_file)

                }
                ##
                if (page_release_linux) {
                    writeLines(text = paste0("runtime_key=", page_release_runtime$runtime_key,
                                             "; release_header_md5=", memory_compile_options$release_header_md5),
                               con = paste0(tools::file_path_sans_ext(model_file), "_page_release_runtime.txt"))
                }
                ##
                compiled[[file_name]] <-  list( stan_source_path = normalizePath(path = model_file, mustWork = TRUE),
                                                original_source_path = normalizePath(path = original_file, mustWork = TRUE),
                                                stan_source_md5 = unname(obj = tools::md5sum(files = model_file)),
                                                library_path = normalizePath(path = model_library, mustWork = TRUE),
                                                library_md5 = unname(obj = tools::md5sum(files = model_library)),
                                                bridgestan_version = as.character(x = utils::packageVersion(pkg = "bridgestan")),
                                                make_args = configuration$make_args,
                                                compile_make_args = compile_make_args,
                                                compile_stanc_args = compile_stanc_args,
                                                memory_release_symbol = if (page_release_linux)
                                                    "nicostan_bs_release_thread_memory_v1" else NULL,
                                                memory_release_method = if (page_release_linux)
                                                    "MADV_DONTNEED" else NULL,
                                                page_release_runtime_key = if (page_release_linux)
                                                    page_release_runtime$runtime_key else NULL)

            }
            ##
            models[[algorithm]] <-  compiled[[file_name]]

        }
        ##
        return(models)

}
##
fn_paper1_run_stan_via_NicoStan <-  function( case,
                                              y,
                                              stan_data,
                                              settings,
                                              runtime
) {

        configuration <-  fn_paper1_stan_via_NicoStan_settings(settings = settings)
        model_specification <-  runtime$stan_models[[case$algorithm]]
        model_cache <-  runtime$stan_via_NicoStan_cache
        if (!is.environment(x = model_cache)) stop("The threaded Stan backend requires an execution-local model cache.")
        ## Chunk size is model data, so a new N or chunk size needs its own initialisation; chains, WCP and repeats do not.
        chunk_size <-  if (case$algorithm == "AD_Stan") "none" else as.character(x = stan_data$chunk_size)
        model_key <-  paste0(model_specification$stan_source_md5, "_N_", case$N, "_chunk_size_", chunk_size)
        model <-  model_cache[[model_key]]
        ##
        if (is.null(x = model)) {

            message(paste0("Initialising Stan evaluator for NicoStan: ", case$algorithm,
                           ", N = ", case$N, ", chunk_size = ", chunk_size))
            ##
            initialised_model <-  NicoStan:::initialise_model( Model_type = "Stan",
                                                               stream = case$base_seed,
                                                               sample_nuisance = TRUE,
                                                               n_nuisance_override = nrow(x = y) * ncol(x = y),
                                                               model_args_list = list(y = y),
                                                               Stan_data_list = stan_data,
                                                               Stan_model_file_path = model_specification$stan_source_path,
                                                               stanc_args = model_specification$compile_stanc_args,
                                                               make_args = if (is.null(x = model_specification$compile_make_args))
                                                                   model_specification$make_args else
                                                                   model_specification$compile_make_args)
            ##
            model_info <-  initialised_model$bs_model$model_info()
            if (!any(grepl(pattern = "STAN_THREADS[[:space:]]*=[[:space:]]*true", x = model_info, ignore.case = TRUE))) {

                stop("Every Stan evaluator used by parallel BayesMVP chains must be compiled with STAN_THREADS=true.")

            }
            ## The current six-test models declare u_raw first, followed by the 43 main coordinates.
            if (initialised_model$n_nuisance != nrow(x = y) * ncol(x = y) ||
                initialised_model$n_params_main != 2 * choose(n = ncol(x = y), k = 2) + 2 * ncol(x = y) + 1) {

                stop("The Stan evaluator dimensions do not match the Paper 1 model.")

            }
            if (!all(startsWith(x = initialised_model$bs_main_param_names[seq_len(length.out = initialised_model$n_nuisance)],
                                 prefix = "u_raw."))) stop("The Paper 1 nuisance block must be the first Stan parameter declaration.")
            ## All four models require beta[2,1,1] >= beta[1,1,1]. Locate their free beta_vec coordinates explicitly.
            ordered_beta_indices <-  match(x = c("beta_vec.1", paste0("beta_vec.", ncol(x = y) + 1L)),
                                            table = initialised_model$bs_main_param_names) - initialised_model$n_nuisance
            if (anyNA(x = ordered_beta_indices) || any(ordered_beta_indices < 1L) ||
                any(ordered_beta_indices > initialised_model$n_params_main)) stop("Could not locate the Stan class-order coordinates.")
            ## Keep the JSON beside the study's evaluator, independent of later package-cache writes.
            json_file <-  file.path(settings$output_dir, "stan_via_NicoStan_models", paste0(model_key, ".json"))
            ## The compiled models now live in the machine cache, so this study folder no longer exists by default.
            dir.create(path = dirname(path = json_file), recursive = TRUE, showWarnings = FALSE)
            if (file.exists(file = json_file)) {

                # if (!identical(unname(tools::md5sum(initialised_model$json_file_path)), unname(tools::md5sum(json_file)))) {
                ## 2026-10-03: compare the data, not the bytes: JSON files written before 29 Sep are compact and newer
                ## ones are pretty-printed, so identical data gave different md5 sums (and stopped the run).
                if (!identical(jsonlite::fromJSON(txt = initialised_model$json_file_path, simplifyVector = FALSE),
                               jsonlite::fromJSON(txt = json_file, simplifyVector = FALSE))) {

                    stop("Existing Stan JSON belongs to different data; use a fresh study output directory.")

                }

            } else if (!file.copy(from = initialised_model$json_file_path, to = json_file, overwrite = FALSE)) {

                stop("Could not preserve the Stan data JSON: ", json_file)

            }
            ##
            model_args <-  initialised_model$Model_args_as_Rcpp_List
            model_args$json_file_path <-  normalizePath(path = json_file, mustWork = TRUE)
            model_args$N <-  as.integer(x = case$N)
            model_args$n_tests <-  as.integer(x = ncol(x = y))
            ## Retain only the sampler inputs, not the large BridgeStan parameter-name lists.
            model <-  list( model_args = model_args,
                             n_nuisance = initialised_model$n_nuisance,
                             n_params_main = initialised_model$n_params_main,
                             ordered_beta_indices = ordered_beta_indices)
            model_cache[[model_key]] <-  model

        }
        ##
        runtime_receipt_after_evaluator_initialisation <-  fn_paper1_stan_runtime_receipt(
            stage = "evaluator_initialisation", require_single_tbb = FALSE)
        ##
        EHMC_args <-  fn_paper1_stan_HMC_args(configuration = configuration)
        ##
        ## Record the L and randomize_tau that NicoStan's converter actually received, outside the clock.
        ##
        checked_trajectory <-  fn_paper1_check_fixed_trajectory( EHMC_args    = EHMC_args,
                                                                 fixed_L      = configuration$fixed_L,
                                                                 package_name = "NicoStan")
        ##
        metric <-  NicoStan:::init_EHMC_Metric_as_Rcpp_List( n_params_main = model$n_params_main,
                                                           n_nuisance = model$n_nuisance,
                                                           metric_shape_main = configuration$metric_shape_main)
        ## Time state allocation, thread setup, the threaded sampler (including C++ model loading), and divergence checking.
        started <-  proc.time()[3]
        initial <-  withr::with_seed(seed = case$seed, code = {

            list( main = matrix(data = stats::runif(n = model$n_params_main * case$n_chains,
                                                    min = -configuration$init_radius, max = configuration$init_radius),
                                 nrow = model$n_params_main, ncol = case$n_chains),
                  nuisance = matrix(data = stats::runif(n = model$n_nuisance * case$n_chains,
                                                        min = -configuration$init_radius, max = configuration$init_radius),
                                     nrow = model$n_nuisance, ncol = case$n_chains))

        })
        ## Sorting these two iid uniform draws samples their conditional distribution under the model's ordering constraint.
        beta_first <-  initial$main[model$ordered_beta_indices[1], ]
        beta_second <-  initial$main[model$ordered_beta_indices[2], ]
        initial$main[model$ordered_beta_indices[1], ] <-  pmin(beta_first, beta_second)
        initial$main[model$ordered_beta_indices[2], ] <-  pmax(beta_first, beta_second)
        ## Stan reduce_sum and the chain workers share this total TBB budget; WCP is not a separate per-chain cap.
        RcppParallel::setThreadOptions(numThreads = case$n_threads)
        output <-  NicoStan:::Rcpp_fn_RcppParallel_EHMC_sampling(  n_threads_R = case$n_chains,
                                                                   n_threads_WCP = case$threads_per_chain,
                                                                   sample_nuisance_R = TRUE,
                                                                   ## The Stan arms keep the FULL nuisance trace on purpose: Stan stores every parameter every
                                                                   ## iteration, so dropping it would make Stan look better than it is.
                                                                   n_nuisance_to_track = model$n_nuisance,
                                                                   seed_R = case$seed,
                                                                   iter_one_by_one = FALSE,
                                                                   n_iter_R = case$n_iter,
                                                                   partitioned_HMC_R = FALSE,
                                                                   diffusion_HMC_R = FALSE,
                                                                   Model_type_R = "Stan",
                                                                   force_autodiff_R = TRUE,
                                                                   force_PartialLog_R = FALSE,
                                                                   multi_attempts_R = FALSE,
                                                                   theta_main_vectors_all_chains_input_from_R = initial$main,
                                                                   theta_us_vectors_all_chains_input_from_R = initial$nuisance,
                                                                   y_Eigen_R = y,
                                                                   Model_args_as_Rcpp_List = model$model_args,
                                                                   EHMC_args_as_Rcpp_List = EHMC_args,
                                                                   EHMC_Metric_as_Rcpp_List = metric,
                                                                   use_disk = FALSE,
                                                                   trace_dir = "/tmp/hmc_traces")
        ##
        divergences <-  sum(unlist(x = output[[2]]))
        ## Match the existing CmdStanR arm: retain divergence counts for inspection rather than changing its failure policy.
        elapsed <-  unname(obj = proc.time()[3] - started)
        ## The second receipt remains outside the timer and rejects any model-load path that introduced another TBB runtime.
        runtime_receipt_after_timed_call <-  fn_paper1_stan_runtime_receipt(
            stage = "timed_call", require_single_tbb = FALSE)
        if (!identical(x = runtime_receipt_after_evaluator_initialisation$mapped_tbb_path,
                       y = runtime_receipt_after_timed_call$mapped_tbb_path) ||
            !identical(x = runtime_receipt_after_evaluator_initialisation$mapped_tbb_md5,
                       y = runtime_receipt_after_timed_call$mapped_tbb_md5)) {

            invisible(NULL)

        }
        ##
        ## ---- Self-check: the Stan arms must really have stored the full nuisance trace (output[[3]] holds it).
        ##
        n_nuisance_tracked <-  model$n_nuisance
        if (!all(vapply(X = output[[3]], FUN = function(trace) nrow(x = trace) == model$n_nuisance, FUN.VALUE = logical(length = 1)))) {

            stop("The Stan arm did not store the full nuisance trace (Stan case ", case$case_id, ").")

        }
        ##
        return(list( elapsed_seconds = elapsed,
                     divergences = divergences,
                     mean_L = as.numeric(x = checked_trajectory$L_main),
                     randomize_tau = checked_trajectory$randomize_tau,
                     n_nuisance_to_track = n_nuisance_tracked,
                     stan_runtime_receipt_status = runtime_receipt_after_timed_call$verification,
                     stan_runtime_tbb_path_after_evaluator_initialisation = runtime_receipt_after_evaluator_initialisation$mapped_tbb_path,
                     stan_runtime_tbb_md5_after_evaluator_initialisation = runtime_receipt_after_evaluator_initialisation$mapped_tbb_md5,
                     stan_runtime_tbb_path_after_timed_call = runtime_receipt_after_timed_call$mapped_tbb_path,
                     stan_runtime_tbb_md5_after_timed_call = runtime_receipt_after_timed_call$mapped_tbb_md5,
                     timing_scope = "NicoStan_BridgeStan_state_allocation_thread_setup_sampling_model_load_divergence_check"))

}
##
## ---- Shared initial state for every BayesMVP cell ---------------------------------------------------------------------------------------
##
fn_paper1_bayesmvp_initial_state <-  function( y,
                                               n_chains,
                                               sampler_settings
) {

        n_tests <-  ncol(x = y)
        ##
        n_corrs <-  2 * choose(n = n_tests, k = 2)
        ##
        n_params_main <-  n_corrs + 2 * n_tests + 1
        ##
        main <-  matrix(data = 0.01, nrow = n_params_main, ncol = n_chains)
        ##
        main[n_corrs + seq_len(length.out = n_tests), ] <-  -1
        ##
        main[n_corrs + n_tests + seq_len(length.out = n_tests), ] <-  1
        ##
        main[n_params_main, ] <-  sampler_settings$init_prevalence_raw
        ##
        nuisance <-  matrix(data = sampler_settings$init_u_value, nrow = nrow(x = y) * n_tests, ncol = n_chains)
        ##
        return(list(main = main,
                    nuisance = nuisance))

}
##
## ---- One case executor, shared by the chunk and scaling analyses ----------------------------------------------------------------------
##
fn_paper1_run_case <-  function( case,
                                 y,
                                 settings,
                                 runtime
) {

        if (case$n_threads != case$n_chains * case$threads_per_chain ||
            case$n_threads > min(settings$total_thread_limit, parallel::detectCores())) stop("Invalid case thread budget.")
        ##
        if (case$algorithm %in% c("MD_BayesMVP_WCP", "AD_Stan_WCP", "Mplus_WCP") && case$threads_per_chain <= 1) {

            stop("Use the non-WCP arm for one thread per chain; WCP arms require at least two.")

        }
        ##
        if (case$N != nrow(x = y) || ncol(x = y) != 6L) stop("Case/data dimensions disagree.")
        ##
        if (grepl(pattern = "^MD_", x = case$algorithm) && case$threads_per_chain > case$num_chunks) {

            stop("BayesMVP WCP threads must not exceed chunk count.")

        }
        ##
        if (grepl(pattern = "^MD_", x = case$algorithm)) {

            n_tests <-  ncol(x = y)
            ##
            n_nuisance <-  nrow(x = y) * n_tests
            ##
            n_params_main <-  2 * choose(n = n_tests, k = 2) + 2 * n_tests + 1
            ##
            sampler <-  settings$bayesmvp
            ##
            ## Nuisance trace kept during the timed call: 0 (none) or "all" (see the runner). NicoStan allows no partial trace.
            if (is.null(x = sampler$n_nuisance_to_track) || length(x = sampler$n_nuisance_to_track) != 1L ||
                !is.numeric(x = sampler$n_nuisance_to_track) || !is.finite(x = sampler$n_nuisance_to_track) || sampler$n_nuisance_to_track < 0) {

                stop("Set paper1_settings$bayesmvp$n_nuisance_to_track (0 = keep no nuisance trace during the timed calls).")

            }
            n_nuisance_tracked <-  if (sampler$n_nuisance_to_track == 0) 0 else n_nuisance
            ##
            EHMC_args <-  BayesMVP:::init_EHMC_args_as_Rcpp_List(diffusion_HMC = sampler$diffusion_HMC)
            ##
            EHMC_args$eps_main <-  sampler$eps_main
            ##
            EHMC_args$tau_main <-  sampler$L_main * sampler$eps_main
            ##
            ## ---- Exactly L_main leapfrog steps every iteration:
            ##
            ## Before this fixed-tau change, these arms did not set randomize_tau; the C++ default was true (NicoStan structures.hpp),
            ## so every iteration drew tau_ii ~ U(0, 2 tau_main), i.e. L between 1 and 2 L_main instead of L_main.
            ## tau_main_ii is set as well because the converter check computes L from it.
            ##
            EHMC_args$tau_main_ii <-  EHMC_args$tau_main
            if (!identical(sampler$randomize_tau, FALSE)) stop("Paper 1 requires bayesmvp$randomize_tau = FALSE for fixed L.")
            EHMC_args$randomize_tau <-  sampler$randomize_tau
            ##
            checked_trajectory <-  fn_paper1_check_fixed_trajectory( EHMC_args    = EHMC_args,
                                                                     fixed_L      = sampler$L_main,
                                                                     package_name = "BayesMVP")
            ##
            metric <-  BayesMVP:::init_EHMC_Metric_as_Rcpp_List( n_params_main = n_params_main,
                                                                 n_nuisance = n_nuisance,
                                                                 metric_shape_main = sampler$metric_shape_main)
            ##
            ## Each N has one fixed dataset in this study; reuse its model arguments across cases and BayesMVP arms.
            model_cache <-  runtime$bayesmvp_model_args
            model_key <-  as.character(x = case$N)
            model_args <-  if (is.environment(x = model_cache)) model_cache[[model_key]] else NULL
            ##
            if (is.null(x = model_args)) {

                model <-  BayesMVP:::initialise_model( Model_type = "LC_MVP",
                                                       stream = case$base_seed,
                                                       sample_nuisance = TRUE,
                                                       n_nuisance_override = NULL,
                                                       model_args_list = list(y = y),
                                                       compile = TRUE,
                                                       force_recompile = FALSE,
                                                       cmdstanr_model_fit_obj = NULL,
                                                       Stan_data_list = NULL,
                                                       Stan_model_file_path = NULL,
                                                       Stan_cpp_user_header = NULL,
                                                       Stan_cpp_flags = NULL,
                                                       stanc_args = NULL)
                ##
                model_args <-  model$Model_args_as_Rcpp_List
                ## Cache the unmodified arguments; chunk count and other case overrides below remain local to this call.
                if (is.environment(x = model_cache)) model_cache[[model_key]] <-  model_args

            }
            ##
            model_args$Model_args_ints[4] <-  case$num_chunks
            ##
            model_args$model_so_file <-  model_args$json_file_path <-  "none"
            ##
            model_args$Model_args_strings[12] <-  "num_diff"
            ## Timed BayesMVP work includes state allocation, thread setup, sampling and the divergence check.
            started <-  proc.time()[3]
            ##
            initial <-  fn_paper1_bayesmvp_initial_state(y = y, n_chains = case$n_chains, sampler_settings = sampler)
            ##
            call_args <-  list( n_threads_R = case$n_chains,
                                n_threads_WCP = case$threads_per_chain,
                                sample_nuisance_R = TRUE,
                                n_nuisance_to_track = n_nuisance_tracked,
                                seed_R = case$seed,
                                iter_one_by_one = FALSE,
                                n_iter_R = case$n_iter,
                                partitioned_HMC_R = sampler$partitioned_HMC,
                                diffusion_HMC_R = sampler$diffusion_HMC,
                                Model_type_R = "LC_MVP",
                                force_autodiff_R = FALSE,
                                force_PartialLog_R = FALSE,
                                multi_attempts_R = FALSE,
                                theta_main_vectors_all_chains_input_from_R = initial$main,
                                theta_us_vectors_all_chains_input_from_R = initial$nuisance,
                                y_Eigen_R = y,
                                Model_args_as_Rcpp_List = model_args,
                                EHMC_args_as_Rcpp_List = EHMC_args,
                                EHMC_Metric_as_Rcpp_List = metric,
                                use_disk = FALSE,
                                trace_dir = "/tmp/hmc_traces")
            ##
            if (case$algorithm == "MD_BayesMVP_multi_process") {

                if (.Platform$OS.type == "windows") stop("The original multiprocess arm requires fork support.")
                ##
                RcppParallel::setThreadOptions(numThreads = 1)
                ##
                divergence_counts <-  parallel::mclapply( X = seq_len(length.out = case$n_chains),
                                                          FUN = function(chain_index) {

                        chain_args <-  call_args
                        ##
                        chain_args$n_threads_R <-  chain_args$n_threads_WCP <-  1
                        ##
                        chain_args$seed_R <-  case$base_seed + chain_index
                        ##
                        chain_args$theta_main_vectors_all_chains_input_from_R <-  initial$main[, chain_index, drop = FALSE]
                        ##
                        chain_args$theta_us_vectors_all_chains_input_from_R <-  initial$nuisance[, chain_index, drop = FALSE]
                        ##
                        output <-  do.call(what = BayesMVP:::Rcpp_fn_RcppParallel_EHMC_sampling, args = chain_args)
                        ##
                        return(sum(unlist(x = output[[2]])))

                },
                                                          mc.cores = case$n_chains)
                ##
                if (any(vapply(X = divergence_counts, FUN = inherits, FUN.VALUE = logical(length = 1), what = "try-error"))) {

                    stop("A BayesMVP multiprocess worker failed.")

                }
                ##
                divergences <-  sum(unlist(x = divergence_counts))

            } else {

                ## Use the full case thread budget, matching PS2 and BayesMVP's sampling wrapper.
                RcppParallel::setThreadOptions(numThreads = case$n_threads)
                ##
                output <-  do.call(what = BayesMVP:::Rcpp_fn_RcppParallel_EHMC_sampling, args = call_args)
                ##
                divergences <-  sum(unlist(x = output[[2]]))
                ##
                ## ---- Self-check: the nuisance trace must really be absent when none was requested (output[[3]] holds it).
                ##
                if (n_nuisance_tracked == 0 && any(vapply(X = output[[3]], FUN = function(trace) length(x = trace) > 0, FUN.VALUE = logical(length = 1)))) {

                    stop("A nuisance trace was stored although n_nuisance_to_track = 0 (BayesMVP case ", case$case_id, ").")

                }

            }
            ##
            if (divergences > 0) stop("Divergences in BayesMVP case ", case$case_id, ".")
            ##
            elapsed <-  unname(obj = proc.time()[3] - started)
            ##
            return(list( elapsed_seconds = elapsed,
                         divergences = divergences,
                         mean_L = as.numeric(x = checked_trajectory$L_main),
                         randomize_tau = checked_trajectory$randomize_tau,
                         n_nuisance_to_track = n_nuisance_tracked,
                         timing_scope = "state_allocation_thread_setup_sampling_divergence_check"))

        }
        ##
        if (grepl(pattern = "^AD_", x = case$algorithm)) {

            stan_data <-  runtime$stan_data[[as.character(x = case$N)]]
            stan_partition <-  fn_paper1_stan_partition_columns(cases = case)
            ##
            if (case$algorithm == "AD_Stan_chunked") stan_data$chunk_size <-  stan_partition$stan_chunk_size
            ##
            if (case$algorithm %in% c("AD_Stan_WCP", "AD_Stan_tape_chunked")) {

                stan_data$chunk_size <-  stan_partition$stan_chunk_size

            }
            ##
            if (fn_paper1_stan_backend(settings = settings) == "NicoStan") {

                return(fn_paper1_run_stan_via_NicoStan(case = case, y = y, stan_data = stan_data,
                                                      settings = settings, runtime = runtime))

            }
            ##
            stan_threads_per_chain <-   if (case$algorithm %in% c("AD_Stan_WCP", 
                                                                  "AD_Stan_tape_chunked")){
                                                case$threads_per_chain 
                                        } else {
                                                NULL 
                                        }
            ##
            started <-  proc.time()[3]
            ##
            fit <-  runtime$stan_models[[case$algorithm]]$sample( data = stan_data,
                                                                  seed = case$seed,
                                                                  chains = case$n_chains,
                                                                  parallel_chains = case$n_chains,
                                                                  threads_per_chain = stan_threads_per_chain,
                                                                  iter_warmup = 0,
                                                                  iter_sampling = case$n_iter,
                                                                  adapt_engaged = FALSE,
                                                                  step_size = settings$stan$step_size,
                                                                  metric = settings$stan$metric,
                                                                  init = settings$stan$init,
                                                                  refresh = 0,
                                                                  save_warmup = FALSE,
                                                                  max_treedepth = settings$stan$max_treedepth)
            ##
            elapsed <-  unname(obj = proc.time()[3] - started)
            ## Diagnostics remain outside the Stan clock. Read named fields rather than relying on array position.
            diagnostics <-  fit$sampler_diagnostics()
            ##
            diagnostic_names <-  dimnames(x = diagnostics)[[3]]
            ##
            divergences <-  if ("divergent__" %in% diagnostic_names) sum(diagnostics[, , "divergent__"]) else NA_real_
            ##
            mean_L <-  if ("n_leapfrog__" %in% diagnostic_names) mean(x = diagnostics[, , "n_leapfrog__"]) else NA_real_
            ##
            return(list(elapsed_seconds = elapsed, divergences = divergences, mean_L = mean_L, timing_scope = "cmdstan_sample_only"))

        }
        ##
        ## Mplus uses its original PS2 model/settings and saves separate input/output files for each case.
        if (!is.function(x = settings$mplus$runner) || !is.function(x = settings$mplus$verifier) || is.null(x = settings$mplus$settings)) {

            stop("Load the Paper 1 Mplus runner, run verifier and settings before selecting Mplus.")

        }
        ##
        mplus_settings <-  settings$mplus$settings
        ##
        mplus_settings$WCP <-  case$algorithm == "Mplus_WCP"
        ##
        mplus_settings$n_threads <-  mplus_settings$n_threads_if_local_HPC_AMD_EPYC <-  case$n_threads
        ##
        mplus_settings$deficit_ratio <-  1
        ##
        mplus_settings$n_chains <-  case$n_chains
        mplus_settings$n_WCP <-  case$threads_per_chain
        mplus_settings$total_thread_limit <-  settings$total_thread_limit
        ##
        mplus_settings$n_thin <-  1
        ##
        mplus_settings$n_iter <-  case$n_iter
        ##
        ## Counts are per-chain iteration requests, not retained posterior draws.
        mplus_settings$n_fb_iter <-  case$n_iter
        ## The untimed warm-up call uses the short-run mode too (BITERATIONS = 1), never the FBITERATIONS long-run mode.
        is_short_call <-  !is.null(case$timing_run_label) &&
            (startsWith(case$timing_run_label, "short_run_") || startsWith(case$timing_run_label, "warm_up_"))
        mplus_settings$iteration_mode <-  if (is_short_call) "BITERATIONS" else case$mplus_iteration_mode
        mplus_settings$biterations_minimum <-  settings$mplus_biterations_minimum
        mplus_settings$bconvergence <-  settings$mplus_bconvergence
        mplus_settings$save_draws <-  settings$mplus_save_draws
        ## An overhead control need not save sampling draws; requested verification still applies to the long run.
        is_overhead_control <-  identical(settings$mplus_short_run_role, "overhead_control") &&
            !is.null(case$timing_run_label) && startsWith(case$timing_run_label, "short_run_")
        ## The untimed warm-up call (BITERATIONS = 1) saves no draws, so it is never verified either.
        is_untimed_warm_up <-  !is.null(case$timing_run_label) && startsWith(case$timing_run_label, "warm_up_")
        mplus_settings$verify_saved_draws <-  settings$mplus_verify_saved_draws && !is_overhead_control && !is_untimed_warm_up
        ##
        mplus_settings$n_superchains <-  min(ceiling(x = case$n_chains / 8), case$n_chains)
        ##
        ## The short and long runs of the two-run timing each keep their own Mplus folder.
        ##
        timing_run_folder_suffix <-  if (is.null(x = case$timing_run_label)) "" else paste0("_", case$timing_run_label)
        mplus_settings$output_dir <-  file.path(settings$output_dir, "Mplus_inputs", case$algorithm,
                                               paste0("N_", case$N, "_chains_", case$n_chains,
                                                      "_WCP_", case$threads_per_chain, "_repeat_", case$run,
                                                      if (!is.null(case$resume_key) && !is.na(case$resume_key)) paste0("_case_", case$resume_key) else "",
                                                      timing_run_folder_suffix))
        ##
        started <-  proc.time()[3]
        ##
        output <-  settings$mplus$runner(y = y, N = case$N, MCMC_seed = case$seed, settings = mplus_settings)
        ##
        elapsed <-  unname(obj = proc.time()[3] - started)
        ##
        if (is.null(x = output)) stop("Mplus returned no result for case ", case$case_id, ".")
        ##
        ## Outside the clock: verify the requested controls; count actual iterations only when explicitly requested.
        ##
        verified_run <-  settings$mplus$verifier( output_file = output$output_file,
                                                  n_fb_iter   = mplus_settings$n_fb_iter,
                                                  n_chains    = case$n_chains,
                                                  iteration_mode = mplus_settings$iteration_mode,
                                                  verify_saved_draws = mplus_settings$verify_saved_draws)
        ##
        return(list(elapsed_seconds = elapsed,
                    divergences = NA_real_,
                    mean_L = NA_real_,
                    mplus_iterations_per_chain = verified_run$iterations_per_chain,
                    mplus_PPPP_run = verified_run$PPPP_run,
                    timing_scope = paste0("Mplus_input_write_engine_run_output_validation_", mplus_settings$iteration_mode)))

}
##
## ---- One execution loop and one results table -----------------------------------------------------------------------------------------
##
## ---- Results files named by device, N and algorithm, following the PS2 convention ------------------------------------------------------
##
fn_paper1_save_results_by_algorithm_and_N <-  function( results,
                                                        output_dir,
                                                        algorithm,
                                                        N
) {

        results_directory <-  file.path(output_dir, "results_by_algorithm_and_N")
        dir.create(path = results_directory, recursive = TRUE, showWarnings = FALSE)
        ##
        selected_results <-  fn_paper1_execution_columns(cases = fn_paper1_stan_partition_columns(cases = results),
                                                         stan_backend = NULL, native_backend = NULL)
        if (!is.null(x = algorithm)) selected_results <-  selected_results[selected_results$algorithm %in% algorithm, , drop = FALSE]
        if (!is.null(x = N)) selected_results <-  selected_results[selected_results$N %in% N, , drop = FALSE]
        ##
        ## Mplus files are split by long-run mode; "" for the other algorithms.
        selected_results$mplus_mode_file_key <-  if ("mplus_iteration_mode" %in% names(x = selected_results))
            ifelse(grepl("^Mplus_", selected_results$algorithm) & !is.na(selected_results$mplus_iteration_mode), selected_results$mplus_iteration_mode, "") else ""
        result_groups <-  unique(x = selected_results[, c("device", "N", "algorithm", "execution_backend", "mplus_mode_file_key"), drop = FALSE])
        saved_files <-  character(length = 0)
        ##
        for (group_index in seq_len(length.out = nrow(x = result_groups))) {

            result_group <-  result_groups[group_index, , drop = FALSE]
            group_results <-  selected_results[selected_results$device == result_group$device &
                                               selected_results$N == result_group$N &
                                               selected_results$algorithm == result_group$algorithm &
                                               selected_results$execution_backend == result_group$execution_backend &
                                               selected_results$mplus_mode_file_key == result_group$mplus_mode_file_key, , drop = FALSE]
            group_results$mplus_mode_file_key <-  NULL
            ##
            file_name <-  paste0(result_group$device, "_paper_1_N_", result_group$N, "_algorithm_", result_group$algorithm)
            if (result_group$execution_backend == "NicoStan_BridgeStan") file_name <-  paste0(file_name, "_backend_NicoStan_BridgeStan")
            if (nzchar(x = result_group$mplus_mode_file_key)) file_name <-  paste0(file_name, "_", result_group$mplus_mode_file_key)
            ##
            ## The timing method is part of every file name; results saved before it existed keep their names.
            ##
            if ("timing_method" %in% names(x = group_results)) {

                file_name <-  paste0(file_name, "_timing_", unique(x = group_results$timing_method)[1])

            }
            results_file <-  file.path(results_directory, paste0(file_name, ".rds"))
            csv_file <-  file.path(results_directory, paste0(file_name, ".csv"))
            ##
            saveRDS(object = group_results, file = results_file)
            ## Keep all rows and fields; only flatten list columns in the readable CSV copy.
            group_results_csv <-  group_results
            for (column in names(x = group_results_csv)[vapply(X = group_results_csv, FUN = is.list, FUN.VALUE = logical(length = 1))]) {

                group_results_csv[[column]] <-  vapply( X = group_results_csv[[column]],
                                                       FUN = paste,
                                                       FUN.VALUE = character(length = 1),
                                                       collapse = ";")

            }
            ##
            utils::write.csv(x = group_results_csv, file = csv_file, row.names = FALSE)
            saved_files <-  c(saved_files, results_file, csv_file)

        }
        ##
        return(invisible(x = saved_files))

}
##
fn_execute_paper1_grid <-  function( grid,
                                     simulated_data,
                                     settings,
                                     runtime,
                                     case_runner,
                                     checkpoint_file,
                                     cache_only = FALSE
) {

        if (anyDuplicated(x = fn_paper1_case_identity(grid = grid))) stop("Duplicate experiment cells must be removed before execution.")
        expected_backends <-  fn_paper1_execution_columns(
            cases = grid[, setdiff(x = names(x = grid), y = "execution_backend"), drop = FALSE],
            stan_backend = fn_paper1_stan_backend(settings = settings),
            native_backend = fn_paper1_native_backend(settings = settings))$execution_backend
        if ("execution_backend" %in% names(x = grid) && !identical(grid$execution_backend, expected_backends)) {

            stop("The grid uses a different execution backend. Regenerate it after changing stan_backend.")

        }
        ## Keep the cache local to this execution so a new study cannot reuse another dataset's model arguments.
        runtime$bayesmvp_model_args <-  new.env(parent = emptyenv())
        runtime$stan_via_NicoStan_cache <-  new.env(parent = emptyenv())
        ##
        results <-  fn_paper1_execution_columns(cases = grid, stan_backend = fn_paper1_stan_backend(settings = settings),
                                                native_backend = fn_paper1_native_backend(settings = settings))
        results <-  fn_paper1_stan_partition_columns(cases = results)
        ##
        results$dataset_md5 <-  NA_character_
        ##
        results$elapsed_seconds <-  results$divergences <-  results$mean_L <-  NA_real_
        ##
        ## ---- Two-run timing columns and the checked trajectory / Mplus iteration records:
        ##
        ## elapsed_seconds is the time every downstream view uses: T_long - S for two-run timing, T for a single run.
        ##
        timing_settings <-  fn_paper1_timing_settings(settings = settings, algorithms = unique(x = grid$algorithm))
        ##
        if (!"timing_method" %in% names(x = grid) || any(grid$timing_method != timing_settings$timing_method)) {

            stop("The grid was built for a different timing_method. Regenerate it after changing settings$timing_method.")

        }
        expected_estimator <-  if (timing_settings$timing_method == "single_run") rep("single_run", nrow(grid)) else
            ifelse(grepl("^Mplus_", grid$algorithm) & identical(settings$mplus_short_run_role, "overhead_control"),
                   "overhead_subtraction", "iteration_difference")
        if (!"timing_estimator" %in% names(grid) || !identical(grid$timing_estimator, expected_estimator)) {
            stop("Regenerate the experiment grid after changing the Mplus short-run role.")
        }
        mplus_rows <- grepl("^Mplus_", grid$algorithm)
        if (any(mplus_rows) &&
            (!"mplus_iteration_mode" %in% names(grid) ||
             anyNA(grid$mplus_iteration_mode[mplus_rows]) ||
             !setequal(unique(grid$mplus_iteration_mode[mplus_rows]), settings$mplus_iteration_mode))) {
            stop("Regenerate the experiment grid after changing the Mplus long-run mode.")
        }
        if (any(mplus_rows) && timing_settings$timing_method == "two_run_difference" &&
            (!"mplus_iteration_mode_short_run" %in% names(grid) ||
             anyNA(grid$mplus_iteration_mode_short_run[mplus_rows]) ||
             any(grid$mplus_iteration_mode_short_run[mplus_rows] != settings$mplus_iteration_mode_short_run))) {
            stop("Regenerate the experiment grid after changing the Mplus short-run mode.")
        }
        ##
        results$elapsed_seconds_short_run <-  results$elapsed_seconds_long_run <-  NA_real_
        results$fixed_cost_seconds <-  results$seconds_per_iteration <-  NA_real_
        results$two_run_timing_flag <-  NA_character_
        results$iteration_multiplier_long_over_short <-  results$n_iter / results$n_iter_short_run
        results$divergences_short_run <-  NA_real_
        results$randomize_tau <-  NA
        results$n_nuisance_to_track <-  NA_real_   ## nuisance trace kept during the timed calls: 0 = none
        results$mplus_iterations_per_chain <-  results$mplus_iterations_per_chain_short_run <-  NA_real_
        results$mplus_PPPP_run <-  NA
        ##
        results$timing_scope <-  NA_character_
        results$execution_process_mode <-  NA_character_
        results$execution_process_pid <-  NA_integer_
        results$execution_process_parent_pid <-  NA_integer_
        results$stan_runtime_receipt_status <-  NA_character_
        results$stan_runtime_tbb_path_after_evaluator_initialisation <-  NA_character_
        results$stan_runtime_tbb_md5_after_evaluator_initialisation <-  NA_character_
        results$stan_runtime_tbb_path_after_timed_call <-  NA_character_
        results$stan_runtime_tbb_md5_after_timed_call <-  NA_character_
        ##
        results$status <-  "pending"
        ##
        results$error <-  NA_character_
        ##
        results$chain_seeds <-  lapply( X = seq_len(length.out = nrow(x = grid)),
                                        FUN = function(row_index) {

                case <-  grid[row_index, ]
                ##
                return(if (case$algorithm == "MD_BayesMVP_multi_process") case$base_seed + seq_len(length.out = case$n_chains) else case$seed)

        })
        results$resume_key <- NA_character_
        results$reused_completed_run <- FALSE
        results$measurement_metadata_file <- NA_character_
        cache_enabled <- !is.null(checkpoint_file) && !is.null(runtime$resume_metadata)
        if (cache_enabled) {
            resume_helpers <- new.env(parent = environment())
            sys.source(file.path(settings$algorithm_study_dir, "paper_1_chunking_and_parallel_scalability",
                                 "R_fns_alg_paper_1_resume.R"), envir = resume_helpers)
            cache_dir <- file.path(dirname(checkpoint_file), "cache", "cases")
            cache_index <- resume_helpers$fn_paper1_resume_cache_index(cache_dir = cache_dir,
                                                                         settings = settings)
        }
        ## Create every named file at the start so unrun algorithms/N values remain visibly pending.
        if (!is.null(x = checkpoint_file)) {

            fn_paper1_save_results_by_algorithm_and_N( results = results,
                                                       output_dir = dirname(path = checkpoint_file),
                                                       algorithm = NULL, N = NULL)

        }
        ## Every unique configuration/repeat occurs once, regardless of how many figures use it.
        for (row_index in seq_len(length.out = nrow(x = grid))) {

            case <-  grid[row_index, , drop = FALSE]
            ##
            dataset_index <-  match(x = case$N, table = simulated_data$N_vec)
            ##
            y <-  simulated_data$y_list[[dataset_index]]
            ##
            results$dataset_md5[row_index] <-  simulated_data$benchmark_data$dataset_md5[[as.character(x = case$N)]]
            if (cache_enabled) {
                key <- resume_helpers$fn_paper1_resume_key(case, settings, results$dataset_md5[row_index], runtime$resume_metadata)
                results$resume_key[row_index] <- key
                case$resume_key <- key
                cached <- resume_helpers$fn_paper1_cache_read(cache_dir = cache_dir,
                                                               key = key,
                                                               index = cache_index,
                                                               preserve_historical_stan_runs = settings$preserve_historical_stan_runs)
                if (!is.null(cached)) {
                    measured_columns <- intersect(names(results), setdiff(names(cached), c(names(grid), "dataset_md5", "resume_key", "reused_completed_run")))
                    results[row_index, measured_columns] <- cached[1, measured_columns, drop = FALSE]
                    ## Reuse the raw measurements, applying the currently selected reporting estimator and flag tolerance.
                    if (timing_settings$timing_method == "two_run_difference") {
                        short <- results$elapsed_seconds_short_run[row_index]
                        long <- results$elapsed_seconds_long_run[row_index]
                        overhead <- identical(expected_estimator[row_index], "overhead_subtraction")
                        actual_iter <- results$mplus_iterations_per_chain[row_index]
                        if (!is.finite(actual_iter) || actual_iter < 1) actual_iter <- case$n_iter
                        ## Match fresh-run normalization, including saved pairs
                        ## from before the one-chain FBITERATIONS count fix.
                        if (overhead && identical(case$mplus_iteration_mode, "FBITERATIONS") && case$n_chains == 1) {
                            actual_iter <- max(actual_iter, 200)
                            results$mplus_iterations_per_chain[row_index] <- actual_iter
                        }
                        per_iter <- (long - short) / if (overhead) actual_iter else (case$n_iter - case$n_iter_short_run)
                        fixed <- if (overhead) short else short - case$n_iter_short_run * per_iter
                        results$elapsed_seconds[row_index] <- if (overhead) per_iter * case$n_iter else long - fixed
                        results$fixed_cost_seconds[row_index] <- fixed
                        results$seconds_per_iteration[row_index] <- per_iter
                        results$two_run_timing_flag[row_index] <- if (long <= short) "long_run_not_slower_than_short_run" else
                            if (fixed < -timing_settings$tolerance_fraction * long) "fixed_cost_negative_beyond_tolerance" else "ok"
                    }
                    results$timing_estimator[row_index] <- expected_estimator[row_index]
                    results$reused_completed_run[row_index] <- TRUE
                    message(paste0("\033[36mSkipping completed case ", row_index, "/", nrow(grid), ": ",
                                   case$algorithm, ", N = ", case$N, ", chains = ", case$n_chains, "\033[0m"))
                    ##
                    ## ---- Refresh the saved results once the last case of this algorithm/N has been restored from the cache:
                    ## (otherwise its results_by_algorithm_and_N file keeps saying "pending" until the whole invocation ends)
                    ##
                    rows_in_same_algorithm_and_N_group <-  which(grid$algorithm == case$algorithm & grid$N == case$N)
                    ##
                    if (!is.null(x = checkpoint_file) && row_index == max(rows_in_same_algorithm_and_N_group)) {

                        saveRDS(object = results, file = checkpoint_file)
                        ##
                        fn_paper1_save_results_by_algorithm_and_N( results = results,
                                                                   output_dir = dirname(path = checkpoint_file),
                                                                   algorithm = case$algorithm,
                                                                   N = case$N)

                    }
                    ##
                    next
                }
            }
            ##
            if (isTRUE(x = cache_only)) {

                stop(paste0("Cache-only assembly found no completed result for ", case$algorithm,
                            ", N = ", case$N, ", chains = ", case$n_chains, ", repeat = ", case$run,
                            ". No sampling was started in this session. Source the runner again to complete it in a child process."),
                     call. = FALSE)

            }
            ##
            stan_partition_label <-  if (!is.na(x = results$stan_chunk_size[row_index])) {

                paste0(", chunk_size=", results$stan_chunk_size[row_index])

            } else ""
            ##
            # message(paste0((  "Case %d/%d: %s, N=%d, chunks=%d, chains=%d, WCP=%d, threads=%d, repeat=%d%s\n",
            #              row_index,
            #              nrow(x = grid),
            #              case$algorithm,
            #              case$N,
            #              
            #              case$num_chunks,
            #              
            #              case$n_chains,
            #              
            #              case$threads_per_chain,
            #              
            #              case$n_threads,
            #              
            #              case$run,
            #              
            #              stan_partition_label))
            message(paste0("Case = ", row_index, "/", nrow(x = grid), "|",
                           " Software/algorithm = ", case$algorithm,"|",
                           " Backend = ", results$execution_backend[row_index], "|",
                           " N = ", case$N, "|",
                           " num_chunks = ", case$num_chunks, "|",
                           " n_chains = ", case$n_chains, "|",
                           " threads_per_chain = ", case$threads_per_chain, "|",
                           " n_threads = ", case$n_threads, "|",
                           " run # = ", case$run, "|",
                           " stan_partition_label = ", stan_partition_label))
            ##
            output <-  tryCatch( expr = {

                ## Raw run times are checked inside; a corrected two-run time may be <= 0 and is then flagged, not dropped.
                isolate_Stan_case <-  fn_paper1_should_isolate_Stan_case( case        = case,
                                                                          settings    = settings,
                                                                          case_runner = case_runner)
                value <-  if (isTRUE(x = isolate_Stan_case)) {

                    fn_paper1_time_Stan_case_in_fresh_R_process( case            = case,
                                                                 y               = y,
                                                                 settings        = settings,
                                                                 runtime         = runtime,
                                                                 timing_settings = timing_settings)

                } else {

                    fn_paper1_time_case_with_timing_method( case            = case,
                                                            y               = y,
                                                            settings        = settings,
                                                            runtime         = runtime,
                                                            case_runner     = case_runner,
                                                            timing_settings = timing_settings)

                }
                if (is.null(x = value$execution_process_mode)) {

                    value$execution_process_mode <-  if (identical(x = case_runner, y = fn_paper1_run_case)) {

                        "current_R_process"

                    } else "custom_case_runner"
                    value$execution_process_pid <-  Sys.getpid()
                    value$execution_process_parent_pid <-  NA_integer_

                }
                ##
                if (length(x = value$elapsed_seconds) != 1L || !is.finite(x = value$elapsed_seconds)) {

                    stop("The case runner returned an invalid elapsed time.")

                }
                ##
                value

            },
                                 error = function(error) error)
            ##
            if (inherits(x = output, what = "error")) {

                results$status[row_index] <-  "failed"
                ##
                results$error[row_index] <-  conditionMessage(c = output)
                ##
                if (!is.null(x = checkpoint_file)) {

                    saveRDS(object = results, file = checkpoint_file)
                    ##
                    fn_paper1_save_results_by_algorithm_and_N( results = results,
                                                               output_dir = dirname(path = checkpoint_file),
                                                               algorithm = case$algorithm,
                                                               N = case$N)

                }
                ##
                stop("Case ", case$case_id, " failed: ", conditionMessage(c = output), call. = FALSE)

            }
            ##
            results$elapsed_seconds[row_index] <-  output$elapsed_seconds
            ##
            results$divergences[row_index] <-  output$divergences
            ##
            results$mean_L[row_index] <-  output$mean_L
            ##
            results$elapsed_seconds_short_run[row_index] <-  output$elapsed_seconds_short_run
            results$elapsed_seconds_long_run[row_index] <-  output$elapsed_seconds_long_run
            results$fixed_cost_seconds[row_index] <-  output$fixed_cost_seconds
            results$seconds_per_iteration[row_index] <-  output$seconds_per_iteration
            results$two_run_timing_flag[row_index] <-  output$two_run_timing_flag
            results$timing_estimator[row_index] <-  output$timing_estimator
            results$divergences_short_run[row_index] <-  output$divergences_short_run
            results$randomize_tau[row_index] <-  fn_paper1_runner_value_or_missing(output, "randomize_tau", NA)
            results$n_nuisance_to_track[row_index] <-  fn_paper1_runner_value_or_missing(output, "n_nuisance_to_track", NA_real_)
            results$mplus_iterations_per_chain[row_index] <-  fn_paper1_runner_value_or_missing(output, "mplus_iterations_per_chain", NA_real_)
            results$mplus_iterations_per_chain_short_run[row_index] <-  output$mplus_iterations_per_chain_short_run
            results$mplus_PPPP_run[row_index] <-  fn_paper1_runner_value_or_missing(output, "mplus_PPPP_run", NA)
            ##
            results$timing_scope[row_index] <-  output$timing_scope
            results$execution_process_mode[row_index] <-  output$execution_process_mode
            results$execution_process_pid[row_index] <-  output$execution_process_pid
            results$execution_process_parent_pid[row_index] <-  output$execution_process_parent_pid
            results$stan_runtime_receipt_status[row_index] <-  fn_paper1_runner_value_or_missing(output, "stan_runtime_receipt_status", NA_character_)
            results$stan_runtime_tbb_path_after_evaluator_initialisation[row_index] <-  fn_paper1_runner_value_or_missing(output, "stan_runtime_tbb_path_after_evaluator_initialisation", NA_character_)
            results$stan_runtime_tbb_md5_after_evaluator_initialisation[row_index] <-  fn_paper1_runner_value_or_missing(output, "stan_runtime_tbb_md5_after_evaluator_initialisation", NA_character_)
            results$stan_runtime_tbb_path_after_timed_call[row_index] <-  fn_paper1_runner_value_or_missing(output, "stan_runtime_tbb_path_after_timed_call", NA_character_)
            results$stan_runtime_tbb_md5_after_timed_call[row_index] <-  fn_paper1_runner_value_or_missing(output, "stan_runtime_tbb_md5_after_timed_call", NA_character_)
            ##
            ## ---- Release this case's sampler outputs now: the timings are recorded in `results`, so the fit
            ## objects (for the Stan arms: the full N x T nuisance trace per iteration per chain, ~4 GB per 180-chain case
            ## at N = 50,000, plus the main trace and the constrained-parameter copies) are garbage. R only collects garbage
            ## when its heap threshold is reached, so without this the child process sits on tens of GB between cases.
            ## This runs outside every timer (the next case starts with its own untimed warm-up call).
            ##
            rm(output)
            invisible(x = gc(full = TRUE))
            ##
            results$status[row_index] <-  "completed"
            if (cache_enabled) {
                results$measurement_metadata_file[row_index] <- if (is.null(runtime$resume_metadata_file)) NA_character_ else runtime$resume_metadata_file
                resume_helpers$fn_paper1_cache_write(cache_dir = cache_dir,
                                                      key = results$resume_key[row_index],
                                                      row = results[row_index, , drop = FALSE],
                                                      metadata = list(source_metadata_file = results$measurement_metadata_file[row_index]),
                                                      preserve_historical_stan_runs = settings$preserve_historical_stan_runs)
            }
            ##
            if (!is.null(x = checkpoint_file)) {

                saveRDS(object = results, file = checkpoint_file)
                ##
                fn_paper1_save_results_by_algorithm_and_N( results = results,
                                                           output_dir = dirname(path = checkpoint_file),
                                                           algorithm = case$algorithm,
                                                           N = case$N)

            }

        }
        ##
        ## An entirely cached invocation still refreshes the current-grid index and named views in the same folder.
        if (!is.null(checkpoint_file)) {
            saveRDS(results, checkpoint_file)
            fn_paper1_save_results_by_algorithm_and_N(results, dirname(checkpoint_file), algorithm = NULL, N = NULL)
        }
        return(results)

}
##
## ---- Both analyses are views of that one results table --------------------------------------------------------------------------------
##
fn_paper1_adjusted_scaling <-  function( configurations ) {

        configurations$normalisation_seconds <-  configurations$adjusted_scaling <-  NA_real_
        configurations$normalisation_configuration_id <-  NA_integer_
        if (!nrow(x = configurations)) return(configurations)
        ## The manuscript scales each implementation by its own minimum time, separately for each N and device.
        label_column <-  if ("Algorithm_label" %in% names(x = configurations)) "Algorithm_label" else "algorithm"
        group_columns <-  c("device", "N", label_column, "execution_backend", "n_iter", "dataset_md5", "timing_scope", "timing_estimator")
        groups <-  unique(x = configurations[, group_columns, drop = FALSE])
        ##
        for (group_index in seq_len(length.out = nrow(x = groups))) {

            keep <-  rep(x = TRUE, times = nrow(x = configurations))
            for (column in group_columns) keep <-  keep & configurations[[column]] == groups[[column]][group_index]
            rows <-  which(x = keep)
            ## A flagged non-positive two-run time can never be the own-minimum-time baseline.
            positive_time_rows <-  rows[configurations$time_mean[rows] > 0]
            if (!length(x = positive_time_rows)) next
            baseline <-  positive_time_rows[which.min(x = configurations$time_mean[positive_time_rows])]
            configurations$normalisation_seconds[rows] <-  configurations$time_mean[baseline]
            configurations$normalisation_configuration_id[rows] <-  configurations$configuration_id[baseline]
            configurations$adjusted_scaling[rows] <-  configurations$chain_rate[rows] * configurations$time_mean[baseline]

        }
        ##
        return(configurations)

}
##
fn_summarise_paper1_benchmark <-  function( results ) {

        if (!"timing_estimator" %in% names(results)) {
            results$timing_estimator <-  if ("timing_method" %in% names(results))
                ifelse(results$timing_method == "two_run_difference", "iteration_difference", "single_run") else "single_run"
        }
        complete <-  results[results$status == "completed", , drop = FALSE]
        complete <-  fn_paper1_execution_columns(cases = fn_paper1_stan_partition_columns(cases = complete),
                                                   stan_backend = NULL, native_backend = NULL)
        ##
        if (!nrow(x = complete)) return(list( configurations = data.frame(),
                                              optimal_combinations = data.frame(),
                                              best_chunks = data.frame(),
                                              wcp_chunk_search = data.frame(),
                                              wcp_best_chunks = data.frame(),
                                              scaling = data.frame(),
                                              chunking = data.frame(),
                                              wcp_matched = data.frame()))
        ##
        columns <-  c("device", "algorithm", "execution_backend",
                      "N", 
                      "num_chunks", "n_threads", "n_chains", "threads_per_chain",
                      "n_iter",
                      "dataset_md5",
                      "timing_scope", "timing_estimator")
        ##
        if ("benchmark_role" %in% names(x = complete)) columns <-  c(columns, "benchmark_role")
        ##
        ## ---- One timing method per analysis:
        ##
        ## Two-run corrected times exclude the fixed per-run cost and single-run times include it, so they are never mixed.
        ##
        if ("timing_method" %in% names(x = complete)) {

            if (length(x = unique(x = complete$timing_method)) != 1L || anyNA(x = complete$timing_method)) {

                stop("These results mix timing methods (", paste(unique(x = complete$timing_method), collapse = ", "),
                     "); summarise each timing method separately.")

            }
            ##
            columns <-  c(columns, "timing_method")

        }
        ##
        configurations <-  unique(x = complete[, columns, drop = FALSE])
        configurations <-  fn_paper1_stan_partition_columns(cases = configurations)
        ##
        summaries <-  lapply( X = seq_len(length.out = nrow(x = configurations)),
                              FUN = function(configuration_index) {

                configuration <-  configurations[configuration_index, , drop = FALSE]
                ##
                keep <-  rep(x = TRUE, times = nrow(x = complete))
                ##
                for (column in columns) keep <-  keep & complete[[column]] == configuration[[column]]
                ##
                rows <-  complete[keep, , drop = FALSE]
                ##
                configuration$n_repeats <-  nrow(x = rows)
                ##
                ## ---- Repeated runs of one configuration are summarised by their arithmetic mean (time_mean and the *_mean columns);
                ##      the medians are kept alongside for reference only, and no table, figure or rate uses them:
                ##
                configuration$time_mean <-  mean(x = rows$elapsed_seconds)
                configuration$time_median <-  stats::median(x = rows$elapsed_seconds)
                ##
                configuration$time_sd <-  stats::sd(x = rows$elapsed_seconds)
                ##
                ## ---- Two-run timing detail per configuration; NA for single-run results:
                ##
                ## A configuration whose median corrected time is not positive keeps its numbers and flags, but has no rate.
                ##
                two_run_rows <-  "two_run_timing_flag" %in% names(x = rows) && any(!is.na(x = rows$two_run_timing_flag))
                configuration$time_short_run_mean <-  if (two_run_rows) mean(x = rows$elapsed_seconds_short_run) else NA_real_
                configuration$time_long_run_mean <-  if (two_run_rows) mean(x = rows$elapsed_seconds_long_run) else NA_real_
                configuration$fixed_cost_mean <-  if (two_run_rows) mean(x = rows$fixed_cost_seconds) else NA_real_
                configuration$seconds_per_iteration_mean <-  if (two_run_rows) mean(x = rows$seconds_per_iteration) else NA_real_
                configuration$time_short_run_median <-  if (two_run_rows) stats::median(x = rows$elapsed_seconds_short_run) else NA_real_
                configuration$time_long_run_median <-  if (two_run_rows) stats::median(x = rows$elapsed_seconds_long_run) else NA_real_
                configuration$fixed_cost_median <-  if (two_run_rows) stats::median(x = rows$fixed_cost_seconds) else NA_real_
                configuration$seconds_per_iteration_median <-  if (two_run_rows) stats::median(x = rows$seconds_per_iteration) else NA_real_
                configuration$n_repeats_two_run_flagged <-  if (two_run_rows) sum(rows$two_run_timing_flag != "ok") else NA_real_
                configuration$time_mean_not_positive <-  configuration$time_mean <= 0
                configuration$time_median_not_positive <-  configuration$time_median <= 0
                ##
                configuration$chain_rate <-  if (configuration$time_mean_not_positive) NA_real_ else configuration$n_chains / configuration$time_mean
                ## Iteration throughput is retained for within-sampler comparisons only.
                configuration$total_iter_per_sec <-  if (configuration$time_mean_not_positive) NA_real_ else
                    configuration$n_chains * configuration$n_iter / configuration$time_mean
                ##
                configuration$max_divergences <-  if (all(is.na(x = rows$divergences))) NA_real_ else max(rows$divergences, na.rm = TRUE)
                ##
                configuration$case_ids <-  list(rows$case_id)
                ##
                return(configuration)

        })
        ##
        configurations <-  do.call(what = rbind, args = summaries)
        configurations$configuration_id <-  seq_len(length.out = nrow(x = configurations))
        ##
        ## First select chunks at fixed chains AND WCP threads, separately within every algorithm.
        best <-  configurations[order(-configurations$chain_rate, configurations$num_chunks), , drop = FALSE]
        ##
        group_columns <-  c("device", "algorithm", "execution_backend", "N", "n_threads", "n_chains", "threads_per_chain",
                             "n_iter", "dataset_md5", "timing_scope", "timing_estimator")
        if ("benchmark_role" %in% names(x = configurations)) group_columns <-  c(group_columns, "benchmark_role")
        ##
        best <-  best[!duplicated(x = best[, group_columns, drop = FALSE]), , drop = FALSE]
        ##
        wcp_chunk_search <-  configurations[configurations$algorithm %in% c("MD_BayesMVP_WCP", "AD_Stan_WCP"), , drop = FALSE]
        ##
        if ("benchmark_role" %in% names(x = wcp_chunk_search)) {

            wcp_chunk_search <-  wcp_chunk_search[wcp_chunk_search$benchmark_role == "main_scaling", , drop = FALSE]

        }
        ##
        wcp_chunk_search$selected_best_chunks <-  wcp_chunk_search$configuration_id %in% best$configuration_id
        wcp_chunk_search$relative_chunk_throughput <-  rep(x = NA_real_, times = nrow(x = wcp_chunk_search))
        ##
        for (row_index in seq_len(length.out = nrow(x = wcp_chunk_search))) {

            same_allocation <-  rep(x = TRUE, times = nrow(x = best))
            for (column in group_columns) {

                same_allocation <-  same_allocation & best[[column]] == wcp_chunk_search[[column]][row_index]

            }
            ##
            wcp_chunk_search$relative_chunk_throughput[row_index] <-
                wcp_chunk_search$chain_rate[row_index] / best$chain_rate[same_allocation]

        }
        ##
        wcp_best_chunks <-  wcp_chunk_search[wcp_chunk_search$selected_best_chunks, , drop = FALSE]
        ##
        ## Optional second-stage summary over WCP choices; keep each algorithm arm separate.
        ## The conditional WCP chunk optima above, not this collapsed view, feed the scaling comparison.
        optimal_combinations <-  best
        optimal_combinations$implementation <-  optimal_combinations$algorithm
        ##
        optimal_combinations <-  optimal_combinations[order(-optimal_combinations$chain_rate,
                                                              optimal_combinations$num_chunks,
                                                              optimal_combinations$threads_per_chain), , drop = FALSE]
        ##
        optimum_groups <-  setdiff(x = group_columns, y = c("n_threads", "threads_per_chain"))
        optimal_combinations <-  optimal_combinations[!duplicated(x = optimal_combinations[, optimum_groups, drop = FALSE]), , drop = FALSE]
        ##
        scaling <-  best
        ##
        if ("benchmark_role" %in% names(x = scaling)) {

            scaling <-  scaling[scaling$benchmark_role == "main_scaling", , drop = FALSE]
            scaling <-  scaling[order(-scaling$chain_rate, scaling$num_chunks, scaling$threads_per_chain), , drop = FALSE]
            ##
            scaling_groups <-  c("device", "algorithm", "N", "n_threads", "n_iter", "dataset_md5", "timing_scope", "timing_estimator")
            scaling <-  scaling[!duplicated(x = scaling[, scaling_groups, drop = FALSE]), , drop = FALSE]

        }
        ##
        scaling$baseline_threads <-  scaling$throughput_ratio <-  NA_real_
        ##
        for (row_index in seq_len(length.out = nrow(x = scaling))) {

            same <-  scaling$device == scaling$device[row_index] & scaling$algorithm == scaling$algorithm[row_index] &
                scaling$N == scaling$N[row_index] & scaling$n_iter == scaling$n_iter[row_index] &
                scaling$dataset_md5 == scaling$dataset_md5[row_index] & scaling$timing_scope == scaling$timing_scope[row_index] &
                scaling$timing_estimator == scaling$timing_estimator[row_index]
            ##
            ## Retain the old fixed-chain interpretation only for results saved before explicit study roles were introduced.
            if (!("benchmark_role" %in% names(x = scaling)) && scaling$algorithm[row_index] == "MD_BayesMVP_WCP") {

                same <-  same & scaling$n_chains == scaling$n_chains[row_index]

            }
            ##
            baseline_index <-  which(x = same)[which.min(x = scaling$n_threads[same])]
            ##
            scaling$baseline_threads[row_index] <-  scaling$n_threads[baseline_index]
            ##
            scaling$throughput_ratio[row_index] <-  scaling$chain_rate[row_index] / scaling$chain_rate[baseline_index]

        }
        ##
        ## Direct chunking comparison at unchanged chain/thread counts, wherever a one-chunk reference exists.
        ##
        chunking <-  configurations
        ##
        chunking$one_chunk_seconds <-  chunking$chunking_speedup <-  NA_real_
        ##
        chunking$baseline_case_ids <-  rep(x = list(integer(length = 0)), times = nrow(x = chunking))
        ##
        for (row_index in seq_len(length.out = nrow(x = chunking))) {

            same <-  configurations$num_chunks == 1
            ##
            for (column in group_columns) same <-  same & configurations[[column]] == chunking[[column]][row_index]
            ##
            baseline <-  which(x = same)
            ##
            if (length(x = baseline) == 1L) {

                chunking$one_chunk_seconds[row_index] <-  configurations$time_mean[baseline]
                ##
                chunking$chunking_speedup[row_index] <-  configurations$time_mean[baseline] / chunking$time_mean[row_index]
                ##
                chunking$baseline_case_ids[[row_index]] <-  configurations$case_ids[[baseline]]

            }

        }
        ##
        ## WCP versus serial at the SAME N, chunk count, chain count, iterations and data; seeds match by repeat.
        ##
        wcp <-  configurations[configurations$algorithm == "MD_BayesMVP_WCP", , drop = FALSE]
        ##
        wcp$serial_seconds <-  wcp$serial_measured_seconds <-  wcp$serial_n_iter <-
            wcp$wcp_speedup <-  wcp$wcp_efficiency <-  rep(x = NA_real_, times = nrow(x = wcp))
        ##
        wcp$baseline_case_ids <-  rep(x = list(integer(length = 0)), times = nrow(x = wcp))
        ##
        for (row_index in seq_len(length.out = nrow(x = wcp))) {

            same <-  configurations$algorithm == "MD_BayesMVP" & configurations$threads_per_chain == 1
            ##
            for (column in c("device", "N", "num_chunks", "n_chains", "dataset_md5", "timing_scope", "timing_estimator")) {

                same <-  same & configurations[[column]] == wcp[[column]][row_index]

            }
            ##
            baseline <-  which(x = same)
            ##
            if (length(x = baseline) == 1L && configurations$time_mean[baseline] > 0 && wcp$time_mean[row_index] > 0) {

                same_iterations <-  configurations$n_iter[baseline] == wcp$n_iter[row_index]
                corrected_pair <-  "timing_method" %in% names(configurations) &&
                    identical(as.character(configurations$timing_method[baseline]), "two_run_difference") &&
                    identical(as.character(wcp$timing_method[row_index]), "two_run_difference")
                if (!same_iterations && !corrected_pair) next
                ## Scale the serial measurement to the WCP budget; unequal budgets require corrected two-run times.
                wcp$serial_n_iter[row_index] <-  configurations$n_iter[baseline]
                wcp$serial_measured_seconds[row_index] <-  configurations$time_mean[baseline]
                wcp$serial_seconds[row_index] <-  configurations$time_mean[baseline] *
                    wcp$n_iter[row_index] / configurations$n_iter[baseline]
                ##
                wcp$wcp_speedup[row_index] <-  wcp$serial_seconds[row_index] / wcp$time_mean[row_index]
                ##
                wcp$wcp_efficiency[row_index] <-  wcp$wcp_speedup[row_index] / wcp$threads_per_chain[row_index]
                ##
                wcp$baseline_case_ids[[row_index]] <-  configurations$case_ids[[baseline]]

            }

        }
        ##
        return(list(configurations = configurations,
                    optimal_combinations = optimal_combinations, 
                    best_chunks = best,
                    wcp_chunk_search = wcp_chunk_search,
                    wcp_best_chunks = wcp_best_chunks,
                    scaling = scaling, 
                    chunking = chunking,
                    wcp_matched = wcp))

}
##
## ---- Public entry point: all arguments and study settings are explicit ---------------------------------------------------------------
##
fn_run_paper1_benchmark <-  function( settings,
                                      dry_run,
                                      Stan_model_obj,
                                      cache_only = FALSE
) {

        if (!identical(settings$package_stack, "NicoStan_BayesMVP")) stop("Set package_stack explicitly to NicoStan_BayesMVP in the runner.")
        if (!identical(settings$stan_reduce_sum_type, "static")) stop("Paper 1 currently supports stan_reduce_sum_type = static only.")
        if (!is.logical(x = settings$run_each_Stan_case_in_fresh_R_process) ||
            length(x = settings$run_each_Stan_case_in_fresh_R_process) != 1L ||
            is.na(x = settings$run_each_Stan_case_in_fresh_R_process)) {

            stop("Set run_each_Stan_case_in_fresh_R_process explicitly to TRUE or FALSE in the Paper 1 runner.")

        }
        settings$stan_backend <-  fn_paper1_stan_backend(settings = settings)
        if (settings$stan_backend == "NicoStan") {

            settings$stan_via_NicoStan <-  fn_paper1_stan_via_NicoStan_settings(settings = settings)
            if (!is.null(x = Stan_model_obj)) stop("Stan_model_obj supplies CmdStanR objects; set it to NULL for the NicoStan backend.")

        } else settings$stan_via_NicoStan <-  NULL
        ##
        grid <-  fn_paper1_benchmark_grid(settings = settings)
        ##
        counts <-  stats::aggregate( x = list(timed_calls = rep(x = 1L, times = nrow(x = grid))),
                                     by = grid[c("benchmark_role", "algorithm", "N")],
                                     FUN = sum)
        ##
        print(x = counts, row.names = FALSE)
        ##
        message(paste0("Total: ", nrow(x = grid), " configuration/repeat pairs; ",
                       nrow(x = grid) * ifelse(settings$timing_method == "two_run_difference", 2, 1), " sampler calls."))
        message(paste0("Stan execution backend: ", settings$stan_backend,
                       if (settings$stan_backend == "NicoStan") paste0("; step size = ", settings$stan_via_NicoStan$step_size,
                           "; fixed L = ", settings$stan_via_NicoStan$fixed_L) else ""))
        ##
        ## ---- Timing method, the iteration counts that will reach each sampler, and the configured run-time budget:
        ##
        timing_settings <-  fn_paper1_timing_settings(settings = settings, algorithms = settings$algorithms)
        timing_iteration_table <-  unique(x = grid[, c("algorithm", "N", "n_iter_short_run", "n_iter"), drop = FALSE])
        names(x = timing_iteration_table)[names(x = timing_iteration_table) == "n_iter"] <-  "n_iter_long_run"
        message(paste0("\033[36m", "Timing method: ", timing_settings$timing_method,
                       if (timing_settings$timing_method == "two_run_difference")
                           " (short run then long run per configuration; reported time = T_long - fixed cost)" else
                           " (one timed run per configuration; fixed per-run cost included)",
                       "\033[0m"))
        print(x = timing_iteration_table, row.names = FALSE)
        message(paste0("Repeats: BayesMVP/Stan = ", settings$n_runs, "; Mplus = ", settings$mplus_n_runs,
                       "; Mplus long mode(s) = ", paste(settings$mplus_iteration_mode, collapse = " + "),
                       "; Mplus short mode = ", settings$mplus_iteration_mode_short_run,
                       "; Mplus short role = ", settings$mplus_short_run_role,
                       "; minimum for BITERATIONS = ", settings$mplus_biterations_minimum,
                       "; BCONVERGENCE = ", settings$mplus_bconvergence,
                       "; save draws = ", settings$mplus_save_draws, "; require saved draws = ", settings$mplus_verify_saved_draws))
        message(paste0("BayesMVP: L = ", settings$bayesmvp$L_main, "; epsilon = ", settings$bayesmvp$eps_main,
                       "; randomize_tau = ", settings$bayesmvp$randomize_tau))
        ##
        ##
        if (isTRUE(x = dry_run)) return(invisible(x = list(grid = grid, counts = counts, settings = settings)))
        ##
        if (any(grid$n_threads > parallel::detectCores())) stop("Selected thread counts exceed this machine's detected logical CPUs.")
        ##
        if (any(grepl(pattern = "^Mplus", x = grid$algorithm)) &&
            (!is.function(x = settings$mplus$runner) || !is.function(x = settings$mplus$verifier) || is.null(x = settings$mplus$settings))) {

            stop("Supply the original compatible Mplus runner/verifier/settings before selecting Mplus; no replacement priors are assumed.")

        }
        ## Check the Mplus executable before starting any of the selected study arms.
        if (any(grepl(pattern = "^Mplus", x = grid$algorithm))) {

            if (!requireNamespace(package = "MplusAutomation", quietly = TRUE)) stop("MplusAutomation is required for Mplus.")
            mplus_command <-  settings$mplus$settings$Mplus_command
            if (is.null(x = mplus_command)) mplus_command <-  MplusAutomation::detectMplus()
            if (!nzchar(x = Sys.which(names = mplus_command)) && !file.exists(file = mplus_command)) {

                stop("Mplus executable not found. Set paper1_settings$mplus$settings$Mplus_command to its path.")

            }
            settings$mplus$settings$Mplus_command <-  mplus_command

        }
        ## Use the new shared sampler and specialised extension for every newly measured native/Stan case.
        if (any(grepl(pattern = "^(MD_|AD_)", x = grid$algorithm))) fn_paper1_load_split_packages(settings = settings)
        if (settings$stan_backend == "NicoStan" && length(x = grep(pattern = "^AD_", x = grid$algorithm))) {

            invisible(x = fn_paper1_stan_HMC_args(configuration = settings$stan_via_NicoStan))

        }
        ##
        ## The stable device directory is reused; completed cases are matched by work/data/build signatures below.
        ##
        data <-  fn_paper1_COVID_data(algorithm_study_dir = settings$algorithm_study_dir, N_vec = settings$N_vec, seed = settings$data_seed)
        ##
        stan_algorithms <-  unique(x = grid$algorithm[grepl(pattern = "^AD_", x = grid$algorithm)])
        ##
        models <-  if (!length(x = stan_algorithms)) list() else if (settings$stan_backend == "NicoStan") {

            fn_paper1_compile_stan_via_NicoStan(settings = settings, algorithms = stan_algorithms)

        } else if (is.null(x = Stan_model_obj)) {

            fn_paper1_compile_stan(settings = settings, algorithms = stan_algorithms)

        } else Stan_model_obj
        ##
        if (any(!stan_algorithms %in% names(x = models))) stop("Missing a selected Stan model.")
        ##
        for (algorithm in if (settings$stan_backend == "cmdstanr")
            intersect(x = stan_algorithms, y = c("AD_Stan_WCP", "AD_Stan_tape_chunked")) else character(length = 0)) {

            expected_stan_file <-  file.path(settings$algorithm_study_dir, "paper_1_chunking_and_parallel_scalability", "stan_models",
                                             "LC_MVP_bin_PartialLog_v5_reduce_sum_static.stan")
            supplied_stan_file <-  models[[algorithm]]$stan_file()
            ## Reject previously supplied dynamic models as well as compiling the distinct static source by default.
            if (is.null(x = supplied_stan_file) ||
                !identical(normalizePath(path = supplied_stan_file), normalizePath(path = expected_stan_file))) {

                stop("Supply the Paper 1 reduce_sum_static model for ", algorithm, "; set Stan_model_obj = NULL to compile it.")

            }
            ## A cached object compiled before the data-field rename still expects the old JSON field.
            supplied_executable <-  models[[algorithm]]$exe_file()
            if (is.null(x = supplied_executable) || !file.exists(file = supplied_executable) ||
                file.info(supplied_executable)$mtime < file.info(expected_stan_file)$mtime) {

                stop("The supplied Stan executable predates the current chunk_size model. Recompile it before running ", algorithm, ".")

            }
            ##
            if (!isTRUE(x = cmdstanr:::model_compile_info(models[[algorithm]]$exe_file())[["STAN_THREADS"]])) {

                stop("Supplied reduce_sum model is not compiled with STAN_THREADS=true.")

            }

        }
        ##
        runtime <-  list( stan_models = models,
                          stan_data = setNames(object = lapply(X = data$y_list, FUN = fn_paper1_stan_data), nm = as.character(x = data$N_vec)))
        ##
        metadata_environment <-  new.env(parent = globalenv())
        ##
        sys.source( file = file.path(settings$algorithm_study_dir, "0_utilities", "shared_functions", "R_fns_benchmark_timing_results.R"),
                    envir = metadata_environment)
        ##
        saved_settings <-  settings
        ##
        saved_settings$mplus$runner <-  if (is.function(x = settings$mplus$runner)) paste( deparse(expr = body(fun = settings$mplus$runner)),
                                                                                           collapse = "\n") else NULL
        saved_settings$mplus$verifier <-  if (is.function(x = settings$mplus$verifier)) paste( deparse(expr = body(fun = settings$mplus$verifier)),
                                                                                               collapse = "\n") else NULL
        ##
        stan_runtime_receipt <-  if (settings$stan_backend == "NicoStan" && length(x = stan_algorithms)) {
            fn_paper1_stan_runtime_receipt(stage = "model_compilation", require_single_tbb = FALSE)
        } else NULL
        stan_load_policy <-  if (settings$stan_backend == "NicoStan" && length(x = stan_algorithms)) {
            fn_paper1_stan_load_policy_fingerprint()
        } else NULL
        ##
        metadata <-  list( settings = saved_settings,
                           data = data$benchmark_data,
                           machine = as.list(Sys.info()[c("nodename", "sysname", "release", "machine")]),
                           ##
                           ## The timing method, every short/long iteration setting and the budget estimate:
                           ##
                           timing = list( timing_method = timing_settings$timing_method,
                                          mplus_short_run_role = settings$mplus_short_run_role,
                                          timing_estimators = unique(grid[, c("algorithm", "timing_estimator")]),
                                          negative_fixed_cost_tolerance_fraction = timing_settings$tolerance_fraction,
                                          iteration_settings = settings[c(unlist(x = fn_paper1_iteration_setting_names(), use.names = FALSE), "bayesmvp_wcp_iterations", "stan_wcp_iterations", "mplus_standard_iterations")],
                                          iterations_by_algorithm_and_N = timing_iteration_table,
                                          stan_fixed_L = if (is.null(x = settings$stan_via_NicoStan)) NA_real_ else settings$stan_via_NicoStan$fixed_L,
                                          bayesmvp_fixed_L = settings$bayesmvp$L_main),
                           runner_functions = lapply( X = c("fn_paper1_benchmark_grid", "fn_paper1_run_case", "fn_execute_paper1_grid",
                                                            "fn_paper1_save_results_by_algorithm_and_N", "fn_summarise_paper1_benchmark",
                                                            "fn_paper1_COVID_data", "fn_paper1_stan_data", "fn_paper1_stan_partition_columns",
                                                            "fn_paper1_bayesmvp_initial_state", "fn_paper1_compile_stan",
                                                            "fn_paper1_stan_backend", "fn_paper1_execution_columns",
                                                            "fn_paper1_distinct_tbb_paths_from_maps", "fn_paper1_stan_runtime_receipt",
                                                            "fn_paper1_stan_load_policy_fingerprint",
                                                            "fn_paper1_stan_via_NicoStan_settings", "fn_paper1_compile_stan_via_NicoStan",
                                                            "fn_paper1_run_stan_via_NicoStan", "fn_paper1_stan_HMC_args",
                                                            "fn_paper1_check_fixed_trajectory", "fn_paper1_native_backend",
                                                            "fn_paper1_load_split_packages", "fn_paper1_adjusted_scaling",
                                                            "fn_paper1_algorithm_family", "fn_paper1_iteration_setting_names",
                                                            "fn_paper1_iterations_for_case", "fn_paper1_timing_settings",
                                                            "fn_paper1_runner_value_or_missing", "fn_paper1_time_case_with_timing_method",
                                                            "fn_paper1_should_isolate_Stan_case",
                                                            "fn_paper1_time_Stan_case_in_fresh_R_process"),
                                                      FUN = function(function_name)
                                                          list( name = function_name,
                                                                definition = deparse(expr = get(x = function_name, mode = "function")))),
                           package_build = metadata_environment$fn_benchmark_package_build_metadata(package_name = "BayesMVP"),
                           NicoStan_build = if ("NicoStan" %in% loadedNamespaces())
                               metadata_environment$fn_benchmark_package_build_metadata(package_name = "NicoStan") else NULL,
                           stan_runtime_receipt = stan_runtime_receipt,
                           stan_load_policy = stan_load_policy,
                           model_builds = if (settings$stan_backend == "NicoStan") models else lapply( X = models,
                                                  FUN = function(model)
                                                      list( executable_path = model$exe_file(),
                                                            executable_md5 = unname(obj = tools::md5sum(files = model$exe_file())),
                                                            stan_source_path = model$stan_file(),
                                                            stan_source_md5 = unname(obj = tools::md5sum(files = model$stan_file())),
                                                            cmdstanr_version = as.character(x = utils::packageVersion(pkg = "cmdstanr")))))
        ##
        if (any(grepl(pattern = "^Mplus", x = grid$algorithm))) {

            mplus_source_file <-  file.path(settings$algorithm_study_dir, "paper_1_chunking_and_parallel_scalability", "R_fns_alg_paper_1_Mplus.R")
            metadata$mplus_build <-  list( source_file = mplus_source_file,
                                           source_md5 = unname(obj = tools::md5sum(files = mplus_source_file)),
                                           source_code = readLines(con = mplus_source_file, warn = FALSE),
                                           MplusAutomation_version = as.character(x = utils::packageVersion(pkg = "MplusAutomation")),
                                           command = if (is.null(x = settings$mplus$settings$Mplus_command))
                                               MplusAutomation::detectMplus() else settings$mplus$settings$Mplus_command)
            mplus_executable <- metadata$mplus_build$command
            if (!file.exists(mplus_executable)) mplus_executable <- unname(Sys.which(mplus_executable))
            metadata$mplus_build$executable_md5 <- unname(tools::md5sum(mplus_executable))

        }
        ##
        dir.create(path = settings$output_dir, recursive = TRUE, showWarnings = FALSE)
        runtime$resume_metadata <- metadata
        metadata_dir <- file.path(settings$output_dir, "cache", "metadata")
        dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
        runtime$resume_metadata_file <- file.path(metadata_dir, paste0(fn_paper1_object_md5(metadata), ".rds"))
        if (!file.exists(runtime$resume_metadata_file)) saveRDS(metadata, runtime$resume_metadata_file)
        dataset_cache_dir <- file.path(settings$output_dir, "cache", "datasets")
        dir.create(dataset_cache_dir, recursive = TRUE, showWarnings = FALSE)
        dataset_cache_file <- file.path(dataset_cache_dir, paste0(fn_paper1_object_md5(data$benchmark_data), ".rds"))
        if (!file.exists(dataset_cache_file)) saveRDS(data, dataset_cache_file)
        ##
        saveRDS(object = metadata, file = file.path(settings$output_dir, "study_metadata.rds"))
        ##
        saveRDS(object = data, file = file.path(settings$output_dir, "COVID19_datasets.rds"))
        ##
        saveRDS(object = grid, file = file.path(settings$output_dir, "experiment_grid.rds"))
        ##
        results <-  fn_execute_paper1_grid( grid = grid,
                                            simulated_data = data,
                                            settings = settings,
                                            runtime = runtime,
                                            case_runner = fn_paper1_run_case,
                                            checkpoint_file = file.path(settings$output_dir, "results.rds"),
                                            cache_only = cache_only)
        ##
        ##
        ## ---- Main analysis on the first Mplus long-run mode only; every mode stays in results.rds / results.csv, and
        ## ---- mplus_iteration_mode_comparison.csv summarises each Mplus configuration separately per mode.
        ##
        analysis_results <-  results[!(grepl("^Mplus_", results$algorithm) & !is.na(results$mplus_iteration_mode) &
                                       results$mplus_iteration_mode != settings$mplus_iteration_mode[1]), , drop = FALSE]
        analysis <-  fn_summarise_paper1_benchmark(results = analysis_results)
        ##
        saveRDS(object = analysis, file = file.path(settings$output_dir, "analysis.rds"))
        ##
        mplus_mode_results <-  results[grepl("^Mplus_", results$algorithm) & !is.na(results$mplus_iteration_mode), , drop = FALSE]
        if (length(x = unique(x = mplus_mode_results$mplus_iteration_mode)) > 1L) {

            mode_summaries <-  lapply( X = split(x = mplus_mode_results, f = mplus_mode_results$mplus_iteration_mode),
                                       FUN = function(mode_results) {

                mode_configurations <-  fn_summarise_paper1_benchmark(results = mode_results)$configurations
                mode_configurations$mplus_iteration_mode <-  unique(x = mode_results$mplus_iteration_mode)
                return(mode_configurations)

            })
            mode_comparison <-  do.call(what = rbind, args = mode_summaries)
            mode_comparison <-  mode_comparison[, vapply(X = mode_comparison, FUN = function(column) !is.list(column), FUN.VALUE = logical(length = 1)), drop = FALSE]
            utils::write.csv(x = mode_comparison, file = file.path(settings$output_dir, "mplus_iteration_mode_comparison.csv"), row.names = FALSE)
            message(paste0("\033[36mMplus FBITERATIONS vs BITERATIONS summary: ", file.path(settings$output_dir, "mplus_iteration_mode_comparison.csv"), "\033[0m"))

        }
        ## CSV views contain the same results, with list columns represented as semicolon-separated values.
        results_csv <-  results
        ##
        results_csv$chain_seeds <-  vapply(X = results$chain_seeds, FUN = paste, FUN.VALUE = character(length = 1), collapse = ";")
        ##
        utils::write.csv(x = results_csv, file = file.path(settings$output_dir, "results.csv"), row.names = FALSE)
        ##
        for (view in names(x = analysis)) {

            view_csv <-  analysis[[view]]
            ##
            for (column in names(x = view_csv)[vapply(X = view_csv, FUN = is.list, FUN.VALUE = logical(length = 1))]) {

                view_csv[[column]] <-  vapply(X = view_csv[[column]], FUN = paste, FUN.VALUE = character(length = 1), collapse = ";")

            }
            ##
            utils::write.csv(x = view_csv, file = file.path(settings$output_dir, paste0(view, ".csv")), row.names = FALSE)

        }
        ##
        return(invisible(x = list(results = results,
                                  analysis = analysis,
                                  metadata = metadata,
                                  output_dir = settings$output_dir)))

}






















