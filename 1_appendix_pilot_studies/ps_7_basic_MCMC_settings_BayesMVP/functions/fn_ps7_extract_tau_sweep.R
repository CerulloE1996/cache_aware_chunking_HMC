##
## -| --------- Extract the tau-scheme sweep from saved ps7 runs -------------------------------------------------
##
## summarize_ps7_results() rebuilds every file name from a settings grid, which means it only finds
## runs whose settings can be stated exactly. This does the opposite: it discovers whatever is on disk and
## recovers the settings FROM THE NAME. That makes it robust to a run list that drifted from the
## driver's option block, and it is the right tool for "what did I actually just produce".
## include_settings / exclude_settings select parsed filenames BEFORE opening any saved fit.
##
## Each RDS carries the full traces (about 36 MB), so files are read ONE AT A TIME, reduced to scalars,
## and dropped before the next. Reading them into a list first will exhaust memory.
##
## The parsed fields mirror the encoders in ps_7_MCMC_settings_BayesMVP_functions.R exactly. Absent
## tag = the original ChESSR behaviour, which is what "encodes to nothing" means there.
##
##
## ---- printed names of the two joint-diffusion integrators (defined ONCE here, so the attribution is changed in one place):
##
R_fn_ps7_integrator_label <-  function( diffusion_HMC_integrator ) {

        integrator_labels <-  c( kick_flow_kick = "KFK (BayesMVP/Cerullo)",
                                flow_kick_flow = "FKF (ADL)")
        return(unname(ifelse(diffusion_HMC_integrator %in% names(integrator_labels),
                             integrator_labels[diffusion_HMC_integrator],
                             diffusion_HMC_integrator)))

}
##
#' fn_ps7_parse_run_name
#' @param x  a saved run's basename.
#' @return one-row tibble of the settings the name encodes.
#' @export
fn_ps7_parse_run_name <-  function( x ) {

        ## Remove only the new suffix while parsing the existing tags; retain it in the file/configuration identity.
        run_file_name <-  x
        ##
        ## ---- Trajectory-jitter tokens, written AFTER "_af1" in this order: "_rb1", "_ts<G|number>", "_ab1".
        ##      Read and strip them FIRST, so the older end-anchored patterns below see exactly the names they were written for.
        ##      Absent = randomize_tau_burnin FALSE / tau_sampling_scale "none" / ADAM bias correction by iteration index.
        ##      "_tbJ" (EXPERIMENTAL tau_adaptation_block axis) is written AFTER "_ab1", still before "_runK";
        ##      absent = tau_adaptation_block "main" (the existing behaviour).
        ##
        {
            trajectory_jitter_suffix_pattern <-  "(_rb1)?(_ts(G|[0-9.][0-9.eE+-]*))?(_ab1)?(_tbJ)?(_run[0-9]+)$"
            trajectory_jitter_suffix <-  regmatches(x = x, m = regexpr(pattern = trajectory_jitter_suffix_pattern, text = x))
            if (length(x = trajectory_jitter_suffix) == 0) trajectory_jitter_suffix <-  ""   ## a name without "_runK"
            randomize_tau_burnin <-  grepl(pattern = "^_rb1", x = trajectory_jitter_suffix)
            tau_sampling_scale_token <-  sub(pattern = "^(_rb1)?(_ts(G|[0-9.][0-9.eE+-]*))?(_ab1)?(_tbJ)?(_run[0-9]+)$",
                                            replacement = "\\3",
                                            x = trajectory_jitter_suffix)
            tau_sampling_scale <-  if (identical(tau_sampling_scale_token, "")) "none" else
                                  if (identical(tau_sampling_scale_token, "G")) "gaussian_matched" else
                                  as.character(as.numeric(sub(pattern = "^\\.", replacement = "0.", x = tau_sampling_scale_token)))
            tau_adam_bias_correction <-  if (grepl(pattern = "_ab1(_tbJ)?_run[0-9]+$", x = x)) "performed_update_counter" else "iteration_index"
            tau_adaptation_block <-  if (grepl(pattern = "_tbJ_run[0-9]+$", x = x)) "joint" else "main"
            x <-  sub(pattern = trajectory_jitter_suffix_pattern, replacement = "\\6", x = x)
        }
        ##
        autodiff_fallback <-  grepl(pattern = "_af1_run[0-9]+$", x = x)
        x <-  sub(pattern = "_af1(_run[0-9]+)$", replacement = "\\1", x = x)
        ##
        grab <-  function( pattern, default = NA_character_ ) {
            m <-  regmatches(x = x, m = regexpr(pattern = pattern, text = x))
            if (length(x = m) == 0) return(default)
            sub(pattern = pattern, replacement = "\\1", x = m)
        }
        ##
        manual_L <-  grab("_mL([0-9.]+)_")
        ## must START with a digit or "." (R_fn_enc_num: 0.5 -> ".5"): "_mt" is also the metric-type prefix ("_mtE_"),
        ## which the old pattern "_mt([0-9.eE+-]+)_" matched as the value "E" (-> NA + a warning for every file):
        manual_tau_value <-  grab("_mt([0-9.][0-9.eE+-]*)_")
        ##
        algorithm_token <-  grab("_ba(ke|ce|cr|cl|sn)_")
        legacy_objective <-  if (grepl(pattern = "_toct_", x = x)) "ChEES_per_tau" else
                             if (grepl(pattern = "_toch", x = x)) "ChEES" else "KE"
        algorithm_map <-  c(ke = "KE", ce = "ChEES", cr = "CHESSR", cl = "CHESSR_log", sn = "SNAPER")
        ##
        ## ---- nuisance mass ("_nuUD_" etc.; codes from R_fn_map_M_typ in ps_7_MCMC_settings_BayesMVP_functions.R):
        ##
        metric_type_nuisance_token <-  grab("_nu(UD|E|U|H)_")
        metric_type_nuisance_map <-  c(UD = "uniform_diag", E = "Empirical", U = "unit", H = "Hessian")
        metric_type_nuisance <-  if (is.na(metric_type_nuisance_token)) NA_character_ else unname(metric_type_nuisance_map[metric_type_nuisance_token])
        burnin_algorithm <-  if (is.na(algorithm_token)) {
            if (legacy_objective == "ChEES_per_tau") "CHESSR_log" else legacy_objective
        } else unname(algorithm_map[algorithm_token])
        ## Retain tau_objective only as a legacy reader/reporting column, never as a second run selector.
        ##
        tibble::tibble( model_type = grab("^ps7_run_(.+)_N[0-9]+_"),
                    N = as.numeric(grab("_N([0-9]+)_")),
                    tau_adaptation_version = as.numeric(grab("_ta([0-9]+)_", default = "1")),
                    burnin_schedule_version = as.numeric(grab("_bs([0-9]+)a", default = "0")),
                    n_adapt = as.numeric(grab("_bs[0-9]+a([0-9]+)m")),
                    metric_adaptation_end_iter = as.numeric(grab("_bs[0-9]+a[0-9]+m([0-9]+)c")),
                    centre_adaptation_end_iter = as.numeric(grab("_bs[0-9]+a[0-9]+m[0-9]+c([0-9]+)_")),
                    ## chains / threads / chunks ("_cb4_s180_wb8_s1_kb25_s25_" = burn-in 4 chains, 180 sampling chains,
                    ## 8 WCP threads in burn-in, 1 in sampling, 25 chunks in burn-in, 25 in sampling):
                    n_chains_burnin = as.numeric(grab("_cb([0-9]+)_")),
                    n_chains_sampling = as.numeric(grab("_cb[0-9]+_s([0-9]+)_")),
                    n_threads_WCP_burnin = as.numeric(grab("_wb([0-9]+)_")),
                    n_threads_WCP_sampling = as.numeric(grab("_wb[0-9]+_s([0-9]+)_")),
                    num_chunks_burnin = as.numeric(grab("_kb([0-9]+)_")),
                    num_chunks_sampling = as.numeric(grab("_kb[0-9]+_s([0-9]+)_")),
                    adapt_delta = as.numeric(sub("^\\.", "0.", grab("_AD([.0-9]+)_"))),
                    test_perm_override = grab("_tp([0-9-]+)_"),
                    ## absent "_ifkf" = the default kick_flow_kick
                    diffusion_HMC_integrator = if (grepl(pattern = "_ifkf", x = x)) "flow_kick_flow" else "kick_flow_kick",
                    ## absent "_toch" = the original kinetic-energy objective
                    burnin_algorithm = burnin_algorithm,
                    tau_objective = if (is.na(algorithm_token)) legacy_objective else burnin_algorithm,
                    ## absent "_tw1" = median across chains with the accept indicator
                    tau_weight_by_p_jump = grepl(pattern = "_tw1", x = x),
                    ## absent "_mL" = tau adapted
                    manual_L = if (is.na(manual_L)) NA_real_ else as.numeric(manual_L),
                    manual_tau_value = if (is.na(manual_tau_value)) NA_real_ else as.numeric(manual_tau_value),
                    n_burnin = as.numeric(grab("_b([0-9]+)_it")),
                    n_iter = as.numeric(grab("_it([0-9]+)_")),
                    learning_rate = as.numeric(sub("^\\.", "0.", grab("_LR([.0-9]+)_AD"))),
                    ## "_Li.15_d" = the learning-rate hold value and its iteration ("d" = default)
                    ## "_LiLR" = the initial LR is the LR itself (written instead of "_Li<LR>_d" only when a name is too long)
                    learning_rate_initial = if (grepl(pattern = "_LiLR_", x = x)) as.numeric(sub("^\\.", "0.", grab("_LR([.0-9]+)_AD"))) else
                                            as.numeric(sub("^\\.", "0.", grab("_Li([.0-9]+)_"))),
                    tau_initial = if (identical(grab("_ti([^_]+)_"), "hpi")) as.character(pi / 2) else grab("_ti([^_]+)_"),
                    metric_estimator = if (grepl(pattern = "_Mp(_tp[0-9-]+)?(_tr[OS])?(_nER)?(_c[F0])?(_cfi[0-9]+)?(_pa[0-9]+)?(_pb[0-9]+)?(_pL[0-9]+)?(_sT)?(_tbbC)?(_JgA)?(_ll1)?_run", x = x)) "pooled" else
                                       if (grepl(pattern = "_Ms(_tp[0-9-]+)?(_tr[OS])?(_nER)?(_c[F0])?(_cfi[0-9]+)?(_pa[0-9]+)?(_pb[0-9]+)?(_pL[0-9]+)?(_sT)?(_tbbC)?(_JgA)?(_ll1)?_run", x = x)) "chain_mean_scaled" else "chain_mean",
                    tau_ramp = if (grepl(pattern = "_trS(_nER)?(_c[F0])?(_cfi[0-9]+)?(_pa[0-9]+)?(_pb[0-9]+)?(_pL[0-9]+)?(_sT)?(_tbbC)?(_JgA)?(_ll1)?_run", x = x)) "staged" else
                               if (grepl(pattern = "_trO(_nER)?(_c[F0])?(_cfi[0-9]+)?(_pa[0-9]+)?(_pb[0-9]+)?(_pL[0-9]+)?(_sT)?(_tbbC)?(_JgA)?(_ll1)?_run", x = x)) "original" else NA_character_,
                    metric_type_nuisance = metric_type_nuisance,
                    ## absent "_nER" = eps re-initialised at the ChEES handover (the default)
                    eps_reinit_at_ChEES_handover = !grepl(pattern = "_nER(_c[F0])?(_cfi[0-9]+)?(_pa[0-9]+)?(_pb[0-9]+)?(_pL[0-9]+)?(_sT)?(_tbbC)?(_JgA)?(_ll1)?_run", x = x),
                    ## absent "_cF"/"_c0" = nuisance centre never frozen (Earlier runs)
                    theta_hat_us_rule = if (grepl(pattern = "_cF(_cfi[0-9]+)?(_pa[0-9]+)?(_pb[0-9]+)?(_pL[0-9]+)?(_sT)?(_tbbC)?(_JgA)?(_ll1)?_run", x = x)) "running_mean_frozen" else
                                        if (grepl(pattern = "_c0(_cfi[0-9]+)?(_pa[0-9]+)?(_pb[0-9]+)?(_pL[0-9]+)?(_sT)?(_tbbC)?(_JgA)?(_ll1)?_run", x = x)) "zero" else "running_mean",
                    ## New schedules record the effective centre endpoint in _bs; legacy overrides use _cfi.
                    theta_hat_us_freeze_iter = as.numeric(grab("_cfi([0-9]+)_",
                        default = grab("_bs[0-9]+a[0-9]+m[0-9]+c([0-9]+)_"))),
                    ## "_clip<k>_int<m>" = end of the initial epsilon-only ramp (clip_iter) and the gap from there to the start of adaptive tau (int):
                    clip_iter = as.numeric(grab("_clip([0-9]+)_int[0-9]+_")),
                    int = as.numeric(grab("_clip[0-9]+_int([0-9]+)_")),
                    ## absent "_pa"/"_pb" = all post-adaptation iterations kept / 125-iteration pre-burnin
                    ## Legacy _pa filenames remain accepted by the patterns above, but are not an active setting.
                    pre_burnin_n_iter = as.numeric(grab("_pb([0-9]+)_")),
                    pre_burnin_L = as.numeric(grab("_pL([0-9]+)_")),
                    ## "_sT" = one tau_ii per burn-in iteration shared by all burn-in chains
                    share_tau_ii_across_chains_in_burnin = grepl(pattern = "_sT(_tbbC)?(_JgA)?(_ll1)?_run", x = x),
                    ## "_tbbC" = burn-in TBB pool of exactly n_chains_burnin threads
                    burnin_TBB_pool_equals_n_chains = grepl(pattern = "_tbbC(_JgA)?(_ll1)?_run", x = x),
                    ## "_JgA" = correlation-transform Jacobian gradient by autodiff (absent = "num_diff", finite differences)
                    J_grad_option = if (grepl(pattern = "_JgA(_ll1)?_run", x = x)) "autodiff" else "num_diff",
                    autodiff_fallback = autodiff_fallback,
                    ## "_ll1" = sampler retained log-likelihood traces; absence keeps the current PS7 default, FALSE.
                    store_log_lik_trace = grepl(pattern = "_ll1_run", x = x),
                    ## "_rb1" = jittered burn-in tau (with "_sT": one tau per iteration shared by all chains); absent = fixed length
                    randomize_tau_burnin = randomize_tau_burnin,
                    ## "_tsG" = tau_sampling_scale "gaussian_matched" applied at the switch to sampling; absent = "none"
                    tau_sampling_scale = tau_sampling_scale,
                    ## "_ab1" = tau ADAM bias correction by PERFORMED updates (fix); absent = by iteration index
                    tau_adam_bias_correction = tau_adam_bias_correction,
                    ## "_tbJ" = tau_adaptation_block "joint" (criterion computed on main AND nuisance parameters concatenated,
                    ## still adapting the single joint tau;, EXPERIMENTAL); absent = "main" (the existing behaviour)
                    tau_adaptation_block = tau_adaptation_block,
                    run = as.numeric(grab("_run([0-9]+)$")),
                    ## the file name encodes EVERY setting, so the name without "_runK" identifies one configuration:
                    configuration = sub(pattern = "_run[0-9]+$", replacement = "", x = run_file_name),
                    file = run_file_name)

}
##
#' Select rows using named setting levels, without opening any saved fit
#' @param settings A table of parsed filename settings or cached run rows.
#' @param include_settings Named list of allowed values. Different settings are combined with AND; values within a setting with OR.
#' @param exclude_settings Named list of disallowed values. Matching ANY exclusion removes that row. NULL/list() means no filter.
#' @return Logical row-selection vector. Unknown names and empty level vectors are errors; NA can be selected explicitly.
fn_ps7_setting_filter_index <-  function( settings,
                                         include_settings = NULL,
                                         exclude_settings = NULL
) {

        selected <-  rep(x = TRUE, times = nrow(x = settings))
        filters <-  list(include_settings = include_settings, exclude_settings = exclude_settings)
        for (filter_name in names(x = filters)) {
            setting_filter <-  filters[[filter_name]]
            if (is.null(x = setting_filter) || identical(x = setting_filter, y = list())) next
            if (!is.list(x = setting_filter) || is.null(x = names(x = setting_filter)) ||
                anyNA(x = names(x = setting_filter)) || any(names(x = setting_filter) == "") ||
                anyDuplicated(x = names(x = setting_filter))) {
                stop(filter_name, " must be a named list, e.g. list(diffusion_HMC_integrator = 'kick_flow_kick').")
            }
            unknown_settings <-  setdiff(x = names(x = setting_filter), y = names(x = settings))
            if (length(x = unknown_settings)) stop("Unknown ", filter_name, " settings: ", paste(unknown_settings, collapse = ", "))
            ##
            for (setting_name in names(x = setting_filter)) {
                selected_levels <-  setting_filter[[setting_name]]
                if (!is.atomic(x = selected_levels) || !length(x = selected_levels) || !is.null(x = dim(x = selected_levels))) {
                    stop(filter_name, "$", setting_name, " must contain at least one level; omit the setting to leave it unrestricted.")
                }
                matches <-  as.character(x = settings[[setting_name]]) %in% as.character(x = selected_levels)
                selected <-  selected & if (filter_name == "include_settings") matches else !matches
            }
        }
        ##
        return(selected)

}
##
#' fn_ps7_extract_tau_sweep
#' @param dir       directory of saved ps7 runs.
#' @param pattern   file-name filter.
#' @param cache     where to save the extracted frame; NULL = do not save.
#' @param verbose   print progress every 25 files.
#' @param include_burnin_history retain the small tau/epsilon/L and cost-objective histories in a list column.
#' @param reconstruct_ta2_superchains explicitly confirm that old ta2 runs used the audited endpoint initializer with ZERO nuisance jitter.
#' @param include_settings Allowed parsed filename levels, e.g. list(diffusion_HMC_integrator = "kick_flow_kick", n_burnin = 250).
#' @param exclude_settings Disallowed parsed filename levels, e.g. list(diffusion_HMC_integrator = "flow_kick_flow").
#'   Both filters run BEFORE readRDS. NULL/list() leaves that filter unrestricted. Exclusions take precedence.
#' @return tibble, one row per saved run, settings + metrics.
#' @export
fn_ps7_extract_tau_sweep <-  function( dir,
                                      pattern = "N2500",
                                      cache = NULL,
                                      verbose = TRUE,
                                      include_burnin_history = FALSE,
                                      reconstruct_ta2_superchains = FALSE,
                                      include_settings = NULL,
                                      exclude_settings = NULL
) {

        files <-  list.files(path = dir, pattern = pattern)
        if (length(x = files) == 0) stop("no files matching '", pattern, "' in ", dir)
        file_settings <-  dplyr::bind_rows(lapply(X = files, FUN = fn_ps7_parse_run_name))
        selected <-  fn_ps7_setting_filter_index( settings = file_settings,
                                                  include_settings = include_settings,
                                                  exclude_settings = exclude_settings)
        files <-  files[selected]
        file_settings <-  file_settings[selected, , drop = FALSE]
        if (!length(x = files)) stop("No files remain after include_settings / exclude_settings; no saved fits were opened.")
        if (verbose) {
            cat("\nFilename filters kept ", sum(selected), " / ", length(selected), " entries; skipped ", sum(!selected),
                " before opening any saved fits.\n", sep = "")
            cat("Extracting ", length(files), " entries across the N and burn-in lengths selected below.\n", sep = "")
            cat("One run file is one fit/seed, not one configuration. Readability is checked below.\n")
            file_count_summary <-  file_settings %>%
                dplyr::group_by(.data$N, .data$n_burnin) %>%
                dplyr::summarise(run_files = dplyr::n(), configurations = dplyr::n_distinct(.data$configuration), .groups = "drop")
            print(x = file_count_summary, n = Inf, width = Inf)
            cat("The regression may subsequently select only one N and burn-in length.\n")
        }
        ## Source the updated package R helper when using live R functions without reinstalling BayesMVP.
        for (helper_name in c("fn_nested_rhat_grouping_from_burnin", "fn_nested_rhat_from_draws_array")) {
            helper_function <-  get0(x = helper_name, mode = "function", inherits = TRUE)
            if (is.null(x = helper_function) && requireNamespace(package = "BayesMVP", quietly = TRUE)) {
                helper_function <-  get0(x = helper_name, envir = asNamespace(ns = "BayesMVP"), inherits = FALSE)
            }
            if (is.null(x = helper_function)) stop("Source BayesMVP/inst/BayesMVP/R/R_fn_create_superchain_ids.R before extracting runs.")
            assign(x = helper_name, value = helper_function)
        }
        ##
        one <-  function(x) if (is.null(x) || length(x) != 1) NA_real_ else as.numeric(x)
        finite_max <-  function(values) if (length(x = values) > 0 && all(is.finite(x = values))) max(values) else NA_real_
        rows <-  vector(mode = "list", length = length(x = files))
        ##
        for (i in seq_along(along.with = files)) {

            settings <-  file_settings[i, , drop = FALSE]
            r <-  try(readRDS(file = file.path(dir, files[i])), silent = TRUE)
            ##
            ## Keep the same columns for unreadable files, so reporting can exclude them without rbind failing.
            readable <-  !inherits(x = r, what = "try-error")
            if (!readable) r <-  NULL
            e <-  r$efficiency_info
            ## Same probabilities as the live post-burn-in print, averaged over ALL sampling iterations and chains.
            ## In joint HMC p_jump_main is the joint proposal's probability; in partitioned HMC it is the main block's.
            sampling_acceptance_probabilities <-  r$sampler_diagnostics$sampling$p_jump_main
            sampling_acceptance_probability <-  if (is.numeric(x = sampling_acceptance_probabilities) &&
                length(x = sampling_acceptance_probabilities) > 0 && all(is.finite(x = sampling_acceptance_probabilities)) &&
                all(sampling_acceptance_probabilities >= 0 & sampling_acceptance_probabilities <= 1)) {
                mean(x = sampling_acceptance_probabilities)
            } else NA_real_
            ##
            ## Prefer the actual initializer mapping saved by new fits. Never infer groups from posterior draws.
            nested_rhat_grouping <-  r$HMC_info$nested_rhat_grouping
            if (is.null(x = nested_rhat_grouping) && readable && isTRUE(x = reconstruct_ta2_superchains) &&
                isTRUE(x = settings$tau_adaptation_version == 2)) {
                nested_rhat_grouping <-  fn_nested_rhat_grouping_from_burnin(
                    n_chains_sampling = r$settings$n_chains_sampling,
                    n_superchains = r$HMC_info$n_superchains,
                    n_chains_burnin = r$settings$n_chains_burnin,
                    nuisance_jitter_scale = 0,
                    source = "reconstructed_audited_ta2_zero_jitter")
            }
            if (is.null(x = nested_rhat_grouping)) {
                nested_rhat_grouping <-  list(status = "unavailable_initial_state_mapping_not_recorded",
                                             source = "not_recorded")
            }
            interest_trace <-  r$trace_gq
            if (length(dim(interest_trace)) == 3L && length(r$diagnostic_parameter_names) > 0L) {
                interest_indices <-  match(x = r$diagnostic_parameter_names, table = r$tibble_gq$parameter)
                if (anyNA(interest_indices)) stop("Saved diagnostic parameter names do not match the generated-quantity summary: ", files[i])
                interest_trace <-  interest_trace[, , interest_indices, drop = FALSE]
            }
            ## Equal-sized groups use the earliest chain indices in each group. ONLY nested R-hat uses this subset.
            max_nested_rhat_interest <-  if (length(x = dim(x = interest_trace)) == 3) finite_max(
                fn_nested_rhat_from_draws_array(draws_array = interest_trace, nested_rhat_grouping = nested_rhat_grouping)) else NA_real_
            max_nested_rhat_main <-  if (length(x = dim(x = r$trace_main)) == 3) finite_max(
                fn_nested_rhat_from_draws_array(draws_array = r$trace_main, nested_rhat_grouping = nested_rhat_grouping)) else NA_real_
            rows[[i]] <-  dplyr::bind_cols(settings,
                               tibble::tibble( readable = readable,
                                           min_ESS = one(r$min_ESS),
                                           max_Rhat = one(r$max_Rhat),
                                           max_nRhat = max_nested_rhat_interest,
                                           max_nRhat_stored = one(r$max_nRhat),
                                           n_superchains_nominal = one(r$HMC_info$n_superchains),
                                           n_superchains_diagnostic = one(nested_rhat_grouping$n_superchains),
                                           n_chains_per_superchain_diagnostic = one(nested_rhat_grouping$n_chains_per_superchain),
                                           n_chains_nested_rhat = one(nested_rhat_grouping$n_chains_used),
                                           n_chains_omitted_nested_rhat = one(nested_rhat_grouping$n_chains_omitted),
                                           nested_rhat_grouping_source = nested_rhat_grouping$source,
                                           nested_rhat_grouping_status = nested_rhat_grouping$status,
                                           ##
                                           ## the same three diagnostics over the MAIN (raw) parameters:
                                           min_ESS_main = one(e$Min_ESS_main),
                                           max_Rhat_main = one(e$Max_rhat_main),
                                           max_nRhat_main = max_nested_rhat_main,
                                           max_nRhat_main_stored = one(e$Max_nested_rhat_main),
                                           ##
                                           ESS_per_grad_samp = one(r$ESS_per_grad_samp),
                                           ESS_per_sec_samp = one(r$ESS_per_sec_samp),
                                           sampling_acceptance_probability = sampling_acceptance_probability,
                                           pct_divs = one(r$divergences$pct_divs),
                                           time_burnin = one(e$time_burnin),
                                           time_sampling = one(e$time_sampling),
                                           time_total = one(e$time_total),
                                           time_summaries = one(e$time_summaries),
                                           L_main_samp = one(e$L_main_during_sampling),
                                           L_main_burnin = one(e$L_main_during_burnin),
                                           eps_main = one(r$HMC_info$eps_main),
                                           tau_main = one(r$HMC_info$tau_main),
                                           grad_evals_per_sec = one(e$grad_evals_per_sec)))
            ##
            ## Optional small adaptation histories, never posterior draw arrays.
            if (isTRUE(x = include_burnin_history)) {
                history <-  if (!readable) NULL else r$efficiency_info[
                    c("tau_main_during_burnin_vec", "eps_main_during_burnin_vec", "L_main_during_burnin_vec",
                      "ChEES_criterion_ema_vec", "ChEES_per_tau_gradient_vec")]
                rows[[i]] <-  dplyr::mutate(.data = rows[[i]], burnin_history = list(history))
            }
            ## the traces dominate the footprint; drop them before the next file
            rm(r, interest_trace, sampling_acceptance_probabilities)
            # if (i %% 10 == 0) gc(verbose = FALSE)   ## removed: explicit collections are not needed here
            if (verbose && i %% 25 == 0) cat("  Read ", i, "/", length(x = files), " matching files\n", sep = "")

        }
        ##
        out <-  dplyr::bind_rows(rows)
        if (verbose) {
            cat("\nExtraction complete: ", nrow(out), " file entries; ", sum(out$readable), " readable fits; ",
                sum(!out$readable), " unreadable entries.\n", sep = "")
            cat("\nNested R-hat grouping (ESS, ordinary R-hat and timings still use ALL sampling chains):\n")
            print(unique(x = out[c("N", "n_superchains_diagnostic", "n_chains_per_superchain_diagnostic",
                                   "n_chains_nested_rhat", "n_chains_omitted_nested_rhat", "nested_rhat_grouping_status")]), n = Inf, width = Inf)
        }
        ##
        ## one readable label per adaptation scheme, matching summarize_ps7_results()
        out$tau_scheme <-  ifelse( test = !is.na(out$manual_tau_value),
                                  yes = paste0("fixed_tau_", out$manual_tau_value),
                                  no = ifelse(test = !is.na(out$manual_L),
                                  yes = paste0("fixed_L_", out$manual_L),
                                  no = paste0(out$tau_objective, ifelse(out$tau_weight_by_p_jump, "_pjump", ""))))
        ##
        if (!is.null(cache)) saveRDS(object = out, file = cache)
        return(out)

}






















