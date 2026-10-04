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
## ---- Every quantity of the naming list, flattened from the time-criterion record of a CHESSR_time / SNAPER_time run into one value
##      each (NA when absent; a vector, such as sampling_timing_probe_L_values, becomes one comma-separated string). The record is NicoStan's
##      model_results$time_criterion, saved in the run as time_criterion (see fn_ps7_time_criterion_record_from_results in
##      ps_7_MCMC_settings_BayesMVP_functions.R):
##        time_criterion_settings                 the settings as NicoStan validated them;
##        time_criterion_sampling_quantities      the sampling quantities the criterion used, each with its source ("user_supplied",
##                                                "previous_run", "sampling_timing_probe", ...), and n_iter_sampling_for_time_criterion;
##        sampling_timing_probe                   the probe result (NULL when the probe did not run);
##        sampling_timing_probe_wall_time, sampling_timing_probe_fraction_of_burnin, time_burnin_without_sampling_timing_probe;
##        burnin                                  the burn-in record: n_iter_burnin, time_per_leapfrog_step_burnin_at_end (the online
##                                                estimate at the end of burn-in and its status), and the at_handover / at_end
##                                                snapshots (the handover is the first tau update, clip_iter + int).
##      Each quantity is read from its place in that record (ps7_time_criterion_record_value_paths, first path found); a record of
##      another layout is searched by name instead (at the top level, then as "<name>" inside "at_handover" / "at_end" for the
##      "_at_handover" / "_at_end" quantities, then anywhere outside "at_handover"). Used by the check after the fit, the console
##      message and fn_ps7_extract_tau_sweep():
##
ps7_time_criterion_record_value_paths <-  list(
        ## the sampling timing probe (settings as run, results, wall time; "sampling_timing_probe_ran" is derived below):
        sampling_timing_probe_status               = list(c("sampling_timing_probe", "status")),
        sampling_timing_probe_L_values_used        = list(c("sampling_timing_probe", "sampling_timing_probe_L_values")),
        sampling_timing_probe_n_iter_per_L_used    = list(c("sampling_timing_probe", "sampling_timing_probe_n_iter_per_L")),
        sampling_timing_probe_n_iter_total         = list(c("sampling_timing_probe", "sampling_timing_probe_n_iter_total")),
        sampling_timing_probe_wall_time            = list("sampling_timing_probe_wall_time", c("sampling_timing_probe", "sampling_timing_probe_wall_time")),
        sampling_timing_probe_fraction_of_burnin   = list("sampling_timing_probe_fraction_of_burnin"),
        time_burnin_without_sampling_timing_probe  = list("time_burnin_without_sampling_timing_probe"),
        ## the sampling quantities the criterion used, and where each came from:
        time_per_leapfrog_step_sampling            = list(c("time_criterion_sampling_quantities", "time_per_leapfrog_step_sampling")),
        time_per_leapfrog_step_sampling_source     = list(c("time_criterion_sampling_quantities", "time_per_leapfrog_step_sampling_source")),
        time_per_iter_overhead_sampling            = list(c("time_criterion_sampling_quantities", "time_per_iter_overhead_sampling")),
        time_per_iter_overhead_sampling_source     = list(c("time_criterion_sampling_quantities", "time_per_iter_overhead_sampling_source")),
        time_per_iter_summaries_sampling           = list(c("time_criterion_sampling_quantities", "time_per_iter_summaries_sampling")),
        time_per_iter_summaries_sampling_source    = list(c("time_criterion_sampling_quantities", "time_per_iter_summaries_sampling_source")),
        sampling_overhead_in_leapfrog_steps        = list(c("time_criterion_sampling_quantities", "sampling_overhead_in_leapfrog_steps")),
        sampling_overhead_in_leapfrog_steps_source = list(c("time_criterion_sampling_quantities", "sampling_overhead_in_leapfrog_steps_source")),
        ## the burn-in leapfrog time (online estimate), at the handover and at the end of burn-in:
        time_per_leapfrog_step_burnin_at_handover  = list(c("burnin", "at_handover", "time_per_leapfrog_step_burnin"), c("at_handover", "time_per_leapfrog_step_burnin")),
        time_per_leapfrog_step_burnin              = list(c("burnin", "time_per_leapfrog_step_burnin_at_end", "time_per_leapfrog_step_burnin"),
                                                          c("time_per_leapfrog_step_burnin_at_end", "time_per_leapfrog_step_burnin")),
        time_per_leapfrog_step_burnin_status       = list(c("burnin", "time_per_leapfrog_step_burnin_at_end", "status"),
                                                          c("time_per_leapfrog_step_burnin_at_end", "status")),
        ## iteration counts, the ESS target and where n_iter_sampling_for_time_criterion came from:
        n_iter_burnin                              = list(c("burnin", "n_iter_burnin")),
        n_iter_sampling_for_time_criterion         = list(c("time_criterion_sampling_quantities", "n_iter_sampling_for_time_criterion")),
        n_iter_sampling_for_time_criterion_source  = list(c("time_criterion_sampling_quantities", "n_iter_sampling_for_time_criterion_source")),
        ESS_target_for_time_criterion              = list(c("time_criterion_sampling_quantities", "ESS_target_for_time_criterion")),
        ESS_per_iter_sampling_expected             = list(c("time_criterion_sampling_quantities", "ESS_per_iter_sampling_expected")),
        ESS_per_iter_sampling_expected_source      = list(c("time_criterion_sampling_quantities", "ESS_per_iter_sampling_expected_source")),
        ## the criterion's step-size-dependent quantities at the handover and at the end of burn-in:
        tau_offset_from_sampling_overhead_at_handover      = list(c("burnin", "at_handover", "tau_offset_from_sampling_overhead")),
        tau_offset_from_sampling_overhead_at_end           = list(c("burnin", "at_end", "tau_offset_from_sampling_overhead")),
        burnin_to_sampling_leapfrog_time_ratio_at_handover = list(c("burnin", "at_handover", "burnin_to_sampling_leapfrog_time_ratio")),
        burnin_to_sampling_leapfrog_time_ratio_at_end      = list(c("burnin", "at_end", "burnin_to_sampling_leapfrog_time_ratio")),
        time_to_target_ESS_tau_penalty_at_handover         = list(c("burnin", "at_handover", "time_to_target_ESS_tau_penalty_at_tau_main")),
        time_to_target_ESS_tau_penalty_at_end              = list(c("burnin", "at_end", "time_to_target_ESS_tau_penalty_at_tau_main")),
        ## the criterion's estimate of ESS_elasticity_wrt_log_tau (the chain mean at the handover and at the last tau update, and
        ## pooled over the second half of the tau updates next to the mean penalty over the same updates: where the criterion
        ## has converged the two agree):
        ESS_elasticity_wrt_log_tau_at_handover             = list(c("burnin", "at_handover", "ESS_elasticity_wrt_log_tau")),
        ESS_elasticity_wrt_log_tau_at_last_update          = list(c("burnin", "at_last_update", "ESS_elasticity_wrt_log_tau")),
        ESS_elasticity_wrt_log_tau_pooled_over_second_half_of_tau_updates   = list(c("burnin", "at_end", "ESS_elasticity_wrt_log_tau_pooled_over_second_half_of_tau_updates")),
        time_to_target_ESS_tau_penalty_mean_over_second_half_of_tau_updates = list(c("burnin", "at_end", "time_to_target_ESS_tau_penalty_mean_over_second_half_of_tau_updates")),
        ## how the probe's summaries time was obtained, and the burn-in iteration at which the criterion was found to have no
        ## interior optimum (NA = never):
        sampling_timing_probe_summaries_method             = list(c("sampling_timing_probe", "time_per_iter_summaries_sampling_method")),
        iteration_of_no_interior_optimum_warning           = list(c("burnin", "iteration_of_no_interior_optimum_warning")))
##
## every name fn_ps7_time_criterion_values() returns (the derived sampling_timing_probe_ran first, the derived
## time_criterion_fallback_to_rate_criterion last):
ps7_time_criterion_record_value_names <-  c("sampling_timing_probe_ran", names(ps7_time_criterion_record_value_paths),
                                           "time_criterion_fallback_to_rate_criterion")
##
fn_ps7_time_criterion_values <-  function(time_criterion_record) {

        if (!is.list(time_criterion_record)) time_criterion_record <-  list()
        fn_value_at_record_path <-  function(value_path) {
            record_value <-  time_criterion_record
            for (record_path_element in value_path) {
                if (!is.list(record_value) || is.null(names(record_value)) || !record_path_element %in% names(record_value)) return(NULL)
                record_value <-  record_value[[record_path_element]]
            }
            return(record_value)
        }
        fn_find_first_value_by_name_anywhere <-  function(record_sub_list, value_name) {
            if (!is.list(record_sub_list)) return(NULL)
            if (!is.null(names(record_sub_list)) && value_name %in% names(record_sub_list) && !is.null(record_sub_list[[value_name]])) return(record_sub_list[[value_name]])
            for (record_element_index in seq_along(record_sub_list)) {
                if (identical(names(record_sub_list)[record_element_index], "at_handover")) next
                found_value <-  fn_find_first_value_by_name_anywhere(record_sub_list = record_sub_list[[record_element_index]], value_name = value_name)
                if (!is.null(found_value)) return(found_value)
            }
            return(NULL)
        }
        fn_find_value_by_name <-  function(value_name) {
            found_value <-  time_criterion_record[[value_name]]
            if (is.null(found_value) && grepl(pattern = "_at_(handover|end)$", x = value_name)) {
                sub_list_name <-  if (grepl(pattern = "_at_handover$", x = value_name)) "at_handover" else "at_end"
                name_inside_sub_list <-  sub(pattern = "_at_(handover|end)$", replacement = "", x = value_name)
                for (record_sub_list in list(time_criterion_record[[sub_list_name]], time_criterion_record$burnin[[sub_list_name]])) {
                    if (is.list(record_sub_list) && is.null(found_value)) found_value <-  record_sub_list[[name_inside_sub_list]]
                }
            }
            if (is.null(found_value)) found_value <-  fn_find_first_value_by_name_anywhere(record_sub_list = time_criterion_record, value_name = value_name)
            return(found_value)
        }
        fn_single_value_or_comma_separated_text <-  function(record_value) {
            if (is.null(record_value) || length(record_value) == 0) return(NA)
            if (is.list(record_value)) return(NA)
            if (length(record_value) > 1) return(paste(record_value, collapse = ","))
            return(record_value)
        }
        ##
        time_criterion_values_flattened <-  list()
        time_criterion_values_flattened$sampling_timing_probe_ran <-  if ("sampling_timing_probe" %in% names(time_criterion_record)) {
              !is.null(time_criterion_record$sampling_timing_probe)
        } else fn_single_value_or_comma_separated_text(fn_find_value_by_name("sampling_timing_probe_ran"))
        for (value_name in names(ps7_time_criterion_record_value_paths)) {
            found_value <-  NULL
            for (value_path in ps7_time_criterion_record_value_paths[[value_name]]) {
                if (is.null(found_value)) found_value <-  fn_value_at_record_path(value_path = value_path)
            }
            if (is.null(found_value)) found_value <-  fn_find_value_by_name(value_name = value_name)
            time_criterion_values_flattened[[value_name]] <-  fn_single_value_or_comma_separated_text(record_value = found_value)
        }
        ## the rate criterion was used when neither sampling_overhead_in_leapfrog_steps nor time_per_leapfrog_step_sampling was available:
        recorded_fallback <-  fn_find_value_by_name("time_criterion_fallback_to_rate_criterion")
        time_criterion_values_flattened$time_criterion_fallback_to_rate_criterion <-  if (!is.null(recorded_fallback)) as.logical(fn_single_value_or_comma_separated_text(recorded_fallback)) else
            if (length(time_criterion_record) == 0) NA else
            isTRUE(grepl(pattern = "^not_available", x = time_criterion_values_flattened$sampling_overhead_in_leapfrog_steps_source)) &&
            !isTRUE(is.finite(suppressWarnings(as.numeric(time_criterion_values_flattened$time_per_leapfrog_step_sampling))))
        return(time_criterion_values_flattened[ps7_time_criterion_record_value_names])

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
        ## ---- Advanced trajectory-length option tokens (written only when non-default; absent = legacy forward / p1 / b2 / uniform):
        ##      "_tg2" = two-ended gradient, "_tc<p>" = cost exponent p, "_jb3" / "_jb4" = ESJD jump power, "_tjH" = Halton jitter.
        ##      These tokens are stripped before the older end-anchored settings patterns are applied.
        ##
        tau_gradient_estimator <-  if (grepl(pattern = "_tg2(_|$)", x = x)) "two_ended" else "forward"
        tau_cost_token <-  regmatches(x = x, m = regexpr(pattern = "_tc([0-9.]+([eE][-+]?[0-9]+)?)(_|$)", text = x))
        tau_cost_exponent <-  if (length(tau_cost_token) == 0) 1 else as.numeric(sub(pattern = "^_tc", replacement = "", x = sub(pattern = "_$", replacement = "", x = tau_cost_token)))
        esjd_jump_power <-  if (grepl(pattern = "_jb4(_|$)", x = x)) 4 else if (grepl(pattern = "_jb3(_|$)", x = x)) 3 else 2
        tau_jitter_burnin <-  if (grepl(pattern = "_tjH(_|$)", x = x)) "halton" else "uniform"
        x <-  sub(pattern = "_tg2", replacement = "", x = x)
        x <-  sub(pattern = "_tc[0-9.]+([eE][-+]?[0-9]+)?", replacement = "", x = x)
        x <-  sub(pattern = "_jb[34]", replacement = "", x = x)
        x <-  sub(pattern = "_tjH", replacement = "", x = x)
        ##
        ## ---- ADAM token "_A<codes>" (fn_ps7_adam_file_name_code in ps_7_MCMC_settings_BayesMVP_functions.R), written INSTEAD of "_ab1"
        ##      and in its place (after "_rb1" / "_ts..", before "_tbJ"). Codes: "i" = tau ADAM bias correction by iteration index;
        ##      eps update "h" beta1, "g" beta2, "k" denominator constant; tau update "t" beta1, "u" beta2, "w" denominator constant;
        ##      "r" = tau learning-rate restart at the metric freeze; "p" = tau_adaptation_scheme "probe_then_average". Names written before the two updates had separate settings may also
        ##      carry "b" / "v" / "s" (beta1 / beta2 / constant of BOTH updates). The whole token is matched, split into its codes, and
        ##      translated back to "_ab1" (nothing with "i"), so every pattern below sees the names it was written for. Absent = the
        ##      values of every earlier run: beta1 0, beta2 0.95, constant 1e-8 for both updates, no restart:
        ##      "q" = tau_adaptation_scheme "fixed_length_probe_then_decay_and_average".
        ##
        {
            # adam_token_pattern <-  "_A([bghikprstuvw][0-9.eE+-]*)+(_tbJ)?(_run[0-9]+)$"
            adam_token_pattern <-  "_A([bghikpqrstuvw][0-9.eE+-]*)+(_tbJ)?(_run[0-9]+)$"
            adam_settings_from_name <-  list(eps_adam_beta1 = 0, eps_adam_beta2 = 0.95, eps_adam_epsilon = 1e-8,
                                             tau_adam_beta1 = 0, tau_adam_beta2 = 0.95, tau_adam_epsilon = 1e-8,
                                             tau_learning_rate_restart_at_metric_end = FALSE, tau_adaptation_scheme = "adam_decay")
            adam_token <-  regmatches(x = x, m = regexpr(pattern = adam_token_pattern, text = x))
            if (length(x = adam_token) == 1) {
                adam_token_body <-  sub(pattern = "^_A", replacement = "", x = sub(pattern = "(_tbJ)?(_run[0-9]+)$", replacement = "", x = adam_token))
                # adam_codes <-  regmatches(x = adam_token_body, m = gregexpr(pattern = "[bghikprstuvw][0-9.eE+-]*", text = adam_token_body))[[1]]
                adam_codes <-  regmatches(x = adam_token_body, m = gregexpr(pattern = "[bghikpqrstuvw][0-9.eE+-]*", text = adam_token_body))[[1]]
                if (!identical(paste(adam_codes, collapse = ""), adam_token_body)) stop("fn_ps7_parse_run_name: unreadable ADAM token in ", run_file_name)
                adam_iteration_index <-  FALSE
                adam_code_targets <-  list(b = c("eps_adam_beta1", "tau_adam_beta1"), v = c("eps_adam_beta2", "tau_adam_beta2"),
                                           s = c("eps_adam_epsilon", "tau_adam_epsilon"),
                                           h = "eps_adam_beta1", g = "eps_adam_beta2", k = "eps_adam_epsilon",
                                           t = "tau_adam_beta1", u = "tau_adam_beta2", w = "tau_adam_epsilon")
                for (adam_code in adam_codes) {
                    adam_code_letter <-  substring(text = adam_code, first = 1, last = 1)
                    adam_code_value <-  substring(text = adam_code, first = 2)
                    if (adam_code_letter == "i") {
                        adam_iteration_index <-  TRUE
                    } else if (adam_code_letter == "p") {
                        adam_settings_from_name$tau_adaptation_scheme <-  "probe_then_average"
                    } else if (adam_code_letter == "q") {
                        adam_settings_from_name$tau_adaptation_scheme <-  "fixed_length_probe_then_decay_and_average"
                    } else if (adam_code_letter == "r") {
                        adam_settings_from_name$tau_learning_rate_restart_at_metric_end <-  TRUE
                    } else {
                        for (adam_setting_name in adam_code_targets[[adam_code_letter]]) {
                            adam_settings_from_name[[adam_setting_name]] <-  as.numeric(sub(pattern = "^\\.", replacement = "0.", x = adam_code_value))
                        }
                    }
                }
                x <-  sub(pattern = adam_token_pattern, replacement = paste0(if (adam_iteration_index) "" else "_ab1", "\\2\\3"), x = x)
            }
        }
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
        ## algorithm_token <-  grab("_ba(ke|ce|cr|cl|sn)_")
        ## algorithm_token <-  grab("_ba(ke|ce|cr|cl|sn|ct|st)_")
        ## "_baej" = ESJD, "_baex" = ESJD_CHESSR (ps_7_MCMC_settings_BayesMVP_functions.R, fn_ps7_encode_burnin_algorithm):
        ## "_baen" = ESJD_SNAPER (ps_7_MCMC_settings_BayesMVP_functions.R, fn_ps7_encode_burnin_algorithm):
        ## algorithm_token <-  grab("_ba(ke|ce|cr|cl|sn|ct|st|ej|ex|en)_")
        ## "_balq" = LQ_ESSR (ps_7_MCMC_settings_BayesMVP_functions.R, fn_ps7_encode_burnin_algorithm):
        algorithm_token <-  grab("_ba(ke|ce|cr|cl|sn|ct|st|ej|ex|en|lq)_")
        ## ---- CHESSR_time / SNAPER_time with non-default time-criterion settings: "_bac<4 hex>" / "_bas<4 hex>", the final "t"
        ##      replaced by the settings hash (fn_ps7_time_criterion_settings_hash in ps_7_MCMC_settings_BayesMVP_functions.R).
        ##      NA = the default settings ("_bact" / "_bast"), or another algorithm:
        time_criterion_hashed_token <-  grab("_ba([cs][0-9a-f]{4})_")
        time_criterion_settings_hash <-  if (is.na(time_criterion_hashed_token)) NA_character_ else substring(text = time_criterion_hashed_token, first = 2)
        if (is.na(algorithm_token) && !is.na(time_criterion_hashed_token)) algorithm_token <-  paste0(substring(text = time_criterion_hashed_token, first = 1, last = 1), "t")
        legacy_objective <-  if (grepl(pattern = "_toct_", x = x)) "ChEES_per_tau" else
                             if (grepl(pattern = "_toch", x = x)) "ChEES" else "KE"
        ## algorithm_map <-  c(ke = "KE", ce = "ChEES", cr = "CHESSR", cl = "CHESSR_log", sn = "SNAPER")
        ## algorithm_map <-  c(ke = "KE", ce = "ChEES", cr = "CHESSR", cl = "CHESSR_log", sn = "SNAPER",
        ##                     ct = "CHESSR_time", st = "SNAPER_time")
        algorithm_map <-  c(ke = "KE", ce = "ChEES", cr = "CHESSR", cl = "CHESSR_log", sn = "SNAPER",
                            ct = "CHESSR_time", st = "SNAPER_time",
                            ej = "ESJD", ex = "ESJD_CHESSR", en = "ESJD_SNAPER", lq = "LQ_ESSR")
        ##
        ## ---- nuisance mass ("_nuUD_" etc.; codes from R_fn_map_M_typ in ps_7_MCMC_settings_BayesMVP_functions.R):
        ##
        metric_type_nuisance_token <-  grab("_nu(UD|E|U|H)_")
        metric_type_nuisance_map <-  c(UD = "uniform_diag", E = "Empirical", U = "unit", H = "Hessian")
        metric_type_nuisance <-  if (is.na(metric_type_nuisance_token)) NA_character_ else unname(metric_type_nuisance_map[metric_type_nuisance_token])
        ##
        ## ---- main-block mass type and shape ("_mtE_", "_msd_" etc., written as "..._mt<type>_nu<type>_ms<shape>_w..."; codes from
        ##      R_fn_map_M_typ and R_fn_map_M_shp in ps_7_MCMC_settings_BayesMVP_functions.R). The type pattern takes letter codes only, so
        ##      a manual tau value "_mt<number>_" (which starts with a digit or ".") never matches it, and "_msc<number>" (M_decay_scale)
        ##      never matches the shape pattern. NA = a name without the token:
        ##
        metric_type_main_token <-  grab("_mt(UD|E|U|H)_")
        metric_type_main_map <-  c(UD = "uniform_diag", E = "Empirical", U = "unit", H = "Hessian")
        metric_type_main <-  if (is.na(metric_type_main_token)) NA_character_ else unname(metric_type_main_map[metric_type_main_token])
        metric_shape_main_token <-  grab("_ms(d|g)_")
        metric_shape_main_map <-  c(d = "dense", g = "diag")
        metric_shape_main <-  if (is.na(metric_shape_main_token)) NA_character_ else unname(metric_shape_main_map[metric_shape_main_token])
        burnin_algorithm <-  if (is.na(algorithm_token)) {
            if (legacy_objective == "ChEES_per_tau") "CHESSR_log" else legacy_objective
        } else unname(algorithm_map[algorithm_token])
        ##
        ## ---- step-size acceptance mean and divergence-triggered tau shrink, from the "_ta" version digit (see
        ##      fn_ps7_tau_adaptation_version_token in ps_7_MCMC_settings_BayesMVP_functions.R): "_ta5" = harmonic / no shrink,
        ##      "_ta6" = arithmetic / no shrink, "_ta7" = harmonic / shrink on, "_ta8" = geometric / no shrink,
        ##      "_ta9" = geometric / shrink on; "_ta4" and every earlier name = arithmetic /
        ##      shrink on. fn_ps7_extract_tau_sweep() replaces both with the saved run's own record when it has one:
        ##
        eps_tau_shrink_version <-  as.numeric(grab("_ta([0-9]+)_", default = "1"))
        # eps_acceptance_mean <-  if (isTRUE(eps_tau_shrink_version %in% c(5, 7))) "harmonic" else "arithmetic"
        # tau_shrink_on_divergence <-  !isTRUE(eps_tau_shrink_version %in% c(5, 6))
        eps_acceptance_mean <-  if (isTRUE(eps_tau_shrink_version %in% c(5, 7))) "harmonic" else
                                if (isTRUE(eps_tau_shrink_version %in% c(8, 9))) "geometric" else "arithmetic"
        tau_shrink_on_divergence <-  !isTRUE(eps_tau_shrink_version %in% c(5, 6, 8))
        ##
        ## ---- pooled metric estimator: window resets and off-diagonal shrinkage, from the code written straight after "_Mp"
        ##      (see fn_ps7_metric_pooled_file_name_code in ps_7_MCMC_settings_BayesMVP_functions.R): <shrinkage><window code>,
        ##      window code "w" = "none", "s" = "stan_style", "r68-135" = resets after iterations 68 and 135; e.g. "_Mp.5w" = "none" / 0.5.
        ##      "_Mp" alone (every run saved before these options existed) = "stan_style" / 0; a "chain_mean" / "chain_mean_scaled"
        ##      name never has the code (the settings are inert there) and also reads as "stan_style" / 0.
        ##      fn_ps7_extract_tau_sweep() replaces both with a pooled run's own saved record when it has one:
        ##
        metric_pooled_code_pattern <-  "_Mp(A|[0-9.][0-9.eE+-]*)(w|s|r[0-9]+(-[0-9]+)*)_"
        metric_pooled_code <-  regmatches(x = x, m = regexpr(pattern = metric_pooled_code_pattern, text = x))
        if (length(x = metric_pooled_code) == 0) {
            metric_pooled_window_resets <-  "stan_style"
            metric_pooled_offdiagonal_shrinkage <-  0
        } else {
            metric_pooled_window_code <-  sub(pattern = metric_pooled_code_pattern, replacement = "\\2", x = metric_pooled_code)
            metric_pooled_window_resets <-  if (identical(metric_pooled_window_code, "w")) "none" else
                                           if (identical(metric_pooled_window_code, "s")) "stan_style" else
                                           sub(pattern = "^r", replacement = "", x = metric_pooled_window_code)
            metric_pooled_shrinkage_code <-  sub(pattern = metric_pooled_code_pattern, replacement = "\\1", x = metric_pooled_code)
            metric_pooled_offdiagonal_shrinkage <-  if (identical(metric_pooled_shrinkage_code, "A")) "adaptive" else
                as.numeric(sub(pattern = "^\\.", replacement = "0.", x = metric_pooled_shrinkage_code))
        }
        ##
        ## ---- per-iteration metric estimator: the off-diagonal shrinkage written straight after "_Mi" (see
        ##      fn_ps7_metric_per_iteration_file_name_code in ps_7_MCMC_settings_BayesMVP_functions.R), e.g. "_Mi1" = 1; "_Mi" alone = 0.
        ##      The window resets are inert for it and stay "stan_style":
        ##
        metric_per_iteration_code_pattern <-  "_Mi(A|[0-9.][0-9.eE+-]*)_"
        metric_per_iteration_code <-  regmatches(x = x, m = regexpr(pattern = metric_per_iteration_code_pattern, text = x))
        if (length(x = metric_per_iteration_code) == 1) {
            metric_per_iteration_shrinkage_code <-  sub(pattern = metric_per_iteration_code_pattern, replacement = "\\1", x = metric_per_iteration_code)
            metric_pooled_offdiagonal_shrinkage <-  if (identical(metric_per_iteration_shrinkage_code, "A")) "adaptive" else
                as.numeric(sub(pattern = "^\\.", replacement = "0.", x = metric_per_iteration_shrinkage_code))
        }
        metric_pooled_offdiagonal_shrinkage_mode <-  if (identical(metric_pooled_offdiagonal_shrinkage, "adaptive")) "adaptive" else "fixed"
        if (identical(metric_pooled_offdiagonal_shrinkage_mode, "adaptive")) metric_pooled_offdiagonal_shrinkage <-  NA_real_
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
                    ## "_tiA" = tau_initial "adaptive" (tau set at the handover from the lambda_max estimates, see tau_initial_lambda_max_* below):
                    ## tau_initial = if (identical(grab("_ti([^_]+)_"), "hpi")) as.character(pi / 2) else grab("_ti([^_]+)_"),
                    # tau_initial = if (identical(grab("_ti([^_]+)_"), "hpi")) as.character(pi / 2) else
                    #               if (identical(grab("_ti([^_]+)_"), "A")) "adaptive" else grab("_ti([^_]+)_"),
                    ## "_tiA" / "_tiAw.5" = adaptive ("w.5" = lambda_max from the last half of [clip_iter, handover]):
                    tau_initial = if (identical(grab("_ti([^_]+)_"), "hpi")) as.character(pi / 2) else
                                  if (grepl(pattern = "^A(all|w[.0-9]+)?$", x = grab("_ti([^_]+)_"))) "adaptive" else grab("_ti([^_]+)_"),
                    tau_initial_moments_window = if (identical(grab("_ti([^_]+)_"), "Aall")) "all" else
                                                 if (grepl(pattern = "^A(w[.0-9]+)?$", x = grab("_ti([^_]+)_"))) "after_clip_iter" else NA_character_,
                    tau_initial_moments_window_fraction = if (grepl(pattern = "^Aw[.0-9]+$", x = grab("_ti([^_]+)_")))
                                                              as.numeric(sub("^\\.", "0.", sub("^Aw", "", grab("_ti([^_]+)_")))) else
                                                          if (identical(grab("_ti([^_]+)_"), "A")) 1 else NA_real_,
                    ## ADAM settings of the eps and the tau updates and the tau learning-rate restart ("_A" token, read above; absent = 0 / 0.95 /
                    ## 1e-8 for both updates and no restart):
                    eps_adam_beta1 = adam_settings_from_name$eps_adam_beta1,
                    eps_adam_beta2 = adam_settings_from_name$eps_adam_beta2,
                    eps_adam_epsilon = adam_settings_from_name$eps_adam_epsilon,
                    tau_adam_beta1 = adam_settings_from_name$tau_adam_beta1,
                    tau_adam_beta2 = adam_settings_from_name$tau_adam_beta2,
                    tau_adam_epsilon = adam_settings_from_name$tau_adam_epsilon,
                    tau_learning_rate_restart_at_metric_end = adam_settings_from_name$tau_learning_rate_restart_at_metric_end,
                    ## "p" code = tau_adaptation_scheme "probe_then_average"; absent = "adam_decay":
                    tau_adaptation_scheme = adam_settings_from_name$tau_adaptation_scheme,
                    ## "_Mp" followed by the optional pooled-estimator code (e.g. "_Mp.5w"; see metric_pooled_code_pattern above):
                    ## metric_estimator = if (grepl(pattern = "_Mp(_tp[0-9-]+)?(_tr[OS])?(_nER)?(_c[F0])?(_cfi[0-9]+)?(_pa[0-9]+)?(_pb[0-9]+)?(_pL[0-9]+)?(_sT)?(_tbbC)?(_JgA)?(_ll1)?_run", x = x)) "pooled" else
                    metric_estimator = if (grepl(pattern = "_Mp((A|[0-9.][0-9.eE+-]*)(w|s|r[0-9]+(-[0-9]+)*))?(_tp[0-9-]+)?(_tr[OS])?(_nER)?(_c[F0])?(_cfi[0-9]+)?(_pa[0-9]+)?(_pb[0-9]+)?(_pL[0-9]+)?(_sT)?(_tbbC)?(_JgA)?(_ll1)?_run", x = x)) "pooled" else
                                       ## if (grepl(pattern = "_Ms(_tp[0-9-]+)?(_tr[OS])?(_nER)?(_c[F0])?(_cfi[0-9]+)?(_pa[0-9]+)?(_pb[0-9]+)?(_pL[0-9]+)?(_sT)?(_tbbC)?(_JgA)?(_ll1)?_run", x = x)) "chain_mean_scaled" else "chain_mean",
                                       if (grepl(pattern = "_Ms(_tp[0-9-]+)?(_tr[OS])?(_nER)?(_c[F0])?(_cfi[0-9]+)?(_pa[0-9]+)?(_pb[0-9]+)?(_pL[0-9]+)?(_sT)?(_tbbC)?(_JgA)?(_ll1)?_run", x = x)) "chain_mean_scaled" else
                                       ## "_Mi" followed by the optional off-diagonal shrinkage (e.g. "_Mi1"; see metric_per_iteration_code_pattern above):
                                       if (grepl(pattern = "_Mi(A|[0-9.][0-9.eE+-]*)?(_tp[0-9-]+)?(_tr[OS])?(_nER)?(_c[F0])?(_cfi[0-9]+)?(_pa[0-9]+)?(_pb[0-9]+)?(_pL[0-9]+)?(_sT)?(_tbbC)?(_JgA)?(_ll1)?_run", x = x)) "per_iteration" else "chain_mean",
                    tau_ramp = if (grepl(pattern = "_trS(_nER)?(_c[F0])?(_cfi[0-9]+)?(_pa[0-9]+)?(_pb[0-9]+)?(_pL[0-9]+)?(_sT)?(_tbbC)?(_JgA)?(_ll1)?_run", x = x)) "staged" else
                               if (grepl(pattern = "_trO(_nER)?(_c[F0])?(_cfi[0-9]+)?(_pa[0-9]+)?(_pb[0-9]+)?(_pL[0-9]+)?(_sT)?(_tbbC)?(_JgA)?(_ll1)?_run", x = x)) "original" else NA_character_,
                    metric_type_nuisance = metric_type_nuisance,
                    ## "_mt<type>" / "_ms<shape>" = main-block mass type ("Empirical", "Hessian", "unit", "uniform_diag") and shape ("dense", "diag")
                    metric_type_main = metric_type_main,
                    metric_shape_main = metric_shape_main,
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
                    ## Advanced trajectory-length fields; absent tokens resolve to the legacy values.
                    tau_gradient_estimator = tau_gradient_estimator,
                    tau_cost_exponent = tau_cost_exponent,
                    esjd_jump_power = esjd_jump_power,
                    tau_jitter_burnin = tau_jitter_burnin,
                    ## "_tsG" = tau_sampling_scale "gaussian_matched" applied at the switch to sampling; absent = "none"
                    tau_sampling_scale = tau_sampling_scale,
                    ## "_ab1" = tau ADAM bias correction by PERFORMED updates (fix); absent = by iteration index
                    tau_adam_bias_correction = tau_adam_bias_correction,
                    ## "_tbJ" = tau_adaptation_block "joint" (criterion computed on main AND nuisance parameters concatenated,
                    ## still adapting the single joint tau;, EXPERIMENTAL); absent = "main" (the existing behaviour)
                    tau_adaptation_block = tau_adaptation_block,
                    ## "_ta5" / "_ta7" = step size adapted on the HARMONIC mean of the per-chain acceptance probabilities;
                    ## "_ta4", "_ta6" and every earlier name = the arithmetic mean; "_ta8" / "_ta9" = the geometric mean
                    eps_acceptance_mean = eps_acceptance_mean,
                    ## "_ta5" / "_ta6" / "_ta8" = no divergence-triggered tau shrink; "_ta4", "_ta7", "_ta9" and every earlier name = the 0.95 shrink
                    tau_shrink_on_divergence = tau_shrink_on_divergence,
                    ## pooled metric estimator (code after "_Mp"): Welford window resets ("none", "stan_style" or the reset
                    ## iterations joined by "-") and off-diagonal shrinkage; "_Mp" alone and every non-pooled name = "stan_style" / 0
                    metric_pooled_window_resets = metric_pooled_window_resets,
                    metric_pooled_offdiagonal_shrinkage = metric_pooled_offdiagonal_shrinkage,
                    metric_pooled_offdiagonal_shrinkage_mode = metric_pooled_offdiagonal_shrinkage_mode,
                    ## "_bac<4 hex>" / "_bas<4 hex>" = CHESSR_time / SNAPER_time with non-default time-criterion settings, the hash
                    ## of their settings text (the settings themselves are added by fn_ps7_extract_tau_sweep() from the saved run);
                    ## NA = the default settings ("_bact" / "_bast") or another algorithm
                    time_criterion_settings_hash = time_criterion_settings_hash,
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
##
## -| --------- Sampling gradient count: expected leapfrog steps per iteration ----------------------------------------------
##
## The sampler draws tau_ii ~ U(0, 2 * tau) per chain and iteration (randomize_tau_sampling = TRUE, NicoStan's default, which the
## PS7 runner does not change) and takes L_ii = max(1, ceiling(tau_ii / eps)) leapfrog steps
## (EHMC_dual_sampler_fns.hpp, "Compute L"). The saved ESS_per_grad_samp divides by the NOMINAL count tau / eps per iteration,
## which omits the rounding up to whole steps (e.g. tau / eps = 2 -> 2.5 steps executed on average, +25%) and so favours
## short trajectories. With a = 2 * tau / eps and f = floor(a):
##     E[ceiling(a * U)] = ( f * (f + 1) / 2 + (a - f) * (f + 1) ) / a,     U ~ U(0, 1);
## without jitter the count is max(1, ceiling(tau / eps)). Early-terminated (divergent) trajectories are not counted separately.
##
fn_ps7_expected_sampling_leapfrog_steps_per_iteration <-  function( tau,
                                                                    eps,
                                                                    randomize_tau_sampling = TRUE) {

        tau_over_eps <-  tau / eps
        if (!isTRUE(is.finite(tau_over_eps)) || tau_over_eps <= 0) return(NA_real_)
        if (!isTRUE(randomize_tau_sampling)) return(max(1, ceiling(tau_over_eps)))
        upper_bound_in_steps <-  2 * tau_over_eps
        whole_steps_below_upper_bound <-  floor(upper_bound_in_steps)
        return(( whole_steps_below_upper_bound * (whole_steps_below_upper_bound + 1) / 2 +
                 (upper_bound_in_steps - whole_steps_below_upper_bound) * (whole_steps_below_upper_bound + 1) ) / upper_bound_in_steps)

}
##
## -| --------- Sampling gradient count: expected gradient evaluations per iteration (path-aware) --------------------------------
##
## Gradient evaluations per completed, non-divergent trajectory of L leapfrog steps differ by native sampling path
## (NicoStan runtime/MCMC; one C++ call per chain covers all n_iter sampling iterations, and the endpoint-reuse flag lives
## for that call, EHMC_single_threaded_samp_fns.hpp):
##     main-only (no nuisance) or standard joint HMC (diffusion off, not partitioned):  L + 1 every iteration
##     joint diffusion, kick_flow_kick:   L + 1 in the first iteration, then L (the endpoint gradient is reused)
##     joint diffusion, flow_kick_flow:   L + 2 in the first iteration, then L + 1
##     partitioned:                        (L_us + 1) nuisance plus (L_main + 1) main, every iteration (no reuse)
## So per iteration this is E[L] + 0 (kick_flow_kick) or E[L] + 1 (all other paths), with E[L] from
## fn_ps7_expected_sampling_leapfrog_steps_per_iteration() at tau = max(tau, eps) (the C++ raises tau to eps before the
## jitter). For a partitioned run, call it once with tau_main / eps_main and once with tau_us / eps_us.
## When n_iter_per_chain_call is given, the one extra evaluation at the start of each chain's sampling call (joint diffusion
## only) is spread over its iterations (+ 1 / n_iter), giving per chain E[L] x n_iter + 1 (kick_flow_kick) or
## (E[L] + 1) x n_iter + 1 (flow_kick_flow). Rejections restore the stored start evaluation, so they keep the reuse; the rare
## exceptions (which clear it) and early-terminated trajectories are not counted separately.
##
fn_ps7_expected_sampling_gradient_evaluations_per_iteration <-  function( tau,
                                                                          eps,
                                                                          randomize_tau_sampling = TRUE,
                                                                          diffusion_HMC,
                                                                          diffusion_HMC_integrator = "kick_flow_kick",
                                                                          partitioned_HMC,
                                                                          has_nuisance,
                                                                          n_iter_per_chain_call = NULL) {

        if (!isTRUE(is.finite(tau)) || !isTRUE(is.finite(eps)) || eps <= 0) return(NA_real_)
        expected_leapfrog_steps <-  fn_ps7_expected_sampling_leapfrog_steps_per_iteration( tau = max(tau, eps),
                                                                                           eps = eps,
                                                                                           randomize_tau_sampling = randomize_tau_sampling)
        diffusion_HMC_integrator_used <-  if (is.null(diffusion_HMC_integrator)) "kick_flow_kick" else diffusion_HMC_integrator
        is_joint_diffusion_path <-  isTRUE(has_nuisance) && isTRUE(diffusion_HMC) && !isTRUE(partitioned_HMC)
        endpoint_evaluations_per_iteration <-  if (is_joint_diffusion_path && identical(diffusion_HMC_integrator_used, "kick_flow_kick")) 0 else 1
        start_of_call_evaluations_per_iteration <-  if (is_joint_diffusion_path && !is.null(n_iter_per_chain_call) && isTRUE(n_iter_per_chain_call > 0)) {
                                                        1 / n_iter_per_chain_call
                                                    } else 0
        return(expected_leapfrog_steps + endpoint_evaluations_per_iteration + start_of_call_evaluations_per_iteration)

}
##
## The expected sampling gradient evaluations per iteration of one saved PS7 run (main parameters' trajectory; path flags from
## HMC_info, else settings; LC_MVP / MVP / LC_MVOP / MVOP / latent_trait runs always sample nuisance parameters):
##
fn_ps7_expected_sampling_gradient_evaluations_for_run <-  function( run_object,
                                                                    include_start_of_call_evaluation = TRUE) {

        saved_value_or_setting <-  function(field_name, default_value) {
              if (!is.null(run_object$HMC_info[[field_name]])) return(run_object$HMC_info[[field_name]])
              if (!is.null(run_object$settings[[field_name]])) return(run_object$settings[[field_name]])
              return(default_value)
        }
        randomize_tau_sampling_saved <-  if (!is.null(run_object$adaptation$randomize_tau_sampling)) {
                                              run_object$adaptation$randomize_tau_sampling
                                         } else if (!is.null(run_object$settings$randomize_tau_sampling)) {
                                              run_object$settings$randomize_tau_sampling
                                         } else TRUE
        model_type_saved <-  if (!is.null(run_object$settings$model_type)) {
                                  run_object$settings$model_type
                             } else sub(pattern = "^ps7_run_(.*)_N[0-9]+_.*$", replacement = "\\1", x = as.character(run_object$run_id)[1])
        return(fn_ps7_expected_sampling_gradient_evaluations_per_iteration( tau = as.numeric(run_object$HMC_info$tau_main),
                                                                            eps = as.numeric(run_object$HMC_info$eps_main),
                                                                            randomize_tau_sampling = randomize_tau_sampling_saved,
                                                                            diffusion_HMC = saved_value_or_setting("diffusion_HMC", TRUE),
                                                                            diffusion_HMC_integrator = saved_value_or_setting("diffusion_HMC_integrator", "kick_flow_kick"),
                                                                            partitioned_HMC = saved_value_or_setting("partitioned_HMC", FALSE),
                                                                            has_nuisance = model_type_saved %in% c("LC_MVP", "MVP", "LC_MVOP", "MVOP", "latent_trait"),
                                                                            n_iter_per_chain_call = if (isTRUE(include_start_of_call_evaluation)) as.numeric(run_object$n_iter) else NULL))

}
##
## The saved (nominal) ESS_per_grad_samp of a run, re-expressed per EXECUTED gradient: nominal x (tau / eps) / expected steps.
## Uses only the saved tau, eps and ESS_per_grad_samp, so no run is repeated:
##
fn_ps7_ESS_per_grad_samp_expected_steps <-  function(run_object) {

        ## computed DIRECTLY from the saved min ESS, sampling iterations, sampling chains, tau and eps, so it never depends on whether the
        ## saved ESS_per_grad_samp was already corrected (runs saved since the correction carry the marker
        ## ESS_per_grad_samp_gradient_count = "executed_leapfrog_steps_expected"; used only when a field below is missing):
        saved_fields_for_direct_computation <-  list(min_ESS = run_object$min_ESS, n_iter = run_object$n_iter, n_chains = run_object$n_chains,
                                                     tau = run_object$HMC_info$tau_main, eps = run_object$HMC_info$eps_main)
        if (all(vapply(X = saved_fields_for_direct_computation, FUN = function(value) length(value) == 1 && isTRUE(is.finite(as.numeric(value))), FUN.VALUE = TRUE))) {
              ## randomize_tau_sampling_saved <-  if (is.null(run_object$settings$randomize_tau_sampling)) TRUE else run_object$settings$randomize_tau_sampling
              ## expected_steps_per_iteration <-  fn_ps7_expected_sampling_leapfrog_steps_per_iteration( tau = as.numeric(saved_fields_for_direct_computation$tau),
              ##                                                                                         eps = as.numeric(saved_fields_for_direct_computation$eps),
              ##                                                                                         randomize_tau_sampling = randomize_tau_sampling_saved)
              ## return(as.numeric(saved_fields_for_direct_computation$min_ESS) /
              ##        (expected_steps_per_iteration * as.numeric(saved_fields_for_direct_computation$n_iter) * as.numeric(saved_fields_for_direct_computation$n_chains)))
              ## expected GRADIENT EVALUATIONS per iteration for the run's native path (fn_ps7_expected_sampling_gradient_evaluations_per_iteration):
              expected_gradient_evaluations_per_iteration <-  fn_ps7_expected_sampling_gradient_evaluations_for_run(run_object = run_object)
              return(as.numeric(saved_fields_for_direct_computation$min_ESS) /
                     (expected_gradient_evaluations_per_iteration * as.numeric(saved_fields_for_direct_computation$n_iter) * as.numeric(saved_fields_for_direct_computation$n_chains)))
        }
        ## if (identical(run_object$ESS_per_grad_samp_gradient_count, "executed_leapfrog_steps_expected")) return(as.numeric(run_object$ESS_per_grad_samp))
        if (identical(run_object$ESS_per_grad_samp_gradient_count, "gradient_evaluations_expected")) return(as.numeric(run_object$ESS_per_grad_samp))
        ##
        tau_sampling <-  run_object$HMC_info$tau_main
        eps_sampling <-  run_object$HMC_info$eps_main
        ESS_per_grad_samp_nominal <-  run_object$ESS_per_grad_samp
        if (is.null(tau_sampling) || is.null(eps_sampling) || is.null(ESS_per_grad_samp_nominal)) return(NA_real_)
        randomize_tau_sampling <-  if (is.null(run_object$settings$randomize_tau_sampling)) TRUE else run_object$settings$randomize_tau_sampling
        expected_steps <-  fn_ps7_expected_sampling_leapfrog_steps_per_iteration( tau = as.numeric(tau_sampling),
                                                                                  eps = as.numeric(eps_sampling),
                                                                                  randomize_tau_sampling = randomize_tau_sampling)
        ## return(as.numeric(ESS_per_grad_samp_nominal) * (as.numeric(tau_sampling) / as.numeric(eps_sampling)) / expected_steps)
        ## the saved value's own denominator per iteration (expected steps for runs with the earlier marker, else nominal tau / eps),
        ## replaced by the expected gradient evaluations (the start-of-call evaluation is left out when n_iter is not saved):
        saved_denominator_per_iteration <-  if (identical(run_object$ESS_per_grad_samp_gradient_count, "executed_leapfrog_steps_expected")) {
                                                 expected_steps
                                            } else as.numeric(tau_sampling) / as.numeric(eps_sampling)
        expected_gradient_evaluations_per_iteration <-  fn_ps7_expected_sampling_gradient_evaluations_for_run( run_object = run_object,
                                                                                                                include_start_of_call_evaluation = isTRUE(is.finite(as.numeric(run_object$n_iter)[1])))
        return(as.numeric(ESS_per_grad_samp_nominal) * saved_denominator_per_iteration / expected_gradient_evaluations_per_iteration)

}
##
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
            ##
            ## ---- step-size acceptance mean and divergence-triggered tau shrink: the saved run's own record when it has one
            ##      (runs saved since these options existed), otherwise the file name, which gives "arithmetic" / TRUE for
            ##      every earlier run. eps_tau_shrink_source says which of the two was used:
            ##
            eps_tau_shrink_recorded <-  readable && !is.null(x = r$eps_acceptance_mean) && !is.null(x = r$tau_shrink_on_divergence)
            if (eps_tau_shrink_recorded) {
                eps_acceptance_mean_saved <-  as.character(x = r$eps_acceptance_mean)
                tau_shrink_on_divergence_saved <-  as.logical(x = r$tau_shrink_on_divergence)
                if (!identical(eps_acceptance_mean_saved, settings$eps_acceptance_mean) ||
                    !identical(tau_shrink_on_divergence_saved, settings$tau_shrink_on_divergence)) {
                    message(paste0("  ", files[i], ": the saved run records eps_acceptance_mean = ", eps_acceptance_mean_saved,
                                   " | tau_shrink_on_divergence = ", tau_shrink_on_divergence_saved, ", but its file name gives ",
                                   settings$eps_acceptance_mean, " | ", settings$tau_shrink_on_divergence, "; the saved record is used."))
                }
                settings$eps_acceptance_mean <-  eps_acceptance_mean_saved
                settings$tau_shrink_on_divergence <-  tau_shrink_on_divergence_saved
            }
            ##
            ## ---- pooled metric estimator (window resets and off-diagonal shrinkage): a POOLED run's own saved record when it has
            ##      one (runs saved since these options existed), otherwise the file name, which gives "stan_style" / 0 for every
            ##      earlier "_Mp" run. Non-pooled runs keep the file name's values (the settings are inert there, although the
            ##      sampler records them). metric_pooled_settings_source says which of the two was used:
            ##
            metric_pooled_settings_recorded <-  readable && identical(settings$metric_estimator, "pooled") &&
                                               !is.null(x = r$metric_pooled_window_resets) && !is.null(x = r$metric_pooled_offdiagonal_shrinkage)
            if (metric_pooled_settings_recorded) {
                metric_pooled_window_resets_saved <-  r$metric_pooled_window_resets
                metric_pooled_window_resets_saved <-  if (is.character(x = metric_pooled_window_resets_saved)) {
                    paste(metric_pooled_window_resets_saved, collapse = "-")
                } else if (length(x = metric_pooled_window_resets_saved) == 0) {
                    "none"
                } else {
                    paste(sub(pattern = "^ +", replacement = "",
                              x = formatC(sort(unique(as.numeric(metric_pooled_window_resets_saved))), format = "f", digits = 0)),
                          collapse = "-")
                }
                metric_pooled_offdiagonal_shrinkage_mode_saved <-  if (identical(as.character(r$metric_pooled_offdiagonal_shrinkage), "adaptive")) "adaptive" else "fixed"
                metric_pooled_offdiagonal_shrinkage_saved <-  if (identical(metric_pooled_offdiagonal_shrinkage_mode_saved, "adaptive")) NA_real_ else
                                                                as.numeric(x = r$metric_pooled_offdiagonal_shrinkage)
                if (!identical(metric_pooled_window_resets_saved, settings$metric_pooled_window_resets) ||
                    !isTRUE(all.equal(metric_pooled_offdiagonal_shrinkage_saved, settings$metric_pooled_offdiagonal_shrinkage))) {
                    message(paste0("  ", files[i], ": the saved run records metric_pooled_window_resets = ", metric_pooled_window_resets_saved,
                                   " | metric_pooled_offdiagonal_shrinkage = ", metric_pooled_offdiagonal_shrinkage_saved,
                                   ", but its file name gives ", settings$metric_pooled_window_resets, " | ",
                                   settings$metric_pooled_offdiagonal_shrinkage, "; the saved record is used."))
                }
                settings$metric_pooled_window_resets <-  metric_pooled_window_resets_saved
                settings$metric_pooled_offdiagonal_shrinkage <-  metric_pooled_offdiagonal_shrinkage_saved
                settings$metric_pooled_offdiagonal_shrinkage_mode <-  metric_pooled_offdiagonal_shrinkage_mode_saved
            }
            ##
            ## ---- CHESSR_time / SNAPER_time ("_bact" / "_bast", or "_bac<4 hex>" / "_bas<4 hex>"): the time-criterion SETTINGS the
            ##      run used, from its saved time_criterion_inputs (the defaults for a "_bact" / "_bast" run without that record),
            ##      and every quantity of the naming list from the sampler's saved time_criterion record (fn_ps7_time_criterion_values).
            ##      All NA for every other algorithm. A settings text whose hash is not the file name's, or a record that disagrees
            ##      with the saved settings, is reported:
            ##
            time_criterion_is_used <-  isTRUE(settings$burnin_algorithm %in% c("CHESSR_time", "SNAPER_time"))
            time_criterion_inputs_saved <-  if (readable && time_criterion_is_used) r$time_criterion_inputs else NULL
            time_criterion_values <-  fn_ps7_time_criterion_values(time_criterion_record = if (readable && time_criterion_is_used) r$time_criterion else NULL)
            time_criterion_settings_saved <-  time_criterion_inputs_saved$effective_time_criterion_settings
            if (time_criterion_is_used && is.null(time_criterion_settings_saved) && is.na(settings$time_criterion_settings_hash)) {
                time_criterion_settings_saved <-  list(run_sampling_timing_probe = TRUE, sampling_timing_probe_L_values = c(2, 8),
                                                      sampling_timing_probe_n_iter_per_L = 10, n_iter_sampling_for_time_criterion_source = "planned_n_iter")
            }
            fn_saved_time_criterion_setting_value <-  function(setting_name, missing_value) {
                time_criterion_value <-  time_criterion_settings_saved[[setting_name]]
                if (is.null(time_criterion_value) || length(time_criterion_value) == 0) return(missing_value)
                if (length(time_criterion_value) > 1) return(paste(time_criterion_value, collapse = ","))
                return(time_criterion_value)
            }
            time_criterion_settings_columns <-  tibble::tibble(
                run_sampling_timing_probe = as.logical(fn_saved_time_criterion_setting_value("run_sampling_timing_probe", NA)),
                sampling_timing_probe_L_values = as.character(fn_saved_time_criterion_setting_value("sampling_timing_probe_L_values", NA_character_)),
                sampling_timing_probe_n_iter_per_L = as.numeric(fn_saved_time_criterion_setting_value("sampling_timing_probe_n_iter_per_L", NA_real_)),
                n_iter_sampling_for_time_criterion_source = as.character(fn_saved_time_criterion_setting_value("n_iter_sampling_for_time_criterion_source", NA_character_)),
                ESS_target_for_time_criterion = as.numeric(fn_saved_time_criterion_setting_value("ESS_target_for_time_criterion", NA_real_)),
                ## one value per run file: several previous runs are listed by their sorted basenames, as in the settings text:
                time_criterion_previous_run_path_basenames = if (is.null(time_criterion_settings_saved$time_criterion_previous_run_path)) NA_character_ else
                                                                 paste(sort(basename(time_criterion_settings_saved$time_criterion_previous_run_path)), collapse = ","),
                time_criterion_settings_text = if (is.null(time_criterion_inputs_saved$time_criterion_settings_text)) NA_character_ else
                                                   as.character(time_criterion_inputs_saved$time_criterion_settings_text),
                ## the parameters whose min ESS the time criterion refers to: the saved PS7 setting, else the set NicoStan recorded, else
                ## its default "diagnostic" (NA for every other algorithm):
                time_criterion_ess_parameter_set = if (!time_criterion_is_used) NA_character_ else
                                                   if (!is.null(time_criterion_settings_saved$time_criterion_ess_parameter_set))
                                                       paste(time_criterion_settings_saved$time_criterion_ess_parameter_set, collapse = ",") else
                                                   if (readable && !is.null(r$time_criterion$time_criterion_settings$time_criterion_ess_parameter_set))
                                                       paste(r$time_criterion$time_criterion_settings$time_criterion_ess_parameter_set, collapse = ",") else "diagnostic")
            if (!is.null(time_criterion_inputs_saved)) {
                saved_time_criterion_settings_hash <-  time_criterion_inputs_saved$time_criterion_settings_hash
                saved_time_criterion_settings_hash <-  if (is.null(saved_time_criterion_settings_hash) || !nzchar(saved_time_criterion_settings_hash)) NA_character_ else
                                                           saved_time_criterion_settings_hash
                if (!identical(saved_time_criterion_settings_hash, settings$time_criterion_settings_hash)) {
                    message(paste0("  ", files[i], ": the saved time-criterion settings hash ", saved_time_criterion_settings_hash, " is not the file name's ",
                                   settings$time_criterion_settings_hash, " (NA = default settings)."))
                }
            }
            ## (NicoStan's n_iter_sampling_for_time_criterion_source is more detailed, e.g. "planned_n_iter_of_this_run" or
            ##  "ESS_target_for_time_criterion / ESS_per_iter_sampling_expected (previous_run)", so it must START with the setting:)
            if (time_criterion_is_used) {
                for (setting_name in c("n_iter_sampling_for_time_criterion_source", "ESS_target_for_time_criterion")) {
                    recorded_value <-  time_criterion_values[[setting_name]]
                    setting_value <-  time_criterion_settings_columns[[setting_name]]
                    values_are_numbers <-  !is.na(suppressWarnings(as.numeric(recorded_value))) && !is.na(suppressWarnings(as.numeric(setting_value)))
                    values_agree <-  startsWith(x = as.character(recorded_value), prefix = as.character(setting_value)) ||
                                    (values_are_numbers && isTRUE(all.equal(as.numeric(recorded_value), as.numeric(setting_value))))
                    if (!is.na(recorded_value) && !is.na(setting_value) && !values_agree) {
                        message(paste0("  ", files[i], ": the saved time-criterion record has ", setting_name, " = ", recorded_value,
                                       ", but the saved settings give ", setting_value, "; the saved settings are kept."))
                    }
                }
            }
            time_criterion_character_value_names <-  c("sampling_timing_probe_status", "sampling_timing_probe_L_values_used", "sampling_timing_probe_summaries_method",
                                                       grep(pattern = "_source$|_status$", x = names(time_criterion_values), value = TRUE))
            time_criterion_logical_value_names <-  c("sampling_timing_probe_ran", "time_criterion_fallback_to_rate_criterion")
            time_criterion_outcome_names <-  setdiff(x = names(time_criterion_values),
                                                    y = c("n_iter_sampling_for_time_criterion_source", "ESS_target_for_time_criterion"))
            time_criterion_outcomes <-  dplyr::bind_cols(time_criterion_settings_columns,
                                                         tibble::as_tibble(lapply(X = stats::setNames(time_criterion_outcome_names, time_criterion_outcome_names),
                                                                                  FUN = function(value_name) {
                time_criterion_value <-  time_criterion_values[[value_name]]
                if (value_name %in% time_criterion_character_value_names) return(as.character(time_criterion_value))
                if (value_name %in% time_criterion_logical_value_names) return(as.logical(time_criterion_value))
                return(suppressWarnings(as.numeric(time_criterion_value)))
            })))
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
                                           eps_tau_shrink_source = if (eps_tau_shrink_recorded) "saved_run" else "file_name",
                                           metric_pooled_settings_source = if (metric_pooled_settings_recorded) "saved_run" else "file_name",
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
                                           ## ESS_per_grad_samp = one(r$ESS_per_grad_samp),
                                           ## per EXECUTED gradient (expected whole leapfrog steps under the sampling jitter); the saved nominal value is kept:
                                           ESS_per_grad_samp = one(fn_ps7_ESS_per_grad_samp_expected_steps(run_object = r)),
                                           ## the nominal (tau / eps) value: as saved for older runs; recomputed for runs saved with executed steps
                                           ## ESS_per_grad_samp_nominal = one(if (identical(r$ESS_per_grad_samp_gradient_count, "executed_leapfrog_steps_expected"))
                                           ##                                     r$ESS_per_grad_samp * fn_ps7_expected_sampling_leapfrog_steps_per_iteration( tau = r$HMC_info$tau_main,
                                           ##                                                                                                                   eps = r$HMC_info$eps_main,
                                           ##                                                                                                                   randomize_tau_sampling = if (is.null(r$settings$randomize_tau_sampling)) TRUE else r$settings$randomize_tau_sampling) /
                                           ##                                     (r$HMC_info$tau_main / r$HMC_info$eps_main) else r$ESS_per_grad_samp),
                                           ## ... and for runs saved with expected gradient evaluations (marker "gradient_evaluations_expected"):
                                           ESS_per_grad_samp_nominal = one(if (identical(r$ESS_per_grad_samp_gradient_count, "gradient_evaluations_expected"))
                                                                               r$ESS_per_grad_samp * fn_ps7_expected_sampling_gradient_evaluations_for_run(run_object = r) /
                                                                               (r$HMC_info$tau_main / r$HMC_info$eps_main) else
                                                                           if (identical(r$ESS_per_grad_samp_gradient_count, "executed_leapfrog_steps_expected"))
                                                                               r$ESS_per_grad_samp * fn_ps7_expected_sampling_leapfrog_steps_per_iteration( tau = r$HMC_info$tau_main,
                                                                                                                                                             eps = r$HMC_info$eps_main,
                                                                                                                                                             randomize_tau_sampling = if (is.null(r$settings$randomize_tau_sampling)) TRUE else r$settings$randomize_tau_sampling) /
                                                                               (r$HMC_info$tau_main / r$HMC_info$eps_main) else r$ESS_per_grad_samp),
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
                                           ## tau_initial "adaptive" (NA for a numeric tau_initial): lambda_max of each block and the tau set at the handover
                                           ## (before the max_tau ceilings), as recorded by the burn-in (saved run element tau_initial_adaptive):
                                           tau_initial_moments_start_iter = one(r$tau_initial_adaptive$moments_start_iter),
                                           tau_initial_moments_n_updates = one(r$tau_initial_adaptive$n_updates),
                                           tau_initial_lambda_max_main = one(r$tau_initial_adaptive$lambda_max_main),
                                           tau_initial_lambda_max_us = one(r$tau_initial_adaptive$lambda_max_us),
                                           tau_initial_handover_tau_main = one(r$tau_initial_adaptive$tau_main_handover),
                                           tau_initial_handover_tau_us = one(r$tau_initial_adaptive$tau_us_handover),
                                           ## elapsed seconds of the estimator inside time_burnin (per-iteration moments, loading the compiled
                                           ## accumulator, the handover eigenvalues), so that they can be reported or subtracted:
                                           tau_initial_seconds_moments = one(r$tau_initial_adaptive$seconds_moments),
                                           tau_initial_seconds_loader = one(r$tau_initial_adaptive$seconds_loader),
                                           tau_initial_seconds_handover = one(r$tau_initial_adaptive$seconds_handover),
                                           grad_evals_per_sec = one(e$grad_evals_per_sec)))
            ## the time-criterion quantities, after every earlier column (NA for every algorithm other than CHESSR_time / SNAPER_time):
            rows[[i]] <-  dplyr::bind_cols(rows[[i]], time_criterion_outcomes)
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







































