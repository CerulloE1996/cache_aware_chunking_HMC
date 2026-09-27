#### =====================================================================================================================================
## fn_ps7_fit_factorial_regression.R
##
## ---- Analyse individual saved fits; sourcing this file never runs MCMC or changes a benchmark -----------------------------------------
##
## Source fn_ps7_extract_tau_sweep.R and fn_ps7_summarise_runs.R first. Those supply the existing setting-name and ESS-target conventions.
## Dependencies: dplyr, tidyr, tibble, and base R. No Stan compilation or additional regression package is required.
##
#' Factorial regression of a positive PS7 efficiency outcome
#'
#' Fits log(outcome) using categorical settings and a fixed block for the repeated run/seed index.
#' Fits an additive model and, where identifiable, a second model with the requested two-factor interactions.
#' Exponentiated additive coefficients are adjusted ratios under an additive log-response assumption.
#' Interaction coefficients exponentiate to RATIOS OF RATIOS, not unconditional factor effects.
#' HC3 intervals allow unequal residual variances; they do not account for arbitrary within-seed correlation,
#' machine-load drift, uncertainty across datasets, multiple testing, or choosing the best configuration afterwards.
#' With only three seeds, treat intervals as exploratory. Exact matching and per-seed summaries accompany the models.
#'
#' @param runs_table Individual fits from fn_ps7_extract_tau_sweep(), or ps7_configuration_summary$runs. NOT configuration averages.
#' @param outcome Positive numeric column, or time_to_target_ESS, ESS_per_1000_grad, ESS_per_sec, burnin_sec, sampling_sec.
#' @param target_min_ESS Same scalar / named-by-N vector / function(N, Model_type) convention as fn_ps7_summarise_runs().
#'   NULL reuses target_min_ESS already in runs_table; otherwise supply the target explicitly.
#' @param factor_settings Settings to compare. Numeric settings are categorical, so learning-rate effects need not be linear.
#'   NULL uses names(expected_factor_levels), or the original seven settings when no planned grid is supplied.
#' @param interaction_pairs Named-column pairs to add to the separate interaction model. Default pairs follow factor_settings;
#'   "recommended" uses the three learning-rate/objective pairs plus the four original schedule/metric pairs.
#'   list() requests no interactions. Explicitly supplied pairs must name selected factors. Constant factors are omitted explicitly.
#' @param reference_levels Preferred reference levels; absent references use the first observed level, reported in factor_levels.
#' @param expected_factor_levels Optional complete planned grid, including levels not observed yet. Does not change any sampling settings.
#' @param expected_run_indices Expected repeated seed/run labels, used to count missing fits. Default matches the present three-run study.
#' @param matched_setting One or more binary settings for exact same-settings, same-run comparisons, e.g. tau_ramp.
#'   Known settings not selected as factors or with only one observed level are skipped explicitly; unknown names fail.
#'   When theta_hat_us_rule is modelled, its two freeze columns may be dependent descriptors only if each rule has one
#'   fully recorded, consistent endpoint: zero uses 0; running_mean_frozen uses one positive endpoint. Other schedule guards remain.
#' @param confidence_level Confidence level for approximate HC3 coefficient intervals.
#' @param print_table Print coverage, adjusted additive ratios, and matched comparisons. Full results are returned as tibbles and lm objects.
#' @return List containing models, effects, coverage, missing cells, factor levels, residual diagnostics, matched runs and excluded runs.
fn_ps7_fit_factorial_regression <-  function( runs_table,
                                               outcome = "time_to_target_ESS",
                                               target_min_ESS = NULL,
                                               factor_settings = NULL,
                                               interaction_pairs = list(c("tau_ramp", "tau_initial"),
                                                                        c("tau_ramp", "eps_reinit_at_ChEES_handover"),
                                                                        c("tau_objective", "tau_initial"),
                                                                        c("metric_estimator", "learning_rate")),
                                               reference_levels = list(tau_ramp = "original", tau_initial = "pi",
                                                                       tau_objective = "ChEES_per_tau", metric_estimator = "chain_mean_scaled",
                                                                       eps_reinit_at_ChEES_handover = FALSE),
                                               expected_factor_levels = NULL,
                                               expected_run_indices = c(1, 2, 3),
                                               matched_setting = "tau_ramp",
                                               confidence_level = 0.95,
                                               print_table = TRUE
) {

        requested_factors <-  if (!is.null(factor_settings)) factor_settings else names(expected_factor_levels)
        use_canonical_algorithm <-  "burnin_algorithm" %in% names(runs_table) &&
                                    (!"tau_adaptation_version" %in% names(runs_table) || any(runs_table$tau_adaptation_version >= 3, na.rm = TRUE))
        if ("tau_objective" %in% requested_factors && !"burnin_algorithm" %in% requested_factors) use_canonical_algorithm <-  FALSE
        if ("burnin_algorithm" %in% requested_factors) use_canonical_algorithm <-  TRUE
        algorithm_setting <-  if (use_canonical_algorithm) "burnin_algorithm" else "tau_objective"
        if (algorithm_setting == "burnin_algorithm" && missing(reference_levels)) {
                reference_levels$tau_objective <-  NULL
                reference_levels$burnin_algorithm <-  "CHESSR_log"
        }
        if (is.null(x = factor_settings)) {
            factor_settings <-  if (!is.null(x = expected_factor_levels)) names(x = expected_factor_levels) else {
                c("learning_rate", algorithm_setting, "learning_rate_initial", "tau_initial",
                  "metric_estimator", "tau_ramp", "eps_reinit_at_ChEES_handover")
            }
        }
        automatic_interaction_pairs <-  missing(interaction_pairs) || identical(x = interaction_pairs, y = "recommended")
        if (identical(x = interaction_pairs, y = "recommended")) {
            interaction_pairs <-  list(c("learning_rate", "tau_objective"),
                                        c("learning_rate", "learning_rate_initial"),
                                        c("tau_objective", "learning_rate_initial"),
                                        c("tau_ramp", "tau_initial"),
                                        c("tau_ramp", "eps_reinit_at_ChEES_handover"),
                                        c("tau_objective", "tau_initial"),
                                        c("metric_estimator", "learning_rate"))
        }
        if (automatic_interaction_pairs && algorithm_setting == "burnin_algorithm") {
                interaction_pairs <-  lapply(X = interaction_pairs, FUN = function(pair) {
                        replace(x = pair, list = pair == "tau_objective", values = "burnin_algorithm")
                })
        }
        if (!is.list(x = interaction_pairs)) stop("interaction_pairs must be 'recommended', list(), or a list of two-setting pairs.")
        ##
        required_columns <-  c("file", "configuration", "run", "readable", "model_type", "N", "n_burnin", factor_settings)
        missing_columns <-  setdiff(x = required_columns, y = names(x = runs_table))
        if (length(x = missing_columns)) stop("Use the individual runs table. Missing: ", paste(missing_columns, collapse = ", "))
        if (!is.data.frame(x = runs_table) || !nrow(x = runs_table)) stop("runs_table must contain individual fits.")
        if (!is.character(x = outcome) || length(x = outcome) != 1 || is.na(x = outcome)) stop("outcome must name one metric.")
        if (!length(x = factor_settings) || anyNA(x = factor_settings) || anyDuplicated(x = factor_settings)) {
            stop("factor_settings must be nonempty, nonmissing and unique.")
        }
        if (!is.list(x = reference_levels) || anyDuplicated(x = names(x = reference_levels)) ||
            any(lengths(x = reference_levels) != 1) || anyNA(x = unlist(x = reference_levels))) stop("Supply scalar, named reference_levels.")
        if (!is.numeric(x = confidence_level) || length(x = confidence_level) != 1 || !is.finite(x = confidence_level) ||
            confidence_level <= 0 || confidence_level >= 1) stop("confidence_level must be between zero and one.")
        if (!length(x = expected_run_indices) || anyNA(x = expected_run_indices) || anyDuplicated(x = expected_run_indices)) {
            stop("expected_run_indices must contain distinct, nonmissing run indices.")
        }
        if (!exists(x = "fn_ps7_parse_run_name", mode = "function")) stop("Source functions/fn_ps7_extract_tau_sweep.R first.")
        runs_table <-  tibble::as_tibble(x = runs_table)
        ##
        ## ---- Keep N, model, burn-in length and every unmodelled parsed setting fixed -----------------------------------------------------
        ##
        setting_names <-  setdiff(x = names(x = fn_ps7_parse_run_name(x = runs_table$file[1])), y = c("file", "configuration", "run"))
        omitted_default_interaction_pairs <-  character()
        if (automatic_interaction_pairs) {
            selected_pairs <-  vapply(interaction_pairs, function(pair) all(pair %in% factor_settings), FUN.VALUE = FALSE)
            omitted_default_interaction_pairs <-  vapply(interaction_pairs[!selected_pairs], paste, collapse = ":", FUN.VALUE = "")
            interaction_pairs <-  interaction_pairs[selected_pairs]
        }
        valid_pairs <-  vapply(X = interaction_pairs, FUN = function(pair) length(x = pair) == 2 &&
                                  !anyNA(x = pair) && !anyDuplicated(x = pair) && all(pair %in% factor_settings), FUN.VALUE = FALSE)
        if (!all(valid_pairs)) stop("Each interaction_pairs element must name two distinct factor_settings.")
        missing_settings <-  setdiff(x = setting_names, y = names(x = runs_table))
        if (length(x = missing_settings)) stop("Re-extract the saved runs; missing setting columns: ", paste(missing_settings, collapse = ", "))
        if (length(x = setdiff(x = factor_settings, y = setting_names))) stop("factor_settings must name parsed input settings, not outcomes.")
        fixed_setting_names <-  setdiff(x = setting_names, y = factor_settings)
        dependent_settings <-  tibble::tibble()
        centre_endpoint_columns <-  c("centre_adaptation_end_iter", "theta_hat_us_freeze_iter")
        if ("theta_hat_us_rule" %in% factor_settings && all(centre_endpoint_columns %in% fixed_setting_names)) {
            centre_schedule_mapping <-  runs_table %>%
                dplyr::distinct(dplyr::across(.cols = dplyr::all_of(x = c("theta_hat_us_rule", centre_endpoint_columns))))
            centre_rules <-  centre_schedule_mapping$theta_hat_us_rule
            centre_endpoints <-  centre_schedule_mapping$centre_adaptation_end_iter
            freeze_endpoints <-  centre_schedule_mapping$theta_hat_us_freeze_iter
            if (!anyNA(x = centre_rules) && !anyDuplicated(x = centre_rules) &&
                all(centre_rules %in% c("zero", "running_mean_frozen")) &&
                all(is.finite(centre_endpoints)) && all(is.finite(freeze_endpoints)) &&
                all(centre_endpoints == freeze_endpoints) && all(centre_endpoints == floor(x = centre_endpoints)) &&
                all(centre_endpoints[centre_rules == "zero"] == 0) &&
                all(centre_endpoints[centre_rules == "running_mean_frozen"] > 0)) {
                ## These columns describe the modelled centre treatment, not separately identifiable predictor effects.
                ## More than one freeze endpoint WITHIN a rule is not exempted: it still needs filtering or explicit modelling.
                dependent_settings <-  centre_schedule_mapping
                fixed_setting_names <-  setdiff(x = fixed_setting_names, y = centre_endpoint_columns)
            }
        }
        ## The parser retains the legacy objective label as a descriptor of the same algorithm.
        if (all(c("burnin_algorithm", "tau_objective") %in% names(runs_table))) {
                if ("burnin_algorithm" %in% factor_settings) fixed_setting_names <-  setdiff(fixed_setting_names, "tau_objective")
                if ("tau_objective" %in% factor_settings) fixed_setting_names <-  setdiff(fixed_setting_names, "burnin_algorithm")
        }
        ## learning_rate_initial = NULL in the runner means the initial LR IS the LR, so it varies only with learning_rate:
        if (all(c("learning_rate", "learning_rate_initial") %in% names(runs_table)) &&
            isTRUE(all(abs(runs_table$learning_rate_initial - runs_table$learning_rate) < 1e-12))) {
                fixed_setting_names <-  setdiff(fixed_setting_names, "learning_rate_initial")
        }
        varying_fixed_settings <-  fixed_setting_names[vapply(X = runs_table[fixed_setting_names],
            FUN = function(values) dplyr::n_distinct(values) > 1, FUN.VALUE = FALSE)]
        if (length(x = varying_fixed_settings)) {
            stop("Filter to one stratum or explicitly model these varying settings first: ", paste(varying_fixed_settings, collapse = ", "),
                 ". Do not silently pool schedules, chain layouts, models or sample sizes.")
        }
        if (anyNA(x = runs_table[c("run", factor_settings)])) stop("Run indices and factor settings must be recorded for every row.")
        if (any(!runs_table$run %in% expected_run_indices)) stop("Observed run indices are absent from expected_run_indices.")
        duplicate_cells <-  runs_table %>%
            dplyr::count(dplyr::across(.cols = dplyr::all_of(x = c(factor_settings, "run"))), name = "n_records") %>%
            dplyr::filter(.data$n_records > 1)
        if (nrow(x = duplicate_cells)) stop("Duplicate setting/run cells. Filter to one build/dataset/source cohort before analysing.")
        ##
        ## ---- Recalculate the requested target consistently: burn-in + scaled sampling AND summaries ------------------------------------
        ##
        if (!is.null(x = target_min_ESS)) {
            if (!exists(x = "R_fn_ps7_resolve_target_min_ESS", mode = "function")) stop("Source functions/fn_ps7_summarise_runs.R first.")
            runs_table$target_min_ESS <-  R_fn_ps7_resolve_target_min_ESS(N_of_each_run = runs_table$N,
                                                                         model_type_of_each_run = runs_table$model_type,
                                                                         target_min_ESS = target_min_ESS)
        }
        if (identical(x = outcome, y = "time_to_target_ESS")) {
            timing_columns <-  c("target_min_ESS", "time_burnin", "time_sampling", "time_summaries", "min_ESS")
            if (!all(timing_columns %in% names(x = runs_table))) stop("Supply target_min_ESS and a cache with sampling AND summary timings.")
            valid_timing <-  with(data = runs_table, expr = is.finite(target_min_ESS) & target_min_ESS > 0 &
                                     is.finite(time_burnin) & time_burnin >= 0 & is.finite(time_sampling) & time_sampling >= 0 &
                                     is.finite(time_summaries) & time_summaries >= 0 & is.finite(min_ESS) & min_ESS > 0)
            runs_table <-  runs_table %>%
                dplyr::mutate(time_to_target_ESS = dplyr::if_else(condition = valid_timing,
                    true = .data$time_burnin + (.data$time_sampling + .data$time_summaries) * .data$target_min_ESS / .data$min_ESS,
                    false = NA_real_))
        }
        outcome_columns <-  c(ESS_per_1000_grad = "ESS_per_grad_samp", ESS_per_sec = "ESS_per_sec_samp",
                              burnin_sec = "time_burnin", sampling_sec = "time_sampling")
        outcome_column <-  if (outcome %in% names(x = outcome_columns)) unname(obj = outcome_columns[outcome]) else outcome
        if (!outcome_column %in% names(x = runs_table) || !is.numeric(x = runs_table[[outcome_column]])) stop("Unknown/non-numeric outcome: ", outcome)
        runs_table <-  runs_table %>%
            dplyr::mutate(outcome_value = .data[[outcome_column]] * if (outcome == "ESS_per_1000_grad") 1000 else 1,
                          exclusion_reason = dplyr::case_when(
                              is.na(.data$readable) | !.data$readable ~ "unreadable fit",
                              !is.finite(.data$outcome_value) ~ "missing/non-finite outcome (or timing/ESS input)",
                              .data$outcome_value <= 0 ~ "nonpositive outcome cannot be log-transformed",
                              TRUE ~ NA_character_))
        model_data <-  runs_table %>% dplyr::filter(is.na(.data$exclusion_reason))
        excluded_runs <-  runs_table %>% dplyr::filter(!is.na(.data$exclusion_reason))
        if (!nrow(x = model_data)) stop("No readable positive outcomes remain for this regression.")
        ## Diagnostics are retained, not used to filter the efficiency results. Missing fits are never assigned a fictitious runtime.
        ##
        ## ---- Coverage of the intended Cartesian product, including incomplete repetitions ----------------------------------------------
        ##
        if (is.null(x = expected_factor_levels)) {
            expected_factor_levels <-  lapply(X = runs_table[factor_settings], FUN = unique)
        }
        if (!is.list(x = expected_factor_levels) || anyDuplicated(x = names(x = expected_factor_levels)) ||
            !setequal(x = names(x = expected_factor_levels), y = factor_settings) || any(lengths(x = expected_factor_levels) == 0)) {
            stop("expected_factor_levels must be a named list covering exactly factor_settings, with nonempty levels.")
        }
        expected_factor_levels <-  lapply(X = expected_factor_levels[factor_settings], FUN = function(values) unique(as.character(x = values)))
        for (setting_name in factor_settings) {
            if (anyNA(x = expected_factor_levels[[setting_name]]) ||
                !all(as.character(x = runs_table[[setting_name]]) %in% expected_factor_levels[[setting_name]])) {
                stop("Unexpected/missing planned levels for ", setting_name)
            }
        }
        expected_cells <-  do.call(what = tidyr::expand_grid, args = expected_factor_levels)
        ##
        ## ---- tau_sampling_scale is inert under jittered burn-in (the runner collapses it to "none" there), so the cells
        ##      randomize_tau_burnin = TRUE x tau_sampling_scale != "none" never exist by design: not "missing" fits.
        ##
        if (all(c("randomize_tau_burnin", "tau_sampling_scale") %in% names(x = expected_cells))) {
            structurally_absent_cells <-  expected_cells$randomize_tau_burnin == "TRUE" & expected_cells$tau_sampling_scale != "none"
            expected_cells <-  expected_cells[!structurally_absent_cells, , drop = FALSE]
        }
        ##
        ## ---- tau_adaptation_block ("_tbJ";, EXPERIMENTAL) is inert whenever tau is pinned (the runner collapses
        ##      it to "main" there), so the cells manual_L/manual_tau_value != NA x tau_adaptation_block = "joint" never
        ##      exist by design: not "missing" fits (mirrors the tau_sampling_scale collapse immediately above).
        ##
        if (all(c("manual_L", "tau_adaptation_block") %in% names(x = expected_cells))) {
            structurally_absent_cells <-  expected_cells$manual_L != "NA" & expected_cells$tau_adaptation_block == "joint"
            expected_cells <-  expected_cells[!structurally_absent_cells, , drop = FALSE]
        }
        if (all(c("manual_tau_value", "tau_adaptation_block") %in% names(x = expected_cells))) {
            structurally_absent_cells <-  expected_cells$manual_tau_value != "NA" & expected_cells$tau_adaptation_block == "joint"
            expected_cells <-  expected_cells[!structurally_absent_cells, , drop = FALSE]
        }
        observed_cells <-  runs_table %>%
            dplyr::mutate(dplyr::across(.cols = dplyr::all_of(x = factor_settings), .fns = as.character)) %>%
            dplyr::group_by(dplyr::across(.cols = dplyr::all_of(x = factor_settings))) %>%
            dplyr::summarise(n_runs_present = dplyr::n(), n_runs_analysed = sum(is.na(.data$exclusion_reason)), .groups = "drop")
        cell_counts <-  expected_cells %>%
            dplyr::left_join(y = observed_cells, by = factor_settings) %>%
            dplyr::mutate(n_runs_present = dplyr::coalesce(.data$n_runs_present, 0),
                          n_runs_analysed = dplyr::coalesce(.data$n_runs_analysed, 0),
                          n_runs_expected = length(x = expected_run_indices),
                          n_runs_missing = .data$n_runs_expected - .data$n_runs_present)
        coverage <-  tibble::tibble(outcome = outcome, configurations_present = nrow(x = observed_cells),
                                    configurations_expected = nrow(x = expected_cells), fits_present = nrow(x = runs_table),
                                    fits_expected = nrow(x = expected_cells) * length(x = expected_run_indices),
                                    fits_analysed = nrow(x = model_data), fits_excluded = nrow(x = excluded_runs))
        if (any(cell_counts$n_runs_missing > 0) || nrow(x = excluded_runs)) {
            warning("Incomplete grid/outcomes: effects depend on the fitted model in unobserved cells. Inspect cell_counts and excluded_runs.",
                    call. = FALSE)
        }
        ##
        ## ---- Set categorical reference levels explicitly; do not assume a linear learning-rate response ---------------------------------
        ##
        factor_level_rows <-  list()
        for (setting_name in factor_settings) {
            observed_levels <-  as.character(x = sort(x = unique(x = model_data[[setting_name]])))
            preferred_reference <-  as.character(x = reference_levels[[setting_name]])
            if (length(x = preferred_reference) && preferred_reference %in% observed_levels) {
                observed_levels <-  c(preferred_reference, setdiff(x = observed_levels, y = preferred_reference))
            }
            model_data[[setting_name]] <-  factor(x = as.character(x = model_data[[setting_name]]), levels = observed_levels)
            factor_level_rows[[setting_name]] <-  tibble::tibble(setting = setting_name, reference = observed_levels[1],
                                                                 levels = paste(observed_levels, collapse = ", "),
                                                                 n_levels = length(x = observed_levels))
        }
        factor_levels <-  dplyr::bind_rows(factor_level_rows)
        active_settings <-  factor_levels$setting[factor_levels$n_levels > 1]
        model_data <-  model_data %>% dplyr::mutate(seed_block = factor(x = .data$run), log_outcome = log(x = .data$outcome_value))
        additive_terms <-  c(active_settings, if (dplyr::n_distinct(model_data$run) > 1) "seed_block")
        if (!length(x = additive_terms)) stop("No varying settings or seed blocks to model.")
        ##
        ## ---- Fit OLS and calculate heteroskedasticity-consistent HC3 covariance without extra dependencies -------------------------------
        ##
        fn_fit_model <-  function(model_terms) {

                model_formula <-  stats::reformulate(termlabels = model_terms, response = "log_outcome")
                contrast_names <-  c(active_settings, if ("seed_block" %in% model_terms) "seed_block")
                contrasts <-  stats::setNames(object = rep(x = list("contr.treatment"), times = length(x = contrast_names)), nm = contrast_names)
                fitted_model <-  stats::lm(formula = model_formula, data = model_data, contrasts = contrasts, na.action = stats::na.fail,
                                          x = TRUE, y = TRUE)
                model_matrix <-  stats::model.matrix(object = fitted_model)
                if (fitted_model$rank < ncol(x = model_matrix) || stats::df.residual(object = fitted_model) < 1) {
                    return(list(model = fitted_model, status = "not identifiable: aliased terms or no residual degrees of freedom"))
                }
                leverage <-  stats::hatvalues(model = fitted_model)
                if (any(1 - leverage < 1e-8)) return(list(model = fitted_model, status = "HC3 unavailable: a fitted row has leverage approximately one"))
                residuals <-  stats::residuals(object = fitted_model)
                inverse_crossproduct <-  summary(object = fitted_model)$cov.unscaled
                inverse_crossproduct <-  inverse_crossproduct[colnames(x = model_matrix), colnames(x = model_matrix), drop = FALSE]
                weighted_matrix <-  model_matrix * as.numeric(x = residuals / (1 - leverage))
                covariance_HC3 <-  inverse_crossproduct %*% crossprod(x = weighted_matrix) %*% inverse_crossproduct
                coefficient_estimates <-  stats::coef(object = fitted_model)
                standard_errors <-  sqrt(x = pmax(diag(x = covariance_HC3), 0))
                critical_value <-  stats::qt(p = (1 + confidence_level) / 2, df = stats::df.residual(object = fitted_model))
                coefficients <-  tibble::tibble(term = names(x = coefficient_estimates), log_ratio = unname(obj = coefficient_estimates),
                                                standard_error_HC3 = unname(obj = standard_errors)) %>%
                    dplyr::mutate(ratio = exp(x = .data$log_ratio),
                                  ratio_lower = exp(x = .data$log_ratio - critical_value * .data$standard_error_HC3),
                                  ratio_upper = exp(x = .data$log_ratio + critical_value * .data$standard_error_HC3),
                                  percent_change = 100 * (.data$ratio - 1),
                                  interpretation = dplyr::case_when(
                                      .data$term == "(Intercept)" ~ "baseline geometric outcome, not a ratio",
                                      startsWith(x = .data$term, prefix = "seed_block") ~ "seed block adjustment",
                                      grepl(pattern = ":", x = .data$term) ~ "ratio of ratios",
                                      TRUE ~ "level / reference (conditional when interactions are included)"),
                                  dplyr::across(.cols = dplyr::all_of(x = c("ratio", "ratio_lower", "ratio_upper", "percent_change")),
                                                .fns = function(values) dplyr::if_else(condition = .data$term == "(Intercept)",
                                                                                     true = NA_real_, false = values)))
                residual_diagnostics <-  model_data %>%
                    dplyr::select(dplyr::any_of(x = c("file", "run", factor_settings, "outcome_value", "sampling_acceptance_probability",
                                                       "pct_divs", "max_Rhat", "max_nRhat"))) %>%
                    dplyr::mutate(fitted_log_outcome = stats::fitted(object = fitted_model), residual_log_outcome = residuals,
                                  leverage = leverage, cooks_distance = stats::cooks.distance(model = fitted_model))
                return(list(model = fitted_model, status = "fitted", coefficients = coefficients, covariance_HC3 = covariance_HC3,
                            baseline_geometric_outcome = exp(x = unname(obj = coefficient_estimates["(Intercept)"])),
                            residual_diagnostics = residual_diagnostics))

        }
        additive <-  fn_fit_model(model_terms = additive_terms)
        if (additive$status != "fitted") stop("Additive model ", additive$status, ". Reduce/filter the requested design.")
        ## Every two-factor cell must be observed before estimating its full interaction. This is necessary, not sufficient:
        ## the joint model's matrix rank is still checked below. Unsupported terms are listed explicitly, never extrapolated silently.
        interaction_coverage <-  dplyr::bind_rows(tibble::tibble(term = character(), cells_observed = integer(),
                                                                 cells_expected = numeric(), included = logical(), reason = character()),
                                                 lapply(X = interaction_pairs, FUN = function(pair) {
            observed_pair_cells <-  nrow(x = dplyr::distinct(.data = dplyr::select(.data = model_data, dplyr::all_of(x = pair))))
            expected_pair_cells <-  prod(vapply(X = model_data[pair], FUN = nlevels, FUN.VALUE = 0))
            pair_varies <-  all(pair %in% active_settings)
            tibble::tibble(term = paste(pair, collapse = ":"), cells_observed = observed_pair_cells, cells_expected = expected_pair_cells,
                          included = pair_varies && observed_pair_cells == expected_pair_cells,
                          reason = if (!pair_varies) "one or both settings constant in analysed fits" else
                              if (observed_pair_cells < expected_pair_cells) "missing two-factor cells" else "included")
        }))
        interaction_terms <-  unique(x = interaction_coverage$term[interaction_coverage$included])
        if (nrow(x = interaction_coverage) && any(interaction_coverage$cells_observed < interaction_coverage$cells_expected)) {
            warning("Not yet estimating interactions with missing two-factor cells: ",
                    paste(interaction_coverage$term[interaction_coverage$cells_observed < interaction_coverage$cells_expected], collapse = ", "),
                    ". See interaction_coverage; they enter automatically once those cells are available.", call. = FALSE)
        }
        interactions <-  if (length(x = interaction_terms)) fn_fit_model(model_terms = c(additive_terms, interaction_terms)) else {
            list(model = NULL, status = "not requested or the relevant settings are constant")
        }
        if (length(x = interaction_terms) && interactions$status != "fitted") {
            warning("Interaction model ", interactions$status, ". No interaction effect table is reported; do not interpret aliased coefficients.",
                    call. = FALSE)
        }
        ## An unavailable interaction model still returns a printable coefficient table, never NULL.
        if (is.null(interactions$coefficients)) interactions$coefficients <-  additive$coefficients[0, , drop = FALSE]
        ##
        ## ---- All pairwise setting contrasts, at every combination of their fitted interaction partners ---------------------------------
        ##
        additive$conditional_contrasts <-  fn_ps7_all_regression_contrasts(
            fitted_result = additive, factor_settings = active_settings, confidence_level = confidence_level)
        interactions$conditional_contrasts <-  fn_ps7_all_regression_contrasts(
            fitted_result = interactions, factor_settings = active_settings, confidence_level = confidence_level)
        ##
        ## ---- Exact same-seed pairing; report each seed separately instead of pretending to have hundreds of independent seeds -------------
        ##
        matched_run_tables <-  list()
        if (!is.null(x = matched_setting) && (!is.character(x = matched_setting) || anyNA(x = matched_setting) ||
            anyDuplicated(x = matched_setting) || !all(matched_setting %in% setting_names))) {
            stop("matched_setting must be NULL or distinct parsed setting names.")
        }
        skipped_matched_settings <-  tibble::tibble(setting = character(), reason = character())
        for (comparison_setting in matched_setting) {
            if (!comparison_setting %in% factor_settings || nlevels(x = model_data[[comparison_setting]]) != 2) {
                skipped_matched_settings <-  dplyr::bind_rows(skipped_matched_settings,
                    tibble::tibble(setting = comparison_setting,
                                   reason = if (!comparison_setting %in% factor_settings) "not selected as a regression factor" else
                                       "requires exactly two observed levels"))
                next
            }
            comparison_levels <-  levels(x = model_data[[comparison_setting]])
            matching_columns <-  c(setdiff(x = factor_settings, y = comparison_setting), "run")
            reference_runs <-  model_data %>% dplyr::filter(.data[[comparison_setting]] == comparison_levels[1]) %>%
                dplyr::select(dplyr::all_of(x = matching_columns), reference_outcome = "outcome_value", reference_file = "file")
            alternative_runs <-  model_data %>% dplyr::filter(.data[[comparison_setting]] == comparison_levels[2]) %>%
                dplyr::select(dplyr::all_of(x = matching_columns), alternative_outcome = "outcome_value", alternative_file = "file")
            matched_run_tables[[comparison_setting]] <-  dplyr::inner_join(x = reference_runs, y = alternative_runs, by = matching_columns) %>%
                dplyr::mutate(matched_setting = comparison_setting,
                              comparison = paste(comparison_levels[2], "/", comparison_levels[1]),
                              ratio = .data$alternative_outcome / .data$reference_outcome, log_ratio = log(x = .data$ratio))
        }
        matched_runs <-  dplyr::bind_rows(matched_run_tables)
        matched_by_seed <-  tibble::tibble(matched_setting = character(), run = numeric(), comparison = character(),
                                           n_pairs = integer(), geometric_mean_ratio = numeric())
        matched_summary <-  tibble::tibble(matched_setting = character(), comparison = character(), n_pairs = integer(),
                                           n_seeds = integer(), geometric_mean_ratio = numeric(), n_alternative_lower = integer())
        if (nrow(x = matched_runs)) {
                matched_by_seed <-  matched_runs %>% dplyr::group_by(.data$matched_setting, .data$run, .data$comparison) %>%
                    dplyr::summarise(n_pairs = dplyr::n(), geometric_mean_ratio = exp(x = mean(x = .data$log_ratio)), .groups = "drop")
                matched_summary <-  matched_runs %>% dplyr::group_by(.data$matched_setting, .data$comparison) %>%
                    dplyr::summarise(n_pairs = dplyr::n(), n_seeds = dplyr::n_distinct(.data$run),
                                     geometric_mean_ratio = exp(x = mean(x = .data$log_ratio)),
                                     n_alternative_lower = sum(.data$ratio < 1), .groups = "drop")
        }
        if (isTRUE(x = print_table)) {
            cat("\n==== PS7 FACTORIAL REGRESSION: ", outcome, " ====\n", sep = "")
            print(x = coverage, width = Inf)
            if (length(omitted_default_interaction_pairs)) {
                cat("Automatic interactions skipped because their factors were not selected: ",
                    paste(omitted_default_interaction_pairs, collapse = ", "), ".\n", sep = "")
            }
            if (nrow(skipped_matched_settings)) {
                cat("\nRequested matched comparisons skipped:\n")
                print(x = skipped_matched_settings, n = Inf, width = Inf)
            }
            cat("Categorical settings; run index is a fixed seed block. Diagnostics do not filter efficiency results.\n")
            print(x = factor_levels, n = Inf, width = Inf)
            if (nrow(x = dependent_settings)) {
                cat("\nCentre endpoints describe the centre-rule contrast below; they are not separate regression terms.\n")
                print(x = dependent_settings, n = Inf, width = Inf)
            }
            cat("\nAdditive-model ratios: level / reference. For times, below 1 is faster; for ESS rates, above 1 is better.\n")
            cat("Intervals are approximate ", confidence_level * 100, "% HC3 intervals, not seed-clustered or selection-adjusted.\n", sep = "")
            print(x = additive$coefficients %>%
                      dplyr::filter(.data$term != "(Intercept)", !startsWith(x = .data$term, prefix = "seed_block")) %>%
                      dplyr::select(dplyr::all_of(x = c("term", "ratio", "ratio_lower", "ratio_upper", "percent_change"))), n = Inf, width = Inf)
            cat("\nInteraction model: ", interactions$status, ". Full coefficients and residual diagnostics are returned.\n", sep = "")
            if (nrow(x = interaction_coverage)) print(x = interaction_coverage, n = Inf, width = Inf)
            if (nrow(x = interactions$conditional_contrasts)) {
                cat("\nAll pairwise setting contrasts from the interaction model (full HC3 covariance):\n")
                cat("Conditioning covers every combination of fitted interaction partners; none means no fitted interaction for that setting.\n")
                cat("Other predictors held at reference levels. Seed blocks are adjustment terms, not settings compared here.\n")
                cat("Intervals are pointwise, not multiplicity-adjusted. Observed support does not imply exact matching on other settings.\n")
                print(x = interactions$conditional_contrasts %>%
                          dplyr::select(dplyr::all_of(x = c("comparison_setting", "comparison", "at", "ratio", "ratio_lower", "ratio_upper", "percent_change", "observed_support"))),
                      n = Inf, width = Inf)
            } else {
                cat("\nAll pairwise setting contrasts from the additive model (no fitted interaction model):\n")
                print(x = additive$conditional_contrasts %>%
                          dplyr::select(dplyr::all_of(x = c("comparison_setting", "comparison", "ratio", "ratio_lower", "ratio_upper", "percent_change"))),
                      n = Inf, width = Inf)
            }
            if (nrow(x = matched_by_seed)) {
                cat("\nExact same-seed comparisons, matched on all other factors (observed pairs only):\n")
                print(x = matched_by_seed, n = Inf, width = Inf)
            }
            if (nrow(x = excluded_runs)) cat("\nExcluded outcomes are retained in $excluded_runs with reasons.\n")
        }
        return(list(outcome = outcome, coverage = coverage, factor_levels = factor_levels, dependent_settings = dependent_settings,
                    constant_settings = dplyr::slice(.data = dplyr::select(.data = runs_table, dplyr::all_of(x = fixed_setting_names)), 1),
                    cell_counts = cell_counts, missing_cells = dplyr::filter(.data = cell_counts, .data$n_runs_missing > 0),
                    model_data = model_data, excluded_runs = excluded_runs, additive = additive, interactions = interactions,
                    interaction_coverage = interaction_coverage,
                    omitted_default_interaction_pairs = omitted_default_interaction_pairs,
                    skipped_matched_settings = skipped_matched_settings,
                    matched_runs = matched_runs, matched_by_seed = matched_by_seed, matched_summary = matched_summary,
                    confidence_level = confidence_level))

}
##
## ---- Contrast two factor levels at specified settings; remaining predictors use their fitted reference levels -------------------------
##
fn_ps7_regression_contrast <-  function( fitted_result,
                                         comparison_setting,
                                         reference,
                                         alternative,
                                         at = list(),
                                         confidence_level = 0.95
) {
        if (!identical(fitted_result$status, "fitted")) stop("The requested model must be identifiable with HC3 covariance available.")
        if (!is.numeric(confidence_level) || length(confidence_level) != 1 || !is.finite(confidence_level) ||
            confidence_level <= 0 || confidence_level >= 1) stop("confidence_level must be between zero and one.")
        fitted_model <-  fitted_result$model
        factor_levels <-  fitted_model$xlevels
        if (length(comparison_setting) != 1 || !comparison_setting %in% names(factor_levels)) stop("comparison_setting must name a fitted factor.")
        if (length(reference) != 1 || length(alternative) != 1 ||
            !all(c(reference, alternative) %in% factor_levels[[comparison_setting]]) || reference == alternative) {
            stop("reference and alternative must be distinct observed levels of comparison_setting.")
        }
        if (!is.list(at) || (length(at) && (is.null(names(at)) || anyDuplicated(names(at)) ||
            any(!names(at) %in% names(factor_levels)) || comparison_setting %in% names(at)))) {
            stop("at must be a named list of other fitted factors.")
        }
        prediction_rows <-  tibble::as_tibble(stats::model.frame(fitted_model)[rep(1, 2), , drop = FALSE])
        for (setting_name in names(factor_levels)) {
            selected_level <-  if (setting_name %in% names(at)) as.character(at[[setting_name]]) else factor_levels[[setting_name]][1]
            if (length(selected_level) != 1 || is.na(selected_level) || !selected_level %in% factor_levels[[setting_name]]) {
                stop("Unknown or nonscalar level for ", setting_name, ".")
            }
            prediction_rows[[setting_name]] <-  factor(rep(selected_level, 2), levels = factor_levels[[setting_name]])
        }
        prediction_rows[[comparison_setting]] <-  factor(c(reference, alternative), levels = factor_levels[[comparison_setting]])
        prediction_matrix <-  stats::model.matrix(stats::delete.response(stats::terms(fitted_model)),
                                                   data = prediction_rows, contrasts.arg = fitted_model$contrasts,
                                                   xlev = factor_levels)
        coefficient_names <-  names(stats::coef(fitted_model))
        contrast_vector <-  prediction_matrix[2, coefficient_names] - prediction_matrix[1, coefficient_names]
        covariance_HC3 <-  fitted_result$covariance_HC3[coefficient_names, coefficient_names, drop = FALSE]
        log_ratio <-  sum(contrast_vector * stats::coef(fitted_model))
        standard_error_HC3 <-  sqrt(max(0, as.numeric(t(contrast_vector) %*% covariance_HC3 %*% contrast_vector)))
        critical_value <-  stats::qt((1 + confidence_level) / 2, df = stats::df.residual(fitted_model))
        remaining_settings <-  setdiff(names(factor_levels), c(comparison_setting, names(at)))
        return(tibble::tibble(
            comparison = paste(alternative, "/", reference),
            at = if (length(at)) paste(paste(names(at), unlist(at), sep = " = "), collapse = " | ") else "none",
            other_reference_settings = paste(vapply(remaining_settings, function(setting_name) {
                paste(setting_name, factor_levels[[setting_name]][1], sep = " = ")
            }, FUN.VALUE = ""), collapse = " | "),
            log_ratio = log_ratio, standard_error_HC3 = standard_error_HC3,
            ratio = exp(log_ratio), ratio_lower = exp(log_ratio - critical_value * standard_error_HC3),
            ratio_upper = exp(log_ratio + critical_value * standard_error_HC3),
            percent_change = 100 * (exp(log_ratio) - 1), confidence_level = confidence_level))
}
##
## ---- Generate all unordered level pairs for each setting, including non-reference comparisons ----------------------------------------
##
fn_ps7_all_regression_contrasts <-  function( fitted_result,
                                             factor_settings,
                                             confidence_level = 0.95
) {
        if (!identical(fitted_result$status, "fitted")) return(tibble::tibble())
        fitted_model <-  fitted_result$model
        model_frame <-  stats::model.frame(fitted_model)
        term_factors <-  attr(stats::terms(fitted_model), "factors")
        contrast_tables <-  list()
        for (comparison_setting in factor_settings) {
            comparison_levels <-  fitted_model$xlevels[[comparison_setting]]
            if (length(comparison_levels) < 2) next
            involving_terms <-  which(term_factors[comparison_setting, ] > 0)
            partner_settings <-  setdiff(rownames(term_factors)[
                rowSums(term_factors[, involving_terms, drop = FALSE] > 0) > 0], comparison_setting)
            conditioning_grid <-  if (length(partner_settings)) {
                do.call(tidyr::expand_grid, fitted_model$xlevels[partner_settings])
            } else tibble::tibble(.unconditional = TRUE)
            level_pairs <-  utils::combn(comparison_levels, 2, simplify = FALSE)
            for (level_pair in level_pairs) {
                for (conditioning_index in seq_len(nrow(conditioning_grid))) {
                    conditioning_values <-  if (length(partner_settings)) as.list(conditioning_grid[conditioning_index, ]) else list()
                    contrast_table <-  fn_ps7_regression_contrast(
                        fitted_result = fitted_result,
                        comparison_setting = comparison_setting,
                        reference = level_pair[1],
                        alternative = level_pair[2],
                        at = conditioning_values,
                        confidence_level = confidence_level)
                    supporting_rows <-  rep(TRUE, nrow(model_frame))
                    for (partner_setting in partner_settings) {
                        supporting_rows <-  supporting_rows &
                            as.character(model_frame[[partner_setting]]) == conditioning_values[[partner_setting]]
                    }
                    n_reference <-  sum(supporting_rows & as.character(model_frame[[comparison_setting]]) == level_pair[1])
                    n_alternative <-  sum(supporting_rows & as.character(model_frame[[comparison_setting]]) == level_pair[2])
                    contrast_tables[[length(contrast_tables) + 1]] <-  contrast_table %>%
                        dplyr::mutate(comparison_setting = comparison_setting,
                                      reference_level = level_pair[1], alternative_level = level_pair[2],
                                      n_reference = n_reference, n_alternative = n_alternative,
                                      observed_support = n_reference > 0 && n_alternative > 0, .before = 1)
                }
            }
        }
        return(dplyr::bind_rows(contrast_tables))
}






















