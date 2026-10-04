
## fn_ps7_summarise_runs.R
##
## Summarise saved ps7 runs: one row per CONFIGURATION (or per selected group), sorted by any metric.
##
## A configuration is the run's file name without "_runK": the file-name builder encodes every setting, so runs
## that differ in ANY setting - even one fn_ps7_parse_run_name() does not know about - are never pooled.
## Works straight from a directory, so it replaces summarize_ps7_results(), which had to rebuild every file
## name from a settings grid (and found nothing whenever one grid field was off).
##
## Every number is a mean (or median, for R-hats) over ALL readable runs - nothing is blanked out. Convergence is
## shown next to the numbers (nRhat_ok / Rhat_ok counts), never used to hide them.
##
## Metrics per run:
##   min_ESS, max_Rhat, max_nRhat          - over the parameters of interest (Se / Sp / prev)
##   min_ESS_main, max_Rhat_main, ...      - over the MAIN (raw) parameters
##   time_to_target                        - burn-in + (sampling + summaries) scaled by target ESS / observed ESS
##   Rhat_floor                            - sqrt(1 + 2 * n_chains_sampling / min_ESS): the split R-hat that PERFECTLY
##                                           mixing chains of this length would give. With many short chains R-hat can
##                                           never go below it, so R-hat near its floor means "chains too short for
##                                           R-hat", not "not converged" - nested R-hat is the meaningful one there.
##
## Usage:
##   fn_ps7_summarise_runs(runs_table)                                             ## defaults
##   fn_ps7_summarise_runs(runs_table, sort_by = "time_to_target_ESS")             ## fastest to the target first
##   fn_ps7_summarise_runs(runs_table, sort_by = c("nRhat_ok", "ESS_grad"))        ## tie-breaks, left to right
##   fn_ps7_summarise_runs(runs_table, target_min_ESS = 2000)                      ## one target for every N
##   fn_ps7_summarise_runs(runs_table, target_min_ESS = c("2500" = 2500, "10000" = 1000))   ## per N
##   fn_ps7_summarise_runs(runs_table, group_by = c("diffusion_HMC_integrator", "learning_rate"))  ## pool the rest
##   fn_ps7_summarise_runs(runs_table[runs_table$diffusion_HMC_integrator == "kick_flow_kick", ])   ## subset first
##   fn_ps7_print_metric_names()                                                   ## every sortable / printable metric
##


##
## ---- every metric the summary can sort by or print. short_name is passed to sort_by / columns_to_print
##      (the full column_name works too); bigger_is_better sets the default sort direction:
##
R_fn_ps7_metric_specification <-  function() {

        metric_specification <-  tibble::tibble(
              short_name       = c("runs", "nRhat_ok", "Rhat_ok",
                                   "ESS_grad", "ESS_sec", "min_ESS", "min_ESS_main",
                                   "burnin_sec", "sampling_sec", "target_ESS", "time_to_target_ESS",
                                   "Rhat", "Rhat_floor", "nRhat", "Rhat_main", "nRhat_main",
                                   "divs_pct", "eps", "L", "grads_per_s", "accept_prob"),
              column_name      = c("n_runs", "n_nRhat_ok", "n_Rhat_ok",
                                   "mean_ESS_per_grad_x1000", "mean_ESS_per_sec", "mean_min_ESS", "mean_min_ESS_main",
                                   "mean_time_burnin", "mean_time_sampling", "target_min_ESS", "mean_time_to_target_min_ESS",
                                   "median_max_Rhat", "median_Rhat_floor", "median_max_nRhat", "median_max_Rhat_main", "median_max_nRhat_main",
                                   "mean_pct_divs", "mean_eps", "mean_L", "mean_grad_evals_per_sec", "mean_sampling_acceptance_probability"),
              bigger_is_better = c(TRUE, TRUE, TRUE,
                                   TRUE, TRUE, TRUE, TRUE,
                                   FALSE, FALSE, FALSE, FALSE,
                                   FALSE, FALSE, FALSE, FALSE, FALSE,
                                   FALSE, TRUE, FALSE, TRUE, TRUE),
              number_format    = c("%.0f", "", "",
                                   "", "%.0f", "%.0f", "%.0f",
                                   "%.1f", "%.1f", "%.0f", "",
                                   "%.3f", "%.3f", "%.3f", "%.3f", "%.3f",
                                   "%.2f", "%.3f", "%.1f", "%.0f", "%.3f"),
              description      = c("runs in the row",
                                   "runs with nested R-hat (Se/Sp/prev) <= max_nRhat_ok and divergences <= max_pct_divs_ok",
                                   "runs with split R-hat (Se/Sp/prev) <= max_Rhat_ok and divergences <= max_pct_divs_ok",
                                   "min ESS (Se/Sp/prev) per 1000 gradients, mean [min-max] over runs",
                                   "min ESS (Se/Sp/prev) per second of sampling",
                                   "min ESS over Se/Sp/prev",
                                   "min ESS over the MAIN (raw) parameters",
                                   "burn-in time (sec)",
                                   "sampling time (sec)",
                                   "target min ESS used for time_to_target_ESS",
                                   "time (sec): burn-in + (sampling + summaries) * target ESS / observed ESS, mean [min-max] over runs",
                                   "median of the max split R-hat (Se/Sp/prev)",
                                   "median R-hat floor from chain length alone",
                                   "median of the max nested R-hat (Se/Sp/prev)",
                                   "median of the max split R-hat (MAIN parameters)",
                                   "median of the max nested R-hat (MAIN parameters)",
                                   "per cent divergent transitions",
                                   "step size (sampling)",
                                   "leapfrog steps per iteration (sampling)",
                                   "gradient evaluations per second (sampling)",
                                   "mean post-burn-in main/joint acceptance probability (0-1), over runs with recorded valid probabilities"))
        ##
        return(metric_specification)

}


##
## ---- friendlier spellings accepted by sort_by / columns_to_print:
##
R_fn_ps7_metric_name_aliases <-  function() {

        return(c( ESS_per_grad       = "ESS_grad",
                  ESS_per_sec        = "ESS_sec",
                  time_to_target     = "time_to_target_ESS",
                  to_target_s        = "time_to_target_ESS",
                  time_burnin        = "burnin_sec",
                  burnin_time        = "burnin_sec",
                  burnin_s           = "burnin_sec",
                  time_sampling      = "sampling_sec",
                  sampling_time      = "sampling_sec",
                  sampling_s         = "sampling_sec",
                  pct_divs           = "divs_pct",
                  divergences        = "divs_pct",
                  max_Rhat           = "Rhat",
                  max_nRhat          = "nRhat",
                  nested_Rhat        = "nRhat",
                  grad_evals_per_sec = "grads_per_s",
                  sampling_acceptance_probability = "accept_prob",
                  sampling_acceptance = "accept_prob"))

}


##
## ---- print every metric name available for sorting / printing:
##
fn_ps7_print_metric_names <-  function() {

        metric_specification <-  R_fn_ps7_metric_specification()
        print(x = dplyr::select(.data = metric_specification, short_name, column_name, bigger_is_better, description), n = Inf, width = Inf)
        metric_name_aliases <-  R_fn_ps7_metric_name_aliases()
        cat("\nalso accepted:", paste0(names(metric_name_aliases), " (= ", metric_name_aliases, ")", collapse = ", "), "\n")
        return(invisible(metric_specification))

}


##
## ---- turn user-facing metric names (short names, aliases or full column names) into rows of the specification:
##
R_fn_ps7_resolve_metric_names <-  function( requested_metric_names,
                                           argument_name) {

        metric_specification <-  R_fn_ps7_metric_specification()
        metric_name_aliases  <-  R_fn_ps7_metric_name_aliases()
        ##
        resolved_short_names <-  sapply(requested_metric_names, function(requested_metric_name) {
              if (requested_metric_name %in% names(metric_name_aliases))         return(unname(metric_name_aliases[requested_metric_name]))
              if (requested_metric_name %in% metric_specification$short_name)   return(requested_metric_name)
              if (requested_metric_name %in% metric_specification$column_name)  return(metric_specification$short_name[metric_specification$column_name == requested_metric_name])
              return(NA_character_)
        })
        ##
        if (anyNA(resolved_short_names)) {
              stop(argument_name, ": unknown metric(s) ", paste(requested_metric_names[is.na(resolved_short_names)], collapse = ", "),
                   ". Valid: ", paste(metric_specification$short_name, collapse = ", "),
                   " (or run fn_ps7_print_metric_names()).")
        }
        ##
        return(dplyr::slice(.data = metric_specification, match(resolved_short_names, metric_specification$short_name)))

}


##
## ---- target min ESS for each run. target_min_ESS may be:
##        NULL                                    -> the ps3 targets per (N, model type), via R_fn_ps7_get_target_min_ESS
##        one number                              -> that target for every N
##        a named vector, e.g. c("2500" = 2500, "10000" = 1000)  -> per N (names are the N values)
##        a function(N, Model_type)               -> a custom rule
##      Stops if any N in the table ends up without a target (a silent NA would blank the time-to-target column):
##
R_fn_ps7_resolve_target_min_ESS <-  function( N_of_each_run,
                                             model_type_of_each_run,
                                             target_min_ESS) {

        if (is.null(target_min_ESS)) {
              target_min_ESS_of_each_run <-  mapply(function(N_of_run, model_type_of_run) R_fn_ps7_get_target_min_ESS(N = N_of_run, Model_type = model_type_of_run),
                                                   N_of_each_run,
                                                   model_type_of_each_run)
        } else if (is.function(target_min_ESS)) {
              target_min_ESS_of_each_run <-  mapply(function(N_of_run, model_type_of_run) target_min_ESS(N = N_of_run, Model_type = model_type_of_run),
                                                   N_of_each_run,
                                                   model_type_of_each_run)
        } else if (is.numeric(target_min_ESS) && length(target_min_ESS) == 1 && is.null(names(target_min_ESS))) {
              target_min_ESS_of_each_run <-  rep(target_min_ESS, length(N_of_each_run))
        } else if (is.numeric(target_min_ESS) && !is.null(names(target_min_ESS))) {
              target_min_ESS_of_each_run <-  unname(target_min_ESS[as.character(N_of_each_run)])
        } else {
              stop("target_min_ESS must be NULL (ps3 targets), one number, a vector named by N (e.g. c(\"2500\" = 2500, \"10000\" = 1000)) or a function(N, Model_type).")
        }
        ##
        target_min_ESS_of_each_run <-  as.numeric(unlist(target_min_ESS_of_each_run))
        if (length(target_min_ESS_of_each_run) != length(N_of_each_run)) stop("target_min_ESS must resolve to one target per run.")
        if (any(is.finite(target_min_ESS_of_each_run) & target_min_ESS_of_each_run <= 0)) stop("target_min_ESS must be positive.")
        N_without_target <-  unique(N_of_each_run[!is.finite(target_min_ESS_of_each_run)])
        if (length(N_without_target) > 0) {
              stop("target_min_ESS: no target for N = ", paste(N_without_target, collapse = ", "),
                   ". Pass one number, or a vector named by N that covers every N in the table.")
        }
        ##
        return(target_min_ESS_of_each_run)

}


##
## ---- the setting columns fn_ps7_parse_run_name() produces (everything except the run index, the key and the file name):
##
R_fn_ps7_setting_column_names <-  function(runs_table) {

        parsed_column_names <-  names(fn_ps7_parse_run_name(x = runs_table$file[1]))
        return(setdiff(parsed_column_names, c("run", "configuration", "file")))

}


##
## ---- summarise a runs table (from fn_ps7_extract_tau_sweep) into one row per configuration (or per group_by group):
##
#' @param runs_table        one row per run, from fn_ps7_extract_tau_sweep(). Subset it first to look at part of a sweep.
#' @param target_min_ESS    NULL = ps3 targets per N; one number = same for every N; vector named by N; or function(N, Model_type).
#' @param sort_by           metric(s) to sort by, left to right = tie-breaks. Short names ("ESS_grad", "time_to_target_ESS",
#'                          "burnin_s", "ESS_sec", "nRhat", "min_ESS_main", ...), aliases ("time_to_target", "time_burnin")
#'                          or full column names. fn_ps7_print_metric_names() lists them all.
#' @param sort_decreasing   NULL = best first (bigger first for ESS-type metrics, smaller first for times / R-hats / divs);
#'                          or TRUE / FALSE, one per sort_by entry.
#' @param sort_within_N     TRUE = rows grouped by N first (times and ESS targets are not comparable across N).
#' @param group_by          NULL = one row per configuration, labelled by every setting that varies. Or setting column
#'                          names, e.g. c("diffusion_HMC_integrator", "learning_rate"): one row per combination of those,
#'                          POOLING runs over every other setting (the pooled settings are printed).
#' @param columns_to_print  metrics shown in the printed table (the returned table always has all of them).
#' @param top_n             print only the best top_n rows (per N when sort_within_N = TRUE). Inf = all.
#' @param max_Rhat_ok, max_nRhat_ok, max_pct_divs_ok   thresholds for the Rhat_ok / nRhat_ok counts only.
#' @return invisibly, list(configurations, runs, label_settings, pooled_settings, constant_settings).
fn_ps7_summarise_runs <-  function( runs_table,
                                   target_min_ESS    = NULL,
                                   sort_by           = "ESS_grad",
                                   sort_decreasing   = NULL,
                                   sort_within_N     = TRUE,
                                   group_by          = NULL,
                                   columns_to_print  = c("runs", "nRhat_ok", "ESS_grad", "ESS_sec", "min_ESS", "min_ESS_main",
                                                         "burnin_sec", "time_to_target_ESS", "Rhat", "Rhat_floor", "nRhat",
                                                         "divs_pct", "eps", "L", "accept_prob"),
                                   top_n             = Inf,
                                   max_Rhat_ok       = 1.05,
                                   max_nRhat_ok      = 1.05,
                                   max_pct_divs_ok   = 1.0,
                                   print_table       = TRUE) {

        if (nrow(runs_table) == 0) stop("fn_ps7_summarise_runs: runs_table has no rows.")
        runs_table <-  tibble::as_tibble(x = runs_table)
        ##
        ## ---- a table made by an OLDER fn_ps7_extract_tau_sweep() (e.g. still loaded in the R session) lacks these columns:
        ##
        {
              required_column_names <-  c("configuration", "model_type", "N", "n_chains_burnin", "num_chunks_burnin", "file",
                                         "min_ESS_main", "max_Rhat_main", "max_nRhat_main")
              missing_column_names  <-  setdiff(required_column_names, names(runs_table))
              if (length(missing_column_names) > 0) {
                    stop("fn_ps7_summarise_runs: runs_table is missing column(s) ", paste(missing_column_names, collapse = ", "),
                         " - it was made by an older fn_ps7_extract_tau_sweep(). Re-source functions/fn_ps7_extract_tau_sweep.R",
                         " and re-run fn_ps7_extract_tau_sweep().")
              }
        }
        ##
        if (!("sampling_acceptance_probability" %in% names(x = runs_table))) {
              runs_table$sampling_acceptance_probability <-  NA_real_
              warning("Sampling acceptance is absent from this cached runs table. Re-source fn_ps7_extract_tau_sweep.R and ",
                      "re-extract the saved runs to populate accept_prob; no model fits are needed.", call. = FALSE)
        }
        ##
        ## ---- check the metric names up front, so a typo fails before any work:
        ##
        {
              sort_by_specification          <-  R_fn_ps7_resolve_metric_names(sort_by, "sort_by")
              columns_to_print_specification <-  R_fn_ps7_resolve_metric_names(columns_to_print, "columns_to_print")
              if (is.null(sort_decreasing)) sort_decreasing <-  sort_by_specification$bigger_is_better
              if (length(sort_decreasing) != nrow(sort_by_specification)) stop("sort_decreasing must be NULL or have one TRUE/FALSE per sort_by entry.")
        }
        ##
        ## ---- per run: target min ESS, time to it, R-hat floor, and the R-hat / nested R-hat pass flags:
        ##
        {
              runs_table$target_min_ESS <-  R_fn_ps7_resolve_target_min_ESS( N_of_each_run          = runs_table$N,
                                                                            model_type_of_each_run = runs_table$model_type,
                                                                            target_min_ESS         = target_min_ESS)
              if (!("time_summaries" %in% names(runs_table))) runs_table$time_summaries <-  NA_real_
              ##
              valid_target_timing <-  is.finite(runs_table$time_burnin) & runs_table$time_burnin >= 0 &
                                      is.finite(runs_table$time_sampling) & runs_table$time_sampling >= 0 &
                                      is.finite(runs_table$time_summaries) & runs_table$time_summaries >= 0 &
                                      is.finite(runs_table$min_ESS) & runs_table$min_ESS > 0
              ##
              runs_table$time_to_target_min_ESS <-  NA_real_
              runs_table$time_to_target_min_ESS[valid_target_timing] <-  runs_table$time_burnin[valid_target_timing] +
                  (runs_table$time_sampling[valid_target_timing] + runs_table$time_summaries[valid_target_timing]) *
                  runs_table$target_min_ESS[valid_target_timing] / runs_table$min_ESS[valid_target_timing]
              ##
              if (any(runs_table$readable & !valid_target_timing, na.rm = TRUE)) {
                    warning("Missing/invalid timing or ESS: time_to_target_ESS is NA for affected configurations. ",
                            "Re-extract saved runs if the cached table lacks time_summaries; missing timings are not treated as zero.")
              }
              runs_table$Rhat_floor_from_chain_length <-  sqrt(1 + (2 * runs_table$n_chains_sampling) / runs_table$min_ESS)
              ##
              divergences_ok <-  !is.finite(runs_table$pct_divs) | (runs_table$pct_divs <= max_pct_divs_ok)
              runs_table$Rhat_ok  <-  runs_table$readable & is.finite(runs_table$max_Rhat)  & (runs_table$max_Rhat  <= max_Rhat_ok)  & divergences_ok
              runs_table$nRhat_ok <-  runs_table$readable & is.finite(runs_table$max_nRhat) & (runs_table$max_nRhat <= max_nRhat_ok) & divergences_ok
        }
        ##
        ## ---- which settings vary, which label the rows, and which (if group_by is given) get pooled over:
        ##
        {
              setting_column_names <-  R_fn_ps7_setting_column_names(runs_table)
              n_distinct_values_per_setting <-  sapply(setting_column_names, function(setting_name) length(unique(as.character(runs_table[[setting_name]]))))
              varying_setting_names  <-  setting_column_names[n_distinct_values_per_setting > 1]
              constant_setting_names <-  setting_column_names[n_distinct_values_per_setting == 1]
              ##
              if (is.null(group_by)) {
                    label_setting_names <-  varying_setting_names
                    group_key_of_each_run <-  runs_table$configuration        ## the file-name key: never pools distinct settings
              } else {
                    unknown_group_by_names <-  setdiff(group_by, setting_column_names)
                    if (length(unknown_group_by_names) > 0) {
                          stop("group_by: unknown setting(s) ", paste(unknown_group_by_names, collapse = ", "),
                               ". Valid: ", paste(setting_column_names, collapse = ", "))
                    }
                    label_setting_names <-  unique(c(if (sort_within_N && "N" %in% varying_setting_names) "N", group_by))
                    group_key_of_each_run <-  do.call(paste, c(lapply(runs_table[label_setting_names], as.character), sep = " | "))
              }
              pooled_setting_names <-  setdiff(varying_setting_names, label_setting_names)
        }
        ##
        ## ---- one row per group:
        ##
        {
              fn_mean_or_NA   <-  function(values) if (sum(is.finite(values)) > 0) mean(values[is.finite(values)])   else NA_real_
              fn_median_or_NA <-  function(values) if (sum(is.finite(values)) > 0) median(values[is.finite(values)]) else NA_real_
              fn_min_or_NA    <-  function(values) if (sum(is.finite(values)) > 0) min(values[is.finite(values)])    else NA_real_
              fn_max_or_NA    <-  function(values) if (sum(is.finite(values)) > 0) max(values[is.finite(values)])    else NA_real_
              ##
              configuration_rows <-  lapply(split(runs_table, group_key_of_each_run), function(runs_of_group) {
                    readable_runs <-  dplyr::filter(.data = runs_of_group, .data$readable)
                    ## Do not compare a subset's target-time mean against all runs' burn-in mean.
                    complete_target_timings <-  nrow(readable_runs) > 0 && all(is.finite(readable_runs$time_to_target_min_ESS))
                    tibble::tibble( !!!dplyr::select(.data = dplyr::slice(.data = runs_of_group, 1), dplyr::all_of(label_setting_names)),
                                n_runs                      = nrow(runs_of_group),
                                n_nRhat_ok                  = sum(runs_of_group$nRhat_ok),
                                n_Rhat_ok                   = sum(runs_of_group$Rhat_ok),
                                mean_ESS_per_grad_x1000     = 1000 * fn_mean_or_NA(readable_runs$ESS_per_grad_samp),
                                min_ESS_per_grad_x1000      = 1000 * fn_min_or_NA(readable_runs$ESS_per_grad_samp),
                                max_ESS_per_grad_x1000      = 1000 * fn_max_or_NA(readable_runs$ESS_per_grad_samp),
                                mean_ESS_per_sec            = fn_mean_or_NA(readable_runs$ESS_per_sec_samp),
                                mean_min_ESS                = fn_mean_or_NA(readable_runs$min_ESS),
                                mean_min_ESS_main           = fn_mean_or_NA(readable_runs$min_ESS_main),
                                mean_time_burnin            = fn_mean_or_NA(readable_runs$time_burnin),
                                mean_time_sampling          = fn_mean_or_NA(readable_runs$time_sampling),
                                target_min_ESS              = fn_mean_or_NA(readable_runs$target_min_ESS),
                                mean_time_to_target_min_ESS = if (complete_target_timings) mean(readable_runs$time_to_target_min_ESS) else NA_real_,
                                min_time_to_target_min_ESS  = if (complete_target_timings) min(readable_runs$time_to_target_min_ESS) else NA_real_,
                                max_time_to_target_min_ESS  = if (complete_target_timings) max(readable_runs$time_to_target_min_ESS) else NA_real_,
                                median_max_Rhat             = fn_median_or_NA(readable_runs$max_Rhat),
                                median_Rhat_floor           = fn_median_or_NA(readable_runs$Rhat_floor_from_chain_length),
                                median_max_nRhat            = fn_median_or_NA(readable_runs$max_nRhat),
                                median_max_Rhat_main        = fn_median_or_NA(readable_runs$max_Rhat_main),
                                median_max_nRhat_main       = fn_median_or_NA(readable_runs$max_nRhat_main),
                                mean_pct_divs               = fn_mean_or_NA(readable_runs$pct_divs),
                                mean_eps                    = fn_mean_or_NA(readable_runs$eps_main),
                                mean_L                      = fn_mean_or_NA(readable_runs$L_main_samp),
                                mean_grad_evals_per_sec     = fn_mean_or_NA(readable_runs$grad_evals_per_sec),
                                mean_sampling_acceptance_probability = fn_mean_or_NA(readable_runs$sampling_acceptance_probability),
                                group_key                   = group_key_of_each_run[match(runs_of_group$file[1], runs_table$file)])
              })
              configuration_table <-  dplyr::bind_rows(configuration_rows)
        }
        ##
        ## ---- sort: N first (if asked and N varies), then each sort_by metric in turn, NA always last:
        ##
        {
              sort_keys <-  list()
              if (sort_within_N && "N" %in% names(configuration_table)) sort_keys[["N"]] <-  as.numeric(as.character(configuration_table$N))
              for (sort_index in seq_len(nrow(sort_by_specification))) {
                    sort_column_values <-  configuration_table[[sort_by_specification$column_name[sort_index]]]
                    sort_keys[[sort_by_specification$column_name[sort_index]]] <-  if (sort_decreasing[sort_index]) -sort_column_values else sort_column_values
              }
              configuration_table <-  dplyr::arrange(.data = configuration_table, !!!unname(sort_keys))
              ##
              ## rank within N (1 = best on the FIRST sort_by metric; ties share a rank):
              ##
              first_sort_values <-  configuration_table[[sort_by_specification$column_name[1]]]
              if (sort_decreasing[1]) first_sort_values <-  -first_sort_values
              N_of_each_row <-  if ("N" %in% names(configuration_table)) as.character(configuration_table$N) else rep("all", nrow(configuration_table))
              configuration_table$rank <-  ave(first_sort_values, N_of_each_row, FUN = function(values) rank(values, na.last = "keep", ties.method = "min"))
        }
        ##
        ## ---- print:
        ##
        if (print_table) {
              constant_settings_text <-  paste(paste0(constant_setting_names, " = ", sapply(constant_setting_names, function(setting_name) as.character(runs_table[[setting_name]][1]))),
                                              collapse = " | ")
              target_min_ESS_by_N <-  dplyr::distinct(.data = runs_table, model_type, N, target_min_ESS)
              ##
              cat("\n================ PS7 RUNS: ONE ROW PER ", if (is.null(group_by)) "CONFIGURATION" else "GROUP", " ================\n", sep = "")
              cat(paste0(formatC(as.integer(nrow(runs_table)), format = "d"), " runs in ", formatC(as.integer(nrow(configuration_table)), format = "d"), " rows. Every number is over ALL readable runs (means; medians for R-hats).\n"))
              cat("sorted by: ", paste0(sort_by_specification$short_name, ifelse(sort_decreasing, " (desc)", " (asc)"), collapse = ", "),
                  if (sort_within_N) " - within each N" else "", "\n", sep = "")
              cat("target min ESS: ", paste0(target_min_ESS_by_N$model_type, " N = ", target_min_ESS_by_N$N, " -> ", target_min_ESS_by_N$target_min_ESS, collapse = " | "), "\n", sep = "")
              cat("\nsettings shared by ALL runs:\n", strwrap(constant_settings_text, width = 120, prefix = "  "), sep = "\n")
              cat("\nsettings labelling the rows: ", paste(label_setting_names, collapse = ", "), "\n", sep = "")
              if (length(pooled_setting_names) > 0) cat("settings POOLED over (they vary but are not in group_by): ", paste(pooled_setting_names, collapse = ", "), "\n", sep = "")
              cat("\n")
              ##
              ## ---- compact printed view (the returned table keeps the full names and every column):
              ##
              short_setting_labels <-  c( n_chains_burnin                      = "chains_b",
                                         n_chains_sampling                    = "chains_s",
                                         n_threads_WCP_burnin                 = "WCP_b",
                                         n_threads_WCP_sampling               = "WCP_s",
                                         num_chunks_burnin                    = "chunks_b",
                                         num_chunks_sampling                  = "chunks_s",
                                         diffusion_HMC_integrator             = "integ",
                                         learning_rate                        = "LR",
                                         learning_rate_initial                = "LR_init",
                                         burnin_algorithm                     = "algorithm",
                                         tau_objective                        = "legacy_obj",
                                         tau_weight_by_p_jump                 = "pjump",
                                         metric_estimator                     = "metric",
                                         metric_type_nuisance                 = "nu_mass",
                                         eps_reinit_at_ChEES_handover         = "eps_reinit",
                                         theta_hat_us_rule                    = "centre",
                                         theta_hat_us_freeze_iter             = "centre_freeze",
                                         burnin_schedule_version              = "schedule_v",
                                         n_adapt                              = "fixed_from",
                                         metric_adaptation_end_iter            = "metric_end",
                                         centre_adaptation_end_iter            = "centre_end",
                                         pre_burnin_n_iter                    = "pre_n",
                                         pre_burnin_L                         = "pre_L",
                                         share_tau_ii_across_chains_in_burnin = "shared_tau",
                                         burnin_TBB_pool_equals_n_chains      = "TBB_chains",
                                         J_grad_option                        = "J_grad",
                                         autodiff_fallback                    = "AD_fallback",
                                         store_log_lik_trace                  = "store_log_lik",
                                         randomize_tau_burnin                 = "jitter_burnin",
                                         tau_sampling_scale                   = "tau_samp_scale",
                                         tau_adam_bias_correction             = "adam_bias_corr",
                                         tau_adaptation_block                 = "tau_block",
                                         eps_acceptance_mean                  = "eps_accept_mean",
                                         tau_shrink_on_divergence             = "tau_div_shrink",
                                         metric_pooled_window_resets          = "pooled_resets",
                                         metric_pooled_offdiagonal_shrinkage  = "pooled_offdiag_shrink",
                                         test_perm_override                   = "test_order")
              ##
              rows_to_print <-  if (is.finite(top_n)) configuration_table$rank <= top_n & !is.na(configuration_table$rank) else rep(TRUE, nrow(configuration_table))
              table_to_print <-  dplyr::filter(.data = configuration_table, rows_to_print)
              ##
              printed_view <-  dplyr::select(.data = table_to_print, rank, dplyr::all_of(label_setting_names))
              if ("diffusion_HMC_integrator" %in% names(printed_view)) {
                    printed_view$diffusion_HMC_integrator <-  R_fn_ps7_integrator_label(printed_view$diffusion_HMC_integrator)
              }
              names(printed_view) <-  ifelse(names(printed_view) %in% names(short_setting_labels), short_setting_labels[names(printed_view)], names(printed_view))
              ##
              fn_format_or_dash <-  function(values, number_format) {
                    number_format_digits <-  switch(number_format,
                                                     "%.0f" = 0L,
                                                     "%.1f" = 1L,
                                                     "%.2f" = 2L,
                                                     "%.3f" = 3L,
                                                     stop("Unsupported metric number format: ", number_format))
                    ifelse(is.finite(values), sub("^ +", "", formatC(values, format = "f", digits = number_format_digits)), "-")
              }
              for (print_index in seq_len(nrow(columns_to_print_specification))) {
                    metric_short_name  <-  columns_to_print_specification$short_name[print_index]
                    metric_column_name <-  columns_to_print_specification$column_name[print_index]
                    printed_view[[metric_short_name]] <-  switch(metric_short_name,
                          nRhat_ok    = paste0(table_to_print$n_nRhat_ok, "/", table_to_print$n_runs),
                          Rhat_ok     = paste0(table_to_print$n_Rhat_ok,  "/", table_to_print$n_runs),
                          ESS_grad    = ifelse(is.finite(table_to_print$mean_ESS_per_grad_x1000),
                                               paste0(sub("^ +", "", formatC(table_to_print$mean_ESS_per_grad_x1000, format = "f", digits = 1)), " [", sub("^ +", "", formatC(table_to_print$min_ESS_per_grad_x1000, format = "f", digits = 1)), "-", sub("^ +", "", formatC(table_to_print$max_ESS_per_grad_x1000, format = "f", digits = 1)), "]"), "-"),
                          time_to_target_ESS = ifelse(is.finite(table_to_print$mean_time_to_target_min_ESS),
                                               paste0(sub("^ +", "", formatC(table_to_print$mean_time_to_target_min_ESS, format = "f", digits = 1)), " [", sub("^ +", "", formatC(table_to_print$min_time_to_target_min_ESS, format = "f", digits = 1)), "-", sub("^ +", "", formatC(table_to_print$max_time_to_target_min_ESS, format = "f", digits = 1)), "]"), "-"),
                          fn_format_or_dash(table_to_print[[metric_column_name]], columns_to_print_specification$number_format[print_index]))
              }
              print(x = printed_view, n = Inf, width = Inf)
              if (sum(!rows_to_print) > 0) cat(paste0("\n  (", formatC(as.integer(sum(!rows_to_print)), format = "d"), " more rows not printed: top_n = ", as.character(top_n), ")\n"))
              ##
              cat("\n")
              for (print_index in seq_len(nrow(columns_to_print_specification))) {
                    cat(paste0("  ", formatC(as.character(columns_to_print_specification$short_name[print_index]), format = "s", width = 12, flag = "-"), " ", as.character(columns_to_print_specification$description[print_index]), "\n"))
              }
              cat("  (sort / print any other metric by name - fn_ps7_print_metric_names() lists them)\n")
        }
        ##
        return(invisible(list( configurations     = configuration_table,
                               runs               = runs_table,
                               label_settings     = label_setting_names,
                               pooled_settings    = pooled_setting_names,
                               constant_settings  = runs_table[1, constant_setting_names, drop = FALSE])))

}


##
## ---- read a directory of saved ps7 runs and summarise it (the replacement for summarize_ps7_results):
##
fn_ps7_summarise_directory <-  function( dir,
                                        pattern           = "^ps7_run_",
                                        target_min_ESS    = NULL,
                                        sort_by           = "ESS_grad",
                                        sort_decreasing   = NULL,
                                        sort_within_N     = TRUE,
                                        group_by          = NULL,
                                        columns_to_print  = c("runs", "nRhat_ok", "ESS_grad", "ESS_sec", "min_ESS", "min_ESS_main",
                                                              "burnin_sec", "time_to_target_ESS", "Rhat", "Rhat_floor", "nRhat",
                                                              "divs_pct", "eps", "L", "accept_prob"),
                                        top_n             = Inf,
                                        max_Rhat_ok       = 1.05,
                                        max_nRhat_ok      = 1.05,
                                        max_pct_divs_ok   = 1.0,
                                        cache             = NULL,
                                        verbose           = TRUE,
                                        reconstruct_ta2_superchains = FALSE) {

        runs_table <-  fn_ps7_extract_tau_sweep( dir                    = dir,
                                                pattern                = pattern,
                                                cache                  = cache,
                                                verbose                = verbose,
                                                include_burnin_history = FALSE,
                                                reconstruct_ta2_superchains = reconstruct_ta2_superchains)
        ##
        return(fn_ps7_summarise_runs( runs_table        = runs_table,
                                      target_min_ESS    = target_min_ESS,
                                      sort_by           = sort_by,
                                      sort_decreasing   = sort_decreasing,
                                      sort_within_N     = sort_within_N,
                                      group_by          = group_by,
                                      columns_to_print  = columns_to_print,
                                      top_n             = top_n,
                                      max_Rhat_ok       = max_Rhat_ok,
                                      max_nRhat_ok      = max_nRhat_ok,
                                      max_pct_divs_ok   = max_pct_divs_ok))

}






















