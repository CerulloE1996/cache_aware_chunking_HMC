##
## -| --------- Read the tau-scheme sweep -----------------------------------------------------------------------
##
## Consumes the frame from fn_ps7_extract_tau_sweep() (one row per saved RUN) and answers, in order:
##
##   0. STRATA      - what is actually in the directory. Runs from an earlier sweep at a different
##                    tau_initial / metric_estimator sit in the same folder; pooling them would
##                    compare settings that differ in more than the tau scheme, so they are reported
##                    and then excluded rather than silently averaged in.
##   1. L CURVE     - ESS per gradient against the PINNED L, per integrator, with the adapted schemes
##                    marked on the same scale. The question the adaptation cannot answer itself.
##   2. ADAPTIVE    - rate/log-rate/SNAPER comparisons matched on other settings, plus the run-to-run SPREAD of
##                    the adapted L. Spread alone does not establish a random walk; use the saved
##                    within-run tau and epsilon histories to investigate that mechanism.
##   3. BY LR       - the L curve split by learning rate, because an optimum that moves with LR would
##                    be hidden by the average over LRs.
##
## Runs that did not converge are excluded everywhere and listed with a reason - a stuck run reports
## min_ESS = Inf and would otherwise rank first on every speed metric.
##
#' @param d            data frame from fn_ps7_extract_tau_sweep().
#' @param tau_initial_keep  which tau_initial stratum to analyse; NULL = the most common one.
#' @param max_nRhat_ok      nested-Rhat limit.
#' @param max_pct_divs_ok   per cent divergent transitions limit.
#' @return invisibly, list(curve, tau_curve, adaptive, paired, by_lr, dropped).
#' @export
fn_ps7_report_tau_sweep <-  function( d,
                                     tau_initial_keep = NULL,
                                     max_nRhat_ok = 1.05,
                                     max_pct_divs_ok = 1.0) {

        d <-  tibble::as_tibble(x = d)
        if (!"manual_tau_value" %in% names(x = d)) d$manual_tau_value <-  NA_real_
        fmt <-  function(x, k = 3) sub("^ +", "", formatC(x = x, format = "f", digits = k))
        ##
        ## ---------------------------------------------------------------- 0. strata
        cat("\n================ 0. WHAT IS IN THE DIRECTORY ================\n")
        strata <-  d %>%
            dplyr::filter(dplyr::if_all(c(tau_initial, metric_estimator, n_burnin), ~ !is.na(.x))) %>%
            dplyr::count(tau_initial, metric_estimator, n_burnin, name = "Freq") %>%
            dplyr::arrange(.data$n_burnin, .data$metric_estimator, .data$tau_initial)
        print(x = strata, n = Inf, width = Inf)
        ##
        if (is.null(tau_initial_keep)) {
            tau_initial_keep <-  names(x = sort(x = table(d$tau_initial), decreasing = TRUE))[1]
        }
        d <-  dplyr::filter(.data = d, .data$tau_initial == .env$tau_initial_keep)
        ##
        ## ---- Pin the OTHER outer axes too.
        ##
        ## The KE scheme encodes to nothing in the file name, so its cell is also filled by every
        ## earlier sweep that happened to share this tau_initial - runs at a different
        ## metric_estimator or learning_rate_initial. Leaving them in would compare ChEES against a
        ## KE group drawn from a wider set of settings, and the run-to-run spread of the adapted L -
        ## the headline number here - would be inflated for KE by the extra settings rather than by
        ## the objective. Pin both axes to the level the ChEES runs used, so every scheme is compared
        ## on identical ground.
        ##
        ## Runs saved with "_ta4" encode EVERY algorithm explicitly ("_bake" for KE), so no KE cell can be filled by an
        ## earlier sweep. Pinning is then not needed, and would silently drop every other swept learning_rate_initial /
        ## metric_estimator level for ALL algorithms as soon as ChEES is among them.
        ##
        legacy_unencoded_KE_runs_present <-  !"tau_adaptation_version" %in% names(x = d) || any(d$tau_adaptation_version < 4)
        for (fld in c("metric_estimator", "learning_rate_initial")) {
            if (!legacy_unencoded_KE_runs_present) next
            if (!fld %in% names(x = d)) next
            chees_rows <-  d %>%
                dplyr::filter(.data$tau_objective == "ChEES", is.na(.data$manual_L), is.na(.data$manual_tau_value)) %>%
                dplyr::pull(dplyr::all_of(fld))
            if (length(x = chees_rows) == 0) next
            level <-  names(x = sort(x = table(chees_rows), decreasing = TRUE))[1]
            d <-  dplyr::filter(.data = d, as.character(.data[[fld]]) == .env$level)
            cat("  pinned ", fld, " = ", level, "\n", sep = "")
        }
        cat("\nanalysing tau_initial = '", tau_initial_keep, "' only ",
            "(mixing strata would compare settings that differ in more than the tau scheme)\n", sep = "")
        ##
        ## ---------------------------------------------------------------- convergence filter
        bad_read <-  !d$readable
        bad_ess  <-  !is.finite(d$ESS_per_grad_samp) | !is.finite(d$min_ESS)
        bad_rhat <-  is.finite(d$max_nRhat) & d$max_nRhat > max_nRhat_ok
        bad_divs <-  is.finite(d$pct_divs) & d$pct_divs > max_pct_divs_ok
        keep <-  !(bad_read | bad_ess | bad_rhat | bad_divs)
        ##
        if (any(!keep)) {
            reason <-  ifelse(bad_read, "unreadable",
                      ifelse(bad_ess, "non-finite ESS",
                      ifelse(bad_rhat, "nested Rhat too high", "too many divergences")))
            cat("\nEXCLUDED ", sum(!keep), " of ", length(x = keep), " runs:\n", sep = "")
            print(table(scheme = d$tau_scheme[!keep], reason = reason[!keep]))
        } else {
            cat("\nall ", nrow(x = d), " runs in this stratum converged\n", sep = "")
        }
        dropped <-  dplyr::filter(.data = d, !keep)
        d <-  dplyr::filter(.data = d, keep)
        if (nrow(x = d) == 0) { cat("nothing left to report.\n"); return(invisible(NULL)) }
        ##
        agg <-  function(x, by, FUN = mean) tapply(X = x, INDEX = by, FUN = FUN, simplify = TRUE)
        ##
        ## ---------------------------------------------------------------- 1. the L curve
        pinned <-  dplyr::filter(.data = d, !is.na(.data$manual_L))
        adaptive <-  dplyr::filter(.data = d, is.na(.data$manual_L), is.na(.data$manual_tau_value))
        curve <-  NULL
        tau_curve <-  NULL
        pinned_tau <-  dplyr::filter(.data = d, !is.na(.data$manual_tau_value))
        if (nrow(x = pinned_tau) > 0) {
            cat("\n================ FIXED-TAU CURVE, BY INTEGRATOR AND LEARNING RATE ================\n")
            tau_curve <-  pinned_tau %>%
                dplyr::filter(dplyr::if_all(c(manual_tau_value, diffusion_HMC_integrator, learning_rate), ~ !is.na(.x))) %>%
                dplyr::group_by(.data$manual_tau_value, .data$diffusion_HMC_integrator, .data$learning_rate) %>%
                dplyr::summarise(dplyr::across(c(ESS_per_grad_samp, ESS_per_sec_samp, L_main_samp, eps_main), mean), .groups = "drop") %>%
                dplyr::arrange(.data$learning_rate, .data$diffusion_HMC_integrator, .data$manual_tau_value)
            print(x = tau_curve, n = Inf, width = Inf)
        }
        ##
        if (nrow(x = pinned) > 0) {
            cat("\n================ 1. FIXED-L CURVE: ESS per gradient ================\n")
            cat("  (mean over learning rates and runs; n = runs behind each cell)\n\n")
            rows <-  list()
            for (L in sort(x = unique(x = pinned$manual_L))) {
                e <-  list(L = L)
                for (intg in sort(x = unique(x = pinned$diffusion_HMC_integrator))) {
                    s <-  pinned[pinned$manual_L == L & pinned$diffusion_HMC_integrator == intg, ]
                    tag <-  if (intg == "kick_flow_kick") "KFK" else "FKF"
                    e[[paste0(tag, "_ESSgrad")]] <-  mean(x = s$ESS_per_grad_samp) * 1000
                    e[[paste0(tag, "_se")]] <-  stats::sd(x = s$ESS_per_grad_samp) / sqrt(nrow(s)) * 1000
                    e[[paste0(tag, "_n")]] <-  nrow(x = s)
                }
                rows[[length(rows) + 1]] <-  tibble::as_tibble(x = e)
            }
            curve <-  dplyr::bind_rows(rows)
            print(x = curve, n = Inf, width = Inf)
            cat("\n  (ESS per 1000 gradients; se = standard error over the runs in that cell)\n")
            ##
            for (tag in c("KFK", "FKF")) {
                col <-  paste0(tag, "_ESSgrad")
                if (!col %in% names(x = curve)) next
                ## an L at which EVERY run of this integrator failed to converge has no mean at all;
                ## exclude those cells from the peak rather than letting NaN propagate
                ok <-  is.finite(curve[[col]])
                if (!any(ok)) { cat("\n  ", tag, ": no converged runs at any L\n", sep = ""); next }
                best <-  curve$L[ok][which.max(x = curve[[col]][ok])]
                failed_at <-  curve$L[!ok]
                cat("\n  ", tag, ": best L = ", best,
                    "   peak ", fmt(max(curve[[col]][ok]), 2),
                    " vs ", fmt(min(curve[[col]][ok]), 2), " at the worst converged L",
                    if (length(x = failed_at) > 0)
                        paste0("\n       NO CONVERGED RUNS AT L = ", paste(failed_at, collapse = ", ")) else "",
                    "\n", sep = "")
            }
        }
        ##
        ## ---------------------------------------------------------------- 2. adaptive schemes
        paired <-  NULL
        if (nrow(x = adaptive) > 0) {
            cat("\n================ 2. ADAPTIVE SCHEMES ================\n\n")
            ##
            ## ---- split by scheme x integrator AND every other setting that varies among these runs (burn-in length,
            ##      chains, WCP, learning rates, ...), so a row never pools different configurations:
            ##
            other_setting_names <-  setdiff(names(fn_ps7_parse_run_name(x = adaptive$file[1])),
                                           c("run", "configuration", "file", "tau_objective", "burnin_algorithm", "tau_weight_by_p_jump",
                                             "manual_L", "manual_tau_value", "diffusion_HMC_integrator"))
            other_varying_setting_names <-  other_setting_names[sapply(other_setting_names, function(setting_name) length(unique(as.character(adaptive[[setting_name]]))) > 1)]
            fn_group_key <-  function(runs) do.call(paste, c(list(runs$tau_scheme, runs$diffusion_HMC_integrator),
                                                              lapply(other_varying_setting_names, function(setting_name) as.character(runs[[setting_name]])),
                                                              sep = "|"))
            adaptive_group_key <-  fn_group_key(adaptive)
            rows <-  list()
            for (group_key_value in sort(x = unique(x = adaptive_group_key))) {
                    s <-  adaptive[adaptive_group_key == group_key_value, , drop = FALSE]
                    rows[[length(rows) + 1]] <-  tibble::tibble(
                        scheme = s$tau_scheme[1],
                        integrator = R_fn_ps7_integrator_label(s$diffusion_HMC_integrator[1]),
                        !!!dplyr::select(.data = dplyr::slice(.data = s, 1), dplyr::all_of(other_varying_setting_names)),
                        n = nrow(x = s),
                        ESSgrad = mean(x = s$ESS_per_grad_samp) * 1000,
                        L_mean = mean(x = s$L_main_samp),
                        L_sd = stats::sd(x = s$L_main_samp),
                        L_cv_pct = 100 * stats::sd(x = s$L_main_samp) / mean(x = s$L_main_samp),
                        pct_divs = mean(x = s$pct_divs, na.rm = TRUE))
            }
            adaptive_tab <-  dplyr::bind_rows(rows)
            print(x = adaptive_tab, n = Inf, width = Inf)
            cat("\n  L_cv_pct is the run-to-run spread of the ADAPTED trajectory length.\n",
                "  This spread also mixes learning-rate effects and is not a test of a random walk.\n",
                "  Inspect the saved tau/epsilon histories within matched settings before making\n",
                "  claims about adaptation convergence.\n", sep = "")
            ##
            ## Prefer the original rate criterion when comparing the newly requested variants.
            baseline_candidates <-  c("CHESSR_pjump", "CHESSR", "KE_pjump", "KE")
            baseline_scheme <-  baseline_candidates[baseline_candidates %in% adaptive$tau_scheme][1L]
            if (!is.na(baseline_scheme)) {
                message(paste0("Adaptive schemes relative to ", baseline_scheme, ", matched on all other varying settings."))
                ## (was integrator x learning_rate only, which paired runs that differed in other settings)
                key <-  function(s) do.call(paste, c(list(s$diffusion_HMC_integrator, s$learning_rate),
                                                    lapply(other_varying_setting_names, function(setting_name) as.character(s[[setting_name]])),
                                                    sep = "|"))
                rows <-  list()
                for (alt in setdiff(x = sort(x = unique(x = adaptive$tau_scheme)), y = baseline_scheme)) {
                    b <-  adaptive[adaptive$tau_scheme == baseline_scheme, ]
                    a <-  adaptive[adaptive$tau_scheme == alt, ]
                    bm <-  tapply(b$ESS_per_grad_samp, key(b), mean)
                    am <-  tapply(a$ESS_per_grad_samp, key(a), mean)
                    shared <-  intersect(x = names(bm), y = names(am))
                    if (length(x = shared) == 0) next
                    rows[[length(rows) + 1]] <-  tibble::tibble(
                        scheme = alt, baseline = baseline_scheme, n_pairs = length(x = shared),
                        pct_change = 100 * mean(x = (am[shared] - bm[shared]) / bm[shared]),
                        wins = paste0(sum(am[shared] > bm[shared]), "/", length(x = shared)))
                }
                if (length(x = rows) > 0) {
                    paired <-  dplyr::bind_rows(rows)
                    print(x = paired, n = Inf, width = Inf)
                    message(paste0("pct_change = mean per cent change in ESS per gradient relative to ", baseline_scheme, "."))
                }
            }
        }
        ##
        ## ---------------------------------------------------------------- 3. does the optimum move with LR?
        by_lr <-  NULL
        if (nrow(x = pinned) > 0 && length(x = unique(x = pinned$learning_rate)) > 1) {
            cat("\n================ 3. BEST L, SPLIT BY LEARNING RATE ================\n")
            cat("  (if these disagree, the averaged curve above is hiding an interaction)\n\n")
            rows <-  list()
            for (lr in sort(x = unique(x = pinned$learning_rate))) {
                for (intg in sort(x = unique(x = pinned$diffusion_HMC_integrator))) {
                    s <-  pinned[pinned$learning_rate == lr & pinned$diffusion_HMC_integrator == intg, ]
                    if (nrow(x = s) == 0) next
                    m <-  tapply(s$ESS_per_grad_samp, s$manual_L, mean)
                    rows[[length(rows) + 1]] <-  tibble::tibble(
                        learning_rate = lr,
                        integrator = R_fn_ps7_integrator_label(intg),
                        best_L = as.numeric(names(m)[which.max(x = m)]),
                        peak_ESSgrad = max(m) * 1000)
                }
            }
            by_lr <-  dplyr::bind_rows(rows)
            print(x = by_lr, n = Inf, width = Inf)
        }
        ##
        invisible(list(curve = curve, tau_curve = tau_curve, adaptive = if (exists("adaptive_tab")) adaptive_tab else NULL,
                       paired = paired, by_lr = by_lr, dropped = dropped))

}






















