
## ps_7_MCMC_settings_BayesMVP_functions.R


## ==========================================================================================================================================
## Pilot Study 7: BayesMVP settings
## ==========================================================================================================================================
##
##   1. R_fn_ps7_num_chunks()    -- NEW helper. Single source of truth for the
##      deterministic N -> num_chunks mapping (extracted from the old inline
##      if/else chain in run_ps7_models(); values IDENTICAL). Needed because the
##      resume check must know num_chunks BEFORE the run loop (it's in the saved
##      filename via _n_chnk_b/_n_chnk_s). NB: the old code called
##      find_num_chunks_MVP() first and then ALWAYS overrode it via the if/else
##      chain (which has a catch-all else), so dropping that call is
##      behaviour-identical.
##
##   2. R_fn_file_name_string()  -- output IDENTICAL to before for all existing
##      calls; only change is `if (grouping == TRUE)` -> `if (isTRUE(grouping))`
##      so the NULL default no longer errors ("argument is of length zero").
##
##   3. run_ps7_models()         -- FIXED RESUME. The old existing-run check built
##      its filename with the RETIRED format (no _n_pops / _n_WCP_* / _n_chnk_*,
##      included true_p, used _M_dcy_pow/_M_dcy_scl instead of _M_pow/_M_scl), so
##      it NEVER matched what R_fn_file_name_string() saves -> every run was
##      always re-done. The check now builds the name via R_fn_file_name_string()
##      itself, using values (n_pops, num_chunks, WCP thread counts, overridden
##      priors) that are hoisted ABOVE the check and then REUSED inside the run
##      loop -- so check-name == save-name by construction, forever.
##      Two further protective changes: (a) the post-sampling metric-diagnostics
##      block (generalised-eigenvalue check on M_inv_dense_main) is wrapped in
##      try() -- with metric_shape_main = "diag" (the ordinal config)
##      M_inv_dense_main may not exist, and a hard error there would throw away a
##      COMPLETED run before it gets saved; (b) large dead commented-out
##      exploration blocks are stripped (live code is verbatim; toggle-comments
##      for priors/settings are kept -- the original file remains the archive for
##      the stripped notes).
##
##   4. summarize_ps7_results()  -- FIXED LOADER. Its inline filename was missing
##      "_n_pops" (which R_fn_file_name_string() writes between _N and _diff_HMC),
##      so it could never find the saved run files (-> silently empty tibble).
##      It now builds the name via R_fn_file_name_string() too (one definition,
##      three users). New argument `n_pops` (default: 3 for LC_MVOP, else 2 --
##      i.e. the ordinal 3-pop and binary COVID 2-pop DGMs). NB: the saved format
##      writes model_args_list$num_chunks into BOTH _n_chnk_b and _n_chnk_s, so
##      num_chunks_sampling is retained in the signature but ignored --
##      num_chunks_burnin is what must match the files on disk.
##
## Unchanged (NOT in this file): get_target_ESS(), summarize_combos(),
## make_combos(), summarize_ps7_by_config(), generate_bayesmvp_latex_table(),
## generate_all_bayesmvp_tables().
##
##
## ---- 1. NEW helper: deterministic N -> num_chunks mapping  --------------------------------------------------------------------------------
##
R_fn_ps7_num_chunks <-  function( N, 
                                  n_tests
) {
  
    n_obs <-  N*n_tests
  
    if (parallel::detectCores() > 16) {  ## local HPC
      
          if (n_obs == 60000) { ## N=10,000 for binary
            num_chunks <-  25
          } else if (n_obs == 15000) { ## N=2500 for binary, N=5000 for ordinal
            num_chunks <-  10
          } else if (n_obs == 3000) {  ## N=500 for binary, N=1000 for ordinal
            num_chunks <-  1    ## 2
          } else if (n_obs == 150000) { ## N=25,000 for binary
            num_chunks <-  50
          } else if (n_obs == 300000) { ## N=50,000 for binary
            num_chunks <-  125 ## 100
          # } else if (n_obs == 15000) {  ## ----------------------------------- ordinal-only
          #   num_chunks <-  10
          # } else if (n_obs == 1000) {  ## ----------------------------------- ordinal-only
          #   num_chunks <-  1    ## 2
          } else if (n_obs == 750) {   ## ----------------------------------- ordinal-only
            num_chunks <-  1    ## 2#
          } else {
            num_chunks <-  1    ## 2
          }
      
    } else { ## Laptop
      
          if (n_obs == 60000) {
            num_chunks <-  50
          } else if (n_obs == 15000) {
            num_chunks <-  10
          } else if (n_obs == 3000) {
            num_chunks <-  1
          } else if (n_obs == 150000) {
            num_chunks <-  125
          } else if (n_obs == 300000) {
            num_chunks <-  250
          # } else if (n_obs == 5000) {  ## ----------------------------------- ordinal-only
          #   num_chunks <-  10
          # } else if (n_obs == 1000) {  ## ----------------------------------- ordinal-only
          #   num_chunks <-  1    ## 2
          } else if (n_obs == 750) {   ## ----------------------------------- ordinal-only
            num_chunks <-  1    ## 2#
          } else {
            num_chunks <-  1    ## 2
          }
      
    }

        return(num_chunks)

}
##
## ---- 2. Filename builder (ONE definition -- used by run_ps7_models() check + save AND summarize_ps7_results()): --------------------------
##


##
## ---- NAME TAGS SHORTENED (only the tag TEXT; every value is unchanged, so no setting is lost):
##   _n_pops->_np  _diff_HMC->_dH  _pt_HMC->_pH  _n_ch_b->_cb  _n_WCP_b->_wb  _n_chnk_b->_kb
##   _n_b->_b  _n_i->_it  _ratio_M->_rM  _M_typ_main->_mt  _nuis->_nu  _M_shp->_ms  _width->_w
##   _M_dcy->_md  pi_LKJ->_LKJ  pi_p->_pp  _m_att->_ma  _M_pow->_mp  _M_scl->_msc
##   _ORIG_COLS->_OC  _M_pooled->_Mp  _ord_grp1->_og1
##
## This frees ~60 bytes of the 247-byte NAME_MAX budget so the integrator AND both warm-start
## fields fit. FILES SAVED UNDER THE OLD TAGS WILL NO LONGER BE MATCHED by the resume check or
## by summarize_ps7_results() - rename or re-run them.
##
## ---- Value encoders (injective -- nothing is lost, only spelled shorter):
##
R_fn_enc_lgl <-  function(x) if (isTRUE(as.logical(x))) "1" else "0"
##
R_fn_enc_num <-  function(x) sub("^(-?)0\\.", "\\1.", as.character(x))   ## 0.05 -> .05
##
## NULL / NA-safe version. paste0("_x", NULL) silently returns "_x", which would fuse two
## fields together and make two different runs share a file name; "d" = the sampler's default,
## "off" = the feature switched off:
R_fn_enc_opt <-  function(x, off_is_null = TRUE) {
        if (is.null(x)) return(if (off_is_null) "d" else "off")
        if (length(x) != 1L || is.na(x)) return("off")
        return(R_fn_enc_num(x))
}
##
## ---- tau_initial levels: a number (pi, 2*pi, ...) or "adaptive" (NicoStan sets tau at the handover to (pi/2) * sqrt(lambda_max),
##      lambda_max from the burn-in draws' covariance in metric coordinates; "_tiA" in the file name).
##      c(pi, 2*pi, "adaptive") is a CHARACTER vector, and as.character(pi) keeps 15 significant digits ("3.14159265358979",
##      which is not pi). A numeric level given as text is therefore read back as a number and, when it lies within 1e-12
##      (relative) of a multiple of pi/8, replaced by that exact multiple, so the sampler and the file-name builder receive
##      exactly the pi and 2*pi that the numeric vector c(pi, 2*pi) gives them. A numeric level is returned unchanged.
##
fn_ps7_tau_initial_value <-  function(tau_initial) {

        if (is.null(tau_initial)) return(NULL)
        if (length(tau_initial) != 1) {
            stop("fn_ps7_tau_initial_value: one tau_initial level at a time; got ", length(tau_initial), ": ",
                 paste(as.character(tau_initial), collapse = ", "), ".")
        }
        if (is.factor(tau_initial)) tau_initial <-  as.character(tau_initial)
        if (is.numeric(tau_initial)) return(tau_initial)
        if (is.na(tau_initial)) return(NA_real_)
        if (identical(tau_initial, "adaptive")) return("adaptive")
        ##
        tau_initial_number <-  suppressWarnings(as.numeric(tau_initial))
        if (!is.finite(tau_initial_number) || tau_initial_number <= 0) {
            stop("tau_initial must be a positive finite number or \"adaptive\"; got '", tau_initial, "'.")
        }
        pi_eighths <-  round(tau_initial_number / (pi / 8))
        if (pi_eighths >= 1 && abs(tau_initial_number - pi_eighths * (pi / 8)) <= 1e-12 * tau_initial_number) {
            return((pi_eighths / 8) * pi)
        }
        return(tau_initial_number)

}
##
## ---- Encoders for the fields added AFTER the original naming scheme.
##
## The basename budget is 247 bytes and the existing scheme already uses most of it, so these
## emit NOTHING when the feature sits at the behaviour every earlier ps7 run had (hold off /
## warm start off / the default integrator). Two consequences, both wanted: a run that predates
## the field keeps EXACTLY the name it has today, so old results still resolve; and the budget
## is only spent when a setting is actually in play. Still injective - absent means "the old
## behaviour", which is one specific configuration, not "unknown".
##
R_fn_enc_hold <-  function(tag, value, iter) {
        if (is.null(value) || length(value) != 1L || is.na(value)) return("")   ## feature off
        paste0(tag, R_fn_enc_num(value), "_", if (is.null(iter)) "d" else R_fn_enc_num(iter))
}
##
R_fn_enc_dHMC_int <-  function(x) {
        if (is.null(x) || identical(as.character(x), "kick_flow_kick")) return("")  ## the default
        paste0("_i", R_fn_enc_map(x, R_fn_map_dHMC_int, field = "diffusion_HMC_integrator"))
}
##
## ---- Trajectory-length adaptation fields (added after the runs on N = 2500 were launched).
##
## All three follow the same rule as the integrator above: emit NOTHING at the behaviour every
## earlier run had, so those files keep exactly the names they have on disk today and the resume
## check and summarizer still find them. Absent therefore means the ORIGINAL ChESSR adaptation:
## kinetic-energy objective, median across chains, accept indicator, tau adapted (not pinned).
##
R_fn_enc_tau_obj <-  function(x) {
        if (is.null(x) || identical(as.character(x), "KE")) return("")   ## the original objective
        paste0("_to", R_fn_enc_map(x, R_fn_map_tau_obj, field = "burnin_algorithm"))
}
##
## New development runs record the sole algorithm selector explicitly.
fn_ps7_normalise_burnin_algorithm <-  function(burnin_algorithm) {
        ## aliases <-  c(ke = "KE", chees = "ChEES", chessr = "CHESSR", cheesr = "CHESSR",
                      ## chessr_log = "CHESSR_log", cheesr_log = "CHESSR_log", chees_per_tau = "CHESSR_log", snaper = "SNAPER")
        ## aliases <-  c(ke = "KE", chees = "ChEES", chessr = "CHESSR", cheesr = "CHESSR",
                      ## chessr_log = "CHESSR_log", cheesr_log = "CHESSR_log", chees_per_tau = "CHESSR_log", snaper = "SNAPER",
                      ## chessr_time = "CHESSR_time", cheesr_time = "CHESSR_time", snaper_time = "SNAPER_time")
        aliases <-  c(ke = "KE", chees = "ChEES", chessr = "CHESSR", cheesr = "CHESSR",
                      chessr_log = "CHESSR_log", cheesr_log = "CHESSR_log", chees_per_tau = "CHESSR_log", snaper = "SNAPER",
                      chessr_time = "CHESSR_time", cheesr_time = "CHESSR_time", snaper_time = "SNAPER_time",
                      esjd = "ESJD", esjd_chessr = "ESJD_CHESSR", esjd_cheesr = "ESJD_CHESSR",
                      esjd_snaper = "ESJD_SNAPER",
                      lq_essr = "LQ_ESSR", lqessr = "LQ_ESSR")
        keys <-  tolower(trimws(as.character(burnin_algorithm)))
        if (anyNA(keys) || any(!keys %in% names(aliases))) {
                ## stop("burnin_algorithm must contain KE, ChEES, CHESSR, CHESSR_log or SNAPER.")
                ## stop("burnin_algorithm must contain KE, ChEES, CHESSR, CHESSR_log, SNAPER, CHESSR_time or SNAPER_time.")
                ## stop("burnin_algorithm must contain KE, ChEES, CHESSR, CHESSR_log, SNAPER, CHESSR_time, SNAPER_time, ESJD, ESJD_CHESSR or ESJD_SNAPER.")
                stop(paste0("burnin_algorithm must contain KE, ChEES, CHESSR, CHESSR_log, SNAPER, CHESSR_time, SNAPER_time, ESJD, ESJD_CHESSR, ",
                            "ESJD_SNAPER or LQ_ESSR."))
        }
        return(unname(aliases[keys]))
}
##
## "_bact" = CHESSR_time and "_bast" = SNAPER_time (the wall-clock time-to-target-ESS variants of CHESSR and SNAPER; same byte
## count as "_bacr" / "_basn"). With non-default time-criterion settings the final "t" becomes 4 hexadecimal digits,
## "_bac<4 hex>" / "_bas<4 hex>" (see fn_ps7_time_criterion_settings_hash below):
##
## "_baej" = ESJD (expected squared jumped distance per unit trajectory length) and "_baex" = ESJD_CHESSR (the geometric mean of the
## ESJD and CHESSR criteria); the same byte count again, since many names are at the 249-byte limit. Neither can collide with the
## hashed "_bac<4 hex>" / "_bas<4 hex>" tokens (third letter "e", not "c" or "s"):
## "_baen" = ESJD_SNAPER (the geometric mean of the ESJD and SNAPER criteria); this token is distinct from the existing tokens and
## from the hashed time-criterion tokens:
## "_balq" = LQ_ESSR (the soft minimum of the linear and quadratic lag-one ESS bounds per unit trajectory length); third letter "l",
## so it cannot collide with "_bacl" (CHESSR_log) or the hashed "_bac<4 hex>" / "_bas<4 hex>" tokens:
##
## fn_ps7_encode_burnin_algorithm <-  function(burnin_algorithm) {
fn_ps7_encode_burnin_algorithm <-  function(burnin_algorithm,
                                           time_criterion_settings_hash = "") {
        algorithm <-  fn_ps7_normalise_burnin_algorithm(burnin_algorithm = burnin_algorithm)
        if (length(algorithm) != 1L) stop("A filename requires one burnin_algorithm.")
        ## non-default time-criterion settings: "_bac<4 hex>" / "_bas<4 hex>" in place of "_bact" / "_bast":
        if (nzchar(time_criterion_settings_hash) && algorithm %in% c("CHESSR_time", "SNAPER_time")) {
                return(paste0("_ba", if (identical(algorithm, "CHESSR_time")) "c" else "s", time_criterion_settings_hash))
        }
        ## tokens <-  c(KE = "ke", ChEES = "ce", CHESSR = "cr", CHESSR_log = "cl", SNAPER = "sn")
        ## tokens <-  c(KE = "ke", ChEES = "ce", CHESSR = "cr", CHESSR_log = "cl", SNAPER = "sn",
                     ## CHESSR_time = "ct", SNAPER_time = "st")
        ## tokens <-  c(KE = "ke", ChEES = "ce", CHESSR = "cr", CHESSR_log = "cl", SNAPER = "sn",
                     ## CHESSR_time = "ct", SNAPER_time = "st",
                     ## ESJD = "ej", ESJD_CHESSR = "ex", ESJD_SNAPER = "en")
        tokens <-  c(KE = "ke", ChEES = "ce", CHESSR = "cr", CHESSR_log = "cl", SNAPER = "sn",
                     CHESSR_time = "ct", SNAPER_time = "st",
                     ESJD = "ej", ESJD_CHESSR = "ex", ESJD_SNAPER = "en", LQ_ESSR = "lq")
        return(paste0("_ba", unname(tokens[algorithm])))
}
##
##
## ---- Wall-clock time-to-target-ESS burn-in criteria ("CHESSR_time" / "SNAPER_time"; "_bact" / "_bast") ------------------------------
##
## Objective. The wall-clock time to reach ESS_target_for_time_criterion (the target ESS over all sampling chains) is
##
##   time_to_target_ESS(tau) = n_iter_burnin * (time_per_iter_overhead_burnin + n_leapfrog_steps_per_iter * time_per_leapfrog_step_burnin)
##                           + (ESS_target_for_time_criterion / ESS_per_iter_sampling(tau))
##                             * (time_per_iter_overhead_sampling + time_per_iter_summaries_sampling
##                                + n_leapfrog_steps_per_iter * time_per_leapfrog_step_sampling),
##
## with n_leapfrog_steps_per_iter = tau / eps. Every time_* quantity is a wall-clock time in seconds:
##   time_per_leapfrog_step_burnin     one leapfrog step (one gradient evaluation for all burn-in chains, run in parallel) in the
##                                     burn-in configuration (n_chains_burnin, n_threads_WCP_burnin, num_chunks_burnin); measured
##                                     online by NicoStan during burn-in, from the wall time of each burn-in iteration after the clip
##                                     phase (a running least-squares fit on the leapfrog steps; its status is saved with the run);
##   time_per_leapfrog_step_sampling   the same in the sampling configuration (n_chains_sampling, n_threads_WCP_sampling,
##                                     num_chunks_sampling); the slope of the sampling timing probe fit;
##   time_per_iter_overhead_sampling   the fixed non-leapfrog time of one sampling iteration; the intercept of the probe fit;
##   time_per_iter_summaries_sampling  the summaries routine per sampling iteration (one draw from every sampling chain);
##   time_per_iter_overhead_burnin     the fixed non-leapfrog time of one burn-in iteration. It does not depend on tau, so it
##                                     drops out of the derivative below and is not measured.
## n_iter_burnin is the number of burn-in iterations run at the adapted tau: iterations clip_iter_tau + 1 to n_burnin, i.e.
## settings$n_burnin (the "_b" token) - clip_iter_tau (clip_iter + int). The earlier burn-in iterations (one leapfrog step before
## clip_iter, then the ramp to tau_initial up to the handover) and the pre-burnin (pre_burnin_L steps) take the same time whatever
## tau is chosen, so they drop out of the derivative below; NicoStan saves the count with its definition in the run. ESS_per_iter_sampling(tau)
## is the ESS per sampling iteration over all sampling chains.
##
## Setting d time_to_target_ESS / d log(tau) = 0 gives the target
##
##   ESS_elasticity_wrt_log_tau = time_to_target_ESS_tau_penalty
##                              = (1 + burnin_to_sampling_leapfrog_time_ratio) * tau / (tau + tau_offset_from_sampling_overhead),
##
##   ESS_elasticity_wrt_log_tau             = d log(ESS_per_iter_sampling) / d log(tau), the quantity the ChEES / SNAPER criterion
##                                            estimates (Sountsov and Hoffman 2022, eq. 2 and 5; Hoffman, Radul and Sountsov 2021);
##   sampling_overhead_in_leapfrog_steps    = (time_per_iter_overhead_sampling + time_per_iter_summaries_sampling)
##                                            / time_per_leapfrog_step_sampling;
##   tau_offset_from_sampling_overhead      = eps * sampling_overhead_in_leapfrog_steps;
##   burnin_to_sampling_leapfrog_time_ratio = (n_iter_burnin * time_per_leapfrog_step_burnin)
##                                            / (n_iter_sampling_for_time_criterion * time_per_leapfrog_step_sampling).
##
## CHESSR / SNAPER ascend ESS_elasticity_wrt_log_tau - 1 with ADAM on log(tau); CHESSR_time / SNAPER_time ascend
## ESS_elasticity_wrt_log_tau - time_to_target_ESS_tau_penalty. With burnin_to_sampling_leapfrog_time_ratio = 0 and
## sampling_overhead_in_leapfrog_steps = 0 they are CHESSR / SNAPER exactly; sampling_overhead_in_leapfrog_steps -> infinity
## gives ChEES.
##
## n_iter_sampling_for_time_criterion (the sampling iterations the criterion assumes) comes from
## n_iter_sampling_for_time_criterion_source:
##   "planned_n_iter"                 the run's planned n_iter (settings$n_iter, the "_it" token; the default);
##   "ESS_target_for_time_criterion"  ESS_target_for_time_criterion / ESS_per_iter_sampling_expected.
##                                    ESS_target_for_time_criterion = NULL takes the ps3 target min ESS for this N and model
##                                    (R_fn_ps7_get_target_min_ESS); ESS_per_iter_sampling_expected is user-supplied, or taken
##                                    from time_criterion_previous_run_path (min_ESS / n_iter of the saved run(s)).
## Which min ESS: time_criterion_ess_parameter_set chooses the parameters whose minimum ESS the target and
## ESS_per_iter_sampling_expected refer to: "diagnostic" (the default; the model's diagnostic parameters, i.e. the generated
## quantities of its Stan skeleton file: for LC_MVP the accuracy parameters p, Se_baseline, Sp_baseline and Fp_baseline, the
## parameters of the ps3 target min ESS and of the saved min_ESS), "main" (the parameters block), or parameter names / base names
## (e.g. c("p", "Se_baseline", "Sp_baseline")). It is used only when ESS_per_iter_sampling_expected is read from
## time_criterion_previous_run_path (the "ESS_target_for_time_criterion" source without a user value), and is inert (absent from
## the file name, the sampler call and the saved settings) otherwise; NicoStan records the set it applied in every
## CHESSR_time / SNAPER_time run.
##
## time_per_leapfrog_step_sampling, time_per_iter_overhead_sampling and time_per_iter_summaries_sampling, in order of precedence:
##   1. user-supplied numbers (the runner settings of the same names);
##   2. time_criterion_previous_run_path, one or more saved PS7 runs, read by NicoStan's fn_time_criterion_quantities_from_saved_run()
##      (one run: time_per_iter_overhead_sampling = 0, so time_per_leapfrog_step_sampling is an upper bound; runs with at least
##      two leapfrog step counts: least squares across the runs);
##   3. the sampling timing probe (run_sampling_timing_probe = TRUE, the default): NicoStan runs sampling_timing_probe_n_iter_per_L
##      iterations at each of sampling_timing_probe_L_values (default 10 at each of 2 and 8 leapfrog steps; plus a 3-iteration call
##      at each number of leapfrog steps, which removes the per-call set-up) in the sampling configuration, after the pre-burn-in
##      and before the burn-in, fits time_per_iter_sampling = time_per_iter_overhead_sampling + n_leapfrog_steps_per_iter *
##      time_per_leapfrog_step_sampling, times the summaries routine
##      on the probe draws, and discards the probe draws. Its wall time (sampling_timing_probe_wall_time) is included in time_burnin.
## The probe is skipped when all three are known from 1 or 2. Without any of them (probe off and nothing supplied) NicoStan falls
## back to burnin_to_sampling_leapfrog_time_ratio = 0 and sampling_overhead_in_leapfrog_steps = 0, i.e. plain CHESSR / SNAPER,
## with a warning.
##
## File names. Every setting at its default writes NOTHING beyond "_bact" / "_bast", and every setting is inert (absent from the
## name, the sampler call and the saved settings) for the other algorithms and whenever tau is pinned (manual_L, manual_tau_value
## or partitioned_HMC). Any EFFECTIVE non-default value replaces the final "t" of the algorithm token by 4 hexadecimal digits,
## "_bac<4 hex>" (CHESSR_time) or "_bas<4 hex>" (SNAPER_time), e.g. "_bac3f9a": the first 4 digits of the md5 of the
## time-criterion settings text (fn_ps7_time_criterion_settings_text), which lists every effective non-default setting as
## "name=value" (run_sampling_timing_probe = FALSE, sampling_timing_probe_L_values / sampling_timing_probe_n_iter_per_L, the
## ESS-target source and its resolved
## target, a non-default time_criterion_ess_parameter_set, user-supplied numbers at 15 significant digits, and the sorted BASENAMES of time_criterion_previous_run_path).
## The hash sits inside the algorithm token because many names are already within a few bytes of the 249-byte limit: its 3 extra
## bytes fit every saved CHESSR / SNAPER configuration after the "_LiLR" shortening, whereas readable tokens appended at the end
## (8 bytes for an ESS target alone) would not.
## Collisions: the text is saved in the run (time_criterion_inputs$time_criterion_settings_text), and run_ps7_models() stops
## before resuming a saved run whose text differs from the current one, so two settings can never share a saved run.
## fn_ps7_parse_run_name() reads the algorithm and the hash back (time_criterion_settings_hash); fn_ps7_extract_tau_sweep() adds
## the settings themselves from the saved run.
##
ps7_time_criterion_algorithms <-  c("CHESSR_time", "SNAPER_time")
##
ps7_time_criterion_default_settings <-  list(run_sampling_timing_probe = TRUE,
                                             sampling_timing_probe_L_values = c(2, 8),
                                             sampling_timing_probe_n_iter_per_L = 10,
                                             n_iter_sampling_for_time_criterion_source = "planned_n_iter",
                                             time_criterion_ess_parameter_set = "diagnostic")
##
## Every settings field of the time criteria (fn_ps7_effective_time_criterion_settings() removes all of them where they are inert):
ps7_time_criterion_setting_names <-  c("run_sampling_timing_probe",
                                       "sampling_timing_probe_L_values",
                                       "sampling_timing_probe_n_iter_per_L",
                                       "n_iter_sampling_for_time_criterion_source",
                                       "ESS_target_for_time_criterion",
                                       "ESS_per_iter_sampling_expected",
                                       "time_criterion_previous_run_path",
                                       "time_per_leapfrog_step_sampling",
                                       "time_per_iter_overhead_sampling",
                                       "time_per_iter_summaries_sampling",
                                       "time_criterion_ess_parameter_set")
##
## The three sampling times that make the probe unnecessary once all of them are supplied:
ps7_sampling_timing_quantity_names <-  c("time_per_leapfrog_step_sampling",
                                         "time_per_iter_overhead_sampling",
                                         "time_per_iter_summaries_sampling")
##
fn_ps7_is_time_criterion_algorithm <-  function(burnin_algorithm) {
        !is.null(burnin_algorithm) && length(burnin_algorithm) == 1 && !is.na(burnin_algorithm) &&
            burnin_algorithm %in% ps7_time_criterion_algorithms
}
##
## ---- Validation: NULL (or NA) = not supplied / the default; anything else must be one valid value (sampling_timing_probe_L_values:
##      at least two distinct whole numbers of leapfrog steps >= 1):
##
fn_ps7_validate_time_criterion_settings <-  function(settings) {

        fn_supplied_number <-  function(supplied_value, setting_name, must_be_positive) {
            if (is.null(supplied_value)) return(NULL)
            if (length(supplied_value) != 1) stop(paste0("settings$", setting_name, " must be NULL or ONE number per run."))
            if (is.na(supplied_value)) return(NULL)
            if (!is.numeric(supplied_value) || !is.finite(supplied_value) || supplied_value < 0 || (must_be_positive && supplied_value == 0)) {
                stop(paste0("settings$", setting_name, " must be NULL or one finite number ", if (must_be_positive) "> 0" else ">= 0", "."))
            }
            return(as.numeric(supplied_value))
        }
        ##
        run_sampling_timing_probe_requested <-  BayesMVP::if_null_then_set_to(settings$run_sampling_timing_probe,
                                                                              ps7_time_criterion_default_settings$run_sampling_timing_probe)
        sampling_timing_probe_L_values_requested <-  BayesMVP::if_null_then_set_to(settings$sampling_timing_probe_L_values,
                                                                                   ps7_time_criterion_default_settings$sampling_timing_probe_L_values)
        sampling_timing_probe_n_iter_per_L_requested <-  BayesMVP::if_null_then_set_to(settings$sampling_timing_probe_n_iter_per_L,
                                                                                       ps7_time_criterion_default_settings$sampling_timing_probe_n_iter_per_L)
        n_iter_sampling_for_time_criterion_source_requested <-  BayesMVP::if_null_then_set_to(settings$n_iter_sampling_for_time_criterion_source,
                                                                                              ps7_time_criterion_default_settings$n_iter_sampling_for_time_criterion_source)
        ##
        if (!is.logical(run_sampling_timing_probe_requested) || length(run_sampling_timing_probe_requested) != 1 ||
            is.na(run_sampling_timing_probe_requested)) {
            stop("settings$run_sampling_timing_probe must be one TRUE or FALSE.")
        }
        if (!is.numeric(sampling_timing_probe_L_values_requested) || anyNA(sampling_timing_probe_L_values_requested) ||
            any(!is.finite(sampling_timing_probe_L_values_requested)) || any(sampling_timing_probe_L_values_requested < 1) ||
            any(sampling_timing_probe_L_values_requested != floor(sampling_timing_probe_L_values_requested)) ||
            length(unique(sampling_timing_probe_L_values_requested)) < 2) {
            stop("settings$sampling_timing_probe_L_values must hold at least two different whole numbers of leapfrog steps >= 1, e.g. c(2, 8).")
        }
        ## (NicoStan times a short call of 3 iterations at each number of leapfrog steps as well, to remove the per-call set-up, so the long call needs more:)
        if (!is.numeric(sampling_timing_probe_n_iter_per_L_requested) || length(sampling_timing_probe_n_iter_per_L_requested) != 1 ||
            !is.finite(sampling_timing_probe_n_iter_per_L_requested) || sampling_timing_probe_n_iter_per_L_requested < 4 ||
            sampling_timing_probe_n_iter_per_L_requested != floor(sampling_timing_probe_n_iter_per_L_requested)) {
            stop("settings$sampling_timing_probe_n_iter_per_L must be one whole number >= 4.")
        }
        if (!is.character(n_iter_sampling_for_time_criterion_source_requested) ||
            length(n_iter_sampling_for_time_criterion_source_requested) != 1 ||
            !n_iter_sampling_for_time_criterion_source_requested %in% c("planned_n_iter", "ESS_target_for_time_criterion")) {
            stop("settings$n_iter_sampling_for_time_criterion_source must be 'planned_n_iter' or 'ESS_target_for_time_criterion'.")
        }
        time_criterion_previous_run_path_requested <-  settings$time_criterion_previous_run_path
        if (!is.null(time_criterion_previous_run_path_requested) && length(time_criterion_previous_run_path_requested) == 1 &&
            is.na(time_criterion_previous_run_path_requested)) {
            time_criterion_previous_run_path_requested <-  NULL
        }
        if (!is.null(time_criterion_previous_run_path_requested) &&
            (!is.character(time_criterion_previous_run_path_requested) || length(time_criterion_previous_run_path_requested) < 1 ||
             anyNA(time_criterion_previous_run_path_requested) || any(!nzchar(time_criterion_previous_run_path_requested)))) {
            stop("settings$time_criterion_previous_run_path must be NULL or one or more paths of saved PS7 runs.")
        }
        ## "diagnostic" (default), "main", or one or more parameter names / base names (NULL or NA = the default):
        time_criterion_ess_parameter_set_requested <-  settings$time_criterion_ess_parameter_set
        if (is.null(time_criterion_ess_parameter_set_requested) ||
            (length(time_criterion_ess_parameter_set_requested) == 1 && is.na(time_criterion_ess_parameter_set_requested))) {
            time_criterion_ess_parameter_set_requested <-  ps7_time_criterion_default_settings$time_criterion_ess_parameter_set
        }
        if (!is.character(time_criterion_ess_parameter_set_requested) || length(time_criterion_ess_parameter_set_requested) < 1 ||
            anyNA(time_criterion_ess_parameter_set_requested) || any(!nzchar(trimws(time_criterion_ess_parameter_set_requested))) ||
            (length(time_criterion_ess_parameter_set_requested) > 1 && any(time_criterion_ess_parameter_set_requested %in% c("diagnostic", "main")))) {
            stop("settings$time_criterion_ess_parameter_set must be 'diagnostic', 'main' or one or more parameter names / base names, e.g. c('p', 'Se_baseline', 'Sp_baseline').")
        }
        ##
        return(list(run_sampling_timing_probe = run_sampling_timing_probe_requested,
                    ## sorted and without repeats, as NicoStan runs them, so that one probe always gives one file name:
                    sampling_timing_probe_L_values = sort(unique(as.numeric(sampling_timing_probe_L_values_requested))),
                    sampling_timing_probe_n_iter_per_L = as.numeric(sampling_timing_probe_n_iter_per_L_requested),
                    n_iter_sampling_for_time_criterion_source = n_iter_sampling_for_time_criterion_source_requested,
                    ESS_target_for_time_criterion = fn_supplied_number(settings$ESS_target_for_time_criterion, "ESS_target_for_time_criterion", TRUE),
                    ESS_per_iter_sampling_expected = fn_supplied_number(settings$ESS_per_iter_sampling_expected, "ESS_per_iter_sampling_expected", TRUE),
                    time_criterion_previous_run_path = time_criterion_previous_run_path_requested,
                    time_per_leapfrog_step_sampling = fn_supplied_number(settings$time_per_leapfrog_step_sampling, "time_per_leapfrog_step_sampling", TRUE),
                    time_per_iter_overhead_sampling = fn_supplied_number(settings$time_per_iter_overhead_sampling, "time_per_iter_overhead_sampling", FALSE),
                    time_per_iter_summaries_sampling = fn_supplied_number(settings$time_per_iter_summaries_sampling, "time_per_iter_summaries_sampling", FALSE),
                    time_criterion_ess_parameter_set = trimws(time_criterion_ess_parameter_set_requested)))

}
##
## ---- The time-criterion settings that actually reach the sampler: NULL when they are inert (any other algorithm, or tau
##      pinned), otherwise the validated values with every inert element removed (NULL):
##        - the previous run is kept only while it supplies something (a sampling time the user did not give, or
##          ESS_per_iter_sampling_expected under the "ESS_target_for_time_criterion" source);
##        - the probe runs only while a sampling time is still unknown; sampling_timing_probe_L_values and
##          sampling_timing_probe_n_iter_per_L are kept only while it runs;
##        - ESS_target_for_time_criterion and ESS_per_iter_sampling_expected only under the "ESS_target_for_time_criterion"
##          source, where a NULL ESS_target_for_time_criterion becomes the ps3 target min ESS (R_fn_ps7_get_target_min_ESS);
##        - time_criterion_ess_parameter_set only while ESS_per_iter_sampling_expected is read from the previous run(s) (the
##          "ESS_target_for_time_criterion" source with no user-supplied ESS_per_iter_sampling_expected).
##      BOTH the file-name builder and run_ps7_models() use this, so the name, the sampler arguments and the saved settings
##      agree (same design as fn_ps7_effective_tau_adaptation_block):
##
fn_ps7_effective_time_criterion_settings <-  function(settings,
                                                     N,
                                                     Model_type) {

        ## (compared as lower-case text, so a legacy tau_objective such as "ChEES_per_tau" never reaches the normaliser here):
        burnin_algorithm_key <-  if (is.null(settings$burnin_algorithm) || length(settings$burnin_algorithm) != 1) "" else
                                    tolower(trimws(as.character(settings$burnin_algorithm)))
        tau_is_pinned <-  (!is.null(settings$manual_L) && length(settings$manual_L) == 1 && !is.na(settings$manual_L)) ||
                         (!is.null(settings$manual_tau_value) && length(settings$manual_tau_value) == 1 && !is.na(settings$manual_tau_value))
        ## partitioned_HMC = TRUE: run_ps7_models() pins tau (manual_tau = TRUE), so the criterion is never computed there either:
        if (tau_is_pinned || isTRUE(settings$partitioned_HMC) ||
            !isTRUE(burnin_algorithm_key %in% c("chessr_time", "cheesr_time", "snaper_time"))) return(NULL)
        ##
        validated_time_criterion_settings <-  fn_ps7_validate_time_criterion_settings(settings = settings)
        effective_time_criterion_settings <-  validated_time_criterion_settings
        n_iter_sampling_uses_ESS_target <-  identical(validated_time_criterion_settings$n_iter_sampling_for_time_criterion_source, "ESS_target_for_time_criterion")
        ##
        if (!n_iter_sampling_uses_ESS_target) {
            effective_time_criterion_settings["ESS_target_for_time_criterion"] <-  list(NULL)
            effective_time_criterion_settings["ESS_per_iter_sampling_expected"] <-  list(NULL)
        } else if (is.null(validated_time_criterion_settings$ESS_target_for_time_criterion)) {
            ESS_target_for_this_N <-  R_fn_ps7_get_target_min_ESS(N = N, Model_type = Model_type)
            if (!is.finite(ESS_target_for_this_N) || ESS_target_for_this_N <= 0) {
                stop(paste0("n_iter_sampling_for_time_criterion_source = 'ESS_target_for_time_criterion' but there is no ps3 target min ESS for ",
                            Model_type, " at N = ", N, "; set ESS_target_for_time_criterion."))
            }
            effective_time_criterion_settings$ESS_target_for_time_criterion <-  ESS_target_for_this_N
        }
        ##
        sampling_times_user_supplied <-  !vapply(X = effective_time_criterion_settings[ps7_sampling_timing_quantity_names],
                                                FUN = is.null, FUN.VALUE = TRUE)
        previous_run_supplies_something <-  !is.null(validated_time_criterion_settings$time_criterion_previous_run_path) &&
                                           (!all(sampling_times_user_supplied) ||
                                            (n_iter_sampling_uses_ESS_target && is.null(effective_time_criterion_settings$ESS_per_iter_sampling_expected)))
        if (!previous_run_supplies_something) effective_time_criterion_settings["time_criterion_previous_run_path"] <-  list(NULL)
        ## the ESS parameter set matters only for ESS_per_iter_sampling_expected read from the previous run(s):
        ESS_per_iter_sampling_expected_from_previous_run <-  n_iter_sampling_uses_ESS_target &&
                                                            is.null(effective_time_criterion_settings$ESS_per_iter_sampling_expected) &&
                                                            !is.null(effective_time_criterion_settings$time_criterion_previous_run_path)
        if (!ESS_per_iter_sampling_expected_from_previous_run) effective_time_criterion_settings["time_criterion_ess_parameter_set"] <-  list(NULL)
        ##
        all_sampling_times_known_without_probe <-  all(sampling_times_user_supplied) ||
                                                  !is.null(effective_time_criterion_settings$time_criterion_previous_run_path)
        if (all_sampling_times_known_without_probe) effective_time_criterion_settings$run_sampling_timing_probe <-  FALSE
        if (!isTRUE(effective_time_criterion_settings$run_sampling_timing_probe)) {
            effective_time_criterion_settings["sampling_timing_probe_L_values"] <-  list(NULL)
            effective_time_criterion_settings["sampling_timing_probe_n_iter_per_L"] <-  list(NULL)
        }
        ##
        return(effective_time_criterion_settings)

}
##
## ---- The time-criterion settings text behind the "_bac<4 hex>" / "_bas<4 hex>" token: every EFFECTIVE non-default setting as
##      "name=value", in a fixed order, joined by ";" (NULL at the defaults and for every other algorithm). Numbers are written at
##      15 significant digits; the previous run by its BASENAME (the saved run's own name, which encodes its settings and seed):
##
fn_ps7_time_criterion_settings_text <-  function(effective_time_criterion_settings) {

        if (is.null(effective_time_criterion_settings)) return(NULL)
        fn_number_as_text_15_significant_digits <-  function(number_value) {
            paste(sub(pattern = "^ +", replacement = "", x = formatC(number_value, digits = 15, format = "g")), collapse = ",")
        }
        time_criterion_settings_text_parts <-  character(0)
        if (!isTRUE(effective_time_criterion_settings$run_sampling_timing_probe)) {
            time_criterion_settings_text_parts <-  c(time_criterion_settings_text_parts, "run_sampling_timing_probe=FALSE")
        } else {
            if (!identical(effective_time_criterion_settings$sampling_timing_probe_L_values,
                           ps7_time_criterion_default_settings$sampling_timing_probe_L_values)) {
                time_criterion_settings_text_parts <-  c(time_criterion_settings_text_parts,
                                                         paste0("sampling_timing_probe_L_values=",
                                                                fn_number_as_text_15_significant_digits(effective_time_criterion_settings$sampling_timing_probe_L_values)))
            }
            if (!identical(effective_time_criterion_settings$sampling_timing_probe_n_iter_per_L,
                           ps7_time_criterion_default_settings$sampling_timing_probe_n_iter_per_L)) {
                time_criterion_settings_text_parts <-  c(time_criterion_settings_text_parts,
                                                         paste0("sampling_timing_probe_n_iter_per_L=",
                                                                fn_number_as_text_15_significant_digits(effective_time_criterion_settings$sampling_timing_probe_n_iter_per_L)))
            }
        }
        if (identical(effective_time_criterion_settings$n_iter_sampling_for_time_criterion_source, "ESS_target_for_time_criterion")) {
            time_criterion_settings_text_parts <-  c(time_criterion_settings_text_parts,
                                                     "n_iter_sampling_for_time_criterion_source=ESS_target_for_time_criterion",
                                                     paste0("ESS_target_for_time_criterion=",
                                                            fn_number_as_text_15_significant_digits(effective_time_criterion_settings$ESS_target_for_time_criterion)))
        }
        if (!is.null(effective_time_criterion_settings$time_criterion_ess_parameter_set) &&
            !identical(effective_time_criterion_settings$time_criterion_ess_parameter_set,
                       ps7_time_criterion_default_settings$time_criterion_ess_parameter_set)) {
            time_criterion_settings_text_parts <-  c(time_criterion_settings_text_parts,
                                                     paste0("time_criterion_ess_parameter_set=",
                                                            paste(effective_time_criterion_settings$time_criterion_ess_parameter_set, collapse = ",")))
        }
        for (supplied_number_name in c("ESS_per_iter_sampling_expected", ps7_sampling_timing_quantity_names)) {
            supplied_value <-  effective_time_criterion_settings[[supplied_number_name]]
            if (!is.null(supplied_value)) {
                time_criterion_settings_text_parts <-  c(time_criterion_settings_text_parts,
                                                         paste0(supplied_number_name, "=", fn_number_as_text_15_significant_digits(supplied_value)))
            }
        }
        if (!is.null(effective_time_criterion_settings$time_criterion_previous_run_path)) {
            time_criterion_settings_text_parts <-  c(time_criterion_settings_text_parts,
                                                     paste0("time_criterion_previous_run=",
                                                            paste(sort(basename(effective_time_criterion_settings$time_criterion_previous_run_path)), collapse = ",")))
        }
        if (length(time_criterion_settings_text_parts) == 0) return(NULL)
        return(paste(time_criterion_settings_text_parts, collapse = ";"))

}
##
## ---- "" at the defaults (and for every other algorithm), otherwise the 4 hexadecimal digits that replace the final "t":
##
fn_ps7_time_criterion_settings_hash <-  function(effective_time_criterion_settings) {

        time_criterion_settings_text <-  fn_ps7_time_criterion_settings_text(effective_time_criterion_settings = effective_time_criterion_settings)
        if (is.null(time_criterion_settings_text)) return("")
        return(substr(digest::digest(object = time_criterion_settings_text, algo = "md5", serialize = FALSE), 1, 4))

}
##
##
## ---- The time-criterion record NicoStan returns for a CHESSR_time / SNAPER_time run: model_results$time_criterion (the settings
##      as validated, the resolved sampling quantities and their sources, the previous-run quantities, the sampling timing probe
##      result and wall time, and the burn-in record with the at_handover / at_last_update / at_end snapshots); else
##      burnin_object$time_criterion (the burn-in record alone); else a list called time_criterion_record. NULL for every other
##      algorithm.
##
fn_ps7_time_criterion_record_from_results <-  function(model_results) {

        for (candidate_record in list(model_results$time_criterion, model_results$burnin_object$time_criterion,
                                      model_results$time_criterion_record)) {
            if (is.list(candidate_record)) return(candidate_record)
        }
        return(NULL)

}
##
## ---- fn_ps7_time_criterion_values() (defined in fn_ps7_extract_tau_sweep.R, its main user) flattens this record into one value
##      per quantity of the naming list.
##
R_fn_enc_tau_wt <-  function(x) {
        if (!isTRUE(as.logical(x))) return("")   ## the original median-of-chains aggregation
        return("_tw1")
}
##
## A pinned trajectory length. "_mL9" = manual_tau with L = 9 leapfrog steps, adaptation OFF.
R_fn_enc_manual_L <-  function(x) {
        if (is.null(x) || length(x) != 1 || is.na(x)) return("")   ## adaptation ON, as before
        paste0("_mL", R_fn_enc_num(x))
}
##
## Absolute integration time, distinct from the existing fixed-L encoding.
R_fn_enc_manual_tau <-  function(x) {
        if (is.null(x)) return("")
        if (length(x) != 1) stop("manual_tau_value must be ONE value per run.")
        if (is.na(x)) return("")
        if (!is.numeric(x) || !is.finite(x) || x <= 0) stop("manual_tau_value must be positive and finite.")
        paste0("_mt", R_fn_enc_num(x))
}
##
R_fn_enc_map <-  function(x, map, field = "setting") {
        ## A filename encodes ONE run, so every setting reaching here must be a single value. A vector
        ## used to reach `if (is.na(out))` and die with "the condition has length > 1", which says
        ## nothing about which setting was wrong - and silently pasting both codes together would be
        ## worse, since it would build a filename that matches nothing.
        if (length(x) != 1L) {
            stop("R_fn_enc_map: '", field, "' must be ONE value per run, got ", length(x), ": ",
                 paste(as.character(x), collapse = ", "),
                 ". If this came from summarize_ps7_results(), the option lists are swept one combination at a",
                 " time - re-source ps_7_MCMC_settings_BayesMVP_functions.R, since an older copy of that",
                 " function passed the whole list through to the filename builder.")
        }
        x   <-  as.character(x)
        out <-  unname(map[x])
        if (is.na(out)) out <-  gsub("[^A-Za-z0-9]", "", x)   ## unknown value -> keep verbatim
        return(out)
}
##
R_fn_map_M_typ <-  c(Hessian = "H", 
                    Empirical = "E",
                    unit = "U", 
                    uniform_diag = "UD")
##
R_fn_map_M_shp <-  c(dense = "d", 
                    diag = "g")
##
R_fn_map_dHMC_int <-  c(kick_flow_kick = "kfk",    ## the sampler's default ordering
                       flow_kick_flow = "fkf")    ## Alenlov-Doucet-Lindsten
##
R_fn_map_tau_obj <-  c(KE = "ke",       ## original: change in kinetic energy (never encoded; see above)
                      ChEES = "ch", CHESSR_log = "ct", ChEES_per_tau = "ct")    ## Hoffman, Radul and Sountsov (2021), position-based
##
## ---- Collapse the axes that a PINNED trajectory length makes inert.
##
## With manual_L set, manual_tau is TRUE, and three of the swept settings then have no effect
## whatsoever on the sampler:
##
##   burnin_algorithm        - the adaptation is off, so nothing reads the objective.
##   tau_weight_by_p_jump - likewise.
##   tau_initial          - every use of it in R_fn_init_and_run_burnin_CHESS sits inside
##                          `if (manual_tau == FALSE)`, so under a pinned L it is dead.
##
## All three are still part of the FILE NAME, so leaving them crossed would run the identical
## sampler several times over and save the results under different names. That is not just wasted
## compute: the duplicates would come back as separate "configurations", so the spread across
## configurations - which is exactly what indicates whether a trajectory length is reliably better -
## would be computed over copies of the same run and look far tighter than it is.
##
## tau_initial is pinned to a HARD-CODED pi rather than to whatever the caller happened to list
## first, so the name a fixed-L run is saved under does not depend on the order or the contents of
## the tau_initial vector, and resuming works across sessions that set it differently.
##
## BOTH the runner (ps_7_MCMC_settings_BayesMVP.R) and the reader (summarize_ps7_results) call THIS
## function, so they cannot drift apart about which file names exist.
##
## ---- Restricting WHICH objectives get the p_jump weighting swept.
##
## weight_p_jump_only_for names the burnin_algorithm values for which tau_weight_by_p_jump = TRUE is
## worth running; every other objective is forced back to FALSE and de-duplicated. NULL (the
## default) sweeps the weighting for all of them, which is the plain cross.
##
## This controls an experimental weighting comparison only. KE uses the squared
## kinetic-energy-change criterion. The trajectory algorithm itself is selected
## solely by burnin_algorithm and is preserved in the run metadata and filename.
##
#' @param sampler_combinations   data frame with columns manual_L, burnin_algorithm,
#'                               tau_weight_by_p_jump and tau_initial.
#' @param weight_p_jump_only_for character vector of burnin_algorithm values for which
#'                               tau_weight_by_p_jump = TRUE is kept; NULL = keep it for all.
#' @return the same data frame with inert axes collapsed on pinned-L rows, de-duplicated.
#' @export
##
## ---- ADAM axes of the grid (eps / tau ADAM settings, the tau learning-rate restart and the tau adaptation scheme): each row is set to its
##      EFFECTIVE values (fn_ps7_effective_adam_settings: the tau settings are inert when tau is pinned, the restart with
##      "probe_then_average"), then the rows are de-duplicated, so the same run is never fitted twice under one name:
##
fn_ps7_collapse_inert_adam_axes <-  function(sampler_combinations) {

        adam_axis_names <-  intersect(x = c("eps_adam_beta1", "eps_adam_beta2", "eps_adam_epsilon", "tau_adam_beta1", "tau_adam_beta2", "tau_adam_epsilon",
                                            "tau_learning_rate_restart_at_metric_end", "tau_adaptation_scheme"),
                                      y = names(x = sampler_combinations))
        if (!length(x = adam_axis_names)) return(sampler_combinations)
        number_of_rows_before <-  nrow(x = sampler_combinations)
        for (row_index in seq_len(length.out = number_of_rows_before)) {
            effective_adam_settings <-  fn_ps7_effective_adam_settings(settings = as.list(sampler_combinations[row_index, , drop = FALSE]))
            for (adam_axis_name in adam_axis_names) sampler_combinations[[adam_axis_name]][row_index] <-  effective_adam_settings[[adam_axis_name]]
        }
        sampler_combinations <-  unique(x = sampler_combinations)
        rownames(x = sampler_combinations) <-  NULL
        number_of_rows_dropped <-  number_of_rows_before - nrow(x = sampler_combinations)
        if (number_of_rows_dropped > 0) {
            message(paste0("collapsed ", number_of_rows_dropped, " sampler combination(s): tau ADAM settings do nothing when tau is pinned, ",
                           "and the tau learning-rate restart does nothing with tau_adaptation_scheme = 'probe_then_average'."))
        }
        return(sampler_combinations)

}
##
fn_ps7_collapse_inert_tau_axes <-  function(sampler_combinations,
                                           weight_p_jump_only_for = NULL) {

        ##
        ## ---- ADAM axes (see fn_ps7_collapse_inert_adam_axes):
        ##
        sampler_combinations <-  fn_ps7_collapse_inert_adam_axes(sampler_combinations = sampler_combinations)
        ##
        ## ---- tau_sampling_scale is inert under jittered burn-in or a pinned tau (see fn_ps7_collapse_inert_tau_sampling_scale):
        ##
        sampler_combinations <-  fn_ps7_collapse_inert_tau_sampling_scale(sampler_combinations = sampler_combinations)
        ##
        ## ---- tau_adaptation_block is inert whenever tau is pinned (see fn_ps7_collapse_inert_tau_adaptation_block):
        ##
        sampler_combinations <-  fn_ps7_collapse_inert_tau_adaptation_block(sampler_combinations = sampler_combinations)
        ##
        if (!"manual_L" %in% names(x = sampler_combinations)) return(sampler_combinations)
        ##
        n_before <-  nrow(x = sampler_combinations)
        ##
        ## ---- restrict the p_jump sweep to the named objectives (before the pinned-L collapse, so
        ##      the two reductions compose and the row count is reported once, at the end)
        if (!is.null(weight_p_jump_only_for) &&
            all(c("burnin_algorithm", "tau_weight_by_p_jump") %in% names(x = sampler_combinations))) {
            ##
            ## unknown <-  setdiff(x = weight_p_jump_only_for, y = c("KE", "ChEES", "CHESSR", "CHESSR_log", "SNAPER"))
            ## unknown <-  setdiff(x = weight_p_jump_only_for, y = c("KE", "ChEES", "CHESSR", "CHESSR_log", "SNAPER", "CHESSR_time", "SNAPER_time"))
            unknown <-  setdiff(x = weight_p_jump_only_for, y = c("KE", "ChEES", "CHESSR", "CHESSR_log", "SNAPER", "CHESSR_time", "SNAPER_time",
                                                                  "ESJD", "ESJD_CHESSR", "ESJD_SNAPER", "LQ_ESSR"))
            if (length(x = unknown) > 0) {
                ## stop("weight_p_jump_only_for must contain only 'KE', 'ChEES', 'CHESSR', 'CHESSR_log' or 'SNAPER'; got: ",
                ## stop("weight_p_jump_only_for must contain only 'KE', 'ChEES', 'CHESSR', 'CHESSR_log', 'SNAPER', 'CHESSR_time' or 'SNAPER_time'; got: ",
                stop("weight_p_jump_only_for must contain only 'KE', 'ChEES', 'CHESSR', 'CHESSR_log', 'SNAPER', 'CHESSR_time', 'SNAPER_time', ",
                     "'ESJD', 'ESJD_CHESSR' or 'ESJD_SNAPER'; got: ",
                     paste(unknown, collapse = ", "))
            }
            not_weighted <-  !(sampler_combinations$burnin_algorithm %in% weight_p_jump_only_for)
            sampler_combinations$tau_weight_by_p_jump[not_weighted] <-  FALSE
        }
        ##
        if (!"manual_tau_value" %in% names(x = sampler_combinations)) sampler_combinations$manual_tau_value <-  NA_real_
        if (any(!is.na(sampler_combinations$manual_tau_value) &
                (!is.finite(sampler_combinations$manual_tau_value) | sampler_combinations$manual_tau_value <= 0))) {
            stop("manual_tau_value must contain NA or positive finite integration times.")
        }
        if (any(!is.na(sampler_combinations$manual_L) & !is.na(sampler_combinations$manual_tau_value))) {
            stop("Choose a fixed-L OR an absolute-tau sweep; set the other control to NA.")
        }
        pinned_rows <-  !is.na(sampler_combinations$manual_L) | !is.na(sampler_combinations$manual_tau_value)
        if (!any(pinned_rows)) {
            sampler_combinations <-  unique(x = sampler_combinations)
            n_dropped <-  n_before - nrow(x = sampler_combinations)
            if (n_dropped > 0) {
                message("dropped ", n_dropped, " sampler combination(s): tau_weight_by_p_jump is only ",
                        "swept for ", paste(weight_p_jump_only_for, collapse = "/"), ".")
            }
            return(sampler_combinations)
        }
        ##
        if ("burnin_algorithm" %in% names(x = sampler_combinations))        sampler_combinations$burnin_algorithm[pinned_rows] <-  "KE"
        if ("tau_weight_by_p_jump" %in% names(x = sampler_combinations)) sampler_combinations$tau_weight_by_p_jump[pinned_rows] <-  FALSE
        if ("tau_initial" %in% names(x = sampler_combinations))          sampler_combinations$tau_initial[pinned_rows] <-  pi
        ## a pinned tau is re-asserted after the burn-in ramp every iteration, so the ramp is inert as well:
        ## (NA = a run saved before the ramp option existed: no "_tr" token, so it must stay NA)
        if ("tau_ramp" %in% names(x = sampler_combinations))             sampler_combinations$tau_ramp[pinned_rows & !is.na(sampler_combinations$tau_ramp)] <-  "original"
        ##
        sampler_combinations <-  unique(x = sampler_combinations)
        ##
        n_dropped <-  n_before - nrow(x = sampler_combinations)
        if (n_dropped > 0) {
            message("collapsed ", n_dropped, " sampler combination(s): burnin_algorithm / ",
                    "tau_weight_by_p_jump / tau_initial do nothing when manual_L or manual_tau_value is set",
                    if (!is.null(weight_p_jump_only_for))
                        paste0(", and tau_weight_by_p_jump is only swept for ",
                               paste(weight_p_jump_only_for, collapse = "/")) else "", ".")
        }
        ##
        return(sampler_combinations)

}
##
##
## ---- tau_sampling_scale / randomize_tau_burnin / tau_adam_bias_correction: validation, effective values, collapse:
##
## tau_sampling_scale (NicoStan; an experimental heuristic derived from the unit-Gaussian
## analysis of Hoffman, Radul and Sountsov 2021 - NOT a published method; see NicoStan docs/adaptation-notes.md)
## multiplies the adapted tau once at the switch to sampling. NicoStan applies it ONLY when tau was adapted
## (manual_tau = FALSE) with a FIXED length (randomize_tau_burnin = FALSE) and sampling is randomised. So it is
## INERT when:
##   - randomize_tau_burnin = TRUE: tau_bar was already adapted under jitter (with share_tau_ii_across_chains_in_burnin
##     = TRUE this is ONE jittered tau per iteration shared by all chains - the published jitter scheme, whose mean
##     cost equals the fixed tau);
##   - manual_L / manual_tau_value is pinned, or partitioned_HMC = TRUE (ps7 then fixes tau with manual_tau = TRUE).
## Inert rows are collapsed to "none" and de-duplicated, as fn_ps7_collapse_inert_tau_axes() does for the pinned-L axes,
## so the same physical run is never saved under two names.
##
## tau_adam_bias_correction: "performed_update_counter" (default; the  NicoStan fix - ADAM bias correction counts
## the tau updates actually performed, "_ab1" in the file name) or "iteration_index" (every earlier run: the iteration
## index was used; no token, so those files keep their names). It is not a sampler argument: NicoStan records the
## scheme it used in model_results$tau_adam_bias_correction (NULL for builds older than the fix), which is checked.
##
## tau_adaptation_block (NicoStan - EXPERIMENTAL test arm, added to re-test the claim
## that adapting the trajectory length on the main block only stops a large nuisance block dominating the adaptation;
## NOT published): "main" (default; the trajectory-length criterion is computed on the main parameters only, the
## existing behaviour, no token) or "joint" (the criterion is computed on the main AND nuisance parameters
## concatenated, still adapting the single joint tau; "_tbJ" in the file name). It IS a sampler argument (passed
## next to randomize_tau_burnin / tau_sampling_scale / burnin_algorithm), and NicoStan also reports the EFFECTIVE
## value in model_results$tau_adaptation_block (NULL for builds older than this option, treated as "main"), which is
## checked. Inert whenever tau is pinned (manual_L or manual_tau_value set): no criterion is computed at all, so it
## collapses to "main" and de-duplicates, exactly as tau_sampling_scale collapses (see
## fn_ps7_collapse_inert_tau_adaptation_block() below).
##
fn_ps7_validate_tau_jitter_settings <-  function(settings) {

        tau_sampling_scale_requested <-  if_null_then_set_to(settings$tau_sampling_scale, "none")
        randomize_tau_burnin_requested <-  if_null_then_set_to(settings$randomize_tau_burnin, FALSE)
        tau_adam_bias_correction_requested <-  if_null_then_set_to(settings$tau_adam_bias_correction, "performed_update_counter")
        tau_adaptation_block_requested <-  if_null_then_set_to(settings$tau_adaptation_block, "main")
        ##
        tau_sampling_scale_is_valid_string <-  is.character(tau_sampling_scale_requested) &&
                                              length(tau_sampling_scale_requested) == 1 &&
                                              tau_sampling_scale_requested %in% c("none", "gaussian_matched")
        tau_sampling_scale_is_valid_number <-  is.numeric(tau_sampling_scale_requested) &&
                                              length(tau_sampling_scale_requested) == 1 &&
                                              is.finite(tau_sampling_scale_requested) &&
                                              tau_sampling_scale_requested > 0
        if (!tau_sampling_scale_is_valid_string && !tau_sampling_scale_is_valid_number) {
            stop(paste0("settings$tau_sampling_scale must be 'none', 'gaussian_matched' or one positive number; got: ",
                        paste(as.character(tau_sampling_scale_requested), collapse = ", ")))
        }
        if (!is.logical(randomize_tau_burnin_requested) || length(randomize_tau_burnin_requested) != 1 ||
            is.na(randomize_tau_burnin_requested)) {
            stop("settings$randomize_tau_burnin must be one TRUE or FALSE.")
        }
        if (!is.character(tau_adam_bias_correction_requested) || length(tau_adam_bias_correction_requested) != 1 ||
            !tau_adam_bias_correction_requested %in% c("performed_update_counter", "iteration_index")) {
            stop("settings$tau_adam_bias_correction must be 'performed_update_counter' or 'iteration_index'.")
        }
        if (!is.character(tau_adaptation_block_requested) || length(tau_adaptation_block_requested) != 1 ||
            is.na(tau_adaptation_block_requested) || !tau_adaptation_block_requested %in% c("main", "joint")) {
            stop("settings$tau_adaptation_block must be 'main' or 'joint'.")
        }
        ##
        return(list(tau_sampling_scale = tau_sampling_scale_requested,
                    randomize_tau_burnin = randomize_tau_burnin_requested,
                    tau_adam_bias_correction = tau_adam_bias_correction_requested,
                    tau_adaptation_block = tau_adaptation_block_requested))

}
##
## The tau_sampling_scale NicoStan will actually apply for these settings ("none" whenever it is inert, see above).
## BOTH the file-name builder and run_ps7_models() use this, so the name, the sampler argument and the self-check agree.
##
fn_ps7_effective_tau_sampling_scale <-  function(settings) {

        validated_tau_jitter_settings <-  fn_ps7_validate_tau_jitter_settings(settings = settings)
        tau_is_pinned <-  (!is.null(settings$manual_L) && length(settings$manual_L) == 1 && !is.na(settings$manual_L)) ||
                         (!is.null(settings$manual_tau_value) && length(settings$manual_tau_value) == 1 && !is.na(settings$manual_tau_value))
        tau_sampling_scale_is_inert <-  isTRUE(validated_tau_jitter_settings$randomize_tau_burnin) ||
                                       tau_is_pinned ||
                                       isTRUE(settings$partitioned_HMC)
        if (tau_sampling_scale_is_inert) return("none")
        return(validated_tau_jitter_settings$tau_sampling_scale)

}
##
## Unit-Gaussian "gaussian_matched" factor, recomputed HERE independently of NicoStan (base R only), so the factor the
## sampler reports is checked against the derivation rather than against itself. Same formulas as NicoStan's
## fn_tau_sampling_scale_gaussian_factor(): KE / ChEES 0.7151, CHESSR / SNAPER 0.7678, CHESSR_log 0.6738.
## CHESSR_time / SNAPER_time use the CHESSR / SNAPER factor 0.7678. It is the exact unit-Gaussian factor when
## burnin_to_sampling_leapfrog_time_ratio = 0 and sampling_overhead_in_leapfrog_steps = 0 (where they ARE CHESSR / SNAPER). With a
## positive tau_offset_from_sampling_overhead their optimum lies between the rate optimum and the ChEES optimum (whose factor is
## 0.7151), and a positive burnin_to_sampling_leapfrog_time_ratio moves it towards shorter tau, so 0.7678 is an approximation there.
##
fn_ps7_expected_gaussian_tau_sampling_factor <-  function(burnin_algorithm) {

        fn_first_maximiser <-  function(criterion_function) {
                optimize(f = criterion_function,
                         interval = c(0.3, 2.0),
                         maximum = TRUE,
                         tol = 1e-10)$maximum
        }
        fn_expected_chees_jittered <-  function(tau_bar) 0.5 - sin(4 * tau_bar) / (8 * tau_bar)
        fn_chees_rate_fixed <-  function(tau_fixed) sin(tau_fixed)^2 / tau_fixed
        fn_expected_per_chain_rate_jittered <-  function(tau_bar) {
                integrate(f = function(tau_value) ifelse(tau_value == 0, 0, sin(tau_value)^2 / tau_value),
                          lower = 0,
                          upper = 2 * tau_bar,
                          rel.tol = 1e-10)$value / (2 * tau_bar)
        }
        ##
        if (burnin_algorithm %in% c("KE", "ChEES")) {
                return(fn_first_maximiser(fn_expected_chees_jittered) / (pi / 2))
        }
        ## if (burnin_algorithm %in% c("CHESSR", "SNAPER")) {
        if (burnin_algorithm %in% c("CHESSR", "SNAPER", "CHESSR_time", "SNAPER_time")) {
                return(fn_first_maximiser(fn_expected_per_chain_rate_jittered) / fn_first_maximiser(fn_chees_rate_fixed))
        }
        if (burnin_algorithm == "CHESSR_log") {
                return(fn_first_maximiser(function(tau_bar) fn_expected_chees_jittered(tau_bar) / tau_bar) /
                       fn_first_maximiser(fn_chees_rate_fixed))
        }
        ##
        ## ---- ESJD (squared jump 2 (1 - cos(t))) 0.7678 and ESJD_CHESSR (geometric mean of the ESJD and CHESSR rates) 0.8083, both
        ##      per-chain means of C_i / tau_i; the fixed-length ESJD optimum (2.3311) lies beyond 2.0, so ESJD takes (0.3, 4.0), which
        ##      holds its first local maximum (fixed and jittered) and no other:
        ##      ESJD_SNAPER has the same unit-Gaussian factor as ESJD_CHESSR because the unit-Gaussian SNAPER rate is the per-chain rate above:
        ##
        fn_squared_jump_rate_fixed <-  function(tau_fixed) 2 * (1 - cos(tau_fixed)) / tau_fixed
        fn_expected_per_chain_squared_jump_rate_jittered <-  function(tau_bar) {
                stats::integrate(f = function(tau_value) ifelse(tau_value == 0, 0, 2 * (1 - cos(tau_value)) / tau_value),
                                 lower = 0,
                                 upper = 2 * tau_bar,
                                 rel.tol = 1e-10)$value / (2 * tau_bar)
        }
        fn_first_maximiser_up_to_four <-  function(criterion_function) {
                stats::optimize(f = criterion_function,
                                interval = c(0.3, 4.0),
                                maximum = TRUE,
                                tol = 1e-10)$maximum
        }
        if (burnin_algorithm == "ESJD") {
                return(fn_first_maximiser_up_to_four(fn_expected_per_chain_squared_jump_rate_jittered) /
                       fn_first_maximiser_up_to_four(fn_squared_jump_rate_fixed))
        }
        if (burnin_algorithm %in% c("ESJD_CHESSR", "ESJD_SNAPER")) {
                return(fn_first_maximiser(function(tau_bar) sqrt(fn_expected_per_chain_squared_jump_rate_jittered(tau_bar) * fn_expected_per_chain_rate_jittered(tau_bar))) /
                       fn_first_maximiser(function(tau_fixed) sqrt(fn_squared_jump_rate_fixed(tau_fixed) * fn_chees_rate_fixed(tau_fixed))))
        }
        ## LQ_ESSR: min(g(a), g(b)) / tau on a unit Gaussian, a = E[cos(t)], b = E[cos(t)^2], g(r) = (1 - r) / (1 + r)
        ## (the same factor as NicoStan's fn_tau_sampling_scale_gaussian_factor):
        if (burnin_algorithm == "LQ_ESSR") {
                fn_ESS_fraction <- function(lag_one_autocorrelation) (1 - lag_one_autocorrelation) / (1 + lag_one_autocorrelation)
                return(fn_first_maximiser(function(tau_bar) min(fn_ESS_fraction(sin(2 * tau_bar) / (2 * tau_bar)),
                                                                fn_ESS_fraction(0.5 + sin(4 * tau_bar) / (8 * tau_bar))) / tau_bar) /
                       fn_first_maximiser(function(tau_fixed) min(fn_ESS_fraction(cos(tau_fixed)),
                                                                  fn_ESS_fraction(cos(tau_fixed)^2)) / tau_fixed))
        }
        stop(paste0("No gaussian_matched factor for burnin_algorithm = '", burnin_algorithm, "'."))

}
##
## Collapse tau_sampling_scale to "none" on the rows where it is inert, then de-duplicate. Called at the start of
## fn_ps7_collapse_inert_tau_axes(), so the runner and every reader that uses that helper collapse it identically.
##
fn_ps7_collapse_inert_tau_sampling_scale <-  function(sampler_combinations) {

        if (!"tau_sampling_scale" %in% names(x = sampler_combinations)) return(sampler_combinations)
        ##
        number_of_rows_before <-  nrow(x = sampler_combinations)
        for (row_index in seq_len(length.out = number_of_rows_before)) {
            row_settings <-  as.list(sampler_combinations[row_index, , drop = FALSE])
            sampler_combinations$tau_sampling_scale[row_index] <-  fn_ps7_effective_tau_sampling_scale(settings = row_settings)
        }
        sampler_combinations <-  unique(x = sampler_combinations)
        rownames(x = sampler_combinations) <-  NULL
        ##
        number_of_rows_dropped <-  number_of_rows_before - nrow(x = sampler_combinations)
        if (number_of_rows_dropped > 0) {
            message(paste0("collapsed ", number_of_rows_dropped, " sampler combination(s): tau_sampling_scale does nothing when ",
                           "randomize_tau_burnin = TRUE (tau_bar already adapted under jitter) or tau is pinned."))
        }
        ##
        return(sampler_combinations)

}
##
## The tau_adaptation_block NicoStan will actually apply for these settings ("main" whenever it is inert, i.e. tau is
## pinned - no trajectory-length criterion is computed at all). BOTH the file-name builder and run_ps7_models() use
## this, so the name, the sampler argument and the self-check agree (same design as fn_ps7_effective_tau_sampling_scale).
##
fn_ps7_effective_tau_adaptation_block <-  function(settings) {

        validated_tau_jitter_settings <-  fn_ps7_validate_tau_jitter_settings(settings = settings)
        tau_is_pinned <-  (!is.null(settings$manual_L) && length(settings$manual_L) == 1 && !is.na(settings$manual_L)) ||
                         (!is.null(settings$manual_tau_value) && length(settings$manual_tau_value) == 1 && !is.na(settings$manual_tau_value))
        if (tau_is_pinned) return("main")
        return(validated_tau_jitter_settings$tau_adaptation_block)

}
##
## Collapse tau_adaptation_block to "main" on the rows where it is inert (tau pinned), then de-duplicate. Called at
## the start of fn_ps7_collapse_inert_tau_axes(), so the runner and every reader that uses that helper collapse it
## identically - exactly as fn_ps7_collapse_inert_tau_sampling_scale() does for tau_sampling_scale.
##
fn_ps7_collapse_inert_tau_adaptation_block <-  function(sampler_combinations) {

        if (!"tau_adaptation_block" %in% names(x = sampler_combinations)) return(sampler_combinations)
        ##
        number_of_rows_before <-  nrow(x = sampler_combinations)
        for (row_index in seq_len(length.out = number_of_rows_before)) {
            row_settings <-  as.list(sampler_combinations[row_index, , drop = FALSE])
            sampler_combinations$tau_adaptation_block[row_index] <-  fn_ps7_effective_tau_adaptation_block(settings = row_settings)
        }
        sampler_combinations <-  unique(x = sampler_combinations)
        rownames(x = sampler_combinations) <-  NULL
        ##
        number_of_rows_dropped <-  number_of_rows_before - nrow(x = sampler_combinations)
        if (number_of_rows_dropped > 0) {
            message(paste0("collapsed ", number_of_rows_dropped, " sampler combination(s): tau_adaptation_block does nothing when ",
                           "tau is pinned (manual_L or manual_tau_value set) - no trajectory-length criterion is computed."))
        }
        ##
        return(sampler_combinations)

}
##
##
## ---- ADAM settings of the burn-in adaptation: validation, effective values and the "_A" file-name code ---------------------------------
##
## The step-size (eps) update and the trajectory-length (tau) update each have their own ADAM settings (runner settings block, ps7_settings):
##   eps_adam_beta1, eps_adam_beta2, eps_adam_epsilon   first-moment decay, second-moment decay and denominator constant of the eps update
##                                                      (NicoStan options NicoStan_eps_adam_beta1 / _beta2 / _epsilon);
##   tau_adam_beta1, tau_adam_beta2, tau_adam_epsilon   the same for the tau update (NicoStan_tau_adam_beta1 / _beta2 / _epsilon);
##   tau_learning_rate_restart_at_metric_end            TRUE = the tau learning rate restarts at the full LR at the tau update of iteration
##                                                      metric_adaptation_end_iter and decays to LR^2 again by the last tau update
##                                                      (NicoStan_tau_learning_rate_restart_at_metric_end).
## The four tau settings are inert when tau is pinned (manual_L, manual_tau_value, or partitioned_HMC = TRUE, where ps7 fixes tau), so
## they collapse to their defaults there. A NULL setting (e.g. a settings list saved before these settings existed) means the value of
## every earlier run: beta1 0, beta2 0.95, denominator constant 1e-8, no restart.
##
## File names. At those defaults nothing changes ("_ab1" as before). Otherwise the token "_A<codes>" is written INSTEAD of "_ab1", in
## its place (after "_rb1" / "_ts..", before "_tbJ"): many names are already at the 249-byte limit, and "_A.." replaces the 4 bytes of
## "_ab1". Codes, in this order, each only when its value is not the default:
##   "i"         tau ADAM bias correction by iteration index (tau_adam_bias_correction = "iteration_index"); without "i" the "_A"
##               token means the performed-update counter, exactly as "_ab1";
##   "h<value>"  eps_adam_beta1;   "g<value>"  eps_adam_beta2;   "k<value>"  eps_adam_epsilon;
##   "t<value>"  tau_adam_beta1;   "u<value>"  tau_adam_beta2;   "w<value>"  tau_adam_epsilon;   "r"  the tau learning-rate restart;
##   "p"         tau_adaptation_scheme = "probe_then_average" (NicoStan_tau_adaptation_scheme; default "adam_decay", no code): the probe from
##               tau_handover / 8, then ADAM with constant LR and the averaged tau for sampling. The restart has no effect with it, so it
##               collapses to FALSE there. Inert (collapses to "adam_decay") when tau is pinned.
##   "q"         tau_adaptation_scheme = "fixed_length_probe_then_decay_and_average": as "p", with fixed-length trajectories during the
##               probe (also in the jittered arm), the top rung at tau_handover (not x 4) and the LR decaying from LR to LR^2 after the
##               probe. The restart collapses to FALSE and the scheme to "adam_decay" as for "p".
## Values are written as R_fn_enc_num() writes them (0.9 -> ".9", 1e-06 -> "1e-06"); no code letter is "e" or "E", so none can occur in
## a value. e.g. "_At.95" = tau beta1 0.95; "_Ah.9t.95" = eps beta1 0.9 and tau beta1 0.95; "_At.9u.99r" = tau beta1 0.9, tau beta2
## 0.99 and the restart. fn_ps7_parse_run_name() reads the codes back (and the codes "b" / "v" / "s", which set both updates, of names
## written before the two updates had separate settings) and translates "_A.." to "_ab1" (or nothing, with "i") before its other patterns.
##
fn_ps7_effective_adam_settings <-  function(settings) {

        fn_value_or_default <-  function(value, default_value) {
                if (is.null(value) || (length(value) == 1 && is.na(value))) return(default_value)
                return(value)
        }
        fn_check_decay <-  function(value, setting_name) {
                if (!is.numeric(value) || length(value) != 1 || !is.finite(value) || value < 0 || value >= 1) {
                    stop(paste0("settings$", setting_name, " must be one number in [0, 1); got: ", paste(as.character(value), collapse = ", ")))
                }
                return(as.numeric(value))
        }
        fn_check_constant <-  function(value, setting_name) {
                if (!is.numeric(value) || length(value) != 1 || !is.finite(value) || value <= 0) {
                    stop(paste0("settings$", setting_name, " must be one positive number; got: ", paste(as.character(value), collapse = ", ")))
                }
                return(as.numeric(value))
        }
        ##
        eps_adam_beta1   <-  fn_check_decay(fn_value_or_default(settings$eps_adam_beta1, 0),         "eps_adam_beta1")
        eps_adam_beta2   <-  fn_check_decay(fn_value_or_default(settings$eps_adam_beta2, 0.95),      "eps_adam_beta2")
        eps_adam_epsilon <-  fn_check_constant(fn_value_or_default(settings$eps_adam_epsilon, 1e-8), "eps_adam_epsilon")
        tau_adam_beta1   <-  fn_check_decay(fn_value_or_default(settings$tau_adam_beta1, 0),         "tau_adam_beta1")
        tau_adam_beta2   <-  fn_check_decay(fn_value_or_default(settings$tau_adam_beta2, 0.95),      "tau_adam_beta2")
        tau_adam_epsilon <-  fn_check_constant(fn_value_or_default(settings$tau_adam_epsilon, 1e-8), "tau_adam_epsilon")
        tau_learning_rate_restart_at_metric_end <-  fn_value_or_default(settings$tau_learning_rate_restart_at_metric_end, FALSE)
        if (!is.logical(tau_learning_rate_restart_at_metric_end) || length(tau_learning_rate_restart_at_metric_end) != 1 ||
            is.na(tau_learning_rate_restart_at_metric_end)) {
            stop("settings$tau_learning_rate_restart_at_metric_end must be one TRUE or FALSE.")
        }
        ##
        tau_adaptation_scheme <-  fn_value_or_default(settings$tau_adaptation_scheme, "adam_decay")
        # if (!is.character(tau_adaptation_scheme) || length(tau_adaptation_scheme) != 1 || !tau_adaptation_scheme %in% c("adam_decay", "probe_then_average")) {
        #     stop("settings$tau_adaptation_scheme must be 'adam_decay' or 'probe_then_average'.")
        if (!is.character(tau_adaptation_scheme) || length(tau_adaptation_scheme) != 1 ||
            !tau_adaptation_scheme %in% c("adam_decay", "probe_then_average", "fixed_length_probe_then_decay_and_average")) {
            stop("settings$tau_adaptation_scheme must be 'adam_decay', 'probe_then_average' or 'fixed_length_probe_then_decay_and_average'.")
        }
        ## the learning-rate restart has no effect with the constant learning rate of "probe_then_average":
        ## (nor with "fixed_length_probe_then_decay_and_average", whose post-probe LR schedule replaces it)
        # if (identical(tau_adaptation_scheme, "probe_then_average")) tau_learning_rate_restart_at_metric_end <-  FALSE
        if (tau_adaptation_scheme %in% c("probe_then_average", "fixed_length_probe_then_decay_and_average")) tau_learning_rate_restart_at_metric_end <-  FALSE
        ##
        tau_is_pinned <-  (!is.null(settings$manual_L) && length(settings$manual_L) == 1 && !is.na(settings$manual_L)) ||
                         (!is.null(settings$manual_tau_value) && length(settings$manual_tau_value) == 1 && !is.na(settings$manual_tau_value)) ||
                         isTRUE(settings$partitioned_HMC)
        if (tau_is_pinned) {
            tau_adam_beta1 <-  0
            tau_adam_beta2 <-  0.95
            tau_adam_epsilon <-  1e-8
            tau_learning_rate_restart_at_metric_end <-  FALSE
            tau_adaptation_scheme <-  "adam_decay"
        }
        ##
        return(list(eps_adam_beta1 = eps_adam_beta1,
                    eps_adam_beta2 = eps_adam_beta2,
                    eps_adam_epsilon = eps_adam_epsilon,
                    tau_adam_beta1 = tau_adam_beta1,
                    tau_adam_beta2 = tau_adam_beta2,
                    tau_adam_epsilon = tau_adam_epsilon,
                    tau_learning_rate_restart_at_metric_end = tau_learning_rate_restart_at_metric_end,
                    tau_adaptation_scheme = tau_adaptation_scheme))

}
##
## The codes of the "_A" token (see above), without the "_A" and the "i"; "" at the defaults of every earlier run:
##
fn_ps7_adam_file_name_code <-  function(settings) {

        effective_adam_settings <-  fn_ps7_effective_adam_settings(settings = settings)
        return(paste0(if (effective_adam_settings$eps_adam_beta1 != 0)      paste0("h", R_fn_enc_num(effective_adam_settings$eps_adam_beta1)) else "",
                      if (effective_adam_settings$eps_adam_beta2 != 0.95)   paste0("g", R_fn_enc_num(effective_adam_settings$eps_adam_beta2)) else "",
                      if (effective_adam_settings$eps_adam_epsilon != 1e-8) paste0("k", R_fn_enc_num(effective_adam_settings$eps_adam_epsilon)) else "",
                      if (effective_adam_settings$tau_adam_beta1 != 0)      paste0("t", R_fn_enc_num(effective_adam_settings$tau_adam_beta1)) else "",
                      if (effective_adam_settings$tau_adam_beta2 != 0.95)   paste0("u", R_fn_enc_num(effective_adam_settings$tau_adam_beta2)) else "",
                      if (effective_adam_settings$tau_adam_epsilon != 1e-8) paste0("w", R_fn_enc_num(effective_adam_settings$tau_adam_epsilon)) else "",
                      if (isTRUE(effective_adam_settings$tau_learning_rate_restart_at_metric_end)) "r" else "",
                      # if (identical(effective_adam_settings$tau_adaptation_scheme, "probe_then_average")) "p" else ""))
                      if (identical(effective_adam_settings$tau_adaptation_scheme, "probe_then_average")) "p" else "",
                      if (identical(effective_adam_settings$tau_adaptation_scheme, "fixed_length_probe_then_decay_and_average")) "q" else ""))

}
##
##
## ---- eps_acceptance_mean / tau_shrink_on_divergence: validation and the file-name version token -------------------------------------
##
## eps_acceptance_mean (NicoStan burn-in): how the per-chain Metropolis acceptance probabilities of one burn-in iteration are
## averaged before the ADAM step-size update (gradient p_jump - adapt_delta on log eps):
##   "harmonic"   (default) K / sum_k (1 / p_jump_k), the harmonic mean of Hoffman, Radul and Sountsov (2021, Algorithm 1)
##                and Sountsov and Hoffman (2022, Appendix, eq. 15), so that every chain sharing the step size can make progress;
##   "arithmetic" mean(p_jump_k), the rule every earlier ps7 run used;
##   "geometric"  exp(mean(log(p_jump_k))), each p_jump_k clamped to [1e-8, 1] first; between the harmonic and the arithmetic mean.
## tau_shrink_on_divergence (NicoStan burn-in): whether tau_main and tau_us are multiplied by 0.95 at an adapting iteration
## (clip_iter < ii < n_adapt, tau not pinned) with sum(div_main) > 1 over the burn-in chains, i.e. at least two divergent chains:
##   FALSE (default) no divergence-triggered shrink (the trajectory length follows the adaptation criterion only);
##   TRUE            the shrink every earlier ps7 run had (NicoStan's own factor and chain threshold are left at their defaults).
## Both ARE sampler arguments (passed next to tau_adaptation_block), and NicoStan reports the values it applied in
## model_results (NULL for builds older than these options, treated as "arithmetic" / TRUE), which is checked after every fit.
##
## File names. The pair is written as the "_ta" version digit, so it costs NO bytes (many names are already at the 249-byte
## limit, and any extra token would push the N = 10,000 names with "_tsG" / "_tbJ" over it even after the "_LiLR" shortening).
## Every earlier development ("_ta4") run used "arithmetic" with the shrink on, so that pair keeps "_ta4" and exactly the name
## it has on disk today (those runs still resume under the old behaviour). Any other pair gets its own digit, so a run of the
## new default can never be resumed from, or saved over, an earlier "_ta4" file:
##   arithmetic + shrink on  -> "_ta4"   (every earlier run)
##   harmonic   + shrink off -> "_ta5"   (the new default)
##   arithmetic + shrink off -> "_ta6"
##   harmonic   + shrink on  -> "_ta7"
##   geometric  + shrink off -> "_ta8"
##   geometric  + shrink on  -> "_ta9"
## fn_ps7_parse_run_name() reads the pair back from the digit ("_ta4" and every earlier version = arithmetic / shrink on).
## A legacy settings list (no burnin_algorithm, "_ta2" names) has no version to carry the pair, so there a missing value
## means the earlier behaviour ("arithmetic" / TRUE) and any other pair stops in R_fn_file_name_string().
##
fn_ps7_validate_eps_tau_shrink_settings <-  function(settings) {

        legacy_settings <-  is.null(settings$burnin_algorithm)
        eps_acceptance_mean_requested <-  BayesMVP::if_null_then_set_to(settings$eps_acceptance_mean, if (legacy_settings) "arithmetic" else "harmonic")
        tau_shrink_on_divergence_requested <-  BayesMVP::if_null_then_set_to(settings$tau_shrink_on_divergence, legacy_settings)
        ##
        # if (!is.character(eps_acceptance_mean_requested) || length(eps_acceptance_mean_requested) != 1 ||
        #     is.na(eps_acceptance_mean_requested) || !eps_acceptance_mean_requested %in% c("harmonic", "arithmetic")) {
        #     stop("settings$eps_acceptance_mean must be 'harmonic' or 'arithmetic'.")
        # }
        if (!is.character(eps_acceptance_mean_requested) || length(eps_acceptance_mean_requested) != 1 ||
            is.na(eps_acceptance_mean_requested) || !eps_acceptance_mean_requested %in% c("harmonic", "arithmetic", "geometric")) {
            stop("settings$eps_acceptance_mean must be 'harmonic', 'arithmetic' or 'geometric'.")
        }
        if (!is.logical(tau_shrink_on_divergence_requested) || length(tau_shrink_on_divergence_requested) != 1 ||
            is.na(tau_shrink_on_divergence_requested)) {
            stop("settings$tau_shrink_on_divergence must be one TRUE or FALSE.")
        }
        ##
        return(list(eps_acceptance_mean = eps_acceptance_mean_requested,
                    tau_shrink_on_divergence = tau_shrink_on_divergence_requested))

}
##
## "_ta4" / "_ta5" / "_ta6" / "_ta7" / "_ta8" / "_ta9" for the six pairs (see above):
##
fn_ps7_tau_adaptation_version_token <-  function(settings) {

        validated_eps_tau_shrink_settings <-  fn_ps7_validate_eps_tau_shrink_settings(settings = settings)
        eps_is_harmonic <-  identical(validated_eps_tau_shrink_settings$eps_acceptance_mean, "harmonic")
        eps_is_geometric <-  identical(validated_eps_tau_shrink_settings$eps_acceptance_mean, "geometric")
        tau_shrink_is_on <-  isTRUE(validated_eps_tau_shrink_settings$tau_shrink_on_divergence)
        if ( eps_is_geometric && !tau_shrink_is_on) return("_ta8")
        if ( eps_is_geometric &&  tau_shrink_is_on) return("_ta9")
        if (!eps_is_harmonic &&  tau_shrink_is_on) return("_ta4")
        if ( eps_is_harmonic && !tau_shrink_is_on) return("_ta5")
        if (!eps_is_harmonic && !tau_shrink_is_on) return("_ta6")
        return("_ta7")

}
##
##
## ---- metric_pooled_window_resets / metric_pooled_offdiagonal_shrinkage: validation and the file-name code ---------------------------
##
## Both are NicoStan burn-in settings of the POOLED metric estimator (metric_estimator = "pooled") only; "chain_mean" and
## "chain_mean_scaled" never read them (their covariance is shrunk every iteration by Rcpp_shrink_matrix(shrinkage_factor)).
## "per_iteration" reads metric_pooled_offdiagonal_shrinkage (its dense main covariance proposal is shrunk in the same way) but not
## metric_pooled_window_resets (nothing is accumulated across iterations); its file-name code is fn_ps7_metric_per_iteration_file_name_code.
##
## metric_pooled_window_resets: iterations at which the pooled Welford accumulators are reset. These accumulators hold the
## pooled covariance of the dense main metric AND the pooled variances of the nuisance diagonal metric (metric_type_nuisance
## "Empirical" or "uniform_diag"), so the resets change both:
##   "stan_style"  (default) round(c(0.30, 0.60) * n_adapt), the resets every earlier "_Mp" run had (the final metric then uses
##                 only the last window, e.g. iterations 136 to 180 of a 250-iteration burn-in);
##   "none"        no reset: ONE window from the metric start to metric_adaptation_end_iter, so the final metric uses every
##                 pooled draw in that range;
##   c(r1, r2, ..) whole numbers >= 0; the accumulators are reset at the start of iteration r + 1 for each value r.
## metric_pooled_offdiagonal_shrinkage: s in [0, 1], applied ONCE to the pooled main covariance proposal, after the Stan-style
## regularisation and before the symmetrisation, as (1 - s) * cov + s * diag(diag(cov)) (the nuisance metric is diagonal):
##   0 (default) = every off-diagonal unshrunk, as in every earlier "_Mp" run; 1 = diagonal.
## Both ARE sampler arguments (passed next to metric_estimator), and NicoStan reports the values it applied in model_results
## (NULL for builds older than these options, treated as "stan_style" / 0), which is checked after every fit.
##
## A missing value means the earlier behaviour ("stan_style" / 0), the same as NicoStan's own defaults, for every settings
## list; other values are set explicitly in the runner. A numeric vector of reset iterations is returned sorted and without
## duplicates (the same resets), so that one set of resets always has one file name.
##
fn_ps7_validate_metric_pooled_settings <-  function(settings) {

        metric_pooled_window_resets_requested <-  if (is.null(settings$metric_pooled_window_resets)) {
            "stan_style"
        } else settings$metric_pooled_window_resets
        ##
        metric_pooled_offdiagonal_shrinkage_requested <-  if (is.null(settings$metric_pooled_offdiagonal_shrinkage)) {
            0
        } else settings$metric_pooled_offdiagonal_shrinkage
        ##
        ## ---- window resets: "none", "stan_style" or whole numbers >= 0 (a zero-length vector = "none"):
        ##
        if (is.character(metric_pooled_window_resets_requested)) {
            if (length(metric_pooled_window_resets_requested) != 1 || is.na(metric_pooled_window_resets_requested) ||
                !metric_pooled_window_resets_requested %in% c("none", "stan_style")) {
                stop(paste0("settings$metric_pooled_window_resets must be 'none', 'stan_style' or a numeric vector of whole numbers >= 0; got: ",
                            paste(as.character(metric_pooled_window_resets_requested), collapse = ", ")))
            }
        } else if (is.numeric(metric_pooled_window_resets_requested) && !is.factor(metric_pooled_window_resets_requested)) {
            if (length(metric_pooled_window_resets_requested) == 0) {
                metric_pooled_window_resets_requested <-  "none"
            } else {
                if (any(!is.finite(metric_pooled_window_resets_requested)) || any(metric_pooled_window_resets_requested < 0) ||
                    any(metric_pooled_window_resets_requested != round(metric_pooled_window_resets_requested))) {
                    stop(paste0("settings$metric_pooled_window_resets: a numeric value must hold whole numbers >= 0 (iterations); got: ",
                                paste(as.character(metric_pooled_window_resets_requested), collapse = ", ")))
                }
                metric_pooled_window_resets_requested <-  sort(unique(as.numeric(metric_pooled_window_resets_requested)))
            }
        } else {
            stop("settings$metric_pooled_window_resets must be 'none', 'stan_style' or a numeric vector of whole numbers >= 0.")
        }
        ##
        ## ---- off-diagonal shrinkage: one number in [0, 1]:
        ##
        adaptive_shrinkage_requested <-  identical(metric_pooled_offdiagonal_shrinkage_requested, "adaptive")
        if (is.character(metric_pooled_offdiagonal_shrinkage_requested) && !adaptive_shrinkage_requested &&
            length(metric_pooled_offdiagonal_shrinkage_requested) == 1 && !is.na(metric_pooled_offdiagonal_shrinkage_requested)) {
            metric_pooled_offdiagonal_shrinkage_requested <-  suppressWarnings(as.numeric(metric_pooled_offdiagonal_shrinkage_requested))
        }
        if (!adaptive_shrinkage_requested &&
            (!is.numeric(metric_pooled_offdiagonal_shrinkage_requested) || is.factor(metric_pooled_offdiagonal_shrinkage_requested) ||
            length(metric_pooled_offdiagonal_shrinkage_requested) != 1 || !is.finite(metric_pooled_offdiagonal_shrinkage_requested) ||
            metric_pooled_offdiagonal_shrinkage_requested < 0 || metric_pooled_offdiagonal_shrinkage_requested > 1)) {
            stop(paste0("settings$metric_pooled_offdiagonal_shrinkage must be one number in [0, 1] or 'adaptive'; got: ",
                        paste(as.character(metric_pooled_offdiagonal_shrinkage_requested), collapse = ", ")))
        }
        ##
        return(list(metric_pooled_window_resets = metric_pooled_window_resets_requested,
                    metric_pooled_offdiagonal_shrinkage = if (adaptive_shrinkage_requested) "adaptive" else
                                                        as.numeric(metric_pooled_offdiagonal_shrinkage_requested)))

}
##
## One printable label per set of window resets: "none", "stan_style", or the reset iterations joined by "-" (e.g. "68-135").
## The same label is the metric_pooled_window_resets column of fn_ps7_parse_run_name(). A zero-length vector is "none":
##
fn_ps7_metric_pooled_window_resets_label <-  function(metric_pooled_window_resets) {

        if (is.character(metric_pooled_window_resets)) return(paste(metric_pooled_window_resets, collapse = "-"))
        if (length(metric_pooled_window_resets) == 0) return("none")
        return(paste(sub(pattern = "^ +", replacement = "", x = formatC(sort(unique(as.numeric(metric_pooled_window_resets))), format = "f", digits = 0)),
                     collapse = "-"))

}
##
## File-name code appended to "_Mp" (pooled estimator only). "" for "stan_style" resets with off-diagonal shrinkage 0, the pair
## of every earlier "_Mp" run, which therefore keeps exactly the name it has on disk (and still resumes under that behaviour).
## Any other pair is written as <shrinkage><window code>, with no separator and no extra underscore:
##   shrinkage   R_fn_enc_num(s), e.g. ".5", "0", ".25", "1";
##   window code "w" = "none" (one whole window), "s" = "stan_style", "r" + the reset iterations joined by "-" (e.g. "r68-135").
## e.g. "_Mp.5w" = "none" / 0.5, "_Mp0w" = no resets and no shrinkage, "_Mp.5s" = Stan-style resets with 0.5.
## "_Mp.5w" costs 3 bytes: many names are within a few bytes of the 249-byte limit, so a separate token for each
## setting (at least 9 bytes) would push the longest saved "_Mp" configurations over it.
## fn_ps7_parse_run_name() reads both settings back ("_Mp" alone, and every non-pooled name, = "stan_style" / 0).
##
fn_ps7_metric_pooled_file_name_code <-  function(settings) {

        if (!identical(settings$metric_estimator, "pooled")) return("")
        ##
        validated_metric_pooled_settings <-  fn_ps7_validate_metric_pooled_settings(settings = settings)
        metric_pooled_window_resets_validated <-  validated_metric_pooled_settings$metric_pooled_window_resets
        metric_pooled_offdiagonal_shrinkage_validated <-  validated_metric_pooled_settings$metric_pooled_offdiagonal_shrinkage
        ##
        if (identical(metric_pooled_window_resets_validated, "stan_style") &&
            is.numeric(metric_pooled_offdiagonal_shrinkage_validated) && metric_pooled_offdiagonal_shrinkage_validated == 0) return("")
        ##
        metric_pooled_window_code <-  if (identical(metric_pooled_window_resets_validated, "none")) "w" else
                                     if (identical(metric_pooled_window_resets_validated, "stan_style")) "s" else
                                     paste0("r", fn_ps7_metric_pooled_window_resets_label(metric_pooled_window_resets = metric_pooled_window_resets_validated))
        ##
        shrinkage_file_name_code <-  if (identical(metric_pooled_offdiagonal_shrinkage_validated, "adaptive")) "A" else
                                       R_fn_enc_num(metric_pooled_offdiagonal_shrinkage_validated)
        return(paste0(shrinkage_file_name_code, metric_pooled_window_code))

}
##
## File-name code appended to "_Mi" (per-iteration estimator only): the off-diagonal shrinkage, which the per-iteration estimator
## applies to its dense main covariance proposal as the pooled estimator does. "" for 0, otherwise R_fn_enc_num(s), e.g. "_Mi1" = 1
## (diagonal), "_Mi.5" = 0.5. The window resets are inert for it and are not written. fn_ps7_parse_run_name() reads the shrinkage back.
##
fn_ps7_metric_per_iteration_file_name_code <-  function(settings) {

        if (!identical(settings$metric_estimator, "per_iteration")) return("")
        ##
        validated_metric_pooled_settings <-  fn_ps7_validate_metric_pooled_settings(settings = settings)
        metric_pooled_offdiagonal_shrinkage_validated <-  validated_metric_pooled_settings$metric_pooled_offdiagonal_shrinkage
        ##
        if (is.numeric(metric_pooled_offdiagonal_shrinkage_validated) && metric_pooled_offdiagonal_shrinkage_validated == 0) return("")
        ##
        return(if (identical(metric_pooled_offdiagonal_shrinkage_validated, "adaptive")) "A" else
                   R_fn_enc_num(metric_pooled_offdiagonal_shrinkage_validated))

}
##
R_fn_map_M_dcy <-  c(inverse = "inv",
                    exponential = "exp",
                    constant = "con",
                    none = "non")
##
##
## Resolve before BOTH resume lookup and fitting. Missing schedule keeps old callers' filenames unchanged.
fn_ps7_resolve_burnin_schedule <-  function(settings) {
        if (is.null(x = settings$burnin_schedule)) {
            settings$resolved_burnin_schedule <-  NULL
            if (is.null(x = settings$metric_adaptation_end_iter) && is.null(x = settings$n_adapt)) return(settings)
            ## Explicit cutoffs still need a token, even for an older caller that did not select a schedule.
            settings$burnin_schedule <-  "legacy"
        }
        resolved_schedule <-  fn_burnin_adaptation_schedule(
            n_burnin = settings$n_burnin,
            n_adapt = settings$n_adapt,
            burnin_schedule = settings$burnin_schedule,
            metric_adaptation_end_iter = settings$metric_adaptation_end_iter,
            theta_hat_us_freeze_iter = settings$theta_hat_us_freeze_iter,
            theta_hat_us_rule = if (is.null(x = settings$theta_hat_us_rule)) "running_mean_frozen" else settings$theta_hat_us_rule,
            clip_iter = settings$clip_iter,
            clip_iter_tau = if (is.null(x = settings$int) || is.null(x = settings$clip_iter)) NULL else settings$clip_iter + settings$int)
        settings$n_adapt <-  resolved_schedule$n_adapt
        settings$metric_adaptation_end_iter <-  resolved_schedule$metric_adaptation_end_iter
        settings$theta_hat_us_freeze_iter <-  resolved_schedule$theta_hat_us_freeze_iter
        if (!is.null(x = resolved_schedule$clip_iter)) settings$clip_iter <-  resolved_schedule$clip_iter
        if (!is.null(x = resolved_schedule$clip_iter_tau)) settings$int <-  resolved_schedule$clip_iter_tau - settings$clip_iter
        settings$resolved_burnin_schedule <-  resolved_schedule
        settings
}
##

##
## ---- Recorded trajectory-criterion options and non-default filename tokens -----------------------------------------------------------
##
fn_ps7_recorded_trajectory_options <-  function(settings) {

        recorded_legacy_values <-  list(tau_gradient_estimator = "forward",
                                       tau_cost_exponent = 1,
                                       esjd_jump_power = 2,
                                       tau_jitter_burnin = "uniform",
                                       ## LQ_ESSR only: parameter families the criterion monitors (NULL = all rows; "_io" file token):
                                       interest_only = NULL)
        options <-  lapply(names(recorded_legacy_values), function(option_name) {
                if (is.null(settings[[option_name]])) recorded_legacy_values[[option_name]] else settings[[option_name]]
        })
        names(options) <-  names(recorded_legacy_values)
        if (!is.character(options$tau_gradient_estimator) || length(options$tau_gradient_estimator) != 1 ||
            is.na(options$tau_gradient_estimator) || !options$tau_gradient_estimator %in% c("forward", "two_ended")) {
                stop("tau_gradient_estimator must be 'forward' or 'two_ended'.")
        }
        if (!is.numeric(options$tau_cost_exponent) || length(options$tau_cost_exponent) != 1 ||
            !is.finite(options$tau_cost_exponent) || options$tau_cost_exponent < 0 || options$tau_cost_exponent > 1.5) {
                stop("tau_cost_exponent must be a finite number in [0, 1.5].")
        }
        if (!is.numeric(options$esjd_jump_power) || length(options$esjd_jump_power) != 1 ||
            !is.finite(options$esjd_jump_power) || !options$esjd_jump_power %in% c(2, 3, 4)) {
                stop("esjd_jump_power must be 2, 3 or 4.")
        }
        if (!is.character(options$tau_jitter_burnin) || length(options$tau_jitter_burnin) != 1 ||
            is.na(options$tau_jitter_burnin) || !options$tau_jitter_burnin %in% c("uniform", "halton")) {
                stop("tau_jitter_burnin must be 'uniform' or 'halton'.")
        }
        algorithm <-  if (is.null(settings$burnin_algorithm)) "KE" else fn_ps7_normalise_burnin_algorithm(settings$burnin_algorithm)
        if (options$tau_cost_exponent != 1 && algorithm %in% c("KE", "ChEES", "CHESSR_log", "CHESSR_time", "SNAPER_time")) {
                stop("This criterion requires the legacy tau_cost_exponent value 1.")
        }
        if (options$esjd_jump_power != 2 && !algorithm %in% c("ESJD", "ESJD_CHESSR", "ESJD_SNAPER")) {
                stop("esjd_jump_power is applied only by the ESJD family; other criteria require 2.")
        }
        if (algorithm == "KE" && options$tau_gradient_estimator != "forward") stop("KE requires the forward trajectory-criterion estimator.")
        if (!is.null(options$interest_only)) {
                if (!is.character(options$interest_only) || length(options$interest_only) < 1 || anyNA(options$interest_only)) {
                        stop("interest_only must be NULL or a character vector of parameter family names, e.g. c(\"beta\", \"p_raw\").")
                }
                if (algorithm != "LQ_ESSR") stop("interest_only is used only by burnin_algorithm = 'LQ_ESSR'.")
                ## the third choice of the criterion's parameters (main = the 44 main parameters, joint = every parameter):
                if (identical(settings$tau_adaptation_block, "joint")) stop("interest_only is its own choice of parameters; run it with tau_adaptation_block = 'main'.")
        }
        if (options$tau_jitter_burnin == "halton" && !isTRUE(settings$randomize_tau_burnin)) {
                stop("Halton burn-in jitter requires randomize_tau_burnin = TRUE.")
        }
        options$tau_cost_exponent <-  as.numeric(options$tau_cost_exponent)
        options$esjd_jump_power <-  as.numeric(options$esjd_jump_power)
        return(options)

}




fn_ps7_trajectory_option_arguments <-  function(settings) {

        options <-  fn_ps7_recorded_trajectory_options(settings)
        # changed <-  c(options$tau_gradient_estimator != "forward",
        #               options$tau_cost_exponent != 1,
        #               options$esjd_jump_power != 2,
        #               options$tau_jitter_burnin != "uniform")
        # return(options[changed])
        ##
        ## ---- all four options are passed every time: R_fn_sample_model is internal and has no defaults, so a run with the
        ##      default values must still pass them (the file name carries only the non-default ones):
        return(options)

}




fn_ps7_trajectory_option_file_name_code <-  function(settings) {

        options <-  fn_ps7_recorded_trajectory_options(settings)
        for (precision in seq_len(17)) {
              cost_exponent_code <-  trimws(formatC(options$tau_cost_exponent, format = "g", digits = precision))
              if (identical(as.numeric(cost_exponent_code), as.numeric(options$tau_cost_exponent))) break
        }
        cost_exponent_code <-  sub(pattern = "^0\\.", replacement = ".", x = cost_exponent_code)
        cost_exponent_code <-  sub(pattern = "e([-+])0+([0-9]+)$", replacement = "e\\1\\2", x = cost_exponent_code)
        return(paste0(if (options$tau_gradient_estimator != "forward") "_tg2" else "",
                      if (options$tau_cost_exponent != 1) paste0("_tc", cost_exponent_code) else "",
                      if (options$esjd_jump_power != 2) paste0("_jb", options$esjd_jump_power) else "",
                      if (options$tau_jitter_burnin != "uniform") "_tjH" else "",
                      ## interest_only: "_io" + the first letter of each family, sorted (c("beta", "p_raw") -> "_iobp"):
                      if (!is.null(options$interest_only)) paste0("_io", paste(substr(sort(options$interest_only), 1, 1), collapse = "")) else ""))

}



R_fn_file_name_string <-  function(Model_type,
                                  N,
                                  settings,
                                  model_args_list,
                                  prior_LKJ_nd,
                                  prior_LKJ_d,
                                  prior_prev_a,
                                  prior_prev_b,
                                  grouping = NULL
) {

        development_algorithm <-  !is.null(settings$burnin_algorithm)
        if (development_algorithm) {
                settings$burnin_algorithm <-  fn_ps7_normalise_burnin_algorithm(burnin_algorithm = settings$burnin_algorithm)
        } else {
                settings$burnin_algorithm <-  if (is.null(settings$tau_objective)) "KE" else settings$tau_objective
        }
        ##
        if (!is.null(settings$burnin_post_adapt_iter)) {
            stop("burnin_post_adapt_iter has been removed; filenames now describe the full n_burnin length.")
        }
        ##
        ## ---- A pinned trajectory length makes burnin_algorithm, tau_weight_by_p_jump and tau_initial
        ## inert (every use of them in the sampler sits behind `manual_tau == FALSE`). Collapsing them
        ## HERE as well as in fn_ps7_collapse_inert_tau_axes() makes the name correct by construction:
        ## the same physical run gets the same name no matter which caller built the settings, so a
        ## direct call cannot mint a name like "_toch_mL12" that nothing will ever be saved under.
        ## The name stays honest, because it encodes what the SAMPLER does, not what was passed.
        invisible(R_fn_enc_manual_tau(settings$manual_tau_value))   ## validate before scalar conditions below
        if ((!is.null(settings$manual_L) && length(settings$manual_L) == 1 && !is.na(settings$manual_L)) ||
            (!is.null(settings$manual_tau_value) && length(settings$manual_tau_value) == 1 && !is.na(settings$manual_tau_value))) {
            if (!is.null(settings$manual_L) && !is.na(settings$manual_L) &&
                !is.null(settings$manual_tau_value) && !is.na(settings$manual_tau_value)) {
                stop("manual_L and manual_tau_value cannot both be fixed.")
            }
            settings$burnin_algorithm        <-  "KE"
            settings$tau_weight_by_p_jump <-  FALSE
            settings$tau_initial          <-  pi
            if (!is.null(settings$tau_ramp)) settings$tau_ramp <-  "original"
        }
        ##
        ## ---- A filename describes ONE run. If a swept option list (metric_estimator, tau_initial,
        ## learning_rate_initial, diffusion_HMC_integrator, ...) arrives here whole instead of one
        ## value at a time, the name silently becomes wrong or the encoders die with an unhelpful
        ## "the condition has length > 1". Name the offending field instead.
        ##
        scalar_settings <-  c("burnin_algorithm", "output_dir", "diffusion_HMC", "partitioned_HMC", "diffusion_HMC_integrator",
                             "n_chains_burnin", "n_chains_sampling", "n_threads_WCP_burnin", "n_threads_WCP_sampling",
                             "num_chunks_burnin", "num_chunks_sampling", "n_burnin", "n_iter", "learning_rate",
                             "adapt_delta", "tau_initial", "clip_iter", "int", "learning_rate_initial",
                             "ratio_M_main", "ratio_M_nuisance", "metric_type_main", "metric_type_nuisance",
                             "metric_shape_main", "int_width", "M_decay_type", "M_decay_power", "M_decay_scale",
                             "multi_attempts", "metric_estimator", "manual_tau_value", "tau_ramp",
                             "eps_reinit_at_ChEES_handover", "theta_hat_us_rule", "pre_burnin_n_iter",
                             "pre_burnin_L", "share_tau_ii_across_chains_in_burnin", "burnin_TBB_pool_equals_n_chains",
                             "J_grad_option", "autodiff_fallback", "theta_hat_us_freeze_iter", "store_log_lik_trace",
                             "burnin_schedule", "metric_adaptation_end_iter", "n_adapt",
                             ## "tau_sampling_scale", "randomize_tau_burnin", "tau_adam_bias_correction", "tau_adaptation_block")
                             "tau_sampling_scale", "randomize_tau_burnin", "tau_adam_bias_correction", "tau_adaptation_block",
                             ## (metric_pooled_window_resets is not here: one run may reset at several iterations; it is checked
                             ##  by fn_ps7_validate_metric_pooled_settings below)
                             ## "eps_acceptance_mean", "tau_shrink_on_divergence")
                             ## "eps_acceptance_mean", "tau_shrink_on_divergence", "metric_pooled_offdiagonal_shrinkage")
                             "eps_acceptance_mean", "tau_shrink_on_divergence", "metric_pooled_offdiagonal_shrinkage",
                             ## ADAM settings of the eps and the tau updates ("_A" code, fn_ps7_adam_file_name_code):
                             "eps_adam_beta1", "eps_adam_beta2", "eps_adam_epsilon", "tau_adam_beta1", "tau_adam_beta2", "tau_adam_epsilon",
                             "tau_learning_rate_restart_at_metric_end", "tau_adaptation_scheme",
                             "tau_gradient_estimator", "tau_cost_exponent", "esjd_jump_power", "tau_jitter_burnin",
                             ## time-criterion settings of CHESSR_time / SNAPER_time (sampling_timing_probe_L_values,
                             ## time_criterion_previous_run_path and time_criterion_ess_parameter_set may hold several values;
                             ## fn_ps7_validate_time_criterion_settings checks them):
                             "run_sampling_timing_probe", "sampling_timing_probe_n_iter_per_L", "n_iter_sampling_for_time_criterion_source",
                             "ESS_target_for_time_criterion", "ESS_per_iter_sampling_expected",
                             "time_per_leapfrog_step_sampling", "time_per_iter_overhead_sampling", "time_per_iter_summaries_sampling")
        for (field in scalar_settings) {
            value <-  settings[[field]]
            if (!is.null(value) && length(value) != 1L) {
                stop("R_fn_file_name_string: settings$", field, " must be ONE value per run, got ", length(value),
                     ": ", paste(as.character(value), collapse = ", "),
                     ".\n  summarize_ps7_results() sweeps the option lists one combination at a time; each option",
                     " reached this from the summary section, re-source",
                     " functions/ps_7_MCMC_settings_BayesMVP_functions.R - an older copy of that function passed",
                     " the whole list straight through to this builder.")
            }
        }
        ##
        ## ---- tau_initial: a number or "adaptive"; a numeric level that arrives as text (a character grid column) is restored exactly:
        if (!is.null(settings$tau_initial)) settings$tau_initial <-  fn_ps7_tau_initial_value(settings$tau_initial)
        ##
        if (!is.null(settings$store_log_lik_trace) &&
            (!is.logical(settings$store_log_lik_trace) || is.na(settings$store_log_lik_trace))) {
            stop("settings$store_log_lik_trace must be TRUE or FALSE (NULL/missing keeps the PS7 default, FALSE).")
        }
        ##
        if (!is.null(settings$autodiff_fallback) &&
            (!is.logical(settings$autodiff_fallback) || is.na(settings$autodiff_fallback))) {
            stop("settings$autodiff_fallback must be a single TRUE or FALSE.")
        }
        ##
        ## pooled metric estimator settings (inert for "chain_mean" / "chain_mean_scaled", but a bad value still stops here):
        invisible(fn_ps7_validate_metric_pooled_settings(settings = settings))
        ##
        settings <-  fn_ps7_resolve_burnin_schedule(settings = settings)
        theta_hat_us_freeze_iter <-  settings$theta_hat_us_freeze_iter
        if (!is.null(theta_hat_us_freeze_iter) &&
            (!is.numeric(theta_hat_us_freeze_iter) || !is.finite(theta_hat_us_freeze_iter) ||
             theta_hat_us_freeze_iter < 1 || theta_hat_us_freeze_iter > .Machine$integer.max ||
             theta_hat_us_freeze_iter != floor(theta_hat_us_freeze_iter))) {
            stop("settings$theta_hat_us_freeze_iter must be NULL or one positive integer iteration.")
        }
        ##
        file_name_string <-  file.path(settings$output_dir, paste0("ps7_run",
                                                                      ##
                                                                      "_", Model_type,
                                                                      ##
                                                                      "_N", N,
                                                                      ##
                                                                      "_np", model_args_list$n_pops,
                                                                      ## Corrected ramp, tau optimiser clock and squared-KE objective.
                                                                      ## Never resume a pre-fix run under the corrected adaptation.
                                                                      ## "_ta5" / "_ta6" / "_ta7" = the step-size acceptance mean and divergence tau shrink pairs
                                                                      ## (see fn_ps7_tau_adaptation_version_token); "_ta4" = the pair every earlier run used.
                                                                      ## if (development_algorithm) "_ta4" else "_ta2",
                                                                      if (development_algorithm) fn_ps7_tau_adaptation_version_token(settings = settings) else "_ta2",
                                                                      settings$resolved_burnin_schedule$filename_token,
                                                                      ##
                                                                      "_dH", R_fn_enc_lgl(settings$diffusion_HMC),
                                                                      "_pH", R_fn_enc_lgl(settings$partitioned_HMC),
                                                                      ##
                                                                      ## joint-diffusion integrator; "" when it is the default kick_flow_kick,
                                                                      ## which is what every run predating this field used:
                                                                      R_fn_enc_dHMC_int(settings$diffusion_HMC_integrator),
                                                                      ##
                                                                      ## trajectory-length adaptation; all "" at the original ChESSR settings:
                                                                      ## if (development_algorithm) fn_ps7_encode_burnin_algorithm(settings$burnin_algorithm) else
                                                                      ## (non-default time-criterion settings of CHESSR_time / SNAPER_time: "_bac<4 hex>" / "_bas<4 hex>")
                                                                      if (development_algorithm) fn_ps7_encode_burnin_algorithm(burnin_algorithm = settings$burnin_algorithm,
                                                                                                                                 time_criterion_settings_hash = fn_ps7_time_criterion_settings_hash(
                                                                                                                                     fn_ps7_effective_time_criterion_settings(settings = settings, N = N, Model_type = Model_type))) else
                                                                          R_fn_enc_tau_obj(settings$burnin_algorithm),
                                                                      R_fn_enc_tau_wt(settings$tau_weight_by_p_jump),
                                                                      R_fn_enc_manual_L(settings$manual_L),
                                                                      R_fn_enc_manual_tau(settings$manual_tau_value),
                                                                      ##
                                                                      "_cb", settings$n_chains_burnin,
                                                                      "_s", settings$n_chains_sampling, ## ----
                                                                      ##
                                                                      "_wb", settings$n_threads_WCP_burnin,
                                                                      "_s", settings$n_threads_WCP_sampling, ## ----
                                                                      ##
                                                                      "_kb", settings$num_chunks_burnin, ## ----
                                                                      "_s", settings$num_chunks_sampling,      ## ----
                                                                      ##
                                                                      "_b", settings$n_burnin,
                                                                      "_it", settings$n_iter,
                                                                      "_LR", R_fn_enc_num(settings$learning_rate),
                                                                      "_AD", R_fn_enc_num(settings$adapt_delta),
                                                                      ## "_tiA" = tau_initial "adaptive" (tau set at the handover from the eigenvalue estimate):
                                                                      ## "_ti", if (identical(settings$tau_initial, pi)) "pi" else
                                                                      ## "_ti", if (identical(settings$tau_initial, "adaptive")) "A" else
                                                                      ## "_tiAw.5" = adaptive with lambda_max from the last half of [clip_iter, handover] (no "w" = all of it):
                                                                      ## "_tiAall" = all iterations from 1 to the handover; historical window codes are unchanged.
                                                                      "_ti", if (identical(settings$tau_initial, "adaptive"))
                                                                                  if (identical(settings$tau_initial_moments_window_fraction, "all")) "Aall" else
                                                                                  paste0("A", if (!is.null(settings$tau_initial_moments_window_fraction) &&
                                                                                                    settings$tau_initial_moments_window_fraction != 1)
                                                                                                    paste0("w", R_fn_enc_num(settings$tau_initial_moments_window_fraction)) else "") else
                                                                          if (identical(settings$tau_initial, pi)) "pi" else
                                                                          if (identical(settings$tau_initial, 2*pi)) "2pi" else
                                                                          ## New schedule names need room for their boundary token; do not rename legacy pi/2 files.
                                                                          if (!is.null(settings$resolved_burnin_schedule) && identical(settings$tau_initial, pi/2)) "hpi" else
                                                                              R_fn_enc_opt(settings$tau_initial),
                                                                      "_clip", settings$clip_iter,
                                                                      "_int", settings$int,
                                                                      ##
                                                                      ## ---- learning-rate hold and step-size warm start. These CHANGE THE POSTERIOR
                                                                      ## adaptation, so they belong in the name; n_refresh deliberately does NOT, since
                                                                      ## it only controls how chatty the burnin is and would invalidate saved runs.
                                                                      ## "d" = sampler default, "off" = feature disabled.
                                                                      R_fn_enc_hold("_Li", settings$learning_rate_initial, settings$learning_rate_initial_iter),
                                                                      R_fn_enc_hold("_Ei", settings$eps_initial, settings$eps_initial_iter),
                                                                      "_rM", R_fn_enc_num(settings$ratio_M_main), "_", R_fn_enc_num(settings$ratio_M_nuisance),
                                                                      ##
                                                                      # "_M_typ", R_fn_enc_map(settings$metric_type_main,  R_fn_map_M_typ),
                                                                      "_mt", R_fn_enc_map(settings$metric_type_main,  R_fn_map_M_typ),
                                                                      "_nu", R_fn_enc_map(settings$metric_type_nuisance,  R_fn_map_M_typ),
                                                                      ##
                                                                      "_ms", R_fn_enc_map(settings$metric_shape_main, R_fn_map_M_shp),
                                                                      "_w", settings$int_width,
                                                                      ##
                                                                      "_md", R_fn_enc_map(settings$M_decay_type, R_fn_map_M_dcy),
                                                                      ##
                                                                      "_LKJ", prior_LKJ_nd, "_", prior_LKJ_d,
                                                                      "_pp", R_fn_enc_num(prior_prev_a), "_", prior_prev_b,
                                                                      ##
                                                                      "_ma", R_fn_enc_lgl(settings$multi_attempts)
        ))
        ##
        if (settings$M_decay_type == "inverse") {
          file_name_string <-  paste0(file_name_string,
                                     "_mp", R_fn_enc_num(settings$M_decay_power),
                                     "_msc", R_fn_enc_num(settings$M_decay_scale))
        }
        ##
        if (isTRUE(grouping)) {
          file_name_string <-  paste0(file_name_string, "_og1")
        }
        ##
        if (Model_type == "LC_MVP") { 
          if (settings$reorder_cols_MVP == FALSE) file_name_string <-  paste0(file_name_string, "_OC")
        }
        ##
        if (settings$metric_estimator == "pooled") { 
          ## "_Mp" + the pooled-estimator code (fn_ps7_metric_pooled_file_name_code): nothing for "stan_style" resets with
          ## off-diagonal shrinkage 0 (every earlier run, and any settings list without the two fields, whose name is therefore
          ## unchanged), e.g. "_Mp.5w" for "none" / 0.5:
          ## file_name_string <-  paste0(file_name_string, "_Mp")
          file_name_string <-  paste0(file_name_string, "_Mp", fn_ps7_metric_pooled_file_name_code(settings = settings))
        } else if (settings$metric_estimator == "chain_mean_scaled") {
          file_name_string <-  paste0(file_name_string, "_Ms")     ## chain-mean variance x n_chains_burnin
        }
        ## "_Mi" = "per_iteration" (the cross-chain (co)variance of each iteration's draws, blended into the metric with ratio_M), followed
        ## by the off-diagonal shrinkage when it is not 0 (fn_ps7_metric_per_iteration_file_name_code), e.g. "_Mi1". Any other value would
        ## be saved under a "chain_mean" name (no code), so it stops here:
        if (identical(settings$metric_estimator, "per_iteration")) {
          file_name_string <-  paste0(file_name_string, "_Mi", fn_ps7_metric_per_iteration_file_name_code(settings = settings))
        } else if (!settings$metric_estimator %in% c("pooled", "chain_mean", "chain_mean_scaled")) {
          stop(paste0("R_fn_file_name_string: metric_estimator must be 'pooled', 'chain_mean', 'chain_mean_scaled' or 'per_iteration'; got: ",
                      paste(as.character(settings$metric_estimator), collapse = ", ")))
        }
        ##
        ## A fixed (user-supplied) test order, e.g. "_tp546132"; absent = order estimated by the pre-burnin:
        if (!is.null(settings$test_perm_override)) {
          ## digits run together for < 10 tests (e.g. "_tp546132"); separated by "-" otherwise, so it stays unambiguous:
          test_perm_separator <-  if (max(settings$test_perm_override) < 10) "" else "-"
          file_name_string <-  paste0(file_name_string, "_tp", paste(settings$test_perm_override, collapse = test_perm_separator))
        }
        ##
        ## Burn-in tau ramp: "_trO" = original, "_trS" = staged. Absent = a run saved before the option existed
        ## (those used the staged ramp from  onwards and the original ramp before that):
        if (!is.null(settings$tau_ramp)) {
          file_name_string <-  paste0(file_name_string, "_tr", if (identical(settings$tau_ramp, "staged")) "S" else "O")
        }
        ##
        ## "_nER" = eps NOT re-initialised at the ChEES handover. Absent = re-initialised (the default, and every
        ## run saved before the option existed):
        if (isFALSE(settings$eps_reinit_at_ChEES_handover)) {
          file_name_string <-  paste0(file_name_string, "_nER")
        }
        ##
        ## Nuisance Gaussian centre: "_cF" = running mean frozen before eps finishes adapting (the default from
        ##, also when NULL), "_c0" = fixed at zero. Absent = "running_mean" (never frozen),
        ## which is what every run saved before the option existed used:
        theta_hat_us_rule_used <-  if (is.null(settings$theta_hat_us_rule)) "running_mean_frozen" else settings$theta_hat_us_rule
        if (identical(theta_hat_us_rule_used, "running_mean_frozen")) file_name_string <-  paste0(file_name_string, "_cF")
        if (identical(theta_hat_us_rule_used, "zero"))                file_name_string <-  paste0(file_name_string, "_c0")
        ## NULL keeps the existing name; every explicit freeze schedule has a distinct tag.
        if (!is.null(theta_hat_us_freeze_iter) && is.null(settings$resolved_burnin_schedule)) {
          file_name_string <-  paste0(file_name_string, "_cfi", paste0(sub("^ +", "", formatC(theta_hat_us_freeze_iter, format = "f", digits = 0))))
        }
        ##
        ## Burn-in speed-ups: "_pa<k>" = only k non-adapting iterations after n_adapt; "_pb<n>" = n-iteration pre-burnin.
        ## Absent = unchanged (all n_burnin - n_adapt of them; 125-iteration pre-burnin):
        ## Historical _pa notes above are retained for reference only; new filenames never generate that token.
        if (!is.null(settings$pre_burnin_n_iter))      file_name_string <-  paste0(file_name_string, "_pb", settings$pre_burnin_n_iter)
        if (!is.null(settings$pre_burnin_L))           file_name_string <-  paste0(file_name_string, "_pL", settings$pre_burnin_L)
        ## "_sT" = one tau_ii per burn-in iteration shared by all burn-in chains (absent = per-chain, as before):
        if (isTRUE(settings$share_tau_ii_across_chains_in_burnin)) file_name_string <-  paste0(file_name_string, "_sT")
        ## "_tbbC" = burn-in TBB pool of exactly n_chains_burnin threads (absent = n_threads_WCP_burnin x n_chains_burnin, as before):
        if (isTRUE(settings$burnin_TBB_pool_equals_n_chains))      file_name_string <-  paste0(file_name_string, "_tbbC")
        ## "_JgA" = correlation-transform Jacobian gradient by autodiff (absent = "num_diff", finite differences, as before):
        if (identical(settings$J_grad_option, "autodiff"))          file_name_string <-  paste0(file_name_string, "_JgA")
        ## "_ll1" = sampling stores the log-likelihood trace; FALSE retains the current PS7 filename.
        if (isTRUE(settings$store_log_lik_trace))                  file_name_string <-  paste0(file_name_string, "_ll1")
        ## "_af1" = optional same-position autodiff fallback. Absent/FALSE preserves the existing filename.
        if (isTRUE(settings$autodiff_fallback) && (is.null(settings$multi_attempts) || isTRUE(settings$multi_attempts)) &&
            Model_type %in% c("LC_MVP", "MVP", "LC_MVOP", "MVOP", "latent_trait")) {
          file_name_string <-  paste0(file_name_string, "_af1")
        }
        ##
        ## ---- Trajectory-jitter tokens, appended LAST (after "_af1") so every older name - which has none of them - still
        ##      parses unchanged (fn_ps7_parse_run_name strips them first). Absent = the behaviour before each option existed:
        ##        "_rb1"  = randomize_tau_burnin = TRUE (with "_sT": ONE jittered tau per iteration shared by all burn-in chains);
        ##        "_tsG"  = tau_sampling_scale = "gaussian_matched" EFFECTIVELY applied (a number x gives "_ts" + x, e.g. "_ts.75");
        ##                  never written when it is inert (see fn_ps7_effective_tau_sampling_scale);
        ##        "_ab1"  = tau ADAM bias correction counts PERFORMED updates (the  NicoStan fix). Absent = the iteration
        ##                  index, as in every earlier run - so pre-fix files are never resumed as post-fix runs.
        ##        "_tbJ"  = tau_adaptation_block = "joint" EFFECTIVELY applied (the trajectory-length criterion computed on the
        ##                  main AND nuisance parameters concatenated;, EXPERIMENTAL). Written
        ##                  AFTER "_ab1"; never written when it is inert (see fn_ps7_effective_tau_adaptation_block). Absent = "main"
        ##                  (the criterion computed on the main parameters only), the existing behaviour.
        ##
        {
          validated_tau_jitter_settings <-  fn_ps7_validate_tau_jitter_settings(settings = settings)
          tau_sampling_scale_effective <-  fn_ps7_effective_tau_sampling_scale(settings = settings)
          tau_adaptation_block_effective <-  fn_ps7_effective_tau_adaptation_block(settings = settings)
          ##
          if (isTRUE(validated_tau_jitter_settings$randomize_tau_burnin)) {
            file_name_string <-  paste0(file_name_string, "_rb1")
          }
          if (identical(tau_sampling_scale_effective, "gaussian_matched")) {
            file_name_string <-  paste0(file_name_string, "_tsG")
          } else if (is.numeric(tau_sampling_scale_effective)) {
            file_name_string <-  paste0(file_name_string, "_ts", R_fn_enc_num(tau_sampling_scale_effective))
          }
          ## ---- ADAM settings: "_A<codes>" INSTEAD of "_ab1" (in its place) when any of them is not the default of every earlier
          ##      run; "" at those defaults, so those names keep "_ab1" (see fn_ps7_adam_file_name_code):
          adam_file_name_code <-  fn_ps7_adam_file_name_code(settings = settings)
          # if (identical(validated_tau_jitter_settings$tau_adam_bias_correction, "performed_update_counter")) {
          if (identical(validated_tau_jitter_settings$tau_adam_bias_correction, "performed_update_counter") && !nzchar(adam_file_name_code)) {
            file_name_string <-  paste0(file_name_string, "_ab1")
          }
          if (nzchar(adam_file_name_code)) {
            file_name_string <-  paste0(file_name_string, "_A",
                                       if (identical(validated_tau_jitter_settings$tau_adam_bias_correction, "iteration_index")) "i" else "",
                                       adam_file_name_code)
          }
          if (identical(tau_adaptation_block_effective, "joint")) {
            file_name_string <-  paste0(file_name_string, "_tbJ")
          }
        }
        ##
        ## ---- Step-size acceptance mean / divergence-triggered tau shrink: carried by the "_ta" digit above (no bytes added).
        ##      A legacy "_ta2" name has no version to carry them, so only the pair of every earlier run is accepted there:
        ##
        if (!development_algorithm &&
            !identical(fn_ps7_tau_adaptation_version_token(settings = utils::modifyList(x = settings, val = list(burnin_algorithm = NULL))), "_ta4")) {
          stop("R_fn_file_name_string: a legacy (\"_ta2\") name cannot encode eps_acceptance_mean / tau_shrink_on_divergence; ",
               "set burnin_algorithm, or eps_acceptance_mean = 'arithmetic' with tau_shrink_on_divergence = TRUE.")
        }
        ##
        ## ---- Guard: NAME_MAX is 255 BYTES per path component. Reserve 6 for the "_runNN"
        ##      suffix. Fails HERE (before the resume check) rather than at saveRDS after a
        ##      run has already completed and is about to be thrown away:
        ##
        file_name_string <-  paste0(file_name_string, fn_ps7_trajectory_option_file_name_code(settings = settings))
        nb <-  nchar(basename(file_name_string), type = "bytes")
        ##
        ## ---- If the name is too long and the initial LR is just the LR (learning_rate_initial = NULL in the runner),
        ##      "_Li<LR>_d" is written as "_LiLR" (saves 5 bytes). Names that already fit (i.e., every run saved so far) are unchanged:
        ##
        if ((nb > 249) && is.null(settings$learning_rate_initial_iter) && (length(settings$learning_rate_initial) == 1) &&
            isTRUE(abs(as.numeric(settings$learning_rate_initial) - as.numeric(settings$learning_rate)) < 1e-12)) {
          file_name_string <-  sub(pattern = paste0("_Li", R_fn_enc_num(settings$learning_rate_initial), "_d_"), replacement = "_LiLR_",
                                   x = file_name_string, fixed = TRUE)
          nb <-  nchar(basename(file_name_string), type = "bytes")
        }
        ##
        ## ---- Still too long (e.g. an "_A" ADAM code on a name already at the limit): three lossless short forms, applied together and
        ##      only here, so every name that already fits (i.e. every run saved so far) is unchanged. fn_ps7_parse_run_name() reads
        ##      none of these segments:
        ##        "_rM<r>_<r>_"             -> "_rM<r>_"     when ratio_M_main = ratio_M_nuisance    (e.g. "_rM.9_.9_" -> "_rM.9_");
        ##        "_LKJ<a>_<a>_"            -> "_LKJ<a>_"    when prior_LKJ_nd = prior_LKJ_d          (e.g. "_LKJ4_4_"  -> "_LKJ4_");
        ##        "_md<inv|exp|con|non>_"   -> "_md<i|e|c|n>_"  (M_decay_type)                          (e.g. "_mdinv_"   -> "_mdi_").
        ##
        if (nb > 249) {
          if (isTRUE(settings$ratio_M_main == settings$ratio_M_nuisance)) {
            file_name_string <-  sub(pattern = paste0("_rM", R_fn_enc_num(settings$ratio_M_main), "_", R_fn_enc_num(settings$ratio_M_nuisance), "_"),
                                     replacement = paste0("_rM", R_fn_enc_num(settings$ratio_M_main), "_"),
                                     x = file_name_string, fixed = TRUE)
          }
          if (isTRUE(prior_LKJ_nd == prior_LKJ_d)) {
            file_name_string <-  sub(pattern = paste0("_LKJ", prior_LKJ_nd, "_", prior_LKJ_d, "_"),
                                     replacement = paste0("_LKJ", prior_LKJ_nd, "_"),
                                     x = file_name_string, fixed = TRUE)
          }
          M_decay_type_code <-  R_fn_enc_map(settings$M_decay_type, R_fn_map_M_dcy)
          M_decay_type_short_code <-  unname(c(inv = "i", exp = "e", con = "c", non = "n")[M_decay_type_code])
          if (length(M_decay_type_short_code) == 1 && !is.na(M_decay_type_short_code)) {
            file_name_string <-  sub(pattern = paste0("_md", M_decay_type_code, "_"),
                                     replacement = paste0("_md", M_decay_type_short_code, "_"),
                                     x = file_name_string, fixed = TRUE)
          }
          nb <-  nchar(basename(file_name_string), type = "bytes")
        }
        ##
        ## ---- Still too long after those three (e.g. "_tc.75_jb3" with "_LR.0125"): one more lossless short form, again applied only
        ##      here, so every name that already fits is unchanged. fn_ps7_parse_run_name() does not read this segment either:
        ##        "_dH<d>_pH<p>_"           -> "_dH<d>pH<p>_"  (diffusion_HMC / partitioned_HMC flags)  (e.g. "_dH1_pH0_" -> "_dH1pH0_").
        ##
        if (nb > 249) {
          file_name_string <-  sub(pattern = "_dH([01])_pH([01])_", replacement = "_dH\\1pH\\2_", x = file_name_string)
          nb <-  nchar(basename(file_name_string), type = "bytes")
        }
        if (nb > 249) {
          stop(paste0("R_fn_file_name_string: basename is ", formatC(as.integer(nb), format = "d"), " bytes (max 249 + run suffix).\n  ", as.character(basename(file_name_string))))
        }
        ##
        return(file_name_string)
        
}
##
## ---- 3. Function to run BayesMVP models (FIXED resume)  -----------------------------------------------------------------------------------
##
##
## ---- target min ESS per model type and N (the ps3 targets used for the "time to min ESS" benchmarks); NA = no target:
##
R_fn_ps7_get_target_min_ESS <-  function( N,
                                         Model_type) {
  
        if (Model_type == "LC_MVP") {
              if (N == 500)   return(7000)
              if (N == 2500)  return(2500)
              if (N == 10000) return(1000)
              if (N == 50000) return(1000)
        } else if (Model_type == "LC_MVOP") {
              if (N == 250)   return(15000)
              if (N == 1000)  return(9000)
              if (N == 5000)  return(3000)
        }
        return(NA_real_)
  
}


##
## ---- Cache of initialise_model() objects, keyed by model type, N, MCMC seed and the model arguments: the initialisation
##      is the same for every sampler configuration at a given N and seed, so run_ps7_models() does it once per seed and
##      reuses it for every later fit. fn_ps7_clear_init_object_cache() empties the cache (e.g. after changing the data).
##      The same cache also holds the objects passed to R_fn_sample_model() (keys starting "sampler_"), built from the final
##      model_args_list so that the sampler reuses them instead of initialising the model again.
##
ps7_init_object_cache <-  new.env(parent = emptyenv())
##
fn_ps7_clear_init_object_cache <-  function() {
        rm(list = ls(envir = ps7_init_object_cache, all.names = TRUE), envir = ps7_init_object_cache)
        invisible(NULL)
}
##
## ---- Per-fit wall-clock profiling of run_ps7_models() (settings$profile_wrapper_steps; NULL/TRUE = print, FALSE = silent).
##      fn_ps7_stamp() returns the current wall-clock time in seconds; the stamps are stored in a named numeric vector
##      per fit and only printed (never saved), so they cannot change what is fitted or saved.
##      fn_ps7_format_seconds() gives "NA" for a missing / non-finite value instead of failing inside paste0().
##
fn_ps7_stamp <-  function() {
        as.numeric(Sys.time())
}
##
fn_ps7_format_seconds <-  function(x) {
        if (length(x) == 1 && is.numeric(x) && is.finite(x)) formatC(x, format = "f", digits = 2) else "NA"
}
##
## ----------------------------------------------------------------------------------------------------------------------------------
##
run_ps7_models <-  function( Model_type,
                            N,
                            y,
                            settings,
                            vect_type,
                            runs_override = NULL,
                            use_disk = use_disk,
                            use_disk_path = "/tmp/hmc_traces",
                            use_disk_path_post_hoc_dir = "/tmp/constrain_traces",
                            ##
                            true_prev,
                            prior_LKJ_nd,
                            prior_LKJ_d,
                            prior_prev_a,
                            prior_prev_b,
                            ##
                            pop,
                            ##
                            LC_MVOP_grouping,
                            ##
                            dry_run = FALSE
) {
  
        require(BayesMVP)
        if (is.null(settings$burnin_algorithm)) {
            stop("New PS7 runs require settings$burnin_algorithm. 
                 Legacy tau_objective settings remain readable through the filename/result readers.")
        }
        settings$vect_type <-  vect_type
        ##
        ## ---- tau_initial: a number or "adaptive" (see fn_ps7_tau_initial_value); the file name, the sampler call and the
        ##      requested/reported check below all read this one value:
        if (length(settings$tau_initial) == 1) settings$tau_initial <-  fn_ps7_tau_initial_value(settings$tau_initial)
        ##
        cat(paste0("\n========== PS7 Phase 1: Running models for N = ", formatC(as.integer(N), format = "d"), " ==========\n"))
        if (!is.null(settings$burnin_post_adapt_iter)) {
            stop("burnin_post_adapt_iter has been removed. Set n_burnin to the total burn-in length; old _pa results remain readable.")
        }
        ## Older callers omitted this setting because PS7 hard-coded FALSE in the sampler call.
        if (is.null(settings$store_log_lik_trace)) settings$store_log_lik_trace <-  FALSE
        ##
        ## ---- ADAM settings of the NicoStan burn-in, separately for the step-size (eps) and the trajectory-length (tau) update: set in the
        ##      runner's settings block (ps7_settings) and read here, before the resume check builds the file name. The same values go to
        ##      the NicoStan options set just before the sampler call, to the file name ("_A" code) and to the saved settings:
        ##
        adam_setting_names <-  c("eps_adam_beta1", "eps_adam_beta2", "eps_adam_epsilon", 
                                 "tau_adam_beta1", "tau_adam_beta2", "tau_adam_epsilon",
                                 "tau_learning_rate_restart_at_metric_end")
        ## a runner session started before tau_adaptation_scheme existed runs the scheme it always ran:
        if (is.null(settings$tau_adaptation_scheme)) settings$tau_adaptation_scheme <-  "adam_decay"
        if (!all(adam_setting_names %in% names(settings)) || any(vapply(settings[adam_setting_names],
                                                                        function(value) length(value) != 1 || is.na(value), FALSE))) {
            stop("ADAM settings missing (", 
                 paste(adam_setting_names, collapse = ", "), 
                 "): re-source the runner (they are set in its settings block).")
        }
        ps7_eps_adam_beta1   <-  settings$eps_adam_beta1
        ps7_eps_adam_beta2   <-  settings$eps_adam_beta2
        ps7_eps_adam_epsilon <-  settings$eps_adam_epsilon
        ps7_tau_adam_beta1   <-  settings$tau_adam_beta1
        ps7_tau_adam_beta2   <-  settings$tau_adam_beta2
        ps7_tau_adam_epsilon <-  settings$tau_adam_epsilon
        ps7_tau_learning_rate_restart_at_metric_end <-  settings$tau_learning_rate_restart_at_metric_end
        ps7_tau_adaptation_scheme <-  settings$tau_adaptation_scheme
        adam_settings_requested <-  fn_ps7_effective_adam_settings(settings = settings)
        message(BayesMVP:::colourise(paste0("ADAM: eps beta1 = ", adam_settings_requested$eps_adam_beta1, ", beta2 = ", 
                                            adam_settings_requested$eps_adam_beta2,
                                            ", denominator constant = ", adam_settings_requested$eps_adam_epsilon,
                                            " | tau beta1 = ", adam_settings_requested$tau_adam_beta1, ", beta2 = ", 
                                            adam_settings_requested$tau_adam_beta2,
                                            ", denominator constant = ", adam_settings_requested$tau_adam_epsilon,
                                            " | tau learning-rate restart at the metric freeze = ",
                                            adam_settings_requested$tau_learning_rate_restart_at_metric_end,
                                            " | tau adaptation scheme = ", adam_settings_requested$tau_adaptation_scheme,
                                            " | file-name code: ", { adam_code_for_message <-  fn_ps7_adam_file_name_code(settings = settings)
                                                                     if (nzchar(adam_code_for_message)) paste0("_A", adam_code_for_message) else "none (_ab1)" }), "cyan"))
        ## the NicoStan ADAM options are set before every sampler call and restored when run_ps7_models() returns:
        NicoStan_adam_option_names <-  c("NicoStan_eps_adam_beta1", "NicoStan_eps_adam_beta2", "NicoStan_eps_adam_epsilon",
                                         "NicoStan_tau_adam_beta1", "NicoStan_tau_adam_beta2", "NicoStan_tau_adam_epsilon",
                                         "NicoStan_tau_learning_rate_restart_at_metric_end", "NicoStan_tau_adaptation_scheme")
        previous_NicoStan_adam_options <-  stats::setNames(lapply(NicoStan_adam_option_names, getOption), NicoStan_adam_option_names)
        on.exit(expr = options(previous_NicoStan_adam_options), add = TRUE)
        ## Snapshot once for resume, fitting and saving. An explicit saved setting takes precedence over the session option.
        if (is.null(settings$autodiff_fallback)) {
            settings$autodiff_fallback <-  getOption(x = "BayesMVP_autodiff_fallback", default = FALSE)
        }
        if (!is.logical(settings$autodiff_fallback) || length(settings$autodiff_fallback) != 1 || is.na(settings$autodiff_fallback)) {
            stop("settings$autodiff_fallback / BayesMVP_autodiff_fallback must be a single TRUE or FALSE.")
        }
        settings$autodiff_fallback <-  isTRUE(settings$autodiff_fallback) &&
            (is.null(settings$multi_attempts) || isTRUE(settings$multi_attempts)) &&
            Model_type %in% c("LC_MVP", "MVP", "LC_MVOP", "MVOP", "latent_trait")
        previous_autodiff_fallback_option <-  options(BayesMVP_autodiff_fallback = settings$autodiff_fallback)
        on.exit(expr = options(previous_autodiff_fallback_option), add = TRUE)
        cat("Same-position autodiff fallback = ", settings$autodiff_fallback,
            if (settings$autodiff_fallback) " (filename: _af1)" else " (existing filenames unchanged)", "\n", sep = "")
        ##
        settings <-  fn_ps7_resolve_burnin_schedule(settings = settings)
        ##
        # n_threads_WCP_burnin   <-  settings$n_threads_WCP_burnin
        # n_threads_WCP_sampling <-  settings$n_threads_WCP_sampling
        ##
        n_tests <-  ncol(y)
        ##
        if (Model_type == "LC_MVOP") {
          N <-  nrow(y)
          grouping <-  LC_MVOP_grouping
        } else {
          N <-  N
          grouping <-  FALSE
        }
        ##
        n_class <-  2
        ##
        ## ============================================================================================
        ## HOISTED deterministic quantities (NEW): these fully determine the saved filename, so they
        ## are computed ONCE here, used by the resume check below, and REUSED inside the run loop --
        ## guaranteeing check-name == save-name.
        ## ============================================================================================
        ##
        ## ---- n_pops (as set in the priors section of the run loop):
        ##
        if (Model_type == "LC_MVOP") {
          n_pops_used <-  3
        } else {
          n_pops_used <-  length(unique(pop))
        }
        ##
        ## ---- Prior overrides (these OVERRIDE the function arguments -- as before, just hoisted):
        ##
        if (Model_type == "LC_MVOP") {
          prior_prev_a <-  1.5
          prior_prev_b <-  10
        } else {
          prior_prev_a <-  2.5
          prior_prev_b <-  10
        }
        ##
        if (Model_type == "LC_MVOP") {
          prior_LKJ_nd <-  4.0 ## -------------------------------------
          prior_LKJ_d  <-  4.0 ## -------------------------------------
        } else {
          prior_LKJ_nd <-  4.0 ## -------------------------------------
          prior_LKJ_d  <-  4.0 ## -------------------------------------
        }
        ##
        ## ---- num_chunks + WCP thread counts (single source of truth; in the saved filename):
        ##
        num_chunks_used <-  R_fn_ps7_num_chunks(N, n_tests)
        ## Chunk counts set explicitly in settings take precedence over the N-based defaults below
        ## (e.g. num_chunks_burnin = num_chunks_sampling = 25 reproduces what every run before
        ##  actually used at N = 10000, whatever its filename said):
        num_chunks_burnin_user   <-  settings$num_chunks_burnin
        num_chunks_sampling_user <-  settings$num_chunks_sampling
        ##
        print(paste("num_chunks_used = ", num_chunks_used))
        ##
        settings$num_chunks_burnin   <-  num_chunks_used
        settings$num_chunks_sampling <-  num_chunks_used
        ##
        # n_threads_WCP_burnin ∈ {6, 12, 25}
        ##
        if (parallel::detectCores() > 16) { ## Local-HPC
          settings$n_threads_WCP_burnin   <-  min(parallel::detectCores()/settings$n_chains_burnin, num_chunks_used)
          settings$n_threads_WCP_sampling <-  1 ## num_chunks_used
        } else {  ## Laptop
          settings$n_threads_WCP_burnin   <-  min(parallel::detectCores()/settings$n_chains_burnin, num_chunks_used)
          settings$n_threads_WCP_sampling <-  1 ## num_chunks_used
        }
        # ## --- PS1 burnin results:
        # HPC
        # N	      4 chains                   	8 chains                   	  16 chains
        # 500	    8 WCP,  10 chunks (1.12x)	  4 WCP,  4 chunks (1.33x)	    4 WCP,  4 chunks (1.24x)    [best = 8 chains, 4 WCP,  4 chunks]
        # 2500	  8 WCP,  10 chunks (2.36x)	  8 WCP,  25 chunks (2.18x)  	  4 WCP,  4 chunks (2.14x)    [best = 4 chains, 8 WCP,  10 chunks]
        # 10,000	16 WCP, 100 chunks (4.23x)	8 WCP,  100 chunks (3.37x)	  4 WCP,  100 chunks (2.74x)  [best = 4 chains, 16 WCP, 100 chunks]
        # 50,000	24 WCP, 500 chunks (4.42x)	12 WCP, 500 chunks (2.87x)	  11 WCP, 100 chunks (2.50x)  [best = 4 chains, 24 WCP, 500 chunks]
        # ##
        if (N == 50000) { 
            settings$n_chains_burnin <-  4 ## ---- from PS1 (burnin part)
            ##
            settings$n_threads_WCP_burnin <-  24   ## ---- from PS1 (burnin part)
            settings$num_chunks_burnin    <-  500  ## ---- from PS1 (burnin part)
            settings$num_chunks_sampling  <-  100  ## ---- from PS1 (sampling part)
        } else if (N == 10000) {
            settings$n_chains_burnin <-  4 ## ---- from PS1 (burnin part)
            ##
            settings$n_threads_WCP_burnin <-  16  ## ---- from PS1 (burnin part)
            settings$num_chunks_burnin    <-  100 ## ---- from PS1 (burnin part)
            settings$num_chunks_sampling  <-  25  ## ---- from PS1 (sampling part)
        } else if (N == 5000) { 
              settings$n_threads_WCP_burnin <-  8 ## --------------------------------------------------------------------------------------
              settings$num_chunks_burnin    <-  10 ## -------------------------------------------------------------------------------------
              settings$num_chunks_sampling  <-  10 ## -------------------------------------------------------------------------------------
        } else if (N == 2500) { 
              settings$n_chains_burnin <-  4 ## -------------------------------------------
              ##
              settings$n_threads_WCP_burnin <-  8  ## ---- from PS1 (burnin part)
              settings$num_chunks_burnin    <-  10 ## ---- from PS1 (burnin part)
              settings$num_chunks_sampling  <-  10
        } else if (N == 500) { 
              settings$n_threads_WCP_burnin <-  2 ## ---- from PS1 (burnin part)
              settings$num_chunks_burnin    <-  2 ## ---- from PS1 (burnin part)
        }
        ##
        ## ---- Burn-in chain count override (2026-09-30): environment variable PS7_N_CHAINS_BURNIN_OVERRIDE (unset = the settings above).
        ##      Used to test whether more burn-in chains change the tau adaptation; the WCP threads per chain are reduced so that the
        ##      burn-in stays within the machine's threads (180 on the local-HPC, 16 on the laptop):
        ps7_n_chains_burnin_override <-  Sys.getenv("PS7_N_CHAINS_BURNIN_OVERRIDE", unset = "")
        if (nzchar(ps7_n_chains_burnin_override)) {
              settings$n_chains_burnin <-  as.numeric(ps7_n_chains_burnin_override)
              ps7_total_burnin_threads <-  if (parallel::detectCores() > 16) 180 else 16
              settings$n_threads_WCP_burnin <-  max(1, min(floor(ps7_total_burnin_threads / settings$n_chains_burnin), settings$num_chunks_burnin))
              message(BayesMVP:::colourise(paste0("PS7_N_CHAINS_BURNIN_OVERRIDE: n_chains_burnin = ", settings$n_chains_burnin,
                                                  ", n_threads_WCP_burnin = ", settings$n_threads_WCP_burnin,
                                                  " (", settings$n_chains_burnin * settings$n_threads_WCP_burnin, " threads)"), "cyan"))
        }
        ##
        if (!is.null(num_chunks_burnin_user))   settings$num_chunks_burnin   <-  num_chunks_burnin_user
        if (!is.null(num_chunks_sampling_user)) settings$num_chunks_sampling <-  num_chunks_sampling_user
        ## each within-chain thread needs at least one chunk, so a user-set chunk count caps the thread counts:
        if (!is.null(num_chunks_burnin_user) && (settings$n_threads_WCP_burnin > settings$num_chunks_burnin)) {
          message("num_chunks_burnin = ", settings$num_chunks_burnin, " < n_threads_WCP_burnin = ", settings$n_threads_WCP_burnin,
                  ": lowering n_threads_WCP_burnin to ", settings$num_chunks_burnin)
          settings$n_threads_WCP_burnin <-  settings$num_chunks_burnin
        }
        if (!is.null(num_chunks_sampling_user) && (settings$n_threads_WCP_sampling > settings$num_chunks_sampling)) {
          settings$n_threads_WCP_sampling <-  settings$num_chunks_sampling
        }
        ##
        message(" ------------------------------------------------------------------------------ ")
        print(paste0("settings$num_chunks_burnin = ", settings$num_chunks_burnin))
        print(paste0("settings$num_chunks_sampling = ", settings$num_chunks_sampling))
        print(paste0("settings$n_threads_WCP_burnin = ", settings$n_threads_WCP_burnin))
        print(paste0("settings$n_threads_WCP_sampling = ", settings$n_threads_WCP_sampling))
        ##
        print(paste("total # threads for burnin = ", settings$n_threads_WCP_burnin * settings$n_chains_burnin  ))
        ##
        ## ============================================================================================
        ##
        ## ---- Storage for all runs:
        ##
        all_runs_results <-  list()
        ##
        run_i <-  1 ## debug
        ##
        ## ---- Check for existing runs:
        ##
        ## FIX: the old inline filename here was the RETIRED format and NEVER matched what
        ## R_fn_file_name_string() saves, so resume never worked. Now the check uses
        ## R_fn_file_name_string() itself, fed with the hoisted values above.
        ##
        existing_runs <-  c()
        ##
        ## ---- Wall-clock profiling switch (NULL/TRUE = print one timing line per fit, FALSE = silent) and the start of the
        ##      file-name / resume step (done once per run_ps7_models() call, so it is reported separately from the per-fit steps):
        ##
        profile_wrapper_steps <-  is.null(settings$profile_wrapper_steps) || isTRUE(settings$profile_wrapper_steps)
        ps7_resume_time_start <-  fn_ps7_stamp()
        ##
        print(paste("clip_iter = ", settings$clip_iter))
        print(paste("int = ", settings$int))
        ##
        ## ---- time-criterion settings of CHESSR_time / SNAPER_time ("_bact" / "_bast"): keep only the EFFECTIVE values in the
        ##      settings (see fn_ps7_effective_time_criterion_settings), so the file name, the sampler arguments and the saved settings
        ##      agree, and every time-criterion field is absent from the settings of every other algorithm (whose saved settings
        ##      are therefore unchanged):
        ##
        {
            effective_time_criterion_settings <-  fn_ps7_effective_time_criterion_settings(settings = settings, N = N, Model_type = Model_type)
            for (time_criterion_setting_name in ps7_time_criterion_setting_names) {
                settings[[time_criterion_setting_name]] <-  NULL   ## removes the element (no named NULL is left behind)
            }
            if (!is.null(effective_time_criterion_settings)) {
                for (time_criterion_setting_name in names(effective_time_criterion_settings)) {
                    if (!is.null(effective_time_criterion_settings[[time_criterion_setting_name]])) {
                        settings[[time_criterion_setting_name]] <-  effective_time_criterion_settings[[time_criterion_setting_name]]
                    }
                }
                if (!is.null(settings$time_criterion_previous_run_path) && !all(file.exists(settings$time_criterion_previous_run_path))) {
                    stop(paste0("time_criterion_previous_run_path: no saved run at ",
                                paste(settings$time_criterion_previous_run_path[!file.exists(settings$time_criterion_previous_run_path)],
                                      collapse = ", ")))
                }
                ##
                ## ---- the previous runs must be of this model, N and sampling configuration: their sampling times are used as this run's.
                ##      (n_iter need not match: every quantity taken from them is per sampling iteration. A setting that is NULL here, e.g.
                ##      num_chunks_sampling left to the sampler, is not compared.)
                ##
                if (!is.null(settings$time_criterion_previous_run_path)) {
                    previous_run_names_parsed <-  dplyr::bind_rows(lapply(X = basename(settings$time_criterion_previous_run_path), 
                                                                          FUN = fn_ps7_parse_run_name))
                    current_run_configuration <-  list(model_type             = Model_type,
                                                       N                      = N,
                                                       n_chains_sampling      = settings$n_chains_sampling,
                                                       n_threads_WCP_sampling = settings$n_threads_WCP_sampling,
                                                       num_chunks_sampling    = settings$num_chunks_sampling)
                    for (configuration_name in names(current_run_configuration)) {
                        if (is.null(current_run_configuration[[configuration_name]])) next
                        previous_run_values <-  previous_run_names_parsed[[configuration_name]]
                        current_value <-  current_run_configuration[[configuration_name]]
                        values_differ <-  if (is.numeric(previous_run_values)) {
                            is.na(previous_run_values) | previous_run_values != as.numeric(current_value)
                        } else {
                            is.na(previous_run_values) | previous_run_values != as.character(current_value)
                        }
                        if (any(values_differ)) {
                            stop(paste0("time_criterion_previous_run_path: ", configuration_name, " of this run is ", current_value, ", but ",
                                        paste(basename(settings$time_criterion_previous_run_path)[values_differ], collapse = ", "),
                                        " has ", paste(previous_run_values[values_differ], collapse = ", "),
                                        ". Use saved runs of the same model, N and sampling configuration."))
                        }
                    }
                }
            }
        }
        ##
        file_name_string_check <-  R_fn_file_name_string( Model_type = Model_type,
                                                         N = N,
                                                         settings = settings,
                                                         model_args_list = list( n_pops = n_pops_used,
                                                                                 num_chunks = settings$num_chunks_burnin),
                                                         prior_LKJ_nd = prior_LKJ_nd,
                                                         prior_LKJ_d = prior_LKJ_d,
                                                         prior_prev_a = prior_prev_a,
                                                         prior_prev_b = prior_prev_b,
                                                         grouping = grouping)
        ##
        ## ---- Dry run: return the expected run files and whether each already exists, BEFORE any readRDS() or fit
        ##      (same names as the resume loop below, one per run_i in 1:settings$n_runs):
        ##
        if (isTRUE(dry_run)) {
            dry_run_files <-  paste0(file_name_string_check, "_run", 1:settings$n_runs)
            return(invisible(list(dry_run = TRUE,
                                  N = N,
                                  n_burnin = settings$n_burnin,
                                  run_files = dry_run_files,
                                  existing = file.exists(dry_run_files))))
        }
        ##
        for (run_i in 1:settings$n_runs) {
          
              run_file <-  paste0(file_name_string_check, "_run", run_i)
              ##
              if (file.exists(run_file)) {
                existing_runs <-  c(existing_runs, run_i)
                all_runs_results[[run_i]] <-  readRDS(run_file)
                ## CHESSR_time / SNAPER_time: the 4-digit settings hash of "_bac<4 hex>" / "_bas<4 hex>" could in principle be shared
                ## by two different settings texts, so a saved run is resumed only when its own text is the current one:
                if (!is.null(effective_time_criterion_settings) &&
                    !identical(all_runs_results[[run_i]]$time_criterion_inputs$time_criterion_settings_text,
                               fn_ps7_time_criterion_settings_text(effective_time_criterion_settings = effective_time_criterion_settings))) {
                    stop(paste0("The saved run ", basename(run_file), " was fitted with the time-criterion settings '",
                                paste(all_runs_results[[run_i]]$time_criterion_inputs$time_criterion_settings_text, collapse = ""),
                                "', not the current '", paste(fn_ps7_time_criterion_settings_text(effective_time_criterion_settings), collapse = ""),
                                "' (the same 4-digit hash); move that run elsewhere before fitting these settings."))
                }
                if (isTRUE(settings$debug_burnin_timing) && is.null(all_runs_results[[run_i]]$sampler_diagnostics$burnin_profile)) {
                    warning("Loaded existing run without burn-in timing diagnostics: ", run_file,
                            ". Use a separate output directory for new profiled runs.", call. = FALSE)
                }
                cat(paste0("  Loaded existing run ", formatC(as.integer(run_i), format = "d"), "\n"))
              }
          
        }
        ##
        runs_to_do <-  setdiff(1:settings$n_runs, existing_runs)
        ##
        ps7_resume_time_seconds <-  fn_ps7_stamp() - ps7_resume_time_start  ## file name + file.exists() + readRDS() of existing runs
        ##
        if (length(runs_to_do) == 0) {
          cat(paste0("All ", formatC(as.integer(settings$n_runs), format = "d"), " runs complete for N=", formatC(as.integer(N), format = "d"), ", n_burnin=", formatC(as.integer(settings$n_burnin), format = "d"), "\n"))
        } else {
          cat(paste0("Running: ", as.character(paste(runs_to_do, collapse = ", ")), "\n"))
        }
        ##
        ## ---- Run models:;
        ##
        if (!(is.null(runs_override))) {
          runs_to_do <-  runs_override
        }
        ##
        run_i <-  1
        ##
        for (run_i in runs_to_do) {
          
          ps7_step_times <-  c(fit_start = fn_ps7_stamp())  ## wall-clock stamps of this fit's steps (printed after the save)
          ##
          MCMC_seed <-  run_i * 1000  # DIFFERENT MCMC seed each run
          ##
          set.seed(MCMC_seed, kind = "L'Ecuyer-CMRG")
          ##
          cat(paste0("\n--- Run ", formatC(as.integer(run_i), format = "d"), "/", formatC(as.integer(settings$n_runs), format = "d"), " (MCMC seed = ", formatC(as.integer(MCMC_seed), format = "d"), ") ---\n"))
          ##
          ## ---- Initialize model:
          ##
          model_args_list <-  list(y = y,
                                  n_tests = n_tests,
                                  n_class = 2,
                                  N = N)
          ##
          ## Use the installed provider and shared NicoStan implementation selected by the driver.
          ## No per-run source()/sourceCpp() overrides may bypass that package build.
          ##
          n_coeffs_MVP  <-  1 ## intercept-only
          n_coeffs_MVOP <-  3
          ##
          if (Model_type == "LC_MVOP") {
            
                model_args_list$X <-  list()
                ##
                for (c in 1:2) {
                  model_args_list$X[[c]] <-  list()
                  for (t in 1:3) {
                    model_args_list$X[[c]][[t]] <-  matrix(data = 1, nrow = N, ncol = 2)
                    model_args_list$X[[c]][[t]][1:N, 2] <-  pop
                  }
                }
                ## Extract the language column (col 2) from any X matrix:
                lang_vec <-  model_args_list$X[[1]][[1]][, 2]
                
                ## Expand to dummies:
                lang_dummies <-  expand_categorical_to_dummies( lang_vec,
                                                               ref_level = 1L, ## ----  Note: 1 = Afrikaans, 2 = Xhosa, 3 = Zulu
                                                               prefix = "lang")
                ## Result: lang_2 and lang_3 columns (lang_1 = reference)
                
                ## Replace col 2 with the two dummies, keeping intercept:
                for (c in 1:n_class) {
                  for (t in 1:n_tests) {
                    X_old <-  model_args_list$X[[c]][[t]]
                    if (n_coeffs_MVOP == 3) model_args_list$X[[c]][[t]] <-  cbind(X_old[, 1, drop = FALSE], lang_dummies)
                    else                    model_args_list$X[[c]][[t]] <-  lang_dummies ## cbind(X_old[, 1, drop = FALSE], lang_dummies)
                  }
                }
            
          } else {
            
                ## (binary: default X from initialise_model -- see original file for the
                ##  commented-out covariate-X experiment that used to live here)
            
          }
          ##
          if (Model_type == "LC_MVOP") {
            if (LC_MVOP_grouping == TRUE) { 
                 model_args_list$n_cat_per_ord_test <-  c(10, 10)
            } else { 
                 model_args_list$n_cat_per_ord_test <-  c(28, 31)
            }
          }
          ##
          ## ---- Initialise once per (model, N, seed, model arguments) and reuse: the initialisation does not depend on the
          ##      sampler configuration, so with n_runs seeds it is done n_runs times per N rather than once per fit.
          ##      initialise_model() returns a plain list, so the cached object cannot be modified by a fit.
          ##
          init_object_cache_key <-  paste0(Model_type, "_N", N, "_seed", MCMC_seed, "_", digest::digest(model_args_list))
          ##
          if (exists(x = init_object_cache_key, envir = ps7_init_object_cache, inherits = FALSE)) {
                init_object <-  get(x = init_object_cache_key, envir = ps7_init_object_cache, inherits = FALSE)
                cat(BayesMVP:::colourise(paste0("  Reusing the cached initialisation for seed ", MCMC_seed), "cyan"), "\n")
          } else {
                init_object <-  BayesMVP::initialise_model( Model_type = Model_type,
                                                  stream = MCMC_seed,
                                                  sample_nuisance = TRUE,
                                                  n_nuisance_override = NULL,
                                                  model_args_list = model_args_list,
                                                  compile = TRUE,
                                                  force_recompile = FALSE,
                                                  cmdstanr_model_fit_obj = NULL,
                                                  Stan_data_list = NULL,
                                                  Stan_model_file_path = NULL,
                                                  Stan_cpp_user_header = NULL,
                                                  Stan_cpp_flags = NULL,
                                                  stanc_args = NULL)
                assign(x = init_object_cache_key, value = init_object, envir = ps7_init_object_cache)
          }
          ##
          ps7_step_times["init"] <-  fn_ps7_stamp()  ## model_args_list + cache key digest + cache hit or initialise_model()
          ##
          model_args_list <-  init_object$model_args_list
          ##
          ## ---- Set priors (values HOISTED above -- see top of function):
          ##
          model_args_list$n_pops <-  n_pops_used
          model_args_list$pop <-  pop
          ##
          model_args_list$prior_prev_a <-  rep(prior_prev_a, model_args_list$n_pops)
          model_args_list$prior_prev_b <-  rep(prior_prev_b, model_args_list$n_pops)
          ##
          model_args_list$prior_prev_a
          model_args_list$prior_prev_b
          ##
          if (Model_type == "LC_MVP") {
            model_args_list$prior_coeffs_mean_mat <-  rep(list(matrix(0.0, nrow = n_coeffs_MVP, ncol = n_tests)), 2)
            model_args_list$prior_coeffs_sd_mat   <-  rep(list(matrix(1.0, nrow = n_coeffs_MVP, ncol = n_tests)), 2)
          } else {
            model_args_list$prior_coeffs_mean_mat <-  rep(list(matrix(0.0, nrow = n_coeffs_MVOP, ncol = n_tests)), 2)
            model_args_list$prior_coeffs_sd_mat   <-  rep(list(matrix(1.0, nrow = n_coeffs_MVOP, ncol = n_tests)), 2)
          }
          ##
          if (Model_type == "LC_MVOP") {
            
                for (c in 1:2) {
                  model_args_list$prior_coeffs_mean_mat[[c]][,] <-  0.0
                  model_args_list$prior_coeffs_sd_mat[[c]][,]   <-  0.50
                }
                ##
                model_args_list$prior_coeffs_mean_mat[[1]][1, 1]  <-  -1.5 ## ---- for ref test
                model_args_list$prior_coeffs_mean_mat[[2]][1, 1]  <-  +1.0 ## +0.5 ## D+ ## ---- for ref test
                ##
                model_args_list$prior_coeffs_sd_mat[[1]][1, 1] <-  0.25 ## 1.0 ## 0.50 ## 0.25 ## D- ## ---- for ref test
                model_args_list$prior_coeffs_sd_mat[[2]][1, 1] <-  0.25 ## 1.0 ## 0.50 ## 0.25 ## D+ ## ---- for ref test
                ##
                #  model_args_list$prior_coeffs_sd_mat[[1]][1, 2:3]   <-  1.0 ## 2.0
                #  model_args_list$prior_coeffs_sd_mat[[2]][1, 2:3]   <-  1.0 ## 2.0
                ##
                model_args_list$prior_coeffs_sd_mat[[1]][2:3, 1] <-  0.50 ## ---- for ref test
                model_args_list$prior_coeffs_sd_mat[[2]][2:3, 1] <-  0.50 ## ---- for ref test
                ##
                # model_args_list$prior_coeffs_mean_mat[[1]][1, 2:3] <-  -1.25 ## ----------------
                # model_args_list$prior_coeffs_mean_mat[[2]][1, 2:3] <-  +1.25 ## ----------------
                # model_args_list$prior_coeffs_sd_mat[[1]][1, 2:3]   <-  0.5 ## ----------------
                # model_args_list$prior_coeffs_sd_mat[[2]][1, 2:3]   <-  0.5 ## ----------------
                ##
                model_args_list$prior_coeffs_sd_mat[[1]][1, 2:3]   <-  1.0 ## ----------------
                model_args_list$prior_coeffs_sd_mat[[2]][1, 2:3]   <-  1.0 ## ----------------
                ##
                # model_args_list$prior_coeffs_sd_mat[[1]][1, 1] <-  0.01 ## 1.0 ## 0.50 ## 0.25
                # model_args_list$prior_coeffs_sd_mat[[2]][1, 1] <-  0.01 ## 1.0 ## 0.50 ## 0.25 ## D+
                ##
                model_args_list$prior_coeffs_mean_mat
                model_args_list$prior_coeffs_sd_mat
            
          }  else {
                
                ## All tests: Sp roughly 97-99%
                ## beta_nd ~ N(-2.0, 0.5) -> median Sp ~ 97.7%, 95% CI ~ (84%, 99.9%)
                model_args_list$prior_coeffs_mean_mat[[1]][1, ] <-  -2.0
                model_args_list$prior_coeffs_sd_mat[[1]][1, ]   <-  0.5
                ##
                ## Lab tests (euroimmun, roche): Se roughly 90-98%
                ## beta_d ~ N(1.5, 0.5) -> median Se ~ 93%, 95% CI ~ (69%, 99.4%)
                model_args_list$prior_coeffs_mean_mat[[2]][1, 1:2] <-  1.5  ## euroimmun, roche
                model_args_list$prior_coeffs_sd_mat[[2]][1, 1:2]   <-  0.5
                ##
                ## Lateral flow (abc19, surescreen, orientgene, biomerica): Se roughly 70-90%
                ## beta_d ~ N(1.0, 0.5) -> median Se ~ 84%, 95% CI ~ (50%, 97.7%)
                model_args_list$prior_coeffs_mean_mat[[2]][1, 3:6] <-  1.0
                model_args_list$prior_coeffs_sd_mat[[2]][1, 3:6]   <-  0.5
            
          }
          ##
          model_args_list$prior_coeffs_mean_mat
          model_args_list$prior_coeffs_sd_mat
          ##
          ## ---- Correlation priors (LKJ values HOISTED above -- see top of function):
          ##
          corr_prior_sets <-  get_correlation_prior_sets()
          corr_prior_settings <-  corr_prior_sets$LC_MVP[[paste0("DGM_", settings$DGM)]][[settings$prior_name]]
          ##
          model_args_list$prior_LKJ <-  matrix(c(prior_LKJ_nd, prior_LKJ_d), ncol = 1) ## --------------------------
          ##
          model_args_list$lkj_cholesky_eta <-  model_args_list$prior_LKJ
          model_args_list$corr_force_positive <-  corr_prior_settings$corr_force_positive
          ##
          if (Model_type == "LC_MVOP") {
            # model_args_list$corr_force_positive <-  TRUE
            model_args_list$corr_force_positive <-  FALSE
          } else {
            model_args_list$corr_force_positive <-  FALSE
          }
          # model_args_list$corr_force_positive <-  TRUE ## --------------------------
          ##
          if (isTRUE(model_args_list$corr_force_positive)) {
            for (c in 1:n_class) model_args_list$lb_corr[[c]] <-  matrix(0.0, ncol = n_tests, nrow = n_tests)
          } else {
            for (c in 1:n_class) model_args_list$lb_corr[[c]] <-  matrix(-1.0, ncol = n_tests, nrow = n_tests)
          }
          for (c in 1:n_class) model_args_list$ub_corr[[c]] <-  matrix(+1.0, ncol = n_tests, nrow = n_tests)
          ##
          if (!is.null(corr_prior_settings$lb_corr)) {
            if (Model_type == "LC_MVOP") {
              model_args_list$lb_corr <-  matrix(-1.0, ncol = n_tests, nrow = n_tests)
              model_args_list$ub_corr <-  matrix(+1.0, ncol = n_tests, nrow = n_tests)
            } else {
              model_args_list$lb_corr <-  corr_prior_settings$lb_corr
              model_args_list$ub_corr <-  corr_prior_settings$ub_corr
            }
          }
          ##
          model_args_list$n_covariates_per_outcome_mat
          ##
          ## ---- Inits:
          ##
          init_lists_per_chain <-  list()
          ##
          for (kk in 1:settings$n_chains_burnin) {
                
                prev_init_LC_MVP  <-  0.10
                ##
                # prev_init_LC_MVOP <-  0.05
                prev_init_LC_MVOP <-  c(0.20, 0.10, 0.025)   ## per-pop (Afrikaans, Xhosa, Zulu)
                ##
                # random_draws <-  rnorm(n = n_tests*N, mean = 0.0, sd = 0.10)
                random_draws <-  rnorm(n = n_tests*N, mean = 0.0, sd = 0.01)
                ##
                n_corrs <-  n_tests * (n_tests - 1) / 2
                ##
                n_thr_total <-  sum(model_args_list$n_thr_per_ord_test)
                ##
                if (Model_type == "LC_MVOP") {
                  
                      # ##
                      # ## ---- Cutpoint inits: spread each test's cutpoints across the probability
                      # ##      scale (equal-mass categories) instead of packing them all into a
                      # ##      narrow band. Offset by the class intercept so the cutpoints sit
                      # ##      where that class's latent scale actually is.
                      # ##
                      # ##      NB: C_raw_vec is a RAGGED FLAT vector of length n_thr_total, blocked
                      # ##      by ordinal test in slot order (test 1's thresholds, then test 2's, ...)
                      # ##      -- must match calculate_start_indices() in the C++/Stan side.
                      # ##
                      # C_raw_init_list <-  list()
                      # ##
                      # for (c in 1:n_class) {
                      #   
                      #       beta_c_intercept <-  init_lists_per_chain_beta_intercept <-  ifelse(c == 1, -1.0, +1.0)
                      #       ##
                      #       blocks <-  list()
                      #       ##
                      #       for (t_ord in seq_along(model_args_list$n_thr_per_ord_test)) {
                      #         
                      #             n_thr_t <-  model_args_list$n_thr_per_ord_test[t_ord]
                      #             K_t     <-  n_thr_t + 1
                      #             ##
                      #             ## ---- Equal-mass cutpoints on the latent (probit) scale, shifted to
                      #             ##      sit around this class's intercept:
                      #             ##
                      #             C_t <-  qnorm(seq_len(n_thr_t) / K_t) + beta_c_intercept
                      #             ##
                      #             ## ---- Invert the log-difference parameterisation:
                      #             ##        C[1] = C_raw[1];  C[k] = C[k-1] + exp(C_raw[k])
                      #             ##
                      #             if (n_thr_t == 1) {
                      #               blocks[[t_ord]] <-  C_t[1]
                      #             } else {
                      #               blocks[[t_ord]] <-  c(C_t[1], log(diff(C_t)))
                      #             }
                      #         
                      #       }
                      #       ##
                      #       C_raw_init_list[[c]] <-  unlist(blocks)
                      #   
                      # }
                      # ##
                      # stopifnot(length(C_raw_init_list[[1]]) == n_thr_total)
                      ##
                      ## ---- Invert the bounding transform: theta = atanh(2*(C_raw - L)/(U - L) - 1)
                      ##
                      # model_args_list$C_raw_lower <-  -7.5
                      # model_args_list$C_raw_upper <-  +2.5
                      ##
                      # model_args_list$C_raw_lower <-  -7.5
                      # model_args_list$C_raw_upper <-  +1.25
                      ##
                      model_args_list$C_raw_lower <-  -10.0
                      model_args_list$C_raw_upper <-  +2.5
                      ##
                      # model_args_list$C_raw_lower <-  -15.0
                      # model_args_list$C_raw_upper <-  +2.5
                      ##
                      to_unc <-  function(x) {
                            z <-  2 * (x - model_args_list$C_raw_lower) / (model_args_list$C_raw_upper - model_args_list$C_raw_lower) - 1
                            atanh(pmin(pmax(z, -0.999), 0.999))   ## guard: atanh(+-1) = Inf
                      }
                      ##
                      # C_raw_init_list[[c]] <-  to_unc(unlist(blocks))
                      to_unc(-2.0)
                            
                      init_lists_per_chain[[kk]] <-  list(
                        u_raw = matrix(random_draws, ncol = n_tests, nrow = N),
                        Omega_unconstrained_vec = list(rep(0.001, n_corrs), rep(0.001, n_corrs)),
                        beta = list(matrix(-1.0, nrow = n_coeffs_MVOP, ncol = n_tests),
                                    matrix(+1.0, nrow = n_coeffs_MVOP, ncol = n_tests)),
                        # p_raw = array(atanh(2 * prev_init_LC_MVOP - 1)),
                        p_raw = array((atanh(2 * prev_init_LC_MVOP - 1))),
                        C_unc_vec = rep(list(rep(to_unc(-2.0), n_thr_total)), n_class)
                        # C_raw_vec = C_raw_init_list
                      )
                      
                      # K <-  model_args_list$n_cat_per_ord_test[t]
                      # C_init <-  qnorm(seq_len(K - 1) / K)
                      # C_raw_init <-  c(C_init[1], log(diff(C_init)))
                      
                      for (c in 1:2) {
                        init_lists_per_chain[[kk]]$beta[[c]][c(2, 3), ] <-  0.0
                      }
                  
                } else if (Model_type == "LC_MVP") {
                  
                      init_lists_per_chain[[kk]] <-  list(
                        u_raw = matrix(random_draws, ncol = n_tests, nrow = N),
                        Omega_unconstrained_vec = list(rep(0.001, n_corrs), rep(0.001, n_corrs)),
                        beta = list(matrix(-1.0, nrow = n_coeffs_MVP, ncol = n_tests),
                                    matrix(+1.0, nrow = n_coeffs_MVP, ncol = n_tests)),
                        # p_raw = array(atanh(2 * prev_init_LC_MVP - 1))
                        p_raw = array(rep(atanh(2 * prev_init_LC_MVP - 1)), model_args_list$n_pops)
                      )
                      
                      try({
                        for (c in 1:2) {
                          init_lists_per_chain[[kk]]$beta[[1]][1, ] <-  -2.0
                          ##
                          init_lists_per_chain[[kk]]$beta[[2]][1, 1:2] <-  +1.5
                          init_lists_per_chain[[kk]]$beta[[2]][1, 3:6] <-  +1.0
                        }
                      })
                  
                }
            
          }
          init_lists_per_chain <-  BayesMVP:::resize_init_list( init_lists_per_chain = init_lists_per_chain,
                                                                  n_chains_new = settings$n_chains_burnin)
          ##
          # model_args_list$overflow_threshold  <-  +5.0
          # model_args_list$underflow_threshold <-  -5.0
          ##
          model_args_list$overflow_threshold  <-  +7.5
          model_args_list$underflow_threshold <-  -7.5
          ##
          # model_args_list$overflow_threshold  <-  +100
          # model_args_list$underflow_threshold <-  -100
          ##
          # model_args_list$overflow_threshold  <-  +0.001
          # model_args_list$underflow_threshold <-  -0.001
          ##
          ## ---- num_chunks / J_grad_option / WCP threads (HOISTED -- see R_fn_ps7_num_chunks()):
          ##
          model_args_list$num_chunks <-  settings$num_chunks_burnin
          ##
          ##
          ## ---- Jacobian gradient of the correlation (LDL / Pinkney) transform: "num_diff" = forward finite differences
          ##      (the default, and Earlier runs), "autodiff" = Stan reverse mode (slower; "_JgA" in the
          ##      file name). Any other string leaves that Jacobian at ZERO in the C++ (silently wrong gradients), so it
          ##      is checked here and read back from the model after the run:
          ##
          J_grad_option_requested <-  if (is.null(settings$J_grad_option)) "num_diff" else settings$J_grad_option
          if (!(J_grad_option_requested %in% c("num_diff",
                                               "autodiff"))) {
            stop("settings$J_grad_option must be 'num_diff' or 'autodiff'; got: ", J_grad_option_requested)
          }
          model_args_list$J_grad_option <-  J_grad_option_requested
          print(paste("J_grad_option = ", model_args_list$J_grad_option))
          ##
          ## (settings$n_threads_WCP_burnin / _sampling already set at top of function,
          ##  BEFORE the resume check; use each phase's resolved chunk count below.)
          ##
          print(paste("num_chunks = ", model_args_list$num_chunks))
          ##
          if (settings$num_chunks_burnin < settings$n_threads_WCP_burnin) {
            stop(" num_chunks_burnin must be >= n_threads_WCP_burnin!")
          }
          if (settings$num_chunks_sampling < settings$n_threads_WCP_sampling) {
            stop(" num_chunks_sampling must be >= n_threads_WCP_sampling!")
          }
          ##
          model_args_list$nuisance_transformation
          ##
          clip_iter_tau <-  settings$clip_iter + settings$int ## + round(clip_iter/2)
          ##
          int_width <-  settings$int_width
          ##
          if (!is.null(settings$n_adapt)) {
            n_adapt <-  settings$n_adapt
          } else if (settings$n_burnin == 1000) {
            n_adapt <-  900
          } else {
            n_adapt <-  settings$n_burnin - round(settings$n_burnin/10)
          }
          ##
          if (settings$partitioned_HMC == TRUE) {
            manual_tau <-  TRUE
            tau_if_manual <-  c(3.0, 1.50)
          } else {
            # manual_tau <-  TRUE ; tau_if_manual <-  c(3.0, 3.0)
            # manual_tau <-  TRUE ; tau_if_manual <-  c(1.50, 1.50)
            ##
            manual_tau <-  FALSE
            tau_if_manual <-  c(3.0, 3.0)
          }
          ##
          ## ---- FIXED-L SWEEP.
          ##
          ## settings$manual_L is a leapfrog-step count. When set, the adaptation is switched off and
          ## tau is pinned to L * eps at every burnin iteration, with eps still dual-averaging to
          ## adapt_delta. That is the experiment the tau adaptation cannot run on itself: does the
          ## adapted L (about 9 for kick_flow_kick, about 23 for flow_kick_flow) actually sit near
          ## the ESS-per-gradient optimum for this model?
          ##
          ## NULL / NA leaves every existing run untouched, adaptation and all.
          tau_if_manual_in_L_units <-  FALSE
          ##
          if (!is.null(settings$manual_L) && length(settings$manual_L) == 1 && !is.na(settings$manual_L)) {
                if (!is.numeric(settings$manual_L) || !is.finite(settings$manual_L) || settings$manual_L < 1) {
                    stop("settings$manual_L must be a single leapfrog-step count >= 1; got: ", settings$manual_L)
                }
                manual_tau <-  TRUE
                tau_if_manual <-  c(settings$manual_L, settings$manual_L)
                tau_if_manual_in_L_units <-  TRUE
          }
          ##
          ## Absolute tau stays fixed as epsilon adapts; do not convert it to L units.
          invisible(R_fn_enc_manual_tau(settings$manual_tau_value))
          if (!is.null(settings$manual_tau_value) && !is.na(settings$manual_tau_value)) {
                if (length(settings$manual_tau_value) != 1 || !is.numeric(settings$manual_tau_value) ||
                    !is.finite(settings$manual_tau_value) || settings$manual_tau_value <= 0) {
                    stop("settings$manual_tau_value must be a single positive finite integration time.")
                }
                if (isTRUE(tau_if_manual_in_L_units)) stop("manual_L and manual_tau_value cannot both be fixed.")
                manual_tau <-  TRUE
                tau_if_manual <-  rep(x = settings$manual_tau_value, times = 2)
                tau_if_manual_in_L_units <-  FALSE
          }
          ##
          ## ---- Trajectory-length objective and acceptance weighting. NULL -> the settings every
          ##      ps7 run predating these fields used, so old results stay reproducible.
          burnin_algorithm_used <-  fn_ps7_normalise_burnin_algorithm(
              burnin_algorithm = if (is.null(settings$burnin_algorithm)) "KE" else settings$burnin_algorithm)
          tau_weight_by_p_jump_used <-  isTRUE(settings$tau_weight_by_p_jump)
          ##
          ## if (!burnin_algorithm_used %in% c("KE", "ChEES", "CHESSR", "CHESSR_log", "SNAPER")) {
              ## stop("settings$burnin_algorithm must be 'KE', 'ChEES', 'CHESSR', 'CHESSR_log' or 'SNAPER'; got: ", burnin_algorithm_used)
          ## if (!burnin_algorithm_used %in% c("KE", "ChEES", "CHESSR", "CHESSR_log", "SNAPER", "CHESSR_time", "SNAPER_time")) {
              ## stop("settings$burnin_algorithm must be 'KE', 'ChEES', 'CHESSR', 'CHESSR_log', 'SNAPER', 'CHESSR_time' or 'SNAPER_time'; got: ", burnin_algorithm_used)
          ## if (!burnin_algorithm_used %in% c("KE", "ChEES", "CHESSR", "CHESSR_log", "SNAPER", "CHESSR_time", "SNAPER_time", "ESJD", "ESJD_CHESSR", "ESJD_SNAPER")) {
          if (!burnin_algorithm_used %in% c("KE", "ChEES", "CHESSR", "CHESSR_log", "SNAPER", "CHESSR_time", "SNAPER_time", "ESJD", "ESJD_CHESSR",
                                            "ESJD_SNAPER", "LQ_ESSR")) {
              stop("settings$burnin_algorithm must be 'KE', 'ChEES', 'CHESSR', 'CHESSR_log', 'SNAPER', 'CHESSR_time', 'SNAPER_time', 'ESJD', 'ESJD_CHESSR', ",
                   "'ESJD_SNAPER' or 'LQ_ESSR'; got: ", burnin_algorithm_used)
          }
          #
          # model_args_list$vect_type <-  "Stan"
          ##
          if (Model_type == "LC_MVOP") {
            
              ##
              # model_args_list$prior_dirichlet_alpha
              # ##
              # dirichlet_priors_PHQ_9 <-  MetaOrdDTA:::generate_dual_ordinal_dirichlet_priors( n_cat = 28,
              #                                                                                nd_peak = 0.15,
              #                                                                                d_peak = 0.40,
              #                                                                                nd_max_alpha = 3,
              #                                                                                d_max_alpha  = 3,
              #                                                                                nd_min_alpha = 1,
              #                                                                                d_min_alpha  = 1,
              #                                                                                smoothness = 0.10)
              # ##
              # dirichlet_priors_CES_D <-  MetaOrdDTA:::generate_dual_ordinal_dirichlet_priors( n_cat = 31,
              #                                                                                nd_peak = 0.15,
              #                                                                                d_peak = 0.50,
              #                                                                                nd_max_alpha = 3,
              #                                                                                d_max_alpha  = 3,
              #                                                                                nd_min_alpha = 1,
              #                                                                                d_min_alpha  = 1,
              #                                                                                smoothness = 0.10)
              # dirichlet_priors_PHQ_9$diseased
              # dirichlet_priors_PHQ_9$diseased[c(11, 12, 13)] ## peaks at categories 11-13
              # ##
              # dirichlet_priors_PHQ_9$non_diseased
              # dirichlet_priors_PHQ_9$non_diseased[c(2:6)]
              # ##
              # dirichlet_priors_CES_D$diseased[c(15:17)] ## peaks at categories 15-17
              # ##
              # dirichlet_priors_CES_D$non_diseased[c(3:7)]
              # ##
              # model_args_list$prior_dirichlet_alpha[[1]][1:28, 1] <-  dirichlet_priors_PHQ_9$non_diseased ## PHQ-9, D-
              # model_args_list$prior_dirichlet_alpha[[2]][1:28, 1] <-  dirichlet_priors_PHQ_9$diseased     ## PHQ-9, D+
              # ##
              # model_args_list$prior_dirichlet_alpha[[1]][1:31, 2] <-  dirichlet_priors_CES_D$non_diseased ## CES-D, D-
              # model_args_list$prior_dirichlet_alpha[[2]][1:31, 2] <-  dirichlet_priors_CES_D$diseased     ## CES-D, D+
              # ##
              # ##
              # dirichlet_priors_PHQ_9_neutral <-  generate_neutral_ordinal_dirichlet_prior(   n_cat = 28,
              #                                                                               peak = 0.25,
              #                                                                               max_alpha = 5, 
              #                                                                               smoothness = 0.25)$alpha
              # ##
              # dirichlet_priors_CES_D_neutral <-  generate_neutral_ordinal_dirichlet_prior(  n_cat = 31,
              #                                                                              peak = 0.25,
              #                                                                              max_alpha = 5, 
              #                                                                              smoothness = 0.25)$alpha
              # ##
              # model_args_list$prior_dirichlet_alpha[[1]][1:28, 1] <-  dirichlet_priors_PHQ_9_neutral
              # model_args_list$prior_dirichlet_alpha[[2]][1:28, 1] <-  dirichlet_priors_PHQ_9_neutral
              # model_args_list$prior_dirichlet_alpha[[1]][1:31, 2] <-  dirichlet_priors_CES_D_neutral
              # model_args_list$prior_dirichlet_alpha[[2]][1:31, 2] <-  dirichlet_priors_CES_D_neutral
              
          }
          ##
          if (Model_type == "LC_MVOP") settings$reorder_cols_MVP <-  NULL
          ##
          burnin_timing_arguments <-  list()
          if (!is.null(settings$debug_burnin_timing)) {
              if (!is.logical(settings$debug_burnin_timing) || length(settings$debug_burnin_timing) != 1 ||
                  is.na(settings$debug_burnin_timing)) stop("settings$debug_burnin_timing must be TRUE or FALSE.")
              if (settings$debug_burnin_timing) {
                  if (!"debug_burnin_timing" %in% names(formals(BayesMVP::R_fn_sample_model))) {
                      stop("This BayesMVP build does not support debug_burnin_timing; rebuild and reinstall the updated package.")
                  }
                  burnin_timing_arguments <-  list(debug_burnin_timing = TRUE)
              }
          }
          ##
          ## ---- run_in_fresh_R_process: NULL = leave the NicoStan default (a fresh child R process per fit); FALSE = fit in
          ##      this process, which the development-tree sourcing and the cached initialise_model() objects both need (a
          ##      child process loads the installed packages and re-initialises the model):
          ##
          fresh_process_arguments <-  list()
          if (!is.null(settings$run_in_fresh_R_process)) {
              if (!is.logical(settings$run_in_fresh_R_process) || length(settings$run_in_fresh_R_process) != 1 ||
                  is.na(settings$run_in_fresh_R_process)) stop("settings$run_in_fresh_R_process must be TRUE, FALSE or NULL.")
              if (!"run_in_fresh_R_process" %in% names(formals(BayesMVP::R_fn_sample_model))) {
                  stop("This NicoStan / BayesMVP build has no run_in_fresh_R_process argument; reinstall NicoStan and BayesMVP together, restart R.")
              }
              fresh_process_arguments <-  list(run_in_fresh_R_process = settings$run_in_fresh_R_process)
          }
          ##
          ## ---- trajectory jitter: requested values, and the tau_sampling_scale NicoStan will actually apply ("none" when inert).
          ##      The SAME helpers build the file name, so the name, the sampler arguments and the checks below agree:
          ##
          {
              validated_tau_jitter_settings <-  fn_ps7_validate_tau_jitter_settings(settings = settings)
              randomize_tau_burnin_requested <-  validated_tau_jitter_settings$randomize_tau_burnin
              tau_sampling_scale_requested <-  fn_ps7_effective_tau_sampling_scale(settings = settings)
              tau_adam_bias_correction_requested <-  validated_tau_jitter_settings$tau_adam_bias_correction
              ## tau_adaptation_block ("_tbJ";, EXPERIMENTAL): the effective value NicoStan will actually apply
              ## ("main" when inert, i.e. tau is pinned), computed by the SAME helper the file name and collapse use:
              tau_adaptation_block_requested <-  fn_ps7_effective_tau_adaptation_block(settings = settings)
              if (!all(c("tau_sampling_scale", "randomize_tau_burnin") %in% names(formals(BayesMVP::R_fn_sample_model)))) {
                  stop("This NicoStan / BayesMVP build has no tau_sampling_scale / randomize_tau_burnin argument; reinstall NicoStan and BayesMVP together, restart R.")
              }
              if (!"tau_adaptation_block" %in% names(formals(BayesMVP::R_fn_sample_model))) {
                  stop("This NicoStan / BayesMVP build has no tau_adaptation_block argument; reinstall NicoStan and BayesMVP together, restart R.")
              }
              ## step-size acceptance mean and divergence-triggered tau shrink ("_ta" digit), checked BEFORE the fit,
              ## so an older build stops here instead of after a completed run:
              validated_eps_tau_shrink_settings <-  fn_ps7_validate_eps_tau_shrink_settings(settings = settings)
              eps_acceptance_mean_requested <-  validated_eps_tau_shrink_settings$eps_acceptance_mean
              tau_shrink_on_divergence_requested <-  validated_eps_tau_shrink_settings$tau_shrink_on_divergence
              if (!all(c("eps_acceptance_mean", "tau_shrink_on_divergence") %in% names(formals(BayesMVP::R_fn_sample_model)))) {
                  stop("This NicoStan / BayesMVP build has no eps_acceptance_mean / tau_shrink_on_divergence argument; reinstall NicoStan and BayesMVP together, restart R.")
              }
              ## pooled metric estimator: window resets and off-diagonal shrinkage (the "_Mp" code), checked BEFORE the fit, so an
              ## older build stops here instead of after a completed run. Passed for every estimator; only "pooled" reads them:
              validated_metric_pooled_settings <-  fn_ps7_validate_metric_pooled_settings(settings = settings)
              metric_pooled_window_resets_requested <-  validated_metric_pooled_settings$metric_pooled_window_resets
              metric_pooled_offdiagonal_shrinkage_requested <-  validated_metric_pooled_settings$metric_pooled_offdiagonal_shrinkage
              if (!all(c("metric_pooled_window_resets", "metric_pooled_offdiagonal_shrinkage") %in% names(formals(BayesMVP::R_fn_sample_model)))) {
                  stop("This NicoStan / BayesMVP build has no metric_pooled_window_resets / metric_pooled_offdiagonal_shrinkage argument; reinstall NicoStan and BayesMVP together, restart R.")
              }
              ## effective starting tau and learning-rate hold value the sampler should report back:
              tau_initial_requested <-  if_null_then_set_to(settings$tau_initial,
                                                           if (identical(settings$metric_estimator, "pooled")) pi else 2 * pi)
              learning_rate_initial_requested <-  if_null_then_set_to(settings$learning_rate_initial, settings$learning_rate)
          }
          ##
          ## ---- CHESSR_time / SNAPER_time ("_bact" / "_bast"): the time-criterion settings handed to NicoStan as ONE argument,
          ##      time_criterion_settings (a named list, see NicoStan's fn_default_time_criterion_settings()), checked BEFORE the fit so
          ##      an older build stops here. It is passed ONLY for these two algorithms (and only when tau is adapted), so the sampler
          ##      call of every other algorithm is unchanged. NicoStan applies the precedence for each sampling time (user-supplied
          ##      number > the previous run(s) > the sampling timing probe) and reads time_criterion_previous_run_path itself.
          ##      n_iter_sampling_for_time_criterion: "planned_n_iter" passes nothing (NicoStan then uses the planned n_iter of the
          ##      run); "ESS_target_for_time_criterion" passes the target, and ESS_per_iter_sampling_expected when it is user-supplied
          ##      (otherwise NicoStan takes it from the previous run(s), which must then be given):
          ##
          time_criterion_is_used <-  fn_ps7_is_time_criterion_algorithm(burnin_algorithm = burnin_algorithm_used) && !isTRUE(manual_tau)
          time_criterion_arguments <-  list()
          time_criterion_inputs <-  NULL
          if (time_criterion_is_used) {
              if (!exists(x = "fn_ps7_time_criterion_values", mode = "function")) {
                  stop("Source functions/fn_ps7_extract_tau_sweep.R before fitting CHESSR_time / SNAPER_time (it defines fn_ps7_time_criterion_values).")
              }
              if (!"time_criterion_settings" %in% names(formals(BayesMVP::R_fn_sample_model))) {
                  stop(paste0("This NicoStan / BayesMVP build has no time_criterion_settings argument (needed for burnin_algorithm = '",
                              burnin_algorithm_used, "'); reinstall NicoStan and BayesMVP together, restart R."))
              }
              ##
              n_iter_sampling_uses_ESS_target <-  identical(settings$n_iter_sampling_for_time_criterion_source, "ESS_target_for_time_criterion")
              if (n_iter_sampling_uses_ESS_target && is.null(settings$ESS_per_iter_sampling_expected) && is.null(settings$time_criterion_previous_run_path)) {
                  stop(paste0("n_iter_sampling_for_time_criterion_source = 'ESS_target_for_time_criterion' needs ESS_per_iter_sampling_expected: ",
                              "set it, or set time_criterion_previous_run_path to one or more saved runs."))
              }
              sampling_time_sources <-  vapply(X = ps7_sampling_timing_quantity_names, FUN = function(sampling_time_name) {
                  if (!is.null(settings[[sampling_time_name]])) return("user_supplied")
                  if (!is.null(settings$time_criterion_previous_run_path)) return("previous_run")
                  if (isTRUE(settings$run_sampling_timing_probe)) return("sampling_timing_probe")
                  return("none")
              }, FUN.VALUE = character(1))
              if (any(sampling_time_sources == "none")) {
                  warning(paste0("burnin_algorithm = '", burnin_algorithm_used, "' with run_sampling_timing_probe = FALSE and no value for ",
                                 paste(names(sampling_time_sources)[sampling_time_sources == "none"], collapse = ", "),
                                 ": NicoStan falls back to burnin_to_sampling_leapfrog_time_ratio = 0 and sampling_overhead_in_leapfrog_steps = 0, ",
                                 "i.e. the ", sub(pattern = "_time$", replacement = "", x = burnin_algorithm_used), " criterion."), call. = FALSE)
              }
              ##
              time_criterion_settings_for_sampler <-  list(run_sampling_timing_probe = isTRUE(settings$run_sampling_timing_probe),
                                                           sampling_timing_probe_L_values = settings$sampling_timing_probe_L_values,
                                                           sampling_timing_probe_n_iter_per_L = settings$sampling_timing_probe_n_iter_per_L,
                                                           time_criterion_previous_run_path = settings$time_criterion_previous_run_path,
                                                           time_per_leapfrog_step_sampling = settings$time_per_leapfrog_step_sampling,
                                                           time_per_iter_overhead_sampling = settings$time_per_iter_overhead_sampling,
                                                           time_per_iter_summaries_sampling = settings$time_per_iter_summaries_sampling,
                                                           ESS_target_for_time_criterion = if (n_iter_sampling_uses_ESS_target) settings$ESS_target_for_time_criterion else NULL,
                                                           ESS_per_iter_sampling_expected = if (n_iter_sampling_uses_ESS_target) settings$ESS_per_iter_sampling_expected else NULL,
                                                           ## NULL (inert) = NicoStan's default, "diagnostic":
                                                           time_criterion_ess_parameter_set = settings$time_criterion_ess_parameter_set)
              ## NULL = not given (the NicoStan default applies), so it is left out of the list:
              time_criterion_settings_for_sampler <-  time_criterion_settings_for_sampler[!vapply(X = time_criterion_settings_for_sampler,
                                                                                                   FUN = is.null, FUN.VALUE = TRUE)]
              time_criterion_arguments <-  list(time_criterion_settings = time_criterion_settings_for_sampler)
              ##
              ## what PS7 handed to the sampler, and where each sampling time is expected to come from (saved with the run):
              time_criterion_inputs <-  list(effective_time_criterion_settings = settings[intersect(x = ps7_time_criterion_setting_names, y = names(settings))],
                                             time_criterion_settings_text = fn_ps7_time_criterion_settings_text(effective_time_criterion_settings = effective_time_criterion_settings),
                                             time_criterion_settings_hash = fn_ps7_time_criterion_settings_hash(effective_time_criterion_settings = effective_time_criterion_settings),
                                             time_criterion_settings_for_sampler = time_criterion_settings_for_sampler,
                                             sampling_time_sources_expected = as.list(sampling_time_sources),
                                             n_iter_sampling_for_time_criterion_source = settings$n_iter_sampling_for_time_criterion_source)
              ##
              message(BayesMVP:::colourise(paste0("Time criterion (", burnin_algorithm_used, "): run_sampling_timing_probe = ",
                                                  isTRUE(settings$run_sampling_timing_probe),
                                                  if (isTRUE(settings$run_sampling_timing_probe))
                                                      paste0(" (", settings$sampling_timing_probe_n_iter_per_L, " iterations at each number of leapfrog steps in sampling_timing_probe_L_values = ",
                                                             paste(settings$sampling_timing_probe_L_values, collapse = ", "), ")") else "",
                                                  " | sampling times from: ", paste(names(sampling_time_sources), "=", sampling_time_sources, collapse = ", "),
                                                  " | n_iter_sampling_for_time_criterion from: ", settings$n_iter_sampling_for_time_criterion_source,
                                                  if (n_iter_sampling_uses_ESS_target) paste0(" (ESS_target_for_time_criterion = ", settings$ESS_target_for_time_criterion, ")") else
                                                      paste0(" (", settings$n_iter, ")")), "cyan"))
          }
          ps7_step_times["prep"] <-  fn_ps7_stamp()  ## priors, correlation prior sets, per-chain inits, argument checks
          ##
          ## ---- Initialisation object handed to the sampler: R_fn_sample_model() keeps the init_object it is given (and skips
          ##      its own initialise_model()) only when that object was built from the SAME final model_args_list (e.g. the priors,
          ##      LKJ values, lb_corr / ub_corr, corr_force_positive, overflow thresholds, num_chunks = num_chunks_burnin and
          ##      J_grad_option set above), the same stream and the same sample_nuisance / n_nuisance_override as the call below.
          ##      The object from the first cache (built before those settings were added) is only used for its model_args_list,
          ##      so a second one is built here from the final arguments, once per key (model type, N, MCMC seed, final
          ##      model_args_list, sample_nuisance, n_nuisance_override), and reused by every later fit with the same key.
          ##      sampler_init_sample_nuisance and sampler_init_n_nuisance_override are the values passed to R_fn_sample_model()
          ##      below; if they ever differ, the sampler re-initialises the model itself (as before), so nothing fitted changes:
          ##
          sampler_init_sample_nuisance      <-  TRUE
          sampler_init_n_nuisance_override  <-  0
          ##
          sampler_init_object_cache_key <-  paste0("sampler_", Model_type, "_N", N, "_seed", MCMC_seed, "_",
                                                   digest::digest(list( model_args_list     = model_args_list,
                                                                        sample_nuisance     = sampler_init_sample_nuisance,
                                                                        n_nuisance_override = sampler_init_n_nuisance_override)))
          ##
          if (exists(x = sampler_init_object_cache_key, envir = ps7_init_object_cache, inherits = FALSE)) {
                init_object <-  get(x = sampler_init_object_cache_key, envir = ps7_init_object_cache, inherits = FALSE)
                cat(BayesMVP:::colourise(paste0("  Reusing the cached sampler initialisation for seed ", MCMC_seed,
                                                " (num_chunks = ", model_args_list$num_chunks,
                                                ", J_grad_option = ", model_args_list$J_grad_option, ")"), "cyan"), "\n")
          } else {
                init_object <-  BayesMVP::initialise_model( Model_type = Model_type,
                                                  stream = MCMC_seed,
                                                  sample_nuisance = sampler_init_sample_nuisance,
                                                  n_nuisance_override = sampler_init_n_nuisance_override,
                                                  model_args_list = model_args_list,
                                                  compile = TRUE,
                                                  force_recompile = FALSE,
                                                  cmdstanr_model_fit_obj = NULL,
                                                  Stan_data_list = NULL,
                                                  Stan_model_file_path = NULL,
                                                  Stan_cpp_user_header = NULL,
                                                  Stan_cpp_flags = NULL,
                                                  stanc_args = NULL)
                assign(x = sampler_init_object_cache_key, value = init_object, envir = ps7_init_object_cache)
          }
          ##
          ps7_step_times["sampler_init"] <-  fn_ps7_stamp()  ## sampler cache key digest + cache hit or initialise_model() on the final model_args_list
          ##
          ## Call the installed provider so shared NicoStan functions retain their namespace bindings.
          ## tau_initial = "adaptive": the part of [clip_iter, handover] that feeds lambda_max (NicoStan reads it as an option; 1 = all of it):
          options(NicoStan_tau_initial_moments_window_fraction = if (is.null(settings$tau_initial_moments_window_fraction)) 1 else
                                                                   settings$tau_initial_moments_window_fraction)
          ## ADAM settings of the eps and the tau updates (read at the start of run_ps7_models(); every option is set explicitly, so the
          ## sampler arguments beta1_adam / beta2_adam / eps_adam below are not used by either update):
          options(NicoStan_eps_adam_beta1   = adam_settings_requested$eps_adam_beta1,
                  NicoStan_eps_adam_beta2   = adam_settings_requested$eps_adam_beta2,
                  NicoStan_eps_adam_epsilon = adam_settings_requested$eps_adam_epsilon,
                  NicoStan_tau_adam_beta1   = adam_settings_requested$tau_adam_beta1,
                  NicoStan_tau_adam_beta2   = adam_settings_requested$tau_adam_beta2,
                  NicoStan_tau_adam_epsilon = adam_settings_requested$tau_adam_epsilon,
                  NicoStan_tau_learning_rate_restart_at_metric_end = adam_settings_requested$tau_learning_rate_restart_at_metric_end,
                  NicoStan_tau_adaptation_scheme = adam_settings_requested$tau_adaptation_scheme)
          model_results <-  do.call(what = BayesMVP::R_fn_sample_model, args = c(list(debug = FALSE,
                                                stream = MCMC_seed,
                                                ##
                                                init_object = init_object,
                                                ##
                                                n_chains_burnin = settings$n_chains_burnin,
                                                init_lists_per_chain = init_lists_per_chain,
                                                ##
                                                parallel_method = "RcppParallel",
                                                # parallel_method = "OpenMP",
                                                ##
                                                Stan_data_list = NULL,
                                                model_args_list = model_args_list,
                                                ##
                                                reorder_cols_MVP = settings$reorder_cols_MVP,
                                                ##
                                                sample_nuisance = TRUE,
                                                n_nuisance_override = 0,
                                                ##
                                                seed = MCMC_seed,
                                                ##
                                                n_burnin = settings$n_burnin,
                                                n_adapt = n_adapt,
                                                gap = NULL,
                                                ##
                                                n_chains_sampling = settings$n_chains_sampling,
                                                n_superchains = settings$n_superchains, ## -----------------
                                                n_iter = settings$n_iter, ## -----------------
                                                ##
                                                adapt_delta = settings$adapt_delta,
                                                learning_rate = settings$learning_rate,
                                                ##
                                                tau_mult = NULL,
                                                tau_initial = settings$tau_initial,
                                                ##
                                                manual_tau = manual_tau,
                                                tau_if_manual = tau_if_manual,
                                                tau_if_manual_in_L_units = tau_if_manual_in_L_units,
                                                ##
                                                ## "KE" = the kinetic-energy-change statistic every earlier ps7 run used;
                                                ## "ChEES" = the position-based criterion of Hoffman et al. (2021):
                                                burnin_algorithm = burnin_algorithm_used,
                                                tau_weight_by_p_jump = tau_weight_by_p_jump_used,
                                                ## burn-in tau ramp: NULL -> "original" (earlier); "staged" = the later one
                                                tau_ramp = settings$tau_ramp,
                                                ## eps at the ChEES handover: NULL/TRUE -> re-initialised (as before); FALSE -> kept
                                                eps_reinit_at_ChEES_handover = settings$eps_reinit_at_ChEES_handover,
                                                ## nuisance Gaussian centre: NULL -> "running_mean_frozen" (the fix); "running_mean" = previously
                                                theta_hat_us_rule = settings$theta_hat_us_rule,
                                                theta_hat_us_freeze_iter = settings$theta_hat_us_freeze_iter,
                                                burnin_schedule = if (is.null(settings$burnin_schedule)) "legacy" else settings$burnin_schedule,
                                                metric_adaptation_end_iter = settings$metric_adaptation_end_iter,
                                                ## burn-in speed-ups (NULL = unchanged)
                                                pre_burnin_n_iter = settings$pre_burnin_n_iter,
                                                pre_burnin_L = settings$pre_burnin_L,
                                                share_tau_ii_across_chains_in_burnin = settings$share_tau_ii_across_chains_in_burnin,
                                                ## burn-in tau jitter ("_rb1"): with share_tau_ii_across_chains_in_burnin = TRUE, ONE
                                                ## tau_ii ~ U(0, 2 tau) per iteration shared by all burn-in chains:
                                                randomize_tau_burnin = randomize_tau_burnin_requested,
                                                ## one-off rescaling of tau at the switch to sampling ("_tsG"); "none" when inert:
                                                tau_sampling_scale = tau_sampling_scale_requested,
                                                ## which block(s) the trajectory-length criterion is computed on ("_tbJ" = "joint";
                                                ##, EXPERIMENTAL); "main" when inert (tau pinned):
                                                tau_adaptation_block = tau_adaptation_block_requested,
                                                ## burn-in step size adapted on the "harmonic" (default), "arithmetic" or "geometric" mean of the
                                                ## per-chain acceptance probabilities ("_ta" digit, see fn_ps7_tau_adaptation_version_token):
                                                eps_acceptance_mean = eps_acceptance_mean_requested,
                                                ## divergence-triggered 0.95 tau shrink: FALSE (default) or TRUE ("_ta" digit):
                                                tau_shrink_on_divergence = tau_shrink_on_divergence_requested,
                                                burnin_TBB_pool_equals_n_chains = settings$burnin_TBB_pool_equals_n_chains,
                                                ##
                                                ## Sampler storage is controlled independently of summary export below.
                                                store_log_lik_trace = settings$store_log_lik_trace,
                                                ##
                                                ##
                                                diffusion_HMC   = settings$diffusion_HMC, ## -----------------
                                                partitioned_HMC = settings$partitioned_HMC,
                                                ##
                                                ## joint-diffusion integrator ordering; NULL -> package default:
                                                diffusion_HMC_integrator = if (is.null(settings$diffusion_HMC_integrator))
                                                                               "kick_flow_kick" else settings$diffusion_HMC_integrator,
                                                ##
                                                clip_iter = settings$clip_iter,
                                                clip_iter_tau = clip_iter_tau,
                                                ##
                                                ## burnin print INTERVAL (every n_refresh iterations); NULL -> n_burnin / 20:
                                                n_refresh = settings$n_refresh,
                                                ##
                                                ## hold the ADAM learning rate HIGH for the first
                                                ## learning_rate_initial_iter iterations, then drop to
                                                ## settings$learning_rate (NULL / NA = no hold):
                                                learning_rate_initial = settings$learning_rate_initial, ## -------------------------
                                                learning_rate_initial_iter = settings$learning_rate_initial_iter, ## ------------------------
                                                ##
                                                ## the analogous warm start on the step-size itself (NULL = off):
                                                eps_initial = settings$eps_initial,
                                                eps_initial_iter = settings$eps_initial_iter,
                                                ##
                                                # use_proposed = TRUE, ## -----------------
                                                use_proposed = FALSE, ## -----------------
                                                ##
                                                # beta1_adam = 0.00,
                                                # beta1_adam = 0.90,
                                                ## unused: the NicoStan_eps_adam_* / NicoStan_tau_adam_* options set above give both updates their settings
                                                beta1_adam = 0,
                                                ##
                                                # beta2_adam = 0.95,
                                                # eps_adam = 1e-8,
                                                beta2_adam = 0.95,
                                                eps_adam = 1e-8,
                                                ##
                                                force_autodiff   = FALSE,
                                                force_PartialLog = FALSE,
                                                multi_attempts   = settings$multi_attempts, ## ---------------------------------------------------
                                                ##
                                                force_autodiff_for_metric = FALSE, ## -----------------
                                                force_PartialLog_for_metric = FALSE,
                                                force_multi_attempts_for_metric = TRUE,
                                                ##
                                                vect_type = vect_type,
                                                Phi_type = "Phi",
                                                inv_Phi_type = "inv_Phi",
                                                ##
                                                metric_type_main = settings$metric_type_main, ## "Empirical",
                                                metric_shape_main = settings$metric_shape_main, ##  "dense",
                                                ratio_M_main = settings$ratio_M_main,
                                                interval_width_main = settings$int_width, ## max(1, floor(0.5*settings$n_burnin/200)),
                                                ##
                                                M_decay_type  = settings$M_decay_type,
                                                M_decay_power = settings$M_decay_power,
                                                M_decay_scale = settings$M_decay_scale,
                                                ##
                                                metric_type_nuisance = settings$metric_type_nuisance,
                                                # metric_type_nuisance = "Empirical",
                                                # metric_type_nuisance = "uniform_diag",
                                                # metric_type_nuisance = ifelse(settings$diffusion_HMC == TRUE, "unit", "uniform_diag"),
                                                ##
                                                metric_shape_nuisance = "diag",
                                                ratio_M_nuisance = settings$ratio_M_nuisance,
                                                interval_width_nuisance = settings$int_width, ## max(1, floor(0.5*settings$n_burnin/200)),
                                                ##
                                                max_tau_main     = 50,
                                                max_tau_nuisance = 50,
                                                ##
                                                max_eps_main     = 1.00,
                                                max_eps_nuisance = 1.00,
                                                ##
                                                max_L = 1024,
                                                n_nuisance_to_track = settings$n_nuisance_to_track,
                                                ##
                                                use_disk = use_disk,
                                                use_disk_path = use_disk_path,
                                                ##
                                                n_threads_WCP_burnin = settings$n_threads_WCP_burnin,
                                                n_threads_WCP_sampling = settings$n_threads_WCP_sampling,
                                                num_chunks_burnin = settings$num_chunks_burnin,
                                                num_chunks_sampling = settings$num_chunks_sampling,
                                                ##
                                                metric_estimator = settings$metric_estimator,
                                                ## pooled estimator only ("_Mp" code; inert for "chain_mean" / "chain_mean_scaled"):
                                                ## Welford window resets ("none" / "stan_style" / iterations) and off-diagonal shrinkage:
                                                metric_pooled_window_resets = metric_pooled_window_resets_requested,
                                                metric_pooled_offdiagonal_shrinkage = metric_pooled_offdiagonal_shrinkage_requested,
                                                ##
                                                ## test_perm_override = settings$test_perm_override), burnin_timing_arguments, fresh_process_arguments))
                                                test_perm_override = settings$test_perm_override), burnin_timing_arguments, fresh_process_arguments,
                                                ## time-criterion arguments of CHESSR_time / SNAPER_time; empty for every other algorithm:
                                                time_criterion_arguments,
                                                fn_ps7_trajectory_option_arguments(settings = settings)))
          ##
          ps7_step_times["sampler_call"] <-  fn_ps7_stamp()  ## the whole R_fn_sample_model() call (initialisations + pre-burn-in + burn-in + sampling)
          trajectory_option_arguments <-  fn_ps7_trajectory_option_arguments(settings = settings)
          for (option_name in names(trajectory_option_arguments)) {
              requested_value <-  trajectory_option_arguments[[option_name]]
              reported_value <-  model_results[[option_name]]
              matches_request <-  if (is.numeric(requested_value)) {
                  length(reported_value) == 1 && is.numeric(reported_value) &&
                  isTRUE(all.equal(as.numeric(reported_value), as.numeric(requested_value), tolerance = 1e-12))
              } else identical(reported_value, requested_value)
              if (!matches_request) {
                  stop(option_name, ": requested ", paste(requested_value, collapse = ", "),
                       " but the sampler reported ", paste(reported_value, collapse = ", "), ".")
              }
          }
          ##
          ## ---- the J_grad_option the C++ actually ran with (row "J_grad_option" of Model_args_strings) must be the one
          ##      requested - otherwise this run would be saved under the wrong file name:
          ##
          J_grad_option_used <-  as.character(model_results$init_object$Model_args_as_Rcpp_List$Model_args_strings["J_grad_option", 1])
          if (!identical(J_grad_option_used, J_grad_option_requested)) {
            stop("J_grad_option: requested '", J_grad_option_requested, "' but the model ran with '", J_grad_option_used, "'.")
          }
          ##
          if ((isTRUE(settings$autodiff_fallback) || !is.null(model_results$autodiff_fallback)) &&
              !identical(model_results$autodiff_fallback, settings$autodiff_fallback)) {
            stop("Autodiff fallback differs from the requested setting. Re-source the updated sampler R functions and load the rebuilt package.")
          }
          ##
          ## ---- the trajectory-length algorithm and p_jump weighting the sampler actually ran with must be the ones
          ##      requested - both are in the file name ("_ba.." and "_tw1"), so a mismatch would save this run under the wrong name:
          ##
          if (!identical(model_results$burnin_algorithm, burnin_algorithm_used)) {
            stop("burnin_algorithm: requested '", burnin_algorithm_used, "' but the sampler ran with '",
                 paste(model_results$burnin_algorithm, collapse = ", "), "'.")
          }
          if (!identical(isTRUE(model_results$tau_weight_by_p_jump), tau_weight_by_p_jump_used)) {
            stop("tau_weight_by_p_jump: requested ", tau_weight_by_p_jump_used, " but the sampler ran with ",
                 paste(model_results$tau_weight_by_p_jump, collapse = ", "), ".")
          }
          message(paste0("Effective burnin_algorithm = ", model_results$burnin_algorithm,
                         " | trajectory_criterion = ", model_results$trajectory_criterion,
                         " | tau_weight_by_p_jump = ", model_results$tau_weight_by_p_jump))
          ##
          ## ---- the nuisance mass and the centre of the nuisance rotation the sampler actually ran with must be the ones
          ##      requested - both are in the file name ("_nu.." and "_cF"/"_c0"). model_results$theta_hat_us_rule is NULL
          ##      for a NicoStan build older, so that stops too (reinstall NicoStan and BayesMVP together):
          ##
          metric_type_nuisance_requested <-  if (is.null(settings$metric_type_nuisance)) "Empirical" else settings$metric_type_nuisance
          theta_hat_us_rule_requested    <-  if (is.null(settings$theta_hat_us_rule)) "running_mean_frozen" else settings$theta_hat_us_rule
          if (!identical(model_results$metric_type_nuisance, metric_type_nuisance_requested)) {
            stop("metric_type_nuisance: requested '", metric_type_nuisance_requested, "' but the sampler ran with '",
                 paste(model_results$metric_type_nuisance, collapse = ", "), "'.")
          }
          if (!identical(model_results$theta_hat_us_rule, theta_hat_us_rule_requested)) {
            stop("theta_hat_us_rule: requested '", theta_hat_us_rule_requested, "' but the sampler ran with '",
                 paste(model_results$theta_hat_us_rule, collapse = ", "),
                 "' (empty = NicoStan build older than 2026-09-22; reinstall NicoStan and BayesMVP together).")
          }
          message(paste0("Effective metric_type_nuisance = ", model_results$metric_type_nuisance,
                         " | theta_hat_us_rule = ", model_results$theta_hat_us_rule))
          ##
          ## ---- trajectory jitter, tau rescaling, ADAM bias correction, starting tau and learning-rate hold: the values the
          ##      sampler actually ran with must be the ones requested. The first three are in the file name ("_rb1", "_tsG",
          ##      "_ab1"), tau_initial and learning_rate_initial too ("_ti..", "_Li.."), so any mismatch would save this run
          ##      under the wrong name. NULL = a NicoStan build older (reinstall NicoStan and BayesMVP together):
          ##
          {
              fn_numbers_match <-  function(value_reported, value_requested) {
                  if (is.null(value_reported) || length(value_reported) != 1 || length(value_requested) != 1) return(FALSE)
                  if (is.na(value_requested)) return(is.na(value_reported))
                  return(!is.na(value_reported) && isTRUE(all.equal(as.numeric(value_reported), as.numeric(value_requested), tolerance = 1e-12)))
              }
              ##
              if (!identical(model_results$randomize_tau_burnin, randomize_tau_burnin_requested)) {
                  stop("randomize_tau_burnin: requested ", randomize_tau_burnin_requested, " but the sampler ran with ",
                       paste(model_results$randomize_tau_burnin, collapse = ", "), ".")
              }
              tau_sampling_scale_effective_expected <-  if (is.numeric(tau_sampling_scale_requested))
                  as.character(tau_sampling_scale_requested) else tau_sampling_scale_requested
              if (!identical(model_results$tau_sampling_scale_effective, tau_sampling_scale_effective_expected)) {
                  stop("tau_sampling_scale: expected the sampler to apply '", tau_sampling_scale_effective_expected, "' but it reports '",
                       paste(model_results$tau_sampling_scale_effective, collapse = ", "), "' (not applied because: ",
                       paste(model_results$tau_sampling_scale_not_applied_reason, collapse = ", "), ").")
              }
              tau_sampling_scale_factor_expected <-  if (identical(tau_sampling_scale_requested, "none")) 1 else
                                                    if (identical(tau_sampling_scale_requested, "gaussian_matched"))
                                                        fn_ps7_expected_gaussian_tau_sampling_factor(burnin_algorithm = burnin_algorithm_used) else
                                                    tau_sampling_scale_requested
              if (is.null(model_results$tau_sampling_scale_factor) ||
                  !isTRUE(all.equal(model_results$tau_sampling_scale_factor, tau_sampling_scale_factor_expected, tolerance = 1e-6))) {
                  stop("tau_sampling_scale: expected factor ", signif(tau_sampling_scale_factor_expected, 7), " but the sampler reports ",
                       paste(model_results$tau_sampling_scale_factor, collapse = ", "), ".")
              }
              if (!isTRUE(all.equal(model_results$tau_main_after_sampling_scale,
                                    tau_sampling_scale_factor_expected * model_results$tau_main_before_sampling_scale,
                                    tolerance = 1e-10))) {
                  stop("tau_sampling_scale: tau_main went from ", model_results$tau_main_before_sampling_scale, " to ",
                       model_results$tau_main_after_sampling_scale, ", not by the factor ", signif(tau_sampling_scale_factor_expected, 7),
                       " (capped at max_tau_main?).")
              }
              ## the summaries compute the sampling L / gradient counts from this list, so it must hold the SAMPLING tau:
              if (!isTRUE(all.equal(model_results$burnin_object$EHMC_args_as_Rcpp_List$tau_main,
                                    model_results$tau_main_after_sampling_scale,
                                    tolerance = 1e-12))) {
                  stop("tau_sampling_scale: burnin_object$EHMC_args_as_Rcpp_List$tau_main does not hold the sampling tau.")
              }
              ##
              tau_adam_bias_correction_reported <-  if (is.null(model_results$tau_adam_bias_correction)) "iteration_index" else
                                                   model_results$tau_adam_bias_correction
              if (!identical(tau_adam_bias_correction_reported, tau_adam_bias_correction_requested)) {
                  stop("tau_adam_bias_correction: requested '", tau_adam_bias_correction_requested, "' but the sampler reports '",
                       tau_adam_bias_correction_reported, "' (iteration_index = field absent, i.e. a NicoStan build older than the fix).")
              }
              ##
              ## eps_acceptance_mean / tau_shrink_on_divergence: the values NicoStan reports it applied (model_results, else the
              ## burn-in object). NULL = a build without these options, treated as "arithmetic" / TRUE (the earlier behaviour):
              eps_acceptance_mean_reported <-  BayesMVP::if_null_then_set_to(model_results$eps_acceptance_mean,
                                                                             BayesMVP::if_null_then_set_to(model_results$burnin_object$eps_acceptance_mean, "arithmetic"))
              tau_shrink_on_divergence_reported <-  BayesMVP::if_null_then_set_to(model_results$tau_shrink_on_divergence,
                                                                                  BayesMVP::if_null_then_set_to(model_results$burnin_object$tau_shrink_on_divergence, TRUE))
              if (!identical(eps_acceptance_mean_reported, eps_acceptance_mean_requested)) {
                  stop("eps_acceptance_mean: requested '", eps_acceptance_mean_requested, "' but the sampler reports '",
                       eps_acceptance_mean_reported, "' ('arithmetic' = field absent, i.e. a NicoStan build without eps_acceptance_mean).")
              }
              if (!identical(as.logical(tau_shrink_on_divergence_reported), tau_shrink_on_divergence_requested)) {
                  stop("tau_shrink_on_divergence: requested ", tau_shrink_on_divergence_requested, " but the sampler reports ",
                       tau_shrink_on_divergence_reported, " (TRUE = field absent, i.e. a NicoStan build without tau_shrink_on_divergence).")
              }
              ##
              message(BayesMVP:::colourise(paste0("Effective eps_acceptance_mean = ", eps_acceptance_mean_reported,
                                                  " | tau_shrink_on_divergence = ", tau_shrink_on_divergence_reported), "cyan"))
              ##
              ## pooled metric estimator: the values NicoStan reports it applied (model_results, else the burn-in object). NULL = a
              ## build without these options, treated as "stan_style" / 0 (the earlier behaviour):
              metric_pooled_window_resets_reported <-  BayesMVP::if_null_then_set_to(model_results$metric_pooled_window_resets,
                                                                                     BayesMVP::if_null_then_set_to(model_results$burnin_object$metric_pooled_window_resets, "stan_style"))
              metric_pooled_offdiagonal_shrinkage_reported <-  BayesMVP::if_null_then_set_to(model_results$metric_pooled_offdiagonal_shrinkage,
                                                                                             BayesMVP::if_null_then_set_to(model_results$burnin_object$metric_pooled_offdiagonal_shrinkage, 0))
              metric_pooled_window_reset_iterations_reported <-  BayesMVP::if_null_then_set_to(model_results$metric_pooled_window_reset_iterations,
                                                                                               model_results$burnin_object$metric_pooled_window_reset_iterations)
              if (!identical(fn_ps7_metric_pooled_window_resets_label(metric_pooled_window_resets = metric_pooled_window_resets_reported),
                             fn_ps7_metric_pooled_window_resets_label(metric_pooled_window_resets = metric_pooled_window_resets_requested))) {
                  stop("metric_pooled_window_resets: requested '", fn_ps7_metric_pooled_window_resets_label(metric_pooled_window_resets = metric_pooled_window_resets_requested),
                       "' but the sampler reports '", fn_ps7_metric_pooled_window_resets_label(metric_pooled_window_resets = metric_pooled_window_resets_reported),
                       "' ('stan_style' = field absent, i.e. a NicoStan build without metric_pooled_window_resets).")
              }
              metric_pooled_offdiagonal_shrinkage_reported <-  fn_ps7_validate_metric_pooled_settings(
                  settings = list(metric_pooled_window_resets = metric_pooled_window_resets_reported,
                                  metric_pooled_offdiagonal_shrinkage = metric_pooled_offdiagonal_shrinkage_reported))$metric_pooled_offdiagonal_shrinkage
              if (!identical(metric_pooled_offdiagonal_shrinkage_reported, metric_pooled_offdiagonal_shrinkage_requested)) {
                  stop("metric_pooled_offdiagonal_shrinkage: requested ", metric_pooled_offdiagonal_shrinkage_requested, " but the sampler reports ",
                       paste(metric_pooled_offdiagonal_shrinkage_reported, collapse = ", "),
                       " (0 = field absent, i.e. a NicoStan build without metric_pooled_offdiagonal_shrinkage).")
              }
              ##
              message(BayesMVP:::colourise(paste0("Effective metric_pooled_window_resets = ",
                                                  fn_ps7_metric_pooled_window_resets_label(metric_pooled_window_resets = metric_pooled_window_resets_reported),
                                                  if (!is.null(metric_pooled_window_reset_iterations_reported))
                                                      paste0(" (reset iterations: ", fn_ps7_metric_pooled_window_resets_label(
                                                          metric_pooled_window_resets = as.numeric(metric_pooled_window_reset_iterations_reported)), ")") else "",
                                                  " | metric_pooled_offdiagonal_shrinkage = ", metric_pooled_offdiagonal_shrinkage_reported,
                                                  ## if (!identical(settings$metric_estimator, "pooled")) " (inert: metric_estimator is not 'pooled')" else ""), "cyan"))
                                                  if (identical(settings$metric_estimator, "per_iteration")) " (window resets inert: metric_estimator is 'per_iteration')" else
                                                  if (!identical(settings$metric_estimator, "pooled")) " (inert: metric_estimator is not 'pooled')" else ""), "cyan"))
              ##
              ## CHESSR_time / SNAPER_time ("_bact" / "_bast"): the sampler must return its time-criterion record, with
              ## n_iter_sampling_for_time_criterion from the source requested (the planned n_iter, or ESS_target_for_time_criterion /
              ## ESS_per_iter_sampling_expected) and every user-supplied value as given; a probe that did not run or did not fit, or a
              ## fall-back to the rate criterion, gives a warning. Only read for these two algorithms:
              time_criterion_record <-  NULL
              if (time_criterion_is_used) {
                  time_criterion_record <-  fn_ps7_time_criterion_record_from_results(model_results = model_results)
                  if (is.null(time_criterion_record)) {
                      stop("burnin_algorithm = '", burnin_algorithm_used, "' but the sampler returned no time-criterion record ",
                           "(a NicoStan build without CHESSR_time / SNAPER_time); reinstall NicoStan and BayesMVP together, restart R.")
                  }
                  time_criterion_values <-  fn_ps7_time_criterion_values(time_criterion_record = time_criterion_record)
                  fn_as_number <-  function(recorded_value) suppressWarnings(as.numeric(recorded_value))
                  n_iter_sampling_for_time_criterion_reported <-  fn_as_number(time_criterion_values$n_iter_sampling_for_time_criterion)
                  n_iter_sampling_for_time_criterion_expected <-  if (identical(settings$n_iter_sampling_for_time_criterion_source, "ESS_target_for_time_criterion"))
                      settings$ESS_target_for_time_criterion / fn_as_number(time_criterion_values$ESS_per_iter_sampling_expected) else settings$n_iter
                  if (!isTRUE(startsWith(x = as.character(time_criterion_values$n_iter_sampling_for_time_criterion_source),
                                         prefix = settings$n_iter_sampling_for_time_criterion_source))) {
                      stop("n_iter_sampling_for_time_criterion_source: requested '", settings$n_iter_sampling_for_time_criterion_source,
                           "' but the sampler reports '", time_criterion_values$n_iter_sampling_for_time_criterion_source, "'.")
                  }
                  if (!isTRUE(all.equal(n_iter_sampling_for_time_criterion_reported, n_iter_sampling_for_time_criterion_expected, tolerance = 1e-9))) {
                      stop("n_iter_sampling_for_time_criterion: expected ", n_iter_sampling_for_time_criterion_expected, " (",
                           settings$n_iter_sampling_for_time_criterion_source, ") but the sampler reports ",
                           time_criterion_values$n_iter_sampling_for_time_criterion, ".")
                  }
                  for (user_supplied_name in c(ps7_sampling_timing_quantity_names, "ESS_per_iter_sampling_expected", "ESS_target_for_time_criterion")) {
                      if (!is.null(time_criterion_settings_for_sampler[[user_supplied_name]]) &&
                          !isTRUE(all.equal(fn_as_number(time_criterion_values[[user_supplied_name]]), time_criterion_settings_for_sampler[[user_supplied_name]],
                                            tolerance = 1e-12))) {
                          stop(user_supplied_name, ": supplied ", time_criterion_settings_for_sampler[[user_supplied_name]], " but the sampler used ",
                               time_criterion_values[[user_supplied_name]], ".")
                      }
                  }
                  if (isTRUE(settings$run_sampling_timing_probe) &&
                      (!isTRUE(as.logical(time_criterion_values$sampling_timing_probe_ran)) || !identical(time_criterion_values$sampling_timing_probe_status, "fitted"))) {
                      warning(paste0("run_sampling_timing_probe = TRUE but the sampler reports that the probe ",
                                     if (isTRUE(as.logical(time_criterion_values$sampling_timing_probe_ran))) "ran without a fit" else "did not run",
                                     " (status '", time_criterion_values$sampling_timing_probe_status, "')."), call. = FALSE)
                  }
                  if (isTRUE(as.logical(time_criterion_values$time_criterion_fallback_to_rate_criterion))) {
                      warning(paste0("burnin_algorithm = '", burnin_algorithm_used, "': no sampling times were available, so the sampler used the ",
                                     sub(pattern = "_time$", replacement = "", x = burnin_algorithm_used), " criterion (burnin_to_sampling_leapfrog_time_ratio = 0, ",
                                     "sampling_overhead_in_leapfrog_steps = 0)."), call. = FALSE)
                  }
                  fn_text_4_significant_digits <-  function(recorded_value) {
                      if (is.numeric(recorded_value) && length(recorded_value) == 1 && is.finite(recorded_value)) as.character(signif(recorded_value, 4)) else
                          as.character(recorded_value)
                  }
                  message(BayesMVP:::colourise(paste0("Effective time criterion (", burnin_algorithm_used, "):",
                                                      " time_per_leapfrog_step_sampling = ", fn_text_4_significant_digits(time_criterion_values$time_per_leapfrog_step_sampling),
                                                      " s (", time_criterion_values$time_per_leapfrog_step_sampling_source, ")",
                                                      ", time_per_iter_overhead_sampling = ", fn_text_4_significant_digits(time_criterion_values$time_per_iter_overhead_sampling),
                                                      " s (", time_criterion_values$time_per_iter_overhead_sampling_source, ")",
                                                      ", time_per_iter_summaries_sampling = ", fn_text_4_significant_digits(time_criterion_values$time_per_iter_summaries_sampling),
                                                      " s (", time_criterion_values$time_per_iter_summaries_sampling_source, ")",
                                                      " | sampling_overhead_in_leapfrog_steps = ", fn_text_4_significant_digits(time_criterion_values$sampling_overhead_in_leapfrog_steps),
                                                      " | probe: ", time_criterion_values$sampling_timing_probe_status,
                                                      ", wall time ", fn_text_4_significant_digits(time_criterion_values$sampling_timing_probe_wall_time), " s",
                                                      " (", fn_text_4_significant_digits(100 * fn_as_number(time_criterion_values$sampling_timing_probe_fraction_of_burnin)), "% of the burn-in)"), "cyan"))
                  message(BayesMVP:::colourise(paste0("  time_per_leapfrog_step_burnin = ", fn_text_4_significant_digits(time_criterion_values$time_per_leapfrog_step_burnin_at_handover),
                                                      " s at the handover, ", fn_text_4_significant_digits(time_criterion_values$time_per_leapfrog_step_burnin), " s at the end (",
                                                      time_criterion_values$time_per_leapfrog_step_burnin_status, ")",
                                                      " | n_iter_burnin = ", time_criterion_values$n_iter_burnin,
                                                      ", n_iter_sampling_for_time_criterion = ", fn_text_4_significant_digits(n_iter_sampling_for_time_criterion_reported),
                                                      " | at the handover / at the end:",
                                                      " tau_offset_from_sampling_overhead = ", fn_text_4_significant_digits(time_criterion_values$tau_offset_from_sampling_overhead_at_handover),
                                                      " / ", fn_text_4_significant_digits(time_criterion_values$tau_offset_from_sampling_overhead_at_end),
                                                      ", burnin_to_sampling_leapfrog_time_ratio = ", fn_text_4_significant_digits(time_criterion_values$burnin_to_sampling_leapfrog_time_ratio_at_handover),
                                                      " / ", fn_text_4_significant_digits(time_criterion_values$burnin_to_sampling_leapfrog_time_ratio_at_end),
                                                      ", time_to_target_ESS_tau_penalty = ", fn_text_4_significant_digits(time_criterion_values$time_to_target_ESS_tau_penalty_at_handover),
                                                      " / ", fn_text_4_significant_digits(time_criterion_values$time_to_target_ESS_tau_penalty_at_end)), "cyan"))
                  message(BayesMVP:::colourise(paste0("  over the second half of the tau updates: ESS_elasticity_wrt_log_tau (pooled estimate) = ",
                                                      fn_text_4_significant_digits(time_criterion_values$ESS_elasticity_wrt_log_tau_pooled_over_second_half_of_tau_updates),
                                                      ", mean time_to_target_ESS_tau_penalty = ",
                                                      fn_text_4_significant_digits(time_criterion_values$time_to_target_ESS_tau_penalty_mean_over_second_half_of_tau_updates),
                                                      " | probe summaries: ", time_criterion_values$sampling_timing_probe_summaries_method), "cyan"))
              }
              ##
              ## tau_adaptation_block ("_tbJ";, EXPERIMENTAL): NULL = a NicoStan build without this option, treated
              ## as "main" (exactly like the tau_adam_bias_correction check above):
              tau_adaptation_block_reported <-  if (is.null(model_results$tau_adaptation_block)) "main" else
                                               model_results$tau_adaptation_block
              if (!identical(tau_adaptation_block_reported, tau_adaptation_block_requested)) {
                  stop("tau_adaptation_block: requested '", tau_adaptation_block_requested, "' but the sampler reports '",
                       tau_adaptation_block_reported, "' (NULL/'main' = a NicoStan build without tau_adaptation_block).")
              }
              ##
              ## ADAM ("_A" code): the settings the burn-in reports the eps and the tau updates used (burnin_object, "<setting>_used"). NULL =
              ## a NicoStan build without that option, which used the sampler arguments 0 / 0.95 / 1e-8: accepted only when that was requested.
              ## The expected tau beta1 is ps7_tau_adam_beta1 (0 when tau is pinned, where it is inert):
              adam_settings_default_values <-  c(eps_adam_beta1 = 0, eps_adam_beta2 = 0.95, eps_adam_epsilon = 1e-8,
                                                 tau_adam_beta1 = 0, tau_adam_beta2 = 0.95, tau_adam_epsilon = 1e-8)
              adam_settings_reported <-  list()
              for (adam_setting_name in names(adam_settings_default_values)) {
                  adam_setting_reported <-  model_results$burnin_object[[paste0(adam_setting_name, "_used")]]
                  if (is.null(adam_setting_reported)) adam_setting_reported <-  unname(adam_settings_default_values[adam_setting_name])
                  if (!fn_numbers_match(adam_setting_reported, adam_settings_requested[[adam_setting_name]])) {
                      stop(adam_setting_name, ": expected ", adam_settings_requested[[adam_setting_name]], " but the burn-in reports ", adam_setting_reported,
                           " (", adam_settings_default_values[adam_setting_name], " = a NicoStan build without that option; re-source the NicoStan R files).")
                  }
                  adam_settings_reported[[adam_setting_name]] <-  as.numeric(adam_setting_reported)
              }
              tau_adam_beta1_reported <-  adam_settings_reported$tau_adam_beta1
              tau_learning_rate_restart_reported <-  if (is.null(model_results$burnin_object$tau_learning_rate_restart_at_metric_end)) FALSE else
                                                    as.logical(model_results$burnin_object$tau_learning_rate_restart_at_metric_end)
              if (!identical(tau_learning_rate_restart_reported, adam_settings_requested$tau_learning_rate_restart_at_metric_end)) {
                  stop("tau learning-rate restart at the metric freeze: requested ", adam_settings_requested$tau_learning_rate_restart_at_metric_end,
                       " but the burn-in reports ", tau_learning_rate_restart_reported,
                       " (FALSE = a NicoStan build without NicoStan_tau_learning_rate_restart_at_metric_end; re-source the NicoStan R files).")
              }
              tau_adaptation_scheme_reported <-  if (is.null(model_results$burnin_object$tau_adaptation_scheme)) "adam_decay" else
                                                model_results$burnin_object$tau_adaptation_scheme
              if (!identical(tau_adaptation_scheme_reported, adam_settings_requested$tau_adaptation_scheme)) {
                  stop("tau adaptation scheme: requested ", adam_settings_requested$tau_adaptation_scheme, " but the burn-in reports ", tau_adaptation_scheme_reported,
                       " ('adam_decay' = a NicoStan build without NicoStan_tau_adaptation_scheme; re-source the NicoStan R files).")
              }
              message(BayesMVP:::colourise(paste0("Effective ADAM: eps beta1 = ", adam_settings_reported$eps_adam_beta1, ", beta2 = ", adam_settings_reported$eps_adam_beta2,
                                                  ", denominator constant = ", adam_settings_reported$eps_adam_epsilon,
                                                  " | tau beta1 = ", adam_settings_reported$tau_adam_beta1, ", beta2 = ", adam_settings_reported$tau_adam_beta2,
                                                  ", denominator constant = ", adam_settings_reported$tau_adam_epsilon,
                                                  " | tau learning-rate restart = ", tau_learning_rate_restart_reported,
                                                  " | tau adaptation scheme = ", tau_adaptation_scheme_reported,
                                                  # if (identical(tau_adaptation_scheme_reported, "probe_then_average"))
                                                  if (tau_adaptation_scheme_reported %in% c("probe_then_average", "fixed_length_probe_then_decay_and_average"))
                                                      paste0(" (probe ended at iteration ", model_results$burnin_object$tau_probe_end_iteration,
                                                             ", averaged tau ", signif(model_results$burnin_object$tau_main_averaged, 5), ")") else "",
                                                  if (isTRUE(tau_learning_rate_restart_reported))
                                                      paste0(" (from iteration ", model_results$burnin_object$tau_learning_rate_restart_iteration, ")") else ""), "cyan"))
              ##
              ## tau_initial "adaptive" is reported back as the string itself; a numeric level as the number it ran with:
              ## if (!fn_numbers_match(model_results$tau_initial, tau_initial_requested)) {
              if (!(if (identical(tau_initial_requested, "adaptive")) identical(model_results$tau_initial, "adaptive") else
                        fn_numbers_match(model_results$tau_initial, tau_initial_requested))) {
                  stop("tau_initial: requested ", tau_initial_requested, " but the sampler ran with ",
                       paste(model_results$tau_initial, collapse = ", "), ".")
              }
              if (!fn_numbers_match(model_results$learning_rate_initial, learning_rate_initial_requested)) {
                  stop("learning_rate_initial: requested ", learning_rate_initial_requested, " but the sampler ran with ",
                       paste(model_results$learning_rate_initial, collapse = ", "), ".")
              }
              message(paste0("Effective randomize_tau_burnin = ", model_results$randomize_tau_burnin,
                             " | tau_sampling_scale = ", model_results$tau_sampling_scale_effective,
                             " (factor ", signif(model_results$tau_sampling_scale_factor, 6),
                             "; tau_main ", signif(model_results$tau_main_before_sampling_scale, 6),
                             " -> ", signif(model_results$tau_main_after_sampling_scale, 6), ")",
                             " | tau_adam_bias_correction = ", tau_adam_bias_correction_reported,
                             " | tau_adaptation_block = ", tau_adaptation_block_reported,
                             ## " | tau_initial = ", signif(model_results$tau_initial, 7),
                             " | tau_initial = ", if (is.character(model_results$tau_initial)) model_results$tau_initial else
                                                      signif(model_results$tau_initial, 7),
                             " | learning_rate_initial = ", model_results$learning_rate_initial))
          }
          ##
          ps7_step_times["post_fit_checks"] <-  fn_ps7_stamp()  ## sampler self-checks and messages after R_fn_sample_model()
          ##
          model_fit_object <-  BayesMVP::create_summary_and_traces(  model_results = model_results,
                                                          ##
                                                          compute_main_params = TRUE,
                                                          compute_transformed_parameters = TRUE,
                                                          compute_generated_quantities = TRUE,
                                                          ##
                                                          save_log_lik_trace = FALSE,
                                                          # save_log_lik_trace = TRUE,
                                                          ##
                                                          save_nuisance_trace = FALSE,
                                                          compute_nested_rhat = TRUE,
                                                          ##
                                                          n_superchains = NULL,
                                                          save_trace_tibbles = FALSE,
                                                          ##
                                                          use_disk = use_disk,
                                                          use_disk_path = use_disk_path,
                                                          use_disk_path_post_hoc_dir = use_disk_path_post_hoc_dir)
            ##
            ps7_step_times["summaries"] <-  fn_ps7_stamp()  ## BayesMVP::create_summary_and_traces()
            ##
            ## ---- Extract:
            ##
            tibble_main <-  model_fit_object$summaries$summary_tibbles$summary_tibble_main_params
            tibble_tp   <-  model_fit_object$summaries$summary_tibbles$summary_tibble_transformed_parameters
            tibble_gq   <-  model_fit_object$summaries$summary_tibbles$summary_tibble_generated_quantities
            ##
            ## Ordinal skeletons pad the rectangular baseline arrays beyond each test's real thresholds.
            ## Keep the full summaries/traces, but exclude those deterministic padding cells from diagnostics.
            tibble_gq_for_diagnostics <-  tibble_gq
            if (Model_type == "LC_MVOP") {
                ordinal_rows <-  grepl(pattern = "^(Se|Sp|Fp)_baseline_ord\\[", x = tibble_gq$parameter)
                ordinal_names <-  tibble_gq$parameter[ordinal_rows]
                ordinal_test <-  as.integer(sub(pattern = "^.*\\[([0-9]+),([0-9]+)\\]$", replacement = "\\1", x = ordinal_names))
                ordinal_threshold <-  as.integer(sub(pattern = "^.*\\[([0-9]+),([0-9]+)\\]$", replacement = "\\2", x = ordinal_names))
                real_thresholds <-  model_args_list$n_thr_per_ord_test
                if (is.null(real_thresholds)) real_thresholds <-  model_args_list$n_cat_per_ord_test - 1L
                stopifnot(!anyNA(ordinal_test), !anyNA(ordinal_threshold), all(ordinal_test >= 1L),
                          all(ordinal_test <= length(real_thresholds)))
                retain_rows <-  rep(TRUE, nrow(tibble_gq))
                retain_rows[ordinal_rows] <-  ordinal_threshold <= real_thresholds[ordinal_test]
                tibble_gq_for_diagnostics <-  tibble_gq[retain_rows, , drop = FALSE]
            }
            baseline_families <-  sub(pattern = "\\[.*$", replacement = "", x = tibble_gq_for_diagnostics$parameter)
            tibble_Se <-  tibble_gq_for_diagnostics[baseline_families %in% c("Se_baseline", "Se_baseline_bin", "Se_baseline_ord"), , drop = FALSE]
            tibble_Sp <-  tibble_gq_for_diagnostics[baseline_families %in% c("Sp_baseline", "Sp_baseline_bin", "Sp_baseline_ord"), , drop = FALSE]
            tibble_prev <-  BayesMVP::extract_params_from_tibble_batch(tibble = tibble_gq, param_strings_vec = "p")
            ##
            trace_main <-  model_fit_object$traces$traces_as_arrays$trace_params_main
            ##
            # str(trace_main)
            # ##
            # # trace_main: n_params x n_iter x n_chains (or however the list is laid out)
            # tr <-  trace_main   # 250 x 180 x 44
            # ##
            # BayesMVP:::Rcpp_compute_MCMC_diagnostics(list(tr[, , 1]), diagnostic = "split_ESS_rank", n_threads = 1)$diagnostics
            # BayesMVP:::Rcpp_compute_MCMC_diagnostics(list(tr[, , 1]), diagnostic = "split_rhat_rank", n_threads = 1)$diagnostics
            # ##
            # acc <-  apply(tr, 2, function(m) mean(rowSums(abs(diff(m))) > 0)); summary(acc)
#           ##
#           which(apply(tr, 2, function(m) min(c(apply(m[1:125, ], 2, sd), apply(m[126:250, ], 2, sd)))) == 0)
#           # 
#           # # max over params of within-chain sd, per chain -> frozen chains are exactly 0
#           # chain_sd <-  apply(tr, 2, function(m) max(apply(m, 2, sd)))   # m is 250 x 44
#           # which(chain_sd == 0)
#           # 
#           # # Stan's check is isApproxToConstant at 1e-12 relative, so also catch "barely moved":
#           # chain_rel <-  apply(tr, 2, function(m) max(apply(m, 2, function(v) diff(range(v)) / max(abs(v), 1e-300))))
#           # which(chain_rel < 1e-10)
#           # 
#           # # any non-finite draw anywhere (the other branch of the same check):
#           # which(apply(tr, 2, function(m) any(!is.finite(m))))
#           # 
#           # # divergences per chain (if divs is 250 x 180):
#           colSums(divs)
#           ##
#           bad <-  c(125, 145, 166)
            # for (k in bad) {
            #   m <-  tr[, k, ]
            #   cat("chain", k,
            #       "| params constant in 1st half:", sum(apply(m[1:125, ], 2, sd) == 0),
            #       "| in 2nd half:", sum(apply(m[126:250, ], 2, sd) == 0), "of 44\n")
            # }
            # 
            # ## divergences for those chains (sampling_object[[2]] is the list of 1 x n_iter div traces)
            # divs_per_chain <-  sapply(model_results$sampling_object[[2]], sum)
            # divs_per_chain[bad]
            # 
            # ## where they're stuck vs everyone else (last draw)
            # round(cbind(stuck = t(tr[250, bad, ]), median_all = apply(tr[250, , ], 2, median)), 3)
            ##
          trace_tp   <-  model_fit_object$traces$traces_as_arrays$trace_transformed_params
          trace_gq   <-  model_fit_object$traces$traces_as_arrays$trace_generated_quantities
          ##
          ## ---- Compare the metric with UNCONSTRAINED main draws, in the metric's own coordinate system.
          ## A short high-dimensional validation fit may not have enough draws for this diagnostic.
          metric_inverse <-  model_results$burnin_object$EHMC_Metric_as_Rcpp_List$M_inv_dense_main
          draws_mat <-  do.call(what = rbind, args = lapply(X = model_results$sampling_object[[1]], FUN = t))
          if (is.matrix(metric_inverse) && nrow(draws_mat) > ncol(draws_mat) &&
              identical(dim(metric_inverse), c(ncol(draws_mat), ncol(draws_mat)))) {
              posterior_covariance <-  stats::cov(x = draws_mat)
              generalised_eigenvalues <-  eigen(x = solve(a = metric_inverse, b = posterior_covariance), only.values = TRUE)$values
              message(paste0("Metric/posterior covariance eigenvalue ratio = ",
                             max(Re(generalised_eigenvalues)) / min(Re(generalised_eigenvalues))))
          }
          ##
          ## Preserve undefined diagnostics instead of removing them and reporting a misleading best-case result.
          max_Rhat <-  if (length(x = tibble_gq_for_diagnostics$Rhat) > 0 && !anyNA(x = tibble_gq_for_diagnostics$Rhat)) max(tibble_gq_for_diagnostics$Rhat) else NA_real_
          max_nRhat <-  if (length(x = tibble_gq_for_diagnostics$n_Rhat) > 0 && !anyNA(x = tibble_gq_for_diagnostics$n_Rhat)) max(tibble_gq_for_diagnostics$n_Rhat) else NA_real_
          ##
          min_ESS <-  if (length(x = tibble_gq_for_diagnostics$n_eff) > 0 && all(is.finite(x = tibble_gq_for_diagnostics$n_eff))) min(tibble_gq_for_diagnostics$n_eff) else NA_real_
          ## the tail ESS over the same parameters (NicoStan's n_eff_tail; NA for a summary made before that column existed):
          min_ESS_tail <-  if (length(x = tibble_gq_for_diagnostics$n_eff_tail) > 0 && all(is.finite(x = tibble_gq_for_diagnostics$n_eff_tail))) min(tibble_gq_for_diagnostics$n_eff_tail) else NA_real_
          ## the ESS of the centred squared draws over the same parameters (NicoStan's n_eff_sd, the posterior SD's ESS; Enzo, 4 Oct 2026:
          ## a separate endpoint of equal standing; NA for a summary made before that column existed):
          min_ESS_sd <-  if (length(x = tibble_gq_for_diagnostics$n_eff_sd) > 0 && all(is.finite(x = tibble_gq_for_diagnostics$n_eff_sd))) min(tibble_gq_for_diagnostics$n_eff_sd) else NA_real_
          ##
          time_total  <-  model_fit_object$summaries$efficiency_info$time_total
          time_burnin <-  model_fit_object$summaries$efficiency_info$time_burnin
          ##
          time_sampling <-  model_fit_object$summaries$efficiency_info$time_sampling
          time_summaries <-  model_fit_object$summaries$efficiency_info$time_summaries
          ##
          summary_tibbles <-  model_fit_object$summaries$summary_tibbles
          efficiency_info <-  model_fit_object$summaries$efficiency_info
          divergences <-  model_fit_object$summaries$divergences
          HMC_info <-  model_fit_object$summaries$HMC_info
          ## pooled metric estimator settings, as the sampler reported them (checked above), also kept in HMC_info
          ## (NicoStan's own entries are kept when it already records them):
          if (is.null(HMC_info$metric_pooled_window_resets))          HMC_info$metric_pooled_window_resets <-  metric_pooled_window_resets_reported
          if (is.null(HMC_info$metric_pooled_offdiagonal_shrinkage))  HMC_info$metric_pooled_offdiagonal_shrinkage <-  as.numeric(metric_pooled_offdiagonal_shrinkage_reported)
          ##
          ESS_per_sec_samp <-  min_ESS / time_sampling
          ESS_per_sec_total <-  min_ESS / time_total
          ##
          L_main_during_sampling <-  (HMC_info$tau_main / HMC_info$eps_main)
          ## n_grad_evals_sampling_main <-  L_main_during_sampling * settings$n_iter * settings$n_chains_sampling
          ## executed gradients: whole leapfrog steps under the sampling jitter (fn_ps7_expected_sampling_leapfrog_steps_per_iteration,
          ## fn_ps7_extract_tau_sweep.R), not the nominal tau / eps:
          ## n_grad_evals_sampling_main <-  fn_ps7_expected_sampling_leapfrog_steps_per_iteration( tau = HMC_info$tau_main,
          ##                                                                                       eps = HMC_info$eps_main,
          ##                                                                                       randomize_tau_sampling = if (is.null(settings$randomize_tau_sampling)) TRUE else settings$randomize_tau_sampling) *
          ##                               settings$n_iter * settings$n_chains_sampling
          ## gradient EVALUATIONS for the native sampling path (fn_ps7_expected_sampling_gradient_evaluations_per_iteration,
          ## fn_ps7_extract_tau_sweep.R): expected leapfrog steps + 1 endpoint evaluation per iteration, except joint diffusion
          ## kick_flow_kick (endpoint reused), plus one start-of-call evaluation per chain on the joint diffusion paths.
          ## Path flags as the summary reported them (else the PS7 settings); these models always sample nuisance parameters:
          diffusion_HMC_sampling <-  if (!is.null(HMC_info$diffusion_HMC)) HMC_info$diffusion_HMC else settings$diffusion_HMC
          diffusion_HMC_integrator_sampling <-  if (!is.null(HMC_info$diffusion_HMC_integrator)) HMC_info$diffusion_HMC_integrator else settings$diffusion_HMC_integrator
          partitioned_HMC_sampling <-  if (!is.null(HMC_info$partitioned_HMC)) HMC_info$partitioned_HMC else settings$partitioned_HMC
          has_nuisance_sampling <-  Model_type %in% c("LC_MVP", "MVP", "LC_MVOP", "MVOP", "latent_trait")
          randomize_tau_sampling_used <-  if (is.null(settings$randomize_tau_sampling)) TRUE else settings$randomize_tau_sampling
          n_grad_evals_sampling_main <-  fn_ps7_expected_sampling_gradient_evaluations_per_iteration( tau = HMC_info$tau_main,
                                                                                                      eps = HMC_info$eps_main,
                                                                                                      randomize_tau_sampling = randomize_tau_sampling_used,
                                                                                                      diffusion_HMC = diffusion_HMC_sampling,
                                                                                                      diffusion_HMC_integrator = diffusion_HMC_integrator_sampling,
                                                                                                      partitioned_HMC = partitioned_HMC_sampling,
                                                                                                      has_nuisance = has_nuisance_sampling,
                                                                                                      n_iter_per_chain_call = settings$n_iter) *
                                        settings$n_iter * settings$n_chains_sampling
          Min_ess_per_grad_main_samp <-  min_ESS / n_grad_evals_sampling_main
          Min_ess_tail_per_grad_main_samp <-  min_ESS_tail / n_grad_evals_sampling_main   ## the same gradient count as the bulk value
          Min_ess_sd_per_grad_main_samp <-  min_ESS_sd / n_grad_evals_sampling_main       ## the same gradient count as the bulk value
          ##
          L_us_during_sampling <-  (HMC_info$tau_us / HMC_info$eps_us)
          ## n_grad_evals_sampling_us <-  L_us_during_sampling  * settings$n_iter * settings$n_chains_sampling
          ## nuisance gradient evaluations of a partitioned run: expected leapfrog steps at tau_us / eps_us + 1 per iteration
          ## (not used below: the weights give the nuisance gradient weight 0 for the joint, non-partitioned PS7 runs):
          n_grad_evals_sampling_us <-  fn_ps7_expected_sampling_gradient_evaluations_per_iteration( tau = HMC_info$tau_us,
                                                                                                    eps = HMC_info$eps_us,
                                                                                                    randomize_tau_sampling = randomize_tau_sampling_used,
                                                                                                    diffusion_HMC = diffusion_HMC_sampling,
                                                                                                    diffusion_HMC_integrator = diffusion_HMC_integrator_sampling,
                                                                                                    partitioned_HMC = TRUE,
                                                                                                    has_nuisance = has_nuisance_sampling) *
                                      settings$n_iter * settings$n_chains_sampling
          Min_ess_per_grad_us_samp <-  min_ESS / n_grad_evals_sampling_us
          ##
          ## (non-partitioned: "main" grad = "all" grad, so use main only)
          weight_nuisance_grad <-  0.00
          weight_main_grad <-  1.00
          ##
          Min_ess_per_grad_samp_weighted <-  Min_ess_per_grad_main_samp
          ##
          print(paste("ESS/sec = ", round(ESS_per_sec_samp, 3)))
          print(paste("ESS/grad = ", round(Min_ess_per_grad_samp_weighted*1000, 3)))
          ##
          # cat(sprintf("  min_ESS = %d, max_Rhat = %.4f\n", round(min_ESS), max_Rhat))
          print(paste("min_ESS = ", round(min_ESS), "max_Rhat =" , max_Rhat))
          ##
          # if (is.infinite(min_ESS)) { 
          #   beepr::beep()
          #   stop("Inf!!!")
          # }
          ##
          ## ---- Retain sampler evidence for rejection/sticking investigations; do not alter or drop any draws.
          ## Acceptance probabilities are iteration x chain matrices, not realised accept/reject indicators.
          ## Joint HMC uses p_jump_main for the whole joint proposal; nuisance has zero columns when not sampled.
          ## Older installed builds do not return this field. Keep that absence explicit rather than inventing values.
          ##
          sampling_diagnostics <-  model_results$sampling_object$sampling_diagnostics
          if (is.null(x = sampling_diagnostics)) {
            warning("Sampling acceptance traces are unavailable from this installed BayesMVP build; rebuild the package to capture them.",
                    call. = FALSE)
          } else {
            ## ---- post-burnin (sampling) acceptance: p_jump_main is n_iter x n_chains_sampling. If it sits well below
            ##      adapt_delta, eps came out of the burn-in too big; if it climbs from the start to the end, the chains
            ##      were still settling after the handover:
            try({
              p_jump_sampling_matrix <-  as.matrix(sampling_diagnostics$p_jump_main)
              n_sampling_iterations_recorded <-  nrow(p_jump_sampling_matrix)
              first_sampling_iterations <-  seq_len(min(10, n_sampling_iterations_recorded))
              second_half_sampling_iterations <-  seq(from = floor(n_sampling_iterations_recorded / 2) + 1, to = n_sampling_iterations_recorded)
              cat(paste0("post-burnin mean acceptance = ", sub("^ +", "", formatC(mean(p_jump_sampling_matrix, na.rm = TRUE), format = "f", digits = 3)), " (target adapt_delta = ", sub("^ +", "", formatC(settings$adapt_delta, format = "f", digits = 2)), ") | iterations 1-", formatC(as.integer(length(first_sampling_iterations)), format = "d"), ": ", sub("^ +", "", formatC(mean(p_jump_sampling_matrix[first_sampling_iterations, , drop = FALSE], na.rm = TRUE), format = "f", digits = 3)), " | second half: ", sub("^ +", "", formatC(mean(p_jump_sampling_matrix[second_half_sampling_iterations, , drop = FALSE], na.rm = TRUE), format = "f", digits = 3)), " | eps = ", sub("^ +", "", formatC(HMC_info$eps_main, format = "f", digits = 3)), ", L = ", sub("^ +", "", formatC(L_main_during_sampling, format = "f", digits = 1)), "\n"))
            })
          }
          ##
          ## ---- time to the target min ESS (ps3 targets; min ESS over the parameters of interest): the burn-in plus the
          ##      sampling and summary times scaled by target ESS / observed ESS; the burn-in time is NEVER scaled:
          ##
          try({
            target_min_ESS <-  R_fn_ps7_get_target_min_ESS( N          = N,
                                                           Model_type = Model_type)
            if (is.finite(target_min_ESS) && target_min_ESS > 0 && is.finite(min_ESS) && min_ESS > 0 &&
                is.finite(time_burnin) && time_burnin >= 0 && is.finite(time_sampling) && time_sampling >= 0 &&
                is.finite(time_summaries) && time_summaries >= 0) {
              sampling_time_to_target_min_ESS <-  time_sampling * target_min_ESS / min_ESS
              summary_time_to_target_min_ESS <-  time_summaries * target_min_ESS / min_ESS
              ## the same for the posterior SD's ESS (min_ESS_sd), a separate endpoint reported alongside:
              if (is.finite(min_ESS_sd) && min_ESS_sd > 0) {
                  time_to_target_min_ESS_sd <-  time_burnin + (time_sampling + time_summaries) * target_min_ESS / min_ESS_sd
                  cat(BayesMVP:::colourise(paste0("time to min ESS(theta^2) (same target; SD's ESS ", round(min_ESS_sd), "): ",
                                                  formatC(time_to_target_min_ESS_sd, format = "f", digits = 1), " s"), "cyan"), "\n")
              }
              cat(BayesMVP:::colourise(paste0("time to min ESS (target ", sub("^ +", "", formatC(target_min_ESS, format = "f", digits = 0)), " at N = ", sub("^ +", "", formatC(N, format = "f", digits = 0)), "): burn-in ", sub("^ +", "", formatC(time_burnin, format = "f", digits = 1)), " s + sampling ", sub("^ +", "", formatC(sampling_time_to_target_min_ESS, format = "f", digits = 1)), " s + summaries ", sub("^ +", "", formatC(summary_time_to_target_min_ESS, format = "f", digits = 1)), " s = ", sub("^ +", "", formatC(time_burnin + sampling_time_to_target_min_ESS + summary_time_to_target_min_ESS, format = "f", digits = 1)), " s | (this run: min ESS ", sub("^ +", "", formatC(min_ESS, format = "f", digits = 0)), " after ", sub("^ +", "", formatC(time_sampling, format = "f", digits = 1)), " s of sampling)"), "blue"), "\n")
            }
          })
          ##
          all_runs_results[[run_i]] <-  list(
            run_id = run_i,
            ##
            MCMC_seed = MCMC_seed,
            ##
            Se_mean = tibble_Se$mean,
            Sp_mean = tibble_Sp$mean,
            prev_mean = tibble_prev$mean,
            ##
            tibble_main = tibble_main,
            tibble_tp = tibble_tp,
            tibble_gq = tibble_gq,
            diagnostic_parameter_names = tibble_gq_for_diagnostics$parameter,
            ##
            trace_main = trace_main,
            ## trace_tp is NOT saved: it was ~76% of every saved run file (50 x 180 x 302 doubles, ~22 MB serialised,
            ## ~0.55 s of gzip per fit), nothing reads it (the readers use tibble_tp), and it can be rebuilt from
            ## trace_main when needed:
            # trace_tp = trace_tp,
            trace_gq = trace_gq,
            ##
            max_Rhat = max_Rhat,
            max_nRhat = max_nRhat,
            ##
            min_ESS = min_ESS,
            min_ESS_tail = min_ESS_tail,
            min_ESS_sd = min_ESS_sd,
            ##
            ESS_per_sec_total = ESS_per_sec_total,
            ESS_per_sec_samp = ESS_per_sec_samp,
            ESS_per_grad_samp = Min_ess_per_grad_samp_weighted,
            ESS_tail_per_grad_samp = Min_ess_tail_per_grad_main_samp,
            ESS_sd_per_grad_samp = Min_ess_sd_per_grad_main_samp,
            ## the gradient count behind ESS_per_grad_samp (readers rescale only runs without this marker, i.e. saved with the nominal tau / eps count):
            ## ESS_per_grad_samp_gradient_count = "executed_leapfrog_steps_expected",
            ## expected gradient evaluations for the native sampling path (fn_ps7_expected_sampling_gradient_evaluations_per_iteration):
            ESS_per_grad_samp_gradient_count = "gradient_evaluations_expected",
            ##
            n_chains = settings$n_chains_sampling,
            n_iter = settings$n_iter,
            ##
            summary_tibbles = summary_tibbles,
            efficiency_info = efficiency_info,
            divergences = divergences,
            HMC_info = HMC_info,
            ## the run's settings, with the pooled metric estimator settings as validated above (and passed to the sampler):
            ## settings = settings,
            settings = utils::modifyList(x = settings,
                                         val = list(metric_pooled_window_resets = metric_pooled_window_resets_requested,
                                                    metric_pooled_offdiagonal_shrinkage = metric_pooled_offdiagonal_shrinkage_requested)),
            burnin_algorithm = burnin_algorithm_used,
            trajectory_adaptation_version = "metric_main_v4",
            vect_type = vect_type,
            trajectory_criterion = model_results$trajectory_criterion,
            adaptation = model_fit_object$adaptation,
            package_versions = c(NicoStan = as.character(utils::packageVersion("NicoStan")),
                                 BayesMVP = as.character(utils::packageVersion("BayesMVP"))),
            store_log_lik_trace = settings$store_log_lik_trace,
            theta_hat_us_freeze_iter = settings$theta_hat_us_freeze_iter,
            manual_tau_value = if (is.null(settings$manual_tau_value)) NA_real_ else settings$manual_tau_value,
            J_grad_option = J_grad_option_used,
            autodiff_fallback = settings$autodiff_fallback,
            ##
            ## trajectory jitter / tau rescaling / ADAM bias correction, as the sampler reported them (all checked above):
            randomize_tau_burnin = model_results$randomize_tau_burnin,
            tau_sampling_scale = model_results$tau_sampling_scale,
            tau_sampling_scale_effective = model_results$tau_sampling_scale_effective,
            tau_sampling_scale_factor = model_results$tau_sampling_scale_factor,
            tau_main_before_sampling_scale = model_results$tau_main_before_sampling_scale,
            tau_main_after_sampling_scale = model_results$tau_main_after_sampling_scale,
            tau_adam_bias_correction = tau_adam_bias_correction_reported,
            ## step-size acceptance mean and divergence-triggered tau shrink, as the sampler reported them (checked above);
            ## absent in runs saved before these options existed, which used "arithmetic" / TRUE:
            eps_acceptance_mean = eps_acceptance_mean_reported,
            tau_shrink_on_divergence = as.logical(tau_shrink_on_divergence_reported),
            ## pooled metric estimator ("_Mp" code), as the sampler reported them (checked above): the window resets as requested
            ## ("none" / "stan_style" / iterations), the reset iterations the burn-in used, and the off-diagonal shrinkage. Absent
            ## in runs saved before these options existed, whose "_Mp" runs used "stan_style" / 0. Recorded for every estimator,
            ## but only "pooled" reads them:
            metric_pooled_window_resets = metric_pooled_window_resets_reported,
            metric_pooled_window_reset_iterations = metric_pooled_window_reset_iterations_reported,
            metric_pooled_offdiagonal_shrinkage = if (identical(metric_pooled_offdiagonal_shrinkage_reported, "adaptive")) "adaptive" else
                                                    as.numeric(metric_pooled_offdiagonal_shrinkage_reported),
            metric_pooled_offdiagonal_shrinkage_mode = if (identical(metric_pooled_offdiagonal_shrinkage_reported, "adaptive")) "adaptive" else "fixed",
            ## trajectory-length criterion block ("_tbJ";, EXPERIMENTAL), as the sampler reported it (checked above):
            tau_adaptation_block = tau_adaptation_block_reported,
            ## ADAM ("_A" code; also in settings above): the settings of the eps and the tau updates and the tau learning-rate restart at
            ## the metric freeze, as the burn-in reported them (checked above), and the restart iteration:
            eps_adam_beta1 = adam_settings_reported$eps_adam_beta1,
            eps_adam_beta2 = adam_settings_reported$eps_adam_beta2,
            eps_adam_epsilon = adam_settings_reported$eps_adam_epsilon,
            tau_adam_beta1 = adam_settings_reported$tau_adam_beta1,
            tau_adam_beta2 = adam_settings_reported$tau_adam_beta2,
            tau_adam_epsilon = adam_settings_reported$tau_adam_epsilon,
            tau_learning_rate_restart_at_metric_end = tau_learning_rate_restart_reported,
            tau_learning_rate_restart_iteration = model_results$burnin_object$tau_learning_rate_restart_iteration,
            ## tau adaptation scheme ("p" code) and, for "probe_then_average", the probe end iteration, the probe tau path and the averaged tau:
            tau_adaptation_scheme = tau_adaptation_scheme_reported,
            tau_probe_end_iteration = model_results$burnin_object$tau_probe_end_iteration,
            tau_probe_tau_path_vec = model_results$burnin_object$tau_probe_tau_path_vec,
            tau_main_averaged = model_results$burnin_object$tau_main_averaged,
            tau_initial_effective = model_results$tau_initial,
            ## tau_initial "adaptive": lambda_max estimates and the tau set at the handover (NULL for a numeric tau_initial):
            tau_initial_adaptive = model_results$burnin_object$tau_initial_adaptive,
            learning_rate_initial_effective = model_results$learning_rate_initial,
            ##
            sampler_diagnostics = list(
                burnin_profile = model_results$burnin_profile,
                burnin_schedule = model_results$burnin_schedule,
                sampling = sampling_diagnostics,
                divergence_traces = model_results$sampling_object[[2]],  ## list of 1 x iteration matrices, one per chain
                burnin_metric = model_results$burnin_object$EHMC_Metric_as_Rcpp_List,
                metric_adaptive_shrinkage = model_results$burnin_object$metric_adaptive_shrinkage,
                metric_adaptive_shrinkage_history = model_results$burnin_object$metric_adaptive_shrinkage_history,
                ## burn-in record for metric diagnosis: burn-in chains' main draws (iteration x parameter x chain), main metric
                ## variances used at each iteration, nuisance metric variance quantiles at each iteration (NULL for runs before it existed):
                burnin_trace_main_all_chains = model_results$burnin_object$burnin_trace_main_all_chains,
                burnin_metric_main_variance_history = model_results$burnin_object$burnin_metric_main_variance_history,
                burnin_metric_nuisance_variance_quantiles_history = model_results$burnin_object$burnin_metric_nuisance_variance_quantiles_history,
                burnin_integration_args = model_results$burnin_object$EHMC_args_as_Rcpp_List),
            ##
            time_total = time_total,
            time_burnin = time_burnin
          )
          ##
          ## CHESSR_time / SNAPER_time only ("_bact" / "_bast"): what PS7 handed to the sampler (time_criterion_inputs: the effective
          ## settings, the settings text and hash of the file name, the time_criterion_settings list and the expected source of each
          ## sampling time), and the sampler's time-criterion record as returned (time_criterion, checked above: the probe settings,
          ## results and wall time, the sampling quantities and their sources, the online burn-in leapfrog time, the iteration counts,
          ## and tau_offset_from_sampling_overhead / burnin_to_sampling_leapfrog_time_ratio at the handover and at the end). Added only
          ## for these two algorithms, so the saved run of every other algorithm keeps exactly its earlier structure:
          if (time_criterion_is_used) {
            all_runs_results[[run_i]]$time_criterion_inputs <-  time_criterion_inputs
            all_runs_results[[run_i]]$time_criterion <-  time_criterion_record
          }
          ##
          # gc(verbose = FALSE)   ## removed: a full collection after every fit cost several seconds of wall time per fit
          ##
          ps7_step_times["extraction"] <-  fn_ps7_stamp()  ## tibbles, diagnostics, cov/eigen, prints and the run list
          ##
          ## ---- Save raw results (SAME builder as the resume check above => names always match):
          ##
          file_name_string <-  R_fn_file_name_string( Model_type = Model_type,
                                                     N = N,
                                                     settings = settings,
                                                     model_args_list = model_args_list,
                                                     prior_LKJ_nd = prior_LKJ_nd,
                                                     prior_LKJ_d = prior_LKJ_d,
                                                     prior_prev_a = prior_prev_a,
                                                     prior_prev_b = prior_prev_b,
                                                     grouping = grouping)
          ##
          if (!identical(file_name_string, file_name_string_check)) {
            warning(paste0("run_ps7_models: save-name != check-name?!\n  save : ", file_name_string,
                           "\n  check: ", file_name_string_check))
          }
          ##
          run_file <-  paste0(file_name_string, "_run", run_i)
          ##
          ps7_step_times["file_name"] <-  fn_ps7_stamp()  ## R_fn_file_name_string() for the save
          ##
          saveRDS(all_runs_results[[run_i]], run_file)
          ##
          ps7_step_times["save"] <-  fn_ps7_stamp()  ## saveRDS() (default gzip compression)
          ##
          ## ---- One-line wall-clock breakdown of this fit (printed only; NOT stored in the saved run, because the save time is
          ##      only known after saveRDS()). time_total is the sampler's sum of its timed sections (burn-in + sampling +
          ##      summaries), so "untimed" = wall time of the R_fn_sample_model() call + create_summary_and_traces() minus
          ##      time_total, i.e. the work inside those two calls that no sampler timer covers:
          ##
          if (profile_wrapper_steps) {
            try({
              ps7_step_seconds <-  diff(ps7_step_times)  ## named by the step that ENDS at each stamp
              ps7_fit_wall_seconds <-  ps7_step_times[["save"]] - ps7_step_times[["fit_start"]]
              ps7_time_total <-  if (is.null(time_total)) NA_real_ else as.numeric(time_total)[1]
              ps7_untimed_seconds <-  ps7_step_seconds[["sampler_call"]] + ps7_step_seconds[["summaries"]] - ps7_time_total
              ps7_other_seconds <-  ps7_fit_wall_seconds - sum(ps7_step_seconds)
              ##
              message(BayesMVP:::colourise(paste0("  [ps7 wrapper timing] run ", run_i,
                                                  ": fit wall ", fn_ps7_format_seconds(ps7_fit_wall_seconds), " s",
                                                  " | sampler time_total ", fn_ps7_format_seconds(ps7_time_total), " s",
                                                  " (burn-in ", fn_ps7_format_seconds(time_burnin),
                                                  " + sampling ", fn_ps7_format_seconds(time_sampling),
                                                  " + summaries ", fn_ps7_format_seconds(time_summaries), ")",
                                                  " | untimed in sampler call + summaries ", fn_ps7_format_seconds(ps7_untimed_seconds), " s",
                                                  " | init ", fn_ps7_format_seconds(ps7_step_seconds[["init"]]),
                                                  " | prep ", fn_ps7_format_seconds(ps7_step_seconds[["prep"]]),
                                                  " | sampler init ", fn_ps7_format_seconds(ps7_step_seconds[["sampler_init"]]),
                                                  " | sampler call ", fn_ps7_format_seconds(ps7_step_seconds[["sampler_call"]]),
                                                  " | post-fit checks ", fn_ps7_format_seconds(ps7_step_seconds[["post_fit_checks"]]),
                                                  " | summaries ", fn_ps7_format_seconds(ps7_step_seconds[["summaries"]]),
                                                  " | extraction ", fn_ps7_format_seconds(ps7_step_seconds[["extraction"]]),
                                                  " | file name ", fn_ps7_format_seconds(ps7_step_seconds[["file_name"]]),
                                                  " | save ", fn_ps7_format_seconds(ps7_step_seconds[["save"]]),
                                                  " | other ", fn_ps7_format_seconds(ps7_other_seconds),
                                                  " | name+resume (once per call) ", fn_ps7_format_seconds(ps7_resume_time_seconds)), "cyan"))
            })
          }
          
        }
        
        return(all_runs_results)
  
}
##
## ---- 4. Function to summarize BayesMVP pilot study 7 results (FIXED loader filename): -----------------------------------------------------
##
summarize_ps7_results <-  function(output_dir,
                                  ##
                                  Model_type,
                                  ##
                                  N_val,
                                  ##
                                  target_ESS_val,
                                  ##
                                  n_pops = ifelse(Model_type == "LC_MVOP", 3, 2),
                                  ##
                                  n_chains_burnin = 8,
                                  n_chains_sampling = 64,
                                  ##
                                  n_iter = 1000,
                                  ##
                                  adapt_delta = 0.80,
                                  ##
                                  n_burnin_vec = c(250, 500, 1000),
                                  clip_iter_vec = c(50, 100, 100),
                                  int_vec = c(75, 100, 300),
                                  ##
                                  ratio_M_main_vec,
                                  ratio_M_nuisance_vec,
                                  ##
                                  int_width_vec = c(1, 1, 1),
                                  ##
                                  learning_rate_vec = c(0.0125, 0.025, 0.0375, 0.05, 0.0625, 0.075),
                                  ##
                                  metric_type_main,
                                  metric_shape_main,
                                  ##
                                  metric_type_nuisance, ## 
                                  ##
                                  diffusion_HMC,
                                  partitioned_HMC,
                                  ##
                                  ## must MIRROR the run settings, because the file name is rebuilt here to FIND
                                  ## the saved runs - a stub field that is missing or wrong means no file matches
                                  ## and every run looks un-done (the same trap as the old missing "_n_pops"):
                                  diffusion_HMC_integrator = "kick_flow_kick",
                                  ##
                                  learning_rate_initial = 0.15,
                                  learning_rate_initial_iter = NULL,
                                  ##
                                  eps_initial = NULL,
                                  eps_initial_iter = NULL,
                                  ##
                                  M_decay_type,
                                  M_decay_power = NULL,
                                  ##
                                  n_runs = 10,
                                  ##
                                  prior_LKJ_nd,
                                  prior_LKJ_d,
                                  prior_prev_a,
                                  prior_prev_b,
                                  ##
                                  grouping,
                                  ##
                                  n_threads_WCP_burnin,
                                  n_threads_WCP_sampling,
                                  ##
                                  num_chunks_burnin,
                                  num_chunks_sampling,
                                  ##
                                  reorder_cols_MVP,
                                  multi_attempts,
                                  ##
                                  metric_estimator,
                                  tau_initial = c(pi, 2*pi),
                                  ##
                                  ## ---- Trajectory-length adaptation axes. Like every other argument here these must
                                  ##      MIRROR the run settings, since the file name is rebuilt to FIND the saved runs.
                                  ##      The defaults are the original ChESSR configuration, which encodes to nothing, so
                                  ##      calling this function exactly as before still resolves every existing file.
                                  ##
                                  ## "KE" = the kinetic-energy-change statistic; "ChEES" = position-based (Hoffman 2021).
                                  ## Pass c("KE", "ChEES", "CHESSR", "CHESSR_log", "SNAPER") to compare all three in one data frame.
                                  ## "CHESSR_time" / "SNAPER_time" ("_bact" / "_bast") are accepted too, at their default time-criterion
                                  ## settings (this reader has no time-criterion arguments, so it never rebuilds a "_bac<4 hex>" /
                                  ## "_bas<4 hex>" name; fn_ps7_extract_tau_sweep() finds every such run from the names on disk).
                                  ## "ESJD" / "ESJD_CHESSR" ("_baej" / "_baex") are accepted as well.
                                  ## "ESJD_SNAPER" ("_baen") is accepted as well.
                                  burnin_algorithm = "KE",
                                  ##
                                  ## FALSE = median across chains with the accept indicator (original);
                                  ## TRUE = acceptance-weighted mean using p_jump.
                                  tau_weight_by_p_jump = FALSE,
                                  ##
                                  ## Pinned leapfrog-step counts, for the fixed-L sweep. NA = tau adapted (original).
                                  ## e.g. manual_L = c(NA, 4, 6, 9, 14, 23, 40).
                                  manual_L = NA_real_,
                                  manual_tau_value = NA_real_,
                                  ##
                                  ## Which objectives get tau_weight_by_p_jump = TRUE swept for them.
                                  ## NULL = all of them (the plain cross). MUST mirror the runner.
                                  tau_weight_by_p_jump_for = NULL,
                                  ##
                                  ## Burn-in tau ramp: "original", "staged", or c("original", "staged") to compare both.
                                  ## NULL = runs saved before the option existed (no "_tr" token in their file names).
                                  tau_ramp = NULL,
                                  ##
                                  ## The fixed test order the runs used, e.g. c(5, 4, 6, 1, 3, 2) ("_tp546132" in the
                                  ## file name). NULL = the order was estimated by the pre-burnin (no "_tp" token).
                                  test_perm_override = NULL,
                                  ##
                                  ## TRUE = eps re-initialised at the ChEES handover (no token; every run before the option existed),
                                  ## FALSE = kept ("_nER"), c(TRUE, FALSE) = compare both.
                                  eps_reinit_at_ChEES_handover = TRUE,
                                  ##
                                  ## "running_mean_frozen" ("_cF"), "zero" ("_c0") or "running_mean" (no token: every run saved
                                  ## previously). A vector compares several.
                                  theta_hat_us_rule = "running_mean_frozen",
                                  ##
                                  ## burn-in speed-ups, as set in the runner ("_pa<k>" / "_pb<n>"); NULL = not used
                                  ## Historical _pa reference above is retired; this interface only exposes pre-burn-in controls.
                                  pre_burnin_n_iter = NULL,
                                  pre_burnin_L = NULL,
                                  share_tau_ii_across_chains_in_burnin = FALSE,
                                  burnin_TBB_pool_equals_n_chains = FALSE,
                                  ##
                                  ## "num_diff" (no token, Earlier runs) or "autodiff" ("_JgA"):
                                  J_grad_option = "num_diff",
                                  ## NULL = automatic schedule; a vector (with optional NA for auto) mirrors the runner.
                                  theta_hat_us_freeze_iter = NULL,
                                  burnin_schedule = NULL,
                                  metric_adaptation_end_iter = NULL,
                                  ## FALSE = current PS7 default; TRUE = "_ll1". A vector compares both.
                                  store_log_lik_trace = FALSE,
                                  ## FALSE keeps existing filenames; TRUE finds "_af1" runs, c(FALSE, TRUE) reads both.
                                  ## autodiff_fallback = FALSE
                                  autodiff_fallback = FALSE,
                                  ## pooled metric estimator (the "_Mp" code): NULL = "stan_style" / 0, which finds every run saved
                                  ## before these options existed (plain "_Mp"); "none" with 0.5 finds the "_Mp.5w" runs:
                                  metric_pooled_window_resets = NULL,
                                  metric_pooled_offdiagonal_shrinkage = NULL

) {
  
        require(dplyr)
        require(tibble)
        
        ## Storage for all results
        all_rows <-  list()
        row_idx <-  1
        loaded_configurations <-  character()
        ##
        ## Each of these arguments accepts one value or a vector, just like the runner's top-level controls.
        sampler_combinations <-  expand.grid( metric_estimator = metric_estimator,
                                              tau_initial = tau_initial,
                                              learning_rate_initial = learning_rate_initial,
                                              diffusion_HMC_integrator = diffusion_HMC_integrator,
                                              burnin_algorithm = burnin_algorithm,
                                              tau_weight_by_p_jump = tau_weight_by_p_jump,
                                              manual_L = manual_L,
                                              manual_tau_value = manual_tau_value,
                                              tau_ramp = if (is.null(tau_ramp)) NA_character_ else tau_ramp,
                                              eps_reinit_at_ChEES_handover = eps_reinit_at_ChEES_handover,
                                              theta_hat_us_rule = theta_hat_us_rule,
                                              theta_hat_us_freeze_iter = if (is.null(theta_hat_us_freeze_iter))
                                                  NA_real_ else theta_hat_us_freeze_iter,
                                              store_log_lik_trace = store_log_lik_trace,
                                              autodiff_fallback = autodiff_fallback,
                                              KEEP.OUT.ATTRS = FALSE,
                                              stringsAsFactors = FALSE)
        ## A pinned trajectory length switches the adaptation off, which makes burnin_algorithm,
        ## tau_weight_by_p_jump AND tau_initial inert. The shared helper collapses them so the reader
        ## rebuilds exactly the names the runner saved:
        sampler_combinations <-  fn_ps7_collapse_inert_tau_axes(sampler_combinations = sampler_combinations,
                                                               weight_p_jump_only_for = tau_weight_by_p_jump_for)
        ##
        if (nrow(x = sampler_combinations) == 0) stop("Select at least one value for each of the sampler options.")
        ## if (any(!sampler_combinations$metric_estimator %in% c("pooled", "chain_mean", "chain_mean_scaled"))) {
        if (any(!sampler_combinations$metric_estimator %in% c("pooled", "chain_mean", "chain_mean_scaled", "per_iteration"))) {
            ## stop("metric_estimator must contain only 'pooled', 'chain_mean' and/or 'chain_mean_scaled'.")
            stop("metric_estimator must contain only 'pooled', 'chain_mean', 'chain_mean_scaled' and/or 'per_iteration'.")
        }
        if (any(!sampler_combinations$diffusion_HMC_integrator %in% c("kick_flow_kick", "flow_kick_flow"))) {
            stop("diffusion_HMC_integrator must contain only 'kick_flow_kick' and/or 'flow_kick_flow'.")
        }
        ## if (any(!sampler_combinations$burnin_algorithm %in% c("KE", "ChEES", "CHESSR", "CHESSR_log", "SNAPER"))) {
            ## stop("burnin_algorithm must contain only 'KE', 'ChEES', 'CHESSR', 'CHESSR_log' or 'SNAPER'.")
        ## if (any(!sampler_combinations$burnin_algorithm %in% c("KE", "ChEES", "CHESSR", "CHESSR_log", "SNAPER", "CHESSR_time", "SNAPER_time"))) {
            ## stop("burnin_algorithm must contain only 'KE', 'ChEES', 'CHESSR', 'CHESSR_log', 'SNAPER', 'CHESSR_time' or 'SNAPER_time'.")
        if (any(!sampler_combinations$burnin_algorithm %in% c("KE", "ChEES", "CHESSR", "CHESSR_log", "SNAPER", "CHESSR_time", "SNAPER_time",
                                                              "ESJD", "ESJD_CHESSR", "ESJD_SNAPER", "LQ_ESSR"))) {
            ## stop("burnin_algorithm must contain only 'KE', 'ChEES', 'CHESSR', 'CHESSR_log', 'SNAPER', 'CHESSR_time', 'SNAPER_time', 'ESJD', 'ESJD_CHESSR' or 'ESJD_SNAPER'.")
            stop(paste0("burnin_algorithm must contain only 'KE', 'ChEES', 'CHESSR', 'CHESSR_log', 'SNAPER', 'CHESSR_time', 'SNAPER_time', 'ESJD', ",
                        "'ESJD_CHESSR', 'ESJD_SNAPER' or 'LQ_ESSR'."))
        }
        if (any(!is.na(sampler_combinations$tau_ramp) & !(sampler_combinations$tau_ramp %in% c("original", "staged")))) {
            stop("tau_ramp must be NULL, or contain only 'original' and/or 'staged'.")
        }
        if (!is.logical(sampler_combinations$eps_reinit_at_ChEES_handover) || any(is.na(sampler_combinations$eps_reinit_at_ChEES_handover))) {
            stop("eps_reinit_at_ChEES_handover must contain only TRUE and/or FALSE.")
        }
        if (any(!is.logical(sampler_combinations$tau_weight_by_p_jump)) ||
            any(is.na(sampler_combinations$tau_weight_by_p_jump))) {
            stop("tau_weight_by_p_jump must contain only TRUE and/or FALSE.")
        }
        if (any(!is.na(sampler_combinations$manual_L) &
                (!is.finite(sampler_combinations$manual_L) | sampler_combinations$manual_L < 1))) {
            stop("manual_L must contain NA (adapt tau) and/or leapfrog-step counts >= 1.")
        }
        if ("adaptive" %in% tau_initial) {
            stop("summarize_ps7_results() reads numeric tau_initial levels only; for tau_initial = \"adaptive\" use fn_ps7_extract_tau_sweep() and fn_ps7_summarise_runs().")
        }
        if (!is.numeric(tau_initial) || any(!is.finite(tau_initial) | tau_initial <= 0)) {
            stop("tau_initial must contain positive finite numbers.")
        }
        if (length(target_ESS_val) != 1 || !is.numeric(target_ESS_val) || !is.finite(target_ESS_val) || target_ESS_val <= 0) {
            stop("target_ESS_val must be a single positive finite number.")
        }
        if (any(lengths(list(clip_iter_vec, int_vec, ratio_M_main_vec, ratio_M_nuisance_vec, int_width_vec)) != length(n_burnin_vec))) {
            stop("Burn-in configuration vectors must have the same length as n_burnin_vec; they are paired by position.")
        }
        ##
        for (combination_index in seq_len(length.out = nrow(x = sampler_combinations))) {
          metric_estimator <-  sampler_combinations$metric_estimator[combination_index]
          tau_initial <-  sampler_combinations$tau_initial[combination_index]
          learning_rate_initial <-  sampler_combinations$learning_rate_initial[combination_index]
          diffusion_HMC_integrator <-  sampler_combinations$diffusion_HMC_integrator[combination_index]
          burnin_algorithm <-  sampler_combinations$burnin_algorithm[combination_index]
          tau_weight_by_p_jump <-  sampler_combinations$tau_weight_by_p_jump[combination_index]
          manual_L <-  sampler_combinations$manual_L[combination_index]
          manual_tau_value <-  sampler_combinations$manual_tau_value[combination_index]
          tau_ramp <-  sampler_combinations$tau_ramp[combination_index]
          eps_reinit_at_ChEES_handover <-  sampler_combinations$eps_reinit_at_ChEES_handover[combination_index]
          theta_hat_us_rule <-  sampler_combinations$theta_hat_us_rule[combination_index]
          theta_hat_us_freeze_iter <-  sampler_combinations$theta_hat_us_freeze_iter[combination_index]
          if (is.na(theta_hat_us_freeze_iter)) theta_hat_us_freeze_iter <-  NULL
          store_log_lik_trace <-  sampler_combinations$store_log_lik_trace[combination_index]

        ## Loop through all configurations using INDEX
        for (ii in seq_along(n_burnin_vec)) {
          
          n_burnin <-  n_burnin_vec[ii]
          clip_iter <-  clip_iter_vec[ii]
          int <-  int_vec[ii]
          ratio_M_main <-  ratio_M_main_vec[ii]
          ratio_M_nuisance <-  ratio_M_nuisance_vec[ii]
          int_width <-  int_width_vec[ii]
          
          if (n_burnin == 1000) {
            n_adapt <-  900
          } else {
            n_adapt <-  n_burnin - round(n_burnin/10)
          }
          ##
          # M_decay_scale <-  n_adapt / 5 ## -----------------
          M_decay_scale <-  n_adapt / 100 ## -----------------
          
          try({
            
            for (learning_rate in learning_rate_vec) {
              
              ## ---- Build filename via the SHARED builder (FIX: old inline copy here was
              ##      missing "_n_pops" so it NEVER matched the saved files):
              ##
              settings_stub <-  list( output_dir = output_dir,
                                     ##
                                     tau_initial = tau_initial,
                                     ##
                                     diffusion_HMC = diffusion_HMC,
                                     partitioned_HMC = partitioned_HMC,
                                     diffusion_HMC_integrator = diffusion_HMC_integrator,
                                     ##
                                     ## must mirror the runner, or the rebuilt name matches nothing:
                                     burnin_algorithm = burnin_algorithm,
                                     tau_weight_by_p_jump = tau_weight_by_p_jump,
                                     manual_L = manual_L,
                                     manual_tau_value = manual_tau_value,
                                     ##
                                     learning_rate_initial = learning_rate_initial,
                                     learning_rate_initial_iter = learning_rate_initial_iter,
                                     eps_initial = eps_initial,
                                     eps_initial_iter = eps_initial_iter,
                                     ##
                                     n_chains_burnin = n_chains_burnin,
                                     n_chains_sampling = n_chains_sampling,
                                     ##
                                     n_threads_WCP_burnin = n_threads_WCP_burnin,
                                     n_threads_WCP_sampling = n_threads_WCP_sampling,
                                     ##
                                     n_burnin = n_burnin,
                                     n_iter = n_iter,
                                     learning_rate = learning_rate,
                                     adapt_delta = adapt_delta,
                                     clip_iter = clip_iter,
                                     int = int,
                                     ratio_M_main = ratio_M_main,
                                     ratio_M_nuisance = ratio_M_nuisance,
                                     metric_type_main = metric_type_main,
                                     metric_type_nuisance = metric_type_nuisance, ## 
                                     metric_shape_main = metric_shape_main,
                                     int_width = int_width,
                                     ##
                                     M_decay_type = M_decay_type,
                                     M_decay_power = M_decay_power,
                                     M_decay_scale = M_decay_scale,
                                     ##
                                     reorder_cols_MVP = reorder_cols_MVP,
                                     multi_attempts = multi_attempts,
                                     ##
                                     metric_estimator = metric_estimator,
                                     ##
                                     num_chunks_burnin = num_chunks_burnin, ##
                                     num_chunks_sampling = num_chunks_sampling, ##
                                     ##
                                     ## "_tr" / "_tp" tokens - must mirror the runner, or runs saved with them are never found:
                                     tau_ramp = if (is.na(tau_ramp)) NULL else tau_ramp,
                                     test_perm_override = test_perm_override,
                                     eps_reinit_at_ChEES_handover = eps_reinit_at_ChEES_handover,
                                     theta_hat_us_rule = theta_hat_us_rule,
                                     theta_hat_us_freeze_iter = theta_hat_us_freeze_iter,
                                     burnin_schedule = burnin_schedule,
                                     metric_adaptation_end_iter = metric_adaptation_end_iter,
                                     pre_burnin_n_iter = pre_burnin_n_iter,
                                     pre_burnin_L = pre_burnin_L,
                                     share_tau_ii_across_chains_in_burnin = share_tau_ii_across_chains_in_burnin,
                                     burnin_TBB_pool_equals_n_chains = burnin_TBB_pool_equals_n_chains,
                                     J_grad_option = J_grad_option,
                                     store_log_lik_trace = store_log_lik_trace,
                                     ## autodiff_fallback = sampler_combinations$autodiff_fallback[combination_index])
                                     autodiff_fallback = sampler_combinations$autodiff_fallback[combination_index],
                                     metric_pooled_window_resets = metric_pooled_window_resets,
                                     metric_pooled_offdiagonal_shrinkage = metric_pooled_offdiagonal_shrinkage)
              ##
              model_args_stub <-  list( n_pops = n_pops,
                                       num_chunks = NULL)
              ##
              file_name_string <-  R_fn_file_name_string( Model_type = Model_type,
                                                         N = N_val,
                                                         settings = settings_stub,
                                                         model_args_list = model_args_stub,
                                                         prior_LKJ_nd = prior_LKJ_nd,
                                                         prior_LKJ_d = prior_LKJ_d,
                                                         prior_prev_a = prior_prev_a,
                                                         prior_prev_b = prior_prev_b,
                                                         grouping = grouping)
              ##
              print(paste(file_name_string))
              ##
              ## Repeated entries in the paired burn-in vectors must not count the same saved runs twice.
              if (file_name_string %in% loaded_configurations) next
              loaded_configurations <-  c(loaded_configurations, file_name_string)
              ##
              if (file.exists(file_name_string)) {
                all_runs <-  readRDS(file_name_string)
              } else {
                ## Load individual run files
                all_runs <-  list()
                for (run_i in 1:n_runs) {
                  run_file <-  paste0(file_name_string, "_run", run_i)
                  if (file.exists(run_file)) {
                    all_runs[[run_i]] <-  readRDS(run_file)
                  }
                }
              }
              
              if (length(all_runs) == 0) next
              
              for (run_i in seq_along(all_runs)) {
                
                try({
                  
                  r <-  all_runs[[run_i]]
                  if (is.null(r)) next
                  
                  eff  <-  r$efficiency_info
                  hmc  <-  r$HMC_info
                  divs <-  r$divergences
                  
                  row <-  tibble(
                    N = N_val,
                    metric_estimator = metric_estimator,
                    tau_initial = tau_initial,
                    learning_rate_initial = learning_rate_initial,
                    diffusion_HMC_integrator = diffusion_HMC_integrator,
                    ##
                    ## ---- Trajectory-length adaptation. Rows loaded from files saved before these
                    ##      fields existed get the original ChESSR settings, which is what they ran:
                    ##      the kinetic-energy objective, median across chains, tau adapted.
                    burnin_algorithm = burnin_algorithm,
                    tau_weight_by_p_jump = tau_weight_by_p_jump,
                    manual_L = manual_L,
                    manual_tau_value = manual_tau_value,
                    ##
                    ## NA = saved before these options existed (staged ramp from, original before;
                    ## test order estimated by the pre-burnin):
                    tau_ramp = tau_ramp,
                    test_perm_override = if (is.null(test_perm_override)) NA_character_ else paste(test_perm_override, collapse = "-"),
                    eps_reinit_at_ChEES_handover = eps_reinit_at_ChEES_handover,
                    theta_hat_us_rule = theta_hat_us_rule,
                    theta_hat_us_freeze_iter = if (is.null(theta_hat_us_freeze_iter)) NA_real_ else theta_hat_us_freeze_iter,
                    store_log_lik_trace = store_log_lik_trace,
                    ##
                    ## One readable label for the whole adaptation scheme, for grouping and plots:
                    tau_scheme = if (!is.na(manual_tau_value)) paste0("fixed_tau_", manual_tau_value)
                                 else if (!is.na(manual_L)) paste0("fixed_L_", manual_L)
                                 else paste0(burnin_algorithm, if (isTRUE(tau_weight_by_p_jump)) "_pjump" else ""),
                    learning_rate_initial_iter = if (is.null(learning_rate_initial_iter)) NA_real_ else learning_rate_initial_iter,
                    target_ESS = target_ESS_val,
                    n_burnin = n_burnin,
                    clip_iter = clip_iter,
                    int = int,
                    ##
                    ratio_M_main = ratio_M_main,
                    ratio_M_nuisance = ratio_M_nuisance,
                    ##
                    int_width = int_width,
                    learning_rate = learning_rate,
                    run = run_i,
                    n_chains_burnin = hmc$n_chains_burnin,
                    n_chains_sampling = hmc$n_chains_sampling,
                    n_iter = eff$n_iter,
                    ##
                    time_total = eff$time_total,
                    time_sampling = eff$time_sampling,
                    time_burnin = eff$time_burnin,
                    time_summaries = eff$time_summaries,
                    time_total_wo_summaries = eff$time_total_wo_summaries,
                    ##
                    max_Rhat = eff$Max_rhat_main,
                    max_nRhat = eff$Max_nested_rhat_main,
                    ##
                    min_ESS = r$min_ESS, ## eff$Min_ESS_main,
                    ## the posterior SD's ESS (NicoStan n_eff_sd; NA for runs saved before 4 Oct 2026) and its rate per 1000 gradients
                    ## with the SAME gradient count as ESS_per_grad_samp, plus its time to the target ESS (burn-in never scaled):
                    min_ESS_sd = if (is.null(r$min_ESS_sd)) NA_real_ else r$min_ESS_sd,
                    ESS_sd_per_grad_samp = if (is.null(r$min_ESS_sd) || !is.finite(r$min_ESS_sd) || !is.finite(r$min_ESS) || r$min_ESS <= 0) NA_real_ else
                                               fn_ps7_ESS_per_grad_samp_expected_steps(run_object = r) * 1000 * r$min_ESS_sd / r$min_ESS,
                    est_time_to_target_ESS_sd_mins = if (!is.null(r$min_ESS_sd) && is.finite(r$min_ESS_sd) && r$min_ESS_sd > 0)
                        (time_burnin + (target_ESS_val/r$min_ESS_sd)*(time_sampling + time_summaries))/60 else NA_real_,
                    ##
                    ESS_per_sec_total = r$ESS_per_sec_total, ## eff$ESS_per_sec_total,
                    ESS_per_sec_samp = r$ESS_per_sec_samp, ## eff$ESS_per_sec_samp,
                    ##
                    n_divs = divs$n_divs,
                    pct_divs = divs$pct_divs,
                    ##
                    ## ESS_per_grad_samp = r$ESS_per_grad_samp*1000,  ## legacy display column: min ESS per 1000 gradients
                    ## min_ESS_per_grad_sampling = r$ESS_per_grad_samp,  ## unscaled: min ESS per ONE gradient
                    ## per EXECUTED gradient (expected whole leapfrog steps under the sampling jitter; see fn_ps7_extract_tau_sweep.R):
                    ESS_per_grad_samp = fn_ps7_ESS_per_grad_samp_expected_steps(run_object = r)*1000,  ## legacy display column: min ESS per 1000 gradients
                    min_ESS_per_grad_sampling = fn_ps7_ESS_per_grad_samp_expected_steps(run_object = r),  ## unscaled: min ESS per ONE gradient
                    min_ESS_per_sec_sampling = r$ESS_per_sec_samp,
                    grad_evals_per_sec = eff$grad_evals_per_sec,
                    ##
                    L_main_burnin = eff$L_main_during_burnin,
                    L_main_samp = eff$L_main_during_sampling,
                    L_us_burnin = eff$L_us_during_burnin,
                    L_us_samp = eff$L_us_during_sampling,
                    ##
                    tau_main = hmc$tau_main,
                    eps_main = hmc$eps_main,
                    tau_us = hmc$tau_us,
                    eps_us = hmc$eps_us,
                    ##
                    adapt_delta = hmc$adapt_delta,
                    ##
                    sampling_time_to_Min_ESS = eff$sampling_time_to_Min_ESS,
                    total_time_to_1000_ESS_wo_summaries = eff$total_time_to_1000_ESS_wo_summaries,
                    total_time_to_10000_ESS_wo_summaries = eff$total_time_to_10000_ESS_wo_summaries,
                    ##
                    ## Linear extrapolations from observed min ESS; invalid ESS cannot give a usable time estimate.
                    sampling_time_to_target_ESS_mins = if (is.finite(min_ESS) && min_ESS > 0)
                        (target_ESS_val/min_ESS)*time_sampling/60 else NA_real_,
                    est_time_to_target_ESS_wo_summaries_mins = time_burnin/60 + sampling_time_to_target_ESS_mins,
                    ## Preserve the existing convention for the estimate INCLUDING scaled summary time.
                    est_time_to_target_ESS_mins = if (is.finite(min_ESS) && min_ESS > 0)
                        (time_burnin + (target_ESS_val/min_ESS)*(time_sampling + time_summaries))/60 else NA_real_
                  )
                  
                  all_rows[[row_idx]] <-  row
                  row_idx <-  row_idx + 1
                  
                })
              }
            }
          })
        }
        } ## sampler combinations
        
        ## Combine all rows
        if (length(all_rows) > 0) {
          results_df <-  bind_rows(all_rows)
          results_df <-  results_df %>% arrange(N, n_burnin, learning_rate, metric_estimator, tau_initial,
                                                learning_rate_initial, diffusion_HMC_integrator,
                                                clip_iter, int, ratio_M_main, ratio_M_nuisance, int_width, run)
        } else {
          warning("No matching PS7 runs found. Check the selected options, chunks, WCP threads and other filename settings.")
          results_df <-  tibble()
        }
        
        return(results_df)
  
}


# get_target_ESS <-  function(N) { 
#   
#     if (N_val == 500)          {  target_ESS_val <-  10000
#     } else if (N_val == 2500)  {  target_ESS_val <-  5000
#     } else if (N_val == 10000) {  target_ESS_val <-  5000
#     }
#     
#     return(target_ESS_val)
#   
# }
# 
# 
# ##
# ## ---- Functions to run BayesMVP models  -------------------------------------------------------------------------------------------------------
# ##
# run_ps7_models <-  function( Model_type,
#                             N,
#                             y,
#                             settings, 
#                             vect_type,
#                             runs_override = NULL,
#                             use_disk = use_disk,
#                             use_disk_path = "/tmp/hmc_traces",
#                             use_disk_path_post_hoc_dir = "/tmp/constrain_traces",
#                             ##
#                             true_prev,
#                             prior_LKJ_nd,
#                             prior_LKJ_d,
#                             prior_prev_a,
#                             prior_prev_b,
#                             ##
#                             pop,
#                             ##
#                             LC_MVOP_grouping
# ) {
#   
#         require(BayesMVP)
#         ##
#         cat(sprintf("\n========== PS7 Phase 1: Running models for N = %d ==========\n", N))
#         ##
#         n_threads_WCP_burnin   <-  settings$n_threads_WCP_burnin
#         n_threads_WCP_sampling <-  settings$n_threads_WCP_sampling
#         ##
#         n_tests <-  ncol(y)
#         ##
#         # Model_type <-  "LC_MVP" ## ----------------------------------
#         ##
#         if (Model_type == "LC_MVOP") { 
#           N <-  nrow(y)
#           grouping <-  LC_MVOP_grouping
#         } else {
#           N <-  N
#           grouping <-  FALSE
#         }
#         # ##
#         n_class <-  2
#         ##
#         ## ---- Storage for all runs:
#         ##
#         all_runs_results <-  list()
#         ##
#         run_i <-  1 ## debug
#         ##
#         ## ---- Check for existing runs:
#         ##
#         existing_runs <-  c()
#         ##
#         print(paste("clip_iter = ", settings$clip_iter))
#         print(paste("int = ", settings$int))
#         ##
#         for (run_i in 1:settings$n_runs) {
#                   {
#                         file_name_string <-  file.path(settings$output_dir, paste0("ps7_run",
#                                                                                    ##
#                                                                                    "_Model", Model_type, ## ----
#                                                                                    ##
#                                                                                    "_N", N,
#                                                                                    ##
#                                                                                    "_diff_HMC", settings$diffusion_HMC, 
#                                                                                    "_pt_HMC", settings$partitioned_HMC, 
#                                                                                    ##
#                                                                                    "_n_ch_b", settings$n_chains_burnin,
#                                                                                    "_n_ch_s", settings$n_chains_sampling,
#                                                                                    "_n_b", settings$n_burnin,
#                                                                                    "_n_i", settings$n_iter,
#                                                                                    "_LR", settings$learning_rate,
#                                                                                    "_AD", settings$adapt_delta,
#                                                                                    "_clip", settings$clip_iter,
#                                                                                    "_int", settings$int,
#                                                                                    "_ratio_M", settings$ratio_M_main, "_",
#                                                                                    settings$ratio_M_nuisance,
#                                                                                    "_M_typ", settings$metric_type_main,
#                                                                                    "_M_shp", settings$metric_shape_main,
#                                                                                    "_width", settings$int_width,
#                                                                                    ##
#                                                                                    "_M_dcy", settings$M_decay_type,
#                                                                                    ##
#                                                                                    "true_p", true_prev,
#                                                                                    "pi_LKJ", prior_LKJ_nd, "_", prior_LKJ_d,
#                                                                                    "pi_p", prior_prev_a, "_", prior_prev_b
#                                                                                  ))
#                     
#                         if (settings$M_decay_type == "inverse") { 
#                           file_name_string <-  paste0(file_name_string,        
#                                                      "_M_dcy_pow", settings$M_decay_power,  
#                                                      "_M_dcy_scl", settings$M_decay_scale  
#                           )
#                         }
#                         if (grouping == TRUE) { 
#                           file_name_string <-  paste0(file_name_string,        
#                                                      "_ord_grp", grouping)
#                         }
#                         ##
#                         run_file <-  paste0(file_name_string, "_run", run_i)
#                 }
#                 ##
#                 if (file.exists(run_file)) {
#                   existing_runs <-  c(existing_runs, run_i)
#                   all_runs_results[[run_i]] <-  readRDS(run_file)
#                   cat(sprintf("  Loaded existing run %d\n", run_i))
#                 }
#           
#         }
#         ##
#         runs_to_do <-  setdiff(1:settings$n_runs, existing_runs)
#         ##
#         if (length(runs_to_do) == 0) {
#           cat(sprintf("All %d runs complete for N=%d, n_burnin=%d\n",
#                       settings$n_runs, N, settings$n_burnin))
#         } else {
#           cat(sprintf("Running: %s\n", paste(runs_to_do, collapse = ", ")))
#         }
#         ##
#         ## ---- Run models:;
#         ##
#         if (!(is.null(runs_override))) { 
#           runs_to_do <-  runs_override
#         }
#         ##
#         run_i <-  1
#         ##
#         # runs_to_do <-  runs_to_do[-1] ## ------------------------------
#         ##
#         for (run_i in runs_to_do) {
#                 
#                 MCMC_seed <-  run_i * 1000  # DIFFERENT MCMC seed each run
#                 ##
#                 set.seed(MCMC_seed, kind = "L'Ecuyer-CMRG")
#                 ##
#                 cat(sprintf("\n--- Run %d/%d (MCMC seed = %d) ---\n", run_i, settings$n_runs, MCMC_seed))
#                 ##
#                 ## ---- Initialize model:
#                 ##
#                 # if (Model_type == "LC_MVOP") {
#                 #   y[2,] <-  y[2,] + 1
#                 #   y[3,] <-  y[3,] + 1
#                 # }
#                 model_args_list <-  list(y = y, 
#                                         n_tests = n_tests, 
#                                         n_class = 2, 
#                                         N = N)
#                 ##
#                 {
#                     source(file.path("/home/enzocerullo/Documents/Work/PhD_work/R_packages/BayesMVP/inst/BayesMVP/R", 
#                                      "R_fns_misc_BayesMVP.R"))
#                     source(file.path("/home/enzocerullo/Documents/Work/PhD_work/R_packages/BayesMVP/inst/BayesMVP/R", 
#                                      "R_fn_initialise_model.R"))
#                     source(file.path("/home/enzocerullo/Documents/Work/PhD_work/R_packages/BayesMVP/inst/BayesMVP/R", 
#                                      "R_fns_init_hard_coded_models.R"))
#                     source(file.path("/home/enzocerullo/Documents/Work/PhD_work/R_packages/BayesMVP/inst/BayesMVP/R", 
#                                      "R_fns_post_burnin_prep_inits.R"))
#                     source(file.path("/home/enzocerullo/Documents/Work/PhD_work/R_packages/BayesMVP/inst/BayesMVP/R", 
#                                      "R_fns_ChESSR_HMC.R"))
#                     source(file.path("/home/enzocerullo/Documents/Work/PhD_work/R_packages/BayesMVP/inst/BayesMVP/R", 
#                                      "R_fn_burnin_helper_fns.R"))
#                     source(file.path("/home/enzocerullo/Documents/Work/PhD_work/R_packages/BayesMVP/inst/BayesMVP/R", 
#                                      "R_fn_init_and_run_burnin_CHESS.R"))
#                     source(file.path("/home/enzocerullo/Documents/Work/PhD_work/R_packages/BayesMVP/inst/BayesMVP/R", 
#                                      "R_fn_sample.R"))
#                     source(file.path("/home/enzocerullo/Documents/Work/PhD_work/R_packages/BayesMVP/inst/BayesMVP/R", 
#                                      "R_fn_create_summary_and_traces.R"))
#                 }
#                 ##
#                 # stream = MCMC_seed
#                 # sample_nuisance = TRUE
#                 # n_nuisance_override = NULL
#                 # model_args_list = model_args_list
#                 # compile = TRUE
#                 # force_recompile = FALSE
#                 # cmdstanr_model_fit_obj = NULL
#                 # Stan_data_list = NULL
#                 # Stan_model_file_path = NULL
#                 # Stan_cpp_user_header = NULL
#                 # Stan_cpp_flags = NULL
#                 # stanc_args = NULL
#                 ###
#                 n_coeffs_MVP  <-  1 ## intercept-only
#                 n_coeffs_MVOP <-  3
#                 ##
#                 if (Model_type == "LC_MVOP") {
#                       
#                       model_args_list$X <-  list()
#                       ##
#                       for (c in 1:2) {
#                         model_args_list$X[[c]] <-  list()
#                         for (t in 1:3) { 
#                           model_args_list$X[[c]][[t]] <-  matrix(data = 1, nrow = N, ncol = 2)
#                           model_args_list$X[[c]][[t]][1:N, 2] <-  pop
#                         }
#                       }
#                       ## Extract the language column (col 2) from any X matrix:
#                       lang_vec <-  model_args_list$X[[1]][[1]][, 2]
#                       
#                       ## Expand to dummies:
#                       lang_dummies <-  expand_categorical_to_dummies( lang_vec,
#                                                                      ref_level = 1L, ## ----  Note: 1 = Afrikaans, 2 = Xhosa, 3 = Zulu
#                                                                      prefix = "lang")
#                       ## Result: lang_2 and lang_3 columns (lang_1 = reference)
#                       
#                       ## Replace col 2 with the two dummies, keeping intercept:
#                       for (c in 1:n_class) {
#                         for (t in 1:n_tests) {
#                           X_old <-  model_args_list$X[[c]][[t]]
#                           if (n_coeffs_MVOP == 3) model_args_list$X[[c]][[t]] <-  cbind(X_old[, 1, drop = FALSE], lang_dummies)
#                           else                    model_args_list$X[[c]][[t]] <-  lang_dummies ## cbind(X_old[, 1, drop = FALSE], lang_dummies)
#                         }
#                       }
#                       
#                 } else { 
#                   
#                       # model_args_list$X <-  list()
#                       # ##
#                       # for (c in 1:2) {
#                       #   model_args_list$X[[c]] <-  list()
#                       #   for (t in 1:6) {
#                       #     model_args_list$X[[c]][[t]] <-  matrix(data = 1, nrow = N, ncol = 2)
#                       #     model_args_list$X[[c]][[t]][1:N, 2] <-  rcat(n = N, p = c(0.1, 0.1, 0.1))
#                       #     model_args_list$X[[c]][[t]][1:N, 1] <-  rep(0.0, N)
#                       #   }
#                       # }
#                       # ## Extract the language column (col 2) from any X matrix:
#                       # lang_vec <-  model_args_list$X[[1]][[1]][, 2]
#                       # 
#                       # ## Expand to dummies:
#                       # lang_dummies <-  expand_categorical_to_dummies( lang_vec,
#                       #                                                ref_level = 1L, ## ----  Note: 1 = Afrikaans, 2 = Xhosa, 3 = Zulu
#                       #                                                prefix = "lang")
#                       # ## Result: lang_2 and lang_3 columns (lang_1 = reference)
#                       # 
#                       # ## Replace col 2 with the two dummies, keeping intercept:
#                       # for (c in 1:n_class) {
#                       #   for (t in 1:n_tests) {
#                       #     X_old <-  model_args_list$X[[c]][[t]]
#                       #     if (n_coeffs_MVP == 3) model_args_list$X[[c]][[t]] <-  cbind(X_old[, 1, drop = FALSE], lang_dummies)
#                       #     else                   model_args_list$X[[c]][[t]] <-  lang_dummies ## cbind(X_old[, 1, drop = FALSE], lang_dummies)
#                       #   }
#                       # }
#                     
#                 }
#                 ##
#                 # model_args_list$X[[c]][[t]][1:50, ]
#                 ##
#                 init_object <-  initialise_model( Model_type = Model_type,
#                                                   stream = MCMC_seed,
#                                                   sample_nuisance = TRUE,
#                                                   n_nuisance_override = NULL,
#                                                   model_args_list = model_args_list,
#                                                   compile = TRUE,
#                                                   force_recompile = FALSE,
#                                                   cmdstanr_model_fit_obj = NULL,
#                                                   Stan_data_list = NULL,
#                                                   Stan_model_file_path = NULL,
#                                                   Stan_cpp_user_header = NULL,
#                                                   Stan_cpp_flags = NULL,
#                                                   stanc_args = NULL)
#                 ##
#                 # tryCatch(
#                 #   # the existing initialise_model(...) call here
#                 #   initialise_model( Model_type = Model_type,
#                 #                     stream = MCMC_seed,
#                 #                     sample_nuisance = TRUE,
#                 #                     n_nuisance_override = NULL,
#                 #                     model_args_list = model_args_list,
#                 #                     compile = TRUE,
#                 #                     force_recompile = FALSE,
#                 #                     cmdstanr_model_fit_obj = NULL,
#                 #                     Stan_data_list = NULL,
#                 #                     Stan_model_file_path = NULL,
#                 #                     Stan_cpp_user_header = NULL,
#                 #                     Stan_cpp_flags = NULL,
#                 #                     stanc_args = NULL),
#                 #   error = function(e) {
#                 #     writeLines(conditionMessage(e), "/tmp/bs_error.txt")
#                 #     message("Full error written to /tmp/bs_error.txt")
#                 #   }
#                 # )
#                 ##
#                 model_args_list <-  init_object$model_args_list
#                 ##
#                 # model_args_list$overflow_threshold  <-  +100
#                 # model_args_list$underflow_threshold <-  -100
#                 ##
#                 # model_args_list$vect_type <-  "Stan"
#                 ##
#                 ## ---- Set priors:
#                 ##
#                 # prior_LKJ_d <-  4.0 ; prior_LKJ_nd <-  16.0 ## --------------------------
#                 # # prior_LKJ_d <-  3.0 ; prior_LKJ_nd <-  14.0
#                 # # prior_LKJ_d <-  2.0 ; prior_LKJ_nd <-  12.0
#                 # # prior_LKJ_d <-  1.5 ; prior_LKJ_nd <-  10.0
#                 # ##
#                 # if (true_prev == 0.10) { 
#                 #   prior_prev_a <-  2.5 ;  prior_prev_b <-  15
#                 # } else if (true_prev == 0.20) { 
#                 #   # prev_priors <-  c(5, 10)
#                 #   prior_prev_a <-  3 ;  prior_prev_b <-  15
#                 # } else { 
#                 #   prior_prev_a <-  2 ;  prior_prev_b <-  2
#                 # }
#                 ##
#                 if (Model_type == "LC_MVOP") {
#                   
#                       prior_prev_a <-  1.5
#                       prior_prev_b <-  10
#                       ##
#                       model_args_list$n_pops <-  3
#                       model_args_list$pop <-  pop
#                       ##
#                       model_args_list$prior_prev_a <-  rep(prior_prev_a, model_args_list$n_pops)
#                       model_args_list$prior_prev_b <-  rep(prior_prev_b, model_args_list$n_pops)
#                   
#                 } else {
#                       
#                       prior_prev_a <-  2.5
#                       prior_prev_b <-  10
#                       ##
#                       model_args_list$n_pops <-  length(unique(pop))
#                       model_args_list$pop <-  pop
#                       ##
#                       model_args_list$prior_prev_a <-  rep(prior_prev_a, model_args_list$n_pops)
#                       model_args_list$prior_prev_b <-  rep(prior_prev_b, model_args_list$n_pops)
#                   
#                 }
#                 ##
#                 model_args_list$prior_prev_a
#                 model_args_list$prior_prev_b
#                 ##
#                 # prev_priors <-  c(prior_prev_a, prior_prev_b)
#                 ##
#                 # base_prior <-  create_base_prior( DGM = settings$DGM,
#                 #                                  n_tests = n_tests,
#                 #                                  prev_prior_params = prev_priors)
#                 # ##
#                 # model_args_list$prev_prior_a <-  base_prior$prev_prior_a
#                 # model_args_list$prev_prior_b <-  base_prior$prev_prior_b
#                 ##
#                 if (Model_type == "LC_MVP") {
#                     # model_args_list$prior_coeffs_mean_mat <-  base_prior$prior_coeffs_mean_mat
#                     # model_args_list$prior_coeffs_sd_mat   <-  base_prior$prior_coeffs_sd_mat
#                     model_args_list$prior_coeffs_mean_mat <-  rep(list(matrix(0.0, nrow = n_coeffs_MVP, ncol = n_tests)), 2)
#                     model_args_list$prior_coeffs_sd_mat   <-  rep(list(matrix(1.0, nrow = n_coeffs_MVP, ncol = n_tests)), 2)
#                 } else { 
#                     model_args_list$prior_coeffs_mean_mat <-  rep(list(matrix(0.0, nrow = n_coeffs_MVOP, ncol = n_tests)), 2)
#                     model_args_list$prior_coeffs_sd_mat   <-  rep(list(matrix(1.0, nrow = n_coeffs_MVOP, ncol = n_tests)), 2)
#                 }
#                 # ##
#                 # if (settings$DGM %in% c(2, 3)) {
#                 #   # model_args_list$prior_coeffs_mean_mat[[1]][1, 1] <-  -2.33
#                 #   model_args_list$prior_coeffs_mean_mat[[1]][1, 1] <-  -1.7
#                 #   model_args_list$prior_coeffs_mean_mat[[2]][1, 1] <-  +0.385 ## D+
#                 #   ##
#                 #   # model_args_list$prior_coeffs_sd_mat[[1]][1, 1] <-  0.50 ## 0.375
#                 #   # model_args_list$prior_coeffs_sd_mat[[2]][1, 1] <-  0.45 ## D+ ## 0.25
#                 #   ##
#                 #   # model_args_list$prior_coeffs_sd_mat[[1]][1, 1] <-  0.375 ## 0.375
#                 #   model_args_list$prior_coeffs_sd_mat[[1]][1, 1] <-  0.30 ## 0.375
#                 #   model_args_list$prior_coeffs_sd_mat[[2]][1, 1] <-  0.25 ## D+ ## 0.25
#                 #   ##
#                 # } else if (settings$DGM == 5) { 
#                 #   model_args_list$prior_coeffs_mean_mat[[1]][1, 1]  <-  -1.30
#                 #   model_args_list$prior_coeffs_mean_mat[[2]][1, 1]  <-  +1.30
#                 #   ##
#                 #   model_args_list$prior_coeffs_sd_mat[[1]][1, 1] <-  0.50
#                 #   model_args_list$prior_coeffs_sd_mat[[2]][1, 1] <-  0.50
#                 # }
#                 if (Model_type == "LC_MVOP") {
#                   
#                         for (c in 1:2) {
#                           model_args_list$prior_coeffs_sd_mat[[c]][,] <-  0.50
#                           ##
#                           model_args_list$prior_coeffs_mean_mat[[c]][1, c(2:n_tests)] <-  0.0
#                           model_args_list$prior_coeffs_sd_mat[[c]][1, c(2:n_tests)]   <-  2.0 ## 1.0 ##  1.0 ## 2.0
#                           ##
#                           # model_args_list$prior_coeffs_sd_mat[[c]][2, 1] <-  0.25
#                           try({  
#                             model_args_list$prior_coeffs_sd_mat[[c]][2:3, 1] <-  0.50 ## 0.25 ## for Ref test (BINARY)
#                           })
#                         }
#                         ##
#                         model_args_list$prior_coeffs_mean_mat[[1]][1, 1]  <-  -1.5
#                         model_args_list$prior_coeffs_mean_mat[[2]][1, 1]  <-  +0.5 ## D+
#                         ##
#                         model_args_list$prior_coeffs_sd_mat[[1]][1, 1] <-  0.25 ## 1.0 ## 0.50 ## 0.25
#                         model_args_list$prior_coeffs_sd_mat[[2]][1, 1] <-  0.25 ## 1.0 ## 0.50 ## 0.25 ## D+
#                         # model_args_list$prior_coeffs_sd_mat[[1]][1, 1] <-  0.01 ## 1.0 ## 0.50 ## 0.25
#                         # model_args_list$prior_coeffs_sd_mat[[2]][1, 1] <-  0.01 ## 1.0 ## 0.50 ## 0.25 ## D+
#                         ##
#                         model_args_list$prior_coeffs_mean_mat
#                         model_args_list$prior_coeffs_sd_mat
#                     
#                 }  else {
#                   
#                         # model_args_list$prior_coeffs_mean_mat[[1]][1, 1]  <-  -1.75
#                         # model_args_list$prior_coeffs_mean_mat[[2]][1, 1]  <-  +0.385 ## +1.75 ## D+
#                         # ##
#                         # model_args_list$prior_coeffs_sd_mat[[1]][1, 1] <-  0.35
#                         # model_args_list$prior_coeffs_sd_mat[[2]][1, 1] <-  0.45 ## 0.35 ## D+
#                         ##
#                         ## All tests: Sp roughly 97-99%
#                         ## beta_nd ~ N(-2.0, 0.5) -> median Sp ~ 97.7%, 95% CI ~ (84%, 99.9%)
#                         model_args_list$prior_coeffs_mean_mat[[1]][1, ] <-  -2.0
#                         model_args_list$prior_coeffs_sd_mat[[1]][1, ]   <-  0.5
#                         ##
#                         ## Lab tests (euroimmun, roche): Se roughly 90-98%
#                         ## beta_d ~ N(1.5, 0.5) -> median Se ~ 93%, 95% CI ~ (69%, 99.4%)
#                         model_args_list$prior_coeffs_mean_mat[[2]][1, 1:2] <-  1.5  ## euroimmun, roche
#                         model_args_list$prior_coeffs_sd_mat[[2]][1, 1:2]   <-  0.5
#                         ##
#                         ## Lateral flow (abc19, surescreen, orientgene, biomerica): Se roughly 70-90%
#                         ## beta_d ~ N(1.0, 0.5) -> median Se ~ 84%, 95% CI ~ (50%, 97.7%)
#                         model_args_list$prior_coeffs_mean_mat[[2]][1, 3:6] <-  1.0
#                         model_args_list$prior_coeffs_sd_mat[[2]][1, 3:6]   <-  0.5
#                         
#                 }
#                 ##
#                 model_args_list$prior_coeffs_mean_mat
#                 model_args_list$prior_coeffs_sd_mat
#                 ##
#                 ## ---- Correlation priors:
#                 ##
#                 corr_prior_sets <-  get_correlation_prior_sets()
#                 corr_prior_settings <-  corr_prior_sets$LC_MVP[[paste0("DGM_", settings$DGM)]][[settings$prior_name]]
#                 ##
#                 # model_args_list$prior_LKJ <-  matrix(c(corr_prior_settings$eta_nd, corr_prior_settings$eta_d), ncol = 1)
#                 ##
#                 if (Model_type == "LC_MVOP") { 
#                     prior_LKJ_nd <-  4 ## -------------------------------------
#                     prior_LKJ_d  <-  4 ## -------------------------------------
#                     # prior_LKJ_nd <-  3 ## -------------------------------------
#                     # prior_LKJ_d  <-  3 ## -------------------------------------
#                     # prior_LKJ_nd <-  2 ## -------------------------------------
#                     # prior_LKJ_d  <-  2 ## -------------------------------------
#                     # prior_LKJ_nd <-  1 ## -------------------------------------
#                     # prior_LKJ_d  <-  1 ## -------------------------------------
#                 } else { 
#                   prior_LKJ_nd <-  4 ## -------------------------------------
#                   prior_LKJ_d  <-  4 ## -------------------------------------
#                   # prior_LKJ_nds <-  3 ## -------------------------------------
#                   # prior_LKJ_d  <-  3 ## -------------------------------------
#                   # prior_LKJ_nd <-  2 ## -------------------------------------
#                   # prior_LKJ_d  <-  2 ## -------------------------------------
#                 }
#                 ##
#                 model_args_list$prior_LKJ <-  matrix(c(prior_LKJ_nd, prior_LKJ_d), ncol = 1) ## --------------------------
#                 ##
#                 model_args_list$lkj_cholesky_eta <-  model_args_list$prior_LKJ
#                 model_args_list$corr_force_positive <-  corr_prior_settings$corr_force_positive
#                 ##
#                 if (Model_type == "LC_MVOP") { 
#                    # model_args_list$corr_force_positive <-  TRUE
#                    model_args_list$corr_force_positive <-  FALSE
#                 } else {
#                    model_args_list$corr_force_positive <-  FALSE
#                 }
#                 # model_args_list$corr_force_positive <-  TRUE ## --------------------------
#                 ##
#                 if (isTRUE(model_args_list$corr_force_positive)) {
#                   for (c in 1:n_class) model_args_list$lb_corr[[c]] <-  matrix(0.0, ncol = n_tests, nrow = n_tests)
#                 } else { 
#                   for (c in 1:n_class) model_args_list$lb_corr[[c]] <-  matrix(-1.0, ncol = n_tests, nrow = n_tests)
#                 }
#                 for (c in 1:n_class) model_args_list$ub_corr[[c]] <-  matrix(+1.0, ncol = n_tests, nrow = n_tests)
#                 ##
#                 if (!is.null(corr_prior_settings$lb_corr)) {
#                   if (Model_type == "LC_MVOP") {
#                       model_args_list$lb_corr <-  matrix(-1.0, ncol = n_tests, nrow = n_tests)
#                       model_args_list$ub_corr <-  matrix(+1.0, ncol = n_tests, nrow = n_tests)
#                   } else {
#                       model_args_list$lb_corr <-  corr_prior_settings$lb_corr
#                       model_args_list$ub_corr <-  corr_prior_settings$ub_corr
#                   }
#                 }
#                 ##
#                 model_args_list$n_covariates_per_outcome_mat
#                 ##
#                 # model_args_list$n_covariates_max
#                 ##
#                 ## ---- Inits:
#                 ##
#                 init_lists_per_chain <-  list()
#                 ##
#                 for (kk in 1:settings$n_chains_burnin) {
#                   
#                       prev_init_LC_MVP  <-  0.10
#                       ##
#                       prev_init_LC_MVOP <-  0.05
#                       ##
#                       # random_draws <-  rnorm(n = n_tests*N, mean = 0.0, sd = 0.10)
#                       random_draws <-  rnorm(n = n_tests*N, mean = 0.0, sd = 0.01)
#                       ##
#                       n_corrs <-  n_tests * (n_tests - 1) / 2
#                       ##
#                       n_thr_total <-  sum(model_args_list$n_thr_per_ord_test)
#                       ##
#                       if (Model_type == "LC_MVOP") {
#                         
#                               init_lists_per_chain[[kk]] <-  list(
#                                 u_raw = matrix(random_draws, ncol = n_tests, nrow = N),
#                                 Omega_unconstrained_vec = list(rep(0.001, n_corrs), rep(0.001, n_corrs)),
#                                 beta = list(matrix(-1.0, nrow = n_coeffs_MVOP, ncol = n_tests), 
#                                             matrix(+1.0, nrow = n_coeffs_MVOP, ncol = n_tests)),
#                                 # p_raw = array(atanh(2 * prev_init_LC_MVOP - 1)),
#                                 p_raw = array(rep(atanh(2 * prev_init_LC_MVOP - 1)), model_args_list$n_pops),
#                                 C_raw_vec = rep(list(rep(-3.0, n_thr_total)), n_class)
#                               )
#                               
#                               for (c in 1:2) {
#                                   init_lists_per_chain[[kk]]$beta[[c]][c(2, 3), ] <-  0.0
#                               }
#                           
#                       } else if (Model_type == "LC_MVP") {
#                             
#                               init_lists_per_chain[[kk]] <-  list(
#                                 u_raw = matrix(random_draws, ncol = n_tests, nrow = N),
#                                 Omega_unconstrained_vec = list(rep(0.001, n_corrs), rep(0.001, n_corrs)),
#                                 beta = list(matrix(-1.0, nrow = n_coeffs_MVP, ncol = n_tests), 
#                                             matrix(+1.0, nrow = n_coeffs_MVP, ncol = n_tests)),
#                                 # p_raw = array(atanh(2 * prev_init_LC_MVP - 1))
#                                 p_raw = array(rep(atanh(2 * prev_init_LC_MVP - 1)), model_args_list$n_pops)
#                               )
#                               
#                               try({
#                                 for (c in 1:2) {
#                                   init_lists_per_chain[[kk]]$beta[[1]][1, ] <-  -2.0
#                                   ##
#                                   init_lists_per_chain[[kk]]$beta[[2]][1, 1:2] <-  +1.5
#                                   init_lists_per_chain[[kk]]$beta[[2]][1, 3:6] <-  +1.0
#                                 }
#                               })
#                               
#                       }
#                       
#                       # init_lists_per_chain[[kk]]$beta[[1]][1, 1] <-  -2.0
#                   
#                 }
#                 init_lists_per_chain <-  BayesMVP:::resize_init_list( init_lists_per_chain = init_lists_per_chain, 
#                                                                         n_chains_new = settings$n_chains_burnin)
#                 ##
#                 # model_args_list$num_chunks <-  ifelse(N > 2500, 6, 1)
#                 # model_args_list$J_grad_option <-  "autodiff"
#                 ##
#                 model_args_list$num_chunks <-  BayesMVP:::find_num_chunks_MVP( N = N,
#                                                                                  n_tests = n_tests)
#                 if (N == 10000) { 
#                     model_args_list$J_grad_option <-  "num_diff"
#                     # model_args_list$J_grad_option <-  "autodiff"
#                     ##
#                     # model_args_list$num_chunks <-  10
#                     model_args_list$num_chunks <-  25 ## 20 ## 20 ## 50 ## 10
#                     # model_args_list$num_chunks <-  20
#                     # model_args_list$num_chunks <-  40
#                     # model_args_list$num_chunks <-  50
#                 } else if (N == 2500) { 
#                     model_args_list$J_grad_option <-  "num_diff"
#                     # model_args_list$J_grad_option <-  "autodiff"
#                     ##
#                     # model_args_list$num_chunks <-  20
#                     model_args_list$num_chunks <-  10
#                     # model_args_list$num_chunks <-  8
#                     # model_args_list$num_chunks <-  5 ## 4
#                     # model_args_list$num_chunks <-  4
#                     # model_args_list$num_chunks <-  1
#                     # model_args_list$num_chunks <-  25
#                     # model_args_list$num_chunks <-  8
#                     # model_args_list$num_chunks <-  3
#                 } else if (N == 500) { 
#                     model_args_list$J_grad_option <-  "num_diff"
#                     # model_args_list$J_grad_option <-  "autodiff"
#                     ##
#                     model_args_list$num_chunks <-  1
#                     # model_args_list$num_chunks <-  2
#                 } else if (N == 25000) {
#                     model_args_list$J_grad_option <-  "num_diff"
#                     ##
#                     # model_args_list$num_chunks <-  50
#                     model_args_list$num_chunks <-  100
#                 } else if (N == 50000) {
#                     model_args_list$J_grad_option <-  "num_diff"
#                     ##
#                     # model_args_list$num_chunks <-  50
#                     # model_args_list$num_chunks <-  100
#                     model_args_list$num_chunks <-  200
#                 } else if (N == 5000) { ## ----------------------------------- ordinal-only
#                     model_args_list$num_chunks <-  5
#                     ##
#                     # model_args_list$J_grad_option <-  "autodiff"
#                     model_args_list$J_grad_option <-  "num_diff"
#                 } else if (N == 1000) { ## ----------------------------------- ordinal-only
#                     model_args_list$num_chunks <-  1
#                     # model_args_list$num_chunks <-  2
#                     ##
#                     # model_args_list$J_grad_option <-  "autodiff"
#                     model_args_list$J_grad_option <-  "num_diff"
#                 } else if (N == 250) {  ## ----------------------------------- ordinal-only
#                     model_args_list$num_chunks <-  1
#                     # model_args_list$num_chunks <-  2
#                     ##
#                     # model_args_list$J_grad_option <-  "autodiff"
#                     model_args_list$J_grad_option <-  "num_diff"
#                 } else {
#                   model_args_list$num_chunks <-  1
#                   # model_args_list$num_chunks <-  2
#                   ##
#                   # model_args_list$J_grad_option <-  "autodiff"
#                   model_args_list$J_grad_option <-  "num_diff"
#                 }
#                 ##
#                 # if (Model_type == "LC_MVP") {
#                 #   
#                 #     if (N == 10000) { 
#                 #       settings$n_threads_WCP_burnin <-  25
#                 #     } else if (N == 2500) { 
#                 #       settings$n_threads_WCP_burnin <-  10
#                 #     } else if (N == 500) { 
#                 #       settings$n_threads_WCP_burnin <-  1
#                 #     }
#                 #   
#                 # } else if (Model_type == "LC_MVOP") { 
#                 #   
#                 # }
#                 # ##
#                 if (parallel::detectCores() > 16) { ## Local-HPC
#                   settings$n_threads_WCP_burnin   <-  model_args_list$num_chunks
#                   settings$n_threads_WCP_sampling <-  1 ## model_args_list$num_chunks
#                 } else {  ## Laptop
#                   settings$n_threads_WCP_burnin   <-  model_args_list$num_chunks
#                   settings$n_threads_WCP_sampling <-  1 ## model_args_list$num_chunks
#                 }
#                 ##
#                 # model_args_list$J_grad_option <-  "autodiff" ## ----------------------------------- debug
#                 ##
#                 # model_args_list$num_chunks <-  1 ## ----------------------------------- debug
#                 # model_args_list$num_chunks <-  5 ## ----------------------------------- debug
#                 # model_args_list$num_chunks <-  2 ## ----------------------------------- debug
#                 # model_args_list$J_grad_option <-  "autodiff"  ## ----------------------------------- debug
#                 ##
#                 print(paste("num_chunks = ", model_args_list$num_chunks))
#                 ##
#                 if (model_args_list$num_chunks < settings$n_threads_WCP_burnin) {
#                   stop(" n_chunks must be ≥ n_threads_WCP_burnin!")
#                 }
#                 if (model_args_list$num_chunks < settings$n_threads_WCP_sampling) {
#                   stop(" n_chunks must be ≥ n_threads_WCP_sampling!")
#                 }
#                 ##
#                 model_args_list$nuisance_transformation
#                 ##
#                 clip_iter_tau <-  settings$clip_iter + settings$int ## + round(clip_iter/2)
#                 # gap <-  clip_iter + clip_iter_tau
#                 ##
#                 Rcpp::sourceCpp("~/Documents/Work/PhD_work/R_packages/BayesMVP/inst/BayesMVP/src/burnin_prep_fns.cpp")
#                 ##
#                 int_width <-  settings$int_width
#                 ##
#                 if (settings$n_burnin == 1000) {
#                   n_adapt <-  900
#                 } else {
#                   n_adapt <-  settings$n_burnin - round(settings$n_burnin/10)
#                 }
#                 ##
#                 # str(model_args_list)
#                 ##
#                 # {
#                 #     debug = FALSE
#                 #     stream = MCMC_seed
#                 #     ##
#                 #     init_object = init_object
#                 #     ##
#                 #     n_chains_burnin = settings$n_chains_burnin
#                 #     init_lists_per_chain = init_lists_per_chain
#                 #     ##
#                 #     parallel_method = "RcppParallel"
#                 #     ##
#                 #     Stan_data_list = NULL
#                 #     model_args_list = model_args_list
#                 #     ##
#                 #     sample_nuisance = TRUE
#                 #     n_nuisance_override = NULL
#                 #     ##
#                 #     seed = MCMC_seed
#                 #     ##
#                 #     n_burnin = settings$n_burnin
#                 #     n_adapt = n_adapt
#                 #     gap = NULL
#                 #     ##
#                 #     n_chains_sampling = settings$n_chains_sampling
#                 #     n_superchains = settings$n_chains_sampling
#                 #     n_iter = settings$n_iter
#                 #     ##
#                 #     adapt_delta = settings$adapt_delta
#                 #     learning_rate = settings$learning_rate
#                 #     ##
#                 #     tau_mult = 1.60
#                 #     tau_initial = 2*pi ## bookmark: this (2*pi) is used in a Python algorithm (source pending verification)
#                 #     ##
#                 #     manual_tau = FALSE
#                 #     # tau_if_manual = c(6.0, 3.0)
#                 #     tau_if_manual = c(3.0, 9.0)
#                 #     ##
#                 #     burnin_algorithm = "ChESSR"
#                 #     ##
#                 #     diffusion_HMC   = settings$diffusion_HMC ## -----------------
#                 #     partitioned_HMC = settings$partitioned_HMC
#                 #     ##
#                 #     clip_iter = settings$clip_iter
#                 #     clip_iter_tau = clip_iter_tau
#                 #     ##
#                 #     n_refresh = 100
#                 #     use_proposed = TRUE ## -----------------
#                 #     ##
#                 #     beta1_adam = 0.00
#                 #     beta2_adam = 0.95
#                 #     eps_adam = 1e-8
#                 #     ##
#                 #     force_autodiff = FALSE
#                 #     force_PartialLog = FALSE
#                 #     multi_attempts = FALSE
#                 #     ## multi_attempts = FALSE
#                 #     ##
#                 #     force_autodiff_for_metric = TRUE
#                 #     force_PartialLog_for_metric = TRUE
#                 #     force_multi_attempts_for_metric = FALSE
#                 #     ##
#                 #     vect_type = vect_type
#                 #     Phi_type = "Phi"
#                 #     inv_Phi_type = "inv_Phi"
#                 #     ##
#                 #     # metric_type_main = "Hessian"
#                 #     # metric_shape_main = "dense"
#                 #     # ratio_M_main = 0.25
#                 #     # interval_width_main = 25
#                 #     ##
#                 #     metric_type_main = settings$metric_type_main ## "Empirical"
#                 #     metric_shape_main = settings$metric_shape_main ##  "dense"
#                 #     ratio_M_main = settings$ratio_M_main
#                 #     interval_width_main = settings$int_width ## max(1, floor(0.5*settings$n_burnin/200))
#                 #     ##
#                 #     M_decay_type  = settings$M_decay_type
#                 #     M_decay_power = settings$M_decay_power
#                 #     M_decay_scale = settings$M_decay_scale
#                 #     ##
#                 #     # metric_type_main = "Empirical"
#                 #     # metric_shape_main = "diag"
#                 #     # ratio_M_main = 0.25
#                 #     # interval_width_main = 25
#                 #     ##
#                 #     # metric_type_nuisance = "Empirical"
#                 #     # metric_type_nuisance = "uniform_diag"
#                 #     metric_type_nuisance = ifelse(settings$diffusion_HMC == TRUE, "unit", "uniform_diag")
#                 #     ##
#                 #     metric_shape_nuisance = "diag"
#                 #     ##
#                 #     ratio_M_nuisance = settings$ratio_M_nuisance
#                 #     interval_width_nuisance = settings$int_width ## max(1, floor(0.5*settings$n_burnin/200))
#                 #     ##
#                 #     max_tau_main = 50
#                 #     max_tau_nuisance = 50
#                 #     ##
#                 #     max_eps_main = 1.00
#                 #     max_eps_nuisance = 100.00
#                 #     ##
#                 #     max_L = 1024
#                 #     n_nuisance_to_track = NULL
#                 #     ##
#                 #     use_disk = use_disk
#                 #     use_disk_path = use_disk_path
#                 # 
#                 # manual_tau <-  TRUE
#                 #     ##
#                 #     # if (N == 500)   tau_if_manual <-  3.0
#                 #     # if (N == 2500)  tau_if_manual <-  6.0
#                 #     # if (N == 10000) tau_if_manual <-  7.5
#                 #     # if (N == 25000) tau_if_manual <-  10.0
#                 #     # if (N == 50000) tau_if_manual <-  10.0
#                 # }
#                 ##
#                 if (settings$partitioned_HMC == TRUE) {
#                   manual_tau <-  TRUE ; tau_if_manual <-  c(3.0, 1.50)
#                 } else { 
#                   # manual_tau <-  TRUE ; tau_if_manual <-  c(3.0, 3.0)
#                   # manual_tau <-  TRUE ; tau_if_manual <-  c(1.50, 1.50)
#                   ##
#                   manual_tau <-  FALSE ; tau_if_manual <-  c(3.0, 3.0)
#                 }
#                 ##
#                 model_results <-  BayesMVP::R_fn_sample_model(
#                   debug = FALSE,
#                   stream = MCMC_seed,
#                   ##
#                   init_object = init_object,
#                   ##
#                   n_chains_burnin = settings$n_chains_burnin,
#                   init_lists_per_chain = init_lists_per_chain,
#                   ##
#                   parallel_method = "RcppParallel",
#                   # parallel_method = "OpenMP",
#                   ##
#                   Stan_data_list = NULL,
#                   model_args_list = model_args_list,
#                   ##
#                   sample_nuisance = TRUE,
#                   n_nuisance_override = NULL,
#                   ##
#                   seed = MCMC_seed,
#                   ##
#                   n_burnin = settings$n_burnin,
#                   n_adapt = n_adapt,
#                   gap = NULL,
#                   ##
#                   n_chains_sampling = settings$n_chains_sampling,
#                   n_superchains = settings$n_superchains, ## -----------------
#                   n_iter = settings$n_iter, ## -----------------
#                   ##
#                   adapt_delta = settings$adapt_delta,
#                   learning_rate = settings$learning_rate,
#                   ##
#                   tau_mult = 1.60,
#                   tau_initial = 2*pi, ## bookmark: this (2*pi) is used in a Python algorithm (source pending verification)
#                   ##
#                   manual_tau = manual_tau,
#                   tau_if_manual = tau_if_manual,
#                   ##
#                   burnin_algorithm = "ChESSR",
#                   ##
#                   diffusion_HMC   = settings$diffusion_HMC, ## -----------------
#                   partitioned_HMC = settings$partitioned_HMC,
#                   ##
#                   clip_iter = settings$clip_iter,
#                   clip_iter_tau = clip_iter_tau,
#                   ##
#                   n_refresh = 100,
#                   # use_proposed = TRUE, ## -----------------
#                   use_proposed = FALSE, ## -----------------
#                   ##
#                   beta1_adam = 0.00,
#                   beta2_adam = 0.95,
#                   eps_adam = 1e-8,
#                   ##
#                   force_autodiff = FALSE,
#                   force_PartialLog = FALSE,
#                   multi_attempts = FALSE, ## -----------------
#                   ##
#                   force_autodiff_for_metric = FALSE, ## -----------------
#                   force_PartialLog_for_metric = FALSE,
#                   force_multi_attempts_for_metric = TRUE,
#                   ##
#                   vect_type = vect_type,
#                   Phi_type = "Phi",
#                   inv_Phi_type = "inv_Phi",
#                   ##
#                   metric_type_main = settings$metric_type_main, ## "Empirical",
#                   metric_shape_main = settings$metric_shape_main, ##  "dense",
#                   ratio_M_main = settings$ratio_M_main,
#                   interval_width_main = settings$int_width, ## max(1, floor(0.5*settings$n_burnin/200)),
#                   ##
#                   M_decay_type  = settings$M_decay_type,
#                   M_decay_power = settings$M_decay_power,
#                   M_decay_scale = settings$M_decay_scale,
#                   ##
#                   metric_type_nuisance = "Empirical",
#                   # metric_type_nuisance = "uniform_diag",
#                   # metric_type_nuisance = ifelse(settings$diffusion_HMC == TRUE, "unit", "uniform_diag"),
#                   # metric_type_nuisance = "Empirical",
#                   ##
#                   metric_shape_nuisance = "diag",
#                   ratio_M_nuisance = settings$ratio_M_nuisance,
#                   interval_width_nuisance = settings$int_width, ## max(1, floor(0.5*settings$n_burnin/200)),
#                   ##
#                   max_tau_main = 50,
#                   max_tau_nuisance = 50,
#                   ##
#                   max_eps_main = 1.00,
#                   max_eps_nuisance = 100.00,
#                   ##
#                   max_L = 1024,
#                   n_nuisance_to_track = 1,
#                   ##
#                   use_disk = use_disk,
#                   use_disk_path = use_disk_path,
#                   ##
#                   n_threads_WCP_burnin = settings$n_threads_WCP_burnin,
#                   n_threads_WCP_sampling = settings$n_threads_WCP_sampling)
#                 # ##
#                 # ## ---- Create summaries:
#                 # ##
#                 # {
#                 #   compute_main_params = TRUE
#                 #   compute_transformed_parameters = TRUE
#                 #   compute_generated_quantities = TRUE
#                 #   ##
#                 #   save_log_lik_trace = FALSE
#                 #   save_nuisance_trace = FALSE
#                 #   compute_nested_rhat = TRUE
#                 #   ##
#                 #   n_superchains = NULL
#                 #   save_trace_tibbles = FALSE
#                 #   ##
#                 #   n_iter_to_store = NULL
#                 # 
#                 ##
#                 # Stan_data_list <-  model_results$init_object$Stan_data_list
#                 ##
#                 # str(model_results$sampling_object) ## [[1]]
#                 ##
#                 model_fit_object <-  BayesMVP::create_summary_and_traces( 
#                   model_results = model_results,
#                   ##
#                   compute_main_params = TRUE,
#                   compute_transformed_parameters = TRUE,
#                   compute_generated_quantities = TRUE,
#                   ##
#                   save_log_lik_trace = FALSE,
#                   # save_log_lik_trace = TRUE,
#                   ##
#                   save_nuisance_trace = FALSE,
#                   compute_nested_rhat = TRUE,
#                   ##
#                   n_superchains = NULL,
#                   save_trace_tibbles = FALSE,
#                   ##
#                   use_disk = use_disk,
#                   use_disk_path = use_disk_path,
#                   use_disk_path_post_hoc_dir = use_disk_path_post_hoc_dir)
#                 ##
#                 # str(model_fit_object$traces$log_lik_trace)
#                 # # log_lik_trace_default_corrs <-  model_fit_object$traces$log_lik_trace
#                 # # log_lik_trace_force_pos_corrs <-  model_fit_object$traces$log_lik_trace
#                 # require(loo)
#                 # # loo(log_lik_trace_default_corrs, log_lik_trace_force_pos_corrs)
#                 # # ?loo
#                 # log_lik_trace_default_corrs_array <-  log_lik_trace_force_pos_corrs_array <-  array(dim = c(settings$n_iter, settings$n_chains_sampling, N))
#                 # for (c in 1:settings$n_chains_sampling) {
#                 #    for (n in 1:N) {
#                 #      log_lik_trace_default_corrs_array[, c, n] <-  log_lik_trace_default_corrs[[c]][n, ]
#                 #      log_lik_trace_force_pos_corrs_array[, c, n] <-  log_lik_trace_force_pos_corrs[[c]][n, ]
#                 #    }
#                 # }
#                 # # loo_default_corrs   <-  loo(log_lik_trace_default_corrs_array)
#                 # loo_force_pos_corrs <-  loo(log_lik_trace_force_pos_corrs_array)
#                 # ##
#                 # loo_compare(loo_default_corrs, loo_force_pos_corrs)
#                 ##
#                 # > tibble_Se
#                 # # A tibble: 6 × 10
#                 # param_group parameter       mean      sd `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
#                 # <chr>       <chr>          <dbl>   <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl>  <dbl>
#                 # 1 Se_baseline Se_baseline[1] 0.964 0.0120   0.937 0.965   0.984 14050  1.01   1.01
#                 # 2 Se_baseline Se_baseline[2] 0.990 0.00770  0.970 0.991   0.999 13133  1.01   1.01
#                 # 3 Se_baseline Se_baseline[3] 0.879 0.0187   0.840 0.879   0.913 26332  1.00   1.00
#                 # 4 Se_baseline Se_baseline[4] 0.934 0.0148   0.902 0.934   0.959 18638  1.01   1.01
#                 # 5 Se_baseline Se_baseline[5] 0.961 0.0116   0.935 0.962   0.981 15992  1.01   1.01
#                 # 6 Se_baseline Se_baseline[6] 0.939 0.0142   0.908 0.940   0.963 23152  1.01   1.00
#                 # > tibble_Sp
#                 # # A tibble: 6 × 10
#                 # param_group parameter       mean      sd `2.5%` `50%` `97.5%` n_eff  Rhat n_Rhat
#                 # <chr>       <chr>          <dbl>   <dbl>  <dbl> <dbl>   <dbl> <dbl> <dbl>  <dbl>
#                 # 1 Sp_baseline Sp_baseline[1] 0.992 0.00198  0.988 0.993   0.996 23722  1.01   1.00
#                 # 2 Sp_baseline Sp_baseline[2] 0.993 0.00207  0.988 0.993   0.996 13001  1.01   1.01
#                 # 3 Sp_baseline Sp_baseline[3] 0.989 0.00219  0.985 0.990   0.993 44955  1.00   1.00
#                 # 4 Sp_baseline Sp_baseline[4] 0.975 0.00334  0.968 0.975   0.981 47727  1.00   1.00
#                 # 5 Sp_baseline Sp_baseline[5] 0.973 0.00348  0.966 0.974   0.980 36549  1.00   1.00
#                 # 6 Sp_baseline Sp_baseline[6] 0.920 0.00578  0.908 0.920   0.931 48566  1.00   1.00
#                 # > tibble_prev
#                 # # A tibble: 2 × 10
#                 # param_group parameter   mean      sd `2.5%`  `50%` `97.5%` n_eff  Rhat n_Rhat
#                 # <chr>       <chr>      <dbl>   <dbl>  <dbl>  <dbl>   <dbl> <dbl> <dbl>  <dbl>
#                 # 1 p           p[1]      0.0784 0.00799 0.0636 0.0782  0.0945 59730  1.00   1.00
#                 # 2 p           p[2]      0.175  0.0102  0.156  0.175   0.195  53322  1.00   1.00
#                 ##
#                 ## ---- Extract:
#                 ##
#                 tibble_main <-  model_fit_object$summaries$summary_tibbles$summary_tibble_main_params
#                 tibble_tp   <-  model_fit_object$summaries$summary_tibbles$summary_tibble_transformed_parameters
#                 tibble_gq   <-  model_fit_object$summaries$summary_tibbles$summary_tibble_generated_quantities
#                 ##
#                 # tibble_tp %>% print(n = 100)
#                 ##
#                 tibble_Se   <-  BayesMVP::extract_params_from_tibble_batch(tibble = tibble_gq, param_strings_vec = "Se_baseline")
#                 tibble_Sp   <-  BayesMVP::extract_params_from_tibble_batch(tibble = tibble_gq, param_strings_vec = "Sp_baseline")
#                 tibble_prev <-  BayesMVP::extract_params_from_tibble_batch(tibble = tibble_gq, param_strings_vec = "p")
#                 ##
#                 trace_main <-  model_fit_object$traces$traces_as_arrays$trace_params_main
#                 trace_tp   <-  model_fit_object$traces$traces_as_arrays$trace_transformed_params
#                 trace_gq   <-  model_fit_object$traces$traces_as_arrays$trace_generated_quantities
#                 ##
#                 draws_mat <-  matrix(trace_main, nrow = prod(dim(trace_main)[1:2]), ncol = dim(trace_main)[3])
#                 S_post <-  cov(draws_mat)
#                 ge <-  eigen(solve(model_results$burnin_object$EHMC_Metric_as_Rcpp_List$M_inv_dense_main, S_post))$values
#                 print(paste("max(Re(ge)) / min(Re(ge)) = ", max(Re(ge)) / min(Re(ge))))
#                 ##
#                 EHMC_Metric_as_Rcpp_List <-  model_results$burnin_object$EHMC_Metric_as_Rcpp_List
#                 M_inv <-  EHMC_Metric_as_Rcpp_List$M_inv_dense_main
#                 ev <-  eigen(M_inv, symmetric = TRUE, only.values = TRUE)$values
#                 summary(ev); min(ev)                          # near-zero eigenvalues?
#                 summary(diag(M_inv) / diag(S_post))           # scale: is M_inv ~ posterior var, or miles off?
#                 ##
#                 max(abs(EHMC_Metric_as_Rcpp_List$M_dense_main %*% M_inv - diag(nrow(M_inv))))
#                 ##
#                 # ##
#                 # ## after sampling, for the main params (using the trace):
#                 # 
#                 # ## Flatten iter × chain per parameter (valid because param is the LAST dim):
#                 # draws_mat <-  matrix(trace_main, nrow = 500 * 180, ncol = dim(trace_main)[3])
#                 # post_var_main <-  matrixStats::colVars(draws_mat)   # pooled posterior variance per param
#                 # 
#                 # ## The decisive ratio:
#                 # M_inv <-  model_results$burnin_object$EHMC_Metric_as_Rcpp_List$M_inv_main_vec   # or burnin_object$EHMC_Metric_as_Rcpp_List
#                 # stopifnot(length(M_inv) == dim(trace_main)[3])
#                 # ratio <-  M_inv / post_var_main
#                 # summary(ratio)
#                 # round(cbind(M_inv = head(M_inv, 10), post_var = head(post_var_main, 10), ratio = head(ratio, 10)), 6)
#                 # 
#                 # 
#                 # per_chain_var  <-  apply(trace_main, c(2, 3), var)    # 180 x 43
#                 # per_chain_mean <-  apply(trace_main, c(2, 3), mean)   # 180 x 43
#                 # W <-  colMeans(per_chain_var)                          # within-chain variance
#                 # B <-  apply(per_chain_mean, 2, var)                    # between-chain variance of means
#                 # summary(B / W)
#                 # 
#                 # worst <-  which.max(B/W)   # recompute B/W on this 1500-iter trace first
#                 # hist(per_chain_mean[, worst], breaks = 40)
#                 # ## and the top-3 worst, and specifically beta[2,1,2] and p_raw[1]
#                 # 
#                 # cor(log(ratio), log(B/W))   # strongly negative → metric miscalibration IS the mixing bottleneck
#                 # 
#                 # 
#                 # tibble(data.frame(param = dimnames(trace_main)[[3]], ratio = ratio)) %>% dplyr::arrange(ratio) %>% print(n=100)
#                 # 
#                 # ## pooled draws, whitened by full covariance:
#                 # draws_mat <-  matrix(trace_main, nrow = 500*180, ncol = 43)
#                 # S <-  cov(draws_mat)
#                 # kappa(S, exact = TRUE)                      # condition number: huge → dense M has headroom
#                 # evd <-  eigen(S)
#                 # white <-  scale(draws_mat, center = TRUE, scale = FALSE) %*% evd$vectors %*% diag(1/sqrt(evd$values))
#                 # ## nonlinearity check: pairwise plots of the worst actors in whitened space
#                 # # pairs(white[sample(nrow(white), 2000), c(<indices of beta[2,1,2], p_raw[1], worst Omegas>)])
#                 # 
#                 # ## ---- Whitened pairs plot: worst mixers + worst ratios ----
#                 # 
#                 # param_names <-  dimnames(trace_main)[[3]]
#                 # 
#                 # ## Pooled draws and whitening:
#                 # draws_mat <-  matrix(trace_main, nrow = 500 * 180, ncol = length(param_names))
#                 # colnames(draws_mat) <-  param_names   # BEFORE whitening, so we can track columns
#                 # 
#                 # S <-  cov(draws_mat)
#                 # cat("condition number kappa(S) =", kappa(S, exact = TRUE), "\n")
#                 # 
#                 # evd   <-  eigen(S, symmetric = TRUE)
#                 # white <-  sweep(draws_mat, 2, colMeans(draws_mat)) %*% evd$vectors %*% diag(1 / sqrt(evd$values))
#                 # 
#                 # ## NOTE: whitening mixes all params together - column j of `white` is the j-th
#                 # ## PRINCIPAL DIRECTION, not parameter j. Two complementary views:
#                 # 
#                 # ## ---- View 1: raw-space pairs of the worst actors (interpretable axes) ----
#                 # worst_params <-  c("beta[2,1,2]", "p_raw[1]",
#                 #                   "Omega_unconstrained_vec[2,1]",
#                 #                   "Omega_unconstrained_vec[2,12]",
#                 #                   "Omega_unconstrained_vec[2,7]",
#                 #                   "beta[2,1,3]")
#                 # idx <-  match(worst_params, param_names)
#                 # stopifnot(!anyNA(idx))   # catches any name mismatch (spacing etc.)
#                 # 
#                 # set.seed(1)
#                 # sub <-  sample(nrow(draws_mat), 2000)
#                 # 
#                 # pairs(draws_mat[sub, idx], pch = ".", col = rgb(0, 0, 0, 0.35),
#                 #       labels = worst_params, main = "Raw space: worst mixers")
#                 # 
#                 # ## ---- View 2: whitened principal directions with the longest autocorrelation ----
#                 # ## Find which whitened directions mix worst (per-chain B/W in whitened space):
#                 # white_arr <-  array(white, dim = c(500, 180, ncol(white)))
#                 # pc_mean <-  apply(white_arr, c(2, 3), mean)
#                 # pc_var  <-  apply(white_arr, c(2, 3), var)
#                 # BW_white <-  apply(pc_mean, 2, var) / colMeans(pc_var)
#                 # worst_dirs <-  order(BW_white, decreasing = TRUE)[1:6]
#                 # cat("worst whitened directions:", worst_dirs, "| B/W:", round(BW_white[worst_dirs], 2), "\n")
#                 # 
#                 # pairs(white[sub, worst_dirs], pch = ".", col = rgb(0, 0, 0, 0.35),
#                 #       labels = paste0("PC", worst_dirs),
#                 #       main = "Whitened space: slowest directions")
#                 # 
#                 # ## Which params load on those slow directions (interpretation aid):
#                 # for (d in worst_dirs[1:3]) {
#                 #   loadings <-  evd$vectors[, d]
#                 #   top <-  order(abs(loadings), decreasing = TRUE)[1:5]
#                 #   cat("PC", d, ":", paste0(param_names[top], " (", round(loadings[top], 2), ")", collapse = ", "), "\n")
#                 # }
#                 # 
#                 # 
#                 # ## How wrong is the adapted dense metric, as a matrix?
#                 # S_post <-  cov(draws_mat)                                   # ground truth (linear part)
#                 # M_inv_adapted <-  model_results$burnin_object$EHMC_Metric_as_Rcpp_List$M_inv_dense_main # what the sampler used
#                 # ## generalized eigenvalues: if M_inv were perfect, all ~ equal
#                 # ge <-  eigen(solve(M_inv_adapted, S_post))$values
#                 # summary(Re(ge)); max(Re(ge))/min(Re(ge))                   # spread = residual condition number under adapted M
#                 # 
#                 # 
#                 # 
#                 # S_post_pd <-  BayesMVP:::Rcpp_near_PD(S_post)
#                 # 
#                 # 
#                 # 
#                 # first_half  <-  apply(trace_main[1:750, , worst], 2, mean)
#                 # second_half <-  apply(trace_main[751:1500, , worst], 2, mean)
#                 # plot(first_half, second_half, asp = 1); abline(0, 1)
#                 # ## outlier chains specifically:
#                 # outliers <-  which(per_chain_mean[, worst] < -0.45)
#                 # matplot(trace_main[, outliers, worst], type = "l", lty = 1)
#                 # ## and: are the SAME chains outliers across parameters?
#                 # rank_mat <-  apply(per_chain_mean, 2, rank)
#                 # sort(rowMeans(rank_mat))[1:10]   # consistently extreme chains?
#                 # 
#                 # 
#                 # d <-  second_half - first_half
#                 # mean(d); t.test(d)                 # population drift ≠ 0?
#                 # mean(d < 0)                        # fraction of chains that moved down
#                 # ## same thing across ALL params at once:
#                 # D <-  apply(trace_main[751:1500, , ], c(2,3), mean) - apply(trace_main[1:750, , ], c(2,3), mean)
#                 # colMeans(D < 0)                    # per-param fraction below diagonal
#                 # 
#                 # matplot(trace_main[, sample(180, 15), worst], type = "l", lty = 1)   # random chains
#                 # plot(colMeans(trace_main[, , worst]), type = "l")                    # ensemble mean vs iteration
#                 # 
#                 
#                 
#                 
#                 
#                 ##
#                 max_Rhat  <-  max(tibble_gq$Rhat, na.rm = TRUE)
#                 max_nRhat <-  max(tibble_gq$n_Rhat, na.rm = TRUE)
#                 ##
#                 min_ESS  <-  min(tibble_gq$n_eff, na.rm = TRUE)
#                 ##
#                 time_total  <-  model_fit_object$summaries$efficiency_info$time_total
#                 time_burnin <-  model_fit_object$summaries$efficiency_info$time_burnin
#                 ##
#                 time_sampling <-  model_fit_object$summaries$efficiency_info$time_sampling
#                 time_summaries <-  model_fit_object$summaries$efficiency_info$time_summaries
#                 ##
#                 ##
#                 summary_tibbles <-  model_fit_object$summaries$summary_tibbles
#                 efficiency_info <-  model_fit_object$summaries$efficiency_info
#                 divergences <-  model_fit_object$summaries$divergences
#                 HMC_info <-  model_fit_object$summaries$HMC_info
#                 ##
#                 ESS_per_sec_samp <-  min_ESS / time_sampling
#                 ESS_per_sec_total <-  min_ESS / time_total
#                 ##
#                 # try({ 
#                   L_main_during_sampling <-  (HMC_info$tau_main / HMC_info$eps_main)
#                   n_grad_evals_sampling_main <-  L_main_during_sampling * settings$n_iter * settings$n_chains_sampling
#                   Min_ess_per_grad_main_samp <-  min_ESS / n_grad_evals_sampling_main
#                 # }, silent = TRUE)
#                 # try({
#                   L_us_during_sampling <-  (HMC_info$tau_us / HMC_info$eps_us)
#                   n_grad_evals_sampling_us <-  L_us_during_sampling  * settings$n_iter * settings$n_chains_sampling
#                   Min_ess_per_grad_us_samp <-  min_ESS / n_grad_evals_sampling_us
#                 # }, silent = TRUE)
#                 # if (partitioned_HMC == TRUE) { ## i.e. if nuisance are sampledseperately 
#                 #   if (Model_type == "Stan") {
#                 #     weight_nuisance_grad <-  1.0
#                 #     weight_main_grad <-  1.0 
#                 #   } else { 
#                 #     weight_nuisance_grad <-  0.3333333
#                 #     weight_main_grad <-  0.6666667 ## main grad takes ~ 2x as long to compute as nuisance grad 
#                 #   }
#                 #   Min_ess_per_grad_samp_weighted <-  (weight_nuisance_grad * Min_ess_per_grad_us_samp + weight_main_grad * Min_ess_per_grad_main_samp) /
#                 #     (weight_nuisance_grad + weight_main_grad)
#                 # } else if (partitioned_HMC == FALSE) {  # if not partitioned, grad isnt seperate and "main" grad = "all" grad so use weight_main_grad only !!
#                   weight_nuisance_grad <-  0.00
#                   weight_main_grad <-  1.00
#                   # Min_ess_per_grad_samp_weighted <-  (weight_nuisance_grad * Min_ess_per_grad_us_samp + weight_main_grad * Min_ess_per_grad_main_samp) / 
#                   #                                   (weight_nuisance_grad + weight_main_grad)
#                   # ##
#                   # Min_ess_per_grad_samp_weighted*1000
#                   ##
#                   Min_ess_per_grad_samp_weighted <-  Min_ess_per_grad_main_samp
#                 # }
#                 ##
#                 print(paste("ESS/sec = ", round(ESS_per_sec_samp, 3)))
#                 print(paste("ESS/grad = ", round(Min_ess_per_grad_samp_weighted*1000, 3)))
#                 ##
#                 # ##
#                 # if (N <= 2500) {
#                 #   print(paste("ESS/grad = ", Min_ess_per_grad_samp_weighted*1000))
#                 #   if ((Min_ess_per_grad_samp_weighted*1000) < 0.25) {
#                 #     stop("ESS/grad too low")
#                 #   }
#                 # }
#                 # if (N == 10000) {
#                 #   print(paste("ESS/grad = ", Min_ess_per_grad_samp_weighted*1000))
#                 #   if ((Min_ess_per_grad_samp_weighted*1000) < 0.025) {
#                 #     stop("ESS/grad too low")
#                 #   }
#                 # }
#                 # ##
#                 cat(sprintf("  min_ESS = %d, max_Rhat = %.4f\n", round(min_ESS), max_Rhat))
#                 ##
#                 all_runs_results[[run_i]] <-  list(
#                   run_id = run_i,
#                   ##
#                   MCMC_seed = MCMC_seed,
#                   ##
#                   Se_mean = tibble_Se$mean,
#                   Sp_mean = tibble_Sp$mean,
#                   prev_mean = tibble_prev$mean,
#                   ##
#                   tibble_main = tibble_main,
#                   tibble_tp = tibble_tp,
#                   tibble_gq = tibble_gq,
#                   ##
#                   trace_main = trace_main,
#                   trace_tp = trace_tp,
#                   trace_gq = trace_gq,
#                   ##
#                   max_Rhat = max_Rhat,
#                   max_nRhat = max_nRhat,
#                   ##
#                   min_ESS = min_ESS,
#                   ##
#                   ESS_per_sec_total = ESS_per_sec_total,
#                   ESS_per_sec_samp = ESS_per_sec_samp,
#                   ESS_per_grad_samp = Min_ess_per_grad_samp_weighted,
#                   ##
#                   n_chains = settings$n_chains,
#                   n_iter = settings$n_iter,
#                   ##
#                   summary_tibbles = summary_tibbles,
#                   efficiency_info = efficiency_info,
#                   divergences = divergences,
#                   HMC_info = HMC_info,
#                   ##
#                   time_total = time_total,
#                   time_burnin = time_burnin
#                 )
#                 ##
#                 gc(verbose = FALSE)
#                 ##
#                 ## ---- Save raw results:
#                 ##
#                 file_name_string <-  R_fn_file_name_string( Model_type = Model_type,
#                                                            N = N,
#                                                            settings = settings,
#                                                            model_args_list = model_args_list,
#                                                            prior_LKJ_nd = prior_LKJ_nd,
#                                                            prior_LKJ_d = prior_LKJ_d,
#                                                            prior_prev_a = prior_prev_a,
#                                                            prior_prev_b = prior_prev_b,
#                                                            grouping = grouping)
#                 ##
#                 run_file <-  paste0(file_name_string, "_run", run_i)
#                 ##
#                 saveRDS(all_runs_results[[run_i]], run_file)
#                 ##
#                 # beepr::beep(sound = 6)
#         }
#         
#         {
#           # for (run_i in 1:settings$n_runs) {
#           #     all_runs_results[[run_i]]$trace_gq <-  NULL
#           # }
#           # ##
#           # ## Save raw results:
#           # ##
#           # saveRDS(all_runs_results, 
#           #         file.path( settings$output_dir, paste0("ps7_raw_runs",
#           #                                                "_N", N,
#           #                                                "_n_ch_burn", settings$n_chains_burnin,
#           #                                                "_n_ch_samp", settings$n_chains_sampling,
#           #                                                "_n_burn", settings$n_burnin,
#           #                                                "_n_iter", settings$n_iter,
#           #                                                "_LR", settings$learning_rate,
#           #                                                "_AD", settings$adapt_delta,
#           #                                                "_clip", settings$clip_iter,
#           #                                                "_int", settings$int,
#           #                                                "_ratio_M", settings$ratio_M_main, "_", settings$ratio_M_nuisance,
#           #                                                "_M_type", settings$metric_type_main,
#           #                                                "_M_shape", settings$metric_shape_main,
#           #                                                "_width", settings$int_width
#           #                                                )))
#                                                         
#                               
#         }
#         
#         return(all_runs_results)
#         
# }
# 
# 
# 
# 
# R_fn_file_name_string <-  function(Model_type, 
#                                   N, 
#                                   settings, 
#                                   model_args_list,
#                                   prior_LKJ_nd, 
#                                   prior_LKJ_d,
#                                   prior_prev_a, 
#                                   prior_prev_b,
#                                   grouping = NULL
# ) {
#   
#             file_name_string <-  file.path(settings$output_dir, paste0("ps7_run",
#                                                                           ##
#                                                                           "_", Model_type,
#                                                                           ##
#                                                                           "_N", N,
#                                                                           ##
#                                                                           "_n_pops", model_args_list$n_pops, ##
#                                                                           ##
#                                                                           "_diff_HMC", settings$diffusion_HMC, 
#                                                                           "_pt_HMC", settings$partitioned_HMC, 
#                                                                           ##
#                                                                           "_n_ch_b", settings$n_chains_burnin,
#                                                                           "_n_ch_s", settings$n_chains_sampling,
#                                                                           ##
#                                                                           "_n_WCP_b", settings$n_threads_WCP_burnin,
#                                                                           "_n_WCP_s", settings$n_threads_WCP_sampling,
#                                                                           ##
#                                                                           "_n_chnk_b", model_args_list$num_chunks,
#                                                                           "_n_chnk_s", model_args_list$num_chunks,
#                                                                           ##
#                                                                           "_n_b", settings$n_burnin,
#                                                                           "_n_i", settings$n_iter,
#                                                                           "_LR", settings$learning_rate,
#                                                                           "_AD", settings$adapt_delta,
#                                                                           "_clip", settings$clip_iter,
#                                                                           "_int", settings$int,
#                                                                           "_ratio_M", settings$ratio_M_main, "_", settings$ratio_M_nuisance,
#                                                                           "_M_typ", settings$metric_type_main,
#                                                                           "_M_shp", settings$metric_shape_main,
#                                                                           "_width", settings$int_width,
#                                                                           ##
#                                                                           "_M_dcy", settings$M_decay_type,
#                                                                           ##
#                                                                           # "true_p", true_prev,
#                                                                           ##
#                                                                           "pi_LKJ", prior_LKJ_nd, "_", prior_LKJ_d,
#                                                                           "pi_p", prior_prev_a, "_", prior_prev_b
#             ))
#             ##
#             if (settings$M_decay_type == "inverse") { 
#               file_name_string <-  paste0(file_name_string,        
#                                          "_M_pow", settings$M_decay_power, ## new
#                                          "_M_scl", settings$M_decay_scale ## new
#               )
#             }
#             if (grouping == TRUE) { 
#               file_name_string <-  paste0(file_name_string,        
#                                          "_ord_grp", grouping)
#             }
#             
#             return(file_name_string)
# 
# }
# 
# 
# ##
# ## ---- Function to summarize BayesMVP pilot study 7 results: -----------------------------------------------------------------------------------------
# ##
# summarize_ps7_results <-  function(output_dir,
#                                   ##
#                                   Model_type,
#                                   ##
#                                   N_val,
#                                   ##
#                                   n_chains_burnin = 8,
#                                   n_chains_sampling = 64,
#                                   ##
#                                   n_iter = 1000,
#                                   ##
#                                   adapt_delta = 0.80,
#                                   ##
#                                   n_burnin_vec = c(250, 500, 1000),
#                                   clip_iter_vec = c(50, 100, 100),
#                                   int_vec = c(75, 100, 300),
#                                   ##
#                                   ratio_M_main_vec,
#                                   ratio_M_nuisance_vec,
#                                   ##
#                                   int_width_vec = c(1, 1, 1),
#                                   ##
#                                   learning_rate_vec = c(0.0125, 0.025, 0.0375, 0.05, 0.0625, 0.075),
#                                   ##
#                                   metric_type_main,
#                                   metric_shape_main,
#                                   ##
#                                   diffusion_HMC,
#                                   partitioned_HMC,
#                                   ##
#                                   M_decay_type,
#                                   M_decay_power = NULL,
#                                   # M_decay_scale = NULL,
#                                   ##
#                                   n_runs = 10,
#                                   ##
#                                   prior_LKJ_nd,
#                                   prior_LKJ_d,
#                                   prior_prev_a,
#                                   prior_prev_b,
#                                   ##
#                                   grouping,
#                                   ##
#                                   n_threads_WCP_burnin,
#                                   n_threads_WCP_sampling,
#                                   ##
#                                   num_chunks_burnin,
#                                   num_chunks_sampling
#                                   
# ) {
#   
#         require(dplyr)
#         require(tibble)
#         
#         ## Storage for all results
#         all_rows <-  list()
#         row_idx <-  1
#         
#         ## Loop through all configurations using INDEX
#         for (ii in seq_along(n_burnin_vec)) {
#           
#           n_burnin <-  n_burnin_vec[ii]
#           clip_iter <-  clip_iter_vec[ii]
#           int <-  int_vec[ii]
#           ratio_M_main <-  ratio_M_main_vec[ii]
#           ratio_M_nuisance <-  ratio_M_nuisance_vec[ii]
#           int_width <-  int_width_vec[ii]
#           
#           if (n_burnin == 1000) { 
#             n_adapt <-  900
#           } else { 
#             n_adapt <-  n_burnin - round(n_burnin/10)
#           }
#           ##
#           # M_decay_scale <-  n_adapt / 5 ## -----------------
#           M_decay_scale <-  n_adapt / 100 ## -----------------
#           
#           try({
#             
#             for (learning_rate in learning_rate_vec) {
#               
#               # cat(sprintf("Trying: ii=%d, n_burnin=%d, clip=%d, int=%d, ratio_M=%s, ratio_M=%s, width=%d, LR=%s\n", 
#               #             ii, n_burnin, clip_iter, int, ratio_M_main, ratio_M_nuisance, int_width, learning_rate))
#               
#               ## Build filename
#               file_name_string <-  file.path(output_dir, paste0(
#                 "ps7_run",
#                 ##
#                 "_", Model_type,
#                 ##
#                 "_N", N_val,
#                 ##
#                 "_diff_HMC", diffusion_HMC, 
#                 "_pt_HMC", partitioned_HMC, 
#                 ##
#                 "_n_ch_b", n_chains_burnin,
#                 "_n_ch_s", n_chains_sampling,
#                 ##
#                 "_n_WCP_b", n_threads_WCP_burnin,
#                 "_n_WCP_s", n_threads_WCP_sampling,
#                 ##
#                 "_n_chnk_b", num_chunks_burnin,
#                 "_n_chnk_s", num_chunks_sampling,
#                 ##
#                 "_n_b", n_burnin,
#                 "_n_i", n_iter,
#                 "_LR", learning_rate,
#                 "_AD", adapt_delta,
#                 "_clip", clip_iter,
#                 "_int", int,
#                 "_ratio_M", ratio_M_main, "_", ratio_M_nuisance,
#                 "_M_typ", metric_type_main,
#                 "_M_shp", metric_shape_main,
#                 "_width", int_width,
#                 ##
#                 "_M_dcy", M_decay_type,
#                 ##
#                 "pi_LKJ", prior_LKJ_nd, "_", prior_LKJ_d,
#                 "pi_p", prior_prev_a, "_", prior_prev_b
#               ))
#               ##
#               if (M_decay_type == "inverse") { 
#                 file_name_string <-  paste0(file_name_string,        
#                                            "_M_pow", M_decay_power, ## new
#                                            "_M_scl", M_decay_scale ## new
#                 )
#               }
#               ##
#               if (grouping == TRUE) { 
#                 file_name_string <-  paste0( file_name_string,        
#                                          "_ord_grp", grouping)
#               }
#               print(paste(file_name_string))
#               ##
#               if (file.exists(file_name_string)) {
#                 all_runs <-  readRDS(file_name_string)
#               } else {
#                 ## Load individual run files
#                 all_runs <-  list()
#                 for (run_i in 1:n_runs) {
#                   run_file <-  paste0(file_name_string, "_run", run_i)
#                   if (file.exists(run_file)) {
#                     all_runs[[run_i]] <-  readRDS(run_file)
#                   }
#                 }
#               }
#               
#               if (length(all_runs) == 0) next
#               
#               for (run_i in seq_along(all_runs)) {
#                 
#                 try({
#                   
#                   r <-  all_runs[[run_i]]
#                   if (is.null(r)) next
#                   
#                   eff <-  r$efficiency_info
#                   hmc <-  r$HMC_info
#                   divs <-  r$divergences
#                   
#                   row <-  tibble(
#                     N = N_val,
#                     n_burnin = n_burnin,
#                     clip_iter = clip_iter,
#                     int = int,
#                     ##
#                     ratio_M_main = ratio_M_main,
#                     ratio_M_nuisance = ratio_M_nuisance,
#                     ##
#                     int_width = int_width,
#                     learning_rate = learning_rate,
#                     run = run_i,
#                     n_chains_burnin = hmc$n_chains_burnin,
#                     n_chains_sampling = hmc$n_chains_sampling,
#                     n_iter = eff$n_iter,
#                     ##
#                     time_total = eff$time_total,
#                     time_sampling = eff$time_sampling,
#                     time_burnin = eff$time_burnin,
#                     time_summaries = eff$time_summaries,
#                     time_total_wo_summaries = eff$time_total_wo_summaries,
#                     ##
#                     max_Rhat = eff$Max_rhat_main,
#                     max_nRhat = eff$Max_nested_rhat_main,
#                     ##
#                     min_ESS = r$min_ESS, ## eff$Min_ESS_main,
#                     ##
#                     ESS_per_sec_total = r$ESS_per_sec_total, ## eff$ESS_per_sec_total,
#                     ESS_per_sec_samp = r$ESS_per_sec_samp, ## eff$ESS_per_sec_samp,
#                     ##
#                     n_divs = divs$n_divs,
#                     pct_divs = divs$pct_divs,
#                     ##
#                     ESS_per_grad_samp = r$ESS_per_grad_samp*1000,  ## 1000 * eff$Min_ess_per_grad_samp_weighted,
#                     grad_evals_per_sec = eff$grad_evals_per_sec,
#                     ##
#                     L_main_burnin = eff$L_main_during_burnin,
#                     L_main_samp = eff$L_main_during_sampling,
#                     L_us_burnin = eff$L_us_during_burnin,
#                     L_us_samp = eff$L_us_during_sampling,
#                     ##
#                     tau_main = hmc$tau_main,
#                     eps_main = hmc$eps_main,
#                     tau_us = hmc$tau_us,
#                     eps_us = hmc$eps_us,
#                     ##
#                     adapt_delta = hmc$adapt_delta,
#                     ##
#                     sampling_time_to_Min_ESS = eff$sampling_time_to_Min_ESS,
#                     total_time_to_1000_ESS_wo_summaries = eff$total_time_to_1000_ESS_wo_summaries,
#                     total_time_to_10000_ESS_wo_summaries = eff$total_time_to_10000_ESS_wo_summaries
#                   )
#                   
#                   all_rows[[row_idx]] <-  row
#                   row_idx <-  row_idx + 1
#                   
#                 })
#               }
#             }
#           })
#         }
#         
#         ## Combine all rows
#         if (length(all_rows) > 0) {
#           results_df <-  bind_rows(all_rows)
#           results_df <-  results_df %>% arrange(N, n_burnin, clip_iter, int, ratio_M_main, ratio_M_nuisance, int_width, learning_rate, run)
#         } else {
#           results_df <-  tibble()
#         }
#         
#         return(results_df)
#         
# }
# 
# 
# 
# 
# ##
# ## ---- Function to summarize specific parameter combinations: -----------------------------------------------------------------------------------------
# ##
summarize_combos <-  function(df,
                             combos_df,
                             group_vars = c("clip_iter", "int", "ratio_M_main", "ratio_M_nuisance", "int_width"),
                             summary_vars = c("ESS_per_grad_samp", "ESS_per_sec_samp", "min_ESS", "pct_divs")) {

        require(dplyr)

        ## Filter df to only rows matching combos_df
        df_filtered <-  df %>%
          dplyr::inner_join(unique(x = combos_df), by = names(combos_df))

        ## Keep the sweep axes separate even when an older caller supplies only the original grouping fields.
        group_vars <-  unique(x = c(group_vars, intersect(x = c("N", "n_burnin", "learning_rate", "metric_estimator", "tau_initial",
                                                               "learning_rate_initial", "diffusion_HMC_integrator", "target_ESS",
                                                               ## see summarize_ps7_by_config: never average across tau schemes
                                                               "tau_scheme"),
                                                        y = names(x = df_filtered))))

        ## Summarize by the grouping variables
        summary_df <-  df_filtered %>%
          dplyr::group_by(across(all_of(group_vars))) %>%
          dplyr::summarise(
            n = dplyr::n(),
            dplyr::across(all_of(summary_vars),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        # median = ~median(.x, na.rm = TRUE),
                        # g_mean = ~exp(mean(log(.x), na.rm = TRUE)),
                        sd = ~sd(.x, na.rm = TRUE),
                        min = ~min(.x, na.rm = TRUE),
                        max = ~max(.x, na.rm = TRUE))),
            .groups = "drop"
          )

        return(summary_df)

}
# ##
# ## ---- Helper to create combos_df easily: -------------------------------------------------------------------------------------------------------------
# ##
# make_combos <-  function(...) {
#   
#         ## Takes named vectors and creates all combinations, OR
#         ## takes a list of specific combos
#         args <-  list(...)
#         
#         ## If first arg is a list of lists, treat as explicit combos
#         if (length(args) == 1 && is.list(args[[1]]) && is.list(args[[1]][[1]])) {
#           return(bind_rows(args[[1]]))
#         }
#         
#         ## Otherwise expand grid
#         expand.grid(args, stringsAsFactors = FALSE) %>% as_tibble()
#   
# }
# 
# 
# 
# ##
# ## ---- Function to create aggregated summary by configuration: ----------------------------------------------------------------------------------------
# ##
summarize_ps7_by_config <-  function(results_df) {

        if (nrow(x = results_df) == 0) return(results_df)
        ## Paired burn-in settings and all sweep axes identify a configuration; average only its replicate runs.
        configuration_columns <-  c("N", "n_burnin", "learning_rate", "metric_estimator", "tau_initial",
                                   "learning_rate_initial", "diffusion_HMC_integrator", "learning_rate_initial_iter",
                                   "clip_iter", "int", "ratio_M_main", "ratio_M_nuisance", "int_width",
                                   "n_chains_burnin", "n_chains_sampling", "n_iter", "adapt_delta", "target_ESS",
                                   ##
                                   ## Trajectory-length adaptation. These MUST be here: two runs that differ only in
                                   ## the tau objective are different configurations, and averaging them together
                                   ## would silently hide the very comparison they were run for. any_of() below means
                                   ## a results_df that predates these columns is unaffected.
                                   "burnin_algorithm", "tau_weight_by_p_jump", "manual_L", "manual_tau_value", "tau_scheme",
                                   "theta_hat_us_rule", "theta_hat_us_freeze_iter", "store_log_lik_trace")
        agg_df <-  results_df %>%
          dplyr::group_by(dplyr::across(dplyr::any_of(configuration_columns))) %>%
          dplyr::summarise(
            n_runs = length(N),
            ##
            time_total_mean = mean(time_total, na.rm = TRUE),
            time_total_sd = sd(time_total, na.rm = TRUE),
            ##
            time_sampling_mean = mean(time_sampling, na.rm = TRUE),
            time_sampling_sd = sd(time_sampling, na.rm = TRUE),
            ##
            time_burnin_mean = mean(time_burnin, na.rm = TRUE),
            time_burnin_sd = sd(time_burnin, na.rm = TRUE),
            ##
            time_summaries_mean = mean(time_summaries, na.rm = TRUE),
            time_summaries_sd = sd(time_summaries, na.rm = TRUE),
            ##
            max_Rhat_mean = mean(max_Rhat, na.rm = TRUE),
            max_nRhat_mean = mean(max_nRhat, na.rm = TRUE),
            ##
            min_ESS_mean = mean(min_ESS, na.rm = TRUE),
            min_ESS_sd = sd(min_ESS, na.rm = TRUE),
            ##
            ESS_per_sec_total_mean = mean(ESS_per_sec_total, na.rm = TRUE),
            ESS_per_sec_samp_mean = mean(ESS_per_sec_samp, na.rm = TRUE),
            ##
            n_divs_mean = mean(n_divs, na.rm = TRUE),
            pct_divs_mean = mean(pct_divs, na.rm = TRUE),
            ##
            ESS_per_grad_samp_mean = mean(ESS_per_grad_samp, na.rm = TRUE),
            ## Explicit units for the comparison tables; legacy ESS_per_grad_samp remains per 1000 gradients.
            min_ESS_per_grad_sampling_mean = mean(min_ESS_per_grad_sampling, na.rm = TRUE),
            min_ESS_per_grad_sampling_sd = sd(min_ESS_per_grad_sampling, na.rm = TRUE),
            min_ESS_per_sec_sampling_mean = mean(min_ESS_per_sec_sampling, na.rm = TRUE),
            min_ESS_per_sec_sampling_sd = sd(min_ESS_per_sec_sampling, na.rm = TRUE),
            sampling_time_to_target_ESS_mins_mean = mean(sampling_time_to_target_ESS_mins, na.rm = TRUE),
            est_time_to_target_ESS_wo_summaries_mins_mean = mean(est_time_to_target_ESS_wo_summaries_mins, na.rm = TRUE),
            est_time_to_target_ESS_wo_summaries_mins_sd = sd(est_time_to_target_ESS_wo_summaries_mins, na.rm = TRUE),
            est_time_to_target_ESS_mins_mean = mean(est_time_to_target_ESS_mins, na.rm = TRUE),
            est_time_to_target_ESS_mins_sd = sd(est_time_to_target_ESS_mins, na.rm = TRUE),
            grad_evals_per_sec_mean = mean(grad_evals_per_sec, na.rm = TRUE),
            ##
            L_main_burnin_mean = mean(L_main_burnin, na.rm = TRUE),
            L_main_samp_mean = mean(L_main_samp, na.rm = TRUE),
            ##
            tau_main_mean = mean(tau_main, na.rm = TRUE),
            eps_main_mean = mean(eps_main, na.rm = TRUE),
            ##
            total_time_to_1000_ESS_mean = mean(total_time_to_1000_ESS_wo_summaries, na.rm = TRUE),
            total_time_to_10000_ESS_mean = mean(total_time_to_10000_ESS_wo_summaries, na.rm = TRUE),
            ##
            .groups = "drop"
          ) %>%
          dplyr::arrange(N, n_burnin, learning_rate)

        return(agg_df)

}
##
## -| --------- Reading a ps7 sweep quickly: fn_ps7_report() ---------------------------------------------------------------------------
##
## A 64 x 42 tibble carries the information but hides the answer. This prints the SAME information
## in four blocks, none of which repeats anything:
##
##   1. CONSTANTS - the settings identical in every row (N, n_burnin, clip_iter, ... ), once.
##   2. RANKING   - the configurations ordered by one metric, with short labels, top_n of them.
##   3. MARGINALS - per swept axis, what each level does on average with everything else varying.
##   4. PAIRED    - for each TWO-level axis, hold every other axis fixed and count how often each
##                  level wins. This is the decisive view when axes interact: a level can lose on
##                  the marginal mean (dragged down by one bad partner setting) yet win nearly
##                  every matched pair, or vice versa.
##
## Rows whose metrics are not finite (e.g. ESS/grad = Inf when a run reports zero gradients) are
## EXCLUDED from ranking, marginals and pairs, and listed on their own - never silently ranked first.
##
## Everything is taken from the data frame passed in; nothing is read from the global environment.
## Works with whichever metric columns are present, so it accepts the output of either
## summarize_ps7_by_config() or summarize_ps7_results() + grouping.
##
#' @param summary_df  one row per configuration (from summarize_ps7_by_config()).
#' @param sweep_axes  the setting columns that vary; NULL = detect them (non-metric, >1 value).
#' @param rank_by     metric column to rank on; NULL = first available time-to-target column.
#' @param top_n       how many configurations to print in the ranking (Inf = all).
#' @param digits      printing precision.
#' @return invisibly, a list(constants, ranking, marginals, paired, dropped).
#' @export
fn_ps7_report <-  function( summary_df,
                           sweep_axes = NULL,
                           rank_by = NULL,
                           top_n = 12,
                           digits = 3,
                           max_Rhat_ok = 1.05,      ## classic split-Rhat
                           max_nRhat_ok = 1.05,     ## nested Rhat (the meaningful one for many short chains)
                           max_pct_divs_ok = 1.0    ## per cent divergent transitions
) {

    stopifnot(is.data.frame(x = summary_df))
    summary_df <-  tibble::as_tibble(x = summary_df)
    if (nrow(x = summary_df) == 0) { cat("fn_ps7_report: nothing to report (0 rows).\n"); return(invisible(NULL)) }
    ##
    ## ---- which columns are RESULTS and which are SETTINGS
    metric_columns <-  grep( pattern = "_(mean|sd|min|max|median|g_mean)$", x = names(x = summary_df), value = TRUE)
    bookkeeping <-  intersect( x = c("n", "n_runs", "n_missing", "target_ESS", "run"), y = names(x = summary_df))
    setting_columns <-  setdiff( x = names(x = summary_df), y = c(metric_columns, bookkeeping))
    ##
    n_levels <-  vapply( X = summary_df[setting_columns], FUN = function(column) length(x = unique(x = column)), FUN.VALUE = 0)
    constant_columns <-  setting_columns[n_levels == 1]
    if (is.null(sweep_axes)) sweep_axes <-  setting_columns[n_levels > 1]
    sweep_axes <-  intersect( x = sweep_axes, y = names(x = summary_df))
    ##
    ## ---- the metrics worth showing, in the order they should be read. direction: -1 = lower is better.
    metric_catalogue <-  list( list( column = "est_time_to_target_ESS_wo_summaries_mins_mean", label = "mins_to_ESS", direction = -1),
                               list( column = "est_time_to_target_ESS_mins_mean",             label = "mins_to_ESS", direction = -1),
                               list( column = "min_ESS_per_grad_sampling_mean",               label = "ESS_per_grad", direction = +1),
                               list( column = "ESS_per_grad_samp_mean",                       label = "ESS_per_grad", direction = +1),
                               list( column = "min_ESS_per_sec_sampling_mean",                label = "ESS_per_sec", direction = +1),
                               list( column = "ESS_per_sec_samp_mean",                        label = "ESS_per_sec", direction = +1),
                               list( column = "L_main_samp_mean",                             label = "L_main", direction = 0),
                               list( column = "max_Rhat_mean",                                label = "max_Rhat", direction = -1),
                               list( column = "pct_divs_mean",                                label = "pct_divs", direction = -1))
    metrics <-  Filter( f = function(m) m$column %in% names(x = summary_df), x = metric_catalogue)
    ##
    seen_labels <-  character()   ## keep only the FIRST available column per label (the two summarisers differ)
    metrics <-  Filter( f = function(m) { keep <-  !(m$label %in% seen_labels); seen_labels <<- c(seen_labels, m$label); keep }, x = metrics)
    if (length(x = metrics) == 0) stop("fn_ps7_report: none of the expected metric columns are present.")
    ##
    metric_columns_used <-  vapply( X = metrics, FUN = function(m) m$column, FUN.VALUE = "")
    metric_labels <-  vapply( X = metrics, FUN = function(m) m$label, FUN.VALUE = "")
    directions <-  vapply( X = metrics, FUN = function(m) m$direction, FUN.VALUE = 0)
    ##
    if (is.null(rank_by)) rank_by <-  metric_columns_used[[1]]
    if (!rank_by %in% names(x = summary_df)) stop("fn_ps7_report: rank_by column not present: ", rank_by)
    rank_direction <-  if (rank_by %in% metric_columns_used) directions[[match(x = rank_by, table = metric_columns_used)]] else -1
    ##
    ## ---- EXCLUDE configurations that did not converge, and say why.
    ## A stuck run reports every n_eff / Rhat as NaN; min()/max() over the survivors then give
    ## min_ESS = +Inf and max_Rhat = -Inf, so the broken run looks infinitely fast and would rank
    ## FIRST on every speed metric. Non-finite metrics and bad Rhat / divergences are therefore
    ## excluded from the ranking, the marginals AND the pairs, and listed separately with a reason.
    fn_fails <-  function( column, limit) {
        if (!column %in% names(x = summary_df)) return(rep( x = FALSE, times = nrow(x = summary_df)))
        value <-  summary_df[[column]]
        is.finite(value) & value > limit
    }
    finite_row <-  Reduce( f = `&`, x = lapply( X = summary_df[metric_columns_used], FUN = is.finite))
    fails_rhat <-  fn_fails( column = "max_Rhat_mean",  limit = max_Rhat_ok) |
                   fn_fails( column = "max_nRhat_mean", limit = max_nRhat_ok)
    fails_divs <-  fn_fails( column = "pct_divs_mean",  limit = max_pct_divs_ok)
    ##
    reason <-  ifelse( test = !finite_row, yes = "non-finite metric",
                       no = ifelse( test = fails_rhat, yes = "Rhat too high",
                                    no = ifelse( test = fails_divs, yes = "too many divergences", no = "")))
    keep <-  finite_row & !fails_rhat & !fails_divs
    dropped <-  summary_df[!keep, , drop = FALSE]
    dropped_reason <-  reason[!keep]
    clean_df <-  summary_df[keep, , drop = FALSE]
    if (nrow(x = clean_df) == 0) {
        cat("fn_ps7_report: every configuration was excluded (", paste( unique(x = dropped_reason), collapse = "; "), ").\n", sep = "")
        return(invisible(x = list( dropped = dropped, dropped_reason = dropped_reason)))
    }
    ##
    ## ---- short labels, so the ranking fits on one line
    abbreviations <-  c( chain_mean = "cm", pooled = "pool", kick_flow_kick = "kfk", flow_kick_flow = "fkf",
                         Empirical = "emp", Hessian = "hess", uniform_diag = "unif_d", unit = "unit",
                         dense = "dense", diag = "diag", inverse = "inv", exponential = "exp", constant = "con")
    fn_short <-  function( values, column_name) {
        if (is.character(x = values) || is.factor(x = values)) {
            out <-  unname(obj = abbreviations[as.character(x = values)])
            return(ifelse( test = is.na(x = out), yes = as.character(x = values), no = out))
        }
        if (column_name == "tau_initial") {
            return(ifelse( test = abs(values - pi) < 1e-8, yes = "pi",
                           no = ifelse( test = abs(values - 2*pi) < 1e-8, yes = "2pi", no = format(x = values, digits = 4))))
        }
        format(x = values, trim = TRUE)
    }
    ##
    ## ---- 1. constants
    cat("\n================ ps7 sweep:", nrow(x = summary_df), "configurations,", length(x = sweep_axes), "axes ================\n")
    if (length(x = constant_columns) > 0) {
        cat("\n-- same in every configuration --\n")
        cat(paste0("   ", paste( constant_columns, "=", vapply( X = constant_columns,
                                                                FUN = function(column) as.character(x = summary_df[[column]][1]),
                                                                FUN.VALUE = ""), collapse = " | ")), "\n")
    }
    for (column in bookkeeping) {
        if (length(x = unique(x = summary_df[[column]])) == 1) cat("   ", column, "=", summary_df[[column]][1], "\n")
    }
    ##
    ## ---- 2. ranking
    ranking <-  clean_df
    ranking <-  dplyr::arrange(.data = ranking, dplyr::desc(.env$rank_direction * .data[[rank_by]]))
    ranking_out <-  tibble::tibble( rank = seq_len(length.out = nrow(x = ranking)))
    for (axis in sweep_axes) ranking_out[[axis]] <-  fn_short( values = ranking[[axis]], column_name = axis)
    for (index in seq_along(along.with = metric_columns_used)) {
        ranking_out[[metric_labels[index]]] <-  round( x = ranking[[metric_columns_used[index]]], digits = digits)
    }
    cat("\n-- ranked by", rank_by, if (rank_direction < 0) "(lower is better)" else "(higher is better)", "--\n")
    print( x = utils::head(x = ranking_out, n = min(top_n, nrow(x = ranking_out))), n = Inf, width = Inf)
    if (nrow(x = ranking_out) > top_n) cat("   ... ", nrow(x = ranking_out) - top_n, " more (top_n = ", top_n, ")\n", sep = "")
    ##
    ## ---- 3. marginals
    cat("\n-- marginal means per axis (every other setting varying) --\n")
    marginals <-  list()
    for (axis in sweep_axes) {
        levels_here <-  unique(x = clean_df[[axis]])
        rows <-  lapply( X = levels_here, FUN = function(level) {
            subset_df <-  clean_df[clean_df[[axis]] == level, , drop = FALSE]
            row <-  tibble::tibble( axis = axis, level = fn_short( values = level, column_name = axis), n_configs = nrow(x = subset_df))
            for (index in seq_along(along.with = metric_columns_used)) {
                row[[metric_labels[index]]] <-  round( x = mean(subset_df[[metric_columns_used[index]]]), digits = digits)
            }
            row
        })
        marginals[[axis]] <-  dplyr::bind_rows(rows)
    }
    marginal_table <-  dplyr::bind_rows(marginals)
    print( x = marginal_table, n = Inf, width = Inf)
    ##
    ## ---- 4. paired comparison for every two-level axis
    paired <-  list()
    two_level_axes <-  sweep_axes[vapply( X = sweep_axes, FUN = function(axis) length(x = unique(x = clean_df[[axis]])) == 2, FUN.VALUE = TRUE)]
    if (length(x = two_level_axes) > 0) {
        cat("\n-- paired comparisons (all other axes held fixed; one matched pair per combination) --\n")
        for (axis in two_level_axes) {
            other_axes <-  setdiff( x = sweep_axes, y = axis)
            key <-  if (length(x = other_axes) == 0) rep( x = "all", times = nrow(x = clean_df)) else
                        do.call( what = paste, args = c(clean_df[other_axes], sep = "|"))
            levels_here <-  sort(x = unique(x = clean_df[[axis]]))
            rows <-  lapply( X = seq_along(along.with = metric_columns_used), FUN = function(index) {
                column <-  metric_columns_used[index]
                direction <-  directions[index]
                wins_a <-  0; wins_b <-  0; pairs <-  0; ratios <-  numeric(0)
                for (group in unique(x = key)) {
                    pair_df <-  clean_df[key == group, , drop = FALSE]
                    value_a <-  pair_df[[column]][pair_df[[axis]] == levels_here[1]]
                    value_b <-  pair_df[[column]][pair_df[[axis]] == levels_here[2]]
                    if (length(x = value_a) != 1 || length(x = value_b) != 1) next
                    pairs <-  pairs + 1
                    if (direction != 0) {
                        if (direction * (value_a - value_b) > 0) wins_a <-  wins_a + 1 else if (value_a != value_b) wins_b <-  wins_b + 1
                    }
                    if (value_b != 0) ratios <-  c(ratios, value_a / value_b)
                }
                tibble::tibble( axis = axis,
                            metric = metric_labels[index],
                            pairs = pairs,
                            wins_A = if (direction == 0) NA_integer_ else wins_a,
                            wins_B = if (direction == 0) NA_integer_ else wins_b,
                            A = fn_short( values = levels_here[1], column_name = axis),
                            B = fn_short( values = levels_here[2], column_name = axis),
                            median_A_over_B = round( x = stats::median(x = ratios), digits = digits))
            })
            paired[[axis]] <-  dplyr::bind_rows(rows)
        }
        paired_table <-  dplyr::bind_rows(paired)
        print( x = paired_table, n = Inf, width = Inf)
        cat("   wins_A / wins_B count matched pairs where that level is better on that metric;",
            "median_A_over_B is the typical ratio A:B.\n")
    } else {
        paired_table <-  NULL
    }
    ##
    ## ---- excluded rows
    if (nrow(x = dropped) > 0) {
        cat("\n-- EXCLUDED: ", nrow(x = dropped), " configuration(s) that did not converge (not ranked, not averaged, not paired) --\n", sep = "")
        excluded_columns <-  c( lapply( X = stats::setNames( object = sweep_axes, nm = sweep_axes),
                                        FUN = function(axis) fn_short( values = dropped[[axis]], column_name = axis)),
                                list( why = dropped_reason),
                                stats::setNames( object = lapply( X = metric_columns_used, FUN = function(column) dropped[[column]]),
                                                 nm = metric_labels))
        print( x = tibble::as_tibble(x = excluded_columns), n = Inf, width = Inf)
    }
    cat("\n")
    ##
    invisible(x = list( constants = constant_columns, ranking = ranking_out, marginals = marginal_table,
                        paired = paired_table, dropped = dropped, rank_by = rank_by))

}
# # ##
# # ## ---- Function to print a nice comparison table: -----------------------------------------------------------------------------------------
# # ##
# # print_ps7_comparison <-  function(agg_df) {
# #   
# #         cat("\n========== PS7 Configuration Comparison ==========\n\n")
# #         
# #         comparison_df <-  agg_df %>%
# #           mutate(
# #             config = paste0("n_burn=", n_burnin, ", LR=", learning_rate)
# #           ) %>%
# #           select(
# #             N, config, n_runs,
# #             min_ESS_mean, 
# #             ESS_per_sec_samp_mean,
# #             ESS_per_grad_samp_mean,
# #             pct_divs_mean,
# #             time_total_mean,
# #             L_main_samp_mean,
# #             eps_main_mean
# #           )
# #         
# #         print(comparison_df, n = Inf)
# #         
# #         ## Find best config by ESS/sec
# #         best_by_ESS_sec <-  agg_df %>%
# #           group_by(N) %>%
# #           slice_max(ESS_per_sec_samp_mean, n = 1) %>%
# #           mutate(config = paste0("n_burn=", n_burnin, ", LR=", learning_rate))
# #         
# #         cat("\n\n========== Best Configuration by ESS/sec (sampling) ==========\n")
# #         print(best_by_ESS_sec %>% select(N, config, ESS_per_sec_samp_mean, min_ESS_mean, pct_divs_mean))
# #         
# #         ## Find best config by ESS/grad
# #         best_by_ESS_grad <-  agg_df %>%
# #           group_by(N) %>%
# #           slice_max(ESS_per_grad_samp_mean, n = 1) %>%
# #           mutate(config = paste0("n_burn=", n_burnin, ", LR=", learning_rate))
# #         
# #         cat("\n\n========== Best Configuration by ESS/grad (sampling) ==========\n")
# #         print(best_by_ESS_grad %>% select(N, config, ESS_per_grad_samp_mean, min_ESS_mean, pct_divs_mean))
# #         
# #         invisible(comparison_df)
# #   
# # }






##
## ---- Function for BayesMVP latex tables (Pilot Study 7) - with summary option: -----------------------------------------------------------------
##
generate_bayesmvp_latex_table <-  function(df,
                                          target_ESS,
                                          N_val = NULL,
                                          n_burnin_val = NULL,
                                          summary_table = FALSE,
                                          ##
                                          dp_time_burnin = 2,
                                          dp_time_samp = 2,
                                          dp_ESS_per_sec = 1,
                                          dp_ESS_per_grad = 2,
                                          dp_ESS_per_grad_sd = 2,
                                          dp_est_time_to_target = 2,
                                          dp_Rhat = 3,
                                          dp_nRhat = 3,
                                          dp_L = 1
) {
  
        require(dplyr)
        
        ## Filter by N if provided
        if (!is.null(N_val)) {
          df <-  df[df$N == N_val, ]
        } else {
          N_val <-  df$N[1]
        }
        
        ## Filter by n_burnin (only if not summary table)
        if (!summary_table && !is.null(n_burnin_val)) {
          df <-  df[df$n_burnin == n_burnin_val, ]
        }
        
        if (nrow(df) == 0) {
          warning(paste0("No data for N=", N_val))
          return(invisible(NULL))
        }
        
        df <-  df %>% mutate( target_ESS = target_ESS,
                             time_total_mins = round(time_total/60, 3),
                             est_time_to_target_ESS_mins = (time_burnin +
                                                              (target_ESS/min_ESS)*time_sampling + 
                                                              (target_ESS/min_ESS)*time_summaries)/60)
        
        ## Helper functions
        fmt_comma <-  function(x) {
          formatC(round(x), format = "d", big.mark = ",")
        }
        
        fmt_dec <-  function(x, digits = 2) {
          sub("^ +", "", formatC(x, format = "f", digits = digits))
        }
        
        ## Start building LaTeX
        latex_out <-  c()
        
        ## Table header
        latex_out <-  c(latex_out, "%%%%")
        latex_out <-  c(latex_out, "%%%%")
        latex_out <-  c(latex_out, "\\begin{table}[H]")
        latex_out <-  c(latex_out, "\\centering")
        latex_out <-  c(latex_out, "\\footnotesize")
        latex_out <-  c(latex_out, "\\begin{threeparttable}")
        
        ## Caption and label
        latex_out <-  c(latex_out, "\\caption{\\begin{tabular}[t]{@{}l@{}}")
        if (summary_table) {
          latex_out <-  c(latex_out, paste0("Summary of pilot study runs using \\textbf{BayesMVP} ($N = ", fmt_comma(N_val), "$). \\\\"))
          latex_out <-  c(latex_out, "On local HPC (using 64 chains for sampling, 8 for burnin).")
          latex_out <-  c(latex_out, "\\end{tabular}}")
          latex_out <-  c(latex_out, paste0("\\label{table_appendix_BayesMVP_pilot_study_summary_N_", N_val, "}"))
        } else {
          latex_out <-  c(latex_out, paste0("Results of pilot study runs using \\textbf{BayesMVP} ($N = ", fmt_comma(N_val), "$, $N_{burn} = ", n_burnin_val, "$). \\\\"))
          latex_out <-  c(latex_out, "On local HPC (using 64 chains for sampling, 8 for burnin).")
          latex_out <-  c(latex_out, "\\end{tabular}}")
          latex_out <-  c(latex_out, paste0("\\label{table_appendix_BayesMVP_pilot_study_N_", N_val, "_nburn_", n_burnin_val, "}"))
        }
        
        ## Begin tabular
        latex_out <-  c(latex_out, "\\begin{tabular}{llllllllll}")
        latex_out <-  c(latex_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        latex_out <-  c(latex_out, "\\toprule")
        latex_out <-  c(latex_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        ## Column headers
        if (summary_table) {
          latex_out <-  c(latex_out, "LR &")
        } else {
          latex_out <-  c(latex_out, "Seed &")
        }
        latex_out <-  c(latex_out, "\\begin{tabular}[c]{@{}l@{}} $ L_{burn}/ $ \\\\ $ L_{samp} $ \\end{tabular} &")
        latex_out <-  c(latex_out, "$ \\hat{R} / n\\hat{R} $ &")
        latex_out <-  c(latex_out, "$ ESS $ &")
        latex_out <-  c(latex_out, "\\begin{tabular}[c]{@{}l@{}} Burnin \\\\ time \\\\ (secs) \\end{tabular} &")
        latex_out <-  c(latex_out, "\\begin{tabular}[c]{@{}l@{}} Sampling \\\\ time \\\\ (secs) \\end{tabular} &")
        latex_out <-  c(latex_out, "\\begin{tabular}[c]{@{}l@{}} $ESS$ \\\\ /sec \\\\  (sampling) \\end{tabular} &")
        latex_out <-  c(latex_out, "\\begin{tabular}[c]{@{}l@{}} $ESS$ \\\\ /grad \\\\  (sampling) \\end{tabular} &")
        latex_out <-  c(latex_out, "\\begin{tabular}[c]{@{}l@{}} $\\approx$ time\\tnote{*} \\\\to \\\\ $\\approx$ target\\\\ $ESS$ \\\\ (mins) \\end{tabular} &")
        latex_out <-  c(latex_out, "\\begin{tabular}[c]{@{}l@{}} $\\approx$ $N_{iter}$ \\\\ needed \\\\ for \\\\ $\\approx$ target\\tnote{*} \\\\ $ESS$ \\end{tabular}")
        latex_out <-  c(latex_out, "\\\\ \\midrule")
        latex_out <-  c(latex_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        
        ##
        ## ---- SUMMARY TABLE: loop through n_burnin, then learning_rate ----
        ##
        if (summary_table) {
          
          n_burnin_values <-  sort(unique(df$n_burnin), decreasing = TRUE)
          
          for (n_burn in n_burnin_values) {
            
            df_burnin <-  df[df$n_burnin == n_burn, ]
            learning_rate_values <-  sort(unique(df_burnin$learning_rate))
            
            ## Group header for n_burnin
            latex_out <-  c(latex_out, paste0("\\multicolumn{10}{l}{$ \\mathbf{ N_{burn} = ", n_burn, " } $}"))
            latex_out <-  c(latex_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            latex_out <-  c(latex_out, "\\\\ \\hline")
            latex_out <-  c(latex_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            
            ## Loop through learning rates - show only averages
            for (lr in learning_rate_values) {
              
              df_group <-  df_burnin[df_burnin$learning_rate == lr, ]
              
              ## Calculate averages
              mean_L_burnin <-  mean(df_group$L_main_burnin, na.rm = TRUE)
              sd_L_burnin <-  sd(df_group$L_main_burnin, na.rm = TRUE)
              mean_L_samp <-  mean(df_group$L_main_samp, na.rm = TRUE)
              sd_L_samp <-  sd(df_group$L_main_samp, na.rm = TRUE)
              mean_Rhat <-  mean(df_group$max_Rhat, na.rm = TRUE)
              mean_nRhat <-  mean(df_group$max_nRhat, na.rm = TRUE)
              mean_ESS <-  mean(df_group$min_ESS, na.rm = TRUE)
              sd_ESS <-  sd(df_group$min_ESS, na.rm = TRUE)
              mean_time_burnin <-  mean(df_group$time_burnin, na.rm = TRUE)
              sd_time_burnin <-  sd(df_group$time_burnin, na.rm = TRUE)
              mean_time_samp <-  mean(df_group$time_sampling, na.rm = TRUE)
              sd_time_samp <-  sd(df_group$time_sampling, na.rm = TRUE)
              mean_ESS_per_sec <-  mean(df_group$ESS_per_sec_samp, na.rm = TRUE)
              sd_ESS_per_sec <-  sd(df_group$ESS_per_sec_samp, na.rm = TRUE)
              mean_ESS_per_grad <-  mean(df_group$ESS_per_grad_samp, na.rm = TRUE)
              sd_ESS_per_grad <-  sd(df_group$ESS_per_grad_samp, na.rm = TRUE)
              
              time_to_target_vec <-  df_group$est_time_to_target_ESS_mins
              mean_time_to_target <-  mean(time_to_target_vec, na.rm = TRUE)
              sd_time_to_target <-  sd(time_to_target_vec, na.rm = TRUE)
              
              n_iter_needed_vec <-  df_group$n_iter * (target_ESS / df_group$min_ESS)
              mean_n_iter_needed <-  mean(n_iter_needed_vec, na.rm = TRUE)
              sd_n_iter_needed <-  sd(n_iter_needed_vec, na.rm = TRUE)
              
              ## Summary row
              row_str <-  paste0(
                "$ ", lr, " $  & ",
                "\\begin{tabular}[c]{@{}l@{}} $ ", fmt_dec(mean_L_burnin, dp_L), " (",
                fmt_dec(sd_L_burnin, dp_L), ") $ \\\\ $ ", 
                fmt_dec(mean_L_samp, dp_L), " (", fmt_dec(sd_L_samp, dp_L), ") $ \\end{tabular} & ",
                "$ ", fmt_dec(mean_Rhat, dp_Rhat), " / ", fmt_dec(mean_nRhat, dp_nRhat), " $ & ",
                "\\begin{tabular}[c]{@{}l@{}} $ ", fmt_comma(round(mean_ESS)), " $ \\\\ $ (", 
                fmt_comma(round(sd_ESS)), ") $ \\end{tabular} & ",
                "\\begin{tabular}[c]{@{}l@{}} $ ", fmt_dec(mean_time_burnin, dp_time_burnin), " $ \\\\ $ (",
                fmt_dec(sd_time_burnin, dp_time_burnin), ") $ \\end{tabular} & ",
                "\\begin{tabular}[c]{@{}l@{}} $ ", fmt_dec(mean_time_samp, dp_time_samp), " $ \\\\ $ (",
                fmt_dec(sd_time_samp, dp_time_samp), ") $ \\end{tabular} & ",
                "\\begin{tabular}[c]{@{}l@{}} $ ", fmt_dec(mean_ESS_per_sec, dp_ESS_per_sec), " $ \\\\ $ (",
                fmt_dec(sd_ESS_per_sec, dp_ESS_per_sec), ") $ \\end{tabular} & ",
                "\\begin{tabular}[c]{@{}l@{}} $ ", fmt_dec(mean_ESS_per_grad, dp_ESS_per_grad), " $ \\\\ $ (", 
                fmt_dec(sd_ESS_per_grad, dp_ESS_per_grad_sd), ") $ \\end{tabular} & ",
                "\\begin{tabular}[c]{@{}l@{}} $ ", fmt_dec(mean_time_to_target, dp_est_time_to_target), " $ \\\\ $ (",
                fmt_dec(sd_time_to_target, dp_est_time_to_target), ") $ \\end{tabular} & ",
                "\\begin{tabular}[c]{@{}l@{}} $ ", fmt_comma(round(mean_n_iter_needed)), " $ \\\\ $ (", fmt_comma(round(sd_n_iter_needed)), ") $ \\end{tabular}"
              )
              
              latex_out <-  c(latex_out, row_str)
              latex_out <-  c(latex_out, "\\\\")
            }
            
            latex_out <-  c(latex_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            latex_out <-  c(latex_out, "\\bottomrule")
            latex_out <-  c(latex_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%")
          }
          
        } else {
          ##
          ## ---- FULL TABLE: original behavior - individual runs + averages ----
          ##
          learning_rate_values <-  sort(unique(df$learning_rate))
          
          for (lr in learning_rate_values) {
            
            df_group <-  df[df$learning_rate == lr, ]
            
            ## Group header
            latex_out <-  c(latex_out, paste0("\\multicolumn{10}{l}{$ \\mathbf{ LR = ", lr, " } $}"))
            latex_out <-  c(latex_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            latex_out <-  c(latex_out, "\\\\ \\hline")
            latex_out <-  c(latex_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            
            ## Data rows
            for (i in 1:nrow(df_group)) {
              
              row <-  df_group[i, ]
              
              n_iter_needed <-  round(row$n_iter * (target_ESS / row$min_ESS))
              
              row_str <-  paste0(
                "$ ", row$run, " $  &  ",
                "$ ", fmt_dec(row$L_main_burnin, dp_L), "/", fmt_dec(row$L_main_samp, dp_L), " $ & ",
                "$ ", fmt_dec(row$max_Rhat, dp_Rhat), " / ", fmt_dec(row$max_nRhat, dp_nRhat), " $  &  ",
                "$ ", fmt_comma(row$min_ESS), " $ & ",
                "$ ", fmt_dec(row$time_burnin, dp_time_burnin), " $ & ",
                "$ ", fmt_dec(row$time_sampling, dp_time_samp), " $ & ",
                "$ ", fmt_dec(row$ESS_per_sec_samp, dp_ESS_per_sec), " $ & ",
                "$ ", fmt_dec(row$ESS_per_grad_samp, dp_ESS_per_grad), " $ & ",
                "$ ", fmt_dec(row$est_time_to_target_ESS_mins, dp_est_time_to_target), " $ & ",
                "$ ", fmt_comma(n_iter_needed), " $"
              )
              
              latex_out <-  c(latex_out, row_str)
              latex_out <-  c(latex_out, "\\\\")
            }
            
            ## Calculate averages
            mean_L_burnin <-  mean(df_group$L_main_burnin, na.rm = TRUE)
            sd_L_burnin <-  sd(df_group$L_main_burnin, na.rm = TRUE)
            mean_L_samp <-  mean(df_group$L_main_samp, na.rm = TRUE)
            sd_L_samp <-  sd(df_group$L_main_samp, na.rm = TRUE)
            mean_Rhat <-  mean(df_group$max_Rhat, na.rm = TRUE)
            mean_nRhat <-  mean(df_group$max_nRhat, na.rm = TRUE)
            mean_ESS <-  mean(df_group$min_ESS, na.rm = TRUE)
            sd_ESS <-  sd(df_group$min_ESS, na.rm = TRUE)
            mean_time_burnin <-  mean(df_group$time_burnin, na.rm = TRUE)
            sd_time_burnin <-  sd(df_group$time_burnin, na.rm = TRUE)
            mean_time_samp <-  mean(df_group$time_sampling, na.rm = TRUE)
            sd_time_samp <-  sd(df_group$time_sampling, na.rm = TRUE)
            mean_ESS_per_sec <-  mean(df_group$ESS_per_sec_samp, na.rm = TRUE)
            sd_ESS_per_sec <-  sd(df_group$ESS_per_sec_samp, na.rm = TRUE)
            mean_ESS_per_grad <-  mean(df_group$ESS_per_grad_samp, na.rm = TRUE)
            sd_ESS_per_grad <-  sd(df_group$ESS_per_grad_samp, na.rm = TRUE)
            
            time_to_target_vec <-  df_group$est_time_to_target_ESS_mins
            mean_time_to_target <-  mean(time_to_target_vec, na.rm = TRUE)
            sd_time_to_target <-  sd(time_to_target_vec, na.rm = TRUE)
            
            n_iter_needed_vec <-  df_group$n_iter * (target_ESS / df_group$min_ESS)
            mean_n_iter_needed <-  mean(n_iter_needed_vec, na.rm = TRUE)
            sd_n_iter_needed <-  sd(n_iter_needed_vec, na.rm = TRUE)
            
            latex_out <-  c(latex_out, "\\midrule")
            latex_out <-  c(latex_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            
            avg_str <-  paste0(
              "\\textbf{Avg}.  & ",
              "\\begin{tabular}[c]{@{}l@{}} $ \\textbf{", fmt_dec(mean_L_burnin, dp_L), 
              " (", fmt_dec(sd_L_burnin, dp_L), ") } $ \\\\ $ \\textbf{ ", 
              fmt_dec(mean_L_samp, dp_L), " (", fmt_dec(sd_L_samp, dp_L), ") } $ \\end{tabular} & ",
              "$ \\textbf{", fmt_dec(mean_Rhat, dp_Rhat), "} / \\textbf{", fmt_dec(mean_nRhat, dp_nRhat), "} $ & ",
              "\\begin{tabular}[c]{@{}l@{}} $ \\textbf{", fmt_comma(round(mean_ESS)), "} $ \\\\ $ \\textbf{(",
              fmt_comma(round(sd_ESS)), ")} $ \\end{tabular} & ",
              "\\begin{tabular}[c]{@{}l@{}} $ \\textbf{", fmt_dec(mean_time_burnin, dp_time_burnin), "} $ \\\\ $ \\textbf{(", 
              fmt_dec(sd_time_burnin, dp_time_burnin), ")} $ \\end{tabular} & ",
              "\\begin{tabular}[c]{@{}l@{}} $ \\textbf{", fmt_dec(mean_time_samp, dp_time_samp), "} $ \\\\ $ \\textbf{(",
              fmt_dec(sd_time_samp, dp_time_samp), ")} $ \\end{tabular} & ",
              "\\begin{tabular}[c]{@{}l@{}} $ \\textbf{", fmt_dec(mean_ESS_per_sec, dp_ESS_per_sec), "} $ \\\\ $ \\textbf{(",
              fmt_dec(sd_ESS_per_sec, dp_ESS_per_sec), ")} $ \\end{tabular} & ",
              "\\begin{tabular}[c]{@{}l@{}} $ \\textbf{", fmt_dec(mean_ESS_per_grad, dp_ESS_per_grad), "} $ \\\\ $ \\textbf{(",
              fmt_dec(sd_ESS_per_grad, dp_ESS_per_grad_sd), ")} $ \\end{tabular} & ",
              "\\begin{tabular}[c]{@{}l@{}} $ \\textbf{", fmt_dec(mean_time_to_target, dp_est_time_to_target), "} $ \\\\ $ \\textbf{(",
              fmt_dec(sd_time_to_target, dp_est_time_to_target), ")} $ \\end{tabular} & ",
              "\\begin{tabular}[c]{@{}l@{}} $ \\textbf{", fmt_comma(round(mean_n_iter_needed)), "} $ \\\\ $ \\textbf{(", 
              fmt_comma(round(sd_n_iter_needed)), ")} $ \\end{tabular}"
            )
            
            latex_out <-  c(latex_out, avg_str)
            latex_out <-  c(latex_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%")
            latex_out <-  c(latex_out, "\\\\ \\bottomrule")
            latex_out <-  c(latex_out, "%%%%%%%%%%%%%%%%%%%%%%%%%%%")
          }
        }
        
        ## Close tabular
        latex_out <-  c(latex_out, "\\end{tabular}")
        
        ## Table notes
        latex_out <-  c(latex_out, "\\begin{tablenotes}")
        latex_out <-  c(latex_out, "\\footnotesize")
        latex_out <-  c(latex_out, paste0("\\item[*] Target ESS for $N = ", fmt_comma(N_val),
                                         "$ is (from pilot study \\#3) equal to $", fmt_comma(target_ESS), "$."))
        latex_out <-  c(latex_out, "\\end{tablenotes}")
        
        ## Close threeparttable and table
        latex_out <-  c(latex_out, "\\end{threeparttable}")
        latex_out <-  c(latex_out, "\\end{table}")
        latex_out <-  c(latex_out, "%%%%")
        latex_out <-  c(latex_out, "%%%%")
        
        ## Print with cat
        cat(paste(latex_out, collapse = "\n"))
        
        ## Return invisibly
        invisible(paste(latex_out, collapse = "\n"))
  
}

## Usage:
## Full tables (one per n_burnin):
# generate_bayesmvp_latex_table(df = results_BayesMVP_df,
#                               target_ESS = 5000,
#                               N_val = 10000,
#                               n_burnin_val = 250,
#                               summary_table = FALSE)

## Summary table (all n_burnin in one table):
# generate_bayesmvp_latex_table(df = results_BayesMVP_df,
#                               target_ESS = 5000,
#                               N_val = 10000,
#                               summary_table = TRUE)
##
## ---- Wrapper function to generate all 3 tables for a given N: ---------------------------------------------------------------------------
##
##
## ---- Wrapper function to generate tables for a given N: ---------------------------------------------------------------------------------
##
generate_all_bayesmvp_tables <-  function(df,
                                         target_ESS,
                                         N_val,
                                         n_burnin_vec = c(1000, 500, 250),
                                         summary_table = FALSE,
                                         ...
) {
        
        all_tables <-  list()
        
        if (summary_table) {
          ## Summary table: just ONE table with all n_burnin values
          cat("\n\n")
          cat(paste0("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
          cat(paste0("%% N = ", N_val, " (Summary Table)\n"))
          cat(paste0("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
          cat("\n")
          
          table_str <-  generate_bayesmvp_latex_table(
            df = df,
            target_ESS = target_ESS,
            N_val = N_val,
            summary_table = TRUE,
            ...
          )
          
          all_tables[["summary"]] <-  table_str
          
        } else {
          ## Full tables: one per n_burnin
          for (n_burnin in n_burnin_vec) {
            
            cat("\n\n")
            cat(paste0("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
            cat(paste0("%% N = ", N_val, ", n_burnin = ", n_burnin, "\n"))
            cat(paste0("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n"))
            cat("\n")
            
            table_str <-  generate_bayesmvp_latex_table(
              df = df,
              target_ESS = target_ESS,
              N_val = N_val,
              n_burnin_val = n_burnin,
              summary_table = FALSE,
              ...
            )
            
            all_tables[[as.character(n_burnin)]] <-  table_str
          }
        }
        
        invisible(all_tables)
  
}




##
## -| --------- Comparing trajectory-length schemes: fn_ps7_compare_tau_schemes() -------------------------------------
##
## The original ChESSR adaptation scores a trajectory by its change in KINETIC ENERGY and ascends
## that in log(tau) with ADAM. Kinetic energy is a function of the state alone and is invariant to
## the momentum flip, so for the leapfrog-plus-flip involution Phi and acceptance probability alpha,
## detailed balance gives E[alpha * h(Phi z)] = E[alpha * h(z)] for ANY h, hence
##
##     E[ alpha * (KE_prop - KE_0) ] = 0        exactly, for every tau and every eps
##
## once the chains are stationary. The objective therefore has no optimum: during the transient it
## pushes tau up (chains outside the typical set convert potential energy into kinetic energy), and
## after convergence it is zero-mean noise that ADAM normalises to a step of roughly +/- the learning
## rate, so log(tau) random-walks. ChEES escapes this by scoring the SQUARED change in distance from
## the running mean rather than a linear difference of a state function.
##
## This function puts the schemes side by side. It does three things, in order of how decisive they
## are:
##
##   1. OVERALL   - each scheme's mean on every metric, with the run-to-run spread. The spread on
##                  L_main is the interesting column: a random walk in log(tau) shows up as a large
##                  sd of the adapted L across replicate runs, and that is the direct evidence for
##                  the argument above.
##   2. PAIRED    - baseline scheme vs each alternative, matched on every OTHER swept axis, so
##                  nothing is confounded. Reports the win count and the mean relative change.
##   3. FIXED-L   - if fixed_L_* rows are present, the L -> metric curve, plus where the ADAPTED
##                  schemes landed. This is the only view that says whether the adapted L is near
##                  the ESS-per-gradient optimum, which the adaptation cannot determine independently.
##
## Rows that failed to converge are excluded on the same terms as fn_ps7_report and listed.
##
#' @param summary_df    per-configuration summary (summarize_ps7_by_config()) containing tau_scheme.
#' @param baseline      the scheme every other one is compared against; "KE" is the original ChESSR.
#' @param metrics       metric columns to report; NULL = auto-detect the usual ps7 ones.
#' @param digits        printing precision.
#' @param max_Rhat_ok   convergence limits, as in fn_ps7_report().
#' @param max_nRhat_ok  .
#' @param max_pct_divs_ok .
#' @return invisibly, list(overall, paired, fixed_L, dropped).
#' @export
fn_ps7_compare_tau_schemes <-  function( summary_df,
                                        baseline = "KE",
                                        metrics = NULL,
                                        digits = 4,
                                        max_Rhat_ok = 1.05,
                                        max_nRhat_ok = 1.05,
                                        max_pct_divs_ok = 1.0
) {

    stopifnot(is.data.frame(x = summary_df))
    summary_df <-  tibble::as_tibble(x = summary_df)
    ##
    if (!"tau_scheme" %in% names(x = summary_df)) {
        stop("fn_ps7_compare_tau_schemes: no 'tau_scheme' column. Re-run summarize_ps7_results() after",
             " re-sourcing ps_7_MCMC_settings_BayesMVP_functions.R - an older copy did not emit it.")
    }
    if (nrow(x = summary_df) == 0) {
        cat("fn_ps7_compare_tau_schemes: nothing to compare (0 rows).\n")
        return(invisible(NULL))
    }
    ##
    ## ---- metrics to report, higher-is-better flagged
    metric_catalogue <-  list( list( column = "min_ESS_per_grad_sampling_mean",               label = "ESS_per_grad",  direction = +1),
                              list( column = "ESS_per_grad_samp_mean",                       label = "ESS_per_grad",  direction = +1),
                              list( column = "min_ESS_per_sec_sampling_mean",                label = "ESS_per_sec",   direction = +1),
                              list( column = "ESS_per_sec_samp_mean",                        label = "ESS_per_sec",   direction = +1),
                              list( column = "est_time_to_target_ESS_wo_summaries_mins_mean", label = "mins_to_ESS",   direction = -1),
                              list( column = "est_time_to_target_ESS_mins_mean",              label = "mins_to_ESS",   direction = -1),
                              list( column = "time_burnin_mean",                              label = "time_burnin",   direction = -1),
                              list( column = "pct_divs_mean",                                 label = "pct_divs",      direction = -1),
                              list( column = "L_main_samp_mean",                              label = "L_main",        direction = 0))
    ##
    available <-  Filter( f = function(m) m$column %in% names(x = summary_df), x = metric_catalogue)
    seen <-  character()
    available <-  Filter( f = function(m) { keep <-  !(m$label %in% seen); seen <<- c(seen, m$label); keep }, x = available)
    if (!is.null(metrics)) available <-  Filter( f = function(m) m$column %in% metrics, x = available)
    if (length(x = available) == 0) stop("fn_ps7_compare_tau_schemes: none of the expected metric columns are present.")
    ##
    metric_columns <-  vapply( X = available, FUN = function(m) m$column, FUN.VALUE = "")
    metric_labels  <-  vapply( X = available, FUN = function(m) m$label,  FUN.VALUE = "")
    directions     <-  vapply( X = available, FUN = function(m) m$direction, FUN.VALUE = 0)
    ##
    ## ---- drop non-converged configurations, exactly as fn_ps7_report does
    fn_fails <-  function( column, limit) {
        if (!column %in% names(x = summary_df)) return(rep( x = FALSE, times = nrow(x = summary_df)))
        value <-  summary_df[[column]]
        is.finite(value) & value > limit
    }
    finite_row <-  Reduce( f = `&`, x = lapply( X = summary_df[metric_columns], FUN = is.finite))
    fails_rhat <-  fn_fails( column = "max_Rhat_mean", limit = max_Rhat_ok) |
                  fn_fails( column = "max_nRhat_mean", limit = max_nRhat_ok)
    fails_divs <-  fn_fails( column = "pct_divs_mean", limit = max_pct_divs_ok)
    ##
    keep <-  finite_row & !fails_rhat & !fails_divs
    dropped <-  summary_df[!keep, , drop = FALSE]
    kept <-  summary_df[keep, , drop = FALSE]
    ##
    if (nrow(x = dropped) > 0) {
        cat("\n---- EXCLUDED (did not converge): ", nrow(x = dropped), " configuration(s) ----\n", sep = "")
        excluded_reason <-  ifelse( test = !finite_row[!keep], yes = "non-finite metric",
                            no = ifelse( test = fails_rhat[!keep], yes = "Rhat too high", no = "too many divergences"))
        print(x = tibble::tibble(tau_scheme = dropped$tau_scheme, reason = excluded_reason), n = Inf, width = Inf)
    }
    if (nrow(x = kept) == 0) { cat("\nNothing left after the convergence filter.\n"); return(invisible(NULL)) }
    ##
    schemes <-  sort(x = unique(x = kept$tau_scheme))
    cat("\n================ TRAJECTORY-LENGTH SCHEMES PRESENT ================\n")
    cat(paste0("  ", schemes, "  (", as.integer(table(kept$tau_scheme)[schemes]), " configuration(s))", collapse = "\n"), "\n", sep = "")
    ##
    ## ================================================================================ 1. OVERALL
    overall_rows <-  list()
    for (scheme in schemes) {
        rows_s <-  kept[kept$tau_scheme == scheme, , drop = FALSE]
        entry <-  tibble::tibble( tau_scheme = scheme, n_configs = nrow(x = rows_s))
        for (mm in seq_along(along.with = metric_columns)) {
            values <-  rows_s[[metric_columns[mm]]]
            entry[[metric_labels[mm]]] <-  round( x = mean( x = values, na.rm = TRUE), digits = digits)
            entry[[paste0(metric_labels[mm], "_sd")]] <-  round( x = stats::sd( x = values, na.rm = TRUE), digits = digits)
        }
        overall_rows[[length(overall_rows) + 1]] <-  entry
    }
    overall <-  dplyr::bind_rows(overall_rows)
    cat("\n================ 1. OVERALL (mean over configurations, sd across them) ================\n")
    print(overall, n = Inf, width = Inf)
    if ("L_main_sd" %in% names(x = overall)) {
        cat("\n  NOTE: L_main_sd is the spread of trajectory lengths across configurations.\n",
            "        It mixes setting differences and run variability, so it cannot establish a random walk.\n",
            "        Inspect the saved within-run tau/epsilon histories at matched settings.\n", sep = "")
    }
    ##
    ## ================================================================================ 2. PAIRED
    ## Match on every column that is NOT a metric, NOT bookkeeping and NOT part of the scheme label.
    metric_like <-  grep( pattern = "_(mean|sd|min|max|median|g_mean)$", x = names(x = kept), value = TRUE)
    bookkeeping <-  intersect( x = c("n", "n_runs", "n_missing", "run"), y = names(x = kept))
    scheme_columns <-  intersect( x = c("tau_scheme", "burnin_algorithm", "tau_weight_by_p_jump", "manual_L", "manual_tau_value"),
                                 y = names(x = kept))
    match_columns <-  setdiff( x = names(x = kept), y = c(metric_like, bookkeeping, scheme_columns))
    match_columns <-  match_columns[vapply( X = kept[match_columns], FUN = function(column) length(x = unique(x = column)) > 1, FUN.VALUE = TRUE)]
    ##
    paired <-  NULL
    if (!baseline %in% schemes) {
        cat("\n================ 2. PAIRED ================\n")
        cat("  baseline scheme '", baseline, "' is not present, so no paired comparison.\n", sep = "")
    } else {
        cat("\n================ 2. PAIRED vs '", baseline, "' (matched on: ",
            if (length(x = match_columns) == 0) "nothing else varies" else paste(match_columns, collapse = ", "), ") ================\n", sep = "")
        ##
        fn_key <-  function(rows) {
            if (length(x = match_columns) == 0) return(rep( x = "all", times = nrow(x = rows)))
            do.call( what = paste, args = c(lapply( X = match_columns, FUN = function(column) as.character(rows[[column]])), list(sep = "|")))
        }
        kept$.match_key <-  fn_key(rows = kept)
        base_rows <-  kept[kept$tau_scheme == baseline, , drop = FALSE]
        ##
        paired_rows <-  list()
        for (scheme in setdiff(x = schemes, y = baseline)) {
            alt_rows <-  kept[kept$tau_scheme == scheme, , drop = FALSE]
            shared <-  intersect( x = base_rows$.match_key, y = alt_rows$.match_key)
            if (length(x = shared) == 0) {
                cat("  ", scheme, ": no matched pairs with the baseline.\n", sep = "")
                next
            }
            entry <-  tibble::tibble( tau_scheme = scheme, n_pairs = length(x = shared))
            for (mm in seq_along(along.with = metric_columns)) {
                base_values <-  base_rows[[metric_columns[mm]]][match( x = shared, table = base_rows$.match_key)]
                alt_values  <-  alt_rows[[metric_columns[mm]]][match( x = shared, table = alt_rows$.match_key)]
                ##
                ## relative change is only meaningful where the baseline is nonzero
                usable <-  is.finite(base_values) & is.finite(alt_values) & base_values != 0
                rel <-  if (any(usable)) 100 * mean(x = (alt_values[usable] - base_values[usable]) / abs(base_values[usable])) else NA_real_
                entry[[paste0(metric_labels[mm], "_pct")]] <-  round( x = rel, digits = 1)
                ##
                if (directions[mm] != 0) {
                    wins <-  sum( (alt_values - base_values) * directions[mm] > 0, na.rm = TRUE)
                    entry[[paste0(metric_labels[mm], "_wins")]] <-  paste0(wins, "/", sum(is.finite(base_values) & is.finite(alt_values)))
                }
            }
            paired_rows[[length(paired_rows) + 1]] <-  entry
        }
        if (length(x = paired_rows) > 0) {
            paired <-  dplyr::bind_rows(paired_rows)
            print(paired, n = Inf, width = Inf)
            cat("\n  *_pct  = mean per cent change vs the baseline (SIGNED: positive means the metric went up,\n",
                "           which is better for ESS_per_grad / ESS_per_sec and worse for mins_to_ESS / pct_divs).\n",
                "  *_wins = matched pairs in which the alternative beat the baseline.\n", sep = "")
        }
        kept$.match_key <-  NULL
    }
    ##
    ## ================================================================================ 3. FIXED-L
    fixed_L <-  NULL
    if ("manual_L" %in% names(x = kept) && any(!is.na(kept$manual_L))) {
        cat("\n================ 3. FIXED-L SWEEP ================\n")
        pinned <-  kept[!is.na(kept$manual_L), , drop = FALSE]
        rows_L <-  list()
        for (L_value in sort(x = unique(x = pinned$manual_L))) {
            rows_s <-  pinned[pinned$manual_L == L_value, , drop = FALSE]
            entry <-  tibble::tibble( L = L_value, n_configs = nrow(x = rows_s))
            for (mm in seq_along(along.with = metric_columns)) {
                entry[[metric_labels[mm]]] <-  round( x = mean( x = rows_s[[metric_columns[mm]]], na.rm = TRUE), digits = digits)
            }
            rows_L[[length(rows_L) + 1]] <-  entry
        }
        fixed_L <-  dplyr::bind_rows(rows_L)
        print(fixed_L, n = Inf, width = Inf)
        ##
        ## where the ADAPTED schemes landed, so the two can be read against each other
        adaptive_rows <-  is.na(kept$manual_L)
        if ("manual_tau_value" %in% names(x = kept)) adaptive_rows <-  adaptive_rows & is.na(kept$manual_tau_value)
        adapted <-  kept[adaptive_rows, , drop = FALSE]
        if (nrow(x = adapted) > 0 && "L_main_samp_mean" %in% names(x = adapted)) {
            cat("\n  Adapted schemes, for comparison (L reached vs the swept grid above):\n")
            for (scheme in sort(x = unique(x = adapted$tau_scheme))) {
                rows_s <-  adapted[adapted$tau_scheme == scheme, , drop = FALSE]
                cat(paste0("    ", formatC(as.character(scheme), format = "s", width = 16, flag = "-"), " L = ", sub("^ +", "", formatC(mean(x = rows_s$L_main_samp_mean, na.rm = TRUE), format = "f", digits = 2)), " (sd ", sub("^ +", "", formatC(stats::sd(x = rows_s$L_main_samp_mean, na.rm = TRUE), format = "f", digits = 2)), " over ", formatC(as.integer(nrow(x = rows_s)), format = "d"), " configs)\n"))
            }
        }
        ##
        if ("ESS_per_grad" %in% names(x = fixed_L)) {
            best <-  fixed_L$L[which.max(x = fixed_L$ESS_per_grad)]
            cat("\n  Best ESS per gradient on the swept grid: L = ", best, ".\n",
                "  If an adapted scheme sits far from this, its trajectory length is not being chosen well,\n",
                "  and the integrator comparison at the adapted L is confounded by that.\n", sep = "")
        }
    }
    ##
    invisible(list( overall = overall, paired = paired, fixed_L = fixed_L, dropped = dropped))

}































