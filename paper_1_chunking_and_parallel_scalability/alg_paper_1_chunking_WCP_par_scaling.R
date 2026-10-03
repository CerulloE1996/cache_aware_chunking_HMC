##
## ======================================================================================================================================
## alg_paper_1_chunking_WCP_par_scaling.R
##
## RUN THIS FILE in RStudio for the unified Paper 1 sampling study.
## Edit the settings below, save, then click Source. With run_benchmark = TRUE, Source starts the benchmark.
## This replaces the old PS1 sampling and PS2 main parallel-scaling entry points for new study runs.
## The separate ps_1_burnin_optimizing_N_chunks_and_WCP_threads.R study remains separate.
##
# rm(list = ls())
# .rs.restartR()
##
## ---- Run controls and paths ----------------------------------------------------------------------------------------------------------
##
{
      run_benchmark <-  TRUE                         ## FALSE = preview the grid only; TRUE = preview, then run it.
      ##
      make_figures_and_tables <-  TRUE               ## Export from saved results after the benchmark finishes.
      ##
      device <-  if (parallel::detectCores() > 16) "HPC" else "Laptop"
      algorithm_study_dir <-  path.expand("~/Documents/Work/PhD_work/Alg_paper_analysis")
      ##
      source( file = file.path(algorithm_study_dir, 
                               "paper_1_chunking_and_parallel_scalability", 
                               "R_fns_alg_paper_1_chunking_WCP_par_scaling.R"),
              local = TRUE)
      source( file = file.path(algorithm_study_dir, 
                               "paper_1_chunking_and_parallel_scalability",
                               "R_fns_alg_paper_1_figures_tables.R"),
              local = TRUE)
      ##
      source( file = file.path(algorithm_study_dir,
                               "paper_1_chunking_and_parallel_scalability",
                               "R_fns_alg_paper_1_Mplus.R"),
              local = TRUE)
      ##
      ## ---- Every study setting is assigned in this runner; the helper supplies no defaults:
      ##
      paper1_settings <-  list( algorithm_study_dir = algorithm_study_dir,
                                  device = device,
                                  package_stack = "NicoStan_BayesMVP", 
                                  package_library = NULL,  ## NULL uses the current R libraries; set a path for an isolated build.
                                  run_each_Stan_case_in_fresh_R_process = FALSE,
                                  stan_reduce_sum_type = "static",
                                  mplus = list(runner = fn_paper1_run_Mplus,
                                               verifier = fn_paper1_verify_Mplus_run))
      ## NULL compiles the selected CmdStanR models; NicoStan always requires NULL here:
      Stan_model_obj <-  NULL
}
##
## ---- Selected implementations, sample sizes and repeats ------------------------------------------------------------------------------
##
{
      paper1_settings$N_vec <-  c(500, 2500, 10000, 50000)
      ##
      ## EVERY algorithm uses the short/long timing method selected below.
      ## n_runs counts repetitions of the WHOLE PAIR: 2 means short, long, short, long per configuration.
      paper1_settings$n_runs       <-  4
      paper1_settings$mplus_n_runs <-  3  ## Retain the existing Mplus repetitions; add repeat 4 only for BayesMVP/Stan.
      paper1_settings$preserve_historical_stan_runs <-  c(3, 4)
      ##
      paper1_settings$data_seed <-  123
      ##
      paper1_settings$algorithms <-  c(  "AD_Stan", ## ---- done (HPC + Laptop)
                                         "AD_Stan_chunked", ## ---- done (HPC + Laptop)
                                         "AD_Stan_tape_chunked", ## ---- done (HPC + Laptop)
                                         "AD_Stan_WCP",  ## ---- done (HPC + Laptop),
                                         ##
                                         "Mplus_WCP",  ## ---- done (HPC + Laptop) ---- re-running w/ BITER as well
                                         "Mplus_standard", ## ---- done (HPC + Laptop) ---- re-running w/ BITER as well
                                         ##
                                         "MD_BayesMVP",  ## ---- done (HPC + Laptop)
                                         "MD_BayesMVP_WCP")  ## ---- done (HPC + Laptop)
                                   
      ##
      ##
      ## ---- Child process of the "fresh R process per algorithm" mode below: measure only the algorithm it was started for,
      ##      and leave the figures and tables to the parent session.
      ##
      if (nzchar(Sys.getenv("PAPER1_CHILD_ALGORITHM"))) {
            paper1_settings$algorithms <-  Sys.getenv("PAPER1_CHILD_ALGORITHM")
            make_figures_and_tables <-  FALSE
            message(paste0("\033[36mPaper 1 child process: measuring ", paper1_settings$algorithms, " only.\033[0m"))
      }
      paper1_settings$N_by_algorithm$MD_BayesMVP <-  paper1_settings$N_vec
      ##
      paper1_settings$N_by_algorithm$MD_BayesMVP_WCP <-  paper1_settings$N_vec
      ##
      paper1_settings$N_by_algorithm$AD_Stan <-  paper1_settings$N_vec
      ##
      paper1_settings$N_by_algorithm$AD_Stan_chunked <-  paper1_settings$N_vec
      ##
      paper1_settings$N_by_algorithm$AD_Stan_tape_chunked <-  paper1_settings$N_vec
      ##
      paper1_settings$N_by_algorithm$AD_Stan_WCP <-  paper1_settings$N_vec
      ##
      paper1_settings$N_by_algorithm$Mplus_standard <-  paper1_settings$N_vec
      paper1_settings$N_by_algorithm$Mplus_WCP <-  paper1_settings$N_vec
      ## Optional multi-process arm, if added to algorithms above.
      paper1_settings$N_by_algorithm$MD_BayesMVP_multi_process <-  paper1_settings$N_vec
}
##
## ---- Chunk counts: use the same chunk grid for BayesMVP with and without WCP ----------------------------------------------------------
##
{
      bayesmvp_chunk_candidates <-  list( "500"   = c(1, 2, 4, 10),
                                          ##
                                          "2500"  = c(1, 2, 4, 10, 25),
                                          ##
                                          "10000" = c(1, 4, 10, 25, 50, 100),
                                          ##
                                          ## "50000" = c(1, 4, 10, 25, 50, 100, 250, 500))
                                          "50000"    = c(1, 4, 10, 25, 50, 100, 250, 500, 1000))
      ##
      paper1_settings$chunks_by_algorithm$MD_BayesMVP     <- bayesmvp_chunk_candidates
      paper1_settings$chunks_by_algorithm$MD_BayesMVP_WCP <- bayesmvp_chunk_candidates
      paper1_settings$chunks_by_algorithm$MD_BayesMVP_multi_process <- bayesmvp_chunk_candidates
      ##
      ## ---- 2026-09-28: Stan chunk counts extended past the L2 line per active thread (the searches at N = 10,000 and 50,000
      ##      selected the largest tested count on both machines); the previous grids are kept above each line.
      ##
      stan_chunk_candidates <-  list( "500"   = c(2, 4, 10),
                                      ##
                                      ## "2500"  = c(4, 10, 25, 50),
                                      "2500"     = c(4, 10, 25, 50, 100, 250),
                                      ##
                                      ## "10000" = c(4, 10, 25, 50, 100),
                                      "10000"    = c(4, 10, 25, 50, 100, 250, 500, 1000),
                                      ##
                                      ## "50000" = c(10, 25, 50, 100, 250, 500))
                                      "50000"    = c(10, 25, 50, 100, 250, 500, 1000, 2000, 5000))
      ##
      paper1_settings$chunks_by_algorithm$AD_Stan_chunked      <- stan_chunk_candidates
      paper1_settings$chunks_by_algorithm$AD_Stan_tape_chunked <- stan_chunk_candidates
      paper1_settings$chunks_by_algorithm$AD_Stan_WCP          <- stan_chunk_candidates
}
##
## ---- PS2 chain counts with sparse joint chunk / WCP search ----------------------------------------------------------------------------
##
## Chunks and WCP are crossed independently for BayesMVP and for Stan reduce_sum.
## Retain PS2 chain counts; total threads = chains x WCP must fit within the device limit.
## The chunking-only arms use one thread per chain; WCP candidate lists contain only values greater than one.
##
if (device == "HPC") {

    paper1_settings$n_threads_vec <-  c(1, 2, 4, 8, 16, 32, 64, 96, 180)
    paper1_settings$total_thread_limit <- 180
    ##
    ## ---- Extra ONE-CHAIN chunk counts for the matched-serial references (MD_BayesMVP and AD_Stan_tape_chunked only):
    ##      the WCP-only configurations use N_chunks = N_threads_per_chain (6 and 8 with 16 chains at small N; 2 at 2500 for Stan),
    ##      which the one-chain chunk grids below do not contain, so their parallel efficiency had no serial reference.
    paper1_settings$serial_reference_extra_chunks <-  list( "500"  = c(6, 8),
                                                            "2500" = c(2, 6, 8))
    ##
    sampling_chain_counts <-  c(1, 2, 4, 8, 16, 32, 64, 96, 180)
    # WCP_chain_counts      <-  c(4, 8, 16)
    ## 2026-10-03: narrow WCP up to 90 x 2 (see below):
    WCP_chain_counts      <-  c(4, 8, 16, 24, 32, 45, 48, 64, 90)
    ##
    ## n_WCP candidates, keyed by N and then sampling chain count:
    ##
    WCP_candidates_for_4_chains_small_N   <- c(2, 4, 8)
    WCP_candidates_for_8_chains_small_N   <- c(2, 4, 8)
    WCP_candidates_for_16_chains_small_N  <- c(2, 4, 6, 8)
    ##
    WCP_candidates_for_4_chains_big_N  <- c(2, 4, 8, 16, 24, 32, 44)
    WCP_candidates_for_8_chains_big_N  <- c(2, 4, 8, 12, 16, 22)
    WCP_candidates_for_16_chains_big_N <- c(2, 4, 6, 8, 11)
    ##
    bayesmvp_WCP_candidates <-  list(
        "500"   = list( "4" = WCP_candidates_for_4_chains_small_N,
                        "8" = WCP_candidates_for_8_chains_small_N,
                       "16" = WCP_candidates_for_16_chains_small_N),
                       # "24" = WCP_candidates_for_24_chains),
        ##
        "2500"  = list( "4" = WCP_candidates_for_4_chains_small_N,
                        "8" = WCP_candidates_for_8_chains_small_N,
                       "16" = WCP_candidates_for_16_chains_small_N),
                       # "24" = WCP_candidates_for_24_chains),
        ##
        "10000" = list( "4" = WCP_candidates_for_4_chains_big_N,
                        "8" = WCP_candidates_for_8_chains_big_N,
                       "16" = WCP_candidates_for_16_chains_big_N),
                       # "24" = WCP_candidates_for_24_chains),
        ##
        "50000" = list( "4" = WCP_candidates_for_4_chains_big_N,
                        "8" = WCP_candidates_for_8_chains_big_N,
                       "16" = WCP_candidates_for_16_chains_big_N))
                       # "24" = WCP_candidates_for_24_chains))
    ##
    ## ---- 2026-10-03: narrow WCP with more chains (up to 90 chains x 2 threads), so that WCP + chunking can fill
    ##      64-180 threads with N_threads/chain = 2 or 4, as on the laptop (where 8 x 2 fills all 16 threads);
    ##      the grid above stops at 16 chains, which forced 6-11 threads per chain at 96-176 threads.
    ##      Added for every N; Stan and Mplus copy these candidates below.
    ##
    narrow_WCP_candidates <-  list( "24" = c(4),       ## 24 x 4 =  96 threads
                                    "32" = c(2, 4),    ## 32 x 2 =  64, 32 x 4 = 128
                                    "45" = c(4),       ## 45 x 4 = 180
                                    "48" = c(2),       ## 48 x 2 =  96
                                    "64" = c(2),       ## 64 x 2 = 128
                                    "90" = c(2))       ## 90 x 2 = 180
    ##
    for (N_key in names(bayesmvp_WCP_candidates)) {
          bayesmvp_WCP_candidates[[N_key]] <-  c(bayesmvp_WCP_candidates[[N_key]], narrow_WCP_candidates)
    }
    ##
    stan_WCP_candidates <- bayesmvp_WCP_candidates
    ##
    # 4*WCP_candidates_for_4_chains_big_N
    # 8*WCP_candidates_for_8_chains_big_N
    # 16*WCP_candidates_for_16_chains_big_N
    # # 24*WCP_candidates_for_24_chains
    # # # 45*c(2, 4)

} else {

    paper1_settings$n_threads_vec <- c(1, 2, 4, 8, 16)
    paper1_settings$total_thread_limit <- 16
    ##
    ## ---- Extra ONE-CHAIN chunk counts for the matched-serial references (see the HPC block above): the laptop WCP-only
    ##      configurations use 2, 4 or 8 chunks, and 8 (and 2 for Stan at N = 2500) are missing from the one-chain chunk grids.
    paper1_settings$serial_reference_extra_chunks <-  list( "500"  = c(8),
                                                            "2500" = c(2, 8))
    ##
    sampling_chain_counts <-  c(1, 2, 4, 8, 16)
    ##
    # WCP_chain_counts      <-  c(2, 4)
    WCP_chain_counts      <-  c(2, 4, 8)
    ##
    bayesmvp_WCP_candidates <-  list(
        "500"   = list( "2" = c(2, 4, 8),
                        "4" = c(2, 4),
                        "8" = c(2)),
        ##
        "2500"  = list( "2" = c(2, 4, 8),
                        "4" = c(2, 4),
                        "8" = c(2)),
        ##
        "10000" = list( "2" = c(2, 4, 8),
                        "4" = c(2, 4),
                        "8" = c(2)),
        ##
        "50000" = list( "2" = c(2, 4, 8),
                        "4" = c(2, 4),
                        "8" = c(2)))
    ##
    stan_WCP_candidates <-  list(
        "500"   = list( "2" = c(2, 4, 8),
                        "4" = c(2, 4),
                        "8" = c(2)),
        ##
        "2500"  = list( "2" = c(2, 4, 8),
                        "4" = c(2, 4),
                        "8" = c(2)),
        ##
        "10000" = list( "2" = c(2, 4, 8),
                        "4" = c(2, 4),
                        "8" = c(2)),
        ##
        "50000" = list( "2" = c(2, 4, 8),
                        "4" = c(2, 4),
                        "8" = c(2)))
    # ##
    # 2*c(2, 4, 8)
    # 4*c(2, 4)

}
##
## Start Mplus with the same per-N/per-chain WCP grid; this copy can be edited independently.
##
mplus_WCP_candidates <-  stan_WCP_candidates
##
## 2026-10-03: Mplus gets only the two-thread narrow allocations (32, 48, 64 and 90 chains x 2 threads;
## HPC only):
##
for (N_key in names(mplus_WCP_candidates)) {
      mplus_WCP_candidates[[N_key]][c("24", "45")] <-  NULL
      if (!is.null(mplus_WCP_candidates[[N_key]][["32"]])) mplus_WCP_candidates[[N_key]][["32"]] <-  c(2)
}
##
## Chain counts come from PS2 sampling, independently of the selected WCP candidates.
##
{
      paper1_settings$n_chains_by_algorithm$MD_BayesMVP               <- sampling_chain_counts
      paper1_settings$n_chains_by_algorithm$MD_BayesMVP_multi_process <- sampling_chain_counts
      paper1_settings$n_chains_by_algorithm$MD_BayesMVP_WCP           <- WCP_chain_counts
      ##
      paper1_settings$n_chains_by_algorithm$AD_Stan              <- sampling_chain_counts
      paper1_settings$n_chains_by_algorithm$AD_Stan_chunked      <- sampling_chain_counts
      paper1_settings$n_chains_by_algorithm$AD_Stan_tape_chunked <- sampling_chain_counts
      paper1_settings$n_chains_by_algorithm$AD_Stan_WCP          <- WCP_chain_counts
      ## 
      paper1_settings$n_chains_by_algorithm$Mplus_standard <- sampling_chain_counts
      # paper1_settings$n_chains_by_algorithm$Mplus_WCP      <- WCP_chain_counts
      ## 2026-10-03: no 24 or 45 chains for Mplus (two-thread narrow allocations only; see above):
      paper1_settings$n_chains_by_algorithm$Mplus_WCP      <- setdiff(WCP_chain_counts, c(24, 45))
      ##
      paper1_settings$threads_per_chain_by_algorithm$MD_BayesMVP_WCP <- bayesmvp_WCP_candidates
      paper1_settings$threads_per_chain_by_algorithm$AD_Stan_WCP     <- stan_WCP_candidates
      paper1_settings$threads_per_chain_by_algorithm$Mplus_WCP       <- mplus_WCP_candidates
}
##
## Each implementation, N and chain count has its own WCP candidates; rows are never pooled across chain counts.
## BayesMVP/Stan exclude WCP > chunks; all arms exclude chains x WCP above the thread limit.
## No extra WCP=1 cases are added; the existing chunking-only arms provide those measurements.
##
## ---- Existing sampling settings and the agreed proportional BayesMVP iteration budget ------------------------------------------------
##
## These are sampling settings; the separate burn-in settings have not been imported.
##
{
      ##
      ## ---- Timing method:
      ##
      ## "two_run_difference": every configuration runs a SHORT run and then a LONG run, back to back, same seed/config.
      ## With T(n) = S + n t: t = (T_long - T_short) / (n_long - n_short), fixed per-run cost S = T_short - n_short t, and
      ## the reported time is T_long - S (what a long run costs when start-up is negligible). "single_run" = old behaviour.
      ## The *_iterations settings below are the LONG runs; the *_short_run settings only measure the fixed cost.
      ##
      paper1_settings$timing_method <-  "two_run_difference"
      ## Mplus uses the short call as an approximate overhead control: reported time = T_long - T_short.
      ## This does not assume the short call completes its requested iterations or isolate the PPPP stage.
      ## "iteration_pair" restores the older slope estimate when both iteration counts are established.
      paper1_settings$mplus_short_run_role <-  "overhead_control"
      ##
      ## ---- Untimed warm-up call before every timed pair:
      ##
      ## At N = 50,000 with 64+ chains the timed call's start-up (state allocation, thread setup, model loading) is 4-11 s
      ## and the FIRST call of a pair is colder than the second, so T_long < T_short happened and the difference went
      ## negative. One untimed 1-iteration call before the short run makes both timed calls warm. Costs one start-up
      ## per pair (about 15-20 min over the grid, almost all at N = 50,000). Cached rows are reused unchanged.
      ##
      paper1_settings$untimed_warm_up_run_before_timing <-  TRUE
      ##
      ## ---- Final long-run counts: N = 500 / 2,500 / 10,000 / 50,000.
      ##      BayesMVP: 400 / 80 / 80 / 50; standard Stan arms: 200 / 40 / 20 / 10. Short runs remain one iteration.
      ##      The 50-iteration, no-trace BayesMVP diagnostic at N = 50,000 had three-repeat CVs of about 0.3-3.7% across
      ##      the tested configurations; this does not imply that every configuration agrees within 3%.
      ##      Changed iteration counts produce new cache keys. Matching completed cases are reused when their other
      ##      measured settings/data are unchanged; reruns are not restricted to WCP arms.
      ##
      paper1_settings$bayesmvp_iterations           <-  c("500" = 400, "2500" = 80, "10000" = 80,  "50000" = 50) ## old: 200 / 40 / 10 / 10
      paper1_settings$bayesmvp_iterations_short_run <-  c("500" = 1,   "2500" = 1,  "10000" = 1,   "50000" = 1)
      ##
      ## ---- BayesMVP WCP arm: separate LONG-run counts, 2,000 / 400 / 200 / 100.
      ##      These are 5x / 5x / 2.5x / 2x the current standard BayesMVP counts, respectively.
      ##      Longer WCP calls give more timing margin than the earlier short calls. The grid includes 4/8/16 chains
      ##      on HPC and 2/4/8 on Laptop. Short runs use bayesmvp_iterations_short_run.
      ##
      paper1_settings$bayesmvp_wcp_iterations <-  c("500" = 2000, "2500" = 400, "10000" = 200,  "50000" = 100) ## old: 2000 / 400 / 100 / 50
      ##
      total_L_bayesmvp <- 10*(paper1_settings$bayesmvp_iterations + paper1_settings$bayesmvp_iterations_short_run)
      message(paste0("total_L_bayesmvp = ", paste(total_L_bayesmvp, collapse = " ")))
      ##
      paper1_settings$stan_iterations           <-  c("500" = 200, "2500" = 40, "10000" = 20,   "50000" = 10) ## old: 200 / 40 / 10 / 4
      paper1_settings$stan_iterations_short_run <-  c("500" = 1,   "2500" = 1,  "10000" = 1,    "50000" = 1)
      ##
      ## ---- Stan WCP arm: separate LONG-run counts, 800 / 200 / 80 / 20. 
      ##      These are 4x / 5x / 4x / 2x the current standard Stan counts; they differ from the BayesMVP WCP counts.
      ##      Saved HPC rates project median corrected durations of about 2.1 / 1.5 / 1.7 / 2.4 seconds at these counts;
      ##      these are projections, not per-configuration minimums. The report accounts for differing iteration budgets.
      ##      Earlier Stan-WCP runs mistakenly used standard stan_iterations; changed-count cases receive new cache keys.
      ##
      paper1_settings$stan_wcp_iterations <- c("500" = 800, "2500" = 200, "10000" = 80, "50000" = 20) ## old: 200 / 40 / 10 / 4 
      ##
      total_L_stan <- 1*(paper1_settings$stan_iterations + paper1_settings$stan_iterations_short_run)
      message(paste0("total_L_stan = ", paste(total_L_stan, collapse = " ")))
      ##
      ## Mplus: choose FBITERATIONS or BITERATIONS for the long run; the short overhead control must use BITERATIONS.
      ## Only the FBITERATIONS count must be a multiple of 100. The short BITERATIONS call remains one requested iteration.
      ## Saved-draw verification, when enabled, applies to the long run; the overhead control need not save draws.
      ##
      # paper1_settings$mplus_iteration_mode           <-  "FBITERATIONS"  ## Long run: choose "FBITERATIONS" or "BITERATIONS".
      # paper1_settings$mplus_iteration_mode           <-  "BITERATIONS"  ## Long run: choose "FBITERATIONS" or "BITERATIONS".
      paper1_settings$mplus_iteration_mode <- c("BITERATIONS", "FBITERATIONS")
      ##
      paper1_settings$mplus_iteration_mode_short_run <-  "BITERATIONS"   ## Note: "BITERATIONS" is mandatory for the short run.
      ##
      paper1_settings$mplus_biterations_minimum <-  0  ## Used only in BITERATIONS mode.
      paper1_settings$mplus_bconvergence        <-  0  ## Mplus ignores BCONVERGENCE in FBITERATIONS mode.
      paper1_settings$mplus_save_draws         <- TRUE
      paper1_settings$mplus_verify_saved_draws <- FALSE
      ##
      paper1_settings$mplus_WCP_iterations       <-  c("500" = 10000, "2500" = 2000, "10000" = 500, "50000" = 100)
      paper1_settings$mplus_iterations_short_run <-  c("500" = 1,     "2500" = 1,    "10000" = 1,   "50000" = 1)
      ##
      ## ---- Mplus_standard: its own, much SMALLER long-run counts; mplus_WCP_iterations above is used by Mplus_WCP only (it was the shared mplus_iterations Previously).
      ##      Mirrors the original PS1 script (ps_1_optimizing_N_chunks_old.R): Mplus_standard ran FBITERATIONS = 5 x N_iter and
      ##      Mplus_WCP 30 x N_iter (N_iter = 400 / 80 / 8 at N = 500 / 2,500 / 25,000), because Mplus_standard runs many
      ##      one-thread chains, each far slower per iteration than a multi-thread Mplus_WCP chain.
      ##      2000 / 400 = the original 5 x N_iter at N = 500 / 2,500; 100 / 20 follow the same 1/N scaling at N = 10,000 / 50,000.
      ##      FBITERATIONS mode needs multiples of 100 (checked in fn_paper1_iterations_for_case), so N = 50,000 would need 100 there.
      ##      Changed-count Mplus_standard cases receive new cache keys (n_iter is part of the key); Mplus_WCP keys are unchanged.
      ##
      paper1_settings$mplus_standard_iterations  <-  c("500" = 2000,  "2500" = 400,  "10000" = 200, "50000" = 100)
      ##
      ## ---- 2026-10-03: long-run counts for the narrow-WCP allocations with many chains (24-90 chains; see the
      ##      HPC grid above). They do 2-6x the work per iteration of the 4-16-chain WCP runs that the WCP counts
      ##      above were sized for, so fewer iterations still give long runs of ~2 s or more (more with more
      ##      chains);
      ##      at N <= 2,500 BayesMVP keeps its WCP counts, because shorter long runs there would last < 1 s.
      ##      WCP-only cases (N_chunks = N_threads/chain: the slowest, memory-bound ones) use half at N >= 10,000.
      ##      The 4-16-chain cases keep the counts above, so their saved runs are reused.
      ##
      paper1_settings$wcp_many_chains_minimum <-  24
      paper1_settings$wcp_many_chains_iterations <-
            list( MD_BayesMVP_WCP = c("500" = 2000, "2500" = 400, "10000" = 100, "50000" = 25),
                  AD_Stan_WCP     = c("500" = 400,  "2500" = 100, "10000" = 20,  "50000" = 5),
                  Mplus_WCP       = c("500" = 1000, "2500" = 200, "10000" = 100, "50000" = 100))
      paper1_settings$wcp_many_chains_iterations_WCP_only <-
            list( MD_BayesMVP_WCP = c("500" = 2000, "2500" = 400, "10000" = 50,  "50000" = 12),
                  AD_Stan_WCP     = c("500" = 400,  "2500" = 100, "10000" = 10,  "50000" = 3))
      ##
      ## Flag (never drop) pairs whose fixed cost is below -10% of the long-run time (repeat noise is ~2-8%).
      ##
      paper1_settings$two_run_negative_fixed_cost_tolerance_fraction <-  0.10
      ##
      ## ---- Grid cells to leave OUT of the exported report (the runs are kept; dropped rows are counted in the report READ_ME).
      ##      The HPC WCP grid ran threads_per_chain = 6 with 8 chains at N = 500 but not at N = 2,500:
      ##
      paper1_settings$report_excluded_cases <-  list( list(device = "HPC", N = 500, n_chains = 8, threads_per_chain = 6))
      ##
      ## ---- Nuisance trace kept during the TIMED sampling calls. Previously this runner stored the WHOLE
      ##      nuisance block (N x T values per chain per iteration) in EVERY arm, whereas the legacy ps1/ps2 runners stored 10 for
      ##      BayesMVP. Allocating that trace and copying it back to R grows with iterations x chains (at N = 50,000, 96 chains,
      ##      50 iterations: ~11 GB and ~13 s of a 27 s call; ~22 s of 41 s with 180 chains), so the two-run difference counted
      ##      it as sampling time, inflating every per-iteration time and wiping out the SMT gain at large N (found by
      ##      alg_paper_1_SMT_diagnostic.R). Configured trace policy:
      ##        - BayesMVP arms: n_nuisance_to_track = 0 (its users do not need the nuisance trace; set in paper1_settings$bayesmvp
      ##          below, which is part of the BayesMVP arms' resume key, so those arms re-run).
      ##        - Stan arms: keep the FULL nuisance trace, because Stan always stores every parameter every iteration; dropping it
      ##          would make Stan look better than it is. Fixed inside the Stan case runner (not a setting), so the Stan cache is reused.
      ##
      ##
      paper1_settings$bayesmvp <-  list( eps_main = 0.00001,
                                         L_main = 10,
                                         randomize_tau   = FALSE,
                                         diffusion_HMC   = FALSE,
                                         partitioned_HMC = FALSE,
                                         metric_shape_main = "dense",
                                         init_u_value = 0.01,
                                         init_prevalence_raw = -0.6931472,
                                         n_nuisance_to_track = 0)   ## 0 = no nuisance trace in the timed calls (see above)
      ##
      ##
      ## ---- Stan trajectory length for the two-run study:
      ##
      ## These values are read by fn_paper1_stan_via_NicoStan_settings(), then checked in C++.
      ## tau = fixed_L * step_size; neither L nor tau is randomised with randomize_tau = FALSE.
      paper1_settings$stan_via_NicoStan <-  list( step_size = 0.00001,
                                                   fixed_L = 1,
                                                   randomize_tau = FALSE,
                                                   metric_shape_main = "diag",
                                                   init_radius = 0.01)
      ##
      ## Separate controls for the optional CmdStanR backend; they do not set NicoStan's path length.
      paper1_settings$stan <-  list( step_size = 0.00001,
                                     max_treedepth = 1, ## L = 2^1 - 1 = 2 - 1 = 1.
                                     metric = "diag_e",
                                     init = 0.01)
      ##
      paper1_settings$stan_backend <-  "NicoStan"
      ## Add missing WCP-only baselines with target n_chunks = n_WCP; overlapping cases are measured once.
      paper1_settings$include_NicoStan_WCP_only <-  TRUE
      ## Same for BayesMVP, so the within-BayesMVP comparison has a WCP-only optimum at every WCP count and N.
      paper1_settings$include_BayesMVP_WCP_only <-  TRUE
      ##
      ## Mplus priors retained from the original PS2 study; executable NULL means auto-detect.
      paper1_settings$mplus$settings <-  list(  prior_IW_nd = 24,
                                                prior_IW_d = 9,
                                                prior_prev_alpha = 5,
                                                prior_prev_beta = 10,
                                                prior_beta_mean_test1_nd = -2.33,
                                                prior_beta_mean_test1_d = 0.385,
                                                prior_beta_sd_test1_nd = 0.5,
                                                prior_beta_sd_test1_d = 0.45,
                                                Mplus_command = NULL)
      ##
      paper1_settings$force_recompile <-  FALSE
}
##
## ---- Compiled Stan models: one cache folder per machine, reused by every study --------------------------------------------------------
##
## Models are recompiled when the .stan source, compiler arguments, BridgeStan version or memory-release runtime changes.
##
paper1_settings$stan_model_cache_dir <-  file.path( algorithm_study_dir,
                                                    "paper_1_chunking_and_parallel_scalability",
                                                    "stan_model_cache",
                                                    Sys.info()[["nodename"]])
##
## ---- Stable output directory: later invocations resume matching completed runs in this folder -----------------------------------------
##
paper1_settings$output_dir <-  file.path( algorithm_study_dir,
                                          "paper_1_chunking_and_parallel_scalability",
                                          "paper_1_computational_outputs", device)
##
## ---- Reload the helper files right before running  ------------------------------------------------------------------------
##
## Running only part of this file (e.g. selecting from the settings blocks down) skips the source() calls in the top block, so R
## would keep whatever helper functions it loaded earlier, even if the files on disk have changed since. Sourcing them again here
## makes every run use the current helper code. The helper files only define functions, so this is quick and has no side effects.
##
{
      for (paper1_helper_file_name in c("R_fns_alg_paper_1_chunking_WCP_par_scaling.R",
                                        "R_fns_alg_paper_1_figures_tables.R",
                                        "R_fns_alg_paper_1_Mplus.R")) {

            source( file  = file.path(algorithm_study_dir,
                                      "paper_1_chunking_and_parallel_scalability",
                                      paper1_helper_file_name),
                    local = TRUE)

      }
      ##
      ## The top block stored the Mplus functions inside the settings, so point them at the freshly loaded versions too.
      ##
      paper1_settings$mplus$runner   <-  fn_paper1_run_Mplus
      paper1_settings$mplus$verifier <-  fn_paper1_verify_Mplus_run
}
##
## ---- Preview the complete grid, then execute it once ----------------------------------------------------------------------------------
##
paper1_plan <-  fn_run_paper1_benchmark( settings = paper1_settings,
                                         dry_run = TRUE,
                                         Stan_model_obj = Stan_model_obj)
##
message(paste0("\033[36mOutput directory: ", paper1_settings$output_dir, "\033[0m"))
##
##
## ---- Run each algorithm in its own fresh R process  ----------------------------------------------------------------------
##
## The updated Linux model libraries release unused autodiff arena pages after each completed sampling chain. In the tested
## N = 50,000, 180-chain case, the peak remains about 186 GiB, while retained RAM falls to about 17.5 GiB after removing the fit
## output and running gc(). The peak is similar to the matched CmdStan test; it is needed while those chains are sampling.
## Keeping one child process per algorithm also releases any remaining worker/helper-thread memory when that algorithm ends.
## Process startup is outside the timers, but memory allocation and release can affect the timed calls. The short/long method
## did not cancel that cost exactly in the paired diagnostic. Each child reuses completed cases; existing repeats remain saved.
## Afterwards this session assembles results.rds from the cache and makes the figures and tables.
##
run_each_algorithm_in_fresh_R_process <-  TRUE
##
if (isTRUE(x = run_benchmark) && isTRUE(x = run_each_algorithm_in_fresh_R_process) && !nzchar(Sys.getenv("PAPER1_CHILD_ALGORITHM"))) {

      runner_file <-  file.path(algorithm_study_dir, "paper_1_chunking_and_parallel_scalability", "alg_paper_1_chunking_WCP_par_scaling.R")
      ##
      ## ---- Each child runs a snapshot copy of this runner: Rscript reads its file as it goes, so saving an
      ## ---- edit to the runner while a child was running made that child read a mix of old and new text (parse error
      ## ---- "unexpected ')'" after the AD_Stan_tape_chunked cases had finished).
      ##
      runner_snapshot_file <-  file.path(tempdir(), paste0("alg_paper_1_runner_snapshot_", format(x = Sys.time(), format = "%Y%m%d_%H%M%S"), ".R"))
      if (!file.copy(from = runner_file, to = runner_snapshot_file, overwrite = TRUE)) stop("Could not copy the runner to ", runner_snapshot_file)
      message(paste0("\033[36mChildren run a snapshot of the runner: ", runner_snapshot_file, "\033[0m"))
      Rscript_path <-  file.path(R.home("bin"), "Rscript")
      ##
      for (child_algorithm in paper1_settings$algorithms) {

            message(paste0("\033[36m==== Starting a fresh R process for ", child_algorithm, " ====\033[0m"))
            child_exit_status <-  system2( command = Rscript_path,
                                           args    = c("--vanilla", shQuote(runner_snapshot_file)),
                                           env     = paste0("PAPER1_CHILD_ALGORITHM=", child_algorithm),
                                           stdout  = "",
                                           stderr  = "")
            if (!identical(as.integer(child_exit_status), 0L)) {
                  stop(paste0("The R process for ", child_algorithm, " failed (exit status ", child_exit_status, "). ",
                              "Its finished cases are cached; fix the error and Source again to continue."))
            }
            message(paste0("\033[36m==== ", child_algorithm, " finished; its R process has exited and freed its memory ====\033[0m"))

      }
      message("\033[36mAll algorithms measured in fresh R processes. This session now assembles results.rds from the cache (no sampling).\033[0m")

}
##
if (isTRUE(x = run_benchmark)) {

      ## Check reporting dependencies before starting the expensive sampling work.
      if (isTRUE(x = make_figures_and_tables)) fn_paper1_report_dependencies()
      ##
      paper1_results <-  fn_run_paper1_benchmark( settings = paper1_settings,
                                                  dry_run = FALSE,
                                                  Stan_model_obj = Stan_model_obj,
                                                  cache_only = isTRUE(x = run_each_algorithm_in_fresh_R_process) &&
                                                               !nzchar(Sys.getenv("PAPER1_CHILD_ALGORITHM")))
      ##
      paper1_analysis <-  paper1_results$analysis
      ##
      message(paste0("\033[36mStudy finished. Results and analysis saved in: ", paper1_results$output_dir, "\033[0m"))
      ##
      if (isTRUE(x = make_figures_and_tables)) {
  
          paper1_report <-  fn_paper1_export_manuscript( study_output_dirs = paper1_results$output_dir,
                                                     manuscript_dir = file.path(dirname(algorithm_study_dir), "Alg_papers_LaTeX",
                                                                                  "paper_1_chunking_and_parallel_scalability"),
                                                         helper_dir = file.path(algorithm_study_dir, "paper_1_chunking_and_parallel_scalability"),
                                                         excluded_cases = paper1_settings$report_excluded_cases)
  
      }

} else {

      message("Preview only. Set run_benchmark <- TRUE at the top of this runner and Source it to start the study.")

}
##
## ---- Inspect after completion; these commands do not rerun any fits -------------------------------------------------------------------
##
## View(paper1_plan$grid)
## View(paper1_results$results)
## View(paper1_analysis$best_chunks)
## View(paper1_analysis$optimal_combinations)
## View(paper1_analysis$chunking)
## View(paper1_analysis$wcp_matched)
## View(paper1_analysis$scaling)






















