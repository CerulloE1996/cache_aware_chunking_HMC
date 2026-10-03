##
## ==============================================================================================================
## analyse_overnight_counters.R
##
## Summarises the hardware-counter passes of 3 Oct 2026 on the local-HPC (AMD EPYC 9654, Zen 4), N = 50,000,
## 180 chains with one thread each (the E4 settings), for Mplus standard (BITERATIONS) and the four R workloads:
##   1. top-down pass (pmc_stat_topdown, PMC_EVENT_SET=topdown): dispatch slots lost to the back end
##      (6 slots per cycle), and the share of those due to loads not completing (memory-bound part);
##   2. floating-point pass (PMC_EVENT_SET=fp): floating-point operations retired (one per instruction, whatever
##      its width), per cycle and as a share of all retired operations;
##   3. floating-point width pass (PMC_EVENT_SET=fpwidth): floating-point uops by width (scalar, 128-, 256- and
##      512-bit; these include the loads and stores of floating-point registers);
##   4. AVX-512 vs AVX2-only builds of BayesMVP (same source), same chunked cases (two-run difference timing).
## R workloads use the window between snapshots 2 and 3 (the long sampling run); Mplus uses its whole run.
##
{
      mechanism_dir <-  path.expand(paste0("~/Documents/Work/PhD_work/Alg_paper_analysis/",
                                           "paper_1_chunking_and_parallel_scalability/mechanism_study"))
      topdown_dir   <-  file.path(mechanism_dir, "topdown_study_2026_10_03")
      avx_dir       <-  file.path(mechanism_dir, "overnight_2026_10_03", "avx_test")
      output_dir    <-  file.path(mechanism_dir, "overnight_2026_10_03")
      ##
      workloads <-  c("Mplus_standard",
                      "MD_BayesMVP_chunks1",
                      "MD_BayesMVP_chunking_500",
                      "AD_Stan_chunks1",
                      "AD_Stan_tape_chunked_250")
      event_names <-  list( topdown = c("cycles", "ex_ret_ops", "no_dispatch_backend", "no_dispatch_frontend",
                                        "no_retire_not_complete", "no_retire_load_not_complete"),
                            fp      = c("cycles", "instructions", "ex_ret_ops", "fp_ops_retired_all",
                                        "no_dispatch_smt_contention", "no_dispatch_backend"),
                            fpwidth = c("cycles", "fp_uops_scalar", "fp_uops_128", "fp_uops_256", "fp_uops_512",
                                        "fp_ops_retired_all"))
      n_threads <-  180
      n_cores   <-  96
}
##
## ---- Counter window of one workload in one pass ----------------------------------------------------------
##
fn_counter_window <-  function( counts_table,
                                workload,
                                columns
) {

        if (workload == "Mplus_standard") {
              rows <-  counts_table[counts_table$label == workload, , drop = FALSE]
              if (nrow(rows) == 0) return(NULL)
              return(list(seconds = rows$seconds[nrow(rows)], values = unlist(rows[nrow(rows), columns])))
        }
        snapshot_2 <-  counts_table[counts_table$label == paste0(workload, ":snapshot_2"), , drop = FALSE]
        snapshot_3 <-  counts_table[counts_table$label == paste0(workload, ":snapshot_3"), , drop = FALSE]
        if (nrow(snapshot_2) == 0 || nrow(snapshot_3) == 0) return(NULL)
        ## the last complete pass of each workload is used (earlier, interrupted passes are ignored):
        snapshot_2 <-  snapshot_2[nrow(snapshot_2), , drop = FALSE]
        snapshot_3 <-  snapshot_3[nrow(snapshot_3), , drop = FALSE]
        list(seconds = snapshot_3$seconds - snapshot_2$seconds,
             values  = unlist(snapshot_3[1, columns]) - unlist(snapshot_2[1, columns]))

}
##
fn_read_pass <-  function(event_set) {

        counts_file <-  file.path(topdown_dir, paste0("counts_", event_set, ".csv"))
        if (!file.exists(counts_file)) return(NULL)
        counts_table <-  utils::read.csv(file = counts_file, header = FALSE, stringsAsFactors = FALSE)
        names(counts_table) <-  c("label", "seconds", "exit_status", "max_rss_kb", event_names[[event_set]])
        windows <-  lapply(X = workloads, FUN = function(workload) {
              fn_counter_window(counts_table = counts_table, workload = workload, columns = event_names[[event_set]])
        })
        names(windows) <-  workloads
        windows[!vapply(X = windows, FUN = is.null, FUN.VALUE = logical(1))]

}
##
## ---- 1. Top-down ----------------------------------------------------------------------------------------------
##
{
      topdown_windows <-  fn_read_pass(event_set = "topdown")
      topdown_summary <-  do.call(what = rbind, args = lapply(X = names(topdown_windows), FUN = function(workload) {
            v <-  topdown_windows[[workload]]$values
            slots <-  6 * v[["cycles"]]
            backend_share <-  v[["no_dispatch_backend"]] / slots
            memory_share_of_backend <-  v[["no_retire_load_not_complete"]] / v[["no_retire_not_complete"]]
            data.frame( workload                         = workload,
                        ops_per_cycle                    = round(v[["ex_ret_ops"]] / v[["cycles"]], 3),
                        backend_bound_percent            = round(100 * backend_share, 1),
                        memory_bound_percent             = round(100 * backend_share * memory_share_of_backend, 1),
                        core_bound_percent               = round(100 * backend_share *
                                                                 (1 - memory_share_of_backend), 1),
                        frontend_bound_percent           = round(100 * v[["no_dispatch_frontend"]] / slots, 1),
                        stringsAsFactors                 = FALSE)
      }))
      message(paste0("\033[36m", "Top-down (share of dispatch slots; local-HPC, N = 50,000, 180 chains):", "\033[0m"))
      print(topdown_summary)
}
##
## ---- 2. Floating-point operations ------------------------------------------------------------------------------
##
{
      fp_windows <-  fn_read_pass(event_set = "fp")
      fp_summary <-  do.call(what = rbind, args = lapply(X = names(fp_windows), FUN = function(workload) {
            v <-  fp_windows[[workload]]$values
            fp_per_thread_cycle <-  v[["fp_ops_retired_all"]] / v[["cycles"]]
            data.frame( workload                          = workload,
                        instructions_per_cycle            = round(v[["instructions"]] / v[["cycles"]], 3),
                        fp_ops_per_thread_cycle           = round(fp_per_thread_cycle, 3),
                        ## both SMT threads of a core share its floating-point pipes (180 threads on 96 cores):
                        fp_ops_per_core_cycle             = round(fp_per_thread_cycle * n_threads / n_cores, 3),
                        fp_share_of_ops_percent           = round(100 * v[["fp_ops_retired_all"]] / v[["ex_ret_ops"]], 1),
                        fp_ops_per_second_billions        = round(v[["fp_ops_retired_all"]] /
                                                                  fp_windows[[workload]]$seconds / 1e9, 1),
                        smt_contention_percent            = round(100 * v[["no_dispatch_smt_contention"]] /
                                                                  (6 * v[["cycles"]]), 1),
                        stringsAsFactors                  = FALSE)
      }))
      message(paste0("\033[36m", "Floating-point operations retired (one per instruction):", "\033[0m"))
      print(fp_summary)
}
##
## ---- 3. Floating-point width --------------------------------------------------------------------------------------
##
{
      fpwidth_windows <-  fn_read_pass(event_set = "fpwidth")
      fpwidth_summary <-  do.call(what = rbind, args = lapply(X = names(fpwidth_windows), FUN = function(workload) {
            v <-  fpwidth_windows[[workload]]$values
            width_uops <-  c(scalar = v[["fp_uops_scalar"]], bits_128 = v[["fp_uops_128"]],
                             bits_256 = v[["fp_uops_256"]], bits_512 = v[["fp_uops_512"]])
            data.frame( workload                    = workload,
                        scalar_percent              = round(100 * width_uops[["scalar"]]   / sum(width_uops), 1),
                        bits_128_percent            = round(100 * width_uops[["bits_128"]] / sum(width_uops), 1),
                        bits_256_percent            = round(100 * width_uops[["bits_256"]] / sum(width_uops), 1),
                        bits_512_percent            = round(100 * width_uops[["bits_512"]] / sum(width_uops), 1),
                        fp_uops_per_core_cycle      = round(sum(width_uops) / v[["cycles"]] * n_threads / n_cores, 3),
                        stringsAsFactors            = FALSE)
      }))
      message(paste0("\033[36m", "Floating-point uops by width (incl. loads/stores of FP registers):", "\033[0m"))
      print(fpwidth_summary)
}
##
## ---- 4. AVX-512 vs AVX2-only builds of BayesMVP ----------------------------------------------------------------------
##
{
      avx_summary <-  NULL
      avx_times_file <-  file.path(avx_dir, "times.csv")
      if (file.exists(avx_times_file)) {
            avx_times <-  utils::read.csv(file = avx_times_file, stringsAsFactors = FALSE)
            avx_times$build <-  sub(pattern = "_.*$", replacement = "", x = avx_times$label)
            avx_times$seconds_per_iteration <-  (avx_times$elapsed_seconds_long_run - avx_times$elapsed_seconds_short_run) /
                                                (avx_times$n_iter_long - avx_times$n_iter_short)
            avx_times$chain_iterations_per_second <-  avx_times$n_chains / avx_times$seconds_per_iteration
            avx_summary <-  stats::aggregate( chain_iterations_per_second ~ build + n_chains + num_chunks,
                                              data = avx_times, FUN = mean)
            avx_wide <-  stats::reshape( data = avx_summary, idvar = c("n_chains", "num_chunks"),
                                         timevar = "build", direction = "wide")
            avx_wide$AVX512_over_AVX2 <-  round(avx_wide$chain_iterations_per_second.avx512 /
                                                avx_wide$chain_iterations_per_second.avx2, 3)
            message(paste0("\033[36m", "BayesMVP, AVX-512 vs AVX2-only build (chain-iterations per second; ",
                           "mean of the repeats):", "\033[0m"))
            print(avx_times[, c("label", "n_chains", "num_chunks", "seconds_per_iteration",
                                "chain_iterations_per_second")])
            print(avx_wide)
            avx_summary <-  avx_wide
      }
}
##
## ---- Save ------------------------------------------------------------------------------------------------------
##
{
      utils::write.csv(x = topdown_summary, file = file.path(output_dir, "summary_topdown.csv"), row.names = FALSE)
      utils::write.csv(x = fp_summary,      file = file.path(output_dir, "summary_fp.csv"),      row.names = FALSE)
      if (!is.null(fpwidth_summary)) {
            utils::write.csv(x = fpwidth_summary, file = file.path(output_dir, "summary_fpwidth.csv"), row.names = FALSE)
      }
      if (!is.null(avx_summary)) {
            utils::write.csv(x = avx_summary, file = file.path(output_dir, "summary_avx.csv"), row.names = FALSE)
      }
      message(paste0("\033[32m", "Saved the summaries in ", output_dir, "\033[0m"))
}























