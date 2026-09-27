##
## =====================================================================================================================================
## analyse_mechanism_study.R
##
## Summarises the mechanism experiments (run_mechanism_experiments.sh). For each case, pmc_stat took three counter snapshots around
## two sampling calls in the same R process (a short run and a long run, after an untimed warm-up call), so
##
##     per-iteration count = (long-run window - short-run window) / (n_iter_long - n_iter_short),
##
## which removes the fixed per-call costs in the same way as the study's two-run timing. Counts cover every thread of the process
## (all chains). "Per chain-iteration" divides by the number of chains. Fills are L1 data-cache fills, by where the data came from.
##
{
      mechanism_dir <-  path.expand("~/Documents/Work/PhD_work/Alg_paper_analysis/paper_1_chunking_and_parallel_scalability/mechanism_study")
      devices <-  c("HPC", "Laptop")
      counter_names <-  c("cycles", "instructions", "fills_local_L2", "fills_local_L3", "fills_other_CCX", "fills_DRAM")
}
##
## ---- Per-case summary ------------------------------------------------------------------------------------------------------------------
##
fn_mechanism_summary <-  function( counts_file,
                                   times_file
) {

        counts <-  utils::read.csv(file = counts_file, stringsAsFactors = FALSE)
        times <-  utils::read.csv(file = times_file, stringsAsFactors = FALSE)
        ##
        snapshot_rows <-  grepl(pattern = ":snapshot_[0-9]+$", x = counts$label)
        snapshots <-  counts[snapshot_rows, , drop = FALSE]
        snapshots$case_label <-  sub(pattern = ":snapshot_[0-9]+$", replacement = "", x = snapshots$label)
        snapshots$snapshot <-  as.numeric(sub(pattern = "^.*:snapshot_", replacement = "", x = snapshots$label))
        ##
        summary_rows <-  lapply(X = seq_len(nrow(times)), FUN = function(row_index) {

                case_row <-  times[row_index, , drop = FALSE]
                case_snapshots <-  snapshots[snapshots$case_label == case_row$label, , drop = FALSE]
                case_snapshots <-  case_snapshots[order(case_snapshots$snapshot), , drop = FALSE]
                if (nrow(case_snapshots) < 3) return(NULL)
                ## If a case was repeated, use its last three snapshots:
                case_snapshots <-  utils::tail(case_snapshots, 3)
                ##
                short_window <-  unlist(case_snapshots[2, counter_names]) - unlist(case_snapshots[1, counter_names])
                long_window <-  unlist(case_snapshots[3, counter_names]) - unlist(case_snapshots[2, counter_names])
                n_iter_difference <-  case_row$n_iter_long - case_row$n_iter_short
                per_iteration <-  (long_window - short_window) / n_iter_difference
                per_chain_iteration <-  per_iteration / case_row$n_chains
                ##
                seconds_per_iteration <-  (case_row$elapsed_seconds_long_run - case_row$elapsed_seconds_short_run) / n_iter_difference
                total_fills <-  sum(per_iteration[c("fills_local_L2", "fills_local_L3", "fills_other_CCX", "fills_DRAM")])
                ##
                data.frame( case_row,
                            seconds_per_iteration          = seconds_per_iteration,
                            chain_iterations_per_second    = case_row$n_chains / seconds_per_iteration,
                            instructions_per_cycle         = per_iteration[["instructions"]] / per_iteration[["cycles"]],
                            cycles_per_chain_iteration     = per_chain_iteration[["cycles"]],
                            instructions_per_chain_iter    = per_chain_iteration[["instructions"]],
                            L2_fills_per_chain_iter        = per_chain_iteration[["fills_local_L2"]],
                            L3_fills_per_chain_iter        = per_chain_iteration[["fills_local_L3"]],
                            other_CCX_fills_per_chain_iter = per_chain_iteration[["fills_other_CCX"]],
                            DRAM_fills_per_chain_iter      = per_chain_iteration[["fills_DRAM"]],
                            DRAM_MB_per_chain_iteration    = 64 * per_chain_iteration[["fills_DRAM"]] / 1e6,
                            share_of_fills_from_L3         = per_iteration[["fills_local_L3"]] / total_fills,
                            share_of_fills_from_DRAM       = per_iteration[["fills_DRAM"]] / total_fills,
                            DRAM_read_GB_per_second        = 64 * per_iteration[["fills_DRAM"]] / seconds_per_iteration / 1e9,
                            stringsAsFactors               = FALSE)

        })
        ##
        do.call(what = rbind, args = summary_rows)

}
##
{
      mechanism_summary <-  NULL
      for (device in devices) {
            counts_file <-  file.path(mechanism_dir, paste0("mechanism_counts_", device, ".csv"))
            times_file <-  file.path(mechanism_dir, paste0("mechanism_times_", device, ".csv"))
            if (!file.exists(counts_file) || !file.exists(times_file)) next
            mechanism_summary <-  rbind(mechanism_summary, fn_mechanism_summary(counts_file = counts_file, times_file = times_file))
      }
      mechanism_summary$experiment <-  sub(pattern = "_.*$", replacement = "", x = mechanism_summary$label)
      mechanism_summary$placement <-  sub(pattern = "^.*_tpc[0-9]+_", replacement = "", x = mechanism_summary$label)
      utils::write.csv(x = mechanism_summary, file = file.path(mechanism_dir, "mechanism_summary.csv"), row.names = FALSE)
      message(paste0("\033[36mMechanism summary: ", nrow(mechanism_summary), " cases written to mechanism_summary.csv\033[0m"))
}

##
## ---- Figures -------------------------------------------------------------------------------------------------------------------------
##
{
      algorithm_labels <-  c( MD_BayesMVP = "BayesMVP (manual gradients)", AD_Stan = "Stan (autodiff)", AD_Stan_tape_chunked = "Stan (autodiff)")
      ##
      E1 <-  mechanism_summary[mechanism_summary$experiment == "E1", , drop = FALSE]
      E1$implementation <-  unname(algorithm_labels[E1$algorithm])
      E1$chunk_label <-  factor( x = paste0(E1$num_chunks, " chunk", ifelse(E1$num_chunks == 1, "", "s")),
                                 levels = paste0(sort(unique(E1$num_chunks)), " chunk", ifelse(sort(unique(E1$num_chunks)) == 1, "", "s")))
      E1$N_label <-  factor(x = paste0("N = ", format(E1$N, big.mark = ",")), levels = c("N = 10,000", "N = 50,000"))
      ##
      for (device in unique(E1$device)) {

            E1_device <-  E1[E1$device == device, , drop = FALSE]
            ##
            plot_throughput <-  ggplot2::ggplot( data = E1_device,
                                                 mapping = ggplot2::aes(x = n_chains, y = chain_iterations_per_second, colour = chunk_label, group = chunk_label)) +
                                ggplot2::geom_line() +
                                ggplot2::geom_point() +
                                ggplot2::scale_x_log10(breaks = sort(unique(E1_device$n_chains))) +
                                ggplot2::scale_y_log10() +
                                ggplot2::facet_grid(rows = ggplot2::vars(N_label), cols = ggplot2::vars(implementation), scales = "free_y") +
                                ggplot2::labs( x = "Number of chains (one thread each)", y = "Throughput (chain-iterations / second)",
                                               colour = "Chunks", title = paste0(device, ": throughput")) +
                                ggplot2::theme_bw()
            ##
            plot_IPC <-  ggplot2::ggplot( data = E1_device,
                                          mapping = ggplot2::aes(x = n_chains, y = instructions_per_cycle, colour = chunk_label, group = chunk_label)) +
                         ggplot2::geom_line() +
                         ggplot2::geom_point() +
                         ggplot2::scale_x_log10(breaks = sort(unique(E1_device$n_chains))) +
                         ggplot2::facet_grid(rows = ggplot2::vars(N_label), cols = ggplot2::vars(implementation)) +
                         ggplot2::labs( x = "Number of chains (one thread each)", y = "Instructions per cycle",
                                        colour = "Chunks", title = paste0(device, ": instructions per cycle")) +
                         ggplot2::theme_bw()
            ##
            plot_DRAM <-  ggplot2::ggplot( data = E1_device,
                                           mapping = ggplot2::aes(x = n_chains, y = DRAM_MB_per_chain_iteration, colour = chunk_label, group = chunk_label)) +
                          ggplot2::geom_line() +
                          ggplot2::geom_point() +
                          ggplot2::scale_x_log10(breaks = sort(unique(E1_device$n_chains))) +
                          ggplot2::scale_y_log10() +
                          ggplot2::facet_grid(rows = ggplot2::vars(N_label), cols = ggplot2::vars(implementation), scales = "free_y") +
                          ggplot2::labs( x = "Number of chains (one thread each)", y = "Data loaded from DRAM per chain-iteration (MB)",
                                         colour = "Chunks", title = paste0(device, ": DRAM traffic per chain-iteration")) +
                          ggplot2::theme_bw()
            ##
            ggplot2::ggsave(filename = file.path(mechanism_dir, paste0("figure_E1_throughput_", device, ".png")), plot = plot_throughput, width = 9, height = 6, dpi = 150)
            ggplot2::ggsave(filename = file.path(mechanism_dir, paste0("figure_E1_IPC_", device, ".png")), plot = plot_IPC, width = 9, height = 6, dpi = 150)
            ggplot2::ggsave(filename = file.path(mechanism_dir, paste0("figure_E1_DRAM_", device, ".png")), plot = plot_DRAM, width = 9, height = 6, dpi = 150)

      }
      message("\033[36mFigures written to the mechanism_study folder.\033[0m")
}

##
## ---- DRAM traffic at the memory controllers (HPC, E5; umc_stat). Background traffic from other programs (~3 GB/s) is not subtracted.
##
{
      umc_file <-  file.path(mechanism_dir, "mechanism_umc_HPC.csv")
      if (file.exists(umc_file)) {

            umc <-  utils::read.csv(file = umc_file, header = FALSE, stringsAsFactors = FALSE,
                                    col.names = c("label", "seconds", "status", "DRAM_read_bytes", "DRAM_write_bytes"))
            umc <-  umc[grepl(pattern = ":snapshot_[0-9]+$", x = umc$label), , drop = FALSE]
            umc$case_label <-  sub(pattern = ":snapshot_[0-9]+$", replacement = "", x = umc$label)
            umc$snapshot <-  as.numeric(sub(pattern = "^.*:snapshot_", replacement = "", x = umc$label))
            times_HPC <-  utils::read.csv(file = file.path(mechanism_dir, "mechanism_times_HPC.csv"), stringsAsFactors = FALSE)
            ##
            bandwidth <-  do.call(what = rbind, args = lapply(X = unique(umc$case_label), FUN = function(case_label) {
                  snapshots_case <-  utils::tail(umc[umc$case_label == case_label, , drop = FALSE][order(umc$snapshot[umc$case_label == case_label]), ], 3)
                  case_row <-  times_HPC[times_HPC$label == case_label, , drop = FALSE]
                  if (nrow(snapshots_case) < 3 || nrow(case_row) != 1) return(NULL)
                  n_iter_difference <-  case_row$n_iter_long - case_row$n_iter_short
                  seconds_per_iteration <-  (case_row$elapsed_seconds_long_run - case_row$elapsed_seconds_short_run) / n_iter_difference
                  window_difference <-  function(column) ((snapshots_case[[column]][3] - snapshots_case[[column]][2]) -
                                                          (snapshots_case[[column]][2] - snapshots_case[[column]][1])) / n_iter_difference
                  data.frame( case_row[, c("label", "algorithm", "N", "num_chunks", "n_chains")],
                              seconds_per_iteration        = seconds_per_iteration,
                              DRAM_read_GB_per_second      = window_difference("DRAM_read_bytes") / seconds_per_iteration / 1e9,
                              DRAM_write_GB_per_second     = window_difference("DRAM_write_bytes") / seconds_per_iteration / 1e9,
                              DRAM_read_GB_per_chain_iter  = window_difference("DRAM_read_bytes") / case_row$n_chains / 1e9,
                              stringsAsFactors             = FALSE)
            }))
            bandwidth$DRAM_total_GB_per_second <-  bandwidth$DRAM_read_GB_per_second + bandwidth$DRAM_write_GB_per_second
            bandwidth$configuration <-  paste0(ifelse(bandwidth$algorithm == "MD_BayesMVP", "BayesMVP", ifelse(bandwidth$algorithm == "AD_Stan", "Stan", "Stan tape-chunked")),
                                               ", ", bandwidth$num_chunks, " chunk", ifelse(bandwidth$num_chunks == 1, "", "s"))
            utils::write.csv(x = bandwidth, file = file.path(mechanism_dir, "mechanism_DRAM_bandwidth_HPC.csv"), row.names = FALSE)
            ##
            plot_bandwidth <-  ggplot2::ggplot( data = bandwidth,
                                                mapping = ggplot2::aes(x = n_chains, y = DRAM_total_GB_per_second, colour = configuration, group = configuration)) +
                               ggplot2::geom_hline(yintercept = 460.8, linetype = "dashed", colour = "grey40") +
                               ggplot2::annotate(geom = "text", x = 8, y = 475, hjust = 0, size = 3,
                                                 label = "Theoretical peak, 12 x DDR5-4800 channels (460.8 GB/s)") +
                               ggplot2::geom_line() +
                               ggplot2::geom_point() +
                               ggplot2::scale_x_log10(breaks = c(8, 48, 96, 180)) +
                               ggplot2::labs( x = "Number of chains (one thread each)", y = "DRAM traffic, reads + writes (GB/s)",
                                              colour = NULL, title = "HPC, N = 50,000: DRAM traffic measured at the memory controllers") +
                               ggplot2::theme_bw()
            ggplot2::ggsave(filename = file.path(mechanism_dir, "figure_E5_DRAM_bandwidth_HPC.png"), plot = plot_bandwidth, width = 9, height = 5, dpi = 150)

      }
}























