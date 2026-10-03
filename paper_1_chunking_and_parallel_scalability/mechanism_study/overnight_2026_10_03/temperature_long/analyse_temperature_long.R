##
## ==============================================================================================================
## analyse_temperature_long.R
##
## Summarises the long CPU temperature runs (run_overnight_measurements.sh, part 1) on the local-HPC
## (AMD EPYC 9654), N = 50,000, 180 chains with one thread each (the E4 settings). Each workload ran for ~600 s
## of sampling after the CPU had cooled to within 2 C (Tctl) of the idle baseline; temp_power_log recorded Tctl,
## the twelve per-CCD temperatures, the CPU package energy (RAPL) and the pages swapped in/out every second.
##
## Per workload: plateau = mean over the last 120 s of the run; peak = maximum over the run; rise = plateau minus
## the idle baseline (the 300 s idle period at the start). Package power = energy difference / time over the same
## 120 s. Counts (pmc_stat, umc_stat) cover the loaded period (snapshot 1 to 2) of each R workload and the whole
## Mplus run. Seconds per iteration (from the sampling calls) are compared with those of E4.
##
{
      temperature_dir <-  path.expand(paste0("~/Documents/Work/PhD_work/Alg_paper_analysis/",
                                             "paper_1_chunking_and_parallel_scalability/mechanism_study/",
                                             "overnight_2026_10_03/temperature_long"))
      workloads <-  c("Mplus_standard",
                      "MD_BayesMVP_chunking_500",
                      "AD_Stan_tape_chunked_250",
                      "MD_BayesMVP_chunks1",
                      "AD_Stan_chunks1")
      workload_display_names <-  c(Mplus_standard           = "Mplus standard (BITERATIONS)",
                                   MD_BayesMVP_chunking_500 = "MD_BayesMVP_chunking (N_chunks = 500)",
                                   AD_Stan_tape_chunked_250 = "AD_Stan_tape_chunked (N_chunks = 250)",
                                   MD_BayesMVP_chunks1      = "MD_BayesMVP (N_chunks = 1)",
                                   AD_Stan_chunks1          = "AD_Stan (N_chunks = 1)")
      ## E4 seconds per iteration at the same settings (mechanism_times_HPC.csv, two-run difference):
      E4_seconds_per_iteration <-  c(MD_BayesMVP_chunking_500 = 0.527,
                                     AD_Stan_tape_chunked_250 = 1.818,
                                     MD_BayesMVP_chunks1      = 6.730,
                                     AD_Stan_chunks1          = 4.474)
      plateau_seconds <-  120
      counter_names <-  c("cycles", "instructions",
                          "fills_local_L2", "fills_local_L3", "fills_other_CCX", "fills_DRAM")
}
##
## ---- Read the logs ---------------------------------------------------------------------------------------
##
{
      temperatures <-  utils::read.csv(file = file.path(temperature_dir, "temperatures.csv"),
                                       stringsAsFactors = FALSE)
      ccd_columns <-  grep(pattern = "^Tccd", x = names(temperatures), value = TRUE)
      temperatures$Tccd_mean <-  rowMeans(temperatures[, ccd_columns])
      temperatures$Tccd_max  <-  apply(X = temperatures[, ccd_columns], MARGIN = 1, FUN = max)
      ##
      counts <-  utils::read.csv(file = file.path(temperature_dir, "counts.csv"),
                                 header = FALSE, stringsAsFactors = FALSE)
      names(counts) <-  c("label", "seconds", "exit_status", "max_rss_kb", counter_names)
      umc <-  utils::read.csv(file = file.path(temperature_dir, "umc.csv"),
                              header = FALSE, stringsAsFactors = FALSE)
      names(umc) <-  c("label", "seconds", "exit_status", "DRAM_read_bytes", "DRAM_write_bytes")
      ## calls.csv exists once the first R workload has finished:
      calls <-  data.frame(label = character(0), n_iter = numeric(0), elapsed_seconds_sampler = numeric(0))
      if (file.exists(file.path(temperature_dir, "calls.csv"))) {
            calls <-  utils::read.csv(file = file.path(temperature_dir, "calls.csv"), stringsAsFactors = FALSE)
      }
}
##
## ---- Counter window: loaded period (snapshot 1 to 2) for R, whole run for Mplus ------------------------------
##
fn_counter_window <-  function( table,
                                workload,
                                columns
) {

        if (workload == "Mplus_standard") {
              end_row <-  table[table$label == workload, , drop = FALSE]
              return(list(seconds = end_row$seconds[1], values = unlist(end_row[1, columns])))
        }
        snapshot_1 <-  table[table$label == paste0(workload, ":snapshot_1"), , drop = FALSE]
        snapshot_2 <-  table[table$label == paste0(workload, ":snapshot_2"), , drop = FALSE]
        list(seconds = snapshot_2$seconds[1] - snapshot_1$seconds[1],
             values  = unlist(snapshot_2[1, columns]) - unlist(snapshot_1[1, columns]))

}
##
## ---- Per-workload summary ----------------------------------------------------------------------------------
##
{
      idle_rows <-  temperatures[temperatures$phase == "idle_baseline", , drop = FALSE]
      idle_Tctl <-  mean(idle_rows$Tctl)
      idle_Tccd_mean <-  mean(idle_rows$Tccd_mean)
      idle_power_watts <-  diff(range(idle_rows$package_energy_joules)) / diff(range(idle_rows$unix_time))
      ##
      summary_rows <-  lapply(X = workloads, FUN = function(workload) {

            run_rows <-  temperatures[temperatures$phase == workload, , drop = FALSE]
            if (nrow(run_rows) == 0) return(NULL)
            plateau_start <-  max(run_rows$unix_time) - plateau_seconds
            plateau_rows <-  run_rows[run_rows$unix_time >= plateau_start, , drop = FALSE]
            plateau_power_watts <-  diff(range(plateau_rows$package_energy_joules)) /
                                    diff(range(plateau_rows$unix_time))
            ## Power while sampling: the R workloads pause for ~5-8 s between sampling calls (R-side setup and
            ## conversion), so the median of the per-second power over the run is used (sampling fills most seconds):
            per_second_power <-  diff(run_rows$package_energy_joules) / diff(run_rows$unix_time)
            sampling_power_watts <-  stats::median(per_second_power[-(1:30)])
            ## Tctl in the seconds at (near) the sampling power, over the plateau window:
            plateau_power_per_second <-  c(NA, diff(plateau_rows$package_energy_joules) / diff(plateau_rows$unix_time))
            sampling_seconds <-  !is.na(plateau_power_per_second) & plateau_power_per_second >= 0.9 * sampling_power_watts
            sampling_Tctl <-  mean(plateau_rows$Tctl[sampling_seconds])
            ## Slope of Tctl over the plateau window (C per minute; close to 0 once the plateau is reached):
            ## (time centred first: raw unix times are collinear with the intercept at lm's tolerance)
            plateau_rows$minutes_from_start <-  (plateau_rows$unix_time - min(plateau_rows$unix_time)) / 60
            plateau_slope <-  stats::coef(stats::lm(Tctl ~ minutes_from_start, data = plateau_rows))[[2]]
            ##
            core_window <-  fn_counter_window(table = counts, workload = workload, columns = counter_names)
            dram_window <-  fn_counter_window(table = umc, workload = workload,
                                              columns = c("DRAM_read_bytes", "DRAM_write_bytes"))
            workload_calls <-  calls[calls$label == workload, , drop = FALSE]
            seconds_per_iteration <-  NA
            if (nrow(workload_calls) > 0) {
                  seconds_per_iteration <-  sum(workload_calls$elapsed_seconds_sampler) /
                                            sum(workload_calls$n_iter)
            }
            ##
            data.frame( workload                    = workload,
                        run_seconds                 = round(diff(range(run_rows$unix_time)), 0),
                        idle_Tctl                   = round(idle_Tctl, 1),
                        plateau_Tctl                = round(mean(plateau_rows$Tctl), 1),
                        peak_Tctl                   = round(max(run_rows$Tctl), 1),
                        rise_Tctl                   = round(mean(plateau_rows$Tctl) - idle_Tctl, 1),
                        plateau_slope_C_per_minute  = round(plateau_slope, 2),
                        plateau_Tccd_mean           = round(mean(plateau_rows$Tccd_mean), 1),
                        rise_Tccd_mean              = round(mean(plateau_rows$Tccd_mean) - idle_Tccd_mean, 1),
                        peak_Tccd                   = round(max(run_rows$Tccd_max), 1),
                        idle_power_watts            = round(idle_power_watts, 0),
                        plateau_power_watts         = round(plateau_power_watts, 0),
                        sampling_power_watts        = round(sampling_power_watts, 0),
                        sampling_Tctl               = round(sampling_Tctl, 1),
                        pages_swapped_in            = diff(range(run_rows$pswpin)),
                        pages_swapped_out           = diff(range(run_rows$pswpout)),
                        instructions_per_cycle      = round(core_window$values[["instructions"]] /
                                                            core_window$values[["cycles"]], 3),
                        DRAM_GB_per_second          = round((dram_window$values[["DRAM_read_bytes"]] +
                                                             dram_window$values[["DRAM_write_bytes"]]) /
                                                            dram_window$seconds / 1e9, 1),
                        seconds_per_iteration       = round(seconds_per_iteration, 3),
                        E4_seconds_per_iteration    = if (workload %in% names(E4_seconds_per_iteration)) {
                                                            E4_seconds_per_iteration[[workload]]
                                                      } else NA,
                        stringsAsFactors            = FALSE)

      })
      temperature_long_summary <-  do.call(what = rbind, args = summary_rows)
      utils::write.csv(x = temperature_long_summary, row.names = FALSE,
                       file = file.path(temperature_dir, "temperature_long_summary.csv"))
      message(paste0("\033[36m", "Long temperature runs (local-HPC, N = 50,000, 180 chains):", "\033[0m"))
      print(temperature_long_summary)
}
##
## ---- Figure: Tctl and package power over time, workloads shaded -----------------------------------------------
##
{
      plot_data <-  temperatures
      plot_data$minutes <-  (plot_data$unix_time - min(plot_data$unix_time)) / 60
      plot_data$segment <-  cumsum(c(1, diff(plot_data$unix_time) > 5))
      plot_data$power_watts <-  c(NA, diff(plot_data$package_energy_joules) / diff(plot_data$unix_time))
      run_bounds <-  do.call(what = rbind, args = lapply(X = workloads, FUN = function(workload) {
            rows <-  plot_data[plot_data$phase == workload, , drop = FALSE]
            if (nrow(rows) == 0) return(NULL)
            data.frame(workload = workload, start = min(rows$minutes), end = max(rows$minutes),
                       stringsAsFactors = FALSE)
      }))
      run_bounds$workload <-  factor(workload_display_names[run_bounds$workload],
                                     levels = workload_display_names)
      long_data <-  rbind(data.frame(minutes = plot_data$minutes, value = plot_data$Tctl,
                                     segment = plot_data$segment, panel = "Tctl (°C)"),
                          data.frame(minutes = plot_data$minutes, value = plot_data$power_watts,
                                     segment = plot_data$segment, panel = "CPU package power (W)"))
      ##
      temperature_plot <-  ggplot2::ggplot() +
                           ggplot2::geom_rect(data = run_bounds,
                                              mapping = ggplot2::aes(xmin = start, xmax = end,
                                                                     ymin = -Inf, ymax = Inf, fill = workload),
                                              alpha = 0.25) +
                           ggplot2::geom_line(data = long_data,
                                              mapping = ggplot2::aes(x = minutes, y = value, group = segment)) +
                           ggplot2::facet_wrap(facets = ggplot2::vars(panel), ncol = 1, scales = "free_y") +
                           ggplot2::labs(x = "Time (minutes)", y = NULL, fill = "Workload",
                                         title = paste0("local-HPC, N = 50,000, 180 chains: CPU temperature ",
                                                        "and power per workload (600 s each)")) +
                           ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2)) +
                           ggplot2::theme_bw(base_size = 12) +
                           ggplot2::theme(legend.position = "bottom")
      ggplot2::ggsave(filename = file.path(temperature_dir, "temperature_long_HPC.png"),
                      plot = temperature_plot, width = 11, height = 7.5, dpi = 150)
      message(paste0("\033[32m", "Saved temperature_long_summary.csv and temperature_long_HPC.png", "\033[0m"))
}
























