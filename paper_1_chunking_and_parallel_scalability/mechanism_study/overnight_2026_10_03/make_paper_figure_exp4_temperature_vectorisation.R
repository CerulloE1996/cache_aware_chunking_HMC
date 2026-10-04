##
## ==============================================================================================================
## make_paper_figure_exp4_temperature_vectorisation.R
##
## Figure for E4, design (v) (local-HPC, N = 50,000, N_chains = 180, N_threads/chain = 1):
##   (A) CPU control temperature (Tctl; 30 s centred moving mean) and (B) CPU package power (RAPL; 31 s
##       running median, i.e., without the ~5-8 s pauses between the sampling calls), against the minutes since
##       the start of each ~600 s full-load run (per-second logs in temperature_long), with the continuous run of
##       MD_BayesMVP_chunking (one sampling call; temperature_continuous) dashed;
##   (C) time per iteration of MD_BayesMVP_chunking with BayesMVP's AVX2 maths functions, and with
##       stan::math's scalar functions, relative to BayesMVP's AVX-512 maths functions (same AVX-512 build;
##       two-run difference timing as in E4, mean of the two repeats; stan_math_test, 360 W power limit).
## Usage: Rscript make_paper_figure_exp4_temperature_vectorisation.R
##
{
      overnight_dir <-  path.expand(paste0("~/Documents/Work/PhD_work/Alg_paper_analysis/",
                                           "paper_1_chunking_and_parallel_scalability/mechanism_study/",
                                           "overnight_2026_10_03"))
      figure_file <-  file.path(overnight_dir, "Figure_exp4_temperature_power_vectorisation.png")
      moving_mean_seconds <-  30
      running_median_seconds <-  31
      ##
      ## Runs (phase names in the logs), in legend order:
      run_names <-  c("MD_BayesMVP_chunking_500",
                      "MD_BayesMVP_chunking_500_continuous",
                      "AD_Stan_tape_chunked_250",
                      "MD_BayesMVP_chunks1",
                      "AD_Stan_chunks1",
                      "Mplus_standard")
      run_legend_labels <-  c(
            expression(paste("MD_BayesMVP_chunking (", N[chunks], " = 500)")),
            expression(paste("MD_BayesMVP_chunking (", N[chunks], " = 500; one continuous sampling call)")),
            expression(paste("AD_Stan_tape_chunked (", N[chunks], " = 250)")),
            expression(paste("MD_BayesMVP (", N[chunks], " = 1)")),
            expression(paste("AD_Stan (", N[chunks], " = 1)")),
            "Mplus (BITERATIONS)")
      run_colours <-  c(MD_BayesMVP_chunking_500            = "#D55E00",
                        MD_BayesMVP_chunking_500_continuous = "#D55E00",
                        AD_Stan_tape_chunked_250            = "#0072B2",
                        MD_BayesMVP_chunks1                 = "#E69F00",
                        AD_Stan_chunks1                     = "#56B4E9",
                        Mplus_standard                      = "#009E73")
      run_linetypes <-  c(MD_BayesMVP_chunking_500            = "solid",
                          MD_BayesMVP_chunking_500_continuous = "22",
                          AD_Stan_tape_chunked_250            = "solid",
                          MD_BayesMVP_chunks1                 = "solid",
                          AD_Stan_chunks1                     = "solid",
                          Mplus_standard                      = "solid")
      ##
      ## Maths-function variants (the vect_type of the timed runs), in legend order:
      variant_names <-  c("AVX512", "AVX2", "Stan")
      variant_legend_labels <-  c("BayesMVP's AVX-512 functions (reference)",
                                  "BayesMVP's AVX2 functions",
                                  "stan::math's scalar functions")
      variant_colours <-  c(AVX512 = "#4D4D4D", AVX2 = "#999999", Stan = "#CC79A7")
}
##
## ---- Per-second traces of one run, aligned at the start of the run -----------------------------------------
##
fn_moving_mean <-  function( values,
                             window ) {

        return(as.numeric(stats::filter(x = values, filter = rep(1 / window, window), sides = 2)))

}


fn_running_median <-  function( values,
                                window ) {

        ## (the first value has no power: carried from the second)
        values[1] <-  values[2]
        return(as.numeric(stats::runmed(x = values, k = window, endrule = "median")))

}


fn_run_trace <-  function( temperatures,
                           run_name ) {

        run_rows <-  temperatures[temperatures$phase == run_name, , drop = FALSE]
        stopifnot(nrow(run_rows) > 100)
        per_second_power <-  c(NA, diff(run_rows$package_energy_joules) / diff(run_rows$unix_time))
        return(data.frame( run           = run_name,
                           minutes       = (run_rows$unix_time - min(run_rows$unix_time)) / 60,
                           Tctl          = fn_moving_mean(run_rows$Tctl, moving_mean_seconds),
                           package_power = fn_running_median(per_second_power, running_median_seconds),
                           stringsAsFactors = FALSE))

}
##
## ---- Read the temperature/power logs -------------------------------------------------------------------------
##
{
      temperatures_long <-  utils::read.csv(file = file.path(overnight_dir, "temperature_long", "temperatures.csv"),
                                            stringsAsFactors = FALSE)
      temperatures_continuous <-  utils::read.csv(file = file.path(overnight_dir, "temperature_continuous",
                                                                   "temperatures.csv"),
                                                  stringsAsFactors = FALSE)
      traces <-  do.call(what = rbind, args = lapply(X = run_names, FUN = function(run_name) {
            if (grepl(pattern = "_continuous$", x = run_name)) {
                  return(fn_run_trace(temperatures = temperatures_continuous, run_name = run_name))
            }
            return(fn_run_trace(temperatures = temperatures_long, run_name = run_name))
      }))
      traces$run <-  factor(x = traces$run, levels = run_names)
      ##
      idle_rows <-  temperatures_long[temperatures_long$phase == "idle_baseline", , drop = FALSE]
      idle_Tctl <-  mean(idle_rows$Tctl)
      idle_power_watts <-  diff(range(idle_rows$package_energy_joules)) / diff(range(idle_rows$unix_time))
      message(paste0("\033[36m", "Idle: Tctl ", round(idle_Tctl, 1), " C, package power ",
                     round(idle_power_watts, 0), " W", "\033[0m"))
}
##
## ---- Maths-function variants: chain-iterations per second, relative time per iteration --------------------
##
{
      times <-  rbind(utils::read.csv(file = file.path(overnight_dir, "stan_math_test", "times.csv"),
                                      stringsAsFactors = FALSE),
                      utils::read.csv(file = file.path(overnight_dir, "stan_math_test", "times_180.csv"),
                                      stringsAsFactors = FALSE))
      times$variant <-  sub(pattern = "_.*$", replacement = "", x = times$label)
      ## Two-run difference: the long run minus the short run removes the fixed (start-up) cost of each call
      times$chain_iterations_per_second <-  times$n_chains * (times$n_iter_long - times$n_iter_short) /
                                            (times$elapsed_seconds_long_run - times$elapsed_seconds_short_run)
      variant_rates <-  stats::aggregate(chain_iterations_per_second ~ variant + n_chains + num_chunks,
                                         data = times, FUN = mean)
      reference_rates <-  variant_rates[variant_rates$variant == "AVX512", c("n_chains", "chain_iterations_per_second")]
      names(reference_rates)[2] <-  "reference_chain_iterations_per_second"
      variant_rates <-  merge(x = variant_rates, y = reference_rates, by = "n_chains")
      variant_rates$relative_time_per_iteration <-  variant_rates$reference_chain_iterations_per_second /
                                                    variant_rates$chain_iterations_per_second
      variant_rates$variant <-  factor(x = variant_rates$variant, levels = variant_names)
      variant_rates$n_chains_label <-  factor(x = variant_rates$n_chains, levels = sort(unique(variant_rates$n_chains)))
      variant_rates <-  variant_rates[order(variant_rates$n_chains, variant_rates$variant), ]
      print(variant_rates, row.names = FALSE)
      ##
      n_chains_levels <-  levels(variant_rates$n_chains_label)
      chunks_per_level <-  vapply(X = n_chains_levels, FUN = function(level) {
            unique(variant_rates$num_chunks[variant_rates$n_chains_label == level])
      }, FUN.VALUE = numeric(1))
      n_chains_axis_labels <-  lapply(X = seq_along(n_chains_levels), FUN = function(index) {
            bquote(atop(N[chains] == .(as.numeric(n_chains_levels[index])),
                        N[chunks] == .(as.numeric(chunks_per_level[index]))))
      })
}
##
## ---- Panels ------------------------------------------------------------------------------------------------
##
{
      run_colour_scale <-  ggplot2::scale_colour_manual(name = NULL, values = run_colours, breaks = run_names,
                                                        labels = run_legend_labels)
      run_linetype_scale <-  ggplot2::scale_linetype_manual(name = NULL, values = run_linetypes, breaks = run_names,
                                                            labels = run_legend_labels)
      ## time_axis <-  ggplot2::scale_x_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10.6))
      ## 4 Oct 2026: the upper limit 10.6 min silently dropped the last seconds of the two longest runs
      ## (MD_BayesMVP, N_chunks = 1: 10.69 min; AD_Stan: 10.62 min); the limit now covers every logged second:
      time_axis <-  ggplot2::scale_x_continuous(breaks = seq(0, 10, by = 2), limits = c(0, 10.75))
      ##
      plot_temperature <-  ggplot2::ggplot(data = traces,
                                           mapping = ggplot2::aes(x = minutes, y = Tctl,
                                                                  colour = run, linetype = run)) +
                           ggplot2::geom_hline(yintercept = idle_Tctl, colour = "grey45", linetype = "dotted") +
                           ggplot2::annotate(geom = "text", x = 10.6, y = idle_Tctl - 1.2, hjust = 1, size = 3,
                                             colour = "grey35",
                                             label = paste0("idle (", formatC(idle_Tctl, format = "f", digits = 1),
                                                            "°C)")) +
                           ggplot2::geom_line(linewidth = 0.55, na.rm = TRUE) +
                           run_colour_scale + run_linetype_scale + time_axis +
                           ggplot2::labs(x = NULL, y = "Tctl (°C)",
                                         title = "CPU control temperature (Tctl)") +
                           ggplot2::theme_bw() +
                           ggplot2::theme(legend.position = "bottom", legend.key.width = grid::unit(1.2, "cm"),
                                          plot.title = ggplot2::element_text(size = 11))
      ##
      plot_power <-  ggplot2::ggplot(data = traces,
                                     mapping = ggplot2::aes(x = minutes, y = package_power,
                                                            colour = run, linetype = run)) +
                     ggplot2::geom_hline(yintercept = idle_power_watts, colour = "grey45", linetype = "dotted") +
                     ggplot2::annotate(geom = "text", x = 5, y = idle_power_watts + 9, hjust = 0.5, size = 3,
                                       colour = "grey35", label = paste0("idle (", round(idle_power_watts, 0), " W)")) +
                     ggplot2::geom_line(linewidth = 0.55, na.rm = TRUE) +
                     run_colour_scale + run_linetype_scale + time_axis +
                     ggplot2::labs(x = "Minutes since the start of the run", y = "Package power (W)",
                                   title = "CPU package power") +
                     ggplot2::theme_bw() +
                     ggplot2::theme(legend.position = "bottom", legend.key.width = grid::unit(1.2, "cm"),
                                    plot.title = ggplot2::element_text(size = 11))
      ##
      plot_variants <-  ggplot2::ggplot(data = variant_rates,
                                        mapping = ggplot2::aes(x = n_chains_label, y = relative_time_per_iteration,
                                                               fill = variant)) +
                        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.75) +
                        ggplot2::geom_text(mapping = ggplot2::aes(label = paste0(formatC(relative_time_per_iteration,
                                                                                         format = "f", digits = 2),
                                                                                 "×")),
                                           position = ggplot2::position_dodge(width = 0.8), vjust = -0.35,
                                           size = 3) +
                        ggplot2::scale_fill_manual(name = NULL, values = variant_colours, breaks = variant_names,
                                                   labels = variant_legend_labels) +
                        ggplot2::scale_x_discrete(labels = n_chains_axis_labels) +
                        ggplot2::scale_y_continuous(limits = c(0, 2.4), breaks = seq(0, 2, by = 0.5),
                                                    expand = ggplot2::expansion(mult = c(0, 0.02))) +
                        ggplot2::labs(x = NULL, y = "Time per iteration\n(relative to AVX-512)",
                                      title = "Maths functions in MD_BayesMVP_chunking") +
                        ggplot2::theme_bw() +
                        ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(size = 11),
                                       panel.grid.major.x = ggplot2::element_blank())
}
##
## ---- Save ---------------------------------------------------------------------------------------------------
##
{
      ## (wrap_plots: the patchwork operators need the package attached)
      combined_plot <-  patchwork::wrap_plots(plot_temperature, plot_power, plot_variants, ncol = 1,
                                              heights = c(1, 1, 0.75), guides = "collect") +
                        patchwork::plot_annotation(tag_levels = "A") &
                        ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")
      ggplot2::ggsave(filename = figure_file, plot = combined_plot, width = 9.5, height = 10.5, dpi = 200)
      message(paste0("\033[32m", "Saved ", figure_file, "\033[0m"))
}
























