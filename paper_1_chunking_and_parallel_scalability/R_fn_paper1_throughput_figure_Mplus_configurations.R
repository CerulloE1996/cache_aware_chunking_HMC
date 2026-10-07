

##
## ---- Paper 1, E1 Part V: within-Mplus throughput figure (Mplus_standard vs. Mplus_WCP, at each N, on both machines,
##      under both iteration modes), drawn from the saved Mplus timing results (no refitting), in the same layout as the
##      E1 Parts I-IV throughput figures (R_fn_paper1_throughput_figures_model_names.R):
##
##        - total throughput = N_chains x N_iter / time (column total_iter_per_sec), against N_threads (log2 scale);
##        - for Mplus_WCP, the fastest allocation (N_chains x N_threads/chain) at each N_threads, i.e. the allocation
##          used in E2 (bold in table paper1_mplus_modes);
##        - colour = configuration (the colours of Mplus_standard and Mplus_WCP in the E2 scaling figure);
##          line type = iteration mode (BITERATIONS solid, FBITERATIONS dashed);
##        - circled = the highest throughput of each configuration and mode at each N.
##
##
paper_1_dir <-  path.expand("~/Documents/Work/PhD_work/Alg_paper_analysis/paper_1_chunking_and_parallel_scalability")
##
Mplus_modes_csv_path <-  file.path( paper_1_dir, "paper_1_computational_outputs",
                                    "manuscript_outputs_final_both_devices_narrow_WCP_2026_10_03",
                                    "data", "mplus_iteration_mode_comparison.csv")
##
output_dir <-  path.expand(paste0("~/Documents/Work/PhD_work/Alg_papers_LaTeX/paper_1_v47_2026_10_04/",
                                  "Mplus_E1_figure_2026_10_07"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
##
message(paste0("\033[36m", "Reading Mplus timing results: ", Mplus_modes_csv_path, "\033[0m"))
Mplus_results <-  utils::read.csv(file = Mplus_modes_csv_path, stringsAsFactors = FALSE)
##
##
## ---- the fastest allocation of each configuration at each N_threads (for Mplus_standard there is only one,
##      N_chains = N_threads):
Mplus_results <-  Mplus_results[order( Mplus_results$device,
                                       Mplus_results$mplus_iteration_mode,
                                       Mplus_results$algorithm,
                                       Mplus_results$N,
                                       Mplus_results$n_threads,
                                       -Mplus_results$total_iter_per_sec), , drop = FALSE]
fastest_allocation_keys <-  c("device", "mplus_iteration_mode", "algorithm", "N", "n_threads")
Mplus_fastest_allocation_at_each_N_threads <-  Mplus_results[!duplicated(Mplus_results[, fastest_allocation_keys]), ,
                                                             drop = FALSE]
##
Mplus_fastest_allocation_at_each_N_threads$model_name <-  paste0( Mplus_fastest_allocation_at_each_N_threads$algorithm,
                                                                  " (",
                                                                  Mplus_fastest_allocation_at_each_N_threads$mplus_iteration_mode,
                                                                  ")")
Mplus_fastest_allocation_at_each_N_threads$N_num   <-  as.numeric(Mplus_fastest_allocation_at_each_N_threads$N)
Mplus_fastest_allocation_at_each_N_threads$N_label <-  factor( paste0( "N = ",
                                                                       formatC( Mplus_fastest_allocation_at_each_N_threads$N_num,
                                                                                format = "d",
                                                                                big.mark = ",")),
                                                               levels = paste0( "N = ",
                                                                                formatC( c(500, 2500, 10000, 50000),
                                                                                         format = "d",
                                                                                         big.mark = ",")))
Mplus_fastest_allocation_at_each_N_threads$allocation <-  paste0( Mplus_fastest_allocation_at_each_N_threads$n_chains,
                                                                  " x ",
                                                                  Mplus_fastest_allocation_at_each_N_threads$threads_per_chain)
##
model_order <-  c( "Mplus_standard (BITERATIONS)",
                   "Mplus_WCP (BITERATIONS)",
                   "Mplus_standard (FBITERATIONS)",
                   "Mplus_WCP (FBITERATIONS)")
##
##
## ---- colours: those of Mplus_standard and Mplus_WCP in the E2 scaling figure (figure 23), i.e. hue_pal() over the
##      sorted names of its eight base configurations (as in R_fns_alg_paper_1_presentation_templates.R):
E2_base_model_names <-  sort(c( "AD_Stan", "AD_Stan_tape_chunked", "AD_Stan_WCP_chunking",
                                "MD_BayesMVP", "MD_BayesMVP_chunking", "MD_BayesMVP_WCP_chunking",
                                "Mplus_standard", "Mplus_WCP"))
E2_colours <-  stats::setNames(scales::hue_pal()(length(E2_base_model_names)), E2_base_model_names)
# colour_values <-  stats::setNames( unname(E2_colours[c("Mplus_standard", "Mplus_WCP", "Mplus_standard", "Mplus_WCP")]),
#                                    model_order)
## ---- (the two E2 Mplus colours are adjacent hues - purple and pink - and hard to tell apart in a two-configuration
##      figure; hence the colours of the E1 Parts I-IV throughput figures are used instead: the baseline colour,
##      hue_pal()(4)[1], for Mplus_standard, and the WCP colour, hue_pal()(4)[3], for Mplus_WCP)
four_colours  <-  scales::hue_pal()(4)
colour_values <-  stats::setNames( four_colours[c(1, 3, 1, 3)],
                                   model_order)
## (line type = iteration mode: BITERATIONS solid, FBITERATIONS dashed)
linetype_values <-  stats::setNames( ifelse(grepl("FBITERATIONS", model_order), "22", "solid"),
                                     model_order)
##
##
fn_plot_Mplus_throughput_by_configuration_and_iteration_mode <-  function( df,
                                                                          output_file,
                                                                          model_order,
                                                                          colour_values,
                                                                          linetype_values,
                                                                          legend_nrow = 2
) {

        df$model_name <-  factor(df$model_name, levels = model_order)
        ##
        plot_list <-  list()
        ##
        for (dev in c("HPC", "Laptop")) {

                df_dev <-  df[df$device == dev, , drop = FALSE]
                if (nrow(df_dev) == 0) next
                ##
                ## every measured thread count gets a tick; on the local-HPC the 176-thread (Mplus_WCP) and 180-thread
                ## points share one combined tick, as they are too close to label separately on a log scale:
                if (dev == "HPC") {
                    x_breaks <-  c(1, 2, 4, 8, 16, 32, 64, 96, 128, 178)
                    x_labels <-  c("1", "2", "4", "8", "16", "32", "64", "96", "128", "176/180")
                    if (!any(df_dev$n_threads == 176)) {
                          x_breaks[length(x_breaks)] <-  180
                          x_labels[length(x_labels)] <-  "180"
                    }
                } else {
                    x_breaks <-  c(1, 2, 4, 8, 16)
                    x_labels <-  c("1", "2", "4", "8", "16")
                }
                ##
                best_points <-  df_dev[order(df_dev$N_num, df_dev$model_name, -df_dev$total_iter_per_sec, df_dev$n_threads), ,
                                       drop = FALSE]
                best_points <-  best_points[!duplicated(best_points[, c("N_num", "model_name")]), , drop = FALSE]
                ##
                p <-  ggplot2::ggplot( df_dev,
                                       ggplot2::aes( x        = n_threads,
                                                     y        = total_iter_per_sec,
                                                     colour   = model_name,
                                                     linetype = model_name,
                                                     group    = model_name)) +
                    ggplot2::geom_line(linewidth = 2) +
                    ggplot2::geom_point(size = 5) +
                    ggplot2::geom_point( data = best_points, shape = 21, fill = "white", size = 7, stroke = 1.2,
                                         show.legend = FALSE) +
                    ggplot2::theme_bw(base_size = 28) +
                    ggplot2::theme( legend.position  = ifelse(dev == "Laptop", "bottom", "none"),
                                    legend.text      = ggplot2::element_text(family = "mono"),
                                    legend.key.width = ggplot2::unit(3, "lines"),
                                    axis.text.x      = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                    ggplot2::scale_colour_manual(values = colour_values[model_order], limits = model_order, drop = FALSE) +
                    ggplot2::scale_linetype_manual(values = linetype_values[model_order], limits = model_order, drop = FALSE) +
                    ggplot2::guides( colour   = ggplot2::guide_legend(title = NULL, nrow = legend_nrow),
                                     linetype = ggplot2::guide_legend(title = NULL, nrow = legend_nrow)) +
                    ## (the multiplication sign as a character: plotmath's %*% is drawn as a centred dot by this device)
                    ggplot2::ylab(expression("Total iterations/second (" * N[chains] * " × " * N[iter] * "/time)")) +
                    ggplot2::xlab(expression(N[threads]~"(log"[2]~"scale)")) +
                    ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels, trans = "log2") +
                    ggplot2::facet_wrap(~ N_label, scales = "free") +
                    ggplot2::ggtitle(ifelse(dev == "HPC", "local-HPC", "Laptop"))
                ##
                plot_list[[dev]] <-  p

        }
        ##
        combined_plot <-  patchwork::wrap_plots(plot_list[["HPC"]], plot_list[["Laptop"]], ncol = 1)
        ggplot2::ggsave(output_file, combined_plot, width = 16, height = 24, dpi = 100)
        message(paste0("\033[32m", "Saved: ", output_file, "\033[0m"))
        invisible(combined_plot)

}
##
##
fn_plot_Mplus_throughput_by_configuration_and_iteration_mode(
        df              = Mplus_fastest_allocation_at_each_N_threads,
        output_file     = file.path(output_dir, "Figure_ps2_Mplus_configurations_total_throughput.png"),
        model_order     = model_order,
        colour_values   = colour_values,
        linetype_values = linetype_values,
        legend_nrow     = 2)
##
##
## ---- the plotted values (one row per point), for checking the figure and the text against the data:
plotted_values_columns <-  c( "device", "N", "mplus_iteration_mode", "algorithm", "n_threads", "allocation",
                              "n_chains", "threads_per_chain", "n_iter", "seconds_per_iteration_mean",
                              "total_iter_per_sec")
plotted_values <-  Mplus_fastest_allocation_at_each_N_threads[, plotted_values_columns, drop = FALSE]
plotted_values_path <-  file.path(output_dir, "Figure_ps2_Mplus_configurations_total_throughput_values.csv")
utils::write.csv(plotted_values, file = plotted_values_path, row.names = FALSE)
message(paste0("\033[32m", "Saved: ", plotted_values_path, "\033[0m"))
























