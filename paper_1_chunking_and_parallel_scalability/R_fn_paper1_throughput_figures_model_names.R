

##
## ---- Paper 1, E1 Parts I-IV: within-sampler throughput figures, re-drawn from the saved extended-chunk-grid presentation views
##      (no refitting), with the configurations named by their model names (as in the text), the local-HPC title, an
##      N_threads (log2 scale) axis with ticks at every measured thread count, and WCP-only drawn dashed so that it stays
##      visible where the chunk search selects it (i.e. where it coincides with chunking + WCP).
##
##        Part I:   MD_BayesMVP, MD_BayesMVP_chunking                                      (Figure_ps2_BayesMVP_chunking_...)
##        Part II:  MD_BayesMVP, MD_BayesMVP_chunking, MD_BayesMVP_WCP, MD_BayesMVP_WCP_chunking (Figure_ps2_BayesMVP_variants_...)
##        Part III: AD_Stan, AD_Stan_chunked, AD_Stan_tape_chunked                          (Figure_ps2_Stan_implementations_...)
##        Part IV:  AD_Stan, AD_Stan_tape_chunked, AD_Stan_WCP, AD_Stan_WCP_chunking         (Figure_ps2_Stan_configurations_...)
##
##      The arm colours follow the earlier figures (hue_pal(4) over the four configurations; container chunking in blue).
##
##
paper_1_dir <-  path.expand("~/Documents/Work/PhD_work/Alg_paper_analysis/paper_1_chunking_and_parallel_scalability")
# views_path  <-  file.path(paper_1_dir, "paper_1_computational_outputs", "manuscript_outputs_final_both_devices_extended_chunk_grid",
#                           "data", "presentation_views.rds")
## 2026-10-03: the export with the narrow-WCP local-HPC runs (up to 90 chains x 2 threads):
views_path  <-  file.path(paper_1_dir, "paper_1_computational_outputs",
                          "manuscript_outputs_final_both_devices_narrow_WCP_2026_10_03", "data", "presentation_views.rds")
# output_dir  <-  path.expand("~/Documents/Work/PhD_work/Alg_papers_LaTeX/paper_1_v42_2026_10_03/figures_v42")
output_dir  <-  path.expand("~/Documents/Work/PhD_work/Alg_papers_LaTeX/paper_1_v42_2026_10_03/figures_v45")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
##
message(paste0("\033[36m", "Reading presentation views: ", views_path, "\033[0m"))
views <-  readRDS(file = views_path)
##
##
fn_plot_throughput_by_model_name <-  function( df,
                                               output_file,
                                               model_order,
                                               colour_values,
                                               legend_nrow = 2
) {

        df$model_name <-  factor(df$model_name, levels = model_order)
        ##
        ## WCP-only dashed, every other configuration solid:
        linetype_values <-  stats::setNames(ifelse(grepl("WCP-only", model_order), "22", "solid"), model_order)
        ##
        plot_list <-  list()
        ##
        for (dev in c("HPC", "Laptop")) {

                df_dev <-  df[df$device == dev, , drop = FALSE]
                if (nrow(df_dev) == 0) next
                ##
                ## 2026-10-03: the 176-thread WCP points were the stand-ins for the missing 180-thread WCP cells
                ## (the WCP grid stopped at 16 chains); now that 180-thread WCP allocations exist, a 176 point is
                ## dropped wherever the same configuration also has a 180 point (the data are unchanged):
                has_180 <-  paste(df_dev$N_num, df_dev$model_name)[df_dev$n_threads == 180]
                drop_176 <-  df_dev$n_threads == 176 & paste(df_dev$N_num, df_dev$model_name) %in% has_180
                df_dev <-  df_dev[!drop_176, , drop = FALSE]
                ##
                ## every measured thread count gets a tick; on the local-HPC the 176-thread (WCP) and 180-thread points share one
                ## combined tick, as they are too close to label separately on a log scale:
                if (dev == "HPC") {
                    x_breaks <-  c(1, 2, 4, 8, 16, 32, 64, 96, 128, 178)
                    x_labels <-  c("1", "2", "4", "8", "16", "32", "64", "96", "128", "176/180")
                    ## 2026-10-03: with no 176-thread points left, the last tick is simply 180:
                    if (!any(df_dev$n_threads == 176)) {
                          x_breaks[length(x_breaks)] <-  180
                          x_labels[length(x_labels)] <-  "180"
                    }
                } else {
                    x_breaks <-  c(1, 2, 4, 8, 16)
                    x_labels <-  c("1", "2", "4", "8", "16")
                }
                ##
                best_points <-  df_dev[order(df_dev$N_num, df_dev$model_name, -df_dev$total_iter_per_sec, df_dev$n_threads), , drop = FALSE]
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
                    ggplot2::geom_point(data = best_points, shape = 21, fill = "white", size = 7, stroke = 1.2, show.legend = FALSE) +
                    ggplot2::theme_bw(base_size = 28) +
                    ggplot2::theme( legend.position = ifelse(dev == "Laptop", "bottom", "none"),
                                    legend.text     = ggplot2::element_text(family = "mono"),
                                    legend.key.width = ggplot2::unit(3, "lines"),
                                    axis.text.x     = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                    ggplot2::scale_colour_manual(values = colour_values[model_order], limits = model_order, drop = FALSE) +
                    ggplot2::scale_linetype_manual(values = linetype_values, limits = model_order, drop = FALSE) +
                    ggplot2::guides( colour   = ggplot2::guide_legend(title = NULL, nrow = legend_nrow),
                                     linetype = ggplot2::guide_legend(title = NULL, nrow = legend_nrow)) +
                    ## ggplot2::ylab(expression(Total~iterations / "second"~(N[chains] %*% N[iter] / time))) +
                    ## (the multiplication sign as a character: plotmath's %*% is drawn as a centred dot by this device)
                    ggplot2::ylab(expression("Total iterations/second (" * N[chains] * " \u00D7 " * N[iter] * "/time)")) +
                    ## ggplot2::xlab(expression(log[2](N[threads]~total))) +
                    ggplot2::xlab(expression(N[threads]~"(log"[2]~"scale)")) +
                    ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels, trans = "log2") +
                    ggplot2::facet_wrap(~ N_label, scales = "free") +
                    ## ggplot2::ggtitle(ifelse(dev == "HPC", "Local HPC", "Laptop"))
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
four_colours <-  scales::hue_pal()(4)
##
## ---- Stan (Parts III and IV):
Stan_names <-  c( AD_Stan              = "AD_Stan",
                  AD_Stan_chunked      = "AD_Stan_chunked",
                  AD_Stan_tape_chunked = "AD_Stan_tape_chunked",
                  WCP_only             = "AD_Stan_WCP (WCP-only)",
                  WCP_chunking         = "AD_Stan_WCP_chunking")
Stan_colours <-  stats::setNames( c(four_colours[1], scales::hue_pal()(3)[3], four_colours[2], four_colours[4], four_colours[3]),
                                  unname(Stan_names))
##
Stan_implementations <-  views$stan[views$stan$algorithm %in% c("AD_Stan", "AD_Stan_chunked", "AD_Stan_tape_chunked"), , drop = FALSE]
Stan_implementations$model_name <-  unname(Stan_names[as.character(Stan_implementations$algorithm)])
fn_plot_throughput_by_model_name( df = Stan_implementations,
                                  output_file = file.path(output_dir, "Figure_ps2_Stan_implementations_total_throughput.png"),
                                  model_order = unname(Stan_names[c("AD_Stan", "AD_Stan_chunked", "AD_Stan_tape_chunked")]),
                                  colour_values = Stan_colours, legend_nrow = 1)
##
Stan_configurations <-  views$NicoStan_by_budget
Stan_mode_to_name <-  c( "NicoStan"                  = Stan_names[["AD_Stan"]],
                         "NicoStan-chunking"         = Stan_names[["AD_Stan_tape_chunked"]],
                         "NicoStan-WCP"              = Stan_names[["WCP_only"]],
                         "NicoStan-chunking_and_WCP" = Stan_names[["WCP_chunking"]])
Stan_configurations$model_name <-  unname(Stan_mode_to_name[as.character(Stan_configurations$comparison_mode)])
fn_plot_throughput_by_model_name( df = Stan_configurations,
                                  output_file = file.path(output_dir, "Figure_ps2_Stan_configurations_total_throughput.png"),
                                  model_order = unname(Stan_names[c("AD_Stan", "AD_Stan_tape_chunked", "WCP_only", "WCP_chunking")]),
                                  colour_values = Stan_colours, legend_nrow = 2)
##
## ---- BayesMVP (Parts I and II):
BayesMVP_mode_to_name <-  c( "BayesMVP"                  = "MD_BayesMVP",
                             "BayesMVP-chunking"         = "MD_BayesMVP_chunking",
                             "BayesMVP-WCP"              = "MD_BayesMVP_WCP (WCP-only)",
                             "BayesMVP-chunking_and_WCP" = "MD_BayesMVP_WCP_chunking")
BayesMVP_colours <-  stats::setNames(four_colours, c("MD_BayesMVP", "MD_BayesMVP_chunking", "MD_BayesMVP_WCP_chunking",
                                                     "MD_BayesMVP_WCP (WCP-only)"))
BayesMVP_all <-  views$BayesMVP_by_budget
BayesMVP_all$model_name <-  unname(BayesMVP_mode_to_name[as.character(BayesMVP_all$comparison_mode)])
fn_plot_throughput_by_model_name( df = BayesMVP_all[BayesMVP_all$model_name %in% c("MD_BayesMVP", "MD_BayesMVP_chunking"), , drop = FALSE],
                                  output_file = file.path(output_dir, "Figure_ps2_BayesMVP_chunking_total_throughput.png"),
                                  model_order = c("MD_BayesMVP", "MD_BayesMVP_chunking"),
                                  colour_values = BayesMVP_colours, legend_nrow = 1)
fn_plot_throughput_by_model_name( df = BayesMVP_all,
                                  output_file = file.path(output_dir, "Figure_ps2_BayesMVP_variants_total_throughput.png"),
                                  model_order = c("MD_BayesMVP", "MD_BayesMVP_chunking", "MD_BayesMVP_WCP (WCP-only)",
                                                  "MD_BayesMVP_WCP_chunking"),
                                  colour_values = BayesMVP_colours, legend_nrow = 2)
























