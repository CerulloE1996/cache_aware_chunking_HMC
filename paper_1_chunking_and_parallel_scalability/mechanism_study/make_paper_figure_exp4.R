##
## ---- Paper 1, experiment 4: DRAM traffic measured at the memory controllers (HPC, N = 50,000), for the manuscript.
##
{
      mechanism_dir <-  path.expand("~/Documents/Work/PhD_work/Alg_paper_analysis/paper_1_chunking_and_parallel_scalability/mechanism_study")
      bandwidth <-  utils::read.csv(file = file.path(mechanism_dir, "mechanism_DRAM_bandwidth_HPC.csv"), stringsAsFactors = FALSE)
      ##
      bandwidth$implementation <-  ifelse(bandwidth$algorithm == "MD_BayesMVP", "NicoStan+BayesMVP", "Stan model (NicoStan)")
      bandwidth$configuration <-  ifelse( bandwidth$algorithm == "MD_BayesMVP",
                                          paste0("NicoStan+BayesMVP, ", bandwidth$num_chunks, " chunk", ifelse(bandwidth$num_chunks == 1, "", "s")),
                                          ifelse( bandwidth$algorithm == "AD_Stan",
                                                  "Stan model (NicoStan), 1 chunk",
                                                  paste0("Stan model (NicoStan) + tape chunking, ", bandwidth$num_chunks, " chunks")))
      configuration_order <-  c( "NicoStan+BayesMVP, 1 chunk", "NicoStan+BayesMVP, 25 chunks", "NicoStan+BayesMVP, 100 chunks", "NicoStan+BayesMVP, 500 chunks",
                                 "Stan model (NicoStan), 1 chunk", "Stan model (NicoStan) + tape chunking, 10 chunks",
                                 "Stan model (NicoStan) + tape chunking, 250 chunks")
      bandwidth$configuration <-  factor(x = bandwidth$configuration, levels = configuration_order)
      ##
      plot_bandwidth <-  ggplot2::ggplot( data = bandwidth,
                                          mapping = ggplot2::aes(x = n_chains, y = DRAM_total_GB_per_second,
                                                                 colour = configuration, linetype = implementation, group = configuration)) +
                         ggplot2::geom_hline(yintercept = 460.8, linetype = "dashed", colour = "grey40") +
                         ggplot2::annotate(geom = "text", x = 8, y = 445, hjust = 0, size = 3.2,
                                           label = "Theoretical peak: 460.8 GB/s (12 channels, DDR5-4800)") +
                         ggplot2::geom_line() +
                         ggplot2::geom_point() +
                         ggplot2::scale_x_log10(breaks = c(8, 48, 96, 180)) +
                         ggplot2::scale_y_continuous(limits = c(0, 480), breaks = seq(0, 450, by = 50)) +
                         ggplot2::scale_linetype_manual(values = c("NicoStan+BayesMVP" = "solid", "Stan model (NicoStan)" = "longdash")) +
                         ## Line-type key (solid = NicoStan+BayesMVP, long-dashed = Stan model) above the colour key; the colour key lines use the same line types:
                         ggplot2::guides( linetype = ggplot2::guide_legend(order = 1, override.aes = list(shape = NA)),
                                          colour = ggplot2::guide_legend(order = 2, override.aes = list(linetype = c(rep("solid", 4), rep("longdash", 3))))) +
                         ggplot2::labs( x = expression(paste(N[chains], " (one thread per chain; log scale)")),
                                        y = "DRAM traffic, reads + writes (GB/s)",
                                        linetype = "Line type",
                                        colour = "Configuration") +
                         ggplot2::theme_bw() +
                         ggplot2::theme(legend.position = "right", legend.key.width = grid::unit(1.5, "cm"))
      ##
      ggplot2::ggsave( filename = file.path(mechanism_dir, "Figure_exp4_DRAM_bandwidth_HPC.png"),
                       plot = plot_bandwidth, width = 9.5, height = 5, dpi = 200)
}
























