##
## =====================================================================================================================================
## R_fns_alg_paper_1_presentation_templates.R
##
## Presentation functions extracted from the existing PS1/PS2 analysis code; no benchmark or legacy file-reading code is executed.
## Called privately by R_fns_alg_paper_1_figures_tables.R. The originals remain available in the legacy files.
##
fn_paper1_presentation_templates <-  function() {

        require(ggplot2)
        ##
        require(dplyr)
        ##
        require(tidyr)
        ##
        require(patchwork)
        ##
        ## ---- Number labels: thousands separators only from 10,000 upwards (e.g. 2500, 10,000, 50,000):
        ##
        fn_paper1_format_number_commas_from_10000 <-  function( numbers ) {

                vapply(X = numbers,
                       FUN = function(one_number) {
                             if (is.finite(one_number) && abs(one_number) >= 10000) format(x = one_number, big.mark = ",", scientific = FALSE, trim = TRUE)
                             else format(x = one_number, scientific = FALSE, trim = TRUE)
                       },
                       FUN.VALUE = character(1),
                       USE.NAMES = FALSE)

        }
        ##
        t_phys_HPC <-  96
        ##
        t_smt_HPC <-  180
        ##
        t_phys_Laptop <-  8
        ##
        t_smt_Laptop <-  16

        ## Missing coverage is shown explicitly without manufacturing observations.
        empty_panel <-  function(label) {

                ggplot() + annotate("text", x = 0, y = 0, label = label, size = 8) +
                    theme_void()

        }

        ## Use the complete input's named palette in every device panel, including absent levels.
        shared_colour_scale <-  function(values) {

                labels <-  sort(unique(as.character(values[!is.na(values)])))
                ##
                colours <-  if (length(labels)) scales::hue_pal()(length(labels)) else character()
                ##
                scale_colour_manual(values = setNames(colours, labels), limits = labels, drop = FALSE)

        }

        ## ---- From
        ## legacy/ps_1_optimizing_N_chunks_and_N_threads/ps_1_optimizing_N_chunks.R:369
        R_fn_plot_ps1_N_chunks_ggplot_1 <-  function( df_both,
                                                      save_plot = TRUE,
                                                      output_path,
                                                      n_threads_for_HPC = 96,
                                                      n_threads_for_Laptop = 8
        ) {

                ##
                make_half <-  function( N_range,
                                        hide_legend) {

                        df_f <-  df_both %>%
                            filter( as.numeric(as.character(N)) %in% N_range,
                                    n_threads == ifelse(device == "HPC", n_threads_for_HPC, n_threads_for_Laptop)) %>%
                            dplyr::mutate(Device = ifelse( device == "HPC",
                                                    paste0("HPC (",    n_threads_for_HPC,    " threads)"),
                                                    paste0("Laptop (", n_threads_for_Laptop, " threads)")))
                        ##
                        if (nrow(df_f) == 0) {

                            return(empty_panel(paste0( "No matching results supplied for N = ",
                                                       paste(N_range, collapse = ", "))))

                        }
                        ##
                        p <-  ggplot(df_f, aes(x = N_chunks, y = time_avg, group = Device)) +
                            geom_point(size = 5) +
                            geom_errorbar( linewidth = 1,
                                           width = 0.02,
                                           aes( ymin = time_avg - time_SD,        ## +/- 1 SD (caption says "standard deviation")
                                                ymax = time_avg + time_SD,
                                                linetype = Device)) +
                            geom_line(linewidth = 1, aes(linetype = Device)) +
                            theme_bw(base_size = 24) +
                            ylab("Time (sec.)") + xlab(expression(N[chunks])) +
                            facet_wrap(~ N_label, scales = "free") +
                            theme( axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                                   legend.position = if (hide_legend) "none" else "bottom")
                        ##
                        p

                }
                ##
                has_bottom <-  any(as.numeric(as.character(df_both$N)) %in% c(10000, 50000) &
                                   df_both$n_threads == ifelse( df_both$device == "HPC",
                                                                n_threads_for_HPC,
                                                                n_threads_for_Laptop))
                ##
                p_top    <-  make_half(c(500, 2500),     hide_legend = has_bottom) + xlab(" ")
                ##
                p_bottom <-  make_half(c(10000, 50000),  hide_legend = FALSE)
                ##
                combined <-  p_top + p_bottom + plot_layout(ncol = 1)
                ##
                print(combined)
                ##
                if (save_plot) {

                    ggsave( file.path( output_path,
                                       paste0( "Figure_N_chunks_pilot_study_plot_1_",
                                               "n_threads_HPC_",
                                               n_threads_for_HPC,
                                               "_",
                                               "Laptop_",
                                               n_threads_for_Laptop,
                                               ".png")),
                            combined,
                            width = 16,
                            height = 16,
                            dpi = 100)

                }
                ##
                invisible(combined)

        }

        ## ---- From
        ## legacy/ps_1_optimizing_N_chunks_and_N_threads/ps_1_optimizing_N_chunks.R:432
        R_fn_plot_ps1_efficiency <-  function( df_HPC,
                                               df_Laptop,
                                               threads_HPC    = c(64, t_phys_HPC, t_smt_HPC),
                                               threads_Laptop = c(4,  t_phys_Laptop, t_smt_Laptop),
                                               save_plot = TRUE,
                                               output_path
        ) {

                make_dev_plot <-  function( df,
                                            threads_keep,
                                            device_name) {

                        df <-  df %>% filter(n_threads %in% threads_keep)
                        ##
                        if (nrow(df) == 0) {

                            return(empty_panel(paste0("No matching ", device_name, " results supplied")))

                        }
                        ##
                        ggplot( df,
                                aes(x = N_chunks, y = Efficiency, colour = N_threads, group = N_threads)) +
                            geom_point(size = 5) +
                            geom_line(linewidth = 1) +
                            theme_bw(base_size = 28) +
                            theme( legend.position = "bottom",
                                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                            ylab(expression(Efficiency~(N[threads] / time))) +
                            xlab(expression(N[chunks])) +
                            facet_wrap(~ N_label, scales = "free")

                }
                ##
                combined <-  make_dev_plot(df_HPC, threads_HPC, "HPC") +
                    make_dev_plot(df_Laptop, threads_Laptop, "Laptop") +
                    plot_layout(ncol = 1)
                ##
                print(combined)
                ##
                if (save_plot) {

                    ggsave( file.path(output_path, "Figure_N_chunks_pilot_study_plot_3.png"),
                            combined,
                            width = 16,
                            height = 24,
                            dpi = 100)

                }
                ##
                invisible(combined)

        }

        ## ---- Paper 1 WCP chunk-search views -------------------------------------------------------------------------------------------
        fn_plot_paper1_WCP_chunk_search <-  function( chunk_search,
                                                       best_chunks,
                                                       output_path,
                                                       file_prefix
        ) {

                required_columns <-  c("device", "algorithm", "N", "configuration_id", "n_chains", "threads_per_chain",
                                       "num_chunks", "chain_rate", "relative_chunk_throughput", "Algorithm_label",
                                       "selected_best_chunks")
                missing_columns <-  setdiff(x = required_columns, y = names(x = chunk_search))
                if (length(x = missing_columns)) stop("WCP chunk-search data are missing: ", paste(missing_columns, collapse = ", "))
                ##
                if (length(x = unique(x = chunk_search$device)) != 1L || length(x = unique(x = chunk_search$algorithm)) != 1L ||
                    length(x = unique(x = chunk_search$N)) != 1L) stop("Supply one device, algorithm and N for each WCP chunk-search plot.")
                ##
                if (!is.logical(x = chunk_search$selected_best_chunks) || anyNA(x = chunk_search$selected_best_chunks) ||
                    !setequal(x = best_chunks$configuration_id, y = chunk_search$configuration_id[chunk_search$selected_best_chunks])) {

                    stop("Selected WCP rows disagree with the conditional chunk optima.")

                }
                ##
                chunk_search$WCP_label <-  factor(x = chunk_search$threads_per_chain,
                                                  levels = sort(x = unique(x = chunk_search$threads_per_chain)))
                chunk_search$chunk_label <-  factor(x = chunk_search$num_chunks, levels = sort(x = unique(x = chunk_search$num_chunks)))
                ##
                plot_title <-  paste0(chunk_search$device[1], ": ", chunk_search$Algorithm_label[1],
                                       ", N = ", fn_paper1_format_number_commas_from_10000(chunk_search$N[1]))
                ##
                ## Keep all measured chunks; the outline marks each fixed chain/WCP allocation's optimum.
                wcp_axis_label <-  if ("execution_backend" %in% names(x = chunk_search) &&
                                        any(chunk_search$execution_backend == "NicoStan_BridgeStan")) {

                    "WCP budget / chain (shared pool)"

                } else "WCP threads / chain"
                ##
                chunk_search_plot <-  ggplot( data = chunk_search,
                                               mapping = aes(x = chunk_label, y = chain_rate, colour = WCP_label, group = WCP_label)) +
                    geom_line(linewidth = 0.8) + geom_point(size = 3) +
                    geom_point(data = chunk_search[chunk_search$selected_best_chunks, , drop = FALSE],
                               shape = 21, fill = "white", size = 5, stroke = 1.2) +
                    theme_bw(base_size = 20) +
                    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
                    labs( x = "Chunks", y = "Within-method efficiency (chains / time)", colour = wcp_axis_label,
                           title = plot_title, subtitle = "All measured chunks; outlined points maximise throughput at each fixed chain/WCP count") +
                    facet_wrap(facets = ~ n_chains, scales = "free_y", labeller = label_both)
                ##
                ## Keep every WCP choice after selecting chunks and label the selected chunk count.
                best_chunks_plot <-  ggplot(data = best_chunks, mapping = aes(x = threads_per_chain, y = chain_rate)) +
                    geom_line(mapping = aes(group = n_chains), linewidth = 0.8) + geom_point(size = 4) +
                    geom_text(mapping = aes(label = num_chunks), vjust = -0.8, size = 4.5) +
                    scale_x_continuous(breaks = sort(x = unique(x = best_chunks$threads_per_chain))) +
                    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
                    theme_bw(base_size = 20) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                    labs( x = wcp_axis_label, y = "Within-method efficiency (chains / time)", title = plot_title,
                           subtitle = "After conditional chunk optimisation; point labels give the selected number of chunks") +
                    facet_wrap(facets = ~ n_chains, scales = "free_y", labeller = label_both)
                ##
                plot_paths <-  file.path(output_path, paste0(file_prefix, c("_chunk_search.png", "_optimal_chunks_by_WCP.png")))
                plot_height <-  4 + 4 * ceiling(x = length(x = unique(x = chunk_search$n_chains)) / 3)
                ggsave(filename = plot_paths[1], plot = chunk_search_plot, width = 16, height = plot_height, dpi = 150)
                ggsave(filename = plot_paths[2], plot = best_chunks_plot, width = 16, height = plot_height, dpi = 150)
                ##
                return(invisible(x = plot_paths))

        }

        ## ---- Paper 1 configuration longtable -----------------------------------------------------------------------------------------
        fn_make_paper1_configuration_table <-  function( configurations,
                                                          caption,
                                                          label,
                                                          output_file,
                                                          include_algorithm = FALSE
        ) {

                fn_escape_tex <-  function(text) {

                        replacements <-  c("\\" = "\\textbackslash{}", "#" = "\\#", "$" = "\\$", "%" = "\\%", "&" = "\\&",
                                           "_" = "\\_", "{" = "\\{", "}" = "\\}", "~" = "\\textasciitilde{}", "^" = "\\textasciicircum{}")
                        text <-  ifelse(test = is.na(x = text), yes = "-", no = as.character(x = text))
                        ## Replace original characters once so inserted LaTeX commands are not escaped again.
                        return(vapply( X = strsplit(x = text, split = "", fixed = TRUE),
                                       FUN = function(characters) {

                                               replace <-  characters %in% names(x = replacements)
                                               characters[replace] <-  replacements[characters[replace]]
                                               return(paste(characters, collapse = ""))

                                       },
                                       FUN.VALUE = character(length = 1)))

                }
                ##
                fn_format_column <-  function(column, digits = 4) {

                        if (!column %in% names(x = configurations)) return(rep(x = "-", times = nrow(x = configurations)))
                        column_values <-  configurations[[column]]
                        ## Mplus stores a placeholder in num_chunks; it has no independent chunk setting in this study.
                        if (column == "num_chunks") column_values[grepl(pattern = "^Mplus_", x = configurations$algorithm)] <-  NA_real_
                        ##
                        return(ifelse(test = is.na(x = column_values), yes = "-",
                                      no = formatC(x = column_values, digits = digits, format = "fg")))

                }
                ##
                rate_column <-  if (include_algorithm) "adjusted_scaling" else "chain_rate"
                rate_heading <-  if (include_algorithm) "Adjusted scaling" else "Chains / time"
                table_columns <-  lapply( X = c("n_chains", "threads_per_chain", "num_chunks", "stan_chunk_size", "n_threads",
                                                "time_mean", rate_column), FUN = fn_format_column)
                if (include_algorithm) table_columns <-  c(list(fn_escape_tex(text = configurations$Algorithm_label)), table_columns)
                ##
                header <-  paste0(if (include_algorithm) "Algorithm & " else "",
                                   "Chains & WCP / chain & Chunks & Stan chunk size & Total threads & Mean seconds & ", rate_heading, " \\\\")
                column_count <-  length(x = table_columns)
                column_spec <-  paste0(if (include_algorithm) "l" else "", "rrrrrrr")
                ##
                lines <-  c("\\begingroup", "\\small", paste0("\\begin{longtable}{", column_spec, "}"),
                            paste0("\\caption{", fn_escape_tex(text = caption), "}\\label{", label, "}\\\\"),
                            "\\toprule", header, "\\midrule", "\\endfirsthead", "\\toprule", header, "\\midrule", "\\endhead",
                            "\\midrule", paste0("\\multicolumn{", column_count, "}{r}{Continued on next page}\\\\"),
                            "\\endfoot", "\\bottomrule", "\\endlastfoot")
                ##
                rows <-  vapply( X = seq_len(length.out = nrow(x = configurations)),
                                  FUN = function(row_index) {

                                          fields <-  vapply(X = table_columns, FUN = function(column) column[row_index],
                                                            FUN.VALUE = character(length = 1))
                                          return(paste0(paste(fields, collapse = " & "), " \\\\"))

                                  },
                                  FUN.VALUE = character(length = 1))
                ##
                writeLines(text = c(lines, rows, "\\end{longtable}", "\\endgroup"), con = output_file)
                return(invisible(x = output_file))

        }

        ## ---- From
        ## legacy/ps_1_optimizing_N_chunks_and_N_threads/ps_1_optimizing_N_chunks.R:473
        get_best_chunks <-  function(df) {

                df %>%
                    group_by(device, N, n_threads) %>%
                    dplyr::arrange(N_chunks_num, .by_group = TRUE) %>%
                    slice_max(Efficiency, n = 1, with_ties = FALSE) %>%    ## first row among ties = smallest N_chunks
                    ungroup() %>%
                    select(device, N, n_threads, best_N_chunks = N_chunks_num, time_avg, Efficiency)

        }

        ## ---- From
        ## legacy/ps_1_optimizing_N_chunks_and_N_threads/ps_1_optimizing_N_chunks.R:491
        get_smt_benefit <-  function( df,
                                      t_phys,
                                      t_smt
        ) {

                wide <-  df %>%
                    filter(n_threads %in% c(t_phys, t_smt)) %>%
                    dplyr::mutate(which_t = ifelse(n_threads == t_phys, "eff_physical", "eff_smt")) %>%
                    select(N, N_chunks_num, which_t, Efficiency) %>%
                    tidyr::pivot_wider(names_from = which_t, values_from = Efficiency)
                ##
                for (column in c("eff_physical", "eff_smt")) {

                    if (!column %in% names(wide)) wide[[column]] <-  rep(NA_real_, nrow(wide))

                }
                ##
                df %>% distinct(N, N_chunks_num) %>%
                    left_join(wide, by = c("N", "N_chunks_num")) %>%
                    dplyr::mutate( smt_better = eff_smt > eff_physical,
                            gain_pct   = 100 * (eff_smt / eff_physical - 1)) %>%
                    dplyr::arrange(N, N_chunks_num)

        }

        ## ---- From
        ## legacy/ps_1_optimizing_N_chunks_and_N_threads/ps_1_optimizing_N_chunks.R:511
        get_smt_threshold <-  function( smt_benefit,
                                        device_name) {

                smt_benefit %>%
                    group_by(N) %>%
                    dplyr::summarise( N_chunks_sel = if (any(smt_better, na.rm = TRUE))
                               min(N_chunks_num[which(smt_better)]) else NA_real_,
                               .groups = "drop") %>%
                    left_join(smt_benefit, by = "N") %>%
                    filter(is.na(N_chunks_sel) | N_chunks_num == N_chunks_sel) %>%
                    distinct(N, .keep_all = TRUE) %>%
                    dplyr::mutate( device = device_name,
                            across( c(eff_physical, eff_smt, gain_pct),
                                    ~ ifelse(is.na(N_chunks_sel), NA_real_, .x))) %>%
                    select(device, N, N_chunks = N_chunks_sel, eff_physical, eff_smt, gain_pct)

        }

        ## ---- From
        ## legacy/ps_1_optimizing_N_chunks_and_N_threads/ps_1_optimizing_N_chunks.R:527
        get_smt_best <-  function( smt_benefit,
                                   device_name) {

                smt_benefit %>%
                    group_by(N) %>%
                    dplyr::arrange(desc(eff_smt), N_chunks_num, .by_group = TRUE) %>%
                    slice_head(n = 1) %>%
                    ungroup() %>%
                    dplyr::mutate( device = device_name,
                            N_chunks_num = ifelse(is.na(eff_smt), NA_real_, N_chunks_num)) %>%
                    select(device, N, N_chunks = N_chunks_num, eff_physical, eff_smt, gain_pct)

        }

        ## ---- From
        ## legacy/ps_1_optimizing_N_chunks_and_N_threads/ps_1_optimizing_N_chunks.R:549
        fmt3 <-  function(x) {

                out <-  formatC(x, digits = 3, format = "fg", flag = "#")
                ##
                out <-  trimws(sub("\\.$", "", out))
                ##
                out[is.na(x)] <-  "---"
                ##
                out

        }

        ## ---- From
        ## legacy/ps_1_optimizing_N_chunks_and_N_threads/ps_1_optimizing_N_chunks.R:550
        fmtN <-  function(N)   paste0("$", fn_paper1_format_number_commas_from_10000(as.numeric(as.character(N))), "$")

        ## ---- From
        ## legacy/ps_1_optimizing_N_chunks_and_N_threads/ps_1_optimizing_N_chunks.R:552
        emit_tex <-  function( lines,
                               file_path) {

                cat(lines, sep = "\n")
                ##
                cat("\n")
                ##
                writeLines(lines, file_path)
                ##
                message("Wrote: ", file_path)

        }

        ## ---- From
        ## legacy/ps_1_optimizing_N_chunks_and_N_threads/ps_1_optimizing_N_chunks.R:560
        make_ps1_best_chunks_table_tex <-  function( best_chunks_df,
                                                     device_name,
                                                     thread_vec,
                                                     caption,
                                                     label,
                                                     file_path
        ) {

                wide <-  best_chunks_df %>%
                    dplyr::mutate(N_num = as.numeric(as.character(N))) %>%
                    select(N_num, n_threads, best_N_chunks) %>%
                    tidyr::pivot_wider(names_from = n_threads, values_from = best_N_chunks) %>%
                    dplyr::arrange(N_num)
                ##
                lines <-  c(
                             "\\begin{table}[H]", "\\centering",
                             paste0("\\caption{", caption, "}"),
                             paste0("\\label{", label, "}"),
                             paste0("\\begin{tabular}{l", paste(rep("c", length(thread_vec)), collapse = ""), "}"),
                             "\\hline",
                             paste0("       & \\multicolumn{", length(thread_vec), "}{c}{$N_{\\text{threads}}$} \\\\"),
                             paste0("$N$          & ", paste(thread_vec, collapse = " & "), " \\\\ \\hline")
                )
                ##
                for (i in seq_len(nrow(wide))) {

                    vals <-  sapply( as.character(thread_vec),
                                     function(tt) if (tt %in% names(wide) && !is.na(wide[[tt]][i])) wide[[tt]][i] else "---")
                    ##
                    lines <-  c(lines, paste0(fmtN(wide$N_num[i]), " & ", paste(vals, collapse = " & "), " \\\\"))

                }
                ##
                lines <-  c(lines, "\\hline", "\\end{tabular}", "\\end{table}")
                ##
                emit_tex(lines, file_path)

        }

        ## ---- From
        ## legacy/ps_1_optimizing_N_chunks_and_N_threads/ps_1_optimizing_N_chunks.R:607
        make_ps1_smt_table_tex <-  function( smt_df,
                                             chunks_col_header,
                                             caption,
                                             label,
                                             file_path
        ) {

                lines <-  c(
                             "\\begin{table}[H]", "\\centering",
                             paste0("\\caption{", caption, "}"),
                             paste0("\\label{", label, "}"),
                             "\\begin{tabular}{llrccc}",
                             "\\hline",
                             paste0( "\\textbf{Device} & $N$ & \\textbf{",
                                     chunks_col_header,
                                     "} & \\textbf{Eff.\\ (physical)} & \\textbf{Eff.\\ (SMT)} & \\textbf{Gain (\\%)} \\\\ \\hline")
                )
                ##
                for (dev in c("HPC", "Laptop")) {

                    block <-  smt_df %>% filter(device == dev) %>%
                        dplyr::mutate(N_num = as.numeric(as.character(N))) %>% dplyr::arrange(N_num)
                    ##
                    for (i in seq_len(nrow(block))) {

                        if (is.na(block$N_chunks[i])) {

                            lines <-  c(lines, paste0(dev, " & ", fmtN(block$N_num[i]), " & --- & --- & --- & --- \\\\"))

                        } else {

                            lines <-  c(lines, paste0( dev,
                                                       "    & ",
                                                       fmtN(block$N_num[i]),
                                                       " & ",
                                                       block$N_chunks[i],
                                                       " & ",
                                                       fmt3(block$eff_physical[i]),
                                                       " & ",
                                                       fmt3(block$eff_smt[i]),
                                                       " & ",
                                                       if (is.na(block$gain_pct[i])) "---" else formatC( x = block$gain_pct[i],
                                                                                                         format = "f",
                                                                                                         digits = 1,
                                                                                                         decimal.mark = "."),
                                                       " \\\\"))

                        }

                    }
                    ##
                    lines <-  c(lines, "\\hline")

                }
                ##
                lines <-  c(lines, "\\end{tabular}", "\\end{table}")
                ##
                emit_tex(lines, file_path)

        }

        ## ---- From
        ## legacy/ps_2_parallel_scaling_vs_Mplus_Stan/functions/R_fns_ps2.R:533
        R_fn_plot_Stan_variants_throughput <-  function( stan_df,
                                                         save_plot = TRUE,
                                                         output_path,
                                                         file_prefix = "Figure_ps2_Stan_variants",
                                                         highlight_best = FALSE
        ) {

                colour_scale <-  shared_colour_scale(stan_df$Stan_variant_label)
                ##
                plot_list <-  list()
                ##
                for (dev in c("HPC", "Laptop")) {

                    ##
                    df_dev <-  stan_df %>% filter(device == dev)
                    ##
                    if (nrow(df_dev) == 0) next
                    ##
                    ## Laptop breaks start at 1: the grid now includes the 1-chain serial baseline on both devices.
                    x_breaks <-  if (dev == "HPC") c(1, 2, 4, 8, 16, 32, 64, 96, 180) else c(1, 2, 4, 8, 16)
                    ##
                    p <-  ggplot( df_dev,
                                  aes( x = n_threads,
                                       y = total_iter_per_sec,
                                       colour = Stan_variant_label,
                                       group  = Stan_variant_label)) +
                        geom_point(size = 5) +
                        geom_line(linewidth = 2) +
                        theme_bw(base_size = 28) +
                        theme( legend.position = ifelse(dev == "Laptop", "bottom", "none"),
                               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                        colour_scale +
                        guides(colour = guide_legend(title = NULL, nrow = 2)) +
                        ylab(expression(Total~iterations / "second"~(N[chains] %*% N[iter] / time))) +
                        xlab(expression(log[2](N[threads]~total))) +
                        scale_x_continuous(breaks = x_breaks, trans = "log2") +
                        facet_wrap(~ N_label, scales = "free") +
                        ggtitle(ifelse(dev == "HPC", "Local HPC", "Laptop"))
                    if (highlight_best) {

                        best_points <-  df_dev %>% group_by(N_num, Stan_variant_label) %>%
                            dplyr::arrange(desc(total_iter_per_sec), n_threads, N_chunks_num, .by_group = TRUE) %>%
                            slice_head(n = 1) %>% ungroup()
                        p <-  p + geom_point(data = best_points, shape = 21, fill = "white", size = 7, stroke = 1.2)

                    }
                    ##
                    plot_list[[dev]] <-  p

                }
                ##
                if (length(plot_list) == 0) {

                    message("No matching observations available for this plot; no file written.")
                    ##
                    return(invisible(NULL))

                }
                ##
                if (length(plot_list) == 2) {

                    combined_plot <-  plot_list[["HPC"]] + plot_list[["Laptop"]] + plot_layout(ncol = 1)

                } else {

                    combined_plot <-  plot_list[[1]] + theme(legend.position = "bottom")

                }
                ##
                print(combined_plot)
                ##
                if (save_plot) {

                    ggsave( file.path(output_path, paste0(file_prefix, "_total_throughput.png")),
                            combined_plot,
                            width = 16,
                            height = ifelse(length(plot_list) == 2, 24, 12),
                            dpi = 100)

                }
                ##
                return(combined_plot)

        }

        ## ---- From
        ## legacy/ps_2_parallel_scaling_vs_Mplus_Stan/functions/R_fns_ps2.R:1056
        R_fn_plot_Stan_variants_relative <-  function( stan_df,
                                                       save_plot = TRUE,
                                                       output_path,
                                                       file_prefix = "Figure_ps2_Stan_variants",
                                                       baseline_label = "Stan",
                                                       baseline_mode = NULL
        ) {

                ## baseline_mode: when given, the baseline rows are those whose Stan_variant_label equals it
                ## (e.g. "BayesMVP" for the within-BayesMVP comparison); otherwise the plain AD_Stan rows, as before.
                is_baseline_row <-  if (is.null(baseline_mode)) stan_df$algorithm == "AD_Stan" else
                    stan_df$Stan_variant_label == baseline_mode
                ##
                baseline_df <-  stan_df[is_baseline_row, , drop = FALSE] %>%
                    select(device, N, n_threads, baseline_total_iter_per_sec = total_iter_per_sec)
                ##
                rel_df <-  stan_df[!is_baseline_row, , drop = FALSE] %>%
                    left_join(baseline_df, by = c("device", "N", "n_threads")) %>%
                    dplyr::mutate(rel_throughput = total_iter_per_sec / baseline_total_iter_per_sec) %>%
                    filter(!is.na(rel_throughput))
                ##
                colour_scale <-  shared_colour_scale(rel_df$Stan_variant_label)
                ##
                plot_list <-  list()
                ##
                for (dev in c("HPC", "Laptop")) {

                    ##
                    df_dev <-  rel_df %>% filter(device == dev)
                    ##
                    if (nrow(df_dev) == 0) next
                    ##
                    ## Laptop breaks start at 1: the grid now includes the 1-chain serial baseline on both devices.
                    x_breaks <-  if (dev == "HPC") c(1, 2, 4, 8, 16, 32, 64, 96, 180) else c(1, 2, 4, 8, 16)
                    ##
                    p <-  ggplot( df_dev,
                                  aes( x = n_threads,
                                       y = rel_throughput,
                                       colour = Stan_variant_label,
                                       group  = Stan_variant_label)) +
                        geom_hline(yintercept = 1.0, linetype = "dashed", colour = "grey40", linewidth = 1) +
                        geom_point(size = 5) +
                        geom_line(linewidth = 2) +
                        theme_bw(base_size = 28) +
                        theme( legend.position = ifelse(dev == "Laptop", "bottom", "none"),
                               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                        colour_scale +
                        guides(colour = guide_legend(title = NULL, nrow = 2)) +
                        ylab(paste0("Total throughput relative to ", baseline_label, " baseline")) +
                        xlab(expression(log[2](N[threads]~total))) +
                        scale_x_continuous(breaks = x_breaks, trans = "log2") +
                        facet_wrap(~ N_label, scales = "free") +
                        ggtitle(ifelse(dev == "HPC", "Local HPC", "Laptop"))
                    ##
                    plot_list[[dev]] <-  p

                }
                ##
                if (length(plot_list) == 0) {

                    message("No matching observations available for this plot; no file written.")
                    ##
                    return(invisible(NULL))

                }
                ##
                if (length(plot_list) == 2) {

                    combined_plot <-  plot_list[["HPC"]] + plot_list[["Laptop"]] + plot_layout(ncol = 1)

                } else {

                    combined_plot <-  plot_list[[1]] + theme(legend.position = "bottom")

                }
                ##
                print(combined_plot)
                ##
                if (save_plot) {

                    ggsave( file.path(output_path, paste0(file_prefix, "_relative_to_baseline.png")),
                            combined_plot,
                            width = 16,
                            height = ifelse(length(plot_list) == 2, 24, 12),
                            dpi = 100)

                }
                ##
                return(combined_plot)

        }

        ## ---- From
        ## legacy/ps_2_parallel_scaling_vs_Mplus_Stan/functions/R_fns_ps2.R:1211
        ##
        ## metric = "serial_speedup": the speed-up S = N_chains x T0_eq / T of equation eq:paper1_serial_efficiency, with
        ## exactly the one-chain, one-thread reference of the ps2 tables (get_serial_speedup() below); ideal scaling is S = N_threads.
        ## values_file (optional): CSV of the plotted values.
        ##
        R_fn_plot_ps2_scaling <-  function( ps2_df,
                                            ## metric      = c("adjusted", "normalised"),
                                            metric      = c("adjusted", "normalised", "serial_speedup"),
                                            save_plot   = TRUE,
                                            ## output_path
                                            output_path,
                                            values_file = NULL
        ) {

                ##
                metric <-  match.arg(metric)
                ##
                if (nrow(ps2_df) == 0) {

                    message("No matching observations available for this plot; no file written.")
                    ##
                    return(invisible(NULL))

                }
                ##
                df_plot <-  ps2_df %>%
                    group_by(device, N_num, N_label, Algorithm_label) %>%
                    dplyr::mutate( work_rate = chain_rate,
                            y_val = if (metric == "normalised") {

                    work_rate / work_rate[n_threads == min(n_threads)][1]     ## all arms start at 1

                } else {

                    ## old manuscript metric; a flagged non-positive two-run time is never the minimum 
                    min(time_avg[time_avg > 0], na.rm = TRUE) * work_rate

                }) %>%
                    ungroup()
                ##
                ## ---- Speed-up over the one-chain, one-thread reference of the ps2 tables (replaces y_val above):
                ##
                if (metric == "serial_speedup") {

                    df_plot <-  get_serial_speedup(ps2_df, t_serial = 1)
                    df_plot$y_val <-  df_plot$serial_speedup_S
                    ##
                    no_reference <-  !is.finite(df_plot$y_val)
                    if (any(no_reference)) {

                        message(paste0("\033[36mR_fn_plot_ps2_scaling: ", sum(no_reference),
                                       " selected configuration(s) have no one-thread reference and are not plotted: ",
                                       paste(unique(paste0(df_plot$device[no_reference], ", N = ", df_plot$N_num[no_reference], ", ",
                                                           df_plot$Algorithm_label[no_reference])), collapse = "; "), "\033[0m"))

                    }
                    df_plot <-  df_plot[!no_reference, , drop = FALSE]
                    ##
                    if (!is.null(values_file)) {

                        value_columns <-  c("device", "N_num", "Algorithm_label", "algorithm", "n_threads", "n_chains", "threads_per_chain",
                                            "num_chunks", "n_iter", "time_avg", "serial_time_equivalent", "serial_speedup_S", "parallel_efficiency_E")
                        utils::write.csv(x = df_plot[, value_columns, drop = FALSE], file = values_file, row.names = FALSE)

                    }

                }
                ##
                y_lab <-  if (metric == "normalised") {

                    "Scaling relative to own lowest-thread config"

                } else if (metric == "serial_speedup") {

                    expression("Speed-up"~~S == N[chains] %*% T[0]^{eq} / T)

                } else {

                        expression(Adj.~ratio~(N[chains] / sec) %*% min[time])

                }
                ##
                colour_scale <-  shared_colour_scale(ps2_df$Algorithm_label)
                ##
                plot_list <-  list()
                ##
                for (dev in c("HPC", "Laptop")) {

                    ##
                    df_dev <-  df_plot %>% filter(device == dev)
                    ##
                    if (nrow(df_dev) == 0) next
                    ##
                    ## Show useful budget labels and endpoints; retain all measured points.
                    tick_candidates <-  if (dev == "HPC") c(1, 2, 4, 8, 16, 32, 64, 96, 128, 180) else c(1, 2, 4, 8, 16)
                    ##
                    x_breaks <-  sort(unique(c(range(df_dev$n_threads),
                                               intersect(tick_candidates, df_dev$n_threads))))
                    ##
                    ## Each normalised reference starts from that algorithm's own thread baseline.
                    ## Adjusted retains the historical shared reference from the best base value.
                    if (metric == "normalised") {

                        ref_line <-  df_dev %>%
                            distinct(device, N_num, N_label, Algorithm_label, n_threads) %>%
                            group_by(device, N_num, N_label, Algorithm_label) %>%
                            dplyr::mutate(perfect = n_threads / min(n_threads)) %>%
                            ungroup()

                    } else {

                        ref_data <-  df_dev %>%
                            group_by(device, N_num, N_label, Algorithm_label) %>%
                            filter(n_threads == min(n_threads)) %>%
                            dplyr::summarise(reference_slope = max(y_val / n_threads), .groups = "drop") %>%
                            group_by(device, N_num, N_label) %>%
                            dplyr::summarise(reference_slope = max(reference_slope), .groups = "drop")
                        ##
                        ref_line <-  df_dev %>%
                            distinct(device, N_num, N_label, n_threads) %>%
                            left_join(ref_data, by = c("device", "N_num", "N_label")) %>%
                            dplyr::mutate(perfect = reference_slope * n_threads)

                    }
                    ##
                    ## The speed-up reference is shared by every arm: ideal scaling is S = N_threads, from 1 thread.
                    if (metric == "serial_speedup") {

                        ref_line <-  dplyr::distinct(df_dev, device, N_num, N_label, n_threads)
                        ref_line$perfect <-  ref_line$n_threads

                    }
                    ##
                    reference_layer <-  if (metric == "normalised") {

                        geom_line( data = ref_line,
                                   aes( x = n_threads,
                                        y = perfect,
                                        colour = Algorithm_label,
                                        group = Algorithm_label),
                                   inherit.aes = FALSE,
                                   show.legend = FALSE,
                                   linewidth = 1,
                                   linetype = "dashed")

                    } else {

                        geom_line( data = ref_line,
                                   aes(x = n_threads, y = perfect),
                                   inherit.aes = FALSE,
                                   linewidth = 1,
                                   linetype = "dashed",
                                   colour = "grey40")

                    }
                    ##
                    p <-  ggplot( df_dev,
                                  aes( x = n_threads,
                                       y = y_val,
                                       colour = Algorithm_label,
                                       group  = Algorithm_label)) +
                        geom_point(size = 5) +
                        geom_line(linewidth = 2) +
                        reference_layer +
                        theme_bw(base_size = 28) +
                        theme( legend.position = ifelse(dev == "Laptop", "bottom", "none"),
                               legend.text = element_text(size = 20),
                               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                        colour_scale +
                        guides(colour = guide_legend(title = NULL, ncol = 2)) +
                        ylab(y_lab) +
                        xlab(expression(log[2](N[threads]~total))) +
                        scale_x_continuous(breaks = x_breaks, trans = "log2") +
                        facet_wrap(~ N_label, scales = "free") +
                        ggtitle(ifelse(dev == "HPC", "Local HPC", "Laptop"))
                    ##
                    plot_list[[dev]] <-  p

                }
                ##
                if (length(plot_list) == 0) {

                    message("No matching observations available for this plot; no file written.")
                    ##
                    return(invisible(NULL))

                }
                ##
                if (length(plot_list) == 2) {

                    combined_plot <-  plot_list[["HPC"]] + plot_list[["Laptop"]] + plot_layout(ncol = 1)

                } else {

                    combined_plot <-  plot_list[[1]] + theme(legend.position = "bottom")

                }
                ##
                print(combined_plot)
                ##
                if (save_plot) {

                    ##
                    ## ---- The file name records the metric, so the old "adjusted" figure is never overwritten:
                    ##      (metric = "serial_speedup" keeps the manuscript file name of figure:ps2_parallel_scalability_plot_1_adj_scalability)
                    ##
                    scaling_figure_filename <-  if (metric == "normalised") {

                        "Figure_ps2_plot_2_norm_scalability.png"

                    } else {

                        "Figure_ps2_plot_2_adj_scalability.png"

                    }
                    ##
                    ggsave( file.path(output_path, scaling_figure_filename),
                            combined_plot,
                            width = 16,
                            height = ifelse(length(plot_list) == 2, 24, 12),
                            dpi = 100)

                }
                ##
                return(combined_plot)

        }

        ## ---- From
        ## legacy/ps_2_parallel_scaling_vs_Mplus_Stan/functions/R_fns_ps2.R:1331
        get_val <-  function( df,
                              dev,
                              N_val,
                              alg_label_df,
                              threads,
                              col
        ) {

                v <-  df %>% filter( device == dev,
                                     N_num == N_val,
                                     Algorithm_label == alg_label_df,
                                     n_threads == threads) %>%
                    pull(all_of(col))
                ##
                if (length(v) == 0) NA_real_ else v[1]

        }
        ##
        ## ---- WCP arms in the ps2 tables: they have no 1-chain run of their own (WCP needs >= 2 threads per
        ##      chain), they run MORE long-run iterations than their serial counterparts, and on the HPC their grids stop at
        ##      176 threads (no 180-thread cell). Previously every one of their cells was therefore "---".
        ##      get_serial_time_equivalent(): the arm's own value at the serial budget if it exists; otherwise the counterpart
        ##      arm's serial time scaled PER ITERATION to the WCP arm's iteration count at the target budget.
        ##      get_stand_in_budget(): the requested budget if the arm measured it; otherwise its largest measured budget above
        ##      the lower budget (rendered with an asterisk, and said in the captions).
        ##
        wcp_counterpart_label <-  c( "BayesMVP + chunking + WCP"                   = "BayesMVP + chunking",
                                     "Stan model (NicoStan) + tape chunking + WCP" = "Stan model (NicoStan) + tape chunking",
                                     "Mplus + WCP"                                 = "Mplus")
        ##
        get_serial_time_equivalent <-  function( df,
                                                 dev,
                                                 N_val,
                                                 alg_label_df,
                                                 t_serial,
                                                 threads_target
        ) {

                own_time <-  get_val(df, dev, N_val, alg_label_df, t_serial, "time_avg")
                if (!is.na(own_time) || !alg_label_df %in% names(wcp_counterpart_label)) return(own_time)
                ##
                counterpart_label <-  wcp_counterpart_label[[alg_label_df]]
                counterpart_time  <-  get_val(df, dev, N_val, counterpart_label, t_serial, "time_avg")
                counterpart_iter  <-  get_val(df, dev, N_val, counterpart_label, t_serial, "N_iter")
                target_iter       <-  get_val(df, dev, N_val, alg_label_df, threads_target, "N_iter")
                ##
                if (is.na(counterpart_time) || is.na(counterpart_iter) || is.na(target_iter) || counterpart_iter <= 0) return(NA_real_)
                return(counterpart_time * target_iter / counterpart_iter)

        }
        ##
        get_stand_in_budget <-  function( df,
                                          dev,
                                          N_val,
                                          alg_label_df,
                                          threads,
                                          threads_lower
        ) {

                if (!is.na(get_val(df, dev, N_val, alg_label_df, threads, "time_avg"))) return(threads)
                ##
                available <-  df %>% filter( device == dev,
                                             N_num == N_val,
                                             Algorithm_label == alg_label_df,
                                             n_threads > threads_lower,
                                             n_threads < threads) %>%
                    pull(n_threads)
                ##
                if (length(available) == 0) return(NA_real_)
                return(max(available))

        }
        ##
        ## ---- get_serial_speedup(): S = N_chains x T0_eq / T and E = S / N_threads (equation eq:paper1_serial_efficiency) for every
        ##      selected configuration, using the same one-chain, one-thread reference as the ps2 tables: the arm's own selected
        ##      configuration at t_serial threads, or, for a WCP arm, the counterpart arm's one at t_serial threads scaled per
        ##      iteration to the WCP arm's N_iter (get_serial_time_equivalent()). A missing reference stays NA.
        ##
        get_serial_speedup <-  function( df,
                                         t_serial = 1
        ) {

                df <-  as.data.frame(df)
                df$serial_time_equivalent <-  NA_real_
                ##
                for (row_index in seq_len(nrow(df))) {

                    df$serial_time_equivalent[row_index] <-  get_serial_time_equivalent( df,
                                                                                          df$device[row_index],
                                                                                          df$N_num[row_index],
                                                                                          df$Algorithm_label[row_index],
                                                                                          t_serial,
                                                                                          df$n_threads[row_index])

                }
                ##
                valid_time <-  is.finite(df$time_avg) & df$time_avg > 0 & is.finite(df$serial_time_equivalent) & df$serial_time_equivalent > 0
                df$serial_speedup_S <-  ifelse(valid_time, df$N_chains * df$serial_time_equivalent / df$time_avg, NA_real_)
                df$parallel_efficiency_E <-  df$serial_speedup_S / df$n_threads
                ##
                return(df)

        }

        ## ---- From
        ## legacy/ps_2_parallel_scaling_vs_Mplus_Stan/functions/R_fns_ps2.R:1352
        make_ratio_table_tex <-  function( df,
                                           dev,
                                           t_lo,
                                           t_hi,
                                           N_vals,
                                           caption,
                                           label,
                                           file_path,
                                           configuration_summary = FALSE,
                                           row_labels = NULL
        ) {

                ## Reuse the manuscript's grouped-N, three-column blocks for within-sampler configuration summaries.
                if (configuration_summary) {

                    selected <-  df[df$device == dev, , drop = FALSE]
                    if (is.null(row_labels)) row_labels <-  unique(selected$Algorithm_label)
                    row_labels <-  row_labels[row_labels %in% selected$Algorithm_label]
                    lines <-  c("\\begin{table}[H]", "\\centering", "\\footnotesize",
                                paste0("\\caption{\\footnotesize{", caption, "}}"), paste0("\\label{", label, "}"),
                                "\\resizebox{\\linewidth}{!}{%",
                                paste0("\\begin{tabular}{l ", paste(rep("ccc", length(N_vals)), collapse = " "), "}"),
                                "\\toprule",
                                paste0(" & ", paste(paste0("\\multicolumn{3}{c}{$N = ",
                                    fn_paper1_format_number_commas_from_10000(N_vals), "$}"), collapse = " & "), " \\\\"),
                                paste(paste0("\\cmidrule(lr){", 2 + 3 * (seq_along(N_vals) - 1), "-",
                                    4 + 3 * (seq_along(N_vals) - 1), "}"), collapse = " "),
                                paste0("Configuration & ", paste(rep("Chunks/WCP & Chains/threads & Iter./sec.",
                                    length(N_vals)), collapse = " & "), " \\\\"), "\\midrule")
                    ##
                    for (row_label in row_labels) {

                        cells <-  character()
                        for (N_value in N_vals) {

                            row <-  selected[selected$Algorithm_label == row_label & selected$N_num == N_value, , drop = FALSE]
                            if (!nrow(row)) {

                                cells <-  c(cells, "---", "---", "---")

                            } else {

                                row <-  row[order(-row$total_iter_per_sec, row$n_threads, row$num_chunks), , drop = FALSE][1, ]
                                rate <-  fmt3(row$total_iter_per_sec)
                                if (row$total_iter_per_sec == max(selected$total_iter_per_sec[selected$N_num == N_value])) {

                                    rate <-  paste0("$\\mathbf{", rate, "}$")

                                }
                                cells <-  c(cells, paste0(row$num_chunks, "/", row$threads_per_chain),
                                            paste0(row$n_chains, "/", row$n_threads), rate)

                            }

                        }
                        escaped_label <-  gsub(pattern = "_", replacement = "\\_", x = row_label, fixed = TRUE)
                        lines <-  c(lines, paste0(escaped_label, " & ", paste(cells, collapse = " & "), " \\\\"))

                    }
                    lines <-  c(lines, "\\bottomrule", "\\end{tabular}", "}", "\\end{table}")
                    emit_tex(lines, file_path)
                    return(invisible(file_path))

                }
                ##
                ##
                M <-  function() matrix(NA_real_, length(manus_rows), length(N_vals), dimnames = list(manus_rows, N_vals))
                ##
                cell_time_lo <-  M()
                ##
                cell_time_hi <-  M()
                ##
                cell_chains_lo <-  M()
                ##
                cell_chains_hi <-  M()
                ##
                for (a in manus_rows) for (j in seq_along(N_vals)) {

                    ## WCP arms: counterpart's serial run scaled per iteration when the arm has no own value at t_lo.
                    cell_time_lo[a, j] <-  get_serial_time_equivalent(df, dev, N_vals[j], manus_to_df_label[[a]], t_lo, t_hi)
                    ##
                    cell_time_hi[a, j] <-  get_val(df, dev, N_vals[j], manus_to_df_label[[a]], t_hi, "time_avg")
                    ##
                    own_chains_lo <-  get_val(df, dev, N_vals[j], manus_to_df_label[[a]], t_lo, "N_chains")
                    cell_chains_lo[a, j] <-  if (is.na(own_chains_lo)) 1 else own_chains_lo   ## the serial equivalent is one chain
                    ##
                    cell_chains_hi[a, j] <-  get_val(df, dev, N_vals[j], manus_to_df_label[[a]], t_hi, "N_chains")

                }
                ##
                ## ---- Ratio per THREAD: (time / chains x threads) at the higher budget over the same at the lower budget.
                ##      For the one-chain-per-thread arms this is exactly time_hi / time_lo as before; for the WCP arms (fewer chains
                ##      than threads) the raw time ratio would reward using 96 threads for 16 chains, so the ratio is scaled by
                ##      threads / chains and 1 still means perfect throughput scaling.
                ##
                cell_ratio <-  (cell_time_hi / cell_chains_hi * t_hi) / (cell_time_lo / cell_chains_lo * t_lo)
                ##
                print(cell_ratio)
                ##
                lines <-  c(
                             "\\begin{table}[H]", "\\centering", "\\footnotesize",
                             paste0("\\caption{\\footnotesize{", caption, "}}"),
                             paste0("\\label{", label, "}"),
                             paste0("\\begin{tabular}{l ", paste(rep("ccc", length(N_vals)), collapse = " "), "}"),
                             "\\toprule",
                             paste0( " & ",
                                     paste( paste0( "\\multicolumn{3}{c}{$N = ",
                                                    fn_paper1_format_number_commas_from_10000(N_vals), "$}"),
                                            collapse = " & "),
                                     " \\\\"),
                             paste( paste0( "\\cmidrule(lr){",
                                            2 + 3 * (seq_along(N_vals) - 1), "-",
                                            4 + 3 * (seq_along(N_vals) - 1), "}"),
                                    collapse = " "),
                             paste0( "Algorithm & ",
                                     paste(rep(paste0(t_lo, "t & ", t_hi, "t & Ratio"), length(N_vals)), collapse = " & "),
                                     " \\\\"),
                             "\\midrule"
                )
                ##
                for (a in manus_rows) {

                    cells <-  c()
                    ##
                    for (j in seq_along(N_vals)) {

                        if (is.na(cell_ratio[a, j])) {

                            cells <-  c(cells, "---", "---", "---")

                        } else {

                            is_best   <-  (cell_ratio[a, j] == min(cell_ratio[, j], na.rm = TRUE))
                            ##
                            ratio_str <-  if (is_best) paste0("$\\mathbf{", fmt3(cell_ratio[a, j]), "\\times}$")
                                else              paste0("$", fmt3(cell_ratio[a, j]), "\\times$")
                            ##
                            cells <-  c(cells, fmt3(cell_time_lo[a, j]), fmt3(cell_time_hi[a, j]), ratio_str)

                        }

                    }
                    ##
                    lines <-  c(lines, paste0(a, " & ", paste(cells, collapse = " & "), " \\\\"))

                }
                ##
                lines <-  c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
                ##
                writeLines(lines, file_path)
                ##
                message("Wrote: ", file_path)
                ##
                return(cat(lines, sep = "\n"))

        }

        ## ---- From
        ## legacy/ps_2_parallel_scaling_vs_Mplus_Stan/functions/R_fns_ps2.R:1418
        make_smt_gain_table_tex <-  function( df,
                                              dev,
                                              t_phys,
                                              t_smt,
                                              N_vals,
                                              caption,
                                              label,
                                              file_path
        ) {

                ##
                M <-  function() matrix(NA_real_, length(manus_rows), length(N_vals), dimnames = list(manus_rows, N_vals))
                ##
                cell_time_p <-  M()
                ##
                cell_time_s <-  M()
                ##
                cell_gain <-  M()
                ##
                cell_stand_in <-  matrix(FALSE, length(manus_rows), length(N_vals), dimnames = list(manus_rows, N_vals))
                ##
                for (a in manus_rows) for (j in seq_along(N_vals)) {

                    ## WCP arms with no cell at t_smt use their largest measured budget above t_phys (asterisk):
                    smt_budget <-  get_stand_in_budget(df, dev, N_vals[j], manus_to_df_label[[a]], t_smt, t_phys)
                    cell_stand_in[a, j] <-  !is.na(smt_budget) && smt_budget != t_smt
                    ##
                    tp <-  get_val(df, dev, N_vals[j], manus_to_df_label[[a]], t_phys, "time_avg")
                    ##
                    ts <-  if (is.na(smt_budget)) NA_real_ else get_val(df, dev, N_vals[j], manus_to_df_label[[a]], smt_budget, "time_avg")
                    ##
                    cp <-  get_val(df, dev, N_vals[j], manus_to_df_label[[a]], t_phys, "N_chains")
                    ##
                    cs <-  if (is.na(smt_budget)) NA_real_ else get_val(df, dev, N_vals[j], manus_to_df_label[[a]], smt_budget, "N_chains")
                    ##
                    cell_time_p[a, j] <-  tp
                    ##
                    cell_time_s[a, j] <-  ts
                    ##
                    if (!is.na(tp) && !is.na(ts)) cell_gain[a, j] <-  100 * ((cs / ts) / (cp / tp) - 1)

                }
                ##
                print(cell_time_p)
                ##
                lines <-  c(
                             "\\begin{table}[H]", "\\centering", "\\footnotesize",
                             paste0("\\caption{\\footnotesize{", caption, "}}"),
                             paste0("\\label{", label, "}"),
                             paste0("\\begin{tabular}{l ", paste(rep("ccc", length(N_vals)), collapse = " "), "}"),
                             "\\toprule",
                             paste0( " & ",
                                     paste( paste0( "\\multicolumn{3}{c}{$N = ",
                                                    fn_paper1_format_number_commas_from_10000(N_vals), "$}"),
                                            collapse = " & "),
                                     " \\\\"),
                             paste( paste0( "\\cmidrule(lr){",
                                            2 + 3 * (seq_along(N_vals) - 1), "-",
                                            4 + 3 * (seq_along(N_vals) - 1), "}"),
                                    collapse = " "),
                             paste0( "Algorithm & ",
                                     paste(rep(paste0(t_phys, "t & ", t_smt, "t & Gain"), length(N_vals)), collapse = " & "),
                                     " \\\\"),
                             "\\midrule"
                )
                ##
                for (a in manus_rows) {

                    cells <-  c()
                    ##
                    for (j in seq_along(N_vals)) {

                        if (is.na(cell_gain[a, j])) {

                            cells <-  c(cells, "---", "---", "---")

                        } else {

                            is_best  <-  (cell_gain[a, j] == max(cell_gain[, j], na.rm = TRUE))
                            ##
                            gain_str <-  if (is_best) paste0("$\\mathbf{", formatC( x = cell_gain[a, j],
                                                                                    format = "f",
                                                                                    digits = 1,
                                                                                    flag = "+",
                                                                                    decimal.mark = "."), "\\%}$")
                                else              paste0("$", formatC( x = cell_gain[a, j],
                                                                       format = "f",
                                                                       digits = 1,
                                                                       flag = "+",
                                                                       decimal.mark = "."), "\\%$")
                            ##
                            cells <-  c(cells, fmt3(cell_time_p[a, j]),
                                        paste0(fmt3(cell_time_s[a, j]), if (cell_stand_in[a, j]) "$^{*}$" else ""), gain_str)

                        }

                    }
                    ##
                    lines <-  c(lines, paste0(a, " & ", paste(cells, collapse = " & "), " \\\\"))

                }
                ##
                lines <-  c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
                ##
                writeLines(lines, file_path)
                ##
                message("Wrote: ", file_path)
                ##
                return(cat(lines, sep = "\n"))

        }
        ##
        ## ---- SMT gain against the 1-chain serial baseline:
        ##
        ## For each algorithm and N: speed-up at t_phys and at t_smt threads relative to the serial (t_serial = 1 chain) run,
        ## speed-up(t) = chain_rate(t) / chain_rate(t_serial) with chain_rate = N_chains / time; their difference (the extra
        ## serial-chain equivalents that SMT buys); that difference per added logical thread (t_smt - t_phys); and the usual
        ## SMT gain in percent, 100 (speed-up(t_smt) / speed-up(t_phys) - 1), which equals the t_phys-vs-t_smt ratio exactly.
        ## Rows come from the same best-per-budget selection as the other ps2 tables, so chunks and chain/WCP splits may
        ## differ between budgets (stated in the caption).
        ##
        make_smt_serial_baseline_table_tex <-  function( df,
                                                         dev,
                                                         t_serial,
                                                         t_phys,
                                                         t_smt,
                                                         N_vals,
                                                         caption,
                                                         label,
                                                         file_path
        ) {

                M <-  function() matrix(NA_real_, length(manus_rows), length(N_vals), dimnames = list(manus_rows, N_vals))
                ##
                cell_speedup_phys <-  M()
                cell_speedup_smt  <-  M()
                cell_difference   <-  M()
                cell_per_thread   <-  M()
                cell_gain         <-  M()
                ##
                chain_rate_at <-  function(N_val, alg, threads) {

                        time_at   <-  get_val(df, dev, N_val, manus_to_df_label[[alg]], threads, "time_avg")
                        chains_at <-  get_val(df, dev, N_val, manus_to_df_label[[alg]], threads, "N_chains")
                        if (is.na(time_at) || is.na(chains_at) || time_at <= 0) NA_real_ else chains_at / time_at

                }
                ##
                cell_stand_in <-  matrix(FALSE, length(manus_rows), length(N_vals), dimnames = list(manus_rows, N_vals))
                ##
                for (a in manus_rows) for (j in seq_along(N_vals)) {

                    ## WCP arms: serial baseline = counterpart's 1-chain run scaled per iteration; largest measured
                    ## budget stands in for t_smt where the arm has no such cell (asterisk).
                    serial_time <-  get_serial_time_equivalent(df, dev, N_vals[j], manus_to_df_label[[a]], t_serial, t_phys)
                    rate_serial <-  if (is.na(serial_time) || serial_time <= 0) NA_real_ else 1 / serial_time
                    rate_phys   <-  chain_rate_at(N_vals[j], a, t_phys)
                    smt_budget  <-  get_stand_in_budget(df, dev, N_vals[j], manus_to_df_label[[a]], t_smt, t_phys)
                    cell_stand_in[a, j] <-  !is.na(smt_budget) && smt_budget != t_smt
                    rate_smt    <-  if (is.na(smt_budget)) NA_real_ else chain_rate_at(N_vals[j], a, smt_budget)
                    ##
                    if (!is.na(rate_serial) && !is.na(rate_phys)) cell_speedup_phys[a, j] <-  rate_phys / rate_serial
                    if (!is.na(rate_serial) && !is.na(rate_smt))  cell_speedup_smt[a, j]  <-  rate_smt / rate_serial
                    ##
                    if (!is.na(cell_speedup_phys[a, j]) && !is.na(cell_speedup_smt[a, j])) {

                        cell_difference[a, j] <-  cell_speedup_smt[a, j] - cell_speedup_phys[a, j]
                        cell_per_thread[a, j] <-  cell_difference[a, j] / (smt_budget - t_phys)
                        cell_gain[a, j]       <-  100 * (cell_speedup_smt[a, j] / cell_speedup_phys[a, j] - 1)

                    }

                }
                ##
                lines <-  c(
                             "\\begin{table}[H]", "\\centering", "\\footnotesize",
                             paste0("\\caption{\\footnotesize{", caption, "}}"),
                             paste0("\\label{", label, "}"),
                             "\\resizebox{\\linewidth}{!}{%",
                             paste0("\\begin{tabular}{l ", paste(rep("ccccc", length(N_vals)), collapse = " "), "}"),
                             "\\toprule",
                             paste0( " & ",
                                     paste( paste0( "\\multicolumn{5}{c}{$N = ",
                                                    fn_paper1_format_number_commas_from_10000(N_vals), "$}"),
                                            collapse = " & "),
                                     " \\\\"),
                             paste( paste0( "\\cmidrule(lr){",
                                            2 + 5 * (seq_along(N_vals) - 1), "-",
                                            6 + 5 * (seq_along(N_vals) - 1), "}"),
                                    collapse = " "),
                             paste0( "Algorithm & ",
                                     paste(rep(paste0("$S_{", t_phys, "}$ & $S_{", t_smt, "}$ & $S_{", t_smt,
                                                      "} - S_{", t_phys, "}$ & per thread & Gain"), length(N_vals)), collapse = " & "),
                                     " \\\\"),
                             "\\midrule"
                )
                ##
                for (a in manus_rows) {

                    cells <-  c()
                    ##
                    for (j in seq_along(N_vals)) {

                        if (is.na(cell_gain[a, j])) {

                            cells <-  c(cells, "---", "---", "---", "---", "---")

                        } else {

                            is_best  <-  (cell_gain[a, j] == max(cell_gain[, j], na.rm = TRUE))
                            ##
                            gain_str <-  if (is_best) paste0("$\\mathbf{", formatC( x = cell_gain[a, j],
                                                                                    format = "f",
                                                                                    digits = 1,
                                                                                    flag = "+",
                                                                                    decimal.mark = "."), "\\%}$")
                                else              paste0("$", formatC( x = cell_gain[a, j],
                                                                       format = "f",
                                                                       digits = 1,
                                                                       flag = "+",
                                                                       decimal.mark = "."), "\\%$")
                            ##
                            cells <-  c(cells,
                                        paste0("$", formatC(x = cell_speedup_phys[a, j], format = "f", digits = 1, decimal.mark = "."), "\\times$"),
                                        paste0(paste0("$", formatC( x = cell_speedup_smt[a, j],
                                                                    format = "f",
                                                                    digits = 1,
                                                                    decimal.mark = "."), "\\times$"), if (cell_stand_in[a, j]) "$^{*}$" else ""),
                                        paste0("$", formatC( x = cell_difference[a, j],
                                                             format = "f",
                                                             digits = 1,
                                                             flag = "+",
                                                             decimal.mark = "."), "$"),
                                        paste0("$", formatC(x = cell_per_thread[a, j], format = "f", digits = 2, decimal.mark = "."), "$"),
                                        gain_str)

                        }

                    }
                    ##
                    lines <-  c(lines, paste0(a, " & ", paste(cells, collapse = " & "), " \\\\"))

                }
                ##
                lines <-  c(lines, "\\bottomrule", "\\end{tabular}", "}", "\\end{table}")
                ##
                writeLines(lines, file_path)
                ##
                message("Wrote: ", file_path)
                ##
                return(invisible(lines))

        }

        ## ---- From legacy/ps_2_parallel_scaling_vs_Mplus_Stan/ps_2_parallel_scaling_vs_Mplus_Stan_v2.R:930
        R_fn_plot_mpmt_throughput <-  function( df,
                                                save_plot = TRUE,
                                                output_path) {

                colour_scale <-  shared_colour_scale(df$mode)
                ##
                plot_list <-  list()
                ##
                for (dev in c("HPC", "Laptop")) {

                    df_dev <-  df %>% filter(device == dev)
                    ##
                    if (nrow(df_dev) == 0) next
                    ##
                    x_breaks <-  sort(unique(df_dev$n_threads))
                    ##
                    p <-  ggplot( df_dev,
                                  aes( x = n_threads,
                                       y = chain_rate,
                                       colour = mode,
                                       linetype = chunk_label,
                                       group = interaction(mode, chunk_label))) +
                        geom_point(size = 4) +
                        geom_line(linewidth = 1.5) +
                        geom_errorbar( width = 0.05,
                                       linewidth = 0.8,
                                       aes( ymin = chain_rate - chain_rate_SD,
                                            ymax = chain_rate + chain_rate_SD)) +
                        theme_bw(base_size = 24) +
                        theme( legend.position = ifelse(dev == "Laptop", "bottom", "none"),
                               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                        colour_scale +
                        guides( colour = guide_legend(title = NULL, nrow = 2),
                                linetype = guide_legend(title = NULL, nrow = 2)) +
                        ylab("Within-method efficiency (chains / time)") +
                        xlab(expression(log[2](N[threads]))) +
                        scale_x_continuous(breaks = x_breaks, trans = "log2") +
                        scale_y_continuous(trans = "log2") +
                        facet_wrap(~ N_label, scales = "free") +
                        ggtitle(ifelse(dev == "HPC", "Local HPC", "Laptop"))
                    ##
                    plot_list[[dev]] <-  p

                }
                ##
                if (length(plot_list) == 0) {

                    message("No matching observations available for this plot; no file written.")
                    ##
                    return(invisible(NULL))

                }
                ##
                combined <-  if (length(plot_list) == 2) plot_list[["HPC"]] + plot_list[["Laptop"]] + plot_layout(ncol = 1)
                    else plot_list[[1]] + theme(legend.position = "bottom")
                ##
                print(combined)
                ##
                if (save_plot) {

                    ggsave( file.path(output_path, "Figure_BayesMVP_MP_vs_MT_throughput.png"),
                            combined,
                            width = 16,
                            height = ifelse(length(plot_list) == 2, 24, 12),
                            dpi = 100)

                }
                ##
                invisible(combined)

        }

        ## ---- From legacy/ps_2_parallel_scaling_vs_Mplus_Stan/ps_2_parallel_scaling_vs_Mplus_Stan_v2.R:976
        R_fn_plot_mpmt_scaling <-  function( df,
                                             save_plot = TRUE,
                                             output_path) {

                colour_scale <-  shared_colour_scale(df$mode)
                ##
                plot_list <-  list()
                ##
                for (dev in c("HPC", "Laptop")) {

                    df_dev <-  df %>% filter(device == dev)
                    ##
                    if (nrow(df_dev) == 0) next
                    ##
                    x_breaks <-  sort(unique(df_dev$n_threads))
                    ##
                    ref_line <-  df_dev %>%
                        distinct(N, N_label, mode, chunk_label, n_threads) %>%
                        group_by(N, N_label, mode, chunk_label) %>%
                        dplyr::mutate(perfect = n_threads / min(n_threads)) %>%
                        ungroup()
                    ##
                    p <-  ggplot( df_dev,
                                  aes( x = n_threads,
                                       y = norm_scaling,
                                       colour = mode,
                                       linetype = chunk_label,
                                       group = interaction(mode, chunk_label))) +
                        geom_point(size = 4) +
                        geom_line(linewidth = 1.5) +
                        geom_line( data = ref_line,
                                   aes( x = n_threads,
                                        y = perfect,
                                        colour = mode,
                                        group = interaction(mode, chunk_label)),
                                   inherit.aes = FALSE,
                                   show.legend = FALSE,
                                   linewidth = 1,
                                   linetype = "dashed") +
                        theme_bw(base_size = 24) +
                        theme( legend.position = ifelse(dev == "Laptop", "bottom", "none"),
                               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                        colour_scale +
                        guides( colour = guide_legend(title = NULL, nrow = 2),
                                linetype = guide_legend(title = NULL, nrow = 2)) +
                        ylab("Scaling relative to own lowest-thread config") +
                        xlab(expression(log[2](N[threads]))) +
                        scale_x_continuous(breaks = x_breaks, trans = "log2") +
                        scale_y_continuous(trans = "log2") +
                        facet_wrap(~ N_label, scales = "free") +
                        ggtitle(ifelse(dev == "HPC", "Local HPC", "Laptop"))
                    ##
                    plot_list[[dev]] <-  p

                }
                ##
                if (length(plot_list) == 0) {

                    message("No matching observations available for this plot; no file written.")
                    ##
                    return(invisible(NULL))

                }
                ##
                combined <-  if (length(plot_list) == 2) plot_list[["HPC"]] + plot_list[["Laptop"]] + plot_layout(ncol = 1)
                    else plot_list[[1]] + theme(legend.position = "bottom")
                ##
                print(combined)
                ##
                if (save_plot) {

                    ggsave( file.path(output_path, "Figure_BayesMVP_MP_vs_MT_scaling.png"),
                            combined,
                            width = 16,
                            height = ifelse(length(plot_list) == 2, 24, 12),
                            dpi = 100)

                }
                ##
                invisible(combined)

        }

        ##
        return(environment())

}






















