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

        ## ---- Cache-capacity lines on the chunk-count figures (links the E4 cache/RAM mechanism to E1/E2) -----------------------------
        ##
        ## Cache sizes (bytes) of the two machines: L1d and L2 per physical core; S_L3 = size of one L3 cache, N_L3 = number of L3
        ## caches, C_L3 = physical cores sharing each L3 cache and s = hardware threads per core (as in equation eq:paper1_auto_n_chunks).
        ##
        paper1_cache_specs <-  list( HPC    = list( physical_cores = 96,                   ## AMD EPYC 9654
                                                    L1d_per_core   = 32 * 1024,
                                                    L2_per_core    = 1024 * 1024,
                                                    S_L3           = 32 * 1024 * 1024,     ## per CCD
                                                    N_L3           = 12,
                                                    C_L3           = 8,
                                                    s              = 2),
                                     Laptop = list( physical_cores = 8,                    ## AMD Ryzen 7 5800H
                                                    L1d_per_core   = 32 * 1024,
                                                    L2_per_core    = 512 * 1024,
                                                    S_L3           = 16 * 1024 * 1024,
                                                    N_L3           = 1,
                                                    C_L3           = 8,
                                                    s              = 2))
        ##
        ## Working set of lp_grad per individual (row) for NicoStan+BayesMVP (binary LC-MVP, 2 classes, T = 6): 31 T + 15 doubles = 1,608 bytes.
        paper1_bytes_per_row_BayesMVP <-  8 * (31 * 6 + 15)
        ##
        ## Working set per individual (row) of one reduce_sum_static partial sum of the Stan model used by AD_Stan_tape_chunked,
        ## AD_Stan_WCP and AD_Stan_WCP_chunking (LC_MVP_bin_PartialLog_v5_reduce_sum_static.stan, 2 classes, T = 6), measured from
        ## single-threaded gradient evaluations (Stan Math, stanc 2.36.0; exactly linear in the rows of the partial sum):
        ## nested autodiff tape (arena) 15,936 + autodiff stacks 2,360 + heap temporaries 804 + data rows read 52 = 19,152 bytes.
        paper1_bytes_per_row_Stan <-  15936 + 2360 + 804 + 52

        ## ---- Cache capacity per active thread, for a total number of active threads (chains x WCP threads per chain):
        ##      L3 per active thread = S_L3 / max(C_L3, A_L3), A_L3 = min(s C_L3, ceiling(N_threads / N_L3)) (equation eq:paper1_auto_n_chunks);
        ##      SMT is in use when the active threads exceed the physical cores, and then L1d and L2 per thread are divided by s.
        ## ---- 2026-10-04: L3 per active thread = S_L3 / A_L3 (no max with C_L3), as in equations
        ##      eq:paper1_auto_n_chunks_active_threads and eq:paper1_auto_n_chunks_target: with fewer active
        ##      threads than cores on one L3 cache, each active thread has more than the per-core share (e.g.
        ##      32 MB with one thread per local-HPC CCD). Unchanged at 96-180 threads (local-HPC) and 8-16
        ##      threads (laptop).
        ##
        fn_paper1_cache_per_thread <-  function( device,
                                                 n_threads_total,
                                                 cache_specs = paper1_cache_specs
        ) {

                spec <-  cache_specs[[device]]
                ##
                if (is.null(spec)) stop(paste0("fn_paper1_cache_per_thread: no cache specification for device = '", device, "'."))
                ##
                SMT <-  n_threads_total > spec$physical_cores
                ##
                A_L3 <-  pmin(spec$s * spec$C_L3, ceiling(n_threads_total / spec$N_L3))
                ##
                data.frame( device           = device,
                            n_threads_total  = n_threads_total,
                            SMT              = SMT,
                            # L3               = spec$S_L3 / pmax(spec$C_L3, A_L3),
                            ## (2026-10-04: no cap at the per-core share; see the note above this function)
                            L3               = spec$S_L3 / A_L3,
                            L2               = spec$L2_per_core  / ifelse(SMT, spec$s, 1),
                            L1               = spec$L1d_per_core / ifelse(SMT, spec$s, 1),
                            stringsAsFactors = FALSE)

        }

        ## ---- Chunk count c* = bytes_per_row x N / (cache per active thread) at which one chunk first fits in each cache level.
        ##      Only thresholds strictly inside chunk_range (the tested chunk counts of a panel) are returned; to the left of the
        ##      L3 line, one chunk exceeds the L3 cache per active thread (and spills to RAM).
        ##
        fn_paper1_cache_threshold_lines <-  function( N,
                                                      chunk_range,
                                                      device,
                                                      n_threads_total,
                                                      bytes_per_row = paper1_bytes_per_row_BayesMVP,
                                                      cache_levels  = c("L3", "L2", "L1"),
                                                      cache_specs   = paper1_cache_specs
        ) {

                per_thread <-  fn_paper1_cache_per_thread( device          = device,
                                                           n_threads_total = n_threads_total,
                                                           cache_specs     = cache_specs)
                ##
                thresholds <-  do.call(rbind, lapply(cache_levels, function(cache_level) {

                        data.frame( N                      = N,
                                    device                 = device,
                                    n_threads_total        = per_thread$n_threads_total,
                                    SMT                    = per_thread$SMT,
                                    level                  = cache_level,
                                    cache_bytes_per_thread = per_thread[[cache_level]],
                                    c_star                 = bytes_per_row * N / per_thread[[cache_level]],
                                    stringsAsFactors       = FALSE)

                }))
                ##
                inside <-  thresholds$c_star > min(chunk_range) & thresholds$c_star < max(chunk_range)
                ##
                thresholds <-  thresholds[inside, , drop = FALSE]
                ##
                rownames(thresholds) <-  NULL
                ##
                thresholds

        }

        ## ---- Rows for the cache-line labels at the top of a panel: a label moves down one row whenever it would overlap a label
        ##      already placed in that row (label widths estimated from the number of characters, in log10 units of the chunk axis):
        ##
        fn_paper1_cache_label_rows <-  function( x,
                                                 labels,
                                                 chunk_range,
                                                 panel_width_in,
                                                 label_size = 4.3,
                                                 label_x = NULL
        ) {

                log10_span <-  1.1 * diff(log10(range(chunk_range)))
                ##
                ## ~0.0215 inches per character per mm of text size, plus the label padding and a small gap:
                half_width <-  0.5 * (nchar(labels) * 0.0215 * label_size + 0.15) / panel_width_in * log10_span
                ##
                ## label centres (by default on the lines themselves; see fn_paper1_cache_label_x for labels kept inside the panel):
                if (is.null(label_x)) label_x <-  x
                ##
                rows <-  rep(NA_real_, length(x))
                ##
                for (i in order(x)) {

                    row <-  0
                    ##
                    repeat {

                        in_row <-  which(rows == row)
                        ##
                        if (!any(abs(log10(label_x[i]) - log10(label_x[in_row])) < half_width[i] + half_width[in_row])) break
                        ##
                        row <-  row + 1

                    }
                    ##
                    rows[i] <-  row

                }
                ##
                rows

        }

        ## ---- Centres for the cache-line labels that keep each label inside the panel: a label centred on its line stays there,
        ##      and a label that would run past either end of the log10 chunk axis (default 5% expansion on each side) moves inwards
        ##      just far enough to fit (label widths estimated as in fn_paper1_cache_label_rows):
        ##
        fn_paper1_cache_label_x <-  function( x,
                                              labels,
                                              chunk_range,
                                              panel_width_in,
                                              label_size = 4.3
        ) {

                log10_range <-  log10(range(chunk_range))
                ##
                log10_span <-  1.1 * diff(log10_range)
                ##
                half_width <-  0.5 * (nchar(labels) * 0.0215 * label_size + 0.15) / panel_width_in * log10_span
                ##
                panel_left  <-  log10_range[1] - 0.05 * diff(log10_range)
                panel_right <-  log10_range[2] + 0.05 * diff(log10_range)
                ##
                label_x <-  x
                ##
                ## (a label wider than the whole panel stays on its line)
                fits_in_panel <-  2 * half_width < panel_right - panel_left
                ##
                past_left  <-  fits_in_panel & log10(x) < panel_left + half_width
                past_right <-  fits_in_panel & log10(x) > panel_right - half_width
                ##
                label_x[past_left]  <-  10^(panel_left + half_width[past_left])
                label_x[past_right] <-  10^(panel_right - half_width[past_right])
                ##
                label_x

        }

        ## ---- Breaks for a log10 chunk axis with free facet scales: the tested chunk counts of the panel whose tested range best
        ##      matches the panel limits (the breaks function is called once per panel, with that panel's limits):
        ##
        fn_paper1_log10_chunk_breaks <-  function( chunks_by_panel ) {

                panel_log10_ranges <-  t(vapply(chunks_by_panel, function(chunks) log10(range(chunks)), numeric(2)))
                ##
                function(limits) {

                        if (any(!is.finite(limits))) return(sort(unique(unlist(chunks_by_panel))))
                        ##
                        distance <-  abs(panel_log10_ranges[, 1] - log10(limits[1])) + abs(panel_log10_ranges[, 2] - log10(limits[2]))
                        ##
                        closest <-  which(abs(distance - min(distance)) < 1e-9)
                        ##
                        chunks <-  sort(unique(unlist(chunks_by_panel[closest])))
                        ##
                        chunks[chunks >= limits[1] & chunks <= limits[2]]

                }

        }

        ## ---- Labels for a log10 chunk axis: every break in priority_chunks is labelled; any other break is labelled only if it is at
        ##      least min_gap_log10 (log10 units) from every labelled break, so neighbouring tick labels do not overlap. Breaks in
        ##      secondary_chunks are considered before the remaining breaks:
        ##
        fn_paper1_log10_chunk_labels <-  function( breaks,
                                                   priority_chunks,
                                                   min_gap_log10,
                                                   secondary_chunks = NULL
        ) {

                keep <-  breaks %in% priority_chunks
                ##
                for (i in which(!keep & breaks %in% secondary_chunks)) {

                    if (all(abs(log10(breaks[i]) - log10(breaks[keep])) >= min_gap_log10)) keep[i] <-  TRUE

                }
                ##
                for (i in which(!keep)) {

                    if (all(abs(log10(breaks[i]) - log10(breaks[keep])) >= min_gap_log10)) keep[i] <-  TRUE

                }
                ##
                chunk_labels <-  fn_paper1_format_number_commas_from_10000(breaks)
                ##
                chunk_labels[!keep] <-  ""
                ##
                ## two labelled breaks closer than min_gap_log10 (both in priority_chunks): the label of the second is padded with
                ## trailing spaces, which moves it further from the axis (rotated tick labels), so the two labels do not overlap:
                kept <-  which(keep)
                ##
                for (j in seq_along(kept)[-1]) {

                    if (abs(log10(breaks[kept[j]]) - log10(breaks[kept[j - 1]])) < min_gap_log10 && !grepl(" $", chunk_labels[kept[j - 1]])) {

                        chunk_labels[kept[j]] <-  paste0(chunk_labels[kept[j]], strrep(" ", 2 * nchar(chunk_labels[kept[j - 1]]) + 1))

                    }

                }
                ##
                chunk_labels

        }

        ## ---- Layers for the cache-capacity lines and their labels. With colour_column/linetype_column the lines take the colour
        ##      and line type of the data line they belong to; otherwise they are drawn in grey:
        ##
        fn_paper1_cache_line_layers <-  function( thresholds,
                                                  colour_column   = NULL,
                                                  linetype_column = NULL,
                                                  label_size      = 4.3
        ) {

                if (is.null(thresholds) || nrow(thresholds) == 0) return(NULL)
                ##
                ## grey labels are centred on label_x when the thresholds carry it (labels kept inside the panel), otherwise on the line:
                label_x_column <-  if ("label_x" %in% names(thresholds)) "label_x" else "c_star"
                ##
                ## labels in the paper's notation (plotmath, e.g. N["threads/chain"]) are drawn with parse = TRUE when the thresholds carry them:
                label_column <-  if ("label_expression" %in% names(thresholds)) "label_expression" else "label"
                ##
                if (is.null(colour_column)) {

                    list( ggplot2::geom_vline( data = thresholds,
                                               mapping = ggplot2::aes(xintercept = c_star),
                                               colour = "grey45",
                                               linetype = "22",
                                               linewidth = 0.6,
                                               alpha = 0.8),
                          ggplot2::geom_label( data = thresholds,
                                               ## mapping = ggplot2::aes(x = .data[[label_x_column]], y = Inf, label = label, vjust = vjust),
                                               mapping = ggplot2::aes(x = .data[[label_x_column]], y = Inf, label = .data[[label_column]], vjust = vjust),
                                               parse = identical(label_column, "label_expression"),
                                               inherit.aes = FALSE,
                                               colour = "grey25",
                                               size = label_size,
                                               label.size = 0,
                                               fill = "white",
                                               label.padding = ggplot2::unit(0.08, "lines")))

                } else {

                    list( ggplot2::geom_vline( data = thresholds,
                                               mapping = ggplot2::aes( xintercept = c_star,
                                                                       colour     = .data[[colour_column]],
                                                                       linetype   = .data[[linetype_column]]),
                                               linewidth = 0.6,
                                               alpha = 0.55,
                                               show.legend = FALSE),
                          ggplot2::geom_label( data = thresholds,
                                               mapping = ggplot2::aes( x      = c_star,
                                                                       y      = Inf,
                                                                       label  = label,
                                                                       colour = .data[[colour_column]],
                                                                       vjust  = vjust),
                                               inherit.aes = FALSE,
                                               size = label_size,
                                               label.size = 0,
                                               fill = "white",
                                               label.padding = ggplot2::unit(0.08, "lines"),
                                               show.legend = FALSE))

                }

        }

        ## ---- Thread counts on the E2 scaling figure after which SMT is in use and the L3 cache per active thread falls below the
        ##      L3 cache per core (equation eq:paper1_auto_n_chunks); changes at the same thread count share one marker:
        ##
        fn_paper1_thread_markers <-  function( device,
                                               n_threads_max,
                                               cache_specs = paper1_cache_specs
        ) {

                spec <-  cache_specs[[device]]
                ##
                per_thread <-  fn_paper1_cache_per_thread( device          = device,
                                                           n_threads_total = seq_len(n_threads_max),
                                                           cache_specs     = cache_specs)
                ##
                L3_per_core <-  spec$S_L3 / spec$C_L3
                ##
              # markers <-  data.frame( device          = device,
              #                         after_n_threads = c( max(per_thread$n_threads_total[!per_thread$SMT]),
              #                                              max(per_thread$n_threads_total[per_thread$L3 >= L3_per_core])),
              #                         change          = c( "SMT",
              #                                              paste0("L3 per thread < ", L3_per_core / (1024 * 1024), " MB")),
              #                         stringsAsFactors = FALSE)
                ##
                ## ---- 2026-10-06: the marker also names the L2 cache per thread, which falls below the L2 cache
                ##      per core once SMT is in use (1 MB local-HPC, 512 KB laptop); the two cache sizes share one
                ##      line, "per thread: L2 < ..., L3 < ...", so the label stays on two lines:
                fn_format_cache_bytes <-  function(bytes) {

                        if (bytes >= 1024 * 1024) return(paste0(bytes / (1024 * 1024), " MB"))
                        ##
                        return(paste0(bytes / 1024, " KB"))

                }
                ##
                n_threads_total <-  per_thread$n_threads_total
                L2_per_core     <-  spec$L2_per_core
                ##
                L2_per_thread_change <-  paste0("L2 < ", fn_format_cache_bytes(L2_per_core))
                L3_per_thread_change <-  paste0("L3 < ", fn_format_cache_bytes(L3_per_core))
                ##
                markers <-  data.frame( device           = device,
                                        after_n_threads  = c( max(n_threads_total[!per_thread$SMT]),
                                                              max(n_threads_total[per_thread$L2 >= L2_per_core]),
                                                              max(n_threads_total[per_thread$L3 >= L3_per_core])),
                                        change           = c("SMT", L2_per_thread_change, L3_per_thread_change),
                                        per_thread_cache = c(FALSE, TRUE, TRUE),
                                        stringsAsFactors = FALSE)
                ##
                markers <-  markers[markers$after_n_threads < n_threads_max, , drop = FALSE]
                ##
                if (nrow(markers) == 0) return(data.frame(device = character(), after_n_threads = numeric(), label = character()))
                ##
              # markers <-  stats::aggregate(change ~ device + after_n_threads, data = markers, FUN = function(changes) paste(changes, collapse = ";\n"))
                ##
                ## ---- 2026-10-06: per thread count, "SMT" first, then the cache sizes per thread on one line:
                fn_combine_changes_at_one_thread_count <-  function(rows) {

                        cache_changes <-  rows$change[rows$per_thread_cache]
                        ##
                        if (length(cache_changes)) {
                            cache_changes <-  paste0("per thread: ", paste(cache_changes, collapse = ", "))
                        }
                        ##
                        changes <-  c(rows$change[!rows$per_thread_cache], cache_changes)
                        ##
                        data.frame( device          = rows$device[1],
                                    after_n_threads = rows$after_n_threads[1],
                                    change          = paste(changes, collapse = ";\n"),
                                    stringsAsFactors = FALSE)

                }
                ##
                markers <-  do.call(what = rbind, args = lapply(X   = split(markers, markers$after_n_threads),
                                                                FUN = fn_combine_changes_at_one_thread_count))
                ##
                rownames(markers) <-  NULL
                ##
                markers$label <-  paste0("> ", markers$after_n_threads, " threads: ", markers$change)
                ##
                markers[order(markers$after_n_threads), c("device", "after_n_threads", "label"), drop = FALSE]

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

        ## ---- Both thread settings on ONE plot: colour = SMT (physical cores only vs. SMT threads),
        ##      line type = device (local-HPC solid, laptop dashed), as in R_fn_plot_ps1_N_chunks_ggplot_1 above:
        ##
        R_fn_plot_ps1_N_chunks_ggplot_SMT_combined <-  function( df_both,
                                                                 save_plot = TRUE,
                                                                 output_path,
                                                                 n_threads_for_HPC_no_SMT = 96,
                                                                 n_threads_for_Laptop_no_SMT = 8,
                                                                 n_threads_for_HPC_SMT = 180,
                                                                 n_threads_for_Laptop_SMT = 16,
                                                                 show_cache_lines = FALSE,
                                                                 bytes_per_row = paper1_bytes_per_row_BayesMVP,
                                                                 thresholds_file = NULL,
                                                                 panel_width_in = 7.2
        ) {

                ##
                SMT_levels <-  c( paste0("Without SMT (HPC: ", n_threads_for_HPC_no_SMT, " threads, laptop: ", n_threads_for_Laptop_no_SMT, " threads)"),
                                  paste0("With SMT (HPC: ",    n_threads_for_HPC_SMT,    " threads, laptop: ", n_threads_for_Laptop_SMT,    " threads)"))
                ##
                SMT_colours <-  stats::setNames(c("#0072B2", "#D55E00"), SMT_levels)
                ##
                ## (legend text in the paper's notation; the level strings above stay the keys of the colour scale)
                SMT_display_labels <-  c( bquote("Without SMT (local-HPC: " * {N[threads] == .(n_threads_for_HPC_no_SMT)} * "; laptop: " *
                                                 {N[threads] == .(n_threads_for_Laptop_no_SMT)} * ")"),
                                          bquote("With SMT (local-HPC: " * {N[threads] == .(n_threads_for_HPC_SMT)} * "; laptop: " *
                                                 {N[threads] == .(n_threads_for_Laptop_SMT)} * ")"))
                ##
                keep_rows <-  (df_both$device == "HPC"    & df_both$n_threads %in% c(n_threads_for_HPC_no_SMT,    n_threads_for_HPC_SMT)) |
                              (df_both$device == "Laptop" & df_both$n_threads %in% c(n_threads_for_Laptop_no_SMT, n_threads_for_Laptop_SMT))
                ##
                df_SMT <-  df_both[keep_rows, , drop = FALSE]
                ##
                df_SMT$Device <-  factor( ifelse(df_SMT$device == "HPC", "HPC", "Laptop"),
                                          levels = c("HPC", "Laptop"))
                ##
                df_SMT$SMT <-  factor( ifelse( df_SMT$n_threads %in% c(n_threads_for_HPC_SMT, n_threads_for_Laptop_SMT),
                                               SMT_levels[2],
                                               SMT_levels[1]),
                                       levels = SMT_levels)
                ##
                ## ---- Optional cache-capacity lines (show_cache_lines = TRUE): for each device and SMT setting, the chunk count at which
                ##      one chunk first fits in the L3, L2 or L1 cache per active thread, drawn within each panel's tested chunk range in
                ##      the colour and line type of the data line it belongs to (left of the L3 line, one chunk spills to RAM):
                ##
                cache_thresholds <-  NULL
                ##
                if (show_cache_lines) {

                        df_SMT$N_chunks_value <-  as.numeric(as.character(df_SMT$N_chunks))
                        ##
                        cache_cases <-  data.frame( Device    = c("HPC", "HPC", "Laptop", "Laptop"),
                                                    n_threads = c( n_threads_for_HPC_no_SMT, n_threads_for_HPC_SMT,
                                                                   n_threads_for_Laptop_no_SMT, n_threads_for_Laptop_SMT),
                                                    SMT_level = SMT_levels[c(1, 2, 1, 2)],
                                                    stringsAsFactors = FALSE)
                        ##
                        cache_thresholds <-  do.call(rbind, lapply(sort(unique(as.numeric(as.character(df_SMT$N)))), function(N_value) {

                                df_N <-  df_SMT[as.numeric(as.character(df_SMT$N)) == N_value, , drop = FALSE]
                                ##
                                thresholds_N <-  do.call(rbind, lapply(seq_len(nrow(cache_cases)), function(i) {

                                        ## only device/thread settings with data in this panel
                                        if (!any(df_N$Device == cache_cases$Device[i] & df_N$n_threads == cache_cases$n_threads[i])) return(NULL)
                                        ##
                                        thresholds_i <-  fn_paper1_cache_threshold_lines( N               = N_value,
                                                                                          chunk_range     = range(df_N$N_chunks_value),
                                                                                          device          = cache_cases$Device[i],
                                                                                          n_threads_total = cache_cases$n_threads[i],
                                                                                          bytes_per_row   = bytes_per_row)
                                        ##
                                        if (nrow(thresholds_i)) thresholds_i$SMT_level <-  cache_cases$SMT_level[i]
                                        ##
                                        thresholds_i

                                }))
                                ##
                                if (is.null(thresholds_N) || nrow(thresholds_N) == 0) return(NULL)
                                ##
                                thresholds_N$N_label <-  as.character(df_N$N_label[1])
                                ##
                                thresholds_N$label <-  paste0(thresholds_N$level, ifelse(thresholds_N$device == "HPC", " HPC", " laptop"))
                                ##
                                thresholds_N$label_row <-  fn_paper1_cache_label_rows( x              = thresholds_N$c_star,
                                                                                       labels         = thresholds_N$label,
                                                                                       chunk_range    = range(df_N$N_chunks_value),
                                                                                       panel_width_in = panel_width_in)
                                ##
                                thresholds_N

                        }))
                        ##
                        if (!is.null(cache_thresholds) && nrow(cache_thresholds)) {

                            cache_thresholds <-  cache_thresholds[order(cache_thresholds$N, cache_thresholds$c_star), , drop = FALSE]
                            ##
                            cache_thresholds$N_label <-  factor(cache_thresholds$N_label, levels = levels(factor(df_SMT$N_label)))
                            cache_thresholds$Device  <-  factor(cache_thresholds$device, levels = c("HPC", "Laptop"))
                            cache_thresholds$SMT     <-  factor(cache_thresholds$SMT_level, levels = SMT_levels)
                            cache_thresholds$vjust   <-  1.3 + 1.35 * cache_thresholds$label_row
                            ##
                            if (!is.null(thresholds_file)) {

                                utils::write.csv( x = cache_thresholds[, c("N", "device", "n_threads_total", "SMT_level", "level",
                                                                           "cache_bytes_per_thread", "c_star", "label", "label_row")],
                                                  file = thresholds_file,
                                                  row.names = FALSE)
                                ##
                                message(paste0("\033[36m", "Cache-capacity thresholds written to: ", thresholds_file, "\033[0m"))

                            }

                        }

                }
                ##
                make_half <-  function( N_range,
                                        hide_legend) {

                        df_f <-  df_SMT[as.numeric(as.character(df_SMT$N)) %in% N_range, , drop = FALSE]
                        ##
                        if (nrow(df_f) == 0) {

                            return(empty_panel(paste0( "No matching results supplied for N = ",
                                                       paste(N_range, collapse = ", "))))

                        }
                        ##
                        cache_thresholds_half <-  if (is.null(cache_thresholds)) NULL else
                            cache_thresholds[cache_thresholds$N %in% N_range, , drop = FALSE]
                        ##
                        cache_line_layers <-  fn_paper1_cache_line_layers( thresholds      = cache_thresholds_half,
                                                                           colour_column   = "SMT",
                                                                           linetype_column = "Device")
                        ##
                        p <-  ggplot2::ggplot(df_f, ggplot2::aes( x = N_chunks,
                                                                  y = time_avg,
                                                                  colour = SMT,
                                                                  linetype = Device,
                                                                  group = interaction(Device, SMT))) +
                            cache_line_layers +
                            ggplot2::geom_point(size = 5) +
                            ggplot2::geom_errorbar( linewidth = 1,
                                                    width = 0.02,
                                                    ggplot2::aes( ymin = time_avg - time_SD,        ## +/- 1 SD (caption says "standard deviation")
                                                                  ymax = time_avg + time_SD)) +
                            ggplot2::geom_line(linewidth = 1) +
                            ## ggplot2::scale_colour_manual(values = SMT_colours, limits = SMT_levels, drop = FALSE) +
                            ggplot2::scale_colour_manual(values = SMT_colours, limits = SMT_levels, drop = FALSE,
                                                         labels = do.call(what = expression, args = SMT_display_labels)) +
                            ggplot2::scale_linetype_manual(values = c(HPC = "solid", Laptop = "22"), drop = FALSE, labels = c(HPC = "local-HPC", Laptop = "laptop")) +
                            ggplot2::theme_bw(base_size = 24) +
                            ggplot2::ylab("Time (sec.)") + ggplot2::xlab(expression(N[chunks])) +
                            ggplot2::labs(colour = NULL, linetype = "Device") +
                            ggplot2::facet_wrap(~ N_label, scales = "free") +
                            ggplot2::guides( colour = ggplot2::guide_legend(order = 1, ncol = 1),
                                             linetype = ggplot2::guide_legend(order = 2, override.aes = list(colour = "black"))) +
                            ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                                            legend.position = if (hide_legend) "none" else "bottom",
                                            legend.box = "vertical",
                                            legend.key.width = ggplot2::unit(3, "lines"))
                        ##
                        ## ---- With cache lines: log10 chunk axis (breaks at each panel's tested chunk counts) and head room for the labels:
                        if (show_cache_lines) {

                            n_label_rows <-  if (is.null(cache_thresholds_half) || nrow(cache_thresholds_half) == 0) 0 else
                                max(cache_thresholds_half$label_row) + 1
                            ##
                            p <-  p +
                                ggplot2::aes(x = N_chunks_value) +
                                ggplot2::scale_x_log10(breaks = fn_paper1_log10_chunk_breaks(split(df_f$N_chunks_value, df_f$N))) +
                                ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.05 + 0.07 * n_label_rows))) +
                                ggplot2::xlab(expression(N[chunks]~"(log scale)"))

                        }
                        ##
                        p

                }
                ##
                has_bottom <-  any(as.numeric(as.character(df_SMT$N)) %in% c(10000, 50000))
                ##
                p_top    <-  make_half(c(500, 2500),     hide_legend = has_bottom) + ggplot2::xlab(" ")
                ##
                p_bottom <-  make_half(c(10000, 50000),  hide_legend = FALSE)
                ##
                combined <-  p_top + p_bottom + patchwork::plot_layout(ncol = 1)
                ##
                print(combined)
                ##
                if (save_plot) {

                    ggplot2::ggsave( file.path( output_path,
                                                if (show_cache_lines) "Figure_N_chunks_pilot_study_plot_1_n_threads_SMT_vs_no_SMT_cache_lines.png" else
                                                    "Figure_N_chunks_pilot_study_plot_1_n_threads_SMT_vs_no_SMT.png"),
                                     combined,
                                     width = 16,
                                     height = 17,
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

        ## ---- Efficiency for both devices on ONE set of four panels: colour = thread level, line type = device
        ##      (local-HPC solid, laptop dashed); log y axis, since the HPC values are about 15 times the laptop's:
        ##
        R_fn_plot_ps1_efficiency_combined <-  function( df_both,
                                                        threads_HPC    = c(64, t_phys_HPC, t_smt_HPC),
                                                        threads_Laptop = c(4,  t_phys_Laptop, t_smt_Laptop),
                                                        save_plot = TRUE,
                                                        output_path,
                                                        show_cache_lines = FALSE,
                                                        bytes_per_row = paper1_bytes_per_row_BayesMVP,
                                                        thresholds_file = NULL,
                                                        panel_width_in = 7.2
        ) {

                ##
                level_labels <-  c( paste0("HPC: ", threads_HPC[1], " threads, laptop: ", threads_Laptop[1], " threads"),
                                    paste0("HPC: ", threads_HPC[2], " threads, laptop: ", threads_Laptop[2], " threads (without SMT)"),
                                    paste0("HPC: ", threads_HPC[3], " threads, laptop: ", threads_Laptop[3], " threads (with SMT)"))
                ##
                level_colours <-  stats::setNames(c("#009E73", "#0072B2", "#D55E00"), level_labels)
                ##
                ## (legend text in the paper's notation; the level strings above stay the keys of the colour scale)
                level_display_labels <-  lapply(X = 1:3, FUN = function(i)
                                             bquote("local-HPC: " * {N[threads] == .(threads_HPC[i])} * "; laptop: " *
                                                    {N[threads] == .(threads_Laptop[i])} * .(c("", " (without SMT)", " (with SMT)")[i])))
                ##
                keep_rows <-  (df_both$device == "HPC"    & df_both$n_threads %in% threads_HPC) |
                              (df_both$device == "Laptop" & df_both$n_threads %in% threads_Laptop)
                ##
                df_eff <-  df_both[keep_rows, , drop = FALSE]
                ##
                if (nrow(df_eff) == 0) return(empty_panel("No matching efficiency results supplied"))
                ##
                df_eff$Device <-  factor( ifelse(df_eff$device == "HPC", "HPC", "Laptop"),
                                          levels = c("HPC", "Laptop"))
                ##
                df_eff$Threads <-  factor( level_labels[ifelse( df_eff$device == "HPC",
                                                                match(df_eff$n_threads, threads_HPC),
                                                                match(df_eff$n_threads, threads_Laptop))],
                                           levels = level_labels)
                ##
                ## ---- Optional cache-capacity lines (show_cache_lines = TRUE): for each device and thread level, the chunk count at which
                ##      one chunk first fits in the L3, L2 or L1 cache per active thread, within each panel's tested chunk range. Thread
                ##      levels of one device with identical thresholds (e.g. 64 and 96 HPC threads, both without SMT) share one line, in
                ##      the colour of the highest of those levels and labelled with all of their thread counts:
                ##
                cache_thresholds <-  NULL
                ##
                y_expand <-  ggplot2::waiver()
                ##
                if (show_cache_lines) {

                        df_eff$N_chunks_value <-  as.numeric(as.character(df_eff$N_chunks))
                        ##
                        cache_thresholds <-  do.call(rbind, lapply(sort(unique(as.numeric(as.character(df_eff$N)))), function(N_value) {

                                df_N <-  df_eff[as.numeric(as.character(df_eff$N)) == N_value, , drop = FALSE]
                                ##
                                thresholds_N <-  do.call(rbind, lapply(c("HPC", "Laptop"), function(dev) {

                                        threads_dev <-  if (dev == "HPC") threads_HPC else threads_Laptop
                                        ##
                                        threads_present <-  threads_dev[threads_dev %in% df_N$n_threads[df_N$Device == dev]]
                                        ##
                                        if (length(threads_present) == 0) return(NULL)
                                        ##
                                        thresholds_dev <-  fn_paper1_cache_threshold_lines( N               = N_value,
                                                                                            chunk_range     = range(df_N$N_chunks_value),
                                                                                            device          = dev,
                                                                                            n_threads_total = threads_present,
                                                                                            bytes_per_row   = bytes_per_row)
                                        ##
                                        if (nrow(thresholds_dev) == 0) return(NULL)
                                        ##
                                        ## one line per distinct threshold of this device:
                                        thresholds_dev$key <-  paste(thresholds_dev$level, signif(thresholds_dev$c_star, 10))
                                        ##
                                        do.call(rbind, lapply(split(thresholds_dev, thresholds_dev$key), function(group) {

                                                group <-  group[order(group$n_threads_total), , drop = FALSE]
                                                ##
                                                highest <-  max(group$n_threads_total)
                                                ##
                                                data.frame( N                      = N_value,
                                                            device                 = dev,
                                                            n_threads_total        = paste(group$n_threads_total, collapse = "/"),
                                                            SMT                    = group$SMT[nrow(group)],
                                                            level                  = group$level[1],
                                                            cache_bytes_per_thread = group$cache_bytes_per_thread[1],
                                                            c_star                 = group$c_star[1],
                                                            Threads_level          = level_labels[match(highest, threads_dev)],
                                                            label                  = paste0( group$level[1],
                                                                                             ## ifelse(dev == "HPC", " HPC ", " laptop "),
                                                                                             ifelse(dev == "HPC", " local-HPC ", " laptop "),
                                                                                             paste(group$n_threads_total, collapse = "/")),
                                                            stringsAsFactors       = FALSE)

                                        }))

                                }))
                                ##
                                if (is.null(thresholds_N) || nrow(thresholds_N) == 0) return(NULL)
                                ##
                                thresholds_N$N_label <-  as.character(df_N$N_label[1])
                                ##
                                thresholds_N$label_row <-  fn_paper1_cache_label_rows( x              = thresholds_N$c_star,
                                                                                       labels         = thresholds_N$label,
                                                                                       chunk_range    = range(df_N$N_chunks_value),
                                                                                       panel_width_in = panel_width_in)
                                ##
                                thresholds_N

                        }))
                        ##
                        if (!is.null(cache_thresholds) && nrow(cache_thresholds)) {

                            cache_thresholds <-  cache_thresholds[order(cache_thresholds$N, cache_thresholds$c_star), , drop = FALSE]
                            ##
                            cache_thresholds$N_label <-  factor(cache_thresholds$N_label, levels = levels(factor(df_eff$N_label)))
                            cache_thresholds$Device  <-  factor(cache_thresholds$device, levels = c("HPC", "Laptop"))
                            cache_thresholds$Threads <-  factor(cache_thresholds$Threads_level, levels = level_labels)
                            cache_thresholds$vjust   <-  1.3 + 1.35 * cache_thresholds$label_row
                            ##
                            y_expand <-  ggplot2::expansion(mult = c(0.05, 0.05 + 0.07 * (max(cache_thresholds$label_row) + 1)))
                            ##
                            if (!is.null(thresholds_file)) {

                                utils::write.csv( x = cache_thresholds[, c("N", "device", "n_threads_total", "SMT", "level", "cache_bytes_per_thread",
                                                                           "c_star", "Threads_level", "label", "label_row")],
                                                  file = thresholds_file,
                                                  row.names = FALSE)
                                ##
                                message(paste0("\033[36m", "Cache-capacity thresholds written to: ", thresholds_file, "\033[0m"))

                            }

                        }

                }
                ##
                cache_line_layers <-  fn_paper1_cache_line_layers( thresholds      = cache_thresholds,
                                                                   colour_column   = "Threads",
                                                                   linetype_column = "Device")
                ##
                combined <-  ggplot2::ggplot(df_eff, ggplot2::aes( x = N_chunks,
                                                                   y = Efficiency,
                                                                   colour = Threads,
                                                                   linetype = Device,
                                                                   group = interaction(Device, Threads))) +
                    cache_line_layers +
                    ggplot2::geom_point(size = 5) +
                    ggplot2::geom_line(linewidth = 1) +
                    ggplot2::scale_y_log10(expand = y_expand) +
                    ## ggplot2::scale_colour_manual(values = level_colours, limits = level_labels, drop = FALSE) +
                    ggplot2::scale_colour_manual(values = level_colours, limits = level_labels, drop = FALSE,
                                                 labels = do.call(what = expression, args = level_display_labels)) +
                    ggplot2::scale_linetype_manual(values = c(HPC = "solid", Laptop = "22"), drop = FALSE, labels = c(HPC = "local-HPC", Laptop = "laptop")) +
                    ggplot2::theme_bw(base_size = 24) +
                    ggplot2::ylab(expression(Efficiency~(N[threads] / time)~"(log scale)")) +
                    ggplot2::xlab(expression(N[chunks])) +
                    ggplot2::labs(colour = NULL, linetype = "Device") +
                    ggplot2::facet_wrap(~ N_label, scales = "free") +
                    ggplot2::guides( colour = ggplot2::guide_legend(order = 1, ncol = 1),
                                     linetype = ggplot2::guide_legend(order = 2, override.aes = list(colour = "black"))) +
                    ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                                    legend.position = "bottom",
                                    legend.box = "vertical",
                                    legend.key.width = ggplot2::unit(3, "lines"))
                ##
                ## ---- With cache lines: log10 chunk axis, with breaks at each panel's tested chunk counts:
                if (show_cache_lines) {

                    combined <-  combined +
                        ggplot2::aes(x = N_chunks_value) +
                        ggplot2::scale_x_log10(breaks = fn_paper1_log10_chunk_breaks(split(df_eff$N_chunks_value, df_eff$N))) +
                        ggplot2::xlab(expression(N[chunks]~"(log scale)"))

                }
                ##
                print(combined)
                ##
                if (save_plot) {

                    ggplot2::ggsave( file.path( output_path,
                                                if (show_cache_lines) "Figure_N_chunks_pilot_study_plot_3_both_devices_cache_lines.png" else
                                                    "Figure_N_chunks_pilot_study_plot_3_both_devices.png"),
                                     combined,
                                     width = 16,
                                     height = 17,
                                     dpi = 100)

                }
                ##
                invisible(combined)

        }

        ## ---- Paper 1 WCP chunk-search views -------------------------------------------------------------------------------------------
        fn_plot_paper1_WCP_chunk_search <-  function( chunk_search,
                                                       best_chunks,
                                                       output_path,
                                                       file_prefix,
                                                       show_cache_lines = FALSE,
                                                       bytes_per_row = paper1_bytes_per_row_BayesMVP,
                                                       thresholds_file = NULL,
                                                       panel_width_in = 4.5
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
                ## plot_title <-  paste0(chunk_search$device[1], ": ", chunk_search$Algorithm_label[1],
                ##                        ", N = ", fn_paper1_format_number_commas_from_10000(chunk_search$N[1]))
                ##
                ## The device as named in the paper, and the two compared configurations by their model names
                ## (WCP-only, i.e. N_chunks = N_threads/chain, vs. chunking + WCP):
                device_label <-  if (identical(chunk_search$device[1], "HPC")) "local-HPC" else chunk_search$device[1]
                WCP_only_model_name     <-  chunk_search$algorithm[1]
                WCP_chunking_model_name <-  paste0(chunk_search$algorithm[1], "_chunking")
                ##
                plot_title <-  paste0(device_label, ", N = ", fn_paper1_format_number_commas_from_10000(chunk_search$N[1]), ": ",
                                       WCP_only_model_name, " vs. ", WCP_chunking_model_name)
                ##
                ## Keep all measured chunks; the outline marks each fixed chain/WCP allocation's optimum.
                ## wcp_axis_label <-  if ("execution_backend" %in% names(x = chunk_search) &&
                ##                         any(chunk_search$execution_backend == "NicoStan_BridgeStan")) {
                ##
                ##     "WCP budget / chain (shared pool)"
                ##
                ## } else "WCP threads / chain"
                ## N_threads/chain, as in the papers (the WCP threads of each chain), for both implementations:
                wcp_axis_label <-  quote(N["threads/chain"])
                ##
                ## ---- Optional cache-capacity lines (show_cache_lines = TRUE): in each chain-count panel the active threads are
                ##      chains x WCP threads per chain. WCP counts with identical thresholds share one grey line, labelled with the cache
                ##      level and the WCP counts it applies to; a line is drawn only where it lies inside the tested chunk range of at
                ##      least one of those WCP series (left of an L3 line, one chunk exceeds the L3 cache per active thread):
                ##
                cache_thresholds <-  NULL
                ##
                if (show_cache_lines) {

                        fn_WCP_range_text <-  function( WCP_group,
                                                        WCP_panel) {

                                if (length(WCP_group) == 1)                  return(paste0("WCP ", WCP_group))
                                if (length(WCP_group) == length(WCP_panel)) return("all WCP")
                                ## (WCP counts that are not neighbours in the panel are listed, not given as a range)
                                if (any(diff(match(sort(WCP_group), WCP_panel)) != 1)) return(paste0("WCP ", paste(sort(WCP_group), collapse = "/")))
                                if (min(WCP_group) == min(WCP_panel))       return(paste0("WCP \u2264 ", max(WCP_group)))
                                if (max(WCP_group) == max(WCP_panel))       return(paste0("WCP \u2265 ", min(WCP_group)))
                                ##
                                paste0("WCP ", min(WCP_group), "-", max(WCP_group))

                        }
                        ##
                        ## The same label in the paper's notation (plotmath, drawn with parse = TRUE), e.g. "L3 (N_threads/chain <= 16)":
                        fn_WCP_range_expression <-  function( level,
                                                              WCP_group,
                                                              WCP_panel) {

                                N_threads_per_chain <-  'N["threads/chain"]'
                                ##
                                inner <-  if (length(WCP_group) == 1) paste0(N_threads_per_chain, "==", WCP_group) else
                                          if (length(WCP_group) == length(WCP_panel)) paste0('"all"~', N_threads_per_chain) else
                                          if (any(diff(match(sort(WCP_group), WCP_panel)) != 1)) {
                                              paste0(N_threads_per_chain, '=="', paste(sort(WCP_group), collapse = "/"), '"')
                                          } else
                                          if (min(WCP_group) == min(WCP_panel)) paste0(N_threads_per_chain, "<=", max(WCP_group)) else
                                          if (max(WCP_group) == max(WCP_panel)) paste0(N_threads_per_chain, ">=", min(WCP_group)) else
                                          paste0(min(WCP_group), "<=", N_threads_per_chain, "<=", max(WCP_group))
                                ##
                                paste0(level, "~(", inner, ")")

                        }
                        ##
                        cache_thresholds <-  do.call(rbind, lapply(sort(unique(chunk_search$n_chains)), function(n_chains_value) {

                                panel <-  chunk_search[chunk_search$n_chains == n_chains_value, , drop = FALSE]
                                ##
                                WCP_panel <-  sort(unique(panel$threads_per_chain))
                                ##
                                per_WCP <-  do.call(rbind, lapply(WCP_panel, function(WCP) {

                                        series_chunks <-  panel$num_chunks[panel$threads_per_chain == WCP]
                                        ##
                                        thresholds_WCP <-  fn_paper1_cache_threshold_lines( N               = chunk_search$N[1],
                                                                                            chunk_range     = c(0, Inf),
                                                                                            device          = chunk_search$device[1],
                                                                                            n_threads_total = n_chains_value * WCP,
                                                                                            bytes_per_row   = bytes_per_row)
                                        ##
                                        thresholds_WCP$threads_per_chain <-  WCP
                                        ##
                                        thresholds_WCP$inside_series <-  thresholds_WCP$c_star > min(series_chunks) & thresholds_WCP$c_star < max(series_chunks)
                                        ##
                                        thresholds_WCP

                                }))
                                ##
                                per_WCP$key <-  paste(per_WCP$level, signif(per_WCP$c_star, 10))
                                ##
                                thresholds_panel <-  do.call(rbind, lapply(split(per_WCP, per_WCP$key), function(group) {

                                        if (!any(group$inside_series)) return(NULL)
                                        ##
                                        ## the line and its label apply only to the WCP series whose tested chunk range contains it:
                                        group <-  group[group$inside_series, , drop = FALSE]
                                        ##
                                        data.frame( N                      = group$N[1],
                                                    device                 = group$device[1],
                                                    n_chains               = n_chains_value,
                                                    threads_per_chain      = paste(group$threads_per_chain, collapse = "/"),
                                                    n_threads_total        = paste(group$n_threads_total, collapse = "/"),
                                                    SMT                    = paste(unique(group$SMT), collapse = "/"),
                                                    level                  = group$level[1],
                                                    cache_bytes_per_thread = group$cache_bytes_per_thread[1],
                                                    c_star                 = group$c_star[1],
                                                    label                  = paste0(group$level[1], " (", fn_WCP_range_text(group$threads_per_chain, WCP_panel), ")"),
                                                    label_expression       = fn_WCP_range_expression( level     = group$level[1],
                                                                                                      WCP_group = group$threads_per_chain,
                                                                                                      WCP_panel = WCP_panel),
                                                    stringsAsFactors       = FALSE)

                                }))
                                ##
                                if (is.null(thresholds_panel) || nrow(thresholds_panel) == 0) return(NULL)
                                ##
                                ## labels near either end of the chunk axis move inwards so they are not cut off by the panel edge:
                                ## (the drawn labels say N_threads/chain instead of WCP, so their width is estimated from that text)
                                label_width_text <-  sub(pattern = "WCP", replacement = "N threads/chain", x = thresholds_panel$label)
                                ##
                                thresholds_panel$label_x <-  fn_paper1_cache_label_x( x              = thresholds_panel$c_star,
                                                                                      ## labels         = thresholds_panel$label,
                                                                                      labels         = label_width_text,
                                                                                      chunk_range    = range(chunk_search$num_chunks),
                                                                                      panel_width_in = panel_width_in,
                                                                                      label_size     = 3.6)
                                ##
                                thresholds_panel$label_row <-  fn_paper1_cache_label_rows( x              = thresholds_panel$c_star,
                                                                                           ## labels         = thresholds_panel$label,
                                                                                           labels         = label_width_text,
                                                                                           chunk_range    = range(chunk_search$num_chunks),
                                                                                           panel_width_in = panel_width_in,
                                                                                           label_size     = 3.6,
                                                                                           label_x        = thresholds_panel$label_x)
                                ##
                                thresholds_panel

                        }))
                        ##
                        if (!is.null(cache_thresholds) && nrow(cache_thresholds)) {

                            cache_thresholds <-  cache_thresholds[order(cache_thresholds$n_chains, cache_thresholds$c_star), , drop = FALSE]
                            ##
                            cache_thresholds$vjust <-  1.3 + 1.35 * cache_thresholds$label_row
                            ##
                            if (!is.null(thresholds_file)) {

                                utils::write.csv( x = cache_thresholds[, c("N", "device", "n_chains", "threads_per_chain", "n_threads_total", "SMT", "level",
                                                                           "cache_bytes_per_thread", "c_star", "label", "label_row")],
                                                  file = thresholds_file,
                                                  row.names = FALSE)
                                ##
                                message(paste0("\033[36m", "Cache-capacity thresholds written to: ", thresholds_file, "\033[0m"))

                            }

                        } else if (!is.null(thresholds_file)) {

                            ## no threshold within the tested chunk range: an empty table records this
                            utils::write.csv( x = data.frame( N = numeric(), device = character(), n_chains = numeric(), threads_per_chain = character(),
                                                              n_threads_total = character(), SMT = character(), level = character(),
                                                              cache_bytes_per_thread = numeric(), c_star = numeric(), label = character(), label_row = numeric()),
                                              file = thresholds_file,
                                              row.names = FALSE)

                        }

                }
                ##
                cache_line_layers <-  fn_paper1_cache_line_layers(thresholds = cache_thresholds, label_size = 3.6)
                ##
                ## chunk_search_plot <-  ggplot( data = chunk_search,
                ##                                mapping = aes(x = chunk_label, y = chain_rate, colour = WCP_label, group = WCP_label)) +
                ##     cache_line_layers +
                ##     geom_line(linewidth = 0.8) + geom_point(size = 3) +
                ##     geom_point(data = chunk_search[chunk_search$selected_best_chunks, , drop = FALSE],
                ##                shape = 21, fill = "white", size = 5, stroke = 1.2) +
                ##     theme_bw(base_size = 20) +
                ##     theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
                ##     labs( x = "Chunks", y = "Within-method efficiency (chains / time)", colour = wcp_axis_label,
                ##            title = plot_title, subtitle = "All measured chunks; outlined points maximise throughput at each fixed chain/WCP count") +
                ##     facet_wrap(facets = ~ n_chains, scales = "free_y", labeller = label_both)
                ##
                ## ---- WCP-only (N_chunks = N_threads/chain) is marked on every series with a filled black triangle, and the fastest
                ##      N_chunks of each N_chains x N_threads/chain allocation (chunking + WCP) with an open circle; the legend names both:
                WCP_only_rows <-  chunk_search$num_chunks == chunk_search$threads_per_chain
                ##
                missing_WCP_only <-  setdiff( x = unique(paste(chunk_search$n_chains, chunk_search$threads_per_chain)),
                                              y = paste(chunk_search$n_chains, chunk_search$threads_per_chain)[WCP_only_rows])
                if (length(missing_WCP_only)) message(paste0("\033[31m", "No WCP-only point (N_chunks = N_threads/chain) for N_chains x N_threads/chain: ",
                                                             paste(sub(" ", " x ", missing_WCP_only), collapse = ", "), "\033[0m"))
                ##
                ## (the multiplication sign as a character: plotmath's %*% is drawn as a centred dot by this device)
                marker_labels <-  do.call(expression, list( bquote("fastest " * N[chunks] * " of each " * N[chains] * " \u00D7 " * N["threads/chain"] *
                                                                   " (" * .(WCP_chunking_model_name) * ")"),
                                                            bquote(.(WCP_only_model_name) * " (WCP-only: " * N[chunks] == N["threads/chain"] * ")")))
                ##
                chunk_search_plot <-  ggplot( data = chunk_search,
                                               mapping = aes(x = chunk_label, y = chain_rate, colour = WCP_label, group = WCP_label)) +
                    cache_line_layers +
                    geom_line(linewidth = 0.8) + geom_point(size = 3) +
                    geom_point( data = chunk_search[chunk_search$selected_best_chunks, , drop = FALSE],
                                mapping = aes(shape = "fastest"), fill = "white", size = 5, stroke = 1.2) +
                    ## geom_point( data = chunk_search[WCP_only_rows, , drop = FALSE],
                    ##             mapping = aes(shape = "WCP_only"), colour = "black", fill = "black", size = 3.4) +
                    ## (each WCP-only triangle is filled with the colour of its N_threads/chain series, with a black outline)
                    geom_point( data = chunk_search[WCP_only_rows, , drop = FALSE],
                                mapping = aes(shape = "WCP_only", fill = WCP_label), colour = "black", size = 3.8, stroke = 1) +
                    scale_fill_hue(guide = "none") +
                    scale_shape_manual( name   = NULL,
                                        values = c(fastest = 21, WCP_only = 24),
                                        breaks = c("fastest", "WCP_only"),
                                        labels = marker_labels) +
                    guides( colour = guide_legend(order = 1, nrow = 2),
                            shape  = guide_legend(order = 2, ncol = 1,
                                                  override.aes = list(colour = "black", fill = c("white", "grey55"), size = c(5, 3.8)))) +
                    theme_bw(base_size = 20) +
                    theme( legend.position = "bottom", legend.box = "vertical", axis.text.x = element_text(angle = 45, hjust = 1)) +
                    labs( x = quote(N[chunks]), y = "Within-method efficiency (chains/second)", colour = wcp_axis_label,
                           title = plot_title) +
                    facet_wrap(facets = ~ n_chains, scales = "free_y", labeller = label_bquote(cols = N[chains] == .(n_chains)))
                ##
                ## ---- With cache lines: log10 chunk axis with breaks at the tested chunk counts (the chunk counts of the smallest WCP
                ##      series first; any other tested count only where its tick label does not overlap a neighbouring one), and head
                ##      room for the line labels:
                if (show_cache_lines) {

                    chunk_breaks <-  sort(unique(chunk_search$num_chunks))
                    ##
                    ## (the chunk counts of the outlined optima are always labelled)
                    chunk_labels <-  fn_paper1_log10_chunk_labels( breaks           = chunk_breaks,
                                                                   priority_chunks  = chunk_search$num_chunks[chunk_search$selected_best_chunks],
                                                                   secondary_chunks = chunk_search$num_chunks[chunk_search$threads_per_chain ==
                                                                                                               min(chunk_search$threads_per_chain)],
                                                                   min_gap_log10    = 0.22 * 1.1 * diff(log10(range(chunk_breaks))) / panel_width_in)
                    ##
                    n_label_rows <-  if (is.null(cache_thresholds) || nrow(cache_thresholds) == 0) 0 else max(cache_thresholds$label_row) + 1
                    ##
                    levels_drawn <-  intersect(c("L3", "L2", "L1"), cache_thresholds$level)
                    ##
                    ## ---- 2026-10-04: head room for the line labels in EACH N_chains panel, in proportion to
                    ##      that panel's own rows of labels (an invisible point above its data), instead of one
                    ##      expansion for every panel: with one L3 line per N_threads/chain below full load, a
                    ##      panel can have nine rows of labels. Label row r ends (1.3 + 1.35 r) label heights
                    ##      (~0.21 inches) below the panel top; the panel height is estimated from the figure
                    ##      height used in ggsave below.
                    headroom_layer <-  NULL
                    ##
                    if (n_label_rows > 0) {

                        n_panel_rows    <-  ceiling(length(unique(chunk_search$n_chains)) / 3)
                        figure_height   <-  4 + 4 * n_panel_rows + min(2, 0.5 * n_label_rows)
                        panel_height_in <-  (figure_height - 5.15 - 0.45 * n_panel_rows) / n_panel_rows
                        ##
                        panel_ids   <-  as.character(sort(unique(chunk_search$n_chains)))
                        panel_rows  <-  tapply(cache_thresholds$label_row, cache_thresholds$n_chains, max) + 1
                        rows_p      <-  ifelse(panel_ids %in% names(panel_rows), panel_rows[panel_ids], 0)
                        label_share <-  ifelse( rows_p > 0,
                                                (1.3 + 1.35 * (rows_p - 1)) * 0.21 / panel_height_in + 0.05,
                                                0)
                        label_share <-  pmin(0.8, label_share)
                        ##
                        y_max <-  tapply(chunk_search$chain_rate, chunk_search$n_chains, max)[panel_ids]
                        y_min <-  tapply(chunk_search$chain_rate, chunk_search$n_chains, min)[panel_ids]
                        ##
                        headroom <-  data.frame( n_chains   = as.numeric(panel_ids),
                                                 num_chunks = min(chunk_search$num_chunks),
                                                 chain_rate = as.numeric(y_max + (y_max - y_min) *
                                                                         1.05 * label_share / (1 - label_share)))
                        ##
                        headroom_layer <-  ggplot2::geom_blank( data = headroom,
                                                                mapping = ggplot2::aes(x = num_chunks,
                                                                                       y = chain_rate),
                                                                inherit.aes = FALSE)

                    }
                    ##
                    chunk_search_plot <-  chunk_search_plot +
                        headroom_layer +
                        ggplot2::aes(x = num_chunks) +
                        ggplot2::scale_x_log10( breaks = chunk_breaks[nzchar(chunk_labels)],
                                                labels = chunk_labels[nzchar(chunk_labels)]) +
## (the shared head room for the line labels before 2026-10-04:)
# ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.05 + 0.08 * n_label_rows))) +
                        ## (2026-10-04: the per-panel head room above replaces the shared expansion)
                        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
                        ggplot2::theme( axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
                                        panel.grid.minor.x = ggplot2::element_blank()) +
                        ## ggplot2::labs( x = "Chunks (log scale)",
                        ggplot2::labs( x = quote(N[chunks] ~ "(log scale)"),
                                       subtitle = if (length(levels_drawn) == 0) {

                                           ## (when every series' L2 threshold lies left of its smallest tested chunk count, every tested
                                           ##  chunk already fits in the L2 cache per active thread)
                                           all_chunks_fit_L2 <-  all(vapply(split(chunk_search, list(chunk_search$n_chains, chunk_search$threads_per_chain), drop = TRUE),
                                                                            function(series) {

                                                                                    L2_per_thread <-  fn_paper1_cache_per_thread( device          = series$device[1],
                                                                                                                                  n_threads_total = series$n_chains[1] * series$threads_per_chain[1])$L2
                                                                                    ##
                                                                                    bytes_per_row * series$N[1] / L2_per_thread <= min(series$num_chunks)

                                                                            }, logical(1)))
                                           ##
                                           ## paste0( chunk_search_plot$labels$subtitle,
                                           ##         "\nNo cache-capacity threshold (L3/L2/L1 per active thread) lies within the tested chunk range",
                                           ##         if (all_chunks_fit_L2) "\n(every tested chunk already fits in the L2 cache per active thread)" else "")
                                           paste0( "No cache-capacity threshold (L3/L2/L1 per active thread) lies within the tested range",
                                                   if (all_chunks_fit_L2) "\n(every tested chunk already fits in the L2 cache per active thread)" else "")

                                       } else {

                                           ## paste0( chunk_search_plot$labels$subtitle,
                                           ##         "\nGrey lines: chunk count at which one chunk first fits in the ",
                                           ##         paste(levels_drawn, collapse = "/"), " cache per active thread")
                                           bquote("Grey lines: " * N[chunks] * " at which one chunk first fits in the " *
                                                  .(paste(levels_drawn, collapse = "/")) * " cache per active thread")

                                       })

                }
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
                ##
                ## ---- With cache lines, only the chunk-search view is written (under its own file name; the view without lines is kept):
                if (show_cache_lines) {

                    plot_paths <-  file.path(output_path, paste0(file_prefix, "_chunk_search_cache_lines.png"))
                    ##
                    ## half an inch taller per row of line labels above the data (at most two inches):
                    ggplot2::ggsave(filename = plot_paths, plot = chunk_search_plot, width = 16, height = plot_height + min(2, 0.5 * n_label_rows), dpi = 150)
                    ##
                    return(invisible(x = plot_paths))

                }
                ##
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
                                                         highlight_best = FALSE,
                                                         legend_nrow = 2
        ) {

                ## legend_nrow: rows of the shared legend (2 by default; the five-arm Stan variants figure uses 5, i.e. one column,
                ## so that its long labels are not cut off at the figure edge).
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
                        guides(colour = guide_legend(title = NULL, nrow = legend_nrow)) +
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
                                            values_file = NULL,
                                            show_markers = FALSE,
                                            markers_file = NULL
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
                ## ---- the legend names each configuration by its model name, as in the text (Mplus with its iteration mode):
                ps2_model_names <-  c( "BayesMVP (1 chunk)"                          = "MD_BayesMVP",
                                       "BayesMVP + chunking"                         = "MD_BayesMVP_chunking",
                                       "BayesMVP + chunking + WCP"                   = "MD_BayesMVP_WCP_chunking",
                                     # "Mplus"                                       = "Mplus (BITERATIONS)",
                                     # "Mplus + WCP"                                 = "Mplus + WCP (BITERATIONS)",
                                       ## ---- 2026-10-06: the Mplus entries are named by their model
                                       ##      names too; the iteration mode (BITERATIONS) is too much
                                       ##      detail for a figure legend:
                                       "Mplus"                                       = "Mplus_standard",
                                       "Mplus + WCP"                                 = "Mplus_WCP",
                                       "Stan model (NicoStan)"                       = "AD_Stan",
                                       "Stan model (NicoStan) + tape chunking"       = "AD_Stan_tape_chunked",
                                       "Stan model (NicoStan) + tape chunking + WCP" = "AD_Stan_WCP_chunking")
                ps2_df$Algorithm_label <-  ifelse(as.character(ps2_df$Algorithm_label) %in% names(ps2_model_names),
                                                  unname(ps2_model_names[as.character(ps2_df$Algorithm_label)]),
                                                  as.character(ps2_df$Algorithm_label))
                ## (the plotted data were built from ps2_df above, so they are relabelled in the same way)
                df_plot$Algorithm_label <-  ifelse(as.character(df_plot$Algorithm_label) %in% names(ps2_model_names),
                                                   unname(ps2_model_names[as.character(df_plot$Algorithm_label)]),
                                                   as.character(df_plot$Algorithm_label))
                ##
                colour_scale <-  shared_colour_scale(ps2_df$Algorithm_label)
                ##
                plot_list <-  list()
                ##
                thread_markers_all <-  NULL
                ##
                for (dev in c("HPC", "Laptop")) {

                    ##
                    df_dev <-  df_plot %>% filter(device == dev)
                    ##
                    if (nrow(df_dev) == 0) next
                    ##
                    ## ---- 2026-10-03: a 176-thread WCP point was the stand-in for a missing 180-thread WCP cell
                    ##      (the WCP grid stopped at 16 chains); it is dropped wherever the same configuration
                    ##      also has a 180-thread point (the data are unchanged):
                    ## ---- 2026-10-04: the 176-thread WCP points are measured allocations in their own right
                    ##      (16 x 11, 8 x 22 and 4 x 44 at N >= 10,000), i.e. the fastest configuration at
                    ##      N_threads = 176, so they are plotted again (the drop below is kept, commented out):
                    # configuration_keys <-  paste(df_dev$N_num, df_dev$Algorithm_label)
                    # has_180 <-  configuration_keys[df_dev$n_threads == 180]
                    # drop_176 <-  df_dev$n_threads == 176 & configuration_keys %in% has_180
                    # df_dev <-  df_dev[!drop_176, , drop = FALSE]
                    ##
                    ## ---- 2026-10-06: N_threads = 176 vs. 180 is not a real difference, so the
                    ##      176-thread points are no longer plotted (dropped from the plotted data
                    ##      only; the values and markers CSVs are unchanged). Hence the "176/180"
                    ##      tick below is not triggered and no open points are drawn:
                    df_dev <-  df_dev[df_dev$n_threads != 176, , drop = FALSE]
                    ##
                    ## Show useful budget labels and endpoints; retain all measured points.
                    tick_candidates <-  if (dev == "HPC") c(1, 2, 4, 8, 16, 32, 64, 96, 128, 180) else c(1, 2, 4, 8, 16)
                    ##
                    x_breaks <-  sort(unique(c(range(df_dev$n_threads),
                                               intersect(tick_candidates, df_dev$n_threads))))
                    x_labels <-  as.character(x_breaks)
                    ##
                    ## ---- the 176-thread (WCP) and 180-thread points are too close to label separately on a log scale: one tick,
                    ##      "176/180", midway between them:
                    if (all(c(176, 180) %in% df_dev$n_threads)) {
                        keep_break <-  !(x_breaks %in% c(176, 180))
                        x_breaks   <-  c(x_breaks[keep_break], 178)
                        x_labels   <-  c(x_labels[keep_break], "176/180")
                    }
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
                    ## ---- Optional light markers (show_markers = TRUE): a dotted vertical line at the thread count after which SMT is in
                    ##      use, and after which the L3 cache per active thread falls below the L3 cache per core (4 MB HPC, 2 MB laptop;
                    ##      equation eq:paper1_auto_n_chunks). Changes at the same thread count share one line and label:
                    ##
                    marker_layers <-  NULL
                    ##
                    if (show_markers) {

                        thread_markers <-  fn_paper1_thread_markers(device = dev, n_threads_max = max(df_dev$n_threads))
                        ##
                        if (nrow(thread_markers)) {

                            thread_markers$vjust <-  1.15 + 2.7 * (seq_len(nrow(thread_markers)) - 1)
                            ##
                            marker_layers <-  list( ggplot2::geom_vline( data = thread_markers,
                                                                         mapping = ggplot2::aes(xintercept = after_n_threads),
                                                                         colour = "grey55",
                                                                         linetype = "dotted",
                                                                         linewidth = 1.1),
                                                    ggplot2::geom_label( data = thread_markers,
                                                                         mapping = ggplot2::aes(x = after_n_threads, y = Inf, label = label, vjust = vjust),
                                                                         inherit.aes = FALSE,
                                                                         hjust = 1.04,
                                                                         colour = "grey30",
                                                                         size = 5.5,
                                                                         lineheight = 0.95,
                                                                         label.size = 0,
                                                                         fill = "white",
                                                                         label.padding = ggplot2::unit(0.1, "lines")))
                            ##
                            thread_markers_all <-  rbind(thread_markers_all, thread_markers[, c("device", "after_n_threads", "label")])

                        }

                    }
                    ##
                    ## ---- 2026-10-04: at N_threads = 176, only the WCP allocations with
                    ##      N_threads/chain = 11, 22 or 44 were run (16 x 11, 8 x 22 and 4 x 44,
                    ##      which cannot use 180 threads), so each 176-thread point is the fastest
                    ##      of these few allocations only. Joined to the 90 x 2 point at 180 threads,
                    ##      it drew a spurious dip at the shared "176/180" tick; the 176-thread points
                    ##      are therefore drawn as open points, not joined to the lines
                    ##      (no point is dropped, and the plotted values CSV is unchanged):
                    is_176_point <-  df_dev$n_threads == 176
                    df_line      <-  df_dev[!is_176_point, , drop = FALSE]
                    df_176       <-  df_dev[is_176_point, , drop = FALSE]
                    ##
                    open_176_layer <-  NULL
                    if (nrow(df_176)) {

                        open_176_layer <-  ggplot2::geom_point( data        = df_176,
                                                                size        = 5,
                                                                shape       = 21,
                                                                fill        = "white",
                                                                stroke      = 2,
                                                                show.legend = FALSE)

                    }
                    ##
                    p <-  ggplot( df_dev,
                                  aes( x = n_threads,
                                       y = y_val,
                                       colour = Algorithm_label,
                                       group  = Algorithm_label)) +
                        marker_layers +
                        # geom_point(size = 5) +
                        # geom_line(linewidth = 2) +
                        ggplot2::geom_point(data = df_line, size = 5) +
                        ggplot2::geom_line(data = df_line, linewidth = 2) +
                        open_176_layer +
                        reference_layer +
                        theme_bw(base_size = 28) +
                        theme( legend.position = ifelse(dev == "Laptop", "bottom", "none"),
                               legend.text = element_text(size = 20),
                               axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
                        colour_scale +
                        guides(colour = guide_legend(title = NULL, ncol = 2)) +
                        ylab(y_lab) +
                        ## xlab(expression(log[2](N[threads]~total))) +
                        ## scale_x_continuous(breaks = x_breaks, trans = "log2") +
                        xlab(expression(N[threads]~"(log"[2]~"scale)")) +
                        scale_x_continuous(breaks = x_breaks, labels = x_labels, trans = "log2") +
                        theme(legend.text = element_text(size = 20, family = "mono")) +
                        facet_wrap(~ N_label, scales = "free") +
                        ## ggtitle(ifelse(dev == "HPC", "Local HPC", "Laptop"))
                        ggtitle(ifelse(dev == "HPC", "local-HPC", "Laptop"))
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
                if (!is.null(markers_file) && !is.null(thread_markers_all)) {

                    thread_markers_all$label <-  gsub("\n", " ", thread_markers_all$label, fixed = TRUE)
                    ##
                    utils::write.csv(x = thread_markers_all, file = markers_file, row.names = FALSE)

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
                    ## the markers figure has its own file name, so the figure without markers is kept:
                    if (show_markers) scaling_figure_filename <-  sub("\\.png$", "_markers.png", scaling_figure_filename)
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






















