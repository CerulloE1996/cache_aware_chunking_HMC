


##
## ---- Paper 1, E1 Parts II and IV: tables of the selected N_chunks of the chunking + WCP configurations
##      (MD_BayesMVP_WCP_chunking and AD_Stan_WCP_chunking), for each device, allocation (N_chains x N_threads/chain) and N,
##      each with the throughput relative to WCP-only (MD_BayesMVP_WCP / AD_Stan_WCP, i.e. N_chunks = N_threads/chain) at the
##      same allocation. Bold marks the fastest chunking + WCP allocation at each N (the configuration used in E2), and a
##      dagger marks the fastest WCP-only allocation.
##
##      Input:  the chunk-search CSVs of the extended chunk grid (one per device, algorithm and N).
##      Output: one .tex table per algorithm (same labels as the tables they replace).
##
## ---- 2026-10-04: the bracket is now the throughput relative to the FASTEST WCP-only allocation on the same
##      device at the same N (one baseline per device and N), not relative to WCP-only at the same allocation
##      (lines 7-8 above): with the same-allocation baseline, a slow WCP-only allocation inflated its bracket, so
##      the largest bracket was not the fastest configuration. The same-allocation ratio is still written to the
##      *_values.csv files (ratio_vs_WCP_only). Cells are "N_chunks (ratio x)", with a section mark on the
##      largest N_chunks tested at that N, as in Main.tex.
##
##
# data_dir   <-  "paper_1_computational_outputs/manuscript_outputs_final_both_devices_extended_chunk_grid/data"
## 2026-10-03: the export with the narrow-WCP local-HPC runs (up to 90 chains x 2 threads):
data_dir   <-  "paper_1_computational_outputs/manuscript_outputs_final_both_devices_narrow_WCP_2026_10_03/data"
# output_dir <-  path.expand("~/Paper1_upload_WCP_only_tables/Files/Supplement/assets/WCP_selected_chunks/tables")
# output_dir <-  path.expand("~/Documents/Work/PhD_work/Alg_papers_LaTeX/paper_1_v42_2026_10_03/tables_v45")
## 2026-10-04: the re-based tables go to their own folder (the tables_v45 outputs above are kept unchanged):
output_dir <-  path.expand(paste0("~/Documents/Work/PhD_work/Alg_papers_LaTeX/paper_1_v46_2026_10_03/",
                                  "tables_2026_10_04_rebased"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
##
N_values <-  c(500, 2500, 10000, 50000)
##
table_settings <-  list(
    AD_Stan_WCP     = list( label         = "S:table:WCP_selected_chunks_Stan",
                            caption_size  = "\\scriptfootnotesize",
                            column_header = "Stan model ($N$)",
                            file_name     = "table_WCP_selected_chunks_Stan_vs_WCP_only.tex"),
    MD_BayesMVP_WCP = list( label         = "S:table:WCP_selected_chunks",
                            caption_size  = "\\footnotesize",
                            column_header = "BayesMVP ($N$)",
                            file_name     = "table_WCP_selected_chunks_BayesMVP_vs_WCP_only.tex")
)
##
fn_texttt <-  function(model_name) paste0("\\texttt{", gsub("_", "\\\\_", model_name), "}")
##
##
for (algorithm in names(table_settings)) {

        settings <-  table_settings[[algorithm]]
        ##
        cells <-  do.call(rbind, lapply(c("HPC", "Laptop"), function(device) do.call(rbind, lapply(N_values, function(N) {

                chunk_search <-  utils::read.csv( file.path(data_dir, paste0(device, "_paper_1_N_", N, "_algorithm_", algorithm,
                                                                             "_chunk_search.csv")),
                                                  stringsAsFactors = FALSE)
                ##
                allocations <-  unique(chunk_search[, c("n_chains", "threads_per_chain")])
                ##
                per_allocation <-  do.call(rbind, lapply(seq_len(nrow(allocations)), function(i) {

                        rows <-  chunk_search[chunk_search$n_chains == allocations$n_chains[i] &
                                              chunk_search$threads_per_chain == allocations$threads_per_chain[i], , drop = FALSE]
                        ##
                        selected <-  rows[rows$selected_best_chunks, , drop = FALSE]
                        WCP_only <-  rows[rows$num_chunks == rows$threads_per_chain, , drop = FALSE]
                        ##
                        if (nrow(selected) != 1) stop(paste0("Expected one selected N_chunks for ", algorithm, ", ", device, ", N = ", N,
                                                             ", ", allocations$n_chains[i], " x ", allocations$threads_per_chain[i]))
                        if (nrow(WCP_only) != 1) stop(paste0("No WCP-only row (N_chunks = N_threads/chain) for ", algorithm, ", ", device,
                                                             ", N = ", N, ", ", allocations$n_chains[i], " x ", allocations$threads_per_chain[i]))
                        ##
                        data.frame( device                      = device,
                                    N                           = N,
                                    n_chains                    = allocations$n_chains[i],
                                    threads_per_chain           = allocations$threads_per_chain[i],
                                    selected_num_chunks         = selected$num_chunks,
                                    selected_total_iter_per_sec = selected$total_iter_per_sec,
                                    WCP_only_total_iter_per_sec = WCP_only$total_iter_per_sec,
                                    ratio_vs_WCP_only           = selected$chain_rate / WCP_only$chain_rate,
                                    selected_is_WCP_only        = selected$num_chunks == selected$threads_per_chain,
                                    ## 2026-10-04: near-tie flag - the second-best N_chunks of this allocation is within 1%:
                                    near_tie                    = nrow(rows) > 1 &&
                                                                  sort(rows$total_iter_per_sec, decreasing = TRUE)[2] >=
                                                                  0.99 * max(rows$total_iter_per_sec),
                                    stringsAsFactors            = FALSE)

                }))
                ##
                ## the fastest chunking + WCP allocation (used in E2) and the fastest WCP-only allocation, at this device and N:
                per_allocation$fastest_chunking_WCP <-  per_allocation$selected_total_iter_per_sec == max(per_allocation$selected_total_iter_per_sec)
                per_allocation$fastest_WCP_only     <-  per_allocation$WCP_only_total_iter_per_sec == max(per_allocation$WCP_only_total_iter_per_sec)
                ##
                ## ---- 2026-10-04: one baseline per device and N - the fastest WCP-only allocation (the dagger
                ##      row); total iterations/second are comparable across allocations here, because N_iter is
                ##      the same for every allocation of this arm at this N:
                per_allocation$ratio_vs_fastest_WCP_only <-  per_allocation$selected_total_iter_per_sec /
                                                             max(per_allocation$WCP_only_total_iter_per_sec)
                ## the largest N_chunks tested at this device and N (section mark in the table):
                per_allocation$largest_tested_num_chunks <-  per_allocation$selected_num_chunks ==
                                                             max(chunk_search$num_chunks)
                ##
                per_allocation

        }))))
        ##
        ## ---- Cell text: "N_chunks (ratio vs. WCP-only)"; "(1)" when WCP-only itself was selected:
        # cells$cell <-  paste0( cells$selected_num_chunks, " (",
        #                        ifelse( cells$selected_is_WCP_only, "1",
        #                                formatC(cells$ratio_vs_WCP_only, format = "f", digits = 2)), ")")
        ## ---- 2026-10-04: cell text "N_chunks (ratio vs. the fastest WCP-only allocation, with the times
        ##      sign)", with a section mark on the largest N_chunks tested at that N:
        cells$cell <-  paste0( cells$selected_num_chunks,
                               ifelse(cells$largest_tested_num_chunks, "$^{\\S}$", ""),
                               ifelse(cells$near_tie, "$^{*}$", ""),
                               " (", formatC(cells$ratio_vs_fastest_WCP_only, format = "f", digits = 2),
                               "$\\times$)")
        cells$cell <-  ifelse(cells$fastest_chunking_WCP, paste0("\\textbf{", cells$cell, "}"), cells$cell)
        ## cells$cell <-  ifelse(cells$fastest_WCP_only, paste0(cells$cell, "$^{\\dagger}$"), cells$cell)
        ## (the text dagger keeps every table row within the 114-character source width)
        cells$cell <-  ifelse(cells$fastest_WCP_only, paste0(cells$cell, "\\dag"), cells$cell)
        ##
        utils::write.csv( x = cells,
                          file = file.path(output_dir, sub("\\.tex$", "_values.csv", settings$file_name)),
                          row.names = FALSE)
        ##
        ## ---- Rows: device, then N_chains, then N_threads/chain (as in the tables they replace):
        allocation_rows <-  unique(cells[, c("device", "n_chains", "threads_per_chain")])
        allocation_rows <-  allocation_rows[order(match(allocation_rows$device, c("HPC", "Laptop")),
                                                  allocation_rows$n_chains, allocation_rows$threads_per_chain), ]
        ##
        table_rows <-  vapply(seq_len(nrow(allocation_rows)), function(i) {

                row_cells <-  vapply(N_values, function(N) {

                        cell <-  cells$cell[cells$device == allocation_rows$device[i] & cells$N == N &
                                            cells$n_chains == allocation_rows$n_chains[i] &
                                            cells$threads_per_chain == allocation_rows$threads_per_chain[i]]
                        ##
                        if (length(cell) == 0) "-" else cell

                }, character(1))
                ##
                first_of_device <-  i == 1 || allocation_rows$device[i] != allocation_rows$device[i - 1]
                device_text     <-  if (!first_of_device) "" else if (allocation_rows$device[i] == "HPC") "local-HPC" else "Laptop"
                ##
                paste0( if (first_of_device && i > 1) "\\midrule\n" else "",
                        device_text, " & $", allocation_rows$n_chains[i], " \\times ", allocation_rows$threads_per_chain[i], "$ & ",
                        paste(row_cells, collapse = " & "), " \\\\")

        }, character(1))
        ##
        WCP_only_name     <-  fn_texttt(algorithm)
        WCP_chunking_name <-  fn_texttt(paste0(algorithm, "_chunking"))
        ##
        tex <-  c( "%%%%",
                   "%%%%",
                   "\\begin{table}[H]",
                   "%%%%",
                   "\\centering",
                   "\\small",
                   "%%%%",
                   "\\caption{",
                   paste0(settings$caption_size, "{"),
                   paste0("        Selected $N_{\\text{chunks}}$ for ", WCP_chunking_name, ","),
                   "        for each device, allocation ($N_{\\text{chains}} \\times N_{\\text{threads/chain}}$),",
                   "        and $N$ (i.e., the $N_{\\text{chunks}}$, with the lowest mean run time for that allocation;",
                   "        see section \\ref{section:paper1_chunk_wcp_selection_design}).",
                   "        %%",
## (the eight caption lines before 2026-10-04, with the same-allocation bracket:)
# paste0("        The number in brackets is the throughput of this configuration relative to ", WCP_only_name),
# "        (i.e., WCP-only, with $N_{\\text{chunks}} = N_{\\text{threads/chain}}$) at the same allocation;",
# "        (1) means that WCP-only itself was selected.",
# "        %%",
# paste0("        Bold marks the fastest ", WCP_chunking_name, " allocation at each $N$"),
# "        (i.e., the configuration which we used in E2;",
# "        see section \\ref{section:methods:experiment_2_cross_algorithm_software}),",
# paste0("        and \\dag{} marks the fastest ", WCP_only_name, " allocation."),
                   ## ---- 2026-10-04: the re-based bracket (one baseline per device and N), as in Main.tex:
                   "        The number in brackets is the throughput of this configuration relative to",
                   paste0("        the fastest ", WCP_only_name, " allocation on that device at that $N$"),
                   "        (\\dag{}; i.e., WCP-only, with $N_{\\text{chunks}} = N_{\\text{threads/chain}}$),",
                   "        so that, for each device, the largest number in each column",
                   "        is the fastest configuration (bold);",
                   "        an $N_{\\text{chunks}}$ equal to $N_{\\text{threads/chain}}$",
                   "        means that WCP-only itself was selected.",
                   "        %%",
                   paste0("        Bold marks the fastest ", WCP_chunking_name, " allocation at each $N$,"),
                   paste0("        and \\dag{} marks the fastest ", WCP_only_name, " allocation;"),
                   "        in E2 (see section \\ref{section:methods:experiment_2_cross_algorithm_software}),",
                   "        we used the fastest allocation at each $N_{\\text{threads}}$",
                   "        (see section \\ref{section:paper1_chunk_wcp_selection_design}).",
                   "        %%",
                   "        Hyphens mark allocations which were not tested at that $N$.",
                   "        $^{\\S}$the largest $N_{\\text{chunks}}$ tested for that $N$;",
                   "        $^{*}$the second-best $N_{\\text{chunks}}$ for that allocation was within $1\\%$ (i.e., essentially tied).",
                   "}}",
                   "%%%%",
                   paste0("\\label{", settings$label, "}"),
                   "%%%%",
                   "\\PaperOneFitTable{",
                   "\\begin{tabular}{llrrrr}",
                   "%%%%",
                   "\\toprule",
                   "%%%%",
                   paste0(" & & \\multicolumn{4}{c}{", settings$column_header, "} \\\\"),
                   "\\cmidrule(lr){3-6}",
                   "Device & $N_{\\text{chains}} \\times N_{\\text{threads/chain}}$ & $500$ & $2500$ & $10,000$ & $50,000$ \\\\",
                   "%%%%",
                   "\\midrule",
                   "%%%%",
                   table_rows,
                   "%%%%",
                   "\\bottomrule",
                   "%%%%",
                   "\\end{tabular}",
                   "}",
                   "%%%%",
                   "\\end{table}",
                   "%%%%",
                   "%%%%")
        ##
        writeLines(text = tex, con = file.path(output_dir, settings$file_name))
        ##
        message(paste0("\033[36m", "Written: ", file.path(output_dir, settings$file_name), "\033[0m"))

}
























