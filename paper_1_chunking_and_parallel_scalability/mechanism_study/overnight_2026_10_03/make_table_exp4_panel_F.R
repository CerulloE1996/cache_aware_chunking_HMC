##
## ==============================================================================================================
## make_table_exp4_panel_F.R
##
## LaTeX rows of panel F of the E4 profiling table (table:exp4_profiling): the long runs of 3 Oct 2026 on the
## local-HPC (N = 50,000, 180 chains with one thread each, ~600 s each from an idle CPU; design (v)).
## Reads the summaries written by analyse_temperature_long.R (temperature, power, DRAM traffic) and
## analyse_overnight_counters.R (top-down and floating-point counters), and prints the rows to paste into Main.tex.
##
{
      overnight_dir <-  path.expand(paste0("~/Documents/Work/PhD_work/Alg_paper_analysis/",
                                           "paper_1_chunking_and_parallel_scalability/mechanism_study/",
                                           "overnight_2026_10_03"))
      temperature_summary <-  utils::read.csv(file = file.path(overnight_dir, "temperature_long",
                                                               "temperature_long_summary.csv"),
                                              stringsAsFactors = FALSE)
      topdown_summary <-  utils::read.csv(file = file.path(overnight_dir, "summary_topdown.csv"), stringsAsFactors = FALSE)
      fp_summary      <-  utils::read.csv(file = file.path(overnight_dir, "summary_fp.csv"),      stringsAsFactors = FALSE)
      ##
      row_names <-  c(Mplus_standard           = "Mplus standard (BITERATIONS)",
                      MD_BayesMVP_chunks1      = "\\texttt{MD\\_BayesMVP}, $N_{\\text{chunks}} = 1$",
                      MD_BayesMVP_chunking_500 = "\\texttt{MD\\_BayesMVP\\_chunking}, $N_{\\text{chunks}} = 500$",
                      AD_Stan_chunks1          = "\\texttt{AD\\_Stan}, $N_{\\text{chunks}} = 1$",
                      AD_Stan_tape_chunked_250 = "\\texttt{AD\\_Stan\\_tape\\_chunked}, $N_{\\text{chunks}} = 250$")
}
##
fn_value <-  function( table,
                       workload,
                       column,
                       digits
) {

        value <-  table[table$workload == workload, column]
        if (length(value) == 0 || is.na(value[1])) return("-")
        formatC(value[1], format = "f", digits = digits)

}
##
{
      latex_rows <-  vapply(X = names(row_names), FUN.VALUE = character(1), FUN = function(workload) {
            paste0(row_names[[workload]], " & 180 & ",
                   fn_value(temperature_summary, workload, "plateau_Tctl",          1), " & ",
                   fn_value(temperature_summary, workload, "sampling_power_watts",  0), " & ",
                   fn_value(topdown_summary,     workload, "memory_bound_percent",  0), " & ",
                   fn_value(temperature_summary, workload, "DRAM_GB_per_second",    0), " & ",
                   fn_value(fp_summary,          workload, "fp_ops_per_core_cycle", 2), " \\\\")
      })
      idle_row <-  paste0("Idle & - & ", formatC(temperature_summary$idle_Tctl[1], format = "f", digits = 1),
                          " & ", formatC(temperature_summary$idle_power_watts[1], format = "f", digits = 0),
                          " & - & - & - \\\\")
      message(paste0("\033[36m", "Panel F rows (paste into Main.tex):", "\033[0m"))
      cat(c(idle_row, latex_rows), sep = "\n")
      writeLines(text = c(idle_row, latex_rows), con = file.path(overnight_dir, "table_exp4_panel_F_rows.tex"))
      message(paste0("\033[32m", "Saved table_exp4_panel_F_rows.tex", "\033[0m"))
}























