##
## =====================================================================================================================================
## make_table_exp4_profiling.R
##
## Reproduces Paper 1's table "table:exp4_profiling" (experiment 4, hardware-counter profiling; panels A, B and C) and every
## experiment-4 value quoted in the text of the paper, using only the published CSVs in the mechanism_study folder, and checks each
## reproduced value against the value printed in the paper.
##
## Case labels in the CSVs use the run-time experiment codes of make_mechanism_cases.py. These map onto the designs of the paper
## (section "Plan for experiment 4") as follows:
##
##     E1 = design (i)   : N_chunks x N_chains, one thread per chain, unpinned (local-HPC and laptop; N = 10,000 and 50,000)
##     E2 = design (ii)  : 8 chains packed onto one CCD ("oneCCD", CPUs 0-7) vs spread over eight CCDs ("spread8CCD",
##                         CPUs 0,8,...,56); local-HPC only, N = 50,000
##     E3 = design (iii) : 1 chain x 8 threads, WCP-only (N_chunks = 8) or chunking + WCP; local-HPC "oneCCD" / "spread8CCD",
##                         laptop unpinned
##     E4 = design (iv)  : SMT. Local-HPC: 16 chains on 8 cores of one CCD ("SMT_oneCCD", CPUs 0-7,96-103);
##                         laptop: 8 chains on its 8 physical cores ("physical_cores_only", CPUs 0,2,...,14)
##     E5 = the DRAM-traffic runs of the key design-(i) cases (local-HPC, N = 50,000; cases_HPC_E5.csv), run under both
##          umc_stat (memory-controller counters) and pmc_stat (core counters)
##
## Definitions (as stated in the paper, and as used in analyse_mechanism_study.R). Each case recorded three counter snapshots, after
## the untimed warm-up call, the short sampling call and the long sampling call, so for any counter
##
##     per-iteration value = ((snapshot_3 - snapshot_2) - (snapshot_2 - snapshot_1)) / (n_iter_long - n_iter_short),
##
## and the time per iteration is the paper's two-run estimate, (T_long - T_short) / (n_iter_long - n_iter_short). Hence:
##
##     throughput (chain-iterations per second)          = N_chains / seconds per iteration
##     IPC (instructions per cycle)                      = instructions per iteration / cycles per iteration (all threads)
##     DRAM traffic (GB/s)                               = (read + write bytes per iteration, summed over the 12 memory
##                                                         controllers) / seconds per iteration / 1e9
##     DRAM reads per chain-iteration (GB)               = read bytes per iteration at the memory controllers / N_chains / 1e9
##     data fetched from RAM per chain-iteration (MB)    = 64 x (L1 data-cache fills from DRAM per iteration) / N_chains / 1e6
##
## Outputs (written to the mechanism_study folder):
##
##     table_exp4_profiling_reproduced.csv   one row per table cell: reproduced value (unrounded and as printed), value in the
##                                           paper, match, and the case label it was computed from
##     table_exp4_profiling_reproduced.tex   the table's tabular, in the same layout as the paper, filled with the reproduced values
##
## Matching rules. Table cells, and values printed without "~", match when the reproduced value rounded to the printed number of
## decimals gives the printed text exactly. Values printed with "~" (approximate) match when the reproduced value rounds to the
## printed value at the printed precision, or lies within 5% of it. Qualitative statements are checked as TRUE/FALSE conditions.
##
##
## ---- Folder of the published CSVs ---------------------------------------------------------------------------------------------------
##
{
      ## Default: the analysis folder. When the folder of this script can be found (Rscript, or source()) and holds the CSVs, that
      ## folder is used instead, so the script also runs from a copy of the GitHub repository.
      mechanism_dir <-  path.expand("~/Documents/Work/PhD_work/Alg_paper_analysis/paper_1_chunking_and_parallel_scalability/mechanism_study")
      ##
      script_file <-  sub(pattern = "^--file=", replacement = "", x = grep(pattern = "^--file=", x = commandArgs(trailingOnly = FALSE), value = TRUE))
      if (length(script_file) != 1) script_file <-  tryCatch(expr = sys.frame(1)$ofile, error = function(e) NULL)
      if (length(script_file) == 1 && file.exists(file.path(dirname(normalizePath(script_file)), "mechanism_times_HPC.csv"))) {
            mechanism_dir <-  dirname(normalizePath(script_file))
      }
}

##
## ---- Console colours ---------------------------------------------------------------------------------------------------------------
##
{
      ## NicoStan's colourise() (cyan = information, green = match, red = mismatch) when NicoStan is installed; plain text otherwise.
      use_NicoStan_colours <-  FALSE
      invisible(utils::capture.output(
            use_NicoStan_colours <-  suppressPackageStartupMessages(suppressWarnings(requireNamespace("NicoStan", quietly = TRUE)))
      ))
      ##
      fn_colour <-  function( text,
                              colour
      ) {

              if (use_NicoStan_colours) return(NicoStan::colourise(text = text, fg = colour))
              text

      }
      ##
      message(fn_colour(paste0("Reproducing table:exp4_profiling and the experiment-4 values in the text, from the CSVs in:\n    ", mechanism_dir), "cyan"))
}

##
## ---- Values printed in the paper (Main.tex): TABLE table:exp4_profiling ------------------------------------------------------------
##
## Typed from the paper's tabular, one line per table row. These are the values under audit; they are not used in any calculation.
##
{
      paper_panel_A <-  utils::read.table( text = "
            implementation                         | chains | Throughput | IPC  | DRAM traffic (GB/s) | DRAM reads per chain-iteration (GB)
            BayesMVP, 1 chunk                      | 8      | 31.0       | 0.86 | 113                 | 2.25
            BayesMVP, 1 chunk                      | 96     | 35.3       | 0.09 | 346                 | 7.07
            BayesMVP, 1 chunk                      | 180    | 26.9       | 0.03 | 357                 | 9.35
            BayesMVP, 100 chunks                   | 8      | 46.1       | 1.63 | 12                  | 0.21
            BayesMVP, 100 chunks                   | 96     | 292.7      | 1.50 | 110                 | 0.25
            BayesMVP, 100 chunks                   | 180    | 355.7      | 0.95 | 159                 | 0.29
            Stan model, 1 chunk                    | 8      | 25.3       | 2.34 | 150                 | 3.44
            Stan model, 1 chunk                    | 96     | 42.8       | 0.37 | 321                 | 4.76
            Stan model, 1 chunk                    | 180    | 39.9       | 0.18 | 325                 | 5.41
            Stan model + tape chunking, 250 chunks | 8      | 22.2       | 3.11 | 17                  | 0.60
            Stan model + tape chunking, 250 chunks | 96     | 128.1      | 2.84 | 105                 | 0.52
            Stan model + tape chunking, 250 chunks | 180    | 99.0       | 0.97 | 260                 | 1.57
            ", sep = "|", header = TRUE, colClasses = "character", strip.white = TRUE, check.names = FALSE)
      ##
      paper_panel_B <-  utils::read.table( text = "
            implementation                         | chains | One CCD | Eight CCDs
            BayesMVP, 1 chunk                      | 8      | 5.7     | 30.5
            BayesMVP, 500 chunks                   | 8      | 44.8    | 44.2
            Stan model, 1 chunk                    | 8      | 8.9     | 27.1
            Stan model + tape chunking, 250 chunks | 8      | 21.7    | 22.7
            ", sep = "|", header = TRUE, colClasses = "character", strip.white = TRUE, check.names = FALSE)
      ##
      paper_panel_C <-  utils::read.table( text = "
            implementation                    | 8 chains x 1 thread (chunked) | 1 chain x 8 threads (chunking + WCP) | 1 chain x 8 threads (WCP-only)
            BayesMVP, local-HPC, one CCD      | 44.8                          | 33.0                                 | 15.9
            BayesMVP, local-HPC, eight CCDs   | 44.2                          | 22.6                                 | 24.0
            Stan model, local-HPC, one CCD    | 21.7                          | 16.5                                 | 8.8
            Stan model, local-HPC, eight CCDs | 22.7                          | 13.3                                 | 13.6
            BayesMVP, laptop                  | 21.9                          | 14.5                                 | 6.1
            Stan model, laptop                | 9.7                           | 6.1                                  | 4.1
            ", sep = "|", header = TRUE, colClasses = "character", strip.white = TRUE, check.names = FALSE)
}

##
## ---- Values printed in the paper (Main.tex): EXPERIMENT-4 VALUES QUOTED IN THE TEXT ------------------------------------------------
##
## One line per quoted value: identifier, location (LaTeX section label), the quoted wording, the value as printed, the matching rule
## ("rounded", "signed" = rounded percentage change with its sign, "approximate", "logical", "text"), the printed number of decimals,
## and the source ("CSV" = reproduced from the published CSVs; "derived" = arithmetic on constants stated in the paper, e.g. the
## 1,608 bytes per individual, the cache sizes and the memory speed, with no CSV input except the design's chain counts).
##
{
      fn_paper <-  function( id,
                             location,
                             quoted,
                             value_in_paper,
                             rule,
                             digits = 0,
                             source = "CSV"
      ) {

              data.frame( id = id, location = location, quoted = quoted, value_in_paper = value_in_paper, rule = rule, digits = digits,
                          source = source, stringsAsFactors = FALSE)

      }
      ##
      plan_section <-  "section:ps1_ps2_plan_for_experiment_4"
      supp_section <-  "S:section:supp_profiling_commands_running"
      results_section <-  "section:ps1_ps2_results_experiment_4"
      auto_section <-  "section:results:experiment_1:auto_N_chunks_BayesMVP"
      discussion_section <-  "section:ps1_ps2_discussion"
      conclusion_section <-  "section:ps1_ps2_discussion_conclusion"
      ##
      paper_text <-  do.call(what = rbind, args = list(
            ##
            ## ---- Plan for experiment 4 (design and measurement)
            fn_paper("plan_idle_background_GB_per_s",       plan_section, "every value includes ~4 GB/s of idle background traffic",            "4",                 "approximate", 0),
            fn_paper("plan_n_iter_short",                   plan_section, "short runs of 2 iterations",                                          "2",                 "text"),
            fn_paper("plan_n_iter_long_BayesMVP_N50000",    plan_section, "long runs of 20 [...] for NicoStan+BayesMVP [...] at N = 50,000",     "20",                "text"),
            fn_paper("plan_n_iter_long_BayesMVP_N10000",    plan_section, "long runs of [...] 50 for NicoStan+BayesMVP [...] at N = 10,000",     "50",                "text"),
            fn_paper("plan_n_iter_long_Stan_N50000",        plan_section, "10 [...] for the Stan model, at N = 50,000",                          "10",                "text"),
            fn_paper("plan_n_iter_long_Stan_N10000",        plan_section, "[...] 30 for the Stan model, at N = 10,000",                          "30",                "text"),
            fn_paper("plan_design_i_N_values",              plan_section, "design (i): at N = 10,000 and N = 50,000",                            "10000, 50000",      "text"),
            fn_paper("plan_design_i_chains_HPC",            plan_section, "design (i): N_chains in {1, 8, 48, 96, 180} on the local-HPC",        "1, 8, 48, 96, 180", "text"),
            fn_paper("plan_design_i_chains_Laptop",         plan_section, "design (i): N_chains in {1, 4, 8, 16} on the laptop",                 "1, 4, 8, 16",       "text"),
            ##
            ## ---- Supplement: running the profiling cases
            fn_paper("supp_cases_HPC",                      supp_section, "128 cases on the local-HPC",                                          "128",               "text"),
            fn_paper("supp_cases_Laptop",                   supp_section, "85 on the laptop",                                                    "85",                "text"),
            fn_paper("supp_cases_HPC_DRAM",                 supp_section, "the 28 DRAM-traffic cases of design (i)",                             "28",                "text"),
            ##
            ## ---- E4 results: memory bandwidth without chunking
            fn_paper("res_plateau_lower_GB_per_s",          results_section, "DRAM traffic plateaus at ~305-360 GB/s from 48 chains onwards (lower)", "305",       "approximate", 0),
            fn_paper("res_plateau_upper_GB_per_s",          results_section, "DRAM traffic plateaus at ~305-360 GB/s from 48 chains onwards (upper)", "360",       "approximate", 0),
            fn_paper("res_peak_GB_per_s",                   results_section, "theoretical peak of 460.8 GB/s",                                   "460.8",             "rounded",     1, "derived"),
            fn_paper("res_IPC_BayesMVP_1chunk_8",           results_section, "IPC collapses [...] from 0.86 [...] for BayesMVP (8 chains)",      "0.86",              "rounded",     2),
            fn_paper("res_IPC_BayesMVP_1chunk_180",         results_section, "IPC collapses [...] to 0.03 for BayesMVP (180 chains)",            "0.03",              "rounded",     2),
            fn_paper("res_IPC_Stan_1chunk_8",               results_section, "[IPC] 2.34 [...] for the Stan model (8 chains)",                   "2.34",              "rounded",     2),
            fn_paper("res_IPC_Stan_1chunk_180",             results_section, "[IPC] to 0.18 for the Stan model (180 chains)",                    "0.18",              "rounded",     2),
            fn_paper("res_throughput_BayesMVP_1chunk_48",   results_section, "47.1 vs. 26.9 chain-iterations per second for BayesMVP at 48 [...]", "47.1",           "rounded",     1),
            fn_paper("res_throughput_BayesMVP_1chunk_180",  results_section, "47.1 vs. 26.9 [...] at [...] 180 chains",                          "26.9",              "rounded",     1),
            ##
            ## ---- E4 results: DRAM reads and IPC with chunking
            fn_paper("res_reads_BayesMVP_1chunk_96",        results_section, "BayesMVP with 100 chunks, from 7.07 [...] GB at 96 chains",        "7.07",              "rounded",     2),
            fn_paper("res_reads_BayesMVP_100chunks_96",     results_section, "BayesMVP with 100 chunks, [...] to 0.25 GB at 96 chains",          "0.25",              "rounded",     2),
            fn_paper("res_reads_BayesMVP_1chunk_180",       results_section, "from 9.35 [...] GB at 180 chains",                                 "9.35",              "rounded",     2),
            fn_paper("res_reads_BayesMVP_100chunks_180",    results_section, "[...] to 0.29 GB at 180 chains",                                   "0.29",              "rounded",     2),
            fn_paper("res_reads_Stan_1chunk_96",            results_section, "Stan model with tape chunking using 250 chunks, from 4.76 [...] GB", "4.76",            "rounded",     2),
            fn_paper("res_reads_Stan_tape250_96",           results_section, "Stan model with tape chunking using 250 chunks, [...] to 0.52 GB", "0.52",              "rounded",     2),
            fn_paper("res_IPC_BayesMVP_100chunks_96",       results_section, "the IPC remains high (1.50 for BayesMVP [...] at 96 chains)",      "1.50",              "rounded",     2),
            fn_paper("res_IPC_Stan_tape250_96",             results_section, "the IPC remains high ([...] 2.84 for the Stan model, at 96 chains)", "2.84",            "rounded",     2),
            fn_paper("res_throughput_100chunks_increases",  results_section, "BayesMVP's throughput keeps increasing up to 180 threads",         "TRUE",              "logical"),
            fn_paper("res_one_chain_speedup_BayesMVP",      results_section, "with one chain, chunking gave only a ~1.5x speed-up for BayesMVP",  "1.5",              "approximate", 1),
            fn_paper("res_one_chain_Stan_tape250",          results_section, "tape chunking made the Stan model slightly slower (3.0 vs. 3.5 [...])", "3.0",          "rounded",     1),
            fn_paper("res_one_chain_Stan_1chunk",           results_section, "tape chunking made the Stan model slightly slower (3.0 vs. 3.5 [...])", "3.5",          "rounded",     1),
            ##
            ## ---- E4 results: packing eight chains onto one CCD (design (ii); table panel B)
            fn_paper("res_packing_slowdown_BayesMVP",       results_section, "packing made the chains ~5.3x (BayesMVP) [...] slower",            "5.3",               "approximate", 1),
            fn_paper("res_packing_slowdown_Stan",           results_section, "and ~3.0x (Stan model) slower",                                    "3.0",               "approximate", 1),
            fn_paper("res_packed_BayesMVP_500",             results_section, "44.8 vs. 44.2 [...] for BayesMVP with 500 chunks (one CCD)",       "44.8",              "rounded",     1),
            fn_paper("res_spread_BayesMVP_500",             results_section, "44.8 vs. 44.2 [...] for BayesMVP with 500 chunks (eight CCDs)",    "44.2",              "rounded",     1),
            fn_paper("res_packed_Stan_tape250",             results_section, "21.7 vs. 22.7 for the Stan model with 250 chunks (one CCD)",       "21.7",              "rounded",     1),
            fn_paper("res_spread_Stan_tape250",             results_section, "21.7 vs. 22.7 for the Stan model with 250 chunks (eight CCDs)",    "22.7",              "rounded",     1),
            ##
            ## ---- E4 results: working set per chunk (arithmetic on the stated 1,608 bytes per individual)
            fn_paper("res_working_set_KB_per_individual",   results_section, "gradient working set is ~1.6 KB per individual",                   "1.6",               "approximate", 1, "derived"),
            fn_paper("res_working_set_MB_1chunk",           results_section, "~80 MB with 1 chunk",                                              "80",                "approximate", 0, "derived"),
            fn_paper("res_working_set_MB_25chunks",         results_section, "~3.2 MB with 25 chunks",                                           "3.2",               "approximate", 1, "derived"),
            fn_paper("res_working_set_MB_100chunks",        results_section, "~0.80 [MB] with 100 [chunks]",                                     "0.80",              "approximate", 2, "derived"),
            fn_paper("res_working_set_MB_500chunks",        results_section, "~0.16 MB with 500 chunks",                                         "0.16",              "approximate", 2, "derived"),
            fn_paper("res_25chunks_fit_L3_not_L2",          results_section, "25 chunks: fits in the L3, but not the L2, cache",                 "TRUE",              "logical",     0, "derived"),
            fn_paper("res_100_500chunks_fit_L2",            results_section, "100 and 500 chunks: fit in each core's 1 MB L2 cache",             "TRUE",              "logical",     0, "derived"),
            ##
            ## ---- E4 results: 25 chunks under SMT
            fn_paper("res_25chunks_most_traffic_removed",   results_section, "With 25 chunks, most of the DRAM traffic was already removed at up to 96 chains", "TRUE", "logical"),
            fn_paper("res_DRAM_traffic_25chunks_96",        results_section, "the DRAM traffic rose from 158 [...] GB/s",                        "158",               "rounded",     0),
            fn_paper("res_DRAM_traffic_25chunks_180",       results_section, "the DRAM traffic rose [...] to 306 GB/s",                          "306",               "rounded",     0),
            fn_paper("res_IPC_25chunks_96",                 results_section, "the IPC dropped from 1.01 [...]",                                  "1.01",              "rounded",     2),
            fn_paper("res_IPC_25chunks_180",                results_section, "the IPC dropped [...] to 0.37",                                    "0.37",              "rounded",     2),
            fn_paper("res_SMT_gain_100chunks",              results_section, "from 96 to 180 chains, chunked BayesMVP gained 22% - 24% (100 chunks)", "+22",          "signed",      0),
            fn_paper("res_SMT_gain_500chunks",              results_section, "from 96 to 180 chains, chunked BayesMVP gained 22% - 24% (500 chunks)", "+24",          "signed",      0),
            fn_paper("res_SMT_gain_1chunk",                 results_section, "whereas unchunked BayesMVP lost 24%",                              "-24",               "signed",      0),
            ##
            ## ---- E4 results: two chains per physical core (design (iv)), local-HPC
            fn_paper("res_SMT_HPC_change_500",              results_section, "increased BayesMVP's throughput by 15% with 500 chunks",           "+15",               "signed",      0),
            fn_paper("res_SMT_HPC_16chains_500",            results_section, "(51.6 vs. 44.8 chain-iterations per second)",                      "51.6",              "rounded",     1),
            fn_paper("res_SMT_HPC_8chains_500",             results_section, "(51.6 vs. 44.8 chain-iterations per second)",                      "44.8",              "rounded",     1),
            fn_paper("res_SMT_HPC_change_25",               results_section, "but reduced it by 16% with 25 chunks",                             "-16",               "signed",      0),
            fn_paper("res_SMT_HPC_16chains_25",             results_section, "(32.0 vs. 37.9)",                                                  "32.0",              "rounded",     1),
            fn_paper("res_SMT_HPC_8chains_25",              results_section, "(32.0 vs. 37.9)",                                                  "37.9",              "rounded",     1),
            fn_paper("res_SMT_HPC_change_1",                results_section, "and by 19% with no chunking",                                      "-19",               "signed",      0),
            fn_paper("res_SMT_HPC_16chains_1",              results_section, "(4.6 vs. 5.7)",                                                    "4.6",               "rounded",     1),
            fn_paper("res_SMT_HPC_8chains_1",               results_section, "(4.6 vs. 5.7)",                                                    "5.7",               "rounded",     1),
            fn_paper("res_SMT_HPC_change_Stan_tape250",     results_section, "for the Stan model with 250 tape chunks, it was reduced by 12%",   "-12",               "signed",      0),
            fn_paper("res_SMT_HPC_16chains_Stan_tape250",   results_section, "(19.2 vs. 21.7)",                                                  "19.2",              "rounded",     1),
            fn_paper("res_SMT_HPC_8chains_Stan_tape250",    results_section, "(19.2 vs. 21.7)",                                                  "21.7",              "rounded",     1),
            ##
            ## ---- E4 results: two chains per physical core, laptop (16 chains unpinned, design (i), vs 8 chains pinned, design (iv))
            fn_paper("res_SMT_Laptop_change_500",           results_section, "the corresponding changes were +10% (23.9 vs. 21.7)",              "+10",               "signed",      0),
            fn_paper("res_SMT_Laptop_16chains_500",         results_section, "(23.9 vs. 21.7)",                                                  "23.9",              "rounded",     1),
            fn_paper("res_SMT_Laptop_8chains_500",          results_section, "(23.9 vs. 21.7)",                                                  "21.7",              "rounded",     1),
            fn_paper("res_SMT_Laptop_change_25",            results_section, "-32% (9.7 vs. 14.3)",                                              "-32",               "signed",      0),
            fn_paper("res_SMT_Laptop_16chains_25",          results_section, "(9.7 vs. 14.3)",                                                   "9.7",               "rounded",     1),
            fn_paper("res_SMT_Laptop_8chains_25",           results_section, "(9.7 vs. 14.3)",                                                   "14.3",              "rounded",     1),
            fn_paper("res_SMT_Laptop_change_1",             results_section, "-10% (2.20 vs. 2.45)",                                             "-10",               "signed",      0),
            fn_paper("res_SMT_Laptop_16chains_1",           results_section, "(2.20 vs. 2.45)",                                                  "2.20",              "rounded",     2),
            fn_paper("res_SMT_Laptop_8chains_1",            results_section, "(2.20 vs. 2.45)",                                                  "2.45",              "rounded",     2),
            fn_paper("res_SMT_Laptop_change_Stan_tape250",  results_section, "-28% (5.9 vs. 8.1)",                                               "-28",               "signed",      0),
            fn_paper("res_SMT_Laptop_16chains_Stan_tape250", results_section, "(5.9 vs. 8.1)",                                                   "5.9",               "rounded",     1),
            fn_paper("res_SMT_Laptop_8chains_Stan_tape250", results_section, "(5.9 vs. 8.1)",                                                    "8.1",               "rounded",     1),
            ##
            ## ---- E4 results: eight chunked chains vs one chain with 8 threads (design (iii); table panel C)
            fn_paper("res_WCP_oneCCD_joint_BayesMVP",       results_section, "on one CCD [...] 1.36x (BayesMVP) [...] faster than chunking + WCP", "1.36",            "rounded",     2),
            fn_paper("res_WCP_oneCCD_joint_Stan",           results_section, "on one CCD [...] 1.32x (Stan model) faster than chunking + WCP",   "1.32",              "rounded",     2),
            fn_paper("res_WCP_oneCCD_WCPonly_BayesMVP",     results_section, "and 2.8x [...] faster than WCP-only (BayesMVP, one CCD)",          "2.8",               "rounded",     1),
            fn_paper("res_WCP_oneCCD_WCPonly_Stan",         results_section, "and [...] 2.5x faster than WCP-only (Stan model, one CCD)",        "2.5",               "rounded",     1),
            fn_paper("res_WCP_eightCCD_joint_BayesMVP",     results_section, "over eight CCDs, 2.0x [...] faster than chunking + WCP (BayesMVP)", "2.0",              "rounded",     1),
            fn_paper("res_WCP_eightCCD_joint_Stan",         results_section, "over eight CCDs, [...] 1.7x faster than chunking + WCP (Stan model)", "1.7",            "rounded",     1),
            fn_paper("res_WCP_eightCCD_WCPonly_BayesMVP",   results_section, "and 1.8x [...] faster than WCP-only (BayesMVP, eight CCDs)",       "1.8",               "rounded",     1),
            fn_paper("res_WCP_eightCCD_WCPonly_Stan",       results_section, "and [...] 1.7x faster than WCP-only (Stan model, eight CCDs)",     "1.7",               "rounded",     1),
            fn_paper("res_WCP_Laptop_lower",                results_section, "and on the laptop (a single die), 1.5-3.6x faster (lower)",        "1.5",               "rounded",     1),
            fn_paper("res_WCP_Laptop_upper",                results_section, "and on the laptop (a single die), 1.5-3.6x faster (upper)",        "3.6",               "rounded",     1),
            fn_paper("res_WCP_only_faster_when_spread",     results_section, "WCP-only [...] was faster when spread over eight CCDs",            "TRUE",              "logical"),
            fn_paper("res_joint_slower_spread_lower",       results_section, "whereas chunking + WCP was 19% - 31% slower (lower)",              "19",                "rounded",     0),
            fn_paper("res_joint_slower_spread_upper",       results_section, "whereas chunking + WCP was 19% - 31% slower (upper)",              "31",                "rounded",     0),
            ##
            ## ---- E1 section on the automatic N_chunks rule: experiment-4 values (design (i), core counters)
            fn_paper("auto_working_set_MB_25chunks",        auto_section, "the 25-chunk working set (~3 MB)",                                    "3",                 "approximate", 0, "derived"),
            fn_paper("auto_25chunks_exceeds_L3_cases",      auto_section, "exceeded the L3 cache per active thread (i.e., at 180 chains on the HPC, and at 8 and 16 chains on the laptop)",
                                                                                                                                                 "HPC 180; Laptop 8; Laptop 16", "text", 0, "derived"),
            fn_paper("auto_RAM_ratio_lower",                auto_section, "the data fetched from RAM per chain-iteration [...] was 3-6 times that with 100 chunks (lower)", "3", "rounded", 0),
            fn_paper("auto_RAM_ratio_upper",                auto_section, "the data fetched from RAM per chain-iteration [...] was 3-6 times that with 100 chunks (upper)", "6", "rounded", 0),
            fn_paper("auto_throughput_ratio_lower",         auto_section, "and the throughput was only 0.40-0.66 of it (lower)",                 "0.40",              "rounded",     2),
            fn_paper("auto_throughput_ratio_upper",         auto_section, "and the throughput was only 0.40-0.66 of it (upper)",                 "0.66",              "rounded",     2),
            fn_paper("auto_RAM_MB_HPC180_25chunks",         auto_section, "e.g., 585 vs. 97 MB per chain-iteration [...] at 180 chains on the HPC", "585",            "rounded",     0),
            fn_paper("auto_RAM_MB_HPC180_100chunks",        auto_section, "e.g., 585 vs. 97 MB per chain-iteration [...] at 180 chains on the HPC", "97",             "rounded",     0),
            fn_paper("auto_throughput_ratio_HPC180",        auto_section, "and 0.56 of the throughput, at 180 chains on the HPC",                "0.56",              "rounded",     2),
            fn_paper("auto_L3_fraction_100chunks_96",       auto_section, "100 chunks (0.19 and 0.36 of the L3 cache per active thread at 96 [...] chains)", "0.19",  "rounded",     2, "derived"),
            fn_paper("auto_L3_fraction_100chunks_180",      auto_section, "100 chunks (0.19 and 0.36 of the L3 cache per active thread at [...] 180 chains)", "0.36", "rounded",     2, "derived"),
            fn_paper("auto_RAM_100_close_to_500",           auto_section, "the data fetched from RAM per chain-iteration was close to (or below) that with 500 chunks", "TRUE", "logical"),
            fn_paper("auto_one_chain_25_vs_100_percent",    auto_section, "with one chain on the HPC, 25 chunks gave 12% lower throughput than 100 chunks", "12",     "rounded",     0),
            fn_paper("auto_25chunks_L3_fraction",           auto_section, "despite using only a tenth of the CCD's L3 cache",                    "0.1",               "approximate", 1, "derived"),
            fn_paper("auto_f_range_lower",                  auto_section, "f = 1/4, which lies within the range supported by experiment 4 (~0.2-0.36) (lower)", "0.2", "approximate", 1, "derived"),
            fn_paper("auto_f_range_upper",                  auto_section, "f = 1/4, which lies within the range supported by experiment 4 (~0.2-0.36) (upper)", "0.36", "approximate", 2, "derived"),
            ##
            ## ---- Discussion and conclusion
            fn_paper("disc_plateau_lower_GB_per_s",         discussion_section, "the measured DRAM traffic plateaus at ~305-360 GB/s (lower)",   "305",               "approximate", 0),
            fn_paper("disc_plateau_upper_GB_per_s",         discussion_section, "the measured DRAM traffic plateaus at ~305-360 GB/s (upper)",   "360",               "approximate", 0),
            fn_paper("disc_peak_GB_per_s",                  discussion_section, "below but reasonably close to its theoretical peak of 460.8 GB/s", "460.8",          "rounded",     1, "derived"),
            fn_paper("disc_BayesMVP_reads_reduction",       discussion_section, "reducing the DRAM reads per chain-iteration by up to ~32x for BayesMVP", "32",       "approximate", 0),
            fn_paper("disc_Stan_plateau_lower_GB_per_s",    discussion_section, "the Stan model also saturates the memory bandwidth (~305-325 GB/s at 48 to 180 chains) (lower)", "305", "approximate", 0),
            fn_paper("disc_Stan_plateau_upper_GB_per_s",    discussion_section, "the Stan model also saturates the memory bandwidth (~305-325 GB/s at 48 to 180 chains) (upper)", "325", "approximate", 0),
            fn_paper("disc_Stan_tape_traffic_reduction_96", discussion_section, "tape chunking reduced its DRAM traffic per chain-iteration by ~9x at 96 chains", "9", "approximate", 0),
            fn_paper("concl_BayesMVP_reads_reduction",      conclusion_section, "reducing DRAM reads per chain-iteration by up to ~32x",         "32",                "approximate", 0)
      ))
}

##
## ---- Constants stated in the paper (used only for the "derived" values) ---------------------------------------------------------------
##
{
      bytes_per_individual <-  1608                              ## B_row: working set per individual (binary LC-MVP, T = 6)
      N_profiling <-  50000                                      ## N used for the table and the quoted working-set sizes
      L3_bytes <-  c(HPC = 2^25, Laptop = 2^24)                  ## S_L3: one L3 cache (32 MB per CCD on the local-HPC; 16 MB on the laptop)
      number_of_L3 <-  c(HPC = 12, Laptop = 1)                   ## N_L3: number of L3 caches
      L2_bytes_HPC <-  2^20                                      ## 1 MB L2 cache per core on the local-HPC
      HPC_memory_channels <-  12                                 ## 12 x DDR5-4800 channels, 8 bytes per transfer
      HPC_transfers_per_second <-  4800e6
      bytes_per_transfer <-  8
}

##
## ---- Read the published CSVs -------------------------------------------------------------------------------------------------------
##
{
      fn_read <-  function(file_name, ...) utils::read.csv(file = file.path(mechanism_dir, file_name), stringsAsFactors = FALSE, ...)
      ##
      times <-  rbind(fn_read("mechanism_times_HPC.csv"), fn_read("mechanism_times_Laptop.csv"))
      ##
      counts_HPC <-  fn_read("mechanism_counts_HPC.csv")
      counts_Laptop <-  fn_read("mechanism_counts_Laptop.csv")
      counts <-  rbind(data.frame(device = "HPC", counts_HPC), data.frame(device = "Laptop", counts_Laptop))
      ##
      umc <-  fn_read("mechanism_umc_HPC.csv", header = FALSE, col.names = c("label", "seconds", "status", "DRAM_read_bytes", "DRAM_write_bytes"))
      cal_umc <-  fn_read("cal_umc.csv", header = FALSE, col.names = c("label", "seconds", "status", "DRAM_read_bytes", "DRAM_write_bytes"))
      ##
      cases_HPC <-  fn_read("cases_HPC.csv")
      cases_Laptop <-  fn_read("cases_Laptop.csv")
      cases_HPC_DRAM <-  fn_read("cases_HPC_E5.csv")
      ##
      summary_published <-  fn_read("mechanism_summary.csv")
      bandwidth_published <-  fn_read("mechanism_DRAM_bandwidth_HPC.csv")
      ##
      message(fn_colour(paste0("Read ", nrow(times), " timed cases (", sum(times$device == "HPC"), " local-HPC, ", sum(times$device == "Laptop"),
                               " laptop), ", nrow(counts), " core-counter rows and ", nrow(umc), " memory-controller rows."), "cyan"))
}

##
## ---- Per-case measures, computed from the raw counters and timings -------------------------------------------------------------------
##
{
      core_counters <-  c("cycles", "instructions", "fills_local_L2", "fills_local_L3", "fills_other_CCX", "fills_DRAM")
      ##
      ## Snapshot rows of one case, in snapshot order (NULL when the case has no snapshots in that file):
      fn_snapshots <-  function( counter_rows,
                                 case_label
      ) {

              snapshot_rows <-  counter_rows[grepl(pattern = ":snapshot_[0-9]+$", x = counter_rows$label), , drop = FALSE]
              snapshot_rows <-  snapshot_rows[sub(pattern = ":snapshot_[0-9]+$", replacement = "", x = snapshot_rows$label) == case_label, , drop = FALSE]
              if (nrow(snapshot_rows) == 0) return(NULL)
              if (nrow(snapshot_rows) != 3) stop(paste0("Case ", case_label, " has ", nrow(snapshot_rows), " counter snapshots (3 expected)."))
              snapshot_rows[order(as.numeric(sub(pattern = "^.*:snapshot_", replacement = "", x = snapshot_rows$label))), , drop = FALSE]

      }
      ##
      ## Per-iteration value of one counter: (long-run window - short-run window) / (n_iter_long - n_iter_short):
      fn_per_iteration <-  function( snapshots,
                                     counter,
                                     n_iter_difference
      ) {

              ((snapshots[[counter]][3] - snapshots[[counter]][2]) - (snapshots[[counter]][2] - snapshots[[counter]][1])) / n_iter_difference

      }
      ##
      case_measures <-  do.call(what = rbind, args = lapply(X = seq_len(nrow(times)), FUN = function(row_index) {

              case_row <-  times[row_index, , drop = FALSE]
              n_iter_difference <-  case_row$n_iter_long - case_row$n_iter_short
              seconds_per_iteration <-  (case_row$elapsed_seconds_long_run - case_row$elapsed_seconds_short_run) / n_iter_difference
              ##
              core_snapshots <-  fn_snapshots(counter_rows = counts[counts$device == case_row$device, , drop = FALSE], case_label = case_row$label)
              if (is.null(core_snapshots)) stop(paste0("Case ", case_row$label, " (", case_row$device, ") has no core-counter snapshots."))
              core <-  sapply(X = core_counters, FUN = function(counter) fn_per_iteration(core_snapshots, counter, n_iter_difference))
              ##
              umc_snapshots <-  if (case_row$device == "HPC") fn_snapshots(counter_rows = umc, case_label = case_row$label) else NULL
              read_bytes <-  if (is.null(umc_snapshots)) NA else fn_per_iteration(umc_snapshots, "DRAM_read_bytes", n_iter_difference)
              write_bytes <-  if (is.null(umc_snapshots)) NA else fn_per_iteration(umc_snapshots, "DRAM_write_bytes", n_iter_difference)
              ##
              data.frame( device                              = case_row$device,
                          label                               = case_row$label,
                          experiment                          = sub(pattern = "_.*$", replacement = "", x = case_row$label),
                          algorithm                           = case_row$algorithm,
                          N                                   = case_row$N,
                          num_chunks                          = case_row$num_chunks,
                          n_chains                            = case_row$n_chains,
                          threads_per_chain                   = case_row$threads_per_chain,
                          placement                           = sub(pattern = "^.*_tpc[0-9]+_", replacement = "", x = case_row$label),
                          n_iter_short                        = case_row$n_iter_short,
                          n_iter_long                         = case_row$n_iter_long,
                          seconds_per_iteration               = seconds_per_iteration,
                          throughput                          = case_row$n_chains / seconds_per_iteration,
                          IPC                                 = core[["instructions"]] / core[["cycles"]],
                          RAM_fills_MB_per_chain_iteration    = 64 * core[["fills_DRAM"]] / case_row$n_chains / 1e6,
                          DRAM_total_GB_per_second            = (read_bytes + write_bytes) / seconds_per_iteration / 1e9,
                          DRAM_read_GB_per_chain_iteration    = read_bytes / case_row$n_chains / 1e9,
                          DRAM_total_GB_per_chain_iteration   = (read_bytes + write_bytes) / case_row$n_chains / 1e9,
                          stringsAsFactors                    = FALSE)

      }))
      ##
      message(fn_colour(paste0("Computed the measures of ", nrow(case_measures), " cases (", sum(!is.na(case_measures$DRAM_total_GB_per_second)),
                               " with memory-controller counters)."), "cyan"))
}

##
## ---- Consistency with the published summaries and case lists ---------------------------------------------------------------------------
##
## The recomputed measures should equal those in mechanism_summary.csv and mechanism_DRAM_bandwidth_HPC.csv (written by
## analyse_mechanism_study.R from the same raw CSVs), and every case in the case lists should have a timing row with the same settings.
##
{
      fn_max_relative_difference <-  function(x, y) max(abs(x - y) / pmax(abs(y), 1e-12))
      ##
      summary_keys <-  paste(summary_published$device, summary_published$label)
      measure_keys <-  paste(case_measures$device, case_measures$label)
      summary_matched <-  summary_published[match(measure_keys, summary_keys), , drop = FALSE]
      bandwidth_rows <-  case_measures[!is.na(case_measures$DRAM_total_GB_per_second), , drop = FALSE]
      bandwidth_matched <-  bandwidth_published[match(bandwidth_rows$label, bandwidth_published$label), , drop = FALSE]
      ##
      consistency <-  c( "throughput vs mechanism_summary.csv"                          = fn_max_relative_difference(case_measures$throughput, summary_matched$chain_iterations_per_second),
                         "IPC vs mechanism_summary.csv"                                 = fn_max_relative_difference(case_measures$IPC, summary_matched$instructions_per_cycle),
                         "RAM fills per chain-iteration vs mechanism_summary.csv"       = fn_max_relative_difference(case_measures$RAM_fills_MB_per_chain_iteration, summary_matched$DRAM_MB_per_chain_iteration),
                         "DRAM traffic vs mechanism_DRAM_bandwidth_HPC.csv"             = fn_max_relative_difference(bandwidth_rows$DRAM_total_GB_per_second, bandwidth_matched$DRAM_total_GB_per_second),
                         "DRAM reads per chain-iteration vs mechanism_DRAM_bandwidth_HPC.csv" = fn_max_relative_difference(bandwidth_rows$DRAM_read_GB_per_chain_iteration, bandwidth_matched$DRAM_read_GB_per_chain_iter))
      ##
      message(fn_colour("\n---- Consistency of the recomputed measures with the published summary CSVs (maximum relative difference) ----", "cyan"))
      for (check_name in names(consistency)) {
            consistent <-  is.finite(consistency[[check_name]]) && consistency[[check_name]] < 1e-9
            message(fn_colour(paste0("    ", check_name, ": ", formatC(consistency[[check_name]], format = "e", digits = 1),
                                     if (consistent) "  (consistent)" else "  (DIFFERENT)"), if (consistent) "green" else "red"))
      }
      ##
      ## Case lists vs timing rows (label and settings):
      fn_cases_complete <-  function( case_list,
                                      device
      ) {

              timed <-  times[times$device == device, , drop = FALSE]
              timed <-  timed[match(case_list$label, timed$label), , drop = FALSE]
              all(!is.na(timed$label)) &&
              all(timed$algorithm == case_list$algorithm) && all(timed$N == case_list$N) && all(timed$num_chunks == case_list$chunks) &&
              all(timed$n_chains == case_list$chains) && all(timed$threads_per_chain == case_list$threads_per_chain) &&
              all(timed$n_iter_long == case_list$n_iter)

      }
      cases_complete <-  c( "cases_HPC.csv"    = fn_cases_complete(cases_HPC, "HPC"),
                            "cases_HPC_E5.csv" = fn_cases_complete(cases_HPC_DRAM, "HPC"),
                            "cases_Laptop.csv" = fn_cases_complete(cases_Laptop, "Laptop"))
      for (case_file in names(cases_complete)) {
            message(fn_colour(paste0("    every case in ", case_file, " has a timing row with the same settings: ", cases_complete[[case_file]]),
                              if (cases_complete[[case_file]]) "green" else "red"))
      }
}

##
## ---- Case selection -------------------------------------------------------------------------------------------------------------------
##
{
      ## Returns the single case with these settings; stops when there is not exactly one.
      fn_case <-  function( device,
                            experiment,
                            algorithm,
                            chunks,
                            chains,
                            placement         = "unpinned",
                            threads_per_chain = 1,
                            N                 = N_profiling
      ) {

              selected <-  case_measures[ case_measures$device == device & case_measures$experiment == experiment & case_measures$algorithm == algorithm &
                                          case_measures$N == N & case_measures$num_chunks == chunks & case_measures$n_chains == chains &
                                          case_measures$threads_per_chain == threads_per_chain & case_measures$placement == placement, , drop = FALSE]
              if (nrow(selected) != 1) {
                    stop(paste0("Expected one case, found ", nrow(selected), ": ", device, ", ", experiment, ", ", algorithm, ", N = ", N, ", ",
                                chunks, " chunks, ", chains, " chains, ", threads_per_chain, " thread(s) per chain, ", placement))
              }
              selected

      }
}

##
## ---- TABLE, panel A ---------------------------------------------------------------------------------------------------------------
##
## Selection: local-HPC, N = 50,000, one thread per chain, unpinned, from the DRAM-traffic runs of design (i) (label prefix E5,
## cases_HPC_E5.csv), since only these runs have memory-controller counters; throughput and IPC are taken from the same runs, so all
## four columns of a row describe one run. Rows: BayesMVP (MD_BayesMVP) with 1 and 100 chunks, the Stan model (AD_Stan, 1 chunk) and
## the Stan model with tape chunking (AD_Stan_tape_chunked, 250 chunks); N_chains = 8, 96 and 180.
## Columns: throughput (1 decimal), IPC (2 decimals), DRAM traffic = reads + writes at the memory controllers in GB/s (0 decimals),
## DRAM reads at the memory controllers per chain-iteration in GB (2 decimals).
##
{
      panel_A_rows <-  data.frame( implementation = c("BayesMVP, 1 chunk", "BayesMVP, 100 chunks", "Stan model, 1 chunk", "Stan model + tape chunking, 250 chunks"),
                                   algorithm      = c("MD_BayesMVP", "MD_BayesMVP", "AD_Stan", "AD_Stan_tape_chunked"),
                                   chunks         = c(1, 100, 1, 250),
                                   stringsAsFactors = FALSE)
      panel_A_columns <-  data.frame( column = c("Throughput", "IPC", "DRAM traffic (GB/s)", "DRAM reads per chain-iteration (GB)"),
                                      metric = c("throughput", "IPC", "DRAM_total_GB_per_second", "DRAM_read_GB_per_chain_iteration"),
                                      digits = c(1, 2, 0, 2),
                                      stringsAsFactors = FALSE)
      panel_A_chains <-  c(8, 96, 180)
}

##
## ---- TABLE, panel B ---------------------------------------------------------------------------------------------------------------
##
## Selection: design (ii) (label prefix E2), local-HPC, N = 50,000, 8 chains, one thread per chain; "One CCD" = placement oneCCD
## (CPUs 0-7), "Eight CCDs" = placement spread8CCD (CPUs 0,8,16,...,56). Rows: BayesMVP with 1 and 500 chunks, the Stan model
## (1 chunk) and the Stan model with tape chunking (250 chunks). Values: throughput (1 decimal).
##
{
      panel_B_rows <-  data.frame( implementation = c("BayesMVP, 1 chunk", "BayesMVP, 500 chunks", "Stan model, 1 chunk", "Stan model + tape chunking, 250 chunks"),
                                   algorithm      = c("MD_BayesMVP", "MD_BayesMVP", "AD_Stan", "AD_Stan_tape_chunked"),
                                   chunks         = c(1, 500, 1, 250),
                                   stringsAsFactors = FALSE)
      panel_B_columns <-  data.frame( column    = c("One CCD", "Eight CCDs"),
                                      placement = c("oneCCD", "spread8CCD"),
                                      stringsAsFactors = FALSE)
}

##
## ---- TABLE, panel C ---------------------------------------------------------------------------------------------------------------
##
## Selection, N = 50,000, throughput (1 decimal):
##     "8 chains x 1 thread (chunked)": local-HPC, design (ii) (E2), 8 chains on one CCD (oneCCD) or eight CCDs (spread8CCD), with
##         BayesMVP 500 chunks / Stan model tape chunking 250 chunks; laptop, design (i) (E1, unpinned), 8 chains, with BayesMVP
##         100 chunks / Stan model tape chunking 500 chunks (the chunk counts stated in the table caption).
##     "1 chain x 8 threads (chunking + WCP)": design (iii) (E3), MD_BayesMVP_WCP with 200 chunks / AD_Stan_WCP with 250 chunks,
##         1 chain x 8 threads; local-HPC oneCCD or spread8CCD (as the row), laptop unpinned.
##     "1 chain x 8 threads (WCP-only)": design (iii) (E3), the same algorithms with N_chunks = 8 (one chunk per thread).
##
{
      panel_C_rows <-  utils::read.table( text = "
            implementation                    | device | placement  | chunked_experiment | chunked_algorithm    | chunked_chunks | WCP_algorithm   | joint_chunks
            BayesMVP, local-HPC, one CCD      | HPC    | oneCCD     | E2                 | MD_BayesMVP          | 500            | MD_BayesMVP_WCP | 200
            BayesMVP, local-HPC, eight CCDs   | HPC    | spread8CCD | E2                 | MD_BayesMVP          | 500            | MD_BayesMVP_WCP | 200
            Stan model, local-HPC, one CCD    | HPC    | oneCCD     | E2                 | AD_Stan_tape_chunked | 250            | AD_Stan_WCP     | 250
            Stan model, local-HPC, eight CCDs | HPC    | spread8CCD | E2                 | AD_Stan_tape_chunked | 250            | AD_Stan_WCP     | 250
            BayesMVP, laptop                  | Laptop | unpinned   | E1                 | MD_BayesMVP          | 100            | MD_BayesMVP_WCP | 200
            Stan model, laptop                | Laptop | unpinned   | E1                 | AD_Stan_tape_chunked | 500            | AD_Stan_WCP     | 250
            ", sep = "|", header = TRUE, strip.white = TRUE, stringsAsFactors = FALSE)
      WCP_only_chunks <-  8
}

##
## ---- Table cells: reproduce, round as printed and compare with the paper -----------------------------------------------------------
##
{
      fn_cell <-  function( panel,
                            implementation,
                            column,
                            digits,
                            value_in_paper,
                            case,
                            metric
      ) {

              if (length(value_in_paper) != 1) stop(paste0("No single paper value for panel ", panel, ", ", implementation, ", ", column, "."))
              value <-  case[[metric]]
              printed <-  unname(formatC(value, format = "f", digits = digits))
              data.frame( panel             = panel,
                          implementation    = implementation,
                          device            = case$device,
                          chains            = case$n_chains,
                          column            = column,
                          value_unrounded   = value,
                          value_as_printed  = printed,
                          value_in_paper    = value_in_paper,
                          match             = isTRUE(printed == value_in_paper),
                          chunks            = case$num_chunks,
                          threads_per_chain = case$threads_per_chain,
                          placement         = case$placement,
                          source_case_label = case$label,
                          stringsAsFactors  = FALSE)

      }
      ##
      table_cells <-  list()
      ##
      ## Panel A:
      for (row_index in seq_len(nrow(panel_A_rows))) {
            for (chains in panel_A_chains) {
                  case <-  fn_case( device = "HPC", experiment = "E5", algorithm = panel_A_rows$algorithm[row_index],
                                    chunks = panel_A_rows$chunks[row_index], chains = chains)
                  paper_row <-  paper_panel_A[paper_panel_A$implementation == panel_A_rows$implementation[row_index] & paper_panel_A$chains == chains, , drop = FALSE]
                  for (column_index in seq_len(nrow(panel_A_columns))) {
                        table_cells[[length(table_cells) + 1]] <-  fn_cell( panel = "A", implementation = panel_A_rows$implementation[row_index],
                                                                            column = panel_A_columns$column[column_index], digits = panel_A_columns$digits[column_index],
                                                                            value_in_paper = paper_row[[panel_A_columns$column[column_index]]],
                                                                            case = case, metric = panel_A_columns$metric[column_index])
                  }
            }
      }
      ##
      ## Panel B:
      for (row_index in seq_len(nrow(panel_B_rows))) {
            paper_row <-  paper_panel_B[paper_panel_B$implementation == panel_B_rows$implementation[row_index], , drop = FALSE]
            for (column_index in seq_len(nrow(panel_B_columns))) {
                  case <-  fn_case( device = "HPC", experiment = "E2", algorithm = panel_B_rows$algorithm[row_index],
                                    chunks = panel_B_rows$chunks[row_index], chains = 8, placement = panel_B_columns$placement[column_index])
                  table_cells[[length(table_cells) + 1]] <-  fn_cell( panel = "B", implementation = panel_B_rows$implementation[row_index],
                                                                      column = panel_B_columns$column[column_index], digits = 1,
                                                                      value_in_paper = paper_row[[panel_B_columns$column[column_index]]],
                                                                      case = case, metric = "throughput")
            }
      }
      ##
      ## Panel C:
      for (row_index in seq_len(nrow(panel_C_rows))) {
            selection <-  panel_C_rows[row_index, , drop = FALSE]
            paper_row <-  paper_panel_C[paper_panel_C$implementation == selection$implementation, , drop = FALSE]
            panel_C_cases <-  list( "8 chains x 1 thread (chunked)"        = fn_case( device = selection$device, experiment = selection$chunked_experiment,
                                                                                      algorithm = selection$chunked_algorithm, chunks = selection$chunked_chunks,
                                                                                      chains = 8, placement = selection$placement),
                                    "1 chain x 8 threads (chunking + WCP)" = fn_case( device = selection$device, experiment = "E3", algorithm = selection$WCP_algorithm,
                                                                                      chunks = selection$joint_chunks, chains = 1, threads_per_chain = 8,
                                                                                      placement = selection$placement),
                                    "1 chain x 8 threads (WCP-only)"       = fn_case( device = selection$device, experiment = "E3", algorithm = selection$WCP_algorithm,
                                                                                      chunks = WCP_only_chunks, chains = 1, threads_per_chain = 8,
                                                                                      placement = selection$placement))
            for (column in names(panel_C_cases)) {
                  table_cells[[length(table_cells) + 1]] <-  fn_cell( panel = "C", implementation = selection$implementation, column = column, digits = 1,
                                                                      value_in_paper = paper_row[[column]], case = panel_C_cases[[column]], metric = "throughput")
            }
      }
      ##
      table_cells <-  do.call(what = rbind, args = table_cells)
      ##
      message(fn_colour(paste0("\n---- Table table:exp4_profiling: ", nrow(table_cells), " cells ----"), "cyan"))
      for (cell_index in seq_len(nrow(table_cells))) {
            cell <-  table_cells[cell_index, , drop = FALSE]
            message(fn_colour(paste0( if (cell$match) "    MATCH     " else "    MISMATCH  ",
                                      "panel ", cell$panel, " | ", cell$implementation, " | ", cell$device, " | ", cell$chains, " chain(s) | ", cell$column,
                                      ": reproduced ", trimws(formatC(cell$value_unrounded, digits = 6, format = "fg")), " -> ", cell$value_as_printed,
                                      "; paper ", cell$value_in_paper, "  [", cell$source_case_label, "]"),
                              if (cell$match) "green" else "red"))
      }
}

##
## ---- Values quoted in the text: reproduce --------------------------------------------------------------------------------------------
##
## Each entry gives the reproduced value (before rounding), the case labels it was computed from, and any detail needed to audit it.
##
{
      fn_reproduced <-  function( value,
                                  sources = character(0),
                                  detail  = ""
      ) {

              list(value = value, sources = paste(sources, collapse = "; "), detail = detail)

      }
      fn_ratio_detail <-  function(numerator, denominator) paste0(trimws(formatC(numerator, digits = 6, format = "fg")), " / ", trimws(formatC(denominator, digits = 6, format = "fg")))
      fn_percent_change <-  function(new, old) 100 * (new / old - 1)
      ##
      ## ---- Design, case lists and iteration counts
      MD_rows <-  grepl(pattern = "^MD", x = times$algorithm)
      design_i <-  case_measures[case_measures$experiment == "E1", , drop = FALSE]
      idle_row <-  cal_umc[cal_umc$label == "idle_3s", , drop = FALSE]
      idle_GB_per_second <-  (idle_row$DRAM_read_bytes + idle_row$DRAM_write_bytes) / idle_row$seconds / 1e9
      ##
      ## ---- Unchunked plateau (DRAM traffic, E5, 48 to 180 chains)
      plateau_cases <-  case_measures[ case_measures$experiment == "E5" & case_measures$num_chunks == 1 & case_measures$n_chains >= 48 &
                                       case_measures$algorithm %in% c("MD_BayesMVP", "AD_Stan"), , drop = FALSE]
      Stan_plateau_cases <-  plateau_cases[plateau_cases$algorithm == "AD_Stan", , drop = FALSE]
      peak_GB_per_second <-  HPC_memory_channels * HPC_transfers_per_second * bytes_per_transfer / 1e9
      ##
      ## ---- Panel A cases (E5) used in the text
      E5 <-  function(algorithm, chunks, chains) fn_case(device = "HPC", experiment = "E5", algorithm = algorithm, chunks = chunks, chains = chains)
      B1_8 <-  E5("MD_BayesMVP", 1, 8);             B1_48 <-  E5("MD_BayesMVP", 1, 48)
      B1_96 <-  E5("MD_BayesMVP", 1, 96);           B1_180 <-  E5("MD_BayesMVP", 1, 180)
      B25_8 <-  E5("MD_BayesMVP", 25, 8);           B25_48 <-  E5("MD_BayesMVP", 25, 48)
      B25_96 <-  E5("MD_BayesMVP", 25, 96);         B25_180 <-  E5("MD_BayesMVP", 25, 180)
      B100_8 <-  E5("MD_BayesMVP", 100, 8);         B100_48 <-  E5("MD_BayesMVP", 100, 48)
      B100_96 <-  E5("MD_BayesMVP", 100, 96);       B100_180 <-  E5("MD_BayesMVP", 100, 180)
      B500_96 <-  E5("MD_BayesMVP", 500, 96);       B500_180 <-  E5("MD_BayesMVP", 500, 180)
      S1_8 <-  E5("AD_Stan", 1, 8);                 S1_96 <-  E5("AD_Stan", 1, 96);                  S1_180 <-  E5("AD_Stan", 1, 180)
      T250_96 <-  E5("AD_Stan_tape_chunked", 250, 96)
      ##
      ## BayesMVP DRAM reads per chain-iteration, 1 chunk vs 100 chunks (the table's chunked configuration), at every E5 chain count:
      E5_chains <-  sort(unique(case_measures$n_chains[case_measures$experiment == "E5"]))
      reads_reduction_100 <-  sapply(X = E5_chains, FUN = function(chains) {
            E5("MD_BayesMVP", 1, chains)$DRAM_read_GB_per_chain_iteration / E5("MD_BayesMVP", 100, chains)$DRAM_read_GB_per_chain_iteration
      })
      ## The same over every chunked BayesMVP configuration of the DRAM-traffic runs (25, 100 and 500 chunks), for information:
      E5_BayesMVP_chunked <-  case_measures[case_measures$experiment == "E5" & case_measures$algorithm == "MD_BayesMVP" & case_measures$num_chunks > 1, , drop = FALSE]
      reads_reduction_all <-  sapply(X = seq_len(nrow(E5_BayesMVP_chunked)), FUN = function(row_index) {
            E5("MD_BayesMVP", 1, E5_BayesMVP_chunked$n_chains[row_index])$DRAM_read_GB_per_chain_iteration / E5_BayesMVP_chunked$DRAM_read_GB_per_chain_iteration[row_index]
      })
      reads_reduction_all_max <-  E5_BayesMVP_chunked[which.max(reads_reduction_all), , drop = FALSE]
      ##
      ## ---- One chain (design (i), E1)
      E1_HPC <-  function(algorithm, chunks, chains) fn_case(device = "HPC", experiment = "E1", algorithm = algorithm, chunks = chunks, chains = chains)
      one_chain_B1 <-  E1_HPC("MD_BayesMVP", 1, 1)
      one_chain_B100 <-  E1_HPC("MD_BayesMVP", 100, 1)
      one_chain_B25 <-  E1_HPC("MD_BayesMVP", 25, 1)
      one_chain_speedups <-  sapply(X = c(4, 25, 100, 500), FUN = function(chunks) E1_HPC("MD_BayesMVP", chunks, 1)$throughput / one_chain_B1$throughput)
      one_chain_S1 <-  E1_HPC("AD_Stan", 1, 1)
      one_chain_T250 <-  E1_HPC("AD_Stan_tape_chunked", 250, 1)
      ##
      ## ---- Design (ii) (E2) and design (iv) (E4)
      E2 <-  function(algorithm, chunks, placement) fn_case(device = "HPC", experiment = "E2", algorithm = algorithm, chunks = chunks, chains = 8, placement = placement)
      E4_HPC <-  function(algorithm, chunks) fn_case(device = "HPC", experiment = "E4", algorithm = algorithm, chunks = chunks, chains = 16, placement = "SMT_oneCCD")
      E1_Laptop <-  function(algorithm, chunks, chains) fn_case(device = "Laptop", experiment = "E1", algorithm = algorithm, chunks = chunks, chains = chains)
      E4_Laptop <-  function(algorithm, chunks) fn_case(device = "Laptop", experiment = "E4", algorithm = algorithm, chunks = chunks, chains = 8, placement = "physical_cores_only")
      ##
      SMT_configurations <-  data.frame( key = c("500", "25", "1", "Stan_tape250"), algorithm = c("MD_BayesMVP", "MD_BayesMVP", "MD_BayesMVP", "AD_Stan_tape_chunked"),
                                         chunks = c(500, 25, 1, 250), stringsAsFactors = FALSE)
      SMT_HPC <-  lapply(X = seq_len(nrow(SMT_configurations)), FUN = function(i) {
            list(new = E4_HPC(SMT_configurations$algorithm[i], SMT_configurations$chunks[i]), old = E2(SMT_configurations$algorithm[i], SMT_configurations$chunks[i], "oneCCD"))
      })
      SMT_Laptop <-  lapply(X = seq_len(nrow(SMT_configurations)), FUN = function(i) {
            list(new = E1_Laptop(SMT_configurations$algorithm[i], SMT_configurations$chunks[i], 16), old = E4_Laptop(SMT_configurations$algorithm[i], SMT_configurations$chunks[i]))
      })
      names(SMT_HPC) <-  SMT_configurations$key
      names(SMT_Laptop) <-  SMT_configurations$key
      ##
      ## ---- Panel C ratios (from the reproduced table cells)
      fn_panel_C <-  function(implementation, column) table_cells[table_cells$panel == "C" & table_cells$implementation == implementation & table_cells$column == column, , drop = FALSE]
      chunked_column <-  "8 chains x 1 thread (chunked)"
      joint_column <-  "1 chain x 8 threads (chunking + WCP)"
      WCP_only_column <-  "1 chain x 8 threads (WCP-only)"
      fn_panel_C_ratio <-  function(implementation, column) fn_panel_C(implementation, chunked_column)$value_unrounded / fn_panel_C(implementation, column)$value_unrounded
      fn_panel_C_sources <-  function(implementation, column) c(fn_panel_C(implementation, chunked_column)$source_case_label, fn_panel_C(implementation, column)$source_case_label)
      laptop_ratios <-  c( fn_panel_C_ratio("BayesMVP, laptop", joint_column), fn_panel_C_ratio("BayesMVP, laptop", WCP_only_column),
                           fn_panel_C_ratio("Stan model, laptop", joint_column), fn_panel_C_ratio("Stan model, laptop", WCP_only_column))
      joint_spread_slowdown <-  c( BayesMVP = 100 * (1 - fn_panel_C("BayesMVP, local-HPC, eight CCDs", joint_column)$value_unrounded / fn_panel_C("BayesMVP, local-HPC, one CCD", joint_column)$value_unrounded),
                                   Stan     = 100 * (1 - fn_panel_C("Stan model, local-HPC, eight CCDs", joint_column)$value_unrounded / fn_panel_C("Stan model, local-HPC, one CCD", joint_column)$value_unrounded))
      WCP_only_faster_spread <-  c( BayesMVP = fn_panel_C("BayesMVP, local-HPC, eight CCDs", WCP_only_column)$value_unrounded > fn_panel_C("BayesMVP, local-HPC, one CCD", WCP_only_column)$value_unrounded,
                                    Stan     = fn_panel_C("Stan model, local-HPC, eight CCDs", WCP_only_column)$value_unrounded > fn_panel_C("Stan model, local-HPC, one CCD", WCP_only_column)$value_unrounded)
      ##
      ## ---- Working sets and the L3 cache per active thread (threads spread evenly over the L3 caches)
      working_set_bytes <-  function(chunks) bytes_per_individual * N_profiling / chunks
      L3_per_active_thread <-  function(device, chains) L3_bytes[[device]] / ceiling(chains / number_of_L3[[device]])
      exceeded_cases <-  unique(design_i[design_i$N == N_profiling & design_i$n_chains > 1, c("device", "n_chains")])
      exceeded_cases <-  exceeded_cases[order(exceeded_cases$device, exceeded_cases$n_chains), , drop = FALSE]
      exceeded_cases <-  exceeded_cases[ mapply(FUN = function(device, chains) working_set_bytes(25) > L3_per_active_thread(device, chains),
                                                exceeded_cases$device, exceeded_cases$n_chains), , drop = FALSE]
      ##
      ## ---- 25 vs 100 chunks where the 25-chunk working set exceeds the L3 cache per active thread (design (i), E1, core counters)
      auto_pairs <-  lapply(X = seq_len(nrow(exceeded_cases)), FUN = function(i) {
            list( chunks_25  = fn_case(device = exceeded_cases$device[i], experiment = "E1", algorithm = "MD_BayesMVP", chunks = 25, chains = exceeded_cases$n_chains[i]),
                  chunks_100 = fn_case(device = exceeded_cases$device[i], experiment = "E1", algorithm = "MD_BayesMVP", chunks = 100, chains = exceeded_cases$n_chains[i]))
      })
      auto_RAM_ratios <-  sapply(X = auto_pairs, FUN = function(pair) pair$chunks_25$RAM_fills_MB_per_chain_iteration / pair$chunks_100$RAM_fills_MB_per_chain_iteration)
      auto_throughput_ratios <-  sapply(X = auto_pairs, FUN = function(pair) pair$chunks_25$throughput / pair$chunks_100$throughput)
      auto_sources <-  unlist(lapply(X = auto_pairs, FUN = function(pair) c(pair$chunks_25$label, pair$chunks_100$label)))
      auto_case_names <-  paste(exceeded_cases$device, exceeded_cases$n_chains, "chains")
      auto_HPC180 <-  list(chunks_25 = E1_HPC("MD_BayesMVP", 25, 180), chunks_100 = E1_HPC("MD_BayesMVP", 100, 180))
      auto_100_vs_500 <-  sapply(X = c(96, 180), FUN = function(chains) {
            E1_HPC("MD_BayesMVP", 100, chains)$RAM_fills_MB_per_chain_iteration / E1_HPC("MD_BayesMVP", 500, chains)$RAM_fills_MB_per_chain_iteration
      })
      close_to_threshold <-  1.05   ## "close to (or below)": at most 5% above the 500-chunk value
      ##
      ## ---- Stan model: DRAM traffic (reads + writes) per chain-iteration, 1 chunk vs tape chunking with 250 chunks, at 96 chains
      Stan_traffic_reduction_96 <-  S1_96$DRAM_total_GB_per_chain_iteration / T250_96$DRAM_total_GB_per_chain_iteration
      Stan_reads_reduction_96 <-  S1_96$DRAM_read_GB_per_chain_iteration / T250_96$DRAM_read_GB_per_chain_iteration
      ##
      ## ---- Fraction of the DRAM reads per chain-iteration removed by 25 chunks (vs 1 chunk), E5, at up to 96 chains
      removed_25 <-  c( "8"  = 1 - B25_8$DRAM_read_GB_per_chain_iteration / B1_8$DRAM_read_GB_per_chain_iteration,
                        "48" = 1 - B25_48$DRAM_read_GB_per_chain_iteration / B1_48$DRAM_read_GB_per_chain_iteration,
                        "96" = 1 - B25_96$DRAM_read_GB_per_chain_iteration / B1_96$DRAM_read_GB_per_chain_iteration)
      ##
      reproduced_text <-  list(
            ##
            ## Plan for experiment 4:
            plan_idle_background_GB_per_s       = fn_reproduced(idle_GB_per_second, "cal_umc.csv: idle_3s", paste0("(reads + writes) / seconds = (", idle_row$DRAM_read_bytes, " + ", idle_row$DRAM_write_bytes, ") / ", idle_row$seconds)),
            plan_n_iter_short                   = fn_reproduced(paste(sort(unique(times$n_iter_short)), collapse = ", "), "mechanism_times_HPC.csv; mechanism_times_Laptop.csv"),
            plan_n_iter_long_BayesMVP_N50000    = fn_reproduced(paste(sort(unique(times$n_iter_long[MD_rows & times$N == 50000])), collapse = ", "), "all MD_BayesMVP* cases, N = 50,000"),
            plan_n_iter_long_BayesMVP_N10000    = fn_reproduced(paste(sort(unique(times$n_iter_long[MD_rows & times$N == 10000])), collapse = ", "), "all MD_BayesMVP* cases, N = 10,000"),
            plan_n_iter_long_Stan_N50000        = fn_reproduced(paste(sort(unique(times$n_iter_long[!MD_rows & times$N == 50000])), collapse = ", "), "all AD_Stan* cases, N = 50,000"),
            plan_n_iter_long_Stan_N10000        = fn_reproduced(paste(sort(unique(times$n_iter_long[!MD_rows & times$N == 10000])), collapse = ", "), "all AD_Stan* cases, N = 10,000"),
            plan_design_i_N_values              = fn_reproduced(paste(sort(unique(design_i$N)), collapse = ", "), "E1 cases"),
            plan_design_i_chains_HPC            = fn_reproduced(paste(sort(unique(design_i$n_chains[design_i$device == "HPC"])), collapse = ", "), "E1 cases, HPC"),
            plan_design_i_chains_Laptop         = fn_reproduced(paste(sort(unique(design_i$n_chains[design_i$device == "Laptop"])), collapse = ", "), "E1 cases, Laptop"),
            ##
            ## Supplement:
            supp_cases_HPC                      = fn_reproduced(as.character(nrow(cases_HPC)), "cases_HPC.csv"),
            supp_cases_Laptop                   = fn_reproduced(as.character(nrow(cases_Laptop)), "cases_Laptop.csv"),
            supp_cases_HPC_DRAM                 = fn_reproduced(as.character(nrow(cases_HPC_DRAM)), "cases_HPC_E5.csv"),
            ##
            ## E4 results: memory bandwidth without chunking:
            res_plateau_lower_GB_per_s          = fn_reproduced(min(plateau_cases$DRAM_total_GB_per_second), plateau_cases$label[which.min(plateau_cases$DRAM_total_GB_per_second)], "minimum over BayesMVP and the Stan model, 1 chunk, 48/96/180 chains"),
            res_plateau_upper_GB_per_s          = fn_reproduced(max(plateau_cases$DRAM_total_GB_per_second), plateau_cases$label[which.max(plateau_cases$DRAM_total_GB_per_second)], "maximum over BayesMVP and the Stan model, 1 chunk, 48/96/180 chains"),
            res_peak_GB_per_s                   = fn_reproduced(peak_GB_per_second, detail = "12 channels x 4800 MT/s x 8 bytes"),
            res_IPC_BayesMVP_1chunk_8           = fn_reproduced(B1_8$IPC, B1_8$label),
            res_IPC_BayesMVP_1chunk_180         = fn_reproduced(B1_180$IPC, B1_180$label),
            res_IPC_Stan_1chunk_8               = fn_reproduced(S1_8$IPC, S1_8$label),
            res_IPC_Stan_1chunk_180             = fn_reproduced(S1_180$IPC, S1_180$label),
            res_throughput_BayesMVP_1chunk_48   = fn_reproduced(B1_48$throughput, B1_48$label),
            res_throughput_BayesMVP_1chunk_180  = fn_reproduced(B1_180$throughput, B1_180$label),
            ##
            ## E4 results: DRAM reads and IPC with chunking:
            res_reads_BayesMVP_1chunk_96        = fn_reproduced(B1_96$DRAM_read_GB_per_chain_iteration, B1_96$label),
            res_reads_BayesMVP_100chunks_96     = fn_reproduced(B100_96$DRAM_read_GB_per_chain_iteration, B100_96$label),
            res_reads_BayesMVP_1chunk_180       = fn_reproduced(B1_180$DRAM_read_GB_per_chain_iteration, B1_180$label),
            res_reads_BayesMVP_100chunks_180    = fn_reproduced(B100_180$DRAM_read_GB_per_chain_iteration, B100_180$label),
            res_reads_Stan_1chunk_96            = fn_reproduced(S1_96$DRAM_read_GB_per_chain_iteration, S1_96$label),
            res_reads_Stan_tape250_96           = fn_reproduced(T250_96$DRAM_read_GB_per_chain_iteration, T250_96$label),
            res_IPC_BayesMVP_100chunks_96       = fn_reproduced(B100_96$IPC, B100_96$label),
            res_IPC_Stan_tape250_96             = fn_reproduced(T250_96$IPC, T250_96$label),
            res_throughput_100chunks_increases  = fn_reproduced(all(diff(c(B100_8$throughput, B100_48$throughput, B100_96$throughput, B100_180$throughput)) > 0),
                                                                c(B100_8$label, B100_48$label, B100_96$label, B100_180$label),
                                                                paste0("8/48/96/180 chains: ", paste(formatC(c(B100_8$throughput, B100_48$throughput, B100_96$throughput, B100_180$throughput), format = "f", digits = 1), collapse = ", "))),
            res_one_chain_speedup_BayesMVP      = fn_reproduced(one_chain_B100$throughput / one_chain_B1$throughput, c(one_chain_B100$label, one_chain_B1$label),
                                                                paste0("100 chunks (the table's chunked configuration); speed-ups with 4/25/100/500 chunks: ", paste(formatC(one_chain_speedups, format = "f", digits = 2), collapse = ", "))),
            res_one_chain_Stan_tape250          = fn_reproduced(one_chain_T250$throughput, one_chain_T250$label),
            res_one_chain_Stan_1chunk           = fn_reproduced(one_chain_S1$throughput, one_chain_S1$label),
            ##
            ## E4 results: packing (design (ii)):
            res_packing_slowdown_BayesMVP       = fn_reproduced(E2("MD_BayesMVP", 1, "spread8CCD")$throughput / E2("MD_BayesMVP", 1, "oneCCD")$throughput,
                                                                c(E2("MD_BayesMVP", 1, "spread8CCD")$label, E2("MD_BayesMVP", 1, "oneCCD")$label), "eight CCDs / one CCD"),
            res_packing_slowdown_Stan           = fn_reproduced(E2("AD_Stan", 1, "spread8CCD")$throughput / E2("AD_Stan", 1, "oneCCD")$throughput,
                                                                c(E2("AD_Stan", 1, "spread8CCD")$label, E2("AD_Stan", 1, "oneCCD")$label), "eight CCDs / one CCD"),
            res_packed_BayesMVP_500             = fn_reproduced(E2("MD_BayesMVP", 500, "oneCCD")$throughput, E2("MD_BayesMVP", 500, "oneCCD")$label),
            res_spread_BayesMVP_500             = fn_reproduced(E2("MD_BayesMVP", 500, "spread8CCD")$throughput, E2("MD_BayesMVP", 500, "spread8CCD")$label),
            res_packed_Stan_tape250             = fn_reproduced(E2("AD_Stan_tape_chunked", 250, "oneCCD")$throughput, E2("AD_Stan_tape_chunked", 250, "oneCCD")$label),
            res_spread_Stan_tape250             = fn_reproduced(E2("AD_Stan_tape_chunked", 250, "spread8CCD")$throughput, E2("AD_Stan_tape_chunked", 250, "spread8CCD")$label),
            ##
            ## E4 results: working set per chunk:
            res_working_set_KB_per_individual   = fn_reproduced(bytes_per_individual / 1e3, detail = "1,608 bytes / 1000"),
            res_working_set_MB_1chunk           = fn_reproduced(working_set_bytes(1) / 1e6, detail = "1,608 bytes x 50,000 / 1 chunk / 1e6"),
            res_working_set_MB_25chunks         = fn_reproduced(working_set_bytes(25) / 1e6, detail = "1,608 bytes x 50,000 / 25 chunks / 1e6"),
            res_working_set_MB_100chunks        = fn_reproduced(working_set_bytes(100) / 1e6, detail = "1,608 bytes x 50,000 / 100 chunks / 1e6"),
            res_working_set_MB_500chunks        = fn_reproduced(working_set_bytes(500) / 1e6, detail = "1,608 bytes x 50,000 / 500 chunks / 1e6"),
            res_25chunks_fit_L3_not_L2          = fn_reproduced(working_set_bytes(25) < L3_bytes[["HPC"]] && working_set_bytes(25) > L2_bytes_HPC,
                                                                detail = paste0(working_set_bytes(25), " bytes vs L3 ", L3_bytes[["HPC"]], " and L2 ", L2_bytes_HPC, " bytes")),
            res_100_500chunks_fit_L2            = fn_reproduced(working_set_bytes(100) < L2_bytes_HPC && working_set_bytes(500) < L2_bytes_HPC,
                                                                detail = paste0(working_set_bytes(100), " and ", working_set_bytes(500), " bytes vs L2 ", L2_bytes_HPC, " bytes")),
            ##
            ## E4 results: 25 chunks under SMT, and the SMT gains from 96 to 180 chains (E5):
            res_25chunks_most_traffic_removed   = fn_reproduced(all(removed_25 > 0.5), c(B25_8$label, B25_48$label, B25_96$label, B1_8$label, B1_48$label, B1_96$label),
                                                                paste0("share of the DRAM reads per chain-iteration removed at 8/48/96 chains: ", paste(formatC(removed_25, format = "f", digits = 2), collapse = ", "))),
            res_DRAM_traffic_25chunks_96        = fn_reproduced(B25_96$DRAM_total_GB_per_second, B25_96$label),
            res_DRAM_traffic_25chunks_180       = fn_reproduced(B25_180$DRAM_total_GB_per_second, B25_180$label),
            res_IPC_25chunks_96                 = fn_reproduced(B25_96$IPC, B25_96$label),
            res_IPC_25chunks_180                = fn_reproduced(B25_180$IPC, B25_180$label),
            res_SMT_gain_100chunks              = fn_reproduced(fn_percent_change(B100_180$throughput, B100_96$throughput), c(B100_180$label, B100_96$label), fn_ratio_detail(B100_180$throughput, B100_96$throughput)),
            res_SMT_gain_500chunks              = fn_reproduced(fn_percent_change(B500_180$throughput, B500_96$throughput), c(B500_180$label, B500_96$label), fn_ratio_detail(B500_180$throughput, B500_96$throughput)),
            res_SMT_gain_1chunk                 = fn_reproduced(fn_percent_change(B1_180$throughput, B1_96$throughput), c(B1_180$label, B1_96$label), fn_ratio_detail(B1_180$throughput, B1_96$throughput))
      )
      ##
      ## E4 results: two chains per physical core (design (iv)), local-HPC and laptop. "16chains" = two chains per core (local-HPC:
      ## E4 SMT_oneCCD; laptop: E1 16 chains, unpinned); "8chains" = one chain per core (local-HPC: E2 oneCCD; laptop: E4 physical_cores_only):
      for (key in SMT_configurations$key) {
            for (device in c("HPC", "Laptop")) {
                  pair <-  if (device == "HPC") SMT_HPC[[key]] else SMT_Laptop[[key]]
                  reproduced_text[[paste0("res_SMT_", device, "_change_", key)]] <-  fn_reproduced( fn_percent_change(pair$new$throughput, pair$old$throughput),
                                                                                                    c(pair$new$label, pair$old$label),
                                                                                                    fn_ratio_detail(pair$new$throughput, pair$old$throughput))
                  reproduced_text[[paste0("res_SMT_", device, "_16chains_", key)]] <-  fn_reproduced(pair$new$throughput, pair$new$label)
                  reproduced_text[[paste0("res_SMT_", device, "_8chains_", key)]] <-  fn_reproduced(pair$old$throughput, pair$old$label)
            }
      }
      ##
      reproduced_text <-  c(reproduced_text, list(
            ##
            ## E4 results: eight chunked chains vs one chain with 8 threads (panel C):
            res_WCP_oneCCD_joint_BayesMVP       = fn_reproduced(fn_panel_C_ratio("BayesMVP, local-HPC, one CCD", joint_column), fn_panel_C_sources("BayesMVP, local-HPC, one CCD", joint_column)),
            res_WCP_oneCCD_joint_Stan           = fn_reproduced(fn_panel_C_ratio("Stan model, local-HPC, one CCD", joint_column), fn_panel_C_sources("Stan model, local-HPC, one CCD", joint_column)),
            res_WCP_oneCCD_WCPonly_BayesMVP     = fn_reproduced(fn_panel_C_ratio("BayesMVP, local-HPC, one CCD", WCP_only_column), fn_panel_C_sources("BayesMVP, local-HPC, one CCD", WCP_only_column)),
            res_WCP_oneCCD_WCPonly_Stan         = fn_reproduced(fn_panel_C_ratio("Stan model, local-HPC, one CCD", WCP_only_column), fn_panel_C_sources("Stan model, local-HPC, one CCD", WCP_only_column)),
            res_WCP_eightCCD_joint_BayesMVP     = fn_reproduced(fn_panel_C_ratio("BayesMVP, local-HPC, eight CCDs", joint_column), fn_panel_C_sources("BayesMVP, local-HPC, eight CCDs", joint_column)),
            res_WCP_eightCCD_joint_Stan         = fn_reproduced(fn_panel_C_ratio("Stan model, local-HPC, eight CCDs", joint_column), fn_panel_C_sources("Stan model, local-HPC, eight CCDs", joint_column)),
            res_WCP_eightCCD_WCPonly_BayesMVP   = fn_reproduced(fn_panel_C_ratio("BayesMVP, local-HPC, eight CCDs", WCP_only_column), fn_panel_C_sources("BayesMVP, local-HPC, eight CCDs", WCP_only_column)),
            res_WCP_eightCCD_WCPonly_Stan       = fn_reproduced(fn_panel_C_ratio("Stan model, local-HPC, eight CCDs", WCP_only_column), fn_panel_C_sources("Stan model, local-HPC, eight CCDs", WCP_only_column)),
            res_WCP_Laptop_lower                = fn_reproduced(min(laptop_ratios), "panel C laptop rows", paste0("laptop ratios (BayesMVP vs chunking + WCP, vs WCP-only; Stan model, the same): ", paste(formatC(laptop_ratios, format = "f", digits = 3), collapse = ", "))),
            res_WCP_Laptop_upper                = fn_reproduced(max(laptop_ratios), "panel C laptop rows", paste0("laptop ratios (BayesMVP vs chunking + WCP, vs WCP-only; Stan model, the same): ", paste(formatC(laptop_ratios, format = "f", digits = 3), collapse = ", "))),
            res_WCP_only_faster_when_spread     = fn_reproduced(all(WCP_only_faster_spread), "panel C, WCP-only column, local-HPC rows",
                                                                paste0("BayesMVP ", fn_panel_C("BayesMVP, local-HPC, eight CCDs", WCP_only_column)$value_as_printed, " vs ", fn_panel_C("BayesMVP, local-HPC, one CCD", WCP_only_column)$value_as_printed,
                                                                       "; Stan model ", fn_panel_C("Stan model, local-HPC, eight CCDs", WCP_only_column)$value_as_printed, " vs ", fn_panel_C("Stan model, local-HPC, one CCD", WCP_only_column)$value_as_printed)),
            res_joint_slower_spread_lower       = fn_reproduced(min(joint_spread_slowdown), "panel C, chunking + WCP column, local-HPC rows", paste0("BayesMVP ", formatC(joint_spread_slowdown[["BayesMVP"]], format = "f", digits = 2), "%, Stan model ", formatC(joint_spread_slowdown[["Stan"]], format = "f", digits = 2), "%")),
            res_joint_slower_spread_upper       = fn_reproduced(max(joint_spread_slowdown), "panel C, chunking + WCP column, local-HPC rows", paste0("BayesMVP ", formatC(joint_spread_slowdown[["BayesMVP"]], format = "f", digits = 2), "%, Stan model ", formatC(joint_spread_slowdown[["Stan"]], format = "f", digits = 2), "%")),
            ##
            ## E1 section on the automatic N_chunks rule:
            auto_working_set_MB_25chunks        = fn_reproduced(working_set_bytes(25) / 1e6, detail = "1,608 bytes x 50,000 / 25 chunks / 1e6"),
            auto_25chunks_exceeds_L3_cases      = fn_reproduced(paste(paste(exceeded_cases$device, exceeded_cases$n_chains), collapse = "; "), "design (i) chain counts at N = 50,000",
                                                                "25-chunk working set > L3 cache / ceiling(N_chains / N_L3)"),
            auto_RAM_ratio_lower                = fn_reproduced(min(auto_RAM_ratios), auto_sources, paste0("RAM fills per chain-iteration, 25 / 100 chunks: ", paste(paste0(auto_case_names, " ", formatC(auto_RAM_ratios, format = "f", digits = 2)), collapse = ", "))),
            auto_RAM_ratio_upper                = fn_reproduced(max(auto_RAM_ratios), auto_sources, paste0("RAM fills per chain-iteration, 25 / 100 chunks: ", paste(paste0(auto_case_names, " ", formatC(auto_RAM_ratios, format = "f", digits = 2)), collapse = ", "))),
            auto_throughput_ratio_lower         = fn_reproduced(min(auto_throughput_ratios), auto_sources, paste0("throughput, 25 / 100 chunks: ", paste(paste0(auto_case_names, " ", formatC(auto_throughput_ratios, format = "f", digits = 3)), collapse = ", "))),
            auto_throughput_ratio_upper         = fn_reproduced(max(auto_throughput_ratios), auto_sources, paste0("throughput, 25 / 100 chunks: ", paste(paste0(auto_case_names, " ", formatC(auto_throughput_ratios, format = "f", digits = 3)), collapse = ", "))),
            auto_RAM_MB_HPC180_25chunks         = fn_reproduced(auto_HPC180$chunks_25$RAM_fills_MB_per_chain_iteration, auto_HPC180$chunks_25$label),
            auto_RAM_MB_HPC180_100chunks        = fn_reproduced(auto_HPC180$chunks_100$RAM_fills_MB_per_chain_iteration, auto_HPC180$chunks_100$label),
            auto_throughput_ratio_HPC180        = fn_reproduced(auto_HPC180$chunks_25$throughput / auto_HPC180$chunks_100$throughput, c(auto_HPC180$chunks_25$label, auto_HPC180$chunks_100$label)),
            auto_L3_fraction_100chunks_96       = fn_reproduced(working_set_bytes(100) / L3_per_active_thread("HPC", 96), detail = "0.804 MB / (2^25 bytes / 8 threads per CCD)"),
            auto_L3_fraction_100chunks_180      = fn_reproduced(working_set_bytes(100) / L3_per_active_thread("HPC", 180), detail = "0.804 MB / (2^25 bytes / 15 threads per CCD)"),
            auto_RAM_100_close_to_500           = fn_reproduced(all(auto_100_vs_500 <= close_to_threshold), c(E1_HPC("MD_BayesMVP", 100, 96)$label, E1_HPC("MD_BayesMVP", 500, 96)$label, E1_HPC("MD_BayesMVP", 100, 180)$label, E1_HPC("MD_BayesMVP", 500, 180)$label),
                                                                paste0("RAM fills per chain-iteration, 100 / 500 chunks at 96 and 180 chains: ", paste(formatC(auto_100_vs_500, format = "f", digits = 3), collapse = ", "), " (close = at most ", close_to_threshold, ")")),
            auto_one_chain_25_vs_100_percent    = fn_reproduced(100 * (1 - one_chain_B25$throughput / one_chain_B100$throughput), c(one_chain_B25$label, one_chain_B100$label), fn_ratio_detail(one_chain_B25$throughput, one_chain_B100$throughput)),
            auto_25chunks_L3_fraction           = fn_reproduced(working_set_bytes(25) / L3_bytes[["HPC"]], detail = "3.216 MB / 2^25 bytes"),
            auto_f_range_lower                  = fn_reproduced(working_set_bytes(100) / L3_per_active_thread("HPC", 96), detail = "as auto_L3_fraction_100chunks_96"),
            auto_f_range_upper                  = fn_reproduced(working_set_bytes(100) / L3_per_active_thread("HPC", 180), detail = "as auto_L3_fraction_100chunks_180"),
            ##
            ## Discussion and conclusion:
            disc_plateau_lower_GB_per_s         = fn_reproduced(min(plateau_cases$DRAM_total_GB_per_second), plateau_cases$label[which.min(plateau_cases$DRAM_total_GB_per_second)], "as res_plateau_lower_GB_per_s"),
            disc_plateau_upper_GB_per_s         = fn_reproduced(max(plateau_cases$DRAM_total_GB_per_second), plateau_cases$label[which.max(plateau_cases$DRAM_total_GB_per_second)], "as res_plateau_upper_GB_per_s"),
            disc_peak_GB_per_s                  = fn_reproduced(peak_GB_per_second, detail = "12 channels x 4800 MT/s x 8 bytes"),
            disc_BayesMVP_reads_reduction       = fn_reproduced(max(reads_reduction_100), paste0("E5 MD_BayesMVP, 1 vs 100 chunks, ", E5_chains[which.max(reads_reduction_100)], " chains"),
                                                                paste0("1 chunk / 100 chunks at 8/48/96/180 chains: ", paste(formatC(reads_reduction_100, format = "f", digits = 2), collapse = ", "))),
            disc_Stan_plateau_lower_GB_per_s    = fn_reproduced(min(Stan_plateau_cases$DRAM_total_GB_per_second), Stan_plateau_cases$label[which.min(Stan_plateau_cases$DRAM_total_GB_per_second)]),
            disc_Stan_plateau_upper_GB_per_s    = fn_reproduced(max(Stan_plateau_cases$DRAM_total_GB_per_second), Stan_plateau_cases$label[which.max(Stan_plateau_cases$DRAM_total_GB_per_second)]),
            disc_Stan_tape_traffic_reduction_96 = fn_reproduced(Stan_traffic_reduction_96, c(S1_96$label, T250_96$label),
                                                                paste0("reads + writes per chain-iteration; reads only: ", formatC(Stan_reads_reduction_96, format = "f", digits = 2))),
            concl_BayesMVP_reads_reduction      = fn_reproduced(max(reads_reduction_100), paste0("E5 MD_BayesMVP, 1 vs 100 chunks, ", E5_chains[which.max(reads_reduction_100)], " chains"), "as disc_BayesMVP_reads_reduction")
      ))
}

##
## ---- Values quoted in the text: compare with the paper -------------------------------------------------------------------------------
##
{
      missing_ids <-  setdiff(paper_text$id, names(reproduced_text))
      if (length(missing_ids) > 0) stop(paste0("No reproduced value for: ", paste(missing_ids, collapse = ", ")))
      ##
      text_checks <-  do.call(what = rbind, args = lapply(X = seq_len(nrow(paper_text)), FUN = function(row_index) {

              claim <-  paper_text[row_index, , drop = FALSE]
              reproduced <-  reproduced_text[[claim$id]]
              value <-  reproduced$value
              ##
              if (claim$rule %in% c("rounded", "approximate")) {
                    printed <-  unname(formatC(value, format = "f", digits = claim$digits))
                    match <-  isTRUE(printed == claim$value_in_paper) ||
                              (claim$rule == "approximate" && abs(value / as.numeric(claim$value_in_paper) - 1) <= 0.05)
              } else if (claim$rule == "signed") {
                    printed <-  unname(formatC(value, format = "f", digits = claim$digits, flag = "+"))
                    match <-  isTRUE(printed == claim$value_in_paper)
              } else {
                    printed <-  unname(as.character(value))
                    match <-  isTRUE(printed == claim$value_in_paper)
              }
              ##
              data.frame( claim,
                          value_unrounded  = if (is.numeric(value)) value else NA,
                          value_as_printed = printed,
                          match            = match,
                          sources          = reproduced$sources,
                          detail           = reproduced$detail,
                          stringsAsFactors = FALSE)

      }))
      ##
      message(fn_colour(paste0("\n---- Experiment-4 values quoted in the text: ", nrow(text_checks), " values ----"), "cyan"))
      for (location in unique(text_checks$location)) {
            message(fn_colour(paste0("  ", location), "cyan"))
            location_checks <-  text_checks[text_checks$location == location, , drop = FALSE]
            for (check_index in seq_len(nrow(location_checks))) {
                  check <-  location_checks[check_index, , drop = FALSE]
                  unrounded <-  if (is.na(check$value_unrounded)) "" else paste0(trimws(formatC(check$value_unrounded, digits = 6, format = "fg")), " -> ")
                  message(fn_colour(paste0( if (check$match) "    MATCH     " else "    MISMATCH  ",
                                            "\"", check$quoted, "\": reproduced ", unrounded, check$value_as_printed, "; paper ", check$value_in_paper,
                                            " (", check$rule, if (check$source == "derived") ", derived from stated constants" else "", ")",
                                            if (nzchar(check$detail)) paste0("\n                  ", check$detail) else "",
                                            if (nzchar(check$sources)) paste0("\n                  [", check$sources, "]") else ""),
                                    if (check$match) "green" else "red"))
            }
      }
      ##
      ## Related values, printed for information only (not counted):
      message(fn_colour(paste0( "\n  For information (not counted): over every chunked BayesMVP configuration of the DRAM-traffic runs (25, 100 and 500 chunks), ",
                                "the largest reduction in DRAM reads per chain-iteration is ", formatC(max(reads_reduction_all), format = "f", digits = 1), "x (",
                                reads_reduction_all_max$num_chunks, " chunks, ", reads_reduction_all_max$n_chains, " chains); with 100 chunks it is ",
                                formatC(max(reads_reduction_100), format = "f", digits = 1), "x."), "cyan"))
}

##
## ---- Outputs: CSV of the table cells, and the tabular in the paper's layout ------------------------------------------------------------
##
{
      utils::write.csv(x = table_cells, file = file.path(mechanism_dir, "table_exp4_profiling_reproduced.csv"), row.names = FALSE)
      ##
      fn_printed <-  function(panel, implementation, column, chains = NULL) {
            rows <-  table_cells$panel == panel & table_cells$implementation == implementation & table_cells$column == column
            if (!is.null(chains)) rows <-  rows & table_cells$chains == chains
            table_cells$value_as_printed[rows]
      }
      fn_row_label <-  function(implementation) {
            if (implementation == "Stan model + tape chunking, 250 chunks") return(r"(\makecell[l]{Stan model + tape \\ chunking, 250 chunks})")
            implementation
      }
      fn_pad <-  function(text, width) formatC(text, width = -width)
      ##
      tex_lines <-  c( r"(\begin{tabular}{llrrrr})",
                       r"(\hline)",
                       r"(\multicolumn{6}{l}{\textit{Panel A}} \\)",
                       r"(Implementation & $N_{\text{chains}}$ & Throughput & IPC & \makecell[r]{DRAM traffic \\ (GB/s)} & \makecell[r]{DRAM reads per \\ chain-iteration (GB)} \\)",
                       r"(\hline)")
      for (implementation in panel_A_rows$implementation) {
            for (chains in panel_A_chains) {
                  first_column <-  if (chains == panel_A_chains[1]) fn_row_label(implementation) else ""
                  values <-  sapply(X = panel_A_columns$column, FUN = function(column) fn_printed("A", implementation, column, chains))
                  tex_lines <-  c(tex_lines, paste0( fn_pad(first_column, 24), " & ", fn_pad(chains, 3), " & ", fn_pad(values[1], 5), " & ", fn_pad(values[2], 4),
                                                     " & ", fn_pad(values[3], 3), " & ", values[4], r"( \\)"))
            }
      }
      tex_lines <-  c( tex_lines,
                       "%%%%",
                       r"(\hline)",
                       "%%%%",
                       r"(\multicolumn{6}{l}{\textit{Panel B}} \\)",
                       r"( & & \makecell[r]{One CCD} & \makecell[r]{Eight CCDs} & & \\)",
                       "%%%%",
                       r"(\hline)",
                       "%%%%")
      for (implementation in panel_B_rows$implementation) {
            tex_lines <-  c(tex_lines, paste0( fn_pad(fn_row_label(implementation), 38), " & 8 & ", fn_pad(fn_printed("B", implementation, "One CCD"), 4),
                                               " & ", fn_printed("B", implementation, "Eight CCDs"), r"( & & \\)"))
      }
      tex_lines <-  c( tex_lines,
                       "%%%%",
                       r"(\hline)",
                       "%%%%",
                       r"(\multicolumn{6}{l}{\textit{Panel C}} \\)",
                       r"( & & \makecell[r]{8 chains $\times$ \\ 1 thread (chunked)} & \makecell[r]{1 chain $\times$ 8 threads \\ (chunking + WCP)} & \makecell[r]{1 chain $\times$ 8 threads \\ (WCP-only)} & \\)",
                       "%%%%",
                       r"(\hline)",
                       "%%%%")
      for (implementation in panel_C_rows$implementation) {
            tex_lines <-  c(tex_lines, paste0( fn_pad(implementation, 33), " & & ", fn_pad(fn_printed("C", implementation, chunked_column), 4), " & ",
                                               fn_pad(fn_printed("C", implementation, joint_column), 4), " & ", fn_pad(fn_printed("C", implementation, WCP_only_column), 4), r"( & \\)"))
      }
      tex_lines <-  c( tex_lines,
                       "%%%%",
                       r"(\hline)",
                       "%%%%",
                       r"(\end{tabular})")
      writeLines(text = tex_lines, con = file.path(mechanism_dir, "table_exp4_profiling_reproduced.tex"))
      ##
      message(fn_colour(paste0("\nWritten: ", file.path(mechanism_dir, "table_exp4_profiling_reproduced.csv"), "\n         ",
                               file.path(mechanism_dir, "table_exp4_profiling_reproduced.tex")), "cyan"))
}

##
## ---- Summary ---------------------------------------------------------------------------------------------------------------------------
##
{
      fn_count_line <-  function(label, matches) {
            message(fn_colour(paste0( "    ", label, ": ", sum(matches), " of ", length(matches), " match",
                                      if (all(matches)) "" else paste0(", ", sum(!matches), " MISMATCH")), if (all(matches)) "green" else "red"))
      }
      ##
      message(fn_colour("\n---- Summary ----", "cyan"))
      fn_count_line("Table table:exp4_profiling (cells)", table_cells$match)
      fn_count_line("Values in the text, reproduced from the CSVs", text_checks$match[text_checks$source == "CSV"])
      fn_count_line("Values in the text, derived from constants stated in the paper", text_checks$match[text_checks$source == "derived"])
      fn_count_line("All values", c(table_cells$match, text_checks$match))
}























