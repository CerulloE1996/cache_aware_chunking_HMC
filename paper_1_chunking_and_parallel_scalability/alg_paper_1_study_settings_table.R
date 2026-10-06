##
## ===============================================================================================================
## alg_paper_1_study_settings_table.R
##
## Builds the Paper 1 study-settings table (E1 Parts I-V and E2) from the settings in the runner
## alg_paper_1_chunking_WCP_par_scaling.R, for both devices (local-HPC and laptop).
## The runner is PARSED, not sourced: only its settings expressions are evaluated (with source() and the Mplus
## runner functions replaced by no-op stubs), so no benchmark is started. The benchmark grid itself is then built with
## the helper fn_paper1_benchmark_grid(), so every cell of the table comes from the cases the runner would measure.
## Output: report_outputs/tables/table_paper1_study_settings.tex
##
{
      algorithm_study_dir <-  path.expand("~/Documents/Work/PhD_work/Alg_paper_analysis")
      paper1_dir          <-  file.path(algorithm_study_dir, "paper_1_chunking_and_parallel_scalability")
      runner_file         <-  file.path(paper1_dir, "alg_paper_1_chunking_WCP_par_scaling.R")
      helper_file         <-  file.path(paper1_dir, "R_fns_alg_paper_1_chunking_WCP_par_scaling.R")
      output_dir          <-  file.path(paper1_dir, "report_outputs", "tables")
      output_file         <-  file.path(output_dir, "table_paper1_study_settings.tex")
}
##
## ---- Evaluate the runner's settings expressions for one device ------------------------------------------------
##
fn_paper1_runner_settings_for_device <-  function( runner_file,
                                                   device_wanted
) {

        runner_expressions <-  parse(file = runner_file, keep.source = FALSE)
        settings_env <-  new.env(parent = globalenv())
        ##
        ## ---- No-op stubs: the runner's source() calls and Mplus function handles are not needed for the settings.
        settings_env$source <-  function(...) invisible(x = NULL)
        settings_env$fn_paper1_run_Mplus        <-  function(...) NULL
        settings_env$fn_paper1_verify_Mplus_run <-  function(...) NULL
        ##
        for (runner_expression in as.list(x = runner_expressions)) {

              expression_text <-  paste0(deparse(expr = runner_expression), collapse = "\n")
              ## Stop before the first call that would build or run the benchmark.
              if (grepl(pattern = "fn_run_paper1_benchmark", x = expression_text, fixed = TRUE)) break
              ##
              suppressMessages(expr = eval(expr = runner_expression, envir = settings_env))
              ##
              ## The runner picks the device from the core count; replace it with the wanted device after it is assigned.
              if (exists(x = "device", envir = settings_env, inherits = FALSE)) {
                    settings_env$device <-  device_wanted
                    if (exists(x = "paper1_settings", envir = settings_env, inherits = FALSE)) {
                          settings_env$paper1_settings$device <-  device_wanted
                    }
              }

        }
        ##
        return(settings_env$paper1_settings)

}
##
## ---- Build both devices' settings and benchmark grids (the helper file only defines functions) ----------------
##
{
      helper_env <-  new.env(parent = globalenv())
      source(file = helper_file, local = helper_env)
      ##
      settings_by_device <-  list()
      grid_by_device     <-  list()
      ##
      for (device_name in c("HPC", "Laptop")) {

            message(paste0("\033[36mReading the runner settings for device = ", device_name, "\033[0m"))
            settings_by_device[[device_name]] <-  fn_paper1_runner_settings_for_device( runner_file   = runner_file,
                                                                                        device_wanted = device_name)
            grid_by_device[[device_name]] <-  suppressMessages(expr = helper_env$fn_paper1_benchmark_grid(
                                                                       settings = settings_by_device[[device_name]]))
            message(paste0("\033[36m  grid rows = ", nrow(x = grid_by_device[[device_name]]), "\033[0m"))

      }
}
##
## ---- Formatting helpers ---------------------------------------------------------------------------------------
##
fn_format_count <-  function( x ) {

        ## 10000 -> "10,000" (as in the paper); numbers below 10,000 are written without a comma.
        return(ifelse(test = x >= 10000, yes = formatC(x = x, format = "d", big.mark = ","), no = as.character(x = x)))

}
##
fn_format_set <-  function( x ) {

        x <-  sort(x = unique(x = x))
        if (length(x = x) == 1) return(paste0("$", fn_format_count(x), "$"))
        return(paste0("$\\{", paste0(fn_format_count(x), collapse = ",\\; "), "\\}$"))

}
##
fn_format_set_wrapped <-  function( x,
                                    values_per_line = 5
) {

        ## Long candidate sets are split over two or more lines so that the column stays narrow.
        x <-  sort(x = unique(x = x))
        if (length(x = x) <= values_per_line + 1) return(fn_format_set(x = x))
        line_groups <-  split(x = fn_format_count(x), f = ceiling(x = seq_along(along.with = x) / values_per_line))
        line_texts <-  sapply(X = line_groups, FUN = paste0, collapse = ",\\; ")
        n_lines <-  length(x = line_texts)
        line_texts <-  paste0(ifelse(test = seq_len(length.out = n_lines) == 1, yes = "\\{", no = ""),
                              line_texts,
                              ifelse(test = seq_len(length.out = n_lines) == n_lines, yes = "\\}", no = ","))
        return(paste0("\\makecell[l]{", paste0("$", line_texts, "$", collapse = "\\\\"), "}"))

}
##
fn_format_power_of_ten <-  function( x ) {

        exponent <-  log10(x = x)
        if (abs(x = exponent - round(x = exponent)) < 1e-12) return(paste0("$10^{", round(x = exponent), "}$"))
        return(paste0("$", format(x = x, scientific = TRUE), "$"))

}
##
fn_unique_or_stop <-  function( x,
                                what
) {

        x <-  unique(x = x)
        if (length(x = x) != 1) stop(paste0("Expected one value for ", what, ", found: ", paste0(x, collapse = ", ")))
        return(x)

}
##
## ---- Per-N iteration and chunk cells from the grids -----------------------------------------------------------
##
{
      many_chains_minimum <-  fn_unique_or_stop( x = c(settings_by_device$HPC$wcp_many_chains_minimum,
                                                       settings_by_device$Laptop$wcp_many_chains_minimum),
                                                 what = "wcp_many_chains_minimum")
      N_values <-  fn_unique_or_stop(x = list(settings_by_device$HPC$N_vec), what = "N_vec")[[1]]
      ##
      all_grid <-  do.call(what = rbind,
                           args = lapply(X = grid_by_device,
                                         FUN = function(g) g[, c("device", "algorithm", "N", "num_chunks", "n_chains",
                                                                 "threads_per_chain", "n_iter", "n_iter_short_run",
                                                                 "run", "mplus_iteration_mode",
                                                                 "mplus_iteration_mode_short_run")]))
      all_grid$WCP_only <-  all_grid$algorithm %in% c("MD_BayesMVP_WCP", "AD_Stan_WCP") &
                            all_grid$num_chunks == all_grid$threads_per_chain
      ##
      check_list <-  character(0)
      ##
      ## N_iter of one arm (set of algorithms) at one N, over both devices; must be unique.
      fn_n_iter_cell <-  function( algorithms,
                                   N_value,
                                   chain_filter,
                                   WCP_only_filter = NULL,
                                   label
      ) {

              rows <-  all_grid[all_grid$algorithm %in% algorithms & all_grid$N == N_value & chain_filter(all_grid$n_chains), ]
              if (!is.null(x = WCP_only_filter)) rows <-  rows[rows$WCP_only == WCP_only_filter, ]
              if (!nrow(x = rows)) return(NA)
              value <-  fn_unique_or_stop(x = rows$n_iter, what = paste0(label, " at N = ", N_value))
              check_list <<-  c(check_list, paste0(label, ", N = ", N_value, ": N_iter = ", value,
                                                   " (", nrow(x = rows), " grid rows, devices: ",
                                                   paste0(unique(x = rows$device), collapse = "+"), ")"))
              return(value)

      }
      ##
      few_chains  <-  function(n_chains) n_chains <  many_chains_minimum
      many_chains <-  function(n_chains) n_chains >= many_chains_minimum
      any_chains  <-  function(n_chains) rep(x = TRUE, times = length(x = n_chains))
      ##
      implementations <-  list(
            list( label         = "\\makecell[l]{NicoStan+\\\\BayesMVP}",
                  plain_label   = "NicoStan+BayesMVP",
                  standard      = c("MD_BayesMVP"),
                  WCP           = "MD_BayesMVP_WCP",
                  chunked       = c("MD_BayesMVP", "MD_BayesMVP_WCP"),
                  serial_extra  = "MD_BayesMVP"),
            list( label         = "\\makecell[l]{Stan model\\\\(via NicoStan)}",
                  plain_label   = "Stan model (via NicoStan)",
                  standard      = c("AD_Stan", "AD_Stan_chunked", "AD_Stan_tape_chunked"),
                  WCP           = "AD_Stan_WCP",
                  chunked       = c("AD_Stan_chunked", "AD_Stan_tape_chunked", "AD_Stan_WCP"),
                  serial_extra  = "AD_Stan_tape_chunked"),
            list( label         = "Mplus",
                  plain_label   = "Mplus",
                  standard      = c("Mplus_standard"),
                  WCP           = "Mplus_WCP",
                  chunked       = NULL,
                  serial_extra  = NULL))
      ##
      panel_A_lines <-  character(0)
      ##
      for (N_value in N_values) {

            N_key <-  as.character(x = N_value)
            for (implementation_index in seq_along(along.with = implementations)) {

                  imp <-  implementations[[implementation_index]]
                  ##
                  ## ---- Candidate N_chunks: the runner's chunk grids (identical on both devices; checked against the grid).
                  chunk_cell <-  "--"
                  serial_cell <-  "--"
                  if (!is.null(x = imp$chunked)) {

                        chunk_sets <-  lapply(X = imp$chunked, FUN = function(a) {
                              lapply(X = settings_by_device, FUN = function(s) s$chunks_by_algorithm[[a]][[N_key]]) })
                        chunk_set <-  fn_unique_or_stop(x = unlist(x = lapply(X = chunk_sets, FUN = function(l) {
                                          lapply(X = l, FUN = function(v) paste0(sort(x = v), collapse = ","))
                                      })), what = paste0(imp$plain_label, " chunk grid at N = ", N_key))
                        chunk_values <-  as.numeric(x = strsplit(x = chunk_set, split = ",")[[1]])
                        ## Check: the chunked standard arm's one-chain rows contain exactly these chunk counts (plus extras).
                        chunk_cell <-  fn_format_set_wrapped(x = chunk_values)
                        check_list <-  c(check_list, paste0(imp$plain_label, ", N = ", N_key, ": candidate N_chunks = ",
                                                            paste0(chunk_values, collapse = ", ")))
                        ##
                        ## ---- Extra one-chain chunk counts (matched-serial references), per device from the grid.
                        serial_by_device <-  sapply(X = names(grid_by_device), FUN = function(d) {
                              g <-  grid_by_device[[d]]
                              extra <-  setdiff(x = g$num_chunks[g$algorithm == imp$serial_extra & g$N == N_value],
                                                y = chunk_values)
                              if (length(x = extra)) paste0(fn_format_count(sort(x = extra)), collapse = ", ") else "--"
                        })
                        if (any(serial_by_device != "--")) {
                              serial_cell <-  paste0(serial_by_device[["HPC"]], " / ", serial_by_device[["Laptop"]])
                        }
                        check_list <-  c(check_list, paste0(imp$plain_label, ", N = ", N_key,
                                                            ": extra one-chain N_chunks (HPC / laptop) = ", serial_cell))

                  }
                  ##
                  standard_iter <-  fn_n_iter_cell( algorithms = imp$standard, N_value = N_value, chain_filter = any_chains,
                                                    label = paste0(imp$plain_label, " standard"))
                  WCP_few_iter  <-  fn_n_iter_cell( algorithms = imp$WCP, N_value = N_value, chain_filter = few_chains,
                                                    label = paste0(imp$plain_label, " WCP, N_chains < ", many_chains_minimum))
                  ##
                  if (is.null(x = imp$chunked)) {

                        WCP_many_iter <-  fn_n_iter_cell( algorithms = imp$WCP, N_value = N_value, chain_filter = many_chains,
                                                          label = paste0(imp$plain_label, " WCP, N_chains >= ", many_chains_minimum))
                        WCP_many_cell <-  fn_format_count(WCP_many_iter)

                  } else {

                        WCP_many_chunked <-  fn_n_iter_cell( algorithms = imp$WCP, N_value = N_value,
                                                             chain_filter = many_chains, WCP_only_filter = FALSE,
                                                             label = paste0(imp$plain_label, " chunking + WCP, N_chains >= ",
                                                                            many_chains_minimum))
                        WCP_many_only    <-  fn_n_iter_cell( algorithms = imp$WCP, N_value = N_value,
                                                             chain_filter = many_chains, WCP_only_filter = TRUE,
                                                             label = paste0(imp$plain_label, " WCP-only, N_chains >= ",
                                                                            many_chains_minimum))
                        WCP_many_cell <-  fn_format_count(WCP_many_chunked)
                        if (!is.na(x = WCP_many_only) && WCP_many_only != WCP_many_chunked) {
                              WCP_many_cell <-  paste0(WCP_many_cell, " (", fn_format_count(WCP_many_only), ")")
                        }
                        ## WCP-only with N_chains < many_chains_minimum must equal the chunking + WCP count (no brackets there).
                        WCP_few_only <-  fn_n_iter_cell( algorithms = imp$WCP, N_value = N_value, chain_filter = few_chains,
                                                         WCP_only_filter = TRUE,
                                                         label = paste0(imp$plain_label, " WCP-only, N_chains < ",
                                                                        many_chains_minimum))
                        if (WCP_few_only != WCP_few_iter) stop("WCP-only N_iter differs at N_chains < many-chain minimum.")

                  }
                  ##
                  N_cell <-  if (implementation_index == 1) {
                        paste0("\\multirow{", length(x = implementations), "}{*}{$", fn_format_count(N_value), "$}")
                  } else ""
                  ##
                  panel_A_lines <-  c(panel_A_lines,
                                      paste0(N_cell, " & ", imp$label, " & ", chunk_cell, " & ", serial_cell, " & ",
                                             fn_format_count(standard_iter), " & ", fn_format_count(WCP_few_iter),
                                             " & ", WCP_many_cell, " \\\\"))

            }
            panel_A_lines <-  c(panel_A_lines, if (N_value != utils::tail(x = N_values, n = 1)) "\\hline" else NULL)

      }
}
##
## ---- Per-device settings (chain grids, WCP allocations) from the grids ----------------------------------------
##
{
      fn_allocation_text <-  function( g,
                                       algorithm,
                                       N_value
      ) {

              rows <-  unique(x = g[g$algorithm == algorithm & g$N == N_value, c("n_chains", "threads_per_chain")])
              chain_values <-  sort(x = unique(x = rows$n_chains))
              pieces <-  sapply(X = chain_values, FUN = function(n_chains) {
                    paste0(n_chains, " \\times ", fn_format_set(x = rows$threads_per_chain[rows$n_chains == n_chains]))
              })
              pieces <-  gsub(pattern = "\\$", replacement = "", x = pieces)
              return(pieces)

      }
      ##
      fn_group_N_by_allocations <-  function( g,
                                              algorithm
      ) {

              allocation_by_N <-  lapply(X = N_values, FUN = function(N_value) fn_allocation_text(g, algorithm, N_value))
              keys <-  sapply(X = allocation_by_N, FUN = paste0, collapse = "; ")
              groups <-  split(x = N_values, f = factor(x = keys, levels = unique(x = keys)))
              return(lapply(X = names(groups), FUN = function(k) list(N = groups[[k]],
                                                                      pieces = allocation_by_N[[match(x = k, table = keys)]])))

      }
      ##
      ## Stan and BayesMVP WCP allocations must be the same (the paper states the candidates are shared); Mplus is a subset.
      for (device_name in names(grid_by_device)) {
            g <-  grid_by_device[[device_name]]
            for (N_value in N_values) {
                  a_B <-  fn_allocation_text(g, "MD_BayesMVP_WCP", N_value)
                  a_S <-  fn_allocation_text(g, "AD_Stan_WCP", N_value)
                  if (!identical(x = a_B, y = a_S)) stop(paste0("BayesMVP and Stan WCP allocations differ: ", device_name,
                                                                ", N = ", N_value))
            }
      }
      ##
      fn_N_range_text <-  function( N_group ) {

              if (identical(x = N_group, y = N_values)) return("all $N$")
              return(paste0("$N \\in ", fn_format_set(x = N_group) |> gsub(pattern = "\\$", replacement = ""), "$"))

      }
      ##
      fn_allocation_cell <-  function( device_name ) {

              groups <-  fn_group_N_by_allocations(g = grid_by_device[[device_name]], algorithm = "MD_BayesMVP_WCP")
              lines <-  character(0)
              for (grp in groups) {
                    ## Greedy wrap at about 60 characters of LaTeX per line keeps the cell narrow.
                    line_index <-  integer(0)
                    current_line <-  1
                    current_width <-  0
                    for (piece in grp$pieces) {
                          if (current_width > 0 && current_width + nchar(x = piece) > 60) {
                                current_line <-  current_line + 1
                                current_width <-  0
                          }
                          line_index <-  c(line_index, current_line)
                          current_width <-  current_width + nchar(x = piece)
                    }
                    piece_lines <-  split(x = grp$pieces, f = line_index)
                    piece_lines <-  sapply(X = piece_lines, FUN = function(p) paste0("$", paste0(p, collapse = ";\\; "), "$"))
                    lines <-  c(lines, paste0(fn_N_range_text(grp$N), ":"), piece_lines)
                    check_list <<-  c(check_list, paste0(device_name, " WCP allocations (BayesMVP = Stan), N = ",
                                                         paste0(grp$N, collapse = ", "), ": ",
                                                         paste0(grp$pieces, collapse = "; ")))
              }
              return(paste0("\\makecell[l]{", paste0(lines, collapse = "\\\\"), "}"))

      }
      ##
      fn_mplus_allocation_cell <-  function( device_name ) {

              g <-  grid_by_device[[device_name]]
              omitted_by_N <-  lapply(X = N_values, FUN = function(N_value) {
                    B <-  unique(x = g[g$algorithm == "MD_BayesMVP_WCP" & g$N == N_value, c("n_chains", "threads_per_chain")])
                    M <-  unique(x = g[g$algorithm == "Mplus_WCP"       & g$N == N_value, c("n_chains", "threads_per_chain")])
                    B_key <-  paste0(B$n_chains, " \\times ", B$threads_per_chain)
                    M_key <-  paste0(M$n_chains, " \\times ", M$threads_per_chain)
                    if (length(x = setdiff(x = M_key, y = B_key))) stop("Mplus allocation not in the BayesMVP grid.")
                    setdiff(x = B_key, y = M_key)
              })
              omitted <-  fn_unique_or_stop(x = sapply(X = omitted_by_N, FUN = paste0, collapse = ";\\; "),
                                            what = paste0("Mplus omitted allocations, ", device_name))
              check_list <<-  c(check_list, paste0(device_name, " Mplus_WCP allocations = BayesMVP ones without: ",
                                                   if (nzchar(x = omitted)) omitted else "none"))
              if (!nzchar(x = omitted)) return("as above")
              return(paste0("\\makecell[l]{as above, without\\\\$", omitted, "$}"))

      }
      ##
      fn_chain_cell <-  function( device_name ) {

              g <-  grid_by_device[[device_name]]
              standard_algorithms <-  c("MD_BayesMVP", "AD_Stan", "AD_Stan_chunked", "AD_Stan_tape_chunked", "Mplus_standard")
              chain_sets <-  sapply(X = standard_algorithms, FUN = function(a) {
                    paste0(sort(x = unique(x = g$n_chains[g$algorithm == a])), collapse = ",")
              })
              chain_set <-  fn_unique_or_stop(x = chain_sets, what = paste0("standard-arm chain grid, ", device_name))
              tpc <-  fn_unique_or_stop(x = g$threads_per_chain[g$algorithm %in% standard_algorithms],
                                        what = "standard-arm threads per chain")
              if (tpc != 1) stop("Standard arms must use one thread per chain.")
              chain_values <-  as.numeric(x = strsplit(x = chain_set, split = ",")[[1]])
              check_list <<-  c(check_list, paste0(device_name, " standard-arm N_chains = ", chain_set,
                                                   " (N_threads/chain = 1)"))
              return(fn_format_set(x = chain_values))

      }
      ##
      chain_cells      <-  sapply(X = names(grid_by_device), FUN = fn_chain_cell)
      allocation_cells <-  sapply(X = names(grid_by_device), FUN = fn_allocation_cell)
      mplus_cells      <-  sapply(X = names(grid_by_device), FUN = fn_mplus_allocation_cell)
}
##
## ---- Device-independent settings ------------------------------------------------------------------------------
##
{
      s_HPC <-  settings_by_device$HPC
      s_Lap <-  settings_by_device$Laptop
      for (setting_name in c("bayesmvp", "stan_via_NicoStan", "n_runs", "mplus_n_runs", "mplus_iteration_mode",
                             "mplus_iteration_mode_short_run", "untimed_warm_up_run_before_timing", "timing_method",
                             "stan_backend")) {
            if (!identical(x = s_HPC[[setting_name]], y = s_Lap[[setting_name]])) stop(paste0(setting_name, " differs by device."))
      }
      if (!identical(x = s_HPC$stan_backend, y = "NicoStan")) stop("The table describes the NicoStan Stan backend only.")
      ##
      overhead_iterations <-  unique(x = c(all_grid$n_iter_short_run))
      overhead_iterations <-  fn_unique_or_stop(x = overhead_iterations, what = "overhead-run iterations")
      ##
      L_bayesmvp <-  s_HPC$bayesmvp$L_main
      L_stan     <-  s_HPC$stan_via_NicoStan$fixed_L
      eps_bayesmvp <-  s_HPC$bayesmvp$eps_main
      eps_stan     <-  s_HPC$stan_via_NicoStan$step_size
      metric_name <-  c(dense = "dense", diag = "diagonal")
      metric_bayesmvp <-  metric_name[[s_HPC$bayesmvp$metric_shape_main]]
      metric_stan     <-  metric_name[[s_HPC$stan_via_NicoStan$metric_shape_main]]
      ##
      repeats_BS    <-  fn_unique_or_stop(x = tapply(X = all_grid$run, INDEX = !grepl(pattern = "^Mplus_", x = all_grid$algorithm),
                                                     FUN = max)[["TRUE"]], what = "repeats")
      repeats_Mplus <-  max(all_grid$run[grepl(pattern = "^Mplus_", x = all_grid$algorithm)])
      mplus_modes   <-  sort(x = unique(x = all_grid$mplus_iteration_mode[!is.na(x = all_grid$mplus_iteration_mode)]))
      mplus_mode_short <-  fn_unique_or_stop(x = all_grid$mplus_iteration_mode_short_run[
                                                  !is.na(x = all_grid$mplus_iteration_mode_short_run)],
                                             what = "Mplus overhead-run mode")
      ##
      check_list <-  c(check_list,
                       paste0("L: BayesMVP L_main = ", L_bayesmvp, ", Stan fixed_L = ", L_stan),
                       paste0("step size: BayesMVP eps_main = ", eps_bayesmvp, ", Stan step_size = ", eps_stan),
                       paste0("metric: BayesMVP ", s_HPC$bayesmvp$metric_shape_main, ", Stan ",
                              s_HPC$stan_via_NicoStan$metric_shape_main),
                       paste0("randomize_tau: BayesMVP ", s_HPC$bayesmvp$randomize_tau, ", Stan ",
                              s_HPC$stan_via_NicoStan$randomize_tau, "; diffusion_HMC ", s_HPC$bayesmvp$diffusion_HMC,
                              "; partitioned_HMC ", s_HPC$bayesmvp$partitioned_HMC),
                       paste0("repeats (grid max run): BayesMVP/Stan ", repeats_BS, " (n_runs = ", s_HPC$n_runs,
                              "), Mplus ", repeats_Mplus, " (mplus_n_runs = ", s_HPC$mplus_n_runs, ")"),
                       paste0("overhead-run iterations (grid n_iter_short_run) = ", overhead_iterations),
                       paste0("untimed warm-up before timing = ", s_HPC$untimed_warm_up_run_before_timing,
                              " (1 iteration, set in fn_paper1_time_case_with_timing_method)"),
                       paste0("Mplus timed-run modes = ", paste0(mplus_modes, collapse = ", "),
                              "; overhead-run mode = ", mplus_mode_short),
                       paste0("Stan backend = ", s_HPC$stan_backend, "; timing method = ", s_HPC$timing_method))
      if (repeats_BS != s_HPC$n_runs || repeats_Mplus != s_HPC$mplus_n_runs) stop("Repeat counts in grid differ from runner.")
      if (!isTRUE(x = s_HPC$untimed_warm_up_run_before_timing)) stop("Caption assumes the untimed warm-up call is used.")
}
##
## ---- Assemble the LaTeX table ---------------------------------------------------------------------------------
##
{
      mc <-  function(n, text, align = "l") paste0("\\multicolumn{", n, "}{", align, "}{", text, "}")
      ##
      general_rows <-  c(
            paste0(mc(3, "\\makecell[l]{$N_{\\text{chains}}$ of the non-WCP arms\\\\($N_{\\text{threads/chain}} = 1$)}"), " & ",
                   mc(2, chain_cells[["HPC"]]), " & ", mc(2, chain_cells[["Laptop"]]), " \\\\"),
            paste0(mc(3, paste0("\\makecell[l]{WCP allocations, $N_{\\text{chains}} \\times N_{\\text{threads/chain}}$\\\\",
                                "(NicoStan+BayesMVP and Stan model)}")), " & ",
                   mc(2, allocation_cells[["HPC"]]), " & ", mc(2, allocation_cells[["Laptop"]]), " \\\\"),
            paste0(mc(3, "WCP allocations (\\texttt{Mplus\\_WCP})"), " & ",
                   mc(2, mplus_cells[["HPC"]]), " & ", mc(2, mplus_cells[["Laptop"]]), " \\\\"),
            "\\hline",
            paste0(mc(3, "$L$ (NicoStan+BayesMVP / Stan model)"), " & ",
                   mc(4, paste0("$", L_bayesmvp, "$ / $", L_stan, "$ (fixed; no randomisation)")), " \\\\"),
            paste0(mc(3, "Step size (NicoStan+BayesMVP / Stan model)"), " & ",
                   mc(4, paste0(fn_format_power_of_ten(eps_bayesmvp), " / ", fn_format_power_of_ten(eps_stan))), " \\\\"),
            paste0(mc(3, "Metric (NicoStan+BayesMVP / Stan model)"), " & ",
                   mc(4, paste0(metric_bayesmvp, " / ", metric_stan)), " \\\\"),
            paste0(mc(3, "\\makecell[l]{Repeats (NicoStan+BayesMVP\\\\and Stan model / Mplus)}"), " & ",
                   mc(4, paste0("$", repeats_BS, "$ / $", repeats_Mplus, "$")), " \\\\"),
            paste0(mc(3, "Mplus timed-run mode"), " & ",
                   mc(4, paste0("\\makecell[l]{", paste0(mplus_modes, collapse = " and "), "\\\\(overhead run: ",
                                       mplus_mode_short, ")}")),
                   " \\\\"))
      ##
      caption_text <-  c(
            "\\scriptfootnotesize{",
            "Settings of the sampling-phase study (E1 Parts I-V and E2), per dataset size $N$.",
            "$N_{\\text{iter}}$ is the number of iterations per chain of the timed run ($n_l$ in equation",
            "\\ref{eq:paper1_two_run_timing});",
            paste0("each timed run was paired with an overhead run of $", overhead_iterations, "$ iteration,"),
            "and each pair was preceded by an untimed warm-up call of one iteration.",
            "Standard refers to the non-WCP arms",
            "(\\texttt{MD\\_BayesMVP} and \\texttt{MD\\_BayesMVP\\_chunking};",
            "\\texttt{AD\\_Stan}, \\texttt{AD\\_Stan\\_chunked} and \\texttt{AD\\_Stan\\_tape\\_chunked};",
            "\\texttt{Mplus\\_standard}),",
            "and WCP to both WCP-only (\\texttt{MD\\_BayesMVP\\_WCP}, \\texttt{AD\\_Stan\\_WCP})",
            "and chunking + WCP (\\texttt{MD\\_BayesMVP\\_WCP\\_chunking}, \\texttt{AD\\_Stan\\_WCP\\_chunking}),",
            paste0("with separate $N_{\\text{iter}}$ for $N_{\\text{chains}} < ", many_chains_minimum,
                   "$ and $N_{\\text{chains}} \\ge ", many_chains_minimum, "$ (local-HPC only);"),
            paste0("for $N_{\\text{chains}} \\ge ", many_chains_minimum, "$, the timed run used $n_l$ iterations,"),
            paste0("and the times were scaled to the $N_{\\text{iter}}$ for $N_{\\text{chains}} < ",
                   many_chains_minimum, "$"),
            "(see section \\ref{section:paper1_chunk_wcp_selection_design});",
            "values in brackets are the WCP-only $n_l$, where these differ.",
            "Candidate $N_{\\text{chunks}}$ apply to every chunked arm",
            "(for chunking + WCP, restricted to $N_{\\text{chunks}} \\ge N_{\\text{threads/chain}}$;",
            "WCP-only uses $N_{\\text{chunks}} = N_{\\text{threads/chain}}$).",
            "Extra $N_{\\text{chunks}}$ (local-HPC / laptop) were run with $N_{\\text{chains}} = 1$ only,",
            "as references for \\texttt{MD\\_BayesMVP} and \\texttt{AD\\_Stan\\_tape\\_chunked}",
            "(see section \\ref{section:paper1_timing_scaling_definitions}).",
            "WCP allocations are listed as $N_{\\text{chains}} \\times \\{N_{\\text{threads/chain}}\\}$.",
            "}")
      ##
      table_lines <-  c(
            "\\begin{table}[H]",
            "\\centering",
            "\\caption{",
            caption_text,
            "}",
            "\\label{table:ps2_parallel_scalability_algorithm_setup_parameters}",
            "\\scriptfootnotesize",
            "\\setlength{\\tabcolsep}{4pt}",
            "\\PaperOneFitTable{",
            "\\begin{tabular}{lllrrrr}",
            "\\hline",
            paste0("$N$ & Implementation & Candidate $N_{\\text{chunks}}$ & ",
                   "\\makecell[l]{Extra $N_{\\text{chunks}}$\\\\(local-HPC / laptop)} & ",
                   "\\makecell[r]{$N_{\\text{iter}}$\\\\standard} & ",
                   "\\makecell[r]{$N_{\\text{iter}}$ WCP\\\\$N_{\\text{chains}} < ", many_chains_minimum, "$} & ",
                   "\\makecell[r]{$n_l$ WCP\\\\$N_{\\text{chains}} \\ge ", many_chains_minimum, "$} \\\\"),
            "\\hline",
            panel_A_lines,
            "\\hline",
            paste0(mc(3, ""), " & ", mc(2, "local-HPC"), " & ", mc(2, "laptop"), " \\\\"),
            "\\hline",
            general_rows,
            "\\hline",
            "\\end{tabular}",
            "}",
            "\\end{table}")
      ##
      if (!dir.exists(paths = output_dir)) dir.create(path = output_dir, recursive = TRUE)
      writeLines(text = table_lines, con = output_file)
      message(paste0("\033[32mWrote ", output_file, "\033[0m"))
      ##
      message("\033[36mCheck list (every value read from the runner's settings and benchmark grid):\033[0m")
      for (check_line in check_list) message(paste0("  ", check_line))
}























