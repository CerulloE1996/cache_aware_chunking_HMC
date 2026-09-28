##
## ---- Stan working-set measurement: data and parameter JSON files (no sampling) ----------------------------------------
##
ws_dir <-  getwd()   ## run from working_set/stan
p1_dir <-  normalizePath(file.path(getwd(), "..", ".."))   ## paper_1_chunking_and_parallel_scalability
##
runner_fns <-  new.env(parent = globalenv())
sys.source(file = file.path(p1_dir, "R_fns_alg_paper_1_chunking_WCP_par_scaling.R"), envir = runner_fns)
##
datasets <-  readRDS(file = file.path(p1_dir, "paper_1_computational_outputs", "HPC", "COVID19_datasets.rds"))
y_full <-  datasets$y_list[[which(datasets$N_vec == 50000)]]
Sigma_nd <-  datasets$Sigma_nd_true_observed_list[[which(datasets$N_vec == 50000)]]
Sigma_d <-  datasets$Sigma_d_true_observed_list[[which(datasets$N_vec == 50000)]]
Se <-  pmin(datasets$Se_true_observed_list[[which(datasets$N_vec == 50000)]]$binary, 0.999)
Sp <-  pmin(datasets$Sp_true_observed_list[[which(datasets$N_vec == 50000)]]$binary, 0.999)
prev <-  datasets$prev_true_observed_list[[which(datasets$N_vec == 50000)]]$overall
##
cases <-  data.frame( N = c(500, 1000, 2500, 5000, 10000, 20000, 20000, 20000),
                      chunk_size = c(500, 1000, 2500, 5000, 10000, 20000, 5000, 1000))
##
for (case_index in seq_len(nrow(cases))) {

    N <-  cases$N[case_index]
    stan_data <-  runner_fns$fn_paper1_stan_data(y = y_full[1:N, , drop = FALSE])
    stan_data$chunk_size <-  as.integer(cases$chunk_size[case_index])
    cmdstanr::write_stan_json(data = stan_data,
                              file = file.path(ws_dir, "json", paste0("data_N", N, "_cs", cases$chunk_size[case_index], ".json")),
                              always_decimal = FALSE)
    ##
    ## ---- point A: u ~ U(0, 1), true Se/Sp/prevalence, first-column correlations from the simulation truth
    set.seed(1)
    u <-  matrix(stats::runif(N * 6), nrow = N)
    inits_A <-  list( u_raw = atanh(2 * u - 1),
                      col_one_raw = rbind(atanh(Sigma_nd[2:6, 1]), atanh(Sigma_d[2:6, 1])),
                      off_raw = matrix(0, nrow = 2, ncol = 10),
                      beta_vec = c(stats::qnorm(1 - Sp), stats::qnorm(Se)),
                      p_raw = array(atanh(2 * prev - 1), dim = 1))
    cmdstanr::write_stan_json(data = inits_A, file = file.path(ws_dir, "json", paste0("init_A_N", N, ".json")))
    ##
    ## ---- point B: all unconstrained parameters at zero
    inits_B <-  list( u_raw = matrix(0, nrow = N, ncol = 6),
                      col_one_raw = matrix(0, nrow = 2, ncol = 5),
                      off_raw = matrix(0, nrow = 2, ncol = 10),
                      beta_vec = rep(0, 12),
                      p_raw = array(0, dim = 1))
    cmdstanr::write_stan_json(data = inits_B, file = file.path(ws_dir, "json", paste0("init_B_N", N, ".json")))

}
message(paste0("Wrote ", nrow(cases), " data files and ", 2 * nrow(cases), " init files"))
