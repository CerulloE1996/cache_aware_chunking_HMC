## -| ------------------ R function to simulate a single dataset for (mixed binary+ordinal) LC_MVOP  ----------------------------
##
## This is the single-study analogue of R_fn_simulate_binary_LC_MVP_data,
## extended to handle mixed binary + ordinal outcomes with multi-population prevalence.
##
## Key differences from binary LC-MVP version:
##   - First n_binary_tests columns are binary (0/1)
##   - Remaining n_ordinal_tests columns are ordinal (0, 1, ..., K_t - 1)
##   - Cutpoints (thresholds) define ordinal category boundaries
##   - Cutpoints can be class-specific (diseased vs non-diseased)
##   - Multi-population prevalence: each observation belongs to a population with its own prev
##
## Key differences from NMA ordinal version:
##   - Single study, NOT meta-analysis (no between-study heterogeneity / random effects)
##   - Full correlation structure between ALL tests (MVP framework)
##   - Uses GHK-style latent variable generation (multivariate normal → categorise)
##



## -| --------- Helper: convert simulation output to BayesMVP model_args_list format  ------------------------------------------

#' convert_sim_to_model_args_list
#' @export
convert_sim_to_model_args_list <-  function( sim_output,
                                            dataset_index = 1L
) {
  
        y   <-  sim_output$y_list[[dataset_index]]
        X   <-  sim_output$X_list[[dataset_index]]
        pop <-  sim_output$pop_list[[dataset_index]]
        N   <-  nrow(y)
        
        dgm <-  sim_output$DGM_info
        
        model_args_list <-  list(
          y = y,
          X = X,
          N = N,
          n_pops = dgm$n_pops,
          pop = pop,
          prior_prev_a = rep(1.0, dgm$n_pops),
          prior_prev_b = rep(1.0, dgm$n_pops)
        )
        
        return(model_args_list)
  
}




## -| --------- Helper: convert simulation output to BayesMVP model_args_list format  ------------------------------------------

#' convert_sim_to_model_args_list_COVID19
#' @export
convert_sim_to_model_args_list_COVID19 <-  function( sim_output,
                                                    dataset_index = 1L
) {
  
  y   <-  sim_output$y_list[[dataset_index]]
  X   <-  sim_output$X_list[[dataset_index]]
  pop <-  sim_output$pop_list[[dataset_index]]
  N   <-  nrow(y)
  
  dgm <-  sim_output$DGM_info
  
  model_args_list <-  list(
    y = y,
    X = X,
    N = N,
    n_pops = dgm$n_pops,
    pop = pop,
    prior_prev_a = rep(1.0, dgm$n_pops),
    prior_prev_b = rep(1.0, dgm$n_pops)
  )
  
  return(model_args_list)
  
}



## -| ------------------ R function to simulate a single dataset for (pure binary) LC_MVP  ----------------------------
##
## Based on true values from an LC-MVP model fitted to the COVID-19 antibody dataset.
## 6 binary tests: euroimmun, roche, abc19, surescreen, orientgene, biomerica
## No reference test (no gold standard).
##
## 2 populations: Police & Fire (reference) and HCW (dummy covariate).
## X matrix: [intercept, dummy_HCW] per observation.
##

#' R_fn_simulate_binary_LC_MVP_data_COVID_19
#' @export
R_fn_simulate_binary_LC_MVP_data_COVID_19 <-  function( study_type = "algorithm_benchmarking",
                                                      seed,
                                                      true_prev_per_pop = c(0.0785, 0.175),  ## (Police & Fire, HCW)
                                                      pop_proportions = NULL,
                                                      corr_force_positive_DGM = FALSE
) {
  
          if (study_type == "sim_study") {
            
                N_vec <-  c( 300, 
                            3000)
                
          } else if (study_type == "algorithm_benchmarking") {
            
                N_vec <-  c( 500,
                            2500,
                            10000)
            
          } else if (study_type == "algorithm_parallel_scaling_tests") {
            
                N_vec <-  c( 500,
                            2500,
                            # 5000,
                            10000, 
                            ## 25000, 
                            50000)
            
          }
          ##
          ## ---- Test configuration:
          ##
          n_tests <-  6L
          test_names <-  c("euroimmun", "roche", "abc19", "surescreen", "orientgene", "biomerica")
          ##
          ## ---- Population setup (2 populations):
          ##
          n_pops <-  2L
          pop_names <-  c("Police_and_Fire", "HCW")
          if (is.null(pop_proportions)) {
            pop_proportions <-  rep(1.0 / n_pops, n_pops)
          }
          stopifnot(length(pop_proportions) == n_pops)
          stopifnot(abs(sum(pop_proportions) - 1.0) < 1e-10)
          stopifnot(length(true_prev_per_pop) == n_pops)
          ##
          ## ---- Covariate setup (intercept-only):
          ##
          n_covariates <-  1L
          ##
          X_per_pop <-  matrix(1, nrow = n_pops, ncol = 1)
          ##
          ## ---- Set the seed:
          ##
          set.seed(seed, kind = "L'Ecuyer-CMRG")
          ##
          ## ---- Make storage lists:
          ##
          y_list <-  list()
          X_list <-  list()
          pop_list <-  list()
          Sigma_nd_true_observed_list <-  Sigma_d_true_observed_list <-  list()
          prev_true_observed_list <-  Se_true_observed_list <-  Sp_true_observed_list <-  list()
          ##
          ## ====================================================================================
          ## Define coefficients: beta[class][k, test]
          ##   class 1 = non-diseased, class 2 = diseased
          ##   k = 1: intercept, k = 2: dummy_HCW
          ##
          ## In the probit model:
          ##   P(Y=1 | D-) = Phi(beta_nd[1, t])  →  Sp = 1 - Phi(beta_nd[1, t]) = Phi(-beta_nd[1, t])
          ##   P(Y=1 | D+) = Phi(beta_d[1, t])   →  Se = Phi(beta_d[1, t])
          ##
          ## So for high Sp: beta_nd should be negative (large negative → high Sp)
          ##    for high Se: beta_d should be positive (large positive → high Se)
          ## ====================================================================================
          ##
          beta_array <-  list()
          beta_array[[1]] <-  matrix(NA, nrow = n_covariates, ncol = n_tests)
          beta_array[[2]] <-  matrix(NA, nrow = n_covariates, ncol = n_tests)
          ##
          ## ---- Class 1 (non-diseased) coefficients:
          ##                              euroimmun  roche    abc19  surescreen  orientgene  biomerica
          beta_array[[1]][1, ] <-  c(        -2.44,     -2.47,     -2.32,     -1.96,        -1.94,        -1.41)   ## intercepts
          ##
          ## ---- Class 2 (diseased) coefficients:
          ##                              euroimmun  roche    abc19  surescreen  orientgene  biomerica
          beta_array[[2]][1, ] <-  c(        1.82,     2.40,     1.17,     1.51,        1.78,        1.55)   ## intercepts
          ##
          ## ====================================================================================
          ## Define correlation structures:
          ## ====================================================================================
          ##
          ## ---- Non-diseased class correlation matrix (6x6):
          ##
          Omega_nd <-  structure(c(1, -0.139, 0.361, 0.341, 0.449, -0.266, -0.139, 1, 
                                  -0.059, -0.123, -0.145, 0.241, 0.361, -0.059, 1, 0.447, 0.495, 
                                  0.233, 0.341, -0.123, 0.447, 1, 0.792, 0.025, 0.449, -0.145, 
                                  0.495, 0.792, 1, -0.025, -0.266, 0.241, 0.233, 0.025, -0.025, 
                                  1), dim = c(6L, 6L))
          diag(Omega_nd) <-  1.0
          ## TODO: Fill in lower triangle from fitted model, e.g.:
          ## Omega_nd[2, 1] <-  xxx  ## roche-euroimmun
          ## Omega_nd[3, 1] <-  xxx  ## abc19-euroimmun
          ## Omega_nd[3, 2] <-  xxx  ## abc19-roche
          ## Omega_nd[4, 1] <-  xxx  ## surescreen-euroimmun
          ## Omega_nd[4, 2] <-  xxx  ## surescreen-roche
          ## Omega_nd[4, 3] <-  xxx  ## surescreen-abc19
          ## Omega_nd[5, 1] <-  xxx  ## orientgene-euroimmun
          ## Omega_nd[5, 2] <-  xxx  ## orientgene-roche
          ## Omega_nd[5, 3] <-  xxx  ## orientgene-abc19
          ## Omega_nd[5, 4] <-  xxx  ## orientgene-surescreen
          ## Omega_nd[6, 1] <-  xxx  ## biomerica-euroimmun
          ## Omega_nd[6, 2] <-  xxx  ## biomerica-roche
          ## Omega_nd[6, 3] <-  xxx  ## biomerica-abc19
          ## Omega_nd[6, 4] <-  xxx  ## biomerica-surescreen
          ## Omega_nd[6, 5] <-  xxx  ## biomerica-orientgene
          ## Omega_nd <-  Omega_nd + t(Omega_nd) - diag(diag(Omega_nd))  ## symmetrise
          ##
          ## ---- Diseased class correlation matrix (6x6):
          ##
          Omega_d <-  structure(c(1, 0.046, 0.298, 0.414, 0.437, -0.068, 0.046, 1, 
                                 0.166, 0.161, 0.121, 0.468, 0.298, 0.166, 1, 0.686, 0.665, 0.374, 
                                 0.414, 0.161, 0.686, 1, 0.7, 0.313, 0.437, 0.121, 0.665, 0.7, 
                                 1, 0.244, -0.068, 0.468, 0.374, 0.313, 0.244, 1), dim = c(6L, 
                                                                                           6L))
          diag(Omega_d) <-  1.0
          ## TODO: Fill in lower triangle from fitted model (same pattern as above)
          ## Omega_d <-  Omega_d + t(Omega_d) - diag(diag(Omega_d))  ## symmetrise
          ##
          require(matrixcalc)
          ## stopifnot(is.positive.definite(Omega_nd))
          ## stopifnot(is.positive.definite(Omega_d))
          ##
          ii_dataset <-  0
          ##
          for (N in N_vec) {
            
                ii_dataset <-  ii_dataset + 1
                ##
                ## ====================================================================
                ## Cholesky decompositions
                ## ====================================================================
                L_Omega_d  <-  t(chol(Omega_d))
                L_Omega_nd <-  t(chol(Omega_nd))
                ##
                ## ====================================================================
                ## Assign observations to populations
                ## ====================================================================
                pop <-  sample(1:n_pops, size = N, replace = TRUE, prob = pop_proportions)
                ##
                ## ====================================================================
                ## Build X matrix (N x n_covariates)
                ## ====================================================================
                X_mat <-  X_per_pop[pop, , drop = FALSE]
                ##
                ## ====================================================================
                ## Generate disease indicators (population-specific prevalence)
                ## ====================================================================
                d_ind <-  integer(N)
                for (g in 1:n_pops) {
                  idx_g <-  which(pop == g)
                  d_ind[idx_g] <-  rbinom(n = length(idx_g), size = 1, prob = true_prev_per_pop[g])
                }
                ##
                n_pos <-  sum(d_ind)
                n_neg <-  N - n_pos
                ##
                ## ====================================================================
                ## Compute observation-specific means: mu_i[t] = X_i %*% beta[class, , t]
                ## ====================================================================
                mu_mat <-  matrix(NA, nrow = N, ncol = n_tests)
                ##
                for (i in 1:N) {
                  cls <-  ifelse(d_ind[i] == 1, 2, 1)
                  for (t in 1:n_tests) {
                    mu_mat[i, t] <-  sum(X_mat[i, ] * beta_array[[cls]][, t])
                  }
                }
                ##
                ## ====================================================================
                ## Generate latent continuous values (multivariate normal)
                ## ====================================================================
                pos_idx <-  which(d_ind == 1)
                neg_idx <-  which(d_ind == 0)
                ##
                latent_results <-  matrix(NA, nrow = N, ncol = n_tests)
                ##
                ## Non-diseased:
                if (n_neg > 0) {
                  Z_nd <-  LaplacesDemon::rmvn(n = n_neg, mu = rep(0, n_tests), Sigma = Omega_nd)
                  latent_results[neg_idx, ] <-  Z_nd + mu_mat[neg_idx, , drop = FALSE]
                }
                ## Diseased:
                if (n_pos > 0) {
                  Z_d <-  LaplacesDemon::rmvn(n = n_pos, mu = rep(0, n_tests), Sigma = Omega_d)
                  latent_results[pos_idx, ] <-  Z_d + mu_mat[pos_idx, , drop = FALSE]
                }
                ##
                ## ====================================================================
                ## Dichotomise to get observed binary outcomes y (threshold at 0)
                ## ====================================================================
                y <-  matrix(NA_real_, nrow = N, ncol = n_tests)
                ##
                for (t in 1:n_tests) {
                  y[, t] <-  ifelse(latent_results[, t] > 0, 1, 0)
                }
                ##
                ## ====================================================================
                ## Compute observed summary quantities
                ## ====================================================================
                prev_observed_per_pop <-  numeric(n_pops)
                for (g in 1:n_pops) {
                  idx_g <-  which(pop == g)
                  prev_observed_per_pop[g] <-  round(mean(d_ind[idx_g]), 4)
                }
                prev_observed_overall <-  round(mean(d_ind), 4)
                ##
                Se_observed <-  numeric(n_tests)
                Sp_observed <-  numeric(n_tests)
                for (t in 1:n_tests) {
                  if (n_pos > 0) Se_observed[t] <-  round(mean(y[pos_idx, t]), 4)
                  if (n_neg > 0) Sp_observed[t] <-  round(1 - mean(y[neg_idx, t]), 4)
                }
                ##
                Sigma_nd_true_observed <-  cor(latent_results[neg_idx, , drop = FALSE])
                Sigma_d_true_observed  <-  cor(latent_results[pos_idx, , drop = FALSE])
                ##
                ## ---- Print summary:
                ##
                cat("\n====================", "| N =", N, "====================\n")
                cat("n_pops =", n_pops, "\n")
                cat("n_covariates =", n_covariates, "\n")
                cat("Observed prev per pop =", prev_observed_per_pop, "\n")
                cat("Observed prev overall =", prev_observed_overall, "\n")
                for (t in 1:n_tests) {
                  cat(paste0("  ", test_names[t], ": Se = ", sub("^ +", "", sub("^ +", "", formatC(Se_observed[t], format = "f", digits = 3))),
                             ", Sp = ", sub("^ +", "", sub("^ +", "", formatC(Sp_observed[t], format = "f", digits = 3))), "\n"))
                }
                ##
                ## ---- Build per-test X matrices (format expected by BayesMVP):
                ## X[[class]][[test]] is an (N x n_covariates) matrix
                X_for_BayesMVP <-  list()
                for (cls in 1:2) {
                  X_for_BayesMVP[[cls]] <-  list()
                  for (t in 1:n_tests) {
                    X_for_BayesMVP[[cls]][[t]] <-  X_mat
                  }
                }
                ##
                ## ---- Store:
                ##
                y_list[[ii_dataset]] <-  y
                X_list[[ii_dataset]] <-  X_for_BayesMVP
                pop_list[[ii_dataset]] <-  pop
                Sigma_nd_true_observed_list[[ii_dataset]] <-  Sigma_nd_true_observed
                Sigma_d_true_observed_list[[ii_dataset]] <-  Sigma_d_true_observed
                prev_true_observed_list[[ii_dataset]] <-  list(per_pop = prev_observed_per_pop,
                                                              overall = prev_observed_overall)
                Se_true_observed_list[[ii_dataset]] <-  list(binary = Se_observed)
                Sp_true_observed_list[[ii_dataset]] <-  list(binary = Sp_observed)
            
          } ## end N loop
          
          ## ====================================================================
          ## Return
          ## ====================================================================
          return(list(
            y_list = y_list,
            X_list = X_list,
            pop_list = pop_list,
            N_vec = N_vec,
            n_tests = n_tests,
            n_binary_tests = n_tests,
            n_ordinal_tests = 0L,
            n_cat_per_ord_test = integer(0),
            n_thr_per_ord_test = integer(0),
            n_covariates = n_covariates,
            test_names = test_names,
            Sigma_nd_true_observed_list = Sigma_nd_true_observed_list,
            Sigma_d_true_observed_list  = Sigma_d_true_observed_list,
            prev_true_observed_list = prev_true_observed_list,
            Se_true_observed_list = Se_true_observed_list,
            Sp_true_observed_list = Sp_true_observed_list,
            DGM_info = list(
              n_tests = n_tests,
              n_binary_tests = n_tests,
              n_ordinal_tests = 0L,
              n_cat_per_ord_test = integer(0),
              n_thr_per_ord_test = integer(0),
              n_covariates = n_covariates,
              beta_array = beta_array,
              Omega_d  = Omega_d,
              Omega_nd = Omega_nd,
              true_prev_per_pop = true_prev_per_pop,
              pop_proportions = pop_proportions,
              n_pops = n_pops,
              pop_names = pop_names,
              test_names = test_names,
              seed = seed
            )
          ))
  
}







# 
# 
# 
#  
# R_fn_simulate_mixed_binary_ordinal_LC_MVOP_data <-  function( study_type = "algorithm_benchmarking",
#                                                              seed,
#                                                              DGM,
#                                                              true_prev_per_pop = c(0.273, 0.117, 0.0324),
#                                                              pop_proportions = NULL,  ## relative size of each population (must sum to 1)
#                                                              grouping
# ) {
#   
#         if (study_type == "sim_study") {
#           
#               N_vec <-  c(300, 
#                          3000)
#               
#         } else if (study_type == "algorithm_benchmarking") {
#           
#               N_vec <-  c(500,
#                          2500,
#                          10000,
#                          25000,
#                          50000)
#               
#         } else if (study_type == "algorithm_parallel_scaling_tests") {
#           
#               N_vec <-  c(500,
#                          2500,
#                          10000,
#                          25000,
#                          50000)
#               
#         }
#         ##
#         ## ---- Test configuration:
#         ##
#         ## 1 binary test + 2 ordinal tests = 3 tests total
#         n_binary_tests  <-  1L
#         n_ordinal_tests <-  2L
#         n_tests <-  n_binary_tests + n_ordinal_tests
#         ##
#         ## ---- Ordinal test configuration:
#         ##
#         ## Test 2 (ordinal): e.g. 7 categories → 6 thresholds  (like GAD-7 w/ fewer cats, or PHQ-2 etc.)
#         ## Test 3 (ordinal): e.g. 21 categories → 20 thresholds (like HADS-A or similar)
#         ## Adjust these to match the real data:
#         n_cat_per_ord_test <-  c(28L, 31L)   ## number of categories for each ordinal test
#         n_thr_per_ord_test <-  n_cat_per_ord_test - 1L  ## number of thresholds = categories - 1
#         ##
#         ## ---- Multi-population setup:
#         n_pops <-  length(true_prev_per_pop)
#         if (is.null(pop_proportions)) {
#           ## Default: equal population sizes
#           pop_proportions <-  rep(1.0 / n_pops, n_pops)
#         }
#         stopifnot(length(pop_proportions) == n_pops)
#         stopifnot(abs(sum(pop_proportions) - 1.0) < 1e-10)
#         ##
#         ## ---- Set the seed (keep OUTSIDE the N loop):
#         ##
#         set.seed(seed, kind = "L'Ecuyer-CMRG")
#         ##
#         ## ---- Make storage lists:
#         ##
#         y_list <-  list()
#         pop_list <-  list()
#         Sigma_nd_true_observed_list <-  Sigma_d_true_observed_list <-  list()
#         prev_true_observed_list <-  Se_true_observed_list <-  Sp_true_observed_list <-  list()
#         true_estimates_observed_list <-  list()
#         ##
#         ## ====================================================================================
#         ## Define Se/Sp for the BINARY test(s)
#         ## ====================================================================================
#         ## Only 1 binary test here:
#         Se_binary <-  0.846
#         Sp_binary <-  0.921
#         Fp_binary <-  1 - Sp_binary
#         ##
#         ## ====================================================================================
#         ## Define mean accuracy (location) for ordinal tests
#         ## ====================================================================================
#         ## beta[class, intercept, test]  - intercept-only model so just 1 covariate (intercept)
#         ## Class 1 = non-diseased, Class 2 = diseased
#         ## These are Phi^{-1}(Fp) and Phi^{-1}(Se) equivalents on probit scale
#         ##
#         ## For ordinal tests, beta defines the LOCATION shift of the latent normal,
#         ## analogous to qnorm(Fp) and qnorm(Se) in the binary case.
#         ##
#         beta_ord <-  array(NA, dim = c(2, n_ordinal_tests))  ## [class, ordinal_test]
#         ##
#         ## Ordinal test 1 (test index 2 overall): moderate accuracy
#         beta_ord[1, 1] <-  -1.556 ## non-diseased location (class 1)
#         beta_ord[2, 1] <-  +0.167 ## diseased location (class 2)
#         ##
#         ## Ordinal test 2 (test index 3 overall): slightly different accuracy  
#         beta_ord[1, 2] <-  -1.253  ## non-diseased location
#         beta_ord[2, 2] <-  +1.253  ## diseased location
#         ##
#         ## ====================================================================================
#         ## Define cutpoints for ordinal tests:
#         ## ====================================================================================
#         ## Cutpoints are class-specific: C_vec[[class]][[ordinal_test]] = vector of thresholds
#         ## These define category boundaries on the probit scale.
#         ##
#         ## Using equally-spaced cutpoints as a simple DGM (these can be replaced with 
#         ## induced-Dirichlet-generated ones for more realistic spacing).
#         ##
#         C_vec <-  list()  ## C_vec[[class]][[ord_test_index]]
#         
#         for (c in 1:2) {
#             
#             C_vec[[c]] <-  list()
#             
#             for (t_ord in 1:n_ordinal_tests) {
#               
#                 n_thr_t <-  n_thr_per_ord_test[t_ord]
#                 ## Generate roughly equally-spaced cutpoints spanning a reasonable range
#                 ## Adjust range based on number of thresholds
#                 range_lower <-  -3.0
#                 range_upper <-  3.0
#                 C_vec[[c]][[t_ord]] <-  seq(from = range_lower, 
#                                            to = range_upper, 
#                                            length.out = n_thr_t)
#               
#             }
#           
#         }
#         ##
#         ## Class = 1, ord_test = 1:
#         ##
#         C_vec[[1]][[1]] <-  c(  -2.721, -2.295, -1.946, -1.67, -1.407, -1.18, -0.974, -0.799, 
#                                -0.644, -0.465, -0.285, -0.138, -0.024, 0.043, 0.268, 0.327, 
#                                0.489, 0.636, 0.876, 0.989, 1.112, 1.233, 1.402, 1.541, 1.709, 
#                                1.918, 2.272)
#         ##
#         ## Class = 2, ord_test = 1:
#         ##
#         C_vec[[2]][[1]] <-  c(-2.049, -1.711, -1.486, -1.325, -1.184, -0.916, -0.792, -0.684, 
#                              -0.552, -0.485, -0.397, -0.3, -0.237, -0.144, -0.028, 0.131, 
#                              0.301, 0.393, 0.628, 0.753, 0.853, 1.11, 1.188, 1.426, 1.651, 
#                              1.809, 2.033)
#         ##
#         ## Class = 1, ord_test = 2:
#         ##
#         C_vec[[1]][[1]] <-  c(-2.402, -2.136, -1.953, -1.682, -1.524, -1.315, -1.16, -1.006, 
#                              -0.848, -0.68, -0.552, -0.42, -0.251, -0.115, 0.016, 0.108, 0.187, 
#                              0.318, 0.453, 0.556, 0.791, 0.977, 1.106, 1.177, 1.373, 1.485, 
#                              1.601, 1.801, 1.991, 2.321)
#         ##
#         ## Class = 2, ord_test = 2:
#         ##
#         C_vec[[2]][[1]] <-  c(-2.072, -1.682, -1.446, -1.255, -1.114, -0.985, -0.866, -0.73, 
#                              -0.609, -0.527, -0.446, -0.37, -0.3, -0.213, -0.119, -0.022, 
#                              0.065, 0.152, 0.254, 0.355, 0.486, 0.617, 0.706, 0.758, 0.927, 
#                              1.179, 1.366, 1.672, 1.9, 2.091)
#         ##
#         ## ====================================================================================
#         ## Define correlation structures:
#         ## ====================================================================================
#         Omega_CI <-  diag(n_tests)
#         ##
#         ## ---- "Varied" structure (rank-1 + diagonal, like latent trait model):
#         ##
#         LT_bs <-  c(0.40, 0.75, 1.10)  ## length = n_tests
#         Sigma_varied <-  diag(n_tests) + t(t(LT_bs)) %*% (t(LT_bs))
#         Omega_varied <-  cov2cor(Sigma_varied)
#         Omega_varied <-  round(Omega_varied, 3)
#         ##
#         LT_bs_half <-  0.5 * LT_bs
#         Sigma_varied_half <-  diag(n_tests) + t(t(LT_bs_half)) %*% (t(LT_bs_half))
#         Omega_varied_half <-  cov2cor(Sigma_varied_half)
#         Omega_varied_half <-  round(Omega_varied_half, 3)
#         ##
#         ## ---- "Highly varied" (heterogeneous) structure - only from LC-MVP model:
#         ##
#         Sigma_highly_varied <-  matrix(c(1,     0.25,  0.25,
#                                         0.25,  1,     0.55,
#                                         0.25,  0.55,  1),
#                                       n_tests, n_tests)
#         Omega_highly_varied <-  cov2cor(Sigma_highly_varied)
#         ##
#         Sigma_highly_varied_half <-  0.5 * Sigma_highly_varied
#         diag(Sigma_highly_varied_half) <-  rep(1, n_tests)
#         Omega_highly_varied_half <-  cov2cor(Sigma_highly_varied_half)
#         ##
#         Sigma_highly_varied_quarter <-  0.25 * Sigma_highly_varied
#         diag(Sigma_highly_varied_quarter) <-  rep(1, n_tests)
#         ##
#         ## ---- Using REAL model-estimates from model fitted to REAL (Baron et al.) depression data:
#         ##
#         # Sigma_Baron_real <-  matrix(c(1,     0.25,  0.25,
#         #                                 0.25,  1,     0.55,
#         #                                 0.25,  0.55,  1),
#         #                               n_tests, n_tests)
#         # Omega_Baron_real <-  cov2cor(Sigma_Baron_real)
#         ##
#         Sigma_Baron_real_nd <-  structure(c( 1, 0.659, 0.709, 
#                                             0.659, 1, 0.672,
#                                             0.709, 0.672, 1),
#                                           dim = c(n_tests, n_tests))
#         Omega_Baron_real_nd <-  cov2cor(Sigma_Baron_real_nd)
#         ##
#         Sigma_Baron_real_d <-  structure(c(1, 0.19, 0.189, 
#                                           0.19, 1, 0.418, 
#                                           0.189, 0.418, 1),
#                                         dim = c(n_tests, n_tests))
#         Omega_Baron_real_d <-  cov2cor(Sigma_Baron_real_d)
#         ##
#         ## ---- Check positive definiteness:
#         ##
#         require(matrixcalc)
#         stopifnot(is.positive.definite(Omega_varied))
#         stopifnot(is.positive.definite(Omega_varied_half))
#         ##
#         stopifnot(is.positive.definite(Omega_highly_varied))
#         stopifnot(is.positive.definite(Omega_highly_varied_half))
#         ##
#         stopifnot(is.positive.definite(Omega_Baron_real_nd))
#         stopifnot(is.positive.definite(Sigma_Baron_real_d))
#         
#         ii_dataset <-  0  ## counter
#         
#         for (N in N_vec) {
#           
#           ii_dataset <-  ii_dataset + 1
#           
#           ## ====================================================================
#           ## Set up DGM
#           ## ====================================================================
#           if (DGM == 1) {
#             
#                 ## DGM 1: Conditional independence
#                 Omega_d  <-  Omega_CI
#                 Omega_nd <-  Omega_CI
#             
#           } else if (DGM == 2) {
#             
#                 ## DGM 2: CD in both groups; varied structure (latent-trait-like)
#                 Omega_d  <-  Omega_varied
#                 Omega_nd <-  Omega_varied_half
#             
#           } else if (DGM == 3) {
#             
#                 ## DGM 3: CD in both groups; highly varied structure (LC-MVP-like)
#                 Omega_d  <-  Omega_highly_varied_half
#                 Omega_nd <-  Sigma_highly_varied_quarter
#                 diag(Omega_nd) <-  rep(1, n_tests)
#                 
#           } else if (DGM == 4) {
#             
#                 Omega_d  <-  Omega_Baron_real_nd
#                 Omega_nd <-  Omega_Baron_real_d
#                 diag(Omega_nd) <-  rep(1, n_tests)
#             
#           } else {
#             stop("DGM must be 1, 2, or 3")
#           }
#           
#           ## ---- PD check:
#           L_Omega_d  <-  t(chol(Omega_d))
#           L_Omega_nd <-  t(chol(Omega_nd))
#           
#           ## ====================================================================
#           ## Assign observations to populations
#           ## ====================================================================
#           pop <-  sample(1:n_pops, size = N, replace = TRUE, prob = pop_proportions)
#           
#           ## ====================================================================
#           ## Generate disease indicators (population-specific prevalence)
#           ## ====================================================================
#           d_ind <-  integer(N)
#           for (g in 1:n_pops) {
#             idx_g <-  which(pop == g)
#             d_ind[idx_g] <-  rbinom(n = length(idx_g), size = 1, prob = true_prev_per_pop[g])
#           }
#           
#           n_pos <-  sum(d_ind)
#           n_neg <-  N - n_pos
#           
#           ## ====================================================================
#           ## Build mean vectors for each class
#           ## ====================================================================
#           ## mu_nd[t] and mu_d[t] for each test t:
#           ##   - Binary test:  mu_nd = qnorm(Fp), mu_d = qnorm(Se)
#           ##   - Ordinal test: mu_nd = beta_ord[1, t_ord], mu_d = beta_ord[2, t_ord]
#           ##
#           mu_nd <-  numeric(n_tests)
#           mu_d  <-  numeric(n_tests)
#           
#           ## Binary test(s):
#           for (t in 1:n_binary_tests) {
#             mu_nd[t] <-  qnorm(Fp_binary)   ## = qnorm(1 - Sp)
#             mu_d[t]  <-  qnorm(Se_binary)
#           }
#           
#           ## Ordinal tests:
#           for (t_ord in 1:n_ordinal_tests) {
#             t <-  n_binary_tests + t_ord  ## overall test index
#             mu_nd[t] <-  beta_ord[1, t_ord]
#             mu_d[t]  <-  beta_ord[2, t_ord]
#           }
#           
#           ## ====================================================================
#           ## Generate latent continuous values (multivariate normal)
#           ## ====================================================================
#           pos_idx <-  which(d_ind == 1)
#           neg_idx <-  which(d_ind == 0)
#           
#           latent_results <-  matrix(NA, nrow = N, ncol = n_tests)
#           
#           if (n_neg > 0) {
#             latent_results[neg_idx, ] <-  LaplacesDemon::rmvn(n = n_neg, mu = mu_nd, Sigma = Omega_nd)
#           }
#           if (n_pos > 0) {
#             latent_results[pos_idx, ] <-  LaplacesDemon::rmvn(n = n_pos, mu = mu_d,  Sigma = Omega_d)
#           }
#           
#           ## ====================================================================
#           ## Dichotomise/categorise to get observed outcomes y
#           ## ====================================================================
#           y <-  matrix(NA_real_, nrow = N, ncol = n_tests)
#           
#           ## ---- Binary test(s): threshold at 0
#           for (t in 1:n_binary_tests) {
#             y[, t] <-  ifelse(latent_results[, t] > 0, 1, 0)
#           }
#           
#           ## ---- Ordinal tests: use cutpoints to assign categories
#           for (t_ord in 1:n_ordinal_tests) {
#             
#             t <-  n_binary_tests + t_ord  ## overall test index
#             n_thr_t <-  n_thr_per_ord_test[t_ord]
#             n_cat_t <-  n_cat_per_ord_test[t_ord]
#             
#             for (i in 1:N) {
#               
#                   ## Select class-specific cutpoints
#                   if (d_ind[i] == 1) {
#                     C_t <-  C_vec[[2]][[t_ord]]  ## diseased
#                   } else {
#                     C_t <-  C_vec[[1]][[t_ord]]  ## non-diseased
#                   }
#                   
#                   ## Assign category based on latent value
#                   latent_val <-  latent_results[i, t]
#                   
#                   if (latent_val <= C_t[1]) {
#                     y[i, t] <-  1  ## first category (1-indexed)
#                   } else if (latent_val > C_t[n_thr_t]) {
#                     y[i, t] <-  n_cat_t  ## last category
#                   } else {
#                     for (k in 2:n_thr_t) {
#                       if (latent_val > C_t[k - 1] && latent_val <= C_t[k]) {
#                         y[i, t] <-  k
#                         break
#                       }
#                     }
#                   }
#                   
#             }
#             
#           }
#           
#           ## ====================================================================
#           ## Compute observed summary quantities
#           ## ====================================================================
#           
#           ## ---- Observed prevalence per population:
#           prev_observed_per_pop <-  numeric(n_pops)
#           for (g in 1:n_pops) {
#                 idx_g <-  which(pop == g)
#                 prev_observed_per_pop[g] <-  round(mean(d_ind[idx_g]), 4)
#           }
#           prev_observed_overall <-  round(mean(d_ind), 4)
#           
#           ## ---- Observed Se/Sp for binary test(s):
#           Se_observed <-  numeric(n_binary_tests)
#           Sp_observed <-  numeric(n_binary_tests)
#           
#           for (t in 1:n_binary_tests) {
#                 if (n_pos > 0) Se_observed[t] <-  round(mean(y[pos_idx, t]), 4)
#                 if (n_neg > 0) Sp_observed[t] <-  round(1 - mean(y[neg_idx, t]), 4)
#           }
#           
#           ## ---- Observed Se/Sp at each threshold for ordinal tests:
#           Se_ord_observed <-  list()
#           Sp_ord_observed <-  list()
#           
#           for (t_ord in 1:n_ordinal_tests) {
#             
#                 t <-  n_binary_tests + t_ord
#                 n_cat_t <-  n_cat_per_ord_test[t_ord]
#                 Se_at_thr <-  numeric(n_cat_t - 1)
#                 Sp_at_thr <-  numeric(n_cat_t - 1)
#                 
#                 for (k in 1:(n_cat_t - 1)) {
#                   ## "positive" = score >= k
#                   if (n_pos > 0) Se_at_thr[k] <-  mean(y[pos_idx, t] >= k)
#                   if (n_neg > 0) Sp_at_thr[k] <-  1 - mean(y[neg_idx, t] >= k)
#                 }
#                 
#                 Se_ord_observed[[t_ord]] <-  round(Se_at_thr, 4)
#                 Sp_ord_observed[[t_ord]] <-  round(Sp_at_thr, 4)
#                 
#           }
#           
#           ## ---- Observed correlations (on latent scale):
#           Sigma_nd_true_observed <-  cor(latent_results[neg_idx, , drop = FALSE])
#           Sigma_d_true_observed  <-  cor(latent_results[pos_idx, , drop = FALSE])
#           
#           ## ---- Print summary:
#           cat("\n==================== DGM", DGM, "| N =", N, "====================\n")
#           cat("n_pops =", n_pops, "\n")
#           cat("Observed prev per pop =", prev_observed_per_pop, "\n")
#           cat("Observed prev overall =", prev_observed_overall, "\n")
#           cat("Binary Se =", Se_observed, "\n")
#           cat("Binary Sp =", Sp_observed, "\n")
#           
#           for (t_ord in 1:n_ordinal_tests) {
#             cat(sprintf("Ordinal test %d: Se range [%.3f, %.3f], Sp range [%.3f, %.3f]\n",
#                         t_ord,
#                         min(Se_ord_observed[[t_ord]]), max(Se_ord_observed[[t_ord]]),
#                         min(Sp_ord_observed[[t_ord]]), max(Sp_ord_observed[[t_ord]])))
#           }
#           
#           ## ---- Store:
#           y_list[[ii_dataset]] <-  y
#           pop_list[[ii_dataset]] <-  pop
#           Sigma_nd_true_observed_list[[ii_dataset]] <-  Sigma_nd_true_observed
#           Sigma_d_true_observed_list[[ii_dataset]] <-  Sigma_d_true_observed
#           prev_true_observed_list[[ii_dataset]] <-  list(
#             per_pop = prev_observed_per_pop,
#             overall = prev_observed_overall
#           )
#           Se_true_observed_list[[ii_dataset]] <-  list(
#             binary = Se_observed,
#             ordinal = Se_ord_observed
#           )
#           Sp_true_observed_list[[ii_dataset]] <-  list(
#             binary = Sp_observed,
#             ordinal = Sp_ord_observed
#           )
#           
#         } ## end N loop
#         
#         ## ====================================================================
#         ## Return
#         ## ====================================================================
#         return(list(
#           y_list = y_list,
#           pop_list = pop_list,
#           ##
#           N_vec = N_vec,
#           n_tests = n_tests,
#           n_binary_tests = n_binary_tests,
#           n_ordinal_tests = n_ordinal_tests,
#           n_cat_per_ord_test = n_cat_per_ord_test,
#           n_thr_per_ord_test = n_thr_per_ord_test,
#           ##
#           Sigma_nd_true_observed_list = Sigma_nd_true_observed_list,
#           Sigma_d_true_observed_list  = Sigma_d_true_observed_list,
#           ##
#           prev_true_observed_list = prev_true_observed_list,
#           Se_true_observed_list = Se_true_observed_list,
#           Sp_true_observed_list = Sp_true_observed_list,
#           ##
#           DGM_info = list(
#             DGM = DGM,
#             n_tests = n_tests,
#             n_binary_tests = n_binary_tests,
#             n_ordinal_tests = n_ordinal_tests,
#             n_cat_per_ord_test = n_cat_per_ord_test,
#             n_thr_per_ord_test = n_thr_per_ord_test,
#             ##
#             Se_binary = Se_binary,
#             Sp_binary = Sp_binary,
#             beta_ord = beta_ord,
#             C_vec = C_vec,
#             ##
#             Omega_d  = Omega_d,
#             Omega_nd = Omega_nd,
#             ##
#             true_prev_per_pop = true_prev_per_pop,
#             pop_proportions = pop_proportions,
#             n_pops = n_pops,
#             seed = seed
#           )
#         ))
#   
# }
# 
# 
# 
# 
# 
# 






















