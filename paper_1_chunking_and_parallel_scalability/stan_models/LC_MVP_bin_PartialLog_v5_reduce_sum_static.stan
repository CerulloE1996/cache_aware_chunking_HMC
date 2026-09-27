


////
//// LC_MVP_bin_PartialLog_v5_reduce_sum_static.stan
////
//// Paper 1 deterministic reduce_sum_static variant of LC_MVP_bin_PartialLog_v5.
//// Implements ONLY the (fully_vectorised == 1) && (handle_numerical_issues == 1) path.
//// The flags are still declared in data for drop-in compatibility with the existing
//// Stan_data_list, but they are ignored (add "chunk_size" to the data list).
////
//// KEY STRUCTURAL CHANGES vs baseline:
////   (1) u_raw is now  array[N] row_vector[n_tests]  (was matrix[N, n_tests]) so that
////       reduce_sum can SLICE the nuisance parameters directly - this is the efficient
////       pattern (each worker only receives/copies its own slice of u_raw, instead of
////       every worker copying the full N x n_tests block). An R matrix init of dim
////       (N, n_tests) still works unchanged via JSON (maps to array-of-rows).
////   (2) The likelihood lives in a plain (non-_lp) partial-sum function, so the
////       tanh Jacobian for u is added to the RETURNED value instead of via target +=
////       (reduce_sum partial functions cannot be _lp). Total log-density is identical.
////   (3) log_lik is NOT a transformed parameter here (it would require evaluating the
////       whole likelihood a second time outside reduce_sum). target is incremented
////       directly. lp__ at a given (mapped) parameter point matches the baseline.
////
//// USAGE (cmdstanr):
////   mod <- cmdstan_model("LC_MVP_bin_PartialLog_v5_reduce_sum_static.stan",
////                        cpp_options = list(stan_threads = TRUE))
////   fit <- mod$sample(data = ..., threads_per_chain = K, ...)
////
////   - chunk_size is an upper bound on observations per partial sum. Static partitioning
////     is deterministic, but need not create exactly ceil(N / chunk_size) slices.
////   - threads_per_chain controls available WCP threads independently of chunk_size.
////     Use one thread per chain for the tape-chunking arm, and more than one for WCP.
////

functions {

      vector lb_ub_lp (vector y, real lb, real ub) {
            int N = num_elements(y);
            vector[N] tanh_y = tanh(y);
            target +=  - log(2)  +  log( (ub - lb) * (1 - square(tanh_y))) ;
            return lb +  (ub - lb) *  0.5 * (1 + tanh_y) ;
      }

      real lb_ub_lp (real y, real lb, real ub) {
             real tanh_y = tanh(y);
             target +=  - log(2)  +  log( (ub - lb) * (1 - square(tanh_y))) ;
             return lb +  (ub - lb) *  0.5 * (1 + tanh_y) ;
      }

      matrix cholesky_corr_constrain_outer_lp( vector col_one_raw,
                                               vector off_raw,
                                               real lb,
                                               real ub) {

            int K = num_elements(col_one_raw) + 1;
            vector[K - 1] z = lb_ub_lp(col_one_raw, lb, ub);
            matrix[K, K] L = diag_matrix(rep_vector(1, K));
            vector[K] D;
            D[1] = 1;
            L[2:K, 1] = z[1:K - 1];
            D[2] = 1 - L[2, 1]^2;
            int cnt = 1;

            for (i in 3:K) {
               D[i] = 1 - L[i, 1]^2;
               L[i, 2:i - 1] = rep_row_vector(1 - L[i, 1]^2, i - 2);
               real l_ij_old = L[i, 2];
              for (j in 2:i - 1) {
                real b1 = dot_product(L[j, 1:(j - 1)], D[1:j - 1]' .* L[i, 1:(j - 1)]);

                  // how to derive the bounds
                  // we know that the correlation value C is bound by
                  // b1 - Ljj * Lij_old <= C <= b1 + Ljj * Lij_old
                  // Now we want our bounds to be enforced too so
                  // max(lb, b1 - Ljj * Lij_old) <= C <= min(ub, b1 + Ljj * Lij_old)
                  // We have the Lij_new = (C - b1) / Ljj
                  // To get the bounds on Lij_new is
                  // (bound - b1) / Ljj

                  real low = max({-sqrt(l_ij_old * D[j]), lb - b1});
                  real up = min({sqrt(l_ij_old * D[j]), ub - b1});
                  if (is_nan(low) || is_nan(up) || is_inf(low) || is_inf(up) || low >= up)
                      reject("empty or nonfinite correlation interval");

                  real x = lb_ub_lp(off_raw[cnt], low, up);
                  L[i, j] = x / D[j];

                  target += -0.5 * log(D[j]);

                  l_ij_old -= D[j] * square(L[i, j]);
                  if (is_nan(l_ij_old) || is_inf(l_ij_old) || l_ij_old <= 0)
                      reject("nonpositive or nonfinite LDL remainder");

                  cnt += 1;
                }
                D[i] = l_ij_old;
              }

              if (is_nan(sum(D)) || is_inf(sum(D)) || min(D) <= 0)
                  reject("nonpositive or nonfinite LDL pivot");
              return diag_post_multiply(L, sqrt(D));

      }

      // need to add citation to this (slight modification from a HP. calculators forum post)
      real inv_Phi_approx_from_prob(real p) {
            return 5.494448514153059 *  sinh(0.33333333333333331483 * asinh( 0.34176618822627863 * logit(p)  )) ;
      }

      // need to add citation to this (slight modification from a HP. calculators forum post)
      vector inv_Phi_approx_from_prob(vector p) {
            return 5.494448514153059 *  sinh(0.33333333333333331483 * asinh( 0.34176618822627863 * logit(p)  )) ;
      }

      // need to add citation to this  (slight modification from a HP. calculators forum post)
      real inv_Phi_approx_from_logit_prob(real logit_p) {
            return 5.494448514153059 *  sinh(0.33333333333333331483 * asinh( 0.34176618822627863 * logit_p  )) ;
      }

      // need to add citation to this (slight modification from a HP. calculators forum post)
      vector inv_Phi_approx_from_logit_prob(vector logit_p) {
            return 5.494448514153059 *  sinh(0.33333333333333331483 * asinh( 0.34176618822627863 * logit_p  )) ;
      }

      vector rowwise_sum(matrix M) {      // M is (N x T) matrix
            return M * rep_vector(1.0, cols(M));
      }

      vector rowwise_max(matrix M) {      // M is (N x T) matrix
            int N =  rows(M);
            vector[N] rowwise_maxes;
            for (n in 1:N) {
              rowwise_maxes[n] = max(M[n, ]);
            }
            return rowwise_maxes;
      }

      vector log_sum_exp_2d(matrix array_2d_to_lse) {
            int N = rows(array_2d_to_lse);
            matrix[N, 2] rowwise_maxes_2d_array;
            rowwise_maxes_2d_array[, 1] =  rowwise_max(array_2d_to_lse);
            rowwise_maxes_2d_array[, 2] =  rowwise_maxes_2d_array[, 1];
            return  rowwise_maxes_2d_array[, 1] + log(rowwise_sum(exp((array_2d_to_lse  -  rowwise_maxes_2d_array))));
      }

      ////
      //// reduce_sum partial-sum function.
      //// Slices u_raw by observation. Body = the (fully_vectorised == 1) &&
      //// (handle_numerical_issues == 1) likelihood path from the baseline model,
      //// with N -> M (slice size) and all row-indexed objects taken as start:end.
      //// Returns: (tanh Jacobian for the u slice) + (sum of per-obs log-lik).
      ////
      real partial_log_lik(array[] row_vector u_raw_slice,
                           int start,
                           int end,
                           data matrix y,
                           data array[] matrix X_nd,
                           data array[] matrix X_d,
                           data array[,] int n_covs_per_outcome,
                           data array[] int pop,
                           array[,,] real beta,
                           array[] matrix L_Omega,
                           matrix L_Omega_diag_recip,
                           vector prev,
                           int n_tests,
                           int n_class,
                           int n_covariates_max,
                           int Phi_type,
                           data real overflow_threshold,
                           data real underflow_threshold) {

            int M = end - start + 1;
            real out = 0.0;

            //// ---- u transform + Jacobian (plain-function equivalent of
            //// ---- lb_ub_lp(u_raw, 0, 1); Jacobian added to the return value):
            matrix[M, n_tests] tanh_u;
            for (m in 1:M)  tanh_u[m, ] = tanh(u_raw_slice[m]);
            matrix[M, n_tests] u = 0.5 * (1.0 + tanh_u);
            out += - M * n_tests * log2() + sum(log1m(square(tanh_u)));

            //// ---- slice-local copies of data (data -> no autodiff cost):
            matrix[M, n_tests] y_s = y[start:end, ];
            array[M] int pop_s = pop[start:end];

            matrix[M, n_class] log_prev;
            for (m in 1:M) {
              log_prev[m, 1] = log1m(prev[pop_s[m]]);
              log_prev[m, 2] = log(prev[pop_s[m]]);
            }

            //// ---- likelihood (hni == 1, fully vectorised path; N -> M):
            matrix[M, n_tests] Z_std_norm;
            vector[M] Bound_Z;
            matrix[M, n_class] lp;
            vector[M] inc;
            matrix[M, n_tests] y1;

            for (c in 1:n_class) {

                  inc = rep_vector(0.0, M);

                  for (t in 1:n_tests) {

                             if (n_covariates_max > 1) {
                                  vector[M] Xbeta;
                                  if (c == 1)   Xbeta =  X_nd[t, 1:n_covs_per_outcome[c,t], start:end]' *   to_vector(beta[c, t, 1:n_covs_per_outcome[c, t]]);
                                  if (c == 2)   Xbeta =  X_d[t,  1:n_covs_per_outcome[c,t], start:end]' *   to_vector(beta[c, t, 1:n_covs_per_outcome[c, t]]);
                                  Bound_Z  = - (Xbeta + inc  )  *  L_Omega_diag_recip[c, t] ; // use as marker for potential overflow
                             } else {
                                  Bound_Z  = - (beta[c, t, 1] + inc  )  *  L_Omega_diag_recip[c, t] ; // use as marker for potential overflow
                             }
               {

                    int num_OK_index = 0 ;
                    int num_Bound_Z_overflows_and_y_eq_1 = 0 ;
                    int num_Bound_Z_underflows_and_y_eq_0 = 0 ;

                    for (m in 1:M) {
                           if       ( (Bound_Z[m]  >  overflow_threshold)    &&  (y_s[m, t] == 1) )      num_Bound_Z_overflows_and_y_eq_1  += 1;
                           else if  ( (Bound_Z[m]  <  underflow_threshold)   &&  (y_s[m, t] == 0) )      num_Bound_Z_underflows_and_y_eq_0 += 1;
                           else   num_OK_index += 1;
                    }

                    if (num_OK_index == M)  { // carry on as normal as no * problematic * overflows/underflows

                           if (Phi_type == 2) {
                                 //// Phi_approx (Phi_type == 2) - done on the log-scale, since 1 - Phi_approx(Bz) rounds to 0 for large Bz (i.e., Bz > ~7.11 when y == 1),
                                 //// which gives log(0), a NaN gradient and Z = +Inf.
                                 //// Note that Phi_approx(x) = inv_logit(poly(x)) with poly(x) = 0.07056 x^3 + 1.5976 x, hence: log(Phi_approx(Bz)) = log_inv_logit(poly) and
                                 //// log(1 - Phi_approx(Bz)) = log_inv_logit(-poly), for every Bz.
                                 //// sign_y = +1 (y == 1) or -1 (y == 0) and poly_signed = sign_y * poly, so y1 = log_inv_logit(-poly_signed).
                                 //// Z = Phi_approx^{-1}(q) (via logit(q)). For y == 1: q = Phi + (1 - Phi) u and 1 - q = (1 - Phi)(1 - u); y == 0 is the mirror
                                 //// image (poly -> -poly, u -> 1 - u, logit q -> -logit q), so with v = u (y == 1) or v = 1 - u (y == 0):
                                 ////     logit(q) = sign_y * [ log_sum_exp(log_inv_logit(poly_signed), log_inv_logit(-poly_signed) + log(v)) - (log(1 - v) + log_inv_logit(-poly_signed)) ].
                                 //// log(v) and log(1 - v) come from log(u) and log1m(u) directly (so 1 - u is never computed).
                                 vector[M] yy = y_s[, t];   //// 0/1 binary outcomes of these rows (data)
                                 vector[M] sign_y = yy + (yy - 1.0);
                                 vector[M] poly_signed = sign_y .* (0.07056 * square(Bound_Z) .* Bound_Z + 1.5976 * Bound_Z);
                                 vector[M] log_Phi_signed = log_inv_logit(poly_signed);
                                 vector[M] log_1m_Phi_signed = log_inv_logit(-poly_signed);
                                 vector[M] log_u = log(u[, t]);
                                 vector[M] log_1m_u = log1m(u[, t]);
                                 vector[M] log_v;
                                 vector[M] log_1m_v;
                                 for (ii in 1:M) {
                                       if (yy[ii] == 1) { log_v[ii] = log_u[ii];    log_1m_v[ii] = log_1m_u[ii]; }
                                       else             { log_v[ii] = log_1m_u[ii]; log_1m_v[ii] = log_u[ii];    }
                                 }
                                 matrix[M, 2] tmp;
                                 tmp[, 1] = log_Phi_signed;
                                 tmp[, 2] = log_1m_Phi_signed + log_v;
                                 Z_std_norm[, t] = inv_Phi_approx_from_logit_prob(sign_y .* (log_sum_exp_2d(tmp) - (log_1m_v + log_1m_Phi_signed)));
                                 y1[, t] = log_1m_Phi_signed;
                           } else if (Phi_type == 0) {
                                 //// Phi_type == 0: exact Phi, on the log-scale, using 1 - Phi(x) = Phi(-x) (so 1 - Phi(x) is never computed by subtraction).
                                 //// Same function as the tail branches below, so the log-lik is continuous at the thresholds (see Phi_exact_binary_Z_and_log_lik).
                                 matrix[M, 2] Phi_exact_Z_and_log_lik = Phi_exact_binary_Z_and_log_lik(Bound_Z, y_s[, t], u[, t]);
                                 Z_std_norm[, t] = Phi_exact_Z_and_log_lik[, 1];
                                 y1[, t] = Phi_exact_Z_and_log_lik[, 2];
                           } else {
                                 //// Phi_type == 1 (i.e., "Phi" in BayesMVP): exact Phi, vectorised, again using 1 - Phi(x) = Phi(-x);
                                 //// continuous with the (exact) tail branches at the thresholds (see Phi_exact_binary_Z_and_log_lik_vectorised).
                                 matrix[M, 2] Phi_exact_Z_and_log_lik = Phi_exact_binary_Z_and_log_lik_vectorised(Bound_Z, y_s[, t], u[, t]);
                                 Z_std_norm[, t] = Phi_exact_Z_and_log_lik[, 1];
                                 y1[, t] = Phi_exact_Z_and_log_lik[, 2];
                           }

                    } else if (num_OK_index < M)  {

                                  int indicator_OK_empty = 0;
                                  if (num_OK_index < 1)  {
                                    num_OK_index = 1;
                                    indicator_OK_empty = 1;
                                  }
                                  ////
                                  int indicator_overflows_and_y_eq_1_empty = 0;
                                  if ( num_Bound_Z_overflows_and_y_eq_1  < 1)  {
                                    num_Bound_Z_overflows_and_y_eq_1  = 1;
                                    indicator_overflows_and_y_eq_1_empty = 1;
                                  }
                                  ////
                                  int indicator_underflows_and_y_eq_0_empty = 0;
                                  if (num_Bound_Z_underflows_and_y_eq_0 < 1)  {
                                    num_Bound_Z_underflows_and_y_eq_0  = 1;
                                    indicator_underflows_and_y_eq_0_empty = 1;
                                  }
                                  ////
                                  array[num_OK_index] int OK_index;
                                  array[num_Bound_Z_overflows_and_y_eq_1] int overflows_and_y_eq_1_index;
                                  array[num_Bound_Z_underflows_and_y_eq_0] int underflows_and_y_eq_0_index;
                                  int counter_1  = 1;
                                  int counter_2  = 1;
                                  int counter_3  = 1;
                                  ////
                                  for (m in 1:M) {
                                         if  (    (Bound_Z[m]  >  overflow_threshold)    &&  (y_s[m, t] == 1) ) {
                                             overflows_and_y_eq_1_index[counter_1] = m;
                                             counter_1 += 1;
                                         } else if  ( (Bound_Z[m]  <  underflow_threshold)   &&  (y_s[m, t] == 0) )  {
                                             underflows_and_y_eq_0_index[counter_2] = m;
                                             counter_2 += 1;
                                         } else {
                                             OK_index[counter_3] = m;
                                             counter_3 += 1;
                                        }
                                  }
                                  ////
                                  if (indicator_OK_empty == 0) {

                                           array[num_OK_index] int index = OK_index;
                                           int local_size = num_OK_index;

                                          if (Phi_type == 2) {
                                             //// Phi_approx (Phi_type == 2) - done on the log-scale, since 1 - Phi_approx(Bz) rounds to 0 for large Bz (i.e., Bz > ~7.11 when y == 1),
                                             //// which gives log(0), a NaN gradient and Z = +Inf.
                                             //// Note that Phi_approx(x) = inv_logit(poly(x)) with poly(x) = 0.07056 x^3 + 1.5976 x, hence: log(Phi_approx(Bz)) = log_inv_logit(poly) and
                                             //// log(1 - Phi_approx(Bz)) = log_inv_logit(-poly), for every Bz.
                                             //// sign_y = +1 (y == 1) or -1 (y == 0) and poly_signed = sign_y * poly, so y1 = log_inv_logit(-poly_signed).
                                             //// Z = Phi_approx^{-1}(q) (via logit(q)). For y == 1: q = Phi + (1 - Phi) u and 1 - q = (1 - Phi)(1 - u); y == 0 is the mirror
                                             //// image (poly -> -poly, u -> 1 - u, logit q -> -logit q), so with v = u (y == 1) or v = 1 - u (y == 0):
                                             ////     logit(q) = sign_y * [ log_sum_exp(log_inv_logit(poly_signed), log_inv_logit(-poly_signed) + log(v)) - (log(1 - v) + log_inv_logit(-poly_signed)) ].
                                             //// log(v) and log(1 - v) come from log(u) and log1m(u) directly (so 1 - u is never computed).
                                             vector[local_size] yy = y_s[index, t];   //// 0/1 binary outcomes of these rows (data)
                                             vector[local_size] sign_y = yy + (yy - 1.0);
                                             vector[local_size] poly_signed = sign_y .* (0.07056 * square(Bound_Z[index]) .* Bound_Z[index] + 1.5976 * Bound_Z[index]);
                                             vector[local_size] log_Phi_signed = log_inv_logit(poly_signed);
                                             vector[local_size] log_1m_Phi_signed = log_inv_logit(-poly_signed);
                                             vector[local_size] log_u = log(u[index, t]);
                                             vector[local_size] log_1m_u = log1m(u[index, t]);
                                             vector[local_size] log_v;
                                             vector[local_size] log_1m_v;
                                             for (ii in 1:local_size) {
                                                   if (yy[ii] == 1) { log_v[ii] = log_u[ii];    log_1m_v[ii] = log_1m_u[ii]; }
                                                   else             { log_v[ii] = log_1m_u[ii]; log_1m_v[ii] = log_u[ii];    }
                                             }
                                             matrix[local_size, 2] tmp;
                                             tmp[, 1] = log_Phi_signed;
                                             tmp[, 2] = log_1m_Phi_signed + log_v;
                                             Z_std_norm[index, t] = inv_Phi_approx_from_logit_prob(sign_y .* (log_sum_exp_2d(tmp) - (log_1m_v + log_1m_Phi_signed)));
                                             y1[index, t] = log_1m_Phi_signed;
                                          } else if (Phi_type == 0) {
                                             //// Phi_type == 0: exact Phi, on the log-scale, using 1 - Phi(x) = Phi(-x) (so 1 - Phi(x) is never computed by subtraction).
                                             //// Same function as the tail branches below, so the log-lik is continuous at the thresholds (see Phi_exact_binary_Z_and_log_lik).
                                             matrix[local_size, 2] Phi_exact_Z_and_log_lik = Phi_exact_binary_Z_and_log_lik(Bound_Z[index], y_s[index,t], u[index, t]);
                                             Z_std_norm[index, t] = Phi_exact_Z_and_log_lik[, 1];
                                             y1[index, t] = Phi_exact_Z_and_log_lik[, 2];
                                          } else {
                                             //// Phi_type == 1 (i.e., "Phi" in BayesMVP): exact Phi, vectorised, again using 1 - Phi(x) = Phi(-x);
                                             //// continuous with the (exact) tail branches at the thresholds (see Phi_exact_binary_Z_and_log_lik_vectorised).
                                             matrix[local_size, 2] Phi_exact_Z_and_log_lik = Phi_exact_binary_Z_and_log_lik_vectorised(Bound_Z[index], y_s[index, t], u[index, t]);
                                             Z_std_norm[index, t] = Phi_exact_Z_and_log_lik[, 1];
                                             y1[index, t] = Phi_exact_Z_and_log_lik[, 2];
                                          }

                                   }
                                   if (indicator_underflows_and_y_eq_0_empty ==  0) { /// underflow + y == 0

                                              array[num_Bound_Z_underflows_and_y_eq_0] int index = underflows_and_y_eq_0_index;
                                              int local_size = num_Bound_Z_underflows_and_y_eq_0;

                                              if (Phi_type != 2) {
                                                    //// Exact Phi (Phi_type == 0 or 1) for the lower tail (y == 0 and Bound_Z < underflow_threshold) - same function as the Phi_type == 0
                                                    //// branch above, so the log-lik, its gradient and Z are all continuous at underflow_threshold. Phi_type == 2 uses the Phi_approx tail below.
                                                    matrix[local_size, 2] Phi_exact_Z_and_log_lik = Phi_exact_binary_Z_and_log_lik(Bound_Z[index], rep_vector(0.0, local_size), u[index, t]);
                                                    Z_std_norm[index, t] = Phi_exact_Z_and_log_lik[, 1];
                                                    y1[index, t] = Phi_exact_Z_and_log_lik[, 2];
                                              } else {
                                              vector[local_size] log_Bound_U_Phi_Bound_Z =  log_inv_logit( 0.07056 * square(Bound_Z[index]) .* Bound_Z[index]  + 1.5976 * Bound_Z[index] );
                                              vector[local_size] log_Phi_Z = log(u[index, t]) +  log_Bound_U_Phi_Bound_Z ;
                                              vector[local_size] log_1m_Phi_Z =   log1m_exp(log(u[index, t])  + log_Bound_U_Phi_Bound_Z);
                                              vector[local_size] logit_Phi_Z = log_Phi_Z - log_1m_Phi_Z;
                                              Z_std_norm[index, t] = inv_Phi_approx_from_logit_prob(logit_Phi_Z);
                                              y1[index, t]  =  log_Bound_U_Phi_Bound_Z ;
                                              }

                                   }
                                   if (indicator_overflows_and_y_eq_1_empty == 0) {  //// overflow + y == 1

                                             array[num_Bound_Z_overflows_and_y_eq_1] int index = overflows_and_y_eq_1_index;
                                             int local_size = num_Bound_Z_overflows_and_y_eq_1;

                                             if (Phi_type != 2) {
                                                   //// Exact Phi (Phi_type == 0 or 1) for the upper tail (y == 1 and Bound_Z > overflow_threshold; log(1 - Phi) = log(Phi(-Bound_Z))) - same function
                                                   //// as the Phi_type == 0 branch above, so the log-lik, its gradient and Z are all continuous at overflow_threshold. Phi_type == 2 uses the Phi_approx tail below.
                                                   matrix[local_size, 2] Phi_exact_Z_and_log_lik = Phi_exact_binary_Z_and_log_lik(Bound_Z[index], rep_vector(1.0, local_size), u[index, t]);
                                                   Z_std_norm[index, t] = Phi_exact_Z_and_log_lik[, 1];
                                                   y1[index, t] = Phi_exact_Z_and_log_lik[, 2];
                                             } else {
                                             vector[local_size] log_Bound_U_Phi_Bound_Z_1m =  log_inv_logit( - 0.07056 * square(Bound_Z[index]) .* Bound_Z[index]  - 1.5976 * Bound_Z[index] );
                                             {
                                               matrix[num_Bound_Z_overflows_and_y_eq_1, 2] tmp_array_2d_to_lse;
                                               tmp_array_2d_to_lse[, 1] = log_Bound_U_Phi_Bound_Z_1m + log(u[index, t]);
                                               vector[local_size] log_Bound_U_Phi_Bound_Z = log1m_exp(log_Bound_U_Phi_Bound_Z_1m);
                                               tmp_array_2d_to_lse[, 2] =  log_Bound_U_Phi_Bound_Z;
                                               vector[local_size] log_Phi_Z = log_sum_exp_2d(tmp_array_2d_to_lse);

                                               vector[local_size] log_1m_Phi_Z  =   log1m(u[index, t])  + log_Bound_U_Phi_Bound_Z_1m;
                                               vector[local_size] logit_Phi_Z = log_Phi_Z - log_1m_Phi_Z;
                                               Z_std_norm[index, t] = inv_Phi_approx_from_logit_prob(logit_Phi_Z);
                                             }
                                             y1[index, t]  =  log_Bound_U_Phi_Bound_Z_1m ;
                                             }

                                   }
                            }

                         }

                            if (t < n_tests)   inc = block(Z_std_norm, 1, 1, M, t) * to_vector(head(L_Omega[c, t + 1, ], t))   ;

                    }  // end of t loop

                    lp[, c] = to_vector(rowwise_sum(y1[, 1:n_tests]))  +  to_vector(log_prev[, c])  ;

              } // end of c loop

              out += sum(log_sum_exp_2d(lp));

              return out;

      }

      //////////////////////////////////////////////////////////////////////////
      //// ---- Exact normal GHK helpers (Phi_type == 0; also used for the tail branches and per-subject loop of Phi_type == 1):
      ////
      //// These use the exact standard normal CDF in every branch (ordinary AND tail), so the log-lik, its gradient and Z are all
      //// continuous at the thresholds. Note that using the Phi_approx tails instead would make the log-lik jump at +/-7.5
      //// by 10.67 nats (log Phi(-7.5) = -31.08 vs. log Phi_approx(-7.5) = -41.75).
      //// The upper side is obtained by reflection (i.e., 1 - Phi(x) = Phi(-x)) on the log-scale, since computing 1 - Phi(x) by subtraction
      //// loses accuracy near +7.5 (and can give Z = +Inf). Z is taken from whichever of log(q) / log(1 - q) is the smaller side.
      //// log_Phi_stable, inv_Phi_from_log_lower and normal_from_log_uniform are copied from the NicoStan four_class models.
      //// The vectorised ordinary branches of Phi_type == 1 use Phi_exact_binary_Z_and_log_lik_vectorised / Phi_exact_interval_Z_and_log_lik_vectorised (below).
      //// Phi_type == 2 never calls these functions.
      ////
      real log_Phi_stable(real x) {
                //// Differentiate the SAME expression used for the log-CDF value. This avoids
                //// the separate approximate derivative in the bundled std_normal_lcdf routine.
                if (x < -10.0) {
                    //// Laplace continued fraction for the Mills ratio Phi(-t) / phi(t).
                    //// Twenty terms agree to double precision in this t > 10 branch.
                    real t = -x;
                    real r = 0.0;
                    for (k in 1:20) r = (21.0 - k) / (t + r);
                    return -0.5 * square(x) - 0.5 * log(2.0 * pi()) - log(t + r);
                }
                if (x <= 0.0) return log(erfc(-x / sqrt(2.0))) - log(2.0);
                return log1m(0.5 * erfc(x / sqrt(2.0)));
      }
      real inv_Phi_from_log_lower(real log_p) {
                //// Called only for a probability <= 0.5. Exponentiating is safe down to -700.
                if (log_p > -700.0) return inv_Phi(exp(log_p));
                real z = -sqrt(-2.0 * log_p);
                for (i in 1:5) {
                    real lp = log_Phi_stable(z);
                    z -= (lp - log_p) * exp(lp - std_normal_lpdf(z));
                }
                return z;
      }
      real normal_from_log_uniform(real log_u, real log1m_u) {
                if (log_u <= -log(2.0)) return inv_Phi_from_log_lower(log_u);
                return -inv_Phi_from_log_lower(log1m_u);
      }
      ////
      //// Binary GHK step (exact normal): y == 1 -> Z truncated to (Bound_Z, Inf), y == 0 -> Z truncated to (-Inf, Bound_Z).
      //// Returns [Z_std_norm, log-likelihood contribution].
      ////
      vector Phi_exact_binary_step(real Bound_Z, real y_obs, real u_obs) {
                vector[2] Z_and_log_lik;
                real log_q;
                real log_1m_q;
                if (y_obs == 1) {
                    //// log(1 - Phi(Bound_Z)) = log Phi(-Bound_Z);  q = Phi(Bound_Z) + (1 - Phi(Bound_Z)) u,  1 - q = (1 - Phi(Bound_Z)) (1 - u).
                    Z_and_log_lik[2] = log_Phi_stable(-Bound_Z);
                    log_1m_q = log1m(u_obs) + Z_and_log_lik[2];
                    log_q = log1m_exp(log_1m_q);
                } else {
                    //// log Phi(Bound_Z);  q = Phi(Bound_Z) u.
                    Z_and_log_lik[2] = log_Phi_stable(Bound_Z);
                    log_q = log(u_obs) + Z_and_log_lik[2];
                    log_1m_q = log1m_exp(log_q);
                }
                Z_and_log_lik[1] = normal_from_log_uniform(log_q, log_1m_q);
                return Z_and_log_lik;
      }
      ////
      //// Interval GHK step (exact normal): Z truncated to (Bound_Z_lo, Bound_Z_hi). Returns [Z_std_norm, log-likelihood contribution].
      //// Bound_Z_lo > 0 uses the reflected (right-side) form, so neither Phi(hi) - Phi(lo) nor Phi(lo) + prob * u is formed near 1.
      ////
      vector Phi_exact_interval_step(real Bound_Z_lo, real Bound_Z_hi, real u_obs) {
                vector[2] Z_and_log_lik;
                real log_q;
                real log_1m_q;
                if (Bound_Z_lo > 0) {
                    //// P = (1 - Phi(lo)) - (1 - Phi(hi)),  1 - q = (1 - u)(1 - Phi(lo)) + u (1 - Phi(hi))  (sum of positive terms).
                    real log_1m_Phi_lo = log_Phi_stable(-Bound_Z_lo);
                    real log_1m_Phi_hi = log_Phi_stable(-Bound_Z_hi);
                    Z_and_log_lik[2] = log_diff_exp(log_1m_Phi_lo, log_1m_Phi_hi);
                    log_1m_q = log_sum_exp(log_1m_Phi_lo + log1m(u_obs), log_1m_Phi_hi + log(u_obs));
                    log_q = log1m_exp(log_1m_q);
                } else {
                    //// P = Phi(hi) - Phi(lo),  q = (1 - u) Phi(lo) + u Phi(hi)  (sum of positive terms; 1 - q >= (1 - u) / 2 here).
                    real log_Phi_lo = log_Phi_stable(Bound_Z_lo);
                    real log_Phi_hi = log_Phi_stable(Bound_Z_hi);
                    Z_and_log_lik[2] = log_diff_exp(log_Phi_hi, log_Phi_lo);
                    log_q = log_sum_exp(log_Phi_lo + log1m(u_obs), log_Phi_hi + log(u_obs));
                    log_1m_q = log1m_exp(log_q);
                }
                Z_and_log_lik[1] = normal_from_log_uniform(log_q, log_1m_q);
                return Z_and_log_lik;
      }
      ////
      //// Vector wrappers: return an [n, 2] matrix, column 1 = Z_std_norm, column 2 = log-likelihood contribution y1.
      ////
      matrix Phi_exact_binary_Z_and_log_lik(vector Bound_Z, vector y_vec, vector u_vec) {
                int n_obs = num_elements(Bound_Z);
                matrix[n_obs, 2] Z_and_log_lik;
                for (obs_index in 1:n_obs) {
                    Z_and_log_lik[obs_index, ] = Phi_exact_binary_step(Bound_Z[obs_index], y_vec[obs_index], u_vec[obs_index])';
                }
                return Z_and_log_lik;
      }
      matrix Phi_exact_interval_Z_and_log_lik(vector Bound_Z_lo, vector Bound_Z_hi, vector u_vec) {
                int n_obs = num_elements(Bound_Z_lo);
                matrix[n_obs, 2] Z_and_log_lik;
                for (obs_index in 1:n_obs) {
                    Z_and_log_lik[obs_index, ] = Phi_exact_interval_step(Bound_Z_lo[obs_index], Bound_Z_hi[obs_index], u_vec[obs_index])';
                }
                return Z_and_log_lik;
      }

      //////////////////////////////////////////////////////////////////////////
      //// ---- VECTORISED exact normal GHK helpers for Phi_type == 1:
      ////
      //// Phi_type == 1 is BayesMVP's exact "Phi" (i.e., Phi_type = "Phi" in R). Its tail branches and per-subject loop use the exact
      //// functions above (same as Phi_type == 0), and its ordinary branches use the two vectorised functions below.
      //// Note that using the Phi_approx tails here would make the log-lik jump at +/-7.5 (by up to 1.26 nats).
      ////
      //// Both use the reflection 1 - Phi(x) = Phi(-x), so every CDF value that is formed is a LOWER-side value Phi(w) = 0.5 * erfc(-w / sqrt(2))
      //// (erfc keeps full relative accuracy for w < 0, unlike 0.5 * (1 + erf)), and 1 - Phi(x) is never formed by subtraction.
      ////   Binary, sign_y = +1 (y == 1) or -1 (y == 0): likelihood Phi(w) with w = -sign_y * Bound_Z, and
      ////     y == 0:  q     = u Phi(Bound_Z)          ->  Z =  inv_Phi(u Phi(w))
      ////     y == 1:  1 - q = (1 - u) Phi(-Bound_Z)   ->  Z = -inv_Phi((1 - u) Phi(w))       (Phi^{-1}(q) = -Phi^{-1}(1 - q))
      ////   Interval (lo, hi): when lo > 0 the interval is reflected to (a, b) = (-hi, -lo) with u -> 1 - u and Z -> -Z, otherwise (a, b) = (lo, hi).
      ////     After reflection a <= 0, P = Phi(b) - Phi(a) and q = (1 - v) Phi(a) + v Phi(b) is a sum of positive terms.
      //// Only the sign / weight selection is an element loop (on bound VALUES, no autodiff); erfc, log and inv_Phi are vectorised calls.
      ////
      //// Precondition (checked in transformed data): these routines are ONLY called from the ordinary branches, where w (binary) and b
      //// (interval) are >= min(underflow_threshold, -overflow_threshold) >= -35, so Phi(w), Phi(b) >= 1e-268 and never underflow (a may be
      //// very negative: Phi(a) = 0 then, with a zero, not NaN, derivative). At the thresholds log(0.5 * erfc(-w / sqrt(2))) equals
      //// log_Phi_stable(w) = log(erfc(-w / sqrt(2))) - log(2) to rounding, so the Phi_type == 1 target is continuous there.
      //// Phi_type == 0 and Phi_type == 2 never call these two functions.
      ////
      matrix Phi_exact_binary_Z_and_log_lik_vectorised(vector Bound_Z, vector y_vec, vector u_vec) {
                int n_obs = num_elements(Bound_Z);
                vector[n_obs] sign_y = y_vec + (y_vec - 1.0);
                //// Phi(w), w = -sign_y * Bound_Z:  0.5 * erfc(-w / sqrt(2)) = 0.5 * erfc(sign_y * Bound_Z / sqrt(2)).
                vector[n_obs] Phi_w = 0.5 * erfc((sign_y .* Bound_Z) / sqrt(2.0));
                //// v = u (y == 0) or 1 - u (y == 1); 1 - u is exact for u >= 0.5 and has relative error <= 1.1e-16 otherwise.
                vector[n_obs] v_vec = (1.0 - y_vec) .* u_vec + y_vec .* (1.0 - u_vec);
                matrix[n_obs, 2] Z_and_log_lik;
                Z_and_log_lik[, 1] = -sign_y .* inv_Phi(v_vec .* Phi_w);
                Z_and_log_lik[, 2] = log(Phi_w);
                return Z_and_log_lik;
      }
      matrix Phi_exact_interval_Z_and_log_lik_vectorised(vector Bound_Z_lo, vector Bound_Z_hi, vector u_vec) {
                int n_obs = num_elements(Bound_Z_lo);
                //// reflect_weight = 1 when Bound_Z_lo > 0 (reflected right-side form), 0 otherwise; keep_weight = 1 - reflect_weight.
                array[n_obs] int reflect_indicator;
                for (obs_index in 1:n_obs) reflect_indicator[obs_index] = (Bound_Z_lo[obs_index] > 0);
                vector[n_obs] reflect_weight = to_vector(reflect_indicator);
                vector[n_obs] keep_weight = 1.0 - reflect_weight;
                vector[n_obs] reflect_sign = keep_weight - reflect_weight;
                vector[n_obs] Bound_a = keep_weight .* Bound_Z_lo - reflect_weight .* Bound_Z_hi;
                vector[n_obs] Bound_b = keep_weight .* Bound_Z_hi - reflect_weight .* Bound_Z_lo;
                vector[n_obs] v_vec    = keep_weight .* u_vec + reflect_weight .* (1.0 - u_vec);
                vector[n_obs] v_vec_1m = keep_weight .* (1.0 - u_vec) + reflect_weight .* u_vec;
                vector[n_obs] Phi_a = 0.5 * erfc(-Bound_a / sqrt(2.0));
                vector[n_obs] Phi_b = 0.5 * erfc(-Bound_b / sqrt(2.0));
                matrix[n_obs, 2] Z_and_log_lik;
                Z_and_log_lik[, 1] = reflect_sign .* inv_Phi(v_vec_1m .* Phi_a + v_vec .* Phi_b);
                Z_and_log_lik[, 2] = log(Phi_b - Phi_a);
                return Z_and_log_lik;
      }

}



data {
      int<lower=1> N;
      int<lower=2> n_tests;
      matrix<lower=0>[N, n_tests]   y;  //////// data
      int<lower=2> n_class;
      int<lower=1> n_pops;
      array[N] int pop;
      ////
      int n_covariates_max_nd;
      int n_covariates_max_d;
      int n_covariates_max;
      array[n_tests] matrix[N, n_covariates_max_nd] X_nd; /////// covariate array (can have  DIFFERENT NUMBERS of covariates for each  outcome - fill rest of array with 999999 if they vary between outcomes)
      array[n_tests] matrix[N, n_covariates_max_d]  X_d; /////// covariate array (can have  DIFFERENT NUMBERS of covariates for each  outcome - fill rest of array with 999999 if they vary between outcomes)
      array[n_class, n_tests] int n_covs_per_outcome;
      ////
      int corr_force_positive;
      int<lower=0, upper=(n_tests * (n_tests - 1)) %/% 2> known_num;
      real overflow_threshold;
      real underflow_threshold;
      ///// priors
      int prior_only;
      array[n_class] matrix[n_covariates_max, n_tests] prior_beta_mean;
      array[n_class] matrix<lower=0>[n_covariates_max, n_tests] prior_beta_sd;
      matrix<lower=0>[n_class, 1] prior_LKJ; // NOTE: Some Stan vector's written as mtx. w/ 1 col to avoid issues w/ custom C++ fns
      matrix<lower=0>[n_pops, 1] prior_p_alpha;
      matrix<lower=0>[n_pops, 1] prior_p_beta;
      ///// other
      int Phi_type;
      int handle_numerical_issues; // ACCEPTED BUT IGNORED (hni == 1 path is hard-wired)
      int fully_vectorised;        // ACCEPTED BUT IGNORED (fully_vectorised == 1 path is hard-wired)
      ////
      int<lower=1> chunk_size;      // Static partial-sum size upper bound (observations), independent of WCP threads
}

transformed data {
      ////
      //// ---- Checks on Phi_type and the thresholds:
      //// Phi_type must be 0 (exact, scalar), 1 (BayesMVP's exact "Phi", vectorised ordinary branch) or 2 (Phi_approx).
      //// The Phi_type == 1 vectorised ordinary branch needs -35 <= underflow_threshold and overflow_threshold <= 35
      //// (see Phi_exact_binary_Z_and_log_lik_vectorised).
      ////
      if (Phi_type < 0 || Phi_type > 2) reject("Phi_type must be 0, 1 or 2; got Phi_type = ", Phi_type);
      if (Phi_type == 1 && (underflow_threshold < -35.0 || overflow_threshold > 35.0))
          reject("Phi_type = 1 needs -35 <= underflow_threshold and overflow_threshold <= 35; got underflow_threshold = ", underflow_threshold,
                 ", overflow_threshold = ", overflow_threshold);

      int k_choose_2 = (n_tests * (n_tests - 1)) / 2;
      int km1_choose_2 = ((n_tests - 1) * (n_tests - 2)) / 2;

      int n_covariates_total_nd =    (sum( (n_covs_per_outcome[1,])));
      int n_covariates_total_d =     (sum( (n_covs_per_outcome[2,])));
      int n_covariates_total =       n_covariates_total_nd + n_covariates_total_d;

      real s = 1 / 1.702;
      real a = 0.07056;
      real b = 1.5976;
      real a_times_3 = 3.0 * 0.07056;
      real<lower=-1, upper=1> lb;
      real<lower=lb, upper=1> ub = 1.0;

      if (corr_force_positive == 1)  lb = 0;
      else lb = -1.0;

}

parameters {

       array[N] row_vector[n_tests] u_raw; // CHANGED (was matrix[N, n_tests]): array-of-rows so reduce_sum slices the nuisance parameters directly
       array[n_class] vector[n_tests - 1] col_one_raw;
       array[n_class] vector[km1_choose_2 - known_num] off_raw;
       vector[n_covariates_total] beta_vec;
       vector[n_pops]  p_raw;

}

transformed parameters {

     array[n_class, n_tests, n_covariates_max] real beta;
     vector<lower=0, upper=1>[n_pops]   prev = lb_ub_lp(p_raw, 0.0, 1.0);
     array[n_class] matrix[n_tests, n_tests] Omega;
     array[n_class] matrix[n_tests, n_tests] L_Omega;
     matrix[n_class, n_tests] L_Omega_diag_recip;

      {
            int counter = 1;
            for (c in 1 : n_class) {
                      for (t in 1:n_tests) {
                        for (k in 1:n_covs_per_outcome[c, t]) {
                           beta[c, t, k] = beta_vec[counter];
                           counter += 1;
                        }
                      }
                    L_Omega[c, :  ] =   cholesky_corr_constrain_outer_lp( to_vector(col_one_raw[c, :]), to_vector(off_raw[c, :]), lb, ub);
                    Omega[c, :  ] = multiply_lower_tri_self_transpose(L_Omega[c, :]);
                    L_Omega_diag_recip[c, ] = to_row_vector(1.0 ./ diagonal(L_Omega[c, :  ]));
            }
      }

}

model {

              for (c in 1 : n_class) {
                  for (t in 1 : n_tests) {
                       for (k in 1 : n_covs_per_outcome[c, t]) {
                         beta[c, t, k] ~ normal(prior_beta_mean[c, k, t], prior_beta_sd[c, k, t]);
                      }
                  }
                   target += lkj_corr_cholesky_lpdf(L_Omega[c,,] | prior_LKJ[c, 1]) ;
              }

              for (g in 1 : n_pops) {
                prev[g] ~ beta(prior_p_alpha[g, 1], prior_p_beta[g, 1]);
              }

              if (prior_only == 0) {
                    //// Deterministic partitioning; chunk_size bounds slice size, not the number of threads.
                    target += reduce_sum_static(partial_log_lik,
                                         u_raw,
                                         chunk_size,
                                         y,
                                         X_nd,
                                         X_d,
                                         n_covs_per_outcome,
                                         pop,
                                         beta,
                                         L_Omega,
                                         L_Omega_diag_recip,
                                         prev,
                                         n_tests,
                                         n_class,
                                         n_covariates_max,
                                         Phi_type,
                                         overflow_threshold,
                                         underflow_threshold);
              }

}

generated quantities {

    vector[n_tests] Se_bin;
    vector[n_tests] Sp_bin;
    vector[n_tests] Fp_bin;
    vector<lower=0, upper=1>[n_pops] p = prev;

   for (c in 1:n_class) {

      for (t in 1:n_tests) { // for binary tests

         if (n_class == 2) { // summary Se and Sp only calculated if n_class = 2 (i.e. the "standard" # of classes for DTA)
              Se_bin[t]  =        Phi(   beta[2, t, 1]   );
              Sp_bin[t]  =    1 - Phi(   beta[1, t, 1]   );
              Fp_bin[t]  =    1 - Sp_bin[t];
        }
        else {
          Se_bin[t] = 999;
          Sp_bin[t] = 999;
          Fp_bin[t] = 999;
        }
    }

}

}

