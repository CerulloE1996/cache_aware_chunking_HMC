////
//// Stan working-set harness: evaluates the log-density gradient of the Paper 1
//// reduce_sum_static LC-MVP model exactly as BridgeStan's bs_model::log_density_gradient
//// does (stan::math::gradient on log_prob_propto_jacobian), single-threaded, and records
//// the reverse-mode arena, AD-stack and heap use of every partial-sum (chunk) call.
////
#include <malloc.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <new>
#include <string>
#include <vector>

#include "model_rs_probed.hpp"
#include <stan/io/json/json_data.hpp>
#include <tbb/global_control.h>

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "usage: ws_harness data.json init.json label [n_reps]\n";
    return 1;
  }
  const std::string data_path = argv[1], init_path = argv[2], label = argv[3];
  const int n_reps = (argc > 4) ? std::atoi(argv[4]) : 3;

  //// Single thread for all TBB work (reduce_sum_static blocks run sequentially on this thread).
  tbb::global_control single_thread(tbb::global_control::max_allowed_parallelism, 1);

  //// Replace the main arena by one with a single 8 GB first block (lazily committed), so that
  //// arena use can be read as a plain pointer difference (no block switches).
  {
    auto& ma = stan::math::ChainableStack::instance_->memalloc_;
    ma.~stack_alloc();
    new (&ma) stan::math::stack_alloc(static_cast<size_t>(8) << 30);
  }

  std::ifstream data_stream(data_path);
  stan::json::json_data data_context(data_stream);
  model_rs_model_namespace::model_rs_model model(data_context, 0, &std::cerr);

  std::ifstream init_stream(init_path);
  stan::json::json_data init_context(init_stream);
  std::vector<int> params_i;
  std::vector<double> params_r_vec;
  model.transform_inits(init_context, params_i, params_r_vec, &std::cerr);
  Eigen::VectorXd params_r = Eigen::Map<Eigen::VectorXd>(params_r_vec.data(), params_r_vec.size());
  const long long P = params_r.size();

  ws_probe::state().chunks.reserve(100000);
  std::vector<double> grad(P);

  for (int rep = 0; rep < n_reps; ++rep) {
    ws_probe::state().chunks.clear();
    ws_probe::state().active = true;
    double lp = 0.0;
    ws_probe::stack_state s_begin, s_forward;
    {
      //// identical to stan::math::gradient(), with two probes added
      stan::math::nested_rev_autodiff nested;
      s_begin = ws_probe::snapshot();
      Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> x_var(params_r);
      stan::math::var fx_var = model.log_prob_propto_jacobian(x_var, &std::cerr);
      s_forward = ws_probe::snapshot();
      lp = fx_var.val();
      stan::math::grad(fx_var.vi_);
      for (long long i = 0; i < P; ++i) grad[i] = x_var.coeff(i).adj();
    }
    ws_probe::state().active = false;

    const auto& st = ws_probe::state();
    const long long outer_to_first_enter = st.chunks.empty() ? 0 : (st.chunks.front().enter.arena - st.outer.arena);
    double gnorm = 0.0;
    for (double g : grad) gnorm += g * g;
    std::printf("RUN label=%s rep=%d P=%lld lp=%.10g grad_l2=%.10g n_chunks=%zu "
                "arena_total_forward=%lld arena_begin_to_outer=%lld arena_outer_to_first_enter=%lld "
                "varstack_total_forward=%zu varnochain_total_forward=%zu heap_begin_to_forward=%lld\n",
                label.c_str(), rep, P, lp, std::sqrt(gnorm), st.chunks.size(),
                static_cast<long long>(s_forward.arena - s_begin.arena),
                static_cast<long long>(st.outer.arena - s_begin.arena),
                outer_to_first_enter,
                s_forward.var_stack - s_begin.var_stack,
                s_forward.var_nochain_stack - s_begin.var_nochain_stack,
                s_forward.heap - s_begin.heap);
    for (size_t k = 0; k < st.chunks.size(); ++k) {
      const auto& c = st.chunks[k];
      std::printf("CHUNK label=%s rep=%d k=%zu rows=%d enter_offset=%lld arena_fn=%lld var_stack_fn=%zu "
                  "var_nochain_fn=%zu var_alloc_fn=%zu heap_at_enter_minus_outer=%lld heap_peak_fn=%lld heap_net_fn=%lld "
                  "var_nochain_at_enter_minus_outer=%zu\n",
                  label.c_str(), rep, k, c.rows,
                  static_cast<long long>(c.enter.arena - st.outer.arena),
                  static_cast<long long>(c.exit.arena - c.enter.arena),
                  c.exit.var_stack - c.enter.var_stack,
                  c.exit.var_nochain_stack - c.enter.var_nochain_stack,
                  c.exit.var_alloc_stack - c.enter.var_alloc_stack,
                  c.enter.heap - st.outer.heap,
                  c.heap_peak_forward - c.enter.heap,
                  c.exit.heap - c.enter.heap,
                  c.enter.var_nochain_stack - st.outer.var_nochain_stack);
    }
  }
  return 0;
}
