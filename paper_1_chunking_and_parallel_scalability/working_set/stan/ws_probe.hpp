#pragma once
////
//// Working-set probes for the reduce_sum_static partial-sum functor.
//// Heap counters are maintained by the malloc interposition in ws_harness.cpp.
////
#include <stan/math/rev/core/chainablestack.hpp>
#include <cstddef>
#include <vector>

extern "C" long long ws_heap_current();
extern "C" long long ws_heap_peak();
extern "C" void ws_heap_reset_peak();

namespace ws_probe {

struct stack_state {
  const char* arena = nullptr;
  long long heap = 0;
  size_t var_stack = 0;
  size_t var_nochain_stack = 0;
  size_t var_alloc_stack = 0;
};

inline stack_state snapshot() {
  auto* s = stan::math::ChainableStack::instance_;
  stack_state out;
  //// alloc(0) returns the next free arena address without advancing it.
  out.arena = static_cast<const char*>(s->memalloc_.alloc(0));
  out.heap = ws_heap_current();
  out.var_stack = s->var_stack_.size();
  out.var_nochain_stack = s->var_nochain_stack_.size();
  out.var_alloc_stack = s->var_alloc_stack_.size();
  return out;
}

struct chunk_record {
  int rows = 0;
  stack_state enter;
  stack_state exit;
  long long heap_peak_forward = 0;
};

struct probe_state {
  bool active = false;
  stack_state outer;
  std::vector<chunk_record> chunks;
};

inline probe_state& state() {
  static probe_state st;
  return st;
}

inline void on_outer() {
  if (!state().active) return;
  state().outer = snapshot();
}

inline void on_enter(int rows) {
  if (!state().active) return;
  state().chunks.emplace_back();  //// capacity is reserved by the harness
  chunk_record& r = state().chunks.back();
  r.rows = rows;
  r.enter = snapshot();
  ws_heap_reset_peak();
}

inline void on_exit() {
  if (!state().active) return;
  chunk_record& r = state().chunks.back();
  r.exit = snapshot();
  r.heap_peak_forward = ws_heap_peak();
}

}  // namespace ws_probe
