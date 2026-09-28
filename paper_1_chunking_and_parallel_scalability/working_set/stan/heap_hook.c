#define _GNU_SOURCE
#include <malloc.h>
#include <stddef.h>
//// ---- malloc interposition (heap current and peak, in usable bytes) -------------------------
void* __libc_malloc(size_t);
void __libc_free(void*);
void* __libc_realloc(void*, size_t);
void* __libc_calloc(size_t, size_t);
void* __libc_memalign(size_t, size_t);
static long long g_heap_current = 0;
static long long g_heap_peak = 0;
static inline void heap_add(void* p) {
  if (!p) return;
  long long n = (long long)(malloc_usable_size(p));
  long long c = __atomic_add_fetch(&g_heap_current, n, __ATOMIC_RELAXED);
  long long pk = __atomic_load_n(&g_heap_peak, __ATOMIC_RELAXED);
  while (c > pk && !__atomic_compare_exchange_n(&g_heap_peak, &pk, c, 1, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) {}
}
static inline void heap_sub(void* p) {
  if (!p) return;
  long long n = (long long)(malloc_usable_size(p));
  __atomic_sub_fetch(&g_heap_current, n, __ATOMIC_RELAXED);
}
void* malloc(size_t n) { void* p = __libc_malloc(n); heap_add(p); return p; }
void free(void* p) { heap_sub(p); __libc_free(p); }
void* calloc(size_t a, size_t b) { void* p = __libc_calloc(a, b); heap_add(p); return p; }
void* realloc(void* p, size_t n) {
  if (!p) return malloc(n);
  if (n == 0) { free(p); return NULL; }
  long long old = (long long)(malloc_usable_size(p));
  void* q = __libc_realloc(p, n);
  if (q) { __atomic_sub_fetch(&g_heap_current, old, __ATOMIC_RELAXED); heap_add(q); }
  return q;
}
void* memalign(size_t al, size_t n) { void* p = __libc_memalign(al, n); heap_add(p); return p; }
void* aligned_alloc(size_t al, size_t n) { void* p = __libc_memalign(al, n); heap_add(p); return p; }
int posix_memalign(void** pp, size_t al, size_t n) {
  void* p = __libc_memalign(al, n);
  if (!p) return 12;
  heap_add(p); *pp = p; return 0;
}
long long ws_heap_current() { return __atomic_load_n(&g_heap_current, __ATOMIC_RELAXED); }
long long ws_heap_peak() { return __atomic_load_n(&g_heap_peak, __ATOMIC_RELAXED); }
void ws_heap_reset_peak() { __atomic_store_n(&g_heap_peak, __atomic_load_n(&g_heap_current, __ATOMIC_RELAXED), __ATOMIC_RELAXED); }
