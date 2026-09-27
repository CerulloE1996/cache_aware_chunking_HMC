/*
 * calibrate_cache_levels: random pointer chase over a working set of a given size, to check that the fill counters in
 * pmc_stat attribute loads to the expected level (small working set -> L2, medium -> L3, large -> DRAM).
 * Usage: calibrate_cache_levels <working_set_kib> <n_steps>
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int main(int argc, char **argv) {
    if (argc < 3) return 2;
    size_t bytes = (size_t)atoll(argv[1]) * 1024, n = bytes / 64;   /* one pointer per 64-byte cache line */
    long long steps = atoll(argv[2]);
    size_t *lines = aligned_alloc(64, n * 64), *order = malloc(n * sizeof(size_t));
    for (size_t i = 0; i < n; ++i) order[i] = i;
    srand(1);
    for (size_t i = n - 1; i > 0; --i) { size_t j = (size_t)rand() % (i + 1); size_t t = order[i]; order[i] = order[j]; order[j] = t; }
    for (size_t i = 0; i < n; ++i) lines[order[i] * 8] = order[(i + 1) % n] * 8;
    size_t p = 0;
    for (long long s = 0; s < steps; ++s) p = lines[p];
    printf("%zu\n", p);
    return 0;
}
