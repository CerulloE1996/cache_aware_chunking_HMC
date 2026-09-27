#include <stdio.h>
#include <stdlib.h>
int main(void) { size_t n = (size_t)1 << 27; double *a = malloc(n * sizeof(double)); for (size_t i = 0; i < n; ++i) a[i] = (double)i;
  double s = 0; for (int r = 0; r < 5; ++r) for (size_t i = 0; i < n; ++i) s += a[i]; printf("%f\n", s); return 0; }
