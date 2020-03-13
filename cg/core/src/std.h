
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>

#define print_vector(V, n) \
  { printf("%s=[%lf", #V, V[0]); for(int __i=1; __i<n; __i++) printf(",%lf", V[__i]); printf("]\n"); }

#define print_matrix(M, m, n) \
  { printf("%s=\n", #M); for(int __i=0; __i<m; __i++) { \
    for(int __j=0; __j<n; __j++) printf(" %10.2lf", M[__i*n + __j]); printf("\n");}};
