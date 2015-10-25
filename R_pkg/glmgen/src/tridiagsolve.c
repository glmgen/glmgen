

/*
 * From
 * https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm
 */

#include "tridiagsolve.h"
#include "stdlib.h"

void tridiagsolve(int N, const double *a, 
  const double *b, const double *c, double *x, double *cprime) 
{
  
  int in;
    
  cprime[0] = c[0] / b[0];
  x[0] = x[0] / b[0];
  
  /* loop from 1 to N - 1 inclusive */
  for (in = 1; in < N; in++) {
      double m = 1.0 / (b[in] - a[in-1] * cprime[in - 1]);
      cprime[in] = c[in] * m;
      x[in] = (x[in] - a[in-1] * x[in - 1]) * m;
  }
  
  /* loop from N - 2 to 0 inclusive, safely testing loop end condition */
  for (in = N - 1; in-- > 0; )
      x[in] = x[in] - cprime[in] * x[in + 1];
}
