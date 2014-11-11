#include "tf.h"

/* This function computes b = D*a,
 * for D the difference operator of order k, defined
 * over the inputs x
 * See "The Falling Factorial Basis" paper for pseudocode
 * of how these can be implemented directly in O(nk) time
 * Note: b should have size n, not n-k. The result is in the first n-k elements.
 */
void tf_dx(double *x, int n, int k,double *a, double *b)
{
  int i, j;

  for(i=0; i < n; i++) b[i] = a[i];

  if( k < 1 || k >= n )
    return;

  for(i=0; i < k; ++i)
  {
    if( i != 0 )
    {
      /* b[i:n-1] = b[i:n-1] ./ ( x[i:n-1] - x[0:n-1-i] ) */
      for(j=i; j < n; ++j)
      {
        b[j] = b[j] / ( x[j] - x[j-i]);
      }
    }

    /* b[i+1:n-1] = -b[i:n-2] + b[i+1:n-1] */
    for(j=n-1; j >= i+1; --j)
    {
      b[j] = b[j] - b[j-1];
    }
  }

  double fact = glmgen_factorial(k-1);
  for(i=0; i < n; ++i)
  {
    b[i] *= fact;
  }

  /* move the solution to the beginning of the array */
  memmove(b, b+k, (n-k)*sizeof(double));
}
