#include <tf.h>

/* NEW: These functions compute b = D*a or b = D^T*a,
 * for D the difference operator of order k, defined
 * over the inputs x
 * See "The Falling Factorial Basis" paper for pseudocode
 * of how these can be implemented directly in O(nk) time */
void tf_dx(double *x, int n, int k,double *a, double *b)
{
    int i=0, j=0;

    memcpy(b, a, n);
    for(i=0; i < k; ++i)
    {
      if( i != 0 )
      {
        /* b[i:n-1] = b[i:n-1] ./ ( x[i:n-1] - x[0:n-1-i] ) */
        for(j=i; j < n; ++j)
        {
          b[j] = b[j] / ( x[j] - x[j-i]);
        }

        /* this was diff */
        /* b[i+1:n-1] = -b[i:n-2] + b[i+1:n-1] */
        for(j=n-1; j >= i+1; --j)
        {
          b[j] = b[j] - b[j-1];
        }
      }
    }
    /* TODO mult(b, factorial(k-1)); */
}
