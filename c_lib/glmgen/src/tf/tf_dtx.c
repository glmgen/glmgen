#include "tf.h"

void tf_dtx(double *x, int n, int k, double *a, double *b)
{
  int i;
  int j;

  memcpy(b, a, n*sizeof(double));
  
  if( k < 1 || k >= n )
    return;

  for(i=k; i > 0; --i)
  {

    /* b[0:n-i] = D' * b[0:n-i-1] for 1 <= i < n */
    b[n-i] = b[n-i-1];
    for(j=n-i-1; j > 0; --j)
    {
      b[j] = b[j-1] - b[j];
    }
    b[0] = -b[0];

    if( i != 1 )
    {
      for(j=0; j <= n-i; ++j)
      {
        b[j] = b[j] / ( x[j+i-1] - x[j] );
      }
    }
  }

  double fact = glmgen_factorial(k-1);
  for(i=0; i < n; ++i)
  {
    b[i] *= fact;
  }
}
