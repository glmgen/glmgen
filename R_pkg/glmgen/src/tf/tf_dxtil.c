#include "tf.h"

/* b = \tilde{D}^{(x,k)} a = k (\Delta_k)^{-1} D^{(x,k)} a
*/
void tf_dxtil(double *x, int n, int k,double *a, double *b)
{
  int i;

  tf_dx(x, n, k, a, b);

  if( k > 0 )
    for(i=0; i < n-k; i++)
    {
      b[i] = b[i] * k/( x[k+i] - x[i] );
    }
}
