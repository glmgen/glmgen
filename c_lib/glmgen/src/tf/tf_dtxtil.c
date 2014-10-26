#include "tf.h"

/* b = (\tilde{D}^{(x,k)})^T a = k ( D^{(x,k)}^T \Delta_k)^{-1} a
*/
void tf_dtxtil(double *x, int n, int k,double *a, double *b)
{
  int i;

  if( k > 0 )
    for(i=0; i < n-k; i++)
    {
      a[i] = a[i] * k/( x[k+i] - x[i] );
    }
  tf_dtx(x, n, k, a, b);

}
