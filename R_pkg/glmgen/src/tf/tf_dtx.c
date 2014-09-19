#include <tf.h>

void tf_dtx(double *x, int n, int k, double *a, double *b)
{
  int i=0;
  int j=0;

  memcpy(b, a, n);

  for(i=k; i > 0; --i)
  {
    /* This was difft; inlining for efficency & to help
       the number of functions down */
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
  /* TODO mult(b, factorial(k-1)); */
}
