#include "tf.h"

/* b(x) = log( 1 + exp(x) ).
 * Avoid computing exp(x) when x >> 0 */
double logi_b(double x)
{
  return x <= 0 ? log( 1 + exp(x) ) : x + log(1+exp(-x));
}

/* b1(x) = b'(x), the first derivative */
double logi_b1(double x)
{
  return x > 0 ? 1 / (1 + exp(-x) ) : exp(x)/ (1 + exp(x) );
}

/* b2(x) = b''(x), the second derivative */
double logi_b2(double x)
{
  x = -fabs(x);
  return exp(x-2*log(1+exp(x)));
}

double pois_b(double x)
{
  return exp(x);
}

double pois_b1(double x)
{
  return exp(x);
}

double pois_b2(double x)
{
  return exp(x);
}
