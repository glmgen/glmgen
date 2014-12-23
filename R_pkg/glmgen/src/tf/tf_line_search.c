#include "tf.h"
#include "math.h"

/* backtracking linesearch. See Jason Lee et. al (2012)
 * f = g + h, g differentiable, h not necessarily so
 * Find t s.t t = gamma^s, x+ = x + td,
 * f(x+) - f(x) \leq alpha t theta
 * where theta = <\nabla g, d> + h(x+d) - h(x)
 * g(x) = W * (-<y,x> + b(x))
 * h(x) = lam |Dx|
 */
double tf_line_search(double * y, double * x, double * w, int n, int k, double lam,
    func_RtoR b, func_RtoR b1,
    double * beta, double * d,
    double alpha, double gamma, int max_iter,
    int * iter, double * Db, double * Dd)
{
  int i, it;
  double norm_Db, norm_Dbn;
  double grad_term, pen_term;
  double pen_diff;
  double ip_yd;
  double theta;
  double t;
  double descent;

  tf_dx(x, n, k+1, beta, Db);
  tf_dx(x, n, k+1, d, Dd);

  /* Compute theta */
  theta = 0;
  norm_Db = 0;
  norm_Dbn = 0;
  grad_term = 0;
  t = 1;
  for(i = 0; i < n; i++)
  {
    theta += w[i] * (-y[i] + b1(beta[i])) * d[i];
    norm_Db += fabs(Db[i]);
    norm_Dbn += fabs(Db[i] + t * Dd[i]);
    grad_term += w[i] * (-y[i] + b1(beta[i])) * d[i];
  }
  theta += lam * ( norm_Dbn - norm_Db );

  ip_yd = 0;
  for(i = 0; i < n; i++)
  {
    ip_yd += w[i] * y[i] * d[i];
  }

  t = 1;
  for(it = 0; it < max_iter; it++)
  {
    /* Compute descent: f(beta+) - f(beta) */
    descent = -t * ip_yd;
    norm_Dbn = 0;
    for(i = 0; i < n; i++)
    {
      descent += w[i] * (b(beta[i] + t * d[i]) - b(beta[i]));
      norm_Dbn += fabs(Db[i] + t * Dd[i]);
    }
    descent += lam * (norm_Dbn - norm_Db);
    pen_term = lam * (norm_Dbn - norm_Db);

    double bound = alpha * t * theta;

    /* Check if the descent is sufficient */
    if (descent <= bound) break;
    else t = t * gamma;
  }

  *iter = it;
  return t;
}

