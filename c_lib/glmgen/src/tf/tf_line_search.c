#include "tf.h"
#include <math.h>

double l1_norm(double * x, int n);
void sum(double * x, double * y, double * z, int n);

/* backtracking linesearch. See Jason Lee et. al (2012)
 * f = g + h, g differentiable, h not necessarily so
 * Find t s.t t = gamma^s, x+ = x + td, 
 * f(x+) - f(x) \leq alpha t theta 
 * where theta = <\nabla g, d> + h(x+d) - h(x)
 * g(x) = -<y,x> + b(x)
 * h(x) = lam |Dx|*/
double tf_line_search(double * y, double * x, int n, int k, 
    double lam, 
    func_RtoR b, func_RtoR b1, 
    double * beta, double * d, 
    double alpha, double gamma, int max_iter,
    int * iter)
{
  int i, it;
  double norm_Db, norm_Dbn;
  double ip_yd;
  double theta;
  double t;
  double descent;

  double * Db, * Dd, *Dbn; /* Dbn : Dbnew */
  
  Db = (double*)malloc(n*sizeof(double));
  Dd = (double*)malloc(n*sizeof(double));
  Dbn = (double*)malloc(n*sizeof(double));

  tf_dx(x, n, k+1, beta, Db);
  tf_dx(x, n, k+1, d, Dd);

  /* Compute theta */
  theta = 0;
  for(i = 0; i < n; i++)
  {
    theta += (-y[i] + b1(beta[i])) * d[i];    
  }
  norm_Db = l1_norm(Db, n);
  sum(Dd, Db, Dbn, n);
  theta += lam * ( l1_norm(Dbn, n) - norm_Db );

  ip_yd = 0;
  for(i =0; i < n; i++)
  {
    ip_yd += y[i] * d[i];
  }

  t = 1;
  for(it = 0; it < max_iter; it++)
  {
    /* Compute descent: f(beta+) - f(beta) */
    descent = -t * ip_yd;
    norm_Dbn = 0;
    for(i = 0; i < n; i++)
    {
      descent += b(beta[i] + t * d[i]) - b(beta[i]);
      norm_Dbn += fabs(Db[i] + t * Dd[i]);
    }
    descent += lam * (norm_Dbn - norm_Db);
    
    /* Check if the descent is sufficient */
    if(descent <= alpha * t * theta)
    {
      *iter = it;
      return t;
    }
    else
      t = t * gamma;
  }

  *iter = it;
  return t;
}


double l1_norm(double * x, int n)
{
  int i;
  double s = 0;
  for(i = 0; i < n; i++)
  {
    s += fabs( x[i] );
  }
  return s;

}

/* z = x + y */
void sum(double * x, double * y, double * z, int n)
{
  int i;
  for(i = 0; i < n; i++)
  {
    z[i] = x[i] + y[i];
  }
}
