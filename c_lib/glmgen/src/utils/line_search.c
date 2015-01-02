/****************************************************************************
 * Copyright (C) 2014 by Taylor Arnold, Veeranjaneyulu Sadhanala,           *
 *                       Ryan Tibshirani                                    *
 *                                                                          *
 * This file is part of the glmgen library / package.                       *
 *                                                                          *
 *   glmgen is free software: you can redistribute it and/or modify it      *
 *   under the terms of the GNU Lesser General Public License as published  *
 *   by the Free Software Foundation, either version 2 of the License, or   *
 *   (at your option) any later version.                                    *
 *                                                                          *
 *   glmgen is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU Lesser General Public License for more details.                    *
 *                                                                          *
 *   You should have received a copy of the GNU Lesser General Public       *
 *   License along with glmgen. If not, see <http://www.gnu.org/licenses/>. *
 ****************************************************************************/

/**
 * @file line_search.c
 * @author Taylor Arnold, Veeranjaneyulu Sadhanala, Ryan Tibshirani
 * @date 2014-12-23
 * @brief Generic line search algorithm.
 *
 * Here.
 */

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
double line_search(double * y, double * x, double * w, int n, int k, double lam,
    func_RtoR b, func_RtoR b1,
    double * beta, double * d,
    double alpha, double gamma, int max_iter,
    int * iter, double * Db, double * Dd)
{
  int i;
  int it;
  double norm_Db;
  double norm_Dbn;
  double grad_term;
  double ip_yd;
  double theta;
  double t;
  double descent;
  double bound;

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
    bound = alpha * t * theta;

    /* Check if the descent is sufficient */
    if (descent <= bound) break;
    else t = t * gamma;
  }

  *iter = it;
  return t;
}

