/****************************************************************************
 * Copyright (C) 2014 by Taylor Arnold, Ryan Tibshirani, Veerun Sadhanala   *
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
 * @file utils.c
 * @author Taylor Arnold, Ryan Tibshirani, Veerun Sadhanala
 * @date 2014-12-23
 * @brief Small generic utility functions for the glmgen library.
 */

#include "utils.h"
#include <math.h>

double glmgen_factorial(int n)
{
  int i=0;
  double x=1;

  for(i=2; i<=n; i++)
  {
    x *= i;
  }
  return x;
}

int is_nan(double x) {
  return (x != x);
}

int count_nans(double * x, int n) {
  int i;
  int num_nans = 0;

  for(i=0; i < n; i++) {
    num_nans += is_nan(x[i]);
  }
  return num_nans;
}

int has_nan(double * x, int n) {
  int i;
  for(i=0; i < n; i++) {
    if(is_nan(x[i]))
      return 1;
  }
  return 0;
}

double l1norm(double * x, int n) {
  int i;
  double s;
  s = 0;
  for(i=0; i < n; i++)
    s += fabs(x[i]);

  return s;
}

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

/* Sets A = w*A */
void diag_times_sparse (const cs * A, double * w)
{
  int i;
  int j;

  for (j = 0; j < A->n; j++)
  {
    for (i = A->p[j] ; i < A->p[j+1] ; i++)
    {
      A->x[i] *= w[ A->i[i] ];
    }
  }

}

/* Calculates A*b + D (identity), for a square matrix
   A and scalar b */
cs * scalar_plus_diag (const cs * A, double b, double *D)
{
  int i;
  int j;
  cs * B;

  B = cs_spalloc(A->m, A->n, A->nzmax, 1, 0);

  for (j = 0; j < A->n; j++)
  {
    B->p[j] = A->p[j];
    for (i = A->p[j] ; i < A->p[j+1] ; i++)
    {
      if(A->i[i] == j)
      {
        B->x[i] = b * A->x[i] + D[j];
      } else {
        B->x[i] = b * A->x[i];
      }
      B->i[i] = A->i[i];
    }
  }
  B->p[j] = A->p[j];

  return B;
}
