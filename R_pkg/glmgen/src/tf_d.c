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
 * @file tf_d.c
 * @author Taylor Arnold, Veeranjaneyulu Sadhanala, Ryan Tibshirani
 * @date 2014-12-24
 * @brief Utility functions for working with the trend filtering penalty matrix.
 * The penalty matrix D in the trend filtering problem has many nice properties
 * which distinquish it from arbitrary sparse of banded matricies. Much of the
 * speed and stability of the trend filtering algorithms come from taking advantage
 * of these properties. Helper functions for doing so are collected here.
 */

#include "tf.h"

/**
 * @brief Creates the penalty matrix of order k.
 * Returns the matrix Dk as a suite sparse style matrix.
 *
 * @param n                    number of observations
 * @param k                    order of the trendfilter
 * @param x                    locations of the responses
 * @return pointer to a csparse matrix
 * @see tf_calc_dktil
 */
cs * tf_calc_dk (int n, int k, const double * x)
{
  long int i;

  int tk = 1; /* "this k" - will iterate until ts = k */

  cs * D1;
  cs * D1_cp;
  cs * Dk;
  cs * Dk_cp;
  cs * delta_k;
  cs * delta_k_cp;
  cs * D1_x_delta;
  cs * Dk_next;
  cs * T;
  cs * eye;

  /* Deal with k=0 separately */
  if(k == 0)
  {
    T = cs_spalloc (n, n, n, 1, 1) ;
    for (i = 0 ; i < n; i++) cs_entry (T, i, i, 1);
    eye = cs_compress (T);
    cs_spfree (T);
    return eye;
  }

  /* Contruct one 'full D1', which persists throughout
     and another copy as Dk */
  D1 = cs_spalloc(n-tk, n, (n-tk)*2, 1, 1);
  Dk = cs_spalloc(n-tk, n, (n-tk)*2, 1, 1);
  D1->nz = (n-tk)*2;
  Dk->nz = (n-tk)*2;
  for (i = 0; i < (n-tk)*2; i++)
  {
    D1->p[i] = (i+1) / 2;
    Dk->p[i] = D1->p[i];
    D1->i[i] = i / 2;
    Dk->i[i] = D1->i[i];
    D1->x[i] = -1 + 2*(i % 2);
    Dk->x[i] = D1->x[i];
  }

  /* Create a column compressed version of Dk, and
     delete the old copy */
  Dk_cp = cs_compress(Dk);
  cs_spfree(Dk);

  for (tk = 1; tk < k; tk++)
  {
    /* 'reduce' the virtual size of D1 to: (n-tk-1) x (n-tk),
       compress into compressed column, saving as D1_cp */
    D1->nz = (n-tk-1)*2;
    D1->m = n-tk-1;
    D1->n = n-tk;
    D1_cp = cs_compress(D1);

    /* Construct diagonal matrix of differences: */
    delta_k = cs_spalloc(n-tk, n-tk, (n-tk), 1, 1);
    for(i = 0; i < n - tk; i++)
    {
      delta_k->p[i] = i;
      delta_k->i[i] = i;
      delta_k->x[i] = tk / (x[tk + i] - x[i]);
    }
    delta_k->nz = n-tk;
    delta_k_cp = cs_compress(delta_k);
    D1_x_delta = cs_multiply(D1_cp, delta_k_cp);

    /* Execute the matrix multiplication */
    Dk_next = cs_multiply(D1_x_delta, Dk_cp);

    /* Free temporary cs matricies created in each loop */
    cs_spfree(D1_cp);
    cs_spfree(delta_k);
    cs_spfree(delta_k_cp);
    cs_spfree(D1_x_delta);
    cs_spfree(Dk_cp);
    Dk_cp = Dk_next;
  }

  cs_spfree(D1);
  return Dk_cp;
}

/**
 * @brief Creates the penalty matrix D tilde of order k.
 * Returns the matrix Dk premultipied by a diagonal
 * matrix of weights.
 *
 * @param n                    number of observations
 * @param k                    order of the trendfilter
 * @param x                    locations of the responses
 * @return pointer to a csparse matrix
 * @see tf_calc_dktil
 */
cs * tf_calc_dktil (int n, int k, const double * x)
{
  cs * delta_k;
  cs * delta_k_cp;
  cs * Dk;
  cs * Dktil;

  int i;

  Dk = tf_calc_dk(n, k, x);

  /* Deal with k=0 separately */
  if(k == 0)
    return Dk;

  /* Construct diagonal matrix of differences: */
  delta_k = cs_spalloc(n-k, n-k, (n-k), 1, 1);
  for(i = 0; i < n - k; i++)
  {
    delta_k->p[i] = i;
    delta_k->i[i] = i;
    delta_k->x[i] = k / (x[k + i] - x[i]);
  }
  delta_k->nz = n-k;
  delta_k_cp = cs_compress(delta_k);
  Dktil = cs_multiply(delta_k_cp, Dk);

  cs_spfree(Dk);
  cs_spfree(delta_k);
  cs_spfree(delta_k_cp);

  return Dktil;
}

/**
 * @brief Multiplies a vector by D, without having to explictly
 * construct or use the matrix D. In symbols, Da = b.
 *
 * @param x                    locations of the responses
 * @param n                    number of observations
 * @param k                    order of the trendfilter
 * @param a                    the input vector to multiply
 * @param b                    allocated space for the output
 * @return void
 * @see tf_dxtil
 */
void tf_dx(double *x, int n, int k,double *a, double *b)
{
  int i;
  int j;
  double fact;

  for(i=0; i < n; i++) b[i] = a[i];

  if( k < 1 || k >= n )
    return;

  for(i=0; i < k; ++i)
  {
    if( i != 0 )
    {
      /* b[i:n-1] = b[i:n-1] ./ ( x[i:n-1] - x[0:n-1-i] ) */
      for(j=i; j < n; ++j)
      {
        b[j] = b[j] / ( x[j] - x[j-i]);
      }
    }

    /* b[i+1:n-1] = -b[i:n-2] + b[i+1:n-1] */
    for(j=n-1; j >= i+1; --j)
    {
      b[j] = b[j] - b[j-1];
    }
  }

  fact = glmgen_factorial(k-1);
  for(i=0; i < n; ++i)
  {
    b[i] *= fact;
  }

  /* move the solution to the beginning of the array */
  memmove(b, b+k, (n-k)*sizeof(double));
}

/**
 * @brief Multiplies a vector by D tilde, without having to
 * explictly construct or use the matrix D.
 *
 * @param x                    locations of the responses
 * @param n                    number of observations
 * @param k                    order of the trendfilter
 * @param a                    the input vector to multiply
 * @param b                    allocated space for the output
 * @return void
 * @see tf_dx
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

/**
 * @brief Multiplies a vector by D transpose, without having
 * to explictly construct or use the matrix D.
 *
 * @param x                    locations of the responses
 * @param n                    number of observations
 * @param k                    order of the trendfilter
 * @param a                    the input vector to multiply
 * @param b                    allocated space for the output
 * @return void
 * @see tf_dtxtil
 */
void tf_dtx(double *x, int n, int k, double *a, double *b)
{
  int i;
  int j;
  double fact;

  for(i=0; i < n-k; i++) b[i] = a[i];

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

  fact = glmgen_factorial(k-1);
  for(i=0; i < n; ++i)
  {
    b[i] *= fact;
  }
}

/**
 * @brief Multiplies a vector by D tilde transpose, without having
 * to explictly construct or use the matrix D.
 *
 * @param x                    locations of the responses
 * @param n                    number of observations
 * @param k                    order of the trendfilter
 * @param a                    the input vector to multiply
 * @param b                    allocated space for the output
 * @return void
 * @see tf_dtx
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
