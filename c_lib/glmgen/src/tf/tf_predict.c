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
 * @file tf_predict.c
 * @author Taylor Arnold, Ryan Tibshirani, Veerun Sadhanala
 * @date 2014-12-23
 * @brief Prediction algorithms for trendfiltering.
 *
 * Routines take a beta vector, and therefore are agnostic to
 * the algorithm used for calculating the fit.
 */

 #include "tf.h"

/**
 * @brief Calculate maximum lambda value.
 * Will return the largest lambda value for which the penalty term is zero
 * for Gaussian losses. Will only be approximate for other loss functions.
 * @param beta                 the beta vector for the prediction; length n-k
 * @param x                    the original positions used in the fit
 * @param n                    number of observations
 * @param k                    order of the fit
 * @param family               family of the fit
 * @param x0                   the new positions to predict at
 * @param n0                   the number of observations in x0
 * @param pred                 allocated space for the predicted values
 * @param zero_tol             tolerance for the fitting algorithm; default is 1e-11
 * @return  void
 * @note The results will not be valid unless the values in x0 are within the range of the original x inputs.
 */
void tf_predict(double * beta, double * x, int n, int k, int family,
                double * x0, int n0, double * pred, double zero_tol)
{
  int i;
  double f;

  tf_predict_gauss(beta, x, n, k, x0, n0, pred, zero_tol);

  switch (family)
  {
    case FAMILY_GAUSSIAN:
      break;

    case FAMILY_LOGISTIC:
      for (i = 0; i < n0; i++)
      {
        f = logi_b1(pred[i]);
        pred[i] = f;
      }
      break;

    case FAMILY_POISSON:
      for (i = 0; i < n0; i++)
      {
        f = pois_b1(pred[i]);
        pred[i] = f;
      }
      break;
  }
}

/**
 * @brief Calculate maximum lambda value.
 * Lower level function for predicting from a Gaussian loss function.
 * @param beta                 the beta vector for the prediction; length n-k
 * @param x                    the original positions used in the fit
 * @param n                    number of observations
 * @param k                    order of the fit
 * @param x0                   the new positions to predict at
 * @param n0                   the number of observations in x0
 * @param pred                 allocated space for the predicted values
 * @param zero_tol             tolerance for the fitting algorithm; default is 1e-11
 * @return  void
 * @see tf_predict
 */
void tf_predict_gauss(double * beta, double * x, int n, int k,
		      double * x0, int n0, double * pred, double zero_tol)
{
  int i;
  int j;
  int l;
  double * phi;
  double * theta;
  double k_fac;
  double h;

  if(n0 <= 0) return;

  /* Compute phi (polynomial coefficients) */
  phi = (double *)malloc((k+1)*sizeof(double));
  poly_coefs(x,k,beta,phi);

  /* Compute theta (falling fact coefficients) */
  theta = (double *)malloc((n)*sizeof(double));
  tf_dx(x,n,k+1,beta,theta);
  k_fac = glmgen_factorial(k);
  for(i=0; i<n-k-1;i++)
    theta[i] /= k_fac;

  /* Threshold small values */
  for (i=0; i<n-k-1; i++) if (fabs(theta[i])<zero_tol) theta[i]=0;

  /* Compute the predictions at each new point x0 */
  for (j=0; j<n0; j++) {
    pred[j] = 0;

    /* Loop over x points, polynomial basis */
    for (i=0; i<k+1; i++) {
      h = 1;
      l=0;
      for (l=0; l<i; l++) {
        h *= (x0[j]-x[l]);
      }
      pred[j] += phi[i]*h;
    }

    /* Loop over x points, falling fact basis */
    for (i=0; i<n-k-1; i++) {
      /* If the current x0 is too small, then break */
      if (x0[j]<=x[i+k]) break;

      /* Otherwise check the ith coef, and if it is nonzero,
       * compute the contribution of the ith basis function */
      if (theta[i]!=0) {
	      h = 1;
	      for (l=0; l<k; l++) {
	        h *= (x0[j]-x[i+l+1]);
	      }
	      pred[j] += theta[i]*h;
      }
    }
  }

  free(phi);
  free(theta);
}

/**
 * @brief Calculate polynomial coefficents of the fit.
 * Helper function to convert the beta vector into a polynomial.
 * @param x                    the original positions used in the fit
 * @param k                    order of the fit
 * @param beta                 the beta vector for the prediction; length n-k
 * @param phi                  allocated memory of length k+1
 * @return  void
 */
void poly_coefs(double *x, int k, double *beta, double *phi)
{
  int j;
  int ell;

  memcpy(phi,beta,(k+1)*sizeof(double));

  for(j=1; j <= k; j++)
  {
    for(ell = k; ell >= j; ell--)
    {
      phi[ell] = (phi[ell] - phi[ell-1]) / ( x[ell] - x[ell-j] );
    }
  }
}

