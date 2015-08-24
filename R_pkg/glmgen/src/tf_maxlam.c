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
 * @file tf_maxlam.c
 * @author Taylor Arnold, Veeranjaneyulu Sadhanala, Ryan Tibshirani
 * @date 2014-12-23
 * @brief Calculate the maximum lambda value for a trend filtering problem.
 * Will return the largest lambda value for which the penalty term is zero
 * for Gaussian losses. Will only be approximate for other loss functions.
 */

#include "tf.h"

/**
 * @brief Calculate maximum lambda value.
 * Will return the largest lambda value for which the penalty term is zero
 * for Gaussian losses. Will only be approximate for other loss functions.
 * @param n                    number of observations
 * @param y                    a vector of responses
 * @param Dt_qr                QR decomposition of the Dt matrix
 * @param w                    vector of sample weights; must be filled
 * @return  Returns the maximum value of lambda.
 */
double tf_maxlam (int n, double * y, gqr * Dt_qr, double * w)
{

  int i;
  int tmp; /*n-k-1*/
  double maxlam;
  double * y_work;

  y_work = (double *) malloc(n * sizeof(double));

  for(i = 0; i < n; i++) y_work[i] = sqrt(w[i]) * y[i];

  glmgen_qrsol(Dt_qr, y_work);

  tmp = Dt_qr->n;
  maxlam = 0;
  for(i = 0; i < tmp; i++) maxlam = MAX(maxlam, fabs(y_work[i]));

  free(y_work);
  return maxlam;
}
