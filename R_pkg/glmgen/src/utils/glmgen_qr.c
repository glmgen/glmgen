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
 * @file glmgen_qr.c
 * @author Taylor Arnold, Veeranjaneyulu Sadhanala, Ryan Tibshirani
 * @date 2014-12-23
 * @brief Construct a single qr decomposition object.
 */

#include "utils.h"
#include "cs.h"

/* Calculates symbolic and numeric qr for a given matrix A */
gqr * glmgen_qr (const cs * A)
{
  gqr * B = malloc(1 * sizeof(cs));
  B->m = A->m;
  B->n = A->n;
  B->S = cs_sqr (1, A, 1) ;
  B->N = cs_qr (A, (B->S)) ;
  B->W = cs_calloc ((B->S) ? ((B->S))->m2 : 1, sizeof (double)) ;
  return (B);
}

/* Free a constructed gqr struct */
csi glmgen_gqr_free (gqr * A)
{
  cs_sfree(A->S);
  cs_nfree(A->N);
  free(A->W);
  free(A);
  return(0);
}

/* Takes gqr structure from matrix A and solves Ax=b, overwriting the input b */
csi glmgen_qrsol (gqr * A, double * b)
{
  csi k;

  if(A->m < A->n) return(1); /* we only deal with m >= n case here */
  cs_ipvec (A->S->pinv, b, A->W, A->m) ;   /* x(0:m-1) = b(p(0:m-1) */
  for (k = 0 ; k < A->n ; k++)       /* apply Householder refl. to x */
  {
    cs_happly (A->N->L, k, A->N->B [k], A->W) ;
  }
  cs_usolve (A->N->U, A->W) ;           /* x = R\x */
  cs_ipvec (A->S->q, A->W, b, A->n) ;      /* b(q(0:n-1)) = x(0:n-1) */

  return (1) ;
}
