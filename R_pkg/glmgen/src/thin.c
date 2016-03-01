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
 * @file thin.c
 * @author Taylor Arnold, Veeranjaneyulu Sadhanala, Ryan Tibshirani
 * @date 2014-12-23
 * @brief Code for thinning large or ill-conditioned data sets.
 */

#include "utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void thin_old(double* x, double* y, double* w, int n, int k,
	      double** xt, double** yt, double** wt, int* nt_ptr, double x_cond)
{
  double r;
  double delta;

  r = x[n-1] - x[0];
  delta = r * pow( n*x_cond, -1./(k+1) );

  thin(x,y,w,n,k,xt,yt,wt,nt_ptr,delta);
}

void thin(double* x, double* y, double* w, int n, int k,
	  double** xt, double** yt, double** wt, int* nt_ptr, double tol)
{
  int i,j, jj;
  int m;  /* number of intervals */
  int nt; /* number of intervals with at least one point */
  double r;
  double mindx;
  int * intvl;
  int intvl_xj;
  int lo;
  int hi;
  int cur_intvl;

  r = x[n-1] - x[0];

  // Do not thin if minimum separation of x (mindx) is > tol
  mindx = r;
  for(i = 0; i < n-1; i++)
    mindx = MIN(x[i+1] - x[i], mindx);

  if( mindx > tol ) return;

  *xt = *yt = *wt = NULL;	
  m = (int) MAX(1, floor(r/tol));	
  tol = r/m;
  intvl = (int*)malloc( n * sizeof(int) );

  nt = 0;
  for(j = 0; j < n; j++)
    {
      intvl_xj = (int) floor( (x[j]-x[0]) / tol ) + 1;
      intvl[j] = MAX(1, MIN(intvl_xj, m));
      if( j == 0 || intvl[j] != intvl[j-1] ) nt++;
    }

  *nt_ptr = nt;

  *xt = (double*)malloc( nt * sizeof(double) );
  *yt = (double*)malloc( nt * sizeof(double) );
  *wt = (double*)malloc( nt * sizeof(double) );

  lo = 0;
  hi = 0;
  i = 0; /* range 0:nt-1 */
  cur_intvl = 1; /* range 1:m */

  for(j = 0; j < n; j++)
    {
      if( intvl[j] > cur_intvl ) /* crossed the current interval */
	{
	  hi = j-1;
	  (*xt)[i] = x[0] + (cur_intvl - 0.5) * tol;

	  (*wt)[i] = (*yt)[i] = 0.;
	  for(jj = lo; jj <= hi; jj++)
	    {
	      (*wt)[i] += w[jj];
	      (*yt)[i] += w[jj] * y[jj];
	    }
	  (*yt)[i] /= (*wt)[i];

	  i++;
	  cur_intvl = intvl[j];
	  lo = j;
	}
      if( i >= nt - 1 ) /* last interval */
	{
	  i = nt - 1;
	  hi = n-1;
	  (*xt)[i] = x[0] + (m - 0.5) * tol;

	  (*wt)[i] = (*yt)[i] = 0.;
	  for( jj = lo; jj <= hi; jj++)
	    {
	      (*wt)[i] += w[jj];
	      (*yt)[i] += w[jj] * y[jj];
	    }
	  (*yt)[i] /= (*wt)[i];

	  break;
	}
    }

  free(intvl);
}
