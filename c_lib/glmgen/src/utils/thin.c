#include "utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void thin( double* x, double* y, double* w, int n, int k, 
  double** xt, double** yt, double** wt, int* nt_ptr, double x_cond)
{
  int i,j, jj;
  int m;  // number of intervals
  int nt; // number of intervals with at least one point
  double r, delta;
  double mindx;
  int * intvl;
  int intvl_xj;
  int lo, hi;
  int cur_intvl;

  r = x[n-1] - x[0];
  delta = r * pow( n*x_cond, -1./(k+1) );
 
  mindx = r;
  for(i = 0; i < n-1; i++) 
    mindx = MIN(x[i+1] - x[i], mindx);

  *xt = *yt = *wt = NULL;
  
  if( mindx >= delta ) return;

  m = (int) MIN( floor(r/delta), 5*n );
  delta = r/m;

  if( m <= 1 ) return; // Not thinning as it merges all points into one

  intvl = (int*)malloc( n * sizeof(int) );

  nt = 0;
  for(j = 0; j < n; j++)
  {
    intvl_xj = (int) floor( (x[j] -x[0]) / delta ) + 1;
    intvl[j] = MAX(1, MIN(intvl_xj, m));

    if( j == 0 || intvl[j] != intvl[j-1] ) nt++;
  }

  *nt_ptr = nt;

  *xt = (double*)malloc( nt * sizeof(double) );
  *yt = (double*)malloc( nt * sizeof(double) );
  *wt = (double*)malloc( nt * sizeof(double) );

  lo = 0;
  hi = 0;
  i = 0; // range 0:nt-1
  cur_intvl = 1; // range 1:m

  for(j = 0; j < n; j++)
  {
    if( intvl[j] > cur_intvl ) //crossed the current interval
    {
      hi = j-1;
      (*xt)[i] = x[0] + (cur_intvl - 0.5) * delta;
      
      (*wt)[i] = (*yt)[i] = 0.;
      for( jj = lo; jj <= hi; jj++) 
      {
        (*wt)[i] += w[jj];
        (*yt)[i] += w[jj] * y[jj];
      }
      (*yt)[i] /= (*wt)[i];

      i++;
      cur_intvl = intvl[j];
      lo = j;
    }
    if( i >= nt - 1 ) // last interval
    {
      i = nt - 1;
      hi = n-1;
      (*xt)[i] = x[0] + (m - 0.5) * delta;
      
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
