/**
  Optimizers for problems dealing with weighted TV-L1 norm regularization.

  @author Álvaro Barbero Jiménez
  @author Suvrit Sra
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <wchar.h>
#include <time.h>

#include "proxTV.h"

/*  tautString_TV1

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator

    min_x 0.5 ||x-y||^2 + lambda sum_i |x_i - x_(i-1)| .

    To do so a Taut String algorithm is used to solve its equivalent problem.

Inputs:
- y: reference signal.
- lambda: penalty parameter.
- x: array in which to store the solution.
- n: length of array y (and x).
 */
int tautString_TV1(double *y,double lambda,double *x,int n) {
  /* Minorant and minorant slopes */
  double mn, mx;
  /* Relative height of majorant and minorant slopes at the current points w.r.t. the tube center */
  double mnHeight, mxHeight;
  /* Last break points of minorant and majorant */
  int mnBreak, mxBreak;
  /* Last break point of taut string */
  int lastBreak;
  /* Auxiliary variables */
  int i, j;
  /* Helpful constants */
  const double minuslambda = -lambda;
  const double lambda2 = 2*lambda;
  const double minuslambda2 = 2*minuslambda;

#define CANCEL(txt,info) \
  printf("tautString_TV1: %s\n",txt); \
  return 0;

  /* Starting point */
  mnHeight = mxHeight = 0;
  mn = minuslambda + y[0];
  mx = lambda + y[0];
  lastBreak = -1;
  mnBreak = mxBreak = 0;

#ifdef DEBUG
  fprintf(DEBUG_FILE,"Starting taut-string with length=%d and penalty=%lf\n",n,lambda); fflush(DEBUG_FILE);
#endif

  /* Proceed along string */
  i = 0;
  while ( i < n ) {
    /* Loop over all points except the last one, that needs special care */
    while ( i < n-1 ) {
#ifdef DEBUG
      fprintf(DEBUG_FILE,"i = %d, mx = %.3f, mn = %.3f, mxHeight = %.3f, mnHeight = %.3f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,mx,mn,mxHeight,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
#endif

      /* Update height of minorant slope w.r.t. tube center */
      /* This takes into account both the slope of the minorant and the change in the tube center */
      mnHeight += mn - y[i];

      /* Check for ceiling violation: tube ceiling at current point is below proyection of minorant slope */
      /* Majorant is r + lambda (except for last point), which is computed on the fly */
      if ( lambda < mnHeight ) {
#ifdef DEBUG
        fprintf(DEBUG_FILE,"CEILING VIOLATION i = %d, mxVal = %f, mnHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,lambda,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
#endif
        /* Break segment at last minorant breaking point */
        i = mnBreak + 1;
        /* Build valid segment up to this point using the minorant slope */
        for ( j = lastBreak+1 ; j <= mnBreak ; j++ )
          x[j] = mn;
        /* Start new segment after the break */
        lastBreak = mnBreak;
        /* Build first point of new segment, which can be done in closed form */
        mn = y[i];
        mx = lambda2+y[i];
        mxHeight = lambda;
        mnHeight = minuslambda;
        mnBreak = mxBreak = i;
        i++;
        continue;
      }

      /* Update height of minorant slope w.r.t. tube center */
      /* This takes into account both the slope of the minorant and the change in the tube center */
      mxHeight += mx - y[i];

      /* Check for minorant violation: minorant at current point is above proyection of majorant slope */
      /* Minorant is r - lambda (except for last point), which is computed on the fly */
      if ( minuslambda > mxHeight ) {
#ifdef DEBUG
        fprintf(DEBUG_FILE,"FLOOR VIOLATION i = %d, mnVal = %f, mxHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,minuslambda,mxHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
#endif
        /* If violated, break segment at last majorant breaking point */
        i = mxBreak + 1;
        /* Build valid segment up to this point using the majorant slope */
        for ( j = lastBreak+1 ; j <= mxBreak ; j++ )
          x[j] = mx;
        /* Start new segment after the break*/
        lastBreak = mxBreak;
        /* Build first point of new segment, which can be done in closed form */
        mx = y[i];
        mn = minuslambda2+y[i];
        mxHeight = lambda;
        mnHeight = minuslambda;
        mnBreak = mxBreak = i;
        i++;
        continue;
      }

      /* No violations at this point */

      /* Check if proyected majorant height is above ceiling */
      if ( mxHeight >= lambda ) {
        /* Update majorant slope */
        mx += ( lambda - mxHeight ) / ( i - lastBreak );
        /* Get correct majorant height (we are touching it!) */
        mxHeight = lambda;
        /* This is a possible majorant breaking point */
        mxBreak = i;
      }

      /* Check if proyected minorant height is under actual minorant */
      if ( mnHeight <= minuslambda ) {
        /* Update minorant slope */
        mn += ( minuslambda - mnHeight ) / ( i - lastBreak );
        /* Compute correct minorant height (we are touching it!) */
        mnHeight = minuslambda;
        /* This is a possible minorant breaking point */
        mnBreak = i;
      }

#ifdef DEBUG
      fprintf(DEBUG_FILE,"i = %d, mx = %.3f, mn = %.3f, mxHeight = %.3f, mnHeight = %.3f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,mx,mn,mxHeight,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
#endif

      /* At this point: no violations, so keep up building current segment */
      i++;
    }

    /* Special case i == n-1 (last point) */
    /* We try to validate the last segment, and if we can, we are finished */
    /* The code is essentially the same as the one for the general case,
       the only different being that here the tube ceiling and floor are both 0 */

    /* Update height of minorant slope w.r.t. tube center */
    /* This takes into account both the slope of the minorant and the change in the tube center */
    mnHeight += mn - y[i];

    /* Check for ceiling violation: tube ceiling at current point is below proyection of minorant slope */
    /* Majorant is 0 at this point */
    if ( IS_POSITIVE(mnHeight) ) { // 0 < mnHeight
#ifdef DEBUG
      fprintf(DEBUG_FILE,"ENDING CEILING VIOLATION i = %d, mxVal = %f, mnHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,0.0,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
#endif
      /* Break segment at last minorant breaking point */
      i = mnBreak + 1;
      /* Build valid segment up to this point using the minorant slope */
      for ( j = lastBreak+1 ; j <= mnBreak ; j++ )
        x[j] = mn;
      /* Start new segment after the break */
      lastBreak = mnBreak;
      /* Go back to main loop, starting a new segment */
      /* We do not precompute the first point of the new segment here, as it might be n-1 and this leads to issues */
      mn = y[i];
      mx = lambda2+y[i];
      mxHeight = mnHeight = minuslambda;
      mnBreak = mxBreak = i;
      continue;
    }

    /* Update height of minorant slope w.r.t. tube center */
    /* This takes into account both the slope of the minorant and the change in the tube center */
    mxHeight += mx - y[i];

    /* Check for minorant violation: minorant at current point is above proyection of majorant slope */
    /* Minorant is 0 at this point */
    if ( IS_NEGATIVE(mxHeight) ) { // 0 > mxHeight
#ifdef DEBUG
      fprintf(DEBUG_FILE,"ENDING FLOOR VIOLATION i = %d, mnVal = %f, mxHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,0.0,mxHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
#endif
      /* If violated, break segment at last majorant breaking point */
      i = mxBreak + 1;
      /* Build valid segment up to this point using the majorant slope */
      for ( j = lastBreak+1 ; j <= mxBreak ; j++ )
        x[j] = mx;
      /* Start new segment after the break*/
      lastBreak = mxBreak;
      /* Go back to main loop, starting a new segment */
      /* We do not precompute the first point of the new segment here, as it might be n-1 and this leads to issues */
      mx = y[i];
      mn = minuslambda2+y[i];
      mxHeight = mnHeight = lambda;
      mnBreak = mxBreak = i;
      continue;
    }

    /* No violations at this point */

    /* Check if proyected minorant height is under actual minorant */
    if ( mnHeight <= 0 ) {
      /* Update minorant slope */
      mn += ( - mnHeight ) / ( i - lastBreak );
    }

#ifdef DEBUG
    fprintf(DEBUG_FILE,"i = %d, mx = %.3f, mn = %.3f, mxHeight = %.3f, mnHeight = %.3f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,mx,mn,mxHeight,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
#endif

    /* At this point: we are finished validating last segment! */
    i++;
  }

  /* Build last valid segment */
  for ( i = lastBreak+1 ; i < n ; i++ )
    x[i] = mn;

  /* Return */
  return 1;

#undef CANCEL
}


/*  tautString_TV1_Weighted

    Given a reference signal y and a penalty parameter lambda, solves the proximity operator

    min_x 0.5 ||x-y||^2 + sum_i \lambda_i |x_i - x_(i-1)| .

    To do so a Taut String algorithm is used to solve its equivalent problem.

Inputs:
- y: reference signal.
- lambda: penalties vector.
- x: array in which to store the solution.
- n: length of array y (and x).
 */
int tautString_TV1_Weighted(double *y,double *lambda,double *x,int n) {
  /* Minorant and minorant slopes */
  double mn, mx;
  /* Relative height of majorant and minorant slopes at the current points w.r.t. the tube center */
  double mnHeight, mxHeight;
  /* Last break points of minorant and majorant */
  int mnBreak, mxBreak;
  /* Last break point of taut string */
  int lastBreak;
  /* Auxiliary variables */
  int i, j;

#define CANCEL(txt,info) \
  printf("tautString_TV1_Weighted: %s\n",txt); \
  return 0;

#ifdef DEBUG
  fprintf(DEBUG_FILE, "tautString_TV1_Weighted start\n");
  fprintf(DEBUG_FILE, "y = ["); for(i=0;i<n&&i<DEBUG_N;i++) fprintf(DEBUG_FILE, " %.3f", y[i]); fprintf(DEBUG_FILE, "]\n");
  fprintf(DEBUG_FILE, "lambda = ["); for(i=0;i<n-1&&i<DEBUG_N;i++) fprintf(DEBUG_FILE, " %.3f", lambda[i]); fprintf(DEBUG_FILE, "]\n");
#endif

  /* Starting point */
  mnHeight = mxHeight = 0;
  mn = -lambda[0] + y[0];
  mx = lambda[0] + y[0];
  lastBreak = -1;
  mnBreak = mxBreak = 0;

  /* Proceed along string */
  i = 0;
  while ( i < n ) {
    /* Loop over all points except the last one, that needs special care */
    while ( i < n-1 ) {
#ifdef DEBUG
      fprintf(DEBUG_FILE,"i = %d, mx = %.3f, mn = %.3f, mxHeight = %.3f, mnHeight = %.3f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,mx,mn,mxHeight,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
#endif

      /* Update height of minorant slope w.r.t. tube center */
      /* This takes into account both the slope of the minorant and the change in the tube center */
      mnHeight += mn - y[i];

      /* Check for ceiling violation: tube ceiling at current point is below proyection of minorant slope */
      /* Majorant is r + lambda (except for last point), which is computed on the fly */
      if ( lambda[i] < mnHeight ) {
#ifdef DEBUG
        fprintf(DEBUG_FILE,"CEILING VIOLATION i = %d, mxVal = %f, mnHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,lambda[i],mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
#endif
        /* Break segment at last minorant breaking point */
        i = mnBreak + 1;
        /* Build valid segment up to this point using the minorant slope */
        for ( j = lastBreak+1 ; j <= mnBreak ; j++ )
          x[j] = mn;
        /* Start new segment after the break */
        lastBreak = mnBreak;
        /* Build first point of new segment, which can be done in closed form */
        mn = y[i] + lambda[i-1] - lambda[i];  // When computing the slopes we need to account for the differences in lambdas
        mx = y[i] + lambda[i-1] + lambda[i]; // Here too
        mxHeight = lambda[i];
        mnHeight = -lambda[i];
        mnBreak = mxBreak = i;
        i++;
        continue;
      }

      /* Update height of minorant slope w.r.t. tube center */
      /* This takes into account both the slope of the minorant and the change in the tube center */
      mxHeight += mx - y[i];

      /* Check for minorant violation: minorant at current point is above proyection of majorant slope */
      /* Minorant is r - lambda (except for last point), which is computed on the fly */
      if ( -lambda[i] > mxHeight ) {
#ifdef DEBUG
        fprintf(DEBUG_FILE,"FLOOR VIOLATION i = %d, mnVal = %f, mxHeight = %f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,-lambda[i],mxHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
#endif
        /* If violated, break segment at last majorant breaking point */
        i = mxBreak + 1;
        /* Build valid segment up to this point using the majorant slope */
        for ( j = lastBreak+1 ; j <= mxBreak ; j++ )
          x[j] = mx;
        /* Start new segment after the break*/
        lastBreak = mxBreak;
        /* Build first point of new segment, which can be done in closed form */
        mx = y[i] - lambda[i-1] + lambda[i]; // When computing the slopes we need to account for the differences in lambdas
        mn = y[i] - lambda[i-1] - lambda[i]; // Here too
        mxHeight = lambda[i];
        mnHeight = -lambda[i];
        mnBreak = mxBreak = i;
        i++;
        continue;
      }

      /* No violations at this point */

      /* Check if proyected majorant height is above ceiling */
      if ( mxHeight >= lambda[i] ) {
        /* Update majorant slope */
        mx += ( lambda[i] - mxHeight ) / ( i - lastBreak );
        /* Get correct majorant height (we are touching it!) */
        mxHeight = lambda[i];
        /* This is a possible majorant breaking point */
        mxBreak = i;
      }

      /* Check if proyected minorant height is under actual minorant */
      if ( mnHeight <= -lambda[i] ) {
        /* Update minorant slope */
        mn += ( -lambda[i] - mnHeight ) / ( i - lastBreak );
        /* Compute correct minorant height (we are touching it!) */
        mnHeight = -lambda[i];
        /* This is a possible minorant breaking point */
        mnBreak = i;
      }

#ifdef DEBUG
      fprintf(DEBUG_FILE,"i = %d, mx = %.3f, mn = %.3f, mxHeight = %.3f, mnHeight = %.3f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,mx,mn,mxHeight,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
#endif

      /* At this point: no violations, so keep up building current segment */
      i++;
    }

    /* Special case i == n-1 (last point) */
    /* We try to validate the last segment, and if we can, we are finished */
    /* The code is essentially the same as the one for the general case,
       the only different being that here the tube ceiling and floor are both 0 */

    /* Update height of minorant slope w.r.t. tube center */
    /* This takes into account both the slope of the minorant and the change in the tube center */
    mnHeight += mn - y[i];

    /* Check for ceiling violation: tube ceiling at current point is below proyection of minorant slope */
    /* Majorant is 0 at this point */
    if ( IS_POSITIVE(mnHeight) ) { // 0 < mnHeight
#ifdef DEBUG
      fprintf(DEBUG_FILE,"ENDING CEILING VIOLATION i = %d, mxVal = %f, mnHeight = %g, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,0.0,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
#endif
      /* Break segment at last minorant breaking point */
      i = mnBreak + 1;
      /* Build valid segment up to this point using the minorant slope */
      for ( j = lastBreak+1 ; j <= mnBreak ; j++ )
        x[j] = mn;
      /* Start new segment after the break */
      lastBreak = mnBreak;
      /* Go back to main loop, starting a new segment */
      /* We do not precompute the first point of the new segment here, as it might be n-1 and this leads to issues */
      mn = y[i] + lambda[i-1] - (i==n-1 ? 0 : lambda[i]);  // When computing the slopes we need to account for the differences in lambdas
      mx = y[i] + lambda[i-1] + (i==n-1 ? 0 : lambda[i]); // Here too
      mxHeight = mnHeight = -lambda[i-1];
      mnBreak = mxBreak = i;
      continue;
    }

    /* Update height of minorant slope w.r.t. tube center */
    /* This takes into account both the slope of the minorant and the change in the tube center */
    mxHeight += mx - y[i];

    /* Check for minorant violation: minorant at current point is above proyection of majorant slope */
    /* Minorant is 0 at this point */
    if ( IS_NEGATIVE(mxHeight) ) { // 0 > mxHeight
#ifdef DEBUG
      fprintf(DEBUG_FILE,"ENDING FLOOR VIOLATION i = %d, mnVal = %f, mxHeight = %g, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,0.0,mxHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
#endif
      /* If violated, break segment at last majorant breaking point */
      i = mxBreak + 1;
      /* Build valid segment up to this point using the majorant slope */
      for ( j = lastBreak+1 ; j <= mxBreak ; j++ )
        x[j] = mx;
      /* Start new segment after the break*/
      lastBreak = mxBreak;
      /* Go back to main loop, starting a new segment */
      /* We do not precompute the first point of the new segment here, as it might be n-1 and this leads to issues */
      mx = y[i] - lambda[i-1] + (i==n-1 ? 0 : lambda[i]); // When computing the slopes we need to account for the differences in lambdas
      mn = y[i] - lambda[i-1] - (i==n-1 ? 0 : lambda[i]); // Here too
      mxHeight = mnHeight = lambda[i-1];
      mnBreak = mxBreak = i;
      continue;
    }

    /* No violations at this point */

    /* Check if proyected minorant height is under actual minorant */
    if ( mnHeight <= 0 ) {
      /* Update minorant slope */
      mn += ( - mnHeight ) / ( i - lastBreak );
    }

#ifdef DEBUG
    fprintf(DEBUG_FILE,"i = %d, mx = %.3f, mn = %.3f, mxHeight = %.3f, mnHeight = %.3f, mxBreak = %d, mnBreak = %d, lastBreak = %d\n",i,mx,mn,mxHeight,mnHeight,mxBreak,mnBreak,lastBreak); fflush(DEBUG_FILE);
#endif

    /* At this point: we are finished validating last segment! */
    i++;
  }

  /* Build last valid segment */
  for ( i = lastBreak+1 ; i < n ; i++ )
    x[i] = mn;

  /* Return */
  return 1;

#undef CANCEL
}
