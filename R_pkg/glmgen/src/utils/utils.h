#ifndef UTILS_H
#define UTILS_H

#include "cs.h"

/* --- custom structure to hold the out of the qr decomposition ------- */
typedef struct gcs_qr  /* holds numeric and symbol qr together */
{
  csi m ;         /* number of rows */
  csi n ;         /* number of columns */
  css *S ;        /* symbolic qr output */
  csn *N ;        /* numeric qr output */
  double * W ;    /* pre-allocated working space for the solver */
} gqr ;

/* Utility functions for solving a linear system with a
   pre-computed sparse qr object (struct gcs_qr / gqr) */
gqr * glmgen_qr (const cs * A);
csi glmgen_qrsol (gqr * B, double * b);
csi glmgen_gqr_free (gqr * A);

void soft_thresh(int n, double *y, double lam, double *beta);
cs * scalar_plus_eye (const cs * A, double b);

#endif
