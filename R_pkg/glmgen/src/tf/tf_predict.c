#include "tf.h"

void tf_predict_gauss(double * beta, double * x, int n, int k,
		      double * x0, int n0, double * pred, double zero_tol)
{
  if(n0 <= 0) return;
  int i=0, j=0;
  /* Compute phi (polynomial coefficients) */
  double *phi = (double *)malloc((k+1)*sizeof(double));
  poly_coefs(x,n,k,beta,phi);

  /* Compute theta (falling fact coefficients) */
  double *theta = (double *)malloc((n)*sizeof(double));
  tf_dx(x,n,k+1,beta,theta);
  double k_fac = glmgen_factorial(k);
  for(i=0; i<n-k-1;i++)
    theta[i] /= k_fac;

  /* Threshold small values */
  for (i=0; i<n-k-1; i++) if (fabs(theta[i])<zero_tol) theta[i]=0;

  /* Compute the predictions at each new point x0 */
  double h;
  for (j=0; j<n0; j++) {
    pred[j] = 0;

    /* Loop over x points, polynomial basis */
    for (i=0; i<k+1; i++) {
      h = 1;
      int l=0;
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
	int l=0;
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

void poly_coefs(double *x, int n, int k, double *beta, double *phi)
{
  memcpy(phi,beta,(k+1)*sizeof(double));
  int j, ell;

  for(j=1; j <= k; j++)
  {
    for(ell = k; ell >= j; ell--)
    {
      phi[ell] = (phi[ell] - phi[ell-1]) / ( x[ell] - x[ell-j] );
    }
  }
}

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
