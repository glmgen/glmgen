#include "tf.h"

double tf_maxlam (int len, double * y, gqr * Dt_qr, double * w, int family)
{
  /* This is exact for the Gaussian case, but only approximate
     for logistic or Poisson losses, which is OK, since we are
     only tasked with finding an interesting range for lambdas. */

  int i;
  int tmp; /*n-k-1*/
  double maxlam;
  double * y_work;

  y_work = (double *) malloc(len * sizeof(double));
    
  for(i = 0; i < len; i++) y_work[i] = sqrt(w[i]) * y[i];

  glmgen_qrsol(Dt_qr, y_work);
  
  tmp = Dt_qr->n;
  maxlam = 0;
  for(i = 0; i < tmp; i++) maxlam = MAX(maxlam, fabs(y_work[i]));

  free(y_work);
  return maxlam;
}
