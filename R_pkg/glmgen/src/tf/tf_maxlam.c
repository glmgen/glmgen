#include "tf.h"

double tf_maxlam (int len, double * y, gqr * Dt_qr, int family)
{
  /* This is exact for the Gaussian case, but only approximate
     for logistic or Poisson losses, which is OK, since we are
     only tasked with finding an interesting range for lambdas. */

  int i;
  double maxlam;
  double * y_work;

  y_work = (double *) malloc(len * sizeof(double));
    
  switch (family)
  {
    case FAMILY_GAUSSIAN:
      for (i = 0; i < len; i++) y_work[i] = y[i];
      break;

    case FAMILY_LOGISTIC:
      for (i = 0; i < len; i++) y_work[i] = y[i];
      break;

    case FAMILY_POISSON:
      for (i = 0; i < len; i++) y_work[i] = y[i];
      break;
  }

  glmgen_qrsol(Dt_qr, y_work);

  maxlam = 0;
  for(i = 0; i < len; i++) maxlam = MAX(maxlam, fabs(y_work[i]));

  free(y_work);
  return maxlam;
}
