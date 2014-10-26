#include "tf.h"

double ts_maxlam (int len, double * y, gqr * Dt_qr, int family)
{
  int i;
  double * y_work;
  double maxlam;
  y_work = (double *) malloc(len * sizeof(double));

  switch (family)
  {
    case FAMILY_GAUSSIAN:
      for (i = 0; i < len; i++) y_work[i] = y[i];
      break;

    case FAMILY_LOGISTIC:
      for (i = 0; i < len; i++) y_work[i] = y[i] - 0.5;
      break;

    case FAMILY_POISSON:
      for (i = 0; i < len; i++) y_work[i] = y[i] - 1;
      break;
  }

  glmgen_qrsol(Dt_qr, y_work);

  maxlam = 0;
  for(i = 0; i < len; i++) maxlam = MAX(maxlam, fabs(y_work[i]));

  free(y_work);
  return maxlam;
}
