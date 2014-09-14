#include "tf.h"

double ts_maxlam (int len, double * Dy, gqr * DDt_qr, int family)
{
  int i;
  double * Dy_work;
  double maxlam;
  Dy_work = (double *) malloc(len * sizeof(double));

  switch (family)
  {
    case FAMILY_GAUSSIAN:
      for (i = 0; i < len; i++) Dy_work[i] = Dy[i];
      break;

    case FAMILY_LOGISTIC:
      for (i = 0; i < len; i++) Dy_work[i] = Dy[i] - 0.5;
      break;

    case FAMILY_POISSON:
      for (i = 0; i < len; i++) Dy_work[i] = Dy[i] - 1;
      break;
  }

  glmgen_qrsol(DDt_qr, Dy_work);

  maxlam = 0;
  for(i = 0; i < len; i++) maxlam = MAX(maxlam, fabs(Dy_work[i]));

  free(Dy_work);
  return maxlam;
}
