#include <stdio.h>

#include "cs.h"
#include "utils.h"
#include "tf.h"

int main()
{
  int n;
  int k;
  int i;
  double * x;
  double * y;
  double * b;
  double * Dy;

  /* Input parameters; will usually be given by parent function to tf_admm */
  n = 6;
  k = 0;

  x = (double *) malloc(n * sizeof(double));
  y = (double *) malloc(n * sizeof(double));
  b = (double *) malloc((n-k-2) * sizeof(double));
  Dy = (double *) malloc((n-k-2) * sizeof(double));

  for (i = 0; i < n; i++) x[i] = i;
  for (i = 0; i < n; i++) y[i] = i + (i > 2)*3;

  /* Allocate D matrices that we want to test */
  cs * D;
  cs * Dt;
  cs * DtD;
  cs * DDt;
  gqr * DDt_qr;
  cs * kernel_mat;

  /* Calculate D_(k+1), D_(k+1)^T * D_(k+1), D_(k+1) * D_(k+1)^t, and Dy */
  D = tf_calc_dk(n, k+1, x);
  Dt = cs_transpose(D, 1);
  DtD = cs_multiply(Dt,D);
  DDt = cs_multiply(D,Dt);
  cs_gaxpy (D, y, Dy);

  /* Print them! */
  printf("\n---------- D_(k+1) -----------------------------\n");
  cs_print(D, 0);
  printf("\n---------- D_(k+1)^T ---------------------------\n");
  cs_print(Dt, 0);
  printf("\n---------- D_(k+1)^T * D_(k+1) -----------------\n");
  cs_print(DtD, 0);
  printf("\n---------- D_(k+1) * D_(k+1)^t -----------------\n");
  cs_print(DDt, 0);
  printf("\n---------- y (for comparison) ------------------\n");
  for (i = 0; i < n; i++) printf("%f\n", y[i]);
  printf("\n---------- D * y -------------------------------\n");
  for (i = 0; i < n-k-2; i++) printf("%f\n", Dy[i]);

  /* Compute the sparse qr of DDt, and solve DDt b = Dy */
  DDt_qr = glmgen_qr(DDt);
  for (i = 0; i < n-k-2; i++) b[i] = Dy[i];
  glmgen_qrsol(DDt_qr, b);

  /* Print the solution */
  printf("\n---------- b solution from -> DDt b = Dy -------\n");
  for (i = 0; i < n-k-2; i++) printf("%f\n", b[i]);

  /* Compute the kernel matrix for rho=0.5, and print it */
  kernel_mat = scalar_plus_eye(DDt, 0.5);
  printf("\n---------- 0.5 * D_(k+1) * D_(k+1)^t + I_(k) ---\n");
  cs_print(kernel_mat, 0);

  /* Free the csparse matricies */
  cs_spfree(D);
  cs_spfree(Dt);
  cs_spfree(DtD);
  cs_spfree(DDt);
  glmgen_gqr_free(DDt_qr);

  /* Free the allocated positions 'x', responses 'y', and inner product 'Dy' */
  free(x);
  free(y);
  free(b);
  free(Dy);
}
