#include <stdio.h>

#include "cs.h"
#include "utils.h"
#include "tf.h"

#define PI 3.14159265

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
  k = 3;

  x = (double *) malloc(n * sizeof(double));
  y = (double *) malloc(n * sizeof(double));
  b = (double *) malloc(n * sizeof(double));
  Dy = (double *) calloc(n, sizeof(double));

  for (i = 0; i < n; i++) x[i] = i/(double)(n-1);
  for (i = 0; i < n; i++) y[i] = i + (i > 2)*3;
  for (i = 0; i < n; i++) y[i] = sin(x[i] * 3.*PI);

  /* Allocate D matrices that we want to test */
  cs * D;
  cs * Dt;
  cs * DtD;
  cs * DDt;
  gqr * Dt_qr;
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
  printf("\n---------- y (for comparison) ------------------\n");
  for (i = 0; i < n; i++) printf("%f\n", y[i]);
  /*printf("\n---------- D_(k+1)^T * D_(k+1) -----------------\n");
  cs_print(DtD, 0);
  printf("\n---------- D_(k+1) * D_(k+1)^t -----------------\n");
  cs_print(DDt, 0);
  printf("\n---------- D * y -------------------------------\n");
  for (i = 0; i < n-k-1; i++) printf("%f\n", Dy[i]);
*/

  /* Compute the sparse qr of DDt, and solve DDt b = Dy */
  Dt_qr = glmgen_qr(Dt);
  for (i = 0; i < n; i++) b[i] = y[i];
  glmgen_qrsol(Dt_qr, b);

  /* Print the solution */
  printf("\n---------- b solution from -> Dt b = y -------\n");
  for (i = 0; i < n; i++) printf("%f\n", b[i]);

  /* Compute the kernel matrix for rho=0.5, and print it */
  kernel_mat = scalar_plus_eye(DDt, 0.5);
/*  printf("\n---------- 0.5 * D_(k+1) * D_(k+1)^t + I_(k) ---\n");
  cs_print(kernel_mat, 0);
*/
  /* Free the csparse matricies */
  cs_spfree(D);
  cs_spfree(Dt);
  cs_spfree(DtD);
  cs_spfree(DDt);
  cs_spfree(kernel_mat);
  glmgen_gqr_free(Dt_qr);


  /* Free the allocated positions 'x', responses 'y', and inner product 'Dy' */
  free(x);
  free(y);
  free(b);
  free(Dy);

  return 0;
}
