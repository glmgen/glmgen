#include "tf.h"

/* Creates the Dk matrix for a given n,k, and x */
cs * tf_calc_dk (int n, int k, const double * x)
{
  long int i;
  int tk = 1; /* "this k" (actually, k+1) - will iterate until ts = k+1 */

  cs * D1;
  cs * D1_cp;
  cs * Dk;
  cs * Dk_cp;
  cs * delta_k;
  cs * delta_k_cp;
  cs * D1_x_delta;
  cs * Dk_next;

  /* Contruct one 'full D1', which persists throughout
     and another copy as Dk */
  D1 = cs_spalloc(n-tk, n, (n-tk)*2, 1, 1);
  Dk = cs_spalloc(n-tk, n, (n-tk)*2, 1, 1);
  D1->nz = (n-tk)*2;
  Dk->nz = (n-tk)*2;
  for (i = 0; i < (n-tk)*2; i++)
  {
     D1->p[i] = (i+1) / 2;
     Dk->p[i] = D1->p[i];
     D1->i[i] = i / 2;
     Dk->i[i] = D1->i[i];
     D1->x[i] = -1 + 2*(i % 2);
     Dk->x[i] = D1->x[i];
  }

  /* Create a column compressed version of Dk, and
     delete the old copy */
  Dk_cp = cs_compress(Dk);
  cs_spfree(Dk);

  for (tk = 1; tk < k+1; tk++)
  {
    /* 'reduce' the virtual size of D1 to: (n-tk-1) x (n-tk),
       compress into compressed column, saving as D1_cp */
    D1->nz = (n-tk-1)*2;
    D1->m = n-tk-1;
    D1->n = n-tk;
    D1_cp = cs_compress(D1);

    /* Construct diagonal matrix of differences: */
    delta_k = cs_spalloc(n-tk, n-tk, (n-tk), 1, 1);
    for(i = 0; i < n - tk; i++)
    {
      delta_k->p[i] = i;
      delta_k->i[i] = i;
      delta_k->x[i] = tk / (x[tk + i] - x[i]);
    }
    delta_k->nz = n-tk;
    delta_k_cp = cs_compress(delta_k);
    D1_x_delta = cs_multiply(D1_cp, delta_k_cp);

    /* Execute the matrix multiplication */
    Dk_next = cs_multiply(D1_x_delta, Dk_cp);

    /* Free temporary cs matricies created in each loop */
    cs_spfree(D1_cp);
    cs_spfree(delta_k);
    cs_spfree(delta_k_cp);
    cs_spfree(D1_x_delta);
    cs_spfree(Dk_cp);
    Dk_cp = Dk_next;
  }

  cs_spfree(D1);
  return Dk_cp;
}
