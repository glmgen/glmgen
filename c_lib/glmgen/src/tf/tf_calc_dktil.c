#include "tf.h"

/* Creates \tilde{D}k = k. (delta_k)^-1 Dk */
cs * tf_calc_dktil (int n, int k, const double * x)
{
  cs * delta_k;  
  cs * delta_k_cp;
  cs * Dk;
  cs * Dktil;
  
  int i;

  Dk = tf_calc_dk(n, k, x);

  /* Construct diagonal matrix of differences: */
  delta_k = cs_spalloc(n-k, n-k, (n-k), 1, 1);
  for(i = 0; i < n - k; i++)
  {
    delta_k->p[i] = i;
    delta_k->i[i] = i;
    delta_k->x[i] = k / (x[k + i] - x[i]);
  }
  delta_k->nz = n-k;
  delta_k_cp = cs_compress(delta_k);
  Dktil = cs_multiply(delta_k_cp, Dk);

  cs_spfree(Dk);
  cs_spfree(delta_k);
  cs_spfree(delta_k_cp);

  return Dktil;
}
