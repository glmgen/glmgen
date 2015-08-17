#include "tf.h"
#include "utils.h"

/* Df will be filled with (# of nonzeros in D^{x,k+1} beta) + k + 1.
   Buffer should be long enough to hold D^{x,k+1} beta. */

double tf_obj(double *x, double *y, double *w, int n, int k, double lambda, 
	      int family, double *beta, int *df, double *buf) {

  double loss, pen;	
  int i, dfval;
  loss = 0; pen = 0;
  dfval = k+1;

  switch (family)
    {
    case FAMILY_GAUSSIAN:
      for (i=0; i<n; i++) loss += w[i]*(y[i]-beta[i])*(y[i]-beta[i])/2;
      break;
    
    case FAMILY_LOGISTIC:
      for (i=0; i<n; i++) loss += w[i]*(-y[i]*beta[i] + logi_b(beta[i]));
      break;
      
    case FAMILY_POISSON:
      for (i=0; i<n; i++) loss += w[i]*(-y[i]*beta[i] + pois_b(beta[i])); 
      break;
      
    default:
      return 0;
    }
  
  tf_dx(x,n,k+1,beta,buf); /* IMPORTANT: use k+1 here! */
  for (i=0; i<n-k-1; i++) {
    pen += fabs(buf[i]);
    if (fabs(buf[i]) > 1e-4) dfval += 1;
  }

  *df = dfval;
  double obj = loss+lambda*pen;
  return obj;
}

