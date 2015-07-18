#include "tf.h"

/* Buffer should be long enough to hold D^{x,k+1} beta. */
double tf_obj(double *y, double *x, double *w, int n, int k, double lambda, 
		double *beta, double *buf) {

	double loss, pen;	
	int i;

	loss = 0; pen = 0;
  for (i=0; i<n; i++) loss += w[i]*(y[i]-beta[i])*(y[i]-beta[i]);
  tf_dx(x,n,k+1,beta,buf); /* IMPORTANT: use k+1 here! */
  for (i=0; i<n-k-1; i++) pen += fabs(buf[i]);
  double obj = loss/2+lambda*pen;
	return obj;

}

