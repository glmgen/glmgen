#include "tf.h"
#include "utils.h"

/*   Buffer buf should be big enough to hold D^{x,k+1} beta. */

double tf_obj(double *x, double *y, double *w, int n, int k, double lambda, 
    int family, double *beta, double *buf) {

  switch (family)
  {
    case FAMILY_GAUSSIAN:
      return tf_obj_gauss(x,y,w,n,k,lambda,beta,buf);

    case FAMILY_LOGISTIC:
      return tf_obj_glm(x,y,w,n,k,lambda,&logi_b,beta,buf);

    case FAMILY_POISSON:
      return tf_obj_glm(x,y,w,n,k,lambda,&pois_b,beta,buf);
  }

  return 0;
}

double tf_obj_gauss(double *x, double *y, double *w, int n, int k, double lambda, 
    double *beta, double *buf) {

  int i;
  double loss, pen;

  loss = 0;
  for (i=0; i<n; i++) loss += w[i]*(y[i]-beta[i])*(y[i]-beta[i])/2;

  pen = 0;
  tf_dx(x,n,k+1,beta,buf); /* IMPORTANT: use k+1 here! */
  for (i=0; i<n-k-1; i++) pen += fabs(buf[i]);
  return loss+lambda*pen;
}

double tf_obj_glm(double *x, double *y, double *w, int n, int k, double lambda, 
    func_RtoR b, double *beta, double *buf) {
  int i;
  double loss, pen;

  loss = 0;
  for (i=0; i<n; i++) loss += w[i] * (-y[i]*beta[i] + b(beta[i]));

  pen = 0;
  tf_dx(x,n,k+1,beta,buf); /* IMPORTANT: use k+1 here! */
  for (i=0; i<n-k-1; i++) pen += fabs(buf[i]);
  return loss+lambda*pen;
}
