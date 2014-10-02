#include "tf.h"
#include "math.h"

void tf_admm_glm (double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj,
       double rho, double obj_tol,
       gqr * sparseQR,
       func_RtoR b, func_RtoR b1, func_RtoR b2)
{
  
  double * d = (double*)malloc(n*sizeof(double)); /* line search direction */
  double * yt = (double*)malloc(n*sizeof(double));/* intermediate y: ytilde */
  double * H = (double*)malloc(n*sizeof(double)); /* Hessian */
  double * obj_admm = (double*)malloc(max_iter*sizeof(double)); 
  double * z = (double*)malloc(n*sizeof(double));
  
  cs * Dk;
  cs * Dkt;
  cs * DktDk;

  cs * kernmat;
  gqr * kernmat_qr;
  
  double pobj, loss, pen;
  
  int verb = 1; 
  if (verb) printf("Iteration\tObjective");

  Dk = tf_calc_dk(n, k-1, x);
  Dkt = cs_transpose(Dk, 1);
  DktDk = cs_multiply(Dkt,Dk);

  int iter;  
  /* One prox-Newton step per iteration */
  for (iter=1; iter<=max_iter; (iter)++)
  {    
    /* set yt with weights and H */
    int i;
    for(i=0; i<n; i++)
    {
      H[i] = b2(beta[i]);
      yt[i] = H[i]*beta[i] + y[i] - b1(beta[i]);      
    }
    
    /* QR factorization */
    kernmat = scalar_plus_diag(DktDk, rho, H);
    kernmat_qr = glmgen_qr(kernmat);  

    /* prox-Newton step */
    /* Set max_iter, tol appropriately */
    double admm_tol = obj_tol * 1e3;

    tf_admm_gauss (yt, x, H, n, k,
       max_iter, lam,
       d, alpha, u,
       obj_admm, rho, admm_tol,
       kernmat_qr);

    cs_spfree(kernmat);
    glmgen_gqr_free(kernmat_qr);
    
    for(i=0; i<n; i++)
    {
      d[i] = d[i] - beta[i];
    }
    
    /* bt linesearch. Also see Dinh et. al ICML 2013*/
    /* TODO: f( x + gamma^s ) \leq f(x) + c gamma^s < f', d> */
    /* double t = linesearch_bt(y, x, H, n, k, d, beta, c, gamma, tol); */
    double t = 1;
    
    for(i=0; i<n; i++)
    {
      beta[i] = beta[i] + t * d[i];      
    }    

    /* Compute objective, if we are told to
     * Compute loss */
    loss = 0;
    for (i=0; i<n; i++)
    {
        loss += -y[i]*beta[i] + b(beta[i]);
    }
    /* Compute penalty */
    tf_dx(x,n,k+1,beta,z); /* IMPORTANT: use k+1 here! */
    pen = 0;
    for (i=0; i<n-k-1; i++)
    {
      pen += fabs(z[i]);
    }
    pobj = loss+lam*pen;
    obj[(iter)-1] = pobj;

    if (verb) printf("%i\t%0.5e\n",iter,pobj);

    if(iter > 1)
    {
      if( pobj < obj_tol || (obj[(iter)-1] - pobj) / pobj < obj_tol )
      {
        break;
      }
    }
  }
  
  /* free */  
  free(d);
  free(yt);
  free(H);
  free(obj_admm);
  free(z);
  cs_spfree(Dk);
  cs_spfree(Dkt);
  cs_spfree(DktDk);
}
