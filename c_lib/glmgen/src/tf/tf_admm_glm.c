#include "tf.h"
#include "math.h"

void tf_admm_glm (double * y, double * x, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj, int * iter,
       double rho, double obj_tol,
       func_RtoR b, func_RtoR b1, func_RtoR b2)
{
  
  double * d = (double*)malloc(n*sizeof(double)); /* line search direction */
  double * yt = (double*)malloc(n*sizeof(double));/* intermediate y: ytilde */
  double * H = (double*)malloc(n*sizeof(double)); /* Hessian */
  double * z = (double*)malloc(n*sizeof(double));
  double * obj_admm;
  
  cs * Dk;
  cs * Dkt;
  cs * DktDk;

  cs * kernmat;
  gqr * kernmat_qr;
  
  double pobj, loss, pen;
  int it;
  double admm_tol;
  int max_iter_admm;

  /* Set max_iter, tol appropriately */
  admm_tol = obj_tol * 1e3;
  max_iter_admm = max_iter;

  obj_admm = (double*)malloc(max_iter_admm*sizeof(double)); 
  
  int verb = 1; 
  if (verb) printf("Iteration\tObjective");

  Dk = tf_calc_dk(n, k, x);
  Dkt = cs_transpose(Dk, 1);
  DktDk = cs_multiply(Dkt,Dk);
 
  
  /* One prox-Newton step per iteration */
  for (it=0; it < max_iter; it++)
  {    
    /* set yt with weights and H */
    int i;
    for(i=0; i<n; i++)
    {
      H[i] = b2(beta[i]);
      /*
      if( fabs(H[i]) > 1e-12 )
        yt[i] = beta[i] + (y[i] - b1(beta[i])) / H[i];
      else
        yt[i] = 0;
      */
      /* set yt = Wyt */
      yt[i] = H[i]*beta[i] + y[i] - b1(beta[i]);

    }
    
    /*for(i=0; i<n; i++) printf("beta=%0.2e\tH=%0.2e\n", beta[i], H[i]);*/
    
    /* QR factorization */
    kernmat = scalar_plus_diag(DktDk, rho, H);
    kernmat_qr = glmgen_qr(kernmat);  

    /* prox-Newton step */
    int iter_admm = 0;
    tf_admm_gauss (yt, x, H, n, k,
       max_iter_admm, lam,
       d, alpha, u,
       obj_admm, &iter_admm, rho, admm_tol,
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
    obj[it] = pobj;

    if (verb) printf("GLM \t%i\t%0.5e\n",it,pobj);

    if(it > 0)
    {
      if( fabs(pobj) < obj_tol || fabs((obj[it-1] - pobj) / pobj) < obj_tol )
      {
        break;
      }
    }
  }
  
  *iter = it;

  /* free */  
  free(d);
  free(yt);
  free(H);
  free(z);
  free(obj_admm);
  cs_spfree(Dk);
  cs_spfree(Dkt);
  cs_spfree(DktDk);
}
