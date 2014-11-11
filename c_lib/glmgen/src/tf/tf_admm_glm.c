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
  
  /* Buffer for line search. Dbn : Dbnew */
  double * Db = (double*)malloc(n*sizeof(double));
  double * Dd = (double*)malloc(n*sizeof(double));
  double * Dbn = (double*)malloc(n*sizeof(double));
  int * iter_ls = (int*)malloc(sizeof(int));
  
  cs * Dk;
  cs * Dkt;
  cs * DktDk;

  cs * kernmat;
  gqr * kernmat_qr;
  
  double pobj, loss, pen;
  double t; /* stepsize */
  int it;
  double admm_tol;
  int max_iter_admm;

  int vary_rho;

  vary_rho = 0; /* TODO pass this as an argument */

  /* Set max_iter, tol appropriately */
  admm_tol = obj_tol * 1e3;
  max_iter_admm = max_iter;

  obj_admm = (double*)malloc(max_iter_admm*sizeof(double)); 
  
  int verb = 0; 
  if (verb) printf("Iteration\tObjective\tLoss\t\tPenalty\n");

  Dk = tf_calc_dktil(n, k, x);
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
      /* IMPORTANT: set yt = Wyt */
      yt[i] = H[i]*beta[i] + y[i] - b1(beta[i]);

    }
    
    if( verb ){
      /* for(i=0; i<n; i++) printf("beta=%0.2e\tH=%0.2e\n", beta[i], H[i]); */
    }
    
    /* QR factorization */
    kernmat = scalar_plus_diag(DktDk, rho, H);
    kernmat_qr = glmgen_qr(kernmat);  

    /* prox-Newton step */
    int iter_admm = 0;
    tf_admm_gauss (yt, x, H, n, k,
       max_iter_admm, lam,
       d, alpha, u,
       obj_admm, &iter_admm, rho, vary_rho, admm_tol,
       kernmat_qr, DktDk);

    cs_spfree(kernmat);
    glmgen_gqr_free(kernmat_qr);
    
    for(i=0; i<n; i++)
    {
      d[i] = d[i] - beta[i];
    }
   
    /* TODO: take the bt line search parameters as inputs */
    double alpha_ls = 0.4;
    double gamma = 0.5;
    int max_iter_ls = 50;

    t = tf_line_search(y, x, n, k, lam, b, b1, beta, d, alpha_ls, gamma, max_iter_ls, iter_ls, Db, Dd, Dbn); 

    /* if(verb) printf("Stepsize t=%.2e,\titers=%d\n", t, *iter_ls); */
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

    if (verb) printf("GLM \t%i\t%0.3e\t%0.3e\t%0.3e\n",it,pobj, loss, lam*pen);

    if(it > 0)
    {
      if( fabs(pobj - obj[it-1]) < fabs(pobj) * obj_tol )
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
  free(iter_ls);
  free(Db);
  free(Dd);
  free(Dbn);
  free(obj_admm);
  cs_spfree(Dk);
  cs_spfree(Dkt);
  cs_spfree(DktDk);
}
