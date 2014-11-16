#include "tf.h"
#include "math.h"

void tf_admm_glm (double * y, double * x, double * w, int n, int k,
       int max_iter, double lam,
       double * beta, double * alpha, double * u,
       double * obj, int * iter,
       double rho, double obj_tol, cs * DktDk,
       func_RtoR b, func_RtoR b1, func_RtoR b2)
{
  
  double * d = (double*)malloc(n*sizeof(double)); /* line search direction */
  double * yt = (double*)malloc(n*sizeof(double));/* working response: ytilde */
  double * H = (double*)malloc(n*sizeof(double)); /* weighted Hessian */
  double * z = (double*)malloc(n*sizeof(double));
  double * obj_admm;
  
  /* Buffers for line search */
  double * Db = (double*)malloc(n*sizeof(double));
  double * Dd = (double*)malloc(n*sizeof(double));
  int * iter_ls = (int*)malloc(sizeof(int));
  
  double pobj, loss, pen;
  double t; /* stepsize */
  int it;
  double admm_tol;
  int max_iter_admm;

  /* Set ADMM parameters appropriately */
  admm_tol = obj_tol * 10;
  max_iter_admm = ADMM_MAX_ITER;

  obj_admm = (double*)malloc(max_iter_admm*sizeof(double)); 
  
  int verb = 1; 
  if (verb) printf("Iteration\tObjective\tLoss\t\tPenalty\n");
 
  /* One Prox Newton step per iteration */
  for (it=0; it < max_iter; it++)
  {    
    /* Define weighted Hessian, and working response */
    int i;
    for(i=0; i<n; i++)
    {
      H[i] = w[i] * b2(beta[i]);
      if (fabs(H[i])>WEIGHT_SMALL)
      {
	yt[i] = beta[i] - (y[i]-b1(beta[i]))/H[i];
      }
      else 
      {
	yt[i] = beta[i] - (y[i]-b1(beta[i]));
      }
      /* IMPORTANT: set yt = Wyt */
      //yt[i] = H[i]*beta[i] + y[i] - b1(beta[i]);
    }
    
    if( verb ){
      /* for(i=0; i<n; i++) printf("beta=%0.2e\tH=%0.2e\n", beta[i], H[i]); */
    }
    
    /* Prox Newton step */
    int iter_admm = 0;
    tf_admm_gauss (yt, x, H, n, k,
       max_iter_admm, lam,
       d, alpha, u,
       obj_admm, &iter_admm, rho, admm_tol,
       DktDk);

    for(i=0; i<n; i++)
    {
      d[i] = d[i] - beta[i];
    }
   
    /* TODO: take the bt line search parameters as inputs */
    double alpha_ls = 0.5;
    double gamma = 0.8;
    int max_iter_ls = 50;

    t = tf_line_search(y, x, w, n, k, lam, b, b1, beta, d, alpha_ls, gamma, max_iter_ls, iter_ls, Db, Dd); 
    //t = 0.01;

    if(verb) printf("Stepsize t=%.2e,\titers=%d\n", t, *iter_ls); 
    for(i=0; i<n; i++)
    {
      beta[i] = beta[i] + t * d[i];
    }    

    /* Compute objective */
    /* Compute loss */
    loss = 0;
    for (i=0; i<n; i++)
    {
      loss += w[i] * (-y[i]*beta[i] + b(beta[i]));
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
  free(obj_admm);
}
