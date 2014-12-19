#include "btree.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* id of the node into which i was merged into */
static int tree_id(int* mt, int n, int i)
{
  int j;
  for(j = i; j < n; )
  {
    if( mt[j] == j )
      return j;
    j = mt[j];
  }
  return -1;
}

static int tree_id_fwd(int* mt, int n, int i)
{
  int j;
  for(j = i+1; j < n; )
  {
    if( mt[j] == j )
      return j;
    j = mt[j];
  }
  return n-1;
}

static int tree_id_bwd(int* mt, int n, int i)
{
  int j;
  for(j = i-1; j >= 0; j-- )
  {
    if( mt[j] == j )
      return j;
  }
  return -1;
}
void presmooth( double* x, double* y, double* w, int n, int k, 
  double** xt, double** yt, double** wt, int* nt_ptr, double x_cond)
{
  int i;
  double * dx;
  int * mt;/* merge_trail: mt[i] is the id into which i is merged into */
  double * wtmp;
  double * ytmp;
  int j, jn, jp;
  double dxj;
  int nt;
  btnode * t; 
  btnode * min_val_node;
  btnode * min_id_node;

  /* n+1 because we want to record the displacement of the first and last points also */ 
  dx = (double*)malloc( (n+1) * sizeof(double) );
  mt = (int*)malloc( n * sizeof(int) );
  wtmp = (double*)malloc( n * sizeof(double) );
  ytmp = (double*)malloc( n * sizeof(double) );

  t = min_val_node = NULL;

  printf("Presmoothing initialization\n");

  /* Initialization */
  dx[0] = dx[n] = 0;
  for(i=0; i<n-1; i++) dx[i+1] = x[i+1] - x[i];
  for(i=0; i<n; i++) wtmp[i] = w[i];
  for(i=0; i<n; i++) ytmp[i] = y[i];
  for(i=0; i<n; i++) mt[i] = i;
  for(i=0; i<n-1; i++) bt_insert(&t, i, dx[i+1]);

  printf("Presmoothing initialization done\n");
  double range;

  nt = n;
  while(1) {
    if( t == NULL ) break;
    // if( nt % 1000 == 0) printf("%d ", n-nt+1);
    bt_find_min(t, &min_val_node);

    dxj = min_val_node->key;
    bt_find_min(min_val_node->ids, &min_id_node);
    j = (int) min_id_node->key;
    
    range = x[n-1]-x[0]-dx[0]-dx[n];
    if( nt * pow( range / dxj, k+1) < x_cond ) break;
    //if(dxj >= delta) break;
    
    // jn = -1 when j is the first node in the current tree. update dx[0]
    // jp = n-1 when j is the last node. update dx[n]
    // Note: mt[n-1] = n-1 always
    jn = tree_id_bwd(mt, n, j);
    jp = tree_id_fwd(mt, n, j);

    /* Merge j and jp into jp */
    bt_delete(&t, j, dx[j+1]);
    if(jn != -1)  bt_delete(&t, jn, dx[jn+1]);
    if(jp != n-1) bt_delete(&t, jp, dx[jp+1]);

    dx[jn+1] += dx[j+1]/2.;
    dx[jp+1] += dx[j+1]/2.;
    
    if(jn != -1)  bt_insert(&t, jn, dx[jn+1]);
    if(jp != n-1) bt_insert(&t, jp, dx[jp+1]);

    /* Average y */
    ytmp[jp] = (wtmp[j]*ytmp[j] + wtmp[jp]*ytmp[jp])/(wtmp[j] + wtmp[jp]);
    
    wtmp[jp] += wtmp[j];
    mt[j] = jp;
    wtmp[j] = 0;
    nt--;
  }

  //printf("nt from #merges = %d\n", nt);

  int* ids = (int*)malloc( nt * sizeof(int) );
  int ntt = 0;
  for(i=0; i < n; i++) {  
    if(mt[i] == i) {
      if( ntt >= nt ) { // Sanity check
        printf("ERROR: nt=%d but ntt >%d\n", nt, ntt);
        break;
      }
      ids[ntt] = i;
      ntt++;
    }
  }

  /* Fill-up outputs */
  *nt_ptr = nt;

  // build responses xt, yt, wt
  *xt = (double*)malloc( nt * sizeof(double) );
  *yt = (double*)malloc( nt * sizeof(double) );
  *wt = (double*)malloc( nt * sizeof(double) );

  // Get the weights at the points remaining after merging
  for(i=0; i < nt; i++) (*wt)[i] = wtmp[ ids[i] ];

  (*xt)[0] = x[0] + dx[0];
  for(i = 1; i < nt; i++) (*xt)[i] = (*xt)[i-1] + dx[ ids[i-1] + 1 ];
  for(i = 0; i < nt; i++) (*yt)[i] = ytmp[ ids[i] ];

  free(dx);
  free(mt);
  free(wtmp);
  free(ytmp);
  free(ids);
  bt_free(t);
}

