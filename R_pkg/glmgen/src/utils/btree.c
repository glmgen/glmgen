#include "btree.h"
#include <stdio.h>
#include <stdlib.h>

#define TRUE 1
#define FALSE 0

/* inserts a new node in a binary search tree */
void bt_insert ( btnode **bt, int id, double val)
{
  if ( *bt == NULL )
  {
    *bt = (btnode*)malloc (sizeof(btnode));

    ( *bt ) -> leftchild = NULL;
    ( *bt ) -> ids = create_node(id);
    ( *bt ) -> val = val;
    ( *bt ) -> rightchild = NULL ;
  }
  else/* search the node to which new node will be attached */
  {
    /* if new val is less, traverse to left */
    if ( val < ( *bt ) -> val )
      bt_insert ( &( ( *bt ) -> leftchild ), id, val ) ;
    else if( val > ( *bt ) -> val ) /* else traverse to right */
      bt_insert ( &( ( *bt ) -> rightchild ), id, val ) ;
    else
      insert_node( &(( *bt ) -> ids), id );
  }
}
/* There should not be a tree node with empty id list after any insert/delete */
/* deletes a node from the binary search tree */
void bt_delete ( btnode **bt, int id, double val )
{
  int found = FALSE;
  btnode *parent, *x, *xsucc ;

  /* if tree is empty */
  if ( *bt == NULL )
  {
    printf ( "\nTree is empty" ) ;
    return ;
  }

  parent = x = NULL ;

  /* call to search function to find the node to be deleted */

  bt_search( bt, id, val, &parent, &x, &found );

  if ( found == FALSE )
    return ;

  delete_node(&(x->ids), id);

  /* if x still has some ids, do not delete it */
  if( !isempty(x -> ids) )
    return;
  
  /* if the node to be deleted has two children */
  if ( x -> leftchild != NULL && x -> rightchild != NULL )
  {
    parent = x ;
    xsucc = x -> rightchild ;

    while ( xsucc -> leftchild != NULL )
    {
      parent = xsucc ;
      xsucc = xsucc -> leftchild ;
    }

    x -> val = xsucc -> val;
    x -> ids = xsucc -> ids; 
    x = xsucc ; /* delete xsucc now */
  }

  /* if the node to be deleted has no child */
  if ( x -> leftchild == NULL && x -> rightchild == NULL )
  {
    if( parent == NULL )
      *bt = NULL;
    else if ( parent -> rightchild == x )
      parent -> rightchild = NULL ;
    else
      parent -> leftchild = NULL ;

    free ( x ) ;
    return ;
  }

  /* if the node to be deleted has only rightchild */
  if ( x -> leftchild == NULL && x -> rightchild != NULL )
  {
    if( parent == NULL )
      (*bt ) = x -> rightchild;    
    else if ( parent -> leftchild == x )
      parent -> leftchild = x -> rightchild ;
    else
      parent -> rightchild = x -> rightchild ;

    free ( x ) ;
    return ;
  }

  /* if the node to be deleted has only left child */
  if ( x -> leftchild != NULL && x -> rightchild == NULL )
  {
    if( parent == NULL )
      (*bt ) = x -> leftchild;    
    else if ( parent -> leftchild == x )
      parent -> leftchild = x -> leftchild ;
    else
      parent -> rightchild = x -> leftchild ;

    free ( x ) ;
    return ;
  }
}

/*returns the address of the node to be deleted, address of its parent and
 *    whether the node is found or not */
void bt_search( btnode **bt, int id, double val,
    btnode **par, btnode **x, int *found )
{
  btnode *q;

  q = *bt;
  *found = FALSE;
  *par = NULL;

  while ( q != NULL )
  {
    if ( q -> val > val ) {
      *par = q;
      q = q -> leftchild ;
    }
    else if( q -> val < val ) {
      *par = q;
      q = q -> rightchild ;
    }
    else
    {
      *found = TRUE;
      /* delete the id */
      /* delete_node(&(q->ids), id); */
      *x = q;
      return;
    }
  }
}

/* traverse a binary search tree in a LDR (Left-id-Right) fashion */
void bt_inorder ( btnode *bt )
{
  if ( bt != NULL )
  {
    bt_inorder ( bt -> leftchild ) ;

    printf("%g", bt -> val ) ;
    printf("["); display( bt -> ids ); printf("]  ");

    bt_inorder ( bt -> rightchild ) ;
  }
}

void bt_find_min ( btnode* bt, btnode** x )
{
  *x = bt;
  if( bt == NULL ) 
    return;
  
  while( (*x) -> leftchild != NULL )
    *x = (*x) -> leftchild;
}

int bt_num_ids( btnode *bt )
{  
  if( bt != NULL )
    return length(bt->ids) + bt_num_ids(bt->leftchild) + bt_num_ids(bt->rightchild);
  else
    return 0;
}

void bt_free( btnode *bt )
{
  if ( bt != NULL )
  {
    bt_free( bt -> leftchild );
    bt_free( bt -> rightchild );
    ll_free( bt -> ids );
    free( bt );
  }
}
