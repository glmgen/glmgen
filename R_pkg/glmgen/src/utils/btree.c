#include "btree.h"
#include <stdio.h>
#include <stdlib.h>

#define TRUE 1
#define FALSE 0

/* inserts a new node in a binary search tree */

void bt_insert (btnode **bt, int id, double val)
{
  if ( *bt == NULL )
  {
    *bt = (btnode*)malloc (sizeof(btnode));

    (*bt)->leftchild = NULL;
    btnode* ids = NULL;
    bt_insert_inner(&ids, id);
    (*bt)->ids = ids;
    (*bt)->key = val;
    (*bt)->node_lvl = FIRST;
    (*bt)->rightchild = NULL;
  }
  else/* search the node to which new node will be attached */
  {
    /* if new val is less, traverse to left */
    if ( val < (*bt)->key )
      bt_insert ( &( ( *bt ) -> leftchild ), id, val ) ;
    else if( val > (*bt)->key ) /* else traverse to right */
      bt_insert ( &( ( *bt ) -> rightchild ), id, val ) ;
    else
      bt_insert_inner( &(( *bt ) -> ids), id );
  }
}
void bt_insert_inner (btnode **bt, int id)
{
  if ( *bt == NULL )
  {
    *bt = (btnode*)malloc (sizeof(btnode));

    (*bt)->leftchild = NULL;
    (*bt)->ids = NULL;
    (*bt)->key = id;
    (*bt)->node_lvl = SECOND;
    (*bt)->rightchild = NULL;
  }
  else /* search the node to which new node will be attached */
  {
    if ( id < (*bt)->key )
      bt_insert_inner ( &( ( *bt ) -> leftchild ), id );
    else if( id > (*bt)->key ) 
      bt_insert_inner ( &( (*bt)->rightchild ), id );
    /* do nothing if id already exists */
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
    return ;

  parent = x = NULL ;

  /* call to search function to find the node to be deleted */
  bt_search( bt, val, &parent, &x, &found );

  if ( found == FALSE )
    return ;

  bt_delete_inner(&(x->ids), id);

  /* if x still has some ids, do not delete it */
  if(x->ids != NULL)
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

    x -> key = xsucc -> key;
    x -> ids = xsucc -> ids; 
    x = xsucc ; /* delete xsucc now */
  }

  bt_delete_found_node(bt, &parent, x);

}

void bt_delete_inner(btnode **bt, int id)
{
  int found = FALSE;
  btnode *parent, *x, *xsucc ;

  if(*bt == NULL)
    return ;

  parent = x = NULL ;

  /* call to search function to find the node to be deleted */
  bt_search( bt, id, &parent, &x, &found );

  if ( found == FALSE )
    return ;

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

    x -> key = xsucc -> key;
    x = xsucc ; /* delete xsucc now */
  }
  
  bt_delete_found_node(bt, &parent, x);
}

void bt_delete_found_node(btnode **bt, btnode **parent, btnode *x)
{  
  /* if the node to be deleted has no child */
  if ( x -> leftchild == NULL && x -> rightchild == NULL )
  {
    if( *parent == NULL )
      *bt = NULL;
    else if ( (*parent) -> rightchild == x )
      (*parent) -> rightchild = NULL ;
    else
      (*parent) -> leftchild = NULL ;

    free ( x ) ;
    return ;
  }

  /* if the node to be deleted has only rightchild */
  if ( x -> leftchild == NULL && x -> rightchild != NULL )
  {
    if( (*parent) == NULL )
      (*bt ) = x -> rightchild;    
    else if ( (*parent) -> leftchild == x )
      (*parent) -> leftchild = x -> rightchild ;
    else
      (*parent) -> rightchild = x -> rightchild ;

    free ( x ) ;
    return ;
  }

  /* if the node to be deleted has only left child */
  if ( x -> leftchild != NULL && x -> rightchild == NULL )
  {
    if( (*parent) == NULL )
      (*bt ) = x -> leftchild;    
    else if ( (*parent) -> leftchild == x )
      (*parent) -> leftchild = x -> leftchild ;
    else
      (*parent) -> rightchild = x -> leftchild ;

    free ( x ) ;
    return ;
  }
}
/*returns the address of the node to be deleted, address of its parent and
 *    whether the node is found or not */
void bt_search( btnode **bt, double val,
    btnode **par, btnode **x, int *found )
{
  btnode *q;

  q = *bt;
  *found = FALSE;
  *par = NULL;
  
  while ( q != NULL )
  {
    if ( q -> key > val ) {
      *par = q;
      q = q -> leftchild ;
    }
    else if( q -> key < val ) {
      *par = q;
      q = q -> rightchild ;
    }
    else
    {
      *found = TRUE;
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

    printf(" %g", bt -> key );
    if( bt->node_lvl == FIRST)
    {
      printf("["); bt_inorder( bt -> ids ); printf(" ] ");
    }

    bt_inorder ( bt -> rightchild ) ;
  }
}

void bt_find_min(btnode* bt, btnode** x)
{
  *x = bt;
  if( bt == NULL ) 
    return;

  while( (*x) -> leftchild != NULL )
    *x = (*x) -> leftchild;
}

void bt_find_min_twice(btnode* bt, btnode** x)
{
  btnode* xtmp;

  bt_find_min(bt, &xtmp);
  if( xtmp == NULL )
    return;

  bt_find_min(xtmp->ids, x);
}

void bt_free( btnode *bt )
{
  if ( bt != NULL )
  {
    bt_free( bt->leftchild );
    bt_free( bt->rightchild );
    if( bt->node_lvl == FIRST )
      bt_free( bt->ids );
    free( bt );
  }
}
