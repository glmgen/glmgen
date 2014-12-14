#ifndef TF_UTILS_BTREE_H
#define TF_UTILS_BTREE_H
#include "linkedlist.h"

struct btreenode
{
  struct btreenode *leftchild ;
  linkedlist* ids;
  double val;
  struct btreenode *rightchild ;
};

typedef struct btreenode btnode;

void bt_insert (btnode** bt, int id, double val);
void bt_delete (btnode** bt, int id, double val);
void bt_search (btnode** bt, int id, double val, 
    btnode **parent, btnode** x, int *found);
void bt_find_min (btnode* bt, btnode** x);
void bt_inorder (btnode* bt);
int bt_num_ids (btnode *bt);
void bt_free (btnode *bt);

#endif

