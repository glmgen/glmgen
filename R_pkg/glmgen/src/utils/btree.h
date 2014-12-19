#ifndef TF_UTILS_BTREE_H
#define TF_UTILS_BTREE_H

typedef enum { FIRST, SECOND } bt_node_lvl;

struct btreenode
{
  struct btreenode *leftchild;
  double key;
  struct btreenode *ids;
  bt_node_lvl node_lvl;
  struct btreenode *rightchild;
};

typedef struct btreenode btnode;

void bt_insert (btnode** bt, int id, double val);
void bt_insert_inner (btnode **bt, int id);
void bt_delete (btnode** bt, int id, double val);
void bt_delete_inner (btnode **bt, int id);
void bt_delete_found_node(btnode **bt, btnode **parent, btnode *x);
void bt_search (btnode **bt, double val,
    btnode **par, btnode **x, int *found);
void bt_find_min (btnode* bt, btnode** x);
void bt_find_min_twice (btnode* bt, btnode** x);
void bt_inorder (btnode* bt);
void bt_free (btnode *bt);

#endif

