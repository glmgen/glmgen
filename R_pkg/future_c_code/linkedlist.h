#ifndef TF_UTILS_LINKEDLIST_H
#define TF_UTILS_LINKEDLIST_H

struct llnode
{
  int value;
  struct llnode *next;
};

typedef struct llnode llnode;
typedef llnode linkedlist;

llnode* create_node(int val);
void insert_node(linkedlist** ll, int val);
void delete_node(linkedlist** ll, int key);
int isempty(linkedlist* ll);
int length(linkedlist* ll);
void display(linkedlist* ll);
void ll_free(linkedlist* ll);
/*
llnode* create_node(int);
void insert_node_first();
void insert_node_last();
void insert_node_pos();
void sorted_ascend();
void delete_pos();
void search();
void update_val();
void display();
void rev_display(llnode *);
*/

#endif
