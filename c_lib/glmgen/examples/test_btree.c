#include "btree.h"
#include <stdio.h>
#include <stdlib.h>

static void test_insert_delete()
{
  struct btreenode *bt;
  int i;

  bt = NULL;

  double values[] = {50, 30, 70, 10, 30, 50};
  for(i=0; i < 6; i++) {
    printf("Inserting (%d,%g)\n", i, values[i]);
    bt_insert(&bt, i, values[i]);
    bt_inorder(bt);
    printf("\n");
  }

  int ids[] = {1, 1, 0, 3, 2, 2, 5, 1, 4, 0 }; 
  double del_values[] = {100, 30, 50, 10, 70, 70, 50, 30, 30, 50};
  for(i=0; i < 10; i++) {
    printf("Deleting (%d,%g)\n", ids[i], del_values[i]);
    bt_delete(&bt, ids[i], del_values[i]);
    bt_inorder(bt);
    printf("\n\n");
  }
  bt_free(bt);

}
int main()
{
  test_insert_delete();
  
  return 0;
}
