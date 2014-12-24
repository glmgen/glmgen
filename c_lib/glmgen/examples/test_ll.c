#include "linkedlist.h"
#include <stdio.h>

void test_insert()
{
  printf("Testing linkedlist insert_node\n");

  int i;
  
  linkedlist* ll = NULL;
  for(i=0; i < 5; i++) {
    printf("Inserting %d\n", i);
    insert_node(&ll, i);
    printf("List: "); display(ll); printf("\n");
    printf("Size: %d\n", length(ll)); 

    printf("\n");
    
  }

  ll_free(ll);
}

void test_delete()
{
  printf("Testing linkedlist delete_node\n");
  
  int i;
  linkedlist* ll = NULL;

  delete_node(&ll, 2);

  for(i=0; i < 5; i++) insert_node(&ll, i);

  int vals[] = { 6, 2, 1, 3, 0, 4, 1};
  for(i=0; i < 7; i++) {
    delete_node(&ll, vals[i]); 
    printf("List: "); display(ll); printf("\n");
    printf("Size: %d\n\n", length(ll)); 
  }
  insert_node(&ll, 3); 
  printf("List: "); display(ll); printf("\n");
  printf("Size: %d\n\n", length(ll)); 

  ll_free(ll);

}

int main()
{
  test_insert();
  test_delete();

  return 0;
}
