#include "linkedlist.h"
#include <stdio.h>
#include <stdlib.h>

llnode* create_node(int val)
{
  llnode* newnode;
  newnode = (llnode *)malloc(sizeof(llnode));
  if (newnode == NULL)
  {
    printf("\nMemory was not allocated");
    return 0;
  }
  else
  {
    newnode->value = val;
    newnode->next = NULL;
    return newnode;
  }
}

void insert_node(linkedlist** ll, int val)
{
  llnode* newnode;
  llnode* temp;
  newnode = create_node(val);
  if (*ll == NULL)
  {
    *ll = newnode;
  }
  else
  {
    temp = *ll;
    *ll = newnode;
    (*ll)->next = temp;
  }
}

void delete_node(linkedlist** ll, int key)
{
  llnode* ptr;
  llnode* prev;
  
  int found = 0;

  if (*ll == NULL)
  {    
    printf("Cannot delete from empty list\n");
    return;
  }
  for (ptr = *ll;ptr != NULL;ptr = ptr->next)
  {
    if (ptr->value == key)
    {
      found = 1;
      if(ptr == *ll)
        *ll = ptr->next;    
      else
        prev->next = ptr->next;      

      free(ptr);      
      break;
    }
    else
    {
      prev = ptr;
    }
  }
  if (found == 0)
  {
    printf("\nElement %d not found in list\n", key);
  }
}

int isempty(linkedlist* ll)
{
  return (ll == NULL);
}
int length(linkedlist* ll)
{
  int len = 0;
  while(ll != NULL)
  {
    len++;
    ll = ll->next;
  }
  return len;
}
void display(linkedlist* ll)
{
  llnode* ptr;
  if (ll == NULL)
  {
    printf("[]\n");
  }
  else
  {
    printf("[");
    for (ptr = ll;ptr != NULL;ptr = ptr->next)
    {    
      printf("%d ", ptr->value);
    }
    printf("]");
  }
}

void ll_free(linkedlist* ll)
{
  if(ll == NULL)
    return;

  ll_free(ll->next);
  free(ll);
}
