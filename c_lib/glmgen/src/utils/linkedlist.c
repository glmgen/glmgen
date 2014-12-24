/****************************************************************************
 * Copyright (C) 2014 by Taylor Arnold, Ryan Tibshirani, Veerun Sadhanala   *
 *                                                                          *
 * This file is part of the glmgen library / package.                       *
 *                                                                          *
 *   glmgen is free software: you can redistribute it and/or modify it      *
 *   under the terms of the GNU Lesser General Public License as published  *
 *   by the Free Software Foundation, either version 2 of the License, or   *
 *   (at your option) any later version.                                    *
 *                                                                          *
 *   glmgen is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU Lesser General Public License for more details.                    *
 *                                                                          *
 *   You should have received a copy of the GNU Lesser General Public       *
 *   License along with glmgen. If not, see <http://www.gnu.org/licenses/>. *
 ****************************************************************************/

/**
 * @file linkedlist.c
 * @author Taylor Arnold, Ryan Tibshirani, Veerun Sadhanala
 * @date 2014-12-23
 * @brief Construct and interact with a linked list.
 */

#include "utils.h"
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
  prev = NULL;

  if (*ll == NULL)
  {
    printf("Cannot delete from empty list\n");
    return;
  }
  for (ptr = *ll; ptr != NULL; ptr = ptr->next)
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
