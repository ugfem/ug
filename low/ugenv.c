// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*	                                                                        */
/* File:      ugenv.c                                                        */
/*                                                                            */
/* Purpose:   implements the ug environment mechanism                        */
/*                                                                            */
/* Author:      Peter Bastian                                                 */
/*              Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen    */
/*              Universitaet Heidelberg                                        */
/*              Im Neuenheimer Feld 368                                        */
/*              6900 Heidelberg                                                */
/*                                                                            */
/* History:   18.02.92 begin, ug version 2.0                                */
/*                                                                            */
/* Revision:  08.09.95                                                      */
/*                                                                            */
/****************************************************************************/

/****************************************************************************/
/*                                                                            */
/* include files                                                            */
/*              system include files                                            */
/*              application include files                                     */
/*                                                                            */
/****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "compiler.h"
#include "heaps.h"
#include "general.h"
#include "misc.h"
#include "ugenv.h"

/****************************************************************************/
/*                                                                            */
/* defines in the following order                                            */
/*                                                                            */
/*          compile time constants defining static data size (i.e. arrays)    */
/*          other constants                                                    */
/*          macros                                                            */
/*                                                                            */
/****************************************************************************/

#define MAXENVPATH        32                /* max depth of the environment tree*/

/****************************************************************************/
/*                                                                            */
/* definition of variables global to this source file only (static!)        */
/*                                                                            */
/****************************************************************************/

static HEAP *envHeap=NULL;                /* heap used for the environment    */
static ENVDIR *path[MAXENVPATH];        /* path to current directory        */
static int pathIndex;                    /* entry to path array                */

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*D
   InitUgEnv - Initialize the Environment and the heap

   SYNOPSIS:
   INT InitUgEnv (INT heapSize);

   PARAMETERS:
   .  heapSize - size of the heap in bytes

   DESCRIPTION:
   This function initializes the Environment and the heap.

   RETURN VALUE:
   INT
   .n    0 if OK
   .n    '__LINE__' if not enough memory, error in initializing or allocating a
      heap structure
   D*/
/****************************************************************************/

INT InitUgEnv (INT heapSize)
{
  void *buffer;
  ENVDIR *root;

  /* allocate memory from system */
  if ((buffer=malloc(heapSize))==NULL) return(__LINE__);

  /* initialize heap structure */
  if ((envHeap=NewHeap(GENERAL_HEAP,heapSize,buffer))==NULL) return(__LINE__);

  /* allocate root directory */
  if ((root=GetMem(envHeap,sizeof(ENVDIR),0))==NULL) return(__LINE__);
  root->type = ROOT_DIR;
  root->next = root->previous = root->down = NULL;
  strcpy(root->name,"root");

  /* set path[0] */
  pathIndex = 0;
  path[0] = root;

  /* return ok */
  return(0);
}


/****************************************************************************/
/*D
   ChangeEnvDir    - Change environment directory

   SYNOPSIS:
   ENVDIR *ChangeEnvDir (const char *s);

   PARAMETERS:
   .  s - pointer to char (const)

   DESCRIPTION:
   This function changes environment directory.

   RETURN VALUE:
   ENVDIR *
   .n      pointer to new environment directory
   .n      NULL if directory not found.
   D*/
/****************************************************************************/

ENVDIR *ChangeEnvDir (const char *s)
{
  ENVDIR *newPath[MAXENVPATH];
  ENVDIR *theDir;
  INT i,k,len;
  char token[NAMESIZE];

  if (s==NULL) return (NULL);

  /* avoid trivial case */
  if ((len=strlen(s))==0) return(NULL);

  /* check max possible length */
  if (len>=MAXENVPATH*NAMESIZE)
    return (NULL);

  /* look at first character */
  if (strncmp(s,DIRSEP,1)==0)
  {
    /* start from root directory */
    newPath[0] = path[0];
    i = 0;
  }
  else
  {
    /* start from current directory */
    for (k=0; k<=pathIndex; k++) newPath[k] = path[k];
    i = pathIndex;
  }

  /* loop through input string */
  do
  {
    if ((s = strntok(s,DIRSEP,NAMELEN,token))==NULL)
      return (NULL);
    if (*token=='\0')
      continue;

    /* process token */
    if (strcmp(token,"..")==0)
    {
      /* one level up */
      if (i>0) i--;
    }
    else
    {
      /* find subdirectory */
      if (i+1>=MAXENVPATH) return(NULL);                  /* path too long    */
      theDir = (ENVDIR *) newPath[i]->down;              /* search next level*/
      while (theDir!=NULL)
      {
        if (theDir->type%2==1)                            /* is type odd ?    */
        {
          if (strcmp(token,theDir->name)==0)              /* name equal ?     */
            break;                                        /* directory found    */
        }
        theDir = (ENVDIR *) theDir->next;
      }
      if (theDir==NULL) return(NULL);                   /* not found error    */
      newPath[++i] = theDir;                              /* extend path        */
    }
  }
  while (*s!='\0');

  /* path found, copy to current path */
  for (k=0; k<=i; k++) path[k] = newPath[k];              /* copy path        */
  pathIndex = i;

  /* return current directory */
  return(path[pathIndex]);
}


/****************************************************************************/
/*D
   GetCurrentDir - Get current environment directory

   SYNOPSIS:
   ENVDIR *GetCurrentDir ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function gets current environment directory.

   RETURN VALUE:
   ENVDIR *
   .n       pointer to current environment directory
   D*/
/****************************************************************************/

ENVDIR *GetCurrentDir ()
{
  return(path[pathIndex]);
}

/****************************************************************************/
/*D
   GetPathName - Assemble pathname of current directory

   SYNOPSIS:
   void GetPathName (char *s);

   PARAMETERS:
   .  s - pointer to buffer for the string

   DESCRIPTION:
   This function assembles pathname of current directory.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void GetPathName (char *s)
{
  int i;

  strcpy(s,DIRSEP);
  for (i=1; i<=pathIndex; i++)
  {
    strcat(s,path[i]->name); strcat(s,DIRSEP);
  }
}

/****************************************************************************/
/*D
   MakeEnvItem - Allocate a new environment item in the current directory

   SYNOPSIS:
   ENVITEM *MakeEnvItem (const char *name, const INT type, const INT size);

   PARAMETERS:
   .  name - name of the new item
   .  type - type of item (5 types possible)
   .  size - size of item in bytes

   DESCRIPTION:
   This function allocates a new environment item in the current directory.

   RETURN VALUE:
   ENVITEM *
   .n      pointer to new item
   .n      NULL if name too long or already in use, type equal to 'ROOT_DIR',
             not enough memory
   .n           in the environment heap or other error occurs.
   D*/
/****************************************************************************/

ENVITEM *MakeEnvItem (const char *name, const INT type, const INT size)
{
  ENVITEM *newItem,*anItem,*lastItem;
  ENVDIR *currentDir;

  /* check if name too long */
  if (strlen(name)+1>NAMESIZE) return(NULL);

  /* check if name not already used in this directory */
  currentDir = path[pathIndex];
  anItem = lastItem = currentDir->down;
  while (anItem!=NULL)
  {
    if ((anItem->v.type==type)&&(strcmp(anItem->v.name,name)==0))
      return(NULL);
    lastItem = anItem;
    anItem = anItem->v.next;
  }

  /* allocate memory from environment heap */
  switch (type)
  {
  case ROOT_DIR :
    return(NULL);

  default :
    if (type%2==0)
    {
      /* new variable */
      newItem=(ENVITEM *) GetMem(envHeap,size,0);
      if (newItem==NULL) return(NULL);
      memset(newItem,0,size);
    }
    else
    {
      /* new directory */
      if (pathIndex+1>=MAXENVPATH) return(NULL);
      newItem=(ENVITEM *) GetMem(envHeap,size,0);
      if (newItem==NULL) return(NULL);
      memset(newItem,0,size);
      newItem->d.down = NULL;
    }
    break;
  }

  /* initialize item */
  newItem->v.type = type;
  newItem->v.locked = 1;
  strcpy(newItem->v.name,name);

  /* put in list */
  if (lastItem==NULL)
  {
    /* directory is empty */
    currentDir->down = newItem;
    newItem->v.next = newItem->v.previous = NULL;
  }
  else
  {
    /* append to last Item */
    lastItem->v.next = newItem;
    newItem->v.previous = lastItem;
    newItem->v.next = NULL;
  }

  /* return pointer to new item */
  return(newItem);
}


/****************************************************************************/
/*D
   RemoveEnvItem - Deallocate an environment item in the current directory

   SYNOPSIS:
   INT RemoveEnvItem (ENVITEM *theItem);

   PARAMETERS:
   .  theItem - pointer to item

   DESCRIPTION:
   This function deallocates an environment item in the current directory.
   Directories may only be deleted if they are empty.

   RETURN VALUE:
   INT
   .n    0 if OK
   .n    1 if item not found in current directory
   .n    2 if attempt is done to delete non empty directory
   .n    3 if attempt is done to delete locked item.
   D*/
/****************************************************************************/

INT RemoveEnvItem (ENVITEM *theItem)
{
  ENVITEM *anItem;
  ENVDIR *currentDir;

  /* check if item is in current directory */
  currentDir = path[pathIndex];
  anItem = currentDir->down;
  while (anItem!=NULL)
  {
    if (anItem==theItem) break;
    anItem = anItem->v.next;
  }
  if (anItem==NULL) return(1);

  /* check if item is allowed to be deleted */
  if (theItem->v.locked) return(3);
  if ((theItem->v.type%2==1)&&(theItem->d.down!=NULL)) return(2);

  /* remove item from double linked list */
  if (theItem->v.previous==NULL)
    currentDir->down = theItem->v.next;
  else
    theItem->v.previous->v.next = theItem->v.next;
  if (theItem->v.next!=NULL)
    theItem->v.next->v.previous = theItem->v.previous;

  /* deallocate memory */
  DisposeMem(envHeap,theItem);

  /* return ok */
  return(0);
}

/****************************************************************************/
/*D
   MoveEnvItem - move an environment item in the tree structure

   SYNOPSIS:
   INT MoveEnvItem (ENVITEM *item, ENVDIR *old, ENVDIR *new)

   PARAMETERS:
   .  item - pointer to item
   .  old - directory containing the item
   .  new - directory to which item has to be moved (if NULL move to root)

   DESCRIPTION:
   This function moves an environment item (or subtree) in the tree structure.

   RETURN VALUE:
   INT
   .n    0 if OK
   D*/
/****************************************************************************/

INT MoveEnvItem (ENVITEM *item, ENVDIR *old, ENVDIR *new)
{
  ENVITEM *anItem;

  if (new==NULL)
    new = path[0];

  for (anItem=ENVDIR_DOWN(old); anItem!=NULL; anItem=NEXT_ENVITEM(anItem))
    if (anItem==item) break;
  if (anItem==NULL) return(1);

  /* remove from old directory */
  if (PREV_ENVITEM(item)!=NULL)
    NEXT_ENVITEM(PREV_ENVITEM(item)) = NEXT_ENVITEM(item);
  else
    ENVDIR_DOWN(old) = NEXT_ENVITEM(item);
  if (NEXT_ENVITEM(item)!=NULL)
    PREV_ENVITEM(NEXT_ENVITEM(item)) = PREV_ENVITEM(item);

  /* insert in new directory */
  PREV_ENVITEM(item) = NULL;
  NEXT_ENVITEM(item) = ENVDIR_DOWN(new);
  ENVDIR_DOWN(new) = item;

  return (0);
}

/****************************************************************************/
/*D
   SearchTree - Search for a given name in the tree

   SYNOPSIS:
   static ENVITEM *SearchTree (const char *name, INT type, INT dirtype);

   PARAMETERS:
   .  name - name of item searched for
   .  type - one of the five basic types
   .  dirtype - directory type to scan or SEARCHALL

   DESCRIPTION:
   This function searches a given name in the tree, sets current directory
   to directory of the item found and returns the item.

   RETURN VALUE:
   ENVITEM *
   .n      pointer to the item in the tree
   .n      NULL if name not found or where not OK.
   D*/
/****************************************************************************/

static ENVITEM *SearchTree (const char *name, INT type, INT dirtype)
{
  ENVDIR *currentDir;
  ENVITEM *theItem,*result;

  currentDir = path[pathIndex];

  /* first loop in current directory */
  theItem = currentDir->down;
  while (theItem!=NULL)
  {
    if (theItem->v.type == type)
      if (strcmp(theItem->v.name,name)==0)
        return(theItem);
    theItem = theItem->v.next;
  }

  /* try recursive search */
  theItem = currentDir->down;
  while (theItem!=NULL)
  {
    if (theItem->v.type%2 == 1)
    {
      /* this is a directory */
      if ((theItem->d.type==dirtype)||(dirtype==SEARCHALL))
      {
        path[++pathIndex] = (ENVDIR *) theItem;
        if ((result=SearchTree(name,type,dirtype))!=NULL)
          return(result);
        pathIndex--;
      }
    }
    theItem = theItem->v.next;
  }

  /* return not found */
  return(NULL);
}

/****************************************************************************/
/*D
   SearchEnv - Search for a given name in the tree

   SYNOPSIS:
   ENVITEM *SearchEnv (const char *name, const char *where, INT type, INT dirtype);

   PARAMETERS:
   .  name - name of item searched for
   .  type - one of the five basic types
   .  where - path to directory where recursive search starts
   .  dirtype - directory type to scan or SEARCHALL

   DESCRIPTION:
   This function searches a given name in the tree, sets current directory
   to directory of the item found and returns the item.

   RETURN VALUE:
   ENVITEM *
   .n     pointer to item in the tree
   .n     NULL if name not found or where not OK.
   D*/
/****************************************************************************/

ENVITEM *SearchEnv (const char *name, const char *where, INT type, INT dirtype)
{
  /* check if search directory is changed */
  if (strcmp(where,".")!=0)
    if (ChangeEnvDir(where)==NULL) return(NULL);

  /* recursive search */
  return(SearchTree(name,type,dirtype));
}


/****************************************************************************/
/*D
   AllocEnvMemory - Allocate memory from environment heap

   SYNOPSIS:
   void *AllocEnvMemory (INT size);

   PARAMETERS:
   .  size - number of bytes to be allocated

   DESCRIPTION:
   This function allocates memory from environment heap.

   RETURN VALUE:
   void *
   .n       pointer to allocated memory
   D*/
/****************************************************************************/

void *AllocEnvMemory (INT size)
{
  return(GetMem(envHeap,size,0));
}

/****************************************************************************/
/*D
   FreeEnvMemory - Deallocate memory from environment heap

   SYNOPSIS:
   void FreeEnvMemory (void *buffer);

   PARAMETERS:
   .  buffer - pointer to buffer previously allocated

   DESCRIPTION:
   This function deallocates memory from environment heap.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void FreeEnvMemory (void *buffer)
{
  DisposeMem(envHeap,buffer);
}

/****************************************************************************/
/*D
   EnvHeapInfo - Print size and used of environment heap to string

   SYNOPSIS:
   void EnvHeapInfo (char *s);

   PARAMETERS:
   .  s - string to print on

   DESCRIPTION:
   This function prints size and used of environment heap to string.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void EnvHeapInfo (char *s)
{
  sprintf(s,"   size: %ld\n   used: %ld\n",HeapSize(envHeap),HeapUsed(envHeap));
}

/****************************************************************************/
/*D
   GetNewEnvDirID - Return a unique 'ID' for a new 'ENVDIR' type

   SYNOPSIS:
   INT GetNewEnvDirID ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function returns a unique 'ID' for a new 'ENVDIR' type.

   RETURN VALUE:
   INT
   .n    'ID' for the new 'ENVDIR'
   D*/
/****************************************************************************/

INT GetNewEnvDirID ()
{
  /* NB: ENVDIRs have odd types and start with 1 */
  static theNewEnvDirID = 1;

  theNewEnvDirID += 2;

  return (theNewEnvDirID);
}

/****************************************************************************/
/*D
   GetNewEnvVarID - Return a unique 'ID' for a new 'ENVVAR' type

   SYNOPSIS:
   INT GetNewEnvVarID ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function returns a unique 'ID' for a new 'ENVVAR' type.

   RETURN VALUE:
   INT
   .n    'ID' for the new 'ENVVAR'
   D*/
/****************************************************************************/

INT GetNewEnvVarID ()
{
  /* NB: ENVVARs have even types and start with 2 */
  static theNewEnvVarID = 0;

  theNewEnvVarID += 2;

  return (theNewEnvVarID);
}
