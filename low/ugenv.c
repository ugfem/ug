// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ugenv.c                                                       */
/*                                                                          */
/* Purpose:   implements the ug environment mechanism                       */
/*                                                                          */
/* Author:      Peter Bastian                                               */
/*              Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen */
/*              Universitaet Heidelberg                                     */
/*              Im Neuenheimer Feld 368                                     */
/*              6900 Heidelberg                                             */
/*                                                                          */
/* History:   18.02.92 begin, ug version 2.0                                */
/*                                                                          */
/* Revision:  08.09.95                                                      */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*              system include files                                        */
/*              application include files                                   */
/*                                                                          */
/****************************************************************************/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "ugtypes.h"
#include "heaps.h"
#include "general.h"
#include "misc.h"
#include "ugenv.h"
/** \todo this is a hierarchy conflict, remove. (VR) */
#include "ugdevices.h"

USING_UG_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*          compile time constants defining static data size (i.e. arrays)  */
/*          other constants                                                 */
/*          macros                                                          */
/*                                                                          */
/****************************************************************************/

#define MAXENVPATH        32            /* max depth of the environment tree*/

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/** \brief Path to current directory

   We only need the first entry to be zero-initialized. */
static ENVDIR *path[MAXENVPATH] = {NULL};
static int pathIndex;                   /* entry to path array              */

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/** \brief Initialize the Environment and the heap

 * @param   heapSize - size of the heap in bytes

   This function initializes the Environment and the heap.

   @return <ul>
 *   <li>     0 if OK
 *   <li>     '__LINE__' if not enough memory, error in initializing or allocating a
      heap structure
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX InitUgEnv ()
{
  ENVDIR *root;

  /* Environment heap already initialized? */
  if (path[0])
    return 0;

  /* allocate root directory */
  if ((root=(ENVDIR*)malloc(sizeof(ENVDIR)))==NULL) return(__LINE__);

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
/** \brief Change environment directory

   This function changes environment directory.

   @return <ul>
 *   <li>       pointer to new environment directory </li>
 *   <li>       NULL if directory not found. </li>
   </ul>
 */
/****************************************************************************/

ENVDIR * NS_PREFIX ChangeEnvDir (const char *s)
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
/** \brief GetCurrentDir - Get current environment directory

   @return <ul>
   <li>        pointer to current environment directory
   </ul>
 */
/****************************************************************************/

ENVDIR * NS_PREFIX GetCurrentDir ()
{
  return(path[pathIndex]);
}

/****************************************************************************/
/** \brief Assemble pathname of current directory

   @param   s - pointer to buffer for the string

 */
/****************************************************************************/

void NS_PREFIX  GetPathName (char *s)
{
  int i;

  strcpy(s,DIRSEP);
  for (i=1; i<=pathIndex; i++)
  {
    strcat(s,path[i]->name); strcat(s,DIRSEP);
  }
}

/****************************************************************************/
/** \brief Allocate a new environment item in the current directory
 *
 * @param   name - name of the new item
 * @param   type - type of item (5 types possible)
 * @param   size - size of item in bytes
 *
 * This function allocates a new environment item in the current directory.
 *
 *  @return <ul>
 *   <li>       pointer to new item
 *   <li>       NULL if name too long or already in use, type equal to 'ROOT_DIR',
             not enough memory
 *   <li>            in the environment heap or other error occurs.
   </ul>
 */
/****************************************************************************/

ENVITEM * NS_PREFIX MakeEnvItem (const char *name, const INT type, const INT size)
{
  ENVITEM *newItem,*anItem,*lastItem;
  ENVDIR *currentDir;

  /* check if name too long */
  if (strlen(name)+1>NAMESIZE)
    return(NULL);

  /* check if name not already used in this directory */
  currentDir = path[pathIndex];
  anItem = lastItem = currentDir->down;

#ifdef Debug
  while (anItem!=NULL)
  {
    if ((anItem->v.type==type)&&(strcmp(anItem->v.name,name)==0))
    {
      UserWriteF("MakeEnvItem(): item '%s' already defined\n",anItem->v.name);
      return(NULL);
    }
    lastItem = anItem;
    anItem = anItem->v.next;
  }
#endif

  /* allocate memory from environment heap */
  switch (type)
  {
  case ROOT_DIR :
    return(NULL);

  default :
    if (type%2==0)
    {
      /* new variable */
      newItem=(ENVITEM*) malloc(size);

      if (newItem==NULL)
      {
        UserWriteF("MakeEnvItem(): envHeap out of memory\n");
        return(NULL);
      }
      memset(newItem,0,size);
    }
    else
    {
      /* new directory */
      if (pathIndex+1>=MAXENVPATH) return(NULL);
      newItem=(ENVITEM*) malloc(size);

      if (newItem==NULL)
      {
        UserWriteF("MakeEnvItem(): envHeap out of memory\n");
        return(NULL);
      }
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
    /* insert before first item */
    newItem->v.previous = NULL;
    currentDir->down->v.previous = newItem;
    newItem->v.next = currentDir->down;
    currentDir->down = newItem;
  }

  /* return pointer to new item */
  return(newItem);
}


/****************************************************************************/
/** \brief Deallocate an environment item in the current directory

   @param   theItem - pointer to item

   This function deallocates an environment item in the current directory.
   Directories may only be deleted if they are empty.

   @return <ul>
 *   <li>     0 if OK
 *   <li>     1 if item not found in current directory
 *   <li>     2 if attempt is done to delete non empty directory
 *   <li>     3 if attempt is done to delete locked item.
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX RemoveEnvItem (ENVITEM *theItem)
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
  free(theItem);

  /* return ok */
  return(0);
}

#ifndef __T3E__
static
#endif
INT RemoveEnvDirContent (ENVITEM *theItem)
{
  ENVITEM *Item,*Next;

  for (Item = theItem; Item != NULL; Item = Next) {
    Next = NEXT_ENVITEM(Item);
    if (IS_ENVDIR(Item))
      RemoveEnvDirContent(ENVITEM_DOWN(Item));
    free(Item);
  }

  return(0);
}

/****************************************************************************/
/** \brief Deallocate an environment directory

 * @param   theItem - pointer to item

   This function deallocates an environment directory.

   @return <ul>
   <li>     0 if OK
   <li>     1 if item not found in current directory
   <li>     2 if theItem is not a directory
   <li>     3 if attempt is done to delete locked item
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX RemoveEnvDir (ENVITEM *theItem)
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
  if (!IS_ENVDIR(theItem)) return(2);
  if (theItem->v.locked) return(3);
  RemoveEnvDirContent(ENVITEM_DOWN(theItem));

  /* remove item from double linked list */
  if (theItem->v.previous==NULL)
    currentDir->down = theItem->v.next;
  else
    theItem->v.previous->v.next = theItem->v.next;
  if (theItem->v.next!=NULL)
    theItem->v.next->v.previous = theItem->v.previous;

  /* deallocate memory */
  free(theItem);

  return(0);
}

/****************************************************************************/
/** \brief Move an environment item in the tree structure

 * @param   item - pointer to item
 * @param   old - directory containing the item
 * @param   new - directory to which item has to be moved (if NULL move to root)

   This function moves an environment item (or subtree) in the tree structure.

   @return <ul>
 *   <li>     0 if OK
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX MoveEnvItem (ENVITEM *item, ENVDIR *oldDir, ENVDIR *newDir)
{
  ENVITEM *anItem;

  if (newDir==NULL)
    newDir = path[0];

  for (anItem=ENVDIR_DOWN(oldDir); anItem!=NULL; anItem=NEXT_ENVITEM(anItem))
    if (anItem==item) break;
  if (anItem==NULL) return(1);

  /* remove from old directory */
  if (PREV_ENVITEM(item)!=NULL)
    NEXT_ENVITEM(PREV_ENVITEM(item)) = NEXT_ENVITEM(item);
  else
    ENVDIR_DOWN(oldDir) = NEXT_ENVITEM(item);
  if (NEXT_ENVITEM(item)!=NULL)
    PREV_ENVITEM(NEXT_ENVITEM(item)) = PREV_ENVITEM(item);

  /* insert in new directory */
  PREV_ENVITEM(item) = NULL;
  NEXT_ENVITEM(item) = ENVDIR_DOWN(newDir);
  ENVDIR_DOWN(newDir) = item;

  return (0);
}

/****************************************************************************/
/** \brief Search for a given name in the tree

 * @param   name - name of item searched for
 * @param   type - one of the five basic types
 * @param   dirtype - directory type to scan or SEARCHALL

   This function searches a given name in the tree, sets current directory
   to directory of the item found and returns the item.

   @return <ul>
 *   <li>       pointer to the item in the tree
 *   <li>       NULL if name not found or where not OK.
   </ul>
 */
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
/** \brief  Search for a given name in the tree
 *
 * @param   name - name of item searched for
 * @param   type - one of the five basic types
 * @param   where - path to directory where recursive search starts
 * @param   dirtype - directory type to scan or SEARCHALL
 *
 * This function searches a given name in the tree, sets current directory
 * to directory of the item found and returns the item.
 *
 * @return <ul>
   <li> pointer to item in the tree
   <li> NULL if name not found or where not OK.
   </ul>
 */
/****************************************************************************/

ENVITEM * NS_PREFIX SearchEnv (const char *name, const char *where, INT type, INT dirtype)
{
  /* check if search directory is changed */
  if (strcmp(where,".")!=0)
    if (ChangeEnvDir(where)==NULL) return(NULL);

  /* recursive search */
  return(SearchTree(name,type,dirtype));
}


/****************************************************************************/
/** \brief For backward compatibility: Allocate memory from environment heap

 * @param   size - number of bytes to be allocated

   \deprecated Simply a wrapper for malloc()

   @return <ul>
   Pointer to allocated memory
 */
/****************************************************************************/

void * NS_PREFIX AllocEnvMemory (INT size)
{
  return malloc(size);
}

/****************************************************************************/
/** \brief For backward compatibility: Deallocate memory from environment heap

 * @param   buffer - pointer to buffer previously allocated

   \deprecated Simply a wrapper for free()
 */
/****************************************************************************/

void NS_PREFIX FreeEnvMemory (void *buffer)
{
  free(buffer);
}

/****************************************************************************/
/** \brief For backward compatibility: Print size and used of environment heap to string

 * @param   s - string to print on

   This function used to print size and used percentage of the environment heap
   to a string.  Ever since the operating system heap is used this information
   is not available anymore.
 */
/****************************************************************************/

void NS_PREFIX EnvHeapInfo (char *s)
{
  sprintf(s,"no heap information available\n");
}

/****************************************************************************/
/** \brief Return a unique 'ID' for a new 'ENVDIR' type

   This function returns a unique 'ID' for a new 'ENVDIR' type.

   @return 'ID' for the new 'ENVDIR'
 */
/****************************************************************************/

INT NS_PREFIX GetNewEnvDirID (void)
{
  /* NB: ENVDIRs have odd types and start with 1 */
  static INT theNewEnvDirID = 1;

  theNewEnvDirID += 2;

  return (theNewEnvDirID);
}

/****************************************************************************/
/** \brief Return a unique 'ID' for a new 'ENVVAR' type

   This function returns a unique 'ID' for a new 'ENVVAR' type.

   @return 'ID' for the new 'ENVVAR'
 */
/****************************************************************************/

INT NS_PREFIX GetNewEnvVarID ()
{
  /* NB: ENVVARs have even types and start with 2 */
  static INT theNewEnvVarID = 0;

  theNewEnvVarID += 2;

  return (theNewEnvVarID);
}


INT NS_PREFIX ExitUgEnv()
{
  RemoveEnvDirContent((ENVITEM*)path[0]);
  path[0] = NULL;
  return 0;
}
