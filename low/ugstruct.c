// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugstruct.c													*/
/*																			*/
/* Purpose:   structure administration for ug 2.0							*/
/*																			*/
/* Author:	  Peter Bastian, Nicolas Neuss									*/
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  email: ug@ica3.uni-stuttgart.de			                        */
/*																			*/
/* History:   18.02.92 begin, ug version 2.0								*/
/*			  10.05.92 begin, ug interpreter by Nicolas Neuss				*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#ifdef __MPW32__
#pragma segment cmd
#endif

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "compiler.h"
#include "general.h"
#include "heaps.h"
#include "misc.h"
#include "ugenv.h"

#include "ugstruct.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define qErr -1

#define MAXENVPATH              32              /* max depth of the environment tree		*/
#define STRUCTSEP ":"
#define STRUCTSEPC ':'

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static ENVDIR *path[MAXENVPATH];
static INT pathIndex;

static INT theStringDirID;                      /* env type for String dirs                     */
static INT theStringVarID;                      /* env type for String vars                     */

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*D
   ChangeStructDir - change current structure directory

   SYNOPSIS:
   ENVDIR *ChangeStructDir (const char *name);

   PARAMETERS:
   .  name - name to which the current structure directory shall be changed

   DESCRIPTION:
   Changes the current structure directory to the one given by name. The name is
   interpreted as starting from the previous end of the path. It is possible to
   use ".." to go to the parent directory and ":" at the beginning
   to start from the root directory.

   RETURN VALUE:
   ENVDIR *
   .n The new directory or NULL, if it was not found.
   D*/
/****************************************************************************/

ENVDIR *ChangeStructDir (const char *name)
{
  ENVDIR *newPath[MAXENVPATH];
  ENVDIR *theDir;
  char token[NAMESIZE];
  const char *s;
  size_t len;
  INT i,k;

  if (name==NULL) return(NULL);

  /* avoid trivial case */
  if ((len=strlen(name))==0) return(NULL);

  /* check max possible length */
  if (len>=MAXENVPATH*NAMESIZE)
    return (NULL);

  /* look at first character */
  if (*name==STRUCTSEPC)
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
  s = name;
  while (*s!='\0')
  {
    if ((s = strntok(s,STRUCTSEP,NAMELEN,token))==NULL)
      return (NULL);
    if (*token=='\0')
      break;

    /* process token (find dir) */
    if (strcmp(token,"..")==0)
    {
      /* one level up */
      if (i>0) i--;
    }
    else
    {
      /* find subdirectory */
      if (i+1>=MAXENVPATH)
      {
        return(NULL);                                                                           /* path too long	*/
      }
      theDir = (ENVDIR *) newPath[i]->down;                             /* search next level*/
      while (theDir!=NULL)
      {
        if (theDir->type%2==1)                                                          /* is type odd ?	*/
        {
          if (strcmp(token,theDir->name)==0)                                    /* name equal ?         */
            break;                                                                                      /* directory found	*/
        }
        theDir = (ENVDIR *) theDir->next;
      }
      if (theDir==NULL)
        return(NULL);                           /* not found */

      newPath[++i] = theDir;                                                            /* extend path		*/
    }
  }

  /* path found, copy to current path */
  for (k=0; k<=i; k++) path[k] = newPath[k];                            /* copy path		*/
  pathIndex = i;

  /* return current directory */
  return(path[pathIndex]);
}

/****************************************************************************/
/*D
   FindStructDir - finds the directory with the given name or its parent

   SYNOPSIS:
   ENVDIR *FindStructDir (const char *name, char **lastnameHnd);

   PARAMETERS:
   .  name - name to be found (from current structure directory)
   .  lastnameHnd - if not NULL it is set to a pointer to the last name

   DESCRIPTION:
   This function searchs (starting from the current structure directory)
   for the structure directory given by name and returns its adress.
   It does not change the current structure directory! If 'lastnameHnd!=NULL'
   it also sets '*lastNameHnd' to be a pointer to the last name of the directory.

   SEE ALSO:
   strntok

   RETURN VALUE:
   ENVDIR *
   .n The directory, or NULL if not found.
   D*/
/****************************************************************************/

static char token[NAMESIZE];                    /* must be global for lastnameHnd ! */
static char nexttoken[NAMESIZE];        /* must be global for lastnameHnd ! */

ENVDIR *FindStructDir (const char *name, char **lastnameHnd)
{
  ENVDIR *newPath[MAXENVPATH];
  ENVDIR *theDir;
  const char *s;
  size_t len;
  INT i,k;

  if (name==NULL) return(NULL);

  /* avoid trivial case */
  if ((len=strlen(name))==0) return(NULL);

  /* check max possible length */
  if (len>=MAXENVPATH*NAMESIZE)
    return (NULL);

  /* look at first character */
  if (*name==STRUCTSEPC)
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

  if ((s = strntok(name,STRUCTSEP,NAMELEN,token))==NULL)
    return (NULL);

  if (*s=='\0')
  {
    if (lastnameHnd!=NULL)
      *lastnameHnd=token;

    return(newPath[i]);
  }

  /* loop through input string */
  while (*s!='\0')
  {
    /* process token (find dir) */
    if (strcmp(token,"..")==0)
    {
      /* one level up */
      if (i>0) i--;
    }
    else
    {
      /* find subdirectory */
      if (i+1>=MAXENVPATH)
      {
        return(NULL);                                                                           /* path too long	*/
      }
      theDir = (ENVDIR *) newPath[i]->down;                                     /* search next level*/
      while (theDir!=NULL)
      {
        if (theDir->type%2==1)                                                          /* is type odd ?	*/
        {
          if (strcmp(token,theDir->name)==0)                                    /* name equal ?         */
            break;                                                                                      /* directory found	*/
        }
        theDir = (ENVDIR *) theDir->next;
      }
      if (theDir==NULL)
        return(NULL);                           /* not found */

      newPath[++i] = theDir;                                                            /* extend path		*/
    }

    if ((s = strntok(s,STRUCTSEP,NAMELEN,nexttoken))==NULL)
      return (NULL);

    if (*nexttoken=='\0')
      break;

    if (lastnameHnd!=NULL)
      if (*s!=STRUCTSEPC) break;

    /* copy next token */
    strcpy(token,nexttoken);
  }

  if (lastnameHnd!=NULL)
    *lastnameHnd=nexttoken;

  return(newPath[i]);
}

/****************************************************************************/
/*D
   FindStringVar - searches a string variable inside a directory/structure

   SYNOPSIS:
   STRVAR *FindStringVar (const ENVDIR *where, const char *name);

   PARAMETERS:
   .  where - the directory to be searched
   .  name - the name to be searched

   DESCRIPTION:
   See above.

   RETURN VALUE:
   STRVAR *
   .n The address of the variable or NULL, if not found.
   D*/
/****************************************************************************/

STRVAR *FindStringVar (const ENVDIR *where, const char *name)
{
  STRVAR *myVar;

  /* find variable */
  myVar = (STRVAR *) where->down;       /* search directory */
  while (myVar!=NULL)
  {
    if (myVar->v.type==theStringVarID)
    {
      if (strcmp(name,myVar->v.name)==0)
        return(myVar);                          /* variable found  */
    }
    myVar = (STRVAR *) myVar->v.next;
  }

  return(NULL);         /* not found */
}

/****************************************************************************/
/*D
   FindStructure - searches a directory/structure inside a directory/structure

   SYNOPSIS:
   ENVDIR *FindStructure (const ENVDIR *where, const char *name)

   PARAMETERS:
   .  where - the directory to be searched
   .  name - the name to be searched

   RETURN VALUE:
   ENVDIR *
   .n The address of the structure or NULL, if not found.
   D*/
/****************************************************************************/

ENVDIR *FindStructure (const ENVDIR *where, const char *name)
{
  ENVDIR *theDir;

  /* find variable */
  theDir = (ENVDIR *) where->down;              /* search directory */
  while (theDir!=NULL)
  {
    if (theDir->type==theStringDirID)
    {
      if (strcmp(name,theDir->name)==0)
        return(theDir);                         /* structure found	*/
    }
    theDir = (ENVDIR *) theDir->next;
  }

  return(NULL);         /* not found */
}

/****************************************************************************/
/*D
    GetStringVar - Get string contained in the variable 'name'

    SYNOPSIS:
    char *GetStringVar (const char *name);

    PARAMETERS:
   .   name - name of the string variable (from the current structure directory)

    DESCRIPTION:
    This function returns the string contained in the variable 'name'.

    RETURN VALUE:
    char *
   .n  pointer to the string
   .n  NULL if not found
   D*/
/****************************************************************************/

char *GetStringVar (const char *name)
{
  char *lastname;
  ENVDIR *theDir;
  STRVAR *myVar;

  if ((theDir=FindStructDir(name,&lastname))==NULL)
    return(NULL);               /* structure directory not found */

  if ((myVar=FindStringVar(theDir,lastname))!=NULL)
    return(myVar->s);
  else
    return(NULL);
}

/****************************************************************************/
/*
   GetStringValue - Get double value of the string contained in the variable 'name'

   SYNOPSIS:
   INT GetStringValue (const char *name, double *value);

   PARAMETERS:
   .  name - name of the string variable
   .  value - address where the result is stored

   DESCRIPTION:
   Same as GetStringValueDouble. Should not be used and may be skipped
   in a future version.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
 */
/****************************************************************************/

INT GetStringValue (const char *name, double *value)
{
  char *lastname;
  ENVDIR *theDir;
  STRVAR *myVar;

  if ((theDir=FindStructDir(name,&lastname))==NULL)
    return(1);                  /* structure directory not found */

  if ((myVar=FindStringVar(theDir,lastname))==NULL)
    return(1);

  if (sscanf(myVar->s,"%lf",value)!=1)
    return (1);

  return (0);
}

/****************************************************************************/
/*D
   GetStringValueDouble	- Get double value of a variable

   SYNOPSIS:
   INT GetStringValueDouble (const char *name, double *value);

   PARAMETERS:
   .  name - name of the string variable
   .  value - address where the result is stored

   DESCRIPTION:
   This function evaluates a string variable holding a floating point number
   and returns the result.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT GetStringValueDouble (const char *name, double *value)
{
  char *lastname;
  ENVDIR *theDir;
  STRVAR *myVar;
  double val;

  if ((theDir=FindStructDir(name,&lastname))==NULL)
    return(1);                  /* structure directory not found */

  if ((myVar=FindStringVar(theDir,lastname))==NULL)
    return(1);

  if (sscanf(myVar->s,"%lf",&val)!=1)
    return (1);

  *value = val;
  return (0);
}

/****************************************************************************/
/*D
   GetStringValueInt - Get integer value of a variable

   SYNOPSIS:
   INT GetStringValueInt (const char *name, int *value);

   PARAMETERS:
   .  name - name of the string variable
   .  value - place result here

   DESCRIPTION:
   This function evaluates a string variable holding an integer number
   and returns the result.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT GetStringValueInt (const char *name, int *value)
{
  char *lastname;
  ENVDIR *theDir;
  STRVAR *myVar;
  int val;

  if ((theDir=FindStructDir(name,&lastname))==NULL)
    return(1);                  /* structure directory not found */

  if ((myVar=FindStringVar(theDir,lastname))==NULL)
    return(1);

  if (sscanf(myVar->s,"%d",&val)!=1)
    return (1);

  *value = val;
  return (0);
}

/****************************************************************************/
/*D
   GetStringDOUBLEInRange - Get the double value of 'name' together with a range check

   SYNOPSIS:
   INT GetStringDOUBLEInRange (const char *name, DOUBLE min,
   DOUBLE max, DOUBLE *value);

   PARAMETERS:
   .  name - memory address of the string variable
   .  min - left endpoint of interval
   .  max - right endpoint of interval
   .  value - address where the value is stored

   DESCRIPTION:
   This function gets the double value of the string variable 'name'
   and checks wether it is lying in the specified interval.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 variable not found
   .n     2 no legal value
   .n     3/4 outside of range
   D*/
/****************************************************************************/

INT GetStringDOUBLEInRange (const char *name, DOUBLE min, DOUBLE max, DOUBLE *value)
{
  char *lastname;
  ENVDIR *theDir;
  STRVAR *myVar;
  double val;

  if ((theDir=FindStructDir(name,&lastname))==NULL) return(1);
  if ((myVar=FindStringVar(theDir,lastname))==NULL) return(1);
  if (sscanf(myVar->s,"%lf",&val)!=1) return(2);
  if (val<min) return(3);
  if (val>max) return(4);
  *value = (DOUBLE) val;

  return (0);
}

/****************************************************************************/
/*D
   GetStringINTInRange - Get the integer value of 'name' together with a range check

   SYNOPSIS:
   INT GetStringINTInRange (const char *name, INT min, INT max, INT *value);

   PARAMETERS:
   .  name - memory address of the string variable
   .  min - left endpoint of interval
   .  max - right endpoint of interval
   .  value - address where the value is stored

   DESCRIPTION:
   This function gets the integer value of the string variable 'name'
   and checks wether it is lying in the specified interval.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 variable not found
   .n     2 no legal value
   .n     3/4 outside of range
   D*/
/****************************************************************************/

INT GetStringINTInRange (const char *name, INT min, INT max, INT *value)
{
  char *lastname;
  ENVDIR *theDir;
  STRVAR *myVar;
  int val;

  if ((theDir=FindStructDir(name,&lastname))==NULL) return(1);
  if ((myVar=FindStringVar(theDir,lastname))==NULL) return(1);
  if (sscanf(myVar->s,"%d",&val)!=1) return(2);
  if (val<min) return(3);
  if (val>max) return(4);
  *value = (INT) val;

  return (0);
}

/****************************************************************************/
/*D
   GetCurrentStructDir - get current structure directory

   SYNOPSIS:
   ENVDIR *GetCurrentStructDir ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function returns the current structure directory.

   RETURN VALUE:
   ENVDIR *
   .n     pointer to the directory
   .n     NULL
   D*/
/****************************************************************************/

ENVDIR *GetCurrentStructDir ()
{
  return(path[pathIndex]);
}

/****************************************************************************/
/*D
   GetStructPathName - Assemble pathname of current structure directory

   SYNOPSIS:
   INT GetStructPathName (char *s, int n);

   PARAMETERS:
   .  s - pointer to buffer for the string
   .  n - length of buffer

   DESCRIPTION:
   This function assembles the pathname of the current structure directory.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT GetStructPathName (char *s, int n)
{
  int i,l;

  /* check length */
  l=2;
  for (i=1; i<=pathIndex; i++)
    l+=strlen(path[i]->name)+1;

  if (l>n)
    return(1);

  strcpy(s,STRUCTSEP);
  for (i=1; i<=pathIndex; i++)
  {
    strcat(s,path[i]->name); strcat(s,STRUCTSEP);
  }
  return(0);
}

/****************************************************************************/
/*D
   MakeStructItem - Allocate a new environment item in the current directory

   SYNOPSIS:
   ENVITEM *MakeStructItem (ENVDIR *where, const char *name, INT type,
   INT size);

   PARAMETERS:
   .  where - directory in which the new one will be allocated
   .  name - name of the new item
   .  type - type of the new item
   .  size - total size if user defined, string size for string variable

   DESCRIPTION:
   This function allocates a new environment item in the current directory.
   It is some reduced form of 'MakeEnvItem' with the additional feature,
   that it can start from a given directory.

   SEE ALSO:
   MakeEnvItem

   RETURN VALUE:
   ENVITEM *
   .n      pointer to
   .n      NULL if not enough memory in the environment heap or other error condition occurs
   D*/
/****************************************************************************/

ENVITEM *MakeStructItem (ENVDIR *where, const char *name, INT type, INT size)
{
  ENVDIR *currentDir;
  ENVITEM *newItem,*anItem,*lastItem;
  INT newsize;

  if (where!=NULL)
    currentDir=where;
  else
    currentDir=path[pathIndex];

  /* check if name too long */
  if (strlen(name)>NAMELEN) return(NULL);

  /* check if name not already used in this directory */
  anItem = lastItem = currentDir->down;
  while (anItem!=NULL)
  {
    if ((anItem->v.type==type)&&(strcmp(anItem->v.name,name)==0))
      return(NULL);

    lastItem = anItem;
    anItem = anItem->v.next;
  }

  /* allocate memory from environment heap */
  if (type==theStringVarID)
  {
    newsize = (1+size/32)*32;
    newItem = (ENVITEM *) AllocEnvMemory(sizeof(STRVAR)+newsize);
    if (newItem==NULL) return(NULL);
    ((STRVAR *) newItem)->length = newsize;
  }
  else if (type==theStringDirID)
  {
    if (pathIndex+1>=MAXENVPATH) return(NULL);
    newItem=(ENVITEM *) AllocEnvMemory(size);
    if (newItem==NULL) return(NULL);
    newItem->d.down = NULL;
  }
  else
    return(NULL);

  /* initialize item */
  newItem->v.type = type;
  newItem->v.locked = NO;
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
   MakeStruct - Allocate a new structure

   SYNOPSIS:
   INT MakeStruct (const char *name);

   PARAMETERS:
   .  name - name of the new directory

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT MakeStruct (const char *name)
{
  char *lastname;
  ENVDIR *theDir,*theStruct;

  if ((theDir=FindStructDir(name,&lastname))==NULL) return(1);

  if ((theStruct=FindStructure(theDir,lastname))!=NULL) return(0);
  if (MakeStructItem(theDir,lastname,theStringDirID,sizeof(ENVDIR))==NULL) return(2);

  return (0);
}

/****************************************************************************/
/*D
   DeleteStruct	- Delete an existing structure

   SYNOPSIS:
   INT DeleteStruct (char *name);

   PARAMETERS:
   .  name - structure name

   SEE ALSO:
   FindStructDir, FindStructure, CheckIfInStructPath, CheckStructTree, RemoveStructTree

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 structure directory not found
   .n     2 structure does not exist
   .n     3 structure is inside structure path
   .n     4 structure contains locked objects
   .n     5 structure could not be removed
   D*/
/****************************************************************************/

INT DeleteStruct (char *name)
{
  char *lastname;
  ENVDIR *theDir,*theStruct;

  if ((theDir=FindStructDir(name,&lastname))==NULL)
    return(1);                  /* structure directory not found */

  if ((theStruct=FindStructure(theDir,lastname))==NULL)
    return(2);                  /* structure does not exist */

  if (CheckIfInStructPath(theStruct))
    return(3);                  /* structure is inside structure path */

  if (CheckStructTree(theStruct))
    return(4);                  /* structure contains locked objects */

  if (RemoveStructTree(theDir,theStruct)!=0)
    return(5);                  /* structure could not be removed */

  return (0);
}

/****************************************************************************/
/*D
   RemoveStringVar - remove a string variable fro the environment tree

   SYNOPSIS:
   INT RemoveStringVar (ENVDIR *homeDir, STRVAR *theVar)

   PARAMETERS:
   .  homeDir - theVar is located in this directory
   .  theVar - string var to remove

   DESCRIPTION:
   This function removes a string variable from the environment tree.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 homeDir==NULL
   .n     2 theVar==NULL
   .n     3 theVar is a directory
   D*/
/****************************************************************************/

INT RemoveStringVar (ENVDIR *homeDir, STRVAR *theVar)
{
  if (homeDir==NULL) return (1);
  if (theVar==NULL) return (2);
  if (theVar->v.type%2!=0) return (3);

  /* remove item from double linked list */
  if (theVar->v.previous==NULL)
    homeDir->down = theVar->v.next;
  else
    theVar->v.previous->v.next = theVar->v.next;

  if (theVar->v.next!=NULL)
    theVar->v.next->v.previous = theVar->v.previous;

  /* deallocate memory */
  FreeEnvMemory(theVar);

  return (0);
}

/****************************************************************************/
/*D
   DeleteVariable	- Delete an existing string variable

   SYNOPSIS:
   INT DeleteVariable (char *name);

   PARAMETERS:
   .  name - variable name

   SEE ALSO:
   FindStructDir, FindStructure, CheckIfInStructPath, CheckStructTree, RemoveStructTree

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 variable directory not found
   .n     2 variable does not exist
   .n     4 variable contains locked objects
   .n     5 variable could not be removed
   D*/
/****************************************************************************/

INT DeleteVariable (char *name)
{
  char *lastname;
  ENVDIR *theDir;
  STRVAR *myVar;

  if ((theDir=FindStructDir(name,&lastname))==NULL)
    return(1);                  /* structure directory not found */

  if ((myVar=FindStringVar(theDir,lastname))==NULL)
    return(2);                  /* variable does not exist */

  if (ENVITEM_LOCKED(myVar))
    return(4);                  /* variable is locked */

  if (RemoveStructTree(theDir,(ENVDIR *) myVar)!=0)
    return(5);                  /* variable could not be removed */

  return (0);
}

/****************************************************************************/
/*D
   SetStringVar - Set a string variable to a given string

   SYNOPSIS:
   INT SetStringVar (const char *name, char *sval);

   PARAMETERS:
   .  name - variable name
   .  sval - address of the string

   DESCRIPTION:
   This function searches a string variable and sets it to the given string.
   If the string variable does not yet exist it is created, if it is
   too short for the string, it is removed and newly created.

   SEE ALSO:
   FindStructDir, FindStringVar, RemoveStringVar, MakeStructItem

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 structure directory not found
   .n     2 could not allocate variable
   D*/
/****************************************************************************/

INT SetStringVar (const char *name, char *sval)
{
  char *lastname;
  ENVDIR *theDir;
  STRVAR *myVar;

  if ((theDir=FindStructDir(name,&lastname))==NULL)
    return(1);                  /* structure directory not found */

  myVar=FindStringVar(theDir,lastname);

  if ((myVar!=NULL) && (myVar->length<=strlen(sval)))
  {
    RemoveStringVar(theDir, myVar);
    myVar=NULL;
  }

  if (myVar==NULL)
  {
    myVar = (STRVAR *) MakeStructItem(theDir,lastname,theStringVarID,strlen(sval));
    if (myVar==NULL)
      return(2);                        /* could not allocate variable */
  }

  strcpy(myVar->s,sval);

  return(0);
}

/****************************************************************************/
/*D
   SetStringVarNotify - Set a string variable to a given string and notify if changed

   SYNOPSIS:
   INT SetStringVarNotify (const char *name, const char *sval);

   PARAMETERS:
   .  name - variable name
   .  sval - address of the string

   DESCRIPTION:
   This function searches a string variable and sets it to the given string.
   If the string variable does not yet exist it is created, if it is
   too short for the string, it is removed and newly created.

   SEE ALSO:
   FindStructDir, FindStringVar, RemoveStringVar, MakeStructItem

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 structure directory not found
   .n     2 could not allocate variable
   D*/
/****************************************************************************/

INT SetStringVarNotify (const char *name, const char *sval)
{
  char *lastname;
  ENVDIR *theDir;
  STRVAR *myVar;
  int notify = SV_NOT_CHANGED;

  if ((theDir=FindStructDir(name,&lastname))==NULL)
    return SV_ERROR;                    /* structure directory not found */

  myVar=FindStringVar(theDir,lastname);

  if ((myVar!=NULL) && (myVar->length<=strlen(sval)))
  {
    RemoveStringVar(theDir, myVar);
    myVar=NULL;
  }

  if (myVar==NULL)
  {
    myVar = (STRVAR *) MakeStructItem(theDir,lastname,theStringVarID,strlen(sval));
    if (myVar==NULL)
      return SV_ERROR;                          /* could not allocate variable */
    notify = SV_CREATED;
  }

  if (notify==SV_NOT_CHANGED && strcmp(myVar->s,sval))
    notify = SV_CHANGED;
  strcpy(myVar->s,sval);

  return notify;
}

/****************************************************************************/
/*D
   SetStringValue - Set a variable to a real value

   SYNOPSIS:
   INT SetStringValue (const char *name, double value);

   PARAMETERS:
   .  name - variable name
   .  value - value

   DESCRIPTION:
   This function sets a variable to a real value.

   SEE ALSO:
   SetStringVar

   RETURN VALUE:
   INT
   .n     see 'SetStringVar'
   D*/
/****************************************************************************/

INT SetStringValue (const char *name, double value)
{
  char buffer[30];

  sprintf(buffer,"%-.14lg",value);
  return(SetStringVar(name,buffer));
}


/****************************************************************************/
/*D
   CheckStructTree - Check if LOCKED-Flag is set somewhere

   SYNOPSIS:
   INT CheckStructTree (const ENVDIR *theDir);

   PARAMETERS:
   .  theDir - directory in which struct is located

   DESCRIPTION:
   This function checks (recursively) if the LOCKED-Flag is set somewhere
   inside the directory tree starting at 'theDir'.
   Up to now the locking is not yet used.

   RETURN VALUE:
   INT
   .n     1 if LOCKED-Flag set somewhere
   .n     0 else.
   D*/
/****************************************************************************/

INT CheckStructTree (const ENVDIR *theDir)
{
  INT status;

  if (ENVITEM_LOCKED(theDir))
    return(1);

  if (theDir->type%2!=0)
    for (theDir=(ENVDIR *) theDir->down; theDir!=NULL; theDir=(ENVDIR *) theDir->next)
      if ((status=CheckStructTree(theDir))!=0)
        return(status);

  return(0);
}

/****************************************************************************/
/*D
   CheckIfInStructPath - Searches '*theDir' inside the current structure path

   SYNOPSIS:
   INT CheckIfInStructPath (const ENVDIR *theDir);

   PARAMETERS:
   .  theDir - structure

   DESCRIPTION:
   See above.

   RETURN VALUE:
   INT
   .n     1 if inside path
   .n     0 if not inside path.
   D*/
/****************************************************************************/

INT CheckIfInStructPath (const ENVDIR *theDir)
{
  int i;

  for (i=0; i<=pathIndex; i++)
    if (theDir==path[i])
      return(1);

  return(0);
}

/****************************************************************************/
/*D
   RemoveStructTree - Removes the whole structure tree

   SYNOPSIS:
   INT RemoveStructTree (ENVDIR *homeDir, ENVDIR *theDir);

   PARAMETERS:
   .  homeDir - home directory
   .  theDir - structure

   DESCRIPTION:
   This function removes recursively the whole structure tree.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT RemoveStructTree (ENVDIR *homeDir, ENVDIR *theDir)
{
  ENVDIR *theDir2;

  if (theDir->type%2!=0)
    for (theDir2=(ENVDIR *) theDir->down; theDir2!=NULL; theDir2=(ENVDIR *) theDir2->next)
      RemoveStructTree(theDir,theDir2);

  /* remove item from double linked list */
  if (theDir->previous==NULL)
    homeDir->down = theDir->next;
  else
    theDir->previous->d.next = theDir->next;

  if (theDir->next!=NULL)
    theDir->next->v.previous = theDir->previous;

  /* deallocate memory */
  FreeEnvMemory(theDir);

  return(0);
}

/****************************************************************************/
/*D
   DirOut - (recursively) print directory contents

   SYNOPSIS:
   static INT DirOut (const ENVDIR *theDir, char *buffer, int bufLen, int ropt);

   PARAMETERS:
   .  theDir - structure
   .  buffer - address of buffer
   .  bufLen - length of buffer
   .  ropt - if set, the contents are printed recursively

   DESCRIPTION:
   This function (recursively) prints the contents of theDir into the buffer.
   It will return the value 4 if the buffer is full and it is not yet ready.
   In this case it has to be called another time. If it is ready it returns
   the value 0.

   RETURN VALUE:
   INT
   .n    0 all is done
   .n    1 buffer too small
   .n    2 theDir is not a string directory
   .n    3 false type encountered (should never occur)
   .n    4 not yet ready with output (this is not an error!)
   .n    5 directory tree too deeply nested (should never occur)
   .n    6 this should also never occur
   D*/
/****************************************************************************/

static INT DirOut (const ENVDIR *theDir, char *buffer, int bufLen, int ropt)
{
  int bufPos, i, n;
  static ENVITEM *theItem;
  static int status,pathPos;
  static const ENVDIR *varPath[MAXENVPATH];
  static char *str;

  bufPos=0;             /* the buffer is initially empty */

  if (bufLen<MAXENVPATH+NAMESIZE+10)
    return(1);                  /* buffer too small */

  if (theDir!=NULL)
  {
    /* a first call */
    if (ENVITEM_TYPE(theDir)!=theStringDirID)
      return(2);                        /* theDir is not a string directory */
    pathPos=0;
    varPath[0]=theDir;
    theItem=ENVITEM_DOWN(theDir);
    status=0;                   /* means that for theItem all has to be done */
  }

  do
  {
    switch (status)
    {
    case 0 :
      while (theItem==NULL)
      {
        if (pathPos==0)
        {
          buffer[bufPos]=(char) 0;
          return(0);                                    /* ready with output */
        }

        /* we have to print a terminating "\t...\t}\n\0" and to go back in our path */
        if (bufLen-bufPos<pathPos+2)
        {
          buffer[bufPos]=(char) 0;
          return(4);                                    /* not yet ready with output */
        }

        for (i=0; i<pathPos-1; i++)
          buffer[bufPos++]='\t';

        buffer[bufPos++]='}'; buffer[bufPos++]='\n';

        theItem=NEXT_ENVITEM(varPath[pathPos--]);
      }

      if ((ENVITEM_TYPE(theItem)!=theStringDirID)&&(ENVITEM_TYPE(theItem)!=theStringVarID))
        return(3);                              /* false type encountered (should not occur) */

      if (bufLen-bufPos<pathPos+1)
      {
        buffer[bufPos]=(char) 0;
        return(4);                              /* not yet ready with output */
      }

      for (i=0; i<pathPos; i++)
        buffer[bufPos++]='\t';

      status=1;

    case 1 :
      n=strlen(ENVITEM_NAME(theItem));
      if (bufLen-bufPos<n+7)                    /* "<name> = {}\n\0" */
      {
        buffer[bufPos]=(char) 0;
        return(4);                              /* not yet ready with output */
      }

      strcpy(buffer+bufPos,ENVITEM_NAME(theItem)); bufPos+=n;
      strcpy(buffer+bufPos," = "); bufPos+=3;
      status=2;

    case 2 :
      if (ENVITEM_TYPE(theItem)==theStringDirID)
      {
        if ((ropt==0)||(ENVITEM_DOWN(theItem)==NULL))
        {
          strcpy(buffer+bufPos,"{}\n"); bufPos+=3;
          theItem=NEXT_ENVITEM(theItem);
        }
        else
        {
          buffer[bufPos++]='{'; buffer[bufPos++]='\n';
          if (pathPos==MAXENVPATH-1)
            return(5);                                          /* directory tree too deeply nested: should not occur! */

          varPath[++pathPos]=(ENVDIR *) theItem;
          theItem=ENVITEM_DOWN(theItem);
        }
        break;
      }
      else
      {
        str=((STRVAR *)theItem)->s;
        status=3;
      }

    case 3 :
      strncpy(buffer+bufPos,str,bufLen-bufPos-2);
      n=strlen(str);
      if (bufLen-bufPos-2<n)
      {
        str+=bufLen-bufPos-2;
        buffer[bufLen-2]=(char) 0;
        return(4);                              /* not yet ready with output */
      }
      else
      {
        bufPos+=n;
        buffer[bufPos++]='\n';
      }
      theItem=NEXT_ENVITEM(theItem);
    }

    status=0;
  } while (TRUE);

  return(6);            /* this should never be reached */
}

/****************************************************************************/
/*D
   VarOut - prints variable = content into a buffer

   SYNOPSIS:
   static INT VarOut (const STRVAR *StrVar, char *buffer, int bufLen,);

   PARAMETERS:
   .  StrVar - string variable
   .  buffer - address of buffer
   .  bufLen - length of buffer

   DESCRIPTION:
   This function prints the contents of StrVar into the buffer.
   It will return the value 4 if the buffer is full and it is not yet ready.
   In this case it has to be called another time. If it is ready it returns
   the value 0.

   RETURN VALUE:
   INT
   .n    0 all is done
   .n    1 buffer too small
   .n    4 not yet ready with output (this is not an error!)
   D*/
/****************************************************************************/

static INT VarOut (const STRVAR *StrVar, char *buffer, int bufLen)
{
  static const char *str;

  if (bufLen<MAXENVPATH+NAMESIZE+10)
    return(1);                  /* buffer too small */

  if (StrVar!=NULL)
  {
    strcpy(buffer,ENVITEM_NAME(StrVar));
    buffer+=strlen(ENVITEM_NAME(StrVar)); bufLen-=strlen(ENVITEM_NAME(StrVar));
    strcpy(buffer," = ");
    buffer+=3; bufLen-=3;
    str=StrVar->s;
  }

  if (strlen(str)+2<bufLen)
  {
    strcpy(buffer, str);
    strcat(buffer, "\n");
    return(0);
  }
  else
  {
    strncpy(buffer, str, bufLen-1);
    str+=bufLen-1;
    buffer[bufLen-1]=(char) 0;
    return(4);
  }
}

/****************************************************************************/
/*D
   PrintStructContents - (recursively) print structure contents

   SYNOPSIS:
   INT PrintStructContents (const char *name, char *buffer, int bufLen, int ropt);

   PARAMETERS:
   .  name - name of structure(dir)
   .  buffer - address of buffer
   .  bufLen - length of buffer
   .  ropt - if set, the contents are printed recursively

   DESCRIPTION:
   This function prints the contents of the variable or structure 'name'.
   More or less it provides a call for the functions DirOut and VarOut.
   You should take care with error handling when using these function.

   SEE ALSO:
   DirOut, VarOut, FindStructDir, FindStringVar, FindStructure

   RETURN VALUE:
   INT
   .n    0 all is done
   .n    1 buffer too small
   .n    2, 3 (should not occur)
   .n    4 not yet ready with output (this is not an error!)
   .n    5, 6 (should not occur)
   .n    7 structure path not found
   .n    8 structure not found
   D*/
/****************************************************************************/

INT PrintStructContents (const char *name, char *buffer, int bufLen, int ropt)
{
  static ENVDIR *theDir;
  static STRVAR *StrVar;
  static int status;
  char *lastname;
  int ret;

  buffer[0]=(char) 0;

  if (name!=NULL)
  {
    /* get variable and/or structure with that name */
    if (strcmp(name,":")==0)
    {
      StrVar=NULL; theDir=path[0];
    }
    else
    {
      if ((theDir=FindStructDir(name,&lastname))==NULL)
        return(7);                              /* structure path not found */
      StrVar=FindStringVar(theDir,lastname);
      theDir=FindStructure(theDir,lastname);
    }
    status=0;
  }

  if (status==0)
    status=(StrVar==NULL) ? 2 : 1;

  if (status==1)
  {
    ret=VarOut(StrVar, buffer, bufLen);
    if ((ret!=0)&&(ret!=4))
      return(ret);
    if (ret==0)
      status=2;
    else
      StrVar=NULL;
    return(4);                  /* also return if a variable was printed to clear the buffer */
  }

  if (status==2)
    status=(theDir==NULL) ? 4 : 3;

  if (status==3)
  {
    ret=DirOut(theDir, buffer, bufLen, ropt);
    if ((ret!=0)&&(ret!=4))
      return(ret);
    if (ret==4)
    {
      theDir=NULL;
      return(4);
    }
  }

  return(0);
}

/****************************************************************************/
/*D
   PrintCurrentStructContents - (recursively) print contents of current structure directory

   SYNOPSIS:
   INT PrintCurrentStructContents (int flag, char *buffer, int bufLen, int ropt);

   PARAMETERS:
   .  flag - 1 for the first call, 0 if others are necessary
   .  buffer - address of buffer
   .  bufLen - length of buffer
   .  ropt - if set, the contents are printed recursively

   DESCRIPTION:
   This function prints the contents of the of the current working structure.
   More or less it calls the function DirOut.

   SEE ALSO:
   DirOut, PrintStructContents, FindStructDir, FindStringVar, FindStructure

   RETURN VALUE:
   INT
   .n    0 all is done
   .n    1 buffer too small
   .n    2, 3 (should not occur)
   .n    4 not yet ready with output (this is not an error!)
   .n    5, 6 (should not occur)
   D*/
/****************************************************************************/

INT PrintCurrentStructContents (int flag, char *buffer, int bufLen, int ropt)
{
  if (flag)
    return (DirOut(path[pathIndex], buffer, bufLen, ropt));
  else
    return (DirOut(NULL, buffer, bufLen, ropt));
}

/****************************************************************************/
/*D
   InitUgStruct	- Initialize ugstruct

   SYNOPSIS:
   INT InitUgStruct ();

   PARAMETERS:
   .  none

   DESCRIPTION:
   This function creates the '/Strings'-directory, gets an ID for
   string directories and sets the current structure path to '/Strings'.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if not enough memory.
   D*/
/****************************************************************************/

INT InitUgStruct ()
{
  ENVDIR *theDir;

  /* install the /Strings directory */
  if (ChangeEnvDir("/")==NULL)
    return(__LINE__);

  theStringDirID = GetNewEnvDirID();
  if (MakeEnvItem("Strings",theStringDirID,sizeof(ENVDIR))==NULL)
    return(__LINE__);
  theStringVarID = GetNewEnvVarID();

  /* set structure path to "/strings" */
  if ((theDir=ChangeEnvDir("/Strings"))==NULL)
    return(__LINE__);

  path[pathIndex=0]=theDir;

  /* return ok */
  return(0);
}
