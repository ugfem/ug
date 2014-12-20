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

#include <config.h>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "ugtypes.h"
#include "general.h"
#include "heaps.h"
#include "misc.h"
#include "ugenv.h"

#include "ugstruct.h"


USING_UG_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define qErr -1

#define MAXENVPATH              32              /* max depth of the environment tree */
#define STRUCTSEP ":"
#define STRUCTSEPC ':'

/****************************************************************************/
/*                                                                          */
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
/** \brief Change current structure directory

   \param name - name to which the current structure directory shall be changed

   Changes the current structure directory to the one given by name. The name is
   interpreted as starting from the previous end of the path. It is possible to
   use ".." to go to the parent directory and ":" at the beginning
   to start from the root directory.

   \return
   The new directory or NULL, if it was not found.
 */
/****************************************************************************/

ENVDIR * NS_PREFIX ChangeStructDir (const char *name)
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

static char token[NAMESIZE];                    /* must be global for lastnameHnd ! */
static char nexttoken[NAMESIZE];        /* must be global for lastnameHnd ! */

/****************************************************************************/
/** \brief Finds the directory with the given name or its parent

   \param name - name to be found (from current structure directory)
   \param lastnameHnd - if not NULL it is set to a pointer to the last name

   This function searchs (starting from the current structure directory)
   for the structure directory given by name and returns its adress.
   It does not change the current structure directory! If 'lastnameHnd!=NULL'
   it also sets '*lastNameHnd' to be a pointer to the last name of the directory.

   \sa
   strntok

   \return
   The directory, or NULL if not found.
 */
/****************************************************************************/

ENVDIR * NS_PREFIX FindStructDir (const char *name, const char **lastnameHnd)
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
/** \brief Searches a string variable inside a directory/structure

   \param where - the directory to be searched
   \param name - the name to be searched

   \return
   The address of the variable or NULL, if not found.
 */
/****************************************************************************/

STRVAR * NS_PREFIX FindStringVar (const ENVDIR *where, const char *name)
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
/** \brief Searches a directory/structure inside a directory/structure

   \param where - the directory to be searched (if NULL then search is in root)
   \param name - the name to be searched

   \return
   The address of the structure or NULL, if not found.
 */
/****************************************************************************/

ENVDIR * NS_PREFIX FindStructure (const ENVDIR *where, const char *name)
{
  ENVDIR *theDir;

  if (where==NULL) where=path[0];

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
/** \brief Get string contained in the variable 'name'

   \param  name - name of the string variable (from the current structure directory)

    This function returns the string contained in the variable 'name'.

    \return
   Pointer to the string, or NULL if not found
 */
/****************************************************************************/

char * NS_PREFIX GetStringVar (const char *name)
{
  const char *lastname;
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
/** \brief Get double value of the string contained in the variable 'name'

   \param name - name of the string variable
   \param value - address where the result is stored

   Same as GetStringValueDouble.

   \todo Should not be used and may be skipped in a future version.

   \return <ul>
   <li> 0 if ok </li>
   <li> 1 if error occured </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX GetStringValue (const char *name, double *value)
{
  const char *lastname;
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
/** \brief Get double value of a variable

   \param name - name of the string variable
   \param value - address where the result is stored

   This function evaluates a string variable holding a floating point number
   and returns the result.

   \return <ul>
   <li> 0 if ok </li>
   <li> 1 if error occured </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX GetStringValueDouble (const char *name, double *value)
{
  const char *lastname;
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
/** \brief Get integer value of a variable

   \param name - name of the string variable
   \param value - place result here

   This function evaluates a string variable holding an integer number
   and returns the result.

   \return <ul>
   <li> 0 if ok </li>
   <li> 1 if error occured </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX GetStringValueInt (const char *name, int *value)
{
  const char *lastname;
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
/** \brief Get the double value of 'name' together with a range check

   \param name - memory address of the string variable
   \param min - left endpoint of interval
   \param max - right endpoint of interval
   \param value - address where the value is stored

   This function gets the double value of the string variable 'name'
   and checks wether it is lying in the specified interval.

   \return <ul>
   <li>  0 if ok </li>
   <li>  1 variable not found </li>
   <li>  2 no legal value </li>
   <li>  3/4 outside of range </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX GetStringDOUBLEInRange (const char *name, DOUBLE min, DOUBLE max, DOUBLE *value)
{
  const char *lastname;
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
/** \brief Get the integer value of 'name' together with a range check

   \param name - memory address of the string variable
   \param min - left endpoint of interval
   \param max - right endpoint of interval
   \param value - address where the value is stored

   This function gets the integer value of the string variable 'name'
   and checks wether it is lying in the specified interval.

   \return <ul>
   <li>  0 if ok </li>
   <li>  1 variable not found </li>
   <li>  2 no legal value </li>
   <li>  3/4 outside of range </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX GetStringINTInRange (const char *name, INT min, INT max, INT *value)
{
  const char *lastname;
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
/** \brief Get current structure directory

   This function returns the current structure directory.

   \return <ul>
   <li>  pointer to the directory </li>
   <li>  NULL </li>
   </ul>
 */
/****************************************************************************/

ENVDIR * NS_PREFIX GetCurrentStructDir ()
{
  return(path[pathIndex]);
}

/****************************************************************************/
/** \brief Assemble pathname of current structure directory

   \param s - pointer to buffer for the string
   \param n - length of buffer

   This function assembles the pathname of the current structure directory.

   \return
   <li>  0 if ok </li>
   <li>  1 if error occured. </li>
 */
/****************************************************************************/

INT NS_PREFIX GetStructPathName (char *s, int n)
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
/** \brief Allocate a new environment item in the current directory

   \param where - directory in which the new one will be allocated
   \param name - name of the new item
   \param type - type of the new item
   \param size - total size if user defined, string size for string variable

   This function allocates a new environment item in the current directory.
   It is some reduced form of 'MakeEnvItem' with the additional feature,
   that it can start from a given directory.

   \sa
   MakeEnvItem

   \return <ul>
   <li>   pointer to </li>
   <li>   NULL if not enough memory in the environment heap or other error condition occurs </li>
   </ul>
 */
/****************************************************************************/

ENVITEM * NS_PREFIX MakeStructItem (ENVDIR *where, const char *name, INT type, INT size)
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
/** \brief Allocate a new structure

   \param name - name of the new directory

   \return <ul>
   <li>  0 if ok </li>
   <li>  1 if error occured </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX MakeStruct (const char *name)
{
  const char *lastname;
  ENVDIR *theDir,*theStruct;

  if ((theDir=FindStructDir(name,&lastname))==NULL) return(1);

  if ((theStruct=FindStructure(theDir,lastname))!=NULL) return(0);
  if (MakeStructItem(theDir,lastname,theStringDirID,sizeof(ENVDIR))==NULL) return(2);

  return (0);
}

/****************************************************************************/
/** \brief Delete an existing structure

   \param name - structure name

   SEE ALSO:
   FindStructDir, FindStructure, CheckIfInStructPath, CheckStructTree, RemoveStructTree

   \return <ul>
   <li>  0 if ok </li>
   <li>  1 structure directory not found </li>
   <li>  2 structure does not exist </li>
   <li>  3 structure is inside structure path </li>
   <li>  4 structure contains locked objects </li>
   <li>  5 structure could not be removed </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX DeleteStruct (const char *name)
{
  const char *lastname;
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
/** \brief Remove a string variable fro the environment tree

   \param homeDir - theVar is located in this directory
   \param theVar - string var to remove

   This function removes a string variable from the environment tree.

   \return <ul>
   <li>  0 if ok </li>
   <li>  1 homeDir==NULL </li>
   <li>  2 theVar==NULL </li>
   <li>  3 theVar is a directory </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX RemoveStringVar (ENVDIR *homeDir, STRVAR *theVar)
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
/** \brief Delete an existing string variable

   \param name - variable name

   SEE ALSO:
   FindStructDir, FindStructure, CheckIfInStructPath, CheckStructTree, RemoveStructTree

   \return <ul>
   <li>  0 if ok </li>
   <li>  1 variable directory not found </li>
   <li>  2 variable does not exist </li>
   <li>  4 variable contains locked objects </li>
   <li>  5 variable could not be removed </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX DeleteVariable (const char *name)
{
  const char *lastname;
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
/** \brief Set a string variable to a given string

   \param name - variable name
   \param sval - address of the string

   This function searches a string variable and sets it to the given string.
   If the string variable does not yet exist it is created, if it is
   too short for the string, it is removed and newly created.

   \sa
   FindStructDir, FindStringVar, RemoveStringVar, MakeStructItem

   \return <ul>
   <li>  0 if ok </li>
   <li>  1 structure directory not found </li>
   <li>  2 could not allocate variable </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX SetStringVar (const char *name, const char *sval)
{
  const char *lastname;
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
/** \brief Set a string variable to a given string

   \param name - variable name
   \param sval - address of the string
   \param n - length of string

   This function searches a string variable and sets it to the given string.
   If the string variable does not yet exist it is created, if it is
   too short for the string, it is removed and newly created.

   The difference to 'SetStringVar' is that the string is not terminated
   by '\0'.

   \sa
   SetStringVar, FindStructDir, FindStringVar, RemoveStringVar, MakeStructItem

   \return <ul>
   <li>  0 if ok </li>
   <li>  1 structure directory not found </li>
   <li>  2 could not allocate variable </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX SetnStringVar (const char *name, const char *sval, int n)
{
  const char *lastname;
  ENVDIR *theDir;
  STRVAR *myVar;

  if ((theDir=FindStructDir(name,&lastname))==NULL)
    return(1);                  /* structure directory not found */

  myVar=FindStringVar(theDir,lastname);

  if ((myVar!=NULL) && (myVar->length<=n))
  {
    RemoveStringVar(theDir, myVar);
    myVar=NULL;
  }

  if (myVar==NULL)
  {
    myVar = (STRVAR *) MakeStructItem(theDir,lastname,theStringVarID,n);
    if (myVar==NULL)
      return(2);                        /* could not allocate variable */
  }

  strncpy(myVar->s,sval,(size_t) n);
  myVar->s[n] = '\0';

  return(0);
}
/****************************************************************************/
/** \brief Set a string variable to a given string and notify if changed

   \param name - variable name
   \param sval - address of the string

   This function searches a string variable and sets it to the given string.
   If the string variable does not yet exist it is created, if it is
   too short for the string, it is removed and newly created.

   \sa
   FindStructDir, FindStringVar, RemoveStringVar, MakeStructItem

   \return <ul>
   <li>  0 if ok </li>
   <li>  1 structure directory not found </li>
   <li>  2 could not allocate variable </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX SetStringVarNotify (const char *name, const char *sval)
{
  const char *lastname;
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
/** \brief Set a variable to a real value

   \param name - variable name
   \param value - value

   This function sets a variable to a real value.

   \sa
   SetStringVar

   \return
    see 'SetStringVar'
 */
/****************************************************************************/

INT NS_PREFIX SetStringValue (const char *name, double value)
{
  char buffer[30];

  sprintf(buffer,"%-.14g",value);
  return(SetStringVar(name,buffer));
}


/****************************************************************************/
/** \brief Check if LOCKED-Flag is set somewhere

   \param theDir - directory in which struct is located

   This function checks (recursively) if the LOCKED-Flag is set somewhere
   inside the directory tree starting at 'theDir'.
   Up to now the locking is not yet used.

   \return <ul>
   <li>  1 if LOCKED-Flag set somewhere </li>
   <li>  0 else </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX CheckStructTree (const ENVDIR *theDir)
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
/** \brief Searches '*theDir' inside the current structure path

   \param theDir - structure

   \return <ul>
   <li>  1 if inside path </li>
   <li>  0 if not inside path </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX CheckIfInStructPath (const ENVDIR *theDir)
{
  int i;

  for (i=0; i<=pathIndex; i++)
    if (theDir==path[i])
      return(1);

  return(0);
}

/****************************************************************************/
/** \brief Removes the whole structure tree

   \param homeDir - home directory
   \param theDir - structure

   This function removes recursively the whole structure tree.

   \return <ul>
   <li>  0 if ok </li>
   <li>  1 if error occured </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX RemoveStructTree (ENVDIR *homeDir, ENVDIR *theDir)
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
/** \brief (recursively) print directory contents

   \param theDir - structure
   \param buffer - address of buffer
   \param bufLen - length of buffer
   \param ropt - if set, the contents are printed recursively

   This function (recursively) prints the contents of theDir into the buffer.
   It will return the value 4 if the buffer is full and it is not yet ready.
   In this case it has to be called another time. If it is ready it returns
   the value 0.

   \return <ul>
   <li> 0 all is done </li>
   <li> 1 buffer too small </li>
   <li> 2 theDir is not a string directory </li>
   <li> 3 false type encountered (should never occur) </li>
   <li> 4 not yet ready with output (this is not an error!) </li>
   <li> 5 directory tree too deeply nested (should never occur) </li>
   <li> 6 this should also never occur </li>
   </ul>
 */
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
  } while (true);

  return(6);            /* this should never be reached */
}

/****************************************************************************/
/** \brief Prints variable = content into a buffer

   \param StrVar - string variable
   \param buffer - address of buffer
   \param bufLen - length of buffer

   This function prints the contents of StrVar into the buffer.
   It will return the value 4 if the buffer is full and it is not yet ready.
   In this case it has to be called another time. If it is ready it returns
   the value 0.

   \return <ul>
   <li> 0 all is done </li>
   <li> 1 buffer too small </li>
   <li> 4 not yet ready with output (this is not an error!) </li>
   </ul>
 */
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
/** \brief (recursively) print structure contents

   \param name - name of structure(dir)
   \param buffer - address of buffer
   \param bufLen - length of buffer
   \param ropt - if set, the contents are printed recursively

   This function prints the contents of the variable or structure 'name'.
   More or less it provides a call for the functions DirOut and VarOut.
   You should take care with error handling when using these function.

   \sa
   DirOut, VarOut, FindStructDir, FindStringVar, FindStructure

   \return <ul>
   <li> 0 all is done </li>
   <li> 1 buffer too small </li>
   <li> 2, 3 (should not occur) </li>
   <li> 4 not yet ready with output (this is not an error!) </li>
   <li> 5, 6 (should not occur) </li>
   <li> 7 structure path not found </li>
   <li> 8 structure not found </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX PrintStructContents (const char *name, char *buffer, int bufLen, int ropt)
{
  static ENVDIR *theDir;
  static STRVAR *StrVar;
  static int status;
  const char *lastname;
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
/** \brief (recursively) print contents of current structure directory

   \param flag - 1 for the first call, 0 if others are necessary
   \param buffer - address of buffer
   \param bufLen - length of buffer
   \param ropt - if set, the contents are printed recursively

   This function prints the contents of the of the current working structure.
   More or less it calls the function DirOut.

   \sa
   DirOut, PrintStructContents, FindStructDir, FindStringVar, FindStructure

   \return <ul>
   <li> 0 all is done </li>
   <li> 1 buffer too small </li>
   <li> 2, 3 (should not occur) </li>
   <li> 4 not yet ready with output (this is not an error!) </li>
   <li> 5, 6 (should not occur) </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX PrintCurrentStructContents (int flag, char *buffer, int bufLen, int ropt)
{
  if (flag)
    return (DirOut(path[pathIndex], buffer, bufLen, ropt));
  else
    return (DirOut(NULL, buffer, bufLen, ropt));
}

/****************************************************************************/
/** \brief Initialize ugstruct

   This function creates the '/Strings'-directory, gets an ID for
   string directories and sets the current structure path to '/Strings'.

   \return <ul>
   <li>  0 if ok </li>
   <li>  1 if not enough memory </li>
   </ul>
 */
/****************************************************************************/

INT NS_PREFIX InitUgStruct ()
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
