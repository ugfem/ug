// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  numproc.c	                                                                                                */
/*																			*/
/* Purpose:   definition of the basic num proc type                                 */
/*																			*/
/* Author:	  Peter Bastian                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   December 12, 1996                                                                         */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "general.h"
#include "gm.h"
#include "ugenv.h"
#include "devices.h"
#include "np.h"
#include "numproc.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define MAXCLASSES              20

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static INT ClassDirID;
static INT ObjectDirID;
static INT ClassVarID;
static INT ObjectVarID;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*D
   CreateClass - create constructor for a new class

   SYNOPSIS:
   INT CreateClass (char *classname, INT size, ConstructorProcPtr Construct);

   PARAMETERS:
   .  classname - name for new class. Must conform to class naming rules.
   .  size - size of objects of this class.
   .  Construct - pointer to constructor function for objects of this class.

   DESCRIPTION:
   This function is used by a class to bring itself into existence
   at initialization time.
   It creates a new class of numerical procedure (TOFKANT - The
   objects formerly known as numproc types :-) ).

   RETURN VALUE:
   INT
   .n    0 no error.
   .n     1 error, usually memory overflow in environment.
   D*/
/****************************************************************************/

INT CreateClass (const char *classname, INT size, ConstructorProcPtr Construct)
{
  NP_CONSTRUCTOR *constructor;

  /* create class directory under root if necessary */
  if (ChangeEnvDir("/") == NULL) return (1);
  if (ChangeEnvDir("NumProcClasses") == NULL) {
    MakeEnvItem("NumProcClasses",ClassDirID,sizeof(ENVDIR));
    if (ChangeEnvDir("NumProcClasses") == NULL) return (1);
  }
  /* allocate memory for constructor structure */
  constructor = (NP_CONSTRUCTOR *) MakeEnvItem (classname,ClassVarID,sizeof(NP_CONSTRUCTOR));
  if (constructor==NULL) return(1);

  /* fill constructor structure */
  constructor->size = size;
  constructor->Construct = Construct;

  /* return OK */
  return(0);
}
/****************************************************************************/
/*D
   GetConstructor - return constructor for given class name

   SYNOPSIS:
   NP_CONSTRUCTOR *GetConstructor (const char *classname)

   PARAMETERS:
   .  classname - final class name of an existing class.

   DESCRIPTION:
   The constructor for the given class name is returned.

   RETURN VALUE:
   NP_CONSTRUCTOR *
   .n    constructor
   .n     NULL if not found
   D*/
/****************************************************************************/

NP_CONSTRUCTOR *GetConstructor (const char *classname)
{
  ENVITEM *item;
  INT m,i;

  item = (ENVITEM *) ChangeEnvDir("/NumProcClasses");
  if (item == NULL) return (NULL);
  for (item=ENVITEM_DOWN(item); item!=NULL; item=NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == ClassVarID)
    {
      /* check classname */
      m = strlen(ENVITEM_NAME(item));
      for (i=m-1; i>=0; i--)
        if (ENVITEM_NAME(item)[i]=='.') break;
      if (strcmp(ENVITEM_NAME(item)+(i+1),classname)==0) break;
    }
  return ((NP_CONSTRUCTOR *)item);
}

/****************************************************************************/
/*D
   CreateObject - create a numproc object of a given numproc class

   SYNOPSIS:
   INT CreateObject (MULTIGRID *theMG, char *objectname, char *classname);

   PARAMETERS:
   .  theMG - numproc objects are now associated with a multigrid
   .  objectname - name for new object.
   .  classname - final class name of an existing class.

   DESCRIPTION:
   This function is used by the npcreate command to create instances
   of a class at runtime.
   It creates a new numproc object of a given class. Memory is
   allocated in the environment under the multigrid item. The size is given
   by the constructor structure. Then the constructor function for that class
   is called with the new object.
   CAUTION - name is classname.objectname. Consider NAMESIZE in ugenv.h

   RETURN VALUE:
   INT
   .n    0 no error.
   .n     __LINE__ error, usually memory overflow in environment or something not found.
   D*/
/****************************************************************************/

INT CreateObject (MULTIGRID *theMG, const char *objectname, const char *classname)
{
  NP_CONSTRUCTOR *constructor;
  NP_BASE *object;
  char name[NAMESIZE];

  /* first find constructor */
  constructor = GetConstructor(classname);
  if (constructor==NULL) {
    PrintErrorMessage('E',"CreateObject","cannot find specified class");
    return(__LINE__);
  }

  /* create objects directory in multigrid if necessary */
  if (ChangeEnvDir("/Multigrids") == NULL) return (__LINE__);
  if (ChangeEnvDir(ENVITEM_NAME(theMG)) == NULL) return (__LINE__);
  if (ChangeEnvDir("Objects") == NULL) {
    MakeEnvItem("Objects",ObjectDirID,sizeof(ENVDIR));
    if (ChangeEnvDir("Objects") == NULL) return (__LINE__);
  }
  /* allocate object */
  if (strlen(objectname)+strlen(ENVITEM_NAME(constructor))+2>NAMESIZE)
    return(__LINE__);
  sprintf(name,"%s.%s",ENVITEM_NAME(constructor),objectname);
  object = (NP_BASE *) MakeEnvItem(name,ObjectVarID,constructor->size);
  if (object==NULL) return(__LINE__);

  /* initialize object with constructor */
  object->mg = theMG;
  object->status =  NP_NOT_INIT;
  object->Init = NULL;
  object->Display = NULL;
  object->Execute = NULL;
  if ((*constructor->Construct)(object)!=0) return(__LINE__);

  /* return OK */
  return(0);
}


/****************************************************************************/
/*D
   GetNumProcByName - find a numproc object

   SYNOPSIS:
   NP_BASE *GetNumProcByName (const MULTIGRID *theMG, char *object_name,
   char *abstract_class_name);

   PARAMETERS:
   .  theMG - numproc objects are now associated with a multigrid
   .  object_name - name of object.
   .  abstract_class_name - leading part of an existing class name.

   DESCRIPTION:
   This function is called by numproc objects to find other numproc objects
   of some abstract base class.
   In order to be found the leading characters of the object name must
   coincide with abstract_class_name and the trailing part of the objects
   name must coincide with object_name.

   RETURN VALUE:
   INT
   .n    NULL pointer    object not found
   .n     else pointer to object
   D*/
/****************************************************************************/

NP_BASE *GetNumProcByName (const MULTIGRID *theMG, const char *object_name, const char *abstract_class_name)
{
  ENVITEM *item;
  INT n,m,i;

  /* OK, it is not nice here to access internas of the ENVITEM
     data structure. It would have been better to extent the
     functionality of ugenv instead. But it has been done
     at other places, so why bother her... See
     GetFirstMultigrid, GetFirstVector, etc.
   */

  /* change to object directory of multigrid */
  if (ChangeEnvDir("/Multigrids") == NULL) return (NULL);
  if (ChangeEnvDir(ENVITEM_NAME(theMG)) == NULL) return (NULL);
  item = (ENVITEM *)ChangeEnvDir("Objects");
  if (item == NULL) return (NULL);

  /* now scan through list of objects */
  n = strlen(abstract_class_name);
  for (item=ENVITEM_DOWN(item); item!=NULL; item=NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == ObjectVarID)
    {
      /* check abstract class name */
      if (strncmp(ENVITEM_NAME(item),abstract_class_name,n)!=0) continue;

      /* check object name */
      m = strlen(ENVITEM_NAME(item));
      for (i=m-1; i>=0; i--)
        if (ENVITEM_NAME(item)[i]=='.') break;
      if (strcmp(ENVITEM_NAME(item)+(i+1),object_name)==0) return((NP_BASE *) item);
    }

  /* not found */
  return(NULL);
}

/****************************************************************************/
/*D
   MGListNPClasses - display classes using general print format function

   SYNOPSIS:
   INT MGListNPClasses (const MULTIGRID *theMG, PrintfProcPtr PrintF)

   PARAMETERS:
   .  theMG - numproc objects are now associated with a multigrid
   .  PrintF - general print format function

   DESCRIPTION:
   This function lists all abstract classes of numerical procedures enroled for
   a certain multigrid.

   RETURN VALUE:
   INT
   .n    0 no error.
   .n     __LINE__ error, usually memory overflow in environment or something not found.
   D*/
/****************************************************************************/

INT MGListNPClasses (const MULTIGRID *theMG, PrintfProcPtr PrintF)
{
  ENVITEM *item;
  INT n,i;
  char *p,classlist[MAXCLASSES][NAMESIZE];

  /* change to object directory of multigrid */
  if (ChangeEnvDir("/Multigrids") == NULL) return (__LINE__);
  if (ChangeEnvDir(ENVITEM_NAME(theMG)) == NULL) return (__LINE__);
  item = (ENVITEM *)ChangeEnvDir("Objects");
  if (item == NULL) return (__LINE__);

  /* now scan through list of objects */
  n = 0;
  for (item=ENVITEM_DOWN(item); item!=NULL; item=NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == ObjectVarID)
    {
      if (n>=MAXCLASSES)
        return (__LINE__);

      strcpy(classlist[n],ENVITEM_NAME(item));
      p = strchr(classlist[n],'.');
      ASSERT(p!=NULL);

      *p = '\0';

      for (i=0; i<n; i++)
        if (strcmp(classlist[n],classlist[i])==0)
          break;
      if ((n==0) || (i>=n))
        n++;
    }

  /* list class names */
  for (i=0; i<n; i++)
    PrintF("%s\n",classlist[i]);

  return(0);
}

/****************************************************************************/
/*D
   MGListNPsOfClass - display num proc contents of all num procs of given class

   SYNOPSIS:
   INT MGListNPsOfClass (const MULTIGRID *theMG, const char *ClassName, PrintfProcPtr PrintF)

   PARAMETERS:
   .  theMG - numproc objects are now associated with a multigrid
   .  ClassName - list num procs of this class
   .  PrintF - general print format function

   DESCRIPTION:
   This function displays the num proc contents of all num procs of a given class.

   RETURN VALUE:
   INT
   .n    0 no error.
   .n     __LINE__ error, usually memory overflow in environment or something not found.
   D*/
/****************************************************************************/

INT MGListNPsOfClass (const MULTIGRID *theMG, const char *ClassName, PrintfProcPtr PrintF)
{
  ENVITEM *item;
  INT n;

  /* change to object directory of multigrid */
  if (ChangeEnvDir("/Multigrids") == NULL) return (__LINE__);
  if (ChangeEnvDir(ENVITEM_NAME(theMG)) == NULL) return (__LINE__);
  item = (ENVITEM *)ChangeEnvDir("Objects");
  if (item == NULL) return (__LINE__);

  /* now scan through list of objects */
  n = strlen(ClassName);
  for (item=ENVITEM_DOWN(item); item!=NULL; item=NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == ObjectVarID)
      if (strncmp(ENVITEM_NAME(item),ClassName,n)==0)
      {
        PrintF("\n=================================================\n"
               "contents of num proc '%s':\n",ENVITEM_NAME(item));
        if (((*((NP_BASE*)item)->Display)((NP_BASE*)item)))
          return (__LINE__);
      }

  return(0);
}

/****************************************************************************/
/*D
   MGListAllNPs - display num proc contents of all num procs

   SYNOPSIS:
   INT MGListAllNPs (const MULTIGRID *theMG, PrintfProcPtr PrintF)

   PARAMETERS:
   .  theMG - numproc objects are now associated with a multigrid
   .  PrintF - general print format function

   DESCRIPTION:
   This function displays the num proc contents of all num procs.

   RETURN VALUE:
   INT
   .n    0 no error.
   .n     __LINE__ error, usually memory overflow in environment or something not found.
   D*/
/****************************************************************************/

INT MGListAllNPs (const MULTIGRID *theMG, PrintfProcPtr PrintF)
{
  ENVITEM *item;

  /* change to object directory of multigrid */
  if (ChangeEnvDir("/Multigrids") == NULL) return (__LINE__);
  if (ChangeEnvDir(ENVITEM_NAME(theMG)) == NULL) return (__LINE__);
  item = (ENVITEM *)ChangeEnvDir("Objects");
  if (item == NULL) return (__LINE__);

  /* now scan through list of objects */
  for (item=ENVITEM_DOWN(item); item!=NULL; item=NEXT_ENVITEM(item))
    if (ENVITEM_TYPE(item) == ObjectVarID)
    {
      PrintF("\n=================================================\n"
             "contents of num proc '%s':\n",ENVITEM_NAME(item));
      if (((*((NP_BASE*)item)->Display)((NP_BASE*)item)))
        return (__LINE__);
    }
  return(0);
}

/****************************************************************************/
/*
   InitNumProcManager - Init this file

   SYNOPSIS:
   INT InitNumProcManager (void);

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function inits this file.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    __LINE__ if error occured.
 */
/****************************************************************************/

INT InitNumProcManager ()
{
  ClassDirID  = GetNewEnvDirID();
  ObjectDirID = GetNewEnvDirID();
  ClassVarID  = GetNewEnvVarID();
  ObjectVarID = GetNewEnvVarID();

  return (0);
}
