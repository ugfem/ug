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
RCSID("$Header$",UG_RCS_STRING)

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

INT CreateClass (char *classname, INT size, ConstructorProcPtr Construct)
{
  NP_CONSTRUCTOR *constructor;

  /* create class directory under root if necessary */
  if (ChangeEnvDir("/") == NULL) return (1);
  if (ChangeEnvDir("NumProcClasses") == NULL)
    MakeEnvItem("NumProcClasses",ClassDirID,sizeof(ENVDIR));
  if (ChangeEnvDir("NumProcClasses") == NULL) return (1);

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
   CreateObject - create a numproc object of a given numproc class

   SYNOPSIS:
   INT CreateObject (MULTIGRID *theMG, char *objectname, char *classname);

   PARAMETERS:
   .  theMG - numproc objects are now associated with a multigrid
   .  objectname - name for new object.
   .  classname - name of an existing class.

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
   .n     1 error, usually memory overflow in environment or something not found.
   D*/
/****************************************************************************/

INT CreateObject (MULTIGRID *theMG, char *objectname, char *classname)
{
  NP_CONSTRUCTOR *constructor;
  NP_BASE *object;
  char name[NAMESIZE];

  /* first find constructor */
  if (ChangeEnvDir("/") == NULL) return (1);
  if (ChangeEnvDir("NumProcClasses") == NULL) return(1);
  constructor = (NP_CONSTRUCTOR *) SearchTree(classname,ClassVarID,ClassDirID);
  if (constructor==NULL) return(1);

  /* create objects directory in multigrid if necessary */
  if (ChangeEnvDir("/Multigrids") == NULL) return (1);
  if (ChangeEnvDir(ENVITEM_NAME(theMG)) == NULL) return (1);
  if (ChangeEnvDir("Objects") == NULL)
    MakeEnvItem("Objects",ObjectDirID,sizeof(ENVDIR));
  if (ChangeEnvDir("Objects") == NULL) return (1);

  /* allocate object */
  if (strlen(objectname)+strlen(classname)+2>NAMESIZE) return(1);
  sprintf(name,"%s.%s",classname,objectname);
  object = (NP_BASE *) MakeEnvItem(name,ObjectVarID,constructor->size);
  if (object==NULL) return(1);

  /* initialize object with constructor */
  if ((*constructor->Construct)(object)!=0) return(1);

  /* return OK */
  return(0);
}


/****************************************************************************/
/*D
   GetNumProcByName - find a numproc object

   SYNOPSIS:
   NP_BASE *GetNumProcByName (MULTIGRID *theMG, char *object_name,
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

NP_BASE *GetNumProcByName (MULTIGRID *theMG, char *object_name, char *abstract_class_name)
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
