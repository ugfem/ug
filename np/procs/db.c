// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  db.c                                                                                                          */
/*																			*/
/* Purpose:   data base interface                                                                       */
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   Sep 11, 1997 begin                                                                */
/*																			*/
/* Remarks:   not finished!                                                                     */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <string.h>

#include "devices.h"
#include "ugenv.h"

#include "scan.h"
#include "numproc.h"
#include "np.h"
#include "ugm.h"
#include "general.h"
#include "fileopen.h"
#include "ugstruct.h"

#include "db.h"

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

typedef struct
{
  NP_ORDERED_LIST db;

  char name[NAMELEN];
  INT n;
  DOUBLE *list;
  DOUBLE regular_step;

} NP_LIST;

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

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   list - numproc list

   DESCRIPTION:
   Reads a file with double values.

   'npinit <name> $f <file>'

   .  <name> - num proc name
   .  $f~<file> - file name
   D*/
/****************************************************************************/

static int cmp_real (const void *p, const void *q)
{
  DOUBLE v1,v2;

  v1 = *((DOUBLE*)p);
  v2 = *((DOUBLE*)q);
  if (v1<v2) return (-1);
  if (v2<v1) return (1);
  return (0);
}

static INT List_PreProcess (NP_ORDERED_LIST *theNP, INT *result)
{
  NP_LIST *np;
  int i,n;
  char buffer[NAMESIZE];

  np = (NP_LIST *)theNP;
  n = np->n;
  np->list = (DOUBLE*) GetMemoryForObject(NP_MG(np),n*sizeof(DOUBLE),-1);
  for (i=0; i<n; i++)
  {
    sprintf(buffer,"%s%d",np->name,(int)i);
    if (GetStringValue(buffer,np->list+i)) return (1);
  }

  /* sort list */
  qsort((void *)np->list,np->n,sizeof(DOUBLE),cmp_real);

  return(0);
}

static INT List_PostProcess (NP_ORDERED_LIST *theNP, INT *result)
{
  NP_LIST *np;

  np = (NP_LIST *)theNP;
  if (PutFreelistMemory(MGHEAP(NP_MG(np)),np->list,np->n*sizeof(DOUBLE)))
    NP_RETURN(1,*result);

  return(0);
}

static INT List_GetListEntry_Index (NP_ORDERED_LIST *theNP, INT n, DOUBLE *Entry, INT *result)
{
  NP_LIST *np;

  np = (NP_LIST *)theNP;
  if (n>=np->n)
  {
    *Entry = 0.0;
    *result = 0;
  }
  else
  {
    *Entry = np->list[n];
    *result = 1;
  }

  return(0);
}

static INT List_GetListEntry_NextHigherEntry (NP_ORDERED_LIST *theNP, DOUBLE value, DOUBLE *Entry, INT *result)
{
  NP_LIST *np;
  INT i;
  DOUBLE tp;

  np = (NP_LIST *)theNP;
  *result = 0;
  for (i=0; i<np->n; i++)
    if (np->list[i]>value)
    {
      *Entry = np->list[i];
      *result = 1;
      break;
    }
  if (np->regular_step>0.0)
  {
    tp = np->regular_step*(floor(value/np->regular_step)+1.0);
    if (*result==1)
    {
      *Entry = MIN(*Entry,tp);
    }
    else
    {
      *result = 1;
      *Entry = tp;
    }
  }

  return(0);
}

INT List_Init (NP_BASE *theNP, INT argc, char **argv)
{
  NP_LIST *np;

  np = (NP_LIST *)theNP;
  if (ReadArgvINT("n",&(np->n),argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvChar ("L",np->name,argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvDOUBLE("s",&(np->regular_step),argc,argv))
  {
    np->regular_step = -1.0;
  }

  return (NP_ACTIVE);
}

INT List_Display (NP_BASE *theNP)
{
  NP_LIST *np;
  INT i;
  char buffer[16];

  np = (NP_LIST *)theNP;
  UserWriteF(DISPLAY_NP_FORMAT_SI,"n",(int)np->n);
  if (np->list!=NULL)
    for (i=0; i<np->n; i++)
    {
      sprintf(buffer,"List[%d]",(int)i);
      UserWriteF(DISPLAY_NP_FORMAT_SF,buffer,np->list[i]);
    }

  return (0);
}

/****************************************************************************/
static INT List_Construct (NP_BASE *theNP)
{
  NP_ORDERED_LIST *np;

  theNP->Init = List_Init;
  theNP->Display = List_Display;
  theNP->Execute = NULL;

  np = (NP_ORDERED_LIST *)theNP;
  np->PreProcess                                          = List_PreProcess;
  np->GetListEntry_Index                          = List_GetListEntry_Index;
  np->GetListEntry_NextHigherEntry        = List_GetListEntry_NextHigherEntry;
  np->PostProcess                                         = List_PostProcess;

  return(0);
}

/****************************************************************************/
/*D
   InitDb - Enrol data base

   SYNOPSIS:
   INT InitDb (void);

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function creates the numproc 'list'.
   It is called in initnp.c.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT InitDb (void)
{
  if (MakeStruct(":DB")!=0) return (__LINE__);

  if (CreateClass(ORDERED_LIST_CLASS_NAME ".list",sizeof(NP_LIST),List_Construct))
    return (__LINE__);

  return (0);
}
