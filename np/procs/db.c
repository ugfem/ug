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
  NP_DATA_BASE db;

  INT n;
  DOUBLE *list;

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

INT DB_Init (NP_BASE *theNP, INT argc, char **argv)
{
  NP_DATA_BASE *np;

  np = (NP_DATA_BASE *)theNP;
  if (ReadArgvChar("f",np->fname,argc,argv))
    return (NP_NOT_ACTIVE);

  return (NP_ACTIVE);
}

INT DB_Display (NP_BASE *theNP)
{
  NP_DATA_BASE *np;

  np = (NP_DATA_BASE *)theNP;
  UserWriteF(DISPLAY_NP_FORMAT_SS,"f",np->fname);

  return (0);
}

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

static INT ListPreProcess (NP_DATA_BASE *theNP, INT *result)
{
  NP_LIST *np;
  FILE *stream;
  INT i;
  int n;
  double a;

  np = (NP_LIST *)theNP;
        #ifdef ModelP
  if (me == master) {
        #endif
  stream = fileopen(theNP->fname,"r");
  if (stream == NULL)
    NP_RETURN(1,*result);
  if (fscanf(stream,"%d",&n) != 1) {
    fclose(stream);
    NP_RETURN(1,*result);
  }
        #ifdef ModelP
}
Broadcast(&n,sizeof(int));
        #endif
  np->n = n;
  if (SetStringValue(":DB:size",(DOUBLE) n))
    NP_RETURN(1,*result);
  np->list = (DOUBLE*) GetMemoryForObject(NP_MG(np),n*sizeof(DOUBLE),-1);
        #ifdef ModelP
  if (me == master) {
        #endif
  for (i=0; i<n; i++) {
    if (fscanf(stream,"%lf",&a) != 1) {
      fclose(stream);
      NP_RETURN(1,*result);
    }
    np->list[i] = a;
  }
  fclose(stream);
        #ifdef ModelP
}
Broadcast(np->list,n*sizeof(DOUBLE));
        #endif

  return(0);
}

static INT ListPostProcess (NP_DATA_BASE *theNP, INT *result)
{
  NP_LIST *np;

  np = (NP_LIST *)theNP;
  if (PutFreelistMemory(MGHEAP(NP_MG(np)),np->list,np->n*sizeof(DOUBLE)))
    NP_RETURN(1,*result);

  return(0);
}

INT ListSize (NP_DATA_BASE *theNP, INT *n, INT *result)
{
  NP_LIST *np;

  np = (NP_LIST *)theNP;
  *n = np->n;

  return(0);
}

INT ListData (NP_DATA_BASE *theNP, INT i, void **data, INT *result)
{
  NP_LIST *np;

  np = (NP_LIST *)theNP;
  data[0] = (void *) (np->list+i);

  return(0);
}

static INT List_Construct (NP_BASE *theNP)
{
  NP_DATA_BASE *np;

  theNP->Init = DB_Init;
  theNP->Display = DB_Display;
  theNP->Execute = NULL;

  np = (NP_DATA_BASE *)theNP;
  np->PreProcess = ListPreProcess;
  np->GetSize = ListSize;
  np->SetSize = NULL;
  np->GetData = ListData;
  np->SetData = NULL;
  np->PostProcess = ListPostProcess;

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

  if (CreateClass(DATA_BASE_CLASS_NAME ".list",sizeof(NP_LIST),List_Construct))
    return (__LINE__);

  return (0);
}
