// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  gridcons.c													*/
/*																			*/
/* Purpose:   basic functions for managing consistency of distributed grids */
/*																			*/
/* Author:	  Stefan Lang, Klaus Birken										*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: birken@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   960906 kb  begin                                                                                          */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#ifdef ModelP

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdlib.h>

#include "debug.h"
#include "parallel.h"
#include "general.h"


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

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)


/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/*
        for all PrioMaster-nodes with remote copies, set exactly one
        to PrioMaster, the other copies to PrioBorder in order to establish
        the BorderNodeIF. this is done for one grid.
 */


static int dddif_ComputeNodeBorderPrios (DDD_OBJ obj)
{
  NODE    *node  = (NODE *)obj;
  int     *plist = DDD_InfoProcList(PARHDR(node));
  int i, min_proc = procs;

  /*
          minimum processor number will get Master-node,
          all others get Border-nodes
   */
  for(i=0; plist[i]>=0; i+=2)
  {
    if (plist[i+1]==PrioMaster && plist[i]<min_proc)
      min_proc = plist[i];
  }

  if (min_proc==procs)
    return;

  if (me!=min_proc)
    DDD_PrioritySet(PARHDR(node), PrioBorder);
}

static int dddif_ComputeVectorBorderPrios (DDD_OBJ obj)
{
  VECTOR  *vector  = (VECTOR *)obj;
  int     *plist = DDD_InfoProcList(PARHDR(vector));
  int i, min_proc = procs;

  /*
          minimum processor number will get Master-node,
          all others get Border-nodes
   */
  for(i=0; plist[i]>=0; i+=2)
  {
    if (plist[i+1]==PrioMaster && plist[i]<min_proc)
      min_proc = plist[i];
  }

  if (min_proc==procs)
    return;

  if (me!=min_proc)
    DDD_PrioritySet(PARHDR(vector), PrioBorder);
}


void dddif_SetBorderPriorities (GRID *theGrid)
{
  DDD_XferBegin();
  DDD_IFAExecLocal(BorderNodeSymmIF, GLEVEL(theGrid),
                   dddif_ComputeNodeBorderPrios);
  DDD_IFAExecLocal(BorderVectorSymmIF, GLEVEL(theGrid),
                   dddif_ComputeVectorBorderPrios);
  DDD_XferEnd();
}


/****************************************************************************/

#endif
