// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      prio.c                                                        */
/*                                                                          */
/* Purpose:   priority management for ddd                                   */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   94/04/25 kb  begin                                            */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

/* standard C library */
#include <stdlib.h>
#include <stdio.h>

#include "dddi.h"


/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/




/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* definition of static variables                                           */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)



/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/



void DDD_PrioritySet (DDD_HDR hdr, DDD_PRIO prio)
{
  TYPE_DESC *desc =  &(theTypeDefs[OBJ_TYPE(hdr)]);

  /* check input parameters */
  if (prio<0 || prio>=MAX_PRIO)
  {
    sprintf(cBuffer,
            "priority must be less than %d in DDD_PrioritySet", MAX_PRIO);
    DDD_PrintError('E', 2305, cBuffer);
    exit(1);
  }


#       ifdef LogObjects
  sprintf(cBuffer, "%4d: LOG DDD_PrioritySet %08x old=%d new=%d\n",
          me, OBJ_GID(hdr), OBJ_PRIO(hdr), prio);
  DDD_PrintDebug(cBuffer);
#       endif

  if (XferActive())
  {
    /* we are during Xfer, therefore initiate PrioChange operation */
    DDD_XferPrioChange(hdr, prio);
  }
  else
  {
    if (! HAS_COUPLING(hdr))
    {
      /* just one local object, we can simply change its priority */
      OBJ_PRIO(hdr) = prio;
    }
    else
    {
      /* distributed object will get inconsistent here. issue warning. */
      if (DDD_GetOption(OPT_WARNING_PRIOCHANGE)==OPT_ON)
      {
        sprintf(cBuffer,
                "creating inconsistency for gid=%08x in DDD_PrioritySet",
                OBJ_GID(hdr));
        DDD_PrintError('W', 2300, cBuffer);
      }

      /* change priority, nevertheless */
      OBJ_PRIO(hdr) = prio;
    }
  }
}



#ifdef F_FRONTEND
void DDD_InfoPriority (DDD_TYPE *type, DDD_OBJ *obj, DDD_PRIO *prio)
{
  *prio = (DDD_PRIO) (OBJ_PRIO(OBJ2HDR(*obj,&(theTypeDefs[*type]))));
}

#endif



/****************************************************************************/

/*
        compute the result of merging two priorities.

        this function merges two priorities p1 and p2 for
        objct of DDD_TYPE desc. the default merge operation
        is the maximum operation. if a merge-matrix has been
        specified, this matrix will be used.

        on return, *pres will contain the resulting priority.
        the return value is:

                PRIO_UNKNOWN   can't decide which priority wins.
                PRIO_FIRST     first priority wins.
                PRIO_SECOND    second priority wins.
                PRIO_ERROR     an error has occurred.
 */

int PriorityMerge (TYPE_DESC *desc, DDD_PRIO p1, DDD_PRIO p2, DDD_PRIO *pres)
{
  if (desc->prioMatrix == NULL)
  {
    if (p1>p2)
    {
      *pres = p1;
      return(PRIO_FIRST);
    }
    if (p2>p1)
    {
      *pres = p2;
      return(PRIO_SECOND);
    }

    *pres = p1;
    return(PRIO_UNKNOWN);
  }
  else
  {
    printf("%4d: allgemeiner prio-merge noch nicht implementiert.\n", me);
    exit(1);
  }

  return(PRIO_ERROR);
}



/****************************************************************************/
