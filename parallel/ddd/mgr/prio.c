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
