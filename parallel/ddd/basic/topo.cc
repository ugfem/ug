// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      topo.c                                                        */
/*                                                                          */
/* Purpose:   maintains communication structure for data-dependent          */
/*            communication topology                                        */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   93/11/30 kb  begin                                            */
/*            95/10/05 kb  added channel disconnection to TopoExit          */
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
#include <config.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "memmgr.h"

#include "dddi.h"

#include "namespace.h"

USING_UG_NAMESPACES

/* PPIF namespace: */
USING_PPIF_NAMESPACE

START_UGDIM_NAMESPACE

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
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

VChannelPtr *theTopology;


/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)


static DDD_PROC    *theProcArray;
static int         *theProcFlags;


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


/* TODO memory usage is O(P) in current implementation! */

void ddd_TopoInit (void)
{
  int i;

  /* get one channel pointer for each partner */
  theTopology = (VChannelPtr *) AllocFix(procs*sizeof(VChannelPtr));
  if (theTopology==NULL)
  {
    DDD_PrintError('E', 1500, STR_NOMEM " in TopoInit");
    return;
  }

  /* initialize channel topology */
  for(i=0; i<procs; i++)
    theTopology[i] = NULL;


  /* get proc array with maxsize = 2 * number of procs */
  theProcArray = (DDD_PROC *) AllocFix(2 * procs*sizeof(DDD_PROC));
  if (theProcArray==NULL)
  {
    DDD_PrintError('E', 1510, STR_NOMEM " in TopoInit");
    return;
  }

  theProcFlags = (int *) AllocFix(2 * procs*sizeof(int));
  if (theProcFlags==NULL)
  {
    DDD_PrintError('E', 1511, STR_NOMEM " in TopoInit");
    return;
  }
}


void ddd_TopoExit (void)
{
  int i;

  FreeFix(theProcArray);
  FreeFix(theProcFlags);


  /* disconnect channels */
  for(i=0; i<procs; i++)
  {
    if (theTopology[i]!=NULL)
    {
      DiscASync(theTopology[i]);
      while (InfoADisc(theTopology[i])!=1)
        ;
    }
  }

  FreeFix(theTopology);
}


/****************************************************************************/


DDD_PROC *DDD_ProcArray (void)
{
  return theProcArray;
}


RETCODE DDD_GetChannels (int nPartners)
{
  int i, nConn;

  if (nPartners>2*(procs-1))
  {
    DDD_PrintError('E', 1520, "topology error in DDD_GetChannels");
    RET_ON_ERROR;
  }

  nConn = 0;
  for(i=0; i<nPartners; i++)
  {
    if (theTopology[theProcArray[i]]==NULL)
    {
      VChannelPtr vc = ConnASync(theProcArray[i], VC_TOPO);

      if (vc==NULL)
      {
        sprintf(cBuffer,
                "can't connect to proc=%d in DDD_GetChannels",
                theProcArray[i]);
        DDD_PrintError('E', 1521, cBuffer);
        RET_ON_ERROR;
      }

      theTopology[theProcArray[i]] = vc;
      nConn++;

      theProcFlags[i] = TRUE;
    }
    else
    {
      theProcFlags[i] = FALSE;
    }
  }


  while (nConn>0)
  {
    for(i=0; i<nPartners; i++)
    {
      if (theProcFlags[i])
      {
        int ret = InfoAConn(theTopology[theProcArray[i]]);
        if (ret==-1)
        {
          sprintf(cBuffer,
                  "PPIF's InfoAConn() failed for connect to proc=%d"
                  " in DDD_GetChannels",
                  theProcArray[i]);
          DDD_PrintError('E', 1530, cBuffer);
          RET_ON_ERROR;
        }


        if (ret==1)
        {
          theProcFlags[i] = FALSE;
          nConn--;
        }
      }
    }
  }
  /* TODO free unused channels? */

  RET_ON_OK;
}


void DDD_DisplayTopo (void)
{
  int p, i;
  char buf[20];

  DDD_SyncAll();

  if (me==0)
  {
    sprintf(cBuffer, "      ");
    for(p=0; p<procs; p++)
    {
      sprintf(buf, "%2d", p);
      strcat(cBuffer, buf);
    }
    strcat(cBuffer,"\n");
    DDD_PrintLine(cBuffer); fflush(stdout);
  }

  for(p=0; p<procs; p++)
  {
    Synchronize();
    if (p==me)
    {
      sprintf(cBuffer, "%4d: ", me);
      for(i=0; i<procs; i++)
      {
        if (theTopology[i]!=NULL)
        {
          strcat(cBuffer,"<>");
        }
        else
        {
          if (i==p)
            strcat(cBuffer,"--");
          else
            strcat(cBuffer,"  ");
        }
      }
      strcat(cBuffer,"\n");
      DDD_PrintLine(cBuffer);
      DDD_Flush();
    }
  }

  DDD_SyncAll();
}


/****************************************************************************/

END_UGDIM_NAMESPACE
