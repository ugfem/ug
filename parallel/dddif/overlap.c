// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  overlap.c														*/
/*																			*/
/* Purpose:   management of grid overlap during adaption                    */
/*																			*/
/* Author:	  Stefan Lang                                                                           */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*																			*/
/* History:   970204 sl begin                                               */
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

/* standard C library */
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* low module */
#include "compiler.h"
#include "debug.h"
#include "heaps.h"
#include "misc.h"
#include "general.h"

/* dev module */
#include "devices.h"

/* gm module */
#include "algebra.h"
#include "evm.h"
#include "gm.h"
#include "refine.h"
#include "rm.h"
#include "ugm.h"

/* parallel modules */
#include "ppif.h"
#include "ddd.h"
#include "parallel.h"
#include "identify.h"
#include "pargm.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* undefine if overlap should be only updated where needed */
/* This does not work since the connection of the overlap needs the
   fatherelements on both sides (ghost and master sons) and this is
   not ensured.
   #define UPDATE_FULLOVERLAP
 */

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
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*
   UpdateElementOverlap -

   SYNOPSIS:
   static INT UpdateElementOverlap (ELEMENT *theElement);

   PARAMETERS:
   .  theElement

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

static INT UpdateElementOverlap (ELEMENT *theElement)
{
  INT i,s,prio;
  INT SonsOfSide,SonSides[MAX_SONS];
  ELEMENT *theNeighbor,*theSon;
  ELEMENT *SonList[MAX_SONS];

  /* yellow_class specific code:                                */
  /* update need to be done for all elements with THEFLAG set,  */
  /* execpt for yellow copies, since their neighbor need not be */
  /* refined (s.l. 971029)                                      */
#ifndef UPDATE_FULLOVERLAP
  if (!THEFLAG(theElement) && REFINECLASS(theElement)!=YELLOW_CLASS) return(GM_OK);
#endif
  /*
          if (!THEFLAG(theElement)) return(GM_OK);
   */

  for (i=0; i<SIDES_OF_ELEM(theElement); i++)
  {
    theNeighbor = NBELEM(theElement,i);
    if (theNeighbor == NULL) continue;

    prio = EPRIO(theNeighbor);
    if (!IS_REFINED(theNeighbor) || !EHGHOSTPRIO(prio)) continue;

    /* yellow_class specific code:                                     */
    /* this is the special situation an update of the element overlap  */
    /* is needed, since the yellow element has now gotten a new yellow */
    /* neighbor (s.l. 971029)                                          */
    /* sending of yellow copies is now done in each situation. To send */
    /* a yellow copy only if needed, THEFLAG(theNeighbor) must be set  */
    /* properly in AdaptGrid() (980114 s.l.)                          */
                #ifndef UPDATE_FULLOVERLAP
    if ((REFINECLASS(theElement)==YELLOW_CLASS && !THEFLAG(theElement)) &&
        !THEFLAG(theNeighbor)) continue;
                #endif

    PRINTDEBUG(gm,1,("%d: EID=%d side=%d NbID=%d " "NbPARTITION=%d\n",me,
                     ID(theElement),i,ID(theNeighbor), EPROCPRIO(theNeighbor,PrioMaster)))

    Get_Sons_of_ElementSide(theElement,i,&SonsOfSide,
                            SonList,SonSides,1,0);
    PRINTDEBUG(gm,1,("%d: SonsOfSide=%d\n",me,SonsOfSide))

    for (s=0; s<SonsOfSide; s++)
    {
      theSon = SonList[s];
      ASSERT(theSon != NULL);

      PRINTDEBUG(gm,1,("%d: Sending Son=%08x/%x SonID=%d "
                       "SonLevel=%d to dest=%d\n", me,EGID(theSon),theSon,
                       ID(theSon),LEVEL(theSon), EPROCPRIO(theNeighbor,PrioMaster)))

      HEAPFAULT(theNeighbor);

      if (EPROCPRIO(theNeighbor,PrioMaster)>=procs) break;

      XFERECOPYX(theSon,EPROCPRIO(theNeighbor,PrioMaster),PrioHGhost,
                 (OBJT(theSon)==BEOBJ) ? BND_SIZE_TAG(TAG(theSon)) :
                 INNER_SIZE_TAG(TAG(theSon)));
      /* send son to all elements where theNeighbor is master, vghost or vhghost */
      if (0)
      {
        INT *proclist = EPROCLIST(theNeighbor);
        proclist += 2;
        while (*proclist != -1)
        {
          if (!EHGHOSTPRIO(*(proclist+1)))
          {
            XFERECOPYX(theSon,*proclist,PrioHGhost,
                       (OBJT(theSon)==BEOBJ) ? BND_SIZE_TAG(TAG(theSon)) :
                       INNER_SIZE_TAG(TAG(theSon)));
          }
          proclist += 2;
        }
      }
    }
  }

  return(GM_OK);
}

/****************************************************************************/
/*
   UpdateGridOverlap -

   SYNOPSIS:
   INT UpdateGridOverlap (GRID *theGrid);

   PARAMETERS:
   .  theGrid

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

INT UpdateGridOverlap (GRID *theGrid)
{
  ELEMENT *theElement;

  for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    if (IS_REFINED(theElement))
      UpdateElementOverlap(theElement);
  }

  return(GM_OK);
}


/****************************************************************************/
/*
   UpdateMultiGridOverlap -

   SYNOPSIS:
   static INT UpdateMultiGridOverlap (MULTIGRID *theMG, INT FromLevel);

   PARAMETERS:
   .  theMG
   .  FromLevel

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

static INT UpdateMultiGridOverlap (MULTIGRID *theMG, INT FromLevel)
{
  INT l;
  GRID    *theGrid;

  ddd_HandlerInit(HSET_REFINE);

  for (l=FromLevel; l<TOPLEVEL(theMG); l++)
  {
    theGrid = GRID_ON_LEVEL(theMG,l);
    UpdateGridOverlap(theGrid);
  }

  return(GM_OK);
}


/****************************************************************************/
/*
   DropUsedFlags -

   SYNOPSIS:
   static INT DropUsedFlags (GRID *theGrid);

   PARAMETERS:
   .  theGrid

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

static INT DropUsedFlags (GRID *theGrid)
{
  ELEMENT *theElement;

  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    if (USED(theElement) == 1)
    {
      REFINE_ELEMENT_LIST(1,theElement,"drop mark");

      ASSERT(EFATHER(theElement)!=NULL);

      /* this father has to be connected */
      SETUSED(EFATHER(theElement),1);
      SETUSED(theElement,0);
    }
  }

  return(GM_OK);
}


/****************************************************************************/
/*
   ConnectGridOverlap -

   SYNOPSIS:
   INT	ConnectGridOverlap (GRID *theGrid);

   PARAMETERS:
   .  theGrid

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

INT     ConnectGridOverlap (GRID *theGrid)
{
  INT i,j,Sons_of_Side,prio;
  INT SonSides[MAX_SIDE_NODES];
  ELEMENT *theElement;
  ELEMENT *theNeighbor;
  ELEMENT *theSon;
  ELEMENT *Sons_of_Side_List[MAX_SONS];

  for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    prio = EPRIO(theElement);

    /* connect only FROM hgost copies */
    if (!IS_REFINED(theElement) || !EHGHOSTPRIO(prio)) continue;

    PRINTDEBUG(gm,1,("%d: Connecting e=%08x/%x ID=%d eLevel=%d\n",
                     me,DDD_InfoGlobalId(PARHDRE(theElement)),
                     theElement,ID(theElement),
                     LEVEL(theElement)));

    for (i=0; i<SIDES_OF_ELEM(theElement); i++)
    {
      if (OBJT(theElement)==BEOBJ
          && SIDE_ON_BND(theElement,i)
          && !INNER_BOUNDARY(theElement,i)) continue;

      theNeighbor = NBELEM(theElement,i);
      if (theNeighbor == NULL) continue;

      prio = EPRIO(theNeighbor);
      /* overlap situation hasn't changed */
      if (!THEFLAG(theElement) && !THEFLAG(theNeighbor)) continue;

      /* connect only TO master copies */
      if (!IS_REFINED(theNeighbor) || !MASTERPRIO(prio)) continue;

      if (Get_Sons_of_ElementSide(theElement,i,&Sons_of_Side,
                                  Sons_of_Side_List,SonSides,1,0)!=GM_OK) RETURN(GM_FATAL);

      IFDEBUG(gm,1)
      UserWriteF(PFMT "                 side=%d NSONS=%d Sons_of_Side=%d:\n",
                 me,i,NSONS(theElement),Sons_of_Side);
      for (j=0; j<Sons_of_Side; j++)
        UserWriteF(PFMT "            son=%08x/%x sonside=%d\n",
                   me,EGID(Sons_of_Side_List[j]),
                   Sons_of_Side_List[j],SonSides[j]);
      printf("%d:         connecting ghostelements:\n",me);
      ENDDEBUG

      /* the ioflag=1 is needed, since not all sended ghosts are needed! */
      if (Connect_Sons_of_ElementSide(theGrid,theElement,i,
                                      Sons_of_Side,Sons_of_Side_List,SonSides,1)!=GM_OK)
        RETURN(GM_FATAL);
    }

    /* yellow_class specific code:                             */
    /* check whether is a valid ghost, which as in minimum one */
    /* master element as neighbor                              */
    /* TODO: move this functionality to ComputeCopies          */
    /* then disposing of theSon can be done in AdaptGrid      */
    /* and the extra Xfer env around ConnectGridOverlap()      */
    /* can be deleted (s.l. 971029)                            */
    {
      ELEMENT *SonList[MAX_SONS];

      GetAllSons(theElement,SonList);
      for (i=0; SonList[i]!=NULL; i++)
      {
        INT ok = 0;
        theSon = SonList[i];
        if (!EHGHOST(theSon)) continue;
        for (j=0; j<SIDES_OF_ELEM(theSon); j++)
        {
          if (NBELEM(theSon,j)!=NULL && EMASTER(NBELEM(theSon,j))) ok = 1;
        }
        if (!ok)
        {
          if (ECLASS(theSon) == YELLOW_CLASS)
          {
            UserWriteF(PFMT "ConnectGridOverlap(): disposing useless yellow ghost  e=" EID_FMTX
                       "f=" EID_FMTX "this ghost is useless!\n",
                       me,EID_PRTX(theSon),EID_PRTX(theElement));
            DisposeElement(UPGRID(theGrid),theSon,TRUE);
          }
          else
          {
            UserWriteF(PFMT "ConnectGridOverlap(): ERROR e=" EID_FMTX
                       "f=" EID_FMTX "this ghost is useless!\n",
                       me,EID_PRTX(theSon),EID_PRTX(theElement));

            /* TODO: better do this
               assert(0); */
          }
        }
      }
    }
  }

  return(GM_OK);
}



/****************************************************************************/
/*
   ConnectMultiGridOverlap -

   SYNOPSIS:
   static INT	ConnectMultiGridOverlap (MULTIGRID *theMG, INT FromLevel);

   PARAMETERS:
   .  theMG
   .  FromLevel

   DESCRIPTION:

   RETURN VALUE:
   INT
 */
/****************************************************************************/

static INT      ConnectMultiGridOverlap (MULTIGRID *theMG, INT FromLevel)
{
  INT l;
  GRID *theGrid;

  /* drop used marks to fathers */
  for (l=FromLevel+1; l<=TOPLEVEL(theMG); l++)
  {
    theGrid = GRID_ON_LEVEL(theMG,l);
    if (DropUsedFlags(theGrid)) RETURN(GM_FATAL);
  }

  /* connect sons of elements with used flag set */
  for (l=FromLevel; l<TOPLEVEL(theMG); l++)
  {

    theGrid = GRID_ON_LEVEL(theMG,l);
    if (ConnectGridOverlap(theGrid)) RETURN(GM_FATAL);
  }

  return(GM_OK);
}

#endif
