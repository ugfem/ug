// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  lbrcb.c														*/
/*																			*/
/* Purpose:   simple static load balancing scheme for testing initial		*/
/*            grid distribution												*/
/*																			*/
/* Author:	  Klaus Birken                                                                          */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: birken@ica3.uni-stuttgart.de							*/
/*																			*/
/* History:   940416 kb  begin                                                                          */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#ifndef __LBRCB_H__
#define __LBRCB_H__

#ifdef ModelP

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/


#include "parallel.h"
#include "evm.h"



/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/


#define SMALL_COORD         1.0E-5      /* resolution when comparing COORDs */



/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct {
  ELEMENT *elem;
  COORD center[DIM];
} LB_INFO;


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

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/



static int sort_rcb_x (const void *e1, const void *e2)
{
  LB_INFO *t1, *t2;

  t1 = (LB_INFO *)e1;
  t2 = (LB_INFO *)e2;

  if (t1->center[0] < t2->center[0] -SMALL_COORD) return(-1);
  if (t1->center[0] > t2->center[0] +SMALL_COORD) return(1);

  /* x coordinates are considered to be equal, compare y now */
  if (t1->center[1] < t2->center[1] -SMALL_COORD) return(-1);
  if (t1->center[1] > t2->center[1] +SMALL_COORD) return(1);

  return(0);
}


static int sort_rcb_y (const void *e1, const void *e2)
{
  LB_INFO *t1, *t2;

  t1 = (LB_INFO *)e1;
  t2 = (LB_INFO *)e2;

  if (t1->center[1] < t2->center[1] -SMALL_COORD) return(-1);
  if (t1->center[1] > t2->center[1] +SMALL_COORD) return(1);

  /* y coordinates are considered to be equal, compare x now */
  if (t1->center[0] < t2->center[0] -SMALL_COORD) return(-1);
  if (t1->center[0] > t2->center[0] +SMALL_COORD) return(1);

  return(0);
}



/****************************************************************************/
/*                                                                          */
/* Function:  theRCB                                                        */
/*                                                                          */
/* Purpose:   simple load balancing algorithm,                              */
/*            balances all local triangles using                            */
/*            a 'recursive coordinate bisection' scheme,                    */
/*                                                                          */
/* Input:     theItems:  LB_INFO array                                      */
/*            nItems:    length of array                                    */
/*            px,py:     bottom left position in 2D processor array         */
/*            dx,dy:     size of 2D processor array                         */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

static void theRCB (LB_INFO *theItems, int nItems, int px, int py, int dx, int dy)
{
  int i, part0, part1, ni0, ni1;

  if (nItems==0)
    return;

  if ((dx<=1)&&(dy<=1))
  {
    for(i=0; i<nItems; i++)
    {
      int dest = py*DimX+px;
      PARTITION(theItems[i].elem) = dest;
    }
    return;
  }

  if (dx>=dy)
  {
    if (nItems>1) qsort(theItems, nItems, sizeof(LB_INFO), sort_rcb_x);

    part0 = dx/2;
    part1 = dx-part0;

    ni0 = (int)(((double)part0)/((double)(dx))*((double)nItems));
    ni1 = nItems-ni0;

    theRCB(theItems,     ni0, px,       py, part0, dy);
    theRCB(theItems+ni0, ni1, px+part0, py, part1, dy);

  }
  else
  {
    if (nItems>1) qsort(theItems, nItems, sizeof(LB_INFO), sort_rcb_y);

    part0 = dy/2;
    part1 = dy-part0;

    ni0 = (int)(((double)part0)/((double)(dy))*((double)nItems));
    ni1 = nItems-ni0;

    theRCB(theItems,     ni0, px, py      , dx, part0);
    theRCB(theItems+ni0, ni1, px, py+part0, dx, part1);
  }
}




/****************************************************************************/
/*                                                                          */
/* Function:  XferElemsAndOverlap                                           */
/*                                                                          */
/* Purpose:   send elements to other procs, keep overlapping region of one  */
/*            element, maintain correct priorities at interfaces.           */
/*                                                                          */
/*            the destination procs have been computed by theRCB function   */
/*            and put into the elements' PARTITION-entries.                 */
/*                                                                          */
/*            NOTE: this routine handles only the local grid case, i.e.,    */
/*                  distributed grids cannot be redistributed               */
/*                                                                          */
/* Input:     -                                                             */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

static void XferElemsAndOverlap (GRID *theGrid)
{
  ELEMENT *elem, *nb;
  int i;


  /* create element copies */
  for(elem=FIRSTELEMENT(theGrid); elem!=NULL; elem=SUCCE(elem))
  {
    /* create element copy */
    DDD_XferCopyObjX(PARHDRE(elem),
                     PARTITION(elem),
                     0,
                     (OBJT(elem)==BEOBJ) ? BND_SIZE(TAG(elem)) : INNER_SIZE(TAG(elem))
                     );

    if (PARTITION(elem)!=me)
    {
      /*DDD_XferDeleteObj(PARHDRE(elem)); */
      /*
                              for(i=0; i<3; i++)
                              {
                                      elem->node[i]->elemCnt--;

                                      if (elem->node[i]->triCnt==0)
                                      {
                                              DDD_XferDeleteObj(PARHDR(elem->node[i]));
                                      }
                              }
       */
    }
  }
}




/****************************************************************************/

static void CenterOfMass (ELEMENT *e, COORD *pos)
{
  int i;

  V_DIM_CLEAR(pos)

  for(i=0; i<CORNERS_OF_ELEM(e); i++)
  {
    V_DIM_LINCOMB(1.0,pos,1.0,CVECT(MYVERTEX(CORNER(e,i))),pos)
  }

  V_DIM_SCALE(1.0/(float)CORNERS_OF_ELEM(e),pos)
}



/****************************************************************************/


int BalanceGrid (MULTIGRID *theMG)
{
  HEAP *theHeap = theMG->theHeap;
  GRID *theGrid = GRID_ON_LEVEL(theMG,CURRENTLEVEL(theMG));
  LB_INFO *lbinfo;
  ELEMENT *e;
  int i;

  /* distributed grids cannot be redistributed by this function */

  if (me==master)
  {
    Mark(theHeap,FROM_TOP);
    lbinfo = (LB_INFO *)
             GetMem(theHeap, NT(theGrid)*sizeof(LB_INFO), FROM_TOP);

    if (lbinfo==NULL)
    {
      Release(theHeap,FROM_TOP);
      UserWrite("ERROR: could not allocate memory from the MGHeap\n");
      return (1);
    }


    /* construct LB_INFO list */
    for (i=0, e=FIRSTELEMENT(theGrid); e!=NULL; i++, e=SUCCE(e))
    {
      lbinfo[i].elem = e;
      CenterOfMass(e, lbinfo[i].center);
    }


    /* apply coordinate bisection strategy */
    theRCB(lbinfo, NT(theGrid), 0, 0, DimX, DimY);

    for (e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
    {
      UserWriteF("elem %08x has dest=%d\n",
                 DDD_InfoGlobalId(PARHDRE(e)), PARTITION(e));
    }
  }


  /* start physical transfer */
  DDD_XferBegin();
  if (me==master && NT(theGrid)>0) XferElemsAndOverlap(theGrid);
  DDD_XferEnd();


  if (me==master)
  {
    Release(theHeap,FROM_TOP);
  }

  DDD_ConsCheck();

  return 0;
}




/****************************************************************************/

#endif  /* ModelP */
#endif
