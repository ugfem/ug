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
#include "general.h"
#include "ugm.h"
#include "devices.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/


#define SMALL_DOUBLE         1.0E-5      /* resolution when comparing DOUBLEs */



/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct {
  ELEMENT *elem;
  DOUBLE center[DIM];
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

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


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

  if (t1->center[0] < t2->center[0] -SMALL_DOUBLE) return(-1);
  if (t1->center[0] > t2->center[0] +SMALL_DOUBLE) return(1);

  /* x coordinates are considered to be equal, compare y now */
  if (t1->center[1] < t2->center[1] -SMALL_DOUBLE) return(-1);
  if (t1->center[1] > t2->center[1] +SMALL_DOUBLE) return(1);

        #ifdef __THREEDIM__
  /* x and y coordinates are considered to be equal, compare y now */
  if (t1->center[2] < t2->center[2] -SMALL_DOUBLE) return(-1);
  if (t1->center[2] > t2->center[2] +SMALL_DOUBLE) return(1);
        #endif

  return(0);
}


static int sort_rcb_y (const void *e1, const void *e2)
{
  LB_INFO *t1, *t2;

  t1 = (LB_INFO *)e1;
  t2 = (LB_INFO *)e2;

  if (t1->center[1] < t2->center[1] -SMALL_DOUBLE) return(-1);
  if (t1->center[1] > t2->center[1] +SMALL_DOUBLE) return(1);

  /* y coordinates are considered to be equal, compare x now */
  if (t1->center[0] < t2->center[0] -SMALL_DOUBLE) return(-1);
  if (t1->center[0] > t2->center[0] +SMALL_DOUBLE) return(1);

        #ifdef __THREEDIM__
  /* y and x coordinates are considered to be equal, compare x now */
  if (t1->center[2] < t2->center[2] -SMALL_DOUBLE) return(-1);
  if (t1->center[2] > t2->center[2] +SMALL_DOUBLE) return(1);
        #endif

  return(0);
}


#ifdef __THREEDIM__
static int sort_rcb_z (const void *e1, const void *e2)
{
  LB_INFO *t1, *t2;

  t1 = (LB_INFO *)e1;
  t2 = (LB_INFO *)e2;

  if (t1->center[2] < t2->center[2] -SMALL_DOUBLE) return(-1);
  if (t1->center[2] > t2->center[2] +SMALL_DOUBLE) return(1);

  /* z coordinates are considered to be equal, compare x now */
  if (t1->center[1] < t2->center[1] -SMALL_DOUBLE) return(-1);
  if (t1->center[1] > t2->center[1] +SMALL_DOUBLE) return(1);

  /* z and y coordinates are considered to be equal, compare x now */
  if (t1->center[0] < t2->center[0] -SMALL_DOUBLE) return(-1);
  if (t1->center[0] > t2->center[0] +SMALL_DOUBLE) return(1);

  return(0);
}
#endif

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
/*            dim:       sort dimension 0=x, 1=y, 2=z                       */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

static void theRCB (LB_INFO *theItems, int nItems, int px, int py, int dx, int dy, int dim)
{
  int i, part0, part1, ni0, ni1;
  static int (*sort_function)(const void *e1, const void *e2);

  /* determine sort function */
  switch (dim) {
  case 0 :
    sort_function = sort_rcb_x;
    break;
  case 1 :
    sort_function = sort_rcb_y;
    break;
                #ifdef __THREEDIM__
  case 2 :
    sort_function = sort_rcb_z;
    break;
                #endif
  default :
    printf("%d: theRCB(): ERROR no valid sort dimension specified\n",me);
    break;
  }

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
    if (nItems>1) qsort(theItems, nItems, sizeof(LB_INFO), sort_function);

    part0 = dx/2;
    part1 = dx-part0;

    ni0 = (int)(((double)part0)/((double)(dx))*((double)nItems));
    ni1 = nItems-ni0;

    theRCB(theItems,     ni0, px,       py, part0, dy,(dim+1)%DIM);
    theRCB(theItems+ni0, ni1, px+part0, py, part1, dy,(dim+1)%DIM);

  }
  else
  {
    if (nItems>1) qsort(theItems, nItems, sizeof(LB_INFO), sort_function);

    part0 = dy/2;
    part1 = dy-part0;

    ni0 = (int)(((double)part0)/((double)(dy))*((double)nItems));
    ni1 = nItems-ni0;

    theRCB(theItems,     ni0, px, py      , dx, part0,(dim+1)%DIM);
    theRCB(theItems+ni0, ni1, px, py+part0, dx, part1,(dim+1)%DIM);
  }
}





/****************************************************************************/

static void CenterOfMass (ELEMENT *e, DOUBLE *pos)
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


static void InheritPartition (ELEMENT *e)
{
  int i;

  for(i=0; i<SONS_OF_ELEM(e); i++)
  {
    ELEMENT *son = SON(e,i);
    if (son==NULL) break;

    PARTITION(son) = PARTITION(e);
    InheritPartition(son);
  }
}



int BalanceGridRCB (MULTIGRID *theMG, int level)
{
  HEAP *theHeap = theMG->theHeap;
  GRID *theGrid = GRID_ON_LEVEL(theMG,level);       /* balance grid of level */
  LB_INFO *lbinfo;
  ELEMENT *e;
  int i, son;

  /* distributed grids cannot be redistributed by this function */

  if (me==master)
  {
    if (NT(theGrid) == 0)
    {
      UserWriteF("WARNING in BalanceGridRCB: no elements in grid\n");
      return (1);
    }

    Mark(theHeap,FROM_TOP);
    lbinfo = (LB_INFO *)
             GetMem(theHeap, NT(theGrid)*sizeof(LB_INFO), FROM_TOP);

    if (lbinfo==NULL)
    {
      Release(theHeap,FROM_TOP);
      UserWrite("ERROR in BalanceGridRCB: could not allocate memory from the MGHeap\n");
      return (1);
    }


    /* construct LB_INFO list */
    for (i=0, e=FIRSTELEMENT(theGrid); e!=NULL; i++, e=SUCCE(e))
    {
      lbinfo[i].elem = e;
      CenterOfMass(e, lbinfo[i].center);
    }


    /* apply coordinate bisection strategy */
    theRCB(lbinfo, NT(theGrid), 0, 0, DimX, DimY, 0);

    IFDEBUG(dddif,1)
    for (e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
    {
      UserWriteF("elem %08x has dest=%d\n",
                 DDD_InfoGlobalId(PARHDRE(e)), PARTITION(e));
    }
    ENDDEBUG

    Release(theHeap,FROM_TOP);
  }

  return 0;
}


/****************************************************************************/

#endif  /* ModelP */
