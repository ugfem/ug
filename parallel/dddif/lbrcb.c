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

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)


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

        #ifdef __THREEDIM__
  /* x and y coordinates are considered to be equal, compare y now */
  if (t1->center[2] < t2->center[2] -SMALL_COORD) return(-1);
  if (t1->center[2] > t2->center[2] +SMALL_COORD) return(1);
        #endif

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

        #ifdef __THREEDIM__
  /* y and x coordinates are considered to be equal, compare x now */
  if (t1->center[2] < t2->center[2] -SMALL_COORD) return(-1);
  if (t1->center[2] > t2->center[2] +SMALL_COORD) return(1);
        #endif

  return(0);
}


#ifdef __THREEDIM__
static int sort_rcb_z (const void *e1, const void *e2)
{
  LB_INFO *t1, *t2;

  t1 = (LB_INFO *)e1;
  t2 = (LB_INFO *)e2;

  if (t1->center[2] < t2->center[2] -SMALL_COORD) return(-1);
  if (t1->center[2] > t2->center[2] +SMALL_COORD) return(1);

  /* z coordinates are considered to be equal, compare x now */
  if (t1->center[1] < t2->center[1] -SMALL_COORD) return(-1);
  if (t1->center[1] > t2->center[1] +SMALL_COORD) return(1);

  /* z and y coordinates are considered to be equal, compare x now */
  if (t1->center[0] < t2->center[0] -SMALL_COORD) return(-1);
  if (t1->center[0] > t2->center[0] +SMALL_COORD) return(1);

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
/*                                                                          */
/* Function:  XferGridWithOverlap                                           */
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

static void XferGridWithOverlap (GRID *theGrid)
{
  ELEMENT *elem;
  NODE *node;

  /* by default, no node stores its vector */
  for(node=FIRSTNODE(theGrid); node!=NULL; node=SUCCN(node))
  {
    SETXFERNODE(node,DEL_VECTOR);
  }


  for(elem=FIRSTELEMENT(theGrid); elem!=NULL; elem=SUCCE(elem))
  {
    int has_local_nb = FALSE;
    int j;

    /* create Master copy */
    DDD_XferCopyObjX(PARHDRE(elem),
                     PARTITION(elem),
                     PrioMaster,
                     (OBJT(elem)==BEOBJ) ? BND_SIZE_TAG(TAG(elem)) : INNER_SIZE_TAG(TAG(elem))
                     );


    /* create 1-overlapping of elements */
    for(j=0; j<SIDES_OF_ELEM(elem); j++)
    {
      ELEMENT *nb = NBELEM(elem,j);

      if (nb!=NULL)
      {
        if (PARTITION(elem)!=PARTITION(nb))
        {
          /* create Ghost copy */
          DDD_XferCopyObjX(PARHDRE(elem),
                           PARTITION(nb),
                           PrioGhost,
                           (OBJT(elem)==BEOBJ) ?
                           BND_SIZE_TAG(TAG(elem)) :
                           INNER_SIZE_TAG(TAG(elem))
                           );
        }

        /* remember any local neighbour element */
        if (PARTITION(nb)==me)
          has_local_nb = TRUE;
      }
    }

    /* consider elements on master-proc */
    if (PARTITION(elem)!=me)
    {
      if (has_local_nb)
      {
        /* element is needed as Ghost copy */
        DDD_PrioritySet(PARHDRE(elem), PrioGhost);
      }
      else
      {
        /* element isn't needed */
        DDD_XferDeleteObj(PARHDRE(elem));
      }
    }
    else
    {
      /* element will be local master */
      int n;
      for(n=0; n<CORNERS_OF_ELEM(elem); n++)
      {
        SETXFERNODE(CORNER(elem,n),KEEP_VECTOR);
      }
    }


    /*
                    {
                            for(i=0; i<3; i++)
                            {
                                    elem->node[i]->elemCnt--;

                                    if (elem->node[i]->triCnt==0)
                                    {
                                            DDD_XferDeleteObj(PARHDR(elem->node[i]));
                                    }
                            }
                    }
     */
  }

  /* evaluate _VECTOR flag for nodes */
  /*
          for(node=FIRSTNODE(theGrid); node!=NULL; node=SUCCN(node))
          {
                  if (XFERNODE(node)==DEL_VECTOR)
                          DDD_XferDeleteObj(PARHDR(NVECTOR(node)));
          }
   */
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



int BalanceGrid (MULTIGRID *theMG)
{
  HEAP *theHeap = theMG->theHeap;
  GRID *theGrid = GRID_ON_LEVEL(theMG,0);       /* balance coarse grid */
  LB_INFO *lbinfo;
  ELEMENT *e;
  int i, son;

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
    theRCB(lbinfo, NT(theGrid), 0, 0, DimX, DimY, 0);

    for (e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
    {
      UserWriteF("elem %08x has dest=%d\n",
                 DDD_InfoGlobalId(PARHDRE(e)), PARTITION(e));
    }
  }


  /* send son elements to father element */
  for (e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
  {
    InheritPartition(e);
  }



  /* start physical transfer */
  ddd_HandlerInit(HSET_XFER);
  DDD_XferBegin();
  if (me==master)
  {
    /* send all grids */
    int g;
    for(g=TOPLEVEL(theMG); g>=0; g--)
    {
      GRID *grid = GRID_ON_LEVEL(theMG,g);
      if (NT(grid)>0) XferGridWithOverlap(grid);
    }
  }

  DDD_XferEnd();


  /* set priorities of border nodes */
  /* TODO this is an extra communication. eventually integrate this
              with grid distribution phase. */
  {
    int g;
    for(g=TOPLEVEL(theMG); g>=0; g--)
    {
      GRID *grid = GRID_ON_LEVEL(theMG,g);
      dddif_SetBorderPriorities(grid);
    }
  }


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
