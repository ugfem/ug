// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  gridcons.c													*/
/*																			*/
/* Purpose:   functions for managing consistency of distributed grids       */
/*																			*/
/* Author:	  Stefan Lang, Klaus Birken										*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
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

#include <config.h>
#include <stdlib.h>

#include "debug.h"
#include "parallel.h"
#include "general.h"
#include "gm.h"
#include "refine.h"
#include "ugm.h"
#include "evm.h"
#include "shapes.h"
#include "ugdevices.h"
#include "namespace.h"

/* UG namespaces: */
USING_UG_NAMESPACES

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#ifdef DDD_PRIO_ENV
#define DDD_XferBegin   DDD_PrioBegin
#define DDD_XferEnd             DDD_PrioEnd
#endif


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

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*
   ConstructConsistentGridLevel -

   SYNOPSIS:
   void ConstructConsistentGridLevel (GRID *theGrid);

   PARAMETERS:
   .  theGrid

   DESCRIPTION:
   Provide a consistent grid level. Do not call for a single grid level, but
   always for the whole multigrid from bottom to top, since otherwise
   consistency is not ensured!

   RETURN VALUE:
   void
 */
/****************************************************************************/

static void ConstructConsistentGridLevel (GRID *theGrid)
{
  INT i,j,k,l,m,o;
  DOUBLE fac,*local;
  ELEMENT *theElement,*theFather,*theNb;
  NODE    *theNode;
  EDGE    *theEdge;
  VERTEX  *theVertex;

#ifdef __PERIODIC_BOUNDARY__
  SetPerVecVOBJECT(theGrid);
#endif

  /* this is the simplest fix for VFATHER zombies  */
  /* just reset all VFATHER pointers and set them  */
  /* only by master nodes of this or upper levels. */
  /* A more complicated fix would be to set the    */
  /* priorities of the vertices correctly.         */
  /* (980126 s.l.)                                 */
  for (theVertex = PFIRSTVERTEX(theGrid); theVertex != NULL;
       theVertex = SUCCV(theVertex)) {
    VFATHER(theVertex) = NULL;
    /*
                if (VXGHOST(theVertex))
                        VFATHER(theVertex) = NULL;
                    else if (OBJT(theVertex) == BVOBJ)
                        if (MOVED(theVertex))
                            {
                                INT n;
                                    DOUBLE *x[MAX_CORNERS_OF_ELEM];

                                    theElement = VFATHER(theVertex);
                                    if (theElement == NULL) continue;
                                    HEAPFAULT(theElement);
                                    CORNER_COORDINATES(theElement,n,x);
                                    UG_GlobalToLocal(n,(const DOUBLE **)x,
                                                                     CVECT(theVertex),LCVECT(theVertex));
                            }
     */
  }

  /* reconstruct VFATHER pointers and                */
  /* make ghost neighborships symmetric (only for 3d)*/
  /*
   * Also: If there are SideVectors, set the VCOUNT field correctly (i.e., the number
   * of elements that reference a given SideVector).  This information has been transferred as-is
   * during load balancing, but may be wrong on new ghost elements.
   */
  for (theElement = PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
  {
    /* This is the SideVector part */
#ifdef __THREEDIM__
    if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
        SETVCOUNT(SVECTOR(theElement,i), (NBELEM(theElement,i) ? 2 : 1));
#endif

    /* Here comes the rest, all the way to the end of the element loop */

    /* TODO: delete now done in ElementObjMkCons()
       #ifdef __THREEDIM__
                    if (EVGHOST(theElement))
                    {
                            ELEMENT *NbElement;

                            for (i=0; i<SIDES_OF_ELEM(theElement); i++)
                            {
                                    NbElement = NBELEM(theElement,i);
                                    for (j=0; j<SIDES_OF_ELEM(NbElement); j++)
                                    {
                                            if (NBELEM(NbElement,j) == theElement) break;
                                    }
                                    if (j>=SIDES_OF_ELEM(NbElement))
                                            SET_NBELEM(theElement,i,NULL);
                            }
                    }
       #endif
     */

    theFather = EFATHER(theElement);

    /* no reconstruction of VFATHER possible */
    if (theFather == NULL) continue;

    for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
    {
      theNode = CORNER(theElement,i);
      if (CORNERTYPE(theNode)) continue;

      theVertex = MYVERTEX(theNode);

      /* this is too few for arbitrary load balancing, since
              VFATHER pointer may have changed (970828 s.l.)
                              if (VFATHER(theVertex)==NULL || EPRIO(VFATHER(theVertex))==PrioHGhost)
       */
      /* this is too few for arbitrary load balancing, since
              VFATHER pointer may be already a zombie pointer (980126 s.l.)
                              if (VFATHER(theVertex)==NULL || EPRIO(theFather)!=PrioHGhost)
       */
      {
        switch (NTYPE(theNode))
        {
        case (MID_NODE) :
        {
          INT co0,co1;

          for (j=0; j<EDGES_OF_ELEM(theFather); j++)
          {
            theEdge = GetEdge(CORNER(theFather,CORNER_OF_EDGE(theFather,j,0)),
                              CORNER(theFather,CORNER_OF_EDGE(theFather,j,1)));
            if (MIDNODE(theEdge) == theNode) break;
          }
          /* here should be an assertion, but not in each situation the
             midnode pointer is set (970829 s.l.)
                                                          ASSERT(j<EDGES_OF_ELEM(theFather));
           */
          if (j>=EDGES_OF_ELEM(theFather))
          {
            for (j=0; j<EDGES_OF_ELEM(theFather); j++)
            {
              theEdge = GetEdge(CORNER(theFather,CORNER_OF_EDGE(theFather,j,0)),
                                CORNER(theFather,CORNER_OF_EDGE(theFather,j,1)));
                                                                #ifdef Debug
              if (theEdge->midnode != NULL)
                PRINTDEBUG(dddif,1,
                           (PFMT " ConstructConsistentGrid(): elem=" EID_FMTX
                            " i=%d n1=" ID_FMTX " n2=" ID_FMTX " midnode= " ID_FMTX  "\n",
                            me,theFather,EID_PRTX(theFather),j,
                            ID_PRTX(NBNODE(LINK0(theEdge))),
                            ID_PRTX(NBNODE(LINK1(theEdge))),
                            ID_PRTX(theEdge->midnode)))
                                                                #endif
            }



            PRINTDEBUG(dddif,1,
                       ("ConstructConsistentGrid(): WARN "
                        " theNode= " ID_FMTX
                        " vertex= " VID_FMTX
                        " recalculation of VFATHER impossible\n",
                        ID_PRTX(NBNODE(LINK0(theEdge))),
                        VID_PRTX(theVertex)));
            /* if it couldn't  be recalculated reset it */
            VFATHER(theVertex) = NULL;
            break;
          }


          /* reconstruct local coordinates of vertex */
          co0 = CORNER_OF_EDGE(theFather,j,0);
          co1 = CORNER_OF_EDGE(theFather,j,1);

          /* local coordinates have to be local towards pe */


          V_DIM_LINCOMB(0.5, LOCAL_COORD_OF_ELEM(theFather,co0),
                        0.5, LOCAL_COORD_OF_ELEM(theFather,co1),
                        LCVECT(theVertex));
          SETONEDGE(theVertex,j);

          break;
        }

                                        #ifdef __THREEDIM__
        case (SIDE_NODE) :
          /* always compute new coords for this case! */
          if (TAG(theFather) == PYRAMID) k=0;
          else
            k =  GetSideIDFromScratch(theElement,theNode);
          ASSERT(k < SIDES_OF_ELEM(theFather));

          SETONSIDE(theVertex,k);

          m = CORNERS_OF_SIDE(theFather,k);
          local = LCVECT(theVertex);
          fac = 1.0 / m;
          V_DIM_CLEAR(local);
          for (o=0; o<m; o++)
          {
            l = CORNER_OF_SIDE(theFather,k,o);
            V_DIM_LINCOMB(1.0,local,1.0,
                          LOCAL_COORD_OF_ELEM(theFather,l),local);
          }
          V_DIM_SCALE(fac,local);

          theNb = NBELEM(theFather,k);
          if (theNb != NULL)
          {
            for (j=0; j<SIDES_OF_ELEM(theNb); j++)
            {
              if (NBELEM(theNb,j) == theFather) break;
            }
            ASSERT(j < SIDES_OF_ELEM(theNb));
            SETONNBSIDE(theVertex,j);
          }
          else SETONNBSIDE(theVertex,MAX_SIDES_OF_ELEM);
          break;

                                        #endif
        case (CENTER_NODE) :
        case (LEVEL_0_NODE) :
          /* nothing to do */
          break;

        case (CORNER_NODE) :
        default :
          assert(0);
          break;
        }
                                #ifdef Debug
        if (theFather != NULL) HEAPFAULT(theFather);
                                #endif
        VFATHER(theVertex) = theFather;

        if (OBJT(theVertex) == BVOBJ)
          if (MOVED(theVertex)) {
            INT n;
            DOUBLE *x[MAX_CORNERS_OF_ELEM];

            CORNER_COORDINATES(theFather,n,x);
            UG_GlobalToLocal(n,(const DOUBLE **)x,
                             CVECT(theVertex),LCVECT(theVertex));
          }
      }
    }
  }
}

/****************************************************************************/
/*
   ConstructConsistentGrid -

   SYNOPSIS:
   void ConstructConsistentGrid (GRID *theGrid);

   PARAMETERS:
   .  theGrid

   DESCRIPTION:
   Provide a consistent grid level. Do not call for a single grid level, but
   always for the whole multigrid from bottom to top, since otherwise
   consistency is not ensured!

   RETURN VALUE:
   void
 */
/****************************************************************************/

void NS_DIM_PREFIX ConstructConsistentGrid (GRID *theGrid)
{
  /* the setting of the priorities has to be done in two waves after */
  /* completion of the grid transfer, since                          */
  /* - decisions about vghost prio can only be done if all sons are  */
  /*   available in SetGhostObjectPriorities()                       */
  /* - setting of the border priorities can only be done if all      */
  /*   ghost objects have their proper priority                      */

  DDD_XferBegin();
  SetGhostObjectPriorities(theGrid);
  DDD_XferEnd();

  DDD_XferBegin();
  SetBorderPriorities(theGrid);
  DDD_XferEnd();

  ConstructConsistentGridLevel(theGrid);
}

/****************************************************************************/
/*
   ConstructConsistentMultiGrid -

   SYNOPSIS:
   void ConstructConsistentMultiGrid (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG

   DESCRIPTION:
   Provide a consistent multigrid.

   RETURN VALUE:
   void
 */
/****************************************************************************/

void NS_DIM_PREFIX ConstructConsistentMultiGrid (MULTIGRID *theMG)
{
  INT l;

  /* this is done in three waves:   */
  /* 1. set priorities of objects   */
  /* 2. set border priorities       */
  /* 3. repair grid inconsistencies */

  /* 1. set priorities of objects   */
  DDD_XferBegin();
  for (l=0; l<=TOPLEVEL(theMG); l++)
  {
    GRID *theGrid = GRID_ON_LEVEL(theMG,l);
    SetGhostObjectPriorities(theGrid);
  }
  DDD_XferEnd();

  /* 2. set border priorities       */
  DDD_XferBegin();
  for (l=0; l<=TOPLEVEL(theMG); l++)
  {
    GRID *theGrid = GRID_ON_LEVEL(theMG,l);
    SetBorderPriorities(theGrid);
  }
  DDD_XferEnd();

  /* 3. repair grid inconsistencies */
  for (l=0; l<=TOPLEVEL(theMG); l++)
  {
    GRID *theGrid = GRID_ON_LEVEL(theMG,l);
    ConstructConsistentGridLevel(theGrid);
  }

}

#endif
