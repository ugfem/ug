// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  lb.c													        */
/*																			*/
/* Purpose:   a collect of simple and special load balancing functionality  */
/*																			*/
/* Author:	  Stefan Lang, Klaus Birken, Christian Wieners	                        */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*																			*/
/* History:   960906 kb  begin                                                                                          */
/*            980204 sl  renamed test.c to lb.c                             */
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

#include <stdio.h>

#include "parallel.h"
#include "general.h"
#include "ugm.h"
#include "ugdevices.h"

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
    TransferGridComplete-

   SYNOPSIS:
   static int TransferGridComplete (MULTIGRID *theMG, INT level);

   PARAMETERS:
   .  theMG
   .  level

   DESCRIPTION:

   RETURN VALUE:
   int
 */
/****************************************************************************/

static int TransferGridComplete (MULTIGRID *theMG, INT level)
{
  ELEMENT *e;
  GRID *theGrid = GRID_ON_LEVEL(theMG,level);

  if (theGrid==NULL)
  {
    UserWriteF(PFMT "TransferGridComplete(): no grid on level=%d\n",me,level);
    return(0);
  }

  /* assign elements of level 0 */
  if (me == master) {
    for (e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
      PARTITION(e) = 1;
  }

  IFDEBUG(dddif,1);
  for (e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
  {
    UserWriteF("elem %08x has dest=%d\n",
               DDD_InfoGlobalId(PARHDRE(e)), PARTITION(e));
  }
  ENDDEBUG

  return(0);
}


/****************************************************************************/
/*
   TransferGridToMaster -

   SYNOPSIS:
   static int TransferGridToMaster (MULTIGRID *theMG, INT fl, INT tl);

   PARAMETERS:
   .  theMG
   .  fl
   .  tl

   DESCRIPTION:

   RETURN VALUE:
   int
 */
/****************************************************************************/

static int TransferGridToMaster (MULTIGRID *theMG, INT fl, INT tl)
{
  ELEMENT *e;
  GRID *theGrid;

  /* send all levels to master */
  if (me!=master)
  {
    int l;

    for (l=fl; l<=tl; l++) {

      theGrid = GRID_ON_LEVEL(theMG,l);

      /* create element copies */
      for(e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
      {
        PARTITION(e) = 0;
      }
    }
  }

  return(0);
}


/****************************************************************************/
/*
   CollectElementsNearSegment -

   SYNOPSIS:
   static int CollectElementsNearSegment(MULTIGRID *theMG, int level, int part, int dest);

   PARAMETERS:
   .  theMG
   .  level
   .  part
   .  dest

   DESCRIPTION:

   RETURN VALUE:
   int
 */
/****************************************************************************/

static int CollectElementsNearSegment(MULTIGRID *theMG,
                                      int fl, int tl, int part, int dest)
{
  ELEMENT *theElement;
  INT dompart,side,sid,nbsid,level;

  for (level=fl; level<=tl; level ++)
    for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,level));
         theElement!=NULL; theElement=SUCCE(theElement))
      if (OBJT(theElement) == BEOBJ)
        for (side=0; side<SIDES_OF_ELEM(theElement); side++) {
          if (INNER_SIDE(theElement,side))
            continue;
          BNDS_BndSDesc(ELEM_BNDS(theElement,side),
                        &sid,&nbsid,&dompart);
          if (part == dompart)
            PARTITION(theElement) = dest;
        }

  return(0);
}

/****************************************************************************/
/*
   CreateDD -

   SYNOPSIS:
   static int CreateDD(GRID *theGrid, int hor_boxes, int vert_boxes );

   PARAMETERS:
   .  theGrid
   .  hor_boxes
   .  vert_boxes

   DESCRIPTION:

   RETURN VALUE:
   int
 */
/****************************************************************************/

static int PartitionElementsForDD(GRID *theGrid, int hor_boxes, int vert_boxes )
{
  ELEMENT *theElement;
  INT i, nrcorners;
  DOUBLE *coord, xmax, ymax;

  for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    nrcorners = CORNERS_OF_ELEM(theElement);
    ASSERT(nrcorners==4 || nrcorners==3 );              /* works only for triangle and quadrilateral grids */

    /* calculate the coordinates xmax, ymax of the element */
    xmax = ymax = 0.0;
    for( i=0; i<nrcorners; i++ )
    {
      coord = CVECT(MYVERTEX(CORNER(theElement,i)));
      xmax = MAX(xmax,coord[0]);
      ymax = MAX(ymax,coord[1]);
    }
    ASSERT(xmax>0.0 && xmax<1.00001);                   /* this function works only for the unit square! */
    ASSERT(ymax>0.0 && ymax<1.00001);
    /*printf( PFMT "element coord %g %g %d %d\n", me, xmax, ymax, hor_boxes, vert_boxes );*/

    /* the according subdomain is determined by the upper right corner
       the idea: for each dimension: calculate from the element coord the position in the
                                                     corresponding array dimension
                             then map the array coordinates onto a PE number */
    PARTITION(theElement) = (int)((ymax-0.00001)*vert_boxes) * hor_boxes + (int)((xmax-0.00001)*hor_boxes);
    /*printf( PFMT "element partition %d\n", me, PARTITION(theElement) );*/
  }

  return(0);
}

static int CreateDD(MULTIGRID *theMG, INT level, int hor_boxes, int vert_boxes )
/* the final call to TransferGridFromLevel() must be done by the current caller */
{
  INT elements;
  int hor, vert;
  GRID *theGrid;

  theGrid = GRID_ON_LEVEL(theMG,level);
  if( hor_boxes*vert_boxes >= 4 )
  {
    elements = NT(theGrid);
    elements = UG_GlobalMaxINT(elements);

    if( elements > 20000 )
    {
      /* we have a case with too heavy load for DDD; thus subdivide */
      /* find a coarser auxiliary partition */
      hor = hor_boxes;
      vert = vert_boxes;

      if( hor%2 == 0 )
        hor /= 2;                               /* hor_boxes can be halfened */
      else if( vert%2 == 0 )
        vert /= 2;                              /* vert_boxes can be halfened */
      else
        assert(0);                              /* further consideration for subdividing necessary */

      CreateDD(theMG, level, hor, vert );                       /* simplify the problem */
      TransferGridFromLevel(theMG,level);
    }
  }

  PartitionElementsForDD(theGrid, hor_boxes, vert_boxes );
}

/****************************************************************************/
/*
   lbs -  interface for simple or special load balancing functionality

   SYNOPSIS:
   void lbs (char *argv, MULTIGRID *theMG);

   PARAMETERS:
   .  argv
   .  theMG

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void lbs (char *argv, MULTIGRID *theMG)
{
  int n,mode,param,fromlevel,tolevel,part,hor_boxes,vert_boxes,dest;

  mode = param = fromlevel = tolevel = 0;

  n = sscanf(argv,"%d %d %d",&param,&fromlevel,&tolevel);
  UserWriteF(PFMT "lbs() param=%d",me,param);
  if (n > 1)
    UserWriteF(" fromlevel=%d",fromlevel);
  if (n > 2)
    UserWriteF(" tolevel=%d",tolevel);
  UserWriteF("\n");

  /* param>100 is used as switch for DDD xfer statistics */
  if (param>=100)
    mode = param-100;
  else
    mode = param;

  /* switch DDD infos on */
  if (param>=100)
    DDD_SetOption(OPT_INFO_XFER, XFER_SHOW_MEMUSAGE);

  /* TODO: remove the InitCurrMG call completely; InitCurrMG may be called
           only once and this is already done in new/open command */
  /*InitCurrMG(theMG);*/
  switch (mode)
  {
  /* transfer vectors of coarsest amg grid to master */
  case (-1) :
    AMGAgglomerate(theMG);
    return;

  /* dies balanciert ein GRID mit RCB */
  case (0) :
    BalanceGridRCB(theMG,0);
    fromlevel = 0;
    break;

  /* dies verschickt ein GRID komplett */
  case (1) :
    TransferGridComplete(theMG,fromlevel);
    break;

  /* dies verschickt ein verteiltes GRID zum master */
  case (2) :
    TransferGridToMaster(theMG,fromlevel,tolevel);
    fromlevel = 0;
    break;

  /* dies balanciert ein GRID mit RCB ab fromlevel */
  case (3) :
    if (fromlevel>=0 && fromlevel<=TOPLEVEL(theMG))
    {
      BalanceGridRCB(theMG,fromlevel);
    }
    else
    {
      UserWriteF(PFMT "lbs(): gridlevel=%d not "
                 "existent!\n",me,fromlevel);
    }
    break;

  /* dies balanciert ein GRID mit RCB ab fromlevel */
  case (4) :
    if (fromlevel>=0 && fromlevel<=TOPLEVEL(theMG) ||
        tolevel>=0 && tolevel<=TOPLEVEL(theMG)     ||
        tolevel < fromlevel)
    {
      int j;

      for (j=fromlevel; j<=tolevel; j++)
        BalanceGridRCB(theMG,j);
      /*
                                      TransferGrid(theMG,fromlevel,tolevel);
       */
    }
    else
    {
      UserWriteF(PFMT "lbs(): ERROR fromlevel=%d "
                 "tolevel=%d\n",me,fromlevel,tolevel);
    }
    break;

  case (5) :
    n = sscanf(argv,"%d %d %d %d %d",
               &param,&part,&dest,&fromlevel,&tolevel);
    if (n < 5) tolevel = TOPLEVEL(theMG);
    if (n < 4) fromlevel = CURRENTLEVEL(theMG);
    if (n < 3) break;
    CollectElementsNearSegment(theMG,fromlevel,tolevel,part,dest);
    UserWriteF(PFMT "lbs() collect from part %d to proc %d\n",
               me,part,dest);
    break;

  /* dies erzeugt eine regelmaessige Domain Decomposition */
  case (6) :
    if (sscanf(argv,"%d %d %d",&param,&hor_boxes,&vert_boxes) != 3) break;
    ASSERT(hor_boxes*vert_boxes == procs );
    fromlevel = TOPLEVEL(theMG);
    CreateDD(theMG,fromlevel,hor_boxes,vert_boxes);
    break;

  default : break;
  }

  TransferGridFromLevel(theMG,fromlevel);

  /* switch DDD infos off */
  if (param>=100)
    DDD_SetOption(OPT_INFO_XFER, XFER_SHOW_NONE);
}

#endif /* ModelP */
