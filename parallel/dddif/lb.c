// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef ModelP


#include <stdio.h>


#include "parallel.h"
#include "general.h"
#include "ugm.h"
#include "devices.h"

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


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

static int CreateDD(GRID *theGrid, int hor_boxes, int vert_boxes )
{
  ELEMENT *theElement;
  INT i;
  DOUBLE *coord, xmax, ymax;

  for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL;
       theElement=SUCCE(theElement))
  {
    ASSERT(CORNERS_OF_ELEM(theElement)==4);             /* works only for quadrilateral grids */

    /* calculate the coordinates xmax, ymax of the element */
    xmax = ymax = 0.0;
    for( i=0; i<4; i++ )
    {
      coord = CVECT(MYVERTEX(CORNER(theElement,i)));
      xmax = MAX(xmax,coord[0]);
      ymax = MAX(ymax,coord[1]);
    }
    printf( PFMT "element coord %g %g %d %d\n", me, xmax, ymax, hor_boxes, vert_boxes );

    /* the according subdomain is determined by the upper right corner */
    PARTITION(theElement) = (int)(ymax*vert_boxes - 0.5) * hor_boxes + (int)(xmax*hor_boxes - 0.5);
    printf( PFMT "element partition %d\n", me, PARTITION(theElement) );
  }

  return(0);
}


/****************************************************************************/
/*
   ddd_test -

   SYNOPSIS:
   void ddd_test (char *argv, MULTIGRID *theMG);

   PARAMETERS:
   .  argv
   .  theMG

   DESCRIPTION:

   RETURN VALUE:
   void
 */
/****************************************************************************/

void ddd_test (char *argv, MULTIGRID *theMG)
{
  int n,mode,param,fromlevel,tolevel,part,hor_boxes,vert_boxes,dest;

  mode = param = fromlevel = tolevel = 0;

  n = sscanf(argv,"%d %d %d",&param,&fromlevel,&tolevel);
  UserWriteF(PFMT "ddd_test() param=%d",me,param);
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

  InitCurrMG(theMG);
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
      UserWriteF(PFMT "ddd_test(): gridlevel=%d not "
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
      UserWriteF(PFMT "ddd_test(): ERROR fromlevel=%d "
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
    UserWriteF(PFMT "ddd_test() collect from part %d to proc %d\n",
               me,part,dest);
    break;

  /* dies erzeugt eine regelmaessige Domain Decomposition */
  case (6) :
    if (sscanf(argv,"%d %d %d",&param,&hor_boxes,&vert_boxes) != 3) break;
    ASSERT(hor_boxes*vert_boxes == procs );
    CreateDD(GRID_ON_LEVEL(theMG,TOPLEVEL(theMG)),hor_boxes,vert_boxes);
    break;

  default : break;
  }

  TransferGridFromLevel(theMG,fromlevel);

  /* switch DDD infos off */
  if (param>=100)
    DDD_SetOption(OPT_INFO_XFER, XFER_SHOW_NONE);
}

#endif /* ModelP */
