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


static int TransferGridToMaster (MULTIGRID *theMG)
{
  ELEMENT *e;
  GRID *theGrid;

  /* send all levels to master */
  if (me!=master)
  {
    int l;

    for (l=0; l<=TOPLEVEL(theMG); l++) {

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

static int CollectElementsNearSegment(MULTIGRID *theMG,
                                      int level, int part, int p)
{
  GRID *theGrid = GRID_ON_LEVEL(theMG,level);
  ELEMENT *theElement;
  INT dompart,side,sid,nbsid;

  for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL;
       theElement=SUCCE(theElement))
    if (OBJT(theElement) == BEOBJ)
      for (side=0; side<SIDES_OF_ELEM(theElement); side++) {
        if (INNER_SIDE(theElement,side))
          continue;
        BNDS_BndSDesc(ELEM_BNDS(theElement,side),&sid,&nbsid,&dompart);
        if (part == dompart)
          PARTITION(theElement) = 0;
      }

  return(0);
}

void ddd_test (char *argv, MULTIGRID *theMG)
{
  int mode,param,fromlevel,tolevel,part;

  mode = param = fromlevel = tolevel = 0;

  sscanf(argv,"%d %d %d",&param,&fromlevel,&tolevel);
  UserWriteF(PFMT "ddd_test() param=%d fromlevel=%d tolevel=%d\n",
             me,param,fromlevel,tolevel);

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
    TransferGridToMaster(theMG);
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

  /* dies balanciert ein GRID mit RCB */
  case (5) :
    if (sscanf(argv,"%d %d",&param,&part) != 2) break;
    fromlevel = CURRENTLEVEL(theMG);
    CollectElementsNearSegment(theMG,fromlevel,part,0);
    break;

  default : break;
  }

  TransferGridFromLevel(theMG,fromlevel);

  /* switch DDD infos off */
  if (param>=100)
    DDD_SetOption(OPT_INFO_XFER, XFER_SHOW_NONE);
}

#endif /* ModelP */
