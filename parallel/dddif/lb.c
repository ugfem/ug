// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef ModelP


#include <stdio.h>


#include "parallel.h"
#include "general.h"
#include "ugm.h"


/* RCS string */
RCSID("$Header$",UG_RCS_STRING)


static int TransferGridComplete (MULTIGRID *theMG)
{
  ELEMENT *e;
  GRID *theGrid = GRID_ON_LEVEL(theMG,0);

  /* assign elements of level 0 */
  if (me == master) {
    for (e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
      PARTITION(e) = 1;
  }

  for (e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
  {
    UserWriteF("elem %08x has dest=%d\n",
               DDD_InfoGlobalId(PARHDRE(e)), PARTITION(e));
  }

  /* start physical transfer */
  ddd_HandlerInit(HSET_XFER);
  DDD_XferBegin();

  if (me==master) {

    TransferGridFromCoarse(theMG);
    if (0) {
      /* create element copies */
      for(e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
      {
        /* create element copy */
        DDD_XferCopyObjX(PARHDRE(e),
                         PARTITION(e),
                         PrioMaster,
                         (OBJT(e)==BEOBJ) ? BND_SIZE_TAG(TAG(e)) : INNER_SIZE_TAG(TAG(e))
                         );

        /* delete local copy */
        DDD_XferDeleteObj(PARHDRE(e));
      }
    }
  }

  DDD_XferEnd();

  DDD_ConsCheck();

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

  TransferGridFromCoarse(theMG);

  return(0);
}

void ddd_test (int param, MULTIGRID *theMG)
{
  int mode;

  /* param>100 is used as switch for DDD xfer statistics */
  if (param>=100)
    mode = param-100;
  else
    mode = param;

  /* switch DDD infos on */
  if (param>=100)
    DDD_SetOption(OPT_INFO_XFER, XFER_SHOW_MEMUSAGE);

  InitCurrMG(theMG);
  switch (mode) {
  /* dies balanciert ein GRID mit RCB */
  case (0) :
    BalanceGridRCB(theMG);
    TransferGridFromCoarse(theMG);
    break;

  /* dies verschickt ein GRID komplett */
  case (1) :
    TransferGridComplete(theMG);
    break;

  /* dies verschickt ein verteiltes GRID zum master */
  case (2) :
    TransferGridToMaster(theMG);
    break;

  default : break;
  }

  /* switch DDD infos off */
  if (param>=100)
    DDD_SetOption(OPT_INFO_XFER, XFER_SHOW_NONE);
}

#endif /* ModelP */
