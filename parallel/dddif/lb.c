// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef ModelP


#include <stdio.h>


#include "parallel.h"


int TransferGrid (MULTIGRID *theMG)
{
  ELEMENT *e;
  GRID *theGrid = GRID_ON_LEVEL(theMG,CURRENTLEVEL(theMG));

  /* assign elements */
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
  DDD_XferBegin();

  if (me==master) {

    /* create element copies */
    for(e=FIRSTELEMENT(theGrid); e!=NULL; e=SUCCE(e))
    {
      /* create element copy */
      DDD_XferCopyObjX(PARHDRE(e),
                       PARTITION(e),
                       0,
                       (OBJT(e)==BEOBJ) ? BND_SIZE(TAG(e)) : INNER_SIZE(TAG(e))
                       );
    }
  }

  DDD_XferEnd();

  DDD_ConsCheck();

  return(0);
}


void ddd_test (int mode, MULTIGRID *theMG)
{
  InitCurrMG(theMG);
  switch (mode) {
  /* dies balanciert ein GRID mit RCB */
  case (0) : BalanceGrid(theMG); break;

  /* dies verschickt ein GRID komplett */
  case (1) : TransferGrid(theMG); break;

  default : break;
  }
}

#endif /* ModelP */
