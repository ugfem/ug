// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef ModelP


#include <stdio.h>


#include "gm.h"
#include "parallel.h"



int theDebugProc = 0;



void ddd_pstat (int cmd)
{
  switch (cmd)
  {
  case 'c' :
    DDD_ConsCheck();
    break;

  case 's' :
    if (me!=theDebugProc) return;
    DDD_Status();
    break;

  case 'i' :
    if (me!=theDebugProc) return;
    DDD_DisplayIF();
    break;

  case 'l' :
    if (me!=theDebugProc) return;
    DDD_ListLocalObjects();
  }

  if (me==theDebugProc) UserWrite("\n");
}

#endif /* ModelP */
