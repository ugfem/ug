// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifdef ModelP


#include <stdio.h>


#include "parallel.h"


void ddd_DisplayContext (void)
{
  int i, last=-1;
  char sep[2] = "";
  char buf[20];

  /* only master should display context */
  if (me!=master)
    return;

  UserWrite("current context: (");
  for(i=0; i<procs+1; i++)
  {
    if (i>=procs || !CONTEXT(i))
    {
      if (last+1==i-1)
      {
        sprintf(buf, "%s%d", sep, last+1); UserWrite(buf);
        sep[0] = ',';
      }
      if (last+1<i-1)
      {
        sprintf(buf, "%s%d-%d", sep, last+1, i-1); UserWrite(buf);
        sep[0] = ',';
      }
      last = i;
    }
  }
  UserWrite(")\n");
}



void ddd_pstat (int cmd)
{
  switch (cmd)
  {
  case 'c' :
    DDD_ConsCheck();
    UserWrite("\n");
    break;

  case 's' :
    SYNC_CONTEXT;
    DDD_Status();
    UserWrite("\n");
    SYNC_END;
    break;

  case 'i' :
    SYNC_CONTEXT;
    DDD_IFDisplay();
    UserWrite("\n");
    SYNC_END;
    break;

  case 'l' :
    SYNC_CONTEXT;
    DDD_ListLocalObjects();
    UserWrite("\n");
    SYNC_END;
  }
}

#endif /* ModelP */
