// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  std_parallel.c				                                                                */
/*																			*/
/* Purpose:   parallel part of standard ug domain description                           */
/*																			*/
/* Author:	  Klaus Birken / Christian Wieners								*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   Sep 12 1996 ug version 3.4                                                                */
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

/* standard C library */
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* low modules */
#include "compiler.h"
#include "heaps.h"
#include "general.h"
#include "debug.h"

/* parallel modules */
#include "parallel.h"

/* domain module */
#include "domain.h"

#include "std_domain.h"

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
/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

void DomInitParallel (INT TypeBndP, INT TypeBndS)
{}

void DomHandlerInit (INT handlerSet)
{}

void BElementXferBndS (BNDS **bnds, int n, int proc, int prio)
{
  INT size,i,size0;

  size = CEIL(sizeof(INT));
  for (i=0; i<n; i++)
    if (bnds[i] != NULL)
    {
      size0 = BND_SIZE(bnds[i]);
      size += CEIL(size0) + CEIL(sizeof(INT));

      PRINTDEBUG(dom,1,("Xfer  me %x %d pid %d n %d size %d\n",
                        me,bnds[i],BND_PATCH_ID(bnds[i]),
                        BND_N(bnds[i]),BND_SIZE(bnds[i])));

    }

  DDD_XferAddData(size,DDD_USER_DATA);
}

void BElementGatherBndS (BNDS **bnds, int n, int cnt, char *data)
{
  INT size,i;

  for (i=0; i<n; i++)
    if (bnds[i] != NULL)
    {

      PRINTDEBUG(dom,1,("Gather %d  me %d %x pid %d n %d size %d\n",i,
                        me,bnds[i],BND_PATCH_ID(bnds[i]),
                        BND_N(bnds[i]),BND_SIZE(bnds[i])));


      size = BND_SIZE(bnds[i]);
      memcpy(data,&i,sizeof(INT));
      data += CEIL(sizeof(INT));
      memcpy(data,bnds[i],size);
      data += CEIL(size);
    }
  i = -1;
  memcpy(data,&i,sizeof(INT));
}

void BElementScatterBndS (BNDS **bnds, int n, int cnt, char *data)
{
  INT size,i;
  BNDS *bs;

  memcpy(&i,data,sizeof(INT));
  while (i != -1)
  {
    data += CEIL(sizeof(INT));
    bs = (BNDS *) data;
    size = BND_SIZE(bs);

    PRINTDEBUG(dom,1,("Scatter %d me %d\n",i,size));

    if (bnds[i] == NULL)
    {
      bs = (BNDS *) memmgr_AllocOMEM((size_t)size,TypeBndS,0,0);
      memcpy(bs,data,size);
      bnds[i] = bs;
    }
    data += CEIL(size);
    memcpy(&i,data,sizeof(INT));
  }
}

void BVertexXferBndP (BNDP *bndp, int proc, int prio)
{
  INT size;

  size = BND_SIZE(bndp);

  PRINTDEBUG(dom,1,("  me %x %d pid %d n %d size %d\n",
                    me,bndp,BND_PATCH_ID(bndp),BND_N(bndp),BND_SIZE(bndp)));

  DDD_XferAddData(size,DDD_USER_DATA);
}

void BVertexGatherBndP (BNDP *bndp, int cnt, char *data)
{
  PRINTDEBUG(dom,1,("  me %d pid %d n %d size %d cnt %d\n",
                    me,BND_PATCH_ID(bndp),
                    BND_N(bndp),BND_SIZE(bndp),cnt));

  ASSERT(cnt == BND_SIZE(bndp));

  memcpy(data,bndp,cnt);
}


void BVertexScatterBndP (BNDP **bndp, int cnt, char *data)
{
  if (*bndp == NULL)
  {
    *bndp = (BNDS *) memmgr_AllocOMEM((size_t)cnt,TypeBndP,0,0);
    memcpy(*bndp,data,cnt);
    PRINTDEBUG(dom,1,("  me %d pid %d n %d size %d cnt %d\n",
                      me,BND_PATCH_ID(*bndp),
                      BND_N(*bndp),BND_SIZE(*bndp),cnt));
  }
}


#endif /* ModelP */
