// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  lgm_parallel.c				                                                                */
/*																			*/
/* Purpose:   parallel part of lgm ug domain description                                        */
/*																			*/
/* Author:	  Stefan Lang / Klaus Birken / Christian Wieners				*/
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

#include "lgm_domain.h"

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

extern LGM_LINE **LinePtrArray;

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

void BElementXferBndS (BNDS **Bnds, int n, int proc, int prio)
{
  LGM_BNDS **bnds = (LGM_BNDS **)Bnds;
  INT size,i,size0;

  PRINTDEBUG(dom,1,(PFMT "BElementXferBndS(): bnds=%x n=%d proc=%d prio=%x\n",
                    me,bnds,n,proc,prio));

  size = sizeof(INT);
  for (i=0; i<n; i++)
    if (bnds[i] != NULL)
    {
      size0 = sizeof(*bnds[i]);
      size += CEIL(size0) + CEIL(sizeof(INT));

      PRINTDEBUG(dom,2,(PFMT "BElementXferBndS(): Xfer i=%d bnds=%x "
                        "lineid=%d size=%d\n",
                        me,i,bnds[i],LGM_LINE_ID(LGM_BNDS_LINE(bnds[i])),size));
    }

  DDD_XferAddData(size,DDD_DOMAIN_DATA);
}

void BElementGatherBndS (BNDS **Bnds, int n, int cnt, char *data)
{
  LGM_BNDS **bnds = (LGM_BNDS **)Bnds;
  INT size,i;

  PRINTDEBUG(dom,1,(PFMT "BElementGatherBndS(): bnds=%x n=%d cnt=%d "
                    "data=%x\n",me,bnds,n,cnt,data));

  for (i=0; i<n; i++)
    if (bnds[i] != NULL)
    {

      PRINTDEBUG(dom,2,(PFMT "BElementGatherBndS(): i=%d bnds=%x "
                        "lineid=%d size=%d\n",
                        me,i,bnds[i],LGM_LINE_ID(LGM_BNDS_LINE(bnds[i])),
                        LGM_BNDS_SIZE(bnds[i])));

      size = sizeof(*bnds[i]);

      memcpy(data,&i,sizeof(INT));
      data += CEIL(sizeof(INT));

      memcpy(data,bnds[i],size);

      /* substiute line pointer by line id */
      PRINTDEBUG(dom,3,(PFMT "BElementGatherBndS(): substituting "
                        "lineptr=%x by lineid=%d\n",me,LGM_BNDS_LINE((LGM_BNDS *)data),
                        LGM_LINE_ID(LGM_BNDS_LINE(bnds[i]))));

      LGM_BNDS_LINE((LGM_BNDS *)data) =
        (LGM_LINE *) LGM_LINE_ID(LGM_BNDS_LINE(bnds[i]));
      data += CEIL(size);
    }

  i = -1;
  memcpy(data,&i,sizeof(INT));
}

void BElementScatterBndS (BNDS **Bnds, int n, int cnt, char *data)
{
  LGM_BNDS **bnds = (LGM_BNDS **)Bnds;
  INT size,i;
  LGM_BNDS *bs;

  PRINTDEBUG(dom,1,(PFMT "BElementScatterBndS(): bnds=%x n=%d cnt=%d data=%x\n",
                    me,bnds,n,cnt,data));

  memcpy(&i,data,sizeof(INT));
  while (i != -1)
  {
    data += CEIL(sizeof(INT));
    size = sizeof(*bs);

    if (bnds[i] == NULL)
    {
      PRINTDEBUG(dom,2,(PFMT "BElementScatterBndS(): scatter i=%d "
                        "size=%d\n",me,i,size));

      bs = (LGM_BNDS *) memmgr_AllocOMEM((size_t)size,TypeBndS,0,0);
      ASSERT(bs!=NULL);

      memcpy(bs,data,size);
      bnds[i] = bs;

      /* substitute line id by its local pointer */
      PRINTDEBUG(dom,3,(PFMT "BElementScatterBndS(): substituting "
                        "lineid=%d by lineptr=%x\n",me,
                        (INT) LGM_BNDS_LINE((LGM_BNDS *)data),
                        LINE_ID_2_LINE((INT) LGM_BNDS_LINE((LGM_BNDS *)data))));

      LGM_BNDS_LINE(bs) =
        LINE_ID_2_LINE((INT) LGM_BNDS_LINE((LGM_BNDS *)data));
    }
    else
      PRINTDEBUG(dom,2,(PFMT "BElementScatterBndS(): ignoring i=%d "
                        "size=%d\n",me,i,size));

    data += CEIL(size);
    memcpy(&i,data,sizeof(INT));
  }
}

void BVertexXferBndP (BNDP *Bndp, int proc, int prio)
{
  LGM_BNDP *bndp = (LGM_BNDP *)Bndp;
  INT i,size;
  INT n;

  PRINTDEBUG(dom,1,(PFMT "BVertexXferBndP(): bndp=%x proc=%d "
                    "prio=%d\n",me,bndp,proc,prio));

  n               = LGM_BNDP_N(bndp);
  size    = CEIL(sizeof(INT))+CEIL(n*sizeof(LGM_BNDP_PLINE));

  PRINTDEBUG(dom,1,(PFMT "BVertexXferBndP(): bndp=%x n=%d size=%d\n"
                    ,me,bndp,n,size));

  DDD_XferAddData(size,DDD_DOMAIN_DATA);
}

void BVertexGatherBndP (BNDP *Bndp, int cnt, char *data)
{
  LGM_BNDP *bndp = (LGM_BNDP *)Bndp;
  INT i;
  INT n;
  INT size;
  char *buffer = data;

  PRINTDEBUG(dom,1,(PFMT "BVertexGatherBndP(): bndp=%x cnt=%d data=%x\n",
                    me,bndp,cnt,data));

  n               = LGM_BNDP_N(bndp);
  size    = CEIL(n*sizeof(LGM_BNDP_PLINE));

  PRINTDEBUG(dom,1,(PFMT "BVertexGatherBndP(): bndp=%x n=%d size=%d\n"
                    ,me,bndp,n,size));

  memcpy(data,&n,sizeof(INT));
  data += CEIL(sizeof(INT));
  memcpy(data,&LGM_BNDP_LINES(bndp,0),size);

  /* substitute the line pointer by their ids */
  for (i=0; i<n; i++)
  {
    PRINTDEBUG(dom,3,(PFMT "BVertexGatherBndP(): substituting i=%d "
                      "lineptr=%x by lineid=%d local=%11.4E\n",
                      me,i,LGM_BNDP_LINE(bndp,i),
                      LGM_LINE_ID(LGM_BNDP_LINE(bndp,i)),
                      LGM_BNDP_LOCAL(bndp,i)));
    fflush(stdout);

    LGM_BNDP_LINE_GLINE(((LGM_BNDP_PLINE*)data)[0]) =
      (LGM_LINE *)LGM_LINE_ID(LGM_BNDP_LINE(bndp,i));

    PRINTDEBUG(dom,3,(PFMT "BVertexGatherBndP(): inside buffer: "
                      "line=%x local=%11.4E\n",me,
                      LGM_BNDP_LINE_GLINE(((LGM_BNDP_PLINE*)data)[0]),
                      LGM_BNDP_LINE_LOCAL(((LGM_BNDP_PLINE*)data)[0])));

    data += sizeof(LGM_BNDP_PLINE);
  }
}

void BVertexScatterBndP (BNDP **Bndp, int cnt, char *data)
{
  LGM_BNDP **bndp = (LGM_BNDP **)Bndp;
  INT i;
  char *buffer = data;

  PRINTDEBUG(dom,1,(PFMT "BVertexScatterBndP(): bndp=%x cnt=%d data=%x\n",
                    me,*bndp,cnt,data));

  if (*bndp == NULL)
  {
    INT n           = LGM_BNDP_N((LGM_BNDP *)data);
    INT size        = CEIL(n*sizeof(LGM_BNDP_PLINE));

    *bndp = (LGM_BNDP *) memmgr_AllocOMEM(sizeof(LGM_BNDP)+(n-1)*
                                          sizeof(struct lgm_bndp_line),TypeBndP,0,0);
    ASSERT(*bndp!=NULL);

    PRINTDEBUG(dom,2,(PFMT "BVertexScatterBndP(): scatter bndp=%x n=%d "
                      "size %d\n",me,*bndp,n,size))

    memcpy(&LGM_BNDP_N(*bndp),&n,sizeof(INT));
    data += CEIL(sizeof(INT));
    memcpy(&LGM_BNDP_LINES(*bndp,0),data,size);

    /* substitute the line ids by their pointers */
    for (i=0; i<n; i++)
    {
      /* substitute line id by its local pointer */
      PRINTDEBUG(dom,3,(PFMT "BVertexScatterBndP(): substituting i=%d "
                        "lineid=%d by lineptr=%x local=%11.4E\n",
                        me,i,(INT) LGM_BNDP_LINE_GLINE(((LGM_BNDP_PLINE *)data)[0]),
                        LINE_ID_2_LINE((INT)
                                       LGM_BNDP_LINE_GLINE(((LGM_BNDP_PLINE *)data)[0])),
                        LGM_BNDP_LINE_LOCAL(((LGM_BNDP_PLINE *)data)[0])));

      LGM_BNDP_LINE(*bndp,i) =
        LINE_ID_2_LINE((INT)
                       LGM_BNDP_LINE_GLINE(((LGM_BNDP_PLINE *)data)[0]));

      PRINTDEBUG(dom,3,(PFMT "BVertexScatterBndP(): outside buffer: "
                        "line=%x local=%11.4E\n",me,
                        LGM_BNDP_LINE(*((LGM_BNDP **)Bndp),i),
                        LGM_BNDP_LOCAL(*((LGM_BNDP **)Bndp),i)));

      data += sizeof(LGM_BNDP_PLINE);
    }
  }
  else
  {
    PRINTDEBUG(dom,2,(PFMT "BVertexScatterBndP(): ignoring bndp=%x",
                      me,*bndp))
  }
}

#endif /* ModelP */
