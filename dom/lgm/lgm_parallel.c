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

#if (LGM_DIM==2)
#define LGM_OBJECT                              LGM_LINE
#define LGM_OBJECT_ID                   LGM_LINE_ID
#define LGM_OBJECT_ID_2_OBJECT  LGM_LINE_ID_2_LINE
#define LGM_BNDS_OBJECT                 LGM_BNDS_LINE
#define LGM_BNDP_OBJECT                 LGM_BNDP_LINE
#define LGM_BNDP_POBJECT                LGM_BNDP_PLINE
#define LGM_OBJECT_ID_2_OBJECT  LGM_LINE_ID_2_LINE
#define LGM_BNDP_OBJECTS                LGM_BNDP_LINES
#define LGM_BNDP_OBJECT_GOBJECT LGM_BNDP_LINE_GLINE
#define LGM_BNDP_OBJECT_LOCAL   LGM_BNDP_LINE_LOCAL
#endif

#if (LGM_DIM==3)
#define LGM_OBJECT                              LGM_SURFACE
#define LGM_OBJECT_ID                   LGM_SURFACE_ID
#define LGM_OBJECT_ID_2_OBJECT  LGM_SURFACE_ID_2_SURFACE
#define LGM_BNDS_OBJECT                 LGM_BNDS_SURFACE
#define LGM_BNDP_OBJECT                 LGM_BNDP_SURFACE
#define LGM_BNDP_POBJECT                LGM_BNDP_PSURFACE
#define LGM_BNDP_OBJECTS                LGM_BNDP_SURFACES
#define LGM_BNDP_OBJECT_GOBJECT LGM_BNDP_SURFACE_GSURFACE
#define LGM_BNDP_OBJECT_LOCAL   LGM_BNDP_SURFACE_LOCAL
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

#if (LGM_DIM==2)
extern LGM_LINE                        **LinePtrArray;
#endif

#if (LGM_DIM==3)
extern LGM_SURFACE                     **SurfacePtrArray;
#endif

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


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

  /* to store end tag */
  size = sizeof(INT);

  for (i=0; i<n; i++)
    if (bnds[i] != NULL)
    {
      size0 = sizeof(*bnds[i]);
      size += CEIL(size0) + CEIL(sizeof(INT));

      PRINTDEBUG(dom,2,(PFMT "BElementXferBndS(): Xfer i=%d bnds=%x "
                        "line/surfaceid=%d size=%d\n",
                        me,i,bnds[i],LGM_OBJECT_ID(LGM_BNDS_OBJECT(bnds[i])),size));
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
                        "line/surfaceid=%d size=%d\n",
                        me,i,bnds[i],LGM_OBJECT_ID(LGM_BNDS_OBJECT(bnds[i])),
                        sizeof(LGM_BNDS)));
      size = sizeof(*bnds[i]);

      memcpy(data,&i,sizeof(INT));
      data += CEIL(sizeof(INT));

      memcpy(data,bnds[i],size);

      /* substiute line pointer by line id */
      PRINTDEBUG(dom,3,(PFMT "BElementGatherBndS(): substituting "
                        "line/surfaceptr=%x by line/surfaceid=%d\n",me,LGM_BNDS_OBJECT((LGM_BNDS *)data),
                        LGM_OBJECT_ID(LGM_BNDS_OBJECT(bnds[i]))));

      LGM_BNDS_OBJECT((LGM_BNDS *)data) =
        (LGM_OBJECT *) LGM_OBJECT_ID(LGM_BNDS_OBJECT(bnds[i]));

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
                        "line/surfaceid=%d by line/surfaceptr=%x\n",me,
                        (INT) LGM_BNDS_OBJECT((LGM_BNDS *)data),
                        LGM_OBJECT_ID_2_OBJECT((INT) LGM_BNDS_OBJECT((LGM_BNDS *)data))));

      LGM_BNDS_OBJECT(bs) =
        LGM_OBJECT_ID_2_OBJECT((INT) LGM_BNDS_OBJECT((LGM_BNDS *)data));
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

  size    = CEIL(sizeof(INT))+CEIL(n*sizeof(LGM_BNDP_POBJECT));

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
  size    = CEIL(n*sizeof(LGM_BNDP_POBJECT));

  PRINTDEBUG(dom,1,(PFMT "BVertexGatherBndP(): bndp=%x n=%d size=%d\n"
                    ,me,bndp,n,size));

  memcpy(data,&n,sizeof(INT));
  data += CEIL(sizeof(INT));
  memcpy(data,&LGM_BNDP_OBJECTS(bndp,0),size);

  /* substitute the line/surface pointer by their ids */
  for (i=0; i<n; i++)
  {
    PRINTDEBUG(dom,3,(PFMT "BVertexGatherBndP(): substituting i=%d "
                      "line/surfaceptr=%x by line/surfaceid=%d local=%11.4E\n",
                      me,i,LGM_BNDP_OBJECT(bndp,i),
                      LGM_OBJECT_ID(LGM_BNDP_OBJECT(bndp,i)),
                      LGM_BNDP_LOCAL(bndp,i)));
    fflush(stdout);

    LGM_BNDP_OBJECT_GOBJECT(((LGM_BNDP_POBJECT *)data)[0]) =
      (LGM_OBJECT *)LGM_OBJECT_ID(LGM_BNDP_OBJECT(bndp,i));

    PRINTDEBUG(dom,3,(PFMT "BVertexGatherBndP(): inside buffer: "
                      "line=%x local=%11.4E\n",me,
                      LGM_BNDP_OBJECT_GOBJECT(((LGM_BNDP_POBJECT *)data)[0]),
                      LGM_BNDP_OBJECT_LOCAL(((LGM_BNDP_POBJECT *)data)[0])));

    data += sizeof(LGM_BNDP_POBJECT);
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
    INT size        = CEIL(n*sizeof(LGM_BNDP_POBJECT));

    *bndp = (LGM_BNDP *) memmgr_AllocOMEM(sizeof(LGM_BNDP)+(n-1)*
                                          sizeof(LGM_BNDP_POBJECT),TypeBndP,0,0);
    ASSERT(*bndp!=NULL);

    PRINTDEBUG(dom,2,(PFMT "BVertexScatterBndP(): scatter bndp=%x n=%d "
                      "size %d\n",me,*bndp,n,size))

    memcpy(&LGM_BNDP_N(*bndp),&n,sizeof(INT));
    data += CEIL(sizeof(INT));
    memcpy(&LGM_BNDP_OBJECTS(*bndp,0),data,size);

    /* substitute the line ids by their pointers */
    for (i=0; i<n; i++)
    {
      /* substitute line id by its local pointer */
      PRINTDEBUG(dom,3,(PFMT "BVertexScatterBndP(): substituting i=%d "
                        "lineid=%d by lineptr=%x local=%11.4E\n",
                        me,i,(INT) LGM_BNDP_OBJECT_GOBJECT(((LGM_BNDP_POBJECT *)data)[0]),
                        LGM_OBJECT_ID_2_OBJECT((INT)
                                               LGM_BNDP_OBJECT_GOBJECT(((LGM_BNDP_POBJECT *)data)[0])),
                        LGM_BNDP_OBJECT_LOCAL(((LGM_BNDP_POBJECT *)data)[0])));

      LGM_BNDP_OBJECT(*bndp,i) =
        LGM_OBJECT_ID_2_OBJECT((INT)
                               LGM_BNDP_OBJECT_GOBJECT(((LGM_BNDP_POBJECT *)data)[0]));

      PRINTDEBUG(dom,3,(PFMT "BVertexScatterBndP(): outside buffer: "
                        "line=%x local=%11.4E\n",me,
                        LGM_BNDP_OBJECT(*((LGM_BNDP **)Bndp),i),
                        LGM_BNDP_LOCAL(*((LGM_BNDP **)Bndp),i)));

      data += sizeof(LGM_BNDP_POBJECT);
    }
  }
  else
  {
    PRINTDEBUG(dom,2,(PFMT "BVertexScatterBndP(): ignoring bndp=%x",
                      me,*bndp))
  }
}

#endif /* ModelP */
