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
#include "memmgr.h"

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

static INT BVP_type = BVP_STANDARD;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

void SetBVPType(INT type)
{
  BVP_type = type;
}

void DomInitParallel (INT TypeBndP, INT TypeBndS)
{}

void DomHandlerInit (INT handlerSet)
{}

static void M_BElementXferBndS (BNDS **bnds, int n, int proc, int prio)
{
  INT size,i,size0;

  size = CEIL(sizeof(INT));
  for (i=0; i<n; i++)
    if (bnds[i] != NULL)
    {
      size0 = M_BNDS_SIZE(bnds[i]);
      size += CEIL(size0) + CEIL(sizeof(INT));

      PRINTDEBUG(dom,2,("BElementXferBndS(): Xfer  me %x "
                        "%d pid %d n %d size %d\n",
                        me,bnds[i],BND_PATCH_ID(bnds[i]),
                        (int)(((M_BNDS *)(bnds[i]))->n),size0));

    }

  DDD_XferAddData(size,DDD_DOMAIN_DATA);
}

static void M_BElementGatherBndS (BNDS **bnds, int n, int cnt, char *data)
{
  INT size,i;

  for (i=0; i<n; i++)
    if (bnds[i] != NULL)
    {
      size = M_BNDS_SIZE(bnds[i]);

      PRINTDEBUG(dom,2,("BElementGatherBndS(): %d  "
                        "me %d %x pid %d n %d size %d\n",i,
                        me,bnds[i],BND_PATCH_ID(bnds[i]),
                        (int)(((M_BNDS *)(bnds[i]))->n),size));


      memcpy(data,&i,sizeof(INT));
      data += CEIL(sizeof(INT));
      memcpy(data,bnds[i],size);
      data += CEIL(size);
    }
  i = -1;
  memcpy(data,&i,sizeof(INT));
}

static void M_BElementScatterBndS (BNDS **bnds, int n, int cnt, char *data)
{
  INT size,i;
  BNDS *bs;

  memcpy(&i,data,sizeof(INT));
  while (i != -1)
  {
    data += CEIL(sizeof(INT));
    bs = (BNDS *) data;
    size = M_BNDS_SIZE(bs);

    PRINTDEBUG(dom,1,("BElementScatterBndS(): %d me %d\n",i,size));

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

static void M_BVertexXferBndP (BNDP *bndp, int proc, int prio)
{
  INT size;

  size = sizeof(M_BNDP);

  PRINTDEBUG(dom,2,("BVertexXferBndP():  me %x %d pid %d n %d size %d\n",
                    me,bndp,BND_PATCH_ID(bndp),BND_N(bndp),BND_SIZE(bndp)));

  DDD_XferAddData(size,DDD_DOMAIN_DATA);
}

static void M_BVertexGatherBndP (BNDP *bndp, int cnt, char *data)
{
  PRINTDEBUG(dom,2,("BVertexGatherBnd():  me %d pid %d "
                    "n %d size %d cnt %d\n",
                    me,BND_PATCH_ID(bndp),
                    BND_N(bndp),BND_SIZE(bndp),cnt));

  ASSERT(cnt == sizeof(M_BNDP));

  memcpy(data,bndp,cnt);
}

static void M_BVertexScatterBndP (BNDP **bndp, int cnt, char *data)
{
  if (*bndp == NULL)
  {
    *bndp = (BNDP *) memmgr_AllocOMEM((size_t)cnt,TypeBndP,0,0);
    memcpy(*bndp,data,cnt);
    PRINTDEBUG(dom,2,("BVertexScatterBndP():  me %d pid "
                      "%d n %d size %d cnt %d\n",
                      me,BND_PATCH_ID(*bndp),
                      BND_N(*bndp),BND_SIZE(*bndp),cnt));
  }
}

void BElementXferBndS (BNDS **bnds, int n, int proc, int prio)
{
  INT size,i,size0;

  if (BVP_type == BVP_MARC) {
    M_BElementXferBndS(bnds,n,proc,prio);
    return;
  }
  size = CEIL(sizeof(INT));
  for (i=0; i<n; i++)
    if (bnds[i] != NULL)
    {
      size0 = BND_SIZE(bnds[i]);
      size += CEIL(size0) + CEIL(sizeof(INT));

      PRINTDEBUG(dom,1,("BElementXferBndS(): Xfer  me %x "
                        "%d pid %d n %d size %d\n",
                        me,bnds[i],BND_PATCH_ID(bnds[i]),
                        BND_N(bnds[i]),BND_SIZE(bnds[i])));

    }

  DDD_XferAddData(size,DDD_DOMAIN_DATA);
}

void BElementGatherBndS (BNDS **bnds, int n, int cnt, char *data)
{
  INT size,i;

  if (BVP_type == BVP_MARC) {
    M_BElementGatherBndS(bnds,n,cnt,data);
    return;
  }
  for (i=0; i<n; i++)
    if (bnds[i] != NULL)
    {

      PRINTDEBUG(dom,1,("BElementGatherBndS(): %d  "
                        "me %d %x pid %d n %d size %d\n",i,
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

  if (BVP_type == BVP_MARC) {
    M_BElementScatterBndS(bnds,n,cnt,data);
    return;
  }
  memcpy(&i,data,sizeof(INT));
  while (i != -1)
  {
    data += CEIL(sizeof(INT));
    bs = (BNDS *) data;
    size = BND_SIZE(bs);

    PRINTDEBUG(dom,1,("BElementScatterBndS(): %d me %d\n",i,size));

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

  if (BVP_type == BVP_MARC) {
    M_BVertexXferBndP(bndp,proc,prio);
    return;
  }
  size = BND_SIZE(bndp);

  PRINTDEBUG(dom,1,("BVertexXferBndP():  me %x %d pid %d n %d size %d\n",
                    me,bndp,BND_PATCH_ID(bndp),BND_N(bndp),BND_SIZE(bndp)));

  DDD_XferAddData(size,DDD_DOMAIN_DATA);
}

void BVertexGatherBndP (BNDP *bndp, int cnt, char *data)
{
  if (BVP_type == BVP_MARC) {
    M_BVertexGatherBndP(bndp,cnt,data);
    return;
  }
  PRINTDEBUG(dom,1,("BVertexGatherBnd():  me %d pid %d "
                    "n %d size %d cnt %d\n",
                    me,BND_PATCH_ID(bndp),
                    BND_N(bndp),BND_SIZE(bndp),cnt));

  ASSERT(cnt == BND_SIZE(bndp));

  memcpy(data,bndp,cnt);
}


void BVertexScatterBndP (BNDP **bndp, int cnt, char *data)
{
  if (BVP_type == BVP_MARC) {
    M_BVertexScatterBndP(bndp,cnt,data);
    return;
  }
  if (*bndp == NULL)
  {
    *bndp = (BNDS *) memmgr_AllocOMEM((size_t)cnt,TypeBndP,0,0);
    memcpy(*bndp,data,cnt);
    PRINTDEBUG(dom,1,("BVertexScatterBndP():  me %d pid "
                      "%d n %d size %d cnt %d\n",
                      me,BND_PATCH_ID(*bndp),
                      BND_N(*bndp),BND_SIZE(*bndp),cnt));
  }
}
#endif /* ModelP */
