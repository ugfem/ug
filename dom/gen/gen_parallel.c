// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  gen_parallel.c				                                                                */
/*																			*/
/* Purpose:   parallel part of general ug domain description                            */
/*																			*/
/* Author:	  Klaus Birken / Christian Wieners								*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   Dec 8 1999 begin                                                                          */
/*																			*/
/* Remarks:   adopted from std_parallel.c                                                               */
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

#include "gen.h"

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

static INT BVP_type = BVP_GENERAL;

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

void BElementXferBndS (BNDS **bnds, int n, int proc, int prio)
{
  INT size,i,size0;

  size = CEIL(sizeof(INT));
  for (i=0; i<n; i++)
    if (bnds[i] != NULL) {
      size0 = sizeof(BS);
      size += CEIL(size0) + CEIL(sizeof(INT));
    }

  DDD_XferAddData(size,DDD_DOMAIN_DATA);
}

void BElementGatherBndS (BNDS **bnds, int n, int cnt, char *data)
{
  INT size,i;

  for (i=0; i<n; i++)
    if (bnds[i] != NULL) {
      size = sizeof(BS);
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
  while (i != -1) {
    data += CEIL(sizeof(INT));
    bs = (BNDS *) data;
    size = sizeof(BS);
    if (bnds[i] == NULL) {
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
  BP *bp = (BP *)bndp;
  INT size;

  size = sizeof(BP)+(bp->n-1)*sizeof(int);

  DDD_XferAddData(size,DDD_DOMAIN_DATA);
}

void BVertexGatherBndP (BNDP *bndp, int cnt, char *data)
{
  memcpy(data,bndp,cnt);

  ASSERT(cnt == sizeof(BP)+(((BP *)bndp)->n-1)*sizeof(int));
}

void BVertexScatterBndP (BNDP **bndp, int cnt, char *data)
{
  if (*bndp == NULL) {
    *bndp = (BNDP *) memmgr_AllocOMEM((size_t)cnt,TypeBndP,0,0);
    memcpy(*bndp,data,cnt);
  }
}
#endif /* ModelP */
