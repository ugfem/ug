// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      reduct.c                                                      */
/*                                                                          */
/* Purpose:   standard parallel routines not supported by ddd               */
/*            reduction operations (GlobalSum, GlobalMax etc)               */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   940128 kb  begin                                              */
/*            960902 kb  copied from fedemo, adapted                        */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

/* standard C library */
#ifdef __SR2201__
#include <stdlib.h>
#include <stdio.h>
#endif

#include "general.h"
#include "ppif.h"
#include "compiler.h"
#include "memmgr.h"


/****************************************************************************/
/*                                                                          */
/* macros                                                                   */
/*                                                                          */
/****************************************************************************/


#ifndef MAX
#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#endif


/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


/* some useful functions by Peter Bastian, from ugp/ug/ugcom.c */

INT UG_GlobalMaxINT (INT i)
{
  int l;
  INT n;

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,&n,sizeof(INT));
    i = MAX(i,n);
  }
  Concentrate(&i,sizeof(INT));
  Broadcast(&i,sizeof(INT));
  return(i);
}

INT UG_GlobalMinINT (INT i)
{
  int l;
  INT n;

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,&n,sizeof(INT));
    i = MIN(i,n);
  }
  Concentrate(&i,sizeof(INT));
  Broadcast(&i,sizeof(INT));
  return(i);
}

DOUBLE UG_GlobalMaxDOUBLE (DOUBLE i)
{
  int l;
  DOUBLE n;

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,&n,sizeof(DOUBLE));
    i = MAX(i,n);
  }
  Concentrate(&i,sizeof(DOUBLE));
  Broadcast(&i,sizeof(DOUBLE));
  return(i);
}

DOUBLE UG_GlobalMinDOUBLE (DOUBLE i)
{
  int l;
  DOUBLE n;

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,&n,sizeof(DOUBLE));
    i = MIN(i,n);
  }
  Concentrate(&i,sizeof(DOUBLE));
  Broadcast(&i,sizeof(DOUBLE));
  return(i);
}

DOUBLE UG_GlobalSumDOUBLE (DOUBLE x)
{
  int l;
  DOUBLE y;

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,&y,sizeof(DOUBLE));
    x += y;
  }
  Concentrate(&x,sizeof(DOUBLE));
  Broadcast(&x,sizeof(DOUBLE));
  return(x);
}

INT UG_GlobalSumINT (INT x)
{
  int l;
  INT y;

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,&y,sizeof(INT));
    x += y;
  }
  Concentrate(&x,sizeof(INT));
  Broadcast(&x,sizeof(INT));
  return(x);
}

void UG_GlobalSumNDOUBLE (INT n, DOUBLE *xs)
{
  int l, i, size=sizeof(DOUBLE)*n;
  DOUBLE *ys;

  ys = (DOUBLE *)memmgr_AllocTMEM(size);

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,ys,size);

    /* execute reduction operation */
    for(i=0; i<n; i++)
      xs[i] += ys[i];
  }
  Concentrate(xs,size);
  Broadcast(xs,size);

  memmgr_FreeTMEM(ys);
}
