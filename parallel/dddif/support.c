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
#if defined(__SR2201__) || defined(__CC__)
#include <stdlib.h>
#include <stdio.h>
#endif

#include <stddef.h>

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

/****************************************************************************/
/*D
   UG_GlobalMaxINT - get maximum for INT value

   SYNOPSIS:
   INT UG_GlobalMaxINT (INT i)

   PARAMETERS:
   .  i - calculate maximum for this value

   DESCRIPTION:
   This function calculates the maximum of i over all processors.

   RETURN VALUE:
   The maximum value of i over all processors.

   D*/
/****************************************************************************/

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

/****************************************************************************/
/*D
   UG_GlobalMinINT - get global minimum for INT value

   SYNOPSIS:
   INT UG_GlobalMaxINT (INT i)

   PARAMETERS:
   .  i - calculate minimum for this value

   DESCRIPTION:
   This function calculates the minimum of i over all processors.

   RETURN VALUE:
   The minimum value of i over all processors

   D*/
/****************************************************************************/

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

/****************************************************************************/
/*D
   UG_GlobalSumINT - get global sum for INT value

   SYNOPSIS:
   INT UG_GlobalSumINT (INT i)

   PARAMETERS:
   .  i - calculate sum for this variable

   DESCRIPTION:
   This function calculates the sum of i over all processors.

   RETURN VALUE:
   The sum of i over all processors

   D*/
/****************************************************************************/

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

/****************************************************************************/
/*D
   UG_GlobalMaxNINT - get maximum for n integer values

   SYNOPSIS:
   void UG_GlobalMaxNINT (INT n, INT *x)

   PARAMETERS:
   .  n - number of elements in array x to be used
   .  x - array of size n

   DESCRIPTION:
   This function calculates the maximum of x[i] over all processors for all
   i from 0 to n-1. x is overwritten with the maximum value.

   RETURN VALUE:
   none

   D*/
/****************************************************************************/

void UG_GlobalMaxNINT (INT n, INT *x)
{
  int i,l,size;
  INT *y;

  size = sizeof(INT)*n;
  y = (INT *)memmgr_AllocTMEM(size, TMEM_STD);

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,y,size);
    for (i=0; i<n; i++)
      x[i] = MAX(x[i],y[i]);
  }
  Concentrate(x,size);
  Broadcast(x,size);

  memmgr_FreeTMEM(y, TMEM_STD);

  return;
}

/****************************************************************************/
/*D
   UG_GlobalMinNINT - get minimum for n integer values

   SYNOPSIS:
   void UG_GlobalMinNINT (INT n, INT *x)

   PARAMETERS:
   .  n - number of elements in array x to be used
   .  x - array of size n

   DESCRIPTION:
   This function calculates the minimum of x[i] over all processors for all
   i from 0 to n-1. x is overwritten with the minimum value.

   RETURN VALUE:
   none

   D*/
/****************************************************************************/

void UG_GlobalMinNINT (INT n, INT *x)
{
  int i,l,size;
  INT *y;

  size = sizeof(INT)*n;
  y = (INT *)memmgr_AllocTMEM(size, TMEM_STD);

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,y,size);
    for (i=0; i<n; i++)
      x[i] = MIN(x[i],y[i]);
  }
  Concentrate(x,size);
  Broadcast(x,size);

  memmgr_FreeTMEM(y, TMEM_STD);

  return;
}

/****************************************************************************/
/*D
   UG_GlobalSumNINT - calculate global sum for n integer values

   SYNOPSIS:
   void UG_GlobalSumNINT (INT n, INT *x)

   PARAMETERS:
   .  n - number of elements in array x to be used
   .  x - array of size n

   DESCRIPTION:
   This function calculates the sum of x[i] over all processors for each
   i from 0 to n-1. x is overwritten with the result.

   RETURN VALUE:
   none

   D*/
/****************************************************************************/

void UG_GlobalSumNINT (INT n, INT *xs)
{
  int l, i, size=sizeof(INT)*n;
  INT *ys;

  ys = (INT *)memmgr_AllocTMEM(size, TMEM_STD);

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,ys,size);

    for(i=0; i<n; i++)
      xs[i] += ys[i];
  }
  Concentrate(xs,size);
  Broadcast(xs,size);

  memmgr_FreeTMEM(ys, TMEM_STD);
}

/****************************************************************************/
/*D
   UG_GlobalMaxDOUBLE - get global maximum for DOUBLE value

   SYNOPSIS:
   DOUBLE UG_GlobalMaxDOUBLE (DOUBLE x)

   PARAMETERS:
   .  x - calculate maximum for this value

   DESCRIPTION:
   This function calculates the maximum of x over all processors.

   RETURN VALUE:
   The maximum value of x over all processors.

   D*/
/****************************************************************************/

DOUBLE UG_GlobalMaxDOUBLE (DOUBLE x)
{
  int l;
  DOUBLE n;

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,&n,sizeof(DOUBLE));
    x = MAX(x,n);
  }
  Concentrate(&x,sizeof(DOUBLE));
  Broadcast(&x,sizeof(DOUBLE));
  return(x);
}

/****************************************************************************/
/*D
   UG_GlobalMinDOUBLE - get global minimum for DOUBLE value

   SYNOPSIS:
   DOUBLE UG_GlobalMinDOUBLE (DOUBLE x)

   PARAMETERS:
   .  x - calculate minimum for this value

   DESCRIPTION:
   This function calculates the minimum of x over all processors.

   RETURN VALUE:
   The minimum value of x over all processors.

   D*/
/****************************************************************************/

DOUBLE UG_GlobalMinDOUBLE (DOUBLE x)
{
  int l;
  DOUBLE y;

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,&y,sizeof(DOUBLE));
    x = MIN(x,y);
  }
  Concentrate(&x,sizeof(DOUBLE));
  Broadcast(&x,sizeof(DOUBLE));
  return(x);
}

/****************************************************************************/
/*D
   UG_GlobalSumDOUBLE - get global sum for DOUBLE value

   SYNOPSIS:
   DOUBLE UG_GlobalSumDOUBLE (DOUBLE i)

   PARAMETERS:
   .  i - calculate sum for this variable

   DESCRIPTION:
   This function calculates the sum of i over all processors.

   RETURN VALUE:
   The sum of i over all processors

   D*/
/****************************************************************************/

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

/****************************************************************************/
/*D
   UG_GlobalMaxNDOUBLE - get maximum over n integer values

   SYNOPSIS:
   void UG_GlobalMaxNDOUBLE (INT n, DOUBLE *x)

   PARAMETERS:
   .  n - number of elements in array x to be used
   .  x - array of size n

   DESCRIPTION:
   This function calculates the maximum of x[i] over all processors for all
   i from 0 to n-1. x is overwritten with the maximum value.

   RETURN VALUE:
   none

   D*/
/****************************************************************************/

void UG_GlobalMaxNDOUBLE (INT n, DOUBLE *x)
{
  int i,l,size;
  DOUBLE *y;

  size = sizeof(DOUBLE)*n;
  y = (DOUBLE *)memmgr_AllocTMEM(size, TMEM_STD);

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,y,size);
    for (i=0; i<n; i++)
      x[i] = MAX(x[i],y[i]);
  }
  Concentrate(x,size);
  Broadcast(x,size);

  memmgr_FreeTMEM(y, TMEM_STD);

  return;
}

/****************************************************************************/
/*D
   UG_GlobalMinNDOUBLE - get minimum over n integer values

   SYNOPSIS:
   void UG_GlobalMinNDOUBLE (INT n, DOUBLE *x)

   PARAMETERS:
   .  n - number of elements in array x to be used
   .  x - array of size n

   DESCRIPTION:
   This function calculates the minimum of x[i] over all processors for all
   i from 0 to n-1. x is overwritten with the minimum value.

   RETURN VALUE:
   none

   D*/
/****************************************************************************/

void UG_GlobalMinNDOUBLE (INT n, DOUBLE *x)
{
  int i,l,size;
  DOUBLE *y;

  size = sizeof(DOUBLE)*n;
  y = (DOUBLE *)memmgr_AllocTMEM(size, TMEM_STD);

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,y,size);
    for (i=0; i<n; i++)
      x[i] = MIN(x[i],y[i]);
  }
  Concentrate(x,size);
  Broadcast(x,size);

  memmgr_FreeTMEM(y, TMEM_STD);

  return;
}

/****************************************************************************/
/*D
   UG_GlobalSumNDOUBLE - calculate global sum for n DOUBLE values

   SYNOPSIS:
   void UG_GlobalSumNDOUBLE (INT n, DOUBLE *x)

   PARAMETERS:
   .  n - number of elements in array x to be used
   .  x - array of size n

   DESCRIPTION:
   This function calculates the sum of x[i] over all processors for each
   i from 0 to n-1. x is overwritten with the result.

   RETURN VALUE:
   none

   D*/
/****************************************************************************/

void UG_GlobalSumNDOUBLE (INT n, DOUBLE *x)
{
  int l, i, size=sizeof(DOUBLE)*n;
  DOUBLE *y;

  y = (DOUBLE *)memmgr_AllocTMEM(size, TMEM_STD);

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,y,size);
    for(i=0; i<n; i++)
      x[i] += y[i];
  }
  Concentrate(x,size);
  Broadcast(x,size);

  memmgr_FreeTMEM(y, TMEM_STD);
}
