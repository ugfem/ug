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
/*
   #include <stdlib.h>
   #include <stdio.h>
 */

#include "ppif.h"
#include "compiler.h"


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

/*

   int ddd_GlobalMinInt (int i)
   {
        int l,n;

        for (l=degree-1; l>=0; l--)
        {
                GetConcentrate(l,&n,sizeof(int));
                i = MIN(i,n);
        }
        Concentrate(&i,sizeof(int));
        Broadcast(&i,sizeof(int));
        return(i);
   }



   int ddd_GlobalSumInt (int x)
   {
        int l;
        int y;

        for (l=degree-1; l>=0; l--)
        {
                GetConcentrate(l,&y,sizeof(int));
                x += y;
        }
        Concentrate(&x,sizeof(int));
        Broadcast(&x,sizeof(int));
        return(x);
   }


 */

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
