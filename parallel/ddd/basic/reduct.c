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


int ddd_GlobalMaxInt (int i)
{
  int l,n;

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,&n,sizeof(int));
    i = MAX(i,n);
  }
  Concentrate(&i,sizeof(int));
  Broadcast(&i,sizeof(int));
  return(i);
}

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



/*
   double ddd_GlobalMaxDouble (double i)
   {
    int l;
        double n;

    for (l=degree-1; l>=0; l--)
    {
        GetConcentrate(l,&n,sizeof(double));
        i = MAX(i,n);
    }
    Concentrate(&i,sizeof(double));
    Broadcast(&i,sizeof(double));
    return(i);
   }

   double ddd_GlobalMinDouble (double i)
   {
    int l;
        double n;

    for (l=degree-1; l>=0; l--)
    {
        GetConcentrate(l,&n,sizeof(double));
        i = MIN(i,n);
    }
    Concentrate(&i,sizeof(double));
    Broadcast(&i,sizeof(double));
    return(i);
   }

   double ddd_GlobalSumDouble (double x)
   {
        int l;
        double y;

        for (l=degree-1; l>=0; l--)
        {
                GetConcentrate(l,&y,sizeof(double));
                x += y;
        }
        Concentrate(&x,sizeof(double));
        Broadcast(&x,sizeof(double));
        return(x);
   }
 */
