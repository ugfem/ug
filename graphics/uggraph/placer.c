// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      placer.c                                                      */
/*                                                                          */
/* Purpose:   automatic placement for pictures inside a window.             */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: birken@ica3.uni-stuttgart.de                        */
/*                                                                          */
/* History:   971107 kb  begin                                              */
/*                                                                          */
/* Remarks:   The algorithm for placing pictures has two levels:            */
/*            The outer level applies an evolutionary strategy (i.e., a     */
/*            simplified genetic algorithm) for optimization of a given     */
/*            order. The order consists of the pictures which must be       */
/*            placed and 'newline' codes. The inner level places the        */
/*            pictures in a given order into the window. It starts at       */
/*            the lower left and proceeds to the right. Whenever a 'newline'*/
/*            code is found, the input focus moves one step up and to the   */
/*            left. During that procedure the Place-algorithm maintains     */
/*            a staircase-like border.                                      */
/*                                                                          */
/****************************************************************************/

#define INSIDE_UG


/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifdef INSIDE_UG
#include "general.h"
#include "wpm.h"
#ifndef __COMPILER__
#include "compiler.h"
#endif
#else
typedef int INT;
typedef double DOUBLE;
#endif


/****************************************************************************/
/*                                                                          */
/* constants, macros, data types                                            */
/*                                                                          */
/****************************************************************************/


/* this is the relative weight of input-order optimization vs. placement */
#define INPUT_ORDER_WEIGHT 0.02


#define MAX_RECTS   128
#define MAX_ORDER   (MAX_RECTS*2)
#define MAX_CORNERS (MAX_RECTS*2)


#define FTOI(f)   ((int)((f)+((f)<0.0 ? -.5 : .5)))



/****************************************************************************/


typedef struct
{
  /* order in input array */
  int input_order;

  /* given aspect ratio (y/x) */
  double aspect_ratio;

  /* given relative size */
  double rel_size;

  /* computed win pos and size */
  double llx, lly, sx, sy;
} PRect;

typedef PRect *PRectPtr;



typedef struct {
  double x, y;
} Corner;



/****************************************************************************/
/*                                                                          */
/* variables used in this source file only                                  */
/*                                                                          */
/****************************************************************************/


/* RCS string */
#ifdef INSIDE_UG
static char RCS_ID("$Header$",UG_RCS_STRING);
#endif



/****************************************************************************/
/*                                                                          */
/* subroutines                                                              */
/*                                                                          */
/****************************************************************************/



static void InitPic (PRect *rect, int i, DOUBLE ar, DOUBLE rs)
{
  assert(ar>0.0);

  /* input values */
  rect->aspect_ratio = ar;
  rect->rel_size     = rs;
  rect->input_order  = i;

  /* dont place here */
  rect->llx = 0.0;
  rect->lly = 0.0;

  /* we scale to: sx*sy==rel_size */
  rect->sx  = sqrt(rs/ar);
  rect->sy  = ar * rect->sx;
}


static void InitWin (PRect *rect, INT *ll, INT *ur)
{
  rect->llx = (double)ll[0];
  rect->lly = (double)ll[1];
  rect->sx  = (double)(ur[0]-ll[0]);
  rect->sy  = (double)(ur[1]-ll[1]);

  rect->aspect_ratio = rect->sy / rect->sx;
  rect->rel_size = rect->sx*rect->sy;
}



/****************************************************************************/


static void ScaleIntoWindow (PRect *win, int n, PRect *rect)
{
  DOUBLE xmax=0.0, ymax=0.0;
  DOUBLE scale_factor;
  int i;

  for(i=0; i<n; i++)
  {
    if (rect[i].llx + rect[i].sx  > xmax)
      xmax = rect[i].llx + rect[i].sx;
    if (rect[i].lly + rect[i].sy  > ymax)
      ymax = rect[i].lly + rect[i].sy;
  }

  /* compute scale_factor */
  if (win->sx/xmax < win->sy/ymax)
  {
    scale_factor = win->sx/xmax;
  }
  else
  {
    scale_factor = win->sy/ymax;
  }

  for(i=0; i<n; i++)
  {
    /* scale size */
    rect[i].sx *= scale_factor;
    rect[i].sy *= scale_factor;

    rect[i].llx = (rect[i].llx*scale_factor) + win->llx;
    rect[i].lly = (rect[i].lly*scale_factor) + win->lly;
  }
}


/****************************************************************************/

static int InitOrder (PRectPtr *order, int n, PRect *rect)
{
  int i;
  for(i=0; i<n; i++)
  {
    order[i]   = &(rect[i]);
    order[i+n] = NULL;               /* NULL ptr means 'newline' in win */
  }
  return(2*n);
}


static void CopyOrder (int n, PRectPtr *order1, PRectPtr *order2)
{
  int i;
  for(i=0; i<n; i++)
    order1[i] = order2[i];
}



/****************************************************************************/

static int p1_last, p2_last;

static void Mutate (PRectPtr *order, int n)
{
  PRectPtr tmp;
  int p1, p2;

  p2 = p1 = ((int)random())%n;
  while(p2==p1)
    p2 = ((int)random())%n;

  p1_last = p1;
  p2_last = p2;

  tmp = order[p1];
  order[p1] = order[p2];
  order[p2] = tmp;
}


static void MutateBack (PRectPtr *order, int n)
{
  PRectPtr tmp;

  tmp = order[p1_last];
  order[p1_last] = order[p2_last];
  order[p2_last] = tmp;
}



/****************************************************************************/


static double Place (PRectPtr *order, int n, PRect *win)
{
  Corner corner[MAX_CORNERS];
  int i, c, cn, ci,  last_idx;
  DOUBLE input_order_dist;
  DOUBLE xmax=0.0, ymax=0.0;
  DOUBLE area_sum=0.0;
  DOUBLE scale_factor, eval_func;

  /* we start at the lower left */
  corner[0].x = 0.0;
  corner[0].y = 0.0;
  c = 0; cn = 1;

  /* the corner array stores the non-konvex corners of the
     staircase from upper left to lower right. (cn-1) is the
     number of steps in that staircase. */

  /* place rects according to order */
  for(i=0; i<n; i++)
  {
    if (order[i]!=NULL)
    {
      DOUBLE ynew, xnew;

      /* it's a rect, place it */
      order[i]->llx = corner[c].x;
      order[i]->lly = corner[c].y;

      xnew = order[i]->llx + order[i]->sx;
      ynew = order[i]->lly + order[i]->sy;

      /* handle corners to the back (=upper left) */
      if (corner[c].x==0.0 || ynew<corner[c-1].y)
      {
        /* create new corner */
        for(ci=cn; ci>c; ci--)
          corner[ci] = corner[ci-1];
        cn++; assert(cn<=MAX_CORNERS);

        /* move old corner up */
        corner[c].y = ynew;
        c++;
      }
      else
      {
        /* move previous corner up */
        corner[c-1].y = ynew;

        while (c-2>=0 && ynew>=corner[c-2].y)
        {
          /* merge corners c-1 and c-2 */
          corner[c-2].y = corner[c-1].y;

          /* eat c-1 corner */
          for(ci=c-1; ci<cn-1; ci++)
            corner[ci] = corner[ci+1];
          cn--;
          c--;
        }
      }

      /* move current x-corner */
      corner[c].x = xnew;

      /* handle corners in direction of lower right */
      while (c+1<cn && xnew>=corner[c+1].x)
      {
        /* merge current with next (move new corner down) */
        corner[c].y = corner[c+1].y;

        /* eat c+1 corner */
        for(ci=c+1; ci<cn-1; ci++)
          corner[ci] = corner[ci+1];
        cn--;
      }

      /* compute evaluation data */
      area_sum += order[i]->sx*order[i]->sy;
      if (xnew  > xmax) xmax = xnew;
      if (ynew  > ymax) ymax = ynew;
    }
    else
    {
      /* it's a newline, move one step of stair up&left */
      if (c>0) c--;
    }
  }


  /* nested optimization for input_order */
  input_order_dist = 0.0;
  last_idx = -1;
  for(i=0; i<n; i++)
  {
    if (order[i]!=NULL)
    {
      int curr_idx = order[i]->input_order;

      if (last_idx!=-1)
      {
        input_order_dist += (DOUBLE)(last_idx-curr_idx);
      }
      last_idx = curr_idx;
    }
  }
  input_order_dist = input_order_dist/(DOUBLE)(n*n/2);


  /* evaluate solution */
  if (win->sx/xmax < win->sy/ymax)
  {
    scale_factor = win->sx/xmax;
  }
  else
  {
    scale_factor = win->sy/ymax;
  }
  scale_factor = scale_factor*scale_factor;

  /* eval-function is percentage of empty area in window */
  eval_func = 1.0 - area_sum*scale_factor/(win->sx*win->sy);

  /* add optimization for input-order */
  eval_func += (input_order_dist*INPUT_ORDER_WEIGHT);


  /*
     printf("win=%lf a_sum=%lf  iod=%lf eval=%lf\n",
     win->sx*win->sy, area_sum*scale_factor,
     input_order_dist*INPUT_ORDER_WEIGHT,
     eval_func);
   */


  /* return evaluation function */
  return(eval_func);
}


/****************************************************************************/


#ifdef INSIDE_UG
INT PlacePictures (PLACEMENT_TASK *task, PLACEMENT_REAL *real)
{
  int n = task->n;
#else
int placer_Exec (INT *win_geom, int n, DOUBLE *rect_defs, INT *geom)
{
#endif

  PRect win, rect[MAX_RECTS];
  PRectPtr order[MAX_ORDER], best_order[MAX_ORDER];
  int order_len;
  int i, k, iterations;
  DOUBLE threshold, sol, sol_last, sol_best, step;

  assert(n<MAX_RECTS);

  /* transform input data */
#ifdef INSIDE_UG
  InitWin(&win, task->winLL, task->winUR);
  for(i=0; i<n; i++)
  {
    InitPic(rect+i, i, task->aspect_ratio[i], task->rel_size[i]);
  }

#else
  for(i=0; i<n; i++)
  {
    InitPic(rect+i, i, rect_defs[0], rect_defs[1]);
    rect_defs+=2;
  }
  InitWin(&win, win_geom, win_geom+2);
#endif


  /* init random generator */
  srandom(1);

  order_len = InitOrder(order, n, rect);
  CopyOrder(MAX_ORDER, best_order, order);

  /* init optimization data */
  sol_last = sol_best = Place(order, order_len, &win);
  iterations = 400*n;
  threshold = sol_last/20.0;
  step = threshold / iterations;

  /* do evolution strategy */
  for(i=0; i<iterations; i++)
  {
    /* do one optimization step */
    Mutate(order, order_len);
    sol = Place(order, order_len, &win);
    if (sol - sol_last < threshold)
    {
      sol_last = sol;
      if (sol_last < sol_best)
      {
        sol_best = sol_last;
        CopyOrder(MAX_ORDER, best_order, order);

        printf("placer: new best in iteration %d: sol=%f\n",
               i, sol_best);
      }
    }
    else
    {
      MutateBack(order, order_len);
    }
    threshold -= step;
  }

  /* reconstruct best order */
  Place(best_order, order_len, &win);


  /* for debugging only
     for(i=0; i<order_len; i++)
     {
          if (best_order[i]!=NULL)
                  printf("best_order.pic %d\n", best_order[i]->input_order);
     }
   */

  /* scale sizes and lower-lefts */
  ScaleIntoWindow(&win, n, rect);


  /* fill return array */
  for(i=0, k=0; i<n; i++, k+=4)
  {
    DOUBLE dx=0.0, dy=0.0;
    DOUBLE mirror_y = win.sy+2*win.lly;

    /* mirror y-coord */
#ifdef INSIDE_UG
    /* output relative to window LL */
    dx = -win.llx;
    dy = -win.lly;
    real->picLL[i][0] = FTOI(dx+rect[i].llx);
    real->picLL[i][1] = FTOI(dy+mirror_y-(rect[i].lly+rect[i].sy));
    real->picUR[i][0] = FTOI(dx+rect[i].llx+rect[i].sx);
    real->picUR[i][1] = FTOI(dy+mirror_y-rect[i].lly);
#else
    geom[k]   = FTOI(rect[i].llx);
    geom[k+1] = FTOI(mirror_y-(rect[i].lly+rect[i].sy));
    geom[k+2] = FTOI(rect[i].llx+rect[i].sx);
    geom[k+3] = FTOI(mirror_y-rect[i].lly);
#endif
  }


  /* return success */
        #ifdef INSIDE_UG
  return(0);
        #else
  return(1);
        #endif
}



/****************************************************************************/
