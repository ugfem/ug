// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  graph.h														*/
/*																			*/
/* Purpose:   low level plot routines and clip management					*/
/*																			*/
/* Author:	  Klaus Johannsen												*/
/*			  Institut fuer Computeranwendungen                                                     */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: ug@ica3.uni-stuttgart.de                                                    */
/*																			*/
/* History:   8.12.94 begin, ug3-version									*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "graph.h"
#include "wpm.h"
#include "devices.h"
#include "misc.h"
#include "general.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define INDEXCHAR               '/'
#define MOVECHAR                '|'
#define REL_CHARWIDTH   0.7
#define REL_INDEXSIZE   0.7

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

static OUTPUTDEVICE *CurrentOutputDevice;               /* current output device	*/
static COORD_POINT CurrCursor;                                  /* current cursor position	*/
static SHORT CurrTextSize;                                              /* current text size		*/

static DOUBLE currClipRegionMaxX;               /* corner of ViewPort having the	*/
static DOUBLE currClipRegionMaxY;               /*largest values for each component */
static DOUBLE currClipRegionMinX;               /* corner of ViewPort having the	*/
static DOUBLE currClipRegionMinY;               /*smallest values for each component*/

static COORD_POINT currClipRegionCorner[4];     /* corners of the view port */

static char buffer[256];                                                /* general purpose text buff*/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* Function:  SubtractCoordPoint											*/
/*																			*/
/* Purpose:   Calculate linearcombination:									*/
/*			  result =	a - b												*/
/*																			*/
/*																			*/
/* output:	  INT 0: ok                                                                                                     */
/*			  INT 1: error													*/
/*																			*/
/****************************************************************************/

static INT SubtractCoordPoint (COORD_POINT a, COORD_POINT b, COORD_POINT *result)
{
  (*result).x = a.x - b.x;
  (*result).y = a.y - b.y;

  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:  EuklidNormCoordPoint											*/
/*																			*/
/* Purpose:   Calculate norm of vector a									*/
/*																			*/
/* input:	  DOUBLE *a: input vector (a[0],a[1],a[2])						*/
/*			  DOUBLE *result: output scalar result[0]						*/
/*																			*/
/* output:	  INT 0: ok                                                                                                     */
/*			  INT 1: error													*/
/*																			*/
/****************************************************************************/

static INT EuklidNormCoordPoint (COORD_POINT a, DOUBLE *result)
{
  DOUBLE sum;

  sum = a.x*a.x + a.y*a.y;
  *result = (DOUBLE)sqrt( (float)sum );

  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:  ScalarProductCoordPoint										*/
/*																			*/
/* Purpose:   Calculate linearcombination:									*/
/*			  result =	a - b												*/
/*																			*/
/*																			*/
/* output:	  INT 0: ok                                                                                                     */
/*			  INT 1: error													*/
/*																			*/
/****************************************************************************/

static INT ScalarProductCoordPoint (COORD_POINT a, COORD_POINT b, DOUBLE *result)
{
  *result = a.x*b.x + a.y*b.y;

  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:  LinCombCoordPoint                                                                                         */
/*																			*/
/* Purpose:   Calculate linearcombination:									*/
/*			  result =	lambdaa * a + lambdab * b							*/
/*																			*/
/*																			*/
/* output:	  INT 0: ok                                                                                                     */
/*			  INT 1: error													*/
/*																			*/
/****************************************************************************/

static INT LinCombCoordPoint (DOUBLE lambdaa, COORD_POINT a, DOUBLE lambdab, COORD_POINT b, COORD_POINT *result)
{
  (*result).x = lambdaa*a.x + lambdab*b.x;
  (*result).y = lambdaa*a.y + lambdab*b.y;

  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:  PrepareGraph													*/
/*																			*/
/* Purpose:   set current view and its port                                                             */
/*																			*/
/* Input:	  VIEW *theView: set graph context for that view				*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

INT PrepareGraph (const PICTURE *thePicture)
{
  /* set current output device */
  CurrentOutputDevice = PIC_OUTPUTDEV(thePicture);

  /* set position of currClipRegion */
  currClipRegionMaxX = MAX(PIC_GUR(thePicture)[0],PIC_GLL(thePicture)[0]);
  currClipRegionMaxY = MAX(PIC_GUR(thePicture)[1],PIC_GLL(thePicture)[1]);
  currClipRegionMinX = MIN(PIC_GUR(thePicture)[0],PIC_GLL(thePicture)[0]);
  currClipRegionMinY = MIN(PIC_GUR(thePicture)[1],PIC_GLL(thePicture)[1]);

  currClipRegionCorner[0].x = currClipRegionMinX;
  currClipRegionCorner[1].x = currClipRegionMaxX;
  currClipRegionCorner[2].x = currClipRegionMaxX;
  currClipRegionCorner[3].x = currClipRegionMinX;
  currClipRegionCorner[0].y = currClipRegionMaxY;
  currClipRegionCorner[1].y = currClipRegionMaxY;
  currClipRegionCorner[2].y = currClipRegionMinY;
  currClipRegionCorner[3].y = currClipRegionMinY;

  /* activate IF Window */
  if ((*CurrentOutputDevice->ActivateOutput)(UGW_IFWINDOW(PIC_UGW(thePicture)))) return (1);

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  PrepareGraphWindow											*/
/*																			*/
/* Purpose:   set current view and its port                                                             */
/*																			*/
/* Input:	  VIEW *theView: set graph context for that view				*/
/*																			*/
/* Output:	  INT 0: ok                                                                                                     */
/*				  1: error													*/
/*																			*/
/****************************************************************************/

INT PrepareGraphWindow (const UGWINDOW *theWindow)
{
  /* set current output device */
  CurrentOutputDevice = UGW_OUTPUTDEV(theWindow);

  /* set position of currClipRegion */
  currClipRegionMaxX = MAX(UGW_LUR(theWindow)[0],UGW_LLL(theWindow)[0]);
  currClipRegionMaxY = MAX(UGW_LUR(theWindow)[1],UGW_LLL(theWindow)[1]);
  currClipRegionMinX = MIN(UGW_LUR(theWindow)[0],UGW_LLL(theWindow)[0]);
  currClipRegionMinY = MIN(UGW_LUR(theWindow)[1],UGW_LLL(theWindow)[1]);

  currClipRegionCorner[0].x = currClipRegionMinX;
  currClipRegionCorner[1].x = currClipRegionMaxX;
  currClipRegionCorner[2].x = currClipRegionMaxX;
  currClipRegionCorner[3].x = currClipRegionMinX;
  currClipRegionCorner[0].y = currClipRegionMaxY;
  currClipRegionCorner[1].y = currClipRegionMaxY;
  currClipRegionCorner[2].y = currClipRegionMinY;
  currClipRegionCorner[3].y = currClipRegionMinY;

  /* activate IF Window */
  if ((*CurrentOutputDevice->ActivateOutput)(UGW_IFWINDOW(theWindow))) return (1);

  return (0);
}

/****************************************************************************/
/*																			*/
/* Function:  ClipPoint                                                                                                         */
/*																			*/
/* Purpose:   clip point against current view port							*/
/*																			*/
/* In/Output: COORD_POINT in: input point									*/
/*			  out: output point                                                                                     */
/*			  INT *reject: 1 if point has to be rejected					*/
/*						   0 if point is in view port						*/
/*																			*/
/* RETURNCODE: INT 0: ok													*/
/*			   INT 1: an error occurred                                                                     */
/*																			*/
/****************************************************************************/

static INT ClipPoint (COORD_POINT in, SHORT_POINT *out, INT *reject)
{
  INT flags;

  *reject = 1;
  flags = 0;
  if (in.y>currClipRegionMaxY) flags |= 1;
  if (in.x>currClipRegionMaxX) flags |= 2;
  if (in.y<currClipRegionMinY) flags |= 4;
  if (in.x<currClipRegionMinX) flags |= 8;
  if (!flags)
  {
    COPY_SC_TO_SH(in,*out);
    *reject = 0;
  }
  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:  ClipLine														*/
/*																			*/
/* Purpose:   clip line against current view port							*/
/*																			*/
/* In/Output: COORD_POINT in1: first input point							*/
/*			  COORD_POINT in2: second input point							*/
/*			  SHORT_POINT *out1: first output point                                                 */
/*			  SHORT_POINT *out2: second output point						*/
/*			  INT *reject: line has to be rejected							*/
/*			  INT *side1 :-1  if point 1 is in the interior                                 */
/*						   0  if point 1 is replaced by one on south bnd	*/
/*						   1  - - - - - - - - - - - - - - - -  east  - -	*/
/*						   2  - - - - - - - - - - - - - - - -  north - -	*/
/*						   3  - - - - - - - - - - - - - - - -  west  - -	*/
/*																			*/
/*			  INT *side2 : same as for point 1								*/
/*																			*/
/* RETURNCODE: INT 0: ok													*/
/*			   INT 1: an error occurred                                                                     */
/*																			*/
/* Remark: the points in1, in2 remain unchanged during the proceedure		*/
/*																			*/
/****************************************************************************/

#define SMALL_DOUBLE            1E-30

static INT ClipLine (COORD_POINT in1, COORD_POINT in2,
                     SHORT_POINT *out1,SHORT_POINT *out2,
                     INT *reject,INT *side1,INT *side2   )
{
  INT flags1,flags2,flags3,flags4;
  DOUBLE dx,dy,slope,invslope;

  flags1 = 0;
  if (in1.y>currClipRegionMaxY) flags1 |= 1;                    /* flag sceme :                                                                                 */
  if (in1.x>currClipRegionMaxX) flags1 |= 2;                    /*														*/
  if (in1.y<currClipRegionMinY) flags1 |= 4;                    /*	1100  |  0100  |  0110			12	|	4  |   6	*/
  if (in1.x<currClipRegionMinX) flags1 |= 8;                    /*	------+--------+-------                 ----+------+-----	*/
  flags2 = 0;                                                                           /*	1000  |  0000  |  0010	  or	 8	|	0  |   2	*/
  if (in2.y>currClipRegionMaxY) flags2 |= 1;                    /*	------+--------+-------                 ----+------+-----	*/
  if (in2.x>currClipRegionMaxX) flags2 |= 2;                    /*	1001  |  0001  |  0011			 9	|	1  |   3	*/
  if (in2.y<currClipRegionMinY) flags2 |= 4;                    /*														*/
  if (in2.x<currClipRegionMinX) flags2 |= 8;                    /*	(3)     2	  (2)									*/
  *side1 = -1;                                                                          /*		---------			0,1,2,3 id of side			*/
  *side2 = -1;                                                                          /*	3  |		 |	1									*/
  /*		---------			(0),(1),..id of corner		*/
  /* reject or accept the line */                               /*	(0)     0	  (1)									*/
  if ((flags1 & flags2) != 0)
  {
    *reject = 1;
    return(0);
  }
  *reject = 0;
  if ( flags1==0 && flags2==0)
  {
    COPY_SC_TO_SH(in1,*out1);
    COPY_SC_TO_SH(in2,*out2);
    return(0);
  }

  /* get differences dx, dy */
  dx = in1.x - in2.x;
  dy = in1.y - in2.y;
  flags3 = (ABS(dx) < SMALL_DOUBLE);
  flags3 |= ((ABS(dy) < SMALL_DOUBLE)<<1);
  switch(flags3)
  {
  case 0 :                      /* line is not parallel to one of the axis	 */
    slope    = dy/dx;
    invslope = dx/dy;

    /* clip to south boundary of the view port */
    flags4 = (in1.y>currClipRegionMaxY);
    flags4 |= ((in2.y>currClipRegionMaxY)<<1);
    switch(flags4)
    {
    case 0 :
      break;
    case 1 :
      *side1 = 0;
      in1.x += (currClipRegionMaxY - in1.y)*invslope;
      in1.y = currClipRegionMaxY;
      break;
    case 2 :
      *side2 = 0;
      in2.x += (currClipRegionMaxY - in2.y)*invslope;
      in2.y = currClipRegionMaxY;
      break;
    case 3 :
      *reject = 1;
      return(0);
    default :
      return(1);
    }

    /* clip to east boundary of the view port */
    flags4 = in1.x>currClipRegionMaxX;
    flags4 |= ((in2.x>currClipRegionMaxX)<<1);
    switch(flags4)
    {
    case 0 :
      break;
    case 1 :
      *side1 = 1;
      in1.y += (currClipRegionMaxX - in1.x)*slope;
      in1.x = currClipRegionMaxX;
      break;
    case 2 :
      *side2 = 1;
      in2.y += (currClipRegionMaxX - in2.x)*slope;
      in2.x = currClipRegionMaxX;
      break;
    case 3 :
      *reject = 1;
      return(0);
    default :
      return(1);
    }

    /* clip to north boundary of the view port */
    flags4 = in1.y<currClipRegionMinY;
    flags4 |= ((in2.y<currClipRegionMinY)<<1);
    switch(flags4)
    {
    case 0 :
      break;
    case 1 :
      *side1 = 2;
      in1.x += (currClipRegionMinY - in1.y)*invslope;
      in1.y = currClipRegionMinY;
      break;
    case 2 :
      *side2 = 2;
      in2.x += (currClipRegionMinY - in2.y)*invslope;
      in2.y = currClipRegionMinY;
      break;
    case 3 :
      *reject = 1;
      return(0);
    default :
      return(1);
    }

    /* clip to west boundary of the view port */
    flags4 = in1.x<currClipRegionMinX;
    flags4 |= ((in2.x<currClipRegionMinX)<<1);
    switch(flags4)
    {
    case 0 :
      break;
    case 1 :
      *side1 = 3;
      in1.y += (currClipRegionMinX - in1.x)*slope;
      in1.x = currClipRegionMinX;
      break;
    case 2 :
      *side2 = 3;
      in2.y += (currClipRegionMinX - in2.x)*slope;
      in2.x = currClipRegionMinX;
      break;
    case 3 :
      *reject = 1;
      return(0);
    default :
      return(1);

    }
    break;

  case 1 :                      /* line is parallel to Y axis				 */
    switch(flags1)
    {
    case 0 :
      break;
    case 1 :
      in1.y = currClipRegionMaxY;
      *side1 = 0;
      break;
    case 4 :
      in1.y = currClipRegionMinY;
      *side1 = 2;
      break;
    default :
      return(1);
    }
    switch(flags2)
    {
    case 0 :
      break;
    case 1 :
      in2.y = currClipRegionMaxY;
      *side2 = 0;
      break;
    case 4 :
      in2.y = currClipRegionMinY;
      *side2 = 2;
      break;
    default :
      return(1);
    }
    break;

  case 2 :                      /* line is parallel to X axis				 */
    switch(flags1)
    {
    case 0 :
      break;
    case 2 :
      in1.x = currClipRegionMaxX;
      *side1 = 1;
      break;
    case 8 :
      in1.x = currClipRegionMinX;
      *side1 = 3;
      break;
    default :
      return(1);
    }
    switch(flags2)
    {
    case 0 :
      break;
    case 2 :
      in2.x = currClipRegionMaxX;
      *side2 = 1;
      break;
    case 8 :
      in2.x = currClipRegionMinX;
      *side2 = 3;
      break;
    default :
      return(1);
    }
    break;

  default :
    return(1);
  }
  COPY_SC_TO_SH(in1,*out1);
  COPY_SC_TO_SH(in2,*out2);
  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:  Clipon												                */
/*																			*/
/* Purpose:   clip line against current view port							*/
/*																			*/
/* Input: COORD_POINT *in, INT nin,                                                                             */
/*		  SHORT_POINT *out, INT *nout										*/
/*																			*/
/* Output: INT 0: ok														*/
/*		   INT 1: an error occurred                                                                             */
/*																			*/
/* Remark: the input polygon remains unchanged during the proceedure		*/
/*																			*/
/****************************************************************************/

#define CLOCKWISE                       0
#define COUNTERCLOCKWISE        1
#define PIXEL_DIFFERENCE        10

static void FillUp (INT a, INT b, INT orientation, SHORT_POINT *out, INT *nptr)
{
  INT k;

  if (orientation == COUNTERCLOCKWISE)
    for (k=(a+1)%4; k!=(b+1)%4; k=(k+1)%4)
    {
      out[*nptr].x = (short)currClipRegionCorner[k].x;
      out[*nptr].y = (short)currClipRegionCorner[k].y;
      (*nptr)++;
    }
  else
    for (k=a; k!=b; k=((k>0) ? (k-1) : (3)))
    {
      out[*nptr].x = (short)currClipRegionCorner[k].x;
      out[*nptr].y = (short)currClipRegionCorner[k].y;
      (*nptr)++;
    }
}

static void FillUpTot (SHORT_POINT *out, INT *nptr)
{
  INT k;

  for (k=0; k<4; k++)
  {
    COPY_SC_TO_SH(currClipRegionCorner[k],out[k]);
  }
  *nptr = 4;
}

static INT ClipPolygon (COORD_POINT *in, INT nin,
                        SHORT_POINT *out, INT *nout)
{
  INT i,FillupStart,side1,side2,FirstSide,reject,flag,left,leftmark,right,rightmark,orientation;
  SHORT_POINT out1,out2;
  COORD_POINT point[3];
  DOUBLE norm, lambda, ScalarPrd;

  /* initializations */
  *nout = 0;

  /* check if polygon is degenerated */
  if (nin<3) return (0);

  /* decide whether polygon is left or right handed */
  left = right = orientation = 0;
  for (i=0; i<nin; i++)
  {
    if (((in[i].x-in[(i-1+nin)%nin].x)*(in[(i+1)%(nin)].y-in[i].y)) >=
        ((in[i].y-in[(i-1+nin)%nin].y)*(in[(i+1)%(nin)].x-in[i].x)))
      left++;
    else
      leftmark = i;

    if (((in[i].x-in[(i-1+nin)%nin].x)*(in[(i+1)%(nin)].y-in[i].y)) <=
        ((in[i].y-in[(i-1+nin)%nin].y)*(in[(i+1)%(nin)].x-in[i].x)))
      right++;
    else
      rightmark = i;
  }
  if (left == nin)
    orientation = CLOCKWISE;
  else if (right == nin)
    orientation = COUNTERCLOCKWISE;
  else if (left == nin-1)
  {
    SubtractCoordPoint(in[leftmark],in[(leftmark-1+nin)%nin],&(point[0]));
    SubtractCoordPoint(in[(leftmark+1)%nin],in[(leftmark-1+nin)%nin],&(point[1]));
    EuklidNormCoordPoint(point[1],&norm);
    ScalarProductCoordPoint(point[0],point[1],&ScalarPrd);
    if (norm < SMALL_C)
      lambda = 1.0;
    else
      lambda = ScalarPrd/norm/norm;
    LinCombCoordPoint(1.0, in[(leftmark-1+nin)%nin], lambda, point[1], &(point[2]));
    SubtractCoordPoint(in[leftmark],point[2],&(point[0]));
    EuklidNormCoordPoint(point[0],&norm);
    in[leftmark] = point[2];
    orientation = CLOCKWISE;
  }
  else if (right == nin-1)
  {
    SubtractCoordPoint(in[rightmark],in[(rightmark-1+nin)%nin],&(point[0]));
    SubtractCoordPoint(in[(rightmark+1)%nin],in[(rightmark-1+nin)%nin],&(point[1]));
    EuklidNormCoordPoint(point[1],&norm);
    ScalarProductCoordPoint(point[0],point[1],&ScalarPrd);
    if (norm < SMALL_C)
      lambda = 1.0;
    else
      lambda = ScalarPrd/norm/norm;
    LinCombCoordPoint(1.0, in[(rightmark-1+nin)%nin], lambda, point[1], &(point[2]));
    SubtractCoordPoint(in[rightmark],point[2],&(point[0]));
    EuklidNormCoordPoint(point[0],&norm);
    in[rightmark] = point[2];
    orientation = COUNTERCLOCKWISE;
  }
  else
    return(1);


  /* the main loop */
  FillupStart = -1;
  flag = 0;
  for(i=0; i<nin; i++)
  {
    ClipLine(in[i],in[(i+1)%(nin)],&out1,&out2,&reject,&side1,&side2);
    if (!reject)
    {
      if (flag == 0)
      {
        FirstSide = side1;
        flag = 1;
      }
      if (FillupStart != -1)
        FillUp (FillupStart,side1,orientation,out,nout);
      FillupStart = side2;
      out[(*nout)++] = out1;
      if (side2 != -1)
        out[(*nout)++] = out2;
    }
  }
  if (flag && (FirstSide!=-1))
    FillUp (FillupStart,FirstSide,orientation,out,nout);
  if (!flag)
  {
    /* test if zero vector is in the interior of the polygon */
    if (orientation==1)
    {
      left  = 0;
      for (i=0; i<nin; i++)
        left += (in[(i+1)%(nin)].x*in[i].y >= in[i].x*in[(i+1)%(nin)].y);
      if (left == nin)
        FillUpTot(out,nout);
    }
    else
    {
      right  = 0;
      for (i=0; i<nin; i++)
        right += (in[(i+1)%(nin)].x*in[i].y <= in[i].x*in[(i+1)%(nin)].y);
      if (right == nin)
        FillUpTot(out,nout);
    }
  }
  return (0);
}

/****************************************************************************/
/*D
   UgMove - Move cursor to screen point 'in'

   SYNOPSIS:
   void UgMove (COORD_POINT in)

   PARAMETERS:
   .  in -

   DESCRIPTION:
   This function moves cursor to screen point 'in'.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgMove (COORD_POINT in)
{
        #ifdef ModelP
  if (me != master)
    return;
        #endif

  CurrCursor      = in;
}

/****************************************************************************/
/*D
   UgDraw -  draw line from current cursor position to point

   SYNOPSIS:
   void UgDraw (COORD_POINT point);

   PARAMETERS:
   .  point -

   DESCRIPTION:
   This function draws line from current cursor position to point
   and puts cursor to the end point (not clipping end point).

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgDraw (COORD_POINT point)
{
  SHORT_POINT out1,out2;
  INT reject,dummy;

        #ifdef ModelP
  if (me != master)
    return;
        #endif

  if (ClipLine (CurrCursor,point,&out1,&out2,&reject,&dummy,&dummy)) return;
  if (!reject)
  {
    (*CurrentOutputDevice->Move)(out1);
    (*CurrentOutputDevice->Draw)(out2);
  }
  CurrCursor = point;
}

/****************************************************************************/
/*D
   UgLine - Draw line from point1 to point2

   SYNOPSIS:
   void UgLine (COORD_POINT point1, COORD_POINT point2);

   PARAMETERS:
   .  point1 -
   .  point2 -

   DESCRIPTION:
   This function draws line from point1 to point2.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgLine (COORD_POINT point1, COORD_POINT point2)
{
  SHORT_POINT out1,out2;
  INT reject,dummy;

        #ifdef ModelP
  if (me != master)
    return;
        #endif

  if (ClipLine (point1,point2,&out1,&out2,&reject,&dummy,&dummy)) return;
  if (!reject)
  {
    (*CurrentOutputDevice->Move)(out1);
    (*CurrentOutputDevice->Draw)(out2);
  }
}

/****************************************************************************/
/*D
   UgStyledLine - Draw line from point1 to point2 with the specified pattern

   SYNOPSIS:
   void UgStyledLine (COORD_POINT point1, COORD_POINT point2, DOUBLE dash_length, DOUBLE space_length );

   PARAMETERS:
   .  point1 -
   .  point2 -
   .  dash_length - length of the small line segments in pixel coordinates
   .  space_length - length of the gap betwenn small line segments in pixel coordinates

   DESCRIPTION:
   Along the line from point1 to point2 there will be drawn a line segment
   of length 'dash_length' followed by a gap of length 'space_length' then
   again a linesegment and so on. The lengths are adjusted a little bit thus
   the whole line ends again with a dash.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgStyledLine (COORD_POINT point1, COORD_POINT point2, DOUBLE dash_length, DOUBLE space_length )
{
  SHORT_POINT out1,out2, end;
  INT reject,dummy;
  register double x1, y1, x2, y2;
  double dx_dash, dy_dash, dx_space, dy_space;
  COORD_POINT temp;

        #ifdef ModelP
  if (me != master)
    return;
        #endif

  /* to avoid rounding errors if the same line is drawn twice, one time p1->p2
     and the other time p2->p1, order the points in a general way */
  if ( point1.x > point2.x )
  {
    temp = point1;
    point1 = point2;
    point2 =  temp;
  }

  if (ClipLine (point1,point2,&out1,&out2,&reject,&dummy,&dummy)) return;
  if (reject)
    return;

  /* adjust the dash- and spacelength to fit exactly into the line;
     the line will start and end with a dash */

  /* misuse of x1 and y1 as the slopes and x2 as linelength; y2 as temp */
  x1 = (double)(out2.x - out1.x); y1 = (double)(out2.y - out1.y);
  x2 = sqrt( x1*x1 + y1 * y1 );
  if( fabs(x2) < 1e-20 )
  {             /* linelength very small */
    (*CurrentOutputDevice->Move)(out1);
    (*CurrentOutputDevice->Draw)(out2);
    return;
  }

  /* adjust the dash- and spacelength to fit exactly into the line;
     the line will start and end with a dash! */
  dummy = x2 / (dash_length + space_length) + 0.5;              /* how many pairs of dash+space fit completely into the line? */
  y2 = x2 / ( (dummy+1)*dash_length + dummy*space_length );             /* scaling factor */
  dash_length *= y2;
  space_length *= y2;

  /* increments for dashes and spaces */
  dx_dash = x1 * dash_length / x2; dy_dash = y1 * dash_length / x2;
  dx_space = x1 * space_length / x2; dy_space = y1 * space_length / x2;

  x1 = out1.x; y1 = out1.y;
  end = out2;                   /* save end point */
  out2 = out1;          /* reset */

  while ( out2.x != end.x || out2.y != end.y )
  {
    x2 = x1 + dx_dash; y2 = y1 + dy_dash;
    out2.x = (short)(x2 + 0.5 ); out2.y = (short)(y2 + 0.5 );             /* rounding */

    (*CurrentOutputDevice->Move)(out1);
    (*CurrentOutputDevice->Draw)(out2);

    x1 = x2 + dx_space; y1 = y2 + dy_space;
    out1.x = (short)(x1 + 0.5 ); out1.y = (short)(y1 + 0.5 );             /* rounding */
  }
}

/****************************************************************************/
/*D
   UgInverseLine - Draw line from point1 to point2

   SYNOPSIS:
   void UgInverseLine (COORD_POINT point1, COORD_POINT point2);

   PARAMETERS:
   .  point1 -
   .  point2 -

   DESCRIPTION:
   This function draws line from point1 to point2.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgInverseLine (COORD_POINT point1, COORD_POINT point2)
{
  SHORT_POINT out[2];
  INT reject,dummy;

        #ifdef ModelP
  if (me != master)
    return;
        #endif

  if (ClipLine (point1,point2,&(out[0]),&(out[1]),&reject,&dummy,&dummy)) return;
  if (!reject) (*CurrentOutputDevice->InversePolyline)(out,2);
}

/****************************************************************************/
/*D
   UgPolyLine - Draw polyline from points[0] to points[n-1]

   SYNOPSIS:
   void UgPolyLine (COORD_POINT *points, INT n);

   PARAMETERS:
   .  points - list of points
   .  n -  nb.of points

   DESCRIPTION:
   This function draws polyline from points[0] to points[n-1].

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgPolyLine (COORD_POINT *points, INT n)
{
  SHORT_POINT out1,out2;
  INT reject,dummy,k;

        #ifdef ModelP
  if (me != master)
    return;
        #endif

  for (k=1; k<n; k++)
    if (ClipLine ( points[k-1],points[k],&out1,&out2,&reject,&dummy,&dummy))
      return;
    else if (!reject)
    {
      (*CurrentOutputDevice->Move)(out1);
      (*CurrentOutputDevice->Draw)(out2);
    }
}

/****************************************************************************/
/*D
   UgPolygon - Draw polygon with edge points[0] to points[n-1]

   SYNOPSIS:
   void UgPolygon (COORD_POINT *points, INT n);

   PARAMETERS:
   .  points - list of edge points
   .  n - nb.of points

   DESCRIPTION:
   This function draws polygon with edge points[0] to points[n-1].

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgPolygon (COORD_POINT *points, INT n)
{
  INT nout;
  SHORT_POINT out[MAX_POINTS_OF_POLY];

        #ifdef ModelP
  if (me != master)
    return;
        #endif

  if (ClipPolygon(points, n, out, &nout)) return;
  if (nout<2) return;
  (*CurrentOutputDevice->Polygon)(out, nout);
}
/****************************************************************************/
/*D
   UgInversePolygon - Invert polygon with edge points[0] to points[n-1]

   SYNOPSIS:
   void UgInversePolygon (COORD_POINT *points, INT n);

   PARAMETERS:
   .  points - list of edge points
   .  n - nb.of points

   DESCRIPTION:
   This function inverts polygon with edge points[0] to points[n-1].

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgInversePolygon (COORD_POINT *points, INT n)
{
  INT nout;
  SHORT_POINT out[MAX_POINTS_OF_POLY];

        #ifdef ModelP
  if (me != master)
    return;
        #endif

  if (ClipPolygon(points, n, out, &nout)) return;
  if (nout<2) return;
  (*CurrentOutputDevice->InversePolygon)(out, nout);
}

/****************************************************************************/
/*D
   UgErasePolygon - Draw polygon with edge points[0] to points[n-1]

   SYNOPSIS:
   void UgErasePolygon (COORD_POINT *points, INT n);

   PARAMETERS:
   .  points - list of edge points
   .  n -	nb.of points

   DESCRIPTION:
   This function draws polygon with edge points[0] to points[n-1].

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgErasePolygon (COORD_POINT *points, INT n)
{
  INT nout;
  SHORT_POINT out[MAX_POINTS_OF_POLY];

        #ifdef ModelP
  if (me != master)
    return;
        #endif

  if (ClipPolygon(points, n, out, &nout)) return;
  if (nout<2) return;
  (*CurrentOutputDevice->ErasePolygon)(out, nout);
  return;
}

/****************************************************************************/
/*D
   UgPolymark - Draw n marks

   SYNOPSIS:
   void UgPolymark (COORD_POINT *points, INT n);

   PARAMETERS:
   .  points - nb. of marks
   .  n - nb. of marks

   DESCRIPTION:
   This function draw n marks.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgPolymark (COORD_POINT *points, INT n)
{
  INT k,reject;
  SHORT_POINT out;

        #ifdef ModelP
  if (me != master)
    return;
        #endif

  for (k=0; k<n; k++)
  {
    ClipPoint (points[k],&out,&reject);
    if (!reject)
      (*CurrentOutputDevice->Polymark)(1,&out);
  }
}
/****************************************************************************/
/*D
   UgText -  Draw text s

   SYNOPSIS:
   void UgText (const char *s, INT mode);

   PARAMETERS:
   .  s -
   .  mode -

   DESCRIPTION:
   This function draws text s.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgText (const char *s, INT mode)
{
  INT reject;
  SHORT_POINT out;
  char *p,*next,*move;
  short baseline;
  short TextSize;

        #ifdef ModelP
  if (me != master)
    return;
        #endif

  ClipPoint(CurrCursor,&out,&reject);
  if (reject)
    return;

  switch (mode)
  {
  case TEXT_INDEXED :
    strcpy(buffer,s);
    p = buffer;
    baseline = out.y;
    TextSize = CurrTextSize;

    next = strchr(p,INDEXCHAR);
    if (next!=NULL)
      *next = '\0';
    move = strchr(p,MOVECHAR);
    if (move!=NULL)
      *move = '\0';

    (*CurrentOutputDevice->Move)(out);
    (*CurrentOutputDevice->Text)(p,TEXT_REGULAR);

    while (next!=NULL)
    {
      /* shift cursor to end of last string */
      if (move!=NULL)
        out.x += REL_CHARWIDTH * CurrTextSize * strlen(p) * CurrentOutputDevice->signx;

      p = next;
      next = strchr(++p,INDEXCHAR);
      if (next!=NULL)
        *next = '\0';
      move = strchr(p,MOVECHAR);
      if (move!=NULL)
        *move = '\0';

      switch (*p)
      {
      case 'N' :
        UgSetTextSize(TextSize);
        out.y  = baseline;
        break;
      case 'H' :
        UgSetTextSize(REL_INDEXSIZE*TextSize);
        out.y  = baseline + 0.5*TextSize*CurrentOutputDevice->signy;
        break;
      case 'T' :
        UgSetTextSize(REL_INDEXSIZE*TextSize);
        out.y  = baseline - 0.5*TextSize*CurrentOutputDevice->signy;
        break;
      }
      p++;
      (*CurrentOutputDevice->Move)(out);
      (*CurrentOutputDevice->Text)(p,TEXT_REGULAR);
    }
    return;

  case TEXT_REGULAR :
  case TEXT_INVERSE :
    (*CurrentOutputDevice->Move)(out);
    (*CurrentOutputDevice->Text)(s,mode);
  }
}

/****************************************************************************/
/*D
   UgCenteredText -  Draw text s

   SYNOPSIS:
   void UgCenteredText (COORD_POINT point, const char *s, INT mode);

   PARAMETERS:
   .  point -
   .  s -
   .  mode -

   DESCRIPTION:
   This function draws text s.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgCenteredText (COORD_POINT point, const char *s, INT mode)
{
  INT reject;
  SHORT_POINT out;

        #ifdef ModelP
  if (me != master)
    return;
        #endif

  ClipPoint(point,&out,&reject);
  if (!reject)
  {
    (*CurrentOutputDevice->CenteredText)(out,s,mode);
  }
}

/****************************************************************************/
/*D
   UgSetColor - Set color to colorIndex

   SYNOPSIS:
   void UgSetColor (long colorIndex);

   PARAMETERS:
   .  colorIndex -

   DESCRIPTION:
   This function  sets color to	colorIndex.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgSetColor (long colorIndex)
{
        #ifdef ModelP
  if (me != master)
    return;
        #endif

  (*CurrentOutputDevice->SetColor)(colorIndex);
}

/****************************************************************************/
/*D
   UgSetMarker - Set marker to Index

   SYNOPSIS:
   void UgSetMarker (short index);

   PARAMETERS:
   .  index -

   DESCRIPTION:
   This function sets marker to Index.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgSetMarker (short index)
{
        #ifdef ModelP
  if (me != master)
    return;
        #endif

  (*CurrentOutputDevice->SetMarker)(index);
}

/****************************************************************************/
/*D
   UgSetMarkerSize -  Set marker size to Index

   SYNOPSIS:
   void UgSetMarkerSize (short Index);

   PARAMETERS:
   .  Index -

   DESCRIPTION:
   This function sets marker size to Index.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgSetMarkerSize (short Index)
{
        #ifdef ModelP
  if (me != master)
    return;
        #endif

  (*CurrentOutputDevice->SetMarkerSize)(Index);
}

/****************************************************************************/
/*D
   UgSetTextSize - Set text size to size

   SYNOPSIS:
   void UgSetTextSize (short size)

   PARAMETERS:
   .  size

   DESCRIPTION:
   This function sets text size to size.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgSetTextSize (short size)
{
        #ifdef ModelP
  if (me != master)
    return;
        #endif

  CurrTextSize = size;
  (*CurrentOutputDevice->SetTextSize)(size);
}

/****************************************************************************/
/*D
   UgSetLineWidth - Set line width to width

   SYNOPSIS:
   void UgSetLineWidth (short width);

   PARAMETERS:
   .  width -

   DESCRIPTION:
   This function sets line width to width.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgSetLineWidth (short width)
{
        #ifdef ModelP
  if (me != master)
    return;
        #endif
  (*CurrentOutputDevice->SetLineWidth)(width);
}

/****************************************************************************/
/*D
   UgClearViewPort -  Set line width to width

   SYNOPSIS:
   void UgClearViewPort (void)

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function sets line width to width.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgClearViewPort (void)
{
        #ifdef ModelP
  if (me != master)
    return;
        #endif

  (*CurrentOutputDevice->ClearViewPort)();
}

/****************************************************************************/
/*D
   UgFlush -  Flush the machines graphics buffer

   SYNOPSIS:
   void UgFlush (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function flushes the machines graphics buffer.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgFlush (void)
{
        #ifdef ModelP
  if (me != master)
    return;
        #endif

  (*CurrentOutputDevice->Flush)();
}

/****************************************************************************/
/*D
   UgWait - wait for a time specified in (parts of) seconds

   SYNOPSIS:
   INT UgWait (DOUBLE wait)

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function waits for a time specified in (parts of) seconds.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

void UgWait (DOUBLE wait)
{
  time_t end,time,delta;

  delta = wait*CLOCKS_PER_SEC;
  end   = clock() + delta;
  while ((time=clock())<end)
    if ((end>2*delta) && (time<delta))
      break;                                    /* after wrap around */

  return;
}
