// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  evm.c                                                                                                                 */
/*																			*/
/* Purpose:   elementary vector manipulations								*/
/*																			*/
/* Author:	  Klaus Johannsen												*/
/*			  Institut fuer Computeranwendungen                                                     */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: ug@ica3.uni-stuttgart.de						    */
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

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stddef.h>

#include "compiler.h"
#include "misc.h"
#include "evm.h"
#include "general.h"
#include "ugdevices.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define OneSixth 0.166666666666666667

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

#ifdef __TWODIM__
const DOUBLE unit_vec[DIM][DIM]={{1,0},{0,1}};
#endif

#ifdef __THREEDIM__
const DOUBLE unit_vec[DIM][DIM]={{1,0,0},{0,1,0},{0,0,1}};
#endif

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of macros											*/
/*																			*/
/****************************************************************************/

#define MIN_DETERMINANT                                 0.0001*SMALL_C

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****																	 ****/
/****		general routines											 ****/
/****																	 ****/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
/*D
   ClipRectangleAgainstRectangle - Clip a rectangle against a rectangle

   SYNOPSIS:
   INT ClipRectangleAgainstRectangle (const DOUBLE *r1min, const DOUBLE *r1max,
   DOUBLE *r2min, DOUBLE *r2max);

   PARAMETERS:
   .  r1min - lower left corner of rectangle 1
   .  r1max - upper right corner of rectangle 1
   .  r2min - lower left corner of rectangle 2
   .  r2max - upper right corner of rectangle 2

   DESCRIPTION:
   This function clips the rectangle given by r2 against the rectangle r1,
   i.e. r2 is modified such that it is completely inside of r1.


   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if window is collapsed.
   D*/
/****************************************************************************/

INT ClipRectangleAgainstRectangle (const DOUBLE *r1min, const DOUBLE *r1max, DOUBLE *r2min, DOUBLE *r2max)
{
  if (r2min[0] < r1min[0]) r2min[0] = r1min[0];
  if (r2min[1] < r1min[1]) r2min[1] = r1min[1];
  if (r2max[0] > r1max[0]) r2max[0] = r1max[0];
  if (r2max[1] > r1max[1]) r2max[1] = r1max[1];

  if (r2min[0] >= r2max[0] || r2min[1] >= r2max[1])
    return(1);
  return(0);
}

/****************************************************************************/
/*D
   CheckRectagleIntersection - Check if two rectangles intersect

   SYNOPSIS:
   INT CheckRectagleIntersection (const DOUBLE *r1min, const DOUBLE *r1max,
   const DOUBLE *r2min, const DOUBLE *r2max);

   PARAMETERS:
   .  r1min - lower left corner of rectangle 1
   .  r1max - upper right corner of rectangle 1
   .  r2min - lower left corner of rectangle 2
   .  r2max - upper right corner of rectangle 2

   DESCRIPTION:
   This function clips a rectangle against a rectangle.

   RETURN VALUE:
   INT
   .n    0 if no intersection
   .n    1 if intersection.
   D*/
/****************************************************************************/

INT CheckRectagleIntersection (const DOUBLE *r1min, const DOUBLE *r1max, const DOUBLE *r2min, const DOUBLE *r2max)
{
  if (r1max[0] < r2min[0]) return(0);
  if (r2max[0] < r1min[0]) return(0);
  if (r1max[1] < r2min[1]) return(0);
  if (r2max[1] < r1min[1]) return(0);

  return(1);
}

/****************************************************************************/
/*D
   CheckRectangle - Check if rectangle has a minimum size

   SYNOPSIS:
   INT CheckRectangle (const DOUBLE *rmin, const DOUBLE *rmax,
   const DOUBLE minsize);

   PARAMETERS:
   .  rmin - lower left corner
   .  rmax - upper right corner
   .  minsize - minimal size of rect in x and y direction

   DESCRIPTION:
   This function checks if a rectangle has at least size 'minsize' in `x` and
   `y`-direction.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if rectangle is smaller than 'minsize' in some direction
   D*/
/****************************************************************************/

INT CheckRectangle (const DOUBLE *rmin, const DOUBLE *rmax, const DOUBLE minsize)
{
  if (rmax[0] <= rmin[0]+minsize) return(1);
  if (rmax[1] <= rmin[1]+minsize) return(1);
  return(0);
}

/****************************************************************************/
/*D
   PointInTriangle - Decide if Point lies in the triangle of Points

   SYNOPSIS:
   INT PointInTriangle (const COORD_POINT *Points,
   const COORD_POINT Point);

   PARAMETERS:
   .  Points - Array of three 'COORD_POINT' structures
   .  Point - Point to check

   STRUCTURES:

   .vb
   struct coord_point
   {
    DOUBLE x;
    DOUBLE y;
   };

   typedef struct coord_point COORD_POINT;
   .ve

   DESCRIPTION:
   This function decides if 'Point' lies in the triangle given by 'Points'.

   RETURN VALUE:
   INT
   .n   0 if it lies not in the triangle
   .n   1 if it lies in the triangle.
   D*/
/****************************************************************************/

INT PointInTriangle (const COORD_POINT *Points, const COORD_POINT Point)
{
  DOUBLE M[9], Inverse[9], rhs[3], lambda[3];

  /* invert a 3x3 system */
  M[0]=Points[0].x, M[3]=Points[1].x, M[6] =Points[2].x;
  M[1]=Points[0].y, M[4]=Points[1].y, M[7] =Points[2].y;
  M[2]=1.0,                 M[5]=1.0,             M[8] =1.0;

  if (M3_Invert(Inverse, M)) return (0);
  rhs[0] = Point.x;
  rhs[1] = Point.y;
  rhs[2] = 1.0;
  M3_TIMES_M3(Inverse,rhs,lambda);

  /* decide if Point lies in the interior of Points */
  if (lambda[0]>=0.0 && lambda[1]>=0.0 && lambda[2]>=0.0)
    return (1);
  return (0);
}

/****************************************************************************/
/*D
   PointInPolygon - Decide if Point lies in the polygon of Points

   SYNOPSIS:
   INT PointInPolygon (const COORD_POINT *Points, INT n,
   const COORD_POINT Point);

   PARAMETERS:
   .  Points - polygon given by array of 'COORD_POINT' structures
   .  n - number of corners
   .  Point - Point in question

   STRUCTURES:

   .vb
   struct coord_point
   {
    DOUBLE x;
    DOUBLE y;
   };
   .ve

   DESCRIPTION:
   This function decides if 'Point' lies in the polygon of 'Points'.

   The number of corners of the polygon must be less than or equal
   than 4 in the current implementation!

   RETURN VALUE:
   INT
   .n     0 when lies not in the polygon
   .n     1 when lies in the polygon.
   D*/
/****************************************************************************/

#define POLYMAX         8

INT PointInPolygon (const COORD_POINT *Points, INT n, COORD_POINT Point)
{
  DOUBLE D[POLYMAX] ,tau[POLYMAX],xa,ya,xe,ye;
  int i, left, right;

  assert (n<=POLYMAX);
  if (n<=2) return (0);

  xa = Points[0].x;
  ya = Points[0].y;
  for (i=1; i<=n; i++)
  {
    xe = Points[i%n].x;
    ye = Points[i%n].y;
    D[i-1] = (xe-xa)*(xe-xa)+(ye-ya)*(ye-ya);
    tau[i-1] = (-(ye-ya)*(Point.x-xa)+(xe-xa)*(Point.y-ya));
    xa = xe;
    ya = ye;
  }
  left = right = 0;
  for (i=0; i<n; i++)
  {
    if (tau[i]>=0.0) left++;
    if (tau[i]<=0.0) right++;
    /*	if (tau[i]>=D[i]*SMALL_C) left++;
            if (tau[i]<=-D[i]*SMALL_C) right++;		*/
  }
  if (left==n || right==n)
    return(1);
  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:  PointInPolygonC												*/
/*																			*/
/* Purpose:   decide if Point lies in the polygon of Points with DOUBLE-desc	*/
/*																			*/
/* input:	  const DOUBLE_VECTOR_2D *Points: polygon						*/
/*			  INT n: number of corners										*/
/*			  const DOUBLE_VECTOR_2D Point									*/
/*																			*/
/* return:	  INT 0: lies not in the polygon								*/
/*				  1: lies in the polygon									*/
/*																			*/
/****************************************************************************/

INT PointInPolygonC (const DOUBLE_VECTOR_2D *Points, INT n, const DOUBLE_VECTOR_2D Point)
{
  DOUBLE tau[POLYMAX],xa,ya,xe,ye;
  int i, left, right;

  assert (n<=POLYMAX);
  if (n<=2) return (0);

  xa = Points[0][0];
  ya = Points[0][1];
  for (i=1; i<=n; i++)
  {
    xe = Points[i%n][0];
    ye = Points[i%n][1];
    tau[i-1] = (-(ye-ya)*(Point[0]-xa)+(xe-xa)*(Point[1]-ya));
    xa = xe;
    ya = ye;
  }
  left = right = 0;
  for (i=0; i<n; i++)
  {
    if (tau[i]>=0.0) left++;
    if (tau[i]<=0.0) right++;
  }
  if (left==n || right==n)
    return(1);
  return(0);
}

/****************************************************************************/
/*																			*/
/* Function:  PolyArea														*/
/*																			*/
/* Purpose:   determine area of polygon								                */
/*																			*/
/* input:	  INT n: nb of corners of polygon								*/
/*			  DOUBLE_VECTOR_2D *Polygon: polygon								*/
/*																			*/
/* output:	  DOUBLE *Area: area												*/
/*																			*/
/* return:	  INT 0: ok														*/
/*				  1: error													*/
/*																			*/
/****************************************************************************/

INT PolyArea (INT n, DOUBLE_VECTOR_2D *Polygon, DOUBLE *Area)
{
  INT i;
  DOUBLE c;
  DOUBLE_VECTOR_2D a, b;


  *Area = 0.0;
  if (n<3) return (0);
  for (i=1; i<n-1; i++)
  {
    V2_SUBTRACT(Polygon[i],Polygon[0],a)
    V2_SUBTRACT(Polygon[i+1],Polygon[0],b)
    V2_VECTOR_PRODUCT(a,b,c)
      (*Area) += ABS(c);
  }
  (*Area) *= 0.5;

  return (0);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****																	 ****/
/****		2D routines                                                                                              ****/
/****																	 ****/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
/*D
   M2_Invert - Calculate inverse of a 2x2 DOUBLE matrix

   SYNOPSIS:
   INT M2_Invert (DOUBLE *Inverse, const DOUBLE *Matrix);

   PARAMETERS:
   .  Inverse - inverse of matrix
   .  Matrix - matrix

   DESCRIPTION:
   This function  calculates inverse of a 2x2 DOUBLE matrix.
   The entries of the matrices are given in a linear array with the
   following order -

   .vb
 | Matrix[0] Matrix[1] |
 | Matrix[2] Matrix[3] |
   .ve

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if matrix is nearly singular.
   D*/
/****************************************************************************/

INT M2_Invert (DOUBLE *Inverse, const DOUBLE *Matrix)
{
  DOUBLE det;

  det = Matrix[0]*Matrix[3]-Matrix[1]*Matrix[2];
  if (ABS(det)<SMALL_C) return (1);
  Inverse[0] = Matrix[3]/det;
  Inverse[1] = -Matrix[1]/det;
  Inverse[2] = -Matrix[2]/det;
  Inverse[3] = Matrix[0]/det;

  return (0);
}

/****************************************************************************/
/*D
   M3_Invert - Calculate inverse of a 3x3 DOUBLE matrix

   SYNOPSIS:
   INT M3_Invert (DOUBLE *Inverse, const DOUBLE *Matrix);

   PARAMETERS:
   .  Inverse - inverse of matrix
   .  Matrix - matrix

   DESCRIPTION:
   This function calculates inverse of a 3x3 DOUBLE matrix.
   The entries of the matrices are given in a linear array with the
   following order -

   .vb
 | Matrix[0] Matrix[1] Matrix[2]|
 | Matrix[3] Matrix[4] Matrix[5]|
 | Matrix[6] Matrix[7] Matrix[8]|
   .ve

   RETURN VALUE:
   INT
   .n    0 when ok
   .n    1 when matrix is nearly singular.
   D*/
/****************************************************************************/

INT M3_Invert (DOUBLE *Inverse, const DOUBLE *Matrix)
{
  DOUBLE determinant,invdet;
  INT i,i1,i2, j,j1,j2;

  for (i=0; i<3; i++)
  {
    i1 = (i+1)%3;
    i2 = (i+2)%3;
    for ( j=0; j<3; j++)
    {
      j1 = (j+1)%3;
      j2 = (j+2)%3;
      Inverse[j+3*i] = Matrix[i1+3*j1]*Matrix[i2+3*j2] - Matrix[i1+3*j2]*Matrix[i2+3*j1];
    }
  }
  determinant = Inverse[0+3*0]*Matrix[0+3*0] + Inverse[0+3*1]*Matrix[1+3*0] + Inverse[0+3*2]*Matrix[2+3*0];

  /* check the determinant */
  if (fabs(determinant) > MIN_DETERMINANT)
  {
    invdet = 1.0/determinant;
    for (i=0; i<3; i++)
      for (j=0; j<3; j++)
        Inverse[i+3*j] *= invdet;
    return(0);
  }

  return(1);
}

/****************************************************************************/
/*D
   V2_Normalize	- Normalize a 2D vector

   SYNOPSIS:
   INT V2_Normalize (DOUBLE *a);

   PARAMETERS:
   .  a - input 2D vector (a[0],a[1])

   DESCRIPTION:
   This function normalizes the 2D vector a.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured
   .n    2 if vector is nearly 0.
   D*/
/****************************************************************************/

INT V2_Normalize (DOUBLE *a)
{
  DOUBLE norm;

  V2_EUKLIDNORM(a,norm);
  if (norm < SMALL_C) return(2);
  V2_SCALE(1.0/norm,a);

  return(0);
}

/****************************************************************************/
/*D
   V2_Rotate - Rotate vector by angle

   SYNOPSIS:
   INT V2_Rotate (DOUBLE *vector, DOUBLE alpha);

   PARAMETERS:
   .  vector - 2D vector
   .  alpha - angle in radiant

   DESCRIPTION:
   This function rotates vector 'vector' by angle 'alpha'.
   The vector is turned in mathmatical pos. sense.

   RETURN VALUE:
   INT
   .n      0
   D*/
/****************************************************************************/

INT V2_Rotate (DOUBLE *vector, DOUBLE alpha)
{
  DOUBLE help[2];
  DOUBLE calpha, salpha;

  /* rotate vector */
  help[0] = -vector[1]; help[1] = vector[0];
  calpha = (DOUBLE)cos((double)alpha);
  salpha = (DOUBLE)sin((double)alpha);
  V2_LINCOMB(salpha,help,calpha,vector,vector);

  return (0);
}

/****************************************************************************/
/*D
   V2_IntersectLineSegments - compute intersection of two line segments

   SYNOPSIS:
   INT V2_IntersectLineSegments (const DOUBLE_VECTOR a0, const DOUBLE_VECTOR a1, const DOUBLE_VECTOR b0, const DOUBLE_VECTOR b1, DOUBLE *lambda)

   PARAMETERS:
   .  a0		- begin of line segment a
   .  a1		- end   of line segment a
   .  b0		- begin of line segment b
   .  b1		- end   of line segment b
   .  lambda	- result: parameter of intersection on segment a

   DESCRIPTION:
   This function returns positive number if vector 2 is "left" of vector 1, i.e.
   the third component of the vector product of (x1,y1,0) and (x2,y2,0).

   RETURN VALUE:
   INT
   .n    0: segments intersect
   .n    1<<2: segments are (nearly) parallel
   .n    bit 0 and bit 1 say whether intersection does not lie on segment a and b resp.
   D*/
/****************************************************************************/

INT V2_IntersectLineSegments (const DOUBLE_VECTOR a0, const DOUBLE_VECTOR a1, const DOUBLE_VECTOR b0, const DOUBLE_VECTOR b1, DOUBLE *lambda)
{
  DOUBLE_VECTOR ta,tb,coeff,r,M[DIM],MI[DIM];
  DOUBLE det;
  INT res;

  /* we search the cutting point of line a0+c0*ta with b0-c1*tb by solving the system

                                    T
          (ta[0]  ta[1])    (c0)   (b0[0]-a0[0])
          (			 )    (  ) = (			 )
          (tb[0]  tb[1])    (c1)   (b0[1]-a0[1])
   */

  V2_SUBTRACT(a1,a0,ta);                                                                        /* vector from a0 to a1 */
  V2_SUBTRACT(b0,b1,tb);                                                                        /* vector from b1 to b0 */
  V2_COPY(ta,M[0]);                                                                                     /* transposed coefficient matrix for cut of lines */
  V2_COPY(tb,M[1]);

  M2_INVERT(M,MI,det);                                                                          /* inverse */
  if (ABS(det)<SMALL_D)
    /* lines are parallel */
    return (1<<2);

  V2_SUBTRACT(b0,a0,r);                                                                         /* right hand side */
  MT2_TIMES_V2(MI,r,coeff);                                                                     /* solve for coefficients */
  *lambda = coeff[0];                                                                                   /* local param on side */

  res = 0;
  if (!((-SMALL_C<coeff[0]) && (coeff[0]<1.0+SMALL_C)))         /* local param of segment a not in [0,1]? */
    res |= 1<<0;
  if (!((-SMALL_C<coeff[1]) && (coeff[1]<1.0+SMALL_C)))         /* local param of segment b not in [0,1]? */
    res |= 1<<1;

  return (res);
}

/****************************************************************************/
/*D
   vp - Return positive number if vector 2 is "left" of vector 1

   SYNOPSIS:
   DOUBLE vp (const DOUBLE x1, const DOUBLE y1, const DOUBLE x2, const DOUBLE y2);

   PARAMETERS:
   .  x1,y1 - coordinates of a 2D vector
   .  x2,y2 - coordinates of a 2D vector

   DESCRIPTION:
   This function returns positive number if vector 2 is "left" of vector 1, i.e.
   the third component of the vector product of (x1,y1,0) and (x2,y2,0).

   RETURN VALUE:
   DOUBLE
   D*/
/****************************************************************************/

DOUBLE vp (const DOUBLE x1, const DOUBLE y1, const DOUBLE x2, const DOUBLE y2)
{
  DOUBLE l1,l2;

  l1 = sqrt(x1*x1+y1*y1);
  l2 = sqrt(x2*x2+y2*y2);
  if ((l1<SMALL_D)||(l2<SMALL_D))
    return(0.0);
  else
    return((x1*y2-y1*x2)/(l1*l2));
}

/****************************************************************************/
/*D
   tarea - Compute area of a triangle

   SYNOPSIS:
   DOUBLE tarea (DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2);

   PARAMETERS:
   .  x0,y0 - coordinates of first point
   .  x1,y1 - coordinates of second point
   .  x2,y2 - coordinates of third point

   DESCRIPTION:
   This function computes the area of a triangle.

   RETURN VALUE:
   DOUBLE  area
   D*/
/****************************************************************************/

DOUBLE tarea (DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2)
{
  return(0.5*fabs((y1-y0)*(x2-x0)-(x1-x0)*(y2-y0)));
}
/****************************************************************************/
/*D
   qarea - Compute area of a convex quadrilateral

   SYNOPSIS:
   DOUBLE qarea (DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2,
   DOUBLE x3,DOUBLE y3);

   PARAMETERS:
   .  x0,y0 - coordinates of first point
   .  x1,y1 - coordinates of second point
   .  x2,y2 - coordinates of third point
   .  x3,y3 - coordinates of fourth point

   DESCRIPTION:
   This function computes the area of a convex quadrilateral.

   RETURN VALUE:
   DOUBLE area
   D*/
/****************************************************************************/
DOUBLE qarea (DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2,DOUBLE x3,DOUBLE y3)
{
  return( 0.5*fabs( (y3-y1)*(x2-x0)-(x3-x1)*(y2-y0) ) );
}
/****************************************************************************/
/*D
   c_tarea - Compute area of triangle

   SYNOPSIS:
   DOUBLE c_tarea (const DOUBLE *x0, const DOUBLE *x1, const DOUBLE *x2)

   PARAMETERS:
   .  x0 - Array with coordinates of first point
   .  x1 - Array with coordinates of second point
   .  x2 - Array with coordinates of third point

   DESCRIPTION:
   This function computes the area of a triangle.

   RETURN VALUE:
   DOUBLE area
   D*/
/****************************************************************************/
DOUBLE c_tarea (const DOUBLE *x0, const DOUBLE *x1, const DOUBLE *x2)
{
  return(0.5*fabs((x1[_Y_]-x0[_Y_])*(x2[_X_]-x0[_X_])-(x1[_X_]-x0[_X_])*(x2[_Y_]-x0[_Y_])));
}
/****************************************************************************/
/*D
   c_qarea - Compute area of a convex quadrilateral

   SYNOPSIS:
   DOUBLE c_qarea (const DOUBLE *x0, const DOUBLE *x1, const DOUBLE *x2, const DOUBLE *x3);

   PARAMETERS:
   .  x0 - Array with coordinates of first point
   .  x1 - Array with coordinates of second point
   .  x2 - Array with coordinates of third point
   .  x3 - Array with coordinates of fourth point

   DESCRIPTION:
   This function computes the area of a convex quadrilateral.

   RETURN VALUE:
   DOUBLE area
   D*/
/****************************************************************************/
DOUBLE c_qarea (const DOUBLE *x0, const DOUBLE *x1, const DOUBLE *x2, const DOUBLE *x3)
{
  return( 0.5*fabs( (x3[_Y_]-x1[_Y_])*(x2[_X_]-x0[_X_])-(x3[_X_]-x1[_X_])*(x2[_Y_]-x0[_Y_]) ) );
}

/****************************************************************************/
/*D
   ctarea - Compute area of element wrt cylinder metric

   SYNOPSIS:
   DOUBLE ctarea (DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2)

   PARAMETERS:
   .  x0,y0 - coordinates of first point
   .  x1,y1 - coordinates of second point
   .  x2,y2 - coordinates of third point

   DESCRIPTION:
   This function computes area of element wrt cylinder metric.

   RETURN VALUE:
   DOUBLE
   D*/
/****************************************************************************/

DOUBLE ctarea (DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2)
{
  return((y0+y1+y2) * fabs((y1-y0)*(x2-x0)-(x1-x0)*(y2-y0)) / 6);
}
/****************************************************************************/
/*D
   cqarea - Compute area of element wrt cylinder metric

   SYNOPSIS:
   DOUBLE cqarea (DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2,
   DOUBLE x3,DOUBLE y3);

   PARAMETERS:
   .  x0,y0 - coordinates of first point
   .  x1,y1 - coordinates of second point
   .  x2,y2 - coordinates of third point
   .  x3,y3 - coordinates of fourth point

   DESCRIPTION:
   This function computes area of element wrt cylinder metric.

   RETURN VALUE:
   DOUBLE
   D*/
/****************************************************************************/
DOUBLE cqarea (DOUBLE x0,DOUBLE y0,DOUBLE x1,DOUBLE y1,DOUBLE x2,DOUBLE y2,DOUBLE x3,DOUBLE y3)
{
  return(
           ((y0+y1+y2) * fabs((y1-y0)*(x2-x0)-(x1-x0)*(y2-y0)) +
            (y0+y2+y3) * fabs((y2-y0)*(x3-x0)-(x2-x0)*(y3-y0))) / 6
           );
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****																	 ****/
/****		3D routines                                                                                              ****/
/****																	 ****/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
/*D
   V3_Normalize	- Normalize vector a 3D vector

   SYNOPSIS:
   INT V3_Normalize (DOUBLE *a);

   PARAMETERS:
   .  a - 3D vector

   DESCRIPTION:
   This function normalizes vector a.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    2 if vector is nearly 0.
   D*/
/****************************************************************************/

INT V3_Normalize (DOUBLE *a)
{
  DOUBLE norm;

  V3_EUKLIDNORM(a,norm);
  if (norm < SMALL_C) return(2);
  V3_SCALE(1.0/norm,a);

  return(0);
}

/****************************************************************************/
/*D
   V3_NormVectorProduct - Calculate norm of vector product  a x b

   SYNOPSIS:
   INT V3_NormVectorProduct (const DOUBLE *a, const DOUBLE *b, DOUBLE *result);

   PARAMETERS:
   .  a - input vector (a[0],a[1],a[2])
   .  b - input vector (b[0],b[1],b[2])
   .  result - output scalar result[0]

   DESCRIPTION:
   This function calculates norm of vector product a x b.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT V3_NormVectorProduct (const DOUBLE *a, const DOUBLE *b, DOUBLE *result)
{
  DOUBLE VectorPrd[3];

  V3_VECTOR_PRODUCT(a,b,VectorPrd);
  V3_EUKLIDNORM(VectorPrd,*result);

  return(0);
}

/****************************************************************************/
/*D
   V3_Rotate - Rotate vector around axis by a given angle

   SYNOPSIS:
   INT V3_Rotate (DOUBLE *vector, const DOUBLE *axis, DOUBLE alpha);

   PARAMETERS:
   .  vector - vector to rotate
   .  axis - axis around whoch to rotate
   .  alpha - angle in radiant

   DESCRIPTION:
   This function rotates a given vector around an axis by angle (looking from the
   top). The vector is turned in mathmatical pos. sense.

   RETURN VALUE:
   INT
   .n    0 if o.k.
   .n    1 if axis is nearly zero.
   D*/
/****************************************************************************/

INT V3_Rotate (DOUBLE *vector, const DOUBLE *axis, DOUBLE alpha)
{
  DOUBLE RotationAxis[3], help[3];
  DOUBLE scalarprd, calpha, salpha;

  /* normalize axis */
  V3_COPY(axis,RotationAxis);
  if (V3_Normalize(RotationAxis)) return(1);

  /* rotate vector */
  calpha = (DOUBLE)cos((double)alpha);
  salpha = (DOUBLE)sin((double)alpha);
  V3_SCALAR_PRODUCT(RotationAxis,vector,scalarprd);
  V3_VECTOR_PRODUCT(RotationAxis,vector,help);
  V3_LINCOMB(calpha,vector,salpha,help,help);
  V3_LINCOMB(1.0,help,(1.0-calpha)*scalarprd,RotationAxis,vector);

  return (0);
}

/****************************************************************************/
/*D
   V3_Angle - Calculate angle between two vectors

   SYNOPSIS:
   INT V3_Angle (const DOUBLE *a, const DOUBLE *b, DOUBLE *result);

   PARAMETERS:
   .  a - first vector
   .  b - second vector
   .  result - places result here

   DESCRIPTION:
   This function calculates angle between two vectors.

   RETURN VALUE:
   INT
   .n    0 if o.k.
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT V3_Angle (const DOUBLE *a, const DOUBLE *b, DOUBLE *result)
{
  DOUBLE c, sc, n1, n2;

  V3_EUKLIDNORM(a,n1)
  V3_EUKLIDNORM(b,n2)
  c = n1*n2;
  if (ABS(c)<SMALL_C)
  {
    *result = 0.0;
    return (1);
  }
  V3_SCALAR_PRODUCT(a,b,sc)
  c = sc/c;
  if (c>=1.0)
    *result = 0.0;
  else if (c<=-1.0)
    *result = PI;
  else
    *result = (DOUBLE)acos((double)c);

  return (0);
}

/****************************************************************************/
/*D
   V3_Orthogonalize - Orthgonalize a vector w.r.t. to another vector.

   SYNOPSIS:
   INT V3_Orthogonalize (const DOUBLE *a, const DOUBLE *b, DOUBLE *r);

   PARAMETERS:
   .  a - vector to orthogonalize
   .  b - vector where 'a' is orthogonalized to
   .  r - resulting vector

   DESCRIPTION:
   This function orthgonalizes vector 'a' w.r.t. to 'b'.

   RETURN VALUE:
   INT
   .n    0 if o.k.
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT V3_Orthogonalize (const DOUBLE *a, const DOUBLE *b, DOUBLE *r)
{
  DOUBLE normb, scprd;

  V3_EUKLIDNORM(b,normb)
  if (normb < SMALL_C)
    V3_COPY(a,r)
    else
    {
      V3_SCALAR_PRODUCT(a,b,scprd)
      V3_LINCOMB(1.0,a,-scprd/normb/normb,b,r)
    }

  return (0);
}

/****************************************************************************/
/*D
   V3_Project - Project a vector onto another vector.

   SYNOPSIS:
   INT V3_Project (const DOUBLE *a, const DOUBLE *b, DOUBLE *r);

   PARAMETERS:
   .  a - vector to project
   .  b - vector onto project
   .  r - resulting vector

   DESCRIPTION:
   This function projects vector 'a' onto 'b' store in 'r'.

   RETURN VALUE:
   INT
   .n    0 if o.k.
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT V3_Project (const DOUBLE *a, const DOUBLE *b, DOUBLE *r)
{
  DOUBLE normb, scprd;

  normb = b[0]*b[0]+b[1]*b[1]+b[2]*b[2];
  if (normb==0.0)
    return (1);
  else
  {
    V3_SCALAR_PRODUCT(a,b,scprd)
    scprd /= normb;
    V3_COPY(b,r)
    V3_SCALE(scprd,r)
  }

  return (0);
}

/****************************************************************************/
/*D
   M4_Invert - Invert a 4X4 Matrix

   SYNOPSIS:
   INT M4_Invert (DOUBLE *Inverse, const DOUBLE *Matrix);

   PARAMETERS:
   .  Inverse - output of inverted matrix
   .  Matrix - input matrix to be inverted

   DESCRIPTION:
   This function inverts a 4X4 Matrix.
   The entries of the matrices are given in a linear array with the
   following order -

   .vb
 | Matrix[0]  Matrix[1]  Matrix[2]  Matrix[3]|
 | Matrix[4]  Matrix[5]  Matrix[6]  Matrix[7]|
 | Matrix[8]  Matrix[9]  Matrix[10] Matrix[11]|
 | Matrix[12] Matrix[13] Matrix[14] Matrix[15]|
   .ve

   RETURN VALUE:
   INT
   .n  0 when ok
   .n  1 when matrix is singular or an error occurred.
   D*/
/****************************************************************************/

INT M4_Invert (DOUBLE *Inverse, const DOUBLE *Matrix)
{
  DOUBLE d,dinv;
  INT i,i1,i2,i3, j,j1,j2,j3,sign;

  sign = 0;                     /* no matter which value!!! */

  /* determine submatrices */
  for ( i=0; i<4; i++ )
  {
    i1 = (i+1) & 3;
    i2 = (i1+1) & 3;
    i3 = (i2+1) & 3;
    for ( j=0; j<4; j++ )
    {
      j1 = (j+1) & 3;
      j2 = (j1+1) & 3;
      j3 = (j2+1) & 3;
      Inverse[j+4*i] =   Matrix[i1+4*j1] * (  Matrix[i2+4*j2] * Matrix[i3+4*j3]
                                              - Matrix[i2+4*j3] * Matrix[i3+4*j2] )
                       + Matrix[i1+4*j2] * (  Matrix[i2+4*j3] * Matrix[i3+4*j1]
                                              - Matrix[i2+4*j1] * Matrix[i3+4*j3] )
                       + Matrix[i1+4*j3] * (  Matrix[i2+4*j1] * Matrix[i3+4*j2]
                                              - Matrix[i2+4*j2] * Matrix[i3+4*j1] );

      if (sign) Inverse[j+4*i] *= -1.0;
      sign = !sign;
    }
    sign = !sign;
  }

  /* determine determinant */
  d = Inverse[0+4*0] * Matrix[0+4*0] + Inverse[0+4*1] * Matrix[1+4*0]
      + Inverse[0+4*2] * Matrix[2+4*0] + Inverse[0+4*3] * Matrix[3+4*0];

  /* check determinant and determine inverse */
  if (ABS(d) > MIN_DETERMINANT)
  {
    dinv = 1.0/d;
    for ( i=0; i<4; i++ )
      for (j=0; j<4; j++ ) Inverse[i+4*j] *= dinv;
    return(0);
  }
  return(1);
}

/****************************************************************************/
/*D
   QuadraticFittedMin - determines the minimum-position of a function y=f(x)

   SYNOPSIS:
   INT QuadraticFittedMin (DOUBLE *x, DOUBLE *y, INT n, DOUBLE *minx);

   PARAMETERS:
   .  x1,y1 - coordinates of a 2D vector
   .  x2,y2 - coordinates of a 2D vector

   DESCRIPTION:
   This function determines the minimum-position of a function y=f(x)
   fitted by n points

   RETURN VALUE:
   INT  0: o.k.
            1: over or underflow of n
                2: no quadratic funct with minimum can be fitted
   D*/
/****************************************************************************/

INT QuadraticFittedMin (DOUBLE *x, DOUBLE *y, INT n, DOUBLE *minx)
{
  INT i,j,k;
  DOUBLE mat[50][3],rhs[3],qm[9],qmi[9],coeff[3];

  if (n>50 || n<3) return (1);

  for (i=0; i<n; i++)
  {
    mat[i][0] = 1.0;
    mat[i][1] = x[i];
    mat[i][2] = x[i]*x[i];
  }
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      qm[i+3*j] = 0.0;
      for (k=0; k<n; k++)
        qm[i+3*j] += mat[k][i]*mat[k][j];
    }
    rhs[i] = 0.0;
    for (j=0; j<n; j++)
      rhs[i] += mat[j][i]*y[j];
  }
  if (M3_Invert (qmi,qm)) return (2);
  M3_TIMES_V3(qmi,rhs,coeff);
  if (coeff[2]<=0.0) return (2);
  *minx = -0.5*coeff[1]/coeff[2];

  return (0);
}

/* volume computations, orientation is same as in general element definition ! */
DOUBLE V_te (const DOUBLE *x0, const DOUBLE *x1,
             const DOUBLE *x2, const DOUBLE *x3)
{
  DOUBLE_VECTOR a, b, h, n;

  V3_SUBTRACT(x1,x0,a);
  V3_SUBTRACT(x2,x0,b);
  V3_SUBTRACT(x3,x0,h);
  V3_VECTOR_PRODUCT(a,b,n);

  return(OneSixth*V3_SCAL_PROD(n,h));
}

DOUBLE V_py (const DOUBLE *x0, const DOUBLE *x1, const DOUBLE *x2,
             const DOUBLE *x3, const DOUBLE *x4)
{
  DOUBLE_VECTOR a,b,h,n;

  V3_SUBTRACT(x2,x0,a);
  V3_SUBTRACT(x3,x1,b);
  V3_SUBTRACT(x4,x0,h);
  V3_VECTOR_PRODUCT(a,b,n);

  return(OneSixth*V3_SCAL_PROD(n,h));
}

DOUBLE V_pr (const DOUBLE *x0, const DOUBLE *x1, const DOUBLE *x2,
             const DOUBLE *x3, const DOUBLE *x4, const DOUBLE *x5)
{
  DOUBLE_VECTOR a,b,c,d,e,m,n;

  V3_SUBTRACT(x4,x0,a);
  V3_SUBTRACT(x1,x3,b);
  V3_SUBTRACT(x1,x0,c);
  V3_SUBTRACT(x2,x0,d);
  V3_SUBTRACT(x5,x0,e);
  a[0] = x4[0]-x0[0]; a[1] = x4[1]-x0[1]; a[2] = x4[2]-x0[2];
  b[0] = x1[0]-x3[0]; b[1] = x1[1]-x3[1]; b[2] = x1[2]-x3[2];
  c[0] = x1[0]-x0[0]; c[1] = x1[1]-x0[1]; c[2] = x1[2]-x0[2];
  d[0] = x2[0]-x0[0]; d[1] = x2[1]-x0[1]; d[2] = x2[2]-x0[2];
  e[0] = x5[0]-x0[0]; e[1] = x5[1]-x0[1]; e[2] = x5[2]-x0[2];

  V3_VECTOR_PRODUCT(a,b,m);
  V3_VECTOR_PRODUCT(c,d,n);
  V3_ADD(n,m,n);

  return(OneSixth*V3_SCAL_PROD(n,e));
}

DOUBLE V_he (const DOUBLE *x0, const DOUBLE *x1, const DOUBLE *x2, const DOUBLE *x3,
             const DOUBLE *x4, const DOUBLE *x5, const DOUBLE *x6, const DOUBLE *x7)
{
  return(V_pr(x0,x1,x2,x4,x5,x6)+V_pr(x0,x2,x3,x4,x6,x7));
}

DOUBLE GeneralElementVolume (INT tag, DOUBLE *x_co[])
{
  switch (tag)
  {
#               ifdef __TWODIM__
  case TRIANGLE :
    return(c_tarea (x_co[0],x_co[1],x_co[2]));

  case QUADRILATERAL :
    return(c_qarea (x_co[0],x_co[1],x_co[2],x_co[3]));
#               endif

#               ifdef __THREEDIM__
  case TETRAHEDRON :
    return(V_te(x_co[0],x_co[1],x_co[2],x_co[3]));

  case PYRAMID :
    return (V_py(x_co[0],x_co[1],x_co[2],x_co[3],x_co[4]));

  case PRISM :
    return (V_pr(x_co[0],x_co[1],x_co[2],x_co[3],x_co[4],x_co[5]));

  case HEXAHEDRON :
    return(V_he(x_co[0],x_co[1],x_co[2],x_co[3],x_co[4],x_co[5],x_co[6],x_co[7]));
#                       endif

  default :
    PrintErrorMessage('E',"GeneralElementVolume","unknown element");
    return(0.0);
  }
}

DOUBLE ElementVolume (const ELEMENT *elem)
{
  DOUBLE *x_co[MAX_CORNERS_OF_ELEM];
  INT i;

  for (i=0; i<CORNERS_OF_ELEM(elem); i++)
    x_co[i] = CVECT(MYVERTEX(CORNER(elem,i)));

  return (GeneralElementVolume(TAG(elem),x_co));
}

/****************************************************************************/
/*D
   EXDecomposeMatrixFLOAT - LU decompose a band matrix (FLOAT numbers)

   SYNOPSIS:
   INT EXDecomposeMatrixFLOAT (FLOAT *Mat, INT bw, INT n);

   PARAMETERS:
   .  Mat - pointer to FLOAT array containing the bandmatrix
   .  bw - bandwidth
   .  n - number of rows (==columns) of the matrix

   DESCRIPTION:
   This function calculates the Gauss decomposition of the given band matrix;
   the L and U factors are stored instead of the band matrix.

   RETURN VALUE:
   INT  0: o.k.
        1: main diagonal element to small (0.0); can not divide

   SEE ALSO:
   EXDecomposeMatrixDOUBLE, EXCopyMatrixFLOAT, EXApplyLUFLOAT
   D*/
/****************************************************************************/

INT EXDecomposeMatrixFLOAT (FLOAT *Mat, INT bw, INT n)
{
  INT i,j,k;
  FLOAT f,d;

  for (i=0; i<n-1; i++)
  {
    d = EX_MAT(Mat,bw,i,i);
    if (d==0.0) return (1);
    for (j=i+1; j<=MIN(i+bw,n-1); j++)
    {
      f = EX_MAT(Mat,bw,j,i)/d;
      EX_MAT(Mat,bw,j,i) = f;
      for (k=i+1; k<=MIN(i+bw,n-1); k++)
        EX_MAT(Mat,bw,j,k) -= f*EX_MAT(Mat,bw,i,k);
    }
  }
  return (0);
}

/****************************************************************************/
/*D
   EXDecomposeMatrixDOUBLE - LU decompose a band matrix (DOUBLE numbers)

   SYNOPSIS:
   INT EXDecomposeMatrixDOUBLE (DOUBLE *Mat, INT bw, INT n);

   PARAMETERS:
   .  Mat - pointer to DOUBLE array containing the bandmatrix
   .  bw - bandwidth
   .  n - number of rows (==columns) of the matrix

   DESCRIPTION:
   This function calculates the Gauss decomposition of the given band matrix;
   the L and U factors are stored instead of the band matrix.

   RETURN VALUE:
   INT  0: o.k.
        1: main diagonal element to small (0.0); can not divide

   SEE ALSO:
   EXDecomposeMatrixFLOAT, EXCopyMatrixDOUBLE, EXApplyLUDOUBLE
   D*/
/****************************************************************************/

INT EXDecomposeMatrixDOUBLE (DOUBLE *Mat, INT bw, INT n)
{
  INT i,j,k;
  DOUBLE f,d;

  for (i=0; i<n-1; i++)
  {
    d = EX_MAT(Mat,bw,i,i);
    if (d==0.0) REP_ERR_RETURN (1);
    for (j=i+1; j<=MIN(i+bw,n-1); j++)
    {
      f = EX_MAT(Mat,bw,j,i)/d;
      EX_MAT(Mat,bw,j,i) = f;
      for (k=i+1; k<=MIN(i+bw,n-1); k++)
        EX_MAT(Mat,bw,j,k) -= f*EX_MAT(Mat,bw,i,k);
    }
  }
  return (0);
}

/****************************************************************************/
/*D
   EXApplyLUFLOAT - applies a LU decomposed band matrix (FLOAT numbers)

   SYNOPSIS:
   INT EXApplyLUFLOAT (FLOAT *Mat, INT bw, INT n, DOUBLE *Vec);

   PARAMETERS:
   .  Mat - pointer to FLOAT array containing the bandmatrix
   .  bw - bandwidth
   .  n - number of rows (==columns) of the matrix
   .  Vec - pointer to DOUBLE array containing the vector

   DESCRIPTION:
   This function solves for the LU decomposed band matrix 'Mat' the equation
   L*U x = f.
   f is provided in 'Vec' and the result x is returned again in 'Vec'.

   Note: 'MAt' contains FLOATs whereas 'Vec' contains DOUBLEs.

   RETURN VALUE:
   INT  0: o.k.

   SEE ALSO:
   EXApplyLUDOUBLE, EXCopyMatrixFLOAT, EXDecomposeMatrixFLOAT
   D*/
/****************************************************************************/

INT EXApplyLUFLOAT (FLOAT *Mat, INT bw, INT n, DOUBLE *Vec)
{
  INT i,j;

  /* invert lower */
  for (i=1; i<n; i++)
    for (j=MAX(i-bw,0); j<i; j++)
      Vec[i] -= EX_MAT(Mat,bw,i,j)*Vec[j];

  /* invert upper */
  for (i=n-1; i>=0; i--)
  {
    for (j=i+1; j<=MIN(i+bw,n-1); j++)
      Vec[i] -= EX_MAT(Mat,bw,i,j)*Vec[j];
    Vec[i] /= EX_MAT(Mat,bw,i,i);
  }
  return (0);
}

/****************************************************************************/
/*D
   EXApplyLUDOUBLE - applies a LU decomposed band matrix (DOUBLE numbers)

   SYNOPSIS:
   INT EXApplyLUDOUBLE (DOUBLE *Mat, INT bw, INT n, DOUBLE *Vec);

   PARAMETERS:
   .  Mat - pointer to DOUBLE array containing the bandmatrix
   .  bw - bandwidth
   .  n - number of rows (==columns) of the matrix
   .  Vec - pointer to DOUBLE array containing the vector

   DESCRIPTION:
   This function solves for the LU decomposed band matrix 'Mat' the equation
   L*U x = f.
   f is provided in 'Vec' and the result x is returned again in 'Vec'.

   RETURN VALUE:
   INT  0: o.k.

   SEE ALSO:
   EXApplyLUFLOAT, EXCopyMatrixDOUBLE, EXDecomposeMatrixDOUBLE
   D*/
/****************************************************************************/

INT EXApplyLUDOUBLE (DOUBLE *Mat, INT bw, INT n, DOUBLE *Vec)
{
  INT i,j;

  /* invert lower */
  for (i=1; i<n; i++)
    for (j=MAX(i-bw,0); j<i; j++)
      Vec[i] -= EX_MAT(Mat,bw,i,j)*Vec[j];

  /* invert upper */
  for (i=n-1; i>=0; i--)
  {
    for (j=i+1; j<=MIN(i+bw,n-1); j++)
      Vec[i] -= EX_MAT(Mat,bw,i,j)*Vec[j];
    Vec[i] /= EX_MAT(Mat,bw,i,i);
  }
  return (0);
}

/****************************************************************************/
/*D
   LineISTriangle3D -  gives intersection-point of a line with a triangle in 3D

   SYNOPSIS:
   INT LineISTriangle3D (c1, c2, c3, p1, p2, lambda);

   PARAMETERS:
   .  c1,c2,c3 - corners of the triangle
   .  p1,p2 - endpoints of the line
   .  lambda - local variable of intersection-point on line

   DESCRIPTION:
   This function  gives intersection-point of a line with a triangle
   in 3D if existing

   RETURN VALUE:
   INT  0: no intersection
        1: intersection

   SEE ALSO:
   D*/
/****************************************************************************/

INT LineISTriangle3D (const DOUBLE *c1, const DOUBLE *c2, const DOUBLE *c3, const DOUBLE *p1, const DOUBLE *p2, DOUBLE *lamda)
{
  DOUBLE M[9],Inv[9],sol[3],rhs[3];

  M[0]=c1[0]-c3[0];   M[3]=c2[0]-c3[0];   M[6]=p1[0]-p2[0];   rhs[0]=p1[0]-c3[0];
  M[1]=c1[1]-c3[1];   M[4]=c2[1]-c3[1];   M[7]=p1[1]-p2[1];   rhs[1]=p1[1]-c3[1];
  M[2]=c1[2]-c3[2];   M[5]=c2[2]-c3[2];   M[8]=p1[2]-p2[2];   rhs[2]=p1[2]-c3[2];
  if (M3_Invert(Inv,M)) return (0);
  M3_TIMES_V3(Inv,rhs,sol);
  if (sol[0]<0.0 || sol[1]<0.0 || sol[0]+sol[1]>1.0 || sol[2]<0.0 || sol[2]>1.0) return (0);
  *lamda=sol[2];

  return (1);
}
