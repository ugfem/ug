// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  shapes.c														*/
/*																			*/
/* Purpose:   shape functions for triangles and quadrilaterals				*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: bastian@iwr1.iwr.uni-heidelberg.de					*/
/*																			*/
/* History:   08.04.92 begin, ug version 2.0								*/
/*			  20.11.94 moved shapes.c from ug to cd folder					*/
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

#include "compiler.h"
#include "misc.h"
#include "gm.h"
#include "evm.h"
#include "shapes2d.h"


#ifdef __THREEDIM__
#error this source file is for 2D ONLY
#endif

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

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*D
   N - Shape function

   SYNOPSIS:
   DOUBLE N (int n, int i, DOUBLE s, DOUBLE t);

   PARAMETERS:
   .  n - number of sides (3 for triangle, 4 for quadrangle)
   .  i - corner number (corner number [0..n-1])
   .  s,t - local DOUBLEinates

   DESCRIPTION:
   This function finds the linear shape functions Ni(s,t) to approximate the
   solution u in the integration point ip for triangles and quadrilaterals

   .n   uip(s,t) = SUM Ni(s,t)*ui

   where the sum runs over all nodes of the element to which the considered
   ip belongs. The shape function is defined as

   - for all elements who do not have the node i as a corner
   .n   Ni = 0

   - for the elements
   .n   Ni(node i) = 1
   .n   Ni(node k) = 0, if k is not equal i.

   RETURN VALUE:
   DOUBLE
   .n          value
   D*/
/****************************************************************************/

DOUBLE N (int n, int i, DOUBLE s, DOUBLE t)
{
  if (n==3)
  {
    switch (i)
    {
    case 0 : return(1-s-t);
    case 1 : return(s);
    case 2 : return(t);
    }
  }
  else if (n==4)
  {
    switch (i)
    {
    case 0 : return(0.25*(1-s)*(1-t));
    case 1 : return(0.25*(1+s)*(1-t));
    case 2 : return(0.25*(1+s)*(1+t));
    case 3 : return(0.25*(1-s)*(1+t));
    }
  }

  /* ERROR: i<0 || i>=n */
  return (-1.0);
}

/****************************************************************************/
/*D
   dNds	-  Partial derivative of shape function

   SYNOPSIS:
   DOUBLE dNds (int n, int i, DOUBLE s, DOUBLE t);

   PARAMETERS:
   .  n - number of sides (3 for triangle, 4 for quadrangle)
   .  i - corner number [0..n-1]
   .  s - local DOUBLEinates
   .  t - local DOUBLEinates

   DESCRIPTION:
   This function calculates the partial derivative of the shape function Ni(s,t)
   with respect to s.

   RETURN VALUE:
   DOUBLE
   .n          value
   D*/
/****************************************************************************/

DOUBLE dNds (int n, int i, DOUBLE s, DOUBLE t)
{
  if (n==3)
  {
    switch (i)
    {
    case 0 : return(-1);
    case 1 : return(1);
    case 2 : return(0);
    }
  }
  else if (n==4)
  {
    switch (i)
    {
    case 0 : return(-0.25*(1-t));
    case 1 : return(0.25*(1-t));
    case 2 : return(0.25*(1+t));
    case 3 : return(-0.25*(1+t));
    }
  }

  /* ERROR: i<0 || i>=n */
  return (-1.0);
}

/****************************************************************************/
/*D
   dNdt	- Partial derivative of shape function

   SYNOPSIS:
   DOUBLE dNdt (int n, int i, DOUBLE s, DOUBLE t);

   PARAMETERS:
   .  n - number of sides (for triangle, 4 for quadrangle)
   .  i - corner number [0..n-1]
   .  s - local DOUBLEinates
   .  t - local DOUBLEinates

   DESCRIPTION:
   This function calculates the partial derivative of the shape function Ni(s,t)
   with respect to t.

   RETURN VALUE:
   DOUBLE
   .n          value
   D*/
/****************************************************************************/

DOUBLE dNdt (int n, int i, DOUBLE s, DOUBLE t)
{
  if (n==3)
  {
    switch (i)
    {
    case 0 : return(-1);
    case 1 : return(0);
    case 2 : return(1);
    }
  }
  else if (n==4)
  {
    switch (i)
    {
    case 0 : return(-0.25*(1-s));
    case 1 : return(-0.25*(1+s));
    case 2 : return(0.25*(1+s));
    case 3 : return(0.25*(1-s));
    }
  }

  /* ERROR: i<0 || i>=n */
  return (-1.0);
}

/****************************************************************************/
/*D
   Derivatives - Compute partial derivatives of the shape functions and
   determinant of coordinate transformation

   SYNOPSIS:
   int Derivatives (int n, const DOUBLE *px, const DOUBLE *py, DOUBLE ips,
   DOUBLE ipt, DOUBLE *dNdx, DOUBLE *dNdy, DOUBLE *DetJ);

   PARAMETERS:
   .  n - number of sides (3 for triangle, 4 for quadrangle)
   .  px - x,y DOUBLEinates of corners
   .  py - x,y DOUBLEinates of corners
   .  ips - location to compute derivatives in s,t DOUBLE
   .  ipt - location to compute derivatives in s,t DOUBLE
   .  dNdx - output array for derivatives
   .  dNdy - output array for derivatives
   .  DetJ - determinant of coordinate transformation

   DESCRIPTION:
   This function computes the partial derivatives of the shape
   functions with respect to x,y coordinates at a given
   point (ips,ipt) in s,t DOUBLEinates.
   It is used to compute derivatives	of variables
   in the interior of an element at specific locations. Furthermore
   the determinant of coordinate transformation is calculated.

   RETURN VALUE:
   int
   .n    0 if ok
   .n    1 if determinant of coordinate transformation too small.
   D*/
/****************************************************************************/

int Derivatives (int n, const DOUBLE *px, const DOUBLE *py, DOUBLE ips, DOUBLE ipt, DOUBLE *dNdx, DOUBLE *dNdy, DOUBLE *DetJ)
{
  DOUBLE dydt,dyds,dxdt,dxds,detJ;
  int j;

  dydt = 0.0; for (j=0; j<n; j++) dydt += dNdt(n,j,ips,ipt)*py[j];
  dyds = 0.0; for (j=0; j<n; j++) dyds += dNds(n,j,ips,ipt)*py[j];
  dxdt = 0.0; for (j=0; j<n; j++) dxdt += dNdt(n,j,ips,ipt)*px[j];
  dxds = 0.0; for (j=0; j<n; j++) dxds += dNds(n,j,ips,ipt)*px[j];
  detJ = dxds*dydt-dyds*dxdt;
  if (fabs(detJ)<=SMALL_D) return(1);
  for (j=0; j<n; j++)
  {
    dNdx[j] = ( dydt*dNds(n,j,ips,ipt)-dyds*dNdt(n,j,ips,ipt))/detJ;
    dNdy[j] = (-dxdt*dNds(n,j,ips,ipt)+dxds*dNdt(n,j,ips,ipt))/detJ;
  }
  *DetJ = detJ;
  return(0);
}
/****************************************************************************/
/*D
   Gradients - Compute gradients of the shape functions and
   determinant of coordinate transformation

   SYNOPSIS:
   INT Gradients (INT n, const COORD **theCorners, DOUBLE ips, DOUBLE ipt,
   DOUBLE_VECTOR Gradient[MAX_CORNERS_OF_ELEM], DOUBLE *DetJ);

   PARAMETERS:
   .  n - number of sides (3 for triangle, 4 for quadrangle)
   .  theCorners - coordinates of the element corners
   .  ips - location to compute derivatives in s,t DOUBLE
   .  ipt - location to compute derivatives in s,t DOUBLE
   .  Gradient - output array for gradients
   .  DetJ - determinant of coordinate transformation

   DESCRIPTION:
   This function computes gradients of the shape functions with respect to
   x,y Coordinates at a given point (ips,ipt) and the determinant of
   coordinate transformation.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if determinant of coordinate transformation too small.
   D*/
/****************************************************************************/

INT Gradients (INT n, const COORD **theCorners, DOUBLE ips, DOUBLE ipt, DOUBLE_VECTOR Gradient[MAX_CORNERS_OF_ELEM], DOUBLE *DetJ)
{
  DOUBLE dydt,dyds,dxdt,dxds,detJ;
  int j;

  dydt = 0.0; for (j=0; j<n; j++) dydt += dNdt(n,j,ips,ipt)*theCorners[j][1];
  dyds = 0.0; for (j=0; j<n; j++) dyds += dNds(n,j,ips,ipt)*theCorners[j][1];
  dxdt = 0.0; for (j=0; j<n; j++) dxdt += dNdt(n,j,ips,ipt)*theCorners[j][0];
  dxds = 0.0; for (j=0; j<n; j++) dxds += dNds(n,j,ips,ipt)*theCorners[j][0];
  detJ = dxds*dydt-dyds*dxdt;
  if (fabs(detJ)<=SMALL_D) return(1);
  for (j=0; j<n; j++)
  {
    Gradient[j][_X_] = ( dydt*dNds(n,j,ips,ipt)-dyds*dNdt(n,j,ips,ipt))/detJ;
    Gradient[j][_Y_] = (-dxdt*dNds(n,j,ips,ipt)+dxds*dNdt(n,j,ips,ipt))/detJ;
  }
  *DetJ = detJ;
  return(0);
}
/************************************************************/
/*D
   COORD x_dot_y - Computes the scalar product in 2D

   SYNOPSIS:
   static COORD x_dot_y (const COORD x[2], const COORD y[2]);

   PARAMETERS:
   .  x[2] - coordinates of a point
   .  y[2] - coordinates of a point

   DESCRIPTION:
   This function computes the scalar product of two coordinate
   vectors in two dimensions.

   RETURN VALUE:
   COORD (x[0]*y[0]+x[1]*y[1])
   D*/
/***********************************************************/
static COORD x_dot_y (const COORD x[2], const COORD y[2])
{
  return (x[0]*y[0]+x[1]*y[1]);
}

/************************************************************/
/*D
   x_dot_normal_to_y - dot product in 2D

   SYNOPSIS:
   static COORD x_dot_normal_to_y (const COORD x[2], const COORD y[2]);

   PARAMETERS:
   .  x[2] - coordinates of a point
   .  y[2] - coordinates of a point

   DESCRIPTION:
   This function computes the dot product of two coordinate vectors
   in two dimensions.

   RETURN VALUE:
   COORD (x[0]*y[1]-x[1]*y[0])
   D*/
/***********************************************************/
static COORD x_dot_normal_to_y (const COORD x[2], const COORD y[2])
{
  return (x[0]*y[1]-x[1]*y[0]);
}

/************************************************************/
/*D
   set_x_plus_ay - Adds two coordinate vectors

   SYNOPSIS:
   static INT set_x_plus_ay (COORD res[2], const COORD x[2], const COORD s,
   const COORD y[2]);

   PARAMETERS:
   .  res[2] - result of adding two coordinate vectors
   .  x[2] - coordinates of a point
   .  s - scalar for multiplication
   .  y[2] - coordinates of a point

   DESCRIPTION:
   This function computes the sum of two coordinate vectors where
   one vector before is multplied with a scalar.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if operation failed.
   D*/
/***********************************************************/
static INT set_x_plus_ay (COORD res[2], const COORD x[2], const COORD s, const COORD y[2])
{
  res[0] = x[0]+s*y[0];
  res[1] = x[1]+s*y[1];

  return (0);
}

/****************************************************************************/
/*D
   LocalToGlobal2d - Transform local coordinates to global

   SYNOPSIS:
   INT LocalToGlobal2d (INT n, const COORD **Corners, const COORD_VECTOR EvalPoint,
   COORD_VECTOR GlobalCoord);

   PARAMETERS:
   .  n - corner number (3 or 4)
   .  Corners - coordinates of corners
   .  EvalPoint - local coordinates
   .  GlobalCoord - resulting global coordinates

   DESCRIPTION:
   This function computes the shape functions in an evaluated point and transforms
   the local coordinates to global in triangular and quadrilateral elements.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT LocalToGlobal2d (INT n, const COORD **Corners, const COORD_VECTOR EvalPoint, COORD_VECTOR GlobalCoord)
{
  COORD N0,N1,N2,N3;

  if (n==3)
  {
    N0 = N(3,0,EvalPoint[0],EvalPoint[1]);
    N1 = N(3,1,EvalPoint[0],EvalPoint[1]);
    N2 = N(3,2,EvalPoint[0],EvalPoint[1]);
    V2_LINCOMB(N0,Corners[0],N1,Corners[1],GlobalCoord)
    V2_LINCOMB(1.0,GlobalCoord,N2,Corners[2],GlobalCoord)
    return (0);
  }
  if (n==4)
  {
    N0 = N(4,0,EvalPoint[0],EvalPoint[1]);
    N1 = N(4,1,EvalPoint[0],EvalPoint[1]);
    N2 = N(4,2,EvalPoint[0],EvalPoint[1]);
    N3 = N(4,3,EvalPoint[0],EvalPoint[1]);
    V2_LINCOMB(N0,Corners[0],N1,Corners[1],GlobalCoord)
    V2_LINCOMB(1.0,GlobalCoord,N2,Corners[2],GlobalCoord)
    V2_LINCOMB(1.0,GlobalCoord,N3,Corners[3],GlobalCoord)
    return (0);
  }
  return (1);
}

/****************************************************************************/
/*D
   L2GDerivative2d - Derivative of LocalToGlobal2d

   SYNOPSIS:
   INT L2GDerivative2d (INT n, const COORD **Corners,
   const COORD_VECTOR EvalPoint, COORD *Derivative);

   PARAMETERS:
   .  n - corner number
   .  Corners - coordinates of corners
   .  EvalPoint - local coordinates
   .  Derivative - df1/ds, df2/ds, df1/dt, df2/dt

   DESCRIPTION:
   This function calculates the derivates of the shape functions in an
   evaluated point and transforms the local coordinates of the derivates to
   global in triangular and quadrilateral elements.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.

   SEE ALSO:
   LocalToGlobal2d
   D*/
/****************************************************************************/

INT L2GDerivative2d (INT n, const COORD **Corners, const COORD_VECTOR EvalPoint, COORD *Derivative)
{
  COORD dNds0,dNds1,dNds2,dNds3;
  COORD dNdt0,dNdt1,dNdt2,dNdt3;

  if (n==3)
  {
    dNds0 = dNds(3,0,EvalPoint[0],EvalPoint[1]);
    dNds1 = dNds(3,1,EvalPoint[0],EvalPoint[1]);
    dNds2 = dNds(3,2,EvalPoint[0],EvalPoint[1]);
    dNdt0 = dNdt(3,0,EvalPoint[0],EvalPoint[1]);
    dNdt1 = dNdt(3,1,EvalPoint[0],EvalPoint[1]);
    dNdt2 = dNdt(3,2,EvalPoint[0],EvalPoint[1]);
    V2_LINCOMB(dNds0,Corners[0],dNds1,Corners[1],Derivative)
    V2_LINCOMB(1.0,Derivative,dNds2,Corners[2],Derivative)
    V2_LINCOMB(dNdt0,Corners[0],dNdt1,Corners[1],Derivative+2)
    V2_LINCOMB(1.0,Derivative+2,dNdt2,Corners[2],Derivative+2)
    return (0);
  }
  if (n==4)
  {
    dNds0 = dNds(4,0,EvalPoint[0],EvalPoint[1]);
    dNds1 = dNds(4,1,EvalPoint[0],EvalPoint[1]);
    dNds2 = dNds(4,2,EvalPoint[0],EvalPoint[1]);
    dNds3 = dNds(4,3,EvalPoint[0],EvalPoint[1]);
    dNdt0 = dNdt(4,0,EvalPoint[0],EvalPoint[1]);
    dNdt1 = dNdt(4,1,EvalPoint[0],EvalPoint[1]);
    dNdt2 = dNdt(4,2,EvalPoint[0],EvalPoint[1]);
    dNdt3 = dNdt(4,3,EvalPoint[0],EvalPoint[1]);
    V2_LINCOMB(dNds0,Corners[0],dNds1,Corners[1],Derivative)
    V2_LINCOMB(1.0,Derivative,dNds2,Corners[2],Derivative)
    V2_LINCOMB(1.0,Derivative,dNds3,Corners[3],Derivative)
    V2_LINCOMB(dNdt0,Corners[0],dNdt1,Corners[1],Derivative+2)
    V2_LINCOMB(1.0,Derivative+2,dNdt2,Corners[2],Derivative+2)
    V2_LINCOMB(1.0,Derivative+2,dNdt3,Corners[3],Derivative+2)
    return (0);
  }
  return (1);
}

/****************************************************************************/
/*D
   GlobalToLocal2d - Transforms global coordinates to local

   SYNOPSIS:
   INT GlobalToLocal2d (INT n, const COORD **Corners,
   const COORD_VECTOR EvalPoint, COORD_VECTOR LocalCoord);

   PARAMETERS:
   .  n - corner number
   .  Corners - coordinates of the element corners
   .  EvalPoint - global coordinates
   .  LocalCoord - resulting local coordinates

   DESCRIPTION:
   This function transforms global coordinates to local in triangular and
   quadrilateral element. The global point has to be inside the triangle/quadrilateral.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if error occured.
   D*/
/****************************************************************************/

INT GlobalToLocal2d (INT n, const COORD **Corners, const COORD_VECTOR EvalPoint, COORD_VECTOR LocalCoord)
{
  COORD xs[2],xt[2],xst[2],x0[2],aux1[2],aux2[2];
  COORD xtdnxst,xsdnxst,a,b,c,t1,t2;
  COORD t1x,t2x,t3x,t1y,t2y,t3y,D;


  if (n==3)
  {
    /* its a triangle */
    t1x = Corners[1][0]-Corners[0][0];
    t2x = Corners[2][0]-Corners[0][0];
    t3x = Corners[0][0];
    t1y = Corners[1][1]-Corners[0][1];
    t2y = Corners[2][1]-Corners[0][1];
    t3y = Corners[0][1];

    D = t1x*t2y-t2x*t1y;

    if (D<0.0)
      return(1);

    LocalCoord[0] = (t2y*(EvalPoint[0]-t3x)-t2x*(EvalPoint[1]-t3y))/D;
    LocalCoord[1] = (-t1y*(EvalPoint[0]-t3x)+t1x*(EvalPoint[1]-t3y))/D;

    return (0);
  }

  if (n!=4)
    return (2);

  /* coeefficients of the local coordinates s,t in the bilinear form:
     EvalPoint = sum_i N_i(s,t) Corners_i */

  /* vector from center of mass to EvalPoint */
  x0[0]  = EvalPoint[0] - 0.25 * ( Corners[0][0]+Corners[1][0]+Corners[2][0]+Corners[3][0]);
  x0[1]  = EvalPoint[1] - 0.25 * ( Corners[0][1]+Corners[1][1]+Corners[2][1]+Corners[3][1]);

  /* vector from mid of side 3 to side 1 */
  xs[0]  = 0.25 * (-Corners[0][0]+Corners[1][0]+Corners[2][0]-Corners[3][0]);
  xs[1]  = 0.25 * (-Corners[0][1]+Corners[1][1]+Corners[2][1]-Corners[3][1]);

  /* vector from mid of side 2 to side 0 */
  xt[0]  = 0.25 * (-Corners[0][0]-Corners[1][0]+Corners[2][0]+Corners[3][0]);
  xt[1]  = 0.25 * (-Corners[0][1]-Corners[1][1]+Corners[2][1]+Corners[3][1]);

  /* vector from mid of diagonal 1 to diagonal 0 */
  xst[0] = 0.25 * ( Corners[0][0]-Corners[1][0]+Corners[2][0]-Corners[3][0]);
  xst[1] = 0.25 * ( Corners[0][1]-Corners[1][1]+Corners[2][1]-Corners[3][1]);


  /* scalarproducts with normal vectors (if vanishing ==> vectors parallel) */
  xtdnxst = x_dot_normal_to_y(xt,xst);
  xsdnxst = x_dot_normal_to_y(xs,xst);

  /* NB: xs || xt only in the degenerate case */

  /* case 1: xst = (0,0)-vector --- this is the case of a parallelogram */
  if (x_dot_y(xst,xst)<SMALL_C)
  {
    a = x_dot_normal_to_y(xs,xt);

    LocalCoord[0] =  x_dot_normal_to_y(x0,xt) / a;
    LocalCoord[1] = -x_dot_normal_to_y(x0,xs) / a;                      /* xddny = -nxdy */

    return (0);
  }

  /* case 2: xs || xst --- this is the case of a trapezoid, side 0 parallel side 2 */
  if (fabs(xsdnxst)<SMALL_C)
  {
    LocalCoord[1] = x_dot_normal_to_y(x0,xst) / xtdnxst;

    set_x_plus_ay(aux1,x0,-LocalCoord[1],xt);
    set_x_plus_ay(aux2,xs,LocalCoord[1],xst);

    LocalCoord[0] = x_dot_y(aux1,aux2) / x_dot_y(aux2,aux2);

    return (0);
  }

  /* case 3: xt || xst --- this is the case of a trapezoid, side 1 parallel side 3 */
  if (fabs(xtdnxst)<SMALL_C)
  {
    LocalCoord[0] = x_dot_normal_to_y(x0,xst) / xsdnxst;

    set_x_plus_ay(aux1,x0,-LocalCoord[0],xs);
    set_x_plus_ay(aux2,xt,LocalCoord[0],xst);

    LocalCoord[1] = x_dot_y(aux1,aux2) / x_dot_y(aux2,aux2);

    return (0);
  }

  /* the general case */
  a = 0.5*(x_dot_normal_to_y(xt,xs) - x_dot_normal_to_y(x0,xst)) / xtdnxst;
  b = -x_dot_normal_to_y(x0,xs) / xtdnxst;
  c = a*a - b;

  if (c<0.0)
    return (3);
  c = sqrt(c);

  t1 = -a + c;
  t2 = -a - c;

  if ((-1.1<=t1) && (t1<=1.1))
  {
    LocalCoord[1] = t1;

    set_x_plus_ay(aux1,x0,-t1,xt);
    LocalCoord[0] = x_dot_normal_to_y(aux1,xst) / xsdnxst;
  }
  else if ((-1.1<=t2) && (t2<=1.1))
  {
    LocalCoord[1] = t2;

    set_x_plus_ay(aux1,x0,-t2,xt);
    LocalCoord[0] = x_dot_normal_to_y(aux1,xst) / xsdnxst;
  }
  else
    return (4);

  if ((-1.1<=LocalCoord[0]) && (LocalCoord[0]<=1.1))
    return (0);
  else
    return (5);
}

/* in the specialGlobalToLocal2d the global point has not to be inside the triangle/quadrilateral */
/**************************************************************/
/*D
   specialGlobalToLocal2d - Transforms global coordinates to local

   SYNOPSIS:
   INT specialGlobalToLocal2d (INT n, const COORD **Corners,
   const COORD_VECTOR EvalPoint, COORD_VECTOR LocalCoord);

   PARAMETERS:
   .  n - corner number
   .  Corners - coordinates of corners
   .  EvalPoint - global coordinates
   .  LocalCoord - resulting local coordinates

   DESCRIPTION:
   This function transforms global coordinates to local in triangular and
   quadrilateral element. The global point has not to be inside the triangle/quadrilateral.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error occured.
   D*/
/***************************************************************/
INT specialGlobalToLocal2d (INT n, const COORD **Corners, const COORD_VECTOR EvalPoint, COORD_VECTOR LocalCoord)
{
  COORD xs[2],xt[2],xst[2],x0[2],aux1[2],aux2[2];
  COORD xtdnxst,xsdnxst,a,b,c,t1,t2,s1,s2;
  COORD t1x,t2x,t3x,t1y,t2y,t3y,D;


  if (n==3)
  {
    /* its a triangle */
    t1x = Corners[1][0]-Corners[0][0];
    t2x = Corners[2][0]-Corners[0][0];
    t3x = Corners[0][0];
    t1y = Corners[1][1]-Corners[0][1];
    t2y = Corners[2][1]-Corners[0][1];
    t3y = Corners[0][1];

    D = t1x*t2y-t2x*t1y;

    if (D<0.0)
      return(1);

    LocalCoord[0] = (t2y*(EvalPoint[0]-t3x)-t2x*(EvalPoint[1]-t3y))/D;
    LocalCoord[1] = (-t1y*(EvalPoint[0]-t3x)+t1x*(EvalPoint[1]-t3y))/D;

    return (0);
  }

  if (n!=4)
    return (2);

  /* coeefficients of the local coordinates s,t in the bilinear form:
     EvalPoint = sum_i N_i(s,t) Corners_i */

  /* center of mass */
  x0[0]  = EvalPoint[0] - 0.25 * ( Corners[0][0]+Corners[1][0]+Corners[2][0]+Corners[3][0]);
  x0[1]  = EvalPoint[1] - 0.25 * ( Corners[0][1]+Corners[1][1]+Corners[2][1]+Corners[3][1]);

  /* vector from mid of side 3 to side 1 */
  xs[0]  = 0.25 * (-Corners[0][0]+Corners[1][0]+Corners[2][0]-Corners[3][0]);
  xs[1]  = 0.25 * (-Corners[0][1]+Corners[1][1]+Corners[2][1]-Corners[3][1]);

  /* vector from mid of side 2 to side 0 */
  xt[0]  = 0.25 * (-Corners[0][0]-Corners[1][0]+Corners[2][0]+Corners[3][0]);
  xt[1]  = 0.25 * (-Corners[0][1]-Corners[1][1]+Corners[2][1]+Corners[3][1]);

  /* vector from mid of diagonal 1 to diagonal 0 */
  xst[0] = 0.25 * ( Corners[0][0]-Corners[1][0]+Corners[2][0]-Corners[3][0]);
  xst[1] = 0.25 * ( Corners[0][1]-Corners[1][1]+Corners[2][1]-Corners[3][1]);


  /* scalarproducts with normal vectors (if vanishing ==> vectors parallel) */
  xtdnxst = x_dot_normal_to_y(xt,xst);
  xsdnxst = x_dot_normal_to_y(xs,xst);

  /* NB: xs || xt only in the degenerate case */

  /* case 1: xst = (0,0)-vector --- this is the case of a parallelogram */
  if (x_dot_y(xst,xst)<SMALL_C)
  {
    a = x_dot_normal_to_y(xs,xt);
    LocalCoord[0] =  x_dot_normal_to_y(x0,xt) / a;
    LocalCoord[1] = -x_dot_normal_to_y(x0,xs) / a;                      /* xddny = -nxdy */

    return (0);
  }

  /* case 2: xs || xst --- this is the case of a trapezoid, side 0 parallel side 2 */
  if (fabs(xsdnxst)<SMALL_C)
  {
    LocalCoord[1] = x_dot_normal_to_y(x0,xst) / xtdnxst;

    set_x_plus_ay(aux1,x0,-LocalCoord[1],xt);
    set_x_plus_ay(aux2,xs,LocalCoord[1],xst);

    LocalCoord[0] = x_dot_y(aux1,aux2) / x_dot_y(aux2,aux2);

    return (0);
  }

  /* case 3: xt || xst --- this is the case of a trapezoid, side 1 parallel side 3 */
  if (fabs(xtdnxst)<SMALL_C)
  {
    LocalCoord[0] = x_dot_normal_to_y(x0,xst) / xsdnxst;

    set_x_plus_ay(aux1,x0,-LocalCoord[0],xs);
    set_x_plus_ay(aux2,xt,LocalCoord[0],xst);

    LocalCoord[1] = x_dot_y(aux1,aux2) / x_dot_y(aux2,aux2);

    return (0);
  }

  /* the general case */
  a = 0.5*(x_dot_normal_to_y(xt,xs) - x_dot_normal_to_y(x0,xst)) / xtdnxst;
  b = -x_dot_normal_to_y(x0,xs) / xtdnxst;
  c = a*a - b;

  if (c<0.0)
    return (3);
  c = sqrt(c);

  /* first solution */
  t1 = -a + c;
  set_x_plus_ay(aux1,x0,-t1,xt);
  s1 = x_dot_normal_to_y(aux1,xst) / xsdnxst;

  /* second solution */
  t2 = -a - c;
  set_x_plus_ay(aux1,x0,-t2,xt);
  s2 = x_dot_normal_to_y(aux1,xst) / xsdnxst;

  /* take the one next to the origin */
  if ((s1*s1+t1*t1) < (s2*s2+t2*t2))
  {
    LocalCoord[0] = s1;
    LocalCoord[1] = t1;
  }
  else
  {
    LocalCoord[0] = s2;
    LocalCoord[1] = t2;
  }

  return (0);
}
