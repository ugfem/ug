// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/****************************************************************************/
/*																			*/
/* File:	  shapes.c														*/
/*																			*/
/* Purpose:   shape functions for triangles and quadrilaterals				*/
/*																			*/
/* Author:	  Peter Bastian 												*/
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: ug@ica3.uni-stuttgart.de					        */
/*																			*/
/* History:   08.04.92 begin, ug version 2.0								*/
/*			  20.11.94 moved shapes.c from ug to cd folder					*/
/*																			*/
/* Remarks: 																*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files 									*/
/*																			*/
/****************************************************************************/

#include <math.h>
#include <assert.h>

#include "compiler.h"
#include "misc.h"
#include "gm.h"
#include "evm.h"
#include "shapes.h"

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

#ifdef __THREEDIM__
COORD_VECTOR HexRefCoord[MAX_CORNERS_OF_ELEM] = {{0.0,0.0,0.0},{1.0,0.0,0.0},
						 {1.0,1.0,0.0},{0.0,1.0,0.0},
						 {0.0,0.0,1.0},{1.0,0.0,1.0},
						 {1.0,1.0,1.0},{0.0,1.0,1.0}};

/* Coordinates of the central points in reference hexahedron of inner surfaces 
   where integration is done in finite-volume discretization.        
   Used in "AssembleElementHEX_FV".  */ 
COORD_VECTOR CenterOfIntergrSurf[MAX_EDGES_OF_ELEM] = {
  {0.5,0.25,0.25},{0.75,0.5,0.25},{0.5,0.75,0.25},{0.25,0.5,0.25},
  {0.5,0.25,0.75},{0.75,0.5,0.75},{0.5,0.75,0.75},{0.25,0.5,0.75},
  {0.75,0.25,0.5},{0.25,0.25,0.5},{0.75,0.75,0.5},{0.25,0.75,0.5}};	   

COORD_VECTOR TransfCoeff[MAX_CORNERS_OF_ELEM]; 
#endif

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* local midpoints */
static COORD_VECTOR_2D LMP_Triangle 		= {0.3333333333333, 0.3333333333333};
static COORD_VECTOR_2D LMP_Quadrilateral	= {0.5, 0.5};
static COORD_VECTOR_3D LMP_Tetrahedron		= {0.25, 0.25, 0.25};



/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*D
   GN - General Shape function for nodes

   SYNOPSIS:
   DOUBLE GN (INT n, INT i, COORD local);

   PARAMETERS:
.  n - number of corners of the element
.  i - corner number (corner number [0..n-1])
.  local - local COORDinates

   DESCRIPTION:
   This function finds the linear/bilinear shape functions Ni(local) to approximate the
   solution u in the integration point ip for triangles/quadrilaterals/tetrahedra
   
.n   uip(local) = SUM Ni(local)*ui

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

DOUBLE GN (INT n, INT i, COORD *local)
{
#ifdef __TWODIM__
	switch (n)
	{
		case 3:
			switch (i)
			{
				case 0 : return((DOUBLE)(1-local[0]-local[1]));
				case 1 : return((DOUBLE)local[0]);
				case 2 : return((DOUBLE)local[1]);
			}
		case 4:
			switch (i)
			{
				case 0 : return((DOUBLE)(0.25*(1-local[0])*(1-local[1])));
				case 1 : return((DOUBLE)(0.25*(1+local[0])*(1-local[1])));
				case 2 : return((DOUBLE)(0.25*(1+local[0])*(1+local[1])));
				case 3 : return((DOUBLE)(0.25*(1-local[0])*(1+local[1])));
			}
		default:
			return (-1.0);
	}
#endif
	
#ifdef __THREEDIM__
	switch (n)
	{
		case 4:
			switch (i)
			{
				case 0 : return((DOUBLE)(1.0-local[0]-local[1]-local[2]));
				case 1 : return((DOUBLE)local[0]);
				case 2 : return((DOUBLE)local[1]);
				case 3 : return((DOUBLE)local[2]);
			}
		default:
			return (-1.0);
	}
#endif
}

/****************************************************************************/
/*D
   LMP - local midpoint

   SYNOPSIS:
   COORD *LMidP (INT n);

   PARAMETERS:
.  n - number of corners of the element

   DESCRIPTION:
   This function gives the local coordinates of the midpoint of an element
      
   RETURN VALUE:
   COORD *
.n          local
D*/   
/****************************************************************************/

COORD *LMP (INT n)
{
#ifdef __TWODIM__
	switch (n)
	{
		case 3:	return (LMP_Triangle);
		case 4:	return (LMP_Quadrilateral);
	}
#endif
	
#ifdef __THREEDIM__
	switch (n)
	{
		case 4:	return (LMP_Tetrahedron);
	}
#endif
}

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

#ifdef __TWODIM__
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
#endif

#ifdef __THREEDIM__						 
DOUBLE N (const INT i, const COORD *LocalCoord)
{
	switch (i)
	{
		case 0 : return((DOUBLE)(1.0-LocalCoord[0]-LocalCoord[1]-LocalCoord[2]));
		case 1 : return((DOUBLE)LocalCoord[0]);
		case 2 : return((DOUBLE)LocalCoord[1]);
		case 3 : return((DOUBLE)LocalCoord[2]);
	}
	
	/* error */
	assert(FALSE);
	return (-1);
}
#endif

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

#ifdef __TWODIM__
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
#endif

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

#ifdef __TWODIM__
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
#endif

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

#ifdef __TWODIM__
INT Derivatives (INT n, const DOUBLE *px, const DOUBLE *py, DOUBLE ips, DOUBLE ipt, DOUBLE *dNdx, DOUBLE *dNdy, DOUBLE *DetJ)
{
	DOUBLE dydt,dyds,dxdt,dxds,detJ;
	INT j;
	
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
#endif

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

#ifdef __TWODIM__
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
#endif

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

#ifdef __TWODIM__
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
#endif

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

#ifdef __TWODIM__
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
#endif

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

static COORD x_dot_y (const COORD x[2], const COORD y[2])
{
	return (x[0]*y[0]+x[1]*y[1]);
}

static INT set_x_plus_ay (COORD res[2], const COORD x[2], const COORD s, const COORD y[2])
{
	res[0] = x[0]+s*y[0];
	res[1] = x[1]+s*y[1];
	
	return (0);
}

static COORD x_dot_normal_to_y (const COORD x[2], const COORD y[2])
{
	return (x[0]*y[1]-x[1]*y[0]);
}

#ifdef __TWODIM__
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
		LocalCoord[1] = -x_dot_normal_to_y(x0,xs) / a;		/* xddny = -nxdy */
		
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
#endif

/*****************************************************************/
/*                                                               */
/*    NiHex: shape functions in the reference hexahedron         */						
/*                                                               */
/*****************************************************************/

#ifdef __THREEDIM__						
DOUBLE NiHex (int i, COORD *alpha)
{
	switch (i)
	{
		case (0): return((DOUBLE)(-(-1.0+alpha[0])*(-1.0+alpha[1])*(-1.0+alpha[2])));
		case (1): return((DOUBLE)(alpha[0]*(alpha[1]-1.0)*(alpha[2]-1.0)));
		case (2): return((DOUBLE)(-alpha[0]*alpha[1]*(alpha[2]-1.0)));
		case (3): return((DOUBLE)((alpha[0]-1.0)*alpha[1]*(alpha[2]-1.0)));
		case (4): return((DOUBLE)((alpha[0]-1.0)*(alpha[1]-1.0)*alpha[2]));
		case (5): return((DOUBLE)(-alpha[0]*(alpha[1]-1.0)*alpha[2]));
		case (6): return((DOUBLE)(alpha[0]*alpha[1]*alpha[2]));
		case (7): return((DOUBLE)(-(alpha[0]-1.0)*alpha[1]*alpha[2]));
	}
}
#endif

/*****************************************************************/
/*                                                               */
/* NiHexDer: derivatives of the shape functions in ref.hex.      */
/*                                                               */
/*****************************************************************/

#ifdef __THREEDIM__
DOUBLE NiHexDer (int i, int j, COORD *alpha)
{
	switch (i)
	{
		case (0):
			switch (j)
			{
				case (0): return((DOUBLE)(-(alpha[1]-1.0)*(alpha[2]-1.0)));
				case (1): return((DOUBLE)(-(alpha[0]-1.0)*(alpha[2]-1.0)));
				case (2): return((DOUBLE)(-(alpha[0]-1.0)*(alpha[1]-1.0)));
			}
		case (1):
			switch (j)
			{
				case (0): return((DOUBLE)((alpha[1]-1.0)*(alpha[2]-1.0)));
				case (1): return((DOUBLE)(alpha[0]*(alpha[2]-1.0)));
				case (2): return((DOUBLE)(alpha[0]*(alpha[1]-1.0)));
			}
		case (2):
			switch (j)
			{
				case (0): return((DOUBLE)(-alpha[1]*(alpha[2]-1.0)));
				case (1): return((DOUBLE)(-alpha[0]*(alpha[2]-1.0)));
				case (2): return((DOUBLE)(-alpha[0]*alpha[1]));			
			}
		case (3):
			switch (j)
			{
				case (0): return((DOUBLE)(alpha[1]*(alpha[2]-1.0)));
				case (1): return((DOUBLE)((alpha[0]-1.0)*(alpha[2]-1.0)));
				case (2): return((DOUBLE)((alpha[0]-1.0)*alpha[1]));
			}
		case (4):
			switch (j)
			{
				case (0): return((DOUBLE)((alpha[1]-1.0)*alpha[2]));
				case (1): return((DOUBLE)((alpha[0]-1.0)*alpha[2]));
				case (2): return((DOUBLE)((alpha[0]-1.0)*(alpha[1]-1.0)));
			}
		case (5):
			switch (j)
			{
				case (0): return((DOUBLE)(-(alpha[1]-1.0)*alpha[2]));
				case (1): return((DOUBLE)(-alpha[0]*alpha[2]));
				case (2): return((DOUBLE)(-alpha[0]*(alpha[1]-1.0)));
			}
		case (6):
			switch (j)
			{
				case (0): return((DOUBLE)(alpha[1]*alpha[2]));
				case (1): return((DOUBLE)(alpha[0]*alpha[2]));
				case (2): return((DOUBLE)(alpha[0]*alpha[1]));
			}
		case (7):
			switch (j)
			{
				case (0): return((DOUBLE)(-alpha[1]*alpha[2]));
				case (1): return((DOUBLE)(-(alpha[0]-1.0)*alpha[2]));
				case (2): return((DOUBLE)(-(alpha[0]-1.0)*alpha[1]));
			}
	}
}
#endif
						
/*****************************************************************************/
/*                                                                           */
/* TransformCoefficients: computes coefficients of trasformation             */
/*                        of theElement to unit cube :                       */
/*  			  x = a_0 + a_1*alpha + a_2*beta + a_3*gamma + ...           */
/*                        y = ...                                            */
/*                        z = ...                                            */
/*                        (x,y,z)-global coordinates of hexahedron vertecies */
/*                        (alpha,beta,gamma) - local.                        */
/*	                            										     */
/* Input:                 pointer to theElement                              */
/* Output:                extern TransCoeff                                  */
/*****************************************************************************/

#ifdef __THREEDIM__						 
INT TransformCoefficients (ELEMENT *theElement)
{
	int i;
	COORD *x[MAX_CORNERS_OF_ELEM];
	
	for (i=0; i<CORNERS_OF_ELEM(theElement); ++i)
		x[i] = CVECT(MYVERTEX(CORNER(theElement,i)));
		
	V3_COPY(x[0],TransfCoeff[0]);
	
	V3_SUBTRACT(x[1],x[0],TransfCoeff[1]);
	V3_SUBTRACT(x[3],x[0],TransfCoeff[2]);
	V3_SUBTRACT(x[4],x[0],TransfCoeff[3]);
	
	V3_SUBTRACT(x[2],x[1],TransfCoeff[4]);
	V3_SUBTRACT(x[7],x[3],TransfCoeff[5]);
	V3_ADD(TransfCoeff[4],TransfCoeff[5],TransfCoeff[7]);
	
	V3_SUBTRACT(TransfCoeff[4],TransfCoeff[2],TransfCoeff[4]);
	V3_SUBTRACT(x[5],x[1],TransfCoeff[6]);
	V3_SUBTRACT(TransfCoeff[6],TransfCoeff[3],TransfCoeff[6]);	
	V3_SUBTRACT(TransfCoeff[5],TransfCoeff[3],TransfCoeff[5]);
	
	V3_SUBTRACT(TransfCoeff[3],TransfCoeff[7],TransfCoeff[7]);
	V3_SUBTRACT(TransfCoeff[7],x[5],TransfCoeff[7]);
	V3_ADD(TransfCoeff[7],x[6],TransfCoeff[7]);
	
	return (0);
} 
#endif

/*************************************************************************/
/*                                                                       */
/*  JacobiMatrix: computes Jacobi matrix of trasformation                */
/*                  x = Fx(alpha[0],alpha[1],alpha[2]);                  */
/*                  y = Fy(alpha[0],alpha[1],alpha[2]);                  */ 
/*                  z = Fz(alpha[0],alpha[1],alpha[2]);                  */ 
/*                                                                       */
/*  Input:  local coordinates alpha[0],alpha[1],alpha[2].                */
/*  Output: Jacobi matrix  Matrix[9]                                     */
/*                                                                       */
/*************************************************************************/

#ifdef __THREEDIM__						 
INT JacobiMatrix (COORD_VECTOR alpha, COORD *Matrix)
{
	COORD grad[3];
	
	   GradientsOfTransormation (alpha, grad, 0);
	   Matrix[0] = grad[0]; Matrix[3] = grad[1]; Matrix[6] = grad[2];
	   GradientsOfTransormation (alpha, grad, 1);
	   Matrix[1] = grad[0]; Matrix[4] = grad[1]; Matrix[7] = grad[2];
	   GradientsOfTransormation (alpha, grad, 2);
	   Matrix[2] = grad[0]; Matrix[5] = grad[1]; Matrix[8] = grad[2];
	   
	return (0);
}

DOUBLE JacobiMatrixDet (COORD *Matrix)
{
	COORD Inverse[9];
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
	return ((DOUBLE)(Inverse[0]*Matrix[0] + Inverse[3]*Matrix[1] + Inverse[6]*Matrix[2]));	
}

INT GradientsOfTransormation (COORD_VECTOR alpha, COORD_VECTOR grad, int i)
{
	grad[0] = TransfCoeff[1][i] + TransfCoeff[4][i]*alpha[1] + TransfCoeff[6][i]*alpha[2] + TransfCoeff[7][i]*alpha[1]*alpha[2];
	grad[1] = TransfCoeff[2][i] + TransfCoeff[4][i]*alpha[0] + TransfCoeff[5][i]*alpha[2] + TransfCoeff[7][i]*alpha[0]*alpha[2]; 
	grad[2] = TransfCoeff[3][i] + TransfCoeff[5][i]*alpha[1] + TransfCoeff[6][i]*alpha[0] + TransfCoeff[7][i]*alpha[0]*alpha[1]; 

	return (0);
}
	
INT LocalToGlobalHEX (COORD_VECTOR Local, COORD_VECTOR Global)
{
	int i;
	
	for (i=0; i<3; ++i)
		Global[i] = TransfCoeff[0][i] + TransfCoeff[1][i]*Local[0] + 
		                                TransfCoeff[2][i]*Local[1] + 
		                                TransfCoeff[3][i]*Local[2] +
			                        TransfCoeff[4][i]*Local[0]*Local[1] +
			                        TransfCoeff[5][i]*Local[2]*Local[1] +
			                        TransfCoeff[6][i]*Local[2]*Local[0] +
			                        TransfCoeff[7][i]*Local[2]*Local[0]*Local[1];
	return (0);			                        
}	
#endif
	
/*****************************************************************************/
/*                                                                           */
/* GlobalToLocalHEX: finds local coordinates in reference hexahedron         */
/*            Input: Global[3] physical coordinates in global system;        */
/*           Output: Local[3] coordinates in reference hexahedron;           */
/*                   at the beginning Local[3] contains initial guess        */
/* Before calling this program 	TransformCoefficients (ELEMENT *theElement)  */
/* must be called                                                            */
/*                                                                           */
/*****************************************************************************/

#ifdef __THREEDIM__						 	
INT GlobalToLocalHEX (COORD_VECTOR Global, COORD_VECTOR Local)
{
	COORD jm[9],jmi[9];
	COORD tmp[3], rhs[3];
	COORD diff=0.0;
	int iter=0;
	
newiter:
	++iter;
	LocalToGlobalHEX (Local, tmp);
	V3_EUKLIDNORM_OF_DIFF(Global,tmp,diff);
	if ((diff<=1e-6)) return (0);
	if ((iter>50)) return (1);
	
	JacobiMatrix (Local, jm);
	
	M3_TIMES_V3(jm,Local,rhs);
	V3_SUBTRACT(tmp,rhs,tmp);
	V3_SUBTRACT(Global,tmp,rhs);
	
	M3_Invert (jmi,jm);
	M3_TIMES_V3(jmi,rhs,Local);
	goto newiter;
}
#endif

/************************************************************************/
/*                                                                      */
/* HexaDerivatives: computes the derivatives of the shape functions     */
/*                  in global coordinates (x,y,z) in the point          */
/*                  (alpha[0],alpha[1],alpha[2]) in local coordinates.  */
/*                                                                      */
/* Input:local coordinates in the refr.hexah. alpha[0],alpha[1],alpha[2]*/
/* Output: theGradient[MAX_CORNERS_OF_ELEM][3] - gradients.             */
/*                                                                      */
/************************************************************************/

#ifdef __THREEDIM__						 
INT HexaDerivatives (COORD *alpha, COORD_VECTOR theGradient[MAX_CORNERS_OF_ELEM])
{
	int i,j;
	COORD JMat[9], JMatInv[9], r[3];
	
	JacobiMatrix (alpha,JMat);
	
	/* Transpose Jacobi matrix */
	r[0] = JMat[3]; JMat[3] = JMat[1]; JMat[1] = r[0];
	r[0] = JMat[6]; JMat[6] = JMat[2]; JMat[2] = r[0];
	r[0] = JMat[7]; JMat[7] = JMat[5]; JMat[5] = r[0];
	
	M3_Invert (JMatInv,JMat); 
	
	for (i=0; i<8; ++i)    /* 8 = number of hexa corners */
	{
		for (j=0; j<3; ++j)
			r[j] = (COORD) NiHexDer (i, j, alpha);
			
		M3_TIMES_V3(JMatInv,r,theGradient[i]);		
	}	
}

INT FV_HexInfo (ELEMENT *theElement, COORD **theCorners, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], COORD_VECTOR GIP[MAX_EDGES_OF_ELEM])
{
	int i,j;
	COORD CenterOfMassElement[3], CetrerOfMassSides[MAX_SIDES_OF_ELEM][3], emp[MAX_EDGES_OF_ELEM][3];
	COORD a[3], b[3], c[3], sp;
	
	/* compute center of mass of the element */
	V3_CLEAR(CenterOfMassElement)
	for (i=0; i<CORNERS_OF_ELEM(theElement); ++i) 
		V3_ADD(theCorners[i],CenterOfMassElement,CenterOfMassElement)
	V3_SCALE(0.125,CenterOfMassElement)
	
	/* compute centers of mass of the sides */
	for (i=0; i<SIDES_OF_ELEM(theElement); ++i)
	{
		V3_CLEAR(CetrerOfMassSides[i])
		for (j=0; j<CORNERS_OF_SIDE(theElement,i); ++j)
			V3_ADD(theCorners[CORNER_OF_SIDE(theElement,i,j)],CetrerOfMassSides[i],CetrerOfMassSides[i])
		V3_SCALE(0.25,CetrerOfMassSides[i])
	} 
	
	for (i=0; i<EDGES_OF_ELEM(theElement); ++i)
	{
		V3_LINCOMB(0.5,theCorners[CORNER_OF_EDGE(theElement,i,0)],0.5,theCorners[CORNER_OF_EDGE(theElement,i,1)],emp[i])
		V3_SUBTRACT(CenterOfMassElement,emp[i],a)
		V3_SUBTRACT(CetrerOfMassSides[SIDE_WITH_EDGE(theElement,i,0)],emp[i],b)
		V3_SUBTRACT(CetrerOfMassSides[SIDE_WITH_EDGE(theElement,i,1)],emp[i],c)
		
		V3_VECTOR_PRODUCT(a,b,Area[i])
		V3_VECTOR_PRODUCT(c,a,b)
		V3_LINCOMB (0.5,Area[i],0.5,b,Area[i])
		
		V3_SUBTRACT (theCorners[CORNER_OF_EDGE(theElement,i,1)],theCorners[CORNER_OF_EDGE(theElement,i,0)],a)
		V3_SCALAR_PRODUCT(Area[i],a,sp)
		if (sp<0.0){
			Area[i][0]=-Area[i][0];  Area[i][1]=-Area[i][1];  Area[i][2]=-Area[i][2];
		}
	}
	
	/*  compute GIP */
	/* ...          */
	
	return (0);
}
#endif

/***********************************************************************/
/*                                                                     */
/*    HexaVolume: computes hexahedron volume by computing volumes      */
/*                of 5 tetrahedrons which form it.                     */
/*                                                                     */
/***********************************************************************/

#ifdef __THREEDIM__						 
INT HexaVolume (COORD **theCorners, COORD *volume)
{
	COORD *Tcorners[4]; /* 4 = corners of tetrahedron */
	COORD Tvolume;
	
	*volume = 0.0;
	
	Tcorners[0] = theCorners[0]; Tcorners[2] = theCorners[3];
	Tcorners[1] = theCorners[1]; Tcorners[3] = theCorners[4];
	TetraVolume(Tcorners,&Tvolume);
	*volume += Tvolume;

	Tcorners[0] = theCorners[5]; Tcorners[2] = theCorners[4];
	Tcorners[1] = theCorners[6]; Tcorners[3] = theCorners[1];
	TetraVolume(Tcorners,&Tvolume);
	*volume += Tvolume;

	Tcorners[0] = theCorners[2]; Tcorners[2] = theCorners[3];
	Tcorners[1] = theCorners[1]; Tcorners[3] = theCorners[6];
	TetraVolume(Tcorners,&Tvolume);
	*volume += Tvolume;

	Tcorners[0] = theCorners[7]; Tcorners[2] = theCorners[4];
	Tcorners[1] = theCorners[3]; Tcorners[3] = theCorners[6];
	TetraVolume(Tcorners,&Tvolume);
	*volume += Tvolume;

	Tcorners[0] = theCorners[3]; Tcorners[2] = theCorners[1];
	Tcorners[1] = theCorners[6]; Tcorners[3] = theCorners[4];
	TetraVolume(Tcorners,&Tvolume);
	*volume += Tvolume;
	
	return (0);
}
#endif

/***************************************************************************/
/*                                                                         */
/*     TransformGlobalToLocal3D: computes transformation from global to    */
/*                               local coordinates.                        */
/*                                                                         */
/*     No difference between tetrahedrons ans hexahedrons!                 */
/*                                                                         */
/***************************************************************************/

#ifdef __THREEDIM__						 
INT TransformGlobalToLocal3D(ELEMENT *theElement, COORD_VECTOR Global, COORD_VECTOR Local)
{
    INT i;
    COORD *x[MAX_CORNERS_OF_ELEM];

    if(TAG(theElement)==TETRAHEDRON) {
        for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
            x[i] = CVECT(MYVERTEX(CORNER(theElement,i)));
        return(GlobalToLocal3d(x,Global,Local));
    }
    if(TAG(theElement)==HEXAHEDRON){
        TransformCoefficients (theElement);
        return(GlobalToLocalHEX(Global,Local));
    }
}
#endif

/****************************************************************************/
/*D
   GlobalToLocal3d - Transform global coordinates to local

   SYNOPSIS:
   INT GlobalToLocal3d (const COORD **Corners, const COORD *EvalPoint, COORD *LocalCoord);

   PARAMETERS:
.  Corners - coordinates of corners 
.  EvalPoint - global coordinates
.  LocalCoord - local coordinates 

   DESCRIPTION:
   This function transforms global coordinates to local in an evaluated point
   in 3D.

   RETURN VALUE:
   INT
.n    0 if ok
.n    1 if error occured.
D*/
/****************************************************************************/

#ifdef __THREEDIM__						 
INT GlobalToLocal3d (const COORD **Corners, const COORD *EvalPoint, COORD *LocalCoord)
{
	COORD_VECTOR a;
	COORD M[9],I[9];

	V3_SUBTRACT(EvalPoint,Corners[0],a)
	V3_SUBTRACT(Corners[1],Corners[0],M)
	V3_SUBTRACT(Corners[2],Corners[0],M+3)
	V3_SUBTRACT(Corners[3],Corners[0],M+6)
	if (M3_Invert(I,M)) return (1);
	M3_TIMES_V3(I,a,LocalCoord)

	return(0);
}
#endif

/****************************************************************************/
/*D
   TetraSideNormals - Calculate inner normals of tetrahedra

   SYNOPSIS:
   INT TetraSideNormals (ELEMENT *theElement, COORD **theCorners, COORD_VECTOR theNormals[MAX_SIDES_OF_ELEM]);

   PARAMETERS:
.  theCorners - list of pointers to phys corner vectors
.  theNormals - normals of tetrahedra

   DESCRIPTION:
   This function calculates the inner normals on the sides of a tetrahedron.	

   RETURN VALUE:
   INT
.n    0 if ok
.n    1 if error occured.
D*/
/****************************************************************************/

#ifdef __THREEDIM__						 
INT TetraSideNormals (ELEMENT *theElement, COORD **theCorners, COORD_VECTOR theNormals[MAX_SIDES_OF_ELEM])
{
   	ELEMENT e;
	COORD_VECTOR a, b;
	COORD h;
	INT j,k;

	/* TODO: changed MAX_CORNERS_OF_ELEM to 4 and subsequently*/
	SETTAG(&e,4);
	for (j=0; j<4; j++)
	{
		k = SIDE_OPP_TO_CORNER(&e,j);
		V3_SUBTRACT(theCorners[(j+1)%4],theCorners[(j+2)%4],a)
		V3_SUBTRACT(theCorners[(j+1)%4],theCorners[(j+3)%4],b)
		V3_VECTOR_PRODUCT(a,b,theNormals[k])
		V3_Normalize(theNormals[k]);
		V3_SUBTRACT(theCorners[j],theCorners[(j+1)%4],a)
		V3_SCALAR_PRODUCT(theNormals[k],a,h);
		if (ABS(h)<SMALL_C) return (1);
		if (h<0.0)
			V3_SCALE(-1.0,theNormals[k]);
	}

	return (0);

}
#endif

/****************************************************************************/
/*D
   TetMaxSideAngle - Calculate maximal side angle of Tetrahedron

   SYNOPSIS:
   INT TetMaxSideAngle (ELEMENT *theElement, COORD **theCorners, COORD *MaxAngle);

   PARAMETERS:
.  theCorners - list of pointers to phys corner vectors
.  MaxAngle - max side angle

   DESCRIPTION:
   This function calculates the maximal side angle of the tetrahedron.

   RETURN VALUE:
   INT
.n    0 if ok
.n    1 if error occured.
D*/
/****************************************************************************/

#ifdef __THREEDIM__						 
INT TetMaxSideAngle (ELEMENT *theElement, const COORD **theCorners, COORD *MaxAngle)
{
	COORD_VECTOR theNormal[MAX_SIDES_OF_ELEM];
	COORD max,help;
	INT i;

	if (TetraSideNormals (theElement,theCorners,theNormal)) return (1);
	max = -1.0;
	for (i=0; i<EDGES_OF_ELEM(theElement); i++)
	{
		V3_SCALAR_PRODUCT(theNormal[SIDE_WITH_EDGE(theElement,i,0)],theNormal[SIDE_WITH_EDGE(theElement,i,1)],help)
		max = MAX(help,max);
	}
	max = MIN(max,1.0);
	*MaxAngle = 180.0/PI*acos(-max);
	
	return (0);
}
#endif

/****************************************************************************/
/*D
   TetAngleAndLength - Calculates side angle and length of edge of Tetrahedron

   SYNOPSIS:
   INT TetAngleAndLength (ELEMENT *theElement, COORD **theCorners, COORD *Angle, COORD *Length);

   PARAMETERS:
.  theCorners - list of pointers to phys corner vectors
.  Angle - side angle
.  Length - sidelength

   DESCRIPTION:
   This function calculates the side angle of a tetrahedron and the lengths of 
   the edges belonging to this side.

   RETURN VALUE:
   INT
.n    0 if ok
.n    1 if error occured.
D*/
/****************************************************************************/

#ifdef __THREEDIM__						 
INT TetAngleAndLength (ELEMENT *theElement, const COORD **theCorners, COORD *Angle, COORD *Length)
{
	COORD_VECTOR theNormals[MAX_SIDES_OF_ELEM],theEdge[MAX_EDGES_OF_ELEM];
	COORD h;
	INT j,k;

	for (j=0; j<EDGES_OF_ELEM(theElement); j++)
	{
		V3_SUBTRACT(theCorners[CORNER_OF_EDGE(theElement,j,1)],theCorners[CORNER_OF_EDGE(theElement,j,0)],theEdge[j])
		V3_EUKLIDNORM(theEdge[j],Length[j])
	}
	for (j=0; j<SIDES_OF_ELEM(theElement); j++)
	{
		V3_VECTOR_PRODUCT(theEdge[EDGE_OF_SIDE(theElement,j,0)],theEdge[EDGE_OF_SIDE(theElement,j,1)],theNormals[j])
		V3_Normalize(theNormals[j]);
		k = EDGE_OF_CORNER(theElement,CORNER_OPP_TO_SIDE(theElement,j),0);
		V3_SCALAR_PRODUCT(theNormals[j],theEdge[k],h)
		if (ABS(h)<SMALL_C) return (1);
		if ( (h<0.0 && CORNER_OF_EDGE(theElement,k,1)==CORNER_OPP_TO_SIDE(theElement,j)) ||
			 (h>0.0 && CORNER_OF_EDGE(theElement,k,0)==CORNER_OPP_TO_SIDE(theElement,j))	 )
			V3_SCALE(-1.0,theNormals[j]);
	}
	for (j=0; j<EDGES_OF_ELEM(theElement); j++)
	{
		V3_SCALAR_PRODUCT(theNormals[SIDE_WITH_EDGE(theElement,j,0)],theNormals[SIDE_WITH_EDGE(theElement,j,1)],Angle[j])
		Angle[j] = MAX(Angle[j],-1.0);
		Angle[j] = MIN(Angle[j], 1.0);
		Angle[j] = (COORD)acos((double)Angle[j]);
	}
	
	return (0);
}
#endif

/****************************************************************************/
/*D
   TetraDerivative - Calculates gradient of shape function for tetrahedron

   SYNOPSIS:
   INT TetraDerivative (ELEMENT *theElement, COORD **theCorners, COORD_VECTOR theGradient[MAX_CORNERS_OF_ELEM])

   PARAMETERS:
.  theCorners - list of pointers to phys corner vectors
.  theGradient - gradients of shape functions

   DESCRIPTION:
   This function calculates the gradient of shape function for tetrahedron.

   RETURN VALUE:
   INT
.n    0 if ok
.n    1 if error occured.
D*/
/****************************************************************************/

#ifdef __THREEDIM__						 
INT TetraDerivative (ELEMENT *theElement, const COORD **theCorners, COORD_VECTOR theGradient[MAX_CORNERS_OF_ELEM])
{
	COORD_VECTOR a, b;
	COORD h;
	INT j;

	for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
	{
		V3_SUBTRACT(theCorners[(j+1)%4],theCorners[(j+2)%4],a)
		V3_SUBTRACT(theCorners[(j+1)%4],theCorners[(j+3)%4],b)
		V3_VECTOR_PRODUCT(a,b,theGradient[j])
		V3_Normalize(theGradient[j]);
		V3_SUBTRACT(theCorners[j],theCorners[(j+1)%4],a)
		V3_SCALAR_PRODUCT(theGradient[j],a,h)
		if (ABS(h)<SMALL_C) return (1);
		V3_SCALE(1/h,theGradient[j])
	}

	return (0);
}
#endif

/****************************************************************************/
/*D
   TetraVolume - Calculate volume of tetrahedron

   SYNOPSIS:
   INT TetraVolume (COORD **theCorners, COORD *volume);

   PARAMETERS:
.  theCorners - list of pointers to phys corner vectors
.  volume - volume of tetrahedron

   DESCRIPTION:
   This function calculates the volume of a tetrahedron.	

   RETURN VALUE:
   INT
.n    0 if ok
.n    1 if error occured.
D*/
/****************************************************************************/

#ifdef __THREEDIM__						 
INT TetraVolume (const COORD **theCorners, COORD *volume)
{
	COORD_VECTOR a, b, n;

	/* compute volume of tetrahedron */
	V3_SUBTRACT(theCorners[0],theCorners[1],a)
	V3_SUBTRACT(theCorners[0],theCorners[2],b)
	V3_VECTOR_PRODUCT(a,b,n)
	V3_SUBTRACT(theCorners[0],theCorners[3],a)
	V3_SCALAR_PRODUCT(n,a,volume[0])
	*volume = ABS(*volume)/6;
	
	return (0); 
}
#endif

/****************************************************************************/
/*D
   FV_TetInfo - Calculate control volume subsurfaces and global integration points

   SYNOPSIS:
   INT FV_TetInfo (const COORD **theCorners, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], 
   COORD_VECTOR GIP[MAX_EDGES_OF_ELEM]);

   PARAMETERS:
.  theCorners - coordinates of the element corners
.  Area - area of subsurface
.  GIP - global integration point

   DESCRIPTION:
   This function calculates the subsurfaces (SCV) of the control volume (CV) in 
   tetrahedron elements. The global integration points are computed in the 
   center of mass on these subsurfaces.

   RETURN VALUE:
   INT
.n    0 if ok
.n    1 if error occured.
D*/
/****************************************************************************/

#ifdef __THREEDIM__						 
INT FV_TetInfo (const COORD **theCorners, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], COORD_VECTOR GIP[MAX_EDGES_OF_ELEM])
{
   	ELEMENT e;
	COORD_VECTOR emp[MAX_EDGES_OF_ELEM], diff, a, b;
	COORD sp;
	INT i;

	/* TODO: changed MAX_EDGES_OF_ELEM to 6 */
	SETTAG(&e,4);
	for (i=0; i<6; i++)
	{
		V3_LINCOMB(0.5,theCorners[CORNER_OF_EDGE(&e,i,0)],0.5,theCorners[CORNER_OF_EDGE(&e,i,1)],emp[i])
		V3_SUBTRACT (theCorners[CORNER_OF_OPPEDGE(&e,i,0)],emp[i],a)
		V3_SUBTRACT (theCorners[CORNER_OF_OPPEDGE(&e,i,1)],emp[i],b)
		V3_VECTOR_PRODUCT(a,b,Area[i])
		V3_SUBTRACT (theCorners[CORNER_OF_EDGE(&e,i,1)],theCorners[CORNER_OF_EDGE(&e,i,0)],diff)
		V3_SCALAR_PRODUCT(Area[i],diff,sp)
		if (sp>0.0)
			V3_SCALE(1/12.0,Area[i])
		else
			V3_SCALE(-1/12.0,Area[i])
	}
	/* TODO: changed MAX_EDGES_OF_ELEM to 6 */
	for (i=0; i<6; i++)
		V3_LINCOMB(17.0/24.0,emp[i],7.0/24.0,emp[OPPOSITE_EDGE(&e,i)],GIP[i])

	return (0);
}
#endif

/****************************************************************************/
/*																			*/
/* Function:  FV_AliTetInfo													*/
/*																			*/
/* Purpose:   calc gradient of shape function i 							*/
/*																			*/
/* Input:	  int i: corner number [0..3]									*/
/*			  COORD **Corners: list of ptrs to phys corner vectors			*/
/*																			*/
/* Output:	  COORD **theGradient: list of ptrs to gradients				*/
/*																			*/
/* Return:	  INT 0: ok 													*/
/*				  1: error													*/
/*																			*/
/****************************************************************************/

#ifdef __THREEDIM__						 
static INT FindCrossParam3D (COORD_VECTOR p1, COORD_VECTOR p2, COORD_VECTOR p3, COORD_VECTOR p4, DOUBLE_VECTOR v, COORD_VECTOR param)
{
	COORD M[9], I[9];
	
	V3_SUBTRACT(p1,p2,M)
	V3_SUBTRACT(p4,p3,M+3)
	V3_COPY(v,M+6)
	if (M3_Invert(I,M))
		return (1);
	V3_SUBTRACT(p1,p3,M)
	M3_TIMES_V3(I,M,param)
	
	return (0);
}

static INT MirrorAtPlane (const COORD *in, const COORD *pp, const COORD *pn, COORD *out)
{
	COORD_VECTOR a;
	
	V3_SUBTRACT(pp,in,a)
	if (V3_Project(a,pn,out))
	V3_LINCOMB(1.0,in,2.0,out,out)
	
	return (0);
}

#define V3_TRI_CM(a,b,c,e)	 {(e)[0] = 0.333333333333*((a)[0]+(b)[0]+(c)[0]);\
					     	  (e)[1] = 0.333333333333*((a)[1]+(b)[1]+(c)[1]);\
							  (e)[2] = 0.333333333333*((a)[2]+(b)[2]+(c)[2]);}
#define V3_QUA_CM(a,b,c,d,e) {(e)[0] = 0.25*((a)[0]+(b)[0]+(c)[0]+(d)[0]);\
							  (e)[1] = 0.25*((a)[1]+(b)[1]+(c)[1]+(d)[1]);\
							  (e)[2] = 0.25*((a)[2]+(b)[2]+(c)[2]+(d)[2]);}

static INT CornerOfSideAndEdge[4][6] = {
   { 2, 0, 1,-1,-1,-1},
   {-1, 3,-1,-1, 2, 1},
   {-1,-1, 3, 2,-1, 0},
   { 3,-1,-1, 1, 0,-1}
};

static INT SideOfCorners[4][4][4] = { 
                  {{-1,-1,-1,-1}, {-1,-1, 0, 3}, {-1, 0,-1, 2}, {-1, 3, 2,-1}},
				  {{-1,-1, 0, 3}, {-1,-1,-1,-1}, { 0,-1,-1, 1}, { 3,-1, 1,-1}},
				  {{-1, 0,-1, 2}, { 0,-1,-1, 1}, {-1,-1,-1,-1}, { 2, 1,-1,-1}},
				  {{-1, 3, 2,-1}, { 3,-1, 1,-1}, { 2, 1,-1,-1}, {-1,-1,-1,-1}}
};

/* TODO: change to macros of the general element concept */

static COORD_VECTOR TRefCoord[4] = 
{{0.0,0.0,0.0},{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};

static INT  CornerOfSide[4][3] = {{0,2,1}, {1,2,3}, {2,0,3}, {3,0,1}};
/* gives the position of specified corner on specified side or -1 if it does not ly on that side */
INT  CornerOfSideInv[4][4]  = {{0,2,1,-1}, {-1,0,1,2}, {1,-1,0,2}, {1,2,-1,0}};

/* The indices of the corner opposite to three corners. */
static INT  CornerOppCorners[4][4][4] =
				{ {{-1,-1,-1,-1}, {-1,-1, 3, 2}, {-1, 3,-1, 1}, {-1, 2, 1,-1}}, 
				  {{-1,-1, 3, 2}, {-1,-1,-1,-1}, { 3,-1,-1, 0}, { 2,-1, 0,-1}},
				  {{-1, 3,-1, 1}, { 3,-1,-1, 0}, {-1,-1,-1,-1}, { 1, 0,-1,-1}},
				  {{-1, 2, 1,-1}, { 2,-1, 0,-1}, { 1, 0,-1,-1}, {-1,-1,-1,-1}} };

/* The indices of the corners of each edge and its opposite edge */
/* dont change order !											 */
static INT  CornerOfEdge[6][2]	 = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};
static INT  CornerOfOppEdge[6][2] = {{2,3},{0,3},{3,1},{1,2},{2,0},{0,1}};

/* the indices of the edges between two corners */
static INT  EdgeWithCorners[4][4] = {{-1,0,2,3},{0,-1,1,4},
												  {2,1,-1,5},{3,4,5,-1}};

/* the indices of the sides around each edge */
static INT  SideWithEdge[6][MAX_SIDES_OF_EDGE] = {{0,3},{0,1},{0,2},{2,3},{1,3},{1,2}};

/* the index of the side having to specified edges */
static INT  SideOf2Edges[6][6] = {{-1,0,0,3,3,-1},{0,-1,0,-1,1,1},{0,0,-1,2,-1,2},{3,-1,2,-1,3,2},{3,1,-1,3,-1,1},{-1,1,2,2,1,-1}};

/* the index of the edge having to specified sides */
static INT  EdgeOf2Sides[4][4] = {{-1,1,2,0},{1,-1,5,4},{2,5,-1,3},{0,4,3,-1}};

/* the indices of the edges of each side */
static INT  EdgeOfSide[4][3] = {{0,1,2},{1,5,4},{2,3,5},{0,4,3}};
static INT  CondensedEdgeOfSide[4] = {0x07,0x32,0x2C,0x19};

/* the indices of opposite corners for each side */
static INT OppositeCorner[4] = {3,0,1,2};

/* the indices of opposite sides for each corner */
static INT OppositeSide[4] = {1,2,3,0};

INT FV_AliTetInfo (const COORD **CornerPoints, COORD_VECTOR Area[6], DOUBLE_VECTOR conv, COORD_VECTOR GIP[6], COORD_VECTOR LIP[6])
{
	COORD sp, alpha, check[2], M[9], Inv[9];
	COORD_VECTOR a, b, c, d, e, cm, normal, param, EdgeMidPoints[6], SideMidPoints[4];
	INT i, help, noutflow, ninflow, outflow[4], inflow[4], OpEdge[3], GrEdge[3], side[3], OpCorner, corner[3], inverted, First;
	INT BackEdge, FrontEdge, BackCorner[2], FrontCorner[2], EdgeF0B0, EdgeF0B1, EdgeF1B0, EdgeF1B1, flags, changed;
	
	/* reset areas and integrationpoints */
	for (i=0; i<6; i++)
	{
		V3_CLEAR(Area[i])
		V3_CLEAR(GIP[i])
		V3_COPY(LIP[i],LIP[i])
	}
	
	/* edge mid points */
	for (i=0; i<6; i++)	V3_LINCOMB(0.5,CornerPoints[CornerOfEdge[i][0]],0.5,CornerPoints[CornerOfEdge[i][1]],EdgeMidPoints[i])
			
	/* side mid points */
	for (i=0; i<4; i++)
	{
		V3_ADD(CornerPoints[CornerOfSide[i][0]],CornerPoints[CornerOfSide[i][1]],a)
		V3_ADD(CornerPoints[CornerOfSide[i][2]],a,SideMidPoints[i])
		V3_SCALE(0.33333333333333333,SideMidPoints[i])
	}

	/* in/outflow sides */
	noutflow = ninflow = 0;
	for (i=0; i<4; i++)
	{
		V3_SUBTRACT(CornerPoints[CornerOfSide[i][1]],CornerPoints[CornerOfSide[i][0]],a)
		V3_SUBTRACT(CornerPoints[CornerOfSide[i][2]],CornerPoints[CornerOfSide[i][0]],b)
		V3_VECTOR_PRODUCT(a,b,normal)
		V3_SUBTRACT(CornerPoints[OppositeCorner[i]],CornerPoints[CornerOfSide[i][0]],a)
		V3_SCALAR_PRODUCT(a,normal,sp)
		if (sp>0.0) V3_SCALE(-1.0,normal)
		V3_SCALAR_PRODUCT(conv,normal,sp)
		if (sp>=0.0) 	outflow[noutflow++] = i;
		if (sp<0.0)		inflow[ninflow++] = i;
	}
	
	/* check */
	check[0]=check[1]=0.0;
	for (i=0; i<ninflow; i++)
	{
		V3_SUBTRACT(CornerPoints[CornerOfSide[inflow[i]][1]],CornerPoints[CornerOfSide[inflow[i]][0]],a)
		V3_SUBTRACT(CornerPoints[CornerOfSide[inflow[i]][2]],CornerPoints[CornerOfSide[inflow[i]][0]],b)
		V3_VECTOR_PRODUCT(a,b,normal)
		V3_SCALAR_PRODUCT(conv,normal,sp)
		check[0] += ABS(0.5*sp);
	}
	for (i=0; i<noutflow; i++)
	{
		V3_SUBTRACT(CornerPoints[CornerOfSide[outflow[i]][1]],CornerPoints[CornerOfSide[outflow[i]][0]],a)
		V3_SUBTRACT(CornerPoints[CornerOfSide[outflow[i]][2]],CornerPoints[CornerOfSide[outflow[i]][0]],b)
		V3_VECTOR_PRODUCT(a,b,normal)
		V3_SCALAR_PRODUCT(conv,normal,sp)
		check[1] += ABS(0.5*sp);
	}
	assert((check[0]-check[1])<=1e-6*(check[0]+check[1]));

	/* change if inflow > outflow */
	inverted = 0;
	if (ninflow>noutflow)
	{
		inverted = 1;
		for (i=0; i<ninflow; i++)	{help = inflow[i]; inflow[i] = outflow[i]; outflow[i] = help;}
		help = ninflow; ninflow = noutflow; noutflow = help;
	}
	if ( ninflow+noutflow!=4 || ninflow<1 || ninflow>2) 
		return (1);
		
	/* handle different cases */
	switch (ninflow)
	{
		case 2:
			/* two in-, two outflow */
			BackEdge = EdgeOf2Sides[inflow[0]][inflow[1]];
			FrontEdge = EdgeOf2Sides[outflow[0]][outflow[1]];
			BackCorner[0] = CornerOfSideAndEdge[outflow[0]][FrontEdge];
			BackCorner[1] = CornerOfSideAndEdge[outflow[1]][FrontEdge];
			FrontCorner[0] = CornerOfSideAndEdge[inflow[0]][BackEdge];
			FrontCorner[1] = CornerOfSideAndEdge[inflow[1]][BackEdge];
			EdgeF0B0 = EdgeWithCorners[FrontCorner[0]][BackCorner[0]];
			EdgeF0B1 = EdgeWithCorners[FrontCorner[0]][BackCorner[1]];
			EdgeF1B0 = EdgeWithCorners[FrontCorner[1]][BackCorner[0]];
			EdgeF1B1 = EdgeWithCorners[FrontCorner[1]][BackCorner[1]];
			
			if (FindCrossParam3D(CornerPoints[FrontCorner[0]],CornerPoints[FrontCorner[1]],CornerPoints[BackCorner[0]],CornerPoints[BackCorner[1]],conv,param)) return (1);
			if (param[0]<0.0 || param[0]>1.0 || param[1]<0.0 || param[1]>1.0)
				i=i;
			changed=0;
			if (param[0]<0.5)
			{
				param[0]=1.0-param[0];
				help=inflow[0]; inflow[0]=inflow[1]; inflow[1]=help;
				changed=1;
			}
			if (param[1]<0.5)
			{
				param[1]=1.0-param[1];
				help=outflow[0]; outflow[0]=outflow[1]; outflow[1]=help;
				changed=1;
			}
			if (changed)
			{
				BackCorner[0] = CornerOfSideAndEdge[outflow[0]][FrontEdge];
				BackCorner[1] = CornerOfSideAndEdge[outflow[1]][FrontEdge];
				FrontCorner[0] = CornerOfSideAndEdge[inflow[0]][BackEdge];
				FrontCorner[1] = CornerOfSideAndEdge[inflow[1]][BackEdge];
				EdgeF0B0 = EdgeWithCorners[FrontCorner[0]][BackCorner[0]];
				EdgeF0B1 = EdgeWithCorners[FrontCorner[0]][BackCorner[1]];
				EdgeF1B0 = EdgeWithCorners[FrontCorner[1]][BackCorner[0]];
				EdgeF1B1 = EdgeWithCorners[FrontCorner[1]][BackCorner[1]];
			}
			flags = (param[0]>0.5) | ((param[1]>0.5)<<1);
			switch (flags)
			{
				case 0:
					/* totaly symmetric */
					V3_SUBTRACT(CornerPoints[BackCorner[1]],CornerPoints[BackCorner[0]],a)
					V3_SUBTRACT(CornerPoints[FrontCorner[1]],CornerPoints[FrontCorner[0]],b)
					V3_VECTOR_PRODUCT(a,b,c)
					V3_SCALE(0.041666666666666666,c) V3_COPY(c,d)
					V3_SCALAR_PRODUCT(conv,c,sp)
					if (sp>0.0) V3_SCALE(-1.0,d)
					else 		{V3_SCALE(-1.0,c); sp = -sp;}
					
					V3_QUA_CM(SideMidPoints[inflow[0]],SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF0B0],EdgeMidPoints[BackEdge],GIP[EdgeF0B0])
					if (CornerOfEdge[EdgeF0B0][0]==FrontCorner[0])		V3_COPY(d,Area[EdgeF0B0])
					else												V3_COPY(c,Area[EdgeF0B0])
			
					V3_QUA_CM(SideMidPoints[inflow[0]],SideMidPoints[outflow[1]],EdgeMidPoints[EdgeF0B1],EdgeMidPoints[BackEdge],GIP[EdgeF0B1])
					if (CornerOfEdge[EdgeF0B1][0]==FrontCorner[0])		V3_COPY(d,Area[EdgeF0B1])
					else												V3_COPY(c,Area[EdgeF0B1])
					
					V3_QUA_CM(SideMidPoints[inflow[1]],SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF1B0],EdgeMidPoints[BackEdge],GIP[EdgeF1B0])
					if (CornerOfEdge[EdgeF1B0][0]==FrontCorner[1])		V3_COPY(d,Area[EdgeF1B0])
					else												V3_COPY(c,Area[EdgeF1B0])

					V3_QUA_CM(SideMidPoints[inflow[1]],SideMidPoints[outflow[1]],EdgeMidPoints[EdgeF1B1],EdgeMidPoints[BackEdge],GIP[EdgeF1B1])
					if (CornerOfEdge[EdgeF1B1][0]==FrontCorner[1])		V3_COPY(d,Area[EdgeF1B1])
					else												V3_COPY(c,Area[EdgeF1B1])
					break;
				case 1:
					/* FrontEdgeMidPoint in 'inflow[0]' ,symmetric in 'outflow' */
					if (param[0]<=0.75)
					{
						/* no coupling of F0 and F1 */
						
						/* entries F1B0 and F1B1 */
						V3_SUBTRACT(EdgeMidPoints[FrontEdge],EdgeMidPoints[EdgeF1B0],a)
						V3_SUBTRACT(SideMidPoints[inflow[1]],SideMidPoints[outflow[0]],b)
						V3_VECTOR_PRODUCT(a,b,c)
						V3_SCALE(0.5,c) V3_COPY(c,d)
						V3_SCALAR_PRODUCT(conv,c,sp)
						if (sp>0.0) V3_SCALE(-1.0,d)
						else 		{V3_SCALE(-1.0,c); sp = -sp;}
						
						V3_QUA_CM(SideMidPoints[inflow[1]],SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF1B0],EdgeMidPoints[FrontEdge],GIP[EdgeF1B0])
						if (CornerOfEdge[EdgeF1B0][0]==FrontCorner[1])	V3_COPY(d,Area[EdgeF1B0])
						else											V3_COPY(c,Area[EdgeF1B0])
						
						V3_QUA_CM(SideMidPoints[inflow[1]],SideMidPoints[outflow[1]],EdgeMidPoints[EdgeF1B1],EdgeMidPoints[FrontEdge],GIP[EdgeF1B1])
						if (CornerOfEdge[EdgeF1B1][0]==FrontCorner[1])	V3_COPY(d,Area[EdgeF1B1])
						else											V3_COPY(c,Area[EdgeF1B1])
						
						/* entries F0B0 and F0B1 */
						V3_SUBTRACT(EdgeMidPoints[FrontEdge],EdgeMidPoints[EdgeF0B0],a)
						V3_SUBTRACT(SideMidPoints[inflow[0]],SideMidPoints[outflow[0]],b)
						V3_VECTOR_PRODUCT(a,b,c)
						V3_SCALE(0.5,c) V3_COPY(c,d)
						V3_SCALAR_PRODUCT(conv,c,sp)
						if (sp>0.0) V3_SCALE(-1.0,d)
						else 		{V3_SCALE(-1.0,c); sp = -sp;}
						
						V3_QUA_CM(SideMidPoints[inflow[0]],SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF0B0],EdgeMidPoints[FrontEdge],GIP[EdgeF0B0])
						if (CornerOfEdge[EdgeF0B0][0]==FrontCorner[0])	V3_COPY(d,Area[EdgeF0B0])
						else											V3_COPY(c,Area[EdgeF0B0])

						V3_QUA_CM(SideMidPoints[inflow[0]],SideMidPoints[outflow[1]],EdgeMidPoints[EdgeF0B1],EdgeMidPoints[FrontEdge],GIP[EdgeF0B1])
						if (CornerOfEdge[EdgeF0B1][0]==FrontCorner[0])	V3_COPY(d,Area[EdgeF0B1])
						else											V3_COPY(c,Area[EdgeF0B1])
					}
					else
					{
						/* coupling of F0 and F1 */
						
						/* entries for F0F1, F0B0, F0B1 */
						if (FindCrossParam3D(EdgeMidPoints[FrontEdge],SideMidPoints[outflow[0]],SideMidPoints[inflow[0]],EdgeMidPoints[EdgeF0B0],conv,param)) return (1);
						V3_SUBTRACT(EdgeMidPoints[FrontEdge],SideMidPoints[outflow[0]],a)
						V3_SUBTRACT(SideMidPoints[inflow[0]],EdgeMidPoints[EdgeF0B0],b)
						V3_VECTOR_PRODUCT(a,b,c) V3_SCALAR_PRODUCT(conv,c,sp) if (sp<0.0) {V3_SCALE(-1.0,c); sp = -sp;}

						V3_LINCOMB(0.5,SideMidPoints[inflow[0]],0.5,EdgeMidPoints[FrontEdge],GIP[FrontEdge])
						V3_COPY(c,Area[FrontEdge])
						if (CornerOfEdge[FrontEdge][0]==FrontCorner[1]) alpha = -param[0]*param[1];
						else											alpha = param[0]*param[1];
						V3_SCALE(alpha,Area[FrontEdge])

						V3_LINCOMB(1.0-param[0],EdgeMidPoints[FrontEdge],param[0],SideMidPoints[outflow[0]],a)
						V3_TRI_CM(SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF0B0],a,GIP[EdgeF0B0])
						V3_COPY(c,Area[EdgeF0B0])
						if (CornerOfEdge[EdgeF0B0][0]==FrontCorner[0])  alpha = 0.5*(1.0-param[0])*(param[1]-1.0);
						else											alpha = 0.5*(1.0-param[0])*(1.0-param[1]);
						V3_SCALE(alpha,Area[EdgeF0B0])
						V3_COPY(c,Area[EdgeF0B1])

							
						V3_SUBTRACT(CornerPoints[FrontCorner[1]],CornerPoints[FrontCorner[0]],a)
						V3_SUBTRACT(EdgeMidPoints[BackEdge],CornerPoints[FrontCorner[0]],b)
						V3_VECTOR_PRODUCT(a,b,d)
						if (MirrorAtPlane(GIP[EdgeF0B0],CornerPoints[FrontCorner[0]],d,GIP[EdgeF0B1])) return (1);
						if (CornerOfEdge[EdgeF0B1][0]==FrontCorner[0])  alpha = 0.5*(1.0-param[0])*(param[1]-1.0);
						else											alpha = 0.5*(1.0-param[0])*(1.0-param[1]);
						V3_SCALE(alpha,Area[EdgeF0B1])
						
						/* entries for F1B0, F1B1 */
						V3_COPY(c,d) alpha=0.5*(1.0-param[0])*param[1]; V3_SCALE(alpha,d)

						V3_LINCOMB(1.0-param[0],EdgeMidPoints[FrontEdge],param[0],SideMidPoints[outflow[0]],a)
						V3_TRI_CM(SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF0B0],a,GIP[EdgeF1B0])
						V3_SCALE(sp*alpha,GIP[EdgeF1B0])
					
						V3_SUBTRACT(SideMidPoints[inflow[0]],EdgeMidPoints[EdgeF1B0],a)
						V3_SUBTRACT(SideMidPoints[inflow[1]],SideMidPoints[outflow[0]],b)
						V3_VECTOR_PRODUCT(a,b,c) V3_SCALAR_PRODUCT(conv,c,sp) if (sp<0.0) {V3_SCALE(-0.5,c);sp=-sp;} else V3_SCALE(0.5,c)
						V3_ADD(c,d,c) V3_COPY(c,d) V3_SCALE(-1.0,d)
							
						
						V3_QUA_CM(SideMidPoints[outflow[0]],SideMidPoints[inflow[1]],EdgeMidPoints[EdgeF1B0],SideMidPoints[inflow[0]],cm)
						V3_LINCOMB(1.0,GIP[EdgeF1B0],0.5*sp,cm,GIP[EdgeF1B0])
						if (CornerOfEdge[EdgeF1B0][0]==FrontCorner[1])	V3_COPY(d,Area[EdgeF1B0])
						else											V3_COPY(c,Area[EdgeF1B0])
						V3_SCALAR_PRODUCT(conv,Area[EdgeF1B0],sp) sp = ABS(sp); if (sp>0.0) V3_SCALE(1.0/sp,GIP[EdgeF1B0])

						V3_SUBTRACT(CornerPoints[FrontCorner[1]],CornerPoints[FrontCorner[0]],a)
						V3_SUBTRACT(EdgeMidPoints[BackEdge],CornerPoints[FrontCorner[0]],b)
						V3_VECTOR_PRODUCT(a,b,e)
						if (MirrorAtPlane(GIP[EdgeF1B0],CornerPoints[FrontCorner[0]],e,GIP[EdgeF1B1])) return (1);
						if (CornerOfEdge[EdgeF1B1][0]==FrontCorner[1])	V3_COPY(d,Area[EdgeF1B1])
						else											V3_COPY(c,Area[EdgeF1B1])
							
						/* the local upwind intergration points */
						V3_COPY(TRefCoord[FrontCorner[0]],LIP[FrontEdge])
					}
					break;
				case 2:
					/* FrontEdgeMidPoint in 'outflow[1]' ,symmetric in 'inflow' */
					if (param[1]<=0.75)
					{
						/* no coupling of B0 and B1 */
						
						/* entries for F0B1 and F1B1 */
						V3_SUBTRACT(EdgeMidPoints[BackEdge],EdgeMidPoints[EdgeF0B1],a)
						V3_SUBTRACT(SideMidPoints[inflow[0]],SideMidPoints[outflow[1]],b)
						V3_VECTOR_PRODUCT(a,b,c)
						V3_SCALE(0.5,c) V3_COPY(c,d)
						V3_SCALAR_PRODUCT(conv,c,sp)
						if (sp>0.0) V3_SCALE(-1.0,d)
						else 		V3_SCALE(-1.0,c)
						
						V3_QUA_CM(SideMidPoints[inflow[0]],SideMidPoints[outflow[1]],EdgeMidPoints[EdgeF0B1],EdgeMidPoints[BackEdge],GIP[EdgeF0B1])
						if (CornerOfEdge[EdgeF0B1][0]==FrontCorner[0])	V3_COPY(d,Area[EdgeF0B1])
						else											V3_COPY(c,Area[EdgeF0B1])
						
						V3_QUA_CM(SideMidPoints[inflow[1]],SideMidPoints[outflow[1]],EdgeMidPoints[EdgeF1B1],EdgeMidPoints[BackEdge],GIP[EdgeF1B1])
						if (CornerOfEdge[EdgeF1B1][0]==FrontCorner[1])	V3_COPY(d,Area[EdgeF1B1])
						else											V3_COPY(c,Area[EdgeF1B1])
						
						/* entries for F0B0 and F1B0 */
						V3_SUBTRACT(EdgeMidPoints[BackEdge],EdgeMidPoints[EdgeF0B0],a)
						V3_SUBTRACT(SideMidPoints[inflow[0]],SideMidPoints[outflow[0]],b)
						V3_VECTOR_PRODUCT(a,b,c)
						V3_SCALE(0.5,c) V3_COPY(c,d)
						V3_SCALAR_PRODUCT(conv,c,sp)
						if (sp>0.0) V3_SCALE(-1.0,d)
						else 		V3_SCALE(-1.0,c)
						
						V3_QUA_CM(SideMidPoints[inflow[0]],SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF0B0],EdgeMidPoints[BackEdge],GIP[EdgeF0B0])
						if (CornerOfEdge[EdgeF0B0][0]==FrontCorner[0])	V3_COPY(d,Area[EdgeF0B0])
						else											V3_COPY(c,Area[EdgeF0B0])
						
						V3_QUA_CM(SideMidPoints[inflow[1]],SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF1B0],EdgeMidPoints[BackEdge],GIP[EdgeF1B0])
						if (CornerOfEdge[EdgeF1B0][0]==FrontCorner[1])	V3_COPY(d,Area[EdgeF1B0])
						else											V3_COPY(c,Area[EdgeF1B0])
					}
					else
					{
						/* coupling of B0 and B1 */
						
						/* entries for B0B1, F0B0 and F1B0 */
						if (FindCrossParam3D(EdgeMidPoints[BackEdge],SideMidPoints[inflow[0]],SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF0B0],conv,param)) return (1);
						V3_SUBTRACT(EdgeMidPoints[BackEdge],SideMidPoints[inflow[0]],a)
						V3_SUBTRACT(SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF0B0],b)
						V3_VECTOR_PRODUCT(a,b,c) V3_SCALAR_PRODUCT(conv,c,sp) if (sp<0.0) {V3_SCALE(-1.0,c); sp = -sp;}

						V3_LINCOMB(0.5,SideMidPoints[outflow[0]],0.5,EdgeMidPoints[BackEdge],GIP[BackEdge])
						V3_COPY(c,Area[BackEdge])
						if (CornerOfEdge[BackEdge][0]==BackCorner[0]) alpha = -param[0]*param[1];
						else										  alpha = param[0]*param[1];
						V3_SCALE(alpha,Area[BackEdge])
							
						V3_LINCOMB(1.0-param[0],EdgeMidPoints[BackEdge],param[0],SideMidPoints[inflow[0]],a)
						V3_TRI_CM(SideMidPoints[inflow[0]],EdgeMidPoints[EdgeF0B0],a,GIP[EdgeF0B0])
						V3_COPY(c,Area[EdgeF0B0])
						if (CornerOfEdge[EdgeF0B0][0]==FrontCorner[0])  alpha = 0.5*(1.0-param[0])*(param[1]-1.0);
						else											alpha = 0.5*(1.0-param[0])*(1.0-param[1]);
						V3_SCALE(alpha,Area[EdgeF0B0])
							
						V3_SUBTRACT(CornerPoints[BackCorner[1]],CornerPoints[BackCorner[0]],a)
						V3_SUBTRACT(EdgeMidPoints[FrontEdge],CornerPoints[BackCorner[0]],b)
						V3_VECTOR_PRODUCT(a,b,d)
						if (MirrorAtPlane(GIP[EdgeF0B0],CornerPoints[BackCorner[0]],d,GIP[EdgeF1B0])) return (1);
						V3_COPY(c,Area[EdgeF1B0])
						if (CornerOfEdge[EdgeF1B0][0]==FrontCorner[1])  alpha = 0.5*(1.0-param[0])*(param[1]-1.0);
						else											alpha = 0.5*(1.0-param[0])*(1.0-param[1]);
						V3_SCALE(alpha,Area[EdgeF1B0])
						
						/* entries for F1B1 and F0B1 */
						V3_COPY(c,d) alpha=0.5*(1.0-param[0])*param[1]; V3_SCALE(alpha,d)

						V3_LINCOMB(1.0-param[0],EdgeMidPoints[BackEdge],param[0],SideMidPoints[inflow[0]],a)
						V3_TRI_CM(SideMidPoints[inflow[0]],EdgeMidPoints[BackEdge],a,GIP[EdgeF0B1])
						V3_SCALE(sp*alpha,GIP[EdgeF0B1])
					
						V3_SUBTRACT(SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF0B1],a)
						V3_SUBTRACT(SideMidPoints[inflow[0]],SideMidPoints[outflow[1]],b)
						V3_VECTOR_PRODUCT(a,b,c) V3_SCALAR_PRODUCT(conv,c,sp) if (sp<0.0) {V3_SCALE(-0.5,c);sp=-sp;} else V3_SCALE(0.5,c)
						V3_ADD(c,d,c) V3_COPY(c,d) V3_SCALE(-1.0,d)
							
						V3_QUA_CM(SideMidPoints[outflow[1]],SideMidPoints[inflow[0]],EdgeMidPoints[EdgeF0B1],SideMidPoints[outflow[0]],cm)
						V3_LINCOMB(1.0,GIP[EdgeF0B1],0.5*sp,cm,GIP[EdgeF0B1])
						if (CornerOfEdge[EdgeF0B1][0]==FrontCorner[0])	V3_COPY(d,Area[EdgeF0B1])
						else											V3_COPY(c,Area[EdgeF0B1])
							
						V3_SUBTRACT(CornerPoints[BackCorner[1]],CornerPoints[BackCorner[0]],a)
						V3_SUBTRACT(EdgeMidPoints[FrontEdge],CornerPoints[BackCorner[0]],b)
						V3_VECTOR_PRODUCT(a,b,e)
						if (MirrorAtPlane(GIP[EdgeF0B1],CornerPoints[BackCorner[0]],e,GIP[EdgeF1B1])) return (1);
						if (CornerOfEdge[EdgeF1B1][0]==FrontCorner[1])	V3_COPY(d,Area[EdgeF1B1])
						else											V3_COPY(c,Area[EdgeF1B1])
							
						/* the local upwind intergration points */
						V3_COPY(TRefCoord[BackCorner[1]],LIP[BackEdge])
					}
					break;
				case 3:
					/* FrontEdgeMidPoint in 'inflow[0]' and 'outflow[1]' */
					if (FindCrossParam3D(SideMidPoints[outflow[0]],EdgeMidPoints[FrontEdge],EdgeMidPoints[BackEdge],SideMidPoints[inflow[0]],conv,param)) return (1);
					if (param[0]>=0.0 )
					{
						/* no coupling of B0 and B1 */
						V3_LINCOMB(1.0-param[0],SideMidPoints[outflow[0]],param[0],EdgeMidPoints[FrontEdge],a)

						/* coupling F0-B1 */
						V3_SUBTRACT(a,EdgeMidPoints[EdgeF0B1],b)
						V3_SUBTRACT(SideMidPoints[inflow[0]],EdgeMidPoints[FrontEdge],c)
						V3_VECTOR_PRODUCT(b,c,d)
						V3_SCALAR_PRODUCT(conv,d,sp) sp = ABS(sp);
						V3_QUA_CM(a,EdgeMidPoints[EdgeF0B1],SideMidPoints[inflow[0]],EdgeMidPoints[FrontEdge],GIP[EdgeF0B1])
						V3_SCALE(0.5*sp,GIP[EdgeF0B1])
						V3_SUBTRACT(EdgeMidPoints[EdgeF0B1],EdgeMidPoints[FrontEdge],b)
						V3_SUBTRACT(SideMidPoints[outflow[1]],EdgeMidPoints[FrontEdge],c)
						V3_VECTOR_PRODUCT(b,c,e)
						V3_SCALAR_PRODUCT(conv,e,sp) sp = ABS(sp);
						V3_TRI_CM(EdgeMidPoints[EdgeF0B1],EdgeMidPoints[FrontEdge],SideMidPoints[outflow[1]],cm)
						V3_LINCOMB(1.0,GIP[EdgeF0B1],0.5*sp,cm,GIP[EdgeF0B1])
						V3_ADD(d,e,d)
						V3_SCALAR_PRODUCT(conv,d,sp) if (sp<0.0) V3_SCALE(-1.0,d)
						if (CornerOfEdge[EdgeF0B1][0]==FrontCorner[0])	V3_SCALE(-1.0,d)
						V3_SCALE(0.5,d) V3_COPY(d,Area[EdgeF0B1])
						V3_SCALAR_PRODUCT(conv,d,sp) sp = ABS(sp); if (sp>0.0) V3_SCALE(1.0/sp,GIP[EdgeF0B1])

						/* coupling F1-B1 */
						V3_SUBTRACT(a,SideMidPoints[inflow[1]],b)
						V3_SUBTRACT(EdgeMidPoints[FrontEdge],EdgeMidPoints[BackEdge],c)
						V3_VECTOR_PRODUCT(b,c,d)
						V3_SCALAR_PRODUCT(conv,d,sp) sp = ABS(sp);
						V3_QUA_CM(a,SideMidPoints[inflow[1]],EdgeMidPoints[FrontEdge],EdgeMidPoints[BackEdge],GIP[EdgeF1B1])
						V3_SCALE(0.5*sp,GIP[EdgeF1B1])
						V3_SUBTRACT(EdgeMidPoints[EdgeF1B1],EdgeMidPoints[FrontEdge],b)
						V3_SUBTRACT(SideMidPoints[inflow[1]],SideMidPoints[outflow[1]],c)
						V3_VECTOR_PRODUCT(b,c,e)
						V3_SCALAR_PRODUCT(conv,e,sp) sp = ABS(sp);
						V3_QUA_CM(SideMidPoints[inflow[1]],SideMidPoints[outflow[1]],EdgeMidPoints[EdgeF1B1],EdgeMidPoints[FrontEdge],cm)
						V3_LINCOMB(1.0,GIP[EdgeF1B1],0.5*sp,cm,GIP[EdgeF1B1])
						V3_ADD(d,e,d)
						V3_SCALAR_PRODUCT(conv,d,sp) if (sp<0.0) V3_SCALE(-1.0,d)
						if (CornerOfEdge[EdgeF1B1][0]==FrontCorner[1])	V3_SCALE(-1.0,d)
						V3_SCALE(0.5,d) V3_COPY(d,Area[EdgeF1B1])
						V3_SCALAR_PRODUCT(conv,d,sp) sp = ABS(sp); if (sp>0.0) V3_SCALE(1.0/sp,GIP[EdgeF1B1])

						/* coupling F1-B0 */
						V3_SUBTRACT(a,EdgeMidPoints[EdgeF1B0],b)
						V3_SUBTRACT(EdgeMidPoints[BackEdge],SideMidPoints[outflow[0]],c)
						V3_VECTOR_PRODUCT(b,c,d)
						V3_SCALAR_PRODUCT(conv,d,sp) sp = ABS(sp);
						V3_QUA_CM(a,EdgeMidPoints[EdgeF1B0],EdgeMidPoints[BackEdge],SideMidPoints[outflow[0]],GIP[EdgeF1B0])
						V3_SCALE(0.5*sp,GIP[EdgeF1B0])
						V3_SUBTRACT(SideMidPoints[inflow[1]],EdgeMidPoints[BackEdge],b)
						V3_SUBTRACT(EdgeMidPoints[EdgeF1B0],EdgeMidPoints[BackEdge],c)
						V3_VECTOR_PRODUCT(b,c,e)
						V3_SCALAR_PRODUCT(conv,e,sp) sp = ABS(sp);
						V3_TRI_CM(SideMidPoints[inflow[1]],EdgeMidPoints[BackEdge],EdgeMidPoints[EdgeF1B0],cm)
						V3_LINCOMB(1.0,GIP[EdgeF1B0],0.5*sp,cm,GIP[EdgeF1B0])
						V3_ADD(d,e,d)
						V3_SCALAR_PRODUCT(conv,d,sp) if (sp<0.0) V3_SCALE(-1.0,d)
						if (CornerOfEdge[EdgeF1B0][0]==FrontCorner[1])	V3_SCALE(-1.0,d)
						V3_SCALE(0.5,d) V3_COPY(d,Area[EdgeF1B0])
						V3_SCALAR_PRODUCT(conv,d,sp) sp = ABS(sp); if (sp>0.0) V3_SCALE(1.0/sp,GIP[EdgeF1B0])

						/* coupling F0-B0 */
						V3_QUA_CM(a,EdgeMidPoints[EdgeF0B0],SideMidPoints[outflow[0]],SideMidPoints[inflow[0]],GIP[EdgeF0B0])
						V3_SUBTRACT(a,EdgeMidPoints[EdgeF0B0],b)
						V3_SUBTRACT(SideMidPoints[outflow[0]],SideMidPoints[inflow[0]],c)
						V3_VECTOR_PRODUCT(b,c,d)
						V3_SCALAR_PRODUCT(conv,d,sp) if (sp<0.0) V3_SCALE(-1.0,d)
						if (CornerOfEdge[EdgeF0B0][0]==FrontCorner[0])	V3_SCALE(-1.0,d)
						V3_SCALE(0.5,d) V3_COPY(d,Area[EdgeF0B0])
					}
					else
					{
						/* coupling of B0 and B1 */
						
						/* entries for B0B1(the half) and F0B0 */
						if (FindCrossParam3D(EdgeMidPoints[BackEdge],SideMidPoints[inflow[0]],SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF0B0],conv,param)) return (1);
						V3_SUBTRACT(EdgeMidPoints[BackEdge],SideMidPoints[inflow[0]],a)
						V3_SUBTRACT(SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF0B0],b)
						V3_VECTOR_PRODUCT(a,b,c) V3_SCALAR_PRODUCT(conv,c,sp) if (sp<0.0) {V3_SCALE(-1.0,c); sp = -sp;}
							
						V3_LINCOMB(1.0-param[0],EdgeMidPoints[BackEdge],param[0],SideMidPoints[inflow[0]],GIP[BackEdge])
						V3_COPY(c,Area[BackEdge])
						if (CornerOfEdge[BackEdge][0]==BackCorner[0]) alpha = -0.5*param[0]*param[1];
						else										  alpha = 0.5*param[0]*param[1];
						V3_SCALE(alpha,Area[BackEdge])
							
						V3_TRI_CM(GIP[BackEdge],SideMidPoints[inflow[0]],EdgeMidPoints[EdgeF0B0],GIP[EdgeF0B0])
						V3_COPY(c,Area[EdgeF0B0])
						if (CornerOfEdge[EdgeF0B0][0]==FrontCorner[0])  alpha = 0.5*(1.0-param[0])*(param[1]-1.0);
						else											alpha = 0.5*(1.0-param[0])*(1.0-param[1]);
						V3_SCALE(alpha,Area[EdgeF0B0])

						/* entry for F0B1 */
						V3_LINCOMB(1.0-param[0],EdgeMidPoints[BackEdge],param[0],SideMidPoints[inflow[0]],e)
						V3_SUBTRACT(EdgeMidPoints[EdgeF0B1],e,a)
						V3_SUBTRACT(EdgeMidPoints[FrontEdge],SideMidPoints[inflow[0]],b)
						V3_VECTOR_PRODUCT(a,b,c)
						V3_SCALAR_PRODUCT(conv,c,sp) sp = ABS(sp);
						V3_QUA_CM(EdgeMidPoints[EdgeF0B1],e,EdgeMidPoints[FrontEdge],SideMidPoints[inflow[0]],GIP[EdgeF0B1])
						V3_SCALE(0.5*sp,GIP[EdgeF0B1])
						V3_SUBTRACT(EdgeMidPoints[EdgeF0B1],EdgeMidPoints[FrontEdge],a)
						V3_SUBTRACT(SideMidPoints[outflow[1]],EdgeMidPoints[FrontEdge],b)
						V3_VECTOR_PRODUCT(a,b,d) V3_ADD(c,d,c) 
						V3_SCALAR_PRODUCT(conv,d,sp) sp = ABS(sp);
						V3_TRI_CM(EdgeMidPoints[EdgeF0B1],EdgeMidPoints[FrontEdge],SideMidPoints[outflow[1]],cm)
						V3_LINCOMB(1.0,GIP[EdgeF0B1],0.5*sp,cm,GIP[EdgeF0B1])
						V3_SUBTRACT(EdgeMidPoints[FrontEdge],e,a)
						V3_SUBTRACT(SideMidPoints[outflow[0]],e,b)
						V3_VECTOR_PRODUCT(a,b,d) V3_ADD(c,d,c) 
						V3_SCALAR_PRODUCT(conv,d,sp) sp = ABS(sp);
						V3_TRI_CM(EdgeMidPoints[FrontEdge],e,SideMidPoints[outflow[0]],cm)
						V3_LINCOMB(1.0,GIP[EdgeF0B1],0.5*sp,cm,GIP[EdgeF0B1])
    					V3_SCALAR_PRODUCT(conv,c,sp) if (sp<0.0) V3_SCALE(-0.5,c) else V3_SCALE(0.5,c)
						if (CornerOfEdge[EdgeF0B1][0]==FrontCorner[0])  V3_SCALE(-1.0,c)
						V3_COPY(c,Area[EdgeF0B1])
						V3_SCALAR_PRODUCT(conv,c,sp) sp = ABS(sp); if (sp>0.0) V3_SCALE(1.0/sp,GIP[EdgeF0B1])
							
						/* entries for B0B1 (the half) and F1B0 */
						if (FindCrossParam3D(EdgeMidPoints[BackEdge],SideMidPoints[inflow[1]],SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF1B0],conv,param)) return (1);
						V3_SUBTRACT(EdgeMidPoints[BackEdge],SideMidPoints[inflow[1]],a)
						V3_SUBTRACT(SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF1B0],b)
						V3_VECTOR_PRODUCT(a,b,c) V3_SCALAR_PRODUCT(conv,c,sp) if (sp<0.0) V3_SCALE(-1.0,c)
						V3_COPY(c,d)
						
						V3_LINCOMB(1.0-param[0],EdgeMidPoints[BackEdge],param[0],SideMidPoints[inflow[1]],e)
						V3_ADD(GIP[BackEdge],e,GIP[BackEdge])
						if (CornerOfEdge[BackEdge][0]==BackCorner[0]) alpha = -0.5*param[0]*param[1];
						else										  alpha = 0.5*param[0]*param[1];
						V3_SCALE(alpha,d) V3_ADD(d,Area[BackEdge],Area[BackEdge])
						V3_ADD(EdgeMidPoints[BackEdge],GIP[BackEdge],GIP[BackEdge]) V3_ADD(SideMidPoints[outflow[0]],GIP[BackEdge],GIP[BackEdge])
						V3_SCALE(0.25,GIP[BackEdge])
							
						V3_TRI_CM(e,SideMidPoints[inflow[1]],EdgeMidPoints[EdgeF1B0],GIP[EdgeF1B0])
						V3_COPY(c,Area[EdgeF1B0])
						if (CornerOfEdge[EdgeF1B0][0]==FrontCorner[1])  alpha = 0.5*(1.0-param[0])*(param[1]-1.0);
						else											alpha = 0.5*(1.0-param[0])*(1.0-param[1]);
						V3_SCALE(alpha,Area[EdgeF1B0])

						/* entry for F1B1 */	
						V3_LINCOMB(1.0-param[0],EdgeMidPoints[BackEdge],param[0],SideMidPoints[inflow[1]],e)
						V3_SUBTRACT(EdgeMidPoints[EdgeF1B1],e,a)
						V3_SUBTRACT(SideMidPoints[inflow[1]],SideMidPoints[outflow[0]],b)
						V3_VECTOR_PRODUCT(a,b,c)
						V3_SCALAR_PRODUCT(conv,c,sp) sp = ABS(sp);
						V3_QUA_CM(EdgeMidPoints[EdgeF1B1],e,SideMidPoints[inflow[1]],SideMidPoints[outflow[0]],GIP[EdgeF1B1])
						V3_SCALE(0.5*sp,GIP[EdgeF1B1])
						V3_SUBTRACT(SideMidPoints[outflow[1]],SideMidPoints[outflow[0]],a)
						V3_SUBTRACT(EdgeMidPoints[EdgeF1B1],EdgeMidPoints[FrontEdge],b)
						V3_VECTOR_PRODUCT(a,b,d) V3_ADD(c,d,c) V3_SCALAR_PRODUCT(conv,c,sp) if (sp<0.0) V3_SCALE(-0.5,c) else V3_SCALE(0.5,c)
						if (CornerOfEdge[EdgeF1B1][0]==FrontCorner[1])  V3_SCALE(-1.0,c)
						V3_COPY(c,Area[EdgeF1B1])
						V3_SCALAR_PRODUCT(conv,d,sp) sp = ABS(sp);
						V3_QUA_CM(SideMidPoints[outflow[1]],SideMidPoints[outflow[0]],EdgeMidPoints[EdgeF1B1],EdgeMidPoints[FrontEdge],cm)
						V3_LINCOMB(1.0,GIP[EdgeF1B1],0.5*sp,cm,GIP[EdgeF1B1])
						V3_SCALAR_PRODUCT(conv,c,sp) sp = ABS(sp); if (sp>0.0) V3_SCALE(1.0/sp,GIP[EdgeF1B1])
							
						/* the local upwind intergration points */
						V3_COPY(TRefCoord[BackCorner[1]],LIP[BackEdge])
					}
					break;
				default :
					return (1);
			}
						
			/* the local upwind intergration points */
			V3_COPY(TRefCoord[BackCorner[0]],LIP[EdgeF0B0])
			V3_COPY(TRefCoord[BackCorner[0]],LIP[EdgeF1B0])
			V3_COPY(TRefCoord[BackCorner[1]],LIP[EdgeF0B1])
			V3_COPY(TRefCoord[BackCorner[1]],LIP[EdgeF1B1])
			break;
		case 1:
			/* one in-, three outflow */
			OpCorner = OppositeCorner[inflow[0]];
			for (i=0; i<3; i++)
			{
				corner[i] = CornerOfSide[inflow[0]][i];
				OpEdge[i] = EdgeWithCorners[corner[i]][OpCorner];
				GrEdge[i] = EdgeWithCorners[corner[i]][CornerOfSide[inflow[0]][(i+1)%3]];
				side[i]   = SideOfCorners[corner[i]][CornerOfSide[inflow[0]][(i+1)%3]][OpCorner];
			}
			
			/* calculate the three quadrilaterals */
			V3_SUBTRACT(CornerPoints[corner[1]],CornerPoints[corner[0]],a)	/* a = P1 - P0 	*/
			V3_SUBTRACT(CornerPoints[corner[2]],CornerPoints[corner[0]],b)	/* b = P2 - P0 	*/
			V3_VECTOR_PRODUCT(a,b,c)
			V3_SCALE(0.02777777777777777777,c)
			for (i=0; i<3; i++)
			{
				V3_SUBTRACT(CornerPoints[corner[(i+2)%3]],CornerPoints[corner[(i+1)%3]],a)
				V3_SUBTRACT(CornerPoints[OpCorner],CornerPoints[corner[(i+1)%3]],b)
				V3_VECTOR_PRODUCT(a,b,d)
				V3_SCALE(0.08333333333333333333,d)
				V3_ADD(c,d,Area[OpEdge[i]])
				
				V3_SCALAR_PRODUCT(conv,Area[OpEdge[i]],sp) sp = ABS(sp);
				V3_QUA_CM(EdgeMidPoints[OpEdge[i]],SideMidPoints[side[i]],SideMidPoints[inflow[0]],SideMidPoints[side[(i+2)%3]],GIP[OpEdge[i]])
					/*V3_SCALE(sp,GIP[OpEdge[i]])*/
			}
			
			/* calculate triangles from crosspoints */
			for (i=0; i<3; i++)
			{
				if (FindCrossParam3D(SideMidPoints[side[i]],EdgeMidPoints[OpEdge[i]],SideMidPoints[inflow[0]],EdgeMidPoints[GrEdge[i]],conv,param)) return (1);
				if (param[0]>0.0)
				{
					/* subtract first triangle to corrisponding quadrilateral */
					V3_SUBTRACT(EdgeMidPoints[OpEdge[i]],SideMidPoints[inflow[0]],a)
					V3_SUBTRACT(SideMidPoints[side[i]],SideMidPoints[inflow[0]],b)
					V3_VECTOR_PRODUCT(a,b,c)
					if (param[0]>=1.0)
						V3_SCALE(0.5,c)
					else
						V3_SCALE(0.5*param[0],c)
					V3_SUBTRACT(Area[OpEdge[i]],c,Area[OpEdge[i]])
					V3_ADD(Area[OpEdge[(i+1)%3]],c,Area[OpEdge[(i+1)%3]])

					V3_LINCOMB(1.0-param[0],SideMidPoints[side[i]],param[0],EdgeMidPoints[OpEdge[i]],e)/*
					V3_TRI_CM(e,SideMidPoints[side[i]],SideMidPoints[inflow[0]],cm)
					V3_SCALAR_PRODUCT(conv,c,sp) sp = ABS(sp);
					V3_SCALE(sp,cm)
					V3_SUBTRACT(GIP[OpEdge[i]],cm,GIP[OpEdge[i]])
					V3_ADD(GIP[OpEdge[(i+1)%3]],cm,GIP[OpEdge[(i+1)%3]])*/		/* still to normalize !! */
							
					/* connections of P0 and P1 */
					V3_SUBTRACT(EdgeMidPoints[OpEdge[i]],EdgeMidPoints[GrEdge[i]],a)
					V3_SUBTRACT(SideMidPoints[side[i]],EdgeMidPoints[GrEdge[i]],b)
					V3_VECTOR_PRODUCT(a,b,Area[GrEdge[i]])
					if (param[0]>=1.0)
						V3_SCALE(0.5,Area[GrEdge[i]])
					else
						V3_SCALE(0.5*param[0],Area[GrEdge[i]])
					V3_TRI_CM(e,SideMidPoints[side[i]],EdgeMidPoints[GrEdge[i]],GIP[GrEdge[i]])					

					First = (CornerOfEdge[GrEdge[i]][0]==corner[i]);
					V3_SCALAR_PRODUCT(conv,Area[GrEdge[i]],sp)
					if ((!inverted && !First && sp<0.0) ||	
					    (!inverted && First && sp>0.0) ||	
				 	    (inverted && First && sp<0.0) ||	
				   	    (inverted && !First && sp>0.0))	
						V3_SCALE(-1.0,Area[GrEdge[i]])
					V3_COPY(TRefCoord[corner[(i+1)%3]],LIP[GrEdge[i]])
				}
				else if (param[0]<0.0)
				{
					if (FindCrossParam3D(SideMidPoints[side[i]],EdgeMidPoints[OpEdge[(i+1)%3]],SideMidPoints[inflow[0]],EdgeMidPoints[GrEdge[i]],conv,param)) return (1);
					if (param[0]>0.0)
					{
						V3_SUBTRACT(SideMidPoints[side[i]],SideMidPoints[inflow[0]],a)
						V3_SUBTRACT(EdgeMidPoints[OpEdge[(i+1)%3]],SideMidPoints[inflow[0]],b)
						V3_VECTOR_PRODUCT(a,b,c)
						if (param[0]>=1.0)
							V3_SCALE(0.5,c)
						else
							V3_SCALE(0.5*param[0],c)
						V3_ADD(Area[OpEdge[i]],c,Area[OpEdge[i]])
						V3_SUBTRACT(Area[OpEdge[(i+1)%3]],c,Area[OpEdge[(i+1)%3]])

						V3_LINCOMB(1.0-param[0],SideMidPoints[side[i]],param[0],EdgeMidPoints[OpEdge[(i+1)%3]],e)/*
						V3_TRI_CM(e,SideMidPoints[side[i]],SideMidPoints[inflow[0]],cm)
						V3_SCALAR_PRODUCT(conv,c,sp) sp = ABS(sp);
						V3_SCALE(sp,cm)
						V3_ADD(GIP[OpEdge[i]],cm,GIP[OpEdge[i]])
						V3_SUBTRACT(GIP[OpEdge[(i+1)%3]],cm,GIP[OpEdge[(i+1)%3]])*/	/* still to normalize !! */
							
						/* connections of P0 and P1 */
						V3_SUBTRACT(EdgeMidPoints[OpEdge[(i+1)%3]],EdgeMidPoints[GrEdge[i]],a)
						V3_SUBTRACT(SideMidPoints[side[i]],EdgeMidPoints[GrEdge[i]],b)
						V3_VECTOR_PRODUCT(a,b,Area[GrEdge[i]])
						if (param[0]>=1.0)
							V3_SCALE(0.5,Area[GrEdge[i]])
						else
							V3_SCALE(0.5*param[0],Area[GrEdge[i]])
						V3_TRI_CM(e,SideMidPoints[side[i]],EdgeMidPoints[GrEdge[i]],GIP[GrEdge[i]])					

						First = (CornerOfEdge[GrEdge[i]][0]==corner[i]);
						V3_SCALAR_PRODUCT(conv,Area[GrEdge[i]],sp)
						if ((!inverted && !First && sp>0.0) ||	
					 	   (!inverted && First && sp<0.0) ||	
				 	  	   (inverted && First && sp>0.0) ||	
				   	   	   (inverted && !First && sp<0.0))	
							V3_SCALE(-1.0,Area[GrEdge[i]])
						V3_COPY(TRefCoord[corner[i]],LIP[GrEdge[i]])
					}
				}
				else
					V3_LINCOMB(0.5,SideMidPoints[side[i]],0.5,EdgeMidPoints[GrEdge[i]],GIP[GrEdge[i]])
			}

			/*for (i=0; i<3; i++)
			{
				V3_SCALAR_PRODUCT(conv,Area[OpEdge[i]],sp) sp = ABS(sp);
				if (sp>0.0)	V3_SCALE(1.0/sp,GIP[OpEdge[i]])
			}*/
			
			/* turn to the right direction and set LIPs */			
			for (i=0; i<3; i++)
			{
				/* the main connections */
				First = (CornerOfEdge[OpEdge[i]][1]==OpCorner);
				V3_SCALAR_PRODUCT(conv,Area[OpEdge[i]],sp)
				if ((!inverted && !First && sp>0.0) ||	
				    (!inverted && First && sp<0.0) ||	
				    (inverted && First && sp>0.0) ||	
				    (inverted && !First && sp<0.0))	
					V3_SCALE(-1.0,Area[OpEdge[i]])
				
				if (inverted)
					V3_COPY(TRefCoord[OpCorner],LIP[OpEdge[i]])
				else
					V3_COPY(TRefCoord[corner[i]],LIP[OpEdge[i]])
			}
			break;			
		default:
			return (1);		
	}
	
	/* local IP's from global IP's */
	V3_SUBTRACT(CornerPoints[1],CornerPoints[0],M)
	V3_SUBTRACT(CornerPoints[2],CornerPoints[0],M+3)
	V3_SUBTRACT(CornerPoints[3],CornerPoints[0],M+6);
	if (M3_Invert(Inv,M)) return (1);
	for (i=0; i<6; i++)
	{ 
		V3_SUBTRACT(GIP[i],CornerPoints[0],a)
		M3_TIMES_V3(Inv,a,LIP[i])
	}
	
	return (0);
}
#endif

/****************************************************************************/
/*																			*/
/* Function:  FV_TetInfo_for_conv											*/
/*																			*/
/* Purpose:   calc gradient of shape function i 							*/
/*																			*/
/* Input:	  int i: corner number [0..3]									*/
/*			  COORD **Corners: list of ptrs to phys corner vectors			*/
/*																			*/
/* Output:	  COORD **theGradient: list of ptrs to gradients				*/
/*																			*/
/* Return:	  INT 0: ok 													*/
/*				  1: error													*/
/*																			*/
/****************************************************************************/

#ifdef __THREEDIM__						 
INT FV_TetInfo_for_conv (ELEMENT *theElement, const COORD **CornerPoints, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], COORD_VECTOR GIP[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM], COORD_VECTOR conv)
{
	COORD sp, spn, spz, alpha1, alpha2;
	COORD_VECTOR a, b, c, normal;
	COORD_VECTOR EdgeMidPoints[6], SideMidPoints[4];
	INT i, j, help, noutflow, ninflow, outflow[4], inflow[4], edge, inverted, side;
	
	/* reset areas and integrationpoints */
	for (i=0; i<EDGES_OF_ELEM(theElement); i++)
	{
		V3_CLEAR(Area[i])
		V3_CLEAR(GIP[i])
		V3_CLEAR(LUIP[i])
	}
	
	/* edge mid points */
	for (i=0; i<6; i++)	V3_LINCOMB(0.5,CornerPoints[CORNER_OF_EDGE(theElement,i,0)],0.5,CornerPoints[CORNER_OF_EDGE(theElement,i,1)],EdgeMidPoints[i])
			
	/* side mid points */
	for (i=0; i<SIDES_OF_ELEM(theElement); i++)
	{
		V3_ADD(CornerPoints[CORNER_OF_SIDE(theElement,i,0)],CornerPoints[CORNER_OF_SIDE(theElement,i,1)],a)
		V3_ADD(CornerPoints[CORNER_OF_SIDE(theElement,i,2)],a,a)
		V3_SCALE(0.33333333333333333,SideMidPoints[i])
	}

	/* in/outflow sides */
	noutflow = ninflow = 0;
	for (i=0; i<4; i++)
	{
		V3_SUBTRACT(CornerPoints[CORNER_OF_SIDE(theElement,i,1)],CornerPoints[CORNER_OF_SIDE(theElement,i,0)],a)
		V3_SUBTRACT(CornerPoints[CORNER_OF_SIDE(theElement,i,2)],CornerPoints[CORNER_OF_SIDE(theElement,i,0)],b)
		V3_VECTOR_PRODUCT(a,b,normal)
		V3_SUBTRACT(CornerPoints[CORNER_OPP_TO_SIDE(theElement,i)],CornerPoints[CORNER_OF_SIDE(theElement,i,0)],a)
		V3_SCALAR_PRODUCT(a,normal,sp)
		if (sp>0.0) V3_SCALE(-1.0,normal)
		V3_SCALAR_PRODUCT(conv,normal,sp)
		if (sp>=0.0) 	outflow[noutflow++] = i;
		if (sp<0.0)		inflow[ninflow++] = i;
	}

	/* change if inflow > outflow */
	inverted=0;
	if (ninflow>noutflow)
	{
		inverted=1;
		for (i=0; i<ninflow; i++)	{help = inflow[i]; inflow[i] = outflow[i]; outflow[i] = help;}
		help = ninflow; ninflow = noutflow; noutflow = help;
	}
	if ( ninflow+noutflow!=4 || ninflow<1 || ninflow>2) 
		return (1);
		
	/* handle different cases */
	switch (ninflow)
	{
		case 2:
		  break;
		case 1:
			/* one in-, three outflow */
			for (i=0; i<3; i++)	
			{
				/* P0 = corners of element-side[i]  */
				edge = EDGE_OF_CORNER(theElement,CORNER_OF_SIDE(theElement,inflow[0],i),CORNER_OPP_TO_SIDE(theElement,inflow[0]));								/* edge: P0 ~ P3	*/
				V3_SUBTRACT(SideMidPoints[inflow[0]],EdgeMidPoints[edge],a)													/* P012 - P3		*/
				V3_SUBTRACT(CornerPoints[CORNER_OF_SIDE(theElement,inflow[0],(i+2)%3)],CornerPoints[CORNER_OF_SIDE(theElement,inflow[0],(i+1)%3)],b)/* P2 - P1 			*/
				V3_VECTOR_PRODUCT(a,b,Area[edge])
				V3_SCALE(0.1666666666,Area[edge])

				/* P2 = corners of element-side[i]  */
				/* calc alpha1 */
				edge = EDGE_OF_CORNER(theElement,CORNER_OF_SIDE(theElement,inflow[0],(i+1)%3),CORNER_OPP_TO_SIDE(theElement,inflow[0]));						/* edge: P0 ~ P3 	*/
				V3_SUBTRACT(EdgeMidPoints[edge],CornerPoints[CORNER_OF_SIDE(theElement,inflow[0],(i+2)%3)],a)							/* P03 - P1 		*/
				V3_VECTOR_PRODUCT(conv,a,c)
				edge = EDGE_OF_CORNER(theElement,CORNER_OF_SIDE(theElement,inflow[0],(i+1)%3),CORNER_OF_SIDE(theElement,inflow[0],(i+2)%3));					/* edge: P0 ~ P1 */
				V3_SUBTRACT(EdgeMidPoints[edge],CornerPoints[CORNER_OF_SIDE(theElement,inflow[0],i)],b)									/* P01 - P2 */
				V3_SCALAR_PRODUCT(b,c,spn)
				if (spn==0.0) continue;
				V3_SUBTRACT(CornerPoints[CORNER_OPP_TO_SIDE(theElement,inflow[0])],CornerPoints[CORNER_OF_SIDE(theElement,inflow[0],i)],b)	   /* P3 - P2 */
				V3_SCALAR_PRODUCT(b,c,spz)
				alpha1 = spz/spn;
				assert (alpha1>=0.0);
				
				/* calc alpha2 */
				edge = EDGE_OF_CORNER(theElement,CORNER_OF_SIDE(theElement,inflow[0],(i+2)%3),CORNER_OPP_TO_SIDE(theElement,inflow[0]));						/* edge: P1 ~ P3 	*/
				V3_SUBTRACT(EdgeMidPoints[edge],CornerPoints[CORNER_OF_SIDE(theElement,inflow[0],(i+1)%3)],a)							/* P13 - P0 		*/
				V3_VECTOR_PRODUCT(conv,a,c)
				edge = EDGE_OF_CORNER(theElement,CORNER_OF_SIDE(theElement,inflow[0],(i+1)%3),CORNER_OF_SIDE(theElement,inflow[0],(i+2)%3));					/* edge: P0 ~ P1 	*/
				V3_SUBTRACT(EdgeMidPoints[edge],CornerPoints[CORNER_OF_SIDE(theElement,inflow[0],i)],b)									/* P01 - P2 		*/
				V3_SCALAR_PRODUCT(b,c,spn)
				assert (spn!=0.0);
				V3_SUBTRACT(CornerPoints[CORNER_OPP_TO_SIDE(theElement,inflow[0])],CornerPoints[CORNER_OF_SIDE(theElement,inflow[0],i)],b)				/* P3 - P2 			*/
				V3_SCALAR_PRODUCT(b,c,spz)
				alpha2 = spz/spn;
				assert (alpha2>=0.0);
				
				/* take the right triangle */
				V3_SUBTRACT(CornerPoints[CORNER_OF_SIDE(theElement,inflow[0],i)],CornerPoints[CORNER_OPP_TO_SIDE(theElement,inflow[0])],a)				/* P2 - P3 			*/
				edge = EDGE_OF_CORNER(theElement,CORNER_OF_SIDE(theElement,inflow[0],(i+1)%3),CORNER_OF_SIDE(theElement,inflow[0],(i+2)%3));					/* edge: P0 ~ P1 	*/
				V3_SUBTRACT(EdgeMidPoints[edge],CornerPoints[CORNER_OF_SIDE(theElement,inflow[0],i)],b)									/* P01 - P2 		*/
				V3_VECTOR_PRODUCT(a,b,c)
				  if (alpha1<alpha2)		j=1;
				else if (alpha2<alpha1) j=0;
				else
				{
					side = SIDE_OPP_TO_CORNER(theElement,CORNER_OF_SIDE(theElement,inflow[0],i));	   				/* side 013 		*/
					V3_SUBTRACT(SideMidPoints[side],CornerPoints[CORNER_OF_SIDE(theElement,inflow[0],(i+2)%3)],a)						/* P013 - P1 		*/
					V3_EUKLIDNORM(a,alpha1)
					V3_SUBTRACT(SideMidPoints[side],CornerPoints[CORNER_OF_SIDE(theElement,inflow[0],(i+1)%3)],a)						/* P013 - P0 		*/
					V3_EUKLIDNORM(a,alpha2)
					if (alpha1<alpha2)		j=1;
					else if (alpha2<alpha1) j=0;
					else continue;
				}
				if (j)
				{
					V3_SCALE(alpha1/18.0,c)
					edge = EDGE_OF_CORNER(theElement,CORNER_OF_SIDE(theElement,inflow[0],i),CORNER_OPP_TO_SIDE(theElement,inflow[0]));							/* edge: P0 ~ P3 	*/
					V3_SUBTRACT(Area[edge],c,Area[edge])			
					edge = EDGE_OF_CORNER(theElement,CORNER_OF_SIDE(theElement,inflow[0],(i+2)%3),CORNER_OPP_TO_SIDE(theElement,inflow[0]));					/* edge: P1 ~ P3 	*/
					V3_ADD(Area[edge],c,Area[edge])			
				}
				else
				{
					V3_SCALE(alpha2/18.0,c)
					edge = EDGE_OF_CORNER(theElement,CORNER_OF_SIDE(theElement,inflow[0],i),CORNER_OPP_TO_SIDE(theElement,inflow[0]));							/* edge: P0 ~ P3 	*/
					V3_ADD(Area[edge],c,Area[edge])			
					edge = EDGE_OF_CORNER(theElement,CORNER_OF_SIDE(theElement,inflow[0],(i+2)%3),CORNER_OPP_TO_SIDE(theElement,inflow[0]));					/* edge: P1 ~ P3 	*/
					V3_SUBTRACT(Area[edge],c,Area[edge])			
				 }
			}
			/* orientate the areas */
			for (i=0; i<6; i++)	
			{
				V3_SUBTRACT(CornerPoints[CORNER_OF_EDGE(theElement,i,1)],CornerPoints[CORNER_OF_EDGE(theElement,i,0)],a)
				V3_SCALAR_PRODUCT(a,Area[i],sp)
				if (sp<0.0)
					V3_SCALE(-1.0,Area[i])
			}
			
			/* set LUIP */
			for (i=0; i<3; i++)	
			{			
				edge = EDGE_OF_CORNER(theElement,CORNER_OF_SIDE(theElement,inflow[0],i),CORNER_OPP_TO_SIDE(theElement,inflow[0]));
				if (inverted)	V3_COPY(TRefCoord[CORNER_OPP_TO_SIDE(theElement,inflow[0])],LUIP[edge])
				else			V3_COPY(TRefCoord[CORNER_OF_SIDE(theElement,inflow[0],i)],LUIP[edge])
			 }
			break;
		default:
			return (1);		
    }
	
	return (0);
}
#endif

/****************************************************************************/
/*D
   Side_TetInfo	- Calculate subsurfaces and integration points on elementside of tetrahedron 

   SYNOPSIS:
   INT Side_TetInfo (COORD **theCorners, INT side, COORD_VECTOR Area, 
   COORD_VECTOR GIP[3]);

   PARAMETERS:
.  theCorners - list of pointers to phys corner vectors
.  side - side of tetrahedron
.  Area - area of subsurface on an element side
.  Gip[3] - surface integration points

   DESCRIPTION:
   The tetrahedrons are subdevided into four subcontrol volumes, defined by the
   edge midpoints, the center of mass of the sides and the center of mass
   of the volume for the discretization. In natural one elementside decomposes
   in three subsurfaces. The area and the global integration points of these
   subsurfaces are calculated in this function.

   RETURN VALUE:
   INT
.n    0 if ok
.n    1 if error occured.
D*/
/****************************************************************************/

#ifdef __THREEDIM__						 
INT Side_TetInfo (COORD **theCorners, INT side, COORD_VECTOR Area, COORD_VECTOR GIP[3])
{
   	ELEMENT e;
	COORD_VECTOR a,b,c;
	COORD scalarprd;
	INT i0, i1, i2, i3;
	
	SETTAG(&e,4);
	i0 = CORNER_OF_SIDE(&e,side,0);
	i1 = CORNER_OF_SIDE(&e,side,1);
	i2 = CORNER_OF_SIDE(&e,side,2);
	i3 = CORNER_OPP_TO_SIDE(&e,side);

	/* get Area */
	V3_SUBTRACT(theCorners[i1],theCorners[i0],a)
	V3_SUBTRACT(theCorners[i2],theCorners[i0],b)	
	V3_SUBTRACT(theCorners[i3],theCorners[i0],c)
	V3_VECTOR_PRODUCT(a,b,Area)
	V3_SCALAR_PRODUCT(c,Area,scalarprd)
	if (scalarprd > 0.0)
		V3_SCALE(-1.0/6.0,Area)
	else
		V3_SCALE( 1.0/6.0,Area)
	
	/* get three surface integration points */
	V3_LINCOMB(14.0/24.0,theCorners[i0],5.0/24.0,theCorners[i1],a)
	V3_LINCOMB(5.0/24.0,theCorners[i2],1.0,a,GIP[0])
	
	V3_LINCOMB(14.0/24.0,theCorners[i1],5.0/24.0,theCorners[i2],a)
	V3_LINCOMB(5.0/24.0,theCorners[i0],1.0,a,GIP[1])
	
	V3_LINCOMB(14.0/24.0,theCorners[i2],5.0/24.0,theCorners[i0],a)
	V3_LINCOMB(5.0/24.0,theCorners[i1],1.0,a,GIP[2])
	
	return (0);
}
#endif

/****************************************************************************/
/*D
   GSUIP - Calculate upwind integration point

   SYNOPSIS:
   INT GSUIP (COORD **theCorners, COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], 
   DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);

   PARAMETERS:
.  theCorners - coordinates of the element corners
.  LIP - integration point
.  conv - velocity
.  LUIP - upwind integration point

   DESCRIPTION:
   This function calculates the upwind integration points LUIP on the surface
   using the `standard upwind discretization`, where the LUIP is defined downstream
   at the node on the intersection of the local flow vector `conv`, defined in the
   integration point IP, with the element boundary.

   RETURN VALUE:
   INT
.n    0 if ok
.n    1 if error occured.	
D*/			
/****************************************************************************/

#ifdef __THREEDIM__						 
INT GetSkewedUIP (const COORD_VECTOR *theCorners, const COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], const DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM])
{
	COORD_VECTOR lconv;
	COORD alpha;
	COORD M[9],I[9];
	INT flags, i;

	V3_SUBTRACT(theCorners[1],theCorners[0],M)
	V3_SUBTRACT(theCorners[2],theCorners[0],M+3)
	V3_SUBTRACT(theCorners[3],theCorners[0],M+6)
	if (M3_Invert(I,M)) return (1);
	for (i=0; i<6; i++)
	{
		M3_TIMES_V3 (I,conv[i],lconv)
		flags = (ABS(lconv[0])<SMALL_C);
		flags |= ((ABS(lconv[1])<SMALL_C)<<1);
		flags |= ((ABS(lconv[2])<SMALL_C)<<2);
		switch(flags)
		{
			case 0:
				if (lconv[0]>0.0)
				{
					alpha = LIP[i][0]/lconv[0];
					LUIP[i][1]=LIP[i][1]-alpha*lconv[1];
					LUIP[i][2]=LIP[i][2]-alpha*lconv[2];
					if (LUIP[i][1]>=0.0 && LUIP[i][2]>=0.0 && LUIP[i][1]+LUIP[i][2]<=1.0)
					{
						LUIP[i][0]=0.0;
						break;
					}
				}
				if (lconv[1]>0.0)
				{
					alpha = LIP[i][1]/lconv[1];
					LUIP[i][0]=LIP[i][0]-alpha*lconv[0];
					LUIP[i][2]=LIP[i][2]-alpha*lconv[2];
					if (LUIP[i][0]>=0.0 && LUIP[i][2]>=0.0 && LUIP[i][0]+LUIP[i][2]<=1.0)
					{
						LUIP[i][1]=0.0;
						break;
					}
				}
				if (lconv[2]>0.0)
				{
					alpha = LIP[i][2]/lconv[2];
					LUIP[i][0]=LIP[i][0]-alpha*lconv[0];
					LUIP[i][1]=LIP[i][1]-alpha*lconv[1];
					if (LUIP[i][0]>=0.0 && LUIP[i][1]>=0.0 && LUIP[i][0]+LUIP[i][1]<=1.0)
					{
						LUIP[i][2]=0.0;
						break;
					}
				}
				alpha = (LIP[i][0]+LIP[i][1]+LIP[i][2]-1.0)/(lconv[0]+lconv[1]+lconv[2]);
				LUIP[i][0]=LIP[i][0]-alpha*lconv[0];
				LUIP[i][1]=LIP[i][1]-alpha*lconv[1];
				LUIP[i][2]=LIP[i][2]-alpha*lconv[2];
				break;
			case 1:
				LUIP[i][0]=LIP[i][0];
				if (lconv[1]>0.0)
				{
					alpha = LIP[i][1]/lconv[1];
					LUIP[i][2]=LIP[i][2]-alpha*lconv[2];
					if (LUIP[i][2]>=0.0 && LUIP[i][0]+LUIP[i][2]<=1.0)
					{
						LUIP[i][1]=0.0;
						break;
					}
				}
				if (lconv[2]>0.0)
				{
					alpha = LIP[i][2]/lconv[2];
					LUIP[i][1]=LIP[i][1]-alpha*lconv[1];
					if (LUIP[i][1]>=0.0 && LUIP[i][0]+LUIP[i][1]<=1.0)
					{
						LUIP[i][2]=0.0;
						break;
					}
				}	
				alpha = (LIP[i][0]+LIP[i][1]+LIP[i][2]-1.0)/(lconv[1]+lconv[2]);
				LUIP[i][1]=LIP[i][1]-alpha*lconv[1];
				LUIP[i][2]=LIP[i][2]-alpha*lconv[2];
				break;
			case 2:
				LUIP[i][1]=LIP[i][1];
				if (lconv[0]>0.0)
				{
					alpha = LIP[i][0]/lconv[0];
					LUIP[i][2]=LIP[i][2]-alpha*lconv[2];
					if (LUIP[i][2]>=0.0 && LUIP[i][1]+LUIP[i][2]<=1.0)
					{
						LUIP[i][0]=0.0;
						break;
					}
				}
				if (lconv[2]>0.0)
				{
					alpha = LIP[i][2]/lconv[2];
					LUIP[i][0]=LIP[i][0]-alpha*lconv[0];
					if (LUIP[i][0]>=0.0 && LUIP[i][0]+LUIP[i][1]<=1.0)
					{
						LUIP[i][2]=0.0;
						break;
					}
				}	
				alpha = (LIP[i][0]+LIP[i][1]+LIP[i][2]-1.0)/(lconv[0]+lconv[2]);
				LUIP[i][0]=LIP[i][0]-alpha*lconv[0];
				LUIP[i][2]=LIP[i][2]-alpha*lconv[2];
				break;
			case 3:
				LUIP[i][0]=LIP[i][0];
				LUIP[i][1]=LIP[i][1];
				if (lconv[2]>0.0)
					LUIP[i][2]=0.0;
				else
					LUIP[i][2]=1.0-LIP[i][0]-LIP[i][1];
				break;
			case 4:
				LUIP[i][2]=LIP[i][2];
				if (lconv[0]>0.0)
				{
					alpha = LIP[i][0]/lconv[0];
					LUIP[i][1]=LIP[i][1]-alpha*lconv[1];
					if (LUIP[i][1]>=0.0 && LUIP[i][1]+LUIP[i][2]<=1.0)
					{
						LUIP[i][0]=0.0;
						break;
					}
				}
				if (lconv[1]>0.0)
				{
					alpha = LIP[i][1]/lconv[1];
					LUIP[i][0]=LIP[i][0]-alpha*lconv[0];
					if (LUIP[i][0]>=0.0 && LUIP[i][0]+LUIP[i][2]<=1.0)
					{
						LUIP[i][1]=0.0;
						break;
					}
				}
				alpha = (LIP[i][0]+LIP[i][1]+LIP[i][2]-1.0)/(lconv[0]+lconv[1]);
				LUIP[i][0]=LIP[i][0]-alpha*lconv[0];
				LUIP[i][1]=LIP[i][1]-alpha*lconv[1];
				break;
			case 5:
				LUIP[i][0]=LIP[i][0];
				if (lconv[1]>0.0)
					LUIP[i][1]=0.0;
				else
					LUIP[i][1]=1.0-LIP[i][0]-LIP[i][2];
				LUIP[i][2]=LIP[i][2];
				break;
			case 6:
				if (lconv[0]>0.0)
					LUIP[i][0]=0.0;
				else
					LUIP[i][0]=1.0-LIP[i][1]-LIP[i][2];
				LUIP[i][1]=LIP[i][1];
				LUIP[i][2]=LIP[i][2];
				break;
			case 7:
				LUIP[i][0]=LIP[i][0];
				LUIP[i][1]=LIP[i][1];
				LUIP[i][2]=LIP[i][2];
				break;
		}
	/*	if ((LUIP[i][0]<-SMALL_C) || (LUIP[i][1]<-SMALL_C) || (LUIP[i][2]<-SMALL_C) || (LUIP[i][0]+LUIP[i][1]+LUIP[i][2]>1.0+SMALL_C))
			return (1);*/
	}
	
	return(0);
}
#endif

/****************************************************************************/
/*D
   GFUIP - Calculate upwind integration point

   SYNOPSIS:
   INT GFUIP (COORD **theCorners, COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], 
   DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);

   PARAMETERS:
.  theCorners - coordinates of the element corners
.  LIP - integration point
.  conv - velocity in the integration point
.  LUIP - upwind integration point

   DESCRIPTION:
   This function calculates the upwind integration points LUIP using the
   `fast upwind discretization`. There the LUIP is defined at the corner with
   the maximum distance upstream to the intersection of the local flow vector
   `conv`, defined in the integration point IP, with the element boundary.
   
   RETURN VALUE:
   INT
.n    0 if ok
.n    1 if error occured.													
D*/
/****************************************************************************/

#ifdef __THREEDIM__						 
INT GFUIP (const COORD **theCorners, const COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM])
{
	COORD_VECTOR lconv;
	COORD sp, min;
	COORD M[9],I[9];
	INT j, i, ip[MAX_CORNERS_OF_ELEM], n;

	V3_SUBTRACT(theCorners[1],theCorners[0],M)
	V3_SUBTRACT(theCorners[2],theCorners[0],M+3)
	V3_SUBTRACT(theCorners[3],theCorners[0],M+6)
	if (M3_Invert(I,M)) return (1);
	for (i=0; i<6; i++)
	{
		M3_TIMES_V3(I,conv[i],lconv)
		min = MAX_C; n = 0;
		for (j=0; j<4; j++)
		{
			V3_SCALAR_PRODUCT(lconv,TRefCoord[j],sp)
			if (min==sp)
				ip[n++] = j;
			if (min>sp)
			{
				n = 0;
				min =sp;
				ip[n++] = j;
			}
		}
		assert(n>0);
		V3_CLEAR(LUIP[i])
		for (j=0; j<n; j++)
			V3_ADD(TRefCoord[ip[j]],LUIP[i],LUIP[i])
		V3_SCALE(1.0/(COORD)n,LUIP[i])
	}
	
	return(0);
}
#endif

/****************************************************************************/
/*D
   GCUIP - Calculate upwind integration point

   SYNOPSIS:
   INT GCUIP (COORD **theCorners, COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], 
   DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);

   PARAMETERS:
.  theCorners - coordinates of the element corners
.  LIP - integration point
.  conv - velocity
.  LUIP - upwind integration point


   DESCRIPTION:
   This function calculates the upwind integration point LUIP using the
   `corner upwind discretization`, where the LUIP is defined downstream at the
   corner next to the intersection-point of the local flow vector `conv`, defined
   in the integration point IP, with the element boundary.

   RETURN VALUE:
   INT
.n    0 if ok
.n    1 if error occured.
D*/
/****************************************************************************/

#ifdef __THREEDIM__						 
INT GCUIP (const COORD **theCorners, const COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM])
{
	COORD_VECTOR a, lconv, SUIP;
	COORD alpha, sp, min;
	COORD M[9],I[9];
	INT flags, i, j, k; 

	V3_SUBTRACT(theCorners[1],theCorners[0],M)
	V3_SUBTRACT(theCorners[2],theCorners[0],M+3)
	V3_SUBTRACT(theCorners[3],theCorners[0],M+6)
	if (M3_Invert(I,M)) return (1);
	for (i=0; i<6; i++)
	{
		V3_SUBTRACT(conv[i],theCorners[0],a)
		M3_TIMES_V3 (I,a,lconv)
		flags = (ABS(lconv[0])<SMALL_C);
		flags |= ((ABS(lconv[1])<SMALL_C)<<1);
		flags |= ((ABS(lconv[2])<SMALL_C)<<2);
		switch(flags)
		{
			case 0:
				if (lconv[0]>0.0)
				{
					alpha = LIP[i][0]/lconv[0];
					SUIP[1]=LIP[i][1]-alpha*lconv[1];
					SUIP[2]=LIP[i][2]-alpha*lconv[2];
					if (SUIP[1]>=0.0 && SUIP[2]>=0.0 && SUIP[1]+SUIP[2]<=1.0)
					{
						SUIP[0]=0.0;
						break;
					}
				}
				if (lconv[1]>0.0)
				{
					alpha = LIP[i][1]/lconv[1];
					SUIP[0]=LIP[i][0]-alpha*lconv[0];
					SUIP[2]=LIP[i][2]-alpha*lconv[2];
					if (SUIP[0]>=0.0 && SUIP[2]>=0.0 && SUIP[0]+SUIP[2]<=1.0)
					{
						SUIP[1]=0.0;
						break;
					}
				}
				if (lconv[2]>0.0)
				{
					alpha = LIP[i][2]/lconv[2];
					SUIP[0]=LIP[i][0]-alpha*lconv[0];
					SUIP[1]=LIP[i][1]-alpha*lconv[1];
					if (SUIP[0]>=0.0 && SUIP[1]>=0.0 && SUIP[0]+SUIP[1]<=1.0)
					{
						SUIP[2]=0.0;
						break;
					}
				}
				alpha = (LIP[i][0]+LIP[i][1]+LIP[i][2]-1.0)/(lconv[0]+lconv[1]+lconv[2]);
				SUIP[0]=LIP[i][0]-alpha*lconv[0];
				SUIP[1]=LIP[i][1]-alpha*lconv[1];
				SUIP[2]=LIP[i][2]-alpha*lconv[2];
				break;
			case 1:
				SUIP[0]=LIP[i][0];
				if (lconv[1]>0.0)
				{
					alpha = LIP[i][1]/lconv[1];
					SUIP[2]=LIP[i][2]-alpha*lconv[2];
					if (SUIP[2]>=0.0 && SUIP[0]+SUIP[2]<=1.0)
					{
						SUIP[1]=0.0;
						break;
					}
				}
				if (lconv[2]>0.0)
				{
					alpha = LIP[i][2]/lconv[2];
					SUIP[1]=LIP[i][1]-alpha*lconv[1];
					if (SUIP[1]>=0.0 && SUIP[0]+SUIP[1]<=1.0)
					{
						SUIP[2]=0.0;
						break;
					}
				}	
				alpha = (LIP[i][0]+LIP[i][1]+LIP[i][2]-1.0)/(lconv[1]+lconv[2]);
				SUIP[1]=LIP[i][1]-alpha*lconv[1];
				SUIP[2]=LIP[i][2]-alpha*lconv[2];
				break;
			case 2:
				SUIP[1]=LIP[i][1];
				if (lconv[0]>0.0)
				{
					alpha = LIP[i][0]/lconv[0];
					SUIP[2]=LIP[i][2]-alpha*lconv[2];
					if (SUIP[2]>=0.0 && SUIP[1]+SUIP[2]<=1.0)
					{
						SUIP[0]=0.0;
						break;
					}
				}
				if (lconv[2]>0.0)
				{
					alpha = LIP[i][2]/lconv[2];
					SUIP[0]=LIP[i][0]-alpha*lconv[0];
					if (SUIP[0]>=0.0 && SUIP[0]+SUIP[1]<=1.0)
					{
						SUIP[2]=0.0;
						break;
					}
				}	
				alpha = (LIP[i][0]+LIP[i][1]+LIP[i][2]-1.0)/(lconv[0]+lconv[2]);
				SUIP[0]=LIP[i][0]-alpha*lconv[0];
				SUIP[2]=LIP[i][2]-alpha*lconv[2];
				break;
			case 3:
				SUIP[0]=LIP[i][0];
				SUIP[1]=LIP[i][1];
				if (lconv[2]>0.0)
					SUIP[2]=0.0;
				else
					SUIP[2]=1.0-LIP[i][0]-LIP[i][1];
				break;
			case 4:
				SUIP[2]=LIP[i][2];
				if (lconv[0]>0.0)
				{
					alpha = LIP[i][0]/lconv[0];
					SUIP[1]=LIP[i][1]-alpha*lconv[1];
					if (SUIP[1]>=0.0 && SUIP[1]+SUIP[2]<=1.0)
					{
						SUIP[0]=0.0;
						break;
					}
				}
				if (lconv[1]>0.0)
				{
					alpha = LIP[i][1]/lconv[1];
					SUIP[0]=LIP[i][0]-alpha*lconv[0];
					if (SUIP[0]>=0.0 && SUIP[0]+SUIP[2]<=1.0)
					{
						SUIP[1]=0.0;
						break;
					}
				}
				alpha = (LIP[i][0]+LIP[i][1]+LIP[i][2]-1.0)/(lconv[0]+lconv[1]);
				SUIP[0]=LIP[i][0]-alpha*lconv[0];
				SUIP[1]=LIP[i][1]-alpha*lconv[1];
				break;
			case 5:
				SUIP[0]=LIP[i][0];
				if (lconv[1]>0.0)
					SUIP[1]=0.0;
				else
					SUIP[1]=1.0-LIP[i][0]-LIP[i][2];
				SUIP[2]=LIP[i][2];
				break;
			case 6:
				if (lconv[0]>0.0)
					SUIP[0]=0.0;
				else
					SUIP[0]=1.0-LIP[i][1]-LIP[i][2];
				SUIP[1]=LIP[i][1];
				SUIP[2]=LIP[i][2];
				break;
			case 7:
				SUIP[0]=LIP[i][0];
				SUIP[1]=LIP[i][1];
				SUIP[2]=LIP[i][2];
				break;
		}
		
		/* now find corner */
		min = MAX_C;
		for (j=0; j<4; j++)
		{
			V3_EUKLIDNORM_OF_DIFF(SUIP,TRefCoord[j],sp)
			if (min>sp)
			{
				k=j;
				min=sp;
			}
		}
		assert(k>=0 && k<4);
		LUIP[i][0] = TRefCoord[k][0];
		LUIP[i][1] = TRefCoord[k][1];
		LUIP[i][2] = TRefCoord[k][2];
	}
	
	return(0);
}
#endif

/****************************************************************************/
/*D
   COPYIP - Copy integration points

   SYNOPSIS:
   INT COPYIP (COORD **theCorners, COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], 
   DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);

   PARAMETERS:
.  theCorners - coordinates of the element corners
.  LIP - integration point
.  conv - velocity 
.  LUIP - upwind integration point
  
   DESCRIPTION:
   This function copies integration points on the upwind integration points.

   RETURN VALUE:
   INT
.n    0 if ok
.n    1 if error occured.
D*/
/****************************************************************************/

#ifdef __THREEDIM__						 
INT COPYIP (const COORD **theCorners, const COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM])
{
	INT i;
	
	for (i=0; i<6; i++)
		V3_COPY(LIP[i],LUIP[i])
	
	return (0);
}
#endif
