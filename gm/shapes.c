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

#define SMALL_DET      1e-50
#define SMALL_DIFF     1e-4
#define MAX_ITER       20 

/* some useful abbreviations */
#define Xi  ((DOUBLE)ip_local[0])
#define Eta ((DOUBLE)ip_local[1])
#define Mu  ((DOUBLE)ip_local[2])

#define J11 J[0][0]
#define J12 J[0][1]
#define J13 J[0][2]
#define J21 J[1][0]
#define J22 J[1][1]
#define J23 J[1][2]
#define J31 J[2][0]
#define J32 J[2][1]
#define J33 J[2][2]

#define X0  co_global[0][0]
#define X1  co_global[1][0]
#define X2  co_global[2][0]
#define X3  co_global[3][0]
#define X4  co_global[4][0]
#define X5  co_global[5][0]
#define X6  co_global[6][0]
#define X7  co_global[7][0]

#define Y0  co_global[0][1]
#define Y1  co_global[1][1]
#define Y2  co_global[2][1]
#define Y3  co_global[3][1]
#define Y4  co_global[4][1]
#define Y5  co_global[5][1]
#define Y6  co_global[6][1]
#define Y7  co_global[7][1]

#define Z0  co_global[0][2]
#define Z1  co_global[1][2]
#define Z2  co_global[2][2]
#define Z3  co_global[3][2]
#define Z4  co_global[4][2]
#define Z5  co_global[5][2]
#define Z6  co_global[6][2]
#define Z7  co_global[7][2]

#define U0  nodal_values[0]
#define U1  nodal_values[1]
#define U2  nodal_values[2]
#define U3  nodal_values[3]
#define U4  nodal_values[4]
#define U5  nodal_values[5]
#define U6  nodal_values[6]
#define U7  nodal_values[7]

#define SQ(x) (x)*(x)

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

/* local midpoints */
static COORD_VECTOR_2D LMP_Triangle 		= {0.33333333333, 0.33333333333};
static COORD_VECTOR_2D LMP_Quadrilateral	= {0.5, 0.5};
static COORD_VECTOR_3D LMP_Tetrahedron		= {0.25, 0.25, 0.25};
static COORD_VECTOR_3D LMP_Pyramid   		= {0.5, 0.5, 0.33333333333333333};
static COORD_VECTOR_3D LMP_Prism    		= {0.33333333333,0.3333333333,0.5};
static COORD_VECTOR_3D LMP_Hexahedron		= {0.5, 0.5, 0.5};

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*D
   LocalCornerCoordinates - Return corners of reference element

   SYNOPSIS:
   INT LocalCornerCoordinates (INT dim, INT tag, INT corner, DOUBLE *result);

   PARAMETERS:
.  dim - space dimension, 2D can elements can be used in 3D
.  tag - identifier for element type
.  corner - corner for which the coordinates should be returned
.  result - vector to store result

   DESCRIPTION:
   This routine returns the coordinates of the corners of the reference element.
   The function is currently implemented in two dimensions for triangles (tag=3)
   and quadrilaterals (tag=4).

   RETURN VALUES:
   0 when o.k.

   1 if an error occured.
D*/
/****************************************************************************/

INT LocalCornerCoordinates (INT dim, INT tag, INT corner, DOUBLE *result)
{
	if (dim==2) 
	{
		if (tag==TRIANGLE)
			switch (corner) {
				case 0: result[0] = 0.0; result[1] = 0.0; return(0);
				case 1: result[0] = 1.0; result[1] = 0.0; return(0);
				case 2: result[0] = 0.0; result[1] = 1.0; return(0);
				default: return(1);
			}
		if (tag==QUADRILATERAL)
			switch (corner) {
				case 0: result[0] = 0.0; result[1] = 0.0; return(0);
				case 1: result[0] = 1.0; result[1] = 0.0; return(0);
				case 2: result[0] = 1.0; result[1] = 1.0; return(0);
				case 3: result[0] = 0.0; result[1] = 1.0; return(0);
				default: return(1);
			}
	}
	if (dim==3)
	{
		switch (tag) {
			case TETRAHEDRON:
				switch (corner) {
					case 0: result[0] = 0.0; result[1] = 0.0; result[2] = 0.0; return(0);
					case 1: result[0] = 1.0; result[1] = 0.0; result[2] = 0.0; return(0);
					case 2: result[0] = 0.0; result[1] = 1.0; result[2] = 0.0; return(0);
					case 3: result[0] = 0.0; result[1] = 0.0; result[2] = 1.0; return(0);
					default: return(1);
				}
			case PYRAMID:
				switch (corner) {
					case 0: result[0] = 0.0; result[1] = 0.0; result[2] = 0.0; return(0);
					case 1: result[0] = 1.0; result[1] = 0.0; result[2] = 0.0; return(0);
					case 2: result[0] = 1.0; result[1] = 1.0; result[2] = 0.0; return(0);
					case 3: result[0] = 0.0; result[1] = 1.0; result[2] = 0.0; return(0);
					case 4: result[0] = 0.0; result[1] = 0.0; result[2] = 1.0; return(0);
					default: return(1);
				}
			case PRISM:
				switch (corner) {
					case 0: result[0] = 0.0; result[1] = 0.0; result[2] = 0.0; return(0);
					case 1: result[0] = 1.0; result[1] = 0.0; result[2] = 0.0; return(0);
					case 2: result[0] = 0.0; result[1] = 1.0; result[2] = 0.0; return(0);
					case 4: result[0] = 0.0; result[1] = 0.0; result[2] = 1.0; return(0);
					case 5: result[0] = 1.0; result[1] = 0.0; result[2] = 1.0; return(0);
					case 6: result[0] = 0.0; result[1] = 1.0; result[2] = 1.0; return(0);
					default: return(1);
				}
			case HEXAHEDRON:
				switch (corner) {
					case 0: result[0] = 0.0; result[1] = 0.0; result[2] = 0.0; return(0);
					case 1: result[0] = 1.0; result[1] = 0.0; result[2] = 0.0; return(0);
					case 2: result[0] = 1.0; result[1] = 1.0; result[2] = 0.0; return(0);
					case 3: result[0] = 0.0; result[1] = 1.0; result[2] = 0.0; return(0);
					case 4: result[0] = 0.0; result[1] = 0.0; result[2] = 1.0; return(0);
					case 5: result[0] = 1.0; result[1] = 0.0; result[2] = 1.0; return(0);
					case 6: result[0] = 1.0; result[1] = 1.0; result[2] = 1.0; return(0);
					case 7: result[0] = 0.0; result[1] = 1.0; result[2] = 1.0; return(0);
					default: return(1);
				}
			default: return(1);
		}
	}				

	return (1);

}

/****************************************************************************/
/*D
   InterpolateFEFunction - Interpolate a finite element function 

   SYNOPSIS:
   INT InterpolateFEFunction (INT dim, INT tag, DOUBLE ip_local[DIM],
        DOUBLE nodal_values[MAX_CORNERS_OF_ELEM], DOUBLE *result);

   PARAMETERS:
.  dim - space dimension, 2D can elements can be used in 3D
.  tag - identifier for element type
.  ip_local - point in coordinates of the reference element where FE function should be evaluated
.  nodal_values - values of finite element function in corners
.  result - where to store result

   DESCRIPTION:
   This function interpolates a finite element function given by its nodal
   values. Currently standard linear elements for triangles (tag=3) and
   bilinear elements for quadrilaterals (tag=4) are implemented.

   RETURN VALUE:
   0 when o.k.

   1 if an error occured.
D*/
/****************************************************************************/

INT InterpolateFEFunction (INT dim, INT tag, DOUBLE ip_local[DIM],
        DOUBLE nodal_values[MAX_CORNERS_OF_ELEM], DOUBLE *result)
{
	if (dim==2)
	{
		if (tag==TRIANGLE) {
			*result = U0 + Xi*(U1-U0) + Eta*(U2-U0);
			return (0);
		}
		if (tag==QUADRILATERAL) {
			*result = U0 + Xi*(U1-U0) + Eta*(U3-U0) + Xi*Eta*(U0-U1+U2-U3);
			return (0);
		}
	}
	if ((dim==3)&&(DIM==3))
	{
		switch (tag) {
			case TETRAHEDRON:
				*result = (1-Xi-Eta-Mu) * nodal_values[0]
				        +      Xi       * nodal_values[1]
				        +      Eta      * nodal_values[2]
				        +      Mu       * nodal_values[3];
				return(0);
			case PYRAMID:
				if (Xi > Eta)
				  {
					*result = ((1.0-Xi)*(1.0-Eta)-Mu*(1.0-Eta))*nodal_values[0]
					  +(Xi*(1.0-Eta)-Mu*Eta)*nodal_values[1]
						+(Xi*Eta+Mu*Eta)*nodal_values[2]
						  +((1.0-Xi)*Eta-Mu*Eta)*nodal_values[3]
							+Mu*nodal_values[4];
    				return (0);
				  }
				else
				  {
					*result = ((1.0-Xi)*(1.0-Eta)-Mu*(1.0-Xi))*nodal_values[0]
					  +(Xi*(1.0-Eta)-Mu*Xi)*nodal_values[1]
						+(Xi*Eta+Mu*Xi)*nodal_values[2]
						  +((1.0-Xi)*Eta-Mu*Xi)*nodal_values[3]
							+Mu*nodal_values[4];
    				return (0);
				  }
			case PRISM:
				*result = (1-Xi-Eta)*(1-Mu) * nodal_values[0]
				        + Xi*(1-Mu)         * nodal_values[1]
				        + Eta*(1-Mu)        * nodal_values[2]
				        + (1-Xi-Eta)*Mu     * nodal_values[3]
				        + Xi*Mu             * nodal_values[4]
				        + Eta*Mu            * nodal_values[5];
				return(0);
			case HEXAHEDRON:
				*result = (1-Xi)*(1-Eta)*(1-Mu) * nodal_values[0]
				        +  (Xi) *(1-Eta)*(1-Mu) * nodal_values[1]
				        +  (Xi) * (Eta) *(1-Mu) * nodal_values[2]
				        + (1-Xi)* (Eta) *(1-Mu) * nodal_values[3]
						+ (1-Xi)*(1-Eta)* (Mu)  * nodal_values[4]
				        +  (Xi) *(1-Eta)* (Mu)  * nodal_values[5]
				        +  (Xi) * (Eta) * (Mu)  * nodal_values[6]
				        + (1-Xi)* (Eta) * (Mu)  * nodal_values[7];
				return(0);
			default: return(1);
		}
	}				
	
	return (1);
}

/****************************************************************************/
/*D
   LinearTrafo - return whether coordinate trafo is linear or not

   SYNOPSIS:
   INT LinearTrafo (INT dim, INT tag);

   PARAMETERS:
.  dim - space dimension, 2D can elements can be used in 3D
.  tag - identifier for element type

   DESCRIPTION:
   This function returns TRUE when the coordinate transformation from the
   reference element to the general element is linear. In that case the
   Jacobian need only be evaluated once per element.

   RETURN VALUES:
   TRUE when transformation is linear

   FALSE when transformation is nonlinear
D*/
/****************************************************************************/

INT LinearTrafo (INT dim, INT tag)
{
	if (dim==2)
	{
		if (tag==TRIANGLE) return(1);
		if (tag==QUADRILATERAL) return(0);
	}
	if (dim==3)
	{
		if (tag==TETRAHEDRON) return(1);
		if (tag==PYRAMID) return(0);
		if (tag==PRISM) return(0);
		if (tag==HEXAHEDRON) return(0);
	}
		
	return(0);
}

/****************************************************************************/
/*D
   JacobianInverse - compute inverse of jacobian

   SYNOPSIS:
   INT JacobianInverse (INT dim, INT tag, DOUBLE co_global[MAX_CORNERS_OF_ELEM][DIM],
       DOUBLE ip_local[DIM], DOUBLE Jinv[DIM][DIM], DOUBLE *detJ);

   PARAMETERS:
.  dim - space dimension, 2D can elements can be used in 3D
.  tag - identifier for element type
.  co_global - global coordintes of corners of the element
.  ip_local - point in coordinates of the reference element where Jacobian should be evaluated
.  Jinv - where to store resulting inverse of Jacobian
.  detJ - where to store determinant of Jacobian

   DESCRIPTION:
   Let T be the coordinate transformation from the reference element to a
   general element.
   This function computes the determinant of the Jacobian of T and
   the inverse of the Jacobian of T at a given point in reference coordinates. 
   Currently triangles (tag=3) and
   quadrilaterals (tag=4) are implemented.

   RETURN VALUES:
   0 when o.k.

   1 if an error occured (e.g. determinant is very small)

   SEE ALSO:
   GradientFEFunction, LinearTrafo
D*/
/****************************************************************************/

INT JacobianInverse (INT dim, INT tag, DOUBLE co_global[MAX_CORNERS_OF_ELEM][DIM],
       DOUBLE ip_local[DIM], DOUBLE Jinv[DIM][DIM], DOUBLE *detJ)
{
	DOUBLE det;
	DOUBLE x[MAX_CORNERS_OF_ELEM], y[MAX_CORNERS_OF_ELEM];
	INT i;
	DOUBLE J[DIM][DIM],xi,eta;
	DOUBLE sub0123, sub0145, sub0347, sub10234567;

	if (dim==2)
	{
		for (i=0; i<tag; i++) {
			x[i] = co_global[i][0]; y[i] = co_global[i][1];
		}
		xi = ip_local[0]; eta = ip_local[1];
	
		if (tag==TRIANGLE) {
	    	J[0][0] = x[1]-x[0]; J[0][1] = y[1]-y[0];
	    	J[1][0] = x[2]-x[0]; J[1][1] = y[2]-y[0];
		}
	
		if (tag==QUADRILATERAL) {
			J[0][0] = -(1-eta)*x[0] + (1-eta)*x[1] + eta*x[2] - eta*x[3];
			J[0][1] = -(1-eta)*y[0] + (1-eta)*y[1] + eta*y[2] - eta*y[3];
			J[1][0] = -(1-xi)*x[0] - xi*x[1] + xi*x[2] + (1-xi)*x[3];
			J[1][1] = -(1-xi)*y[0] - xi*y[1] + xi*y[2] + (1-xi)*y[3];
		}
	
	    det = J[0][0]*J[1][1] - J[0][1]*J[1][0];
	 	if (fabs(det)<=1.0E-15) return (1);
	
		/* fill inverse */
	    Jinv[0][0] =  J[1][1]/det; Jinv[0][1] = -J[0][1]/det;
	    Jinv[1][0] = -J[1][0]/det; Jinv[1][1] =  J[0][0]/det;
	
		/* return determinant */
		*detJ = det;
		return(0);
	}
	if ((dim==3)&&(DIM==3))
	{
		/* fill jacobian */
		switch (tag) {
			case TETRAHEDRON:
				J11 = X1-X0;	J12 = Y1-Y0;	J13 = Z1-Z0;
				J21 = X2-X0;	J22 = Y2-Y0;	J23 = Z2-Z0;
				J31 = X3-X0;	J32 = Y3-Y0;	J33 = Z3-Z0;
				break;
			case PYRAMID:
				if (Xi > Eta)
				  {
					J11= -(1.0-Eta)*X0 +(1.0-Eta)*X1+Eta*X2    -Eta*X3;
					J21= (Mu-1.0+Xi)*X0-(Xi+Mu)*X1  +(Xi+Mu)*X2+(1.0-Xi-Mu)*X3;
					J31= -(1.0-Eta)*X0 -Eta*X1;     +Eta*X2    -Eta*X3     +X4;
					J12= -(1.0-Eta)*Y0 +(1.0-Eta)*Y1+Eta*Y2    -Eta*Y3;
					J22= (Mu-1.0+Xi)*Y0-(Xi+Mu)*Y1  +(Xi+Mu)*Y2+(1.0-Xi-Mu)*Y3;
					J32= -(1.0-Eta)*Y0 -Eta*Y1;     +Eta*Y2    -Eta*Y3     +Y4;
					J13= -(1.0-Eta)*Z0 +(1.0-Eta)*Z1+Eta*Z2    -Eta*Z3;
					J23= (Mu-1.0+Xi)*Z0-(Xi+Mu)*Z1  +(Xi+Mu)*Z2+(1.0-Xi-Mu)*Z3;
					J33= -(1.0-Eta)*Z0 -Eta*Z1;     +Eta*Z2    -Eta*Z3     +Z4;
					break;
				  }
				else
				  {
					J11=(Mu+Eta-1.0)*X0+(1.0-Eta-Mu)*X1+(Eta+Mu)*X2-(Eta+Mu)*X3;
					J21=-(1.0-Xi)*X0   -Xi*X1          +Xi*X2      +(1.0-Xi)*X3;
					J31=-(1.0-Xi)*X0   -Xi*X1          +Xi*X2      -Xi*X3   +X4;
					J12=(Mu+Eta-1.0)*Y0+(1.0-Eta-Mu)*Y1+(Eta+Mu)*Y2-(Eta+Mu)*Y3;
					J22=-(1.0-Xi)*Y0   -Xi*Y1          +Xi*Y2      +(1.0-Xi)*Y3;
					J32=-(1.0-Xi)*Y0   -Xi*Y1          +Xi*Y2      -Xi*Y3   +Y4;
					J13=(Mu+Eta-1.0)*Z0+(1.0-Eta-Mu)*Z1+(Eta+Mu)*Z2-(Eta+Mu)*Z3;
					J23=-(1.0-Xi)*Z0   -Xi*Z1          +Xi*Z2      +(1.0-Xi)*Z3;
					J33=-(1.0-Xi)*Z0   -Xi*Z1          +Xi*Z2      -Xi*Z3   +Z4;
					break;
				  }
			  case PRISM:
				J11=-(1.0-Mu)*X0+(1.0-Mu)*X1     -Mu*X3+Mu*X4;
				J21=-(1.0-Mu)*X0+(1.0-Mu)*X2     -Mu*X3+Mu*X4;
				J31=-(1.0-Xi-Eta)*X0-Xi*X1-Eta*X2+(1.0-Xi-Eta)*X3+Xi*X4+Eta*X4;
				J12=-(1.0-Mu)*Y0+(1.0-Mu)*Y1     -Mu*Y3+Mu*Y4;
				J22=-(1.0-Mu)*Y0+(1.0-Mu)*Y2     -Mu*Y3+Mu*Y4;
				J32=-(1.0-Xi-Eta)*Y0-Xi*Y1-Eta*Y2+(1.0-Xi-Eta)*Y3+Xi*Y4+Eta*Y4;
				J13=-(1.0-Mu)*Z0+(1.0-Mu)*Z1     -Mu*Z3+Mu*Z4;
				J23=-(1.0-Mu)*Z0+(1.0-Mu)*Z2     -Mu*Z3+Mu*Z4;
				J33=-(1.0-Xi-Eta)*Z0-Xi*Z1-Eta*Z2+(1.0-Xi-Eta)*Z3+Xi*Z4+Eta*Z4;
				break;
			case HEXAHEDRON:
				sub0123 = X0-X1+X2-X3;
				sub0145 = X0-X1-X4+X5;
				sub0347 = X0-X3-X4+X7;
				sub10234567 = X1-X0-X2+X3+X4-X5+X6-X7;
				J11 = X1-X0 + Eta*(sub0123) + Mu *(sub0145) + Eta*Mu*(sub10234567);
				J21 = X3-X0 + Xi *(sub0123) + Mu *(sub0347) + Xi *Mu*(sub10234567);
				J31 = X4-X0 + Xi *(sub0145) + Eta*(sub0347) + Xi*Eta*(sub10234567);

				sub0123 = Y0-Y1+Y2-Y3;
				sub0145 = Y0-Y1-Y4+Y5;
				sub0347 = Y0-Y3-Y4+Y7;
				sub10234567 = Y1-Y0-Y2+Y3+Y4-Y5+Y6-Y7;
				J12 = Y1-Y0 + Eta*(sub0123) + Mu *(sub0145) + Eta*Mu*(sub10234567);
				J22 = Y3-Y0 + Xi *(sub0123) + Mu *(sub0347) + Xi *Mu*(sub10234567);
				J32 = Y4-Y0 + Xi *(sub0145) + Eta*(sub0347) + Xi*Eta*(sub10234567);

				sub0123 = Z0-Z1+Z2-Z3;
				sub0145 = Z0-Z1-Z4+Z5;
				sub0347 = Z0-Z3-Z4+Z7;
				sub10234567 = Z1-Z0-Z2+Z3+Z4-Z5+Z6-Z7;
				J13 = Z1-Z0 + Eta*(sub0123) + Mu *(sub0145) + Eta*Mu*(sub10234567);
				J23 = Z3-Z0 + Xi *(sub0123) + Mu *(sub0347) + Xi *Mu*(sub10234567);
				J33 = Z4-Z0 + Xi *(sub0145) + Eta*(sub0347) + Xi*Eta*(sub10234567);
				break;
			default: return(1);
		}
		
		/* compute determinant */
		det = J11*J22*J33 + J12*J23*J31 + J13*J21*J32
		    - J13*J22*J31 - J11*J23*J32 - J12*J21*J33;
	 	if (fabs(det)<=1.0E-15) return (1);
		*detJ = det;

		/* invert jacobian */
		Jinv[0][0] =  (J22*J33-J23*J32)/det; Jinv[0][1] = -(J12*J33-J13*J32)/det; Jinv[0][2] =  (J12*J23-J13*J22)/det;	    	
		Jinv[1][0] = -(J21*J33-J23*J31)/det; Jinv[1][1] =  (J11*J33-J13*J31)/det; Jinv[1][2] = -(J11*J23-J13*J21)/det;	    	
		Jinv[2][0] =  (J21*J32-J22*J31)/det; Jinv[2][1] = -(J11*J32-J12*J31)/det; Jinv[2][2] =  (J11*J22-J12*J21)/det;	    	

		return(0);
	}

	return(0);
}

/****************************************************************************/
/*D
   GradientFEFunction - compute gradient of finite element function 

   SYNOPSIS:
   INT GradientFEFunction (INT dim, INT tag, DOUBLE ip_local[DIM], 
       DOUBLE Jinv[DIM][DIM], DOUBLE nodal_values[MAX_CORNERS_OF_ELEM],
       DOUBLE result[DIM]);

   PARAMETERS:
.  dim - space dimension, 2D can elements can be used in 3D
.  tag - identifier for element type
.  ip_local - point in coordinates of the reference element where gradient should be evaluated
.  Jinv - inverse of Jacobian
.  nodal_values - values of finite element function in corners
.  result - where to store result

   DESCRIPTION:
   This function computes the gradient of a finite element function in a point
   given in reference coordinates.

   RETURN VALUES:
   0 when o.k.

   1 if an error occured.

   SEE ALSO:
   JacobianInverse
D*/
/****************************************************************************/

INT GradientFEFunction (INT dim, INT tag, DOUBLE ip_local[DIM], 
       DOUBLE Jinv[DIM][DIM], DOUBLE nodal_values[MAX_CORNERS_OF_ELEM],
       DOUBLE result[DIM])
{
	DOUBLE GradRef[DIM],sub1,sub2,sub3,sub4;
	INT i;

	if (dim==2)
	{
		/* Gradients of basis function on reference element */
		if (tag==TRIANGLE) {
			GradRef[0] = U1-U0; 
			GradRef[1] = U2-U0; 
		}
		if (tag==QUADRILATERAL) {
			sub1 = U0-U1+U2-U3;
			GradRef[0] = U1-U0 + Eta*sub1; 
			GradRef[1] = U3-U0 + Xi *sub1; 
		}
	
		/* transform gradient to (x,y) coordinates */
		result[0] = Jinv[0][0]*GradRef[0] + Jinv[0][1]*GradRef[1];
		result[1] = Jinv[1][0]*GradRef[0] + Jinv[1][1]*GradRef[1];
		return (0);
	}
	if ((dim==3)&&(DIM==3))
	{
		switch (tag) {
			case TETRAHEDRON:
				GradRef[0] = U1-U0; 
				GradRef[1] = U2-U0; 
				GradRef[2] = U3-U0; 
				break;
			case PYRAMID:
				sub1 = U0-U1+U2-U3;
				if (Xi > Eta)
				  {
					GradRef[0] = U1-U0 + Eta*sub1; 
					GradRef[1] = U3-U0 + (Xi + Mu)*sub1; 
					GradRef[2] = U4-U0 + Eta*sub1;
					break;
				  }
				else
				  {
					GradRef[0] = U1-U0 + (Eta + Mu)*sub1; 
					GradRef[1] = U3-U0 + Xi *sub1; 
					GradRef[2] = U4-U0 + Xi *sub1; 
					break;
				  }
			case PRISM:
				sub1 = U0-U1+U3-U4;
				sub2 = U0-U2+U3-U5;
				GradRef[0] = U1-U0 + Mu * sub1; 
				GradRef[1] = U3-U0 + Mu * sub2; 
				GradRef[2] = U4-U0 + Xi * sub1 + Eta * sub2; 
				break;
			case HEXAHEDRON:
				sub1 = U0-U1+U2-U3;
				sub2 = U0-U1-U4+U5;
				sub3 = U0-U3-U4+U7;
				sub4 = U1-U0-U2+U3+U4-U5+U6-U7;
				GradRef[0] = U1-U0 + Eta*sub1 +Mu *sub2 +Eta*Mu*sub4; 
				GradRef[1] = U3-U0 + Xi *sub1 +Mu *sub3 +Xi*Mu *sub4; 
				GradRef[2] = U4-U0 + Xi *sub2 +Eta*sub3 +Xi*Eta*sub4; 
				break;
			default: return(1);
		}

		/* transform gradient to (x,y) coordinates */
		result[0] = Jinv[0][0]*GradRef[0] + Jinv[0][1]*GradRef[1] + Jinv[0][2]*GradRef[2];
		result[1] = Jinv[1][0]*GradRef[0] + Jinv[1][1]*GradRef[1] + Jinv[1][2]*GradRef[2];
		result[2] = Jinv[2][0]*GradRef[0] + Jinv[2][1]*GradRef[1] + Jinv[2][2]*GradRef[2];
		return (0);
	}
	
	return (0);	
}

/****************************************************************************/
/*D
   SurfaceElement - compute surface element for surface integral of first kind

   SYNOPSIS:
   INT SurfaceElement (INT dim, INT nc, DOUBLE co_global[MAX_CORNERS_OF_ELEM][DIM],
       DOUBLE ip_local[DIM], DOUBLE *result)

   PARAMETERS:
.  dim - space dimension of mapped element (i.e. 2 and 3 are allowed)
.  nc - number of corners
.  co_global - global coordintes of corners of the element
.  ip_local - point in coordinates of the reference element where surfel should be evaluated

   DESCRIPTION:
   The reference elements in 1D and 2D are mapped to an arbitrary element in the
   next higher dimension, e.g. the unit quadrilateral to a general quadrilateral
   in 3D space. This function computes the surface integral for the
   surface integral of the first kind (see Bronstein 20. Auflage Kap. 3.1.12.2).

   RETURN VALUES:
   0 is ok else error

   SEE ALSO:

D*/
/****************************************************************************/

INT SurfaceElement (INT dim, INT nc, 
					DOUBLE co_global[MAX_CORNERS_OF_ELEM][DIM],
					DOUBLE ip_local[DIM], DOUBLE *result)
{
	DOUBLE E,G,F;
	
	if (dim==2)
	{
		*result = sqrt((X1-X0)*(X1-X0)+(Y1-Y0)*(Y1-Y0));
		return(0);
	}
	if (dim==3)
	{
		switch (nc)
		{
			case 3:
				E = (X1-X0)*(X1-X0) + (Y1-Y0)*(Y1-Y0) + (Z1-Z0)*(Z1-Z0);
				G = (X2-X0)*(X2-X0) + (Y2-Y0)*(Y2-Y0) + (Z2-Z0)*(Z2-Z0);
				F = (X1-X0)*(X2-X0) + (Y1-Y0)*(Y2-Y0) + (Z1-Z0)*(Z2-Z0);
				*result = sqrt(E*G-F*F);
				return(0);
				
			case 4:
				E = SQ((X0-X1)*(1-Eta)+(X2-X3)*Eta) + SQ((Y0-Y1)*(1-Eta)+(Y2-Y3)*Eta) + SQ((Z0-Z1)*(1-Eta)+(Z2-Z3)*Eta);
				G = SQ((X3-X0)*(1-Xi )+(X2-X1)*Xi ) + SQ((Y3-Y0)*(1-Xi )+(Y2-Y1)*Xi ) + SQ((Z3-Z0)*(1-Xi )+(Z2-Z1)*Xi );
				F = ((X0-X1)*(1-Eta)+(X2-X3)*Eta)*((X3-X0)*(1-Xi )+(X2-X1)*Xi )
				  + ((Y0-Y1)*(1-Eta)+(Y2-Y3)*Eta)*((Y3-Y0)*(1-Xi )+(Y2-Y1)*Xi )
				  + ((Z0-Z1)*(1-Eta)+(Z2-Z3)*Eta)*((Z3-Z0)*(1-Xi )+(Z2-Z1)*Xi );
				*result = sqrt(E*G-F*F);
				return(0);
		}
	}
	return(1);			
}

/****************************************************************************/
/*D
   GN - General Shape function for nodes

   SYNOPSIS:
   DOUBLE GN (INT n, INT i, COORD ip_local);

   PARAMETERS:
.  n - number of corners of the element
.  i - corner number (corner number [0..n-1])
.  local - local COORDinates

   DESCRIPTION:
   This function finds the value of the shape function i for the reference 
   elements at the given local coordinate.

   The shape functions fullfill
.n   Ni(node i) = 1
.n   Ni(node k) = 0, if k is not equal i.
   
   RETURN VALUE:
   DOUBLE
.n          value
D*/   
/****************************************************************************/

DOUBLE GN (INT n, INT i, COORD *ip_local)
{
    #ifdef __TWODIM__
    switch (n)
	  {
	  case 3:
		switch (i)
		  {
		  case 0 : return(1.0-Xi-Eta);
		  case 1 : return(Xi);
		  case 2 : return(Eta);
		  }
	  case 4:
		switch (i)
		  {
			case 0 : return((1.0-Xi)*(1.0-Eta));
			case 1 : return(Xi*(1.0-Eta));
			case 2 : return(Xi*Eta);
			case 3 : return((1.0-Xi)*Eta);
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
		  case 0 : return(1.0-Xi-Eta-Mu);
		  case 1 : return(Xi);
		  case 2 : return(Eta);
		  case 3 : return(Mu);
			}
	  case 5:
		switch (i)
		  {
		  case 0 : 
			if (Xi > Eta)
			  return((1.0-Xi)*(1.0-Eta) - Mu*(1.0-Eta));
			else
			  return((1.0-Xi)*(1.0-Eta) - Mu*(1.0-Xi));
		  case 1 : 
			if (Xi > Eta)
			  return(Xi*(1.0-Eta)       - Mu*Eta);
			else
			  return(Xi*(1.0-Eta)       - Mu*Xi);
		  case 2 : 
			if (Xi > Eta)
			  return(Xi*Eta             + Mu*Eta);
			else
			  return(Xi*Eta             + Mu*Xi);
		  case 3 : 			  
			if (Xi > Eta)
			  return((1.0-Xi)*Eta       - Mu*Eta);
			else
			  return((1.0-Xi)*Eta       - Mu*Xi);
		  case 4 : return(Mu);
		  }
	  case 6:
		switch (i)
		  {
		  case 0 : return((1.0-Xi-Eta)*(1.0-Mu));
		  case 1 : return(Xi*(1.0-Mu));
		  case 2 : return(Eta*(1.0-Mu));
		  case 3 : return((1.0-Xi-Eta)*Mu);
		  case 4 : return(Xi*Mu);
		  case 5 : return(Eta*Mu);
		  }
	  case 8:
		switch (i)
		  {
		  case 0 : return((1.0-Xi)*(1.0-Eta)*(1.0-Mu));
		  case 1 : return(Xi*(1.0-Eta)*(1.0-Mu));
		  case 2 : return(Xi*Eta*(1.0-Mu));
		  case 3 : return((1.0-Xi)*Eta*(1.0-Mu));
		  case 4 : return((1.0-Xi)*(1.0-Eta)*Mu);
		  case 5 : return(Xi*(1.0-Eta)*Mu);
		  case 6 : return(Xi*Eta*Mu);
		  case 7 : return((1.0-Xi)*Eta*Mu);
		  }
	  default:
		return (-1.0);
	  }
#endif
}

/****************************************************************************/
/*D
   GNs - General Shape function for nodes

   SYNOPSIS:
   INT GNs (INT n, COORD *ip_local, DOUBLE *result)

   PARAMETERS:
.  n - number of corners of the element
.  local - local COORDinates
.  result - vector of values

   DESCRIPTION:
   This function finds the value of the shape functions for the reference 
   elements at the given local coordinate.

   The shape functions fullfill
.n   Ni(node i) = 1
.n   Ni(node k) = 0, if k is not equal i.
   
   RETURN VALUE:
   INT
.n    0 if ok 
.n    1 if determinant of coordinate transformation too small.
D*/   
/****************************************************************************/

INT GNs (INT n, COORD *ip_local, DOUBLE *result)
{
    #ifdef __TWODIM__
    switch (n)
	  {
	  case 3:
		result[0] = 1.0-Xi-Eta;
		result[1] = Xi;
		result[2] = Eta;
		return (0);
	  case 4:
		result[0] = (1.0-Xi)*(1.0-Eta);
		result[1] = Xi*(1.0-Eta);
		result[2] = Xi*Eta;
		result[3] = (1.0-Xi)*Eta;
		return (0);
	  default:
		return (1);
	  }
    #endif
	
    #ifdef __THREEDIM__
	switch (n)
	  {
	  case 4:
		result[0] = 1.0-Xi-Eta-Mu;
		result[1] = Xi;
		result[2] = Eta;
		result[3] = Mu;
	  case 5:
		if (Xi > Eta)
		  {
			result[0] = (1.0-Xi)*(1.0-Eta) - Mu*(1.0-Eta);
			result[1] = Xi*(1.0-Eta)       - Mu*Eta;
			result[2] = Xi*Eta             + Mu*Eta;
			result[3] = (1.0-Xi)*Eta       - Mu*Eta;
			result[4] = Mu;
			return (0);
		  }
		else
		  {
			result[0] = (1.0-Xi)*(1.0-Eta) - Mu*(1.0-Xi);
			result[1] = Xi*(1.0-Eta)       - Mu*Xi;
			result[2] = Xi*Eta             + Mu*Xi;
			result[3] = (1.0-Xi)*Eta       - Mu*Xi;
			result[4] = Mu;
			return (0);
		  }
	  case 6:
		result[0] = (1.0-Xi-Eta)*(1.0-Mu);
		result[1] = Xi*(1.0-Mu);
		result[2] = Eta*(1.0-Mu);
		result[3] = (1.0-Xi-Eta)*Mu;
		result[4] = Xi*Mu;
		result[5] = Eta*Mu;
		return (0);
	  case 8:
		result[0] = (1.0-Xi)*(1.0-Eta)*(1.0-Mu);
		result[1] = Xi*(1.0-Eta)*(1.0-Mu);
		result[2] = Xi*Eta*(1.0-Mu);
		result[3] = (1.0-Xi)*Eta*(1.0-Mu);
		result[4] = (1.0-Xi)*(1.0-Eta)*Mu;
		result[5] = Xi*(1.0-Eta)*Mu;
		result[6] = Xi*Eta*Mu;
		result[7] = (1.0-Xi)*Eta*Mu;
		return (0);
	  default:
		return (1);
	  }
#endif
}

/****************************************************************************/
/*D
   D_GN - General Shape function for nodes

   SYNOPSIS:
   INT D_GN (INT n, INT i, COORD *ip_local, DOUBLE *derivative);

   PARAMETERS:
.  n - number of corners of the element
.  i - corner number (corner number [0..n-1])
.  ip_local - local COORDinates
.  derivative - derivative

   DESCRIPTION:
   This function computes the derivative of the shape functions GN.
   
   RETURN VALUE:
   INT
.n    0 if ok 
.n    1 if determinant of coordinate transformation too small.
D*/   
/****************************************************************************/

INT D_GN (INT n, INT i, COORD *ip_local, DOUBLE *derivative)
{
    #ifdef __TWODIM__
	switch (n)
	  {
	  case 3:
		switch (i) 
		  {
		  case 0: 
			derivative[0] = -1.0;
			derivative[1] = -1.0;
			return(0);
		  case 1: 
			derivative[0] = 1.0;
			derivative[1] = 0.0;
			return(0);
		  case 2: 
			derivative[0] = 0.0;
			derivative[1] = 1.0;
			return(0);
		  }
	  case 4:
		switch (i) 
		  {
		  case 0: 
			derivative[0] = -1.0+Eta;
			derivative[1] = -1.0+Xi;
			return(0);
		  case 1: 
			derivative[0] = 1.0-Eta;
			derivative[1] = -Xi;
			return(0);
		  case 2: 
			derivative[0] = -Eta;
			derivative[1] = 1.0-Xi;
			return(0);
		  case 3: 
			derivative[0] = Eta;
			derivative[1] = Xi;
			return(0);
		  }
	  }
    #endif

    #ifdef __THREEDIM__
	switch (n)
	  {
	  case 4:
		switch (i) 
		  {
		  case 0: 
			derivative[0] = -1.0;
			derivative[1] = -1.0;
			derivative[2] = -1.0;
			return(0);
		  case 1: 
			derivative[0] = 1.0;
			derivative[1] = 0.0;
			derivative[2] = 0.0;
			return(0);
		  case 2: 
			derivative[0] = 0.0;
			derivative[1] = 1.0;
			derivative[2] = 0.0;
			return(0);
		  case 3: 
			derivative[0] = 0.0;
			derivative[1] = 0.0;
			derivative[2] = 1.0;
			return(0);
		  }
	  case 5:
		switch (i) 
		  {
		  case 0: 
			if (Xi > Eta)
			  {
				derivative[0] = -(1.0-Eta);
				derivative[1] = -(1.0-Xi) + Mu;
				derivative[2] = -(1.0-Eta);
				return(0);
			  }
			else
			  {
				derivative[0] = -(1.0-Eta) + Mu;
				derivative[1] = -(1.0-Xi);
				derivative[2] = -(1.0-Xi);
				return(0);
			  }
		  case 1: 
			if (Xi > Eta)
			  {
				derivative[0] = (1.0-Eta);
				derivative[1] = -Xi - Mu;
				derivative[2] = -Eta;
				return(0);
			  }
			else
			  {
				derivative[0] = (1.0-Eta) - Mu;
				derivative[1] = -Xi;
				derivative[2] = -Xi;
				return(0);
			  }
		  case 2:
			if (Xi > Eta)
			  {
				derivative[0] = Eta;
				derivative[1] = Xi + Mu;
				derivative[2] = Eta;
				return(0);
			  }
			else
			  {
				derivative[0] = Eta + Mu;
				derivative[1] = Xi;
				derivative[2] = Xi;
				return(0);
			  } 
		  case 3: 
			if (Xi > Eta)
			  {
				derivative[0] = -Eta;
				derivative[1] = 1.0-Xi - Mu;
				derivative[2] = -Eta;
				return(0);
			  }
			else
			  {
				derivative[0] = -Eta - Mu;
				derivative[1] = 1.0-Xi;
				derivative[2] = -Xi;
				return(0);
			  } 
          case 4:
			derivative[0] = 0.0;
			derivative[1] = 0.0;
			derivative[2] = 1.0;
			return(0);
		  }
	  case 6:
		switch (i) 
		  {
		  case 0: 
			derivative[0] = -(1.0-Mu);
			derivative[1] = -(1.0-Mu);
			derivative[2] = -(1.0-Xi-Eta);
			return(0);
		  case 1: 
			derivative[0] = (1.0-Mu);
			derivative[1] = 0.0;
			derivative[2] = -Xi;
			return(0);
		  case 2: 
			derivative[0] = 0.0;
			derivative[1] = (1.0-Mu);
			derivative[2] = -Eta;
			return(0);
		  case 3: 
			derivative[0] = -Mu;
			derivative[1] = -Mu;
			derivative[2] = 1.0-Xi-Eta;
            return(0);
          case 4:
			derivative[0] = Mu;
			derivative[1] = 0.0;
			derivative[2] = Xi;
			return(0);
          case 5:
			derivative[0] = 0.0;
			derivative[1] = Mu;
			derivative[2] = Eta;
			return(0);
		  }
	  case 8:
		switch (i) 
		  {
		  case 0: 
			derivative[0] = -(1.0-Eta)*(1.0-Mu);
			derivative[1] = -(1.0-Xi)*(1.0-Mu);
			derivative[2] = -(1.0-Xi)*(1.0-Eta);
			return(0);
		  case 1: 
			derivative[0] = (1.0-Eta)*(1.0-Mu);
			derivative[1] = -(Xi)*(1.0-Mu);
			derivative[2] = -(Xi)*(1.0-Eta);
			return(0);
		  case 2: 
			derivative[0] = (Eta)*(1.0-Mu);
			derivative[1] = (Xi)*(1.0-Mu);
			derivative[2] = -Xi*Eta;
			return(0);
		  case 3: 
			derivative[0] = -(Eta)*(1.0-Mu);
			derivative[1] = (1.0-Xi)*(1.0-Mu);
			derivative[2] = -(1.0-Xi)*(Eta);
            return(0);
          case 4:
			derivative[0] = -(1.0-Eta)*(-Mu);
			derivative[1] = -(1.0-Xi)*(-Mu);
			derivative[2] = (1.0-Xi)*(1.0-Eta);
			return(0);
		  case 5: 
			derivative[0] = (1.0-Eta)*Mu;
			derivative[1] = -(Xi)*Mu;
			derivative[2] = (Xi)*(1.0-Eta);
			return(0);
		  case 6: 
			derivative[0] = (Eta)*Mu;
			derivative[1] = (Xi)*Mu;
			derivative[2] = Xi*Eta;
			return(0);
		  case 7: 
			derivative[0] = -(Eta)*Mu;
			derivative[1] = (1.0-Xi)*Mu;
			derivative[2] = (1.0-Xi)*Eta;
            return(0);
		  }
	  }
    #endif

    return(1);
}

/****************************************************************************/
/*D
   LMP - local midpoint

   SYNOPSIS:
   COORD *LMP (INT n);

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
		case 5:	return (LMP_Pyramid);
		case 6:	return (LMP_Prism);
		case 8:	return (LMP_Hexahedron);
	}
#endif
}

/****************************************************************************/
/*D
   GlobalToLocal - Transform global coordinates to local

   SYNOPSIS:
   INT GlobalToLocal (INT n, const COORD **Corners, const COORD *EvalPoint, 
   COORD *LocalCoord);

   PARAMETERS:
.  n - number of corners
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

INT GlobalToLocal (INT n, const COORD **Corners, 
				   const COORD *EvalPoint, COORD *LocalCoord)
{
	COORD_VECTOR tmp,diff,M[DIM],IM[DIM];
	DOUBLE s,IMdet;
	INT i;
	
	V_DIM_SUBTRACT(EvalPoint,Corners[0],diff);
	if (n == DIM+1)
	  {
		TRANSFORMATION(DIM+1,Corners,LocalCoord,M);
		M_DIM_INVERT(M,IM,IMdet);
		MM_TIMES_V_DIM(IM,diff,LocalCoord);	
		return(0);
	  }
	V_DIM_CLEAR(LocalCoord);
	TRANSFORMATION(n,Corners,LocalCoord,M);
	M_DIM_INVERT(M,IM,IMdet);
	MM_TIMES_V_DIM(IM,diff,LocalCoord);	
	for (i=0; i<MAX_ITER; i++)
	  {
		LOCAL_TO_GLOBAL (n,Corners,LocalCoord,tmp);
		V_DIM_SUBTRACT(tmp,EvalPoint,diff);
		V_DIM_EUKLIDNORM(diff,s);
		if (s <= SMALL_DIFF) 
			return (0);
		TRANSFORMATION(n,Corners,LocalCoord,M);
		M_DIM_INVERT(M,IM,IMdet);
		MM_TIMES_V_DIM(IM,diff,tmp);
		V_DIM_SUBTRACT(LocalCoord,tmp,LocalCoord);
	  }

	return(1);
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
			case 0 : return(t-1);
			case 1 : return(1-t);
			case 2 : return(t);
			case 3 : return(-t);
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
			case 0 : return(s-1);
			case 1 : return(-s);
			case 2 : return(s);
			case 3 : return(1-s);
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
    if (fabs(detJ)<=SMALL_DET) return(1);
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
   L2GDerivative2d - Derivative of LOCAL_TO_GLOBAL

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

	if (TetraSideNormals (theElement,(COORD **)theCorners,theNormal)) return (1);
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
static INT FindCrossParam3D (COORD *p1, COORD *p2, COORD *p3, COORD *p4, DOUBLE_VECTOR v, COORD *param)
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
			
			if (FindCrossParam3D((COORD *)CornerPoints[FrontCorner[0]],
								 (COORD *)CornerPoints[FrontCorner[1]],
								 (COORD *)CornerPoints[BackCorner[0]],
								 (COORD *)CornerPoints[BackCorner[1]],
								 conv,param)) return (1);
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
