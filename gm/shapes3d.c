// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/****************************************************************************/
/*																			*/
/* File:	  shapes3d.c														*/
/*																			*/
/* Purpose:   shape functions for triangles and quadrangles 				*/
/*																			*/
/* Author:	  Peter Bastian 												*/
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: bastian@iwr1.iwr.uni-heidelberg.de					*/
/*																			*/
/* History:   08.04.92 begin, ug version 2.0								*/
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

#include "evm.h"
#include "shapes3d.h"
#include "simplex.h"

#ifdef __TWODIM__
#error this source file is for 3D ONLY
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
   DOUBLE N (INT i, COORD *LocalCoord);

   PARAMETERS:
.  i - corner number [0..3]
.  LocalCoord - local coordinates
  
   DESCRIPTION:
   This function finds the linear shape functions Ni(s,t) to approximate the
   solution u in the integration point ip for tetrahedrons
   
.n   uip(s,t) = SUM Ni(s,t)*ui

   where the sum runs over all nodes of the element to which the considered
   ip belongs. 

   RETURN VALUE:
   DOUBLE
D*/
/****************************************************************************/

DOUBLE N (INT i, COORD *LocalCoord)
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

/****************************************************************************/
/*D
   GlobalToLocal3d - Transform global coordinates to local

   SYNOPSIS:
   INT GlobalToLocal3d (COORD **Corners,COORD *EvalPoint, 
   COORD *LocalCoord);

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

INT GlobalToLocal3d (COORD **Corners,COORD *EvalPoint, COORD *LocalCoord)
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

/****************************************************************************/
/*D
   TetraSideNormals - Calculate inner normals of tetrahedra

   SYNOPSIS:
   INT TetraSideNormals (COORD **theCorners, COORD_VECTOR theNormals[MAX_SIDES_OF_ELEM]);

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

INT TetraSideNormals (COORD **theCorners, COORD_VECTOR theNormals[MAX_SIDES_OF_ELEM])
{
	COORD_VECTOR a, b;
	COORD h;
	INT j,k;

	for (j=0; j<MAX_CORNERS_OF_ELEM; j++)
	{
		k = OppositeSide[j];
		V3_SUBTRACT(theCorners[(j+1)%MAX_CORNERS_OF_ELEM],theCorners[(j+2)%MAX_CORNERS_OF_ELEM],a)
		V3_SUBTRACT(theCorners[(j+1)%MAX_CORNERS_OF_ELEM],theCorners[(j+3)%MAX_CORNERS_OF_ELEM],b)
		V3_VECTOR_PRODUCT(a,b,theNormals[k])
		V3_Normalize(theNormals[k]);
		V3_SUBTRACT(theCorners[j],theCorners[(j+1)%MAX_CORNERS_OF_ELEM],a)
		V3_SCALAR_PRODUCT(theNormals[k],a,h);
		if (ABS(h)<SMALL_C) return (1);
		if (h<0.0)
			V3_SCALE(-1.0,theNormals[k]);
	}

	return (0);

}

/****************************************************************************/
/*D
   TetMaxSideAngle - Calculate maximal side angle of Tetrahedron

   SYNOPSIS:
   INT TetMaxSideAngle (COORD **theCorners, COORD *MaxAngle);

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

INT TetMaxSideAngle (COORD **theCorners, COORD *MaxAngle)
{
	COORD_VECTOR theNormal[MAX_SIDES_OF_ELEM];
	COORD max,help;
	INT i;

	if (TetraSideNormals (theCorners,theNormal)) return (1);
	max = -1.0;
	for (i=0; i<MAX_EDGES_OF_ELEM; i++)
	{
		V3_SCALAR_PRODUCT(theNormal[SideWithEdge[i][0]],theNormal[SideWithEdge[i][1]],help)
		max = MAX(help,max);
	}
	max = MIN(max,1.0);
	*MaxAngle = 180.0/PI*acos(-max);
	
	return (0);
}

/****************************************************************************/
/*D
   TetAngleAndLength - Calculates side angle and length of edge of Tetrahedron

   SYNOPSIS:
   INT TetAngleAndLength (COORD **theCorners, COORD *Angle, COORD *Length);

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

INT TetAngleAndLength (COORD **theCorners, COORD *Angle, COORD *Length)
{
	COORD_VECTOR theNormals[MAX_SIDES_OF_ELEM],theEdge[MAX_EDGES_OF_ELEM];
	COORD h;
	INT j,k;

	for (j=0; j<MAX_EDGES_OF_ELEM; j++)
	{
		V3_SUBTRACT(theCorners[CornerOfEdge[j][1]],theCorners[CornerOfEdge[j][0]],theEdge[j])
		V3_EUKLIDNORM(theEdge[j],Length[j])
	}
	for (j=0; j<MAX_SIDES_OF_ELEM; j++)
	{
		V3_VECTOR_PRODUCT(theEdge[EdgeOfSide[j][0]],theEdge[EdgeOfSide[j][1]],theNormals[j])
		V3_Normalize(theNormals[j]);
		k = EdgesOfCorner[OppositeCorner[j]][0];
		V3_SCALAR_PRODUCT(theNormals[j],theEdge[k],h)
		if (ABS(h)<SMALL_C) return (1);
		if ( (h<0.0 && CornerOfEdge[k][1]==OppositeCorner[j]) ||
			 (h>0.0 && CornerOfEdge[k][0]==OppositeCorner[j])	 )
			V3_SCALE(-1.0,theNormals[j]);
	}
	for (j=0; j<MAX_EDGES_OF_ELEM; j++)
	{
		V3_SCALAR_PRODUCT(theNormals[SideWithEdge[j][0]],theNormals[SideWithEdge[j][1]],Angle[j])
		Angle[j] = MAX(Angle[j],-1.0);
		Angle[j] = MIN(Angle[j], 1.0);
		Angle[j] = (COORD)acos((double)Angle[j]);
	}
	
	return (0);
}

/****************************************************************************/
/*D
   TetraDerivative - Calculates gradient of shape function for tetrahedron

   SYNOPSIS:
   INT TetraDerivative (COORD **theCorners, COORD_VECTOR theGradient[MAX_CORNERS_OF_ELEM])

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

INT TetraDerivative (const COORD **theCorners, COORD_VECTOR theGradient[MAX_CORNERS_OF_ELEM])
{
	COORD_VECTOR a, b;
	COORD h;
	INT j;

	for (j=0; j<MAX_CORNERS_OF_ELEM; j++)
	{
		V3_SUBTRACT(theCorners[(j+1)%MAX_CORNERS_OF_ELEM],theCorners[(j+2)%MAX_CORNERS_OF_ELEM],a)
		V3_SUBTRACT(theCorners[(j+1)%MAX_CORNERS_OF_ELEM],theCorners[(j+3)%MAX_CORNERS_OF_ELEM],b)
		V3_VECTOR_PRODUCT(a,b,theGradient[j])
		V3_Normalize(theGradient[j]);
		V3_SUBTRACT(theCorners[j],theCorners[(j+1)%MAX_CORNERS_OF_ELEM],a)
		V3_SCALAR_PRODUCT(theGradient[j],a,h)
		if (ABS(h)<SMALL_C) return (1);
		V3_SCALE(1/h,theGradient[j])
	}

	return (0);
}

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

INT TetraVolume (COORD **theCorners, COORD *volume)
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

/****************************************************************************/
/*D
   FV_TetInfo - Calculate control volume subsurfaces and global integration points

   SYNOPSIS:
   INT FV_TetInfo (COORD **theCorners, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], 
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

INT FV_TetInfo (COORD **theCorners, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], COORD_VECTOR GIP[MAX_EDGES_OF_ELEM])
{
	COORD_VECTOR emp[MAX_EDGES_OF_ELEM], diff, a, b;
	COORD sp;
	INT i;

	for (i=0; i<MAX_EDGES_OF_ELEM; i++)
	{
		V3_LINCOMB(0.5,theCorners[CornerOfEdge[i][0]],0.5,theCorners[CornerOfEdge[i][1]],emp[i])
		V3_SUBTRACT (theCorners[CornerOfOppEdge[i][0]],emp[i],a)
		V3_SUBTRACT (theCorners[CornerOfOppEdge[i][1]],emp[i],b)
		V3_VECTOR_PRODUCT(a,b,Area[i])
		V3_SUBTRACT (theCorners[CornerOfEdge[i][1]],theCorners[CornerOfEdge[i][0]],diff)
		V3_SCALAR_PRODUCT(Area[i],diff,sp)
		if (sp>0.0)
			V3_SCALE(1/12.0,Area[i])
		else
			V3_SCALE(-1/12.0,Area[i])
	}
	for (i=0; i<MAX_EDGES_OF_ELEM; i++)
		V3_LINCOMB(17.0/24.0,emp[i],7.0/24.0,emp[OppositeEdge[i]],GIP[i])

	return (0);
}

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

INT FV_TetInfo_for_conv (COORD **CornerPoints, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], COORD_VECTOR GIP[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM], COORD_VECTOR conv)
{
	COORD sp, spn, spz, alpha1, alpha2;
	COORD_VECTOR a, b, c, normal, ex, ey;
	COORD_VECTOR EdgeMidPoints[6], SideMidPoints[4];
	INT i, j, k, help, noutflow, ninflow, outflow[4], inflow[4], edge, inverted, side;
	
	/* reset areas and integrationpoints */
	for (i=0; i<MAX_EDGES_OF_ELEM; i++)
	{
		V3_CLEAR(Area[i])
		V3_CLEAR(GIP[i])
		V3_CLEAR(LUIP[i])
	}
	
	/* edge mid points */
	for (i=0; i<6; i++)	V3_LINCOMB(0.5,CornerPoints[CornerOfEdge[i][0]],0.5,CornerPoints[CornerOfEdge[i][1]],EdgeMidPoints[i])
			
	/* side mid points */
	for (i=0; i<MAX_SIDES_OF_ELEM; i++)
	{
		V3_ADD(CornerPoints[CornerOfSide[i][0]],CornerPoints[CornerOfSide[i][1]],a)
		V3_ADD(CornerPoints[CornerOfSide[i][2]],a,a)
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
				edge = EdgeWithCorners[CornerOfSide[inflow[0]][i]][OppositeCorner[inflow[0]]];								/* edge: P0 ~ P3	*/
				V3_SUBTRACT(SideMidPoints[inflow[0]],EdgeMidPoints[edge],a)													/* P012 - P3		*/
				V3_SUBTRACT(CornerPoints[CornerOfSide[inflow[0]][(i+2)%3]],CornerPoints[CornerOfSide[inflow[0]][(i+1)%3]],b)/* P2 - P1 			*/
				V3_VECTOR_PRODUCT(a,b,Area[edge])
				V3_SCALE(0.1666666666,Area[edge])

				/* P2 = corners of element-side[i]  */
				/* calc alpha1 */
				edge = EdgeWithCorners[CornerOfSide[inflow[0]][(i+1)%3]][OppositeCorner[inflow[0]]];						/* edge: P0 ~ P3 	*/
				V3_SUBTRACT(EdgeMidPoints[edge],CornerPoints[CornerOfSide[inflow[0]][(i+2)%3]],a)							/* P03 - P1 		*/
				V3_VECTOR_PRODUCT(conv,a,c)
				edge = EdgeWithCorners[CornerOfSide[inflow[0]][(i+1)%3]][CornerOfSide[inflow[0]][(i+2)%3]];					/* edge: P0 ~ P1 */
				V3_SUBTRACT(EdgeMidPoints[edge],CornerPoints[CornerOfSide[inflow[0]][i]],b)									/* P01 - P2 */
				V3_SCALAR_PRODUCT(b,c,spn)
				if (spn==0.0) continue;
				V3_SUBTRACT(CornerPoints[OppositeCorner[inflow[0]]],CornerPoints[CornerOfSide[inflow[0]][i]],b)				/* P3 - P2 */
				V3_SCALAR_PRODUCT(b,c,spz)
				alpha1 = spz/spn;
				assert (alpha1>=0.0);
				
				/* calc alpha2 */
				edge = EdgeWithCorners[CornerOfSide[inflow[0]][(i+2)%3]][OppositeCorner[inflow[0]]];						/* edge: P1 ~ P3 	*/
				V3_SUBTRACT(EdgeMidPoints[edge],CornerPoints[CornerOfSide[inflow[0]][(i+1)%3]],a)							/* P13 - P0 		*/
				V3_VECTOR_PRODUCT(conv,a,c)
				edge = EdgeWithCorners[CornerOfSide[inflow[0]][(i+1)%3]][CornerOfSide[inflow[0]][(i+2)%3]];					/* edge: P0 ~ P1 	*/
				V3_SUBTRACT(EdgeMidPoints[edge],CornerPoints[CornerOfSide[inflow[0]][i]],b)									/* P01 - P2 		*/
				V3_SCALAR_PRODUCT(b,c,spn)
				assert (spn!=0.0);
				V3_SUBTRACT(CornerPoints[OppositeCorner[inflow[0]]],CornerPoints[CornerOfSide[inflow[0]][i]],b)				/* P3 - P2 			*/
				V3_SCALAR_PRODUCT(b,c,spz)
				alpha2 = spz/spn;
				assert (alpha2>=0.0);
				
				/* take the right triangle */
				V3_SUBTRACT(CornerPoints[CornerOfSide[inflow[0]][i]],CornerPoints[OppositeCorner[inflow[0]]],a)				/* P2 - P3 			*/
				edge = EdgeWithCorners[CornerOfSide[inflow[0]][(i+1)%3]][CornerOfSide[inflow[0]][(i+2)%3]];					/* edge: P0 ~ P1 	*/
				V3_SUBTRACT(EdgeMidPoints[edge],CornerPoints[CornerOfSide[inflow[0]][i]],b)									/* P01 - P2 		*/
				V3_VECTOR_PRODUCT(a,b,c)
				if (alpha1<alpha2)		j=1;
				else if (alpha2<alpha1) j=0;
				else
				{
					side = OppositeSide[CornerOfSide[inflow[0]][i]];														/* side 013 		*/
					V3_SUBTRACT(SideMidPoints[side],CornerPoints[CornerOfSide[inflow[0]][(i+2)%3]],a)						/* P013 - P1 		*/
					V3_EUKLIDNORM(a,alpha1)
					V3_SUBTRACT(SideMidPoints[side],CornerPoints[CornerOfSide[inflow[0]][(i+1)%3]],a)						/* P013 - P0 		*/
					V3_EUKLIDNORM(a,alpha2)
					if (alpha1<alpha2)		j=1;
					else if (alpha2<alpha1) j=0;
					else continue;
				}
				if (j)
				{
					V3_SCALE(alpha1/18.0,c)
					edge = EdgeWithCorners[CornerOfSide[inflow[0]][i]][OppositeCorner[inflow[0]]];							/* edge: P0 ~ P3 	*/
					V3_SUBTRACT(Area[edge],c,Area[edge])			
					edge = EdgeWithCorners[CornerOfSide[inflow[0]][(i+2)%3]][OppositeCorner[inflow[0]]];					/* edge: P1 ~ P3 	*/
					V3_ADD(Area[edge],c,Area[edge])			
				}
				else
				{
					V3_SCALE(alpha2/18.0,c)
					edge = EdgeWithCorners[CornerOfSide[inflow[0]][i]][OppositeCorner[inflow[0]]];							/* edge: P0 ~ P3 	*/
					V3_ADD(Area[edge],c,Area[edge])			
					edge = EdgeWithCorners[CornerOfSide[inflow[0]][(i+2)%3]][OppositeCorner[inflow[0]]];					/* edge: P1 ~ P3 	*/
					V3_SUBTRACT(Area[edge],c,Area[edge])			
				}
			}
			
			/* orientate the areas */
			for (i=0; i<6; i++)	
			{
				V3_SUBTRACT(CornerPoints[CornerOfEdge[i][1]],CornerPoints[CornerOfEdge[i][0]],a)
				V3_SCALAR_PRODUCT(a,Area[i],sp)
				if (sp<0.0)
					V3_SCALE(-1.0,Area[i])
			}
			
			/* set LUIP */
			for (i=0; i<3; i++)	
			{			
				edge = EdgeWithCorners[CornerOfSide[inflow[0]][i]][OppositeCorner[inflow[0]]];
				if (inverted)	V3_COPY(TRefCoord[OppositeCorner[inflow[0]]],LUIP[edge])
				else			V3_COPY(TRefCoord[CornerOfSide[inflow[0]][i]],LUIP[edge])
			}
			break;
		default:
			return (1);		
	}
	
	return (0);
}

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

INT Side_TetInfo (COORD **theCorners, INT side, COORD_VECTOR Area, COORD_VECTOR GIP[3])
{
	COORD_VECTOR a,b,c;
	COORD scalarprd;
	INT i0, i1, i2, i3;
	
	i0 = CornerOfSide[side][0];
	i1 = CornerOfSide[side][1];
	i2 = CornerOfSide[side][2];
	i3 = OppositeCorner[side];

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

INT GSUIP (COORD **theCorners, COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM])
{
	COORD_VECTOR a, lconv;
	COORD alpha;
	COORD M[9],I[9];
	INT flags, i;

	V3_SUBTRACT(theCorners[1],theCorners[0],M)
	V3_SUBTRACT(theCorners[2],theCorners[0],M+3)
	V3_SUBTRACT(theCorners[3],theCorners[0],M+6)
	if (M3_Invert(I,M)) return (1);
	for (i=0; i<MAX_EDGES_OF_ELEM; i++)
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

INT GFUIP (COORD **theCorners, COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM])
{
	COORD_VECTOR lconv;
	COORD sp, min;
	COORD M[9],I[9];
	INT j, i, ip[MAX_CORNERS_OF_ELEM], n;

	V3_SUBTRACT(theCorners[1],theCorners[0],M)
	V3_SUBTRACT(theCorners[2],theCorners[0],M+3)
	V3_SUBTRACT(theCorners[3],theCorners[0],M+6)
	if (M3_Invert(I,M)) return (1);
	for (i=0; i<MAX_EDGES_OF_ELEM; i++)
	{
		M3_TIMES_V3(I,conv[i],lconv)
		min = MAX_C; n = 0;
		for (j=0; j<MAX_CORNERS_OF_ELEM; j++)
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

INT GCUIP (COORD **theCorners, COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM])
{
	COORD_VECTOR a, lconv, SUIP;
	COORD alpha, sp, min;
	COORD M[9],I[9];
	INT flags, i, j, k; 

	V3_SUBTRACT(theCorners[1],theCorners[0],M)
	V3_SUBTRACT(theCorners[2],theCorners[0],M+3)
	V3_SUBTRACT(theCorners[3],theCorners[0],M+6)
	if (M3_Invert(I,M)) return (1);
	for (i=0; i<MAX_EDGES_OF_ELEM; i++)
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
		for (j=0; j<MAX_CORNERS_OF_ELEM; j++)
		{
			V3_EUKLIDNORM_OF_DIFF(SUIP,TRefCoord[j],sp)
			if (min>sp)
			{
				k=j;
				min=sp;
			}
		}
		assert(k>=0 && k<MAX_CORNERS_OF_ELEM);
		LUIP[i][0] = TRefCoord[k][0];
		LUIP[i][1] = TRefCoord[k][1];
		LUIP[i][2] = TRefCoord[k][2];
	}
	
	return(0);
}

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

INT COPYIP (COORD **theCorners, COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM])
{
	INT i;
	
	for (i=0; i<MAX_EDGES_OF_ELEM; i++)
		V3_COPY(LIP[i],LUIP[i])
	
	return (0);
}

