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
/*			  internet: ug@ica3.uni-stuttgart.de					*/
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

static INT MirrorAtPlane (COORD *in, COORD *pp, COORD *pn, COORD *out)
{
	COORD_VECTOR a;
	
	V3_SUBTRACT(pp,in,a)
	if (V3_Project(a,pn,out))
	V3_LINCOMB(1.0,in,2.0,out,out)
	
	return (0);
}

#define V3_TRI_CM(a,b,c,e)			   {(e)[0] = 0.333333333333*((a)[0]+(b)[0]+(c)[0]);\
										(e)[1] = 0.333333333333*((a)[1]+(b)[1]+(c)[1]);\
										(e)[2] = 0.333333333333*((a)[2]+(b)[2]+(c)[2]);}
#define V3_QUA_CM(a,b,c,d,e)		   {(e)[0] = 0.25*((a)[0]+(b)[0]+(c)[0]+(d)[0]);\
										(e)[1] = 0.25*((a)[1]+(b)[1]+(c)[1]+(d)[1]);\
										(e)[2] = 0.25*((a)[2]+(b)[2]+(c)[2]+(d)[2]);}

INT FV_AliTetInfo (COORD **CornerPoints, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv, COORD_VECTOR GIP[MAX_EDGES_OF_ELEM], COORD_VECTOR LIP[MAX_EDGES_OF_ELEM])
{
	COORD sp, alpha, check[2], M[9], Inv[9];
	COORD_VECTOR a, b, c, d, e, cm, normal, param, EdgeMidPoints[6], SideMidPoints[4];
	INT i, help, noutflow, ninflow, outflow[4], inflow[4], OpEdge[3], GrEdge[3], side[3], OpCorner, corner[3], inverted, First;
	INT BackEdge, FrontEdge, BackCorner[2], FrontCorner[2], EdgeF0B0, EdgeF0B1, EdgeF1B0, EdgeF1B1, flags, changed;
	
	/* reset areas and integrationpoints */
	for (i=0; i<MAX_EDGES_OF_ELEM; i++)
	{
		V3_CLEAR(Area[i])
		V3_CLEAR(GIP[i])
		V3_COPY(LIP[i],LIP[i])
	}
	
	/* edge mid points */
	for (i=0; i<6; i++)	V3_LINCOMB(0.5,CornerPoints[CornerOfEdge[i][0]],0.5,CornerPoints[CornerOfEdge[i][1]],EdgeMidPoints[i])
			
	/* side mid points */
	for (i=0; i<MAX_SIDES_OF_ELEM; i++)
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

INT GetSkewedUIP (COORD_VECTOR *theCorners, COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM])
{
	COORD_VECTOR lconv;
	COORD alpha;
	COORD M[9],I[9];
	INT flags, i;

	V3_SUBTRACT(theCorners[1],theCorners[0],M)
	V3_SUBTRACT(theCorners[2],theCorners[0],M+3)
	V3_SUBTRACT(theCorners[3],theCorners[0],M+6)
	if (M3_Invert(I,M)) return (1);
	for (i=0; i<MAX_EDGES_OF_ELEM; i++)
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

