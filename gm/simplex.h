// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

/***************************************************************************
*																		   *
*	File:		 simplex.h												   *
*																		   *
*	Purpose:	 definition and algorithms for simplices header file	   *
*																		   *
*	Author:          Juergen Bey											   *
*				 Mathematisches Institut								   *
*				 Auf der Morgenstelle 10								   *
*				 7400 Tuebingen                                                                                    *
*				 email : juergen@tue-num2.mathematik.uni-tuebingen.de	   *
*																		   *
*	History:	 10.12.92  begin										   *
*																		   *
*	Remarks:															   *
*																		   *
***************************************************************************/

#ifndef __SIMPLEX__
#define __SIMPLEX__

/***************************************************************************
*																		   *
*								 Includes								   *
*																		   *
***************************************************************************/

#ifdef __MISC__
#include "misc.h"
#endif

#ifndef __GM__
#include "gm.h"
#endif

/***************************************************************************
*																		   *
*								 Constants								   *
*																		   *
***************************************************************************/


/***************************************************************************
*																		   *
*								  Macros								   *
*																		   *
***************************************************************************/

#define CORNERS_OF_TETRASIDE    3
#define EDGES_OF_TETRASIDE              3

/***************************************************************************
*																		   *
*							  Data structures							   *
*																		   *
***************************************************************************/


/***************************************************************************
*																		   *
*							  Global variables							   *
*																		   *
***************************************************************************/

extern COORD TRefCoord[MAX_CORNERS_OF_ELEM][DIM];
extern COORD_VECTOR LIP[MAX_EDGES_OF_ELEM];
extern COORD_VECTOR LSIP[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_SIDE];

extern INT CornerOfSide[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_SIDE];
extern INT SignOfSide[2][MAX_SIDES_OF_ELEM];
extern INT CornerOfEdge[MAX_EDGES_OF_ELEM][CORNERS_OF_EDGE];
extern INT CornerOfOppEdge[MAX_EDGES_OF_ELEM][CORNERS_OF_EDGE];
extern INT EdgeWithCorners[MAX_CORNERS_OF_ELEM][MAX_CORNERS_OF_ELEM];
extern INT SideWithEdge[MAX_EDGES_OF_ELEM][MAX_SIDES_OF_EDGE];
extern INT EdgeOfSide[MAX_SIDES_OF_ELEM][EDGES_OF_TETRASIDE];
extern INT OppositeCorner[MAX_SIDES_OF_ELEM];
extern INT OppositeSide[MAX_SIDES_OF_ELEM];
extern INT CornerOppCorners[MAX_CORNERS_OF_ELEM][MAX_CORNERS_OF_ELEM][MAX_CORNERS_OF_ELEM];
extern INT CornerIndex[24][MAX_CORNERS_OF_ELEM];
extern INT SideCornerIndex[MAX_SIDES_OF_ELEM][24][MAX_CORNERS_OF_SIDE];
extern INT CondensedEdgeOfSide[MAX_SIDES_OF_ELEM];
extern INT TriSectionEdge[64][2];
extern INT CornerOfSideInv[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_ELEM];
extern INT MidNodeOfEdge[MAX_EDGES_OF_ELEM];
extern INT OppositeEdge[MAX_EDGES_OF_ELEM];
extern INT EdgesOfCorner[MAX_CORNERS_OF_ELEM][3];
extern INT SideEdgesOfEdge[MAX_EDGES_OF_ELEM][2][2];
extern INT NextEdgeToSIP[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_SIDE][2];
extern INT SideOf2Edges[MAX_EDGES_OF_ELEM][MAX_EDGES_OF_ELEM];
extern INT EdgeOf2Sides[MAX_SIDES_OF_ELEM][MAX_SIDES_OF_ELEM];
/***************************************************************************
*																		   *
*							  Extern functions							   *
*																		   *
***************************************************************************/


#endif

/***************************************************************************
*																		   *
*								 The End								   *
*																		   *
***************************************************************************/
