// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/***************************************************************************
*																		   *
*	File:		 simplex.c												   *
*																		   *
*	Purpose:	 definition and algorithms for simlices                                    *
*																		   *
*	Author:          Juergen Bey											   *
*				 Mathematisches Institut								   *
*				 Auf der Morgenstelle 10								   *
*				 7400 Tuebingen                                                                                    *
*			     email: ug@ica3.uni-stuttgart.de					       *
*																		   *
*	History:	 10.12.92  begin										   *
*																		   *
*	Remarks:															   *
*																		   *
***************************************************************************/




/***************************************************************************
*																		   *
*								 Includes								   *
*																		   *
***************************************************************************/

#include  "simplex.h"
#include  "gm.h"


#ifdef __TWODIM__
#error this source file is for 3D ONLY
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

/***************************************************************************
*																		   *
*								 Tetrahedras							   *
*																		   *
***************************************************************************/

/* definition of a reference tetrahedron for local coordinates */
COORD_VECTOR TRefCoord[MAX_CORNERS_OF_ELEM] = {{0.0,0.0,0.0},{1.0,0.0,0.0},
                                               {0.0,1.0,0.0},{0.0,0.0,1.0}};

/* local itegration points in reference tet */
COORD_VECTOR LIP[MAX_EDGES_OF_ELEM] = {{17.0/48.0,7.0/48.0,7.0/48.0},
                                       {17.0/48.0,17.0/48.0,7.0/48.0},
                                       {7.0/48.0,17.0/48.0,7.0/48.0},
                                       {7.0/48.0,7.0/48.0,17.0/48.0},
                                       {17.0/48.0,7.0/48.0,17.0/48.0},
                                       {7.0/48.0,17.0/48.0,17.0/48.0}};

/* local surface itegration points in reference tet */
COORD_VECTOR LSIP[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_SIDE] =        {{{ 5.0/24.0, 5.0/24.0,              0.0},{14.0/24.0, 5.0/24.0,              0.0},{ 5.0/24.0,14.0/24.0,              0.0}},
                                                                    {{14.0/24.0, 5.0/24.0, 5.0/24.0},{ 5.0/24.0,14.0/24.0, 5.0/24.0},{ 5.0/24.0, 5.0/24.0,14.0/24.0}},
                                                                    {{              0.0,14.0/24.0, 5.0/24.0},{              0.0, 5.0/24.0,14.0/24.0},{              0.0, 5.0/24.0, 5.0/24.0}},
                                                                    {{ 5.0/24.0,      0.0,14.0/24.0},{ 5.0/24.0,      0.0, 5.0/24.0},{14.0/24.0,      0.0, 5.0/24.0}}};

/* get next two edges of surface itegration point */
INT NextEdgeToSIP[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_SIDE][2] = {{{0,2},{0,1},{1,2}},
                                                                {{1,4},{1,5},{4,5}},
                                                                {{2,5},{3,5},{2,3}},
                                                                {{3,4},{0,3},{0,4}}};

/* The indices of the corners of each side. */
INT CornerOfSide[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_SIDE] = {{0,2,1}, {1,2,3}, {2,0,3}, {3,0,1}};
/* gives the position of specified corner on specified side or -1 if it does not ly on that side */
INT CornerOfSideInv[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_ELEM]  = {{0,2,1,-1}, {-1,0,1,2}, {1,-1,0,2}, {1,2,-1,0}};

/* The indices of the corner opposite to three corners. */
INT CornerOppCorners[MAX_CORNERS_OF_ELEM][MAX_CORNERS_OF_ELEM][MAX_CORNERS_OF_ELEM] =
{ {{-1,-1,-1,-1}, {-1,-1, 3, 2}, {-1, 3,-1, 1}, {-1, 2, 1,-1}},
  {{-1,-1, 3, 2}, {-1,-1,-1,-1}, { 3,-1,-1, 0}, { 2,-1, 0,-1}},
  {{-1, 3,-1, 1}, { 3,-1,-1, 0}, {-1,-1,-1,-1}, { 1, 0,-1,-1}},
  {{-1, 2, 1,-1}, { 2,-1, 0,-1}, { 1, 0,-1,-1}, {-1,-1,-1,-1}} };


/* The indices of the corners of each edge and its opposite edge */
/* dont change order !											 */
INT CornerOfEdge[MAX_EDGES_OF_ELEM][CORNERS_OF_EDGE]    = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};
INT CornerOfOppEdge[MAX_EDGES_OF_ELEM][CORNERS_OF_EDGE] = {{2,3},{0,3},{3,1},{1,2},{2,0},{0,1}};


/* the indices of the edges between two corners */
INT EdgeWithCorners[MAX_CORNERS_OF_ELEM][MAX_CORNERS_OF_ELEM] = {{-1,0,2,3},{0,-1,1,4},
                                                                 {2,1,-1,5},{3,4,5,-1}};

/* the indices of the sides around each edge */
INT SideWithEdge[MAX_EDGES_OF_ELEM][MAX_SIDES_OF_EDGE] = {{0,3},{0,1},{0,2},{2,3},{1,3},{1,2}};

/* the index of the side having to specified edges */
INT SideOf2Edges[MAX_EDGES_OF_ELEM][MAX_EDGES_OF_ELEM] = {{-1,0,0,3,3,-1},{0,-1,0,-1,1,1},{0,0,-1,2,-1,2},{3,-1,2,-1,3,2},{3,1,-1,3,-1,1},{-1,1,2,2,1,-1}};

/* the index of the edge having to specified sides */
INT EdgeOf2Sides[MAX_SIDES_OF_ELEM][MAX_SIDES_OF_ELEM] = {{-1,1,2,0},{1,-1,5,4},{2,5,-1,3},{0,4,3,-1}};

/* the indices of the edges of each side */
INT EdgeOfSide[MAX_SIDES_OF_ELEM][EDGES_OF_TETRASIDE] = {{0,1,2},{1,5,4},{2,3,5},{0,4,3}};
INT CondensedEdgeOfSide[MAX_SIDES_OF_ELEM] = {0x07,0x32,0x2C,0x19};


/* the indices of opposite corners for each side */
INT OppositeCorner[MAX_SIDES_OF_ELEM] = {3,0,1,2};


/* the indices of opposite sides for each corner */
INT OppositeSide[MAX_CORNERS_OF_ELEM] = {1,2,3,0};

/* determine from NodeOrder (0,24( the order of nodes of a tetrahedron :			*/
/*	   Index=0 means the one with the largest Z-component in Cutting ref system     */
/*	   Index=3 means the one with the smallest Z-component in Cutting ref system	*/
INT CornerIndex[24][MAX_CORNERS_OF_ELEM] = {{0,1,2,3},{0,1,3,2},{0,2,1,3},{0,2,3,1},{0,3,1,2},{0,3,2,1},
                                            {1,0,2,3},{1,0,3,2},{1,2,0,3},{1,2,3,0},{1,3,0,2},{1,3,2,0},
                                            {2,0,1,3},{2,0,3,1},{2,1,0,3},{2,1,3,0},{2,3,0,1},{2,3,1,0},
                                            {3,0,1,2},{3,0,2,1},{3,1,0,2},{3,1,2,0},{3,2,0,1},{3,2,1,0}  };

/* determine from NodeOrder (0,24( the order of nodes of a side of a tetrahedron :	*/
/*	   Index=0 means the one with the largest Z-component in Cutting ref system     */
/*	   Index=2 means the one with the smallest Z-component in Cutting ref system	*/
INT SideCornerIndex[MAX_SIDES_OF_ELEM][24][MAX_CORNERS_OF_SIDE] =
{
  { {0,1,2},{0,1,2},{0,2,1},{0,2,1},{0,1,2},{0,2,1},
    {1,0,2},{1,0,2},{1,2,0},{1,2,0},{1,0,2},{1,2,0},
    {2,0,1},{2,0,1},{2,1,0},{2,1,0},{2,0,1},{2,1,0},
    {0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0} },
  { {1,2,3},{1,3,2},{2,1,3},{2,3,1},{3,1,2},{3,2,1},
    {1,2,3},{1,3,2},{1,2,3},{1,2,3},{1,3,2},{1,3,2},
    {2,1,3},{2,3,1},{2,1,3},{2,1,3},{2,3,1},{2,3,1},
    {3,1,2},{3,2,1},{3,1,2},{3,1,2},{3,2,1},{3,2,1} },
  { {0,2,3},{0,3,2},{0,2,3},{0,2,3},{0,3,2},{0,3,2},
    {0,2,3},{0,3,2},{2,0,3},{2,3,0},{3,0,2},{3,2,0},
    {2,0,3},{2,0,3},{2,0,3},{2,3,0},{2,3,0},{2,3,0},
    {3,0,2},{3,0,2},{3,0,2},{3,2,0},{3,2,0},{3,2,0} },
  { {0,1,3},{0,1,3},{0,1,3},{0,3,1},{0,3,1},{0,3,1},
    {1,0,3},{1,0,3},{1,0,3},{1,3,0},{1,3,0},{1,3,0},
    {0,1,3},{0,3,1},{1,0,3},{1,3,0},{3,0,1},{3,1,0},
    {3,0,1},{3,0,1},{3,1,0},{3,1,0},{3,0,1},{3,1,0} }
};

/* determine number of edge from reduced (i.e. restricted to one side) edgepattern */
/* if there are two edges marked for bisection, if not deliver -1. If the edge-    */
/* is not reduced (i.e. marked edges lying on more than one side) deliver -2	   */
INT TriSectionEdge[64][2] = {  {-1,-1},{-1,-1},{-1,-1},{ 1, 0},{-1,-1},{ 0, 2},{ 2, 1},{-1,-1},
                               {-1,-1},{ 3, 0},{-2,-2},{-2,-2},{ 2, 3},{-2,-2},{-2,-2},{-2,-2},
                               {-1,-1},{ 0, 4},{ 4, 1},{-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2},
                               { 4, 3},{-1,-1},{-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2},
                               {-1,-1},{-2,-2},{ 1, 5},{-2,-2},{ 5, 2},{-2,-2},{-2,-2},{-2,-2},
                               { 3, 5},{-2,-2},{-2,-2},{-2,-2},{-1,-1},{-2,-2},{-2,-2},{-2,-2},
                               { 5, 4},{-2,-2},{-1,-1},{-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2},
                               {-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2}      };

/* midnode of edge */
INT MidNodeOfEdge[MAX_EDGES_OF_ELEM] = {4,5,6,7,8,9};

/* the unique disjoint edge */
INT OppositeEdge[MAX_EDGES_OF_ELEM] = {5,3,4,1,2,0};

/* two times [second entry] two [third entry] edges for spec. edge [first entry] which bild up a side */
INT SideEdgesOfEdge[MAX_EDGES_OF_ELEM][2][2] = {{{1,2},{3,4}},
                                                {{0,2},{4,5}},
                                                {{0,1},{3,5}},
                                                {{2,5},{0,4}},
                                                {{1,5},{0,3}},
                                                {{1,4},{2,3}}};

/* edges with a given corner */
INT EdgesOfCorner[MAX_CORNERS_OF_ELEM][3] = {{0,2,3},{0,1,4},{1,2,5},{3,4,5}};

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* data for CVS */
static char rcsid[] = "$Header$";

/***************************************************************************
*																		   *
*							 Forward functions							   *
*																		   *
***************************************************************************/


/***************************************************************************
*																		   *
*							  Implementation							   *
*																		   *
***************************************************************************/

/***************************************************************************
*																		   *
*								 The End								   *
*																		   *
***************************************************************************/
