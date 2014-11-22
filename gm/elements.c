// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      elements.c                                                    */
/*																			*/
/* Purpose:   implements a general element concept							*/
/*																			*/
/* Author:	  Peter Bastian, Stefan Lang                                                                    */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de							    */
/*																			*/
/* History:   24.03.95 begin, ug version 3.0                                */
/*			  18.03.96 ug3.1                                                                */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#include <config.h>
#include <assert.h>

#include "general.h"
#include "ugdevices.h"

#include "gm.h"
#include "ugm.h"

#ifdef ModelP
#include "parallel.h"
#endif

#include "elements.h"

USING_UG_NAMESPACE
USING_UGDIM_NAMESPACE

/****************************************************************************/
/*													*/
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

INT NS_DIM_PREFIX n_offset[TAGS];
INT NS_DIM_PREFIX father_offset[TAGS];
INT NS_DIM_PREFIX sons_offset[TAGS];
INT NS_DIM_PREFIX nb_offset[TAGS];
INT NS_DIM_PREFIX evector_offset[TAGS];
INT NS_DIM_PREFIX svector_offset[TAGS];
INT NS_DIM_PREFIX side_offset[TAGS];
INT NS_DIM_PREFIX data_offset[TAGS];

GENERAL_ELEMENT * NS_DIM_PREFIX element_descriptors[TAGS];
GENERAL_ELEMENT * NS_DIM_PREFIX reference_descriptors[MAX_CORNERS_OF_ELEM+1];
INT NS_DIM_PREFIX reference2tag[MAX_CORNERS_OF_ELEM+1];

#ifndef ModelP
static INT nOBJT, OBJT4Elements[MAXOBJECTS];
#endif

/****************************************************************************/
/*																			*/
/* definition of local variables											*/
/*																			*/
/****************************************************************************/

#ifdef __TWODIM__
static GENERAL_ELEMENT def_triangle = {
  3,                                                                                    /* tag							*/
  4,                                                                                    /* max number of sons			*/
  3,                                                                                    /* number of sides				*/
  3,                                                                                    /* number of corners			*/
  {{0.0,0.0},{1.0,0.0},{0.0,1.0}},                      /* local coordinates			*/
  3,                                                                                    /* number of edges				*/
  {1,1,1,-1},                                                                   /* edges for each side	(2D!)	*/
  {2,2,2,-1},                                                                   /* corners for each side		*/
  2,                                                                                    /* an edge has 2 corners		*/
  {{0,-1,-1},{1,-1,-1},{2,-1,-1},{-1,-1,-1}},       /* number of edge j of side i   */
  {{0,1,-1},{1,2,-1},{2,0,-1},{-1,-1,-1}},          /* number of corner j of side i */
  {{0,1},{1,2},{2,0},{-1,-1},{-1,-1},{-1,-1}}       /* number of corner j of edge i */
} ;

static GENERAL_ELEMENT def_quadrilateral = {
  4,                                                                                    /* tag							*/
  4,                                                                                    /* max number of sons			*/
  4,                                                                                    /* number of sides				*/
  4,                                                                                    /* number of corners			*/
  {{0.0,0.0},{1.0,0.0},{1.0,1.0},
   {0.0,1.0}},                                                          /* local coordinates			*/
  4,                                                                                    /* number of edges				*/
  {1,1,1,1},                                                                    /* edges for each side	(2D!)	*/
  {2,2,2,2},                                                                    /* corners for each side		*/
  2,                                                                                    /* an edge has 2 corners		*/
  {{0,-1,-1},{1,-1,-1},{2,-1,-1},{3,-1,-1}},       /* number of edge j of side i   */
  {{0,1,-1},{1,2,-1},{2,3,-1},{3,0,-1}},           /* number of corner j of side i */
  {{0,1},{1,2},{2,3},{3,0},{-1,-1},{-1,-1}}        /* number of corner j of edge i */
} ;
#endif

#ifdef __THREEDIM__
static GENERAL_ELEMENT def_tetrahedron = {
  4,                                                                                    /* tag							*/
  12,                                                                                   /* max number of sons			*/
  4,                                                                                    /* number of sides				*/
  4,                                                                                    /* number of corners			*/
  {{0.0,0.0,0.0},{1.0,0.0,0.0},
   {0.0,1.0,0.0},{0.0,0.0,1.0}},                        /* local coordinates			*/
  6,                                                                                    /* number of edges				*/
  {3,3,3,3,-1,-1},                                                      /* edges for each side	(2D!)	*/
  {3,3,3,3,-1,-1},                                                      /* corners for each side		*/
  2,                                                                                    /* an edge has 2 corners		*/
  {{2,1,0,-1},{1,5,4,-1},{3,5,2,-1},{0,4,3,-1}},       /* number of edge j of side i   */
  {{0,2,1,-1},{1,2,3,-1},{0,3,2,-1},{0,1,3,-1}},       /* number of corner j of side i */
  {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3} }        /* number of corner j of edge i */
} ;

static GENERAL_ELEMENT def_pyramid = {
  5,                                                                                    /* tag							*/
  0,                                                                                    /* max number of sons			*/
  5,                                                                                    /* number of sides				*/
  5,                                                                                    /* number of corners			*/
  {{0.0,0.0,0.0},{1.0,0.0,0.0},{1.0,1.0,0.0},
   {0.0,1.0,0.0},{0.0,0.0,1.0}},                        /* local coordinates			*/
  8,                                                                                    /* number of edges				*/
  {4,3,3,3,3,-1},                                                       /* edges for each side	(2D!)	*/
  {4,3,3,3,3,-1},                                                       /* corners for each side		*/
  2,                                                                                    /* an edge has 2 corners		*/
  {{3,2,1,0},{0,5,4,-1},{1,6,5,-1},                     /* number of edge j of side i   */
   {2,7,6,-1},{3,4,7,-1}},
  {{0,3,2,1},{0,1,4,-1},{1,2,4,-1},                     /* number of corner j of side i */
   {2,3,4,-1},{3,0,4,-1}},
  {{0,1},{1,2},{2,3},{3,0},{0,4},{1,4},         /* number of corner j of edge i */
   {2,4},{3,4}}
} ;

static GENERAL_ELEMENT def_prism = {
  6,                                                                                    /* tag							*/
  0,                                                                                    /* max number of sons			*/
  5,                                                                                    /* number of sides				*/
  6,                                                                                    /* number of corners			*/
  {{0.0,0.0,0.0},{1.0,0.0,0.0},{0.0,1.0,0.0},
   {0.0,0.0,1.0},{1.0,0.0,1.0},{0.0,1.0,1.0}},       /* local coordinates		*/
  9,                                                                                    /* number of edges				*/
  {3,4,4,4,3,-1},                                                       /* edges for each side	(2D!)	*/
  {3,4,4,4,3,-1},                                                       /* corners for each side		*/
  2,                                                                                    /* an edge has 2 corners		*/
  {{2,1,0,-1},{0,4,6,3},{1,5,7,4},                      /* number of edge j of side i   */
   {2,3,8,5},{6,7,8,-1}},
  {{0,2,1,-1},{0,1,4,3},{1,2,5,4},                      /* number of corner j of side i */
   {2,0,3,5},{3,4,5,-1}},
  {{0,1},{1,2},{2,0},{0,3},{1,4},{2,5},         /* number of corner j of edge i */
   {3,4},{4,5},{5,3}}
} ;

static GENERAL_ELEMENT def_hexahedron = {
  7,                                                                                    /* tag							*/
  30,                                                                                   /* max number of sons			*/
  6,                                                                                    /* number of sides				*/
  8,                                                                                    /* number of corners			*/
  {{0.0,0.0,0.0},{1.0,0.0,0.0},
   {1.0,1.0,0.0},{0.0,1.0,0.0},
   {0.0,0.0,1.0},{1.0,0.0,1.0},
   {1.0,1.0,1.0},{0.0,1.0,1.0}},                        /* local coordinates			*/
  12,                                                                                   /* number of edges				*/
  {4,4,4,4,4,4},                                                        /* edges for each side	(2D!)	*/
  {4,4,4,4,4,4},                                                        /* corners for each side		*/
  2,                                                                                    /* an edge has 2 corners		*/
  {{3,2,1,0},{0,5,8,4},{1,6,9,5},                       /* number of edge j of side i   */
   {2,7,10,6},{3,4,11,7},{8,9,10,11}},
  {{0,3,2,1},{0,1,5,4},{1,2,6,5},                       /* number of corner j of side i */
   {2,3,7,6},{3,0,4,7},{4,5,6,7}},
  {{0,1},{1,2},{2,3},{3,0},{0,4},{1,5},         /* number of corner j of edge i */
   {2,6},{3,7},{4,5},{5,6},{6,7},{7,4}}
} ;
#endif

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/** \brief Compute index fields for a given element type

   \param el pointer to an element description

   STRUCTURES:
   \verbatim
   typedef struct {
    INT tag;                                // element type to be defined

    // the following parameters determine size of refs array in element
    INT max_sons_of_elem;                   // max number of sons for this type
    INT sides_of_elem;                      // how many sides ?
    INT corners_of_elem;                    // how many corners ?

    // more size parameters
    INT edges_of_elem;                      // how many edges ?
    INT edges_of_side[MAX_SIDES_OF_ELEM];   // number of edges for each side
    INT corners_of_side[MAX_SIDES_OF_ELEM]; // number of corners for each side
    INT corners_of_edge;                    // is always 2 !

    // index computations
    // Within each element sides, edges, corners are numbered in some way.
    // Within each side the edges and corners are numbered, within the edge the
    // corners are numbered. The following arrays map the local numbers within
    // the side or edge to the numbering within the element.
    INT edge_of_side[MAX_SIDES_OF_ELEM][MAX_EDGES_OF_SIDE];
    INT corner_of_side[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_SIDE];
    INT corner_of_edge[MAX_EDGES_OF_ELEM][MAX_CORNERS_OF_EDGE];

    // the following parameters are derived from data above
    INT mapped_inner_objt;                  // tag to objt mapping for free list
    INT mapped_bnd_objt;                    // tag to objt mapping for free list
    INT inner_size, bnd_size;               // size in bytes used for alloc
    INT edge_with_corners[MAX_CORNERS_OF_ELEM][MAX_CORNERS_OF_ELEM];
    INT side_with_edge[MAX_EDGES_OF_ELEM][MAX_SIDES_OF_EDGE];
    INT corner_of_side_inv[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_ELEM];
    INT edges_of_corner[MAX_CORNERS_OF_ELEM][MAX_EDGES_OF_ELEM];
    INT corner_of_oppedge[MAX_EDGES_OF_ELEM][MAX_CORNERS_OF_EDGE];
    INT corner_opp_to_side[MAX_SIDES_OF_ELEM];
    INT opposite_edge[MAX_EDGES_OF_ELEM];
    INT side_opp_to_corner[MAX_CORNERS_OF_ELEM];
    INT edge_of_corner[MAX_CORNERS_OF_ELEM][MAX_EDGES_OF_ELEM];
    INT edge_of_two_sides[MAX_SIDES_OF_ELEM][MAX_SIDES_OF_ELEM];

   } GENERAL_ELEMENT;
   \endverbatim


   This function processes a topology description of an element type and computes
   index mappings. Currently descriptions for triangles,
   quadrilaterals and tetrahedra are included. Hexahedral elements have been implemented
   in a prototype version.

   CAUTION: The above data structure is filled UP TO appropriate sizes for memory allocation
   as well as offsets in the 'refs' array of the 'generic_element'! For complete filling
   you will have to call 'ProcessElementDescription'

   Only the following components of the GENERAL_ELEMENT structure must be provided.
   All other components are derived from the given information.

   . tag - New tag for the elememt which will be delivered by the 'TAG' macro.
   . max_sons_of_elem - Max number of sons allowed for that element type.
   . sides_of_elem - Number of sides for that element type.
   . corners_of_elem - Number of corners for that element type.
   . edges_of_elem - Number of edges for that element type.
   . edges_of_side - Number of edges for each side.
   . corners_of_side - Number of corners of each side.
   . corners_of_edge - is always 2.
   . edge_of_side[s][e] - The edges are numbered in the element and in each side of the
   element. This array provides a mapping that tells you the number of edge 'e' of side
   's' with respect to the numbering in the element.
   . corner_of_side[s][c] - The corners edges are numbered in the element and in each side of the
   element. This array provides a mapping that tells you the number of corner 'c' of side
   's' with respect to the numbering in the element.
   . corner_of_edge[e][c] - Tells you the number of corner 'c' in edge 'e' with respect
   to the numbering in the element.

   SEE ALSO:

   'ELEMENT', 'ProcessElementDescription'.

   \return <ul>
   <li> GM_OK if ok </li>
   <li> GM_ERROR if error occured </li>
   <ul>
 */
/****************************************************************************/

static INT PreProcessElementDescription (GENERAL_ELEMENT *el)
{
  INT tag;
  INT i,j,k,l,n,from,to;
        #ifdef __THREEDIM__
  INT m,n1,n2;
        #endif

  tag = el->tag;

  /* derive additional index fields */

  /* edge_with_corners(i,j) : number of edge between corners i and j, -1 if no such edge ex. */
  for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    for (j=0; j<MAX_CORNERS_OF_ELEM; j++) el->edge_with_corners[i][j] = -1;
  for (i=0; i<el->edges_of_elem; i++)
  {
    el->edge_with_corners[el->corner_of_edge[i][0]][el->corner_of_edge[i][1]] = i;
    el->edge_with_corners[el->corner_of_edge[i][1]][el->corner_of_edge[i][0]] = i;
  }

  /* side_with_edge(i,j) : edge i is an edge of side side_with_edge(i,j) */
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
    for (j=0; j<MAX_SIDES_OF_EDGE; j++) el->side_with_edge[i][j] = -1;
  for (k=0; k<el->edges_of_elem; k++)
  {
    from = el->corner_of_edge[k][0];
    to   = el->corner_of_edge[k][1];

    for (i=0; i<el->sides_of_elem; i++) {
      n = el->corners_of_side[i];
      for (j=0; j<n; j++)
      {
        if ((el->corner_of_side[i][j]==from)&&(el->corner_of_side[i][(j+1)%n]==to))
          el->side_with_edge[k][1] = i;
        if ((el->corner_of_side[i][j]==to)&&(el->corner_of_side[i][(j+1)%n]==from))
          el->side_with_edge[k][0] = i;
      }
    }
  }

  /* corner_of_side_inv(i,j) : j is number of a corner in the element. Then this
     array returns the local number of this corner within side i or -1 if side i
     does not contain this corner. */
  for (i=0; i<MAX_SIDES_OF_ELEM; i++)
    for (j=0; j<MAX_CORNERS_OF_ELEM; j++) el->corner_of_side_inv[i][j] = -1;
  for (i=0; i<el->sides_of_elem; i++)
    for (j=0; j<el->corners_of_side[i]; j++)
    {
      n = el->corner_of_side[i][j];
      el->corner_of_side_inv[i][n] = j;
    }

  /* edges_of_corner(i,j) : i is a number of a corner within the element, then
     edges_of_corner(i,j) gives the number of an edge adjacent to corner i or -1 */
  for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    for (j=0; j<MAX_EDGES_OF_ELEM; j++) el->edges_of_corner[i][j] = -1;
  for (i=0; i<el->edges_of_elem; i++)
    for (j=0; j<el->corners_of_edge; j++)
    {
      n = el->corner_of_edge[i][j];
      for (k=0; k<MAX_EDGES_OF_ELEM; k++)
        if (el->edges_of_corner[n][k]<0)
        {
          el->edges_of_corner[n][k] = i;
          break;
        }
    }


  /* fields not valid for all elements */

  /* corner_of_oppedge(i,j) */
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
    for (j=0; j<MAX_CORNERS_OF_EDGE; j++)
      el->corner_of_oppedge[i][j] = -1;

  /* corner_opp_to_side(i) */
  for (i=0; i<MAX_SIDES_OF_ELEM; i++)
    el->corner_opp_to_side[i] = -1;

  /* opposite_edge(i) */
  for (i=0; i<MAX_EDGES_OF_ELEM; i++)
    el->opposite_edge[i] = -1;

  /* side_opp_to_corner(i) */
  for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    el->side_opp_to_corner[i] = -1;

  /* edge_of_corner(i,j) */
  for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    for (j=0; j<MAX_EDGES_OF_ELEM; j++)
      el->edge_of_corner[i][j] = -1;

#ifdef __TWODIM__
  switch (tag)
  {
  case TRIANGLE :
    /* corner_of_oppedge(i,j) */
    /* is not defined!		  */

    /* corner_opp_to_side(i)  */
    /* is not defined!		  */

    /* opposite_edge(i)		  */
    /* is not defined!		  */

    /* side_opp_to_corner(i)  */
    /* is not defined!		  */

    /* edge_of_corner(i,j)	  */
    for (i=0; i<el->edges_of_elem; i++) {
      for (j=0; j<el->corners_of_edge; j++) {
        if (el->corner_of_edge[i][j] >=0) {
          for (k=0; k<el->edges_of_elem; k++)
            if (el->edge_of_corner[el->corner_of_edge[i][j]][k] < 0)
              break;
          assert(k<el->edges_of_elem);
          el->edge_of_corner[el->corner_of_edge[i][j]][k] = i;
        }
      }
    }

    break;

  case QUADRILATERAL :
    /* corner_of_oppedge(i,j) */
    for (i=0; i<el->edges_of_elem; i++) {
      for (j=0; j<el->edges_of_elem; j++) {
        n=1;
        for (k=0; k<el->corners_of_edge; k++)
          for (l=0; l<el->corners_of_edge; l++)
            if (el->corner_of_edge[i][k]==el->corner_of_edge[j][l])
              n=0;
        if (n) {
          el->corner_of_oppedge[i][0] = el->corner_of_edge[j][0];
          el->corner_of_oppedge[i][1] = el->corner_of_edge[j][1];
          break;
        }
      }
      assert(j<el->edges_of_elem);
    }

    /* corner_opp_to_side(i)  */
    /* is not defined!		  */

    /* opposite_edge(i)		  */
    for (i=0; i<el->edges_of_elem; i++) {
      n = 0;
      for (j=0; j<el->corners_of_edge; j++) {
        for (k=0; k<el->edges_of_elem; k++) {
          if (el->edges_of_corner[el->corner_of_edge[i][j]][k] >= 0)
            n |= (0x1<<(el->edges_of_corner[el->corner_of_edge[i][j]][k]));
        }
      }
      for (j=0; j<el->edges_of_elem; j++)
        if (((n>>j) & 0x1) == 0)
          break;
      assert(j<el->edges_of_elem);

      el->opposite_edge[i] = j;
    }

    /* side_opp_to_corner(i)  */
    /* is not defined!		  */

    /* edge_of_corner(i,j)	  */
    for (i=0; i<el->edges_of_elem; i++) {
      for (j=0; j<el->corners_of_edge; j++) {
        if (el->corner_of_edge[i][j] >=0) {
          for (k=0; k<el->edges_of_elem; k++)
            if (el->edge_of_corner[el->corner_of_edge[i][j]][k] < 0)
              break;
          assert(k<el->edges_of_elem);
          el->edge_of_corner[el->corner_of_edge[i][j]][k] = i;
        }
      }
    }

    break;
  }
#endif

#ifdef __THREEDIM__
  /* edge_of_two_sides(i,j) */
  for (i=0; i<MAX_SIDES_OF_ELEM; i++)
    for (j=0; j<MAX_SIDES_OF_ELEM; j++)
      el->edge_of_two_sides[i][j] = -1;

  switch (tag)
  {
  case TETRAHEDRON :

    /* corner_of_oppedge(i,j) */
    for (i=0; i<el->edges_of_elem; i++) {
      for (j=0; j<el->edges_of_elem; j++) {
        n=1;
        for (k=0; k<el->corners_of_edge; k++)
          for (l=0; l<el->corners_of_edge; l++)
            if (el->corner_of_edge[i][k]==el->corner_of_edge[j][l])
              n=0;
        if (n) {
          el->corner_of_oppedge[i][0] = el->corner_of_edge[j][0];
          el->corner_of_oppedge[i][1] = el->corner_of_edge[j][1];
          break;
        }
      }
      assert(j<el->edges_of_elem);
    }

    /* corner_opp_to_side(i) */
    for (i=0; i<el->sides_of_elem; i++) {
      n = 0;
      for (j=0; j<el->corners_of_side[i]; j++) {
        n |= (0x1<<(el->corner_of_side[i][j]));
      }
      for (j=0; j<el->corners_of_elem; j++) {
        if (((n>>j) & 0x1) == 0)
          break;
      }
      assert(j<el->corners_of_elem);
      el->corner_opp_to_side[i] = j;
    }

    /* opposite_edge(i)		  */
    for (i=0; i<el->edges_of_elem; i++) {
      n = 0;
      for (j=0; j<el->corners_of_edge; j++) {
        for (k=0; k<el->edges_of_elem; k++) {
          if (el->edges_of_corner[el->corner_of_edge[i][j]][k] >= 0)
            n |= (0x1<<(el->edges_of_corner[el->corner_of_edge[i][j]][k]));
        }
      }
      for (j=0; j<el->edges_of_elem; j++)
        if (((n>>j) & 0x1) == 0)
          break;
      assert(j<el->edges_of_elem);

      el->opposite_edge[i] = j;
    }

    /* side_opp_to_corner(i)  */
    for (i=0; i<el->corners_of_elem; i++) {
      for (j=0; j<el->sides_of_elem; j++) {
        n = 0;
        for (k=0; k<el->corners_of_side[j]; k++)
          n |= (0x1<<(el->corner_of_side[j][k]));
        if (((n>>i) & 0x1) == 0) {
          el->side_opp_to_corner[i] = j;
          break;
        }
      }
      assert(j<el->sides_of_elem);
    }

    /* edge_of_corner(i,j)	  */
    for (i=0; i<el->edges_of_elem; i++) {
      for (j=0; j<el->corners_of_edge; j++) {
        if (el->corner_of_edge[i][j] >=0) {
          for (k=0; k<el->edges_of_elem; k++)
            if (el->edge_of_corner[el->corner_of_edge[i][j]][k] < 0)
              break;
          assert(k<el->edges_of_elem);
          el->edge_of_corner[el->corner_of_edge[i][j]][k] = i;
        }
      }
    }

    break;

  case PYRAMID :

    /* corner_of_oppedge(i,j) */
    /* is not defined!		  */

    /* corner_opp_to_side(i) */
    for (i=0; i<el->sides_of_elem; i++) {
      if (el->corners_of_side[i] == 4) {
        n = 0;
        for (j=0; j<el->corners_of_side[i]; j++) {
          n |= (0x1<<(el->corner_of_side[i][j]));
        }
        for (j=0; j<el->corners_of_elem; j++) {
          if (((n>>j) & 0x1) == 0)
            break;
        }
        assert(j<el->corners_of_elem);
        el->corner_opp_to_side[i] = j;
      }
    }

    /* opposite_edge(i)		  */
    /* is not defined!		  */

    /* side_opp_to_corner(i)  */
    for (i=0; i<el->corners_of_elem; i++) {
      for (j=0; j<el->sides_of_elem; j++) {
        n = 0;
        for (k=0; k<el->corners_of_side[j]; k++)
          n |= (0x1<<(el->corner_of_side[j][k]));
        if (((n>>i) & 0x1) == 0) {
          el->side_opp_to_corner[i] = j;
          break;
        }
      }
      assert(j<el->sides_of_elem);
    }

    /* edge_of_corner(i,j)	  */
    for (i=0; i<el->edges_of_elem; i++) {
      for (j=0; j<el->corners_of_edge; j++) {
        if (el->corner_of_edge[i][j] >=0) {
          for (k=0; k<el->edges_of_elem; k++)
            if (el->edge_of_corner[el->corner_of_edge[i][j]][k] < 0)
              break;
          assert(k<el->edges_of_elem);
          el->edge_of_corner[el->corner_of_edge[i][j]][k] = i;
        }
      }
    }

    break;

  case PRISM :

    /* corner_of_oppedge(i,j) */
    /* is not defined!		  */

    /* corner_opp_to_side(i)  */
    /* is not defined!		  */

    /* opposite_edge(i)		  */
    /* is not defined!		  */

    /* side_opp_to_corner(i)  */
    /* is not defined!		  */

    /* edge_of_corner(i,j)	  */
    for (i=0; i<el->edges_of_elem; i++) {
      for (j=0; j<el->corners_of_edge; j++) {
        if (el->corner_of_edge[i][j] >=0) {
          for (k=0; k<el->edges_of_elem; k++)
            if (el->edge_of_corner[el->corner_of_edge[i][j]][k] < 0)
              break;
          assert(k<el->edges_of_elem);
          el->edge_of_corner[el->corner_of_edge[i][j]][k] = i;
        }
      }
    }

    break;

  case HEXAHEDRON :

    /* corner_of_oppedge(i,j) */
    for (i=0; i<el->edges_of_elem; i++) {
      n = 0;
      for (j=0; j<el->corners_of_edge; j++) {
        n1 = el->corner_of_edge[i][j];
        for (k=0; k<el->edges_of_elem; k++) {
          if (el->edges_of_corner[n1][k] >= 0) {
            n |= (0x1<<(el->edges_of_corner[n1][k]));
            for (l=0; l<el->corners_of_edge; l++) {
              n2 = el->corner_of_edge[el->edges_of_corner[n1][k]][l];
              if (n2 != n1) {
                for (m=0; m<el->edges_of_elem; m++)
                  if (el->edges_of_corner[n2][m] >= 0)
                    n |= (0x1<<(el->edges_of_corner[n2][m]));
              }
            }
          }
        }
      }
      for (k=0; k<el->edges_of_elem; k++)
        if (((n>>k) & 0x1) == 0)
          break;
      assert(k<el->edges_of_elem);

      el->corner_of_oppedge[i][0] = el->corner_of_edge[k][0];
      el->corner_of_oppedge[i][1] = el->corner_of_edge[k][1];
    }

    /* corner_opp_to_side(i)  */
    /* is not defined!		  */

    /* opposite_edge(i)		  */
    for (i=0; i<el->edges_of_elem; i++) {
      n = 0;
      for (j=0; j<el->corners_of_edge; j++) {
        n1 = el->corner_of_edge[i][j];
        for (k=0; k<el->edges_of_elem; k++) {
          if (el->edges_of_corner[n1][k] >= 0) {
            n |= (0x1<<(el->edges_of_corner[n1][k]));
            for (l=0; l<el->corners_of_edge; l++) {
              n2 = el->corner_of_edge[el->edges_of_corner[n1][k]][l];
              if (n2 != n1) {
                for (m=0; m<el->edges_of_elem; m++) {
                  if (el->edges_of_corner[n2][m] >= 0)
                    n |= (0x1<<(el->edges_of_corner[n2][m]));
                }
              }
            }
          }
        }
      }
      for (k=0; k<el->edges_of_elem; k++)
        if (((n>>k) & 0x1) == 0)
          break;
      assert(k<el->edges_of_elem);

      el->opposite_edge[i] = k;
    }

    /* side_opp_to_corner(i)  */
    /* is not defined!		  */

    /* edge_of_corner(i,j)	  */
    for (i=0; i<el->edges_of_elem; i++) {
      for (j=0; j<el->corners_of_edge; j++) {
        if (el->corner_of_edge[i][j] >=0) {
          for (k=0; k<el->edges_of_elem; k++)
            if (el->edge_of_corner[el->corner_of_edge[i][j]][k] < 0)
              break;
          assert(k<el->edges_of_elem);
          el->edge_of_corner[el->corner_of_edge[i][j]][k] = i;
        }
      }
    }

    break;
  }

  for (i=0; i<el->sides_of_elem; i++)
    for (j=0; j<el->sides_of_elem; j++)
      for (k=0; k<el->edges_of_side[i]; k++)
        for (l=0; l<el->edges_of_side[j]; l++)
          if (el->edge_of_side[i][k] == el->edge_of_side[j][l])
          {
            assert( (i==j) || (el->edge_of_two_sides[i][j]==-1) ||
                    (el->edge_of_two_sides[i][j]==el->edge_of_side[i][k]) );
            el->edge_of_two_sides[i][j] = el->edge_of_side[i][k];
          }
#endif

  /* make description globally available */
  element_descriptors[tag] = el;
  reference_descriptors[el->corners_of_elem] = el;
  reference2tag[el->corners_of_elem] = tag;

  return(GM_OK);
}

/****************************************************************************/
/** \brief Compute offsets and size for a given element type

   \param theMG multigrid for format dependent pointer offsets in elements
   \param el    pointer to an element description

   STRUCTURES:
   \verbatim
   typedef struct {
    INT tag;                                // element type to be defined

    // the following parameters determine size of refs array in element
    INT max_sons_of_elem;                   // max number of sons for this type
    INT sides_of_elem;                      // how many sides ?
    INT corners_of_elem;                    // how many corners ?

    // more size parameters
    INT edges_of_elem;                      // how many edges ?
    INT edges_of_side[MAX_SIDES_OF_ELEM];   // number of edges for each side
    INT corners_of_side[MAX_SIDES_OF_ELEM]; // number of corners for each side
    INT corners_of_edge;                    // is always 2 !

    // index computations
    // Within each element sides, edges, corners are numbered in some way.
    // Within each side the edges and corners are numbered, within the edge the
    // corners are numbered. The following arrays map the local numbers within
    // the side or edge to the numbering within the element.
    INT edge_of_side[MAX_SIDES_OF_ELEM][MAX_EDGES_OF_SIDE];
    INT corner_of_side[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_SIDE];
    INT corner_of_edge[MAX_EDGES_OF_ELEM][MAX_CORNERS_OF_EDGE];

    // the following parameters are derived from data above
    INT mapped_inner_objt;                  // tag to objt mapping for free list
    INT mapped_bnd_objt;                    // tag to objt mapping for free list
    INT inner_size, bnd_size;               // size in bytes used for alloc
    INT edge_with_corners[MAX_CORNERS_OF_ELEM][MAX_CORNERS_OF_ELEM];
    INT side_with_edge[MAX_EDGES_OF_ELEM][MAX_SIDES_OF_EDGE];
    INT corner_of_side_inv[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_ELEM];
    INT edges_of_corner[MAX_CORNERS_OF_ELEM][MAX_EDGES_OF_ELEM];
        INT corner_of_oppedge[MAX_EDGES_OF_ELEM][MAX_CORNERS_OF_EDGE];
        INT corner_opp_to_side[MAX_SIDES_OF_ELEM];
        INT opposite_edge[MAX_EDGES_OF_ELEM];
        INT side_opp_to_corner[MAX_CORNERS_OF_ELEM];
        INT edge_of_corner[MAX_CORNERS_OF_ELEM][MAX_EDGES_OF_ELEM];
        INT edge_of_two_sides[MAX_SIDES_OF_ELEM][MAX_SIDES_OF_ELEM];

   } GENERAL_ELEMENT;
   \endverbatim

   This function processes a topology description of an element type and computes
   the appropriate sizes for memory allocation and offsets in the 'refs' array of the
   'generic_element'. All other data are fixed and do not depend on the multigrid or format.
   Before calling this function 'PreProcessElementDescription' has to be called once
   during initialization.

   SEE ALSO:

   'ELEMENT', 'PreProcessElementDescription'.

   \return <ul>
   <li> GM_OK if ok </li>
   <li> GM_ERROR if error occured. </li>
   </ul>
 */
/****************************************************************************/

static INT ProcessElementDescription (MULTIGRID *theMG, GENERAL_ELEMENT *el)
{
  INT p_count, tag;

  tag = el->tag;
  p_count = 0;

  /* the corners */
  n_offset[tag] = p_count; p_count += el->corners_of_elem;

  /* the father */
  father_offset[tag] = p_count; p_count++;

  /* the sons */
  sons_offset[tag] = 0;
  /*
     #ifdef __TWODIM__
                  sons_offset[tag] = p_count; p_count += el->max_sons_of_elem;
     #endif
     #ifdef __THREEDIM__
                  sons_offset[tag] = p_count; p_count++;
     #endif
   */
  /* for 2D/3D one son pointer is stored in serial case     */
  /* for 2D/3D two son pointer are stored in parallel case: */
  /*    one to master elements, other to (h/v) ghosts       */
        #ifdef ModelP
  sons_offset[tag] = p_count; p_count+=2;
        #else
  sons_offset[tag] = p_count; p_count++;
        #endif

  /* the neighbors */
  nb_offset[tag] = p_count; p_count += el->sides_of_elem;

  /* element vector */
  evector_offset[tag] = 0;
  if (VEC_DEF_IN_OBJ_OF_MG(theMG,ELEMVEC))
  {
    evector_offset[tag] = p_count;
    p_count++;
  }

  /* side vector */
  svector_offset[tag] = 0;
        #ifdef __THREEDIM__
  if (VEC_DEF_IN_OBJ_OF_MG(theMG,SIDEVEC))
  {
    svector_offset[tag] = p_count;
    p_count += el->sides_of_elem;
  }
        #endif

  if (EDATA_DEF_IN_MG(theMG)) {
    data_offset[tag] = p_count;
    p_count++;
  }

  /* so far for an inner element */
  el->inner_size = sizeof(struct generic_element) + (p_count-1)*sizeof(void *);

  /* the element sides */
  side_offset[tag] = p_count; p_count += el->sides_of_elem;

  /* now the size of an element on the boundary */
  el->bnd_size = sizeof(struct generic_element) + (p_count-1)*sizeof(void *);

  /* get a free object id for free list */
  /** \todo OBJT is always allocated when this functions is called but never released
                   this will probably cause problems when several mgs are open: switching between
                   them will lead to an overflow of the UsedOBJT variable in ugm.c
                   possible remedy: store element OBJT in mg and release when it is closed. Also
                   don't reallocate them for a given mg */
  el->mapped_inner_objt = GetFreeOBJT();
  if (el->mapped_inner_objt < 0) return(GM_ERROR);

#ifndef ModelP
  if (nOBJT>=MAXOBJECTS-1) return(GM_ERROR);
  OBJT4Elements[nOBJT++]=el->mapped_inner_objt;
#endif

  el->mapped_bnd_objt = GetFreeOBJT();
  if (el->mapped_bnd_objt < 0) return(GM_ERROR);

#ifndef ModelP
  OBJT4Elements[nOBJT++]=el->mapped_bnd_objt;
  if (nOBJT>=MAXOBJECTS-1) return(GM_ERROR);
#endif

  return(GM_OK);
}

/****************************************************************************/
/** \brief Pre-initialize general element data up to multigrid dependent stuff

   This function pre-initializes the general element data up to multigrid dependent stuff.

   \return <ul>
   <li> GM_OK if ok </li>
   <li> GM_ERROR if error occured. </li>
   </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX PreInitElementTypes (void)
{
  INT err;

#ifdef __TWODIM__
  err = PreProcessElementDescription(&def_triangle);
  if (err!=GM_OK) return(err);
  err = PreProcessElementDescription(&def_quadrilateral);
  if (err!=GM_OK) return(err);
#endif

#ifdef __THREEDIM__
  err = PreProcessElementDescription(&def_tetrahedron);
  if (err!=GM_OK) return(err);
  err = PreProcessElementDescription(&def_pyramid);
  if (err!=GM_OK) return(err);
  err = PreProcessElementDescription(&def_prism);
  if (err!=GM_OK) return(err);
  err = PreProcessElementDescription(&def_hexahedron);
  if (err!=GM_OK) return(err);
#endif

  return (0);
}

/****************************************************************************/
/** \brief Initialize topological information for element types

   This function initializes topological information for element types and
   is called once during startup. Add your initialization of a new element
   type here.

   \return <ul>
   <li> GM_OK if ok </li>
   <li> GM_ERROR if error occured </li>
   </ul>
 */
/****************************************************************************/

INT NS_DIM_PREFIX InitElementTypes (MULTIGRID *theMG)
{
  INT i,err;

  if (theMG==NULL)
    return(GM_ERROR);

#ifndef ModelP
  /* release allocated OBJTs */
  for (i=0; i<nOBJT; i++)
    if (ReleaseOBJT (OBJT4Elements[i]))
      return (GM_ERROR);
  nOBJT=0;
#endif

#ifdef __TWODIM__
  err = ProcessElementDescription(theMG,&def_triangle);
  if (err!=GM_OK) return(err);
  err = ProcessElementDescription(theMG,&def_quadrilateral);
  if (err!=GM_OK) return(err);
#endif

#ifdef __THREEDIM__
  err = ProcessElementDescription(theMG,&def_tetrahedron);
  if (err!=GM_OK) return(err);
  err = ProcessElementDescription(theMG,&def_pyramid);
  if (err!=GM_OK) return(err);
  err = ProcessElementDescription(theMG,&def_prism);
  if (err!=GM_OK) return(err);
  err = ProcessElementDescription(theMG,&def_hexahedron);
  if (err!=GM_OK) return(err);
#endif

#ifdef ModelP
  InitCurrMG(theMG);
#endif

  return(GM_OK);
}
