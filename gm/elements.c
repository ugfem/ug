// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  elements.c													*/
/*																			*/
/* Purpose:   implements a general element concept							*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de							*/
/*																			*/
/* History:   24.03.95 begin, ug version 3.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/* define this to exclude extern definition of global arrays */
#define __COMPILE_EL__

#include "devices.h"

#include "switch.h"
#include "gm.h"
#include "ugm.h"

#include "elements.h"

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

INT n_offset[TAGS];
INT father_offset[TAGS];
INT sons_offset[TAGS];
INT nb_offset[TAGS];
INT evector_offset[TAGS];
INT svector_offset[TAGS];
INT data_offset[TAGS];
INT side_offset[TAGS];

GENERAL_ELEMENT *element_descriptors[TAGS];

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
  3,                                                                                    /* number of edges				*/
  {1,1,1,0},                                                                    /* edges for each side	(2D!)	*/
  {2,2,2,0},                                                                    /* corners for each side		*/
  2,                                                                                    /* an edge has 2 corners		*/
  { {0,0,0}, {1,0,0}, {2,0,0}, {0,0,0} },       /* number of edge j of side i   */
  { {0,1,0}, {1,2,0}, {2,0,0}, {0,0,0} },       /* number of corner j of side i */
  { {0,1},{1,2},{2,0},{0,0},{0,0},{0,0} }       /* number of corner j of edge i */
} ;

static GENERAL_ELEMENT def_quadrilateral = {
  4,                                                                                    /* tag							*/
  4,                                                                                    /* max number of sons			*/
  4,                                                                                    /* number of sides				*/
  4,                                                                                    /* number of corners			*/
  4,                                                                                    /* number of edges				*/
  {1,1,1,1},                                                                    /* edges for each side	(2D!)	*/
  {2,2,2,2},                                                                    /* corners for each side		*/
  2,                                                                                    /* an edge has 2 corners		*/
  { {0,0,0}, {1,0,0}, {2,0,0}, {3,0,0} },       /* number of edge j of side i   */
  { {0,1,0}, {1,2,0}, {2,3,0}, {3,0,0} },       /* number of corner j of side i */
  { {0,1},{1,2},{2,3},{3,0},{0,0},{0,0} }       /* number of corner j of edge i */
} ;
#endif

#ifdef __THREEDIM__
static GENERAL_ELEMENT def_tetrahedron = {
  4,                                                                                    /* tag							*/
  12,                                                                                   /* max number of sons			*/
  4,                                                                                    /* number of sides				*/
  4,                                                                                    /* number of corners			*/
  6,                                                                                    /* number of edges				*/
  {3,3,3,3},                                                                    /* edges for each side	(2D!)	*/
  {3,3,3,3},                                                                    /* corners for each side		*/
  2,                                                                                    /* an edge has 2 corners		*/
  { {0,1,2}, {1,5,4}, {2,3,5}, {0,4,3} },       /* number of edge j of side i   */
  { {0,2,1}, {1,2,3}, {2,0,3}, {3,0,1} },       /* number of corner j of side i */
  { {0,1},{1,2},{0,2},{0,3},{1,3},{2,3} }       /* number of corner j of edge i */
} ;
#endif

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*D
   ProcessElementDescription - compute offsets and size for a given element type

   SYNOPSIS:
   static INT ProcessElementDescription (GENERAL_ELEMENT *el);

   PARAMETERS:
   .  el - pointer to an element description

   STRUCTURES:

   .vb
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

   } GENERAL_ELEMENT;
   .ve

   DESCRIPTION:
   This function processes a topology description of an element type and computes
   the appropriated sizes for memory allocation, offsets in the 'refs' array of the
   'generic_element' and index mappings. Currently descriptions for triangles,
   quadrilaterals and tetrahedra are included. Hexahedral elements have been implemented
   in a prototype version.

   `Only the following components of the 'GENERAL_ELEMENT' structure must be provided.
   All other components are derived from the given information.`

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

   'ELEMENT'.

   RETURN VALUE:
   INT
   .n    GM_OK if ok
   .n    GM_ERROR if error occured.
   D*/
/****************************************************************************/

static INT ProcessElementDescription (GENERAL_ELEMENT *el)
{
  INT p_count, tag;
  INT i,j,k,n;

  tag = el->tag;
  p_count = 0;

  /* the corners */
  n_offset[tag] = p_count; p_count += el->corners_of_elem;

  /* the father */
  father_offset[tag] = p_count; p_count++;

  /* the sons */
  sons_offset[tag] = 0;
        #ifdef __TWODIM__
  sons_offset[tag] = p_count; p_count += el->max_sons_of_elem;
        #endif
        #ifdef __THREEDIM__
  sons_offset[tag] = p_count; p_count++;
        #endif

  /* the neighbors */
  nb_offset[tag] = p_count; p_count += el->sides_of_elem;

  /* element vector */
  evector_offset[tag] = 0;
        #ifdef __ELEMDATA__
  evector_offset[tag] = p_count; p_count++;
        #endif

  /* side vector */
  svector_offset[tag] = 0;
        #ifdef __SIDEDATA__
  svector_offset[tag] = p_count; p_count += el->sides_of_elem;
        #endif

  /* data vector in version 2.3 mode */
  data_offset[tag] = 0;
        #ifdef __version23__
  data_offset[tag] = p_count; p_count++;
        #endif

  /* so far for an inner element */
  el->inner_size = sizeof(struct generic_element) + (p_count-1)*sizeof(void *);

  /* the element sides */
  side_offset[tag] = p_count; p_count += el->sides_of_elem;

  /* now the size of an element on the boundary */
  el->bnd_size = sizeof(struct generic_element) + (p_count-1)*sizeof(void *);

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
  for (i=0; i<el->sides_of_elem; i++)
    for (j=0; j<el->edges_of_side[i]; j++)
    {
      n = el->edge_of_side[i][j];
      for (k=0; k<MAX_SIDES_OF_EDGE; k++)
        if (el->side_with_edge[n][k]<0)
        {
          el->side_with_edge[n][k] = i;
          break;
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

  /* make description globally available */
  element_descriptors[tag] = el;

  /* get a free object id for free list */
  el->mapped_inner_objt = GetFreeOBJT();
  el->mapped_bnd_objt = GetFreeOBJT();

  return(GM_OK);
}


/****************************************************************************/
/*D
   InitElementTypes - Initialize topological information for element types

   SYNOPSIS:
   INT InitElementTypes (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function initializes topological information for element types and
   is called once during startup. Add your initialization of a new element
   type here.

   RETURN VALUE:
   INT
   .n   GM_OK if ok
   .n   GM_ERROR if error occured.
   D*/
/****************************************************************************/

INT InitElementTypes (void)
{
  INT err;

#ifdef __TWODIM__
  err = ProcessElementDescription(&def_triangle);
  if (err!=GM_OK) return(err);
  err = ProcessElementDescription(&def_quadrilateral);
  if (err!=GM_OK) return(err);
#endif

#ifdef __THREEDIM__
  err = ProcessElementDescription(&def_tetrahedron);
  if (err!=GM_OK) return(err);
#endif

  return(GM_OK);
}
