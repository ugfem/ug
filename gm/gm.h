// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  gm.h															*/
/*																			*/
/* Purpose:   grid manager header file (the heart of ug)					*/
/*																			*/
/* Author:	  Peter Bastian, Klaus Johannsen								*/
/*																			*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/*			  blockvector data structure:									*/
/*			  Christian Wrobel                                                                              */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de					            */
/*																			*/
/* History:   09.03.92 begin, ug version 2.0  (as ugtypes2.h)				*/
/*			  13.12.94 begin, ug version 3.0								*/
/*			  27.09.95 blockvector implemented (Christian Wrobel)			*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/*																			*/
/****************************************************************************/

/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __GM__
#define __GM__

#include <limits.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __HEAPS__
#include "heaps.h"
#endif

#ifndef __UGENV__
#include "ugenv.h"
#endif

#ifndef __MISC__
#include "misc.h"
#endif

#ifndef __DEBUG__
#include "debug.h"
#endif

#ifndef __DOMAIN__
#include "domain.h"
#endif

#ifndef __PARGM_H__
#include "pargm.h"
#endif

/****************************************************************************/
/*																			*/
/* consistency of commandline defines										*/
/*																			*/
/****************************************************************************/

#if (!defined _2) && (!defined _3)
#error ****	define dimension _2 or _3		****
#endif

#ifdef _2
        #ifdef _3
        #error ****	define EITHER dimension _2 OR _3	   ****
        #endif
#define two
#endif

#ifdef _3
#define three
#endif

#ifdef two
#ifdef three
#error ****	define at most dimension two OR three		****
#endif
#endif

#ifndef two
#ifndef three
#error ****	define at least dimension two OR three		****
#endif
#endif

#ifdef two
#ifdef Sideon
#error ****   two dimensional case cannot have side data	****
#endif
#endif

/****************************************************************************/
/*																			*/
/* derive additional switches from commandline specified basic switches		*/
/*																			*/
/****************************************************************************/

#ifdef two
#define __TWODIM__
#endif

#ifdef three
#define __THREEDIM__
#endif

#ifdef ModelP
#define MODEL "PARALLEL"
#else
#define MODEL "SEQUENTIAL"
#endif

#undef __GRAPE_TRUE__
#ifdef _GRAPE
#define __GRAPE_TRUE__
#endif

#ifdef _GRAPE
#define GRAPE_SUPPORT "ON"
#else
#define GRAPE_SUPPORT "OFF"
#endif

#ifdef _NETGEN
#define NETGEN_SUPPORT "ON"
#else
#define NETGEN_SUPPORT "OFF"
#endif

#ifdef Debug
#define DEBUG_MODE "ON"
#else
#define DEBUG_MODE "OFF"
#endif

/****************************************************************************/
/*																			*/
/* "hard" switches for interpolation matrix and block-vectors				*/
/*																			*/
/****************************************************************************/

/* if interpolation matrix is stored */
#define __INTERPOLATION_MATRIX__

/* if block vector descriptors are used*/
#define __BLOCK_VECTOR_DESC__

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* necessary for most C runtime libraries */
#undef DOMAIN

/* some size parameters */
#define DIM_MAX                                 3       /* maximal space dimension				*/
#define DIM_OF_BND_MAX                  2       /* maximal dimension of boundary surface*/
#define MAXLEVEL                                32      /* maximum depth of triangulation		*/
#define MAXOBJECTS                              32      /* use 5 bits for object identification */
#define MAXSELECTION               100  /* max number of elements in selection	*/
#define TAGS                                     8      /* number of different element types	*/

/* some size macros for allocation purposes */
#define MAX_SIDES_OF_ELEM               6                       /* max number of sides of an elem	*/
#define MAX_EDGES_OF_ELEM               12                      /* max number of edges of an element*/
#define MAX_CORNERS_OF_ELEM             8                       /* max number of corners of an eleme*/
#define MAX_EDGES_OF_SIDE               4           /* max number of edges of a side	*/
#define MAX_EDGES_OF_CORNER             4                       /* max number of edges meeting in co*/
#define MAX_CORNERS_OF_SIDE     4           /* max number of corners of a side  */
#define MAX_CORNERS_OF_EDGE             2                       /* an edge has always two corners.. */
#define MAX_SIDES_OF_EDGE               2                       /* two sides have one edge in common*/
#define MAX_SONS                        30          /* max number of sons of an element */
#define MAX_SIDE_NODES                  9                       /* max number of nodes on elem side */

#define MAX_SIDES_TOUCHING              10                      /* max #fine sides touching a coarse*/

#define MAX_ELEM_VECTORS                (MAX_CORNERS_OF_ELEM+MAX_EDGES_OF_ELEM+1+MAX_SIDES_OF_ELEM)

#define MAX_NDOF_MOD_32         3           /* max number of doubles in a vector or matrix mod 32 */

/****************************************************************************/
/*																			*/
/* switch-define dependent defines											*/
/*																			*/
/****************************************************************************/

#ifdef __TWODIM__
        #define DIM                                             2                       /* space dimension								*/
        #define DIM_OF_BND                                      1                       /* dimension of boundary surface				*/
#endif

#ifdef __THREEDIM__
        #define DIM                                             3                       /* space dimension								*/
        #define DIM_OF_BND                                      2                       /* dimension of boundary surface				*/
#endif

/****************************************************************************/
/*																			*/
/* defines for algebra														*/
/*																			*/
/****************************************************************************/

#define MAXVOBJECTS                                             4                       /* four different data types					*/
#define MAXVECTORS                                              4                       /* max number of abstract vector types			*/
#if (MAXVECTORS<MAXVOBJECTS)
        #error *** MAXVECTORS must not be smaller than MAXVOBJECTS ***
#endif
#define NOVTYPE                                                 -1                      /* to indicate type not defined					*/
#define MAXDOMPARTS                                             4                       /* max number of geometric domain parts			*/

#define BITWISE_TYPE(t) (1<<(t))                                        /* transforms type into bitpattern				*/

/* derived sizes for algebra */
#define MAXMATRICES             ((MAXVECTORS*(MAXVECTORS+1))/2) /* max number of diff. matrix types			*/
#define MAXCONNECTIONS  (MAXMATRICES + MAXVECTORS)              /* max number of diff. connections              */

/* defines for vectors */
#define NODEVEC                                                 0                       /* vector associated to a node					*/
#define EDGEVEC                                                 1                       /* vector associated to an edge                                 */
#define ELEMVEC                                                 2                       /* vector associated to an element				*/
#define SIDEVEC                                                 3                       /* vector associated to an elementside			*/

/* some constants for abstract vector type names */
#define FROM_VTNAME                                             '0'
#define TO_VTNAME                                               'z'
#define MAXVTNAMES                                              (1+TO_VTNAME-FROM_VTNAME)

/* constants for blockvector description (BVD) */
#define NO_BLOCKVECTOR ((BLOCKNUMBER) ~0)        /* number for "there is no blockvector"; largest number of type BLOCKNUMBER */
#define MAX_BV_NUMBER (NO_BLOCKVECTOR - 1)      /* largest admissible blockvector number */
#define MAX_BV_LEVEL UCHAR_MAX                          /* largest admissible blockvector level number */
#define BVD_MAX_ENTRIES (sizeof(BVD_ENTRY_TYPE)*CHAR_BIT)       /* maximum number
                                                                   of entries in a BVD; NOTE: the actual available
                                                                   number of entries depends on the range of each entry */

/* constants for BLOCKVECTOR */
#define BVDOWNTYPEVECTOR        0       /* symbolic value for BVDOWNTYPE */
#define BVDOWNTYPEBV            1       /* symbolic value for BVDOWNTYPE */
#define BVDOWNTYPEDIAG          2       /* symbolic value for BVDOWNTYPE */
#define BV1DTV                          0       /* symbolic value for BVTVTYPE */
#define BV2DTV                          1       /* symbolic value for BVTVTYPE */
#define BVNOORIENTATION         0       /* no special orientation for BVORIENTATION */
#define BVHORIZONTAL            1       /* vectors form a horizontal line for BVORIENTATION */
#define BVVERTICAL                      2       /* vectors form a vertical line for BVORIENTATION */


/****************************************************************************/
/*																			*/
/* various defines															*/
/*																			*/
/****************************************************************************/

/* result codes of user supplied functions	0 = OK as usual */
#define OUT_OF_RANGE                    1       /* coordinate out of range				*/
#define CANNOT_INIT_PROBLEM     1       /* configProblem could not init problem */

/* use of GSTATUS (for grids), use power of 2 */
#define GSTATUS_BDF                             1
#define GSTATUS_INTERPOLATE             2
#define GSTATUS_ASSEMBLED               4

/* selection mode */
#define nodeSelection                   1                       /* objects selected are nodes			*/
#define elementSelection                2                       /* objects selected are elements		*/
#define vectorSelection                 3                       /* objects selected are vectors			*/

/* possible values for rule in MarkForRefinement */
#define NO_REFINEMENT           0
#define COPY                            1
#define RED                                     2
#define BLUE                            3
#define COARSE                          4
#ifdef __TWODIM__
#define BISECTION_1             5
#define BISECTION_2_Q           6
#define BISECTION_2_T1          7
#define BISECTION_2_T2          8
#define BISECTION_3             9
#endif
#ifdef __THREEDIM__
#define TETRA_RED_HEX           5
#endif

/* values for element class */
#define NO_CLASS                0
#define YELLOW_CLASS    1
#define GREEN_CLASS             2
#define RED_CLASS               3
#define SWITCH_CLASS    4

/* values for node types (relative to the father element of the vertex) */
#define CORNER_NODE             0
#define MID_NODE        1
#define SIDE_NODE       2
#define CENTER_NODE     3
#define LEVEL_0_NODE    4

/* macros for the control word management									*/
#define MAX_CONTROL_WORDS       20              /* maximum number of control words		*/
#define MAX_CONTROL_ENTRIES 100         /* max number of entries				*/

/* macros for the multigrid user data space management						*/
#define OFFSET_IN_MGUD(id)              (GetMGUDBlockDescriptor(id)->offset)
#define IS_MGUDBLOCK_DEF(id)    (GetMGUDBlockDescriptor(id)!=NULL)

/* REMARK: TOPNODE no more available since 970411
   because of problems in parallelisation
   to use it in serial version uncomment define
   #define TOPNODE(p)		(p)->iv.topnode
 */

/****************************************************************************/
/*																			*/
/* format definition data structure                                                                             */
/*																			*/
/****************************************************************************/

/*----------- general typedefs ---------------------------------------------*/

typedef DOUBLE DOUBLE_VECTOR[DIM];
typedef DOUBLE DOUBLE_VECTOR_2D[2];
typedef DOUBLE DOUBLE_VECTOR_3D[3];


/*----------- typedef for functions ----------------------------------------*/

typedef INT (*ConversionProcPtr)(               /* print user data --> string		*/
  void *,                                                               /* pointer to user data				*/
  const char *,                                                 /* prefix for each line				*/
  char *                                                                /* resulting string					*/
  );
typedef INT (*TaggedConversionProcPtr)( /* tagged print user data --> string*/
  INT,                                                                  /* tag for data identification		*/
  void *,                                                               /* pointer to user data				*/
  const char *,                                                 /* prefix for each line				*/
  char *                                                                /* resulting string					*/
  );


/*----------- definition of structs ----------------------------------------*/

struct format {

  /* fields for enironment variable */
  ENVDIR d;

  /* variables of format */
  INT sVertex;                                   /* size of vertex user data struc. in bytes*/
  INT sMultiGrid;                        /* size of mg user data structure in bytes	*/
  INT VectorSizes[MAXVECTORS];       /* number of doubles in vectors                    */
  char VTypeNames[MAXVECTORS];       /* a single char for abstract type name    */
  INT MatrixSizes[MAXMATRICES];      /* number of doubles in matrices			*/
#ifdef __INTERPOLATION_MATRIX__
  INT IMatrixSizes[MAXMATRICES];      /* number of doubles in matrices	        */
#endif
  INT ConnectionDepth[MAXMATRICES];      /* depth of connection for matrices    */
  INT elementdata;
  INT nodeelementlist;
  INT nodedata;

  /* print user data to string */
  ConversionProcPtr PrintVertex;
  ConversionProcPtr PrintGrid;
  ConversionProcPtr PrintMultigrid;
  TaggedConversionProcPtr PrintVector;          /* tag indicates VTYPE			*/
  TaggedConversionProcPtr PrintMatrix;          /* tag indicates MTP			*/

  /* table connecting parts, objects and types */
  INT po2t[MAXDOMPARTS][MAXVOBJECTS];                   /* (part,obj) --> vtype,		*/
  /* -1 if not defined			*/

  /* derived components */
  INT MaxConnectionDepth;               /* maximal connection depth                     */
  INT NeighborhoodDepth;                /* geometrical depth corresponding		        */
  /* algebraic con with depth 1			        */
  /* both derived from ConnectionDepth		*/
  INT t2p[MAXVECTORS];                  /* type --> part, bitwise, not unique		*/
  INT t2o[MAXVECTORS];                  /* type --> object, bitwise, not unique		*/
  /* both derived from po2t					*/
  char t2n[MAXVECTORS];                 /* type --> type name						*/
  INT n2t[MAXVTNAMES];                  /* type name --> type						*/
  INT OTypeUsed[MAXVOBJECTS];           /* 0 if vector not needed for geom object	*/
  INT MaxPart;                                  /* largest part used						*/
  /* both derived from po2t					*/
  INT MaxType;                                  /* largest type used						*/
  /* derived from VectorSizes					*/
} ;

typedef struct {
  int tp;                                               /* abstract type is described here	                */
  /* description only complete with po2t info	*/
  char name;                                            /* a single char as name of abstract type	*/
  int size;                                             /* data size in bytes						*/
} VectorDescriptor ;

typedef struct {
  int from;                                             /* This connection goes from position from	*/
  int to;                                               /* to position to							*/
  int size;                                             /* with size bytes per connection			*/
  int isize;                                            /* size of interpolation matrices			*/
  int depth;                                            /* connect with depth in dual graph             */
} MatrixDescriptor ;

/****************************************************************************/
/*																			*/
/* matrix/vector/blockvector data structure									*/
/*																			*/
/****************************************************************************/

/* data structure for BlockvectorDescription */
typedef unsigned INT BVD_ENTRY_TYPE;    /* memory providing storage for level numbers */
typedef unsigned INT BLOCKNUMBER;       /* valid numbers are 0..MAX_BV_NUMBER */
typedef unsigned char BLOCKLEVEL;       /* valid levels are 0..MAX_BV_LEVEL */

struct blockvector_description_format           /* describes how a struct of type
                                                                                           blockvector_description is to
                                                                                           be interpreted				*/
{
  INT bits;                                                                       /* bits per blocknumber entry */
  BLOCKLEVEL max_level;                                                   /* max. number of entries	*/
  BVD_ENTRY_TYPE level_mask[BVD_MAX_ENTRIES];
  /* level_mask[i] = mask entries for levels 0..i		*/
  BVD_ENTRY_TYPE neg_digit_mask[BVD_MAX_ENTRIES];
  /* neg_digit_mask[i] = masks out entry for level i	*/
};
typedef struct blockvector_description_format BV_DESC_FORMAT;

struct blockvector_description  /* describes the position of a blockvector	*/
{                                                               /* in a hierarchy of blockvectors			*/
  BVD_ENTRY_TYPE entry;                 /* sequence of block levels	according to	*/
  /* a certain blockvector_description_format	*/
  BLOCKLEVEL current;                           /* levels 0..current-1 currently valid		*/
  BLOCKLEVEL read;                              /* level read is next to be read			*/
};
typedef struct blockvector_description BV_DESC;

struct vector {

  unsigned INT control;                         /* object identification, various flags */
  union geom_object *object;                    /* associated object					*/

        #ifdef ModelP
  DDD_HEADER ddd;
        #endif

  struct vector *pred,*succ;                    /* double linked list of vectors		*/

  unsigned INT index;                           /* ordering of unknowns                                 */
  unsigned INT skip;                                    /* used bitwise to skip unknowns		*/
  struct matrix *start;                         /* implements matrix					*/

    #ifdef __BLOCK_VECTOR_DESC__
  BV_DESC block_descr;                          /* membership to the blockvector levels	*/
        #endif

    #ifdef __INTERPOLATION_MATRIX__
  struct matrix *istart;            /* implements interpolation matrix      */
    #endif

  /* user data */
  DOUBLE value[1];                                      /* array of doubles                                     */
};
typedef struct vector VECTOR;

struct matrix {
  unsigned INT control;                         /* object identification, various flags */

  struct matrix *next;                          /* row list                                                     */
  struct vector *vect;                          /* destination vector					*/

  /* user data */
  DOUBLE value[1];                                      /* array of doubles                                     */
};
typedef struct matrix MATRIX;
typedef struct matrix CONNECTION;

struct blockvector
{
  unsigned INT control;                         /* object identification, various flags	  */

  BLOCKNUMBER number;                                   /* logical blockvectornumber			  */
  struct blockvector *pred,*succ;       /* double linked list of blockvectors	  */
  VECTOR *first_vec;                                    /* start vector of this blockvector       */
  VECTOR *last_vec;                                     /* last vector of this blockvector		  */
  INT vec_number;                                       /* number of covered VECTORs			  */
  void *user_data;                                      /* pointer to any data					  */

  struct blockvector *first_son;        /* start of blockvector list on next level*/
  struct blockvector *last_son;         /* end of blockvector list on next level  */
};
typedef struct blockvector BLOCKVECTOR;

/****************************************************************************/
/*																			*/
/* unstructured grid data structures										*/
/*																			*/
/****************************************************************************/

/*----------- typedef for functions ----------------------------------------*/


/*----------- definition of structs ----------------------------------------*/

struct ivertex {                                        /* inner vertex structure				*/

  /* variables */
  unsigned INT control;                         /* object identification, various flags */
  INT id;                                                       /* unique id used for load/store		*/
  DOUBLE x[DIM];                                        /* vertex position						*/
  DOUBLE xi[DIM];                                       /* local coordinates in father element	*/

        #ifdef ModelP
  DDD_HEADER ddd;
        #endif

  /* pointers */
  union vertex *pred,*succ;                     /* double linked list of vertices		*/
  void *data;                                           /* associated user data structure		*/
  union element *father;                        /* father element						*/
        #ifdef TOPNODE
  /* REMARK: TOPNODE no more available since 970411
          because of problems in parallelisation
          to use it in serial version uncomment define TOPNODE */
  struct node *topnode;                         /* highest node where defect is valid	*/
        #endif
} ;

struct bvertex {                                        /* boundary vertex structure			*/

  /* variables */
  unsigned INT control;                         /* object identification, various flags */
  INT id;                                                       /* unique id used for load/store		*/
  DOUBLE x[DIM];                                        /* vertex position						*/
  DOUBLE xi[DIM];                                       /* local coordinates in father element	*/

        #ifdef ModelP
  DDD_HEADER ddd;
        #endif

  /* pointers */
  union vertex *pred,*succ;                     /* double linked list of vertices		*/
  void *data;                                           /* associated user data structure		*/
  union element *father;                        /* father element						*/
        #ifdef TOPNODE
  /* REMARK: TOPNODE no more available since 970411
          because of problems in parallelisation
          to use it in serial version uncomment define TOPNODE */
  struct node *topnode;                         /* highest node where defect is valid	*/
        #endif

  BNDP *bndp;                                       /* pointer boundary point decriptor		*/
} ;

union vertex {                                          /* only used to define pointer to vertex*/
  struct ivertex iv;
  struct bvertex bv;
} ;

struct elementlist {
  union element *el;
  struct elementlist *next;
} ;

struct node {                                           /* level dependent part of a vertex     */

  /* variables */
  unsigned INT control;                         /* object identification, various flags */
  INT id;                                                       /* unique id used for load/store		*/

        #ifdef ModelP
  DDD_HEADER ddd;
        #endif

  /* pointers */
  struct node *pred,*succ;                      /* double linked list of nodes per level*/
  struct link *start;                           /* list of links						*/
  union geom_object *father;                    /* node or edge on coarser level
                                                             (NULL if none) */
  struct node *son;                                     /* node on finer level (NULL if none)	*/
  union vertex *myvertex;                       /* corresponding vertex structure		*/

  /* WARNING: the allocation of the vector pointer depends on the format      */

  /* associated vector */
  VECTOR *vector;                                       /* associated vector					*/

  /* WARNING: the allocation of the data pointer depends on the format        */
  void *data;                                       /* associated data pointer              */
} ;

struct link {

  /* variables */
  unsigned INT control;                         /* object identification, various flags */
  struct link *next;                                    /* ptr to next link                                     */
  struct node *nbnode;                          /* ptr to neighbor node                                 */
} ;


struct edge {                                           /* undirected edge of the grid graph	*/

  /* variables */
  struct link links[2];                         /* two links							*/

        #if defined(ModelP) && defined(__THREEDIM__)
  DDD_HEADER ddd;
        #endif

  struct node *midnode;                         /* pointer to mid node on next finer gri*/

  /* WARNING: the allocation of the vector pointer depends on the format      */

  /* associated vector */
  VECTOR *vector;                                       /* associated vector					*/
} ;

struct generic_element {            /* no difference between inner and bndel*/

  /* variables */
  unsigned INT control;             /* object identification, various flags */
  INT id;                           /* unique id used for load/store        */
  unsigned INT flag;                /* additional flags for elements        */
  INT property;                         /* to store NodeOrder for hexahedrons   */

        #ifdef ModelP
  DDD_HEADER ddd;
  INT ptmp;                                             /* stores parition information			*/
  INT ptmp1;                                            /* stores cluster pointer				*/
  INT ptmp2;                                            /* stores number of descendents			*/
        #endif

  /* pointers */
  union element *pred, *succ;       /* double linked list of elements       */
  void *refs[1];                                        /* variable length array managed by ug  */
} ;

struct triangle {

  /* variables */
  unsigned INT control;                         /* object identification, various flags */
  INT id;                                                       /* unique id used for load/store		*/
  unsigned INT flag;                            /* additional flags for elements		*/
  INT property;                                 /* we need more bits ...				*/

        #ifdef ModelP
  DDD_HEADER ddd;
  INT ptmp;                                             /* stores parition information			*/
  INT ptmp1;                                            /* stores cluster pointer				*/
  INT ptmp2;                                            /* stores number of descendents			*/
        #endif

  /* pointers */
  union element *pred, *succ;           /* doubly linked list of elements		*/
  struct node *n[3];                                    /* corners of that element				*/
  union element *father;                        /* father element on coarser grid		*/
        #ifdef ModelP
  union element *sons[2];                       /* element tree                                                 */
        #else
  union element *sons[1];                       /* element tree                                                 */
        #endif
  union element *nb[3];                         /* dual graph							*/

  /* WARNING: the allocation of the vector pointer depends on the format      */
  /* void *ptr[4] would be possible too:                                      */
  /* if there are no element vectors, the sides will be ptr[0],ptr[1],ptr[2]  */
  /* Use the macros to find the correct address!                              */

  /* associated vector */
  VECTOR *vector;                                       /* associated vector					*/

  BNDS *bnds[3];                        /* only on bnd, NULL if interior side	*/

  /* WARNING: the allocation of the data pointer depends on the format        */
  void *data;                                       /* associated data pointer              */
} ;

struct quadrilateral {

  /* variables */
  unsigned INT control;                         /* object identification, various flags */
  INT id;                                                       /* unique id used for load/store		*/
  unsigned INT flag;                            /* additional flags for elements		*/
  INT property;                                 /* we need more bits ...				*/

        #ifdef ModelP
  DDD_HEADER ddd;
  INT ptmp;                                             /* stores parition information			*/
  INT ptmp1;                                            /* stores cluster pointer				*/
  INT ptmp2;                                            /* stores number of descendents			*/
        #endif

  /* pointers */
  union element *pred, *succ;           /* doubly linked list of elements		*/
  struct node *n[4];                                    /* corners of that element				*/
  union element *father;                        /* father element on coarser grid		*/
        #ifdef ModelP
  union element *sons[2];                       /* element tree                                                 */
        #else
  union element *sons[1];                       /* element tree                                                 */
        #endif
  union element *nb[4];                         /* dual graph							*/

  /* WARNING: the allocation of the vector pointer depends on the format      */
  /* void *ptr[5] would be possible too:                                      */
  /* if there are no element vectors, the sides will be ptr[0],ptr[1], ..     */
  /* Use the macros to find the correct address!                              */

  /* associated vector */
  VECTOR *vector;                                       /* associated vector					*/

  BNDS *bnds[4];                        /* only on bnd, NULL if interior side	*/

  /* WARNING: the allocation of the data pointer depends on the format        */
  void *data;                                       /* associated data pointer              */
} ;

struct tetrahedron {

  /* variables */
  unsigned INT control;                         /* object identification, various flags */
  INT id;                                                       /* unique id used for load/store		*/
  unsigned INT flag;                            /* additional flags for elements		*/
  INT property;                                 /* we need more bits ...				*/

        #ifdef ModelP
  DDD_HEADER ddd;
  INT ptmp;
        #endif

  /* pointers */
  union element *pred, *succ;           /* doubly linked list of elements		*/
  struct node *n[4];                                    /* corners of that element				*/
  union element *father;                        /* father element on coarser grid		*/
        #ifdef ModelP
  union element *sons[2];                       /* element tree                                                 */
        #else
  union element *sons[1];                       /* element tree                                                 */
        #endif
  union element *nb[4];                         /* dual graph							*/

  /* WARNING: the allocation of the vector pointer depends on the format      */
  /* void *ptr[9] would be possible too:                                      */
  /* if there are no element vectors, the sides will be ptr[0],ptr[1], ..     */
  /* Use the macros to find the correct address!                              */

  /* associated vector */
  VECTOR *vector;                                       /* associated vector					*/

  /* associated vector */
  VECTOR *sidevector[4];                        /* associated vectors for sides			*/

  BNDS *bnds[4];                        /* only on bnd, NULL if interior side	*/

  /* WARNING: the allocation of the data pointer depends on the format        */
  void *data;                                       /* associated data pointer              */
} ;

struct pyramid {

  /* variables */
  unsigned INT control;                         /* object identification, various flags */
  INT id;                                                       /* unique id used for load/store		*/
  unsigned INT flag;                            /* additional flags for elements		*/
  INT property;                                 /* we need more bits ...				*/

        #ifdef ModelP
  DDD_HEADER ddd;
  INT ptmp;
        #endif

  /* pointers */
  union element *pred, *succ;           /* doubly linked list of elements		*/
  struct node *n[5];                                    /* corners of that element				*/
  union element *father;                        /* father element on coarser grid		*/
        #ifdef ModelP
  union element *sons[2];                       /* element tree                                                 */
        #else
  union element *sons[1];                       /* element tree                                                 */
        #endif
  union element *nb[5];                         /* dual graph							*/

  /* WARNING: the allocation of the vector pointer depends on the format      */
  /* void *ptr[11] would be possible too:                                     */
  /* if there are no element vectors, the sides will be ptr[0],ptr[1], ..     */
  /* Use the macros to find the correct address!                              */

  /* associated vector */
  VECTOR *vector;                                       /* associated vector					*/

  /* associated vector */
  VECTOR *sidevector[5];                        /* associated vectors for sides			*/

  BNDS *bnds[5];                        /* only on bnd, NULL if interior side	*/

  /* WARNING: the allocation of the data pointer depends on the format        */
  void *data;                                       /* associated data pointer              */
} ;

struct prism {

  /* variables */
  unsigned INT control;                         /* object identification, various flags */
  INT id;                                                       /* unique id used for load/store		*/
  unsigned INT flag;                            /* additional flags for elements		*/
  INT property;                                 /* we need more bits ...				*/

        #ifdef ModelP
  DDD_HEADER ddd;
  INT ptmp;
        #endif

  /* pointers */
  union element *pred, *succ;           /* doubly linked list of elements		*/
  struct node *n[6];                                    /* corners of that element				*/
  union element *father;                        /* father element on coarser grid		*/
        #ifdef ModelP
  union element *sons[2];                       /* element tree                                                 */
        #else
  union element *sons[1];                       /* element tree                                                 */
        #endif
  union element *nb[5];                         /* dual graph							*/

  /* WARNING: the allocation of the vector pointer depends on the format      */
  /* void *ptr[11] would be possible too:                                     */
  /* if there are no element vectors, the sides will be ptr[0],ptr[1], ..     */
  /* Use the macros to find the correct address!                              */

  /* associated vector */
  VECTOR *vector;                                       /* associated vector					*/

  /* associated vector */
  VECTOR *sidevector[5];                        /* associated vectors for sides			*/

  BNDS *bnds[5];                        /* only on bnd, NULL if interior side	*/

  /* WARNING: the allocation of the data pointer depends on the format        */
  void *data;                                       /* associated data pointer              */
};

struct hexahedron {

  /* variables */
  unsigned INT control;                         /* object identification, various flags */
  INT id;                                                       /* unique id used for load/store		*/
  unsigned INT flag;                            /* additional flags for elements		*/
  INT property;                                 /* we need more bits ...				*/

        #ifdef ModelP
  DDD_HEADER ddd;
  INT ptmp;
        #endif

  /* pointers */
  union element *pred, *succ;           /* doubly linked list of elements		*/
  struct node *n[8];                                    /* corners of that element				*/
  union element *father;                        /* father element on coarser grid		*/
        #ifdef ModelP
  union element *sons[2];                       /* element tree                                                 */
        #else
  union element *sons[1];                       /* element tree                                                 */
        #endif
  union element *nb[6];                         /* dual graph							*/

  /* WARNING: the allocation of the vector pointer depends on the format      */
  /* void *ptr[13] would be possible too:                                     */
  /* if there are no element vectors, the sides will be ptr[0],ptr[1], ..     */
  /* Use the macros to find the correct address!                              */

  /* associated vector */
  VECTOR *vector;                                       /* associated vector					*/

  /* associated vector */
  VECTOR *sidevector[6];                        /* associated vectors for sides			*/

  BNDS *bnds[6];                        /* only on bnd, NULL if interior side	*/

  /* WARNING: the allocation of the data pointer depends on the format        */
  void *data;                                       /* associated data pointer              */
} ;

union element {
  struct generic_element ge;
    #ifdef __TWODIM__
  struct triangle tr;
  struct quadrilateral qu;
        #endif
    #ifdef __THREEDIM__
  struct tetrahedron te;
  struct pyramid py;
  struct prism pr;
  struct hexahedron he;
        #endif
} ;

union geom_object {                                             /* objects that can hold a vector		*/
  struct node nd;
  struct edge ed;
  union element el;
} ;

union selection_object {                                        /* objects than can be selected			*/
  struct node nd;
  union element el;
  struct vector ve;
} ;

typedef struct
{
  unsigned INT VecReserv[MAXVECTORS][MAX_NDOF_MOD_32];
  unsigned INT MatReserv[MAXMATRICES][MAX_NDOF_MOD_32];
  unsigned INT VecConsistentStatus[MAXMATRICES][MAX_NDOF_MOD_32];
  unsigned INT VecCollectStatus[MAXMATRICES][MAX_NDOF_MOD_32];
} DATA_STATUS;

struct grid {

  /* variables */
  unsigned INT control;                         /* object identification, various flags */
  INT status;                                           /* possible values see defines above	*/
  INT level;                                                    /* level of that grid					*/
  INT nVert;                                                    /* number of vertices					*/
  INT nNode;                                                    /* number of nodes						*/
  INT nElem;                                                    /* number of elements					*/
  INT nEdge;                                                    /* number of edges						*/
  INT nVector;                                          /* number of vectors					*/
  INT nCon;                                                     /* number of Connections				*/
#ifdef __INTERPOLATION_MATRIX__
  INT nIMat;                        /* number of interpolation matrices     */
#endif
  DATA_STATUS data_status;          /* memory management for vectors|matrix */
                                    /* status for consistent and collect    */
  /* pointers */
  union  element *elements[ELEMENT_LISTPARTS];       /* pointer to first element*/
  union  element *lastelement[ELEMENT_LISTPARTS];      /*pointer to last element*/
  union  vertex *vertices[VERTEX_LISTPARTS];            /* pointer to first vertex	*/
  union  vertex *lastvertex[VERTEX_LISTPARTS];      /* pointer to last vertex	*/
  struct node *firstNode[NODE_LISTPARTS];       /* pointer to first node		*/
  struct node *lastNode[NODE_LISTPARTS];        /* pointer to last node                 */
  VECTOR *firstVector[VECTOR_LISTPARTS];        /* pointer to first vector		*/
  VECTOR *lastVector[VECTOR_LISTPARTS];         /* pointer to last vector		*/
  BLOCKVECTOR *firstblockvector;        /* pointer to the first blockvector		*/
  BLOCKVECTOR *lastblockvector;         /* pointer to the last blockvector		*/
  struct grid *coarser, *finer;         /* coarser and finer grids				*/
  struct multigrid *mg;                         /* corresponding multigrid structure	*/
} ;

struct multigrid {

  /* env item */
  ENVDIR v;

  /* variables */
  INT status;                                           /* possible values, see above			*/
  INT magic_cookie;                                     /* used for identification			    */
  INT vertIdCounter;                                    /* count objects in that multigrid		*/
  INT nodeIdCounter;                                    /* count objects in that multigrid		*/
  INT elemIdCounter;                                    /* count objects in that multigrid		*/
  INT topLevel;                                         /* depth of the element tree			*/
  INT currentLevel;                                     /* level we are working on				*/
  INT fullrefineLevel;                          /* last level with complete surface     */
  INT bottomLevel;                                      /* bottom level for AMG                 */
  BVP *theBVP;                                          /* pointer to BndValProblem				*/
  BVP_DESC theBVPD;                                     /* description of BVP-properties		*/
  struct format *theFormat;                     /* pointer to format definition                 */
  HEAP *theHeap;                                        /* associated heap structure			*/
  INT nProperty;                                        /* max nb of properties used in elements*/

  DATA_STATUS data_status;          /* memory management for vectors|matrix */
                                    /* status for consistent and collect    */
  /* pointers */
  struct grid *amggrids[MAXLEVEL];      /* pointers to the grids				*/
  struct grid *grids[MAXLEVEL];         /* pointers to the grids				*/

  /* NodeElementPointerArray used for an O(n) InsertElement               */
  union element ***ndelemptrarray;                      /* pointer to the node element blocks   */

  /* selection */
  INT NbOfSelections;                           /* number of selected objects			*/
  INT SelectionMode;                                    /* selectionmode (see above)			*/
  union selection_object *Selection[MAXSELECTION];       /* pointer to selec obj*/

  /* user data */
  void *GenData;                                        /* general user data space				*/
  HEAP *UserHeap;                                       /* user heap							*/
  void *genpurp;                                        /* general purpose pointer				*/

  /* i/o handing */
  INT saved;                                                    /* 1 if multigrid saved					*/
  char filename[NAMESIZE];                      /* filename if saved					*/

  INT CoarseGridFixed;                          /* coarse grid complete					*/
} ;

/****************************************************************************/
/*																			*/
/*					typedef for structs                                                                     */
/*																			*/
/****************************************************************************/

/* geometrical part */
typedef struct format FORMAT;

typedef union  vertex VERTEX;
typedef struct elementlist ELEMENTLIST;
typedef struct node NODE;
typedef union  element ELEMENT;
typedef struct link LINK;
typedef struct edge EDGE;
typedef union  geom_object GEOM_OBJECT;
typedef union  selection_object SELECTION_OBJECT;
typedef struct grid GRID;
typedef struct multigrid MULTIGRID;

/****************************************************************************/
/*																			*/
/*					structs for evaluation functions						*/
/*																			*/
/****************************************************************************/

/*----------- typedef for functions ----------------------------------------*/

typedef INT (*PreprocessingProcPtr)(const char *, MULTIGRID *);
typedef DOUBLE (*ElementEvalProcPtr)(const ELEMENT *,const DOUBLE **,DOUBLE *);
typedef void (*ElementVectorProcPtr)(const ELEMENT *,const DOUBLE **,DOUBLE *,DOUBLE *);
typedef DOUBLE (*MatrixEvalProcPtr)(const MATRIX *);

/*----------- definition of structs ----------------------------------------*/

struct elementvalues {

  /* fields for enironment list variable */
  ENVVAR v;

  PreprocessingProcPtr PreprocessProc;                  /* prepare eval values					*/
  ElementEvalProcPtr EvalProc;                                  /* pointer to corresponding function	*/
} ;

struct elementvector {

  /* fields for enironment list variable */
  ENVVAR v;

  PreprocessingProcPtr PreprocessProc;                  /* prepare eval values					*/
  ElementVectorProcPtr EvalProc;                                /* pointer to corresponding function	*/
  int dimension;                                                                /* dimension of result vector			*/
} ;

struct matrixvalues {

  /* fields for enironment list variable */
  ENVVAR v;

  PreprocessingProcPtr PreprocessProc;                  /* prepare eval values					*/
  MatrixEvalProcPtr EvalProc;                                   /* pointer to corresponding function	*/
} ;

typedef struct elementvalues EVALUES ;
typedef struct elementvector EVECTOR ;
typedef struct matrixvalues MVALUES ;

/****************************************************************************/
/*																			*/
/* algebraic dependency for vector ordering                                                             */
/*																			*/
/****************************************************************************/

typedef INT (*DependencyProcPtr)(GRID *, const char *);

struct AlgebraicDependency {

  /* fields for enironment list variable */
  ENVVAR v;

  DependencyProcPtr DependencyProc;             /* pointer to dependency function			*/
} ;

typedef struct AlgebraicDependency ALG_DEP;

/****************************************************************************/
/*																			*/
/* find cut for vector ordering				                                                                */
/*																			*/
/****************************************************************************/

typedef VECTOR *(*FindCutProcPtr)(GRID *, VECTOR *, INT *);

typedef struct {

  /* fields for enironment list variable */
  ENVVAR v;

  FindCutProcPtr FindCutProc;           /* pointer to find cut function				*/

} FIND_CUT;

/****************************************************************************/
/*																			*/
/*					dynamic management of control words                                     */
/*																			*/
/****************************************************************************/

/* status of control word */
#define CW_FREE                                         0
#define CW_USED                                         1

/* status of control entry */
#define CE_FREE                                         0
#define CE_USED                                         1
#define CE_LOCKED                                       1

/* initializer macros for control entry and word predefines */
#define CW_INIT(used,cw,objs)                           {used, STR(cw), cw ## CW, cw ## OFFSET,objs}
#define CW_INIT_UNUSED                                          {CW_FREE,0,0,0}
#define CE_INIT(mode,cw,ce,objs)                        {mode, STR(ce), cw ## CW, ce ## CE, ce ## SHIFT, ce ## LEN, objs}
#define CE_INIT_UNUSED                                          {CE_FREE, 0, 0, 0, 0, 0, 0}

/* description of a control word */
typedef struct {
  INT used;                                                             /* this struct is used				*/
  char *name;                                                           /* name string						*/
  unsigned INT offset_in_object ;               /* where in object is it ?			*/
  INT objt_used;                                                /* bitwise object ID				*/
  unsigned INT used_mask ;                              /* used bits						*/
} CONTROL_WORD;

/* manage part of a control word */
typedef struct {
  INT used;                                                             /* this struct is used				*/
  char *name;                                                           /* name string						*/
  INT control_word ;                                            /* pointer to corresponding controlw*/
  INT offset_in_word;                                   /* shift in control word			*/
  INT length;                                                   /* number of bits used				*/
  INT objt_used;                                                /* bitwise object ID				*/
  unsigned INT offset_in_object;                /* copy from control word (faster)	*/
  unsigned INT mask;                                            /* 1 where bits are used			*/
  unsigned INT xor_mask;                                /* 0 where bits are used			*/
} CONTROL_ENTRY;

#ifndef __COMPILE_CW__
extern CONTROL_WORD
  control_words[MAX_CONTROL_WORDS];             /* global array with descriptions	*/
extern CONTROL_ENTRY
  control_entries[MAX_CONTROL_ENTRIES];      /* predefined control words		*/
#endif

/* general query macros */

/* dynamic control words */
#undef _DEBUG_CW_
#if (defined _DEBUG_CW_) && \
  !(defined __COMPILE_CW__)                             /* to avoid infinite recursion during ReadCW */

/* map cw read/write to functions */
        #define CW_READ(p,ce)                                   ReadCW(p,ce)
        #define CW_READ_STATIC(p,s,t)                   ReadCW(p,s ## CE)
        #define CW_WRITE(p,ce,n)                                WriteCW(p,ce,n)
        #define CW_WRITE_STATIC(p,s,t,n)                WriteCW(p,s ## CE,n)

#else   /* _DEBUG_CW_ */

        #define ControlWord(p,ce)  (((unsigned INT *)(p))[control_entries[ce].offset_in_object])

        #ifndef __T3E__
        #define CW_READ(p,ce)      ((ControlWord(p,ce) & control_entries[ce].mask)>>control_entries[ce].offset_in_word)
        #endif

/* very special hack */
        #ifdef __T3E__
        #define CW_READ(p,ce)      ((int)((ControlWord(p,ce) & control_entries[ce].mask)>>control_entries[ce].offset_in_word) )
        #endif

        #define CW_WRITE(p,ce,n)   ControlWord(p,ce) = (ControlWord(p,ce)&control_entries[ce].xor_mask)|(((n)<<control_entries[ce].offset_in_word)&control_entries[ce].mask)

/* static control words */
        #define StaticControlWord(p,t)            (((unsigned INT *)(p))[t ## OFFSET])
        #define StaticControlWordMask(s)          ((POW2(s ## LEN) - 1) << s ## SHIFT)

        #ifndef __T3E__
        #define CW_READ_STATIC(p,s,t)                                                \
  ((StaticControlWord(p,t) &  StaticControlWordMask(s)) >> s ## SHIFT)
        #endif

/* very special hack */
        #ifdef __T3E__
        #define CW_READ_STATIC(p,s,t)                                                \
  ((int)     ((StaticControlWord(p,t) &  StaticControlWordMask(s)) >> s ## SHIFT))
        #endif

        #define CW_WRITE_STATIC(p,s,t,n)                                             \
  StaticControlWord(p,t) =                                           \
    (StaticControlWord(p,t) &  (~StaticControlWordMask(s))) |          \
    (((n) << s ## SHIFT) & StaticControlWordMask(s))

#endif  /* _DEBUG_CW_ */

/* enumeration list of all control words of gm.h */
enum GM_CW {
  VECTOR_CW,
  MATRIX_CW,
  BLOCKVECTOR_CW,
  VERTEX_CW,
  NODE_CW,
  LINK_CW,
  EDGE_CW,
  ELEMENT_CW,
  FLAG_CW,
  PROPERTY_CW,
  GRID_CW,
  GRID_STATUS_CW,
  MULTIGRID_STATUS_CW,

  GM_N_CW
};

/* enumeration list of all control entry of gm.h */
enum GM_CE {
  VTYPE_CE,
  VOTYPE_CE,
  VPART_CE,
  VCOUNT_CE,
  VECTORSIDE_CE,
  VCLASS_CE,
  VDATATYPE_CE,
  VNCLASS_CE,
  VNEW_CE,
  VCCUT_CE,
  VCCOARSE_CE,
  NEW_DEFECT_CE,
  FINE_GRID_DOF_CE,
  MOFFSET_CE,
  MROOTTYPE_CE,
  MDESTTYPE_CE,
  MDIAG_CE,
  MSIZE_CE,
  MNEW_CE,
  CEXTRA_CE,
  MDOWN_CE,
  MUP_CE,
  BVDOWNTYPE_CE,
  BVLEVEL_CE,
  BVTVTYPE_CE,
  BVORIENTATION_CE,
  OBJ_CE,
  USED_CE,
  TAG_CE,
  LEVEL_CE,
  THEFLAG_CE,
  MOVE_CE,
  MOVED_CE,
  ONEDGE_CE,
  ONSIDE_CE,
  ONNBSIDE_CE,
  NOOFNODE_CE,
  NSUBDOM_CE,
  MODIFIED_CE,
  NTYPE_CE,
  NPROP_CE,
  LOFFSET_CE,
  NO_OF_ELEM_CE,
  AUXEDGE_CE,
  EDGENEW_CE,
  EDSUBDOM_CE,
  ECLASS_CE,
  NSONS_CE,
  NEWEL_CE,
  SUBDOMAIN_CE,
  NODEORD_CE,
  PROP_CE,
#ifdef ModelP
  XFERVECTOR_CE,
  XFERMATX_CE,
#endif

  GM_N_CE
};

/****************************************************************************/
/*																			*/
/* Macro definitions for algebra structures                                                             */
/*																			*/
/*																			*/
/* Use of the control word:                                                                                             */
/*																			*/
/* macro name|bits	|V|M|use												*/
/*																			*/
/* all objects:                                                                                                                         */
/*																			*/
/* vectors:                                                                                                                             */
/* VOTYPE	 |0 - 1 |*| | node-,edge-,side- or elemvector					*/
/* VCFLAG	 |3		|*| | flag for general use								*/
/* VCUSED	 |4		|*| | flag for general use								*/
/* VCOUNT	 |5-6	|*| |                                                                                                   */
/* VECTORSIDE|7 - 9 |*| | nb of side the side vect corr. to (in object elem)*/
/* VCLASS	 |11-12 |*| | class of v. (3: if corr. to red/green elem)		*/
/*					  (2: if corr. to first algebraic nb.)					*/
/*					  (1: if corr. to second algebraic nb.)                                 */
/* VDATATYPE |13-16 |*| | data type used bitwise							*/
/* VNCLASS	 |17-18 |*| | type of elem on finer grid the v. lies geom. in:	*/
/*							0: no elem on finer grid						*/
/*							1: elem of 'second alg. nbhood' only			*/
/*							2: elem of 'first alg. nbhood' only                     */
/*							3: red or green elem							*/
/* VNEW          |19	|*| | 1 if vector is new								*/
/* VCNEW	 |20	|*| | 1 if vector has a new connection					*/
/* VCCUT	 |26	|*| |                                                                                                   */
/* VTYPE	 |27-28 |*| | abstract vector type								*/
/* VPART	 |29-30 |*| | domain part										*/
/* VCCOARSE  |31    |*| | indicate algebraic part of VECTOR-MATRIX graph	*/
/*																			*/
/* matrices:																*/
/* MOFFSET	 |0     | |*| 0 if first matrix in connection else 1			*/
/* MROOTTYPE |1 - 2 | |*| VTYPE of root vector								*/
/* MDESTTYPE |3 - 4 | |*| VTYPE of destination vector						*/
/* MDIAG	 |5     | |*| 1 if diagonal matrix element						*/
/* MTYPE	 |6 -11 | |*| one of the (different sized) 10+4 matricees		*/
/* MUSED	 |12	| |*| general purpose flag								*/
/* MSIZE	 |13-27 | |*| size of the matrix in bytes						*/
/* MNEW          |28	| |*| 1 if matrix/connection is new                                     */
/* CEXTRA	 |29	| |*| 1 if is extra connection							*/
/*																			*/
/* Use of the control word in 'BLOCKVECTOR':								*/
/* BVDOWNTYPE 0	 BVDOWNTYPEVECTOR if the down component points to a vector,	*/
/*				 BVDOWNTYPEBV if it points to a further blockvector (son)	*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* general macros															*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* macros for VECTORs														*/
/*																			*/
/****************************************************************************/

/* control word offset */
#define VECTOR_OFFSET                           0

/* predefined control word entries */
#define VOTYPE_SHIFT                            0
#define VOTYPE_LEN                                      2
#define VOTYPE(p)                                       CW_READ_STATIC(p,VOTYPE_,VECTOR_)
#define SETVOTYPE(p,n)                          CW_WRITE_STATIC(p,VOTYPE_,VECTOR_,n)
#if (MAXVOBJECTS > POW2(VOTYPE_LEN))
        #error  *** VOTYPE_LEN too small ***
#endif

#define VTYPE_SHIFT                             2
#define VTYPE_LEN                                       2
#define VTYPE(p)                                        CW_READ_STATIC(p,VTYPE_,VECTOR_)
#define SETVTYPE(p,n)                           CW_WRITE_STATIC(p,VTYPE_,VECTOR_,n)
#if (MAXVTYPES > POW2(VTYPE_LEN))
        #error  *** VTYPE_LEN too small ***
#endif

#define VDATATYPE_SHIFT                         4
#define VDATATYPE_LEN                           4
#define VDATATYPE(p)                            CW_READ_STATIC(p,VDATATYPE_,VECTOR_)
#define SETVDATATYPE(p,n)                       CW_WRITE_STATIC(p,VDATATYPE_,VECTOR_,n)
#if (MAXVTYPES > VDATATYPE_LEN)
        #error  *** VDATATYPE_LEN too small ***
#endif

#define VCLASS_SHIFT                            8
#define VCLASS_LEN                                      2
#define VCLASS(p)                                       CW_READ_STATIC(p,VCLASS_,VECTOR_)
#define SETVCLASS(p,n)                          CW_WRITE_STATIC(p,VCLASS_,VECTOR_,n)

#define VNCLASS_SHIFT                           10
#define VNCLASS_LEN                             2
#define VNCLASS(p)                                      CW_READ_STATIC(p,VNCLASS_,VECTOR_)
#define SETVNCLASS(p,n)                         CW_WRITE_STATIC(p,VNCLASS_,VECTOR_,n)

#define VNEW_SHIFT                                      12
#define VNEW_LEN                                        1
#define VNEW(p)                                         CW_READ_STATIC(p,VNEW_,VECTOR_)
#define SETVNEW(p,n)                            CW_WRITE_STATIC(p,VNEW_,VECTOR_,n)

#define VCCUT_SHIFT                             13
#define VCCUT_LEN                                       1
#define VCCUT(p)                                        CW_READ_STATIC(p,VCCUT_,VECTOR_)
#define SETVCCUT(p,n)                           CW_WRITE_STATIC(p,VCCUT_,VECTOR_,n)

#define VCOUNT_SHIFT                            14
#define VCOUNT_LEN                                      2
#define VCOUNT(p)                                       CW_READ_STATIC(p,VCOUNT_,VECTOR_)
#define SETVCOUNT(p,n)                          CW_WRITE_STATIC(p,VCOUNT_,VECTOR_,n)

#define VECTORSIDE_SHIFT                        16
#define VECTORSIDE_LEN                          3
#define VECTORSIDE(p)                           CW_READ_STATIC(p,VECTORSIDE_,VECTOR_)
#define SETVECTORSIDE(p,n)                      CW_WRITE_STATIC(p,VECTORSIDE_,VECTOR_,n)

#define VCCOARSE_SHIFT                          19
#define VCCOARSE_LEN                            1
#define VCCOARSE(p)                                     CW_READ_STATIC(p,VCCOARSE_,VECTOR_)
#define SETVCCOARSE(p,n)                        CW_WRITE_STATIC(p,VCCOARSE_,VECTOR_,n)

#define FINE_GRID_DOF_SHIFT             20
#define FINE_GRID_DOF_LEN                       1
#define FINE_GRID_DOF(p)                        CW_READ_STATIC(p,FINE_GRID_DOF_,VECTOR_)
#define SETFINE_GRID_DOF(p,n)           CW_WRITE_STATIC(p,FINE_GRID_DOF_,VECTOR_,n)

#define NEW_DEFECT_SHIFT                        21
#define NEW_DEFECT_LEN                          1
#define NEW_DEFECT(p)                           CW_READ_STATIC(p,NEW_DEFECT_,VECTOR_)
#define SETNEW_DEFECT(p,n)                      CW_WRITE_STATIC(p,NEW_DEFECT_,VECTOR_,n)

#ifdef ModelP
#define XFERVECTOR_SHIFT                        20
#define XFERVECTOR_LEN                          2
#define XFERVECTOR(p)                           CW_READ(p,XFERVECTOR_CE)
#define SETXFERVECTOR(p,n)                      CW_WRITE(p,XFERVECTOR_CE,n)
#endif /* ModelP */

#define VPART_SHIFT                             22
#define VPART_LEN                                       2
#define VPART(p)                                        CW_READ_STATIC(p,VPART_,VECTOR_)
#define SETVPART(p,n)                           CW_WRITE_STATIC(p,VPART_,VECTOR_,n)
#if (MAXDOMPARTS > POW2(VPART_LEN))
        #error  *** VPART_LEN too small ***
#endif

#define VCFLAG(p)                                       THEFLAG(p)
#define SETVCFLAG(p,n)                          SETTHEFLAG(p,n)

#define VCUSED(p)                                       USED(p)
#define SETVCUSED(p,n)                          SETUSED(p,n)

#define VOBJECT(v)                                      ((v)->object)
#ifdef ModelP
#define PPREDVC(p,v)                            (((v)==PRIO_FIRSTVECTOR(p,PrioMaster)) ? \
                                                 PRIO_LASTVECTOR(p,PrioBorder) : (v)->pred)
#else
#define PPREDVC(p,v)                            ((v)->pred)
#endif
#define PREDVC(v)                                       ((v)->pred)
#define SUCCVC(v)                                       ((v)->succ)
#define VINDEX(v)                                       ((v)->index)
#define V_IN_DATATYPE(v,dt)                     (VDATATYPE(v) & (dt))
#define VSKIPME(v,n)                            ((((v)->skip)>>n) & 1)
#define VVECSKIP(v,n)                           ((((v)->skip)>>n) & 15)
#define VFULLSKIP(v,n)                          (VVECSKIP(v,n)==15)
#define SETVSKIPME(v,n)                         (((v)->skip=n))
#define VECSKIP(v)                                      ((v)->skip)
#define VECSKIPBIT(v,n)                         (((v)->skip) & (1<<n))
#define SETVECSKIPBIT(v,n)                      (v)->skip = ((v)->skip & (~(1<<n))) | (1<<n)
#define VSTART(v)                                       ((v)->start)
#ifdef __INTERPOLATION_MATRIX__
#define VISTART(v)                                      ((v)->istart)
#endif
#define VVALUE(v,n)                             ((v)->value[n])
#define VVALUEPTR(v,n)                          (&((v)->value[n]))
#define VMYNODE(v)                                      ((NODE*)((v)->object))
#define VMYEDGE(v)                                      ((EDGE*)((v)->object))
#define VMYELEMENT(v)                           ((ELEMENT*)((v)->object))
#define VUP(p)                                          LoWrd(VINDEX(p))
#define SETVUP(p,n)                             SetLoWrd(VINDEX(p),n)
#define VDOWN(p)                                        HiWrd(VINDEX(p))
#define SETVDOWN(p,n)                           SetHiWrd(VINDEX(p),n)
#ifdef __BLOCK_VECTOR_DESC__
#define VBVD(v)                                         ((v)->block_descr)
#define VMATCH(v,bvd,bvdf)                      BVD_IS_SUB_BLOCK( &(v)->block_descr, bvd, bvdf )
#endif

/* user for nodes, edges and elements */
#define CAST_NVECTOR(p)                         NVECTOR(p)
#define CAST_EDVECTOR(p)                        EDVECTOR(p)
#define CAST_SVECTOR(p,i)                       SVECTOR(p,i)
#define CAST_EVECTOR(p)                         EVECTOR(p)


/****************************************************************************/
/*																			*/
/* macros for MATRIXs														*/
/*																			*/
/****************************************************************************/

/* control word offset */
#define MATRIX_OFFSET                           0

#define MOFFSET_SHIFT                           0
#define MOFFSET_LEN                             1
#define MOFFSET(p)                                      CW_READ_STATIC(p,MOFFSET_,MATRIX_)
#define SETMOFFSET(p,n)                         CW_WRITE_STATIC(p,MOFFSET_,MATRIX_,n)

#define MROOTTYPE_SHIFT                         1
#define MROOTTYPE_LEN                           2
#define MROOTTYPE(p)                            CW_READ_STATIC(p,MROOTTYPE_,MATRIX_)
#define SETMROOTTYPE(p,n)                       CW_WRITE_STATIC(p,MROOTTYPE_,MATRIX_,n)
#if (MAXVTYPES > POW2(MROOTTYPE_LEN))
        #error  *** MROOTTYPE_LEN too small ***
#endif

#define MDESTTYPE_SHIFT                         3
#define MDESTTYPE_LEN                           2
#define MDESTTYPE(p)                            CW_READ_STATIC(p,MDESTTYPE_,MATRIX_)
#define SETMDESTTYPE(p,n)                       CW_WRITE_STATIC(p,MDESTTYPE_,MATRIX_,n)
#if (MAXVTYPES > POW2(MDESTTYPE_LEN))
        #error  *** MDESTTYPE_LEN too small ***
#endif

#define MDIAG_SHIFT                             5
#define MDIAG_LEN                                       1
#define MDIAG(p)                                        CW_READ_STATIC(p,MDIAG_,MATRIX_)
#define SETMDIAG(p,n)                           CW_WRITE_STATIC(p,MDIAG_,MATRIX_,n)

#define MNEW_SHIFT                                      6
#define MNEW_LEN                                        1
#define MNEW(p)                                         CW_READ_STATIC(p,MNEW_,MATRIX_)
#define SETMNEW(p,n)                            CW_WRITE_STATIC(p,MNEW_,MATRIX_,n)

#define CEXTRA_SHIFT                            7
#define CEXTRA_LEN                                      1
#define CEXTRA(p)                                       CW_READ_STATIC(p,CEXTRA_,MATRIX_)
#define SETCEXTRA(p,n)                          CW_WRITE_STATIC(p,CEXTRA_,MATRIX_,n)

#define MDOWN_SHIFT                             8
#define MDOWN_LEN                                       1
#define MDOWN(p)                                        CW_READ_STATIC(p,MDOWN_,MATRIX_)
#define SETMDOWN(p,n)                           CW_WRITE_STATIC(p,MDOWN_,MATRIX_,n)

#define MUP_SHIFT                                       9
#define MUP_LEN                                         1
#define MUP(p)                                          CW_READ_STATIC(p,MUP_,MATRIX_)
#define SETMUP(p,n)                             CW_WRITE_STATIC(p,MUP_,MATRIX_,n)

#define MSIZE_SHIFT                             10
#define MSIZE_LEN                                       15
#define MSIZEMAX                                        (POW2(MSIZE_LEN)-1)
#define MSIZE(p)                                        (CW_READ(p,MSIZE_CE)+sizeof(MATRIX)-sizeof(DOUBLE))
#define SETMSIZE(p,n)                           CW_WRITE(p,MSIZE_CE,(n-sizeof(MATRIX)+sizeof(DOUBLE)))

#define MTYPE(p)                                        (MatrixType[MROOTTYPE(p)][MDESTTYPE(p)])

#define MUSED(p)                                        USED(p)
#define SETMUSED(p,n)               SETUSED(p,n)

#ifdef ModelP
#define XFERMATX_SHIFT                          25
#define XFERMATX_LEN                            2
#define XFERMATX(p)                             CW_READ(p,XFERMATX_CE)
#define SETXFERMATX(p,n)                        CW_WRITE(p,XFERMATX_CE,n)
#endif

#define MINC(m)                                         ((MATRIX*)(((char *)(m))+MSIZE(m)))
#define MDEC(m)                                         ((MATRIX*)(((char *)(m))-MSIZE(m)))
#define MNEXT(m)                                        ((m)->next)
#define MDEST(m)                                        ((m)->vect)
#define MADJ(m)                                         ((MDIAG(m)) ? (m) : ((MOFFSET(m)) ? (MDEC(m)) : (MINC(m))))
#define MROOT(m)                                        MDEST(MADJ(m))
#define MMYCON(m)                                       ((MOFFSET(m)) ? (MDEC(m)) : (m))
#define MVALUE(m,n)                             ((m)->value[n])
#define MVALUEPTR(m,n)                          (&((m)->value[n]))
#define MDESTINDEX(m)                           ((m)->vect->index)
#define MSTRONG(p)                                      (MDOWN(p) && MUP(p))

/****************************************************************************/
/*																			*/
/* macros for CONNECTIONs													*/
/*																			*/
/****************************************************************************/

#define CMATRIX0(m)                             (m)
#define CMATRIX1(m)                             ((MDIAG(m)) ? (NULL) : (MINC(m)))
#define SETCUSED(c,n)                           {SETMUSED(CMATRIX0(c),n); SETMUSED(MADJ(CMATRIX0(c)),n);}

/****************************************************************************/
/*																			*/
/* macros for struct blockvector_description (BV_DESC)						*/
/*																			*/
/****************************************************************************/

/* access to members of struct blockvector_description (BV_DESC) */
#define BVD_NR_ENTRIES(bvd)                                     ((bvd)->current)

/* macros for blockvectordescription */

#define BVD_INIT(bvd)                                           ((bvd)->current=0)

/* sequential access operations on struct blockvector_description  (BV_DESC) */
#define BVD_PUSH_ENTRY(bvd,bnr,bvdf)            PushEntry( (bvd), (bnr), (bvdf) )
#define BVD_DISCARD_LAST_ENTRY(bvd)                     {assert(BVD_NR_ENTRIES(bvd)>0);BVD_NR_ENTRIES(bvd)--;}
#define BVD_INC_LAST_ENTRY(bvd,incr,bvdf)       ((bvd)->entry = (((((bvd)->entry >> ((bvdf)->bits * (BVD_NR_ENTRIES(bvd)-1))) + (incr)) & ((bvdf)->level_mask[0])) << ((bvdf)->bits * (BVD_NR_ENTRIES(bvd)-1))) | ((bvd)->entry & (bvdf)->neg_digit_mask[BVD_NR_ENTRIES(bvd)-1]))
#define BVD_DEC_LAST_ENTRY(bvd,decr,bvdf)       ((bvd)->entry = (((((bvd)->entry >> ((bvdf)->bits * (BVD_NR_ENTRIES(bvd)-1))) - (decr)) & ((bvdf)->level_mask[0])) << ((bvdf)->bits * (BVD_NR_ENTRIES(bvd)-1))) | ((bvd)->entry & (bvdf)->neg_digit_mask[BVD_NR_ENTRIES(bvd)-1]))
#define BVD_INIT_SEQ_READ(bvd)                          ((bvd)->read = 0)
#define BVD_READ_NEXT_ENTRY(bvd,bvdf)           ( ((bvd)->read<BVD_NR_ENTRIES(bvd)) ? BVD_GET_ENTRY((bvd),(bvd)->read++,(bvdf)) : NO_BLOCKVECTOR )

/* random access operations on struct blockvector_description  (BV_DESC) */
#define BVD_SET_ENTRY(bvd,level,bnr,bvdf)       ( (bvd)->entry = ( ((bvd)->entry & (bvdf)->neg_digit_mask[(level)]) | ( (bnr) << ( (bvdf)->bits*(level) ) ) ) )
#define BVD_GET_ENTRY(bvd,level,bvdf)           ( ((bvd)->entry >> ((bvdf)->bits * (level))) & (bvdf)->level_mask[0] )

#define BVD_IS_SUB_BLOCK(bvd_a,bvd_b,bvdf)      ( (BVD_NR_ENTRIES(bvd_a) >= BVD_NR_ENTRIES(bvd_b)) && (((bvd_a)->entry & (((bvdf)->level_mask[BVD_NR_ENTRIES(bvd_b)-1]))) == ((((bvd_b)->entry & (bvdf)->level_mask[BVD_NR_ENTRIES(bvd_b)-1])))))

/****************************************************************************/
/*																			*/
/* macros for BLOCKVECTOR													*/
/*																			*/
/****************************************************************************/

/* control word offset */
#define BLOCKVECTOR_OFFSET                              0

#define BVDOWNTYPE_SHIFT                                0
#define BVDOWNTYPE_LEN                                  2
#define BVDOWNTYPE(bv)                                  CW_READ_STATIC(bv,BVDOWNTYPE_,BLOCKVECTOR_)
#define SETBVDOWNTYPE(bv,n)                     CW_WRITE_STATIC(bv,BVDOWNTYPE_,BLOCKVECTOR_,n)

#define BVLEVEL_SHIFT                                   2
#define BVLEVEL_LEN                                     4
#define BVLEVEL(bv)                                             CW_READ_STATIC(bv,BVLEVEL_,BLOCKVECTOR_)
#define SETBVLEVEL(bv,n)                                CW_WRITE_STATIC(bv,BVLEVEL_,BLOCKVECTOR_,n)

#define BVTVTYPE_SHIFT                                  6
#define BVTVTYPE_LEN                                    1
#define BVTVTYPE(bv)                                    CW_READ_STATIC(bv,BVTVTYPE_,BLOCKVECTOR_)
#define SETBVTVTYPE(bv,n)                               CW_WRITE_STATIC(bv,BVTVTYPE_,BLOCKVECTOR_,n)

#define BVORIENTATION_SHIFT                             7
#define BVORIENTATION_LEN                               2
#define BVORIENTATION(bv)                               CW_READ_STATIC(bv,BVORIENTATION_,BLOCKVECTOR_)
#define SETBVORIENTATION(bv,n)                  CW_WRITE_STATIC(bv,BVORIENTATION_,BLOCKVECTOR_,n)

/* access to members of struct blockvector */
#define BVNUMBER(bv)                                    ((bv)->number)
#define BVUSERDATA(bv)                                  ((bv)->user_data)
#define BVPRED(bv)                                              ((bv)->pred)
#define BVSUCC(bv)                                              ((bv)->succ)
#define BVFIRSTVECTOR(bv)                               ((bv)->first_vec)
#define BVLASTVECTOR(bv)                                ((bv)->last_vec)
#define BVENDVECTOR(bv)                                 (BVSUCC(BVLASTVECTOR(bv)))
#define BVNUMBEROFVECTORS(bv)                   ((bv)->vec_number)
#define BVDOWNVECTOR(bv)                                ((bv)->first_vec)
#define BVDOWNBV(bv)                                    ((bv)->first_son)
#define BVDOWNBVLAST(bv)                                ((bv)->last_son)
#define BVDOWNBVEND(bv)                                 (BVSUCC(BVDOWNBVLAST(bv)))

#define BV_GEN_F                                                0
#define BV_GEN_L                                                1
#define BV_GEN_C                                                2
#define SETBV_GC(bv,flc,cyc)                    (BVNUMBER(bv)=3*(cyc)+(flc))
#define BV_IS_GC(bv,gen,cyc)                    (BVNUMBER(bv)%3==(gen) && BVNUMBER(bv)/3==(cyc))
#define BV_GEN(bv)                                              (BVNUMBER(bv)%3)                        /* gen.: FLC    */
#define BV_CYC(bv)                                              (BVNUMBER(bv)/3)                        /* nb. of cycle */
#define SET_ORD_GCL(n,flc,cyc,line)             {SetLoWrd(n,3*(cyc)+(flc)); SetHiWrd(n,line);}
#define ORD_GEN(n)                                              (LoWrd(n)%3)            /* gen.: FLC    */
#define ORD_CYC(n)                                              (LoWrd(n)/3)            /* nb. of cycle */
#define ORD_LIN(n)                                              HiWrd(n)

/* operations on struct block */
#define BV_IS_EMPTY(bv)                                 (BVNUMBEROFVECTORS(bv)==0)
#define BV_IS_LEAF_BV(bv)                               (BVDOWNTYPE(bv)==BVDOWNTYPEVECTOR)
#define BV_IS_DIAG_BV(bv)                               (BVDOWNTYPE(bv)==BVDOWNTYPEDIAG)

/****************************************************************************/
/*																			*/
/* Macro definitions for geometric objects									*/
/*																			*/
/*																			*/
/* Use of the control word:                                                                                             */
/*																			*/
/* macro name|bits	|V|N|L|E|V|M|  use										*/
/*							 C											*/
/* all objects:                                                                                                                         */
/* TAG		 |18-20 | | | |*| | |general purpose tag field					*/
/* LEVEL	 |21-25 |*|*| |*| | |level of a node/element (imp. for copies)	*/
/* THEFLAG	 |26	|*|*|*|*| | |general purp.,  leave them as you found 'em*/
/* USED          |27	|*|*|*|*| | |object visited, leave them as you found 'em*/
/* OBJT          |28-31 |*|*|*|*| | |object type identification                                 */

/*																			*/
/* vertices:																*/
/* MOVED	 |0     |*| | | | | |boundary vertex not lying on edge midpoint */
/* MOVE          |1-2	|*| | | | | |vertex can be moved on a 0(1,2,3) dim subsp*/
/* ONEDGE	 |3 - 6 |*| | | | | |no. of edge in father element				*/
/* ONSIDE	 |3 - 5 |*| | | | | |no. of side in father element				*/
/* ONNBSIDE	 |6 - 8 |*| | | | | |no. of side in the neigbor of the father   */
/* NOOFNODE	 |9 -13 |*| | | | | |???									    */
/*																			*/
/* nodes:																	*/
/* NSUBDOM	 |0-3	| |*| | | | |subdomain id                                       */
/* MODIFIED  |6         | |*| | | | |1 if node must be assembled				*/
/* N_OUTFLOW |0-7	|														*/
/* N_INFLOW  |8-15	|														*/
/*																			*/
/* links and edges:                                                                                                             */
/* LOFFSET	 |0     | | |*| | | |position of link in links array			*/
/* EDGENEW	 |1     | | |*| | | |status of edge								*/
/* NOOFELEM  |2-8	| | |*| | | |nb. of elem. the edge is part of			*/
/* AUXEDGE	 |9		|														*/
/* EDSUBDOM  |12-17 | | | |*| | |subdomain of edge if inner edge, 0 else	*/
/*																			*/
/* elements:																*/
/* ECLASS	 |8-9	| | | |*| | |element class from enumeration type		*/
/* NSONS	 |10-13 | | | |*| | |number of sons                                                     */
/* NEWEL	 |14	| | | |*| | |element just created						*/
/* VSIDES	 |11-14 | | | |*| | |viewable sides                                                     */
/* NORDER	 |15-19 | | | |*| | |view position order of the nodes			*/
/* CUTMODE	 |26-27 | | | |*| | |elem intersects cutplane or...                     */
/*																			*/
/****************************************************************************/

/* object identification */
enum GM_OBJECTS {

  MGOBJ,                                                /* multigrid object	                        */
  IVOBJ,                                                /* inner vertex                                         */
  BVOBJ,                                                /* boundary vertex					*/
  IEOBJ,                                                /* inner element					*/
  BEOBJ,                                                /* boundary element                             */
  EDOBJ,                                                /* edge object						*/
  NDOBJ,                                                /* node object						*/
  GROBJ,                                                /* grid object						*/

  /* object numbers for algebra */
  VEOBJ,                                                /* vector object					*/
  MAOBJ,                                                /* matrix object					*/
  BLOCKVOBJ,                        /* blockvector object               */

  NPREDEFOBJ,                                           /* no of predefined objects             */

  NOOBJ = -1
};
#define LIOBJ           EDOBJ           /* link and edge are identified		        */
#define COOBJ           MAOBJ           /* connection and matrix are identified		*/

/****************************************************************************/
/*																			*/
/* general macros															*/
/*																			*/
/****************************************************************************/

/* control word offset */
#define GENERAL_CW                                      NODE_CW         /* any of the geom objects	*/
#define GENERAL_OFFSET                          0

#define OBJ_SHIFT                                       28
#define OBJ_LEN                                         4
#define OBJT(p)                                         CW_READ_STATIC(p,OBJ_,GENERAL_)
#define SETOBJT(p,n)                            CW_WRITE_STATIC(p,OBJ_,GENERAL_,n)
#define OBJT_MAX                                        (POW2(OBJ_LEN)-1)

#define USED_SHIFT                                      27
#define USED_LEN                                        1
#define USED(p)                                         CW_READ_STATIC(p,USED_,GENERAL_)
#define SETUSED(p,n)                            CW_WRITE_STATIC(p,USED_,GENERAL_,n)

#define THEFLAG_SHIFT                           26
#define THEFLAG_LEN                             1
#define THEFLAG(p)                                      CW_READ_STATIC(p,THEFLAG_,GENERAL_)
#define SETTHEFLAG(p,n)                         CW_WRITE_STATIC(p,THEFLAG_,GENERAL_,n)

#define LEVEL_SHIFT                             21
#define LEVEL_LEN                                       5
#define LEVEL(p)                                        CW_READ_STATIC(p,LEVEL_,GENERAL_)
#define SETLEVEL(p,n)                           CW_WRITE_STATIC(p,LEVEL_,GENERAL_,n)

#define TAG_SHIFT                                       18
#define TAG_LEN                                         3
#define TAG(p)                                          CW_READ_STATIC(p,TAG_,GENERAL_)
#define SETTAG(p,n)                             CW_WRITE_STATIC(p,TAG_,GENERAL_,n)

#define CTRL(p)         (*((unsigned INT *)(p)))
#define ID(p)           (((INT *)(p))[1])

/****************************************************************************/
/*																			*/
/* macros for vertices														*/
/*																			*/
/****************************************************************************/

/* control word offset */
#define VERTEX_OFFSET                           0

#define MOVE_SHIFT                                      1
#define MOVE_LEN                                        2
#define MOVE(p)                                         CW_READ_STATIC(p,MOVE_,VERTEX_)
#define SETMOVE(p,n)                            CW_WRITE_STATIC(p,MOVE_,VERTEX_,n)

#define MOVED_SHIFT                             0
#define MOVED_LEN                                       1
#define MOVED(p)                                        CW_READ_STATIC(p,MOVED_,VERTEX_)
#define SETMOVED(p,n)                           CW_WRITE_STATIC(p,MOVED_,VERTEX_,n)

#define ONEDGE_SHIFT                            3
#define ONEDGE_LEN                                      4
#define ONEDGE(p)                                       CW_READ_STATIC(p,ONEDGE_,VERTEX_)
#define SETONEDGE(p,n)                          CW_WRITE_STATIC(p,ONEDGE_,VERTEX_,n)

/* the following two overlap with ONEDGE */
#define ONSIDE_SHIFT                            3
#define ONSIDE_LEN                                      3
#define ONSIDE(p)                                       CW_READ_STATIC(p,ONSIDE_,VERTEX_)
#define SETONSIDE(p,n)                          CW_WRITE_STATIC(p,ONSIDE_,VERTEX_,n)

#define ONNBSIDE_SHIFT                          6
#define ONNBSIDE_LEN                            3
#define ONNBSIDE(p)                                     CW_READ_STATIC(p,ONNBSIDE_,VERTEX_)
#define SETONNBSIDE(p,n)                        CW_WRITE_STATIC(p,ONNBSIDE_,VERTEX_,n)

#define NOOFNODE_SHIFT                          9
#define NOOFNODE_LEN                            5
#define NOOFNODEMAX                                     POW2(NOOFNODE_LEN)
#if (MAXLEVEL > NOOFNODEMAX)
#error  ****  set NOOFNODEMAX/_LEN appropriate to MAXLEVEL: 2^NOOFNODE_LEN = NOOFNODEMAX >= MAXLEVEL ****
#endif
#define NOOFNODE(p)                                     CW_READ_STATIC(p,NOOFNODE_,VERTEX_)
#define SETNOOFNODE(p,n)                        CW_WRITE_STATIC(p,NOOFNODE_,VERTEX_,n)
#define INCNOOFNODE(p)                          SETNOOFNODE(p,NOOFNODE(p)+1)
#define DECNOOFNODE(p)                          SETNOOFNODE(p,NOOFNODE(p)-1)

#define PREDV(p)                (p)->iv.pred
#define SUCCV(p)                (p)->iv.succ
#define CVECT(p)                (p)->iv.x
#define XC(p)                   (p)->iv.x[0]
#define YC(p)                   (p)->iv.x[1]
#define ZC(p)                   (p)->iv.x[2]
#define LCVECT(p)               (p)->iv.xi
#define XI(p)                   (p)->iv.xi[0]
#define ETA(p)                  (p)->iv.xi[1]
#define NU(p)                   (p)->iv.xi[2]
#define VDATA(p)                (p)->iv.data
#define VFATHER(p)              (p)->iv.father

/* for boundary vertices */
#define V_BNDP(p)               (p)->bv.bndp

/* parallel macros */
#ifdef ModelP
#define PARHDRV(p)                      (&((p)->iv.ddd))
#endif /* ModelP */

/****************************************************************************/
/*																			*/
/* macros for nodes                                                                                                             */
/*																			*/
/****************************************************************************/

/* control word offset */
#define NODE_OFFSET                                     0

#define NTYPE_SHIFT                                     0
#define NTYPE_LEN                                       3
#define NTYPE(p)                                        CW_READ_STATIC(p,NTYPE_,NODE_)
#define SETNTYPE(p,n)                           CW_WRITE_STATIC(p,NTYPE_,NODE_,n)

#define NSUBDOM_SHIFT                           3
#define NSUBDOM_LEN                                     6
#define NSUBDOM(p)                                      CW_READ_STATIC(p,NSUBDOM_,NODE_)
#define SETNSUBDOM(p,n)                         CW_WRITE_STATIC(p,NSUBDOM_,NODE_,n)

#define NPROP_SHIFT                 11
#define NPROP_LEN                   4
#define NPROP(p)                    CW_READ_STATIC(p,NPROP_,NODE_)
#define SETNPROP(p,n)               CW_WRITE_STATIC(p,NPROP_,NODE_,n)

#define MODIFIED_SHIFT                          15
#define MODIFIED_LEN                            1
#define MODIFIED(p)                             CW_READ_STATIC(p,MODIFIED_,NODE_)
#define SETMODIFIED(p,n)                        CW_WRITE_STATIC(p,MODIFIED_,NODE_,n)

#define PREDN(p)        (p)->pred
#define SUCCN(p)        (p)->succ
#define START(p)        (p)->start

#define NFATHER(p)                      (p)->father
#define SETNFATHER(p,n)         (p)->father = n
/*
   #define NFATHER(p)			((NTYPE(p) == CORNER_NODE) ? (p)->father : NULL)
   #define NFATHEREDGE(p)		((NTYPE(p) == MID_NODE) ? (EDGE *)(p)->father : NULL)
   #define SETNFATHEREDGE(p,e)	(p)->father = (NODE *) e
 */
#define CORNERTYPE(p)   (NTYPE(p) == CORNER_NODE)
#define MIDTYPE(p)              (NTYPE(p) == MID_NODE)
#define SIDETYPE(p)             (NTYPE(p) == SIDE_NODE)
#define CENTERTYPE(p)   (NTYPE(p) == CENTER_NODE)

#define SONNODE(p)      (p)->son
#define MYVERTEX(p) (p)->myvertex
#define NDATA(p)        (p)->data
#define NVECTOR(p)      (p)->vector

#define NODE_ELEMENT_LIST(p)    ((ELEMENTLIST *)(p)->data)
#define ELEMENT_PTR(p)                  ((p)->el)

/****************************************************************************/
/*																			*/
/* macros for links                                                                                                             */
/*																			*/
/****************************************************************************/

/* CAUTION: the controlword of LINK0 and its edge are identical (AVOID overlapping of flags) */

/* control word offset */
#define LINK_OFFSET                             0

#define LOFFSET_SHIFT                           0
#define LOFFSET_LEN                             1
#define LOFFSET(p)                                      CW_READ(p,LOFFSET_CE)
#define SETLOFFSET(p,n)                         CW_WRITE(p,LOFFSET_CE,n)

#define NBNODE(p)       (p)->nbnode
#define NEXT(p)         (p)->next
#define LDATA(p)        (p)->matelem
#define MATELEM(p)      (p)->matelem  /* can be used for node and link */

#define MYEDGE(p)       ((EDGE *)((p)-LOFFSET(p)))
#define REVERSE(p)      ((p)+(1-LOFFSET(p)*2))

/****************************************************************************/
/*																			*/
/* macros for edges                                                                                                             */
/*																			*/
/****************************************************************************/

/* control word offset */
#define EDGE_OFFSET                             0

#define NO_OF_ELEM_SHIFT                        2
#define NO_OF_ELEM_LEN                          7
#define NO_OF_ELEM_MAX                          128
#define NO_OF_ELEM(p)                           CW_READ(p,NO_OF_ELEM_CE)
#define SET_NO_OF_ELEM(p,n)             CW_WRITE(p,NO_OF_ELEM_CE,n)
#define INC_NO_OF_ELEM(p)                       SET_NO_OF_ELEM(p,NO_OF_ELEM(p)+1)
#define DEC_NO_OF_ELEM(p)                       SET_NO_OF_ELEM(p,NO_OF_ELEM(p)-1)

#define AUXEDGE_SHIFT                           9
#define AUXEDGE_LEN                             1
#define AUXEDGE(p)                                      CW_READ(p,AUXEDGE_CE)
#define SETAUXEDGE(p,n)                         CW_WRITE(p,AUXEDGE_CE,n)

#define EDGENEW_SHIFT                           1
#define EDGENEW_LEN                             1
#define EDGENEW(p)                                      CW_READ(p,EDGENEW_CE)
#define SETEDGENEW(p,n)                         CW_WRITE(p,EDGENEW_CE,n)

/* boundary edges will be indicated by a subdomain id of 0 */
#define EDSUBDOM_SHIFT                          12
#define EDSUBDOM_LEN                            6
#define EDSUBDOM(p)                                     CW_READ(p,EDSUBDOM_CE)
#define SETEDSUBDOM(p,n)                        CW_WRITE(p,EDSUBDOM_CE,n)

#define LINK0(p)        (&((p)->links[0]))
#define LINK1(p)        (&((p)->links[1]))
#define MIDNODE(p)      ((p)->midnode)
#define EDDATA(p)       ((p)->data)
#define EDVECTOR(p) ((p)->vector)

/****************************************************************************/
/*																			*/
/* macros for elements														*/
/*																			*/
/****************************************************************************/

/* TAG values */
#define TRIANGLE                3
#define QUADRILATERAL   4
#define TETRAHEDRON     4
#define PYRAMID                 5
#define PRISM           6
#define HEXAHEDRON              7

/* control word offsets */
#define ELEMENT_OFFSET                                  0
#define FLAG_OFFSET                                             2
#define PROPERTY_OFFSET                                 3

/* macros for control word */
#define ECLASS_SHIFT                                    8
#define ECLASS_LEN                                              2
#define ECLASS(p)                                               CW_READ(p,ECLASS_CE)
#define SETECLASS(p,n)                                  CW_WRITE(p,ECLASS_CE,n)

#define NSONS_SHIFT                                     10
#define NSONS_LEN                                               5
#define NSONS(p)                                                CW_READ(p,NSONS_CE)
#define SETNSONS(p,n)                                   CW_WRITE(p,NSONS_CE,n)

#define NEWEL_SHIFT                                     17
#define NEWEL_LEN                                               1
#define NEWEL(p)                                                CW_READ(p,NEWEL_CE)
#define SETNEWEL(p,n)                                   CW_WRITE(p,NEWEL_CE,n)

/* macros for flag word                           */
/* are obviously all for internal use */

/* the property field */
#define SUBDOMAIN_SHIFT                 24
#define SUBDOMAIN_LEN                   6
#define SUBDOMAIN(p)                    CW_READ(p,SUBDOMAIN_CE)
#define SETSUBDOMAIN(p,n)               CW_WRITE(p,SUBDOMAIN_CE,n)

#define NODEORD_SHIFT                   0
#define NODEORD_LEN                     24
#define NODEORD(p)                      CW_READ(p,NODEORD_CE)
#define SETNODEORD(p,n)                 CW_WRITE(p,NODEORD_CE,n)

#define PROP_SHIFT                      30
#define PROP_LEN                        2
#define PROP(p)                         CW_READ(p,PROP_CE)
#define SETPROP(p,n)                    CW_WRITE(p,PROP_CE,n)

/* parallel macros */
#ifdef ModelP
#define PARTITION(p)                    ((p)->ge.ptmp)
#define PARHDRE(p)                      (&((p)->ge.ddd))
#endif


/*******************************/
/* the general element concept */
/*******************************/

/* this structure contains all topological properties of an element and more .. */
typedef struct {
  INT tag;                                                                      /* element type to be defined       */

  /* the following parameters determine size of refs array in element */
  INT max_sons_of_elem;                                         /* max number of sons for this type */
  INT sides_of_elem;                                                    /* how many sides ?					*/
  INT corners_of_elem;                                          /* how many corners ?				*/

  /* local geometric description of the element */
  DOUBLE_VECTOR local_corner[MAX_CORNERS_OF_ELEM];                      /* local coordinates of the corners of the element */

  /* more size parameters */
  INT edges_of_elem;                                                    /* how many edges ?					*/
  INT edges_of_side[MAX_SIDES_OF_ELEM];         /* number of edges for each side	*/
  INT corners_of_side[MAX_SIDES_OF_ELEM];       /* number of corners for each side  */
  INT corners_of_edge;                                          /* is always 2 !					*/

  /* index computations */
  /* Within each element sides, edges, corners are numbered in some way.      */
  /* Within each side the edges and corners are numbered, within the edge the */
  /* corners are numbered. The following arrays map the local numbers within  */
  /* the side or edge to the numbering within the element.					*/
  INT edge_of_side[MAX_SIDES_OF_ELEM][MAX_EDGES_OF_SIDE];
  INT corner_of_side[MAX_SIDES_OF_ELEM][MAX_CORNERS_OF_SIDE];
  INT corner_of_edge[MAX_EDGES_OF_ELEM][MAX_CORNERS_OF_EDGE];

  /* the following parameters are derived from data above */
  INT mapped_inner_objt;                                        /* tag to objt mapping for free list*/
  INT mapped_bnd_objt;                                          /* tag to objt mapping for free list*/
  INT inner_size, bnd_size;                                     /* size in bytes used for alloc     */
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

  /* ... the refinement rules should be placed here later */
} GENERAL_ELEMENT;

#ifndef __COMPILE_EL__
/* these are the offsets into the variable length pointer array of the element */
extern INT n_offset[TAGS];
extern INT father_offset[TAGS];
extern INT sons_offset[TAGS];
extern INT nb_offset[TAGS];
extern INT evector_offset[TAGS];
extern INT svector_offset[TAGS];
extern INT side_offset[TAGS];
extern INT data_offset[TAGS];

/* the element descriptions are also globally available, these are pointers ! */
extern GENERAL_ELEMENT *element_descriptors[TAGS], *reference_descriptors[MAX_CORNERS_OF_ELEM+1];
#endif

/****************************************************************************/
/*																			*/
/* macros for element descriptors											*/
/*																			*/
/****************************************************************************/

/* macros to access element descriptors by element pointers		*/

#define SIDES_OF_ELEM(p)                (element_descriptors[TAG(p)]->sides_of_elem)
#define EDGES_OF_ELEM(p)                (element_descriptors[TAG(p)]->edges_of_elem)
#define CORNERS_OF_ELEM(p)              (element_descriptors[TAG(p)]->corners_of_elem)
#define LOCAL_COORD_OF_ELEM(p,c)        (element_descriptors[TAG(p)]->local_corner[(c)])


#define SONS_OF_ELEM(p)                 (element_descriptors[TAG(p)]->max_sons_of_elem) /* this is the number of pointers ! */

#define EDGES_OF_SIDE(p,i)              (element_descriptors[TAG(p)]->edges_of_side[(i)])
#define CORNERS_OF_SIDE(p,i)            (element_descriptors[TAG(p)]->corners_of_side[(i)])

#define CORNERS_OF_EDGE                         2

#define EDGE_OF_SIDE(p,s,e)             (element_descriptors[TAG(p)]->edge_of_side[(s)][(e)])
#define EDGE_OF_TWO_SIDES(p,s,t)        (element_descriptors[TAG(p)]->edge_of_two_sides[(s)][(t)])
#define CORNER_OF_SIDE(p,s,c)           (element_descriptors[TAG(p)]->corner_of_side[(s)][(c)])
#define CORNER_OF_EDGE(p,e,c)           (element_descriptors[TAG(p)]->corner_of_edge[(e)][(c)])

#define EDGE_WITH_CORNERS(p,c0,c1)      (element_descriptors[TAG(p)]->edge_with_corners[(c0)][(c1)])
#define SIDE_WITH_EDGE(p,e,k)           (element_descriptors[TAG(p)]->side_with_edge[(e)][(k)])
#define CORNER_OF_SIDE_INV(p,s,c)       (element_descriptors[TAG(p)]->corner_of_side_inv[(s)][(c)])
#define EDGES_OF_CORNER(p,c,k)          (element_descriptors[TAG(p)]->edges_of_corner[(c)][(k)])
#define CORNER_OF_OPPEDGE(p,e,c)        (element_descriptors[TAG(p)]->corner_of_oppedge[(e)][(c)])
#define CORNER_OPP_TO_SIDE(p,s)         (element_descriptors[TAG(p)]->corner_opp_to_side[(s)])
#define OPPOSITE_EDGE(p,e)                  (element_descriptors[TAG(p)]->opposite_edge[(e)])
#define SIDE_OPP_TO_CORNER(p,c)         (element_descriptors[TAG(p)]->side_opp_to_corner[(c)])
#define EDGE_OF_CORNER(p,c,e)           (element_descriptors[TAG(p)]->edge_of_corner[(c)][(e)])

#define CTRL2(p)        ((p)->ge.flag)
#define FLAG(p)                 (p)->ge.flag
#define SUCCE(p)                (p)->ge.succ
#define PREDE(p)                (p)->ge.pred

#define CORNER(p,i)     ((NODE *) (p)->ge.refs[n_offset[TAG(p)]+(i)])
#define EFATHER(p)              ((ELEMENT *) (p)->ge.refs[father_offset[TAG(p)]])
#define SON(p,i)                ((ELEMENT *) (p)->ge.refs[sons_offset[TAG(p)]+(i)])
#define NBELEM(p,i)     ((ELEMENT *) (p)->ge.refs[nb_offset[TAG(p)]+(i)])
#define ELEM_BNDS(p,i)  ((BNDS *) (p)->ge.refs[side_offset[TAG(p)]+(i)])
#define EVECTOR(p)              ((VECTOR *) (p)->ge.refs[evector_offset[TAG(p)]])
#define SVECTOR(p,i)    ((VECTOR *) (p)->ge.refs[svector_offset[TAG(p)]+(i)])
#define EDATA(p)            ((void *) (p)->ge.refs[data_offset[TAG(p)]])
#define SIDE_ON_BND(p,i) (ELEM_BNDS(p,i) != NULL)
#define INNER_SIDE(p,i)  (ELEM_BNDS(p,i) == NULL)
#define INNER_BOUNDARY(p,i) (InnerBoundary(p,i))
/* TODO: replace by function call */

#ifdef __TWODIM__
#define EDGE_ON_BND(p,i) (ELEM_BNDS(p,i) != NULL)
#endif

#ifdef __THREEDIM__
#define EDGE_ON_BND(p,i) (SIDE_ON_BND(p,SIDE_WITH_EDGE(p,i,0)) || \
                          SIDE_ON_BND(p,SIDE_WITH_EDGE(p,i,1)))
#endif

/* use the following macros to assign values, since definition  */
/* above is no proper lvalue.									*/
#define SET_CORNER(p,i,q)       (p)->ge.refs[n_offset[TAG(p)]+(i)] = q
#define SET_EFATHER(p,q)        (p)->ge.refs[father_offset[TAG(p)]] = q
#define SET_SON(p,i,q)          (p)->ge.refs[sons_offset[TAG(p)]+(i)] = q
#define SET_NBELEM(p,i,q)       (p)->ge.refs[nb_offset[TAG(p)]+(i)] = q
#define VOID_NBELEM(p,i)        (p)->ge.refs[nb_offset[TAG(p)]+(i)]
#define SET_BNDS(p,i,q)         (p)->ge.refs[side_offset[TAG(p)]+(i)] = q
#define SET_EVECTOR(p,q)        (p)->ge.refs[evector_offset[TAG(p)]] = q
#define SET_SVECTOR(p,i,q)      (p)->ge.refs[svector_offset[TAG(p)]+(i)] = q
#define SET_EDATA(p,q)      (p)->ge.refs[data_offset[TAG(p)]] = q

#define SideBndCond(t,side,l,v,type)  BNDS_BndCond(ELEM_BNDS(t,side),l,NULL,v,type)
#define Vertex_BndCond(p,w,i,v,t)     BNDP_BndCond(V_BNDP(p),w,i,NULL,v,t)


/* macros to access corner pointers directly */

#define CORNER_OF_EDGE_PTR(e,i,j)               (CORNER(e,CORNER_OF_EDGE(e,i,j)))
#define CORNER_OF_SIDE_PTR(e,i,j)               (CORNER(e,CORNER_OF_SIDE(e,i,j)))


/* macros to access element descriptors by element tags	*/

#define INNER_SIZE_TAG(t)                       (element_descriptors[t]->inner_size)
#define BND_SIZE_TAG(t)                         (element_descriptors[t]->bnd_size)
#define MAPPED_INNER_OBJT_TAG(t)                (element_descriptors[t]->mapped_inner_objt)
#define MAPPED_BND_OBJT_TAG(t)                  (element_descriptors[t]->mapped_bnd_objt)

#define SIDES_OF_TAG(t)                         (element_descriptors[t]->sides_of_elem)
#define EDGES_OF_TAG(t)                         (element_descriptors[t]->edges_of_elem)
#define CORNERS_OF_TAG(t)                       (element_descriptors[t]->corners_of_elem)
#define LOCAL_COORD_OF_TAG(t,c)                 (element_descriptors[t]->local_corner[(c)])

#define SONS_OF_TAG(t)                                  (element_descriptors[t]->max_sons_of_elem) /* this is the number of pointers ! */

#define EDGES_OF_SIDE_TAG(t,i)                  (element_descriptors[t]->edges_of_side[(i)])
#define CORNERS_OF_SIDE_TAG(t,i)                (element_descriptors[t]->corners_of_side[(i)])

#define EDGE_OF_SIDE_TAG(t,s,e)                 (element_descriptors[t]->edge_of_side[(s)][(e)])
#define EDGE_OF_TWO_SIDES_TAG(t,s,u)    (element_descriptors[t]->edge_of_two_sides[(s)][(u)])
#define CORNER_OF_SIDE_TAG(t,s,c)               (element_descriptors[t]->corner_of_side[(s)][(c)])
#define CORNER_OF_EDGE_TAG(t,e,c)               (element_descriptors[t]->corner_of_edge[(e)][(c)])

#define EDGE_WITH_CORNERS_TAG(t,c0,c1)  (element_descriptors[t]->edge_with_corners[(c0)][(c1)])
#define SIDE_WITH_EDGE_TAG(t,e,k)               (element_descriptors[t]->side_with_edge[(e)][(k)])
#define CORNER_OF_SIDE_INV_TAG(t,s,c)   (element_descriptors[t]->corner_of_side_inv[(s)][(c)])
#define EDGES_OF_CORNER_TAG(t,c,k)              (element_descriptors[t]->edges_of_corner[(c)][(k)])
#define CORNER_OF_OPPEDGE_TAG(t,e,c)    (element_descriptors[t]->corner_of_oppedge[(e)][(c)])
#define CORNER_OPP_TO_SIDE_TAG(t,s)             (element_descriptors[t]->corner_opp_to_side[(s)])
#define OPPOSITE_EDGE_TAG(t,e)              (element_descriptors[t]->opposite_edge[(e)])
#define SIDE_OPP_TO_CORNER_TAG(t,c)             (element_descriptors[t]->side_opp_to_corner[(c)])
#define EDGE_OF_CORNER_TAG(t,c,e)               (element_descriptors[t]->edge_of_corner[(c)][(e)])


/* macros to access reference descriptors by number of element corners	*/

#define SIDES_OF_REF(n)                   (reference_descriptors[n]->sides_of_elem)
#define EDGES_OF_REF(n)                   (reference_descriptors[n]->edges_of_elem)
#define CORNERS_OF_REF(n)                 (reference_descriptors[n]->corners_of_elem)
#define LOCAL_COORD_OF_REF(n,c)           (reference_descriptors[n]->local_corner[(c)])
#define EDGES_OF_SIDE_REF(n,i)            (reference_descriptors[n]->edges_of_side[(i)])
#define CORNERS_OF_SIDE_REF(n,i)          (reference_descriptors[n]->corners_of_side[(i)])
#define EDGE_OF_SIDE_REF(n,s,e)           (reference_descriptors[n]->edge_of_side[(s)][(e)])
#define EDGE_OF_TWO_SIDES_REF(n,s,t)  (reference_descriptors[n]->edge_of_two_sides[(s)][(t)])
#define CORNER_OF_SIDE_REF(n,s,c)         (reference_descriptors[n]->corner_of_side[(s)][(c)])
#define CORNER_OF_EDGE_REF(n,e,c)         (reference_descriptors[n]->corner_of_edge[(e)][(c)])
#define EDGE_WITH_CORNERS_REF(n,c0,c1) (reference_descriptors[n]->edge_with_corners[(c0)][(c1)])
#define SIDE_WITH_EDGE_REF(n,e,k)         (reference_descriptors[n]->side_with_edge[(e)][(k)])
#define CORNER_OF_SIDE_INV_REF(n,s,c) (reference_descriptors[n]->corner_of_side_inv[(s)][(c)])
#define EDGES_OF_CORNER_REF(n,c,k)        (reference_descriptors[n]->edges_of_corner[(c)][(k)])
#define CORNER_OF_OPPEDGE_REF(n,e,c)  (reference_descriptors[n]->corner_of_oppedge[(e)][(c)])
#define CORNER_OPP_TO_SIDE_REF(n,s)       (reference_descriptors[n]->corner_opp_to_side[(s)])
#define OPPOSITE_EDGE_REF(n,e)            (reference_descriptors[n]->opposite_edge[(e)])
#define SIDE_OPP_TO_CORNER_REF(n,c)       (reference_descriptors[n]->side_opp_to_corner[(c)])
#define EDGE_OF_CORNER_REF(n,c,e)         (reference_descriptors[n]->edge_of_corner[(c)][(e)])

/****************************************************************************/
/*																			*/
/* macros for grids                                                                                                             */
/*																			*/
/****************************************************************************/

/* control word offset */
#define GRID_OFFSET                                             0
#define GRID_STATUS_OFFSET                              1

#define GLEVEL(p)                                       ((p)->level)
#define SETGLOBALGSTATUS(p)             ((p)->status=~0)
#define GSTATUS(p,n)                            ((p)->status&(n))
#define RESETGSTATUS(p,n)                       ((p)->status&=~(n))

#ifdef ModelP
#define PFIRSTELEMENT(p)                                ((LISTPART_FIRSTELEMENT(p,0)!=NULL) ?\
                                                         (LISTPART_FIRSTELEMENT(p,0)) : (FIRSTELEMENT(p)))
#define PRIO_FIRSTELEMENT(p,prio)               ((p)->elements[PRIO2LISTPART(ELEMENT_LIST,prio)])
#define LISTPART_FIRSTELEMENT(p,part)   ((p)->elements[part])
#define FIRSTELEMENT(p)                                 ((p)->elements[PRIO2LISTPART(ELEMENT_LIST,PrioMaster)])

#define PLASTELEMENT(p)                                 LASTELEMENT(p)
#define PRIO_LASTELEMENT(p,prio)                ((p)->lastelement[PRIO2LISTPART(ELEMENT_LIST,prio)])
#define LISTPART_LASTELEMENT(p,part)    ((p)->lastelement[part])
#define LASTELEMENT(p)                                  ((p)->lastelement[PRIO2LISTPART(ELEMENT_LIST,PrioMaster)])
#else
#define FIRSTELEMENT(p)         ((p)->elements[0])
#define PFIRSTELEMENT(p)        FIRSTELEMENT(p)
#define LASTELEMENT(p)          ((p)->lastelement[0])
#define PLASTELEMENT(p)         LASTELMENT(p)
#endif

#ifdef ModelP
#define PFIRSTVERTEX(p)                                 ((LISTPART_FIRSTVERTEX(p,0)!=NULL) ?\
                                                         (LISTPART_FIRSTVERTEX(p,0)) :\
                                                         ((LISTPART_FIRSTVERTEX(p,1)!=NULL) ?\
                                                          (LISTPART_FIRSTVERTEX(p,1)) : (FIRSTVERTEX(p))))
#define PRIO_FIRSTVERTEX(p,prio)                ((p)->vertices[PRIO2LISTPART(VERTEX_LIST,prio)])
#define LISTPART_FIRSTVERTEX(p,part)    ((p)->vertices[part])
#define FIRSTVERTEX(p)                                  (((p)->vertices[PRIO2LISTPART(VERTEX_LIST,\
                                                                                      PrioBorder)]!=NULL) ?\
                                                         (p)->vertices[PRIO2LISTPART(VERTEX_LIST,PrioBorder)] :\
                                                         (p)->vertices[PRIO2LISTPART(VERTEX_LIST,PrioMaster)])
#define SFIRSTVERTEX(p)                                  (p)->vertices[PRIO2LISTPART(VERTEX_LIST,PrioMaster)]

#define PLASTVERTEX(p)                                  LASTVERTEX(p)
#define PRIO_LASTVERTEX(p,prio)                 ((p)->lastvertex[PRIO2LISTPART(VERTEX_LIST,prio)])
#define LISTPART_LASTVERTEX(p,part)     ((p)->lastvertex[part])
#define LASTVERTEX(p)                                   ((p)->lastvertex[PRIO2LISTPART(VERTEX_LIST,PrioMaster)])
#else
#define FIRSTVERTEX(p)          ((p)->vertices[0])
#define PFIRSTVERTEX(p)         FIRSTVERTEX(p)
#define SFIRSTVERTEX(p)         FIRSTVERTEX(p)
#define LASTVERTEX(p)           ((p)->lastvertex[0])
#define PLASTVERTEX(p)          LASTVERTEX(p)
#endif

#define FIRSTELEMSIDE(p)        ((p)->sides)

#ifdef ModelP
#define PFIRSTNODE(p)                                   ((LISTPART_FIRSTNODE(p,0)!=NULL) ?\
                                                         (LISTPART_FIRSTNODE(p,0)) :\
                                                         ((LISTPART_FIRSTNODE(p,1)!=NULL) ?\
                                                          (LISTPART_FIRSTNODE(p,1)) : (FIRSTNODE(p))))
#define PRIO_FIRSTNODE(p,prio)                  ((p)->firstNode[PRIO2LISTPART(NODE_LIST,prio)])
#define LISTPART_FIRSTNODE(p,part)              ((p)->firstNode[part])
#define FIRSTNODE(p)                                    (((p)->firstNode[PRIO2LISTPART(NODE_LIST,\
                                                                                       PrioBorder)]!=NULL) ?\
                                                         (p)->firstNode[PRIO2LISTPART(NODE_LIST,PrioBorder)] :\
                                                         (p)->firstNode[PRIO2LISTPART(NODE_LIST,PrioMaster)])
#define SFIRSTNODE(p)                                    (p)->firstNode[PRIO2LISTPART(NODE_LIST,PrioMaster)]

#define PLASTNODE(p)                                    LASTNODE(p)
#define PRIO_LASTNODE(p,prio)                   ((p)->lastNode[PRIO2LISTPART(NODE_LIST,prio)])
#define LISTPART_LASTNODE(p,part)               ((p)->lastNode[part])
#define LASTNODE(p)                                     ((p)->lastNode[PRIO2LISTPART(NODE_LIST,PrioMaster)])
#else
#define FIRSTNODE(p)            ((p)->firstNode[0])
#define PFIRSTNODE(p)           FIRSTNODE(p)
#define SFIRSTNODE(p)           FIRSTNODE(p)
#define LASTNODE(p)             ((p)->lastNode[0])
#define PLASTNODE(p)            LASTNODE(p)
#endif

#ifdef ModelP
#define PFIRSTVECTOR(p)                                 ((LISTPART_FIRSTVECTOR(p,0)!=NULL) ?\
                                                         (LISTPART_FIRSTVECTOR(p,0)) :\
                                                         ((LISTPART_FIRSTVECTOR(p,1)!=NULL) ?\
                                                          (LISTPART_FIRSTVECTOR(p,1)) : (FIRSTVECTOR(p))))
#define PRIO_FIRSTVECTOR(p,prio)                ((p)->firstVector[PRIO2LISTPART(VECTOR_LIST,prio)])
#define LISTPART_FIRSTVECTOR(p,part)    ((p)->firstVector[part])
#define FIRSTVECTOR(p)                                  (((p)->firstVector[PRIO2LISTPART(VECTOR_LIST,\
                                                                                         PrioBorder)]!=NULL) ?\
                                                         (p)->firstVector[PRIO2LISTPART(VECTOR_LIST,PrioBorder)] :\
                                                         (p)->firstVector[PRIO2LISTPART(VECTOR_LIST,PrioMaster)])
#define SFIRSTVECTOR(p)                                  (p)->firstVector[PRIO2LISTPART(VECTOR_LIST,PrioMaster)]

#define PLASTVECTOR(p)                                  LASTVECTOR(p)
#define PRIO_LASTVECTOR(p,prio)                 ((p)->lastVector[PRIO2LISTPART(VECTOR_LIST,prio)])
#define LISTPART_LASTVECTOR(p,part)     ((p)->lastVector[part])
#define LASTVECTOR(p)                                   ((p)->lastVector[PRIO2LISTPART(VECTOR_LIST,PrioMaster)])
#else
#define FIRSTVECTOR(p)          ((p)->firstVector[0])
#define PFIRSTVECTOR(p)         FIRSTVECTOR(p)
#define SFIRSTVECTOR(p)         FIRSTVECTOR(p)
#define LASTVECTOR(p)           ((p)->lastVector[0])
#define PLASTVECTOR(p)          LASTVECTOR(p)
#endif

#define GFIRSTBV(p)             ((p)->firstblockvector)
#define GLASTBV(p)                      ((p)->lastblockvector)
#define UPGRID(p)                       ((p)->finer)
#define DOWNGRID(p)             ((p)->coarser)
#define MYMG(p)                         ((p)->mg)
#define NV(p)                           ((p)->nVert)
#define NN(p)                           ((p)->nNode)
#define NT(p)                           ((p)->nElem)
#define NE(p)                           ((p)->nEdge)
#define NS(p)                           ((p)->nSide)
#define NVEC(p)                         ((p)->nVector)
#define NC(p)                           ((p)->nCon)
#define VEC_DEF_IN_OBJ_OF_GRID(p,tp)     ((p)->mg->theFormat->OTypeUsed[(tp)]>0)
#define NIMAT(p)                        ((p)->nIMat)
#define NELIST_DEF_IN_GRID(p)  ((p)->mg->theFormat->nodeelementlist)
#define EDATA_DEF_IN_GRID(p)   ((p)->mg->theFormat->elementdata)
#define NDATA_DEF_IN_GRID(p)   ((p)->mg->theFormat->nodedata)

/****************************************************************************/
/*																			*/
/* macros for multigrids													*/
/*																			*/
/****************************************************************************/

/* control word offset */
#define MULTIGRID_STATUS_OFFSET           ((sizeof(ENVDIR))/sizeof(unsigned INT))

#define MGSTATUS(p)                     ((p)->status)
#define RESETMGSTATUS(p)                {(p)->status=0; (p)->magic_cookie = (int)time(NULL); (p)->saved=0;}
#define MG_MAGIC_COOKIE(p)              ((p)->magic_cookie)
#define VIDCNT(p)                       ((p)->vertIdCounter)
#define NIDCNT(p)                       ((p)->nodeIdCounter)
#define EIDCNT(p)                       ((p)->elemIdCounter)
#define TOPLEVEL(p)                     ((p)->topLevel)
#define BOTTOMLEVEL(p)                  ((p)->bottomLevel)
#define CURRENTLEVEL(p)                 ((p)->currentLevel)
#define FULLREFINELEVEL(p)              ((p)->fullrefineLevel)
#define MGFORMAT(p)                     ((p)->theFormat)
#define DATAFORMAT(p)                   ((p)->theFormat)
#define MG_BVP(p)                               ((p)->theBVP)
#define MG_BVPD(p)                              (&((p)->theBVPD))
#define MGBNDSEGDESC(p,i)               (&((p)->segments[i]))
#define MGVERTEX(p,k)                   ((p)->corners[k])
#define MGNOOFCORNERS(p)                ((p)->numOfCorners)
#define MGHEAP(p)                               ((p)->theHeap)
#define MG_NPROPERTY(p)                 ((p)->nProperty)
#define GRID_ON_LEVEL(p,i)              ((p)->grids[i])
/* macros for the NodeElementsBlockArray . . .  */
#define ELEMS_OF_NODE_MAX               75
#define NDELEM_BLKS_MAX                 100
#define NO_NODES_OF_BLK                 1000
#define MGNDELEMPTRARRAY(p)             ((p)->ndelemptrarray)
#define MGNDELEMBLK(p,i)                (*(((p)->ndelemptrarray)+i))
#define MGNDELEMOFFS(i,o)               (i*ELEMS_OF_NODE_MAX+o)
#define MGNDELEMBLKENTRY(p,b,i) (*((*(((p)->ndelemptrarray)+b))+i))
/* . . . macros for the NodeElementsBlockArray  */
#define SELECTIONSIZE(p)                ((p)->NbOfSelections)
#define SELECTIONMODE(p)                ((p)->SelectionMode)
#define SELECTIONOBJECT(p,i)    ((p)->Selection[(((i)<MAXSELECTION) ? (i) : (MAXSELECTION-1))])
#define MGNAME(p)                               ((p)->v.name)
#define MG_USER_HEAP(p)                 ((p)->UserHeap)
#define GEN_MGUD(p)                     ((p)->GenData)
#define GEN_MGUD_ADR(p,o)               ((void *)(((char *)((p)->GenData))+(o)))
#define VEC_DEF_IN_OBJ_OF_MG(p,tp)       ((p)->theFormat->OTypeUsed[(tp)]>0)
#define NELIST_DEF_IN_MG(p)     ((p)->theFormat->nodeelementlist)
#define EDATA_DEF_IN_MG(p)      ((p)->theFormat->elementdata)
#define NDATA_DEF_IN_MG(p)      ((p)->theFormat->nodedata)
#define MG_GENPURP(p)                   ((p)->genpurp)
#define MG_SAVED(p)                             ((p)->saved)
#define MG_FILENAME(p)                  ((p)->filename)
#define MG_COARSE_FIXED(p)              ((p)->CoarseGridFixed)

/****************************************************************************/
/*																			*/
/* macros for formats														*/
/*																			*/
/****************************************************************************/

#define FMT_ELEM_DATA(f)                                ((f)->elementdata)
#define FMT_NODE_DATA(f)                                ((f)->nodedata)
#define FMT_NODE_ELEM_LIST(f)                   ((f)->nodeelementlist)
#define FMT_S_VERTEX(f)                                 ((f)->sVertex)
#define FMT_S_MG(f)                                             ((f)->sMultiGrid)
#define FMT_S_VEC_TP(f,t)                               ((f)->VectorSizes[t])
#define FMT_VTYPE_NAME(f,t)                             ((f)->VTypeNames[t])
#define FMT_S_MAT_TP(f,t)                               ((f)->MatrixSizes[t])
#define FMT_S_MATPTR(f)                                 ((f)->MatrixSizes)
#define FMT_S_IMAT_TP(f,t)                              ((f)->IMatrixSizes[t])
#define FMT_CONN_DEPTH_TP(f,t)                  ((f)->ConnectionDepth[t])
#define FMT_CONN_DEPTH_PTR(f)                   ((f)->ConnectionDepth)
#define FMT_CONN_DEPTH_MAX(f)                   ((f)->MaxConnectionDepth)
#define FMT_NB_DEPTH(f)                                 ((f)->NeighborhoodDepth)
#define FMT_PR_VERTEX(f)                                ((f)->PrintVertex)
#define FMT_PR_GRID(f)                                  ((f)->PrintGrid)
#define FMT_PR_MG(f)                                    ((f)->PrintMultigrid)
#define FMT_PR_VEC(f)                                   ((f)->PrintVector)
#define FMT_PR_MAT(f)                                   ((f)->PrintMatrix)
#define FMT_PO2T(f,p,o)                                 ((f)->po2t[p][o])
#define FMT_T2P(f,t)                                    ((f)->t2p[t])
#define FMT_TYPE_IN_PART(f,t,o)                 ((f)->t2p[o] & (1<<o))
#define FMT_T2O(f,o)                                    ((f)->t2o[o])
#define FMT_TYPE_USES_OBJ(f,t,o)                ((f)->t2o[o] & (1<<o))
#define FMT_USES_OBJ(f,o)                               ((f)->OTypeUsed[o])
#define FMT_MAX_PART(f)                                 ((f)->MaxPart)
#define FMT_MAX_TYPE(f)                                 ((f)->MaxType)

#define FMT_N2T(f,c)                                    (((c)<FROM_VTNAME) ? NOVTYPE : ((c)>TO_VTNAME) ? NOVTYPE : (f)->n2t[(c)-FROM_VTNAME])
#define FMT_SET_N2T(f,c,t)                              ((f)->n2t[(c)-FROM_VTNAME] = t)
#define FMT_T2N(f,t)                                    (((f)->t2n[t]))

/* constants for USED flags of objects */
#define MG_ELEMUSED     1
#define MG_NODEUSED     2
#define MG_EDGEUSED     4
#define MG_VERTEXUSED   8
#define MG_VECTORUSED   16
#define MG_MATRIXUSED   32

/****************************************************************************/
/*																			*/
/* declaration of exported global variables									*/
/*																			*/
/****************************************************************************/

/* predefined blockvector description formats */
extern const BV_DESC_FORMAT DH_bvdf;            /* bvdf for domain halfening	*/
extern const BV_DESC_FORMAT one_level_bvdf;     /* bvdf for only 1 blocklevel	*/
extern const BV_DESC_FORMAT two_level_bvdf;     /* bvdf for 2 blocklevels		*/
extern const BV_DESC_FORMAT three_level_bvdf;   /* bvdf for 3 blocklevels	*/

/****************************************************************************/
/*																			*/
/* interface functions for module grid manager								*/
/*																			*/
/****************************************************************************/

/* return values for functions returning an INT. The usual rule is: 0 ok, >0 error */
#define GM_OK                                           0
#define GM_ERROR                                        1
#define GM_FILEOPEN_ERROR                       2
#define GM_RULE_WITH_ORIENTATION        3
#define GM_RULE_WITHOUT_ORIENTATION 4
#define GM_OUT_OF_MEM                           5
#define GM_OUT_OF_RANGE                         6
#define GM_NOT_FOUND                            7
#define GM_INCONSISTANCY                        8
#define GM_COARSE_NOT_FIXED                     9
#define GM_FATAL                                        999

/* some constants passed as parameters */
#define GM_KEEP_BOUNDARY_NODES          0
#define GM_MOVE_BOUNDARY_NODES          1
#define GM_REFINE_TRULY_LOCAL           2
#define GM_COPY_ALL                             3
#define GM_REFINE_NOT_CLOSED            4
#define GM_REFINE_PARALLEL                      0
#define GM_REFINE_SEQUENTIAL            1
#define GM_REFINE_NOHEAPTEST            0
#define GM_REFINE_HEAPTEST                      1
#define GM_FCFCLL                                       1
#define GM_FFCCLL                                       2
#define GM_FFLLCC                                       3
#define GM_FFLCLC                                       4
#define GM_CCFFLL                                       5
#define GM_LOV_BEGIN                            1
#define GM_LOV_END                                      2
#define GM_GEN_FIRST                            0
#define GM_GEN_LAST                                     1
#define GM_GEN_CUT                                      2
#define GM_ALL_LEVELS                           1
#define GM_CURRENT_LEVEL                        2
#define GM_ORDER_IN_COLS                        0
#define GM_ORDER_IN_ROWS                        1
#define GM_PUT_AT_BEGIN                         1               /* put skip vectors at begin of the list	*/
#define GM_PUT_AT_END                           2               /* put skip vectors at end of the list		*/
#define GM_TAKE_SKIP                            (1<<0)
#define GM_TAKE_NONSKIP                         (1<<1)

/* get/set current multigrid, loop through multigrids */
MULTIGRID               *MakeMGItem                             (const char *name);
MULTIGRID               *GetMultigrid                           (const char *name);
MULTIGRID               *GetFirstMultigrid                      (void);
MULTIGRID               *GetNextMultigrid                       (const MULTIGRID *theMG);

/* format definition */
FORMAT                   *GetFormat                             (const char *name);
FORMAT                   *GetFirstFormat                        (void);
FORMAT                   *GetNextFormat                         (FORMAT * fmt);
INT                               ChangeToFormatDir                     (const char *name);
INT                               DeleteFormat                          (const char *name);
FORMAT                   *CreateFormat (char *name, INT sVertex, INT sMultiGrid,
                                        ConversionProcPtr PrintVertex,
                                        ConversionProcPtr PrintGrid,
                                        ConversionProcPtr PrintMultigrid,
                                        TaggedConversionProcPtr PrintVector,
                                        TaggedConversionProcPtr PrintMatrix,
                                        INT nvDesc, VectorDescriptor *vDesc,
                                        INT nmDesc, MatrixDescriptor *mDesc,
                                        SHORT ImatTypes[],
                                        INT po2t[MAXDOMPARTS][MAXVOBJECTS],
                                        INT nodeelementlist, INT edata, INT ndata);

/* create, saving and disposing a multigrid structure */
MULTIGRID *CreateMultiGrid (char *MultigridName, char *BndValProblem,
                            char *format, MEM heapSize, INT optimizedIE);
MULTIGRID   *OpenMGFromDataFile(MULTIGRID *theMG, INT number, char *type, char *DataFileName, INT Force);
MULTIGRID       *LoadMultiGrid  (char *MultigridName, char *name, char *type,
                                 char *BndValProblem, char *format,
                                 unsigned long heapSize,INT force,INT optimizedIE, INT autosave);
INT             SaveMultiGrid (MULTIGRID *theMG, char *name, char *type, char *comment, INT autosave);
INT         DisposeGrid             (GRID *theGrid);
INT             DisposeMultiGrid                (MULTIGRID *theMG);
INT         DisposeAMGLevel         (MULTIGRID *theMG);
INT         DisposeAMGLevels        (MULTIGRID *theMG);
#ifdef __TWODIM__
INT                     SaveCnomGridAndValues (MULTIGRID *theMG, char *FileName, char *plotprocName, char *tagName);
#endif

/* coarse grid manipulations */
NODE        *InsertInnerNode            (GRID *theGrid, DOUBLE *pos);
NODE        *InsertBoundaryNode     (GRID *theGrid, BNDP *bndp);

INT             DeleteNodeWithID                (GRID *theGrid, INT id);
INT             DeleteNode                              (GRID *theGrid, NODE *theNode);
ELEMENT     *InsertElementFromIDs       (GRID *theGrid, INT n, INT  *idList);
ELEMENT     *InsertElement                      (GRID *theGrid, INT n, NODE **NodeList, ELEMENT **ElemList, INT *NbgSdList);
INT         InsertMesh              (MULTIGRID *theMG, MESH *theMesh);
INT             DeleteElementWithID     (MULTIGRID *theMG, INT id);
INT             DeleteElement                   (MULTIGRID *theMG, ELEMENT *theElement);

/* refinement */
INT             EstimateHere                    (ELEMENT *theElement);
INT         MarkForRefinement       (ELEMENT *theElement, INT rule, void *data);
INT             GetRefinementMark               (ELEMENT *theElement, INT *rule, void *data);
INT             GetRefinementMarkType   (ELEMENT *theElement);
INT             RefineMultiGrid                 (MULTIGRID *theMG, INT flag, INT seq, INT mgtest);
INT         TestRefineInfo          (MULTIGRID *theMG);
INT         SetRefineInfo           (MULTIGRID *theMG);
INT             ClearMarksOnLevel               (GRID *theGrid, INT ClearType);


NODE            *GetFineNodeOnEdge              (const ELEMENT *theElement, INT side);
/*INT			GetFineSidesTouchingCoarseSide (const ELEMENT *theElement, INT side, INT *nfine, ELEMENT *Elements[MAX_SIDES_TOUCHING], INT Sides[MAX_SIDES_TOUCHING]);*/

/* moving nodes */
INT         GetMidNodeParam         (NODE * theNode, DOUBLE *lambda);
INT         GetCenterNodeParam      (NODE * theNode, DOUBLE *lambda);
#ifdef __THREEDIM__
NODE            GetSideNodeParam        (NODE * theNode, DOUBLE *lambda);
INT                     GetSideIDFromScratch    (ELEMENT *theElement, NODE *theNode);
#endif
INT         MoveMidNode             (MULTIGRID *theMG, NODE *theNode, DOUBLE lambda);
INT         MoveCenterNode          (MULTIGRID *theMG, NODE *theNode, DOUBLE *lambda);
#ifdef __THREEDIM__
INT         MoveSideNode             (MULTIGRID *theMG, NODE *theNode, DOUBLE *lambda);
#endif
INT         MoveNode                (MULTIGRID *theMG, NODE *theNode, DOUBLE *newPos);
INT             SmoothMultiGrid                 (MULTIGRID *theMG, INT niter, INT bdryFlag);
INT         SmoothGrid              (GRID *theGrid, const DOUBLE LimitLocDis, INT *MoveInfo, const INT ForceLevelSet, const INT bnd_num, const INT *bnd);
INT         SmoothGridReset         (GRID *theGrid, INT *MoveInfo);

/* handling struct blockvector_description_format (BV_DESC_FORMAT) */
INT InitBVDF                                            ( BV_DESC_FORMAT *bvdf, BLOCKNUMBER max_blocks );

/* handling struct blockvector_description (BV_DESC) */
INT PushEntry                                           ( BV_DESC *bvd, BLOCKNUMBER bnr, const BV_DESC_FORMAT *bvdf );

/* functions to create a BLOCKVECTOR structure for a regular rectangular grid */
INT CreateBVStripe2D                            ( GRID *grid, INT vectors, INT vectors_per_stripe );
INT CreateBVStripe3D                            ( GRID *grid, INT inner_vectors, INT stripes_per_plane, INT vectors_per_stripe );
INT CreateBVDomainHalfening                     ( GRID *grid, INT side, INT leaf_size );

/* general functions for BLOCKVECTOR */
INT CreateBlockvector                           ( GRID *theGrid, BLOCKVECTOR **BVHandle );
INT CreateBlockvector_l0                        ( GRID *theGrid, BLOCKVECTOR **BVHandle, BLOCKVECTOR *insertBV, INT after);
INT DisposeBlockvector                          ( GRID *theGrid, BLOCKVECTOR *bv );
void FreeAllBV                                          ( GRID *grid );
void SetLevelnumberBV                           ( BLOCKVECTOR *bv, INT level );

/* algebraic connections */
CONNECTION      *CreateExtraConnection  (GRID *theGrid, VECTOR *from, VECTOR *to);
INT             DisposeExtraConnections (GRID *theGrid);
INT             DisposeConnectionsInGrid (GRID *theGrid);
MATRIX          *GetMatrix                              (const VECTOR *FromVector, const VECTOR *ToVector);
MATRIX      *GetOrderedMatrix       (const VECTOR *FromVector, const VECTOR *ToVector);
CONNECTION      *GetConnection                  (const VECTOR *FromVector, const VECTOR *ToVector);
#ifdef __INTERPOLATION_MATRIX__
MATRIX      *GetIMatrix             (VECTOR *FineVector, VECTOR *CoarseVector);
MATRIX      *CreateIMatrix          (GRID *theGrid, VECTOR *fvec, VECTOR *cvec);
INT                     DisposeIMatrixList              (GRID *theGrid, VECTOR *theVector);
INT             DisposeIMatricesInGrid  (GRID *theGrid);
#endif
INT         GetAllVectorsOfElement  (GRID *theGrid, ELEMENT *theElement,
                                     VECTOR **vec);

/* searching */
NODE            *FindNodeFromId                 (GRID *theGrid, INT id);
NODE            *FindNodeFromPosition   (GRID *theGrid, DOUBLE *pos, DOUBLE *tol);
VECTOR          *FindVectorFromPosition (GRID *theGrid, DOUBLE *pos, DOUBLE *tol);
ELEMENT         *FindElementFromId              (GRID *theGrid, INT id);
ELEMENT         *FindElementFromPosition(GRID *theGrid, DOUBLE *pos);
ELEMENT     *FindElementOnSurface   (MULTIGRID *theMG, DOUBLE *global);
ELEMENT     *NeighbourElement       (ELEMENT *t, INT side);
INT          InnerBoundary          (ELEMENT *t, INT side);
BLOCKVECTOR *FindBV                                     (const GRID *grid, const BV_DESC *bvd, const BV_DESC_FORMAT *bvdf );

/* list */
void            ListMultiGridHeader     (const INT longformat);
void            ListMultiGrid                   (MULTIGRID *theMG, const INT isCurrent, const INT longformat);
INT         MultiGridStatus         (MULTIGRID *theMG, INT gridflag, INT greenflag, INT lbflag, INT verbose);
void            ListGrids                               (const MULTIGRID *theMG);
void            ListNode                                (MULTIGRID *theMG, NODE *theNode,               INT dataopt, INT bopt, INT nbopt, INT vopt);
void            ListNodeSelection               (MULTIGRID *theMG,                                              INT dataopt, INT bopt, INT nbopt, INT vopt);
void            ListNodeRange                   (MULTIGRID *theMG, INT from, INT to,    INT dataopt, INT bopt, INT nbopt, INT vopt);
void            ListElement                     (MULTIGRID *theMG, ELEMENT *theElement, INT dataopt, INT bopt, INT nbopt, INT vopt);
void            ListElementSelection    (MULTIGRID *theMG,                                              INT dataopt, INT bopt, INT nbopt, INT vopt);
void            ListElementRange                (MULTIGRID *theMG, INT from, INT to,    INT dataopt, INT bopt, INT nbopt, INT vopt, INT lopt);
void            ListVector                              (MULTIGRID *theMG, VECTOR *theVector,   INT matrixopt, INT dataopt);
void            ListVectorSelection     (MULTIGRID *theMG,                                              INT matrixopt, INT dataopt);
void            ListVectorOfElementSelection(MULTIGRID *theMG,                                  INT matrixopt, INT dataopt);
void            ListVectorRange                 (MULTIGRID *theMG,                      INT fl, INT tl, INT fromV, INT toV, INT matrixopt, INT dataopt);

/* query */
LINK            *GetLink                                (NODE *from, NODE *to);
EDGE            *GetSonEdge                             (EDGE *theEdge);
#ifdef __THREEDIM__
EDGE            *FatherEdge                             (NODE **SideNodes, INT ncorners, NODE **Nodes, EDGE *theEdge);
#endif
EDGE            *GetEdge                                (NODE *from, NODE *to);
INT             GetSons                                 (ELEMENT *theElement, ELEMENT *SonList[MAX_SONS]);
#ifdef ModelP
INT             GetAllSons                              (ELEMENT *theElement, ELEMENT *SonList[MAX_SONS]);
#endif
INT             VectorPosition                  (const VECTOR *theVector, DOUBLE *position);
INT             VectorInElement                 (ELEMENT *theElement, VECTOR *theVector);
INT             MinMaxAngle                     (ELEMENT *theElement, DOUBLE *amin, DOUBLE *amax);

/* check */
#ifndef ModelP
INT                     CheckGrid                               (GRID *theGrid, INT checkgeom, INT checkalgebra, INT checklists);
#else
INT                     CheckGrid                               (GRID *theGrid, INT checkgeom, INT checkalgebra, INT checklists, INT checkif);
#endif
INT                     CheckLists                              (GRID *theGrid);

/* selection */
void            ClearSelection                  (MULTIGRID *theMG);
INT             AddNodeToSelection              (MULTIGRID *theMG, NODE *theNode);
INT             IsNodeSelected                  (MULTIGRID *theMG, NODE *theNode);
INT             AddElementToSelection   (MULTIGRID *theMG, ELEMENT *theElement);
INT             IsElementSelected               (MULTIGRID *theMG, ELEMENT *theElement);
INT             AddVectorToSelection    (MULTIGRID *theMG, VECTOR *theVector);
INT             IsVectorSelected                (MULTIGRID *theMG, VECTOR *theVector);
INT             RemoveNodeFromSelection (MULTIGRID *theMG, NODE *theNode);
INT             RemoveElementFromSelection(MULTIGRID *theMG, ELEMENT *theElement);
INT             RemoveVectorFromSelection(MULTIGRID *theMG, VECTOR *theVector);

/* multigrid user data space management (using the heaps.c block heap management) */
INT             AllocateControlEntry    (INT cw_id, INT length, INT *ce_id);
INT             FreeControlEntry                (INT ce_id);
void            ListCWofObject                  (const void *obj, INT offset);
void            ListAllCWsOfObject              (const void *obj);
void            ListAllCWsOfAllObjectTypes (PrintfProcPtr myprintf);
unsigned INT ReadCW                                     (const void *obj, INT ce);
void            WriteCW                                 (void *obj, INT ce, INT n);
void            ResetCEstatistics               (void);
void            PrintCEstatistics               (void);
INT             DefineMGUDBlock                 (BLOCK_ID id, MEM size);
INT             FreeMGUDBlock                   (BLOCK_ID id);
BLOCK_DESC      *GetMGUDBlockDescriptor (BLOCK_ID id);

/* ordering of degrees of freedom */
ALG_DEP         *CreateAlgebraicDependency              (char *name, DependencyProcPtr DependencyProc);
FIND_CUT        *CreateFindCutProc                              (char *name, FindCutProcPtr FindCutProc);
INT                     LexOrderVectorsInGrid                   (GRID *theGrid, const INT *order, const INT *sign, INT which, INT SpecSkipVecs, INT AlsoOrderMatrices);
INT             OrderVectors                                    (MULTIGRID *theMG, INT levels, INT mode, INT PutSkipFirst, INT SkipPat, const char *dependency, const char *dep_options, const char *findcut);
INT                     ShellOrderVectors                               (GRID *theGrid, VECTOR *seed);
INT                     PrepareForLineorderVectors              (GRID *theGrid);
INT                     MarkBeginEndForLineorderVectors (ELEMENT *elem, INT dt, INT ot, const INT *mark);
INT                     LineOrderVectors                                (MULTIGRID *theMG, INT levels, const char *dependency, const char *dep_options, const char *findcut, INT verboselevel);
INT                     RevertVecOrder                                  (GRID *theGrid);

/* functions for evaluation-fct management */
INT              InitEvalProc                                                           (void);
EVALUES         *CreateElementValueEvalProc                             (const char *name, PreprocessingProcPtr PreProc, ElementEvalProcPtr EvalProc);
EVECTOR         *CreateElementVectorEvalProc                            (const char *name, PreprocessingProcPtr PreProc, ElementVectorProcPtr EvalProc, INT d);
MVALUES         *CreateMatrixValueEvalProc                                      (const char *name, PreprocessingProcPtr PreProc, MatrixEvalProcPtr EvalProc);
EVALUES         *CreateElementValueEvalProcFromCoeffProc        (const char *name, CoeffProcPtr CoeffProc);
EVECTOR         *CreateElementVectorEvalProcFromCoeffProc       (const char *name, CoeffProcPtr CoeffProc, INT d);
EVALUES         *GetElementValueEvalProc                                        (const char *name);
EVECTOR         *GetElementVectorEvalProc                                       (const char *name);
MVALUES         *GetMatrixValueEvalProc                                         (const char *name);
EVALUES         *GetFirstElementValueEvalProc                           (void);
EVALUES         *GetNextElementValueEvalProc                            (EVALUES *EvalProc);
EVECTOR         *GetFirstElementVectorEvalProc                          (void);
EVECTOR         *GetNextElementVectorEvalProc                           (EVECTOR *EvecProc);

/* miscellaneous */
INT             RenumberMultiGrid                                       (MULTIGRID *theMG, INT *nboe, INT *nioe, INT *nbov, INT *niov, INT *foid, INT *non);
INT                     OrderNodesInGrid                                        (GRID *theGrid, const INT *order, const INT *sign, INT AlsoOrderLinks);
INT             PutAtEndOfList                                          (GRID *theGrid, INT cnt, ELEMENT **elemList);
INT         MGSetVectorClasses                              (MULTIGRID *theMG);
INT         SetEdgeSubdomainFromElements        (GRID *theGrid);
INT         SetSubdomainIDfromBndInfo           (MULTIGRID *theMG);
INT         FixCoarseGrid                       (MULTIGRID *theMG);
INT                     ClearMultiGridUsedFlags                         (MULTIGRID *theMG, INT FromLevel, INT ToLevel, INT mask);

#endif
