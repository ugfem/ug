// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  dif.h															*/
/*																			*/
/* Purpose:   header file for domain interface                                                          */
/*																			*/
/* Author:	  Klaus Johannsen                                                                                               */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: klaus@ica3.uni-stuttgart.de							*/
/*																			*/
/* History:   04.06.96 begin												*/
/*																			*/
/* Remarks:                                                                                                                             */
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

#ifndef __LGM_TR__
#define __LGM_TR__

#ifndef __DOMAIN__
#include "domain.h"
#endif

#include "heaps.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#ifdef _2
#define LGM_TRDIM                                                       2
#endif

#ifdef _3
#define LGM_TRDIM                                                       3
#endif

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

#if (LGM_TRDIM==2)
/****************************************************************************/
/*																			*/
/*	2D structures															*/
/*																			*/
/****************************************************************************/

struct lgm_domain_info {

  /* parameters */
  char Name[128];                                               /* name of domain								*/
  char ProblemName[128];                                /* problems name of domain						*/
  int Dimension;                                                /* dimension of domain							*/
  int Convex;                                                           /* 0 (no) or 1 (yes)							*/
  float Radius, MidPoint[3];                            /* sphere of which domain is a subset			*/
  int nSubDomain;                                               /* nb. of subdomains							*/
  int nPolyline;                                                /* nb. of lines									*/
  int nPoint;                                                           /* nb. of points								*/
};

struct lgm_sizes {

  int *Subdom_nLine;                                            /* nb. of lines of each subdomain				*/
  int *Polyline_nPoint;                                 /* nb. of points of each polyline				*/
};

struct lgm_subdomain_info {

  char Unit[128];                                               /* Unit-identification							*/
  int *LineNumber;                                              /* ids of the lines                                                     */
};

struct lgm_line_info {

  int left, right;                                              /* subdomain on left and right side				*/
  int *point;                                                           /* global ids of the points                                     */
};

struct lgm_point_info {

  float position[LGM_TRDIM];                            /* position of corner							*/
};

struct lgm_mesh_info {

  int dummy;
};

typedef struct lgm_domain_info LGM_DOMAIN_INFO;
typedef struct lgm_sizes LGM_SIZES;
typedef struct lgm_subdomain_info LGM_SUBDOMAIN_INFO;
typedef struct lgm_line_info LGM_LINE_INFO;
typedef struct lgm_point_info LGM_POINT_INFO;
typedef struct lgm_mesh_info LGM_MESH_INFO;

#endif


#if (LGM_TRDIM==3)
/****************************************************************************/
/*																			*/
/*	3D structures															*/
/*																			*/
/****************************************************************************/

struct lgm_domain_info {

  char Name[128];                                               /* name of domain								*/
  char ProblemName[128];                                /* problems name of domain						*/
  int Dimension;                                                /* dimension of domain							*/
  int Convex;                                                           /* 0 (no) or 1 (yes)							*/
  float Radius, MidPoint[3];                            /* sphere of which domain is a subset			*/
  int nSubDomain;                                               /* nb. of subdomains							*/
  int nSurface;                                                 /* nb. of surfaces								*/
  int nPolyline;                                                        /* nb. of lines								*/
  int nPoint;                                                           /* nb. of points								*/
};

struct lgm_sizes {

  int *Subdom_nSurf;                                            /* nb. of surfaces of each subdomain			*/
  int *Surf_nPolyline;                                  /* nb. of polylines of each surface				*/
  int *Surf_nTriangle;                                  /* nb. of triangles of each surface				*/
  int *Surf_nPoint;                                             /* nb. of points of each surface				*/
  int *Polyline_nPoint;                                 /* nb. of points of each line					*/
};

struct lgm_subdomain_info {

  char Unit[128];                                               /* Unit-identification							*/
  int *SurfaceNumber;                                           /* ids of the surfaces                                                  */
};

struct lgm_triangle_info {

  int corner[3];                                                /* local (w.r.t. the surface) ids of the corners*/
  int neighbor[3];                                              /* local ids of the neighborsfor that surface	*/
};

struct lgm_surface_info {

  int left, right;                                              /* subdomain on left and right side				*/
  int nTriangles;                                               /* nb. of triangles representing the surface    */
  int nPoint;                                                           /* nb. of points used for this surface			*/
  int nLine;                                                            /* nb. of lines used for this surface			*/
  /*	LGM_SURFACE_DATA *SurfaceData;*/		/* data for surface								*/
  struct lgm_triangle_info *Triangle;           /* ptr to list of triangle_info					*/
  int *point;                                                           /* ptr to array of corner (global) id's			*/
  int *line;                                                            /* ptr to array of (global) line id's			*/
  int **point_list;                                             /* ptr to array for finding neighbour triangles */
  int length;                                                           /* length of point_list array					*/
};

struct lgm_line_info {

  int *point;                                                           /* global ids of the points                                     */
};

struct lgm_point_info {

  float position[LGM_TRDIM];                            /* position of corner							*/
};

struct lgm_mesh_info {

  int nBndP;                         /* nb. of boundary points              */
  int *BndP_nSurf;                   /* nb. of surfaces per bound. point    */
  int **BndP_SurfID;                 /* id of each surface                  */
  int **BndP_Cor_TriaID;                 /* id of corr. triangle of each surface*/
  float ***BndP_lcoord;              /* local coord of BndP on each surface */
  float **BndPosition;                           /* list of boundary points	            */
  int nInnP;                         /* nb. of inner nodes                  */
  float **InnPosition;               /* positions of inner nodes            */
  int nSubDomains;                   /* nb. of subdomains                   */
  int *nSides;                       /* nb. of boundary sides per subdomain */
  int **Side_corners;                /* nb. of side corners                 */
  int ***Side_corner_ids;                /* corner ids                          */
  int *nElements;                    /* nb. of element corners              */
  int **Element_corners;             /* nb. of element corners              */
  int ***Element_corner_ids;         /* nb. of side corners                 */
  int ***nbElements;                 /* nb. of side corners                 */
};

typedef struct lgm_domain_info LGM_DOMAIN_INFO;
typedef struct lgm_sizes LGM_SIZES;
typedef struct lgm_subdomain_info LGM_SUBDOMAIN_INFO;
typedef struct lgm_triangle_info LGM_TRIANGLE_INFO;
typedef struct lgm_surface_info LGM_SURFACE_INFO;
typedef struct lgm_line_info LGM_LINE_INFO;
typedef struct lgm_point_info LGM_POINT_INFO;
typedef struct lgm_mesh_info LGM_MESH_INFO;

#endif

/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

int LGM_ReadDomain                      (HEAP *theHeap, char *filename, LGM_DOMAIN_INFO *domain_info, INT MarkKey);
int LGM_ReadSizes                       (LGM_SIZES *lgm_sizes);
int LGM_ReadSubDomain           (int i, LGM_SUBDOMAIN_INFO *subdom_info);
int LGM_ReadLines                       (int i, LGM_LINE_INFO *line_info);
int LGM_ReadPoints                      (LGM_POINT_INFO *lgm_point_info);

#if (LGM_TRDIM==3)
int LGM_ReadSurface             (int i, LGM_SURFACE_INFO *surface_info);
#endif

#if (LGM_TRDIM==3)
int HGM_ReadSurface             (int i, LGM_SURFACE_INFO *surface_info);
#endif

FILE *LGM_WriteOpenFile         (char* name);
INT InitLGMTransfer             (void);

#endif
