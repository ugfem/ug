// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  connectuggrape.c												*/
/*																			*/
/* Purpose:   Give Grape the information to Display ug data					*/
/*																			*/
/* Author:	  Klaus Johannsen				                                                                */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/*			  Monika Wierse/Martin Metscher									*/
/*            Universitaet Freiburg											*/
/*            Institut fuer Angewandte Mathematik							*/
/*            Hermann--Herder--Str. 10										*/
/*            D-79104 Freiburg												*/
/*																			*/
/* History:   27.04.96 begin, ug version 3.1								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

/* ug includes */
#include "defs.h"
#include "gm.h"
#include "evm.h"
#include "general.h"

/* COARSE is defined in gm.h and grape.h: we do not need it at all, so: */
#ifdef COARSE
#undef COARSE
#endif

/* GRAPE includes */
#include "grape.h"
#if (DIM == 3)

#define HELEMENT HELEMENT3D
#define HMESH HMESH3D
#define MESH MESH3D
#define HELEMENT_DESCRIPTION HELEMENT3D_DESCRIPTION
#define ELEMENTG ELEMENT3D
#define F_HDATA F_HDATA3D
#define F_HEL_INFO F_HEL_INFO3D

#elif (DIM == 2)

#define HELEMENT HELEMENT2D
#define HMESH HMESH2D
#define MESH MESH2D
#define HELEMENT_DESCRIPTION HELEMENT2D_DESCRIPTION
#define ELEMENTG ELEMENT2D
#define F_HDATA F_HDATA2D
#define F_HEL_INFO F_HEL_INFO2D

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

/* debug option */
/* #define TEST  */
/* #define DEBUG */

/*    print_element(&((d)->hel));\*/
#define FILL_ELEMENT(d,s)       {int i,k;\
                                 (d)->hel.eindex=ID(s);\
                                 (d)->hel.parent = WOL((d)->previous);\
                                 for(i=0; i<CORNERS_OF_ELEM(s); i++)\
                                 {\
                                   (d)->vindex[i]=ID(MYVERTEX(CORNER(s,i)));\
                                   ASSURE((d)->vindex[i] < (d)->hel.mesh->max_vindex," vertex index too large",return (NULL));\
                                   for(k=0; k<DIM; k++)\
                                     (d)->vertex[i][k] = CVECT(MYVERTEX(CORNER(s,i)))[k];\
                                 }\
                                 (d)->hel.user_data         = (s);\
                                 ASSURE((d)->hel.user_data, "FILL_ELEMENT: cannot allocate memory",return (NULL));\
                                 ASSURE((d)->hel.eindex < (d)->hel.mesh->max_eindex," element index too large",return (NULL));\
}

#define IL(el) ((DL_HELEM *) (el))
#define WOL(el) ((HELEMENT *) (el))
#define HMESH_MG(p)             ((MULTIGRID*)((p)->user_data))
#define ELEMENT_ELEM(p) ((ELEMENT*)((p)->user_data))
#define ELEM_CVECT(p,i)         (CVECT(MYVERTEX(CORNER(p,i)))

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct hel_in_list
{
  HELEMENT hel ;                                                                /* H ELEMENT of GRAPE		*/

  /* explicit information in hel_in_list */
  struct  hel_in_list *previous,*next ;                 /* double linked list		*/
  float vcoord[MAX_CORNERS_OF_ELEM][DIM];               /* coordinates of corners	*/
  float *vertex[MAX_CORNERS_OF_ELEM];                           /* ptr to corner vectors	*/
  int vindex[MAX_CORNERS_OF_ELEM];                              /* VERTEX ids of corners	*/

  /* management of son-elements */
  ELEMENT *actSons[MAX_SONS];                                           /* actual son				*/
  INT nextSonIndex;                                                             /* index of next son		*/

  /* for complete hirarchical procedures */
  double minmax[2];                                                             /* min/max value on ELEM	*/

} DL_HELEM;

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

/* ug: */
static EVALUES *ElemEval;
static EVECTOR *ElemVec;
static PreprocessingProcPtr EvalPreProc;
static ElementEvalProcPtr EvalPlotProc;
static ElementVectorProcPtr EvecPlotProc;
static int FirstEvec;

/* grape: */
static MANAGER *mgr;
static int nCompEntries;
static DL_HELEM *last_dl_helem = NULL;                                                          /* last DL_HELEMs of double linked list */

static HELEMENT *simplex_neighbour(HELEMENT *, int, int, double *, float *);
static void save_actuell_status(void);
static void add_methods(void) ;
static int build_list_of_data(int *);
static int  simplex_boundary( HELEMENT *, int);
static int  simplex_world_to_coord( HELEMENT *, float *, double *);
static void simplex_coord_to_world( HELEMENT *, double *, float *);
static int  simplex_check_inside( HELEMENT *, double *);

static HELEMENT *simplex_first_child(HELEMENT *);
static HELEMENT *simplex_next_child(HELEMENT *);
static HELEMENT *simplex_first_macro( HMESH *);
static ELEMENTG *simplex_next_element(ELEMENTG *el) ;

#if (DIM == 3)
/* vertex indices of the polygons for a tetrahedron   nach ug: in element.c */
static int t_v0[3] = {0,2,1},   t_v1[3] = {1,2,3};
static int t_v2[3] = {2,0,3},   t_v3[3] = {3,0,1};

/* polygon adjacencies  for a tetrahedron  nach ug */
static int t_p0[3] = {2,1,3},   t_p1[3] = {0,2,3};
static int t_p2[3] = {0,3,1},   t_p3[3] = {2,0,1};


/* local coordinates of the vertices for a tetrahedron nach ug */
static double t_c0[4] = {0.,0.,0.,1.}, t_c1[4] = {1.,0.,0.,0.};
static double t_c2[4] = {0.,1.,0.,0.}, t_c3[4] = {0.,0.,1.,0.};

static int simplex_polygon_length[4] = {3, 3, 3, 3};
static int    *simplex_vertex[4] =       {t_v0,t_v1,t_v2,t_v3};
static int    *simplex_next_polygon[4] = {t_p0,t_p1,t_p2,t_p3};
static double *simplex_coord[4] =        {t_c0,t_c1,t_c2,t_c3};

#elif (DIM == 2)

/* vertex indices of the polygons for a triangles   nach ug: dort Feld CornerOfSide */
static int t_v0[2] = {0,1},   t_v1[2] = {1,2}, t_v2[2] = {2,0};

/* polygon adjacencies  for a triangles  nach ug */
static int t_p0[2] = {2,1},   t_p1[2] = {0,2}, t_p2[2] = {0,3};


/* local coordinates of the vertices for a triangles nach ug */
static double t_c0[3] = {0.,0.,1.}, t_c1[3] = {1.,0.,0.};
static double t_c2[3] = {0.,1.,0.};

static int simplex_polygon_length[3] = {3, 3, 3};
static int    *simplex_vertex[3] =       {t_v0,t_v1,t_v2};
static int    *simplex_next_polygon[3] = {t_p0,t_p1,t_p2};
/*static double *simplex_coord[3] =        {t_c0,t_c1,t_c2};*/

/* test */
static double t_c[3][3] = {{0.,0.,1.}, {1.,0.,0.},{0.,1.,0.}};
static double *simplex_coord[3] =        {t_c[0],t_c[1],t_c[2]};

#endif

#if (DIM == 2)

static HELEMENT2D_DESCRIPTION simplex_description =
{
  /* configuration */
  (DIM+1),                                                      /* n_points								*/
  (DIM+1),                                                      /* n_polygons                                                   */
  simplex_coord,                                        /* barycentr. coords of corners			*/
  1,                                                            /* param. degree: conforming elems (?)	*/

  /* function ptrs */
  simplex_world_to_coord,                       /* Global2Local                                                 */
  simplex_coord_to_world,                       /* Local2Global							*/
  simplex_check_inside,                         /* PointInElement and more				*/
  simplex_neighbour,                                    /* find a neighbor ...					*/
  simplex_boundary                                      /* identify boundary side				*/
};

#elif (DIM == 3)

static HELEMENT3D_DESCRIPTION simplex_description =
{
  /* configuration */
  (DIM+1),                                                      /* n_points								*/
  (DIM+1),                                                      /* n_polygons                                                   */
  simplex_polygon_length,                       /* length of each polygon				*/
  simplex_vertex,                                       /* corners of polygons					*/
  simplex_next_polygon,                         /* neighbor polys of poly				*/
  DIM+1,                                                        /* dim of barycenric coords				*/
  simplex_coord,                                        /* barycentr. coords of corners			*/
  1,                                                            /* param. degree: conforming elems (?)	*/

  /* function ptrs */
  simplex_world_to_coord,                       /* Global2Local                                                 */
  simplex_coord_to_world,                       /* Local2Global							*/
  simplex_check_inside,                         /* PointInElement and more				*/
  simplex_neighbour,                                    /* find a neighbor ...					*/
  simplex_boundary                                      /* identify boundary side				*/
};

#endif

/*	static VINHERIT VinheritLookUp[3][2][4] =
        { { {{1,hi0,hd1},{1,hi2,hd1},{1,hi3,hd1},{2,hi01,hd01}},
            {{1,hi1,hd1},{1,hi3,hd1},{1,hi2,hd1},{2,hi01,hd01}} },
          { {{1,hi0,hd1},{1,hi2,hd1},{1,hi3,hd1},{2,hi01,hd01}},
            {{1,hi1,hd1},{1,hi2,hd1},{1,hi3,hd1},{2,hi01,hd01}}  },
          { {{1,hi0,hd1},{1,hi2,hd1},{1,hi3,hd1},{2,hi01,hd01}},
            {{1,hi1,hd1},{1,hi2,hd1},{1,hi3,hd1},{2,hi01,hd01}}  } };
 */

static int selectmode = 0;
static F_HDATA    *data_f_data_struct = NULL;
static F_HEL_INFO *data_f_hel_info_struct = NULL;

static BUTTON *buttfswitch;


/* RCS string */
RCSID("$Header$",UG_RCS_STRING)

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   get_dl_helem - get a DL_HELEM structure

   SYNOPSIS:
   static DL_HELEM *get_dl_helem (void);

   PARAMETERS:
   .  void - nothing

   DESCRIPTION:
   The function get a DL_HELEM structure from free list or allocates a
   new one

   RETURN VALUE:
   DL_HELEM *
      dl_helem - ptr to new DL_HELEM

   SEE ALSO:
   D*/
/****************************************************************************/

static DL_HELEM *get_dl_helem (void)
{
  INT i;
  DL_HELEM  *dl_helem;

  if (last_dl_helem!=NULL)
  {
    dl_helem = last_dl_helem;
    last_dl_helem = last_dl_helem->next;

    /* erweitern falls next  von dl_helem != NULL */
  }
  else
  {
    dl_helem = (DL_HELEM  *)mem_alloc(sizeof(DL_HELEM ));
    for(i=0; i<MAX_CORNERS_OF_ELEM; i++)
      dl_helem->vertex[i]                     = dl_helem->vcoord[i];
    dl_helem->hel.vertex                            = dl_helem->vertex;
    dl_helem->hel.vindex                            = dl_helem->vindex;
    dl_helem->hel.descr                             = &simplex_description;
    dl_helem->hel.size_of_user_data         = sizeof(ELEMENT);
    dl_helem->next                                  = NULL;
  }

  return (dl_helem);
}

/****************************************************************************/
/*D
   free_dl_helem - free DL_HELEM structure

   SYNOPSIS:
   void free_dl_helem (DL_HELEM *dl_helem);

   PARAMETERS:
   .  dl_helem - structre to free

   DESCRIPTION:
   The function frees a DL_HELEM structure from free list

   RETURN VALUE:

   SEE ALSO:
   D*/
/****************************************************************************/

static void free_dl_helem (DL_HELEM *dl_helem)
{
  dl_helem->next = last_dl_helem;
  last_dl_helem = dl_helem;

  /* erweitern falls next  von dl_helem != NULL */
}

/****************************************************************************/
/*
   SYNOPSIS:

   PARAMETERS:

   DESCRIPTION:

   RETURN VALUE:

   SEE ALSO:
 */
/****************************************************************************/

void print_element (HELEMENT *el)
{
  int i,j ;

  printf("info zu Element %d: \n",el->eindex);

  /*
     for (i=0; i< MAX_CORNERS_OF_ELEM; i++)
          printf("vertices %d \n",el->vindex[i]);

     #if (DIM == 3)
     for (i=0; i< MAX_CORNERS_OF_ELEM; i++)
          printf(" %f %f %f\n",el->vertex[i][0],el->vertex[i][1],el->vertex[i][2]);
     #elif (DIM == 2)
     for (i=0; i< MAX_CORNERS_OF_ELEM; i++)
          printf(" %f %f \n",el->vertex[i][0],el->vertex[i][1]);
     #endif
     for (i=0; i< MAX_CORNERS_OF_ELEM; i++)
          printf(" Nachbar an %d. Seite %d\n",i,ID(NBELEM(ELEMENT_ELEM(el),i)));
   */

  return;
}

/****************************************************************************/
/*D
   simplex_check_inside - check if point is inside element

   SYNOPSIS:
   int simplex_check_inside (HELEMENT *el, double *coord);

   PARAMETERS:
   .  el - HELEMENT
   .  coord - barycentr. coordinates

   DESCRIPTION:
   The function checks if point inside element, if not gives nb of
   polygon abrought which it is to find

   RETURN VALUE:
   int
      -1 - point lies inside
          >=0 - nb of polygon

   SEE ALSO:
   D*/
/****************************************************************************/

static int simplex_check_inside (HELEMENT *el, double *coord)
{
  int i;
  int iret = -1;
  double cmin = -1.0E-5;


  for (i=0; i<el->descr->number_of_vertices; i++)
    if (coord[i] < cmin)
    {
      cmin = coord[i];
      iret = i;
    }

  return(iret);
}

/****************************************************************************/
/*D
   simplex_boundary - check if side is on boundary

   SYNOPSIS:
   int simplex_boundary( HELEMENT *el, int np);

   PARAMETERS:
   .  el - HELEMENT
   .  np - nb of polygon

   DESCRIPTION:
   The function checks an element polygon lies on the boundary

   RETURN VALUE:
   int
      0 - not on the boundary
          !=0 - on the boundary

   SEE ALSO:
   D*/
/****************************************************************************/

static int simplex_boundary (HELEMENT *el, int np)
{
  ASSURE(el,"simplex_boundary: no element",return (0));

  if (OBJT(ELEMENT_ELEM(el))==IEOBJ) return (0);
  return (SIDE(ELEMENT_ELEM(el),np)!=NULL);
}

/****************************************************************************/
/*D
   simplex_free_element - free ELEMENTG from list

   SYNOPSIS:
   static void simplex_free_element (ELEMENTG *el);

   PARAMETERS:
   .  el - ELEMENTG

   DESCRIPTION:
   The function frees the DL_HELEM associated with the ELEMENTG

   RETURN VALUE:
   void

   SEE ALSO:
   D*/
/****************************************************************************/

static void simplex_free_element (ELEMENTG *el)
{
  free_dl_helem(IL((HELEMENT*)el));
  return;
}

/****************************************************************************/
/*D
   simplex_free_element - free HELEMENT from list

   SYNOPSIS:
   static void simplex_free_element (HELEMENT *el);

   PARAMETERS:
   .  el - HELEMENT

   DESCRIPTION:
   The function frees the DL_HELEM associated with the HELEMENT

   RETURN VALUE:
   void

   SEE ALSO:
   D*/
/****************************************************************************/

static HELEMENT *simplex_neighbour (HELEMENT *el, int np, int flag, double *coord, float *xyz)
{
  ELEMENT *elem;

  ASSURE(el,"simplex_neighbour: no element",return (NULL));

  elem = NBELEM(ELEMENT_ELEM(el),np);
  if (elem==NULL)
  {
    simplex_free_element((ELEMENTG*)el);
    return(NULL);
  }
  FILL_ELEMENT(IL(el),elem);

  /* Parents ??? */
  return(el);
}

/****************************************************************************/
/*D
   simplex_world_to_coord - Global2Local

   SYNOPSIS:
   static int simplex_world_to_coord (HELEMENT3D *el, float *xyz, double *coord);

   PARAMETERS:
   .  el - HELEMENT
   .  xyz - global
   .  coord - local

   DESCRIPTION:
   Global2Local

   RETURN VALUE:
   int
      -1 - if ok
          !=-1 - else

   SEE ALSO:
   D*/
/****************************************************************************/

#if (DIM == 3)

static int simplex_world_to_coord (HELEMENT3D *el, float *xyz, double *coord)
{
  double a[3][3],b[3] ;
  int i,j;
  int *vx;
  float **v ;

  ASSURE(el, "simplex_coord_to_world: no element", return (0));
  v = el->vertex ;

  /* compute the barycentric coordinates of Point */
  /* by solving the linear system:
     ( x0-x3  x1-x3  x2-x3 )  ( coord[0] )  =   ( x - x3)
     ( y0-y3  y1-y3  y2-y3 )  ( coord[1] )  =   ( y - y3)
     ( z0-z3  z1-z3  z2-z3 )  ( coord[2] )  =   ( z - z3)

     coord[3] : = 1 - coord[0] - coord[1] - coord[2]
   */

  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      a[j][i] = v[i][j] - v[3][j];

  for (j=0; j<3; j++)
    b[j] = xyz[j] - v[3][j];


  if(!g_solve3(a,b,coord))
  {
    printf("simplex_coord_to_world: ERROR\n");
    return(-2);
  }
  coord[3] = 1. - coord[0] - coord[1] - coord[2];

  return (simplex_check_inside(el,coord));
}

#elif (DIM == 2)

static int simplex_world_to_coord( HELEMENT2D *el,float *xyz, double *coord)
{
  double a[2][2],b[2] ;
  int i,j;
  int *vx;
  float **v ;

  ASSURE(el, "simplex_coord_to_world: no element", return (0));
  v = el->vertex ;

  /* compute the barycentric coordinates of Point */
  /* by solving the linear system:
     ( x0-x3  x1-x3   )  ( coord[0] )  =   ( x - x2)
     ( y0-y3  y1-y3   )  ( coord[1] )  =   ( y - y2)
     coord[2] : = 1 - coord[0] - coord[1]
   */

  for (i=0; i<2; i++)
    for (j=0; j<2; j++)
      a[j][i] =  v[i][j] - v[3][j];

  for (j=0; j<2; j++)
    b[j] = xyz[j] - v[3][j];


  if(!g_solve2(a,b,coord))
  {
    printf("simplex_coord_to_world: ERROR\n");
    return(-2);
  }
  coord[3] = 1. - coord[0] - coord[1];

  return (simplex_check_inside(el,coord));
}

#endif

/****************************************************************************/
/*D
   simplex_coord_to_world - Local2Global

   SYNOPSIS:
   static void simplex_coord_to_world (HELEMENT *el, double *coord, float *xyz);

   PARAMETERS:
   .  el - HELEMENT
   .  coord - local
   .  xyz - global

   DESCRIPTION:
   Local2Global

   RETURN VALUE:
   void

   SEE ALSO:
   D*/
/****************************************************************************/

static void simplex_coord_to_world (HELEMENT *el, double *coord, float *xyz)
{
  int i,j;
  float **v ;

  ASSURE(el,"simplex_coord_to_world: no element",return );

  v = el->vertex ;

  for (j=0; j<DIM; j++) xyz[j]=coord[0]*v[0][j];
  for (i=1; i<DIM+1; i++)
    for (j = 0; j < DIM; j++)
      xyz[j] += coord[i] * v[i][j];

  return;
}

/****************************************************************************/
/*D
   simplex_copy_element - copies element

   SYNOPSIS:
   static ELEMENTG *simplex_copy_element (ELEMENTG *e1);

   PARAMETERS:
   .  e1 - ELEMENTG

   DESCRIPTION:
   Copies ELEMENTG (ELEMENT2D resp. ELEMENT3D)

   RETURN VALUE:
   ELEMENTG
   .n    NULL if error
   .n    ptr to copy if ok

   SEE ALSO:
   D*/
/****************************************************************************/

static ELEMENTG *simplex_copy_element (ELEMENTG *e1)
{
  ELEMENTG *e2;
  int i,nop;

  nop = e1->mesh->max_number_of_vertices;
  e2 = (ELEMENTG*)mem_alloc(sizeof(ELEMENTG));
  memcpy(e2,e1,sizeof(ELEMENTG));
  e2->vertex = (float**)mem_alloc(nop*sizeof(float*));
  e2->vindex = (int *)mem_alloc(nop*sizeof(int));
  e2->user_data = (char *)mem_alloc(e2->size_of_user_data);
  memcpy(e2->vertex,e1->vertex,nop*sizeof(float*));
  memcpy(e2->user_data,e1->user_data,e2->size_of_user_data);
  memcpy(e2->vindex,e1->vindex,nop*sizeof(int));

  /* not clear if it works: check if used */
  return (NULL);
  /*return (e2);*/
}

/****************************************************************************/
/*D
   simplex_first_macro - gives first macro element

   SYNOPSIS:
   static HELEMENT *simplex_first_macro (HMESH *mesh);

   PARAMETERS:
   .  mesh - HMESH from GRAPE

   DESCRIPTION:
   gives first macro element constructed by first element on coarsest
   level of MULTIGRID from ug

   RETURN VALUE:
   HELEMENT
   .n    NULL if error
   .n    element if ok

   SEE ALSO:
   D*/
/****************************************************************************/

static HELEMENT *simplex_first_macro (HMESH *mesh)
{
  DL_HELEM *dl_helem;
  int i,j ;

  ASSURE(mesh,"simplex_macro_element: no mesh",return (NULL));

  /* Start of level list */
  dl_helem = get_dl_helem() ;
  dl_helem->hel.mesh      = (HMESH *)mesh;
  dl_helem->hel.vinh      = NULL ;
  dl_helem->previous      = NULL ;

  FILL_ELEMENT(dl_helem,FIRSTELEMENT(GRID_ON_LEVEL(HMESH_MG(mesh),0)));

#ifdef DEBUG
  printf(" in first_macro \n");
  print_element(WOL(dl_helem));
#endif

  return (WOL(dl_helem));
}

/****************************************************************************/
/*D
   simplex_next_macro - gives next macro element

   SYNOPSIS:
   static HELEMENT *simplex_next_macro(HELEMENT *el);

   PARAMETERS:
   .  el - HELEMENT

   DESCRIPTION:
   gives next macro element constructed by successor element on coarsest
   level of MULTIGRID from ug

   RETURN VALUE:
   HELEMENT
   .n    NULL if error
   .n    element if ok

   SEE ALSO:
   D*/
/****************************************************************************/

static HELEMENT *simplex_next_macro(HELEMENT *el)
{
  int i,j ;
  ASSURE(el,"simplex_next_macro: no element", return (NULL));

  if (SUCCE(ELEMENT_ELEM(el))==NULL)
  {
    free_dl_helem(IL(el));
    return(NULL) ;
  }
  FILL_ELEMENT(IL(el),SUCCE(ELEMENT_ELEM(el)));

#ifdef DEBUG
  printf(" in next_macro \n");
  print_element(el) ;
#endif

  el->vinh = NULL ;
  IL(el)->previous = NULL;

  return(el) ;
}

/****************************************************************************/
/*D
   simplex_first_element - gives first element

   SYNOPSIS:
   static ELEMENTG *simplex_first_element (HMESH *mesh);

   PARAMETERS:
   .  mesh - HMESH from GRAPE

   DESCRIPTION:
   gives first element constructed by first leaf element of MULTIGRID
   from ug

   RETURN VALUE:
   HELEMENT
   .n    NULL if error
   .n    element if ok

   SEE ALSO:
   D*/
/****************************************************************************/

static ELEMENTG *simplex_first_element (MESH *mesh)
{
  DL_HELEM *helil,*hel;

  ASSURE(mesh,"simplex_first_element: no mesh",return (NULL));
  helil = IL(simplex_first_macro((HMESH*)mesh));
  while ((hel = IL(simplex_first_child(WOL(helil)))) != NULL)
    helil = hel ;

#ifdef DEBUG
  printf("in First element\n") ;
  print_element(WOL(helil)) ;
#endif

  return((ELEMENTG *)helil) ;
}

/****************************************************************************/
/*D
   simplex_next_element - gives next element

   SYNOPSIS:
   static ELEMENTG *simplex_next_element (HELEMENT *el);

   PARAMETERS:
   .  el - HELEMENT

   DESCRIPTION:
   gives next macro element constructed by successing leaf element of
   MULTIGRID from ug

   RETURN VALUE:
   HELEMENT
   .n    NULL if error
   .n    element if ok

   SEE ALSO:
   D*/
/****************************************************************************/

static ELEMENTG *simplex_next_element (ELEMENTG *el)
{
  int level ;
  HELEMENT *ell ;

  ASSURE(el,"simplex_next_element: no element",return (NULL));


#ifdef DEBUG
  print_element((HELEMENT *)el) ;
  printf(" in next element \n") ;
#endif

  level = LEVEL(ELEMENT_ELEM(el));
  if(level == 0)
  {
    if((el = (ELEMENTG *)simplex_next_macro((HELEMENT *)el)) == NULL)
      return(NULL);
  }
  else if(IL(el)->nextSonIndex < 0)
  {
    ell = WOL(IL(el)->previous);
    IL(ell)->next = NULL ;
    free_dl_helem(IL(el)) ;
    return (simplex_next_element((ELEMENTG *)ell));
  }
  else
    el =  (ELEMENTG *)simplex_next_child((HELEMENT *)el) ;

  if( simplex_first_child((HELEMENT *)el) == NULL)                      /* we are on the finest level  */
    return((ELEMENTG*)el) ;
  else
  {
    /* Going  to the level_of_interest or the finest existing */
    while( (ell = simplex_first_child((HELEMENT *)el)) != NULL )
      el = (ELEMENTG *)ell ;
    return((ELEMENTG*)el);
  }
}

/****************************************************************************/
/*D
   simplex_first_child - gives first child element

   SYNOPSIS:
   static HELEMENT *simplex_first_child (HELEMENT *el);

   PARAMETERS:
   .  el - HELEMENT

   DESCRIPTION:
   gives first child element from element

   RETURN VALUE:
   HELEMENT
   .n    NULL if no child
   .n    element if there is a child

   SEE ALSO:
   D*/
/****************************************************************************/

static HELEMENT *simplex_first_child (HELEMENT *el)
{
  int level ;
  ELEMENT *elem;
  DL_HELEM *child;

  ASSURE(el,"simplex_first_child: no mesh",return (NULL));
  elem = ELEMENT_ELEM(el);
  if (NSONS(elem)==0 || (level =  LEVEL(elem)+1) > ((HMESH *)el->mesh)->level_of_interest)
    return(NULL);

  child = get_dl_helem() ;
  child->nextSonIndex = NSONS(elem)-2;       /* index of son to process next in 'actSons'*/
  /* will be decreased by simplex_next_child    */
  /* stop at -1	            */

  child->hel.mesh  = el->mesh;
  child->hel.vinh  = NULL ;          /* "Andern */
  child->previous  = IL(el);

  if (GetSons(elem,child->actSons))
  {
    printf("get sons errors\n");
    assert(0);
  }
  FILL_ELEMENT(child,child->actSons[NSONS(elem)-1]);

#ifdef DEBUG
  print_element(WOL(child));
#endif

  return (WOL(child));
}

/****************************************************************************/
/*D
   simplex_next_child - gives next child element

   SYNOPSIS:
   static HELEMENT *simplex_next_child (HELEMENT *el);

   PARAMETERS:
   .  el - HELEMENT

   DESCRIPTION:
   gives next child element from element

   RETURN VALUE:
   HELEMENT
   .n    NULL if last child
   .n    element if there is a next child

   SEE ALSO:
   D*/
/****************************************************************************/

static HELEMENT *simplex_next_child (HELEMENT *el)
{
  int i;

  ASSURE(el,"simplex_next_child: no mesh",return (NULL));

  if (LEVEL(ELEMENT_ELEM(el)) == 0)
    return(simplex_next_macro(el));

  if( IL(el)->nextSonIndex < 0)
  {
    free_dl_helem(IL(el)) ;
    return(NULL) ;
  }
  FILL_ELEMENT(IL(el),IL(el)->actSons[IL(el)->nextSonIndex]);
  IL(el)->nextSonIndex--;

#ifdef DEBUG
  print_element(el) ;
#endif

  return(el) ;
}

/****************************************************************************/
/*D
   data_f_hel_info - some information about data

   SYNOPSIS:
   static void data_f_hel_info (HELEMENT *el, F_HEL_INFO *f_hel_info);

   PARAMETERS:
   .  el - HELEMENT
   .  f_hel_info - F_HEL_INFO structure

   DESCRIPTION:
   ???

   RETURN VALUE:
   void

   SEE ALSO:
   D*/
/****************************************************************************/

static void data_f_hel_info (HELEMENT *el, F_HEL_INFO *f_hel_info)
{
  if (data_f_hel_info_struct)
    f_hel_info->polynomial_degree = (data_f_hel_info_struct+selectmode)->polynomial_degree;

  return;
}

/****************************************************************************/
/*D
   data_f_scalar - scalar data from actual eval proc

   SYNOPSIS:
   static void data_f_scalar (HELEMENT *el, int j, double *coord, double *val);

   PARAMETERS:
   .  el - HELEMENT
   .  j - nb of corner
   .  coord - local coordinates
   .  val - scalar value

   DESCRIPTION:
   gives scalar value on local coordinate of the element, or if local
   coordinates==NULL value at corner j

   RETURN VALUE:
   void

   SEE ALSO:
   D*/
/****************************************************************************/

static void data_f_scalar (HELEMENT *el, int j, double *coord, double *val)
{
  int dim,i ;
  COORD_VECTOR LocCoord;

  dim = el->mesh->f_data->dimension_of_value;
  if (dim!=1)
  {
    printf(" dim != 1 \n") ;
    assert(0);
  }
#if (DIM == 3)
  if (coord==NULL)
    V3_COPY(el->descr->coord[j],LocCoord)
    else
      V3_COPY(coord,LocCoord)

#elif (DIM == 2)
  if (coord==NULL)
    V2_COPY(el->descr->coord[j],LocCoord)
    else
      V2_COPY(coord,LocCoord)

#endif
      *val = (*EvalPlotProc)(ELEMENT_ELEM(el),el->vertex,LocCoord);

#ifdef DEBUG
  printf(" values in data_F %d %f \n",el->eindex,val[0]) ;
#endif

  return;
}

/****************************************************************************/
/*D
   data_f_vector - vector data from actual eval proc

   SYNOPSIS:
   static void data_f_vector (HELEMENT *el, int j, double *coord, double *val);

   PARAMETERS:
   .  el - HELEMENT
   .  j - nb of corner
   .  coord - local coordinates
   .  val - scalar value

   DESCRIPTION:
   gives vector value on local coordinate of the element, or if local
   coordinates==NULL value at corner j

   RETURN VALUE:
   void

   SEE ALSO:
   D*/
/****************************************************************************/

static void data_f_vector (HELEMENT *el, int j, double *coord, double *val)
{
  int dim,i ;
  COORD_VECTOR LocCoord;

  dim = el->mesh->f_data->dimension_of_value;

#if (DIM == 3)
  if (coord==NULL)
    V3_COPY(el->descr->coord[j],LocCoord)
    else
      V3_COPY(coord,LocCoord)

#elif (DIM == 2)
  if (coord==NIL)
    V2_COPY(el->descr->coord[j],LocCoord)
    else
      V2_COPY(coord,LocCoord)

#endif

      (*EvecPlotProc)(ELEMENT_ELEM(el),el->vertex,LocCoord,val);

#ifdef DEBUG
  printf(" values in data_F %d %f \n",el->eindex,val[0]) ;
#endif

  return;
}

/****************************************************************************/
/*D
   simplex_get_bounds - give bounds

   SYNOPSIS:
   static void simplex_get_bounds (HELEMENT *el, double *min, double *max);

   PARAMETERS:
   .  el - HELEMENT
   .  min - lower bound
   .  max - upper bound

   DESCRIPTION:
   gives lower and upper bound on element (here: trivial)

   RETURN VALUE:
   void

   SEE ALSO:
   D*/
/****************************************************************************/

static void simplex_get_bounds (HELEMENT *el, double *min, double *max)
{
  *min = -100. ;
  *max = 100. ;
}

/****************************************************************************/
/*D
   simplex_get_estimate - gives estimate

   SYNOPSIS:
   static double simplex_get_estimate (HELEMENT *hel, double *VertexEst);

   PARAMETERS:
   .  el - HELEMENT
   .  VertexEst - estimate value

   DESCRIPTION:
   estimate value >=0, <=1(for example: approximation of second derivative),
   here trivial

   RETURN VALUE:
   void

   SEE ALSO:
   D*/
/****************************************************************************/

static double simplex_get_estimate (HELEMENT *hel, double *VertexEst)
{
  return(1.);
}

/****************************************************************************/
/*D
   CallGrape - jump into grape

   SYNOPSIS:
   int CallGrape (MULTIGRID *theMG);

   PARAMETERS:
   .  theMG - MULTIGRID from ug

   DESCRIPTION:
   jump into grape

   RETURN VALUE:
   int
   .n    0 if ok
   .n    1 if error

   SEE ALSO:
   D*/
/****************************************************************************/

int CheckElementsInMG (MULTIGRID *theMG)
{
  INT i, ret;
  GRID *theGrid;
  ELEMENT *theElement;

  ret=0;
  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
#if (DIM == 3)
      if (TAG(theElement)!=TETRAHEDRON)
#elif (DIM == 2)
      if (TAG(theElement)!=TRIANGLE)
#endif
        ret = 1;
  }

  if (ret)
    UserWrite("ERROR: GRAPE cannot process all elements\n");

  return (ret);
}

int CallGrape (MULTIGRID *theMG)
{
  static HMESH *mesh;
  int i,j;
  GRAPHICDEVICE *devg ;
  CONTROLDEVICE *devc ;
  SCENE *s ;
  char *name ;

  /* check if elements of theMG can be processed */
  if (CheckElementsInMG(theMG)) return (1);

  mgr = (MANAGER *)GRAPE(Manager,"get-stdmgr") ();

  if(mesh ==  NULL)
  {

#if (DIM == 3)
    mesh = (HMESH *)GRAPE(HMesh3d,"new-instance") ("ug-data-in-mesh");

#elif (DIM == 2)
    mesh = (HMESH *)GRAPE(HMesh2d,"new-instance") ("ug-data-in-mesh");

#endif

    ASSURE(mesh, "cannot create new instance", END_METHOD(NULL));

    mesh->first_element = simplex_first_element;
    mesh->next_element = simplex_next_element;
    mesh->copy_element = simplex_copy_element;
    mesh->free_element = simplex_free_element;
    mesh->first_child = simplex_first_child;
    mesh->next_child = simplex_next_child;
    mesh->first_macro = simplex_first_macro;
    mesh->next_macro = simplex_next_macro;
    mesh->max_number_of_vertices = MAX_CORNERS_OF_ELEM;
    mesh->max_dimension_of_coord = DIM+1;

    mesh->f_data = (F_HDATA *)mem_alloc(sizeof(F_HDATA));
    ASSURE(mesh->f_data, "CallGrape : can't alloc memory for f_data",return (NULL));
    mesh->f_data->name                = "data to hierarhical grid";
    mesh->f_data->dimension_of_value  = 1;
    mesh->f_data->continuous_data     = 1;
    mesh->f_data->f                   = data_f_scalar;
    mesh->f_data->f_el_info           = data_f_hel_info;
    mesh->f_data->size_of_user_data   = 0;
    mesh->f_data->user_data           = NULL;
    mesh->f_data->last                = NULL;
    mesh->f_data->next                = NULL;
    mesh->f_data->get_bounds          = simplex_get_bounds ;
    mesh->f_data->get_estimate        = NULL;     /*simplex_get_estimate ;*/
    mesh->f_data->est_bound           = 0.1;
    mesh->size_of_user_data           = sizeof(MULTIGRID);
    mesh->user_data                   = (void*)theMG;

#if (DIM == 2)
    mesh->dimension_of_world = 2 ;
#endif

    data_f_hel_info_struct = (F_HEL_INFO *)mem_alloc(sizeof(F_HEL_INFO));
    data_f_hel_info_struct->polynomial_degree = 1;
    ASSURE(build_list_of_data(&FirstEvec),"problems to build list of function buttons",END_METHOD(NULL));

    /* setzen der ersten skalaren Funktion */
    ElemEval = GetElementValueEvalProc(ENVITEM_NAME(GetFirstElementValueEvalProc()));
    if (ElemEval == NULL) return (1);
    if ((EvalPreProc = ElemEval->PreprocessProc)!=NULL)
      if ((*EvalPreProc)(NULL,theMG)) return (1);
    if ((EvalPlotProc = ElemEval->EvalProc)==NULL) return (1);

    add_methods();
  }
  mesh->max_eindex = theMG->elemIdCounter;
  mesh->max_vindex = theMG->vertIdCounter;
  mesh->max_level = theMG->topLevel ;
  mesh->level_of_interest = theMG->currentLevel ;

  printf ("currentLevel %d maxeindex %d \n",theMG->currentLevel, mesh->max_eindex);

  s = (SCENE *)GRAPE(Scene,"new-instance") ("scene");
  ASSIGN(s->object , (INSTANCE *)mesh ) ;
  GRAPE(s,"hmesh-interactive-send") ();
  GRAPE(mgr,"handle") (s);
  save_actuell_status();

  /* Abschiessen des Graphikfensters und des Kontrolfensters: es tut so nicht ! */

  /* devg = (GRAPHICDEVICE *)GRAPE(GraphicDevice,"get-stddev")();
     GRAPE(devg,"delete")();
     devc= (CONTROLDEVICE *)GRAPE(ControlDevice,"get-stddev")();
     GRAPE(devc,"delete")(); */

  return(0);
}

/****************************************************************************/
/*D
   button_change_function - changes method by pressing a button

   SYNOPSIS:
   BUTTON *button_change_function(void);

   PARAMETERS:

   DESCRIPTION:
   changes a method by pressing a buttom

   RETURN VALUE:
   BUTTON *
   .n    button pressed

   SEE ALSO:
   D*/
/****************************************************************************/

BUTTON *button_change_function(void)
{
  MANAGER *mgr;
  BUTTON *button;
  INSTANCE *inst;
  HMESH *mesh ;
  SCENE *sc ;

  button = (BUTTON *)START_METHOD(G_CLASS | G_INSTANCE);
  mgr = (MANAGER *)GRAPE(Manager,"get-stdmgr") ();
  inst = (INSTANCE *)GRAPE(mgr,"get-current-object") ();
  ASSURE(g_has_superclass(inst, Scene),"at the moment only for scene",
         END_METHOD(NULL));
  sc = (SCENE *)inst;
  mesh = (HMESH *)sc->object;

  if(button->name[0] == 's')
  {
    /* skalare Gr"osse */
    mesh->f_data->f = data_f_scalar;
    mesh->f_data->dimension_of_value = 1 ;
    ElemEval = GetElementValueEvalProc(button->name+2);

    ASSURE(ElemEval, "problems with scalar functions", END_METHOD(NULL));
    if ((EvalPreProc = ElemEval->PreprocessProc)!=NULL)
      ASSURE((*EvalPreProc)(button->name+2,mesh->user_data)!=1, "problems with the functions", END_METHOD(button));
    ASSURE((EvalPlotProc = ElemEval->EvalProc), "problems with the functions", END_METHOD(button));
  }
  else
  {
    mesh->f_data->f = data_f_vector;

    ElemVec = GetElementVectorEvalProc(button->name+2);
    ASSURE(ElemVec, "problems with vector functions", END_METHOD(NULL));
    if ((EvalPreProc = ElemVec->PreprocessProc)!=NULL)
      ASSURE((*EvalPreProc)(button->name+2,mesh->user_data)!=1, "problems with the functions",END_METHOD(button));
    ASSURE((EvecPlotProc = ElemVec->EvalProc), "problems with the functions",
           END_METHOD(button));
    mesh->f_data->dimension_of_value = ElemVec->dimension;
  }
  END_METHOD(button);
}

/****************************************************************************/
/*D
   add_methods - add method

   SYNOPSIS:
   static void add_methods (void);

   PARAMETERS:

   DESCRIPTION:

   RETURN VALUE:
   void

   SEE ALSO:
   D*/
/****************************************************************************/

static void add_methods (void)
{
  PROJECT *pro;

  extern BUTTON *button_change_function(void);
  GRAPE(Button,"add-method") ("change-function",button_change_function);

  return;
}

/****************************************************************************/
/*D
   build_list_of_data - build list of data

   SYNOPSIS:
   static int build_list_of_data (int *theFirstEvec);

   PARAMETERS:

   DESCRIPTION:
   set element evaluation functionx and the necessary lists for the scalar
   and vectoriell case


   RETURN VALUE:
   int
   .n    1 if ok

   SEE ALSO:
   D*/
/****************************************************************************/


static BUTTON *bs[20],*bv[20] ;
static int build_list_of_data (int *theFirstEvec)
{
  GROUP *group ;
  RADIO *radio ;
  int i,lengthlist ;
  EVALUES *theEval;
  EVECTOR *theEvec;
  char string[20] ;

  mgr = (MANAGER *)GRAPE(Manager,"get-stdmgr") ();
  radio = (RADIO *)GRAPE(Radio, "new-instance") ("");

  /* Einh"angen der Listen zur Ansteuerung verschiedener L"osungsfunktionen auf dem Gitter */

  /* Liste f"ur skalare Gr"ossen */
  i = 0;
  for (theEval=GetFirstElementValueEvalProc(); theEval!=NULL; theEval=GetNextElementValueEvalProc(theEval))
  {
    strcpy(string,"s_");
    strcat(string,ENVITEM_NAME(theEval));
    bs[i]=(BUTTON *)new_item(Button, I_Name, string,I_Method, "change-function", I_Self,I_Size, 5.5, 1.0, I_Pos, 0.0, i+1.0, I_End);
    GRAPE(radio, "add-inter") (bs[i]);
    i++ ;
  }

  (*theFirstEvec) = i;
  for (theEvec=GetFirstElementVectorEvalProc(); theEvec!=NULL; theEvec=GetNextElementVectorEvalProc(theEvec))
  {
    strcpy(string,"v_");
    strcat(string,ENVITEM_NAME(theEvec));
    bv[i]=(BUTTON *)new_item(Button, I_Name, string,I_Method, "change-function", I_Self,I_Size, 5.5, 1.0, I_Pos, 6.0, i+1.0, I_End);
    GRAPE(radio, "add-inter") (bv[i]);
    i ++;
  }
  lengthlist = (*theFirstEvec < (i-(*theFirstEvec)) ) ? (i-(*theFirstEvec)) : (*theFirstEvec);
  ASSURE(lengthlist < 19,"not enough buttons for list of functions",return (0)) ;

  GRAPE(radio,"set-pref-size") (12.0,(double)lengthlist);
  /*	group = (GROUP *)(	new_item(Group, I_Border, bfBorder|bfTitle,I_Size, 12.,
                                          (double)lengthlist+1.3,I_Pos, 0.0, (double)lengthlist+1.3, I_Item, radio,
                                          I_Item,new_item(StaticText,
                                          I_Text, "scalar data <-----", I_Size, 5.75, 1.,
                                          I_Pos, 0., 1.0, I_End),
                                          I_Item,new_item(StaticText,
                                          I_Text, "---->vectoriell data:", I_Size, 5.75, 1.,
                                          I_Pos, 5.75, 1.0, I_End),I_End)); */

  /*
     group = (GROUP *)(	new_item(Group, I_Border, bfBorder|bfTitle,
                                          I_Size, 12., 3.+1.3, I_Pos, 0.0, 3+1.3, I_Item, radio,
                                          I_Item,new_item(StaticText,
                                          I_Text, "scalar data <-:", I_Size, 5.75, 1.,
                                          I_Pos, 0., 1.0, I_End),
                                          I_Item,new_item(StaticText,
                                          I_Text, "->vectoriell data:", I_Size, 5.75, 1.,
                                          I_Pos, 5.75, 1.0, I_End),I_End));
   */

  GRAPE(radio,"set-pref-pos") (0.,3.0);
  /*GRAPE(group,"set-label")("functions");*/
  GRAPE(mgr, "add-inter") (radio);

  return(1) ;
}

/****************************************************************************/
/*D
   save_actuell_status - save actual status

   SYNOPSIS:
   static void save_actuell_status (void)

   PARAMETERS:
   .  void

   DESCRIPTION:
   save actual status of the grape.

   RETURN VALUE:
   void
   D*/
/****************************************************************************/

static void save_actuell_status (void)
{
  GRAPHICDEVICE *devg ;
  int ver = 1;
  if(devg = (GRAPHICDEVICE *) GRAPE (GraphicDevice,"get-stddev") ())
  {
    XDR *xdrp;
    if(xdrp = g_xdr_open_file("default.st", omWrite))
    {
      if (!g_xdr_version(xdrp, &ver) || !g_xdr_graphicdevice(xdrp, devg))
        ALERT(0, "Error writing instance.", ;);
      g_xdr_close_file(xdrp);
    }
    else
      ALERT(0, "Error opening file.", ;);
  }
  else g_infobox("write st","No GraphicDevice for writing.");

  return;
}

/****************************************************************************/
/*D
   InitGrape - Init grape

   SYNOPSIS:
   void usleep (unsigned long time);

   PARAMETERS:
   .  time

   DESCRIPTION:
   usleep dummy.

   RETURN VALUE:
   D*/
/****************************************************************************/

void usleep (unsigned long time)
{
  return;
}

/****************************************************************************/
/*D
   InitGrape - Init grape

   SYNOPSIS:
   INT InitGrape (void)

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function initializes the grape.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error
   D*/
/****************************************************************************/

INT InitGrape (void)
{
  return (0);
}
