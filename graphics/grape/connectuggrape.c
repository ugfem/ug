// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*
 * File:	  connectuggrape.c
 *
 * Purpose:   Give Grape the information to Display ug data
 *
 * Author:	  Monika Wierse
 *                           Universitaet Freiburg
 *                           Institut fuer Angewandte Mathematik
 *                           Hermann--Herder--Str. 10
 *                           D-79104 Freiburg
 *                Klaus Johannsen
 *			  Institut fuer Computeranwendungen III
 *			  Universitaet Stuttgart
 *			  Pfaffenwaldring 27
 *			  70550 Stuttgart
 *                         phone: 0049-(0)711-685-7003
 *			  fax  : 0049-(0)711-685-7000
 *
 * History:   9.10.95
 *
 * Erledigen: manager- Steruerung? Listen von Funktionen
 * ncompEntries f"uellen, tuen traces?, copy, neighbour f"ullen
 ****************************************************************************/


#include "grape.h"
#include "defs.h"
#include "gm.h"
#include "evm.h"
#include "simplex.h"


#if (DIM == 3)

#define HELEMENT HELEMENT3D
#define HMESH HMESH3D
#define HELEMENT_DESCRIPTION HELEMENT3D_DESCRIPTION
#define ELEMENTG ELEMENT3D
#define F_HDATA F_HDATA3D
#define F_HEL_INFO F_HEL_INFO3D

#elif (DIM == 2)

#define HELEMENT HELEMENT2D
#define HMESH HMESH2D
#define HELEMENT_DESCRIPTION HELEMENT2D_DESCRIPTION
#define ELEMENTG ELEMENT2D
#define F_HDATA F_HDATA2D
#define F_HEL_INFO F_HEL_INFO2D

#endif

/* #define TEST
 #define DEBUG */

/***************************** TYPEDEFS **********************************************/
/*                                                                                   */
/*************************************************************************************/

typedef struct hel_in_list
{
  HELEMENT hel ;
  struct  hel_in_list *previous,*next ;
  float vcoord[MAX_CORNERS_OF_ELEM][DIM];
  float *vertex[MAX_CORNERS_OF_ELEM];
  int vindex[MAX_CORNERS_OF_ELEM];
  double minmax[2];
  ELEMENT *actSons[MAX_SONS];
  INT nextSonIndex;
} HEL_IN_LIST;


/********************************** MACROS   ************************************/
/*                                                                              */
/********************************************************************************/
/*    print_element(&((d)->hel));\*/
#define FILL_ELEMENT(d,s)       {int i,k;\
                                 (d)->hel.eindex=ID(s);\
                                 for(i=0; i<CORNERS_OF_ELEM(s); i++)\
                                 {\
                                   (d)->vindex[i]=ID(MYVERTEX(CORNER(s,i)));\
                                   ASSURE((d)->vindex[i] < (d)->hel.mesh->max_vindex," vertex index too large",return (NULL));\
                                   for(k=0; k<DIM; k++)\
                                     (d)->vertex[i][k] = CVECT(MYVERTEX(CORNER(s,i)))[k];\
                                 }\
                                 (d)->hel.parent = WOL((d)->previous);\
                                 (d)->hel.user_data         = (s);\
                                 ASSURE((d)->hel.user_data, "FILL_ELEMENT: cannot allocate memory",return (NULL));\
                                 ASSURE((d)->hel.eindex < (d)->hel.mesh->max_eindex," element index too large",return (NULL));\
}

#define IL(el) ((HEL_IN_LIST *) (el))
#define WOL(el) ((HELEMENT *) (el))
#define HMESH_MG(p)             ((MULTIGRID*)((p)->user_data))
#define ELEMENT_ELEM(p) ((ELEMENT*)((p)->user_data))
#define ELEM_CVECT(p,i)         (CVECT(MYVERTEX(CORNER(p,i)))

/********************************** statische Variablen **************************/
/*                                                                               */
/*********************************************************************************/

/* ug: */

static EVALUES *ElemEval;
static PreprocessingProcPtr EvalPreProc;
static ElementEvalProcPtr EvalPlotProc;

/* grape: */
static MANAGER *mgr;
static int nCompEntries;
static HEL_IN_LIST *helfree = NULL;

static HELEMENT *tetra_neighbour(HELEMENT *, int, int,
                                 double *, float *);
static int  tetra_boundary( HELEMENT *, int);
static int  tetra_world_to_coord( HELEMENT *, float *, double *);
static void tetra_coord_to_world( HELEMENT *, double *, float *);
static int  tetra_check_inside( HELEMENT *, double *);

static HELEMENT *tetra_first_child(HELEMENT *);
static HELEMENT *tetra_next_child(HELEMENT *);
static HELEMENT *tetra_first_macro( HMESH *);
static ELEMENTG *tetra_next_element(HELEMENT *el) ;

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

static int tetra_polygon_length[4] = {3, 3, 3, 3};
static int    *tetra_vertex[4] =       {t_v0,t_v1,t_v2,t_v3};
static int    *tetra_next_polygon[4] = {t_p0,t_p1,t_p2,t_p3};
static double *tetra_coord[4] =        {t_c0,t_c1,t_c2,t_c3};

#elif (DIM == 2)
/* vertex indices of the polygons for a triangles   nach ug: dort Feld CornerOfSide */
static int t_v0[2] = {0,1},   t_v1[2] = {1,2}, t_v2[2] = {2,0};

/* polygon adjacencies  for a triangles  nach ug */
static int t_p0[2] = {2,1},   t_p1[2] = {0,2}, t_p2[2] = {0,3};


/* local coordinates of the vertices for a triangles nach ug */
static double t_c0[3] = {0.,0.,1.}, t_c1[3] = {1.,0.,0.};
static double t_c2[3] = {0.,1.,0.};

static int tetra_polygon_length[3] = {3, 3, 3};
static int    *tetra_vertex[3] =       {t_v0,t_v1,t_v2};
static int    *tetra_next_polygon[3] = {t_p0,t_p1,t_p2};
static double *tetra_coord[3] =        {t_c0,t_c1,t_c2};
#endif

static HELEMENT_DESCRIPTION tetra_description =
{
  (DIM+1) , (DIM+1) ,  /* n_points, n_polygons */
  tetra_polygon_length,
  tetra_vertex,
  tetra_next_polygon,
  (DIM+1),  /* dim coord */
  tetra_coord,
  1,   /* param. degree */
  tetra_world_to_coord, tetra_coord_to_world, tetra_check_inside,
  tetra_neighbour,  tetra_boundary
};



/*static VINHERIT VinheritLookUp[3][2][4] =
        { { {{1,hi0,hd1},{1,hi2,hd1},{1,hi3,hd1},{2,hi01,hd01}},
            {{1,hi1,hd1},{1,hi3,hd1},{1,hi2,hd1},{2,hi01,hd01}} },
          { {{1,hi0,hd1},{1,hi2,hd1},{1,hi3,hd1},{2,hi01,hd01}},
            {{1,hi1,hd1},{1,hi2,hd1},{1,hi3,hd1},{2,hi01,hd01}}  },
          { {{1,hi0,hd1},{1,hi2,hd1},{1,hi3,hd1},{2,hi01,hd01}},
            {{1,hi1,hd1},{1,hi2,hd1},{1,hi3,hd1},{2,hi01,hd01}}  } };
 */

/******************* routines to handle the list of HELEMENT3D s *************/
/*                                                                           */
/*****************************************************************************/

static HEL_IN_LIST *get_hel_in_list()

{
  int i;
  HEL_IN_LIST  *hel_in_list;

  if(helfree) {
    hel_in_list = helfree;
    helfree = helfree->next;

    /* erweitern falls next  von hel_in_list != NULL */

  }
  else {
    hel_in_list = (HEL_IN_LIST  *)mem_alloc(sizeof(HEL_IN_LIST ));
    for(i=0; i<MAX_CORNERS_OF_ELEM; i++)
      hel_in_list->vertex[i]     =  hel_in_list->vcoord[i];
    hel_in_list->hel.vertex       =  hel_in_list->vertex;
    hel_in_list->hel.vindex       =  hel_in_list->vindex;
    hel_in_list->hel.descr        = &tetra_description;
    hel_in_list->hel.size_of_user_data = sizeof(ELEMENT);
    hel_in_list->next            = NULL;
  }
  return(hel_in_list);
}

static void free_hel_in_list( HEL_IN_LIST  *hel_in_list)
{
  hel_in_list->next = helfree;
  helfree = hel_in_list;

  /* erweitern falls next  von hel_in_list != NULL */
}


/****************************************************************************/
/*  Standard-Funktionen auf tetraeder                                       */
/****************************************************************************/

void print_element(HELEMENT *el)
{
  int i,j ;

  printf("info zu Element %d: \n",el->eindex) ;
  /*    for (i=0; i< MAX_CORNERS_OF_ELEM; i++)
     printf(" vertices %d \n",el->vindex[i]) ;
     #if (DIM == 3)
     for (i=0; i< MAX_CORNERS_OF_ELEM; i++)
     printf(" %f %f %f\n",el->vertex[i][0],el->vertex[i][1],el->vertex[i][2]);
     #elif (DIM == 2)
     for (i=0; i< MAX_CORNERS_OF_ELEM; i++)
     printf(" %f %f \n",el->vertex[i][0],el->vertex[i][1]);
     #endif
     for (i=0; i< MAX_CORNERS_OF_ELEM; i++)
     printf(" Nachbar an %d. Seite %d\n",i,ID(NBELEM(ELEMENT_ELEM(el),i))) ; */

}


static int tetra_check_inside( HELEMENT *el,  double *coord)
{
  int i;
  int iret = -1;
  double cmin = -1.0E-5;

  for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
    if (coord[i] < cmin) {
      cmin = coord[i];
      iret = i;
    }
  return(iret);
}

static int tetra_boundary( HELEMENT *el,  int np)

{
  ASSURE(el,       "tetra_boundary: no element",
         return (0));

  if (OBJT(ELEMENT_ELEM(el))==IEOBJ)
    return (0);
  return (SIDE(ELEMENT_ELEM(el),np)!=NULL);
}


static void tetra_free_element(HELEMENT *el)
{

  free_hel_in_list(IL(el));
  return;
}

static HELEMENT *tetra_neighbour(HELEMENT *el,  int np, int flag,
                                 double *coord,  float *xyz)
{
  ELEMENT *elem;

  ASSURE(el,       "tetra_neighbour: no element",
         return (NULL));

  elem = NBELEM(ELEMENT_ELEM(el),np);
  if (elem==NULL) {
    tetra_free_element(el);
    return(NULL);
  }
  FILL_ELEMENT(IL(el),elem);
  /* Parents ??? */
  return(el);
}

#if (DIM == 3)
static int tetra_world_to_coord( HELEMENT3D *el,float *xyz, double *coord)

{
  double a[3][3],b[3] ;
  int i,j;
  int *vx;
  float **v ;

  ASSURE(el, "tetra_coord_to_world: no element", return (0));
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
      a[j][i] =  v[i][j] - v[3][j];

  for (j=0; j<3; j++)
    b[j] = xyz[j] - v[3][j];


  if(!g_solve3(a,b,coord)) {
    printf("tetra_coord_to_world: ERROR\n");
    return(-2);
  }
  coord[3] = 1. - coord[0] - coord[1] - coord[2];


  return(tetra_check_inside(el,coord));
}
#elif (DIM == 2)
static int tetra_world_to_coord( HELEMENT2D *el,float *xyz, double *coord)

{
  double a[2][2],b[2] ;
  int i,j;
  int *vx;
  float **v ;

  ASSURE(el, "tetra_coord_to_world: no element", return (0));
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


  if(!g_solve2(a,b,coord)) {
    printf("tetra_coord_to_world: ERROR\n");
    return(-2);
  }
  coord[3] = 1. - coord[0] - coord[1];
  return(tetra_check_inside(el,coord));
}
#endif

static void tetra_coord_to_world( HELEMENT *el,
                                  double *coord, float *xyz)
{
  int i,j;
  float **v ;

  ASSURE(el,       "tetra_coord_to_world: no element",
         return );

  v = el->vertex ;

  for (j = 0; j < DIM; j++) xyz[j] = coord[0] * v[0][j];
  for (i = 1; i < DIM+1; i++) {
    for (j = 0; j < DIM; j++) xyz[j] += coord[i] * v[i][j];
  }

  return;
}

static ELEMENTG *tetra_copy_element( ELEMENTG *e1)
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
  printf(" not implemented yet \n") ;
}


/**********************************************************************************

   tetra_first_macro( HMESH *mesh): get the first element on level 0

**********************************************************************************/

static HELEMENT *tetra_first_macro( HMESH *mesh)
{
  HEL_IN_LIST *hel_in_list;
  int i,j ;

  ASSURE(mesh,       "tetra_macro_element: no mesh",
         return (NULL));

  /* Start of level list */

  hel_in_list = get_hel_in_list() ;
  hel_in_list->hel.mesh      = (HMESH *)mesh;
  hel_in_list->hel.vinh      = NULL ;
  hel_in_list->previous      = NULL ;

  FILL_ELEMENT(hel_in_list,FIRSTELEMENT(GRID_ON_LEVEL(HMESH_MG(mesh),0)));

#ifdef DEBUG
  printf(" in first_macro \n") ;
  print_element(WOL(hel_in_list)) ;
#endif

  return(WOL(hel_in_list)) ;
}

/**********************************************************************************

   tetra_next_macro(HELEMENT *el): running through the macro grid

 ************************************************************************+********/

static HELEMENT *tetra_next_macro(HELEMENT *el)
{
  int i,j ;
  ASSURE(el,       "tetra_next_macro: no element", return (NULL));

  if( SUCCE(ELEMENT_ELEM(el)) == NULL) {
    free_hel_in_list(IL(el)) ;
    return(NULL) ;
  }
  FILL_ELEMENT(IL(el),SUCCE(ELEMENT_ELEM(el)));

#ifdef DEBUG
  printf(" in next_macro \n") ;
  print_element(el) ;
#endif
  el->vinh = NULL ;
  IL(el)->previous = NULL ;
  return(el) ;
}

/***************************************************************************************************

   tetra_first_element( HMESH *mesh): going to the first child on the level_of_interest
                                     of the first element on the macro grid

 ****************************************************************************************************/

static ELEMENTG *tetra_first_element( HMESH *mesh)
{
  HEL_IN_LIST *helil,*hel;
  ASSURE(mesh,       "tetra_first_element: no mesh",return (NULL));
  helil = IL(tetra_first_macro(mesh));
  while( (hel = IL(tetra_first_child(WOL(helil)))) != NULL )
    helil = hel ;

  printf(" in First element \n") ;

#ifdef DEBUG
  print_element(WOL(helil)) ;
#endif
  return((ELEMENTG *)helil) ;

}

/*********************************************************************************************

   tetra_next_element(HELEMENT *el): running through the level_of_interest

*********************************************************************************************/


static ELEMENTG *tetra_next_element(HELEMENT *el)
{
  int level,i ;
  HELEMENT *ell ;

  ASSURE(el,       "tetra_next_element: no element",
         return (NULL));

  print_element(el) ;

#ifdef DEBUG
  printf(" in next element \n") ;
#endif

  level = LEVEL(ELEMENT_ELEM(el));

  if( tetra_next_child(el) == NULL) {
    if( level == 0 )   /* memory is given free in next_child */
      return(NULL);
    else{   /* going back to upper level */
      (IL(el)->previous)->next = NULL ;
      ell = WOL(IL(el)->previous) ;
      free_hel_in_list(IL(el)) ;
      return(tetra_next_element(ell));
    }
  }
  /* el is now set to  next_child, this is done in tetra_next_child */

  if( tetra_first_child(el) == NULL) /* we are on the finest level  */
    return((ELEMENTG*)el) ;
  else{   /* Going  to the level_of_interest or the finest existing */
    while( (ell = tetra_first_child(el)) != NULL )
      el = ell ;
    return((ELEMENTG*)el);
  }
}


/***********************************************************************************************

   tetra_first_child(HELEMENT *el):  gives back the first child of an element or NULL, if there                                                        is no one

 *************************************************************************************************/

static HELEMENT *tetra_first_child(HELEMENT *el)

{

  int level ;
  ELEMENT *elem;
  HEL_IN_LIST *child;

  ASSURE(el,       "tetra_first_child: no mesh",
         return (NULL));
  elem = ELEMENT_ELEM(el);
  if (NSONS(elem)==0 || (level =  LEVEL(elem)+1) > ((HMESH *)el->mesh)->level_of_interest)
    return(NULL);

  child = get_hel_in_list() ;
  child->nextSonIndex = NSONS(elem)-2; /* index of son to process next in 'actSons'*/
                                       /* will be decreased by tetra_next_child    */
                                       /* stop at -1	            */

  child->hel.mesh  = el->mesh;
  child->hel.vinh  = NULL ;    /* "Andern */
  child->previous  = IL(el);

  if(GetSons(elem,child->actSons)) {
    printf(" get sons errors \n") ;
    exit(1) ;
  }
  FILL_ELEMENT(child,child->actSons[NSONS(elem)-1]);
#ifdef DEBUG
  print_element(WOL(child)) ;
#endif

  return(WOL(child)) ;
}




/***************************************************************************************************

   tetra_next_child(HELEMENT *el):   going through the list of children of one element ;
                                  if there is no more child NULL will be returned;
                                  on macro grid the same as next_macro

 **************************************************************************************************/


static HELEMENT *tetra_next_child(HELEMENT *el)
{
  int i;

  ASSURE(el,       "tetra_next_child: no mesh",
         return (NULL));

  if( LEVEL(ELEMENT_ELEM(el)) == 0 )
    return(tetra_next_macro(el)) ;

  if( IL(el)->nextSonIndex < 0) {
    free_hel_in_list(IL(el)) ;
    return(NULL) ;
  }
  FILL_ELEMENT(IL(el),IL(el)->actSons[IL(el)->nextSonIndex]);
  IL(el)->nextSonIndex--;
#ifdef DEBUG
  print_element(el) ;
#endif
  return(el) ;
}





/*******************************function on data *********************************************/
/*			                                                                    */
/************************************************************************************************/

static int selectmode = 0;
static F_HDATA    *data_f_data_struct = NULL;
static F_HEL_INFO *data_f_hel_info_struct = NULL;

static void data_f_hel_info(HELEMENT *el, F_HEL_INFO *f_hel_info)

{
  if (data_f_hel_info_struct)
    f_hel_info->polynomial_degree =
      (data_f_hel_info_struct+selectmode)->polynomial_degree;
  return;
}


static void data_f( HELEMENT *el, int j,
                    double *coord, double *val)
{
  int dim,i ;
  COORD_VECTOR LocCoord;

  dim = data_f_data_struct[selectmode].dimension_of_value;
  if (dim!=1) {
    printf(" dim != 1 \n") ;
    exit(1) ;
  }
#if (DIM == 3)
  if (coord==NULL)
    V3_COPY(TRefCoord[j],LocCoord)
    else
      V3_COPY(coord,LocCoord)
#elif (DIM == 2)
  if (coord==NULL)
    V2_COPY(TRefCoord[j],LocCoord)
    else
      V2_COPY(coord,LocCoord)
#endif

      if (dim==1)
        val[0] = (*EvalPlotProc)(ELEMENT_ELEM(el),NULL,LocCoord);
      else
        val[0] = 0.0;

#ifdef DEBUG
  printf(" values in data_F %d %f \n",el->eindex,val[0]) ;
#endif
  return;

}


static void tetra_get_bounds(HELEMENT *el,double *min,double *max)
{
  *min = -100. ;
  *max = 100. ;
}


static double tetra_get_estimate(HELEMENT *hel,double *VertexEst)
{
  return(1.);
}

int CallGrape(MULTIGRID *theMG)

{
  static HMESH *mesh ;
  int i ;

  mgr = (MANAGER *)GRAPE(Manager,"get-stdmgr") ();

  if(mesh ==  NULL) {

#if (DIM == 3)
    mesh = (HMESH *)GRAPE(HMesh3d,"new-instance") ("ug-data-in-mesh");
#elif (DIM == 2)
    mesh = (HMESH *)GRAPE(HMesh2d,"new-instance") ("ug-data-in-mesh");
#endif

    ASSURE(mesh, "cannot create new instance", END_METHOD(NULL));

    mesh->first_element = tetra_first_element;
    mesh->next_element = tetra_next_element;
    mesh->copy_element = tetra_copy_element;
    mesh->free_element = tetra_free_element;
    mesh->first_child = tetra_first_child;
    mesh->next_child = tetra_next_child;
    mesh->first_macro = tetra_first_macro;
    mesh->next_macro = tetra_next_macro;
    mesh->max_number_of_vertices = MAX_CORNERS_OF_ELEM;
    mesh->max_dimension_of_coord = DIM+1;

    mesh->f_data = (F_HDATA *)mem_alloc(sizeof(F_HDATA));
    ASSURE(mesh->f_data, "CallGrape : can't alloc memory for f_data",
           return (NULL));
    mesh->f_data->name                = "data to hierarhical grid";
    mesh->f_data->dimension_of_value  = 1;
    mesh->f_data->continuous_data     = 1;
    mesh->f_data->f                   = data_f;
    mesh->f_data->f_el_info           = data_f_hel_info;
    mesh->f_data->size_of_user_data   = 0;
    mesh->f_data->user_data           = NULL;
    mesh->f_data->last                = NULL;
    mesh->f_data->next                = NULL;
    mesh->f_data->get_bounds          = tetra_get_bounds ;
    mesh->f_data->get_estimate        = tetra_get_estimate ;
    mesh->f_data->est_bound           = 0.1;
    mesh->size_of_user_data           = sizeof(MULTIGRID);
    mesh->user_data                   = (void*)theMG;

    /* set element evaluation function */
    ElemEval = GetElementValueEvalProc("x");
    if (ElemEval == NULL) return (1);
    if ((EvalPreProc = ElemEval->PreprocessProc)!=NULL)
      if((*EvalPreProc)(NULL,theMG))
        return (1);
    if ((EvalPlotProc = ElemEval->EvalProc)==NULL) return (1);

    data_f_hel_info_struct = (F_HEL_INFO *)mem_alloc(sizeof(F_HEL_INFO));
    data_f_hel_info_struct[0].polynomial_degree = 1;
    data_f_data_struct = (F_HDATA *)mem_alloc(sizeof(F_HDATA));
    data_f_data_struct[0].name                 = mem_alloc(100);
    strcpy(data_f_data_struct[0].name,"density") ;  /* Typ der Groesse */
    data_f_data_struct[0].continuous_data      = 1;
    data_f_data_struct[0].dimension_of_value   = 1 ;
    add_methods() ;
  }
  mesh->max_eindex = theMG->elemIdCounter;
  mesh->max_vindex = theMG->vertIdCounter;
  mesh->max_level = theMG->topLevel ;
  mesh->level_of_interest = theMG->currentLevel ;
  printf(" currentLevel %d maxeindex %d \n",theMG->currentLevel, mesh->max_eindex ) ;

  GRAPE(mgr,"handle") (mesh);
  return(0) ;
}


/* einziger globaler Zugang zu obigen static Definitionen */


static BUTTON *buttfswitch;

BUTTON *button_fswitch_fdiff()

{
  BUTTON *button;

  button = (BUTTON *)START_METHOD(G_INSTANCE);
  selectmode = (selectmode+1)%nCompEntries;
  GRAPE(button,"set-name") (data_f_data_struct[selectmode].name);
  END_METHOD(button);

}

INT InitGrape(void)
{}

HMESH3D *hmesh3d_change_level_of_int_send()

{
  static RULER *rul_level = NULL;
  MANAGER *mgr ;
  HMESH3D *mesh ;

  mesh = (HMESH3D *)START_METHOD(G_INSTANCE);
  ASSURE(mesh,
         "change-level-of-int-send : can't manage classes",
         END_METHOD(NULL));

  if(rul_level == NULL)
    rul_level = (RULER *)GRAPE(Ruler,"new-instance") (&mesh->level_of_interest,
                                                      "level-of-interest",RULER_INTEGER);

  if(mesh->level_of_interest > mesh->max_level)
    mesh->level_of_interest = mesh->max_level ;

  printf("mesh->level_of_interest %d mesh->max_level %d \n",
         mesh->level_of_interest,mesh->max_level) ;

  mgr = (MANAGER *)GRAPE(Manager,"get-stdmgr") ();
  GRAPE(mgr,"add-inter") (rul_level);

  END_METHOD(mesh);
}
/*
   HMESH2D *hmesh2d_change_level_of_int_send()

   {static RULER *rul_level = NULL;
   MANAGER *mgr ;
   HMESH2D *mesh ;

   mesh = (HMESH2D *)START_METHOD(G_INSTANCE);
   ASSURE(mesh,
         "change-level-of-int-send : can't manage classes",
         END_METHOD(NULL));

   if(rul_level == NULL)
    rul_level = (RULER *)GRAPE(Ruler,"new-instance")(&mesh->level_of_interest,
        "level-of-interest",RULER_INTEGER);

   if(mesh->level_of_interest > mesh->max_level)
   mesh->level_of_interest = mesh->max_level ;

   printf("mesh->level_of_interest %d mesh->max_level %d \n",
        mesh->level_of_interest,mesh->max_level) ;

   mgr = (MANAGER *)GRAPE(Manager,"get-stdmgr")();
   GRAPE(mgr,"add-inter")(rul_level);

   END_METHOD(mesh);
   }
 */

add_methods()
{
  HMESH3D *hmesh3d_change_level_of_int_send() ;
  /*  HMESH2D *hmesh2d_change_level_of_int_send() ; */

  GRAPE(HMesh3d,"add-method")
    ("change-level-of-int-send",hmesh3d_change_level_of_int_send);
  /*  GRAPE(HMesh2d,"add-method")
      ("change-level-of-int-send",hmesh2d_change_level_of_int_send);  */
}
