// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  mgio.h														*/
/*																			*/
/* Purpose:   header file for mgio.c		                                                                */
/*																			*/
/* Author:	  Klaus Johannsen                                                                                               */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   13.11.96 begin												*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/* switch */
#define __MGIO_USE_IN_UG__

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __MGIO__
#define __MGIO__

#include <stdio.h>


/****************************************************************************/
/*																			*/
/* configuration of interface                                                                                           */
/*																			*/
/****************************************************************************/

#define __MGIO_USE_IN_UG__
#define MGIO_DIM                        3

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#ifdef __MGIO_USE_IN_UG__

        #include "heaps.h"
        #include "gm.h"
        #include "rm.h"
        #include "domain.h"

        #define MGIO_BNDP                                                       BNDP
        #define MGIO_BNDS                                                       BNDS

        #undef MGIO_DIM
        #define MGIO_DIM                                                        DIM
        #define MGIO_MAX_SONS_OF_ELEM                           MAX_SONS
        #define MGIO_MAX_EDGES_OF_ELEM                          MAX_EDGES_OF_ELEM
        #define MGIO_MAX_CORNERS_OF_ELEM                        MAX_CORNERS_OF_ELEM
        #define MGIO_MAX_SIDES_OF_ELEM                          MAX_SIDES_OF_ELEM
        #define MGIO_MAX_NEW_CORNERS                            MAX_NEW_CORNERS_DIM
        #define MGIO_MAX_CORNERS_OF_SIDE            MAX_CORNERS_OF_SIDE

        #define MGIO_MAXLEVEL                                           MAXLEVEL
        #define MGIO_TAGS                                                       TAGS

#else

        #define MGIO_MAXLEVEL                                           32
        #define MGIO_TAGS                                                       8

        #if (MGIO_DIM==2)
                #define MGIO_MAX_SONS_OF_ELEM                   4
                #define MGIO_MAX_EDGES_OF_ELEM                  4
                #define MGIO_MAX_CORNERS_OF_ELEM                4
                #define MGIO_MAX_SIDES_OF_ELEM                  4
                #define MGIO_MAX_NEW_CORNERS                    5
                #define MGIO_MAX_CORNERS_OF_SIDE        2
        #endif
        #if (MGIO_DIM==3)
                #define MGIO_MAX_SONS_OF_ELEM                   12
                #define MGIO_MAX_EDGES_OF_ELEM                  12
                #define MGIO_MAX_CORNERS_OF_ELEM                8
                #define MGIO_MAX_SIDES_OF_ELEM                  6
                #define MGIO_MAX_NEW_CORNERS                    19
                #define MGIO_MAX_CORNERS_OF_SIDE        4
        #endif

#endif

/* defines for MGIO_MG_GENERAL */
#define MGIO_NAMELEN                            128
#define MGIO_NODEVECTOR                         1
#define MGIO_EDGEVECTOR                         2
#define MGIO_SIDEVECTOR                         4
#define MGIO_ELEMVECTOR                         8

/* macros for MGIO_MG_GENERAL */
#define MGIO_SET_VECTYPE(p,i)           (p)->VectorTypes=(i)
#define MGIO_ADD_VECTYPE(p,i)           (p)->VectorTypes|=(i)
#define MGIO_IS_VECTYPE(p,i)            (p)->VectorTypes==(i)
#define MGIO_CONT_VECTYPE(p,i)          (p)->VectorTypes&(i)

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

struct mgio_mg_general {

  /* information about the file */
  int mode;                                     /* macros see above								*/

  /* number of objects */
  int nLevel;                                   /* nb of levels of the mg						*/
  int nNode;                                    /* nb of nodes, i.e. corners of elements		*/
  int nPoint;                                   /* nb of points, i.e. diff node locations in mg	*/
  int nElement;                         /* nb of elements in mg							*/

  /* information on geometry */
  int dim;                                                              /* dimension						*/
  char DomainName[MGIO_NAMELEN];                /* name of domain in ug				*/

  /* information on algebraic structure */
  char Formatname[MGIO_NAMELEN];                /* name of format used in ug		*/
  int VectorTypes;                                              /* macros see above					*/
};

struct mgio_ge_general {

  int nGenElement;                      /* nb of general elements used					*/
};

struct mgio_ge_element {

  int tag;                                              /* identification of general element in ug	*/
  int nCorner;                                  /* nb of corners							*/
  int nEdge;                                            /* nb of edges								*/
  int nSide;                                            /* nb of sides								*/
  int CornerOfEdge[MGIO_MAX_EDGES_OF_ELEM][2];          /* corners of edge		*/
  int CornerOfSide[MGIO_MAX_SIDES_OF_ELEM][4];          /* corners of side		*/
};

struct mgio_rr_general {

  int nRules;                                                                           /* nb of rules used							*/
};

struct mgio_sondata {
  SHORT tag;                                            /* which element type is the son                */
  SHORT corners[MGIO_MAX_CORNERS_OF_ELEM];              /* corners of the son                           */
  SHORT nb[MGIO_MAX_SIDES_OF_ELEM];                     /* neighbors of this son                        */
  /* < 20 if neighbor has same father           */
  /* >= 20 if neighbor has other father		*/
  INT path;                                             /* path used in GetSons() for ug                        */
};

struct mgio_rr_rule {

  int class;                                            /* class of rule:3bits for COPY, IREG, REG      */
  int nsons;                                            /* number of sons rule creates                  */
  int pattern[MGIO_MAX_NEW_CORNERS];                    /* stores which edges are refined               */
  int sonandnode[MGIO_MAX_NEW_CORNERS][2];              /* for each new node the number of the son      */
  /* and the local node number of the node      */
  struct mgio_sondata sons[MGIO_MAX_SONS_OF_ELEM];       /* for all new sons					*/
};

struct mgio_cg_general {

  int nPoint;                                                           /* nb of points on coarse grid, sum of next two		*/
  int nBndPoint;                                                /* nb of bnd points on coarse grid					*/
  int nInnerPoint;                                              /* nb of inner points on coarse grid				*/
  int nElement;                                                 /* nb of elements on coarse grid, sum of next two	*/
  int nBndElement;                                              /* nb of bnd elements on coarse grid				*/
  int nInnerElement;                                            /* nb of inner elements on coarse grid				*/
};

struct mgio_cg_point {

  double position[MGIO_DIM];                            /* position of the point							*/
};

struct mgio_movedcorner {

  int id;                                                               /* local id of moved corner							*/
  double position[MGIO_DIM];                            /* position of the point							*/
};

struct mgio_cg_element {

  int ge;                                                                               /* id of general element					*/
  int cornerid[MGIO_MAX_CORNERS_OF_ELEM];               /* ids of corners							*/
  int nbid[MGIO_MAX_SIDES_OF_ELEM];                             /* ids of neighbor elements                             */
  int nhe;                                                                              /* nb of he_elements of this element		*/
};

struct mgio_refinement {

  int refrule;                                                                  /* id of refinement rule					*/
  int nnewcorners;                                                              /* nb of new corners on next level			*/
  int newcornerid[MGIO_MAX_CORNERS_OF_ELEM+MGIO_MAX_NEW_CORNERS];      /* ids of new corners    */
  int nmoved;                                                                           /* nmoved new corners moved					*/
  struct mgio_movedcorner mvcorner[MGIO_MAX_NEW_CORNERS];       /* array of moved corners		*/
};

struct mgio_bd_general {

  int nBndP;                                                                            /* n BNDP in mg, only ug					*/
  int nBndS;                                                                            /* n BNDS in mg, only ug					*/
};

typedef struct mgio_mg_general MGIO_MG_GENERAL;
typedef struct mgio_ge_general MGIO_GE_GENERAL;
typedef struct mgio_ge_element MGIO_GE_ELEMENT;
typedef struct mgio_rr_general MGIO_RR_GENERAL;
typedef struct mgio_rr_rule MGIO_RR_RULE;
typedef struct mgio_cg_general MGIO_CG_GENERAL;
typedef struct mgio_cg_point MGIO_CG_POINT;
typedef struct mgio_cg_element MGIO_CG_ELEMENT;
typedef struct mgio_refinement MGIO_REFINEMENT;
typedef struct mgio_bd_general MGIO_BD_GENERAL;

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

/* read functions */
int             Read_OpenFile           (char *filename);
int             Read_MG_General         (MGIO_MG_GENERAL *mg_general);
int             Read_GE_General         (MGIO_GE_GENERAL *ge_general);
int             Read_GE_Elements        (int n, MGIO_GE_ELEMENT *ge_element);
int             Read_RR_General         (MGIO_RR_GENERAL *rr_general);
int             Read_RR_Rules           (int n, MGIO_RR_RULE    *rr_rules);
int     Read_CG_General         (MGIO_CG_GENERAL *cg_general);
int             Read_CG_Points          (int n, MGIO_CG_POINT   *cg_point);
int             Read_CG_Elements        (int n, MGIO_CG_ELEMENT *cg_element);
int     Read_Refinement         (int n, MGIO_REFINEMENT *refinement);
int             Read_BD_General         (MGIO_BD_GENERAL *bd_general);

/* write functions */
int             Write_OpenFile          (char *filename);
int             Write_MG_General        (MGIO_MG_GENERAL *mg_general);
int             Write_GE_General        (MGIO_GE_GENERAL *ge_general);
int             Write_GE_Elements       (int n, MGIO_GE_ELEMENT *ge_element);
int             Write_RR_General        (MGIO_RR_GENERAL *rr_general);
int             Write_RR_Rules          (int n, MGIO_RR_RULE    *rr_rules);
int     Write_CG_General        (MGIO_CG_GENERAL *cg_general);
int             Write_CG_Points         (int n, MGIO_CG_POINT   *cg_point);
int             Write_CG_Elements       (int n, MGIO_CG_ELEMENT *cg_element);
int     Write_Refinement        (int n, MGIO_REFINEMENT *refinement);
int             Write_BD_General        (MGIO_BD_GENERAL *bd_general);

#ifdef __MGIO_USE_IN_UG__
int             Read_PBndDesc           (BVP *theBVP, HEAP *theHeap, int n, BNDP **BndPList);
int             Write_PBndDesc          (int n, BNDP **BndPList);
#endif

/* general functions */
int     CloseFile                       ();
int     MGIO_Init                       ();

#endif
