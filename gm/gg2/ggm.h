// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ggm.h                                                         */
/*                                                                          */
/* Purpose:   header file for gg manager                                                        */
/*                                                                          */
/* Author:    Wolfgang Hoffmann, Henrik Renz-Reichert	                    */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart, Germany										*/
/*			  email: ug@ica3.uni-stuttgart.de							    */
/*																			*/
/* History:   08.03.94 begin, ug version 2.2                                */
/*                15.10.95 implemented in ug31                                  */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __GGM__
#define __GGM__

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __GM__
#include "gm.h"
#endif

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/* object types for front list, front component and independent front list	*/

/* possible orientations of FLs */
#define MATHPOS                 1
#define MATHNEG                 -1

/* macros for grid data */
#define STARTIFL(pg)    (pg)->first
#define LASTIFL(pg)             (pg)->last
#define NIFL(pg)                (pg)->nIndepFrontlist


/* macros for independent front lists */
#define SUCCIFL(ipfl)   ((ipfl)->succifl)
#define PREDIFL(ipfl)   ((ipfl)->predifl)
#define MYGRID(ipfl)    ((ipfl)->myGrid)
#define STARTFL(ipfl)   ((ipfl)->startfl)
#define LASTFL(ipfl)    ((ipfl)->lastfl)
#define NFL(ipfl)               ((ipfl)->nFrontlist)

/* macros for front lists */
#define SUCCFL(pfl)             ((pfl)->succfl)
#define PREDFL(pfl)             ((pfl)->predfl)
/* #define MYGRID(pfl)	((pfl)->myGrid)  <-- s.a. */
#define FLORIENTATION(pfl)      (pfl->orientation)
#define MYIFL(pfl)              ((pfl)->myIFL)
#define STARTFC(pfl)    ((pfl)->startfc)
#define LASTFC(pfl)             ((pfl)->lastfc)
#define NFC(pfl)                ((pfl)->nFrontcomp)

/* macros for front components */
#define SUCCFC(pfc)             ((pfc)->succfc)
#define PREDFC(pfc)             ((pfc)->predfc)
#define MYFL(pfc)               ((pfc)->myFL)
#define FRONTN(pfc)             ((pfc)->frontnode)
#define FCNGB(pfc)              ((pfc)->Neighbor)
#define FCNGBS(pfc)             ((pfc)->NeighborSide)

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/


struct frontcomp {
  unsigned INT control;                 /* object identification, various flags */
  struct frontcomp *succfc,*predfc;      /* double linked list of front comps	*/
  struct frontlist *myFL;                       /* pointer to my front list				*/
  NODE *frontnode;                                      /* ptr to corresponding front node		*/
  ELEMENT *Neighbor;
  int NeighborSide;        /*nfck*/
};



struct frontlist {
  unsigned INT control;                 /* object identification, various flags */
  struct frontlist *succfl,*predfl;      /* double linked list of front lists	*/
  GRID *myGrid;                                         /* pointer to my grid					*/
  struct indepfrontlist *myIFL;         /* pointer to my indep. front list		*/
  INT orientation;                                      /* MATHPOS or MATHNEG					*/
  INT SubdomainID;
  struct frontcomp *startfc;                    /* entry to front the component list	*/
  struct frontcomp *lastfc;                     /* entry to front the component list	*/
  long int nFrontcomp;                          /* # of front components in this list	*/
};

struct indepfrontlist {
  unsigned INT control;                                         /* object identification, various flags                 */
  struct indepfrontlist *succifl,*predifl;              /* double linked list of independent front lists*/
  GRID *myGrid;                                                                 /* pointer to my grid							*/
  struct frontlist *startfl;                                            /* entry to front list							*/
  struct frontlist *lastfl;                                             /* reach end of front list list					*/
  int nFrontlist;                                                               /* number of lists								*/
};

typedef struct frontlist FRONTLIST;
typedef struct frontcomp FRONTCOMP;
typedef struct indepfrontlist INDEPFRONTLIST;

typedef struct {

  INDEPFRONTLIST *first;                /* entry to independent frontlists		   */
  INDEPFRONTLIST *last;                 /* reach the end of independent frontlists */
  int nIndepFrontlist;                  /* number of independent lists		           */
} MG_GGDATA;

typedef struct {
  DOUBLE CheckCos;
  DOUBLE epsi;
  DOUBLE searchconst;
  DOUBLE h_global;
  INT msizecoeffno;
} GG_PARAM;

typedef struct {
  INT doanimate;
  INT doupdate;
  INT dostep;
  INT plotfront;
  INT printelem;
  INT equilateral;
  INT doedge;
  INT doangle;
  INT doEdge;
  INT doAngle;
} GG_ARG;


/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

INT                     SetFlagsfortemporaryGGObjects(INT IflObject,INT FlObject,INT FcObject);
INDEPFRONTLIST  *CreateIndepFrontList   (GRID *theGrid);
FRONTLIST               *CreateFrontList                (INDEPFRONTLIST *theIFL, INT SubdomainID);
FRONTCOMP               *CreateFrontComp                (FRONTLIST *mylist, FRONTCOMP *after, INT ncomp, NODE **NodeHandle);

INT DisposeADVfront                     (GRID *theGrid);
INT DisposeIndepFrontList       (INDEPFRONTLIST *theIFL);
INT DisposeFrontList            (FRONTLIST *theFL);
INT DisposeFrontComp            (FRONTLIST *myList, FRONTCOMP *theFC);

MG_GGDATA *GetMGdataPointer (MULTIGRID *theMG);

#endif
