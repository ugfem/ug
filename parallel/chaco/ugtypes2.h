// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ugtypes2.h                                                    */
/*                                                                          */
/* Purpose:   defines data structure for unstructured grids  2D version     */
/*                                                                          */
/* Author:    Peter Bastian                                                 */
/*            Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen   */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            6900 Heidelberg                                               */
/*                                                                          */
/* History:   06.02.92 begin, ug version 2.0                                */
/*            27.08.92 parallel extension                                   */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __UGTYPES2__
#define __UGTYPES2__

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __UGENV__
#include "ugenv.h"
#endif

#ifndef __HEAPS__
#include "heaps.h"
#endif

#ifndef __PPIF__
#include "ppif.h"
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

/* necessary for most C runtime libraries */
#undef DOMAIN

/* parameters for unstructured grid structure */
#define DIM                     2                                         /* space dimension				*/
#define DIM1            1                                         /* space dimension minus one      */
#define SIDES           4                                         /* max number of sides of an elem */
#define SONS            4                                         /* max number of sons of an elem	*/

#define MAXLEVEL    32                    /* maximum depth of triangulation */
#define MAXOBJECTS      16                                        /* use 4 bits for object identific*/
#define MAXIF           11                                        /* number of interface lists		*/

#define MAXCOUNTER  0x00FFFFFF            /* max number of transferable objs*/

/* environment variable id */
#define DOMAIN_DIR              7
#define BOUNDARY_VAR    8
#define PROBLEM_DIR             9               /* this is a problem directory				*/
#define BNDCOND_VAR             16              /* this item is a boundary condition funct..*/
#define FORMAT_DIR              15
#define SYMBOL_VAR              10

/* boundary segment types */
#define PERIODIC        1
#define NON_PERIODIC    2

/* result codes of user supplied functions  0 = OK as usual */
#define OUT_OF_RANGE            1       /* coordinate out of range					*/
#define CANNOT_INIT_PROBLEM     1       /* configProblem could not init problem     */

/* values for where field in symbols */
#define inNodeData              1
#define inNodeDiag              2
#define inVertex                4
#define inElement               8
#define inLink                  16
#define inEdge                  32

/* values for kind field of symbol, only Real supported */
#define isDOUBLE                1
#define isOther                 2

/* interface lists */
#define IF_vertices                                       0
#define IF_original_to_slave_nodes        1
#define IF_master_to_slave_nodes          2
#define IF_original_and_master_nodes  3
#define IF_slave_to_slave_nodes           4
#define IF_slave_to_original_nodes        5
#define IF_slave_to_master_nodes          6
#define IF_neighbor_elements              7
#define IF_master_to_slave_elements   8
#define IF_slave_to_slave_elements    9
#define IF_slave_to_master_elements   10


/****************************************************************************/
/*                                                                          */
/* domain definition data structures                                        */
/*                                                                          */
/****************************************************************************/

struct domain {

  /* fields for environment directory */
  INT type;                                                       /* set to DOMAIN_DIR                  */
  INT locked;                                                 /* may not be changed or deleted      */
  struct domain *next;                /* linked list of domains             */
  struct domain *previous;            /* linked list of domains             */
  char name[NAMESIZE];                /* real name                          */
  struct boundary_segment *segments;      /* list of boundary segments          */

  /* domain variables */
  COORD xmin,xmax,ymin,ymax;              /* bounding box of this domain        */
  INT numOfSegments;                                      /* number of boundary segments		*/
  INT numOfCorners;                                       /* number of corner points			*/
} ;

typedef INT (*BndSegFuncPtr)(void *,COORD *,COORD *);

struct boundary_segment {

  /* fields for environment directory */
  INT type;                                                       /* set to BOUNDARY_VAR                */
  INT locked;                                                 /* may not be changed or deleted      */
  struct boundary_segment *next;      /* linked list of domains             */
  struct boundary_segment *previous;  /* linked list of domains             */
  char name[NAMESIZE];                /* real name                          */

  /* fields for boundary segment */
  INT left,right;                         /* number of left and right subdomain */
  INT id;                             /* unique id of that segment          */
  INT segType;                        /* segment type, see above            */
  INT from,to;                            /* numbers of start and end vertices  */
  INT resolution;                         /* measure for the curvature          */
  COORD alpha,beta;                   /* parameter interval used alpha>beta */

  BndSegFuncPtr BndSegFunc;                       /* pointer to definition function     */
  void *data;                                                     /* can be used by applic to find data */
} ;

typedef struct domain DOMAIN;
typedef struct boundary_segment BOUNDARY_SEGMENT;


/****************************************************************************/
/*                                                                          */
/* problem data structure                                                                                       */
/*                                                                          */
/****************************************************************************/

typedef INT (*BndCondProcPtr)(void *, DOUBLE *, DOUBLE *, INT *);
typedef INT (*CoeffProcPtr)(DOUBLE *, DOUBLE *);
typedef INT (*ConfigProcPtr)(void);

struct problem {

  /* fields for environment directory */
  INT type;                                             /* set to PROBLEM_DIR			            */
  INT locked;                                           /* may not be changed or deleted            */
  struct problem *next;
  struct problem *previous;         /* double linked list of environment items  */
  char name[NAMESIZE];                  /* name of that item						*/
  struct bndcond *theList;              /* entry to list of user supplied functions */

  /* fields for problem */
  INT problemID;                                /* used to identify problem type			*/
  ConfigProcPtr ConfigProblem;      /* procedure to reinitialize problem        */
  INT numOfCoefficients;                /* # of coefficient functions				*/
  CoeffProcPtr Coefficients[1];      /* coefficient functions					*/
} ;

struct bndcond {

  /* fields for environment variable */
  INT type;                                             /* set to BNDCOND_VAR			            */
  INT locked;                                           /* may not be changed or deleted            */
  struct bndcond *next;                 /* double linked list of environment items  */
  struct bndcond *previous;
  char name[NAMESIZE];                  /* name of that item						*/

  /* fields for boundary condition */
  INT id;                                               /* corresponds to boundary segment id !     */
  BndCondProcPtr BndCond;               /* function defining boundary condition		*/
} ;

typedef struct problem PROBLEM ;
typedef struct bndcond BNDCOND ;



/****************************************************************************/
/*                                                                          */
/* format definition data structure                                                                     */
/*                                                                          */
/****************************************************************************/

typedef INT (*ConversionProcPtr)(void *, char *);
typedef INT (*CombineProcPtr)(void *, void *);

struct format {

  /* fields for enironment variable */
  INT type;                                      /* one of the variable types above             */
  INT locked;                                    /* may not be changed or deleted               */
  struct format *next;       /* linked list of formats                      */
  struct format *previous;   /* linked list of formats                      */
  char name[NAMESIZE];       /* real name                                   */
  struct symbol *symbols;        /* list of symbols for that format				*/

  /* variables of format */
  INT sVertex;               /* size of vertex user data structure in bytes */
  INT sNode;                 /* size of node user data structure in bytes   */
  INT sDiag;                 /* size of diagonal user data structure in byt */
  INT sElement;              /* size of element user data structure in bytes*/
  INT sLink;                 /* size of link user data structure in bytes   */
  INT sEdge;                 /* size of edge user data structure in bytes   */

  ConversionProcPtr SaveVertex;                 /* write user data to string        */
  ConversionProcPtr SaveNode;
  ConversionProcPtr SaveDiag;
  ConversionProcPtr SaveElement;
  ConversionProcPtr SaveLink;
  ConversionProcPtr SaveEdge;

  ConversionProcPtr LoadVertex;                 /* read user data from string       */
  ConversionProcPtr LoadNode;
  ConversionProcPtr LoadDiag;
  ConversionProcPtr LoadElement;
  ConversionProcPtr LoadLink;
  ConversionProcPtr LoadEdge;

  ConversionProcPtr PrintVertex;                /* print user data to string        */
  ConversionProcPtr PrintNode;
  ConversionProcPtr PrintDiag;
  ConversionProcPtr PrintElement;
  ConversionProcPtr PrintLink;
  ConversionProcPtr PrintEdge;

  CombineProcPtr CombineVertex;                 /* combine user data in load transf */
  CombineProcPtr CombineNode;
  CombineProcPtr CombineDiag;
  CombineProcPtr CombineElement;
  CombineProcPtr CombineLink;
  CombineProcPtr CombineEdge;
} ;

typedef struct format FORMAT;

struct symbol {

  /* fields for enironment list variable */
  INT type;                                      /* one of the variable types above             */
  INT locked;                                    /* may not be changed or deleted               */
  struct symbol *next;       /* linked list of formats                      */
  struct symbol *previous;   /* linked list of formats                      */
  char name[NAMESIZE];       /* real name                                   */

  int index;                                     /* index in DOUBLE  array						*/
  int where;                                     /* where this symbol is stored, see above		*/
  int kind;                                      /* type of data associated with symbol			*/
  int offset;                                    /* offset in bytes from begin of user data stru*/
  int len;                                       /* length in bytes								*/
} ;

typedef struct symbol SYMBOL;



/****************************************************************************/
/*                                                                          */
/* unstructured grid data structures                                                                    */
/*                                                                          */
/****************************************************************************/

struct listobj {                                        /* prototype for double linked list elem*/
  unsigned INT control;                         /* object type and various flags		*/
  struct listobj *pred,*succ;                   /* double linked list					*/
  unsigned INT id;                                      /* globally unique id					*/
} ;

struct ivertex {                    /* inner vertex structure               */

  /* variables */
  unsigned INT control;                 /* object identification, various flags */
  union vertex *pred,*succ;             /* double linked list of nodes per level*/
  unsigned INT id;                      /* unique id used for load/store        */
  unsigned INT counter;                         /* needed for load transfer				*/
  struct coupling *couple_me;           /* inter processor coupling                     */

  COORD x[DIM];                         /* vertex position                      */
  COORD xi[DIM];                        /* local coordinates in father element  */

  /* pointers */
  void *data;                           /* associated user data structure       */
  union element *father;                /* father element                       */
  struct node *topnode;                         /* leave node where defect is computed  */
} ;


struct bvertex {                    /* boundary vertex structure            */

  /* variables */
  unsigned INT control;                 /* object identification, various flags */
  union vertex *pred,*succ;             /* double linked list of nodes per level*/
  unsigned INT id;                      /* unique id used for load/store        */
  unsigned INT counter;                         /* needed for load transfer				*/
  struct coupling *couple_me;           /* inter processor coupling                     */

  COORD x[DIM];                         /* vertex position                      */
  COORD xi[DIM];                        /* local coordinates in father element  */

  /* pointers */
  void *data;                           /* associated user data structure       */
  union element *father;                /* father element                       */
  struct node *topnode;                         /* leave node where defect is computed  */

  /* only boundary vertex */
  struct bndsegdesc *segment;           /* pointer to boundary segment          */
  COORD lambda[DIM1];                   /* position on boundary segment         */
  COORD zeta[DIM1];                     /* local coordinates in father edge     */
  /* Def.: l = (1-z)*l0 + z*l1            */
  INT side;                             /* father edge                          */
} ;

union vertex {                      /* only used to define pointer to vertex*/
  struct ivertex iv;
  struct bvertex bv;
} ;

struct node {                       /* level dependent part of a vertex     */

  /* variables */
  unsigned INT control;                 /* object identification, various flags */
  struct node *pred,*succ;              /* double linked list of nodes per level*/
  unsigned INT id;                      /* unique id used for load/store        */
  unsigned INT counter;                         /* needed for load transfer				*/
  struct coupling *couple_me;           /* inter processor coupling                     */

  INT index;                            /* discrete coordinates for ordering    */
  unsigned SHORT vskip;                 /* used bitwise for unkowns in myvertex */
  unsigned SHORT nskip;                 /* used bitwise for unknowns in node    */

  /* pointers */
  struct link *start;                   /* list of neighbors                    */
  struct node *father;                  /* father node or lcoord struct !       */
  struct node *son;                     /* node on finer level (NULL if none)   */
  union vertex *myvertex;               /* corresponding vertex structure       */

  /* node data */
  void *data;                           /* associated user data                 */
  void *diag;                       /* diagonal coefficient (as link data)  */
} ;


struct link {                       /* directed graph structure node/link   */

  /* variables */
  unsigned INT control;                         /* space for some flags                 */

  /* pointers */
  struct node *nbnode;                  /* destination of that link             */
  struct link *next;                    /* linked list of neighbors             */

  /* link data */
  void *diag;                           /* associated user data structure       */
} ;


struct edge {                       /* undirected edge of the grid graph    */

  /* variables */
  struct link links[2];                 /* two directed links                   */

  /* edge data */
  void *data;                           /* associated user data structure       */
} ;


struct ielement {                   /* triangle/quadrangle structure        */

  /* variables */
  unsigned INT control;                 /* object identification, various flags */
  union element *pred, *succ;           /* double linked list of elements       */
  unsigned INT id;                                      /* unique id used for load/store        */
  unsigned INT counter;                         /* needed for load transfer				*/
  struct coupling *couple_me;           /* inter processor coupling                     */
  unsigned INT lbdata1;                         /* *temporary* storage used in load     */
  unsigned INT lbdata2;                         /* balancing, error estimators, ...     */

  /* pointers */
  struct node *n[SIDES];                /* corners of that element              */
  union element *father;                /* father element on coarser grid       */
  union element *sons[SONS];            /* element tree                         */
  union element *nb[SIDES];                     /* dual graph                           */

  /* element data */
  void *data;                           /* associated user data structure       */
} ;

struct belement {                   /* element with one side on the boundary*/

  /* variables */
  unsigned INT control;                 /* object identification, various flags */
  union element *pred, *succ;           /* double linked list of elements       */
  unsigned INT id;                                      /* unique id used for load/store        */
  unsigned INT counter;                         /* needed for load transfer				*/
  struct coupling *couple_me;           /* inter processor coupling                     */
  unsigned INT lbdata1;                         /* *temporary* storage used in load     */
  unsigned INT lbdata2;                         /* balancing, error estimators, ...     */

  /* pointers */
  struct node *n[SIDES];                /* corners of that element              */
  union element *father;                /* father element on coarser grid       */
  union element *sons[SONS];            /* element tree                         */
  union element *nb[SIDES];                     /* dual graph                           */

  /* element data */
  void *data;                           /* associated user data structure       */

  /* only boundary element */
  struct elementside *side[SIDES];      /* NULL if interior side                */
} ;

union element {
  struct ielement ie;
  struct belement be;
} ;

struct elementside {

  /* variables */
  unsigned INT control;                 /* object identification, various flags */
  struct elementside *pred,*succ;       /* double linked list					*/

  struct bndsegdesc *segment;           /* side is part of that segment         */
  COORD lambda[DIM][DIM1];              /* parameter range                      */
} ;

struct bndsegdesc {                                     /* descriptor for one boundary segment  */
  unsigned INT control;                 /* object identification, various flags */
  unsigned INT id;                                      /* unique id used for load/store        */

  BOUNDARY_SEGMENT *theSegment;         /* (1) coordinate definition            */
  BNDCOND *theBoundaryCondition;        /* (2) boundary condition definition    */
} ;

struct grid {

  /* variables */
  unsigned INT control;                 /* object identification, various flags */
  INT level;                                                    /* level of that grid                   */
  INT nVert;                            /* number of vertices                   */
  INT nNode;                            /* number of nodes locally              */
  INT nElem;                                                    /* number of elements locally           */
  INT nEdge;                            /* number of edges locally              */
  INT nSide;                            /* number of element sides locally      */
  INT status;                           /* possible values see defines above    */

  /* pointers */
  union  element *elements;             /* pointer to first element             */
  union  vertex *vertices;                      /* pointer to first vertex				*/
  struct elementside *sides;                    /* pointer to first boundary side		*/
  struct node *firstNode;           /* pointer to first node                */
  struct node *lastNode;            /* pointer to last node                 */
  struct grid *coarser, *finer;         /* coarser and finer grids              */
  struct multigrid *mg;                 /* corresponding multigrid structure    */

  /* parallel extension */
  struct listhead *interfaces[MAXIF];       /* entry points on that level       */
  void *lbdata;                                         /* data for load balancing algorithm    */
} ;


struct multigrid {

  /* variables */
  unsigned INT control;                 /* object identification, various flags */
  INT status;                                                   /* possible values, see above           */
  INT vertIdCounter;                /* count objects in that multigrid      */
  INT nodeIdCounter;                /* count objects in that multigrid      */
  INT elemIdCounter;                /* count objects in that multigrid      */
  INT topLevel;                     /* depth of the element tree            */
  INT currentLevel;                                     /* level we are working on				*/
  DOMAIN *theDomain;                /* pointer to domain definition         */
  FORMAT *theFormat;                /* pointer to format definition         */
  PROBLEM *theProblem;                          /* pointer to problem definition        */
  struct bndsegdesc *segments;          /* array of combined boundary descriptio*/
  INT numOfSegments;                                    /* number of entries in the array above */
  INT numOfCorners;                                     /* number of corners in this domain     */
  Heap *theHeap;                    /* associated heap structure            */

  struct grid *grids[MAXLEVEL];     /* pointers to the grids                */

  INT nfreeObjects[MAXOBJECTS];         /* number of objects in each free list  */
  void *freeObjects[MAXOBJECTS];    /* pointer to allocated but unused objs */
  void *lbdata;                                         /* data for load balancing algorithm    */
  struct mglisthead *interfaces[MAXIF];       /* over all levels !              */
} ;


typedef struct listobj LISTOBJ;
typedef struct listobj    *LISTOBJPTR;
typedef struct listobj   **LISTOBJHANDLE;
typedef union  vertex VERTEX;
typedef struct node NODE;
typedef union  element ELEMENT;
typedef struct elementside ELEMENTSIDE;
typedef struct bndsegdesc BNDSEGDESC;
typedef struct link LINK;
typedef struct edge EDGE;
typedef struct grid GRID;
typedef struct multigrid MULTIGRID;


/****************************************************************************/
/*                                                                          */
/* unstructured grid data structures parallel extension                                         */
/*                                                                          */
/****************************************************************************/

struct coupling {                                       /* object stored on several processors	*/

  unsigned INT control;                         /* has every ug object					*/
  struct coupling *predA,*succA;        /* double linked list of atoms                  */
  struct coupling *predM,*succM;        /* double linked list for message		*/

  void *object;                                         /* pointer to atom stored on this proc  */
  struct listhead *myHead;                      /* the list I am in						*/
  unsigned INT proc;                                    /* proc where atom is stored                    */
  unsigned INT key;                                     /* to sort message						*/
  unsigned INT index;                                   /* position within the sorted message	*/
} ;


struct listhead {                                       /* head for message list (level wise)   */

  unsigned INT control;                         /* has every ug object					*/
  struct listhead *pred,*succ;          /* itself a double linked list			*/
  int listindex;                                        /* in interfaces array					*/

  unsigned INT proc;                                    /* destination processor				*/
  VChannelPtr vc;                                       /* async virtual channel to that proc	*/
  msgid fromId,toId;                                /* unique ids of messages                   */
  char *fromBuffer,*toBuffer;                   /* temporarily allocated buffers		*/
  INT itemSize;                                         /* bufferSize = nItems*itemSize			*/
  INT nItems;                                                   /* number of items in the following list*/
  DOUBLE **vecHandle;                                   /* pointer array for fast gather/scatter*/
  struct coupling *theList;                     /* first element of double linked list	*/
  /* NOTE: This is the entry point for    */
  /* this level, containing nItems.       */
  struct coupling *last;                        /* last coupling record on this level   */
  struct mglisthead *all;               /* listhead in multigrid structure      */
} ;

struct mglisthead {                                     /* head for message list (all levels)   */

  unsigned INT control;                         /* has every ug object					*/
  struct mglisthead *pred,*succ;        /* itself a double linked list			*/
  int listindex;                                        /* in interfaces array					*/

  unsigned INT proc;                                    /* destination processor				*/
  VChannelPtr vc;                                       /* async virtual channel to that proc	*/
  msgid fromId,toId;                                /* unique ids of messages                   */
  char *fromBuffer,*toBuffer;                   /* temporarily allocated buffers		*/
  INT itemSize;                                         /* bufferSize = nItems*itemSize			*/
  INT nItems;                                                   /* number of items in the following list*/
  DOUBLE **vecHandle;                                   /* pointer array for fast gather/scatter*/
  struct coupling *theList;                     /* first element of double linked list	*/
  /* NOTE: This is the entry point for    */
  /* this level, containing nItems.       */
  struct listhead *lh[MAXLEVEL];        /* listheads an all the levels          */
} ;

enum COPY_TYPE {NO_COPY,SLAVE_COPY,ORIGINAL_COPY,MASTER_COPY};
enum COUPLING_USAGE {NBE_COUPLING,ELEMENT_COUPLING,NODE_COUPLING,VERTEX_COUPLING};

typedef struct coupling COUPLING;
typedef struct listhead LISTHEAD;
typedef struct mglisthead MGLISTHEAD;

/****************************************************************************/
/*                                                                          */
/* Macro definitions                                                        */
/*                                                                          */
/*                                                                          */
/* Use of the control word:                                                 */
/*                                                                          */
/* macro name|bits  |V|N|L|E|G|M|use                                        */
/*                                                                          */
/* all objects:                                                             */
/* OBJT      |28-31 | |*|*|*|*|*|object type identification                 */
/* DUMMY     |27    | |*| |*| | |dummy objects                              */
/* TAG       |24-26 | |*|*|*|*|*|general purpose tag field                  */
/* USED      |23    | |*|*|*|*|*|object visited, leave them as you found 'em*/
/* LEVEL     |17-22 | |*| |*| | |level of a node/element (imp. for copies)  */
/* FLAG      |16    | |*|*|*|*|*|general purp.,  leave them as you found 'em*/
/*                                                                          */
/* vertices:                                                                */
/* MOVE      |0-1   |*| | | | | |vertex can be moved on a 0(1,2,3) dim subsp*/
/*                                                                          */
/* nodes:                                                                   */
/* CLASS     |0-2   | |*| | | | |class of node on current level             */
/* NCLASS    |3-5   | |*| | | | |class of node on next level                */
/* MODIFIED  |6     | |*| | |*| |at least one ELEMENT has changed during ref*/
/*                                                                          */
/* links and edges:                                                         */
/* LOFFSET   |0     | | |*| | | |position of link in links array            */
/* REFINED   |1     | | |*| | | |set to 1 if edge will be refined           */
/* EXTRA     |2     | | |*| | | |set to 1 if edge is no triangulation edge  */
/*                                                                          */
/* elements:                                                                */
/* OLDREF    |0-4   | | | |*| | |how element is currently refined           */
/* ECLASS    |5-6   | | | |*| | |element class from enumeration type        */
/* NEWREF    |7-11  | | | |*| | |how element will be refined                */
/* COARSEN   |12    | | | |*| | |true if element will be derefined          */
/* NSONS     |13-15 | | | |*| | |number of sons                             */
/*                                                                          */
/* coupling objects:                                                        */
/* foreign_type     |0-1   | | | | | | |current shared type                 */
/* new_foreign_type |2-3   | | | | | | |shared type after load transfer     */
/* transfer_type    |4-5   | | | | | | |type of load transfer               */
/* nb_cnt           |7-9   | | | | | | |this is nb nb_cnt in object (elem)  */
/* c_usage          |10-11 | | | | | | |hor. or vertical                    */
/*                                                                          */
/****************************************************************************/

/* object numbers */
#define IVOBJ 0                                                 /* inner vertex						*/
#define BVOBJ 1                                                 /* boundary vertex					*/
#define IEOBJ 2                                                 /* inner element					*/
#define BEOBJ 3                                                 /* boundary element					*/
#define EDOBJ 4                                                 /* edge object						*/
#define NDOBJ 5                                                 /* node object						*/
#define ESOBJ 6                                                 /* element side object				*/
#define GROBJ 7                                                 /* grid object						*/
#define MGOBJ 8                                                 /* multigrid object					*/

#define COOBJ 9                                                 /* coupling object					*/
#define LHOBJ 10                                                /* list head object					*/
#define MHOBJ 11                                                /* multigrid list head object		*/

/****************************************************************************/
/*                                                                          */
/* general macros                                                           */
/*                                                                          */
/****************************************************************************/

#define OBJMASK 0xF0000000
#define OBJSHIFT 28
#define OBJT(p) (((*((unsigned INT *)(p))) & OBJMASK)>>OBJSHIFT)
#define SETOBJT(p,n) *((unsigned INT *)(p)) = ((*((unsigned INT *)(p)))&(~OBJMASK))|((n)<<OBJSHIFT)

#define DUMMYMASK 0x08000000
#define DUMMYSHIFT 27
#define DUMMY(p) (((*((unsigned int *)(p))) & DUMMYMASK)>>DUMMYSHIFT)
#define SETDUMMY(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~DUMMYMASK))|(((n)<<DUMMYSHIFT)&DUMMYMASK)

#define TAGMASK 0x07000000
#define TAGSHIFT 24
#define TAG(p) (((*((unsigned int *)(p))) & TAGMASK)>>TAGSHIFT)
#define SETTAG(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~TAGMASK))|(((n)<<TAGSHIFT)&TAGMASK)

#define USEDMASK 0x00800000
#define USEDSHIFT 23
#define USED(p) ((*((unsigned int *)(p))) & USEDMASK)
#define SETUSED(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~USEDMASK))|(((n)<<USEDSHIFT)&USEDMASK)

#define LEVELMASK 0x007E0000
#define LEVELSHIFT 17
#define LEVEL(p) (((*((unsigned int *)(p))) & LEVELMASK)>>LEVELSHIFT)
#define SETLEVEL(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~LEVELMASK))|(((n)<<LEVELSHIFT)&LEVELMASK)

#define FLAGMASK 0x00010000
#define FLAGSHIFT 16
#define FLAG(p) ((*((unsigned int *)(p))) & FLAGMASK)
#define SETFLAG(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~FLAGMASK))|(((n)<<FLAGSHIFT)&FLAGMASK)

#define CTRL(p)         (*((unsigned INT *)(p)))
#define ID(p)           (((struct listobj *)(p))->id)

#define PRED(p) ((LISTOBJPTR)(p))->pred
#define SUCC(p) ((LISTOBJPTR)(p))->succ

#define COUPLE_ME(p)    ((NODE *)p)->couple_me

/****************************************************************************/
/*                                                                          */
/* counter word in vertices, nodes and elements                             */
/*                                                                          */
/* usage:                                                                   */
/*   Bit 30-31   :  current shared object type                              */
/*   Bit 28-29   :  shared object type after load transfer                  */
/*   Bit 24-27   :  son_cnt (only in elements)                              */
/*   Bit 23      :  kind of object table entry                              */
/*   Bit 0 -22   :  relative object count in load transfer for one message  */
/*                                                                          */
/****************************************************************************/

#define COUNTERWORD(p)         (((NODE *)(p))->counter)

#define get_domestic_type(p)     (((((NODE *)(p))->counter)>>30)&3)
#define set_domestic_type(p,n)   ((NODE *)(p))->counter = ((((NODE *)(p))->counter)&0x3FFFFFFF) | (((n)&3)<<30)

#define get_new_domestic_type(p)    (((((NODE *)(p))->counter)&0x30000000)>>28)
#define set_new_domestic_type(p,n)  ((NODE *)(p))->counter = ((((NODE *)(p))->counter)&0xCFFFFFFF) | (((n)&3)<<28)

#define get_son_cnt(p)     (((((NODE *)(p))->counter)&0x0F000000)>>24)
#define set_son_cnt(p,n)   ((NODE *)(p))->counter = ((((NODE *)(p))->counter)&0xF0FFFFFF) | (((n)&15)<<24)

#define get_obj_table_flag(p)     (((((NODE *)(p))->counter)&0x00800000)>>23)
#define set_obj_table_flag(p,n)   ((NODE *)(p))->counter = ((((NODE *)(p))->counter)&0xFF7FFFFF) | (((n)&1)<<23)

#define get_counter(p)     ((((NODE *)(p))->counter)&0x007FFFFF)
#define set_counter(p,n)   ((NODE *)(p))->counter = ((((NODE *)(p))->counter)&0xFF800000) | ((n)&0x007FFFFF)
#define get_dest(p)        ((((NODE *)(p))->counter)&0x007FFFFF)
#define set_dest(p,n)      ((NODE *)(p))->counter = ((((NODE *)(p))->counter)&0xFF800000) | ((n)&0x007FFFFF)

#define SLAVE(p)           (get_domestic_type(p)==SLAVE_COPY)
#define MASTER(p)          (get_domestic_type(p)==MASTER_COPY)
#define ORIGORMASTER(p)    ((get_domestic_type(p)==ORIGINAL_COPY)||(get_domestic_type(p)==MASTER_COPY))

/****************************************************************************/
/*                                                                          */
/* macros for vertices                                                      */
/*                                                                          */
/****************************************************************************/

#define MOVEMASK 0x00000003
#define MOVE(p) ((*((unsigned INT *)(p))) & MOVEMASK)
#define SETMOVE(p,n) *((unsigned INT *)(p)) = ((*((unsigned INT *)(p)))&(~MOVEMASK))|((n)&MOVEMASK)

#define PREDV(p)        (p)->iv.pred
#define SUCCV(p)        (p)->iv.succ
#define CVECT(p)        (p)->iv.x
#define LCVECT(p)       (p)->iv.xi
#define XC(p)           (p)->iv.x[0]
#define YC(p)           (p)->iv.x[1]
#define ZC(p)           (p)->iv.x[2]
#define XI(p)           (p)->iv.xi[0]
#define ETA(p)          (p)->iv.xi[1]
#define NU(p)           (p)->iv.xi[2]
#define VDATA(p)        (p)->iv.data
#define VFATHER(p)      (p)->iv.father
#define TOPNODE(p)      (p)->iv.topnode

#define BSEG(p)         (p)->bv.segment
#define PVECT(p)        (p)->bv.lambda
#define LPVECT(p)       (p)->bv.zeta
#define LAMBDA(p)       (p)->bv.lambda[0]
#define MU(p)           (p)->bv.lambda[1]
#define ZETA(p)         (p)->bv.zeta[0]
#define TAU(p)          (p)->bv.zeta[1]
#define ONSIDE(p)       (p)->bv.side


/****************************************************************************/
/*                                                                          */
/* macros for nodes                                                         */
/*                                                                          */
/****************************************************************************/

#define CLASSMASK 0x00000007
#define CLASS(p) ((*((unsigned INT *)(p))) & CLASSMASK)
#define SETCLASS(p,n) *((unsigned INT *)(p)) = ((*((unsigned INT *)(p)))&(~CLASSMASK))|((n)&CLASSMASK)

#define NCLASSMASK 0x00000038
#define NCLASSSHIFT 3
#define NCLASS(p) (((*((unsigned int *)(p))) & NCLASSMASK)>>NCLASSSHIFT)
#define SETNCLASS(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~NCLASSMASK))|(((n)<<NCLASSSHIFT)&NCLASSMASK)

#define MODIFIEDMASK 0x00000040
#define MODIFIEDSHIFT 6
#define MODIFIED(p) (((*((unsigned int *)(p))) & MODIFIEDMASK)>>MODIFIEDSHIFT)
#define SETMODIFIED(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~MODIFIEDMASK))|(((n)<<MODIFIEDSHIFT)&MODIFIEDMASK)

#define DELAYUPDATEMASK 0x00000080
#define DELAYUPDATESHIFT 7
#define DELAYUPDATE(p) (((*((unsigned int *)(p))) & DELAYUPDATEMASK)>>DELAYUPDATESHIFT)
#define SETDELAYUPDATE(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~DELAYUPDATEMASK))|(((n)<<DELAYUPDATESHIFT)&DELAYUPDATEMASK)

#define BORDERMASK 0x00000100
#define BORDERSHIFT 8
#define BORDER(p) (((*((unsigned int *)(p))) & BORDERMASK)>>BORDERSHIFT)
#define SETBORDER(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~BORDERMASK))|(((n)<<BORDERSHIFT)&BORDERMASK)

#define PREDN(p)        (p)->pred
#define SUCCN(p)        (p)->succ
#define INDEX(p)        (p)->index
#define VSKIP(p)        (p)->vskip
#define NSKIP(p)        (p)->nskip
#define START(p)        (p)->start
#define NFATHER(p)      (p)->father
#define SONNODE(p)      (p)->son
#define MYVERTEX(p) (p)->myvertex
#define NDATA(p)        (p)->data
#define NDIAG(p)        (p)->diag
#define COEFF(p)        (p)->diag

/****************************************************************************/
/*                                                                          */
/* macros for links                                                             */
/*                                                                          */
/****************************************************************************/

#define LOFFSETMASK 0x00000001
#define LOFFSET(p) ((*((unsigned INT *)(p))) & LOFFSETMASK)
#define SETLOFFSET(p,n) *((unsigned INT *)(p)) = ((*((unsigned INT *)(p)))&(~LOFFSETMASK))|((n)&LOFFSETMASK)

#define NBNODE(p)       (p)->nbnode
#define NEXT(p)         (p)->next
#define LDATA(p)        (p)->diag

#define MYEDGE(p)       ((EDGE *)((p)-LOFFSET(p)))
#define REVERSE(p)      ((p)+(1-LOFFSET(p)*2))

/****************************************************************************/
/*                                                                          */
/* macros for edges                                                             */
/*                                                                          */
/****************************************************************************/

#define REFINEDMASK 0x00000002
#define REFINEDSHIFT 1
#define REFINED(p) (((*((unsigned int *)(p))) & REFINEDMASK)>>REFINEDSHIFT)
#define SETREFINED(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~REFINEDMASK))|(((n)<<REFINEDSHIFT)&REFINEDMASK)

#define EXTRAMASK 0x00000004
#define EXTRASHIFT 2
#define EXTRA(p) (((*((unsigned int *)(p))) & EXTRAMASK)>>EXTRASHIFT)
#define SETEXTRA(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~EXTRAMASK))|(((n)<<EXTRASHIFT)&EXTRAMASK)

#define LINK0(p)        (&((p)->links[0]))
#define LINK1(p)        (&((p)->links[1]))
#define NODE0(p)        ((p)->links[0].nbnode)
#define NODE1(p)        ((p)->links[1].nbnode)
#define EDDATA(p)       (p)->data

/****************************************************************************/
/*                                                                          */
/* macros for elements                                                          */
/*                                                                          */
/****************************************************************************/

#define OLDREFMASK 0x0000001F
#define OLDREF(p) (((p)->ie.control) & OLDREFMASK)
#define SETOLDREF(p,n) (p)->ie.control = (((p)->ie.control)&(~OLDREFMASK))|((n)&OLDREFMASK)

#define NEWREFMASK 0x00000F80
#define NEWREFSHIFT 7
#define NEWREF(p) ((((p)->ie.control) & NEWREFMASK)>>NEWREFSHIFT)
#define SETNEWREF(p,n) (p)->ie.control = (((p)->ie.control)&(~NEWREFMASK))|(((n)<<NEWREFSHIFT)&NEWREFMASK)

/* REFINE and REFVAR are not used anymore ! (replaced by OLDREF) */
#define REFINEMASK 0x00000007
#define REFINE(p) ((*((unsigned INT *)(p))) & REFINEMASK)
#define SETREFINE(p,n) *((unsigned INT *)(p)) = ((*((unsigned INT *)(p)))&(~REFINEMASK))|((n)&REFINEMASK)

#define REFVARMASK 0x00000018
#define REFVARSHIFT 3
#define REFVAR(p) (((*((unsigned int *)(p))) & REFVARMASK)>>REFVARSHIFT)
#define SETREFVAR(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~REFVARMASK))|(((n)<<REFVARSHIFT)&REFVARMASK)

#define ECLASSMASK 0x00000060
#define ECLASSSHIFT 5
#define ECLASS(p) (((*((unsigned int *)(p))) & ECLASSMASK)>>ECLASSSHIFT)
#define SETECLASS(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~ECLASSMASK))|(((n)<<ECLASSSHIFT)&ECLASSMASK)

/* MARK and MARKVAR are not used anymore ! (replaced by NEWREF) */
#define MARKMASK 0x00000380
#define MARKSHIFT 7
#define MARK(p) (((*((unsigned int *)(p))) & MARKMASK)>>MARKSHIFT)
#define SETMARK(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~MARKMASK))|(((n)<<MARKSHIFT)&MARKMASK)

#define MARKVARMASK 0x00000C00
#define MARKVARSHIFT 10
#define MARKVAR(p) (((*((unsigned int *)(p))) & MARKVARMASK)>>MARKVARSHIFT)
#define SETMARKVAR(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~MARKVARMASK))|(((n)<<MARKVARSHIFT)&MARKVARMASK)

#define COARSENMASK 0x00001000
#define COARSENSHIFT 12
#define COARSEN(p) (((*((unsigned int *)(p))) & COARSENMASK)>>COARSENSHIFT)
#define SETCOARSEN(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~COARSENMASK))|(((n)<<COARSENSHIFT)&COARSENMASK)

#define NSONSMASK 0x0000E000
#define NSONSSHIFT 13
#define NSONS(p) (((*((unsigned int *)(p))) & NSONSMASK)>>NSONSSHIFT)
#define SETNSONS(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~NSONSMASK))|(((n)<<NSONSSHIFT)&NSONSMASK)

/* TAG values */
#define TRIANGLE        3
#define QUADRANGLE      4
#define TETRAHEDRON     4

#define SUCCE(p)        (p)->ie.succ
#define PREDE(p)        (p)->ie.pred
#define CORNER(p,i) (p)->ie.n[(i)]
#define EFATHER(p)      (p)->ie.father
#define SON(p,i)        (p)->ie.sons[(i)]
#define NBELEM(p,i) (p)->ie.nb[(i)]
#define EDATA(p)        (p)->ie.data
#define SIDE(p,i)       (p)->be.side[(i)]

/****************************************************************************/
/*                                                                          */
/* macros for element sides                                                 */
/*                                                                          */
/****************************************************************************/

#define SUCCS(p)        (p)->succ
#define PREDS(p)        (p)->pred
#define SEGDESC(p)      (p)->segment
#define PARAM(p,i,j) (p)->lambda[i][j]

/****************************************************************************/
/*                                                                          */
/* macros for boundary segment descriptors                                  */
/*                                                                          */
/****************************************************************************/

#define SEGID(p)                (p)->theSegment->id
#define LEFT(p)                 (p)->theSegment->left
#define RIGHT(p)                (p)->theSegment->right
#define SEGTYPE(p)              (p)->theSegment->type
#define FROM(p)                 (p)->theSegment->from
#define TO(p)                   (p)->theSegment->to
#define RES(p)                  (p)->theSegment->resolution
#define ALPHA(p)                (p)->theSegment->alpha
#define BETA(p)                 (p)->theSegment->beta
#define BNDSEGFUNC(p)   (*((p)->theSegment->BndSegFunc))
#define BNDDATA(p)          ((p)->theSegment->data)
#define BNDCONDFUNC(p)  ((p)->theBoundaryCondition->BndCond)

/****************************************************************************/
/*                                                                          */
/* macros for coupling objects                                              */
/*                                                                          */
/****************************************************************************/

#define foreign_type_mask 0x00000003
#define get_foreign_type(p) ((*((unsigned INT *)(p))) & foreign_type_mask)
#define set_foreign_type(p,n) *((unsigned INT *)(p)) = ((*((unsigned INT *)(p)))&(~foreign_type_mask))|((n)&foreign_type_mask)

#define new_foreign_type_mask 0x0000000C
#define new_foreign_type_shift 2
#define get_new_foreign_type(p) (((*((unsigned int *)(p))) & new_foreign_type_mask)>>new_foreign_type_shift)
#define set_new_foreign_type(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~new_foreign_type_mask))|(((n)<<new_foreign_type_shift)&new_foreign_type_mask)

#define transfer_mask 0x00000030
#define transfer_shift 4
#define get_transfer_type(p) (((*((unsigned int *)(p))) & transfer_mask)>>transfer_shift)
#define set_transfer_type(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~transfer_mask))|(((n)<<transfer_shift)&transfer_mask)

#define nb_cnt_mask 0x000000C0
#define nb_cnt_shift 6
#define get_nb_cnt(p) (((*((unsigned int *)(p))) & nb_cnt_mask)>>nb_cnt_shift)
#define set_nb_cnt(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~nb_cnt_mask))|(((n)<<nb_cnt_shift)&nb_cnt_mask)

#define c_usage_mask 0x00000300
#define c_usage_shift 8
#define get_c_usage(p) (((*((unsigned int *)(p))) & c_usage_mask)>>c_usage_shift)
#define set_c_usage(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~c_usage_mask))|(((n)<<c_usage_shift)&c_usage_mask)

#define if_list_mask 0x00007C00
#define if_list_shift 10
#define get_if_list(p) (((*((unsigned int *)(p))) & if_list_mask)>>if_list_shift)
#define set_if_list(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~if_list_mask))|(((n)<<if_list_shift)&if_list_mask)

#define transient_mask 0x00008000
#define transient_shift 15
#define get_transient_flag(p) (((*((unsigned int *)(p))) & transient_mask)>>transient_shift)
#define set_transient_flag(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~transient_mask))|(((n)<<transient_shift)&transient_mask)

#define nb_back_cnt_mask 0x00030000
#define nb_back_cnt_shift 16
#define get_nb_back_cnt(p) (((*((unsigned int *)(p))) & nb_back_cnt_mask)>>nb_back_cnt_shift)
#define set_nb_back_cnt(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~nb_back_cnt_mask))|(((n)<<nb_back_cnt_shift)&nb_back_cnt_mask)

#define join_third_mask 0xFFFC0000
#define join_third_shift 18
#define get_join_third(p) (((*((unsigned int *)(p))) & join_third_mask)>>join_third_shift)
#define set_join_third(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~join_third_mask))|(((n)<<join_third_shift)&join_third_mask)


/****************************************************************************/
/*                                                                          */
/* macros for listhead                                                      */
/*                                                                          */
/****************************************************************************/

#define send_flag_mask 0x00000001
#define get_send_flag(p) ((*((unsigned INT *)(p))) & send_flag_mask)
#define set_send_flag(p,n) *((unsigned INT *)(p)) = ((*((unsigned INT *)(p)))&(~send_flag_mask))|((n)&send_flag_mask)

#define recv_flag_mask 0x00000002
#define recv_flag_shift 1
#define get_recv_flag(p) (((*((unsigned int *)(p))) & recv_flag_mask)>>recv_flag_shift)
#define set_recv_flag(p,n) *((unsigned int *)(p)) = ((*((unsigned int *)(p)))&(~recv_flag_mask))|(((n)<<recv_flag_shift)&recv_flag_mask)

/****************************************************************************/
/*                                                                          */
/* macros for grids                                                             */
/*                                                                          */
/****************************************************************************/

#define GLEVEL(p)               ((p)->level)
#define GSTATUS(p)              ((p)->status)
#define FIRSTELEMENT(p) ((p)->elements)
#define FIRSTVERTEX(p)  ((p)->vertices)
#define FIRSTNODE(p)    ((p)->firstNode)
#define LASTNODE(p)             ((p)->lastNode)
#define UPGRID(p)               ((p)->finer)
#define DOWNGRID(p)             ((p)->coarser)
#define MYMG(p)                 ((p)->mg)
#define NV(p)                   ((p)->nVert)
#define NN(p)                   ((p)->nNode)
#define NT(p)                   ((p)->nElem)
#define NE(p)                   ((p)->nEdge)
#define NS(p)                   ((p)->nSide)
#define INTERFACE(p,n)  ((p)->interfaces[n])

/****************************************************************************/
/*                                                                          */
/* macros for multigrids                                                    */
/*                                                                          */
/****************************************************************************/

#define TOPLEVEL(p)                     ((p)->topLevel)
#define CURRENTLEVEL(p)         ((p)->currentLevel)
#define GRIDLEVEL(p,i)          ((p)->grids[i])
#define MGHEAP(p)                       ((p)->theHeap)

#endif
