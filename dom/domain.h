// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      domain.h                                                      */
/*                                                                          */
/* Purpose:   standard header file template                                 */
/*                                                                          */
/* Author:      klaus Johannsen                                             */
/*              Institut fuer Computeranwendungen III                       */
/*              Universitaet Stuttgart                                      */
/*              Pfaffenwaldring 27                                          */
/*              70550 Stuttgart                                             */
/*              email: ug@ica3.uni-stuttgart.de                             */
/*                                                                          */
/* History:   18.06.96 begin, ug version 3.2                                */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


/* RCS_ID
   $Header$
 */

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __DOMAIN__
#define __DOMAIN__

#include <stdio.h>

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __UGENV__
#include "ugenv.h"
#endif

#ifndef __HEAPS__
#include "heaps.h"
#endif

#ifdef __MWCW__
#include "MWCW.cmdlinedefs"
#endif

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*          compile time constants defining static data size (i.e. arrays)  */
/*          other constants                                                 */
/*          macros                                                          */
/*                                                                          */
/****************************************************************************/

#ifdef _2
#ifdef _3
#error ****    define EITHER dimension _2 OR _3       ****
#endif
#define __TWODIM__
#define DIM 2
#define DIM_OF_BND 1
#endif

#ifdef _3
#define __THREEDIM__
#define DIM 3
#define DIM_OF_BND 2
#endif

#ifndef _2
#ifndef _3
#error ****    define at least dimension two OR three        ****
#endif
#endif

/* offset of additional parameters for boundary condition function call */
enum DOM_IN_PARAMS {

  DOM_GLB_X,                                            /* general bnd cond. x-coordinate               */
  DOM_GLB_Y,                                            /* general bnd cond. y-coordinate               */
#ifdef __THREEDIM__
  DOM_GLB_Z,                                            /* general bnd cond. z-coordinate               */
#endif
  DOM_EVAL_FOR_SD,                              /* evaluate bc for this subdomain		*/
  DOM_N_IN_PARAMS
};

#define DOM_LOC_X       DOM_GLB_X       /* parametrized bnd cond. local x-coord.*/
#ifdef __THREEDIM__
#define DOM_LOC_Y       DOM_GLB_Y       /* parametrized bnd cond. local y-coord.*/
#endif

#define DOM_PARAM_OFFSET        DOM_N_IN_PARAMS

/* subdomain of BC evaluation unknown */
#define DOM_EVAL_SD_UNKNOWN             -1.0

/* boundary types */
#define FIXED         0
#define FREE          1
#define PERIODIC      2
#define NON_PERIODIC  3

/* status for mesh */
#define MESHSTAT_NOTINIT     0
#define MESHSTAT_EMPTY       1
#define MESHSTAT_CNODES      2
#define MESHSTAT_SURFMESH    3
#define MESHSTAT_MESH        4

/* function formats */
typedef INT (*ConfigProcPtr)(INT argc, char **argv);
typedef INT (*CoeffProcPtr)(DOUBLE *, DOUBLE *);
typedef INT (*UserProcPtr)(DOUBLE *, DOUBLE *);

/* macros for BVPDescriptor */
#define BVPD_NAME(d)         ((d)->name)
#define BVPD_MIDPOINT(d)     ((d)->midpoint)
#define BVPD_RADIUS(d)       ((d)->radius)
#define BVPD_CONVEX(d)       ((d)->convex)
#define BVPD_NSUBDOM(d)      ((d)->nSubDomains)
#define BVPD_NPARTS(d)       ((d)->nDomainParts)
#define BVPD_CONFIG(d)       ((d)->ConfigProc)
#define BVPD_S2P_PTR(d)      ((d)->s2p)
#define BVPD_S2P(d,s)        ((d)->s2p[s])
#define BVPD_NCOEFFF(d)      ((d)->numOfCoeffFct)
#define BVPD_NUSERF(d)       ((d)->numOfUserFct)

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

typedef void *BVP;                 /* structure handled by domain module    */
typedef void *BNDS;                /* structure handled by domain module    */
typedef void *BNDP;                /* structure handled by domain module    */

struct BVP_Descriptor
{
  /* general part */
  char name[NAMELEN];                /* name of the BVP                     */

  /* domain part */
  DOUBLE midpoint[DIM];               /* sphere in which the domain lies     */
  DOUBLE radius;
  INT convex;                        /* 1 if domain is convex, 0 if not     */
  INT nSubDomains;                   /* nb. of subdomains,
                                                                                exterior not counted                */
  INT nDomainParts;                                      /* number of parts in the domain		*/
  INT *s2p;                                                      /* pointer to table subbdom --> part	*/

  /* problem part */
  ConfigProcPtr ConfigProc;          /* configuration function              */
  INT numOfCoeffFct;                 /* nb. of coefficient functions        */
  INT numOfUserFct;                  /* nb. of user functions               */
};
typedef struct BVP_Descriptor BVP_DESC;

struct mesh
{
  INT mesh_status;                                       /* status see above					*/
  INT nBndP;                         /* nb. of boundary points              */
  BNDP **theBndPs;                                       /* list of boundary points	            */
  INT nInnP;                         /* nb. of inner nodes                  */
  DOUBLE **Position;                 /* positions of inner nodes            */
  INT nSubDomains;                   /* nb. of subdomains                   */
  INT *nSides;                       /* nb. of boundary sides per subdomain */
  INT **Side_corners;                /* nb. of side corners                 */
  INT **xy_Side;                                         /* triangle_id for prism                */
  INT ***Side_corner_ids;                /* corner ids                          */
  INT *nElements;                    /* nb. of elements per subdomain       */
  INT **Element_corners;             /* nb. of element corners              */
  INT ***Element_corner_ids;         /* nb. of side corners                 */
  INT ***nbElements;                 /* nb. of side corners                 */
  INT **ElemSideOnBnd;               /* used bitwise: sides on bnd for elem */

  /* parallel part */
  char *VertexLevel;                                     /* level of vertex						*/
  /* NULL if all vertex on level 0		*/
  char *VertexPrio;                                      /* priority of vertex					*/
  /* NULL if all vertex are master		*/
  char **ElementLevel;                                   /* level of element in subdomain	*/
  /* NULL if all elements on level 0		*/
  char **ElementPrio;                                    /* priority of element in subdomain	*/
  /* NULL if all elements are master		*/
};
typedef struct mesh MESH;

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* functions for BVP														*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   BVP_GetFirst - Return a pointer to the first BVP

   SYNOPSIS:
   BVP *BVP_GetFirst (void);

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function returns a pointer to the first BVP defined in the domain
   subsystem..

   RETURN VALUE:
   BVP *
   .n   pointer to BVP
   .n   NULL if not found.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BVP_GetNext - Return a pointer to the next BVP

   SYNOPSIS:
   BVP *BVP_GetNext (BVP *theBVP)

   PARAMETERS:
   .  theBVP - BVP structure

   DESCRIPTION:
   This function returns a pointer to the next BVP defined in the domain
   subsystem..

   RETURN VALUE:
   BVP *
   .n   pointer to BVP
   .n   NULL if not found.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BVP_Save - save a BVP

   SYNOPSIS:
   INT BVP_Save (BVP *theBVP, char *name, char argc, char **argv)

   PARAMETERS:
   .  theBVP - BVP structure
   .  name - name of file
   .  argc, argv - command parameters

   DESCRIPTION:
   This function saves a BVP to file named <name>.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BVP_Load - load a BVP

   SYNOPSIS:
   BVP *BVP_Load (char *name, INT argc, char **argv)

   PARAMETERS:
   .  name - name of file
   .  argc, argv - command parameters

   DESCRIPTION:
   This function loads a BVP from file named <name>.

   RETURN VALUE:
   BVP *
   .n   pointer to BVP
   .n   NULL if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BVP_GetByName - get pointer to BVP by name

   SYNOPSIS:
   BVP *BVP_GetByName (char *name)

   PARAMETERS:
   .  name - name of BVP

   DESCRIPTION:
   This function returns the pointer to the BVP specified by its <name>.

   RETURN VALUE:
   BVP *
   .n   pointer to BVP
   .n   NULL if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*
   BVP_Init - initialize a BVP and return a mesh

   SYNOPSIS:
   INT BVP_Init (char *name, HEAP *Heap, MESH *Mesh, INT MarkKey)

   PARAMETERS:
   .  filename - name of file
   .  theHeap - heap
   .  MarkKey - use key for temporary memory allocation (do not Mark/Release)

   DESCRIPTION:
   Function initialize a BVP and returns a mesh.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.

   SEE ALSO:
 */
/****************************************************************************/

/****************************************************************************/
/*D
   BVP_Dispose - dispose a BVP

   SYNOPSIS:
   INT BVP_Dispose (BVP *theBVP)

   PARAMETERS:
   .  theBVP - BVP structure

   DESCRIPTION:
   This function disposes a BVP.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BVP_SetBVPDesc - set BVP-descriptor

   SYNOPSIS:
   INT BVP_SetBVPDesc (BVP *theBVP, BVP_DESC *theBVPDesc)

   PARAMETERS:
   .  theBVP - BVP structure
   .  theBVPDesc - descriptor to set

   DESCRIPTION:
   This function sets the BVP descriptor according to the BVP.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BVP_SetCoeffFct - set coefficient function(s)

   SYNOPSIS:
   INT BVP_SetCoeffFct (BVP *theBVP, INT n, CoeffProcPtr *CoeffFct)

   PARAMETERS:
   .  theBVP - BVP structure
   .  n - nb. of coefficient function or -1 for all

   DESCRIPTION:
   This function one or all coefficient functions.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BVP_SetUserFct - set coefficient function(s)

   SYNOPSIS:
   INT BVP_SetUserFct (BVP *theBVP, INT n, UserProcPtr *UserFct)

   PARAMETERS:
   .  theBVP - BVP structure
   .  n - nb. of user function or -1 for all

   DESCRIPTION:
   This function gives one or all user functions.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BVP_Check - check consistency of BVP

   SYNOPSIS:
   INT BVP_Check (BVP *aBVP)

   PARAMETERS:
   .  aBVP - BVP structure
   .  CheckResult - 0 if ok, 1 if error detected

   DESCRIPTION:
   This function checks consistency of BVP

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* functions called by script commands										*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   BVP_InsertBndP - sets a BNDP from command input

   SYNOPSIS:
   BNDP *BVP_InsertBndP (HEAP *Heap, BVP *theBVP, INT argc, char **argv)

   PARAMETERS:
   .  theBVP - BVP structure
   .  argc, argv - command parameters
   .  theBndP - the BNDP to set

   DESCRIPTION:
   This function sets a BNDP from command input parameters.
   Options are implementation specific.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDP_SaveInsertedBndP - write command to insert this BNDP

   SYNOPSIS:
   INT BNDP_SaveInsertedBndP (BNDP *theBndP, char *data, INT max_data_size)

   PARAMETERS:
   .  theBndP - BNDP structure
   .  data - string to store command
   .  max_data_size - maximal datasize to use

   DESCRIPTION:
   This function writes a command to string which inserts the BNDP.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* functions for BNDP														*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDP_Global - return global coordinates of BNDP

   SYNOPSIS:
   INT BNDP_Global (BNDP *aBndP, DOUBLE *global)

   PARAMETERS:
   .  aBndP - BNDP structure
   .  global - global coordinates

   DESCRIPTION:
   This function returns global coordinates of BNDP

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDP_Move - change global coordinates of free boundary point

   SYNOPSIS:
   INT BNDP_Move (BNDP *aBndP, const DOUBLE global[])

   PARAMETERS:
   .  aBndP - BNDP structure
   .  global - new global coordinates

   DESCRIPTION:
   This function sets global coordinates of a free BNDP. The local coordinates
   will be kept fixed which implies that the topology of the boundary triangulation
   must not be changed!

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDP_BndCond - gets bnd conditions for a BNDP

   SYNOPSIS:
   INT BNDP_BndCond (BNDP *aBndP, INT *n, INT i,
   DOUBLE *in, DOUBLE *value, INT *type)

   PARAMETERS:
   .  aBndP - BNDP structure
   .  i	 - evaluate on patch i
   .  n     - number of BNDS
   .  in    - input vector (if !=NULL has to be allocated with >= DOM_N_IN_PARAMS DOUBLES)
   .  type  - type of bnd cond
   .  value - values

   DESCRIPTION:
   This function gets bnd conditions for a BNDP

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDP_BndPDesc - sets descriptor for BNDP

   SYNOPSIS:
   INT BNDP_BndPDesc (BNDP *theBndP, INT *move, INT *part)

   PARAMETERS:
   .  aBndP - BNDP structure
   .  move  - movable flag (0: no, 1:yes)
   .  part  - domain part

   DESCRIPTION:
   This function sets the descriptor for a BNDP.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDP_BndEDesc - sets descriptor for BNDE (boundary edge)

   SYNOPSIS:
   INT BNDP_BndEDesc (BNDP *aBndP0, BNDP *aBndP1, INT *part)

   PARAMETERS:
   .  aBndP0 - first BNDP
   .  aBndP1 - second BNDP
   .  part	  - domain part ID

   DESCRIPTION:
   This function sets the descriptor for a BNDE.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDP_CreateBndS - creates a BNDS from a nb of BNDPs

   SYNOPSIS:
   BNDS *BNDP_CreateBndS (HEAP *Heap, BNDP **aBndP, INT n)

   PARAMETERS:
   .  Heap  - heap to allocate from
   .  aBndP - ptr to list of BNDP structures
   .  n     - nb of BNDPs

   DESCRIPTION:
   This function creates a BNDS from n BNDPs

   RETURN VALUE:
   BNDS *
   .n   pointer
   .n   NULL if the points describe an inner side
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDP_CreateBndP - sets BNDP from a two of BNDPs

   SYNOPSIS:
   BNDP *BNDP_CreateBndP (HEAP *Heap, BNDP *aBndP0, BNDP *aBndP1, DOUBLE lcoord)

   PARAMETERS:
   .  aBndP0 - first BNDP
   .  aBndP1 - second BNDP
   .  lcoord - local coordinate between P0 and P1 where the BNDP will be created

   DESCRIPTION:
   This function sets a BNDP from two BNDPs

   RETURN VALUE:
   BNDS *
   .n   pointer
   .n   NULL if the points describe an inner point
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDP_Dispose - dispose a BNDP

   SYNOPSIS:
   INT BNDP_Dispose (HEAP *Heap, BNDP *aBndP)

   PARAMETERS:
   .  Heap - heap
   .  aBndP - BNDP

   DESCRIPTION:
   This function disposes a BNDP

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDP_SaveBndP - save a BNDP

   SYNOPSIS:
   INT BNDP_SaveBndP (BNDP *theBndP, FILE *stream)

   PARAMETERS:
   .  theBndP - BNDP
   .  stream - file

   DESCRIPTION:
   This function saves a BNDP on a file.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BVP_LoadBndP - load a BNDP

   SYNOPSIS:
   BNDP *BNDP_LoadBndP (BVP *theBVP, HEAP *Heap, FILE *stream)

   PARAMETERS:
   .  theBVP - BVP structure
   .  Heap   - heap to allocate from
   .  stream - file

   DESCRIPTION:
   This function loads a BNDP with the format given by BVP_SaveBndP and
   allocates it from the given heap.

   RETURN VALUE:
   INT
   .n      0 if ok
   .n      1 if error.
   D*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* functions for BNDS														*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDS_Global - gets global coordinates of local position on BNDS

   SYNOPSIS:
   INT BNDS_Local2Global (BNDS *aBndS, DOUBLE *local, DOUBLE *global)

   PARAMETERS:
   .  aBndS  - BNDS structure
   .  local  - local coordinate on BNDS
   .  global - global coordinate

   DESCRIPTION:
   This function gets global coordinates of local position on BNDS

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDS_BndCond - gets bnd conditions of local position on BNDS

   SYNOPSIS:
   INT BNDS_BndCond (BNDS *aBndS, DOUBLE *local, DOUBLE *in,
   INT *type, DOUBLE *value)

   PARAMETERS:
   .  aBndS - BNDS structure
   .  local - local coordinate on BNDS
   .  in    - input vector (if !=NULL has to be allocated with >= DOM_N_IN_PARAMS DOUBLES)
   .  type  - type of bnd cond
   .  value - values

   DESCRIPTION:
   This function gets bnd conditions of local position on BNDS.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDS_BndSDesc - sets descriptor for BNDS

   SYNOPSIS:
   INT BNDS_BndSDesc (BNDS *aBndS, INT *id, INT *nbid, INT *part)

   PARAMETERS:
   .  aBndS - BNDS structure
   .  id  - subdomain ID of the element with aBndS
   .  nbid  - subdomain ID of the neighbour element (across the boundary)
   .  right - ID of right subdomain
   .  part  - domain part

   DESCRIPTION:
   This function sets the descriptor for a BNDS.

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDS_CreateBndP - create BNDP on BNDS

   SYNOPSIS:
   BNDP *BNDS_CreateBndP (HEAP *Heap, BNDS *aBndS, DOUBLE *local)

   PARAMETERS:
   .  Heap  - heap to allocate from
   .  aBndS - BNDS structure
   .  local - local coordinate on BNDS
   .  size  - size used for aBndP

   DESCRIPTION:
   This function creates a boundary point (BNDP) on a BNDS

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/

/****************************************************************************/
/*D
   BNDS_Dispose - dispose BNDS

   SYNOPSIS:
   INT BNDS_Dispose (HEAP *Heap, BNDS *theBndS)

   PARAMETERS:
   .  Heap - heap
   .  theBndS - BNDS struct

   DESCRIPTION:
   This function disposes a BNDS

   RETURN VALUE:
   INT
   .n   0 if ok
   .n   1 if error.
   D*/
/****************************************************************************/
BVP        *BVP_GetNext           (BVP *theBVP);
INT         BVP_Save              (BVP *theBVP, char *name, char *mgname, HEAP *theHeap, INT argc, char **argv);
BVP        *BVP_Load              (char *name, INT argc, char **argv);
BVP        *BVP_GetByName         (char *name);
BVP        *BVP_Init              (char *name, HEAP *Heap, MESH *Mesh, INT MarkKey);
INT         BVP_Dispose           (BVP *theBVP);
INT         BVP_SetBVPDesc        (BVP *theBVP, BVP_DESC *theBVPDesc);
INT         BVP_SetCoeffFct       (BVP *theBVP, INT n, CoeffProcPtr *CoeffFct);
INT         BVP_SetUserFct        (BVP *theBVP, INT n, UserProcPtr *UserFct);
INT             BVP_Check                         (BVP *aBVP);
BNDP*           BVP_InsertBndP            (HEAP *Heap, BVP *theBVP, INT argc, char **argv);
INT         BNDP_SaveInsertedBndP (BNDP *theBndP, char *data, INT max_data_size);
MESH       *BVP_GenerateMesh      (HEAP *Heap, BVP *aBVP, INT argc, char **argv, INT MarkKey);
INT         BNDP_Global           (BNDP *theBndP, DOUBLE *global);
INT                     BNDP_Move                         (BNDP *aBndP, const DOUBLE global[]);
INT         BNDP_BndCond          (BNDP *theBndP, INT *n, INT i, DOUBLE *in, DOUBLE *value, INT *type);
INT         BNDP_BndPDesc         (BNDP *theBndP, INT *move, INT *part);
INT         BNDP_BndEDesc         (BNDP *theBndP0, BNDP *theBndP1, INT *part);
BNDS*       BNDP_CreateBndS       (HEAP *Heap, BNDP **theBndP, INT n);
BNDP*       BNDP_CreateBndP       (HEAP *Heap, BNDP *theBndP0, BNDP *theBndP1, DOUBLE lcoord);
INT         BNDP_Dispose          (HEAP *Heap, BNDP *theBndP);
INT         BNDP_SaveBndP         (BNDP *theBndP);
BNDP       *BNDP_LoadBndP         (BVP *theBVP, HEAP *Heap);
INT         BNDS_Global           (BNDS *theBndS, DOUBLE *local, DOUBLE *global);
INT         BNDS_BndCond          (BNDS *theBndS, DOUBLE *local, DOUBLE *in, DOUBLE *value, INT *type);
INT         BNDS_BndSDesc         (BNDS *theBndS, INT *id, INT *nbid, INT *part);
BNDP*       BNDS_CreateBndP       (HEAP *Heap, BNDS *theBndS, DOUBLE *local);
BVP        *BVP_GetFirst          (void);
INT         BNDS_Dispose          (HEAP *Heap, BNDS *theBndS);


#ifdef ModelP

void DomInitParallel     (INT TypeBndP, INT TypeBndS);
void DomHandlerInit      (INT handlerSet);

void BElementXferBndS    (BNDS **bnds, int n, int proc, int prio);
void BElementGatherBndS  (BNDS **bnds, int n, int cnt, char *data);
void BElementScatterBndS (BNDS **bnds, int n, int cnt, char *data);

void BVertexXferBndP     (BNDP *bndp, int proc, int prio);
void BVertexGatherBndP   (BNDP *bndp, int cnt, char *data);
void BVertexScatterBndP  (BNDP **bndp, int cnt, char *data);

#endif

/* miscellaneous */
INT         InitDom               (void);

#endif
