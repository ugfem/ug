// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/** \defgroup dom Domains
 * \ingroup ug
 */
/*! \file domain.h
 * \ingroup dom
 */

/** \addtogroup dom
 *
 * @{
 */


/****************************************************************************/
/*                                                                          */
/* File:      domain.h                                                      */
/*                                                                          */
/* Purpose:   standard header file template                                 */
/*                                                                          */
/* Author:    klaus Johannsen                                               */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70550 Stuttgart                                               */
/*            email: ug@ica3.uni-stuttgart.de                               */
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

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*          compile time constants defining static data size (i.e. arrays)  */
/*          other constants                                                 */
/*          macros                                                          */
/*                                                                          */
/****************************************************************************/

#include "compiler.h"
#include "ugenv.h"
#include "heaps.h"

#include "namespace.h"

START_UGDIM_NAMESPACE

/** \brief Offset of additional parameters for boundary condition function call */
enum DOM_IN_PARAMS {

  DOM_GLB_X,              /*!< General bnd cond. x-coordinate               */
  DOM_GLB_Y,              /*!< General bnd cond. y-coordinate               */
  DOM_GLB_Z,              /*!< General bnd cond. z-coordinate               */
  DOM_EVAL_FOR_SD,        /*!< Evaluate bc for this subdomain               */
  DOM_N_IN_PARAMS
};

#define DOM_LOC_X       DOM_GLB_X       /*!< Parametrized bnd cond. local x-coord.*/
#define DOM_LOC_Y       DOM_GLB_Y       /*!< parametrized bnd cond. local y-coord.*/

#define DOM_PARAM_OFFSET        DOM_N_IN_PARAMS

/** \brief Subdomain of BC evaluation unknown */
#define DOM_EVAL_SD_UNKNOWN             -1.0

/** \brief Boundary types */
enum BoundaryType {FIXED, FREE, PERIODIC, NON_PERIODIC};

/** \brief Status for mesh */
enum MeshStatus {MESHSTAT_NOTINIT,
                 MESHSTAT_EMPTY,
                 MESHSTAT_CNODES,
                 MESHSTAT_SURFMESH,
                 MESHSTAT_MESH};

/** @name Function formats */
/*@{*/
typedef INT (*ConfigProcPtr)(INT argc, char **argv);
typedef INT (*CoeffProcPtr)(DOUBLE *, DOUBLE *);
typedef INT (*UserProcPtr)(DOUBLE *, DOUBLE *);
/*@}*/

/** @name Macros for BVPDescriptor */
/*@{*/
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
/*@}*/

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

typedef void *BVP;                 /*!< Structure handled by domain module    */
typedef void *BNDS;                /*!< Structure handled by domain module    */
typedef void *BNDP;                /*!< Structure handled by domain module    */

/** \todo Please doc me! */
struct BVP_Descriptor
{
  /** @name General part */
  /*@{*/
  /** \brief Name of the BVP */
  char name[NS_PREFIX NAMELEN];
  /*@}*/

  /** @name Domain part */
  /*@{*/
  /** \brief Midpoint of a sphere in which the domain lies     */
  DOUBLE midpoint[3];

  /** \brief Radius of a sphere in which the domain lies     */
  DOUBLE radius;

  /** \brief 1 if domain is convex, 0 if not */
  INT convex;

  /** \brief Number of subdomains, exterior not counted                */
  INT nSubDomains;

  /** \brief Number of parts in the domain               */
  INT nDomainParts;

  /** \brief Pointer to table subbdom --> part   */
  INT *s2p;
  /*@}*/

  /** @name Problem part */
  /*@{*/
  /** \brief Configuration function              */
  ConfigProcPtr ConfigProc;

  /** \brief Number of coefficient functions        */
  INT numOfCoeffFct;

  /** \brief Number of user functions               */
  INT numOfUserFct;
  /*@}*/
};
typedef struct BVP_Descriptor BVP_DESC;


/** \todo Please doc me! */
struct mesh
{

  /** \brief Status see above                                    */
  INT mesh_status;

  /** \brief Number of boundary points              */
  INT nBndP;

  /** \brief List of boundary points                 */
  BNDP **theBndPs;

  /** \brief Number of inner nodes                  */
  INT nInnP;

  /** \brief Positions of inner nodes            */
  DOUBLE **Position;

  /** \brief Number of subdomains                   */
  INT nSubDomains;

  /** \brief Number of boundary sides per subdomain */
  INT *nSides;

  /** \brief Number of side corners                 */
  INT **Side_corners;

  /** \brief Triangle_id for prism                */
  INT **xy_Side;

  /** \brief Corner ids                          */
  INT ***Side_corner_ids;

  /** \brief Number of elements per subdomain       */
  INT *nElements;

  /** \brief Number of element corners              */
  INT **Element_corners;

  /** \brief Number of side corners                 */
  INT ***Element_corner_ids;

  /** \brief Number of side corners                 */
  INT ***nbElements;

  /** \brief Used bitwise: sides on bnd for elem */
  INT **ElemSideOnBnd;

  /** @name Parallel part */
  /*@}*/
  /** \brief Level of vertex NULL if all vertex on level 0 */
  char *VertexLevel;

  /** \brief Priority of vertex   NULL if all vertex are master               */
  char *VertexPrio;

  /** \brief Level of element in subdomain   NULL if all elements on level 0             */
  char **ElementLevel;

  /** \brief Priority of element in subdomain  NULL if all elements are master             */
  char **ElementPrio;
  /*@}*/
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
/*                                                                          */
/* functions for BVP                                                        */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* functions called by script commands                                      */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* functions for BNDP                                                       */
/*                                                                          */
/****************************************************************************/









/****************************************************************************/
/*                                                                          */
/* functions for BNDS                                                       */
/*                                                                          */
/****************************************************************************/






/****************************************************************************/
/** \brief Return a pointer to the next BVP

 * @param theBVP - BVP structure

   This function returns a pointer to the next BVP defined in the domain
   subsystem.

 * @return <ul>
 *   <li> pointer to BVP </li>
 *   <li> NULL if not found. </li>
 * </ul>
 */
/****************************************************************************/
BVP        *BVP_GetNext           (BVP *theBVP);

/****************************************************************************/
/** \brief Save a BVP
 *
 * @param theBVP - BVP structure
 * @param name - name of file
 * @param argc, argv - command parameters

   This function saves a BVP to file named \<name\>.

 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error. </li>
 * </ul>
 */
/****************************************************************************/
INT BVP_Save (BVP *theBVP, char *name, char *mgname, NS_PREFIX HEAP *theHeap, INT argc, char **argv);

/****************************************************************************/
/** \brief Load a BVP
 *
 * @param name - name of file
 * @param argc, argv - command parameters

   This function loads a BVP from file named \<name\>.

 * @return <ul>
 *   <li> pointer to BVP </li>
 *   <li> NULL if error. </li>
 * </ul>
 */
/****************************************************************************/
BVP        *BVP_Load              (char *name, INT argc, char **argv);

/****************************************************************************/
/** \brief Get pointer to BVP by name
 *
 * @param name - name of BVP

   This function returns the pointer to the BVP specified by its \<name\>.

 * @return <ul>
 *   <li> pointer to BVP </li>
 *   <li> NULL if error. </li>
 * </ul>
 */
/****************************************************************************/
BVP        *BVP_GetByName         (const char *name);

void Set_Current_BVP(BVP* theBVP);

/****************************************************************************/
/** \brief Initialize a BVP and return a mesh
 *
 * @param filename - name of file
 * @param theHeap - heap
 * @param MarkKey - use key for temporary memory allocation (do not Mark/Release)

   Function initialize a BVP and returns a mesh.

 * @return <ul>
 *   <li>    0 if ok  </li>
 *   <li>    1 if error.       </li>
 * </ul>
 */
/****************************************************************************/
BVP *BVP_Init (char *filename, NS_PREFIX HEAP *theHeap, MESH *Mesh, INT MarkKey);

/****************************************************************************/
/** \brief Dispose a BVP
 *
 * @param theBVP - BVP structure

   This function disposes a BVP.

 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error. </li>
 * </ul>
 */
/****************************************************************************/
INT         BVP_Dispose           (BVP *theBVP);

/****************************************************************************/
/** \brief Set BVP-descriptor
 *
 * @param theBVP - BVP structure
 * @param theBVPDesc - descriptor to set

   This function sets the BVP descriptor according to the BVP.

 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error. </li>
 * </ul>
 */
/****************************************************************************/
INT         BVP_SetBVPDesc        (BVP *theBVP, BVP_DESC *theBVPDesc);


/****************************************************************************/
/** \brief Set coefficient function(s)
 *
 * @param theBVP - BVP structure
 * @param n - Number of coefficient function or -1 for all

   This function one or all coefficient functions.

 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error. </li>
 * </ul>
 */
/****************************************************************************/
INT         BVP_SetCoeffFct       (BVP *theBVP, INT n, CoeffProcPtr *CoeffFct);

/****************************************************************************/
/** \brief Set coefficient function(s)
 *
 * @param theBVP - BVP structure
 * @param n - Number of user function or -1 for all

   This function gives one or all user functions.

 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error. </li>
 * </ul>
 */
/****************************************************************************/
INT         BVP_SetUserFct        (BVP *theBVP, INT n, UserProcPtr *UserFct);

/****************************************************************************/
/** \brief Check consistency of BVP
 *
 * @param aBVP - BVP structure

   This function checks consistency of BVP

 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error. </li>
 * </ul>
 */
/****************************************************************************/
INT             BVP_Check                         (BVP *aBVP);

/****************************************************************************/
/** \brief Sets a BNDP from command input
 *
 * @param theBVP - BVP structure
 * @param argc, argv - command parameters

   This function sets a BNDP from command input parameters.
   Options are implementation specific.

 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error. </li>
 * </ul>
 */
/****************************************************************************/
BNDP* BVP_InsertBndP (NS_PREFIX HEAP *Heap, BVP *theBVP, INT argc, char **argv);

/****************************************************************************/
/** \brief Write command to insert this BNDP
 *
 * @param theBndP - BNDP structure
 * @param data - string to store command
 * @param max_data_size - maximal datasize to use

   This function writes a command to string which inserts the BNDP.

 * @return <ul>
 *   <li> 0 if ok
 *   <li> 1 if error.
 * </ul> */
/****************************************************************************/
INT         BNDP_SaveInsertedBndP (BNDP *theBndP, char *data, INT max_data_size);

MESH       *BVP_GenerateMesh (NS_PREFIX HEAP *Heap, BVP *aBVP, INT argc, char **argv, INT MarkKey);

/****************************************************************************/
/** \brief Return global coordinates of BNDP
 *
 * @param theBndP - BNDP structure
 * @param global - global coordinates

   This function returns global coordinates of BNDP

 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error. </li>
 * </ul> */
/****************************************************************************/
INT         BNDP_Global           (BNDP *theBndP, DOUBLE *global);

/****************************************************************************/
/** \brief Change global coordinates of free boundary point
 *
 * @param aBndP - BNDP structure
 * @param global - new global coordinates

   This function sets global coordinates of a free BNDP. The local coordinates
   will be kept fixed which implies that the topology of the boundary triangulation
   must not be changed!

 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error. </li>
 * </ul> */
/****************************************************************************/
INT         BNDP_Move                         (BNDP *aBndP, const DOUBLE global[]);

/****************************************************************************/
/** \brief Gets boundary conditions for a BNDP
 *
 * @param theBndP - BNDP structure
 * @param i     - evaluate on patch i
 * @param n     - number of BNDS
 * @param in    - input vector (if !=NULL has to be allocated with >= DOM_N_IN_PARAMS DOUBLES)
 * @param type  - type of bnd cond
 * @param value - values

   This function gets bnd conditions for a BNDP

 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error. </li>
 * </ul> */
/****************************************************************************/
INT         BNDP_BndCond          (BNDP *theBndP, INT *n, INT i, DOUBLE *in, DOUBLE *value, INT *type);

/****************************************************************************/
/** \brief Sets descriptor for BNDP
 *
 * @param theBndP - BNDP structure
 * @param move  - movable flag (0: no, 1:yes)
 * @param part  - domain part

   This function sets the descriptor for a BNDP.

 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error. </li>
 * </ul> */
/****************************************************************************/
INT         BNDP_BndPDesc         (BNDP *theBndP, INT *move, INT *part);

/****************************************************************************/
/** \brief Sets descriptor for BNDE (boundary edge)
 *
 * @param theBndP0 - first BNDP
 * @param theBndP1 - second BNDP
 * @param part   - domain part ID

   This function sets the descriptor for a BNDE.

 * @return <ul>
 *   <li> 0 if ok
 *   <li> 1 if error.
 * </ul> */
/****************************************************************************/
INT         BNDP_BndEDesc         (BNDP *theBndP0, BNDP *theBndP1, INT *part);

/****************************************************************************/
/** \brief Creates a BNDS from a nb of BNDPs
 *
 * @param Heap  - heap to allocate from
 * @param theBndP - Pointer to list of BNDP structures
 * @param n     - Number of BNDPs

   This function creates a BNDS from n BNDPs

 * @return <ul>
 *   <li> pointer </li>
 *   <li> NULL if the points describe an inner side </li>
 * </ul> */
/****************************************************************************/
BNDS* BNDP_CreateBndS (NS_PREFIX HEAP *Heap, BNDP **theBndP, INT n);

/****************************************************************************/
/** \brief Sets BNDP from a two of BNDPs
 *
 * @param theBndP0 - first BNDP
 * @param theBndP1 - second BNDP
 * @param lcoord - local coordinate between P0 and P1 where the BNDP will be created

   This function sets a BNDP from two BNDPs

 * @return <ul>
 *   <li> pointer </li>
 *   <li> NULL if the points describe an inner point </li>
 * </ul> */
/****************************************************************************/
BNDP*       BNDP_CreateBndP (NS_PREFIX HEAP *Heap, BNDP *theBndP0, BNDP *theBndP1, DOUBLE lcoord);

/****************************************************************************/
/** \brief Dispose a BNDP
 *
 * @param Heap - heap
 * @param theBndP - BNDP

   This function disposes a BNDP

 * @return <ul>
 *   <li>    0 if ok  </li>
 *   <li>    1 if error.    </li>
 * </ul> */
/****************************************************************************/
INT         BNDP_Dispose          (NS_PREFIX HEAP *Heap, BNDP *theBndP);

/****************************************************************************/
/** \brief Save a BNDP
 *
 * @param theBndP - BNDP
 * @param stream - file

   This function saves a BNDP on a file.

 * @return <ul>
 *   <li>    0 if ok   </li>
 *   <li>    1 if error.     </li>
 * </ul> */
/****************************************************************************/
INT         BNDP_SaveBndP         (BNDP *theBndP);
INT         BNDP_SaveBndP_Ext     (BNDP *theBndP);

/****************************************************************************/
/** \brief Load a BNDP
 *
 * @param theBVP - BVP structure
 * @param Heap   - heap to allocate from
 * @param stream - file

   This function loads a BNDP with the format given by BVP_SaveBndP and
   allocates it from the given heap.

 * @return <ul>
 *   <li>    0 if ok  </li>
 *   <li>    1 if error.      </li>
 * </ul> */
/****************************************************************************/
BNDP       *BNDP_LoadBndP         (BVP *theBVP, NS_PREFIX HEAP *Heap);
BNDP       *BNDP_LoadBndP_Ext     (void);

/****************************************************************************/
/** \brief Gets global coordinates of local position on BNDS
 *
 * @param theBndS  - BNDS structure
 * @param local  - local coordinate on BNDS
 * @param global - global coordinate

   This function gets global coordinates of local position on BNDS

 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error. </li>
 * </ul> */
/****************************************************************************/
INT         BNDS_Global           (BNDS *theBndS, DOUBLE *local, DOUBLE *global);

/****************************************************************************/
/** \brief Gets boundary conditions of local position on BNDS
 *
 * @param theBndS - BNDS structure
 * @param local - local coordinate on BNDS
 * @param in    - input vector (if !=NULL has to be allocated with >= DOM_N_IN_PARAMS DOUBLES)
 * @param type  - type of bnd cond
 * @param value - values

   This function gets bnd conditions of local position on BNDS.

 * @return <ul>
 *   <li> 0 if ok
 *   <li> 1 if error.
 * </ul> */
/****************************************************************************/
INT         BNDS_BndCond          (BNDS *theBndS, DOUBLE *local, DOUBLE *in, DOUBLE *value, INT *type);

/****************************************************************************/
/** \brief Sets descriptor for BNDS
 *
 * @param theBndS - BNDS structure
 * @param id  - subdomain ID of the element with aBndS
 * @param nbid  - subdomain ID of the neighbour element (across the boundary)
 * @param part  - domain part

   This function sets the descriptor for a BNDS.

 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error. </li>
 * </ul>
 */
/****************************************************************************/
INT         BNDS_BndSDesc         (BNDS *theBndS, INT *id, INT *nbid, INT *part);

/****************************************************************************/
/** \brief Create BNDP on BNDS
 *
 * @param Heap  - heap to allocate from
 * @param theBndS - BNDS structure
 * @param local - local coordinate on BNDS

   This function creates a boundary point (BNDP) on a BNDS

 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error. </li>
 * </ul> */
/****************************************************************************/
BNDP*       BNDS_CreateBndP       (NS_PREFIX HEAP *Heap, BNDS *theBndS, DOUBLE *local);

/****************************************************************************/
/** \brief Return a pointer to the first BVP
 *
   This function returns a pointer to the first BVP defined in the domain
   subsystem..

 * @return <ul>
 *   <li> pointer to BVP </li>
 *   <li> NULL if not found. </li>
 * </ul> */
/****************************************************************************/
BVP        *BVP_GetFirst          (void);

/****************************************************************************/
/** \brief Dispose BNDS
 *
 * @param Heap - heap
 * @param theBndS - BNDS struct

   This function disposes a BNDS

 * @return <ul>
 *   <li> 0 if ok </li>
 *   <li> 1 if error. </li>
 * </ul> */
/****************************************************************************/
INT         BNDS_Dispose          (NS_PREFIX HEAP *Heap, BNDS *theBndS);

/****************************************************************************/
/** \brief Gets surface ids for a BNDP
 *
 * @param aBndP - BNDP structure
 * @param i     - evaluate on patch i
 * @param n     - number of BNDS

   This function returns surface ids of the n surfaces on which
   the BNDP resides

 * @return  surface id
 */
/****************************************************************************/
INT         BNDP_SurfaceId        (BNDP *aBndP, INT *n, INT i);


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


END_NAMESPACE

#endif

/** @} */
