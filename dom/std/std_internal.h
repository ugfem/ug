// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef UG_STD_INTERNAL_H
#define UG_STD_INTERNAL_H

#if defined(AUTOTOOLS_BUILD) && ! defined(UGLIB)
#error internal header! Must not be used in Applications!
#endif

#include "domain.h"

#include "std_domain.h"

#include "namespace.h"
#include "dimension.h"

START_UGDIM_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*    compile time constants defining static data size (i.e. arrays)        */
/*    other constants                                                       */
/*    macros                                                                */
/*                                                                          */
/****************************************************************************/

typedef DOUBLE COORD_BND_VECTOR[DIM_OF_BND];

enum {BVP_STANDARD,
      BVP_MARC};

enum PatchType {POINT_PATCH_TYPE,
                LINE_PATCH_TYPE,
                LINEAR_PATCH_TYPE,
                PARAMETRIC_PATCH_TYPE,
                MARC_0_PATCH_TYPE,
                MARC_1_PATCH_TYPE,
                MARC_2_PATCH_TYPE};

/** @name  Macros for DOMAIN_PART_INFO */
/*@{*/
#define DPI_SD2P_PTR(p)                                 ((p)->sd2part)
#define DPI_SD2P(p,sd)                                  ((p)->sd2part[sd])
#define DPI_SG2P_PTR(p)                                 ((p)->sg2part)
#define DPI_SG2P(p,sg)                                  ((p)->sg2part[sg])
#define DPI_LN2P_PTR(p)                                 ((p)->ln2part)
#define DPI_LN2P(p,c0,c1)                               ((p)->ln2part[c0][c1])
#define DPI_PT2P_PTR(p)                                 ((p)->pt2part)
#define DPI_PT2P(p,pt)                                  ((p)->pt2part[pt])
/*@}*/

/** @name Macros for DOMAIN */
/*@{*/
#define DOMAIN_MIDPOINT(p)                              ((p)->MidPoint)
#define DOMAIN_RADIUS(p)                                ((p)->radius)
#define DOMAIN_CONVEX(p)                                ((p)->domConvex)
#define DOMAIN_NSEGMENT(p)                              ((p)->numOfSegments)
#define DOMAIN_NCORNER(p)                               ((p)->numOfCorners)
#define DOMAIN_NPARTS(p)                                ((p)->nParts)
#define DOMAIN_PARTINFO(p)                              ((p)->dpi)
/*@}*/

/** @name Macros for STD_BVP */
/*@{*/
#define MAX_CORNERS_OF_LINEAR_PATCH      DIM


#define STD_BVP_NAME(p)                                 ENVITEM_NAME(p)

#define STD_BVP_DOMAIN(p)                               ((p)->Domain)
#define STD_BVP_PROBLEM(p)                              ((p)->Problem)

#define STD_BVP_MIDPOINT(p)                             ((p)->MidPoint)
#define STD_BVP_RADIUS(p)                               ((p)->radius)
#define STD_BVP_CONVEX(p)                               ((p)->domConvex)

#define STD_BVP_NCORNER(p)                              ((p)->ncorners)
#define STD_BVP_NSIDES(p)                               ((p)->nsides)
#define STD_BVP_SIDEOFFSET(p)                   ((p)->sideoffset)
#define STD_BVP_PATCHES(p)                              ((p)->patches)
#define STD_BVP_PATCH(p,i)                              ((p)->patches[i])
#define STD_BVP_CORNERPATCH(p,i)                ((p)->patches[i])
#define STD_BVP_LINEPATCH(p,i)                  ((p)->patches[i+STD_BVP_NCORNER(p)])
#define STD_BVP_SIDEPATCH(p,i)                  ((p)->patches[i+STD_BVP_SIDEOFFSET(p)])

#define STD_BVP_NSUBDOM(p)                              ((p)->numOfSubdomains)
#define STD_BVP_NDOMPART(p)                             ((p)->nDomainParts)
#define STD_BVP_S2P_PTR(p)                              ((p)->s2p)
#define STD_BVP_S2P(p,s)                                ((p)->s2p[s])
#define STD_BVP_NPATCH(p)                               ((p)->nPatch)
#define STD_BVP_NCOEFFPROC(p)                   ((p)->numOfCoeffFct)
#define STD_BVP_NUSERPROC(p)                    ((p)->numOfUserFct)
#define STD_BVP_CONFIGPROC(p)                   ((p)->ConfigProc)
#define STD_BVP_COEFFPROC(p,i)                  ((CoeffProcPtr)((p)->CU_ProcPtr[i]))
#define STD_BVP_USERPROC(p,i)                   ((UserProcPtr)((p)->CU_ProcPtr[i+STD_BVP_NCOEFFPROC(p)]))

#define GetSTD_BVP(p)                           ((STD_BVP *)(p))
/*@}*/

/****************************************************************************/
/*                                                                          */
/* macros for boundary segments                                             */
/*                                                                          */
/****************************************************************************/

/** @name Macros for boundary segments */
/*@{*/
#define SEG_LEFT(p)                                             ((p)->left)
#define SEG_RIGHT(p)                                    ((p)->right)
#define SEG_ID(p)                                               ((p)->id)
#define SEG_TYPE(p)                                             ((p)->segType)
#define SEG_POINTS(p)                                   ((p)->points)
#define SEG_POINT(p,i)                                  ((p)->points[i])
#define SEG_RESOLUTION(p)                               ((p)->resolution)
#define SEG_ALPHA(p)                                    ((p)->alpha)
#define SEG_FROM(p,i)                                   ((p)->alpha[i])
#define SEG_BETA(p)                                             ((p)->beta)
#define SEG_TO(p,i)                                             ((p)->beta[i])
#define SEG_FUNC(p)                                             ((p)->BndSegFunc)
#define SEG_DATA(p)                                             ((p)->data)
/*@}*/

/****************************************************************************/
/*                                                                          */
/* macros for patches                                                       */
/*                                                                          */
/****************************************************************************/

/** @name Macros for patches */
/*@{*/
#define PATCH_FIXED                             0
#define PATCH_BND_OF_FREE               1
#define PATCH_FREE                              2

#define PATCH_TYPE(p)           (p)->ge.type
#define PATCH_STATE(p)          (p)->ge.state
#define PATCH_IS_FREE(p)                ((p)->ge.state==PATCH_FREE)
#define PATCH_IS_FIXED(p)               ((p)->ge.state==PATCH_FIXED)
#define PATCH_ID(p)             (p)->ge.id
#define POINT_PATCH_N(p)        (p)->po.npatches
#define POINT_PATCH_PID(p,i)    (p)->po.pop[i].patch_id
#define POINT_PATCH_CID(p,i)    (p)->po.pop[i].corner_id
#define LINE_PATCH_C0(p)                ((p)->li.c0)
#define LINE_PATCH_C1(p)                ((p)->li.c1)
#define LINE_PATCH_N(p)         (p)->li.npatches
#define LINE_PATCH_PID(p,i)     (p)->li.lop[i].patch_id
#define LINE_PATCH_CID0(p,i)    (p)->li.lop[i].corner_id[0]
#define LINE_PATCH_CID1(p,i)    (p)->li.lop[i].corner_id[1]
#define PARAM_PATCH_LEFT(p)     (p)->pa.left
#define PARAM_PATCH_RIGHT(p)    (p)->pa.right
#define PARAM_PATCH_POINTS(p,i) (p)->pa.points[i]
#define PARAM_PATCH_RES(p)              ((p)->pa.resolution)
#define PARAM_PATCH_RANGE(p)    (p)->pa.range
#define PARAM_PATCH_BS(p)       (p)->pa.BndSegFunc
#define PARAM_PATCH_BSD(p)      (p)->pa.bs_data
#define PARAM_PATCH_BC(p)       (p)->pa.BndCond
#define PARAM_PATCH_BCD(p)      (p)->pa.bc_data
#define LINEAR_PATCH_LEFT(p)    (p)->lp.left
#define LINEAR_PATCH_RIGHT(p)   (p)->lp.right
#define LINEAR_PATCH_N(p)       (p)->lp.corners
#define LINEAR_PATCH_POINTS(p,i) (p)->lp.points[i]
#define LINEAR_PATCH_POS(p,i)    (p)->lp.pos[i]

#define BND_PATCH_ID(p)         (((BND_PS *)p)->patch_id)
#define BND_DATA(p)             (((BND_PS *)p)->data)
#define BND_N(p)                (((BND_PS *)p)->n)
#define BND_LOCAL(p,i)          (((BND_PS *)p)->local[i])
#define BND_SIZE(p)             ((((BND_PS *)p)->n-1)*sizeof(COORD_BND_VECTOR)+sizeof(BND_PS))
#define M_BNDS_NSIZE(n)         (((n)-1)*sizeof(M_BNDP)+sizeof(M_BNDS))
#define M_BNDS_SIZE(p)          M_BNDS_NSIZE(((M_BNDS *)(p))->n)

#define IF_MARC(p) \
  if (PATCH_TYPE(currBVP->patches[BND_PATCH_ID(p)]) >= MARC_0_PATCH_TYPE)
/*@}*/

/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* domain definition data structures                                        */
/*                                                                          */
/****************************************************************************/


/*----------- definition of structs ----------------------------------------*/


/** \brief Data type describing a domain. */
struct domain {

  /** \brief Fields for environment directory */
  NS_PREFIX ENVDIR d;

  /** \brief A point in the middle of the domain */
  DOUBLE MidPoint[DIM];

  /** \brief Defines sphere around MidPoint containing the domain */
  DOUBLE radius;

  /** \brief Number of boundary segments */
  INT numOfSegments;

  /** \brief Number of corner points */
  INT numOfCorners;

  /** \brief Is the domain convex? */
  INT domConvex;

  /** @name Description of domain parts */
  /*@{*/

  /** \brief Number of parts in the domain */
  INT nParts;

  /** \brief Domain part info */
  const DOMAIN_PART_INFO *dpi;
  /*@}*/
};

/** \brief Data structure defining part of the boundary of a domain */
struct boundary_segment {

  /** \brief Field for environment directory */
  NS_PREFIX ENVVAR v;

  /** @name Fields for boundary segment */
  /*@{*/
  /** \brief Number of left and right subdomain */
  INT left,right;

  /** \brief Unique id of that segment */
  INT id;

  /** \brief Segment type, see above
   *
   * \todo See where???*/
  INT segType;

  /** \brief Numbers of the vertices (ID) */
  INT points[CORNERS_OF_BND_SEG];

  /** \brief Measure for the curvature */
  INT resolution;

  /** \brief Parameter interval used*/
  DOUBLE alpha[DIM_OF_BND],beta[DIM_OF_BND];

  /** \brief Pointer to definition function */
  BndSegFuncPtr BndSegFunc;

  /** \brief Can be used by application to find data */
  void *data;
  /*@}*/
};

/** \brief ???
 *
 * \todo Please doc me!
 */
struct linear_segment {

  /** \brief Field for environment directory */
  NS_PREFIX ENVVAR v;

  /* fields for boundary segment */
  /** \brief  Number of left and right subdomain */
  INT left,right;

  /** \brief  Unique id of that segment                  */
  INT id;

  /** \brief  Number of corners                  */
  INT n;

  /** \brief  Numbers of the vertices (ID)*/
  INT points[MAX_CORNERS_OF_LINEAR_PATCH];

  /** \brief  Coordinates  */
  DOUBLE x[MAX_CORNERS_OF_LINEAR_PATCH][DIM_OF_BND];
};

/****************************************************************************/
/*                                                                          */
/* problem data structure                                                   */
/*                                                                          */
/****************************************************************************/

/*----------- definition of structs ----------------------------------------*/

/** \brief Data type describing a problem. */
struct problem {

  /** \brief Field for environment directory
   *
   * The problem is an environment directory. This directory is a subdirectory
   * of the domain where this problem corresponds to. d also contains the
   * name of the problem.
   */
  NS_PREFIX ENVDIR d;

  /* fields for problem */
  /** \brief Used to identify problem type
   *
   * Problem class identification number. This number is used to determine
   * that the problem description coincides with the pde solved by the
   * problem class library.
   */
  INT problemID;

  /** \brief Procedure to reinitialize problem
   *
   * Pointer to a user definable function that is executed when the reinit
   * command is given in the UG shell.
   */
  ConfigProcPtr ConfigProblem;

  /** \brief Number of coefficient functions
   *
   *  User definable coefficient functions come in two flavours.
   * They are either of type CoeffProcPtr or of type UserProcPtr.
   * numOfCoeffFct and numOfUserFct give the number of functions of each type that
   * make up the problem description.
   */
  INT numOfCoeffFct;

  /** \brief Number of User functions
   *
   * User definable coefficient functions come in two flavours.
   * They are either of type CoeffProcPtr or of type UserProcPtr.
   * numOfCoeffFct and numOfUserFct give the number of functions of each type that
   * make up the problem description.
   */
  INT numOfUserFct;

  /** \brief Coefficient functions
   *
   *  Array that stores the pointers to coefficient and user functions.
   * Since access to this array is provided through macros (see below) the layout
   * is not important. Note that this array is allocated dynamically to the desired length.
   */
  void *CU_ProcPtr[1];
};

/** \brief ???
 *
 * \todo Please doc me!
 */
struct bndcond {

  /** \brief Field for environment variable */
  NS_PREFIX ENVVAR v;

  /* fields for boundary condition */
  /** \brief Corresponds to boundary segment id ! */
  INT id;

  /** \brief Function defining boundary condition */
  BndCondProcPtr BndCond;

  /** \brief Additional data for bnd cond */
  void *data;
};

/****************************************************************************/
/*                                                                          */
/* BoundaryValueProblem data structure                                      */
/*                                                                          */
/****************************************************************************/
/** \brief ???
 *
 * \todo Please doc me!
 */
struct std_BoundaryValueProblem
{
  /** \brief Fields for environment directory */
  NS_PREFIX ENVDIR d;

  /* init */
  INT type;

  /** \brief Domain pointer                      */
  struct domain *Domain;

  /** \brief Problem pointer                     */
  struct problem *Problem;

  /** \brief File name for boundary infos        */
  char bnd_file[NS_PREFIX NAMESIZE];

  /** \brief File name for meshinfos             */
  char mesh_file[NS_PREFIX NAMESIZE];

  /** @name Domain part */
  /*@{*/
  /** \brief Center of a sphere containing the domain */
  DOUBLE MidPoint[DIM];

  /** \brief Radius of a sphere containing the domain */
  DOUBLE radius;

  /** \brief 1 if domain is convex, 0 if not     */
  INT domConvex;

  /** \brief Number of subdomains, exterior not counted                */
  INT numOfSubdomains;

  /** \brief Number of parts in the domain               */
  INT nDomainParts;

  /** \brief Pointer to table subbdom --> part   */
  INT *s2p;
  /*@}*/

  /** @name Boundary decription */
  /*@{*/
  INT ncorners;
  INT nsides;
  INT sideoffset;

  /** \brief list of patches */
  union patch **patches;
  /*@}*/

  /** @name Problem part */
  /*@{*/
  /** \brief Configuration function              */
  ConfigProcPtr ConfigProc;

  /** \brief Number of coefficient functions        */
  INT numOfCoeffFct;

  /** \brief Number of user functions               */
  INT numOfUserFct;

  /** \brief General bnd. cond. (if exists)      */
  BndCondProcPtr GeneralBndCond;

  /** \brief Coefficient functions                           */
  void *CU_ProcPtr[1];
  /*@}*/
};

/****************************************************************************/
/*                                                                          */
/* Patch data structure                                                     */
/*                                                                          */
/****************************************************************************/
/** \brief ???
 *
 * \todo Please doc me!
 */
struct generic_patch {


  /** \brief Patch type */
  enum PatchType type;

  /** \brief Fixed/bnd of free/free */
  INT state;

  /** \brief Unique id used for load/store */
  INT id;
};

/** \brief ???
 *
 * \todo Please doc me!
 */
struct point_on_patch {

  INT patch_id;
  INT corner_id;
};
/** \brief ???
 *
 * \todo Please doc me!
 */
struct point_patch {

  /** \brief Patch type */
  enum PatchType type;

  /** \brief Fixed/bnd of free/free */
  INT state;

  /** \brief Unique id used for load/store */
  INT id;

  /** \brief Number of patches */
  INT npatches;

  /** \brief Reference to surface */
  struct point_on_patch pop[1];
};

#ifdef __THREEDIM__
/** \brief ???
 *
 * \todo Please doc me!
 */
struct line_on_patch {

  INT patch_id;
  INT corner_id[2];
};
/** \brief ???
 *
 * \todo Please doc me!
 */
struct line_patch {

  /** \brief Patch type */
  enum PatchType type;

  /** \brief Fixed/bnd of free/free */
  INT state;

  /** \brief Unique id used for load/store */
  INT id;

  /** \brief Number of patches */
  INT npatches;

  /** \brief Corner 0 of line */
  INT c0;

  /** \brief Corner 1 of line */
  INT c1;

  /** \brief Reference to surface */
  struct line_on_patch lop[1];
};
#endif
/** \brief ???
 *
 * \todo Please doc me!
 */
struct linear_patch {

  /** \brief Patch type */
  enum PatchType type;

  /** \brief Fixed/bnd of free/free */
  INT state;

  /** \brief Unique id used for load/store */
  INT id;

  /** \brief Id of left and right subdomain */
  INT left,right;

  /** \brief Number of corners */
  INT corners;

  /** \brief Ids of points */
  INT points[MAX_CORNERS_OF_LINEAR_PATCH];

  /** \brief Position */
  DOUBLE pos[MAX_CORNERS_OF_LINEAR_PATCH][DIM];
};
/** \brief ???
 *
 * \todo Please doc me!
 */
struct parameter_patch {

  /** \brief Patch type */
  enum PatchType type;

  /** \brief Fixed/bnd of free/free */
  INT state;

  /** \brief Unique id used for load/store */
  INT id;

  /** \brief Id of left and right subdomain */
  INT left,right;

  /** \brief Measure for curvature */
  INT resolution;

  /** \brief Ids of points */
  INT points[CORNERS_OF_BND_SEG];

  /** \brief Parameter range */
  DOUBLE range[2][DIM_OF_BND];

  /** \brief Pointer to definition function */
  BndSegFuncPtr BndSegFunc;

  /** \brief Can be used by applic to find data */
  void *bs_data;

  /** @name Fields for boundary condition */
  /*@{*/
  /** \brief Function defining boundary condition */
  BndCondProcPtr BndCond;

  /** \brief Additional data for bnd cond */
  void *bc_data;
  /*@}*/
};

/** \brief ???
 *
 * \todo Please doc me!
 */
struct marc_0_patch {

  /** \brief Patch type                           */
  enum PatchType type;

  /** \brief Fixed/bnd of free/free               */
  INT state;

  /** \brief Unique id used for load/store        */
  INT id;

  /** \brief Position                             */
  DOUBLE pos[3];
};

/** \brief ???
 *
 * \todo Please doc me!
 */
struct marc_1_patch {

  /** \brief Patch type                           */
  enum PatchType type;

  /** \brief Fixed/bnd of free/free               */
  INT state;

  /** \brief Unique id used for load/store        */
  INT id;

  /** \brief Line between two points              */
  INT p[2];
};

/** \brief ???
 *
 * \todo Please doc me!
 */
struct marc_2_patch {

  /** \brief Patch type                           */
  enum PatchType type;

  /** \brief Fixed/bnd of free/free               */
  INT state;

  /** \brief Unique id used for load/store        */
  INT id;

  /** \brief Bnd cond                             */
  INT c;

  /** \brief Triangle of three points             */
  INT p[3];
};

/** \brief ???
 *
 * \todo Please doc me!
 */
union patch {
  struct generic_patch ge;
  struct point_patch po;
  struct linear_patch lp;
  struct parameter_patch pa;
    #ifdef __THREEDIM__
  struct line_patch li;
        #endif
  struct marc_0_patch m0;
  struct marc_1_patch m1;
  struct marc_2_patch m2;
} ;

/** \brief ???
 *
 * \todo Please doc me!
 */
struct bnd_ps {

  /** \brief Associated patch                     */
  INT patch_id;

  /** \brief E.g. global coordinates, pointers... */
  void *data;

  /** \brief Number of arguments                  */
  INT n;

  /** \brief Parameter range                      */
  COORD_BND_VECTOR local[1];
};
/** \brief ???
 *
 * \todo Please doc me!
 */
struct marc_bndp {

  /** \brief Associated patch                     */
  INT patch_id;

  /** \brief Position                             */
  DOUBLE pos[3];
};

/** \brief ???
 *
 * \todo Please doc me!
 */
struct marc_bnds {

  /** \brief associated patch */
  INT patch_id;

  /** \brief number of corners */
  INT n;

  /** \brief corners */
  struct marc_bndp p[1];
};

/*----------- typedef for structs ------------------------------------------*/

/* typedefs */
#undef DOMAIN
typedef struct domain DOMAIN;
typedef struct linear_segment LINEAR_SEGMENT;
typedef struct boundary_segment BOUNDARY_SEGMENT;
typedef struct problem PROBLEM;
typedef struct bndcond BOUNDARY_CONDITION;

typedef struct std_BoundaryValueProblem STD_BVP;
typedef union patch PATCH;
typedef struct point_patch POINT_PATCH;
typedef struct line_patch LINE_PATCH;
typedef struct linear_patch LINEAR_PATCH;
typedef struct parameter_patch PARAMETER_PATCH;
typedef struct marc_0_patch M0_PATCH;
typedef struct marc_1_patch M1_PATCH;
typedef struct marc_2_patch M2_PATCH;
typedef struct bnd_ps BND_PS;
typedef struct marc_bndp M_BNDP;
typedef struct marc_bnds M_BNDS;

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

# ifdef __THREEDIM__
INT RepairMesh (NS_PREFIX HEAP *Heap, INT MarkKey, MESH *mesh);
INT CheckPrisms (INT *corner, INT n0, INT n1 , INT n2, INT n3);
# endif
void SetBVPType(INT type);

DOMAIN *GetDomain                           (char *name);

LINEAR_SEGMENT *CreateLinearSegment (char *name,
                                     INT left, INT right,INT id,
                                     INT n, INT *point,
                                     DOUBLE x[MAX_CORNERS_OF_LINEAR_PATCH][DIM_OF_BND]);

/* BVP definition */
INT   STD_BVP_Configure           (INT argc, char **argv);

/* scanning of coordinates */
INT   ReadAndPrintArgvPosition    (char *name, INT argc, char **argv, DOUBLE *pos);

END_NAMESPACE

#endif
