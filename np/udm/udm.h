// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  udm.h															*/
/*																			*/
/* Purpose:   user data manager (header file)								*/
/*																			*/
/* Author:	  Peter Bastian													*/
/*																			*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/*																			*/
/* History:   02.12.96 begin, ug version 3.4								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __UDM__
#define __UDM__

#include "compiler.h"
#include "gm.h"

/****************************************************************************/
/*																			*/
/* data structures                                                                                                                      */
/*																			*/
/****************************************************************************/

#define NVECTYPES                               MAXVECTORS
#define NMATTYPES                               (MAXVECTORS*MAXVECTORS)
#define MAX_SINGLE_VEC_COMP              9      /* max nb of comp in one TYPE                   */
#define MAX_SINGLE_MAT_COMP             81      /* max nb of comp in one TYPE		    */
#define MAX_VEC_COMP                     9      /* max nb of comp in one VECDATA_DESC	*/
#define MAX_MAT_COMP                    81      /* max nb of comp in one VECDATA_DESC	*/

#define NVECOFFSETS                             (NVECTYPES+1)
/* for offset component in VECDATA_DESC	*/
#define NMATOFFSETS                             (NMATTYPES+1)

typedef struct {

  /* fields for environment list variable */
  ENVVAR v;

  INT locked;                          /* locked for dynamic allocation         */
  char compNames[MAX_VEC_COMP];    /* names for symbol components           */
  SHORT NCmpInType[NVECTYPES];     /* number of components of a vector      */
                                   /* per type                              */
  SHORT *CmpsInType[NVECTYPES];    /* pointer to SHORT vector containing    */
  /*    the components                     */
  /* redundant (but frequently used) information                          */
  SHORT IsScalar;                  /* TRUE if desc is scalar:               */
                                   /*  same settings in all types           */
  SHORT ScalComp;                  /* location of scalar component          */
  INT ScalTypeMask;                /* mask for used vectypes                */
  SHORT offset[NVECOFFSETS];       /* offsets for VEC_SCALARs               */
  SHORT Components[1];                 /* memory for component mapping	        */

} VECDATA_DESC;

typedef struct {

  ENVVAR v;

  INT locked;                          /* locked for dynamic allocation         */
  char compNames[2*MAX_MAT_COMP];   /* names for symbol components          */
  SHORT RowsInType[NMATTYPES];          /* number of rows of a matrix per type  */
  SHORT ColsInType[NMATTYPES];          /* number of columns of a matrix        */
                                        /* per type                             */
  SHORT *CmpsInType[NMATTYPES];         /* pointer to SHORT vector containing   */
                                        /* the components                       */
  /* redundant (but frequently used) information                          */
  SHORT IsScalar;                       /* TRUE if desc is scalar:              */
  /* same settings in all types           */
  SHORT ScalComp;                       /* location of scalar component         */
  INT ScalRowTypeMask;                  /* mask for used vectypes in rows       */
  INT ScalColTypeMask;                  /* mask for used vectypes in cols       */
  SHORT offset[NMATOFFSETS];            /* offsets for what ever you need it    */
  SHORT Components[1];                  /* memory for component mapping	        */

} MATDATA_DESC;


typedef DOUBLE VEC_SCALAR[MAX_VEC_COMP];

/****************************************************************************/
/*																			*/
/* definition of exported data structures									*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* definition of exported functions											*/
/*																			*/
/****************************************************************************/

VECDATA_DESC *CreateVecDesc (MULTIGRID *theMG, char *name, char *compNames,
                             SHORT *NCmpInType);
MATDATA_DESC *CreateMatDesc (MULTIGRID *theMG, char *name, char *compNames,
                             SHORT *RowsInType, SHORT *ColsInType);

VECDATA_DESC *GetVecDataDescByName (MULTIGRID *theMG, char *name);
MATDATA_DESC *GetMatDataDescByName (MULTIGRID *theMG, char *name);

INT AllocVDfromVD (MULTIGRID *theMG, INT fl, INT tl,
                   VECDATA_DESC *template_desc, VECDATA_DESC **new_desc);
INT AllocMDfromVD (MULTIGRID *theMG, INT fl, INT tl,
                   VECDATA_DESC *x, VECDATA_DESC *y, MATDATA_DESC **new_desc);
INT AllocMDfromMD (MULTIGRID *theMG, INT fl, INT tl,
                   MATDATA_DESC *template_desc, MATDATA_DESC **new_desc);
INT FreeVD        (MULTIGRID *theMG, INT fl, INT tl, VECDATA_DESC *x);
INT FreeMD        (MULTIGRID *theMG, INT fl, INT tl, MATDATA_DESC *A);

INT ConstructVecOffsets (SHORT *NCmpInType, SHORT *offset);
INT ConstructMatOffsets (SHORT *RowsInType, SHORT *ColsInType, SHORT *offset);


/* init user data manager */
INT InitUserDataManager (void);

#endif
