// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  project.c                                                                                                     */
/*																			*/
/* Purpose:   projection into subspaces                                                                 */
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   Nov 5, 1997 begin                                                                 */
/*																			*/
/* Remarks:   not finished!                                                                     */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "devices.h"
#include "ugenv.h"

#include "scan.h"
#include "numproc.h"
#include "np.h"
#include "ugm.h"
#include "general.h"
#include "fileopen.h"
#include "ugstruct.h"

#include "commands.h"
#include "assemble.h"

#include "project.h"
#include "shapes.h"
#include "evm.h"
#include "gm.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  project the solution of an eigenvalue-calculation in a space          */
/*		  where the eigenvalues aren't zero.                                                            */
/*		                                                                                                                                */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct
{
  NP_PROJECT prj;

  NP_NL_ASSEMBLE *Assemble;

  VECDATA_DESC *t;
  VECDATA_DESC *b;

} NP_PRJ;

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

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

INT Project_Init (NP_PROJECT *theNP, INT argc, char **argv)
{
  MULTIGRID *theMG;
  theMG = theNP->base.mg ;

  /* assign x, if calling the generic execute routine */

  theNP->x =  ReadArgvVecDesc(theMG,"x",argc,argv);

  return (NP_ACTIVE);
}

INT Project_Display (NP_PROJECT *theNP)
{

  if ((theNP->x) != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"x",ENVITEM_NAME(theNP->x));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"x","no extern assign");

  return (0);
}

/****************************************************************************/
/*D
   project - projection num proc

   DESCRIPTION:

   'npinit <name> $P <equation>'

   .  <name> - num proc name
   .  $P <equation> - name of the equation e.g. pln : Laplace-Equation
                                             ppn : Plate-Equation
                                             pen : Elasticity-Equation
   D*/
/****************************************************************************/

static INT Prj_Init (NP_BASE *theNP, INT argc, char **argv)
{
  NP_PRJ *np = (NP_PRJ *)theNP;

  np->Assemble = (NP_NL_ASSEMBLE *)
                 ReadArgvNumProc(theNP->mg,"A",NL_ASSEMBLE_CLASS_NAME,argc,argv);

  if (np->Assemble == NULL)
    return(NP_NOT_ACTIVE);

  return(Project_Init(&np->prj,argc,argv));
}

static INT Prj_Display (NP_BASE *theNP)
{
  NP_PRJ *np = (NP_PRJ *)theNP;

  UserWrite("configuration parameters:\n");
  if (np->Assemble != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Assemble",ENVITEM_NAME(np->Assemble));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"Assemble","no Assemble-Routine");

  Project_Display(&np->prj);

  return (0);
}


/* Projection for Laplace-Equation */

static INT ProjectionVectorLN (NP_PROJECT *theNP, INT fl, INT tl,
                               INT index, VECDATA_DESC *t, INT *result)
{
  NP_PRJ    *np    = (NP_PRJ *) theNP;
  MULTIGRID *theMG = theNP->base.mg;

  if (index != 0)
    NP_RETURN(1,result[0]);

  /* assign t to the one-vector */

  if (dset(theMG,fl,tl,ALL_VECTORS,t,1.0) != NUM_OK)
    NP_RETURN(1,result[0]);

  return(0);
}


/*  Old and false Projection for Plate-Equation  */

static INT ProjectionVectorPL1 (NP_PROJECT *theNP, INT fl, INT tl,
                                INT index, VECDATA_DESC *t, INT *result)
{
  NP_PRJ    *np    = (NP_PRJ *) theNP;
  MULTIGRID *theMG = theNP->base.mg;
  VECTOR *v;
  DOUBLE_VECTOR pos;
  INT lev,vtype,ncomp,comp;

  for (lev=fl; lev<=tl; lev++)
    for (v=FIRSTVECTOR(GRID_ON_LEVEL(theMG,lev)); v!=NULL; v=SUCCVC(v)) {
      vtype = VTYPE(v);
      ncomp = VD_NCMPS_IN_TYPE(t,vtype);
      if (ncomp == 0) continue;
      VectorPosition(v,pos);
      comp = VD_CMP_OF_TYPE(t,vtype,0);
      switch(index)
      {
      case 0 :
        VVALUE(v,comp) = 1.0 ;
        continue;
      case 1 :
        VVALUE(v,comp) = pos[0] ;
        continue;
      case 2 :
        VVALUE(v,comp) = pos[1] ;
        continue;
      case 3 :
        VVALUE(v,comp) = pos[0]*pos[0] ;
        continue;
      case 4 :
        VVALUE(v,comp) = pos[0]*pos[1] ;
        continue;
      case 5 :
        VVALUE(v,comp) = pos[1]*pos[1] ;
        continue;
      }
    }
  return(0);
}


/* modified projection for the plate-equation */
/* the node and middle point are differently valued */

static Normale(INT n, INT i, DOUBLE **x, DOUBLE_VECTOR normal)
{
  DOUBLE s ;
  DOUBLE_VECTOR y;

  V2_SUBTRACT(x[(i+1)%n],x[i],y);
  s = 0.0 ;
  s = sqrt(y[0]*y[0] + y[1]* y[1]);
  V2_SCALE(1.0/s,y);
  V2_SUBTRACT(x[(i+1)%n],x[(i+2)%n],normal);
  V2_SCALAR_PRODUCT(normal,y,s);
  V2_LINCOMB(1.0,normal,-s,y,normal);
  s = 0.0 ;
  s = sqrt(normal[0]*normal[0] + normal[1]* normal[1] );
  V2_SCALE(1.0/s,normal);

  return(0);

}


#define ORIENTATION(n,s,pos)         (((pos)[(s)][0]>(pos)[((s)+1)%n][0]) ?        \
                                      1.0 :                                 \
                                      ( ((pos)[(s)][0]==(pos)[((s)+1)%n][0])    \
                                        ?                                    \
                                        ( ((pos)[(s)][1]>(pos)[((s)+1)%n][1])    \
                                                                               ? 1.0 : -1.0 )                      \
                                        : -1.0 ) )


static INT ProjectionVectorPL (NP_PROJECT *theNP, INT fl, INT tl,
                               INT index, VECDATA_DESC *t, INT *result)
{
  NP_PRJ    *np    = (NP_PRJ *) theNP;
  MULTIGRID *theMG = theNP->base.mg;
  ELEMENT *theElement;
  VECTOR *v;
  DOUBLE_VECTOR pos ;
  DOUBLE_VECTOR normal;
  DOUBLE s;
  INT lev,vtype,ncomp,comp,side,n;
  DOUBLE *x_pos[3] ;

  for (lev=fl; lev<=tl; lev++) {
    for (v=FIRSTVECTOR(GRID_ON_LEVEL(theMG,lev)); v!=NULL; v=SUCCVC(v)) {
      if (VTYPE(v) == 1) continue;
      vtype = VTYPE(v);
      ncomp = VD_NCMPS_IN_TYPE(t,vtype);
      comp = VD_CMP_OF_TYPE(t,vtype,0);
      if (ncomp == 0) continue;
      VectorPosition(v,pos);
      switch(index)
      {
      case 0 :
        VVALUE(v,comp) = 1.0 ;
        continue;
      case 1 :
        VVALUE(v,comp) = pos[0] ;
        continue;
      case 2 :
        VVALUE(v,comp) = pos[1] ;
        continue;
      case 3 :
        VVALUE(v,comp) = pos[0]*pos[0] ;
        continue;
      case 4 :
        VVALUE(v,comp) = pos[0]*pos[1] ;
        continue;
      case 5 :
        VVALUE(v,comp) = pos[1]*pos[1] ;
        continue;
      }
    }
    for (theElement=FIRSTELEMENT(GRID_ON_LEVEL(theMG,lev)); theElement!=NULL;
         theElement=SUCCE(theElement))
      for (side=0; side<EDGES_OF_ELEM(theElement); side++) {
        v = EDVECTOR(GetEdge(CORNER(theElement,
                                    CORNER_OF_EDGE(theElement,side,0)),
                             CORNER(theElement,
                                    CORNER_OF_EDGE(theElement,side,1))));
        CORNER_COORDINATES(theElement,n,x_pos);
        VectorPosition(v,pos);
        s = ORIENTATION(n,side,x_pos);
        comp = VD_CMP_OF_TYPE(t,EDGEVEC,0);
        Normale(n,side,x_pos,normal);

        switch(index)
        {
        case 0 :
          VVALUE(v,comp) = 0.0;
          continue;
        case 1 :
          VVALUE(v,comp) = normal[0] * s ;
          continue;
        case 2 :
          VVALUE(v,comp) = normal[1] * s ;
          continue;
        case 3 :
          VVALUE(v,comp) = 2* pos[0]* normal[0] * s ;
          continue;
        case 4 :
          VVALUE(v,comp) = (pos[1]* normal[0] + pos[0]* normal[1]) * s ;
          continue;
        case 5 :
          VVALUE(v,comp) =  2* pos[1]* normal[1] * s ;
          continue;
        }
      }
  }

  return(0);
}



/* Projection for Elasticity-Equation */

static INT ProjectionVectorEL (NP_PROJECT *theNP, INT fl, INT tl,
                               INT index, VECDATA_DESC *t, INT *result)
{
  NP_PRJ    *np    = (NP_PRJ *) theNP;
  MULTIGRID *theMG = theNP->base.mg;
  VECTOR *v;
  DOUBLE_VECTOR pos;
  INT lev,vtype,ncomp,comp;

  for (lev=fl; lev<=tl; lev++)
    for (v=FIRSTVECTOR(GRID_ON_LEVEL(theMG,lev)); v!=NULL; v=SUCCVC(v)) {
      vtype = VTYPE(v);
      ncomp = VD_NCMPS_IN_TYPE(t,vtype);
      if (ncomp == 0) continue;
      ASSERT(ncomp == DIM);
      VectorPosition(v,pos);
      comp = VD_CMP_OF_TYPE(t,vtype,0);

      switch(index)
      {
      case 0 :
        VVALUE(v,comp) = 1 ;
        VVALUE(v,comp+1) = 0 ;
        VVALUE(v,comp+2) = 0 ;
        continue;
      case 1 :
        VVALUE(v,comp) = 0 ;
        VVALUE(v,comp+1) = 1 ;
        VVALUE(v,comp+2) = 0 ;
        continue ;
      case 2 :
        VVALUE(v,comp) = 0 ;
        VVALUE(v,comp+1) = 0 ;
        VVALUE(v,comp+2) = 1 ;
        continue;
      case 3 :
        VVALUE(v,comp) = - pos[1] ;
        VVALUE(v,comp+1) = pos[0] ;
        VVALUE(v,comp+2) = 0 ;
        continue ;
      case 4 :
        VVALUE(v,comp) = pos[2] ;
        VVALUE(v,comp+1) = 0 ;
        VVALUE(v,comp+2) = - pos[0] ;
        continue ;
      case 5 :
        VVALUE(v,comp) = 0 ;
        VVALUE(v,comp+1) = - pos[2] ;
        VVALUE(v,comp+2) = pos[1] ;
        continue ;
      }

    }
  return(0);
}

static INT Projection (NP_PROJECT *theNP, INT fl, INT tl,
                       VECDATA_DESC *x, INT *result)
{
  NP_PRJ    *np    = (NP_PRJ *) theNP;
  MULTIGRID *theMG = theNP->base.mg;
  DOUBLE a0,a1;
  INT i;

  /* assign x, if extern evaluable */

  if ((theNP->x) != NULL)
    x = theNP->x ;

  np->t = NULL;
  np->b = NULL;

  if (AllocVDFromVD(theMG,fl,tl,x,&np->t))
    NP_RETURN(1,result[0]);
  if (AllocVDFromVD(theMG,fl,tl,x,&np->b))
    NP_RETURN(1,result[0]);
  for (i=0; i<theNP->dim; i++) {
    if ((*theNP->ProjectionVector)(theNP,fl,tl,i,np->t,result))
      NP_RETURN(1,result[0]);
    if ((*np->Assemble->NLAssembleDefect)
          (np->Assemble,fl,tl,np->t,np->b,NULL,result))
      RETURN(1);
    if (ddot(theMG,fl,tl,ON_SURFACE,np->t,np->b,&a0) != NUM_OK)
      return(1);
    if (ddot(theMG,fl,tl,ON_SURFACE,x,np->b,&a1) != NUM_OK)
      return(1);
    ASSERT(a0 != 0.0);
    if (daxpy(theMG,fl,tl,ALL_VECTORS,x,-a1/a0,np->t) != NUM_OK)
      return(1);
  }
  FreeVD(theMG,fl,tl,np->t);
  FreeVD(theMG,fl,tl,np->b);

  return(0);
}

static INT PLN_Construct (NP_BASE *theNP)
{
  NP_PROJECT *np;

  theNP->Init = Prj_Init;
  theNP->Display = Prj_Display;
  theNP->Execute = NULL;

  np = (NP_PROJECT *) theNP;
  np->PreProcess = NULL;
  np->Project = Projection ;
  np->ProjectionVector = ProjectionVectorLN ;
  np->PostProcess = NULL;

  np->dim = 1;

  return(0);
}

static INT PPN_Construct (NP_BASE *theNP)
{
  NP_PROJECT *np;

  theNP->Init = Prj_Init;
  theNP->Display = Prj_Display;
  theNP->Execute = NULL;

  np = (NP_PROJECT *) theNP;
  np->PreProcess = NULL;
  np->Project = Projection;
  np->ProjectionVector = ProjectionVectorPL ;
  np->PostProcess = NULL;

  /*	np->dim = 6; */
  np->dim = 3;

  return(0);
}

static INT PEN_Construct (NP_BASE *theNP)
{
  NP_PROJECT *np;

  theNP->Init = Prj_Init;
  theNP->Display = Prj_Display;
  theNP->Execute = NULL;

  np = (NP_PROJECT *) theNP;
  np->PreProcess = NULL;
  np->Project = Projection ;
  np->ProjectionVector = ProjectionVectorEL;
  np->PostProcess = NULL;

  np->dim = 6;

  return(0);
}

/****************************************************************************/
/*D
   InitProject - Enrol project num procs

   SYNOPSIS:
   INT InitProject (void);

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function creates the numproc 'project'.
   It is called in initnp.c.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT InitProject (void)
{

  if (CreateClass(PROJECT_CLASS_NAME ".pln",sizeof(NP_PRJ),
                  PLN_Construct))
    return (__LINE__);
  if (CreateClass(PROJECT_CLASS_NAME ".ppn",sizeof(NP_PRJ),
                  PPN_Construct))
    return (__LINE__);
  if (CreateClass(PROJECT_CLASS_NAME ".pen",sizeof(NP_PRJ),
                  PEN_Construct))
    return (__LINE__);

  return (0);
}
