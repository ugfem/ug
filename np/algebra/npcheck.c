// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      npcheck.c                                                     */
/*                                                                          */
/* Purpose:   check of numerical structures                                                     */
/*                                                                          */
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de			                                */
/*																			*/
/* History:   Juli 1 97 begin                                                                           */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <math.h>

#include "devices.h"
#include "compiler.h"
#include "gm.h"
#include "np.h"
#include "disctools.h"
#include "debug.h"
#include "general.h"

#include "pargm.h"
#ifdef ModelP
#include "parallel.h"
#endif

#include "npcheck.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

#ifdef ModelP
static INT pcheck;
#endif

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

INT CheckSymmetryOfMatrix (GRID *theGrid, MATDATA_DESC *A)
{
  MATRIX *m,*mt;
  VECTOR *v,*w;
  DOUBLE *mptr,*mtptr;
  register SHORT i,j,rcomp,ccomp,*mcomp,*mtcomp,vtype,mtype;

  for (v=FIRSTVECTOR(theGrid); v!=NULL; v=SUCCVC(v)) {
    vtype = VTYPE(v);
    for (m=VSTART(v); m!=NULL; m=MNEXT(m)) {
      mt = MADJ(m);
      w = MDEST(m);
      mtype = MTP(vtype,VTYPE(w));
      rcomp = MD_ROWS_IN_MTYPE(A,mtype);
      if (rcomp == 0) continue;
      ccomp = MD_COLS_IN_MTYPE(A,mtype);
      if (ccomp == 0) continue;
      mcomp = MD_MCMPPTR_OF_MTYPE(A,mtype);
      mptr = MVALUEPTR(m,0);
      mtcomp = MD_MCMPPTR_OF_MTYPE(A,MTP(VTYPE(w),vtype));
      mtptr = MVALUEPTR(m,0);
      for (i=0; i<ccomp; i++)
        for (j=0; j<rcomp; j++)
          if (mptr[mcomp[i*rcomp+j]] != mtptr[mtcomp[j*ccomp+i]])
            return(1);
    }
  }

  return(0);
}

static INT CheckVector (GRID *theGrid, VECTOR *v)
{
  FORMAT *theFormat;
  NODE *theNode;
  VECTOR *w;
  INT nerr = 0;

  /* get format */
  theFormat = MGFORMAT(MYMG(theGrid));
  if (FMT_S_MATPTR(theFormat)[MatrixType[VTYPE(v)][VTYPE(v)]]
      && (!GHOST(v))) {
    if (VSTART(v) == NULL) {
      UserWriteF(PFMT "ERROR: no diagonal matrix vec=" VINDEX_FMTX "\n",
                 me,VINDEX_PRTX(v));
      nerr++;
    }
    else if (!MDIAG(VSTART(v))) {
      UserWriteF(PFMT "ERROR: VSTART no diagonal matrix vec="
                 VINDEX_FMTX "\n",
                 me,VINDEX_PRTX(v));
      nerr++;
    }
  }

  /* check flags locally */
  if (NEW_DEFECT(v) != (VCLASS(v)>=2)) {
    UserWriteF(PFMT "ERROR: classes not match vec="
               VINDEX_FMTX " NEW_DEFECT %d VCLASS %d\n",
               me,VINDEX_PRTX(v),NEW_DEFECT(v),VCLASS(v));
    nerr++;
  }
  if (FINE_GRID_DOF(v) != ((VCLASS(v)>=2)&&(VNCLASS(v)<=1))) {
    UserWriteF(PFMT "ERROR: classes not match vec="
               VINDEX_FMTX " FINE_GRID_DOF %d VNCLASS %d VCLASS %d\n",
               me,VINDEX_PRTX(v),FINE_GRID_DOF(v),VNCLASS(v),VCLASS(v));
    nerr++;
  }
  if (FINE_GRID_DOF(v))
    if (FULLREFINELEVEL(MYMG(theGrid)) > GLEVEL(theGrid)) {
      UserWriteF(PFMT "ERROR: FULLREFINELEVEL too large vec="
                 VINDEX_FMTX " FINE_GRID_DOF %d FULLREFINELEVEL %d\n",
                 me,VINDEX_PRTX(v),FINE_GRID_DOF(v),
                 FULLREFINELEVEL(MYMG(theGrid)));
      nerr++;
    }
  if (VOTYPE(v) == NODEVEC) {
    theNode = (NODE *) VOBJECT(v);
    if (theNode == NULL) {
      if (GLEVEL(theGrid) >= 0) {
        UserWriteF(PFMT "ERROR: nodevector has no NODE vec="
                   VINDEX_FMTX " \n",
                   me,VINDEX_PRTX(v));
        nerr++;
      }
    }
    else {
      if (OBJT(theNode) != NDOBJ) {
        UserWriteF(PFMT "ERROR: nodevector has no NODE object vec="
                   VINDEX_FMTX " OBJT %d\n",
                   me,VINDEX_PRTX(v),OBJT(theNode));
        nerr++;
      }
      if (NTYPE(theNode) == CORNER_NODE) {
        theNode = (NODE *)NFATHER(theNode);
        if (theNode != NULL) {
          w = NVECTOR(theNode);
          if (w == NULL) {
            UserWriteF(PFMT "ERROR:"
                       " cornernode vector has no father vec="
                       VINDEX_FMTX "\n",
                       me,VINDEX_PRTX(v));
            nerr++;
          }
          if (VNCLASS(w) != VCLASS(v)) {
            UserWriteF(PFMT "ERROR:"
                       " VCLASS and VNCLASS not matches vec="
                       VINDEX_FMTX " VCLASS %d father vec "
                       VINDEX_FMTX " VNCLASS %d\n",
                       me,VINDEX_PRTX(v),VCLASS(v),
                       VINDEX_PRTX(w),VNCLASS(w));
            nerr++;
          }
        }
      }
    }
  }
  return(nerr);
}

static INT CheckVectors (GRID *theGrid)
{
  INT nerr = 0;
  VECTOR *theVector;

  for (theVector=PFIRSTVECTOR(theGrid); theVector!=NULL;
       theVector=SUCCVC(theVector)) {
    nerr += CheckVector(theGrid,theVector);
  }
  return(nerr);
}

#ifdef ModelP
static int Gather_VectorFlags (DDD_OBJ obj, void *data)
{
  VECTOR *pv = (VECTOR *)obj;
  INT *idata = (INT *)data;

  idata[0] = VECSKIP(pv);
  idata[1] = VCLASS(pv);
  idata[2] = VNCLASS(pv);
  idata[3] = NEW_DEFECT(pv);
  idata[4] = FINE_GRID_DOF(pv);
  idata[5] = VTYPE(pv);
  idata[6] = VOTYPE(pv);
  idata[7] = VDATATYPE(pv);
  idata[8] = VNEW(pv);
  idata[9] = VECTORSIDE(pv);
  idata[10] = VPART(pv);

  return (0);
}

static int Scatter_VectorFlags (DDD_OBJ obj, void *data)
{
  VECTOR *pv = (VECTOR *)obj;
  INT *idata = (INT *)data;

  if (idata[0] != VECSKIP(pv)) {
    printf(PFMT "ERROR:"
           " VECSKIP not matches vec="
           VINDEX_FMTX " %d master %d\n",
           VINDEX_PRTX(pv),VECSKIP(pv),idata[0]);
    pcheck++;
    ASSERT(0);
  }
  if (idata[1] != VCLASS(pv)) {
    printf(PFMT "ERROR:"
           " VCLASS not matches vec="
           VINDEX_FMTX " %d master %d\n",
           VINDEX_PRTX(pv),VCLASS(pv),idata[1]);
    ASSERT(0);
  }
  if (idata[2] != VNCLASS(pv)) {
    printf(PFMT "ERROR:"
           " VNCLASS not matches vec="
           VINDEX_FMTX " %d master %d\n",
           VINDEX_PRTX(pv),VNCLASS(pv),idata[2]);
    pcheck++;
    ASSERT(0);
  }
  if (idata[3] != NEW_DEFECT(pv)) {
    printf(PFMT "ERROR:"
           " NEW_DEFECT not matches vec="
           VINDEX_FMTX " %d master %d\n",
           VINDEX_PRTX(pv),NEW_DEFECT(pv),idata[3]);
    pcheck++;
    ASSERT(0);
  }
  if (idata[4] != FINE_GRID_DOF(pv)) {
    printf(PFMT "ERROR:"
           " FINE_GRID_DOF not matches vec="
           VINDEX_FMTX " %d master %d\n",
           VINDEX_PRTX(pv),FINE_GRID_DOF(pv),idata[4]);
    pcheck++;
    ASSERT(0);
  }
  if (idata[5] != VTYPE(pv)) {
    printf(PFMT "ERROR:"
           " VTYPE not matches vec="
           VINDEX_FMTX " %d master %d\n",
           VINDEX_PRTX(pv),VTYPE(pv),idata[5]);
    pcheck++;
    ASSERT(0);
  }
  if (idata[6] != VOTYPE(pv)) {
    printf(PFMT "ERROR:"
           " VOTYPE not matches vec="
           VINDEX_FMTX " %d master %d\n",
           VINDEX_PRTX(pv),VOTYPE(pv),idata[6]);
    pcheck++;
    ASSERT(0);
  }
  if (idata[7] != VDATATYPE(pv)) {
    printf(PFMT "ERROR:"
           " VDATATYPE not matches vec="
           VINDEX_FMTX " %d master %d\n",
           VINDEX_PRTX(pv),VDATATYPE(pv),idata[7]);
    pcheck++;
    ASSERT(0);
  }
  if (idata[8] != VNEW(pv)) {
    printf(PFMT "ERROR:"
           " VNEW not matches vec="
           VINDEX_FMTX " %d master %d\n",
           VINDEX_PRTX(pv),VNEW(pv),idata[8]);
    pcheck++;
    ASSERT(0);
  }
  if (idata[9] != VECTORSIDE(pv)) {
    printf(PFMT "ERROR:"
           " VECTORSIDE not matches vec="
           VINDEX_FMTX " %d master %d\n",
           VINDEX_PRTX(pv),VECTORSIDE(pv),idata[9]);
    pcheck++;
    ASSERT(0);
  }
  if (idata[10] != VPART(pv)) {
    printf(PFMT "ERROR:"
           " VPART not matches vec="
           VINDEX_FMTX " %d master %d\n",
           VINDEX_PRTX(pv),VPART(pv),idata[10]);
    pcheck++;
    ASSERT(0);
  }

  return (0);
}
#endif

INT CheckNP (MULTIGRID *theMG, INT argc, char **argv)
{
  MATDATA_DESC *A;
  VECDATA_DESC *x,*y;
  INT i,level,nerr;
  char value[VALUELEN];
  VEC_SCALAR damp;
  DOUBLE nrm,diff;

  if (ReadArgvChar("A",value,argc,argv) == 0) {
    A = GetMatDataDescByName(theMG,value);
    if (A == NULL) {
      UserWriteF("ERROR: no matrix %s in npckeck\n",value);
      return(1);
    }
    if (ReadArgvOption("S",argc,argv)) {
      for (level=theMG->bottomLevel; level<=TOPLEVEL(theMG); level++)
        if (CheckSymmetryOfMatrix(GRID_ON_LEVEL(theMG,level),A))
          UserWriteF("matrix %s not symmetric on level %d\n",
                     ENVITEM_NAME(A),level);
      return(0);
    }
    if (ReadArgvOption("G",argc,argv)) {
      if (ReadArgvChar("x",value,argc,argv)) {
        UserWriteF("ERROR: no vector in npckeck\n");
        return(1);
      }
      x = GetVecDataDescByName(theMG,value);
      if (x == NULL) {
        UserWriteF("ERROR: no vector %s in npckeck\n",value);
        return(1);
      }
      level = CURRENTLEVEL(theMG);
      if (level == BOTTOMLEVEL(theMG)) {
        UserWriteF("ERROR: no GalerkinCheck,"
                   "level %d is bottomlevel\n",level);
        return(1);
      }
      if (AllocVDFromVD(theMG,level-1,level,x,&y))
        return(1);
      dmatset(theMG,level-1,level-1,ALL_VECTORS,A,0.0);
      dset(theMG,level,level,ALL_VECTORS,x,0.0);
      dset(theMG,level-1,level,ALL_VECTORS,y,0.0);
      AssembleGalerkinByMatrix(GRID_ON_LEVEL(theMG,level),A,0);
      for (i=0; i<VD_NCOMP(x); i++) damp[i] = 1.0;
      InterpolateCorrectionByMatrix(GRID_ON_LEVEL(theMG,level),x,x,damp);
      if (dmatmul(theMG,level,level,ALL_VECTORS,y,A,x) != NUM_OK)
        return(1);
      RestrictByMatrix(GRID_ON_LEVEL(theMG,level),y,y,damp);
      IFDEBUG(np,1)
      UserWriteF("x %d\n",level-1);
      PrintVector(GRID_ON_LEVEL(theMG,level-1),x,3,3);
      UserWriteF("x %d\n",level);
      PrintVector(GRID_ON_LEVEL(theMG,level),x,3,3);
      UserWriteF("y %d\n",level);
      PrintVector(GRID_ON_LEVEL(theMG,level),y,3,3);
      UserWriteF("y %d\n",level-1);
      PrintVector(GRID_ON_LEVEL(theMG,level-1),y,3,3);
      ENDDEBUG
      if (dmatmul_minus(theMG,level-1,level-1,ALL_VECTORS,y,A,x)!=NUM_OK)
        return(1);
      IFDEBUG(np,1)
      UserWriteF("y %d\n",level-1);
      PrintVector(GRID_ON_LEVEL(theMG,level-1),y,3,3);
      ENDDEBUG
      dnrm2(theMG,level-1,level-1,ALL_VECTORS,x,&nrm);
      dnrm2(theMG,level-1,level-1,ALL_VECTORS,y,&diff);
      UserWriteF("Galerkin test: nrm(x) = %f nrm(Ax-RAPx) = %f\n",
                 nrm,diff);
      return(0);
    }
  }
  for (level=theMG->bottomLevel; level<=TOPLEVEL(theMG); level++) {
    UserWriteF("[%d: numeric: ",level);
    nerr = CheckVectors(GRID_ON_LEVEL(theMG,level));
        #ifdef ModelP
    nerr = UG_GlobalSumINT(nerr);
        #endif
    if (nerr)
      UserWriteF("ERROR: vector flags not correctly set] ");
    else
      UserWrite("ok] ");
  }
    #ifdef ModelP
  pcheck = 0;
  DDD_IFOneway(VectorVAllIF, IF_FORWARD, 11 * sizeof(INT),
               Gather_VectorFlags, Scatter_VectorFlags);
  pcheck = UG_GlobalSumINT(pcheck);
  if (pcheck == 0)
    UserWriteF("[parallel numeric: ok]");
  else
    UserWriteF("[parallel numeric: %d errors]",pcheck);
        #endif
  UserWrite("\n");

  return(0);
}
