// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  iter_2.c			                                                                                */
/*																			*/
/* Purpose:   iteration num procs, some more                                    */
/*																			*/
/*																			*/
/* Author:	  Klaus Johannsen                                                                                       */
/*            IWR/Technische Simulation                                     */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            69120 Heidelberg                                              */
/*            email: klaus.johannsen@iwr.uni-heidelberg.de                  */
/*                                                                          */
/* History:   27.02.03 begin                                                */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
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

#include "general.h"
#include "debug.h"
#include "dlmgr.h"
#include "gm.h"
#include "ugm.h"
#include "algebra.h"
#include "scan.h"
#include "numproc.h"
#include "np.h"
#include "ugdevices.h"
#include "udm.h"
#include "pcr.h"
#include "debug.h"
#include "fifo.h"
#include "evm.h"
#include "misc.h"
#include "ugstruct.h"

#include "transfer.h"
#include "ls.h"
#include "iter.h"
#include "project.h"
#include "disctools.h"
#include "block.h"

#include "ff_gen.h"
#include "ff.h"
#include "ugblas.h"
#include "order.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
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
  NP_ITER iter;

  VEC_SCALAR damp;
  DOUBLE alpha,Gamma;
  MATDATA_DESC *B;
} NP_SOR_A;

typedef struct
{
  NP_ITER iter;

  VEC_SCALAR damp;
  DOUBLE alpha,Gamma;
  MATDATA_DESC *B;
} NP_ILU_A;

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

static VEC_SCALAR Factor_One;

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   sora - numproc for SOR_A smoother

   DESCRIPTION:
   This numproc executes an SOR_A (auto-damp) smoother.

   .vb
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>]
       $damp <sc double list> $alpha <double> $Gamma <double>;
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $damp~<sc~double~list> - additional damping factors for each component
   .  $alpha~<double> - alpha
   .  $Gamma~<double> - Gamma

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p];'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT SOR_A_CopyMatrix (MULTIGRID *mg, INT level, MATDATA_DESC *B, MATDATA_DESC *A, DOUBLE alpha, DOUBLE Gamma)
{
  INT i,comp0,comp,nc,nr;
  VECTOR *v;
  MATRIX *m;
  DOUBLE sum;

  if (dmatcopy(mg,level,level,ALL_VECTORS,B,A)) return(1);
  nc = MD_COLS_IN_RT_CT(A,NODEVEC,NODEVEC);
  nr = MD_ROWS_IN_RT_CT(A,NODEVEC,NODEVEC);
  assert(nc==nr);
  comp0=MD_MCMP_OF_RT_CT(B,NODEVEC,NODEVEC,0);
  for (v=FIRSTVECTOR(GRID_ON_LEVEL(mg,level)); v!=NULL; v=SUCCVC(v))
  {
    for (i=0; i<nc; i++)
    {
      /* component of i-th diagonal */
      comp=comp0+i*(nc+1);

      /* set diagonal */
      if (VECSKIP(v)) continue;
      sum=0.0;
      for (m=MNEXT(VSTART(v)); m!=NULL; m=MNEXT(m))
        if (!VECSKIP(MDEST(m)))
          sum+=ABS(MVALUE(m,comp)-MVALUE(MADJ(m),comp));
      sum*=0.25*alpha*Gamma;
      if (MVALUE(VSTART(v),comp)<0.0) sum=-sum;
      MVALUE(VSTART(v),comp)+=sum;

      /* set lower triangle */
      for (m=MNEXT(VSTART(v)); m!=NULL; m=MNEXT(m))
        if (!VECSKIP(MDEST(m)))
          if (MDESTINDEX(m)<VINDEX(v))
          {
            MVALUE(m,comp)=0.5*(1.0+alpha)*MVALUE(m,comp)+0.5*(1.0-alpha)*MVALUE(MADJ(m),comp);
            MVALUE(MADJ(m),comp)=0.0;
          }
    }
  }

  return(0);
}

static INT SOR_F_CopyMatrix (MULTIGRID *mg, INT level, MATDATA_DESC *B, MATDATA_DESC *A)
{
  INT comp;
  VECTOR *v;
  MATRIX *m;
  DOUBLE sum;

  if (dmatcopy(mg,level,level,ALL_VECTORS,B,A)) return(1);
  comp=MD_MCMP_OF_RT_CT(B,NODEVEC,NODEVEC,0);
  for (v=FIRSTVECTOR(GRID_ON_LEVEL(mg,level)); v!=NULL; v=SUCCVC(v))
  {
    /* modify diagonal */
    if (VECSKIP(v)) continue;
    sum=0.0;
    for (m=MNEXT(VSTART(v)); m!=NULL; m=MNEXT(m))
      if (!VECSKIP(MDEST(m)))
        sum+=ABS(MVALUE(m,comp));
    if (sum<=ABS(MVALUE(VSTART(v),comp))) continue;
    if (MVALUE(VSTART(v),comp)<0.0) sum=-sum;
    MVALUE(VSTART(v),comp)=sum;
  }

  return(0);
}

static INT SOR_A_Init (NP_BASE *theNP, INT argc , char **argv)
{
  INT i;
  NP_SOR_A *np=(NP_SOR_A *)theNP;

  for (i=0; i<MAX_VEC_COMP; i++) np->damp[i]=1.0;
  sc_read(np->damp,NP_FMT(np),np->iter.b,"damp",argc,argv);
  if (ReadArgvDOUBLE("alpha",&(np->alpha),argc,argv)) np->alpha = 1.5;
  if (ReadArgvDOUBLE("Gamma",&(np->Gamma),argc,argv)) np->Gamma = 1.0;

  return (NPIterInit(&np->iter,argc,argv));
}

static INT SOR_A_Display (NP_BASE *theNP)
{
  NP_SOR_A *np=(NP_SOR_A *)theNP;

  NPIterDisplay(&np->iter);
  UserWrite("configuration parameters:\n");
  if (sc_disp(np->damp,np->iter.b,"damp")) REP_ERR_RETURN (1);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"alpha",np->alpha);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"Gamma",np->Gamma);

  return (0);
}

static INT SOR_A_PreProcess  (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_SOR_A *np=(NP_SOR_A *)theNP;

  if (l_setindex(NP_GRID(theNP,level))) NP_RETURN(1,result[0]);
  np->B=NULL;
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->B)) NP_RETURN(1,result[0]);
  if (np->alpha<0.0 || np->Gamma<0.0)
  {
    if (SOR_F_CopyMatrix(NP_MG(theNP),level,np->B,A)) NP_RETURN(1,result[0]);
  }
  else
  {
    if (SOR_A_CopyMatrix(NP_MG(theNP),level,np->B,A,np->alpha,np->Gamma)) NP_RETURN(1,result[0]);
  }
  *baselevel = level;

  return (0);
}


static INT SOR_A_Smoother (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_SOR_A *np=(NP_SOR_A *)theNP;

  if (l_lsor(NP_GRID(theNP,level),x,np->B,b,Factor_One,NULL)) NP_RETURN(1,result[0]);
  if (dscalx(NP_MG(theNP),level,level,ALL_VECTORS,x,np->damp)) NP_RETURN(1,result[0]);
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)) NP_RETURN(1,result[0]);

  return (0);
}

static INT SOR_A_PostProcess (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_SOR_A *np=(NP_SOR_A *)theNP;

  if (FreeMD(NP_MG(theNP),level,level,np->B)) REP_ERR_RETURN(1);

  return(0);
}

static INT SOR_A_Construct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init=SOR_A_Init;
  theNP->Display=SOR_A_Display;
  theNP->Execute=NPIterExecute;

  np=(NP_ITER *)theNP;
  np->PreProcess=SOR_A_PreProcess;
  np->Iter=SOR_A_Smoother;
  np->PostProcess=SOR_A_PostProcess;

  return(0);
}

/****************************************************************************/
/*D
   ilua - numproc for ILU_A smoother

   DESCRIPTION:
   This numproc executes an ILU_A (auto-damp) smoother.

   .vb
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>]
       $damp <sc double list> $alpha <double> $Gamma <double>;
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $damp~<sc~double~list> - additional damping factors for each component
   .  $alpha~<double> - alpha
   .  $Gamma~<double> - Gamma

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p];'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT ILU_A_CopyMatrix (MULTIGRID *mg, INT level, MATDATA_DESC *B, MATDATA_DESC *A, DOUBLE alpha, DOUBLE Gamma)
{
  INT Acomp,Bcomp;
  VECTOR *v;
  MATRIX *m;
  DOUBLE sum;

  if (dmatcopy(mg,level,level,ALL_VECTORS,B,A)) return(1);
  Acomp=MD_MCMP_OF_RT_CT(A,NODEVEC,NODEVEC,0);
  Bcomp=MD_MCMP_OF_RT_CT(B,NODEVEC,NODEVEC,0);
  for (v=FIRSTVECTOR(GRID_ON_LEVEL(mg,level)); v!=NULL; v=SUCCVC(v))
  {
    /* set diagonal */
    if (VECSKIP(v)) continue;
    sum=0.0;
    for (m=MNEXT(VSTART(v)); m!=NULL; m=MNEXT(m))
      if (!VECSKIP(MDEST(m)))
        sum+=ABS(MVALUE(m,Acomp)-MVALUE(MADJ(m),Acomp));
    sum*=0.25*alpha*Gamma;
    if (MVALUE(VSTART(v),Acomp)<0.0) sum=-sum;
    MVALUE(VSTART(v),Bcomp)+=sum;

    /* set lower and upper triangle */
    for (m=MNEXT(VSTART(v)); m!=NULL; m=MNEXT(m))
      if (!VECSKIP(MDEST(m)))
        if (MDESTINDEX(m)!=VINDEX(v))
        {
          MVALUE(m,Bcomp)=0.5*(1.0+alpha)*MVALUE(m,Acomp)+0.5*(1.0-alpha)*MVALUE(MADJ(m),Acomp);
        }
  }

  return(0);
}

static INT ILU_A_Init (NP_BASE *theNP, INT argc , char **argv)
{
  INT i;
  NP_ILU_A *np=(NP_ILU_A *)theNP;

  for (i=0; i<MAX_VEC_COMP; i++) np->damp[i]=1.0;
  sc_read(np->damp,NP_FMT(np),np->iter.b,"damp",argc,argv);
  if (ReadArgvDOUBLE("alpha",&(np->alpha),argc,argv)) np->alpha = 1.5;
  if (ReadArgvDOUBLE("Gamma",&(np->Gamma),argc,argv)) np->Gamma = 1.0;

  return (NPIterInit(&np->iter,argc,argv));
}

static INT ILU_A_Display (NP_BASE *theNP)
{
  NP_ILU_A *np=(NP_ILU_A *)theNP;

  NPIterDisplay(&np->iter);
  UserWrite("configuration parameters:\n");
  if (sc_disp(np->damp,np->iter.b,"damp")) REP_ERR_RETURN (1);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"alpha",np->alpha);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"Gamma",np->Gamma);

  return (0);
}

static INT ILU_A_PreProcess  (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_ILU_A *np=(NP_ILU_A *)theNP;

  if (l_setindex(NP_GRID(theNP,level))) NP_RETURN(1,result[0]);
  np->B=NULL;
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->B)) NP_RETURN(1,result[0]);
  if (ILU_A_CopyMatrix(NP_MG(theNP),level,np->B,A,np->alpha,np->Gamma)) NP_RETURN(1,result[0]);
  if (l_ilubthdecomp(GRID_ON_LEVEL(NP_MG(theNP),level),np->B,0,NULL,NULL,NULL)!=NUM_OK)
  {
    PrintErrorMessage('E',"ILUAPreProcess","decomposition failed");
    NP_RETURN(1,result[0]);
  }
  *baselevel = level;

  return (0);
}


static INT ILU_A_Smoother (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_ILU_A *np=(NP_ILU_A *)theNP;

  if (l_luiter(NP_GRID(theNP,level),x,np->B,b)!=NUM_OK) NP_RETURN(1,result[0]);
  if (dscalx(NP_MG(theNP),level,level,ALL_VECTORS,x,np->damp)) NP_RETURN(1,result[0]);
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)) NP_RETURN(1,result[0]);

  return (0);
}

static INT ILU_A_PostProcess (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_ILU_A *np=(NP_ILU_A *)theNP;

  if (FreeMD(NP_MG(theNP),level,level,np->B)) REP_ERR_RETURN(1);

  return(0);
}

static INT ILU_A_Construct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init=ILU_A_Init;
  theNP->Display=ILU_A_Display;
  theNP->Execute=NPIterExecute;

  np=(NP_ITER *)theNP;
  np->PreProcess=ILU_A_PreProcess;
  np->Iter=ILU_A_Smoother;
  np->PostProcess=ILU_A_PostProcess;

  return(0);
}

/****************************************************************************/
/*
   InitIter_2	- Init this file

   SYNOPSIS:
   INT InitIter_2 ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function inits this file.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
 */
/****************************************************************************/

INT InitIter_2 ()
{
  INT i;

  for (i=0; i<MAX_VEC_COMP; i++) Factor_One[i] = 1.0;

  if (CreateClass(ITER_CLASS_NAME ".sora",sizeof(NP_SOR_A),SOR_A_Construct)) REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".ilua",sizeof(NP_ILU_A),ILU_A_Construct)) REP_ERR_RETURN (__LINE__);

  return (0);
}
