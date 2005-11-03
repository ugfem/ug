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
#include "blocking.h"

USING_UG_NAMESPACES

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
  INT reg;
  MATDATA_DESC *B;
} NP_SOR_A;

typedef struct
{
  NP_ITER iter;

  VEC_SCALAR damp;
  DOUBLE alpha,Gamma;
  INT reg;
  MATDATA_DESC *B;
} NP_SSOR_A;

typedef struct
{
  NP_ITER iter;

  VEC_SCALAR damp;
  DOUBLE alpha,Gamma;
  INT reg;
  MATDATA_DESC *B;
} NP_ILU_A;

#define OBGS_MODE_NOT_INIT        0
#define OBGS_JAC                  1
#define OBGS_GS                   2
#define OBGS_SGS                  3

typedef struct
{
  NP_ITER iter;

  VEC_SCALAR damp,omega;
  NP_BLOCKING *blocking;
  INT mode,optimizeBand,gnu;

  DOUBLE *Vec[MAXLEVEL],**Mat[MAXLEVEL];
  INT MarkKey[MAXLEVEL];
  INT *bw[MAXLEVEL],*n_Mat[MAXLEVEL];
  BLOCKING_STRUCTUR bs[MAXLEVEL];
} NP_OBGS;

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

static INT AutoDamp_CopyMatrix (MULTIGRID *mg, INT level, MATDATA_DESC *B, MATDATA_DESC *A, DOUBLE alpha, DOUBLE Gamma, INT regularize)
{
  INT i,j,Acomp0,Bcomp0,Acomp,Bcomp,nc,nr;
  DOUBLE scale;
  VECTOR *v;
  MATRIX *m;

  if (dmatcopy(mg,level,level,ALL_VECTORS,B,A)) return(1);
  nc = MD_COLS_IN_RT_CT(A,NODEVEC,NODEVEC);
  nr = MD_ROWS_IN_RT_CT(A,NODEVEC,NODEVEC);
  assert(nc==nr);
  Acomp0=MD_MCMP_OF_RT_CT(A,NODEVEC,NODEVEC,0);
  Bcomp0=MD_MCMP_OF_RT_CT(B,NODEVEC,NODEVEC,0);

  for (v=FIRSTVECTOR(GRID_ON_LEVEL(mg,level)); v!=NULL; v=SUCCVC(v))
  {
    for (i=0; i<nc; i++)
    {
      if (VECSKIPBIT(v,i)) continue;

      /* components of i-th diagonal */
      Acomp=Acomp0+i*(nc+1);
      Bcomp=Bcomp0+i*(nc+1);

      /* set off-diagonals */
      for (m=MNEXT(VSTART(v)); m!=NULL; m=MNEXT(m))
        if (MDESTINDEX(m)!=VINDEX(v))
          if (!VECSKIPBIT(MDEST(m),i))
            MVALUE(m,Bcomp)=0.5*(1.0+alpha)*MVALUE(m,Acomp)+0.5*(1.0-alpha)*MVALUE(MADJ(m),Acomp);

      /* scale scalar sub-problems */
      scale=0.0;
      for (m=MNEXT(VSTART(v)); m!=NULL; m=MNEXT(m))
        if (!VECSKIP(MDEST(m)))
          scale+=ABS(MVALUE(m,Acomp)-MVALUE(MADJ(m),Acomp));
      scale*=0.25*alpha*Gamma;
      scale=1+scale/ABS(MVALUE(VSTART(v),Acomp));
      for (j=0; j<nc; j++)
        MVALUE(VSTART(v),Bcomp0+j+i*nc)*=scale;
    }

    /* set scale according to stability of diagonal block */
    if (regularize)
      switch (nc)
      {
      case 1 :
        /* nothing to do */
        break;
      case 2 :
        scale=ABS(MVALUE(VSTART(v),Bcomp0)*MVALUE(VSTART(v),Bcomp0+3)-MVALUE(VSTART(v),Bcomp0+1)*MVALUE(VSTART(v),Bcomp0+2));
        assert(scale!=0.0);
        scale=(ABS(MVALUE(VSTART(v),Bcomp0)*MVALUE(VSTART(v),Bcomp0+3))+ABS(MVALUE(VSTART(v),Bcomp0+1)*MVALUE(VSTART(v),Bcomp0+2)))/scale;
        for (i=0; i<4; i++) MVALUE(VSTART(v),Bcomp0+i)*=scale;
        break;
      }
  }

  return(0);
}

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

static INT SOR_A_Init (NP_BASE *theNP, INT argc , char **argv)
{
  INT i;
  NP_SOR_A *np=(NP_SOR_A *)theNP;

  for (i=0; i<MAX_VEC_COMP; i++) np->damp[i]=1.0;
  sc_read(np->damp,NP_FMT(np),np->iter.b,"damp",argc,argv);
  if (ReadArgvDOUBLE("alpha",&(np->alpha),argc,argv)) np->alpha = 1.5;
  if (ReadArgvDOUBLE("Gamma",&(np->Gamma),argc,argv)) np->Gamma = 1.0;
  if (np->Gamma<0.0) REP_ERR_RETURN (1);
  if (ReadArgvINT("reg",&(np->reg),argc,argv)) np->reg=1;

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
  UserWriteF(DISPLAY_NP_FORMAT_SI,"reg",(int)np->reg);

  return (0);
}

static INT SOR_A_PreProcess  (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_SOR_A *np=(NP_SOR_A *)theNP;

  if (l_setindex(NP_GRID(theNP,level))) NP_RETURN(1,result[0]);
  np->B=NULL;
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->B)) NP_RETURN(1,result[0]);
  if (AutoDamp_CopyMatrix(NP_MG(theNP),level,np->B,A,np->alpha,np->Gamma,np->reg)) NP_RETURN(1,result[0]);
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
   ssora - numproc for SSOR_A smoother

   DESCRIPTION:
   This numproc executes an SSOR_A (auto-damp) smoother.

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

static INT SSOR_A_Init (NP_BASE *theNP, INT argc , char **argv)
{
  INT i;
  NP_SSOR_A *np=(NP_SSOR_A *)theNP;

  for (i=0; i<MAX_VEC_COMP; i++) np->damp[i]=1.0;
  sc_read(np->damp,NP_FMT(np),np->iter.b,"damp",argc,argv);
  if (ReadArgvDOUBLE("alpha",&(np->alpha),argc,argv)) np->alpha = 1.5;
  if (ReadArgvDOUBLE("Gamma",&(np->Gamma),argc,argv)) np->Gamma = 1.0;
  if (np->Gamma<0.0) REP_ERR_RETURN (1);
  if (ReadArgvINT("reg",&(np->reg),argc,argv)) np->reg=1;

  return (NPIterInit(&np->iter,argc,argv));
}

static INT SSOR_A_Display (NP_BASE *theNP)
{
  NP_SSOR_A *np=(NP_SSOR_A *)theNP;

  NPIterDisplay(&np->iter);
  UserWrite("configuration parameters:\n");
  if (sc_disp(np->damp,np->iter.b,"damp")) REP_ERR_RETURN (1);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"alpha",np->alpha);
  UserWriteF(DISPLAY_NP_FORMAT_SF,"Gamma",np->Gamma);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"reg",(int)np->reg);

  return (0);
}

static INT SSOR_A_PreProcess  (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_SSOR_A *np=(NP_SSOR_A *)theNP;

  if (l_setindex(NP_GRID(theNP,level))) NP_RETURN(1,result[0]);
  np->B=NULL;
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->B)) NP_RETURN(1,result[0]);
  if (AutoDamp_CopyMatrix(NP_MG(theNP),level,np->B,A,np->alpha,np->Gamma,np->reg)) NP_RETURN(1,result[0]);
  *baselevel = level;

  return (0);
}


static INT SSOR_A_Smoother (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_SSOR_A *np=(NP_SSOR_A *)theNP;
  VECDATA_DESC *vtmp=NULL;

  if (AllocVDFromVD(NP_MG(theNP),level,level,x,&vtmp)) NP_RETURN(1,result[0]);

  /* foreward iteration */
  if (l_lsor(NP_GRID(theNP,level),x,np->B,b,Factor_One,NULL)) NP_RETURN(1,result[0]);
  if (dscalx(NP_MG(theNP),level,level,ALL_VECTORS,x,np->damp)) NP_RETURN(1,result[0]);
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)) NP_RETURN(1,result[0]);

  /* backward iteration */
  if (l_usor(NP_GRID(theNP,level),vtmp,np->B,b,Factor_One,NULL)) NP_RETURN(1,result[0]);
  if (dscalx(NP_MG(theNP),level,level,ALL_VECTORS,vtmp,np->damp)) NP_RETURN(1,result[0]);
  if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,vtmp)) NP_RETURN(1,result[0]);

  /* add correction */
  if (dadd(NP_MG(theNP),level,level,ALL_VECTORS,x,vtmp) != NUM_OK) NP_RETURN(1,result[0]);

  if (FreeVD(NP_MG(theNP),level,level,vtmp)) NP_RETURN(1,result[0]);

  return (0);
}

static INT SSOR_A_PostProcess (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_SSOR_A *np=(NP_SSOR_A *)theNP;

  if (FreeMD(NP_MG(theNP),level,level,np->B)) REP_ERR_RETURN(1);

  return(0);
}

static INT SSOR_A_Construct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init=SSOR_A_Init;
  theNP->Display=SSOR_A_Display;
  theNP->Execute=NPIterExecute;

  np=(NP_ITER *)theNP;
  np->PreProcess=SSOR_A_PreProcess;
  np->Iter=SSOR_A_Smoother;
  np->PostProcess=SSOR_A_PostProcess;

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

static INT ILU_A_Init (NP_BASE *theNP, INT argc , char **argv)
{
  INT i;
  NP_ILU_A *np=(NP_ILU_A *)theNP;

  for (i=0; i<MAX_VEC_COMP; i++) np->damp[i]=1.0;
  sc_read(np->damp,NP_FMT(np),np->iter.b,"damp",argc,argv);
  if (ReadArgvDOUBLE("alpha",&(np->alpha),argc,argv)) np->alpha = 1.5;
  if (ReadArgvDOUBLE("Gamma",&(np->Gamma),argc,argv)) np->Gamma = 1.0;
  if (ReadArgvINT("reg",&(np->reg),argc,argv)) np->reg=1;

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
  UserWriteF(DISPLAY_NP_FORMAT_SI,"reg",(int)np->reg);

  return (0);
}

static INT ILU_A_PreProcess  (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  NP_ILU_A *np=(NP_ILU_A *)theNP;

  if (l_setindex(NP_GRID(theNP,level))) NP_RETURN(1,result[0]);
  np->B=NULL;
  if (AllocMDFromMD(NP_MG(theNP),level,level,A,&np->B)) NP_RETURN(1,result[0]);
  if (AutoDamp_CopyMatrix(NP_MG(theNP),level,np->B,A,np->alpha,np->Gamma,np->reg)) NP_RETURN(1,result[0]);
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
/*D
   obgs - numproc for OBGS smoother

   DESCRIPTION:
   This numproc executes an overlapping block Gauss-Seidel (type) smoother.

   .vb
   npinit <name> [$c <cor>] [$b <rhs>] [$A <mat>]
       $damp <sc double list> $omega <double>;
   .ve

   .  $c~<cor> - correction vector
   .  $b~<rhs> - right hand side vector
   .  $A~<mat> - stiffness matrix
   .  $damp~<sc~double~list> - additional damping factors for each component
   .  $omega~<double> - sor/ssor-local damping for each component

   .  <sc~double~list>  - [nd <double  list>] | [ed <double  list>] | [el <double  list>] | [si <double  list>]
   .  <double~list>  - <double> {: <double>}*
   .n     nd = nodedata, ed = edgedata, el =  elemdata, si = sidedata

   'npexecute <name> [$i] [$s] [$p];'

   .  $i - preprocess
   .  $s - smooth
   .  $p - postprocess
   D*/
/****************************************************************************/

static INT OBGS_Init (NP_BASE *theNP, INT argc , char **argv)
{
  INT i;
  NP_OBGS *np=(NP_OBGS *)theNP;
  char buffer[128];

  for (i=0; i<MAX_VEC_COMP; i++) np->damp[i]=1.0;
  sc_read(np->damp,NP_FMT(np),np->iter.b,"damp",argc,argv);
  for (i=0; i<MAX_VEC_COMP; i++) np->omega[i]=1.0;
  sc_read(np->omega,NP_FMT(np),np->iter.b,"omega",argc,argv);
  np->blocking=(NP_BLOCKING*)ReadArgvNumProc(theNP->mg,"B",BLOCKING_CLASS_NAME,argc,argv);
  if (np->blocking==NULL) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvChar("mode",buffer,argc,argv)) strcpy(buffer,"gs");
  np->mode=OBGS_MODE_NOT_INIT;
  if (strcmp(buffer,"jac")==0) np->mode=OBGS_JAC;
  if (strcmp(buffer,"gs")==0) np->mode=OBGS_GS;
  if (strcmp(buffer,"sgs")==0) np->mode=OBGS_SGS;
  if (np->mode==OBGS_MODE_NOT_INIT) REP_ERR_RETURN (NP_NOT_ACTIVE);
  if (ReadArgvINT ("o",&np->optimizeBand,argc,argv)) np->optimizeBand=1;
  if (ReadArgvINT ("gnu",&np->gnu,argc,argv)) np->gnu=0;

  return (NPIterInit(&np->iter,argc,argv));
}

static INT OBGS_Display (NP_BASE *theNP)
{
  NP_OBGS *np=(NP_OBGS *)theNP;

  NPIterDisplay(&np->iter);
  UserWrite("configuration parameters:\n");
  if (sc_disp(np->damp,np->iter.b,"damp")) REP_ERR_RETURN (1);
  if (sc_disp(np->damp,np->iter.b,"omega")) REP_ERR_RETURN (1);
  if (np->blocking!=NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"B",ENVITEM_NAME(np->blocking));
  else UserWriteF(DISPLAY_NP_FORMAT_SS,"B","---");
  if (np->mode==OBGS_MODE_NOT_INIT) UserWriteF(DISPLAY_NP_FORMAT_SS,"mode","---");
  if (np->mode==OBGS_JAC) UserWriteF(DISPLAY_NP_FORMAT_SS,"mode","jac");
  if (np->mode==OBGS_GS) UserWriteF(DISPLAY_NP_FORMAT_SS,"mode","gs");
  if (np->mode==OBGS_SGS) UserWriteF(DISPLAY_NP_FORMAT_SS,"mode","sgs");
  UserWriteF(DISPLAY_NP_FORMAT_SS,"o",BOOL_2_YN(np->optimizeBand));
  UserWriteF(DISPLAY_NP_FORMAT_SS,"gnu",BOOL_2_YN(np->gnu));

  return (0);
}

static HEAP *OBGS_Heap;
static INT OBGS_Key;
static void *OBGS_GetMem (MEM n) {
  return (GetTmpMem(OBGS_Heap,n,OBGS_Key));
}

static INT OBGS_PrintMatrix (DOUBLE *Mat, INT n, INT bw)
{
  INT i,j;
  DOUBLE entry;

  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
    {
      if (ABS(i-j)<=bw) entry=EX_MAT(Mat,bw,i,j);
      else entry=0.0;
      printf("%15.8e ",entry);
    }
    printf("\n");
  }
  printf("\n");

  return(0);
}

static INT OBGS_PreProcess  (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *baselevel, INT *result)
{
  INT i,j,n,bl,max,max_comp_in_type,k,index,bw,rindex,rtype,rcomp,cindex,ctype,ccomp,n_Mat,n_Mat_max,MarkKey,fifo_buffer_size,LocalMarkKey;
  SHORT *comp;
  NP_OBGS *np=(NP_OBGS *)theNP;
  GRID *theGrid=NP_GRID(theNP,level);
  HEAP *theHeap=MGHEAP(NP_MG(theNP));
  VECTOR *theV,*theW,**bv,**bv_local;
  MATRIX *theM;
  DOUBLE *Mat;
  FIFO fifo;
  void *fifo_buffer;
  FILE *file;
  DOUBLE_VECTOR pos;
  char name[16];

  /* memory management */
  if (MarkTmpMem(theHeap,&(np->MarkKey[level]))) REP_ERR_RETURN(1);
  OBGS_Heap=theHeap; OBGS_Key=MarkKey=np->MarkKey[level];

  /* preprocess blocking */
  if (np->blocking==NULL) REP_ERR_RETURN(1);
  if (np->blocking->PreProcess!=NULL)
    if ((*np->blocking->PreProcess)(np->blocking,level,result)) REP_ERR_RETURN(1);

  /* perform blocking */
  if (np->blocking->Blocking!=NULL)
    if ((*np->blocking->Blocking)(np->blocking,OBGS_GetMem,level,A,&(np->bs[level]),result)) REP_ERR_RETURN(1);
  assert(np->bs[level].n!=0);

  /* allocate vectors and matrices */
  np->bw[level]=(INT *)GetTmpMem(theHeap,sizeof(INT)*np->bs[level].n,MarkKey);
  np->n_Mat[level]=(INT *)GetTmpMem(theHeap,sizeof(INT)*np->bs[level].n,MarkKey);
  np->Mat[level]=(DOUBLE **)GetTmpMem(theHeap,sizeof(DOUBLE *)*np->bs[level].n,MarkKey);
  for (theV=FIRSTVECTOR(theGrid),max_comp_in_type=0; theV!=NULL; theV=SUCCVC(theV)) { SETVCUSED(theV,0); SETVACTIVE(theV,0); k=VD_NCMPS_IN_TYPE(x,VTYPE(theV)); max_comp_in_type=MAX(max_comp_in_type,k); }
  if (np->optimizeBand)
  {
    max=0;
    for (bl=0; bl<np->bs[level].n; bl++)
      max=MAX(max,np->bs[level].nb[bl]);
    fifo_buffer_size=max*sizeof(VECTOR*);
    fifo_buffer=(void *)GetTmpMem(theHeap,fifo_buffer_size,MarkKey); assert(fifo_buffer!=NULL);
  }
  for (bl=n_Mat_max=0; bl<np->bs[level].n; bl++)
  {
    /* admin */
    n=np->bs[level].nb[bl]; bv=np->bs[level].vb[bl];
    for (i=0; i<n; i++) SETVCUSED(bv[i],1);

    /* print gnufile iff */
    if (np->gnu)
    {
      sprintf(name,"gnu_%d_%d",level,bl);
      file=fopen(name,"w");
      if (file==NULL) REP_ERR_RETURN(1);
      for (i=0; i<n; i++)
      {
        VectorPosition(bv[i],pos);
                                #ifdef __TWODIM__
        fprintf(file,"%e %e\n",pos[0],pos[1]);
                                #else
        fprintf(file,"%e %e %e\n",pos[0],pos[1],pos[2]);
                                #endif
      }
      fclose(file);
    }

    /* optimize band */
    if (np->optimizeBand)
    {
      /* allocate temporary bv-list */
      if (MarkTmpMem(theHeap,&LocalMarkKey)) REP_ERR_RETURN(1);
      bv_local=(VECTOR **)GetTmpMem(theHeap,sizeof(VECTOR *)*n,LocalMarkKey);
      for (i=0; i<n; i++) SETVACTIVE(bv[i],0);
      for (i=0; i<n; i++) SETVCFLAG(bv[i],0);
      fifo_init(&fifo,fifo_buffer,fifo_buffer_size);
      for(i=0; i<n; )
      {
        for(j=0; j<n; j++) if (!VACTIVE(bv[j])) break;
        fifo_in(&fifo,(void *)bv[j]); SETVCFLAG(bv[j],1);
        while(!fifo_empty(&fifo))
        {
          theV=(VECTOR *)fifo_out(&fifo);
          for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
            if (!VCFLAG(MDEST(theM)) && VCUSED(MDEST(theM)))
            {
              fifo_in(&fifo,(void *)MDEST(theM)); SETVCFLAG(MDEST(theM),1);
            }
        }
        fifo_in(&fifo,(void *)theV); SETVCFLAG(theV,0);
        while(!fifo_empty(&fifo))
        {
          theV=(VECTOR *)fifo_out(&fifo);
          bv_local[i++]=theV; SETVACTIVE(theV,1);
          for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
            if (VCFLAG(MDEST(theM)) && VCUSED(MDEST(theM)))
            {
              fifo_in(&fifo,(void *)MDEST(theM)); SETVCFLAG(MDEST(theM),0);
            }
        }
      }
      for (i=0; i<n; i++) SETVACTIVE(bv[i],0);
      for (i=0; i<n; i++) bv[i]=bv_local[i];

      /* release temporary bv-list */
      ReleaseTmpMem(theHeap,LocalMarkKey);
    }

    /* get bandwidth */
    for (i=1,VINDEX(bv[0])=0; i<n; i++) VINDEX(bv[i])=VINDEX(bv[i-1])+VD_NCMPS_IN_TYPE(x,VTYPE(bv[i-1]));
    for (i=0,max=0; i<n; i++)
    {
      index=VINDEX(bv[i]);
      for (theM=MNEXT(VSTART(bv[i])); theM!=NULL; theM=MNEXT(theM))
        if (VCUSED(MDEST(theM)))
        {
          k=index-MDESTINDEX(theM);
          k=ABS(k);
          max=MAX(max,k);
        }
    }
    bw=max+max_comp_in_type-1; n_Mat=VINDEX(bv[n-1])+VD_NCMPS_IN_TYPE(x,VTYPE(bv[n-1])); n_Mat_max=MAX(n_Mat_max,n_Mat);

    /* get block matrices */
    Mat=(DOUBLE*)GetTmpMem(theHeap,n_Mat*(2*bw+1)*sizeof(DOUBLE),MarkKey); assert(Mat!=NULL);
    for (i=0; i<n_Mat*(2*bw+1); i++) Mat[i]=0.0;

    /* set block matrix */
    for (i=0; i<n; i++)
    {
      rindex = VINDEX(bv[i]);
      rtype = VTYPE(bv[i]);
      rcomp = VD_NCMPS_IN_TYPE(x,rtype);
      for (theM=VSTART(bv[i]); theM!=NULL; theM=MNEXT(theM))
      {
        theW=MDEST(theM);
        if (!VCUSED(theW)) continue;
        cindex = VINDEX(theW);
        ctype = VTYPE(theW);
        ccomp = VD_NCMPS_IN_TYPE(x,ctype);
        comp = MD_MCMPPTR_OF_RT_CT(A,rtype,ctype);
        for (j=0; j<rcomp; j++)
          for (k=0; k<ccomp; k++)
            EX_MAT(Mat,bw,rindex+j,cindex+k) = MVALUE(theM,comp[j*ccomp+k]);
      }
    }

    /* decompose matrix */
    if (EXDecomposeMatrixDOUBLE(Mat,bw,n_Mat)) REP_ERR_RETURN(1);

    /* admin */
    np->bw[level][bl]=bw; np->Mat[level][bl]=Mat; np->n_Mat[level][bl]=n_Mat;
    for (i=0; i<n; i++) SETVCUSED(bv[i],0);
  }

  /* allocate help vector */
  np->Vec[level]=(DOUBLE *)GetTmpMem(theHeap,sizeof(DOUBLE)*n_Mat_max,MarkKey); assert(np->Vec[level]!=NULL);

  *baselevel = level;

  return (0);
}

static INT OBGS_Smoother (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_OBGS *np=(NP_OBGS *)theNP;
  VECTOR **bv,*theW;
  SHORT *comp,*rcomp,*ccomp;
  MATRIX *theM;
  INT i,j,k,n,bl,bw,type,rtype,ctype,nccomp,nrcomp,n_Mat;
  DOUBLE *Mat;
  DOUBLE *vec=np->Vec[level];

  /* reset correction */
  if (dset(NP_MG(theNP),level,level,ALL_VECTORS,x,0.0)) REP_ERR_RETURN(1);

  /* iterate */
  if (np->mode==OBGS_JAC)
  {
    for (bl=0; bl<np->bs[level].n; bl++)
    {
      /* admin */
      n=np->bs[level].nb[bl]; bv=np->bs[level].vb[bl]; Mat=np->Mat[level][bl]; bw=np->bw[level][bl]; n_Mat=np->n_Mat[level][bl];

      /* copy b to vec */
      for (i=j=0; i<n; i++)
      {
        type=VTYPE(bv[i]);
        comp=VD_CMPPTR_OF_TYPE(b,type);
        VINDEX(bv[i])=j;
        for (k=0; k<VD_NCMPS_IN_TYPE(b,type); k++)
          vec[j++]=VVALUE(bv[i],comp[k]);
      }

      /* apply inverse */
      if (EXApplyLUDOUBLE(Mat,bw,n_Mat,vec)) REP_ERR_RETURN(1);

      /* add correction to x */
      for (i=j=0; i<n; i++)
      {
        type=VTYPE(bv[i]);
        comp=VD_CMPPTR_OF_TYPE(x,type);
        for (k=0; k<VD_NCMPS_IN_TYPE(b,type); k++)
        {
          vec[j]*=np->damp[k];
          VVALUE(bv[i],comp[k])+=vec[j++];
        }
      }
    }

    /* update defect */
    if (dscalx(NP_MG(theNP),level,level,ALL_VECTORS,x,np->damp)) NP_RETURN(1,result[0]);
    if (dmatmul_minus(NP_MG(theNP),level,level,ALL_VECTORS,b,A,x)) NP_RETURN(1,result[0]);
  }
  if (np->mode==OBGS_GS || np->mode==OBGS_SGS)
  {
    for (bl=0; bl<np->bs[level].n; bl++)
    {
      /* admin */
      n=np->bs[level].nb[bl]; bv=np->bs[level].vb[bl]; Mat=np->Mat[level][bl]; bw=np->bw[level][bl]; n_Mat=np->n_Mat[level][bl];

      /* copy b to vec */
      for (i=j=0; i<n; i++)
      {
        type=VTYPE(bv[i]);
        comp=VD_CMPPTR_OF_TYPE(b,type);
        VINDEX(bv[i])=j;
        for (k=0; k<VD_NCMPS_IN_TYPE(b,type); k++)
          vec[j++]=VVALUE(bv[i],comp[k]);
      }

      /* apply inverse */
      if (EXApplyLUDOUBLE(Mat,bw,n_Mat,vec)) REP_ERR_RETURN(1);

      /* add correction to x */
      for (i=j=0; i<n; i++)
      {
        type=VTYPE(bv[i]);
        comp=VD_CMPPTR_OF_TYPE(x,type);
        for (k=0; k<VD_NCMPS_IN_TYPE(b,type); k++)
        {
          vec[j]*=np->omega[k];
          VVALUE(bv[i],comp[k])+=vec[j++];
        }
      }

      /* update defect */
      for (i=0; i<n; i++)
      {
        ctype=VTYPE(bv[i]);
        ccomp=VD_CMPPTR_OF_TYPE(x,ctype);
        nccomp=VD_NCMPS_IN_TYPE(x,ctype);
        for (theM=VSTART(bv[i]); theM!=NULL; theM=MNEXT(theM))
        {
          theW=MDEST(theM);
          rtype = VTYPE(theW);
          rcomp = VD_CMPPTR_OF_TYPE(b,rtype);
          nrcomp=VD_NCMPS_IN_TYPE(b,rtype);
          comp = MD_MCMPPTR_OF_RT_CT(A,rtype,ctype);
          for (j=0; j<nccomp; j++)
            for (k=0; k<nrcomp; k++)
              VVALUE(theW,rcomp[k])-=MVALUE(MADJ(theM),comp[k*nccomp+j])*vec[VINDEX(bv[i])+j];
        }
      }
    }
  }
  if (np->mode==OBGS_SGS)
  {
    for (bl=np->bs[level].n-1; bl>=0; bl--)
    {
      /* admin */
      n=np->bs[level].nb[bl]; bv=np->bs[level].vb[bl]; Mat=np->Mat[level][bl]; bw=np->bw[level][bl]; n_Mat=np->n_Mat[level][bl];

      /* copy b to vec */
      for (i=j=0; i<n; i++)
      {
        type=VTYPE(bv[i]);
        comp=VD_CMPPTR_OF_TYPE(b,type);
        VINDEX(bv[i])=j;
        for (k=0; k<VD_NCMPS_IN_TYPE(b,type); k++)
          vec[j++]=VVALUE(bv[i],comp[k]);
      }

      /* apply inverse */
      if (EXApplyLUDOUBLE(Mat,bw,n_Mat,vec)) REP_ERR_RETURN(1);

      /* add correction to x */
      for (i=j=0; i<n; i++)
      {
        type=VTYPE(bv[i]);
        comp=VD_CMPPTR_OF_TYPE(x,type);
        for (k=0; k<VD_NCMPS_IN_TYPE(b,type); k++)
        {
          vec[j]*=np->omega[k];
          VVALUE(bv[i],comp[k])+=vec[j++];
        }
      }

      /* update defect */
      for (i=0; i<n; i++)
      {
        ctype=VTYPE(bv[i]);
        ccomp=VD_CMPPTR_OF_TYPE(x,ctype);
        nccomp=VD_NCMPS_IN_TYPE(x,ctype);
        for (theM=VSTART(bv[i]); theM!=NULL; theM=MNEXT(theM))
        {
          theW=MDEST(theM);
          rtype = VTYPE(theW);
          rcomp = VD_CMPPTR_OF_TYPE(b,rtype);
          nrcomp=VD_NCMPS_IN_TYPE(b,rtype);
          comp = MD_MCMPPTR_OF_RT_CT(A,rtype,ctype);
          for (j=0; j<nccomp; j++)
            for (k=0; k<nrcomp; k++)
              VVALUE(theW,rcomp[k])-=MVALUE(MADJ(theM),comp[k*nccomp+j])*vec[VINDEX(bv[i])+j];
        }
      }
    }
  }

  return (0);
}

static INT OBGS_PostProcess (NP_ITER *theNP, INT level, VECDATA_DESC *x, VECDATA_DESC *b, MATDATA_DESC *A, INT *result)
{
  NP_OBGS *np=(NP_OBGS *)theNP;
  HEAP *theHeap=MGHEAP(NP_MG(theNP));

  if (np->blocking->PostProcess!=NULL)
    if ((*np->blocking->PostProcess)(np->blocking,level,result)) REP_ERR_RETURN(1);

  /* memory management */
  ReleaseTmpMem(theHeap,np->MarkKey[level]);

  return(0);
}

static INT OBGS_Construct (NP_BASE *theNP)
{
  NP_ITER *np;

  theNP->Init=OBGS_Init;
  theNP->Display=OBGS_Display;
  theNP->Execute=NPIterExecute;

  np=(NP_ITER *)theNP;
  np->PreProcess=OBGS_PreProcess;
  np->Iter=OBGS_Smoother;
  np->PostProcess=OBGS_PostProcess;

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

INT NS_DIM_PREFIX InitIter_2 ()
{
  INT i;

  for (i=0; i<MAX_VEC_COMP; i++) Factor_One[i] = 1.0;

  if (CreateClass(ITER_CLASS_NAME ".sora",sizeof(NP_SOR_A),SOR_A_Construct)) REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".ssora",sizeof(NP_SSOR_A),SSOR_A_Construct)) REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".ilua",sizeof(NP_ILU_A),ILU_A_Construct)) REP_ERR_RETURN (__LINE__);
  if (CreateClass(ITER_CLASS_NAME ".obgs",sizeof(NP_OBGS),OBGS_Construct)) REP_ERR_RETURN (__LINE__);

  return (0);
}
