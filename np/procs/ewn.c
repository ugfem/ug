// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ew.c	                                                                                                        */
/*																			*/
/* Purpose:   eigenvalue solver num procs                                       */
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart			                                                                */
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   Januar 7, 1997                                                                            */
/*																			*/
/* Remarks:                                                                                                                             */
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

#include "general.h"

#include "debug.h"
#include "ugstruct.h"
#include "ugdevices.h"
#include "debug.h"
#include "gm.h"
#include "ugblas.h"
#include "disctools.h"
#include "scan.h"
#include "numproc.h"
#include "pcr.h"
#include "formats.h"
#include "np.h"
#include "ugstruct.h"
#include "block.h"

#include "assemble.h"
#include "transfer.h"
#include "ls.h"

#include "ew.h"

#include "quadrature.h"
#include "shapes.h"
#include "evm.h"

#include "project.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define ABS_LIMIT 1e-10
#define VERY_SMALL 1e-10

#define MD_UNDEF        0
#define MD_STD          1
#define MD_0            2

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct
{
  NP_EW_SOLVER ew;

  /* numproc references */
  NP_LINEAR_SOLVER *LS;
  NP_TRANSFER *Transfer;
  NP_PROJECT *Project;

  /* parameters */
  INT maxiter;
  INT baselevel;
  INT display;
  INT reset;
  INT c_n;                                                               /* convergence on c_n vectors      */
  INT c_d;                                                               /* up to c_d digits                */
  INT type;                                  /* usage of mass matrix            */
  DOUBLE shift;                                                  /* shift of center of inverse iter */
  DOUBLE scale;                                                  /* scaling factor for eigenvalues  */

  /* dynamical data */
  VECDATA_DESC *r;                           /* help vector                     */
  VECDATA_DESC *t;                           /* help vector                     */
  VECDATA_DESC *E;                           /* mass matrix, lumped, modified   */
  MATDATA_DESC *M;                           /* matrix                          */
  VECDATA_DESC *e[MAX_NUMBER_EW];            /* eigenvectors                    */

} NP_EWN;

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

static INT SetUnsymmetric (MULTIGRID *mg, INT fl, INT tl, const VECDATA_DESC *x, INT xclass, INT index)
{
  SHORT i,j;
  INT vtype;
  INT lev;

  for (lev=fl; lev<=tl; lev++)
    l_setindex(GRID_ON_LEVEL(mg,lev));
  j = 0;
  for (vtype=0; vtype<NVECTYPES; vtype++)
    if (VD_ISDEF_IN_TYPE(x,vtype))
    {
      SHORT ncomp = VD_NCMPS_IN_TYPE(x,vtype);
      VECTOR *v;
      DOUBLE_VECTOR pos;

      A_VLOOP__TYPE_CLASS(lev,fl,tl,v,mg,vtype,xclass)
      {
        for (i=0; i<ncomp; i++)
          VVALUE(v,VD_CMP_OF_TYPE(x,vtype,i)) = 0.0;
        if (VECSKIP(v) != 0) continue;
        if (j++ < index) continue;
        if (VINDEX(v) % (index+2) == 0) continue;
        VectorPosition(v,pos);
        for (i=0; i<ncomp; i++)
          VVALUE(v,VD_CMP_OF_TYPE(x,vtype,i)) = pos[i]+1.0/(1.0+index*VINDEX(v)*VINDEX(v));
      }
    }
        #ifdef ModelP
  if (a_vector_consistent(mg,fl,tl,x)) return(NUM_ERROR);
    #endif

  return (NUM_OK);
}

static INT dm0add2 (MULTIGRID *mg, INT fl, INT tl, INT mode, VECDATA_DESC *x, MATDATA_DESC *A)
{
  return (1);
}

static INT dmassadd (MULTIGRID *mg, INT fl, INT tl, INT mode, const MATDATA_DESC *A, const VECDATA_DESC *x, INT type)
{
  switch(type)
  {
  case MD_0 :
    return (dm0add(mg,fl,tl,mode,x,A));
  default :
    return(1);
  }
}

static INT EWPreProcess (NP_EW_SOLVER *theNP, INT level, INT nev, VECDATA_DESC **ev, NP_NL_ASSEMBLE *Assemble, INT *result)
{
  NP_EWN *np;
  INT i,bl;

  np = (NP_EWN *) theNP;

  bl = 0;
  for (i=1; i<nev; i++)
    if (AllocVDFromVD(theNP->base.mg,bl,level,ev[0],&ev[i])) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(theNP->base.mg,bl,level,ev[0],&np->r)) NP_RETURN(1,result[0]);
  if (AllocVDFromVD(theNP->base.mg,bl,level,ev[0],&np->t)) NP_RETURN(1,result[0]);
  if (AllocMDFromVD(theNP->base.mg,bl,level,ev[0],ev[0],&np->M)) NP_RETURN(1,result[0]);
  if (np->reset)
    for (i=0; i<nev; i++)
      if (SetUnsymmetric(theNP->base.mg,bl,level,ev[i],EVERY_CLASS,i)) NP_RETURN(1,result[0]);
  if (np->shift!=0.0)
  {
    if (dcopy(theNP->base.mg,bl,level,ALL_VECTORS,np->t,np->E)) return(1);
    if (dscal(theNP->base.mg,bl,level,ALL_VECTORS,np->t,-np->shift)) return(1);
    if (dmassadd(theNP->base.mg,bl,level,ALL_VECTORS,np->M,np->t,np->type)) return(1);
  }
  np->reset = 0;
  if (np->LS->PreProcess != NULL) if ((*np->LS->PreProcess)(np->LS,level,ev[0],np->r,np->M, &np->baselevel,result)) return(1);

  return (0);
}

static INT EWPostProcess (NP_EW_SOLVER *theNP, INT level, INT nev, VECDATA_DESC **ev, NP_NL_ASSEMBLE *Assemble, INT *result)
{
  NP_EWN *np;
  INT i,bl;

  np = (NP_EWN *) theNP;
  bl = 0;

  for (i=1; i<nev; i++)
    if (FreeVD(theNP->base.mg,bl,level,ev[i])) NP_RETURN(1,result[0]);
  if (FreeVD(theNP->base.mg,bl,level,np->r)) NP_RETURN(1,result[0]);
  if (FreeVD(theNP->base.mg,bl,level,np->t)) NP_RETURN(1,result[0]);
  if (FreeMD(theNP->base.mg,bl,level,np->M)) NP_RETURN(1,result[0]);
  for (i=0; i<nev; i++)
    if ((*np->Transfer->ProjectSolution)(np->Transfer,bl,level,ev[i],result)) NP_RETURN(1,result[0]);
  if (np->shift!=0.0)
  {
    if (dcopy(theNP->base.mg,bl,level,ALL_VECTORS,np->t,np->E)) return(1);
    if (dscal(theNP->base.mg,bl,level,ALL_VECTORS,np->t,np->shift)) return(1);
    if (dmassadd(theNP->base.mg,bl,level,ALL_VECTORS,np->M,np->t,np->type)) return(1);
  }
  if (np->LS->PostProcess != NULL)
    if ((*np->LS->PostProcess)(np->LS,level,ev[0],np->r,np->M,result)) NP_RETURN(1,result[0]);

  return (0);
}

static INT EWInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_EW_SOLVER *npew;
  NP_EWN *np;
  INT i,n;
  char *token,*names,buffer[128];

  npew = (NP_EW_SOLVER *) theNP;
  np = (NP_EWN *) theNP;

  np->reset = 1;
  np->LS = (NP_LINEAR_SOLVER *)ReadArgvNumProc(theNP->mg,"L",LINEAR_SOLVER_CLASS_NAME,argc,argv);
  if (np->LS == NULL) return(NP_NOT_ACTIVE);
  np->Transfer = (NP_TRANSFER *)ReadArgvNumProc(theNP->mg,"T",TRANSFER_CLASS_NAME,argc,argv);
  np->Project = (NP_PROJECT *)ReadArgvNumProc(theNP->mg,"P",PROJECT_CLASS_NAME,argc,argv);

  np->M = ReadArgvMatDesc(theNP->mg,"M",argc,argv);
  if (np->M==NULL) return(NP_NOT_ACTIVE);
  np->E = ReadArgvVecDesc(theNP->mg,"E",argc,argv);
  if (np->E==NULL) return(NP_NOT_ACTIVE);
  np->t = ReadArgvVecDesc(theNP->mg,"t",argc,argv);
  np->r = ReadArgvVecDesc(theNP->mg,"r",argc,argv);
  if (ReadArgvINT("m",&(np->maxiter),argc,argv))
    return(NP_NOT_ACTIVE);
  np->display = ReadArgvDisplay(argc,argv);
  np->baselevel = 0;

  n = 0;
  for (i=1; i<argc; i++)
    if (argv[i][0]=='e')
    {
      if (sscanf(argv[i],"e %s",buffer)!=1)
      {
        UserWrite("Missing symbol for eigenvector in init of ew\n");
        return(NP_NOT_ACTIVE);
      }
      names=argv[i];
      names++;
      while ((*names==' ')||(*names=='\t')) names++;
      token = strtok(names," ");
      npew->ev[n] = GetVecDataDescByName(npew->base.mg,token);
      if (npew->ev[n] == NULL) npew->ev[n] = CreateVecDescOfTemplate(npew->base.mg,token,NULL);
      if (npew->ev[n++] == NULL) return(NP_NOT_ACTIVE);
      token = strtok(NULL," ");
      if (token!=NULL)
        if (sscanf(token,"%d",&n) != 1)
        {
          n = 1;
          while (token!=NULL) {
            npew->ev[n] = GetVecDataDescByName(npew->base.mg,token);
            if (npew->ev[n] == NULL)
              npew->ev[n] = CreateVecDescOfTemplate(npew->base.mg,
                                                    token,NULL);
            if (npew->ev[n++] == NULL)
              return(NP_NOT_ACTIVE);
            token = strtok(NULL," ");
          }
        }
    }
  npew->nev = n;
  if (ReadArgvINT("c_n",&np->c_n,argc,argv)) np->c_n=npew->nev;
  if (np->c_n<1 || np->c_n>npew->nev) return(NP_NOT_ACTIVE);
  if (ReadArgvINT("c_d",&np->c_d,argc,argv)) np->c_d=6;
  if (np->c_d<1 || np->c_d>16) return(NP_NOT_ACTIVE);
  if (sc_read(npew->abslimit,NP_FMT(npew),npew->ev[0],"abslimit",argc,argv))
    for (i=0; i<MAX_VEC_COMP; i++)
      npew->abslimit[i] = ABS_LIMIT;
  if (sc_read(npew->reduction,NP_FMT(npew),npew->ev[0],"red",argc,argv))
    return(NP_ACTIVE);
  if (ReadArgvChar("type",buffer,argc,argv)) return(NP_ACTIVE);
  if (strcmp(buffer,"std")==0) np->type=MD_STD;
  else if (strcmp(buffer,"0")==0) np->type=MD_0;
  else np->type=MD_UNDEF;
  if (ReadArgvDOUBLE("scale",&np->scale,argc,argv)) np->scale=1.0;
  if (np->scale<=0.0) return(NP_ACTIVE);
  if (ReadArgvDOUBLE("shift",&np->shift,argc,argv)) np->shift=0.0;
  np->shift/=np->scale;

  return(NP_EXECUTABLE);
}

static INT EWDisplay (NP_BASE *theNP)
{
  INT i;
  NP_EW_SOLVER *npew;
  NP_EWN *np;

  np = (NP_EWN *) theNP;
  npew = (NP_EW_SOLVER *) theNP;

  if (npew->nev > 0) UserWrite("symbolic user data:\n");
  for (i=0; i<npew->nev; i++)
    if (i<10) UserWriteF("ev[%d]            = %-35.32s\n", i,ENVITEM_NAME(npew->ev[i]));
    else UserWriteF("ev[%d]           = %-35.32s\n", i,ENVITEM_NAME(npew->ev[i]));
  UserWrite("\n");
  UserWrite("configuration parameters:\n");
  if (sc_disp(npew->reduction,npew->ev[0],"red")) return (1);
  if (sc_disp(npew->abslimit,npew->ev[0],"abslimit")) return (1);

  UserWriteF(DISPLAY_NP_FORMAT_SI,"m",(int)np->maxiter);
  if (np->LS != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"L",ENVITEM_NAME(np->LS));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"L","---");
  if (np->Transfer != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T",ENVITEM_NAME(np->Transfer));
  else
    UserWriteF(DISPLAY_NP_FORMAT_SS,"T","---");
  if (np->display == PCR_NO_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","NO_DISPLAY");
  else if (np->display == PCR_RED_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","RED_DISPLAY");
  else if (np->display == PCR_FULL_DISPLAY)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"DispMode","FULL_DISPLAY");
  if (np->r != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"r",ENVITEM_NAME(np->r));
  if (np->t != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"t",ENVITEM_NAME(np->t));
  if (np->M != NULL)
    UserWriteF(DISPLAY_NP_FORMAT_SS,"M",ENVITEM_NAME(np->M));

  return(0);
}

static INT EWExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_EWN *np;
  EWRESULT ewresult;
  INT i,result,level,nev,m;

  np = (NP_EWN *) theNP;
  level = CURRENTLEVEL(theNP->mg);

  nev = np->ew.nev;       /* voreingestellte Anzahl EWe sichern */
  if (!ReadArgvINT("m",&m,argc,argv))
    if (0 < m && m < nev)
      np->ew.nev = m;          /* Anzahl EWe reduzieren */
    else
      UserWriteF("EWExecute: $m %d out of range - working with maximum %d EV\n",m,nev);

  np->reset = ReadArgvOption("r",argc,argv);
  if ((*np->ew.PreProcess)(&np->ew,level,np->ew.nev,np->ew.ev,NULL,&result))
  {
    UserWriteF("EWExecute: PreProcess failed, error code %d\n",result);
    return (1);
  }
  if ((*np->ew.Solver)(&np->ew,level,np->ew.nev,np->ew.ev,np->ew.ew,NULL,np->ew.abslimit,np->ew.reduction,&ewresult))
  {
    UserWriteF("EWSolverExecute: Solver failed, error code %d\n",ewresult.error_code);
    return (1);
  }

  if ((*np->ew.PostProcess)(&np->ew,level,np->ew.nev,np->ew.ev,NULL,&result))
  {
    UserWriteF("EWExecute: PostProcess failed, error code %d\n",result);
    return (1);
  }
  np->ew.nev = nev;       /* restore old value */

  return(0);
}

static int EWCompare (DOUBLE **index1, DOUBLE **index2)
{
  DOUBLE a,b;

  a=**index1;
  b=**index2;

  if (ABS(a)>ABS(b)) return (1);
  if (ABS(a)<ABS(b)) return (-1);
  if (ABS(a)==ABS(b) && b<0.0) return (1);
  if (ABS(a)==ABS(b) && a<0.0) return (-1);
  return (0);
}

static INT SmallEWNSolver_Sci (INT nev, DOUBLE A[MAX_NUMBER_EW][MAX_NUMBER_EW], DOUBLE B[MAX_NUMBER_EW][MAX_NUMBER_EW], DOUBLE *re, DOUBLE *im, DOUBLE E[MAX_NUMBER_EW][MAX_NUMBER_EW])
{
  INT i,j;
  FILE *file;

  file=fopen("/tmp/Sci_A","w");
  for (i=0; i<nev; i++)
  {
    for (j=0; j<nev; j++)
      fprintf(file,"%16.16e ",A[i][j]);
    fprintf(file,"\n");
  }
  fclose(file);

  file=fopen("/tmp/Sci_B","w");
  for (i=0; i<nev; i++)
  {
    for (j=0; j<nev; j++)
      fprintf(file,"%16.16e ",B[i][j]);
    fprintf(file,"\n");
  }
  fclose(file);

  system("sci_spec");

  file=fopen("/tmp/Sci_S","r");
  for (i=0; i<nev; i++)
    fscanf(file,"%lf",re+i);
  for (i=0; i<nev; i++)
    fscanf(file,"%lf",im+i);
  fclose(file);

  file=fopen("/tmp/Sci_E","r");
  for (i=0; i<nev; i++)
    for (j=0; j<nev; j++)
      fscanf(file,"%lf",&E[i][j]);
  fclose(file);

  return(0);
}

static DOUBLE Round (DOUBLE v, INT n)
{
  DOUBLE s,sign;

  if (v==0.0) return(0.0);
  if (v>0.0) sign=1.0;
  else { v=-v; sign=-1.0; }

  s=pow(10.0,n+floor(-log10(v)));
  v=floor(v*s+0.5)/s;
  return(sign*v);
}

static INT dmassdot (MULTIGRID *mg, INT fl, INT tl, INT mode, const VECDATA_DESC *x, const VECDATA_DESC *y, INT type)
{
  switch(type)
  {
  case MD_STD :
    return (dpdot(mg,fl,tl,mode,x,y));
  case MD_0 :
    return (dm0dot(mg,fl,tl,mode,x,y));
  default :
    return(1);
  }
}

static INT EWNSolver (NP_EW_SOLVER *theNP, INT level, INT New, VECDATA_DESC **ev, DOUBLE *ew, NP_NL_ASSEMBLE *Assemble, VEC_SCALAR abslimit, VEC_SCALAR reduction, EWRESULT *ewresult)
{
  NP_EWN     *np    = (NP_EWN *) theNP;
  MULTIGRID *theMG = theNP->base.mg;
  INT i,j,k,l,PrintID,iter,done;
  char text[DISPLAY_WIDTH+4],format1[64],format2[64],format3[64],formatr1[64],formatr2[64];
  DOUBLE a[2],rq,s;
  DOUBLE A[MAX_NUMBER_EW][MAX_NUMBER_EW];
  DOUBLE B[MAX_NUMBER_EW][MAX_NUMBER_EW];
  DOUBLE BL[MAX_NUMBER_EW*MAX_NUMBER_EW];
  DOUBLE L[MAX_NUMBER_EW*MAX_NUMBER_EW];
  DOUBLE G[MAX_NUMBER_EW][MAX_NUMBER_EW];
  DOUBLE GL[MAX_NUMBER_EW][MAX_NUMBER_EW];
  DOUBLE E[MAX_NUMBER_EW][MAX_NUMBER_EW];
  DOUBLE ew_re[MAX_NUMBER_EW],ew_im[MAX_NUMBER_EW];
  DOUBLE tmp_re[MAX_NUMBER_EW],tmp_im[MAX_NUMBER_EW];
  DOUBLE old_re[MAX_NUMBER_EW],old_im[MAX_NUMBER_EW];
  DOUBLE* table[MAX_NUMBER_EW];
  INT index[MAX_NUMBER_EW];
  INT bl = 0;

  ewresult->error_code = 0;
  CenterInPattern(text,DISPLAY_WIDTH,ENVITEM_NAME(np),'§',"\n"); UserWrite(text);
  sprintf(text,"%d.%d",np->c_d+1,np->c_d-1);
  strcpy(format1," %-3d  %3d: (% "); strcat(format1,text); strcat(format1,"e, % "); strcat(format1,text); strcat(format1,"e)\n");
  strcpy(format2,"      %3d: (% "); strcat(format2,text); strcat(format2,"e, % "); strcat(format2,text); strcat(format2,"e)\n");
  strcpy(format3,"      %3d: [% "); strcat(format3,text); strcat(format3,"e, % "); strcat(format3,text); strcat(format3,"e]\n");
  strcpy(formatr1," %-3d   res: %3d: (% "); strcat(formatr1,text); strcat(formatr1,"e, % "); strcat(formatr1,text); strcat(formatr1,"e)\n");
  strcpy(formatr2,"            %3d: (% "); strcat(formatr2,text); strcat(formatr2,"e, % "); strcat(formatr2,text); strcat(formatr2,"e)\n");
  for (iter=0; iter<np->maxiter; iter++)
  {
    if (iter > 0)
      for (i=0; i<New; i++)
      {
        if (dcopy(theMG,bl,level,ALL_VECTORS,np->t,ev[i])) NP_RETURN(1,ewresult->error_code);
        if (dmassdot(theMG,bl,level,ALL_VECTORS,np->t,np->E,np->type)) NP_RETURN(1,ewresult->error_code);
        if ((*np->LS->Defect)(np->LS,level,ev[i],np->t,np->M, &ewresult->error_code)) NP_RETURN(1,ewresult->error_code);
        if ((*np->LS->Residuum)(np->LS,0,level,ev[i],np->t,np->M, &ewresult->lresult[i])) NP_RETURN(1,ewresult->error_code);
        if ((*np->LS->Solver)(np->LS,level,ev[i],np->t,np->M, abslimit,reduction, &ewresult->lresult[i])) NP_RETURN(1,ewresult->error_code);
        if (np->Project != NULL)
          if (np->Project->Project(np->Project,bl,level, ev[i],&ewresult->error_code) != NUM_OK) NP_RETURN(1,ewresult->error_code);
      }
    for (i=0; i<New; i++)
    {
      if (ddot(theMG,0,level,ON_SURFACE,ev[i],ev[i],&B[i][i])) NP_RETURN(1,ewresult->error_code);
      if (dscal(theMG,0,level,ALL_VECTORS,ev[i],1/sqrt(B[i][i])) != NUM_OK) NP_RETURN(1,ewresult->error_code);
    }
    if (dscal(theMG,0,level,ALL_VECTORS,np->r,1/sqrt(B[0][0])) != NUM_OK) NP_RETURN(1,ewresult->error_code);
    for (i=0; i<New; i++)
    {
      if (dmatmul (theMG,0,level,ON_SURFACE,np->t,np->M,ev[i]) != NUM_OK) NP_RETURN(1,ewresult->error_code);
      for (j=0; j<New; j++)
        if (ddot(theMG,0,level,ON_SURFACE,np->t,ev[j],&A[i][j])) NP_RETURN(1,ewresult->error_code);
    }
    for (i=0; i<New; i++)
    {
      if (dcopy(theMG,bl,level,ALL_VECTORS,np->t,ev[i])) NP_RETURN(1,ewresult->error_code);
      if (dmassdot(theMG,bl,level,ALL_VECTORS,np->t,np->E,np->type)) NP_RETURN(1,ewresult->error_code);
      for (j=0; j<New; j++)
        if (ddot(theMG,0,level,ON_SURFACE,np->t,ev[j],&B[i][j])) NP_RETURN(1,ewresult->error_code);
    }

    /* Special Eigenvalue problem  G E_i = lambda E_i */
    SmallEWNSolver_Sci(New,A,B,ew_re,ew_im,E);

    IFDEBUG(np,1)
    UserWriteF("A\n"); for (i=0; i<New; i++) { for (j=0; j<New; j++) UserWriteF("%8.4f\t",A[i][j]);UserWriteF("\n"); }
    UserWriteF("B\n"); for (i=0; i<New; i++) { for (j=0; j<New; j++) UserWriteF("%8.4f\t",B[i][j]);UserWriteF("\n"); }
    UserWriteF("E\n"); for (i=0; i<New; i++) { for (j=0; j<New; j++) UserWriteF("%8.4f\t",E[i][j]);UserWriteF("\n"); }
    ENDDEBUG

    for (i=0; i<New; i++)
      if (AllocVDFromVD(theMG,bl,level,ev[0],&np->e[i])) NP_RETURN(1,ewresult->error_code);
    for (i=0; i<New; i++)
    {
      if (dset(theMG,bl,level,ALL_VECTORS,np->e[i],0.0)) NP_RETURN(1,ewresult->error_code);
      for (j=0; j<New; j++)
        if (daxpy(theMG,bl,level,ALL_VECTORS,np->e[i],E[j][i],ev[j]) != NUM_OK) NP_RETURN(1,ewresult->error_code);
    }

    /* sort eigen-values/vectors */
    for (i=0; i<New; i++)
    {
      ew[i]=sqrt(ew_re[i]*ew_re[i]+ew_im[i]*ew_im[i]);
      if (ew_im[i]<0.0) ew[i]*=-1.0;
      tmp_re[i]=ew_re[i]; tmp_im[i]=ew_im[i];
      table[i] = &ew[i];
    }
    qsort(table, New, sizeof(*table),(int (*)(const void *, const void *))EWCompare);
    for (i=0; i<New; i++)
      for (j=0; j<New; j++)
        if (table[i]==&ew[j])
          index[i] = j;

    for (i=0; i<New; i++)
    {
      if (dcopy(theMG,bl,level,ALL_VECTORS,ev[i],np->e[index[i]])) NP_RETURN(1,ewresult->error_code);
      ew_re[i]=tmp_re[index[i]]; ew_im[i]=tmp_im[index[i]];
      if (FreeVD(theMG,bl,level,np->e[index[i]])) NP_RETURN(1,ewresult->error_code);
    }

    /* scale/shift eigenvalues */
    for (i=0; i<New; i++)
    {
      ew_re[i]=np->scale*(ew_re[i]+np->shift);
      ew_im[i]=np->scale*ew_im[i];
    }

    /* display */
    if (np->display > PCR_NO_DISPLAY)
    {
      for (i=0; i<np->ew.nev; i++)
        if (i==0)
          UserWriteF(format1,(int)iter,(int)i,ew_re[i],ew_im[i]);
        else if (i<np->c_n)
          UserWriteF(format2,(int)i,ew_re[i],ew_im[i]);
        else
          UserWriteF(format3,(int)i,ew_re[i],ew_im[i]);
      UserWriteF("\n");
    }

    /* check convergence */
    if (iter>0)
    {
      done=1;
      for (i=0; i<np->c_n; i++)
      {
        if (Round(ew_re[i],np->c_d)!=Round(old_re[i],np->c_d)) done=0;
        if (Round(ew_im[i],np->c_d)!=Round(old_im[i],np->c_d)) done=0;
      }
    }
    else
      done=0;
    for (i=0; i<New; i++)
    {
      old_re[i]=ew_re[i];
      old_im[i]=ew_im[i];
    }
    if (done) break;
  }

  /* print result */
  for (i=0; i<np->c_n; i++)
    if (i==0)
      UserWriteF(formatr1,(int)iter,(int)i,ew_re[i],ew_im[i]);
    else
      UserWriteF(formatr2,(int)i,ew_re[i],ew_im[i]);
  UserWriteF("\n");

  return (0);
}

static INT EWNConstruct (NP_BASE *theNP)
{
  NP_EWN *np;

  theNP->Init = EWInit;
  theNP->Display = EWDisplay;
  theNP->Execute = EWExecute;

  np = (NP_EWN *) theNP;
  np->ew.PreProcess = EWPreProcess;
  np->ew.Rayleigh = NULL;
  np->ew.Solver = EWNSolver;
  np->ew.PostProcess = EWPostProcess;

  return(0);
}

/****************************************************************************/
/*
   InitEW	- Init this file

   SYNOPSIS:
   INT InitEW ();

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function inits this file.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    __LINE__ if error occured.
 */
/****************************************************************************/

INT InitEWN ()
{
  INT i;

  if (CreateClass(EW_SOLVER_CLASS_NAME ".ewn",sizeof(NP_EWN),EWNConstruct)) return (__LINE__);

  return (0);
}
