// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  order.c                                                                                                       */
/*																			*/
/* Purpose:   order vectors                                                                             */
/*																			*/
/* Author:	  Klaus Johannsen                                                                       */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   Nov 18, 1997 begin                                                                */
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

#include "ugdevices.h"
#include "ugenv.h"

#include "scan.h"
#include "numproc.h"
#include "np.h"
#include "ugm.h"
#include "evm.h"
#include "general.h"
#include "fileopen.h"
#include "ugstruct.h"
#include "fifo.h"

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
  NP_ORDER order;

  char config[32];
  INT mode;
  INT ordering[3];
  INT sign[3];
  INT which;
  INT AlsoOrderMatrices;
  INT SpecialTreatSkipVecs;

} NP_ORDER_LEX;

typedef struct
{
  NP_ORDER order;

  INT bandwidth;

} NP_ORDER_BW;

typedef INT (*MatDepProcPtr)(GRID *theGrid, MATDATA_DESC *A, INT comp);
typedef INT (*CutProcPtr)(GRID *theGrid, MATDATA_DESC *A, INT comp, VECTOR **c_list, INT *c_inserted);
typedef INT (*ModePreProcPtr)(GRID *theGrid);
typedef INT (*ModeProcPtr)(GRID *theGrid, VECTOR ***c_list, INT *f_inserted, INT *l_inserted, INT *v_remaining);
typedef INT (*ModePostProcPtr)(GRID *theGrid);

typedef struct
{
  NP_ORDER order;

  INT comp;
  INT ncyc;
  INT ncut;

} NP_ORDER_SO;

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

static void PVP (VECTOR *v)
{
  DOUBLE_VECTOR pos;

  VectorPosition(v,pos);
  printf("%f %f\n",pos[0],pos[1]);

  return;
}

static void PVPI (VECTOR *v, INT f)
{
  DOUBLE_VECTOR pos;
  INT ipos[2];

  VectorPosition(v,pos);
  ipos[0]=(INT)(pos[0]*f+0.5);
  ipos[1]=(INT)(pos[1]*f+0.5);
  printf("%d %d\n",ipos[0],ipos[1]);

  return;
}

static void PVPC (VECTOR *v, char *c)
{
  DOUBLE_VECTOR pos;

  VectorPosition(v,pos);
  printf("%s: %f %f\n",c,pos[0],pos[1]);

  return;
}

static void PVPCI (VECTOR *v, char *c, INT f)
{
  DOUBLE_VECTOR pos;
  INT ipos[2];

  VectorPosition(v,pos);
  ipos[0]=(INT)(pos[0]*f+0.5);
  ipos[1]=(INT)(pos[1]*f+0.5);
  printf("%s: %d %d\n",ipos[0],ipos[1]);

  return;
}

INT ORDER_Init (NP_BASE *theNP, INT argc, char **argv)
{
  NP_ORDER *np;

  np = (NP_ORDER *)theNP;
  np->A = ReadArgvMatDesc(np->base.mg,"A",argc,argv);
  if (np->A == NULL) return(NP_ACTIVE);

  return (NP_EXECUTABLE);
}

INT ORDER_Display (NP_BASE *theNP)
{
  NP_ORDER *np;

  np = (NP_ORDER *)theNP;
  UserWrite("symbolic user data:\n");
  if (np->A != NULL) UserWriteF(DISPLAY_NP_FORMAT_SS,"A",ENVITEM_NAME(np->A));

  return (0);
}

INT NPOrderExecute (NP_BASE *theNP, INT argc , char **argv)
{
  NP_ORDER *np;
  INT result,i,fromlevel,tolevel;

  np = (NP_ORDER *) theNP;
  tolevel = CURRENTLEVEL(theNP->mg);

  if (ReadArgvOption("a",argc,argv)) fromlevel = 0;
  else fromlevel = tolevel;

  for (i=fromlevel; i<=tolevel; i++)
    if ((*(np->Order))(np,i,np->A,&result))
      REP_ERR_RETURN (1);

  return(0);
}

/****************************************************************************/
/*D
   order lex - numproc order lex

   DESCRIPTION:
   Reads a file with double values.

   'npinit <name> $config <xx>'

   .  <name> - num proc name
   .  $f~<file> - file name
   D*/
/****************************************************************************/

static INT OrderLexInit (NP_BASE *theNP, INT argc , char **argv)
{
  INT i;
  NP_ORDER_LEX *np;

  np = (NP_ORDER_LEX *) theNP;
  if (ReadArgvChar ("config",np->config,argc,argv)) return (NP_NOT_ACTIVE);
  if (strlen(np->config)!=DIM) return (NP_NOT_ACTIVE);
  np->which = GM_TAKE_SKIP | GM_TAKE_NONSKIP;
  np->AlsoOrderMatrices = FALSE;
  np->SpecialTreatSkipVecs = FALSE;
  np->mode = OV_CARTES;

  for (i=0; i<strlen(np->config); i++)
    switch(np->config[i])
    {
    case 'r' :
      np->ordering[i] = _X_;
      np->sign[i] = 1;
      break;
    case 'l' :
      np->ordering[i] = _X_;
      np->sign[i] = -1;
      break;
    case 'u' :
      np->ordering[i] = _Y_;
      np->sign[i] = 1;
      break;
    case 'd' :
      np->ordering[i] = _Y_;
      np->sign[i] = -1;
      break;
#ifdef __THREEDIM__
    case 'f' :
      np->ordering[i] = _Z_;
      np->sign[i] = 1;
      break;
    case 'b' :
      np->ordering[i] = _Z_;
      np->sign[i] = -1;
      break;
#endif
    default :
      return (NP_NOT_ACTIVE);
    }

  return (NP_EXECUTABLE);
}

INT OrderLexDisplay (NP_BASE *theNP)
{
  NP_ORDER_LEX *np;

  np = (NP_ORDER_LEX *) theNP;
  UserWriteF(DISPLAY_NP_FORMAT_SS,"config",np->config);

  return (0);
}

static INT OrderLex (NP_ORDER *theNP, INT level, MATDATA_DESC *A, INT *result)
{
  GRID *theGrid;
  NP_ORDER_LEX *np;

  np = (NP_ORDER_LEX *) theNP;
  theGrid = NP_GRID(theNP,level);
  if (LexOrderVectorsInGrid(theGrid,np->mode,np->ordering,np->sign,np->which,np->SpecialTreatSkipVecs,np->AlsoOrderMatrices)!=GM_OK)
    return (1);

  return(0);
}

static INT OrderLex_Construct (NP_BASE *theNP)
{
  NP_ORDER_LEX *np;

  theNP->Init =    OrderLexInit;
  theNP->Display = OrderLexDisplay;
  theNP->Execute = NPOrderExecute;

  np = (NP_ORDER_LEX *) theNP;
  np->order.Order = OrderLex;

  return(0);
}

/****************************************************************************/
/*D
   order bandwidth - numproc order bw

   DESCRIPTION:
   Reads a file with double values.

   'npinit <name> $config <xx>'

   .  <name> - num proc name
   .  $f~<file> - file name
   D*/
/****************************************************************************/

static INT OrderBWInit (NP_BASE *theNP, INT argc , char **argv)
{
  return (NP_EXECUTABLE);
}

INT OrderBWDisplay (NP_BASE *theNP)
{
  NP_ORDER_BW *np;

  np = (NP_ORDER_BW *) theNP;
  UserWriteF(DISPLAY_NP_FORMAT_SI,"bandwidth",np->bandwidth);

  return (0);
}

static INT OrderBW (NP_ORDER *theNP, INT level, MATDATA_DESC *A, INT *result)
{
  GRID *theGrid;
  HEAP *theHeap;
  VECTOR *theV,**vlist;
  MATRIX *theM;
  INT i,n,MarkKey,bw,k,index;
  void *buffer;
  FIFO myfifo;
  NP_ORDER_BW *np;

  np = (NP_ORDER_BW *) theNP;
  theGrid = NP_GRID(theNP,level);
  theHeap = MGHEAP(MYMG(theGrid));
  n = 0; for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV)) n++;

  /* reorder vector-list */
  MarkTmpMem(theHeap,&MarkKey);
  buffer=(void *)GetTmpMem(theHeap,sizeof(VECTOR*)*n,MarkKey);
  vlist = (VECTOR**)GetTmpMem(theHeap,sizeof(VECTOR*)*n,MarkKey);
  fifo_init(&myfifo,buffer,sizeof(VECTOR*)*n);
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    SETVCUSED(theV,0);
  fifo_in(&myfifo,(void *)FIRSTVECTOR(theGrid));
  SETVCUSED(FIRSTVECTOR(theGrid),1);
  while(!fifo_empty(&myfifo))
  {
    theV = (VECTOR *)fifo_out(&myfifo);
    for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
      if (!VCUSED(MDEST(theM)))
      {
        fifo_in(&myfifo,(void *)MDEST(theM));
        SETVCUSED(MDEST(theM),1);
      }
  }
  fifo_in(&myfifo,(void *)theV);
  SETVCUSED(theV,0); i=0;
  while(!fifo_empty(&myfifo))
  {
    theV = (VECTOR *)fifo_out(&myfifo);
    vlist[i++] = theV;
    for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
      if (VCUSED(MDEST(theM)))
      {
        fifo_in(&myfifo,(void *)MDEST(theM));
        SETVCUSED(MDEST(theM),0);
      }
  }
  assert(i==n);
  for (i=0; i<n; i++) GRID_UNLINK_VECTOR(theGrid,vlist[i]);
  for (i=0; i<n; i++) GRID_LINK_VECTOR(theGrid,vlist[i],PRIO(vlist[i]));
  ReleaseTmpMem(theHeap,MarkKey);

  /* determine bandwidth */
  for (theV=FIRSTVECTOR(theGrid),i=0; theV!=NULL; theV=SUCCVC(theV),i++)
    VINDEX(theV) = i;
  bw = 0;
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
  {
    index = VINDEX(theV);
    for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
    {
      k = index-MDESTINDEX(theM);
      k = ABS(k);
      bw = MAX(bw,k);
    }
  }
  np->bandwidth = bw;

  return(0);
}

static INT OrderBW_Construct (NP_BASE *theNP)
{
  NP_ORDER_LEX *np;

  theNP->Init =    OrderBWInit;
  theNP->Display = OrderBWDisplay;
  theNP->Execute = NPOrderExecute;

  np = (NP_ORDER_LEX *) theNP;
  np->order.Order = OrderBW;

  return(0);
}

/****************************************************************************/
/*D
   order lex - numproc order lex

   DESCRIPTION:
   Reads a file with double values.

   'npinit <name> $config <xx>'

   .  <name> - num proc name
   .  $f~<file> - file name
   D*/
/****************************************************************************/

static INT OrderSOInit (NP_BASE *theNP, INT argc , char **argv)
{
  NP_ORDER_SO *np;

  np = (NP_ORDER_SO *) theNP;
  if (ReadArgvINT ("comp",&(np->comp),argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);

  return (ORDER_Init(theNP,argc,argv));
}

INT OrderSODisplay (NP_BASE *theNP)
{
  NP_ORDER_SO *np;

  np = (NP_ORDER_SO *) theNP;
  UserWriteF(DISPLAY_NP_FORMAT_SI,"comp",(int)np->comp);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"ncyc",(int)np->ncyc);
  UserWriteF(DISPLAY_NP_FORMAT_SI,"ncut",(int)np->ncut);

  return (ORDER_Display(theNP));
}

static INT MatrixDep_Adjoint (GRID *theGrid, MATDATA_DESC *A, INT comp)
{
  VECTOR *theV;
  MATRIX *theM;
  DOUBLE diag, offdiag;

  /* set matrix flags */
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
  {
    diag = MVALUE(VSTART(theV),comp);
    for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
    {
      SETMUP(theM,0);
      offdiag = MVALUE(theM,comp) - MVALUE(MADJ(theM),comp);
      if (diag*offdiag > 0.0) SETMDOWN(theM,1);
      else SETMDOWN(theM,0);
    }
  }

  return (0);
}

#define MUPP(m)         MDOWN(MADJ(m))
#define WH_INIT              0
#define WH_CYCLE_INIT        1
#define WH_ACYCLE_INIT       2
#define WH_DOWN              0
#define WH_UP                1

typedef struct {
  VECTOR *v;
  DOUBLE angle;
  INT mode;

} MORDER;

int WH_mcmp (const void* x, const void *y)
{
  MORDER *mx,*my;

  mx=(MORDER*)x;
  my=(MORDER*)y;
  if (mx->angle<my->angle) return (-1);
  if (mx->angle>my->angle) return (1);
  return (0);
}

static INT WH_IsStarVector (VECTOR *v)
{
  INT i,n,ndown,ndu,nud;
  MATRIX *m;
  MORDER mo[30];
  DOUBLE_VECTOR z,p;

  n=ndown=0;
  VectorPosition(v,z);
  for (m=MNEXT(VSTART(v)); m!=NULL; m=MNEXT(m))
  {
    if (MDOWN(m) && !VCUSED(MDEST(m)))
    {
      mo[n].mode=WH_DOWN;
      ndown++;
    }
    else if (MUPP(m))
      mo[n].mode=WH_UP;
    else
      continue;
    VectorPosition(MDEST(m),p);
    p[0]-=z[0];
    p[1]-=z[1];
    mo[n].angle=atan2(p[1],p[0]);
    n++;
  }
  if (ndown<=1) return (0);
  if (n<=3) return (0);
  qsort((void*)mo,n,sizeof(MORDER),WH_mcmp);
  ndu=nud=0;
  for (i=0; i<n; i++)
    if (mo[i].mode==WH_DOWN && mo[(i+1)%n].mode==WH_UP)
      ndu++;
    else if (mo[i].mode==WH_UP && mo[(i+1)%n].mode==WH_DOWN)
      nud++;
  if (nud!=1 || ndu!=1) return (1);
  return (0);
}

static INT WH_PrintStarVectors (GRID *g)
{
  VECTOR *v;

  for (v=FIRSTVECTOR(g); v!=NULL; v=SUCCVC(v))
    if (!VCUSED(v))
      if (WH_IsStarVector(v))
        PVP(v);
  return(1);
}

static VECTOR *WH_next_left (VECTOR *v)
{
  INT i,n,ndown,ndu,nud;
  VECTOR *w,*s;
  MATRIX *m;
  MORDER mo[30];
  DOUBLE_VECTOR z,p;

  n=ndown=0;
  VectorPosition(v,z);
  for (m=MNEXT(VSTART(v)); m!=NULL; m=MNEXT(m))
  {
    if (MDOWN(m) && !VCUSED(MDEST(m)))
    {
      mo[n].mode=WH_DOWN;
      ndown++;
      mo[n].v=s=MDEST(m);
    }
    else if (MUPP(m))
      mo[n].mode=WH_UP;
    else
      continue;
    VectorPosition(MDEST(m),p);
    p[0]-=z[0];
    p[1]-=z[1];
    mo[n].angle=atan2(p[1],p[0]);
    n++;
  }
  if (ndown==0) return (NULL);
  if (ndown==1) return (s);
  qsort((void*)mo,n,sizeof(MORDER),WH_mcmp);
  ndu=nud=0;
  for (i=0; i<n; i++)
    if (mo[i].mode==WH_DOWN && mo[(i+1)%n].mode==WH_UP)
    {
      ndu++;
      w=mo[i].v;
    }
    else if (mo[i].mode==WH_UP && mo[(i+1)%n].mode==WH_DOWN)
      nud++;
  if (nud!=1 || ndu!=1)
  {
    assert(0);
  }
  return (w);
}

static VECTOR *WH_next_right (VECTOR *v)
{
  INT i,n,ndown,ndu,nud;
  VECTOR *w,*s;
  MATRIX *m;
  MORDER mo[30];
  DOUBLE_VECTOR z,p;

  n=ndown=0;
  VectorPosition(v,z);
  for (m=MNEXT(VSTART(v)); m!=NULL; m=MNEXT(m))
  {
    if (MDOWN(m) && !VCUSED(MDEST(m)))
    {
      mo[n].mode=WH_DOWN;
      ndown++;
      mo[n].v=s=MDEST(m);
    }
    else if (MUPP(m))
      mo[n].mode=WH_UP;
    else
      continue;
    VectorPosition(MDEST(m),p);
    p[0]-=z[0];
    p[1]-=z[1];
    mo[n].angle=atan2(p[1],p[0]);
    n++;
  }
  if (ndown==0) return (NULL);
  if (ndown==1) return (s);
  qsort((void*)mo,n,sizeof(MORDER),WH_mcmp);
  ndu=nud=0;
  for (i=0; i<n; i++)
    if (mo[i].mode==WH_DOWN && mo[(i+1)%n].mode==WH_UP)
      ndu++;
    else if (mo[i].mode==WH_UP && mo[(i+1)%n].mode==WH_DOWN)
    {
      nud++;
      w=mo[(i+1)%n].v;
    }
  if (nud!=1 || ndu!=1)
  {
    assert(0);
  }

  return (w);
}

static INT SimpleCut (GRID *g, VECTOR **vlist, INT *ncut)
{
  INT n,nl,nr;
  VECTOR *v0,*v1,*vl,*vr;
  MATRIX *m;
  DOUBLE_VECTOR p0,p1;
  DOUBLE min;

  /* run left */
  for (vl=FIRSTVECTOR(g); vl!=NULL; vl=SUCCVC(vl))
  {
    SETVCFLAG(vl,0);
    SETVCCUT(vl,0);
  }
  for (vl=FIRSTVECTOR(g),nl=0; vl!=NULL; vl=WH_next_left(vl))
    if (!VCFLAG(vl))
      SETVCFLAG(vl,1);
    else if (!VCCUT(vl))
    {
      SETVCCUT(vl,1);
      nl++;
    }
    else
      break;
  assert(vl!=NULL);

  /* run right */
  for (vr=FIRSTVECTOR(g); vr!=NULL; vr=SUCCVC(vr))
  {
    SETVCFLAG(vr,0);
    SETVCCUT(vr,0);
  }
  for (vr=FIRSTVECTOR(g); vr!=NULL; vr=SUCCVC(vr))
  {
    SETVCFLAG(vr,0);
    SETVCCUT(vr,0);
  }
  for (vr=FIRSTVECTOR(g),nr=0; vr!=NULL; vr=WH_next_right(vr))
    if (!VCFLAG(vr))
      SETVCFLAG(vr,1);
    else if (!VCCUT(vr))
    {
      SETVCCUT(vr,1);
      nr++;
    }
    else
      break;
  assert(vr!=NULL);

  if (nl<nr) v0=vl;
  else v0=vr;
  n=0;
  while(1)
  {
    vlist[n++]=v0;
    VectorPosition(v0,p0);
    v1=NULL; min=p0[0];
    for (m=MNEXT(VSTART(v0)); m!=NULL; m=MNEXT(m))
      if (!VCUSED(MDEST(m)))
      {
        VectorPosition(MDEST(m),p1);
        if (p1[0]<min)
        {
          min=p1[0];
          v1=MDEST(m);
        }
      }
    if (v1==NULL) break;
    v0=v1;
  }
  *ncut=n;

  return(0);
}

static void FirstInsertInVList (GRID *g, VECTOR *v, VECTOR **vlist, INT n, INT unlink)
{
  MATRIX *theM;

  vlist[n]=v;
  SETVCUSED(v,1);
  for (theM=MNEXT(VSTART(v)); theM!=NULL; theM=MNEXT(theM))
    if (MDOWN(theM) && !VCUSED(MDEST(theM)))
      SETVUP(MDEST(theM),VUP(MDEST(theM))-1);
  if (unlink)
    GRID_UNLINK_VECTOR(g,v);
  SETVCFLAG(v,0);

  return;
}

static void LastInsertInVList (GRID *g, VECTOR *v, VECTOR **vlist, INT n, INT unlink)
{
  MATRIX *theM;

  vlist[n]=v;
  SETVCUSED(v,1);
  for (theM=MNEXT(VSTART(v)); theM!=NULL; theM=MNEXT(theM))
    if (MUPP(theM) && !VCUSED(MDEST(theM)))
      SETVDOWN(MDEST(theM),VDOWN(MDEST(theM))-1);
  if (unlink)
    GRID_UNLINK_VECTOR(g,v);
  SETVCFLAG(v,0);

  return;
}

static INT OrderSO (NP_ORDER *theNP, INT level, MATDATA_DESC *A, INT *result)
{
  HEAP *theHeap;
  GRID *theGrid;
  NP_ORDER_SO *np;
  VECTOR *theV,**vlist;
  MATRIX *theM,*theM0;
  INT i,cnt,comp,ncut,n,MarkKey,fno,fni,lno,lni;
  DOUBLE pos[2];

  np = (NP_ORDER_SO *) theNP;
  np->ncut = np->ncyc = 0;
  theGrid = NP_GRID(theNP,level);
  A = np->order.A;
  if (A==NULL) return (1);
  comp = MD_MCMP_OF_RT_CT(A,NODEVEC,NODEVEC,np->comp);

  /* set matrix dependencies */
  if (MatrixDep_Adjoint(theGrid,A,comp)) return (1);

  /* count matrix dependencies */
  for (theV=FIRSTVECTOR(theGrid),n=0; theV!=NULL; theV=SUCCVC(theV),n++)
  {
    SETVCUSED(theV,0);
    SETVDOWN(theV,0);
    SETVUP(theV,0);
    for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
      if (MDOWN(theM))
      {
        SETVDOWN(theV,(VDOWN(theV)+1));                         /* number of vectors depending on theV */
      }
      else if (MDOWN(MADJ(theM)))
      {
        SETVUP(theV,(VUP(theV)+1));                                     /* number of vectors theV depends on   */
      }
  }

  /* allocate list for vectors */
  theHeap = MGHEAP(MYMG(theGrid));
  MarkTmpMem(theHeap,&MarkKey);
  vlist = (VECTOR**)GetTmpMem(theHeap,sizeof(VECTOR*)*n,MarkKey);
  assert(vlist!=NULL);

  /* insert first/last set */
  fni = fno = 0; lni = lno = n-1;
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
  {
    if (VUP(theV)==NULL) FirstInsertInVList(theGrid,theV,vlist,fni++,0);
    else if (VDOWN(theV)==NULL) LastInsertInVList(theGrid,theV,vlist,lni--,0);
  }
  for (i=fno; i<fni; i++) GRID_UNLINK_VECTOR(theGrid,vlist[i]);
  for (i=lno; i>lni; i--) GRID_UNLINK_VECTOR(theGrid,vlist[i]);

  /* order */
  for (cnt=0;; cnt++)
  {
    /* process first set */
    for (i=fno; i<fni; i++)
      for (theM0=MNEXT(VSTART(vlist[i])); theM0!=NULL; theM0=MNEXT(theM0))
      {
        theV=MDEST(theM0);
        if (!VCUSED(theV) && MDOWN(theM0) && VUP(theV)==0)
          FirstInsertInVList(theGrid,theV,vlist,fni++,1);
      }
    fno=fni;

    /* process last set */
    for (i=lno; i>lni; i--)
      for (theM0=MNEXT(VSTART(vlist[i])); theM0!=NULL; theM0=MNEXT(theM0))
      {
        theV=MDEST(theM0);
        if (!VCUSED(theV) && MUPP(theM0) && VDOWN(theV)==0)
          LastInsertInVList(theGrid,theV,vlist,lni--,1);
      }
    lno=lni;

    /* break if vlist complete */
    if (fni>lni) break;

    /* cut star vectors */
    if (!cnt)
      for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
        if (!VCUSED(theV) && WH_IsStarVector(theV))
        {
          FirstInsertInVList(theGrid,theV,vlist,fni,0);
          LastInsertInVList(theGrid,theV,vlist,fni++,1);
          np->ncut++;
          SETVCFLAG(theV,1);
        }

    /* insert cut as first set */
    if (SimpleCut(theGrid,vlist+fni,&ncut)) return(1);
    assert(ncut>0); assert(ncut<=lni-fni+1);
    np->ncut+=ncut; np->ncyc++;
    for (i=fni; i<fni+ncut; i++)
    {
      FirstInsertInVList(theGrid,vlist[i],vlist,i,0);
      LastInsertInVList(theGrid,vlist[i],vlist,i,1);
      SETVCFLAG(vlist[i],1);
    }
    fni+=ncut;

    /* insert cut-neighbors as last set */
    for (i=fni-ncut; i<fni; i++)
      for (theM0=MNEXT(VSTART(vlist[i])); theM0!=NULL; theM0=MNEXT(theM0))
      {
        theV=MDEST(theM0);
        if (!VCUSED(theV) && MUPP(theM0) && VDOWN(theV)==0)
          LastInsertInVList(theGrid,theV,vlist,lni--,1);
      }
  }

  /* postprocess */
  for (i=0; i<n; i++)
    GRID_LINK_VECTOR(theGrid,vlist[i],PRIO(vlist[i]));
  ReleaseTmpMem(theHeap,MarkKey);

  return(0);
}

static INT OrderSO_Construct (NP_BASE *theNP)
{
  NP_ORDER_LEX *np;

  theNP->Init =    OrderSOInit;
  theNP->Display = OrderSODisplay;
  theNP->Execute = NPOrderExecute;

  np = (NP_ORDER_LEX *) theNP;
  np->order.Order = OrderSO;

  return(0);
}

/****************************************************************************/
/*D
   InitOrder - Enrol order numprocs

   SYNOPSIS:
   INT InitOrder (void);

   PARAMETERS:
   .  void -

   DESCRIPTION:
   This function creates the numproc classes for ordering.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured.
   D*/
/****************************************************************************/

INT InitOrder (void)
{
  if (CreateClass(ORDER_CLASS_NAME ".lex",sizeof(NP_ORDER_LEX),OrderLex_Construct)) return (__LINE__);
  if (CreateClass(ORDER_CLASS_NAME ".bw",sizeof(NP_ORDER_BW),OrderBW_Construct)) return (__LINE__);
  if (CreateClass(ORDER_CLASS_NAME ".so",sizeof(NP_ORDER_SO),OrderSO_Construct)) return (__LINE__);

  return (0);
}
