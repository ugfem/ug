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

#include "devices.h"
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
  char dep[32];
  char mode[32];
  char cut[32];

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
  for (i=0; i<n; i++) GRID_LINK_VECTOR(theGrid,vlist[i],PrioMaster);
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

  if (ReadArgvINT("comp",&(np->comp),argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvChar ("dep",np->dep,argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvChar ("mode",np->mode,argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);
  if (ReadArgvChar ("cut",np->cut,argc,argv)) REP_ERR_RETURN(NP_NOT_ACTIVE);


  return (ORDER_Init(theNP,argc,argv));
}

INT OrderSODisplay (NP_BASE *theNP)
{
  NP_ORDER_SO *np;

  np = (NP_ORDER_SO *) theNP;

  UserWriteF(DISPLAY_NP_FORMAT_SI,"comp",(int)np->comp);
  UserWriteF(DISPLAY_NP_FORMAT_SS,"dep",np->dep);
  UserWriteF(DISPLAY_NP_FORMAT_SS,"mode",np->mode);
  UserWriteF(DISPLAY_NP_FORMAT_SS,"cut",np->cut);

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
      offdiag = MVALUE(theM,comp) - MVALUE(MADJ(theM),comp);
      if (diag*offdiag >= 0.0) SETMDOWN(theM,1);
      else SETMDOWN(theM,0);
    }
  }

  return (0);
}

/* static variables for mode DFCFCLL */
static VECTOR **DFCFCLL_vlist;
static INT DFCFCLL_MarkKey,DFCFCLL_f_nextout,DFCFCLL_f_nextin,DFCFCLL_l_nextout,DFCFCLL_l_nextin,DFCFCLL_f_inserted,DFCFCLL_l_inserted,DFCFCLL_n;

static INT DFCFCLL_PutFirst (GRID *theGrid, VECTOR *theV)
{
  MATRIX *theM;

  DFCFCLL_vlist[DFCFCLL_f_nextin++] = theV;
  SETVCUSED(theV,1);
  GRID_UNLINK_VECTOR(theGrid,theV);
  for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
    if (MDOWN(theM) && VCUSED(MDEST(theM))==0)
    {
      assert(VUP(MDEST(theM))>0);
      SETVUP(MDEST(theM),(VUP(MDEST(theM))-1));
    }
  DFCFCLL_f_inserted++;

  return (0);
}

static INT DFCFCLL_PutLast (GRID *theGrid, VECTOR *theV)
{
  MATRIX *theM;

  DFCFCLL_vlist[DFCFCLL_l_nextin--] = theV;
  SETVCUSED(theV,1);
  GRID_UNLINK_VECTOR(theGrid,theV);
  for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
    if (MUP(theM) && VCUSED(MDEST(theM))==0)
    {
      assert(VDOWN(MDEST(theM))>0);
      SETVDOWN(MDEST(theM),(VDOWN(MDEST(theM))-1));
    }
  DFCFCLL_l_inserted++;

  return (0);
}

static INT DFCFCLL_ModePre (GRID *theGrid)
{
  HEAP *theHeap;
  VECTOR *theV;
  INT n;

  /* count vectors and reset VCUSED */
  for (theV=FIRSTVECTOR(theGrid),n=0; theV!=NULL; theV=SUCCVC(theV),n++)
    SETVCUSED(theV,0);
  DFCFCLL_n = n;

  /* allocate list for vectors */
  theHeap = MGHEAP(MYMG(theGrid));
  MarkTmpMem(theHeap,&DFCFCLL_MarkKey);
  DFCFCLL_vlist = (VECTOR**)GetTmpMem(theHeap,sizeof(VECTOR*)*n,DFCFCLL_MarkKey);

  /* initialize vlist with first and last */
  DFCFCLL_f_nextout = DFCFCLL_f_nextin = 0;
  DFCFCLL_l_nextout = DFCFCLL_l_nextin = n-1;
  DFCFCLL_f_inserted = DFCFCLL_l_inserted = 0;
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
  {
    if (VUP(theV)==0)
    {
      if (DFCFCLL_PutFirst(theGrid,theV))
        return (1);
      continue;
    }
    if (VDOWN(theV)==0)
    {
      if (DFCFCLL_PutLast(theGrid,theV))
        return (1);
      continue;
    }
  }

  return (0);
}

static INT DFCFCLL_Mode (GRID *theGrid, VECTOR ***c_list, INT *f_inserted, INT *l_inserted, INT *v_remaining)
{
  INT i;
  MATRIX *theM;
  VECTOR *theV;

  /* fill in first */
  for (i=DFCFCLL_f_nextout; i<DFCFCLL_f_nextin; i++)
    for (theM=MNEXT(VSTART(DFCFCLL_vlist[i])); theM!=NULL; theM=MNEXT(theM))
    {
      theV = MDEST(theM);
      if (VCUSED(theV)) continue;
      if (VUP(theV)>0) continue;
      if (DFCFCLL_PutFirst(theGrid,theV))
        return (1);
    }
  *f_inserted = DFCFCLL_f_inserted; DFCFCLL_f_inserted = 0;

  /* fill in last */
  for (i=DFCFCLL_l_nextout; i>DFCFCLL_f_nextin; i--)
    for (theM=MNEXT(VSTART(DFCFCLL_vlist[i])); theM!=NULL; theM=MNEXT(theM))
    {
      theV = MDEST(theM);
      if (VCUSED(theV)) continue;
      if (VDOWN(theV)>0) continue;
      if (DFCFCLL_PutLast(theGrid,theV))
        return (1);
    }
  *l_inserted = DFCFCLL_l_inserted; DFCFCLL_l_inserted = 0;

  /* remaining vectors to order */
  *v_remaining = DFCFCLL_l_nextin - DFCFCLL_f_nextin + 1;

  return (0);
}

static INT DFCFCLL_ModePost (GRID *theGrid)
{
  HEAP *theHeap;
  INT i;

  for (i=0; i<DFCFCLL_n; i++)
    GRID_LINK_VECTOR(theGrid,DFCFCLL_vlist[i],PrioMaster);

  theHeap = MGHEAP(MYMG(theGrid));
  ReleaseTmpMem(theHeap,DFCFCLL_MarkKey);

  return (0);
}

static INT Cut_CenterToBnd (GRID *theGrid, MATDATA_DESC *A, INT comp, VECTOR **c_list, INT *c_inserted)
{
  VECTOR *theV,*theStart;
  MATRIX *theM;
  DOUBLE min,max,asp;
  INT i;
  DOUBLE_VECTOR pos;

  /* try to find center for start, i.e. a vector in the inner with minimum antisymmetric part */
  theStart = NULL; min = 1e30;
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
  {
    assert (VCUSED(theV));
    asp = 0.0;
    for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
    {
      asp += ABS(MVALUE(theM,comp) - MVALUE(MADJ(theM),comp));
      if (VCUSED(MDEST(theM)))
        break;
    }
    if (theM!=NULL) continue;

    if (asp < min)
    {
      min = asp;
      theStart = theV;
    }
  }

  /* if there is no center take the leftmost vector */
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
  {
    if (VectorPosition(theV,pos)) return (1);
    if (pos[0]<min)
    {
      min = pos[0];
      theStart = theV;
    }
  }
  assert(theStart);

  /* run in x-direction and fill in cut-vectors */
  for (i=0; theStart!=NULL; theStart=theV)
  {
    c_list[i++] = theStart;
    if (VectorPosition(theStart,pos)) return (1);
    theV = NULL; max = pos[0];
    for (theM=MNEXT(VSTART(theStart)); theM!=NULL; theM=MNEXT(theM))
    {
      if (VCUSED(MDEST(theM))) continue;
      if (VectorPosition(MDEST(theM),pos)) return (1);
      if (pos[0]>max)
      {
        max = pos[0];
        theV = MDEST(theM);
      }
    }
  }
  *c_inserted = i;

  return (0);
}

static INT OrderSO (NP_ORDER *theNP, INT level, MATDATA_DESC *A, INT *result)
{
  GRID *theGrid;
  NP_ORDER_SO *np;
  VECTOR *theV,**vlist,**c_list;
  MATRIX *theM;
  INT i,comp,n,MarkKey,f_inserted,l_inserted,v_remaining,c_inserted;
  MatDepProcPtr MatDepProc;
  CutProcPtr CutProc;
  ModePreProcPtr ModePreProc;
  ModeProcPtr ModeProc;
  ModePostProcPtr ModePostProc;

  np = (NP_ORDER_SO *) theNP;
  theGrid = NP_GRID(theNP,level);
  A = np->order.A;
  if (A==NULL) return (1);
  comp = MD_MCMP_OF_RT_CT(A,NODEVEC,NODEVEC,np->comp);

  /* configure */
  if (strcmp(np->dep,"adj")==0) MatDepProc  = MatrixDep_Adjoint;
  else return (1);
  if (strcmp(np->mode,"dfcfcll")==0)
  {
    ModePreProc     = DFCFCLL_ModePre;
    ModeProc        = DFCFCLL_Mode;
    ModePostProc    = DFCFCLL_ModePost;
  }
  else return (1);
  if (strcmp(np->cut,"ctb")==0) CutProc     = Cut_CenterToBnd;
  else return (1);

  /* set matrix dependencies */
  if ((*MatDepProc)(theGrid,A,comp)) return (1);

  /* count matrix dependencies */
  for (theV=FIRSTVECTOR(theGrid),n=0; theV!=NULL; theV=SUCCVC(theV),n++)
  {
    SETVDOWN(theV,0);
    SETVUP(theV,0);
    for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM))
      if (MDOWN(theM))
      {
        assert(0);
        SETVDOWN(theV,(VDOWN(theV)+1));                       /* number of vectors depending on theV */
        SETVUP(MDEST(theM),(VUP(MDEST(theM))+1));             /* number of vectors depending on theV */
      }
      else if (MDOWN(MADJ(theM)))
        SETVUP(theV,(VUP(theV)+1));                 /* number of vectors theV depends on   */
  }

  /* order */
  if ((*ModePreProc)(theGrid)) return (1);
  for (i=0;; i++)
  {
    if ((*ModeProc)(theGrid,&c_list,&f_inserted,&l_inserted,&v_remaining)) return (1);
    UserWriteF("%d: FIRST: %d\n",i,f_inserted);
    UserWriteF("   LAST:  %d\n",f_inserted);
    UserWriteF("   REM:   %d\n",v_remaining);

    if (v_remaining==0) break;
    if ((*CutProc)(theGrid,A,comp,c_list,&c_inserted)) return (1);
    UserWriteF("   CUT:   %d\n",c_inserted);
  }
  if ((*ModePostProc)(theGrid)) return (1);

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
  /*	if (CreateClass(ORDER_CLASS_NAME ".so",sizeof(NP_ORDER_SO),OrderSO_Construct))   return (__LINE__); */

  return (0);
}
