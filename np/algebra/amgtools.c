// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:	  amgtools.c													*/
/*                                                                          */
/* Purpose:   tools for algebraic multigrid                                         */
/*                                                                          */
/* Author:	  Nicolas Neuss                                                                                     */
/*			  Institut fuer Angewandte Mathematik                           */
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 294										*/
/*			  69120 Heidelberg												*/
/*			  email: neuss@iwr.uni-heidelberg.de			                        */
/*																			*/
/* History:   1994-1995 in old ug2.0							            */
/*            May 1997  in new ug3.7                                        */
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "algebra.h"
#include "block.h"
#include "cmdline.h"
#include "compiler.h"
#include "debug.h"
#include "devices.h"
#include "evm.h"
#include "general.h"
#include "gm.h"
#include "heaps.h"
#include "misc.h"
#include "np.h"
#include "ugm.h"
#include "algebra.h"
#include "fifo.h"
#ifdef ModelP
#include "pargm.h"
#include "parallel.h"
#endif

#include "amgtools.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define REUSKEN 0
#define WAGNER 1

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);
REP_ERR_FILE;

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/* Some routines marking strong connections  (!change for systems!)         */
/****************************************************************************/

INT UnmarkAll(GRID *theGrid, MATDATA_DESC *A, DOUBLE theta, INT vcomp)
{
  VECTOR *vect;
  MATRIX *mat;

  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    for (mat=VSTART(vect); mat!=NULL; mat=MNEXT(mat))
      SETSTRONG(mat,0);
  return(0);
}

INT MarkAll(GRID *theGrid, MATDATA_DESC *A, DOUBLE theta, INT vcomp)
{
  VECTOR *vect;
  MATRIX *mat;

  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    for (mat=VSTART(vect); mat!=NULL; mat=MNEXT(mat))
      SETSTRONG(mat,1);
  return(0);
}

INT MarkOffDiagWithoutDirichlet(GRID *theGrid, MATDATA_DESC *A, DOUBLE theta, INT vcomp)
{
  VECTOR *vect;
  MATRIX *mat;

  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    if (VECSKIP(vect)==0)
      for (mat=VSTART(vect); mat!=NULL; mat=MNEXT(mat))
        if (VECSKIP(MDEST(mat))==0)
          SETSTRONG(mat,1);
  return(0);
}

INT MarkAbsolute(GRID *theGrid, MATDATA_DESC *A, DOUBLE theta, INT vcomp)
{
  VECTOR *vect;
  MATRIX *mat;
  INT mcomp;

  /* Check that only mtype=0x0 (for other, AMG is unclear anyhow).        */
  /* Define and set the variables scalar,blockN,blockNN, and error. */
  BLOCK_SETUP(A);
  if (error) {
    PrintErrorMessage('E',"MarkAbsolute","not yet for general matrices");
    REP_ERR_RETURN(error);
  }
  mcomp=MD_MCMP_OF_MTYPE(A,0,0);
  if (vcomp>MD_ROWS_IN_MTYPE(A,0)) {
    PrintErrorMessage('E',"MarkAbsolute","vcomp too large");
    REP_ERR_RETURN(error);
  }
  if (vcomp>=0)
    mcomp += vcomp*(MD_COLS_IN_MTYPE(A,0)+1);
  else
  {
    PrintErrorMessage('E',"MarkAbsolute","whole block handling not implemented for this marking");
    REP_ERR_RETURN(error);
  }

  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
  {
    /* if it is a Dirichlet vector there is no strong influence on the vector */
    if (VECSKIP(vect)) continue;

    /* the diagonal is special anyhow, we don't mark it... */
    for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
    {
      /* we also say there is no strong influence from Dirichlet vectors
         so that they wont't be used for coarse grid vectors */
      if (VECSKIP(MDEST(mat))==0)
        if (-MVALUE(mat,mcomp)>=theta)
          SETSTRONG(mat,1);
    }
  }

  return(0);
}

INT MarkRelative(GRID *theGrid, MATDATA_DESC *A, DOUBLE theta, INT vcomp)
{
  VECTOR *vect,*vect2;
  MATRIX *matD,*mat;
  INT mcomp;
  DOUBLE s,threshold;

  /* Check that only mtype=0x0 (for other, AMG is unclear anyhow).        */
  /* Define and set the variables scalar,blockN,blockNN, and error. */
  BLOCK_SETUP(A);
  if (error) {
    PrintErrorMessage('E',"MarkAbsolute","not yet for general matrices");
    REP_ERR_RETURN(error);
  }

  mcomp=MD_MCMP_OF_MTYPE(A,0,0);
  if (vcomp>MD_ROWS_IN_MTYPE(A,0)) {
    PrintErrorMessage('E',"MarkRelative","vcomp too large");
    REP_ERR_RETURN(error);
  }
  if (vcomp>=0)
    mcomp += vcomp*(MD_COLS_IN_MTYPE(A,0)+1);
  else
  {
    PrintErrorMessage('E',"MarkRelative","whole block handling not implemented for this marking");
    REP_ERR_RETURN(error);
  }

  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
  {
    /* if it is a Dirichlet vector there is no strong influence on the vector */
    if (VECSKIP(vect)) continue;

    matD=VSTART(vect);
    /* the diagonal is special anyhow, we don't mark it... */

    /* find maximum of -a_{ij} */
    s=0.0;
    for (mat=MNEXT(matD); mat!=NULL; mat=MNEXT(mat))
    {
      vect2=MDEST(mat);
      if (VECSKIP(vect2)==0)
        if (s < -MVALUE(mat,mcomp))
          s = -MVALUE(mat,mcomp);
    }
    threshold=s*theta;
    for (mat=MNEXT(matD); mat!=NULL; mat=MNEXT(mat))
    {
      /* we also say there is no strong influence from Dirichlet vectors
         so that they wont't be used for coarse grid vectors */
      if (VECSKIP(MDEST(mat))==0)
        if (-MVALUE(mat,mcomp)>=threshold)
          SETSTRONG(mat,1);
    }
  }

  return(0);
}

INT MarkVanek(GRID *theGrid, MATDATA_DESC *A, DOUBLE theta, INT vcomp)
{
  VECTOR *vect,*vect2;
  MATRIX *mat,*matD;
  DOUBLE norm_ii,norm_jj,norm_ij;
  INT mcomp;

  /* Check that only mtype=0x0 (for other, AMG is unclear anyhow).        */
  /* Define and set the variables scalar,blockN,blockNN, and error. */
  BLOCK_SETUP(A);
  if (error) {
    PrintErrorMessage('E',"MarkAbsolute","not yet for general matrices");
    REP_ERR_RETURN(error);
  }
  if (vcomp>MD_ROWS_IN_MTYPE(A,0)) {
    PrintErrorMessage('E',"MarkAbsolute","vcomp too large");
    REP_ERR_RETURN(error);
  }

  mcomp=MD_MCMP_OF_MTYPE(A,0,0);
  if (vcomp>0)
    mcomp+=vcomp*(MD_COLS_IN_MTYPE(A,0)+1);

  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
  {
    /* if it is a Dirichlet vector there is no strong influence on the vector */
    if (VECSKIP(vect)) continue;

    matD=VSTART(vect);
    if (vcomp<0)
    {BLOCK_NORM(&(MVALUE(matD,mcomp)),norm_ii);}
    else
      norm_ii=fabs(MVALUE(matD,mcomp));

    for (mat=MNEXT(matD); mat!=NULL; mat=MNEXT(mat))
    {
      vect2=MDEST(mat);

      /* we also say there is no strong influence from Dirichlet vectors
         so that they wont't be used for coarse grid vectors */
      if (VECSKIP(vect2)!=0) continue;

      matD=VSTART(vect2);
      if (vcomp<0) {
        BLOCK_NORM(&(MVALUE(matD,mcomp)),norm_jj);
        BLOCK_NORM(&(MVALUE(mat ,mcomp)),norm_ij);
      }
      else {
        norm_jj=fabs(MVALUE(matD,mcomp));
        norm_ij=fabs(MVALUE(mat,mcomp));
      }

      if (norm_ij>=theta*sqrt(norm_ii*norm_jj))
        SETSTRONG(mat,1);
    }
  }

  return(0);
}

/****************************************************************************/
/* Some routines marking strong connections  (!change for systems!)         */
/****************************************************************************/

INT SetupInitialList(GRID *theGrid, HEAP *theHeap, AVECTOR **initialSH, AVECTOR **initialEH, INT MarkKey)
{
  VECTOR *vect;
  AVECTOR *avect;

  /* sets index field and copies vectors to initial list */
  *initialSH=*initialEH=NULL;
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
  {
    SETVCUSED(vect,0);
    SETVCCOARSE(vect,0);

    if ((avect=(AVECTOR *) GetTmpMem(theHeap,sizeof(AVECTOR),MarkKey))==NULL)
    {
      PrintErrorMessage('E',"SetupInitialList",
                        "could not allocate avector");
      REP_ERR_RETURN(1);
    }

    CTRL(avect)=0;
    SETAVCOARSE(avect,0);
    SETAVFINE(avect,0);
    SETAVSKIP(avect,0);
    SETAVTESTED(avect,0);
    STRONG_IN(avect)=0;
    STRONG_OUT(avect)=0;

    VECT(avect)=vect;
    VISTART(vect)=(MATRIX *) avect;             /* VISTART field is used for establishing bijection */

    ADDATEND_LIST2(*initialSH,*initialEH,avect);
  }

  return(DONE);
}

INT DistributeInitialList(AVECTOR **La, AVECTOR **Le, AVECTOR **Da, AVECTOR **De, AVECTOR **Ua, AVECTOR **Ue)
{
  INT i;
  AVECTOR *avect;

  /* we sort the avects according to the number of points they influence
     we also put Dirichlet values straight away in the final tested fine grid */
  while ((avect=*La)!=NULL)
  {
    ELIMINATE_LIST2(*La,*Le,avect);
    if (VECSKIP(VECT(avect)) != 0)
    {
      SETAVFINE(avect,1);
      SETAVTESTED(avect,1);
      SETAVSKIP(avect,1);
      ADDATEND_LIST2(*Da,*De,avect);
    }
                #ifdef ModelP
    else if (DDD_InfoNCopies(PARHDR(VECT(avect))) > 0)
    {
      SETAVFINE(avect,1);
      SETAVTESTED(avect,1);
      SETAVSKIP(avect,1);
      ADDATEND_LIST2(*Da,*De,avect);
    }
                #endif
    else
    {
      i=STRONG_OUT(avect);
      ADDATEND_LIST2(Ua[i],Ue[i],avect);
    }
  }
  return(DONE);
}

INT CountStrongNeighbors(AVECTOR *initialS, DOUBLE *avNrOfStrongNbsHnd, INT *maxNeighbors)
{
  int nrOfPoints,sumOfStrongLinks;
  int nrOfNbs,nrOfStrongNbs;
  VECTOR *vect,*vect2;
  AVECTOR *avect,*avect2;
  MATRIX *mat;

  *avNrOfStrongNbsHnd=0.0;
  sumOfStrongLinks=nrOfPoints=*maxNeighbors=0;
  for (avect=initialS; avect!=NULL; avect=avect->succ)
  {
    nrOfPoints++;
    vect=VECT(avect);
    nrOfNbs=0;
    nrOfStrongNbs=0;
    for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
    {
      if (STRONG(mat))
      {
        sumOfStrongLinks++;
        vect2=MDEST(mat);
        avect2=(AVECTOR *) VISTART(vect2);
        STRONG_OUT(avect2)++;
        nrOfStrongNbs++;
      }
      nrOfNbs++;
    }
    if (nrOfNbs>*maxNeighbors) *maxNeighbors=nrOfNbs;
    STRONG_IN(avect)=nrOfStrongNbs;
  }

  *avNrOfStrongNbsHnd=((DOUBLE) sumOfStrongLinks)/((DOUBLE) nrOfPoints);

  return(DONE);
}

/****************************************************************************/
/*D
   l_vectorflags_consistent - make flags consistent

   SYNOPSIS:
   static INT l_vectorflags_consistent (GRID *g);

   PARAMETERS:
   .  g - pointer to grid

   DESCRIPTION:
   This function makes the coarsening flags consistent.

   RETURN VALUE:
   INT
   .n    NUM_OK      if ok
   .n    NUM_ERROR   if error occurrs
   D*/
/****************************************************************************/

#ifdef ModelP
static int Gather_VectorFlags (DDD_OBJ obj, void *data)
{
  VECTOR *pv = (VECTOR *)obj;
  INT *flag = (INT *)data;

  flag[0] = VCCOARSE(pv);

  PRINTDEBUG(np,3,("%d:gather  ind %3d gid %08x  c %d prio %d, attr %d\n",
                   me,VINDEX(pv),
                   DDD_InfoGlobalId(PARHDR(pv)),VCCOARSE(pv),
                   DDD_InfoPriority(PARHDR(pv)),
                   DDD_InfoAttr(PARHDR(pv))));

  return (0);
}

static int Scatter_VectorFlags (DDD_OBJ obj, void *data)
{
  VECTOR *pv = (VECTOR *)obj;
  INT *flag = (INT *)data;

  if (flag[0] == 1)
    SETVCCOARSE(pv,1);

  PRINTDEBUG(np,3,("%d:scatter ind %3d gid %08x  c %d\n",
                   me,VINDEX(pv),DDD_InfoGlobalId(PARHDR(pv)),VCCOARSE(pv)));

  return (0);
}

static INT l_vectorflag_consistent (GRID *g)
{
  PRINTDEBUG(np,3,("%d: l_vectorflag_consistent\n",me));

  DDD_IFAExchange(BorderVectorSymmIF, GRID_ATTR(g), sizeof(INT),
                  Gather_VectorFlags, Scatter_VectorFlags);
  return (NUM_OK);
}
#endif

/****************************************************************************/
/*                                                                          */
/* The following routines are for the coarsening by selection of CG points. */
/* It is essentially the algorithm described by Ruge-Stueben in 1987,       */
/* see also Neuss' thesis, ICA-Preprint 1996-07.                            */
/* CoarsenRugeStueben marks coarse grid points, then calls                  */
/* GenerateNewGrid to generate the new grid together with an                */
/* injection interpolation matrix (from each CG node to its father).        */
/*                                                                          */
/****************************************************************************/

static INT CheckImat(GRID *theGrid, int i)
{
  VECTOR *vect;
  MATRIX *imat;

  UserWriteF("Checking at point %d\n",i);
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
  {
    if (VECSKIP(vect)!=0)
      if (VISTART(vect)!=NULL)
        UserWrite("imat on Dirichlet node");

    for (imat=VISTART(vect); imat!=NULL; imat=MNEXT(imat))
    {
      if (MDEST(imat)==NULL)
      {
        if (VCCOARSE(vect))
          UserWrite("bad imat on coarse vect");
        else
          UserWrite("bad imat on fine vect");
      }
      else
      {
        UserWriteF("From ID(N(v)) %d to ID(N(v)) %d, value: %f",ID(VMYNODE(vect)),
                   ID(VMYNODE(MDEST(imat))),MVALUE(imat,0));
        if (CEXTRA(imat))
          UserWrite("  <extra>\n");
        else
          UserWrite("\n");
      }
    }
  }
  return(0);
}

static INT GenerateNewGrid(GRID *theGrid)
{
  INT noc,nof,m;
  VECTOR *vect,*newVect;
  GRID *newGrid;
  MULTIGRID *theMG;

        #ifdef ModelP
  l_vectorflag_consistent(theGrid);
        #endif

  /* if the new grid is all or empty we're done else generate it */
  nof=noc=0;
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
  {
    if (VCCOARSE(vect))
      noc++;
    else
      nof++;
  }
  m = noc * nof;
    #ifdef ModelP
  PRINTDEBUG(np,3,("%d: noc * nof %d\n",me,m));
        #ifdef Debug
  if (rep_err_count > 0) m = 0;
        #endif
  m = UG_GlobalMinINT(m);
  PRINTDEBUG(np,1,("%d: noc * nof %d\n",me,m));
        #endif
  if (m == 0)
    REP_ERR_RETURN(1);

  theMG=MYMG(theGrid);
  if ((newGrid=CreateNewLevelAMG(theMG))==NULL)
  {
    PrintErrorMessage('E',"GenerateNewGrid",
                      "could not create new amg level");
    REP_ERR_RETURN(1);
  }

    #ifdef ModelP
  DDD_IdentifyBegin();
        #endif

  /* generate vectors of newGrid */
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
  {
    assert(VISTART(vect)==NULL);
    if (VCCOARSE(vect))
    {
      /* generate coarse grid vector at end of list */
      if (CreateVector(newGrid,VOTYPE(vect),VOBJECT(vect),&newVect))
      {
        PrintErrorMessage('E',"GenerateNewGrid",
                          "could not create vector");
        REP_ERR_RETURN(1);
      }
      SETVCLASS(newVect,3);
      SETVNCLASS(newVect,VCLASS(vect));
      SETNEW_DEFECT(newVect,1);
      SETFINE_GRID_DOF(newVect,0);
      SETPRIO(newVect,PRIO(vect));
      VECSKIP(newVect)=VECSKIP(vect);

                        #ifdef ModelP
      if (DDD_InfoPrioCopies(PARHDR(vect)) > 0) {
        int *proclist = DDD_InfoProcList(PARHDR(vect));

        proclist += 2;
        PRINTDEBUG(np,3,("%d: ind %d gid %08x n %d ",me,VINDEX(vect),
                         DDD_InfoGlobalId(PARHDR(vect)),
                         DDD_InfoNCopies(PARHDR(vect))));
        while (*proclist != -1) {
          if (!GHOSTPRIO(proclist[1])) {
            PRINTDEBUG(np,3,("%d: pl %d\n",me,*proclist));
            DDD_IdentifyObject(PARHDR(newVect),*proclist,
                               PARHDR(vect));
          }
          proclist += 2;
        }
        PRINTDEBUG(np,3,(" prio %d, attr %d\n",
                         DDD_InfoPriority(PARHDR(newVect)),
                         DDD_InfoAttr(PARHDR(newVect))));

        PRINTDEBUG(np,3,("\n"));
      }
                        #endif

      /* an interpolation matrix is created ... */
      if (CreateIMatrix(theGrid,vect,newVect) == NULL)
      {
        PrintErrorMessage('E',"GenerateNewGrid",
                          "could not create interpolation matrix");
        REP_ERR_RETURN(1);
      }
      assert(VISTART(vect) != NULL);
      assert(MDEST(VISTART(vect)) != NULL);
      assert(VSTART(newVect) == NULL);
      if (CreateConnection(newGrid,newVect,newVect) == NULL)
      {
        PrintErrorMessage('E',"GenerateNewGrid",
                          "could not create diag matrix");
        REP_ERR_RETURN(1);
      }
      assert(VSTART(newVect) != NULL);
      assert(MDEST(VSTART(newVect)) == newVect);
    }
  }

    #ifdef ModelP
  DDD_IdentifyEnd();
        #endif

  PRINTDEBUG(np,1,("%d: IdentifyEnd %d\n",me,GLEVEL(theGrid)));

  return(DONE);
}

INT GeometricCoarsening(GRID *theGrid)
{
  VECTOR *vect,*cvect;
  NODE *theNode;

  if (GLEVEL(theGrid)<=0)
    REP_ERR_RETURN(1);

  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    SETVCCOARSE(vect,0);

  for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    if (CORNERTYPE(theNode))
    {
      vect = NVECTOR(theNode);
      cvect = NVECTOR((NODE *) NFATHER(theNode));
      SETVCCOARSE(vect,1);
      if (CreateIMatrix(theGrid,vect,cvect) == NULL)
        REP_ERR_RETURN(1);
    }

  return(0);
}

INT CoarsenRugeStueben(GRID *theGrid)
{
  int flag,i,k,maxNeighbors;
  INT error;
  DOUBLE avNosN;
  MULTIGRID *theMG;
  HEAP *theHeap;
  VECTOR *vect,*vect2,*vect3;
  AVECTOR *avect,*avect2,*avect3,*testCoarse;
  AVECTOR *initialS,*initialE,*Ca,*Ce,*Da,*De,*Fa,*Fe,*Ta,*Te;
  AVECTOR *Ua[2*MAXNEIGHBORS+1],*Ue[2*MAXNEIGHBORS+1];
  MATRIX *mat,*mat2,*mat3;
  INT MarkKey;

  theMG=MYMG(theGrid);

  /*	We now allocate some administrative copy of the fine grid, since we don't want
          to destroy any information upon it. This may be avoided for several
          special coarse grid choices. */

  theHeap=MGHEAP(theMG);
  MarkTmpMem(theHeap,&MarkKey);

  if ((error=SetupInitialList(theGrid,theHeap,&initialS,&initialE,MarkKey))!=DONE)
  {
    ReleaseTmpMem(theHeap,MarkKey);
    REP_ERR_RETURN(error);
  }

  if ((error=CountStrongNeighbors(initialS,&avNosN,&maxNeighbors))!=DONE)
  {
    ReleaseTmpMem(theHeap,MarkKey);
    REP_ERR_RETURN(error);
  }

  if (maxNeighbors>MAXNEIGHBORS)
  {
    ReleaseTmpMem(theHeap,MarkKey);
    REP_ERR_RETURN(1);
  }

  Ca=Ce=Fa=Fe=Ta=Te=Da=De=NULL;
  for (i=0; i<=2*maxNeighbors; i++)
    Ua[i]=Ue[i]=NULL;

  if ((error=DistributeInitialList(&initialS,&initialE,&Da,&De,Ua,Ue))!=DONE)
  {
    ReleaseTmpMem(theHeap,MarkKey);
    REP_ERR_RETURN(error);
  }

  /* first part: preliminary C-point choice */
  i=maxNeighbors;
  while (i>=0)
  {
    while ((avect=Ua[i])!=NULL)
    {
      ELIMINATE_LIST2(Ua[i],Ue[i],avect);
      ADDATEND_LIST2(Ca,Ce,avect);
      SETAVCOARSE(avect,1);
      vect=VECT(avect);
      for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
      {
        vect2=MDEST(mat);
        avect2=(AVECTOR *) VISTART(vect2);
        if (AVFINE(avect2)||AVCOARSE(avect2)) continue;
        mat2=MADJ(mat);
        if (mat2==NULL)
        {
          PrintErrorMessage('E',"CoarsenRugeStueben",
                            "G(A) is not symmetric");
          ReleaseTmpMem(theHeap,MarkKey);
          REP_ERR_RETURN(1);
        }
        if (STRONG(mat2))
        {
          k=STRONG_OUT(avect2);
          ELIMINATE_LIST2(Ua[k],Ue[k],avect2)
          ADDATEND_LIST2(Fa,Fe,avect2)
          SETAVFINE(avect2,1);
          for (mat3=MNEXT(VSTART(vect2)); mat3!=NULL; mat3=MNEXT(mat3))
            if (STRONG(mat3))
            {
              vect3=MDEST(mat3);
              avect3=(AVECTOR *) VISTART(vect3);
              if (AVFINE(avect3)||AVCOARSE(avect3)) continue;
              k=STRONG_OUT(avect3);
              ELIMINATE_LIST2(Ua[k],Ue[k],avect3)
              k++;
              if (k>i) i=k;
              STRONG_OUT(avect3)=k;
              ADDATEND_LIST2(Ua[k],Ue[k],avect3)
            }
        }
      }

      for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
        if (STRONG(mat))
        {
          vect2=MDEST(mat);
          avect2=(AVECTOR *) VISTART(vect2);
          if (AVFINE(avect2)||AVCOARSE(avect2)) continue;
          k=STRONG_OUT(avect2);
          ELIMINATE_LIST2(Ua[k],Ue[k],avect2)
          STRONG_OUT(avect2)=--k;
          ADDATEND_LIST2(Ua[k],Ue[k],avect2)
        }

    }
    i--;
  }

  /* second part: final C-point choice
     (tests if all F points i depend only on points j depending on C
          if not, either j  or i is made C-point) */
  while ((avect=Fa)!=NULL)
  {
    ELIMINATE_LIST2(Fa,Fe,avect)
    ADDATEND_LIST2(Ta,Te,avect)
    SETAVTESTED(avect,1);

    vect=VECT(avect);
    for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
    {
      vect2=MDEST(mat);
      avect2=(AVECTOR *) VISTART(vect2);
      if (AVCOARSE(avect2))
        SETVCUSED(vect2,1);                             /* is in coarse neighborhood of vect */
    }

    testCoarse=NULL;
    for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
      if (STRONG(mat))
      {
        vect2=MDEST(mat);
        if (VCUSED(vect2)) continue;

        flag=0;
        for (mat2=MNEXT(VSTART(vect2)); mat2!=NULL; mat2=MNEXT(mat2))
          if (STRONG(mat2))
            if (VCUSED(MDEST(mat2)))
            {
              flag=1; break;
            }

        if (flag==0)
        {
          if (testCoarse==NULL)
          {
            testCoarse=(AVECTOR *) VISTART(vect2);
            SETVCUSED(vect2,1);
          }
          else
          {
            testCoarse=avect;
            break;
          }
        }
      }

    if (testCoarse!=NULL)
    {
      if (AVTESTED(testCoarse))
        ELIMINATE_LIST2(Ta,Te,testCoarse)
        else
          ELIMINATE_LIST2(Fa,Fe,testCoarse)

          ADDATEND_LIST2(Ca,Ce,testCoarse)
          SETAVTESTED(testCoarse,0);
      SETAVFINE(testCoarse,0);
      assert(VECSKIP(VECT(testCoarse)) == 0);
      SETAVCOARSE(testCoarse,1);
    }

    /* set back used flag */
    for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
      SETVCUSED(MDEST(mat),0);
  }

  Fa=Ta; Fe=Te;

  /* copy coarse/fine information to theGrid
     and generate vectors of newGrid */
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
  {
    avect=(AVECTOR *) VISTART(vect);
    if (AVCOARSE(avect))
      SETVCCOARSE(vect,1);
    VISTART(vect)=NULL;             /* ==(AVECTOR *) VISTART(vect) */
  }
  error=GenerateNewGrid(theGrid);
  ReleaseTmpMem(theHeap,MarkKey);
  REP_ERR_RETURN(error);
}

INT CoarsenAverage (GRID *theGrid)
{
  INT error;
  FIFO myfifo;
  void *buffer;
  VECTOR *theV,*theW;
  NODE *theNode;
  MATRIX *theM;
  HEAP *theHeap;
  INT MarkKey;
  INT n,m,d,dmin,dmax;

  n = 0;
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV)) n++;
  dmin = n;
  PRINTDEBUG(np,3,("d "));
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV)) {
    if (VECSKIP(theV) == 0)
      VINDEX(theV) = -1;
    else {
      d = 0;
      for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM)) d++;
      VINDEX(theV) = d;
      dmin = MIN(d,dmin);
    }
  }
  PRINTDEBUG(np,3,("\ndmin %d\n",dmin));
  theHeap = MGHEAP(MYMG(theGrid));
  MarkTmpMem(theHeap,&MarkKey);
  buffer=(void *)GetTmpMem(theHeap,sizeof(VECTOR*)*n,MarkKey);
  if (buffer == NULL) {
    ReleaseTmpMem(theHeap,MarkKey);
    return(1);
  }
  fifo_init(&myfifo,buffer,sizeof(void *) * n);
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    if (VINDEX(theV) == dmin)
      fifo_in(&myfifo,(void *)theV);
  dmin++;
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    if (VINDEX(theV) == dmin)
      fifo_in(&myfifo,(void *)theV);
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    if (VINDEX(theV) > dmin)
      fifo_in(&myfifo,(void *)theV);
    #ifdef ModelP
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    if (VECSKIP(theV) == 0)
      if (DDD_InfoPrioCopies(PARHDR(theV)) > 0)
        fifo_in(&myfifo,(void *)theV);
        #endif
  dmax = 0;
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV)) {
    VINDEX(theV) = n;
    if (VECSKIP(theV) != 0) continue;
    if (VOTYPE(theV) != NODEVEC) continue;
    theNode = (NODE *) VOBJECT(theV);
    if (theNode == NULL) continue;
    if (OBJT(MYVERTEX(theNode)) == BVOBJ) {
      d = 0;
      for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM)) d++;
      VINDEX(theV) = d;
      dmax = MAX(d,dmax);
    }
  }
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    if (VINDEX(theV) == dmax)
      fifo_in(&myfifo,(void *)theV);
  dmax--;
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    if (VINDEX(theV) == dmax)
      fifo_in(&myfifo,(void *)theV);
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV))
    if (VINDEX(theV) < dmax)
      fifo_in(&myfifo,(void *)theV);
  m = 0;
  for (theV=FIRSTVECTOR(theGrid); theV!=NULL; theV=SUCCVC(theV)) {
    VINDEX(theV) = m++;
    SETVCUSED(theV,0);
  }
  while(!fifo_empty(&myfifo)) {
    theV = (VECTOR *)fifo_out(&myfifo);
    assert(VSTART(theV) != NULL);
    if (VCUSED(theV)) {
      if (VCCOARSE(theV) == 0) {
        for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM)) {
          theW = MDEST(theM);
          if (!VCUSED(theW)) break;
        }
        if (theM != NULL) {
          m--;
          SETVCCOARSE(theW,1);
          SETVCUSED(theW,1);
          fifo_in(&myfifo,(void *)theW);
          PRINTDEBUG(np,3,("out %d inc %d\n",
                           VINDEX(theV),VINDEX(theW)));
        }
      }
      else {
        PRINTDEBUG(np,3,("%d:outc %d",me,VINDEX(theV)));
        for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM)) {
          theW = MDEST(theM);
          if (VCUSED(theW)) continue;
          m--;
          SETVCUSED(theW,1);
          SETVCCOARSE(theW,0);
          fifo_in(&myfifo,(void *)theW);
          PRINTDEBUG(np,3,(" f %d",VINDEX(theW)));
        }
        PRINTDEBUG(np,3,("\n"));
      }
    }
    else {
      m--;
      SETVCCOARSE(theV,1);
      SETVCUSED(theV,1);
      PRINTDEBUG(np,3,("%d:c %d",me,VINDEX(theV)));
      for (theM=MNEXT(VSTART(theV)); theM!=NULL; theM=MNEXT(theM)) {
        theW = MDEST(theM);
        if (VCUSED(theW)) continue;
        m--;
        SETVCUSED(theW,1);
        SETVCCOARSE(theW,0);
        fifo_in(&myfifo,(void *)theW);
        PRINTDEBUG(np,3,(" f %d",VINDEX(theW)));
      }
      PRINTDEBUG(np,3,("\n"));
    }
  }
  ReleaseTmpMem(theHeap,MarkKey);
        #ifdef ModelP
  PRINTDEBUG(np,3,("%d: m %d\n",me,m));
  m = UG_GlobalMaxINT(m) - UG_GlobalMinINT(m);
  PRINTDEBUG(np,3,("%d: m %d\n",me,m));
        #endif
  if (m != 0)
    return(1);
  error = GenerateNewGrid(theGrid);
  return(error);
}

/****************************************************************************/
/*D
        CoarsenGreedy -

        DESCRIPTION:

        .vb
        .ve

        SEE ALSO:
        D*/
/****************************************************************************/

INT CoarsenGreedy (GRID *theGrid)
{
  VECTOR *vi,*vj;
  MATRIX *mij;
  INT nCoarse,nFine;
  INT error;

  /* set back used flag */
  for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
    SETVCUSED(vi,0);

  nCoarse = nFine = 0;
  for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
    if (!VCUSED(vi))
    {
      SETVCCOARSE(vi,1);
      SETVCUSED(vi,1);
      (nCoarse)++;
      for (mij=MNEXT(VSTART(vi)); mij!=NULL; mij=MNEXT(mij))
      {
        vj = MDEST(mij);
        if (!VCUSED(vj))
        {
          SETVCCOARSE(vj,0);
          SETVCUSED(vj,1);
          (nFine)++;
        }
      }
    }

  if ( (nFine+nCoarse) != NVEC(theGrid) )
    PrintErrorMessage('W',"CoarsenGreedy","not all vectors labeled!");

  /* label all Dirichlet fine:
      for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
              if ( VECSKIP(vi) )
              {
                      SETVCCOARSE(vi,0);
                      SETVCUSED(vi,1);
              } */

  error = GenerateNewGrid(theGrid);
  REP_ERR_RETURN(error);
}

INT CoarsenGreedyWithBndLoop (GRID *theGrid)
{
  VECTOR *vi,*vj;
  MATRIX *mij;
  INT minNeighbors,nNeighbors;
  DOUBLE vx, vy;
  char buffer[64];
  INT nCoarse,nFine;
  INT error;

  nCoarse = nFine = 0;
  /* set back used flag */
  for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
    SETVCUSED(vi,0);

  /* first loop over boundary vectors to find domain corners
          (supposing that vecs belonging to the corners have the minor number of neihgbors): */
  minNeighbors = NVEC(theGrid);
  for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
  {
    if (OBJT(MYVERTEX(VMYNODE(vi))) != BVOBJ) continue;
    nNeighbors = 0;
    for (mij=MNEXT(VSTART(vi)); mij!=NULL; mij=MNEXT(mij))
      nNeighbors++;
    minNeighbors = MIN(nNeighbors,minNeighbors);
  }

  sprintf(buffer," min no of conns: %d\n", minNeighbors);
  UserWrite(buffer);
  for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
    if (!VCUSED(vi) && (OBJT(MYVERTEX(VMYNODE(vi))) == BVOBJ) )
    {
      nNeighbors = 0;
      for (mij=MNEXT(VSTART(vi)); mij!=NULL; mij=MNEXT(mij))
        nNeighbors++;
      if (nNeighbors == minNeighbors)
      {
        vx = XC(MYVERTEX(VMYNODE(vi)));
        vy = YC(MYVERTEX(VMYNODE(vi)));
        sprintf(buffer," min no of conns at: x: %7.4f   y: %7.4f\n", vx, vy);
        UserWrite(buffer);
        SETVCCOARSE(vi,1);
        SETVCUSED(vi,1);
        (nCoarse)++;
        for (mij=MNEXT(VSTART(vi)); mij!=NULL; mij=MNEXT(mij))
        {
          vj = MDEST(mij);
          if (!VCUSED(vj) && (OBJT(MYVERTEX(VMYNODE(vj))) == BVOBJ) )
          {
            SETVCCOARSE(vj,0);
            SETVCUSED(vj,1);
            (nFine)++;
          }
        }
      }
    }

  /* second loop over boundary vectors: */
  for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
    if (!VCUSED(vi) &&  (OBJT(MYVERTEX(VMYNODE(vi))) == BVOBJ) )
    {
      SETVCCOARSE(vi,1);
      SETVCUSED(vi,1);
      (nCoarse)++;
      for (mij=MNEXT(VSTART(vi)); mij!=NULL; mij=MNEXT(mij))
      {
        vj = MDEST(mij);
        if (!VCUSED(vj) &&  (OBJT(MYVERTEX(VMYNODE(vj))) == BVOBJ) )
        {
          SETVCCOARSE(vj,0);
          SETVCUSED(vj,1);
          (nFine)++;
        }
      }
    }

  /* label the rest */
  for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
  {
    if (!VCUSED(vi))
    {
      SETVCCOARSE(vi,1);
      SETVCUSED(vi,1);
      (nCoarse)++;
      for (mij=MNEXT(VSTART(vi)); mij!=NULL; mij=MNEXT(mij))
      {
        vj = MDEST(mij);
        if (!VCUSED(vj))
        {
          SETVCCOARSE(vj,0);
          SETVCUSED(vj,1);
          (nFine)++;
        }
      }
    }
  }

  if ( (nFine+nCoarse) != NVEC(theGrid) )
    PrintErrorMessage('W',"CoarsenGreedy","not all vectors labeled!");

  /* set back used flag */
  for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
    SETVCUSED(vi,0);

  error = GenerateNewGrid(theGrid);
  REP_ERR_RETURN(error);
}

/****************************************************************************/
/*D
   bfs - do breadth first search labeling

   SYNOPSIS:
   static INT bfs (FIFO *Fifo, VECTOR *theSeedVector,
                                   long INT *nFine,long INT *nCoarse,
                                   long INT *nIsolated)

   PARAMETERS:
   .  Fifo - pointer to fifo data structure.
   .  nFine, nCoarse, nIsolated - number of vectors with according label

   DESCRIPTION:
   A vector vj in N_i (neighborhood of vector vi) is labeled coarse if no
   other vector vk in N_j is labeled coarse, otherwise vj is labeled fine.

   Only strong connections are considered in labeling, e.g. we use a
   reduced graph.

   NOTE: In "bfs" vectors are labeled BEFORE put in fifo! In other words:
   Vectors popped from fifo are always labeled.

   RETURN VALUE:
   INT
   .n    0      if ok
   .n    1      if fifo_in failed

   SEE ALSO:
   D*/
/****************************************************************************/

static INT bfs (FIFO *Fifo, VECTOR *theSeedVector,
                long INT *nFine,long INT *nCoarse,
                long INT *nIsolated)
{
  VECTOR *vi,*vj,*vk;
  MATRIX *mij,*mjk;
  INT no_coarse_neighbor;

  /* label seed vector */
  if (MNEXT(VSTART(theSeedVector))==NULL)
  {
    SETVCCOARSE(theSeedVector,0);
    (*nIsolated)++;
    (*nFine)++;
    return (0);
  }
  else         /* label seed vector coarse; put in queue: */
  {
    SETVCCOARSE(theSeedVector,1);
    (*nCoarse)++;
    if (fifo_in (Fifo, theSeedVector)==1) {
      PrintErrorMessage('E',"bfs", "fifo_in failed");
      UserWriteF(" used: %d, size: %d\n", Fifo->used, Fifo->size);
      return(1);
    }
  }
  SETVCUSED(theSeedVector,1);

  while (!fifo_empty (Fifo))
  {
    vi = fifo_out (Fifo);
#ifdef DebugAMG
    UserWriteF("pop  vector %d (node %d) from fifo\n", VINDEX(vi),ID(VMYNODE(vi)));
#endif
    /* run over all vectors vj adjacent to vi: */
    for (mij=MNEXT(VSTART(vi)); mij!=NULL; mij=MNEXT(mij))
    {
      vj = MDEST(mij);
      if (VCUSED(vj))
        continue;                                       /* j is always labeled (so: don't put in fifo again!)! */

      /* run over all vectors vk adjacent to vj: */
      no_coarse_neighbor = TRUE;
      for (mjk=MNEXT(VSTART(vj)); mjk!=NULL; mjk=MNEXT(mjk))
      {
        vk = MDEST(mjk);
        if ( VCCOARSE(vk) )
          if (STRONG(mjk))
          {
            no_coarse_neighbor = FALSE;
            break;
          }
      }

      if (no_coarse_neighbor)
      {
        SETVCCOARSE(vj,1);
        (*nCoarse)++;
#ifdef DebugAMG
        UserWriteF("    --> vj %d (node %d) coarse\n", VINDEX(vj),ID(VMYNODE(vj)));
#endif
      }
      else
      {
        SETVCCOARSE(vj,0);
        (*nFine)++;
#ifdef DebugAMG
        UserWriteF("    --> vj %d (node %d) fine\n", VINDEX(vj),ID(VMYNODE(vj)));
#endif
      }

      /* put in fifo: */
      SETVCUSED(vj,1);
#ifdef DebugAMG
      UserWriteF("push  vector vj %d (node %d) in fifo\n", VINDEX(vj),ID(VMYNODE(vj)));
#endif
      if (fifo_in (Fifo, vj)==1) {
        PrintErrorMessage('E',"bfs", "fifo_in failed");
        UserWriteF(" used: %d, size: %d\n", Fifo->used, Fifo->size);
        return(1);
      }
    }
  }       /* end while */

  return(0);
}

/****************************************************************************/
/*D
   CoarsenBreadthFirst - Breadth first search traversal(s) through
                         the grid graph

   SYNOPSIS:
   INT CoarsenBreadthFirst (GRID *theGrid)

   PARAMETERS:
   .  theGrid - pointer to grid.

   DESCRIPTION:
   Do breadth first search traversal(s) through the grid graph until
   all vectors are labeled. "CoarsenBreadthFirst" inits a fifo,
   controls the traversal(s) and calls "bfs" where the actual labeling
   is done.  Since the graph need not to be connected we call "bfs"
   nTraversals-times with a seed vector not yet labeled.

   RETURN VALUE:
   INT
   .n    0      if ok
   .n    1      if no memory for fifo or bfs failed
   .n    error  if GenerateNewGrid failed

   SEE ALSO:
   D*/
/****************************************************************************/

INT CoarsenBreadthFirst (GRID *theGrid)
{
  HEAP *theHeap;
  void *buffer;
  FIFO myFifo;
  long INT fifosize;
  INT MarkKey;

  VECTOR *vi, *theSeedVector;

  long INT nFine,nCoarse,nIsolated,nLabeled;
  long INT nFinetotal, nCoarsetotal, nIsolatedtotal, nLabeledtotal;
  INT nTraversals,error;

  nTraversals = nLabeledtotal = nFinetotal = nCoarsetotal = nIsolatedtotal = 0;

  /* set back used flag */
  for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
    SETVCUSED(vi,0);

  /* init fifo: */
  theHeap = MGHEAP(MYMG(theGrid));
  MarkTmpMem(theHeap,&MarkKey);
  fifosize = 2*NVEC(theGrid)*sizeof(VECTOR*);
  buffer=(void *)GetTmpMem(theHeap,fifosize,MarkKey);
  if (buffer == NULL)
  {
    PrintErrorMessage('E',"CoarsenBreadthFirst",
                      "could not get temp mem");
    ReleaseTmpMem(theHeap,MarkKey);
    REP_ERR_RETURN(1);
  }
  fifo_init(&myFifo,buffer,fifosize);

  while (nLabeledtotal < NVEC(theGrid))
  {
    for (theSeedVector=FIRSTVECTOR(theGrid); theSeedVector!=NULL; theSeedVector=SUCCVC(theSeedVector))
      if (!VCUSED(theSeedVector))
        break;
    if (theSeedVector==NULL) break;

    nFine = nCoarse = nIsolated = 0;
    if ((error = bfs (&myFifo,theSeedVector,&nFine,&nCoarse,&nIsolated))!=0)
    {
      PrintErrorMessage('E',"CoarsenBreadthFirst",
                        "bfs failed");
      REP_ERR_RETURN(1);
    }
    nTraversals++;
    nLabeled  = (nFine + nCoarse);             /* + nIsolated); now isolated are also fine! */
    IFDEBUG(np,4)
    UserWriteF("CoarsenBreadthFirst: %d vectors labeled by bfs in traversal %d(fine %d, coarse %d isolated %d)\n",nLabeled,nTraversals,nFine,nCoarse,nIsolated);
    ENDDEBUG

      nFinetotal        += nFine;
    nCoarsetotal      += nCoarse;
    nIsolatedtotal    += nIsolated;
    nLabeledtotal     += nLabeled;
  }

  /* clear fifo: */
  fifo_clear(&myFifo);
  ReleaseTmpMem(theHeap,MarkKey);

  /* label all Dirichlet fine: */
  for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
    if ( VECSKIP(vi) )
    {
      SETVCCOARSE(vi,0);
      SETVCUSED(vi,1);
      nFinetotal++;
    }

  error = GenerateNewGrid(theGrid);
  REP_ERR_RETURN(error);
}

/****************************************************************************/
/* The following routines are for the coarsening by clustering of points.   */
/* It is essentially the algorithm described by Vanek, et al (1994),        */
/* see also Neuss' thesis, ICA-Preprint 1996-07.                            */
/* PutVectInCluster, GenerateClusters are auxiliary functions to            */
/* CoarsenVanek. Observe that grid and IMATRIX-generation are not           */
/* completely seperated, since the clusters are defined by the IMATRIX      */
/* initially.                                                               */
/****************************************************************************/

static INT GenerateClusters(AVECTOR **Ua, AVECTOR **Ue, GRID *theGrid, GRID *newGrid, int minSizeOfCluster)
{
  int i,k,nc;
  VECTOR *vect,*vect2,*newVect;
  AVECTOR *avect,*avect2;
  MATRIX *mat,*mat2,*imat;
  AVECTOR *Ca,*Ce;

  if (minSizeOfCluster<0) minSizeOfCluster=0;

  i=MAXNEIGHBORS;
  while (i>=minSizeOfCluster)
  {
    while ((avect=Ua[i])!=NULL)
    {
      Ca=Ce=NULL;
      ELIMINATE_LIST2(Ua[i],Ue[i],avect);
      ADDATEND_LIST2(Ca,Ce,avect);
      vect=VECT(avect);
      SETVCCOARSE(vect,1);

      nc=1;
      for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
      {
        mat2=MADJ(mat);
        if (STRONG(mat2))
        {
          vect2=MDEST(mat);
          if (VCCOARSE(vect2)) continue;                                /* already belongs to a cluster */
          avect2=(AVECTOR *) VISTART(vect2);
          k=STRONG_OUT(avect2);
          ELIMINATE_LIST2(Ua[k],Ue[k],avect2);
          ADDATEND_LIST2(Ca,Ce,avect2);
          SETVCCOARSE(vect2,1);
          nc++;
        }
      }

      /* generate a coarse grid vector for the cluster */
      if (CreateVector(newGrid,VOTYPE(vect),VOBJECT(vect),&newVect))
      {
        PrintErrorMessage('E',"GenerateClusters",
                          "could not create vector");
        REP_ERR_RETURN(1);
      }
      SETVCLASS(newVect,3);
      SETVNCLASS(newVect,VCLASS(vect));
      SETNEW_DEFECT(newVect,1);
      SETFINE_GRID_DOF(newVect,0);
      VINDEX(newVect)=nc;
      VOBJECT(newVect) = VOBJECT(vect);

      /* generate the cluster as imatrices to newVect */
      for (avect=Ca; avect!=NULL; avect=avect->succ)
      {
        vect=VECT(avect);

        VISTART(vect)=NULL;                             /* !was used for connection
                                                           to avect up to now! */

        /* interpolation matrices are created to vect */
        if ((imat=CreateIMatrix(theGrid,vect,newVect))==NULL)
        {
          PrintErrorMessage('E',"GenerateClusters",
                            "could not create interpolation matrix");
          REP_ERR_RETURN(1);
        }

        /* change the order for the neighbors of the cluster */
        for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
          if (STRONG(mat))
          {
            vect2=MDEST(mat);
            if (VCCOARSE(vect2)) continue;
            avect2=(AVECTOR *) VISTART(vect2);
            k=STRONG_OUT(avect2);
            ELIMINATE_LIST2(Ua[k],Ue[k],avect2);
            STRONG_OUT(avect2)=--k;
            ADDATEND_LIST2(Ua[k],Ue[k],avect2);
          }
      }
    }
    i--;
  }

  return(DONE);
}

INT CoarsenVanek(GRID *theGrid)
{
  int i,k,minSize;
  INT error,maxNeighbors;
  DOUBLE avNosN;
  GRID *newGrid;
  MULTIGRID *theMG;
  HEAP *theHeap;
  VECTOR *vect,*vect2,*newVect,*newVect0;
  AVECTOR *avect,*avect2,*initialS,*initialE;
  AVECTOR *Da,*De,*Ua[2*MAXNEIGHBORS+1],*Ue[2*MAXNEIGHBORS+1];
  MATRIX *mat,*imat;
  INT MarkKey;

  theMG=MYMG(theGrid);
  theHeap=MGHEAP(theMG);
  MarkTmpMem(theHeap,&MarkKey);

  if ((error=SetupInitialList(theGrid,theHeap,&initialS,&initialE,MarkKey))!=DONE)
  {
    ReleaseTmpMem(theHeap,MarkKey);
    REP_ERR_RETURN(error);
  }

  if ((error=CountStrongNeighbors(initialS,&avNosN,&maxNeighbors))!=DONE)
  {
    ReleaseTmpMem(theHeap,MarkKey);
    REP_ERR_RETURN(error);
  }

  if (maxNeighbors>MAXNEIGHBORS)
  {
    PrintErrorMessage('E',"CoarsenVanek","too many neighbors");
    ReleaseTmpMem(theHeap,MarkKey);
    REP_ERR_RETURN(1);
  }

  if ((newGrid=CreateNewLevelAMG(theMG))==NULL)
  {
    PrintErrorMessage('E',"CoarsenCommand","could not create new amg level");
    ReleaseTmpMem(theHeap,MarkKey);
    REP_ERR_RETURN(1);
  }

  Da=De=NULL;
  for (i=0; i<=2*MAXNEIGHBORS; i++)
    Ua[i]=Ue[i]=NULL;

  if ((error=DistributeInitialList(&initialS,&initialE,&Da,&De,Ua,Ue))!=DONE)
  {
    ReleaseTmpMem(theHeap,MarkKey);
    REP_ERR_RETURN(error);
  }

  /* no interpolation to Dirichlet values */
  for (avect=Da; avect!=NULL; avect=avect->succ)
    VISTART(VECT(avect))=NULL;

  /* first step:
     we generate clusters with a size greater than (avNosN+1)*.66-1 */
  if ((error=GenerateClusters(Ua,Ue,theGrid,newGrid,
                              (INT)((avNosN+1.0)*0.66-1.0)))!=DONE)
  {
    ReleaseTmpMem(theHeap,MarkKey);
    REP_ERR_RETURN(error);
  }

  /* second step:
     try to add the rest of the points to an already existing cluster
     if there are several possibilities add it to the smallest cluster */
  i=0;
  while (i<MAXNEIGHBORS)
  {
    avect=Ua[i];
    while (avect!=NULL)
    {
      /* check if the point can be added to a cluster */
      vect=VECT(avect);
      minSize=999; newVect0=NULL;
      for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
        if (STRONG(mat))
        {
          vect2=MDEST(mat);
          if (VCCOARSE(vect2))
          {
            newVect=MDEST(VISTART(vect2));
            if (VINDEX(newVect)<minSize)
            {
              minSize=VINDEX(newVect);
              newVect0=newVect;
            }
          }
        }

      if (newVect0!=NULL)
      {
        SETVCCOARSE(vect,1);

        /* correct the order of all neighbors */
        for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
          if (STRONG(mat))
          {
            vect2=MDEST(mat);
            if (VCCOARSE(vect2)) continue;
            avect2=(AVECTOR *) VISTART(vect2);

            k=STRONG_OUT(avect2);
            ELIMINATE_LIST2(Ua[k],Ue[k],avect2);
            STRONG_OUT(avect2)=--k;
            ADDATEND_LIST2(Ua[k],Ue[k],avect2);
          }

        /* we now eliminate avect
           (but its successor is still valid!) */
        ELIMINATE_LIST2(Ua[i],Ue[i],avect);

        VISTART(vect)=NULL;
        if ((imat=CreateIMatrix(theGrid,vect,newVect0))==NULL)
        {
          PrintErrorMessage('E',"CoarsenVanek",
                            "could not create interpolation matrix");
          ReleaseTmpMem(theHeap,MarkKey);
          REP_ERR_RETURN(1);
        }

        VINDEX(newVect0)++;
      }
      avect=avect->succ;
    }

    i++;
  }

  /* third pass: go down to clusters of size zero */
  if ((error=GenerateClusters(Ua,Ue,theGrid,newGrid,0))!=DONE)
  {
    ReleaseTmpMem(theHeap,MarkKey);
    REP_ERR_RETURN(error);
  }

  /* we now don't need the heap anymore ... */
  ReleaseTmpMem(theHeap,MarkKey);

  return(DONE);
}


/****************************************************************************/
/*                                                                          */
/* Function:  IpRoutines                                                    */
/*                                                                          */
/* Purpose:   Generates a good interpolation/restriction for amg            */
/*            starting from a crude one (injection or piecewise constant)   */
/*            (changes for systems necessary!!)                             */
/*                                                                          */
/****************************************************************************/

static DOUBLE Dist (VECTOR *v, VECTOR *w)
{
  DOUBLE_VECTOR a,b;
  DOUBLE s;

  VectorPosition(v,a);
  VectorPosition(w,b);

  V_DIM_EUKLIDNORM_OF_DIFF(a,b,s);

  return(s);
}

INT IpAverage (GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I)
{
  INT ncomp,i,j,n,nmax;
  DOUBLE s;
  GRID *newGrid;
  VECTOR *vect,*dest,*newVect;
  MATRIX *mat,*imat;

  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    if (VCCOARSE(vect)) {
      assert(VISTART(vect)!=NULL);
      newVect = MDEST(VISTART(vect));
      assert(newVect!=NULL);
      VECSKIP(newVect) = VECSKIP(vect);
    }

  newGrid=theGrid->coarser;
  nmax = 1 + NVEC(theGrid) / NVEC(newGrid);
  PRINTDEBUG(np,3,("%d:Ip nmax %d\n",me,nmax));
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    if (VCCOARSE(vect) == 0) {
      for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
        SETMUSED(mat,0);
      ncomp = MD_COLS_IN_RT_CT(A,VTYPE(vect),VTYPE(vect));
      n = 0;
            #ifdef ModelP
      if (DDD_InfoPrioCopies(PARHDR(vect)) > 0) {
        for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
          if (DDD_InfoPrioCopies(PARHDR(MDEST(mat))) > 0)
            if (VCCOARSE(MDEST(mat))==1) {
              SETMUSED(mat,1);
              n++;
            }
      }
      else
            #endif
      if (VECSKIP(vect)) {
        for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
          if (VECSKIP(MDEST(mat)))
            if (VCCOARSE(MDEST(mat))==1) {
              SETMUSED(mat,1);
              n++;
            }
      }
      if (n == 0)
        for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat)) {
          if (VCCOARSE(MDEST(mat))==1) {
            SETMUSED(mat,1);
            n++;
          }
          if (n >= nmax) break;
        }
      PRINTDEBUG(np,2,("%d:%d ->",VINDEX(vect),n));
      if (n == 0) {
        PrintErrorMessage('E',"IpAverage",
                          "can't build Interpolation matrix (n = 0)");
        REP_ERR_RETURN(1);
      }
      /*
         sum = 0.0;
         for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat)) {
          if (!MUSED(mat)) continue;
          dest = MDEST(mat);
              sum += Dist(vect,dest);
         }
         sum = 1.0 / sum;
       */
      for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat)) {
        if (!MUSED(mat)) continue;
        dest = MDEST(mat);
        PRINTDEBUG(np,3,(" %d",VINDEX(dest)));
        if (MD_COLS_IN_RT_CT(A,VTYPE(vect),VTYPE(dest)) != ncomp) {
          PrintErrorMessage('E',"IpAverage",
                            "can't handle this format");
          REP_ERR_RETURN(1);
        }
        assert(VISTART(dest)!=NULL);
        newVect = MDEST(VISTART(dest));
        assert(newVect!=NULL);
        if ((imat=CreateIMatrix(theGrid,vect,newVect))==NULL) {
          PrintErrorMessage('E',"IpAverage",
                            "could not create interpolation matrix");
          REP_ERR_RETURN(1);
        }
        SETMDIAG(imat,1);
        /* s = sum * Dist(vect,dest); */
        s = 1.0/((DOUBLE) n);
        for (i=0; i<ncomp; i++)
          for (j=0; j<ncomp; j++) {
            if (i == j) MVALUE(imat,i*ncomp+j) = s;
            else MVALUE(imat,i*ncomp+j) = 0.0;
          }
      }
      PRINTDEBUG(np,3,("\n"));
    }
    else {
      ncomp = MD_COLS_IN_RT_CT(A,VTYPE(vect),VTYPE(vect));
      imat = VISTART(vect);
      assert (imat != NULL);
      SETMDIAG(imat,1);
      for (i=0; i<ncomp; i++)
        for (j=0; j<ncomp; j++) {
          if (i == j) MVALUE(imat,i*ncomp+j) = 1.0;
          else MVALUE(imat,i*ncomp+j) = 0.0;
        }
    }
  PRINTDEBUG(np,3,("\n"));

  IFDEBUG(np,2)
  CheckImat(theGrid,2);
  ENDDEBUG

  return(DONE);
}


INT IpRugeStueben(GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I)
{
  INT icomp,mcomp,vskip,nCoarse;
  DOUBLE modDiag[MAX_MAT_COMP],modDiagInv[MAX_MAT_COMP];
  DOUBLE sum[MAX_MAT_COMP],sumInv[MAX_MAT_COMP],factor[MAX_MAT_COMP];
  GRID *newGrid;
  VECTOR *vect,*vect2,*vect3,*newVect;
  MATRIX *mat,*mat2,*imat,*imat2;

  /* Check that only mtype=0x0 (for other, AMG is unclear anyhow).        */
  /* Define and set the variables scalar,blockN,blockNN, and error. */
  BLOCK_SETUP(A);
  if (error) {
    PrintErrorMessage('E',"IpRugeStueben",
                      "not yet for general matrices");
    REP_ERR_RETURN(error);
  }

  icomp = 0;       /* preliminary, later this should be obtained via I */
  mcomp=MD_MCMP_OF_MTYPE(A,0,0);

  newGrid=theGrid->coarser;

  /*    For the fine grid points we have to solve for e_i:
          0 == a_{ii} e_i + \sum_{j \in C_i} a_{ij} e_j +
               \sum_{j \in W_i-C} a_{ij} e_i + \sum_{j \in S_i-C} a_{ij} e_j
          where C_i:coarse, S_i:strong, W_i:weak in N_i
          e_j in the last term is approximated by
               \frac{{\sum_{C_j \cut C_i} a_{jk} e_k}}{\sum{C_j \cut C_i} a_{jk}}
      Please note that the IP-matrices directly connecting
          the coarse grid points to their fathers
          are used here to store intermediate values of a local computation. */

  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    if (VCCOARSE(vect)==0)
    {
      vskip=VECSKIP(vect);
      nCoarse=0;
      if (vskip==0)
      {
        mat=VSTART(vect);
        BLOCK_COPY(&(MVALUE(mat,mcomp)),modDiag);
        for (mat=MNEXT(mat); mat!=NULL; mat=MNEXT(mat))
        {
          vect2=MDEST(mat);
          assert(VCUSED(vect2)==0);
          if (VCCOARSE(vect2)==1)
          {
            SETVCUSED(vect2,1);                                 /* is coarse in neighborhood of vect */
            BLOCK_COPY(&(MVALUE(mat,mcomp)),&(MVALUE(VISTART(vect2),icomp)));
            nCoarse++;
          }
          else
          {
            /* since with weakly connected points the below
               interpolation to coarse grid points
               is not by the coarse grid choice guaranteed to work,
               these are simply lumped to the diagonal */
            if ((STRONG(mat)==0)&&(VECSKIP(vect2)==0))
              BLOCK_ADD1(&(MVALUE(mat,mcomp)),modDiag);
          }
        }

#ifdef DebugAMG
        UserWriteF("NID=%d: s=%f, ",ID(VMYNODE(vect)),s);
#endif

        for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
          if (STRONG(mat))                       /* no strong influence on vskip!=0 nodes! */
          {
            vect2=MDEST(mat);
            if (VCCOARSE(vect2)==0)
            {
              /* compute \sum_{l \in C_i} a_{jl} */
              BLOCK_CLEAR(sum);
              for (mat2=MNEXT(VSTART(vect2)); mat2!=NULL; mat2=MNEXT(mat2))
              {
                vect3=MDEST(mat2);
                if (VCUSED(vect3))
                  BLOCK_ADD1(&(MVALUE(mat2,mcomp)),sum);
              }
#ifdef DebugAMG
              UserWriteF("(%d)%f, ",ID(VMYNODE(vect2)),sum);
#endif
              BLOCK_INVERT(sum,sumInv);
              BLOCK_MUL(&(MVALUE(mat,mcomp)),sumInv,factor);
              for (mat2=MNEXT(VSTART(vect2)); mat2!=NULL; mat2=MNEXT(mat2))
              {
                vect3=MDEST(mat2);
                if (VCUSED(vect3))
                  BLOCK_MUL_ADD(factor, &(MVALUE(mat2,mcomp)),
                                &(MVALUE(VISTART(vect3),icomp)));
              }
            }
          }
#ifdef DebugAMG
        UserWrite("\n");
#endif
        BLOCK_INVERT(modDiag, modDiagInv);
        BLOCK_SCALE1(-1.0, modDiagInv);
      }

      for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
      {
        vect2=MDEST(mat);
        if (VCUSED(vect2))
        {
          SETVCUSED(vect2,0);

          imat=VISTART(vect2);
          newVect=MDEST(imat);
          if ((imat2=CreateIMatrix(theGrid,vect,newVect))==NULL)
          {
            PrintErrorMessage('E',"IpRugeStueben","could not create interpolation matrix");
            REP_ERR_RETURN(1);
          }
          /* in Average: {BLOCK_SCALIDENTITY(1/((DOUBLE) nCoarse), &(MVALUE(imat2,icomp)));} */
          BLOCK_MUL_NNN(modDiagInv, &(MVALUE(imat,icomp)), &(MVALUE(imat2,icomp)));
        }
      }
    }

  /* Finally, we set imat to identity on the direct coarse grid fathers */
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    if (VCCOARSE(vect))
      BLOCK_IDENTITY(&(MVALUE(VISTART(vect),icomp)));

  IFDEBUG(np,4)
  CheckImat(theGrid,2);
  ENDDEBUG

  return(DONE);
}


/****************************************************************************/
/*D
   IpSchurMG - Create interpolation matrices and store rest/prol weights

   SYNOPSIS:
   static INT IpSchurMG(GRID *theGrid, MATDATA_DESC *A, INT RWflag)

   PARAMETERS:
   .  theGrid - pointer to grid.
   .  A - matrix descriptor for stiffness matrix.
   .  RWflag - flag denoting type of lumping: REUSKEN - simple lump; WAGNER - matdep lump


   DESCRIPTION:
   We first compute entries in a so called "lumped" matrix in the following way:
   For each fine vector vi go through its neighborhood N_i and find its
   fine neighbors j, then in the same way find j's coarse neighbors k.
   Replace in the i's equation these fine j's by a mean (arithmetic for "REUSKEN",
   matrix-coefficent weighted for the "WAGNER" case) of the coarse k's.
   This corresponds to new entries appearing in the fine-coarse block of the stiffness
   matrix which are not actually built.

   According to the entries in the lumped matrix we then create the interpolation
   matrices, and compute and store restriction/prolongation weights.

   We build a reduced matrix before lumping:
   Weak fine-fine connections A_ij are lumped (sorry, a second kind of lumping) to A_ii,
   e.g., such a weakly coupled fine vector vj in N_i is NOT substituted by the described
   (arithmetic or matrix dependent) mean of coarse vectors vk in N_j!

   The restriction operator is build due to (~ means "lumped")
   .vb
   restr = [-A_cf*Inv(~A_ff), I].
   .ve
   In components:
   .vb
   [restr]_ji = -[A_cf]_ji*[Inv(~A_ff)]_ii.
   .ve

   The prolongation operator is build due to
   .vb
   prol = [-Inv(~A_ff)*(~A_fc), I]^T.
   .ve
   In components:
   .vb
   [prol]_ij = -[Inv(~A_ff)]_ii*[~A_fc]_ij.
   .ve

   NOTE: restriction contains NO additional fine-coarse connections!

   The RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured (neg. VINDEX if CreateIMatrix failed).
   D*/
/****************************************************************************/

static INT IpSchurMG(GRID *theGrid, MATDATA_DESC *A, INT RWflag)
{
  INT icomp,mcomp,rcomp,nCoarse;
  INT arithmetic_lump,vsmask;
  DOUBLE Sum_jk[MAX_MAT_COMP], Weight_ik[MAX_MAT_COMP];
  DOUBLE InvSum_jk[MAX_MAT_COMP], Inv_ii[MAX_MAT_COMP];
  GRID *newGrid;
  VECTOR *vi,*vj,*vk,*newVect;
  MATRIX *mij,*mjk,*imat,*imat2;

  DOUBLE cut = 1.0e-3;       /* preliminary, perhaps later this should be obtained via np */

  /* Check that only mtype=0x0 (for other, AMG is unclear anyhow).        */
  /* Define and set the variables scalar, blockN, blockNN, and error.     */
  BLOCK_SETUP(A);
  if (error) {
    PrintErrorMessage('E',"IpSchurMG",
                      "not yet for general matrices");
    REP_ERR_RETURN(error);
  }

  icomp = 0;       /* preliminary, later this should be obtained via I?? */
  rcomp = blockNN;
  mcomp=MD_MCMP_OF_MTYPE(A,0,0);

  newGrid=theGrid->coarser;

  /*    For the fine grid points we have to solve for e_i:
          0 == A_{ii} e_i + \sum_{j \in C_i} A_{ij} e_j +
               \sum_{j \in N_i-C_i} A_{ij} e_j
          where N_i: neighborhood of node i, C_i:coarse in N_i
          e_j in the last term is approximated by
               \frac{{\sum_{k \in C_j} e_k}}{\sum{k \in C_j} 1}
          in the case of the Reusken variant and
               \frac{{\sum_{k \in C_j} A_{jk} e_k}}{\sum{k \in C_j} A_{jk}}
          for Wagners variant.
      Please note that the IP-matrices directly connecting
          the coarse grid points to their fathers
          are used here to store intermediate values of a local computation. */

  for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
  {
    SETVCUSED(vi,0);
    if (VCCOARSE(vi))
      BLOCK_CLEAR( &(MVALUE(VISTART(vi),icomp)) );
  }

  vsmask=(1<<blockN)-1;
  for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
    if ((VCCOARSE(vi)==0))
    {
      /* skip, if SKIP flags are set for all components */
      if ((VECSKIP(vi) & vsmask) == vsmask)
        continue;

#ifdef DebugAMG
      UserWriteF("NID=%d: A_ii[%d]=%f, ",ID(VMYNODE(vi)),mcomp,MVALUE(VSTART(vi),mcomp));
#endif

      /* compute lumped fine-coarse block ~A_fc and
         store result in direct Imats */
      for (mij=MNEXT(VSTART(vi)); mij!=NULL; mij=MNEXT(mij))
      {
        vj=MDEST(mij);
        if (VCCOARSE(vj)==1)                         /* original entry in A_fc */
        {BLOCK_ADD1( &(MVALUE(mij,mcomp)), &(MVALUE(VISTART(vj),icomp)) );}
        else
        {
          if (STRONG(mij))
          {
            /* do lumping: compute \sum_{k \in C_j} A_{jk} or nb of coarse neighbors */
            BLOCK_CLEAR(Sum_jk);

            nCoarse = 0;
            for (mjk=MNEXT(VSTART(vj)); mjk!=NULL; mjk=MNEXT(mjk))
            {
              vk=MDEST(mjk);
              if (VCCOARSE(vk))
              {
                nCoarse++;
                if (RWflag==WAGNER)
                  BLOCK_ADD1( &(MVALUE(mjk,mcomp)), Sum_jk );
              }
            }

            if (RWflag==WAGNER)
              BLOCK_INVERT(Sum_jk,InvSum_jk);

            /* use arithmetic average for REUSKEN, if BLOCK_INVERT returned error!=0,
               or if Wagner's cut criterion is fulfilled */
            if ( (RWflag==REUSKEN) || (error!=0) ||
                 ((RWflag==WAGNER) && (fabs(Sum_jk[0]) <= cut*MVALUE(VSTART(vj),mcomp))) )
              arithmetic_lump = 1;
            else
              arithmetic_lump = 0;
#ifdef DebugAMG
            UserWriteF("Sum_jk[0] (%d)%f, ",ID(VMYNODE(vj)),Sum_jk[0]);
#endif

            if ( arithmetic_lump )
            {
              /* Weight_ik = Aff_ij/nCoarse: */
              BLOCK_SCALE( 1.0 /((DOUBLE) nCoarse), &(MVALUE(mij,mcomp)), Weight_ik );
            }
            else
              /* Weight_ik = Aff_ij*InvSum_jk (in this order!): */
              BLOCK_MUL( &(MVALUE(mij,mcomp)), InvSum_jk, Weight_ik );

            /* lumped fine-coarse block ~A_fc*/
            for (mjk=MNEXT(VSTART(vj)); mjk!=NULL; mjk=MNEXT(mjk))
            {
              vk=MDEST(mjk);
              if (VCCOARSE(vk))
              {
                if ( arithmetic_lump )
                {BLOCK_ADD1( Weight_ik, &(MVALUE(VISTART(vk),icomp)) );}
                else
                  BLOCK_MUL_ADD( Weight_ik, &(MVALUE(mjk,mcomp)), &(MVALUE(VISTART(vk),icomp)) );
              }
            }
          }                                  /* end if STRONG(mij) */
          else
          {
            /* lump weak fine-fine connections A_ij on A_ii */
            BLOCK_ADD1( &(MVALUE(mij,mcomp)),&(MVALUE(VSTART(vi),mcomp)) );
          }
        }                                /* end is fine */
      }
#ifdef DebugAMG
      UserWrite("\n");
#endif

      BLOCK_INVERT( &(MVALUE(VSTART(vi),mcomp)), Inv_ii );
      if (error)
      {
        PrintErrorMessage('E',"IpSchurMG","inversion of Aff_ii failed!");
        BLOCK_WRITEOUT( &(MVALUE(VSTART(vi),mcomp)) );
        UserWriteF("    vi %d, level %d\n", VINDEX(vi), GLEVEL(theGrid));
      }
      BLOCK_SCALE1( -1.0, Inv_ii );
      /* now we create the interpolation matrices */
      for (mij=MNEXT(VSTART(vi)); mij!=NULL; mij=MNEXT(mij))
      {
        vj=MDEST(mij);
        if (VCCOARSE(vj)==1)
        {
          imat=VISTART(vj);
          newVect=MDEST(imat);
          /* [prol]_ij = -[Inv(~A_ff)]_ii*[~A_fc]_ij
             [restr]_ji = -[A_cf]_ji*[Inv(~A_ff)]_ii */
          if ((imat2=CreateIMatrix(theGrid,vi,newVect))==NULL)
          {
            PrintErrorMessage('E',"IpSchurMG","could not create interpolation matrix");
            REP_ERR_RETURN(1);
          }
          BLOCK_MUL_ADD_NNT( Inv_ii, &(MVALUE(imat,icomp)), &(MVALUE(imat2,icomp)) );
          BLOCK_MUL_ADD_NNT( &(MVALUE(MADJ(mij),mcomp)), Inv_ii, &(MVALUE(imat2,rcomp)) );
          SETVCUSED(vj,1);
        }
      }

      /* once again to handle additional connections */
      for (mij=MNEXT(VSTART(vi)); mij!=NULL; mij=MNEXT(mij))
        if (STRONG(mij))
        {
          vj=MDEST(mij);
          if ( VCCOARSE(vj)==0 )
          {
            for (mjk=MNEXT(VSTART(vj)); mjk!=NULL; mjk=MNEXT(mjk))
            {
              vk=MDEST(mjk);
              if ( (VCCOARSE(vk)==1) && (VCUSED(vk)==0) )
              {
                imat=VISTART(vk);
                newVect=MDEST(imat);
                if ((imat2=CreateIMatrix(theGrid,vi,newVect))==NULL)
                {
                  PrintErrorMessage('E',"IpSchurMG","could not create interpolation matrix");
                  REP_ERR_RETURN(1);
                }
                BLOCK_MUL_ADD_NNT( Inv_ii, &(MVALUE(imat,icomp)), &(MVALUE(imat2,icomp)) );
                BLOCK_CLEAR( &(MVALUE(imat2,rcomp)) );                                                 /* restriction contains no extra conn's! */
                SETVCUSED(vk,1);
              }
            }
          }
        }

      /* set back flags */
      for (mij=MNEXT(VSTART(vi)); mij!=NULL; mij=MNEXT(mij))
      {
        vj=MDEST(mij);
        if (VCCOARSE(vj)==1)
        {
          if (VCUSED(vj)==1)
          {
            imat=VISTART(vj);
            BLOCK_CLEAR( &(MVALUE(imat,icomp)) );
            SETVCUSED(vj,0);
          }
        }
        else                        /* fine */
        {
          for (mjk=MNEXT(VSTART(vj)); mjk!=NULL; mjk=MNEXT(mjk))
          {
            vk=MDEST(mjk);
            if (VCUSED(vk)==1)
            {
              imat=VISTART(vk);
              BLOCK_CLEAR( &(MVALUE(imat,icomp)) );
              SETVCUSED(vk,0);
            }
          }
        }
      }

    }

  /* Finally, we set imat to identity on the direct coarse grid fathers */
  for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
    if (VCCOARSE(vi))
    {
      BLOCK_IDENTITY( &(MVALUE(VISTART(vi),icomp)) );
      BLOCK_IDENTITY( &(MVALUE(VISTART(vi),rcomp)) );
    }

  IFDEBUG(np,4)
  CheckImat(theGrid,2);
  ENDDEBUG

  return(DONE);
}

INT IpReusken(GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I)
{
  return (IpSchurMG(theGrid,A,REUSKEN));
}

INT IpWagner(GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I)
{
  return (IpSchurMG(theGrid,A,WAGNER));
}

/****************************************************************************/
/* The next routines are for the SchurMG algorithm                          */
/* It is unclear if their effect can be obtained also by a smoothing step   */
/****************************************************************************/

/****************************************************************************/
/*D
   NBFineGridCorrection - Calculate Schur complement fine grid correction

   SYNOPSIS:
   INT NBFineGridCorrection (GRID *theGrid, const VECDATA_DESC *to,
                             const VECDATA_DESC *from)

   PARAMETERS:
   .  theGrid - pointer to grid

   DESCRIPTION:
   .vb
   c_f += Inv(~A_ff)*d_f; c: correction, d: defect (matdep: transformed defect).
   .ve

   For Matrix dependent lumping we suppose that the defect is always transformed
   before the fine grid correction step!

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured (negative vector index if inversion of small block failed).
   D*/
/****************************************************************************/

INT NBFineGridCorrection (GRID *theGrid, const VECDATA_DESC *to,const VECDATA_DESC *from,
                          const MATDATA_DESC *A)
{
  VECTOR *vi;
  DOUBLE Inv_ii[MAX_MAT_COMP];

  INT vsmask;                           /* vec skip mask        */
  register SHORT mcomp;                 /* mat-component(s)     */
  register SHORT tvcomp,fvcomp;         /* #vector-component(s) */

  /* Check that only mtype=0x0 (for other, AMG is unclear anyhow).        */
  /* Define and set the variables scalar, blockN, blockNN, and error. */
  BLOCK_SETUP(A);
  if (error) {
    PrintErrorMessage('E',"NBFineGridCorrection", "not yet for general matrices");
    REP_ERR_RETURN(error);
  }
  mcomp  = MD_MCMP_OF_MTYPE(A,0,0);
  tvcomp = VD_CMP_OF_TYPE(to,0,0);
  fvcomp = VD_CMP_OF_TYPE(from,0,0);

  vsmask=(1<<blockN)-1;
  for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
    if (VCCOARSE(vi)==0)
    {
      if ((VECSKIP(vi) & vsmask) == vsmask)
        continue;
      BLOCK_INVERT( &(MVALUE(VSTART(vi),mcomp)), Inv_ii );
      if ( error )
      {
        PrintErrorMessage('E',"NBFineGridCorrection","inversion of Aff_ii failed!");
        BLOCK_WRITEOUT( &(MVALUE(VSTART(vi),mcomp)) );
        UserWriteF("   vi %d\n", VINDEX(vi));
        return (1);
      }
      else
        BLOCK_MATVECADD( Inv_ii, &(VVALUE(vi,fvcomp)),&(VVALUE(vi,tvcomp)) );
    }
  return (NUM_OK);

}

/****************************************************************************/
/*D
        NBTransformDefect - Transform defect in Schur complement mg

   SYNOPSIS:
   INT NBTransformDefect (GRID *theGrid, const VECDATA_DESC *to,
   const VECDATA_DESC *from, const MATDATA_DESC *Mat)

   PARAMETERS:
   .  theGrid - pointer to grid
   .  to, from - locations in vector data field
   .  Mat - locations in stiffness resp. lumped matrix

   DESCRIPTION:
   In Schur complement mg with matrix dependent lumping the defect in the fine vectors has to be
   transformed according to
   .vb
   (F*d)_f = {2*I - A_ff*Inv(~A_ff)}*d_f,
   .ve
   The defect in the coarse vectors remains unchanged.
   .vb
   (F*d)_c =  d_c.
   .ve
   In components:
   .vb
   [(F*d)_f]_i = 2*[d_f]_i - Sum_j{[A_ff]_ij * [Inv(~A_ff)]_jj * [d_f]_j
   .ve

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    1 if error occured (negative vector index if inversion of small block failed).
   D*/
/****************************************************************************/

INT NBTransformDefect (GRID *theGrid, const VECDATA_DESC *to, const VECDATA_DESC *from,
                       const MATDATA_DESC *A)
{
  VECTOR *vi, *vj;
  MATRIX *mij;

  DOUBLE Inv_jj[MAX_MAT_COMP],Weight_ij[MAX_MAT_COMP];
  DOUBLE sum_j[MAX_VEC_COMP];

  INT vsmask;                           /* vec skip mask       */
  register SHORT mcomp;                 /* mat-component(s)    */
  register SHORT tvcomp,fvcomp;         /* vector-component(s) */

  /* Check that only mtype=0x0 (for other, AMG is unclear anyhow).        */
  /* Define and set the variables scalar, blockN ,blockNN, and error. */
  BLOCK_SETUP(A);
  if (error) {
    PrintErrorMessage('E',"NBTransformDefect", "not yet for general matrices");
    REP_ERR_RETURN(error);
  }
  mcomp  = MD_MCMP_OF_MTYPE(A,0,0);
  tvcomp = VD_CMP_OF_TYPE(to,0,0);
  fvcomp = VD_CMP_OF_TYPE(from,0,0);

  vsmask=(1<<blockN)-1;
  for (vi=FIRSTVECTOR(theGrid); vi!=NULL; vi=SUCCVC(vi))
    if (VCCOARSE(vi)==0)
    {
      if ((VECSKIP(vi) & vsmask) == vsmask)
        continue;

      /*sum_j = 0.0;*/
      BLOCK_VECCLEAR(sum_j)
      for (mij=VSTART(vi); mij!=NULL; mij=MNEXT(mij))
      {
        vj = MDEST(mij);
        if ( VCCOARSE(vj)==0 )
        {
          BLOCK_INVERT( &(MVALUE(VSTART(vj),mcomp)), Inv_jj );
          if ( error )
          {
            PrintErrorMessage('E',"NBTransformDefect","inversion of Aff_jj failed!");
            BLOCK_WRITEOUT( &(MVALUE(VSTART(vj),mcomp)) );
            UserWriteF("    vi %d --> vj %d\n", VINDEX(vi),VINDEX(vj));
            return (NUM_ERROR);
          }
          /* sum_j = Sum_{j} (Matff_ij * InvLumpff_jj * defect_j) */
          BLOCK_MUL( &(MVALUE(mij,mcomp)),Inv_jj,Weight_ij );
          BLOCK_MATVECADD( Weight_ij, &(VVALUE(vj,fvcomp)),sum_j)
        }
      }
      /* ~defect_i = 2.0*defect_i - sum_j */
      BLOCK_VECSCALE( 2.0, &(VVALUE(vi,fvcomp)), &(VVALUE(vi,tvcomp)) );
      BLOCK_VECSUB1( sum_j, &(VVALUE(vi,tvcomp)) );

    }

  return (NUM_OK);
}

/****************************************************************************/
/* Here start the interpolation routines of ClusterAMG type                 */
/****************************************************************************/

INT IpPiecewiseConstant(GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I)
{
  INT i,j,ncomp;
  VECTOR *vect;
  MATRIX *imat;

  /* we set imat to piecewise constant on clusters */
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    if ((imat=VISTART(vect))!=NULL)
    {
      ncomp = MD_COLS_IN_RT_CT(A,VTYPE(vect),VTYPE(vect));
      SETMDIAG(imat,1);
      for (i=0; i<ncomp; i++)
        for (j=0; j<ncomp; j++)
          if (i == j) MVALUE(imat,i*ncomp + j) = 1.0;
          else MVALUE(imat,i*ncomp + j) = 0.0;

    }

  return(DONE);
}

INT IpVanek(GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I)
{
  INT icomp,mcomp;
  DOUBLE sum[MAX_MAT_COMP],factor[MAX_MAT_COMP];
  VECTOR *vect,*vect2,*newVect2;
  MATRIX *mat,*imat,*imat2,*imats;

  /* Check that only mtype=0x0 (for other, AMG is unclear anyhow).        */
  /* Define and set the variables scalar,blockN,blockNN, and error. */
  BLOCK_SETUP(A);
  if (error) {
    PrintErrorMessage('E',"IpVanek",
                      "not yet for general matrices");
    REP_ERR_RETURN(error);
  }

  icomp = 0;       /* preliminary, later this should be obtained via I */
  mcomp=MD_MCMP_OF_MTYPE(A,0,0);

  /* first we set imat to piecewise constant on clusters */
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    if ((imat=VISTART(vect))!=NULL)
      BLOCK_IDENTITY(&(MVALUE(imat,icomp)));

  /* then this prolongation is smoothed */
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
  {
    if (VECSKIP(vect)!=0) continue;

    /* compute in s the diagonal of the filtered matrix */
    mat=VSTART(vect);
    BLOCK_COPY(&(MVALUE(mat,mcomp)), sum);
    for (mat=MNEXT(mat); mat!=NULL; mat=MNEXT(mat))
    {
      vect2=MDEST(mat);
      /* !!!! Vorzeichenfehler in [VMB]  !!!! */
      if ((STRONG(mat)==0)&&(VECSKIP(vect2)==0))
        BLOCK_ADD1(&(MVALUE(mat,mcomp)), sum);
    }
    BLOCK_INVERT(sum, factor);
    BLOCK_SCALE1(-0.666666666, factor);

    imat=VISTART(vect); imats=MNEXT(imat);
    for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
      if (STRONG(mat))
      {
        vect2=MDEST(mat);                         /* can't be Dirichlet, since it is strongly connected */
        newVect2=MDEST(VISTART(vect2));
        if ((imat2=GetIMatrix(vect,newVect2))==NULL)
        {
          if ((imat2=CreateIMatrix(theGrid,vect,newVect2))==NULL)
          {
            PrintErrorMessage('E',"IpVanek","could not create interpolation matrix");
            REP_ERR_RETURN(1);
          }

          /* but the cluster imat should remain the first one */
          MNEXT(imat2)=imats; MNEXT(imat)=imat2; VISTART(vect)=imat; imats=imat2;
          BLOCK_CLEAR(&(MVALUE(imat2,icomp)));
        }
        BLOCK_MUL_ADD_NNT(factor, &(MVALUE(mat,mcomp)), &(MVALUE(imat2,icomp)));
      }
  }

  IFDEBUG(np,4)
  CheckImat(theGrid,2);
  ENDDEBUG

  return(DONE);
}

/****************************************************************************/
/*                                                                          */
/* Function:  FastGalerkinFromInterpolation                                 */
/*                                                                          */
/* Purpose:   Given a prolongation this routine computes a Galerkin         */
/*            coarse grid matrix.                                           */
/*            (changes for systems necessary!!)                             */
/*                                                                          */
/****************************************************************************/


INT FastGalerkinFromInterpolation(GRID *theGrid, MATDATA_DESC *A,
                                  MATDATA_DESC *I, INT type)
{
  INT icomp,mcomp,rcomp;
  DOUBLE R_ikA_kl[MAX_MAT_COMP];
  GRID *coarseGrid;
  register VECTOR *cvi,*cvj;
  VECTOR *vk,*vl;
  register MATRIX *cmij,*plj;
  MATRIX *mkl,*rik;
  INT symmetric, rinj, pinj;

  /* Check that only mtype=0x0 (for other, AMG is unclear anyhow).        */
  /* Define and set the variables scalar,blockN,blockNN, and error. */
  BLOCK_SETUP(A);
  if (error) {
    PrintErrorMessage('E',"FastGalerkinFromInterpolation",
                      "not yet for general matrices, use AssembleGalerkinFromInterpolation");
    REP_ERR_RETURN(error);
  }

  symmetric = rinj = pinj = 0;
  if (type & 1) symmetric = 1;
  if (type & 2) rinj = 1;
  if (type & 4) pinj = 1;
  if (type & 8) rcomp = blockNN;else rcomp = 0;
  icomp = 0;
  mcomp = MD_MCMP_OF_MTYPE(A,0,0);


  coarseGrid=theGrid->coarser;

  /* even if this should not be necessary for newly generated AMG grids */
  for (cvi=FIRSTVECTOR(coarseGrid); cvi!=NULL; cvi=SUCCVC(cvi))
    if (VISTART(cvi)!=NULL)
    {
      PrintErrorMessage('E',"FastGalerkinFromInterpolation","VISTART not empty on coarse grid");
      REP_ERR_RETURN(NUM_ERROR);
    }

  /* this sums up all R_ik*A_kl*P_lj */
  for (vk=FIRSTVECTOR(theGrid); vk!=NULL; vk=SUCCVC(vk))
  {
    if (rinj)
      if (VCCOARSE(vk) == 0) continue;
    for (rik=VISTART(vk); rik!=NULL; rik=MNEXT(rik))
    {
      cvi=MDEST(rik);

      /* to keep access to matrices fast we store the addresses of
              coarse grid matrices cvi->cvj in VISTART(cvj)!!! */
      for (cmij=VSTART(cvi); cmij!=NULL; cmij=MNEXT(cmij))
        VISTART(MDEST(cmij))=cmij;

      for (mkl=VSTART(vk); mkl!=NULL; mkl=MNEXT(mkl))
      {
        vl=MDEST(mkl);
        if (rinj)
        {BLOCK_COPY(&(MVALUE(mkl,mcomp)), R_ikA_kl);}
        else
        {BLOCK_MUL_NNN(&(MVALUE(rik,rcomp)), &(MVALUE(mkl,mcomp)), R_ikA_kl);}

        for (plj=VISTART(vl); plj!=NULL; plj=MNEXT(plj))
        {
          cvj=MDEST(plj);
          if ((cmij=VISTART(cvj))==NULL)
          {
            if ((cmij=CreateExtraConnection(coarseGrid,cvi,cvj))==NULL)
            {
              PrintErrorMessage('E',"FastGalerkinFromInterpolation","could not create stiffness matrix");
              REP_ERR_RETURN(NUM_ERROR);
            }
            BLOCK_CLEAR(&(MVALUE(cmij,mcomp)));
            BLOCK_CLEAR(&(MVALUE(MADJ(cmij),mcomp)));

            /* and keep VISTART up to date! */
            VISTART(cvj)=cmij;
          }
          /*	if (pinj) {MVALUE(cmij,mcomp)+=R_ikA_kl;break;} else */
          BLOCK_MUL_ADD_NTN(R_ikA_kl,&(MVALUE(plj,icomp)),&(MVALUE(cmij,mcomp)));
        }
      }

      /* now it is necessary to clear the VISTART fields again!! */
      for (cmij=VSTART(cvi); cmij!=NULL; cmij=MNEXT(cmij))
        VISTART(MDEST(cmij))=NULL;
    }
  }

        #ifdef ModelP
  if (l_ghostmatrix_collect(coarseGrid,A))
    return(NUM_ERROR);
        #endif

  return(NUM_OK);
}

INT AssembleGalerkinFromInterpolation(GRID *theGrid, MATDATA_DESC *A,
                                      MATDATA_DESC *I, INT type)
{
  return(AssembleGalerkinByMatrix(theGrid,A,type));
}
/****************************************************************************/
/*                                                                          */
/* Function:  SparsenCGMatrix                                               */
/*                                                                          */
/* Purpose:   Sometimes it may be useful to sparsen the computed            */
/*            coarse grid matrix. This routine assumes that the             */
/*            matrices chosen to be kept are marked as strong               */
/*            by the routines from the beginning.                           */
/*                                                                          */
/****************************************************************************/

INT SparsenCGMatrix(GRID *theGrid, MATDATA_DESC *A, INT lumpFlag)
{
  INT mcomp;
  VECTOR *vect;
  MATRIX *mat,*matD,*matS;

  /* Check that only mtype=0x0 (for other, AMG is unclear anyhow).        */
  /* Define and set the variables scalar,blockN,blockNN,mcomp, and error. */
  BLOCK_SETUP(A);
  if (error) {
    PrintErrorMessage('E',"SparsenCGMatrix", "not yet for general matrices");
    REP_ERR_RETURN(error);
  }

  mcomp=MD_MCMP_OF_MTYPE(A,0,0);


  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
  {
    matD=VSTART(vect);
    for (mat=MNEXT(matD); mat!=NULL; mat=matS)
    {
      matS=MNEXT(mat);
      if ((STRONG(mat)==0)&&(STRONG(MADJ(mat))==0))
      {
        if (lumpFlag)
          BLOCK_ADD1( &(MVALUE(mat,mcomp)), &(MVALUE(matD,mcomp)) )

          if (DisposeConnection(theGrid,MMYCON(mat))!=0)
          {
            PrintErrorMessage('E',"SparsenCGMatrix",
                              "could not dispose connection");
            REP_ERR_RETURN(1);
          }
      }
    }
  }

  return(DONE);
}

/****************************************************************************/
/*                                                                          */
/* Function:  ReorderFineGrid                                               */
/*                                                                          */
/* Purpose:   The ordering of the points often is important.                */
/*            eg. RUGESTUEBEN requires COARSEFINE whereas                   */
/*            WAGNER/REUSKEN requires FINECOARSE.                           */
/*                                                                          */
/****************************************************************************/

INT ReorderFineGrid(GRID *theGrid, INT orderType)
{
    #ifndef ModelP
  VECTOR *vect,*CaV,*CeV,*FaV,*FeV,*TaV,*TeV;

  switch (orderType)
  {
  case ASBEFORE : break;
  case COARSEFINE :
  case FINECOARSE :
    TaV=TeV=CaV=CeV=FaV=FeV=NULL;
    while ((vect=FIRSTVECTOR(theGrid))!=NULL)
    {
      ELIMINATE_LIST2(FIRSTVECTOR(theGrid),LASTVECTOR(theGrid),vect);
      if (VECSKIP(vect))
      {
        ADDATEND_LIST2(TaV,TeV,vect);
        continue;
      }

      if (VCCOARSE(vect)==1)
      {
        ADDATEND_LIST2(CaV,CeV,vect);
        continue;
      }
      else
      {
        ADDATEND_LIST2(FaV,FeV,vect);
        continue;
      }
    }

    if (orderType==COARSEFINE)
    {
      APPEND_LIST2(FIRSTVECTOR(theGrid),LASTVECTOR(theGrid),CaV,CeV);
      APPEND_LIST2(FIRSTVECTOR(theGrid),LASTVECTOR(theGrid),FaV,FeV);
      APPEND_LIST2(FIRSTVECTOR(theGrid),LASTVECTOR(theGrid),TaV,TeV);
    }
    else
    {
      APPEND_LIST2(FIRSTVECTOR(theGrid),LASTVECTOR(theGrid),FaV,FeV);
      APPEND_LIST2(FIRSTVECTOR(theGrid),LASTVECTOR(theGrid),CaV,CeV);
      APPEND_LIST2(FIRSTVECTOR(theGrid),LASTVECTOR(theGrid),TaV,TeV);
    }
    break;
  }
        #endif

  return(DONE);
}
