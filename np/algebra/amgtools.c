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

#undef DebugAMG

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

INT UnmarkAll(GRID *theGrid, MATDATA_DESC *A, DOUBLE theta)
{
  VECTOR *vect;
  MATRIX *mat;

  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    for (mat=VSTART(vect); mat!=NULL; mat=MNEXT(mat))
      SETSTRONG(mat,0);
  return(0);
}

INT MarkAll(GRID *theGrid, MATDATA_DESC *A, DOUBLE theta)
{
  VECTOR *vect;
  MATRIX *mat;

  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    for (mat=VSTART(vect); mat!=NULL; mat=MNEXT(mat))
      SETSTRONG(mat,1);
  return(0);
}

INT MarkOffDiagWithoutDirichlet(GRID *theGrid, MATDATA_DESC *A, DOUBLE theta)
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

INT MarkAbsolute(GRID *theGrid, MATDATA_DESC *A, DOUBLE theta)
{
  VECTOR *vect;
  MATRIX *mat;
  INT mcomp;

  if (MD_IS_SCALAR(A)==FALSE)
    REP_ERR_RETURN(1);
  mcomp=MD_SCALCMP(A);

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

INT MarkRelative(GRID *theGrid, MATDATA_DESC *A, DOUBLE theta)
{
  VECTOR *vect,*vect2;
  MATRIX *matD,*mat;
  INT mcomp;
  DOUBLE s,threshold;

  if (MD_IS_SCALAR(A)==FALSE)
    REP_ERR_RETURN(1);
  mcomp=MD_SCALCMP(A);

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

/****************************************************************************/
/* Some routines marking strong connections  (!change for systems!)         */
/****************************************************************************/

INT SetupInitialList(GRID *theGrid, HEAP *theHeap, AVECTOR **initialSH, AVECTOR **initialEH)
{
  VECTOR *vect;
  AVECTOR *avect;

  /* sets index field and copies vectors to initial list */
  *initialSH=*initialEH=NULL;
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
  {
    SETVCUSED(vect,0);
    SETVCCOARSE(vect,0);

    if ((avect=(AVECTOR *) GetMem(theHeap,sizeof(AVECTOR),FROM_TOP))==NULL)
    {
      PrintErrorMessage('E',"SetupInitialList",
                        "could not allocate avector");
      REP_ERR_RETURN(1);
    }

    CTRL(avect)=0;
    SETAVCOARSE(avect,0);
    SETAVFINE(avect,0);
    SETAVTESTED(avect,0);
    STRONG_IN(avect)=0;
    STRONG_OUT(avect)=0;

    VECT(avect)=vect;
    VISTART(vect)=(MATRIX *) avect;             /* VISTART field is used for establishing bijection */

    ADDATEND_LIST2(*initialSH,*initialEH,avect);
  }

  return(DONE);
}

INT DistributeInitialList(AVECTOR **La, AVECTOR **Le, AVECTOR **Ta, AVECTOR **Te, AVECTOR **Ua, AVECTOR **Ue)
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
      ADDATEND_LIST2(*Ta,*Te,avect);
    }
                #ifdef ModelP
    else if (DDD_InfoNCopies(PARHDR(VECT(avect))) > 0)
    {

      printf("border skiped\n");

      SETAVFINE(avect,1);
      SETAVTESTED(avect,1);
      ADDATEND_LIST2(*Ta,*Te,avect);
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
  INT noc,nof;
  VECTOR *vect,*newVect;
  GRID *newGrid;
  MULTIGRID *theMG;
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
  INT noc,nof;
  VECTOR *vect,*newVect;
  GRID *newGrid;
  MULTIGRID *theMG;

  /* if the new grid is all or empty we're done else generate it */
  nof=noc=0;
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
  {
    if (VCCOARSE(vect))
      noc++;
    else
      nof++;
  }
  if (noc*nof==0)
    return(DONE);

  theMG=MYMG(theGrid);
  if ((newGrid=CreateNewLevelAMG(theMG))==NULL)
  {
    PrintErrorMessage('E',"GenerateNewGrid",
                      "could not create new amg level");
    REP_ERR_RETURN(1);
  }

  /* generate vectors of newGrid */
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
  {
    assert(VISTART(vect)==NULL);
    if (VCCOARSE(vect))
    {
      /* generate coarse grid vector at end of list */
      if ((newVect=CreateVector(newGrid,VTYPE(vect),vect->object))==NULL)
      {
        PrintErrorMessage('E',"GenerateNewGrid",
                          "could not create vector");
        REP_ERR_RETURN(1);
      }
      SETVCLASS(newVect,3);
      SETVNCLASS(newVect,VCLASS(vect));

      /* an interpolation matrix is created ... */
      if (CreateIMatrix(theGrid,vect,newVect) == NULL)
      {
        PrintErrorMessage('E',"GenerateNewGrid",
                          "could not create interpolation matrix");
        REP_ERR_RETURN(1);
      }
      assert(VISTART(vect)!=NULL);
      assert(MDEST(VISTART(vect))!=NULL);
    }
  }
  return(DONE);
}


INT CoarsenRugeStueben(GRID *theGrid)
{
  int flag,i,k,maxNeighbors;
  INT error;
  DOUBLE avNosN;
  MULTIGRID *theMG;
  HEAP *theHeap;
  VECTOR *vect,*vect2,*vect3;
  AVECTOR *avect,*avect2,*avect3,*testCoarse,*initialS,*initialE,*Ca,*Ce,*Fa,*Fe,*Ta,*Te;
  AVECTOR *Ua[2*MAXNEIGHBORS+1],*Ue[2*MAXNEIGHBORS+1];
  MATRIX *mat,*mat2,*mat3;

  theMG=MYMG(theGrid);

  /*	We now allocate some administrative copy of the fine grid, since we don't want
          to destroy any information upon it. This may be avoided for several
          special coarse grid choices. */

  theHeap=MGHEAP(theMG);
  Mark(theHeap,FROM_TOP);

  if ((error=SetupInitialList(theGrid,theHeap,&initialS,&initialE))!=DONE)
  {
    Release(theHeap,FROM_TOP);
    REP_ERR_RETURN(error);
  }

  if ((error=CountStrongNeighbors(initialS,&avNosN,&maxNeighbors))!=DONE)
  {
    Release(theHeap,FROM_TOP);
    REP_ERR_RETURN(error);
  }

  if (maxNeighbors>MAXNEIGHBORS)
  {
    Release(theHeap,FROM_TOP);
    REP_ERR_RETURN(1);
  }

  Ca=Ce=Fa=Fe=Ta=Te=NULL;
  for (i=0; i<=2*maxNeighbors; i++)
    Ua[i]=Ue[i]=NULL;

  if ((error=DistributeInitialList(&initialS,&initialE,&Ta,&Te,Ua,Ue))!=DONE)
  {
    Release(theHeap,FROM_TOP);
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
      assert(VECSKIP(VECT(avect)) == 0);
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
          Release(theHeap,FROM_TOP);
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
  Release(theHeap,FROM_TOP);
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
  INT error;
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
      if ((newVect=CreateVector(newGrid,VTYPE(vect),vect->object))==NULL)
      {
        PrintErrorMessage('E',"GenerateClusters","could not create vector");
        REP_ERR_RETURN(1);
      }
      SETVCLASS(newVect,3);
      SETVNCLASS(newVect,CLASS(vect));
      INDEX(newVect)=nc;

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

  theMG=MYMG(theGrid);
  theHeap=MGHEAP(theMG);
  Mark(theHeap,FROM_TOP);

  if ((error=SetupInitialList(theGrid,theHeap,&initialS,&initialE))!=DONE)
  {
    Release(theHeap,FROM_TOP);
    REP_ERR_RETURN(error);
  }

  if ((error=CountStrongNeighbors(initialS,&avNosN,&maxNeighbors))!=DONE)
  {
    Release(theHeap,FROM_TOP);
    REP_ERR_RETURN(error);
  }

  if (maxNeighbors>MAXNEIGHBORS)
  {
    PrintErrorMessage('E',"CoarsenVanek","too many neighbors");
    Release(theHeap,FROM_TOP);
    REP_ERR_RETURN(1);
  }

  if ((newGrid=CreateNewLevelAMG(theMG))==NULL)
  {
    PrintErrorMessage('E',"CoarsenCommand","could not create new amg level");
    Release(theHeap,FROM_TOP);
    REP_ERR_RETURN(1);
  }

  Da=De=NULL;
  for (i=0; i<=2*MAXNEIGHBORS; i++)
    Ua[i]=Ue[i]=NULL;

  if ((error=DistributeInitialList(&initialS,&initialE,&Da,&De,Ua,Ue))!=DONE)
  {
    Release(theHeap,FROM_TOP);
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
    Release(theHeap,FROM_TOP);
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
            if (INDEX(newVect)<minSize)
            {
              minSize=INDEX(newVect);
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
          Release(theHeap,FROM_TOP);
          REP_ERR_RETURN(1);
        }

        INDEX(newVect0)++;
      }
      avect=avect->succ;
    }

    i++;
  }

  /* third pass: go down to clusters of size zero */
  if ((error=GenerateClusters(Ua,Ue,theGrid,newGrid,0))!=DONE)
  {
    Release(theHeap,FROM_TOP);
    REP_ERR_RETURN(error);
  }

  /* we now don't need the heap anymore ... */
  Release(theHeap,FROM_TOP);

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

INT IpAverage (GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I)
{
  INT icomp,mcomp,ncomp,i,j,n;
  DOUBLE s,t,sum,factor;
  GRID *newGrid;
  VECTOR *vect,*dest,*newVect;
  MATRIX *mat,*mat2,*imat;

  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    if (VCCOARSE(vect)) {
      dest = vect;
      assert(VECSKIP(dest)==0);
      assert(VISTART(dest)!=NULL);
      newVect = MDEST(VISTART(dest));
      assert(newVect!=NULL);
    }

  newGrid=theGrid->coarser;
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    if ((VCCOARSE(vect)==0)&&(VECSKIP(vect)==0))
    {
      ncomp = MD_COLS_IN_RT_CT(A,VTYPE(vect),VTYPE(vect));
      n = 0;
      for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
        if (VCCOARSE(MDEST(mat))==1) n++;
      assert(n > 0);
      s = 1.0 / n;
      for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat)) {
        dest = MDEST(mat);
        if (VCCOARSE(dest) == 0) continue;
        if (MD_COLS_IN_RT_CT(A,VTYPE(vect),VTYPE(dest)) != ncomp) {
          PrintErrorMessage('E',"IpAverage","can't handle this format");
          REP_ERR_RETURN(1);
        }
        assert(VISTART(dest)!=NULL);
        newVect = MDEST(VISTART(dest));
        assert(newVect!=NULL);
        if ((imat=CreateIMatrix(theGrid,vect,newVect))==NULL) {
          PrintErrorMessage('E',"IpAverge",
                            "could not create interpolation matrix");
          REP_ERR_RETURN(1);
        }
        SETMDIAG(imat,1);
        MVALUE(imat,0) = s;
        for (i=0; i<ncomp; i++)
          for (j=0; j<ncomp; j++)
            if (i == j) MVALUE(imat,i*ncomp + j) = s;
            else MVALUE(imat,i*ncomp + j) = 0.0;
      }
    }
    else {
      imat = VISTART(vect);
      if (imat == NULL) continue;
      SETMDIAG(imat,1);
      MVALUE(imat,0) = 1.0;
      for (i=0; i<ncomp; i++)
        for (j=0; j<ncomp; j++)
          if (i == j) MVALUE(imat,i*ncomp + j) = 1.0;
          else MVALUE(imat,i*ncomp + j) = 0.0;
    }

  IFDEBUG(np,4)
  CheckImat(theGrid,2);
  ENDDEBUG

  return(DONE);
}


INT IpRugeStueben(GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I)
{
  INT icomp,mcomp;
  DOUBLE s,t,sum,factor;
  GRID *newGrid;
  VECTOR *vect,*vect2,*vect3,*newVect;
  MATRIX *mat,*mat2,*imat;

  if (MD_IS_SCALAR(A)==FALSE)
    REP_ERR_RETURN(1);
  mcomp=MD_SCALCMP(A);

  /* if (MD_IS_SCALAR(I)==FALSE)
          REP_ERR_RETURN(1);
     icomp=MD_SCALCMP(I); */
  icomp = 0;       /* preliminary, later this should be obtained via I */

  newGrid=theGrid->coarser;

  /*    For the fine grid points we have to solve for e_i:
          0 == a_{ii} e_i + \sum_{j \in C_i} a_{ij} e_j +
               \sum_{j \in W_i-C} a_{ij} e_i + \sum_{j \in S_i-C} a_{ij} e_j
          where C_i:coarse, S_i:strong, W_i:weak in N_i
          The last term is approximated by
               \frac{\sum{C_j \cut C_i} a_{jk}}{{\sum_{C_j \cut C_i} a_{jk} e_k}}
      Please note that the IP-matrices directly connecting
          the coarse grid points to their fathers
          are used here to store intermediate values of a local computation. */

  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    if ((VCCOARSE(vect)==0)&&(VECSKIP(vect)==0))
    {
      mat=VSTART(vect);
      s=MVALUE(mat,mcomp);
      for (mat=MNEXT(mat); mat!=NULL; mat=MNEXT(mat))
      {
        vect2=MDEST(mat);
        assert(VCUSED(vect2)==0);
        if (VCCOARSE(vect2)==1)
        {
          SETVCUSED(vect2,1);                           /* is coarse in neighborhood of vect */
          MVALUE(VISTART(vect2),icomp)=MVALUE(mat,mcomp);
        }
        else
        {
          /* since with weakly connected points the below
             interpolation to coarse grid points
             is not by the coarse grid choice guaranteed to work,
             these are simply lumped to the diagonal */
          if ((STRONG(mat)==0)&&(VECSKIP(vect2)==0))
            s+=MVALUE(mat,mcomp);
        }
      }

#ifdef DebugAMG
      UserWriteF("NID=%d: s=%f, ",ID(VMYNODE(vect)),s);
#endif
      /*******unfortunately we have to store s to formulate
         WAGNER for later use in amgc:
         in the previous version this was done in FACTOREDMAT.
         For the new version one
         should  look at this point again.

         switch (ipType)
         {
              case REUSKEN:
              case WAGNER:
                      VVALUE(vect,FACTOREDMAT)=s;
                      break;
              default:
                      break;
         }********/

      for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
        if (STRONG(mat))
        {
          vect2=MDEST(mat); if (VCCOARSE(vect2)==0)
          {
            /* compute \sum_{l \in C_i} a_{jl} */
            sum=0.0;
            for (mat2=MNEXT(VSTART(vect2)); mat2!=NULL; mat2=MNEXT(mat2))
            {
              vect3=MDEST(mat2);
              if (VCUSED(vect3))
                sum+=MVALUE(mat2,mcomp);
            }
#ifdef DebugAMG
            UserWriteF("(%d)%f, ",ID(VMYNODE(vect2)),sum);
#endif
            factor=MVALUE(mat,mcomp)/sum;
            for (mat2=MNEXT(VSTART(vect2)); mat2!=NULL; mat2=MNEXT(mat2))
            {
              vect3=MDEST(mat2);
              if (VCUSED(vect3))
                MVALUE(VISTART(vect3),icomp)+=MVALUE(mat2,mcomp)*factor;
            }
          }
        }
#ifdef DebugAMG
      UserWrite("\n");
#endif
      for (mat=MNEXT(VSTART(vect)); mat!=NULL; mat=MNEXT(mat))
      {
        vect2=MDEST(mat);
        if (VCUSED(vect2))
        {
          SETVCUSED(vect2,0);

          imat=VISTART(vect2);
          t=-MVALUE(imat,icomp)/s;
          newVect=MDEST(imat);
          if ((imat=CreateIMatrix(theGrid,vect,newVect))==NULL)
          {
            PrintErrorMessage('E',"IpRugeStueben","could not create interpolation matrix");
            REP_ERR_RETURN(1);
          }
          MVALUE(imat,icomp)=t;

          /************** for Reusken/Wagner the restriction is different [r= - L_{21} D_{11}^{-1}]
             ... perhaps later ...
             switch (ipType)
             {
                  case REUSKEN:
                  case WAGNER:
                  MVALUE(imat,Rij)=-VALUE(GetMatrix(vect2,vect),mcomp)/s;
                  break;

                  default:
                  MVALUE(matR,Rij)=t;
                  break;
             } ****************/
        }
      }
    }

  /* Finally, we set imat to identity on the direct coarse grid fathers */
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    if (VCCOARSE(vect))
      MVALUE(VISTART(vect),icomp)=1.0;

  IFDEBUG(np,4)
  CheckImat(theGrid,2);
  ENDDEBUG

  return(DONE);
}


INT IpVanek(GRID *theGrid, MATDATA_DESC *A, MATDATA_DESC *I)
{
  INT icomp,mcomp;
  DOUBLE s,factor;
  VECTOR *vect,*vect2,*newVect2;
  MATRIX *mat,*imat,*imat2,*imats;

  if (MD_IS_SCALAR(A)==FALSE)
    REP_ERR_RETURN(1);
  mcomp=MD_SCALCMP(A);

  /* if (MD_IS_SCALAR(I)==FALSE)
          REP_ERR_RETURN(1);
     icomp=MD_SCALCMP(I); */
  icomp = 0;       /* preliminary, later this should be obtained via I */

  /* first we set imat to piecewise constant on clusters */
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
    if ((imat=VISTART(vect))!=NULL)
      MVALUE(imat,icomp)=1.0;

  /* then this prolongation is smoothed */
  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
  {
    if (VECSKIP(vect)!=0) continue;

    /* compute in s the diagonal of the filtered matrix */
    mat=VSTART(vect);
    s=MVALUE(mat,mcomp);
    for (mat=MNEXT(mat); mat!=NULL; mat=MNEXT(mat))
    {
      vect2=MDEST(mat);
      if ((STRONG(mat)==0)&&(VECSKIP(vect2)==0))
        s-=MVALUE(mat,mcomp);
    }
    factor=0.666666666/s;

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
          MVALUE(imat2,icomp)=0.0;
        }
        MVALUE(imat2,icomp)-=factor*MVALUE(mat,mcomp);
      }
  }

  IFDEBUG(np,4)
  CheckImat(theGrid,2);
  ENDDEBUG

  return(DONE);
}

/****************************************************************************/
/*                                                                          */
/* Function:  GalerkinCGMatrixFromInterpolation                             */
/*                                                                          */
/* Purpose:   Given a prolongation this routine computes a Galerkin         */
/*            coarse grid matrix.                                           */
/*            (changes for systems necessary!!)                             */
/*                                                                          */
/****************************************************************************/


INT GalerkinCGMatrixFromInterpolation(GRID *theGrid,
                                      MATDATA_DESC *A, MATDATA_DESC *I)
{
  INT icomp,mcomp,level;
  DOUBLE Akl,IikAkl;
  GRID *newGrid;
  VECTOR *vect,*vect2,*newVect,*newVect2;
  MATRIX *mat,*mat2,*imat,*imat2;

  level = GLEVEL(theGrid) - 1;
  if (AllocMDFromMD(MYMG(theGrid),level,level,A,&A))
    REP_ERR_RETURN(1);
  if (AssembleGalerkinByMatrix(theGrid,A,1))
    REP_ERR_RETURN(1);

  return(DONE);
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

  if (MD_IS_SCALAR(A)==FALSE)
    REP_ERR_RETURN(1);
  mcomp=MD_SCALCMP(A);

  for (vect=FIRSTVECTOR(theGrid); vect!=NULL; vect=SUCCVC(vect))
  {
    matD=VSTART(vect);
    for (mat=MNEXT(matD); mat!=NULL; mat=matS)
    {
      matS=MNEXT(mat);
      if ((STRONG(mat)==0)&&(STRONG(MADJ(mat))==0))
      {
        if (lumpFlag)
          MVALUE(matD,mcomp)+=MVALUE(mat,mcomp);

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
  VECTOR *vect,*CaV,*CeV,*FaV,*FeV,*TaV,*TeV;

    #ifndef ModelP
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
