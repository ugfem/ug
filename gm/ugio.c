// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugio.c														*/
/*																			*/
/* Purpose:   ug's grid input/output                                        */
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: bastian@iwr1.iwr.uni-heidelberg.de					*/
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

#ifdef __MPW32__
#pragma segment ugio
#endif


/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "compiler.h"
#include "fileopen.h"
#include "heaps.h"
#include "defaults.h"

#include "devices.h"

#include "switch.h"
#include "gm.h"
#include "algebra.h"
#include "misc.h"
#include "ugm.h"
#include "ugio.h"

/* include refine because of macros accessed  */
#include "ugrefine.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define BeginOfData ">-UG-DATA-<"

#define SEARCH_OFFSET 2000

/* for file transfer */
#define xOLDREFMASK 0x0000001F
#define xOLDREF(p) (((p)) & xOLDREFMASK)

#define xNEWREFMASK 0x00000F80
#define xNEWREFSHIFT 7
#define xNEWREF(p) ((((p)) & xNEWREFMASK)>>xNEWREFSHIFT)

#define xECLASSMASK 0x00000060
#define xECLASSSHIFT 5
#define xECLASS(p) (((p) & xECLASSMASK)>>xECLASSSHIFT)

#define xCOARSENMASK 0x00001000
#define xCOARSENSHIFT 12
#define xCOARSEN(p) (((p) & xCOARSENMASK)>>xCOARSENSHIFT)

#define xNSONSMASK 0x0000E000
#define xNSONSSHIFT 13
#define xNSONS(p) (((p) & xNSONSMASK)>>xNSONSSHIFT)

#define xMOVEMASK 0x00000003
#define xMOVE(p) ((*((unsigned INT *)(p))) & xMOVEMASK)

/* ug22/23 files */
#define ug22_23ECLASSMASK 0x00001800
#define ug22_23ECLASSSHIFT 11
#define ug22_23ECLASS(cw1) (((cw1) & ug22_23ECLASSMASK)>>ug22_23ECLASSSHIFT)

#define ug22_23NSONSMASK 0x0000E000
#define ug22_23NSONSSHIFT 13
#define ug22_23NSONS(cw1) (((cw1) & ug22_23NSONSMASK)>>ug22_23NSONSSHIFT)

#define ug22_23OLDREFMASK 0x000000FF
#define ug22_23OLDREF(cw2) ((cw2) & ug22_23OLDREFMASK)

#define ug22_23NEWREFMASK 0x0000FF00
#define ug22_23NEWREFSHIFT 8
#define ug22_23NEWREF(cw2) (((cw2) & ug22_23NEWREFMASK)>>ug22_23NEWREFSHIFT)

#define ug22_23COARSENMASK 0x00010000
#define ug22_23COARSENSHIFT 16
#define ug22_23COARSEN(cw2) (((cw2) & ug22_23COARSENMASK)>>ug22_23COARSENSHIFT)

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/



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

static INT gridpaths_set=FALSE;

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* Function:  MGSetVectorClasses											*/
/*																			*/
/* Purpose:   Returns highest vector class of a dof on next level			*/
/*																			*/
/* Input:	  *theElement													*/
/*																			*/
/* Output:	  INT															*/
/*																			*/
/****************************************************************************/

INT MGSetVectorClasses (MULTIGRID *theMG)
{
  INT i;
  GRID *theGrid;
  ELEMENT *theElement;

  /* set vector classes */
  for (i=0; i<=TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    if (ClearVectorClasses(theGrid)) return (1);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    {
      if (ECLASS(theElement)!=REGULAR_CLASS && ECLASS(theElement)!=IRREGULAR_CLASS) continue;
      if (SeedVectorClasses(theElement)) return (1);
    }
    if (PropagateVectorClasses(theGrid)) return (1);
  }

  /* set NextVectorClasses */
  if (ClearNextVectorClasses(GRID_ON_LEVEL(theMG,TOPLEVEL(theMG)))) return (1);
  for (i=0; i<TOPLEVEL(theMG); i++)
  {
    theGrid = GRID_ON_LEVEL(theMG,i);
    if (ClearNextVectorClasses(theGrid)) return (1);
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
    {
      if (NSONS(theElement)==0) continue;
      if (ECLASS(SON(theElement,0))!=REGULAR_CLASS && ECLASS(SON(theElement,0))!=IRREGULAR_CLASS) continue;
      if (SeedNextVectorClasses(theElement)) return (1);
    }
    if (PropagateNextVectorClasses(theGrid)) return (1);
  }

  return (0);
}

/****************************************************************************/
/*D
   SaveMultiGrid - Save complete multigrid structure in a text file

   SYNOPSIS:
   INT SaveMultiGrid (MULTIGRID *theMG, char *name, char *comment);

   PARAMETERS:
   .  theMG - pointer to multigrid
   .  name - name of the text file
   .  comment - to be included at beginning of file

   DESCRIPTION:
   This function saves a complete 'MULTIGRID' to a text file in such a way
   that it can be read again with the 'LoadMultiGrid' function.
   This function should be able to read also grid files from older versions.

   RETURN VALUE:
   INT
   .n    0 if ok
   .n    >0 if error occured.
   D*/
/****************************************************************************/

INT SaveMultiGrid (MULTIGRID *theMG, char *name, char *comment)
{
  FILE *stream;
  FORMAT *theFormat;
  GRID *theGrid;
  VERTEX *theVertex;
  NODE *theNode;
  ELEMENT *theElement;
  LINK *theLink;
  EDGE *theEdge;
  ELEMENTSIDE *theSide;
  VSEGMENT *theVertexSegment;
  int i,j,k,n,m;
  long i1,i2,i3;
  char buffer[1024];

  RenumberMultiGrid(theMG);

  if (gridpaths_set)
    /* this way grids are stored to path[0] */
    stream = FileOpenUsingSearchPaths(name,"w","gridpaths");
  else
    stream = fileopen(name,"w");
  if (stream==NULL)
  {
    PrintErrorMessage('E',"SaveMultiGrid","cannot open file");
    return(GM_FILEOPEN_ERROR);
  }

  j = theMG->topLevel;
  theFormat = theMG->theFormat;

  /* write header */
#ifdef __version23__
  fprintf(stream,">-version UG23-<\n");
#endif
#ifdef __version3__
  fprintf(stream,">-version UG31-<\n");
#endif

  fprintf(stream,"%s\n",comment);
  fprintf(stream,"%s\n",BeginOfData);

  /* save dimension */
  fprintf(stream,"(DIM = %ld)\n",DIM);

  /* start with multigrid information */
  /* NB: keep first two ints written for compatibility mode (former controlword) */
  fprintf(stream,"(MG %d %d %ld %ld %ld %ld %ld %ld \n\"%s\"\n\"%s\"\n\"%s\"\n",
          0,
          0,
          (long) theMG->status,
          (long) theMG->vertIdCounter,
          (long) theMG->nodeIdCounter,
          (long) theMG->elemIdCounter,
          (long) theMG->topLevel,
          (long) theMG->currentLevel,
          theMG->theDomain->d.name,
          theMG->theProblem->d.name,
          theMG->theFormat->d.name);
  fprintf(stream,"%ld",(long)theMG->numOfCorners);
  for (i=0; i<theMG->numOfCorners; i++)
    fprintf(stream," %ld",(long)ID(theMG->corners[i]));
  fprintf(stream,"\n");

  /* now for all levels */
  for (k=0; k<=j; k++)
  {
    theGrid = theMG->grids[k];

    /* write grid info */
    i1 = -1;
    fprintf(stream,"(GR %lx %ld %ld %ld %ld %ld %ld %ld\n",
            (unsigned long) CTRL(theGrid),
            (long) theGrid->level,
            (long) theGrid->nVert,
            (long) theGrid->nNode,
            (long) theGrid->nElem,
            (long) theGrid->nEdge,
            (long) theGrid->status,
            i1);

    /* vertices */
    for (theVertex=FIRSTVERTEX(theGrid); theVertex!=NULL; theVertex=SUCCV(theVertex))
    {
      /* header */
      if (OBJT(theVertex)==IVOBJ)
        fprintf(stream,"(IV ");
      else
        fprintf(stream,"(BV ");

      /* all vertices */
      fprintf(stream,"%lx %ld",
              (unsigned long) CTRL(theVertex),(long) ID(theVertex));
      for (i=0; i<DIM; i++)
        fprintf(stream," %lg",(double) (CVECT(theVertex)[i]));
      for (i=0; i<DIM; i++)
        fprintf(stream," %lg",(double) (LCVECT(theVertex)[i]));
      i1 = -1;
      if (VFATHER(theVertex)!=NULL) i1 = ID(VFATHER(theVertex));
      fprintf(stream," %ld",i1);

      /* data, no save with data available in ug 3.0 ! */
      fprintf(stream," ()");

      /* boundary vertices only */
      if (OBJT(theVertex)==BVOBJ)
      {
        i = 0;
        for (theVertexSegment=VSEG(theVertex); theVertexSegment!= NULL; theVertexSegment=NEXTSEG(theVertexSegment))
          i++;
        fprintf(stream," %ld",i);
        for (theVertexSegment=VSEG(theVertex); theVertexSegment!= NULL; theVertexSegment=NEXTSEG(theVertexSegment))
        {
          fprintf(stream," (VS %lx %ld",CTRL(theVertexSegment),ID(BSEGDESC(theVertexSegment)));
          for (i=0; i<DIM_OF_BND; i++)
            fprintf(stream," %lg",(double) (LAMBDA(theVertexSegment,i)));
          fprintf(stream,")");
        }
      }

      /* trailer */
      fprintf(stream,")\n");
    }

    /* nodes */
    for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
    {
      /* node info */
      i1 = i2 = i3 = -1;
      if (NFATHER(theNode)!=NULL) i1 = ID(NFATHER(theNode));
      if (SONNODE(theNode)!=NULL) i2 = ID(SONNODE(theNode));
      if (MYVERTEX(theNode)!=NULL) i3 = ID(MYVERTEX(theNode));
#ifdef __version23__
      fprintf(stream,"(ND %lx %ld %ld %hu %hu %ld %ld %ld",
              (unsigned long) CTRL(theNode),
              (long) ID(theNode),
              (long) INDEX(theNode),
              (unsigned short) VSKIP(theNode),
              (unsigned short) NSKIP(theNode),
              i1,i2,i3);
#endif
#ifdef __version3__
      fprintf(stream,"(ND %lx %ld %ld %ld %ld %ld",
              (unsigned long) CTRL(theNode),
              (long) ID(theNode),
              (long) INDEX(theNode),
              i1,i2,i3);
#endif
#ifdef __version23__
      /* data */
      fprintf(stream," ()");
      fprintf(stream," ()");
#endif
      /* trailer */
      fprintf(stream,")\n");
    }

    /* edges */
    for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
      for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
        if (ID(NBNODE(theLink))>ID(theNode))
        {
          theEdge = MYEDGE(theLink);

          /* from and to */
          fprintf(stream,"(ED %ld %ld",
                  (long) ID(NBNODE(LINK1(theEdge))),
                  (long) ID(NBNODE(LINK0(theEdge))));

          /* save ID from MidNode if */
#ifdef __THREEDIM__
          {
            i1 = -1;
            if (MIDNODE(theEdge) != NULL) i1 = ID(MIDNODE(theEdge));
            fprintf(stream," %ld",i1);
          }
#endif
#ifdef __version23__
          /* edge data */
          fprintf(stream," ()");

          /* link 0 */
          fprintf(stream," %lx",(long) CTRL(LINK0(theEdge)));
          fprintf(stream," ()");

          /* link 1 */
          fprintf(stream," %lx",(long) CTRL(LINK1(theEdge)));
          fprintf(stream," ()");
#endif

          /* trailer */
          fprintf(stream,")\n");
        }

    /* elements */
    for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
    {
      /* header */
      if (OBJT(theElement)==IEOBJ)
        fprintf(stream,"(IE ");
      else
        fprintf(stream,"(BE ");

      /* all elements */
      fprintf(stream,"%lx %lx %ld",
              (unsigned long) CTRL(theElement),(unsigned long)CTRL2(theElement),(long) ID(theElement));
      for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
        fprintf(stream," %ld",ID(CORNER(theElement,i)));
      i1 = -1;
      if (EFATHER(theElement)!=NULL) i1 = ID(EFATHER(theElement));
      fprintf(stream," %ld",i1);
#ifdef __TWODIM__
      for (i=0; i<NSONS(theElement); i++)
        fprintf(stream," %ld",ID(SON(theElement,i)));
#endif
#ifdef __THREEDIM__
      if (SON(theElement,0)!=NULL) i1 = ID(SON(theElement,0));
      fprintf(stream," %ld",i1);
#endif

      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      {
        i1 = -1;
        if (NBELEM(theElement,i)!=NULL) i1 = ID(NBELEM(theElement,i));
        fprintf(stream," %ld",i1);
      }

#ifdef __version23__
      /* element data */
      fprintf(stream," ()");
#endif

      /* boundary elements only */
      if (OBJT(theElement)==BEOBJ)
      {
        m = 0;
        for (i=0; i<SIDES_OF_ELEM(theElement); i++)
          if (SIDE(theElement,i)!=NULL)
            m++;
        fprintf(stream," %d",m);
        for (i=0; i<SIDES_OF_ELEM(theElement); i++)
        {
          theSide = SIDE(theElement,i);
          if (theSide!=NULL)
          {
            fprintf(stream," (SI %ld %lx %ld",(long)i,(unsigned long) CTRL(theSide),(long) ID(SEGDESC(theSide)));
            for (n=0; n<CORNERS_OF_SIDE(theElement,i); n++)
              for (m=0; m<DIM_OF_BND; m++)
                fprintf(stream," %lg",(double)(PARAM(theSide,n,m)));

            /* trailer */
            fprintf(stream,")");
          }
        }
      }

      /* trailer */
      fprintf(stream,")\n");
    }

    /* trailer */
    fprintf(stream,")\n");

    sprintf(buffer,"[%d] ",k);
    UserWrite(buffer);
  }

  /* trailer */
  fprintf(stream,")\n");
  UserWrite("\n");

  fclose(stream); return(GM_OK);
}

/****************************************************************************/
/*D
   LoadMultiGrid - Load complete multigrid structure from a text file

   SYNOPSIS:
   MULTIGRID *LoadMultiGrid (char *MultigridName, char *FileName,
      char *domain, char *problem, char *format, unsigned long heapSize);

   PARAMETERS:
   .  MultigridName - Name of the new 'MULTIGRID' structure in memory.
   .  FileName - Name of the file to be read.
   .  domain - `Name` of the 'DOMAIN' to be used for the 'MULTIGRID'.
   .  problem - `Name` of the 'PROBLEM' to be used for the 'MULTIGRID'.
   .  format - `Name` of the 'FORMAT' to be used for the 'MULTIGRID'.
   .  heapSize - Size of the heap in bytes that will be allocated for the 'MULTIGRID'.

   DESCRIPTION:
   This function can read grid files produced with the 'SaveMultiGrid' function.

   RETURN VALUE:
   INT
   .n    NULL if an error occured
   .n    else pointer to new 'MULTIGRID'
   D*/
/****************************************************************************/

MULTIGRID *LoadMultiGrid (char *MultigridName, char *FileName, char *domain, char *problem, char *format, unsigned long heapSize)
{
  FILE *stream;
  FORMAT *theFormat;
  GRID *theGrid;
  VERTEX *theVertex,**pv,**Vertices;
  VSEGMENT *vs;
  NODE *theNode,**Nodes;
  ELEMENT *theElement,**Elements;
  LINK *theLink;
  EDGE *theEdge;
  ELEMENTSIDE *theSide;
  HEAP *theHeap,*theUserHeap;
  MULTIGRID *theMG;
  BOUNDARY_SEGMENT *theSegment;
  BOUNDARY_CONDITION *theBndCond;
  DOMAIN *theDomain;
  PROBLEM *theProblem;
  VIRT_HEAP_MGMT *theGenMGUDM;
  INT ds;
  int i,j,j1,k,n,m,o,dim;
  long i1,i2,i3,i4,i5,l,nv,nn,ne,ned,nl,cl;
  unsigned long ul,ul2;
  char buffer[1024],c,c1,c2;
  double val,lambda;
  unsigned short us1,us2;
  long version;
  INT point,onside;
  unsigned INT cw1, cw2;
  INT maxNsubdomain,minNsubdomain,Nsubdomain,err;
  INT *counter;

  /* open file */
  if (gridpaths_set)
    stream = FileOpenUsingSearchPaths(FileName,"r","gridpaths");
  else
    stream = fileopen(FileName,"r");
  if (stream==NULL)
  {
    PrintErrorMessage('E',"LoadMultiGrid","cannot open file");
    return(NULL);
  }

  /* find version of ug file */
  if (fscanf(stream,">-version UG%ld-<",&version)!=1)
  {
    PrintErrorMessage('E',"LoadMultiGrid","I think it's not an ug file");
    fclose(stream); return(NULL);
  }
  sprintf(buffer,"grid file ug version %ld\n",version);
  UserWrite(buffer);

  /* find begin of data */
  strcpy(buffer,BeginOfData);
  n = strlen(buffer); k=0;
  while (1)
  {
    k++;
    if (k > SEARCH_OFFSET) {fclose(stream); return(NULL);}
    for (i=0; i<n; i++)
      if ((c=getc(stream))!=buffer[i])
        break;
    if (i==n) break;
  }

  /* check dimension for version >= 23 */
  if (version>=23)
  {
    fscanf(stream," (DIM =%ld",&dim);
    while ((c=getc(stream))!='\n') ;
  }
  else
  {
    dim = 2;
    PrintErrorMessage('W',"LoadMultiGrid","That is an old ug file, save it to get an up to date version");
  }
  if (dim != DIM)
  {
    PrintErrorMessage('W',"LoadMultiGrid","grid has wrong dimension");
    fclose(stream); return (NULL);
  }
#ifdef __version23__
  if (version>23)
  {
    PrintErrorMessage('W',"LoadMultiGrid","You can read only files with version <=2.3 with this application");
    fclose(stream); return(NULL);
  }
#endif

  /* allocate multigrid envitem */
  theMG = MakeMGItem(MultigridName);
  if (theMG==NULL) return(NULL);

  /* allocate the heap */
  theHeap = NewHeap(SIMPLE_HEAP, heapSize, malloc(heapSize));
  if (theHeap==NULL)
  {
    PrintErrorMessage('E',"LoadMultiGrid","could not allocate heap");
    DisposeMultiGrid(theMG);
    return(NULL);
  }


  /* read multigrid info */
  /* NB: keep first two ints read for compatibility mode (former controlword) */
  fscanf(stream," (MG %d %d %ld %ld %ld %ld %ld %ld",&i,&i,&i1,&nv,&nn,&ne,&nl,&cl);
  theMG->status = (INT) i1;
  theMG->currentLevel = 0;
  theMG->vertIdCounter = 0;
  theMG->nodeIdCounter = 0;
  theMG->elemIdCounter = 0;
  theMG->topLevel = -1;
  SELECTIONSIZE(theMG) = 0;

  /* find domain structure */
  while ((c=getc(stream))!='"') ;
  i = 0;
  while ((c=getc(stream))!='"') buffer[i++] = c;
  buffer[i] = (char) 0;
  if (domain==NULL)
    theDomain = GetDomain(buffer);
  else
    theDomain = GetDomain(domain);
  if (theDomain==NULL)
  {
    PrintErrorMessage('E',"LoadMultiGrid","domain not found");
    DisposeMultiGrid(theMG);
    fclose(stream); return(NULL);
  }

  /* find problem structure */
  while ((c=getc(stream))!='"') ;
  i = 0;
  while ((c=getc(stream))!='"') buffer[i++] = c;
  buffer[i] = (char) 0;
  if (problem==NULL)
    theProblem = GetProblem(ENVITEM_NAME(theDomain),buffer);
  else
    theProblem = GetProblem(ENVITEM_NAME(theDomain),problem);
  if (theProblem==NULL)
  {
    PrintErrorMessage('E',"LoadMultiGrid","problem not found");
    DisposeMultiGrid(theMG);
    fclose(stream); return(NULL);
  }

  /* find format structure */
  while ((c=getc(stream))!='"') ;
  i = 0;
  while ((c=getc(stream))!='"') buffer[i++] = c;
  buffer[i] = (char) 0;
  if (format==NULL)
    theFormat = GetFormat(buffer);
  else
    theFormat = GetFormat(format);
  if (theFormat==NULL)
  {
    PrintErrorMessage('E',"LoadMultiGrid","format not found");
    DisposeMultiGrid(theMG);
    fclose(stream); return(NULL);
  }


  /* allocate user data from that heap */

  /* 1: general user data space */
  theGenMGUDM = GetGenMGUDM();
  if (!theGenMGUDM->locked)
    CalcAndFixTotalSize(theGenMGUDM);
  ds = theGenMGUDM->TotalSize;
  if (ds!=0)
  {
    GEN_MGUD(theMG) = GetMem(theHeap,ds,FROM_BOTTOM);
    if (GEN_MGUD(theMG)==NULL) {DisposeMultiGrid(theMG); fclose(stream); return(NULL);}
    memset(GEN_MGUD(theMG),0,ds);
  }
  else
    GEN_MGUD(theMG) = NULL;

  /* 2: user heap */
  ds = theFormat->sMultiGrid;
  if (ds!=0)
  {
    theUserHeap = NewHeap(SIMPLE_HEAP, ds, GetMem(theHeap,ds,FROM_BOTTOM));
    if (theUserHeap==NULL)
    {
      DisposeMultiGrid(theMG);
      return(NULL);
    }
    MG_USER_HEAP(theMG) = theUserHeap;
  }
  else
    MG_USER_HEAP(theMG) = NULL;

  /* domain, problem and format */
  theMG->theDomain = theDomain;
  theMG->theFormat = theFormat;
  theMG->theProblem = theProblem;
  theMG->theHeap = theHeap;
  for (i=0; i<MAXLEVEL; i++) theMG->grids[i] = NULL;
  for (i=0; i<MAXOBJECTS; i++) theMG->freeObjects[i] = NULL;
  for (i=0; i<MAXVECTORS; i++) theMG->freeVectors[i] = NULL;
  for (i=0; i<MAXCONNECTIONS; i++) theMG->freeConnections[i] = NULL;

  /* allocate boundary descriptors */
  n = theDomain->numOfSegments;
  theMG->numOfSegments = n;
  theMG->segments = (BNDSEGDESC *) GetMem(theHeap,n*sizeof(BNDSEGDESC),FROM_BOTTOM);
  if (theMG->segments==NULL)
  {
    PrintErrorMessage('E',"LoadMultiGrid","no memory for segments");
    DisposeMultiGrid(theMG);
    fclose(stream); return(NULL);
  }
  for (i=0; i<n; i++)
  {
    theMG->segments[i].theSegment = NULL;
    theMG->segments[i].theBoundaryCondition = NULL;
  }

  /* combine boundary coordinates and boundary conditions */
  for (theSegment=GetFirstBoundarySegment(theDomain); theSegment!=NULL; theSegment = GetNextBoundarySegment(theSegment))
  {
    i = theSegment->id;
    if ((i<0)||(i>=n))
    {
      PrintErrorMessage('E',"LoadMultiGrid","segment id out of range");
      DisposeMultiGrid(theMG);
      fclose(stream); return(NULL);
    }
    if (theMG->segments[i].theSegment!=NULL)
    {
      PrintErrorMessage('E',"LoadMultiGrid","multiply defined segment");
      DisposeMultiGrid(theMG);
      fclose(stream); return(NULL);
    }
    theMG->segments[i].theSegment = theSegment;
    theMG->segments[i].id = i;
  }
  for (theBndCond=GetFirstBoundaryCondition(theProblem); theBndCond!=NULL; theBndCond = GetNextBoundaryCondition(theBndCond))
  {
    i = theBndCond->id;
    if ((i<0)||(i>=n))
    {
      PrintErrorMessage('E',"LoadMultiGrid","boundary condition id out of range");
      DisposeMultiGrid(theMG);
      fclose(stream); return(NULL);
    }
    if (theMG->segments[i].theBoundaryCondition!=NULL)
    {
      PrintErrorMessage('E',"LoadMultiGrid","multiply defined boundary condition");
      DisposeMultiGrid(theMG);
      fclose(stream); return(NULL);
    }
    theMG->segments[i].theBoundaryCondition = theBndCond;
  }

  /* check if all pointers are correctly defined */
  maxNsubdomain = -MAX_I;
  minNsubdomain =  MAX_I;
  for (i=0; i<n; i++)
  {
    theSegment = theMG->segments[i].theSegment;

    if (theSegment==NULL)
    {
      PrintErrorMessage('E',"LoadMultiGrid","boundary segment not found");
      DisposeMultiGrid(theMG);
      return(NULL);
    }
    if (theMG->segments[i].theBoundaryCondition==NULL)
    {
      PrintErrorMessage('E',"LoadMultiGrid","boundary conditionnot found");
      DisposeMultiGrid(theMG);
      return(NULL);
    }
    if (theSegment->left==theSegment->right)
    {
      sprintf(buffer,"ERROR: left==right for segment %d",i);
      PrintErrorMessage('E',"LoadMultiGrid",buffer);
      DisposeMultiGrid(theMG);
      return(NULL);
    }
    maxNsubdomain = MAX(maxNsubdomain,MAX(theSegment->left,theSegment->right));
    minNsubdomain = MIN(minNsubdomain,MIN(theSegment->left,theSegment->right));
  }

  /* check subdomains */

  if (minNsubdomain!=0)
  {
    PrintErrorMessage('E',"LoadMultiGrid","ERROR: subdomain IDs must be >= 0");
    DisposeMultiGrid(theMG);
    return(NULL);
  }

  /* get storage for subdomain counters */
  Nsubdomain = maxNsubdomain+1;
  Mark(theHeap,FROM_TOP);
  counter = (INT *) GetMem(theHeap,Nsubdomain*sizeof(INT),FROM_TOP);
  for (i=0; i<Nsubdomain; i++) counter[i] = 0;
  for (i=0; i<n; i++)
  {
    counter[LEFT(theMG->segments+i)]++;
    counter[RIGHT(theMG->segments+i)]++;
  }
  err = 0;
  for (i=0; i<Nsubdomain; i++)
    if (counter[i]<=0)
    {
      err++;
      sprintf(buffer,"ERROR: subdomain ID %d not used\n",i);
      UserWrite(buffer);
    }
  Release(theHeap,FROM_TOP);
  if (err)
  {
    free(theHeap);
    return (NULL);
  }
  theMG->numOfSubdomains = Nsubdomain;

  /* load corner vertices ids */
  n = theDomain->numOfCorners;
  theMG->numOfCorners = n;
  pv = (VERTEX **) GetMem(theHeap,n*sizeof(VERTEX *),FROM_BOTTOM);
  if (pv==NULL)
  {
    PrintErrorMessage('E',"LoadMultiGrid","not enough memory for corners array");
    DisposeMultiGrid(theMG);
    fclose(stream); return(NULL);
  }
  theMG->corners = pv;
  fscanf(stream," %ld",&l);
  if (l!=n)
  {
    PrintErrorMessage('E',"LoadMultiGrid","number of corners in domain and file different");
    DisposeMultiGrid(theMG);
    fclose(stream); return(NULL);
  }
  for (i=0; i<n; i++)
  {
    /* store ids from verticees */
    fscanf(stream," %ld",&l);
    pv[i] = (VERTEX *) l;
  }

  /* allocate the big arrays of pointers */
  Mark(theHeap,FROM_TOP);
  Vertices = (VERTEX **) GetMem(theHeap,nv*sizeof(VERTEX *),FROM_TOP);
  if (Vertices==NULL)
  {
    PrintErrorMessage('E',"LoadMultiGrid","not enough memory for pointer array");
    DisposeMultiGrid(theMG);
    fclose(stream); return(NULL);
  }
  Nodes = (NODE **) GetMem(theHeap,nn*sizeof(NODE *),FROM_TOP);
  if (Nodes==NULL)
  {
    PrintErrorMessage('E',"LoadMultiGrid","not enough memory for pointer array");
    DisposeMultiGrid(theMG);
    fclose(stream); return(NULL);
  }
  Elements = (ELEMENT **) GetMem(theHeap,ne*sizeof(ELEMENT *),FROM_TOP);
  if (Elements==NULL)
  {
    PrintErrorMessage('E',"LoadMultiGrid","not enough memory for pointer array");
    DisposeMultiGrid(theMG);
    fclose(stream); return(NULL);
  }

  /* read level per level */
  for (k=0; k<=nl; k++)
  {
    theGrid = CreateNewLevel(theMG);
    if (theGrid==NULL)
    {
      PrintErrorMessage('E',"LoadMultiGrid","cannot allocate level");
      DisposeMultiGrid(theMG);
      fclose(stream); return(NULL);
    }
    fscanf(stream," (GR %lx %ld %ld %ld %ld %ld %ld %ld",&ul,&i1,&nv,&nn,&ne,&ned,&i2,&i3);
    CTRL(theGrid) = (unsigned INT) ul;
    theGrid->status = (INT) i2;
    /* i3 contains the old activeNodes pointer no longer supported since it has never been used */

    /* vertices */
    theVertex = NULL;
    for (i=0; i<nv; i++)
    {
      fscanf(stream," (%c%c %lx %ld",&c1,&c2,&ul,&i1);
      if (c1=='B')
        theVertex = CreateBoundaryVertex(theGrid,theVertex);
      else
        theVertex = CreateInnerVertex(theGrid,theVertex);
      if (theVertex==NULL)
      {
        PrintErrorMessage('E',"LoadMultiGrid","cannot allocate vertex");
        DisposeMultiGrid(theMG);
        fclose(stream); return(NULL);
      }
      CTRL(theVertex) = (unsigned INT) ul;
      cw1 = CTRL(theVertex);
      Vertices[i1] = theVertex;

      for (j=0; j<DIM; j++)
      {
        fscanf(stream," %lf",&val);
        CVECT(theVertex)[j] = (COORD) val;
      }
      for (j=0; j<DIM; j++)
      {
        fscanf(stream," %lf",&val);
        LCVECT(theVertex)[j] = (COORD) val;
      }

      fscanf(stream," %ld",&l);
      VFATHER(theVertex) = (ELEMENT *) l;

      /* data */
      while ((c=getc(stream))!='(') ;
      j = 0;
      while ((c=getc(stream))!=')') buffer[j++] = c;
      buffer[j] = (char) 0;

      /* adjust MOVE Flag for inner vertex */
      if ((version<23)&&(OBJT(theVertex)==IVOBJ))
      {
        SETMOVED(theVertex,0);
        SETMOVE(theVertex,2);
        SETONEDGE(theVertex,0);
      }

      /* boundary vertices only */
      if (OBJT(theVertex)==BVOBJ)
      {
        if (version>=23)                         /* read vertex segments */
        {
          fscanf(stream," %ld",&l);
          for (j=0; j<l; j++)
          {
            fscanf(stream," (VS %lx %ld",&ul, &i1);
            vs = CreateVertexSegment (theGrid,theVertex);
            if (vs == NULL)
            {
              PrintErrorMessage('E',"LoadMultiGrid","cannot allocate vertex segment");
              DisposeMultiGrid(theMG);
              fclose(stream); return(NULL);
            }
            CTRL(vs) = (unsigned INT) ul;
            BSEGDESC(vs) = &(theMG->segments[i1]);
            for (j1=0; j1<DIM_OF_BND; j1++)
            {
              fscanf(stream," %lf",&val);
              LAMBDA(vs,j1) = (COORD) val;
            }
            fscanf(stream,")");
          }
        }
        else                         /* we must construct vertex segments somehow ... */
        {
          /* if vertex can be moved we have only one segment */
          if (xMOVE(theVertex)!=0)
          {
            vs = CreateVertexSegment (theGrid,theVertex);
            if (vs == NULL)
            {
              PrintErrorMessage('E',"LoadMultiGrid","cannot allocate vertex segment");
              DisposeMultiGrid(theMG);
              fclose(stream); return(NULL);
            }
            fscanf(stream," %ld",&i1);                                     /* segment no. */
            BSEGDESC(vs) = &(theMG->segments[i1]);
            for (j1=0; j1<DIM_OF_BND; j1++)
            {
              fscanf(stream," %lf",&val);
              LAMBDA(vs,j1) = (COORD) val;
            }
            for (j1=0; j1<DIM_OF_BND; j1++)
            {
              fscanf(stream," %lf",&val);                                           /* dummy read */
            }
            fscanf(stream," %ld",&i1);
            onside = i1;
            SETMOVE(theVertex,1);
          }
          else                               /* scan segment list */
          {
            /* find logical number of vertex */
            fscanf(stream," %ld",&i1);                                     /* segment no. */
            theSegment = theMG->segments[i1].theSegment;
            /* old version can only be 2D ! */
            fscanf(stream," %lf",&lambda);
            fscanf(stream," %lf",&val);                                     /* dummy read */
            fscanf(stream," %ld",&i1);
            onside = i1;
            point = -1;
            if (fabs(lambda-theSegment->alpha[0])<1.0E-10)
              point = theSegment->points[0];
            if (fabs(lambda-theSegment->beta[0])<1.0E-10)
              point = theSegment->points[1];
            if (point<0)
            {
              PrintErrorMessage('E',"LoadMultiGrid","cannot fill vseg list");
              DisposeMultiGrid(theMG);
              fclose(stream); return(NULL);
            }
            for (j1=0; j1<theMG->numOfSegments; j1++)
            {
              theSegment = theMG->segments[j1].theSegment;
              if (theSegment->points[0]==point)
              {
                vs = CreateVertexSegment (theGrid,theVertex);
                if (vs == NULL)
                {
                  PrintErrorMessage('E',"LoadMultiGrid","cannot allocate vertex segment");
                  DisposeMultiGrid(theMG);
                  fclose(stream); return(NULL);
                }
                BSEGDESC(vs) = &(theMG->segments[j1]);
                LAMBDA(vs,0) = theSegment->alpha[0];
              }
              if (theSegment->points[1]==point)
              {
                vs = CreateVertexSegment (theGrid,theVertex);
                if (vs == NULL)
                {
                  PrintErrorMessage('E',"LoadMultiGrid","cannot allocate vertex segment");
                  DisposeMultiGrid(theMG);
                  fclose(stream); return(NULL);
                }
                BSEGDESC(vs) = &(theMG->segments[j1]);
                LAMBDA(vs,0) = theSegment->beta[0];
              }
            }
            SETMOVE(theVertex,0);
          }
          SETONEDGE(theVertex,onside);
          SETMOVED(theVertex,0);
        }
      }

      /* trailer */
      fscanf(stream," )");
    }

    /* nodes */
    theNode = NULL;
    for (i=0; i<nn; i++)
    {
      theNode = CreateNode(theGrid,theNode);
      if (theNode==NULL)
      {
        PrintErrorMessage('E',"LoadMultiGrid","cannot allocate node");
        DisposeMultiGrid(theMG);
        fclose(stream); return(NULL);
      }
      if (version<30)
        fscanf(stream," (%c%c %lx %ld %ld %hu %hu %ld %ld %ld",&c1,&c2,&ul,&i1,&i2,&us1,&us2,&i3,&i4,&i5);
      else
        /* no skip fields */
        fscanf(stream," (%c%c %lx %ld %ld %ld %ld %ld",&c1,&c2,&ul,&i1,&i2,&i3,&i4,&i5);
      CTRL(theNode) = (unsigned INT) ul;
      Nodes[i1] = theNode;
      INDEX(theNode) = (INT) i2;
#ifdef __version23__
      VSKIP(theNode) = (unsigned SHORT) us1;
      NSKIP(theNode) = (unsigned SHORT) us2;
#endif
      NFATHER(theNode) = (NODE *) i3;
      SONNODE(theNode) = (NODE *) i4;
      MYVERTEX(theNode) = (VERTEX *) i5;

      if (version<=23)                   /* then we have data fields (possibly empty) */
      {
        /* read data in any case */
        while ((c=getc(stream))!='(') ;
        j = 0;
        while ((c=getc(stream))!=')') buffer[j++] = c;
        buffer[j] = (char) 0;
        /* no dataopt in version 3 */
        while ((c=getc(stream))!='(') ;
        j = 0;
        while ((c=getc(stream))!=')') buffer[j++] = c;
        buffer[j] = (char) 0;
        /* no dataopt in version 3 */
      }

      /* trailer */
      fscanf(stream," )");
    }

    /* edges */
    for (i=0; i<ned; i++)
    {
      /* from and to */
      fscanf(stream," (%c%c %ld %ld",&c1,&c2,&i1,&i2);

      theEdge = CreateEdge(theGrid,Nodes[i1],Nodes[i2]);
      if (theEdge==NULL)
      {
        PrintErrorMessage('E',"LoadMultiGrid","cannot allocate edge");
        DisposeMultiGrid(theMG);
        fclose(stream); return(NULL);
      }

      /* load ID from MidNode if */
                        #ifdef __THREEDIM__
      {
        fscanf(stream," %ld",&i3);
        MIDNODE(theEdge) = (NODE *)i3;
      }
                        #endif

      /* edge data */
      if (version<=23)                   /* then we have data fields (possibly empty) */
      {
        while ((c=getc(stream))!='(') ;
        j = 0;
        while ((c=getc(stream))!=')') buffer[j++] = c;
        buffer[j] = (char) 0;
        /* no dataopt in version 3 */

        /* link 0 */
        fscanf(stream," %lx",&ul);
        CTRL(LINK0(theEdge)) = (unsigned INT) ul;
        while ((c=getc(stream))!='(') ;
        j = 0;
        while ((c=getc(stream))!=')') buffer[j++] = c;
        buffer[j] = (char) 0;
        /* no dataopt in version 3 */

        /* link 1 */
        fscanf(stream," %lx",&ul);
        CTRL(LINK1(theEdge)) = (unsigned INT) ul;
        while ((c=getc(stream))!='(') ;
        j = 0;
        while ((c=getc(stream))!=')') buffer[j++] = c;
        buffer[j] = (char) 0;
        /* no dataopt in version 3 */
      }

      /* trailer */
      fscanf(stream," )");
    }

    /* elements */
    theElement = NULL;
    for (i=0; i<ne; i++)
    {
      if (version<21)
      {
        fscanf(stream," (%c%c %lx %ld",&c1,&c2,&ul,&i1);
        ul2=0;
      }
      if (version==22)
      {
        fscanf(stream," (%c%c %lx %ld %lx",&c1,&c2,&ul,&i1,&ul2);
      }
      if (version>=23)
      {
        fscanf(stream," (%c%c %lx %lx %ld",&c1,&c2,&ul,&ul2,&i1);
      }
      if (c1=='B')
        theElement = CreateBoundaryElement(theGrid,theElement,TAG(&ul));
      else
        theElement = CreateInnerElement(theGrid,theElement,TAG(&ul));
      if (theElement==NULL)
      {
        PrintErrorMessage('E',"LoadMultiGrid","cannot allocate element");
        DisposeMultiGrid(theMG);
        fclose(stream); return(NULL);
      }

      /* all elements */
      CTRL(theElement) = (unsigned INT) ul;
      CTRL2(theElement) = (unsigned INT) ul2;
      Elements[i1] = theElement;
      for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
      {
        fscanf(stream," %ld",&l);
        SET_CORNER(theElement,j,(NODE *) l);
      }
      cw1 = CTRL(theElement); cw2 = CTRL2(theElement);
      if (version<21)
      {
        SETMARK(theElement,xNEWREF(cw1));
        SETREFINE(theElement,xOLDREF(cw1));
        SETCOARSEN(theElement,xCOARSEN(cw1));
        SETECLASS(theElement,xECLASS(cw1));
        SETNSONS(theElement,xNSONS(cw1));
        SETDECOUPLED(theElement,0);
      }
      if ((version>=22)&&(version<=23))
      {
        SETMARK(theElement,ug22_23NEWREF(cw2));
        SETREFINE(theElement,ug22_23OLDREF(cw2));
        SETCOARSEN(theElement,ug22_23COARSEN(cw2));
        SETECLASS(theElement,ug22_23ECLASS(cw1));
        SETNSONS(theElement,ug22_23NSONS(cw1));
      }
      fscanf(stream," %ld",&l);
      SET_EFATHER(theElement,(ELEMENT *) l);
                        #ifdef __TWODIM__
      for (j=0; j<NSONS(theElement); j++)
      {
        fscanf(stream," %ld",&l);
        SET_SON(theElement,j,(ELEMENT *) l);
      }
                        #endif
                        #ifdef __THREEDIM__
      fscanf(stream," %ld",&l);
      SET_SON(theElement,0,(ELEMENT *) l);
                        #endif

      for (j=0; j<SIDES_OF_ELEM(theElement); j++)
      {
        fscanf(stream," %ld",&l);
        SET_NBELEM(theElement,j,(ELEMENT *) l);
      }

      /* element data */
      if (version<=23)
      {
        while ((c=getc(stream))!='(') ;
        j = 0;
        while ((c=getc(stream))!=')') buffer[j++] = c;
        buffer[j] = (char) 0;
        /* no dataopt in version 3 */
      }

      /* boundary elements only */
      if (OBJT(theElement)==BEOBJ)
      {
        fscanf(stream," %d",&o);
        for (j=0; j<o; j++)
        {
          theSide = CreateElementSide(theGrid);
          if (theSide==NULL)
          {
            PrintErrorMessage('E',"LoadMultiGrid","cannot allocate element side");
            DisposeMultiGrid(theMG);
            fclose(stream); return(NULL);
          }
          fscanf(stream," (%c%c %ld %lx %ld",&c1,&c2,&i1,&ul,&i2);
          SET_SIDE(theElement,i1,theSide);
          CTRL(theSide) = (unsigned INT) ul;
          SEGDESC(theSide) = &(theMG->segments[i2]);
          for (n=0; n<CORNERS_OF_SIDE(theElement,i1); n++)
            for (m=0; m<DIM_OF_BND; m++)
            {
              fscanf(stream," %lf",&val);
              PARAM(theSide,n,m) = (COORD) val;
            }

          /* trailer */
          fscanf(stream," )");
        }
      }

      /* trailer */
      fscanf(stream," )");
    }

    /* trailer */
    fscanf(stream," )");

    sprintf(buffer,"[%d] ",k);
    UserWrite(buffer);
  }

  /* trailer */
  fscanf(stream," )");
  UserWrite("\n");

  /* pointers in multigrid */
  theMG->currentLevel = cl;
  for (i=0; i<theMG->numOfCorners; i++)
    pv[i] = Vertices[((INT)pv[i])];

  /* scan through levels */
  for (k=0; k<=nl; k++)
  {
    theGrid = theMG->grids[k];

    /* vertices */
    for (theVertex=FIRSTVERTEX(theGrid); theVertex!=NULL; theVertex=SUCCV(theVertex))
    {
      l = ((long)VFATHER(theVertex));
      if (l>=0)
        VFATHER(theVertex) = Elements[l];
      else
        VFATHER(theVertex) = NULL;
    }

    /* nodes */
    for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
    {
      l = ((long)NFATHER(theNode));
      if (l>=0)
        NFATHER(theNode) = Nodes[l];
      else
        NFATHER(theNode) = NULL;
      l = ((long)SONNODE(theNode));
      if (l>=0)
        SONNODE(theNode) = Nodes[l];
      else
        SONNODE(theNode) = NULL;
      l = ((long)MYVERTEX(theNode));
      if (l>=0)
        MYVERTEX(theNode) = Vertices[l];
      else
        MYVERTEX(theNode) = NULL;
      if (CLASS(theNode)>=2)
        TOPNODE(MYVERTEX(theNode)) = theNode;
    }

    if (DIM == 3)
      for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
        for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
          if (ID(NBNODE(theLink))>ID(theNode))
          {
            theEdge = MYEDGE(theLink);
            /*						l = (long)MIDNODE(theEdge);
                                                            if (l>=0)
                                                                    MIDNODE(theEdge) = Nodes[l];
                                                            else
                                                                    MIDNODE(theEdge) = NULL;*/
          }

    /* elements */
    for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
    {
      for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
      {
        l = ((long)CORNER(theElement,i));
        SET_CORNER(theElement,i,Nodes[l]);
      }
      l = ((long)EFATHER(theElement));
      if (l>=0)
        SET_EFATHER(theElement,Elements[l]);
      else
        SET_EFATHER(theElement,NULL);
      for (i=0; i<NSONS(theElement); i++)
      {
        l = ((long)SON(theElement,i));
        SET_SON(theElement,i,Elements[l]);
      }
      for (i=0; i<SIDES_OF_ELEM(theElement); i++)
      {
        l = ((long)NBELEM(theElement,i));
        if (l>=0)
          SET_NBELEM(theElement,i,Elements[l]);
        else
          SET_NBELEM(theElement,i,NULL);
      }
    }
  }

  /* throw away pointers */
  Release(theHeap,FROM_TOP);

#ifdef __version3__
  /* handle algebra */
  if (MGCreateConnection(theMG))
  {
    UserWrite("could not create connection in multigrid\n");
    DisposeMultiGrid(theMG);
    return(NULL);
  }
  if (MGSetVectorClasses(theMG))
  {
    UserWrite("cannot comput vector classes in multigrid\n");
    DisposeMultiGrid(theMG);
    return(NULL);
  }
#endif

  /* return ok */
  fclose(stream); return(theMG);
}

INT InitUgio ()
{
  /* read gridpaths from defaults file (iff) */
  gridpaths_set = FALSE;
  if (ReadSearchingPaths(DEFAULTSFILENAME,"gridpaths")==0)
    gridpaths_set = TRUE;

  return (0);
}
