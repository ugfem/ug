// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  dataexplorer.c												*/
/*																			*/
/* Purpose:   dataexplorer output			                                                                */
/*																			*/
/* Author:	  Volker Reichenberger / Michael Lampe							*/
/*			  IWR Technische Simulation				                                                */
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  69120 Heidelberg												*/
/*			  email: volker.reichenberger@iwr.uni-heidelberg.de				*/
/*                   michael.lampe@iwr.uni-heidelberg.de                    */
/*																			*/
/* History:   23.03.2000 begin												*/
/*																			*/
/* Remarks:   DX supports only irregular grids made of triangles/           */
/*            quadrilaterals (2D) or tetrahedrons/hexahedrons (3D),         */
/*            respectively -- and you can't mix the two types. So we        */
/*            use degenerate quadrilaterals/hexahedrons to emulate the      */
/*            missing element types. Let's hope this works in all cases!    */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			system include files											*/
/*			application include files										*/
/*																			*/
/****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "ugdevices.h"
#include "enrol.h"
#include "compiler.h"
#include "misc.h"
#include "general.h"
#include "pfile.h"

#include "gm.h"
#include "ugenv.h"
#include "ugm.h"
#include "algebra.h"
#include "cmdint.h"
#include "commands.h"
#include "helpmsg.h"
#include "shapes.h"
#include "cmdline.h"
#include "num.h"
#include "rm.h"

#include "dataexplorer.h"

/****************************************************************************/
/*																			*/
/* defines in the following order										        */
/*																			*/
/*		compile time constants defining static data size (i.e. arrays)		*/
/*		other constants												                */
/*		macros																*/
/*																			*/
/****************************************************************************/

#define MAXVARIABLES    50                      /* max number of eval procs				*/

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* data for CVS	*/
static char RCS_ID("$Header$",UG_RCS_STRING);

/****************************************************************************/
/*																			*/
/* definition of functions													*/
/*																			*/
/****************************************************************************/

#ifdef ModelP
static INT get_offset (INT n)
{
  INT i,subtreesize[50],sum,offset;

  /* get number of objects downtree */
  sum = n;
  for (i=0; i<degree; i++) {
    GetConcentrate(i,subtreesize+i,sizeof(INT));
    sum += subtreesize[i];
  }

  /* get offset */
  if (me==master)
  {
    offset = 0;
  }
  else
  {
    Concentrate(&sum,sizeof(INT));
    GetSpread(&offset,sizeof(INT));
  }

  /* send offsets for downtree nodes */
  sum = offset+n;
  for (i=0; i<degree; i++) {
    Spread(i,&sum,sizeof(INT));
    sum += subtreesize[i];
  }

  return(offset);
}

/*
   static void LocallyUniqueIDs (MULTIGRID *mg)
   {
        VERTEX *vx;
        INT k, nv;

        nv = 0;
        for (k=0; k<=TOPLEVEL(mg); k++)
                for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
                        ID(vx) = nv++;
   }
 */

static void LocallyUniqueIDs (MULTIGRID *mg)
{
  VERTEX *vx;
  ELEMENT *el;
  INT i, k, nv;

  for (k=0; k<=TOPLEVEL(mg); k++)
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
      SETUSED(vx,0);

  nv=0;
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el)) {
      if (!EstimateHere(el)) continue;
      for (i=0; i<CORNERS_OF_ELEM(el); i++) {
        vx = MYVERTEX(CORNER(el,i));
        if (USED(vx)) continue;
        SETUSED(vx,1);
        ID(vx) = nv++;
      }
    }
}

#endif

static double clampf(double x)
{
  if (x < (double)FLT_MIN && x > (double)(-FLT_MIN))
    return 0.0;
  else
    return x;
}

/****************************************************************************/
/*D
   dataexplorer - file output in DataExplorer format

   DESCRIPTION:
   The dataexplorer command writes the grid and multiple vector or scalar
   grid functions in DataExplorer readable format to a file.

   'dataexplorer <filename> [$ns <nep> $s <vd>]* [$nv <nep> $s <vd>]*
                                            [$cs <eep> $s <vd>]* [$cv <eep> $s <vd>]*'

   .  $ns...			- plot function for scalar nodal values
   .  $nv...			- plot function for vector nodal values
   .  $cs...			- plot function for scalar element values
   .  $cv...			- plot function for vector element values

   .  <vd>				- vecdata desc
   .  <nep>			- eval proc (nodal values)
   .  <eep>			- eval proc (element values)

   KEYWORDS:
   graphics, plot, file, output, DataExplorer

   EXAMPLE:
   'dataexplorer NavierStokesSolution.dx $ns nvalue $s psol $nv nvector $s velsol'
   D*/
/****************************************************************************/

#ifdef __MWCW__
#pragma global_optimizer on
#pragma optimization_level 1
#endif

static INT DataExplorerCommand (INT argc, char **argv)
{
  INT i,j,k,v;                                  /* counters etc.							*/
  INT counter;                                  /* for formatting output					*/
  char item[1024],it[256];              /* item buffers								*/
  INT ic=0;                                             /* item length								*/
  VERTEX *vx;                                           /* a vertex pointer							*/
  ELEMENT *el;                                  /* an element pointer						*/

  MULTIGRID *mg;                                /* our multigrid							*/
  HEAP *heap;

  char filename[NAMESIZE];              /* file name for output file				*/
  PFILE *pf;                                            /* the output file pointer					*/

  INT ns;                                               /* number of scalar eval procs				*/
  INT nv;                                               /* number of vector eval procs				*/
  EVALUES *es[MAXVARIABLES];            /* pointers to scalar eval function desc	*/
  char es_name[MAXVARIABLES][NAMESIZE];         /* names for eval functions     */
  EVECTOR *ev[MAXVARIABLES];            /* pointers to vector eval function desc	*/
  char ev_name[MAXVARIABLES][NAMESIZE];         /* names for eval functions     */
  INT ns_cell;                                  /* number of scalar eval procs				*/
  INT nv_cell;                                  /* number of vector eval procs				*/
  EVALUES *es_cell[MAXVARIABLES];       /* pointers to scalar eval function desc*/
  char es_cell_name[MAXVARIABLES][NAMESIZE];            /* names for eval functions	*/
  EVECTOR *ev_cell[MAXVARIABLES];       /* pointers to vector eval function desc*/
  char ev_cell_name[MAXVARIABLES][NAMESIZE];            /* names for eval functions	*/
  char s[NAMESIZE];                             /* name of eval proc						*/

  INT numVertices;                              /* number of data points locally			*/
  INT numElements;                              /* number of elements locally				*/
  INT gnumVertices;                             /* number of data points globally			*/
  INT gnumElements;                             /* number of elements globallay				*/

  PreprocessingProcPtr pre;             /* pointer to prepare function				*/
  ElementEvalProcPtr eval_s;            /* pointer to scalar evaluation function	*/
  ElementVectorProcPtr eval_v;      /* pointer to vector evaluation function	*/

  DOUBLE *CornersCoord[MAX_CORNERS_OF_ELEM];       /* pointers to coordinates   */
  DOUBLE LocalCoord[DIM];               /* is one of the corners local coordinates	*/
  DOUBLE local[DIM];                            /* local coordinate in DOUBLE				*/
  DOUBLE value;                                 /* returned by user eval proc				*/
  DOUBLE x,y,z;                                 /* scalar values							*/
  DOUBLE vval[DIM];                             /* result of vector evaluation function		*/

  time_t ltime;

  INT oe,ov;
  INT blocks;

  INT *Id2Position, key;

  /* get current multigrid	*/
  mg = GetCurrentMultigrid();
  if (mg==NULL)
  {
    PrintErrorMessage('W',"dataexplorer","no multigrid open\n");
    return (OKCODE);
  }

  /* scan options	*/
  ns = nv = ns_cell = nv_cell = 0;

  for(i=1; i<argc; i++)
  {
    if (strncmp(argv[i],"ns",2)==0) {
      if (ns>=MAXVARIABLES)
      {
        PrintErrorMessage('E',"dataexplorer:","too many scalar variables "\
                          "specified\n");
        break;
      }
      sscanf(argv[i],"ns %s", s);
      es[ns] = GetElementValueEvalProc(s);
      if (es[ns]==NULL)
      {
        PrintErrorMessageF('E',"dataexplorer:","could not find scalar "\
                           "eval proc %s\n",s);
        break;
      }
      if (sscanf(argv[i+1],"s %s", s) == 1)
      {
        strcpy(es_name[ns],s);
        i++;
      }
      else
        strcpy(es_name[ns],es[ns]->v.name);
      ns++;
      continue;
    }

    if (strncmp(argv[i],"nv",2)==0) {
      if (nv>=MAXVARIABLES)
      {
        PrintErrorMessage('E',"dataexplorer:","too many vector variables "\
                          "specified\n");
        break;
      }
      sscanf(argv[i],"nv %s", s);
      ev[nv] = GetElementVectorEvalProc(s);
      if (ev[nv]==NULL)
      {
        PrintErrorMessageF('E',"dataexplorer:","could not find vector "\
                           "eval proc %s\n",s);
        break;
      }
      if (sscanf(argv[i+1],"s %s", s) == 1)
      {
        strcpy(ev_name[nv],s);
        i++;
      }
      else
        strcpy(ev_name[nv],ev[nv]->v.name);
      nv++;
      continue;
    }

    if (strncmp(argv[i],"cs",2)==0) {
      if (ns_cell>=MAXVARIABLES)
      {
        PrintErrorMessage('E',"dataexplorer:","too many scalar variables "\
                          "specified\n");
        break;
      }
      sscanf(argv[i],"cs %s", s);
      es_cell[ns_cell] = GetElementValueEvalProc(s);
      if (es_cell[ns_cell]==NULL)
      {
        PrintErrorMessageF('E',"dataexplorer:","could not find scalar "\
                           "eval proc %s\n",s);
        break;
      }
      if (sscanf(argv[i+1],"s %s", s) == 1)
      {
        strcpy(es_cell_name[ns_cell],s);
        i++;
      }
      else
        strcpy(es_cell_name[ns_cell],es_cell[ns_cell]->v.name);
      ns_cell++;
      continue;
    }

    if (strncmp(argv[i],"cv",2)==0) {
      if (nv_cell>=MAXVARIABLES)
      {
        PrintErrorMessage('E',"dataexplorer:","too many vector variables "\
                          "specified\n");
        break;
      }
      sscanf(argv[i],"cv %s", s);
      ev_cell[nv_cell] = GetElementVectorEvalProc(s);
      if (ev_cell[nv_cell]==NULL)
      {
        PrintErrorMessageF('E',"dataexplorer:","could not find vector "\
                           "eval proc %s\n",s);
        break;
      }
      if (sscanf(argv[i+1],"s %s", s) == 1)
      {
        strcpy(ev_cell_name[nv_cell],s);
        i++;
      }
      else
        strcpy(ev_cell_name[nv_cell],ev_cell[nv_cell]->v.name);
      nv_cell++;
      continue;
    }
  }

  if (ns==0 && nv==0 && ns_cell==0 && nv_cell==0)
    UserWrite("dataexplorer: no variables given, printing mesh data only\n");

  /* get file name and open output file	*/
  if (sscanf(argv[0],expandfmt(CONCAT3(" dataexplorer %",NAMELENSTR,"[ -~]")),
             filename)!=1)
  {
    PrintErrorMessage('E',"dataexplorer","could not read name of output file");
    return(PARAMERRORCODE);
  }
  pf = pfile_open(filename);
  if (pf==NULL) {
    PrintErrorMessage('E',"dataexplorer","could not open output file");
    return(PARAMERRORCODE);
  }

  /********************************/
  /* TITLE                                              */
  /********************************/

  time(&ltime);
  sprintf(it,"#\n# DataExplorer file written by UG\n");
  strcpy(item+ic,it); ic+=strlen(it);
  sprintf(it,"# Date: %s#\n",ctime(&ltime));
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;

  sprintf(it,"# Command: ");
  strcpy(item+ic,it); ic+=strlen(it);
  for ( i=0; i<argc; i++ )    {
    sprintf(it," %s",argv[i]);
    strcpy(item+ic,it); ic+=strlen(it);
  }
  sprintf(it,"\n#\n");
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;

  /********************************/
  /* compute sizes				*/
  /********************************/

  /* count vertices */
  /*
          numVertices = 0;
          for (k=0; k<=TOPLEVEL(mg); k++)
                  for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
                          numVertices++;
   */
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
      SETUSED(vx,0);
  numVertices = 0;
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el)) {
      if (!EstimateHere(el)) continue;
      for (i=0; i<CORNERS_OF_ELEM(el); i++) {
        vx = MYVERTEX(CORNER(el,i));
        if (USED(vx)) continue;
        SETUSED(vx,1);
        numVertices++;
      }
    }

  /* count surface elements */
  numElements = 0;
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el)) {
      if (!EstimateHere(el)) continue;
      numElements++;
    }

#ifdef ModelP
  gnumVertices = UG_GlobalSumINT(numVertices);
  gnumElements = UG_GlobalSumINT(numElements);
  LocallyUniqueIDs(mg);
  ov = get_offset(numVertices);
  oe = get_offset(numElements);
#else
  gnumVertices = numVertices;
  gnumElements = numElements;
  ov = oe = 0;
#endif

  /****************************************************************/
  /*	1. write vertex coordinates                                             */
  /****************************************************************/

  sprintf(it,"\n#\n# positions\n#\n");
  strcpy(item+ic,it); ic+=strlen(it);
  sprintf(it,"object 1 class array type float rank 1 shape %d items %d"\
          " data follows\n", DIM, gnumVertices);
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;

  /* unmark vertices */
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
      SETUSED(vx,0);

  /* allocate memory for mapping */
  heap = mg->theHeap;
  MarkTmpMem(heap, &key);
  Id2Position = (INT *)GetTmpMem(heap, numVertices*sizeof(INT), key);
  if (Id2Position == NULL) {
    ReleaseTmpMem(heap, key);
    UserWrite("dataexplorer: OOM\n");
    return CMDERRORCODE;
  }

  /* write vertex coordinates */
  counter=0;
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el)) {
      if (!EstimateHere(el)) continue;
      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        vx = MYVERTEX(CORNER(el,i));
        if (USED(vx)) continue;
        SETUSED(vx,1);

        /* map: vertex id -> position in coordinate list */
        Id2Position[ID(vx)] = counter;

        /* write the thing */
#ifdef __TWODIM__
        sprintf(it,"\t%g\t%g\n", clampf(XC(vx)), clampf(YC(vx)));
#else
        sprintf(it,"\t%g\t%g\t%g\n", clampf(XC(vx)), clampf(YC(vx)),
                clampf(ZC(vx)));
#endif
        pfile_tagged_puts(pf,it,counter+ov);
        counter++;
      }
    }
  pfile_sync(pf);


  /****************************************************************/
  /*	2. write connections										*/
  /****************************************************************/

  sprintf(it,"\n#\n# connections\n#\n");
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;
  sprintf(it,"object 2 class array type int rank 1 shape %d items %d"\
          " data follows\n", 4*(DIM-1), gnumElements);
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;

  counter = 0;
  for (k=0; k<=TOPLEVEL(mg); k++)         {
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (!EstimateHere(el)) continue;

      switch(CORNERS_OF_ELEM(el))
      {
      case 3 :
        sprintf(it,"\t%d\t%d\t%d\t%d\n",
                Id2Position[ID(MYVERTEX(CORNER(el,0)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,1)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,2)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,2)))]+ov);
        break;

      case 4 :
#ifdef __TWODIM__
        sprintf(it,"\t%d\t%d\t%d\t%d\n",
                Id2Position[ID(MYVERTEX(CORNER(el,0)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,1)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,3)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,2)))]+ov);
#else
        sprintf(it,"\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                Id2Position[ID(MYVERTEX(CORNER(el,0)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,1)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,2)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,2)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,3)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,3)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,3)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,3)))]+ov);
#endif
        break;

      case 5 :
        sprintf(it,"\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                Id2Position[ID(MYVERTEX(CORNER(el,0)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,1)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,3)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,2)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,4)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,4)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,4)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,4)))]+ov);
        break;

      case 6 :
        sprintf(it,"\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                Id2Position[ID(MYVERTEX(CORNER(el,0)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,1)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,2)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,2)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,3)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,4)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,5)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,5)))]+ov);
        break;

      case 8 :
        sprintf(it,"\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                Id2Position[ID(MYVERTEX(CORNER(el,0)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,1)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,3)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,2)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,4)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,5)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,7)))]+ov,
                Id2Position[ID(MYVERTEX(CORNER(el,6)))]+ov);
        break;
      }
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_tagged_puts(pf,item,counter+oe); ic=0;
      counter++;
    }
  }
  pfile_sync(pf);

#ifdef __TWODIM__
  sprintf(it,"attribute \"element type\" string \"quads\"\n");
#else
  sprintf(it,"attribute \"element type\" string \"cubes\"\n");
#endif
  strcpy(item+ic,it); ic+=strlen(it);
  sprintf(it,"attribute \"ref\" string \"positions\"\n\n");
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;

  /* delete map array */
  ReleaseTmpMem(heap, key);

  /****************************************************************/
  /*	3. write node data											*/
  /****************************************************************/

  blocks = 1;

  /* write all scalar node data */
  for (v=0; v<ns; v++)
  {
    pre     = es[v]->PreprocessProc;
    if (pre!=NULL) pre(es_name[v],mg);

    sprintf(it,"#\n# data block %d\n#\n", blocks);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"object %d class array type float rank 0"\
            " items %d data follows\n", blocks+2, gnumVertices);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    for (k=0; k<=TOPLEVEL(mg); k++)
      for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
        SETUSED(vx,0);

    counter=0;
    for (k=0; k<=TOPLEVEL(mg); k++) {
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if (!EstimateHere(el)) continue;
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
          CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
        {
          vx = MYVERTEX(CORNER(el,i));
          if (USED(vx)) continue;
          SETUSED(vx,1);

          /* get local coordinate of corner */
          LocalCornerCoordinates(DIM,TAG(el),i,local);
          for (j=0; j<DIM; j++) LocalCoord[j] = local[j];

          /* scalar components */
          eval_s = es[v]->EvalProc;
          value = eval_s(el,(const DOUBLE **)CornersCoord,LocalCoord);
          sprintf(it,"\t%g\n",clampf(value));
          strcpy(item+ic,it); ic+=strlen(it);
          pfile_tagged_puts(pf,item,counter+ov); ic=0;
          counter++;
        }
      }
    }
    pfile_sync(pf);

    sprintf(it,"\nobject \"data%d\" class field\n", blocks);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"positions\" value 1\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"connections\" value 2\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"data\" value %d\n", blocks+2);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    blocks++;
  }

  /* write all vector node data */
  for (v=0; v<nv; v++)
  {
    pre     = ev[v]->PreprocessProc;
    if (pre!=NULL) pre(ev_name[v],mg);

    sprintf(it,"#\n# data block %d\n#\n", blocks);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"object %d class array type float rank 1 shape %d items %d"\
            " data follows\n", blocks+2, DIM, gnumVertices);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    for (k=0; k<=TOPLEVEL(mg); k++)
      for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
        SETUSED(vx,0);
    counter=0;
    for (k=0; k<=TOPLEVEL(mg); k++)         {
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if (!EstimateHere(el)) continue;
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
          CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
        {
          vx = MYVERTEX(CORNER(el,i));
          if (USED(vx)) continue;
          SETUSED(vx,1);

          /* get local coordinate of corner */
          LocalCornerCoordinates(DIM,TAG(el),i,local);
          for (j=0; j<DIM; j++) LocalCoord[j] = local[j];

          /* vector components */
          eval_v = ev[v]->EvalProc;
          eval_v(el,(const DOUBLE **)CornersCoord,LocalCoord,vval);
#ifdef __TWODIM__
          sprintf(it,"\t%g\t%g\n",clampf(vval[0]),clampf(vval[1]));
#else
          sprintf(it,"\t%g\t%g\t%g\n",clampf(vval[0]),clampf(vval[1]),
                  clampf(vval[2]));
#endif
          strcpy(item+ic,it); ic+=strlen(it);
          pfile_tagged_puts(pf,item,counter+ov); ic=0;
          counter++;
        }
      }
    }
    pfile_sync(pf);

    sprintf(it,"\nobject \"data%d\" class field\n", blocks);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"positions\" value 1\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"connections\" value 2\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"data\" value %d\n", blocks+2);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    blocks++;
  }

  /****************************************************************/
  /*	3. element data												*/
  /****************************************************************/

  /* write all scalar element data */
  for (v=0; v<ns_cell; v++)
  {
    pre     = es_cell[v]->PreprocessProc;
    if (pre!=NULL) pre(es_cell_name[v],mg);

    sprintf(it,"#\n# data block %d\n#\n", blocks);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"object %d class array type float rank 0"\
            " items %d data follows\n", blocks+2, gnumElements);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    counter=0;
    for (k=0; k<=TOPLEVEL(mg); k++) {
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if (!EstimateHere(el)) continue;
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
          CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));
        for (j=0; j<DIM; j++) LocalCoord[j] = 0.0;
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
        {
          /* get local coordinate of corner	*/
          LocalCornerCoordinates(DIM,TAG(el),i,local);
          for (j=0; j<DIM; j++) LocalCoord[j] += local[j];
        }
        for (j=0; j<DIM; j++) LocalCoord[j] /= ((DOUBLE)CORNERS_OF_ELEM(el));

        /* scalar component	*/
        eval_s = es_cell[v]->EvalProc;
        value = eval_s(el,(const DOUBLE **)CornersCoord,LocalCoord);
        sprintf(it,"\t%g\n",clampf(value));
        strcpy(item+ic,it); ic+=strlen(it);
        pfile_tagged_puts(pf,item,counter+oe); ic=0;
        counter++;
      }
    }
    pfile_sync(pf);

    sprintf(it,"attribute \"dep\" string \"connections\"\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"\nobject \"data%d\" class field\n", blocks);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"positions\" value 1\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"connections\" value 2\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"data\" value %d\n", blocks+2);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    blocks++;
  }

  /* write all element vector data */
  for (v=0; v<nv_cell; v++)
  {
    pre     = ev_cell[v]->PreprocessProc;
    if (pre!=NULL) pre(ev_cell_name[v],mg);

    sprintf(it,"#\n# data block %d\n#\n", blocks);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"object %d class array type float rank 1 shape %d items %d"\
            " data follows\n", blocks+2, DIM, gnumElements);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    counter=0;
    for (k=0; k<=TOPLEVEL(mg); k++) {
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if (!EstimateHere(el)) continue;
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
          CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));
        for (j=0; j<DIM; j++) LocalCoord[j] = 0.0;
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
        {
          LocalCornerCoordinates(DIM,TAG(el),i,local);
          for (j=0; j<DIM; j++) LocalCoord[j] += local[j];
        }
        for (j=0; j<DIM; j++) LocalCoord[j] /= ((DOUBLE)CORNERS_OF_ELEM(el));
        eval_v = ev_cell[v]->EvalProc;
        eval_v(el,(const DOUBLE **)CornersCoord,LocalCoord,vval);
#ifdef __TWODIM__
        sprintf(it,"\t%g\t%g\n",clampf(vval[0]),clampf(vval[1]));
#else
        sprintf(it,"\t%g\t%g\t%g\n",clampf(vval[0]),clampf(vval[1]),
                clampf(vval[2]));
#endif
        strcpy(item+ic,it); ic+=strlen(it);

        pfile_tagged_puts(pf,item,counter+oe); ic=0;
        counter++;
      }
    }
    pfile_sync(pf);

    sprintf(it,"attribute \"dep\" string \"positions\"\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"\nobject \"data%d\" class field\n", blocks);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"positions\" value 1\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"connections\" value 2\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"data\" value %d", blocks+2);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    blocks++;
  }

  pfile_close(pf);

  return(OKCODE);
}

/****************************************************************************/
/*																			*/
/* Function:  InitDataExplorer												*/
/*																			*/
/* Purpose:   register all formats for porous media library					*/
/*																			*/
/* Input:	 void															*/
/*																			*/
/* Output:	  INT 0: ok														*/
/*				  else line number where error occured						*/
/*																			*/
/****************************************************************************/

INT InitDataExplorer (void)
{
  if (CreateCommand("dataexplorer",DataExplorerCommand)==NULL) return (__LINE__);

  return(0);
}
