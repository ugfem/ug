// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  dataexplorer.c												*/
/*																			*/
/* Purpose:   dataexplorer output			                                                                */
/*																			*/
/* Author:	  Volker Reichenberger											*/
/*			  IWR Technische Simulation				                                                */
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  69120 Heidelberg												*/
/*			  email: volker.reichenberger@iwr.uni-heidelberg.de				*/
/*																			*/
/* History:   23.03.2000 begin												*/
/*																			*/
/* Remarks:   based on avs.c												*/
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

#define MAXVARIABLES    20                      /* max number of eval procs				*/
#define VALUES_PER_LINE 10                      /* number of data values per line		*/

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

/* data for CVS	*/
static char RCS_ID("$Header$",UG_RCS_STRING);


#ifdef ModelP
static INT get_offset (INT n)
{
  INT i,subtreesize[50],sum,offset;

  /* get number of objects downtree	*/
  sum = n;
  for (i=0; i<degree; i++) {
    GetConcentrate(i,subtreesize+i,sizeof(INT));
    sum += subtreesize[i];
  }

  /* get offset	*/
  if (me==master)
  {
    offset = 0;
  }
  else
  {
    Concentrate(&sum,sizeof(INT));
    GetSpread(&offset,sizeof(INT));
  }

  /* send offsets for downtree nodes	*/
  sum = offset+n;
  for (i=0; i<degree; i++) {
    Spread(i,&sum,sizeof(INT));
    sum += subtreesize[i];
  }

  return(offset);
}

static INT GloballyUniqueIDs (MULTIGRID *theMG)
{
  GRID *theGrid;
  VERTEX *theVertex;
  NODE *theNode;
  ELEMENT *theElement;
  INT nv,nn,ne;
  INT ov,on,oe;
  int k,j;

  nv = ne = nn = 0;

  j = theMG->topLevel;
  for (k=0; k<=j; k++)
  {
    theGrid = GRID_ON_LEVEL(theMG,k);

    /* vertices	*/
    for (theVertex=FIRSTVERTEX(theGrid); theVertex!=NULL; theVertex=SUCCV(theVertex))
      ID(theVertex) = nv++;

    /* nodes	*/
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
      ID(theNode) = nn++;

    /* elements	*/
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      ID(theElement) = ne++;
  }

  theMG->vertIdCounter = nv;
  theMG->nodeIdCounter = nn;
  theMG->elemIdCounter = ne;

  ov = get_offset(nv);
  on = get_offset(nn);
  oe = get_offset(ne);

  for (k=0; k<=j; k++)
  {
    theGrid = GRID_ON_LEVEL(theMG,k);

    /* vertices	*/
    for (theVertex=FIRSTVERTEX(theGrid); theVertex!=NULL; theVertex=SUCCV(theVertex))
      ID(theVertex) += ov;

    /* nodes	*/
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
      ID(theNode) += on;

    /* elements	*/
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      ID(theElement) += oe;
  }

  return(0);
}

static INT DataExplorer_GlobalSumINT (INT i)
{
  int l;
  INT n;

  for (l=degree-1; l>=0; l--)
  {
    GetConcentrate(l,&n,sizeof(INT));
    i += n;
  }
  Concentrate(&i,sizeof(INT));
  Broadcast(&i,sizeof(INT));
  return(i);
}
#endif

/****************************************************************************/
/*																			*/
/* definition of functions													*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*D
   dataexplorer - file output in DataExplorer format

   DESCRIPTION:
   The dataexplorer command stores grid and multiple vector or scalar grid functions
   in DataExplorer readable format to file.

   'dataexplorer <filename> [$scale <factor>] [$zcoord <nep> $s <vd>]
                                   [$ns <nep> $s <vd>]* [$nv <nep> $s <vd>]*
                                   [$es <eep> $s <vd>]* [$ev <eep> $s <vd>]*
                                   [$noco]'

   .  $scale~<factor>	- ?
   .  $zcoord...		- plot function to be drawn on z-coordinate (2D grids only)
   .  $ns...			- plot function for scalar nodal values
   .  $nv...			- plot function for vector nodal values
   .  $cs...			- plot function for scalar element values
   .  $cv...			- plot function for vector element values
   .  $noco			- don't include connections in output

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
  char filename[NAMESIZE];              /* file name for output file				*/
  PFILE *pf;                                            /* the output file pointer					*/

  EVALUES *zcoord;                              /* use scalar as z coordinate in 2D only	*/
  char zcoord_name[NAMESIZE];           /* name for eval functions                              */

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
  INT numNodes;                                 /* number of data points locally			*/
  INT numElements;                              /* number of elements locally				*/
        #ifdef __TWODIM__
  INT numTriangles;                             /* number of elements locally				*/
  INT numQuads;                                 /* number of elements locally				*/
        #else
  INT numTets;                                  /* number of tetrahedons locally			*/
  INT numPyramids;                              /* number of pyramids locally				*/
  INT numPrisms;                                /* number of prisms locally					*/
  INT numHexes;                                 /* number of hexahedrons locally			*/
        #endif
  INT gnumNodes;                                /* number of data points globally			*/
  INT gnumElements;                             /* number of elements globallay				*/
        #ifdef __TWODIM__
  INT gnumTriangles;                            /* number of triangles globallay			*/
  INT gnumQuads;                                /* number of quadrilaterals globallay		*/
        #else
  INT gnumTets;                                 /* number of tetrahedons globallay			*/
  INT gnumPyramids;                             /* number of pyramids globallay				*/
  INT gnumPrisms;                               /* number of prisms globallay				*/
  INT gnumHexes;                                /* number of hexahedrons globallay			*/
        #endif
  PreprocessingProcPtr pre;             /* pointer to prepare function				*/
  ElementEvalProcPtr eval_s;            /* pointer to scalar evaluation function	*/
  ElementVectorProcPtr eval_v;      /* pointer to vector evaluation function	*/
  DOUBLE *CornersCoord[MAX_CORNERS_OF_ELEM];       /* pointers to coordinates   */
  DOUBLE LocalCoord[DIM];               /* is one of the corners local coordinates	*/
  DOUBLE local[DIM];                            /* local coordinate in DOUBLE				*/
  DOUBLE value;                                 /* returned by user eval proc				*/
  DOUBLE x,y,z;                                 /* scalar values							*/
  DOUBLE vval[DIM];                             /* result of vector evaluation function		*/
  DOUBLE scale;
  VECDATA_DESC *displacement;
  const SHORT *comp;
  time_t ltime;
  double scaling=1.0;
  INT oe,on,n;
        #ifdef __TWODIM__
  INT oTriangles, oQuads;
        #else
  INT oTets, oPyramids, oPrisms, oHexes;
        #endif
  INT tag;
  INT datablock, nitem;
  INT doConnections = TRUE;


        #ifdef __THREEDIM__
  UserWrite("3D output for data explorer doesn't work yet.\n");
  return 1;
        #endif

  /* get current multigrid	*/
  mg = GetCurrentMultigrid();
  if (mg==NULL)
  {
    PrintErrorMessage('W',"dataexplorer","no multigrid open\n");
    return (OKCODE);
  }

  /* scan options	*/
  ns = nv = ns_cell = nv_cell = 0;
  zcoord = NULL;
  displacement = ReadArgvVecDesc(mg,"displacement",argc,argv);
  if (displacement != NULL) {
    comp = VD_ncmp_cmpptr_of_otype(displacement,NODEVEC,&n);
    if (n < DIM)
      return(CMDERRORCODE);
    if (ReadArgvDOUBLE("scale",&scale,argc,argv))
      scale = 1.0;
  }
  for(i=1; i<argc; i++)
  {
    if (strncmp(argv[i],"scale",5)==0) {
      sscanf(argv[i],"scale %s", s);
      sscanf(s,"%lg",&scaling);
      continue;
    }

    if (strncmp(argv[i],"zcoord",6)==0) {
      sscanf(argv[i],"zcoord %s", s);
      zcoord = GetElementValueEvalProc(s);
      if (zcoord==NULL)
      {
        PrintErrorMessageF('E',"dataexplorer:","could not find scalar eval proc %s\n",s);
        break;
      }
      if (sscanf(argv[i+1],"s %s", s) == 1)
      {
        strcpy(zcoord_name,s);
        i++;
      }
      else
        strcpy(zcoord_name,zcoord->v.name);
      continue;
    }

    if (strncmp(argv[i],"ns",2)==0) {
      if (ns>=MAXVARIABLES)
      {
        PrintErrorMessage('E',"dataexplorer:","too many scalar variables specified\n");
        break;
      }
      sscanf(argv[i],"ns %s", s);
      es[ns] = GetElementValueEvalProc(s);
      if (es[ns]==NULL)
      {
        PrintErrorMessageF('E',"dataexplorer:","could not find scalar eval proc %s\n",s);
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
        PrintErrorMessage('E',"dataexplorer:","too many vector variables specified\n");
        break;
      }
      sscanf(argv[i],"nv %s", s);
      ev[nv] = GetElementVectorEvalProc(s);
      if (ev[nv]==NULL)
      {
        PrintErrorMessageF('E',"dataexplorer:","could not find vector eval proc %s\n",s);
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
        PrintErrorMessage('E',"dataexplorer:","too many scalar variables specified\n");
        break;
      }
      sscanf(argv[i],"cs %s", s);
      es_cell[ns_cell] = GetElementValueEvalProc(s);
      if (es_cell[ns_cell]==NULL)
      {
        PrintErrorMessageF('E',"dataexplorer:","could not find scalar eval proc %s\n",s);
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
        PrintErrorMessage('E',"dataexplorer:","too many vector variables specified\n");
        break;
      }
      sscanf(argv[i],"cv %s", s);
      ev_cell[nv_cell] = GetElementVectorEvalProc(s);
      if (ev_cell[nv_cell]==NULL)
      {
        PrintErrorMessageF('E',"dataexplorer:","could not find vector eval proc %s\n",s);
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

    if (strncmp(argv[i],"noco",4)==0) {
      doConnections = FALSE;
      continue;
    }

  }
  if (ns==0 && nv==0) UserWrite("dataexplorer: no variables given, printing mesh data only\n");

  /* get file name and open output file	*/
  if (sscanf(argv[0],expandfmt(CONCAT3(" dataexplorer %",NAMELENSTR,"[ -~]")),filename)!=1)
  {
    PrintErrorMessage('E',"dataexplorer","could not read name of output file");
    return(PARAMERRORCODE);
  }
  pf = pfile_open(filename);
  if (pf==NULL) return(PARAMERRORCODE);

  /********************************/
  /* TITLE                                              */
  /********************************/

  time(&ltime);
  sprintf(it,"# DataExplorer data produced by UG\n");
  strcpy(item+ic,it); ic+=strlen(it);
  sprintf(it,"# with multigrid: %s\n",mg->v.name);
  strcpy(item+ic,it); ic+=strlen(it);
  sprintf(it,"# date:	  %s\n",ctime(&ltime));
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;

  sprintf(it,"#\n# produced with	");
  strcpy(item+ic,it); ic+=strlen(it);
  for ( i=0; i<argc; i++ )    {
    sprintf(it," %s",argv[i]);
    strcpy(item+ic,it); ic+=strlen(it);
  }
  sprintf(it,"\n");
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;

  /********************************/
  /* compute sizes				*/
  /********************************/

  /* clear USED flag in vertices on all levels	*/
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
      SETUSED(vx,0);

  /* run thru all levels of elements and set index	*/
  numNodes = numElements = 0;
        #ifdef __TWODIM__
  numTriangles = numQuads = 0;
        #else
  numTets = numPyramids = numPrisms = numHexes = 0;
        #endif
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (!EstimateHere(el)) continue;                                  /* process finest level elements only	*/
                        #ifdef __TWODIM__
      if ( CORNERS_OF_ELEM(el)==3 ) numTriangles++;
      else numQuads++;
                        #else
      if ( CORNERS_OF_ELEM(el)==4 ) numTets++;
      else if ( CORNERS_OF_ELEM(el)==5 ) numPyramids++;
      else if ( CORNERS_OF_ELEM(el)==6 ) numPrisms++;
      else numHexes++;
                        #endif
      numElements++;
      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        vx = MYVERTEX(CORNER(el,i));
        if (USED(vx)) continue;                                 /* we have this one already	*/

        ++numNodes;                                                             /* number of data points, begins with 1 !	*/
        SETUSED(vx,1);                                                  /* tag vector as visited	*/
      }
    }

        #ifdef ModelP
  gnumNodes = DataExplorer_GlobalSumINT(numNodes);
  gnumElements = DataExplorer_GlobalSumINT(numElements);
  gnumTriangles = DataExplorer_GlobalSumINT(numTriangles);
  gnumQuads = DataExplorer_GlobalSumINT(numQuads);
  on=get_offset(numNodes);
  oe=get_offset(numElements);
  GloballyUniqueIDs(mg);       /* renumber objects	*/
        #else
  gnumNodes = numNodes;
  gnumElements = numElements;
    #ifdef __TWODIM__
  gnumTriangles = numTriangles;
  gnumQuads = numQuads;
    #endif
  oe=on=0;
        #endif

  /****************************************************************/
  /*	1. write "positions". i.e. vertex coordinates				*/
  /****************************************************************/

  sprintf(it,"# positions in %dD\n",DIM);
  strcpy(item+ic,it); ic+=strlen(it);
  sprintf(it,"object 1 class array type float rank 1 shape %d items %d data follows\n",
          3,gnumNodes);
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;

  if (DIM==2 && zcoord!=NULL)
  {
    pre     = zcoord->PreprocessProc;
    if (pre!=NULL) pre(zcoord_name,mg);
  }

  /* clear USED flag in vertices on all levels	*/
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
      SETUSED(vx,0);

  counter=0;
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (!EstimateHere(el)) continue;                                  /* process finest level elements only	*/
      for (i=0; i<CORNERS_OF_ELEM(el); i++)
        CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));                         /* x,y,z of corners	*/
      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        vx = MYVERTEX(CORNER(el,i));
        if (USED(vx)) continue;                                 /* we have this one already	*/
        SETUSED(vx,1);                                                  /* tag vector as visited	*/

        /* get local coordinate of corner	*/
        LocalCornerCoordinates(DIM,TAG(el),i,local);
        for (j=0; j<DIM; j++) LocalCoord[j] = local[j];

        x = XC(vx);
        y = YC(vx);
        if (DIM==3)
        {
          z = ZC(vx);
        }
        else
        {
          if (zcoord!=NULL)
          {
            eval_s = zcoord->EvalProc;
            z = scaling*eval_s(el,(const DOUBLE **)CornersCoord,LocalCoord);
          }
          else
            z = 0.0;
        }
        if (displacement != NULL) {
          x += scale * VVALUE(NVECTOR(CORNER(el,i)),comp[0]);
          y += scale * VVALUE(NVECTOR(CORNER(el,i)),comp[1]);
                                        #ifdef __THREEDIM__
          z += scale * VVALUE(NVECTOR(CORNER(el,i)),comp[2]);
                                        #endif
        }
        sprintf(it,"\t%g\t%g\t%g\n",x,y,z);
        pfile_tagged_puts(pf,it,counter+on);
        counter++;
      }
    }

  pfile_sync(pf);       /* end of segment	*/

  /* TODO: Write element center positions, if needed by cell data. */

  /****************************************************************/
  /*	2. cell connectivity										*/
  /****************************************************************/
  if ( doConnections==TRUE )      {

    sprintf(it,"# connections\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    counter=0;
    if (DIM==2)     {
      for ( tag=3; tag<=4; tag++ )    {
                #ifdef __TWODIM__
        sprintf(it,"object %d class array type float rank 1 shape %d items %d data follows\n",
                tag,3,(tag==3) ? gnumTriangles : gnumQuads);
                #endif
        strcpy(item+ic,it); ic+=strlen(it);
        pfile_master_puts(pf,item); ic=0;

        for (k=0; k<=TOPLEVEL(mg); k++)         {
          for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
          {
            if (!EstimateHere(el)) continue;                                                    /* process finest level elements only	*/

            /* cell number and material id	*/
            if  (  tag==3 && CORNERS_OF_ELEM(el)==3 )       {
              sprintf(it,"\t%d\t%d\t%d\n",
                      ID(MYVERTEX(CORNER(el,0))),
                      ID(MYVERTEX(CORNER(el,1))),
                      ID(MYVERTEX(CORNER(el,2))));
              strcpy(item+ic,it); ic+=strlen(it);
              pfile_tagged_puts(pf,item,counter+oe); ic=0;
              counter++;
            }

            if  (  tag==4 && CORNERS_OF_ELEM(el)==4 )       {
              sprintf(it,"\t%d\t%d\t%d\t%d\n",
                      ID(MYVERTEX(CORNER(el,0))),
                      ID(MYVERTEX(CORNER(el,1))),
                      ID(MYVERTEX(CORNER(el,2))),
                      ID(MYVERTEX(CORNER(el,3))));
              strcpy(item+ic,it); ic+=strlen(it);
              pfile_tagged_puts(pf,item,counter+oe); ic=0;
              counter++;
            }
          }
        }

        if ( tag==3 )   {
          sprintf(it,"attribute \"element type\" string \"triangles\"\n");
          strcpy(item+ic,it); ic+=strlen(it);
          sprintf(it,"attribute \"ref\" string \"positions\"\n\n");
          strcpy(item+ic,it); ic+=strlen(it);
          pfile_master_puts(pf,item); ic=0;
        }
        if ( tag==4 )   {
          sprintf(it,"attribute \"element type\" string \"cubes \"\n");
          strcpy(item+ic,it); ic+=strlen(it);
          sprintf(it,"attribute \"ref\" string \"positions\"\n\n");
          strcpy(item+ic,it); ic+=strlen(it);
          pfile_master_puts(pf,item); ic=0;
        }

      }
    }

    pfile_sync(pf);             /* end of segment	*/
  }

  /****************************************************************/
  /*	3. node data												*/
  /****************************************************************/
  datablock = 10;
  if (1+ns+nv>0)
  {
    /* execute all prepare functions */
    for (v=0; v<ns; v++)
    {
      pre = es[v]->PreprocessProc;

      /* execute prepare function */
      if (pre!=NULL) pre(es_name[v],mg);
    }

    /* All scalar variables										*/
    for (v=0; v<ns; v++) {
      /* Write comment about what we are about to write		*/
      sprintf(it,"# Data block \"NodeScalar%s \"\n",ENVITEM_NAME(es[v]));
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"object \"NodeScalar%s\" class array type float rank 0 items %d data follows\n",ENVITEM_NAME(es[v]), gnumNodes);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      /* now the scalar data.									*/
      /* four values in each line, but this is arbitrary.		*/
      /* clear USED flag in vertices on all levels			*/
      for (k=0; k<=TOPLEVEL(mg); k++)
        for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
          SETUSED(vx,0);
      counter=0;  nitem = 0;
      for (k=0; k<=TOPLEVEL(mg); k++)     {
        for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
        {
          if (!EstimateHere(el)) continue;                      /* process finest level elements only */
          for (i=0; i<CORNERS_OF_ELEM(el); i++)
            CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));             /* x,y,z of corners */
          for (i=0; i<CORNERS_OF_ELEM(el); i++)
          {
            vx = MYVERTEX(CORNER(el,i));
            if (USED(vx)) continue;                     /* we have this one already */
            SETUSED(vx,1);                              /* tag vector as visited */

            /* get local coordinate of corner */
            LocalCornerCoordinates(DIM,TAG(el),i,local);
            for (j=0; j<DIM; j++) LocalCoord[j] = local[j];

            /* scalar components */
            eval_s = es[v]->EvalProc;
            value = eval_s(el,(const DOUBLE **)CornersCoord,LocalCoord);
            sprintf(it,"\t%g",value);
            strcpy(item+ic,it); ic+=strlen(it);

            if ( ++nitem>4 ) {
              sprintf(it,"\n"); nitem=0;
              strcpy(item+ic,it); ic+=strlen(it);
            }
            pfile_tagged_puts(pf,item,counter+on); ic=0;
            counter++;
          }
        }
      }
      if ( nitem!=0 ) {
        sprintf(it,"\n"); strcpy(item+ic,it); ic+=strlen(it);
        pfile_master_puts(pf,item); ic=0;
      }

      sprintf(it,"attribute \"dep\" string \"positions\"\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      /*  And now the Field	*/
      sprintf(it,"\n# The field\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"object \"FieldNode%s\" class field\n",
              ENVITEM_NAME(es[v]));
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"positions\" value 1\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"connections\" value 3\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"connections\" value 4\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"data\" value \"NodeScalar%s\"\n\n\n",ENVITEM_NAME(es[v]));
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      pfile_sync(pf);       /* end of segment	*/

    }

    for (v=0; v<nv; v++)
    {
      pre = ev[v]->PreprocessProc;

      /* execute prepare function */
      if (pre!=NULL) pre(ev_name[v],mg);
    }

    /* All vector variables										*/
    for (v=0; v<nv; v++) {
      /* Write comment about what we are about to write		*/
      sprintf(it,"# Data block \"NodeVector%s\"\n",ENVITEM_NAME(ev[i]));
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"object \"NodeVector%s\" class array type float rank 1 shape %d items %d data follows\n",
              ENVITEM_NAME(ev[i]), DIM, gnumNodes);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      /* now the vector data.									*/
      /* four values in each line, but this is arbitrary.		*/
      /* clear USED flag in vertices on all levels			*/
      for (k=0; k<=TOPLEVEL(mg); k++)
        for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
          SETUSED(vx,0);
      counter=0;
      for (k=0; k<=TOPLEVEL(mg); k++)     {
        for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
        {
          if (!EstimateHere(el)) continue;                      /* process finest level elements only */
          for (i=0; i<CORNERS_OF_ELEM(el); i++)
            CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));             /* x,y,z of corners */
          for (i=0; i<CORNERS_OF_ELEM(el); i++)
          {
            vx = MYVERTEX(CORNER(el,i));
            if (USED(vx)) continue;                     /* we have this one already */
            SETUSED(vx,1);                              /* tag vector as visited */

            /* get local coordinate of corner */
            LocalCornerCoordinates(DIM,TAG(el),i,local);
            for (j=0; j<DIM; j++) LocalCoord[j] = local[j];

            /* vector components */
            eval_v = ev[v]->EvalProc;
            eval_v(el,(const DOUBLE **)CornersCoord,LocalCoord,vval);
            if (DIM==2)
              sprintf(it,"\t%g\t%g\n",vval[0],vval[1]);
            else
              sprintf(it,"\t%g\t%g\t%g\n",vval[0],vval[1],vval[2]);
            strcpy(item+ic,it); ic+=strlen(it);

            pfile_tagged_puts(pf,item,counter+on); ic=0;
            counter++;
          }
        }
      }
      if ( nitem!=0 ) {
        sprintf(it,"\n"); strcpy(item+ic,it); ic+=strlen(it);
        pfile_master_puts(pf,item); ic=0;
      }
      sprintf(it,"attribute \"dep\" string \"positions\"\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      /*  And now the Field	*/
      sprintf(it,"\n# The field\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"object \"FieldNode%s\" class field\n",
              ENVITEM_NAME(es[v]));
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"positions\" value 1\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"connections\" value 3\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"connections\" value 4\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"data\" value \"NodeScalar%s\"\n\n\n",ENVITEM_NAME(es[v]));
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;
      pfile_sync(pf);       /* end of segment	*/
    }
  }

  /****************************************************************/
  /*	3. element data												*/
  /****************************************************************/

  if (ns_cell+nv_cell>0)
  {
    /* execute prepare functions	*/
    for (v=0; v<ns_cell; v++)
    {
      pre     = es_cell[v]->PreprocessProc;

      /* execute prepare function	*/
      if (pre!=NULL) pre(es_cell_name[v],mg);
    }

    for (v=0; v<ns_cell; v++)
    {
      /* Write comment about what we are about to write		*/
      sprintf(it,"# Data block \"ElementScalar%s \"\n",ENVITEM_NAME(es_cell[v]));
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"object \"ElementScalar%s\" class array type float rank 0 items %d data follows\n",
              ENVITEM_NAME(es_cell[v]), gnumNodes);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      /* now the scalar data                                                                    */
      counter=0; nitem=0;
      for (k=0; k<=TOPLEVEL(mg); k++)     {
        for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
        {
          if (!EstimateHere(el)) continue;                      /* process finest level elements only	*/
          for (i=0; i<CORNERS_OF_ELEM(el); i++)
            CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));             /* x,y,z of corners	*/
          /* compute center in local coordinates	*/
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
          sprintf(it,"\t%g",value);
          strcpy(item+ic,it); ic+=strlen(it);

          if ( ++nitem>3 ) {
            sprintf(it,"\n");
            strcpy(item+ic,it); ic+=strlen(it);
            nitem=0;
          }
          pfile_tagged_puts(pf,item,counter+oe); ic=0;
          counter++;
        }
      }

      if ( nitem!=0 ) {
        sprintf(it,"\n"); strcpy(item+ic,it); ic+=strlen(it);
        pfile_master_puts(pf,item); ic=0;
      }
      sprintf(it,"attribute \"dep\" string \"positions\"\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      /*  And now the Field	*/
      sprintf(it,"\n# The field\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"object \"FieldElement%s\" class field\n",
              ENVITEM_NAME(es[v]));
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"positions\" value 1\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"connections\" value 3\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"connections\" value 4\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"data\" value \"ElementScalar%s\"\n\n\n",ENVITEM_NAME(es_cell[v]));
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      pfile_sync(pf);       /* end of segment	*/
    }

    /* execute prepare functions	*/
    for (v=0; v<nv_cell; v++)
    {
      pre     = ev_cell[v]->PreprocessProc;

      /* execute prepare function	*/
      if (pre!=NULL) pre(ev_cell_name[v],mg);
    }


    for (v=0; v<nv_cell; v++)
    {
      /* Write comment about what we are about to write		*/
      sprintf(it,"# Data block \"ElementScalar%s \"\n",ENVITEM_NAME(es_cell[v]));
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"object \"ElementScalar%s\" class array type float rank 0 items %d data follows\n",
              ENVITEM_NAME(es_cell[v]), gnumNodes);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      /* now the scalar data                                                                    */
      counter=0; nitem=0;
      for (k=0; k<=TOPLEVEL(mg); k++)     {
        for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
        {
          if (!EstimateHere(el)) continue;                      /* process finest level elements only	*/
          for (i=0; i<CORNERS_OF_ELEM(el); i++)
            CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));             /* x,y,z of corners	*/
          /* compute center in local coordinates	*/
          for (j=0; j<DIM; j++) LocalCoord[j] = 0.0;
          for (i=0; i<CORNERS_OF_ELEM(el); i++)
          {
            /* get local coordinate of corner	*/
            LocalCornerCoordinates(DIM,TAG(el),i,local);
            for (j=0; j<DIM; j++) LocalCoord[j] += local[j];
          }
          for (j=0; j<DIM; j++) LocalCoord[j] /= ((DOUBLE)CORNERS_OF_ELEM(el));

          /* scalar component	*/
          eval_v = ev_cell[v]->EvalProc;
          eval_v(el,(const DOUBLE **)CornersCoord,LocalCoord,vval);
          if (DIM==2)
            sprintf(it,"\t%g\t%g\n",vval[0],vval[1]);
          else
            sprintf(it,"\t%g\t%g\t%g\n",vval[0],vval[1],vval[2]);
          strcpy(item+ic,it); ic+=strlen(it);

          pfile_tagged_puts(pf,item,counter+oe); ic=0;
          counter++;
        }
      }

      sprintf(it,"attribute \"dep\" string \"positions\"\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      /*  And now the Field	*/
      sprintf(it,"\n# The field\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"object \"FieldElement%s\" class field\n",
              ENVITEM_NAME(es[v]));
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"positions\" value 1\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"connections\" value 3\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"connections\" value 4\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"\tcomponent \"data\" value \"ElementVector%s\"\n\n\n",ENVITEM_NAME(es_cell[v]));
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      pfile_sync(pf);       /* end of segment	*/
    }

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
