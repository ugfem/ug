// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
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
#include "elements.h"
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
#include "udm.h"

#include "dataexplorer.h"
#include "cw.h"

#ifdef __cplusplus
#ifdef __TWODIM__
using namespace UG2d;
#else
using namespace UG3d;
#endif
#endif

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

static void LocallyUniqueIDs (MULTIGRID *mg)
{
  VERTEX *vx;
  ELEMENT *el;
  INT i, k, nv;

  for (k=0; k<=TOPLEVEL(mg); k++)
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
      SETUSED(vx,0);

  nv=0;
  /* first set master ids */
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el)) {
      if (!EstimateHere(el)) continue;
      for (i=0; i<CORNERS_OF_ELEM(el); i++) {
        vx = MYVERTEX(CORNER(el,i));
        if (USED(vx)) continue;
        if (DDD_InfoPriority(PARHDR(CORNER(el,i)))!=PrioMaster) continue;
        SETUSED(vx,1);
        ID(vx) = nv++;
      }
    }

  /* next set remaining ids */
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

static int Gather_VertexID (DDD_OBJ obj, void *data)
{
  VECTOR *pv = (VECTOR *)obj;

  ((INT *)data)[0] = ID(MYVERTEX((NODE *)VOBJECT(pv)));

  return (NUM_OK);
}
static int Scatter_VertexID (DDD_OBJ obj, void *data)
{
  VECTOR *pv = (VECTOR *)obj;

  ID(MYVERTEX((NODE *)VOBJECT(pv))) = ((INT *)data)[0];

  return (NUM_OK);
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
   grid functions in DataExplorer readable format to a header and a
   data file. The data file is written as a text file (default) or in
   binary format. In order to save disk space it is possible to write
   only the data and not the grid at later time steps (if the grid is
   not changing).

   'dataexplorer <filename> [$ns <nep> $s <vd>]* [$nv <nep> $s <vd>]*
                                            [$cs <eep> $s <vd>]* [$cv <eep> $s <vd>]*
                                                        [$b 0|1|2]* [$bin]* [$fgrid]* [$s [+|-]id]*'

   .  $ns...			- plot function for scalar nodal values
   .  $nv...			- plot function for vector nodal values
   .  $cs...			- plot function for scalar element values
   .  $cv...			- plot function for vector element values
   .  $b...            - write boundary data 0=no (default) | 1=inner | 2=all
   .  $bin...          - write grid and data in binary format
   .  $fgrid...        - if not initial time step grid is not written
   .  $s...             - write only specified subdomains. Positive id means that the
                                        subdomain is written.  Negative id means that the subdomain is
                                        not included.  If you use $s, default is that subdomains are not
                                        written unless specified.  However, if the first argument
                                        to $s is negative, all subdomains are selected for output (except
                                        the specified one).

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
  INT i,j,k,l,v;                                /* counters etc.							*/
  INT counter;                                  /* for formatting output					*/
  INT dat_pos;                      /* indicates offset in data file            */
  INT old_pos;
  char item[1024],it[256];              /* item buffers								*/
  char out_form[256];
  INT ic=0;                                             /* item length								*/
  VERTEX *vx;                                           /* a vertex pointer							*/
  ELEMENT *el;                                  /* an element pointer						*/
  NODE *no1, *no2;                  /* two node pointers                        */
  LINK *li;                         /* a link pointer                           */

  MULTIGRID *mg;                                /* our multigrid							*/
  HEAP *heap;

  char filename[NAMESIZE];              /* file name for header file				*/
  char filename_dat[NAMESIZE];      /* file name for data output file           */
  char filename_grid[NAMESIZE];      /* file name for grid output file          */
  char *c_ptr;
  PFILE *pf;                                            /* the output header file pointer               */
  PFILE *pf_txt;                    /* file pointer for ascii output            */
  PFILE_BIN *pf_bin;                /* file pointer for binary output           */

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
  INT numVerticesTot;                           /* total number of data points locally		*/
  INT numElements;                              /* number of elements locally				*/
  INT gnumVertices;                             /* number of data points globally			*/
  INT gnumVerticesTot;              /* total number of data points globally		*/
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

  INT notOnlyTetra;                 /* flag for tetrahedrons only grids         */
  INT notOnlyTriang;                /* flag for triangles only boundaries       */
  INT nibnd;                        /* number of inner boundary faces/lines     */
  INT gnibnd;                       /* global number of inner bndary faces/lines*/
  INT nobnd;                        /* number of outer boundary faces/lines     */
  INT gnobnd;                       /* global number of outer bndary faces/lines*/
  INT writeBnds=0;                  /* flag: write boundaries? (1=Inner, 2=All) */
  INT binaryOutput=0;               /* flag: write data in binary? (0=No|1=Yes) */
  INT writeGrid=1;                  /* flag: write grid? (0=No|1=Yes)           */
  INT buffer_INT[8];
  FLOAT buffer_FLOAT[8];
  INT usedBuf;
  INT *subdom=NULL;                             /* which subdomains to draw					*/
  INT subdomains=0;                             /* flag: draw only some subdomains			*/
  INT sdkey;
  char inex;
  time_t ltime;

  INT oe,ov;
  INT oibnd,oobnd;
  INT blocks;

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

    /* write data in binary format? */
    if (strncmp(argv[i],"bin",3)==0)
      binaryOutput=1;

    /* write grid or use existing data file as reference? */
    if (strncmp(argv[i],"fgrid",5)==0)
      writeGrid=0;

    if (strncmp(argv[i],"sd",1)==0) {
      if ( sscanf(argv[i],"sd %d",&j) )
      {
        if ( subdom==NULL )
        {
          heap = mg->theHeap;
          MarkTmpMem(heap, &sdkey);
          subdom = (INT *)GetTmpMem(heap, (mg->theBVPD.nSubDomains+1)*sizeof(INT), sdkey);
          if ( subdom==NULL ) {
            ReleaseTmpMem(heap, sdkey);
            PrintErrorMessage('E',"dataexplorer","could not allocate memory");
            return(PARAMERRORCODE);
          }
          /* default is that all subdomains are used */
          for (k=1; k<=mg->theBVPD.nSubDomains; k++) subdom[k]=(j>0) ? 0 : 1;
        }

        k = (j>0) ? j : -j;
        if ( k <= mg->theBVPD.nSubDomains && k > 0) subdom[k] = (j>0) ? 1 : 0;
        else UserWriteF("There is no subdomain %d\n",k);
        subdomains = 1;
      }
    }
  }

  if (ns==0 && nv==0 && ns_cell==0 && nv_cell==0)
    UserWrite("dataexplorer: no variables given, printing mesh data only\n");

  /* are inner and/or outer boundaries to be written? */
  if (ReadArgvINT("b",&writeBnds,argc,argv))
    writeBnds = 0;
  if (writeBnds>2) writeBnds=2;
  if (writeBnds<0) writeBnds=0;

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

  /* on first time step always write grid */
  if (strstr(filename,"0000")!=NULL)
    writeGrid=1;

  if (binaryOutput) {
    strcpy(filename_dat,filename);
    c_ptr=strrchr(filename_dat,'.');
    /*                  if (c_ptr!=NULL) */
    /*                          memset(c_ptr, '\0', 1); */
    strcat(filename_dat, ".bin");
    pf_bin = pfile_open_bin(filename_dat);
    if (pf_bin==NULL) {
      PrintErrorMessage('E',"dataexplorer","could not open data file");
      return(PARAMERRORCODE);
    }
    strcpy(out_form,"binary");
  } else {
    strcpy(filename_dat,filename);
    c_ptr=strrchr(filename_dat,'.');
    /*                  if (c_ptr!=NULL) */
    /*                          memset(c_ptr, '\0', 1); */
    strcat(filename_dat, ".dat");
    pf_txt = pfile_open(filename_dat);
    if (pf_txt==NULL) {
      PrintErrorMessage('E',"dataexplorer","could not open data file");
      return(PARAMERRORCODE);
    }
    strcpy(out_form,"text");
  }
  if (writeGrid)
    strcpy(filename_grid, filename_dat);
  else {
    strcpy(filename_grid, filename);
    c_ptr=strrchr(filename_grid,'.');
    memset(c_ptr, '\0', 1);
    if (binaryOutput)
      strcat(filename_grid, ".0000.bin");
    else
      strcat(filename_grid, ".0000.dat");
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
  numVertices = 0; numVerticesTot=0;
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el)) {
      if (!EstimateHere(el)) continue;
      if ( subdomains && !subdom[SUBDOMAIN(el)] ) continue;
      for (i=0; i<CORNERS_OF_ELEM(el); i++) {
        vx = MYVERTEX(CORNER(el,i));
        if (USED(vx)) continue;
        SETUSED(vx,1);
#ifdef ModelP
        if (DDD_InfoPriority(PARHDR(CORNER(el,i)))==PrioMaster)                         /* count only Master vertices */
#endif
        numVertices++;
        numVerticesTot++;
      }
    }

  if (subdomains)
  {
    /*
     * Total #vertices neccessary for correct data allocation.
     * If only some subdomains are drawn, we still need the total
     * number of vertices.
     */
    for (k=0; k<=TOPLEVEL(mg); k++)
      for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
        SETUSED(vx,0);
    numVerticesTot = 0;
    for (k=0; k<=TOPLEVEL(mg); k++)
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el)) {
        if (!EstimateHere(el)) continue;
        for (i=0; i<CORNERS_OF_ELEM(el); i++) {
          vx = MYVERTEX(CORNER(el,i));
          if (USED(vx)) continue;
          SETUSED(vx,1);
          numVerticesTot++;
        }
      }
  }

  /* count surface elements */
  numElements = 0;
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el)) {
      if (!EstimateHere(el)) continue;
      if ( subdomains && !subdom[SUBDOMAIN(el)] ) continue;
      numElements++;
    }

#ifdef ModelP
  gnumVertices = UG_GlobalSumINT(numVertices);
  gnumVerticesTot = UG_GlobalSumINT(numVerticesTot);
  gnumElements = UG_GlobalSumINT(numElements);
  /*    LocallyUniqueIDs(mg); */
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
  dat_pos = old_pos = 0;
  sprintf(it,"object 1 class array type float rank 1 shape %d items %d %s\ndata file %s,%d\n",
          DIM, gnumVertices, out_form, filename_grid, dat_pos);
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;
  pfile_sync(pf);
  if (binaryOutput)
    dat_pos+=DIM*gnumVertices*sizeof(FLOAT);

  /* unmark vertices */
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
      SETUSED(vx,0);

  /* write vertex coordinates */
  {
    INT count_rest = 0;

    counter=ov;
    for (k=0; k<=TOPLEVEL(mg); k++) {
      if (binaryOutput && !writeGrid) break;
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el)) {
        if (!EstimateHere(el)) continue;
        if ( subdomains && !subdom[SUBDOMAIN(el)] ) continue;

        for (i=0; i<CORNERS_OF_ELEM(el); i++)
        {
          vx = MYVERTEX(CORNER(el,i));
          if (USED(vx)) continue;
          SETUSED(vx,1);

#ifdef ModelP
          if (DDD_InfoPriority(PARHDR(CORNER(el,i)))==PrioMaster) {
#endif
          ID(vx) = counter;
          /* write the thing */
          if (binaryOutput) {
            buffer_FLOAT[0]=clampf(XC(vx));
            buffer_FLOAT[1]=clampf(YC(vx));
            buffer_FLOAT[2]=clampf(ZC(vx));
          } else {
#ifdef __TWODIM__
            sprintf(it,"\t%g\t%g\n", clampf(XC(vx)), clampf(YC(vx)));
#else
            sprintf(it,"\t%g\t%g\t%g\n", clampf(XC(vx)), clampf(YC(vx)),
                    clampf(ZC(vx)));
#endif
          }
          if (binaryOutput)
            pfile_tagged_write_FLOAT(pf_bin, buffer_FLOAT, DIM, counter);
          else {
            if (writeGrid)
              pfile_tagged_puts(pf_txt,it,counter);
            old_pos+=strlen(it);
          }
          counter++;
#ifdef ModelP
        }
        else {
          ID(vx) = ov+numVertices+count_rest;
          count_rest++;
        }
#endif
        }
      }
    }
  }
  if (binaryOutput)
    pfile_sync_bin(pf_bin);
  else {
    pfile_sync(pf_txt);
#ifdef ModelP
    dat_pos = UG_GlobalSumINT(old_pos);
#else
    dat_pos = old_pos;
#endif
    old_pos = 0;
  }

#ifdef ModelP
  {
    /* communicate local IDs */
    for (k=0; k<=TOPLEVEL(mg); k++) {
      DDD_IFAOneway(BorderVectorIF,GRID_ATTR(GRID_ON_LEVEL(mg,k)),IF_BACKWARD,sizeof(INT),Gather_VertexID, Scatter_VertexID);
    }
  }
#endif

  /* check if grid consists of tetrahedrons or triangles (in 2D) only */
  notOnlyTetra=0;
  for (k=0; k<=TOPLEVEL(mg); k++) {
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el)) {
      if (!EstimateHere(el)) continue;
      if ( subdomains && !subdom[SUBDOMAIN(el)] ) continue;
#ifdef __TWODIM__
      if (CORNERS_OF_ELEM(el)!=3) {
        notOnlyTetra = 1;
        break;
      }
#else
      if (CORNERS_OF_ELEM(el)!=4) {
        notOnlyTetra = 1;
        break;
      }
#endif
    }
    if (notOnlyTetra) break;
  }

#ifdef ModelP
  k = UG_GlobalSumINT(notOnlyTetra);
  notOnlyTetra = k;
#endif

  /****************************************************************/
  /*	2. write connections for domain								*/
  /****************************************************************/

  sprintf(it,"\n#\n# connections\n#\n");
  strcpy(item+ic,it); ic+=strlen(it);
  if (notOnlyTetra) {
    sprintf(it,"object 2 class array type int rank 1 shape %d items %d %s\ndata file %s,%d\n",
            4*(DIM-1), gnumElements, out_form, filename_grid, dat_pos);
    if (binaryOutput)
      dat_pos+=4*(DIM-1)*gnumElements*sizeof(INT);
  } else {
    sprintf(it,"object 2 class array type int rank 1 shape %d items %d %s\ndata file %s,%d\n",
            DIM+1, gnumElements, out_form, filename_grid, dat_pos);
    if (binaryOutput)
      dat_pos+=(DIM+1)*gnumElements*sizeof(INT);
  }
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;

  pfile_sync(pf);

  counter = 0;
  for (k=0; k<=TOPLEVEL(mg); k++) {
    if (binaryOutput && !writeGrid) break;
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (!EstimateHere(el)) continue;
      if ( subdomains && !subdom[SUBDOMAIN(el)] ) continue;

      switch(CORNERS_OF_ELEM(el))
      {
      case 3 :
        if (notOnlyTetra) {                               /* in 2D and NOT only triangles */
          if (binaryOutput) {
            buffer_INT[0]=ID(MYVERTEX(CORNER(el,0)));
            buffer_INT[1]=ID(MYVERTEX(CORNER(el,1)));
            buffer_INT[2]=ID(MYVERTEX(CORNER(el,2)));
            buffer_INT[3]=ID(MYVERTEX(CORNER(el,2)));
            usedBuf=4;
          } else
            sprintf(it,"\t%d\t%d\t%d\t%d\n",
                    ID(MYVERTEX(CORNER(el,0))),
                    ID(MYVERTEX(CORNER(el,1))),
                    ID(MYVERTEX(CORNER(el,2))),
                    ID(MYVERTEX(CORNER(el,2))));
        } else {                                          /* in 2D and only triangles */
          if (binaryOutput) {
            buffer_INT[0]=ID(MYVERTEX(CORNER(el,0)));
            buffer_INT[1]=ID(MYVERTEX(CORNER(el,1)));
            buffer_INT[2]=ID(MYVERTEX(CORNER(el,2)));
            usedBuf=3;
          } else
            sprintf(it,"\t%d\t%d\t%d\n",
                    ID(MYVERTEX(CORNER(el,0))),
                    ID(MYVERTEX(CORNER(el,1))),
                    ID(MYVERTEX(CORNER(el,2))));
        }
        break;

      case 4 :
#ifdef __TWODIM__
        if (binaryOutput) {
          buffer_INT[0]=ID(MYVERTEX(CORNER(el,0)));
          buffer_INT[1]=ID(MYVERTEX(CORNER(el,1)));
          buffer_INT[2]=ID(MYVERTEX(CORNER(el,3)));
          buffer_INT[3]=ID(MYVERTEX(CORNER(el,2)));
          usedBuf=4;
        } else
          sprintf(it,"\t%d\t%d\t%d\t%d\n",
                  ID(MYVERTEX(CORNER(el,0))),
                  ID(MYVERTEX(CORNER(el,1))),
                  ID(MYVERTEX(CORNER(el,3))),
                  ID(MYVERTEX(CORNER(el,2))));
#else
        if (notOnlyTetra) {                               /* in 3D and NOT only Tetrahedrons */
          if (binaryOutput) {
            buffer_INT[0]=ID(MYVERTEX(CORNER(el,0)));
            buffer_INT[1]=ID(MYVERTEX(CORNER(el,1)));
            buffer_INT[2]=ID(MYVERTEX(CORNER(el,2)));
            buffer_INT[3]=ID(MYVERTEX(CORNER(el,2)));
            buffer_INT[4]=ID(MYVERTEX(CORNER(el,3)));
            buffer_INT[5]=ID(MYVERTEX(CORNER(el,3)));
            buffer_INT[6]=ID(MYVERTEX(CORNER(el,3)));
            buffer_INT[7]=ID(MYVERTEX(CORNER(el,3)));
            usedBuf=8;
          } else
            sprintf(it,"\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                    ID(MYVERTEX(CORNER(el,0))),
                    ID(MYVERTEX(CORNER(el,1))),
                    ID(MYVERTEX(CORNER(el,2))),
                    ID(MYVERTEX(CORNER(el,2))),
                    ID(MYVERTEX(CORNER(el,3))),
                    ID(MYVERTEX(CORNER(el,3))),
                    ID(MYVERTEX(CORNER(el,3))),
                    ID(MYVERTEX(CORNER(el,3))));
        } else {                                         /* in 3D and only Tetrahedrons in grid */
          /*                                    if ((ID(MYVERTEX(CORNER(el,0))) > numVerticesTot+1) || */
          /*                                            (ID(MYVERTEX(CORNER(el,1))) > numVerticesTot+1) || */
          /*                                            (ID(MYVERTEX(CORNER(el,2))) > numVerticesTot+1) || */
          /*                                        (ID(MYVERTEX(CORNER(el,3))) > numVerticesTot+1)) { */
          /*                                            assert(0); */
          /*                                    } */

          if (binaryOutput) {
            buffer_INT[0]=ID(MYVERTEX(CORNER(el,0)));
            buffer_INT[1]=ID(MYVERTEX(CORNER(el,1)));
            buffer_INT[2]=ID(MYVERTEX(CORNER(el,2)));
            buffer_INT[3]=ID(MYVERTEX(CORNER(el,3)));
            usedBuf=4;
          } else
            sprintf(it,"\t%d\t%d\t%d\t%d\n",
                    ID(MYVERTEX(CORNER(el,0))),
                    ID(MYVERTEX(CORNER(el,1))),
                    ID(MYVERTEX(CORNER(el,2))),
                    ID(MYVERTEX(CORNER(el,3))));
        }
#endif
        break;

      case 5 :
        if (binaryOutput) {
          buffer_INT[0]=ID(MYVERTEX(CORNER(el,0)));
          buffer_INT[1]=ID(MYVERTEX(CORNER(el,1)));
          buffer_INT[2]=ID(MYVERTEX(CORNER(el,3)));
          buffer_INT[3]=ID(MYVERTEX(CORNER(el,2)));
          buffer_INT[4]=ID(MYVERTEX(CORNER(el,4)));
          buffer_INT[5]=ID(MYVERTEX(CORNER(el,4)));
          buffer_INT[6]=ID(MYVERTEX(CORNER(el,4)));
          buffer_INT[7]=ID(MYVERTEX(CORNER(el,4)));
          usedBuf=8;
        } else
          sprintf(it,"\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                  ID(MYVERTEX(CORNER(el,0))),
                  ID(MYVERTEX(CORNER(el,1))),
                  ID(MYVERTEX(CORNER(el,3))),
                  ID(MYVERTEX(CORNER(el,2))),
                  ID(MYVERTEX(CORNER(el,4))),
                  ID(MYVERTEX(CORNER(el,4))),
                  ID(MYVERTEX(CORNER(el,4))),
                  ID(MYVERTEX(CORNER(el,4))));
        break;

      case 6 :
        if (binaryOutput) {
          buffer_INT[0]=ID(MYVERTEX(CORNER(el,0)));
          buffer_INT[1]=ID(MYVERTEX(CORNER(el,1)));
          buffer_INT[2]=ID(MYVERTEX(CORNER(el,2)));
          buffer_INT[3]=ID(MYVERTEX(CORNER(el,2)));
          buffer_INT[4]=ID(MYVERTEX(CORNER(el,3)));
          buffer_INT[5]=ID(MYVERTEX(CORNER(el,4)));
          buffer_INT[6]=ID(MYVERTEX(CORNER(el,5)));
          buffer_INT[7]=ID(MYVERTEX(CORNER(el,5)));
          usedBuf=8;
        } else
          sprintf(it,"\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                  ID(MYVERTEX(CORNER(el,0))),
                  ID(MYVERTEX(CORNER(el,1))),
                  ID(MYVERTEX(CORNER(el,2))),
                  ID(MYVERTEX(CORNER(el,2))),
                  ID(MYVERTEX(CORNER(el,3))),
                  ID(MYVERTEX(CORNER(el,4))),
                  ID(MYVERTEX(CORNER(el,5))),
                  ID(MYVERTEX(CORNER(el,5))));
        break;

      case 8 :
        if (binaryOutput) {
          buffer_INT[0]=ID(MYVERTEX(CORNER(el,0)));
          buffer_INT[1]=ID(MYVERTEX(CORNER(el,1)));
          buffer_INT[2]=ID(MYVERTEX(CORNER(el,3)));
          buffer_INT[3]=ID(MYVERTEX(CORNER(el,2)));
          buffer_INT[4]=ID(MYVERTEX(CORNER(el,4)));
          buffer_INT[5]=ID(MYVERTEX(CORNER(el,5)));
          buffer_INT[6]=ID(MYVERTEX(CORNER(el,7)));
          buffer_INT[7]=ID(MYVERTEX(CORNER(el,6)));
          usedBuf=8;
        } else
          sprintf(it,"\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
                  ID(MYVERTEX(CORNER(el,0))),
                  ID(MYVERTEX(CORNER(el,1))),
                  ID(MYVERTEX(CORNER(el,3))),
                  ID(MYVERTEX(CORNER(el,2))),
                  ID(MYVERTEX(CORNER(el,4))),
                  ID(MYVERTEX(CORNER(el,5))),
                  ID(MYVERTEX(CORNER(el,7))),
                  ID(MYVERTEX(CORNER(el,6))));
        break;
      }
      if (binaryOutput)
        pfile_tagged_write_INT(pf_bin, buffer_INT, usedBuf, counter+oe);
      else {
        if (writeGrid)
          pfile_tagged_puts(pf_txt,it,counter+oe);
        old_pos+=strlen(it);
      }
      counter++;
    }
  }

  if (binaryOutput)
    pfile_sync_bin(pf_bin);
  else {
    pfile_sync(pf_txt);
#ifdef ModelP
    dat_pos += UG_GlobalSumINT(old_pos);
#else
    dat_pos += old_pos;
#endif
    old_pos = 0;
  }

#ifdef __TWODIM__
  if (notOnlyTetra)
    sprintf(it,"attribute \"element type\" string \"quads\"\n");
  else
    sprintf(it,"attribute \"element type\" string \"triangles\"\n");
#else
  if (notOnlyTetra)
    sprintf(it,"attribute \"element type\" string \"cubes\"\n");
  else
    sprintf(it,"attribute \"element type\" string \"tetrahedra\"\n");
#endif
  strcpy(item+ic,it); ic+=strlen(it);
  sprintf(it,"attribute \"ref\" string \"positions\"\n\n");
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;

  /**************************************************/
  /* 2.a write connections for boundaries           */
  /**************************************************/

  if (writeBnds) {
    /* unmark all elements */
    for (k=0; k<=TOPLEVEL(mg); k++)
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
        SETUSED(el,0);

    /* count number of sides on inner/outer boundaries */
    nibnd=0;
    nobnd=0;
    notOnlyTriang=0;
    for (k=0; k<=TOPLEVEL(mg); k++)         {
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if ((!EstimateHere(el)) || (OBJT(el)!=BEOBJ)) continue;

        SETUSED(el,1);
        for (i=0; i<SIDES_OF_ELEM(el); i++) {
#ifndef ModelP
          if (NBELEM(el,i)!=NULL && USED(NBELEM(el,i))) continue;
#else
          if (NBELEM(el,i)!=NULL && !EGHOST(NBELEM(el,i)) && USED(NBELEM(el,i))) continue;
#endif
          if (SIDE_ON_BND(el,i)) {
            if (NBELEM(el,i)==NULL) nobnd++;
            else nibnd++;
#ifdef __THREEDIM
            if (CORNERS_OF_SIDE(el,i)!=3) notOnlyTriang=1;
#endif
          }
        }
      }
    }
#ifdef ModelP
    k = UG_GlobalSumINT(notOnlyTriang);
    notOnlyTetra = k;
    gnibnd = UG_GlobalSumINT(nibnd);
    oibnd = get_offset(nibnd);

    gnobnd = UG_GlobalSumINT(nobnd);
    oobnd = get_offset(nobnd);

#else
    gnibnd = nibnd;
    gnobnd = nobnd;
    oibnd = 0;
    oobnd = 0;
#endif

  }

  if (writeBnds) {
    /* unmark all elements yet again */
    for (k=0; k<=TOPLEVEL(mg); k++)
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
        SETUSED(el,0);

    /* now writing inner boundary faces or lines */
    sprintf(it,"\n#\n# connections for inner boundaries\n#\n");
    strcpy(item+ic,it); ic+=strlen(it);
#ifdef __TWODIM__
    sprintf(it,"object 3 class array type int rank 1 shape %d items %d %s\ndata file %s,%d\n",
            2, gnibnd, out_form, filename_grid, dat_pos);
    if (binaryOutput)
      dat_pos+=2*gnibnd*sizeof(INT);
#else
    if (notOnlyTriang) {
      sprintf(it,"object 3 class array type int rank 1 shape %d items %d %s\ndata file %s,%d\n",
              4, gnibnd, out_form, filename_grid, dat_pos);
      if (binaryOutput)
        dat_pos+=4*gnibnd*sizeof(INT);
    } else {
      sprintf(it,"object 3 class array type int rank 1 shape %d items %d %s\ndata file %s,%d\n",
              3, gnibnd, out_form, filename_grid, dat_pos);
      if (binaryOutput)
        dat_pos+=3*gnibnd*sizeof(INT);
    }
#endif
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;
    pfile_sync(pf);

    counter = 0;
    for (k=0; k<=TOPLEVEL(mg); k++) {
      if (binaryOutput && !writeGrid) break;
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if ((!EstimateHere(el)) || (OBJT(el)!=BEOBJ)) continue;

        SETUSED(el,1);
        for (i=0; i<SIDES_OF_ELEM(el); i++) {
#ifndef ModelP
          if (NBELEM(el,i)!=NULL && USED(NBELEM(el,i))) continue;
#else
          if (NBELEM(el,i)!=NULL && !EGHOST(NBELEM(el,i)) && USED(NBELEM(el,i))) continue;
#endif
          if (SIDE_ON_BND(el,i) && NBELEM(el,i)!=NULL) {
            switch(CORNERS_OF_SIDE(el,i))
            {
            case 2 :                                        /* in 2D sides are edges => only two points */
              if (binaryOutput) {
                buffer_INT[0]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0))));
                buffer_INT[1]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1))));
                usedBuf=2;
              } else
                sprintf(it,"\t%d\t%d\n",
                        ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0)))),
                        ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1)))));
              break;
            case 3 :                                        /* has to be 3D */
              if (notOnlyTriang) {                                              /* not only faces with 3 corners */
                if (binaryOutput) {
                  buffer_INT[0]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0))));
                  buffer_INT[1]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1))));
                  buffer_INT[2]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2))));
                  buffer_INT[3]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2))));
                  usedBuf=4;
                } else
                  sprintf(it,"\t%d\t%d\t%d\t%d\n",
                          ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0)))),
                          ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1)))),
                          ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2)))),
                          ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2)))));
              } else {                                                          /* only triangles as faces */
                if (binaryOutput) {
                  buffer_INT[0]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0))));
                  buffer_INT[1]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1))));
                  buffer_INT[2]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2))));
                  usedBuf=3;
                } else
                  sprintf(it,"\t%d\t%d\t%d\n",
                          ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0)))),
                          ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1)))),
                          ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2)))));
              }
              break;
            case 4 :                                        /* has to be 3D and quads */
              if (binaryOutput) {
                buffer_INT[0]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0))));
                buffer_INT[1]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1))));
                buffer_INT[2]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,3))));
                buffer_INT[3]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2))));
                usedBuf=4;
              } else
                sprintf(it,"\t%d\t%d\t%d\t%d\n",
                        ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0)))),
                        ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1)))),
                        ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,3)))),
                        ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2)))));
              break;
            }
            if (binaryOutput)
              pfile_tagged_write_INT(pf_bin, buffer_INT, usedBuf, counter+oibnd);
            else {
              if (writeGrid)
                pfile_tagged_puts(pf_txt,it,counter+oibnd);
              old_pos+=strlen(it);
            }
            counter++;
          }
        }
      }
    }
    if(binaryOutput)
      pfile_sync_bin(pf_bin);
    else {
      pfile_sync(pf_txt);
#ifdef ModelP
      dat_pos += UG_GlobalSumINT(old_pos);
#else
      dat_pos += old_pos;
#endif
      old_pos = 0;
    }

#ifdef __TWODIM__
    sprintf(it,"attribute \"element type\" string \"lines\"\n");
#else
    if (notOnlyTriang)
      sprintf(it,"attribute \"element type\" string \"quads\"\n");
    else
      sprintf(it,"attribute \"element type\" string \"triangles\"\n");
#endif
    strcpy(item+ic,it); ic+=strlen(it);
    sprintf(it,"attribute \"ref\" string \"positions\"\n\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;
  }

  if (writeBnds==2) {
    /* unmark all elements again */
    for (k=0; k<=TOPLEVEL(mg); k++)
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
        SETUSED(el,0);

    /* now writing outer boundary faces or lines */
    sprintf(it,"\n#\n# connections for outer boundaries\n#\n");
    strcpy(item+ic,it); ic+=strlen(it);
#ifdef __TWODIM__
    sprintf(it,"object 4 class array type int rank 1 shape %d items %d %s\ndata file %s,%d\n",
            2, gnobnd, out_form, filename_grid, dat_pos);
    if (binaryOutput)
      dat_pos+=2*gnobnd*sizeof(INT);
#else
    if (notOnlyTriang) {
      sprintf(it,"object 4 class array type int rank 1 shape %d items %d %s\ndata file %s,%d\n",
              4, gnobnd, out_form, filename_grid, dat_pos);
      if (binaryOutput)
        dat_pos+=4*gnobnd*sizeof(INT);
    } else {
      sprintf(it,"object 4 class array type int rank 1 shape %d items %d %s\ndata file %s,%d\n",
              3, gnobnd, out_form, filename_grid, dat_pos);
      if (binaryOutput)
        dat_pos+=3*gnobnd*sizeof(INT);
    }
#endif
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;
    pfile_sync(pf);

    counter = 0;
    for (k=0; k<=TOPLEVEL(mg); k++) {
      if (binaryOutput && !writeGrid) break;
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if ((!EstimateHere(el)) || (OBJT(el)!=BEOBJ)) continue;

        SETUSED(el,1);
        for (i=0; i<SIDES_OF_ELEM(el); i++) {
          if ( NBELEM(el,i)!=NULL ) continue;
          if (SIDE_ON_BND(el,i) && NBELEM(el,i)==NULL) {
            switch(CORNERS_OF_SIDE(el,i))
            {
            case 2 :                                        /* in 2D sides are edges => only two points */
              if (binaryOutput) {
                buffer_INT[0]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0))));
                buffer_INT[1]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1))));
                usedBuf=2;
              } else
                sprintf(it,"\t%d\t%d\n",
                        ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0)))),
                        ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1)))));
              break;
            case 3 :                                        /* has to be 3D */
              if (notOnlyTriang) {                                              /* not only faces with 3 corners */
                if (binaryOutput) {
                  buffer_INT[0]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0))));
                  buffer_INT[1]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1))));
                  buffer_INT[2]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2))));
                  buffer_INT[3]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2))));
                  usedBuf=4;
                } else
                  sprintf(it,"\t%d\t%d\t%d\t%d\n",
                          ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0)))),
                          ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1)))),
                          ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2)))),
                          ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2)))));
              } else {                                                          /* only triangles as faces */
                if (binaryOutput) {
                  buffer_INT[0]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0))));
                  buffer_INT[1]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1))));
                  buffer_INT[2]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2))));
                  usedBuf=3;
                } else
                  sprintf(it,"\t%d\t%d\t%d\n",
                          ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0)))),
                          ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1)))),
                          ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2)))));
              }
              break;
            case 4 :                                        /* has to be 3D and quads */
              if (binaryOutput) {
                buffer_INT[0]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0))));
                buffer_INT[1]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1))));
                buffer_INT[2]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,3))));
                buffer_INT[3]=ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2))));
                usedBuf=4;
              } else
                sprintf(it,"\t%d\t%d\t%d\t%d\n",
                        ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,0)))),
                        ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,1)))),
                        ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,3)))),
                        ID(MYVERTEX(CORNER(el, CORNER_OF_SIDE(el,i,2)))));
              break;
            }
            if (binaryOutput)
              pfile_tagged_write_INT(pf_bin, buffer_INT, usedBuf, counter+oobnd);
            else {
              if (writeGrid)
                pfile_tagged_puts(pf_txt,it,counter+oobnd);
              old_pos+=strlen(it);
            }
            counter++;
          }
        }
      }
    }
    if (binaryOutput)
      pfile_sync_bin(pf_bin);
    else {
      pfile_sync(pf_txt);
#ifdef ModelP
      dat_pos += UG_GlobalSumINT(old_pos);
#else
      dat_pos += old_pos;
#endif
      old_pos = 0;
    }

#ifdef __TWODIM__
    sprintf(it,"attribute \"element type\" string \"lines\"\n");
#else
    if (notOnlyTriang)
      sprintf(it,"attribute \"element type\" string \"quads\"\n");
    else
      sprintf(it,"attribute \"element type\" string \"triangles\"\n");
#endif
    strcpy(item+ic,it); ic+=strlen(it);
    sprintf(it,"attribute \"ref\" string \"positions\"\n\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;
  }

  /* reset data file position counter if grid has not been written */
  if (!writeGrid)
    dat_pos=0;

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

    sprintf(it,"object %d class array type float rank 0 items %d %s\ndata file %s,%d\n",
            blocks+2+writeBnds, gnumVertices, out_form, filename_dat, dat_pos);
    if (binaryOutput)
      dat_pos+=gnumVertices*sizeof(FLOAT);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;
    pfile_sync(pf);

    for (k=0; k<=TOPLEVEL(mg); k++)
      for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
        SETUSED(vx,0);

    counter=0;
    for (k=0; k<=TOPLEVEL(mg); k++) {
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if (!EstimateHere(el)) continue;
        if ( subdomains && !subdom[SUBDOMAIN(el)] ) continue;

        for (i=0; i<CORNERS_OF_ELEM(el); i++)
          CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
        {
#ifdef ModelP
          if (DDD_InfoPriority(PARHDR(CORNER(el,i)))!=PrioMaster) continue;
#endif
          vx = MYVERTEX(CORNER(el,i));
          if (USED(vx)) continue;
          SETUSED(vx,1);

          /* get local coordinate of corner */
          LocalCornerCoordinates(DIM,TAG(el),i,local);
          for (j=0; j<DIM; j++) LocalCoord[j] = local[j];

          /* scalar components */
          eval_s = es[v]->EvalProc;
          value = eval_s(el,(const DOUBLE **)CornersCoord,LocalCoord);
          if (binaryOutput) {
            buffer_FLOAT[0]=clampf(value);
            pfile_tagged_write_FLOAT(pf_bin, buffer_FLOAT, 1, counter+ov);
          } else {
            sprintf(it,"\t%g\n",clampf(value));
            pfile_tagged_puts(pf_txt,it,counter+ov);
            old_pos+=strlen(it);
          }
          counter++;
        }
      }
    }
    if (binaryOutput)
      pfile_sync_bin(pf_bin);
    else {
      pfile_sync(pf_txt);
#ifdef ModelP
      dat_pos += UG_GlobalSumINT(old_pos);
#else
      dat_pos += old_pos;
#endif
      old_pos = 0;
    }

    /* domain data fields */
    sprintf(it,"\nobject \"data%d\" class field\n", blocks);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"positions\" value 1\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"connections\" value 2\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"data\" value %d\n\n", blocks+2+writeBnds);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    if (writeBnds) {
      /* inner boundaries data fields */
      sprintf(it,"\nobject \"ibnd%d\" class field\n", blocks);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"positions\" value 1\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"connections\" value 3\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"data\" value %d\n", blocks+2+writeBnds);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;
    }

    if (writeBnds==2) {
      /* outer boundaries data fields */
      sprintf(it,"\nobject \"obnd%d\" class field\n", blocks);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"positions\" value 1\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"connections\" value 4\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"data\" value %d\n", blocks+2+writeBnds);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;
    }

    blocks++;
  }

  /* write all vector node data */
  for (v=0; v<nv; v++)
  {
    pre     = ev[v]->PreprocessProc;
    if (pre!=NULL) pre(ev_name[v],mg);

    sprintf(it,"#\n# data block %d\n#\n", blocks);
    strcpy(item+ic,it); ic+=strlen(it);

    sprintf(it,"object %d class array type float rank 1 shape %d items %d %s\ndata file %s,%d\n",
            blocks+2+writeBnds, DIM, gnumVertices, out_form, filename_dat, dat_pos);
    if (binaryOutput)
      dat_pos+=DIM*gnumVertices*sizeof(FLOAT);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;
    pfile_sync(pf);

    for (k=0; k<=TOPLEVEL(mg); k++)
      for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
        SETUSED(vx,0);
    counter=0;
    for (k=0; k<=TOPLEVEL(mg); k++)         {
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if (!EstimateHere(el)) continue;
        if ( subdomains && !subdom[SUBDOMAIN(el)] ) continue;

        for (i=0; i<CORNERS_OF_ELEM(el); i++)
          CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
        {
#ifdef ModelP
          if (DDD_InfoPriority(PARHDR(CORNER(el,i)))!=PrioMaster) continue;
#endif
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
          if (binaryOutput) {
            buffer_FLOAT[0]=clampf(vval[0]);
            buffer_FLOAT[1]=clampf(vval[1]);
          } else
            sprintf(it,"\t%g\t%g\n",clampf(vval[0]),clampf(vval[1]));
#else
          if (binaryOutput) {
            buffer_FLOAT[0]=clampf(vval[0]);
            buffer_FLOAT[1]=clampf(vval[1]);
            buffer_FLOAT[2]=clampf(vval[2]);
          } else
            sprintf(it,"\t%g\t%g\t%g\n",clampf(vval[0]),clampf(vval[1]),
                    clampf(vval[2]));
#endif
          if (binaryOutput)
            pfile_tagged_write_FLOAT(pf_bin, buffer_FLOAT, DIM, counter+ov);
          else {
            pfile_tagged_puts(pf_txt,it,counter+ov);
            old_pos+=strlen(it);
          }
          counter++;
        }
      }
    }
    if (binaryOutput)
      pfile_sync_bin(pf_bin);
    else {
      pfile_sync(pf_txt);
#ifdef ModelP
      dat_pos += UG_GlobalSumINT(old_pos);
#else
      dat_pos += old_pos;
#endif
      old_pos = 0;
    }

    /* domain data fields */
    sprintf(it,"\nobject \"data%d\" class field\n", blocks);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"positions\" value 1\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"connections\" value 2\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"data\" value %d\n", blocks+2+writeBnds);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    if (writeBnds) {
      /* inner boundaries data fields */
      sprintf(it,"\nobject \"ibnd%d\" class field\n", blocks);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"positions\" value 1\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"connections\" value 3\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"data\" value %d\n", blocks+2+writeBnds);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;
    }

    if (writeBnds==2) {
      /* outer boundaries data fields */
      sprintf(it,"\nobject \"obnd%d\" class field\n", blocks);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"positions\" value 1\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"connections\" value 4\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"data\" value %d\n", blocks+2+writeBnds);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;
    }

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

    sprintf(it,"object %d class array type float rank 0 items %d %s\ndata file %s,%d\n",
            blocks+2+writeBnds, gnumElements, out_form, filename_dat, dat_pos);
    if (binaryOutput)
      dat_pos+=gnumElements*sizeof(FLOAT);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;
    pfile_sync(pf);

    counter=0;
    for (k=0; k<=TOPLEVEL(mg); k++) {
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if (!EstimateHere(el)) continue;
        if ( subdomains && !subdom[SUBDOMAIN(el)] ) continue;

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
        if (binaryOutput) {
          buffer_FLOAT[0]=clampf(value);
          pfile_tagged_write_FLOAT(pf_bin, buffer_FLOAT, 1, counter+oe);
        } else {
          sprintf(it,"\t%g\n",clampf(value));
          pfile_tagged_puts(pf_txt,it,counter+oe);
          old_pos+=strlen(it);
        }
        counter++;
      }
    }
    if (binaryOutput)
      pfile_sync_bin(pf_bin);
    else {
      pfile_sync(pf_txt);
#ifdef ModelP
      dat_pos += UG_GlobalSumINT(old_pos);
#else
      dat_pos += old_pos;
#endif
      old_pos = 0;
    }

    sprintf(it,"attribute \"dep\" string \"connections\"\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    /* domain data fields */
    sprintf(it,"\nobject \"data%d\" class field\n", blocks);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"positions\" value 1\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"connections\" value 2\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"data\" value %d\n", blocks+2+writeBnds);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    if (writeBnds) {
      /* inner boundaries data fields */
      sprintf(it,"\nobject \"ibnd%d\" class field\n", blocks);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"positions\" value 1\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"connections\" value 3\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"data\" value %d\n", blocks+2+writeBnds);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;
    }

    if (writeBnds==2) {
      /* outer boundaries data fields */
      sprintf(it,"\nobject \"obnd%d\" class field\n", blocks);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"positions\" value 1\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"connections\" value 4\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"data\" value %d\n", blocks+2+writeBnds);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;
    }

    blocks++;
  }

  /* write all element vector data */
  for (v=0; v<nv_cell; v++)
  {
    pre     = ev_cell[v]->PreprocessProc;
    if (pre!=NULL) pre(ev_cell_name[v],mg);

    sprintf(it,"#\n# data block %d\n#\n", blocks);
    strcpy(item+ic,it); ic+=strlen(it);

    sprintf(it,"object %d class array type float rank 1 shape %d items %d %s\ndata file %s,%d\n",
            blocks+2+writeBnds, DIM, gnumElements, out_form, filename_dat, dat_pos);
    if (binaryOutput)
      dat_pos+=DIM*gnumElements*sizeof(FLOAT);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;
    pfile_sync(pf);

    counter=0;
    for (k=0; k<=TOPLEVEL(mg); k++) {
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if (!EstimateHere(el)) continue;
        if ( subdomains && !subdom[SUBDOMAIN(el)] ) continue;

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
        if (binaryOutput) {
          buffer_FLOAT[0]=clampf(vval[0]);
          buffer_FLOAT[1]=clampf(vval[1]);
        } else
          sprintf(it,"\t%g\t%g\n",clampf(vval[0]),clampf(vval[1]));
#else
        if (binaryOutput) {
          buffer_FLOAT[0]=clampf(vval[0]);
          buffer_FLOAT[1]=clampf(vval[1]);
          buffer_FLOAT[2]=clampf(vval[2]);
        } else
          sprintf(it,"\t%g\t%g\t%g\n",clampf(vval[0]),clampf(vval[1]),
                  clampf(vval[2]));
#endif
        if (binaryOutput)
          pfile_tagged_write_FLOAT(pf_bin, buffer_FLOAT, DIM, counter+ov);
        else {
          pfile_tagged_puts(pf_txt,it,counter+ov);
          old_pos+=strlen(it);
        }
        counter++;
      }
    }
    if (binaryOutput)
      pfile_sync_bin(pf_bin);
    else {
      pfile_sync(pf_txt);
#ifdef ModelP
      dat_pos += UG_GlobalSumINT(old_pos);
#else
      dat_pos += old_pos;
#endif
      old_pos = 0;
    }

    sprintf(it,"attribute \"dep\" string \"connections\"\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    /* domain data fields */
    sprintf(it,"\nobject \"data%d\" class field\n", blocks);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"positions\" value 1\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"connections\" value 2\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    sprintf(it,"component \"data\" value %d\n", blocks+2+writeBnds);
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    if (writeBnds) {
      /* inner boundaries data fields */
      sprintf(it,"\nobject \"ibnd%d\" class field\n", blocks);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"positions\" value 1\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"connections\" value 3\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"data\" value %d\n", blocks+2+writeBnds);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;
    }

    if (writeBnds==2) {
      /* outer boundaries data fields */
      sprintf(it,"\nobject \"obnd%d\" class field\n", blocks);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"positions\" value 1\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"connections\" value 4\n");
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;

      sprintf(it,"component \"data\" value %d\n", blocks+2+writeBnds);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;
    }

    blocks++;
  }

  /* put all data blocks into a group for select module */
  if (blocks > 1) {
    sprintf(it,"\nobject \"default\" class group\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;
    for (i = 1; i < blocks; i++) {
      sprintf(it,"member \"data%d\" value \"data%d\"\n", i, i);
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;
    }
  }

  pfile_close(pf);
  if (binaryOutput)
    pfile_close_bin(pf_bin);
  else
    pfile_close(pf_txt);

  if ( subdomains )
    ReleaseTmpMem(heap, sdkey);

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
