// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      avs.c		                                                    */
/*                                                                          */
/* Purpose:   avs output                                                                                        */
/*                                                                          */
/* Author:	  Peter Bastian                                                                                         */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: peter@ica3.uni-stuttgart.de							*/
/*			  fon: 0049-(0)711-685-7003										*/
/*			  fax: 0049-(0)711-685-7000										*/
/*																			*/
/* History:   07.11.96 begin, ug version 3.4, derived from tecplot.c		*/
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "devices.h"
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

#include "avs.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
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

/* data for CVS */
static char RCS_ID("$Header$",UG_RCS_STRING);


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
    theGrid = theMG->grids[k];

    /* vertices */
    for (theVertex=FIRSTVERTEX(theGrid); theVertex!=NULL; theVertex=SUCCV(theVertex))
      ID(theVertex) = nv++;

    /* nodes */
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
      ID(theNode) = nn++;

    /* elements */
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
    theGrid = theMG->grids[k];

    /* vertices */
    for (theVertex=FIRSTVERTEX(theGrid); theVertex!=NULL; theVertex=SUCCV(theVertex))
      ID(theVertex) += ov;

    /* nodes */
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
      ID(theNode) += on;

    /* elements */
    for (theElement=FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
      ID(theElement) += oe;
  }

  return(0);
}

static INT AVS_GlobalSumINT (INT i)
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

static INT AVSCommand (INT argc, char **argv)
{
  INT i,j,k,v;                                  /* counters etc.							*/
  INT counter;                                  /* for formatting output					*/
  char buf[128];                                /* for messages								*/
  char item[1024],it[256];              /* item buffers								*/
  INT ic=0;                                             /* item length								*/
  VERTEX *vx;                                           /* a vertex pointer							*/
  ELEMENT *el;                                  /* an element pointer						*/

  MULTIGRID *mg;                                /* our multigrid							*/
  char filename[NAMESIZE];      /* file name for output file				*/
  PFILE *pf;                                            /* the output file pointer					*/

  EVALUES *zcoord;                              /* use scalar as z coordinate in 2D only	*/
  char zcoord_name[NAMESIZE];           /* name for eval functions                  */

  INT ns;                                               /* number of scalar eval procs				*/
  INT nv;                                               /* number of vector eval procs				*/
  EVALUES *es[MAXVARIABLES];            /* pointers to scalar eval function desc	*/
  char es_name[MAXVARIABLES][NAMESIZE];         /* names for eval functions     */
  EVECTOR *ev[MAXVARIABLES];            /* pointers to vector eval function desc	*/
  char ev_name[MAXVARIABLES][NAMESIZE];         /* names for eval functions     */
  INT ns_cell;                                  /* number of scalar eval procs				*/
  INT nv_cell;                                  /* number of vector eval procs				*/
  EVALUES *es_cell[MAXVARIABLES];       /* pointers to scalar eval function desc*/
  char es_cell_name[MAXVARIABLES][NAMESIZE];            /* names for eval functions */
  EVECTOR *ev_cell[MAXVARIABLES];       /* pointers to vector eval function desc*/
  char ev_cell_name[MAXVARIABLES][NAMESIZE];            /* names for eval functions */
  char s[NAMESIZE];                             /* name of eval proc						*/
  INT numNodes;                                 /* number of data points locally			*/
  INT numElements;                              /* number of elements locally				*/
  INT gnumNodes;                                /* number of data points globally			*/
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
  DOUBLE scale;
  VECDATA_DESC *displacement;
  const SHORT *comp;
  time_t ltime;
  double scaling=1.0;
  INT oe,on;


  /* get current multigrid */
  mg = GetCurrentMultigrid();
  if (mg==NULL)
  {
    PrintErrorMessage('W',"avs","no multigrid open\n");
    return (OKCODE);
  }

  /* scan options */
  ns = nv = ns_cell = nv_cell = 0;
  zcoord = NULL;
  displacement = ReadArgvVecDesc(mg,"displacement",argc,argv);
  if (displacement != NULL) {
    if (VD_NCMPS_IN_TYPE(displacement,NODEVECTOR) < DIM)
      return(CMDERRORCODE);
    comp = VD_CMPPTR_OF_TYPE(displacement,NODEVECTOR,0);
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
        sprintf(buf,"could not find scalar eval proc %s\n",s);
        PrintErrorMessage('E',"avs:",buf);
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
        PrintErrorMessage('E',"avs:","too many scalar variables specified\n");
        break;
      }
      sscanf(argv[i],"ns %s", s);
      es[ns] = GetElementValueEvalProc(s);
      if (es[ns]==NULL)
      {
        sprintf(buf,"could not find scalar eval proc %s\n",s);
        PrintErrorMessage('E',"avs:",buf);
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
        PrintErrorMessage('E',"avs:","too many vector variables specified\n");
        break;
      }
      sscanf(argv[i],"nv %s", s);
      ev[nv] = GetElementVectorEvalProc(s);
      if (ev[nv]==NULL)
      {
        sprintf(buf,"could not find vector eval proc %s\n",s);
        PrintErrorMessage('E',"avs:",buf);
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
        PrintErrorMessage('E',"avs:","too many scalar variables specified\n");
        break;
      }
      sscanf(argv[i],"cs %s", s);
      es_cell[ns_cell] = GetElementValueEvalProc(s);
      if (es_cell[ns_cell]==NULL)
      {
        sprintf(buf,"could not find scalar eval proc %s\n",s);
        PrintErrorMessage('E',"avs:",buf);
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
        PrintErrorMessage('E',"avs:","too many vector variables specified\n");
        break;
      }
      sscanf(argv[i],"cv %s", s);
      ev_cell[nv_cell] = GetElementVectorEvalProc(s);
      if (ev_cell[nv_cell]==NULL)
      {
        sprintf(buf,"could not find vector eval proc %s\n",s);
        PrintErrorMessage('E',"avs:",buf);
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
  if (ns==0 && nv==0) UserWrite("avs: no variables given, printing mesh data only\n");

  /* get file name and open output file */
  if (sscanf(argv[0],expandfmt(CONCAT3(" avs %",NAMELENSTR,"[ -~]")),filename)!=1)
  {
    PrintErrorMessage('E',"avs","could not read name of output file");
    return(PARAMERRORCODE);
  }
  pf = pfile_open(filename);
  if (pf==NULL) return(PARAMERRORCODE);

  /********************************/
  /* TITLE                                              */
  /********************************/

  time(&ltime);
  sprintf(it,"# AVS data produced by UG\n");
  strcpy(item+ic,it); ic+=strlen(it);
  sprintf(it,"# multigrid: %s\n",mg->v.name);
  strcpy(item+ic,it); ic+=strlen(it);
  sprintf(it,"# date:      %s",ctime(&ltime));
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;

  /********************************/
  /* compute sizes				*/
  /********************************/

  /* clear USED flag in vertices on all levels */
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
      SETUSED(vx,0);

  /* run thru all levels of elements and set index */
  numNodes = numElements = 0;
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (NSONS(el)>0) continue;                        /* process finest level elements only */
      numElements++;                                            /* increase element counter */
      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        vx = MYVERTEX(CORNER(el,i));
        if (USED(vx)) continue;                         /* we have this one already */

        ++numNodes;                                                     /* number of data points, begins with 1 ! */
        SETUSED(vx,1);                                          /* tag vector as visited */
      }
    }

        #ifdef ModelP
  gnumNodes = AVS_GlobalSumINT(numNodes);
  gnumElements = AVS_GlobalSumINT(numElements);
  on=get_offset(numNodes);
  oe=get_offset(numElements);
  GloballyUniqueIDs(mg);       /* renumber objects */
        #else
  gnumNodes = numNodes;
  gnumElements = numElements;
  oe=on=0;
        #endif

  /********************************/
  /* (1) now we are ready to write*/
  /* the sizes to the file		*/
  /********************************/

  sprintf(it,"%d %d %d %d %d\n",
          gnumNodes,                                    /* total number of nodes to write       */
          gnumElements,                         /* total number of cells			*/
          ns+1+nv*DIM,                          /* number of doubles per node		*/
          ns_cell+nv_cell*DIM,          /* doubles per cell is always zero !*/
          0);                                                   /* doubles per model                            */
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;

  /********************************/
  /* (2) next are the coordinates */
  /* of nodes in point block form */
  /* This is pfile segment 0      */
  /********************************/

  if (DIM==2 && zcoord!=NULL)
  {
    pre    = zcoord->PreprocessProc;
    if (pre!=NULL) pre(zcoord_name,mg);
  }

  /* clear USED flag in vertices on all levels */
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
      SETUSED(vx,0);

  counter=0;
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (NSONS(el)>0) continue;                        /* process finest level elements only */
      for (i=0; i<CORNERS_OF_ELEM(el); i++)
        CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));         /* x,y,z of corners */
      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        vx = MYVERTEX(CORNER(el,i));
        if (USED(vx)) continue;                         /* we have this one already */
        SETUSED(vx,1);                                          /* tag vector as visited */

        /* get local coordinate of corner */
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
        sprintf(it,"%d %lg %lg %lg\n",ID(vx),x,y,z);
        pfile_tagged_puts(pf,it,counter+on);
        counter++;
      }
    }

  pfile_sync(pf);       /* end of segment */

  /********************************/
  /* (3) next is the cell               */
  /*     connectivity list		*/
  /********************************/

  counter=0;
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (NSONS(el)>0) continue;                        /* process finest level elements only */

      /* cell number and material id */
                        #ifdef ModelP
      sprintf(it,"%d %d ",ID(el),me);
                        #else
      sprintf(it,"%d %d ",ID(el),1);
                        #endif
      strcpy(item+ic,it); ic+=strlen(it);

      if (DIM==2)
      {
        switch (CORNERS_OF_ELEM(el))
        {
        case 3 :
          sprintf(it,"tri %d %d %d\n",
                  ID(MYVERTEX(CORNER(el,0))),
                  ID(MYVERTEX(CORNER(el,1))),
                  ID(MYVERTEX(CORNER(el,2))));
          break;
        case 4 :
          sprintf(it,"quad %d %d %d %d\n",
                  ID(MYVERTEX(CORNER(el,0))),
                  ID(MYVERTEX(CORNER(el,1))),
                  ID(MYVERTEX(CORNER(el,2))),
                  ID(MYVERTEX(CORNER(el,3))));
          break;
        }
      }
      else                   /* DIM==3 */
      {
        switch (CORNERS_OF_ELEM(el))
        {
        case 4 :
          sprintf(it,"tet %d %d %d %d\n",
                  ID(MYVERTEX(CORNER(el,0))),
                  ID(MYVERTEX(CORNER(el,2))),
                  ID(MYVERTEX(CORNER(el,1))),
                  ID(MYVERTEX(CORNER(el,3))));
          break;
        case 5 :
          sprintf(it,"pyr %d %d %d %d %d\n",
                  ID(MYVERTEX(CORNER(el,4))),
                  ID(MYVERTEX(CORNER(el,0))),
                  ID(MYVERTEX(CORNER(el,1))),
                  ID(MYVERTEX(CORNER(el,2))),
                  ID(MYVERTEX(CORNER(el,3))));
          break;
        case 6 :
          sprintf(it,"prism %d %d %d %d %d %d\n",
                  ID(MYVERTEX(CORNER(el,3))),
                  ID(MYVERTEX(CORNER(el,4))),
                  ID(MYVERTEX(CORNER(el,5))),
                  ID(MYVERTEX(CORNER(el,0))),
                  ID(MYVERTEX(CORNER(el,1))),
                  ID(MYVERTEX(CORNER(el,2))));
          break;
        case 8 :
          sprintf(it,"hex %d %d %d %d %d %d %d %d\n",
                  ID(MYVERTEX(CORNER(el,4))),
                  ID(MYVERTEX(CORNER(el,5))),
                  ID(MYVERTEX(CORNER(el,6))),
                  ID(MYVERTEX(CORNER(el,7))),
                  ID(MYVERTEX(CORNER(el,0))),
                  ID(MYVERTEX(CORNER(el,1))),
                  ID(MYVERTEX(CORNER(el,2))),
                  ID(MYVERTEX(CORNER(el,3))));
          break;
        }
      }
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_tagged_puts(pf,item,counter+oe); ic=0;
      counter++;
    }

  pfile_sync(pf);       /* end of segment */

  /********************************/
  /* node data					*/
  /********************************/

  if (1+ns+nv>0)
  {

    /********************************/
    /* (4) data associated with		*/
    /*     nodes					*/
    /********************************/

    sprintf(it,"%d ",ns+nv+1);                                          /* number of components */
    strcpy(item+ic,it); ic+=strlen(it);
    for (k=0; k<ns; k++) {
      sprintf(it,"1 ");                                                         /* scalar component		*/
      strcpy(item+ic,it); ic+=strlen(it);
    }
    sprintf(it,"1 ");                                                           /* scalar component me	*/
    strcpy(item+ic,it); ic+=strlen(it);
    for (k=0; k<nv; k++) {
      sprintf(it,"%d ",DIM);                                            /* vector component		*/
      strcpy(item+ic,it); ic+=strlen(it);
    }
    sprintf(it,"\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    /********************************/
    /* (5) component labels			*/
    /*                                                  */
    /********************************/

    /* now all scalar variables */
    for (i=0; i<ns; i++) {
      sprintf(it,"%s, \n",ENVITEM_NAME(es[i]));           /* no units */
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;
    }
    sprintf(it,"%s, \n","procid");             /* no units */
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;
    for (i=0; i<nv; i++) {
      sprintf(it,"%s, \n",ENVITEM_NAME(ev[i]));           /* no units */
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;
    }

    /********************************/
    /* (6) all the node data		*/
    /*                                                  */
    /********************************/

    /* CAUTION: AVS supports only a point block format. Therefore
       all prepare functions are called at the beginning.
       PREPARE FUNCTIONS MUST NOT INTERFERE WITH EACH OTHER !
     */

    /* execute all prepare functions */
    for (v=0; v<ns; v++)
    {
      pre    = es[v]->PreprocessProc;

      /* execute prepare function */
      if (pre!=NULL) pre(es_name[v],mg);
    }
    for (v=0; v<nv; v++)
    {
      pre    = ev[v]->PreprocessProc;

      /* execute prepare function */
      if (pre!=NULL) pre(ev_name[v],mg);
    }

    /* now the data in point block format */
    /* clear USED flag in vertices on all levels */
    for (k=0; k<=TOPLEVEL(mg); k++)
      for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
        SETUSED(vx,0);
    counter=0;
    for (k=0; k<=TOPLEVEL(mg); k++)
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if (NSONS(el)>0) continue;                              /* process finest level elements only */
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
          CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));                       /* x,y,z of corners */
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
        {
          vx = MYVERTEX(CORNER(el,i));
          if (USED(vx)) continue;                               /* we have this one already */
          SETUSED(vx,1);                                        /* tag vector as visited */

          /* get local coordinate of corner */
          LocalCornerCoordinates(DIM,TAG(el),i,local);
          for (j=0; j<DIM; j++) LocalCoord[j] = local[j];

          /* write node id */
          sprintf(it,"%d ",ID(MYVERTEX(CORNER(el,i))));
          strcpy(item+ic,it); ic+=strlen(it);

          /* scalar components */
          for (v=0; v<ns; v++)
          {
            eval_s = es[v]->EvalProc;
            value = eval_s(el,(const DOUBLE **)CornersCoord,LocalCoord);
            sprintf(it,"%lg ",value);
            strcpy(item+ic,it); ic+=strlen(it);
          }

                                        #ifdef ModelP
          value = (DOUBLE) me;
                                        #else
          value = 1.0;
                                        #endif
          sprintf(it,"%lg ",value);
          strcpy(item+ic,it); ic+=strlen(it);

          /* vector components */
          for (v=0; v<nv; v++)
          {
            eval_v = ev[v]->EvalProc;
            eval_v(el,(const DOUBLE **)CornersCoord,LocalCoord,vval);
            if (DIM==2)
              sprintf(it,"%lg %lg ",vval[0],vval[1]);
            else
              sprintf(it,"%lg %lg %lg ",vval[0],vval[1],vval[2]);
            strcpy(item+ic,it); ic+=strlen(it);
          }

          sprintf(it,"\n");
          strcpy(item+ic,it); ic+=strlen(it);
          pfile_tagged_puts(pf,item,counter+on); ic=0;
          counter++;
        }
      }
  }

  pfile_sync(pf);       /* end of segment */

  /********************************/
  /* cell data					*/
  /********************************/

  if (ns_cell+nv_cell>0)
  {

    /********************************/
    /* (4) data associated with		*/
    /*     elements					*/
    /********************************/

    sprintf(it,"%d ",ns_cell+nv_cell);                          /* number of components */
    strcpy(item+ic,it); ic+=strlen(it);
    for (k=0; k<ns_cell; k++) {
      sprintf(it,"1 ");                                                         /* scalar component		*/
      strcpy(item+ic,it); ic+=strlen(it);
    }
    for (k=0; k<nv_cell; k++) {
      sprintf(it,"%d ",DIM);                                            /* vector component		*/
      strcpy(item+ic,it); ic+=strlen(it);
    }
    sprintf(it,"\n");
    strcpy(item+ic,it); ic+=strlen(it);
    pfile_master_puts(pf,item); ic=0;

    /********************************/
    /* (5) component labels			*/
    /*                                                  */
    /********************************/

    /* now all scalar variables */
    for (i=0; i<ns_cell; i++) {
      sprintf(it,"%s, \n",ENVITEM_NAME(es_cell[i]));           /* no units */
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;
    }
    for (i=0; i<nv_cell; i++) {
      sprintf(it,"%s, \n",ENVITEM_NAME(ev_cell[i]));           /* no units */
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_master_puts(pf,item); ic=0;
    }

    /********************************/
    /* (6) all the cell data		*/
    /*                                                  */
    /********************************/

    /* CAUTION: AVS supports only a point block format. Therefore
       all prepare functions are called at the beginning.
       PREPARE FUNCTIONS MUST NOT INTERFERE WITH EACH OTHER !
     */

    /* execute all prepare functions */
    for (v=0; v<ns_cell; v++)
    {
      pre    = es_cell[v]->PreprocessProc;

      /* execute prepare function */
      if (pre!=NULL) pre(es_cell_name[v],mg);
    }
    for (v=0; v<nv_cell; v++)
    {
      pre    = ev_cell[v]->PreprocessProc;

      /* execute prepare function */
      if (pre!=NULL) pre(ev_cell_name[v],mg);
    }

    /* now the data in point block format */
    counter=0;
    for (k=0; k<=TOPLEVEL(mg); k++)
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if (NSONS(el)>0) continue;                              /* process finest level elements only */
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
          CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));                       /* x,y,z of corners */
        /* compute center in local coordinates */
        for (j=0; j<DIM; j++) LocalCoord[j] = 0.0;
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
        {
          /* get local coordinate of corner */
          LocalCornerCoordinates(DIM,TAG(el),i,local);
          for (j=0; j<DIM; j++) LocalCoord[j] += local[j];
        }
        for (j=0; j<DIM; j++) LocalCoord[j] /= ((DOUBLE)CORNERS_OF_ELEM(el));

        /* write cell id */
        sprintf(it,"%d ",ID(el));
        strcpy(item+ic,it); ic+=strlen(it);

        /* scalar components */
        for (v=0; v<ns_cell; v++)
        {
          eval_s = es_cell[v]->EvalProc;
          value = eval_s(el,(const DOUBLE **)CornersCoord,LocalCoord);
          sprintf(it,"%lg ",value);
          strcpy(item+ic,it); ic+=strlen(it);
        }

        /* vector components */
        for (v=0; v<nv_cell; v++)
        {
          eval_v = ev_cell[v]->EvalProc;
          eval_v(el,(const DOUBLE **)CornersCoord,LocalCoord,vval);
          if (DIM==2)
            sprintf(it,"%lg %lg ",vval[0],vval[1]);
          else
            sprintf(it,"%lg %lg %lg ",vval[0],vval[1],vval[2]);
          strcpy(item+ic,it); ic+=strlen(it);
        }

        sprintf(it,"\n");
        strcpy(item+ic,it); ic+=strlen(it);
        pfile_tagged_puts(pf,item,counter+oe); ic=0;
        counter++;
      }
  }

  pfile_sync(pf);       /* end of segment */

  /********************************/
  /* IN THIS VERSION:				*/
  /* no cell and model data, we   */
  /* are ready !                                        */
  /********************************/

  pfile_close(pf);

  return(OKCODE);
}



/****************************************************************************/
/*																			*/
/* Function:  InitAVS														*/
/*																			*/
/* Purpose:   register all formats for porous media library					*/
/*																			*/
/* Input:     void															*/
/*																			*/
/* Output:	  INT 0: ok														*/
/*				  else line number where error occured						*/
/*																			*/
/****************************************************************************/

INT InitAVS (void)
{
  if (CreateCommand("avs",AVSCommand)==NULL) return (__LINE__);

  return(0);
}
