// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      teti.c		                                                */
/*                                                                          */
/* Purpose:   grape hierarchical file format for tets						*/
/*                                                                          */
/* Author:	  Peter Bastian                                                                                         */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de					                */
/*																			*/
/* History:   27.11.96 begin, ug version 3.4, derived from avs.c			*/
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

#include "ugstruct.h"
#include "devices.h"
#include "enrol.h"
#include "compiler.h"
#include "misc.h"
#include "general.h"

#include "gm.h"
#include "ugenv.h"
#include "ugm.h"
#include "algebra.h"
#include "cmdint.h"
#include "commands.h"
#include "helpmsg.h"
#include "shapes.h"
#include "cmdline.h"
#include "refine.h"             /* this is not clean ! */

#include "teti.h"

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

/****************************************************************************/
/*																			*/
/* definition of functions													*/
/*																			*/
/****************************************************************************/

static INT HierElement (FILE *stream, ELEMENT *el)
{
  int rule,i;
  ELEMENT *sons[MAX_SONS];

  /* write refined flag */
  if (NSONS(el)==0) {
    fprintf(stream,"0\n");
    return(0);
  }
  else
    fprintf(stream,"1\n");

  /* determine refinement rule */
  switch (REFINE(el))
  {
  case 239 : rule = 1; break;
  case 240 : rule = 2; break;
  case 241 : rule = 3; break;
  default :
    UserWrite("teti: refinement rule not allowed !\n");
    return(1);
  }
  fprintf(stream,"%d\n",rule);

  /* write ids of new nodes */
  GetSons(el,sons);
  switch (rule)
  {
  case 1 :
    fprintf(stream,"%d %d %d %d %d %d\n",
            ID(MYVERTEX(CORNER(sons[0],1))),                       /* 4 */
            ID(MYVERTEX(CORNER(sons[1],1))),                       /* 5 */
            ID(MYVERTEX(CORNER(sons[0],2))),                       /* 6 */
            ID(MYVERTEX(CORNER(sons[0],3))),                       /* 7 */
            ID(MYVERTEX(CORNER(sons[1],2))),                       /* 8 */
            ID(MYVERTEX(CORNER(sons[2],2))));                      /* 9 */
    break;

  case 2 :
    fprintf(stream,"%d %d %d %d %d %d\n",
            ID(MYVERTEX(CORNER(sons[1],0))),                       /* 4 */
            ID(MYVERTEX(CORNER(sons[0],0))),                       /* 5 */
            ID(MYVERTEX(CORNER(sons[0],1))),                       /* 6 */
            ID(MYVERTEX(CORNER(sons[2],0))),                       /* 7 */
            ID(MYVERTEX(CORNER(sons[1],2))),                       /* 8 */
            ID(MYVERTEX(CORNER(sons[0],2))));                      /* 9 */
    break;

  case 3 :
    fprintf(stream,"%d %d %d %d %d %d\n",
            ID(MYVERTEX(CORNER(sons[1],0))),                       /* 4 */
            ID(MYVERTEX(CORNER(sons[1],1))),                       /* 5 */
            ID(MYVERTEX(CORNER(sons[2],2))),                       /* 6 */
            ID(MYVERTEX(CORNER(sons[0],0))),                       /* 7 */
            ID(MYVERTEX(CORNER(sons[0],1))),                       /* 8 */
            ID(MYVERTEX(CORNER(sons[0],2))));                      /* 9 */
    break;
  }

  /* now do sons recursively */
  for (i=0; i<NSONS(el); i++)
    HierElement(stream,sons[i]);

  return(0);
}


static INT TetiCommand (INT argc, char **argv)
{
  INT i,j,k,v;                                  /* counters etc.							*/
  char buf[128];                                /* for messages								*/
  VERTEX *vx;                                           /* a vertex pointer							*/
  ELEMENT *el;                                  /* an element pointer						*/

  MULTIGRID *mg;                                /* our multigrid							*/
  char filename[NAMESIZE];      /* file name for output file				*/
  char dataname[NAMESIZE];      /* file name for data output file			*/
  FILE *stream;                                 /* the output file pointer					*/

  char s[NAMESIZE];                             /* name of eval proc						*/
  INT ns;                                               /* number of scalar eval procs				*/
  EVALUES *es[MAXVARIABLES];            /* pointers to scalar eval function desc	*/
  INT numNodes;                                 /* number of data points					*/
  INT numElements;                              /* number of elements						*/
  PreprocessingProcPtr pre;             /* pointer to prepare function				*/
  ElementEvalProcPtr eval_s;            /* pointer to scalar evaluation function	*/
  DOUBLE *CornersCoord[MAX_CORNERS_OF_ELEM];       /* pointers to coordinates    */
  DOUBLE LocalCoord[DIM];               /* is one of the corners local coordinates	*/
  DOUBLE local[DIM];                            /* local coordinate in DOUBLE				*/
  DOUBLE *values;                               /* returned by user eval proc				*/
  INT frame;                                            /* frame number								*/
  double timeval;                               /* time value								*/

  if (DIM!=3) {
    UserWrite("teti: only for tetrahedra !\n");
    return(0);
  }

  /* get current multigrid */
  mg = GetCurrentMultigrid();
  if (mg==NULL)
  {
    PrintErrorMessage('W',"teti","no multigrid open\n");
    return (OKCODE);
  }

  /* scan options */
  ns = 0;
  for(i=1; i<argc; i++)
  {
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
        PrintErrorMessage('E',"teti:",buf);
        break;
      }
      ns++;
      continue;
    }
  }
  if (ns==0) UserWrite("teti: no variables given, printing mesh data only\n");

  /* get file name and open output file */
  if (sscanf(argv[0]," teti %s %d",filename,&frame)!=2)
  {
    PrintErrorMessage('E',"teti","could not read filename");
    return(PARAMERRORCODE);
  }
  stream = fopen(filename,"w");
  if (stream==NULL) return(PARAMERRORCODE);

  /********************************/
  /* TITLE                                              */
  /********************************/

  /* no title */

  /********************************/
  /* renumber objects				*/
  /********************************/

  RenumberMultiGrid(mg);       /* ids are from level 0 to J according to list */

  fprintf(stream,"%d\n",VIDCNT(mg));
  fprintf(stream,"%d\n",EIDCNT(mg));

  /********************************/
  /* Write coarse grid			*/
  /********************************/

  /* count vertices */
  numNodes = 0;
  for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,0)); vx!=NULL; vx=SUCCV(vx))
    numNodes++;
  fprintf(stream,"%d\n",numNodes);

  /* write vertex coordinates */
  for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,0)); vx!=NULL; vx=SUCCV(vx))
    fprintf(stream,"%lg %lg %lg\n",XC(vx),YC(vx),ZC(vx));

  /* count elements */
  numElements = 0;
  for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,0)); el!=NULL; el=SUCCE(el))
    numElements++;
  fprintf(stream,"%d\n",numElements);

  /* coarse grid elements and neighbors */
  for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,0)); el!=NULL; el=SUCCE(el))
  {
    for (i=0; i<CORNERS_OF_ELEM(el); i++)             /* vertex ids */
      fprintf(stream,"%d ",ID(MYVERTEX(CORNER(el,i))));
    for (i=0; i<SIDES_OF_ELEM(el); i++)              /* neighbor ids */
      if (NBELEM(el,i)==NULL)
        fprintf(stream,"%d ",-1);
      else
        fprintf(stream,"%d ",ID(NBELEM(el,i)));
    fprintf(stream,"\n");
  }

  /********************************/
  /* hierarchy					*/
  /********************************/

  for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,0)); el!=NULL; el=SUCCE(el))
    HierElement(stream,el);

  fclose(stream);

  if (ns<=0) return(0);

  /* open data file */
  sprintf(dataname,"%s.%04d",filename,frame);
  stream = fopen(dataname,"w");
  if (stream==NULL) return(PARAMERRORCODE);

  /* write time */
  if (GetStringValueDouble("TIME",&timeval)!=0)
    timeval = (double) frame;

  fprintf(stream,"%lg\n",timeval);

  /* jumps flag is always zero */
  fprintf(stream,"%d\n",0);

  /* number of points */
  fprintf(stream,"%d\n",VIDCNT(mg));

  /* number of components */
  fprintf(stream,"%d\n",ns);

  /* names of components */
  for (v=0; v<ns; v++)
    fprintf(stream,"%s\n",es[v]->v.name);

  /********************************/
  /* node data					*/
  /********************************/

  if (ns>0)
  {
    /* Strategy: allocate array of size nv and fill it */
    Mark(MGHEAP(mg),FROM_BOTTOM);
    values = GetMem(MGHEAP(mg),VIDCNT(mg)*sizeof(DOUBLE),FROM_BOTTOM);
    if (values==NULL) {
      UserWrite("teti: not enough memory\n");
      Release(MGHEAP(mg),FROM_BOTTOM);
      fclose(stream);
      return(1);
    }

    /* loop over all evaluate functions */
    for (v=0; v<ns; v++)
    {
      /* execute prepare function */
      pre    = es[v]->PreprocessProc;
      eval_s = es[v]->EvalProc;
      if (pre!=NULL) pre("es[v]->v.name",mg);

      /* clear used flags in vertices */
      for (k=0; k<=TOPLEVEL(mg); k++)
        for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
          SETUSED(vx,0);

      /* fill values array with scalar data */
      for (k=0; k<=TOPLEVEL(mg); k++)
        for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
        {
          if (NSONS(el)>0) continue;                                    /* process finest level elements only */
          for (i=0; i<CORNERS_OF_ELEM(el); i++)
            CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));                             /* x,y,z of corners */
          for (i=0; i<CORNERS_OF_ELEM(el); i++)
          {
            vx = MYVERTEX(CORNER(el,i));
            if (USED(vx)) continue;                                     /* we have this one already */
            SETUSED(vx,1);                                                      /* tag vector as visited */

            /* get local coordinate of corner */
            LocalCornerCoordinates(DIM,TAG(el),i,local);
            for (j=0; j<DIM; j++) LocalCoord[j] = local[j];

            /* evaluate plot function */
            values[ID(vx)] = eval_s(el,(const DOUBLE **)CornersCoord,LocalCoord);
          }
        }

      /* write values array */
      for (i=0; i<VIDCNT(mg); i++)
        fprintf(stream,"%lg\n",values[i]);
    }

    /* free memory */
    Release(MGHEAP(mg),FROM_BOTTOM);
  }

  fclose(stream);

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

INT InitTeti (void)
{
  if (CreateCommand("teti",TetiCommand)==NULL) return (__LINE__);

  return(0);
}
