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

#include "gm.h"
#include "ugenv.h"
#include "ugm.h"
#include "algebra.h"
#include "cmdint.h"
#include "commands.h"
#include "helpmsg.h"
#include "num.h"
#include "shapes.h"
#include "cmdline.h"

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
RCSID("$Header$",UG_RCS_STRING)

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
  VERTEX *vx;                                           /* a vertex pointer							*/
  ELEMENT *el;                                  /* an element pointer						*/

  MULTIGRID *mg;                                /* our multigrid							*/
  char filename[NAMESIZE];      /* file name for output file				*/
  FILE *stream;                                 /* the output file pointer					*/

  INT ns;                                               /* number of scalar eval procs				*/
  INT nv;                                               /* number of vector eval procs				*/
  EVALUES *es[MAXVARIABLES];            /* pointers to scalar eval function desc	*/
  EVECTOR *ev[MAXVARIABLES];            /* pointers to vector eval function desc	*/
  INT ns_cell;                                  /* number of scalar eval procs				*/
  INT nv_cell;                                  /* number of vector eval procs				*/
  EVALUES *es_cell[MAXVARIABLES];       /* pointers to scalar eval function desc*/
  EVECTOR *ev_cell[MAXVARIABLES];       /* pointers to vector eval function desc*/
  char s[NAMESIZE];                             /* name of eval proc						*/
  INT numNodes;                                 /* number of data points					*/
  INT numElements;                              /* number of elements						*/
  PreprocessingProcPtr pre;             /* pointer to prepare function				*/
  ElementEvalProcPtr eval_s;            /* pointer to scalar evaluation function	*/
  ElementVectorProcPtr eval_v;      /* pointer to vector evaluation function	*/
  COORD *CornersCoord[MAX_CORNERS_OF_ELEM];       /* pointers to coordinates    */
  COORD LocalCoord[DIM];                /* is one of the corners local coordinates	*/
  DOUBLE local[DIM];                            /* local coordinate in DOUBLE				*/
  DOUBLE value;                                 /* returned by user eval proc				*/
  DOUBLE x,y,z;                                 /* scalar values							*/
  DOUBLE vval[DIM];                             /* result of vector evaluation function		*/
  time_t ltime;


  /* get current multigrid */
  mg = GetCurrentMultigrid();
  if (mg==NULL)
  {
    PrintErrorMessage('W',"avs","no multigrid open\n");
    return (OKCODE);
  }

  /* scan options */
  ns = nv = ns_cell = nv_cell = 0;
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
        PrintErrorMessage('E',"avs:",buf);
        break;
      }
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
      nv_cell++;
      continue;
    }

  }
  if (ns==0 && nv==0) UserWrite("avs: no variables given, printing mesh data only\n");

  /* get file name and open output file */
  if (sscanf(argv[0],expandfmt(CONCAT3(" avs %",NAMELENSTR,"[ -~]")),filename)!=1)
  {
    PrintErrorMessage('E',"avs","could not read name of logfile");
    return(PARAMERRORCODE);
  }
  stream = fopen(filename,"w");
  if (stream==NULL) return(PARAMERRORCODE);

  /********************************/
  /* TITLE                                              */
  /********************************/

  time(&ltime);
  fprintf(stream,"# AVS data produced by UG\n");
  fprintf(stream,"# multigrid: %s\n",mg->v.name);
  fprintf(stream,"# date:      %s",ctime(&ltime));

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

  /********************************/
  /* (1) now we are ready to write*/
  /* the sizes to the file		*/
  /********************************/

  fprintf(stream,"%d %d %d %d %d\n",
          numNodes,                                     /* total number of nodes to write       */
          numElements,                          /* total number of cells			*/
          ns+nv*DIM,                                    /* number of doubles per node		*/
          ns_cell+nv_cell*DIM,          /* doubles per cell is always zero !*/
          0);                                                   /* doubles per model                            */


  /********************************/
  /* (2) next are the coordinates */
  /* of nodes in point block form */
  /********************************/

  /* clear USED flag in vertices on all levels */
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
      SETUSED(vx,0);

  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (NSONS(el)>0) continue;                        /* process finest level elements only */
      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        vx = MYVERTEX(CORNER(el,i));
        if (USED(vx)) continue;                         /* we have this one already */
        SETUSED(vx,1);                                          /* tag vector as visited */
        x = XC(vx);
        y = YC(vx);
        if (DIM==3)
          z = ZC(vx);
        else
          z = 0.0;
        fprintf(stream,"%d %lg %lg %lg\n",ID(vx),x,y,z);
      }
    }

  /********************************/
  /* (3) next is the cell               */
  /*     connectivity list		*/
  /********************************/

  counter=0;
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (NSONS(el)>0) continue;                        /* process finest level elements only */
      counter++;                                                                /* this is the cell number                        */

      /* cell number and material id */
      fprintf(stream,"%d %d ",counter,1 /*LEVEL(el)*/);

      if (DIM==2)
      {
        switch (CORNERS_OF_ELEM(el))
        {
        case 3 :
          fprintf(stream,"tri %d %d %d\n",
                  ID(MYVERTEX(CORNER(el,0))),
                  ID(MYVERTEX(CORNER(el,1))),
                  ID(MYVERTEX(CORNER(el,2))));
          break;
        case 4 :
          fprintf(stream,"quad %d %d %d %d\n",
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
          fprintf(stream,"tet %d %d %d %d\n",
                  ID(MYVERTEX(CORNER(el,0))),
                  ID(MYVERTEX(CORNER(el,2))),
                  ID(MYVERTEX(CORNER(el,1))),
                  ID(MYVERTEX(CORNER(el,3))));
          break;
        case 5 :
          fprintf(stream,"pyr %d %d %d %d %d\n",
                  ID(MYVERTEX(CORNER(el,4))),
                  ID(MYVERTEX(CORNER(el,0))),
                  ID(MYVERTEX(CORNER(el,1))),
                  ID(MYVERTEX(CORNER(el,2))),
                  ID(MYVERTEX(CORNER(el,3))));
          break;
        case 6 :
          fprintf(stream,"prism %d %d %d %d %d %d\n",
                  ID(MYVERTEX(CORNER(el,3))),
                  ID(MYVERTEX(CORNER(el,4))),
                  ID(MYVERTEX(CORNER(el,5))),
                  ID(MYVERTEX(CORNER(el,0))),
                  ID(MYVERTEX(CORNER(el,1))),
                  ID(MYVERTEX(CORNER(el,2))));
          break;
        case 8 :
          fprintf(stream,"hex %d %d %d %d %d %d %d %d\n",
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
    }

  /********************************/
  /* node data					*/
  /********************************/

  if (ns+nv>0)
  {

    /********************************/
    /* (4) data associated with		*/
    /*     nodes					*/
    /********************************/

    fprintf(stream,"%d ",ns+nv);                                        /* number of components */
    for (k=0; k<ns; k++)
      fprintf(stream,"1 ");                                                     /* scalar component		*/
    for (k=0; k<nv; k++)
      fprintf(stream,"%d ",DIM);                                                /* vector component		*/
    fprintf(stream,"\n");

    /********************************/
    /* (5) component labels			*/
    /*                                                  */
    /********************************/

    /* now all scalar variables */
    for (i=0; i<ns; i++)
      fprintf(stream,"%s, \n",es[i]->v.name);           /* no units */
    for (i=0; i<nv; i++)
      fprintf(stream,"%s, \n",ev[i]->v.name);           /* no units */

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
      if (pre!=NULL) pre("es[v]->v.name",mg);
    }
    for (v=0; v<nv; v++)
    {
      pre    = ev[v]->PreprocessProc;

      /* execute prepare function */
      if (pre!=NULL) pre("ev[v]->v.name",mg);
    }

    /* now the data in point block format */
    /* clear USED flag in vertices on all levels */
    for (k=0; k<=TOPLEVEL(mg); k++)
      for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
        SETUSED(vx,0);
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
          fprintf(stream,"%d ",ID(MYVERTEX(CORNER(el,i))));

          /* scalar components */
          for (v=0; v<ns; v++)
          {
            eval_s = es[v]->EvalProc;
            value = eval_s(el,(const COORD **)CornersCoord,LocalCoord);
            fprintf(stream,"%lg ",value);
          }

          /* vector components */
          for (v=0; v<nv; v++)
          {
            eval_v = ev[v]->EvalProc;
            eval_v(el,(const COORD **)CornersCoord,LocalCoord,vval);
            if (DIM==2)
              fprintf(stream,"%lg %lg ",vval[0],vval[1]);
            else
              fprintf(stream,"%lg %lg %lg ",vval[0],vval[1],vval[2]);
          }

          fprintf(stream,"\n");
        }
      }
  }

  /********************************/
  /* cell data					*/
  /********************************/

  if (ns_cell+nv_cell>0)
  {

    /********************************/
    /* (4) data associated with		*/
    /*     nodes					*/
    /********************************/

    fprintf(stream,"%d ",ns_cell+nv_cell);                                      /* number of components */
    for (k=0; k<ns_cell; k++)
      fprintf(stream,"1 ");                                                     /* scalar component		*/
    for (k=0; k<nv_cell; k++)
      fprintf(stream,"%d ",DIM);                                                /* vector component		*/
    fprintf(stream,"\n");

    /********************************/
    /* (5) component labels			*/
    /*                                                  */
    /********************************/

    /* now all scalar variables */
    for (i=0; i<ns_cell; i++)
      fprintf(stream,"%s, \n",es_cell[i]->v.name);           /* no units */
    for (i=0; i<nv_cell; i++)
      fprintf(stream,"%s, \n",ev_cell[i]->v.name);           /* no units */

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
      if (pre!=NULL) pre("es_cell[v]->v.name",mg);
    }
    for (v=0; v<nv_cell; v++)
    {
      pre    = ev_cell[v]->PreprocessProc;

      /* execute prepare function */
      if (pre!=NULL) pre("ev_cell[v]->v.name",mg);
    }

    /* now the data in point block format */
    counter=0;
    for (k=0; k<=TOPLEVEL(mg); k++)
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if (NSONS(el)>0) continue;                              /* process finest level elements only */
        counter++;                                                                      /* this is the cell id				  */
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
        for (j=0; j<DIM; j++) LocalCoord[j] /= ((COORD)CORNERS_OF_ELEM(el));

        /* write cell id */
        fprintf(stream,"%d ",counter);

        /* scalar components */
        for (v=0; v<ns_cell; v++)
        {
          eval_s = es_cell[v]->EvalProc;
          value = eval_s(el,(const COORD **)CornersCoord,LocalCoord);
          fprintf(stream,"%lg ",value);
        }

        /* vector components */
        for (v=0; v<nv_cell; v++)
        {
          eval_v = ev_cell[v]->EvalProc;
          eval_v(el,(const COORD **)CornersCoord,LocalCoord,vval);
          if (DIM==2)
            fprintf(stream,"%lg %lg ",vval[0],vval[1]);
          else
            fprintf(stream,"%lg %lg %lg ",vval[0],vval[1],vval[2]);
        }

        fprintf(stream,"\n");
      }
  }

  /********************************/
  /* IN THIS VERSION:				*/
  /* no cell and model data, we   */
  /* are ready !                                        */
  /********************************/

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

INT InitAVS (void)
{
  if (CreateCommand("avs",AVSCommand)==NULL) return (__LINE__);

  return(0);
}
