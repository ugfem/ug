// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      tecplot.c                                                     */
/*                                                                          */
/* Purpose:   tecplot output                                                                                    */
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
/* History:   29.06.95 begin, ug version 3.0								*/
/*            22.08.96 revised version for 3D,local refinement,user eval pr */
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

#include "ugdevices.h"
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

#include "tecplot.h"

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/

#define MAXVARIABLES    20                      /* max number of eval procs				*/
#define VALUES_PER_LINE 10                      /* number of data values per line		*/

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

/****************************************************************************/
/*D
   tecplot - file output in tecplot format

   DESCRIPTION:
   ...

   KEYWORDS:
   graphics, plot, file, output, tecplot
   D*/
/****************************************************************************/

static INT TecplotCommand (INT argc, char **argv)
{
  INT i,j,k,v;                                  /* counters etc.							*/
  INT counter;                                  /* for formatting output					*/
  VECTOR *vc;                                           /* a vector pointer							*/
  ELEMENT *el;                                  /* an element pointer						*/

  MULTIGRID *mg;                                /* our multigrid							*/
  char filename[NAMESIZE];      /* file name for output file				*/
  FILE *stream;                                 /* the output file pointer					*/

  INT nv;                                               /* number of variables (eval functions)		*/
  EVALUES *ev[MAXVARIABLES];            /* pointers to eval function descriptors	*/
  char ev_name[MAXVARIABLES][NAMESIZE];         /* names for eval functions     */
  char s[NAMESIZE];                             /* name of eval proc						*/
  INT numNodes;                                 /* number of data points					*/
  INT numElements;                              /* number of elements						*/
  PreprocessingProcPtr pre;             /* pointer to prepare function				*/
  ElementEvalProcPtr eval;              /* pointer to evaluation function			*/
  DOUBLE *CornersCoord[MAX_CORNERS_OF_ELEM];       /* pointers to coordinates    */
  DOUBLE LocalCoord[DIM];               /* is one of the corners local coordinates	*/
  DOUBLE local[DIM];                            /* local coordinate in DOUBLE				*/
  DOUBLE value;                                 /* returned by user eval proc				*/
  INT maxCorners;                               /* dimension dependent                                          */

  INT saveGeometry;                             /* save geometry flag						*/


  /* get current multigrid */
  mg = GetCurrentMultigrid();
  if (mg==NULL)
  {
    PrintErrorMessage('W',"tecplot","no multigrid open\n");
    return (OKCODE);
  }

  /* scan options */
  nv = 0; saveGeometry = 0;
  for(i=1; i<argc; i++)
  {
    switch(argv[i][0])
    {
    case 'e' :            /* read eval proc */
      if (nv>=MAXVARIABLES)
      {
        PrintErrorMessage('E',"tecplot","too many variables specified\n");
        break;
      }
      sscanf(argv[i],"e %s", s);
      ev[nv] = GetElementValueEvalProc(s);
      if (ev[nv]==NULL)
      {
        PrintErrorMessageF('E',"tecplot","could not find eval proc %s\n",s);
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
      break;

    case 'g' :
      sscanf(argv[i],"g %d", &saveGeometry);
      if (saveGeometry<0) saveGeometry=0;
      if (saveGeometry>1) saveGeometry=1;
      break;
    }
  }
  if (nv==0) UserWrite("tecplot: no variables given, printing mesh data only\n");

  /* get file name and open output file */
  if (sscanf(argv[0],expandfmt(CONCAT3(" tecplot %",NAMELENSTR,"[ -~]")),filename)!=1)
  {
    PrintErrorMessage('E',"tecplot","could not read name of logfile");
    return(PARAMERRORCODE);
  }
  stream = fopen(filename,"w");
  if (stream==NULL) return(PARAMERRORCODE);

  /********************************/
  /* TITLE                                              */
  /********************************/

  fprintf(stream,"TITLE = \"UG TECPLOT OUTPUT\"\n");
  fprintf(stream,"VARIABLES = \"X\", \"Y\"");
  if (DIM==3) fprintf(stream,", \"Z\"");
  for (i=0; i<nv; i++)
    fprintf(stream,", \"%s\"",ev[i]->v.name);
  fprintf(stream,"\n");

  /********************************/
  /* compute sizes				*/
  /********************************/

  /* clear VCFLAG on all levels */
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (vc=FIRSTVECTOR(GRID_ON_LEVEL(mg,k)); vc!=NULL; vc=SUCCVC(vc))
      SETVCFLAG(vc,0);

  /* run thru all levels of elements and set index */
  numNodes = numElements = 0;
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (NSONS(el)>0) continue;                        /* process finest level elements only */
      numElements++;                                            /* increase element counter */
      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        vc = NVECTOR(CORNER(el,i));
        if (VCFLAG(vc)) continue;                       /* we have this one already */

        VINDEX(vc) = ++numNodes;                        /* number of data points, begins with 1 ! */
        SETVCFLAG(vc,1);                                        /* tag vector as visited */
      }
    }

  /********************************/
  /* write ZONE data				*/
  /* uses FEBLOCK for data		*/
  /* uses QUADRILATERAL in 2D		*/
  /* and BRICK in 3D				*/
  /********************************/

  /* write zone record header */
  if (DIM==2)
    fprintf(stream,"ZONE N=%d, E=%d, F=FEBLOCK, ET=QUADRILATERAL\n",numNodes,numElements);
  if (DIM==3)
    fprintf(stream,"ZONE N=%d, E=%d, F=FEBLOCK, ET=BRICK\n",numNodes,numElements);

  /* write data in FEBLOCK format, i.e. all values of one variable in sequence */
  counter = 1;      /* reset counter needed for fixing the number of values per line */

  /* first the X coordinate */
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (vc=FIRSTVECTOR(GRID_ON_LEVEL(mg,k)); vc!=NULL; vc=SUCCVC(vc))
      SETVCFLAG(vc,0);           /* clear all flags */
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (NSONS(el)>0) continue;                        /* process finest level elements only */
      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        vc = NVECTOR(CORNER(el,i));
        if (VCFLAG(vc)) continue;                       /* we have this one already */
        SETVCFLAG(vc,1);                                        /* tag vector as visited */
        fprintf(stream,"%g ",(double)XC(MYVERTEX(CORNER(el,i))));
        counter++;                                                              /* count values	*/
        if (counter%VALUES_PER_LINE==0)
          fprintf(stream,"\n");
      }
    }

  /* then the Y coordinate */
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (vc=FIRSTVECTOR(GRID_ON_LEVEL(mg,k)); vc!=NULL; vc=SUCCVC(vc))
      SETVCFLAG(vc,0);           /* clear all flags */
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (NSONS(el)>0) continue;                        /* process finest level elements only */
      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        vc = NVECTOR(CORNER(el,i));
        if (VCFLAG(vc)) continue;                       /* we have this one already */
        SETVCFLAG(vc,1);                                        /* tag vector as visited */
        fprintf(stream,"%g ",(double)YC(MYVERTEX(CORNER(el,i))));
        counter++;                                                              /* count values	*/
        if (counter%VALUES_PER_LINE==0)
          fprintf(stream,"\n");
      }
    }

  /* then the Z coordinate */
  if (DIM==3) {
    for (k=0; k<=TOPLEVEL(mg); k++)
      for (vc=FIRSTVECTOR(GRID_ON_LEVEL(mg,k)); vc!=NULL; vc=SUCCVC(vc))
        SETVCFLAG(vc,0);                 /* clear all flags */
    for (k=0; k<=TOPLEVEL(mg); k++)
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if (NSONS(el)>0) continue;                              /* process finest level elements only */
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
        {
          vc = NVECTOR(CORNER(el,i));
          if (VCFLAG(vc)) continue;                             /* we have this one already */
          SETVCFLAG(vc,1);                                              /* tag vector as visited */
          fprintf(stream,"%g ",(double)ZC(MYVERTEX(CORNER(el,i))));
          counter++;                                                                    /* count values	*/
          if (counter%VALUES_PER_LINE==0)
            fprintf(stream,"\n");
        }
      }
  }

  /* now all the user variables */
  for (v=0; v<nv; v++)
  {
    pre =  ev[v]->PreprocessProc;
    eval = ev[v]->EvalProc;

    /* execute prepare function */
    if (pre!=NULL) pre(ev_name[v],mg);

    /* now the data */
    for (k=0; k<=TOPLEVEL(mg); k++)
      for (vc=FIRSTVECTOR(GRID_ON_LEVEL(mg,k)); vc!=NULL; vc=SUCCVC(vc))
        SETVCFLAG(vc,0);                 /* clear all flags */
    for (k=0; k<=TOPLEVEL(mg); k++)
      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
      {
        if (NSONS(el)>0) continue;                              /* process finest level elements only */
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
          CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));                       /* x,y,z of corners */
        for (i=0; i<CORNERS_OF_ELEM(el); i++)
        {
          vc = NVECTOR(CORNER(el,i));
          if (VCFLAG(vc)) continue;                             /* we have this one already */
          SETVCFLAG(vc,1);                                              /* tag vector as visited */

          /* get local coordinate of corner */
          LocalCornerCoordinates(DIM,TAG(el),i,local);
          for (j=0; j<DIM; j++) LocalCoord[j] = local[j];

          /* call eval function */
          value = eval(el,(const DOUBLE **)CornersCoord,LocalCoord);
          fprintf(stream,"%g ",value);
          counter++;                                                                    /* count values	*/
          if (counter%VALUES_PER_LINE==0)
            fprintf(stream,"\n");
        }
      }
  }
  fprintf(stream,"\n");

  /* finally write the connectivity list */
  if (DIM==2) maxCorners=4;else maxCorners=8;
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (NSONS(el)>0) continue;                        /* process finest level elements only */

      /* write indices of corners of element */
      for (i=0; i<CORNERS_OF_ELEM(el); i++)
        fprintf(stream,"%d ",VINDEX(NVECTOR(CORNER(el,i))));

      /* now we have i=#corners (e.g. 3 for triangle) */
      j = VINDEX(NVECTOR(CORNER(el,i-1)));                   /* last index */

      /* fill up with last index */
      for ( ; i<maxCorners; i++)
        fprintf(stream,"%d ",j);

      /* end of line is needed */
      fprintf(stream,"\n");
    }

  /********************************/
  /* GEOMETRY                                   */
  /* we will do this later, since */
  /* domain interface will change */
  /********************************/

  fclose(stream);

  return(OKCODE);
}



/****************************************************************************/
/*																			*/
/* Function:  InitTecplot													*/
/*																			*/
/* Purpose:   register all formats for porous media library					*/
/*																			*/
/* Input:     void															*/
/*																			*/
/* Output:	  INT 0: ok														*/
/*				  else line number where error occured						*/
/*																			*/
/****************************************************************************/

INT InitTecplot (void)
{
  if (CreateCommand("tecplot",TecplotCommand)==NULL) return (__LINE__);

  return(0);
}
