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
    theGrid = GRID_ON_LEVEL(theMG,k);

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
    theGrid = GRID_ON_LEVEL(theMG,k);

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

static INT TPL_GlobalSumINT (INT i)
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
  INT counter;                      /* for formatting output                    */
  char item[1024],it[256];      /* item buffers                             */
  INT ic=0;                     /* item length                              */
  VECTOR *vc;                                           /* a vector pointer							*/
  ELEMENT *el;                                  /* an element pointer						*/

  MULTIGRID *mg;                                /* our multigrid							*/
  char filename[NAMESIZE];      /* file name for output file				*/
  PFILE *pf;                    /* the output file pointer                  */


  INT nv;                                               /* number of variables (eval functions)		*/
  EVALUES *ev[MAXVARIABLES];            /* pointers to eval function descriptors	*/
  char ev_name[MAXVARIABLES][NAMESIZE];         /* names for eval functions     */
  char s[NAMESIZE];                             /* name of eval proc						*/
  INT numNodes;                                 /* number of data points					*/
  INT numElements;                              /* number of elements						*/
  INT gnumNodes;                /* number of data points globally           */
  INT gnumElements;             /* number of elements globallay             */
  PreprocessingProcPtr pre;             /* pointer to prepare function				*/
  ElementEvalProcPtr eval;              /* pointer to evaluation function			*/
  DOUBLE *CornersCoord[MAX_CORNERS_OF_ELEM];       /* pointers to coordinates    */
  DOUBLE LocalCoord[DIM];               /* is one of the corners local coordinates	*/
  DOUBLE local[DIM];                            /* local coordinate in DOUBLE				*/
  DOUBLE value;                                 /* returned by user eval proc				*/
  INT maxCorners;                               /* dimension dependent                                          */
  INT oe,on,n;

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
  pf = pfile_open(filename);
  if (pf==NULL) return(PARAMERRORCODE);

  /********************************/
  /* TITLE                                              */
  /********************************/

  ic = 0;
  sprintf(it,"TITLE = \"UG TECPLOT OUTPUT\"\n");
  strcpy(item+ic,it); ic+=strlen(it);
  sprintf(it,"VARIABLES = \"X\", \"Y\"");
  strcpy(item+ic,it); ic+=strlen(it);
  if (DIM==3) {
    sprintf(it,", \"Z\"");
    strcpy(item+ic,it); ic+=strlen(it);
  }
  for (i=0; i<nv; i++) {
    sprintf(it,", \"%s\"",ev[i]->v.name);
    strcpy(item+ic,it); ic+=strlen(it);
  }
  sprintf(it,"\n");
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;

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
      if (!EstimateHere(el)) continue;                          /* process finest level elements only */
      numElements++;                                            /* increase element counter */
      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        vc = NVECTOR(CORNER(el,i));
        if (VCFLAG(vc)) continue;                       /* we have this one already */

        VINDEX(vc) = ++numNodes;                        /* number of data points, begins with 1 ! */
        SETVCFLAG(vc,1);                                        /* tag vector as visited */
      }
    }

        #ifdef ModelP
  gnumNodes = TPL_GlobalSumINT(numNodes);
  gnumElements = TPL_GlobalSumINT(numElements);
  on=get_offset(numNodes);
  oe=get_offset(numElements);
  GloballyUniqueIDs(mg);   /* renumber objects */
    #else
  gnumNodes = numNodes;
  gnumElements = numElements;
  oe=on=0;
    #endif


  /********************************/
  /* write ZONE data				*/
  /* uses FEPOINT for data		*/
  /* uses QUADRILATERAL in 2D		*/
  /* and BRICK in 3D				*/
  /********************************/

  /* write zone record header */
  if (DIM==2) sprintf(it,"ZONE N=%d, E=%d, F=FEPOINT, ET=QUADRILATERAL\n",gnumNodes,gnumElements);
  if (DIM==3) sprintf(it,"ZONE N=%d, E=%d, F=FEPOINT, ET=BRICK\n",gnumNodes,gnumElements);
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;

  /* write data in FEPOINT format, i.e. all variables of a node per line*/

  for (k=0; k<=TOPLEVEL(mg); k++)
    for (vc=FIRSTVECTOR(GRID_ON_LEVEL(mg,k)); vc!=NULL; vc=SUCCVC(vc))
      SETVCFLAG(vc,0);           /* clear all flags */

  counter=0;
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (!EstimateHere(el)) continue;                  /* process finest level elements only */

      for (i=0; i<CORNERS_OF_ELEM(el); i++)
        CornersCoord[i] = CVECT(MYVERTEX(CORNER(el,i)));                        /* x,y,z of corners */

      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        vc = NVECTOR(CORNER(el,i));
        if (VCFLAG(vc)) continue;                       /* we have this one alre ady */
        SETVCFLAG(vc,1);                                /* tag vector as visited */

        sprintf(it,"%g",(double)XC(MYVERTEX(CORNER(el,i))));
        strcpy(item+ic,it); ic+=strlen(it);
        sprintf(it," %g",(double)YC(MYVERTEX(CORNER(el,i))));
        strcpy(item+ic,it); ic+=strlen(it);
        if (DIM == 3)
        {
          sprintf(it," %g",(double)ZC(MYVERTEX(CORNER(el,i))));
          strcpy(item+ic,it); ic+=strlen(it);
        }

        /* now all the user variables */

        /* get local coordinate of corner */
        LocalCornerCoordinates(DIM,TAG(el),i,local);
        for (j=0; j<DIM; j++) LocalCoord[j] = local[j];

        for (v=0; v<nv; v++)
        {
          pre =  ev[v]->PreprocessProc;
          eval = ev[v]->EvalProc;

          /* execute prepare function */
          /* This is not really equivalent to
             the FEBLOCK-version sinc we call "pre" more
             often than there. D.Werner */

          if (pre!=NULL) pre(ev_name[v],mg);

          /* call eval function */
          value = eval(el,(const DOUBLE **)CornersCoord,LocalCoord);
          sprintf(it," %g",value);
          strcpy(item+ic,it); ic+=strlen(it);
        }
        sprintf(it,"\n");
        strcpy(item+ic,it); ic+=strlen(it);
        pfile_tagged_puts(pf,item,counter+on);
        counter++;
      }
    }
  pfile_sync(pf);       /* end of segment */

  sprintf(it,"\n");
  strcpy(item+ic,it); ic+=strlen(it);
  pfile_master_puts(pf,item); ic=0;

  /* finally write the connectivity list */
  counter=0;
  for (k=0; k<=TOPLEVEL(mg); k++)
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (!EstimateHere(el)) continue;           /* process finest level elements only */

      switch(DIM) {
      case 2 :
        switch(TAG(el)) {
        case TRIANGLE :
          sprintf(it,"%d %d %d %d\n",
                  VINDEX(NVECTOR(CORNER(el,0))),
                  VINDEX(NVECTOR(CORNER(el,1))),
                  VINDEX(NVECTOR(CORNER(el,2))),
                  VINDEX(NVECTOR(CORNER(el,2)))
                  );
          break;
        case QUADRILATERAL :
          sprintf(it,"%d %d %d %d\n",
                  VINDEX(NVECTOR(CORNER(el,0))),
                  VINDEX(NVECTOR(CORNER(el,1))),
                  VINDEX(NVECTOR(CORNER(el,2))),
                  VINDEX(NVECTOR(CORNER(el,3)))
                  );
          break;
        default :
          UserWriteF("tecplot: unknown 2D element type with tag(el) = %d detected. Aborting further processing of command tecplot\n", TAG(el));
          return CMDERRORCODE;
          break;
        }
        break;
      case 3 :
        switch(TAG(el)) {
        case HEXAHEDRON :
          sprintf(it,"%d %d %d %d "
                  "%d %d %d %d\n",
                  VINDEX(NVECTOR(CORNER(el,0))),
                  VINDEX(NVECTOR(CORNER(el,1))),
                  VINDEX(NVECTOR(CORNER(el,2))),
                  VINDEX(NVECTOR(CORNER(el,3))),
                  VINDEX(NVECTOR(CORNER(el,4))),
                  VINDEX(NVECTOR(CORNER(el,5))),
                  VINDEX(NVECTOR(CORNER(el,6))),
                  VINDEX(NVECTOR(CORNER(el,7)))
                  );
          break;
        case TETRAHEDRON :
          sprintf(it,"%d %d %d %d "
                  "%d %d %d %d\n",
                  VINDEX(NVECTOR(CORNER(el,0))),
                  VINDEX(NVECTOR(CORNER(el,1))),
                  VINDEX(NVECTOR(CORNER(el,2))),
                  VINDEX(NVECTOR(CORNER(el,2))),
                  VINDEX(NVECTOR(CORNER(el,3))),
                  VINDEX(NVECTOR(CORNER(el,3))),
                  VINDEX(NVECTOR(CORNER(el,3))),
                  VINDEX(NVECTOR(CORNER(el,3)))
                  );
          break;
        case PYRAMID :
          sprintf(it,"%d %d %d %d "
                  "%d %d %d %d\n",
                  VINDEX(NVECTOR(CORNER(el,0))),
                  VINDEX(NVECTOR(CORNER(el,1))),
                  VINDEX(NVECTOR(CORNER(el,2))),
                  VINDEX(NVECTOR(CORNER(el,3))),
                  VINDEX(NVECTOR(CORNER(el,4))),
                  VINDEX(NVECTOR(CORNER(el,4))),
                  VINDEX(NVECTOR(CORNER(el,4))),
                  VINDEX(NVECTOR(CORNER(el,4)))
                  );
          break;
        case PRISM :
          sprintf(it,"%d %d %d %d "
                  "%d %d %d %d\n",
                  VINDEX(NVECTOR(CORNER(el,0))),
                  VINDEX(NVECTOR(CORNER(el,1))),
                  VINDEX(NVECTOR(CORNER(el,2))),
                  VINDEX(NVECTOR(CORNER(el,2))),
                  VINDEX(NVECTOR(CORNER(el,3))),
                  VINDEX(NVECTOR(CORNER(el,4))),
                  VINDEX(NVECTOR(CORNER(el,5))),
                  VINDEX(NVECTOR(CORNER(el,5)))
                  );
          break;
        default :
          UserWriteF("tecplot: unknown 3D element type with tag(el) = %d detected. Aborting further processing of command tecplot\n", TAG(el));
          return CMDERRORCODE;
          break;
        }
        break;
      }
      strcpy(item+ic,it); ic+=strlen(it);
      pfile_tagged_puts(pf,item,counter+oe); ic=0;
      counter++;

    }

  pfile_sync(pf);       /* end of segment */

  /********************************/
  /* GEOMETRY                                   */
  /* we will do this later, since */
  /* domain interface will change */
  /********************************/

  pfile_close(pf);

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
