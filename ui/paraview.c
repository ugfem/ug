// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      paraview.c                                                    */
/*                                                                          */
/* Purpose:   paraview output file                                          */
/*                                                                          */
/* History:   25.02.2005 begin	                                            */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*									    */
/* include files                                                            */
/*			system include files				    */
/*			application include files			    */
/*								            */
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

#include "paraview.h"
#include "cw.h"


USING_UG_NAMESPACES

/****************************************************************************/
/*								            */
/* defines in the following order					    */
/*									    */
/*	     compile time constants defining static data size (i.e. arrays) */
/*	     other constants						    */
/*	     macros	                                                    */
/*                                                                          */
/****************************************************************************/

#define MAXVARIABLES 50                 /* max number of eval procs */

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)	    */
/*                                                                          */
/****************************************************************************/


/* data for CVS	*/
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/
/*                                                                          */
/* definition of functions						    */
/*			                                                    */
/****************************************************************************/

static INT ParaViewCommand (INT argc, char **argv)
{

  int i, j, k, l, a, coe, numVertices, numElements, offset;

  char item[1024], it[256];   /*item buffers*/
  INT ic = 0;                 /*item lenght*/
  char filename[NAMESIZE];    /*filename for output file*/
  char* c_ptr;
  FILE *outputFile;           /*file pointer for output file*/

  MULTIGRID *mg;              /*our multigrid*/
  ELEMENT *el;                /*an element*/
  VERTEX *vx;                 /*a vertex*/

  EVALUES *theEVal[MAXVARIABLES];          /*pointers to scalar eval function desc*/
  char s[NAMESIZE];
  DOUBLE value;               /*returned by user eval proc*/

  const DOUBLE *x[MAX_CORNERS_OF_ELEM];

  /*************************************************************************/
  /*GET CURRENT MULTIGRID                                                  */
  /*************************************************************************/


  mg = GetCurrentMultigrid();
  if (mg==NULL)
  {
    PrintErrorMessage('W',"dataexplorer","no multigrid open\n");
    return (OKCODE);
  }


  /*************************************************************************/
  /*SCAN OPTIONS                                                           */
  /*************************************************************************/

  int ns = 0;
  for(i=1; i<argc; i++)
  {
    PrintErrorMessageF('E',"paraview:","scan option %d\n",i);
    sscanf(argv[i],"a %s", s);
    PrintErrorMessageF('E',"paraview:","sscanf %d\n",i);
    theEVal[i] = GetElementValueEvalProc(s);
    PrintErrorMessageF('E',"paraview:","theEval[ %d\n",i);
    ns++;
  }

  /*************************************************************************/
  /*GET FILENAME AND OPEN OUTPUTFILE                                       */
  /*************************************************************************/


  if (sscanf(argv[0],expandfmt(CONCAT3(" paraview %",NAMELENSTR,"[ -~]")),
             filename)!=1)
  {
    PrintErrorMessage('E',"paraview","could not read name of output file");
    return(PARAMERRORCODE);
  }

  c_ptr = strrchr(filename, '.');
  if (c_ptr != NULL) memset(c_ptr, '\0', 1);
  strcat(filename, ".vtu");

  outputFile = fopen(filename, "w");

  if(outputFile == NULL)
  {
    PrintErrorMessage('E', "paraview", "could not open output file");
    return(PARAMERRORCODE);
  }

  /**************************************************************************/
  /*TITLE                                                                   */
  /**************************************************************************/


  sprintf(it, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;

  sprintf(it, "\t<UnstructuredGrid>\n");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;

  /**************************************************************************/
  /*WRITE NUMBER OF VERTICES AND ELEMENTS                                   */
  /**************************************************************************/


  for (k=0; k<=TOPLEVEL(mg); k++)
  {
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
    {
      SETUSED(vx,0);
    }
  }

  numVertices = 0;

  for (k=0; k<=TOPLEVEL(mg); k++)
  {
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (!EstimateHere(el)) continue;

      numElements ++;

      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        vx = MYVERTEX(CORNER(el,i));
        if (USED(vx)) continue;
        SETUSED(vx,1);

        ID(vx) = numVertices;

        numVertices++;
      }
    }
  }

  sprintf(it, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", numVertices, numElements);
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;


  /****************************************************************************/
  /*WRITE POINTS                                                              */
  /****************************************************************************/


  sprintf(it, "\t\t\t<Points>\n");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;


  sprintf(it, "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n\t\t\t\t\t");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;


  for (k=0; k<=TOPLEVEL(mg); k++)
  {
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
    {
      SETUSED(vx,0);
    }
  }

  for (k=0; k<=TOPLEVEL(mg); k++)
  {
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (!EstimateHere(el)) continue;

      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        vx = MYVERTEX(CORNER(el,i));
        if (USED(vx)) continue;
        SETUSED(vx,1);

        sprintf(it, "%g %g %g \n", XC(vx), YC(vx), ZC(vx));
        strcpy(item+ic, it);
        ic += strlen(it);
        fputs(it, outputFile);
        ic = 0;
      }
    }
  }

  sprintf(it, "\n\t\t\t\t</DataArray>");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;


  sprintf(it, "\n\t\t\t</Points>\n");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;


  /****************************************************************************/
  /*WRITE CELLS:                                                              */
  /****************************************************************************/


  /****************************************************************************/
  /*WRITE CELL CONNECTIVITY                                                   */
  /****************************************************************************/


  sprintf(it, "\t\t\t<Cells>\n");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;

  sprintf(it, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n\t\t\t\t\t");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;


  for (k=0; k<=TOPLEVEL(mg); k++)
  {
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,k)); vx!=NULL; vx=SUCCV(vx))
    {
      SETUSED(vx,0);
    }
  }


  for (k=0; k<=TOPLEVEL(mg); k++)
  {
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      sprintf(it, "\n");
      fputs(it, outputFile);
      ic = 0;

      if (!EstimateHere(el)) continue;

      for (i=0; i<CORNERS_OF_ELEM(el); i++)
      {
        vx = MYVERTEX(CORNER(el,i));

        sprintf(it, "%d ", ID(vx));
        fputs(it, outputFile);
        ic = 0;
      }
    }
  }

  sprintf(it, "\n\t\t\t\t</DataArray>\n");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;


  /****************************************************************************/
  /*WRITE CELL OFFSETS                                                        */
  /****************************************************************************/


  sprintf(it, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t\t\t\t");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;

  offset = 0;

  for (k=0; k<=TOPLEVEL(mg); k++)
  {
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (!EstimateHere(el)) continue;

      offset += CORNERS_OF_ELEM(el);

      sprintf(it, "%d \n", offset);
      strcpy(item+ic, it);
      ic += strlen(it);
      fputs(it, outputFile);
      ic = 0;

    }
  }


  sprintf(it, "\n\t\t\t\t</DataArray>\n");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;


  /****************************************************************************/
  /*WRITE CELL TYPES                                                          */
  /****************************************************************************/


  sprintf(it, "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n\t\t\t\t\t");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;


  for (k=0; k<=TOPLEVEL(mg); k++)
  {
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el))
    {
      if (!EstimateHere(el)) continue;

      if(TAG(el) == TRIANGLE) sprintf(it, "5 \n");
      else if(TAG(el) == QUADRILATERAL) sprintf(it, "9 \n");
      else if(TAG(el) == TETRAHEDRON) sprintf(it, "10 \n");
      else if(TAG(el) == PYRAMID) sprintf(it, "14 \n");
      else if(TAG(el) == PRISM) sprintf(it, "13 \n");
      else if(TAG(el) == HEXAHEDRON) sprintf(it, "12 \n");
      else sprintf(it, "   TAG NOT FOUND   \n");

      strcpy(item+ic, it);
      ic += strlen(it);
      fputs(it, outputFile);
      ic = 0;
    }
  }


  sprintf(it, "\n\t\t\t\t</DataArray>\n");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;


  sprintf(it, "\t\t\t</Cells>\n");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;


  /****************************************************************************/
  /*WRITE POINT DATA                                                          */
  /****************************************************************************/


  /*   sprintf(it, "\t\t\t<PointData>\n"); */
  /*   strcpy(item+ic, it); */
  /*   ic += strlen(it); */
  /*   fputs(it, outputFile); */
  /*   ic = 0; */

  /*   sprintf(it, "\t\t\t\t<DataArray type=\"Float64\" Name=\"values\" format=\"ascii\">\n\t\t\t\t\t"); */
  /*   strcpy(item+ic, it); */
  /*   ic += strlen(it); */
  /*   fputs(it, outputFile); */
  /*   ic = 0; */


  /* save data from eval procs */

  /*   for (k=0; k<=TOPLEVEL(mg); k++) */
  /*     { */
  /*       PrintErrorMessageF('E',"paraview:","WRITE POINT DATA k %d\n",k); */
  /*       for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,k)); el!=NULL; el=SUCCE(el)) */
  /*    { */
  /*      if (!EstimateHere(el)) continue; */
  /*      coe = CORNERS_OF_ELEM(el); */
  /*      for (l=0; l<coe; l++) */
  /*        { */
  /*          PrintErrorMessageF('E',"paraview:","WRITE POINT DATA l %d\n",l); */
  /*          value = (*(theEVal[0]->EvalProc))(el,x,LOCAL_COORD_OF_ELEM(el,l)); */ /*fehler*/
  /*               PrintErrorMessageF('E',"paraview:","WRITE POINT DATA l %d\n",l); */
  /*          sprintf(it,"\t\t\t\t\%g\n",value); */
  /*          strcpy(item+ic, it); */
  /*          ic += strlen(it); */
  /*          fputs(it, outputFile); */
  /*          ic = 0; */
  /*        } */
  /*    } */
  /*     } */

  /*   sprintf(it, "\n\t\t\t\t</DataArray>\n"); */
  /*   strcpy(item+ic, it); */
  /*   ic += strlen(it); */
  /*   fputs(it, outputFile); */
  /*   ic = 0; */


  /*   sprintf(it, "\t\t\t</PointData>\n"); */
  /*   strcpy(item+ic, it); */
  /*   ic += strlen(it); */
  /*   fputs(it, outputFile); */
  /*   ic = 0; */


  /****************************************************************************/
  /*WRITE END OF FILE AND CLOSE FILE                                          */
  /****************************************************************************/


  sprintf(it, "\t\t</Piece>\n");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;


  sprintf(it, "\t</UnstructuredGrid>\n");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;


  sprintf(it, "</VTKFile>");
  strcpy(item+ic, it);
  ic += strlen(it);
  fputs(it, outputFile);
  ic = 0;


  fclose(outputFile);

  return(0);

}


/****************************************************************************/
/*					                                    */
/* Function:  InitParaView                                                  */
/*	                                                                    */
/* Purpose:                                                                 */
/*	                                                                    */
/* Input:	 void	                                                    */
/*                                                                          */
/* Output:	  INT 0: ok						    */
/*			     else line number where error occured	    */
/*                                                                          */
/****************************************************************************/


INT InitParaView (void)
{
  if (CreateCommand("paraview",ParaViewCommand)==NULL) return (__LINE__);

  return(0);
}
