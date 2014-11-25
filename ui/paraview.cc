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
/* Usage:     paraview <filename> $a <sol>                                  */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*									    */
/* include files                                                            */
/*			system include files				    */
/*			application include files			    */
/*								            */
/****************************************************************************/

#include <config.h>
#include <cstdio>
#include <cstring>
#include <ctype.h>
#include <cmath>
#include <ctime>

#include "ugdevices.h"
#include "enrol.h"
#include "ugtypes.h"
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


/* max number of eval procs */
#define MAXVARIABLES 50


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


#ifdef ModelP

static INT WriteMetaFile(char* filename, VECDATA_DESC* vd)
{

  /*counters*/
  int p, level, counterEquations, numEquations;

  /*file stuff*/
  char buffer[1024];          /*buffer for outputfile*/
  char* c_ptr;                /*for searching patterns in strings, resp. a dot in filename*/
  FILE* metaOutputFile;       /*file pointer for output file*/

  /*ug topology*/
  MULTIGRID* mg;              /*our multigrid*/
  GRID* g;                    /*the grid*/


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
  /*GET FILENAME AND OPEN METAOUTPUTFILE                                       */
  /*************************************************************************/


  c_ptr = strrchr(filename, '.');
  if (c_ptr != NULL) memset(c_ptr, '\0', 1);

  strcat(filename, ".MetaFile.pvtu");

  metaOutputFile = fopen(filename, "w");

  if(metaOutputFile == NULL)
  {
    PrintErrorMessage('E', "paraview", "could not open output file");
    return(PARAMERRORCODE);
  }


  /**************************************************************************/
  /*TITLE                                                                   */
  /**************************************************************************/


  sprintf(buffer, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fputs(buffer, metaOutputFile);

  sprintf(buffer, "\t<PUnstructuredGrid Ghostlevel=\"0\">\n");
  fputs(buffer, metaOutputFile);


  /****************************************************************************/
  /*WRITE PPOINTS                                                              */
  /****************************************************************************/


  sprintf(buffer, "\t\t\t<PPoints>\n");
  fputs(buffer, metaOutputFile);

  sprintf(buffer, "\t\t\t\t<PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"/>\n\t\t\t\t\t");
  fputs(buffer, metaOutputFile);

  sprintf(buffer, "\n\t\t\t</PPoints>\n");
  fputs(buffer, metaOutputFile);


  /****************************************************************************/
  /*WRITE PPOINT DATA                                                          */
  /****************************************************************************/

  /*vector data*/

  sprintf(buffer, "\t\t\t<PPointData>\n");
  fputs(buffer, metaOutputFile);

  for (level=0; level<=TOPLEVEL(mg); level++)
  {
    g = GRID_ON_LEVEL(mg,level);
    numEquations = VD_NCMPS_IN_TYPE(vd,VTYPE(FIRSTVECTOR(g)));

    for (counterEquations=0; counterEquations<numEquations; counterEquations++)
    {
      sprintf(buffer, "\t\t\t\t<PDataArray type=\"Float64\" Name=\"%c\" format=\"ascii\"/>\n\t\t\t\t\t", VM_COMP_NAME(vd,counterEquations));
      fputs(buffer, metaOutputFile);
    }
  }

  sprintf(buffer, "\n\t\t\t</PPointData>\n");
  fputs(buffer, metaOutputFile);


  /****************************************************************************/
  /*WRITE PIECES (SOURCES)                                                    */
  /****************************************************************************/


  c_ptr = strrchr(filename, '.');
  if (c_ptr != NULL) memset(c_ptr, '\0', 1);

  c_ptr = strrchr(filename, '.');
  if (c_ptr != NULL) memset(c_ptr, '\0', 1);

  for(p=0; p<procs; p++)
  {
    sprintf(buffer, "\t\t\t<Piece Source=\"%s.%d.vtu\"/>\n", filename, p);
    fputs(buffer, metaOutputFile);
  }


  /****************************************************************************/
  /*WRITE END OF FILE AND CLOSE FILE                                          */
  /****************************************************************************/


  sprintf(buffer, "\t</PUnstructuredGrid>\n");
  fputs(buffer, metaOutputFile);

  sprintf(buffer, "</VTKFile>");
  fputs(buffer, metaOutputFile);


  fclose(metaOutputFile);


  return(0);
}

#endif


static INT ParaViewCommand (INT argc, char **argv)
{

  /*counters*/
  int i, counterCorners, level, counterEquations, numVertices, numElements, offset, numEquations;

  /*file stuff*/
  char buffer[1024], it[5];   /*buffers for outputfile*/
  char filename[NAMESIZE];    /*filename for output file*/
  char* c_ptr;                /*for searching patterns in strings, resp. a dot in filename*/
  FILE *outputFile;           /*file pointer for output file*/

  /*ug topology*/
  MULTIGRID *mg;              /*our multigrid*/
  GRID *g;                    /*the grid*/
  ELEMENT *el;                /*an element*/
  VERTEX *vx;                 /*a vertex*/

  /*ug vector solutions*/
  VECTOR *v;                  /*vector*/
  VECDATA_DESC *vd;           /*the vector data descriptor*/
  DOUBLE value;               /*solution value*/

  /*special cases*/
  int ugPrismArray[6];        /*for renumbering vertex IDs*/
  int vtkWedgeArray[6];       /*for renumbering vertex IDs*/


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
  /*GET FILENAME AND OPEN OUTPUTFILE                                       */
  /*************************************************************************/


#ifdef ModelP

  if (sscanf(argv[0],expandfmt(CONCAT3(" paraview %",NAMELENSTR,"[ -~]")),
             filename)!=1)
  {
    PrintErrorMessage('E',"paraview","could not read name of output file");
    return(PARAMERRORCODE);
  }

  c_ptr = strrchr(filename, '.');
  if (c_ptr != NULL) memset(c_ptr, '\0', 1);

  sprintf(it, ".%d", me); /*process id*/
  strcat(filename, it);

  strcat(filename, ".vtu");

  outputFile = fopen(filename, "w");

  if(outputFile == NULL)
  {
    PrintErrorMessage('E', "paraview", "could not open output file");
    return(PARAMERRORCODE);
  }

#endif

#ifndef ModelP

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

#endif


  /**************************************************************************/
  /*TITLE                                                                   */
  /**************************************************************************/


  sprintf(buffer, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fputs(buffer, outputFile);

  sprintf(buffer, "\t<UnstructuredGrid Ghostlevel=\"0\">\n");
  fputs(buffer, outputFile);


  /**************************************************************************/
  /*WRITE NUMBER OF VERTICES AND ELEMENTS                                   */
  /**************************************************************************/


  for (level=0; level<=TOPLEVEL(mg); level++)
  {
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,level)); vx!=NULL; vx=SUCCV(vx))
    {
      SETUSED(vx,0);
    }
  }

  numVertices = 0;
  numElements = 0;

  for (level=0; level<=TOPLEVEL(mg); level++)
  {
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,level)); el!=NULL; el=SUCCE(el))
    {
      if (!EstimateHere(el)) continue;

      numElements ++;

      for (counterCorners=0; counterCorners<CORNERS_OF_ELEM(el); counterCorners++)
      {
        vx = MYVERTEX(CORNER(el,counterCorners));
        if (USED(vx)) continue;
        SETUSED(vx,1);

        ID(vx) = numVertices;

        numVertices++;
      }
    }
  }

  sprintf(buffer, "\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", numVertices, numElements);
  fputs(buffer, outputFile);


  /****************************************************************************/
  /*WRITE POINTS                                                              */
  /****************************************************************************/


  sprintf(buffer, "\t\t\t<Points>\n");
  fputs(buffer, outputFile);

  sprintf(buffer, "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n\t\t\t\t\t");
  fputs(buffer, outputFile);


  for (level=0; level<=TOPLEVEL(mg); level++)
  {
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,level)); vx!=NULL; vx=SUCCV(vx))
    {
      SETUSED(vx,0);
    }
  }

  for (level=0; level<=TOPLEVEL(mg); level++)
  {
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,level)); el!=NULL; el=SUCCE(el))
    {
      if (!EstimateHere(el)) continue;

      for (counterCorners=0; counterCorners<CORNERS_OF_ELEM(el); counterCorners++)
      {
        vx = MYVERTEX(CORNER(el,counterCorners));
        if (USED(vx)) continue;
        SETUSED(vx,1);

        sprintf(buffer, "%g %g %g \n", XC(vx), YC(vx), ZC(vx));
        fputs(buffer, outputFile);
      }
    }
  }

  sprintf(buffer, "\n\t\t\t\t</DataArray>");
  fputs(buffer, outputFile);

  sprintf(buffer, "\n\t\t\t</Points>\n");
  fputs(buffer, outputFile);


  /****************************************************************************/
  /*WRITE CELLS:                                                              */
  /****************************************************************************/


  /****************************************************************************/
  /*WRITE CELL CONNECTIVITY                                                   */
  /****************************************************************************/


  sprintf(buffer, "\t\t\t<Cells>\n");
  fputs(buffer, outputFile);

  sprintf(buffer, "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n\t\t\t\t\t");
  fputs(buffer, outputFile);

  for (level=0; level<=TOPLEVEL(mg); level++)
  {
    for (vx=FIRSTVERTEX(GRID_ON_LEVEL(mg,level)); vx!=NULL; vx=SUCCV(vx))
    {
      SETUSED(vx,0);
    }
  }

  for (level=0; level<=TOPLEVEL(mg); level++)
  {
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,level)); el!=NULL; el=SUCCE(el))
    {
      sprintf(buffer, "\n");
      fputs(buffer, outputFile);

      if (!EstimateHere(el)) continue;

      for (counterCorners=0; counterCorners<CORNERS_OF_ELEM(el); counterCorners++)
      {
        vx = MYVERTEX(CORNER(el,counterCorners));

        /*if the element is a prism: reorder vertex IDs: UG(0,1,2,3,4,5) -> VTK(0,2,1,3,5,4) */

        if(TAG(el) == PRISM)
        {
          ugPrismArray[counterCorners] = ID(vx);

          vtkWedgeArray[0] = ugPrismArray[0];
          vtkWedgeArray[1] = ugPrismArray[2];
          vtkWedgeArray[2] = ugPrismArray[1];
          vtkWedgeArray[3] = ugPrismArray[3];
          vtkWedgeArray[4] = ugPrismArray[5];
          vtkWedgeArray[5] = ugPrismArray[4];

          for(i=0; i<6; i++)
          {
            sprintf(buffer, "%d ", vtkWedgeArray[i]);
            fputs(buffer, outputFile);
          }
        }

        else
        {
          sprintf(buffer, "%d ", ID(vx));
          fputs(buffer, outputFile);
        }
      }
    }
  }

  sprintf(buffer, "\n\t\t\t\t</DataArray>\n");
  fputs(buffer, outputFile);


  /****************************************************************************/
  /*WRITE CELL OFFSETS                                                        */
  /****************************************************************************/


  sprintf(buffer, "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n\t\t\t\t\t");
  fputs(buffer, outputFile);

  offset = 0;

  for (level=0; level<=TOPLEVEL(mg); level++)
  {
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,level)); el!=NULL; el=SUCCE(el))
    {
      if (!EstimateHere(el)) continue;

      offset += CORNERS_OF_ELEM(el);

      sprintf(buffer, "%d \n", offset);
      fputs(buffer, outputFile);
    }
  }

  sprintf(buffer, "\n\t\t\t\t</DataArray>\n");
  fputs(buffer, outputFile);


  /****************************************************************************/
  /*WRITE CELL TYPES                                                          */
  /****************************************************************************/


  sprintf(buffer, "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n\t\t\t\t\t");
  fputs(buffer, outputFile);

  for (level=0; level<=TOPLEVEL(mg); level++)
  {
    for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,level)); el!=NULL; el=SUCCE(el))
    {
      if (!EstimateHere(el)) continue;

      if(TAG(el) == TRIANGLE) sprintf(buffer, "5 \n");
      else if(TAG(el) == QUADRILATERAL) sprintf(buffer, "9 \n");
      else if(TAG(el) == TETRAHEDRON) sprintf(buffer, "10 \n");
      else if(TAG(el) == PYRAMID) sprintf(buffer, "14 \n");
      else if(TAG(el) == PRISM) sprintf(buffer, "13 \n");
      else if(TAG(el) == HEXAHEDRON) sprintf(buffer, "12 \n");
      else sprintf(buffer, "   TAG NOT FOUND   \n");

      fputs(buffer, outputFile);
    }
  }

  sprintf(buffer, "\n\t\t\t\t</DataArray>\n");
  fputs(buffer, outputFile);

  sprintf(buffer, "\t\t\t</Cells>\n");
  fputs(buffer, outputFile);


  /****************************************************************************/
  /*WRITE POINT DATA                                                          */
  /****************************************************************************/

  /*vector data*/

  sprintf(buffer, "\t\t\t<PointData>\n");
  fputs(buffer, outputFile);

  if ((vd = ReadArgvVecDesc(mg,"a",argc,argv))==NULL)
  {
    PrintErrorMessage('E',"paraview","could not read vec symbol");
    return (PARAMERRORCODE);
  }

  for (level=0; level<=TOPLEVEL(mg); level++)
  {
    g = GRID_ON_LEVEL(mg,level);
    numEquations = VD_NCMPS_IN_TYPE(vd,VTYPE(FIRSTVECTOR(g)));

    for (counterEquations=0; counterEquations<numEquations; counterEquations++)
    {
      sprintf(buffer, "\t\t\t\t<DataArray type=\"Float64\" Name=\"%c\" format=\"ascii\">\n\t\t\t\t\t", VM_COMP_NAME(vd,counterEquations));
      fputs(buffer, outputFile);

      for (el=FIRSTELEMENT(GRID_ON_LEVEL(mg,level)); el!=NULL; el=SUCCE(el))
      {
        if (!EstimateHere(el)) continue;

        for (counterCorners=0; counterCorners<CORNERS_OF_ELEM(el); counterCorners++)
        {
          v = NVECTOR(CORNER(el,counterCorners));

          value = VVALUE(v,VD_CMP_OF_TYPE(vd,VTYPE(v),counterEquations));
          sprintf(buffer,"%g\n",value);
          fputs(buffer, outputFile);
        }
      }

      sprintf(buffer, "\n\t\t\t\t</DataArray>\n");
      fputs(buffer, outputFile);
    }
  }

  sprintf(buffer, "\t\t\t</PointData>\n");
  fputs(buffer, outputFile);


  /****************************************************************************/
  /*WRITE END OF FILE AND CLOSE FILE                                          */
  /****************************************************************************/


  sprintf(buffer, "\t\t</Piece>\n");
  fputs(buffer, outputFile);

  sprintf(buffer, "\t</UnstructuredGrid>\n");
  fputs(buffer, outputFile);

  sprintf(buffer, "</VTKFile>");
  fputs(buffer, outputFile);


  fclose(outputFile);



  /****************************************************************************/
  /*WRITE METAFILE                                                            */
  /****************************************************************************/


#ifdef ModelP

  if(me == master)
  {
    WriteMetaFile(filename, vd);
  }

#endif


  return(0);

}


/****************************************************************************/
/*					                                    */
/* Function:  InitParaView                                                  */
/*	                                                                    */
/* Purpose:                                                                 */
/*	                                                                    */
/* Input:     void	                                                    */
/*                                                                          */
/* Output:    INT 0: ok						            */
/*            else line number where error occured	                    */
/*                                                                          */
/****************************************************************************/


INT InitParaView (void)
{
  if (CreateCommand("paraview",ParaViewCommand)==NULL) return (__LINE__);

  return(0);
}
