// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugmgs.c														*/
/*																			*/
/* Purpose:   analyse i/o file for mg from ug                               */
/*																			*/
/* Author:	  Klaus Johannsen												*/
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  email: ug@ica3.uni-stuttgart.de					            */
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "compiler.h"
#include "fileopen.h"
#include "heaps.h"
#include "defaults.h"
#include "general.h"
#include "debug.h"
#include "bio.h"
#include "ugstruct.h"

#include "devices.h"

#include "gm.h"
#include "algebra.h"
#include "misc.h"
#include "ugm.h"
#include "ugio.h"
#include "elements.h"
#include "shapes.h"
#include "mgio.h"
#include "dio.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

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

/* RCS string */
/*static char RCS_ID("$Header$",UG_RCS_STRING);*/

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

int PrintMGFileInfo (char *filename, int *magic_cookie)
{
  MGIO_MG_GENERAL mg_general;

  /* open file and check for mg-file */
  if (Read_OpenMGFile (filename))
  {
    printf("ERROR: cannot open file %s\n",filename);
    return (0);
  }
  if (Read_MG_General(&mg_general))
  {
    CloseMGFile();
    return (2);
  }

  /* print info for mg-file */
  printf("\n############# mg-file ###############\n\n");
  switch (mg_general.mode)
  {
  case BIO_DEBUG :
    printf("Mode:           DEBUG\n");
    break;
  case BIO_ASCII :
    printf("Mode:           ASCII\n");
    break;
  case BIO_BIN :
    printf("Mode:           BINARY\n");
    break;
  default :
    CloseMGFile();
    return (0);
  }
  printf("Version:        %s\n",mg_general.version);
  printf("Identification: %s\n",mg_general.ident);
  printf("Magic Cookie:   %d\n",(int)mg_general.magic_cookie);
  printf("\n");
  printf("Dimension:      %d\n",mg_general.dim);
  printf("BVP:            %s\n",mg_general.DomainName);
  printf("MG-Name:        %s\n",mg_general.MultiGridName);
  printf("Format-Name:    %s\n",mg_general.Formatname);
  printf("Heapsize:       %d kByte\n",mg_general.heapsize);
  printf("\n");
  printf("# Level:        %d\n",mg_general.nLevel);
  printf("# Node:         %d\n",mg_general.nNode);
  printf("# Vertex:       %d\n",mg_general.nPoint);
  printf("# Element:      %d\n",mg_general.nElement);
  printf("\n");

  /* set magic_cookie */
  *magic_cookie = mg_general.magic_cookie;

  /* close file */
  if (CloseMGFile ()) return (0);

  return (1);
}

int PrintDataFileInfo (char *filename, int *magic_cookie)
{
  int i,j;
  DIO_GENERAL dio_general;

  if (Read_OpenDTFile (filename))
  {
    printf("ERROR: cannot open file %s\n",filename);
    return (0);
  }

  /* read general information */
  if (Read_DT_General (&dio_general))
  {
    CloseDTFile();
    return (2);
  }

  /* print info for data-file */
  printf("\n############# data-file  #############\n\n");
  switch (dio_general.mode)
  {
  case BIO_DEBUG :
    printf("Mode:           DEBUG\n");
    break;
  case BIO_ASCII :
    printf("Mode:           ASCII\n");
    break;
  case BIO_BIN :
    printf("Mode:           BINARY\n");
    break;
  default :
    CloseDTFile();
    return (0);
  }
  printf("Version:        %s\n",dio_general.version);
  printf("MG File:        %s\n",dio_general.mgfile);
  printf("Time:           %f\n",(float)dio_general.time);
  printf("TimeStep:       %f\n",(float)dio_general.dt);
  printf("Magic Cookie:   %d\n",(int)dio_general.magic_cookie);
  printf("\n");
  printf("# VecData:      %d\n",dio_general.nVD);
  printf("\n");
  printf("    # comp  |  name  |  type\n");
  printf("  ----------+--------+---------\n");
  for (i=0; i<dio_general.nVD; i++)
  {
    printf("    %3d     |  %-4s  |  ",dio_general.VDncomp[i],dio_general.VDname[i]);
    if (dio_general.VDtype[i]==DIO_SCALAR)
      printf("scal\n");
    if (dio_general.VDtype[i]==DIO_VECTOR)
      printf("vect\n");
    if (dio_general.VDtype[i]==DIO_MULTIPLE_SCALAR)
      printf("mscal\n");
  }
  printf("\n");

  /* set magic_cookie */
  *magic_cookie = dio_general.magic_cookie;

  CloseDTFile();
  return (1);
}

int main (int argc, char **argv)
{
  int i,ret,j;
  int mgmc[20], datamc[20], mglist[20], datalist[20], nmg, ndata, mc;

  if (argc<2 || argc>20)
  {
    printf("usage: ugmgs <file1> [<file2> <file3> ...]\n");
    return (0);
  }

  nmg = ndata = 0;
  for (i=1; i<argc; i++)
  {
    /* print mg-file info if */
    ret = PrintMGFileInfo (argv[i],mgmc+nmg);
    if (ret == 0) return (0);
    if (ret == 1)
    {
      mglist[nmg++] = i;
      continue;
    }

    /* print data-file info if */
    ret = PrintDataFileInfo (argv[i],datamc+ndata);
    if (ret == 0) return (1);
    if (ret == 1)
    {
      datalist[ndata++] = i;
      continue;
    }

    return (0);
  }

  if (nmg==0 || ndata==0) return (0);

  /* check consistency */
  printf("\n####################### consistency  #######################\n\n");
  printf("                          |");
  for (j=0; j<ndata; j++)
    printf("%-30s",argv[datalist[j]]);
  printf("\n-------------------------");
  for (j=0; j<ndata; j++)
    printf("------------------------------");
  printf("\n");
  for (i=0; i<nmg; i++)
  {
    printf("%-25s |",argv[mglist[i]]);
    for (j=0; j<ndata; j++)
      if (mgmc[i]==datamc[j])
        printf("       C       ");
      else
        printf("       -       ");
    printf("\n");
  }
  printf("\n");

  return (0);
}
