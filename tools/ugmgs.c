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

#include "switch.h"
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

/****************************************************************************/
/*																			*/
/* Function:  MGSetVectorClasses											*/
/*																			*/
/* Purpose:   Returns highest vector class of a dof on next level			*/
/*																			*/
/* Input:	  *theElement													*/
/*																			*/
/* Output:	  INT															*/
/*																			*/
/****************************************************************************/

int main (int argc, char **argv)
{
  MGIO_MG_GENERAL mg_general;
  DIO_GENERAL dio_general;

  if (argc!=2 && argc!=3)
  {
    printf("usage: ugmgs <mg-file> [<data-file>]\n");
    return (0);
  }

  /* open file */
  if (Read_OpenMGFile (argv[1]))
  {
    printf("ERROR: cannot open file %s\n",argv[1]);
    return (1);
  }
  if (Read_MG_General(&mg_general))
  {
    printf("ERROR: file %s is probably not an i/o-file\n",argv[1]);
    CloseMGFile();
    return (1);
  }

  /* print specification */
  printf("\n############# ug-i/o specification ###########\n\n");
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
    return (1);
  }
  printf("Version:        %s\n",mg_general.version);
  printf("Identification: %s\n",mg_general.ident);
  printf("\n");
  printf("Dimension:      %d\n",mg_general.dim);
  printf("Domain:         %s\n",mg_general.DomainName);
  printf("MG-Name:        %s\n",mg_general.MultiGridName);
  printf("Format-Name:    %s\n",mg_general.Formatname);
  printf("Heapsize:       %d\n",mg_general.heapsize);
  printf("\n");
  printf("# Level:        %d\n",mg_general.nLevel);
  printf("# Node:         %d\n",mg_general.nNode);
  printf("# Vertex:       %d\n",mg_general.nPoint);
  printf("# Element:      %d\n",mg_general.nElement);
  printf("\n");

  /* close file */
  if (CloseMGFile ()) return (1);

  /* return if only mg-file */
  if (argc == 2) return (0);

  /* open file */
  if (Read_OpenDTFile (argv[2]))
  {
    printf("ERROR: cannot open file %s\n",argv[2]);
    return (1);
  }

  /* read general information */
  if (Read_DT_General (&dio_general))
  {
    printf("ERROR: file %s is probably not an data-file\n",argv[2]);
    CloseDTFile();
    return (1);
  }

  /* output */
  printf("--------------- consistency ------------------\n\n");
  if (mg_general.magic_cookie == dio_general.magic_cookie)
    printf("<%s> and <%s> ARE consistent\n",argv[1],argv[2]);
  else
    printf("<%s> and <%s> ARE NOT consistent\n",argv[1],argv[2]);
  printf("\n");

  return (0);
}
