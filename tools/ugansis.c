// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  uganasis.c													*/
/*																			*/
/* Purpose:   converter from ANSIS into GEN format                                              */
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   Nov 5, 1999                                                                               */
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

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define MAX_LEN 128

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

/*  UG_RCS_STRING
   $Header$
 */

#include "gen.h"

int main (int argc, char **argv)
{
  FILE *file;
  char line[MAX_LEN],name[MAX_LEN];
  int nP,nC;
  POINT *P;
  CELL *C;
  int i;

  if (argc<2)     {
    printf("filename required\n");
    return (1);
  }
  strcpy(name,argv[1]);
  strcat(name,".ans");
  file = fopen(name,"r");
  if (file == NULL) {
    printf("cannot open file %s\n",name);
    return (1);
  }
  nP = nC = 0;
  while (!feof(file))
  {
    fgets(line, MAX_LEN, file);
    if (line[0] == 'N') nP++;
    if (line[0] == 'E') nC++;
  }
  printf("nP %d nC %d\n",nP,nC);
  fclose(file);

  P = (POINT *) malloc(nP*sizeof(POINT));
  C = (CELL *) malloc(nC*sizeof(CELL));

  file = fopen(name,"r");
  if (file == NULL) {
    printf("cannot open file %s\n",name);
    return (1);
  }
  nP = nC = 0;
  while (!feof(file))
  {
    fgets(line, MAX_LEN, file);
    if (line[0] == 'N')
    {
      double x[3];
      int id,k;

      sscanf(line,"N,%d,%lf,%lf,%lf",&id,x,x+1,x+2);
      for (k=0; k<3; k++)
        P[nP].x[k] = x[k];
      nP++;
    }
    if (line[0] == 'E') {
      int x[3];
      int id,k;

      sscanf(line,"EN,%d,%d,%d,%d,%d",&id,x,x+1,x+2,x+3);
      for (k=0; k<4; k++)
        C[nC].P[k] = x[k] - 1;
      nC++;

    }
  }
  fclose(file);

  strcpy(name,argv[1]);
  strcat(name,".mesh");
  file = fopen(name,"w");
  fprintf(file,"POINTS:   %d\n",nP);
  for (i=0; i<nP; i++)
    fprintf(file,"%22.15e %22.15e %22.15e\n",
            P[i].x[0],
            P[i].x[1],
            P[i].x[2]);
  fprintf(file,"CELLS:   %d\n",nC);
  for (i=0; i<nC; i++)
    fprintf(file,"%d %d %d %d %d %d\n",
            4,1,C[i].P[0],C[i].P[1],C[i].P[2],C[i].P[3]);
  fclose(file);

  return(0);
}
