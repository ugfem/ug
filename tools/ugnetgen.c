// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugnetgen.c													*/
/*																			*/
/* Purpose:   converter from NETGEN into MARC format                                            */
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

int main (int argc, char **argv)
{
  FILE *file;
  char line[MAX_LEN],name[MAX_LEN];
  int i,k,nv,nf,ne,m;
  int *f,*e;
  double *c;

  if (argc<2)     {
    printf("filename required\n");
    return (1);
  }
  strcpy(name,argv[1]);
  strcat(name,".mesh");
  file = fopen(name,"r");
  if (file == NULL) {
    printf("cannot open file %s\n",name);
    return (1);
  }

  fgets(line, MAX_LEN, file);
  sscanf(line,"%d",&nv);
  printf("nv %d\n",nv);
  c = (double *) malloc(nv*sizeof(double)*3);
  for (i=0; i<nv; i++) {
    fgets(line, MAX_LEN, file);
    sscanf(line,"%lf %lf %lf",c+3*i,c+3*i+1,c+3*i+2);
  }
  for (i=0; i<-nv; i++) {
    k = i;
    printf("%5d%8.5f+0%8.4f+0%8.4f+0\n",
           i+1,c[3*k],c[3*k+1],c[3*k+2]);
  }

  fgets(line, MAX_LEN, file);
  sscanf(line,"%d",&ne);
  printf("ne %d\n",ne);
  e = (int *) malloc(ne*sizeof(int)*4);
  for (i=0; i<ne; i++) {
    fgets(line, MAX_LEN, file);
    sscanf(line,"%d %d %d %d %d",&k,e+4*i,e+4*i+1,e+4*i+2,e+4*i+3);
  }
  for (i=0; i<-ne; i++) {
    printf("%5d%5d%5d%5d%5d%5d\n",i+1,3,
           e[4*i],
           e[4*i+1],
           e[4*i+2],
           e[4*i+3]);
  }

  fgets(line, MAX_LEN, file);
  sscanf(line,"%d",&nf);
  printf("nf %d\n",nf);
  f = (int *) malloc(nf*sizeof(int)*3);
  for (i=0; i<nf; i++) {
    fgets(line, MAX_LEN, file);
    sscanf(line,"%d %d %d %d",&k,f+3*i,f+3*i+1,f+3*i+2);
  }
  for (i=0; i<-nf; i++) {
    printf("%5d%5d%5d%5d%5d\n",i+1,3,
           f[3*i],
           f[3*i+1],
           f[3*i+2]);
  }
  m = 0;
  for (i=0; i<nf*3; i++)
    if (f[i] > m)
      m = f[i];
  printf("m %d\n",m);

  strcpy(name,argv[1]);
  strcat(name,".fem");
  file = fopen(name,"w");
  fprintf(file,"connectivity\n\n");
  for (i=0; i<nf; i++)
    fprintf(file,"%5d%5d%5d%5d%5d\n",i+1,3,
            f[3*i],
            f[3*i+1],
            f[3*i+2]);
  fprintf(file,"coordinates\n%5d%5d\n",3,m);
  for (i=0; i<m; i++)
    fprintf(file,"%5d%8.3f+0%8.3f+0%8.3f+0\n",
            i+1,c[3*i],c[3*i+1],c[3*i+2]);
  fclose(file);

  strcpy(name,argv[1]);
  strcat(name,".feb");
  file = fopen(name,"w");
  fprintf(file,"connectivity\n\n");
  for (i=0; i<ne; i++)
    fprintf(file,"%5d%5d%5d%5d%5d%5d\n",i+1,3,
            e[4*i],
            e[4*i+1],
            e[4*i+2],
            e[4*i+3]);
  fprintf(file,"coordinates\n%5d%5d\n",3,nv);
  for (i=0; i<nv; i++)
    fprintf(file,"%5d%8.3f+0%8.3f+0%8.3f+0\n",
            i+1,c[3*i],c[3*i+1],c[3*i+2]);
  fprintf(file,"end option\nsurface contact\n");
  for (i=0; i<m; i++)
    fprintf(file,"%5d%5d%5d\n",i+1,0,0);
  fprintf(file,"end\n");
  fclose(file);

  return(0);
}
