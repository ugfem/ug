// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugmarc.c														*/
/*																			*/
/* Purpose:   converter into MARC format                                                                */
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   Apr 28, 1999                                                                              */
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

#define ModelP

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

#include "ugdevices.h"

#include "gm.h"
#include "algebra.h"
#include "misc.h"
#include "ugm.h"
#include "ugio.h"
#include "elements.h"
#include "shapes.h"
#include "mgio.h"
#include "dio.h"
#include "parallel.h"


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

static int file_readline (FILE *f, char *key)
{
  int len = strlen(key);
  char line[MAX_LEN];

  do {
    fgets(line, MAX_LEN, f);
    /* printf("%s",line); */
    if (strncmp(line,key,len) == 0) return(0);
  } while (!feof(f));

  return(1);
}

int main (int argc, char **argv)
{
  FILE *file;
  char line[MAX_LEN],name[32];
  int i,j,k,nv,nb,nl,nf,ne,m;
  int *v,*iv,*b,*l,*f,*fc,*e,*ec,*t;
  double *c;

  if (argc<2)     {
    printf("filename required\n");
    return (1);
  }
  strcpy(name,argv[1]);
  strcat(name,".net");
  file = fopen(name,"r");
  if (file == NULL) {
    printf("cannot open file %s\n",name);
    return (1);
  }
  fgets(line, MAX_LEN, file);
  fgets(line, MAX_LEN, file);
  sscanf(line+3,"VertexNumber: %d",&nv);
  fgets(line, MAX_LEN, file);
  sscanf(line+3,"EdgeNumber: %d",&nl);
  fgets(line, MAX_LEN, file);
  sscanf(line+3,"FaceNumber: %d",&nf);
  fgets(line, MAX_LEN, file);
  sscanf(line+3,"ElementNumber: %d",&ne);

  printf("%s: %d vertices %d edges %d faces %d elements\n",
         argv[1],nv,nl,nf,ne);

  v = (int *) malloc(nv*sizeof(int));
  iv = (int *) malloc(nv*sizeof(int));
  b = (int *) malloc(nv*sizeof(int));
  c = (double *) malloc(nv*sizeof(double)*3);
  for (i=0; i<nv; i++)
    b[i] = 0;
  l = (int *) malloc(nl*sizeof(int)*2);
  f = (int *) malloc(nf*sizeof(int)*3);
  fc = (int *) malloc(nf*sizeof(int)*3);
  e = (int *) malloc(ne*sizeof(int)*4);
  ec = (int *) malloc(ne*sizeof(int)*4);
  t = (int *) malloc(ne*sizeof(int)*4);

  file_readline(file,"% Vert");
  for (i=0; i<nv; i++) {
    fgets(line, MAX_LEN, file);
    sscanf(line,"%lf %lf %lf",&(c[3*i]),&(c[3*i+1]),&(c[3*i+2]));
    /* printf("%s: %f %f %f\n",line,c[3*i],c[3*i+1],c[3*i+2]); */
  }
  file_readline(file,"% Edge");
  for (i=0; i<nl; i++) {
    fgets(line, MAX_LEN, file);
    sscanf(line,"%d: %d %d",&k,&(l[2*i]),&(l[2*i+1]));
    /*printf("%s: %d %d\n",line,l[2*i],l[2*i+1]); */
    if (k > 0) b[l[2*i]] = b[l[2*i+1]] = 1;
  }
  file_readline(file,"% Face");
  for (i=0; i<nf; i++) {
    fgets(line, MAX_LEN, file);
    sscanf(line,"%d: %d %d %d",&(fc[i]),&(f[3*i]),&(f[3*i+1]),&(f[3*i+2]));
    /* printf("%s: %d %d\n",line,f[3*i],f[3*i+1],f[3*i+2]); */
  }
  file_readline(file,"% Element");
  for (i=0; i<ne; i++) {
    fgets(line, MAX_LEN, file);
    sscanf(line,"%d: %d %d %d %d",
           &(ec[i]),&(e[4*i]),&(e[4*i+1]),&(e[4*i+2]),&(e[4*i+3]));
    /* printf("%s: %d %d\n",line,e[4*i],e[4*i+1]); */

    t[4*i] = t[4*i+1] = t[4*i+2] = t[4*i+3] = -1;
    for (j=0; j<4; j++)
      for (k=0; k<3; k++)
        for (m=0; m<2; m++) {
          if (l[f[e[4*i+j]*3+k]*2+m] == t[4*i]) continue;
          if (l[f[e[4*i+j]*3+k]*2+m] == t[4*i+1]) continue;
          if (l[f[e[4*i+j]*3+k]*2+m] == t[4*i+2]) continue;
          if (l[f[e[4*i+j]*3+k]*2+m] == t[4*i+3]) continue;
          if (t[4*i] == -1) {
            t[4*i] = l[f[e[4*i+j]*3+k]*2+m];
            continue;
          }
          if (t[4*i+1] == -1) {
            t[4*i+1] = l[f[e[4*i+j]*3+k]*2+m];
            continue;
          }
          if (t[4*i+2] == -1) {
            t[4*i+2] = l[f[e[4*i+j]*3+k]*2+m];
            continue;
          }
          if (t[4*i+3] == -1) {
            t[4*i+3] = l[f[e[4*i+j]*3+k]*2+m];
            continue;
          }
        }
    /* printf("e%d: %d %d %d %d\n",i,t[4*i],t[4*i+1],t[4*i+2],t[4*i+3]); */
  }
  fclose(file);

  j=0;
  for (i=0; i<nv; i++)
    if (b[i] == 1) {
      iv[i] = j;
      v[j++] = i;
    }
  nb = j;
  for (i=0; i<nv; i++)
    if (b[i] == 0)  {
      iv[i] = j;
      v[j++] = i;
    }
  for (i=0; i<nv; i++) {
    /* printf("%d(%d,%d): %f %f %f\n",v[i],iv[i],b[i],c[3*i],c[3*i+1],c[3*i+2]);*/
  }
  strcpy(name,argv[1]);
  strcat(name,".fem");
  file = fopen(name,"w");
  fprintf(file,"connectivity\n\n");
  k=1;
  for (i=0; i<nf; i++) {
    if (fc[i])
      fprintf(file,"%5d%5d%5d%5d%5d\n",k++,fc[i],
              v[l[2*f[3*i]]]+1,
              v[l[2*f[3*i]+1]]+1,
              v[l[2*f[3*i+1]+1]]+1);
  }
  fprintf(file,"coordinates\n%5d%5d\n",3,nb);
  for (i=0; i<nb; i++) {
    k = iv[i];
    fprintf(file,"%5d%8.5f+0%8.4f+0%8.4f+0\n",
            i+1,c[3*k],c[3*k+1],c[3*k+2]);
  }
  fclose(file);

  strcpy(name,argv[1]);
  strcat(name,".feb");
  file = fopen(name,"w");
  fprintf(file,"connectivity\n\n");
  for (i=0; i<ne; i++) {
    fprintf(file,"%5d%5d%5d%5d%5d%5d\n",i+1,3,
            v[t[4*i]]+1,
            v[t[4*i+1]]+1,
            v[t[4*i+2]]+1,
            v[t[4*i+3]]+1);
  }
  fprintf(file,"coordinates\n%5d%5d\n",3,nv);
  for (i=0; i<nv; i++) {
    k = iv[i];
    fprintf(file,"%5d%8.5f+0%8.4f+0%8.4f+0\n",
            i+1,c[3*k],c[3*k+1],c[3*k+2]);
  }
  fprintf(file,"end option\nsurface contact\n");
  for (i=0; i<nb; i++)
    fprintf(file,"%5d%5d%5d\n",i+1,0,0);
  fprintf(file,"end\n");
  fclose(file);

  return(0);
}
