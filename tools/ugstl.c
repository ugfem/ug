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
#include <math.h>
#include <assert.h>

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
#define MAX_FACES   100000
#define MAX_CORNERS 100000

#define V3_SUBTRACT(A,B,C)                         {(C)[0] = (A)[0] - (B)[0];\
                                                    (C)[1] = (A)[1] - (B)[1];\
                                                    (C)[2] = (A)[2] - (B)[2];}

#define V3_SCALE(c,C)                              {(C)[0] = (c)*(C)[0];\
                                                    (C)[1] = (c)*(C)[1];\
                                                    (C)[2] = (c)*(C)[2];}

#define V3_VECTOR_PRODUCT(A,B,C)           {(C)[0] = (A)[1]*(B)[2] - (A)[2]*(B)[1];\
                                            (C)[1] = (A)[2]*(B)[0] - (A)[0]*(B)[2];\
                                            (C)[2] = (A)[0]*(B)[1] - (A)[1]*(B)[0];}

#define V3_EUKLIDNORM(A,b)                              (b) = (sqrt((double)((A)[0]*(A)[0]+(A)[1]*(A)[1]+(A)[2]*(A)[2])));

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

int Marc_Extended = 0;

static int ExpandLine (char *theLine)
{
  int i,j,k,l;

  if (Marc_Extended) {
    i = 76;
    j = 69;
    for (l=0; l<3; l++) {
      theLine[i] = ' ';
      i--;
      theLine[i] = theLine[j];
      i--;
      j--;
      theLine[i] = theLine[j];
      i--;
      j--;
      theLine[i] = 'e';
      i--;
      for (k=0; k<18; k++) {
        theLine[i] = theLine[j];
        i--;
        j--;
      }
    }
    theLine[i] = ' ';
    /*printf("%s",theLine); */
  }
  else {
    i = 41;
    j = 34;
    for (l=0; l<3; l++) {
      theLine[i] = ' ';
      i--;
      theLine[i] = theLine[j];
      i--;
      j--;
      theLine[i] = theLine[j];
      i--;
      j--;
      theLine[i] = 'e';
      i--;
      for (k=0; k<8; k++) {
        theLine[i] = theLine[j];
        i--;
        j--;
      }
    }
    theLine[i] = ' ';
  }

  for (i=5; i<strlen(theLine); i++)
    if (theLine[i] == '-')
      if (theLine[i+1] == 'e') {
        theLine[i] = 'e';
        theLine[i+1] = '-';
      }

  return(0);
}

static int file_readline (FILE *f, char *key)
{
  char theLine[MAX_LEN];

  int len = strlen(key);

  do {
    fgets(theLine, MAX_LEN, f);
    if (strncmp(theLine,key,len) == 0) return(0);
  } while (!feof(f));

  return(1);
}

int PrintNormal (FILE *f, double *a, double *b, double *c)
{
  double d[3],e[3],n[3];
  double l;

  V3_SUBTRACT(b,a,d);
  V3_SUBTRACT(c,a,e);
  V3_VECTOR_PRODUCT(d,e,n);
  V3_EUKLIDNORM(n,l);
  if (l==0.0)
    return(1);
  l = 1.0 / l;
  V3_SCALE(l,n);

  fprintf(f,"  facet normal %8.4fe+000 %8.4fe+000 %8.4fe+000\n",
          n[0],
          n[1],
          n[2]);

  return(0);
}

int main (int argc, char **argv)
{
  FILE *file;
  char line[MAX_LEN],name[MAX_LEN];
  int i,j,k,nv,nf,nc;
  int f[4*MAX_FACES];
  int n[MAX_FACES];
  double c[MAX_CORNERS];

  if (argc<2)     {
    printf("filename required\n");
    return (1);
  }
  strcpy(name,argv[1]);
  strcat(name,".fem");
  file = fopen(name,"r");

  if (file_readline(file,"extended")) {
    Marc_Extended = 0;
    fclose(file);
    file = fopen(name,"r");
  }
  else {
    Marc_Extended = 1;
  }

  if (file_readline(file,"connectivity")) {
    printf("wrong file format\n");
    return(1);
  }
  nf = 0;
  do {
    fgets(line, MAX_LEN, file);
    if (strlen(line) < 3) continue;
    n[nf] = sscanf(line,"%5d%5d%5d%5d%5d%5d",&i,&j,
                   &(f[4*nf]),
                   &(f[4*nf+1]),
                   &(f[4*nf+2]),
                   &(f[4*nf+3])) - 2;
    if (n[nf] <3) break;
    /* printf("%d %d %d %d %d\n",i,nf,f[4*nf],f[4*nf+1],f[4*nf+2]); */
    nf++;
    if (nf > MAX_FACES) {
      printf("increase MAX_FACES\n");
      return(1);
    }
  } while (!feof(file));
  fgets(line, MAX_LEN, file);
  nc = 0;
  do {
    fgets(line, MAX_LEN, file);
    /*printf("%s",line); */
    ExpandLine(line);
    if (sscanf(line,"%d",&k) != 1)
      break;
    if (sscanf(line,"%d %lg %lg %lg",&i,
               &(c[3*k]),
               &(c[3*k+1]),
               &(c[3*k+2])) != 4)
      break;
    nc++;
    /* printf("%d %d %f %f %f\n",k,nc,c[3*k],c[3*k+1],c[3*k+2]);*/
    if (k >= MAX_CORNERS) {
      printf("increase MAX_CORNERS\n");
      return(1);
    }
  } while (!feof(file));
  printf("nf %d nc %d\n",nf,nc);

  strcpy(name,argv[1]);
  strcat(name,".stl");
  file = fopen(name,"w");

  fprintf(file,"solid\n");
  for (i=0; i<nf; i++)
  {
    if (PrintNormal(file,&(c[3*f[4*i]]),&(c[3*f[4*i+1]]),&(c[3*f[4*i+2]])))
    {
      printf("%d %d %d %d\n",i,f[4*i],f[4*i+1],f[4*i+2]);
      printf("c %d %f %f %f\n",f[4*i],
             c[3*f[4*i]],
             c[3*f[4*i]+1],
             c[3*f[4*i]+2]);
      printf("c %d %f %f %f\n",f[4*i+1],
             c[3*f[4*i+1]],
             c[3*f[4*i+1]+1],
             c[3*f[4*i+1]+2]);
      printf("c %d %f %f %f\n",f[4*i+2],
             c[3*f[4*i+2]],
             c[3*f[4*i+2]+1],
             c[3*f[4*i+2]+2]);
      return(1);
    }
    fprintf(file,"    outer loop\n");
    fprintf(file,"      vertex %8.4fe+000 %8.4fe+000 %8.4fe+000\n",
            c[3*f[4*i]],
            c[3*f[4*i]+1],
            c[3*f[4*i]+2]);
    fprintf(file,"      vertex %8.4fe+000 %8.4fe+000 %8.4fe+000\n",
            c[3*f[4*i+1]],
            c[3*f[4*i+1]+1],
            c[3*f[4*i+1]+2]);
    fprintf(file,"      vertex %8.4fe+000 %8.4fe+000 %8.4fe+000\n",
            c[3*f[4*i+2]],
            c[3*f[4*i+2]+1],
            c[3*f[4*i+2]+2]);
    fprintf(file,"    endloop\n");
    fprintf(file,"  endfacet\n");

    if (n[i] == 3) continue;

    if (PrintNormal(file,&(c[3*f[4*i]]),&(c[3*f[4*i+2]]),&(c[3*f[4*i+3]])))
    {
      printf("%d %d %d %d\n",i,f[4*i],f[4*i+2],f[4*i+3]);
      printf("c %d %f %f %f\n",f[4*i],
             c[3*f[4*i]],
             c[3*f[4*i]+2],
             c[3*f[4*i]+3]);
      printf("c %d %f %f %f\n",f[4*i+1],
             c[3*f[4*i+1]],
             c[3*f[4*i+1]+2],
             c[3*f[4*i+1]+3]);
      printf("c %d %f %f %f\n",f[4*i+2],
             c[3*f[4*i+2]],
             c[3*f[4*i+2]+2],
             c[3*f[4*i+2]+3]);
      return(1);
    }
    fprintf(file,"    outer loop\n");
    fprintf(file,"      vertex %8.4fe+000 %8.4fe+000 %8.4fe+000\n",
            c[3*f[4*i]],
            c[3*f[4*i]+1],
            c[3*f[4*i]+2]);
    fprintf(file,"      vertex %8.4fe+000 %8.4fe+000 %8.4fe+000\n",
            c[3*f[4*i+2]],
            c[3*f[4*i+2]+1],
            c[3*f[4*i+2]+2]);
    fprintf(file,"      vertex %8.4fe+000 %8.4fe+000 %8.4fe+000\n",
            c[3*f[4*i+3]],
            c[3*f[4*i+3]+1],
            c[3*f[4*i+3]+2]);
    fprintf(file,"    endloop\n");
    fprintf(file,"  endfacet\n");
  }
  fprintf(file,"endsolid\n");

  fclose(file);

  return(0);
}
