// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugmesh.c                                                                                                      */
/*																			*/
/* Purpose:   converter from mesh into UG format                                                */
/*																			*/
/* Author:	  Christian Wieners                                                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   Dec 1, 1999                                                                               */
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

#define MAX_LEN      128
#define MAX_LINES     16
#define MAX_SEGMENTS   8

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

/* UG_RCS_STRING
   $Header$
 */

typedef struct {

  int id;
  int n;
  int subdomain;
  int property;
  int P[8];
  int F[6];
  int C[6];

} CELL;

typedef struct {

  int id;
  int n;
  int C;
  int side;
  int S;
  int P[4];

} FACE;

typedef struct {

  int id;
  int bnd;
  double x[3];

} POINT;

typedef struct {

  int id;
  int L[MAX_LINES];

} CORNER;

typedef struct {

  int id;
  int S[MAX_SEGMENTS];

} LINE;

void PrintCell (CELL C, FACE *F)
{
  printf("id %d n %d subdomain %d property %d\n",
         C.id,
         C.n,
         C.subdomain,
         C.property);
  printf("   points %d %d %d %d %d %d %d %d\n",
         C.P[0],
         C.P[1],
         C.P[2],
         C.P[3],
         C.P[4],
         C.P[5],
         C.P[6],
         C.P[7]);
  printf("   faces %d %d %d %d %d %d\n",
         C.F[0],
         C.F[1],
         C.F[2],
         C.F[3],
         C.F[4],
         C.F[5]);
  if (F != NULL)
    printf("   segments %d %d %d %d %d %d\n",
           F[C.F[0]].S,
           F[C.F[1]].S,
           F[C.F[2]].S,
           F[C.F[3]].S,
           F[C.F[4]].S,
           F[C.F[5]].S);
  printf("   cells %d %d %d %d %d %d\n",
         C.C[0],
         C.C[1],
         C.C[2],
         C.C[3],
         C.C[4],
         C.C[5]);
}

void PrintFace (FACE F)
{
  printf("id %d n %d cell %d side %d segment %d\n",
         F.id,
         F.n,
         F.C,
         F.side,
         F.S);
  printf("   points %d %d %d %d\n",
         F.P[0],
         F.P[1],
         F.P[2],
         F.P[3]);
}

void PrintPoint (POINT P)
{
  printf("id %d P %f %f %f\n",
         P.id,
         P.x[0],
         P.x[1],
         P.x[2]);
}

void ReadPoint (char *line, POINT *P)
{
  int id;
  double x,y,z;
  int k = sscanf(line,"c %d %lf %lf %lf",&id,&x,&y,&z);
  assert(k == 4);
  P->id = id-1;
  P->x[0] = x;
  P->x[1] = y;
  P->x[2] = z;
}

static int faces[6][4] = {
  {0,3,2,1},{0,1,4,5},{1,2,5,6},{2,3,6,7},{3,0,7,4},{4,5,6,7}
};

void ReadCell (char *line, CELL *C, int *nF, FACE *F)
{
  int id,p,p0,p1,p2,p3,p4,p5,p6,p7,f[6],i;
  int k = sscanf(line,
                 "x %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                 &id,&p,
                 &p0,&p1,&p2,&p3,&p4,&p5,&p6,&p7,
                 &(f[0]),&(f[1]),&(f[2]),&(f[3]),&(f[4]),&(f[5]));
  assert(k == 16);
  C->id = id-1;
  C->property = p;
  C->n = 8;
  C->subdomain = 1;
  C->P[0] = p0-1;
  C->P[1] = p1-1;
  C->P[2] = p2-1;
  C->P[3] = p3-1;
  C->P[4] = p4-1;
  C->P[5] = p5-1;
  C->P[6] = p6-1;
  C->P[7] = p7-1;
  for (i=0; i<6; i++)
  {
    int j;

    F[*nF].id = *nF;
    F[*nF].side = i;
    F[*nF].n = 4;
    F[*nF].C = id-1;
    F[*nF].S = f[i];
    for (j=0; j<4; j++) {
      F[*nF].P[j] = C->P[faces[i][j]];
    }
    C->F[i] = (*nF)++;
    C->C[i] = -1-f[i];
  }
}

static void Sort (int *f, int *ff)
{
  int i;

  ff[0] = f[0];
  for (i=1; i<4; i++)
  {
    int j;

    for (j=0; j<i; j++)
      if (f[i] < ff[j])
      {
        int k;

        for (k=i; k>j; k--)
          ff[k] = ff[k-1];
        ff[k] = f[i];
        break;
      }
    ff[j] = f[i];
  }

  printf("f  %d %d %d %d\nff %d %d %d %d\n",
         f[0],
         f[1],
         f[2],
         f[3],
         ff[0],
         ff[1],
         ff[2],
         ff[3]);
}

static int FaceCompare (FACE **F, FACE **G)
{
  int *f = F[0]->P;
  int *g = G[0]->P;
  int ff[4],gg[4];
  int i;

  Sort(f,ff);
  Sort(g,gg);

  for (i=0; i<4; i++)
    if (ff[i] < gg[i]) return(1);
    else if (ff[i] > gg[i]) return(-1);

  return(0);
}

int main (int argc, char **argv)
{
  FILE *file;
  char line[MAX_LEN],name[MAX_LEN];
  int i,nP,nF,nC;
  CELL *C;
  FACE *F,**FF;
  POINT *P;

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
  nC = 0;
  nP = 0;
  while (!feof(file)) {
    fgets(line, MAX_LEN, file);
    if (line[0] == 'x') nC++;
    if (line[0] == 'c') nP++;
  }
  fclose(file);

  P = (POINT *) malloc(nP*sizeof(POINT));
  C = (CELL *) malloc(nC*sizeof(CELL));
  F = (FACE *) malloc(6*nC*sizeof(FACE));
  FF = (FACE **) malloc(6*nC*sizeof(FACE *));

  file = fopen(name,"r");
  nC = 0;
  nP = 0;
  nF = 0;
  while (!feof(file)) {
    fgets(line, MAX_LEN, file);
    if (line[0] == 'x') ReadCell(line,C+(nC++),&nF,F);
    if (line[0] == 'c') ReadPoint(line,P+(nP++));
  }
  fclose(file);

  for (i=0; i<nF; i++)
    FF[i] = F+i;

  qsort(FF,nF,sizeof(FACE *),
        (int (*)(const void *, const void *))FaceCompare);

  for (i=1; i<nF; i++) {

    PrintFace(FF[i-1][0]);
    PrintFace(FF[i][0]);

    printf("comp %d\n",FaceCompare(&(FF[i-1]),&(FF[i])));

    if (FaceCompare(&(FF[i-1]),&(FF[i])) == 0) {
      C[FF[i-1]->C].C[FF[i-1]->side] = FF[i]->C;
      C[FF[i]->C].C[FF[i]->side] = FF[i-1]->C;
      i++;
    }
  }
  for (i=0; i<nF; i++)
    PrintFace(FF[i][0]);
  for (i=0; i<nP; i++)
    PrintPoint(P[i]);
  for (i=0; i<nC; i++)
    PrintCell(C[i],F);


  /*



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
      fprintf(file,"%5d%8.5f+0%8.4f+0%8.4f+0\n",
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
      fprintf(file,"%5d%8.5f+0%8.4f+0%8.4f+0\n",
                          i+1,c[3*i],c[3*i+1],c[3*i+2]);
     fprintf(file,"end option\nsurface contact\n");
     for (i=0; i<m; i++)
      fprintf(file,"%5d%5d%5d\n",i+1,0,0);
     fprintf(file,"end\n");
     fclose(file);
   */

  return(0);
}
