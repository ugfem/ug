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
  int S[6];
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

void PrintCell (CELL C, POINT *P)
{
  printf("id %d n %d subdomain %d property %d\n",
         C.id,
         C.n,
         C.subdomain,
         C.property);
  printf("   points %d %d %d %d %d %d %d %d\n",
         P[C.P[0]].id,
         P[C.P[1]].id,
         P[C.P[2]].id,
         P[C.P[3]].id,
         P[C.P[4]].id,
         P[C.P[5]].id,
         P[C.P[6]].id,
         P[C.P[7]].id);
  printf("   faces %d %d %d %d %d %d\n",
         C.F[0],
         C.F[1],
         C.F[2],
         C.F[3],
         C.F[4],
         C.F[5]);
  printf("   segments %d %d %d %d %d %d\n",
         C.S[0],
         C.S[1],
         C.S[2],
         C.S[3],
         C.S[4],
         C.S[5]);
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
  printf("id %d bnd %d P %f %f %f\n",
         P.id,
         P.bnd,
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
  P->bnd = 0;
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
    C->S[i] = f[i];
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
}

static int FaceCompare (FACE *F, FACE *G)
{
  int *f = F->P;
  int *g = G->P;
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
  int i,n,nP,nF,nC,nBP;
  CELL *C;
  FACE *F;
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

  if (C[nC-1].id == C[nC-2].id) {
    fprintf(stderr,"add newline at the end of %s\n",name);
    return(1);
  }

  /* sort faces */
  qsort(F,nF,sizeof(FACE),
        (int (*)(const void *, const void *))FaceCompare);

  /* install neighborhoods */
  for (i=1; i<nF; i++)
    if (FaceCompare(&(F[i-1]),&(F[i])) == 0) {
      C[F[i-1].C].C[F[i-1].side] = F[i].C;
      C[F[i].C].C[F[i].side] = F[i-1].C;
      i++;
    }
  /* label boundary */
  for (i=0; i<nC; i++)
  {
    int j;

    for (j=0; j<6; j++)
    {
      int k;

      if (C[i].C[j] >= 0) continue;
      for (k=0; k<4; k++)
        P[C[i].P[faces[j][k]]].bnd = 1;
    }
  }
  /* renumber points */
  n = 0;
  for (i=0; i<nP; i++)
    if (P[i].bnd == 1)
      P[i].id = n++;
  nBP = n;
  for (i=0; i<nP; i++)
    if (P[i].bnd == 0)
      P[i].id = n++;

  if (argc>=3) {
    if ((argv[2][0] == '-') && (argv[2][1] == 'v')) {
      for (i=0; i<nF; i++)
        PrintFace(F[i]);
      for (i=0; i<nP; i++)
        if (P[i].bnd == 1)
          PrintPoint(P[i]);
      for (i=0; i<nP; i++)
        if (P[i].bnd == 0)
          PrintPoint(P[i]);
      for (i=0; i<nC; i++)
        PrintCell(C[i],P);
    }
  }
  strcpy(name,argv[1]);
  strcat(name,".fem");
  file = fopen(name,"w");
  fprintf(file,"connectivity\n\n");
  n = 0;
  for (i=0; i<nC; i++)
  {
    int j;

    for (j=0; j<6; j++)
    {
      int k;

      if (C[i].C[j] >= 0) continue;
      fprintf(file,"%5d",1+n++);
      for (k=0; k<4; k++)
        fprintf(file,"%5d",P[C[i].P[faces[j][k]]].id+1);
      fprintf(file,"\n");
    }
  }
  fprintf(file,"coordinates\n%5d%5d\n",3,nBP);
  n = 0;
  for (i=0; i<nP; i++)
    if (P[i].bnd == 1)
      fprintf(file,"%5d%8.2f+0%8.2f+0%8.2f+0\n",
              1+n++,P[i].x[0],P[i].x[1],P[i].x[2]);
  fclose(file);

  strcpy(name,argv[1]);
  strcat(name,".feb");
  file = fopen(name,"w");
  fprintf(file,"connectivity\n\n");
  n=0;
  for (i=0; i<nC; i++)
  {
    int j;

    fprintf(file,"%5d%5d",1+n++,C[i].n-1);
    for (j=0; j<C[i].n; j++)
      fprintf(file,"%5d",P[C[i].P[j]].id+1);
    fprintf(file,"\n");
  }
  fprintf(file,"coordinates\n%5d%5d\n",3,nP);
  n = 0;
  for (i=0; i<nP; i++)
    if (P[i].bnd == 1)
      fprintf(file,"%5d%8.3f+0%8.3f+0%8.3f+0\n",
              1+n++,P[i].x[0],P[i].x[1],P[i].x[2]);
  for (i=0; i<nP; i++)
    if (P[i].bnd == 0)
      fprintf(file,"%5d%8.2f+0%8.2f+0%8.2f+0\n",
              1+n++,P[i].x[0],P[i].x[1],P[i].x[2]);
  fprintf(file,"end option\nsurface contact\n");
  for (i=0; i<nBP; i++)
    fprintf(file,"%5d%5d%5d\n",i+1,0,0);
  fprintf(file,"end\n");
  fclose(file);

  return(0);
}
