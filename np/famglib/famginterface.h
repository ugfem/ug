// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      famginterface.h												*/
/*																			*/
/* Purpose:   famginterface													*/
/*																			*/
/* Author:    Christian Wagner												*/
/*			  Institut fuer Computeranwendungen  III						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: chris@ica3.uni-stuttgart.de							*/
/*																			*/
/*																			*/
/* History:   November 97 begin, Stuttgart									*/
/*			  August 98 integration into ug (Christian Wrobel)				*/
/*																			*/
/* Remarks:																	*/
/*																			*/
/****************************************************************************/

#ifndef __CMG_INTERFACE__
#define __CMG_INTERFACE__

/* RCS_ID
   $Header$
 */

#include "gm.h"        /* for ug data structure               */

#define CMG_RHS      0
#define CMG_UNKNOWN  1
#define CMG_DEFECT   2
#define CMG_TVA       3
#define CMG_TVB       4
#define CMG_NVECTORS 5

struct CMG_Interface
{
  int n;
  int nl;
  int nv;
  int *start;
  int *index;
  double *entry;
  double *vector[CMG_NVECTORS];
  void **extra;
};

struct CMG_Parameter
{
  int heap;
  int nv;
  int gamma;
  int n1;
  int n2;
  double ilut;
  double cgilut;
  int cgnodes;
  double mincoarse;
  int conloops;
  int type;
  int stv;
  double tol;
  double sigma;
  double omegar;
  double omegal;
  double error1;
  double error2;
  int maxit;
  double alimit;
  double rlimit;
  double divlimit;
  double reduction;
  char solver[10];
  char presmoother[10];
  char postsmoother[10];
  char cgsmoother[10];
};

struct CMG_IndexBitField
{
  unsigned type : 1;
  unsigned id : 31;
};

struct CMG_MatrixPtr
{
  struct CMG_IndexBitField *index;
  double *entry;
  double **adjoined;
  int nc;
};

typedef struct CMG_MatrixPtr CMG_MatrixPtr;

struct CMG_TransferBitField
{
  unsigned f1 : 1;
  unsigned f0 : 31;
};


struct CMG_TransferEntry
{
  double data;
  struct CMG_TransferBitField id;
  struct CMG_TransferEntry *next;
};

typedef struct CMG_TransferEntry CMG_TransferEntry;

int CMGSolveSystem(struct CMG_Interface*,struct CMG_Parameter*);
void **CMG_GetExtraPtr(int level);
int CMG_GetN(int level);
int CMG_GetNF(int level);
CMG_MatrixPtr  *CMG_GetMatrixPtr(int level,int i);
double * CMG_GetVector(int level, int i);
int CMG_GetMaxLevel();
CMG_TransferEntry  *CMG_GetTransferEntry(int level,int i);

#endif
