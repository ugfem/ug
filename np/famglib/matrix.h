// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      matrix.h														*/
/*																			*/
/* Purpose:   cmg matrix classes											*/
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

#ifndef __CMG_MATRIX__
#define __CMG_MATRIX__

#include "misc.h"
#include "transfer.h"
#include "graph.h"

/* RCS_ID
   $Header$
 */

// for other data structures, CMGMatrixPtr may contain only
// one pointer to the matrix entry

struct CMGIndexBitField
{
  unsigned type : 1;
  unsigned id : 31;
};


class CMGMatrixPtr
{
public:
  void SetNC(int);
  void SetIndexPtr(CMGIndexBitField *);
  void SetIndex(int);
  void SetType(int);
  void SetEntryPtr(double*);
  void SetAdjoined(double*);
  void SetAdjoinedPtr(double**);
  void SetData(double);
  int GetNext();
  double GetData();
  double GetAdjData();
  int GetIndex();
  int GetType();
  double **GetAdjoinedPtr();
  double *GetAdjoined();
  double *GetEntry();
private:
  CMGIndexBitField *index;
  double *entry;
  double **adjoined;
  int nc;         // remaining columns/rows
};

// for other data structures e.g. ptr++ or ptr = ptr->next

inline int CMGMatrixPtr::GetNext()
{
  if( nc <= 0) return 0;
  index++;
  entry++;
  adjoined++;
  nc--;
  return 1;
}

// for other data structures e.g. ptr = ptr->data
inline double CMGMatrixPtr::GetData() {
  return *entry;
}
inline int CMGMatrixPtr::GetIndex() {
  return index->id;
}
inline int CMGMatrixPtr::GetType() {
  return index->type;
}
inline double CMGMatrixPtr::GetAdjData() {
  return **adjoined;
}
inline double *CMGMatrixPtr::GetAdjoined() {
  return *adjoined;
}
inline double **CMGMatrixPtr::GetAdjoinedPtr() {
  return adjoined;
}
inline double *CMGMatrixPtr::GetEntry() {
  return entry;
}
inline void CMGMatrixPtr::SetData(double d) {
  (*entry) = d;
}
inline void CMGMatrixPtr::SetType(int typ) {
  index->type = typ;
}
inline void CMGMatrixPtr::SetIndex(int id) {
  index->id = id;
}
inline void CMGMatrixPtr::SetNC(int i) {
  nc = i;
}
inline void CMGMatrixPtr::SetIndexPtr(CMGIndexBitField *ptr) {
  index = ptr;
}
inline void CMGMatrixPtr::SetEntryPtr(double* ptr) {
  entry = ptr;
}
inline void CMGMatrixPtr::SetAdjoinedPtr(double** ptr) {
  adjoined = ptr;
}
inline void CMGMatrixPtr::SetAdjoined(double* ptr) {
  (*adjoined) = ptr;
}


class CMGMatrix
{
public:
  int GetN() const;
  int GetNL() const;
  double GetDiag(int) const;
  void SetN(int);
  void SetNL(int);
  void SetStartPtr(int *);
  void SetIndex(int*);
  void SetEntry(double*);
  void SetAdjoined(double**);
  int *GetStartPtr() const;
  CMGMatrixPtr GetStart(int i);
  int GetType(int i) const;
  void DevideFGDefect(double *unknown, double *defect);
  void VecMinusMatVec(double *defect, double *rhs, double *unknown);
  void Mult(double *vout, double *vin);
  void MultTrans(double *vout, double *vin);
  int Init(int nn);
  int Init2(int nn);
  int OrderColumns();
  int OrderColumns(int *map);
  int ReorderColumns(int *map);
  int OrderColumns2(CMGGraph*);
  int ConstructAdjoined();
  int ConstructAdjoinedB();
  int ConstructAdjoined2(CMGGraph*);
  int ConstructEnd();
  int ConstructEndB();
  int ConstructEnd2();
  int CGMatrix(CMGMatrix *fgmatrix, CMGTransfer *transfer, int *father);
  int TmpMatrix(CMGMatrix *matrix, CMGTransfer *transfer, CMGGraph *graph);
  void MarkUnknowns(CMGGraph *graph);
  int Order(int *mapping);
  int Reorder(int *mapping);
  int GetSmallestIndex();
  void ModifyIndex(int f);
  void ModifyIndex(int *type, int f);
  void RemodifyIndex(int f);
  void RemodifyIndex(int *type, int f);
  void JAC(double *);
  void FGS(double *);
  void BGS(double *);
  void SGS(double *);
private:
  int n;
  int nl;
  int *start;
  int *end;       // for reorder
  int *index;
  double *entry;
  double **adjoined;
};

inline int CMGMatrix::GetN() const {
  return n;
}
inline int CMGMatrix::GetNL() const {
  return nl;
}
inline int *CMGMatrix::GetStartPtr() const {
  return start;
}
inline double CMGMatrix::GetDiag(int i) const {
  return entry[start[i]];
}
inline int CMGMatrix::GetType(int i) const {
  return ((CMGIndexBitField *)(index+start[i]))->type;
}
inline void CMGMatrix::SetN(int nn) {
  n = nn;
}
inline void CMGMatrix::SetNL(int nn) {
  nl = nn;
}
inline void CMGMatrix::SetStartPtr(int* p) {
  start = p;
}
inline void CMGMatrix::SetIndex(int* p) {
  index = p;
}
inline void CMGMatrix::SetEntry(double* p) {
  entry = p;
}
inline void CMGMatrix::SetAdjoined(double** p) {
  adjoined = p;
}

inline CMGMatrixPtr CMGMatrix::GetStart(int i)
{
  CMGMatrixPtr mat;

  int offset = start[i];
  mat.SetEntryPtr(entry+offset);
  mat.SetIndexPtr((CMGIndexBitField *)(index+offset));
  mat.SetAdjoinedPtr(adjoined+offset);
  mat.SetNC(end[i]-offset-1);

  return mat;
}

#endif
