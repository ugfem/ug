// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      matrix.h														*/
/*																			*/
/* Purpose:   famg matrix classes											*/
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

#ifndef __FAMG_MATRIX__
#define __FAMG_MATRIX__

#include "misc.h"
#include "transfer.h"
#include "graph.h"

/* RCS_ID
   $Header$
 */

// for other data structures, FAMGMatrixPtr may contain only
// one pointer to the matrix entry

struct FAMGIndexBitField
{
  unsigned type : 1;
  unsigned id : 31;
};


class FAMGMatrixPtr
{
public:
  void SetNC(int);
  void SetIndexPtr(FAMGIndexBitField *);
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
  FAMGIndexBitField *index;
  double *entry;
  double **adjoined;
  int nc;         // remaining columns/rows
};

// for other data structures e.g. ptr++ or ptr = ptr->next

inline int FAMGMatrixPtr::GetNext()
{
  if( nc <= 0) return 0;
  index++;
  entry++;
  adjoined++;
  nc--;
  return 1;
}

// for other data structures e.g. ptr = ptr->data
inline double FAMGMatrixPtr::GetData() {
  return *entry;
}
inline int FAMGMatrixPtr::GetIndex() {
  return index->id;
}
inline int FAMGMatrixPtr::GetType() {
  return index->type;
}
inline double FAMGMatrixPtr::GetAdjData() {
  return **adjoined;
}
inline double *FAMGMatrixPtr::GetAdjoined() {
  return *adjoined;
}
inline double **FAMGMatrixPtr::GetAdjoinedPtr() {
  return adjoined;
}
inline double *FAMGMatrixPtr::GetEntry() {
  return entry;
}
inline void FAMGMatrixPtr::SetData(double d) {
  (*entry) = d;
}
inline void FAMGMatrixPtr::SetType(int typ) {
  index->type = typ;
}
inline void FAMGMatrixPtr::SetIndex(int id) {
  index->id = id;
}
inline void FAMGMatrixPtr::SetNC(int i) {
  nc = i;
}
inline void FAMGMatrixPtr::SetIndexPtr(FAMGIndexBitField *ptr) {
  index = ptr;
}
inline void FAMGMatrixPtr::SetEntryPtr(double* ptr) {
  entry = ptr;
}
inline void FAMGMatrixPtr::SetAdjoinedPtr(double** ptr) {
  adjoined = ptr;
}
inline void FAMGMatrixPtr::SetAdjoined(double* ptr) {
  (*adjoined) = ptr;
}


class FAMGMatrix
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
  FAMGMatrixPtr GetStart(int i);
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
  int OrderColumns2(FAMGGraph*);
  int ConstructAdjoined();
  int ConstructAdjoinedB();
  int ConstructAdjoined2(FAMGGraph*);
  int ConstructEnd();
  int ConstructEndB();
  int ConstructEnd2();
  int CGMatrix(FAMGMatrix *fgmatrix, FAMGTransfer *transfer, int *father);
  int TmpMatrix(FAMGMatrix *matrix, FAMGTransfer *transfer, FAMGGraph *graph);
  void MarkUnknowns(FAMGGraph *graph);
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

inline int FAMGMatrix::GetN() const {
  return n;
}
inline int FAMGMatrix::GetNL() const {
  return nl;
}
inline int *FAMGMatrix::GetStartPtr() const {
  return start;
}
inline double FAMGMatrix::GetDiag(int i) const {
  return entry[start[i]];
}
inline int FAMGMatrix::GetType(int i) const {
  return ((FAMGIndexBitField *)(index+start[i]))->type;
}
inline void FAMGMatrix::SetN(int nn) {
  n = nn;
}
inline void FAMGMatrix::SetNL(int nn) {
  nl = nn;
}
inline void FAMGMatrix::SetStartPtr(int* p) {
  start = p;
}
inline void FAMGMatrix::SetIndex(int* p) {
  index = p;
}
inline void FAMGMatrix::SetEntry(double* p) {
  entry = p;
}
inline void FAMGMatrix::SetAdjoined(double** p) {
  adjoined = p;
}

inline FAMGMatrixPtr FAMGMatrix::GetStart(int i)
{
  FAMGMatrixPtr mat;

  int offset = start[i];
  mat.SetEntryPtr(entry+offset);
  mat.SetIndexPtr((FAMGIndexBitField *)(index+offset));
  mat.SetAdjoinedPtr(adjoined+offset);
  mat.SetNC(end[i]-offset-1);

  return mat;
}

#endif
