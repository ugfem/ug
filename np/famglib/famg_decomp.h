// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      famg_decomp.h													*/
/*																			*/
/* Purpose:   famg decomp classes											*/
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

#ifndef __FAMG_DECOMP__
#define __FAMG_DECOMP__

#include "famg_algebra.h"

/* RCS_ID
   $Header$
 */

struct FAMGDecompBitField
{
  unsigned f0 : 23;
  unsigned f1 : 1;
};

// Class FAMGDecompEntry

class FAMGDecompEntry
{
public:
  double GetData(void) const;
  int GetId(void) const;
  FAMGDecompEntry *GetNext(void) const;
  FAMGDecompEntry *GetReverse();
  FAMGDecompEntry *GetReverse(int);
  void SetData(double);
  void AddData(double);
  void SetId(int);
  void SetRev(int);
  void SetNext(FAMGDecompEntry *);
  void Init(int);
private:
  double data;
  FAMGDecompBitField id;
  FAMGDecompEntry *next;
};


inline double FAMGDecompEntry::GetData() const {
  return data;
}
inline void FAMGDecompEntry::SetData(double val) {
  data = val;
}
inline void FAMGDecompEntry::AddData(double val) {
  data += val;
}
inline int FAMGDecompEntry::GetId() const {
  return id.f0;
}
inline void FAMGDecompEntry::SetId(int i) {
  id.f0 = i;
}
inline void FAMGDecompEntry::SetRev(int i) {
  id.f1 = i;
}
inline FAMGDecompEntry *FAMGDecompEntry::GetNext() const {
  return next;
}
inline void FAMGDecompEntry::SetNext(FAMGDecompEntry *mat) {
  next = mat;
}
inline FAMGDecompEntry *FAMGDecompEntry::GetReverse()
{
  return this+(1 - 2*id.f1);
}


class FAMGDecompRow
{
public:
  double GetData() const;
  int GetId() const;
  void SetData(double);
  void AddData(double);
  void SetId(int);
  FAMGDecompEntry *GetRight() const;
  FAMGDecompEntry *GetLeft() const;
  FAMGDecompEntry* GetEntry(int) const;
  FAMGDecompEntry* NewEntry(FAMGDecompRow*);
  int SaveEntry(FAMGDecompRow*,double);
  int SaveEntry(FAMGDecompRow*,double,FAMGDecompEntry**);
  void SetRight(FAMGDecompEntry*);
  void MoveLeft(FAMGDecompEntry*);
  void Init(int);
private:
  int id;
  double data;
  FAMGDecompEntry *right, *left;
};

inline double FAMGDecompRow::GetData() const {
  return data;
}
inline int FAMGDecompRow::GetId() const {
  return id;
}
inline FAMGDecompEntry *FAMGDecompRow::GetRight() const {
  return right;
}
inline FAMGDecompEntry *FAMGDecompRow::GetLeft() const {
  return left;
}
inline void FAMGDecompRow::SetId(int i) {
  id = i;
}
inline void FAMGDecompRow::SetData(double val) {
  data = val;
}
inline void FAMGDecompRow::AddData(double val) {
  data += val;
}
inline void FAMGDecompRow::SetRight(FAMGDecompEntry *mat) {
  right = mat;
}

class FAMGDecomp
{
public:
  FAMGDecompRow* GetRow() const;
  FAMGDecompRow* GetRow(int) const;
  void SetN(int);
  int Init(FAMGMatrixAlg *);
  int Construct(int);
  void ILUT(double *);
private:
  int n;
  FAMGDecompRow *row;
};

inline FAMGDecompRow *FAMGDecomp::GetRow() const {
  return row;
}
inline FAMGDecompRow *FAMGDecomp::GetRow(int i) const {
  return row+i;
}
inline void FAMGDecomp::SetN(int i) {
  n = i;
}

#endif
