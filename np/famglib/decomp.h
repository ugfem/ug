// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      decomp.h														*/
/*																			*/
/* Purpose:   cmg decomp classes											*/
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

#ifndef __CMG_DECOMP__
#define __CMG_DECOMP__

#include "matrix.h"

/* RCS_ID
   $Header$
 */

struct CMGDecompBitField
{
  unsigned f0 : 23;
  unsigned f1 : 1;
};

// Class CMGDecompEntry

class CMGDecompEntry
{
public:
  double GetData(void) const;
  int GetId(void) const;
  CMGDecompEntry *GetNext(void) const;
  CMGDecompEntry *GetReverse();
  CMGDecompEntry *GetReverse(int);
  void SetData(double);
  void AddData(double);
  void SetId(int);
  void SetRev(int);
  void SetNext(CMGDecompEntry *);
  void Init(int);
private:
  double data;
  CMGDecompBitField id;
  CMGDecompEntry *next;
};


inline double CMGDecompEntry::GetData() const {
  return data;
}
inline void CMGDecompEntry::SetData(double val) {
  data = val;
}
inline void CMGDecompEntry::AddData(double val) {
  data += val;
}
inline int CMGDecompEntry::GetId() const {
  return id.f0;
}
inline void CMGDecompEntry::SetId(int i) {
  id.f0 = i;
}
inline void CMGDecompEntry::SetRev(int i) {
  id.f1 = i;
}
inline CMGDecompEntry *CMGDecompEntry::GetNext() const {
  return next;
}
inline void CMGDecompEntry::SetNext(CMGDecompEntry *mat) {
  next = mat;
}
inline CMGDecompEntry *CMGDecompEntry::GetReverse()
{
  return this+(1 - 2*id.f1);
}


class CMGDecompRow
{
public:
  double GetData() const;
  int GetId() const;
  void SetData(double);
  void AddData(double);
  void SetId(int);
  CMGDecompEntry *GetRight() const;
  CMGDecompEntry *GetLeft() const;
  CMGDecompEntry* GetEntry(int) const;
  CMGDecompEntry* NewEntry(CMGDecompRow*);
  int SaveEntry(CMGDecompRow*,double);
  int SaveEntry(CMGDecompRow*,double,CMGDecompEntry**);
  void SetRight(CMGDecompEntry*);
  void MoveLeft(CMGDecompEntry*);
  void Init(int);
private:
  int id;
  double data;
  CMGDecompEntry *right, *left;
};

inline double CMGDecompRow::GetData() const {
  return data;
}
inline int CMGDecompRow::GetId() const {
  return id;
}
inline CMGDecompEntry *CMGDecompRow::GetRight() const {
  return right;
}
inline CMGDecompEntry *CMGDecompRow::GetLeft() const {
  return left;
}
inline void CMGDecompRow::SetId(int i) {
  id = i;
}
inline void CMGDecompRow::SetData(double val) {
  data = val;
}
inline void CMGDecompRow::AddData(double val) {
  data += val;
}
inline void CMGDecompRow::SetRight(CMGDecompEntry *mat) {
  right = mat;
}

class CMGDecomp
{
public:
  CMGDecompRow* GetRow() const;
  CMGDecompRow* GetRow(int) const;
  void SetN(int);
  int Init(CMGMatrix *);
  int Construct(int);
  void ILUT(double *);
private:
  int n;
  CMGDecompRow *row;
};

inline CMGDecompRow *CMGDecomp::GetRow() const {
  return row;
}
inline CMGDecompRow *CMGDecomp::GetRow(int i) const {
  return row+i;
}
inline void CMGDecomp::SetN(int i) {
  n = i;
}

#endif
