// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      matrix.h														*/
/*																			*/
/* Purpose:   cmg tranfer classes											*/
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

#ifndef __CMG_TRANSFER__
#define __CMG_TRANSFER__

/* RCS_ID
   $Header$
 */

struct CMGTransferBitField
{
  unsigned f1 : 1;
  unsigned f0 : 31;
};

// Class CMGTransferEntry

class CMGTransferEntry
{
public:
  double GetData(void) const;
  int GetId(void) const;
  CMGTransferEntry *GetNext(void) const;
  CMGTransferEntry *GetReverse();
  CMGTransferEntry *GetReverse(int);
  CMGTransferEntry* GetEntry(int);
  CMGTransferEntry* NewEntry(CMGTransferEntry*);
  int SaveEntry(CMGTransferEntry*,double);
  int SaveEntry(CMGTransferEntry*,double,CMGTransferEntry**);
  void SetData(double);
  void AddData(double);
  void SetId(int);
  void SetRev(int);
  void SetNext(CMGTransferEntry *);
  void Init(int);
private:
  double data;
  CMGTransferBitField id;
  CMGTransferEntry *next;
};


inline double CMGTransferEntry::GetData() const {
  return data;
}
inline void CMGTransferEntry::SetData(double val) {
  data = val;
}
inline void CMGTransferEntry::AddData(double val) {
  data += val;
}
inline int CMGTransferEntry::GetId() const {
  return id.f0;
}
inline void CMGTransferEntry::SetId(int i) {
  id.f0 = i;
}
inline void CMGTransferEntry::SetRev(int i) {
  id.f1 = i;
}
inline CMGTransferEntry *CMGTransferEntry::GetNext() const {
  return next;
}
inline void CMGTransferEntry::SetNext(CMGTransferEntry *mat) {
  next = mat;
}
inline CMGTransferEntry *CMGTransferEntry::GetReverse()
{
  return this+(1 - 2*id.f1);
}




class CMGTransfer
{
public:
  CMGTransferEntry* GetRow() const;
  CMGTransferEntry* GetRow(int) const;
  void SetRow(CMGTransferEntry*);
  int Init(class CMGGrid*);
  int Construct(class CMGGrid*);
  int Reorder(int *);
  int Order(int *);
private:
  int n;
  CMGTransferEntry *row;
};

inline CMGTransferEntry *CMGTransfer::GetRow() const {
  return row;
}
inline CMGTransferEntry *CMGTransfer::GetRow(int i) const {
  return row+i;
}
inline void CMGTransfer::SetRow(CMGTransferEntry *ptr) {
  row = ptr;
}



#endif
