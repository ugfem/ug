// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      famg_matrix.h													*/
/*																			*/
/* Purpose:   famg tranfer classes											*/
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

#ifndef __FAMG_TRANSFER__
#define __FAMG_TRANSFER__

/* RCS_ID
   $Header$
 */

struct FAMGTransferBitField
{
  unsigned f1 : 1;
  unsigned f0 : 31;
};

// Class FAMGTransferEntry

class FAMGTransferEntry
{
public:
  double GetData(void) const;
  int GetId(void) const;
  FAMGTransferEntry *GetNext(void) const;
  FAMGTransferEntry *GetReverse();
  FAMGTransferEntry *GetReverse(int);
  FAMGTransferEntry* GetEntry(int);
  FAMGTransferEntry* NewEntry(FAMGTransferEntry*);
  int SaveEntry(FAMGTransferEntry*,double);
  int SaveEntry(FAMGTransferEntry*,double,FAMGTransferEntry**);
  void SetData(double);
  void AddData(double);
  void SetId(int);
  void SetRev(int);
  void SetNext(FAMGTransferEntry *);
  void Init(int);
private:
  double data;
  FAMGTransferBitField id;
  FAMGTransferEntry *next;
};


inline double FAMGTransferEntry::GetData() const {
  return data;
}
inline void FAMGTransferEntry::SetData(double val) {
  data = val;
}
inline void FAMGTransferEntry::AddData(double val) {
  data += val;
}
inline int FAMGTransferEntry::GetId() const {
  return id.f0;
}
inline void FAMGTransferEntry::SetId(int i) {
  id.f0 = i;
}
inline void FAMGTransferEntry::SetRev(int i) {
  id.f1 = i;
}
inline FAMGTransferEntry *FAMGTransferEntry::GetNext() const {
  return next;
}
inline void FAMGTransferEntry::SetNext(FAMGTransferEntry *mat) {
  next = mat;
}
inline FAMGTransferEntry *FAMGTransferEntry::GetReverse()
{
  return this+(1 - 2*id.f1);
}




class FAMGTransfer
{
public:
  FAMGTransferEntry* GetRow() const;
  FAMGTransferEntry* GetRow(int) const;
  void SetRow(FAMGTransferEntry*);
  int Init(class FAMGGrid*);
  int Construct(class FAMGGrid*);
  int Reorder(int *);
  int Order(int *);
private:
  int n;
  FAMGTransferEntry *row;
};

inline FAMGTransferEntry *FAMGTransfer::GetRow() const {
  return row;
}
inline FAMGTransferEntry *FAMGTransfer::GetRow(int i) const {
  return row+i;
}
inline void FAMGTransfer::SetRow(FAMGTransferEntry *ptr) {
  row = ptr;
}



#endif
