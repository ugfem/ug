// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      famg_transfer.h												*/
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

#ifdef USE_UG_DS
extern "C"
{
#include "gm.h"
}
#include "famg_ugalgebra.h"

#else

#include "famg_arrayalgebra.h"
#endif

#include "famg_misc.h"

// forward declaration
class FAMGGrid;

/* RCS_ID
   $Header$
 */

// Class FAMGTransferEntry

class FAMGTransferEntry
{
public:
  FAMGVectorEntry GetCol(void) const;
  FAMGTransferEntry *GetNext(void) const;
  void SetNext(FAMGTransferEntry *);
  FAMGTransferEntry *GetEntry(const FAMGVectorEntry &cg_vec);
  double GetProlongation(void) const;
  double GetRestriction(void) const;
  void SetProlongation(double);
  void SetRestriction(double);
  //void AddProlongation(double);
  //void AddRestriction(double);
#ifdef USE_UG_DS
  void GetColInVar(FAMGugVectorEntry &ugve) const;
#else
  int GetId() {
    return id;
  }
  void SetId(int newid) {
    id=newid;
  }
#endif

private:
  static const int PROLONGATION_COMP;
  static const int RESTRICTION_COMP;

#ifdef USE_UG_DS
  // if USE_UG_DS, no object of type FAMGTransferEntry may be created; only pointer to
  // FAMGTransferEntry are allowed, since FAMGTransferEntry is only a synonym for a
  // ug transfermatrix entry
  // to check proper use of FAMGTransferEntry you can activate a pure virtual function
  // to ensure that there is no object created, but only pointer.
  // USE ONLY FOR COMPILING CHECKS!!!! DON'T RUN THE PROGRAMM, BECAUSE THE ADDITIONAL
  // POINTER FOR VIRTUAL FUNCTION MECHANISM DESTROYS THE CORRESPONDENCE TO MATRIX-TYPE
  //virtual void check_no_object()=NULL;

  // this entries are only set for debugging purpose
  unsigned int control;         /* object identification, various flags */

  MATRIX *next;                         /* row list                                                     */
  VECTOR *vect;                         /* destination vector					*/

  /* user data */
  DOUBLE value[2];                      /* array of doubles                                     */

#else
  FAMGTransferEntry *next;
  int id;
  double data[2];
#endif
};


#ifdef USE_UG_DS
        #ifdef ONLY_ONE_ALGEBRA_DS
inline FAMGVectorEntry FAMGTransferEntry::GetCol() const {
  return (FAMGVectorEntry)(MDEST((MATRIX*)this));
}
        #else
inline FAMGVectorEntry FAMGTransferEntry::GetCol() const {
  return FAMGVectorEntry(new FAMGugVectorEntryRef(MDEST((MATRIX*)this)));
}
        #endif
inline FAMGTransferEntry *FAMGTransferEntry::GetNext() const {
  return (FAMGTransferEntry*)MNEXT((MATRIX*)this);
}
inline void FAMGTransferEntry::SetNext(FAMGTransferEntry *mat) {
  MNEXT((MATRIX*)this) = (MATRIX*)mat;
}
inline double FAMGTransferEntry::GetProlongation() const {
  return MVALUE((MATRIX*)this,PROLONGATION_COMP);
}
inline double FAMGTransferEntry::GetRestriction() const {
  return MVALUE((MATRIX*)this,RESTRICTION_COMP);
}
inline void FAMGTransferEntry::SetProlongation(double val) {
  MVALUE((MATRIX*)this,PROLONGATION_COMP) = val;
}
inline void FAMGTransferEntry::SetRestriction(double val) {
  MVALUE((MATRIX*)this,RESTRICTION_COMP) = val;
}
// weg inline void FAMGTransferEntry::AddProlongation(double val) {MVALUE((MATRIX*)this,PROLONGATION_COMP) += val;}
// weg inline void FAMGTransferEntry::AddRestriction(double val) {MVALUE((MATRIX*)this,RESTRICTION_COMP) += val;}

// for specialized functions
inline void FAMGTransferEntry::GetColInVar(FAMGugVectorEntry &ugve) const {
  ugve = MDEST((MATRIX*)this);
}
#else
inline FAMGVectorEntry FAMGTransferEntry::GetCol() const {
  return FAMGVectorEntry(new FAMGarrayVectorEntryRef(GetId(id)));
}
inline FAMGTransferEntry *FAMGTransferEntry::GetNext() const {
  return next;
}
inline void FAMGTransferEntry::SetNext(FAMGTransferEntry *mat) {
  next = mat;
}
inline double FAMGTransferEntry::GetProlongation() const {
  return data[PROLONGATION_COMP];
}
inline double FAMGTransferEntry::GetRestriction() const {
  return data[RESTRICTION_COMP];
}
inline void FAMGTransferEntry::SetProlongation(double val) {
  data[PROLONGATION_COMP] = val;
}
inline void FAMGTransferEntry::SetRestriction(double val) {
  data[RESTRICTION_COMP] = val;
}
inline void FAMGTransferEntry::AddProlongation(double val) {
  data[PROLONGATION_COMP] += val;
}
inline void FAMGTransferEntry::AddRestriction(double val) {
  data[RESTRICTION_COMP] += val;
}
#endif


class FAMGTransfer
{
public:
  int Init(class FAMGGrid*);
  FAMGTransferEntry* GetFirstEntry(const FAMGVectorEntry& fg_vec) const;
  FAMGTransferEntry* NewEntry(const FAMGVectorEntry& fg_vec, const FAMGVectorEntry& cg_vec);
  int SetEntries(const FAMGVectorEntry& fg_vec, const FAMGVectorEntry& cg_vec,double prolongation_val, double restriction_val);
  int SetDestinationToCoarse( const FAMGGrid &fg, const FAMGGrid &cg );
#if defined  USE_UG_DS && !defined ONLY_ONE_ALGEBRA_DS
  FAMGTransferEntry* GetFirstEntry(const FAMGugVectorEntry& fg_vec) const;
#else
#endif

private:

#ifdef USE_UG_DS
  GRID *mygrid;
#else
  int n;
  FAMGTransferEntry **row_array;                // array of diagonal elements of the interpolation matrix
#endif
};

#ifdef USE_UG_DS
inline FAMGTransferEntry *FAMGTransfer::GetFirstEntry(const FAMGVectorEntry& fg_vec) const {
  return (FAMGTransferEntry*)VISTART(((FAMGugVectorEntryRef*)(fg_vec.GetPointer()))->myvector());
}
#ifndef ONLY_ONE_ALGEBRA_DS
inline FAMGTransferEntry *FAMGTransfer::GetFirstEntry(const FAMGugVectorEntry& fg_vec) const {
  return (FAMGTransferEntry*)VISTART(fg_vec.myvector());
}
#endif
#else
inline FAMGTransferEntry *FAMGTransfer::GetFirstEntry(const FAMGVectorEntry& fg_vec) const {
  return row_array[((FAMGarrayVectorEntry*)(fg_vec->GetPointer()))->myid()];
}
#endif

#endif
