// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      famg_ugalgebra.h												*/
/*																			*/
/* Purpose:   famg matrix classes											*/
/*																			*/
/* Author:    Christian Wrobel												*/
/*			  Institut fuer Computeranwendungen  III						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: ug@ica3.uni-stuttgart.de							*/
/*																			*/
/*																			*/
/* History:   August 98 begin (Christian Wrobel)							*/
/*																			*/
/* Remarks:																	*/
/*																			*/
/****************************************************************************/

#ifndef __FAMG_UGALGEBRA__
#define __FAMG_UGALGEBRA__

#include "famg_algebra.h"

/* RCS_ID
   $Header$
 */

class FAMGugVectorEntryRef : public FAMGVectorEntryRef
{
  friend class FAMGugVector;

public:
  FAMGugVectorEntryRef() : vp(NULL) {};

  void copy( FAMGugVectorEntryRef* vep ) {
    *this=*vep;
  };
  virtual FAMGVectorEntryRef* clone() {
    FAMGugVectorEntryRef* res = new FAMGugVectorEntryRef(); res->copy(this); return res;
  };

  virtual FAMGVectorEntryRef& operator++ () {vp = SUCCVC(vp); return *this;};
  virtual FAMGVectorEntryRef& operator-- () {vp = PREDVC(vp); return *this;};

private:
  FAMGugVectorEntryRef( VECTOR *vec) : vp(vec) {};
  VECTOR* & myvector() {
    return vp;
  };
  VECTOR *vp;
};

class FAMGugVector : public FAMGVector
{
public:
  FAMGugVector( GRID* grid, int vec_comp ) : mygrid(grid), comp(vec_comp) {};

  virtual double & operator[] (FAMGVectorEntry& ve) {return VVALUE(((FAMGugVectorEntryRef*)ve.getpointer())->myvector(),comp);};
  virtual int is_valid(FAMGVectorEntry& ve) {
    return ((FAMGugVectorEntryRef*)ve.getpointer())->myvector()!=NULL;
  };
  virtual int is_end(FAMGVectorEntry& ve) {
    return ((FAMGugVectorEntryRef*)ve.getpointer())->myvector()==NULL;
  };
  virtual FAMGVectorEntry firstEntry() {
    return FAMGVectorEntry( new FAMGugVectorEntryRef(FIRSTVECTOR(mygrid)));
  };
  virtual const FAMGVectorEntry endEntry() {
    return FAMGVectorEntry( new FAMGugVectorEntryRef(NULL));
  };

private:
  GRID *mygrid;
  int comp;
};

#endif
