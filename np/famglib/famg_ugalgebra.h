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

//
// vector stuff
//

class FAMGugVectorEntryRef : public FAMGVectorEntryRef
{
  friend class FAMGugVector;
  friend class FAMGugMatrix;

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
  FAMGugVectorEntryRef( VECTOR *vec ) : vp(vec) {};
  VECTOR*& myvector() {
    return vp;
  };
  VECTOR *vp;
};

class FAMGugVector : public FAMGVector
{
public:
  FAMGugVector( GRID* grid, int vec_comp ) : mygrid(grid), comp(vec_comp) {};

  virtual double& operator[] ( FAMGVectorEntry& ve ) {return VVALUE(((FAMGugVectorEntryRef*)ve.getpointer())->myvector(),comp);};
  virtual int is_valid( FAMGVectorEntry& ve ) {
    return ((FAMGugVectorEntryRef*)ve.getpointer())->myvector()!=NULL;
  };
  virtual int is_end( FAMGVectorEntry& ve ) {
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

//
// matrix stuff
//

class FAMGugMatrixEntryRef : public FAMGMatrixEntryRef
{
  friend class FAMGugMatrix;

public:
  FAMGugMatrixEntryRef() : matp(NULL) {};

  void copy( FAMGugMatrixEntryRef* mep ) {
    *this=*mep;
  };
  virtual FAMGMatrixEntryRef* clone() {
    FAMGugMatrixEntryRef* res = new FAMGugMatrixEntryRef(); res->copy(this); return res;
  };

  virtual FAMGMatrixEntryRef& operator++ () {matp = MNEXT(matp); return *this;};

private:
  FAMGugMatrixEntryRef( MATRIX *mat) : matp(mat) {};
  MATRIX*& myMatrix() {
    return matp;
  };
  MATRIX *matp;
};

class FAMGugMatrix : public FAMGMatrixAlg
{
public:
  FAMGugMatrix( int mat_comp ) : comp(mat_comp) {};

  virtual double& operator[] ( FAMGMatrixEntry& me ) {return MVALUE(((FAMGugMatrixEntryRef*)me.getpointer())->myMatrix(),comp);};
  virtual int is_valid( FAMGVectorEntry& row_ve, FAMGMatrixEntry& me ) {
    return ((FAMGugMatrixEntryRef*)me.getpointer())->myMatrix()!=NULL;
  };
  virtual int is_end( FAMGVectorEntry& row_ve, FAMGMatrixEntry& me ) {
    return ((FAMGugMatrixEntryRef*)me.getpointer())->myMatrix()==NULL;
  };
  virtual FAMGMatrixEntry firstEntry( FAMGVectorEntry& row_ve ) {
    return FAMGMatrixEntry( new FAMGugMatrixEntryRef(VSTART(((FAMGugVectorEntryRef*)row_ve.getpointer())->myvector())));
  };
  virtual const FAMGMatrixEntry endEntry( FAMGVectorEntry& row_ve ) {
    return FAMGMatrixEntry( new FAMGugMatrixEntryRef(NULL));
  };

private:
  int comp;
};

#endif
