// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      famg_onlyugalgebra.h											*/
/*																			*/
/* Purpose:   famg algera classes fixed for ug data structure				*/
/*																			*/
/* Author:    Christian Wrobel												*/
/*			  Institut fuer Computeranwendungen  III						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: ug@ica3.uni-stuttgart.de							*/
/*																			*/
/*																			*/
/* History:   February 99 begin (Christian Wrobel)							*/
/*																			*/
/* Remarks:	  don't include this file into your own files!					*/
/*			  only for use within famg_algebra.h							*/
/*																			*/
/****************************************************************************/

#ifndef __FAMG_ONLYUGALGEBRA__
#define __FAMG_ONLYUGALGEBRA__

/* RCS_ID
   $Header$
 */

extern "C"
{
#include "gm.h"
#include "udm.h"

#ifdef ModelP
#include "np.h"
#endif
}


//
// vector stuff
//
class FAMGVectorEntry; // forward declaration
typedef FAMGVectorEntry FAMGVectorEntryRef;             // synonym; for compatibility only

class FAMGVectorEntry
{
public:
  FAMGVectorEntry() : vp(NULL) {}
  FAMGVectorEntry( const FAMGVectorEntry & ve ) : vp(ve.vp) {}
  FAMGVectorEntry( VECTOR *vptr ) : vp(vptr) {}
  // weg FAMGVectorEntry( FAMGVectorEntryRef* vep ) : vp(vep) {}
  ~FAMGVectorEntry() {}

  FAMGVectorEntry& operator=( const FAMGVectorEntry & ve ) { vp=ve.vp; return *this;}
  FAMGVectorEntry& operator++ () {vp=SUCCVC(vp); return *this;}                 // prefix
  void operator++ (int) {vp=SUCCVC(vp);}                                                                // postfix
  FAMGVectorEntry& operator-- () {vp=PREDVC(vp); return *this;}                 // prefix
  void operator-- (int) {vp=PREDVC(vp);}                                                                // postfix
  int operator==( const FAMGVectorEntry & ve ) const {return (vp==ve.vp);}

  VECTOR *myvector() const {
    return vp;
  }
  FAMGVectorEntryRef* GetPointer() const {
    return (FAMGVectorEntryRef*)this;
  }
  int GetIndex() const {
    return VINDEX(vp);
  }

private:
  VECTOR *vp;
};

class FAMGGridVector
{
public:
  FAMGGridVector( GRID* grid) : mygrid(grid) {}
  ~FAMGGridVector() {};                 // nothing to do

  int is_valid( const FAMGVectorEntry& ve ) const {
    return ve.myvector()!=NULL;
  }
  int is_end( const FAMGVectorEntry& ve ) const {
    return ve.myvector()==NULL;
  }
  int is_beforefirst( const FAMGVectorEntry& ve ) const {
    return ve.myvector()==NULL;
  }
  FAMGVectorEntry firstEntry() const {
    return (FAMGVectorEntry)PFIRSTVECTOR(mygrid);
  }
  FAMGVectorEntry lastEntry() const {
    return (FAMGVectorEntry)LASTVECTOR(mygrid);
  }
  FAMGVectorEntry endEntry() const {
    return (FAMGVectorEntry)NULL;
  }

  int IsCG( const FAMGVectorEntry& ve ) const {
    return VCCOARSE(ve.myvector());
  }
  int IsFG( const FAMGVectorEntry& ve ) const {
    return !IsCG(ve);
  }
  void SetCG( const FAMGVectorEntry& ve ) const {
    SETVCCOARSE(ve.myvector(),1);
  }
  void SetFG( const FAMGVectorEntry& ve ) const {
    SETVCCOARSE(ve.myvector(),0);
  }

  GRID *GetGrid() const {
    return mygrid;
  }

private:
  GRID *mygrid;
};

class FAMGVector
{
  friend class FAMGVectorIter;
  friend class FAMGVectorRevIter;

public:
  typedef class FAMGVectorEntry VectorEntry;
  typedef class FAMGVectorIter Iterator;
  typedef class FAMGVectorRevIter RevIterator;

  FAMGVector(const FAMGGridVector & gridvec) : mygridvector(gridvec) {}
  ~FAMGVector();
  FAMGVector( const FAMGGridVector & gridvec, VECDATA_DESC *vec_desc ) : mygridvector(gridvec), mydesc(vec_desc) {
    allocatedVD=0; assert(VD_IS_SCALAR(vec_desc)); comp = VD_SCALCMP(vec_desc);
  }
  FAMGVector( const FAMGGridVector & gridvec, const FAMGVector &pattern_vec) : mygridvector(gridvec), comp(pattern_vec.comp), mydesc(pattern_vec.mydesc) {}
  FAMGVector* create_new() const;                       // create copy of my; incl. memory for data but without copying the data
  const FAMGGridVector& GetGridVector() const {
    return mygridvector;
  }

  double& operator[] ( const FAMGVectorEntry & ve ) {return VVALUE(ve.myvector(),GetComp());}
  double operator[] ( const FAMGVectorEntry & ve ) const {return VVALUE(ve.myvector(),GetComp());}
  FAMGVector& operator=( const FAMGVector &v ) {CopyValue(*this,v); return *this;}
  FAMGVector& operator+=( const FAMGVector &v ) {AddValue(*this,v); return *this;}
  FAMGVector& operator-=( const FAMGVector &v ) {SubtractValue(*this,v); return *this;}
  double operator=(double val) {SetValue(*this,val); return val;}
  double operator*(const FAMGVector &v) {return ScalProd(*this,v);}                     // scalar product
  FAMGVector& operator*=(double scale) {Scale(*this,scale); return *this;}

  int is_valid( const FAMGVectorEntry& ve ) const {
    return ve.myvector()!=NULL;
  }
  int is_end( const FAMGVectorEntry& ve ) const {
    return ve.myvector()==NULL;
  }
  int is_beforefirst( const FAMGVectorEntry& ve ) {
    return ve.myvector()==NULL;
  }
  FAMGVectorEntry firstEntry(const FAMGVectorEntry& ve) const {
    return GetGridVector().firstEntry();
  }
  FAMGVectorEntry lastEntry() const {
    return GetGridVector().lastEntry();
  }
  FAMGVectorEntry endEntry() const {
    return (FAMGVectorEntry)NULL;
  }

  int IsCG( const FAMGVectorEntry& ve ) const {
    return VCCOARSE(ve.myvector());
  }
  int IsFG( const FAMGVectorEntry& ve ) const {
    return !IsCG(ve);
  }
  void SetCG( const FAMGVectorEntry& ve ) const {
    SETVCCOARSE(ve.myvector(),1);
  }
  void SetFG( const FAMGVectorEntry& ve ) const {
    SETVCCOARSE(ve.myvector(),0);
  }

  double norm() const {
    return ::norm(*this);
  }
  double sum() const {
    return ::sum(*this);
  }
  void AddScaledVec( double scale, const FAMGVector &source ) {
    ::AddScaledValue( *this, scale, source);
  }
  void VecMinusMatVec( const FAMGVector &rhs, const FAMGMatrixAlg &mat, const FAMGVector &sol ) {
    ::VecMinusMatVec( *this, rhs, mat, sol);
  }
  void MatVec( const FAMGMatrixAlg &mat, const FAMGVector &source ) {
    ::MatVec( *this, mat, source);
  }

  void JacobiSmoother( const FAMGMatrixAlg &mat, const FAMGVector &def ) {
    ::JacobiSmoother( *this, mat, def );
  }
  void dampedJacobiSmoother( const FAMGMatrixAlg &mat, const FAMGVector &def ) {
    ::dampedJacobiSmoother( *this, mat, def );
  }
  void FGSSmoother( const FAMGMatrixAlg &mat, FAMGVector &def ) {
    ::FGSSmoother( *this, mat, def );
  }
  void BGSSmoother( const FAMGMatrixAlg &mat, FAMGVector &def ) {
    ::BGSSmoother( *this, mat, def );
  }
  void SGSSmoother( const FAMGMatrixAlg &mat, FAMGVector &def ) {
    ::SGSSmoother( *this, mat, def );
  }
  void JacobiSmoothFG( const FAMGMatrixAlg &mat, const FAMGVector &def ) {
    ::JacobiSmoothFG( *this, mat, def );
  }

  VECDATA_DESC *GetUgVecDesc () const {
    return mydesc;
  }
protected:
  const FAMGGridVector &mygridvector;

private:
  int GetComp() const {
    return comp;
  }
  VECDATA_DESC *mydesc;
  int comp;
  int allocatedVD;                      // 1 if a new VECDATA_DESC was allocated for this vector
};

class FAMGVectorIter
{
public:
  FAMGVectorIter( const FAMGGridVector & gv ) : mygridvector(gv), current_ve(gv.firstEntry()) {}
  FAMGVectorIter( const FAMGVector & v ) : mygridvector(v.GetGridVector()), current_ve(v.GetGridVector().firstEntry()) {}

  int operator() ( FAMGVectorEntry& ve ) {
    int res = !mygridvector.is_end(ve=current_ve); if(res) ++current_ve;return res;
  }
  void reset() {
    current_ve = mygridvector.firstEntry();
  }
private:
  const FAMGGridVector& mygridvector;
  FAMGVectorEntry current_ve;
};

class FAMGVectorRevIter
{
public:
  FAMGVectorRevIter( const FAMGGridVector & gv ) : mygridvector(gv), current_ve(gv.lastEntry()) {}
  FAMGVectorRevIter( const FAMGVector & v ) : mygridvector(v.mygridvector), current_ve(v.GetGridVector().lastEntry()) {}

  int operator() ( FAMGVectorEntry& ve ) {
    int res = !mygridvector.is_beforefirst(ve=current_ve); if(res) --current_ve;return res;
  }
  void reset() {
    current_ve = mygridvector.lastEntry();
  }
private:
  const FAMGGridVector &mygridvector;
  FAMGVectorEntry current_ve;
};

//
// matrix stuff
//

class FAMGMatrixEntry
{
  friend class FAMGMatrixAlg;
  friend class FAMGMatrixIter;

public:
  FAMGMatrixEntry() : matp(NULL) {}
  FAMGMatrixEntry( const FAMGMatrixEntry & me ) : matp(me.GetMyMatrix()) {}
  ~FAMGMatrixEntry() {}                 // nothing to do

  FAMGMatrixEntry& operator=( const FAMGMatrixEntry & me ) { matp=me.GetMyMatrix(); return *this;}
  FAMGMatrixEntry& operator++ () {matp = MNEXT(matp); return *this;}                                            // prefix
  FAMGMatrixEntry& operator++ (int) {FAMGMatrixEntry& tmp = *this; ++*this; return tmp;}                // postfix

  FAMGVectorEntry dest() const {
    return (FAMGVectorEntry)MDEST(GetMyMatrix());
  }

private:
  FAMGMatrixEntry( MATRIX *mat) : matp(mat) {}
  MATRIX* GetMyMatrix() const {
    return matp;
  }
  MATRIX*& myMatrix() {
    return matp;
  }
  MATRIX *matp;
};

class FAMGMatrixAlg
{
public:
  typedef class FAMGVector Vector;
  typedef class FAMGGridVector GridVector;
  typedef class FAMGMatrixEntry MatrixEntry;
  typedef class FAMGMatrixIter Iterator;

  FAMGMatrixAlg( int nr_vecs, int nr_links ) : n(nr_vecs), nlinks(nr_links) {}
  FAMGMatrixAlg( GRID *grid, MATDATA_DESC *md, int nrVec, int nrLink ) : n(nrVec), nlinks(nrLink), mygrid(grid), matdesc(md), comp(MD_SCALCMP(md)) {
    assert(MD_IS_SCALAR(md));
  }
  FAMGMatrixAlg( GRID *grid, MATDATA_DESC *md ) : n(0), nlinks(0), mygrid(grid), matdesc(md), comp(MD_SCALCMP(md)) {
    assert(MD_IS_SCALAR(md));
  }                                                                                                                                                          // CAUTION: set N and NLinks explicitly
  FAMGMatrixAlg( GRID *grid, const FAMGMatrixAlg &pattern_mat) :  n(0), nlinks(0), mygrid(grid) {
    matdesc = NULL; if(pattern_mat.GetMatDesc()->locked) matdesc = pattern_mat.GetMatDesc();else AllocMDFromMD(MYMG(grid), GLEVEL(grid), GLEVEL(grid), pattern_mat.GetMatDesc(), &matdesc);assert(matdesc!=NULL); assert(MD_IS_SCALAR(matdesc)); comp = MD_SCALCMP(matdesc);
  }
  ~FAMGMatrixAlg() {
    if (FreeMD(MYMG(GetMyGrid()),GLEVEL(GetMyGrid()),GLEVEL(GetMyGrid()),GetMatDesc())) assert(0);
  }

  double operator[] ( const FAMGMatrixEntry & me ) const {return MVALUE(me.GetMyMatrix(),GetComp());}
  double& operator[] ( const FAMGMatrixEntry & me ) {return MVALUE(me.GetMyMatrix(),GetComp());}
  int is_valid( const FAMGVectorEntry& row_ve, const FAMGMatrixEntry& me ) const {
    return me.GetMyMatrix()!=NULL;
  }
  int is_end( const FAMGVectorEntry& row_ve, const FAMGMatrixEntry& me ) const {
    return me.GetMyMatrix()==NULL;
  }
  FAMGMatrixEntry firstEntry( const FAMGVectorEntry& row_ve ) const {
    return (FAMGMatrixEntry)VSTART(row_ve.myvector());
  }
  FAMGMatrixEntry endEntry( const FAMGVectorEntry& row_ve ) const {
    return (FAMGMatrixEntry)NULL;
  }
  double DiagValue( const FAMGVectorEntry& row_ve ) const {
    return MVALUE(VSTART(row_ve.myvector()),GetComp());
  }
  double GetAdjData( const FAMGMatrixEntry& me ) const {
    return MVALUE(MADJ(me.GetMyMatrix()),GetComp());
  }
  void AddEntry(double mval, const FAMGVectorEntry &row, const FAMGVectorEntry &col);
  int ConstructGalerkinMatrix( const FAMGGrid &fg ) {
    return ::ConstructGalerkinMatrix(*this, fg);
  }

  int &GetN() {
    return n;
  }
  int &GetNLinks() {
    return nlinks;
  }
  MATDATA_DESC *GetMatDesc() const {
    return matdesc;
  }

private:
  int n;                        // number of vectors
  int nlinks;                   // number of links (==entries) in the matrix
  GRID *GetMyGrid() {
    return mygrid;
  }
  int GetComp() const {
    return comp;
  }
  int comp;
  GRID *mygrid;
  MATDATA_DESC *matdesc;
};

class FAMGMatrixIter
{
public:
  FAMGMatrixIter( const FAMGMatrixAlg & m, const FAMGVectorEntry & row_ve) : myMatrix(m), myrow_ve(row_ve), current_me(m.firstEntry(row_ve)) {}

  int operator() ( FAMGMatrixEntry& me ) {
    me=current_me; int res = (current_me.GetMyMatrix()!=NULL); if(res) ++current_me;return res;
  }
  void reset() {
    current_me = myMatrix.firstEntry(myrow_ve);
  }

private:
  const FAMGMatrixAlg& myMatrix;
  const FAMGVectorEntry& myrow_ve;
  FAMGMatrixEntry current_me;
};

#endif
