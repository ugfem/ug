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

extern "C"
{
#include "gm.h"
#include "udm.h"

#ifdef ModelP
#include "np.h"
#endif
}

#include "famg_algebra.h"

/* RCS_ID
   $Header$
 */

//
// vector stuff
//

#ifdef ModelP
// auxiliaries for ug VECTOR
#define IS_FAMG_MASTER(vec) (PRIO(vec)==PrioMaster)
#define IS_FAMG_GHOST(vec)  (PRIO(vec)<=PrioBorder)
#endif

class FAMGugVectorEntryRef : public FAMGVectorEntryRef
{
public:
  friend class FAMGugVector;
  friend class FAMGugMatrix;
  friend class FAMGugVectorIter;

  FAMGugVectorEntryRef() : vp(NULL) {}
  FAMGugVectorEntryRef( VECTOR *vec ) : vp(vec) {}

  void copy( const FAMGugVectorEntryRef* vep ) {
    *this=*vep;
  }
  virtual FAMGVectorEntryRef* clone() {
    FAMGugVectorEntryRef* res = new FAMGugVectorEntryRef(); res->copy(this); return res;
  }

  virtual FAMGVectorEntryRef& operator++ () {vp = SUCCVC(vp); return *this;}
  virtual FAMGVectorEntryRef& operator-- () {vp = PREDVC(vp); return *this;}

  VECTOR *myvector() const {
    return vp;
  }
  virtual int GetIndex() const {
    return VINDEX(myvector());
  }

protected:
#ifdef __SGI10__
  virtual long comparable_value() const {
    return (long)myvector();
  }
#else
  virtual size_t comparable_value() const {
    return (size_t)myvector();
  }
#endif
private:
  VECTOR *vp;
};

class FAMGugGridVector : public FAMGGridVector
{
public:
  typedef class FAMGugVectorEntry VectorEntry;
  typedef class FAMGugVectorIter Iterator;

  friend class FAMGugVectorIter;
  friend class FAMGugVectorRevIter;

  FAMGugGridVector( GRID* grid) : mygrid(grid) {}
  //virtual ~FAMGugGridVector()	{};	// nothing to do

  virtual int is_valid( const FAMGVectorEntry& ve ) const {
    return ((FAMGugVectorEntryRef*)ve.GetPointer())->myvector()!=NULL;
  }
  virtual int is_end( const FAMGVectorEntry& ve ) const {
    return ((FAMGugVectorEntryRef*)ve.GetPointer())->myvector()==NULL;
  }
  virtual int is_beforefirst( const FAMGVectorEntry& ve ) const {
    return ((FAMGugVectorEntryRef*)ve.GetPointer())->myvector()==NULL;
  }
  virtual FAMGVectorEntry firstEntry() const {
    return FAMGVectorEntry( new FAMGugVectorEntryRef(PFIRSTVECTOR(mygrid)));
  }
  virtual FAMGVectorEntry lastEntry() const {
    return FAMGVectorEntry( new FAMGugVectorEntryRef(LASTVECTOR(mygrid)));
  }
  virtual FAMGVectorEntry endEntry() const {
    return FAMGVectorEntry( new FAMGugVectorEntryRef(NULL));
  }

  virtual int IsCG( const FAMGVectorEntry& ve ) const {
    return VCCOARSE(((FAMGugVectorEntryRef*)ve.GetPointer())->myvector());
  }
  virtual int IsFG( const FAMGVectorEntry& ve ) const {
    return !IsCG(ve);
  }
  virtual void SetCG( const FAMGVectorEntry& ve ) const {
    SETVCCOARSE(((FAMGugVectorEntryRef*)ve.GetPointer())->myvector(),1);
  }
  virtual void SetFG( const FAMGVectorEntry& ve ) const {
    SETVCCOARSE(((FAMGugVectorEntryRef*)ve.GetPointer())->myvector(),0);
  }

  GRID *GetGrid() const {
    return mygrid;
  }

  // only for specialized functions
  int is_valid( const FAMGugVectorEntry& ve ) const;
  int is_end( const FAMGugVectorEntry& ve ) const;
  int is_beforefirst( const FAMGugVectorEntry& ve ) const;

  int IsCG( const FAMGugVectorEntry& ve ) const;
  int IsFG( const FAMGugVectorEntry& ve ) const;
  void SetCG( const FAMGugVectorEntry& ve ) const;
  void SetFG( const FAMGugVectorEntry& ve ) const;
private:
  GRID *mygrid;
};

class FAMGugVector : public FAMGVector
{
public:
  typedef class FAMGugVectorEntry VectorEntry;
  typedef class FAMGugVectorIter Iterator;
  typedef class FAMGugVectorRevIter RevIterator;

  friend class FAMGugVectorIter;
  friend class FAMGugVectorRevIter;

  FAMGugVector( const FAMGugGridVector & gridvec, VECDATA_DESC *vec_desc ) : FAMGVector(gridvec), mydesc(vec_desc) {
    allocatedVD=0; assert(VD_IS_SCALAR(vec_desc)); comp = VD_SCALCMP(vec_desc);
  }
  FAMGugVector( const FAMGGridVector & gridvec, const FAMGugVector &pattern_vec) : FAMGVector(gridvec), comp(pattern_vec.comp), mydesc(pattern_vec.mydesc) {}
  virtual ~FAMGugVector();

  virtual double& operator[] ( const FAMGVectorEntry & ve ) {return VVALUE(((FAMGugVectorEntryRef*)ve.GetPointer())->myvector(),GetComp());}
  virtual double operator[] ( const FAMGVectorEntry & ve ) const {return VVALUE(((FAMGugVectorEntryRef*)ve.GetPointer())->myvector(),GetComp());}
  virtual FAMGVector& operator=( const FAMGVector &v );
  virtual FAMGVector& operator+=( const FAMGVector &v );
  virtual FAMGVector& operator-=( const FAMGVector &v );
  virtual double operator=(double);
  virtual double operator*( const FAMGVector &v );                      // scalar product
  virtual FAMGVector& operator*=( double scale );

  virtual FAMGVector* create_new() const;

  virtual double norm() const;
  virtual double sum() const;
  virtual void AddScaledVec( double scale, const FAMGVector &source );
  virtual void VecMinusMatVec( const FAMGVector &rhs, const FAMGMatrixAlg &mat, const FAMGVector &sol );
  virtual void MatVec( const FAMGMatrixAlg &mat, const FAMGVector &source );

  virtual void JacobiSmoother( const FAMGMatrixAlg &mat, const FAMGVector &def );
  virtual void dampedJacobiSmoother( const FAMGMatrixAlg &mat, const FAMGVector &def );
  virtual void FGSSmoother( const FAMGMatrixAlg &mat, FAMGVector &def );
  virtual void BGSSmoother( const FAMGMatrixAlg &mat, FAMGVector &def );
  virtual void SGSSmoother( const FAMGMatrixAlg &mat, FAMGVector &def );
  virtual void JacobiSmoothFG( const FAMGMatrixAlg &mat, const FAMGVector &def );

  VECDATA_DESC *GetUgVecDesc () const;
  // only for specialized functions
  double& operator[] ( const FAMGugVectorEntry & ve );
  double operator[] ( const FAMGugVectorEntry & ve ) const;
  int is_valid( const FAMGugVectorEntry& ve ) const;
  int is_end( const FAMGugVectorEntry& ve ) const;
  int is_beforefirst( const FAMGugVectorEntry& ve ) const;
  int IsCG( const FAMGugVectorEntry& ve ) const;
  int IsFG( const FAMGugVectorEntry& ve ) const;
  void SetCG( const FAMGugVectorEntry& ve );
  void SetFG( const FAMGugVectorEntry& ve );

private:
  FAMGugVector(const FAMGugGridVector & gridvec) : FAMGVector(gridvec), mydesc(NULL) {
    allocatedVD=0;
  };                                                                                                                    // only for create_new
  int GetComp() const {
    return comp;
  }
  VECDATA_DESC *mydesc;
  int comp;
  int allocatedVD;                      // 1 if a new VECDATA_DESC was allocated for this vector
};

// specialized classes to increase performance

class FAMGugVectorEntry
{
public:
  FAMGugVectorEntry() : vp(NULL) {}
  FAMGugVectorEntry( VECTOR* vecp ) : vp(vecp) {}

  FAMGugVectorEntry& operator=( const FAMGugVectorEntry & ve ) {vp=ve.myvector(); return *this;}
  FAMGugVectorEntry& operator=( VECTOR *vep ) {vp=vep; return *this;}
  FAMGugVectorEntry& operator++ () {vp=SUCCVC(vp); return *this;}               // prefix
  void operator++ (int) {vp=SUCCVC(vp);}                                                                // postfix
  FAMGugVectorEntry& operator-- () {vp=PREDVC(vp); return *this;}               // prefix
  void operator-- (int) {vp=PREDVC(vp); }                                                               // postfix

  VECTOR* myvector() const {
    return vp;
  }
  int GetIndex() const {
    return VINDEX(myvector());
  }

private:
  VECTOR *vp;
};

class FAMGugVectorIter
{
public:
  FAMGugVectorIter( const FAMGugGridVector & gv ) : first_vp(PFIRSTVECTOR(gv.mygrid)) {
    current_vp=first_vp;
  }
  FAMGugVectorIter( const FAMGugVector & v ) : first_vp(PFIRSTVECTOR(((FAMGugGridVector&)(v.GetGridVector())).mygrid)) {
    current_vp=first_vp;
  }

  int operator() ( FAMGugVectorEntry& ve ) {
    ve=current_vp; int res = (current_vp!=NULL); if(res) current_vp=SUCCVC(current_vp);return res;
  }
  void reset() {
    current_vp = first_vp;
  }
private:
  VECTOR *current_vp, *first_vp;
};

class FAMGugVectorRevIter
{
public:
  FAMGugVectorRevIter( const FAMGugGridVector & gv ) : last_vp(PFIRSTVECTOR(gv.mygrid)) {
    current_vp=last_vp;
  }
  FAMGugVectorRevIter( const FAMGugVector & v ) : last_vp(LASTVECTOR(((FAMGugGridVector&)(v.GetGridVector())).mygrid)) {
    current_vp=last_vp;
  }

  int operator() ( FAMGugVectorEntry& ve ) {
    ve=current_vp; int res = (current_vp!=NULL); if(res) current_vp=PREDVC(current_vp);return res;
  }
  void reset() {
    current_vp = last_vp;
  }
private:
  VECTOR *current_vp, *last_vp;
};



//
// matrix stuff
//

class FAMGugMatrixEntryRef : public FAMGMatrixEntryRef
{
  friend class FAMGugMatrix;

public:
  FAMGugMatrixEntryRef() : matp(NULL) {}

  void copy( FAMGugMatrixEntryRef* mep ) {
    *this=*mep;
  }
  virtual FAMGMatrixEntryRef* clone() {
    FAMGugMatrixEntryRef* res = new FAMGugMatrixEntryRef(); res->copy(this); return res;
  }

  virtual FAMGMatrixEntryRef& operator++ () {matp = MNEXT(matp); return *this;}

  virtual FAMGVectorEntry dest() const {
    return FAMGVectorEntry( new FAMGugVectorEntryRef(MDEST(GetMyMatrix())));
  }

private:
  FAMGugMatrixEntryRef( MATRIX *mat) : matp(mat) {}
  MATRIX* GetMyMatrix() const {
    return matp;
  }
  MATRIX*& myMatrix() {
    return matp;
  }
  MATRIX *matp;
};

class FAMGugMatrix : public FAMGMatrixAlg
{
public:
  typedef class FAMGugVector Vector;
  typedef class FAMGugGridVector GridVector;
  typedef class FAMGugMatrixEntry MatrixEntry;
  typedef class FAMGugMatrixIter Iterator;

  FAMGugMatrix( GRID *grid, MATDATA_DESC *md, int nrVec, int nrLink ) : FAMGMatrixAlg(nrVec,nrLink), mygrid(grid), matdesc(md), comp(MD_SCALCMP(md)) {
    assert(MD_IS_SCALAR(md));
  }
  FAMGugMatrix( GRID *grid, MATDATA_DESC *md ) : FAMGMatrixAlg(0,0), mygrid(grid), matdesc(md), comp(MD_SCALCMP(md)) {
    assert(MD_IS_SCALAR(md));
  }                                                                                                                                                            // CAUTION: set N and NLinks explicitly
  FAMGugMatrix( GRID *grid, const FAMGugMatrix &pattern_mat) : FAMGMatrixAlg(0,0), mygrid(grid) {
    matdesc = NULL; if(pattern_mat.GetMatDesc()->locked) matdesc = pattern_mat.GetMatDesc();else AllocMDFromMD(MYMG(grid), GLEVEL(grid), GLEVEL(grid), pattern_mat.GetMatDesc(), &matdesc);assert(matdesc!=NULL); assert(MD_IS_SCALAR(matdesc)); comp = MD_SCALCMP(matdesc);
  }
  ~FAMGugMatrix() {
    if (FreeMD(MYMG(GetMyGrid()),GLEVEL(GetMyGrid()),GLEVEL(GetMyGrid()),GetMatDesc())) assert(0);
  }

  //void SetNumbers();

  virtual double operator[] ( const FAMGMatrixEntry & me ) const {return MVALUE(((FAMGugMatrixEntryRef*)me.GetPointer())->myMatrix(),GetComp());}
  virtual double& operator[] ( const FAMGMatrixEntry & me ) {return MVALUE(((FAMGugMatrixEntryRef*)me.GetPointer())->myMatrix(),GetComp());}
  virtual int is_valid( const FAMGVectorEntry& row_ve, const FAMGMatrixEntry& me ) const {
    return ((FAMGugMatrixEntryRef*)me.GetPointer())->GetMyMatrix()!=NULL;
  }
  virtual int is_end( const FAMGVectorEntry& row_ve, const FAMGMatrixEntry& me ) const {
    return ((FAMGugMatrixEntryRef*)me.GetPointer())->GetMyMatrix()==NULL;
  }
  virtual FAMGMatrixEntry firstEntry( const FAMGVectorEntry& row_ve ) const {
    return FAMGMatrixEntry( new FAMGugMatrixEntryRef(VSTART(((FAMGugVectorEntryRef*)row_ve.GetPointer())->myvector())));
  }
  virtual FAMGMatrixEntry endEntry( const FAMGVectorEntry& row_ve ) const {
    return FAMGMatrixEntry( new FAMGugMatrixEntryRef(NULL));
  }
  virtual double DiagValue( const FAMGVectorEntry& row_ve ) const {
    return MVALUE(VSTART(((FAMGugVectorEntryRef*)row_ve.GetPointer())->myvector()),GetComp());
  }

  virtual double GetAdjData( const FAMGMatrixEntry& me ) const {
    return MVALUE(MADJ(((FAMGugMatrixEntryRef*)me.GetPointer())->myMatrix()),GetComp());
  }

  virtual int ConstructGalerkinMatrix( const FAMGGrid &fg );

  MATDATA_DESC *GetMatDesc() const {
    return matdesc;
  }

  // only for specialized functions
  double& operator[] ( const FAMGugMatrixEntry & me );
  double operator[] ( const FAMGugMatrixEntry & me ) const;
  double DiagValue( const FAMGugVectorEntry& row_ve ) const;
  void AddEntry(double mval, const FAMGugVectorEntry &row, const FAMGugVectorEntry &col);
private:
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

// specialized classes to increase performance

class FAMGugMatrixEntry
{
public:
  FAMGugMatrixEntry() : matp(NULL) {}
  FAMGugMatrixEntry(MATRIX *mat) : matp(mat) {}

  FAMGugMatrixEntry& operator++ () {matp = MNEXT(matp); return *this;}
  FAMGugMatrixEntry& operator=(MATRIX *mp) {matp = mp; return *this;}
  FAMGugVectorEntry dest() const {
    return FAMGugVectorEntry(MDEST(GetMyMatrix()));
  }

  MATRIX* GetMyMatrix() const {
    return matp;
  }
private:
  MATRIX *matp;
};

class FAMGugMatrixIter
{
public:
  FAMGugMatrixIter( const FAMGugMatrix & m, const FAMGugVectorEntry & row_ve) : current_mp(VSTART(row_ve.myvector())) {}

  int operator() ( FAMGugMatrixEntry& me ) {
    me=current_mp; int res = (current_mp!=NULL); if(res) current_mp=MNEXT(current_mp);return res;
  }

private:
  MATRIX *current_mp;
};


//
// inline implementation for GridVector
//
inline int FAMGugGridVector::is_valid( const FAMGugVectorEntry& ve ) const {
  return ve.myvector()!=NULL;
}
inline int FAMGugGridVector::is_end( const FAMGugVectorEntry& ve ) const {
  return ve.myvector()==NULL;
}
inline int FAMGugGridVector::is_beforefirst( const FAMGugVectorEntry& ve ) const {
  return ve.myvector()==NULL;
}

inline int FAMGugGridVector::IsCG( const FAMGugVectorEntry& ve ) const {
  return VCCOARSE(ve.myvector());
}
inline int FAMGugGridVector::IsFG( const FAMGugVectorEntry& ve ) const {
  return !IsCG(ve);
}
inline void FAMGugGridVector::SetCG( const FAMGugVectorEntry& ve ) const {
  SETVCCOARSE(ve.myvector(),1);
}
inline void FAMGugGridVector::SetFG( const FAMGugVectorEntry& ve ) const {
  SETVCCOARSE(ve.myvector(),0);
}

//
// inline implementation for vector
//
inline double& FAMGugVector::operator[] ( const FAMGugVectorEntry & ve ) {return VVALUE(ve.myvector(),GetComp());}
inline double FAMGugVector::operator[] ( const FAMGugVectorEntry & ve ) const {return VVALUE(ve.myvector(),GetComp());}
inline FAMGVector &FAMGugVector::operator=( const FAMGVector &v ) { CopyValue(*this,(FAMGugVector&)v); return *this;}
inline FAMGVector &FAMGugVector::operator+=( const FAMGVector &v ) { AddValue(*this,(FAMGugVector&)v); return *this;}
inline FAMGVector &FAMGugVector::operator-=( const FAMGVector &v ) { SubtractValue(*this,(FAMGugVector&)v); return *this;}
inline double FAMGugVector::operator=( double val ) { SetValue(*this,val); return val;}
inline double FAMGugVector::operator*( const FAMGVector &v ) { return ScalProd(*this,(FAMGugVector&)v);}
inline FAMGVector & FAMGugVector::operator*=( double scale ) {Scale(*this,scale); return *this;}
inline double FAMGugVector::norm() const {
  return ::norm(*this);
}                                                                // cast is very ugly!
inline double FAMGugVector::sum() const {
  return ::sum(*this);
}                                                              // cast is very ugly!

inline int FAMGugVector::is_valid( const FAMGugVectorEntry& ve ) const {
  return ((FAMGugGridVector&)GetGridVector()).is_valid(ve);
}
inline int FAMGugVector::is_end( const FAMGugVectorEntry& ve ) const {
  return ((FAMGugGridVector&)GetGridVector()).is_end(ve);
}
inline int FAMGugVector::is_beforefirst( const FAMGugVectorEntry& ve ) const {
  return ((FAMGugGridVector&)GetGridVector()).is_beforefirst(ve);
}
inline int FAMGugVector::IsCG( const FAMGugVectorEntry& ve ) const {
  return ((FAMGugGridVector&)GetGridVector()).IsCG(ve);
}
inline int FAMGugVector::IsFG( const FAMGugVectorEntry& ve ) const {
  return ((FAMGugGridVector&)GetGridVector()).IsFG(ve);
}
inline void FAMGugVector::SetCG( const FAMGugVectorEntry& ve ) {
  ((FAMGugGridVector&)GetGridVector()).SetCG(ve);
}
inline void FAMGugVector::SetFG( const FAMGugVectorEntry& ve ) {
  ((FAMGugGridVector&)GetGridVector()).SetFG(ve);
}

inline VECDATA_DESC *FAMGugVector::GetUgVecDesc () const {
  return mydesc;
}

inline void FAMGugVector::AddScaledVec( double scale, const FAMGVector &source )
{
  ::AddScaledValue( *this, scale, (FAMGugVector&)source);
}

inline void FAMGugVector::VecMinusMatVec( const FAMGVector &rhs, const FAMGMatrixAlg &mat, const FAMGVector &sol )
{
  ::VecMinusMatVec( *this, (FAMGugVector&)rhs, (FAMGugMatrix&)mat, (FAMGugVector&)sol);
}

inline void FAMGugVector::MatVec( const FAMGMatrixAlg &mat, const FAMGVector &source )
{
  ::MatVec( *this, (FAMGugMatrix&)mat, (FAMGugVector&)source);
}

inline void FAMGugVector::JacobiSmoother( const FAMGMatrixAlg &mat, const FAMGVector &def )
{
  ::JacobiSmoother( *this, (FAMGugMatrix&)mat, (FAMGugVector&)def );
}

inline void FAMGugVector::dampedJacobiSmoother( const FAMGMatrixAlg &mat, const FAMGVector &def )
{
  ::dampedJacobiSmoother( *this, (FAMGugMatrix&)mat, (FAMGugVector&)def );
}

inline void FAMGugVector::FGSSmoother( const FAMGMatrixAlg &mat, FAMGVector &def )
{
  ::FGSSmoother( *this, (FAMGugMatrix&)mat, (FAMGugVector&)def );
}

inline void FAMGugVector::BGSSmoother( const FAMGMatrixAlg &mat, FAMGVector &def )
{
  ::BGSSmoother( *this, (FAMGugMatrix&)mat, (FAMGugVector&)def );
}

inline void FAMGugVector::SGSSmoother( const FAMGMatrixAlg &mat, FAMGVector &def )
{
  ::SGSSmoother( *this, (FAMGugMatrix&)mat, (FAMGugVector&)def );
}

inline void FAMGugVector::JacobiSmoothFG( const FAMGMatrixAlg &mat, const FAMGVector &def )
{
  ::JacobiSmoothFG( *this, (FAMGugMatrix&)mat, (FAMGugVector&)def );
}

//
// inline implementation for matrix
//

inline double& FAMGugMatrix::operator[] ( const FAMGugMatrixEntry & me ) {return MVALUE(me.GetMyMatrix(),GetComp());}
inline double FAMGugMatrix::operator[] ( const FAMGugMatrixEntry & me ) const {return MVALUE(me.GetMyMatrix(),comp);}
inline double FAMGugMatrix::DiagValue( const FAMGugVectorEntry& row_ve ) const {
  return MVALUE(VSTART(row_ve.myvector()),GetComp());
}
inline int FAMGugMatrix::ConstructGalerkinMatrix( const FAMGGrid &fg )
{
  return ::ConstructGalerkinMatrix(*this, fg);
}
#endif
