// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      famg_algebra.h												*/
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

#ifndef __FAMG_ALGEBRA__
#define __FAMG_ALGEBRA__

/* RCS_ID
   $Header$
 */

// if defined, only one type of alegra data structure is active (ug DS if USE_UG_DS is
// defined otherwise array DS); this avoids the use of virtual classes and is thus
// much faster
#define ONLY_ONE_ALGEBRA_DS

#ifdef ModelP

#define FAMG_INNER_FIRST

// FAMG_SINGLESTEP works only if the given grid has full geometric infos.
// This is true for level 0 grid from UG or if FAMG_FULL_OVERLAP is active.
// The geometric info is used to force the same elimination sequence independent
// of the number of processors (incl. 1 PE).
//#define FAMG_SINGLESTEP
//#define FAMG_SINGLESTEP_FULL_OUTPUT
//#define FAMG_FULL_OVERLAP

#endif

// TODO: remove it:
//#define SIMULATE_HALFENING

//#define PROTOCOLNUMERIC

// TODO: remove
extern "C" int PrintLocal(char *);


// Use RETURN(x) it instead of return(x) if you want to have an error check.
// This is recommended for all functions containing global communications to avoid
// deadlocks (an assert is much easyier to debug).
#if !defined USE_UG_DS && !defined RETURN
#ifdef Debug
#define RETURN(rcode)   {int rc = rcode; assert(!rc); return (rc);}
#else
#define RETURN(rcode)   {return (rc);}
#endif
#endif


// forward declarations
class FAMGGraph;
class FAMGGrid;
class FAMGVectorIter;
class FAMGMatrixAlg;

#ifndef ONLY_ONE_ALGEBRA_DS

//
// vector stuff
//

class FAMGVectorEntryRef        // a single VectorEntryRef must be able to denote all corresponding entries in different vectors
{       // abstract base class
public:
  virtual FAMGVectorEntryRef* clone() = 0;

  virtual FAMGVectorEntryRef& operator++ () = 0;                // prefix
  FAMGVectorEntryRef& operator++ (int) {FAMGVectorEntryRef& tmp = *this; ++*this; return tmp;}              // postfix
  virtual FAMGVectorEntryRef& operator-- () = 0;                // prefix
  FAMGVectorEntryRef& operator-- (int) {FAMGVectorEntryRef& tmp = *this; --*this; return tmp;}              // postfix
  int operator==( const FAMGVectorEntryRef & ve ) const {return (comparable_value()==ve.comparable_value());}

  virtual int GetIndex() const = 0;

protected:
#ifdef __SGI10__
  virtual long comparable_value() const = 0;
#else
  virtual size_t comparable_value() const = 0;
#endif
};

class FAMGVectorEntry
{
public:
  FAMGVectorEntry() : vecentry(NULL) {}
  FAMGVectorEntry( const FAMGVectorEntry & ve ) : vecentry(ve.GetPointer()->clone()) {}
  FAMGVectorEntry( FAMGVectorEntryRef* vep ) : vecentry(vep) {}
  ~FAMGVectorEntry() {
    delete vecentry;
  }

  FAMGVectorEntry& operator=( const FAMGVectorEntry & ve ) {if(this!=&ve) {delete vecentry; vecentry=ve.GetPointer()->clone();} return *this;}
  FAMGVectorEntry& operator++ () {++(*vecentry); return *this;}                 // prefix
  void operator++ (int) {++(*vecentry);}                                                                // postfix
  FAMGVectorEntry& operator-- () {--(*vecentry); return *this;}                 // prefix
  void operator-- (int) {--(*vecentry);}                                                                // postfix
  int operator==( const FAMGVectorEntry & ve ) const {return (*GetPointer()==*ve.GetPointer());}

  FAMGVectorEntryRef* GetPointer() const {
    return vecentry;
  }
  int GetIndex() const {
    return vecentry->GetIndex();
  }

private:
  FAMGVectorEntryRef* vecentry;
};

class FAMGGridVector
{
public:
  typedef class FAMGuVectorEntry VectorEntry;
  typedef class FAMGugVectorIter Iterator;

  virtual ~FAMGGridVector() {};                 // nothing to do

  virtual int is_valid( const FAMGVectorEntry& ve ) const = 0;
  virtual int is_end( const FAMGVectorEntry& ve ) const = 0;
  virtual int is_beforefirst( const FAMGVectorEntry& ve ) const = 0;
  virtual FAMGVectorEntry firstEntry() const = 0;
  virtual FAMGVectorEntry lastEntry() const = 0;
  virtual FAMGVectorEntry endEntry() const = 0;

  virtual int IsCG( const FAMGVectorEntry& ve ) const = 0;
  virtual int IsFG( const FAMGVectorEntry& ve ) const = 0;
  virtual void SetCG( const FAMGVectorEntry& ve ) const = 0;
  virtual void SetFG( const FAMGVectorEntry& ve ) const = 0;
};

class FAMGVector
{
public:
  friend class FAMGVectorIter;
  friend class FAMGVectorRevIter;

  FAMGVector(const FAMGGridVector & gridvec) : mygridvector(gridvec) {}
  virtual ~FAMGVector() {};                     // nothing to do
  virtual FAMGVector* create_new() const = 0;                           // create copy of my; incl. memory for data but without copying the data
  const FAMGGridVector& GetGridVector() const {
    return mygridvector;
  }

  virtual double& operator[] ( const FAMGVectorEntry & ve ) = 0;
  virtual double operator[] ( const FAMGVectorEntry & ve ) const = 0;
  virtual FAMGVector& operator=( const FAMGVector &v ) = 0;
  virtual FAMGVector& operator+=( const FAMGVector &v ) = 0;
  virtual FAMGVector& operator-=( const FAMGVector &v ) = 0;
  virtual double operator=(double val) = 0;
  virtual double operator*(const FAMGVector &v) = 0;                    // scalar product
  virtual FAMGVector& operator*=(double scale) = 0;

  int is_valid( const FAMGVectorEntry& ve ) const {
    return GetGridVector().is_valid(ve);
  }
  int is_end( const FAMGVectorEntry& ve ) const {
    return GetGridVector().is_end(ve);
  }
  int is_beforefirst( const FAMGVectorEntry& ve ) {
    return GetGridVector().is_beforefirst(ve);
  }
  FAMGVectorEntry firstEntry(const FAMGVectorEntry& ve) const {
    return GetGridVector().firstEntry();
  }
  FAMGVectorEntry lastEntry() const {
    return GetGridVector().lastEntry();
  }
  FAMGVectorEntry endEntry() const {
    return GetGridVector().endEntry();
  }

  int IsCG( const FAMGVectorEntry& ve ) const {
    return GetGridVector().IsCG(ve);
  }
  int IsFG( const FAMGVectorEntry& ve ) const {
    return GetGridVector().IsFG(ve);
  }
  void SetCG( const FAMGVectorEntry& ve ) const {
    GetGridVector().SetCG(ve);
  }
  void SetFG( const FAMGVectorEntry& ve ) const {
    GetGridVector().SetFG(ve);
  }

  virtual double norm() const = 0;
  virtual double sum() const = 0;
  virtual void AddScaledVec( double scale, const FAMGVector &source ) = 0;
  virtual void VecMinusMatVec( const FAMGVector &rhs, const FAMGMatrixAlg &mat, const FAMGVector &sol ) = 0;
  virtual void MatVec( const FAMGMatrixAlg &mat, const FAMGVector &source ) = 0;

  virtual void JacobiSmoother( const FAMGMatrixAlg &mat, const FAMGVector &def ) = 0;
  virtual void dampedJacobiSmoother( const FAMGMatrixAlg &mat, const FAMGVector &def ) = 0;
  virtual void FGSSmoother( const FAMGMatrixAlg &mat, FAMGVector &def ) = 0;
  virtual void BGSSmoother( const FAMGMatrixAlg &mat, FAMGVector &def ) = 0;
  virtual void SGSSmoother( const FAMGMatrixAlg &mat, FAMGVector &def ) = 0;
  virtual void JacobiSmoothFG( const FAMGMatrixAlg &mat, const FAMGVector &def ) = 0;
protected:
  const FAMGGridVector &mygridvector;
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

class FAMGMatrixEntryRef
{       // abstract base class
public:
  virtual FAMGMatrixEntryRef* clone() = 0;

  virtual FAMGMatrixEntryRef& operator++ () = 0;                // prefix
  FAMGMatrixEntryRef& operator++ (int) {FAMGMatrixEntryRef& tmp = *this; ++*this; return tmp;}              // postfix

  virtual FAMGVectorEntry dest() const = 0;

  virtual int is_strong() const = 0;
  virtual void set_strong(const int) = 0;
};

class FAMGMatrixEntry
{
public:
  FAMGMatrixEntry() : matentry(NULL) {}
  FAMGMatrixEntry( const FAMGMatrixEntry & me ) : matentry(me.GetPointer()->clone()) {}
  FAMGMatrixEntry( FAMGMatrixEntryRef* mep ) : matentry(mep) {}
  ~FAMGMatrixEntry() {
    delete matentry;
  }

  FAMGMatrixEntry& operator=( const FAMGMatrixEntry & me ) {if(this!=&me) {delete matentry; matentry=me.GetPointer()->clone();} return *this;}
  FAMGMatrixEntry& operator++ () {++(*matentry); return *this;}                 // prefix
  void operator++ (int) {++(*matentry);}                                                                // postfix

  FAMGVectorEntry dest() const {
    return matentry->dest();
  }
  FAMGMatrixEntryRef* GetPointer() const {
    return matentry;
  }

  int is_strong() const {
    return matentry->is_strong();
  }
  void set_strong(const int n) {
    matentry->set_strong(n);
  }

private:
  FAMGMatrixEntryRef* matentry;
};

class FAMGMatrixAlg
{
public:
  FAMGMatrixAlg( int nr_vecs, int nr_links ) : n(nr_vecs), nlinks(nr_links) {}
  virtual double operator[] ( const FAMGMatrixEntry & me ) const = 0;
  virtual double& operator[] ( const FAMGMatrixEntry & me ) = 0;
  virtual int is_valid( const FAMGVectorEntry& row_ve, const FAMGMatrixEntry& me ) const = 0;
  virtual int is_end( const FAMGVectorEntry& row_ve, const FAMGMatrixEntry& me ) const = 0;
  virtual FAMGMatrixEntry firstEntry( const FAMGVectorEntry& row_ve ) const = 0;
  virtual FAMGMatrixEntry endEntry( const FAMGVectorEntry& row_ve ) const  = 0;
  virtual double DiagValue( const FAMGVectorEntry& row_ve ) const  = 0;
  virtual double GetAdjData( const FAMGMatrixEntry& me ) const = 0;

  virtual int ConstructGalerkinMatrix( const FAMGGrid &fg ) = 0;
  virtual void MarkStrongLinks( const FAMGGrid &grid ) = 0;

  int &GetN() {
    return n;
  }
  int &GetNLinks() {
    return nlinks;
  }

private:
  int n;                        // number of vectors
  int nlinks;                   // number of links (==entries) in the matrix
};

class FAMGMatrixIter
{
public:
  FAMGMatrixIter( const FAMGMatrixAlg & m, const FAMGVectorEntry & row_ve) : myMatrix(m), myrow_ve(row_ve), current_me(m.firstEntry(row_ve)) {}

  int operator() ( FAMGMatrixEntry& me ) {
    int res = !myMatrix.is_end(myrow_ve,me=current_me); if(res) ++current_me;return res;
  }
  void reset() {
    current_me = myMatrix.firstEntry(myrow_ve);
  }

private:
  const FAMGMatrixAlg& myMatrix;
  const FAMGVectorEntry& myrow_ve;
  FAMGMatrixEntry current_me;
};

#endif  // ONLY_ONE_ALGEBRA_DS


//
// template functions to profit by special implementations
//		implementation in algebra.C
//

template<class VT>
void SetValue( VT &v, double val );

template<class VT>
void AddValue( VT &dest, const VT &source );

template<class VT>
void AddScaledValue( VT &dest, double scale, const VT &source );

template<class VT>
void SubtractValue( VT &dest, const VT &source );

template<class VT>
void CopyValue( VT &dest, const VT &source );

template<class VT>
double norm( const VT& v );

template<class VT>
double ScalProd( const VT& v, const VT& w );

template<class VT>
double sum( const VT& v );

template<class VT>
void Scale( VT& v, double scale );

template<class VT,class MT>
void VecMinusMatVec( VT &d, const VT &f, const MT &M, const VT &u );

template<class VT,class MT>
void MatVec( VT &dest, const MT &M, const VT &source );

template<class VT,class MT>
void JacobiSmoother( VT &sol, const MT &M, const VT &def );

template<class VT,class MT>
void dampedJacobiSmoother( VT &sol, const MT &M, const VT &def );

template<class VT,class MT>
void FGSSmoother( VT &sol, const MT &M, VT &def );

template<class VT,class MT>
void BGSSmoother( VT &sol, const MT &M, VT &def );

template<class VT,class MT>
void SGSSmoother( VT &sol, const MT &M, VT &def );

template<class VT,class MT>
void JacobiSmoothFG( VT &sol, const MT &M, const VT &def );

template<class MT>
int ConstructGalerkinMatrix( MT &Mcg, const FAMGGrid &fg );

template<class MT>
void MarkStrongLinks(const MT &A, const FAMGGrid &grid);

#ifdef ONLY_ONE_ALGEBRA_DS

#ifdef USE_UG_DS
#include "famg_onlyugalgebra.h"
#else
#error not yet implemented
#endif

#endif  // ONLY_ONE_ALGEBRA_DS

#ifdef __KCC
#include "famg_algebra.C"       // enable automatic template instantiation with explicit inclusion; remove famg_algebra.o from Makefile
#endif

#ifdef __SGI10__
#include "famg_algebra.C"       // enable automatic template instantiation with explicit inclusion; remove famg_algebra.o from Makefile
#endif

#endif
