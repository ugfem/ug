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

class FAMGVectorEntryRef
{       // abstract base class
public:
  virtual FAMGVectorEntryRef* clone() = NULL;

  virtual FAMGVectorEntryRef & operator++ () = NULL;                    // prefix
  FAMGVectorEntryRef & operator++ (int) {FAMGVectorEntryRef& tmp = *this; ++*this; return tmp;};              // postfix
  virtual FAMGVectorEntryRef & operator-- () = NULL;                    // prefix
  FAMGVectorEntryRef & operator-- (int) {FAMGVectorEntryRef& tmp = *this; --*this; return tmp;};              // postfix
};

class FAMGVectorEntry
{
public:
  FAMGVectorEntry() : vecentry(NULL) {};
  FAMGVectorEntry( const FAMGVectorEntry & ve ) : vecentry(ve.getpointer()->clone()) {};
  FAMGVectorEntry( FAMGVectorEntryRef* vep ) : vecentry(vep) {};
  ~FAMGVectorEntry() {
    delete vecentry;
  };

  FAMGVectorEntry& operator=( const FAMGVectorEntry & ve ) {if(this!=&ve) {delete vecentry; vecentry=ve.getpointer()->clone();} return *this;};
  FAMGVectorEntry& operator++ () {++(*vecentry); return *this;};                // prefix
  void operator++ (int) {++(*vecentry); };                                                              // postfix
  FAMGVectorEntry& operator-- () {--(*vecentry); return *this;};                // prefix
  void operator-- (int) {--(*vecentry); };                                                              // postfix

  FAMGVectorEntryRef* getpointer() const {
    return vecentry;
  };

private:
  FAMGVectorEntryRef* vecentry;
};

class FAMGVector
{
public:
  virtual double & operator[] (FAMGVectorEntry& ve) = NULL;
  virtual int is_valid(FAMGVectorEntry& ve) = NULL;
  virtual int is_end(FAMGVectorEntry& ve) = NULL;
  virtual FAMGVectorEntry firstEntry() = NULL;
  virtual const FAMGVectorEntry endEntry() = NULL;

};

class FAMGVectorIter
{
public:
  FAMGVectorIter( FAMGVector& v) : myvector(v), current_ve(v.firstEntry()) {};

  int operator() ( FAMGVectorEntry& ve ) {
    int res = !myvector.is_end(ve=current_ve); if(res) ++current_ve;return res;
  };

private:
  FAMGVector& myvector;
  FAMGVectorEntry current_ve;
};
#endif
