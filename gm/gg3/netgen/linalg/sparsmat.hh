// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_SPARSMAT
#define FILE_SPARSMAT

/**************************************************************************/
/* File:   sparsmat.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 94                                                    */
/**************************************************************************/

/*
   Data type sparse matrix
 */


class LOWER_SparseMatrix;
class UPPER_SparseMatrix;

class SparseMatrix : public BaseMatrix
{
public:
  class colstruct { public: INDEX colnr; double data; };
  class linestruct { public: int size; int maxsize; colstruct * col; };

  ARRAY<linestruct> lins;

public:
  SparseMatrix ();
  SparseMatrix (INDEX h, INDEX w = 0);
  SparseMatrix (const SparseMatrix & m2);
  ~SparseMatrix ();

  virtual void SetSize (INDEX h, INDEX w = 0);
  virtual void SetSymmetric (int sym = 1);
  virtual void ChangeSize (INDEX h, INDEX w = 0);
  void DeleteElements ();

  virtual double & operator() (INDEX i, INDEX j);
  virtual double operator() (INDEX i, INDEX j) const;

  SparseMatrix & operator= (const SparseMatrix & m2);
  SparseMatrix & operator*= (double v);

  virtual void Mult (const BaseVector & v, BaseVector & prod) const;
  virtual void MultTrans (const BaseVector & v, BaseVector & prod) const;
  virtual void Residuum (const BaseVector & x, const BaseVector & b, BaseVector & res) const;
  virtual void ResiduumTrans (const BaseVector & x, const BaseVector & b, BaseVector & res) const;

  virtual BaseMatrix * Copy () const;


  void GSStep (const Vector & b, Vector & x, double dump) const;
  void GSStepBack (const Vector & b, Vector & x, double dump) const;

  void GSStepInner (const Vector & b, Vector & x, double dump,
                    const BitArray & inner) const;
  void GSStepBackInner (const Vector & b, Vector & x, double dump,
                        const BitArray & inner) const;

  void GSStepToInner (const Vector & b, Vector & x, double dump,
                      const BitArray & inner) const;
  void GSStepBackToInner (const Vector & b, Vector & x, double dump,
                          const BitArray & inner) const;

  friend void Transpose (const SparseMatrix & m1, SparseMatrix & m2);

  virtual void Solve (const Vector & b, Vector & x) const;
  friend SparseMatrix operator* (const SparseMatrix & m1, const SparseMatrix & m2);
  SparseMatrix & operator+= (const SparseMatrix & m2);
  SparseMatrix & operator*= (const SparseMatrix & m2);
  virtual ostream & Print (ostream & s) const;

  virtual void AddElementMatrix (const ARRAY<INDEX> & pnum, const BaseMatrix & elemmat);
  void AddRowMatrix (INDEX row, const SparseMatrix & m2);
  double RowTimesVector (INDEX i, const Vector & v) const;
  void AddRowToVector (INDEX i, double s, Vector & v) const;

  int ElementsInLine (INDEX i) const { return lins[i].size; }
  INDEX GetIndex (INDEX i, int nr) const { return lins[i].col[nr-1].colnr; }
  double GetData (INDEX i, int nr) const { return lins[i].col[nr-1].data; }
  double & GetDataRef (INDEX i, int nr) const { return lins[i].col[nr-1].data; }
  void Delete (INDEX i, int nr);

  void DeleteElem (INDEX i, INDEX j);
  double Get(INDEX i, INDEX j) const;
  void Set(INDEX i, INDEX j, double v);
  double & Elem(INDEX i, INDEX j);
  char Used (INDEX i, INDEX j) const;

  void SetLineAllocSize (INDEX i, int j);

private:
  colstruct * NewColStruct (int s);
  void DeleteColStruct (colstruct * cs, int s);
};




#endif
