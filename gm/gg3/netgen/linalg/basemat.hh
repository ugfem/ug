// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_BASEMAT
#define FILE_BASEMAT

/**************************************************************************/
/* File:   basemat.hh                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 94                                                    */
/**************************************************************************/

/*
   Base type for linear operator
 */

class DenseMatrix;

class BaseMatrix
{
protected:
  INDEX height, width;
  BOOL symmetric;
  static double shit;

public:
  BaseMatrix ();
  BaseMatrix (INDEX h, INDEX w = 0);
  virtual ~BaseMatrix () { };

  INDEX Width () const { return width; }
  INDEX Height () const { return height; }
  BOOL Symmetric () const { return symmetric; }

  virtual void SetSize (INDEX h, INDEX w = 0);
  virtual void SetSymmetric (int sym = 1);

  virtual double & operator() (INDEX i, INDEX j);
  virtual double operator() (INDEX i, INDEX j) const;

  friend ostream & operator<<(ostream & s, const BaseMatrix & m);
  virtual ostream & Print (ostream & s) const;

  TempVector operator* (const BaseVector & v) const;

  /*
     friend Mult (const BaseMatrix & m1, const BaseMatrix & m2,
        DenseMatrix & m3);
     friend Add (const BaseMatrix & m1, const BaseMatrix & m2,
        DenseMatrix & m3);
     friend void Transpose (const BaseMatrix & m1, DenseMatrix & m2);
     friend DenseMatrix operator* (const BaseMatrix & m1, const BaseMatrix & m2);
     friend DenseMatrix operator+ (const BaseMatrix & m1, const BaseMatrix & m2);
     friend DenseMatrix Transpose (const BaseMatrix & m1);
   */

  virtual void Mult (const BaseVector & v, BaseVector & prod) const;
  virtual void MultTrans (const BaseVector & v, BaseVector & prod) const;
  virtual void Residuum (const BaseVector & x, const BaseVector & b, BaseVector & res) const;
  virtual void ResiduumTrans (const BaseVector & x, const BaseVector & b, BaseVector & res) const;
  //  virtual double EvaluateBilinearform (const BaseVector & x);

  virtual BaseMatrix * Copy () const;
  virtual BaseVector * CreateVector () const;

  virtual void AddElementMatrix (const ARRAY<INDEX> & /* pnum */,
                                 const BaseMatrix & /* elemmat */) { };
  virtual void MultElementMatrix (const ARRAY<INDEX> & /* pnum */,
                                  const BaseVector & /* x */, BaseVector & /* y */) { };
  virtual void MultTransElementMatrix (const ARRAY<INDEX> & /* pnum */,
                                       const BaseVector & /* x */, BaseVector & /* y */) { };


  virtual void SolveDestroy (const Vector & b, Vector & x);
  void Solve (const Vector & b, Vector & x) const;
  Vector SolveDestroy (const Vector & b) const;
  Vector Solve (const Vector & b) const;
  virtual void LU_Decomposition (DenseMatrix & l, DenseMatrix & u) const;

  Vector Jacobi (const Vector & b, float eps = 1E-6) const;
  void Jacobi (const Vector & b, Vector & x, float eps = 1E-6) const;

  void PrecondRichardson (const BaseMatrix & pre, const BaseVector & b,
                          BaseVector & x, double tau, double eps) const;
};


#endif
