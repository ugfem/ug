// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_DENSEMAT
#define FILE_DENSEMAT

/**************************************************************************/
/* File:   densemat.hh                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 94                                                    */
/**************************************************************************/

/*
   Data type dense matrix
 */



class DenseMatrix : public BaseMatrix
{
protected:
  double * data;

public:
  DenseMatrix ();
  DenseMatrix (INDEX h, INDEX w = 0);
  DenseMatrix (const DenseMatrix & m2);
  DenseMatrix (const BaseMatrix & m2);
  ~DenseMatrix ();

  virtual void SetSize (INDEX h, INDEX w = 0);


  virtual double & operator() (INDEX i, INDEX j);
  virtual double operator() (INDEX i, INDEX j) const;

  DenseMatrix & operator= (const BaseMatrix & m2);
  DenseMatrix & operator= (const DenseMatrix & m2);

  DenseMatrix & operator+= (const DenseMatrix & m2);
  DenseMatrix & operator-= (const DenseMatrix & m2);

  DenseMatrix & operator= (double v);
  DenseMatrix & operator*= (double v);

  virtual void Mult (const BaseVector & v, BaseVector & prod) const;
  virtual void MultTrans (const BaseVector & v, BaseVector & prod) const;
  virtual void Residuum (const BaseVector & x, const BaseVector & b, BaseVector & res) const;
  virtual double EvaluateBilinearform (const BaseVector & x) const;

  virtual BaseMatrix * Copy () const;

  virtual void MultElementMatrix (const ARRAY<INDEX> & pnum,
                                  const BaseVector & x, BaseVector & y);
  virtual void MultTransElementMatrix (const ARRAY<INDEX> & pnum,
                                       const BaseVector & x, BaseVector & y);


  double Det () const;

  friend DenseMatrix operator* (const DenseMatrix & m1, const DenseMatrix & m2);
  friend DenseMatrix operator+ (const DenseMatrix & m1, const DenseMatrix & m2);

  friend void Mult (const DenseMatrix & m1, const DenseMatrix & m2, DenseMatrix & m3);
  friend void CalcInverse (const DenseMatrix & m1, DenseMatrix & m2);
  friend void CalcAAt (const DenseMatrix & a, DenseMatrix & m2);
  friend void CalcAtA (const DenseMatrix & a, DenseMatrix & m2);
  friend void CalcABt (const DenseMatrix & a, const DenseMatrix & b, DenseMatrix & m2);
  friend void CalcAtB (const DenseMatrix & a, const DenseMatrix & b, DenseMatrix & m2);
  virtual void SolveDestroy (const Vector & b, Vector & x);



  double Get(INDEX i, INDEX j) const { return data[(i-1)*width+j-1]; }
  double Get(INDEX i) const { return data[i-1]; }
  void Set(INDEX i, INDEX j, double v) { data[(i-1)*width+j-1] = v; }
  double & Elem(INDEX i, INDEX j) { return data[(i-1)*width+j-1]; }
  const double & ConstElem(INDEX i, INDEX j) const { return data[(i-1)*width+j-1]; }
};



#endif
