// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   linalg.hh                                                    */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/

#ifndef FILE_LINALG
#define FILE_LINALG


class TEMP_VECTOR;

class VECTOR;

class SUPER_VECTOR;

class BITARRAY;



class BASE_VECTOR

{

protected:

  INDEX length;

  static double shit;



public:

  BASE_VECTOR ();

  BASE_VECTOR (INDEX alength);

  virtual ~BASE_VECTOR () { };



  virtual void SetLength (INDEX alength);

  virtual void ChangeLength (INDEX alength);

  INDEX Length () const { return length; }



  virtual BASE_VECTOR & operator= (const BASE_VECTOR & /* v2 */)
  { return *this; }

  virtual BASE_VECTOR & operator= (double /* scal */) { return *this; }



  virtual double & operator() (INDEX /* i */ ) { return shit; }

  virtual double operator() (INDEX /* i */) const { return 0; }



  virtual double SupNorm () const = 0;

  virtual double L2Norm () const = 0;

  virtual double L1Norm () const = 0;

  virtual double Min () const = 0;

  virtual double Max () const = 0;



  virtual BASE_VECTOR & operator+= (const BASE_VECTOR & v2) = 0;

  virtual BASE_VECTOR & operator-= (const BASE_VECTOR & v2) = 0;

  virtual BASE_VECTOR & operator*= (double scal) = 0;

  virtual BASE_VECTOR & operator/= (double scal);



  virtual TEMP_VECTOR operator+ (const BASE_VECTOR & v2) const;

  virtual TEMP_VECTOR operator- (const BASE_VECTOR & v2) const;

  virtual TEMP_VECTOR operator- () const;

  virtual double operator* (const BASE_VECTOR & v2) const = 0;

  TEMP_VECTOR operator* (double scal) const;

  virtual TEMP_VECTOR operator/ (double scal) const;

  friend TEMP_VECTOR operator* (double scal, const BASE_VECTOR & v1);



  virtual BASE_VECTOR & Add (double /* scal */, const BASE_VECTOR & /* v2 */)
  { return *this; }

  virtual BASE_VECTOR & Add (double /* scal */, const BASE_VECTOR & /* v2 */,

                             double /* scal3 */, const BASE_VECTOR & /* v3 */)
  { return *this; }

  virtual BASE_VECTOR & Set (double /* scal */, const BASE_VECTOR & /* v2 */)
  { return *this; }

  virtual BASE_VECTOR & Set (double /* scal */, const BASE_VECTOR & /* v2 */ ,

                             double /* scal3 */, const BASE_VECTOR & /* v3 */)
  { return *this; }



  virtual void SetRandom () { };



  friend ostream & operator<<(ostream & s, const BASE_VECTOR & v);

  virtual ostream & Print (ostream & s) const { return s; }



  virtual BASE_VECTOR * Copy () const;

  virtual VECTOR & CastToVector () { return *(VECTOR*)this; }

  virtual const VECTOR & CastToVector () const { return *(VECTOR*)this; }

  virtual SUPER_VECTOR & CastToSuperVector () { return *(SUPER_VECTOR*)this; }

  virtual const SUPER_VECTOR & CastToSuperVector () const { return *(SUPER_VECTOR*)this; }

};





class TEMP_VECTOR : public BASE_VECTOR

{

  BASE_VECTOR * vec;



public:

  TEMP_VECTOR (BASE_VECTOR & v1) { vec = & v1; }

  ~TEMP_VECTOR ();

  virtual VECTOR & CastToVector ()

  { return vec->CastToVector(); }

  virtual const VECTOR & CastToVector () const

  { return vec->CastToVector(); }

  virtual SUPER_VECTOR & CastToSuperVector ()

  { return vec->CastToSuperVector(); }

  virtual const SUPER_VECTOR & CastToSuperVector () const

  { return vec->CastToSuperVector(); }





  virtual BASE_VECTOR & operator+= (const BASE_VECTOR & /* v2 */)
  { return *this; }

  virtual BASE_VECTOR & operator-= (const BASE_VECTOR & /* v2 */)
  { return *this; }

  virtual BASE_VECTOR & operator*= (double /* scal */) { return *this; }

  virtual double operator* (const BASE_VECTOR & /* v2 */) const { return 0; }



  virtual double SupNorm () const { return vec->SupNorm(); }

  virtual double L2Norm () const { return vec->L2Norm(); }

  virtual double L1Norm () const { return vec->L1Norm(); }

  virtual double Min () const { return vec->Min(); }

  virtual double Max () const { return vec->Max(); }



  virtual void Swap (BASE_VECTOR &) { };

  virtual BASE_VECTOR * Copy () const;

};





class VECTOR : public BASE_VECTOR

{

  double * data;

  static double shit;



public:

  VECTOR ();

  VECTOR (INDEX alength);

  VECTOR (const VECTOR & v2);

  virtual ~VECTOR ();



  virtual void SetLength (INDEX alength);

  virtual void ChangeLength (INDEX alength);



  virtual BASE_VECTOR & operator= (const BASE_VECTOR & v2);

  virtual BASE_VECTOR & operator= (double scal);



  double & operator() (INDEX i);

  double operator() (INDEX i) const;



  virtual double SupNorm () const;

  virtual double L2Norm () const;

  virtual double L1Norm () const;

  virtual double Min () const;

  virtual double Max () const;



  virtual BASE_VECTOR & operator+= (const BASE_VECTOR & v2);

  virtual BASE_VECTOR & operator-= (const BASE_VECTOR & v2);

  virtual BASE_VECTOR & operator*= (double scal);

  //  virtual BASE_VECTOR & operator/= (double scal);



  virtual double operator* (const BASE_VECTOR & v2) const;





  virtual BASE_VECTOR & Add (double scal, const BASE_VECTOR & v2);

  virtual BASE_VECTOR & Add (double scal, const BASE_VECTOR & v2,

                             double scal3, const BASE_VECTOR & v3);

  virtual BASE_VECTOR & Set (double scal, const BASE_VECTOR & v2);

  virtual BASE_VECTOR & Set (double scal , const BASE_VECTOR & v2,

                             double scal3, const BASE_VECTOR & v3);



  virtual void SetRandom ();





  /*

     friend double operator* (const VECTOR & v1, const VECTOR & v2);

     friend VECTOR operator+ (const VECTOR & v1, const VECTOR & v2);

     friend VECTOR operator- (const VECTOR & v1, const VECTOR & v2);

     friend VECTOR operator- (const VECTOR & v1);

     friend VECTOR operator* (double scal, const VECTOR & v1);

   */

  virtual ostream & Print (ostream & s) const;

  virtual BASE_VECTOR * Copy () const;

  virtual void Swap (BASE_VECTOR &);



  const double & Get (INDEX i) const { return data[i-1]; }

  void Set (INDEX i, double v) { data[i-1] = v; }

  double & Elem (INDEX i) { return data[i-1]; }

};









class MATRIX;



class BASE_MATRIX

{

protected:

  INDEX height, width;

  int symmetric;

  static double shit;



public:

  BASE_MATRIX ();

  BASE_MATRIX (INDEX h, INDEX w = 0);

  virtual ~BASE_MATRIX () { };



  INDEX Width () const { return width; }

  INDEX Height () const { return height; }

  int Symmetric () const { return symmetric; }



  virtual void SetSize (INDEX h, INDEX w = 0);

  virtual void SetSymmetric (int sym = 1);



  virtual double & operator() (INDEX i, INDEX j);

  virtual double operator() (INDEX i, INDEX j) const;



  friend ostream & operator<<(ostream & s, const BASE_MATRIX & m);

  virtual ostream & Print (ostream & s) const;



  TEMP_VECTOR operator* (const BASE_VECTOR & v) const;



  friend Mult (const BASE_MATRIX & m1, const BASE_MATRIX & m2,

               MATRIX & m3);

  friend Add (const BASE_MATRIX & m1, const BASE_MATRIX & m2,

              MATRIX & m3);

  friend void Transpose (const BASE_MATRIX & m1, MATRIX & m2);

  friend MATRIX operator* (const BASE_MATRIX & m1, const BASE_MATRIX & m2);

  friend MATRIX operator+ (const BASE_MATRIX & m1, const BASE_MATRIX & m2);

  friend MATRIX Transpose (const BASE_MATRIX & m1);



  virtual void Mult (const BASE_VECTOR & v, BASE_VECTOR & prod) const;

  virtual void MultTrans (const BASE_VECTOR & v, BASE_VECTOR & prod) const;

  virtual void Residuum (const BASE_VECTOR & x, const BASE_VECTOR & b, BASE_VECTOR & res) const;

  //  virtual double EvaluateBilinearform (const BASE_VECTOR & x);



  virtual BASE_MATRIX * Copy () const;

  virtual BASE_VECTOR * CreateVector () const;



  virtual void AddElementMatrix (const ARRAY<INDEX> & /* pnum */,
                                 const BASE_MATRIX & /* elemmat */) { };





  virtual void SolveDestroy (const VECTOR & b, VECTOR & x);

  void Solve (const VECTOR & b, VECTOR & x) const;

  VECTOR SolveDestroy (const VECTOR & b) const;

  VECTOR Solve (const VECTOR & b) const;

  void SolvePI (const VECTOR & b, VECTOR & x) const;

  virtual void LU_Decomposition (MATRIX & l, MATRIX & u) const;

  virtual void Gerschgorin (double & min, double & max) const;



  VECTOR Jacobi (const VECTOR & b, float eps = 1E-6) const;

  void Jacobi (const VECTOR & b, VECTOR & x, float eps = 1E-6) const;



  VECTOR CG (const VECTOR & b, float eps = 1E-6) const;

  void CG (const BASE_VECTOR & b, BASE_VECTOR & x, float eps = 1E-6) const;

  VECTOR PCG (const BASE_MATRIX & pre, const VECTOR & b, float eps = 1E-6) const;

  void PCG (const BASE_MATRIX & pre, const VECTOR & b, VECTOR & x, float eps = 1E-6) const;

  void PrecondRichardson (const BASE_MATRIX & pre, const BASE_VECTOR & b,

                          BASE_VECTOR & x, double tau, double eps) const;



  void LargestEigenValue (double & val, VECTOR & x) const;

  void EigenSystem (ARRAY<double> & val, ARRAY<VECTOR> & x) const;

  void ConstraintEigenSystem (ARRAY<double> & val, ARRAY<VECTOR> & x,

                              const ARRAY<VECTOR> & cons) const;

  double SpectralNorm () const;

};









class MATRIX : public BASE_MATRIX

{

protected:

  double * data;



public:

  MATRIX ();

  MATRIX (INDEX h, INDEX w = 0);

  MATRIX (const MATRIX & m2);

  MATRIX (const BASE_MATRIX & m2);

  ~MATRIX ();



  virtual void SetSize (INDEX h, INDEX w = 0);





  virtual double & operator() (INDEX i, INDEX j);

  virtual double operator() (INDEX i, INDEX j) const;



  MATRIX & operator= (const BASE_MATRIX & m2);

  MATRIX & operator= (const MATRIX & m2);



  MATRIX & operator+= (const MATRIX & m2);

  MATRIX & operator-= (const MATRIX & m2);



  MATRIX & operator= (double v);

  MATRIX & operator*= (double v);



  virtual void Mult (const BASE_VECTOR & v, BASE_VECTOR & prod) const;

  virtual void MultTrans (const BASE_VECTOR & v, BASE_VECTOR & prod) const;

  virtual void Residuum (const BASE_VECTOR & x, const BASE_VECTOR & b, BASE_VECTOR & res) const;

  virtual double EvaluateBilinearform (const BASE_VECTOR & x) const;



  virtual BASE_MATRIX * Copy () const;





  double Det () const;



  friend MATRIX operator* (const MATRIX & m1, const MATRIX & m2);

  friend MATRIX operator+ (const MATRIX & m1, const MATRIX & m2);



  friend void Mult (const MATRIX & m1, const MATRIX & m2, MATRIX & m3);

  friend void CalcInverse (const MATRIX & m1, MATRIX & m2);

  friend void CalcAAt (const MATRIX & a, MATRIX & m2);

  friend void CalcAtA (const MATRIX & a, MATRIX & m2);

  friend void CalcABt (const MATRIX & a, const MATRIX & b, MATRIX & m2);

  friend void CalcAtB (const MATRIX & a, const MATRIX & b, MATRIX & m2);

  virtual void SolveDestroy (const VECTOR & b, VECTOR & x);





  double Get(INDEX i, INDEX j) const { return data[(i-1)*width+j-1]; }

  double Get(INDEX i) const { return data[i-1]; }

  void Set(INDEX i, INDEX j, double v) { data[(i-1)*width+j-1] = v; }

  double & Elem(INDEX i, INDEX j) { return data[(i-1)*width+j-1]; }

  const double & ConstElem(INDEX i, INDEX j) const { return data[(i-1)*width+j-1]; }

};







class LOWER_SPARSE_MATRIX;

class UPPER_SPARSE_MATRIX;



class SPARSE_MATRIX : public BASE_MATRIX

{

protected:

  struct colstruct { INDEX colnr; double data; };

  struct linestruct { int size; int maxsize; colstruct * col; };



  ARRAY<linestruct> lins;



public:

  SPARSE_MATRIX ();

  SPARSE_MATRIX (INDEX h, INDEX w = 0);

  SPARSE_MATRIX (const SPARSE_MATRIX & m2);

  ~SPARSE_MATRIX ();



  virtual void SetSize (INDEX h, INDEX w = 0);

  virtual void SetSymmetric (int sym = 1);

  virtual void ChangeSize (INDEX h, INDEX w = 0);

  void DeleteElements ();



  virtual double & operator() (INDEX i, INDEX j);

  virtual double operator() (INDEX i, INDEX j) const;



  SPARSE_MATRIX & operator= (const SPARSE_MATRIX & m2);

  SPARSE_MATRIX & operator*= (double v);



  virtual void Mult (const BASE_VECTOR & v, BASE_VECTOR & prod) const;

  virtual void MultTrans (const BASE_VECTOR & v, BASE_VECTOR & prod) const;

  virtual void Residuum (const BASE_VECTOR & x, const BASE_VECTOR & b, BASE_VECTOR & res) const;



  virtual BASE_MATRIX * Copy () const;





  void GSStep (const VECTOR & b, VECTOR & x, double dump) const;

  void GSStepBack (const VECTOR & b, VECTOR & x, double dump) const;



  void GSStepInner (const VECTOR & b, VECTOR & x, double dump,

                    const BITARRAY & inner) const;

  void GSStepBackInner (const VECTOR & b, VECTOR & x, double dump,

                        const BITARRAY & inner) const;



  void GSStepToInner (const VECTOR & b, VECTOR & x, double dump,

                      const BITARRAY & inner) const;

  void GSStepBackToInner (const VECTOR & b, VECTOR & x, double dump,

                          const BITARRAY & inner) const;



  void ILU_Decomposition (LOWER_SPARSE_MATRIX & l, UPPER_SPARSE_MATRIX & u) const;

  VECTOR ILU (const VECTOR & v, float eps = 1E-6) const;

  VECTOR ILU_PCG (const VECTOR & v, float eps = 1E-6) const;

  void ILU_PCG (const VECTOR & v, VECTOR & x, const LOWER_SPARSE_MATRIX & l,

                const UPPER_SPARSE_MATRIX & u, float eps = 1E-6) const;

  virtual void Gerschgorin (double & min, double & max) const;

  friend void Transpose (const SPARSE_MATRIX & m1, SPARSE_MATRIX & m2);



  virtual void Solve (const VECTOR & b, VECTOR & x) const;

  friend SPARSE_MATRIX operator* (const SPARSE_MATRIX & m1, const SPARSE_MATRIX & m2);

  SPARSE_MATRIX & operator+= (const SPARSE_MATRIX & m2);

  SPARSE_MATRIX & operator*= (const SPARSE_MATRIX & m2);

  virtual ostream & Print (ostream & s) const;



  virtual void AddElementMatrix (const ARRAY<INDEX> & pnum, const BASE_MATRIX & elemmat);

  void AddRowMatrix (INDEX row, const SPARSE_MATRIX & m2);



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







class LOWER_SPARSE_MATRIX : public SPARSE_MATRIX

{

public:

  LOWER_SPARSE_MATRIX () : SPARSE_MATRIX () { };

  LOWER_SPARSE_MATRIX (INDEX n) : SPARSE_MATRIX (n, n) { };

  ~LOWER_SPARSE_MATRIX () { };



  virtual void Solve (const VECTOR & b, VECTOR & x) const;

};





class UPPER_SPARSE_MATRIX : public SPARSE_MATRIX

{

public:

  UPPER_SPARSE_MATRIX () : SPARSE_MATRIX () { };

  UPPER_SPARSE_MATRIX (INDEX n) : SPARSE_MATRIX (n, n) { };

  ~UPPER_SPARSE_MATRIX () { };



  virtual void Solve (const VECTOR & b, VECTOR & x) const;

};





extern void Lanczos (const BASE_MATRIX & a);

extern void GeneralLanczos (const BASE_MATRIX & a, const BASE_MATRIX & b);

extern void PrecondLanczos (const BASE_MATRIX & a, const BASE_MATRIX & b);

#endif
