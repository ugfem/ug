// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_VECTOR
#define FILE_VECTOR

/**************************************************************************/
/* File:   vector.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 94                                                    */
/**************************************************************************/

/*
   Data types for vectors
 */

class TempVector;
class Vector;
class BlockVector;

class BaseVector
{
protected:
  INDEX length;
  static double shit;

public:
  BaseVector ();
  BaseVector (INDEX alength);
  virtual ~BaseVector () { };

  virtual void SetLength (INDEX alength);
  virtual void ChangeLength (INDEX alength);
  INDEX Length () const { return length; }

  virtual BaseVector & operator= (const BaseVector & /* v2 */) { return *this; }
  virtual BaseVector & operator= (double /* scal */) { return *this; }

  virtual double & operator() (INDEX /* i */) { return shit; }
  virtual double operator() (INDEX /* i */) const { return 0; }

  virtual double SupNorm () const = 0;
  virtual double L2Norm () const = 0;
  virtual double L1Norm () const = 0;
  virtual double Min () const = 0;
  virtual double Max () const = 0;

  virtual BaseVector & operator+= (const BaseVector & v2) = 0;
  virtual BaseVector & operator-= (const BaseVector & v2) = 0;
  virtual BaseVector & operator*= (double scal) = 0;
  virtual BaseVector & operator/= (double scal);

  virtual TempVector operator+ (const BaseVector & v2) const;
  virtual TempVector operator- (const BaseVector & v2) const;
  virtual TempVector operator- () const;
  virtual double operator* (const BaseVector & v2) const = 0;
  TempVector operator* (double scal) const;
  virtual TempVector operator/ (double scal) const;
  friend TempVector operator* (double scal, const BaseVector & v1);

  virtual BaseVector & Add (double /* scal */, const BaseVector & /* v2 */) { return *this; }
  virtual BaseVector & Add (double /* scal */, const BaseVector & /* v2 */,
                            double /* scal3 */, const BaseVector & /* v3 */) { return *this; }
  virtual BaseVector & Set (double /* scal */, const BaseVector & /* v2 */) { return *this; }
  virtual BaseVector & Set (double /* scal */, const BaseVector & /* v2 */,
                            double /* scal3 */, const BaseVector & /* v3 */) { return *this; }

  virtual void SetRandom () { };

  friend ostream & operator<<(ostream & s, const BaseVector & v);
  virtual ostream & Print (ostream & s) const { return s; }

  virtual BaseVector * Copy () const;
  virtual BOOL IsBlockVector () const { return 0; }
  virtual Vector & CastToVector () { return *(Vector*)this; }
  virtual const Vector & CastToVector () const { return *(Vector*)this; }
  virtual BlockVector & CastToSuperVector () { return *(BlockVector*)this; }
  virtual const BlockVector & CastToSuperVector () const { return *(BlockVector*)this; }
};


class TempVector : public BaseVector
{
  BaseVector * vec;

public:
  TempVector (BaseVector & v1) { vec = & v1; }
  ~TempVector ();
  virtual Vector & CastToVector ()
  { return vec->CastToVector(); }
  virtual const Vector & CastToVector () const
  { return vec->CastToVector(); }
  virtual BlockVector & CastToSuperVector ()
  { return vec->CastToSuperVector(); }
  virtual const BlockVector & CastToSuperVector () const
  { return vec->CastToSuperVector(); }


  virtual BaseVector & operator+= (const BaseVector & /* v2 */) { return *this; }
  virtual BaseVector & operator-= (const BaseVector & /* v2 */) { return *this; }
  virtual BaseVector & operator*= (double /* scal */) { return *this; }
  virtual double operator* (const BaseVector & /* v2 */) const { return 0; }

  virtual double SupNorm () const { return vec->SupNorm(); }
  virtual double L2Norm () const { return vec->L2Norm(); }
  virtual double L1Norm () const { return vec->L1Norm(); }
  virtual double Min () const { return vec->Min(); }
  virtual double Max () const { return vec->Max(); }

  virtual void Swap (BaseVector &) { };
  virtual BaseVector * Copy () const;
};


class Vector : public BaseVector
{
  double * data;
  static double shit;

public:
  Vector ();
  Vector (INDEX alength);
  Vector (const Vector & v2);
  virtual ~Vector ();

  virtual void SetLength (INDEX alength);
  virtual void ChangeLength (INDEX alength);

  virtual BaseVector & operator= (const BaseVector & v2);
  virtual BaseVector & operator= (const Vector & v2);
  virtual BaseVector & operator= (double scal);

  double & operator() (INDEX i);
  double operator() (INDEX i) const;

  virtual double SupNorm () const;
  virtual double L2Norm () const;
  virtual double L1Norm () const;
  virtual double Min () const;
  virtual double Max () const;

  virtual BaseVector & operator+= (const BaseVector & v2);
  virtual BaseVector & operator-= (const BaseVector & v2);
  virtual BaseVector & operator*= (double scal);
  //  virtual BaseVector & operator/= (double scal);

  virtual double operator* (const BaseVector & v2) const;


  virtual BaseVector & Add (double scal, const BaseVector & v2);
  virtual BaseVector & Add (double scal, const BaseVector & v2,
                            double scal3, const BaseVector & v3);
  virtual BaseVector & Set (double scal, const BaseVector & v2);
  virtual BaseVector & Set (double scal , const BaseVector & v2,
                            double scal3, const BaseVector & v3);

  virtual void SetRandom ();

  virtual ostream & Print (ostream & s) const;
  virtual BaseVector * Copy () const;
  virtual void Swap (BaseVector &);

  const double & Get (INDEX i) const { return data[i-1]; }
  void Set (INDEX i, double v) { data[i-1] = v; }
  double & Elem (INDEX i) { return data[i-1]; }
};

#endif
