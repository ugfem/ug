// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <fstream.h>
#include <string.h>
#include <math.h>

#include <template.hh>
#include <myadt.hh>
#include <linalg/linalg.hh>


extern ofstream myerr;

double BaseVector :: shit = 0;

// %FD Constructs a vector of length zero
BaseVector :: BaseVector ()
{
  length = 0;
}

// %FD Constructs a vector of given length
BaseVector :: BaseVector (
  INDEX alength    // length of the vector
  )
{
  length = alength;
}

// %FD Sets length of the vector, old vector will be destroyed
void
BaseVector :: SetLength (
  INDEX alength          // new length of the vector
  )
{
  length = alength;
}

// %FD Changes length of the vector, old values stay valid
void
BaseVector :: ChangeLength (
  INDEX alength          // new length of the vector
  )
{
  length = alength;
}



// %FD { Write a vector with the help of the '<<' operator onto a stream }
ostream &    // stream for further use
operator<< (
  ostream & s,              // stream to write vector onto
  const BaseVector & v     // vector to write
  )
{
  return v.Print (s);
}


// %FD{ Divides every component of the vector by the scalar c.
//      The function checks for division by zero }
BaseVector &      // result vector
BaseVector :: operator/= (
  double c         // scalar to divide by
  )
{
  if (c)
    return (*this) *= (1/c);
  else
  {
    myerr << "operator/=: division by zero" << endl;
    return *this;
  }
}


// %FD Creates a copy of the object
BaseVector *      // pointer to the new vector
BaseVector :: Copy () const
{
  cerr << "Base_vector::Copy called" << endl << flush;
  return NULL;
}



TempVector :: ~TempVector ()
{
  delete vec;
}

TempVector BaseVector :: operator+ (const BaseVector & v2) const
{
  return (*Copy()) += v2;
}

TempVector BaseVector :: operator- (const BaseVector & v2) const
{
  return (*Copy()) -= v2;
}

TempVector BaseVector :: operator- () const
{
  return (*Copy()) *= -1;
}


TempVector BaseVector :: operator* (double scal) const
{
  return (*Copy()) *= scal;
}

TempVector BaseVector :: operator/ (double scal) const
{
  return (*Copy()) /= scal;
}


TempVector operator* (double scal, const BaseVector & v1)
{
  return v1 * scal;
}





BaseVector * TempVector :: Copy () const
{
  return vec->Copy();
}










double Vector :: shit = 0;



Vector :: Vector () : BaseVector()
{
  data = NULL;
}

Vector :: Vector (INDEX alength) : BaseVector (alength)
{

  //  (*testout) << "New Vector of size " << alength << endl;
  if (length)
  {
    data = new double[length];
    if (!data)
    {
      length = 0;
      myerr << "Vector not allocated" << endl;
    }
  }
  else
    data = NULL;
}

Vector :: Vector (const Vector & v2)
{
  length = v2.length;

  //  (*testout) << "New Vector (cc) of size " << length << endl;
  if (length)
  {
    data = new double[length];

    if (data)
    {
      memcpy (data, v2.data, length * sizeof (double));
    }
    else
    {
      length = 0;
      myerr << "Vector::Vector : Vector not allocated" << endl;
    }
  }
  else
    data = NULL;
}


Vector :: ~Vector ()
{
  //  (*testout) << "delete Vector of size " << length << endl;
  if (data) delete [] data;
}

void Vector :: SetLength (INDEX alength)
{
  if (length == alength) return;
  if (data) delete [] data;
  data = NULL;
  length = alength;

  if (length == 0) return;
  data = new double[length];
  if (!data)
  {
    length = 0;
    myerr << "Vector::SetLength: Vector not allocated" << endl;
  }
}

void Vector :: ChangeLength (INDEX alength)
{
  if (length == alength) return;

  if (alength == 0)
  {
    delete [] data;
    length = 0;
    return;
  }

  double * olddata = data;

  data = new double[alength];
  if (!data)
  {
    length = 0;
    myerr << "Vector::SetLength: Vector not allocated" << endl;
    delete [] olddata;
  }

  memcpy (data, olddata, min(alength, length));

  delete [] olddata;
  length = alength;
}

double & Vector :: operator() (INDEX i)
{
  if (i >= 1 && i <= length) return Elem(i);
  else myerr << "\nindex " << i << " out of range ("
             << 1 << "," << Length() << ")\n";
  return shit;
}

double Vector :: operator() (INDEX i) const
{
  if (i >= 1 && i <= length) return Get(i);
  else myerr << "\nindex " << i << " out of range ("
             << 1 << "," << Length() << ")\n" << flush;
  return shit;
}



double Vector :: SupNorm () const
{
  double sup = 0;

  for (INDEX i = 1; i <= Length(); i++)
    if (fabs (Get(i)) > sup)
      sup = fabs(Get(i));

  return sup;
}

double Vector :: L2Norm () const
{
  double sum = 0;

  for (INDEX i = 1; i <= Length(); i++)
    sum += Get(i) * Get(i);

  return sqrt (sum);
}

double Vector :: L1Norm () const
{
  double sum = 0;

  for (INDEX i = 1; i <= Length(); i++)
    sum += fabs (Get(i));

  return sum;
}

double Vector :: Max () const
{
  if (!Length()) return 0;
  double m = Get(1);
  for (INDEX i = 2; i <= Length(); i++)
    if (Get(i) > m) m = Get(i);
  return m;
}

double Vector :: Min () const
{
  if (!Length()) return 0;
  double m = Get(1);
  for (INDEX i = 2; i <= Length(); i++)
    if (Get(i) < m) m = Get(i);
  return m;
}


/*
   ostream & operator<<(ostream & s, const Vector & v)
   {
   int w = s.width();
   if (v.Length())
    {
    s.width(0);
    s << '(';
    for (INDEX i = 1; i < v.Length(); i++)
      {
      s.width(w);
      s << v.Get(i) << ",";
      if (i % 8 == 0) s << endl << ' ';
      }
    s.width(w);
    s << v.Get(v.Length()) << ')';
    }
   else
    s << "(Vector not allocated)";

   return s;
   }
 */

ostream & Vector :: Print (ostream & s) const
{
  int w = s.width();
  if (Length())
  {
    s.width(0);
    s << '(';
    for (INDEX i = 1; i < Length(); i++)
    {
      s.width(w);
      s << Get(i) << ",";
      if (i % 8 == 0) s << endl << ' ';
    }
    s.width(w);
    s << Get(Length()) << ')';
  }
  else
    s << "(Vector not allocated)";

  return s;
}



BaseVector & Vector :: operator+= (const BaseVector & v2)
{
  const Vector & hv2 = v2.CastToVector();

  if (Length() == hv2.Length())
    for (INDEX i = 1; i <= Length(); i++)
      Elem (i) += hv2.Get(i);
  else
    myerr << "operator+= illegal dimension" << endl;
  return *this;
}

BaseVector & Vector :: operator-= (const BaseVector & v2)
{
  const Vector & hv2 = v2.CastToVector();

  if (Length() == hv2.Length())
    for (INDEX i = 1; i <= Length(); i++)
      Elem (i) -= hv2.Get(i);
  else
    myerr << "operator-= illegal dimension" << endl;
  return *this;
}

BaseVector & Vector :: operator*= (double c)
{
  for (INDEX i = 1; i <= Length(); i++)
    Elem(i) *= c;
  return *this;
}



BaseVector & Vector :: Add (double scal, const BaseVector & v2)
{
  const Vector & hv2 = v2.CastToVector();

  if (Length() == hv2.Length())
  {
    double * p1 = data;
    double * p2 = hv2.data;

    for (INDEX i = Length(); i > 0; i--)
    {
      (*p1) += scal * (*p2);
      p1++; p2++;
    }
  }
  else
    myerr << "Vector::Add: illegal dimension" << endl;
  return *this;
}

BaseVector & Vector :: Add (double scal, const BaseVector & v2,
                            double scal3, const BaseVector & v3)
{
  const Vector & hv2 = v2.CastToVector();
  const Vector & hv3 = v3.CastToVector();

  if (Length() == hv2.Length())
  {
    double * p1 = data;
    double * p2 = hv2.data;
    double * p3 = hv3.data;

    for (INDEX i = Length(); i > 0; i--)
    {
      (*p1) += (scal * (*p2) + scal3 * (*p3));
      p1++; p2++; p3++;
    }
  }
  else
    myerr << "Vector::Add: illegal dimension" << endl;
  return *this;
}

BaseVector & Vector :: Set (double scal, const BaseVector & v2)
{
  const Vector & hv2 = v2.CastToVector();

  if (Length() == hv2.Length())
  {
    double * p1 = data;
    double * p2 = hv2.data;

    for (INDEX i = Length(); i > 0; i--)
    {
      (*p1) = scal * (*p2);
      p1++; p2++;
    }
  }
  else
    myerr << "Vector::Set: illegal dimension" << endl;
  return *this;
}


BaseVector & Vector :: Set (double scal , const BaseVector & v2,
                            double scal3, const BaseVector & v3)
{
  const Vector & hv2 = v2.CastToVector();
  const Vector & hv3 = v3.CastToVector();


  if (Length() == hv2.Length())
  {
    double * p1 = data;
    double * p2 = hv2.data;
    double * p3 = hv3.data;

    for (INDEX i = Length(); i > 0; i--)
    {
      (*p1) = scal * (*p2) + scal3 * (*p3);
      p1++; p2++; p3++;
    }
  }
  else
    myerr << "Vector::Set: illegal dimension" << endl;
  return *this;
}



double Vector :: operator* (const BaseVector & v2) const
{
  const Vector & hv2 = v2.CastToVector();

  double sum = 0;
  double * p1 = data;
  double * p2 = hv2.data;

  if (Length() == hv2.Length())
  {
    for (INDEX i = Length(); i > 0; i--)
    {
      sum += (*p1) * (*p2);
      p1++; p2++;
    }
  }
  else
    myerr << "Scalarproduct illegal dimension" << endl;
  return sum;
}

void Vector :: SetRandom ()
{
  INDEX i;
  for (i = 1; i <= Length(); i++)
    Elem(i) = 1.0 / double(i);
}


/*
   TempVector Vector :: operator- () const
   {
   Vector & sum = *(Vector*)Copy();

   if (sum.Length () == Length())
    {
    for (INDEX i = 1; i <= Length(); i++)
      sum.Set (i, Get(i));
    }
   else
    myerr << "operator+ (Vector, Vector): sum.Length() not ok" << endl;
   return sum;
   }
 */

BaseVector & Vector::operator= (const Vector & v2)
{
  SetLength (v2.Length());

  if (data == v2.data) return *this;

  if (v2.Length() == Length())
    memcpy (data, v2.data, sizeof (double) * Length());
  else
    myerr << "Vector::operator=: not allocated" << endl;

  return *this;
}

BaseVector & Vector::operator= (const BaseVector & v2)
{
  const Vector & hv2 = v2.CastToVector();

  SetLength (hv2.Length());

  if (data == hv2.data) return *this;

  if (hv2.Length() == Length())
    memcpy (data, hv2.data, sizeof (double) * Length());
  else
    myerr << "Vector::operator=: not allocated" << endl;

  return *this;
}


BaseVector & Vector::operator= (double scal)
{
  if (!Length()) myerr << "Vector::operator= (double) : data not allocated"
                       << endl;

  for (INDEX i = 1; i <= Length(); i++)
    Set (i, scal);

  return *this;
}


BaseVector * Vector :: Copy () const
{
  return new Vector (*this);
}


void Vector :: Swap (BaseVector & v2)
{
  Vector & hv2 = v2.CastToVector();
  swap (length, hv2.length);
  swap (data, hv2.data);
}
