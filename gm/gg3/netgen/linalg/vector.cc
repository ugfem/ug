// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>

#include <fstream.h>

#include <string.h>

#include <math.h>



#include <template.hh>

#include <array.hh>

#include <linalg/linalg.hh>





extern ofstream myerr;







double BASE_VECTOR :: shit = 0;



// %FD Constructs a vector of length zero

BASE_VECTOR :: BASE_VECTOR ()

{

  length = 0;

}



// %FD Constructs a vector of given length

BASE_VECTOR :: BASE_VECTOR (

  INDEX alength    // length of the vector

  )

{

  length = alength;

}



// %FD Sets length of the vector, old vector will be destroyed

void

BASE_VECTOR :: SetLength (

  INDEX alength          // new length of the vector

  )

{

  length = alength;

}



// %FD Changes length of the vector, old values stay valid

void

BASE_VECTOR :: ChangeLength (

  INDEX alength          // new length of the vector

  )

{

  length = alength;

}







// %FD { Write a vector with the help of the '<<' operator onto a stream }

ostream &    // stream for further use

operator<< (

  ostream & s,              // stream to write vector onto

  const BASE_VECTOR & v     // vector to write

  )

{

  return v.Print (s);

}





// %FD{ Divides every component of the vector by the scalar c.

//      The function checks for division by zero }

BASE_VECTOR &      // result vector

BASE_VECTOR :: operator/= (

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

BASE_VECTOR *      // pointer to the new vector

BASE_VECTOR :: Copy () const

{

  cerr << "Base_vector::Copy called" << endl << flush;

  return NULL;

}







TEMP_VECTOR :: ~TEMP_VECTOR ()

{

  delete vec;

}



TEMP_VECTOR BASE_VECTOR :: operator+ (const BASE_VECTOR & v2) const

{

  return (*Copy()) += v2;

}



TEMP_VECTOR BASE_VECTOR :: operator- (const BASE_VECTOR & v2) const

{

  return (*Copy()) -= v2;

}



TEMP_VECTOR BASE_VECTOR :: operator- () const

{

  return (*Copy()) *= -1;

}





TEMP_VECTOR BASE_VECTOR :: operator* (double scal) const

{

  return (*Copy()) *= scal;

}



TEMP_VECTOR BASE_VECTOR :: operator/ (double scal) const

{

  return (*Copy()) /= scal;

}





TEMP_VECTOR operator* (double scal, const BASE_VECTOR & v1)

{

  return v1 * scal;

}











BASE_VECTOR * TEMP_VECTOR :: Copy () const

{

  return vec->Copy();

}





















double VECTOR :: shit = 0;







VECTOR :: VECTOR () : BASE_VECTOR()

{

  data = NULL;

}



VECTOR :: VECTOR (INDEX alength) : BASE_VECTOR (alength)

{

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



VECTOR :: VECTOR (const VECTOR & v2)

{

  length = v2.length;



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

      myerr << "VECTOR::VECTOR : Vector not allocated" << endl;

    }

  }

  else

    data = NULL;

}





VECTOR :: ~VECTOR ()

{

  if (data) delete [] data;

}



void VECTOR :: SetLength (INDEX alength)

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

    myerr << "VECTOR::SetLength: Vector not allocated" << endl;

  }

}



void VECTOR :: ChangeLength (INDEX alength)

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

    myerr << "VECTOR::SetLength: Vector not allocated" << endl;

    delete [] olddata;

  }



  memcpy (data, olddata, min(alength, length));



  delete [] olddata;

  length = alength;

}



double & VECTOR :: operator() (INDEX i)

{

  if (i >= 1 && i <= length) return Elem(i);

  else myerr << "\nindex " << i << " out of range ("

             << 1 << "," << Length() << ")\n";

  return shit;

}



double VECTOR :: operator() (INDEX i) const

{

  if (i >= 1 && i <= length) return Get(i);

  else myerr << "\nindex " << i << " out of range ("

             << 1 << "," << Length() << ")\n" << flush;

  return shit;

}







double VECTOR :: SupNorm () const

{

  double sup = 0;



  for (INDEX i = 1; i <= Length(); i++)

    if (fabs (Get(i)) > sup)

      sup = fabs(Get(i));



  return sup;

}



double VECTOR :: L2Norm () const

{

  double sum = 0;



  for (INDEX i = 1; i <= Length(); i++)

    sum += Get(i) * Get(i);



  return sqrt (sum);

}



double VECTOR :: L1Norm () const

{

  double sum = 0;



  for (INDEX i = 1; i <= Length(); i++)

    sum += fabs (Get(i));



  return sum;

}



double VECTOR :: Max () const

{

  if (!Length()) return 0;

  double m = Get(1);

  for (INDEX i = 2; i <= Length(); i++)

    if (Get(i) > m) m = Get(i);

  return m;

}



double VECTOR :: Min () const

{

  if (!Length()) return 0;

  double m = Get(1);

  for (INDEX i = 2; i <= Length(); i++)

    if (Get(i) < m) m = Get(i);

  return m;

}





/*

   ostream & operator<<(ostream & s, const VECTOR & v)

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



ostream & VECTOR :: Print (ostream & s) const

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







BASE_VECTOR & VECTOR :: operator+= (const BASE_VECTOR & v2)

{

  const VECTOR & hv2 = v2.CastToVector();



  if (Length() == hv2.Length())

    for (INDEX i = 1; i <= Length(); i++)

      Elem (i) += hv2.Get(i);

  else

    myerr << "operator+= illegal dimension" << endl;

  return *this;

}



BASE_VECTOR & VECTOR :: operator-= (const BASE_VECTOR & v2)

{

  const VECTOR & hv2 = v2.CastToVector();



  if (Length() == hv2.Length())

    for (INDEX i = 1; i <= Length(); i++)

      Elem (i) -= hv2.Get(i);

  else

    myerr << "operator-= illegal dimension" << endl;

  return *this;

}



BASE_VECTOR & VECTOR :: operator*= (double c)

{

  for (INDEX i = 1; i <= Length(); i++)

    Elem(i) *= c;

  return *this;

}







BASE_VECTOR & VECTOR :: Add (double scal, const BASE_VECTOR & v2)

{

  const VECTOR & hv2 = v2.CastToVector();



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

    myerr << "VECTOR::Add: illegal dimension" << endl;

  return *this;

}



BASE_VECTOR & VECTOR :: Add (double scal, const BASE_VECTOR & v2,

                             double scal3, const BASE_VECTOR & v3)

{

  const VECTOR & hv2 = v2.CastToVector();

  const VECTOR & hv3 = v3.CastToVector();



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

    myerr << "VECTOR::Add: illegal dimension" << endl;

  return *this;

}



BASE_VECTOR & VECTOR :: Set (double scal, const BASE_VECTOR & v2)

{

  const VECTOR & hv2 = v2.CastToVector();



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

    myerr << "VECTOR::Set: illegal dimension" << endl;

  return *this;

}





BASE_VECTOR & VECTOR :: Set (double scal , const BASE_VECTOR & v2,

                             double scal3, const BASE_VECTOR & v3)

{

  const VECTOR & hv2 = v2.CastToVector();

  const VECTOR & hv3 = v3.CastToVector();





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

    myerr << "VECTOR::Set: illegal dimension" << endl;

  return *this;

}







double VECTOR :: operator* (const BASE_VECTOR & v2) const

{

  const VECTOR & hv2 = v2.CastToVector();



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



void VECTOR :: SetRandom ()

{

  INDEX i;

  for (i = 1; i <= Length(); i++)

    Elem(i) = 1.0 / double(i);

}





/*

   TEMP_VECTOR VECTOR :: operator- () const

   {

   VECTOR & sum = *(VECTOR*)Copy();



   if (sum.Length () == Length())

    {

    for (INDEX i = 1; i <= Length(); i++)

      sum.Set (i, Get(i));

    }

   else

    myerr << "operator+ (VECTOR, VECTOR): sum.Length() not ok" << endl;

   return sum;

   }

 */



BASE_VECTOR & VECTOR::operator= (const BASE_VECTOR & v2)

{

  const VECTOR & hv2 = v2.CastToVector();



  SetLength (hv2.Length());



  if (hv2.Length() == Length())

    memcpy (data, hv2.data, sizeof (double) * Length());

  else

    myerr << "Vector::operator=: not allocated" << endl;



  return *this;

}





BASE_VECTOR & VECTOR::operator= (double scal)

{

  if (!Length()) myerr << "VECTOR::operator= (double) : data not allocated"

                       << endl;



  for (INDEX i = 1; i <= Length(); i++)

    Set (i, scal);



  return *this;

}





BASE_VECTOR * VECTOR :: Copy () const

{

  return new VECTOR (*this);

}





void VECTOR :: Swap (BASE_VECTOR & v2)

{

  VECTOR & hv2 = v2.CastToVector();

  swap (length, hv2.length);

  swap (data, hv2.data);

}
