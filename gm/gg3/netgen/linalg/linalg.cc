// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>

#include <fstream.h>

#include <string.h>

#include <math.h>



#include <template.hh>

#include <array.hh>

#include <linalg/linalg.hh>





ofstream myerr ("error.out");

// ofstream myerr ("NUL");









double BASE_MATRIX :: shit = 0;





BASE_MATRIX :: BASE_MATRIX ()

{

  height = width = 0;

  symmetric = 0;

}



BASE_MATRIX :: BASE_MATRIX (INDEX h, INDEX w)

{

  if (!w) w = h;

  height = h;

  width = w;

  symmetric = 0;

}



void BASE_MATRIX :: SetSize (INDEX h, INDEX w)

{

  if (!w) w = h;

  height = h;

  width = w;

}



void BASE_MATRIX :: SetSymmetric (int sym)

{

  symmetric = sym;

}



double & BASE_MATRIX :: operator() (INDEX, INDEX)

{

  myerr << "BASE_MATRIX: operator() called" << endl;

  return shit;

}



double BASE_MATRIX :: operator() (INDEX, INDEX) const

{

  myerr << "BASE_MATRIX: operator() called" << endl;

  return 0;

}







ostream & operator<<(ostream & s, const BASE_MATRIX & m)

{

  return m.Print (s);

}



ostream & BASE_MATRIX :: Print (ostream & s) const

{

  for (INDEX i = 1; i <= Height(); i++)

  {

    for (INDEX j = 1; j < Width(); j++)

      s << (*this)(i, j) << "  ";

    s << (*this)(i, Width()) << endl;

  }



  return s;

}







TEMP_VECTOR BASE_MATRIX :: operator* (const BASE_VECTOR & v) const

{

  VECTOR * prod = new VECTOR(Height());



  if (Width() != v.Length())

  {

    myerr << "\nMatrix and Vector don't fit" << endl;

  }

  else if (Height() != prod->Length())

  {

    myerr << "Base_Matrix::operator*(VECTOR): prod vector not ok" << endl;

  }

  else

  {

    Mult (v, *prod);

  }



  return *prod;

}







MATRIX operator* (const BASE_MATRIX & m1, const BASE_MATRIX & m2)

{

  MATRIX temp (m1.Height(), m2.Width());

  double sum;



  if (m1.Width() != m2.Height())

  {

    myerr << "BASE_MATRIX :: operator*: Matrix Size does not fit" << endl;

  }

  else if (temp.Height() != m1.Height())

  {

    myerr << "BASE_MATRIX :: operator*: temp not allocated" << endl;

  }

  else

  {

    for (INDEX i = 1; i <= m1.Height(); i++)

      for (INDEX j = 1; j <= m2.Width(); j++)

      {

        sum = 0;

        for (INDEX k = 1; k <= m1.Width(); k++)

          sum += m1(i, k) * m2(k, j);

        temp(i, j) = sum;

      }

  }

  return temp;

}





MATRIX operator+ (const BASE_MATRIX & m1, const BASE_MATRIX & m2)

{

  MATRIX temp (m1.Height(), m1.Width());

  INDEX i, j;



  if (m1.Width() != m2.Width() || m1.Height() != m2.Height())

  {

    myerr << "BASE_MATRIX :: operator+: Matrix Size does not fit" << endl;

  }

  else if (temp.Height() != m1.Height())

  {

    myerr << "BASE_MATRIX :: operator+: temp not allocated" << endl;

  }

  else

  {

    for (i = 1; i <= m1.Height(); i++)

      for (j = 1; j <= m1.Width(); j++)

      {

        temp(i, j) = m1(i, j) + m2(i, j);

      }

  }

  return temp;

}





void BASE_MATRIX :: Mult (const BASE_VECTOR & /* v */,

                          BASE_VECTOR & /* prod */) const

{

  myerr << "BASE_MATRIX :: Mult called" << endl;

}



void BASE_MATRIX :: MultTrans (const BASE_VECTOR & /* v */,

                               BASE_VECTOR & /* prod */) const

{

  myerr << "BASE_MATRIX :: MultTrans called" << endl;

}



void BASE_MATRIX :: Residuum (const BASE_VECTOR & /* x */,

                              const BASE_VECTOR & /* b */, BASE_VECTOR & /* res */) const

{

  myerr << "BASE_MATRIX :: Residuum called" << endl;

}



BASE_MATRIX * BASE_MATRIX :: Copy () const

{

  myerr << "BASE_MATRIX :: Copy called" << endl;

  return NULL;

}





BASE_VECTOR * BASE_MATRIX :: CreateVector () const

{

  return new VECTOR (Height());

}







/*

   void BASE_MATRIX :: Mult (const VECTOR & v, VECTOR & prod) const

   {

   float sum;



   prod.SetLength (Height());



   if (Width() != v.Length())

    {

    myerr << "\nMatrix and Vector don't fit" << endl;

    }

   else if (Height() != prod.Length())

    {

    myerr << "Base_Matrix::operator*(VECTOR): prod vector not ok" << endl;

    }

   else

    {

    for (INDEX i = 1; i <= Height(); i++)

      {

      sum = 0;



      for (INDEX j = 1; j <= Width(); j++)

        sum += (*this)(i,j) * v.Get(j);



      prod.Set (i, sum);

      }

    }

   }





   void BASE_MATRIX :: MultTrans (const VECTOR & v, VECTOR & prod) const

   {

   float sum;



   prod.SetLength (Width());



   if (Height() != v.Length())

    {

    myerr << "\nMatrix and Vector don't fit" << endl;

    }

   else if (Width() != prod.Length())

    {

    myerr << "Base_Matrix::operator*(VECTOR): prod vector not ok" << endl;

    }

   else

    {

    for (INDEX i = 1; i <= Width(); i++)

      {

      sum = 0;



      for (INDEX j = 1; j <= Height(); j++)

        sum += (*this)(j, i) * v.Get(j);



      prod.Set (i, sum);

      }

    }

   }





   void BASE_MATRIX :: Residuum (const VECTOR & x, const VECTOR & b, VECTOR & res) const

   {

   float sum;



   res.SetLength (Height());



   if (Width() != x.Length() || Height() != b.Length())

    {

    myerr << "\nMatrix and Vector don't fit" << endl;

    }

   else if (Height() != res.Length())

    {

    myerr << "Base_Matrix::operator*(VECTOR): prod vector not ok" << endl;

    }

   else

    {

    for (INDEX i = 1; i <= Height(); i++)

      {

      sum = b.Get(i);



      for (INDEX j = 1; j <= Width(); j++)

        sum -= (*this)(i,j) * x.Get(j);



      res.Set (i, sum);

      }

    }

   }

 */















void BASE_MATRIX :: SolveDestroy (const VECTOR & v, VECTOR & sol)

{

  INDEX i, j, k;

  double q;



  if (Width() != Height())

  {

    myerr << "SolveDestroy: Matrix not square";

    return;

  }

  if (Width() != v.Length())

  {

    myerr << "SolveDestroy: Matrix and Vector don't fit";

    return;

  }



  sol = v;

  if (Height() != sol.Length())

  {

    myerr << "SolveDestroy: Solution Vector not ok";

    return;

  }



  for (i = 1; i <= Height(); i++)

  {

    for (j = i+1; j <= Height(); j++)

    {

      q=(*this)(j,i) / (*this)(i,i);

      for (k = i+1; k <= Height(); k++)

      {

        (*this)(j, k) -= q * (*this)(i,k);

      }

      sol.Elem(j) -= q * sol.Get(i);

    }

  }



  for (i = Height(); i >= 1; i--)

  {

    q = sol(i);

    for (j = i+1; j <= Height(); j++)

    {

      q -= (*this)(i,j) * sol.Get(j);

    }

    sol.Set(i, q / (*this)(i,i));

  }

}



void BASE_MATRIX :: Solve (const VECTOR & v, VECTOR & sol) const

{

  BASE_MATRIX * temp = Copy();



  if (temp->Height() != Height())

  {

    myerr << "Solve: Matrix temp not allocated" << endl;

    return;

  }



  temp->SolveDestroy (v, sol);



  delete temp;

}







VECTOR BASE_MATRIX :: Solve (const VECTOR & v) const

{

  VECTOR sol (v.Length());



  if (Width() != Height())

  {

    myerr << "Solve: Matrix not square";

    return v;

  }

  if (Width() != v.Length())

  {

    myerr << "Solve: Matrix and Vector don't fit";

    return v;

  }

  if (Width() != sol.Length())

  {

    myerr << "Solve: Vector sol not allocated" << endl;

  }



  Solve (v, sol);



  return sol;

}





void BASE_MATRIX :: SolvePI (const VECTOR & v, VECTOR & sol) const

{

  /*
     VECTOR sol2(v.Length());



     Solve (v, sol2);

     Solve (*(VECTOR*)&(v - (*this) * sol2), sol);

     sol += sol2;
   */
}















void BASE_MATRIX :: LU_Decomposition (MATRIX & l, MATRIX & u) const

{

  INDEX i, j ,k;

  double sum;

  l.SetSize (Width());

  u.SetSize (Width());



  for (i = 1; i <= Width(); i++)

    for (j = 1; j <= Width(); j++)

      l(i, j) = u(i, j) = 0;



  for (i = 1; i <= Width(); i++)

  {

    for (k = 1; k < i; k++)

    {

      sum = (*this)(i, k);

      for (j = 1; j < k; j++)

        sum -= l(i, j) * u(j, k);

      l(i, k) = sum / u(k, k);

    }

    l(i, i) = 1;



    for (k = i; k <= Width(); k++)

    {

      sum = (*this)(i, k);

      for (j = 1; j < i; j++)

        sum -= l(i, j) * u(j, k);

      u(i, k) = sum;

    }

  }

}







void Transpose (const BASE_MATRIX & m1, MATRIX & m2)

{

  m2.SetSize (m1.Width(), m1.Height());

  INDEX i, j;



  for (i = 1; i <= m1.Height(); i++)

    for (j = 1; j <= m1.Width(); j++)

      m2(j, i) = m1(i, j);

}











void BASE_MATRIX :: Gerschgorin (double & min, double & max) const

{

  min = 0;

  max = 0;

}
