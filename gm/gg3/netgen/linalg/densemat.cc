// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>

#include <fstream.h>

#include <math.h>

#include <string.h>



#include <template.hh>

#include <array.hh>

#include <linalg/linalg.hh>





extern ofstream myerr;



extern ofstream testout;











MATRIX :: MATRIX () : BASE_MATRIX ()

{

  data = NULL;

}



MATRIX :: MATRIX (INDEX h, INDEX w) : BASE_MATRIX (h, w)

{

  if (!w) w = h;

  data = new double[h*w];

  if(!data)

  {

    height = width = 0;

    myerr << "Matrix not allocated" << endl;

  }

}



MATRIX :: MATRIX (const MATRIX & m2)

{

  data = NULL;

  SetSize (m2.Height(), m2.Width());



  if (data)

    memcpy (data, m2.data, sizeof(double) * m2.Height() * m2.Width());

  else

    myerr << "MATRIX::MATIRX(MATRIX&): Matrix not allocated" << endl;

}





MATRIX :: MATRIX (const BASE_MATRIX & m2)

{

  data = NULL;

  height = width = 0;



  *this = m2;



  if (!data)

    myerr << "MATRIX::MATRIX (m2): Matrix not allocated" << endl;

}





MATRIX :: ~MATRIX ()

{

  if (data) delete [] data;

}





void MATRIX :: SetSize (INDEX h, INDEX w)

{

  if (!w) w = h;

  if (height == h && width == w) return;



  height = h;

  width = w;



  if (data) delete[] data;

  data = new double[h*w];



  if(!data)

  {

    height = width = 0;

    myerr << "Matrix::SetSize: Matrix not allocated" << endl;

  }

}





MATRIX & MATRIX :: operator= (const BASE_MATRIX & m2)

{

  INDEX i, j;



  SetSize (m2.Height(), m2.Width());



  if (data)

    for (i = 1; i <= Height(); i++)

      for (j = 1; j <= Width(); j++)

        Set (i, j, m2(i, j));

  else

    myerr << "MATRIX::Operator=: Matrix not allocated" << endl;



  return *this;

}







MATRIX & MATRIX :: operator= (const MATRIX & m2)

{

  SetSize (m2.Height(), m2.Width());



  if (data)

    memcpy (data, m2.data, sizeof(double) * m2.Height() * m2.Width());

  else

    myerr << "MATRIX::Operator=: Matrix not allocated" << endl;



  return *this;

}





MATRIX & MATRIX :: operator+= (const MATRIX & m2)

{

  INDEX i;

  double * p, * q;



  if (Height() != m2.Height() || Width() != m2.Width())

  {

    myerr << "MATRIX::Operator+=: Sizes don't fit" << endl;

    return *this;

  }



  if (data)

  {

    p = data;

    q = m2.data;

    for (i = Width() * Height(); i > 0; i--)

    {

      *p += *q;

      p++;

      q++;

    }

  }

  else

    myerr << "MATRIX::Operator+=: Matrix not allocated" << endl;



  return *this;

}



MATRIX & MATRIX :: operator-= (const MATRIX & m2)

{

  INDEX i;

  double * p, * q;



  if (Height() != m2.Height() || Width() != m2.Width())

  {

    myerr << "MATRIX::Operator-=: Sizes don't fit" << endl;

    return *this;

  }



  if (data)

  {

    p = data;

    q = m2.data;

    for (i = Width() * Height(); i > 0; i--)

    {

      *p -= *q;

      p++;

      q++;

    }

  }

  else

    myerr << "MATRIX::Operator-=: Matrix not allocated" << endl;



  return *this;

}











double & MATRIX :: operator() (INDEX i, INDEX j)

{

  if (i >= 1 && j >= 1 && i <= height && j <= width)

    return Elem(i,j);

  else myerr << "\nindex (" << i << "," << j << ") out of range (1.."

             << height << ",1.." << width << ")\n";

  return shit;

}



double MATRIX :: operator() (INDEX i, INDEX j) const

{

  if (i >= 1 && j >= 1 && i <= height && j <= width)

    return Get(i,j);

  else myerr << "\nindex (" << i << "," << j << ") out of range (1.."

             << height << ",1.." << width << ")\n";

  return shit;

}





MATRIX & MATRIX :: operator= (double v)

{

  INDEX i;

  double * p = data;



  if (data)

    for (i = width*height; i > 0; i--, p++)

      *p = v;



  return *this;

}







MATRIX & MATRIX :: operator*= (double v)

{

  INDEX i;

  double * p = data;



  if (data)

    for (i = width*height; i > 0; i--, p++)

      *p *= v;



  return *this;

}





double MATRIX :: Det () const

{

  if (width != height)

  {

    myerr << "MATRIX :: Det: width != height" << endl;

    return 0;

  }



  switch (width)

  {

  case 1 : return Get(1, 1);

  case 2 : return Get(1) * Get(4) - Get(2) * Get(3);



  case 3 : return Get(1) * Get(5) * Get(9)

           + Get(2) * Get(6) * Get(7)

           + Get(3) * Get(4) * Get(8)

           - Get(1) * Get(6) * Get(8)

           - Get(2) * Get(4) * Get(9)

           - Get(3) * Get(5) * Get(7);

  default :

  {

    myerr << "Matrix :: Det:  general size not implemented (size=" << width << ")" << endl;

    return 0;

  }

  }

}





void CalcInverse (const MATRIX & m1, MATRIX & m2)

{

  double det;



  if (m1.width != m1.height)

  {

    myerr << "CalcInverse: matrix not symmetric" << endl;

    return;

  }

  if (m1.width != m2.width || m1.height != m2.height)

  {

    myerr << "CalcInverse: dim(m2) != dim(m1)" << endl;

    return;

  }



  det = m1.Det();

  if (det == 0)

  {

    myerr << "CalcInverse: Matrix singular" << endl;

    return;

  }



  det = 1/det;

  switch (m1.width)

  {

  case 1 :

  {

    m2.Set(1, 1, det);

    return;

  }

  case 2 :

  {

    m2.Set(1, 1, det * m1.Get(4));

    m2.Set(2, 2, det * m1.Get(1));

    m2.Set(1, 2, - det * m1.Get(2));

    m2.Set(2, 1, - det * m1.Get(3));

    return;

  }

  case 3 :

  {

    m2.Set(1, 1,  det * (m1.Get(5) * m1.Get(9) - m1.Get(6) * m1.Get(8)));

    m2.Set(2, 1, -det * (m1.Get(4) * m1.Get(9) - m1.Get(6) * m1.Get(7)));

    m2.Set(3, 1,  det * (m1.Get(4) * m1.Get(8) - m1.Get(5) * m1.Get(7)));



    m2.Set(1, 2, -det * (m1.Get(2) * m1.Get(9) - m1.Get(3) * m1.Get(8)));

    m2.Set(2, 2,  det * (m1.Get(1) * m1.Get(9) - m1.Get(3) * m1.Get(7)));

    m2.Set(3, 2, -det * (m1.Get(1) * m1.Get(8) - m1.Get(2) * m1.Get(7)));



    m2.Set(1, 3,  det * (m1.Get(2) * m1.Get(6) - m1.Get(3) * m1.Get(5)));

    m2.Set(2, 3, -det * (m1.Get(1) * m1.Get(6) - m1.Get(3) * m1.Get(4)));

    m2.Set(3, 3,  det * (m1.Get(1) * m1.Get(5) - m1.Get(2) * m1.Get(4)));

    return;

  }

  default :

  {

    myerr << "CalcInverse: general Case not implemented (size = " << m1.width << ")" << endl;

    return;

  }

  }

}





void CalcAAt (const MATRIX & a, MATRIX & m2)

{

  INDEX n1 = a.Height();

  INDEX n2 = a.Width();

  INDEX i, j, k;

  double sum;

  const double *p, *q, *p0;



  if (m2.Height() != n1 || m2.Width() != n1)

  {

    myerr << "CalcAAt: sizes don't fit" << endl;

    return;

  }



  for (i = 1; i <= n1; i++)

  {

    sum = 0;

    p = &a.ConstElem(i, 1);

    for (k = 1; k <= n2; k++)

    {

      sum += *p * *p;

      p++;

    }

    m2.Set(i, i, sum);



    p0 = &a.ConstElem(i, 1);

    q = a.data;

    for (j = 1; j < i; j++)

    {

      sum = 0;

      p = p0;



      for (k = 1; k <= n2; k++)

      {

        sum += *p * *q;

        p++;

        q++;

      }

      m2.Set(i, j, sum);

      m2.Set(j, i, sum);

    }

  }

}







void CalcAtA (const MATRIX & a, MATRIX & m2)

{

  INDEX n1 = a.Height();

  INDEX n2 = a.Width();

  INDEX i, j, k;

  double sum;



  if (m2.Height() != n2 || m2.Width() != n2)

  {

    myerr << "CalcAtA: sizes don't fit" << endl;

    return;

  }



  for (i = 1; i <= n2; i++)

    for (j = 1; j <= n2; j++)

    {

      sum = 0;

      for (k = 1; k <= n1; k++)

        sum += a.Get(k, i) * a.Get(k, j);

      m2(i, j) = sum;

    }

}













void CalcABt (const MATRIX & a, const MATRIX & b, MATRIX & m2)

{

  INDEX n1 = a.Height();

  INDEX n2 = a.Width();

  INDEX n3 = b.Height();

  INDEX i, j, k;

  double sum;



  if (m2.Height() != n1 || m2.Width() != n3 || b.Width() != n2)

  {

    myerr << "CalcABt: sizes don't fit" << endl;

    return;

  }



  for (i = 1; i <= n1; i++)

    for (j = 1; j <= n3; j++)

    {

      sum = 0;

      for (k = 1; k <= n2; k++)

        sum += a.Get(i, k) * b.Get(j, k);

      m2(i, j) = sum;

    }

}







void CalcAtB (const MATRIX & a, const MATRIX & b, MATRIX & m2)

{

  INDEX n1 = a.Height();

  INDEX n2 = a.Width();

  INDEX n3 = b.Width();

  INDEX i, j, k;

  double sum;



  if (m2.Height() != n2 || m2.Width() != n3 || b.Height() != n1)

  {

    myerr << "CalcAtB: sizes don't fit" << endl;

    return;

  }



  for (i = 1; i <= n2; i++)

    for (j = 1; j <= n3; j++)

    {

      sum = 0;

      for (k = 1; k <= n1; k++)

        sum += a.Get(k, i) * b.Get(k, j);

      m2(i, j) = sum;

    }

}















MATRIX operator* (const MATRIX & m1, const MATRIX & m2)

{

  MATRIX temp (m1.Height(), m2.Width());



  if (m1.Width() != m2.Height())

  {

    myerr << "MATRIX :: operator*: Matrix Size does not fit" << endl;

  }

  else if (temp.Height() != m1.Height())

  {

    myerr << "MATRIX :: operator*: temp not allocated" << endl;

  }

  else

  {

    Mult (m1, m2, temp);

  }

  return temp;

}





void Mult (const MATRIX & m1, const MATRIX & m2, MATRIX & m3)

{

  INDEX n2;

  double sum;

  double *p1, *p1s, *p1sn, *p1snn, *p2, *p2s, *p2sn, *p3;



  if (m1.Width() != m2.Height())

  {

    myerr << "MATRIX :: operator*: Matrix Size does not fit" << endl;

    myerr << "m1: " << m1.Height() << " x " << m1.Width() << endl;

    myerr << "m2: " << m2.Height() << " x " << m2.Width() << endl;

  }

  else if (m3.Height() != m1.Height())

  {

    myerr << "MATRIX :: operator*: temp not allocated" << endl;

  }

  else

  {

    /*

        for (i = 1; i <= m1.Height(); i++)

          for (j = 1; j <= m2.Width(); j++)

            {

            sum = 0;

            for (k = 1; k <= m1.Width(); k++)

              sum += m1.Get(i, k) * m2.Get(k, j);

            m3.Set(i, j, sum);

            }

     */

    n2 = m2.Width();



    p3 = m3.data;

    p1s = m1.data;

    p2sn = m2.data + m2.Width();

    p1snn = p1s + m1.Width() * m1.Height();



    while (p1s != p1snn)

    {

      p1sn = p1s + m1.Width();

      p2s = m2.data;



      while (p2s != p2sn)

      {

        sum = 0;

        p1 = p1s;

        p2 = p2s;

        p2s++;

        while (p1 != p1sn)

        {

          sum += *p1++ * *p2;

          p2 += n2;

        }

        *p3++ = sum;

      }

      p1s = p1sn;

    }

  }

}









MATRIX operator+ (const MATRIX & m1, const MATRIX & m2)

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

        temp.Set(i, j, m1.Get(i, j) + m2.Get(i, j));

      }

  }

  return temp;

}

















void MATRIX :: Mult (const BASE_VECTOR & bv, BASE_VECTOR & bprod) const

{

  double sum;

  const double * mp, * sp;

  double * dp;

  const VECTOR & v = bv.CastToVector();

  VECTOR & prod = bprod.CastToVector();



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

    mp = data;

    dp = &prod.Elem(1);

    for (INDEX i = 1; i <= Height(); i++)

    {

      sum = 0;

      sp = &v.Get(1);



      for (INDEX j = 1; j <= Width(); j++)

      {

        //        sum += Get(i,j) * v.Get(j);

        sum += *mp * *sp;

        mp++;

        sp++;

      }



      //      prod.Set (i, sum);

      *dp = sum;

      dp++;

    }

  }

}





void MATRIX :: MultTrans (const BASE_VECTOR & bv, BASE_VECTOR & bprod) const

{

  float sum;

  const VECTOR & v = bv.CastToVector();

  VECTOR & prod = bprod.CastToVector();



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

        sum += Get(j, i) * v.Get(j);



      prod.Set (i, sum);

    }

  }

}





void MATRIX :: Residuum (const BASE_VECTOR & bx, const BASE_VECTOR & bb,

                         BASE_VECTOR & bres) const

{

  float sum;

  const VECTOR & x = bx.CastToVector();

  const VECTOR & b = bb.CastToVector();

  VECTOR & res = bres.CastToVector();



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

        sum -= Get(i,j) * x.Get(j);



      res.Set (i, sum);

    }

  }

}



double MATRIX :: EvaluateBilinearform (const BASE_VECTOR & x) const

{

  double sum = 0, hsum;

  const VECTOR & hx = x.CastToVector();

  INDEX i, j;



  if (Width() != hx.Length() || Height() != hx.Length())

  {

    myerr << "Matrix::EvaluateBilinearForm: sizes don't fit" << endl;

  }

  else

  {

    for (i = 1; i <= Height(); i++)

    {

      hsum = 0;

      for (j = 1; j <= Height(); j++)

      {

        hsum += Get(i, j) * hx.Get(j);

      }

      sum += hsum * hx.Get(i);

    }

  }



  //  testout << "sum = " << sum << endl;

  return sum;

}











void MATRIX :: SolveDestroy (const VECTOR & v, VECTOR & sol)

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

      q=Get(j,i) / Get(i,i);

      for (k = i+1; k <= Height(); k++)

      {

        Elem(j, k) -= q * Get(i,k);

      }

      sol.Elem(j) -= q * sol.Get(i);

    }

  }



  for (i = Height(); i >= 1; i--)

  {

    q = sol(i);

    for (j = i+1; j <= Height(); j++)

    {

      q -= Get(i,j) * sol.Get(j);

    }

    sol.Set(i, q / Get(i,i));

  }

}





BASE_MATRIX * MATRIX :: Copy () const

{

  return new MATRIX (*this);

}
