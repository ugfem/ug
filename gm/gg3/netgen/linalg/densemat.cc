// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <fstream.h>
#include <math.h>
#include <string.h>

#include <template.hh>
#include <myadt.hh>
#include <linalg/linalg.hh>


extern ofstream myerr;

DenseMatrix :: DenseMatrix () : BaseMatrix ()
{
  data = NULL;
}

DenseMatrix :: DenseMatrix (INDEX h, INDEX w) : BaseMatrix (h, w)
{
  if (!w) w = h;
  data = new double[h*w];
  if(!data)
  {
    height = width = 0;
    myerr << "Matrix not allocated" << endl;
  }
}

DenseMatrix :: DenseMatrix (const DenseMatrix & m2)
{
  data = NULL;
  SetSize (m2.Height(), m2.Width());

  if (data)
    memcpy (data, m2.data, sizeof(double) * m2.Height() * m2.Width());
  else
    myerr << "DenseMatrix::MATIRX(DenseMatrix&): Matrix not allocated" << endl;
}


DenseMatrix :: DenseMatrix (const BaseMatrix & m2)
{
  data = NULL;
  height = width = 0;

  *this = m2;

  if (!data)
    myerr << "DenseMatrix::DenseMatrix (m2): Matrix not allocated" << endl;
}


DenseMatrix :: ~DenseMatrix ()
{
  if (data) delete [] data;
}


void DenseMatrix :: SetSize (INDEX h, INDEX w)
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


DenseMatrix & DenseMatrix :: operator= (const BaseMatrix & m2)
{
  INDEX i, j;

  SetSize (m2.Height(), m2.Width());

  if (data)
    for (i = 1; i <= Height(); i++)
      for (j = 1; j <= Width(); j++)
        Set (i, j, m2(i, j));
  else
    myerr << "DenseMatrix::Operator=: Matrix not allocated" << endl;

  return *this;
}



DenseMatrix & DenseMatrix :: operator= (const DenseMatrix & m2)
{
  SetSize (m2.Height(), m2.Width());

  if (data)
    memcpy (data, m2.data, sizeof(double) * m2.Height() * m2.Width());
  else
    myerr << "DenseMatrix::Operator=: Matrix not allocated" << endl;

  return *this;
}


DenseMatrix & DenseMatrix :: operator+= (const DenseMatrix & m2)
{
  INDEX i;
  double * p, * q;

  if (Height() != m2.Height() || Width() != m2.Width())
  {
    myerr << "DenseMatrix::Operator+=: Sizes don't fit" << endl;
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
    myerr << "DenseMatrix::Operator+=: Matrix not allocated" << endl;

  return *this;
}

DenseMatrix & DenseMatrix :: operator-= (const DenseMatrix & m2)
{
  INDEX i;
  double * p, * q;

  if (Height() != m2.Height() || Width() != m2.Width())
  {
    myerr << "DenseMatrix::Operator-=: Sizes don't fit" << endl;
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
    myerr << "DenseMatrix::Operator-=: Matrix not allocated" << endl;

  return *this;
}





double & DenseMatrix :: operator() (INDEX i, INDEX j)
{
  if (i >= 1 && j >= 1 && i <= height && j <= width)
    return Elem(i,j);
  else myerr << "\nindex (" << i << "," << j << ") out of range (1.."
             << height << ",1.." << width << ")\n";
  return shit;
}

double DenseMatrix :: operator() (INDEX i, INDEX j) const
{
  if (i >= 1 && j >= 1 && i <= height && j <= width)
    return Get(i,j);
  else myerr << "\nindex (" << i << "," << j << ") out of range (1.."
             << height << ",1.." << width << ")\n";
  return shit;
}


DenseMatrix & DenseMatrix :: operator= (double v)
{
  INDEX i;
  double * p = data;

  if (data)
    for (i = width*height; i > 0; i--, p++)
      *p = v;

  return *this;
}



DenseMatrix & DenseMatrix :: operator*= (double v)
{
  INDEX i;
  double * p = data;

  if (data)
    for (i = width*height; i > 0; i--, p++)
      *p *= v;

  return *this;
}


double DenseMatrix :: Det () const
{
  if (width != height)
  {
    myerr << "DenseMatrix :: Det: width != height" << endl;
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


void CalcInverse (const DenseMatrix & hm1, DenseMatrix & m2)
{
  int i, j, k, n;
  double det, q;
  DenseMatrix m1 = hm1;

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


  if (m1.Width() <= 3)
  {
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
    }
  }

  else
  {
    n = m1.Height();
    m2 = 0;
    for (i = 1; i <= n; i++)
      m2.Elem(i, i) = 1;

    for (i = 1; i <= n; i++)
    {
      q = m1.Get(i, i);
      for (k = 1; k <= n; k++)
      {
        m1.Elem(i, k) /= q;
        m2.Elem(i, k) /= q;
      }

      for (j = i+1; j <= n; j++)
      {
        q = m1.Elem(j, i);
        for (k = 1; k <= n; k++)
        {
          m1.Elem(j, k) -= q * m1.Elem(i, k);
          m2.Elem(j, k) -= q * m2.Elem(i, k);
        }
      }
    }

    for (i = n; i >= 1; i--)
      for (j = 1; j < i; j++)
      {
        q = m1.Elem(j, i);
        for (k = 1; k <= n; k++)
        {
          m1.Elem(j, k) -= q * m1.Elem(i, k);
          m2.Elem(j, k) -= q * m2.Elem(i, k);
        }
      }
  }
}


void CalcAAt (const DenseMatrix & a, DenseMatrix & m2)
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



void CalcAtA (const DenseMatrix & a, DenseMatrix & m2)
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






void CalcABt (const DenseMatrix & a, const DenseMatrix & b, DenseMatrix & m2)
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



void CalcAtB (const DenseMatrix & a, const DenseMatrix & b, DenseMatrix & m2)
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







DenseMatrix operator* (const DenseMatrix & m1, const DenseMatrix & m2)
{
  DenseMatrix temp (m1.Height(), m2.Width());

  if (m1.Width() != m2.Height())
  {
    myerr << "DenseMatrix :: operator*: Matrix Size does not fit" << endl;
  }
  else if (temp.Height() != m1.Height())
  {
    myerr << "DenseMatrix :: operator*: temp not allocated" << endl;
  }
  else
  {
    Mult (m1, m2, temp);
  }
  return temp;
}


void Mult (const DenseMatrix & m1, const DenseMatrix & m2, DenseMatrix & m3)
{
  INDEX n2;
  double sum;
  double *p1, *p1s, *p1sn, *p1snn, *p2, *p2s, *p2sn, *p3;

  if (m1.Width() != m2.Height())
  {
    myerr << "DenseMatrix :: operator*: Matrix Size does not fit" << endl;
    myerr << "m1: " << m1.Height() << " x " << m1.Width() << endl;
    myerr << "m2: " << m2.Height() << " x " << m2.Width() << endl;
  }
  else if (m3.Height() != m1.Height())
  {
    myerr << "DenseMatrix :: operator*: temp not allocated" << endl;
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




DenseMatrix operator+ (const DenseMatrix & m1, const DenseMatrix & m2)
{
  DenseMatrix temp (m1.Height(), m1.Width());
  INDEX i, j;

  if (m1.Width() != m2.Width() || m1.Height() != m2.Height())
  {
    myerr << "BaseMatrix :: operator+: Matrix Size does not fit" << endl;
  }
  else if (temp.Height() != m1.Height())
  {
    myerr << "BaseMatrix :: operator+: temp not allocated" << endl;
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








void DenseMatrix :: Mult (const BaseVector & bv, BaseVector & bprod) const
{
  double sum;
  const double * mp, * sp;
  double * dp;
  const Vector & v = bv.CastToVector();
  Vector & prod = bprod.CastToVector();

  prod.SetLength (Height());

  if (Width() != v.Length())
  {
    myerr << "\nMatrix and Vector don't fit" << endl;
  }
  else if (Height() != prod.Length())
  {
    myerr << "Base_Matrix::operator*(Vector): prod vector not ok" << endl;
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


void DenseMatrix :: MultTrans (const BaseVector & bv, BaseVector & bprod) const
{
  float sum;
  const Vector & v = bv.CastToVector();
  Vector & prod = bprod.CastToVector();

  prod.SetLength (Width());

  if (Height() != v.Length())
  {
    myerr << "\nMatrix and Vector don't fit" << endl;
  }
  else if (Width() != prod.Length())
  {
    myerr << "Base_Matrix::operator*(Vector): prod vector not ok" << endl;
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


void DenseMatrix :: Residuum (const BaseVector & bx, const BaseVector & bb,
                              BaseVector & bres) const
{
  float sum;
  const Vector & x = bx.CastToVector();
  const Vector & b = bb.CastToVector();
  Vector & res = bres.CastToVector();

  res.SetLength (Height());

  if (Width() != x.Length() || Height() != b.Length())
  {
    myerr << "\nMatrix and Vector don't fit" << endl;
  }
  else if (Height() != res.Length())
  {
    myerr << "Base_Matrix::operator*(Vector): prod vector not ok" << endl;
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

double DenseMatrix :: EvaluateBilinearform (const BaseVector & x) const
{
  double sum = 0, hsum;
  const Vector & hx = x.CastToVector();
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


void DenseMatrix :: MultElementMatrix (const ARRAY<INDEX> & pnum,
                                       const BaseVector & x, BaseVector & y)
{
  int i, j;
  const Vector & hx = x.CastToVector();
  Vector & hy = y.CastToVector();

  if (Symmetric())
  {
    for (i = 1; i <= Height(); i++)
    {
      for (j = 1; j < i; j++)
      {
        hy.Elem(pnum.Get(i)) += Get(i, j) * hx.Get(pnum.Get(j));
        hy.Elem(pnum.Get(j)) += Get(i, j) * hx.Get(pnum.Get(i));
      }
      hy.Elem(pnum.Get(j)) += Get(i, i) * hx.Get(pnum.Get(i));
    }
  }
  else
    for (i = 1; i <= Height(); i++)
      for (j = 1; j <= Width(); j++)
        hy.Elem(pnum.Get(i)) += Get(i, j) * hx.Get(pnum.Get(j));

}

void DenseMatrix :: MultTransElementMatrix (const ARRAY<INDEX> & pnum,
                                            const BaseVector & x, BaseVector & y)
{
  int i, j;
  const Vector & hx = x.CastToVector();
  Vector & hy = y.CastToVector();

  if (Symmetric())
    MultElementMatrix (pnum, x, y);
  else
    for (i = 1; i <= Height(); i++)
      for (j = 1; j <= Width(); j++)
        hy.Elem(pnum.Get(i)) += Get(j, i) * hx.Get(pnum.Get(j));
}



void DenseMatrix :: SolveDestroy (const Vector & v, Vector & sol)
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


BaseMatrix * DenseMatrix :: Copy () const
{
  return new DenseMatrix (*this);
}
