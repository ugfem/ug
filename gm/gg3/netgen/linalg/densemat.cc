// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <stdlib.h>
#include <fstream.h>
#include <math.h>
#include <string.h>

#include <template.hh>
#include <myadt.hh>
#include <linalg/linalg.hh>


//extern ofstream myerr;

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
    //    myerr << "Matrix not allocated" << endl;
  }
}

DenseMatrix :: DenseMatrix (INDEX h, INDEX w, const double * d)
  : BaseMatrix (h, w)
{
  int size = h * w;
  int i;

  data = new double[size];
  for (i = 0; i < size; i++)
    data[i] = d[i];
}


DenseMatrix :: DenseMatrix (const DenseMatrix & m2)
{
  data = NULL;
  SetSize (m2.Height(), m2.Width());
  SetSymmetric (m2.Symmetric());

  if (data)
    memcpy (data, m2.data, sizeof(double) * m2.Height() * m2.Width());
  //  else
  //    myerr << "DenseMatrix::MATIRX(DenseMatrix&): Matrix not allocated" << endl;
}


DenseMatrix :: DenseMatrix (const BaseMatrix & m2)
{
  data = NULL;
  height = width = 0;

  *this = m2;

  //  if (!data)
  //    myerr << "DenseMatrix::DenseMatrix (m2): Matrix not allocated" << endl;
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
    //         myerr << "Matrix::SetSize: Matrix not allocated" << endl;
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
  //  else
  //    myerr << "DenseMatrix::Operator=: Matrix not allocated" << endl;

  return *this;
}



DenseMatrix & DenseMatrix :: operator= (const DenseMatrix & m2)
{
  SetSize (m2.Height(), m2.Width());

  if (data)
    memcpy (data, m2.data, sizeof(double) * m2.Height() * m2.Width());
  //  else
  //    myerr << "DenseMatrix::Operator=: Matrix not allocated" << endl;

  return *this;
}


DenseMatrix & DenseMatrix :: operator+= (const DenseMatrix & m2)
{
  INDEX i;
  double * p, * q;

  if (Height() != m2.Height() || Width() != m2.Width())
  {
    //    myerr << "DenseMatrix::Operator+=: Sizes don't fit" << endl;
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
  //  else
  //    myerr << "DenseMatrix::Operator+=: Matrix not allocated" << endl;

  return *this;
}

DenseMatrix & DenseMatrix :: operator-= (const DenseMatrix & m2)
{
  INDEX i;
  double * p, * q;

  if (Height() != m2.Height() || Width() != m2.Width())
  {
    //    myerr << "DenseMatrix::Operator-=: Sizes don't fit" << endl;
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
  //  else
  //    myerr << "DenseMatrix::Operator-=: Matrix not allocated" << endl;

  return *this;
}





double & DenseMatrix :: operator() (INDEX i, INDEX j)
{
  if (i >= 1 && j >= 1 && i <= height && j <= width)
    return Elem(i,j);
  //  else myerr << "DenseMatrix: index (" << i << "," << j << ") out of range (1.."
  //            << height << ",1.." << width << ")\n";
  return shit;
}

double DenseMatrix :: operator() (INDEX i, INDEX j) const
{
  if (i >= 1 && j >= 1 && i <= height && j <= width)
    return Get(i,j);
  //  else myerr << "DenseMatrix: index (" << i << "," << j << ") out of range (1.."
  //            << height << ",1.." << width << ")\n";
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
    //    myerr << "DenseMatrix :: Det: width != height" << endl;
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
    //      myerr << "Matrix :: Det:  general size not implemented (size=" << width << ")" << endl;
    return 0;
  }
  }
}


void CalcInverse (const DenseMatrix & m1, DenseMatrix & m2)
{
  //  int i, j, k, n;
  double det;
  //  DenseMatrix m1 = hm1;

  if (m1.width != m1.height)
  {
    //    myerr << "CalcInverse: matrix not symmetric" << endl;
    return;
  }
  if (m1.width != m2.width || m1.height != m2.height)
  {
    //    myerr << "CalcInverse: dim(m2) != dim(m1)" << endl;
    return;
  }


  if (m1.Width() <= 3)
  {
    det = m1.Det();
    if (det == 0)
    {
      //      myerr << "CalcInverse: Matrix singular" << endl;
      return;
    }

    det = 1e0 / det;
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
    int i, j, k, n;
    n = m1.Height();

    int dots = (n > 200);

    // Cholesky

    double x;
    Vector p(n);

    m2 = m1;
    m2.SetSymmetric();
    if (!m2.Symmetric())
      cerr << "m should be symmetric for Cholesky" << endl;

    for (i = 1; i <= n; i++)
      for (j = 1; j < i; j++)
        m2.Elem(j, i) = m2.Get(i, j);

    for (i = 1; i <= n; i++)
    {
      if (dots && i % 10 == 0)
        cout << "." << flush;

      for (j = i; j <= n; j++)
      {
        x = m2.Get(i, j);

        const double * pik = &m2.Get(i, 1);
        const double * pjk = &m2.Get(j, 1);

        for (k = i-2; k >= 0; --k, ++pik, ++pjk)
          x -= (*pik) * (*pjk);

        // for (k = i-1; k >= 1; --k)
        //   x -= m2.Get(j, k) * m2.Get(i, k);

        if (i == j)
        {
          if (x <= 0)
          {
            cerr << "Matrix indefinite" << endl;
            return;
          }

          p.Elem(i) = 1 / sqrt(x);
        }
        else
        {
          m2.Elem(j, i) = x * p.Get(i);
        }
      }
    }

    for (i = 1; i <= n; i++)
      m2.Elem(i, i) = 1 / p.Get(i);

    // calc L^{-1}, store upper triangle


    //      DenseMatrix hm(n);
    //      hm = m2;

    for (i = 1; i <= n; i++)
    {
      if (dots && i % 10 == 0)
        cout << "+" << flush;

      for (j = i; j <= n; j++)
      {
        x = 0;
        if (j == i) x = 1;

        const double * pjk = &m2.Get(j, i);
        const double * pik = &m2.Get(i, i);
        for (k = i; k < j; k++, ++pjk, ++pik)
          x -= *pik * *pjk;

        //  for (k = i; k < j; k++)
        //  x -= m2.Get(j, k) * m2.Get(i, k);

        m2.Elem(i, j) = x / m2.Get(j, j);
      }
    }

    // calc A^-1 = L^-T * L^-1

    for (i = 1; i <= n; i++)
    {
      if (dots && i % 10 == 0)
        cout << "-" << flush;

      for (j = 1; j <= i; j++)
      {
        x = 0;
        k = i;
        if (j > i) k = j;

        const double * pik = &m2.Get(i, k);
        const double * pjk = &m2.Get(j, k);

        for ( ; k <= n; ++k, ++pik, ++pjk)
          x += *pik * *pjk;
        // for (  ; k <= n; k++)
        //   x += m2.Get(i, k) * m2.Get(j, k);

        m2.Elem(i, j) = x;
      }
    }

    for (i = 1; i <= n; i++)
      for (j = 1; j < i; j++)
        m2.Elem(j, i) = m2.Get(i, j);

    if (dots) cout << endl;

    /*

       Gauss - Jordan - algorithm

       int r, hi;
       double max, hr;


       ARRAY<int> p(n);   // pivot-permutation
       Vector hv(n);


       m2 = m1;

       if (m2.Symmetric())
       for (i = 1; i <= n; i++)
        for (j = 1; j < i; j++)
          m2.Elem(j, i) = m2.Get(i, j);


       // Algorithm of Stoer, Einf. i. d. Num. Math, S 145

       for (j = 1; j <= n; j++)
       p.Set(j, j);

       for (j = 1; j <= n; j++)
       {
        // pivot search

        max = fabs(m2.Get(j, j));
        r = j;

        for (i = j+1; i <= n ;i++)
          if (fabs (m2.Get(i, j)) > max)
            {
              r = i;
              max = fabs (m2.Get(i, j));
            }

        if (max < 1e-20)
          {
            cerr << "Inverse matrix: matrix singular" << endl;
            return;
          }

        r = j;

        // exchange rows
        if (r > j)
          {
            for (k = 1; k <= n; k++)
              {
                hr = m2.Get(j, k);
                m2.Elem(j, k) = m2.Get(r, k);
                m2.Elem(r, k) = hr;
              }
            hi = p.Get(j);
            p.Elem(j) = p.Get(r);
            p.Elem(r) = hi;
          }


        // transformation

        hr = 1 / m2.Get(j, j);
        for (i = 1; i <= n; i++)
          m2.Elem(i, j) *= hr;
        m2.Elem(j, j) = hr;

        for (k = 1; k <= n; k++)
          if (k != j)
            {
              for (i = 1; i <= n; i++)
                if (i != j)
                  m2.Elem(i, k) -= m2.Elem(i, j) * m2.Elem(j, k);
              m2.Elem(j, k) *= -hr;
            }
       }

       // col exchange

       for (i = 1; i <= n; i++)
       {
        for (k = 1; k <= n; k++)
          hv.Elem(p.Get(k)) = m2.Get(i, k);
        for (k = 1; k <= n; k++)
          m2.Elem(i, k) = hv.Get(k);
       }

     */





    /*
       if (m1.Symmetric())
       for (i = 1; i <= n; i++)
        for (j = 1; j < i; j++)
          m1.Elem(j, i) = m1.Get(i, j);

       m2 = 0;

       for (i = 1; i <= n; i++)
       m2.Elem(i, i) = 1;

       for (i = 1; i <= n; i++)
       {
        //	cout << '.' << flush;
       q = m1.Get(i, i);
       for (k = 1; k <= n; k++)
        {
        m1.Elem(i, k) /= q;
        m2.Elem(i, k) /= q;
        }

       for (j = i+1; j <= n; j++)
        {
        q = m1.Elem(j, i);

        double * m1pi = &m1.Elem(i, i);
        double * m1pj = &m1.Elem(j, i);

        for (k = n; k >= i; --k, ++m1pi, ++m1pj)
       *m1pj -= q * (*m1pi);

        double * m2pi = &m2.Elem(i, 1);
        double * m2pj = &m2.Elem(j, 1);

        for (k = i; k > 0; --k, ++m2pi, ++m2pj)
       *m2pj -= q * (*m2pi);

            //        for (k = 1; k <= n; k++)
            //          {
            //          m1.Elem(j, k) -= q * m1.Elem(i, k);
            //          m2.Elem(j, k) -= q * m2.Elem(i, k);
            //          }

        }
       }

       for (i = n; i >= 1; i--)
       {
        //	cout << "+" << flush;
        for (j = 1; j < i; j++)
          {
            q = m1.Elem(j, i);

            double * m2pi = &m2.Elem(i, 1);
            double * m2pj = &m2.Elem(j, 1);

            for (k = n; k > 0; --k, ++m2pi, ++m2pj)
       *m2pj -= q * (*m2pi);


            //	    for (k = 1; k <= n; k++)
            //	      {
            //		m1.Elem(j, k) -= q * m1.Elem(i, k);
            //		m2.Elem(j, k) -= q * m2.Elem(i, k);
            //	      }
          }
       }

       if (m2.Symmetric())
       {
        for (i = 1; i <= n; i++)
          for (j = 1; j < i; j++)
            m2.Elem(i, j) = m2.Elem(j, i);
       }
     */
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
    //    myerr << "CalcAAt: sizes don't fit" << endl;
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




BaseMatrix * DenseMatrix :: InverseMatrix (const BitArray * /* inner */) const
{
  if (Height() != Width())
  {
    //      myerr << "BaseMatrix::InverseMatrix(): Matrix not symmetric" << endl;
    return new DenseMatrix(1);
  }
  else
  {
    if (Symmetric())
    {
      //	  cout << "Invmat not available" << endl;
      BaseMatrix * invmat = NULL;
      return invmat;
    }

    DenseMatrix * invmat = new DenseMatrix (Height());

    CalcInverse (*this, *invmat);
    return invmat;
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
    //    myerr << "CalcAtA: sizes don't fit" << endl;
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
    //    myerr << "CalcABt: sizes don't fit" << endl;
    return;
  }

  double * pm2 = &m2.Elem(1, 1);
  const double * pa1 = &a.Get(1, 1);

  for (i = 1; i <= n1; i++)
  {
    const double * pb = &b.Get(1, 1);
    for (j = 1; j <= n3; j++)
    {
      sum = 0;
      const double * pa = pa1;

      for (k = 1; k <= n2; k++)
      {
        sum += *pa * *pb;
        pa++; pb++;
      }

      *pm2 = sum;
      pm2++;
    }
    pa1 += n2;
  }
}


void CalcAtB (const DenseMatrix & a, const DenseMatrix & b, DenseMatrix & m2)
{
  INDEX n1 = a.Height();
  INDEX n2 = a.Width();
  INDEX n3 = b.Width();
  INDEX i, j, k;

  if (m2.Height() != n2 || m2.Width() != n3 || b.Height() != n1)
  {
    //    myerr << "CalcAtB: sizes don't fit" << endl;
    return;
  }

  for (i = 1; i <= n2 * n3; i++)
    m2.data[i-1] = 0;

  for (i = 1; i <= n1; i++)
    for (j = 1; j <= n2; j++)
    {
      const double va = a.Get(i, j);
      double * pm2 = &m2.Elem(j, 1);
      const double * pb = &b.Get(i, 1);

      for (k = 1; k <= n3; ++k, ++pm2, ++pb)
        *pm2 += va * *pb;
      //	for (k = 1; k <= n3; k++)
      //	  m2.Elem(j, k) += va * b.Get(i, k);
    }
  /*
     for (i = 1; i <= n2; i++)
     for (j = 1; j <= n3; j++)
      {
        sum = 0;
        for (k = 1; k <= n1; k++)
          sum += a.Get(k, i) * b.Get(k, j);
        m2.Elem(i, j) = sum;
      }
   */
}







DenseMatrix operator* (const DenseMatrix & m1, const DenseMatrix & m2)
{
  DenseMatrix temp (m1.Height(), m2.Width());

  if (m1.Width() != m2.Height())
  {
    //    myerr << "DenseMatrix :: operator*: Matrix Size does not fit" << endl;
    ;
  }
  else if (temp.Height() != m1.Height())
  {
    //    myerr << "DenseMatrix :: operator*: temp not allocated" << endl;
    ;
  }
  else
  {
    Mult (m1, m2, temp);
  }
  return temp;
}


void Mult (const DenseMatrix & m1, const DenseMatrix & m2, DenseMatrix & m3)
{
  double sum;
  double *p1, *p1s, *p1sn, *p1snn, *p2, *p2s, *p2sn, *p3;

  if (m1.Width() != m2.Height() || m1.Height() != m3.Height() ||
      m2.Width() != m3.Width() )
  {
    //    myerr << "DenseMatrix :: Mult: Matrix Size does not fit" << endl;
    //    myerr << "m1: " << m1.Height() << " x " << m1.Width() << endl;
    //    myerr << "m2: " << m2.Height() << " x " << m2.Width() << endl;
    //    myerr << "m3: " << m3.Height() << " x " << m3.Width() << endl;
    return;
  }
  else if (m1.Symmetric() || m2.Symmetric() || m3.Symmetric())
  {
    //    myerr << "DenseMatrix :: Mult: not implemented for symmetric matrices" << endl;
    return;
  }
  else
  {
    //      int i, j, k;
    int n1 = m1.Height();
    int n2 = m2.Width();
    int n3 = m1.Width();

    /*
       for (i = n1 * n2-1; i >= 0; --i)
       m3.data[i] = 0;

       const double * pm1 = &m1.Get(1, 1);
       for (i = 1; i <= n1; i++)
       {
        const double * pm2 = &m2.Get(1, 1);
        double * pm3i = &m3.Elem(i, 1);

        for (j = 1; j <= n3; j++)
          {
            const double vm1 = *pm1;
       ++pm1;
            //	      const double vm1 = m1.Get(i, j);
            double * pm3 = pm3i;
            //	      const double * pm2 = &m2.Get(j, 1);

            for (k = 0; k < n2; k++)
              {
       *pm3 += vm1 * *pm2;
       ++pm2;
       ++pm3;
              }

          //	    for (k = 1; k <= n2; k++)
          //	      m3.Elem(i, k) += m1.Get(i, j) * m2.Get(j, k);
          }
       }
     */

    /*
       for (i = 1; i <= n1; i++)
       for (j = 1; j <= n2; j++)
        {
          sum = 0;
          for (k = 1; k <= n3; k++)
            sum += m1.Get(i, k) * m2.Get(k, j);
          m3.Set(i, j, sum);
        }
     */


    /*
       for (i = 1; i <= n1; i++)
       {
        const double pm1i = &m1.Get(i, 1);
        const double pm2j = &m2.Get(1, 1);

        for (j = 1; j <= n2; j++)
          {
            double sum = 0;
            const double * pm1 = pm1i;
            const double * pm2 = pm2j;
            pm2j++;

            for (k = 1; k <= n3; k++)
              {
                sum += *pm1 * *pm2;
       ++pm1;
                pm2 += n2;
              }

            m3.Set (i, j, sum);
          }
       }
     */


    p3 = m3.data;
    p1s = m1.data;
    p2sn = m2.data + n2;
    p1snn = p1s + n1 * n3;

    while (p1s != p1snn)
    {
      p1sn = p1s + n3;
      p2s = m2.data;

      while (p2s != p2sn)
      {
        sum = 0;
        p1 = p1s;
        p2 = p2s;
        p2s++;

        while (p1 != p1sn)
        {
          sum += *p1 * *p2;
          p1++;
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
    //    myerr << "BaseMatrix :: operator+: Matrix Size does not fit" << endl;
    ;
  }
  else if (temp.Height() != m1.Height())
  {
    //    myerr << "BaseMatrix :: operator+: temp not allocated" << endl;
    ;
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




void Transpose (const DenseMatrix & m1, DenseMatrix & m2)
{
  int w = m1.Width();
  int h = m1.Height();
  int i, j;

  m2.SetSize (w, h);

  double * pm2 = &m2.Elem(1, 1);
  for (j = 1; j <= w; j++)
  {
    const double * pm1 = &m1.Get(1, j);
    for (i = 1; i <= h; i++)
    {
      *pm2 = *pm1;
      pm2 ++;
      pm1 += w;
    }
  }
}



void DenseMatrix :: Mult (const BaseVector & bv, BaseVector & bprod) const
{
  double sum, val;
  const double * mp, * sp;
  double * dp;
  const Vector & v = bv.CastToVector();
  Vector & prod = bprod.CastToVector();


  int n = Height();
  int m = Width();

  if (prod.Length() != n)
    prod.SetLength (n);

  if (m != v.Length())
  {
    //    myerr << "\nMatrix and Vector don't fit" << endl;
    ;
  }
  else if (Height() != prod.Length())
  {
    //    myerr << "Base_Matrix::operator*(Vector): prod vector not ok" << endl;
    ;
  }
  else
  {
    if (Symmetric())
    {
      INDEX i, j;


      for (i = 1; i <= n; i++)
      {
        sp = &v.Get(1);
        dp = &prod.Elem(1);
        mp = &Get(i, 1);

        val = v.Get(i);
        sum = Get(i, i) * val;

        for (j = 1; j < i; ++j, ++mp, ++sp, ++dp)
        {
          sum += *mp * *sp;
          *dp += val * *mp;
        }

        prod.Elem(i) = sum;
      }
    }
    else
    {
      mp = data;
      dp = &prod.Elem(1);
      for (INDEX i = 1; i <= n; i++)
      {
        sum = 0;
        sp = &v.Get(1);

        for (INDEX j = 1; j <= m; j++)
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
}


void DenseMatrix :: MultTrans (const BaseVector & bv, BaseVector & bprod) const
{
  const Vector & v = (const Vector &)bv; // .CastToVector();
  Vector & prod = (Vector & )bprod;     // .CastToVector();

  /*
     if (Height() != v.Length())
     {
     //    myerr << "\nMatrix and Vector don't fit" << endl;
     }
     else if (Width() != prod.Length())
     {
     //    myerr << "Base_Matrix::operator*(Vector): prod vector not ok" << endl;
     }
     else
   */
  {
    int i, j;
    int w = Width(), h = Height();
    if (prod.Length() != w)
      prod.SetLength (w);

    const double * pmat = &Get(1, 1);
    const double * pv = &v.Get(1);

    prod = 0;

    for (i = 1; i <= h; i++)
    {
      double val = *pv;
      ++pv;

      double * pprod = &prod.Elem(1);

      for (j = w-1; j >= 0; --j, ++pmat, ++pprod)
      {
        *pprod += val * *pmat;
      }
    }

    /*
       double sum;

       for (i = 1; i <= Width(); i++)
       {
        sum = 0;

        for (INDEX j = 1; j <= Height(); j++)
          sum += Get(j, i) * v.Get(j);

        prod.Set (i, sum);
       }
     */
  }
}


void DenseMatrix :: Residuum (const BaseVector & bx, const BaseVector & bb,
                              BaseVector & bres) const
{
  double sum;
  const Vector & x = bx.CastToVector();
  const Vector & b = bb.CastToVector();
  Vector & res = bres.CastToVector();

  res.SetLength (Height());

  if (Width() != x.Length() || Height() != b.Length())
  {
    //    myerr << "\nMatrix and Vector don't fit" << endl;
    ;
  }
  else if (Height() != res.Length())
  {
    //    myerr << "Base_Matrix::operator*(Vector): prod vector not ok" << endl;
    ;
  }
  else
  {
    int i, j;
    int h = Height();
    int w = Width();
    const double * mp = &Get(1, 1);

    for (i = 1; i <= h; i++)
    {
      sum = b.Get(i);
      const double * xp = &x.Get(1);

      for (j = 1; j <= w; ++j, ++mp, ++xp)
        sum -= *mp * *xp;

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
    //    myerr << "Matrix::EvaluateBilinearForm: sizes don't fit" << endl;
    ;
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
  double q;

  if (Width() != Height())
  {
    //    myerr << "SolveDestroy: Matrix not square";
    return;
  }
  if (Width() != v.Length())
  {
    //    myerr << "SolveDestroy: Matrix and Vector don't fit";
    return;
  }

  sol = v;
  if (Height() != sol.Length())
  {
    //    myerr << "SolveDestroy: Solution Vector not ok";
    return;
  }


  if (0 /* Symmetric() */)
  {

    // Cholesky factorization

    int i, j, k, n;
    n = Height();

    // Cholesky

    double x;
    Vector p(n);

    for (i = 1; i <= n; i++)
      for (j = 1; j < i; j++)
        Elem(j, i) = Get(i, j);

    for (i = 1; i <= n; i++)
    {
      // cout << "." << flush;
      for (j = i; j <= n; j++)
      {
        x = Get(i, j);

        const double * pik = &Get(i, 1);
        const double * pjk = &Get(j, 1);

        for (k = i-2; k >= 0; --k, ++pik, ++pjk)
          x -= (*pik) * (*pjk);

        // for (k = i-1; k >= 1; --k)
        //   x -= Get(j, k) * Get(i, k);

        if (i == j)
        {
          if (x <= 0)
          {
            cerr << "Matrix indefinite" << endl;
            return;
          }

          p.Elem(i) = 1 / sqrt(x);
        }
        else
        {
          Elem(j, i) = x * p.Get(i);
        }
      }
    }

    for (i = 1; i <= n; i++)
      Elem(i, i) = 1 / p.Get(i);

    // A = L L^t
    // L stored in left-lower triangle


    sol = v;

    // Solve L sol = sol

    for (i = 1; i <= n; i++)
    {
      double val = sol.Get(i);

      const double * pij = &Get(i, 1);
      const double * solj = &sol.Get(1);

      for (j = 1; j < i; j++, ++pij, ++solj)
        val -= *pij * *solj;
      //	  for (j = 1; j < i; j++)
      //	    val -= Get(i, j) * sol.Get(j);

      sol.Elem(i) = val / Get(i, i);
    }

    // Solve L^t sol = sol

    for (i = n; i >= 1; i--)
    {
      double val = sol.Get(i) / Get(i, i);
      sol.Elem(i) = val;

      double * solj = &sol.Elem(1);
      const double * pij = &Get(i, 1);

      for (j = 1; j < i; ++j, ++pij, ++solj)
        *solj -= val * *pij;
      //	  for (j = 1; j < i; j++)
      //	    sol.Elem(j) -= Get(i, j) * val;
    }


  }
  else
  {
    //      cout << "gauss" << endl;
    int i, j, k, n = Height();
    for (i = 1; i <= n; i++)
    {
      for (j = i+1; j <= n; j++)
      {
        q = Get(j,i) / Get(i,i);
        if (q)
        {
          const double * pik = &Get(i, i+1);
          double * pjk = &Elem(j, i+1);

          for (k = i+1; k <= n; ++k, ++pik, ++pjk)
            *pjk -= q * *pik;

          //  for (k = i+1; k <= Height(); k++)
          //	Elem(j, k) -= q * Get(i,k);


          sol.Elem(j) -= q * sol.Get(i);
        }
      }
    }

    for (i = n; i >= 1; i--)
    {
      q = sol.Get(i);
      for (j = i+1; j <= n; j++)
        q -= Get(i,j) * sol.Get(j);

      sol.Set(i, q / Get(i,i));
    }
  }
}


BaseMatrix * DenseMatrix :: Copy () const
{
  return new DenseMatrix (*this);
}
