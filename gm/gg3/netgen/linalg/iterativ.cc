// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/************************************************************************/
/*                                                                      */
/* This file is a part of NETGEN                                        */
/*                                                                      */
/* File:   iterativ.cc                                                  */
/* Author: Joachim Schoeberl                                            */
/*                                                                      */
/************************************************************************/



/* Iterative Solvers */

#include <stdlib.h>

#include <fstream.h>

#include <math.h>

#include <limits.h>



#include <template.hh>

#include <array.hh>

#include <bitarray.hh>



#include <linalg/linalg.hh>





ofstream itfile ("iterativ.out");



extern void MyError (char * ch);



extern void IterationDone (int it, double err);







VECTOR BASE_MATRIX :: Jacobi (const VECTOR & b, float eps) const

{

  VECTOR x(b.Length());

  VECTOR r(b.Length());

  INDEX i, n = 0;

  float res = 2 * eps;

  double c = 0;



  for (i = 1; i <= Height(); i++)

    if ((*this)(i, i) > c) c = (*this)(i, i);



  c = 1/c;

  x = 0;



  while (n++ < 10000 && res > eps)

  {

    (*this).Residuum (x, b, r);

    x.Add (c, r);

    res = r.L2Norm() / sqrt (x.Length());

    //    if (! (n % 50)) itfile << "n = " << n << " res = " << res << endl;

  }



  return x;

}



/*

   void BASE_MATRIX :: CG (const BASE_VECTOR & b, BASE_VECTOR & x, float eps) const

   {

   const INDEX len = b.Length();

   VECTOR r(len), p(len), h(len);

   INDEX n = 0;

   double res, lam, hd;



   VECTOR & hx = x.CastToVector ();

   hx.SetLength (Width());



   hx = 0;

   p = r = b;

   eps *= b.L2Norm();

   res = 2 * eps;



   //  itfile << "this = " << (*this) << endl;

   //  itfile << "x = " << x << endl;

   //  itfile << "b = " << b << endl;





   while (n++ < 50000 && res > eps)

    {

    (*this).Mult (p, h);



    hd = h * p;

    if (!hd) break;

    lam = (r * p) / hd;



    hx.Add (lam, p);

    r.Add (-lam, h);

    p.Set (1, r, -(r*h) / hd, p);

   //    p = r - (r*h) / hd * p;

    res = r.L2Norm();

   //    res = lam;



   //    itfile << "n = " << n << " res = " << res << " hd = " << hd << endl;



   //    if (n % 1000 == 0)

   //      cout << "n = " << n << " res = " << res << endl;

   //    itfile << "n = " << n << " res = " << res << "hd = " << hd

   //           << "lam = " << lam << endl;

    }

   itfile << n << " its" << endl;

   itfile << "err = " << res << endl;

   }

 */





void BASE_MATRIX :: CG (const BASE_VECTOR & b, BASE_VECTOR & x, float eps) const

{

  BASE_VECTOR & d = *b.Copy();

  BASE_VECTOR & w = *b.Copy();

  BASE_VECTOR & s = *b.Copy();

  BASE_VECTOR & ks = *b.Copy();



  INDEX n = 0;

  double al, be, wd, wdn, kss;



  x = 0;

  d = b;

  w = d;

  //  pre.Mult (d, w);

  s = w;



  eps = eps*eps;

  wdn = w * d;

  while (n++ < 10000 && wdn > eps)

  {

    (*this).Mult (s, ks);



    wd = wdn;

    kss = ks * s;



    if (!kss) break;



    al = wd / kss;



    x.Add (al, s);

    d.Add (-al, ks);



    //    pre.Mult (d, w);

    w = d;



    wdn = w * d;

    be = wdn / wd;



    s *= be;

    s.Add (1, w);



    //    itfile << n << " " << wdn << endl;

  }

  //  itfile << "CG: " << n << " Iterations" << endl;



  d.~BASE_VECTOR();

  w.~BASE_VECTOR();

  s.~BASE_VECTOR();

  ks.~BASE_VECTOR();

}



















VECTOR BASE_MATRIX :: CG (const VECTOR & b, float eps) const

{

  const INDEX len = b.Length();

  VECTOR r(len), p(len), x(len), h(len);

  INDEX n = 0;

  double res, lam, hd;



  x = 0;

  p = r = b;

  eps *= b.L2Norm();

  res = 2 * eps;



  while (n++ < 1000 && res > eps)

  {

    (*this).Mult (p, h);



    hd = h * p;

    if (!hd) break;

    lam = (r * p) / hd;



    x.Add (lam, p);

    r.Add (-lam, h);

    p = r - (r*h) / hd * p;



    res = r.L2Norm();

    itfile << "n = " << n << " res = " << res << endl;

  }



  return x;

}









void BASE_MATRIX :: PCG (const BASE_MATRIX & pre, const VECTOR & f, VECTOR & u,

                         float eps) const

{

  const INDEX len = f.Length();

  VECTOR d(len), w(len), s(len), ks(len);

  INDEX n = 0;

  double al, be, wd, wdn, kss;



  u = 0;

  d = f;

  pre.Mult (d, w);

  s = w;



  wdn = w * d;

  eps = eps*eps * wdn + 1e-20;



  while (n++ < 10000 && wdn > eps)

  {

    (*this).Mult (s, ks);



    wd = wdn;

    kss = ks * s;



    //    itfile << "wdn = " << wdn << "  kss = " << kss << endl;



    if (!kss) break;



    al = wd / kss;



    u.Add (al, s);

    d.Add (-al, ks);



    pre.Mult (d, w);



    wdn = w * d;

    be = wdn / wd;



    s *= be;

    s.Add (1, w);



    //    cout << n << " " << wdn << endl;

  }

  itfile << "PCG: " << n << " Iterations" << endl;

}









void SPARSE_MATRIX :: GSStep (const VECTOR & b, VECTOR & x, double dump) const

{

  double sum, diag, val;

  INDEX i, j, n;

  colstruct * col;

  const linestruct * lin;



  n = Height();

  if (n != b.Length() || n != x.Length())

  {

    MyError ("GSStep: Sizes don't fit");

    return;

  }



  if (! Symmetric ())

  {

    for (i = 1; i <= n; i++)

    {

      lin = &lins.Get(i);

      sum = b.Get(i);

      col = lin->col;

      diag = 1;



      for (j = lin->size; j > 0; j--, col++)

      {

        sum -= col->data * x.Get(col->colnr);

        if (col->colnr == i) diag = col->data;

      }



      x.Elem(i) += dump / diag * sum;

    }

  }

  else

  {

    for (i = 1; i <= n; i++)

    {

      lin = &lins.Get(i);

      col = lin->col;

      val = x.Get(i);



      for (j = lin->size; j > 1; j--, col++)

        x.Elem(col->colnr) -= col->data * val;



      x.Set(i, b.Get(i));

    }



    for (i = 1; i <= n; i++)

    {

      lin = &lins.Get(i);

      col = lin->col;

      sum = x.Get(i);



      for (j = lin->size; j > 1; j--, col++)

        sum -= col->data * x.Get(col->colnr);



      x.Set(i, sum / col->data);

    }

  }

}



void SPARSE_MATRIX :: GSStepBack (const VECTOR & b, VECTOR & x, double dump) const

{

  double sum, diag, val;

  INDEX i, j, n;

  colstruct * col;

  const linestruct * lin;



  n = Height();



  if (n != b.Length() || n != x.Length())

  {

    MyError ("GSStepBack: Sizes don't fit");

    return;

  }



  if (! Symmetric ())

  {

    for (i = n; i >= 1; i--)

    {

      lin = &lins.Get(i);

      sum = b.Get(i);

      col = lin->col;

      diag = 1;



      for (j = lin->size; j > 0; j--, col++)

      {

        sum -= col->data * x.Get(col->colnr);

        if (col->colnr == i) diag = col->data;

      }



      x.Elem(i) += dump / diag * sum;

    }

  }

  else

  {

    for (i = n; i >= 1; i--)

    {

      lin = &lins.Get(i);

      col = lin->col;

      sum = b.Get(i);



      for (j = lin->size; j > 1; j--, col++)

        sum -= col->data * x.Get(col->colnr);



      x.Set(i, sum);

    }

    for (i = n; i >= 1; i--)

    {

      lin = &lins.Get(i);

      col = lin->col;

      //      x.Elem(i) /= GetData(i, ElementsInLine(i));

      x.Elem(i) /= col[lin->size - 1].data;

      val = x.Get(i);



      for (j = lin->size; j > 1; j--, col++)

        x.Elem(col->colnr) -= col->data * val;

    }

  }

  /*

      {

      for (i = n; i >= 1; i--)

        {

        sum = b.Get(i);

        for (j = 1; j < ElementsInLine(i); j++)

          {

          colnr = GetIndex (i, j);

          sum -= GetData (i, j) * x.Get(colnr);

          }

        x.Set(i, sum);

        }

      for (i = n; i >= 1; i--)

        {

        x.Elem(i) /= GetData(i, ElementsInLine(i));



        for (j = 1; j < ElementsInLine(i); j++)

          {

          colnr = GetIndex (i, j);

          x.Elem(colnr) -= GetData (i, j) * x.Get(i);

          }

        }

      }

   */

}



















void SPARSE_MATRIX :: GSStepInner (const VECTOR & b, VECTOR & x, double dump,

                                   const BITARRAY & inner) const

{

  double sum, diag, val;

  INDEX i, j, n;

  colstruct * col;

  const linestruct * lin;



  n = Height();

  if (n != b.Length() || n != x.Length())

  {

    MyError ("GSStep: Sizes don't fit");

    return;

  }



  if (! Symmetric ())

  {

    for (i = 1; i <= n; i++)

    {

      lin = &lins.Get(i);

      sum = b.Get(i);

      col = lin->col;

      diag = 1;



      for (j = lin->size; j > 0; j--, col++)

      {

        sum -= col->data * x.Get(col->colnr);

        if (col->colnr == i) diag = col->data;

      }



      x.Elem(i) += dump / diag * sum;

    }

  }

  else

  {

    for (i = 1; i <= n; i++)

      if (inner.Test(i))

      {

        lin = &lins.Get(i);

        col = lin->col;

        val = x.Get(i);



        for (j = lin->size; j > 1; j--, col++)

          if (inner.Test(col->colnr))

            x.Elem(col->colnr) -= col->data * val;



        x.Set(i, b.Get(i));

      }





    for (i = 1; i <= n; i++)

      if (inner.Test(i))

      {

        lin = &lins.Get(i);

        col = lin->col;

        sum = x.Get(i);



        for (j = lin->size; j > 1; j--, col++)

          if (inner.Test(col->colnr))

            sum -= col->data * x.Get(col->colnr);



        x.Set(i, sum / col->data);

      }

  }

}



void SPARSE_MATRIX :: GSStepBackInner (const VECTOR & b, VECTOR & x,

                                       double dump, const BITARRAY & inner) const

{

  double sum, diag, val;

  INDEX i, j, n;

  colstruct * col;

  const linestruct * lin;



  n = Height();



  if (n != b.Length() || n != x.Length())

  {

    MyError ("GSStepBack: Sizes don't fit");

    return;

  }



  if (! Symmetric ())

  {

    for (i = n; i >= 1; i--)

    {

      lin = &lins.Get(i);

      sum = b.Get(i);

      col = lin->col;

      diag = 1;



      for (j = lin->size; j > 0; j--, col++)

      {

        sum -= col->data * x.Get(col->colnr);

        if (col->colnr == i) diag = col->data;

      }



      x.Elem(i) += dump / diag * sum;

    }

  }

  else

  {

    for (i = n; i >= 1; i--)

      if (inner.Test(i))

      {

        lin = &lins.Get(i);

        col = lin->col;

        sum = b.Get(i);



        for (j = lin->size; j > 1; j--, col++)

          if (inner.Test(col->colnr))

            sum -= col->data * x.Get(col->colnr);



        x.Set(i, sum);

      }



    for (i = n; i >= 1; i--)

      if (inner.Test(i))

      {

        lin = &lins.Get(i);

        col = lin->col;

        x.Elem(i) /= col[lin->size - 1].data;

        val = x.Get(i);



        for (j = lin->size; j > 1; j--, col++)

          if (inner.Test(col->colnr))

            x.Elem(col->colnr) -= col->data * val;

      }



  }

  /*

      {

      for (i = n; i >= 1; i--)

        {

        sum = b.Get(i);

        for (j = 1; j < ElementsInLine(i); j++)

          {

          colnr = GetIndex (i, j);

          sum -= GetData (i, j) * x.Get(colnr);

          }

        x.Set(i, sum);

        }

      for (i = n; i >= 1; i--)

        {

        x.Elem(i) /= GetData(i, ElementsInLine(i));



        for (j = 1; j < ElementsInLine(i); j++)

          {

          colnr = GetIndex (i, j);

          x.Elem(colnr) -= GetData (i, j) * x.Get(i);

          }

        }

      }

   */

}

























void SPARSE_MATRIX :: GSStepToInner (const VECTOR & b, VECTOR & x, double dump,

                                     const BITARRAY & inner) const

{

  double sum, diag, val;

  INDEX i, j, n;

  colstruct * col;

  const linestruct * lin;



  n = Height();

  if (n != b.Length() || n != x.Length())

  {

    MyError ("GSStep: Sizes don't fit");

    return;

  }



  if (! Symmetric ())

  {

    MyError ("GSStepToInner: non symmetric case not implemented");

  }

  else

  {

    for (i = 1; i <= n; i++)

    //      if (inner.Test(i))

    {

      lin = &lins.Get(i);

      col = lin->col;

      val = x.Get(i);



      for (j = lin->size; j > 1; j--, col++)

        if (inner.Test(col->colnr))

          x.Elem(col->colnr) -= col->data * val;



      if (inner.Test(i))

        x.Set(i, b.Get(i));

    }



    for (i = 1; i <= n; i++)

      if (inner.Test(i))

      {

        lin = &lins.Get(i);

        col = lin->col;

        sum = x.Get(i);



        for (j = lin->size; j > 1; j--, col++)

          //          if (inner.Test(col->colnr))

          sum -= col->data * x.Get(col->colnr);



        x.Set(i, sum / col->data);

      }

  }

}



void SPARSE_MATRIX :: GSStepBackToInner (const VECTOR & b, VECTOR & x,

                                         double dump, const BITARRAY & inner) const

{

  double sum, diag, val;

  INDEX i, j, n;

  colstruct * col;

  const linestruct * lin;

  //  VECTOR hx(x.Length());



  n = Height();



  if (n != b.Length() || n != x.Length())

  {

    MyError ("GSStepBack: Sizes don't fit");

    return;

  }



  if (! Symmetric ())

  {

    MyError ("GSStepToInner: non symmetric case not implemented");

  }

  else

  {

    for (i = n; i >= 1; i--)

      if (inner.Test(i))

      {

        lin = &lins.Get(i);

        col = lin->col;

        val = x.Get(i) / col[lin->size - 1].data;

        x.Elem(i) = val;



        for (j = lin->size; j > 1; j--, col++)

          //          if (inner.Test(col->colnr))

          x.Elem(col->colnr) -= col->data * val;

      }





    for (i = n; i >= 1; i--)

    //      if (inner.Test(i))

    {

      lin = &lins.Get(i);

      col = lin->col;



      sum = b.Get(i);

      if (!inner.Test(i)) sum += x.Get(i);



      for (j = lin->size; j > 1; j--, col++)

        if (inner.Test(col->colnr))

          sum -= col->data * x.Get(col->colnr);



      x.Set(i, sum);

    }

  }

}































VECTOR SPARSE_MATRIX :: ILU (const VECTOR & b, float eps) const

{

  VECTOR x(b.Length()), r(b.Length()), temp(b.Length()), temp2(b.Length());

  INDEX n = 0;

  float res = 2 * eps;



  LOWER_SPARSE_MATRIX l;

  UPPER_SPARSE_MATRIX u;



  ILU_Decomposition (l, u);



  //  x = 0;



  l.Solve (b, temp);

  u.Solve (temp, x);

  Residuum (x, b, r);

  res = r.L2Norm() / sqrt (x.Length());



  while (n++ < 2000 && res > eps)

  {

    l.Solve (r, temp);

    u.Solve (temp, temp2);



    x += temp2;



    Residuum (x, b, r);

    res = r.L2Norm() / sqrt (x.Length());



    //    if (n % 50 == 0) itfile << "n = " << n << " res = " << res << endl;

  }

  return x;

}







void SPARSE_MATRIX :: ILU_PCG (const VECTOR & b, VECTOR & x, const LOWER_SPARSE_MATRIX & l,

                               const UPPER_SPARSE_MATRIX & u, float eps) const

{

  const INDEX len = b.Length();

  VECTOR r(len), p(len), h(len), h2(len), h3(len);

  INDEX n = 0;

  double res = 2 * eps, lam, shp, shh2;



  x = 0;

  r = b;



  l.Solve (r, h2);

  u.Solve (h2, h3);

  //  p = u.Solve (l.Solve (r));

  p = h3;



  while (n++ < 1000 && res > eps)

  {

    (*this).Mult (p, h);



    shp = h * p;

    if (!shp) break;

    lam = (r * p) / shp;



    x.Add (lam, p);

    r.Add (-lam, h);



    l.Solve (r, h3);

    u.Solve (h3, h2);

    //    h2 = u.Solve (l.Solve (r));

    shh2 = h * h2;



    h2.Add (-shh2 / shp, p);

    p = h2;





    res = r.L2Norm() / sqrt (len);

    itfile << "n = " << n << " res = " << res << endl;

  }

  itfile << "ILU_PCG: " << n << " Iterations" << endl;

}









VECTOR SPARSE_MATRIX :: ILU_PCG (const VECTOR & b, float eps) const

{

  VECTOR x(b.Length());



  LOWER_SPARSE_MATRIX l;

  UPPER_SPARSE_MATRIX u;



  ILU_Decomposition (l, u);





  ILU_PCG (b, x, l, u, eps);



  return x;

}





void BASE_MATRIX :: PrecondRichardson (const BASE_MATRIX & pre,

                                       const BASE_VECTOR & b, BASE_VECTOR & x, double tau, double eps) const

{

  double err, err0;

  int it;

  BASE_VECTOR * d = b.Copy();

  BASE_VECTOR * w = b.Copy();



  it = 0;



  do

  {

    it++;



    Residuum (x, b, *d);

    pre.Mult (*d, *w);

    x.Add (tau, *w);



    err = sqrt (*d * *w);

    if (it == 1) err0 = err;

    if (err0 < 1e-30) break;

    IterationDone (it, err/err0);

  }

  while (err > eps * err0);



  delete d;

  delete w;

}
