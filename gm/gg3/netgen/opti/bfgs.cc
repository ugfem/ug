// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/***************************************************************************/
/*                                                                         */
/* Vorlesung Optimierung I, Gfrerer, WS94/95                               */
/* BFGS-Verfahren zur Lösung freier nichtlinearer Optimierungsprobleme     */
/*                                                                         */
/* Programmautor:  Joachim Schöberl                                        */
/* Matrikelnummer: 9155284                                                 */
/*                                                                         */
/***************************************************************************/

#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <values.h>       // MINDOUBLE, MAXDOUBLE

#include <myadt.hh>  // min, max, sqr

#include <linalg/linalg.hh>
#include <opti/opti.hh>

void MultLDLt (const DenseMatrix & l, const Vector & d, const Vector & g, Vector & p)
{
  int i, j, n;
  double val;

  n = l.Height();
  p = g;
  for (i = 1; i <= n; i++)
  {
    val = 0;
    for (j = i; j <= n; j++)
      val += p.Get(j) * l.Get(j, i);
    p.Set(i, val);
  }
  for (i = 1; i <= n; i++)
    p.Elem(i) *= d.Get(i);

  for (i = n; i >= 1; i--)
  {
    val = 0;
    for (j = 1; j <= i; j++)
      val += p.Get(j) * l.Get(i, j);
    p.Set(i, val);
  }
}

void SolveLDLt (const DenseMatrix & l, const Vector & d, const Vector & g, Vector & p)
{
  int i, j, n;
  double val;

  n = l.Height();
  p = g;

  for (i = 1; i <= n; i++)
  {
    val = 0;
    for (j = 1; j < i; j++)
      val += p.Get(j) * l.Get(i, j);
    p.Elem(i) -= val;
  }
  for (i = 1; i <= n; i++)
    p.Elem(i) /= d.Get(i);

  for (i = n; i >= 1; i--)
  {
    val = 0;
    for (j = i+1; j <= n; j++)
      val += p.Get(j) * l.Get(j, i);
    p.Elem(i) -= val;
  }
}

int LDLtUpdate (DenseMatrix & l, Vector & d, double a, const Vector & u)
{
  // Bemerkung: Es wird a aus R erlaubt
  // Rueckgabewert: 0 .. D bleibt positiv definit
  //               1 .. sonst

  int i, j, n;

  n = l.Height();

  Vector v(n);
  double t, told, xi;

  told = 1;
  v = u;

  for (j = 1; j <= n; j++)
  {
    t = told + a * sqr (v.Elem(j)) / d.Get(j);

    if (t <= 0) return 1;

    xi = a * v.Elem(j) / (d.Get(j) * t);

    d.Elem(j) *= t / told;

    for (i = j + 1; i <= n; i++)
    {
      v.Elem(i) -= v.Elem(j) * l.Elem(i, j);
      l.Elem(i, j) += xi * v.Elem(i);
    }

    told = t;
  }

  return 0;
}


void BFGS (
  Vector & x,         // i: Startwert
                      // o: Loesung, falls IFAIL = 0
  const MinFunction & fun
  )

{

  int i, j, n = x.Length();
  long it;
  char a1crit, a3acrit;


  Vector d(n), g(n), p(n), temp(n), bs(n), xneu(n), y(n), s(n);
  DenseMatrix l(n);

  double /* normg, */ alphahat, hd, fold;
  double a1, a2;
  const double mu1 = 0.1, sigma = 0.1, xi1 = 1, xi2 = 10;
  const double tau = 0.1, tau1 = 0.1, tau2 = 0.6;

  Vector typx(x.Length());      // i: typische Groessenordnung der Komponenten
  double f, typf;               // i: typische Groessenordnung der Loesung
  double fmin = -1e5;           // i: untere Schranke fuer Funktionswert
  double eps = 1e-8;            // i: Abbruchschranke fuer relativen Gradienten
  double tauf = 0.1;            // i: Abbruchschranke fuer die relative Aenderung der
                                //    Funktionswerte
  long itmax = 100;             // i: Maximale Iterationszahl
  int ifail;                    // o:  0 .. Erfolg
                                //    -1 .. Unterschreitung von fmin
                                //     1 .. kein Erfolg bei Liniensuche
                                //     2 .. Überschreitung von itmax

  typx = 1;
  typf = 1;


  l = 0;
  for (i = 1; i <= n; i++)
    l.Elem(i, i) = 1;

  f = fun.FuncGrad (x, g);

  it = 0;
  do
  {
    // Restart

    if (it % (5 * n) == 0)
    {
      for (i = 1; i <= n; i++)
        d.Elem(i) = 1;
      for (i = 2; i <= n; i++)
        for (j = 1; j < i; j++)
          l.Elem(i, j) = 0;
    }

    it++;
    if (it > itmax)
    {
      ifail = 2;
      break;
    }


    // Solve with factorized B

    SolveLDLt (l, d, g, p);


    p *= -1;
    y = g;

    fold = f;

    // line search

    alphahat = 1;
    lines (x, xneu, p, f, g, fun, alphahat, fmin,
           mu1, sigma, xi1, xi2, tau, tau1, tau2, ifail);

    if (ifail != 0) break;

    s = xneu - x;
    y *= -1;
    y += g;

    x = xneu;

    // BFGS Update

    MultLDLt (l, d, s, bs);

    a1 = y * s;
    a2 = s * bs;

    if (a1 > 0 && a2 > 0)
    {
      if (LDLtUpdate (l, d, 1 / a1, y) != 0)
      {
        ifail = 1;
        break;
      }

      if (LDLtUpdate (l, d, -1 / a2, bs) != 0)
      {
        ifail = 1;
        break;
      }
    }

    // Calculate stop conditions

    hd = eps * max (typf, fabs (f));
    a1crit = 1;
    for (i = 1; i <= n; i++)
      if ( fabs (g.Elem(i)) * max (typx.Elem(i), fabs (x.Elem(i))) > hd)
        a1crit = 0;


    a3acrit = (fold - f <= tauf * max (typf, fabs (f)));

    //    testout << "g = " << g << endl;
    //    testout << "a1crit, a3crit = " << int(a1crit) << ", " << int(a3acrit) << endl;

    /*
        // Output for tests

        normg = sqrt (g * g);

        testout << "it =" << setw (5) << it
             << " f =" << setw (12) << setprecision (5) << f
             << " |g| =" << setw (12) << setprecision (5) << normg;

        testout << " x = (" << setw (12) << setprecision (5) << x.Elem(1);
        for (i = 2; i <= n; i++)
          testout << "," << setw (12) << setprecision (5) << x.Elem(i);
        testout << ")" << endl;
     */
  }
  while (!a1crit || !a3acrit);

  //  testout << "it = " << it << " g = " << g << " f = " << f << endl;
}
