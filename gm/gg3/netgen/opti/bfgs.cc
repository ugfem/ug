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

#include <template.hh>  // min, max, sqr
#include <array.hh>
#include <linalg/linalg.hh>

extern ofstream testout;

extern void lines (
  VECTOR & x,         // i: Ausgangspunkt der Liniensuche
  VECTOR & xneu,      // o: Loesung der Liniensuche bei Erfolg
  VECTOR & p,         // i: Suchrichtung
  double & f,         // i: Funktionswert an der Stelle x
                      // o: Funktionswert an der Stelle xneu, falls ifail = 0
  VECTOR & g,         // i: Gradient an der Stelle x
                      // o: Gradient an der Stelle xneu, falls ifail = 0

  double (*fun)(const VECTOR & x, VECTOR & g),

  double & alphahat,  // i: Startwert für alpha_hat
                      // o: Loesung falls ifail = 0
  double fmin,        // i: untere Schranke für f
  double mu1,         // i: Parameter mu_1 aus Alg.2.1
  double sigma,       // i: Parameter sigma aus Alg.2.1
  double xi1,         // i: Parameter xi_1 aus Alg.2.1
  double xi2,         // i: Parameter xi_1 aus Alg.2.1
  double tau,         // i: Parameter tau aus Alg.2.1
  double tau1,        // i: Parameter tau_1 aus Alg.2.1
  double tau2,        // i: Parameter tau_2 aus Alg.2.1
  int & ifail);        // o:  0 bei erfolgreicher Liniensuche
                       //    -1 bei Abbruch wegen Unterschreiten von fmin
                       //     1 bei Abbruch, aus sonstigen Gründen


void MultLDLt (const MATRIX & l, const VECTOR & d, const VECTOR & g, VECTOR & p)
{
  int i, j, n;
  double val;

  n = l.Height();
  p = g;
  for (i = 1; i <= n; i++)
  {
    val = 0;
    for (j = i; j <= n; j++)
      val += p(j) * l(j, i);
    p(i) = val;
  }
  for (i = 1; i <= n; i++)
    p(i) *= d(i);

  for (i = n; i >= 1; i--)
  {
    val = 0;
    for (j = 1; j <= i; j++)
      val += p(j) * l(i, j);
    p(i) = val;
  }
}

void SolveLDLt (const MATRIX & l, const VECTOR & d, const VECTOR & g, VECTOR & p)
{
  int i, j, n;
  double val;

  n = l.Height();
  p = g;

  for (i = 1; i <= n; i++)
  {
    val = 0;
    for (j = 1; j < i; j++)
      val += p(j) * l(i, j);
    p(i) -= val;
  }
  for (i = 1; i <= n; i++)
    p(i) /= d(i);

  for (i = n; i >= 1; i--)
  {
    val = 0;
    for (j = i+1; j <= n; j++)
      val += p(j) * l(j, i);
    p(i) -= val;
  }
}

int LDLtUpdate (MATRIX & l, VECTOR & d, double a, const VECTOR & u)
{
  // Bemerkung: Es wird a aus R erlaubt
  // Rueckgabewert: 0 .. D bleibt positiv definit
  //               1 .. sonst

  int i, j, n;

  n = l.Height();

  VECTOR v(n);
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
  VECTOR & x,         // i: Startwert
                      // o: Loesung, falls IFAIL = 0
  double (*fun)(const VECTOR & x, VECTOR & g)
  )

{

  int i, j, n = x.Length();
  long it;
  char a1crit, a3acrit;


  VECTOR d(n), g(n), p(n), temp(n), bs(n), xneu(n), y(n), s(n);
  MATRIX l(n);

  double normg, alphahat, hd, fold;
  double a1, a2;
  const double mu1 = 0.1, sigma = 0.1, xi1 = 1, xi2 = 10;
  const double tau = 0.1, tau1 = 0.1, tau2 = 0.6;

  VECTOR typx(x.Length());      // i: typische Groessenordnung der Komponenten
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

  f = fun (x, g);

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

    /*
        l.Solve (g, temp);
        for (i = 1; i <= n; i++)
          temp.Elem(i) /= d.Elem(i);
        l.SolveTranspose (temp, p);
     */

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
    /*
        for (i = 1; i <= n; i++)
          s.Elem(i) = xneu.Elem(i) - x.Elem(i);
        for (i = 1; i <= n; i++)
          y.Elem(i) = g.Elem(i) - y.Elem(i);
     */
    x = xneu;

    // BFGS Update

    MultLDLt (l, d, s, bs);
    /*
        l.MultTranspose (s, temp);
        for (i = 1; i <= n; i++)
          temp.Elem(i) *= d.Elem(i);
        l.Mult (temp, bs);
     */
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
