// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/***************************************************************************/
/*                                                                         */
/* Problem:        Liniensuche                                             */
/*                                                                         */
/* Programmautor:  Joachim Schöberl                                        */
/* Matrikelnummer: 9155284                                                 */
/*                                                                         */
/* Algorithmus nach:                                                       */
/*                                                                         */
/*   Optimierung I, Gfrerer, WS94/95                                       */
/*   Algorithmus 2.1: Liniensuche Problem (ii)                             */
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

const double eps0 = 1E-15;

// Liniensuche


double MinFunction :: Func (const Vector & /* x */) const
{
  cerr << "Func of MinFunction called" << endl;
  return 0;
}

void MinFunction :: Grad (const Vector & /* x */, Vector & /* g */) const
{
  cerr << "Grad of MinFunction called" << endl;
}

double MinFunction :: FuncGrad (const Vector & /* x */, Vector & /* g */) const
{
  cerr << "Grad of MinFunction called" << endl;
  return 0;
}


void lines (
  Vector & x,         // i: Ausgangspunkt der Liniensuche
  Vector & xneu,      // o: Loesung der Liniensuche bei Erfolg
  Vector & p,         // i: Suchrichtung
  double & f,         // i: Funktionswert an der Stelle x
                      // o: Funktionswert an der Stelle xneu, falls ifail = 0
  Vector & g,         // i: Gradient an der Stelle x
                      // o: Gradient an der Stelle xneu, falls ifail = 0
  const MinFunction & fun,  // function to minimize
  //  double (*fun) (const Vector & x, Vector & g),
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
  int & ifail)        // o:  0 bei erfolgreicher Liniensuche
                      //    -1 bei Abbruch wegen Unterschreiten von fmin
                      //     1 bei Abbruch, aus sonstigen Gründen

{
  double phi0, phi0prime, phi1, phi1prime, phihatprime;
  double alpha1, alpha2, alphaincr, c;
  char flag = 1;
  long it;

  alpha1 = 0;
  alpha2 = 1e50;
  phi0 = phi1 = f;

  phi0prime = g * p;

  if (phi0prime > 0)
  {
    ifail = 1;
    return;
  }

  phi1prime = phi0prime;

  it = 100000l;

  while (it-- >= 0)
  {
    xneu.Set (1, x, alphahat, p);

    f = fun.FuncGrad (xneu, g);

    if (f < fmin)
    {
      ifail = -1;
      break;
    }

    if (alpha2 - alpha1 < eps0 * alpha2)
    {
      ifail = 0;
      break;
    }

    if (f - phi0 > mu1 * alphahat * phi1prime + eps0 * fabs (phi0))

    {
      flag = 0;
      alpha2 = alphahat;

      c = (f - phi1 - phi1prime * (alphahat-alpha1)) / sqr (alphahat-alpha1);

      alphahat = alpha1 - 0.5 * phi1prime / c;

      if (alphahat > alpha2)
        alphahat = alpha1 + 1/(4*c) * ( (sigma+mu1) * phi0prime - 2*phi1prime
                                        + sqrt (sqr(phi1prime - mu1 * phi0prime) -
                                                4 * (phi1 - phi0 - mu1 * alpha1 * phi0prime) * c));

      alphahat = max (alphahat, alpha1 + tau * (alpha2 - alpha1));
      alphahat = min (alphahat, alpha2 - tau * (alpha2 - alpha1));
    }

    else

    {
      //      dfunk (xneu, f, g, funk);
      f = fun.FuncGrad (xneu, g);

      phihatprime = g * p;

      if (phihatprime < sigma * phi0prime * (1 + eps0))

      {
        if (phi1prime < phihatprime)   // Approximationsfunktion ist konvex

          alphaincr = (alphahat - alpha1) * phihatprime /
                      (phi1prime - phihatprime);

        else
          alphaincr = MAXDOUBLE;

        if (flag)
        {
          alphaincr = max (alphaincr, xi1 * (alphahat-alpha1));
          alphaincr = min (alphaincr, xi2 * (alphahat-alpha1));
        }
        else
        {
          alphaincr = max (alphaincr, tau1 * (alpha2 - alphahat));
          alphaincr = min (alphaincr, tau2 * (alpha2 - alphahat));
        }

        alpha1 = alphahat;
        alphahat += alphaincr;
        phi1 = f;
        phi1prime = phihatprime;
      }

      else

      {
        ifail = 0;     // Erfolg !!
        break;
      }
    }

  }

  if (it < 0)
    ifail = 1;
}



void SteepestDescent (Vector & x, const MinFunction & fun)
{
  int it, n = x.Length();
  Vector xnew(n), p(n), g(n), g2(n);
  double val, alphahat;
  int fail;

  val = fun.FuncGrad(x, g);

  alphahat = 1;
  for (it = 0; it < 10; it++)
  {
    p = -g;
    lines (x, xnew, p, val, g, fun, alphahat, -1e5,
           0.1, 0.1, 1, 10, 0.1, 0.1, 0.6, fail);

    x = xnew;
  }
}
