// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef FILE_OPTI
#define FILE_OPTI

/**************************************************************************/
/* File:   opti.hh                                                        */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/*
   Optimisation functions
 */


class MinFunction
{
public:
  virtual double Func (const Vector & x) const;
  virtual void Grad (const Vector & x, Vector & g) const;
  virtual double FuncGrad (const Vector & x, Vector & g) const;
};


extern void BFGS (Vector & x, const MinFunction & fun);
void SteepestDescent (Vector & x, const MinFunction & fun);


extern void lines (
  Vector & x,         // i: Ausgangspunkt der Liniensuche
  Vector & xneu,      // o: Loesung der Liniensuche bei Erfolg
  Vector & p,         // i: Suchrichtung
  double & f,         // i: Funktionswert an der Stelle x
                      // o: Funktionswert an der Stelle xneu, falls ifail = 0
  Vector & g,         // i: Gradient an der Stelle x
                      // o: Gradient an der Stelle xneu, falls ifail = 0

  const MinFunction & fun,  // function to minmize

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





#endif
