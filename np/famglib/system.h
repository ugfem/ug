// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      system.h														*/
/*																			*/
/* Purpose:   cmg system class												*/
/*																			*/
/* Author:    Christian Wagner												*/
/*			  Institut fuer Computeranwendungen  III						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: chris@ica3.uni-stuttgart.de							*/
/*																			*/
/*																			*/
/* History:   November 97 begin, Stuttgart									*/
/*			  August 98 integration into ug (Christian Wrobel)				*/
/*																			*/
/* Remarks:																	*/
/*																			*/
/****************************************************************************/

#ifndef __CMG_SYSTEM__
#define __CMG_SYSTEM__

#include <string.h>
#include "matrix.h"
#include "grid.h"
#include "multigrid.h"

/* RCS_ID
   $Header$
 */

const int CMGMULTIGRIDS=1;

class CMGParameter  // struct makes it easier than class
{
public:
  int Getheap();
  int Getnv();
  int Getgamma();
  int Getn1();
  int Getn2();
  double Getcgilut();
  double Getilut();
  int Getcgnodes();
  double Getmincoarse();
  int Getconloops();
  int Gettype();
  int Getstv();
  double Gettol();
  double Getsigma();
  double Getomegar();
  double Getomegal();
  double Geterror1();
  double Geterror2();
  int Getmaxit();
  double Getalimit();
  double Getrlimit();
  double Getdivlimit();
  double Getreduction();
  char* Getsolver();
  char* Getpresmoother();
  char* Getpostsmoother();
  char* Getcgsmoother();
  int Read();
  void Setheap(int i);
  void Setnv(int i);
  void Setgamma(int i);
  void Setn1(int i);
  void Setn2(int i);
  void Setilut(double d);
  void Setcgilut(double d);
  void Setcgnodes(int i);
  void Setmincoarse(double d);
  void Setconloops(int i);
  void Settype(int i);
  void Setstv(int i);
  void Settol(double d);
  void Setsigma(double d);
  void Setomegar(double d);
  void Setomegal(double d);
  void Seterror1(double d);
  void Seterror2(double d);
  void Setmaxit(int i);
  void Setalimit(double d);
  void Setrlimit(double d);
  void Setdivlimit(double d);
  void Setreduction(double d);
  void Setsolver(char *solver);
  void Setpresmoother(char *presmoother);
  void Setpostsmoother(char *postsmoother);
  void Setcgsmoother(char *cgsmoother);
private:
  int heap;
  int nv;
  int gamma;
  int n1;
  int n2;
  double ilut;
  double cgilut;
  int cgnodes;
  double mincoarse;
  int conloops;
  int type;
  int stv;
  double tol;
  double sigma;
  double omegar;
  double omegal;
  double error1;
  double error2;
  int maxit;
  double alimit;
  double rlimit;
  double divlimit;
  double reduction;
  char solver[10];
  char presmoother[10];
  char postsmoother[10];
  char cgsmoother[10];
};

inline int CMGParameter::Getheap() {
  return heap;
}
inline int CMGParameter::Getnv() {
  return nv;
}
inline int CMGParameter::Getgamma() {
  return gamma;
}
inline int CMGParameter::Getn1() {
  return n1;
}
inline int CMGParameter::Getn2() {
  return n2;
}
inline double CMGParameter::Getilut() {
  return ilut;
}
inline double CMGParameter::Getcgilut() {
  return cgilut;
}
inline int CMGParameter::Getcgnodes() {
  return cgnodes;
}
inline double CMGParameter::Getmincoarse() {
  return mincoarse;
}
inline int CMGParameter::Getconloops() {
  return conloops;
}
inline int CMGParameter::Gettype() {
  return type;
}
inline int CMGParameter::Getstv() {
  return stv;
}
inline double CMGParameter::Gettol() {
  return tol;
}
inline double CMGParameter::Getsigma() {
  return sigma;
}
inline double CMGParameter::Getomegar() {
  return omegar;
}
inline double CMGParameter::Getomegal() {
  return omegal;
}
inline double CMGParameter::Geterror1() {
  return error1;
}
inline double CMGParameter::Geterror2() {
  return error2;
}
inline int CMGParameter::Getmaxit() {
  return maxit;
}
inline double CMGParameter::Getalimit() {
  return alimit;
}
inline double CMGParameter::Getrlimit() {
  return rlimit;
}
inline double CMGParameter::Getdivlimit() {
  return divlimit;
}
inline double CMGParameter::Getreduction() {
  return reduction;
}
inline char* CMGParameter::Getsolver() {
  return solver;
}
inline char* CMGParameter::Getpresmoother() {
  return presmoother;
}
inline char* CMGParameter::Getpostsmoother() {
  return postsmoother;
}
inline char* CMGParameter::Getcgsmoother() {
  return cgsmoother;
}

inline void CMGParameter::Setheap(int i) {
  heap = i;
}
inline void CMGParameter::Setnv(int i) {
  nv = i;
}
inline void CMGParameter::Setgamma(int i) {
  gamma = i;
}
inline void CMGParameter::Setilut(double d) {
  ilut = d;
}
inline void CMGParameter::Setcgilut(double d) {
  cgilut = d;
}
inline void CMGParameter::Setn1(int i) {
  n1 = i;
}
inline void CMGParameter::Setn2(int i) {
  n2 = i;
}
inline void CMGParameter::Setcgnodes(int i) {
  cgnodes = i;
}
inline void CMGParameter::Setmincoarse(double d) {
  mincoarse = d;
}
inline void CMGParameter::Setconloops(int i) {
  conloops = i;
}
inline void CMGParameter::Settype(int i) {
  type = i;
}
inline void CMGParameter::Setstv(int i) {
  stv = i;
}
inline void CMGParameter::Settol(double d) {
  tol = d;
}
inline void CMGParameter::Setsigma(double d) {
  sigma = d;
}
inline void CMGParameter::Setomegar(double d) {
  omegar = d;
}
inline void CMGParameter::Setomegal(double d) {
  omegal = d;
}
inline void CMGParameter::Seterror1(double d) {
  error1 = d;
}
inline void CMGParameter::Seterror2(double d) {
  error2 = d;
}
inline void CMGParameter::Setmaxit(int i) {
  maxit = i;
}
inline void CMGParameter::Setalimit(double d) {
  alimit = d;
}
inline void CMGParameter::Setrlimit(double d) {
  rlimit = d;
}
inline void CMGParameter::Setdivlimit(double d) {
  divlimit = d;
}
inline void CMGParameter::Setreduction(double d) {
  reduction = d;
}
inline void CMGParameter::Setsolver(char *ptr) {
  strcpy(solver,ptr);
}
inline void CMGParameter::Setpresmoother(char *ptr) {
  strcpy(presmoother,ptr);
}
inline void CMGParameter::Setpostsmoother(char *ptr) {
  strcpy(postsmoother,ptr);
}
inline void CMGParameter::Setcgsmoother(char *ptr) {
  strcpy(cgsmoother,ptr);
}


class CMGSystem
{
public:
  CMGSystem();
  CMGMatrix * GetMatrix() const;
  int GetN() const;
  double *GetVector(int i) const;
  void **GetExtra() const;
  CMGMultiGrid *GetMultiGrid(int) const;
  void SetMatrix(CMGMatrix *);
  void SetExtra(void **);
  void SetVector(int, double *);
  void SetN(int);
  int Init();
  CMGMultiGrid *CreateMultiGrid();
  int Solve(double *rhs, double *defect, double *unknown);
  int LinIt();
  int AdTVSolve();
  int BiCGStab();
  int BiCG();
  int Arnoldi(CMGMultiGrid *mg0, double **vec, double *H, double *G, double *Q, double *P, double &q0, int con);
  int ArnoldiTrans(CMGMultiGrid *mg0, double **vec, double *H, double *G, double *Q, double *P, double &q0, int con);
  int UpdateSolution(CMGMultiGrid *mg0, double **vec, double *H, double *Q, double &q0, int con);
  int ComputeEigenVector(CMGMultiGrid *mg0, double **vec, double *G, double *P, int con);
  int ComputeEigenVectorTrans(CMGMultiGrid *mg0, double **vec, double *G, double *P, int con);
  int GMRES();
  int Construct(double *entr, int *index, int *start, int n, int nl, double *tvA, double *tvB, void **extraptr);
  int ConstructSimple(double *entr, int *index, int *start, int n, int nl, void **extraptr);
  int Deconstruct();
  int DeconstructSimple();
private:
  int nmg;
  int n;             // unknowns
  CMGMultiGrid *mg[CMGMULTIGRIDS];
  CMGMatrix *matrix;
  double *vector[CMGMAXVECTORS];
  int *colmap;
  void **extra;
  int (CMGSystem::*SolverPtr)(void);

};

inline CMGMatrix *CMGSystem::GetMatrix() const {
  return matrix;
}
inline int CMGSystem::GetN() const {
  return n;
}
inline double *CMGSystem::GetVector(int i) const {
  return vector[i];
}
inline void **CMGSystem::GetExtra() const {
  return extra;
}
inline CMGMultiGrid *CMGSystem::GetMultiGrid(int i) const {
  return mg[i];
}
inline void CMGSystem::SetN(int i) {
  n = i;
}
inline void CMGSystem::SetMatrix(CMGMatrix *ptr) {
  matrix = ptr;
}
inline void CMGSystem::SetExtra(void **ptr) {
  extra = ptr;
}
inline void CMGSystem::SetVector(int i, double *ptr) {
  vector[i] = ptr;
}

CMGParameter * CMGGetParameter();
void CMGSetParameter(CMGParameter *ptr);

#endif
