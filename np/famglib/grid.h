// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      grid.h														*/
/*																			*/
/* Purpose:   cmg grid classes												*/
/*                                                                          */
/* Author:    Christian Wagner                                                              */
/*            Institut fuer Computeranwendungen  III			            */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27	                                        */
/*            70569 Stuttgart                                               */
/*            internet: chris@ica3.uni-stuttgart.de                         */
/*                                                                                                                  */
/*                                                                          */
/* History:   November 97 begin, Stuttgart                                  */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#ifndef __CMG_GRID__
#define __CMG_GRID__

#include "matrix.h"
#include "decomp.h"
#include "transfer.h"
#include "graph.h"

/* RCS_ID
   $Header$
 */

const int CMGRHS=0;
const int CMGUNKNOWN=1;
const int CMGDEFECT=2;
const int CMGTVA=3;
const int CMGTVB=4;
const int CMGMAXVECTORS=5;

class CMGGrid
{
public:
  void **GetNode() const;
  CMGMatrix *GetMatrix() const;
  CMGMatrix *GetTmpMatrix() const;
  CMGDecomp *GetDecomp() const;
  CMGGraph *GetGraph() const;
  double *GetVector(int) const;
  double **GetVectorPtr();
  CMGTransfer *GetTransfer() const;
  int GetN() const;
  int GetNF() const;
  int* GetFather() const;
  int* GetMap() const;
  void SetTmpMatrix(CMGMatrix *tmp);
  void SetVector(int, double *);

  void Defect();
  void DefectTrans();
  void JACSmooth();
  void SGSSmooth();
  void FGSSmooth();
  void BGSSmooth();
  void ILUTSmooth();
  void DevideFGDefect();
  void Prolongation(const CMGGrid *cg);
  void ProlongationTrans(const CMGGrid *cg);
  void Restriction(CMGGrid *cg) const;
  void RestrictionTrans(CMGGrid *cg) const;
  int ILUTDecomp(int);
  int ConstructTransfer();
  void CopyVector(int, int);
  void AddVector(int, int);
  void SubVector(int, int);
  void SetVector(int, double);
  void MultVector(int source, double factor);
  int AnalyseParents(int i);
  int InitLevel0(const class CMGSystem &);
  int Init(int);
  int Construct(CMGGrid *);
  void Deconstruct();
  int BiCGStab();
  void SmoothTV();
  int OrderVector(int,int*);
  int ReorderVector(int,int*);
  int Order(int*);
  int Reorder();
  int AnalyseNodeSimple(int i, CMGPaList *&palist);


  int AnalyseNode0(int i, CMGPaList *&palist);
  double BestFirst0(CMGPaList *&palist, double mii, double prt, double plt, double ff, double gg, struct CMGSpecialData *sd, int nnb);
  double BestSecond0(CMGPaList *&palist, double mii, double prt, double plt, double ff, double gg, struct CMGSpecialData *sd, int nnb);
  double BestFirst2(CMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, CMGSpecialData sd);
  void FF0(int j, double &ff, double &gg, double *f, double *g);
  void DFF0(int j);
  void JK0(int k, double &hjk, double &ejk, double* h, double *e);
  int LocalJ0(int j, double* h, double *e);
  void DLocalJ0(int j);
  void FJ0(int j, double &fj, double &gj, double *f, double *g);
  void JJ0(int j, double &hjj, double &ejj);

  int AnalyseNode1(int i, CMGPaList *&palist);
  double BestFirst1(CMGPaList *&palist, double mii, double prt, double plt, double ff, double gg, struct CMGSpecialData *sd, int nnb);
  double BestFirst5(CMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, CMGSpecialData sd);
  double BestSecond1(CMGPaList *&palist, double mii, double prt, double plt, double ff, double gg, struct CMGSpecialData *sd, int nnb);
  void FF1(int j, double &ff, double &gg, double *f, double *g);
  void DFF1(int j);
  void JK1(int k, double &hjk, double &ejk, double* h, double *e);
  int LocalJ1(int j, double* h, double *e);
  void DLocalJ1(int j);
  void FJ1(int j, double &fj, double &gj, double *f, double *g);
  void JJ1(int j, double &hjj, double &ejj);

  int AnalyseNode2(int i, CMGPaList *&palist);
  int AnalyseNode3(int i, CMGPaList *&palist);
  int AnalyseNode4(int i, CMGPaList *&palist);
  int AnalyseNode5(int i, CMGPaList *&palist);

  int SetFlagsAndCount(int i, int f);
  int Connected(int i, int z);
  void SetFlags(int i, int f);
  int SaveCoeffs(int i, int np, int *pa, double *coeff, double *coefft);
  int CountLinks(int i);
  int UpdateNBNewCG(int i);
  int UpdateNeighborsCG(int i);
  void Stencil();
  void CGSmooth();
  void PreSmooth();
  void PostSmooth();
  void GetSmoother();
private:
  int n;                     // unknowns
  int nf;
  CMGMatrix *matrix;
  CMGMatrix *tmpmatrix;
  CMGDecomp *decomp;
  CMGTransfer *transfer;
  CMGGraph *graph;
  double *vector[CMGMAXVECTORS];
  int *father;
  int *map;        // actually either map or father
  void  **vertex;
  void (CMGGrid::*CGSmootherPtr)(void);
  void (CMGGrid::*PreSmootherPtr)(void);
  void (CMGGrid::*PostSmootherPtr)(void);
};

inline void **CMGGrid::GetNode() const {
  return vertex;
}
inline CMGMatrix *CMGGrid::GetMatrix() const {
  return matrix;
}
inline CMGMatrix *CMGGrid::GetTmpMatrix() const {
  return tmpmatrix;
}
inline CMGDecomp *CMGGrid::GetDecomp() const {
  return decomp;
}
inline CMGGraph *CMGGrid::GetGraph() const {
  return graph;
}
inline double **CMGGrid::GetVectorPtr()  {
  return vector;
}
inline double *CMGGrid::GetVector(int i) const {
  return vector[i];
}
inline int CMGGrid::GetN() const {
  return n;
}
inline int CMGGrid::GetNF() const {
  return nf;
}
inline int* CMGGrid::GetFather() const {
  return father;
}
inline int* CMGGrid::GetMap() const {
  return map;
}
inline CMGTransfer *CMGGrid::GetTransfer() const {
  return transfer;
}
inline void CMGGrid::SetTmpMatrix(CMGMatrix *tmp) {
  tmpmatrix = tmp;
}
inline void CMGGrid::SetVector(int i, double *p) {
  vector[i] = p;
}

#endif
