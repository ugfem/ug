// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      famg_grid.h													*/
/*																			*/
/* Purpose:   famg grid classes												*/
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

#ifndef __FAMG_GRID__
#define __FAMG_GRID__

#include "famg_matrix.h"
#include "famg_decomp.h"
#include "famg_transfer.h"
#include "famg_graph.h"

/* RCS_ID
   $Header$
 */

const int FAMGRHS=0;
const int FAMGUNKNOWN=1;
const int FAMGDEFECT=2;
const int FAMGTVA=3;
const int FAMGTVB=4;
const int FAMGMAXVECTORS=5;

class FAMGGrid
{
public:
  void **GetNode() const;
  FAMGMatrix *GetMatrix() const;
  FAMGMatrix *GetTmpMatrix() const;
  FAMGDecomp *GetDecomp() const;
  FAMGGraph *GetGraph() const;
  double *GetVector(int) const;
  double **GetVectorPtr();
  FAMGTransfer *GetTransfer() const;
  int GetN() const;
  int GetNF() const;
  int* GetFather() const;
  int* GetMap() const;
  void SetTmpMatrix(FAMGMatrix *tmp);
  void SetVector(int, double *);

  void Defect();
  void DefectTrans();
  void JACSmooth();
  void SGSSmooth();
  void FGSSmooth();
  void BGSSmooth();
  void ILUTSmooth();
  void DevideFGDefect();
  void Prolongation(const FAMGGrid *cg);
  void ProlongationTrans(const FAMGGrid *cg);
  void Restriction(FAMGGrid *cg) const;
  void RestrictionTrans(FAMGGrid *cg) const;
  int ILUTDecomp(int);
  int ConstructTransfer();
  void CopyVector(int, int);
  void AddVector(int, int);
  void SubVector(int, int);
  void SetVector(int, double);
  void MultVector(int source, double factor);
  int AnalyseParents(int i);
  int InitLevel0(const class FAMGSystem &);
  int Init(int);
  int Construct(FAMGGrid *);
  void Deconstruct();
  int BiCGStab();
  void SmoothTV();
  int OrderVector(int,int*);
  int ReorderVector(int,int*);
  int Order(int*);
  int Reorder();
  int AnalyseNodeSimple(int i, FAMGPaList *&palist);


  int AnalyseNode0(int i, FAMGPaList *&palist);
  double BestFirst0(FAMGPaList *&palist, double mii, double prt, double plt, double ff, double gg, struct FAMGSpecialData *sd, int nnb);
  double BestSecond0(FAMGPaList *&palist, double mii, double prt, double plt, double ff, double gg, struct FAMGSpecialData *sd, int nnb);
  double BestFirst2(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData sd);
  void FF0(int j, double &ff, double &gg, double *f, double *g);
  void DFF0(int j);
  void JK0(int k, double &hjk, double &ejk, double* h, double *e);
  int LocalJ0(int j, double* h, double *e);
  void DLocalJ0(int j);
  void FJ0(int j, double &fj, double &gj, double *f, double *g);
  void JJ0(int j, double &hjj, double &ejj);

  int AnalyseNode1(int i, FAMGPaList *&palist);
  double BestFirst1(FAMGPaList *&palist, double mii, double prt, double plt, double ff, double gg, struct FAMGSpecialData *sd, int nnb);
  double BestFirst5(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData sd);
  double BestSecond1(FAMGPaList *&palist, double mii, double prt, double plt, double ff, double gg, struct FAMGSpecialData *sd, int nnb);
  void FF1(int j, double &ff, double &gg, double *f, double *g);
  void DFF1(int j);
  void JK1(int k, double &hjk, double &ejk, double* h, double *e);
  int LocalJ1(int j, double* h, double *e);
  void DLocalJ1(int j);
  void FJ1(int j, double &fj, double &gj, double *f, double *g);
  void JJ1(int j, double &hjj, double &ejj);

  int AnalyseNode2(int i, FAMGPaList *&palist);
  int AnalyseNode3(int i, FAMGPaList *&palist);
  int AnalyseNode4(int i, FAMGPaList *&palist);
  int AnalyseNode5(int i, FAMGPaList *&palist);

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
  FAMGMatrix *matrix;
  FAMGMatrix *tmpmatrix;
  FAMGDecomp *decomp;
  FAMGTransfer *transfer;
  FAMGGraph *graph;
  double *vector[FAMGMAXVECTORS];
  int *father;
  int *map;        // actually either map or father
  void  **vertex;
  void (FAMGGrid::*CGSmootherPtr)(void);
  void (FAMGGrid::*PreSmootherPtr)(void);
  void (FAMGGrid::*PostSmootherPtr)(void);
};

inline void **FAMGGrid::GetNode() const {
  return vertex;
}
inline FAMGMatrix *FAMGGrid::GetMatrix() const {
  return matrix;
}
inline FAMGMatrix *FAMGGrid::GetTmpMatrix() const {
  return tmpmatrix;
}
inline FAMGDecomp *FAMGGrid::GetDecomp() const {
  return decomp;
}
inline FAMGGraph *FAMGGrid::GetGraph() const {
  return graph;
}
inline double **FAMGGrid::GetVectorPtr()  {
  return vector;
}
inline double *FAMGGrid::GetVector(int i) const {
  return vector[i];
}
inline int FAMGGrid::GetN() const {
  return n;
}
inline int FAMGGrid::GetNF() const {
  return nf;
}
inline int* FAMGGrid::GetFather() const {
  return father;
}
inline int* FAMGGrid::GetMap() const {
  return map;
}
inline FAMGTransfer *FAMGGrid::GetTransfer() const {
  return transfer;
}
inline void FAMGGrid::SetTmpMatrix(FAMGMatrix *tmp) {
  tmpmatrix = tmp;
}
inline void FAMGGrid::SetVector(int i, double *p) {
  vector[i] = p;
}

#endif
