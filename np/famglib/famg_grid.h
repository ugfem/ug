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

#include "famg_algebra.h"
#include "famg_decomp.h"
#include "famg_graph.h"
#include "famg_transfer.h"

/* RCS_ID
   $Header$
 */

//additional measurement and output
#define XFERTIMING

const int FAMGRHS=0;
const int FAMGUNKNOWN=1;
const int FAMGDEFECT=2;
const int FAMGTVA=3;
const int FAMGTVB=4;
const int FAMGMAXVECTORS=5;


class FAMGTransfer;

class FAMGGrid
{
public:
  ~FAMGGrid();
  int GetN() const {
    return n;
  };
  FAMGMatrixAlg *GetMatrix() const;
#ifdef FAMG_SPARSE_BLOCK
  FAMGMatrixAlg *GetDiagMatrix() const;
  void SetDiagMatrix(FAMGMatrixAlg *mat);
  int ConstructDiagonalSum();
  int ConstructDiagonalLump();
  int ConstructDiagonalInvLump();
  int ConstructDiagonalTV();
  int ConstructDiagonalInvTV();
  int ConstructDiagonalInverse();
  int ConstructSparseBlockStructure(FAMGGrid *fg);
#endif
  FAMGMatrixAlg *GetTmpMatrix() const;
  FAMGDecomp *GetDecomp() const;
  FAMGGraph *GetGraph() const;
  FAMGGridVector &GetGridVector() const;
  FAMGVector *GetVector(int) const;
  FAMGVector **GetVectorPtr();
  FAMGTransfer *GetTransfer() const;
  int GetNF() const;
  int* GetMap() const;
  void SetMatrix(FAMGMatrixAlg *mat);
  void SetTmpMatrix(FAMGMatrixAlg *tmp);
  void SetVector(int, FAMGVector *);

  void Defect( void ) const;
  void Defect(FAMGVector &defect, FAMGVector &rhs, FAMGVector &solution) const;
  void JACSmooth();
  void SGSSmooth();
  void FGSSmooth();
  void BGSSmooth();
  void ILUTSmooth();
  void JacobiSmoothFG();
  void Prolongation(const FAMGGrid *cg, const FAMGVector &cgsol, FAMGVector &fgsol, FAMGVector &fgdefect, FAMGVector *c);
  void Restriction(FAMGVector &fgsolution, FAMGVector &fgdefect, FAMGVector &cgdefect, FAMGVector &fgglobsolution) const;
  int ConstructTransfer();
  int AnalyseParents(int i);
  int InitLevel0(const class FAMGSystem &);
  int Init(int n, const FAMGGrid &grid_pattern);
  int Construct(FAMGGrid *);
  void Deconstruct( int level);
  int SolveCoarseGrid();
  int BiCGStab();
  void SmoothTV();
  int AnalyseNodeSimple(FAMGNode* nodei, FAMGPaList *&palist);

#ifdef USE_UG_DS
  GRID *GetugGrid() const {
    return mygrid;
  }
#else
  int* GetFather() const;
#endif

#ifdef FAMG_ILU
  int ILUTDecomp(int);
  int OrderVector(int,int*);
  int ReorderVector(int,int*);
  int Order(int*);
  int Reorder();
#endif

#ifdef UG_DRAW
  void **GetNode() const;
#endif


  int AnalyseNode0(const FAMGVectorEntry &veci, FAMGPaList *&palist);
  double BestFirst0(FAMGPaList *&palist, double mii, double prt, double plt, double ff, double gg, struct FAMGSpecialData *sd, int nnb);
  double BestSecond0(FAMGPaList *&palist, double mii, double prt, double plt, double ff, double gg, struct FAMGSpecialData *sd, int nnb);
  double BestFirst2(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData sd);
  void FF0(const FAMGVectorEntry &veci, double &ff, double &gg, double *f, double *g);
  void DFF0(const FAMGVectorEntry &veci);
  void JK0(int k, const FAMGVectorEntry &veck, double &hjk, double &ejk, double* h, double *e);
  int LocalJ0(int j, const FAMGVectorEntry &vecj, double* h, double *e);
  void DLocalJ0(int j, const FAMGVectorEntry &vecj);
  void FJ0(int j, const FAMGVectorEntry &vecj, double &fj, double &gj, double *f, double *g);
  void JJ0(const FAMGVectorEntry &vecj, double &hjj, double &ejj);

  int AnalyseNode1(const FAMGVectorEntry &veci, FAMGPaList *&palist);
  double BestFirst1(FAMGPaList *&palist, double mii, double prt, double plt, double ff, double gg, struct FAMGSpecialData *sd, int nnb);
  double BestFirst5(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData sd);
  double BestSecond1(FAMGPaList *&palist, double mii, double prt, double plt, double ff, double gg, struct FAMGSpecialData *sd, int nnb);
  void FF1(const FAMGVectorEntry &veci, double &ff, double &gg, double *f, double *g);
  void DFF1(const FAMGVectorEntry &veci);
  void JK1(int k, const FAMGVectorEntry &veck, double &hjk, double &ejk, double* h, double *e);
  int LocalJ1(int j, const FAMGVectorEntry &vecj, double* h, double *e);
  void DLocalJ1(int j, const FAMGVectorEntry &vecj);
  void FJ1(int j, const FAMGVectorEntry &vecj, double &fj, double &gj, double *f, double *g);
  void JJ1(const FAMGVectorEntry &vecj, double &hjj, double &ejj);

  int AnalyseNode2(const FAMGVectorEntry &veci, FAMGPaList *&palist);
  int AnalyseNode3(const FAMGVectorEntry &veci, FAMGPaList *&palist);
  int AnalyseNode4(const FAMGVectorEntry &veci, FAMGPaList *&palist);
  int AnalyseNode5(const FAMGVectorEntry &veci, FAMGPaList *&palist);

  int ConstructLocalMatrix(const FAMGVectorEntry &veci, struct FAMGMatrixLocal *localmatrix, double *&tv1, double *&tv2, double &normr, double &norml);
  int AnalyseNode6(const FAMGVectorEntry &veci, FAMGPaList *&palist);
  int GetLocalMinimum1(FAMGPaList *&palist, double *w1, double *w2, double *t1, double *t2, struct FAMGMatrixLocal *localmatrix);
  int GetLocalMinimum2(FAMGPaList *&palist, double *w1, double *w2, double *t1, double *t2, struct FAMGMatrixLocal *localmatrix);

  int ConstructLocalMatrixA(const FAMGVectorEntry &veci, struct FAMGMatrixLocal *localmatrix, double *&tv1, double *&tv2, double &normr, double &norml);
  int AnalyseNode17(const FAMGVectorEntry &veci, FAMGPaList *&palist);
  int GetLocalMinimum1A(FAMGPaList *&palist, double *w1, double *w2, double *t1, double *t2, struct FAMGMatrixLocal *localmatrix);
  int GetLocalMinimum2A(FAMGPaList *&palist, double *w1, double *w2, double *t1, double *t2, struct FAMGMatrixLocal *localmatrix);
  int GetLocalMinimum2B(FAMGPaList *&palist, double *w1, double *w2, double *t1, double *t2, struct FAMGMatrixLocal *localmatrix);

  int SetFlagsAndCount(int i, int f);
  int Connected(int i, int z);
  void SetFlags(int i, int f);
#ifdef FAMG_SPARSE_BLOCK
  int SaveCoeffs(const FAMGVectorEntry& i, int np, const int pa[], double **coeff, double **coefft);
#else
  int SaveCoeffs(const FAMGVectorEntry& i, int np, const int pa[], double coeff[], double coefft[]);
#endif
  int CountLinks(int i);
  int UpdateNBNewCG(int i);
  int UpdateNeighborsCG(int i);
  void Stencil();
  void CGSmooth();
  void PreSmooth();
  void PostSmooth();
  void GetSmoother();
  FAMGMatrixAlg *GetConsMatrix() const;                 // in serial synonym for the stiffness matrix
  // in parallel the (partly) consistent stiffness matrix
  void SetConsMatrix(FAMGMatrixAlg *cm);
#ifdef ModelP
  int CommunicateNodeStatus();
  void ConstructOverlap();
  INT GetNrMasterVectors() const {
    return NVEC_PRIO(GetugGrid(),PrioMaster);
  }
  INT GetNrBorderVectors() const {
    return NVEC_PRIO(GetugGrid(),PrioBorder);
  }
  INT GetNrGhostVectors() const {
    return NVEC_PRIO(GetugGrid(),PrioVGhost)+NVEC_PRIO(GetugGrid(),PrioVHGhost)+NVEC_PRIO(GetugGrid(),PrioHGhost);
  }
  GRID *GetugGrid() {
    return mygrid;
  }
#endif

private:

#ifdef USE_UG_DS
  void SetugGrid(GRID *grid) {
    mygrid=grid;
  }
#else
  void SetFather(int *f) {
    father=f;
  }
#endif

  int n;                                        // number unknowns
  int nf;                                                       // number fine grid unknowns
  FAMGMatrixAlg *matrix;                        // stiffness matrix
  FAMGMatrixAlg *Consmatrix;                    // (partly) consistent stiffness matrix
  FAMGMatrixAlg *tmpmatrix;                     // temp. stiffness matrix for a double-step todo: remove
#ifdef FAMG_SPARSE_BLOCK
  FAMGMatrixAlg *diagmatrix;
#endif
  FAMGTransfer *transfer;                       // transfer matrix
  FAMGGraph *graph;                                     // node graph for elemination
  FAMGGridVector *mygridvector;
  FAMGVector *vector[FAMGMAXVECTORS];

#ifdef USE_UG_DS
  GRID *mygrid;
#else
  int *father;
#endif

#ifdef FAMG_ILU
  FAMGDecomp *decomp;
  int *map;        // actually either map or father
#endif

#ifdef UG_DRAW
  void  **vertex;
#endif
  void (FAMGGrid::*CGSmootherPtr)(void);
  void (FAMGGrid::*PreSmootherPtr)(void);
  void (FAMGGrid::*PostSmootherPtr)(void);
};

inline FAMGMatrixAlg *FAMGGrid::GetMatrix() const {
  return matrix;
}
inline FAMGMatrixAlg *FAMGGrid::GetTmpMatrix() const {
  return tmpmatrix;
}
inline FAMGGraph *FAMGGrid::GetGraph() const {
  return graph;
}
inline FAMGGridVector &FAMGGrid::GetGridVector() const {
  return *mygridvector;
}
inline FAMGVector **FAMGGrid::GetVectorPtr()  {
  return vector;
}
inline FAMGVector *FAMGGrid::GetVector(int i) const {
  return vector[i];
}
inline int FAMGGrid::GetNF() const {
  return nf;
}
inline FAMGTransfer *FAMGGrid::GetTransfer() const {
  return transfer;
}
inline void FAMGGrid::SetMatrix(FAMGMatrixAlg *mat) {
  matrix = mat;
}
inline void FAMGGrid::SetTmpMatrix(FAMGMatrixAlg *tmp) {
  tmpmatrix = tmp;
}
inline void FAMGGrid::SetVector(int i, FAMGVector *p) {
  vector[i] = p;
}
inline void FAMGGrid::SetConsMatrix(FAMGMatrixAlg *cm) {
  Consmatrix = cm;
}
inline FAMGMatrixAlg *FAMGGrid::GetConsMatrix() const {
  return Consmatrix;
}                                                                               // in parallel the (partly) consistent stiffness matrix
inline FAMGGrid::~FAMGGrid() {
  delete GetGraph(); delete &GetGridVector();
  delete GetMatrix(); delete GetConsMatrix();
  for( int i=0; i<FAMGMAXVECTORS; i++) delete GetVector(i);
};
inline void FAMGGrid::Defect(FAMGVector &defect, FAMGVector &rhs, FAMGVector &solution) const
{
  defect.VecMinusMatVec(rhs, *GetMatrix(), solution);
}
inline void FAMGGrid::Defect( void ) const
{       // take default values
  Defect(*GetVector(FAMGDEFECT), *GetVector(FAMGRHS), *GetVector(FAMGUNKNOWN));
}
#ifdef FAMG_SPARSE_BLOCK
inline FAMGMatrixAlg *FAMGGrid::GetDiagMatrix() const {
  return diagmatrix;
}
inline void FAMGGrid::SetDiagMatrix(FAMGMatrixAlg *mat) {
  diagmatrix = mat;
}
#endif

#ifdef USE_UG_DS
#else
inline int* FAMGGrid::GetFather() const {
  return father;
}
#endif

#ifdef FAMG_ILU
inline FAMGDecomp *FAMGGrid::GetDecomp() const {
  return decomp;
}
inline int* FAMGGrid::GetMap() const {
  return map;
}
#endif

#ifdef UG_DRAW
inline void **FAMGGrid::GetNode() const {
  return vertex;
}
#endif

// only for debugging printv printm printim
void prv( int level, int x_nr );
void prim(int level);
void prm(int level, int comp = 0);
void prvGeom( int level, int x_nr );
void primGeom(int level);
void prmGeom(int level, int comp = 0);
void printlist( FAMGGraph *graph);
void printlistGeom( FAMGGraph *graph);
#endif
