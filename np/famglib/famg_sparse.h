// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      famg_sparse.h			                                                                        */
/*																			*/
/* Purpose:   famg sparse point block classes                               */
/*																			*/
/* Author:    Christian Wagner												*/
/*			  IWR Technische Simulation                                                             */
/*			  Universitaet Heidelberg										*/
/*			  INF 368                                                                                               */
/*			  D- 69120 Heidelberg	                                                                        */
/*			  internet: Christian.Wagner@iwr.uni-heidelberg.de	                */
/*																			*/
/*																			*/
/* History:   June 99 begin, Heidelberg                                                                 */
/*			                                                                                */
/*																			*/
/* Remarks:					                                                                                        */
/*																			*/
/****************************************************************************/


#ifndef __FAMG_SPARSE__
#define __FAMG_SPARSE__

// #define FAMG_SPARSE_BLOCK 1

#ifdef FAMG_SPARSE_BLOCK
extern "C"
{
#include "sm.h"
}


class FAMGSparseVector
{
public:
  FAMGSparseVector();
  ~FAMGSparseVector();
  FAMGSparseVector(const FAMGSparseVector &sv);
  FAMGSparseVector(const short);
  FAMGSparseVector(const short, const short *);
  FAMGSparseVector(const short *, const short *, const short);
  FAMGSparseVector &operator=(const FAMGSparseVector &sb);
  void Init(const FAMGSparseVector *sv, int off);
  void Construct(const FAMGSparseVector *sv);
  void Construct(const short, const short *);
  void Product(const class FAMGSparseBlock *sb, const FAMGSparseVector *sv);
  short Get_n() const;
  short *Get_comp() const;
  short Get_comp(const short i) const;
  short Get_maxcomp() const;
  void ScalProdConstruct(const class FAMGSparseBlock *sb);
  void ScalProdConstruct(const class FAMGSparseBlock *sb1, const class FAMGSparseBlock *sb2);
  void ConstructSparseTransfer(const FAMGSparseVector *sv1, const FAMGSparseVector *sv2, const FAMGSparseVector *sv3, const FAMGSparseVector *sv4);
private:
  short n;
  short *comp;
  short maxcomp;
};

inline short FAMGSparseVector::Get_n() const {
  return n;
}
inline short *FAMGSparseVector::Get_comp() const {
  return comp;
}
inline short FAMGSparseVector::Get_comp(const short i) const {
  return comp[i];
}
inline short FAMGSparseVector::Get_maxcomp() const {
  return maxcomp;
}


class FAMGSparseBlock
{
public:
  FAMGSparseBlock();
  ~FAMGSparseBlock();
  FAMGSparseBlock(SPARSE_MATRIX*);
  FAMGSparseBlock(const FAMGSparseBlock &sb);
  FAMGSparseBlock &operator=(const FAMGSparseBlock &sb);
  void Init(SPARSE_MATRIX* ugsm);
  short Get_ne() const;
  short Get_nr() const;
  short Get_nc() const;
  short *Get_start() const;
  short Get_start(const short i) const;
  short *Get_index() const;
  short Get_index(const short i) const;
  short *Get_offset() const;
  short Get_offset(const short i) const;
  short Get_maxoffset() const;
  int Product(const FAMGSparseBlock *sb1, const FAMGSparseBlock *sb2);
  int Product(const FAMGSparseBlock *sb, const FAMGSparseVector *sv);
  int Product(const FAMGSparseVector *sv, const FAMGSparseBlock *sb);
  int Product(const FAMGSparseVector *svl, const FAMGSparseBlock *sb, const FAMGSparseVector *svr);
  int Transposed(const FAMGSparseBlock *sb);
  void FixDiag();
  int CheckStructureforAdd(const FAMGSparseBlock *sb) const;
private:
  short ne;
  short nr;
  short nc;
  short *start;
  short *index;
  short *offset;
  short maxoffset;

  friend void AdaptStructure(FAMGSparseBlock *sb1, FAMGSparseBlock *sb2);
  friend void FAMGTestSparseBlock();
};

inline short FAMGSparseBlock::Get_ne() const {
  return ne;
}
inline short FAMGSparseBlock::Get_nr() const {
  return nr;
}
inline short FAMGSparseBlock::Get_nc() const {
  return nc;
}
inline short *FAMGSparseBlock::Get_start() const {
  return start;
}
inline short FAMGSparseBlock::Get_start(const short i) const {
  return start[i];
}
inline short *FAMGSparseBlock::Get_index() const {
  return index;
}
inline short FAMGSparseBlock::Get_index(const short i) const {
  return index[i];
}
inline short *FAMGSparseBlock::Get_offset() const {
  return offset;
}
inline short FAMGSparseBlock::Get_offset(const short i) const {
  return offset[i];
}
inline short FAMGSparseBlock::Get_maxoffset() const {
  return maxoffset;
}

void SparseBlockVSet(const FAMGSparseVector *sv, double *v, double val);
void SparseBlockMSet(const FAMGSparseBlock *sb, double *a, double val);
void SparseBlockMInvertDiag(const FAMGSparseBlock *sb, double *dest, double *source);
double SparseBlockMNorm(const FAMGSparseBlock *sb,double *a);
void SparseBlockMSetDiag(const FAMGSparseBlock *sb, double *a, double val);
int SparseBlockMMAdd(const FAMGSparseBlock *sbdest, const FAMGSparseBlock *sbsource, double *dest, const double *source);
int SparseBlockMMAdd(const FAMGSparseBlock *sbdest, const FAMGSparseBlock *sbsource, double *dest, const double *source, double factor);
int SparseBlockMMProduct(const FAMGSparseBlock *sp, const FAMGSparseBlock *sb1, const FAMGSparseBlock *sb, double *ap, double *a1, double *a2);
int SparseBlockMMProduct(const FAMGSparseBlock *sd, const FAMGSparseVector *sv, const FAMGSparseBlock *sb, double *dest, double *d, double *a);
int SparseBlockMMProduct(const FAMGSparseBlock *sd, const FAMGSparseBlock *sb, const FAMGSparseVector *sv, double *dest, double *a, double *d);
int SparseBlockMMProduct(const FAMGSparseBlock *sd, const FAMGSparseVector *svl, const FAMGSparseBlock *sb, const FAMGSparseVector *svr, double *dest, double *dl, double *a, double *dr);
int SparseBlockMMAddProduct(const FAMGSparseBlock *sp, const FAMGSparseBlock *sb1, const FAMGSparseBlock *sb2, double *ap, double *a1, double *a2, double factor);
int SparseBlockMVAddProduct(const FAMGSparseVector *dcomp, const FAMGSparseBlock *sb, const FAMGSparseVector *scomp, double *vd, double *a, double *vs, const double factor);
int SparseBlockMVProduct(const FAMGSparseVector *dest, const FAMGSparseBlock *sb, const FAMGSparseVector *source, double *vd, const double *a, const double *vs);
int SparseBlockMVAddProduct(const FAMGSparseVector *dest, const FAMGSparseVector *sv, const FAMGSparseVector *source, double *vd, double *d, double *vs, const double factor);
int SparseBlockMVProduct(const FAMGSparseVector *dest, const FAMGSparseVector *sv, const FAMGSparseVector *source, double *vd, const double *d, const double *vs);
void SparseBlockMCopy(const FAMGSparseBlock *sd, const FAMGSparseBlock *ss, double *dest, double *source, double factor);
void SparseBlockVCopy(const FAMGSparseVector *svd, const FAMGSparseVector *svs, double *dest, double *source, double factor);
void SparseBlockVAdd(const FAMGSparseVector *svd, const FAMGSparseVector *svs1, const FAMGSparseVector *svs2, double *dest, double *source1, double *source2);
void SparseBlockVSub(const FAMGSparseVector *svd, const FAMGSparseVector *svs1, const FAMGSparseVector *svs2, double *dest, double *source1, double *source2);
void SparseBlockVAdd(const FAMGSparseVector *svd, const FAMGSparseVector *svs1, const FAMGSparseVector *svs2, double *dest, double *source1, double *source2, double factor);


void AdaptStructure(FAMGSparseBlock *sb1, FAMGSparseBlock *sb2);
void SparseBlockRowAddScalProd(const FAMGSparseVector *sp, const FAMGSparseBlock *sb, double *scalprod, double *a, double *b);
void SparseBlockRowAddScalProd(const FAMGSparseVector *sp, const FAMGSparseBlock *sb1, const FAMGSparseBlock *sb2, double *scalprod, double *a, double *b);

void SparseBlockMCopyDense(double *decomp, const FAMGSparseBlock *sb, double *matptr);
void SparseBlockDiagApprox(const FAMGSparseBlock *sb, const FAMGSparseBlock *sbd, const FAMGSparseVector *svt, double *mij, const double *miidecomp, const double factor, const double *tv);
void SparseBlockDiagApprox(const FAMGSparseBlock *sb, const FAMGSparseBlock *sbd, const FAMGSparseBlock *sbo, const FAMGSparseVector *svt, double *mij, const double *miidecomp, const double *matij, const double *tv);
void SparseBlockDiagApprox(const FAMGSparseBlock *sb, const FAMGSparseBlock *sbd, const FAMGSparseBlock *sbo, const FAMGSparseVector *svt, double *mij, const double *miidecomp, const double *matij, const double *mjjdecomp, const double *tv);
void SparseBlockDiagApproxT(const FAMGSparseBlock *sb, const FAMGSparseBlock *sbd, const FAMGSparseVector *svt, double *mij, const double *miidecomp, const double factor, const double *tv);
void SparseBlockDiagApproxT(const FAMGSparseBlock *sb, const FAMGSparseBlock *sbd, const FAMGSparseBlock *sbo, const FAMGSparseVector *svt, double *mij, const double *miidecomp, const double *matij, const double *tv);
void SparseBlockDiagApproxT(const FAMGSparseBlock *sb, const FAMGSparseBlock *sbd, const FAMGSparseBlock *sbo, const FAMGSparseVector *svt, double *mij, const double *miidecomp, const double *matij, const double *mjjdecomp, const double *tv);
void SparseBlockGalDiagApprox(const FAMGSparseBlock *sb, const FAMGSparseVector *sbr, const FAMGSparseBlock *sbd, const FAMGSparseVector *sbp, const FAMGSparseVector *svt, double *mij, const double *ris, const double *mss, const double *psj, const double *tv);
void SparseBlockGalDiagApproxT(const FAMGSparseBlock *sb, const FAMGSparseVector *sbr, const FAMGSparseBlock *sbd, const FAMGSparseVector *sbp, const FAMGSparseVector *svt, double *mij, const double *ris, const double *mss, const double *psj, const double *tv);

#endif



#endif
