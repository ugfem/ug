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

/* #define FAMG_SPARSE_BLOCK 1 */

#ifdef FAMG_SPARSE_BLOCK

class FAMGSparseTransfer
{
public:
  int Get_nr() const;
  int *Get_offset() const;
  int Get_offset(int i) const;
  void Construct(const class FAMGSparseBlock &sb);
private:
  int nr;
  int *offset;

};

inline int FAMGSparseTransfer::Get_nr() const {
  return nr;
}
inline int *FAMGSparseTransfer::Get_offset() const {
  return offset;
}
inline int FAMGSparseTransfer::Get_offset(int i) const {
  return offset[i];
}


class FAMGSparseBlock
{
public:
  int Get_ne() const;
  int Get_nid() const;
  int Get_nr() const;
  int Get_nc() const;
  int *Get_start() const;
  int Get_start(const int i) const;
  int *Get_index() const;
  int Get_index(const int i) const;
  int *Get_offset() const;
  int Get_offset(const int i) const;
  void Transposed(const FAMGSparseBlock &sb);
  int Product(const FAMGSparseBlock &sb);
  int CheckCGIdent(const FAMGSparseTransfer &sp, const FAMGSparseTransfer &sr);
private:
  int ne;
  int nid;
  int nr;
  int nc;
  int *start;
  int *index;
  int *offset;

  friend void FAMGTestSparseBlock();
  friend void SparseBlockProduct();
};

inline int FAMGSparseBlock::Get_ne() const {
  return ne;
}
inline int FAMGSparseBlock::Get_nid() const {
  return nid;
}
inline int FAMGSparseBlock::Get_nr() const {
  return nr;
}
inline int FAMGSparseBlock::Get_nc() const {
  return nc;
}
inline int *FAMGSparseBlock::Get_start() const {
  return start;
}
inline int FAMGSparseBlock::Get_start(const int i) const {
  return start[i];
}
inline int *FAMGSparseBlock::Get_index() const {
  return index;
}
inline int FAMGSparseBlock::Get_index(const int i) const {
  return index[i];
}
inline int *FAMGSparseBlock::Get_offset() const {
  return offset;
}
inline int FAMGSparseBlock::Get_offset(const int i) const {
  return offset[i];
}


void FAMGTestSparseBlock();

#endif



#endif
