// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      graph.h														*/
/*																			*/
/* Purpose:   cmg graph classes												*/
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
#ifndef __CMG_GRAPH__
#define __CMG_GRAPH__

#include "matrix.h"

/* RCS_ID
   $Header$
 */

const int CMGMAXPARENTS=2;


class CMGList
{
public:
  int GetData() const;
  void Insert(class CMGNode *);
  void Init(CMGList *,CMGList *,int );
  CMGList* GetPred() const;
  CMGList* GetSucc() const;
  class CMGNode* GetFirst() const;
  class CMGNode* GetLast() const;
  void SetData(int);
  void SetPred(CMGList *);
  void SetSucc(CMGList *);
  void SetFirst(class CMGNode *);
  void SetLast(class CMGNode *);
private:
  int data;
  CMGList *succ, *pred;
  class CMGNode *first, *last;
};

inline int CMGList::GetData() const {
  return data;
}
inline CMGList* CMGList::GetPred() const {
  return pred;
}
inline CMGList* CMGList::GetSucc() const {
  return succ;
}
inline class CMGNode* CMGList::GetFirst() const {
  return first;
}
inline class CMGNode* CMGList::GetLast() const {
  return last;
}
inline void CMGList::SetData(int val) {
  data = val;
}
inline void CMGList::SetPred(CMGList *p) {
  pred = p;
}
inline void CMGList::SetSucc(CMGList *s) {
  succ = s;
}
inline void CMGList::SetFirst(class CMGNode *f) {
  first = f;
}
inline void CMGList::SetLast(class CMGNode *l) {
  last = l;
}


class CMGPaList
{
public:
  int GetNp() const;
  int GetPa(int) const;
  int *GetPa() const;
  int GetNewLinks() const;
  double GetNewCG() const;
  double GetApprox() const;
  double GetCoeff(int) const;
  double GetCoefft(int) const;
  double *GetCoeff() const;
  double *GetCoefft() const;
  CMGPaList* GetNext() const;
  void SetNp(int);
  void SetPa(int, int);
  void SetNext(CMGPaList* );
  void SetNewLinks(int);
  void SetNewCG(double);
  void SetCoeff(int, double);
  void SetCoefft(int, double);
  void Init(CMGPaList *nex, int n, int *p, double *c, double *ct, double error);
  void MarkParents(class CMGGrid *grid);
  double TotalWeight();
private:
  int np;
  int pa[CMGMAXPARENTS];
  double coeff[CMGMAXPARENTS];
  double coefft[CMGMAXPARENTS];
  double approx;
  int newlinks;
  double newcg;
  class CMGPaList *next;
};

inline int CMGPaList::GetNp() const {
  return np;
}
inline int CMGPaList::GetPa(int i) const {
  return pa[i];
}
inline int *CMGPaList::GetPa() const {
  return (int *) pa;
}
inline int CMGPaList::GetNewLinks() const {
  return newlinks;
}
inline double CMGPaList::GetApprox() const {
  return approx;
}
inline double CMGPaList::GetNewCG() const {
  return newcg;
}
inline double CMGPaList::GetCoeff(int i) const {
  return coeff[i];
}
inline double CMGPaList::GetCoefft(int i) const {
  return coefft[i];
}
inline double *CMGPaList::GetCoeff() const {
  return (double *) coeff;
}
inline double *CMGPaList::GetCoefft() const {
  return (double *)coefft;
}
inline CMGPaList* CMGPaList::GetNext() const {
  return next;
}
inline void CMGPaList::SetNp(int v) {
  np = v;
}
inline void CMGPaList::SetPa(int i, int p) {
  pa[i] = p;
}
inline void CMGPaList::SetNewLinks(int v) {
  newlinks = v;
}
inline void CMGPaList::SetNewCG(double v) {
  newcg = v;
}
inline void CMGPaList::SetCoeff(int i, double c) {
  coeff[i] = c;
}
inline void CMGPaList::SetCoefft(int i, double c) {
  coefft[i] = c;
}
inline void CMGPaList::SetNext(CMGPaList* ptr) {
  next = ptr;
}



struct CMGNodeBitField
{
  unsigned f0 : 1;
  unsigned f1 : 1;
  unsigned f2 : 1;
  unsigned nt : 2;
  unsigned ns : 11;
  int lid : 16;
};

class CMGNode
{
public:
  int GetData() const;
  int GetId() const;
  int GetNSons() const;
  int GetLocalId() const;
  double GetLocalNormA() const;
  double GetLocalNormB() const;
  CMGNode* GetPred() const;
  CMGNode* GetSucc() const;
  CMGList* GetList() const;
  CMGPaList* GetPaList() const;
  void SetPred(CMGNode *);
  void SetSucc(CMGNode *);
  void SetList(CMGList *);
  void SetData(int val);
  void SetId(int);
  void SetNSons(int);
  void SetLocalId(int);
  void SetLocalNormA(double);
  void SetLocalNormB(double);
  void SetPaList(CMGPaList*);
  int GetFlag() const;
  int GetFlag1() const;
  int GetFlag2() const;
  int IsCGNode() const;
  int IsFGNode() const;
  void SetFlag(int);
  void SetFlag1(int);
  void SetFlag2(int);
  void MarkCGNode();
  void MarkFGNode();
  void Init(int);
  int UpdateNeighborsFG(class CMGGrid *grid);
  int Eliminate(CMGGrid *grid);
  void MarkBestNeighbor(CMGGrid *grid);
  int CountNewLinks(CMGGrid *grid, class CMGGraph *graph);
  void CountNewCG(class CMGGraph *graph);
  int CountNewCG(CMGGraph *graph, int j);
  void ComputeTotalWeight();
private:
  int data;
  int id;
  class CMGNode *pred, *succ;
  class CMGList *list;
  class CMGPaList *palist;
  CMGNodeBitField control;
};

inline int CMGNode::GetData() const {
  return data;
}
inline int CMGNode::GetId() const {
  return id;
}
inline int CMGNode::GetNSons() const {
  return control.ns;
}
inline int CMGNode::GetLocalId() const {
  return control.lid;
}
inline CMGNode* CMGNode::GetPred() const {
  return pred;
}
inline CMGNode* CMGNode::GetSucc() const {
  return succ;
}
inline CMGList* CMGNode::GetList() const {
  return list;
}
inline CMGPaList* CMGNode::GetPaList() const {
  return palist;
}
inline void CMGNode::SetData(int val) {
  data = val;
}
inline void CMGNode::SetId(int i) {
  id = i;
}
inline void CMGNode::SetNSons(int i) {
  control.ns = i;
}
inline void CMGNode::SetLocalId(int i) {
  control.lid = i;
}
inline void CMGNode::SetList(CMGList *l) {
  list = l;
}
inline void CMGNode::SetPred(CMGNode *p) {
  pred = p;
}
inline void CMGNode::SetSucc(CMGNode *s) {
  succ = s;
}
inline void CMGNode::SetPaList(CMGPaList *ptr) {
  palist = ptr;
}
inline int CMGNode::IsCGNode() const {
  return (control.nt == 2);
}
inline int CMGNode::IsFGNode() const {
  return (control.nt == 1);
}
inline int CMGNode::GetFlag() const {
  return control.f0;
}
inline int CMGNode::GetFlag1() const {
  return control.f1;
}
inline int CMGNode::GetFlag2() const {
  return control.f2;
}
inline void CMGNode::MarkCGNode() {
  control.nt = 2;
}
inline void CMGNode::MarkFGNode() {
  control.nt = 1;
}
inline void CMGNode::SetFlag(int f) {
  control.f0 = f;
}
inline void CMGNode::SetFlag1(int f) {
  control.f1 = f;
}
inline void CMGNode::SetFlag2(int f) {
  control.f2 = f;
}


class CMGGraph
{
public:
  int* GetMap() const;
  class CMGDecompRow* GetRow() const;
  int GetNF() const;
  CMGList* GetFreeList() const;
  void SetFreeList(CMGList*);
  int Insert(CMGNode *);
  void Remove(CMGNode *);
  int InsertH(CMGNode *);
  void RemoveH(CMGNode *);
  void Store(CMGNode *);
  int InsertHelplist();
  CMGNode *GetNode() const;
  CMGNode *GetFirstNode();
  CMGNode *GetFirstNodeH();
  CMGNode *GetLastNode();
  int Init(class CMGGrid*);
  int Construct(class CMGGrid *);
  int Construct2(class CMGGrid *);
  int InitList();
  void MarkFGNode(CMGNode *);
  void MarkCGNode(CMGNode *);
  void ClearPaList(CMGPaList *);
  void ClearPaListRev(CMGPaList *&);
  void CorrectPaList(CMGPaList *&palist, double threshold);
  int SavePaList(CMGPaList *&list, int np, int *pa, double *c, double *ct, double error);
  int EliminateNodes(CMGGrid *gridptr);
  int RemainingNodes(CMGGrid *gridptr);
  void UpdateNSons(CMGPaList *newlist, CMGPaList *oldlist, CMGGrid *grid);
  void InitNSons();
  int OrderSpecial(class CMGMatrix *matrix);
  int OrderILUT(CMGMatrix *matrix);
private:
  int n;
  int nf;
  int nc;
  CMGNode *node;
  CMGList *list;
  CMGNode *helplist;
  CMGPaList *freepalist;
  CMGList *freelist;
  int *map;
  class CMGDecompRow *row;
};

inline int* CMGGraph::GetMap() const {
  return map;
}
inline class CMGDecompRow* CMGGraph::GetRow() const {
  return row;
}
inline int CMGGraph::GetNF() const {
  return nf;
}
inline CMGNode *CMGGraph::GetNode() const {
  return node;
}
inline CMGList *CMGGraph::GetFreeList() const {
  return freelist;
}
inline void CMGGraph::SetFreeList(CMGList *pl) {
  freelist = pl;
}


#endif
