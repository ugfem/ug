/****************************************************************************/
/*																			*/
/* File:      famg_approxSparse.C											*/
/*																			*/
/* Purpose:   famg parents selection process                                */
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
/* Remarks:																*/
/*																			*/
/****************************************************************************/


#include <iostream.h>
#include <strstream.h>
#include <math.h>

#include "famg_algebra.h"
#include "famg_misc.h"
#include "famg_heap.h"
#include "famg_grid.h"
#include "famg_graph.h"
#include "famg_system.h"
#include "famg_sparse.h"

#ifdef FAMG_SPARSE_BLOCK

struct FAMGMatrixLocal
{
    int nn;
    int ns;
    int *index;
    int *start;
    int *end;
    double *mat;
    double *matD;
    double *matT;
    double *vec;
    double *vecT;
    int *lmap;
};

static int time1 = 0;
static int time2 = 0;
static int time3 = 0;
static int time4 = 0;

    // test
static double help[2];

int FAMGGrid::ConstructLocalMatrixA(const FAMGVectorEntry &veci, struct FAMGMatrixLocal *localmatrix, double *&tv1, double *&tv2, double &normr, double &norml)
{  
    const FAMGGraph &graph = *GetGraph();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
	const FAMGMatrixAlg &A = *GetTmpMatrix();
	const FAMGMatrixAlg &D = *GetDiagMatrix();
    FAMGMatrixIter mi_iter(A,veci); 
	FAMGVectorEntry vecj, veck;
	FAMGMatrixEntry matij, matjk, matkk;
    double norm1, norm2;
    int i, j, k, lid, ind, ix;
    int nn, ns, nm, nm1, nm2;
    int  *start, *end, *index, *lmap;
    double *mat, *matD, *matT, *vec, *vecT;

    const double omega = FAMGGetParameter()->Getomegar();

    const FAMGSparseBlock *Dsb = D.GetDiagSparseBlockPtr();
    const FAMGSparseBlock *Asb = A.GetSparseBlockPtr();
    const FAMGSparseBlock *AsbT = A.GetSparseBlockAdjPtr();
    const FAMGSparseBlock *sb1 = graph.Get_sb1Ptr();
    const FAMGSparseBlock *sb1T = graph.Get_sb1TPtr();

    const FAMGSparseVector *tvAsv = tvA.GetSparseVectorPtr();
    const FAMGSparseVector *tvBsv = tvB.GetSparseVectorPtr();
    const FAMGSparseVector *stv1 = graph.Get_stvPtr();
    const FAMGSparseVector *stv2 = graph.Get_stvTPtr();


    short sbsize = sb1->Get_maxoffset()+1;
    short sbTsize = sb1T->Get_maxoffset()+1;
    short vsize = stv1->Get_maxcomp()+1; 
    short vTsize = stv2->Get_maxcomp()+1; 


    double *t1, *t2;
    double *miiinv, *mjjinv, *mkkinv, *mij, *mji, *mjk, *mkj;
    mij = new double[sbsize];
    mji = new double[sbTsize];
    mjk = new double[sbsize];
    mkj = new double[sbTsize];

    // set local IDs
    graph.GetNode(veci.GetIndex())->SetLocalId(0);
    mi_iter(matij);   
    nn = 1; 
    norm1 = norm2 = 0.0;
	while(mi_iter(matij))
    {
        if(matij.is_strong())
        {
            graph.GetNode(matij.dest().GetIndex())->SetLocalId(nn);
            nn++;
        }
    }
    ns = nn;

    nm1 = 0; nm2 = 0;
    mi_iter.reset();
    mi_iter(matij);   
	while(mi_iter(matij))
    {
        if(matij.is_strong())
        {
            vecj = matij.dest();
            FAMGMatrixIter mj_iter(A,vecj);
            mj_iter(matjk);
            while(mj_iter(matjk))
            {
                if(matjk.is_strong())
                {
                    k = matjk.dest().GetIndex();
                    lid = graph.GetNode(k)->GetLocalId();
                    if(lid < 0)
                    {
                        graph.GetNode(k)->SetLocalId(nn);
                        lid = nn;
                        nn++;
                    }
                    // number of entries per row is stored in NSons, 
                    // not very nice though.
                    graph.GetNode(k)->SetNSons(graph.GetNode(k)->GetNSons()+1);
                    if (lid < ns) nm1++;
                    else nm2++;
                    
                }
            }
        }
    }
    nm = ns + nm1 + 2*nm2 + nn -1;


    // allocate start, end, index, mat, matT
    start = (int *) FAMGGetMem((nn+1)*sizeof(int), FAMG_FROM_TOP); 
    end = (int *) FAMGGetMem((nn+1)*sizeof(int), FAMG_FROM_TOP); 
    index = (int *) FAMGGetMem((nm)*sizeof(int), FAMG_FROM_TOP); 
    mat = (double *) FAMGGetMem((nm*sbsize)*sizeof(double), FAMG_FROM_TOP); 
    matD = (double *) FAMGGetMem((nm*sbsize)*sizeof(double), FAMG_FROM_TOP);
    matT = (double *) FAMGGetMem((nm*sbTsize)*sizeof(double), FAMG_FROM_TOP); 
    vec = (double *)  FAMGGetMem((ns*sbsize)*sizeof(double), FAMG_FROM_TOP); 
    vecT = (double *)  FAMGGetMem((ns*sbTsize)*sizeof(double), FAMG_FROM_TOP); 
    tv1 = (double *)  FAMGGetMem((ns*vsize)*sizeof(double), FAMG_FROM_TOP); 
    tv2 = (double *)  FAMGGetMem((ns*vTsize)*sizeof(double), FAMG_FROM_TOP); 
    lmap = (int *)  FAMGGetMem((nn)*sizeof(int), FAMG_FROM_TOP); 


    // copy matrix
    // this function uses a couple of strange tricks
    // i guess it works :-)
    i = 1; 
    lmap[0] = veci.GetIndex();
    start[0] = 0; end[0] = ns; start[1] = ns; 
    start[ns] = ns+nm1+nm2+ns-1; end[ns] = start[ns] + 1;
    index[0] = 0;

    mi_iter.reset(); mi_iter(matij);
    miiinv = D.GetValuePtr(matij);



    // test
    help[0] = help[1] = 1.0;
    LR_Solve(2,miiinv+Dsb->Get_offset(0),help,help); 










    t1 = tvA.GetValuePtr(veci);
    SparseBlockDiagApprox(sb1,Dsb,tvAsv,matD+0*sbsize,miiinv,1.0-omega,t1);
    SparseBlockMSetDiag(sb1,mat+0*sbsize,1.0-omega);
    SparseBlockMSetDiag(sb1T,matT+0*sbTsize,1.0-omega);
    SparseBlockVSet(stv1,tv1+0*vsize,0.0);
    SparseBlockVSet(stv2,tv2+0*vTsize,0.0);

    while(mi_iter(matij))
    {
        ind = start[i];
        vecj = matij.dest();
        t1 = tvA.GetValuePtr(vecj); t2 = tvB.GetValuePtr(vecj);

        SparseBlockDiagApprox(sb1,Dsb,Asb,tvAsv,mij,miiinv, A.GetValuePtr(matij),t1);
        SparseBlockDiagApproxT(sb1T,Dsb,AsbT,tvBsv,mji,miiinv, A.GetAdjValuePtr(matij),t2);
 
        // SparseBlockMMProduct(sb1,Dsb,Asb,mij,miiinv,A.GetValuePtr(matij)); 
        // SparseBlockMMProduct(sb1T,Dsb,AsbT,mji,miiinv,A.GetAdjValuePtr(matij));

        SparseBlockMVAddProduct(stv1,sb1,tvAsv,tv1+0*vsize,mij,t1,-1.0);
        SparseBlockMVAddProduct(stv2,sb1T,tvBsv,tv2+0*vTsize,mji,t2,-1.0);
        norm1 += SparseBlockMNorm(sb1,mij); 
        norm2 += SparseBlockMNorm(sb1T,mji); 
        if(matij.is_strong())
        {
            lmap[i] = vecj.GetIndex();
            SparseBlockMCopy(sb1,sb1,vec+i*sbsize,mij,-1.0);
            SparseBlockMCopy(sb1T,sb1T,vecT+i*sbTsize,mji,-1.0);
            SparseBlockVCopy(stv1,tvAsv,tv1+i*vsize,t1,1.0);
            SparseBlockVCopy(stv2,tvBsv,tv2+i*vTsize,t2,1.0);
            FAMGMatrixIter mj_iter(A,vecj);  
            mj_iter(matjk);
            mjjinv = D.GetValuePtr(matjk);
            SparseBlockMCopy(sb1,sb1,mat+i*sbsize,mij,-omega);
            SparseBlockMCopy(sb1T,sb1T,matT+i*sbTsize,mji,-omega);
            // recompute mij, can be accelerated
            SparseBlockDiagApprox(sb1,Dsb,Asb,tvAsv,mij,miiinv, A.GetValuePtr(matij),mjjinv,t1);
            SparseBlockMCopy(sb1,sb1,matD+i*sbsize,mij,-omega);
            index[i] = i;
            SparseBlockMSetDiag(sb1,mat+ind*sbsize,1.0-omega);
            SparseBlockDiagApprox(sb1,Dsb,tvAsv,matD+ind*sbsize,mjjinv,1.0-omega,t1);
            SparseBlockMSetDiag(sb1T,matT+ind*sbTsize,1.0-omega);
            index[ind] = i;
            ind++;
            while(mj_iter(matjk))
            {
                if(matjk.is_strong())
                {
                    veck = matjk.dest(); 
                    k = veck.GetIndex(); 
                    lid = graph.GetNode(k)->GetLocalId();
                    // todo: not sure about tv
                    t1 = tvA.GetValuePtr(veck); t2 = tvB.GetValuePtr(veck);
                    FAMGMatrixIter mk_iter(A,veck); 
                    mk_iter(matkk);
                    mkkinv = D.GetValuePtr(matkk);
                    SparseBlockDiagApprox(sb1,Dsb,Asb,tvAsv,mjk,mjjinv, A.GetValuePtr(matjk),t1);
                    SparseBlockDiagApproxT(sb1T,Dsb,AsbT,tvBsv,mkj,mjjinv, A.GetAdjValuePtr(matjk),t2);
                    SparseBlockMCopy(sb1,sb1,mat+ind*sbsize,mjk,-omega);
                    SparseBlockMCopy(sb1T,sb1T,matT+ind*sbTsize,mkj,-omega);
                    // recompute mjk, can be accelerated
                    SparseBlockDiagApprox(sb1,Dsb,Asb,tvAsv,mjk,mjjinv, A.GetValuePtr(matjk),mkkinv,t1);
                    SparseBlockMCopy(sb1,sb1,matD+ind*sbsize,mjk,-omega);
                    index[ind] = lid;
                    ind++;
                    // construct rows for the non-neighbor nodes
                    if(lid >= ns)
                    {
                        ix = end[lid];
                        // todo: not sure about tv
                        // t1 = tvA.GetValuePtr(vecj); t2 = tvB.GetValuePtr(vecj);
                        SparseBlockDiagApprox(sb1,Dsb,Asb,tvAsv,mkj,mkkinv, A.GetAdjValuePtr(matjk),t1);
                        SparseBlockDiagApproxT(sb1T,Dsb,AsbT,tvBsv,mjk,mkkinv, A.GetValuePtr(matjk),t2);
                        SparseBlockMCopy(sb1,sb1,mat+ix*sbsize,mkj,-omega);
                        SparseBlockMCopy(sb1T,sb1T,matT+ix*sbTsize,mjk,-omega);
                        // recompute mkj, can be accelerated
                        SparseBlockDiagApprox(sb1,Dsb,Asb,tvAsv,mkj,mkkinv, A.GetAdjValuePtr(matjk),mjjinv,t1);
                        SparseBlockMCopy(sb1,sb1,matD+ix*sbsize,mkj,-omega);
                        // GetAdjValuePtr yields the connection k->j. 
                        // For mat we need this block, for matT we need the transposed. 
                        //  Thus, the two lines above should be correct.
                        index[ix] = i;
                        end[lid] = ix+1;
                        if(start[lid] == ix-1)
                        {
                            // prepare start/end for the next non-neighbor node
                            // The ordering must be the same as in the previous loop
                            if(lid  < nn)
                            {
                                lmap[lid] = k;
                                index[start[lid]] = lid;
                                // t1 = tvA.GetValuePtr(veck);
                                SparseBlockDiagApprox(sb1,Dsb,tvAsv,matD+start[lid]*sbsize,mkkinv,1.0-omega,t1);
                                SparseBlockMSetDiag(sb1,mat+start[lid]*sbsize,1.0-omega);
                                SparseBlockMSetDiag(sb1T,matT+start[lid]*sbTsize,1.0-omega);
                                start[lid+1] = ix + graph.GetNode(k)->GetNSons();
                                end[lid+1] = start[lid+1] + 1;
                            }
                        }
                    }

                }
            }
            end[i] = ind;
            i++;
            start[i] = ind;
        }
    }

    normr = norm1; norml = norm2;

             
    // reset flags and counter
    mi_iter.reset();
    graph.GetNode(veci.GetIndex())->SetLocalId(-1);
    mi_iter(matij);
    while(mi_iter(matij))
    {
        if(matij.is_strong())
        {
            vecj = matij.dest();
            j = vecj.GetIndex();
            graph.GetNode(j)->SetLocalId(-1);
            FAMGMatrixIter mj_iter(A,vecj);  
            mj_iter(matjk);
            while(mj_iter(matjk))
            {
                if(matjk.is_strong())
                {
                    k = matjk.dest().GetIndex();
                    graph.GetNode(k)->SetLocalId(-1);
                    graph.GetNode(k)->SetNSons(0);
                }
            }
        }
    }



    localmatrix->index = index;
    localmatrix->start = start;
    localmatrix->end = end;
    localmatrix->mat = mat;
    localmatrix->matD = matD;
    localmatrix->matT = matT;
    localmatrix->vec = vec;
    localmatrix->vecT = vecT;
    localmatrix->ns = ns;
    localmatrix->nn = nn;
    localmatrix->lmap = lmap;

    // global variables todo: save them somewhere
    

    delete mij; delete mji; delete mjk; delete mkj;


    return 0;

}



int FAMGGrid::ConstructLocalMatrix(const FAMGVectorEntry &veci, struct FAMGMatrixLocal *localmatrix, double *&tv1, double *&tv2, double &normr, double &norml)
{  
    const FAMGGraph &graph = *GetGraph();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
	const FAMGMatrixAlg &A = *GetTmpMatrix();
	const FAMGMatrixAlg &D = *GetDiagMatrix();
    FAMGMatrixIter mi_iter(A,veci); 
	FAMGVectorEntry vecj, veck;
	FAMGMatrixEntry matij, matjk, matkk;
    double norm1, norm2;
    int i, j, k, lid, ind, ix;
    int nn, ns, nm, nm1, nm2;
    int  *start, *end, *index, *lmap;
    double *mat, *matT, *vec, *vecT;

    const double omega = FAMGGetParameter()->Getomegar();

    const FAMGSparseBlock *Dsb = D.GetDiagSparseBlockPtr();
    const FAMGSparseBlock *Asb = A.GetSparseBlockPtr();
    const FAMGSparseBlock *AsbT = A.GetSparseBlockAdjPtr();
    const FAMGSparseBlock *sb1 = graph.Get_sb1Ptr();
    const FAMGSparseBlock *sb1T = graph.Get_sb1TPtr();

    const FAMGSparseVector *tvAsv = tvA.GetSparseVectorPtr();
    const FAMGSparseVector *tvBsv = tvB.GetSparseVectorPtr();
    const FAMGSparseVector *stv1 = graph.Get_stvPtr();
    const FAMGSparseVector *stv2 = graph.Get_stvTPtr();


    short sbsize = sb1->Get_maxoffset()+1;
    short sbTsize = sb1T->Get_maxoffset()+1;
    short vsize = stv1->Get_maxcomp()+1; 
    short vTsize = stv2->Get_maxcomp()+1; 


    double *t1, *t2;
    double *miiinv, *mjjinv, *mkkinv, *mij, *mji, *mjk, *mkj;
    mij = new double[sbsize];
    mji = new double[sbTsize];
    mjk = new double[sbsize];
    mkj = new double[sbTsize];

    // set local IDs
    graph.GetNode(veci.GetIndex())->SetLocalId(0);
    mi_iter(matij);   
    nn = 1; 
    norm1 = norm2 = 0.0;
	while(mi_iter(matij))
    {
        if(matij.is_strong())
        {
            graph.GetNode(matij.dest().GetIndex())->SetLocalId(nn);
            nn++;
        }
    }
    ns = nn;

    nm1 = 0; nm2 = 0;
    mi_iter.reset();
    mi_iter(matij);   
	while(mi_iter(matij))
    {
        if(matij.is_strong())
        {
            vecj = matij.dest();
            FAMGMatrixIter mj_iter(A,vecj);
            mj_iter(matjk);
            while(mj_iter(matjk))
            {
                if(matjk.is_strong())
                {
                    k = matjk.dest().GetIndex();
                    lid = graph.GetNode(k)->GetLocalId();
                    if(lid < 0)
                    {
                        graph.GetNode(k)->SetLocalId(nn);
                        lid = nn;
                        nn++;
                    }
                    // number of entries per row is stored in NSons, 
                    // not very nice though.
                    graph.GetNode(k)->SetNSons(graph.GetNode(k)->GetNSons()+1);
                    if (lid < ns) nm1++;
                    else nm2++;
                    
                }
            }
        }
    }
    nm = ns + nm1 + 2*nm2 + nn -1;


    // allocate start, end, index, mat, matT
    start = (int *) FAMGGetMem((nn+1)*sizeof(int), FAMG_FROM_TOP); 
    end = (int *) FAMGGetMem((nn+1)*sizeof(int), FAMG_FROM_TOP); 
    index = (int *) FAMGGetMem((nm)*sizeof(int), FAMG_FROM_TOP); 
    mat = (double *) FAMGGetMem((nm*sbsize)*sizeof(double), FAMG_FROM_TOP); 
    matT = (double *) FAMGGetMem((nm*sbTsize)*sizeof(double), FAMG_FROM_TOP); 
    vec = (double *)  FAMGGetMem((ns*sbsize)*sizeof(double), FAMG_FROM_TOP); 
    vecT = (double *)  FAMGGetMem((ns*sbTsize)*sizeof(double), FAMG_FROM_TOP); 
    tv1 = (double *)  FAMGGetMem((ns*vsize)*sizeof(double), FAMG_FROM_TOP); 
    tv2 = (double *)  FAMGGetMem((ns*vTsize)*sizeof(double), FAMG_FROM_TOP); 
    lmap = (int *)  FAMGGetMem((nn)*sizeof(int), FAMG_FROM_TOP); 


    // copy matrix
    // this function uses a couple of strange tricks
    // i guess it works :-)
    i = 1; 
    lmap[0] = veci.GetIndex();
    start[0] = 0; end[0] = ns; start[1] = ns; 
    start[ns] = ns+nm1+nm2+ns-1; end[ns] = start[ns] + 1;
    index[0] = 0;

    mi_iter.reset(); mi_iter(matij);
    miiinv = D.GetValuePtr(matij);

    SparseBlockMSetDiag(sb1,mat+0*sbsize,1.0-omega);
    SparseBlockMSetDiag(sb1T,matT+0*sbTsize,1.0-omega);
    SparseBlockVSet(stv1,tv1+0*vsize,0.0);
    SparseBlockVSet(stv2,tv2+0*vTsize,0.0);

    while(mi_iter(matij))
    {
        ind = start[i];
        vecj = matij.dest();
        t1 = tvA.GetValuePtr(vecj); t2 = tvB.GetValuePtr(vecj);

        SparseBlockDiagApprox(sb1,Dsb,Asb,tvAsv,mij,miiinv, A.GetValuePtr(matij),t1);
        SparseBlockDiagApproxT(sb1T,Dsb,AsbT,tvBsv,mji,miiinv, A.GetAdjValuePtr(matij),t2);
 
        // SparseBlockMMProduct(sb1,Dsb,Asb,mij,miiinv,A.GetValuePtr(matij)); 
        // SparseBlockMMProduct(sb1T,Dsb,AsbT,mji,miiinv,A.GetAdjValuePtr(matij));

        SparseBlockMVAddProduct(stv1,sb1,tvAsv,tv1+0*vsize,mij,t1,-1.0);
        SparseBlockMVAddProduct(stv2,sb1T,tvBsv,tv2+0*vTsize,mji,t2,-1.0);
        norm1 += SparseBlockMNorm(sb1,mij); 
        norm2 += SparseBlockMNorm(sb1T,mji); 
        if(matij.is_strong())
        {
            lmap[i] = vecj.GetIndex();
            SparseBlockMCopy(sb1,sb1,vec+i*sbsize,mij,-1.0);
            SparseBlockMCopy(sb1T,sb1T,vecT+i*sbTsize,mji,-1.0);
            SparseBlockVCopy(stv1,tvAsv,tv1+i*vsize,t1,1.0);
            SparseBlockVCopy(stv2,tvBsv,tv2+i*vTsize,t2,1.0);
            FAMGMatrixIter mj_iter(A,vecj);  
            mj_iter(matjk);
            mjjinv = D.GetValuePtr(matjk);
            SparseBlockMCopy(sb1,sb1,mat+i*sbsize,mij,-omega);
            SparseBlockMCopy(sb1T,sb1T,matT+i*sbTsize,mji,-omega);
            index[i] = i;
            SparseBlockMSetDiag(sb1,mat+ind*sbsize,1.0-omega);
            SparseBlockMSetDiag(sb1T,matT+ind*sbTsize,1.0-omega);
            index[ind] = i;
            ind++;
            while(mj_iter(matjk))
            {
                if(matjk.is_strong())
                {
                    veck = matjk.dest(); 
                    k = veck.GetIndex(); 
                    lid = graph.GetNode(k)->GetLocalId();
                    // todo: not sure about tv
                    t1 = tvA.GetValuePtr(veck); t2 = tvB.GetValuePtr(veck);
                    FAMGMatrixIter mk_iter(A,veck); 
                    mk_iter(matkk);
                    mkkinv = D.GetValuePtr(matkk);
                    SparseBlockDiagApprox(sb1,Dsb,Asb,tvAsv,mjk,mjjinv, A.GetValuePtr(matjk),t1);
                    SparseBlockDiagApproxT(sb1T,Dsb,AsbT,tvBsv,mkj,mjjinv, A.GetAdjValuePtr(matjk),t2);
                    // SparseBlockMMProduct(sb1,Dsb,Asb,mjk,mjjinv,A.GetValuePtr(matjk)); 
                    // SparseBlockMMProduct(sb1T,Dsb,AsbT,mkj,mjjinv,A.GetAdjValuePtr(matjk));
                    SparseBlockMCopy(sb1,sb1,mat+ind*sbsize,mjk,-omega);
                    SparseBlockMCopy(sb1T,sb1T,matT+ind*sbTsize,mkj,-omega);
                    index[ind] = lid;
                    ind++;
                    // construct rows for the non-neighbor nodes
                    if(lid >= ns)
                    {
                        ix = end[lid];
                        // todo: not sure about tv
                        t1 = tvA.GetValuePtr(veck); t2 = tvB.GetValuePtr(veck);
                        SparseBlockDiagApprox(sb1,Dsb,Asb,tvAsv,mkj,mkkinv, A.GetAdjValuePtr(matjk),t1);
                        SparseBlockDiagApproxT(sb1T,Dsb,AsbT,tvBsv,mjk,mkkinv, A.GetValuePtr(matjk),t2);
                        // SparseBlockMMProduct(sb1,Dsb,Asb,mkj,mkkinv,A.GetAdjValuePtr(matjk)); 
                        //  SparseBlockMMProduct(sb1T,Dsb,AsbT,mjk,mkkinv,A.GetValuePtr(matjk));
                        SparseBlockMCopy(sb1,sb1,mat+ix*sbsize,mkj,-omega);
                        SparseBlockMCopy(sb1T,sb1T,matT+ix*sbTsize,mjk,-omega);
                        // GetAdjValuePtr yields the connection k->j. 
                        // For mat we need this block, for matT we need the transposed. 
                        //  Thus, the two lines above should be correct.
                        index[ix] = i;
                        end[lid] = ix+1;
                        if(start[lid] == ix-1)
                        {
                            // prepare start/end for the next non-neighbor node
                            // The ordering must be the same as in the previous loop
                            if(lid  < nn)
                            {
                                lmap[lid] = k;
                                index[start[lid]] = lid;
                                SparseBlockMSetDiag(sb1,mat+start[lid]*sbsize,1.0-omega);
                                SparseBlockMSetDiag(sb1T,matT+start[lid]*sbTsize,1.0-omega);
                                start[lid+1] = ix + graph.GetNode(k)->GetNSons();
                                end[lid+1] = start[lid+1] + 1;
                            }
                        }
                    }

                }
            }
            end[i] = ind;
            i++;
            start[i] = ind;
        }
    }

    normr = norm1; norml = norm2;

             
    // reset flags and counter
    mi_iter.reset();
    graph.GetNode(veci.GetIndex())->SetLocalId(-1);
    mi_iter(matij);
    while(mi_iter(matij))
    {
        if(matij.is_strong())
        {
            vecj = matij.dest();
            j = vecj.GetIndex();
            graph.GetNode(j)->SetLocalId(-1);
            FAMGMatrixIter mj_iter(A,vecj);  
            mj_iter(matjk);
            while(mj_iter(matjk))
            {
                if(matjk.is_strong())
                {
                    k = matjk.dest().GetIndex();
                    graph.GetNode(k)->SetLocalId(-1);
                    graph.GetNode(k)->SetNSons(0);
                }
            }
        }
    }



    localmatrix->index = index;
    localmatrix->start = start;
    localmatrix->end = end;
    localmatrix->mat = mat;
    localmatrix->matT = matT;
    localmatrix->vec = vec;
    localmatrix->vecT = vecT;
    localmatrix->ns = ns;
    localmatrix->nn = nn;
    localmatrix->lmap = lmap;

    // global variables todo: save them somewhere
    

    delete mij; delete mji; delete mjk; delete mkj;


    return 0;

}

static int LocalMatMult(const FAMGGraph &graph, struct FAMGMatrixLocal *new_localmatrix, struct FAMGMatrixLocal *localmatrix)
{
    double *v1ik, *v2ik, *v1kj, *v2kj, *factor, *factorT;
    int i, j, k, ind, new_ind, ix;
    int new_nn, new_ns, nn, ns, nm;
    int *new_start, *new_end, *new_index, *new_lmap, *start, *end, *index, *lmap;
    double *new_mat, *new_matT, *mat, *matT, *vec, *vecT, *new_mat2, *new_matT2;

    index = localmatrix->index;
    start = localmatrix->start;
    end = localmatrix->end;
    mat = localmatrix->mat;
    matT = localmatrix->matT;
    vec = localmatrix->vec;
    vecT = localmatrix->vecT;
    ns = localmatrix->ns;
    nn = localmatrix->nn;
    lmap = localmatrix->lmap;

    const FAMGSparseBlock *sb1 = graph.Get_sb1Ptr();
    const FAMGSparseBlock *sb1T = graph.Get_sb1TPtr();
    const FAMGSparseBlock *sb2 = graph.Get_sb2Ptr();
    const FAMGSparseBlock *sb2T = graph.Get_sb2TPtr();
    const FAMGSparseBlock *sb3 = graph.Get_sb3Ptr();
    const FAMGSparseBlock *sb3T = graph.Get_sb3TPtr();


    short sb1size = sb1->Get_maxoffset()+1;
    short sb1Tsize = sb1T->Get_maxoffset()+1;
    short sb2size = sb2->Get_maxoffset()+1;
    short sb2Tsize = sb2T->Get_maxoffset()+1;
    short sb3size = sb3->Get_maxoffset()+1;
    short sb3Tsize = sb3T->Get_maxoffset()+1;

    nm = ns*nn; // just to be sure

    // allocate start, end, index, mat, matT
    new_start = (int *) FAMGGetMem((ns+1)*sizeof(int), FAMG_FROM_TOP); 
    new_end = (int *) FAMGGetMem((ns+1)*sizeof(int), FAMG_FROM_TOP); 
    new_index = (int *) FAMGGetMem((nm)*sizeof(int), FAMG_FROM_TOP); 
    new_mat = (double *) FAMGGetMem((nn*sb3size + (nm-nn)*sb2size)*sizeof(double), FAMG_FROM_TOP); 
    new_matT = (double *) FAMGGetMem((nn*sb3Tsize + (nm-nn)*sb2Tsize)*sizeof(double), FAMG_FROM_TOP); 
    new_lmap = lmap;
    new_nn = nn;
    new_ns = ns;
    
    new_mat2 = new_mat + nn*(sb3size -sb2size); // help pointer
    new_matT2 = new_matT + nn*(sb3Tsize -sb2Tsize); // help pointer

    
    memset((void *)new_mat,0,(nn*sb3size + (nm-nn)*sb2size)*sizeof(double));
    memset((void *)new_matT,0,(nn*sb3Tsize + (nm-nn)*sb2Tsize)*sizeof(double));
 
    FAMGMarkHeap(FAMG_FROM_TOP);
    double *helpvec1 = (double *) FAMGGetMem((nn*sb2size)*sizeof(double), FAMG_FROM_TOP);
    double *helpvec2 = (double *) FAMGGetMem((nn*sb2Tsize)*sizeof(double), FAMG_FROM_TOP);
    short *compute = (short *) FAMGGetMem((nn)*sizeof(short), FAMG_FROM_TOP);

    new_start[0] = 0; new_end[0] = nn; new_start[1] = nn; // never change this !
    for(i = 1; i < ns; i++)
    {
        memset((void *)helpvec1,0,(nn*sb2size)*sizeof(double)); 
        memset((void *)helpvec2,0,(nn*sb2Tsize)*sizeof(double));
        memset((void *)compute,0,(nn)*sizeof(short));

        for(ind = start[i]; ind < end[i]; ind++)
        {
            k = index[ind];
            v1ik = mat+ind*sb1size;
            v2ik = matT+ind*sb1Tsize;

            for(ix = start[k]; ix < end[k]; ix++)
            {
                j = index[ix];
                v1kj = mat+ix*sb1size;
                v2kj = matT+ix*sb1Tsize;

                SparseBlockMMAddProduct(sb2,sb1,sb1,helpvec1+j*sb2size,v1ik,v1kj,1.0); 
                SparseBlockMMAddProduct(sb2T,sb1T,sb1T,helpvec2+j*sb2Tsize,v2ik,v2kj,1.0); 
                compute[j] = 1;
            }
        }
        new_ind = new_start[i];
        for(j = 0; j < nn; j++)
        {
            if(compute[j])
            {
                SparseBlockMCopy(sb2,sb2,new_mat2+new_ind*sb2size,helpvec1+j*sb2size,1.0);
                SparseBlockMCopy(sb2T,sb2T,new_matT2+new_ind*sb2Tsize,helpvec2+j*sb2Tsize,1.0);
                new_index[new_ind] = j;
                new_ind++;
            }
        }
        new_start[i+1] = new_end[i] = new_ind;
    }
                
 

    new_ind = new_start[0];
    for(ind = new_start[0], j = 0; ind < new_end[0]; ind++, j++)
    {
        new_index[ind] = j;
    }
    for(i = 1; i < ns; i++)
    {
        factor = vec+i*sb1size;
        factorT = vecT+i*sb1Tsize;
        for(ind = new_start[i]; ind < new_end[i]; ind++)
        {
            j = new_index[ind];
            SparseBlockMMAddProduct(sb3,sb1,sb2,new_mat+(new_ind+j)*sb3size,factor,new_mat+ind*sb2size,1.0); 
            SparseBlockMMAddProduct(sb3T,sb1T,sb2T,new_matT+(new_ind+j)*sb3size,factorT,new_matT+ind*sb2Tsize,1.0); 
        }
    }

    new_localmatrix->index = new_index;
    new_localmatrix->start = new_start;
    new_localmatrix->end = new_end;
    new_localmatrix->mat = new_mat;
    new_localmatrix->matT = new_matT;
    new_localmatrix->ns = new_ns;
    new_localmatrix->nn = new_nn;
    new_localmatrix->lmap = new_lmap;

    FAMGReleaseHeap(FAMG_FROM_TOP);

    return 0;
                
}

static int LocalMatMultA(const FAMGGraph &graph, struct FAMGMatrixLocal *new_localmatrix, struct FAMGMatrixLocal *localmatrix)
{
    double *v1ik, *v2ik, *v1kj, *v2kj, *factor, *factorT;
    int i, j, k, ind, new_ind, ix;
    int new_nn, new_ns, nn, ns, nm;
    int *new_start, *new_end, *new_index, *new_lmap, *start, *end, *index, *lmap;
    double *new_mat, *new_matT, *mat, *matD, *matT, *vec, *vecT, *new_mat2, *new_matT2;

    index = localmatrix->index;
    start = localmatrix->start;
    end = localmatrix->end;
    mat = localmatrix->mat;
    matD = localmatrix->matD;
    matT = localmatrix->matT;
    vec = localmatrix->vec;
    vecT = localmatrix->vecT;
    ns = localmatrix->ns;
    nn = localmatrix->nn;
    lmap = localmatrix->lmap;

    const FAMGSparseBlock *sb1 = graph.Get_sb1Ptr();
    const FAMGSparseBlock *sb1T = graph.Get_sb1TPtr();
    const FAMGSparseBlock *sb2 = graph.Get_sb2Ptr();
    const FAMGSparseBlock *sb2T = graph.Get_sb2TPtr();
    const FAMGSparseBlock *sb3 = graph.Get_sb3Ptr();
    const FAMGSparseBlock *sb3T = graph.Get_sb3TPtr();


    short sb1size = sb1->Get_maxoffset()+1;
    short sb1Tsize = sb1T->Get_maxoffset()+1;
    short sb2size = sb2->Get_maxoffset()+1;
    short sb2Tsize = sb2T->Get_maxoffset()+1;
    short sb3size = sb3->Get_maxoffset()+1;
    short sb3Tsize = sb3T->Get_maxoffset()+1;

    nm = ns*nn; // just to be sure

    // allocate start, end, index, mat, matT
    new_start = (int *) FAMGGetMem((ns+1)*sizeof(int), FAMG_FROM_TOP); 
    new_end = (int *) FAMGGetMem((ns+1)*sizeof(int), FAMG_FROM_TOP); 
    new_index = (int *) FAMGGetMem((nm)*sizeof(int), FAMG_FROM_TOP); 
    new_mat = (double *) FAMGGetMem((nn*sb3size + (nm-nn)*sb2size)*sizeof(double), FAMG_FROM_TOP); 
    new_matT = (double *) FAMGGetMem((nn*sb3Tsize + (nm-nn)*sb2Tsize)*sizeof(double), FAMG_FROM_TOP); 
    new_lmap = lmap;
    new_nn = nn;
    new_ns = ns;
    
    new_mat2 = new_mat + nn*(sb3size -sb2size); // help pointer
    new_matT2 = new_matT + nn*(sb3Tsize -sb2Tsize); // help pointer

    
    memset((void *)new_mat,0,(nn*sb3size + (nm-nn)*sb2size)*sizeof(double));
    memset((void *)new_matT,0,(nn*sb3Tsize + (nm-nn)*sb2Tsize)*sizeof(double));
 
    FAMGMarkHeap(FAMG_FROM_TOP);
    double *helpvec1 = (double *) FAMGGetMem((nn*sb2size)*sizeof(double), FAMG_FROM_TOP);
    double *helpvec2 = (double *) FAMGGetMem((nn*sb2Tsize)*sizeof(double), FAMG_FROM_TOP);
    short *compute = (short *) FAMGGetMem((nn)*sizeof(short), FAMG_FROM_TOP);

    new_start[0] = 0; new_end[0] = nn; new_start[1] = nn; // never change this !
    for(i = 1; i < ns; i++)
    {
        memset((void *)helpvec1,0,(nn*sb2size)*sizeof(double)); 
        memset((void *)helpvec2,0,(nn*sb2Tsize)*sizeof(double));
        memset((void *)compute,0,(nn)*sizeof(short));

        for(ind = start[i]; ind < end[i]; ind++)
        {
            k = index[ind];
            v1ik = mat+ind*sb1size;
            v2ik = matT+ind*sb1Tsize;

            for(ix = start[k]; ix < end[k]; ix++)
            {
                j = index[ix];
                v1kj = matD+ix*sb1size;
                v2kj = matT+ix*sb1Tsize;
                // v2kj = matD+ix*sb1Tsize;

                SparseBlockMMAddProduct(sb2,sb1,sb1,helpvec1+j*sb2size,v1ik,v1kj,1.0); 
                SparseBlockMMAddProduct(sb2T,sb1T,sb1T,helpvec2+j*sb2Tsize,v2ik,v2kj,1.0); 
                compute[j] = 1;
            }
        }
        new_ind = new_start[i];
        for(j = 0; j < nn; j++)
        {
            if(compute[j])
            {
                SparseBlockMCopy(sb2,sb2,new_mat2+new_ind*sb2size,helpvec1+j*sb2size,1.0);
                SparseBlockMCopy(sb2T,sb2T,new_matT2+new_ind*sb2Tsize,helpvec2+j*sb2Tsize,1.0);
                new_index[new_ind] = j;
                new_ind++;
            }
        }
        new_start[i+1] = new_end[i] = new_ind;
    }
                
 

    new_ind = new_start[0];
    for(ind = new_start[0], j = 0; ind < new_end[0]; ind++, j++)
    {
        new_index[ind] = j;
    }
    for(i = 1; i < ns; i++)
    {
        factor = vec+i*sb1size;
        factorT = vecT+i*sb1Tsize;
        for(ind = new_start[i]; ind < new_end[i]; ind++)
        {
            j = new_index[ind];
            SparseBlockMMAddProduct(sb3,sb1,sb2,new_mat+(new_ind+j)*sb3size,factor,new_mat+ind*sb2size,1.0); 
            SparseBlockMMAddProduct(sb3T,sb1T,sb2T,new_matT+(new_ind+j)*sb3size,factorT,new_matT+ind*sb2Tsize,1.0); 
        }
    }

    new_localmatrix->index = new_index;
    new_localmatrix->start = new_start;
    new_localmatrix->end = new_end;
    new_localmatrix->mat = new_mat;
    new_localmatrix->matT = new_matT;
    new_localmatrix->ns = new_ns;
    new_localmatrix->nn = new_nn;
    new_localmatrix->lmap = new_lmap;

    FAMGReleaseHeap(FAMG_FROM_TOP);

    return 0;
                
}



static int LocalScalarProducts(const FAMGGraph &graph, struct FAMGMatrixLocal *localmatrix, double *&w1, double *&w2)
{
    double *w1ii, *w2ii, *w1ij, *w2ij; 
    int i, j, k, ind;
    int *index, *start, *end, ns, nn;
    double *mat, *matT, *hl1, *hl2;
 
    index = localmatrix->index;
    start = localmatrix->start;
    end = localmatrix->end;
    mat = localmatrix->mat;
    matT = localmatrix->matT;
    ns = localmatrix->ns;
    nn = localmatrix->nn;

    const FAMGSparseBlock *sbi = graph.Get_sb2Ptr();
    const FAMGSparseBlock *sbiT = graph.Get_sb2TPtr();
    const FAMGSparseBlock *sb0 = graph.Get_sb3Ptr();
    const FAMGSparseBlock *sb0T = graph.Get_sb3TPtr();
    const FAMGSparseBlock *sbia = graph.Get_sbiaPtr();
    const FAMGSparseBlock *sbiTa = graph.Get_sbiaTPtr();
    const FAMGSparseBlock *sb0a = graph.Get_sb0aPtr();
    const FAMGSparseBlock *sb0Ta = graph.Get_sb0aTPtr();

    const FAMGSparseVector *spi = graph.Get_spiPtr();
    const FAMGSparseVector *spiT = graph.Get_spiTPtr();
    const FAMGSparseVector *sp0 = graph.Get_sp0Ptr();
    const FAMGSparseVector *sp0T = graph.Get_sp0TPtr();
    const FAMGSparseVector *sp0i = graph.Get_sp0iPtr();
    const FAMGSparseVector *sp0iT = graph.Get_sp0iTPtr();
    
    
    short ncmpi = spi->Get_maxcomp()+1;
    short ncmpiT = spiT->Get_maxcomp()+1;
    short ncmp0 = sp0->Get_maxcomp()+1;
    short ncmp0T = sp0T->Get_maxcomp()+1;
    short ncmp0i = sp0i->Get_maxcomp()+1;
    short ncmp0iT = sp0iT->Get_maxcomp()+1;

    w1 = (double *) FAMGGetMem((ncmp0 + 2*(ns-1)*ncmp0i+(ns-1)*(ns-1)*ncmpi)*sizeof(double), FAMG_FROM_TOP); 
    w2 = (double *) FAMGGetMem((ncmp0T + 2*(ns-1)*ncmp0iT+(ns-1)*(ns-1)*ncmpiT)*sizeof(double), FAMG_FROM_TOP); 
    memset((void *)w1,0,(ncmp0 + 2*(ns-1)*ncmp0i+(ns-1)*(ns-1)*ncmpi)*sizeof(double));
    memset((void *)w2,0,(ncmp0T + 2*(ns-1)*ncmp0iT+(ns-1)*(ns-1)*ncmpiT)*sizeof(double));


    short sbsize = sbi->Get_maxoffset()+1;
    short sb0size = sb0->Get_maxoffset()+1;
    short sbTsize = sbiT->Get_maxoffset()+1;
    short sb0Tsize = sb0T->Get_maxoffset()+1;
   
    FAMGMarkHeap(FAMG_FROM_TOP);
    if(sbsize > sb0size) hl1 = (double *)  FAMGGetMem(nn*sbsize*sizeof(double), FAMG_FROM_TOP); 
    else  hl1 = (double *)  FAMGGetMem(nn*sb0size*sizeof(double), FAMG_FROM_TOP); 
    if(sbTsize > sb0Tsize) hl2 = (double *)  FAMGGetMem(nn*sbTsize*sizeof(double), FAMG_FROM_TOP); 
    else  hl2 = (double *)  FAMGGetMem(nn*sb0Tsize*sizeof(double), FAMG_FROM_TOP); 


    // i = 0 seperatly 
    
    memset((void *)hl1,0,nn*sb0size*sizeof(double));
    memset((void *)hl2,0,nn*sb0Tsize*sizeof(double));
    w1ii = w1;
    w2ii = w2;
    SparseBlockVSet(sp0,w1ii,0.0);
    SparseBlockVSet(sp0T,w2ii,0.0);
    for(ind = start[0]; ind < end[0]; ind++)
    {
        k = index[ind];
        SparseBlockMCopy(sb0,sb0,hl1+k*sb0size, mat+ind*sb0size, 1.0);
        SparseBlockMCopy(sb0T,sb0T,hl2+k*sb0Tsize, matT+ind*sb0Tsize, 1.0);
        SparseBlockRowAddScalProd(sp0,sb0,w1ii,mat+ind*sb0size,mat+ind*sb0size);
        SparseBlockRowAddScalProd(sp0T,sb0T,w2ii,matT+ind*sb0Tsize,matT+ind*sb0Tsize);
    }
    for(j = 1; j < ns; j++) 
    {
        w1ij = w1+ncmp0+(j-1)*ncmp0i;
        w2ij = w2+ncmp0T+(j-1)*ncmp0iT;
        SparseBlockVSet(sp0i,w1ij,0.0);
        SparseBlockVSet(sp0iT,w2ij,0.0);
        for(ind = start[j]; ind < end[j]; ind++)
        {
            k = index[ind];
            // adapt structure !!! 
            SparseBlockRowAddScalProd(sp0i,sb0a,sbia,w1ij,hl1+k*sb0size,mat+ind*sb0size);
            SparseBlockRowAddScalProd(sp0iT,sb0Ta,sbiTa,w2ij,hl2+k*sb0Tsize,matT+ind*sb0size);
        }
    }


    double *mat1 = mat + nn*(sb0size -sbsize); // help pointer
    double *mat1T = matT + nn*(sb0Tsize -sbTsize); // help pointer
    double *w1i0, *w2i0; // help pointer

    for(i = 1; i < ns; i++)
    {
        memset((void *)hl1,0,nn*sbsize*sizeof(double));
        memset((void *)hl2,0,nn*sbTsize*sizeof(double));
        w1i0 = w1+ncmp0+(ns-1)*ncmp0i+(i-1)*(ncmp0i+ncmpi*(ns-1));
        w2i0 = w2+ncmp0T+(ns-1)*ncmp0iT+(i-1)*(ncmp0iT+ncmpiT*(ns-1));
        w1ii = w1i0+ncmp0i+(i-1)*ncmpi;
        w2ii = w2i0+ncmp0iT+(i-1)*ncmpiT;
        SparseBlockVSet(spi,w1ii,0.0);
        SparseBlockVSet(spiT,w2ii,0.0);

        for(ind = start[i]; ind < end[i]; ind++)
        {
            k = index[ind];
            SparseBlockMCopy(sbi,sbi,hl1+k*sbsize, mat1+ind*sbsize, 1.0);
            SparseBlockMCopy(sbiT,sbiT,hl2+k*sbTsize, mat1T+ind*sbTsize, 1.0);
            SparseBlockRowAddScalProd(spi,sbi,w1ii,mat1+ind*sbsize,mat1+ind*sbsize);
            SparseBlockRowAddScalProd(spiT,sbiT,w2ii,mat1T+ind*sbTsize,mat1T+ind*sbTsize);
        }
        for(j = i+1; j < ns; j++) 
        {
            w1ij = w1i0+ncmp0i+(j-1)*ncmpi;
            w2ij = w2i0+ncmp0iT+(j-1)*ncmpiT;
            SparseBlockVSet(spi,w1ij,0.0);
            SparseBlockVSet(spiT,w2ij,0.0);
            for(ind = start[j]; ind < end[j]; ind++)
            {
                k = index[ind];
                SparseBlockRowAddScalProd(spi,sbi,w1ij,mat1+ind*sbsize,hl1+k*sbsize);
                SparseBlockRowAddScalProd(spiT,sbiT,w2ij,mat1T+ind*sbTsize,hl2+k*sbTsize);
            }
        }
    } 

    
    FAMGReleaseHeap(FAMG_FROM_TOP);

    return 0;
           
}

int FAMGGrid::GetLocalMinimum1(FAMGPaList *&palist, double *w1, double *w2, double *t1, double *t2, struct FAMGMatrixLocal *localmatrix)
{
    // try one parent scheme
    
    double *t10, *t1i, *t20, *t2i;
    double *w100, *w1i0, *w10i, *w1ii, *pi;
    double *w200, *w2i0, *w20i, *w2ii, *ri;
    double *normp, *normr;
    double norm, minnorm, norm1, norm2;
    double *coeffA[1], *coeffB[1];
    int i, pa[1], ns, *lmap;
    short ib;

    const double tol = FAMGGetParameter()->Gettol()*FAMGGetParameter()->Gettol();
    minnorm = 1e-5;
    minnorm = FAMGGetParameter()->Geterror1();
    minnorm = minnorm*minnorm;

    ns = localmatrix->ns;
    lmap = localmatrix->lmap;
    
    const FAMGSparseVector *stv = graph->Get_stvPtr();
    const FAMGSparseVector *stvT = graph->Get_stvTPtr();
    const FAMGSparseVector *swi = graph->Get_spiPtr();
    const FAMGSparseVector *swiT = graph->Get_spiTPtr();
    const FAMGSparseVector *sw0 = graph->Get_sp0Ptr();
    const FAMGSparseVector *sw0T = graph->Get_sp0TPtr();
    const FAMGSparseVector *sw0i = graph->Get_sp0iPtr();
    const FAMGSparseVector *sw0iT = graph->Get_sp0iTPtr();
    const FAMGSparseVector *sp = graph->Get_spPtr();
    const FAMGSparseVector *sr = graph->Get_srPtr();
   
    short nb = sp->Get_n();
    short psize = sp->Get_maxcomp()+1;
    short rsize = sr->Get_maxcomp()+1;
    pi = new double[psize];
    ri = new double[rsize];
    normp = new double[psize];
    normr = new double[rsize];
    

    short *pcomputed = new short[psize];
    short *rcomputed = new short[rsize];

    short *pcomp = sp->Get_comp();
    short *rcomp = sr->Get_comp();
    short ncmpi = swi->Get_maxcomp()+1;
    short ncmpiT = swiT->Get_maxcomp()+1;
    short ncmp0 = sw0->Get_maxcomp()+1;
    short ncmp0T = sw0T->Get_maxcomp()+1;
    short ncmp0i = sw0i->Get_maxcomp()+1;
    short ncmp0iT = sw0iT->Get_maxcomp()+1;
    short vsize = stv->Get_maxcomp()+1;
    short vTsize = stvT->Get_maxcomp()+1;
    
       
    t10 = (t1+0*vsize);
    t20 = (t2+0*vTsize);
    w100 = w1;
    w200 = w2;
                 
    for(i = 1; i < ns; i++)
    {
        if( (graph->GetNode(lmap[i]))->IsFGNode()) continue;

        t1i = (t1+i*vsize);
        t2i = (t2+i*vTsize);
            
        w10i = (w1+ncmp0+(i-1)*ncmp0i);
        w20i = (w2+ncmp0T+(i-1)*ncmp0iT);

        w1i0 = w1+ncmp0+(ns-1)*ncmp0i+(i-1)*(ncmp0i+ncmpi*(ns-1));
        w2i0 = w2+ncmp0T+(ns-1)*ncmp0iT+(i-1)*(ncmp0iT+ncmpiT*(ns-1));
        w1ii = (w1i0+ncmp0i+(i-1)*ncmpi);
        w2ii = (w2i0+ncmp0iT+(i-1)*ncmpiT);

        memset((void*)pcomputed,0,psize*sizeof(short));
        memset((void*)rcomputed,0,rsize*sizeof(short));
        for(ib = 0; ib < nb; ib++)
        {
            if((pcomputed[pcomp[ib]]) && (rcomputed[rcomp[ib]])) continue;
            pcomputed[pcomp[ib]] = rcomputed[rcomp[ib]] = 1;

            pi[pcomp[ib]] = t10[stv->Get_comp(ib)]/t1i[stv->Get_comp(ib)];
            normp[pcomp[ib]] = w100[sw0->Get_comp(ib)] + 
                               w1ii[swi->Get_comp(ib)]*pi[pcomp[ib]]*pi[pcomp[ib]]
                               -2.0*w10i[sw0i->Get_comp(ib)]*pi[pcomp[ib]];

            ri[rcomp[ib]] = t20[stvT->Get_comp(ib)]/t2i[stvT->Get_comp(ib)];
            normr[rcomp[ib]] = w200[sw0T->Get_comp(ib)] + 
                               w2ii[swiT->Get_comp(ib)]*ri[rcomp[ib]]*ri[rcomp[ib]]
                               -2.0*w20i[sw0iT->Get_comp(ib)]*ri[rcomp[ib]];
        }
        norm1 = norm2 = 0.0;
        for(ib = 0; ib < nb; ib++)
        {
            // norm1 += normp[pcomp[ib]];
            // norm2 += normr[rcomp[ib]];
            if(Abs(normp[pcomp[ib]]) > norm1) norm1 = Abs(normp[pcomp[ib]]);
            if(Abs(normr[rcomp[ib]]) > norm2) norm2 = Abs(normr[rcomp[ib]]);
        }
            

        norm = norm1*norm2;
        if(tol*norm < minnorm)
        {
            if(norm < minnorm) minnorm = norm;
            pa[0] = lmap[i]; 
            coeffA[0] = pi; 
            coeffB[0] = ri; 
            if(graph->SavePaList(palist,1,pa,coeffA,coeffB,norm)) return 1;
        }
    }
                
    graph->CorrectPaList(palist,minnorm/tol);

    delete pcomputed;
    delete rcomputed;
    delete pi;
    delete ri;
    delete normp;
    delete normr;

    return 0;
}


int FAMGGrid::GetLocalMinimum2(FAMGPaList *&palist, double *w1, double *w2, double *t1, double *t2, struct FAMGMatrixLocal *localmatrix)
{
    // two parents scheme
   
    double *t10, *t1i, *t1j, *t20, *t2i, *t2j;
    double *w100, *w1i0, *w10i, *w10j, *w1j0, *w1ij, *w1ii, *w1jj, *pi, *pj;
    double *w200, *w2i0, *w20i, *w20j, *w2j0, *w2ij, *w2ii, *w2jj, *ri, *rj;
    double zaehler1, nenner1, laenge1, norm1, *normp;
    double zaehler2, nenner2, laenge2, norm2, *normr;
    double norm, minnorm, *coeffA[2], *coeffB[2];
    int i, j, pa[2], *lmap, ns;
    short ib;

    const double tol = FAMGGetParameter()->Gettol()*FAMGGetParameter()->Gettol();
    minnorm = 1.0;
    minnorm = FAMGGetParameter()->Geterror2();
    minnorm = minnorm*minnorm;

    ns = localmatrix->ns;
    lmap = localmatrix->lmap;

    const FAMGSparseVector *stv = graph->Get_stvPtr();
    const FAMGSparseVector *stvT = graph->Get_stvTPtr();
    const FAMGSparseVector *swi = graph->Get_spiPtr();
    const FAMGSparseVector *swiT = graph->Get_spiTPtr();
    const FAMGSparseVector *sw0 = graph->Get_sp0Ptr();
    const FAMGSparseVector *sw0T = graph->Get_sp0TPtr();
    const FAMGSparseVector *sw0i = graph->Get_sp0iPtr();
    const FAMGSparseVector *sw0iT = graph->Get_sp0iTPtr();
    const FAMGSparseVector *sp = graph->Get_spPtr();
    const FAMGSparseVector *sr = graph->Get_srPtr();

    short nb = sp->Get_n();
    short psize = sp->Get_maxcomp()+1;
    short rsize = sr->Get_maxcomp()+1;
    pi = new double[psize];
    pj = new double[psize];
    ri = new double[rsize];
    rj = new double[rsize];
    normp = new double[psize];
    normr = new double[rsize];
    
    short *pcomputed = new short[psize];
    short *rcomputed = new short[rsize];

    short *pcomp = sp->Get_comp();
    short *rcomp = sr->Get_comp();
    short ncmpi = swi->Get_maxcomp()+1;
    short ncmpiT = swiT->Get_maxcomp()+1;
    short ncmp0 = sw0->Get_maxcomp()+1;
    short ncmp0T = sw0T->Get_maxcomp()+1;
    short ncmp0i = sw0i->Get_maxcomp()+1;
    short ncmp0iT = sw0iT->Get_maxcomp()+1;
    short vsize = stv->Get_maxcomp()+1;
    short vTsize = stvT->Get_maxcomp()+1;
    
       
    t10 = (t1+0*vsize);
    t20 = (t2+0*vTsize);
    w100 = w1;
    w200 = w2;

    for(i = 1; i < ns; i++)
    {
        if( (graph->GetNode(lmap[i]))->IsFGNode()) continue;

        t1i = (t1+i*vsize);
        t2i = (t2+i*vTsize);
            
        w10i = (w1+ncmp0+(i-1)*ncmp0i);
        w20i = (w2+ncmp0T+(i-1)*ncmp0iT);

        w1i0 = w1+ncmp0+(ns-1)*ncmp0i+(i-1)*(ncmp0i+ncmpi*(ns-1));
        w2i0 = w2+ncmp0T+(ns-1)*ncmp0iT+(i-1)*(ncmp0iT+ncmpiT*(ns-1));
        w1ii = (w1i0+ncmp0i+(i-1)*ncmpi);
        w2ii = (w2i0+ncmp0iT+(i-1)*ncmpiT);

        for(j = i+1; j < ns; j++)
        {
            if( (graph->GetNode(lmap[j]))->IsFGNode()) continue;

            t1j = (t1+j*vsize);
            t2j = (t2+j*vTsize);

            w10j = (w1+ncmp0+(j-1)*ncmp0i);
            w20j = (w2+ncmp0T+(j-1)*ncmp0iT);

            w1j0 = w1+ncmp0+(ns-1)*ncmp0i+(j-1)*(ncmp0i+ncmpi*(ns-1));
            w2j0 = w2+ncmp0T+(ns-1)*ncmp0iT+(j-1)*(ncmp0iT+ncmpiT*(ns-1));
            w1jj = (w1j0+ncmp0i+(j-1)*ncmpi);
            w2jj = (w2j0+ncmp0iT+(j-1)*ncmpiT);
            w1ij = (w1i0+ncmp0i+(j-1)*ncmpi);
            w2ij = (w2i0+ncmp0iT+(j-1)*ncmpiT);

            memset((void*)pcomputed,0,psize*sizeof(short));
            memset((void*)rcomputed,0,rsize*sizeof(short));
            for(ib = 0; ib < nb; ib++)
            {
                if(pcomputed[pcomp[ib]]) continue;
                pcomputed[pcomp[ib]] = 1;

                nenner1 = t1i[stv->Get_comp(ib)]*t1i[stv->Get_comp(ib)]*w1jj[swi->Get_comp(ib)]
                        - 2.0*t1i[stv->Get_comp(ib)]*t1j[stv->Get_comp(ib)]*w1ij[swi->Get_comp(ib)]
                        + t1j[stv->Get_comp(ib)]*t1j[stv->Get_comp(ib)]*w1ii[swi->Get_comp(ib)];
                  if(nenner1 < -1.0e+20) 
                  {
                      cout << "Warning! GetLocalMinimum: nenner1 too small." << endl << flush;
                      continue;
                  }
            
            
                // construct prolongation
                if( Abs(t1j[stv->Get_comp(ib)]) > Abs(t1i[stv->Get_comp(ib)]))
                {
                    zaehler1 = t1j[stv->Get_comp(ib)]*t1j[stv->Get_comp(ib)]*w10i[sw0i->Get_comp(ib)]
                              -t1j[stv->Get_comp(ib)]*t1i[stv->Get_comp(ib)]*w10j[sw0i->Get_comp(ib)]
                              -t10[stv->Get_comp(ib)]*t1j[stv->Get_comp(ib)]*w1ij[swi->Get_comp(ib)]
                              +t10[stv->Get_comp(ib)]*t1i[stv->Get_comp(ib)]*w1jj[swi->Get_comp(ib)];
                    laenge1 = t1j[stv->Get_comp(ib)]*t1j[stv->Get_comp(ib)]*w100[sw0->Get_comp(ib)]
                             -2.0*t10[stv->Get_comp(ib)]*t1j[stv->Get_comp(ib)]*w10j[sw0i->Get_comp(ib)]
                             +t10[stv->Get_comp(ib)]*t10[stv->Get_comp(ib)]*w1jj[sw0i->Get_comp(ib)];

                    pi[pcomp[ib]] = zaehler1/nenner1;
                    pj[pcomp[ib]] = (t10[stv->Get_comp(ib)] - pi[pcomp[ib]]*t1i[stv->Get_comp(ib)])/t1j[stv->Get_comp(ib)];
                    normp[pcomp[ib]] = (laenge1 - zaehler1*zaehler1/nenner1)/(t1j[stv->Get_comp(ib)]*t1j[stv->Get_comp(ib)]);
                }

                else
                {
                    zaehler1 = t1i[stv->Get_comp(ib)]*t1i[stv->Get_comp(ib)]*w10j[sw0i->Get_comp(ib)]
                              -t1i[stv->Get_comp(ib)]*t1j[stv->Get_comp(ib)]*w10i[sw0i->Get_comp(ib)]
                              -t10[stv->Get_comp(ib)]*t1i[stv->Get_comp(ib)]*w1ij[swi->Get_comp(ib)]
                              +t10[stv->Get_comp(ib)]*t1j[stv->Get_comp(ib)]*w1ii[swi->Get_comp(ib)];
                    laenge1 = t1i[stv->Get_comp(ib)]*t1i[stv->Get_comp(ib)]*w100[sw0->Get_comp(ib)]
                             -2.0*t10[stv->Get_comp(ib)]*t1i[stv->Get_comp(ib)]*w10i[sw0i->Get_comp(ib)]
                             +t10[stv->Get_comp(ib)]*t10[stv->Get_comp(ib)]*w1ii[sw0i->Get_comp(ib)];
                    pj[pcomp[ib]] = zaehler1/nenner1;
                    pi[pcomp[ib]] = (t10[stv->Get_comp(ib)] - pj[pcomp[ib]]*t1j[stv->Get_comp(ib)])/t1i[stv->Get_comp(ib)];
                    normp[pcomp[ib]] = (laenge1 - zaehler1*zaehler1/nenner1)/(t1i[stv->Get_comp(ib)]*t1i[stv->Get_comp(ib)]);
                }

            }

            // same procedure for the restriction

            for(ib = 0; ib < nb; ib++)
            {
                if(rcomputed[rcomp[ib]]) continue;
                rcomputed[rcomp[ib]] = 1;

                nenner2 = t2i[stvT->Get_comp(ib)]*t2i[stvT->Get_comp(ib)]*w2jj[swiT->Get_comp(ib)]
                        - 2.0*t2i[stvT->Get_comp(ib)]*t2j[stvT->Get_comp(ib)]*w2ij[swiT->Get_comp(ib)]
                        + t2j[stvT->Get_comp(ib)]*t2j[stvT->Get_comp(ib)]*w2ii[swiT->Get_comp(ib)];
                 if(nenner2 < -1e+15) 
                 {
                     cout << "Warning! GetLocalMinimum: nenner2 too small." << endl << flush;
                     continue;
                 }
            
            
                // construct restriction
                if( Abs(t2j[stvT->Get_comp(ib)]) > Abs(t2i[stvT->Get_comp(ib)]))
                {
                    zaehler2 = t2j[stvT->Get_comp(ib)]*t2j[stvT->Get_comp(ib)]*w20i[sw0iT->Get_comp(ib)]
                              -t2j[stvT->Get_comp(ib)]*t2i[stvT->Get_comp(ib)]*w20j[sw0iT->Get_comp(ib)]
                              -t20[stvT->Get_comp(ib)]*t2j[stvT->Get_comp(ib)]*w2ij[swiT->Get_comp(ib)]
                              +t20[stvT->Get_comp(ib)]*t2i[stvT->Get_comp(ib)]*w2jj[swiT->Get_comp(ib)];
                    laenge2 = t2j[stvT->Get_comp(ib)]*t2j[stvT->Get_comp(ib)]*w200[sw0T->Get_comp(ib)]
                             -2.0*t20[stvT->Get_comp(ib)]*t2j[stvT->Get_comp(ib)]*w20j[sw0iT->Get_comp(ib)]
                             +t20[stvT->Get_comp(ib)]*t20[stvT->Get_comp(ib)]*w2jj[sw0iT->Get_comp(ib)];

                    ri[rcomp[ib]] = zaehler2/nenner2;
                    rj[rcomp[ib]] = (t20[stvT->Get_comp(ib)] - ri[rcomp[ib]]*t2i[stvT->Get_comp(ib)])/t2j[stvT->Get_comp(ib)];
                    normr[rcomp[ib]] = (laenge2 - zaehler2*zaehler2/nenner2)/(t2j[stvT->Get_comp(ib)]*t2j[stvT->Get_comp(ib)]);
                }

                else
                {
                    zaehler2 = t2i[stvT->Get_comp(ib)]*t2i[stvT->Get_comp(ib)]*w20j[sw0iT->Get_comp(ib)]
                              -t2i[stvT->Get_comp(ib)]*t2j[stvT->Get_comp(ib)]*w20i[sw0iT->Get_comp(ib)]
                              -t20[stvT->Get_comp(ib)]*t2i[stvT->Get_comp(ib)]*w2ij[swiT->Get_comp(ib)]
                              +t20[stvT->Get_comp(ib)]*t2j[stvT->Get_comp(ib)]*w2ii[swiT->Get_comp(ib)];
                    laenge2 = t2i[stvT->Get_comp(ib)]*t2i[stvT->Get_comp(ib)]*w200[sw0T->Get_comp(ib)]
                             -2.0*t20[stvT->Get_comp(ib)]*t2i[stvT->Get_comp(ib)]*w20i[sw0iT->Get_comp(ib)]
                             +t20[stvT->Get_comp(ib)]*t20[stvT->Get_comp(ib)]*w2ii[sw0iT->Get_comp(ib)];
                    rj[rcomp[ib]] = zaehler2/nenner2;
                    ri[rcomp[ib]] = (t20[stvT->Get_comp(ib)] - rj[rcomp[ib]]*t2j[stvT->Get_comp(ib)])/t2i[stvT->Get_comp(ib)];
                    normr[rcomp[ib]] = (laenge2 - zaehler2*zaehler2/nenner2)/(t2i[stvT->Get_comp(ib)]*t2i[stvT->Get_comp(ib)]);
                }
            }



            //pi[0] = pi[1]; pj[0] = pj[1];
            //ri[0] = ri[1]; rj[0] = rj[1];
            //pi[0] = 0.0; pj[0] = 0.0;
            //ri[0] = 0.0; rj[0] = 0.0;

           
            norm1 = norm2 = 0.0;
            for(ib = 0; ib < nb; ib++)
            {
                norm1 += normp[pcomp[ib]];
                norm2 += normr[rcomp[ib]];
                //if(Abs(normp[pcomp[ib]]) > norm1) norm1 = Abs(normp[pcomp[ib]]);
                //if(Abs(normr[rcomp[ib]]) > norm2) norm2 = Abs(normr[rcomp[ib]]);
            }
     
            norm = Abs(norm1*norm2);
            if(tol*norm < minnorm)
            {
                if(norm < minnorm) minnorm = norm;
                pa[0] = lmap[i]; pa[1] = lmap[j]; 
                coeffA[0] = pi; coeffA[1] = pj;
                coeffB[0] = ri; coeffB[1] = rj; 
                if(graph->SavePaList(palist,2,pa,coeffA,coeffB,norm)) return 1;
            }
        }
    }
                
    graph->CorrectPaList(palist,minnorm/tol);

    delete pcomputed;
    delete rcomputed;
    delete pi;
    delete pj;
    delete ri;
    delete rj;
    delete normp;
    delete normr;

    return 0;
}


int FAMGGrid::AnalyseNode6(const FAMGVectorEntry &veci, FAMGPaList *&palist)
{  
    FAMGMatrixLocal localmatrix, new_localmatrix;
    double *tv1, *tv2, *w1, *w2, normr, norml;
    int t, dt;

    palist = NULL;
 
    FAMGMarkHeap(FAMG_FROM_TOP); 

    t = (int) clock();
    ConstructLocalMatrixA(veci,&localmatrix,tv1,tv2,normr, norml); 
    dt = (int) clock() - t;
    time1 += dt;

    if((normr < 1e-25) || (norml < 1e-25))
    {
         if(graph->SavePaList(palist,0,NULL,NULL,NULL,0.0)) return 1;
         FAMGReleaseHeap(FAMG_FROM_TOP);
         return 0;
    }

    t = (int) clock(); 
    if(n < 0) 
        LocalMatMultA(*GetGraph(),&new_localmatrix,&localmatrix);
    else 
        LocalMatMult(*GetGraph(),&new_localmatrix,&localmatrix);
    dt = (int) clock() - t;
    time2 += dt;


    t = (int) clock();
    LocalScalarProducts(*GetGraph(),&new_localmatrix,w1,w2);
    dt = (int) clock() - t;
    time3 += dt;


    t = (int) clock();
    // GetLocalMinimum1(palist,w1,w2,tv1,tv2,&new_localmatrix);
    if(palist == 0) GetLocalMinimum2(palist,w1,w2,tv1,tv2,&new_localmatrix);
    dt = (int) clock() - t;
    time4 += dt;

    FAMGReleaseHeap(FAMG_FROM_TOP);


    if(veci.GetIndex() == n-1) 
    {
        cout << "time1 " << time1 << endl;
        cout << "time2 " << time2 << endl;
        cout << "time3 " << time3 << endl;
        cout << "time4 " << time4 << endl;
        time1 = time2 = time3 = time4 = 0;
    }

     return 0;

}
 
#endif
