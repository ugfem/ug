/****************************************************************************/
/*																			*/
/* File:      famg_approxNew.C													*/
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


struct FAMGMatrixLocal
{
    int nn;
    int ns;
    int *index;
    int *start;
    int *end;
    double *mat;
    double *matT;
    double *vec;
    double *vecT;
    int *lmap;
};

static int time1 = 0;
static int time2 = 0;
static int time3 = 0;
static int time4 = 0;


int FAMGGrid::ConstructLocalMatrix(const FAMGVectorEntry &veci, struct FAMGMatrixLocal *localmatrix, double *&tv1, double *&tv2, double &normr, double &norml)
{  
    const FAMGGraph &graph = *GetGraph();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
	const FAMGMatrixAlg &A = *GetTmpMatrix();
    FAMGMatrixIter mi_iter(A,veci); 
	FAMGVectorEntry vecj, veck;
	FAMGMatrixEntry matij, matjk, matkk;
    double miiinv, mjjinv, mkkinv, mij, mji, mjk, mkj, t1, t2, norm1, norm2;
    int i, j, k, lid, ind, ix;
    int nn, ns, nm, nm1, nm2;
    int  *start, *end, *index, *lmap;
    double *mat, *matT, *vec, *vecT;

    const double omega = FAMGGetParameter()->Getomegar();

    // set local IDs
    graph.GetNode(veci.GetIndex())->SetLocalId(0);
    mi_iter(matij);   
    nn = 1; 
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
    norm1 = norm2 = 0.0;
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
    mat = (double *) FAMGGetMem((nm)*sizeof(double), FAMG_FROM_TOP); 
    matT = (double *) FAMGGetMem((nm)*sizeof(double), FAMG_FROM_TOP); 
    vec = (double *)  FAMGGetMem((ns)*sizeof(double), FAMG_FROM_TOP); 
    vecT = (double *)  FAMGGetMem((ns)*sizeof(double), FAMG_FROM_TOP); 
    tv1 = (double *)  FAMGGetMem((ns)*sizeof(double), FAMG_FROM_TOP); 
    tv2 = (double *)  FAMGGetMem((ns)*sizeof(double), FAMG_FROM_TOP); 
    lmap = (int *)  FAMGGetMem((nn)*sizeof(int), FAMG_FROM_TOP); 


    // copy matrix
    // this function uses a couple of strange tricks
    // i guess it works :-)
    i = 1; 
    lmap[0] = veci.GetIndex();;
    start[0] = 0; end[0] = ns; start[1] = ns; 
    start[ns] = ns+nm1+nm2+ns-1; end[ns] = start[ns] + 1;
    mat[0] = 1.0-omega;
    matT[0] = 1.0-omega;
    index[0] = 0;
    tv1[0] = 0.0; tv2[0] = 0.0;

    mi_iter.reset(); mi_iter(matij);
    miiinv = 1.0/A[matij];
    while(mi_iter(matij))
    {
        ind = start[i];
        vecj = matij.dest();
        mij = A[matij]*miiinv;
        mji = A.GetAdjData(matij)*miiinv;
        t1 = tvA[vecj]; t2 = tvB[vecj];
        tv1[0] -= mij*t1;
        tv2[0] -= mji*t2;
        norm1 += Abs(mij);
        norm2 += Abs(mji);
        if(matij.is_strong())
        {
            lmap[i] = vecj.GetIndex();
            vec[i] = -mij;
            vecT[i] = -mji;
            tv1[i] = t1; 
            tv2[i] = t2;
            FAMGMatrixIter mj_iter(A,vecj);  
            mj_iter(matjk);
            mjjinv = 1.0/A[matjk]; 
            mat[i] = -omega*mij;
            matT[i] = -omega*mji;
            index[i] = i;
            mat[ind] = 1.0-omega;
            matT[ind] = 1.0-omega;
            index[ind] = i;
            ind++;
            while(mj_iter(matjk))
            {
                if(matjk.is_strong())
                {
                    veck = matjk.dest();
                    k = veck.GetIndex();
                    mjk = A[matjk];
                    mkj = A.GetAdjData(matjk);
                    lid = graph.GetNode(k)->GetLocalId();
                    mat[ind] = -omega*mjk*mjjinv;
                    matT[ind] = -omega*mkj*mjjinv;
                    index[ind] = lid;
                    ind++;
                    
                    // construct rows for the non-neighbor nodes
                    if(lid >= ns)
                    {

                        FAMGMatrixIter mk_iter(A,veck); mk_iter(matkk);
                        mkkinv = 1.0/A[matkk];
                        ix = end[lid];
                        mat[ix] = -omega*mkj*mkkinv;
                        matT[ix] = -omega*mjk*mkkinv;
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
                                mat[start[lid]] = 1.0-omega;
                                matT[start[lid]] = 1.0-omega;
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


    return 0;

}


static int LocalMatMult(struct FAMGMatrixLocal *new_localmatrix, struct FAMGMatrixLocal *localmatrix)
{
    double v1ik, v2ik, v1kj, v2kj, factor, factorT;
    int i, j, k, ind, new_ind, ix;
    int new_nn, new_ns, nn, ns, nm;
    int *new_start, *new_end, *new_index, *new_lmap, *start, *end, *index, *lmap;
    double *new_mat, *new_matT, *mat, *matT, *vec, *vecT;

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

    nm = ns*nn; // just to be sure

    // allocate start, end, index, mat, matT
    new_start = (int *) FAMGGetMem((ns+1)*sizeof(int), FAMG_FROM_TOP); 
    new_end = (int *) FAMGGetMem((ns+1)*sizeof(int), FAMG_FROM_TOP); 
    new_index = (int *) FAMGGetMem((nm)*sizeof(int), FAMG_FROM_TOP); 
    new_mat = (double *) FAMGGetMem((nm)*sizeof(double), FAMG_FROM_TOP); 
    new_matT = (double *) FAMGGetMem((nm)*sizeof(double), FAMG_FROM_TOP); 
    new_lmap = lmap;
    new_nn = nn;
    new_ns = ns;

    
    memset((void *)new_mat,0,(nm)*sizeof(double));
    memset((void *)new_matT,0,(nm)*sizeof(double));
 
    FAMGMarkHeap(FAMG_FROM_TOP);
    double *helpvec1 = (double *) FAMGGetMem((nn)*sizeof(double), FAMG_FROM_TOP);
    double *helpvec2 = (double *) FAMGGetMem((nn)*sizeof(double), FAMG_FROM_TOP);

    new_start[0] = 0; new_end[0] = nn;
    new_start[1] = nn;
    for(i = 1; i < ns; i++)
    {
        memset((void *)helpvec1,0,(nn)*sizeof(double));
        memset((void *)helpvec2,0,(nn)*sizeof(double));

        for(ind = start[i]; ind < end[i]; ind++)
        {
            k = index[ind];
            v1ik = mat[ind];
            v2ik = matT[ind];

            for(ix = start[k]; ix < end[k]; ix++)
            {
                j = index[ix];
                v1kj = mat[ix];
                v2kj = matT[ix];

                helpvec1[j] += v1ik*v1kj;
                helpvec2[j] += v2ik*v2kj;
            }
        }
        new_ind = new_start[i];
        for(j = 0; j < nn; j++)
        {
            if((helpvec1[j] != 0.0) || (helpvec2[j] != 0.0))
            {
                new_mat[new_ind] = helpvec1[j];
                new_matT[new_ind] = helpvec2[j];
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
        factor = vec[i];
        factorT = vecT[i];
        for(ind = new_start[i]; ind < new_end[i]; ind++)
        {
            j = new_index[ind];
            new_mat[new_ind+j] += factor*new_mat[ind];
            new_matT[new_ind+j] += factorT*new_matT[ind];
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



static int LocalScalarProducts(struct FAMGMatrixLocal *localmatrix, double *&w1, double *&w2)
{
    double w1ii, w2ii, w1ij, w2ij, v1k, v2k; 
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

    w1 = (double *) FAMGGetMem((ns)*(ns)*sizeof(double), FAMG_FROM_TOP); 
    w2 = (double *) FAMGGetMem((ns)*(ns)*sizeof(double), FAMG_FROM_TOP); 


    FAMGMarkHeap(FAMG_FROM_TOP);
    hl1 = (double *)  FAMGGetMem((nn)*sizeof(double), FAMG_FROM_TOP); 
    hl2 = (double *)  FAMGGetMem((nn)*sizeof(double), FAMG_FROM_TOP); 

    memset((void *)w1,0,(ns)*(ns)*sizeof(double));
    memset((void *)w2,0,(ns)*(ns)*sizeof(double));

 
    for(i = 0; i < ns; i++)
    {
        memset((void *)hl1,0,(nn)*sizeof(double));
        memset((void *)hl2,0,(nn)*sizeof(double));
        w1ii = w2ii = 0.0;
        for(ind = start[i]; ind < end[i]; ind++)
        {
            k = index[ind];
            v1k = mat[ind];
            v2k = matT[ind];
            hl1[k] = v1k;
            hl2[k] = v2k;
            w1ii += v1k*v1k;
            w2ii += v2k*v2k;
        }
        w1[ns*i+i] = w1ii;
        w2[ns*i+i] = w2ii;
        for(j = i+1; j < ns; j++) 
        {
            w1ij = w2ij = 0.0;
            for(ind = start[j]; ind < end[j]; ind++)
            {
                k = index[ind];
                w1ij += mat[ind]*hl1[k];
                w2ij += matT[ind]*hl2[k];
            }
            w1[ns*i+j] = w1ij;
            w2[ns*i+j] = w2ij;
        }
    } 


    FAMGReleaseHeap(FAMG_FROM_TOP);

    return 0;
           
}

int FAMGGrid::GetLocalMinimum1(FAMGPaList *&palist, double *w1, double *w2, double *t1, double *t2, struct FAMGMatrixLocal *localmatrix)
{
    // try one parent scheme
    
    double t10, t1i, t20, t2i;
    double w100, w1i0, w1ii, pi;
    double w200, w2i0, w2ii, ri;
    double w2h00, w2hi0, w2hii;
    double norm, minnorm, norm1, norm2;
    double coeffA[1], coeffB[1];
    int i, pa[1], ns, *lmap;

    const double tol = FAMGGetParameter()->Gettol()*FAMGGetParameter()->Gettol();
    minnorm = 1e-5;
    minnorm = FAMGGetParameter()->Geterror1();
    minnorm = minnorm*minnorm;

    ns = localmatrix->ns;
    lmap = localmatrix->lmap;


    t10 = t1[0]; w100 = w1[0];
    t20 = t2[0]; w200 = w2[0];
    for(i = 1; i < ns; i++)
    {
        t1i = t1[i]; w1i0 = w1[i]; w1ii = w1[ns*i+i];
        t2i = t2[i]; w2i0 = w2[i]; w2ii = w2[ns*i+i];
        pi = t10/t1i;
        norm1 = w100 + w1ii*pi*pi -2.0*w1i0*pi;

        /*
          this part is not yet checked !!
          
          double v200, v20i, v2i0, v2ii, v2ji, hf;
          v200 = matT[start[ns]+ns]; v20i = matT[start[ns]+i-1]; v20j = matT[start[ns]+j-1];
          v2i0 = matT[start[i-1]+ns]; v2ii = matT[start[i-1]+i-1]; v2ij = matT[start[i-1]+j-1];
          v2j0 = matT[start[j-1]+ns]; v2ji = matT[start[j-1]+i-1]; v2jj = matT[start[j-1]+j-1];

          hf = pi*pi  - 1.0;
          w2h00 = w200 + hf*v200*v200 + 2.0*v200*pi*v20i;
          w2hii = w2ii + hf*v2i0*v2i0 + 2.0*v2i0*pi*v2ii;
          w2hi0 = w2i0 + hf*v2i0*v200 + v2i0*pi*v20i + v200*pi*v2ii;
            
        */          
            
        w2h00 = w200;w2hii = w2ii;w2hi0 = w2i0;

        ri = t20/t2i;
        norm2 = w2h00 + w2hii*ri*ri -2.0*w2hi0*ri;

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

    return 0;
}


int FAMGGrid::GetLocalMinimum2(FAMGPaList *&palist, double *w1, double *w2, double *t1, double *t2, struct FAMGMatrixLocal *localmatrix)
{
    // two parents scheme
   
    double t10, t1i, t1j, t20, t2i, t2j;
    double w100, w1i0, w1j0, w1ij, w1ii, w1jj, pi, pj;
    double w200, w2i0, w2j0, w2ij, w2ii, w2jj, ri, rj;
    double w2h00, w2hi0, w2hj0, w2hij, w2hii, w2hjj;
    double zaehler1, nenner1, laenge1, norm1;
    double zaehler2, nenner2, laenge2, norm2;
    double norm, minnorm, coeffA[2], coeffB[2];
    int i, j, pa[2], *lmap, ns;

    const double tol = FAMGGetParameter()->Gettol()*FAMGGetParameter()->Gettol();
    minnorm = 1.0;
    minnorm = FAMGGetParameter()->Geterror2();
    minnorm = minnorm*minnorm;

    ns = localmatrix->ns;
    lmap = localmatrix->lmap;

    t10 = t1[0]; w100 = w1[0];
    t20 = t2[0]; w200 = w2[0];
    for(i = 1; i < ns; i++)
    {
        t1i = t1[i]; w1i0 = w1[i]; w1ii = w1[ns*i+i];
        t2i = t2[i]; w2i0 = w2[i]; w2ii = w2[ns*i+i];
        for(j = i+1; j < ns; j++)
        {
            t1j = t1[j]; w1j0 = w1[j]; w1jj = w1[ns*j+j];
            t2j = t2[j]; w2j0 = w2[j]; w2jj = w2[ns*j+j];
            w1ij = w1[ns*i+j];
            w2ij = w2[ns*i+j];
            
            nenner1 = t1i*t1i*w1jj - 2.0*t1i*t1j*w1ij + t1j*t1j*w1ii;
            //  if(nenner1 < 1e-15) continue;
            
            
            // construct prolongation
            if( Abs(t1j) > Abs(t1i))
            {
                zaehler1 = t1j*t1j*w1i0 - t1j*t1i*w1j0 - t10*t1j*w1ij + t10*t1i*w1jj;
                laenge1 = w100*t1j*t1j - 2.0*t10*t1j*w1j0 + t10*t10*w1jj;
                pi = zaehler1/nenner1;
                pj = (t10 - pi*t1i)/t1j;
                norm1 = (laenge1 - zaehler1*zaehler1/nenner1)/(t1j*t1j);
            }

            else
            {
                zaehler1 = t1i*t1i*w1j0 - t1i*t1j*w1i0 - t10*t1i*w1ij + t10*t1j*w1ii;
                laenge1 = w100*t1i*t1i - 2.0*t10*t1i*w1i0 + t10*t10*w1ii;
                pj = zaehler1/nenner1;
                pi = (t10 - pj*t1j)/t1i;
                norm1 = (laenge1 - zaehler1*zaehler1/nenner1)/(t1i*t1i);
           }

            /*
            double v200, v20i, v20j, v2i0, v2ii, v2ij, v2j0, v2ji, v2jj, hf;

            v200 = matT[start[ns]+ns]; v20i = matT[start[ns]+i-1]; v20j = matT[start[ns]+j-1];
            v2i0 = matT[start[i-1]+ns]; v2ii = matT[start[i-1]+i-1]; v2ij = matT[start[i-1]+j-1];
            v2j0 = matT[start[j-1]+ns]; v2ji = matT[start[j-1]+i-1]; v2jj = matT[start[j-1]+j-1];

            hf = pi*pi + pj*pj - 1.0;
            w2h00 = w200 + hf*v200*v200 + 2.0*v200*(pi*v20i+ pj*v20j);
            w2hii = w2ii + hf*v2i0*v2i0 + 2.0*v2i0*(pi*v2ii + pj*v2ij);
            w2hjj = w2jj + hf*v2j0*v2j0 + 2.0*v2j0*(pi*v2ji + pj*v2jj) ;
            w2hi0 = w2i0 + hf*v2i0*v200 + v2i0*(pi*v20i + pj*v20j) + v200*(pi*v2ii + pj*v2ij);
            w2hj0 = w2j0 + hf*v2j0*v200 + v2j0*(pi*v20i + pj*v20j) + v200*(pi*v2ji + pj*v2jj);
            w2hij = w2ij + hf*v2i0*v2j0 + v2i0*(pi*v2ji + pj*v2jj) + v2j0*(pi*v2ii + pj*v2ij);
            
            */          
            
            w2h00 = w200;w2hii = w2ii;w2hjj = w2jj;w2hi0 = w2i0;w2hj0 = w2j0; w2hij = w2ij;

            // same procedure for the restriction

            nenner2 = t2i*t2i*w2hjj - 2.0*t2i*t2j*w2hij + t2j*t2j*w2hii;
            //  if(nenner2 < 1e-15) continue;

            if( Abs(t2j) > Abs(t2i))
            {
                zaehler2 = t2j*t2j*w2hi0 - t2j*t2i*w2hj0 - t20*t2j*w2hij + t20*t2i*w2hjj;
                laenge2 = w2h00*t2j*t2j - 2.0*t20*t2j*w2hj0 + t20*t20*w2hjj;
                ri = zaehler2/nenner2;
                rj = (t20 - ri*t2i)/t2j;
                norm2 = (laenge2 - zaehler2*zaehler2/nenner2)/(t2j*t2j);
            }

            else
            {
                zaehler2 = t2i*t2i*w2hj0 - t2i*t2j*w2hi0 - t20*t2i*w2hij + t20*t2j*w2hii;
                laenge2 = w2h00*t2i*t2i - 2.0*t20*t2i*w2hi0 + t20*t20*w2hii;
                rj = zaehler2/nenner2;
                ri = (t20 - rj*t2j)/t2i;
                norm2 = (laenge2 - zaehler2*zaehler2/nenner2)/(t2i*t2i);
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
    ConstructLocalMatrix(veci,&localmatrix,tv1,tv2,normr, norml);
    dt = (int) clock() - t;
    time1 += dt;

    if((normr < 1e-25) || (norml < 1e-25))
    {
        if(graph->SavePaList(palist,0,NULL,NULL,NULL,0.0)) return 1;
        FAMGReleaseHeap(FAMG_FROM_TOP);
        return 0;
    }

    t = (int) clock();
    LocalMatMult(&new_localmatrix,&localmatrix);
    dt = (int) clock() - t;
    time2 += dt;


    t = (int) clock();
    LocalScalarProducts(&new_localmatrix,w1,w2);
    dt = (int) clock() - t;
    time3 += dt;


    t = (int) clock();
    GetLocalMinimum1(palist,w1,w2,tv1,tv2,&new_localmatrix);
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
