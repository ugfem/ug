/****************************************************************************/
/*																			*/
/* File:      famg_approx.C													*/
/*																			*/
/* Purpose:   famg graph classes functions									*/
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

#include <iostream.h>
#include <strstream.h>
#include <math.h>

#include "famg_misc.h"
#include "famg_matrix.h"
#include "famg_heap.h"
#include "famg_grid.h"
#include "famg_graph.h"
#include "famg_system.h"

/* RCS_ID
$Header$
*/

struct FAMGSpecialData
{
    int j;
    double hjj, ejj, fj, gj, rj, lj;
};



void FAMGGrid::JJ1(int j, double &hjj, double &ejj)
{
    FAMGMatrixPtr matjz;
    double esum, hsum, mjj, mjz, mzj;

    double omega = FAMGGetParameter()->Getomegar();
    hsum = 0.0;
    esum = 0.0;
    matjz = tmpmatrix->GetStart(j);
    mjj = matjz.GetData();
    while(matjz.GetNext())
    {
        mjz = matjz.GetData();
        mzj = matjz.GetAdjData();
        hsum += mjz*mjz;
        esum += mzj*mzj;
    }

    hjj = (1.0-omega)*(1.0-omega)+omega*omega*hsum/(mjj*mjj);
    ejj = (1.0-omega)*(1.0-omega)+omega*omega*esum/(mjj*mjj);
    
    return;
}
        

void FAMGGrid::FJ1(int j, double &fj, double &gj, double *f, double *g)
{
    FAMGMatrixPtr matjz;
    FAMGNode *node;
    double mjz, mzj, mjj, tf, tg;
    int z, lz;
    
    double omega = FAMGGetParameter()->Getomegar();
    node = graph->GetNode();
    matjz = tmpmatrix->GetStart(j);
    mjj = matjz.GetData();
    lz = (node+j)->GetLocalId();
    fj = (1.0-omega)*f[lz];
    gj = (1.0-omega)*g[lz];
    tf = 0.0;
    tg = 0.0;
    while(matjz.GetNext())
    {
        z = matjz.GetIndex();
        lz = (node+z)->GetLocalId();
        mjz = matjz.GetData();
        mzj = matjz.GetAdjData();
        tf -= mjz*f[lz];
        tg -= mzj*g[lz];
    }
    
    fj += omega*tf/mjj;
    gj += omega*tg/mjj;

    return;
}


int FAMGGrid::LocalJ1(int j, double *h, double *e)
{
    FAMGMatrixPtr matjz;
    FAMGNode *node, *nodej, *nodez;
    double mjj;
    int z, nn;

    double omega = FAMGGetParameter()->Getomegar();
    node = graph->GetNode();
 
    matjz = tmpmatrix->GetStart(j);
    mjj = matjz.GetData();
    nodej = node+j;
    nodej->SetLocalId(0);
    h[0] = 1.0-omega;
    e[0] = 1.0-omega;
    nn = 1;
    while(matjz.GetNext())
    {
        z = matjz.GetIndex();
        nodez = node+z;
        nodez->SetLocalId(nn);
        h[nn] = -omega*matjz.GetData()/mjj;
        e[nn] = -omega*matjz.GetAdjData()/mjj;
        if(nn < 500) nn++;
        else 
        {
            ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
            FAMGWarning(ostr);
        }
    }

    return 0;
}

void FAMGGrid::DLocalJ1(int j)
{
    FAMGMatrixPtr matjz;
    FAMGNode *node, *nodej, *nodez;
    int z;

    node = graph->GetNode();
    nodej = node+j;
    nodej->SetLocalId(-1);
    matjz = tmpmatrix->GetStart(j);
    while(matjz.GetNext())
    {
        z = matjz.GetIndex();
        nodez = node+z;
        nodez->SetLocalId(-1);
    }

   return;
}
   
   
void FAMGGrid::JK1(int k, double &hjk, double &ejk, double *h, double *e)        
{
    FAMGMatrixPtr matkz;
    FAMGNode *node, *nodek, *nodez;
    double mkz, mzk, mkk, th, te;
    int z, lz;
    
    double omega = FAMGGetParameter()->Getomegar();

    node = graph->GetNode();
    nodek = node+k;

    matkz = tmpmatrix->GetStart(k);
    mkk = matkz.GetData();

    hjk = ejk = 0.0;
    lz = nodek->GetLocalId();
    if(lz >= 0)
    {
        hjk += (1.0-omega)*h[lz];
        ejk += (1.0-omega)*e[lz];
    }
    th = 0.0;
    te = 0.0;
    while(matkz.GetNext())
    {
        z = matkz.GetIndex();
        nodez = node+z;
        lz = nodez->GetLocalId();
        if(lz >= 0)
        {
            mkz = matkz.GetData();
            mzk = matkz.GetAdjData();
            th -= mkz*h[lz];
            te -= mzk*e[lz];
        }
    }
    hjk += omega*th/mkk;
    ejk += omega*te/mkk;

    return;
}


void FAMGGrid::FF1(int i, double &ff, double &gg, double *f, double *g)
{
    FAMGMatrixPtr matjz, matij;
    FAMGNode *node, *nodez, *nodej;
    double mjj, mzj, mjz, rj, lj;
    int z, j, lz, nn;

    double omega = FAMGGetParameter()->Getomegar();
    node = graph->GetNode();
    matij = tmpmatrix->GetStart(i);
    nn = 0;
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        matjz = tmpmatrix->GetStart(j);
        mjj = matjz.GetData();
        rj = matij.GetData()/mjj;
        lj = matij.GetAdjData()/mjj;
        
        // most expensive part
        while(matjz.GetNext())
        {
            z = matjz.GetIndex();
            mjz = matjz.GetData();
            mzj = matjz.GetAdjData();
            nodez = node+z;
            lz = nodez->GetLocalId();
            if(lz < 0)
            {
                nodez->SetLocalId(nn);
                f[nn] = -omega*mjz*rj;
                g[nn] = -omega*mzj*lj;
                if(nn < 4000) nn++;
                else 
                {
                    ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
                    FAMGWarning(ostr);
                }
            }
            else
            {
                f[lz] -= omega*mjz*rj;
                g[lz] -= omega*mzj*lj;
            }
        }
    }
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        mjj = tmpmatrix->GetDiag(j);
        rj = matij.GetData();
        lj = matij.GetAdjData();

        nodej = node+j;
        lz = nodej->GetLocalId();
        if(lz < 0)
        {
            nodej->SetLocalId(nn);
            f[nn] = (1.0-omega)*rj;
            g[nn] = (1.0-omega)*lj;
            if(nn < 4000) nn++;
            else 
            {
                ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
                FAMGWarning(ostr);
            }
        }
        else
        {
            f[lz] += (1.0-omega)*rj;
            g[lz] += (1.0-omega)*lj;
        }
    }

    // compute norm
    ff = 0.0; gg = 0.0;
    for(lz = 0; lz < nn; lz++) 
    {
        ff += f[lz]*f[lz];
        gg += g[lz]*g[lz];
    }

    return;
}
        
void FAMGGrid::DFF1(int i)
{
    FAMGMatrixPtr matjz, matij;
    FAMGNode *node, *nodez, *nodej;
    int z, j;

    node = graph->GetNode();
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        nodej = node+j;
        nodej->SetLocalId(-1);
        
        // most expensive part
        matjz = tmpmatrix->GetStart(j);
        while(matjz.GetNext())
        {
            z = matjz.GetIndex();
            nodez = node+z;
            nodez->SetLocalId(-1);
        }
    }

    return;
}
        


double FAMGGrid::BestSecond1(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData *sd, int nnb)
{
    double *tvA, *tvB;
    double aj, ak, bj, bk;
    double hjj, ejj, fj, gj, hkk, ekk, fk, gk, hjk, ejk, h[500], e[500];
    double error, errorA, errorB, min, nenner, tvAj, tvAk, tvBj, tvBk;
    double coeffA[2], coeffB[2];
    int k, np, pa[2], j, z, y;

    double tol = FAMGGetParameter()->Gettol();
    tvA = vector[FAMGTVA];
    tvB = vector[FAMGTVB];

    min = FAMGGetParameter()->Geterror2();
    for(z = 0; z < nnb-1; z++)
    {
        j = sd[z].j;
        FAMGMarkHeap(FAMG_FROM_TOP);
        hjj = (sd[z]).hjj; ejj = sd[z].ejj;
        fj = sd[z].fj; gj = sd[z].gj;
        tvAj = tvA[j];
        tvBj = tvB[j];
            
        if(LocalJ1(j,h,e)) return 1e+20;

        for(y = z+1; y < nnb; y++)
        {
            k = sd[y].j;
            hkk = sd[y].hjj; ekk = sd[y].ejj;
            fk = sd[y].fj; gk = sd[y].gj;
            tvAk = tvA[k];
            tvBk = tvB[k];
 
            JK1(k,hjk,ejk,h,e);

            // compute aj, ak
           if( (Abs(tvAj) < 1e-15) && (Abs(tvAk) < 1e-15))
            {
                if (Abs(rt) > 1e-15) continue;

                 nenner = hjj*hkk-hjk*hjk;
                if (Abs(nenner) < 1e-15) continue;
                ak = (hjj*fk-hjk*fj)/nenner;
                nenner = hjj*hjj+hjk*hjk;
                if (Abs(nenner) < 1e-15) continue;
                aj = (hjj*fj+hjk*fk-(hjj*hjk+hjk*hkk)*ak)/nenner;
            }
            else if(Abs(tvAj) > Abs(tvAk))
            {
                nenner = hjj*tvAk*tvAk + hkk*tvAj*tvAj - 2.0*hjk*tvAj*tvAk;
                if(Abs(nenner) < 1e-15) 
                {
                    aj = rt*tvAj/(tvAj*tvAj + tvAk*tvAk);
                    ak = rt*tvAk/(tvAj*tvAj + tvAk*tvAk);
                }
                else
                {
                    ak = ((hjj*tvAk-hjk*tvAj)*rt + (tvAj*fk-tvAk*fj)*tvAj)/nenner;
                    aj = (rt - ak*tvAk)/tvAj;
                }
            }
            else
            {
                nenner = hkk*tvAj*tvAj + hjj*tvAk*tvAk - 2.0*hjk*tvAk*tvAj;
                if(Abs(nenner) < 1e-15) 
                {
                    aj = rt*tvAj/(tvAj*tvAj + tvAk*tvAk);
                    ak = rt*tvAk/(tvAj*tvAj + tvAk*tvAk);
                }
                else
                {
                    aj = ((hkk*tvAj-hjk*tvAk)*rt + (tvAk*fj-tvAj*fk)*tvAk)/nenner;
                    ak = (rt - aj*tvAj)/tvAk;
                }
            }

            // same for bj, bk

            if( (Abs(tvBj) < 1e-15) && (Abs(tvBk) < 1e-15))
            {
                if (Abs(lt) > 1e-15) continue;

                nenner = ejj*ekk-ejk*ejk;
                if (Abs(nenner) < 1e-15) continue;
                bk = (ejj*gk-ejk*gj)/nenner;
                nenner = ejj*ejj+ejk*ejk;
                if (Abs(nenner) < 1e-15) continue;
                bj = (ejj*gj+ejk*gk-(ejj*ejk+ejk*ekk)*bk)/nenner;
            }
            else if(Abs(tvBj) > Abs(tvBk))
            {
                nenner = ejj*tvBk*tvBk + ekk*tvBj*tvBj - 2.0*ejk*tvBj*tvBk;
                if(Abs(nenner) < 1e-15) 
                {
                    bj = lt*tvBj/(tvBj*tvBj + tvBk*tvBk);
                    bk = lt*tvBk/(tvBj*tvBj + tvBk*tvBk);
                }
                else
                {
                    bk = ((ejj*tvBk-ejk*tvBj)*lt + (tvBj*gk-tvBk*gj)*tvBj)/nenner;
                    bj = (lt - bk*tvBk)/tvBj;
                }
            }
            else
            {
                nenner = ekk*tvBj*tvBj + ejj*tvBk*tvBk - 2.0*ejk*tvBk*tvBj;
                if(Abs(nenner) < 1e-15) 
                {
                    bj = lt*tvBj/(tvBj*tvBj + tvBk*tvBk);
                    bk = lt*tvBk/(tvBj*tvBj + tvBk*tvBk);
                }
                else
                {
                    bj = ((ekk*tvBj-ejk*tvBk)*lt + (tvBk*gj-tvBj*gk)*tvBk)/nenner;
                    bk = (lt - bj*tvBj)/tvBk;
                }
            }

            errorA = ff+aj*aj*hjj+ak*ak*hkk+2.0*ak*aj*hjk-2.0*(aj*fj+ak*fk);
            errorB = gg+bj*bj*ejj+bk*bk*ekk+2.0*bk*bj*ejk-2.0*(bj*gj+bk*gk);
            error = sqrt(errorA*errorB)/(mii*mii);

            if(error*tol < min)
            {
                if (error < min) min = error;
                // k must be first parents !
                np = 2; pa[0] = j; pa[1] = k; 
                coeffA[0] = -aj/mii; coeffA[1] = -ak/mii;
                coeffB[0] = -bj/mii; coeffB[1] = -bk/mii; 
                if(graph->SavePaList(palist,np,pa,coeffA,coeffB,error)) return 1e+20;
            }
        }
        DLocalJ1(j);
        FAMGReleaseHeap(FAMG_FROM_TOP);
    }

    graph->CorrectPaList(palist,min/tol);

    return min;
}


double FAMGGrid::BestFirst1(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData *sd, int nnb)
{
    double *tvA, *tvB, aj, bj, coeffA[1], coeffB[1];
    double errorA, errorB, error, min;
    double hjj, ejj, fj, gj;
    int j, np, pa[1], z;
    
    double tol = FAMGGetParameter()->Gettol();
    tvA = vector[FAMGTVA];
    tvB = vector[FAMGTVB];
 
    min = 1.01*FAMGGetParameter()->Geterror1(); 
    for(z = 0; z < nnb; z++)
    {
        j = sd[z].j;
        aj = rt/tvA[j];
        bj = lt/tvB[j];
            
        hjj = sd[z].hjj; ejj = sd[z].ejj;
        fj = sd[z].fj; gj = sd[z].gj;

        errorA = aj*aj*hjj-2.0*aj*fj+ff;
        errorB = bj*bj*ejj-2.0*bj*gj+gg;
        error = sqrt(errorA*errorB)/(mii*mii);

        if(error*tol < min) // large tolerance
        {
            if (error < min) min = error;
                
            np = 1; pa[0] = j; 
            coeffA[0] = -aj/mii;
            coeffB[0] = -bj/mii; 
            if(graph->SavePaList(palist,np,pa,coeffA,coeffB,error)) return 1e+20;
        }
    }

    graph->CorrectPaList(palist,min/tol); 


    return min;
}

double FAMGGrid::BestFirst5(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData sd)
{
    double *tvA, *tvB, aj, bj, coeffA[1], coeffB[1];
    double errorA, errorB, error, min;
    double hjj, ejj, fj, gj;
    int j, np, pa[1];
    
    double tol = FAMGGetParameter()->Gettol();
    tvA = vector[FAMGTVA];
    tvB = vector[FAMGTVB];
 
    min = 1.01*FAMGGetParameter()->Geterror1(); 
    j = sd.j;
    aj = rt/tvA[j];
    bj = lt/tvB[j];
            
    hjj = sd.hjj; ejj = sd.ejj;
    fj = sd.fj; gj = sd.gj;

    errorA = aj*aj*hjj-2.0*aj*fj+ff;
    errorB = bj*bj*ejj-2.0*bj*gj+gg;
    error = sqrt(errorA*errorB)/(mii*mii);

    if(error*tol < min) 
    {
        np = 1; pa[0] = j; 
        coeffA[0] = -aj/mii;
        coeffB[0] = -bj/mii; 
        if(graph->SavePaList(palist,np,pa,coeffA,coeffB,error)) return 1e+20;
    }

    return error;
}


int FAMGGrid::AnalyseNode3(int i, FAMGPaList *&palist)
{
    FAMGMatrixPtr matij;
    double rj, lj, rt, lt, mii, rmax, lmax, min1, *tvA, *tvB, ff, gg, hjj, ejj, fj, gj, normr, norml, f[4000], g[4000];
    int j, nnb, z;
    FAMGSpecialData *sd;

    tvA = vector[FAMGTVA];
    tvB = vector[FAMGTVB];

    palist = NULL;

    matij = matrix->GetStart(i);
    mii = matij.GetData();
    rt = lt = 0.0;
    normr = norml = 0.0;
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        rj = matij.GetData();
        lj = matij.GetAdjData();
        normr += Abs(rj);
        norml += Abs(lj);
        rt += rj*tvA[j];
        lt += lj*tvB[j];
    }

    normr = normr/Abs(mii);
    norml = norml/Abs(mii);

    if((normr < 1e-15) || (norml < 1e-15))
    {
        if(graph->SavePaList(palist,0,NULL,NULL,NULL,0.0)) return 1;
        return 0;
    }
        
    
    rmax = 0.0; lmax = 0.0; nnb = 0;
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        rj = matij.GetData();
        lj = matij.GetAdjData();
        if (Abs(rj*tvA[j]) > rmax) rmax = Abs(rj*tvA[j]);
        if (Abs(lj*tvB[j]) > lmax) lmax = Abs(lj*tvB[j]);
        nnb++;
    }

    double sigma = FAMGGetParameter()->Getsigma();
    lmax = lmax*sigma; rmax = rmax*sigma;

 
    FAMGMarkHeap(FAMG_FROM_TOP);
    sd = (struct FAMGSpecialData*) FAMGGetMem(nnb*sizeof(struct FAMGSpecialData), FAMG_FROM_TOP);
    if (sd == NULL) return 1;


    FF1(i,ff,gg, f, g);

    z = 0;
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    { 
        j = matij.GetIndex();
        rj = matij.GetData();
        lj = matij.GetAdjData();
        if((Abs(rj*tvA[j]) > rmax) || (Abs(lj*tvB[j]) > lmax))
        {

            JJ1(j,hjj,ejj);
            FJ1(j,fj,gj,f,g);
            sd[z].j = j;
            sd[z].hjj = hjj;
            sd[z].ejj = ejj;
            sd[z].fj = fj;
            sd[z].gj = gj;
            sd[z].rj = rj;
            sd[z].lj = lj;
            z++;
        }
    }
    DFF1(i);


    min1 = BestFirst1(palist,mii,rt,lt,ff,gg,sd,z);
    if(min1 > FAMGGetParameter()->Geterror1()) 
    {
        graph->ClearPaListRev(palist);
        BestSecond1(palist,mii,rt,lt,ff,gg,sd,z);
    } 

    FAMGReleaseHeap(FAMG_FROM_TOP);
    return 0;
}

int FAMGGrid::AnalyseNode4(int i, FAMGPaList *&palist)
{
    FAMGMatrixPtr matij;
    double rj, lj, rt, lt, hr, hl, mii, rmax, lmax, min1, *tvA, *tvB, ff, gg, hjj, ejj, fj, gj, normr, norml, f[4000], g[4000], rmax2, lmax2;
    int j, nnb, z;
    FAMGSpecialData *sd;

    tvA = vector[FAMGTVA];
    tvB = vector[FAMGTVB];
    palist = NULL;

    matij = matrix->GetStart(i);
    mii = matij.GetData();
    rt = lt = 0.0;
    normr = norml = 0.0;
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        rj = matij.GetData();
        lj = matij.GetAdjData();
        normr += Abs(rj);
        norml += Abs(lj);
        rt += rj*tvA[j];
        lt += lj*tvB[j];
    }

    normr = normr/Abs(mii);
    norml = norml/Abs(mii);

    if((normr < 1e-15) || (norml < 1e-15))
    {
        if(graph->SavePaList(palist,0,NULL,NULL,NULL,0.0)) return 1;
        return 0;
    }
                
    rmax2 = rmax = 0.0; lmax2 = lmax = 0.0; nnb = 0;
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        rj = matij.GetData();
        lj = matij.GetAdjData();
        hr = Abs(rj*tvA[j]);
        hl = Abs(lj*tvB[j]);
        if (hr > rmax2) 
        {
            if (hr > rmax) 
            {
                rmax2 = rmax;
                rmax = hr;
            }
            else rmax2 = hr;
        }
        if (hl > lmax2) 
        {
            if (hl > lmax) 
            {
                lmax2 = lmax;
                lmax = hl;
            }
            else lmax2 = hl;
        }
        nnb++;
    }
    
    double sigma = FAMGGetParameter()->Getsigma();
    lmax = lmax2*sigma; rmax = rmax2*sigma;
 
    FAMGMarkHeap(FAMG_FROM_TOP);
    sd = (struct FAMGSpecialData*) FAMGGetMem(nnb*sizeof(struct FAMGSpecialData), FAMG_FROM_TOP);
    if (sd == NULL) return 1;

    FF1(i,ff,gg, f, g);

    z = 0;
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    { 
        j = matij.GetIndex();
        rj = matij.GetData();
        lj = matij.GetAdjData();
        hr = Abs(rj*tvA[j]);
        hl = Abs(lj*tvB[j]);
        if(( hr > rmax) || (hl >  lmax))
        {
            JJ1(j,hjj,ejj);
            FJ1(j,fj,gj,f,g);
            sd[z].j = j;
            sd[z].hjj = hjj;
            sd[z].ejj = ejj;
            sd[z].fj = fj;
            sd[z].gj = gj;
            sd[z].rj = rj;
            sd[z].lj = lj;
            z++;
        }
    }
    DFF1(i);


    min1 = BestFirst1(palist,mii,rt,lt,ff,gg,sd,z);
    if(min1 > FAMGGetParameter()->Geterror1()) 
    {
        graph->ClearPaListRev(palist);
        BestSecond1(palist,mii,rt,lt,ff,gg,sd,z);
    } 

    FAMGReleaseHeap(FAMG_FROM_TOP);
    return 0;
}



int FAMGGrid::AnalyseNode5(int i, FAMGPaList *&palist)
{
    FAMGMatrixPtr matij;
    double rj, lj, rt, lt, mii, min1, min2, *tvA, *tvB, ff, gg, hjj, ejj, fj, gj, normr, norml, f[4000], g[4000], *err1, sigma;
    int j, nnb, z, y, k;
    FAMGSpecialData *sd, *sd2;

    tvA = vector[FAMGTVA];
    tvB = vector[FAMGTVB];
    palist = NULL;

    matij = matrix->GetStart(i);
    mii = matij.GetData();
    rt = lt = 0.0;
    normr = norml = 0.0;
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        rj = matij.GetData();
        lj = matij.GetAdjData();
        normr += Abs(rj);
        norml += Abs(lj);
        rt += rj*tvA[j];
        lt += lj*tvB[j];
    }

    normr = normr/Abs(mii);
    norml = norml/Abs(mii);

    if((normr < 1e-15) || (norml < 1e-15))
    {
        if(graph->SavePaList(palist,0,NULL,NULL,NULL,0.0)) return 1;
        return 0;
    }
        
    nnb = 0;
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    { 
        nnb++;
    }
 
    sigma = FAMGGetParameter()->Getsigma();

    FF1(i,ff,gg, f, g);
   
    FAMGMarkHeap(FAMG_FROM_TOP);
    sd = (struct FAMGSpecialData*) FAMGGetMem(nnb*sizeof(struct FAMGSpecialData), FAMG_FROM_TOP);
    if (sd == NULL) return 1;
    err1= (double*) FAMGGetMem(nnb*sizeof(double), FAMG_FROM_TOP);
    if (err1 == NULL) return 1;


    z = 0;
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    { 
        j = matij.GetIndex();
        JJ1(j,hjj,ejj);
        FJ1(j,fj,gj,f,g);
        sd[z].j = j;
        sd[z].hjj = hjj;
        sd[z].ejj = ejj;
        sd[z].fj = fj;
        sd[z].gj = gj;
        sd[z].rj = rj;
        sd[z].lj = lj;
        z++;
    }
    DFF1(i);


    min1 = min2 = 10000.0*FAMGGetParameter()->Geterror2();
    for(k = 0; k < z; k++)
    {
        err1[k] = BestFirst5(palist,mii,rt,lt,ff,gg,sd[k]);
        if(err1[k] < min2)
        {
            if (err1[k] < min1)
            {
                min2 = min1;
                min1 = err1[k];
            }
            else min2 = err1[k]; 
        }           
    }

    if(min1 < FAMGGetParameter()->Geterror1())
    {
        graph->CorrectPaList(palist,min1/FAMGGetParameter()->Gettol()); 
    }
    else
    {
        sd2 = (struct FAMGSpecialData*) FAMGGetMem(nnb*sizeof(struct FAMGSpecialData), FAMG_FROM_TOP);
        if (sd2 == NULL) return 1;
        y = 0;
        for(k = 0; k < z; k++)
        {
            if(sigma*err1[k] < min2)
            {
                sd2[y] = sd[k];
                y++;
            }
        }
        graph->ClearPaListRev(palist);
        BestSecond1(palist,mii,rt,lt,ff,gg,sd2,y);
    } 

    FAMGReleaseHeap(FAMG_FROM_TOP);
    return 0;
}


void FAMGGrid::JJ0(int j, double &hjj, double &ejj)
{
    FAMGMatrixPtr matjz;
    double esum, hsum, mjj, mjz, mzj, mzz;
    int z;

    const double omegar = FAMGGetParameter()->Getomegar();
    const double omegal = 0.5*FAMGGetParameter()->Getomegal();
    hsum = 0.0;
    esum = 0.0;
    matjz = tmpmatrix->GetStart(j);
    mjj = matjz.GetData();
    while(matjz.GetNext())
    {
        z = matjz.GetIndex();
        mjz = matjz.GetData();
        mzj = matjz.GetAdjData();
        mzz = tmpmatrix->GetDiag(z);
        hsum += mjz*mjz;
        esum += mzj*mzj/(mzz*mzz);
    }

    hjj = (1.0-omegar)*(1.0-omegar) + omegar*omegar*hsum/(mjj*mjj);
    ejj = ((1.0-omegal)*(1.0-omegal) + omegal*omegal*esum)/(mjj*mjj);
    
    return;
}
        

void FAMGGrid::FJ0(int j, double &fj, double &gj, double *f, double *g)
{
    FAMGMatrixPtr matjz;
    FAMGNode *node;
    double mzz, mjz, mzj, mjj, tf, tg;
    int z, lz;
    
    const double omegar = FAMGGetParameter()->Getomegar();
    const double omegal = 0.5*FAMGGetParameter()->Getomegal();
    node = graph->GetNode();
    matjz = tmpmatrix->GetStart(j);
    mjj = matjz.GetData();
    lz = (node+j)->GetLocalId();
    fj = (1.0-omegar)*f[lz];
    gj = (1.0-omegal)*g[lz]/mjj;
    tf = 0.0;
    tg = 0.0;
    while(matjz.GetNext())
    {
        z = matjz.GetIndex();
        lz = (node+z)->GetLocalId();
        mzz = tmpmatrix->GetDiag(z);
        mjz = matjz.GetData();
        mzj = matjz.GetAdjData();
        tf -= mjz*f[lz];
        tg -= mzj*g[lz]/mzz;
    }
    
    fj += omegar*tf/mjj;
    gj += omegal*tg/mjj;

    return;
}


int FAMGGrid::LocalJ0(int j, double *h, double *e)
{
    FAMGMatrixPtr matjz;
    FAMGNode *node, *nodej, *nodez;
    double mjj, mzz;
    int z, nn;

    node = graph->GetNode();

    const double omegar = FAMGGetParameter()->Getomegar();
    const double omegal = 0.5*FAMGGetParameter()->Getomegal();
    matjz = tmpmatrix->GetStart(j);
    mjj = matjz.GetData();
    nodej = node+j;
    nodej->SetLocalId(0);
    h[0] = (1.0-omegar);
    e[0] = (1.0-omegal)/mjj;
    nn = 1;
    while(matjz.GetNext())
    {
        z = matjz.GetIndex();
        mzz = tmpmatrix->GetDiag(z);
        nodez = node+z;
        nodez->SetLocalId(nn);
        h[nn] = -omegar*matjz.GetData()/mjj;
        e[nn] = -omegal*matjz.GetAdjData()/(mzz*mjj);
        if(nn < 500) nn++;
        else 
        {
            ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
            FAMGWarning(ostr);
        }
    }

    return 0;
}

void FAMGGrid::DLocalJ0(int j)
{
    FAMGMatrixPtr matjz;
    FAMGNode *node, *nodej, *nodez;
    int z;

    node = graph->GetNode();

    nodej = node+j;
    nodej->SetLocalId(-1);
    matjz = tmpmatrix->GetStart(j);
    while(matjz.GetNext())
    {
        z = matjz.GetIndex();
        nodez = node+z;
        nodez->SetLocalId(-1);
    }

   return;
}
   
   
void FAMGGrid::JK0(int k, double &hjk, double &ejk, double *h, double *e)        
{
    FAMGMatrixPtr matkz;
    FAMGNode *node, *nodek, *nodez;
    double mkz, mzk, mkk, mzz, th, te;
    int z, lz;
    
    node = graph->GetNode();
    nodek = node+k;

    matkz = tmpmatrix->GetStart(k);
    mkk = matkz.GetData();

    const double omegar = FAMGGetParameter()->Getomegar();
    const double omegal = 0.5*FAMGGetParameter()->Getomegal();
    hjk = ejk = 0.0;
    lz = nodek->GetLocalId();
    if(lz >= 0)
    {
        hjk += (1.0-omegar)*h[lz];
        ejk += (1.0-omegal)*e[lz]/mkk;
    }
    th = 0.0;
    te = 0.0;
    while(matkz.GetNext())
    {
        z = matkz.GetIndex();
        nodez = node+z;
        lz = nodez->GetLocalId();
        if(lz >= 0)
        {
            mkz = matkz.GetData();
            mzk = matkz.GetAdjData();
            mzz = tmpmatrix->GetDiag(z);
            th -= mkz*h[lz];
            te -= mzk*e[lz]/mzz;
        }
    }
    hjk += omegar*th/mkk;
    ejk += omegal*te/mkk;

    return;
}


void FAMGGrid::FF0(int i, double &ff, double &gg, double *f, double *g)
{
    FAMGMatrixPtr matij, matjz;
    FAMGNode *node, *nodez, *nodej;
    double mjj, mzz, mzj, mjz, rj, lj;
    int z, j, lz, nn;

    node = graph->GetNode();
    const double omegar = FAMGGetParameter()->Getomegar();
    const double omegal = 0.5*FAMGGetParameter()->Getomegal();
    nn = 0;
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        matjz = tmpmatrix->GetStart(j);
        mjj = matjz.GetData();
        rj = omegar*matij.GetData()/mjj;
        lj = omegal*matij.GetAdjData()/mjj;
        
        // most expensive part
        while(matjz.GetNext())
        {
            z = matjz.GetIndex();
            mjz = matjz.GetData();
            mzj = matjz.GetAdjData();
            mzz = tmpmatrix->GetDiag(z);
            nodez = node+z;
            lz = nodez->GetLocalId();
            if(lz < 0)
            {
                nodez->SetLocalId(nn);
                f[nn] = -mjz*rj;
                g[nn] = -mzj*lj/mzz;
                if(nn < 4000) nn++;
                else 
                {
                    ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
                    FAMGWarning(ostr);
                }
            }
            else
            {
                f[lz] -= mjz*rj;
                g[lz] -= mzj*lj/mzz;
            }
        }
    }
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        mjj = tmpmatrix->GetDiag(j);
        rj = matij.GetData();
        lj = matij.GetAdjData();

        nodej = node+j;
        lz = nodej->GetLocalId();
        if(lz < 0)
        {
            nodej->SetLocalId(nn);
            f[nn] = (1.0-omegar)*rj;
            g[nn] = (1.0-omegal)*lj/mjj;
            if(nn < 4000) nn++;
            else 
            {
                ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
                FAMGWarning(ostr);
            }
        }
        else
        {
            f[lz] += (1.0-omegar)*rj;
            g[lz] += (1.0-omegal)*lj/mjj;
        }
     }

    // compute norm
    ff = 0.0; gg = 0.0;
    for(lz = 0; lz < nn; lz++) 
    {
        ff += f[lz]*f[lz];
        gg += g[lz]*g[lz];
    }


    return;
}
        
void FAMGGrid::DFF0(int i)
{
    FAMGMatrixPtr matjz, matij;
    FAMGNode *node, *nodez, *nodej;
    int z, j;

    node = graph->GetNode();
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        nodej = node+j;
        nodej->SetLocalId(-1);
        matjz = tmpmatrix->GetStart(j);
        while(matjz.GetNext())
        {
            z = matjz.GetIndex();
            nodez = node+z;
            nodez->SetLocalId(-1);
        }
    }

    return;
}

#ifdef ORIG_TEST
void FAMGGrid::JJ0(int j, double &hjj, double &ejj)
{
    FAMGMatrixPtr matjz;
    double esum, hsum, mjj, mjz, mzj, mzz;
    int z;

    const double omega = FAMGomega();
    const double omega2 = omega*omega;
    hsum = 0.0;
    esum = 0.0;
    matjz = tmpmatrix->GetStart(j);
    mjj = matjz.GetData();
    while(matjz.GetNext())
    {
        z = matjz.GetIndex();
        mjz = matjz.GetData();
        mzj = matjz.GetAdjData();
        mzz = tmpmatrix->GetDiag(z);
        hsum += mjz*mjz;
        esum += mzj*mzj/(mzz*mzz);
    }

    hjj = (1.0-omega)*(1.0-omega) + omega2*hsum/(mjj*mjj);
    ejj = ((omega+omega-omega2)*(omega+omega-omega2) + esum*omega2*omega2)/(mjj*mjj);
    
    return;
}
        

void FAMGGrid::FJ0(int j, double &fj, double &gj, double *f, double *g)
{
    FAMGMatrixPtr matjz;
    FAMGNode *node;
    double mzz, mjz, mzj, mjj, tf, tg;
    int z, lz;
    
    const double omega = FAMGomega();
    const double omega2 = omega*omega;
    node = graph->GetNode();
    matjz = tmpmatrix->GetStart(j);
    mjj = matjz.GetData();
    lz = (node+j)->GetLocalId();
    fj = (1.0-omega)*f[lz];
    gj = (omega+omega-omega2)*g[lz]/mjj;
    tf = 0.0;
    tg = 0.0;
    while(matjz.GetNext())
    {
        z = matjz.GetIndex();
        lz = (node+z)->GetLocalId();
        mzz = tmpmatrix->GetDiag(z);
        mjz = matjz.GetData();
        mzj = matjz.GetAdjData();
        tf -= mjz*f[lz];
        tg -= mzj*g[lz]/mzz;
    }
    
    fj += omega*tf/mjj;
    gj += omega2*tg/mjj;

    return;
}


int FAMGGrid::LocalJ0(int j, double *h, double *e)
{
    FAMGMatrixPtr matjz;
    FAMGNode *node, *nodej, *nodez;
    double mjj, mzz;
    int z, nn;

    node = graph->GetNode();

    const double omega = FAMGomega();
    const double omega2 = omega*omega;
    matjz = tmpmatrix->GetStart(j);
    mjj = matjz.GetData();
    nodej = node+j;
    nodej->SetLocalId(0);
    h[0] = (1.0-omega);
    e[0] = (omega+omega-omega2)/mjj;
    nn = 1;
    while(matjz.GetNext())
    {
        z = matjz.GetIndex();
        mzz = tmpmatrix->GetDiag(z);
        nodez = node+z;
        nodez->SetLocalId(nn);
        h[nn] = -omega*matjz.GetData()/mjj;
        e[nn] = -omega2*matjz.GetAdjData()/(mzz*mjj);
        if(nn < 500) nn++;
        else 
        {
            ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
            FAMGWarning(ostr);
        }
    }

    return 0;
}

void FAMGGrid::DLocalJ0(int j)
{
    FAMGMatrixPtr matjz;
    FAMGNode *node, *nodej, *nodez;
    int z;

    node = graph->GetNode();

    nodej = node+j;
    nodej->SetLocalId(-1);
    matjz = tmpmatrix->GetStart(j);
    while(matjz.GetNext())
    {
        z = matjz.GetIndex();
        nodez = node+z;
        nodez->SetLocalId(-1);
    }

   return;
}
   
   
void FAMGGrid::JK0(int k, double &hjk, double &ejk, double *h, double *e)        
{
    FAMGMatrixPtr matkz;
    FAMGNode *node, *nodek, *nodez;
    double mkz, mzk, mkk, mzz, th, te;
    int z, lz;
    
    node = graph->GetNode();
    nodek = node+k;

    matkz = tmpmatrix->GetStart(k);
    mkk = matkz.GetData();

    const double omega = FAMGomega();
    const double omega2 = omega*omega;
    hjk = ejk = 0.0;
    lz = nodek->GetLocalId();
    if(lz >= 0)
    {
        hjk += (1.0-omega)*h[lz];
        ejk += (omega+omega-omega2)*e[lz]/mkk;
    }
    th = 0.0;
    te = 0.0;
    while(matkz.GetNext())
    {
        z = matkz.GetIndex();
        nodez = node+z;
        lz = nodez->GetLocalId();
        if(lz >= 0)
        {
            mkz = matkz.GetData();
            mzk = matkz.GetAdjData();
            mzz = tmpmatrix->GetDiag(z);
            th -= mkz*h[lz];
            te -= mzk*e[lz]/mzz;
        }
    }
    hjk += omega*th/mkk;
    ejk += omega2*te/mkk;

    return;
}


void FAMGGrid::FF0(int i, double &ff, double &gg, double *f, double *g)
{
    FAMGMatrixPtr matij, matjz;
    FAMGNode *node, *nodez, *nodej;
    double mjj, mzz, mzj, mjz, rj, lj;
    int z, j, lz, nn;

    node = graph->GetNode();
    const double omega = FAMGomega();
    const double omega2 = omega*omega;
    nn = 0;
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        matjz = tmpmatrix->GetStart(j);
        mjj = matjz.GetData();
        rj = omega*matij.GetData()/mjj;
        lj = omega2*matij.GetAdjData()/mjj;
        
        // most expensive part
        while(matjz.GetNext())
        {
            z = matjz.GetIndex();
            mjz = matjz.GetData();
            mzj = matjz.GetAdjData();
            mzz = tmpmatrix->GetDiag(z);
            nodez = node+z;
            lz = nodez->GetLocalId();
            if(lz < 0)
            {
                nodez->SetLocalId(nn);
                f[nn] = -mjz*rj;
                g[nn] = -mzj*lj/mzz;
                if(nn < 4000) nn++;
                else 
                {
                    ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
                    FAMGWarning(ostr);
                }
            }
            else
            {
                f[lz] -= mjz*rj;
                g[lz] -= mzj*lj/mzz;
            }
        }
    }
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        mjj = tmpmatrix->GetDiag(j);
        rj = matij.GetData();
        lj = matij.GetAdjData();

        nodej = node+j;
        lz = nodej->GetLocalId();
        if(lz < 0)
        {
            nodej->SetLocalId(nn);
            f[nn] = (1.0-omega)*rj;
            g[nn] = (omega+omega-omega2)*lj/mjj;
            if(nn < 4000) nn++;
            else 
            {
                ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
                FAMGWarning(ostr);
            }
        }
        else
        {
            f[lz] += (1.0-omega)*rj;
            g[lz] += (omega+omega-omega2)*lj/mjj;
        }
     }

    // compute norm
    ff = 0.0; gg = 0.0;
    for(lz = 0; lz < nn; lz++) 
    {
        ff += f[lz]*f[lz];
        gg += g[lz]*g[lz];
    }


    return;
}
        
void FAMGGrid::DFF0(int i)
{
    FAMGMatrixPtr matjz, matij;
    FAMGNode *node, *nodez, *nodej;
    int z, j;

    node = graph->GetNode();
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        nodej = node+j;
        nodej->SetLocalId(-1);
        matjz = tmpmatrix->GetStart(j);
        while(matjz.GetNext())
        {
            z = matjz.GetIndex();
            nodez = node+z;
            nodez->SetLocalId(-1);
        }
    }

    return;
}
#endif
#ifdef NEW_TEST

void FAMGGrid::JJ0(int j, double &hjj, double &ejj)
{
    FAMGMatrixPtr matjz;
    double esum, hsum, mjj, mjz, mzj, mzz;
    int z;

    const double omega = FAMGomega();
    const double omega2 = omega*omega;
    hsum = 0.0;
    esum = 0.0;
    matjz = tmpmatrix->GetStart(j);
    mjj = matjz.GetData();
    while(matjz.GetNext())
    {
        z = matjz.GetIndex();
        mjz = matjz.GetData();
        mzj = matjz.GetAdjData();
        mzz = tmpmatrix->GetDiag(z);
        hsum += mjz*mjz/(mzz*mzz);
        esum += mzj*mzj/(mzz*mzz);
    }

    hjj = ((omega+omega-omega2)*(omega+omega-omega2) + hsum*omega2*omega2)/(mjj*mjj);
    ejj = ((omega+omega-omega2)*(omega+omega-omega2) + esum*omega2*omega2)/(mjj*mjj);
    
    return;
}
        

void FAMGGrid::FJ0(int j, double &fj, double &gj, double *f, double *g)
{
    FAMGMatrixPtr matjz;
    FAMGNode *node;
    double mzz, mjz, mzj, mjj, tf, tg;
    int z, lz;
    
    const double omega = FAMGomega();
    const double omega2 = omega*omega;
    node = graph->GetNode();
    matjz = tmpmatrix->GetStart(j);
    mjj = matjz.GetData();
    lz = (node+j)->GetLocalId();
    fj = (omega+omega-omega2)*f[lz]/mjj;
    gj = (omega+omega-omega2)*g[lz]/mjj;
    tf = 0.0;
    tg = 0.0;
    while(matjz.GetNext())
    {
        z = matjz.GetIndex();
        lz = (node+z)->GetLocalId();
        mzz = tmpmatrix->GetDiag(z);
        mjz = matjz.GetData();
        mzj = matjz.GetAdjData();
        tf -= mjz*f[lz]/mzz;
        tg -= mzj*g[lz]/mzz;
    }
    
    fj += omega2*tf/mjj;
    gj += omega2*tg/mjj;

    return;
}


int FAMGGrid::LocalJ0(int j, double *h, double *e)
{
    FAMGMatrixPtr matjz;
    FAMGNode *node, *nodej, *nodez;
    double mjj, mzz;
    int z, nn;

    node = graph->GetNode();

    const double omega = FAMGomega();
    const double omega2 = omega*omega;
    matjz = tmpmatrix->GetStart(j);
    mjj = matjz.GetData();
    nodej = node+j;
    nodej->SetLocalId(0);
    h[0] = (omega+omega-omega2)/mjj;
    e[0] = (omega+omega-omega2)/mjj;
    nn = 1;
    while(matjz.GetNext())
    {
        z = matjz.GetIndex();
        mzz = tmpmatrix->GetDiag(z);
        nodez = node+z;
        nodez->SetLocalId(nn);
        h[nn] = -omega2*matjz.GetData()/(mzz*mjj);
        e[nn] = -omega2*matjz.GetAdjData()/(mzz*mjj);
        if(nn < 500) nn++;
        else 
        {
            ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
            FAMGWarning(ostr);
        }
    }

    return 0;
}

void FAMGGrid::DLocalJ0(int j)
{
    FAMGMatrixPtr matjz;
    FAMGNode *node, *nodej, *nodez;
    int z;

    node = graph->GetNode();

    nodej = node+j;
    nodej->SetLocalId(-1);
    matjz = tmpmatrix->GetStart(j);
    while(matjz.GetNext())
    {
        z = matjz.GetIndex();
        nodez = node+z;
        nodez->SetLocalId(-1);
    }

   return;
}
   
   
void FAMGGrid::JK0(int k, double &hjk, double &ejk, double *h, double *e)        
{
    FAMGMatrixPtr matkz;
    FAMGNode *node, *nodek, *nodez;
    double mkz, mzk, mkk, mzz, th, te;
    int z, lz;
    
    node = graph->GetNode();
    nodek = node+k;

    matkz = tmpmatrix->GetStart(k);
    mkk = matkz.GetData();

    const double omega = FAMGomega();
    const double omega2 = omega*omega;
    hjk = ejk = 0.0;
    lz = nodek->GetLocalId();
    if(lz >= 0)
    {
        hjk += (omega+omega-omega2)*h[lz]/mkk;
        ejk += (omega+omega-omega2)*e[lz]/mkk;
    }
    th = 0.0;
    te = 0.0;
    while(matkz.GetNext())
    {
        z = matkz.GetIndex();
        nodez = node+z;
        lz = nodez->GetLocalId();
        if(lz >= 0)
        {
            mkz = matkz.GetData();
            mzk = matkz.GetAdjData();
            mzz = tmpmatrix->GetDiag(z);
            th -= mkz*h[lz]/mzz;
            te -= mzk*e[lz]/mzz;
        }
    }
    hjk += omega2*th/mkk;
    ejk += omega2*te/mkk;

    return;
}


void FAMGGrid::FF0(int i, double &ff, double &gg, double *f, double *g)
{
    FAMGMatrixPtr matij, matjz;
    FAMGNode *node, *nodez, *nodej;
    double mjj, mzz, mzj, mjz, rj, lj;
    int z, j, lz, nn;

    node = graph->GetNode();
    const double omega = FAMGomega();
    const double omega2 = omega*omega;
    nn = 0;
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        matjz = tmpmatrix->GetStart(j);
        mjj = matjz.GetData();
        rj = omega2*matij.GetData()/mjj;
        lj = omega2*matij.GetAdjData()/mjj;
        
        // most expensive part
        while(matjz.GetNext())
        {
            z = matjz.GetIndex();
            mjz = matjz.GetData();
            mzj = matjz.GetAdjData();
            mzz = tmpmatrix->GetDiag(z);
            nodez = node+z;
            lz = nodez->GetLocalId();
            if(lz < 0)
            {
                nodez->SetLocalId(nn);
                f[nn] = -mjz*rj/mzz;
                g[nn] = -mzj*lj/mzz;
                if(nn < 4000) nn++;
                else 
                {
                    ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
                    FAMGWarning(ostr);
                }
            }
            else
            {
                f[lz] -= mjz*rj/mzz;
                g[lz] -= mzj*lj/mzz;
            }
        }
    }
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        mjj = tmpmatrix->GetDiag(j);
        rj = matij.GetData();
        lj = matij.GetAdjData();

        nodej = node+j;
        lz = nodej->GetLocalId();
        if(lz < 0)
        {
            nodej->SetLocalId(nn);
            f[nn] = (omega+omega-omega2)*rj/mjj;
            g[nn] = (omega+omega-omega2)*lj/mjj;
            if(nn < 4000) nn++;
            else 
            {
                ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
                FAMGWarning(ostr);
            }
       }
        else
        {
            f[lz] += (omega+omega-omega2)*rj/mjj;
            g[lz] += (omega+omega-omega2)*lj/mjj;
        }
     }

    // compute norm
    ff = 0.0; gg = 0.0;
    for(lz = 0; lz < nn; lz++) 
    {
        ff += f[lz]*f[lz];
        gg += g[lz]*g[lz];
    }


    return;
}
        
void FAMGGrid::DFF0(int i)
{
    FAMGMatrixPtr matjz, matij;
    FAMGNode *node, *nodez, *nodej;
    int z, j;

    node = graph->GetNode();
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        nodej = node+j;
        nodej->SetLocalId(-1);
        matjz = tmpmatrix->GetStart(j);
        while(matjz.GetNext())
        {
            z = matjz.GetIndex();
            nodez = node+z;
            nodez->SetLocalId(-1);
        }
    }

    return;
}
#endif

double FAMGGrid::BestSecond0(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData *sd, int nnb)
{
    double *tvA, *tvB;
    double aj, ak, bj, bk;
    double hjj, ejj, fj, gj, hkk, ekk, fk, gk, hjk, ejk, h[500], e[500];
    double error, errorA, errorB, min, nenner, tvAj, tvAk, tvBj, tvBk;
    double coeffA[2], coeffB[2];
    int k, np, pa[2], j, z, y;

    double tol = FAMGGetParameter()->Gettol();
    tvA = vector[FAMGTVA];
    tvB = vector[FAMGTVB];
 
    min = FAMGGetParameter()->Geterror2();
    for(z = 0; z < nnb; z++)
    {
        j = sd[z].j;
        FAMGMarkHeap(FAMG_FROM_TOP);
        hjj = (sd[z]).hjj; ejj = sd[z].ejj;
        fj = sd[z].fj; gj = sd[z].gj;
        if(LocalJ0(j,h,e)) return 1e+20;
        tvAj = tvA[j];
        tvBj = tvB[j];

        for(y = z+1; y < nnb; y++)
        {
            k = sd[y].j;
            hkk = sd[y].hjj; ekk = sd[y].ejj;
            fk = sd[y].fj; gk = sd[y].gj;
            tvAk = tvA[k];
            tvBk = tvB[k];
 
            JK0(k,hjk,ejk,h,e);

            // compute aj, ak

            if( (Abs(tvAj) < 1e-15) && (Abs(tvAk) < 1e-15))
            {
                if(Abs(rt) > 1e-15*Abs(mii)) continue;

                nenner = hjj*hkk-hjk*hjk;
                if (Abs(nenner) < 1e-15) continue;
                ak = (hjj*fk-hjk*fj)/nenner;
                aj = (hkk*fj-hjk*fk)/nenner;
            }
            else if(Abs(tvAj) > Abs(tvAk))
            {
                nenner = hjj*tvAk*tvAk + hkk*tvAj*tvAj - 2.0*hjk*tvAj*tvAk;
                if(Abs(nenner) < 1e-15) 
                {
                    aj = rt*tvAj/(tvAj*tvAj + tvAk*tvAk);
                    ak = rt*tvAk/(tvAj*tvAj + tvAk*tvAk);
                }
                else
                {
                    ak = ((hjj*tvAk-hjk*tvAj)*rt + (tvAj*fk-tvAk*fj)*tvAj)/nenner;
                    aj = (rt - ak*tvAk)/tvAj;
                }
            }
            else
            {
                nenner = hkk*tvAj*tvAj + hjj*tvAk*tvAk - 2.0*hjk*tvAk*tvAj;
                if(Abs(nenner) < 1e-15) 
                {
                    aj = rt*tvAj/(tvAj*tvAj + tvAk*tvAk);
                    ak = rt*tvAk/(tvAj*tvAj + tvAk*tvAk);
                }
                else
                {
                    aj = ((hkk*tvAj-hjk*tvAk)*rt + (tvAk*fj-tvAj*fk)*tvAk)/nenner;
                    ak = (rt - aj*tvAj)/tvAk;
                }
            }

            // same for bj, bk

            if( (Abs(tvBj) < 1e-15) && (Abs(tvBk) < 1e-15))
            {
                if(Abs(lt) > 1e-15*Abs(mii)) continue;

                nenner = ejj*ekk-ejk*ejk;
                if (Abs(nenner) < 1e-15) continue;
                bk = (ejj*gk-ejk*gj)/nenner;
                bj = (ekk*gj-ejk*gk)/nenner;
            }
            else if(Abs(tvBj) > Abs(tvBk))
            {
                nenner = ejj*tvBk*tvBk + ekk*tvBj*tvBj - 2.0*ejk*tvBj*tvBk;
                if(Abs(nenner) < 1e-15) 
                {
                    bj = lt*tvBj/(tvBj*tvBj + tvBk*tvBk);
                    bk = lt*tvBk/(tvBj*tvBj + tvBk*tvBk);
                }
                else
                {
                    bk = ((ejj*tvBk-ejk*tvBj)*lt + (tvBj*gk-tvBk*gj)*tvBj)/nenner;
                    bj = (lt - bk*tvBk)/tvBj;
                }
            }
            else
            {
                nenner = ekk*tvBj*tvBj + ejj*tvBk*tvBk - 2.0*ejk*tvBk*tvBj;
                if(Abs(nenner) < 1e-15) 
                {
                    bj = lt*tvBj/(tvBj*tvBj + tvBk*tvBk);
                    bk = lt*tvBk/(tvBj*tvBj + tvBk*tvBk);
                }
                else
                {
                    bj = ((ekk*tvBj-ejk*tvBk)*lt + (tvBk*gj-tvBj*gk)*tvBk)/nenner;
                    bk = (lt - bj*tvBj)/tvBk;
                }
            }

            errorA = ff+aj*aj*hjj+ak*ak*hkk+2.0*ak*aj*hjk-2.0*(aj*fj+ak*fk);
            errorB = gg+bj*bj*ejj+bk*bk*ekk+2.0*bk*bj*ejk-2.0*(bj*gj+bk*gk);
            errorA = sqrt(errorA)/Abs(mii);
            errorB = sqrt(errorB);
            error = errorB*errorA;
            

            if(error*tol < min)
            {
                if (error < min) min = error;
                // k must be first parents !
                np = 2; pa[0] = j; pa[1] = k; 
                coeffA[0] = -aj/mii; coeffA[1] = -ak/mii;
                coeffB[0] = -bj/mii; coeffB[1] = -bk/mii; 
                if(graph->SavePaList(palist,np,pa,coeffA,coeffB,error)) return 1e+20;
            }
        }
        DLocalJ0(j);
        FAMGReleaseHeap(FAMG_FROM_TOP);
    }

    graph->CorrectPaList(palist,min/tol);

    return min;
}

double FAMGGrid::BestFirst0(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData *sd, int nnb)
{
    double *tvA, *tvB, aj, bj, coeffA[1], coeffB[1];
    double errorA, errorB, error, min;
    double hjj, ejj, fj, gj;
    int j, np, pa[1], z;
    
    double tol = FAMGGetParameter()->Gettol();
    tvA = vector[FAMGTVA];
    tvB = vector[FAMGTVB];

    min = 1.01*FAMGGetParameter()->Geterror1(); 
    for(z = 0; z < nnb; z++)
    {
        j = sd[z].j;
        aj = rt/tvA[j];
        bj = lt/tvB[j];

        hjj = sd[z].hjj; ejj = sd[z].ejj;
        fj = sd[z].fj; gj = sd[z].gj;

        errorA = aj*aj*hjj-2.0*aj*fj+ff;
        errorB = bj*bj*ejj-2.0*bj*gj+gg;
        errorA = sqrt(errorA)/Abs(mii);
        errorB = sqrt(errorB);
        error = errorB*errorA;

        if(error*tol < min) // large tolerance
        {
            if (error < min) min = error;
                
            np = 1; pa[0] = j; 
            coeffA[0] = -aj/mii;
            coeffB[0] = -bj/mii; 
            if(graph->SavePaList(palist,np,pa,coeffA,coeffB,error)) return 1e+20;
        }
    }

    graph->CorrectPaList(palist,min/tol); 


    return min;
}
double FAMGGrid::BestFirst2(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData sd)
{
    double *tvA, *tvB, aj, bj, coeffA[1], coeffB[1];
    double errorA, errorB, error, min;
    double hjj, ejj, fj, gj;
    int j, np, pa[1];
    
    double tol = FAMGGetParameter()->Gettol();
    tvA = vector[FAMGTVA];
    tvB = vector[FAMGTVB];
 
    min = 1.01*FAMGGetParameter()->Geterror1(); 
    j = sd.j;
    aj = rt/tvA[j];
    bj = lt/tvB[j];
            
    hjj = sd.hjj; ejj = sd.ejj;
    fj = sd.fj; gj = sd.gj;

    errorA = aj*aj*hjj-2.0*aj*fj+ff;
    errorB = bj*bj*ejj-2.0*bj*gj+gg;
    errorA = sqrt(errorA)/Abs(mii);
    errorB = sqrt(errorB);
   
    error = errorB*errorA;

    if(error*tol < min) 
    {
        np = 1; pa[0] = j; 
        coeffA[0] = -aj/mii;
        coeffB[0] = -bj/mii; 
        if(graph->SavePaList(palist,np,pa,coeffA,coeffB,error))return 1e+20;
    }

    return error;
}

int FAMGGrid::AnalyseNode0(int i, FAMGPaList *&palist)
{
    FAMGMatrixPtr matij;
    double rj, lj, rt, lt, mjj, mii, rmax, lmax, min1, *tvA, *tvB, ff, gg, hjj, ejj, fj, gj, f[4000], g[4000], normr, norml;
    int j, nnb, z;
    FAMGSpecialData *sd;

    tvA = vector[FAMGTVA];
    tvB = vector[FAMGTVB];

    palist = NULL;

    rt = lt = 0.0;
    normr = norml = 0.0;
    matij = matrix->GetStart(i);
    mii = matij.GetData();
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        rj = matij.GetData();;
        lj = matij.GetAdjData();;
        mjj = matrix->GetDiag(j);
        normr += Abs(rj/mjj);
        norml += Abs(lj);
        rt += rj*tvA[j];
        lt += lj*tvB[j];
    }
        
    normr = normr/Abs(mii);

    if((normr < 1e-15) || (norml < 1e-15))
    {
        if(graph->SavePaList(palist,0,NULL,NULL,NULL,0.0)) return 1;
        return 0;
    }

    rmax = 0.0; lmax = 0.0; nnb = 0;
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        rj = matij.GetData();;
        lj = matij.GetAdjData();;
        mjj = tmpmatrix->GetDiag(j);
        if (Abs(rj*tvA[j]) > rmax) rmax = Abs(rj*tvA[j]);
        if (Abs(lj*tvB[j]) > lmax) lmax = Abs(lj*tvB[j]);
        nnb++;
    }
                
    double sigma = FAMGGetParameter()->Getsigma();
    lmax = lmax*sigma; rmax = rmax*sigma;

    FAMGMarkHeap(FAMG_FROM_TOP);
    sd = (struct FAMGSpecialData*) FAMGGetMem(nnb*sizeof(struct FAMGSpecialData), FAMG_FROM_TOP);
    if (sd == NULL) return 1;

    FF0(i,ff,gg, f, g);

    z = 0;
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        rj = matij.GetData();;
        lj = matij.GetAdjData();;
        mjj = tmpmatrix->GetDiag(j);
        if((Abs(rj*tvA[j]) > rmax) || (Abs(lj*tvB[j]) > lmax))
        {

            JJ0(j,hjj,ejj);
            FJ0(j,fj,gj, f, g);
            sd[z].j = j;
            sd[z].hjj = hjj;
            sd[z].ejj = ejj;
            sd[z].fj = fj;
            sd[z].gj = gj;
            sd[z].rj = rj;
            sd[z].lj = lj/mjj;
            z++;
        }
    }
    DFF0(i);


    min1 = BestFirst0(palist,mii,rt,lt,ff,gg,sd,z);
    if(min1 > FAMGGetParameter()->Geterror1()) 
    {
        graph->ClearPaListRev(palist);
        BestSecond0(palist,mii,rt,lt,ff,gg,sd,z);
    } 

    FAMGReleaseHeap(FAMG_FROM_TOP);
    return 0;
}


int FAMGGrid::AnalyseNode1(int i, FAMGPaList *&palist)
{
    FAMGMatrixPtr matij;
    double rj, lj, rt, lt, mii, rmax, lmax, min1, *tvA, *tvB, ff, gg, hjj, ejj, fj, gj, normr, norml, f[4000], g[4000], hr, hl, rmax2, lmax2, mjj;
    int j, nnb, z;
    FAMGSpecialData *sd;

    tvA = vector[FAMGTVA];
    tvB = vector[FAMGTVB];
    palist = NULL;

    matij = matrix->GetStart(i);
    mii = matij.GetData();
    rt = lt = 0.0;
    normr = norml = 0.0;
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        rj = matij.GetData();
        lj = matij.GetAdjData();
        mjj = matrix->GetDiag(j);
        normr += Abs(rj);
        norml += Abs(lj/mjj);
        rt += rj*tvA[j];
        lt += lj*tvB[j];
    }

    normr = normr/Abs(mii);

    if((normr < 1e-15) || (norml < 1e-15))
    {
        if(graph->SavePaList(palist,0,NULL,NULL,NULL,0.0)) return 1;
        return 0;
    }


    rmax = rmax2 = 0.0; lmax = lmax2 = 0.0; nnb = 0;
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        rj = matij.GetData();
        lj = matij.GetAdjData();
        mjj = tmpmatrix->GetDiag(j);
        hr = Abs(rj*tvA[j]);
        hl = Abs(lj*tvB[j]);
        if(hr > rmax2)
        {
            if (hr > rmax)
            {
                rmax2 = rmax;
                rmax = hr;
            }
            else rmax2 = hr;
        }
        if(hl > lmax2)
        {
            if (hl > lmax)
            {
                lmax2 = lmax;
                lmax = hl;
            }
            else lmax2 = hl;
        }
        nnb++;
    }
                
    double sigma = FAMGGetParameter()->Getsigma();
    lmax = lmax2*sigma; rmax = rmax2*sigma;

    FAMGMarkHeap(FAMG_FROM_TOP);
    sd = (struct FAMGSpecialData*) FAMGGetMem(nnb*sizeof(struct FAMGSpecialData), FAMG_FROM_TOP);
    if (sd == NULL) return 1;

    FF0(i,ff,gg, f, g);

    z = 0;
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    { 
        j = matij.GetIndex();
        rj = matij.GetData();
        lj = matij.GetAdjData();
        mjj = tmpmatrix->GetDiag(j);
        hr = Abs(rj*tvA[j]);
        hl = Abs(lj*tvB[j]);
        if((hr > rmax) || (hl > lmax))
        {
            JJ0(j,hjj,ejj);
            FJ0(j,fj,gj, f, g);
            sd[z].j = j;
            sd[z].hjj = hjj;
            sd[z].ejj = ejj;
            sd[z].fj = fj;
            sd[z].gj = gj;
            sd[z].rj = rj;
            sd[z].lj = lj/mjj;
            z++;
        }
    }
    DFF0(i);


    min1 = BestFirst0(palist,mii,rt,lt,ff,gg,sd,z);
    if(min1 > FAMGGetParameter()->Geterror1()) 
    {
        graph->ClearPaListRev(palist);
        BestSecond0(palist,mii,rt,lt,ff,gg,sd,z);
    } 

    FAMGReleaseHeap(FAMG_FROM_TOP);
    return 0;
}



int FAMGGrid::AnalyseNode2(int i, FAMGPaList *&palist)
{
    FAMGMatrixPtr matij;
    double rj, lj, rt, lt, mii, min1, min2, *tvA, *tvB, ff, gg, hjj, ejj, fj, gj, normr, norml, f[4000], g[4000], *err1, sigma, mjj;
    int j, nnb, z, y, k;
    FAMGSpecialData *sd, *sd2;

    tvA = vector[FAMGTVA];
    tvB = vector[FAMGTVB];
    palist = NULL;

    matij = matrix->GetStart(i);
    mii = matij.GetData();
    rt = lt = 0.0;
    normr = norml = 0.0;
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        rj = matij.GetData();
        lj = matij.GetAdjData();
        mjj = matrix->GetDiag(j);
        normr += Abs(rj);
        norml += Abs(lj/mjj);
        rt += rj*tvA[j];
        lt += lj*tvB[j];
    }

    normr = normr/Abs(mii);

    if((normr < 1e-15) || (norml < 1e-15))
    {
        if(graph->SavePaList(palist,0,NULL,NULL,NULL,0.0)) return 1;
        return 0;
    }
        
    nnb = 0;
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    { 
        nnb++;
    }

    FF0(i,ff,gg, f, g);
    
    sigma = FAMGGetParameter()->Getsigma();
 
    FAMGMarkHeap(FAMG_FROM_TOP);
    sd = (struct FAMGSpecialData*) FAMGGetMem(nnb*sizeof(struct FAMGSpecialData), FAMG_FROM_TOP);
    if (sd == NULL) return 1;

    err1= (double*) FAMGGetMem(nnb*sizeof(double), FAMG_FROM_TOP);
    if (err1 == NULL) return 1;

    z = 0;
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    { 
        j = matij.GetIndex();
        rj = matij.GetData();
        lj = matij.GetAdjData();
        mjj = tmpmatrix->GetDiag(j);
        JJ0(j,hjj,ejj);
        FJ0(j,fj,gj,f,g);
        sd[z].j = j;
        sd[z].hjj = hjj;
        sd[z].ejj = ejj;
        sd[z].fj = fj;
        sd[z].gj = gj;
        sd[z].rj = rj;
        sd[z].lj = lj/mjj;
        z++;
    }
    DFF0(i);


    min1 = min2 = 10000.0*FAMGGetParameter()->Geterror2();
    for(k = 0; k < z; k++)
    {
        err1[k] = BestFirst2(palist,mii,rt,lt,ff,gg,sd[k]);
        if(err1[k] < min2)
        {
            if (err1[k] < min1)
            {
                min2 = min1;
                min1 = err1[k];
            }
            else min2 = err1[k]; 
        }           
    }

    if(min1 < FAMGGetParameter()->Geterror1())
    {
        graph->CorrectPaList(palist,min1/FAMGGetParameter()->Gettol()); 
    }
    else
    {
        sd2 = (struct FAMGSpecialData*) FAMGGetMem(nnb*sizeof(struct FAMGSpecialData), FAMG_FROM_TOP);
        if (sd2 == NULL) return 1;
        y = 0;
        for(k = 0; k < z; k++)
        {
            if(sigma*err1[k] < min2)
            {
                sd2[y] = sd[k];
                y++;
            }
        }
        graph->ClearPaListRev(palist);
        BestSecond0(palist,mii,rt,lt,ff,gg,sd2,y);
    } 

    FAMGReleaseHeap(FAMG_FROM_TOP);
    return 0;
}

