/****************************************************************************/
/*																			*/
/* File:      famg_approx.C													*/
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
/* Remarks:																	*/
/*																			*/
/****************************************************************************/

/*     
This is a pretty old version of the famg parents selection process !
A new version will be checked in soon !
*/


#include <iostream.h>
#include <strstream.h>
#include <math.h>

#include "famg_algebra.h"
#include "famg_misc.h"
#include "famg_heap.h"
#include "famg_grid.h"
#include "famg_graph.h"
#include "famg_system.h"

/* RCS_ID
$Header$
*/

#ifndef FAMG_SPARSE_BLOCK

struct FAMGSpecialData
{
    int j;
    double hjj, ejj, fj, gj, rj, lj;
};



void FAMGGrid::JJ1(const FAMGVectorEntry &vecj, double &hjj, double &ejj)
{
	const FAMGMatrixAlg& M = *GetTmpMatrix();
    FAMGMatrixEntry matjz;
    double esum, hsum, mjj, mjz, mzj;
	
    double omega = FAMGGetParameter()->Getomegar();
    hsum = 0.0;
    esum = 0.0;
	FAMGMatrixIter miter(M,vecj);
	miter(matjz);   // get diagonal
    mjj = M[matjz];
    while(miter(matjz))
    {
        mjz = M[matjz];
        mzj = M.GetAdjData(matjz);
        hsum += mjz*mjz;
        esum += mzj*mzj;
    }

    hjj = (1.0-omega)*(1.0-omega)+omega*omega*hsum/(mjj*mjj);
    ejj = (1.0-omega)*(1.0-omega)+omega*omega*esum/(mjj*mjj);
    
    return;
}
        

void FAMGGrid::FJ1(int j, const FAMGVectorEntry &vecj, double &fj, double &gj, double *f, double *g)
{
	const FAMGMatrixAlg& M = *GetTmpMatrix();
    FAMGMatrixEntry matjz;
    const FAMGGraph &graph = *GetGraph();
    double mjz, mzj, mjj, tf, tg;
    int lz;
    
    double omega = FAMGGetParameter()->Getomegar();
	FAMGMatrixIter miter(M,vecj);
	miter(matjz);   // get diagonal
    mjj = M[matjz];
    lz = graph.GetNode(j)->GetLocalId();
    fj = (1.0-omega)*f[lz];
    gj = (1.0-omega)*g[lz];
    tf = 0.0;
    tg = 0.0;
    while(miter(matjz))
    {
        lz = graph.GetNode(matjz.dest().GetIndex())->GetLocalId();
        mjz = M[matjz];
        mzj = M.GetAdjData(matjz);
        tf -= mjz*f[lz];
        tg -= mzj*g[lz];
    }
    
    fj += omega*tf/mjj;
    gj += omega*tg/mjj;

    return;
}


int FAMGGrid::LocalJ1(int j, const FAMGVectorEntry &vecj, double *h, double *e)
{
	const FAMGMatrixAlg& M = *GetTmpMatrix();
    FAMGMatrixEntry matjz;
    const FAMGGraph &graph = *GetGraph();
    double mjj;
    int nn;

    double omega = FAMGGetParameter()->Getomegar();
 
	graph.GetNode(j)->SetLocalId(0);
    h[0] = 1.0-omega;
    e[0] = 1.0-omega;
    nn = 1;
	FAMGMatrixIter miter(M,vecj);
	miter(matjz);   // get diagonal
    mjj = M[matjz];
    while(miter(matjz))
    {
		graph.GetNode(matjz.dest().GetIndex())->SetLocalId(nn);
        h[nn] = -omega*M[matjz]/mjj;
        e[nn] = -omega*M.GetAdjData(matjz)/mjj;
        if(nn < 500) 
			nn++;
        else 
        {
            ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
            FAMGWarning(ostr);
        }
    }

    return 0;
}

void FAMGGrid::DLocalJ1(int j, const FAMGVectorEntry &vecj)
{
	const FAMGMatrixAlg& M = *GetTmpMatrix();
    FAMGMatrixEntry matjz;
    const FAMGGraph &graph = *GetGraph();

	graph.GetNode(j)->SetLocalId(-1);
	FAMGMatrixIter miter(M,vecj);
	miter(matjz);   // skip diagonal
    while(miter(matjz))
		graph.GetNode(matjz.dest().GetIndex())->SetLocalId(-1);

   return;
}
   
   
void FAMGGrid::JK1(int k, const FAMGVectorEntry &veck, double &hjk, double &ejk, double *h, double *e)        
{
	const FAMGMatrixAlg& M = *GetTmpMatrix();
    FAMGMatrixEntry matkz;
    const FAMGGraph &graph = *GetGraph();
    double mkz, mzk, mkk, th, te;
    int lz;
    
    double omega = FAMGGetParameter()->Getomegar();
    hjk = ejk = 0.0;
    lz = graph.GetNode(k)->GetLocalId();
    if(lz >= 0)
    {
        hjk += (1.0-omega)*h[lz];
        ejk += (1.0-omega)*e[lz];
    }
    th = 0.0;
    te = 0.0;

	FAMGMatrixIter miter(M,veck);
	miter(matkz);    // get diagonal
    mkk = M[matkz];
    while(miter(matkz))
    {
        lz = graph.GetNode(matkz.dest().GetIndex())->GetLocalId();
        if(lz >= 0)
        {
            mkz = M[matkz];
            mzk = M.GetAdjData(matkz);
            th -= mkz*h[lz];
            te -= mzk*e[lz];
        }
    }
    hjk += omega*th/mkk;
    ejk += omega*te/mkk;

    return;
}


void FAMGGrid::FF1(const FAMGVectorEntry &veci, double &ff, double &gg, double *f, double *g)
{
	FAMGVectorEntry vecj;
	const FAMGMatrixAlg& M = *GetTmpMatrix();
    FAMGMatrixEntry matij, matjz;
    const FAMGGraph &graph = *GetGraph();
	FAMGNode *nodej, *nodez;
    double mjj, mzj, mjz, rj, lj;
    int lz, nn;

    double omega = FAMGGetParameter()->Getomegar();
    nn = 0;
	FAMGMatrixIter miter(M,veci);
	miter(matij);   // skip diagonal
    while(miter(matij))
    {
        FAMGMatrixIter mjziter(M,matij.dest());
        mjziter(matjz);    // get diagonal
        mjj = M[matjz];
        rj = M[matij]/mjj;
        lj = M.GetAdjData(matij)/mjj;
        
        // most expensive part
        while(mjziter(matjz))
        {
            mjz = M[matjz];
            mzj = M.GetAdjData(matjz);
            nodez = graph.GetNode(matjz.dest().GetIndex());
            lz = nodez->GetLocalId();
            if(lz < 0)
            {
                nodez->SetLocalId(nn);
                f[nn] = -omega*mjz*rj;
                g[nn] = -omega*mzj*lj;
                if(nn < 4000)
					nn++;
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
    miter.reset();
	miter(matij);    // skip diagonal
    while(miter(matij))
    {
		vecj = matij.dest();
        mjj = M.DiagValue(vecj);
        rj = M[matij];
        lj = M.GetAdjData(matij);

		nodej = graph.GetNode(vecj.GetIndex());
        lz = nodej->GetLocalId();
        if(lz < 0)
        {
            nodej->SetLocalId(nn);
            f[nn] = (1.0-omega)*rj;
            g[nn] = (1.0-omega)*lj;
            if(nn < 4000)
				nn++;
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
        
void FAMGGrid::DFF1(const FAMGVectorEntry &veci)
{
	FAMGVectorEntry vecj;
	const FAMGMatrixAlg& M = *GetTmpMatrix();
    FAMGMatrixEntry matij, matjz;
    const FAMGGraph &graph = *GetGraph();

	FAMGMatrixIter miter(M,veci);
	miter(matij);   // skip diagonal
    while(miter(matij))
    {
		vecj = matij.dest();
		graph.GetNode(vecj.GetIndex())->SetLocalId(-1);
        
        // most expensive part
		FAMGMatrixIter mjziter(M,matij.dest());
		mjziter(matjz); // skip diagonal
        while(mjziter(matjz))
			graph.GetNode(matjz.dest().GetIndex())->SetLocalId(-1);
    }

    return;
}
        

double FAMGGrid::BestSecond1(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData *sd, int nnb)
{
    double aj, ak, bj, bk;
    double hjj, ejj, fj, gj, hkk, ekk, fk, gk, hjk, ejk, h[500], e[500];
    double error, errorA, errorB, min, nenner, tvAj, tvAk, tvBj, tvBk;
    double coeffA[2], coeffB[2];
    int k, np, pa[2], j, z, y;
    FAMGGraph &graph = *GetGraph();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
	FAMGVectorEntry vecj, veck;

    double tol = FAMGGetParameter()->Gettol();

    min = FAMGGetParameter()->Geterror2();
    for(z = 0; z < nnb-1; z++)
    {
        j = sd[z].j;
		vecj = graph.GetNode(j)->GetVec();

        FAMGMarkHeap(FAMG_FROM_TOP);
        hjj = (sd[z]).hjj; ejj = sd[z].ejj;
        fj = sd[z].fj; gj = sd[z].gj;
        tvAj = tvA[vecj];
        tvBj = tvB[vecj];
            
        if(LocalJ1(j,vecj,h,e)) 
			return 1e+20;

        for(y = z+1; y < nnb; y++)
        {
            k = sd[y].j;
			veck = graph.GetNode(k)->GetVec();

            hkk = sd[y].hjj; ekk = sd[y].ejj;
            fk = sd[y].fj; gk = sd[y].gj;
            tvAk = tvA[veck];
            tvBk = tvB[veck];
 
            JK1(k,veck,hjk,ejk,h,e);

            // compute aj, ak
			if( (Abs(tvAj) < 1e-15) && (Abs(tvAk) < 1e-15))
            {
                if (Abs(rt) > 1e-15) 
					continue;
				nenner = hjj*hkk-hjk*hjk;
                if (Abs(nenner) < 1e-15) 
					continue;
                ak = (hjj*fk-hjk*fj)/nenner;
                nenner = hjj*hjj+hjk*hjk;
                if (Abs(nenner) < 1e-15) 
					continue;
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
                if (Abs(lt) > 1e-15)
					continue;
                nenner = ejj*ekk-ejk*ejk;
                if (Abs(nenner) < 1e-15)
					continue;
                bk = (ejj*gk-ejk*gj)/nenner;
                nenner = ejj*ejj+ejk*ejk;
                if (Abs(nenner) < 1e-15)
					continue;
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
                if (error < min)
					min = error;
                // k must be first parents !
                np = 2; pa[0] = j; pa[1] = k; 
                coeffA[0] = -aj/mii; coeffA[1] = -ak/mii;
                coeffB[0] = -bj/mii; coeffB[1] = -bk/mii; 
                if(graph.SavePaList(palist,np,pa,coeffA,coeffB,error))
					return 1e+20;
            }
        }
        DLocalJ1(j,vecj);
        FAMGReleaseHeap(FAMG_FROM_TOP);
    }

    graph.CorrectPaList(palist,min/tol);

    return min;
}


double FAMGGrid::BestFirst1(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData *sd, int nnb)
{
    double aj, bj, coeffA[1], coeffB[1];
    double errorA, errorB, error, min;
    double hjj, ejj, fj, gj;
    int j, np, pa[1], z;
	FAMGGraph &graph = *GetGraph();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
	FAMGVectorEntry vecj;
   
    double tol = FAMGGetParameter()->Gettol();
    min = 1.01*FAMGGetParameter()->Geterror1(); 
    for(z = 0; z < nnb; z++)
    {
        j = sd[z].j;
		vecj = graph.GetNode(j)->GetVec();
        aj = rt/tvA[vecj];
        bj = lt/tvB[vecj];
            
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
            if(graph.SavePaList(palist,np,pa,coeffA,coeffB,error)) 
				return 1e+20;
        }
    }

    graph.CorrectPaList(palist,min/tol); 


    return min;
}

double FAMGGrid::BestFirst5(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData sd)
{
    double aj, bj, coeffA[1], coeffB[1];
    double errorA, errorB, error, min;
    double hjj, ejj, fj, gj;
    int j, np, pa[1];
	FAMGGraph &graph = *GetGraph();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
	FAMGVectorEntry vecj;
    
    double tol = FAMGGetParameter()->Gettol();
    min = 1.01*FAMGGetParameter()->Geterror1(); 
    j = sd.j;
	vecj = graph.GetNode(j)->GetVec();
    aj = rt/tvA[vecj];
    bj = lt/tvB[vecj];
            
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
        if(graph.SavePaList(palist,np,pa,coeffA,coeffB,error)) 
			return 1e+20;
    }

    return error;
}


int FAMGGrid::AnalyseNode3(const FAMGVectorEntry &veci, FAMGPaList *&palist)
{
    double rj, lj, rt, lt, mii, rmax, lmax, min1, ff, gg, hjj, ejj, fj, gj, normr, norml, f[4000], g[4000];
    int j, nnb, z;
    FAMGSpecialData *sd;
	FAMGGraph &graph = *GetGraph();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
	const FAMGMatrixAlg &M = *GetConsMatrix();
	FAMGVectorEntry vecj;
	FAMGMatrixEntry matij;

    palist = NULL;

    rt = lt = 0.0;
    normr = norml = 0.0;
    FAMGMatrixIter miter(M,veci);
	miter(matij);
    mii = M[matij];
    while(miter(matij))
    {
		vecj = matij.dest();
        rj = M[matij];
        lj = M.GetAdjData(matij);
        normr += Abs(rj);
        norml += Abs(lj);
        rt += rj*tvA[vecj];
        lt += lj*tvB[vecj];
    }

    normr = normr/Abs(mii);
    norml = norml/Abs(mii);

    if((normr < 1e-15) || (norml < 1e-15))
    {
        if(graph.SavePaList(palist,0,NULL,NULL,NULL,0.0)) return 1;
        return 0;
    }
        
    
    rmax = 0.0; lmax = 0.0; nnb = 0;
	const FAMGMatrixAlg &MTemp = *GetTmpMatrix();
	FAMGMatrixIter mtiter(MTemp,veci);
    mtiter(matij);    // skip diagonal
    while(mtiter(matij))
    {
        vecj = matij.dest();
        rj = MTemp[matij];
        lj = MTemp.GetAdjData(matij);
        if (Abs(rj*tvA[vecj]) > rmax) 
			rmax = Abs(rj*tvA[vecj]);
        if (Abs(lj*tvB[vecj]) > lmax)
			lmax = Abs(lj*tvB[vecj]);
        nnb++;
    }

    double sigma = FAMGGetParameter()->Getsigma();
    lmax = lmax*sigma; rmax = rmax*sigma;

 
    FAMGMarkHeap(FAMG_FROM_TOP);
    sd = (struct FAMGSpecialData*) FAMGGetMem(nnb*sizeof(struct FAMGSpecialData), FAMG_FROM_TOP);
    if (sd == NULL) 
		return 1;

    FF1(veci,ff,gg, f, g);

    z = 0;
	mtiter.reset();
	mtiter(matij);// skip diagonal
    while(mtiter(matij))
    { 
        vecj = matij.dest();
		j = vecj.GetIndex();
        rj = MTemp[matij];
		lj = MTemp.GetAdjData(matij);
        if((Abs(rj*tvA[vecj]) > rmax) || (Abs(lj*tvB[vecj]) > lmax))
        {

            JJ1(vecj,hjj,ejj);
            FJ1(j,vecj,fj,gj,f,g);
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
    DFF1(veci);


    min1 = BestFirst1(palist,mii,rt,lt,ff,gg,sd,z);
    if(min1 > FAMGGetParameter()->Geterror1()) 
    {
        graph.ClearPaListRev(palist);
        BestSecond1(palist,mii,rt,lt,ff,gg,sd,z);
    } 

    FAMGReleaseHeap(FAMG_FROM_TOP);
    return 0;
}

int FAMGGrid::AnalyseNode4(const FAMGVectorEntry &veci, FAMGPaList *&palist)
{
    double rj, lj, rt, lt, hr, hl, mii, rmax, lmax, min1, ff, gg, hjj, ejj, fj, gj, normr, norml, f[4000], g[4000], rmax2, lmax2;
    int j, nnb, z;
    FAMGSpecialData *sd;
	FAMGGraph &graph = *GetGraph();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
	const FAMGMatrixAlg &M = *GetConsMatrix();
	FAMGVectorEntry vecj;
	FAMGMatrixEntry matij;

    palist = NULL;

    rt = lt = 0.0;
    normr = norml = 0.0;
    FAMGMatrixIter miter(M,veci);
	miter(matij);
    mii = M[matij];
    while(miter(matij))
    {
        vecj = matij.dest();
        rj = M[matij];
        lj = M.GetAdjData(matij);
        normr += Abs(rj);
        norml += Abs(lj);
        rt += rj*tvA[vecj];
        lt += lj*tvB[vecj];
    }

    normr = normr/Abs(mii);
    norml = norml/Abs(mii);

    if((normr < 1e-15) || (norml < 1e-15))
    {
        if(graph.SavePaList(palist,0,NULL,NULL,NULL,0.0)) 
			return 1;
        return 0;
    }
                
    rmax2 = rmax = 0.0; lmax2 = lmax = 0.0; nnb = 0;
	const FAMGMatrixAlg &MTemp = *GetTmpMatrix();
	FAMGMatrixIter mtiter(MTemp,veci);
    mtiter(matij);    // skip diagonal
    while(mtiter(matij))
    {
        vecj = matij.dest();
        rj = MTemp[matij];
        lj = MTemp.GetAdjData(matij);
        hr = Abs(rj*tvA[vecj]);
        hl = Abs(lj*tvB[vecj]);
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
    if (sd == NULL) 
		return 1;

    FF1(veci,ff,gg, f, g);

    z = 0;
	mtiter.reset();
	mtiter(matij);// skip diagonal
    while(mtiter(matij))
    { 
        vecj = matij.dest();
		j = vecj.GetIndex();
        rj = MTemp[matij];
		lj = MTemp.GetAdjData(matij);
        hr = Abs(rj*tvA[vecj]);
        hl = Abs(lj*tvB[vecj]);
        if(( hr > rmax) || (hl >  lmax))
        {
            JJ1(vecj,hjj,ejj);
            FJ1(j,vecj,fj,gj,f,g);
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
    DFF1(veci);


    min1 = BestFirst1(palist,mii,rt,lt,ff,gg,sd,z);
    if(min1 > FAMGGetParameter()->Geterror1()) 
    {
        graph.ClearPaListRev(palist);
        BestSecond1(palist,mii,rt,lt,ff,gg,sd,z);
    } 

    FAMGReleaseHeap(FAMG_FROM_TOP);
    return 0;
}



int FAMGGrid::AnalyseNode5(const FAMGVectorEntry &veci, FAMGPaList *&palist)
{
    double rj, lj, rt, lt, mii, min1, min2, ff, gg, hjj, ejj, fj, gj, normr, norml, f[4000], g[4000], *err1, sigma;
    int nnb, j, z, y, k;
    FAMGSpecialData *sd, *sd2;
	FAMGGraph &graph = *GetGraph();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
	const FAMGMatrixAlg &M = *GetConsMatrix();
	FAMGVectorEntry vecj;
	FAMGMatrixEntry matij;

    palist = NULL;

    rt = lt = 0.0;
    normr = norml = 0.0;
    FAMGMatrixIter miter(M,veci);
	miter(matij);
    mii = M[matij];
    while(miter(matij))
    {
        vecj = matij.dest();
        rj = M[matij];
        lj = M.GetAdjData(matij);
        normr += Abs(rj);
        norml += Abs(lj);
        rt += rj*tvA[vecj];
        lt += lj*tvB[vecj];
    }

    normr = normr/Abs(mii);
    norml = norml/Abs(mii);

    if((normr < 1e-15) || (norml < 1e-15))
    {
        if(graph.SavePaList(palist,0,NULL,NULL,NULL,0.0))
			return 1;
        return 0;
    }
        
    nnb = 0;
	const FAMGMatrixAlg &MTemp = *GetTmpMatrix();
	FAMGMatrixIter mtiter(MTemp,veci);
    mtiter(matij);    // skip diagonal
    while(mtiter(matij))
        nnb++;
 
    sigma = FAMGGetParameter()->Getsigma();

    FF1(veci,ff,gg, f, g);
   
    FAMGMarkHeap(FAMG_FROM_TOP);
    sd = (struct FAMGSpecialData*) FAMGGetMem(nnb*sizeof(struct FAMGSpecialData), FAMG_FROM_TOP);
    if (sd == NULL)
		return 1;
    err1= (double*) FAMGGetMem(nnb*sizeof(double), FAMG_FROM_TOP);
    if (err1 == NULL)
		return 1;

    z = 0;
	mtiter.reset();
	mtiter(matij);// skip diagonal
    while(mtiter(matij))
    { 
        vecj = matij.dest();
		j = vecj.GetIndex();
        JJ1(vecj,hjj,ejj);
        FJ1(j,vecj,fj,gj,f,g);
        sd[z].j = j;
        sd[z].hjj = hjj;
        sd[z].ejj = ejj;
        sd[z].fj = fj;
        sd[z].gj = gj;
        sd[z].rj = rj;
        sd[z].lj = lj;
        z++;
    }
    DFF1(veci);


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
            else 
				min2 = err1[k]; 
        }           
    }

    if(min1 < FAMGGetParameter()->Geterror1())
    {
        graph.CorrectPaList(palist,min1/FAMGGetParameter()->Gettol()); 
    }
    else
    {
        sd2 = (struct FAMGSpecialData*) FAMGGetMem(nnb*sizeof(struct FAMGSpecialData), FAMG_FROM_TOP);
        if (sd2 == NULL)
			return 1;
        y = 0;
        for(k = 0; k < z; k++)
        {
            if(sigma*err1[k] < min2)
            {
                sd2[y] = sd[k];
                y++;
            }
        }
        graph.ClearPaListRev(palist);
        BestSecond1(palist,mii,rt,lt,ff,gg,sd2,y);
    } 

    FAMGReleaseHeap(FAMG_FROM_TOP);
    return 0;
}


void FAMGGrid::JJ0(const FAMGVectorEntry &vecj, double &hjj, double &ejj)
{
	const FAMGMatrixAlg& M = *GetTmpMatrix();
    FAMGMatrixEntry matjz;
    double esum, hsum, mjj, mjz, mzj, mzz;

    const double omegar = FAMGGetParameter()->Getomegar();
    const double omegal = 0.5*FAMGGetParameter()->Getomegal();
    hsum = 0.0;
    esum = 0.0;
	FAMGMatrixIter miter(M,vecj);
	miter(matjz);   // get diagonal
    mjj = M[matjz];
    while(miter(matjz))
    {
        mjz = M[matjz];
        mzj = M.GetAdjData(matjz);
        mzz = M.DiagValue(matjz.dest());
        hsum += mjz*mjz;
        esum += mzj*mzj/(mzz*mzz);
    }

    hjj = (1.0-omegar)*(1.0-omegar) + omegar*omegar*hsum/(mjj*mjj);
    ejj = ((1.0-omegal)*(1.0-omegal) + omegal*omegal*esum)/(mjj*mjj);
    
    return;
}
        

void FAMGGrid::FJ0(int j, const FAMGVectorEntry &vecj, double &fj, double &gj, double *f, double *g)
{
	const FAMGMatrixAlg& M = *GetTmpMatrix();
    FAMGMatrixEntry matjz;
    const FAMGGraph &graph = *GetGraph();
    double mzz, mjz, mzj, mjj, tf, tg;
    int lz;
    
    const double omegar = FAMGGetParameter()->Getomegar();
    const double omegal = 0.5*FAMGGetParameter()->Getomegal();
	FAMGMatrixIter miter(M,vecj);
	miter(matjz);   // get diagonal
    mjj = M[matjz];
    lz = graph.GetNode(j)->GetLocalId();
    fj = (1.0-omegar)*f[lz];
    gj = (1.0-omegal)*g[lz]/mjj;
    tf = 0.0;
    tg = 0.0;
    while(miter(matjz))
    {
		lz = graph.GetNode(matjz.dest().GetIndex())->GetLocalId();
        mzz = M.DiagValue(matjz.dest());
        mjz = M[matjz];
        mzj = M.GetAdjData(matjz);
        tf -= mjz*f[lz];
        tg -= mzj*g[lz]/mzz;
    }
    
    fj += omegar*tf/mjj;
    gj += omegal*tg/mjj;

    return;
}


int FAMGGrid::LocalJ0(int j, const FAMGVectorEntry &vecj, double *h, double *e)
{
	const FAMGMatrixAlg& M = *GetTmpMatrix();
    FAMGMatrixEntry matjz;
    const FAMGGraph &graph = *GetGraph();
    double mjj, mzz;
    int nn;

    const double omegar = FAMGGetParameter()->Getomegar();
    const double omegal = 0.5*FAMGGetParameter()->Getomegal();
	graph.GetNode(j)->SetLocalId(0);
	FAMGMatrixIter miter(M,vecj);
	miter(matjz);   // get diagonal
    mjj = M[matjz];
    h[0] = (1.0-omegar);
    e[0] = (1.0-omegal)/mjj;
    nn = 1;
    while(miter(matjz))
    {
		graph.GetNode(matjz.dest().GetIndex())->SetLocalId(nn);
        mzz = M.DiagValue(matjz.dest());
        h[nn] = -omegar*M[matjz]/mjj;
        e[nn] = -omegal*M.GetAdjData(matjz)/(mzz*mjj);
        if(nn < 500) nn++;
        else 
        {
            ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
            FAMGWarning(ostr);
        }
    }

    return 0;
}

void FAMGGrid::DLocalJ0(int j, const FAMGVectorEntry &vecj)
{
	const FAMGMatrixAlg& M = *GetTmpMatrix();
    FAMGMatrixEntry matjz;
    const FAMGGraph &graph = *GetGraph();

	graph.GetNode(j)->SetLocalId(-1);
	FAMGMatrixIter miter(M,vecj);
	miter(matjz);   // skip diagonal
    while(miter(matjz))
		graph.GetNode(matjz.dest().GetIndex())->SetLocalId(-1);

   return;
}


void FAMGGrid::JK0(int k, const FAMGVectorEntry &veck, double &hjk, double &ejk, double *h, double *e)        
{
	const FAMGMatrixAlg& M = *GetTmpMatrix();
    FAMGMatrixEntry matkz;
    const FAMGGraph &graph = *GetGraph();
    double mkz, mzk, mkk, mzz, th, te;
    int lz;

    const double omegar = FAMGGetParameter()->Getomegar();
    const double omegal = 0.5*FAMGGetParameter()->Getomegal();
    hjk = ejk = 0.0;
    lz = graph.GetNode(k)->GetLocalId();

	FAMGMatrixIter miter(M,veck);
	miter(matkz);    // get diagonal
    mkk = M[matkz];

    if(lz >= 0)
    {
        hjk += (1.0-omegar)*h[lz];
        ejk += (1.0-omegal)*e[lz]/mkk;
    }
    th = 0.0;
    te = 0.0;
    while(miter(matkz))
    {
        lz = graph.GetNode(matkz.dest().GetIndex())->GetLocalId();
        if(lz >= 0)
        {
            mkz = M[matkz];
            mzk = M.GetAdjData(matkz);
            mzz = M.DiagValue(matkz.dest());
            th -= mkz*h[lz];
            te -= mzk*e[lz]/mzz;
        }
    }
    hjk += omegar*th/mkk;
    ejk += omegal*te/mkk;

    return;
}


void FAMGGrid::FF0(const FAMGVectorEntry &veci, double &ff, double &gg, double *f, double *g)
{
	FAMGVectorEntry vecj, vecz;
	const FAMGMatrixAlg& M = *GetTmpMatrix();
    FAMGMatrixEntry matij, matjz;
    const FAMGGraph &graph = *GetGraph();
    FAMGNode *nodez, *nodej;
    double mjj, mzz, mzj, mjz, rj, lj;
    int lz, nn;

    const double omegar = FAMGGetParameter()->Getomegar();
    const double omegal = 0.5*FAMGGetParameter()->Getomegal();
    nn = 0;
	FAMGMatrixIter miter(M,veci);
	miter(matij);   // skip diagonal
    while(miter(matij))
    {
        FAMGMatrixIter mjziter(M,matij.dest());
        mjziter(matjz);    // get diagonal
        mjj = M[matjz];
        rj = omegar*M[matij]/mjj;
        lj = omegal*M.GetAdjData(matij)/mjj;
        
        // most expensive part
        while(mjziter(matjz))
        {
            mjz = M[matjz];
            mzj = M.GetAdjData(matjz);
			vecz = matjz.dest();
            mzz = M.DiagValue(vecz);
            nodez = graph.GetNode(vecz.GetIndex());
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
    miter.reset();
	miter(matij);    // skip diagonal
    while(miter(matij))
    {
		vecj = matij.dest();
        mjj = M.DiagValue(vecj);
        rj = M[matij];
        lj = M.GetAdjData(matij);

        nodej = graph.GetNode(vecj.GetIndex());
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
        
void FAMGGrid::DFF0(const FAMGVectorEntry &veci)
{
	FAMGVectorEntry vecj;
	const FAMGMatrixAlg& M = *GetTmpMatrix();
    FAMGMatrixEntry matij, matjz;
    const FAMGGraph &graph = *GetGraph();

	FAMGMatrixIter miter(M,veci);
	miter(matij);   // skip diagonal
    while(miter(matij))
    {
		vecj = matij.dest();
		graph.GetNode(vecj.GetIndex())->SetLocalId(-1);

        
        // most expensive part
		FAMGMatrixIter mjziter(M,matij.dest());
		mjziter(matjz); // skip diagonal
        while(mjziter(matjz))
			graph.GetNode(matjz.dest().GetIndex())->SetLocalId(-1);
    }

    return;
}


double FAMGGrid::BestSecond0(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData *sd, int nnb)
{
    double aj, ak, bj, bk;
    double hjj, ejj, fj, gj, hkk, ekk, fk, gk, hjk, ejk, h[500], e[500];
    double error, errorA, errorB, min, nenner, tvAj, tvAk, tvBj, tvBk;
    double coeffA[2], coeffB[2];
    int k, np, pa[2], j, z, y;
    FAMGGraph &graph = *GetGraph();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
	FAMGVectorEntry vecj, veck;

    double tol = FAMGGetParameter()->Gettol();
 
    min = FAMGGetParameter()->Geterror2();
    for(z = 0; z < nnb; z++)
    {
        j = sd[z].j;
		vecj = graph.GetNode(j)->GetVec();

        FAMGMarkHeap(FAMG_FROM_TOP);
        hjj = (sd[z]).hjj; ejj = sd[z].ejj;
        fj = sd[z].fj; gj = sd[z].gj;
        if(LocalJ0(j,vecj,h,e))
			return 1e+20;
        tvAj = tvA[vecj];
        tvBj = tvB[vecj];

        for(y = z+1; y < nnb; y++)
        {
            k = sd[y].j;
			veck = graph.GetNode(k)->GetVec();

            hkk = sd[y].hjj; ekk = sd[y].ejj;
            fk = sd[y].fj; gk = sd[y].gj;
            tvAk = tvA[veck];
            tvBk = tvB[veck];
 
            JK0(k,veck,hjk,ejk,h,e);

            // compute aj, ak

            if( (Abs(tvAj) < 1e-15) && (Abs(tvAk) < 1e-15))
            {
                if(Abs(rt) > 1e-15*Abs(mii))
					continue;

                nenner = hjj*hkk-hjk*hjk;
                if (Abs(nenner) < 1e-15)
					continue;
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
                if(Abs(lt) > 1e-15*Abs(mii))
					continue;

                nenner = ejj*ekk-ejk*ejk;
                if (Abs(nenner) < 1e-15)
					continue;
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
                if(graph.SavePaList(palist,np,pa,coeffA,coeffB,error))
					return 1e+20;
            }
        }
        DLocalJ0(j,vecj);
        FAMGReleaseHeap(FAMG_FROM_TOP);
    }

    graph.CorrectPaList(palist,min/tol);

    return min;
}

double FAMGGrid::BestFirst0(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData *sd, int nnb)
{
    double aj, bj, coeffA[1], coeffB[1];
    double errorA, errorB, error, min;
    double hjj, ejj, fj, gj;
    int j, np, pa[1], z;
	FAMGGraph &graph = *GetGraph();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
	FAMGVectorEntry vecj;
  
    double tol = FAMGGetParameter()->Gettol();
    min = 1.01*FAMGGetParameter()->Geterror1(); 
    for(z = 0; z < nnb; z++)
    {
        j = sd[z].j;
		vecj = graph.GetNode(j)->GetVec();

		aj = rt/tvA[vecj];
        bj = lt/tvB[vecj];

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
            if(graph.SavePaList(palist,np,pa,coeffA,coeffB,error))
				return 1e+20;
        }
    }

    graph.CorrectPaList(palist,min/tol); 


    return min;
}


double FAMGGrid::BestFirst2(FAMGPaList *&palist, double mii, double rt, double lt, double ff, double gg, FAMGSpecialData sd)
{
    double aj, bj, coeffA[1], coeffB[1];
    double errorA, errorB, error, min;
    double hjj, ejj, fj, gj;
    int j, np, pa[1];
 	FAMGGraph &graph = *GetGraph();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
	FAMGVectorEntry vecj;
   
    double tol = FAMGGetParameter()->Gettol(); 
    min = 1.01*FAMGGetParameter()->Geterror1(); 
    j = sd.j;
    aj = rt/tvA[vecj];
    bj = lt/tvB[vecj];
            
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
        if(graph.SavePaList(palist,np,pa,coeffA,coeffB,error))
			return 1e+20;
    }

    return error;
}

int FAMGGrid::AnalyseNode0(const FAMGVectorEntry &veci, FAMGPaList *&palist)
{
    double rj, lj, rt, lt, mjj, mii, rmax, lmax, min1, ff, gg, hjj, ejj, fj, gj, f[4000], g[4000], normr, norml;
    int j, nnb, z;
    FAMGSpecialData *sd;
	FAMGGraph &graph = *GetGraph();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
	const FAMGMatrixAlg &M = *GetConsMatrix();
	FAMGVectorEntry vecj;
	FAMGMatrixEntry matij;

    palist = NULL;

    rt = lt = 0.0;
    normr = norml = 0.0;
    FAMGMatrixIter miter(M,veci);
	miter(matij);
    mii = M[matij];
    while(miter(matij))
    {
		vecj = matij.dest();
        rj = M[matij];
        lj = M.GetAdjData(matij);
        mjj = M.DiagValue(matij.dest());
        normr += Abs(rj/mjj);
        norml += Abs(lj);
        rt += rj*tvA[vecj];
        lt += lj*tvB[vecj];
    }

    normr = normr/Abs(mii);

    if((normr < 1e-15) || (norml < 1e-15))
    {
        if(graph.SavePaList(palist,0,NULL,NULL,NULL,0.0))
			return 1;
        return 0;
    }

    rmax = 0.0; lmax = 0.0; nnb = 0;
	const FAMGMatrixAlg &MTemp = *GetTmpMatrix();
	FAMGMatrixIter mtiter(MTemp,veci);
    mtiter(matij);    // skip diagonal
    while(mtiter(matij))
    {
        vecj = matij.dest();
        rj = MTemp[matij];
        lj = MTemp.GetAdjData(matij);
        mjj = M.DiagValue(matij.dest());
        if (Abs(rj*tvA[vecj]) > rmax)
			rmax = Abs(rj*tvA[vecj]);
        if (Abs(lj*tvB[vecj]) > lmax)
			lmax = Abs(lj*tvB[vecj]);
        nnb++;
    }
          
    double sigma = FAMGGetParameter()->Getsigma();
    lmax = lmax*sigma; rmax = rmax*sigma;

    FAMGMarkHeap(FAMG_FROM_TOP);
    sd = (struct FAMGSpecialData*) FAMGGetMem(nnb*sizeof(struct FAMGSpecialData), FAMG_FROM_TOP);

	if (sd == NULL)
		return 1;

    FF0(veci,ff,gg, f, g);

    z = 0;
	mtiter.reset();
	mtiter(matij);// skip diagonal
    while(mtiter(matij))
    { 
        vecj = matij.dest();
		j = vecj.GetIndex();
        rj = MTemp[matij];
		lj = MTemp.GetAdjData(matij);
        mjj = M.DiagValue(matij.dest());
        if((Abs(rj*tvA[vecj]) > rmax) || (Abs(lj*tvB[vecj]) > lmax))
        {
            JJ0(vecj,hjj,ejj);
            FJ0(j,vecj,fj,gj, f, g);
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
    DFF0(veci);


    min1 = BestFirst0(palist,mii,rt,lt,ff,gg,sd,z);
    if(min1 > FAMGGetParameter()->Geterror1()) 
    {
        graph.ClearPaListRev(palist);
        BestSecond0(palist,mii,rt,lt,ff,gg,sd,z);
    } 

    FAMGReleaseHeap(FAMG_FROM_TOP);
    return 0;
}


int FAMGGrid::AnalyseNode1(const FAMGVectorEntry &veci, FAMGPaList *&palist)
{
    double rj, lj, rt, lt, mii, rmax, lmax, min1, ff, gg, hjj, ejj, fj, gj, normr, norml, f[4000], g[4000], hr, hl, rmax2, lmax2, mjj;
    int j, nnb, z;
    FAMGSpecialData *sd;
	FAMGGraph &graph = *GetGraph();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
	const FAMGMatrixAlg &M = *GetConsMatrix();
	FAMGVectorEntry vecj;
	FAMGMatrixEntry matij;

    palist = NULL;

    rt = lt = 0.0;
    normr = norml = 0.0;
    FAMGMatrixIter miter(M,veci);
	miter(matij);
    mii = M[matij];
     while(miter(matij))
    {
        vecj = matij.dest();
        rj = M[matij];
        lj = M.GetAdjData(matij);
        mjj = M.DiagValue(matij.dest());
        normr += Abs(rj);
        norml += Abs(lj/mjj);
        rt += rj*tvA[vecj];
        lt += lj*tvB[vecj];
    }

    normr = normr/Abs(mii);

    if((normr < 1e-15) || (norml < 1e-15))
    {
        if(graph.SavePaList(palist,0,NULL,NULL,NULL,0.0))
			return 1;
        return 0;
    }


    rmax = rmax2 = 0.0; lmax = lmax2 = 0.0; nnb = 0;
	const FAMGMatrixAlg &MTemp = *GetTmpMatrix();
	FAMGMatrixIter mtiter(MTemp,veci);
    mtiter(matij);    // skip diagonal
    while(mtiter(matij))
    {
        vecj = matij.dest();
        rj = MTemp[matij];
        lj = MTemp.GetAdjData(matij);
        mjj = M.DiagValue(matij.dest());
        hr = Abs(rj*tvA[vecj]);
        hl = Abs(lj*tvB[vecj]);
        if(hr > rmax2)
        {
            if (hr > rmax)
            {
                rmax2 = rmax;
                rmax = hr;
            }
            else 
				rmax2 = hr;
        }
        if(hl > lmax2)
        {
            if (hl > lmax)
            {
                lmax2 = lmax;
                lmax = hl;
            }
            else 
				lmax2 = hl;
        }
        nnb++;
    }
                
    double sigma = FAMGGetParameter()->Getsigma();
    lmax = lmax2*sigma; rmax = rmax2*sigma;

    FAMGMarkHeap(FAMG_FROM_TOP);
    sd = (struct FAMGSpecialData*) FAMGGetMem(nnb*sizeof(struct FAMGSpecialData), FAMG_FROM_TOP);
    if (sd == NULL)
		return 1;

    FF0(veci,ff,gg, f, g);

    z = 0;
	mtiter.reset();
	mtiter(matij);// skip diagonal
    while(mtiter(matij))
    { 
        vecj = matij.dest();
		j = vecj.GetIndex();
        rj = MTemp[matij];
		lj = MTemp.GetAdjData(matij);
        mjj = M.DiagValue(matij.dest());
        hr = Abs(rj*tvA[vecj]);
        hl = Abs(lj*tvB[vecj]);
        if((hr > rmax) || (hl > lmax))
        {
            JJ0(vecj,hjj,ejj);
            FJ0(j,vecj,fj,gj, f, g);
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
    DFF0(veci);


    min1 = BestFirst0(palist,mii,rt,lt,ff,gg,sd,z);
    if(min1 > FAMGGetParameter()->Geterror1()) 
    {
        graph.ClearPaListRev(palist);
        BestSecond0(palist,mii,rt,lt,ff,gg,sd,z);
    } 

    FAMGReleaseHeap(FAMG_FROM_TOP);
    return 0;
}


int FAMGGrid::AnalyseNode2(const FAMGVectorEntry &veci, FAMGPaList *&palist)
{
    double rj, lj, rt, lt, mii, min1, min2, ff, gg, hjj, ejj, fj, gj, normr, norml, f[4000], g[4000], *err1, sigma, mjj;
    int j, nnb, z, y, k;
    FAMGSpecialData *sd, *sd2;
	FAMGGraph &graph = *GetGraph();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
	const FAMGMatrixAlg &M = *GetConsMatrix();
	FAMGVectorEntry vecj;
	FAMGMatrixEntry matij;

    palist = NULL;

    rt = lt = 0.0;
    normr = norml = 0.0;
    FAMGMatrixIter miter(M,veci);
	miter(matij);
    mii = M[matij];
	while(miter(matij))
    {
        vecj = matij.dest();
        rj = M[matij];
        lj = M.GetAdjData(matij);
        mjj = M.DiagValue(matij.dest());
        normr += Abs(rj);
        norml += Abs(lj/mjj);
        rt += rj*tvA[vecj];
        lt += lj*tvB[vecj];
    }

    normr = normr/Abs(mii);

    if((normr < 1e-15) || (norml < 1e-15))
    {
        if(graph.SavePaList(palist,0,NULL,NULL,NULL,0.0)) 
			return 1;
        return 0;
    }
        
    nnb = 0;
	const FAMGMatrixAlg &MTemp = *GetTmpMatrix();
	FAMGMatrixIter mtiter(MTemp,veci);
    mtiter(matij);    // skip diagonal
    while(mtiter(matij))
       nnb++;

    FF0(veci,ff,gg, f, g);
    
    sigma = FAMGGetParameter()->Getsigma();
 
    FAMGMarkHeap(FAMG_FROM_TOP);
    sd = (struct FAMGSpecialData*) FAMGGetMem(nnb*sizeof(struct FAMGSpecialData), FAMG_FROM_TOP);
    if (sd == NULL)
		return 1;

    err1= (double*) FAMGGetMem(nnb*sizeof(double), FAMG_FROM_TOP);
    if (err1 == NULL)
		return 1;

    z = 0;
	mtiter.reset();
	mtiter(matij);// skip diagonal
    while(mtiter(matij))
    { 
        vecj = matij.dest();
		j = vecj.GetIndex();
        rj = MTemp[matij];
		lj = MTemp.GetAdjData(matij);
        mjj = M.DiagValue(matij.dest());
        JJ0(vecj,hjj,ejj);
        FJ0(j,vecj,fj,gj,f,g);
        sd[z].j = j;
        sd[z].hjj = hjj;
        sd[z].ejj = ejj;
        sd[z].fj = fj;
        sd[z].gj = gj;
        sd[z].rj = rj;
        sd[z].lj = lj/mjj;
        z++;
    }
    DFF0(veci);


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
            else
				min2 = err1[k]; 
        }           
    }

    if(min1 < FAMGGetParameter()->Geterror1())
    {
        graph.CorrectPaList(palist,min1/FAMGGetParameter()->Gettol()); 
    }
    else
    {
        sd2 = (struct FAMGSpecialData*) FAMGGetMem(nnb*sizeof(struct FAMGSpecialData), FAMG_FROM_TOP);
        if (sd2 == NULL)
			return 1;
        y = 0;
        for(k = 0; k < z; k++)
        {
            if(sigma*err1[k] < min2)
            {
                sd2[y] = sd[k];
                y++;
            }
        }
        graph.ClearPaListRev(palist);
        BestSecond0(palist,mii,rt,lt,ff,gg,sd2,y);
    } 

    FAMGReleaseHeap(FAMG_FROM_TOP);
    return 0;
}

#endif
