/****************************************************************************/
/*																			*/
/* File:      grid.C														*/
/*																			*/
/* Purpose:   CMG grid classes functions									*/
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

#include "misc.h"
#include "grid.h"
#include "graph.h"
#include "matrix.h"
#include "decomp.h"
#include "transfer.h"
#include "heap.h"
#include "system.h"


#ifdef UG_DRAW

extern "C"
{
#include "wpm.h"
#include "wop.h"
#include "connectuggrape.h"
#include "uginterface.h"
}

#endif

/* RCS_ID
$Header$
*/

// Class CMGGrid
 
void CMGGrid::DevideFGDefect()
{
    matrix->DevideFGDefect(vector[CMGUNKNOWN],vector[CMGDEFECT]);

    return;
}
       
void CMGGrid::Defect()
{
    matrix->VecMinusMatVec(vector[CMGDEFECT], vector[CMGRHS], vector[CMGUNKNOWN]); 
    return;
}

void CMGGrid::Restriction(CMGGrid *cg) const
{
    CMGTransferEntry *transij, *trans;
    CMGMatrixPtr matik, matjk;
    double *cgdefect, *fgdefect, sum, tij;
    int ic, i, j, k, nc, *father;

    nc = cg->GetN();
    father = cg->GetFather();
    cgdefect = cg->GetVector(CMGRHS);
    fgdefect = vector[CMGDEFECT];
    trans = transfer->GetRow();
 
    for(ic = 0; ic < nc; ic++)
    {
        i = father[ic];
        sum = fgdefect[i];
        for(transij = (trans+i)->GetNext(); transij != NULL; transij = transij->GetNext())
        {
            j = transij->GetId();
            tij = transij->GetData();
            matjk = matrix->GetStart(j);
            while(matjk.GetNext())
            {
                k = matjk.GetIndex();
                if(matjk.GetType()) // FG node
                {
                    sum -= tij*matjk.GetData()*fgdefect[k];
                }
            }
        }
        matik = matrix->GetStart(i);
        while(matik.GetNext())
        {
            k = matik.GetIndex();
            if(matik.GetType())
            {
                sum -= matik.GetData()*fgdefect[k];
            }
        }
            
        cgdefect[ic] = sum;
    }

    return;
}
    
void CMGGrid::Prolongation(const CMGGrid *cg)
{
    CMGTransferEntry *trans, *transij;
    CMGMatrixPtr matij;
    double sum, *cgunknown, *fghelp, *fgunknown, mii;
    int i, j, nc, *father;

    nc = cg->GetN();
    father = cg->GetFather();
    cgunknown = cg->GetVector(CMGUNKNOWN);

    fgunknown = vector[CMGUNKNOWN];
    fghelp = vector[CMGDEFECT];
    trans = transfer->GetRow();
    
    for(i = 0; i < nc; i++)
    {
        fghelp[father[i]] = cgunknown[i];
        fgunknown[father[i]] += cgunknown[i];
    }
    for(i = 0; i < n; i++)
    {
        if (matrix->GetType(i)) // FG node
        {
            for(transij = (trans+i)->GetNext(); transij != NULL; transij = transij->GetNext())
            {
                j = transij->GetId();
                fghelp[i] += transij->GetData() * fghelp[j];
            }
        }
    }
    for(i = 0; i < n; i++)
    {
        matij = matrix->GetStart(i);
        if(matij.GetType())
        {
            sum = 0.0;
            mii = matij.GetData();
            while(matij.GetNext())
            {
                j = matij.GetIndex();
                sum += matij.GetData() * fghelp[j];
            }
            fgunknown[i] -= sum/mii;
        }
    }

    return;
}

void CMGGrid::CGSmooth()
{
    (this->*CGSmootherPtr)();
    return;
}

void CMGGrid::PreSmooth()
{
    (this->*PreSmootherPtr)();
    return;
}

void CMGGrid::PostSmooth()
{
    (this->*PostSmootherPtr)();
    return;
}

void CMGGrid::JACSmooth()
{
    matrix->JAC(vector[CMGDEFECT]);
    return;
}

void CMGGrid::FGSSmooth()
{
    matrix->FGS(vector[CMGDEFECT]);
    return;
}

void CMGGrid::BGSSmooth()
{
    matrix->BGS(vector[CMGDEFECT]);
    return;
}

void CMGGrid::SGSSmooth()
{
    matrix->SGS(vector[CMGDEFECT]);
    return;
}

void CMGGrid::ILUTSmooth()
{
    decomp->ILUT(vector[CMGDEFECT]);
    return;
}

int CMGGrid::BiCGStab()
{
    double rlimit,alimit,reduction,limit,defectnorm,startdefect,oldnorm;
    double *vec[4], rho, oldrho, alpha, beta, omega, nenner;
    int maxit,i;
    ostrstream ostr; 

    // at this point we could call an extern solver e.g.
    // rhs = def
    // def = M^{-1} def
    // u = u + def
    // def = rhs - K u

    const int CMGR = 0;
    const int CMGV = 1;
    const int CMGP = 2;
    const int CMGT = 3;

    CMGMarkHeap(CMG_FROM_BOTTOM);
    vec[CMGR] = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_BOTTOM);
    if(vec[CMGR] == NULL) return 1;
    vec[CMGV] = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_BOTTOM);
    if(vec[CMGV] == NULL) return 1;
    vec[CMGP] = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_BOTTOM);
    if(vec[CMGP] == NULL) return 1;
    vec[CMGT] = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_BOTTOM);
    if(vec[CMGT] == NULL) return 1;

    maxit = 50;
    rlimit = 1e-10;
    alimit = 1e-14;
    reduction = 1e-10;

    
    CMGCopyVector(n,vec[CMGR],vector[CMGDEFECT]);
    defectnorm = CMGNorm(n,vec[CMGR]);
    startdefect = defectnorm;
    limit = rlimit*startdefect;
    // ostr << "cg " << 0 << "\t" <<  startdefect << endl;
    // CMGWrite(ostr);

    CMGCopyVector(n,vec[CMGP],vec[CMGR]);
    
    rho = CMGSum(n,vec[CMGR]); // \tilde{r} = (1,...,1)
    // rho = CMGScalProd(n,vector[CMGRHS],vec[CMGR]); // \tilde{r} = rhs
    if (Abs(rho) < 1e-10*alimit) 
    {
       ostr << __FILE__ << ", line " << __LINE__ << ": rho too small" << endl;
       CMGWarning(ostr);
    }

    for(i = 0; i < maxit; i++)
    {
        CMGCopyVector(n,vector[CMGDEFECT],vec[CMGP]);
        CGSmooth();
        matrix->Mult(vec[CMGV],vector[CMGDEFECT]);

        nenner = CMGSum(n,vec[CMGV]);// \tilde{r} = (1,...,1)
        // nenner = CMGScalProd(n,vector[CMGRHS],vec[CMGV]); // \tilde{r} = rhs
        if (Abs(nenner) < 1e-15*Abs(rho)) 
        {
            ostr << __FILE__ << ", line " << __LINE__ << ": nenner too small" << endl;
            CMGWarning(ostr);
        }
    
        alpha = rho/nenner;
        CMGAddVector(n,vector[CMGUNKNOWN],vector[CMGDEFECT],alpha);
        CMGAddVector(n,vec[CMGR],vec[CMGV],-alpha);

        oldnorm = defectnorm; 
        defectnorm = oldnorm; // in order to avoid a warning !
        defectnorm = CMGNorm(n,vec[CMGR]);
        // ostr << "cg " << i+0.5 << "\t" << defectnorm << "\t" << defectnorm/oldnorm;
        // ostr << "\t" << alpha  << endl;    
        // CMGWrite(ostr);
        if((defectnorm < alimit) || (defectnorm < limit)) break;
        
        CMGCopyVector(n,vector[CMGDEFECT],vec[CMGR]);
        CGSmooth();
        matrix->Mult(vec[CMGT],vector[CMGDEFECT]);

        omega = CMGScalProd(n,vec[CMGT],vec[CMGR])/CMGScalProd(n,vec[CMGT],vec[CMGT]);
        if (Abs(omega) < 1e-15) 
        {
            ostr << __FILE__ << ", line " << __LINE__ << ": omega too small" << endl;
            CMGWarning(ostr);
        }

        CMGAddVector(n,vector[CMGUNKNOWN],vector[CMGDEFECT],omega);
        CMGAddVector(n,vec[CMGR],vec[CMGT],-omega);
        CMGAddVector(n,vec[CMGP],vec[CMGV],-omega);

        oldnorm = defectnorm;
        defectnorm = CMGNorm(n,vec[CMGR]);
        // ostr << "cg " << i+1 << "\t" << defectnorm << "\t" << defectnorm/oldnorm;
        // ostr << "\t" << omega << endl;    
        // CMGWrite(ostr);
        if((defectnorm < alimit) || (defectnorm < limit)) break;

        oldrho = rho;
        rho = CMGSum(n,vec[CMGR]); // \tilde{r} = (1,...,1)
        // rho = CMGScalProd(n,vector[CMGRHS],vec[CMGR]); // \tilde{r} = rhs
        if (Abs(rho) < 1e-10*alimit) 
        {
            ostr << __FILE__ << ", line " << __LINE__ << ": rho too small" << endl;
            CMGWarning(ostr);
        }

        beta = rho*alpha/(oldrho*omega);
        // could be accelerated
        CMGMultVector(n,vec[CMGP],beta);
        CMGAddVector(n,vec[CMGP],vec[CMGR]);
    }
    CMGCopyVector(n,vector[CMGDEFECT],vec[CMGR]);

    // ostr << "cg: " << i+1 << "  " << defectnorm/startdefect << endl << flush;
    // CMGWrite(ostr);

    CMGReleaseHeap(CMG_FROM_BOTTOM);
         
    if ((defectnorm > startdefect*reduction) && (defectnorm > alimit))
    {
            ostr << __FILE__ << ", line " << __LINE__ << ": coarse grid defect reduction: " << defectnorm/startdefect << endl;
            CMGWarning(ostr);
    }

    return 0;
}

void CMGGrid::SmoothTV()
{
    CMGMatrixPtr matij;
    double *tvA, *tvB, sumA, normA, nd, mii;
    int i,j,k;

    const int stv = CMGGetParameter()->Getstv();
    tvA = vector[CMGTVA];
    tvB = vector[CMGTVB];

    for(k = 0; k < stv; k++)
    {
        for(i = 0; i < n; i++)
        {
            sumA = 0.0;
            matij = matrix->GetStart(i);
            mii = matij.GetData();
            while(matij.GetNext())
            {
                j = matij.GetIndex();
                sumA -= matij.GetData()*tvA[j];
            }
            tvA[i] = sumA/mii;
        }

        for(i = n-1; i >= 0; i--)
        {
            sumA = 0.0;
            matij = matrix->GetStart(i);
            mii = matij.GetData();
            while(matij.GetNext())
            {
                j = matij.GetIndex();
                sumA -= matij.GetData()*tvA[j];
            }
            tvA[i] = sumA/mii;
        }
    }
     
    normA = 0.0;
    for(i = 0; i < n; i++) 
    {
        normA += tvA[i]*tvA[i];
    }
    
    normA = sqrt(normA);
    nd = sqrt((double)n);
    for(i = 0; i < n; i++) 
    {
        tvB[i] = tvA[i]=nd*tvA[i]/normA;
     }

    
    return;
}
    

void CMGGrid::CopyVector(int source, int dest)
{
    double *vecsource, *vecdest;
    int i;

    vecsource = vector[source];
    vecdest = vector[dest];
    
    for(i = 0; i < n; i++) 
    {
        vecdest[i] = vecsource[i];
    }

    return;
}

void CMGGrid::AddVector(int source, int dest)
{
    double *vecsource, *vecdest;
    int i;

    vecsource = vector[source];
    vecdest = vector[dest];
    
    for(i = 0; i < n; i++) 
    {
        vecdest[i] += vecsource[i];
    }

    return;
}

void CMGGrid::MultVector(int source, double factor)
{
    double *vecsource;
    int i;

    vecsource = vector[source];
    
    for(i = 0; i < n; i++) 
    {
        vecsource[i] = factor*vecsource[i];
    }

    return;
}

void CMGGrid::SubVector(int source, int dest)
{
    double *vecsource, *vecdest;
    int i;

    vecsource = vector[source];
    vecdest = vector[dest];
    
    for(i = 0; i < n; i++) 
    {
        vecdest[i] -= vecsource[i];
    }

    return;
}


void CMGGrid::SetVector(int dest, double val)
{
    double *vecdest;
    int i;

    vecdest = vector[dest];
    
    for(i = 0; i < n; i++) 
    {
        vecdest[i] = val;
    }

    return;
}


int CMGGrid::ILUTDecomp(int cgilut)
{   
    CMGMarkHeap(CMG_FROM_BOTTOM);
    if(graph == NULL)
    {
        graph = (CMGGraph *) CMGGetMem(sizeof(CMGGraph),CMG_FROM_BOTTOM);
        if(graph == NULL) return 1;
    }
    if(graph->Init(this)) return 1;
    if(graph->OrderILUT(matrix)) return 1;
    Order(graph->GetMap());
    graph = NULL;
    CMGReleaseHeap(CMG_FROM_BOTTOM);

    decomp = (CMGDecomp *) CMGGetMem(sizeof(CMGDecomp),CMG_FROM_TOP);
    if(decomp == NULL) return 1;

    if(decomp->Init(matrix)) return 1;
    if(decomp->Construct(cgilut)) return 1;

    
    return 0;
}


void CMGGrid::Stencil()
{
    int nn, nl;

    nn = matrix->GetN();
    nl = matrix->GetNL();
    ostrstream ostr; 
    ostr << "unknowns: " << nn << "\t";
    ostr << "avg. stencil: " << (double)nl/(double)nn << endl;
    CMGWrite(ostr);

    return;
}
               
void CMGGrid::GetSmoother()
{
    char *cgsmoother = CMGGetParameter()->Getcgsmoother();
    char *presmoother = CMGGetParameter()->Getpresmoother();
    char *postsmoother = CMGGetParameter()->Getpostsmoother();

    CGSmootherPtr = &CMGGrid::ILUTSmooth;
    if(strcmp(cgsmoother,"ilut") == 0)
    {
        CGSmootherPtr = &CMGGrid::ILUTSmooth;
    }
    else if(strcmp(cgsmoother,"fgs") == 0)
    {
        CGSmootherPtr = &CMGGrid::FGSSmooth;
    }
    else if(strcmp(cgsmoother,"bgs") == 0)
    {
        CGSmootherPtr = &CMGGrid::BGSSmooth;
    }
    else if(strcmp(cgsmoother,"sgs") == 0)
    {
        CGSmootherPtr = &CMGGrid::SGSSmooth;
    }
    else if(strcmp(cgsmoother,"jac") == 0)
    {
        CGSmootherPtr = &CMGGrid::JACSmooth;
    }
    else
    {
        ostrstream ostr;
        ostr << __FILE__ << __LINE__ <<  "cgsmoother = ilut" << endl;
        CMGWarning(ostr);
    }

    PreSmootherPtr = &CMGGrid::FGSSmooth;
    if(strcmp(presmoother,"ilut") == 0)
    {
        PreSmootherPtr = &CMGGrid::ILUTSmooth;
    }
    else if(strcmp(presmoother,"fgs") == 0)
    {
        PreSmootherPtr = &CMGGrid::FGSSmooth;
    }
    else if(strcmp(presmoother,"bgs") == 0)
    {
        PreSmootherPtr = &CMGGrid::BGSSmooth;
    }
    else if(strcmp(presmoother,"sgs") == 0)
    {
        PreSmootherPtr = &CMGGrid::SGSSmooth;
    }
    else if(strcmp(presmoother,"jac") == 0)
    {
        PreSmootherPtr = &CMGGrid::JACSmooth;
    }
    else
    {
        ostrstream ostr;
        ostr << __FILE__ << __LINE__ <<  "presmoother = fgs" << endl;
        CMGWarning(ostr);
    }

    PostSmootherPtr = &CMGGrid::BGSSmooth;
    if(strcmp(postsmoother,"ilut") == 0)
    {
        PostSmootherPtr = &CMGGrid::ILUTSmooth;
    }
    else if(strcmp(postsmoother,"fgs") == 0)
    {
        PostSmootherPtr = &CMGGrid::FGSSmooth;
    }
    else if(strcmp(postsmoother,"bgs") == 0)
    {
        PostSmootherPtr = &CMGGrid::BGSSmooth;
    }
    else if(strcmp(postsmoother,"sgs") == 0)
    {
        PostSmootherPtr = &CMGGrid::SGSSmooth;
    }
    else if(strcmp(postsmoother,"jac") == 0)
    {
        PostSmootherPtr = &CMGGrid::JACSmooth;
    }
    else
    {
        ostrstream ostr;
        ostr << __FILE__ << __LINE__ <<  "postsmoother = bgs" << endl;
        CMGWarning(ostr);
    }


    return;
}

int  CMGGrid::InitLevel0(const class CMGSystem &system)
{
    int i;

    n = system.GetN();  
    nf = 0;
    matrix = system.GetMatrix();
    tmpmatrix = matrix;
    decomp = NULL;
    transfer = NULL;
    for(i = 0; i < CMGMAXVECTORS; i++) vector[i] = system.GetVector(i);
    map = (int *) CMGGetMem(n*sizeof(int),CMG_FROM_TOP);
    if(map == NULL) return 1;
    for(i = 0; i < n; i++) map[i] = i;
    father = NULL;
    graph = NULL;
    vertex = system.GetExtra();

    GetSmoother();

    return 0;
}

int CMGGrid::Init(int nn)
{
    int i;

    n = nn;
    nf = 0;
    matrix = (CMGMatrix *) CMGGetMem(sizeof(CMGMatrix),CMG_FROM_TOP);
    if(matrix == NULL) return 1;
    if(matrix->Init(n)) return 1;
    tmpmatrix = matrix;
    transfer = NULL;
    decomp = NULL;

    father = (int *) CMGGetMem(n*sizeof(int),CMG_FROM_TOP);
    if(father == NULL) return 1;

    map = NULL;
    graph = NULL;

    for(i = 0; i < CMGMAXVECTORS; i++) // here we could save some memory
    {
        vector[i] = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_TOP);
        if(vector[i] == NULL) return 1;
    }
        

    vertex = (void **) CMGGetMem(n*sizeof(void *),CMG_FROM_TOP);
    if(vertex == NULL) return 1;
    
    GetSmoother();

    return 0;
}


void CMGGrid::Deconstruct()
{
    tmpmatrix = matrix;
    decomp = NULL;
    father = NULL;

    return;
}    

int CMGGrid::Construct(CMGGrid *fg)
{
    int i, j;

    CMGMatrix *fmatrix = fg->matrix;
    int fn = fg->GetN();
    void  **vertexfg = fg->GetNode();
    double *tvAcg = vector[CMGTVA];
    double *tvAfg = fg->GetVector(CMGTVA);
    double *tvBcg = vector[CMGTVB];
    double *tvBfg = fg->GetVector(CMGTVB);

    j = 0;
    for(i = 0; i < fn; i++)
    {
        if(fmatrix->GetType(i) == 0)
        {
            father[j] = i;
            tvAcg[j] = tvAfg[i];
            tvBcg[j] = tvBfg[i];
            j++;
        }
    }

    if(vertexfg != NULL)
    {
        j = 0;
        for(i = 0; i < fn; i++)
        {
            if(fmatrix->GetType(i) == 0)
            {
                vertex[j] = vertexfg[i];
                j++;
            }
        }
    }

    if(matrix->CGMatrix(fg->GetMatrix(), fg->GetTransfer(), father)) return 1;

    return 0;
}


int CMGGrid::ConstructTransfer()
{
    int i, conloops;

    conloops = CMGGetParameter()->Getconloops();
    transfer = (CMGTransfer *) CMGGetMem(sizeof(CMGTransfer),CMG_FROM_TOP);
    if(transfer == NULL) return 1;
    if (transfer->Init(this)) return 1;

    CMGMarkHeap(CMG_FROM_BOTTOM);

    graph = (CMGGraph *) CMGGetMem(sizeof(CMGGraph),CMG_FROM_BOTTOM);
    if (graph == NULL) {CMGReleaseHeap(CMG_FROM_BOTTOM);  return 1;}

    // test
    // FGSSmoothTV();

    if (graph->Init(this)) { CMGReleaseHeap(CMG_FROM_BOTTOM); return 1;}
    if (graph->Construct(this)) { CMGReleaseHeap(CMG_FROM_BOTTOM); return 1;}
    if (graph->EliminateNodes(this)) { CMGReleaseHeap(CMG_FROM_BOTTOM); return 1;}
    
    for(i = 0; i < conloops; i++)
    {
        CMGMarkHeap(CMG_FROM_BOTTOM);
        tmpmatrix = (CMGMatrix *) CMGGetMem(sizeof(CMGMatrix),CMG_FROM_BOTTOM);
        if(tmpmatrix == NULL) return 1;
        if(tmpmatrix->Init2(n)) return 1;
        if(tmpmatrix->TmpMatrix(matrix,transfer,graph)) return 1;
 
        if (graph->Construct2(this)) { CMGReleaseHeap(CMG_FROM_BOTTOM); return 1;}
        if (graph->EliminateNodes(this)) { CMGReleaseHeap(CMG_FROM_BOTTOM); return 1;}
        CMGReleaseHeap(CMG_FROM_BOTTOM);
    }

    if (graph->RemainingNodes(this)) { CMGReleaseHeap(CMG_FROM_BOTTOM); return 1;}
    
    nf = graph->GetNF();
    matrix->MarkUnknowns(graph);
 
    

#ifdef UG_DRAW
    /* test */
    PICTURE *thePic;
    thePic = GetCurrentPicture();
    if (thePic!=NULL)
    {
        DrawUgPicture(thePic);
        DrawPictureFrame(thePic,WOP_ACTIVE);
    } 
#endif

    CMGReleaseHeap(CMG_FROM_BOTTOM);


    return 0;
}
    
int CMGGrid::OrderVector(int vn, int *mapping)
{
    double *helpvect, *vect;
    int i;

    CMGMarkHeap(CMG_FROM_TOP);
    helpvect = (double *) CMGGetMem(n*sizeof(double), CMG_FROM_TOP);
    if (helpvect == NULL) return 1;
    vect = vector[vn];
    if(vect != NULL)
    {
        for(i = 0; i < n; i++) helpvect[i] = vect[i];
        for(i = 0; i < n; i++) vect[mapping[i]] = helpvect[i];
    }

    CMGReleaseHeap(CMG_FROM_TOP);

    return 0;
}

int CMGGrid::ReorderVector(int vn, int *mapping)
{
    double *helpvect, *vect;
    int i;

    CMGMarkHeap(CMG_FROM_TOP);
    helpvect = (double *) CMGGetMem(n*sizeof(double), CMG_FROM_TOP);
    if (helpvect == NULL) return 1;
    vect = vector[vn];
    if(vect != NULL)
    {
        for(i = 0; i < n; i++) helpvect[i] = vect[i];
        for(i = 0; i < n; i++) vect[i] = helpvect[mapping[i]];
    }

    CMGReleaseHeap(CMG_FROM_TOP);

    return 0;
}

int CMGGrid::Order(int *mapping)
{
    void **helpvertex;
    int *helpfather, i;
    

    /* order father */
    if(father != NULL)
    {  
        CMGMarkHeap(CMG_FROM_TOP);
        helpfather = (int *) CMGGetMem(n*sizeof(int), CMG_FROM_TOP);
        if (helpfather == NULL) return 1;
        for(i = 0; i < n; i++) helpfather[i] = father[i];
        for(i = 0; i < n; i++) father[mapping[i]] = helpfather[i];
        CMGReleaseHeap(CMG_FROM_TOP);
    }

    /* order vertices */
    if(vertex != NULL)
    {
        CMGMarkHeap(CMG_FROM_TOP);
        helpvertex = (void **) CMGGetMem(n*sizeof(void *), CMG_FROM_TOP);
        if (helpvertex == NULL) return 1;
        for(i = 0; i < n; i++) helpvertex[i] = vertex[i];
        for(i = 0; i < n; i++) vertex[mapping[i]] = helpvertex[i];
        CMGReleaseHeap(CMG_FROM_TOP);
    }
        
    /* order test vectors */
    if(OrderVector(CMGTVA,mapping)) return 1;
    if(OrderVector(CMGTVB,mapping)) return 1;
   
    /* order matrix */
    if (matrix->Order(mapping)) return 1;

    /* order transfer */
    if(transfer != NULL)
    {
        if (transfer->Order(mapping)) return 1;
    }

    return 0;
}

int CMGGrid::Reorder()
{
    void **helpvertex;
    int *helpfather, i;
    
    /* reorder father */
    if(father != NULL)
    {  
        CMGMarkHeap(CMG_FROM_TOP);
        helpfather = (int *) CMGGetMem(n*sizeof(int), CMG_FROM_TOP);
        if (helpfather == NULL) return 1;
        for(i = 0; i < n; i++) helpfather[i] = father[i];
        for(i = 0; i < n; i++) father[i] = helpfather[map[i]];
        CMGReleaseHeap(CMG_FROM_TOP);
    }

    /* reorder vertices */
    if(vertex != NULL)
    {
        CMGMarkHeap(CMG_FROM_TOP);
        helpvertex = (void **) CMGGetMem(n*sizeof(void *), CMG_FROM_TOP);
        if (helpvertex == NULL) return 1;
        for(i = 0; i < n; i++) helpvertex[i] = vertex[i];
        for(i = 0; i < n; i++) vertex[i] = helpvertex[map[i]];
        CMGReleaseHeap(CMG_FROM_TOP);
    }
        
    /* reorder test vectors */
    if(ReorderVector(CMGTVA,map)) return 1;
    if(ReorderVector(CMGTVB,map)) return 1;
    
    /* reorder matrix */
    if (matrix->Reorder(map)) return 1;

    /* reorder transfer */ // debug
    if(transfer != NULL)
    {
        if (transfer->Reorder(map)) return 1;
    }

    return 0;
}

    
    


