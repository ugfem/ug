/****************************************************************************/
/*																			*/
/* File:      famg_grid.C													*/
/*																			*/
/* Purpose:   FAMG grid classes functions									*/
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
#include "famg_grid.h"
#include "famg_graph.h"
#include "famg_matrix.h"
#include "famg_decomp.h"
#include "famg_transfer.h"
#include "famg_heap.h"
#include "famg_system.h"


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

// Class FAMGGrid
 
void FAMGGrid::DevideFGDefect()
{
    matrix->DevideFGDefect(vector[FAMGUNKNOWN],vector[FAMGDEFECT]);

    return;
}
       
void FAMGGrid::Defect()
{
    matrix->VecMinusMatVec(vector[FAMGDEFECT], vector[FAMGRHS], vector[FAMGUNKNOWN]); 
    return;
}

void FAMGGrid::Restriction(FAMGGrid *cg) const
{
    FAMGTransferEntry *transij, *trans;
    FAMGMatrixPtr matik, matjk;
    double *cgdefect, *fgdefect, sum, tij;
    int ic, i, j, k, nc, *father;

    nc = cg->GetN();
    father = cg->GetFather();
    cgdefect = cg->GetVector(FAMGRHS);
    fgdefect = vector[FAMGDEFECT];
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
    
void FAMGGrid::Prolongation(const FAMGGrid *cg)
{
    FAMGTransferEntry *trans, *transij;
    FAMGMatrixPtr matij;
    double sum, *cgunknown, *fghelp, *fgunknown, mii;
    int i, j, nc, *father;

    nc = cg->GetN();
    father = cg->GetFather();
    cgunknown = cg->GetVector(FAMGUNKNOWN);

    fgunknown = vector[FAMGUNKNOWN];
    fghelp = vector[FAMGDEFECT];
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

void FAMGGrid::CGSmooth()
{
    (this->*CGSmootherPtr)();
    return;
}

void FAMGGrid::PreSmooth()
{
    (this->*PreSmootherPtr)();
    return;
}

void FAMGGrid::PostSmooth()
{
    (this->*PostSmootherPtr)();
    return;
}

void FAMGGrid::JACSmooth()
{
    matrix->JAC(vector[FAMGDEFECT]);
    return;
}

void FAMGGrid::FGSSmooth()
{
    matrix->FGS(vector[FAMGDEFECT]);
    return;
}

void FAMGGrid::BGSSmooth()
{
    matrix->BGS(vector[FAMGDEFECT]);
    return;
}

void FAMGGrid::SGSSmooth()
{
    matrix->SGS(vector[FAMGDEFECT]);
    return;
}

void FAMGGrid::ILUTSmooth()
{
    decomp->ILUT(vector[FAMGDEFECT]);
    return;
}

int FAMGGrid::BiCGStab()
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

    const int FAMGR = 0;
    const int FAMGV = 1;
    const int FAMGP = 2;
    const int FAMGT = 3;

    FAMGMarkHeap(FAMG_FROM_BOTTOM);
    vec[FAMGR] = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_BOTTOM);
    if(vec[FAMGR] == NULL) return 1;
    vec[FAMGV] = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_BOTTOM);
    if(vec[FAMGV] == NULL) return 1;
    vec[FAMGP] = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_BOTTOM);
    if(vec[FAMGP] == NULL) return 1;
    vec[FAMGT] = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_BOTTOM);
    if(vec[FAMGT] == NULL) return 1;

    maxit = 50;
    rlimit = 1e-10;
    alimit = 1e-14;
    reduction = 1e-10;

    
    FAMGCopyVector(n,vec[FAMGR],vector[FAMGDEFECT]);
    defectnorm = FAMGNorm(n,vec[FAMGR]);
    startdefect = defectnorm;
    limit = rlimit*startdefect;
    // ostr << "cg " << 0 << "\t" <<  startdefect << endl;
    // FAMGWrite(ostr);

    FAMGCopyVector(n,vec[FAMGP],vec[FAMGR]);
    
    rho = FAMGSum(n,vec[FAMGR]); // \tilde{r} = (1,...,1)
    // rho = FAMGScalProd(n,vector[FAMGRHS],vec[FAMGR]); // \tilde{r} = rhs
    if (Abs(rho) < 1e-10*alimit) 
    {
       ostr << __FILE__ << ", line " << __LINE__ << ": rho too small" << endl;
       FAMGWarning(ostr);
    }

    for(i = 0; i < maxit; i++)
    {
        FAMGCopyVector(n,vector[FAMGDEFECT],vec[FAMGP]);
        CGSmooth();
        matrix->Mult(vec[FAMGV],vector[FAMGDEFECT]);

        nenner = FAMGSum(n,vec[FAMGV]);// \tilde{r} = (1,...,1)
        // nenner = FAMGScalProd(n,vector[FAMGRHS],vec[FAMGV]); // \tilde{r} = rhs
        if (Abs(nenner) < 1e-15*Abs(rho)) 
        {
            ostr << __FILE__ << ", line " << __LINE__ << ": nenner too small" << endl;
            FAMGWarning(ostr);
        }
    
        alpha = rho/nenner;
        FAMGAddVector(n,vector[FAMGUNKNOWN],vector[FAMGDEFECT],alpha);
        FAMGAddVector(n,vec[FAMGR],vec[FAMGV],-alpha);

        oldnorm = defectnorm; 
        defectnorm = oldnorm; // in order to avoid a warning !
        defectnorm = FAMGNorm(n,vec[FAMGR]);
        // ostr << "cg " << i+0.5 << "\t" << defectnorm << "\t" << defectnorm/oldnorm;
        // ostr << "\t" << alpha  << endl;    
        // FAMGWrite(ostr);
        if((defectnorm < alimit) || (defectnorm < limit)) break;
        
        FAMGCopyVector(n,vector[FAMGDEFECT],vec[FAMGR]);
        CGSmooth();
        matrix->Mult(vec[FAMGT],vector[FAMGDEFECT]);

        omega = FAMGScalProd(n,vec[FAMGT],vec[FAMGR])/FAMGScalProd(n,vec[FAMGT],vec[FAMGT]);
        if (Abs(omega) < 1e-15) 
        {
            ostr << __FILE__ << ", line " << __LINE__ << ": omega too small" << endl;
            FAMGWarning(ostr);
        }

        FAMGAddVector(n,vector[FAMGUNKNOWN],vector[FAMGDEFECT],omega);
        FAMGAddVector(n,vec[FAMGR],vec[FAMGT],-omega);
        FAMGAddVector(n,vec[FAMGP],vec[FAMGV],-omega);

        oldnorm = defectnorm;
        defectnorm = FAMGNorm(n,vec[FAMGR]);
        // ostr << "cg " << i+1 << "\t" << defectnorm << "\t" << defectnorm/oldnorm;
        // ostr << "\t" << omega << endl;    
        // FAMGWrite(ostr);
        if((defectnorm < alimit) || (defectnorm < limit)) break;

        oldrho = rho;
        rho = FAMGSum(n,vec[FAMGR]); // \tilde{r} = (1,...,1)
        // rho = FAMGScalProd(n,vector[FAMGRHS],vec[FAMGR]); // \tilde{r} = rhs
        if (Abs(rho) < 1e-10*alimit) 
        {
            ostr << __FILE__ << ", line " << __LINE__ << ": rho too small" << endl;
            FAMGWarning(ostr);
        }

        beta = rho*alpha/(oldrho*omega);
        // could be accelerated
        FAMGMultVector(n,vec[FAMGP],beta);
        FAMGAddVector(n,vec[FAMGP],vec[FAMGR]);
    }
    FAMGCopyVector(n,vector[FAMGDEFECT],vec[FAMGR]);

    // ostr << "cg: " << i+1 << "  " << defectnorm/startdefect << endl << flush;
    // FAMGWrite(ostr);

    FAMGReleaseHeap(FAMG_FROM_BOTTOM);
         
    if ((defectnorm > startdefect*reduction) && (defectnorm > alimit))
    {
            ostr << __FILE__ << ", line " << __LINE__ << ": coarse grid defect reduction: " << defectnorm/startdefect << endl;
            FAMGWarning(ostr);
    }

    return 0;
}

void FAMGGrid::SmoothTV()
{
    FAMGMatrixPtr matij;
    double *tvA, *tvB, sumA, normA, nd, mii;
    int i,j,k;

    const int stv = FAMGGetParameter()->Getstv();
    tvA = vector[FAMGTVA];
    tvB = vector[FAMGTVB];

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
    

void FAMGGrid::CopyVector(int source, int dest)
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

void FAMGGrid::AddVector(int source, int dest)
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

void FAMGGrid::MultVector(int source, double factor)
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

void FAMGGrid::SubVector(int source, int dest)
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


void FAMGGrid::SetVector(int dest, double val)
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


int FAMGGrid::ILUTDecomp(int cgilut)
{   
    FAMGMarkHeap(FAMG_FROM_BOTTOM);
    if(graph == NULL)
    {
        graph = (FAMGGraph *) FAMGGetMem(sizeof(FAMGGraph),FAMG_FROM_BOTTOM);
        if(graph == NULL) return 1;
    }
    if(graph->Init(this)) return 1;
    if(graph->OrderILUT(matrix)) return 1;
    Order(graph->GetMap());
    graph = NULL;
    FAMGReleaseHeap(FAMG_FROM_BOTTOM);

    decomp = (FAMGDecomp *) FAMGGetMem(sizeof(FAMGDecomp),FAMG_FROM_TOP);
    if(decomp == NULL) return 1;

    if(decomp->Init(matrix)) return 1;
    if(decomp->Construct(cgilut)) return 1;

    
    return 0;
}


void FAMGGrid::Stencil()
{
    int nn, nl;

    nn = matrix->GetN();
    nl = matrix->GetNL();
    ostrstream ostr; 
    ostr << "unknowns: " << nn << "\t";
    ostr << "avg. stencil: " << (double)nl/(double)nn << endl;
    FAMGWrite(ostr);

    return;
}
               
void FAMGGrid::GetSmoother()
{
    char *cgsmoother = FAMGGetParameter()->Getcgsmoother();
    char *presmoother = FAMGGetParameter()->Getpresmoother();
    char *postsmoother = FAMGGetParameter()->Getpostsmoother();

    CGSmootherPtr = &FAMGGrid::ILUTSmooth;
    if(strcmp(cgsmoother,"ilut") == 0)
    {
        CGSmootherPtr = &FAMGGrid::ILUTSmooth;
    }
    else if(strcmp(cgsmoother,"fgs") == 0)
    {
        CGSmootherPtr = &FAMGGrid::FGSSmooth;
    }
    else if(strcmp(cgsmoother,"bgs") == 0)
    {
        CGSmootherPtr = &FAMGGrid::BGSSmooth;
    }
    else if(strcmp(cgsmoother,"sgs") == 0)
    {
        CGSmootherPtr = &FAMGGrid::SGSSmooth;
    }
    else if(strcmp(cgsmoother,"jac") == 0)
    {
        CGSmootherPtr = &FAMGGrid::JACSmooth;
    }
    else
    {
        ostrstream ostr;
        ostr << __FILE__ << __LINE__ <<  "cgsmoother = ilut" << endl;
        FAMGWarning(ostr);
    }

    PreSmootherPtr = &FAMGGrid::FGSSmooth;
    if(strcmp(presmoother,"ilut") == 0)
    {
        PreSmootherPtr = &FAMGGrid::ILUTSmooth;
    }
    else if(strcmp(presmoother,"fgs") == 0)
    {
        PreSmootherPtr = &FAMGGrid::FGSSmooth;
    }
    else if(strcmp(presmoother,"bgs") == 0)
    {
        PreSmootherPtr = &FAMGGrid::BGSSmooth;
    }
    else if(strcmp(presmoother,"sgs") == 0)
    {
        PreSmootherPtr = &FAMGGrid::SGSSmooth;
    }
    else if(strcmp(presmoother,"jac") == 0)
    {
        PreSmootherPtr = &FAMGGrid::JACSmooth;
    }
    else
    {
        ostrstream ostr;
        ostr << __FILE__ << __LINE__ <<  "presmoother = fgs" << endl;
        FAMGWarning(ostr);
    }

    PostSmootherPtr = &FAMGGrid::BGSSmooth;
    if(strcmp(postsmoother,"ilut") == 0)
    {
        PostSmootherPtr = &FAMGGrid::ILUTSmooth;
    }
    else if(strcmp(postsmoother,"fgs") == 0)
    {
        PostSmootherPtr = &FAMGGrid::FGSSmooth;
    }
    else if(strcmp(postsmoother,"bgs") == 0)
    {
        PostSmootherPtr = &FAMGGrid::BGSSmooth;
    }
    else if(strcmp(postsmoother,"sgs") == 0)
    {
        PostSmootherPtr = &FAMGGrid::SGSSmooth;
    }
    else if(strcmp(postsmoother,"jac") == 0)
    {
        PostSmootherPtr = &FAMGGrid::JACSmooth;
    }
    else
    {
        ostrstream ostr;
        ostr << __FILE__ << __LINE__ <<  "postsmoother = bgs" << endl;
        FAMGWarning(ostr);
    }


    return;
}

int  FAMGGrid::InitLevel0(const class FAMGSystem &system)
{
    int i;

    n = system.GetN();  
    nf = 0;
    matrix = system.GetMatrix();
    tmpmatrix = matrix;
    decomp = NULL;
    transfer = NULL;
    for(i = 0; i < FAMGMAXVECTORS; i++) vector[i] = system.GetVector(i);
    map = (int *) FAMGGetMem(n*sizeof(int),FAMG_FROM_TOP);
    if(map == NULL) return 1;
    for(i = 0; i < n; i++) map[i] = i;
    father = NULL;
    graph = NULL;
    vertex = system.GetExtra();

    GetSmoother();

    return 0;
}

int FAMGGrid::Init(int nn)
{
    int i;

    n = nn;
    nf = 0;
    matrix = (FAMGMatrix *) FAMGGetMem(sizeof(FAMGMatrix),FAMG_FROM_TOP);
    if(matrix == NULL) return 1;
    if(matrix->Init(n)) return 1;
    tmpmatrix = matrix;
    transfer = NULL;
    decomp = NULL;

    father = (int *) FAMGGetMem(n*sizeof(int),FAMG_FROM_TOP);
    if(father == NULL) return 1;

    map = NULL;
    graph = NULL;

    for(i = 0; i < FAMGMAXVECTORS; i++) // here we could save some memory
    {
        vector[i] = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_TOP);
        if(vector[i] == NULL) return 1;
    }
        

    vertex = (void **) FAMGGetMem(n*sizeof(void *),FAMG_FROM_TOP);
    if(vertex == NULL) return 1;
    
    GetSmoother();

    return 0;
}


void FAMGGrid::Deconstruct()
{
    tmpmatrix = matrix;
    decomp = NULL;
    father = NULL;

    return;
}    

int FAMGGrid::Construct(FAMGGrid *fg)
{
    int i, j;

    FAMGMatrix *fmatrix = fg->matrix;
    int fn = fg->GetN();
    void  **vertexfg = fg->GetNode();
    double *tvAcg = vector[FAMGTVA];
    double *tvAfg = fg->GetVector(FAMGTVA);
    double *tvBcg = vector[FAMGTVB];
    double *tvBfg = fg->GetVector(FAMGTVB);

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


int FAMGGrid::ConstructTransfer()
{
    int i, conloops;

    conloops = FAMGGetParameter()->Getconloops();
    transfer = (FAMGTransfer *) FAMGGetMem(sizeof(FAMGTransfer),FAMG_FROM_TOP);
    if(transfer == NULL) return 1;
    if (transfer->Init(this)) return 1;

    FAMGMarkHeap(FAMG_FROM_BOTTOM);

    graph = (FAMGGraph *) FAMGGetMem(sizeof(FAMGGraph),FAMG_FROM_BOTTOM);
    if (graph == NULL) {FAMGReleaseHeap(FAMG_FROM_BOTTOM);  return 1;}

    // test
    // FGSSmoothTV();

    if (graph->Init(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
    if (graph->Construct(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
    if (graph->EliminateNodes(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
    
    for(i = 0; i < conloops; i++)
    {
        FAMGMarkHeap(FAMG_FROM_BOTTOM);
        tmpmatrix = (FAMGMatrix *) FAMGGetMem(sizeof(FAMGMatrix),FAMG_FROM_BOTTOM);
        if(tmpmatrix == NULL) return 1;
        if(tmpmatrix->Init2(n)) return 1;
        if(tmpmatrix->TmpMatrix(matrix,transfer,graph)) return 1;
 
        if (graph->Construct2(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
        if (graph->EliminateNodes(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
        FAMGReleaseHeap(FAMG_FROM_BOTTOM);
    }

    if (graph->RemainingNodes(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
    
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

    FAMGReleaseHeap(FAMG_FROM_BOTTOM);


    return 0;
}
    
int FAMGGrid::OrderVector(int vn, int *mapping)
{
    double *helpvect, *vect;
    int i;

    FAMGMarkHeap(FAMG_FROM_TOP);
    helpvect = (double *) FAMGGetMem(n*sizeof(double), FAMG_FROM_TOP);
    if (helpvect == NULL) return 1;
    vect = vector[vn];
    if(vect != NULL)
    {
        for(i = 0; i < n; i++) helpvect[i] = vect[i];
        for(i = 0; i < n; i++) vect[mapping[i]] = helpvect[i];
    }

    FAMGReleaseHeap(FAMG_FROM_TOP);

    return 0;
}

int FAMGGrid::ReorderVector(int vn, int *mapping)
{
    double *helpvect, *vect;
    int i;

    FAMGMarkHeap(FAMG_FROM_TOP);
    helpvect = (double *) FAMGGetMem(n*sizeof(double), FAMG_FROM_TOP);
    if (helpvect == NULL) return 1;
    vect = vector[vn];
    if(vect != NULL)
    {
        for(i = 0; i < n; i++) helpvect[i] = vect[i];
        for(i = 0; i < n; i++) vect[i] = helpvect[mapping[i]];
    }

    FAMGReleaseHeap(FAMG_FROM_TOP);

    return 0;
}

int FAMGGrid::Order(int *mapping)
{
    void **helpvertex;
    int *helpfather, i;
    

    /* order father */
    if(father != NULL)
    {  
        FAMGMarkHeap(FAMG_FROM_TOP);
        helpfather = (int *) FAMGGetMem(n*sizeof(int), FAMG_FROM_TOP);
        if (helpfather == NULL) return 1;
        for(i = 0; i < n; i++) helpfather[i] = father[i];
        for(i = 0; i < n; i++) father[mapping[i]] = helpfather[i];
        FAMGReleaseHeap(FAMG_FROM_TOP);
    }

    /* order vertices */
    if(vertex != NULL)
    {
        FAMGMarkHeap(FAMG_FROM_TOP);
        helpvertex = (void **) FAMGGetMem(n*sizeof(void *), FAMG_FROM_TOP);
        if (helpvertex == NULL) return 1;
        for(i = 0; i < n; i++) helpvertex[i] = vertex[i];
        for(i = 0; i < n; i++) vertex[mapping[i]] = helpvertex[i];
        FAMGReleaseHeap(FAMG_FROM_TOP);
    }
        
    /* order test vectors */
    if(OrderVector(FAMGTVA,mapping)) return 1;
    if(OrderVector(FAMGTVB,mapping)) return 1;
   
    /* order matrix */
    if (matrix->Order(mapping)) return 1;

    /* order transfer */
    if(transfer != NULL)
    {
        if (transfer->Order(mapping)) return 1;
    }

    return 0;
}

int FAMGGrid::Reorder()
{
    void **helpvertex;
    int *helpfather, i;
    
    /* reorder father */
    if(father != NULL)
    {  
        FAMGMarkHeap(FAMG_FROM_TOP);
        helpfather = (int *) FAMGGetMem(n*sizeof(int), FAMG_FROM_TOP);
        if (helpfather == NULL) return 1;
        for(i = 0; i < n; i++) helpfather[i] = father[i];
        for(i = 0; i < n; i++) father[i] = helpfather[map[i]];
        FAMGReleaseHeap(FAMG_FROM_TOP);
    }

    /* reorder vertices */
    if(vertex != NULL)
    {
        FAMGMarkHeap(FAMG_FROM_TOP);
        helpvertex = (void **) FAMGGetMem(n*sizeof(void *), FAMG_FROM_TOP);
        if (helpvertex == NULL) return 1;
        for(i = 0; i < n; i++) helpvertex[i] = vertex[i];
        for(i = 0; i < n; i++) vertex[i] = helpvertex[map[i]];
        FAMGReleaseHeap(FAMG_FROM_TOP);
    }
        
    /* reorder test vectors */
    if(ReorderVector(FAMGTVA,map)) return 1;
    if(ReorderVector(FAMGTVB,map)) return 1;
    
    /* reorder matrix */
    if (matrix->Reorder(map)) return 1;

    /* reorder transfer */ // debug
    if(transfer != NULL)
    {
        if (transfer->Reorder(map)) return 1;
    }

    return 0;
}

    
    


