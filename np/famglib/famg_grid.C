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
#include "famg_algebra.h"
#include "famg_decomp.h"
#include "famg_transfer.h"
#include "famg_heap.h"
#include "famg_system.h"

#ifdef USE_UG_DS
extern "C"
{
#include "gm.h"
#include "ugm.h"
}
#endif

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
 
void FAMGGrid::Defect() const
{
    GetVector(FAMGDEFECT)->VecMinusMatVec(*GetVector(FAMGRHS), *GetMatrix(), *GetVector(FAMGUNKNOWN)); 
    return;
}

void FAMGGrid::Restriction(FAMGGrid *cg) const
// including smoothing of fine nodes
{
	FAMGVector &fgsolution = *GetVector(FAMGUNKNOWN);
	FAMGVector &fgdefect = *GetVector(FAMGDEFECT);
	FAMGVector &cgdefect = *(cg->GetVector(FAMGDEFECT));
	FAMGVectorEntry fvec;
	FAMGTransferEntry *transfc;
	const FAMGTransfer &transfer = *GetTransfer();
	FAMGTransferEntry *transfg;
	
	// jacobi smoothing for the fine nodes
    fgsolution.JacobiSmoothFG( *GetMatrix(), fgdefect );
	// correct defect
	Defect();

	cgdefect = 0.0;
	
	FAMGVectorIter fiter(GetGridVector());
	while( fiter(fvec) )
	    for(transfg = transfer.GetFirstEntry(fvec); transfg != NULL; transfg = transfg->GetNext())
			cgdefect[transfg->GetCol()] += transfg->GetRestriction() * fgdefect[fvec];
		
    return;
}
    
void FAMGGrid::Prolongation(const FAMGGrid *cg, FAMGVector *c = NULL)
// adds the prolongued solution-update to the fine grid solution
// including smoothing of fine nodes
{
	FAMGVector &fgsol = *GetVector(FAMGUNKNOWN);
	FAMGVector &fgdefect = *GetVector(FAMGDEFECT);
	const FAMGVector &cgsol = *(cg->GetVector(FAMGUNKNOWN));
	FAMGVectorEntry fvec;
	FAMGTransferEntry *transfc;
	const FAMGTransfer &transfer = *GetTransfer();
	FAMGTransferEntry *transfg;
	register double sum;
	
	FAMGVectorIter fiter(GetGridVector());
	while( fiter(fvec) )
	{
		sum = 0.0;
	    for(transfg = transfer.GetFirstEntry(fvec); transfg != NULL; transfg = transfg->GetNext())
        	sum += transfg->GetProlongation() * cgsol[transfg->GetCol()];
		fgsol[fvec] += sum;
	}

    
// prepare defect for jacobi smoothing
	Defect();
	
	if(c == NULL)
    {
        // famg as solver

        // jacobi smoothing for the fine nodes
        fgsol.JacobiSmoothFG( *GetMatrix(), fgdefect );
        
        // correct defect
        Defect();
    }
    else
    {
        // famg as transfer
        (*c) += fgsol;
        fgsol = 0.0;

        // jacobi smoothing for the fine nodes
        fgsol.JacobiSmoothFG( *GetMatrix(), fgdefect );

        // defect is computed in ug (Lmgc)
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
	GetVector(FAMGUNKNOWN)->dampedJacobiSmoother( *GetMatrix(), *GetVector(FAMGDEFECT) );
    return;
}

void FAMGGrid::FGSSmooth()
{
	GetVector(FAMGUNKNOWN)->FGSSmoother( *GetMatrix(), *GetVector(FAMGDEFECT) );
    return;
}

void FAMGGrid::BGSSmooth()
{
	GetVector(FAMGUNKNOWN)->BGSSmoother( *GetMatrix(), *GetVector(FAMGDEFECT) );
    return;
}

void FAMGGrid::SGSSmooth()
{
	GetVector(FAMGUNKNOWN)->SGSSmoother( *GetMatrix(), *GetVector(FAMGDEFECT) );
    return;
}

#ifdef FAMG_ILU
void FAMGGrid::ILUTSmooth()
{
    decomp->ILUT(vector[FAMGDEFECT]);
    return;
}
#endif

int FAMGGrid::SolveCoarseGrid()
{
#ifdef FAMG_BICG
	return BiCGStab();
#endif
	
	FAMGVector &sol = *GetVector(FAMGUNKNOWN);
	FAMGVector &def = *GetVector(FAMGDEFECT);
	FAMGVector &rhs = *GetVector(FAMGRHS);
	
	double mynorm = def.norm();
	double endnorm = mynorm * 1e-8;
	
	while (mynorm>endnorm)
	{
		SGSSmooth();
		Defect();
		mynorm  = def.norm();
	}
	
	return 0;
}

#ifdef FAMG_BICG
int FAMGGrid::BiCGStab()
// specialy apadted version for use as coarse grid solver
{
    double rlimit,alimit,reduction,limit,defectnorm,startdefect,oldnorm;
    double rho, oldrho, alpha, beta, omega, nenner;
    FAMGVector *vec[4];
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

	FAMGVector &defect = *GetVector(FAMGDEFECT);
	FAMGVector &solution = *GetVector(FAMGUNKNOWN);
	
    FAMGMarkHeap(FAMG_FROM_BOTTOM);
    vec[FAMGR] = solution.create_new();
    if(vec[FAMGR] == NULL) return 1;
    vec[FAMGV] = solution.create_new();
    if(vec[FAMGV] == NULL) return 1;
    vec[FAMGP] = solution.create_new();
    if(vec[FAMGP] == NULL) return 1;
    vec[FAMGT] = solution.create_new();
    if(vec[FAMGT] == NULL) return 1;

    maxit = 50;
    rlimit = 1e-10;
    alimit = 1e-14;
    reduction = 1e-10;

    *vec[FAMGR] = defect;
    defectnorm = vec[FAMGR]->norm();
    startdefect = defectnorm;
    limit = rlimit*startdefect;
    // ostr << "cg " << 0 << "\t" <<  startdefect << endl;
    // FAMGWrite(ostr);

    *vec[FAMGP] = *vec[FAMGR];
	
	rho = vec[FAMGR]->sum();
    // rho = FAMGScalProd(n,vector[FAMGRHS],vec[FAMGR]); // \tilde{r} = rhs
    if (Abs(rho) < 1e-10*alimit) 
    {
       ostr << __FILE__ << ", line " << __LINE__ << ": rho too small" << endl;
       FAMGWarning(ostr);
    }

    for(i = 0; i < maxit; i++)
    {
        *vector[FAMGDEFECT] = *vec[FAMGP];
        CGSmooth();
		vec[FAMGV]->MatVec(*GetMatrix(),defect);

        nenner = vec[FAMGV]->sum();// \tilde{r} = (1,...,1)
        // nenner = FAMGScalProd(n,vector[FAMGRHS],vec[FAMGV]); // \tilde{r} = rhs
        if (Abs(nenner) < 1e-15*Abs(rho)) 
        {
            ostr << __FILE__ << ", line " << __LINE__ << ": nenner too small" << endl;
            FAMGWarning(ostr);
        }
    
        alpha = rho/nenner;
		solution.AddScaledVec(alpha,defect);
		vec[FAMGR]->AddScaledVec(-alpha,*vec[FAMGV]);

        oldnorm = defectnorm; 
        defectnorm = oldnorm; // in order to avoid a warning !
        defectnorm = vec[FAMGR]->norm();
        // ostr << "cg " << i+0.5 << "\t" << defectnorm << "\t" << defectnorm/oldnorm;
        // ostr << "\t" << alpha  << endl;    
        // FAMGWrite(ostr);
        if((defectnorm < alimit) || (defectnorm < limit)) break;
        
		defect = *vec[FAMGR];
        CGSmooth();
		vec[FAMGT]->MatVec(*GetMatrix(),defect);

        omega = (*vec[FAMGT] * *vec[FAMGR]) / (*vec[FAMGT] * *vec[FAMGT]);
        if (Abs(omega) < 1e-15) 
        {
            ostr << __FILE__ << ", line " << __LINE__ << ": omega too small" << endl;
            FAMGWarning(ostr);
        }

		solution.AddScaledVec(omega,defect);
		vec[FAMGR]->AddScaledVec(-omega,*vec[FAMGT]);
		vec[FAMGP]->AddScaledVec(-omega,*vec[FAMGV]);

        oldnorm = defectnorm;
        defectnorm = vec[FAMGR]->norm();
        // ostr << "cg " << i+1 << "\t" << defectnorm << "\t" << defectnorm/oldnorm;
        // ostr << "\t" << omega << endl;    
        // FAMGWrite(ostr);
        if((defectnorm < alimit) || (defectnorm < limit)) break;

        oldrho = rho;
        rho = vec[FAMGR]->sum(); // \tilde{r} = (1,...,1)
        // rho = FAMGScalProd(n,vector[FAMGRHS],vec[FAMGR]); // \tilde{r} = rhs
        if (Abs(rho) < 1e-10*alimit) 
        {
            ostr << __FILE__ << ", line " << __LINE__ << ": rho too small" << endl;
            FAMGWarning(ostr);
        }

        beta = rho*alpha/(oldrho*omega);
        // could be accelerated
		*vec[FAMGP] *= beta;
		*vec[FAMGP] += *vec[FAMGR];
    }
	defect = *vec[FAMGR];

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
#endif

void FAMGGrid::SmoothTV()
{
    double sum, normA, mii, scaling;
    FAMGVector &tvA = *vector[FAMGTVA];
    FAMGVector &tvB = *vector[FAMGTVB];
	FAMGMatrixAlg &mat = *matrix;
	FAMGMatrixEntry me;
	FAMGVectorEntry ve;
    int k;

    const int stv = FAMGGetParameter()->Getstv();

    for(k = 0; k < stv; k++)
    {
		FAMGVectorIter tviter(tvA);
		while( tviter(ve) )
		{
			FAMGMatrixIter miter(*matrix,ve);
            sum = 0.0;
			
			miter(me);	// diagonal element
			mii = mat[me];
			while( miter(me) )
				sum -= mat[me]*tvA[me.dest()];
			tvA[ve] = sum/mii;
		}	
		
		FAMGVectorRevIter tvReviter(tvA);
		while( tvReviter(ve) )
		{
			FAMGMatrixIter miter(*matrix,ve);
            sum = 0.0;
			
			miter(me);	// diagonal element
			mii = mat[me];
			while( miter(me) )
				sum -= mat[me]*tvA[me.dest()];
			tvA[ve] = sum/mii;
		}	
    }
     
    normA = tvA.norm();
    scaling = sqrt((double)mat.GetN()) / normA;
	
	FAMGVectorIter tviter(tvA);
	while( tviter(ve) )
	{
		tvB[ve] = tvA[ve] *= scaling;
	}

	return;
}
    

#ifdef FAMG_ILU
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
#endif


void FAMGGrid::Stencil()
{
    int nn, nl;

    nn = matrix->GetN();
    nl = matrix->GetNLinks();
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

    //CGSmootherPtr = &FAMGGrid::ILUTSmooth;
    CGSmootherPtr = &FAMGGrid::JACSmooth;
    if(strcmp(cgsmoother,"ilut") == 0)
    {
		assert(0);
#ifdef FAMG_ILU
        CGSmootherPtr = &FAMGGrid::ILUTSmooth;
#endif
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
		assert(0);
#ifdef FAMG_ILU
        PreSmootherPtr = &FAMGGrid::ILUTSmooth;
#endif
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
		assert(0);
#ifdef FAMG_ILU
        PostSmootherPtr = &FAMGGrid::ILUTSmooth;
#endif
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

int FAMGGrid::InitLevel0(const class FAMGSystem &system)
{
    int i;
	FAMGVector *new_vector;

	n = system.GetN();
    nf = 0;
	
	mygridvector = system.GetGridVector();
	
#ifdef USE_UG_DS
	SetugGrid(system.GetFineGrid());
#else
    SetFather(NULL);
#endif
	
    matrix = system.GetMatrix();
    tmpmatrix = matrix;
	
#ifdef FAMG_ILU
    decomp = NULL;
#endif
	
    transfer = NULL;
    for(i = 0; i < FAMGMAXVECTORS; i++)
		vector[i] = system.GetVector(i);
	
#ifdef FAMG_REORDERmap	
    map = (int *) FAMGGetMem(n*sizeof(int),FAMG_FROM_TOP);
    if(map == NULL) return 1;
    for(i = 0; i < n; i++) map[i] = i;
#endif

    graph = NULL;

    GetSmoother();

    return 0;
}



int FAMGGrid::Init(int nn, const FAMGGrid& grid_pattern)
{
    int i;

    n = nn;
    nf = 0;
	
#ifdef USE_UG_DS
	GRID *new_grid;
	FAMGugVector *new_vector;
	
	new_grid = CreateNewLevelAMG(MYMG(grid_pattern.GetugGrid()));
	if( new_grid == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create coarser grid" << endl;
		FAMGError(ostr);
		assert(0);
	}
	SetugGrid(new_grid);
	
	mygridvector = new FAMGugGridVector(new_grid);
	if( mygridvector == NULL )
	{
		ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create gridvector" << endl;
		FAMGError(ostr);
		assert(0);
	}
	
	matrix = new FAMGugMatrix(new_grid, *(FAMGugMatrix*)grid_pattern.GetMatrix());
    if(matrix == NULL)
		return 1;
#else	
    matrix = (FAMGMatrix *) FAMGGetMem(sizeof(FAMGMatrix),FAMG_FROM_TOP);
    if(matrix == NULL)
		return 1;
    if(matrix->Init(n))
		return 1;
	mygridvector = NULL; // aendern
#endif
	
    tmpmatrix = matrix;
    transfer = NULL;
	
#ifdef FAMG_ILU
    decomp = NULL;
#endif
	
#ifndef USE_UG_DS	
    father = (int *) FAMGGetMem(n*sizeof(int),FAMG_FROM_TOP);
    if(father == NULL) return 1;
#endif
	
#ifdef FAMG_ILU
    map = NULL;
#endif
    graph = NULL;

    for(i = 0; i < FAMGMAXVECTORS; i++) // here we could save some memory
    {
#ifdef USE_UG_DS
		new_vector = new FAMGugVector(GetGridVector(),*(FAMGugVector*)grid_pattern.GetVector(i));
#else
		assert(0);
        //new_vector = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_TOP);
#endif
		if( new_vector == NULL )
		{
			ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": cannot create vector nr. " << i << endl;
			FAMGError(ostr);
			return 0;
		}
        SetVector(i,new_vector);
    }

#ifdef UG_DRAW
    vertex = (void **) FAMGGetMem(n*sizeof(void *),FAMG_FROM_TOP);
    if(vertex == NULL) return 1;
#endif
	
    GetSmoother();

    return 0;
}


void FAMGGrid::Deconstruct()
{
	int i;

    for(i = 0; i < FAMGMAXVECTORS; i++)
    {
#ifdef USE_UG_DS
		delete GetVector(i);
#else
		assert(0);
#endif
    }
	
#ifdef USE_UG_DS
	delete GetMatrix();
#else
	assert(0);
    father = NULL;
    tmpmatrix = matrix;
#endif
	
#ifdef FAMG_ILU
    decomp = NULL;
#endif
	
	delete mygridvector; mygridvector = NULL;
	
    return;
}    

int FAMGGrid::Construct(FAMGGrid *fg)
{
    int i, j;

    FAMGMatrixAlg *fmatrix = fg->GetMatrix();
    int fn = fg->GetN();
#ifdef UG_DRAW
    void  **vertexfg = fg->GetNode();
#endif
	FAMGVector &tvAcg = *GetVector(FAMGTVA);
	const FAMGVector &tvAfg = *(fg->GetVector(FAMGTVA));
	FAMGVector &tvBcg = *GetVector(FAMGTVB);
	const FAMGVector &tvBfg = *(fg->GetVector(FAMGTVB));

	if(fg->GetTransfer()->SetDestinationToCoarse(*fg,*this))
    {
        ostrstream ostr;
        ostr << __FILE__ << __LINE__  << "can not bend transfer entries to coarse grid" << endl;
        FAMGError(ostr);
        assert(0);
    }
	
	const FAMGGridVector &fg_gridvec = fg->GetGridVector();
	FAMGVectorIter viter(fg_gridvec);
	FAMGVectorEntry fg_ve, cg_ve;
	
	// transfer testvectors by trivial injection to coarse grid
    j = 0;
	while( viter(fg_ve) )
	{
		if (fg_gridvec.IsCG(fg_ve))
		{
			cg_ve = GetTransfer()->GetFirstEntry(fg_ve)->GetCol();
			tvAcg[cg_ve] = tvAfg[fg_ve];
			tvBcg[cg_ve] = tvBfg[fg_ve];
			j++;
		}
	}

    if (j != GetN())
    {
        ostrstream ostr;
        ostr << __FILE__ << __LINE__  << "number of coarse grid node doesn't match" << endl;
        FAMGError(ostr);
        assert(0);
    }
        
#ifdef UG_DRAW
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
#endif
	
    if(GetMatrix()->ConstructGalerkinMatrix(*fg)) 
		return 1;

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
	
//////////////////////////////////////////
#ifdef WEG
	// war nur testen drin
	FAMGNode *node;
	FAMGPaList* pa, *pf;
	int p1, p2, m;
	double d;
	
	node = graph->GetNode(10);
	p1=3; p2=17;	
	for(pa = node->GetPaList(); pa!=NULL; pa = pa->GetNext() )
		if( !((pa->GetPa(0)==p1 && pa->GetPa(1)==p2) || (pa->GetPa(0)==p2 && pa->GetPa(1)==p1)))
			pa->SetApprox(99999.0);
		else
		{
			pf = node->GetPaList();
			
			m=pa->GetNp(); pa->SetNp(pf->GetNp()); pf->SetNp(m);
			assert(m==2);
			m=pa->GetPa(0); pa->SetPa(0,pf->GetPa(1)); pf->SetPa(1,m);
			m=pa->GetPa(1); pa->SetPa(1,pf->GetPa(0)); pf->SetPa(0,m);
			d=pa->GetCoeff(0); pa->SetCoeff(0,pf->GetCoeff(1)); pf->SetCoeff(1,d);
			d=pa->GetCoeff(1); pa->SetCoeff(1,pf->GetCoeff(0)); pf->SetCoeff(0,d);
			d=pa->GetCoefft(0); pa->SetCoefft(0,pf->GetCoefft(1)); pf->SetCoefft(1,d);
			d=pa->GetCoefft(1); pa->SetCoefft(1,pf->GetCoefft(0)); pf->SetCoefft(0,d);
			//m=pa->GetPa(0); pa->SetPa(0,pf->GetPa(0)); pf->SetPa(0,m);
			//m=pa->GetPa(1); pa->SetPa(1,pf->GetPa(1)); pf->SetPa(1,m);
			//d=pa->GetCoeff(0); pa->SetCoeff(0,pf->GetCoeff(0)); pf->SetCoeff(0,d);
			//d=pa->GetCoeff(1); pa->SetCoeff(1,pf->GetCoeff(1)); pf->SetCoeff(1,d);
			//d=pa->GetCoefft(0); pa->SetCoefft(0,pf->GetCoefft(0)); pf->SetCoefft(0,d);
			//d=pa->GetCoefft(1); pa->SetCoefft(1,pf->GetCoefft(1)); pf->SetCoefft(1,d);
			d=pa->GetApprox(); pa->SetApprox(pf->GetApprox()); pf->SetApprox(d);
			m=pa->GetNewLinks(); pa->SetNewLinks(pf->GetNewLinks()); pf->SetNewLinks(m);
			d=pa->GetNewCG(); pa->SetNewCG(pf->GetNewCG()); pf->SetNewCG(d);
			
		}
	
	node = graph->GetNode(22);
	p1=15; p2=29;	
	for(pa = node->GetPaList(); pa!=NULL; pa = pa->GetNext() )
		if( !((pa->GetPa(0)==p1 && pa->GetPa(1)==p2) || (pa->GetPa(0)==p2 && pa->GetPa(1)==p1)))
			pa->SetApprox(99999.0);
		else
		{
			pf = node->GetPaList();
			
			m=pa->GetNp(); pa->SetNp(pf->GetNp()); pf->SetNp(m);
			assert(m==2);
			m=pa->GetPa(0); pa->SetPa(0,pf->GetPa(1)); pf->SetPa(1,m);
			m=pa->GetPa(1); pa->SetPa(1,pf->GetPa(0)); pf->SetPa(0,m);
			d=pa->GetCoeff(0); pa->SetCoeff(0,pf->GetCoeff(1)); pf->SetCoeff(1,d);
			d=pa->GetCoeff(1); pa->SetCoeff(1,pf->GetCoeff(0)); pf->SetCoeff(0,d);
			d=pa->GetCoefft(0); pa->SetCoefft(0,pf->GetCoefft(1)); pf->SetCoefft(1,d);
			d=pa->GetCoefft(1); pa->SetCoefft(1,pf->GetCoefft(0)); pf->SetCoefft(0,d);
			//m=pa->GetPa(0); pa->SetPa(0,pf->GetPa(0)); pf->SetPa(0,m);
			//m=pa->GetPa(1); pa->SetPa(1,pf->GetPa(1)); pf->SetPa(1,m);
			//d=pa->GetCoeff(0); pa->SetCoeff(0,pf->GetCoeff(0)); pf->SetCoeff(0,d);
			//d=pa->GetCoeff(1); pa->SetCoeff(1,pf->GetCoeff(1)); pf->SetCoeff(1,d);
			//d=pa->GetCoefft(0); pa->SetCoefft(0,pf->GetCoefft(0)); pf->SetCoefft(0,d);
			//d=pa->GetCoefft(1); pa->SetCoefft(1,pf->GetCoefft(1)); pf->SetCoefft(1,d);
			d=pa->GetApprox(); pa->SetApprox(pf->GetApprox()); pf->SetApprox(d);
			m=pa->GetNewLinks(); pa->SetNewLinks(pf->GetNewLinks()); pf->SetNewLinks(m);
			d=pa->GetNewCG(); pa->SetNewCG(pf->GetNewCG()); pf->SetNewCG(d);
			
		}
	
	node = graph->GetNode(26);
	p1=19; p2=33;	
	for(pa = node->GetPaList(); pa!=NULL; pa = pa->GetNext() )
		if( !((pa->GetPa(0)==p1 && pa->GetPa(1)==p2) || (pa->GetPa(0)==p2 && pa->GetPa(1)==p1)))
			pa->SetApprox(99999.0);
		else
		{
			pf = node->GetPaList();
			
			m=pa->GetNp(); pa->SetNp(pf->GetNp()); pf->SetNp(m);
			assert(m==2);
			m=pa->GetPa(0); pa->SetPa(0,pf->GetPa(1)); pf->SetPa(1,m);
			m=pa->GetPa(1); pa->SetPa(1,pf->GetPa(0)); pf->SetPa(0,m);
			d=pa->GetCoeff(0); pa->SetCoeff(0,pf->GetCoeff(1)); pf->SetCoeff(1,d);
			d=pa->GetCoeff(1); pa->SetCoeff(1,pf->GetCoeff(0)); pf->SetCoeff(0,d);
			d=pa->GetCoefft(0); pa->SetCoefft(0,pf->GetCoefft(1)); pf->SetCoefft(1,d);
			d=pa->GetCoefft(1); pa->SetCoefft(1,pf->GetCoefft(0)); pf->SetCoefft(0,d);
			//m=pa->GetPa(0); pa->SetPa(0,pf->GetPa(0)); pf->SetPa(0,m);
			//m=pa->GetPa(1); pa->SetPa(1,pf->GetPa(1)); pf->SetPa(1,m);
			//d=pa->GetCoeff(0); pa->SetCoeff(0,pf->GetCoeff(0)); pf->SetCoeff(0,d);
			//d=pa->GetCoeff(1); pa->SetCoeff(1,pf->GetCoeff(1)); pf->SetCoeff(1,d);
			//d=pa->GetCoefft(0); pa->SetCoefft(0,pf->GetCoefft(0)); pf->SetCoefft(0,d);
			//d=pa->GetCoefft(1); pa->SetCoefft(1,pf->GetCoefft(1)); pf->SetCoefft(1,d);
			d=pa->GetApprox(); pa->SetApprox(pf->GetApprox()); pf->SetApprox(d);
			m=pa->GetNewLinks(); pa->SetNewLinks(pf->GetNewLinks()); pf->SetNewLinks(m);
			d=pa->GetNewCG(); pa->SetNewCG(pf->GetNewCG()); pf->SetNewCG(d);
			
		}
	
	node = graph->GetNode(38);
	p1=31; p2=45;	
	for(pa = node->GetPaList(); pa!=NULL; pa = pa->GetNext() )
		if( !((pa->GetPa(0)==p1 && pa->GetPa(1)==p2) || (pa->GetPa(0)==p2 && pa->GetPa(1)==p1)))
			pa->SetApprox(99999.0);
		else
		{
			pf = node->GetPaList();
			
			m=pa->GetNp(); pa->SetNp(pf->GetNp()); pf->SetNp(m);
			assert(m==2);
			m=pa->GetPa(0); pa->SetPa(0,pf->GetPa(1)); pf->SetPa(1,m);
			m=pa->GetPa(1); pa->SetPa(1,pf->GetPa(0)); pf->SetPa(0,m);
			d=pa->GetCoeff(0); pa->SetCoeff(0,pf->GetCoeff(1)); pf->SetCoeff(1,d);
			d=pa->GetCoeff(1); pa->SetCoeff(1,pf->GetCoeff(0)); pf->SetCoeff(0,d);
			d=pa->GetCoefft(0); pa->SetCoefft(0,pf->GetCoefft(1)); pf->SetCoefft(1,d);
			d=pa->GetCoefft(1); pa->SetCoefft(1,pf->GetCoefft(0)); pf->SetCoefft(0,d);
			//m=pa->GetPa(0); pa->SetPa(0,pf->GetPa(0)); pf->SetPa(0,m);
			//m=pa->GetPa(1); pa->SetPa(1,pf->GetPa(1)); pf->SetPa(1,m);
			//d=pa->GetCoeff(0); pa->SetCoeff(0,pf->GetCoeff(0)); pf->SetCoeff(0,d);
			//d=pa->GetCoeff(1); pa->SetCoeff(1,pf->GetCoeff(1)); pf->SetCoeff(1,d);
			//d=pa->GetCoefft(0); pa->SetCoefft(0,pf->GetCoefft(0)); pf->SetCoefft(0,d);
			//d=pa->GetCoefft(1); pa->SetCoefft(1,pf->GetCoefft(1)); pf->SetCoefft(1,d);
			d=pa->GetApprox(); pa->SetApprox(pf->GetApprox()); pf->SetApprox(d);
			m=pa->GetNewLinks(); pa->SetNewLinks(pf->GetNewLinks()); pf->SetNewLinks(m);
			d=pa->GetNewCG(); pa->SetNewCG(pf->GetNewCG()); pf->SetNewCG(d);
			
		}
#endif
//////////////////////////////////////////
	
    if (graph->EliminateNodes(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
    
    for(i = 0; i < conloops; i++)
    {
        FAMGMarkHeap(FAMG_FROM_BOTTOM);
#ifdef USE_UG_DS
        assert(0);	// i don't want to have a second matrix
#else
		assert(0);	// todo: change
        tmpmatrix = (FAMGMatrixAlg *) FAMGGetMem(sizeof(FAMGMatrixAlg),FAMG_FROM_BOTTOM);
        if(tmpmatrix == NULL)
			return 1;
        if(tmpmatrix->Init2(n))
			return 1;
        if(tmpmatrix->TmpMatrix(matrix,transfer,graph))
			return 1;
#endif
 
        if (graph->Construct2(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
        if (graph->EliminateNodes(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
        FAMGReleaseHeap(FAMG_FROM_BOTTOM);
    }

    if (graph->RemainingNodes()) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
    
    nf = graph->GetNF();
    // not neceassary any more: matrix->MarkUnknowns(graph);
 
printim(GLEVEL(GetugGrid()));//?????????????????????????????????????????????
  

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

#ifdef FAMG_ILU
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
#endif
