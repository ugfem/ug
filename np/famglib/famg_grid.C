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

#ifdef ModelP
// forward declaration for ug functions
extern "C"{
void VectorXferCopy (DDD_OBJ obj, DDD_PROC proc, DDD_PRIO prio);
void VectorGatherMatX (DDD_OBJ obj, int cnt, DDD_TYPE type_id, char **Data);
void VectorScatterConnX (DDD_OBJ obj, int cnt, DDD_TYPE type_id, char **Data, int newness);
}
#endif

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
    
void FAMGGrid::Prolongation(const FAMGGrid *cg, FAMGVector *c)
// adds the prolongued solution-update to the fine grid solution
// including smoothing of fine nodes
// c is set for FAMG transfer; not used for FAMG solver
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
	
#ifdef ModelP
	// in parallel now only the nodes in the border of the core partition are in the list
	
	VECTOR *vec;
	MATRIX *mat;
	FAMGNode *nodei;
	
	for ( int p=0; p<procs; p++)
	{
		if( me == p )
		{
		    if (graph->EliminateNodes(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
    		if (graph->RemainingNodes()) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
		}
		
		CommunicateNodeStatus();
	
	    for(i = 0; i < n; i++)
    	{
        	nodei = graph->GetNode(i);
			vec = ((FAMGugVectorEntryRef*)(nodei->GetVec().GetPointer()))->myvector();
		
			if( IS_FAMG_GHOST(vec) )
				continue; // only master vectors can be in border the of the core partition
		
			// now only the unmarked nodes must be inserted into the list
			if( !nodei->IsCGNode() && !nodei->IsFGNode() )
				if(graph->InsertNode(this, nodei))
					return 0;
	    }
	}	
#endif
	
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
    
#ifdef ModelP
	// update the ghost and border nodes
	CommunicateNodeStatus();
#endif
	
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

// *****************************************************************************
// *********** parallel extensions *********************************************
// *****************************************************************************

#ifdef ModelP
void FAMGVectorXferCopy (DDD_OBJ obj, DDD_PROC proc, DDD_PRIO prio)
// derived from dddif/handler.c/VectorXferCopy()
{
	INT		nmat	= 0;
	MATRIX	*imat;
	VECTOR	*pv		= (VECTOR *)obj;
	INT		level		= ATTR_TO_GLEVEL(DDD_InfoAttr(PARHDR(pv)));
	GRID		*theGrid	= GRID_ON_LEVEL(dddctrl.currMG,level);
	/* TODO: define this static global                    */
	/* TODO: take size as maximum of possible connections */
	size_t	sizeArray[256];
	INT flag;

	PRINTDEBUG(dddif,1,(PFMT " FAMGVectorXferCopy(): v=" VINDEX_FMTX " proc=%d "
		"prio=%d vtype=%d\n",me,VINDEX_PRTX(pv),proc,prio,VTYPE(pv)))
	
    if (DDD_XferWithAddData()) {
	    for(imat=VISTART(pv); imat!=NULL; imat=MNEXT(imat)) {
		    ASSERT(nmat<256);
			sizeArray[nmat++] = MSIZE(imat);
		}
		PRINTDEBUG(dddif,2,(PFMT " FAMGVectorXferCopy(): v=" VINDEX_FMTX 
							" AddData nmat=%d\n",me,VINDEX_PRTX(pv),nmat))
		DDD_XferAddDataX(nmat,TypeMatrix,sizeArray);
	}
	else
	{
		PRINTDEBUG(dddif,2,(PFMT " FAMGVectorXferCopy(): no AddData\n"))
	}
}

void FAMGVectorGatherIMatX (DDD_OBJ obj, int cnt, DDD_TYPE type_id, char **Data)
// derived from dddif/handler.c/FAMGVectorGatherMatX()
{
	VECTOR	*vec = (VECTOR *)obj;
	MATRIX	*imat;
	INT		nmat = 0;

	PRINTDEBUG(dddif,3,(PFMT " FAMGVectorGatherIMatX(): v=" VINDEX_FMTX 
		" VOBJID=%d cnt=%d type=%d veobj=%d vtype=%d\n",
		me,VINDEX_PRTX(vec),ID(VOBJECT(vec)),cnt,type_id,
		OBJT(vec),VTYPE(vec)))

	if (cnt<=0) return;

	for (imat=VISTART((VECTOR *) vec); imat!=NULL; imat=MNEXT(imat))
	{
		int Size;

		IFDEBUG(dddif,0)
		if (cnt<nmat+1)
		{
			PRINTDEBUG(dddif,0,(PFMT " FAMGVectorGatherIMatX(): v=" VINDEX_FMTX 
				" cnt=%d nmat=%d type=%d veobj=%d\n",
				me,VINDEX_PRTX(vec),cnt,nmat,type_id,OBJT(vec)))
			assert(0);
		}
		ENDDEBUG

		Size = MSIZE(imat);
		memcpy(Data[nmat],imat,Size);

		PRINTDEBUG(dddif,3,(PFMT " FAMGVectorGatherIMatX(): v=" VINDEX_FMTX 
			" mat=%x Size=%d vectoID=" VINDEX_FMTX "\n",
			me,VINDEX_PRTX(vec),imat,Size,VINDEX_PRTX(MDEST(imat))))

		nmat++;
	}
}


void FAMGVectorScatterIMatX (DDD_OBJ obj, int cnt, DDD_TYPE type_id, char **Data, int newness)
// derived from dddif/handler.c/FAMGVectorScatterConnX()
{
	VECTOR		*vec		= (VECTOR *)obj;
	CONNECTION	*first		= NULL,
				*last		= NULL;
	INT			level		= ATTR_TO_GLEVEL(DDD_InfoAttr(PARHDR(vec)));
	GRID		*theGrid	= GRID_ON_LEVEL(dddctrl.currMG,level);
	INT			prio 		= PRIO(vec);
	INT			i;
	INT			nconn		= 0;
	INT			newconn		= 0;

	PRINTDEBUG(dddif,3,(PFMT " FAMGVectorScatterIMatX(): v=" VINDEX_FMTX 
		" cnt=%d type=%d veobj=%d vtype=%d\n",\
		me,VINDEX_PRTX(vec),cnt,type_id,OBJT(vec),VTYPE(vec)))
	
	if (cnt<=0) return;

	for (i=0; i<cnt; i++)
	{
		MATRIX *mcopy = (MATRIX *)Data[i];

		/* reset MNEXT */
		MNEXT(mcopy) = NULL;

		if (MDEST(mcopy)==NULL)
		{
			/* destination vector is not on this processor  */
			/* -> matrix entry is useless, throw away       */
			PRINTDEBUG(dddif,4,(PFMT " FAMGVectorScatterIMatX(): v=" VINDEX_FMTX
				" mat=%x Size=%d, USELESS no dest vector \n",
				me,VINDEX_PRTX(vec),mcopy,MSIZE(mcopy)))
			continue;
		}

		{
			MATRIX *m,*mat=NULL;

			/* does matrix entry already exist? */
			/* TODO not nice, linear search, change this! */
			for (m=VISTART((VECTOR *)vec); 
				 m!=NULL && (mat==NULL); 
				 m=MNEXT(m))
			{
				if (MDEST(m)==MDEST(mcopy)) mat=m;
			}

			/* matrix entry is really new */
			if (mat == NULL)
			{
				/* handle diagonal entry */
				if (MDIAG(mcopy))
				{
					/* matrix diagonal entry, no other vector is involved */
					assert(0);
				}
				/* handle off-diagonal entry */
				else
				{
					/* matrix off-diagonal entry, another vector is involved */
					VECTOR *other = MDEST(mcopy);
					MATRIX *m, *back=NULL, *newm;
	
					/* does connection already exist for other vec? no, can not be here for Imats*/
					{
						MATRIX *otherm;
						CONNECTION *conn = (CONNECTION *)
						#ifndef DYNAMIC_MEMORY_ALLOCMODEL
							GetFreelistMemory(dddctrl.currMG->theHeap,
											  MSIZE(mcopy));
						#else
					  		GetMemoryForObject(dddctrl.currMG,MSIZE(mcopy),MAOBJ);
						#endif

						nconn++; newconn++;

						if (conn==NULL)
						{
							ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << me << ":  FAMGVectorScatterIMatX(): can't get mem for mat=" << &mcopy << endl;
							FAMGError(ostr);
							return;
						}
	

						if (MOFFSET(mcopy))
						{
							assert(0);
						}
						else
						{
							newm = (MATRIX *) conn;
							PRINTDEBUG(dddif,4,(PFMT " FAMGVectorScatterIMatX(): v="
							VINDEX_FMTX " conn=%x newm=%x Size=%d vectoID=" 
							VINDEX_FMTX " GETMEM\n",
							me,VINDEX_PRTX(vec),conn,newm, MSIZE(mcopy),
							VINDEX_PRTX(MDEST(mcopy))))
						}
					}

					/* TODO: define clearly     
					memcpy(newm,mcopy,MSIZE(mcopy)); */
					memset(newm,0,MSIZE(mcopy));
					memcpy(newm,mcopy,sizeof(MATRIX)-sizeof(DOUBLE));

					if (first==NULL) first = newm;
					else MNEXT(last) = newm;
					last = newm;
				}
			}
			/* matrix entry does already exist */
			else
			{
				PRINTDEBUG(dddif,4,(PFMT " FAMGVectorScatterIMatX(): v="
					VINDEX_FMTX " mat=%x Size=%d vectoID=" VINDEX_FMTX 
					" FOUND\n",me,VINDEX_PRTX(vec),mat,MSIZE(mcopy),
					VINDEX_PRTX(MDEST(mcopy))))
			}
		}
	}

	/* enter matrix list at the beginning of existing list for this vector */
	// is already done
	
		#ifdef Debug
/*
		{
		MATRIX *mat;
        PRINTDEBUG(dddif,4,(PFMT " FAMGVectorScatterIMatX():  v="
                VINDEX_FMTX "new matrices:\n",me,VINDEX_PRTX(vec)));
		for (mat=first; mat!=NULL; mat=MNEXT(mat))
		{
        	PRINTDEBUG(dddif,4,(PFMT "     mat=%x dest=" EID_FMTX "\n",me,mat,VINDEX_PRTX(MDEST(mat))));
		}
		}
*/
		#endif
}

void FAMGGrid::CommunicateNodeStatus()
// idea: send a copy of the vec to all its already existing copies because then 
//       the I-matrix entries can be copied too (as the stiffness matrix entries
//		 are handled in the original VectorXfer routine).
{
	INT i, size;
	int *proclist;
	VECTOR *vec;
	FAMGNode *nodei;
	
	// set my own handlers
	DDD_SetHandlerXFERCOPY         (TypeVector, FAMGVectorXferCopy);
	DDD_SetHandlerXFERGATHERX      (TypeVector, FAMGVectorGatherIMatX);
	DDD_SetHandlerXFERSCATTERX     (TypeVector, FAMGVectorScatterIMatX);

	DDD_XferBegin();
	
    for(i = 0; i < n; i++)
   	{
       	nodei = graph->GetNode(i);
		vec = ((FAMGugVectorEntryRef*)(nodei->GetVec().GetPointer()))->myvector();
		size = sizeof(VECTOR)-sizeof(DOUBLE)
					  +FMT_S_VEC_TP(MGFORMAT(dddctrl.currMG),VTYPE(vec));
		proclist = DDD_InfoProcList(PARHDR(vec));
		proclist += 2;	// skip the entry for this instance itself
		while( *proclist != -1 )
		{
			DDD_XferCopyObjX(PARHDR(vec), proclist[0], proclist[1], size);
			proclist += 2;
		}

	}
	
	DDD_XferEnd();
	
	// restore the ug handlers
	DDD_SetHandlerXFERCOPY         (TypeVector, VectorXferCopy);
	DDD_SetHandlerXFERGATHERX      (TypeVector, VectorGatherMatX);
	DDD_SetHandlerXFERSCATTERX     (TypeVector, VectorScatterConnX);
}
#endif
