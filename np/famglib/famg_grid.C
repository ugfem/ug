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

#ifdef ModelP
#include "famg_coloring.h"

//#define EXTENDED_OUTPUT_FOR_DEBUG 

// nodes are only useful for debugging to have geometric positions 
// of the vectors but aren't used for calculations.
// In addition there can arise a dangerous situation:
//		on coarser grids there may be Xfer'ed a vector and it's node.
//		This node takes it's NODEVEC with, but due to the abuse of
//		the nodes on AMG levels this NODEVEC may be on level 0.
// 		Thus an already constructed grid gets changed which 
//		leads to inconsistent data structures! 
//		DON'T DO THIS ANYMORE!
// Don't define FAMG_SEND_NODES
//#define FAMG_SEND_NODES
#endif

#ifdef USE_UG_DS
extern "C"
{
#include "gm.h"
#include "ugm.h"
#include "disctools.h" // for AssembleDirichletBoundary
#include "commands.h" // nur zum testen fuer GetCurrentMultigrid()
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

INT l_force_consistence (GRID *g, const VECDATA_DESC *x);
INT l_vector_collectAll (GRID *g, const VECDATA_DESC *x);

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
	
#ifdef ModelP
	if (l_vector_collect(mygrid,((FAMGugVector&)fgdefect).GetUgVecDesc())!=NUM_OK) 
		assert(0);
#endif
	
	// jacobi smoothing for the fine nodes
    fgsolution.JacobiSmoothFG( *GetConsMatrix(), fgdefect );
#ifdef ModelP
	if (l_vector_consistent(mygrid,((FAMGugVector&)fgsolution).GetUgVecDesc())!=NUM_OK) 
		assert(0);
#endif
	
	// correct defect
	Defect();
#ifdef ModelP
	if (l_vector_collect(mygrid,((FAMGugVector&)fgdefect).GetUgVecDesc())!=NUM_OK) 
		assert(0);
#endif

	cgdefect = 0.0;
	
	FAMGVectorIter fiter(GetGridVector());
	while( fiter(fvec) )
#ifdef ModelP
		if(IS_FAMG_MASTER(((FAMGugVectorEntryRef*)(fvec.GetPointer()))->myvector()))
#endif		
	    	for(transfg = transfer.GetFirstEntry(fvec); transfg != NULL; transfg = transfg->GetNext())
				cgdefect[transfg->GetCol()] += transfg->GetRestriction() * fgdefect[fvec];
		
#ifdef ModelP
	if (l_vector_collectAll(DOWNGRID(mygrid),((FAMGugVector&)cgdefect).GetUgVecDesc())!=NUM_OK) 
		assert(0);
#endif
	
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
	
	
    #ifdef ModelP
	// distribute master values to V(H)Ghosts
	if (l_ghostvector_consistent(DOWNGRID(mygrid),((FAMGugVector&)cgsol).GetUgVecDesc())!= NUM_OK)
		assert(0);
	#endif        
	
#ifdef PROTOCOLNUMERIC
	FAMGVector &cgdefect = *(cg->GetVector(FAMGDEFECT));
	cgdefect = 1.0;
	#ifdef ModelP
	FAMGVectorIter citer(cg->GetGridVector());
	FAMGVectorEntry cvec;
	while( citer(cvec) )	// make an inconsistent vector with the const value 1.0
	{
		if(!IS_FAMG_MASTER(((FAMGugVectorEntryRef*)(cvec.GetPointer()))->myvector()))
			cgdefect[cvec] = 0.0;
	}
	#endif
    cout<<"FAMGGrid::Prolongation: 1*cgsol before prolongation "<<cgdefect*cgsol<<endl;
#endif
	
	FAMGVectorIter fiter(GetGridVector());
	while( fiter(fvec) )
	{
#ifdef ModelP
		if(!IS_FAMG_MASTER(((FAMGugVectorEntryRef*)(fvec.GetPointer()))->myvector()))
			continue;
#endif		
		sum = 0.0;
	    for(transfg = transfer.GetFirstEntry(fvec); transfg != NULL; transfg = transfg->GetNext())
        	sum += transfg->GetProlongation() * cgsol[transfg->GetCol()];
		fgsol[fvec] = sum;
	}
#ifdef ModelP
	if (l_force_consistence(mygrid,((FAMGugVector&)fgsol).GetUgVecDesc())!=NUM_OK) 
		assert(0);
#endif

#ifdef PROTOCOLNUMERIC
    cout<<"FAMGGrid::Prolongation: defect after prolongation "<<fgdefect*fgdefect<<endl;
    cout<<"FAMGGrid::Prolongation: sol after prolongation "<<fgsol*fgsol<<endl;
    cout<<"FAMGGrid::Prolongation: defect*sol after prolongation "<<fgdefect*fgsol<<endl;
#endif
	
	// prepare defect for jacobi smoothing
	Defect();
#ifdef ModelP
	if (l_vector_collect(mygrid,((FAMGugVector&)fgdefect).GetUgVecDesc())!=NUM_OK) 
		assert(0);
#endif
	
#ifdef PROTOCOLNUMERIC
    cout<<"FAMGGrid::Prolongation: defect after defect "<<fgdefect*fgdefect<<endl;
    cout<<"FAMGGrid::Prolongation: sol after defect "<<fgsol*fgsol<<endl;
    cout<<"FAMGGrid::Prolongation: defect*sol after defect "<<fgdefect*fgsol<<endl;
#endif

	if(c == NULL)
    {
        // famg as solver

		#ifdef ModelP
		assert(0);	// "FAMG as solver" not ported to parallel
		#endif

        // jacobi smoothing for the fine nodes
        fgsol.JacobiSmoothFG( *GetConsMatrix(), fgdefect );
        
        // correct defect
        Defect();
    }
    else
    {
        // famg as transfer
        (*c) += fgsol;
        fgsol = 0.0;

#ifdef PROTOCOLNUMERIC
	    cout<<"FAMGGrid::Prolongation: sol after update "<<(*c)*(*c)<<endl;
    	cout<<"FAMGGrid::Prolongation: defect*c after update "<<fgdefect*(*c)<<endl;
#endif

		// jacobi smoothing for the fine nodes
        fgsol.JacobiSmoothFG( *GetConsMatrix(), fgdefect );
#ifdef ModelP
		if (l_vector_consistent(mygrid,((FAMGugVector&)fgsol).GetUgVecDesc())!=NUM_OK) 
			assert(0);
#endif

#ifdef PROTOCOLNUMERIC
	    cout<<"FAMGGrid::Prolongation: defect after smoothing "<<fgdefect*fgdefect<<endl;
    	cout<<"FAMGGrid::Prolongation: sol after smoothing "<<fgsol*fgsol<<endl;
	    cout<<"FAMGGrid::Prolongation: defect*sol after smoothing "<<fgdefect*fgsol<<endl;
#endif

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
	GetVector(FAMGUNKNOWN)->dampedJacobiSmoother( *GetConsMatrix(), *GetVector(FAMGDEFECT) );
    return;
}

void FAMGGrid::FGSSmooth()
{
	GetVector(FAMGUNKNOWN)->FGSSmoother( *GetConsMatrix(), *GetVector(FAMGDEFECT) );
    return;
}

void FAMGGrid::BGSSmooth()
{
	GetVector(FAMGUNKNOWN)->BGSSmoother( *GetConsMatrix(), *GetVector(FAMGDEFECT) );
    return;
}

void FAMGGrid::SGSSmooth()
{
	GetVector(FAMGUNKNOWN)->SGSSmoother( *GetConsMatrix(), *GetVector(FAMGDEFECT) );
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
	FAMGMatrixAlg &mat = *Consmatrix;
	FAMGMatrixEntry me;
	FAMGVectorEntry ve;
    int k;

    const int stv = FAMGGetParameter()->Getstv();

    for(k = 0; k < stv; k++)
    {
		assert(0); // TODO: noch ueberlegen ob matrix Consmatrix richtig ist
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
	k = mat.GetN();
#ifdef ModelP
	k = UG_GlobalSumINT( k );
#endif
    scaling = sqrt((double)k) / normA;
	
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
#ifdef ModelP
	ostr << me << ": ";
#endif
	ostr << "unknowns: " << nn << "\t";
	ostr << "avg. stencil: " << (double)nl/(double)nn;
#ifdef ModelP
	ostr << " master= "<<GetNrMasterVectors()<<" border= "<<GetNrBorderVectors()<<" ghosts= "<<GetNrGhostVectors();
#endif
	ostr << endl;
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
        ostr << __FILE__ << ", line " << __LINE__ << ": cgsmoother = ilut" << endl;
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
        ostr << __FILE__ << ", line " << __LINE__ <<  ": presmoother = fgs" << endl;
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
        ostr << __FILE__ << ", line " <<  __LINE__ <<  "postsmoother = bgs" << endl;
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
	#ifdef ModelP
	GetNrMasterVectors() = -1;		// default
	GetNrBorderVectors() = -1;		// default
	GetNrGhostVectors() = -1;		// default
	#endif
#else
    SetFather(NULL);
#endif
	
    matrix = system.GetMatrix();
    Consmatrix = system.GetConsMatrix();
    tmpmatrix = Consmatrix;
	
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
	if(grid_pattern.GetMatrix()!=grid_pattern.GetConsMatrix())
	{
		Consmatrix = new FAMGugMatrix(new_grid, *(FAMGugMatrix*)grid_pattern.GetConsMatrix());
	    if(Consmatrix == NULL)
			return 1;
	}
	else
		Consmatrix = matrix;
	
	#ifdef ModelP
	GetNrMasterVectors() = -1;		// default
	GetNrBorderVectors() = -1;		// default
	GetNrGhostVectors() = -1;		// default
	#endif
#else	
    matrix = (FAMGMatrix *) FAMGGetMem(sizeof(FAMGMatrix),FAMG_FROM_TOP);
    if(matrix == NULL)
		return 1;
    if(matrix->Init(n))
		return 1;
	mygridvector = NULL; // aendern
#endif
	
    tmpmatrix = Consmatrix;
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
	if( GetConsMatrix()!=GetMatrix() )
	{
		delete GetConsMatrix();
		SetConsMatrix(NULL);
	}
	delete GetMatrix();
	SetMatrix(NULL);
#else
	assert(0);
    father = NULL;
    tmpmatrix = GetConsMatrix();
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
        ostr << __FILE__ << ", line " <<  __LINE__  << ": can not bend transfer entries to coarse grid" << endl;
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
        ostr << __FILE__ << ", line " <<  __LINE__  << ": j="<<j<<" cg->N="<<GetN()<<" number of coarse grid node doesn't match" << endl;
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
	double TotalTime;

#ifdef ModelP
	double GraphColorTime, BorderTime;
	
	// set a consistent copy of the stiffmat
	// in order to be able to construct a consistent interpolationmatrix
	
	MULTIGRID *mg = MYMG(mygrid);
	INT level = GLEVEL(mygrid);
	MATDATA_DESC *A, *ACons;
	
	A = ((FAMGugMatrix*)GetMatrix())->GetMatDesc();
	ACons = ((FAMGugMatrix*)GetConsMatrix())->GetMatDesc();
	
    if (dmatcopy(mg,level,level,ALL_VECTORS,ACons,A) != NUM_OK) 
	   assert(0);
	if (l_matrix_consistent(mygrid,ACons,MAT_CONS) != NUM_OK) 
	   assert(0);
#endif

	TotalTime = CURRENT_TIME_LONG;
	
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
	int BorderCycles;
	
	GraphColorTime = CURRENT_TIME;
	
	if( ConstructColoringGraph(GRID_ATTR(mygrid)) )
	{
		cout << "FAMGGrid::ConstructTransfer(): ERROR in constructing the coloring graph"<<endl<<fflush;
		FAMGReleaseHeap(FAMG_FROM_BOTTOM);
		return 1;
	}
		
	if( ConstructColoring( FAMGGetParameter()->GetColoringMethod() ) )
	{
		cout << "FAMGGrid::ConstructTransfer(): ERROR in coloring the graph"<<endl<<fflush;
		FAMGReleaseHeap(FAMG_FROM_BOTTOM);
		return 1;
	}
		
	GraphColorTime = CURRENT_TIME - GraphColorTime;
	BorderCycles = UG_GlobalMaxINT((int)FAMGMyColor);
	cout <<me<<": level = "<<level<<" my color = "<<FAMGMyColor<<" max. color = "<<BorderCycles<<" ColoringMethod = "<<FAMGGetParameter()->GetColoringMethod()<<endl;
	
	BorderTime = CURRENT_TIME_LONG;
	for ( int color = 0; color <= BorderCycles; color++)
	{
		if( FAMGMyColor == color )
		{
		    if (graph->InsertHelplist()) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
		    if (graph->EliminateNodes(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
    		if (graph->RemainingNodes()) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
		}
		
//prim(GLEVEL(GetugGrid()));//?????????????????????????????????????????????
		CommunicateNodeStatus();
	}	
	BorderTime = CURRENT_TIME_LONG - BorderTime;
	
	// put the undecided nodes from the core partition into the list
    for(i = 0; i < n; i++)
   	{
       	nodei = graph->GetNode(i);
	
		// now only the unmarked nodes must be inserted into the list
		if( nodei->IsUndecidedNode() )
		{
			vec = ((FAMGugVectorEntryRef*)(nodei->GetVec().GetPointer()))->myvector();
			if( IS_FAMG_MASTER(vec) ) // only master vectors can be in border the of the core partition
				if(graph->InsertNode(this, nodei))
				{
					FAMGReleaseHeap(FAMG_FROM_BOTTOM);
					return 1;
				}
		}
    }
#endif

#ifdef SIMULATE_HALFENING	// TODO: remove it
	FAMGNode *nodei;
	VECTOR *vec;
	
	if (graph->InsertHelplist()) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
	if (graph->EliminateNodes(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
	if (graph->RemainingNodes()) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
	// put the undecided nodes from the core partition into the list
    for(i = 0; i < n; i++)
   	{
       	nodei = graph->GetNode(i);
	
		// now only the unmarked nodes must be inserted into the list
		if( nodei->IsUndecidedNode() )
		{
			vec = ((FAMGugVectorEntryRef*)(nodei->GetVec().GetPointer()))->myvector();
				if(graph->InsertNode(this, nodei))
				{
					FAMGReleaseHeap(FAMG_FROM_BOTTOM);
					return 1;
				}
		}
    }
#endif
	
	if( graph->InsertHelplist() ) {FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
    if (graph->EliminateNodes(this)) {FAMGReleaseHeap(FAMG_FROM_BOTTOM); return 1;}
    
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
//prim(GLEVEL(GetugGrid()));//?????????????????????????????????????????????
	CommunicateNodeStatus();
#endif
	
    nf = graph->GetNF();
    // TODO: not neceassary any more: matrix->MarkUnknowns(graph);
 
	TotalTime = CURRENT_TIME_LONG - TotalTime;
	cout
#ifdef ModelP
		<< me << ": "
#endif
		<< "level "<< level 
	    << " FAMG time = " << TotalTime
#ifdef ModelP
	    << " ( graph coloring = " << GraphColorTime
	    << " border = " << BorderTime << " )"
#endif
		<< endl;
		
//prim(GLEVEL(GetugGrid()));//?????????????????????????????????????????????

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


///////////////////////////////////////////////////////////////////////////////
// ConstructOverlap
///////////////////////////////////////////////////////////////////////////////

#ifdef ModelP

static inline void TransferVector( VECTOR *v, DDD_PROC dest_pe, DDD_PRIO dest_prio, int size )
{	
	DDD_XferCopyObjX(PARHDR(v), dest_pe, dest_prio, size);
	
#ifdef FAMG_SEND_NODES
	/* for debugging purpose send the node too */
	/* TODO: be sure that this isn't done for production runs */
	IFDEBUG(np,0)
	GEOM_OBJECT *obj = VOBJECT(v);
	if( obj != NULL )
		if (VOTYPE(v) == NODEVEC)
		{
			PRINTDEBUG(dddif,2,(PFMT " TransferVector(): vec= " VINDEX_FMTX" n=" ID_FMTX " Xfer v=" 
				VID_FMTX "\n",me,VINDEX_PRTX(v),ID_PRTX((NODE*)obj),VID_PRTX(MYVERTEX((NODE*)obj))))
			DDD_XferCopyObj(PARHDR((NODE*)obj), dest_pe, dest_prio);
		}
		else
			assert(0); /* unrecognized vector type; extend code if necessary */	
	ENDDEBUG
#endif
}

static int SendToMaster( DDD_OBJ obj)
{
	VECTOR *vec = (VECTOR *)obj, *w;
	MATRIX *mat;
	DDD_PROC masterPe;
	int *proclist, i, found, size;
	
	if( IS_FAMG_MASTER(vec) )
		return 0;		// we want only border vectors here
	
	PRINTDEBUG(np,1,("%d: SendToMaster: "VINDEX_FMTX"\n",me,VINDEX_PRTX(vec)));
	
	for( mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat) )
		if( IS_FAMG_MASTER(MDEST(mat)) )
			break;	// now vec is in overlap1 without core partition
	
	if( mat==NULL )
		return 0;	// vec has no master neighbor; hence is not in overlap1
	
	proclist = DDD_InfoProcList(PARHDR(vec));
	masterPe = (DDD_PROC)me;	// init with an unpossible value
	for( i=0; proclist[i]!=-1; i+=2 )
		if( proclist[i+1] == PrioMaster )
		{
			masterPe = proclist[i];
			break;
		}
	assert(masterPe!=(DDD_PROC)me);	// a master copy must exist on an other processor

	// It seems to be necessary to send also the vec itself to let new 
	// connection be constructed from both (vec and its neighbor)
	PRINTDEBUG(np,1,("%d: SendToMaster %d: myself "VINDEX_FMTX"\n",me,masterPe,VINDEX_PRTX(vec)));
	size = sizeof(VECTOR)-sizeof(DOUBLE)+FMT_S_VEC_TP(MGFORMAT(dddctrl.currMG),VTYPE(vec));
	TransferVector(vec, masterPe, PrioMaster, size);
	
	for( mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat) )
	{
		w = MDEST(mat);
		
		// search whether w has a copy on masterPe
		proclist = DDD_InfoProcList(PARHDR(w));
		found = FALSE;
		for( i=0; proclist[i]!=-1; i+=2 )
			if( masterPe == proclist[i] )
			{
				found = TRUE;
				break;
			}
		
		if( !found )	// if it has no copy send w to masterPe
		{
			PRINTDEBUG(np,1,("%d: SendToMaster %d:     -> "VINDEX_FMTX"\n",me,masterPe,VINDEX_PRTX(w)));
			size = sizeof(VECTOR)-sizeof(DOUBLE)+FMT_S_VEC_TP(MGFORMAT(dddctrl.currMG),VTYPE(w));
			TransferVector(w, masterPe, PrioBorder, size);
		}
		else
			PRINTDEBUG(np,1,("%d: SendToMaster %d:     is "VINDEX_FMTX"\n",me,masterPe,VINDEX_PRTX(w)));
	}
	
	return 0;
}

static int SendToOverlap1( DDD_OBJ obj)
// every master sends itself to each processor where a neighbor has a border copy
{
	VECTOR *vec = (VECTOR *)obj, *w, *wn;
	MATRIX *mat, *matw;
	int *proclist_vec, *proclist_w, i, found, size;
	
	if( !IS_FAMG_MASTER(vec) )
		return 0;		// we want only border vectors here

	PRINTDEBUG(np,1,("%d: SendToOverlap1: "VINDEX_FMTX"\n",me,VINDEX_PRTX(vec)));

	proclist_vec = DDD_InfoProcList(PARHDR(vec));
	
	for( mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat) )
	{
		w = MDEST(mat);
		
		proclist_w = DDD_InfoProcList(PARHDR(w));
		for( ; proclist_w[0]!=-1; proclist_w += 2 )
		{		
			// send only to all not-master copies
			if( proclist_w[1] == PrioMaster )
				continue;
			
			// don't send info to myself
			if(proclist_w[0]==me)
				continue;
			
			// search whether vec has already a copy on the PE of the neighbor border copy
			found = FALSE;
			for( i=0; proclist_vec[i]!=-1; i+=2 )
				if( proclist_w[0] == proclist_vec[i] )
				{
					found = TRUE;
					break;
				}
		
		//temp weg if( !found )	// if it has no copy send vec
			{
				PRINTDEBUG(np,1,("%d: SendToOverlap1 %d:     -> "VINDEX_FMTX"\n",me,proclist_w[0],VINDEX_PRTX(w)));
				size = sizeof(VECTOR)-sizeof(DOUBLE)+FMT_S_VEC_TP(MGFORMAT(dddctrl.currMG),VTYPE(vec));
				TransferVector(vec, proclist_w[0], PrioBorder, size);
				// s.u. size = sizeof(VECTOR)-sizeof(DOUBLE)+FMT_S_VEC_TP(MGFORMAT(dddctrl.currMG),VTYPE(w));
				// wird unten mit gemacht DDD_XferCopyObjX(PARHDR(w), proclist_w[0], proclist_w[1], size);
			}

			for( matw=VSTART(w); matw!=NULL; matw=MNEXT(matw) )
			{
				PRINTDEBUG(np,1,("%d: SendToOverlap1 %d:     -> "VINDEX_FMTX"\n",me,proclist_w[0],VINDEX_PRTX(wn)));
				wn = MDEST(matw);
                size = sizeof(VECTOR)-sizeof(DOUBLE)+FMT_S_VEC_TP(MGFORMAT(dddctrl.currMG),VTYPE(wn));
				TransferVector(wn, proclist_w[0], PrioBorder, size);
			}
		}
	}

	return 0;
}

void FAMGGrid::ConstructOverlap()
// extend the overlap as far as necessary; at least 2 links deep
// the vectorlist will be renumbered
{
	VECTOR *vec, *mv;
	INT i;
	FAMGMatrixAlg *matrix_tmp;
	
	if(GLEVEL(mygrid)==0)
	{	
		/* change ghosts to border */
		DDD_XferBegin();
		for( vec=PFIRSTVECTOR(mygrid); PRIO(vec) != PrioBorder && PRIO(vec) != PrioMaster; vec = SUCCVC(vec))
		{
			#ifdef FAMG_SEND_NODES
			DDD_XferPrioChange( PARHDR((NODE*)VOBJECT(vec)), PrioBorder ); /* TODO: cancel this line; its only for beauty in checks */
			#endif
			DDD_XferPrioChange( PARHDR(vec), PrioBorder );
		}
		/* elements with old ghostprios cause errors in check; but that's ok */
		DDD_XferEnd();
		/* vectors which just become border have not set the VECSKIP flag correctly! 
		   this flags will be corrected by the subsequent a_vector_vecskip */

ASSERT(!DDD_ConsCheck());	
	}	

	DDD_XferBegin();
		DDD_IFAExecLocal( BorderVectorIF, GRID_ATTR(mygrid), SendToMaster );
	DDD_XferEnd();
	
	DDD_XferBegin();
		DDD_IFAExecLocal( BorderVectorIF, GRID_ATTR(mygrid), SendToOverlap1 );
	DDD_XferEnd();
	

	// count & set number of vectors
	mv = FIRSTVECTOR(mygrid);
	i = 0;
	GetNrMasterVectors()=0;
	GetNrBorderVectors()=0;
	GetNrGhostVectors()=0;
	for( vec=PFIRSTVECTOR(mygrid); vec!=mv; vec=SUCCVC(vec))
	{
		VINDEX(vec) = i++;
		GetNrGhostVectors()++;
	}
	for( ; vec!=NULL; vec=SUCCVC(vec))
	{
		VINDEX(vec) = i++;
		if( IS_FAMG_MASTER(vec) )
			GetNrMasterVectors()++;
		else
			GetNrBorderVectors()++;
	}
	
	n = NVEC(mygrid);
	assert(i==n);	// otherwise the vectorlist became inconsistent
	assert(i==(GetNrMasterVectors()+GetNrBorderVectors()+GetNrGhostVectors()));	// otherwise the vectorlist became inconsistent
	
	// correct the number of vectors in matrices
	matrix_tmp = GetMatrix();
	if( matrix_tmp!=NULL)
		matrix_tmp->GetN() = n;
	matrix_tmp = GetTmpMatrix();
	if( matrix_tmp!=NULL)
		matrix_tmp->GetN() = n;
	matrix_tmp = GetConsMatrix();
	if( matrix_tmp!=NULL)
		matrix_tmp->GetN() = n;
	
	
	if(GLEVEL(mygrid)==0)
	{		
		// do modifications for dirichlet vectors (as coarsegrid solver)
		// communicate vecskip flags and dirichlet values
		if (a_vector_vecskip(MYMG(mygrid),0,0,((FAMGugVector*)GetVector(FAMGUNKNOWN))->GetUgVecDesc()) != NUM_OK)
			abort();
		// set dirichlet modification in the matrix
		if (AssembleDirichletBoundary (mygrid,((FAMGugMatrix*)GetMatrix())->GetMatDesc(),((FAMGugVector*)GetVector(FAMGUNKNOWN))->GetUgVecDesc(),((FAMGugVector*)GetVector(FAMGDEFECT))->GetUgVecDesc()))
			abort();
	}	

}
#endif


///////////////////////////////////////////////////////////////////////////////
// CommunicateNodeStatus
///////////////////////////////////////////////////////////////////////////////

#ifdef ModelPQQQQQQQQQQQQQQQQQQQQQ

/*
	Das Problem bei dieser Methode ist zur Zeit, dass beim Zielknoten VINDEX und 
	VCCOARSE ueberschrieben werden, bevor noch geeignete Ueberpruefungen (VCCOARSE)
	oder Verhinderungen (VINDEX) stattfinden koennen. 
	Wenn VINDEX richtig waere, koennte man uden zugehoerigen Knoten fingen und so 
	das Problem mit VCCOARSE umgehen. VINDEX geht aber verloren.
	Man koennte hinterher den Index (da die Verkettung unveraendert geblieben ist)
	wieder richtig setzen, dann ist es aber zu spaet, weil die geschickten 
	I-Matrizen nur waehrend des Haenfdleraufrufs fuer den einzelnen
	Knoten verfuegbar sind.
	So ist das Vorgehen also nicht praktikabel.
	
	Im derzeit vorhandenen Code fehlen auch noch die Aenderungsarbeiten an 
	der FAMG Datenstruktur (Nodes, Listen).
*/
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


#ifdef ModelP

/*
	2. Version
*/
static FAMGGraph *Communication_Graph = NULL;
static FAMGGrid *Communication_Grid = NULL;
enum FAMG_MSG_TYPE {FAMG_TYPE_NONE, FAMG_TYPE_COARSE, FAMG_TYPE_FINE};

#ifdef EXTENDED_OUTPUT_FOR_DEBUG
static int Gather_NodeStatus (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
#else
static int Gather_NodeStatus (DDD_OBJ obj, void *data)
#endif
{
	VECTOR *vec = (VECTOR *)obj;
	char *buffer = (char*)data;
	FAMGNode *node = Communication_Graph->GetNode(VINDEX(vec));
	FAMG_MSG_TYPE &msgtype = *(FAMG_MSG_TYPE*)buffer;

	if( node->GetFlagNewMarked() )
	{
		if( node->IsFGNode() )
		{
			#ifdef EXTENDED_OUTPUT_FOR_DEBUG
			PRINTDEBUG(np,1,("%d: Gather_NodeStatus: Fine "VINDEX_FMTX" for PE %d with prio %d\n",me,VINDEX_PRTX(vec),proc, prio));
			#else
			PRINTDEBUG(np,1,("%d: Gather_NodeStatus: Fine "VINDEX_FMTX"\n",me,VINDEX_PRTX(vec)));
			#endif
		
			msgtype = FAMG_TYPE_FINE;
			buffer += sizeof(msgtype);
			
			int &NumberEntries = *(int*)buffer;
			NumberEntries = 0;
			
			// prepare array offsets into the buffer
			buffer = (char*)data + CEIL(sizeof(msgtype)+sizeof(int));	// round up to achive alignment
			DDD_GID *pgid = (DDD_GID*)buffer;
			DOUBLE *prolongation = (DOUBLE*)(buffer+CEIL(FAMGMAXPARENTS*sizeof(DDD_GID)));
			DOUBLE *restriction = prolongation+FAMGMAXPARENTS;
			FAMGTransferEntry *trans;
			
			for ( MATRIX *imat=VISTART(vec); imat!= NULL; imat = MNEXT(imat) )
			{
				assert(NumberEntries<FAMGMAXPARENTS);
				pgid[NumberEntries] = DDD_InfoGlobalId(PARHDR(MDEST(imat)));
				
				trans = (FAMGTransferEntry*)imat;	// dirty cast
				
				prolongation[NumberEntries] = trans->GetProlongation();
				restriction[NumberEntries] = trans->GetRestriction();
		    	PRINTDEBUG(np,1,("%d: Gather_NodeStatus:      ->"VINDEX_FMTX" P=%g R=%g\n",me,VINDEX_PRTX(MDEST(imat)),prolongation[NumberEntries],restriction[NumberEntries] ));

				NumberEntries++;
			}
		}
		else
		{
			#ifdef EXTENDED_OUTPUT_FOR_DEBUG
			PRINTDEBUG(np,1,("%d: Gather_NodeStatus: Coarse "VINDEX_FMTX" for PE %d with prio %d\n",me,VINDEX_PRTX(vec),proc,prio));
			#else
			PRINTDEBUG(np,1,("%d: Gather_NodeStatus: Coarse "VINDEX_FMTX"\n",me,VINDEX_PRTX(vec)));
			#endif

			assert(node->IsCGNode());
			msgtype = FAMG_TYPE_COARSE;
		}
		
		// don't reset SetFlagNewMarked here, because this vector may be called for 
		// further interfaces and needs the NewMarked flag	
	}
	else
	{	// the F/C status of the node has not changed
		#ifdef EXTENDED_OUTPUT_FOR_DEBUG
		PRINTDEBUG(np,1,("%d: Gather_NodeStatus: not changed "VINDEX_FMTX" for PE %d with prio %d\n",me,VINDEX_PRTX(vec),proc,prio));
		#else
		PRINTDEBUG(np,1,("%d: Gather_NodeStatus: not changed "VINDEX_FMTX"\n",me,VINDEX_PRTX(vec)));
		#endif
		
		msgtype = FAMG_TYPE_NONE;
	}
	
	return 0;
}

#ifdef EXTENDED_OUTPUT_FOR_DEBUG
static int Scatter_NodeStatus (DDD_OBJ obj, void *data, DDD_PROC proc, DDD_PRIO prio)
#else
static int Scatter_NodeStatus (DDD_OBJ obj, void *data)
#endif
{
	VECTOR *vec = (VECTOR *)obj;
	char *buffer = (char*)data;
	FAMG_MSG_TYPE msgtype = *(FAMG_MSG_TYPE*)buffer;
	buffer += sizeof(FAMG_MSG_TYPE);

	if( msgtype == FAMG_TYPE_NONE )
	{
		#ifdef EXTENDED_OUTPUT_FOR_DEBUG
		PRINTDEBUG(np,1,("%d: Scatter_NodeStatus: not changed from PE %d "VINDEX_FMTX"\n",me,proc,VINDEX_PRTX(vec)));
		#else
		PRINTDEBUG(np,1,("%d: Scatter_NodeStatus: not changed "VINDEX_FMTX"\n",me,VINDEX_PRTX(vec)));
		#endif
	
		return 0;		// no info for me
	}
	
	FAMGNode *node = Communication_Graph->GetNode(VINDEX(vec));
	
#ifdef Debug
	if( !node->IsUndecidedNode() )
	{
		// check that the state of the node is consistent with the message
		if( msgtype == FAMG_TYPE_COARSE )
		{
			#ifdef EXTENDED_OUTPUT_FOR_DEBUG
			PRINTDEBUG(np,1,("%d: Scatter_NodeStatus: Coarse for already coarse "VINDEX_FMTX" from PE %d with prio %d\n",me,VINDEX_PRTX(vec),proc,prio));
			#else
			PRINTDEBUG(np,1,("%d: Scatter_NodeStatus: Coarse for already coarse "VINDEX_FMTX"\n",me,VINDEX_PRTX(vec)));
			#endif
			
			assert(node->IsCGNode());
		}
		else if( msgtype == FAMG_TYPE_FINE )
		{
			#ifdef EXTENDED_OUTPUT_FOR_DEBUG
			PRINTDEBUG(np,1,("%d: Scatter_NodeStatus: Fine for already fine "VINDEX_FMTX" from PE %d with prio %d\n",me,VINDEX_PRTX(vec),proc,prio));
			#else
			PRINTDEBUG(np,1,("%d: Scatter_NodeStatus: Fine for already fine "VINDEX_FMTX"\n",me,VINDEX_PRTX(vec)));
			#endif
		
			assert(node->IsFGNode());
			// check wether the parents are identical
			MATRIX *imat;
			DDD_GID pgid;
			INT i, nnb, np = *(int*)buffer;		// fetch number of parents from buffer
			// prepare array offsets into the buffer
			DDD_GID *buffer_gid  = (DDD_GID *)((char*)data + CEIL(sizeof(msgtype)+sizeof(int))); // round up to achive alignment
			
			for( nnb=0,imat=VISTART(vec); imat!=NULL; nnb++,imat=MNEXT(imat) )
			{
				pgid = DDD_InfoGlobalId(PARHDR(MDEST(imat)));
				// search the gid in the buffer
				for( i=0; i<np; i++ )
					if( pgid == buffer_gid[i] )
						break;
				if(i>=np)
				{
				    PRINTDEBUG(np,1,("%d: Scatter_NodeStatus: ERROR vector "VINDEX_FMTX" has interpolation matrix to vec "VINDEX_FMTX", but it was not in message\n",me,VINDEX_PRTX(vec),VINDEX_PRTX(MDEST(imat))));
				}
			}
			assert(nnb==np);	// else there was a parent node in the buffer, but not in the I-matrix
		}	
		else
		{
			assert(0); // unrecognized message type
			return 1;
		}
	}
	else
	{
#endif
	if( msgtype == FAMG_TYPE_COARSE )
	{
		#ifdef EXTENDED_OUTPUT_FOR_DEBUG
		PRINTDEBUG(np,1,("%d: Scatter_NodeStatus: Coarse "VINDEX_FMTX" from PE %d with prio %d\n",me,VINDEX_PRTX(vec),proc,prio));
		#else
		PRINTDEBUG(np,1,("%d: Scatter_NodeStatus: Coarse "VINDEX_FMTX"\n",me,VINDEX_PRTX(vec)));
		#endif

		// code from PaList::MarkParents
        Communication_Graph->Remove(node);
        Communication_Graph->MarkCGNode(node);
        Communication_Graph->UpdateNSons(NULL,node->GetPaList(),Communication_Grid);
        Communication_Graph->ClearPaList(node->GetPaList());
        node->SetPaList(NULL);
        Communication_Grid->UpdateNeighborsCG(node->GetId());
	}
	else if( msgtype == FAMG_TYPE_FINE )
	{
		DOUBLE prolongations[FAMGMAXPARENTS], restrictions[FAMGMAXPARENTS];
		int np, parents[FAMGMAXPARENTS], bp, pos = 0;
		DDD_GID pgid;
		MATRIX *mat;	
		
		np = *(int*)buffer;		// fetch number of parents from buffer
		
		#ifdef EXTENDED_OUTPUT_FOR_DEBUG
	    PRINTDEBUG(np,1,("%d: Scatter_NodeStatus: Fine "VINDEX_FMTX" %d parents from PE %d with prio %d\n",me,VINDEX_PRTX(vec), np, proc,prio));
		#else
	    PRINTDEBUG(np,1,("%d: Scatter_NodeStatus: Fine "VINDEX_FMTX" %d parents\n",me,VINDEX_PRTX(vec), np));
		#endif
		
		// prepare array offsets into the buffer
		DDD_GID *buffer_gid  = (DDD_GID *)((char*)data + CEIL(sizeof(msgtype)+sizeof(int))); // round up to achive alignment
		DOUBLE *buffer_prolo = (DOUBLE *)((char*)buffer_gid+CEIL(FAMGMAXPARENTS*sizeof(DDD_GID)));
		DOUBLE *buffer_restr = buffer_prolo+FAMGMAXPARENTS;

		for( bp=0; bp<np; bp++)	// process parent nodes
		{
			pgid = buffer_gid[bp];	// fetch GID of a parent node from buffer
			
		    PRINTDEBUG(np,1,("%d: Scatter_NodeStatus:     -> GID %08x P=%g R=%g\n",me,pgid,buffer_prolo[bp],buffer_restr[bp]));
			for( mat = MNEXT(VSTART(vec)); mat!=NULL; mat=MNEXT(mat) )
				if( DDD_InfoGlobalId(PARHDR(MDEST(mat)))==pgid )
					break;
			if( mat==NULL )
			{
				continue; // this parent node is not on this processor
				// is skippping OK  ?????????
			}
			
			assert( pos<FAMGMAXPARENTS );
			
			parents[pos] = VINDEX(MDEST(mat));
			prolongations[pos] = buffer_prolo[bp];	// TODO: avoid the almost unneccessary copying of the DOUBLEs but consider: pos and bp may be different!
			restrictions[pos] = buffer_restr[bp];
		    PRINTDEBUG(np,1,("%d: Scatter_NodeStatus:     -> "VINDEX_FMTX" P=%g R=%g\n",me,VINDEX_PRTX(MDEST(mat)),prolongations[pos],restrictions[pos]));
			pos++;
		}

		// put the values for the parent nodes into the node
		FAMGPaList *palist = NULL;
		if(Communication_Graph->SavePaList(palist,pos,parents,prolongations,restrictions,1.0))
		{
			assert(0);
			return 1;
		}
		Communication_Graph->ClearPaList(node->GetPaList());		
		node->SetPaList(palist);
		
		// process node as fine node
		node->Eliminate(Communication_Grid);
		node->UpdateNeighborsFG(Communication_Grid);
	}
	else
	{
		assert(0); // unrecognized message type
		return 1;
	}
#ifdef Debug
	}	// end of else
#endif
		
	// call Local_ClearNodeFlag afterwards!
	
	return 0;
}

static int Local_ClearNodeFlag (DDD_OBJ obj)
// clear NewMarked Flag in all interface nodes and their "parents"
{
	VECTOR *vec = (VECTOR *)obj;
	FAMGTransferEntry *trans;
	
	FAMGNode *node = Communication_Graph->GetNode(VINDEX(vec));
	
	node->SetFlagNewMarked(0);	// not newly set but only updated
	if( node->IsFGNode() )
		for( trans=Communication_Grid->GetTransfer()->GetFirstEntry(node->GetVec()); trans!=NULL; trans=trans->GetNext())
			Communication_Graph->GetNode(trans->GetCol())->SetFlagNewMarked(0);	// not newly set but only updated

	return 0;
}


static char text_print[100];
static int Local_Print(DDD_OBJ obj)
// nur zum Testen drin
{
	VECTOR *vec = (VECTOR *)obj;
	
	printf("%d Local_Print: %s %x "VINDEX_FMTX"\n",me,text_print,vec,VINDEX_PRTX(vec));
	return 0;
}

int PrintLocal(char *text)
{
	if (DDD_ConsCheck()) assert(0);
	strcpy(text_print,text);
	DDD_IFAExecLocal( BorderVectorSymmIF, GRID_ATTR(GRID_ON_LEVEL(GetCurrentMultigrid(),0)), Local_Print );
	strcpy(text_print,"OuterVectorSymmIF");
	DDD_IFAExecLocal( OuterVectorSymmIF, GRID_ATTR(GRID_ON_LEVEL(GetCurrentMultigrid(),0)), Local_Print );
	return 0;
}


void FAMGGrid::CommunicateNodeStatus()
{
	PRINTDEBUG(np,1,("%d: FAMGGrid::CommunicateNodeStatus\n",me));
	Communication_Graph = GetGraph();	// set global variable to pass the graph to the Handlers
	Communication_Grid = this;		// set global variable to pass the grid to the Handlers
	
	int size = CEIL(sizeof(FAMG_MSG_TYPE) + sizeof(char)) + 
				CEIL(FAMGMAXPARENTS * ( sizeof(DDD_GID)) + 2 * FAMGMAXPARENTS * sizeof(DOUBLE) );
	size = CEIL(size);
		
	#ifdef EXTENDED_OUTPUT_FOR_DEBUG
	DDD_IFAOnewayX( BorderVectorSymmIF, GRID_ATTR(mygrid), IF_FORWARD, size,
					Gather_NodeStatus, Scatter_NodeStatus);
	#else
	DDD_IFAOneway( BorderVectorSymmIF, GRID_ATTR(mygrid), IF_FORWARD, size,
					Gather_NodeStatus, Scatter_NodeStatus);
	#endif
	
	DDD_IFAExecLocal( BorderVectorSymmIF, GRID_ATTR(mygrid), Local_ClearNodeFlag );
}

///////////////////////////////////////////////////////////////////////////////
// l_force_consistence
///////////////////////////////////////////////////////////////////////////////

static VECDATA_DESC *ConsVector = NULL;

static int Gather_VectorComp (DDD_OBJ obj, void *data)
// exactly copy of Gather_VectorComp
{
	VECTOR *pv = (VECTOR *)obj;
	INT i,type;
	const SHORT *Comp;	

	if (VD_IS_SCALAR(ConsVector)) {
		if (VD_SCALTYPEMASK(ConsVector) & VDATATYPE(pv))
		  *((DOUBLE *)data) = VVALUE(pv,VD_SCALCMP(ConsVector));

		return (NUM_OK);
	}
   
	type = VTYPE(pv);
	Comp = VD_CMPPTR_OF_TYPE(ConsVector,type);
	for (i=0; i<VD_NCMPS_IN_TYPE(ConsVector,type); i++)
		((DOUBLE *)data)[i] = VVALUE(pv,Comp[i]);

	return (NUM_OK);
}

static int Scatter_VectorComp (DDD_OBJ obj, void *data)
// exactly copy of Scatter_VectorComp
{
	VECTOR *pv = (VECTOR *)obj;
	INT i,type,vecskip;
	const SHORT *Comp;	

	if (VD_IS_SCALAR(ConsVector)) {
  	    if (VD_SCALTYPEMASK(ConsVector) & VDATATYPE(pv))
		    if (!VECSKIP(pv))
			    VVALUE(pv,VD_SCALCMP(ConsVector)) += *((DOUBLE *)data);

		return (NUM_OK);
	}

	type = VTYPE(pv);
	Comp = VD_CMPPTR_OF_TYPE(ConsVector,type);
	vecskip = VECSKIP(pv);
	if (vecskip == 0)
		for (i=0; i<VD_NCMPS_IN_TYPE(ConsVector,type); i++)
			VVALUE(pv,Comp[i]) += ((DOUBLE *)data)[i]; 
	else
		for (i=0; i<VD_NCMPS_IN_TYPE(ConsVector,type); i++)
			if (!(vecskip & (1<<i)))				
				VVALUE(pv,Comp[i]) += ((DOUBLE *)data)[i]; 

	return (NUM_OK);
}


INT l_force_consistence (GRID *g, const VECDATA_DESC *x)
// copy values from master to border
{
    INT tp,m; 

    ConsVector = (VECDATA_DESC *)x;

	m = 0;
	for (tp=0; tp<NVECTYPES; tp++)
 	    m = MAX(m,VD_NCMPS_IN_TYPE(ConsVector,tp));

	DDD_IFAOneway(BorderVectorIF, GRID_ATTR(g), IF_BACKWARD, m * sizeof(DOUBLE),
				  Gather_VectorComp, Scatter_VectorComp);
	// TODO: perhaps this communication to the ghosts can be omitted
	DDD_IFAOneway(VectorIF, GRID_ATTR(g), IF_FORWARD, m * sizeof(DOUBLE),
				  Gather_VectorComp, Scatter_VectorComp);

	return (NUM_OK);
}

///////////////////////////////////////////////////////////////////////////////
// l_vector_collectAll
///////////////////////////////////////////////////////////////////////////////

static int Gather_VectorCompCollect (DDD_OBJ obj, void *data)
{
	VECTOR *pv = (VECTOR *)obj;
	INT vc,i,type;
	const SHORT *Comp;	
	
	if (VD_IS_SCALAR(ConsVector)) {
		if (VD_SCALTYPEMASK(ConsVector) & VDATATYPE(pv)) {
		    vc = VD_SCALCMP(ConsVector);
			*((DOUBLE *)data) = VVALUE(pv,vc);
			VVALUE(pv,vc) = 0.0;
		}
		return (NUM_OK);
	  }

	type = VTYPE(pv);
	Comp = VD_CMPPTR_OF_TYPE(ConsVector,type);
	for (i=0; i<VD_NCMPS_IN_TYPE(ConsVector,type); i++) {
		((DOUBLE *)data)[i] = VVALUE(pv,Comp[i]);
		VVALUE(pv,Comp[i]) = 0.0;
	}

	return (NUM_OK);
}

INT l_vector_collectAll (GRID *g, const VECDATA_DESC *x)
{
    INT tp,m; 

    ConsVector = (VECDATA_DESC *)x;

	m = 0;
	for (tp=0; tp<NVECTYPES; tp++)
	  m = MAX(m,VD_NCMPS_IN_TYPE(ConsVector,tp));

	// TODO: both communications within 1 call
	DDD_IFAOneway(BorderVectorIF, GRID_ATTR(g), IF_FORWARD, m * sizeof(DOUBLE),
				  Gather_VectorCompCollect, Scatter_VectorComp);
	DDD_IFAOneway(VectorIF, GRID_ATTR(g), IF_BACKWARD, m * sizeof(DOUBLE),
				  Gather_VectorCompCollect, Scatter_VectorComp);

	return (NUM_OK);
}


pl(FAMGGraph *graph)
{
	FAMGList *list;
	FAMGNode *n;
	list = graph->GetList();
	cout << "LISTE: "<<endl;
	while(list!=NULL)
	{
		cout<<"   data" << list->GetData()<<" first="<<list->GetFirst()->GetId()<<" last="<<list->GetLast()->GetId()<<endl<<flush;
		n = list->GetFirst();
		while( n!=NULL)
		{
			cout<<"       "<<n->GetId();
			if(n->GetPred()!=NULL)
				cout<<" pred= "<<n->GetPred()->GetId();
			if(n->GetSucc()!=NULL)
				cout<<" succ= "<<n->GetSucc()->GetId();
			cout<<endl<<flush;
			n = n->GetSucc();
		}
		list = list->GetSucc();
	}
	cout <<"helplist:";
	n=graph->GetHelpList();
	while(n!=NULL)
	{
		cout<<" "<<n->GetId();
		n = n->GetSucc();
	}
	cout<<endl<<flush;

	return 0;
}
#endif
