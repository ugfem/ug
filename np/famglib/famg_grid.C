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

extern "C"
{
#include "udm.h"
}

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
//		DON'T DO THIS for level < 0 !
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

static int OverlapForLevel0;
static int* CopyPEBuffer = NULL;
#endif

INT l_force_consistence (GRID *g, const VECDATA_DESC *x);
INT l_vector_collectAll (GRID *g, const VECDATA_DESC *x);

// Class FAMGGrid


#ifdef FAMG_SPARSE_BLOCK
void FAMGGrid::Restriction(FAMGVector &fgsolution, FAMGVector &fgdefect, FAMGVector &cgdefect) const
// including smoothing of fine nodes
{
	FAMGVectorEntry fvec;
	const FAMGTransfer &transfer = *GetTransfer();
	FAMGTransferEntry *transfg;
	
#ifdef ModelP
	if (l_vector_collect(mygrid,((FAMGugVector&)fgdefect).GetUgVecDesc())!=NUM_OK) 
		assert(0);
#endif
	
	
	// jacobi smoothing for the fine nodes
    fgsolution.JacobiSmoothFG( *GetMatrix(), fgdefect );
    

#ifdef ModelP
	if (l_vector_consistent(mygrid,((FAMGugVector&)fgsolution).GetUgVecDesc())!=NUM_OK) 
		assert(0);
#endif
	
	
	// correct defect
	Defect(fgdefect, fgdefect, fgsolution);
#ifdef ModelP
	if (l_vector_collect(mygrid,((FAMGugVector&)fgdefect).GetUgVecDesc())!=NUM_OK) 
		assert(0);
#endif
	
	
	const FAMGSparseVector *svfg  = fgdefect.GetSparseVectorPtr();
	const FAMGSparseVector *svcg  = cgdefect.GetSparseVectorPtr();
	const FAMGSparseVector *svr  = transfer.Get_sr();

	cgdefect = 0.0;

	FAMGVectorIter fiter(GetGridVector());
	while( fiter(fvec) )
#ifdef ModelP
		if(IS_FAMG_MASTER(((FAMGugVectorEntryRef*)(fvec.GetPointer()))->myvector()))
#endif	
        {
	    	for(transfg = transfer.GetFirstEntry(fvec); transfg != NULL; transfg = transfg->GetNext())
            {
				// cgdefect[transfg->GetCol()] += transfg->GetRestriction() * fgdefect[fvec];
                SparseBlockMVAddProduct(svcg,svr,svfg,
                                        cgdefect.GetValuePtr(transfg->GetCol()),
                                        transfg->GetRestrictionPtr(),
                                        cgdefect.GetValuePtr(fvec),1.0);
            }
        }
		
#ifdef ModelP
	if (l_vector_collectAll(DOWNGRID(mygrid),((FAMGugVector&)cgdefect).GetUgVecDesc())!=NUM_OK) 
		assert(0);
#endif
	
	
    return;
}
    
void FAMGGrid::Prolongation(const FAMGGrid *cg, const FAMGVector &cgsol, FAMGVector &fgsol, FAMGVector &fgdefect, FAMGVector *c)
// adds the prolongued solution-update to the fine grid solution
// including smoothing of fine nodes
// c is set for FAMG transfer; not used for FAMG solver
{
	FAMGVectorEntry fvec;
	const FAMGTransfer &transfer = *GetTransfer();
	FAMGTransferEntry *transfg;
	
	
#ifdef ModelP
	// distribute master values to V(H)Ghosts
	if (l_ghostvector_consistent(DOWNGRID(mygrid),((FAMGugVector&)cgsol).GetUgVecDesc())!= NUM_OK)
		assert(0);
#endif        
	
	
	const FAMGSparseVector *svfg  = fgsol.GetSparseVectorPtr();
	const FAMGSparseVector *svcg  = cgsol.GetSparseVectorPtr();
	const FAMGSparseVector *svp  = transfer.Get_sp();

	FAMGVectorIter fiter(GetGridVector());
	while( fiter(fvec) )
	{
#ifdef ModelP
		if(!IS_FAMG_MASTER(((FAMGugVectorEntryRef*)(fvec.GetPointer()))->myvector()))
			continue;
#endif		
		// fgsol[fvec] = 0.0;
        SparseBlockVSet(svfg,fgsol.GetValuePtr(fvec),0.0);
	    for(transfg = transfer.GetFirstEntry(fvec); transfg != NULL; transfg = transfg->GetNext())
        {
                SparseBlockMVAddProduct(svfg,svp,svcg,
                                        fgsol.GetValuePtr(fvec),
                                        transfg->GetProlongationPtr(),
                                        cgsol.GetValuePtr(transfg->GetCol()),1.0);
                // fgsol[fvec] += transfg->GetProlongation() * cgsol[transfg->GetCol()];
        }
	}
#ifdef ModelP
	if (l_force_consistence(mygrid,((FAMGugVector&)fgsol).GetUgVecDesc())!=NUM_OK) 
		assert(0);
#endif

	
	// prepare defect for jacobi smoothing
	Defect(fgdefect, fgdefect, fgsol);
#ifdef ModelP
	if (l_vector_collect(mygrid,((FAMGugVector&)fgdefect).GetUgVecDesc())!=NUM_OK) 
		assert(0);
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
        Defect(fgdefect, fgdefect, fgsol);
    }
    else
    {
        // famg as transfer
        (*c) += fgsol;
        fgsol = 0.0;

		// jacobi smoothing for the fine nodes
        fgsol.JacobiSmoothFG( *GetMatrix(), fgdefect );
#ifdef ModelP
		if (l_vector_consistent(mygrid,((FAMGugVector&)fgsol).GetUgVecDesc())!=NUM_OK) 
			assert(0);
#endif


		// defect is computed in ug (Lmgc)
    }

	return;
}
#else

void FAMGGrid::Restriction(FAMGVector &fgsolution, FAMGVector &fgdefect, FAMGVector &cgdefect) const
// including smoothing of fine nodes
{
	FAMGVectorEntry fvec;
	const FAMGTransfer &transfer = *GetTransfer();
	FAMGTransferEntry *transfg;
	
#ifdef ModelP
	if (l_vector_collect(mygrid,((FAMGugVector&)fgdefect).GetUgVecDesc())!=NUM_OK) 
		assert(0);
#endif
	
#ifdef PROTOCOLNUMERIC
    cout<<"FAMGGrid::Restriction: defect before JacobiSmoothFG "<<fgdefect*fgdefect<<endl;
    cout<<"FAMGGrid::Restriction: sol before JacobiSmoothFG "<<fgsolution*fgsolution<<endl;
    cout<<"FAMGGrid::Restriction: defect*sol before JacobiSmoothFG "<<fgdefect*fgsolution<<endl;
#endif
	
	// jacobi smoothing for the fine nodes
    fgsolution.JacobiSmoothFG( *GetConsMatrix(), fgdefect );
#ifdef ModelP
	if (l_vector_consistent(mygrid,((FAMGugVector&)fgsolution).GetUgVecDesc())!=NUM_OK) 
		assert(0);
#endif
	
#ifdef PROTOCOLNUMERIC
    cout<<"FAMGGrid::Restriction: defect after JacobiSmoothFG "<<fgdefect*fgdefect<<endl;
    cout<<"FAMGGrid::Restriction: sol after JacobiSmoothFG "<<fgsolution*fgsolution<<endl;
    cout<<"FAMGGrid::Restriction: defect*sol after JacobiSmoothFG "<<fgdefect*fgsolution<<endl;
#endif
	
	// correct defect
	Defect(fgdefect, fgdefect, fgsolution);
#ifdef ModelP
	if (l_vector_collect(mygrid,((FAMGugVector&)fgdefect).GetUgVecDesc())!=NUM_OK) 
		assert(0);
#endif
#ifdef PROTOCOLNUMERIC
    cout<<"FAMGGrid::Restriction: defect after defect "<<fgdefect*fgdefect<<endl;
    cout<<"FAMGGrid::Restriction: sol after defect "<<fgsolution*fgsolution<<endl;
    cout<<"FAMGGrid::Restriction: defect*sol after defect "<<fgdefect*fgsolution<<endl;
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
	
#ifdef PROTOCOLNUMERIC
    cout<<"FAMGGrid::Restriction: defect after gridtransfer "<<cgdefect*cgdefect<<endl;
#endif
	
    return;
}
    
void FAMGGrid::Prolongation(const FAMGGrid *cg, const FAMGVector &cgsol, FAMGVector &fgsol, FAMGVector &fgdefect, FAMGVector *c)
// adds the prolongued solution-update to the fine grid solution
// including smoothing of fine nodes
// c is set for FAMG transfer; not used for FAMG solver
{
	FAMGVectorEntry fvec;
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
	
	if( 0 )	// Test whether a constant vector is prolongated onto a constant vector
	{
		FAMGVectorIter fiterTest(GetGridVector());
		while( fiterTest(fvec) )
		{
#ifdef ModelP
			if(!IS_FAMG_MASTER(((FAMGugVectorEntryRef*)(fvec.GetPointer()))->myvector()))
				continue;
#endif		
			sum = 0.0;
		    for(transfg = transfer.GetFirstEntry(fvec); transfg != NULL; transfg = transfg->GetNext())
        		sum += transfg->GetProlongation() * cgdefect[transfg->GetCol()];
			fgsol[fvec] = sum;
		}
		prv(GLEVEL(GetugGrid()), VD_SCALCMP(fgsol.GetUgVecDesc()) );	// result should be constant 1
		prv(GLEVEL(cg->GetugGrid()), VD_SCALCMP(cgdefect.GetUgVecDesc()) );	// should be constant 1
	}
	
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
	Defect(fgdefect, fgdefect, fgsol);
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
        Defect(fgdefect, fgdefect, fgsol);
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
#endif

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
    if(vec[FAMGR] == NULL) RETURN(1);
    vec[FAMGV] = solution.create_new();
    if(vec[FAMGV] == NULL) RETURN(1);
    vec[FAMGP] = solution.create_new();
    if(vec[FAMGP] == NULL) RETURN(1);
    vec[FAMGT] = solution.create_new();
    if(vec[FAMGT] == NULL) RETURN(1);

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

#ifdef FAMG_SPARSE_BLOCK
void FAMGGrid::SmoothTV()
{

    FAMGGridVector gv = GetGridVector();
    FAMGVector &tvA = *vector[FAMGTVA];
    FAMGVector &tvB = *vector[FAMGTVB];
    FAMGVector &help = *vector[FAMGUNKNOWN];
	FAMGMatrixAlg &M = *matrix;
    FAMGVectorIter viter(gv);
	FAMGMatrixEntry col;
	FAMGVectorEntry row;
	double *tvAptr, *tvBptr, *helpptr, *mptr;
    int k;

	const FAMGSparseVector *svtvA  = tvA.GetSparseVectorPtr();
	const FAMGSparseVector *svtvB  = tvB.GetSparseVectorPtr();
	const FAMGSparseVector *svh  = help.GetSparseVectorPtr();
	const FAMGSparseBlock *sb  = M.GetSparseBlockPtr();
    const FAMGSparseBlock *sbd  = M.GetDiagSparseBlockPtr();

    const int stv = FAMGGetParameter()->Getstv();

    short nr = sbd->Get_nr();
    FAMGSparseVector svsum(nr);
    double *sum = new double[svsum.Get_maxcomp()+1];
    double *update = new double[svsum.Get_maxcomp()+1];
    double *decomp = new double[nr*nr]; 
    short *pivotmap = new short[nr]; 

    for(k = 0; k < stv; k++)
    {

        while(viter(row))
        {
		
            helpptr = help.GetValuePtr(row);
            tvAptr = tvA.GetValuePtr(row);
        
            // diagonal 
            FAMGMatrixIter miter(M,row);
            miter(col);
            mptr = M.GetValuePtr(col);
            SparseBlockVSet(&svsum,sum,0.0);
            SparseBlockMVAddProduct(&svsum,sbd,svtvA,sum,mptr,tvAptr,1.0);
            SparseBlockMCopyDense(decomp,sbd,mptr);
            while(miter(col))
            {
                tvAptr = tvA.GetValuePtr(col.dest());
                mptr = M.GetValuePtr(col);
                SparseBlockMVAddProduct(&svsum,sb,svtvA,sum,mptr,tvAptr,1.0);
            }
            if(LR_Decomp(nr,decomp,pivotmap)) assert(0);
            if(LR_Solve(nr,decomp,pivotmap,update,sum)) assert(0);
            SparseBlockVAdd(svh,&svsum,&svsum,helpptr,update,update,0.0);
        }
        viter.reset();
        while(viter(row))
        {
            helpptr = help.GetValuePtr(row);
            tvAptr = tvA.GetValuePtr(row);
            SparseBlockVAdd(svtvA,svtvA,svh,tvAptr,tvAptr,helpptr,-1.0);
        }
        viter.reset();
           
	}

    // make sure absolute value is large enough for every tv component
    short n = svtvA->Get_n();
    while(viter(row))
    {
        tvAptr = tvA.GetValuePtr(row);
        tvBptr = tvB.GetValuePtr(row);
        for(short i = 0; i < n; i++)
        {
            if(Abs(tvAptr[svtvA->Get_comp(i)]) < 1e-10)
            {
                tvAptr[svtvA->Get_comp(i)] = 1e-10;
            }
            if(Abs(tvBptr[svtvB->Get_comp(i)]) < 1e-10)
            {
                tvBptr[svtvB->Get_comp(i)] = 1e-10;
            }
        }
    }
    viter.reset();



    if(stv > 0)
    {
        while(viter(row))
        {
            tvAptr = tvA.GetValuePtr(row);
            tvBptr = tvB.GetValuePtr(row);
            SparseBlockVAdd(svtvB,svtvA,svtvA,tvBptr,tvAptr,tvAptr,0.0);
            // test
            // SparseBlockVAdd(&svsum,svtvA,svtvA,tvAptr,tvAptr,tvAptr,0.0);
        }
    }

	return;

}
#else
void FAMGGrid::SmoothTV()
{
    double sum, normA, mii, scalingA, normB, scalingB;
    FAMGVector &tvA = *vector[FAMGTVA];
    FAMGVector &tvB = *vector[FAMGTVB];
    FAMGVector &help = *vector[FAMGUNKNOWN];
	FAMGMatrixAlg &mat = *Consmatrix;
	FAMGMatrixEntry me, de;
	FAMGVectorEntry ve;
    int k;

    const int stv = FAMGGetParameter()->Getstv();

    for(k = 0; k < stv; k++)
    {
		// TODO: noch ueberlegen ob matrix Consmatrix richtig ist
		FAMGVectorIter tviter(tvA);
		while( tviter(ve) )
		{
			FAMGMatrixIter miter(*matrix,ve);
            sum = 0.0;
			
			miter(me);	// diagonal element
			mii = mat[me];
			while( miter(me) )
            {
                FAMGMatrixIter diag(*matrix,me.dest());
                diag(de);
				// sum += mat[me]*tvA[me.dest()]/mat[de];
                sum += mat[me]*tvA[me.dest()];
            }
			//help[ve] = 0.15*tvA[ve] - 0.85*sum;
			help[ve] = 0.15*tvA[ve] - 0.85*sum/mii;
		}	
        tviter.reset();
		while( tviter(ve) ) tvA[ve] = help[ve];



		FAMGVectorIter tvTiter(tvB);
		while( tvTiter(ve) )
		{
			FAMGMatrixIter miter(*matrix,ve);
            sum = 0.0;
			
			miter(me);	// diagonal element
			mii = mat[me];
			while( miter(me) )
            {
                FAMGMatrixIter diag(*matrix,me.dest());
                diag(de);
				 sum += mat[me]*tvB[me.dest()];
				//sum += mat[me]*tvB[me.dest()]/mat[de];
            }
			 help[ve] = 0.15*tvB[ve] - 0.85*sum/mii;
			// help[ve] = 0.15*tvB[ve] - 0.85*sum;
		}	
        tvTiter.reset();
		while( tvTiter(ve) ) tvB[ve] = help[ve];

    }
     
    normA = tvA.norm();
    normB = tvB.norm();
	k = mat.GetN();
#ifdef ModelP
	k = UG_GlobalSumINT( k );
#endif
    scalingA = sqrt((double)k) / normA;
    scalingB = sqrt((double)k) / normB;
	
	FAMGVectorIter tviter(tvA);
	while( tviter(ve) )
	{
		tvA[ve] *= scalingA;
	}
	FAMGVectorIter tvTiter(tvB);
	while( tvTiter(ve) )
	{
		tvB[ve] *= scalingB;
	}

	return;
}
#endif    

#ifdef FAMG_ILU
int FAMGGrid::ILUTDecomp(int cgilut)
{
    FAMGMarkHeap(FAMG_FROM_BOTTOM);
    if(graph == NULL)
    {
        graph = (FAMGGraph *) FAMGGetMem(sizeof(FAMGGraph),FAMG_FROM_BOTTOM);
        if(graph == NULL) RETURN(1);
    }
    if(graph->Init(this)) RETURN(1);
    if(graph->OrderILUT(matrix)) RETURN(1);
    Order(graph->GetMap());
    graph = NULL;
    FAMGReleaseHeap(FAMG_FROM_BOTTOM);

    decomp = (FAMGDecomp *) FAMGGetMem(sizeof(FAMGDecomp),FAMG_FROM_TOP);
    if(decomp == NULL) RETURN(1);

    if(decomp->Init(matrix)) RETURN(1);
    if(decomp->Construct(cgilut)) RETURN(1);

    
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
	ostr << "avg. stencil: " << (double)nl/(double)nn;
#ifdef	USE_UG_DS
	ostr << " level " << GLEVEL(GetugGrid());
#endif
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

    n = system.GetN();
    nf = 0;
	mygridvector = system.GetGridVector();
	
#ifdef USE_UG_DS
	SetugGrid(system.GetFineGrid());
#else
    SetFather(NULL);
#endif
	
    matrix = system.GetMatrix();
    Consmatrix = system.GetConsMatrix();
    tmpmatrix = Consmatrix;
#ifdef FAMG_SPARSE_BLOCK
    diagmatrix = system.GetDiagMatrix();
#endif
	
#ifdef FAMG_ILU
    decomp = NULL;
#endif
	
    transfer = NULL;
    for(i = 0; i < FAMGMAXVECTORS; i++)
		vector[i] = system.GetVector(i);
	
#ifdef FAMG_REORDERmap	
    map = (int *) FAMGGetMem(n*sizeof(int),FAMG_FROM_TOP);
    if(map == NULL) RETURN(1);
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
		RETURN(1);
#ifdef FAMG_SPARSE_BLOCK
	diagmatrix = new FAMGugMatrix(new_grid, *(FAMGugMatrix*)grid_pattern.GetDiagMatrix());
    if(diagmatrix == NULL)
		RETURN(1);
#endif
	if(grid_pattern.GetMatrix()!=grid_pattern.GetConsMatrix())
	{
		Consmatrix = new FAMGugMatrix(new_grid, *(FAMGugMatrix*)grid_pattern.GetConsMatrix());
	    if(Consmatrix == NULL)
			RETURN(1);
	}
	else
		Consmatrix = matrix;
#else	
    matrix = (FAMGMatrix *) FAMGGetMem(sizeof(FAMGMatrix),FAMG_FROM_TOP);
    if(matrix == NULL)
		RETURN(1);
    if(matrix->Init(n))
		RETURN(1);
	mygridvector = NULL; // aendern
#endif
	
    tmpmatrix = Consmatrix;
    transfer = NULL;
	
#ifdef FAMG_ILU
    decomp = NULL;
#endif
	
#ifndef USE_UG_DS	
    father = (int *) FAMGGetMem(n*sizeof(int),FAMG_FROM_TOP);
    if(father == NULL) RETURN(1);
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
			RETURN(1);
		}
        SetVector(i,new_vector);
    }

#ifdef UG_DRAW
    vertex = (void **) FAMGGetMem(n*sizeof(void *),FAMG_FROM_TOP);
    if(vertex == NULL) RETURN(1);
#endif
	
    GetSmoother();

    return 0;
}


void FAMGGrid::Deconstruct(int level)
{
	int i;

	if( level == 0 )
		return;	// objects for level 0 have not been allocated by FAMGGrid, thus are not allowed to be freed here
	
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
    int j;

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
#ifdef FAMG_SPARSE_BLOCK
    j = 0;
    const FAMGSparseVector *svAcg = tvAcg.GetSparseVectorPtr(); 
    const FAMGSparseVector *svAfg = tvAfg.GetSparseVectorPtr(); 
    const FAMGSparseVector *svBcg = tvBcg.GetSparseVectorPtr(); 
    const FAMGSparseVector *svBfg = tvBfg.GetSparseVectorPtr(); 
    double *tvAcg_val, *tvAfg_val, *tvBcg_val, *tvBfg_val; 

	while( viter(fg_ve) )
	{
		if (fg_gridvec.IsCG(fg_ve))
		{
			cg_ve = GetTransfer()->GetFirstEntry(fg_ve)->GetCol();
            tvAcg_val = tvAcg.GetValuePtr(cg_ve); 
            tvAfg_val = tvAfg.GetValuePtr(fg_ve); 
            tvBcg_val = tvBcg.GetValuePtr(cg_ve); 
            tvBfg_val = tvBfg.GetValuePtr(fg_ve); 
            SparseBlockVCopy(svAcg, svAfg, tvAcg_val, tvAfg_val, 1.0); 
            SparseBlockVCopy(svBcg, svBfg, tvBcg_val, tvBfg_val, 1.0); 
			j++;
		}
	}
#else
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
#endif

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
        for(int i = 0; i < fn; i++)
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
		RETURN(1);
	
    return 0;
}

#ifdef FAMG_SINGLESTEP
#define MIN_data 0
#define MIN_pe 1
#define xCOORD 2
#define yCOORD 3
#define BufferSizeGlobalBestChoice 4
#define eps 0.00001

void FAMG_GlobalBestChoice( DOUBLE choice[BufferSizeGlobalBestChoice] )
// OLD: Global optimum is: min data (and within them: min pe --- to be unique)
// Global optimum is: min data (and within them: lexicogrphic smallest); MIN_pe currently not used 
// or equivalently: search global node with lex. min. of (data,x,y)
{
	int l;
	DOUBLE n[BufferSizeGlobalBestChoice];

	for (l=degree-1; l>=0; l--)
	{
		GetConcentrate(l,n,BufferSizeGlobalBestChoice*sizeof(DOUBLE));
		// OLD: if( (n[MIN_data] < choice[MIN_data]) || ( (n[MIN_data] == choice[MIN_data]) && (n[MIN_pe] < choice[MIN_pe]) ) )
		if( (n[MIN_data] < choice[MIN_data]-eps) || ( (ABSDIFF(n[MIN_data],choice[MIN_data])<eps) && ((n[xCOORD] < choice[xCOORD]-eps) || ( (ABSDIFF(n[xCOORD],choice[xCOORD])<eps) && (n[yCOORD] < choice[yCOORD]-eps) )) ) )
		{
			choice[MIN_data] = n[MIN_data];
			choice[MIN_pe] = n[MIN_pe];
			choice[xCOORD] = n[xCOORD];
			choice[yCOORD] = n[yCOORD];
		}
	}
	Concentrate(choice,BufferSizeGlobalBestChoice*sizeof(DOUBLE));
	Broadcast(choice,BufferSizeGlobalBestChoice*sizeof(DOUBLE));
	return;
}
#endif // FAMG_SINGLESTEP


int FAMGGrid::ConstructTransfer()
// Do not leave this function premature in ModelP! This can cause deadlocks.
{
    int i, conloops;
	double TotalTime;

#ifdef ModelP
	double GraphColorTime, BorderTime;
	FAMGNode *nodei;
	VECTOR *vec;
	
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

    FAMGMarkHeap(FAMG_FROM_BOTTOM);

    graph = (FAMGGraph *) FAMGGetMem(sizeof(FAMGGraph),FAMG_FROM_BOTTOM);
    if (graph == NULL) {FAMGReleaseHeap(FAMG_FROM_BOTTOM);  RETURN(1);}

    // test
    // FGSSmoothTV();
    
	if( FAMGGetParameter()->Gettype() == 6 )
		GetTmpMatrix()->MarkStrongLinks(*this);	// needed only for type 6 in the moment
	
    if (graph->Init(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}
    if (graph->Construct(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}

    transfer = (FAMGTransfer *) FAMGGetMem(sizeof(FAMGTransfer),FAMG_FROM_TOP);
    if(transfer == NULL) assert(0);
    if (transfer->Init(this)) assert(0);

#ifdef FAMG_SINGLESTEP
	#define DUMMY_DATA MAX_I
	#define DUMMY_PE (procs+1)
	#define DUMMY_COORD MAX_D

	INT finished, step=0;
	DOUBLE choice[BufferSizeGlobalBestChoice];

	#if !defined __TWODIM__
	#error only for 2 dim
	#endif

	GraphColorTime = -999.999;	// init
	BorderTime = -999.999;		// init

	if( level == 0 )
	{
		if (graph->EliminateDirichletNodes(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}
		CommunicateNodeStatus();
	}
	if( graph->InsertHelplist() ) {FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}

	finished = 0;
	//prv(level,0);//prm(level,1);

	while( !finished )
	{
		//cout<<me<<": #"<<step<<endl; printlist(graph);

		// taken from FAMGGrid::EliminateNodes
		if( graph->GetList() == NULL ) 
			nodei = NULL;
		else 
			nodei = graph->GetList()->GetFirst();

		if( nodei==NULL )
		{
			choice[MIN_data] = DUMMY_DATA;
			choice[MIN_pe] = DUMMY_PE;
			choice[xCOORD] = DUMMY_COORD;
			choice[yCOORD] = DUMMY_COORD;
		}
		else
		{
			// search node from the first sublist, which is lexicographic the smallest
			DOUBLE xmin=DUMMY_COORD, ymin=DUMMY_COORD;
			FAMGNode *nod = nodei;

			while( nod != NULL )
			{
				DOUBLE pos[DIM];
				VECTOR *vec = ((FAMGugVectorEntryRef*)(nod->GetVec().GetPointer()))->myvector();
				VectorPosition(vec,pos);

				if( (pos[0] < xmin-eps) || ( (ABSDIFF(pos[0],xmin)<eps) && (pos[1] < ymin-eps) ) )
				{
					nodei = nod;
					xmin = pos[0];
					ymin = pos[1];
				}

				nod = nod->GetSucc();
				//cout<<me<<": "<<step<<"! sc "<<pos[0]<<" "<<pos[1]<<" "<<xmin<<" "<<ymin<<VINDEX(vec)<<endl;
			}


			if(nodei->GetPaList() == NULL)
			{
				// done
				choice[MIN_data] = DUMMY_DATA;
				choice[MIN_pe] = DUMMY_PE;
				choice[xCOORD] = DUMMY_COORD;
				choice[yCOORD] = DUMMY_COORD;
				//cout << me<<": "<<step<<"! choice palist=0 "<<"node "<<nodei->GetId()<<endl;
			}
			else if (nodei->CheckPaList(graph))
			{
				//cout << me<<": "<<step<<"! choice checkpalist "<<"node "<<nodei->GetId()<<endl;
				// CheckPaList is necessary for non structure
				//   symmetric matrices. It is not nice and
				//   should be avoided. 
				graph->Remove(nodei);
				nodei->ComputeTotalWeight();
				if(graph->Insert(nodei)) RETURN(1);
				choice[MIN_data] = DUMMY_DATA;
				choice[MIN_pe] = DUMMY_PE;
				choice[xCOORD] = DUMMY_COORD;
				choice[yCOORD] = DUMMY_COORD;
			}           
			else
			{
				choice[MIN_data] = nodei->GetData();
				choice[MIN_pe] = me;
				choice[xCOORD] = xmin;
				choice[yCOORD] = ymin;
				//cout<<me<<": "<<step<<"! sbuf "<<xmin<<" "<<ymin<< " "<<choice[xCOORD]<<" "<<choice[yCOORD]<<endl;
			}
		}

		// Communicate global best choice
		//cout << me<<": "<<step<<"! choice vor "<< choice[MIN_data] <<" node "<<nodei->GetId()<<endl;
		FAMG_GlobalBestChoice ( choice );
		//cout << me<<": "<<step<<"! choice nach "<< choice[MIN_data]<<" pe "<<choice[MIN_pe] <<endl;

		finished = (ABSDIFF(choice[MIN_pe],DUMMY_PE)<eps); // no PE left which has a node to be eliminated

		if( !finished )
		{
			if( ABSDIFF(me,choice[MIN_pe])<eps )
			{	// Do the elimination
				//cout << me<<": "<<step<<"! el "<<nodei->GetId()<<" w= "<<choice[MIN_data]<<" x= "<<choice[xCOORD]<<" y= "<<choice[yCOORD]<<endl;
				graph->Remove(nodei);
				if(nodei->Eliminate(this)) RETURN(1);
				if(nodei->UpdateNeighborsFG(this)) RETURN(1); 
				if(graph->InsertHelplist()) RETURN(1);
			}


			CommunicateNodeStatus();
			if(graph->InsertHelplist()) RETURN(1);
		}

		step++;
	}

	if(graph->InsertHelplist()) RETURN(1);
	if (graph->RemainingNodes()) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}
	CommunicateNodeStatus();

#else // FAMG_SINGLESTEP

#ifdef ModelP
	if( level == 0 )
	{
		if (graph->EliminateDirichletNodes(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}
		CommunicateNodeStatus();
	}
	
#ifdef FAMG_INNER_FIRST
	if( graph->InsertHelplist() )
		RETURN(1);

	if (graph->EliminateNodes(this))
		RETURN(1);

	// let undecided nodes in the list; perhaps they can be eliminated
	// in the border step; if not they become coarse there
	//if (graph->RemainingNodes()) RETURN(1);

	CommunicateNodeStatus();
	
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
					RETURN(1);
				}
		}
	}
#endif // FAMG_INNER_FIRST

	// in parallel now only the nodes in the border of the core partition are in the list
	
	int NrNbPe;
	FAMGColor MyColor, BorderCycles;
	
	GraphColorTime = CURRENT_TIME;

	if( (NrNbPe=ConstructColoringGraph(GRID_ATTR(mygrid),FAMGGetParameter()->GetColoringMethod())) < 0) 
	{
		cout << "FAMGGrid::ConstructTransfer(): ERROR in constructing the coloring graph"<<endl<<fflush;
		FAMGReleaseHeap(FAMG_FROM_BOTTOM);
		RETURN(1);
	}

	if( ConstructColoring( FAMGGetParameter()->GetColoringMethod(), MyColor, BorderCycles ) )
	{
		cout << "FAMGGrid::ConstructTransfer(): ERROR in Coloring the graph"<<endl<<fflush;
		FAMGReleaseHeap(FAMG_FROM_BOTTOM);
		RETURN(1);
	}
		
	GraphColorTime = CURRENT_TIME - GraphColorTime;
	cout <<me<<": level = "<<level<<" my color = "<<MyColor<<" max. color = "<<BorderCycles<<" ColoringMethod = "<<FAMGGetParameter()->GetColoringMethod()<<" NrNbPe = "<<NrNbPe<<endl;

	BorderTime = CURRENT_TIME_LONG;
	for ( int color = 0; color <= BorderCycles; color++)
	{
		if( MyColor == color )
		{
		    if (graph->InsertHelplist()) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}
		    if (graph->EliminateNodes(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}
    		if (graph->RemainingNodes()) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}
		}
		
//prim(GLEVEL(GetugGrid()));//?????????????????????????????????????????????
		CommunicateNodeStatus();
	}	
	BorderTime = CURRENT_TIME_LONG - BorderTime;

#ifndef FAMG_INNER_FIRST
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
					RETURN(1);
				}
		}
    }
#endif //FAMG_INNER_FIRST

#endif // ModelP

#ifdef SIMULATE_HALFENING	// TODO: remove it
	FAMGNode *nodei;
	VECTOR *vec;
	
	if (graph->InsertHelplist()) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}
	if (graph->EliminateNodes(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}
	if (graph->RemainingNodes()) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}
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
					RETURN(1);
				}
		}
    }
#endif
	
#if !(defined ModelP && defined FAMG_INNER_FIRST) // i.e !ModelP || (ModelP && !FAMG_INNER_FIRST)
	if( graph->InsertHelplist() ) {FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}
    if (graph->EliminateNodes(this)) {FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}
#endif 

    for(i = 0; i < conloops; i++)
    {
        // FAMGMarkHeap(FAMG_FROM_BOTTOM);
#ifdef USE_UG_DS
#ifndef FAMG_SPARSE_BLOCK
        // assert(0);	// i don't want to have a second matrix. [Maybe we do not need one ? :-) chris] 

        if (graph->Construct2(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}
        if( graph->InsertHelplist() ) {FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);} // can not remember why
        if (graph->EliminateNodes(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}

#endif
#else
		
        tmpmatrix = (FAMGMatrixAlg *) FAMGGetMem(sizeof(FAMGMatrixAlg),FAMG_FROM_BOTTOM);
        if(tmpmatrix == NULL)
			RETURN(1);
        if(tmpmatrix->Init2(n))
			RETURN(1);
        if(tmpmatrix->TmpMatrix(matrix,transfer,graph))
			RETURN(1);
#endif
#ifdef FAMG_SPARSE_BLOCK
 
        if (graph->Construct2(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}
        if( graph->InsertHelplist() ) {FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);} // can not remember why
        if (graph->EliminateNodes(this)) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}
#endif
        // FAMGReleaseHeap(FAMG_FROM_BOTTOM);
    }

#if !(defined ModelP && defined FAMG_INNER_FIRST) // i.e !ModelP || (ModelP && !FAMG_INNER_FIRST)
    if (graph->RemainingNodes()) { FAMGReleaseHeap(FAMG_FROM_BOTTOM); RETURN(1);}

#ifdef ModelP
	// update the ghost and border nodes
	//prim(GLEVEL(GetugGrid()));//?????????????????????????????????????????????
	CommunicateNodeStatus();

#endif
#endif // #if !(defined ModelP && defined FAMG_INNER_FIRST)

#endif //FAMG_SINGLESTEP

    nf = graph->GetNF();
    // TODO: not neceassary any more: matrix->MarkUnknowns(graph);
 
	TotalTime = CURRENT_TIME_LONG - TotalTime;
	cout
#ifdef ModelP
		<< me << ": "
		<< "level "<< level 
#endif
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
    if (helpvect == NULL) RETURN(1);
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
    if (helpvect == NULL) RETURN(1);
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
        if (helpfather == NULL) RETURN(1);
        for(i = 0; i < n; i++) helpfather[i] = father[i];
        for(i = 0; i < n; i++) father[mapping[i]] = helpfather[i];
        FAMGReleaseHeap(FAMG_FROM_TOP);
    }

    /* order vertices */
    if(vertex != NULL)
    {
        FAMGMarkHeap(FAMG_FROM_TOP);
        helpvertex = (void **) FAMGGetMem(n*sizeof(void *), FAMG_FROM_TOP);
        if (helpvertex == NULL) RETURN(1);
        for(i = 0; i < n; i++) helpvertex[i] = vertex[i];
        for(i = 0; i < n; i++) vertex[mapping[i]] = helpvertex[i];
        FAMGReleaseHeap(FAMG_FROM_TOP);
    }
        
    /* order test vectors */
    if(OrderVector(FAMGTVA,mapping)) RETURN(1);
    if(OrderVector(FAMGTVB,mapping)) RETURN(1);
   
    /* order matrix */
    if (matrix->Order(mapping)) RETURN(1);

    /* order transfer */
    if(transfer != NULL)
    {
        if (transfer->Order(mapping)) RETURN(1);
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
        if (helpfather == NULL) RETURN(1);
        for(i = 0; i < n; i++) helpfather[i] = father[i];
        for(i = 0; i < n; i++) father[i] = helpfather[map[i]];
        FAMGReleaseHeap(FAMG_FROM_TOP);
    }

    /* reorder vertices */
    if(vertex != NULL)
    {
        FAMGMarkHeap(FAMG_FROM_TOP);
        helpvertex = (void **) FAMGGetMem(n*sizeof(void *), FAMG_FROM_TOP);
        if (helpvertex == NULL) RETURN(1);
        for(i = 0; i < n; i++) helpvertex[i] = vertex[i];
        for(i = 0; i < n; i++) vertex[i] = helpvertex[map[i]];
        FAMGReleaseHeap(FAMG_FROM_TOP);
    }
        
    /* reorder test vectors */
    if(ReorderVector(FAMGTVA,map)) RETURN(1);
    if(ReorderVector(FAMGTVB,map)) RETURN(1);
    
    /* reorder matrix */
    if (matrix->Reorder(map)) RETURN(1);

    /* reorder transfer */ // debug
    if(transfer != NULL)
    {
        if (transfer->Reorder(map)) RETURN(1);
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

void CountOverlap (GRID *g)	// only for testing; TODO: remove 
// destroys VINDEX
{
	VECTOR *vec, *nb;
	MATRIX *mat;
	unsigned int level, i;
	int nrborder = NVEC_PRIO(g,PrioBorder);
	int shelllevel[ 10000 ], inside, outside;
	
	for ( level = 0; level < 10000; level++ )
		shelllevel[level] = 0;
	
	for( vec=PFIRSTVECTOR(g); vec!=NULL; vec=SUCCVC(vec))
		VINDEX(vec) = 0;
	
	for( vec=PFIRSTVECTOR(g); vec!=NULL; vec=SUCCVC(vec))
		if ( PRIO(vec)==PrioMaster )
		{
			for( mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat) )
			{
				nb = MDEST(mat);
				if( PRIO(nb)!=PrioMaster )
				{
					assert(PRIO(nb)==PrioBorder);
					if( VINDEX(nb) == 0 )
					{
						VINDEX(nb) = 1;
						nrborder--;
						shelllevel[1]++;
					}
				}
			}
		}
	
	level = 1;
	while( nrborder > 0 )
	{
		for( vec=PFIRSTVECTOR(g); vec!=NULL; vec=SUCCVC(vec))
			if ( PRIO(vec)==PrioBorder && VINDEX(vec) == level)
			{
				for( mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat) )
				{
					nb = MDEST(mat);
					if( PRIO(nb)!=PrioMaster )
					{
						assert(PRIO(nb)==PrioBorder);
						if( VINDEX(nb) == 0 )
						{
							VINDEX(nb) = level+1;
							nrborder--;
							shelllevel[level+1]++;
						}
					}
				}
			}
		level++;
		assert(level<10000);
	}
	assert( nrborder== 0 );
	
	inside = shelllevel[1]+shelllevel[2];
	outside = 0;
	for( i = 3; i <= level; i++ )
		outside += shelllevel[i];
	printf ("%d: GLEVEL %d %d Bordervec on %d shelllevels in %d out %d:", me, GLEVEL(g), NVEC_PRIO(g,PrioBorder), level, inside, outside);
	for( i = 1; i <= level; i++ )
		printf (" [ %d ]= %d", i, shelllevel[i]);
	printf ("\n");
}

static void TransferVector( VECTOR *v, DDD_PROC dest_pe )
{	
	int size = sizeof(VECTOR)-sizeof(DOUBLE)+FMT_S_VEC_TP(MGFORMAT(dddctrl.currMG),VTYPE(v));
	DDD_XferCopyObjX(PARHDR(v), dest_pe, PrioBorder, size);
	
#ifndef FAMG_SEND_NODES
	if( OverlapForLevel0 )
	{
#endif
		/* for debugging purpose send the node too */
		/* TODO: be sure that this isn't done for production runs */
		IFDEBUG(np,0)
		GEOM_OBJECT *obj = VOBJECT(v);
		if( obj != NULL )
			if (VOTYPE(v) == NODEVEC)
			{
				PRINTDEBUG(dddif,2,(PFMT " TransferVector(): node for vec= " VINDEX_FMTX" n=" ID_FMTX " Xfer v=" 
					VID_FMTX "\n",me,VINDEX_PRTX(v),ID_PRTX((NODE*)obj),VID_PRTX(MYVERTEX((NODE*)obj))))
				DDD_XferCopyObj(PARHDR((NODE*)obj), dest_pe, PrioBorder);
			}
			else
				assert(0); /* unrecognized vector type; extend code if necessary */	
		ENDDEBUG
#ifndef FAMG_SEND_NODES
	}
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
	for( i=2; proclist[i]!=-1; i+=2 )
		if( proclist[i+1] == PrioMaster )
		{
			masterPe = proclist[i];
			break;
		}
	assert(masterPe!=(DDD_PROC)me);	// a master copy must exist on an other processor

	PRINTDEBUG(np,1,("%d: SendToMaster %d: myself "VINDEX_FMTX"\n",me,masterPe,VINDEX_PRTX(vec)));
	
	// It seems to be necessary to send also the vec itself to let new 
	// connection be constructed from both (vec and its neighbor)
	for( mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat) )
	{
		w = MDEST(mat);

		PRINTDEBUG(np,1,("%d: SendToMaster %d:     -> "VINDEX_FMTX"\n",me,masterPe,VINDEX_PRTX(w)));
		TransferVector(w, masterPe);
	}
	
	return 0;
}

#ifdef FAMG_SINGLESTEP
static int SendToOverlap1FULL( DDD_OBJ obj)
// every master sends itself and its neighbors to each processor where it is a border and a neighbor has a master copy
// the reason: on that processor the vector is in overlap1 and must extend with its neighbors the overlap2
{
	VECTOR *vec = (VECTOR *)obj, *w;
	MATRIX *mat;
	int dest_pe;
	
	if( !IS_FAMG_MASTER(vec) )
		return 0;		// we want only master vectors here

	PRINTDEBUG(np,1,("%d: SendToOverlap1: "VINDEX_FMTX"\n",me,VINDEX_PRTX(vec)));

	// traverse all pes where vec has a border/ghost copy
	for( dest_pe=0; dest_pe<procs; dest_pe++ )
	{
		if( me == dest_pe )
			continue;	// skip myself

		// now we found a pe that should receive the neighborhood of vec (incl. vec)
		for( mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat) )
		{
			w = MDEST(mat);
			PRINTDEBUG(np,1,("%d: SendToOverlap1 %d:     -> "VINDEX_FMTX"\n",me,dest_pe,VINDEX_PRTX(w)));
			TransferVector(w, dest_pe);
		}
	}

	return 0;
}
#endif // FAMG_SINGLESTEP

static int SendToOverlap1( DDD_OBJ obj)
// every master sends itself and its neighbors to each processor where it is a border and a neighbor has a master copy
// the reason: on that processor the vector is in overlap1 and must extend with its neighbors the overlap2
{
	VECTOR *vec = (VECTOR *)obj, *w;
	MATRIX *mat;
	int *proclist_vec, *proclist_w, i, size, dest_pe, found, *destPE_vec_ptr;
	
	if( !IS_FAMG_MASTER(vec) )
		return 0;		// we want only master vectors here

	PRINTDEBUG(np,1,("%d: SendToOverlap1: "VINDEX_FMTX"\n",me,VINDEX_PRTX(vec)));

	// construct list of all pes where vec has a copy
	destPE_vec_ptr = CopyPEBuffer;	// init
	proclist_vec = DDD_InfoProcList(PARHDR(vec));
	proclist_vec += 2; 	// skip entry for me
	while( (*destPE_vec_ptr++ = proclist_vec[0]) != -1 )
	{
		proclist_vec += 2;
	}
	
	// traverse all pes where vec has a border/ghost copy
	for( destPE_vec_ptr = CopyPEBuffer; *destPE_vec_ptr!=-1; destPE_vec_ptr++ )
	{
		dest_pe = *destPE_vec_ptr;
		
		// search whether vec has a neighbor w with a master copy on dest_pe 
		for( mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat) )
		{
			w = MDEST(mat);
		
			if( IS_FAMG_MASTER(w) )
				continue;	// if w is master here, he can not be master on an other processor
			
			// look if w has a master copy on dest_pe and thus vec's copy is there in overlap1
			found = 0;
			proclist_w = DDD_InfoProcList(PARHDR(w));
			for( ; proclist_w[0]!=-1; proclist_w += 2 )
				if( proclist_w[0] == dest_pe && proclist_w[1] == PrioMaster )
				{
					found = 1;
					break;
				}

			if( !found )
				continue;	// check the next neighbor
			
			// now we found a pe that should receive the neighborhood of vec (incl. vec)
			for( mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat) )
			{
				w = MDEST(mat);
				PRINTDEBUG(np,1,("%d: SendToOverlap1 %d:     -> "VINDEX_FMTX"\n",me,dest_pe,VINDEX_PRTX(w)));
				TransferVector(w, dest_pe);
			}
			
			break; // check next copy of vec
		}
	}

	return 0;
}

void FAMGGrid::ConstructOverlap()
// extend the overlap as far as necessary; at least 2 links deep
// the vectorlist will be renumbered
{
	VECTOR *vec;
	INT i, mc = MD_SCALCMP(((FAMGugMatrix*)GetMatrix())->GetMatDesc());
	FAMGMatrixAlg *matrix_tmp;
	MATRIX *mat;

	if(GLEVEL(mygrid)==0)
	{
ASSERT(!DDD_ConsCheck());
		OverlapForLevel0 = 1;
		
		// do modifications for dirichlet vectors (as coarsegrid solver)
		// communicate vecskip flags and dirichlet values 
		if (a_vector_vecskip(MYMG(mygrid),0,0,((FAMGugVector*)GetVector(FAMGUNKNOWN))->GetUgVecDesc()) != NUM_OK)
			abort();		
	}
	else
	{
		OverlapForLevel0 = 0;
	}
	
	FAMGMarkHeap(FAMG_FROM_BOTTOM);
	CopyPEBuffer = (int *) FAMGGetMem((procs+1)*sizeof(int),FAMG_FROM_BOTTOM);
	if( CopyPEBuffer==NULL )
	{
		cout << "FAMGGrid::ConstructOverlap ERROR not enough memory for CopyPEBuffer" << endl << fflush;
		abort();
	}
	
	DDD_XferBegin();
    #ifdef DDDOBJMGR
    DDD_ObjMgrBegin();
    #endif
		DDD_IFAExecLocal( OuterVectorSymmIF, GRID_ATTR(mygrid), SendToMaster );
    #ifdef DDDOBJMGR
    DDD_ObjMgrEnd();
    #endif
	DDD_XferEnd();
	
#ifdef FAMG_SINGLESTEP
	int equal = 0, step = 0;
	int minvec, maxvec;

	while( !equal )
	{
		DDD_XferBegin();
    	#ifdef DDDOBJMGR
   		DDD_ObjMgrBegin();
    	#endif
			DDD_IFAExecLocal( OuterVectorSymmIF, GRID_ATTR(mygrid), SendToOverlap1FULL );
    	#ifdef DDDOBJMGR
    	DDD_ObjMgrEnd();
    	#endif
		DDD_XferEnd();

		minvec = maxvec = NVEC(mygrid);
		minvec = UG_GlobalMinINT(minvec);
		maxvec = UG_GlobalMaxINT(maxvec);

		equal = (minvec==maxvec);	

		step++;
		if( me == master )
			cout<<" 0: ConstrOverl "<<step<<" min "<<minvec<<" max "<<maxvec<<endl;
	}
#else
	DDD_XferBegin();
    #ifdef DDDOBJMGR
    DDD_ObjMgrBegin();
    #endif
		DDD_IFAExecLocal( OuterVectorSymmIF, GRID_ATTR(mygrid), SendToOverlap1 );
    #ifdef DDDOBJMGR
    DDD_ObjMgrEnd();
    #endif
	DDD_XferEnd();
#endif // FAMG_SINGLESTEP

	// In some cases ghost vectors are left. They disturb the calculations 
	// and must be removed here.
	// I see 2 alternatives:
	//	A) Delete this ghosts
	//	B) Change ghosts to border
	// B) would be a little faster, because it needs the little cheaper Prio-environment 
	// instead of the Xfer-environment, but it produces in the moment an 
	// inconsistent data structure w.r.t. blasm.c.
	if(OverlapForLevel0)
	{
		DDD_XferBegin();
	  	#ifdef DDDOBJMGR
		DDD_ObjMgrBegin();
	   	#endif

		VECTOR *mastervec = FIRSTVECTOR(mygrid);
		for( vec=PFIRSTVECTOR(mygrid); vec!=mastervec; vec=SUCCVC(vec))
			DDD_XferDeleteObj(PARHDR(vec));
    		#ifdef DDDOBJMGR
		DDD_ObjMgrEnd();
		#endif
		DDD_XferEnd();
	}

	#ifdef DIAGMATWORKSAFTERPRIOCHANGE
	// The following piece of code would be a (little) faster alternative to the above part. 
	// But unfortunately even after adding the diagoanl matrix entry (if needed) 
	// blasm.c doesn't work (can not find a consistent diagonal matrix entry).
	// Thus omit it.
	if(OverlapForLevel0)
	{
		VECTOR *tmpvec, *mastervec = FIRSTVECTOR(mygrid);
		MATRIX *mat;

		DDD_PrioBegin();
		for( vec=PFIRSTVECTOR(mygrid); vec!=mastervec; )
		{
			if( VSTART(vec) == NULL || MDEST(VSTART(vec))!=vec )
			{
				if((mat=CreateConnection(mygrid, vec, vec))==NULL)
				{
					cerr << me <<": ERROR: Cannot create diagonal connection for former ghost\n"<<endl<<fflush;
					abort();
				}
				memset(mat+sizeof(MATRIX)-sizeof(DOUBLE), 0, MSIZE(mat));
			}
			// note: DDD_PrioChange changes immediately the vectorlist; since tmpvec is necessary
			tmpvec = vec;
			vec = SUCCVC(vec);
			DDD_PrioChange(PARHDR(tmpvec),PrioBorder);
		}
		DDD_PrioEnd();
	}
	#endif

	FAMGReleaseHeap(FAMG_FROM_BOTTOM);
	
//CountOverlap (mygrid);	// TODO: remove; only for testing

	// renumber vector list
	i = 0;
	for( vec=PFIRSTVECTOR(mygrid); vec!=NULL; vec=SUCCVC(vec))
		VINDEX(vec) = i++;
	
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
		// For saftey: once again make consistent Dirichlet information; could be skipped perhaps? 
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
	
	if( node->IsUndecidedNode() )
	{
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
					// The missing node must be out of core partition and out of
					// overlap 2 (otherwise it would be present).
					// Skipping it is ok since the solution process requires values only
					// in the core partition and in overlap of depth 1. 
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
				RETURN(1);
			}
			Communication_Graph->ClearPaList(node->GetPaList());		
			node->SetPaList(palist);
		
			// process node as fine node
    	    Communication_Graph->Remove(node);
			node->Eliminate(Communication_Grid);
			node->UpdateNeighborsFG(Communication_Grid);
		}
		else
		{
			abort(); // unrecognized message type
		}
	}
#ifdef Debug
	else
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
			abort(); // unrecognized message type
		}
	}
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

void ppal( FAMGNode *node )
{
	FAMGPaList *palist = node->GetPaList();

	cout << " PaList = ";
	while( palist!=NULL )
	{
		cout << "(";
		for( int i=0; i<palist->GetNp(); i++ )
		{
			if( i!=0 )
				cout << ";";
			cout << palist->GetPa(i);
		}
		cout << ")";
		//cout << "A"<<palist->GetApprox()<<"L"<<palist->GetNewLinks()<<"C"<<palist->GetNewCG();
		//cout << "T"<<palist->TotalWeight();
		palist = palist->GetNext();
	}
	//cout << "\n";	
}


void printlist(FAMGGraph *graph)
{
	FAMGList *list;
	FAMGNode *n;
	list = graph->GetList();
	#ifdef ModelP
	cout << me << ": ";
	#endif
	cout << "LISTE: "<<endl;
	while(list!=NULL)
	{
		#ifdef ModelP
		cout << me << ": ";
		#endif
		cout<<"   data " << list->GetData()<<" first="<<list->GetFirst()->GetId()<<" last="<<list->GetLast()->GetId()<<endl<<flush;
		n = list->GetFirst();
		while( n!=NULL)
		{
			#ifdef ModelP
			cout << me << ": ";
			#endif
			cout<<"       "<<n->GetId();
			if(n->GetPred()!=NULL)
				cout<<" pred= "<<n->GetPred()->GetId();
			if(n->GetSucc()!=NULL)
				cout<<" succ= "<<n->GetSucc()->GetId();
			ppal(n);
			cout<<endl<<flush;
			n = n->GetSucc();
		}
		list = list->GetSucc();
	}
	#ifdef ModelP
	cout << me << ": ";
	#endif
	cout <<"helplist:";
	n=graph->GetHelpList();
	while(n!=NULL)
	{
		cout<<" "<<n->GetId();
		ppal(n);
		n = n->GetSucc();
	}
	cout<<endl<<flush;
}
#endif
