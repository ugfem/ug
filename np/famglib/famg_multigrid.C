/****************************************************************************/
/*																			*/
/* File:      famg_multigrid.C												*/
/*																			*/
/* Purpose:   famg   multigrid classes functions							*/
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
#include <math.h>

#include "famg_misc.h"
#include "famg_grid.h"
#include "famg_multigrid.h"
#include "famg_graph.h"
#include "famg_algebra.h"
#include "famg_heap.h"
#include "famg_system.h"
#include "famg_sparse.h"

/* RCS_ID
$Header$
*/

// data structure local to this file
struct leaveinfo
{
	double coarsefrac;
	int cgnodes;
#ifdef ModelP
	int cgminnodespe;
#endif
};

typedef struct leaveinfo FAMGLeaveInfo;


// Class FAMGMultigrid



int FAMGMultiGrid::Init(const FAMGSystem &system)
{
    FAMGGrid *grid0;
    int i;


    n = 0;
    for(i = 0; i < FAMGMAXGRIDS; i++) grid[i] = NULL;

    grid0 = (FAMGGrid *) FAMGGetMem(sizeof(FAMGGrid),FAMG_FROM_TOP);
    if(grid0 == NULL)
		RETURN(1);

    if( grid0->InitLevel0(system))
		RETURN(1);
    grid[n] = grid0;
    n++;

    return 0;
}    

#ifdef ModelP
static void GlobalLeave (FAMGLeaveInfo *x)
// input: x filled with my data
// output: x filled with global values according to
// calculate GlobalAvg(x->coarsefrac)
//           GlobalSum(x->cgnodes)
//           GlobalMin(x->cgminnodespe)
{
	int l;
	FAMGLeaveInfo n;

	for (l=degree-1; l>=0; l--)
	{
		GetConcentrate(l,&n,sizeof(FAMGLeaveInfo));
		x->coarsefrac += n.coarsefrac;
		x->cgnodes += n.cgnodes;
		x->cgminnodespe = MIN(x->cgminnodespe,n.cgminnodespe);
	}
	Concentrate(x,sizeof(FAMGLeaveInfo));
	Broadcast(x,sizeof(FAMGLeaveInfo));
	x->coarsefrac /= (DOUBLE)procs;

	return;
}
#endif

int FAMGMultiGrid::Construct()
{
    FAMGGrid *g, *cg;
    int level, nnc, nn, ilu, cgilu, leave;
	DOUBLE coarsefrac = 0.0, t, time, cgtime;
	FAMGLeaveInfo myleaveinfo;

    // read parameter
    const int cgnodes = FAMGGetParameter()->Getcgnodes();
#ifdef ModelP
    const int cgminnodespe = FAMGGetParameter()->Getcgminnodespe();
#endif
    const int cglevels = FAMGGetParameter()->Getcglevels();
    const double mincoarse = FAMGGetParameter()->Getmincoarse();
    const int gamma = FAMGGetParameter()->Getgamma();
	if ((strcmp("ilut",FAMGGetParameter()->Getpresmoother()) == 0)
		|| (strcmp("ilut",FAMGGetParameter()->Getpostsmoother()) == 0))  
	{
		ilu = 1; 
	}
	else ilu = 0;
	if (strcmp("ilut",FAMGGetParameter()->Getcgsmoother()) == 0)
	{
		cgilu = 1;
	}
	else cgilu = 0;

	g = grid[0];
#ifdef FAMG_SPARSE_BLOCK
    g->SmoothTV(); 
#else
    g->SmoothTV();
#endif
    FAMGMarkHeap(FAMG_FROM_TOP); // release in Deconstruct
    for(level = 0; level < FAMGMAXGRIDS-1; level++)
    {
		time = CURRENT_TIME;
        // g->SmoothTV();
#ifdef ModelP
		nn  = g->GetNrMasterVectors();	// we are interested only in the master vectors
#else
        nn = g->GetN();
#endif
	
		leave = 0;

		myleaveinfo.coarsefrac = coarsefrac;
		myleaveinfo.cgnodes = nn;
#ifdef ModelP
		myleaveinfo.cgminnodespe = nn;
		GlobalLeave( &myleaveinfo );
#endif

		if( myleaveinfo.coarsefrac > mincoarse )
		{
			leave = 1;
			#ifdef ModelP
			if( me==master )	
			#endif
				cout << "FAMG finished; coarsening rate " << 1.0/myleaveinfo.coarsefrac << " < " << 1.0/mincoarse << endl; 
		}

		if( level >= cglevels )
		{
			leave = 1;
			#ifdef ModelP
			if( me==master )	
			#endif
				cout << "FAMG finished; levels " << level << " >= " << cglevels << endl; 
		}

		if( myleaveinfo.cgnodes <= cgnodes )
		{
			leave = 1;
			#ifdef ModelP
			if( me==master )	
			#endif
				cout << "FAMG finished; cg nodes " << myleaveinfo.cgnodes << " <= " << cgnodes << endl; 
		}

#ifdef ModelP
		if( myleaveinfo.cgminnodespe <= cgminnodespe )
		{
			leave = 1;
			if( me==master )	
				cout << "FAMG finished; min cg nodes per PE " << myleaveinfo.cgminnodespe << " <= " << cgminnodespe << endl; 
		}
#endif

        if (leave)
			break;

        if (gamma < 1) return 0;	// ModelP: simple return because gamma is known to all processors
#ifdef ModelP
//prv(-level,0);
		g->ConstructOverlap();
//prv(-level,0);
		assert(g->GetNrMasterVectors() == nn );
#endif
        g->Stencil();

#ifdef FAMG_ILU
        if(ilu)
        {
            if (g->ILUTDecomp(0)) 
				RETURN(1);
        }
#endif


#ifdef FAMG_SPARSE_BLOCK
        if (g->ConstructDiagonalInverse()) 
			RETURN(1);       
#endif
        if (g->ConstructTransfer()) 
			RETURN(1);       

#ifdef FAMG_SPARSE_BLOCK
        //        if (g->ConstructDiagonalLump()) 
		// 	RETURN(1);       
#endif

		cgtime = CURRENT_TIME;

        nnc = (g->GetN())-(g->GetNF());
        cg = (FAMGGrid *) FAMGGetMem(sizeof(FAMGGrid),FAMG_FROM_TOP);
        if(cg == NULL)
			RETURN(1);
        if(cg->Init(nnc,*g))
			RETURN(1);
        if(cg->Construct(g))
			RETURN(1);
        grid[n] = cg;
        g = cg;
        n++;
//printf("after Galerkin:\n");
//prm(0,0);prm(0,1); prim(0);
//prm(-1,0);
//prv(-level-1,0);

// for debugging: print some consistent and inconsistent matrices:
//GRID *tmpgrid = cg->GetugGrid();
//int tmplevel = GLEVEL(tmpgrid);
//MATDATA_DESC *tmpA = ((FAMGugMatrix*)cg->GetMatrix())->GetMatDesc();
//MATDATA_DESC *tmpACons = ((FAMGugMatrix*)cg->GetConsMatrix())->GetMatDesc();
//prvGeom(tmplevel,0); primGeom(tmplevel+1); prmGeom(tmplevel,MD_SCALCMP(tmpA));
//if (dmatcopy(MYMG(tmpgrid),tmplevel,tmplevel,ALL_VECTORS,tmpACons,tmpA) != NUM_OK) assert(0);
//if (l_matrix_consistent(tmpgrid,tmpACons,MAT_CONS) != NUM_OK) assert(0);
//prvGeom(tmplevel,0); primGeom(tmplevel+1); prmGeom(tmplevel,MD_SCALCMP(tmpACons));

		t = CURRENT_TIME;
		time = t - time;
		cgtime = t - cgtime;

#ifdef ModelP
		cout << me << ": ";
		nnc = cg->GetNrMasterVectors();	// we are interested only in the master vectors
#endif
		coarsefrac = (double)nnc/nn;
		if( nnc <= 0 )
			coarsefrac = 0.000000999;	// dummy
		cout << "amglevel " << -level << " coarsening rate " << 1.0/coarsefrac << " time "<<time<<" cg "<<cgtime<<endl;
    }


    if(level == FAMGMAXGRIDS-1)
    {
        ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": maximum number of levels reached. " << endl;
        FAMGWarning(ostr);
    }
	
    g->Stencil();
#ifdef FAMG_ILU	
    if(cgilu)
    {
        if (g->ILUTDecomp(1)) RETURN(1);
    }
#endif
	
    return 0;
}


int FAMGMultiGrid::Deconstruct()
{
    int i;
    
#ifdef FAMG_REORDERCOLUMN
    if(Reorder()) RETURN(1);
#endif
    for(i = n-1; i >= 0; i--) if(grid[i] != NULL)  grid[i]->Deconstruct(i);
    FAMGReleaseHeap(FAMG_FROM_TOP); // mark in construct
    n = 1;

    return 0;
}
    
int FAMGMultiGrid::Step(int level)
{    
    // Input:  right hand side, initial guess and defect
    // Output: right hand side, new solution and new defect
    
    FAMGGrid *g, *cg;
    int i,n1,n2,gamma;

    n1 = FAMGGetParameter()->Getn1();
    n2 = FAMGGetParameter()->Getn2();
    gamma = FAMGGetParameter()->Getgamma();
	g = grid[level];
	
    if(gamma < 1)
    {

        for(i = 0; i < n1; i++)
        {
            g->PreSmooth();
            // not necessary any more: g->AddVector(FAMGDEFECT,FAMGUNKNOWN);
            g->Defect();
        }        

        for(i = 0; i < n2; i++)
        {
            g->PostSmooth();
            // not necessary any more: g->AddVector(FAMGDEFECT,FAMGUNKNOWN);
            g->Defect();
        }

        return 0;
    }

    if(level == (n-1))
    { 
        if(g->SolveCoarseGrid())
			RETURN(1);
    }
    else
    {
        cg = grid[level+1];

        for(i = 0; i < n1; i++)
        {
            g->PreSmooth();
            // not necessary any more: g->AddVector(FAMGDEFECT,FAMGUNKNOWN);
            g->Defect();
        }        

        // g->DevideFGDefect(); included in the new restriction
        g->Restriction(*(g->GetVector(FAMGUNKNOWN)), *(g->GetVector(FAMGDEFECT)), *(cg->GetVector(FAMGDEFECT)), *(g->GetVector(FAMGUNKNOWN))/*only as dummy!*/ );
        *(cg->GetVector(FAMGUNKNOWN)) = 0.0;
        *(cg->GetVector(FAMGRHS)) = *(cg->GetVector(FAMGDEFECT));
		for(i = 0; i < gamma; i++) 
			if(Step(level+1)) 
				RETURN(1);
        g->Prolongation(cg, *(cg->GetVector(FAMGUNKNOWN)), *(g->GetVector(FAMGUNKNOWN)), *(g->GetVector(FAMGDEFECT)), NULL);
        // g->Defect(); included in the new restriction

        for(i = 0; i < n2; i++)
        {
            g->PostSmooth();
            // not necessary any more: g->AddVector(FAMGDEFECT,FAMGUNKNOWN);
            g->Defect();
        }

    }
	
    return 0;
}


int FAMGMultiGrid::SGSStep(int level)
{    
    // Input:  right hand side, initial guess and defect
    // Output: right hand side, new solution and new defect
    
    FAMGGrid *g;
    g = grid[level];

    g->SGSSmooth();
    // not necessary any more: g->AddVector(FAMGDEFECT,FAMGUNKNOWN);
    g->Defect();
  
    return 0;
}

#ifdef WEG
void FAMGMultiGrid::Mult(double *vout, double *vin)
{
    FAMGMatrix *matrix;

    matrix = grid[0]->GetMatrix();
    matrix->Mult(vout,vin);

    return;
}
#endif

#ifdef FAMG_REORDERCOLUMN
int FAMGMultiGrid::Order()
{
    FAMGGrid *grid0;
    
    grid0 = grid[0];
    if((n > 0) && (grid0 != NULL))
    {
        if(grid0->Order(grid0->GetMap())) RETURN(1);
    }

    return 0;
}

int FAMGMultiGrid::Reorder()
{
    FAMGGrid *grid0;
    
    grid0 = grid[0];
    if((n > 0) && (grid0 != NULL))
    {
        if(grid0->Reorder()) RETURN(1);
    }

    return 0;
}
#endif
