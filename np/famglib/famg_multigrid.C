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

/* RCS_ID
$Header$
*/

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

int FAMGMultiGrid::Construct()
{
    FAMGGrid *g, *cg;
    int level, nnc, nn, ilu, cgilu, leave;
	DOUBLE coarsefrac = 0.0;
#ifdef ModelP
	DOUBLE commBuf [2];
#endif

    // read parameter
    const int cgnodes = FAMGGetParameter()->Getcgnodes();
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
    g->SmoothTV();
    FAMGMarkHeap(FAMG_FROM_TOP); // release in Deconstruct
    for(level = 0; level < FAMGMAXGRIDS-1; level++)
    {
#ifdef ModelP
		nn  = g->GetNrMasterVectors();	// we are interested only in the master vectors
#else
        nn = g->GetN();
#endif
		
#ifdef ModelP
		leave = ( nn <= cgnodes || level>=cglevels) && (gamma > 0);
		commBuf [0] = coarsefrac;
		commBuf[1] = (DOUBLE)leave;
		UG_GlobalMaxNDOUBLE(2,commBuf);
		leave = (commBuf[1]==1.0) || (commBuf[0] > mincoarse) && (gamma > 0);
#else
		leave = ( coarsefrac > mincoarse || nn <= cgnodes || level>=cglevels) && (gamma > 0);
#endif
        if (leave)
			break;

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
        if (gamma < 1) return 0;	// ModelP: simple return because gamma is known to all processors

        if (g->ConstructTransfer()) 
			RETURN(1);       
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

#ifdef ModelP
		cout << me << ": ";
		nnc = cg->GetNrMasterVectors();	// we are interested only in the master vectors
#endif
		coarsefrac = (double)nnc/nn;
		cout << "amglevel " << -level;
		if( nnc > 0 )
		{
			cout << " coarsening rate " << 1.0/coarsefrac << endl;
		}
		else
		{
			cout << " no coarsening" << endl;
		}
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
        g->Restriction(*(g->GetVector(FAMGUNKNOWN)), *(g->GetVector(FAMGDEFECT)), *(cg->GetVector(FAMGDEFECT)));
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
