/****************************************************************************/
/*																			*/
/* File:      multigrid.C													*/
/*																			*/
/* Purpose:   cmg   multigrid classes functions								*/
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

#include "misc.h"
#include "grid.h"
#include "multigrid.h"
#include "graph.h"
#include "matrix.h"
#include "heap.h"
#include "system.h"

/* RCS_ID
$Header$
*/

// Class CMGMultigrid


int CMGMultiGrid::Init(const CMGSystem &system)
{
    CMGGrid *grid0;
    int i;

    n = 0;
    for(i = 0; i < CMGMAXGRIDS; i++) grid[i] = NULL;

    grid0 = (CMGGrid *) CMGGetMem(sizeof(CMGGrid),CMG_FROM_TOP);
    if(grid0 == NULL) return(NULL);

    if( grid0->InitLevel0(system)) return 1;
    grid[n] = grid0;
    n++;

    return 0;
}    

int CMGMultiGrid::Construct()
{
    CMGGrid *g, *cg;
    int level, nnc, nn, ilu, cgilu;

    // read parameter
    const int cgnodes = CMGGetParameter()->Getcgnodes();
    const double mincoarse = CMGGetParameter()->Getmincoarse();
    const int gamma = CMGGetParameter()->Getgamma();
    if ((strcmp("ilut",CMGGetParameter()->Getpresmoother()) == 0)
     || (strcmp("ilut",CMGGetParameter()->Getpostsmoother()) == 0))  
    {
        ilu = 1; 
        
    }
    else ilu = 0;
    if (strcmp("ilut",CMGGetParameter()->Getcgsmoother()) == 0)
    {
        cgilu = 1;
        
    }
    else cgilu = 0;

    g = grid[0];
    g->SmoothTV();
    CMGMarkHeap(CMG_FROM_TOP); // release in Deconstruct
    for(level = 0; level < CMGMAXGRIDS-1; level++)
    {
        g->Stencil();
        nn = g->GetN();
        if ((nn <= cgnodes) && (gamma > 0)) break;
        
        if(ilu)
        {
            if (g->ILUTDecomp(0)) return 1;
        }
        if (gamma < 1) return 0;

        if (g->ConstructTransfer()) return 1;        
        nnc = (g->GetN())-(g->GetNF());
        // if (nnc == 0) return 0; 
        cg = (CMGGrid *) CMGGetMem(sizeof(CMGGrid),CMG_FROM_TOP);
        if(cg == NULL) return(NULL);
        if(cg->Init(nnc)) return 1;
        if(cg->Construct(g)) return 1;
        grid[n] = cg;
        g = cg;
        n++;
        if(nnc > nn*mincoarse) break;
    }

    if(level == CMGMAXGRIDS-1)
    {
        ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": maximum number of levels reached. " << endl;
        CMGWarning(ostr);
    }

    g->Stencil();
    if(cgilu)
    {
        if (g->ILUTDecomp(1)) return 1;
    }

    return 0;
}


int CMGMultiGrid::Deconstruct()
{
    int i;
    
    if(Reorder()) return 1;
    for(i = 0; i < n; i++) if(grid[i] != NULL)  grid[i]->Deconstruct();
    CMGReleaseHeap(CMG_FROM_TOP); // mark in construct
    n = 1;

    return 0;
}
    
int CMGMultiGrid::Step(int level)
{    
    // Input:  right hand side, initial guess and defect
    // Output: right hand side, new solution and new defect
    
    CMGGrid *g, *cg;
    int i,n1,n2,gamma;

    n1 = CMGGetParameter()->Getn1();
    n2 = CMGGetParameter()->Getn2();
    gamma = CMGGetParameter()->Getgamma();
 

    if(gamma < 1)
    {

        g = grid[level];
        for(i = 0; i < n1; i++)
        {
            g->PreSmooth();
            g->AddVector(CMGDEFECT,CMGUNKNOWN);
            g->Defect();
        }        

        for(i = 0; i < n2; i++)
        {
            g->PostSmooth();
            g->AddVector(CMGDEFECT,CMGUNKNOWN);
            g->Defect();
        }

        return 0;
    }

    g = grid[level];
    if(level == (n-1))
    { 
        if(g->BiCGStab()) return 1;
    }
    else
    {
        cg = grid[level+1];

        for(i = 0; i < n1; i++)
        {
            g->PreSmooth();
            g->AddVector(CMGDEFECT,CMGUNKNOWN);
            g->Defect();
        }        

        g->DevideFGDefect();
        g->Restriction(cg);
        cg->SetVector(CMGUNKNOWN,0.0);
        cg->CopyVector(CMGRHS,CMGDEFECT);
        for(i = 0; i < gamma; i++) { if(Step(level+1)) return 1; }
        g->Prolongation(cg);
        g->Defect();

        for(i = 0; i < n2; i++)
        {
            g->PostSmooth();
            g->AddVector(CMGDEFECT,CMGUNKNOWN);
            g->Defect();
        }

    }

    return 0;
}


int CMGMultiGrid::SGSStep(int level)
{    
    // Input:  right hand side, initial guess and defect
    // Output: right hand side, new solution and new defect
    
    CMGGrid *g;
    g = grid[level];

    g->SGSSmooth();
    g->AddVector(CMGDEFECT,CMGUNKNOWN);
    g->Defect();
  
    return 0;
}

void CMGMultiGrid::Mult(double *vout, double *vin)
{
    CMGMatrix *matrix;

    matrix = grid[0]->GetMatrix();
    matrix->Mult(vout,vin);

    return;
}
        
int CMGMultiGrid::Order()
{
    CMGGrid *grid0;
    
    grid0 = grid[0];
    if((n > 0) && (grid0 != NULL))
    {
        if(grid0->Order(grid0->GetMap())) return 1;
    }

    return 0;
}

int CMGMultiGrid::Reorder()
{
    CMGGrid *grid0;
    
    grid0 = grid[0];
    if((n > 0) && (grid0 != NULL))
    {
        if(grid0->Reorder()) return 1;
    }

    return 0;
}

    
 
