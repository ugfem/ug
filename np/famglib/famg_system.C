/****************************************************************************/
/*																			*/
/* File:      famg_system.C													*/
/*																			*/
/* Purpose:   famg system class functions									*/
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
#include <fstream.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "famg_system.h"
#include "famg_heap.h"
#include "famg_grid.h"
#include "famg_misc.h"

static FAMGParameter *famgparaptr;
static int famgfirsti=0; // first index (0 or 1 ?)
const int famgnv = 50; // maximum number of Arnoldi vectors

/* RCS_ID
$Header$
*/

// Class  Parameter
void FAMGSetParameter(FAMGParameter *ptr)
{
    famgparaptr = ptr;
}
    
FAMGParameter * FAMGGetParameter()
{
    return famgparaptr;
}
    

int FAMGSystem::Init()
{
    int i;

    nmg = 0;
    for(i = 0; i < FAMGMULTIGRIDS; i++) mg[i] = NULL; 
    matrix = NULL;
    for(i = 0; i < FAMGMAXVECTORS; i++) vector[i] = NULL;
    
    SolverPtr = NULL;
 
    return 0;
    
}

int FAMGSystem::ConstructSimple( FAMGMatrixAlg* stiffmat, FAMGVector* tvA, FAMGVector* tvB )
{
    FAMGMarkHeap(FAMG_FROM_BOTTOM);
    
    matrix = stiffmat;

#ifdef FAMG_REORDERCOLUMN	
    if(matrix->OrderColumns(colmap)) return 1;    
#endif
	
    char *solver = FAMGGetParameter()->Getsolver();
    SolverPtr = &FAMGSystem::BiCGStab;
    if(strcmp(solver,"bicgstab") == 0)
    {
        SolverPtr = &FAMGSystem::BiCGStab;
    }
    else if(strcmp(solver,"linit") == 0)
    {
        SolverPtr = &FAMGSystem::LinIt;
    }
    else if(strcmp(solver,"gmres") == 0)
    {
        SolverPtr = &FAMGSystem::GMRES;
    }
    else
    {
        ostrstream ostr;
        ostr << __FILE__ << __LINE__ <<  "solver = bicgstab" << endl;
        FAMGWarning(ostr);
    }

    // mg0 supposed to be constructed
    FAMGMultiGrid *mg0 = mg[0];
    if(mg0 == NULL)
		return 1;
#ifdef FAMG_REORDERCOLUMN	
    if(mg0->Order()) return 1;
#endif
	
    return 0;
}


int FAMGSystem::Construct( FAMGGridVector *gridvector, FAMGMatrixAlg* stiffmat, FAMGVector *vectors[FAMGMAXVECTORS] )
{
    FAMGMarkHeap(FAMG_FROM_TOP);
    FAMGMarkHeap(FAMG_FROM_BOTTOM);
	int i;
	
	SetGridVector(gridvector);
    if (mygridvector == NULL)
	{
        ostrstream ostr;
        ostr << __FILE__ << __LINE__ <<  "you must provide a gridvector" << endl;
        FAMGError(ostr);
		return 1;
    }
 
#ifdef USE_UG_DS
	SetFineGrid(((FAMGugGridVector*)gridvector)->GetGrid());
#endif
	
    SetMatrix(stiffmat);
    if (GetMatrix() == NULL)
	{
        ostrstream ostr;
        ostr << __FILE__ << __LINE__ <<  "you must provide a matrix" << endl;
        FAMGError(ostr);
		return 1;
    }
 
#ifdef FAMG_REORDERCOLUMN	
    // only for reorder column
    colmap = (int*) FAMGGetMem(famg_interface->n*sizeof(int),FAMG_FROM_TOP);
    if (colmap == NULL) return 1;
    if(matrix->OrderColumns(colmap)) return 1;
#endif
        
	for ( i=0; i<FAMGMAXVECTORS; i++ )
	{
	    if(vectors[i] == NULL)
    	{
        	ostrstream ostr;
	        ostr << __FILE__ << __LINE__ <<  "you must provide vector " << i << endl;
    	    FAMGError(ostr);
			return 1;
	    }
		vector[i] = vectors[i];
	}
	
    char *solver = FAMGGetParameter()->Getsolver();
    SolverPtr = &FAMGSystem::BiCGStab;
    if(strcmp(solver,"bicgstab") == 0)
    {
        SolverPtr = &FAMGSystem::BiCGStab;
    }
    else if(strcmp(solver,"linit") == 0)
    {
        SolverPtr = &FAMGSystem::LinIt;
    }
    else if(strcmp(solver,"gmres") == 0)
    {
        SolverPtr = &FAMGSystem::GMRES;
    }
    else
    {
        ostrstream ostr;
        ostr << __FILE__ << __LINE__ <<  "solver = bicgstab" << endl;
        FAMGWarning(ostr);
    }

    FAMGMultiGrid *mg0 = CreateMultiGrid();
    if(mg0 == NULL)
		return 1;
    if (mg0->Init(*this))
		return 1;
    if (mg0->Construct())
		return 1;

    return 0;
}


int FAMGSystem::Solve(FAMGVector *rhs, FAMGVector *defect, FAMGVector *unknown)
{
    FAMGVector *d, *u;
    int i, u_alloc = 0, d_alloc = 0;
    int status;

    if (rhs == NULL)
	{
        ostrstream ostr;
        ostr << __FILE__ << __LINE__ <<  "you must provide rhs vector" << endl;
        FAMGError(ostr);
		return 1;
    }
    vector[FAMGRHS] = rhs;
    
    FAMGMarkHeap(FAMG_FROM_TOP);

    d = defect;
    if((d == NULL) || (d == rhs))
    {
        d = rhs->create_new();
        if (d == NULL)
			return 1;
		d_alloc = 1;
    }
    vector[FAMGDEFECT] = d;
 
    u = unknown;
    if((u == NULL) || (u == rhs) || (u == defect) )
    {
        u = rhs->create_new();
        if (u == NULL) 
			return 1;
		u_alloc = 1;
    }
    vector[FAMGUNKNOWN] = u;
    
    FAMGGrid *g = mg[0]->GetGrid(0);
    g->SetVector(FAMGUNKNOWN,u);
    g->SetVector(FAMGDEFECT,d);
    g->SetVector(FAMGRHS,rhs);
#ifdef FAMG_REORDERCOLUMN	
    if(g->OrderVector(FAMGRHS,g->GetMap())) return 1;
#endif
	
	status = (this->*SolverPtr)();
	if(status != 0)
		return status;
       
#ifdef FAMG_REORDERCOLUMN	
	if(g->ReorderVector(FAMGRHS,g->GetMap())) return 1;
	if(g->ReorderVector(FAMGDEFECT,g->GetMap())) return 1;
	if(g->ReorderVector(FAMGUNKNOWN,g->GetMap())) return 1;
#endif

	if(defect != d) 
	{
		if (defect != NULL)
			*defect = *d;
	}
	if(d_alloc)
	{
		delete d;
		d = NULL;
		vector[FAMGDEFECT] = d;
		g->SetVector(FAMGDEFECT,d);
	}
	
	if(unknown != u) 
	{ 
		if (unknown != NULL)
			*unknown = *u;
	}
	if(u_alloc)
	{
		delete u;
		u = NULL;
		vector[FAMGUNKNOWN] = u;
		g->SetVector(FAMGUNKNOWN,u);
	}
    
    FAMGReleaseHeap(FAMG_FROM_TOP);

    return status;
}

int FAMGSystem::Deconstruct()
{
    FAMGMultiGrid *mg0=mg[0];

    if(mg0 != NULL)
		if(mg0->Deconstruct())
			return 1;

    // repair matrix
#ifdef FAMG_REORDERCOLUMN	
    if(matrix->ReorderColumns(colmap)) return 1;
#endif
    // ??? matrix->RemodifyIndex(famgfirsti);
	
    FAMGReleaseHeap(FAMG_FROM_BOTTOM);
    FAMGReleaseHeap(FAMG_FROM_TOP);

    return 0;
}

int FAMGSystem::DeconstructSimple()
{
    FAMGMultiGrid *mg0=mg[0];

#ifdef FAMG_REORDERCOLUMN	
    if(mg0 != NULL) if(mg0->Reorder()) return 1;

    // repair matrix
    if(matrix->ReorderColumns(colmap)) return 1;
#endif
    // ??? matrix->RemodifyIndex(colmap,famgfirsti);
    FAMGReleaseHeap(FAMG_FROM_BOTTOM);
    return 0;
}

FAMGSystem::FAMGSystem()
{
    int i;

    nmg = 0;
    for(i = 0; i < FAMGMULTIGRIDS; i++) mg[i] = NULL; 
    matrix = NULL;
    for(i = 0; i < FAMGMAXVECTORS; i++) vector[i] = NULL;
    SolverPtr = &FAMGSystem::LinIt;
}


FAMGMultiGrid *FAMGSystem::CreateMultiGrid()
{
    FAMGMultiGrid *newmg;
   
    newmg = (FAMGMultiGrid *) FAMGGetMem(sizeof(FAMGMultiGrid),FAMG_FROM_TOP);
    if(newmg == NULL)
		return(NULL);
    mg[nmg] = newmg;
    nmg++;
    
    return newmg;
}

#ifdef WRITE_VECTOR_TEST
static void FAMGWriteVector(int n, double *vec)
{
    ofstream vfile("vec.dat",ios::out);
    if (!vfile)
		return;
    for(int i = 0; i < n; i++)
    {
        vfile << i << "\t" << vec[i] << endl;
    }

    return;
}
#endif


int FAMGSystem::LinIt()
{
    double rlimit,alimit,divlimit,reduction,limit,defectnorm,startdefect,oldnorm;
    int maxit,i;
    ostrstream ostr;

    maxit = famgparaptr->Getmaxit();
    rlimit = famgparaptr->Getrlimit();
    alimit = famgparaptr->Getalimit();
    divlimit =famgparaptr->Getdivlimit();
    reduction = famgparaptr->Getreduction();
	FAMGVector &defect = *GetVector(FAMGDEFECT);
	FAMGVector &rhs = *GetVector(FAMGRHS);
	FAMGVector &solution = *GetVector(FAMGUNKNOWN);
	
    solution = 0.0;
	defect = rhs;
    startdefect = defect.norm();
    limit = rlimit*startdefect;
    defectnorm = startdefect;
    ostr << 0 << "\t" <<  startdefect << endl;
    FAMGWrite(ostr);

    // FAMGMultiGrid * mg0 = CreateMultiGrid();
    // if(mg0 == NULL) return 1;
    // if (mg0->Init(*this)) return 1;
    // if (mg0->Construct()) return 1;
    FAMGMultiGrid *mg0 = mg[0];
   
    for(i = 0; i < maxit; i++)
    {
        if((defectnorm < alimit) || (defectnorm < limit)) 
			break;
        if(mg0->Step(0)) 
		{
			mg0->Deconstruct();
			return 1;
		}
        oldnorm = defectnorm;
        defectnorm = defect.norm();
        ostr << i+1 << "\t" << defectnorm << "\t" << defectnorm/oldnorm << endl;    
        FAMGWrite(ostr);
    
        if((defectnorm < alimit) || (defectnorm < limit)) 
			break;
        if (defectnorm/oldnorm > divlimit)
        {
            ostr << "diverged" << endl;
            FAMGWrite(ostr);
            break;
        }
        
        if((defectnorm < alimit) || (defectnorm < limit)) break;
            
    }

    if (defectnorm < startdefect*reduction)
    {
        // Yeah, problem solved !!!

        return 0;
    }
    
    return 1;
}

int FAMGSystem::AdTVSolve()
{
    double rlimit,alimit,divlimit,reduction,limit,defectnorm,startdefect,oldnorm;
    int maxit,i,j,init;
    ostrstream ostr;

    maxit = famgparaptr->Getmaxit();
    rlimit = famgparaptr->Getrlimit();
    alimit = famgparaptr->Getalimit();
    divlimit =famgparaptr->Getdivlimit();
    reduction = famgparaptr->Getreduction();
    init = 8;
	FAMGVector &defect = *GetVector(FAMGDEFECT);
	FAMGVector &rhs = *GetVector(FAMGRHS);
	FAMGVector &solution = *GetVector(FAMGUNKNOWN);
	FAMGVector &tva = *GetVector(FAMGTVA);
	FAMGVector &tvb = *GetVector(FAMGTVB);

	solution = 0.0;
	defect = rhs;
	startdefect = defect.norm();
    limit = rlimit*startdefect;
    defectnorm = startdefect;
    ostr << 0 << "\t" <<  startdefect << endl;
    FAMGWrite(ostr);

    // FAMGMultiGrid * mg0 = CreateMultiGrid();
    // if(mg0 == NULL) return 1;
    // if (mg0->Init(*this)) return 1;
    // if (mg0->Construct()) return 1;
    FAMGMultiGrid *mg0 = mg[0];
   
    for(i = 0; i < maxit;)
    {
        if((defectnorm < alimit) || (defectnorm < limit))
			break;
        for(j = 0; (j < init) && (i < maxit); j++, i++)
        {
			tva = solution;
			tvb = solution;
            if(mg0->Step(0))
			{
				mg0->Deconstruct();
				return 1;
			
			}
			tva -= solution;
			tvb -= solution;
            oldnorm = defectnorm;
            defectnorm = defect.norm();
            ostr << i+1 << "\t" << defectnorm << "\t" << defectnorm/oldnorm << endl;    
            FAMGWrite(ostr);
            
            if((defectnorm < alimit) || (defectnorm < limit))
				break;
            if (defectnorm/oldnorm > divlimit)
            {
                ostr << "diverged" << endl;
                FAMGWrite(ostr);
                break;
            }


        }
        if (defectnorm/oldnorm > divlimit)
        {
            ostr << "diverged" << endl;
            FAMGWrite(ostr);
            break;
        }

        if((defectnorm < alimit) || (defectnorm < limit))
			break;
        if(i == maxit)
			break;

        if(mg0->Deconstruct())
			return 1;
        if (mg0->Construct())
			return 1;
    }
    
    if (defectnorm < startdefect*reduction)
    {
        // Yeah, problem solved !!!

        return 0;
    }
    
    return 1;
}




// BiCGStab

int FAMGSystem::BiCGStab()
{
	assert(0); // to port
#ifdef FAMG_BICG
    double rlimit,alimit,divlimit,reduction,limit,defectnorm,startdefect,oldnorm;
    double *vec[6], rho, oldrho, alpha, beta, omega, nenner;
    int maxit,i;
    ostrstream ostr;


    const int FAMGX = 0;
    const int FAMGR = 1;
    const int FAMGV = 2;
    const int FAMGP = 3;
    const int FAMGT = 4;
    const int FAMGR0 = 5;

    FAMGMarkHeap(FAMG_FROM_BOTTOM);
    vec[FAMGX] = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_BOTTOM);
    if(vec[FAMGX] == NULL)
		return 1;
    vec[FAMGR] = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_BOTTOM);
    if(vec[FAMGR] == NULL)
		return 1;
    vec[FAMGV] = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_BOTTOM);
    if(vec[FAMGV] == NULL)
		return 1;
    vec[FAMGP] = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_BOTTOM);
    if(vec[FAMGP] == NULL)
		return 1;
    vec[FAMGT] = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_BOTTOM);
    if(vec[FAMGT] == NULL)
		return 1;
    vec[FAMGR0] = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_BOTTOM);
    if(vec[FAMGR0] == NULL)
		return 1;

    maxit = famgparaptr->Getmaxit();
    rlimit = famgparaptr->Getrlimit();
    alimit = famgparaptr->Getalimit();
    divlimit =famgparaptr->Getdivlimit();
    reduction = famgparaptr->Getreduction();


    // construct preconditioner
    // FAMGMultiGrid * mg0 = CreateMultiGrid();
    // if(mg0 == NULL) return 1;
    // if (mg0->Init(*this)) return 1;
    // if (mg0->Construct()) return 1;
    FAMGMultiGrid *mg0 = mg[0];
   
    FAMGSetVector(n,vec[FAMGX],0.0);
    FAMGCopyVector(n,vec[FAMGR0],vector[FAMGRHS]);

    
    defectnorm = FAMGNorm(n,vec[FAMGR0]);
    startdefect = defectnorm;
    limit = rlimit*startdefect;
    ostr << 0 << "\t" <<  startdefect << endl;
    FAMGWrite(ostr);

    FAMGCopyVector(n,vec[FAMGR],vec[FAMGR0]);
    FAMGCopyVector(n,vec[FAMGP],vec[FAMGR]);
    rho = FAMGSum(n,vec[FAMGR]);
    if (Abs(rho) < 1e-10*alimit) 
    {
       ostr << __FILE__ << ", line " << __LINE__ << ": rho too small" << endl;
       FAMGWarning(ostr);
    }

    for(i = 0; i < maxit; i++)
    {

        FAMGCopyVector(n,vector[FAMGRHS],vec[FAMGP]);
        FAMGSetVector(n,vector[FAMGUNKNOWN],0.0);
        FAMGCopyVector(n,vector[FAMGDEFECT],vector[FAMGRHS]);
        if(mg0->Step(0))
		{
			mg0->Deconstruct();
			return 1;
		}
        FAMGSetSubVector(n,vec[FAMGV],vec[FAMGP],vector[FAMGDEFECT]);

        nenner = FAMGSum(n,vec[FAMGV]);
        if (Abs(nenner) < 1e-15*Abs(rho)) 
        {
            ostr << __FILE__ << ", line " << __LINE__ << ": nenner too small" << endl;
            FAMGWarning(ostr);
        }
    
        alpha = rho/nenner;
        FAMGAddVector(n,vec[FAMGX],vector[FAMGUNKNOWN],alpha);
        FAMGAddVector(n,vec[FAMGR],vec[FAMGV],-alpha);

        oldnorm = defectnorm;
        defectnorm = FAMGNorm(n,vec[FAMGR]);
        ostr << i+0.5 << "\t" << defectnorm << "\t" << defectnorm/oldnorm;
        ostr << "\t" << alpha  << endl; 
        FAMGWrite(ostr);
        if((defectnorm < alimit) || (defectnorm < limit)) break;
        if (defectnorm/oldnorm > divlimit)
        {
            ostr << "diverged" << endl;
            FAMGWrite(ostr);
           break;
        }
        
        FAMGCopyVector(n,vector[FAMGRHS],vec[FAMGR]);
        FAMGSetVector(n,vector[FAMGUNKNOWN],0.0);
        FAMGCopyVector(n,vector[FAMGDEFECT],vector[FAMGRHS]);
        if(mg0->Step(0))
		{
			mg0->Deconstruct();
			return 1;
		}
        FAMGSetSubVector(n,vec[FAMGT],vec[FAMGR],vector[FAMGDEFECT]);

        omega = FAMGScalProd(n,vec[FAMGT],vec[FAMGR])/FAMGScalProd(n,vec[FAMGT],vec[FAMGT]);
        if (Abs(omega) < 1e-15) 
        {
            ostr << __FILE__ << ", line " << __LINE__ << ": omega too small" << endl;
            FAMGWarning(ostr);
        }

        FAMGAddVector(n,vec[FAMGX],vector[FAMGUNKNOWN],omega);
        FAMGAddVector(n,vec[FAMGR],vec[FAMGT],-omega);
        FAMGAddVector(n,vec[FAMGP],vec[FAMGV],-omega);

        oldnorm = defectnorm;
        defectnorm = FAMGNorm(n,vec[FAMGR]);
        ostr << i+1 << "\t" << defectnorm << "\t" << defectnorm/oldnorm;
        ostr << "\t" << omega << endl;   
        FAMGWrite(ostr);
        if((defectnorm < alimit) || (defectnorm < limit)) break;
        if (defectnorm/oldnorm > divlimit)
        {
            ostr << "diverged" << endl;
            FAMGWrite(ostr);            
            break;
        }

        oldrho = rho;
        rho = FAMGSum(n,vec[FAMGR]);
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
    FAMGCopyVector(n,vector[FAMGUNKNOWN],vec[FAMGX]);
    FAMGCopyVector(n,vector[FAMGRHS],vec[FAMGR0]);

    FAMGReleaseHeap(FAMG_FROM_BOTTOM);
    
    if (defectnorm < startdefect*reduction)
    {
        // Yeah, problem solved !!!

        return 0;
    }
#endif
	
    return 1;
}

    // GMRES


void FAMGMatVecMult(int n, int nr, double* M, double* h)
{
    double sum, *g;
    int i,j;

    
    FAMGMarkHeap(FAMG_FROM_BOTTOM);
    g = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_BOTTOM);

    for(i = 0; i < n; i++)
    {
        sum = 0.0;
        for(j = 0; j < n; j++)
        {
             sum += M[i*nr+j]*h[j];
        }
        g[i] = sum;
    }
    for(i = 0; i < n; i++) h[i] = g[i];
            
    FAMGReleaseHeap(FAMG_FROM_BOTTOM);

    return;
}

#ifdef FAMG_GMRES
static void UpdateQ(int n, int nr, double *Q, double c, double s)
{
    int row, srow, j;

    row =n*nr;
    srow = (n+1)*nr;
    for(j = 0; j <= n; j++)
    {
        Q[srow+j] = -s*Q[row+j];
        Q[row+j] = c*Q[row+j];
    }
    Q[srow+n+1] = c;
    Q[row+n+1] = s;

    return;
}

static void UpdateH(int j, int nr, double *H, double *h)
{
    int i;

    for(i = 0; i <= j; i++)
    {
        H[i*nr+j] = h[i];
    }
    
    return;
}
        
static void ComputeX(int n, int nv, double *x, double *H, double *Q, double q)
{
    double sum;
    int i,j;

    for(i = 0; i < n; i++)
    {
        x[i] = Q[i*(nv+1)]*q;
    }

    for(i = n-1; i >= 0; i--)
    {
        sum = x[i];
        for(j = i+1; j < n; j++)
        {
            sum -= H[i*nv+j]*x[j];
        }
        x[i] = sum/H[i*nv+i];
    }
    
    return;
}


static double ComputeEV(int n, int nv, double *H, double *G, double *e)
{
    double sum, norm, s[100], sig;
    int i, j, k, row;

    // initial guess
    sig = -1;
    for(i = 0; i < n; i++) 
    {
        e[i] = sig;
        sig = -sig;
    }

    // test

    for(k = 0; k < 100; k++)
    {
        // Invert G^T. The invers is G since G G^T = I
        for(i = 0; i < n; i++)
        {
            sum = 0.0;
            row = i*(nv+1);
            for(j = 0; j < n; j++) 
            {
                sum += G[row+j]*e[j];
            }
            s[i] = sum;
        }
    
        // Invert H. H is upper triangular
        norm = 0.0;
        for(i = n-1; i >= 0; i--)
        {
            sum = s[i];
            row = i*nv;
            for(j = i+1; j < n; j++)
            {
                sum -= H[row+j]*s[j];
            }
            s[i] = sum/H[row+i];
            norm += s[i]*s[i];
        }
        
        norm = sqrt(norm);
        for(i = 0; i < n; i++) e[i] = s[i]/norm;
    }

    return 1.0/norm;
}
#endif

int FAMGSystem::Arnoldi(FAMGMultiGrid *mg0, double **vec, double *H, double *G, double *Q, double *P, double &q0, int con)
{
	assert(0); // to port
#ifdef FAMG_GMRES
	double q, gamma, oldgamma, h[famgnv], c, s, nenner;
     int i,j,nv;
     ostrstream ostr;
           
     nv = famgparaptr->Getnv();  
     for(i = 0; i < (nv+1)*(nv+1); i++) Q[i] = 0.0;

     q = FAMGNorm(n,vector[FAMGDEFECT]);
     FAMGCopyScaledVector(n,vec[0],vector[FAMGDEFECT],1.0/q); 
     q0 = q;
     Q[0] = 1.0;
     oldgamma = q0;

     // Arnoldi process
     for(i = 0; i < nv; i++)
     {
         // def = M^{-1}A v_{i-1}
         FAMGCopyScaledVector(n,vec[i],vector[FAMGDEFECT],1.0/q); 
         FAMGCopyVector(n,vector[FAMGRHS],vec[i]);
         FAMGSetVector(n,vector[FAMGUNKNOWN],0.0);
         FAMGCopyVector(n,vector[FAMGDEFECT],vector[FAMGRHS]);
         if(con)
         {
             if(mg0->Step(0))
			 {
			 	mg0->Deconstruct();
			 	return 1;
			 }
         }
         else
         {
             mg0->SGSStep(0);
         }
         FAMGSetSubVector(n,vector[FAMGDEFECT],vector[FAMGRHS],vector[FAMGDEFECT]);
 
         // modified Gram-Schmitt
         for(j = 0; j <= i; j++)
         {
             h[j] = FAMGScalProd(n,vector[FAMGDEFECT],vec[j]);
             FAMGAddVector(n,vector[FAMGDEFECT],vec[j],-h[j]);
         }

         q = FAMGNorm(n,vector[FAMGDEFECT]);
         if(q < 1e-10) 
         {
             ostr << __FILE__ << ", line " << __LINE__ << ": GMRES breakdownl" << endl;
             FAMGWarning(ostr);
         }
         // => def = v_i*q

         // update H
         FAMGMatVecMult(i+1,nv+1,Q,h);
         UpdateH(i,nv,H,h);

         if(i == nv-1)
         {
             for(j = 0; j < (nv*nv); j++) G[j] = H[j];
             for(j = 0; j < ((nv+1)*(nv+1)); j++) P[j] = Q[j];
         }
           

         // eliminate H_{i+1,i} = q
         nenner = sqrt(h[i]*h[i] + q*q);
         c = h[i]/nenner; s = q/nenner;
         UpdateQ(i,nv+1,Q,c,s);
        
         // only H_{i,i} must be modified at this point
         H[nv*i+i] = nenner;

         // check convergence
         gamma = Abs(q0*Q[(nv+1)*(i+1)]);
         ostr << i+1 << "\t" << gamma << "\t" << gamma/oldgamma << endl;
         FAMGWrite(ostr);
         oldgamma = gamma;

     }
#endif

     return 0;
}

int FAMGSystem::UpdateSolution(FAMGMultiGrid *mg0, double **vec, double *H, double *Q, double &q0, int con)
{
	assert(0); // to port
#ifdef FAMG_GMRES
    double x[famgnv];
    int j, nv;
    
    nv = famgparaptr->Getnv();

    ComputeX(nv,nv,x,H,Q,q0);
    FAMGSetVector(n,vector[FAMGRHS],0.0);
    for(j = 0; j < nv; j++)
    {
        FAMGAddVector(n,vector[FAMGRHS],vec[j],x[j]);
    }
            

    FAMGSetVector(n,vector[FAMGUNKNOWN],0.0);
    FAMGCopyVector(n,vector[FAMGDEFECT],vector[FAMGRHS]);
    
    
    if(con)
    {
        if(mg0->Step(0))
		{
			mg0->Deconstruct(); 
			return 1;
		}
    }
    else
    {
        mg0->SGSStep(0);
    }
#endif

    return 0;
}

int FAMGSystem::ComputeEigenVector(FAMGMultiGrid *mg0, double **vec, double *G, double *P, int con)
{
	assert(0); // to port
#ifdef FAMG_GMRES
    double e[famgnv], ew, norm;
    int j, nv;

    nv = famgparaptr->Getnv();

    // compute eigenvector to smallest eigenvalue
    ew = ComputeEV(nv,nv,G,P,e);
    ostrstream ostr; ostr << "ew = " << ew << endl;
    FAMGWrite(ostr);
    if(Abs(ew) < 0.1)
    {
        FAMGSetVector(n,vector[FAMGRHS],0.0);
        FAMGSetVector(n,vector[FAMGUNKNOWN],0.0);
        for(j = 0; j < nv; j++)
        {
            FAMGAddVector(n,vector[FAMGRHS],vec[j],e[j]);
        }


        FAMGSetVector(n,vector[FAMGUNKNOWN],0.0);
        FAMGCopyVector(n,vector[FAMGDEFECT],vector[FAMGRHS]);
        if(con)
        {
            if(mg0->Step(0))
			{
				mg0->Deconstruct();
				return 1;
			}
        }
        else
        {
            mg0->SGSStep(0);
        }
        norm = FAMGNorm(n,vector[FAMGUNKNOWN]);
        FAMGCopyScaledVector(n,vector[FAMGTVA],vector[FAMGUNKNOWN],1.0/norm);
        FAMGCopyScaledVector(n,vector[FAMGTVB],vector[FAMGUNKNOWN],1.0/norm);

        return -1;
    }
#endif

    return 0;
}
            
int FAMGSystem::GMRES()
{
	assert(0); // to port
#ifdef FAMG_GMRES
    double rlimit,alimit,reduction,limit,defectnorm,startdefect;
    int maxit;
    double H[famgnv*famgnv], G[famgnv*famgnv], Q[(famgnv+1)*(famgnv+1)], P[(famgnv+1)*(famgnv+1)], *vec[famgnv], *sol, *rhs, *def, q0;
    int i, k, newtv, con, nv;
    ostrstream ostr;

    maxit = famgparaptr->Getmaxit();
    rlimit = famgparaptr->Getrlimit();
    alimit = famgparaptr->Getalimit();
    reduction = famgparaptr->Getreduction();
    nv = famgparaptr->Getnv();

    FAMGMarkHeap(FAMG_FROM_BOTTOM);

    sol = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_BOTTOM);
    if(sol == NULL)
		return 1;
    rhs = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_BOTTOM);
    if(rhs == NULL)
		return 1;
    
    // test
    def = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_BOTTOM);
    if(def == NULL)
		return 1;

    for(i = 0; i < nv; i++)
    {
        vec[i] = (double *) FAMGGetMem(n*sizeof(double),FAMG_FROM_BOTTOM);
        if(vec[i] == NULL)
			return 1;
    }
    

    FAMGSetVector(n,vector[FAMGUNKNOWN],0.0);
    FAMGSetVector(n,sol,0.0);
    FAMGCopyVector(n,vector[FAMGDEFECT],vector[FAMGRHS]);

    startdefect = defectnorm = FAMGNorm(n,vector[FAMGDEFECT]);
    limit = rlimit*defectnorm;
    ostr << "start defect: " << defectnorm  << endl;
    FAMGWrite(ostr);


    // FAMGMultiGrid * mg0 = CreateMultiGrid();
    // if(mg0 == NULL) return 1;
    // if (mg0->Init(*this)) return 1;
    // if (mg0->Construct()) return 1;
    FAMGMultiGrid *mg0 = mg[0];

    newtv = 0;
    con = 1;

    for(k = 0; k*nv < maxit; k++)
    {

        // construct preconditioner
        if(newtv)
        {
            if(con) 
            {
                if(mg0->Deconstruct())
					return 1;
            }
            if (mg0->Construct())
				return 1;
            con = 1;
       }

 
        FAMGCopyVector(n,rhs,vector[FAMGRHS]);
        FAMGCopyVector(n,sol,vector[FAMGUNKNOWN]);
        //test
        FAMGCopyVector(n,def,vector[FAMGDEFECT]);

        if(Arnoldi(mg0,vec,H,G,Q,P,q0,con))
			return 1;

        newtv = 0;
        // status = ComputeEigenVector(mg0,vec,G,P,con);
        // if(status > 0) return 1; // error
        // if(status < 0) 
        // {
        //     newtv = 1; // restart with new TV
        // }


        if(UpdateSolution(mg0,vec,H,Q,q0,con))
			return 1;

        FAMGAddVector(n,vector[FAMGUNKNOWN],sol);

        FAMGCopyVector(n,vector[FAMGRHS],rhs);
        mg0->GetGrid(0)->Defect();
        defectnorm = FAMGNorm(n,vector[FAMGDEFECT]);
        ostr << defectnorm << endl;
        FAMGWrite(ostr);
        if( (defectnorm < alimit) || (defectnorm < limit)) break;
    }

    if(!con) 
    {
        ostr  << __FILE__ << __LINE__ << endl;
        FAMGError(ostr);
        return 1;
    }

    FAMGReleaseHeap(FAMG_FROM_BOTTOM);

    if (defectnorm < startdefect*reduction)
    {
        // Yeah, problem solved !!!

        return 0;
    }
#endif

    return 1;
}

