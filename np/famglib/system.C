/****************************************************************************/
/*																			*/
/* File:      system.C														*/
/*																			*/
/* Purpose:   cmg system class functions									*/
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
#include "system.h"
#include "heap.h"
#include "grid.h"
#include "misc.h"

static CMGParameter *cmgparaptr;
static int cmgfirsti=0; // first index (0 or 1 ?)
const int cmgnv = 50; // maximum number of Arnoldi vectors

/* RCS_ID
$Header$
*/

// Class  Parameter
void CMGSetParameter(CMGParameter *ptr)
{
    cmgparaptr = ptr;
}
    
CMGParameter * CMGGetParameter()
{
    return cmgparaptr;
}
    

int CMGSystem::Init()
{
    int i;

    nmg = 0;
    for(i = 0; i < CMGMULTIGRIDS; i++) mg[i] = NULL; 
    matrix = NULL;
    for(i = 0; i < CMGMAXVECTORS; i++) vector[i] = NULL;
    extra = NULL;
    
    SolverPtr = NULL;
 
    return 0;
    
}
int CMGSystem::ConstructSimple(double *entry, int *index, int *start, int nn, int nl, void **extraptr)
{
    CMGMarkHeap(CMG_FROM_BOTTOM);
    
    n = nn;
    matrix->SetN(n);
    matrix->SetNL(nl);
    matrix->SetIndex(index);
    matrix->SetEntry(entry);
    matrix->SetStartPtr(start);

    matrix->ModifyIndex(colmap,cmgfirsti);
    if(matrix->ConstructEndB()) return 1;
    if(matrix->OrderColumns(colmap)) return 1;
    if(matrix->ConstructAdjoinedB()) return 1;
    

    char *solver = CMGGetParameter()->Getsolver();
    SolverPtr = &CMGSystem::BiCGStab;
    if(strcmp(solver,"bicgstab") == 0)
    {
        SolverPtr = &CMGSystem::BiCGStab;
    }
    else if(strcmp(solver,"linit") == 0)
    {
        SolverPtr = &CMGSystem::LinIt;
    }
    else if(strcmp(solver,"gmres") == 0)
    {
        SolverPtr = &CMGSystem::GMRES;
    }
    else
    {
        ostrstream ostr;
        ostr << __FILE__ << __LINE__ <<  "solver = bicgstab" << endl;
        CMGWarning(ostr);
    }

    extra = extraptr;
    
    // mg0 supposed to be constructed
    CMGMultiGrid *mg0 = mg[0];
    if(mg0 == NULL) return 1;
    if(mg0->Order()) return 1;

    return 0;
}


int CMGSystem::Construct(double *entry, int *index, int *start, int nn, int nl, double *tvA, double *tvB, void **extraptr)
{
    int i;

    CMGMarkHeap(CMG_FROM_TOP);
    CMGMarkHeap(CMG_FROM_BOTTOM);

    matrix = (CMGMatrix *) CMGGetMem(sizeof(CMGMatrix),CMG_FROM_TOP);
    if (matrix == NULL) return 1;
    
    // only for reorder column
    colmap = (int*) CMGGetMem(nn*sizeof(int),CMG_FROM_TOP);
    if (colmap == NULL) return 1;
    
    n = nn;
    matrix->SetN(n);
    matrix->SetNL(nl);
    matrix->SetIndex(index);
    matrix->SetEntry(entry);
    matrix->SetStartPtr(start);

    cmgfirsti = matrix->GetSmallestIndex();
    matrix->ModifyIndex(cmgfirsti);
    if(matrix->ConstructEndB()) return 1;
    if(matrix->OrderColumns(colmap)) return 1;
    if(matrix->ConstructAdjoinedB()) return 1;
    
    if(tvA == NULL)
    {
        tvA = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_TOP);
        for(i = 0; i < n; i++) tvA[i] = 1.0;
    }
    if(tvB == NULL)
    {
        tvB = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_TOP);
        for(i = 0; i < n; i++) tvB[i] = 1.0;
    }
    vector[CMGTVA] = tvA;
    vector[CMGTVB] = tvB;

    char *solver = CMGGetParameter()->Getsolver();
    SolverPtr = &CMGSystem::BiCGStab;
    if(strcmp(solver,"bicgstab") == 0)
    {
        SolverPtr = &CMGSystem::BiCGStab;
    }
    else if(strcmp(solver,"linit") == 0)
    {
        SolverPtr = &CMGSystem::LinIt;
    }
    else if(strcmp(solver,"gmres") == 0)
    {
        SolverPtr = &CMGSystem::GMRES;
    }
    else
    {
        ostrstream ostr;
        ostr << __FILE__ << __LINE__ <<  "solver = bicgstab" << endl;
        CMGWarning(ostr);
    }

    extra = extraptr;
    
    CMGMultiGrid *mg0 = CreateMultiGrid();
    if(mg0 == NULL) return 1;
    if (mg0->Init(*this)) return 1;
    if (mg0->Construct()) return 1;

    return 0;
}


int CMGSystem::Solve(double *rhs, double *defect, double *unknown)
{
    double *d, *u;
    int i;
    int status;

    if (rhs == NULL) return 1;
    vector[CMGRHS] = rhs;
    
    CMGMarkHeap(CMG_FROM_TOP);

    d = defect;
    if((d == NULL) || (d == rhs))
    {
        d = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_TOP);
        if (d == NULL) return 1;
    }
    vector[CMGDEFECT] = d;
 
    u = unknown;
    if((u == NULL) || (u == rhs) || (u == defect) )
    {
        u = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_TOP);
        if (u == NULL) return 1;
    }
    vector[CMGUNKNOWN] = u;
    
    CMGGrid *g = mg[0]->GetGrid(0);
    g->SetVector(CMGUNKNOWN,u);
    g->SetVector(CMGDEFECT,d);
    g->SetVector(CMGRHS,rhs);
    if(g->OrderVector(CMGRHS,g->GetMap())) return 1;

   status = (this->*SolverPtr)();
   if(status != NULL) return status;
       
   if(g->ReorderVector(CMGRHS,g->GetMap())) return 1;
   if(g->ReorderVector(CMGDEFECT,g->GetMap())) return 1;
   if(g->ReorderVector(CMGUNKNOWN,g->GetMap())) return 1;

  if(defect != d) 
   {
       if (defect != NULL)
       {
           for(i = 0; i < n; i++) defect[i] = d[i];
       }
       d = NULL;
       vector[CMGDEFECT] = NULL;
       g->SetVector(CMGDEFECT,(double *)NULL);
   }
   if(unknown != u) 
   { 
       if (unknown != NULL)
       {
           for(i = 0; i < n; i++) unknown[i] = u[i];
       }
       u = NULL;
       vector[CMGUNKNOWN] = NULL;
       g->SetVector(CMGUNKNOWN,(double *)NULL);
   }
    
    CMGReleaseHeap(CMG_FROM_TOP);

    return status;
}

int CMGSystem::Deconstruct()
{
    CMGMultiGrid *mg0=mg[0];

    if(mg0 != NULL) if(mg0->Deconstruct()) return 1;

    // repair matrix
    if(matrix->ReorderColumns(colmap)) return 1;
    matrix->RemodifyIndex(cmgfirsti);
    CMGReleaseHeap(CMG_FROM_BOTTOM);
    CMGReleaseHeap(CMG_FROM_TOP);

    return 0;
}

int CMGSystem::DeconstructSimple()
{
    CMGMultiGrid *mg0=mg[0];

    if(mg0 != NULL) if(mg0->Reorder()) return 1;

    // repair matrix
    if(matrix->ReorderColumns(colmap)) return 1;
    matrix->RemodifyIndex(colmap,cmgfirsti);
    CMGReleaseHeap(CMG_FROM_BOTTOM);
    return 0;
}

CMGSystem::CMGSystem()
{
    int i;

    nmg = 0;
    for(i = 0; i < CMGMULTIGRIDS; i++) mg[i] = NULL; 
    matrix = NULL;
    for(i = 0; i < CMGMAXVECTORS; i++) vector[i] = NULL;
    extra = NULL;
    SolverPtr = &CMGSystem::LinIt;
}


CMGMultiGrid *CMGSystem::CreateMultiGrid()
{
    CMGMultiGrid *newmg;
   
    newmg = (CMGMultiGrid *) CMGGetMem(sizeof(CMGMultiGrid),CMG_FROM_TOP);
    if(newmg == NULL) return(NULL);
    mg[nmg] = newmg;
    nmg++;
    
    return newmg;
}

#ifdef WRITE_VECTOR_TEST
static void CMGWriteVector(int n, double *vec)
{
    ofstream vfile("vec.dat",ios::out);
    if (!vfile) return;
    for(int i = 0; i < n; i++)
    {
        vfile << i << "\t" << vec[i] << endl;
    }

    return;
}
#endif


int CMGSystem::LinIt()
{
    double rlimit,alimit,divlimit,reduction,limit,defectnorm,startdefect,oldnorm;
    int maxit,i;
    ostrstream ostr;

    maxit = cmgparaptr->Getmaxit();
    rlimit = cmgparaptr->Getrlimit();
    alimit = cmgparaptr->Getalimit();
    divlimit =cmgparaptr->Getdivlimit();
    reduction = cmgparaptr->Getreduction();

    CMGSetVector(n,vector[CMGUNKNOWN],0.0);
    CMGCopyVector(n,vector[CMGDEFECT],vector[CMGRHS]);
    startdefect = CMGNorm(n,vector[CMGDEFECT]);
    limit = rlimit*startdefect;
    defectnorm = startdefect;
    ostr << 0 << "\t" <<  startdefect << endl;
    CMGWrite(ostr);

    // CMGMultiGrid * mg0 = CreateMultiGrid();
    // if(mg0 == NULL) return 1;
    // if (mg0->Init(*this)) return 1;
    // if (mg0->Construct()) return 1;
    CMGMultiGrid *mg0 = mg[0];
   
    for(i = 0; i < maxit; i++)
    {
        if((defectnorm < alimit) || (defectnorm < limit)) break;
        if(mg0->Step(0)) { mg0->Deconstruct(); return 1;}
        oldnorm = defectnorm;
        defectnorm = CMGNorm(n,vector[CMGDEFECT]);
        ostr << i+1 << "\t" << defectnorm << "\t" << defectnorm/oldnorm << endl;    
        CMGWrite(ostr);
    
        if((defectnorm < alimit) || (defectnorm < limit)) break;
        if (defectnorm/oldnorm > divlimit)
        {
            ostr << "diverged" << endl;
            CMGWrite(ostr);
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

int CMGSystem::AdTVSolve()
{
    double rlimit,alimit,divlimit,reduction,limit,defectnorm,startdefect,oldnorm;
    int maxit,i,j,init;
    ostrstream ostr;

    maxit = cmgparaptr->Getmaxit();
    rlimit = cmgparaptr->Getrlimit();
    alimit = cmgparaptr->Getalimit();
    divlimit =cmgparaptr->Getdivlimit();
    reduction = cmgparaptr->Getreduction();
    init = 8;

    CMGSetVector(n,vector[CMGUNKNOWN],0.0);
    CMGCopyVector(n,vector[CMGDEFECT],vector[CMGRHS]);
    startdefect = CMGNorm(n,vector[CMGDEFECT]);
    limit = rlimit*startdefect;
    defectnorm = startdefect;
    ostr << 0 << "\t" <<  startdefect << endl;
    CMGWrite(ostr);

    // CMGMultiGrid * mg0 = CreateMultiGrid();
    // if(mg0 == NULL) return 1;
    // if (mg0->Init(*this)) return 1;
    // if (mg0->Construct()) return 1;
    CMGMultiGrid *mg0 = mg[0];
   
    for(i = 0; i < maxit;)
    {
        if((defectnorm < alimit) || (defectnorm < limit)) break;
        for(j = 0; (j < init) && (i < maxit); j++, i++)
        {
            CMGCopyVector(n,vector[CMGTVA],vector[CMGUNKNOWN]);
            CMGCopyVector(n,vector[CMGTVB],vector[CMGUNKNOWN]);
            if(mg0->Step(0)) { mg0->Deconstruct(); return 1;}
            CMGSubVector(n,vector[CMGTVA],vector[CMGUNKNOWN]);
            CMGSubVector(n,vector[CMGTVB],vector[CMGUNKNOWN]);
            oldnorm = defectnorm;
            defectnorm = CMGNorm(n,vector[CMGDEFECT]);
            ostr << i+1 << "\t" << defectnorm << "\t" << defectnorm/oldnorm << endl;    
            CMGWrite(ostr);
            
            if((defectnorm < alimit) || (defectnorm < limit)) break;
            if (defectnorm/oldnorm > divlimit)
            {
                ostr << "diverged" << endl;
                CMGWrite(ostr);
                break;
            }


        }
        if (defectnorm/oldnorm > divlimit)
        {
            ostr << "diverged" << endl;
            CMGWrite(ostr);
            break;
        }

        if((defectnorm < alimit) || (defectnorm < limit)) break;
        if(i == maxit) break;

        if(mg0->Deconstruct()) return 1;
        if (mg0->Construct()) return 1;
            
    }
    
    if (defectnorm < startdefect*reduction)
    {
        // Yeah, problem solved !!!

        return 0;
    }
    
    return 1;
}




// BiCGStab

int CMGSystem::BiCGStab()
{
    double rlimit,alimit,divlimit,reduction,limit,defectnorm,startdefect,oldnorm;
    double *vec[6], rho, oldrho, alpha, beta, omega, nenner;
    int maxit,i;
    ostrstream ostr;


    const int CMGX = 0;
    const int CMGR = 1;
    const int CMGV = 2;
    const int CMGP = 3;
    const int CMGT = 4;
    const int CMGR0 = 5;

    CMGMarkHeap(CMG_FROM_BOTTOM);
    vec[CMGX] = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_BOTTOM);
    if(vec[CMGX] == NULL) return 1;
    vec[CMGR] = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_BOTTOM);
    if(vec[CMGR] == NULL) return 1;
    vec[CMGV] = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_BOTTOM);
    if(vec[CMGV] == NULL) return 1;
    vec[CMGP] = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_BOTTOM);
    if(vec[CMGP] == NULL) return 1;
    vec[CMGT] = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_BOTTOM);
    if(vec[CMGT] == NULL) return 1;
    vec[CMGR0] = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_BOTTOM);
    if(vec[CMGR0] == NULL) return 1;

    maxit = cmgparaptr->Getmaxit();
    rlimit = cmgparaptr->Getrlimit();
    alimit = cmgparaptr->Getalimit();
    divlimit =cmgparaptr->Getdivlimit();
    reduction = cmgparaptr->Getreduction();


    // construct preconditioner
    // CMGMultiGrid * mg0 = CreateMultiGrid();
    // if(mg0 == NULL) return 1;
    // if (mg0->Init(*this)) return 1;
    // if (mg0->Construct()) return 1;
    CMGMultiGrid *mg0 = mg[0];
   
    CMGSetVector(n,vec[CMGX],0.0);
    CMGCopyVector(n,vec[CMGR0],vector[CMGRHS]);

    
    defectnorm = CMGNorm(n,vec[CMGR0]);
    startdefect = defectnorm;
    limit = rlimit*startdefect;
    ostr << 0 << "\t" <<  startdefect << endl;
    CMGWrite(ostr);

    CMGCopyVector(n,vec[CMGR],vec[CMGR0]);
    CMGCopyVector(n,vec[CMGP],vec[CMGR]);
    rho = CMGSum(n,vec[CMGR]);
    if (Abs(rho) < 1e-10*alimit) 
    {
       ostr << __FILE__ << ", line " << __LINE__ << ": rho too small" << endl;
       CMGWarning(ostr);
    }

    for(i = 0; i < maxit; i++)
    {

        CMGCopyVector(n,vector[CMGRHS],vec[CMGP]);
        CMGSetVector(n,vector[CMGUNKNOWN],0.0);
        CMGCopyVector(n,vector[CMGDEFECT],vector[CMGRHS]);
        if(mg0->Step(0)) { mg0->Deconstruct(); return 1;}
        CMGSetSubVector(n,vec[CMGV],vec[CMGP],vector[CMGDEFECT]);

        nenner = CMGSum(n,vec[CMGV]);
        if (Abs(nenner) < 1e-15*Abs(rho)) 
        {
            ostr << __FILE__ << ", line " << __LINE__ << ": nenner too small" << endl;
            CMGWarning(ostr);
        }
    
        alpha = rho/nenner;
        CMGAddVector(n,vec[CMGX],vector[CMGUNKNOWN],alpha);
        CMGAddVector(n,vec[CMGR],vec[CMGV],-alpha);

        oldnorm = defectnorm;
        defectnorm = CMGNorm(n,vec[CMGR]);
        ostr << i+0.5 << "\t" << defectnorm << "\t" << defectnorm/oldnorm;
        ostr << "\t" << alpha  << endl; 
        CMGWrite(ostr);
        if((defectnorm < alimit) || (defectnorm < limit)) break;
        if (defectnorm/oldnorm > divlimit)
        {
            ostr << "diverged" << endl;
            CMGWrite(ostr);
           break;
        }
        
        CMGCopyVector(n,vector[CMGRHS],vec[CMGR]);
        CMGSetVector(n,vector[CMGUNKNOWN],0.0);
        CMGCopyVector(n,vector[CMGDEFECT],vector[CMGRHS]);
        if(mg0->Step(0)) { mg0->Deconstruct(); return 1;}
        CMGSetSubVector(n,vec[CMGT],vec[CMGR],vector[CMGDEFECT]);

        omega = CMGScalProd(n,vec[CMGT],vec[CMGR])/CMGScalProd(n,vec[CMGT],vec[CMGT]);
        if (Abs(omega) < 1e-15) 
        {
            ostr << __FILE__ << ", line " << __LINE__ << ": omega too small" << endl;
            CMGWarning(ostr);
        }

        CMGAddVector(n,vec[CMGX],vector[CMGUNKNOWN],omega);
        CMGAddVector(n,vec[CMGR],vec[CMGT],-omega);
        CMGAddVector(n,vec[CMGP],vec[CMGV],-omega);

        oldnorm = defectnorm;
        defectnorm = CMGNorm(n,vec[CMGR]);
        ostr << i+1 << "\t" << defectnorm << "\t" << defectnorm/oldnorm;
        ostr << "\t" << omega << endl;   
        CMGWrite(ostr);
        if((defectnorm < alimit) || (defectnorm < limit)) break;
        if (defectnorm/oldnorm > divlimit)
        {
            ostr << "diverged" << endl;
            CMGWrite(ostr);            
            break;
        }

        oldrho = rho;
        rho = CMGSum(n,vec[CMGR]);
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
    CMGCopyVector(n,vector[CMGUNKNOWN],vec[CMGX]);
    CMGCopyVector(n,vector[CMGRHS],vec[CMGR0]);

    CMGReleaseHeap(CMG_FROM_BOTTOM);
    
    if (defectnorm < startdefect*reduction)
    {
        // Yeah, problem solved !!!

        return 0;
    }

    
    return 1;
}

    // GMRES


void CMGMatVecMult(int n, int nr, double* M, double* h)
{
    double sum, *g;
    int i,j;

    
    CMGMarkHeap(CMG_FROM_BOTTOM);
    g = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_BOTTOM);

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
            
    CMGReleaseHeap(CMG_FROM_BOTTOM);

    return;
}

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
        
int CMGSystem::Arnoldi(CMGMultiGrid *mg0, double **vec, double *H, double *G, double *Q, double *P, double &q0, int con)
{
     double q, gamma, oldgamma, h[cmgnv], c, s, nenner;
     int i,j,nv;
     ostrstream ostr;
           
     nv = cmgparaptr->Getnv();  
     for(i = 0; i < (nv+1)*(nv+1); i++) Q[i] = 0.0;

     q = CMGNorm(n,vector[CMGDEFECT]);
     CMGCopyScaledVector(n,vec[0],vector[CMGDEFECT],1.0/q); 
     q0 = q;
     Q[0] = 1.0;
     oldgamma = q0;

     // Arnoldi process
     for(i = 0; i < nv; i++)
     {
         // def = M^{-1}A v_{i-1}
         CMGCopyScaledVector(n,vec[i],vector[CMGDEFECT],1.0/q); 
         CMGCopyVector(n,vector[CMGRHS],vec[i]);
         CMGSetVector(n,vector[CMGUNKNOWN],0.0);
         CMGCopyVector(n,vector[CMGDEFECT],vector[CMGRHS]);
         if(con)
         {
             if(mg0->Step(0)) { mg0->Deconstruct(); return 1;}
         }
         else
         {
             mg0->SGSStep(0);
         }
         CMGSetSubVector(n,vector[CMGDEFECT],vector[CMGRHS],vector[CMGDEFECT]);
 
         // modified Gram-Schmitt
         for(j = 0; j <= i; j++)
         {
             h[j] = CMGScalProd(n,vector[CMGDEFECT],vec[j]);
             CMGAddVector(n,vector[CMGDEFECT],vec[j],-h[j]);
         }

         q = CMGNorm(n,vector[CMGDEFECT]);
         if(q < 1e-10) 
         {
             ostr << __FILE__ << ", line " << __LINE__ << ": GMRES breakdownl" << endl;
             CMGWarning(ostr);
         }
         // => def = v_i*q

         // update H
         CMGMatVecMult(i+1,nv+1,Q,h);
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
         CMGWrite(ostr);
         oldgamma = gamma;

     }

     return 0;
}

int CMGSystem::UpdateSolution(CMGMultiGrid *mg0, double **vec, double *H, double *Q, double &q0, int con)
{
    double x[cmgnv];
    int j, nv;
    
    nv = cmgparaptr->Getnv();

    ComputeX(nv,nv,x,H,Q,q0);
    CMGSetVector(n,vector[CMGRHS],0.0);
    for(j = 0; j < nv; j++)
    {
        CMGAddVector(n,vector[CMGRHS],vec[j],x[j]);
    }
            

    CMGSetVector(n,vector[CMGUNKNOWN],0.0);
    CMGCopyVector(n,vector[CMGDEFECT],vector[CMGRHS]);
    
    
    if(con)
    {
        if(mg0->Step(0)) { mg0->Deconstruct(); return 1;}
    }
    else
    {
        mg0->SGSStep(0);
    }


    return 0;
}

int CMGSystem::ComputeEigenVector(CMGMultiGrid *mg0, double **vec, double *G, double *P, int con)
{
    double e[cmgnv], ew, norm;
    int j, nv;

    nv = cmgparaptr->Getnv();

    // compute eigenvector to smallest eigenvalue
    ew = ComputeEV(nv,nv,G,P,e);
    ostrstream ostr; ostr << "ew = " << ew << endl;
    CMGWrite(ostr);
    if(Abs(ew) < 0.1)
    {
        CMGSetVector(n,vector[CMGRHS],0.0);
        CMGSetVector(n,vector[CMGUNKNOWN],0.0);
        for(j = 0; j < nv; j++)
        {
            CMGAddVector(n,vector[CMGRHS],vec[j],e[j]);
        }


        CMGSetVector(n,vector[CMGUNKNOWN],0.0);
        CMGCopyVector(n,vector[CMGDEFECT],vector[CMGRHS]);
        if(con)
        {
            if(mg0->Step(0)) { mg0->Deconstruct(); return 1;}
        }
        else
        {
            mg0->SGSStep(0);
        }
        norm = CMGNorm(n,vector[CMGUNKNOWN]);
        CMGCopyScaledVector(n,vector[CMGTVA],vector[CMGUNKNOWN],1.0/norm);
        CMGCopyScaledVector(n,vector[CMGTVB],vector[CMGUNKNOWN],1.0/norm);

        return -1;
    }

    return 0;
}
            
int CMGSystem::GMRES()
{
    double rlimit,alimit,reduction,limit,defectnorm,startdefect;
    int maxit;
    double H[cmgnv*cmgnv], G[cmgnv*cmgnv], Q[(cmgnv+1)*(cmgnv+1)], P[(cmgnv+1)*(cmgnv+1)], *vec[cmgnv], *sol, *rhs, *def, q0;
    int i, k, newtv, con, nv;
    ostrstream ostr;

    maxit = cmgparaptr->Getmaxit();
    rlimit = cmgparaptr->Getrlimit();
    alimit = cmgparaptr->Getalimit();
    reduction = cmgparaptr->Getreduction();
    nv = cmgparaptr->Getnv();

    CMGMarkHeap(CMG_FROM_BOTTOM);

    sol = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_BOTTOM);
    if(sol == NULL) return 1;
    rhs = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_BOTTOM);
    if(rhs == NULL) return 1;
    
    // test
    def = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_BOTTOM);
    if(def == NULL) return 1;

    for(i = 0; i < nv; i++)
    {
        vec[i] = (double *) CMGGetMem(n*sizeof(double),CMG_FROM_BOTTOM);
        if(vec[i] == NULL) return 1;
    }
    

    CMGSetVector(n,vector[CMGUNKNOWN],0.0);
    CMGSetVector(n,sol,0.0);
    CMGCopyVector(n,vector[CMGDEFECT],vector[CMGRHS]);

    startdefect = defectnorm = CMGNorm(n,vector[CMGDEFECT]);
    limit = rlimit*defectnorm;
    ostr << "start defect: " << defectnorm  << endl;
    CMGWrite(ostr);


    // CMGMultiGrid * mg0 = CreateMultiGrid();
    // if(mg0 == NULL) return 1;
    // if (mg0->Init(*this)) return 1;
    // if (mg0->Construct()) return 1;
    CMGMultiGrid *mg0 = mg[0];

    newtv = 0;
    con = 1;

    for(k = 0; k*nv < maxit; k++)
    {

        // construct preconditioner
        if(newtv)
        {
            if(con) 
            {
                if(mg0->Deconstruct()) return 1;
            }
            if (mg0->Construct()) return 1;
            con = 1;
       }

 
        CMGCopyVector(n,rhs,vector[CMGRHS]);
        CMGCopyVector(n,sol,vector[CMGUNKNOWN]);
        //test
        CMGCopyVector(n,def,vector[CMGDEFECT]);

        if(Arnoldi(mg0,vec,H,G,Q,P,q0,con)) return 1;

        newtv = 0;
        // status = ComputeEigenVector(mg0,vec,G,P,con);
        // if(status > 0) return 1; // error
        // if(status < 0) 
        // {
        //     newtv = 1; // restart with new TV
        // }


        if(UpdateSolution(mg0,vec,H,Q,q0,con)) return 1;

        CMGAddVector(n,vector[CMGUNKNOWN],sol);

        CMGCopyVector(n,vector[CMGRHS],rhs);
        mg0->GetGrid(0)->Defect();
        defectnorm = CMGNorm(n,vector[CMGDEFECT]);
        ostr << defectnorm << endl;
        CMGWrite(ostr);
        if( (defectnorm < alimit) || (defectnorm < limit)) break;
    }

    if(!con) 
    {
        ostr  << __FILE__ << __LINE__ << endl;
        CMGError(ostr);
        return 1;
    }

    CMGReleaseHeap(CMG_FROM_BOTTOM);

    if (defectnorm < startdefect*reduction)
    {
        // Yeah, problem solved !!!

        return 0;
    }

    return 1;
}

