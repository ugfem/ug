/****************************************************************************/
/*																			*/
/* File:      famg_uginterface.C											*/
/*																			*/
/* Purpose:   famg - ug interface											*/
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

extern "C"
{
#include "famginterface.h"
}
#include "famg_matrix.h"
#include "famg_system.h"
#include "famg_heap.h"

/* RCS_ID
$Header$
*/

    // actually these should be a system member functions. Maybe it is more
    // flexibel this way.

static FAMGSystem *famgsystemptr;

static void ReadParameter(FAMGParameter *parameter, FAMGParameter *in_parameter)
{
    parameter->Setheap(in_parameter->Getheap());
    parameter->Setgamma(in_parameter->Getgamma());
    parameter->Setn1(in_parameter->Getn1());
    parameter->Setn2(in_parameter->Getn2());
    parameter->Setilut(in_parameter->Getilut());
    parameter->Setcgilut(in_parameter->Getcgilut());
    parameter->Setcgnodes(in_parameter->Getcgnodes());
    parameter->Setconloops(in_parameter->Getconloops());
    parameter->Setmincoarse(in_parameter->Getmincoarse());
    parameter->Settype(in_parameter->Gettype());
    parameter->Setstv(in_parameter->Getstv());
    parameter->Settol(in_parameter->Gettol());
    parameter->Setsigma(in_parameter->Getsigma());
    parameter->Setomegar(in_parameter->Getomegar());
    parameter->Setomegal(in_parameter->Getomegal());
    parameter->Seterror1(in_parameter->Geterror1());
    parameter->Seterror2(in_parameter->Geterror2());
    parameter->Setmaxit(in_parameter->Getmaxit());
    parameter->Setalimit(in_parameter->Getalimit());
    parameter->Setrlimit(in_parameter->Getrlimit());
    parameter->Setdivlimit(in_parameter->Getdivlimit());
    parameter->Setreduction(in_parameter->Getreduction());
    parameter->Setsolver(in_parameter->Getsolver());
    parameter->Setpresmoother(in_parameter->Getpresmoother());
    parameter->Setpostsmoother(in_parameter->Getpostsmoother());
    parameter->Setcgsmoother(in_parameter->Getcgsmoother());
}    


static int FAMGConstructParameter(FAMGParameter *in_parameter)
{
    FAMGParameter *parameter = new FAMGParameter;
    if(parameter == NULL) return 1;
    ReadParameter(parameter, in_parameter); 
    FAMGSetParameter(parameter);

    return 0;
}

static void FAMGDeconstructParameter()
{
    FAMGParameter *parameter = FAMGGetParameter();
    delete parameter;
}

static int FAMGConstruct(double *matrix, int *index, int *start, int n, int nl, double *tvA, double *tvB, void **extra)
{
 

    FAMGHeap *heap = new FAMGHeap (FAMGGetParameter()->Getheap());
    if (heap == NULL) return 1;
    FAMGSetHeap(heap);

    famgsystemptr = (FAMGSystem *) FAMGGetMem(sizeof(FAMGSystem),FAMG_FROM_TOP);
    if(famgsystemptr == NULL) return 1;
    
    famgsystemptr->Init();
    if(famgsystemptr->Construct(matrix,index,start,n,nl,tvA,tvB,extra)) return 1;

    return 0;
}
static int FAMGConstructSimple(double *matrix, int *index, int *start, int n, int nl, void **extra)
{
    if(famgsystemptr->ConstructSimple(matrix,index,start,n,nl,extra)) return 1;

    return 0;
}

static int FAMGSolve(double *rhs, double *defect, double *unknown)
{
    return famgsystemptr->Solve(rhs,defect,unknown);
}

static void FAMGDeconstruct()
{
    famgsystemptr->Deconstruct();
    
    FAMGHeap *heap = FAMGGetHeap();
    delete heap;
}

static void FAMGDeconstructSimple()
{
    famgsystemptr->DeconstructSimple();

    // remove heap, parameter
}

int FAMGSolveSystem(FAMG_Interface *interface, FAMG_Parameter *in_parameter)
{
    int status;


    FAMGConstructParameter((FAMGParameter*)in_parameter);

    FAMGConstruct(interface->entry,interface->index,interface->start,interface->n,interface->nl,interface->vector[FAMG_TVA],interface->vector[FAMG_TVB],interface->extra);

    status = FAMGSolve(interface->vector[FAMG_RHS],interface->vector[FAMG_DEFECT],interface->vector[FAMG_UNKNOWN]);

	FAMGDeconstructSimple();

	FAMGDeconstructParameter();
 
    return status;
}

void **FAMG_GetExtraPtr(int level)
{
    void **ptr;

    ptr = famgsystemptr->GetMultiGrid(0)->GetGrid(level)->GetNode();

    return ptr;
}

int FAMG_GetN(int level)
{
    int n;

    n = famgsystemptr->GetMultiGrid(0)->GetGrid(level)->GetN();
    return n;
}

int FAMG_GetNF(int level)
{
    int nf;

    nf = famgsystemptr->GetMultiGrid(0)->GetGrid(level)->GetNF();
    return nf;
}

FAMG_MatrixPtr  *FAMG_GetMatrixPtr(int level,int i)
{
    FAMGMatrix *matrix;
    FAMGMatrixPtr mat;
    
    matrix = famgsystemptr->GetMultiGrid(0)->GetGrid(level)->GetMatrix();
    mat = matrix->GetStart(i);
    return ((FAMG_MatrixPtr *) &mat);
}

FAMG_TransferEntry  *FAMG_GetTransferEntry(int level,int i)
{
    FAMGTransfer *transfer;
    FAMGTransferEntry *trans;
    
    transfer = famgsystemptr->GetMultiGrid(0)->GetGrid(level)->GetTransfer();
    if(transfer == NULL) return NULL;
    trans = transfer->GetRow(i);
    return ((FAMG_TransferEntry *) trans);
}

double * FAMG_GetVector(int level, int i)
{
    double *vector;
    
    vector = famgsystemptr->GetMultiGrid(0)->GetGrid(level)->GetVector(i);
    return vector;
}


int FAMG_GetMaxLevel()
{
    int n;

    n = famgsystemptr->GetMultiGrid(0)->GetN();
    return n-1;
}
