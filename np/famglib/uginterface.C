/****************************************************************************/
/*																			*/
/* File:      uginterface.C													*/
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
#include "matrix.h"
#include "system.h"
#include "heap.h"

/* RCS_ID
$Header$
*/

    // actually these should be a system member functions. Maybe it is more
    // flexibel this way.

static CMGSystem *cmgsystemptr;

static void ReadParameter(CMGParameter *parameter, CMGParameter *in_parameter)
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


static int CMGConstructParameter(CMGParameter *in_parameter)
{
    CMGParameter *parameter = new CMGParameter;
    if(parameter == NULL) return 1;
    ReadParameter(parameter, in_parameter); 
    CMGSetParameter(parameter);

    return 0;
}

static void CMGDeconstructParameter()
{
    CMGParameter *parameter = CMGGetParameter();
    delete parameter;
}

static int CMGConstruct(double *matrix, int *index, int *start, int n, int nl, double *tvA, double *tvB, void **extra)
{
 

    CMGHeap *heap = new CMGHeap (CMGGetParameter()->Getheap());
    if (heap == NULL) return 1;
    CMGSetHeap(heap);

    cmgsystemptr = (CMGSystem *) CMGGetMem(sizeof(CMGSystem),CMG_FROM_TOP);
    if(cmgsystemptr == NULL) return 1;
    
    cmgsystemptr->Init();
    if(cmgsystemptr->Construct(matrix,index,start,n,nl,tvA,tvB,extra)) return 1;

    return 0;
}
static int CMGConstructSimple(double *matrix, int *index, int *start, int n, int nl, void **extra)
{
    if(cmgsystemptr->ConstructSimple(matrix,index,start,n,nl,extra)) return 1;

    return 0;
}

static int CMGSolve(double *rhs, double *defect, double *unknown)
{
    return cmgsystemptr->Solve(rhs,defect,unknown);
}

static void CMGDeconstruct()
{
    cmgsystemptr->Deconstruct();
    
    CMGHeap *heap = CMGGetHeap();
    delete heap;
}

static void CMGDeconstructSimple()
{
    cmgsystemptr->DeconstructSimple();

    // remove heap, parameter
}

int CMGSolveSystem(CMG_Interface *interface, CMG_Parameter *in_parameter)
{
    int status;


    CMGConstructParameter((CMGParameter*)in_parameter);

    CMGConstruct(interface->entry,interface->index,interface->start,interface->n,interface->nl,interface->vector[CMG_TVA],interface->vector[CMG_TVB],interface->extra);

    status = CMGSolve(interface->vector[CMG_RHS],interface->vector[CMG_DEFECT],interface->vector[CMG_UNKNOWN]);

    CMGDeconstructSimple();

    CMGConstructSimple(interface->entry,interface->index,interface->start,interface->n,interface->nl,interface->extra);

    status = CMGSolve(interface->vector[CMG_RHS],interface->vector[CMG_DEFECT],interface->vector[CMG_UNKNOWN]);

    CMGDeconstruct();
    CMGConstruct(interface->entry,interface->index,interface->start,interface->n,interface->nl,interface->vector[CMG_TVA],interface->vector[CMG_TVB],interface->extra);

    status = CMGSolve(interface->vector[CMG_RHS],interface->vector[CMG_DEFECT],interface->vector[CMG_UNKNOWN]);

    status = CMGSolve(interface->vector[CMG_RHS],interface->vector[CMG_DEFECT],interface->vector[CMG_UNKNOWN]);

    CMGDeconstructSimple();

    CMGConstructSimple(interface->entry,interface->index,interface->start,interface->n,interface->nl,interface->extra);

    status = CMGSolve(interface->vector[CMG_RHS],interface->vector[CMG_DEFECT],interface->vector[CMG_UNKNOWN]);

    status = CMGSolve(interface->vector[CMG_RHS],interface->vector[CMG_DEFECT],interface->vector[CMG_UNKNOWN]);

     status = CMGSolve(interface->vector[CMG_RHS],interface->vector[CMG_DEFECT],interface->vector[CMG_UNKNOWN]);

   CMGDeconstruct();

   CMGDeconstructParameter();
 
    return status;
}

void **CMG_GetExtraPtr(int level)
{
    void **ptr;

    ptr = cmgsystemptr->GetMultiGrid(0)->GetGrid(level)->GetNode();

    return ptr;
}

int CMG_GetN(int level)
{
    int n;

    n = cmgsystemptr->GetMultiGrid(0)->GetGrid(level)->GetN();
    return n;
}

int CMG_GetNF(int level)
{
    int nf;

    nf = cmgsystemptr->GetMultiGrid(0)->GetGrid(level)->GetNF();
    return nf;
}

CMG_MatrixPtr  *CMG_GetMatrixPtr(int level,int i)
{
    CMGMatrix *matrix;
    CMGMatrixPtr mat;
    
    matrix = cmgsystemptr->GetMultiGrid(0)->GetGrid(level)->GetMatrix();
    mat = matrix->GetStart(i);
    return ((CMG_MatrixPtr *) &mat);
}

CMG_TransferEntry  *CMG_GetTransferEntry(int level,int i)
{
    CMGTransfer *transfer;
    CMGTransferEntry *trans;
    
    transfer = cmgsystemptr->GetMultiGrid(0)->GetGrid(level)->GetTransfer();
    if(transfer == NULL) return NULL;
    trans = transfer->GetRow(i);
    return ((CMG_TransferEntry *) trans);
}

double * CMG_GetVector(int level, int i)
{
    double *vector;
    
    vector = cmgsystemptr->GetMultiGrid(0)->GetGrid(level)->GetVector(i);
    return vector;
}


int CMG_GetMaxLevel()
{
    int n;

    n = cmgsystemptr->GetMultiGrid(0)->GetN();
    return n-1;
}
