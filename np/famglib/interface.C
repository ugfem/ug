/****************************************************************************/
/*																			*/
/* File:      interface.C													*/
/*																			*/
/* Purpose:   cmg interface													*/
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

#include "system.h"
#include "heap.h"
#include "interface.h"

/* RCS_ID
$Header$
*/

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

int CMGConstructParameter(CMGParameter *in_parameter)
{
    CMGParameter *parameter = new CMGParameter;
    if(parameter == NULL) return 1;
    ReadParameter(parameter, in_parameter); 
    CMGSetParameter(parameter);

    return 0;
}

int CMGDeconstructParameter()
{
    CMGParameter *parameter = CMGGetParameter();
    delete parameter;

    return 0;
}

int CMGConstruct(double *matrix, int *index, int *start, int n, int nl, double *tvA, double *tvB, void **extra)
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
int CMGPrepare(double *matrix, int *index, int *start, int n, int nl, void **extra)
{
    if(cmgsystemptr->ConstructSimple(matrix,index,start,n,nl,extra)) return 1;

    return 0;
}

int CMGSolve(double *rhs, double *defect, double *unknown)
{
    return cmgsystemptr->Solve(rhs,defect,unknown);
}

int CMGDeconstruct()
{
    if(cmgsystemptr->Deconstruct()) return 1;
    
    CMGHeap *heap = CMGGetHeap();
    delete heap;
    
    return 0;
}

int CMGRepair()
{
    if(cmgsystemptr->DeconstructSimple()) return 1;
    
    return 0;
}
