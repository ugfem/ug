/****************************************************************************/
/*																			*/
/* File:      famg_interface.C												*/
/*																			*/
/* Purpose:   famg interface												*/
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

#include "famg_system.h"
#include "famg_heap.h"
#include "famg_interface.h"

/* RCS_ID
$Header$
*/

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

int FAMGConstructParameter(FAMGParameter *in_parameter)
{
    FAMGParameter *parameter = new FAMGParameter;
    if(parameter == NULL) RETURN(1);
    ReadParameter(parameter, in_parameter); 
    FAMGSetParameter(parameter);

    return 0;
}

int FAMGDeconstructParameter()
{
    FAMGParameter *parameter = FAMGGetParameter();
    delete parameter;

    return 0;
}

int FAMGConstruct(FAMG_Interface* famg_interface)
{
 

    FAMGHeap *heap = new FAMGHeap (FAMGGetParameter()->Getheap());
    if (heap == NULL) RETURN(1);
    FAMGSetHeap(heap);

    famgsystemptr = (FAMGSystem *) FAMGGetMem(sizeof(FAMGSystem),FAMG_FROM_TOP);
    if(famgsystemptr == NULL) RETURN(1);
    
    famgsystemptr->Init();
    if(famgsystemptr->Construct(famg_interface)) RETURN(1);

    return 0;
}

int FAMGPrepare(FAMG_Interface* famg_interface)
{
    if(famgsystemptr->ConstructSimple(famg_interface)) RETURN(1);

    return 0;
}

int FAMGSolve(double *rhs, double *defect, double *unknown)
{
    return famgsystemptr->Solve(rhs,defect,unknown);
}

int FAMGDeconstruct()
{
    if(famgsystemptr->Deconstruct()) RETURN(1);
    
    FAMGHeap *heap = FAMGGetHeap();
    delete heap;
    
    return 0;
}

int FAMGRepair()
{
    if(famgsystemptr->DeconstructSimple()) RETURN(1);
    
    return 0;
}

FAMGSystem *FAMG_GetSystem()
{
	return famgsystemptr;
}

