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

#include <iostream.h>

#include "famg_uginterface.h"
#include "famg_algebra.h"
#include "famg_system.h"
#include "famg_heap.h"

/* RCS_ID
$Header$
*/

    // actually these should be a system member functions. Maybe it is more
    // flexibel this way.

static FAMGSystem *famgsystemptr;

#ifdef WEG
// TODO: weg!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/* nur fuer Tests */
#include "famg_ugalgebra.h"
extern "C"
{
#include "commands.h"
}

void test_algebra(void)
{
	GRID *grid = GRID_ON_LEVEL(GetCurrentMultigrid(),0);
	
	FAMGugVector sol ( grid,0);
	//FAMGVectorEntry ve(sol.firstEntry());
	FAMGVectorEntry ve;
	
	FAMGVectorIter soliter(sol);
	
	FAMGugMatrix M(0);
	FAMGMatrixEntry me;

	while( soliter(ve) )
	{
		cout << sol[ve] << "Matrixeintraege:";
		
		FAMGMatrixIter Miter(M,ve);
		while( Miter(me) )
		{
			cout << M[me] << " ";	
		}
		cout << endl;
	}
	
	printf( "test\n" );
	return;
}
#endif

static void ReadParameter(FAMGParameter *parameter, FAMGParameter_ug *in_parameter)
{

//test_algebra(); assert(0);

    parameter->Setheap(in_parameter->heap);
    parameter->Setgamma(in_parameter->gamma);
    parameter->Setn1(in_parameter->n1);
    parameter->Setn2(in_parameter->n2);
    parameter->Setilut(in_parameter->ilut);
    parameter->Setcgilut(in_parameter->cgilut);
    parameter->Setcgnodes(in_parameter->cgnodes);
    parameter->Setcglevels(in_parameter->cglevels);
    parameter->Setconloops(in_parameter->conloops);
    parameter->Setmincoarse(in_parameter->mincoarse);
    parameter->Settype(in_parameter->type);
    parameter->Setstv(in_parameter->stv);
    parameter->Settol(in_parameter->tol);
    parameter->Setsigma(in_parameter->sigma);
    parameter->Setomegar(in_parameter->omegar);
    parameter->Setomegal(in_parameter->omegal);
    parameter->Seterror1(in_parameter->error1);
    parameter->Seterror2(in_parameter->error2);
    parameter->Setmaxit(in_parameter->maxit);
    parameter->Setalimit(in_parameter->alimit);
    parameter->Setrlimit(in_parameter->rlimit);
    parameter->Setdivlimit(in_parameter->divlimit);
    parameter->Setreduction(in_parameter->reduction);
    parameter->Setsolver(in_parameter->solver);
    parameter->Setpresmoother(in_parameter->presmoother);
    parameter->Setpostsmoother(in_parameter->postsmoother);
    parameter->Setcgsmoother(in_parameter->cgsmoother);
}    


int FAMGConstructParameter(struct FAMGParameter_ug *in_parameter)
{
    FAMGParameter *parameter = new FAMGParameter;
    if(parameter == NULL) return 1;
    ReadParameter((FAMGParameter*)parameter, in_parameter); 
    FAMGSetParameter(parameter);

    return 0;
}

void FAMGDeconstructParameter()
{
    FAMGParameter *parameter = FAMGGetParameter();
    delete parameter;
}

int FAMGConstruct(FAMGGridVector *gridvector, FAMGMatrixAlg *matrix, FAMGMatrixAlg *Consmatrix, FAMGVector *vectors[FAMG_NVECTORS])
{
 

    FAMGHeap *heap = new FAMGHeap (FAMGGetParameter()->Getheap());
    if (heap == NULL) 
	{
		ostrstream ostr; ostr  << __FILE__ << ", line " << __LINE__ << ": can not allocate Heap" << endl;
		FAMGError(ostr);
		assert(0);
	}
    FAMGSetHeap(heap);

    famgsystemptr = (FAMGSystem *) FAMGGetMem(sizeof(FAMGSystem),FAMG_FROM_TOP);
    if(famgsystemptr == NULL)
	{
		ostrstream ostr; ostr  << __FILE__ << ", line " << __LINE__ << ": can not allocate FAMGSystem" << endl;
		FAMGError(ostr);
		assert(0);
	}
    
    famgsystemptr->Init();
    if(famgsystemptr->Construct(gridvector,matrix,Consmatrix,vectors))
	{
		ostrstream ostr; ostr  << __FILE__ << ", line " << __LINE__ << ": can not construct FAMGSystem" << endl;
		FAMGError(ostr);
		assert(0);
	}

    return 0;
}

int FAMGConstructSimple(FAMGMatrixAlg *matrix, FAMGVector *tvA, FAMGVector *tvB)
{
    if(famgsystemptr->ConstructSimple(matrix,tvA,tvB)) return 1;

    return 0;
}

int FAMGSolve(FAMGVector *rhs, FAMGVector *defect, FAMGVector *unknown)
{
    return famgsystemptr->Solve(rhs,defect,unknown);
}

void FAMGDeconstruct()
{
    famgsystemptr->Deconstruct();
    
    FAMGHeap *heap = FAMGGetHeap();
    delete heap;
}

void FAMGDeconstructSimple()
{
    famgsystemptr->DeconstructSimple();

    // remove heap, parameter
}

int FAMGSolveSystem(FAMG_Interface *interface)
{
    int status;


    //verschoben nach Preprocess FAMGConstructParameter((FAMGParameter*)in_parameter);

    //verschoben nach Preprocess FAMGConstruct(interface->entry,interface->index,interface->start,interface->n,interface->nl,interface->vector[FAMG_TVA],interface->vector[FAMG_TVB],interface->extra);

    status = FAMGSolve(interface->vector[FAMG_RHS],interface->vector[FAMG_DEFECT],interface->vector[FAMG_UNKNOWN]);

	FAMGDeconstructSimple();

	FAMGDeconstructParameter();
 
    return status;
}

int FAMG_RestrictDefect( int fine_level, VECDATA_DESC *to, VECDATA_DESC *from, VECDATA_DESC *smooth_sol, VECDATA_DESC *smooth_def )
// fine_level in the famg grid stack
{
	INT i;

	FAMGGrid &fg = *FAMG_GetSystem()->GetMultiGrid(0)->GetGrid(fine_level);
	FAMGGrid &cg = *FAMG_GetSystem()->GetMultiGrid(0)->GetGrid(fine_level+1);
		
	assert(from==smooth_def);
	//assert(VD_IS_SCALAR(from)); has not to be fulfilled
	//assert(VD_IS_SCALAR(to)); has not to be fulfilled
	//assert(VD_SCALCMP(from)==VD_SCALCMP(to)); has not to be fulfilled
	// now it is sure we can use for np->smooth_def, from and to the same component
	
	FAMGugVector fgsol(*(const FAMGugGridVector*)&fg.GetGridVector(),smooth_sol);
	FAMGugVector fgdef(*(const FAMGugGridVector*)&fg.GetGridVector(),from);
	FAMGugVector cgdef(*(const FAMGugGridVector*)&cg.GetGridVector(),to);
	
	// be sure not to use an out-dated pointer
	for (i=0; i<FAMGMAXVECTORS; i++ )
	{
		fg.SetVector(i, NULL);
		cg.SetVector(i, NULL);
	}
	
	fg.SetVector(FAMGUNKNOWN, &fgsol);
	fg.SetVector(FAMGDEFECT, &fgdef);
	fg.SetVector(FAMGRHS, &fgdef);		// rhs and defect are the same for smoothing/defect calculation 
	cg.SetVector(FAMGDEFECT, &cgdef);
	
	fg.Restriction(&cg);
	
	return 0;
}


int FAMG_ProlongCorrection( int fine_level, VECDATA_DESC *to, VECDATA_DESC *from, VECDATA_DESC *smooth_sol, VECDATA_DESC *smooth_def )
// fine_level in the famg grid stack
{
	INT i;

	FAMGGrid &fg = *FAMG_GetSystem()->GetMultiGrid(0)->GetGrid(fine_level);
	FAMGGrid &cg = *FAMG_GetSystem()->GetMultiGrid(0)->GetGrid(fine_level+1);
	
	assert(to==smooth_sol);
	//assert(VD_IS_SCALAR(from)); has not to be fulfilled
	//assert(VD_IS_SCALAR(to)); has not to be fulfilled
	// now it is sure we can use for np->smooth_sol, from and to the same component

	FAMGugVector fgsol(*(const FAMGugGridVector*)&fg.GetGridVector(),to);
	FAMGugVector fgdef(*(const FAMGugGridVector*)&fg.GetGridVector(),smooth_def);
	FAMGugVector c(*(const FAMGugGridVector*)&fg.GetGridVector(),from);
	FAMGugVector cgsol(*(const FAMGugGridVector*)&cg.GetGridVector(),from);

	// be sure not to use an out-dated pointer
	for (i=0; i<FAMGMAXVECTORS; i++ )
	{
		fg.SetVector(i, NULL);
		cg.SetVector(i, NULL);
	}
	
	fg.SetVector(FAMGUNKNOWN, &fgsol);
	fg.SetVector(FAMGDEFECT, &fgdef);
	fg.SetVector(FAMGRHS, &fgdef);		// rhs and defect are the same for smoothing/defect calculation 
	cg.SetVector(FAMGUNKNOWN, &cgsol);
	
	fg.Prolongation(&cg,&c);
	
	return 0;
}


int FAMG_GetN(int level)
{
    int n;

    n = FAMG_GetSystem()->GetMultiGrid(0)->GetGrid(level)->GetN();
    return n;
}

int FAMG_GetNF(int level)
{
    int nf;

    nf = FAMG_GetSystem()->GetMultiGrid(0)->GetGrid(level)->GetNF();
    return nf;
}

FAMGSystem *FAMG_GetSystem()
{
	return famgsystemptr;
}

#ifdef UG_DRAW
FAMG_Matrix* FAMG_GetMatrixPtr(int level,int i)
{
    FAMGMatrix *matrix;
    FAMGMatrixPtr mat;
    
    matrix = famgsystemptr->GetMultiGrid(0)->GetGrid(level)->GetMatrix();
    mat = matrix->GetStart(i);
    return (*((FAMG_MatrixPtr *)&mat));
}

FAMG_TransferEntry* FAMG_GetTransferEntry(int level,FAMGVectorEntry &row)
{
    FAMGTransfer *transfer;
    FAMGTransferEntry *trans;
    
    transfer = famgsystemptr->GetMultiGrid(0)->GetGrid(level)->GetTransfer();
    if(transfer == NULL) return NULL;
    trans = transfer->GetFirstEntry(row);
    return ((FAMG_TransferEntry *) trans);
}

FAMGVector* FAMG_GetVector(int level, int i)
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
#endif
