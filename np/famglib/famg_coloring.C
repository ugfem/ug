/****************************************************************************/
/*																			*/
/* File:      famg_coloring.C												*/
/*																			*/
/* Purpose:   parallel graph coloring functions for FAMG					*/
/*																			*/
/* Author:    Christian Wrobel												*/
/*			  Institut fuer Wissenschaftliches Rechnen						*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  69120 Heidelberg												*/
/*			  internet: Christian.Wrobel@iwr.uni-heidelberg.de				*/
/*																			*/
/*																			*/
/* History:   February 99 begin, Stuttgart									*/
/*																			*/
/* Remarks:																	*/
/*																			*/
/****************************************************************************/

#include <string.h> 		// for memset

#include "gm.h"
#include "ddd.h"
#include "famg_coloring.h"

/* RCS_ID
$Header$
*/

FAMGColor FAMGMyColor;		// the color of this PE

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static DDD_PROC Nb[FAMGColorMaxNb];	// list of all neighbor PE's
static int NrNb = 0;				// number of valid entries in Nb
static int *helpNbPtr= NULL;

static int DetermineNbs( DDD_OBJ obj)
{
	VECTOR *vec = (VECTOR *)obj;
	int *proclist, i;
	
	proclist = DDD_InfoProcList(PARHDR(vec));
	for( i = 2; proclist[i] != -1	; i += 2 )
		helpNbPtr[proclist[i]] = 1;		// perhaps to much neighbors because of counting ghost-copies
}


int ConstructColoringGraph( DDD_ATTR grid_attr)
{
	int i;
	int helpNb[FAMGColorMaxProcs];
	
	assert(FAMGColorMaxProcs>=procs);	// otherwise increase the constant FAMGColorMaxProcs 
	
	// helpNb[] = 0
	memset( helpNb, 0, sizeof(int)*FAMGColorMaxProcs );
	helpNbPtr = helpNb;
	
	// determine the neighboring PE's
	DDD_IFAExecLocal( BorderVectorSymmIF, grid_attr, DetermineNbs );
	
	NrNb = 0;
	for( i = 0; i < FAMGColorMaxProcs; i++ )
		if( helpNb[i] != 0 )
			Nb[NrNb++] = i;
	
	assert(FAMGColorMaxNb>=NrNb);	// otherwise increase the constant FAMGColorMaxNb 

	return 0;
}

int ConstructColoring()
// according to 
//		Robert K. Gjertsen jr., Mark T. Jones and Paul E. Plassmann
//		Parallel Heuristics for Improved, Balanced Graph Colorings
//		to appear in journal of Parallel and Distributed Computing
{
	int i, nbpe;	
	VChannelPtr NbCh[FAMGColorMaxNb];		// communication channels to the neighbor Pe's
	double MyWeight;						// weight of this Pe (see PLF algorithm)
	DDD_PROC SendQueue[FAMGColorMaxNb];		// list of the neighbors that get a color-message from me
	int NrSend;								// number of entries in SendQueue (i.e. neighbors with a smaller weight than me)
	int NrWait;								// number of neighbors with a larger weight than me (from which I receive their color-messages)
	FAMGColor NbColor[FAMGColorMaxNb];		// list of colors of the neighbors with a larger weight than me
	
	
}
