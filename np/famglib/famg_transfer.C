/****************************************************************************/
/*																			*/
/* File:      famg_transfer.C												*/
/*																			*/
/* Purpose:   famg transfer classes functions								*/
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

#include "famg_transfer.h"
#include "famg_heap.h"
#include "famg_misc.h"
#include "famg_grid.h"
#include "famg_algebra.h"

#ifdef USE_UG_DS
extern "C"
{
#include "gm.h"
#include "algebra.h"
}
#endif

/* RCS_ID
$Header$
*/

//
// Class FAMGTransferEntry
//

// define the 2 static constants
const int FAMGTransferEntry::PROLONGATION_COMP;
const int FAMGTransferEntry::RESTRICTION_COMP;

FAMGTransferEntry *FAMGTransferEntry::GetEntry(const FAMGVectorEntry &cg_vec)
{
    FAMGTransferEntry *transij;

    for(transij = this; transij != NULL; transij = transij->GetNext())
        if(transij->GetCol() == cg_vec) 
			return transij;

    return NULL;
}

//
// Class FAMGTransfer
//

int FAMGTransfer::Init(FAMGGrid *grid) 
{
#ifdef USE_UG_DS
	mygrid = grid->GetugGrid();
	return 0;
#else
    FAMGTransferEntry *rowi;
    int i;

    n = grid->GetN();

    row_array = (FAMGTransferEntry*) FAMGGetMem(n*sizeof(FAMGTransferEntry*), FAMG_FROM_TOP);
    if (row_array == NULL) return(1);
	
    for(i = 0,rowi=row_array; i < n; i++)
        *row_array++ = NULL;

    return(0);
#endif
} 

FAMGTransferEntry *FAMGTransfer::NewEntry(const FAMGVectorEntry& fg_vec, const FAMGVectorEntry& cg_vec)
{
#ifdef USE_UG_DS
	MATRIX *imat;
	
	imat = CreateIMatrix (mygrid, 
						   ((FAMGugVectorEntryRef*)(fg_vec.GetPointer()))->myvector(), 
						   ((FAMGugVectorEntryRef*)(cg_vec.GetPointer()))->myvector());
	
	FAMGTransferEntry *transij = (FAMGTransferEntry*)imat; // dirty cast
	if(transij!=NULL)
	{		
		transij->SetProlongation(0.0);
		transij->SetRestriction(0.0);
	}
	assert(MDEST(imat)!=NULL);
	return transij;

#else
    FAMGTransferEntry *ptr;
    int i,j;
	FAMGarrayVectorEntryRef &row_ent = *(FAMGarrayVectorEntryRef*)(fg_vec.GetPointer);
	FAMGarrayVectorEntryRef &col_ent = *(FAMGarrayVectorEntryRef*)(cg_vec.GetPointer);
	
    ptr = (FAMGTransferEntry *) FAMGGetMem(sizeof(FAMGTransferEntry),FAMG_FROM_TOP);
    if(ptr == NULL) return(NULL);

    i = row_ent.GetId();
    j = col_ent.GetId();
	
    ptr->SetProlongation(0.0);
    ptr->SetRestriction(0.0);
    ptr->SetId(j);

	// insert new entry as first element in the list
	if(row_array[i] == NULL)
	    ptr->SetNext(NULL);
	else
	    ptr->SetNext(row_array[i]->next);
	row_array[i] = ptr;
	
    return ptr;
#endif
}

int FAMGTransfer::SetEntries(const FAMGVectorEntry& fg_vec, const FAMGVectorEntry& cg_vec,double prolongation_val, double restriction_val)
{
    FAMGTransferEntry *transij, *row_entry;

    // test 
    if( (Abs(prolongation_val) < 1e-20) && (Abs(restriction_val) < 1e-20) )
		return 0;

	row_entry = GetFirstEntry(fg_vec);
	if(row_entry==NULL)
	{
		transij = NewEntry(fg_vec,cg_vec);
		if(transij==NULL)
	    {
			ostrstream ostr;
			ostr << __FILE__ << __LINE__  << "can not create new transfer matrix entry" << endl;
			FAMGError(ostr);
			assert(0);
		}
	}
	else
	{
		transij = row_entry->GetEntry(cg_vec);
		if(transij==NULL)
	    {	// not found
			transij = NewEntry(fg_vec,cg_vec);
			if(transij==NULL)
		    {
				ostrstream ostr;
				ostr << __FILE__ << __LINE__  << "can not create new transfer matrix entry" << endl;
				FAMGError(ostr);
				assert(0);
			}
		}
	}
	
	transij->SetProlongation(prolongation_val);
	transij->SetRestriction(restriction_val);
	
    return 0;
}

int FAMGTransfer::SetDestinationToCoarse( const FAMGGrid &fg, const FAMGGrid &cg )
// until now, the destination of the tansferentries was a coarse vector in the fine grid;
// now we must change the destination to the corresponding vector in the coarse grid
// if USE_UG_DS is defined, here the coarse grid vectors are allocated
{
	const FAMGGridVector& fg_gridvec = fg.GetGridVector();
	FAMGVectorEntry fg_ve;
	FAMGVectorIter fg_iter(fg_gridvec);
	FAMGTransferEntry *transfc;
	
	int nrVec=0;
	
#ifdef USE_UG_DS
	VECTOR *ugnew_vec, *ugfg_vec;
	MATRIX *imat;
	GRID *ugfg = fg.GetugGrid();
	GRID *ugcg = cg.GetugGrid();
	
    #ifdef ModelP
	DDD_IdentifyBegin();
	DDD_XferBegin();
	#endif

	// first step: create the coarse grid vectors and the transfer matrix to them
	while( fg_iter(fg_ve) )
	{
		if( fg_gridvec.IsCG(fg_ve) )
		{	
			// create coarse vector on coarse grid as in np/algebra/amgtools.c/GenerateNewGrid()
			ugfg_vec = ((FAMGugVectorEntryRef*)fg_ve.GetPointer())->myvector();
			if (CreateVector(ugcg,VOTYPE(ugfg_vec),VOBJECT(ugfg_vec),&ugnew_vec))
		    {
				ostrstream ostr;
				ostr << __FILE__ << __LINE__  << "can not create coarse grid vector" << endl;
				FAMGError(ostr);
				assert(0);
			}
			VINDEX(ugnew_vec) = nrVec++;
			SETVCLASS(ugnew_vec,3);
			SETVNCLASS(ugnew_vec,VCLASS(ugfg_vec));
			SETNEW_DEFECT(ugnew_vec,1);
			SETFINE_GRID_DOF(ugnew_vec,0);
			SETPRIO(ugnew_vec,PRIO(ugfg_vec));
			VECSKIP(ugnew_vec)=VECSKIP(ugfg_vec);
			SETVCCOARSE(ugnew_vec,0);

			#ifdef ModelP
			if (DDD_InfoPrioCopies(PARHDR(ugfg_vec)) > 0) {
			    int *proclist = DDD_InfoProcList(PARHDR(ugfg_vec));

				proclist += 2;
				PRINTDEBUG(np,3,("%d: ind %d gid %08x n %d ",me,VINDEX(ugfg_vec),
								 DDD_InfoGlobalId(PARHDR(ugfg_vec)),
								 DDD_InfoNCopies(PARHDR(ugfg_vec))));
				while (*proclist != -1) {
				    // if (!GHOSTPRIO(proclist[1])) original changed!
					{ 
					    PRINTDEBUG(np,3,("%d: pl %d\n",me,*proclist));
						DDD_IdentifyObject(PARHDR(ugnew_vec),*proclist,
										   PARHDR(ugfg_vec));
					}
					proclist += 2;
				}
				PRINTDEBUG(np,3,(" prio %d, attr %d\n",
								 DDD_InfoPriority(PARHDR(ugnew_vec)),
								 DDD_InfoAttr(PARHDR(ugnew_vec))));

				PRINTDEBUG(np,3,("\n"));
			}
			#endif

			// transfer entry for coarse node must be allocated now
			imat = CreateIMatrix (mygrid, ugfg_vec, ugnew_vec);
			if(imat!=NULL)
			{
				transfc = (FAMGTransferEntry*)imat;	// dirty cast
		
				transfc->SetProlongation(1.0);
				transfc->SetRestriction(1.0);
			}
			else
		    {
				ostrstream ostr;
				ostr << __FILE__ << __LINE__  << "can not create coarse vector transfer entry" << endl;
				FAMGError(ostr);
				assert(0);
			}
			assert(MDEST(imat)!=NULL);
			assert(transfc->GetNext()==NULL);	// only 1 transfer entry possible for coarse node			
		}
	}

    #ifdef ModelP
	DDD_XferEnd();
	DDD_IdentifyEnd();
	#endif
	
	cg.GetMatrix()->GetN() = nrVec;
	
	// second step: reconnect the transfer matrices from fine vectors to their coarse grid parents ,
	// which has been created in the first step
	fg_iter.reset();
	while( fg_iter(fg_ve) )
		if( fg_gridvec.IsFG(fg_ve) )
		{
			ugfg_vec = ((FAMGugVectorEntryRef*)fg_ve.GetPointer())->myvector();
			for(imat = VISTART(ugfg_vec); imat!=NULL; imat = MNEXT(imat))
			{
				assert(VCCOARSE(MDEST(imat))==1);
				MDEST(imat) = MDEST(VISTART(MDEST(imat)));
			}
		}
	
#else
	assert(0); // to implement
#endif
	return 0;
}
