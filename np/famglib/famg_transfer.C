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
#include "dlmgr.h"
	#ifdef ModelP
	#include "parallel.h"
	#include "np.h"
	#endif
}
#endif

/* RCS_ID
$Header$
*/

//
// Class FAMGTransferEntry
//

// define the 2 static constants
const int FAMGTransferEntry::PROLONGATION_COMP = 0;
const int FAMGTransferEntry::RESTRICTION_COMP = 1;

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
#ifdef FAMG_SPARSE_BLOCK
    const FAMGSparseVector *sptmp = grid->GetGraph()->Get_spPtr();
    const FAMGSparseVector *srtmp = grid->GetGraph()->Get_srPtr();

	short MType = MATRIXTYPE(NODEVEC,NODEVEC); // only node vectors !!!
    short availsize = FMT_S_IMAT_TP(MGFORMAT(MYMG(mygrid)),MType);

    short needsize = (sptmp->Get_maxcomp()+1 + srtmp->Get_maxcomp()+1)*sizeof(double);
    if(needsize > availsize)
    {
        ostrstream ostr;
        ostr << __FILE__ << ", line " <<  __LINE__  << "size of transfer matrix too small !" << endl;
        FAMGError(ostr);
        assert(0);
    }
    
    sp.Init(sptmp,0);
    sr.Init(srtmp,sptmp->Get_maxcomp()+1);
        
#endif

	return 0;
#else
    FAMGTransferEntry *rowi;
    int i;

    n = grid->GetN();

    row_array = (FAMGTransferEntry*) FAMGGetMem(n*sizeof(FAMGTransferEntry*), FAMG_FROM_TOP);
    if (row_array == NULL) RETURN(1);
	
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

#ifdef FAMG_SPARSE_BLOCK
int FAMGTransfer::SetEntries(const FAMGVectorEntry& fg_vec, const FAMGVectorEntry& cg_vec, const FAMGSparseVector *sploc, const FAMGSparseVector *srloc, double *prolongation_val, double *restriction_val)
{

    FAMGTransferEntry *transij, *row_entry;

    // test 
    short i; short small = 1;
    for(i = 0; i < sploc->Get_n(); i++)
    {
        if (Abs(prolongation_val[sploc->Get_comp(i)]) > 1e-20) {small = 0; break;}
    }
    if(small)
    {
        for(i = 0; i < srloc->Get_n(); i++)
        {
            if (Abs(restriction_val[srloc->Get_comp(i)]) > 1e-20) {small = 0; break;}
        }
    }
    if(small) return 0; // nothing to save.

 
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
	
	transij->SetTransferEntry(&sp,sploc,prolongation_val);
	transij->SetTransferEntry(&sr,srloc,restriction_val);
	
    return 0;
}


void FAMGTransferEntry::SetTransferEntry(const FAMGSparseVector *st, const FAMGSparseVector *stloc, double *val)
{
    short ncmp = st->Get_n();
    for(short i = 0; i < ncmp; i++)
    {
        MVALUE((MATRIX*)this,st->Get_comp(i)) = val[stloc->Get_comp(i)]; 
    }
    return;
}
    
void FAMGTransferEntry::SetTransferEntry(const FAMGSparseVector *st, double val)
{
    short ncmp = st->Get_n();
    for(short i = 0; i < ncmp; i++)
    {
        MVALUE((MATRIX*)this,st->Get_comp(i)) = 1.0; 
    }
    return;
}
    

#else
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
#endif

#ifdef ModelP
#ifdef XFERTIMING
static int LocalNr;
static int CountInterfaceLengthCB(DDD_OBJ obj)
{
	LocalNr++;
	return 0;
}
#endif
#endif

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

	#ifdef XFERTIMING
	DOUBLE time;
	#endif

	DDD_IdentifyBegin();
	DDD_PrioBegin();	// SETPRIO needs it
	#endif

	// first step: create the coarse grid vectors and the transfer matrix to them
	while( fg_iter(fg_ve) )
	{
		IFDEBUG(np,3)
		ugfg_vec = ((FAMGugVectorEntryRef*)fg_ve.GetPointer())->myvector();
#ifdef ModelP
		PRINTDEBUG(np,3,("%d: ind %d gid %08x prio %d n %d\n",me,VINDEX(ugfg_vec),
						 DDD_InfoGlobalId(PARHDR(ugfg_vec)),
						 DDD_InfoPriority(PARHDR(ugfg_vec)),
						 DDD_InfoNCopies(PARHDR(ugfg_vec))));
#endif
		ENDDEBUG
		
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
			VECSKIP(ugnew_vec)=VECSKIP(ugfg_vec);
			SETVCCOARSE(ugnew_vec,0);
			VSTART(ugnew_vec) = NULL;
			VISTART(ugnew_vec) = NULL;
			SETPRIO(ugnew_vec,PRIO(ugfg_vec));

			#ifdef ModelP
			if (DDD_InfoPrioCopies(PARHDR(ugfg_vec)) > 0) {
			    int *proclist = DDD_InfoProcList(PARHDR(ugfg_vec));

				proclist += 2;
				PRINTDEBUG(np,3,("%d: ind %d gid %08x n %d\n",me,VINDEX(ugfg_vec),
								 DDD_InfoGlobalId(PARHDR(ugfg_vec)),
								 DDD_InfoNCopies(PARHDR(ugfg_vec))));
				while (*proclist != -1) {
				    // if (!GHOSTPRIO(proclist[1])) original changed!
					{ 
					    PRINTDEBUG(np,3,("%d: to pe %d with prio %d\n",me,proclist[0],proclist[1]));
						DDD_IdentifyObject(PARHDR(ugnew_vec),*proclist,
										   PARHDR(ugfg_vec));
					}
					proclist += 2;
				}
				PRINTDEBUG(np,3,("%d: prio new %d, prio old %d, attr %d\n", me,
								 DDD_InfoPriority(PARHDR(ugnew_vec)),
								 DDD_InfoPriority(PARHDR(ugfg_vec)),
								 DDD_InfoAttr(PARHDR(ugnew_vec))));
			}
			#endif

			// transfer entry for coarse node must be allocated now
			imat = CreateIMatrix (mygrid, ugfg_vec, ugnew_vec);
			if(imat!=NULL)
			{
				transfc = (FAMGTransferEntry*)imat;	// dirty cast
#ifdef FAMG_SPARSE_BLOCK		
                transfc->SetTransferEntry(&sp,1.0);
                transfc->SetTransferEntry(&sr,1.0);
#else
				transfc->SetProlongation(1.0);
				transfc->SetRestriction(1.0);
#endif
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
	
	#ifdef XFERTIMING
	time = CURRENT_TIME;
	#endif

	DDD_PrioEnd();
	DDD_IdentifyEnd();	// this constructs distributed objects; afterwards XferEnd to update prio info

	#ifdef XFERTIMING
	LocalNr=0;
	DDD_IFAExecLocal( OuterVectorSymmIF, GRID_ATTR(mygrid)-1, CountInterfaceLengthCB );
	time = CURRENT_TIME - time;
	cout <<me<<": Dest lev="<<GLEVEL(mygrid)<<' '<<time<<' '<<NVEC(mygrid)<<' '<<LocalNr<<endl;
	#endif

#endif
	
	cg.GetMatrix()->GetN() = nrVec;
	cg.GetConsMatrix()->GetN() = nrVec;
	
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
