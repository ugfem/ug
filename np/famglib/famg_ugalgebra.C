/****************************************************************************/
/*																			*/
/* File:      famg_ugalgebra.C												*/
/*																			*/
/* Purpose:   famg ug vector classes functions								*/
/*																			*/
/* Author:    Christian Wrobel												*/
/*			  Institut fuer Computeranwendungen  III						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: ug@ica3.uni-stuttgart.de							*/
/*																			*/
/*																			*/
/* History:   August 98 begin												*/
/*																			*/
/* Remarks:																	*/
/*																			*/
/****************************************************************************/

extern "C" {
#include <math.h>
#include <gm.h>
#include <algebra.h>
}

#include <strstream.h>
#include "famg_ugalgebra.h"
#include "famg_misc.h"
#include "famg_sparse.h"

/* RCS_ID
$Header$
*/
 

#ifdef FAMG_SPARSE_BLOCK
void SetValueSkip(const FAMGugVector &v, double val )
{
	// typename is a new C++ keyword!
	FAMGugVectorIter viter(v);
	FAMGugVectorEntry ve; 
	short ncmp = v.GetSparseVectorPtr()->Get_n();
	short *comp = v.GetSparseVectorPtr()->Get_comp();
    double *vptr;

    VECTOR *vector;

    
	while(viter(ve))
    {
        vector = ve.myvector();
        vptr = v.GetValuePtr(ve);
        for(short i = 0; i < ncmp; i++) 
        {
            if(VSKIPME(vector,i)) vptr[comp[i]] = 0.0;
            else vptr[comp[i]] = val;
        }
    }
}


void SetValueSkip(const FAMGugVector &v, double *val )
{
	// typename is a new C++ keyword!
	FAMGugVectorIter viter(v);
	FAMGugVectorEntry ve; 
	short ncmp = v.GetSparseVectorPtr()->Get_n();
	short *comp = v.GetSparseVectorPtr()->Get_comp();
    double *vptr;

    VECTOR *vector;

    
	while(viter(ve))
    {
        vector = ve.myvector();
        vptr = v.GetValuePtr(ve);
        for(short i = 0; i < ncmp; i++) 
        {
            if(VSKIPME(vector,i)) vptr[comp[i]] = 0.0;
            else vptr[comp[i]] = val[i];
        }
    }
}
#else
void SetValueSkip(FAMGugVector &v, double val )
{
	// typename is a new C++ keyword!
	FAMGugVectorIter viter(v);
	FAMGugVectorEntry ve; 
    VECTOR *vector;
    
	while(viter(ve))
    {
        vector = ve.myvector();	
		if(VSKIPME(vector,0)) v[ve] = 0.0;
        else  v[ve] = val;
    }
}

#endif

#ifdef FAMG_SPARSE_BLOCK
FAMGVector* FAMGugVector::create_new() const
{
	FAMGugGridVector &gridvector = (FAMGugGridVector&)GetGridVector();
	GRID *grid = gridvector.GetGrid();

	FAMGugVector *new_v = new FAMGugVector(gridvector); 
	if(new_v!=NULL)
	{
		new_v->mydesc = NULL;
		if(AllocVDFromVD( MYMG(grid), GLEVEL(grid), GLEVEL(grid), mydesc, &new_v->mydesc)) 
		{
	        ostrstream ostr;
    	    ostr << __FILE__ << __LINE__ <<  "error in AllocVDFromVD" << endl;
        	FAMGError(ostr);
			return NULL;
    	}
		new_v->comp = VD_CMP_OF_TYPE(new_v->mydesc,0,0);
		new_v->allocatedVD = 1;
        new_v->sv.Construct(VD_NCMPS_IN_TYPE(new_v->mydesc,0),VD_CMPPTR_OF_TYPE(new_v->mydesc,0));
	}
	return new_v;
};
#else
FAMGVector* FAMGugVector::create_new() const
{
	FAMGugGridVector &gridvector = (FAMGugGridVector&)GetGridVector();
	GRID *grid = gridvector.GetGrid();

	FAMGugVector *new_v = new FAMGugVector(gridvector); 
	if(new_v!=NULL)
	{
		new_v->mydesc = NULL;
		if(AllocVDFromVD( MYMG(grid), GLEVEL(grid), GLEVEL(grid), mydesc, &new_v->mydesc)) 
		{
	        ostrstream ostr;
    	    ostr << __FILE__ << __LINE__ <<  "error in AllocVDFromVD" << endl;
        	FAMGError(ostr);
			return NULL;
    	}
		new_v->comp = VD_SCALCMP(new_v->mydesc);
		new_v->allocatedVD = 1;
	}
	return new_v;
};
#endif

#ifdef ONLY_ONE_ALGEBRA_DS
FAMGVector::~FAMGVector()
#else
FAMGugVector::~FAMGugVector()
#endif
{
	if( allocatedVD == 1 )
	{
		FAMGugGridVector &gridvector = (FAMGugGridVector&)GetGridVector();
		GRID *grid = gridvector.GetGrid();
		FreeVD(MYMG(grid), GLEVEL(grid), GLEVEL(grid), mydesc);
		mydesc = NULL;
		allocatedVD = 0;
	}
}

#ifdef WEG // wird wahrscheinlich nicht mehr benoetigt
void FAMGugMatrix::SetNumbers()
{
	VECTOR *vec;
	MATRIX *mat;
	register int n=0, nlinks=0;
	
	for( vec=FIRSTVECTOR(GetMyGrid()); vec!=NULL; vec = SUCCVC(vec) )
	{
		n++;
		for( mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat) )
			nlinks++;
	}
	
	SetN(n);
	SetNLinks(nlinks);
}
#endif

void FAMGugMatrix::AddEntry(double mval, const FAMGugVectorEntry &row, const FAMGugVectorEntry &col)
{
	FAMGugMatrixIter miter(*this, row);
	FAMGugMatrixEntry mij;
	VECTOR *ug_col_vec = col.myvector();
	MATRIX *newmat;

	while( miter(mij) )
		if ( mij.dest().myvector() == ug_col_vec )
		{
			(*this)[mij] += mval;
			return;
		}

	// the desired entry is not yet here, thus allocate it
	newmat = CreateConnection(GetMyGrid(),row.myvector(), ug_col_vec);
	if( newmat == NULL )
	{
        ostrstream ostr;
   	    ostr << __FILE__ << __LINE__ <<  "cannot allocate new matrix entry" << endl;
       	FAMGError(ostr);
		assert(0);
   	}
#ifdef DYNAMIC_MEMORY_ALLOCMODEL
	MVALUE(MADJ(newmat),GetComp())=0.0;	// circumvent bug in ug-alloc
#endif
    MVALUE(newmat,GetComp())=mval;


	GetNLinks()++;
}

#ifdef FAMG_SPARSE_BLOCK
void FAMGugMatrix::AddEntry(const FAMGSparseBlock *sbm, const double *mval, const FAMGugVectorEntry &row, const FAMGugVectorEntry &col)
{
	FAMGugMatrixIter miter(*this, row);
	FAMGugMatrixEntry mij;
	VECTOR *ug_col_vec = col.myvector();
	MATRIX *newmat;
    const FAMGSparseBlock *sb, *sbT;
    double *val, *valT;
    if(row == col) sb = GetDiagSparseBlockPtr();
    else  sb = GetSparseBlockPtr();

	while( miter(mij) )
    {
		if ( mij.dest().myvector() == ug_col_vec )
		{
	        val = GetValuePtr(mij);
            SparseBlockMMAdd(sb,sbm,val,mval);
			return;
		}
    }

	// the desired entry is not yet here, thus allocate it
	newmat = CreateConnection(GetMyGrid(),row.myvector(), ug_col_vec);
	if( newmat == NULL )
	{
        ostrstream ostr;
   	    ostr << __FILE__ << __LINE__ <<  "cannot allocate new matrix entry" << endl;
       	FAMGError(ostr);
		assert(0);
   	}

    if(row == col) 
    {
        val = GetValuePtr(newmat);
        SparseBlockMSet(sb,val,0.0);
        SparseBlockMMAdd(sb,sbm,val,mval);
    }
    else 
    {
        sbT = GetSparseBlockAdjPtr();
        val = GetValuePtr(newmat);
        valT = GetAdjValuePtr(newmat);       
        SparseBlockMSet(sb,val,0.0);
        SparseBlockMSet(sbT,valT,0.0);
        SparseBlockMMAdd(sb,sbm,val,mval);
    }



	GetNLinks()++;
}

void FAMGugMatrix::AddEntry(const FAMGSparseBlock *sbm, const double *mval, const FAMGugVectorEntry &row, const FAMGugVectorEntry &col,double factor)
{
	FAMGugMatrixIter miter(*this, row);
	FAMGugMatrixEntry mij;
	VECTOR *ug_col_vec = col.myvector();
	MATRIX *newmat;
    const FAMGSparseBlock *sb, *sbT;
    double *val, *valT;
    if(row == col) sb = GetDiagSparseBlockPtr();
    else  sb = GetSparseBlockPtr();

	while( miter(mij) )
    {
		if ( mij.dest().myvector() == ug_col_vec )
		{
	        val = GetValuePtr(mij);
            SparseBlockMMAdd(sb,sbm,val,mval,factor);
			return;
		}
    }

	// the desired entry is not yet here, thus allocate it
	newmat = CreateConnection(GetMyGrid(),row.myvector(), ug_col_vec);
	if( newmat == NULL )
	{
        ostrstream ostr;
   	    ostr << __FILE__ << __LINE__ <<  "cannot allocate new matrix entry" << endl;
       	FAMGError(ostr);
		assert(0);
   	}

    if(row == col) 
    {
        val = GetValuePtr(newmat);
        SparseBlockMSet(sb,val,0.0);
        SparseBlockMMAdd(sb,sbm,val,mval,factor);
    }
    else 
    {
        sbT = GetSparseBlockAdjPtr();
        val = GetValuePtr(newmat);
        valT = GetAdjValuePtr(newmat);       
        SparseBlockMSet(sb,val,0.0);
        SparseBlockMSet(sbT,valT,0.0);
        SparseBlockMMAdd(sb,sbm,val,mval,factor);
    }



	GetNLinks()++;
}
#endif
