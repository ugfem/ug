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

/* RCS_ID
$Header$
*/

FAMGVector* FAMGugVector::create_new() const
{
	FAMGugGridVector &gridvector = (FAMGugGridVector&)GetGridVector();
	GRID *grid = gridvector.GetGrid();

	FAMGugVector *new_v = new FAMGugVector(gridvector); 
	if(new_v!=NULL)
	{
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

#ifndef ONLY_ONE_ALGEBRA_DS
FAMGugVector::~FAMGugVector()
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
#endif	// ONLY_ONE_ALGEBRA_DS

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
