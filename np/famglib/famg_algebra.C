/****************************************************************************/
/*																			*/
/* File:      famg_algebra.C												*/
/*																			*/
/* Purpose:   general famg vector classes functions							*/
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

#include <iostream.h>
#include <math.h>
#include <assert.h>

#include "famg_algebra.h"
#include "famg_graph.h"
#include "famg_grid.h"

/* RCS_ID
$Header$
*/

void FAMGGridVector::MarkUnknowns(FAMGGraph *graph)
{
	FAMGVectorEntry ve;
	FAMGVectorIter viter(*this);
	
	while(viter(ve))
		if (graph->GetNode(ve)->IsFGNode())	
			SetFG(ve);
		else
			SetCG(ve);
			

    return;
}

//
// template functions to profit by special implementations
//

template<class VT>
void SetValue( VT &v, double val )
{
	// typename is a new C++ keyword!
	typename VT::Iterator viter(v); 
	typename VT::VectorEntry ve; 
	
	while(viter(ve))
		v[ve] = val;
}


template<class VT>
void AddValue( VT &dest, const VT &source )
{
	// typename is a new C++ keyword!
	typename VT::Iterator viter(dest); 
	typename VT::VectorEntry ve; 
	
	while(viter(ve))
		dest[ve] += source[ve];
}


template<class VT>
void AddScaledValue( VT &dest, double scale, const VT &source )
{
	// typename is a new C++ keyword!
	typename VT::Iterator viter(dest); 
	typename VT::VectorEntry ve; 
	
	while(viter(ve))
		dest[ve] += scale * source[ve];
}


template<class VT>
void SubtractValue( VT &dest, const VT &source )
{
	// typename is a new C++ keyword!
	typename VT::Iterator viter(dest); 
	typename VT::VectorEntry ve; 
	
	while(viter(ve))
		dest[ve] -= source[ve];
}


template<class VT>
void CopyValue( VT &dest, const VT &source )
{
	// typename is a new C++ keyword!
	typename VT::Iterator viter(dest); 
	typename VT::VectorEntry ve; 
	
	if( &dest == &source )
		return; // nothing to do
	
	while(viter(ve))
		dest[ve] = source[ve];
}	

template<class VT>
double norm( const VT& v )
{
	// typename is a new C++ keyword!
	typename VT::Iterator viter(v); 
	typename VT::VectorEntry ve; 
	register double val, res=0.0;
	
	while(viter(ve))
	{
		val = v[ve];
		res += val*val; 
	}
	
#ifdef ModelP
	res = UG_GlobalSumDOUBLE( res );
#endif
	
	return sqrt(res);
}


template<class VT>
double ScalProd( const VT& v, const VT& w )
{
	// typename is a new C++ keyword!
	typename VT::Iterator viter(v); 
	typename VT::VectorEntry ve; 
	register double res=0.0;
	
#ifdef ModelP
	assert(0);	// TODO: make it for parallel
#endif
	
	while(viter(ve))
		res += v[ve]*w[ve]; 
	
	return res;
}


template<class VT>
double sum( const VT& v )
{
	// typename is a new C++ keyword!
	typename VT::Iterator viter(v); 
	typename VT::VectorEntry ve; 
	register double res=0.0;
	
	while(viter(ve))
		res += v[ve];
	
	return res;
}

template<class VT>
void Scale( VT& v, double scale )
{
	// typename is a new C++ keyword!
	typename VT::Iterator viter(v); 
	typename VT::VectorEntry ve; 
	
	while(viter(ve))
		v[ve] *= scale;	
}


template<class VT,class MT>
void VecMinusMatVec( VT &d, const VT &f, const MT &M, const VT &u )
{
	typename VT::Iterator viter(d); 
	typename VT::VectorEntry row; 
	typename MT::MatrixEntry col; 
	register double sum;
	
	while(viter(row))
	{
		typename MT::Iterator miter(M,row);
		
		sum = 0.0;
		while(miter(col))
			sum += M[col] * u[col.dest()];
		d[row] = f[row] - sum;
	}
}


template<class VT,class MT>
void MatVec( VT &dest, const MT &M, const VT &source )
{
	typename VT::Iterator viter(dest); 
	typename VT::VectorEntry row; 
	typename MT::MatrixEntry col; 
	register double sum;
	
	while(viter(row))
	{
		typename MT::Iterator miter(M,row);
		
		sum = 0.0;
		while(miter(col))
			sum += M[col] * source[col.dest()];
		dest[row] = sum;
	}
}

       
template<class VT,class MT>
void JacobiSmoothFG( VT &sol, const MT &M, const VT &def )
// changes only the fine unknowns
// result in sol_vec; def_vec: correct defect before call, after call destroyed
{
	typename VT::Iterator viter(sol); 
	typename VT::VectorEntry ve; 

	double d,m,s;
	
	while(viter(ve))
		if( sol.IsFG(ve) )
			sol[ve] += def[ve] / M.DiagValue(ve);
	
	return;
}

template<class VT,class MT>
void JacobiSmoother( VT &sol, const MT &M, const VT &def )
// result in sol_vec; def_vec: correct defect before call, after call destroyed
{
	typename VT::Iterator viter(sol); 
	typename VT::VectorEntry ve; 

	while(viter(ve))
		sol[ve] += def[ve] / M.DiagValue(ve);

	return;
}

template<class VT,class MT>
void dampedJacobiSmoother( VT &sol, const MT &M, const VT &def )
// result in sol_vec; def_vec: correct defect before call, after call destroyed
{
    static const double omega = 2.0/3.0;

	typename VT::Iterator viter(sol); 
	typename VT::VectorEntry ve; 

	while(viter(ve))
		sol[ve] += omega * def[ve] / M.DiagValue(ve);

	return;
}

template<class VT,class MT>
void FGSSmoother( VT &sol, const MT &M, VT &def )
// forward Gauss-Seidel
// result in sol_vec; def_vec: correct defect before call, after call destroyed
{
	typename VT::Iterator viter(def); 
	typename VT::VectorEntry row, col; 
	typename MT::MatrixEntry me; 
	register double sum, diag;
	int row_index;
	
	while(viter(row))
	{
		typename MT::Iterator miter(M,row);
		
		row_index = row.GetIndex();
		sum = def[row];
		miter(me);
		diag = M[me];
		while(miter(me))
		{
			col = me.dest();
			if( col.GetIndex() < row_index )
				sum -= M[me] * def[col];
		}
		def[row] = sum / diag;
	}
	
	sol += def;	// update solution
	
	return;
}

template<class VT,class MT>
void BGSSmoother( VT &sol, const MT &M, VT &def )
// backward Gauss-Seidel
// result in sol_vec; def_vec: correct defect before call, after call destroyed
{
	typename VT::RevIterator vReviter(def); 
	typename VT::VectorEntry row, col; 
	typename MT::MatrixEntry me; 
	register double sum, diag;
	int row_index;
	
	while(vReviter(row))
	{
		typename MT::Iterator miter(M,row);
		
		row_index = row.GetIndex();
		sum = def[row];
		miter(me);
		diag = M[me];
		while(miter(me))
		{
			col = me.dest();
			if( col.GetIndex() > row_index )
				sum -= M[me] * def[col];
		}
		def[row] = sum / diag;
	}
	
	sol += def;	// update solution
	
	return;
}


template<class VT,class MT>
void SGSSmoother( VT &sol, const MT &M, VT &def )
// backward Gauss-Seidel
// result in sol_vec; def_vec: correct defect before call, after call destroyed
{
	typename VT::Iterator viter(def); 
	typename VT::RevIterator vReviter(def); 
	typename VT::VectorEntry row, col; 
	typename MT::MatrixEntry me; 
	register double sum, diag;
	int row_index;

    /* symmetric Gauss-Seidel */

	while(viter(row))
	{
		typename MT::Iterator miter(M,row);
		
		row_index = row.GetIndex();
		sum = def[row];
		miter(me);
		diag = M[me];
		while(miter(me))
		{
			col = me.dest();
			if( col.GetIndex() < row_index )
				sum -= M[me] * def[col];
		}
		def[row] = sum / diag;
	}

	viter.reset();
	while(viter(row))
		def[row] *= M.DiagValue(row);
	
	while(vReviter(row))
	{
		typename MT::Iterator miter(M,row);
		
		row_index = row.GetIndex();
		sum = def[row];
		miter(me);
		diag = M[me];
		while(miter(me))
		{
			col = me.dest();
			if( col.GetIndex() > row_index )
				sum -= M[me] * def[col];
		}
		def[row] = sum / diag;
	}

	sol += def;	// update solution
	
	return;
}


template<class MT>
int ConstructGalerkinMatrix( MT &Mcg, const FAMGGrid &fg )
// this matrix lives on the coarse grid
// calculates Mcg := R * Mfg * P and with indices:
// Mcg_(i,j) := \sum_{s,t} R_(i,s) * Mfg_(s,t) * P_(t,j)
{
	typedef typename MT::Vector VT;
	
	const FAMGTransfer &transfer = *fg.GetTransfer();
	
	const typename MT::GridVector& fg_gridvec = (typename MT::GridVector&)fg.GetGridVector();
	const MT& Mfg = (MT&)*fg.GetMatrix();
	typename MT::MatrixEntry mij, mis;
	typename VT::VectorEntry i_fg, i_cg, j_fg, j_cg, s_fg, s_cg, t_cg;
	FAMGTransferEntry *pjs, *pij, *pst;

	typename VT::Iterator viter(fg_gridvec);
	
	while (viter(i_fg) )
	{
		if (fg_gridvec.IsCG(i_fg) )
		{
			// i is coarse
		
			transfer.GetFirstEntry(i_fg)->GetColInVar(i_cg);
			
			typename MT::Iterator mijiter(Mfg,i_fg);
			while( mijiter(mij) )
			{
				j_fg = mij.dest();
				
				if( fg_gridvec.IsCG(j_fg) )
				{
					transfer.GetFirstEntry(j_fg)->GetColInVar(j_cg);
					Mcg.AddEntry(Mfg[mij], i_cg, j_cg);               // Mcc
				}
				else
				{
					for( pjs=transfer.GetFirstEntry(j_fg); pjs != NULL; pjs = pjs->GetNext())
					{
						pjs->GetColInVar(s_cg);
						Mcg.AddEntry(Mfg[mij]*pjs->GetProlongation(), i_cg, s_cg);      // Mcf*P
					}
				}
			}
		}
		else
		{
			// i is fine

			typename MT::Iterator misiter(Mfg,i_fg);
			while( misiter(mis) )
			{
				s_fg = mis.dest();

				for( pij=transfer.GetFirstEntry(i_fg); pij != NULL; pij = pij->GetNext())
				{
					pij->GetColInVar(j_cg);

					if( fg_gridvec.IsCG(s_fg) )
					{
						transfer.GetFirstEntry(s_fg)->GetColInVar(s_cg);
						// pij is equivalent to rji 
						Mcg.AddEntry(pij->GetRestriction()*Mfg[mis], j_cg, s_cg);          // R*Mfc
					}
					else
					{	// s is fine 
						for( pst=transfer.GetFirstEntry(s_fg); pst != NULL; pst = pst->GetNext())
						{
							pst->GetColInVar(t_cg);
							// pij is equivalent to rji
							Mcg.AddEntry(pij->GetRestriction()*Mfg[mis]*pst->GetProlongation(), j_cg, t_cg);// R*Mff*P
						}
					}
				}
				
			}
		}
	}
}
