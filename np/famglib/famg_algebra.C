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

extern "C"
{
/* ug library */
#include "gm.h"        /* for data structure               */
} // extern "C"

#ifndef __SGI10__
#include "famg_system.h"
#endif
#include "famg_grid.h"
#include "famg_algebra.h"
#include "famg_graph.h"
#include "famg_sparse.h"

/* RCS_ID
$Header$
*/

#ifdef __AIX__
// xlC on AIX doesn't know the key word typename
#define typename
#endif

//
// template functions to profit by special implementations
//

#ifdef FAMG_SPARSE_BLOCK
template<class VT>
void SetValue( VT &v, double val )
{
	// typename is a new C++ keyword!
	typename VT::Iterator viter(v); 
	typename VT::VectorEntry ve; 
	short ncmp = v.GetSparseVectorPtr()->Get_n();
	short *comp = v.GetSparseVectorPtr()->Get_comp();
    double *vptr;

    
	while(viter(ve))
    {
        vptr = v.GetValuePtr(ve);
        for(short i = 0; i < ncmp; i++) vptr[comp[i]] = val;
    }
}


template<class VT>
void AddValue( VT &dest, const VT &source )
{
    // not really tesyed
	// typename is a new C++ keyword!
	typename VT::Iterator viter(dest); 
	typename VT::VectorEntry ve; 
	short ncmp_d = dest.GetSparseVectorPtr()->Get_n();
	short ncmp_s = source.GetSparseVectorPtr()->Get_n();
	short *comp_d = dest.GetSparseVectorPtr()->Get_comp();
 	short *comp_s = source.GetSparseVectorPtr()->Get_comp();
    short ncmp;
    double *vptr_d, *vptr_s;
	
    ncmp = Min(ncmp_d,ncmp_s);

	while(viter(ve))
    {
        vptr_d = dest.GetValuePtr(ve);
        vptr_s = source.GetValuePtr(ve);
		for(short i = 0; i < ncmp; i++) vptr_d[comp_d[i]] += vptr_s[comp_s[i]];
    }
}

template<class VT>
void AddScaledValue( VT &dest, double scale, const VT &source )
{
    // not really tesyed
	// typename is a new C++ keyword!
	typename VT::Iterator viter(dest); 
	typename VT::VectorEntry ve; 
	short ncmp_d = dest.GetSparseVectorPtr()->Get_n();
	short ncmp_s = source.GetSparseVectorPtr()->Get_n();
	short *comp_d = dest.GetSparseVectorPtr()->Get_comp();
 	short *comp_s = source.GetSparseVectorPtr()->Get_comp();
    short ncmp;
    double *vptr_d, *vptr_s;
	
    ncmp = Min(ncmp_d,ncmp_s);

	while(viter(ve))
    {
        vptr_d = dest.GetValuePtr(ve);
        vptr_s = source.GetValuePtr(ve);
		for(short i = 0; i < ncmp; i++) vptr_d[comp_d[i]] += scale*vptr_s[comp_s[i]];
    }
}

template<class VT>
void SubtractValue( VT &dest, const VT &source )
{
    // not really tesyed
	// typename is a new C++ keyword!
	typename VT::Iterator viter(dest); 
	typename VT::VectorEntry ve; 
	short ncmp_d = dest.GetSparseVectorPtr()->Get_n();
	short ncmp_s = source.GetSparseVectorPtr()->Get_n();
	short *comp_d = dest.GetSparseVectorPtr()->Get_comp();
 	short *comp_s = source.GetSparseVectorPtr()->Get_comp();
    short ncmp;
    double *vptr_d, *vptr_s;
	
    ncmp = Min(ncmp_d,ncmp_s);

	while(viter(ve))
    {
        vptr_d = dest.GetValuePtr(ve);
        vptr_s = source.GetValuePtr(ve);
		for(short i = 0; i < ncmp; i++) vptr_d[comp_d[i]] -= vptr_s[comp_s[i]];
    }
}


template<class VT>
void CopyValue( VT &dest, const VT &source )
{
    // not really tesyed
	// typename is a new C++ keyword!
	typename VT::Iterator viter(dest); 
	typename VT::VectorEntry ve; 
	short ncmp_d = dest.GetSparseVectorPtr()->Get_n();
	short ncmp_s = source.GetSparseVectorPtr()->Get_n();
	short *comp_d = dest.GetSparseVectorPtr()->Get_comp();
 	short *comp_s = source.GetSparseVectorPtr()->Get_comp();
    short ncmp;
    double *vptr_d, *vptr_s;
	
    ncmp = Min(ncmp_d,ncmp_s);

	while(viter(ve))
    {
        vptr_d = dest.GetValuePtr(ve);
        vptr_s = source.GetValuePtr(ve);
		for(short i = 0; i < ncmp; i++) vptr_d[comp_d[i]] = vptr_s[comp_s[i]];
    }
}



template<class VT>
double norm( const VT& v )
{
    assert(0); // todo: adapt to sparse matrix
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
    assert(0);// todo: adapt to sparse matrix
	// typename is a new C++ keyword!
	typename VT::Iterator viter(v); 
	typename VT::VectorEntry ve; 
	register double res=0.0;
	
	while(viter(ve))
		res += v[ve]*w[ve]; 
	
#ifdef ModelP
	res = UG_GlobalSumDOUBLE( res );
#endif

	return res;
}


template<class VT>
double sum( const VT& v )
{
    assert(0);// todo: adapt to sparse matrix
	// typename is a new C++ keyword!
	typename VT::Iterator viter(v); 
	typename VT::VectorEntry ve; 
	register double res=0.0;
	
	while(viter(ve))
		res += v[ve];
	
	return res;
}

template<class VT>
void Scale( VT &v, double scale )
{
	// typename is a new C++ keyword!
	typename VT::Iterator viter(v); 
	typename VT::VectorEntry ve; 
	short ncmp = v.GetSparseVectorPtr()->Get_n();
	short *comp = v.GetSparseVectorPtr()->Get_comp();
    double *vptr;
    
	while(viter(ve))
    {
        vptr = v.GetValuePtr(ve);
        for(short i = 0; i < ncmp; i++) vptr[comp[i]] *= scale;
    }
}
template<class VT,class MT>
void VecMinusMatVec( VT &d, const VT &f, const MT &M, const VT &u )
{
	typename VT::Iterator viter(d); 
	typename VT::VectorEntry row; 
	typename MT::MatrixEntry col; 
	double *dptr, *fptr, *uptr, *mptr;
	const FAMGSparseVector *svu  = u.GetSparseVectorPtr();
	const FAMGSparseVector *svf  = f.GetSparseVectorPtr();
	const FAMGSparseVector *svd  = d.GetSparseVectorPtr();
	const FAMGSparseBlock *sb  = M.GetSparseBlockPtr();
    const FAMGSparseBlock *sbd  = M.GetDiagSparseBlockPtr();
    FAMGSparseVector svsum_d, svsum_o;

    svsum_d.Product(sbd,svu);
    svsum_o.Product(sb,svu);

    double *sum_d = new double[svsum_d.Get_maxcomp()+1];
    double *sum_o = new double[svsum_o.Get_maxcomp()+1];
	
	while(viter(row))
	{
		typename MT::Iterator miter(M,row);
		
        dptr = d.GetValuePtr(row);
        fptr = f.GetValuePtr(row);
        
        // diagonal 
        miter(col);
        uptr = u.GetValuePtr(col.dest());
        mptr = M.GetValuePtr(col);
        SparseBlockVSet(&svsum_d,sum_d,0.0);
        SparseBlockVSet(&svsum_o,sum_o,0.0);
        SparseBlockMVAddProduct(&svsum_d,sbd,svu,sum_d,mptr,uptr,1.0);
		while(miter(col))
        {
            uptr = u.GetValuePtr(col.dest());
            mptr = M.GetValuePtr(col);
            SparseBlockMVAddProduct(&svsum_o,sb,svu,sum_o,mptr,uptr,1.0);
        }
        SparseBlockVSub(svd,svf,&svsum_o,dptr,fptr,sum_o);
        SparseBlockVSub(svd,svd,&svsum_d,dptr,dptr,sum_d);
	}

    delete sum_d;
    delete sum_o;
}

template<class VT,class MT>
void JacobiSmoothFG( VT &sol, const MT &A, const VT &def )
// changes only the fine unknowns
// result in sol_vec; def_vec: correct defect before call, after call destroyed
{
	typename VT::Iterator viter(sol); 
	typename VT::VectorEntry ve; 
	const FAMGSparseVector *svsol  = sol.GetSparseVectorPtr();
	const FAMGSparseVector *svdef  = def.GetSparseVectorPtr();
	const FAMGSparseBlock *sb  = A.GetDiagSparseBlockPtr();
	double *solptr, *defptr, *matptr;

    short nr = sb->Get_nr();
    if(nr != sb->Get_nc()) assert(0);
    if(nr != svsol->Get_n()) assert(0);
    if(nr != svdef->Get_n()) assert(0);

    // todo: implement for more general vectors
    for(short i = 1; i < nr; i++)
    {
        if(svsol->Get_comp(i) - svsol->Get_comp(i-1) != 1) assert(0);
        if(svdef->Get_comp(i) - svdef->Get_comp(i-1) != 1) assert(0);
    }

    short sol_off = svsol->Get_comp(0);
    short def_off = svdef->Get_comp(0);

    double *decomp = new double[nr*nr]; 
    short *pivotmap = new short[nr]; 

	while(viter(ve))
    {
		if( sol.IsFG(ve) )
        {
            solptr = sol.GetValuePtr(ve)+sol_off;
            defptr = def.GetValuePtr(ve)+def_off;
            matptr = A.GetDiagValuePtr(ve);
            SparseBlockMCopyDense(decomp,sb,matptr);
            if(LR_Decomp(nr,decomp,pivotmap)) assert(0);
            if(LR_Solve(nr,decomp,pivotmap,solptr,defptr)) assert(0);
        }
    }

    delete decomp;
    delete pivotmap;
	
	return;
}

template<class VT,class MT>
void JacobiSmoothFGSimple( VT &sol, const MT &D, const VT &def )
// changes only the fine unknowns
// result in sol_vec; def_vec: correct defect before call, after call destroyed
{
	typename VT::Iterator viter(sol); 
	typename VT::VectorEntry ve; 
	const FAMGSparseVector *svsol  = sol.GetSparseVectorPtr();
	const FAMGSparseVector *svdef  = def.GetSparseVectorPtr();
	const FAMGSparseBlock *sb  = D.GetDiagSparseBlockPtr();
	double *solptr, *defptr, *matptr;

	while(viter(ve))
    {
		if( sol.IsFG(ve) )
        {
            solptr = sol.GetValuePtr(ve);
            defptr = def.GetValuePtr(ve);
            matptr = D.GetDiagValuePtr(ve);
            SparseBlockMVAddProduct(svsol,sb,svdef,solptr,matptr,defptr,1.0);
        }
    }
	
	return;
}


#else

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
	
	while(viter(ve))
		res += v[ve]*w[ve]; 
	
#ifdef ModelP
	res = UG_GlobalSumDOUBLE( res );
#endif

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
void JacobiSmoothFG( VT &sol, const MT &M, const VT &def )
// changes only the fine unknowns
// result in sol_vec; def_vec: correct defect before call, after call destroyed
{
	typename VT::Iterator viter(sol); 
	typename VT::VectorEntry ve; 

	while(viter(ve))
		if( sol.IsFG(ve) )
			sol[ve] += def[ve] / M.DiagValue(ve);
	
	return;
}

#endif

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

#ifdef FAMG_SPARSE_BLOCK
template<class MT>
void MarkStrongLinks(const MT &A, const FAMGGrid &grid)
{
	typedef typename MT::Vector VT;
	const typename MT::GridVector& gridvec = (typename MT::GridVector&)grid.GetGridVector();
	typename MT::MatrixEntry matij;
	typename VT::VectorEntry vi;
	typename VT::Iterator viter(gridvec);
    
     
	while (viter(vi))
	{
        typename MT::Iterator mij_iter(A,vi);
        mij_iter(matij);
        matij.set_strong(1);
        while( mij_iter(matij) )
        {
            // if(VSKIPME(matij.dest().myvector(),0)) matij.set_strong(0);
            // else matij.set_strong(1);
             matij.set_strong(1);
        }

    }
    
	return;
}
#else
template<class MT>
void MarkStrongLinks(const MT &A, const FAMGGrid &grid)
{
	typedef typename MT::Vector VT;
	const typename MT::GridVector& gridvec = (typename MT::GridVector&)grid.GetGridVector();
	typename MT::MatrixEntry matij;
	typename VT::VectorEntry vi;
	typename VT::Iterator viter(gridvec);

    double rlist[20], llist[20], mij, mji, rmax, lmax;
    int z, y;
    const double sigma = FAMGGetParameter()->Getsigma();
    const int minsl = 2 - 1;

	while (viter(vi))
	{
        for(z = 0; z <= minsl; z++)
        {
            rlist[z] = llist[z] = 0.0;
        }

        typename MT::Iterator mij_iter(A,vi);
        mij_iter(matij); // skip diagonal
        while( mij_iter(matij) )
        {
            mij = Abs(A[matij]);
            mji = Abs(A.GetAdjData(matij));

            for(z = minsl; z >= 0; z--)
            {
                if (mij < rlist[z]) break;
            }
            for(y = minsl; y > z+1; y--)
            { 
                rlist[y] = rlist[y-1];
            }
            if(z+1 <= minsl) rlist[z+1] = mij;

            for(z = minsl; z >= 0; z--)
            {
                if (mji < llist[z]) break;
            }
            for(y = minsl; y > z+1; y--)
            { 
                llist[y] = llist[y-1];
            }
            if(z+1 <= minsl) llist[z+1] = mji;            
        }

        rmax = rlist[minsl]*sigma; 
        lmax = llist[minsl]*sigma; 

        mij_iter.reset();
        mij_iter(matij);
        matij.set_strong(1);
        while( mij_iter(matij) )
        {
            mij = Abs(A[matij]);
            mji = Abs(A.GetAdjData(matij));
            if((mij > rmax) || (mji > lmax)) 
            {
                matij.set_strong(1);
            }
            else matij.set_strong(0);
        }

    }
    
	return;
}
#endif

#ifdef FAMG_SPARSE_BLOCK

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
	const MT& Dfg = (MT&)*fg.GetDiagMatrix();
	const VT &tvA = *fg.GetVector(FAMGTVA);
	const VT &tvB = *fg.GetVector(FAMGTVB);
	typename MT::MatrixEntry mij, mis;
	typename VT::VectorEntry i_fg, i_cg, j_fg, j_cg, s_fg, s_cg, t_cg;
	FAMGTransferEntry *pjs, *pij, *pst;
	typename VT::Iterator viter(fg_gridvec);
	
    
    // cast because GetSparseBlockPtr returns a const FAMGSparseBlock * pointer
    FAMGSparseBlock *cmatsb_d = (FAMGSparseBlock *)Mcg.GetDiagSparseBlockPtr();
    FAMGSparseBlock *cmatsb_o = (FAMGSparseBlock *)Mcg.GetSparseBlockPtr();

    const FAMGSparseBlock *dmatsb = Dfg.GetDiagSparseBlockPtr();
    const FAMGSparseBlock *fmatsb_o = Mfg.GetSparseBlockPtr();
    const FAMGSparseBlock *fmatsb_d = Mfg.GetDiagSparseBlockPtr();
    const FAMGSparseVector *sp = transfer.Get_sp();
    const FAMGSparseVector *sr = transfer.Get_sr();
    const FAMGSparseVector *tvAsv = tvA.GetSparseVectorPtr();
    const FAMGSparseVector *tvBsv = tvB.GetSparseVectorPtr();
    double *tvAptr, *tvBptr; 

    FAMGSparseBlock sb_o_p, sb_r_o, sb_r_o_p, sb_r_d_p, sb_r_dmat_p; // only offdiagonal blocks

    sb_o_p.Product(fmatsb_o,sp);
    sb_r_o.Product(sr,fmatsb_o);
    sb_r_o_p.Product(sr,fmatsb_o,sp);
    // sb_r_dmat_p.Product(sr,dmatsb,sp);
    sb_r_dmat_p = (*fmatsb_o);
    sb_r_d_p.Product(sr,fmatsb_d,sp);
    

    // chech sparse block structure
    if(cmatsb_o->CheckStructureforAdd(fmatsb_o)) return 1;
    if(cmatsb_o->CheckStructureforAdd(&sb_o_p)) return 1;
    if(cmatsb_o->CheckStructureforAdd(&sb_r_o)) return 1;
    if(cmatsb_o->CheckStructureforAdd(&sb_r_o_p)) return 1;
    if(cmatsb_o->CheckStructureforAdd(&sb_r_dmat_p)) return 1;
    if(cmatsb_d->CheckStructureforAdd(fmatsb_d)) return 1;
    if(cmatsb_d->CheckStructureforAdd(&sb_r_d_p)) return 1;
    if(cmatsb_d->CheckStructureforAdd(&sb_o_p)) return 1;
    if(cmatsb_d->CheckStructureforAdd(&sb_r_o)) return 1;
    if(cmatsb_d->CheckStructureforAdd(&sb_r_o_p)) return 1;
    if(cmatsb_d->CheckStructureforAdd(&sb_r_dmat_p)) return 1;


    short maxoffset = sb_o_p.Get_maxoffset();
    maxoffset = Max(maxoffset,sb_r_o.Get_maxoffset());
    maxoffset = Max(maxoffset,sb_r_o_p.Get_maxoffset());
    maxoffset = Max(maxoffset,sb_r_dmat_p.Get_maxoffset());
    maxoffset = Max(maxoffset,sb_r_d_p.Get_maxoffset());

    double *val = new double[maxoffset+1];
    double *diaginv = new double[dmatsb->Get_maxoffset()+1];


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
					// Mcg.AddEntry(Mfg[mij], i_cg, j_cg);               // Mcc
					if(i_cg == j_cg) Mcg.AddEntry(fmatsb_d,Mfg.GetValuePtr(mij), i_cg, j_cg);
                    else Mcg.AddEntry(fmatsb_o,Mfg.GetValuePtr(mij), i_cg, j_cg);    // Mcc
				}
				else
				{
					for( pjs=transfer.GetFirstEntry(j_fg); pjs != NULL; pjs = pjs->GetNext())
					{
						pjs->GetColInVar(s_cg);
                        SparseBlockMMProduct(&sb_o_p,fmatsb_o,sp,val,Mfg.GetValuePtr(mij),pjs->GetProlongationPtr());
                        Mcg.AddEntry(&sb_o_p,val,i_cg, s_cg);

						// Mcg.AddEntry(Mfg[mij]*pjs->GetProlongation(), i_cg, s_cg);      // Mcf*P
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
						// Mcg.AddEntry(pij->GetRestriction()*Mfg[mis], j_cg, s_cg);          // R*Mfc
                         SparseBlockMMProduct(&sb_r_o,sr,fmatsb_o,val,pij->GetRestrictionPtr(),Mfg.GetValuePtr(mis));
                         Mcg.AddEntry(&sb_r_o,val,j_cg, s_cg);
                       
					}
					else
					{
                        // s is fine 
                        if(s_fg == i_fg)
                        {
                            // special treatment for the A_{i,i} to keep block sparsity pattern  
                            for( pst=transfer.GetFirstEntry(s_fg); pst != NULL; pst = pst->GetNext())
                            {
                                pst->GetColInVar(t_cg);
                                // pij is equivalent to rji
                                // Mcg.AddEntry(pij->GetRestriction()*Mfg[mis]*pst->GetProlongation(), j_cg, t_cg);// R*Mff*P
                                SparseBlockMMProduct(&sb_r_d_p,sr,fmatsb_d,sp,val,pij->GetRestrictionPtr(),Mfg.GetValuePtr(mis),pst->GetProlongationPtr());
                                //Mcg.AddEntry(&sb_r_d_p,val,j_cg, j_cg); // lump to diagonal
                                Mcg.AddEntry(&sb_r_d_p,val,t_cg, t_cg); // lump to diagonal
                                
                                // todo: make sure lumping preserves filter condition
                                if(j_cg != t_cg)
                                {
                                    // SparseBlockMInvertDiag(dmatsb, diaginv, Dfg.GetValuePtr(mis));
                                    // SparseBlockMMProduct(&sb_r_dmat_p,sr,dmatsb,sp,val,pij->GetRestrictionPtr(),diaginv,pst->GetProlongationPtr());
                                    tvAptr = tvA.GetValuePtr(t_cg); tvBptr = tvB.GetValuePtr(t_cg);
                                    SparseBlockGalDiagApprox(&sb_r_dmat_p,sr,fmatsb_d,sp,tvAsv,val,pij->GetRestrictionPtr(),Mfg.GetValuePtr(mis),pst->GetProlongationPtr(),tvAptr);
                                    // SparseBlockGalDiagApproxT(&sb_r_dmat_p,sr,fmatsb_d,sp,tvBsv,val,pij->GetRestrictionPtr(),Mfg.GetValuePtr(mis),pst->GetProlongationPtr(),tvBptr);
                                    Mcg.AddEntry(&sb_r_dmat_p,val,j_cg, t_cg); 
                                    // Mcg.AddEntry(&sb_r_dmat_p,val,j_cg, j_cg,-1.0); 
                                    Mcg.AddEntry(&sb_r_dmat_p,val,t_cg, t_cg,-1.0); 
                                }
                                
                            }
						}
                        else
                        {
                            for( pst=transfer.GetFirstEntry(s_fg); pst != NULL; pst = pst->GetNext())
                            {
                                pst->GetColInVar(t_cg);
                                // pij is equivalent to rji
                                // Mcg.AddEntry(pij->GetRestriction()*Mfg[mis]*pst->GetProlongation(), j_cg, t_cg);// R*Mff*P
                                SparseBlockMMProduct(&sb_r_o_p,sr,fmatsb_o,sp,val,pij->GetRestrictionPtr(),Mfg.GetValuePtr(mis),pst->GetProlongationPtr());
                                Mcg.AddEntry(&sb_r_o_p,val,j_cg, t_cg);
                            }
                        }
					}
				}
				
			}
		}
	}

    delete val;
    delete diaginv;

	return 0;
}
#else
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
	return 0;
}
#endif
