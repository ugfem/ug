/****************************************************************************/
/*																			*/
/* File:      famg_sparse.C			    									*/
/*																			*/
/* Purpose:   famg sparse point block functions                             */
/*																			*/
/* Author:    Christian Wagner												*/
/*			  IWR Technische Simulation             						*/
/*			  Universitaet Heidelberg										*/
/*			  INF 368            											*/
/*			  D- 69120 Heidelberg	     									*/
/*			  internet: Christian.Wagner@iwr.uni-heidelberg.de	    		*/
/*																			*/
/*																			*/
/* History:   June 99 begin, Heidelberg     								*/
/*			                                                				*/
/*																			*/
/* Remarks:					    											*/
/*																			*/
/****************************************************************************/
 
#include <string.h>
#include <strstream.h>
#include "famg_sparse.h"
#include "famg_misc.h"
#include "famg_grid.h"
#include "famg_heap.h"

#ifdef FAMG_SPARSE_BLOCK

FAMGSparseVector::FAMGSparseVector()
{
    n = 0;
    comp = NULL;
    maxcomp = -1;
}

FAMGSparseVector::~FAMGSparseVector()
{
    n = 0;
    if(comp != NULL) delete comp;
    comp = NULL;
    maxcomp = -1;
}

FAMGSparseVector::FAMGSparseVector(const FAMGSparseVector &sv)
{
    n = sv.n;
    if(sv.comp == NULL) comp = NULL;
    else
    {
        comp = new short[n];
        for(short i = 0; i < n; i++) comp[i] = sv.comp[i];
    }
    maxcomp = sv.maxcomp;
}

FAMGSparseVector::FAMGSparseVector(const short nv)
{
    n = nv;
    comp = new short[n];
    for(short i = 0; i < n; i++) comp[i] = i;
    maxcomp = n-1;
}

FAMGSparseVector::FAMGSparseVector(const short nv, const short *comps)
{
    n = nv;
    comp = new short[n];
    maxcomp = -1;
    for(int i = 0; i < n; i++) 
    {
        comp[i] = comps[i];
        if (comp[i] > maxcomp) maxcomp = comp[i];
    }
}
    
FAMGSparseVector::FAMGSparseVector(const short *comps, const short *map, const short ncmp )
{
    n = ncmp;
    comp = new short[n];
    maxcomp = -1;
    for(int i = 0; i < n; i++) 
    {
        comp[i] = comps[map[i]];
        if (comp[i] > maxcomp) maxcomp = comp[i];
    }
}
    
FAMGSparseVector &FAMGSparseVector::operator=(const FAMGSparseVector &sv)
{
    n = sv.n;
    maxcomp = sv.maxcomp;    
    if (comp != NULL) delete comp;
    comp = new short[n]; 
    for(short i = 0; i < n; i++) comp[i] = sv.comp[i];

    return *this;
}

void FAMGSparseVector::Construct(const short nv, const short *comps)
{
    n = nv;
    if(comp != NULL) delete comp;
    comp = new short[n];
    maxcomp = -1;
    for(int i = 0; i < n; i++) 
    {
        comp[i] = comps[i];
        if (comp[i] > maxcomp) maxcomp = comp[i];
    }
}

void FAMGSparseVector::Construct(const FAMGSparseVector *sv)
{
    short i,j;
    n = sv->n;
    if(comp != NULL) delete comp;
    comp = new short[n];

    for(i = 0; i < n; i++) comp[i] = -1;
    short cmp = 0;
    for(i = 0; i < n; i++)
    {
        if(comp[i] >= 0) continue;
        comp[i] = cmp;
        for(j = i+1; j < n; j++)
        {
            if(sv->comp[i] == sv->comp[j])
            {
                comp[j] = cmp;
            }
        }
        cmp++;
    }

    maxcomp = cmp-1;

}
void FAMGSparseVector::Init(const FAMGSparseVector *sv, int off)
{
    short i,j;
    n = sv->n;
    comp = new short[n];

    for(i = 0; i < n; i++) comp[i] = -1;
    short cmp = off;
    for(i = 0; i < n; i++)
    {
        if(comp[i] >= 0) continue;
        comp[i] = cmp;
        for(j = i+1; j < n; j++)
        {
            if(sv->comp[i] == sv->comp[j])
            {
                comp[j] = cmp;
            }
        }
        cmp++;
    }

    maxcomp = cmp-1;

}

void FAMGSparseVector::Product(const FAMGSparseBlock *sb, const FAMGSparseVector *sv)
{
    short i,j, ii, jj;
    // not really tested, detailled construction of comp[] maybe not necessary

    if(sb->Get_nc() != sv->Get_n()) assert(0);
    n = sb->Get_nr();

    // construct comp
    if(comp != NULL) delete comp;
    comp = new short[n];

    // map offsets to standard range
    short *hob = new short[sb->Get_ne()];
    short *hov = new short[sv->Get_n()];

    for(i = 0; i < sb->Get_ne(); i++) hob[i] = -1;
    for(i = 0; i < sv->Get_n(); i++) hov[i] = -1;
    short hb = 0;
    for(i = 0; i < sb->Get_ne(); i++)
    {
        if(hob[i] >= 0) continue;
        hob[i] = hb;
        for(j = i+1; j < sb->Get_ne(); j++)
        {
            if(sb->Get_offset(i) == sb->Get_offset(j)) hob[j] = hb;
        }
        hb++;
    }
    short hv = 0;
    for(i = 0; i < sv->Get_n(); i++)
    {
        if(hov[i] >= 0) continue;
        hov[i] = hv;
        for(j = i+1; j < sv->Get_n(); j++)
        {
            if(sv->Get_comp(i) == sv->Get_comp(j)) hov[j] = hv;
        }
        hv++;
    }
        
    short *helpmat = new short[n*hb*hv];
    memset((void*)helpmat, 0, n*hb*hv*sizeof(short));

    short equal, cj, oij, *hm;

    for(i = 0; i < n; i++)
    {
        for(ii = sb->Get_start(i); ii < sb->Get_start(i+1); ii++)
        {
            j = sb->Get_index(ii);
            oij = hob[ii];            
            cj = hov[j];
            helpmat[hb*hv*i+oij*hv+cj]++;
        }
    }
    
    short cmp = 0; 
    for(i = 0; i < n; i++)
    {  
        hm = &helpmat[hb*hv*i];
        if(hm[0] < 0) continue;
        comp[i] = cmp;
        for(j = i+1; j < n; j++)
        {
            // compare help matrices
            equal = 1;
            for(ii = 0; ii < hb; ii++)
            {
                for(jj = 0; jj < hv; jj++)
                {
                    if(hm[ii*hv+jj] != helpmat[hb*hv*j+ii*hv+jj])
                    {
                        equal = 0;
                        break;
                    }
                }
                if(!equal) break;
            }
            if(equal) 
            {
                comp[i] = cmp;
                helpmat[hb*hv*j] = -1;
            }
        }
        cmp++;
    }
    
    maxcomp = cmp-1;

    delete hob;
    delete hov;
    delete helpmat;

    return;

}


void FAMGSparseVector::ScalProdConstruct(const FAMGSparseBlock *sb)
{
    int i, j, jj, cmp, maxoffset, equal;

     n = sb->Get_nr();
     if(comp != NULL) delete comp;
     comp = new short[n];

     // save maxoffeset
     maxoffset = sb->Get_maxoffset()+1;
     short *hv = new short[maxoffset*n];
     memset((void*)hv, 0, maxoffset*n*sizeof(short));
     for(i = 0; i < n; i++)
     {
         for(jj = sb->Get_start(i); jj < sb->Get_start(i+1); jj++)
         {
             hv[i*maxoffset+sb->Get_offset(jj)]++;
         }
     }

     cmp = 0;
     for(i = 0; i < n; i++)
     {
         if(hv[i*maxoffset] < 0) continue;
         comp[i] = cmp;
         for(j = i+1; j < n; j++)
         {
             equal = 1;
             for(jj = 0; jj < maxoffset; jj++)
             {
                 if(hv[i*maxoffset+jj] != hv[j*maxoffset+jj])
                 {
                     equal = 0;
                     break;
                 }
             }
             if(equal)
             {
                 hv[j*maxoffset] = -1;
                 comp[j] = cmp;
             }
         }
         hv[i*maxoffset] = -1;
         cmp++;
     }

     maxcomp = cmp-1;
     delete hv;
     return;
}

void FAMGSparseVector::ScalProdConstruct(const FAMGSparseBlock *sb1, const FAMGSparseBlock *sb2)
{
    // this function only makes sense, if AdaptStructure(sb1,sb2) has been called !
    int i, j, jj, kk, cmp, maxoffset1, maxoffset2, equal;

     n = sb1->Get_nr();
     if(comp != NULL) delete comp;
     comp = new short[n];

     // save maxoffeset
     maxoffset1 = sb1->Get_maxoffset()+1;
     maxoffset2 = sb2->Get_maxoffset()+1;
     short *hm = new short[maxoffset1*maxoffset2*n];
     memset((void*)hm, 0, maxoffset1*maxoffset2*n*sizeof(short));
     for(i = 0; i < n; i++)
     {
         for(jj = sb1->Get_start(i); jj < sb1->Get_start(i+1); jj++)
         {
             // only possible because of AdaptStructure(sb1,sb2)
             hm[i*maxoffset1*maxoffset2+maxoffset1*sb2->Get_offset(jj)+sb1->Get_offset(jj)]++;
         }
     }

     cmp = 0;
     for(i = 0; i < n; i++)
     {
         if(hm[i*maxoffset1*maxoffset2] < 0) continue;
         comp[i] = cmp;
         for(j = i+1; j < n; j++)
         {
             equal = 1;
             for(jj = 0; jj < maxoffset1; jj++)
             {
                 for(kk = 0; kk < maxoffset2; kk++)
                 {
                     if(hm[i*maxoffset1*maxoffset2+kk*maxoffset1+jj] != hm[j*maxoffset1*maxoffset2+kk*maxoffset1+jj])
                     {
                         equal = 0;
                         break;
                     }
                 }
                 if(!equal) break;
             }
             if(equal)
             {
                 hm[j*maxoffset1*maxoffset2] = -1;
                 comp[j] = cmp;
             }
         }
         hm[i*maxoffset1*maxoffset2] = -1;
         cmp++;
     }

     maxcomp = cmp-1;
     delete hm;
     return;
}

void FAMGSparseVector::ConstructSparseTransfer(const FAMGSparseVector *sv1, const FAMGSparseVector *sv2, const FAMGSparseVector *sv3, const FAMGSparseVector *sv4)
{
    short cmp, equal, i, j;

    n = sv1->n; // todo ?: check sv1->n == sv2->n == sv3->n == sv4->n
    if(comp != NULL) delete comp;
    comp = new short[n];

    for(i = 0; i < n; i++) comp[i] = -1;
    cmp = 0;
    for(i = 0; i < n; i++)
    {
        if(comp[i] >= 0) continue;
        comp[i] = cmp;
        equal = 1;
        for(j = i+1; j < n; j++)
        {
            if(sv1->comp[i] != sv1->comp[j]) {equal = 0;}
            else if(sv2->comp[i] != sv2->comp[j]) {equal = 0;}
            else if(sv3->comp[i] != sv3->comp[j]) {equal = 0;}
            else if(sv4->comp[i] != sv4->comp[j]) {equal = 0;}
            
            if(equal)
            {
                comp[j] = cmp;
            }
        }
        cmp++;
    }

    maxcomp = cmp-1;

    return;
}
            

    
 
FAMGSparseBlock::FAMGSparseBlock()
{
    ne = 0;
    nr = 0;
    nc = 0;
    start = NULL;
    index = NULL;
    offset = NULL;
    maxoffset = -1;
}

FAMGSparseBlock::~FAMGSparseBlock()
{
    ne = 0;
    nr = 0;
    nc = 0;
    if (start != NULL) delete start;
    if (index != NULL) delete index;
    if (offset != NULL) delete offset;
    start = NULL;
    index = NULL;
    offset = NULL;
    maxoffset = -1;
}
 
// copy constructor, very important but expensive.
// So, make sure it is only called if really necessary
FAMGSparseBlock::FAMGSparseBlock(const FAMGSparseBlock &sb)
{
    short i;

    ne = sb.ne;
    nr = sb.nr;
    nc = sb.nc;
    maxoffset = sb.maxoffset;
    
    if(sb.start == NULL) start = NULL;
    else 
    {
        start = new short[nr+1]; 
        for(i = 0; i <= nr; i++) start[i] = sb.start[i];
    }
    if(sb.index == NULL) index = NULL;
    else 
    {
        index = new short[ne];
        for(i = 0; i < ne; i++) index[i] = sb.index[i];
    }
    if(sb.offset == NULL) offset = NULL;
    else 
    {
        offset = new short[ne]; 
        for(i = 0; i < ne; i++) offset[i] = sb.offset[i];
    }
    
    return;
}
 

FAMGSparseBlock::FAMGSparseBlock(SPARSE_MATRIX* ugsm)
{
    if(ugsm == NULL)
    {
        ne = 0;
        nr = 0;
        nc = 0;
        start = NULL;
        index = NULL;
        offset = NULL;
        maxoffset = -1;

        return;
    }
        
    ne = (short) ugsm->N;
    nr = (short) ugsm->nrows;
    nc= (short) ugsm->ncols;
    start = new short[nr+1]; 
    index = new short[ne];
    offset = new short[ne]; 

    maxoffset = -1;
    short i; 
    for(i = 0; i <= nr; i++) start[i] = (short) ugsm->row_start[i];
    for(i = 0; i < ne; i++) 
    {
        index[i] = (short) ugsm->col_ind[i];
        offset[i] = (short) ugsm->offset[i];
        if(offset[i] > maxoffset) maxoffset = offset[i];
    }
}

FAMGSparseBlock &FAMGSparseBlock::operator=(const FAMGSparseBlock &sb)
{
    ne = sb.ne;
    nr = sb.nr;
    nc = sb.nc;
    
    if (start != NULL) delete start;
    if (index != NULL) delete index;
    if (offset != NULL) delete offset;

    start = new short[nr+1]; 
    index = new short[ne];
    offset = new short[ne]; 

    short i;
    for(i = 0; i <= nr; i++) start[i] = sb.start[i];
    for(i = 0; i < ne; i++) 
    {
        index[i] = sb.index[i];
        offset[i] = sb.offset[i];
    }
    maxoffset = sb.maxoffset;

    return *this;
}

void FAMGSparseBlock::Init(SPARSE_MATRIX* ugsm)
{
    if (start != NULL) delete start;
    if (index != NULL) delete index;
    if (offset != NULL) delete offset;

    if(ugsm == NULL)
    {
        ne = 0;
        nr = 0;
        nc = 0;
        start = NULL;
        index = NULL;
        offset = NULL;

        return;
    }
        
    ne = (short) ugsm->N;
    nr = (short) ugsm->nrows;
    nc= (short) ugsm->ncols;
    start = new short[nr+1]; 
    index = new short[ne];
    offset = new short[ne]; 

    maxoffset = -1;
    short i; 
    for(i = 0; i <= nr; i++) start[i] = (short) ugsm->row_start[i];
    for(i = 0; i < ne; i++) 
    {
        index[i] = (short) ugsm->col_ind[i];
        offset[i] = (short) ugsm->offset[i];
        if(offset[i] > maxoffset) maxoffset = offset[i];
    }

    return;
}

int FAMGSparseBlock::Transposed(const FAMGSparseBlock *sb)
{
    short i, j, total, delta, off, ii;
    
    if (start != NULL) delete start;
    if (index != NULL) delete index;
    if (offset != NULL) delete offset;

    ne = sb->ne;
    nr = sb->nc;
    nc = sb->nr;
    maxoffset = sb->maxoffset;

    if(sb->start == NULL)
    {
        start = NULL;
        index = NULL;
        offset = NULL;

        return 0;
    }

    start = new short[nr+1];
    index = new short[ne];
    offset = new short[ne];

    for(i = 0; i <= nr; i++) start[i] = 0;
    for(j = 0; j < sb->ne; j++)
    {
        i = sb->index[j];
        start[i]++;
    }
    total = 0;
    for(i = 0; i <= nr; i++)
    {
        delta = start[i];
        start[i] = total;
        total += delta;
    }

    short *end = new short[nr];
    for(i = 0; i < nr; i++) end[i] = start[i];
    for(i = 0; i < sb->nr; i++)
    {
        for(ii = sb->start[i]; ii < sb->start[i+1]; ii++)
        {
            j = sb->index[ii];
            off = sb->offset[ii];
            index[end[j]] = i;
            offset[end[j]] = off;
            end[j]++;
        }
    }

    delete end;

    return 0;
}

int FAMGSparseBlock::Product(const FAMGSparseBlock *sb1, const FAMGSparseBlock *sb2)
{
    short i, j, k, ii, jj, kk;
    
    if (start != NULL) delete start;
    if (index != NULL) delete index;
    if (offset != NULL) delete offset;

    if(sb1->nc != sb2->nr)
    {
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return 1;
    }
    nr = sb1->nr;
    nc = sb2->nc;


    // construct start first
    start = new short[nr+1];
    short *help = new short[nc];

    start[0] = 0;
    for(i = 0; i < nr; i++)
    {
        memset((void *)help, 0, nc*sizeof(short));
        for(ii = sb1->start[i]; ii < sb1->start[i+1]; ii++)
        {
            j = sb1->index[ii];
            for(jj = sb2->start[j]; jj < sb2->start[j+1]; jj++)
            {
                k = sb2->index[jj];
                help[k] = 1;
            }
        }
        start[i+1] = start[i];
        for(k = 0; k < nc; k++)
        {
            if(help[k]) start[i+1]++;
        }
    }

    ne = start[nr];

    // construct index
    index = new short[ne];
    for(i = 0; i < nr; i++)
    {
        memset((void *)help, 0, nc*sizeof(short));
        for(ii = sb1->start[i]; ii < sb1->start[i+1]; ii++)
        {
            j = sb1->index[ii];
            for(jj = sb2->start[j]; jj < sb2->start[j+1]; jj++)
            {
                k = sb2->index[jj];
                help[k] = 1;
            }
        }
        ii = start[i];
        for(k = 0; k < nc; k++)
        {
            if(help[k]) 
            {
                index[ii] = k;
                ii++;
            }
        }
    }

    // construct offsets
    offset = new short[ne];
    // map offsets to standard range
    short *ho1 = new short[sb1->ne];
    short *ho2 = new short[sb2->ne];

    for(i = 0; i < sb1->ne; i++) ho1[i] = -1;
    for(i = 0; i < sb2->ne; i++) ho2[i] = -1;
    short h1 = 0;
    for(i = 0; i < sb1->ne; i++)
    {
        if(ho1[i] >= 0) continue;
        ho1[i] = h1;
        for(j = i+1; j < sb1->ne; j++)
        {
            if(sb1->offset[i] == sb1->offset[j]) ho1[j] = h1;
        }
        h1++;
    }
    short h2 = 0;
    for(i = 0; i < sb2->ne; i++)
    {
        if(ho2[i] >= 0) continue;
        ho2[i] = h2;
        for(j = i+1; j < sb2->ne; j++)
        {
            if(sb2->offset[i] == sb2->offset[j]) ho2[j] = h2;
        }
        h2++;
    }
        
    short *helpmat = new short[ne*h1*h2];
    memset((void*)helpmat, 0, ne*h1*h2*sizeof(short));

    short equal, ojk, oij, *hm;

    for(i = 0; i < nr; i++)
    {
        for(ii = start[i]; ii < start[i+1]; ii++)
        {
            k = index[ii];
            help[k] = ii;
        }
        for(ii = sb1->start[i]; ii < sb1->start[i+1]; ii++)
        {
            j = sb1->index[ii];
            oij = ho1[ii];
            for(jj = sb2->start[j]; jj < sb2->start[j+1]; jj++)
            {
                k = sb2->index[jj];
                ojk = ho2[jj];
                kk = help[k];
                helpmat[h1*h2*kk+oij*h2+ojk]++;
            }
        }
    }
    
    short off= 0; 
    for(ii = 0; ii < ne; ii++)
    {  
        hm = &helpmat[h1*h2*ii];
        if(hm[0] < 0) continue;
        offset[ii] = off;
        for(i = ii+1; i < ne; i++)
        {
            // compare help matrices
            equal = 1;
            for(j = 0; j < h1; j++)
            {
                for(k = 0; k < h2; k++)
                {
                    if(hm[j*h2+k] != helpmat[h1*h2*i+j*h2+k])
                    {
                        equal = 0;
                        break;
                    }
                }
                if(!equal) break;
            }
            if(equal) 
            {
                offset[i] = off;
                helpmat[h1*h2*i] = -1;
            }
        }
        off++;
    }
    
    maxoffset = off-1;

    delete help;
    delete ho1;
    delete ho2;
    delete helpmat;

    return 0;

}

int FAMGSparseBlock::Product(const FAMGSparseBlock *sb, const FAMGSparseVector *sv)
{
    short i, j, k, r, ii, jj, rr, oj, or, cj, cr, off;
    
    if (start != NULL) delete start;
    if (index != NULL) delete index;
    if (offset != NULL) delete offset;

    if(sb->nc != sv->Get_n())
    {
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return 1;
    }
    nr = sb->nr;
    nc = sv->Get_n();
    ne = sb->ne;


    // construct start first
    start = new short[nr+1];
    for(i = 0; i <= nr; i++) start[i] = sb->start[i];

    // construct index
    index = new short[ne];
    for(ii = 0; ii < ne; ii++) index[ii] = sb->index[ii];

    // construct offsets
    offset = new short[ne];
    for(ii = 0; ii < ne; ii++) offset[ii] = -1;

    off = 0;
    for(i = 0; i < nr; i++)
    {
        for(jj = start[i]; jj < start[i+1]; jj++)
        {
            if(offset[jj] >= 0) continue;
            offset[jj] = off;
            j = index[jj];
            oj = sb->offset[jj];
            cj = sv->Get_comp(j);
            // look for equal entries
            for(k = i; k < nr; k++)
            {
                for(rr = start[k]; rr < start[k+1]; rr++)
                {  
                    if(offset[rr] >= 0) continue;
                    r = index[rr];
                    or = sb->offset[rr];
                    cr = sv->Get_comp(r);
                    if((cr == cj) && (or == oj))
                    {
                        offset[rr] = off;
                    }
                }
            }
            off++;
        }
    }

    maxoffset = off-1;

    return 0;
}

int FAMGSparseBlock::Product(const FAMGSparseVector *sv, const FAMGSparseBlock *sb)
{
    short i, k, ii, jj, rr, oj, or, ci, ck, off;
    
    if (start != NULL) delete start;
    if (index != NULL) delete index;
    if (offset != NULL) delete offset;

    if(sb->nr != sv->Get_n())
    {
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return 1;
    }
    nr = sv->Get_n();
    nc = sb->nc;
    ne = sb->ne;


    // construct start first
    start = new short[nr+1];
    for(i = 0; i <= nr; i++) start[i] = sb->start[i];

    // construct index
    index = new short[ne];
    for(ii = 0; ii < ne; ii++) index[ii] = sb->index[ii];

    // construct offsets
    offset = new short[ne];
    for(ii = 0; ii < ne; ii++) offset[ii] = -1;

    off = 0;
    for(i = 0; i < nr; i++)
    {
        ci = sv->Get_comp(i);
        for(jj = start[i]; jj < start[i+1]; jj++)
        {
            if(offset[jj] >= 0) continue;
            offset[jj] = off;
            oj = sb->offset[jj];
            // look for equal entries
            for(k = i; k < nr; k++)
            {
                ck = sv->Get_comp(k);
                if(ck != ci) continue;
                for(rr = start[k]; rr < start[k+1]; rr++)
                {  
                    if(offset[rr] >= 0) continue;
                    or = sb->offset[rr];
                    if(or == oj)
                    {
                        offset[rr] = off;
                    }
                }
            }
            off++;
        }
    }

    maxoffset = off-1;

    return 0;
}

int FAMGSparseBlock::Product(const FAMGSparseVector *svl, const FAMGSparseBlock *sb, const FAMGSparseVector *svr)
{
    short i, j, k, r, ii, jj, rr, oj, or, cj, cr, off, ci, ck;
    
    if (start != NULL) delete start;
    if (index != NULL) delete index;
    if (offset != NULL) delete offset;

    if((sb->nc != svr->Get_n()) || (sb->nr != svl->Get_n())) 
    {
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return 1;
    }
    nr = svl->Get_n();
    nc = svr->Get_n();
    ne = sb->ne;


    // construct start first
    start = new short[nr+1];
    for(i = 0; i <= nr; i++) start[i] = sb->start[i];

    // construct index
    index = new short[ne];
    for(ii = 0; ii < ne; ii++) index[ii] = sb->index[ii];

    // construct offsets
    offset = new short[ne];
    for(ii = 0; ii < ne; ii++) offset[ii] = -1;

    off = 0;
    for(i = 0; i < nr; i++)
    {
        ci = svl->Get_comp(i);
        for(jj = start[i]; jj < start[i+1]; jj++)
        {
            if(offset[jj] >= 0) continue;
            offset[jj] = off;
            j = index[jj];
            oj = sb->offset[jj];
            cj = svr->Get_comp(j);
            // look for equal entries
            for(k = i; k < nr; k++)
            {
                ck = svl->Get_comp(k);
                if(ck != ci) continue;
                for(rr = start[k]; rr < start[k+1]; rr++)
                {  
                    if(offset[rr] >= 0) continue;
                    r = index[rr];
                    or = sb->offset[rr];
                    cr = svr->Get_comp(r);
                    if((cr == cj) && (or == oj))
                    {
                        offset[rr] = off;
                    }
                }
            }
            off++;
        }
    }

    maxoffset = off-1;

    return 0;
}

void FAMGSparseBlock::FixDiag()
{
    // make sure, that the sparse block contains the diagonal
    // all additional entries are identified

    // todo: this function is not yet tested
 
    short i, jj, newentries, new_maxoffset, found, ind;
  
    newentries = 0;
    for(i = 0; i < nc; i++)
    {
        found = 0;
        for(jj  = start[i]; jj < start[i+1]; jj++)
        {
            if (index[jj] == i) 
            { 
                found = 1;
                break;
            }
        }
        if(!found) newentries++;
    } 

    if(newentries == 0) return; // diagonal ok.
    
    short *new_start = new short[nc];
    short *new_index = new short[ne+newentries];
    short *new_offset = new short[ne+newentries];

    new_start[0] = 0; ind = 0;
    for(i = 0; i < nc; i++)
    {
        found = 0;
        for(jj  = start[i]; jj < start[i+1]; jj++)
        {
            if (index[jj] == i) 
            {
                found = 1;
                break;
            }
        } 
        if(!found)
        {
            new_index[ind] = i;
            new_offset[ind] = maxoffset+1;
            ind++;
            new_maxoffset = maxoffset+1;
        }
        for(jj  = start[i]; jj < start[i+1]; jj++)
        {
            new_index[ind] = index[jj];
            new_offset[ind] = offset[jj];
            ind++;
        }
        new_start[i+1] = ind;
    }

    delete start;
    delete index;
    delete offset;

    start = new_start;
    index = new_index;
    offset = new_offset;
    maxoffset = new_maxoffset;
    ne = ne + newentries;

    return;
}
            
void AdaptStructure(FAMGSparseBlock *sb1, FAMGSparseBlock *sb2)
{
    // not really tested yet
    short i, j, jj, kk, new_ne, found, ind, new_maxoffset1, new_maxoffset2;
    short *new_start1, *new_index1, *new_offset1, *new_start2, *new_index2, *new_offset2;

    if(sb1->Get_nr() != sb2->Get_nr())
    {
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
    }


    new_ne = 0;
    for(i = 0; i < sb1->Get_nr(); i++)
    {
        for(jj = sb1->Get_start(i); jj < sb1->Get_start(i+1); jj++)
        {
            j = sb1->Get_index(jj);
            for(kk = sb2->Get_start(i); kk < sb2->Get_start(i+1); kk++)
            {
                if(j == sb2->Get_index(kk)) 
                {
                    new_ne++;
                    break;
                }
            }
        }
    }

    
    new_start1 = new short[sb1->Get_nr()+1]; 
    new_start2 = new short[sb1->Get_nr()+1]; 
    new_index1 = new short[new_ne];
    new_index2 = new short[new_ne];
    new_offset1 = new short[new_ne]; 
    new_offset2 = new short[new_ne]; 

    new_maxoffset1 = new_maxoffset2 = -1;
    ind = 0; new_start1[0] = new_start2[0] = ind;
    for(i = 0; i < sb1->Get_nr(); i++)
    {
        for(jj = sb1->Get_start(i); jj < sb1->Get_start(i+1); jj++)
        {
            j = sb1->Get_index(jj);
            found = 0;
            for(kk = sb2->Get_start(i); kk < sb2->Get_start(i+1); kk++)
            {
                if(j == sb2->Get_index(kk)) 
                {
                    found = 1;
                    break;
                }
            }
            if(found)
            {
                new_index1[ind] = new_index2[ind] = j;
                new_offset1[ind] = sb1->Get_offset(jj);
                new_offset2[ind] = sb2->Get_offset(kk);
                if(new_offset1[ind] > new_maxoffset1) new_maxoffset1 = new_offset1[ind];
                if(new_offset2[ind] > new_maxoffset2) new_maxoffset2 = new_offset2[ind];
                ind++;
            }
        }
        new_start1[i+1] = new_start2[i+1] = ind;
    }
    
    if (sb1->Get_start() != NULL) delete sb1->Get_start();
    if (sb1->Get_index() != NULL) delete sb1->Get_index();
    if (sb1->Get_offset() != NULL) delete sb1->Get_offset();
    if (sb2->Get_start() != NULL) delete sb2->Get_start();
    if (sb2->Get_index() != NULL) delete sb2->Get_index();
    if (sb2->Get_offset() != NULL) delete sb2->Get_offset();
    
    sb1->start = new_start1;
    sb1->index = new_index1;
    sb1->offset = new_offset1;
    sb1->ne = new_ne;
    sb1->maxoffset = new_maxoffset1;
    sb2->start = new_start2;
    sb2->index = new_index2;
    sb2->offset = new_offset2;
    sb2->ne = new_ne;
    sb2->maxoffset = new_maxoffset2;

    return;
}
            



int FAMGSparseBlock::CheckStructureforAdd(const FAMGSparseBlock *sb) const
{
    short i, j, jj, oj, jjs, osj, oss, r, s;
    short sbnc, ncmax, sbnr,  sbmaxoffset;
    short *sbstart, *sbindex, *sboffset;

    // not yet really tested !

    sbnr = sb->nr; 
    sbnc = sb->nc; 
    ncmax = Max(nc,sbnc);
    sbmaxoffset = sb->maxoffset;
    sbstart = sb->start; 
    sbindex = sb->index;  
    sboffset = sb->offset;  
    

    if(nr != sbnr) 
    { 
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return 1;
    }

    short *ho = new short[maxoffset+1];
    short *hos = new short[sbmaxoffset+1];
    short *hv = new short[ncmax];

    memset((void *)ho,0,sizeof(short)*(maxoffset+1));
    for(j = 0; j <= maxoffset; j++) ho[j] = -1;
    for(j = 0; j <= sbmaxoffset; j++) hos[j] = -1;
    for(j = 0; j < ncmax; j++) hv[j] = -1;

    for(i = 0; i < nr; i++)
    {
        for(jj = start[i]; jj < start[i+1]; jj++)
        {
            j = index[jj];
            oj = offset[jj];
            hv[j] = oj;
            if(ho[oj] < 0) 
            {
                ho[oj] = i*nc+j;
            }
            else
            {
                // dest_ij is not computed, because dest_rs is equal
                // ckeck if  source_ij and source_rs are equal too
                s = ho[oj]%nc; r = (ho[oj]-s)/nc;
                osj = -1;
                for(jjs = sbstart[i]; jjs < sbstart[i+1]; jjs++)
                {
                    if(sbindex[jjs] == j) {osj = sboffset[jjs]; break;}
                }
                oss = -1;
                for(jjs = sbstart[r]; jjs < sbstart[r+1]; jjs++)
                {
                    if(sbindex[jjs] == s) {oss = sboffset[jjs]; break;}
                }
                if(osj != oss)
                { 
                    ostrstream ostr; 
                    ostr << "wrong  dimensions";
                    FAMGWrite(ostr);
                    return 1;
                }
            }
        }
        for(jj = sbstart[i]; jj < sbstart[i+1]; jj++)
        {
            j = sbindex[jj];
            osj = sboffset[jj];

            if(hv[j] < 0)
            { 
                ostrstream ostr; 
                ostr << "wrong  dimensions";
                FAMGWrite(ostr);
                return 1;
            }
        }
    }

    delete hos;
    delete ho;
    delete hv;

    return 0;
}


double SparseBlockMNorm(const FAMGSparseBlock *sb,double *a)
{
    short i, jj;
    short nr = sb->Get_nr(); 
    short *start = sb->Get_start();
    short* offset = sb->Get_offset();
    double sum, norm;

    norm = 0.0;
    for(i = 0; i < nr; i++)
    {
        sum = 0.0;
        for(jj = start[i]; jj < start[i+1]; jj++)
        {
            sum += a[offset[jj]]*a[offset[jj]];
        }
        norm += sqrt(sum);
    }

    return norm;
}
    
void SparseBlockVSet(const FAMGSparseVector *sv, double *v, double val)
{
    short i;
    short n = sv->Get_n(); 
    short *comp = sv->Get_comp();
    
    for(i = 0; i < n; i++) v[comp[i]] = val;

    return;
}

void SparseBlockMSet(const FAMGSparseBlock *sb, double *a, double val)
{
    short i;
    short n = sb->Get_ne(); 
    short *offset = sb->Get_offset();
    
    for(i = 0; i < n; i++) a[offset[i]] = val;

    return;
}


void SparseBlockMSetDiag(const FAMGSparseBlock *sb, double *a, double val)
{
    short i, jj;
    short nr = sb->Get_nr(); 
    short *start = sb->Get_start();
    short *index = sb->Get_index(); 
    short* offset = sb->Get_offset();
    
    for(i = 0; i < nr; i++)
    {
        for(jj = start[i]; jj < start[i+1]; jj++)
        {
            if(index[jj] == i) a[offset[jj]] = val;
        }
    }

    return;
}

void SparseBlockMInvertDiag(const FAMGSparseBlock *sb, double *dest, double *source)
{
    short i, jj;
    short nr = sb->Get_nr(); 
    short *start = sb->Get_start();
    short *index = sb->Get_index(); 
    short* offset = sb->Get_offset();
    
    for(i = 0; i < nr; i++)
    {
        for(jj = start[i]; jj < start[i+1]; jj++)
        {
            if(index[jj] == i) dest[offset[jj]] = 1.0/source[offset[jj]];
            else  dest[offset[jj]] = 0.0;
        }
    }

    return;
}

void SparseBlockMCopy(const FAMGSparseBlock *sd, const FAMGSparseBlock *ss, double *dest, double *source, double factor)
{
    short i;
    short ne = sd->Get_ne(); 
    short* soffset = ss->Get_offset();
    short* doffset = sd->Get_offset();
    
    if(ne != ss->Get_ne()) {return;} //error todo: error return, more checks ?

    for(i = 0; i < ne; i++) dest[doffset[i]] = factor*source[soffset[i]];

    return;
}

void SparseBlockVCopy(const FAMGSparseVector *svd, const FAMGSparseVector *svs, double *dest, double *source, double factor)
{
    short i;
    short n = svd->Get_n(); 
    short* scomp = svs->Get_comp();
    short* dcomp = svd->Get_comp();
    
    if(n > svs->Get_n()) {assert(0);} 

    for(i = 0; i < n; i++) dest[dcomp[i]] = factor*source[scomp[i]];

    return;
}

void SparseBlockVAdd(const FAMGSparseVector *svd, const FAMGSparseVector *svs1, const FAMGSparseVector *svs2, double *dest, double *source1, double *source2)
{
    short i;
    short n = svd->Get_n(); 
    short* scomp1 = svs1->Get_comp();
    short* scomp2 = svs2->Get_comp();
    short* dcomp = svd->Get_comp();
    
    if(n > svs1->Get_n()) {assert(0);} 
    if(n > svs2->Get_n()) {assert(0);} 

    for(i = 0; i < n; i++) dest[dcomp[i]] = source1[scomp1[i]] + source2[scomp2[i]];

    return;
}

void SparseBlockVSub(const FAMGSparseVector *svd, const FAMGSparseVector *svs1, const FAMGSparseVector *svs2, double *dest, double *source1, double *source2)
{
    short i;
    short n = svd->Get_n(); 
    short* scomp1 = svs1->Get_comp();
    short* scomp2 = svs2->Get_comp();
    short* dcomp = svd->Get_comp();
    
    if(n > svs1->Get_n()) {assert(0);} 
    if(n > svs2->Get_n()) {assert(0);} 

    for(i = 0; i < n; i++) dest[dcomp[i]] = source1[scomp1[i]] - source2[scomp2[i]];

    return;
}

void SparseBlockVAdd(const FAMGSparseVector *svd, const FAMGSparseVector *svs1, const FAMGSparseVector *svs2, double *dest, double *source1, double *source2, double factor)
{
    short i;
    short n = svd->Get_n(); 
    short* scomp1 = svs1->Get_comp();
    short* scomp2 = svs2->Get_comp();
    short* dcomp = svd->Get_comp();
    
    if(n > svs1->Get_n()) {assert(0);} 
    if(n > svs2->Get_n()) {assert(0);} 

    for(i = 0; i < n; i++) dest[dcomp[i]] = source1[scomp1[i]] + factor*source2[scomp2[i]];

    return;
}


int SparseBlockMMAdd(const FAMGSparseBlock *sbdest, const FAMGSparseBlock *sbsource, double *dest, const double *source)
{
    short i, j, jj, oj, compute;
    short nc, ncsource, ncmax, nr, nrsource, maxoffset;
    short *start, *startsource, *index, *indexsource, *offset, *offsetsource;

    // not yet really tested !

    nr = sbdest->Get_nr(); nrsource = sbsource->Get_nr(); 
    nc = sbdest->Get_nc(); ncsource = sbsource->Get_nc(); ncmax = Max(nc,ncsource);
    maxoffset = sbdest->Get_maxoffset();
    start = sbdest->Get_start(); startsource = sbsource->Get_start(); 
    index = sbdest->Get_index(); indexsource = sbsource->Get_index(); 
    offset = sbdest->Get_offset(); offsetsource = sbsource->Get_offset(); 
    

    if(nr != nrsource) 
    { 
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return 1;
    }

    short *ho = new short[maxoffset+1];
    short *hv = new short[ncmax];

    memset((void *)ho,0,sizeof(short)*(maxoffset+1));
    for(j = 0; j < ncmax; j++) hv[j] = -1;
    for(i = 0; i < nr; i++)
    {
        compute = 0;
        for(jj = start[i]; jj < start[i+1]; jj++)
        {
            j = index[jj];
            oj = offset[jj];
            if(!ho[oj]) // not yet computed
            {
                ho[oj] = 1;
                hv[j] = oj;
                compute = 1;
            }
        }
        if(!compute) continue;
        for(jj = startsource[i]; jj < startsource[i+1]; jj++)
        {
            j = indexsource[jj];
            if(hv[j] < 0) continue;
            dest[hv[j]] += source[offsetsource[jj]];
        }
    }

    delete ho;
    delete hv;

    return 0;
}

int SparseBlockMMAdd(const FAMGSparseBlock *sbdest, const FAMGSparseBlock *sbsource, double *dest, const double *source, double factor)
{
    short i, j, jj, oj, compute;
    short nc, ncsource, ncmax, nr, nrsource, maxoffset;
    short *start, *startsource, *index, *indexsource, *offset, *offsetsource;

    // not yet really tested !

    nr = sbdest->Get_nr(); nrsource = sbsource->Get_nr(); 
    nc = sbdest->Get_nc(); ncsource = sbsource->Get_nc(); ncmax = Max(nc,ncsource);
    maxoffset = sbdest->Get_maxoffset();
    start = sbdest->Get_start(); startsource = sbsource->Get_start(); 
    index = sbdest->Get_index(); indexsource = sbsource->Get_index(); 
    offset = sbdest->Get_offset(); offsetsource = sbsource->Get_offset(); 
    

    if(nr != nrsource) 
    { 
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return 1;
    }

    short *ho = new short[maxoffset+1];
    short *hv = new short[ncmax];

    memset((void *)ho,0,sizeof(short)*(maxoffset+1));
    for(j = 0; j < ncmax; j++) hv[j] = -1;
    for(i = 0; i < nr; i++)
    {
        compute = 0;
        for(jj = start[i]; jj < start[i+1]; jj++)
        {
            j = index[jj];
            oj = offset[jj];
            if(!ho[oj]) // not yet computed
            {
                ho[oj] = 1;
                hv[j] = oj;
                compute = 1;
            }
        }
        if(!compute) continue;
        for(jj = startsource[i]; jj < startsource[i+1]; jj++)
        {
            j = indexsource[jj];
            if(hv[j] < 0) continue;
            dest[hv[j]] += factor*source[offsetsource[jj]];
        }
    }

    delete ho;
    delete hv;

    return 0;
}


int SparseBlockMMProduct(const FAMGSparseBlock *sp, const FAMGSparseBlock *sb1, const FAMGSparseBlock *sb2, double *ap, double *a1, double *a2)
{
    short i, j, k, ii, jj, oj, compute;
    double aij, bjk;
    short nc, nc1, nc2, nr, nr1, nr2, ne, maxoffset;
    short *start, *start1, *start2, *index, *index1, *index2, *offset, *offset1, *offset2;

    nc = sp->Get_nc(); nc1 = sb1->Get_nc(); nc2 = sb2->Get_nc();
    nr = sp->Get_nr(); nr1 = sb1->Get_nr(); nr2 = sb2->Get_nr();
    ne = sp->Get_ne(); maxoffset = sp->Get_maxoffset();
    start = sp->Get_start(); start1 = sb1->Get_start(); start2 = sb2->Get_start();
    index = sp->Get_index(); index1 = sb1->Get_index(); index2 = sb2->Get_index();
    offset = sp->Get_offset(); offset1 = sb1->Get_offset(); offset2 = sb2->Get_offset();
    


    if((nc1 != nr2) || (nr != nr1) || (nc != nc2))
    { 
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return 1;
    }

    short *ho = new short[maxoffset+1];
    short *hv = new short[nc];

    for(i = 0; i < ne; i++) ap[offset[i]] = 0.0;
    memset((void *)ho,0,sizeof(short)*(maxoffset+1));
    for(j = 0; j < nc; j++) hv[j] = -1;
    for(i = 0; i < nr; i++)
    {
        compute = 0;
        for(ii = start[i]; ii < start[i+1]; ii++)
        {
            j = index[ii];
            oj = offset[ii];
            if(!ho[oj]) // not yet computed
            {
                ho[oj] = 1;
                hv[j] = oj;
                compute = 1;
            }
        }
        if(!compute) continue;
        for(ii = start1[i]; ii < start1[i+1]; ii++)
        {
            j = index1[ii];
            aij = a1[offset1[ii]];
            for(jj = start2[j]; jj < start2[j+1]; jj++)
            {
                k = index2[jj];
                bjk = a2[offset2[jj]];
                if(hv[k] > -1) ap[hv[k]] += aij*bjk;
            }
        }

        for(j = 0; j < nc; j++) hv[j] = -1;
    }

    delete ho;
    delete hv;

    return 0;
}

int SparseBlockMMAddProduct(const FAMGSparseBlock *sp, const FAMGSparseBlock *sb1, const FAMGSparseBlock *sb2, double *ap, double *a1, double *a2, double factor)
{
    short i, j, k, ii, jj, oj, compute;
    double aij, bjk;
    short nc, nc1, nc2, nr, nr1, nr2, maxoffset;
    short *start, *start1, *start2, *index, *index1, *index2, *offset, *offset1, *offset2;

    nc = sp->Get_nc(); nc1 = sb1->Get_nc(); nc2 = sb2->Get_nc();
    nr = sp->Get_nr(); nr1 = sb1->Get_nr(); nr2 = sb2->Get_nr();
    maxoffset = sp->Get_maxoffset();
    start = sp->Get_start(); start1 = sb1->Get_start(); start2 = sb2->Get_start();
    index = sp->Get_index(); index1 = sb1->Get_index(); index2 = sb2->Get_index();
    offset = sp->Get_offset(); offset1 = sb1->Get_offset(); offset2 = sb2->Get_offset();
    


    if((nc1 != nr2) || (nr != nr1) || (nc != nc2))
    { 
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return 1;
    }

    short *ho = new short[maxoffset+1];
    short *hv = new short[nc];

    memset((void *)ho,0,sizeof(short)*(maxoffset+1));
    for(j = 0; j < nc; j++) hv[j] = -1;
    for(i = 0; i < nr; i++)
    {
        compute = 0;
        for(ii = start[i]; ii < start[i+1]; ii++)
        {
            j = index[ii];
            oj = offset[ii];
            if(!ho[oj]) // not yet computed
            {
                ho[oj] = 1;
                hv[j] = oj;
                compute = 1;
            }
        }
        if(!compute) continue;
        for(ii = start1[i]; ii < start1[i+1]; ii++)
        {
            j = index1[ii];
            aij = a1[offset1[ii]];
            for(jj = start2[j]; jj < start2[j+1]; jj++)
            {
                k = index2[jj];
                bjk = a2[offset2[jj]];
                if(hv[k] > -1) ap[hv[k]] += factor*aij*bjk;
            }
        }

        for(j = 0; j < nc; j++) hv[j] = -1;
    }

    delete ho;
    delete hv;

    return 0;
}

int SparseBlockMMProduct(const FAMGSparseBlock *sd, const FAMGSparseBlock *sb, const FAMGSparseVector *sv, double *dest, double *a, double *d)
{
    short i, j, ii, oj;
    short ne, nc, nc1, nc2, nr, nr1, nr2, maxoffset;
    short *start, *index, *offset, *offset1, *comp;

    nc = sd->Get_nc(); nc1 = sb->Get_nc(); nc2 = sv->Get_n();
    nr = sd->Get_nr(); nr1 = sb->Get_nr(); nr2 = nc2;
    maxoffset = sd->Get_maxoffset(); ne = sd->Get_ne();
    start = sd->Get_start(); 
    index = sd->Get_index(); 
    offset = sd->Get_offset(); 
    offset1 = sb->Get_offset(); 
    comp = sv->Get_comp();
    


    if((nc1 != nr2) || (nr != nr1) || (nc != nc2))
    { 
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return 1;
    }

    short *computed = new short[maxoffset+1];
    memset((void *)computed,0,sizeof(short)*(maxoffset+1));

    for(i = 0; i < ne; i++) dest[offset[i]] = 0.0; // maybe we can skip this to save time
    for(i = 0; i < nr; i++)
    {
        for(ii = start[i]; ii < start[i+1]; ii++)
        {
            oj = offset[ii];
            j = index[ii];
            if(computed[oj]) continue;
            computed[oj] = 1;
            dest[oj] = a[offset1[ii]]*d[comp[j]];
        }
    }

    delete computed;

    return 0;
}

int SparseBlockMMProduct(const FAMGSparseBlock *sd, const FAMGSparseVector *sv, const FAMGSparseBlock *sb, double *dest, double *d, double *a)
{
    short i, ii, oj;
    double di;
    short ne, nc, nc1, nc2, nr, nr1, nr2, maxoffset;
    short *start, *offset, *offset1, *comp;

    nc = sd->Get_nc(); nc1 = sb->Get_nc(); nc2 = sv->Get_n();
    nr = sd->Get_nr(); nr1 = sb->Get_nr(); nr2 = nc2;
    maxoffset = sd->Get_maxoffset(); ne = sd->Get_ne();
    start = sd->Get_start(); 
    offset = sd->Get_offset(); 
    offset1 = sb->Get_offset(); 
    comp = sv->Get_comp();
    


    if((nc2 != nr1) || (nr != nr2) || (nc != nc1))
    { 
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return 1;
    }

    short *computed = new short[maxoffset+1];
    memset((void *)computed,0,sizeof(short)*(maxoffset+1));

    for(i = 0; i < ne; i++) dest[offset[i]] = 0.0; // maybe we can skip this to save time
    for(i = 0; i < nr; i++)
    {
        di = d[comp[i]];
        for(ii = start[i]; ii < start[i+1]; ii++)
        {
            oj = offset[ii];
            if(computed[oj]) continue;
            computed[oj] = 1;
            dest[oj] = di*a[offset1[ii]];
        }
    }

    delete computed;

    return 0;
}

int SparseBlockMMProduct(const FAMGSparseBlock *sd, const FAMGSparseVector *svl, const FAMGSparseBlock *sb, const FAMGSparseVector *svr, double *dest, double *dl, double *a, double *dr)
{
    short i, j, ii, oj;
    double di;
    short ne, nc, nc1, ncr, ncl, nr, nr1, nrr, nrl, maxoffset;
    short *start, *index, *offset, *offset1, *compr, *compl;

    nc = sd->Get_nc(); nc1 = sb->Get_nc(); ncr = svr->Get_n(); ncl = svl->Get_n();
    nr = sd->Get_nr(); nr1 = sb->Get_nr(); nrr = ncr; nrl = ncl;
    maxoffset = sd->Get_maxoffset(); ne = sd->Get_ne();
    start = sd->Get_start(); 
    index = sd->Get_index(); 
    offset = sd->Get_offset(); 
    offset1 = sb->Get_offset(); 
    compr = svr->Get_comp();
    compl = svl->Get_comp();
    


    if((ncl != nr1) || (nrr != nc1) || (nr != nrl) || (nc != ncr))
    { 
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return 1;
    }

    short *computed = new short[maxoffset+1];
    memset((void *)computed,0,sizeof(short)*(maxoffset+1));

    for(i = 0; i < ne; i++) dest[offset[i]] = 0.0; // maybe we can skip this to save time
    for(i = 0; i < nr; i++)
    {
        di = dl[compl[i]];
        for(ii = start[i]; ii < start[i+1]; ii++)
        {
            j = index[ii];
            oj = offset[ii];
            if(computed[oj]) continue;
            computed[oj] = 1;
            dest[oj] = di*a[offset1[ii]]*dr[compr[j]];
        }
    }

    delete computed;

    return 0;
}


int SparseBlockMVAddProduct(const FAMGSparseVector *dest, const FAMGSparseBlock *sb, const FAMGSparseVector *source, double *vd, double *a, double *vs, const double factor)
{
    short i, jj, nr, *start, *offset, *scomp, *dcomp, *index;
    double sum;

    nr = sb->Get_nr();
    if((nr != dest->Get_n()) || (sb->Get_nc() != source->Get_n()))
    {
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return 1;
    }
    start = sb->Get_start();
    offset = sb->Get_offset();
    index = sb->Get_index();
    scomp = source->Get_comp();
    dcomp = dest->Get_comp();
    
    short *computed = new short[dest->Get_maxcomp()+1];
    memset((void*)computed,0,sizeof(short)*(dest->Get_maxcomp()+1)); 

    for(i = 0; i < nr; i++)
    {
        if(computed[dcomp[i]]) continue;
        computed[dcomp[i]] = 1;
        sum = 0.0;
        for(jj = start[i]; jj < start[i+1]; jj++)
        {
            sum += a[offset[jj]]*vs[scomp[index[jj]]];
        }
        vd[dcomp[i]] = vd[dcomp[i]] + factor * sum;
    }

    delete computed;

    return 0;    
}

int SparseBlockMVProduct(const FAMGSparseVector *dest, const FAMGSparseBlock *sb, const FAMGSparseVector *source, double *vd, const double *a, const double *vs)
{
    short i, jj, nr, *start, *offset, *scomp, *dcomp, *index;
    double sum;

    nr = sb->Get_nr();
    if((nr != dest->Get_n()) || (sb->Get_nc() != source->Get_n()))
    {
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return 1;
    }
    start = sb->Get_start();
    offset = sb->Get_offset();
    index = sb->Get_index();
    scomp = source->Get_comp();
    dcomp = dest->Get_comp();
    
    short *computed = new short[dest->Get_maxcomp()+1];
    memset((void*)computed,0,sizeof(short)*(dest->Get_maxcomp()+1)); 

    for(i = 0; i < nr; i++)
    {
        if(computed[dcomp[i]]) continue;
        computed[dcomp[i]] = 1;
        sum = 0.0;
        for(jj = start[i]; jj < start[i+1]; jj++)
        {
            sum += a[offset[jj]]*vs[scomp[index[jj]]];
        }
        vd[dcomp[i]] = sum;
    }

    delete computed;

    return 0;    
}

int SparseBlockMVAddProduct(const FAMGSparseVector *dest, const FAMGSparseVector *sv, const FAMGSparseVector *source, double *vd, double *d, double *vs, const double factor)
{
    // matrix is diagonal and therefore saved as vector
    short i, n, *comp, *scomp, *dcomp;

    n = sv->Get_n();
    if((n != dest->Get_n()) || (n != source->Get_n()))
    {
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        assert(0);
    }
    comp = sv->Get_comp();
    scomp = source->Get_comp();
    dcomp = dest->Get_comp();
    
    for(i = 0; i < n; i++)
    {
        vd[dcomp[i]] += factor*d[comp[i]]*vs[scomp[i]];
    }

    return 0;    
}

int SparseBlockMVProduct(const FAMGSparseVector *dest, const FAMGSparseVector *sv, const FAMGSparseVector *source, double *vd, const double *d, const double *vs)
{
    // matrix is diagonal and therefore saved as vector
    short i, n, *comp, *scomp, *dcomp;

    n = sv->Get_n();
    if((n != dest->Get_n()) || (n != source->Get_n()))
    {
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        assert(0);
    }
    comp = sv->Get_comp();
    scomp = source->Get_comp();
    dcomp = dest->Get_comp();
    
    for(i = 0; i < n; i++)
    {
        vd[dcomp[i]] = d[comp[i]]*vs[scomp[i]];
    }

    return 0;    
}


void SparseBlockRowAddScalProd(const FAMGSparseVector *sp, const FAMGSparseBlock *sb, double *scalprod, double *a, double *b)
{
    // scalprod is not cleared !
    short k, jj;
    short maxcomp = sp->Get_maxcomp();
    short *comp = sp->Get_comp();
    short n = sp->Get_n();
    short *start = sb->Get_start();
    short *offset = sb->Get_offset();
 
    if(n != sb->Get_nr()) 
    {
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return;
    }

    short *computed = new short[maxcomp+1];
    memset((void*)computed,0,(maxcomp+1)*sizeof(short));
    
    for(k = 0; k < n; k++)
    {
        if(computed[comp[k]]) continue;
        computed[comp[k]] = 1;
        for(jj = start[k]; jj < start[k+1]; jj++)
        {
            scalprod[comp[k]] += a[offset[jj]] * b[offset[jj]];
        }
    }

    delete computed;
    return;
}

void SparseBlockRowAddScalProd(const FAMGSparseVector *sp, const FAMGSparseBlock *sb1, const FAMGSparseBlock *sb2, double *scalprod, double *a, double *b)
{
    // scalprod is not cleared !
    // make sure AdaptStructure(sb1,sb2) has been called !
    short k, jj;
    short maxcomp = sp->Get_maxcomp();
    short *comp = sp->Get_comp();
    short n = sp->Get_n();
    short *start = sb1->Get_start();
    short *offset1 = sb1->Get_offset();
    short *offset2 = sb2->Get_offset();
 
    if(n != sb1->Get_nr()) 
    {
        ostrstream ostr; 
        ostr << "wrong  dimensions";
        FAMGWrite(ostr);
        return;
    }

    short *computed = new short[maxcomp+1];
    memset((void*)computed,0,(maxcomp+1)*sizeof(short));
    
    for(k = 0; k < n; k++)
    {
        if(computed[comp[k]]) continue;
        computed[comp[k]] = 1;
        for(jj = start[k]; jj < start[k+1]; jj++)
        {
            scalprod[comp[k]] += a[offset1[jj]] * b[offset2[jj]];
        }
    }

    delete computed;
    return;
}




void SparseBlockMCopyDense(double *decomp, const FAMGSparseBlock *sb, double *matptr)
{
    short i, jj, off_i;

    short nr = sb->Get_nr();
    if(nr != sb->Get_nc()) assert(0);
    
    short *start = sb->Get_start();
    short *index = sb->Get_index();
    short *offset = sb->Get_offset();

    memset((void*)decomp,0,nr*nr*sizeof(double));

    for(i = 0; i < nr; i++)
    {
        off_i = i*nr;
        for(jj = start[i]; jj < start[i+1]; jj++)
        {
            decomp[off_i + index[jj]] = matptr[offset[jj]];
        }
    }

    return;
}
    
void SparseBlockDiagApprox(const FAMGSparseBlock *sb, const FAMGSparseBlock *sbd, const FAMGSparseVector *svt, double *mij, const double *miidecomp, const double factor, const double *tv)
{
    double aii;
    short i, j, jj;
    short nr = sb->Get_nr();
    FAMGSparseVector svh(nr);
    double *help1 = new double[nr];

    LR_Solve(nr,miidecomp+sbd->Get_offset(0),help1,tv+svt->Get_comp(0)); 
    
    for(i = 0; i < sb->Get_ne(); i++) mij[sb->Get_offset(i)] = 0.0;
    for(i = 0; i < nr; i++)
    {
        for(jj = sb->Get_start(i); jj < sb->Get_start(i+1); jj++)
        {
            j = sb->Get_index(jj);
            if(i == j)
            {
                aii = factor*help1[i]/tv[svt->Get_comp(i)];
                if(Abs(aii) > Abs(mij[sb->Get_offset(jj)]))
                {                 
                    mij[sb->Get_offset(jj)] = aii;
                }
                break;
            }
        }
    }

    delete help1;

    return;
}


void SparseBlockDiagApprox(const FAMGSparseBlock *sb, const FAMGSparseBlock *sbd, const FAMGSparseBlock *sbo, const FAMGSparseVector *svt, double *mij, const double *miidecomp, const double *matij, const double *tv)
{
    double aii;
    short i, j, jj;
    short nr = sb->Get_nr();
    FAMGSparseVector svh(nr);
    double *help = new double[nr];

    SparseBlockMVProduct(&svh,sbo,svt,help,matij,tv);
    LR_Solve(nr,miidecomp+sbd->Get_offset(0),help,help); 
    // help as in and out vector only possible without pivot search !
    
    for(i = 0; i < sb->Get_ne(); i++) mij[sb->Get_offset(i)] = 0.0;
    for(i = 0; i < nr; i++)
    {
        for(jj = sb->Get_start(i); jj < sb->Get_start(i+1); jj++)
        {
            j = sb->Get_index(jj);
            if(i == j)
            {
                aii = help[i]/tv[svt->Get_comp(i)];
                if(Abs(aii) > Abs(mij[sb->Get_offset(jj)]))
                {                 
                    mij[sb->Get_offset(jj)] = aii;
                }
                break;
            }
        }
    }

    delete help;

    return;
}
    
void SparseBlockDiagApprox(const FAMGSparseBlock *sb, const FAMGSparseBlock *sbd, const FAMGSparseBlock *sbo, const FAMGSparseVector *svt, double *mij, const double *miidecomp, const double *matij, const double *mjjdecomp, const double *tv)
{
    double aii;
    short i, j, jj;
    short nr = sb->Get_nr();
    FAMGSparseVector svh(nr);
    double *help1 = new double[nr];
    double *help2 = new double[nr];

    LR_Solve(nr,mjjdecomp+sbd->Get_offset(0),help1,tv+svt->Get_comp(0)); 
    SparseBlockMVProduct(&svh,sbo,&svh,help2,matij,help1);
    LR_Solve(nr,miidecomp+sbd->Get_offset(0),help1,help2); 
    
    for(i = 0; i < sb->Get_ne(); i++) mij[sb->Get_offset(i)] = 0.0;
    for(i = 0; i < nr; i++)
    {
        for(jj = sb->Get_start(i); jj < sb->Get_start(i+1); jj++)
        {
            j = sb->Get_index(jj);
            if(i == j)
            {
                aii = help1[i]/tv[svt->Get_comp(i)];
                if(Abs(aii) > Abs(mij[sb->Get_offset(jj)]))
                {                 
                    mij[sb->Get_offset(jj)] = aii;
                }
                break;
            }
        }
    }

    delete help1;
    delete help2;

    return;
}

void SparseBlockDiagApproxT(const FAMGSparseBlock *sb, const FAMGSparseBlock *sbd, const FAMGSparseVector *svt, double *mij, const double *miidecomp, const double factor, const double *tv)
{
    double aii;
    short i, j, jj;
    short nr = sb->Get_nr();
    FAMGSparseVector svh(nr);
    double *help1 = new double[nr];

    LR_SolveT(nr,miidecomp+sbd->Get_offset(0),help1,tv+svt->Get_comp(0)); 
    
    for(i = 0; i < sb->Get_ne(); i++) mij[sb->Get_offset(i)] = 0.0;
    for(i = 0; i < nr; i++)
    {
        for(jj = sb->Get_start(i); jj < sb->Get_start(i+1); jj++)
        {
            j = sb->Get_index(jj);
            if(i == j)
            {
                aii = factor*help1[i]/tv[svt->Get_comp(i)];
                if(Abs(aii) > Abs(mij[sb->Get_offset(jj)]))
                {                 
                    mij[sb->Get_offset(jj)] = aii;
                }
                break;
            }
        }
    }

    delete help1;

    return;
}


void SparseBlockDiagApproxT(const FAMGSparseBlock *sb, const FAMGSparseBlock *sbd, const FAMGSparseBlock *sbo, const FAMGSparseVector *svt, double *mij, const double *miidecomp, const double *matij, const double *tv)
{
    double aii;
    short i, j, jj;
    short nr = sb->Get_nr();
    FAMGSparseVector svh(nr);
    double *help = new double[nr];

    SparseBlockMVProduct(&svh,sbo,svt,help,matij,tv);
    LR_SolveT(nr,miidecomp+sbd->Get_offset(0),help,help); 
    // help as in and out vector only possible without pivot search !
    
    for(i = 0; i < sb->Get_ne(); i++) mij[sb->Get_offset(i)] = 0.0;
    for(i = 0; i < nr; i++)
    {
        for(jj = sb->Get_start(i); jj < sb->Get_start(i+1); jj++)
        {
            j = sb->Get_index(jj);
            if(i == j)
            {
                aii = help[i]/tv[svt->Get_comp(i)];
                if(Abs(aii) > Abs(mij[sb->Get_offset(jj)]))
                {                 
                    mij[sb->Get_offset(jj)] = aii;
                }
                break;
            }
        }
    }

    delete help;

    return;
}
    
void SparseBlockDiagApproxT(const FAMGSparseBlock *sb, const FAMGSparseBlock *sbd, const FAMGSparseBlock *sbo, const FAMGSparseVector *svt, double *mij, const double *miidecomp, const double *matij, const double *mjjdecomp, const double *tv)
{
    double aii;
    short i, j, jj;
    short nr = sb->Get_nr();
    FAMGSparseVector svh(nr);
    double *help1 = new double[nr];
    double *help2 = new double[nr];

    LR_SolveT(nr,mjjdecomp+sbd->Get_offset(0),help1,tv+svt->Get_comp(0)); 
    SparseBlockMVProduct(&svh,sbo,&svh,help2,matij,help1);
    LR_SolveT(nr,miidecomp+sbd->Get_offset(0),help1,help2); 
    
    for(i = 0; i < sb->Get_ne(); i++) mij[sb->Get_offset(i)] = 0.0;
    for(i = 0; i < nr; i++)
    {
        for(jj = sb->Get_start(i); jj < sb->Get_start(i+1); jj++)
        {
            j = sb->Get_index(jj);
            if(i == j)
            {
                aii = help1[i]/tv[svt->Get_comp(i)];
                if(Abs(aii) > Abs(mij[sb->Get_offset(jj)]))
                {                 
                    mij[sb->Get_offset(jj)] = aii;
                }
                break;
            }
        }
    }

    delete help1;
    delete help2;

    return;
}
    
void SparseBlockGalDiagApprox(const FAMGSparseBlock *sb, const FAMGSparseVector *sbr, const FAMGSparseBlock *sbd, const FAMGSparseVector *sbp, const FAMGSparseVector *svt, double *mij, const double *ris, const double *mss, const double *psj, const double *tv)
{
    double aii;
    short i, j, jj;
    short nr = sb->Get_nr();
    FAMGSparseVector svh(nr);
    double *help1 = new double[nr];
    double *help2 = new double[nr];


    SparseBlockMVProduct(&svh,sbp,svt,help1,psj,tv);
    SparseBlockMVProduct(&svh,sbd,&svh,help2,mss,help1);
    SparseBlockMVProduct(&svh,sbr,&svh,help1,ris,help2);
    

    for(i = 0; i < sb->Get_ne(); i++) mij[sb->Get_offset(i)] = 0.0;
    for(i = 0; i < nr; i++)
    {
        for(jj = sb->Get_start(i); jj < sb->Get_start(i+1); jj++)
        {
            j = sb->Get_index(jj);
            if(i == j)
            {
                aii = help1[i]/tv[svt->Get_comp(i)];
                if(Abs(aii) > Abs(mij[sb->Get_offset(jj)]))
                {                 
                    mij[sb->Get_offset(jj)] = aii;
                }
                break;
            }
        }
    }
    
    delete help1;
    delete help2;

    return;
}

void SparseBlockGalDiagApproxT(const FAMGSparseBlock *sb, const FAMGSparseVector *sbr, const FAMGSparseBlock *sbd, const FAMGSparseVector *sbp, const FAMGSparseVector *svt, double *mij, const double *ris, const double *mss, const double *psj, const double *tv)
{
    double aii;
    short i, j, jj;
    short nr = sb->Get_nr();
    FAMGSparseVector svh(nr);
    FAMGSparseBlock sbdT;
    double *help1 = new double[nr];
    double *help2 = new double[nr];

    sbdT.Transposed(sbd);

    SparseBlockMVProduct(&svh,sbr,svt,help1,ris,tv);
    SparseBlockMVProduct(&svh,&sbdT,&svh,help2,mss,help1);
    SparseBlockMVProduct(&svh,sbp,&svh,help1,psj,help2);
   

    for(i = 0; i < sb->Get_ne(); i++) mij[sb->Get_offset(i)] = 0.0;
    for(i = 0; i < nr; i++)
    {
        for(jj = sb->Get_start(i); jj < sb->Get_start(i+1); jj++)
        {
            j = sb->Get_index(jj);
            if(i == j)
            {
                aii = help1[i]/tv[svt->Get_comp(i)];
                if(Abs(aii) > Abs(mij[sb->Get_offset(jj)]))
                {                 
                    mij[sb->Get_offset(jj)] = aii;
                }
                break;
            }
        }
    }
    
    delete help1;
    delete help2;

    return;
}



// lump the diagonal block
int FAMGGrid::ConstructDiagonalLump()
{
    FAMGGridVector gv = GetGridVector();
    FAMGMatrixAlg &A = *GetMatrix();
    FAMGMatrixAlg &D = *GetDiagMatrix();
    FAMGVectorEntry vi;
    FAMGVectorIter viter(gv);
    double *Dii, *Aii, sum;
    short i, jj;

    const FAMGSparseBlock *Dsb = D.GetDiagSparseBlockPtr();
    const FAMGSparseBlock *Asb = A.GetDiagSparseBlockPtr();


    short *D_diag_offset = new short[Dsb->Get_nr()]; 
 
    if(Dsb->Get_nr() != Asb->Get_nr())
    {
        ostrstream ostr; 
        ostr << "wrong block structure of the diagmatrix";
        FAMGWrite(ostr);
        return(1);
    }
    for(i = 0; i < Dsb->Get_nr(); i++) D_diag_offset[i] = -1;
    for(i = 0; i < Dsb->Get_nr(); i++)
    {
       for(jj = Dsb->Get_start(i); jj < Dsb->Get_start(i+1); jj++)
       {
           if(Dsb->Get_index(jj) == i)
           {
               D_diag_offset[i] = Dsb->Get_offset(jj);
               break;
           }
           else
           {
               ostrstream ostr; 
               ostr << "wrong block structure of the diagmatrix";
               FAMGWrite(ostr);
               return(1);
           }
       }
    }



    while(viter(vi))
    {
        Dii = D.GetDiagValuePtr(vi);
        Aii = A.GetDiagValuePtr(vi);
        for(i = 0; i < Dsb->Get_ne(); i++) Dii[Dsb->Get_offset(i)] = 0.0;
        for(i = 0; i < Dsb->Get_nr(); i++)
        {
            sum = 0.0;
            for(jj = Asb->Get_start(i); jj < Asb->Get_start(i+1); jj++)
            {
                sum += Aii[Asb->Get_offset(jj)];
            }
            if(D_diag_offset[i] < 0) continue;
            Dii[D_diag_offset[i]] = 1.0/sum;
        } 
    }

 
    delete D_diag_offset;

    return 0;

}

// sum of the offdiagonals
int FAMGGrid::ConstructDiagonalSum()
{
    FAMGGridVector gv = GetGridVector();
    FAMGMatrixAlg &A = *GetMatrix();
    FAMGMatrixAlg &D = *GetDiagMatrix();
    FAMGVectorEntry vi;
    FAMGMatrixEntry dii, aij;
    FAMGVectorIter viter(gv);
    double *Dii, *Aij, d;
    short i, j, jj;

    const FAMGSparseBlock *Dsb = D.GetDiagSparseBlockPtr();
    const FAMGSparseBlock *Asb = A.GetSparseBlockPtr();

    short *D_diag_offset = new short[Dsb->Get_nr()]; 
    short *A_diag_offset = new short[Asb->Get_nr()];

    if(Dsb->Get_nr() != Asb->Get_nr())
    {
        ostrstream ostr; 
        ostr << "wrong block structure of the diagmatrix";
        FAMGWrite(ostr);
        return(1);
    }
    for(i = 0; i < Dsb->Get_nr(); i++) D_diag_offset[i] = -1;
    for(i = 0; i < Dsb->Get_nr(); i++)
    {
       for(jj = Dsb->Get_start(i); jj < Dsb->Get_start(i+1); jj++)
       {
           if(Dsb->Get_index(jj) == i)
           {
               D_diag_offset[i] = Dsb->Get_offset(jj);
               break;
           }
           else
           {
               ostrstream ostr; 
               ostr << "wrong block structure of the diagmatrix";
               FAMGWrite(ostr);
               return(1);
           }
       }
    }
    for(i = 0; i < Asb->Get_nr(); i++) A_diag_offset[i] = -1;
    for(i = 0; i < Asb->Get_nr(); i++)
    {
       for(jj = Asb->Get_start(i); jj < Asb->Get_start(i+1); jj++)
       {
           if(Asb->Get_index(jj) == i)
           {
               A_diag_offset[i] = Asb->Get_offset(jj);
               break;
           }
       }
    }


    for(i = 0; i < Dsb->Get_nr(); i++) 
    {
        if(D_diag_offset[i] < 0) 
        {
            if(A_diag_offset[i] >= 0)
            {
                ostrstream ostr; 
                ostr << "wrong block structure of the diagmatrix";
                FAMGWrite(ostr);
                return(1);
            }
            else
            {
                continue;
            }
        }
        for(j = i+1; j < Dsb->Get_nr(); j++)
        {
            if(D_diag_offset[i] == D_diag_offset[j])
            {
                if(A_diag_offset[i] != A_diag_offset[j])
                {
                    ostrstream ostr; 
                    ostr << "wrong block structure of the diagmatrix";
                    FAMGWrite(ostr);
                    return(1);
                }
                else
                {
                    D_diag_offset[j] = -1;
                    A_diag_offset[j] = -1;
                }
            }
        }
    }

    // sum up all off-diagonal entries
    while(viter(vi))
    {
        FAMGMatrixIter ai_iter(A,vi);
        FAMGMatrixIter di_iter(D,vi);
        di_iter(dii);
        Dii = D.GetValuePtr(dii);
        for(i = 0; i < Dsb->Get_ne(); i++) Dii[Dsb->Get_offset(i)] = 0.0;
        ai_iter(aij); // skip diagonal
        while(ai_iter(aij))
        {
            Aij = A.GetValuePtr(aij);
            for(i = 0; i < Dsb->Get_nr(); i++)
            {
                if(D_diag_offset[i] < 0) continue;
                Dii[D_diag_offset[i]] += Aij[A_diag_offset[i]];
            }   
        }
    }

    // compute -D^{-1}
    viter.reset();
    while(viter(vi))
    {
        FAMGMatrixIter di_iter(D,vi);
        di_iter(dii);
        Dii = D.GetValuePtr(dii);
        for(i = 0; i < Dsb->Get_nr(); i++)
        {
            if(D_diag_offset[i] < 0) continue;
            d = Dii[D_diag_offset[i]];
            if(Abs(d) > 1e-10) d = -1.0/d;
            else d = 1.0;
            Dii[D_diag_offset[i]] = d;
        }
    }



    delete D_diag_offset;
    delete A_diag_offset;

    return 0;

}

// lump the diagonal block
int FAMGGrid::ConstructDiagonalInvLump()
{
    FAMGGridVector gv = GetGridVector();
    FAMGMatrixAlg &A = *GetMatrix();
    FAMGMatrixAlg &D = *GetDiagMatrix();
    FAMGVectorEntry vi;
    FAMGVectorIter viter(gv);
    double *Dii, *Aii;
    short i, jj;

    const FAMGSparseBlock *Dsb = D.GetDiagSparseBlockPtr();
    const FAMGSparseBlock *Asb = A.GetDiagSparseBlockPtr();


    short *D_diag_offset = new short[Dsb->Get_nr()]; 
 
    if(Dsb->Get_nr() != Asb->Get_nr())
    {
        ostrstream ostr; 
        ostr << "wrong block structure of the diagmatrix";
        FAMGWrite(ostr);
        return(1);
    }
    for(i = 0; i < Dsb->Get_nr(); i++) D_diag_offset[i] = -1;
    for(i = 0; i < Dsb->Get_nr(); i++)
    {
       for(jj = Dsb->Get_start(i); jj < Dsb->Get_start(i+1); jj++)
       {
           if(Dsb->Get_index(jj) == i)
           {
               D_diag_offset[i] = Dsb->Get_offset(jj);
               break;
           }
           else
           {
               ostrstream ostr; 
               ostr << "wrong block structure of the diagmatrix";
               FAMGWrite(ostr);
               return(1);
           }
       }
    }



    if(Asb->Get_ne() != 4) assert(0);
    double *help = new double[4];
    double *helpinv = new double[4];
    double det;

    while(viter(vi))
    {
        Dii = D.GetDiagValuePtr(vi);
        Aii = A.GetDiagValuePtr(vi);
        for(i = 0; i < Dsb->Get_ne(); i++) Dii[Dsb->Get_offset(i)] = 0.0;
        for(i = 0; i < 4; i++)
        {
            help[i] = Aii[Asb->Get_offset(i)];
        }
        
        det = help[0]*help[3] - help[1]*help[2];
        helpinv[0] = help[3]/det;
        helpinv[1] = -help[1]/det;
        helpinv[2] = -help[2]/det;
        helpinv[3] = help[0]/det;

        Dii[D_diag_offset[0]] = helpinv[0];
        Dii[D_diag_offset[1]] = helpinv[3];
    }

 
    delete help;
    delete helpinv;

    delete D_diag_offset;

    return 0;

}


// use test vector for constructing the diagonal approximation
int FAMGGrid::ConstructDiagonalTV()
{
    FAMGGridVector gv = GetGridVector();
    FAMGMatrixAlg &A = *GetMatrix();
    FAMGMatrixAlg &D = *GetDiagMatrix();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
    FAMGVectorEntry vi;
    FAMGVectorIter viter(gv);
    double *Dii, *Aii, sum, *tA;
    short i, j, jj;

    const FAMGSparseBlock *Dsb = D.GetDiagSparseBlockPtr();
    const FAMGSparseBlock *Asb = A.GetDiagSparseBlockPtr();
    const FAMGSparseVector *tvAsv = tvA.GetSparseVectorPtr();
    const FAMGSparseVector *tvBsv = tvB.GetSparseVectorPtr();


    short *D_diag_offset = new short[Dsb->Get_nr()]; 
 
    if(Dsb->Get_nr() != Asb->Get_nr())
    {
        ostrstream ostr; 
        ostr << "wrong block structure of the diagmatrix";
        FAMGWrite(ostr);
        return(1);
    }
    for(i = 0; i < Dsb->Get_nr(); i++) D_diag_offset[i] = -1;
    for(i = 0; i < Dsb->Get_nr(); i++)
    {
       for(jj = Dsb->Get_start(i); jj < Dsb->Get_start(i+1); jj++)
       {
           if(Dsb->Get_index(jj) == i)
           {
               D_diag_offset[i] = Dsb->Get_offset(jj);
               break;
           }
           else
           {
               ostrstream ostr; 
               ostr << "wrong block structure of the diagmatrix";
               FAMGWrite(ostr);
               return(1);
           }
       }
    }



    while(viter(vi))
    {
        Dii = D.GetDiagValuePtr(vi);
        Aii = A.GetDiagValuePtr(vi);
        tA = tvA.GetValuePtr(vi);
        for(i = 0; i < Dsb->Get_ne(); i++) Dii[Dsb->Get_offset(i)] = 0.0;
        for(i = 0; i < Dsb->Get_nr(); i++)
        {
            sum = 0.0;
            for(jj = Asb->Get_start(i); jj < Asb->Get_start(i+1); jj++)
            {
                j = Asb->Get_index(jj);
                sum += Aii[Asb->Get_offset(jj)]*tA[tvAsv->Get_comp(j)];
;
            }
            if(D_diag_offset[i] < 0) continue;
            if(Abs(sum) > 1e-10*tA[tvAsv->Get_comp(i)]) 
            {
                Dii[D_diag_offset[i]] = tA[tvAsv->Get_comp(i)]/sum;
            }
            else
            {
                Dii[D_diag_offset[i]] = 1e+10;
            }
        } 

    }

 
    delete D_diag_offset;

    return 0;

}

// use test vector for constructing the diagonal approximation
int FAMGGrid::ConstructDiagonalInvTV()
{
    FAMGGridVector gv = GetGridVector();
    FAMGMatrixAlg &A = *GetMatrix();
    FAMGMatrixAlg &D = *GetDiagMatrix();
	const FAMGVector &tvA = *vector[FAMGTVA];
	const FAMGVector &tvB = *vector[FAMGTVB];
    FAMGVectorEntry vi;
    FAMGVectorIter viter(gv);
    double *Dii, *Aii, *tA;
    short i, jj;

    const FAMGSparseBlock *Dsb = D.GetDiagSparseBlockPtr();
    const FAMGSparseBlock *Asb = A.GetDiagSparseBlockPtr();
    const FAMGSparseVector *tvAsv = tvA.GetSparseVectorPtr();
    const FAMGSparseVector *tvBsv = tvB.GetSparseVectorPtr();


    short *D_diag_offset = new short[Dsb->Get_nr()]; 
 
    if(Dsb->Get_nr() != Asb->Get_nr())
    {
        ostrstream ostr; 
        ostr << "wrong block structure of the diagmatrix";
        FAMGWrite(ostr);
        return(1);
    }
    for(i = 0; i < Dsb->Get_nr(); i++) D_diag_offset[i] = -1;
    for(i = 0; i < Dsb->Get_nr(); i++)
    {
       for(jj = Dsb->Get_start(i); jj < Dsb->Get_start(i+1); jj++)
       {
           if(Dsb->Get_index(jj) == i)
           {
               D_diag_offset[i] = Dsb->Get_offset(jj);
               break;
           }
           else
           {
               ostrstream ostr; 
               ostr << "wrong block structure of the diagmatrix";
               FAMGWrite(ostr);
               return(1);
           }
       }
    }

    short nr = Asb->Get_nr();
    double *help = new double[nr];
    double *decomp = new double[nr*nr]; 
    short *pivotmap = new short[nr]; 

    while(viter(vi))
    {
        tA = tvA.GetValuePtr(vi)+tvAsv->Get_comp(0);
        Aii = A.GetDiagValuePtr(vi);
        Dii = D.GetDiagValuePtr(vi);
        SparseBlockMCopyDense(decomp,Asb,Aii);
        if(LR_Decomp(nr,decomp,pivotmap)) assert(0);
        if(LR_Solve(nr,decomp,pivotmap,help,tA)) assert(0);

        for(i = 0; i < Dsb->Get_ne(); i++) Dii[Dsb->Get_offset(i)] = 0.0;
        for(i = 0; i < Dsb->Get_nr(); i++)
        {
            if(tA[i] > 1e-10*Abs(help[i])) 
            {
                Dii[D_diag_offset[i]] = help[i]/tA[i];
            }
            else
            {
                for(jj = Asb->Get_start(i); jj < Asb->Get_start(i+1); jj++)
                {
                    if (Asb->Get_index(jj) == i)
                    {
                        Dii[D_diag_offset[i]] = 1.0/Aii[Asb->Get_offset(jj)];
                        break;
                    }
                }
            }
        } 

    }

 
    delete D_diag_offset;
    delete help;
    delete pivotmap;
    delete decomp;

    return 0;

}

// use test vector for constructing the diagonal approximation
int FAMGGrid::ConstructDiagonalInverse()
{
    FAMGGridVector gv = GetGridVector();
    FAMGMatrixAlg &A = *GetMatrix();
    FAMGMatrixAlg &D = *GetDiagMatrix();
    FAMGVectorEntry vi;
    FAMGVectorIter viter(gv);
    double *Dii, *Aii;

    const FAMGSparseBlock *Dsb = D.GetDiagSparseBlockPtr();
    const FAMGSparseBlock *Asb = A.GetDiagSparseBlockPtr();

    short nr = Dsb->Get_nr();
    if(Asb->Get_nr() != nr) assert(0);
    if(Dsb->Get_ne() != nr*nr) assert(0);

    if(Dsb->Get_offset(0) + Dsb->Get_ne() != Dsb->Get_maxoffset()+1) assert(0);

    while(viter(vi))
    {
        Aii = A.GetDiagValuePtr(vi);
        Dii = D.GetDiagValuePtr(vi)+Dsb->Get_offset(0);
        SparseBlockMCopyDense(Dii,Asb,Aii);
        if(LR_Decomp(nr,Dii)) assert(0);

    }

 
    return 0;

}



#endif 

