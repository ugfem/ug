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

#ifdef FAMG_SPARSE_BLOCK



void FAMGSparseBlock::Transposed(const FAMGSparseBlock &sb)
{
    int i, j, total, delta, off, ii;
    
    ne = sb.ne;
    nid = sb.nid;
    nr = sb.nc;
    nc = sb.nr;

    start = new int[nr+1];
    index = new int[ne];
    offset = new int[ne];

    for(i = 0; i <= nr; i++) start[i] = 0;
    for(j = 0; j < sb.ne; j++)
    {
        i = sb.index[j];
        start[i]++;
    }
    total = 0;
    for(i = 0; i <= nr; i++)
    {
        delta = start[i];
        start[i] = total;
        total += delta;
    }

    int *end = new int[nr];
    for(i = 0; i < nr; i++) end[i] = start[i];
    for(i = 0; i < sb.nr; i++)
    {
        for(ii = sb.start[i]; ii < sb.start[i+1]; ii++)
        {
            j = sb.index[ii];
            off = sb.offset[ii];
            index[end[j]] = i;
            offset[end[j]] = off;
            end[j]++;
        }
    }

    delete end;

    return;
}

int FAMGSparseBlock::Product(const FAMGSparseBlock &sb)
{
    int i, j, k, ii, jj, kk, oij, ojk, off, equal;

    if(sb.nr != sb.nc) return 1;
    nr = nc = sb.nr;

    // construct start first
    start = new int[nr+1];
    int *help = new int[nc];

    start[0] = 0;
    for(i = 0; i < nr; i++)
    {
        memset((void *)help, 0, nc*sizeof(int));
        for(ii = sb.start[i]; ii < sb.start[i+1]; ii++)
        {
            j = sb.index[ii];
            for(jj = sb.start[j]; jj < sb.start[j+1]; jj++)
            {
                k = sb.index[jj];
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
    index = new int[ne];
    for(i = 0; i < nr; i++)
    {
        memset((void *)help, 0, nc*sizeof(int));
        for(ii = sb.start[i]; ii < sb.start[i+1]; ii++)
        {
            j = sb.index[ii];
            for(jj = sb.start[j]; jj < sb.start[j+1]; jj++)
            {
                k = sb.index[jj];
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


    // construct offset


    offset = new int[ne];
    if ((sb.nid > 5) || (sb.nid == sb.ne)) // no identification
    {     
        for(i = 0; i < ne; i++) offset[i] = i;
        
        delete help;
        return 0;
    }

    // look for equal entries

    int *helpmat = new int[ne*sb.nid*sb.nid];
    memset((void*)helpmat, 0, ne*sb.nid*sb.nid*sizeof(int));

    for(i = 0; i < nr; i++)
    {
        for(ii = start[i]; ii < start[i+1]; ii++)
        {
            k = index[ii];
            help[k] = ii;
        }
        for(ii = sb.start[i]; ii < sb.start[i+1]; ii++)
        {
            j = sb.index[ii];
            oij = sb.offset[ii];
            for(jj = sb.start[j]; jj < sb.start[j+1]; jj++)
            {
                k = sb.index[jj];
                ojk = sb.offset[jj];
                kk = help[k];
                helpmat[sb.nid*sb.nid*kk+oij*sb.nid+ojk]++;
            }
        }
    }
    
    int *hm = helpmat;
    off= 0; ii = 0;
    while(hm != NULL)
    {  
        offset[ii] = off;
        for(i = ii+1; i < ne; i++)
        {
            // compare help matrices
            equal = 1;
            for(j = 0; j < sb.nid; j++)
            {
                for(k = 0; k < sb.nid; k++)
                {
                    if(hm[j*sb.nid+k] != helpmat[sb.nid*sb.nid*i+j*sb.nid+k])
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
                helpmat[sb.nid*sb.nid*i] = -1;
            }
        }
        off++;
            
        hm = NULL;
        for(j = ii+1; j < ne; j++)
        {
            if (helpmat[sb.nid*sb.nid*j] > -1)
            {
                hm = &helpmat[sb.nid*sb.nid*j];
                ii = j;
                break;
            }
        }
    }
        
    nid = off; 

    delete help;
    delete helpmat;
 
    return 0;
}
 
int FAMGSparseBlock::CheckCGIdent(const FAMGSparseTransfer &sp, const FAMGSparseTransfer &sr)
{
    int i, j, k, jj, kk, oj, oi, ok;

    int pnr = sp.Get_nr();
    int *poffset = sp.Get_offset();

    int rnr = sr.Get_nr();
    int *roffset = sr.Get_offset();

    if ((nr != rnr) || (nc != pnr))
    {
        ostrstream ostr; 
        ostr << "wrong prolongation/restriction dimensions";
        FAMGWrite(ostr);
        return(1);
    }
        

    // check if prolongation preservs identification structure
    for(jj = 0; jj < ne; jj++)
    {
        j = index[jj];
        oj = offset[jj];
        for(kk = jj+1; kk < ne; kk++)
        {
            k = index[kk];
            ok = offset[kk];
            if((oj == ok) && (poffset[j] != poffset[k]))
            {
                    ostrstream ostr; 
                    ostr << "different CG identification structure: not yet implemented.";
                    FAMGWrite(ostr);
                    return(1);
            }
        }
    }

    // check if restriction preservs identification structure

    for(i = 0; i < nr; i++)
    {
        for(jj = start[i]; jj < start[i+1]; jj++)
        {
            oi = offset[jj];
            for(j = i+1; j < nr; j++)
            {
                for(kk = start[j]; kk < start[j+1]; kk++)
                {
                    oj = offset[kk];
                    if((oi == oj) && (roffset[i] != roffset[j]))
                    {
                        ostrstream ostr; 
                        ostr << "different CG identification structure: not yet implemented.";
                        FAMGWrite(ostr);
                        return(1);
                    }
                }
            }
        }
    }

    return 0;
}
         
                   

    

void FAMGSparseTransfer::Construct(const FAMGSparseBlock &sb)
{
    int i, j, ii, equal, off, sbnid;

    sbnid = sb.Get_nid();
    nr = sb.Get_nr();
    offset = new int[nr];

    int *help = new int[nr*sbnid];    
    memset((void*) help, 0, nr*sbnid*sizeof(int));


    for(i = 0; i < nr; i++)
    {
        for(ii = sb.Get_start(i); ii < sb.Get_start(i+1); ii++)
        {
            help[i*sbnid+sb.Get_offset(ii)]++;
        }
    }
            
    int *hv = help;
    off= 0; ii = 0;
    while(hv != NULL)
    {  
        offset[ii] = off;
        for(i = ii+1; i < nr; i++)
        {
            // compare help matrices
            equal = 1;
            for(j = 0; j < sbnid; j++)
            {
                if(hv[j] != help[sbnid*i+j])
                {
                    equal = 0;
                    break;
                }
            }
            if(equal) 
            {
                offset[i] = off;
                help[sbnid*i] = -1;
            }
        }
        off++;
            
        hv = NULL;
        for(j = ii+1; j < nr; j++)
        {
            if (help[sbnid*j] > -1)
            {
                hv = &help[sbnid*j];
                ii = j;
                break;
            }
        }
    }
        

    delete help;
    
    return;
}


void SparseBlockProduct()
{
    int i, j, k, ii, jj, oj, compute;
    double aij, bjk;
    FAMGSparseBlock sb, sp;

    int *ho = new int[sb.nid];
    int *hv = new int[sb.nc];

    for(i = 0; i < sp.ne; i++) /* C */(sp.offset[i]) = 0.0;
    for(i = 0; i < sp.nid; i++) ho[i] = 0;
    for(j = 0; j < sp.nc; j++) hv[j] = -1;
    for(i = 0; i < sp.nr; i++)
    {
        compute = 0;
        for(ii = sp.start[i]; ii < sp.start[i+1]; ii++)
        {
            j = sp.index[ii];
            oj = sp.offset[ii];
            if(!ho[oj]) // not yet computed
            {
                ho[oj] = 1;
                hv[j] = ii;
                compute = 1;
            }
        }
        if(!compute) continue;
        for(ii = sb.start[i]; ii < sb.start[i+1]; ii++)
        {
            j = sb.index[ii];
            aij = /* A */ (sb.offset[ii]);
            for(jj = sb.start[j]; jj < sb.start[j+1]; jj++)
            {
                k = sb.index[jj];
                bjk = /* B */ (sb.offset[jj]);
                if(hv[k] > -1) /* C */ (sp.offset[hv[k]]) += aij*bjk;
            }
        }

        for(j = 0; j < sb.nc; j++) hv[j] = -1;
    }

    delete ho;
    delete hv;

    return;
}

        

void FAMGTestSparseBlock()
{
    FAMGSparseBlock sb, sbt, sbp;
    FAMGSparseTransfer st;
   
    sb.ne = 5;
    sb.nid = 3;    
    sb.nc = 3;
    sb.nr = 3;

    sb.start = new int[sb.nr+1];
    sb.index = new int[sb.ne];
    sb.offset = new int[sb.ne];

    sb.start[0] = 0; sb.start[1] = 2; sb.start[2] = 4; sb.start[3] = 5;
    sb.index[0] = 0; sb.index[1] = 2; sb.index[2] = 0; sb.index[3] = 1;
    sb.index[4] = 2;
    sb.offset[0] = 0; sb.offset[1] = 1; sb.offset[2] = 1; sb.offset[3] = 0;
    sb.offset[4] = 2;

    sbt.Transposed(sb);
    sbp.Product(sb);
    st.Construct(sb);
    sb.CheckCGIdent(st,st);

    return ;
}

#endif 
