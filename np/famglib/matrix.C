/****************************************************************************/
/*																			*/
/* File:      matrix.C														*/
/*																			*/
/* Purpose:   cmg matrix classes functions									*/
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
#include <math.h>
#include "misc.h"
#include "matrix.h"
#include "heap.h"
#include "grid.h"

/* RCS_ID
$Header$
*/

static double cmgzero = 0.0;


void CMGMatrix::DevideFGDefect(double *unknown, double *defect)
{
    CMGMatrixPtr matij;
    double *u = unknown;
    double *d = defect;

    for(int i = 0; i < n; i++)
    {
        matij = GetStart(i);
        if(matij.GetType()) // FG node
        {
            (*d) = (*d)/matij.GetData();
            (*u) += (*d);
        }
        d++; u++;
    }

    return;
}
 
void CMGMatrix::VecMinusMatVec(double *defect, double *rhs, double *unknown)
{
    CMGMatrixPtr matij;
    double sum;
    double *d = defect;
    double *r = rhs;
    double *u = unknown;
    
    for(int i = 0; i < n; i++)
    {
        matij = GetStart(i);
        sum = matij.GetData()*(*u);
        while(matij.GetNext())
        {
            sum += matij.GetData()*unknown[matij.GetIndex()];
        }
        (*d) = (*r) - sum;
        u++; d++; r++;
    }

    return;
}

void CMGMatrix::JAC(double *vec)
{
    const double omega = 0.6666666;
    double *v = vec;
    for(int i = 0; i < n; i++)
    {
        (*v) = omega*(*v)/GetDiag(i);
        v++;
    }
    return;
}

void CMGMatrix::FGS(double *vec)
{
    CMGMatrixPtr matij;
    double sum, mii;
    int i,j;

    /* forward Gauss-Seidel */

    double *v = vec;
    for(i = 0; i < n; i++)
    {
        sum = *v;
        matij = GetStart(i);
        mii = matij.GetData();
        while(matij.GetNext())
        {
            j = matij.GetIndex();
            if (j < i)
            {
                sum -= matij.GetData()*vec[j];
            }
        }
        (*v) = sum/mii;
        v++;
    }
    return;
}

void CMGMatrix::BGS(double *vec)
{
    CMGMatrixPtr matij;
    double sum, mii;
    int i,j;

    /* backward Gauss-Seidel */

    double *v = vec+(n-1);
    for(i = n-1; i >= 0; i--)
    {
        sum = *v;
        matij = GetStart(i);
        mii = matij.GetData();
        while(matij.GetNext())
        {
            j = matij.GetIndex();
            if (j > i)
            {
                sum -= matij.GetData()*vec[j];
            }
        }
        (*v) = sum/mii;
        v--;
    }
    return;
}

void CMGMatrix::SGS(double *vec)
{
    CMGMatrixPtr matij;
    double sum, mii;
    int i,j;

    /* symmetric Gauss-Seidel */

    double *v = vec;
    for(i = 0; i < n; i++)
    {
        sum = *v;
        matij = GetStart(i);
        mii = matij.GetData();
        while(matij.GetNext())
        {
            j = matij.GetIndex();
            if (j < i)
            {
                sum -= matij.GetData()*vec[j];
            }
        }
        (*v) = sum/mii;
        v++;
    }

    v = vec;
    for(i = 0; i < n; i++)
    {
        *v = (*v)*GetDiag(i);
    }

    v = vec+(n-1);
    for(i = n-1; i >= 0; i--)
    {
        sum = *v;
        matij = GetStart(i);
        mii = matij.GetData();
        while(matij.GetNext())
        {
            j = matij.GetIndex();
            if (j > i)
            {
                sum -= matij.GetData()*vec[j];
            }
        }
        (*v) = sum/mii;
        v--;
    }
    return;
}

int CMGMatrix::GetSmallestIndex()
{
    int j,i,min,nc;

    nc = start[1]; // smallest index is supposed to be in the first row
    min = index[0];
    for(i = 1; i < nc; i++)
    {
        j = index[i];
        if(j < min) min = j;
    }

    return min;
}


void CMGMatrix::ModifyIndex(int f)
{
    CMGIndexBitField bf;
    int *ind;
    
    ind = index;
    for(int i = 0; i < nl; i++)
    {
        bf.id  = (*(ind) - f);
        bf.type = 0;
        (*ind) = *((int*) &bf);
        ind++;
    }

    return;
}
        
void CMGMatrix::RemodifyIndex(int f)
{
    CMGIndexBitField *bf;
    int *ind;
    
    ind = index;
    for(int i = 0; i < nl; i++)
    {
        bf  = (CMGIndexBitField *)ind;
        (*ind) = (bf->id+f);
        ind++;
    }

    return;
}
        
void CMGMatrix::ModifyIndex(int *type, int f)
{
    CMGIndexBitField bf;
    int *ind, ii;
    
    ind = index; 
    for(int i = 0; i < nl; i++)
    {
        ii = *(ind)-f;
        bf.id  = ii;
        bf.type = type[ii];
        (*ind) = *((int*) &bf);
        ind++; 
    }

    return;
}
        
void CMGMatrix::RemodifyIndex(int *type, int f)
{
    CMGIndexBitField *bf;
    int *ind, ii;
    
    ind = index; 
    for(int i = 0; i < nl; i++)
    {
        bf  = (CMGIndexBitField *)ind;
        ii =  bf->id;
        type[ii] = bf->type;
        (*ind) = ii+f;
        ind++; 
    }

    return;
}
        
int CMGMatrix::OrderColumns(int *map)
{
    double hd, *mat, *firstentry;
    int i, j, nc;
    CMGIndexBitField *firstindex, *ind, hi;

    for(i = 0; i < n; i++)
    {
        nc = start[i+1]-start[i];
        firstindex = ind =  (CMGIndexBitField *) index+start[i];
        firstentry = mat = entry+start[i];
        for(j = 0; j < nc; j++)
        {
            if(ind->id == i) break;
            ind++; mat++;
        }
        if(j == nc)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": wrong matrix structure." << endl;
            CMGError(ostr);
            return(1);
        }
        map[i] = j;
        if(j == 0) continue;
        hi = (*firstindex);
        (*firstindex) = (*ind);
        (*ind) = hi;
        hd = (*firstentry);
        (*firstentry) = (*mat);
        (*mat) = hd;
    }

    return 0;
}

int CMGMatrix::OrderColumns()
{
    double hd, *mat, *firstentry;
    int i, j, nc;
    CMGIndexBitField *firstindex, *ind, hi;

    for(i = 0; i < n; i++)
    {
        nc = start[i+1]-start[i];
        firstindex = ind =  (CMGIndexBitField *) index+start[i];
        firstentry = mat = entry+start[i];
        for(j = 0; j < nc; j++)
        {
            if(ind->id == i) break;
            ind++; mat++;
        }
        if(j == nc)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": wrong matrix structure." << endl;
            CMGError(ostr);
            return(1);
        }
        if(j == 0) continue;
        hi = (*firstindex);
        (*firstindex) = (*ind);
        (*ind) = hi;
        hd = (*firstentry);
        (*firstentry) = (*mat);
        (*mat) = hd;
    }

    return 0;
}

int CMGMatrix::ReorderColumns(int *map)
{
    double hd, *mat, *firstentry;
    int i;
    CMGIndexBitField *firstindex, *ind, hi;

    for(i = 0; i < n; i++)
    {
        if(map[i] == 0) continue;
        firstindex = (CMGIndexBitField *) index+start[i];
        firstentry = entry+start[i];
        ind = firstindex + map[i];
        mat = firstentry + map[i];

        hi = (*firstindex);
        (*firstindex) = (*ind);
        (*ind) = hi;
        hd = (*firstentry);
        (*firstentry) = (*mat);
        (*mat) = hd;
    }

    return 0;
}

int CMGMatrix::OrderColumns2(CMGGraph *graph)
{
    CMGNode *node;
    double hd, *mat, *firstentry;
    int i, j, nc;
    CMGIndexBitField *firstindex, *ind, hi;

    node = graph->GetNode();

    for(i = 0; i < n; i++)
    {
        if((node+i)->IsFGNode()) continue;
        nc = start[i+1]-start[i];
        firstindex = ind = (CMGIndexBitField *) index+start[i];
        firstentry = mat = entry+start[i];
        for(j = 0; j < nc; j++)
        {
            if(ind->id == i) break;
            ind++; mat++;
        }
        if(j == nc)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": wrong matrix structure." << endl;
            CMGError(ostr);
            return(1);
        }
        if(j == 0) continue;
        hi = (*firstindex);
        (*firstindex) = (*ind);
        (*ind) = hi;
        hd = (*firstentry);
        (*firstentry) = (*mat);
        (*mat) = hd;
    }

    return 0;
}


int CMGMatrix::ConstructAdjoined()
{
    CMGMatrixPtr matij, matjk;
    int i, j, k, found;

    adjoined = (double **) CMGGetMem(nl*sizeof(double *), CMG_FROM_TOP);
    if (adjoined == NULL) return 1;

    for(i = 0; i < n; i++)
    {
        matij = GetStart(i);
        do
        {
            j = matij.GetIndex();
            if (j == i)
            {
                *(matij.GetAdjoinedPtr()) = matij.GetEntry();
            }
            else 
            {
                matjk = GetStart(j);
                found = 0;
                do              
                {
                    k = matjk.GetIndex();
                    if (k == i)
                    {
                        *(matij.GetAdjoinedPtr()) = matjk.GetEntry();
                        found = 1;
                        break;
                    }
                } while(matjk.GetNext());
                if(!found)
                {
                    *(matij.GetAdjoinedPtr()) = &cmgzero;
                }
            }
        } while(matij.GetNext());

    }

    return 0;
}

int CMGMatrix::ConstructAdjoinedB()
{
    CMGMatrixPtr matij, matjk;
    int i, j, k, found;

    adjoined = (double **) CMGGetMem(nl*sizeof(double *), CMG_FROM_BOTTOM);
    if (adjoined == NULL) return 1;

    for(i = 0; i < n; i++)
    {
        matij = GetStart(i);
        do
        {
            j = matij.GetIndex();
            if (j == i)
            {
                *(matij.GetAdjoinedPtr()) = matij.GetEntry();
            }
            else 
            {
                matjk = GetStart(j);
                found = 0;
                do              
                {
                    k = matjk.GetIndex();
                    if (k == i)
                    {
                        *(matij.GetAdjoinedPtr()) = matjk.GetEntry();
                        found = 1;
                        break;
                    }
                } while(matjk.GetNext());
                if(!found)
                {
                    *(matij.GetAdjoinedPtr()) = &cmgzero;
                }
            }
        } while(matij.GetNext());

    }

    return 0;
}


int CMGMatrix::ConstructAdjoined2(CMGGraph *graph)
{
    CMGMatrixPtr matij, matjk;
    int i, j, k, found;
    CMGNode *node = graph->GetNode();

    adjoined = (double **) CMGGetMem(nl*sizeof(double *), CMG_FROM_BOTTOM);
    if (adjoined == NULL) return 1;

    for(i = 0; i < n; i++)
    {
       if((node+i)->IsFGNode()) continue;
       matij = GetStart(i);
        do
        {
            j = matij.GetIndex();
            if (j == i)
            {
                *(matij.GetAdjoinedPtr()) = matij.GetEntry();
            }
            else
            {
                matjk = GetStart(j);
                found = 0;
                do              
                {
                    k = matjk.GetIndex();
                    if (k == i)
                    {
                        *(matij.GetAdjoinedPtr()) = matjk.GetEntry();
                        found = 1;
                        break;
                    }
                } while(matjk.GetNext());
                if(!found)
                {
                    *(matij.GetAdjoinedPtr()) = &cmgzero;
                    // return(1);
                }
            }
        } while(matij.GetNext());

    }

    return 0;
}


int CMGMatrix::ConstructEnd()
{
    int i;

    end = (int *) CMGGetMem(sizeof(int)*n,CMG_FROM_TOP);
    if(end == NULL) return 1;
     
    for(i = 0; i < n; i++) end[i] = start[i+1];

    return 0;
}
    
int CMGMatrix::ConstructEndB()
{
    int i;

    end = (int *) CMGGetMem(sizeof(int)*n,CMG_FROM_BOTTOM);
    if(end == NULL) return 1;
     
    for(i = 0; i < n; i++) end[i] = start[i+1];

    return 0;
}
    
int CMGMatrix::ConstructEnd2()
{
    int i;

    end = (int *) CMGGetMem(sizeof(int)*n,CMG_FROM_BOTTOM);
    if(end == NULL) return 1;
     
    for(i = 0; i < n; i++) end[i] = start[i+1];

    return 0;
}
    
int CMGMatrix::CGMatrix(CMGMatrix *fgmatrix, CMGTransfer *transfer, int *father)
{    
    CMGMatrixPtr matik, matjk;
    CMGTransferEntry *trans, *transij, *transks;
    double mis, tij, mjk, tks, mik, *firstentry, *mat;
    int i, j, s, k, z, nsaved, ii, offset, fnn, *fathermap, *tmpflag, *firstindex, *ind, found;
    
    fnn = fgmatrix->GetN();
    trans = transfer->GetRow();

    CMGMarkHeap(CMG_FROM_BOTTOM);

    fathermap = (int *) CMGGetMem(sizeof(int)*fnn,CMG_FROM_BOTTOM);
    if (fathermap == NULL) {CMGReleaseHeap(CMG_FROM_BOTTOM);  return 1;}

    tmpflag = (int *) CMGGetMem(sizeof(int)*n,CMG_FROM_BOTTOM);
    if (tmpflag == NULL) {CMGReleaseHeap(CMG_FROM_BOTTOM);  return 1;}

    for(z = 0; z < n; z++)
    {
        fathermap[father[z]] = z;
        tmpflag[z] = 0;
    }

    // count matrix entries  
    start[0] = 0;
    offset = 0;
    for(z = 0; z < n; z++)
    {
        i = father[z];
        matik = fgmatrix->GetStart(i);
        do
        {
            k = matik.GetIndex();
            if (matik.GetType()) // FGnode
            {
                 for(transks = (trans+k)->GetNext(); transks != NULL; transks = transks->GetNext())
                {
                    s = fathermap[transks->GetId()];
                    if(tmpflag[s] == 0) 
                    {
                        tmpflag[s] = 1;
                        offset++;
                    }     
                } 
            }
            else
            {
                s = fathermap[k];
                if(tmpflag[s] == 0) 
                {
                    tmpflag[s] = 1;
                    offset++;
                }     
            }      
        } while(matik.GetNext());
        for(transij = (trans+i)->GetNext(); transij != NULL; transij = transij->GetNext())
        {
            j = transij->GetId();
            matjk = fgmatrix->GetStart(j);
            do
            {
                k = matjk.GetIndex();
                if(matjk.GetType())
                {
                    for(transks = (trans+k)->GetNext(); transks != NULL; transks = transks->GetNext())
                    {
                        s = fathermap[transks->GetId()];
                        if(tmpflag[s] == 0) 
                        {
                            tmpflag[s] = 1;
                            offset++;
                        }     
                    } 
                }
                else
                {
                    s = fathermap[k];
                    if(tmpflag[s] == 0) 
                    {
                        tmpflag[s] = 1;
                        offset++;
                    }     
                }      
            } while(matjk.GetNext());
        }
        start[z+1] = offset;

        // reset flags
        matik = fgmatrix->GetStart(i);
        do
        {
            k = matik.GetIndex();
            if(matik.GetType())
            {
                for(transks = (trans+k)->GetNext(); transks != NULL; transks = transks->GetNext())
                {
                    s = fathermap[transks->GetId()];
                    tmpflag[s] = 0; 
                } 
            }
            else
            {
                s = fathermap[k];
                tmpflag[s] = 0; 
            }      
        } while(matik.GetNext());
        for(transij = (trans+i)->GetNext(); transij != NULL; transij = transij->GetNext())
        {
            j = transij->GetId();
            matjk = fgmatrix->GetStart(j);
            do
            {
                k = matjk.GetIndex();
                if(matjk.GetType())
                {
                    for(transks = (trans+k)->GetNext(); transks != NULL; transks = transks->GetNext())
                    {
                        s = fathermap[transks->GetId()];
                        tmpflag[s] = 0; 
                    } 
                }
                else
                {
                    s = fathermap[k];
                    tmpflag[s] = 0; 
                }      
            } while(matjk.GetNext());
        }
    }

    nl = offset;
    index = (int *) CMGGetMem(nl*sizeof(int), CMG_FROM_TOP);
    if (index == NULL) return 1;
    
    entry = (double *) CMGGetMem(nl*sizeof(double), CMG_FROM_TOP);
    if (entry == NULL) return 1;
    
    // save entries
    for(z = 0; z < n; z++)
    {
        i = father[z];
        nsaved = 0;
        firstindex = index+start[z];
        firstentry = entry+start[z];
        matik = fgmatrix->GetStart(i);
        do
        {
            k = matik.GetIndex();
            mik = matik.GetData();
            if(matik.GetType())
            {
                for(transks = (trans+k)->GetNext(); transks != NULL; transks = transks->GetNext())
                {
                    s = fathermap[transks->GetId()];
                    tks = transks->GetData();
                    mis = mik*tks;
                    // does mis already exist ?
                    ind = firstindex;
                    mat = firstentry;
                    found = 0;
                    for(ii = 0; ii < nsaved; ii++)
                    {
                        if( (*ind) == s) {(*mat) += mis; found = 1; break;}
                        ind++; mat++;
                    }
                    if(!found) // no
                    {
                        (*ind) = s;
                        (*mat) = mis;
                        nsaved++;
                    }
                        
                 } 
            }
            else
            {
                mis = mik;
                s = fathermap[k];
                // does mis already exist ?
                ind = firstindex;
                mat = firstentry;
                found = 0;
                for(ii = 0; ii < nsaved; ii++)
                {
                    if( (*ind) == s) {(*mat) += mis; found = 1; break;}
                    ind++; mat++;
                }
                if(!found) // no
                {
                    (*ind) = s;
                    (*mat) = mis;
                    nsaved++;
                }
             }      
        } while(matik.GetNext());
        for(transij = (trans+i)->GetNext(); transij != NULL; transij = transij->GetNext())
        {
            j = transij->GetId();
            tij = transij->GetData();
            matjk = fgmatrix->GetStart(j);
            do
            {
                k = matjk.GetIndex();
                mjk = matjk.GetData();
                if(matjk.GetType())
                {
                    for(transks = (trans+k)->GetNext(); transks != NULL; transks = transks->GetNext())
                    {
                        s = fathermap[transks->GetId()];
                        tks = transks->GetData();
                        mis = tij*mjk*tks;
                        // does mis already exist ?
                        ind = firstindex;
                        mat = firstentry;
                        found = 0;
                        for(ii = 0; ii < nsaved; ii++)
                        {
                            if( (*ind) == s) {(*mat) += mis; found = 1; break;}
                            ind++; mat++;
                        }
                        if(!found) // no
                        {
                            (*ind) = s;
                            (*mat) = mis;
                            nsaved++;
                        }
                    } 
                }
                else
                {
                    s = fathermap[k];
                    mis = tij*mjk;
                    // does mis already exist ?
                    ind = firstindex;
                    mat = firstentry;
                    found = 0;
                    for(ii = 0; ii < nsaved; ii++)
                    {
                        if( (*ind) == s) {(*mat) += mis; found = 1; break;}
                        ind++; mat++;
                    }
                    if(!found) // no
                    {
                        (*ind) = s;
                        (*mat) = mis;
                        nsaved++;
                    }
                }      
            } while(matjk.GetNext());
        }
        if((start[z]+nsaved) != start[z+1])
        {  
            ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
            CMGError(ostr);
            return(1);
        }
    }

    ModifyIndex(0);
    if(ConstructEnd()) return 1;
    if(OrderColumns()) return 1;
    if(ConstructAdjoined()) return 1;

    CMGReleaseHeap(CMG_FROM_BOTTOM);

    return 0;
}


int CMGMatrix::TmpMatrix(CMGMatrix *matrix, CMGTransfer *transfer, CMGGraph *graph)
{    
    CMGMatrixPtr matik, matjk;
    CMGTransferEntry *trans, *transij, *transks;
    CMGNode *node;
    double mis, tij, mjk, tks, mik, *firstentry, *mat;
    int i, j, s, k, nsaved, ii, offset, *firstindex, *ind, found;
   
    trans = transfer->GetRow();
    node = graph->GetNode();

    // count matrix entries  
    start[0] = 0;
    offset = 0;
    for(i = 0; i < n; i++)
    {
        if((node+i)->IsFGNode()) 
        {
            start[i+1] = offset;
            continue;
        }
        matik = matrix->GetStart(i);
        do
        {
            k = matik.GetIndex();
            if((node+k)->IsFGNode())
            {
                for(transks = (trans+k)->GetNext(); transks != NULL; transks = transks->GetNext())
                {
                    s = transks->GetId();
                    if((node+s)->GetFlag() == 0) 
                    {
                        (node+s)->SetFlag(1);
                        offset++;
                    }     
                } 
            }
            else
            {
                s = k;
                if((node+s)->GetFlag() == 0) 
                {
                    (node+s)->SetFlag(1);
                    offset++;
                }     
            }      
        } while(matik.GetNext());
        for(transij = (trans+i)->GetNext(); transij != NULL; transij = transij->GetNext())
        {
            j = transij->GetId();
            matjk = matrix->GetStart(j);
            do
            {
                k = matjk.GetIndex();
                if((node+k)->IsFGNode())
                {
                    for(transks = (trans+k)->GetNext(); transks != NULL; transks = transks->GetNext())
                    {
                        s = transks->GetId();
                        if((node+s)->GetFlag() == 0) 
                        {
                            (node+s)->SetFlag(1);
                            offset++;
                        }     
                    } 
                }
                else
                {
                    s = k;
                    if((node+s)->GetFlag() == 0) 
                    {
                        (node+s)->SetFlag(1);
                        offset++;
                    }     
                }      
            } while(matjk.GetNext());
        }
        start[i+1] = offset;

        // reset flags
        matik = matrix->GetStart(i);
        do
        {
            k = matik.GetIndex();
            if((node+k)->IsFGNode())
            {
                for(transks = (trans+k)->GetNext(); transks != NULL; transks = transks->GetNext())
                {
                    s = transks->GetId();
                    (node+s)->SetFlag(0);
                } 
            }
            else
            {
                s = k;
                (node+s)->SetFlag(0); 
            }      
        } while(matik.GetNext());
        for(transij = (trans+i)->GetNext(); transij != NULL; transij = transij->GetNext())
        {
            j = transij->GetId();
            matjk = matrix->GetStart(j);
            do
            {
                k = matjk.GetIndex();
                if((node+k)->IsFGNode())
                {
                    for(transks = (trans+k)->GetNext(); transks != NULL; transks = transks->GetNext())
                    {
                        s = transks->GetId();
                        (node+s)->SetFlag(0);
                    } 
                }
                else
                {
                    s = k;
                    (node+s)->SetFlag(0); 
                }      
            } while(matjk.GetNext());
        }
    }

    nl = offset;
    index = (int *) CMGGetMem(nl*sizeof(int), CMG_FROM_BOTTOM);
    if (index == NULL) return 1;
    
    entry = (double *) CMGGetMem(nl*sizeof(double), CMG_FROM_BOTTOM);
    if (entry == NULL) return 1;
    
    // save entries
    for(i = 0; i < n; i++)
    {        
        if((node+i)->IsFGNode()) continue;
        nsaved = 0;
        firstindex = index+start[i];
        firstentry = entry+start[i];
        matik = matrix->GetStart(i);
        do
        {
            k = matik.GetIndex();
            mik = matik.GetData();
            if((node+k)->IsFGNode())
            {
                for(transks = (trans+k)->GetNext(); transks != NULL; transks = transks->GetNext())
                {
                    s = transks->GetId();
                    tks = transks->GetData();
                    mis = mik*tks;
                    // does mis already exist ?
                    ind = firstindex;
                    mat = firstentry;
                    found = 0;
                    for(ii = 0; ii < nsaved; ii++)
                    {
                        if( (*ind) == s) {(*mat) += mis; found = 1; break;}
                        ind++; mat++;
                    }
                    if(!found) // no
                    {
                        (*ind) = s;
                        (*mat) = mis;
                        nsaved++;
                    }
                        
                 } 
            }
            else
            {
                mis = mik;
                s = k;
                // does mis already exist ?
                ind = firstindex;
                mat = firstentry;
                found = 0;
                for(ii = 0; ii < nsaved; ii++)
                {
                    if( (*ind) == s) {(*mat) += mis; found = 1; break;}
                    ind++; mat++;
                }
                if(!found) // no
                {
                    (*ind) = s;
                    (*mat) = mis;
                    nsaved++;
                }
             }      
        } while(matik.GetNext());
        for(transij = (trans+i)->GetNext(); transij != NULL; transij = transij->GetNext())
        {
            j = transij->GetId();
            tij = transij->GetData();
            matjk = matrix->GetStart(j);
            do
            {
                k = matjk.GetIndex();
                mjk = matjk.GetData();
                if((node+k)->IsFGNode())
                {
                    for(transks = (trans+k)->GetNext(); transks != NULL; transks = transks->GetNext())
                    {
                        s = transks->GetId();
                        tks = transks->GetData();
                        mis = tij*mjk*tks;
                        // does mis already exist ?
                        ind = firstindex;
                        mat = firstentry;
                        found = 0;
                        for(ii = 0; ii < nsaved; ii++)
                        {
                            if( (*ind) == s) {(*mat) += mis; found = 1; break;}
                            ind++; mat++;
                        }
                        if(!found) // no
                        {
                            (*ind) = s;
                            (*mat) = mis;
                            nsaved++;
                        }
                    } 
                }
                else
                {
                    s = k;
                    mis = tij*mjk;
                    // does mis already exist ?
                    ind = firstindex;
                    mat = firstentry;
                    found = 0;
                    for(ii = 0; ii < nsaved; ii++)
                    {
                        if( (*ind) == s) {(*mat) += mis; found = 1; break;}
                        ind++; mat++;
                    }
                    if(!found) // no
                    {
                        (*ind) = s;
                        (*mat) = mis;
                        nsaved++;
                    }
                }      
            } while(matjk.GetNext());
        }
        if((start[i]+nsaved) != start[i+1])
        {  
            ostrstream ostr; ostr << __FILE__ << __LINE__ << endl;
            CMGError(ostr);
            return(1);
        }
    }

    ModifyIndex(0);
    if(ConstructEnd2()) return 1;
    if(OrderColumns2(graph)) return 1;
    if(ConstructAdjoined2(graph)) return 1;

    return 0;
}

void CMGMatrix::MarkUnknowns(CMGGraph *graph)
{
    CMGIndexBitField *bf;
    int id;
    CMGNode *node = graph->GetNode();
    int *ind = index;
    
    for(int i = 0; i < nl; i++)
    {
        bf = (CMGIndexBitField *) ind;
        id = bf->id;
        if((node+id)->IsFGNode()) bf->type = 1;
        ind++;
    }

    return;
}


int CMGMatrix::Init(int nn)
{
    n = nn;
    start = (int *) CMGGetMem((nn+1)*sizeof(int), CMG_FROM_TOP);
    if (start == NULL) return 1;

    nl = 0;
    entry = NULL;
    index = NULL; 
    adjoined = NULL;

    return 0;
}

int CMGMatrix::Init2(int nn) // memory from BOTTOM
{
    n = nn;
    start = (int *) CMGGetMem((nn+1)*sizeof(int), CMG_FROM_BOTTOM);
    if (start == NULL) return 1;

    nl = 0;
    entry = NULL;
    index = NULL; 
    adjoined = NULL;

    return 0;
}
        
void CMGMatrix::Mult(double *vout, double *vin)
{
    CMGMatrixPtr matij;
    double sum;
    int i;

    for(i = 0; i < n; i++)
    {
        matij = GetStart(i);
        sum = matij.GetData()*vin[i];
        while(matij.GetNext())
        {
            sum += matij.GetData()*vin[matij.GetIndex()];
        }
        vout[i] = sum;
    }

    return;
}

void CMGMatrix::MultTrans(double *vout, double *vin)
{
    CMGMatrixPtr matij;
    double sum;
    int i;

    for(i = 0; i < n; i++)
    {
        matij = GetStart(i);
        sum = matij.GetData()*vin[i];
        while(matij.GetNext())
        {
            sum += matij.GetAdjData()*vin[matij.GetIndex()];
        }
        vout[i] = sum;
    }

    return;
}


int CMGMatrix::Order(int *mapping)
{
    int i, *helparray, *ind;

    CMGMarkHeap(CMG_FROM_TOP);
    helparray = (int *) CMGGetMem(n*sizeof(int), CMG_FROM_TOP);
    if (helparray == NULL) return 1;
    for(i = 0; i < n; i++) helparray[i] = start[i];
    for(i = 0; i < n; i++) start[mapping[i]] = helparray[i];
    for(i = 0; i < n; i++) helparray[i] = end[i];
    for(i = 0; i < n; i++) end[mapping[i]] = helparray[i];
    CMGReleaseHeap(CMG_FROM_TOP);

    ind = index;
    CMGIndexBitField *bf;
    for(i = 0; i < nl; i++)
    {
        bf = (CMGIndexBitField *)ind;
        bf->id = mapping[bf->id];
        ind++;
    }
    
    return 0;
}

int CMGMatrix::Reorder(int *mapping)
{
    int i, *helparray, *ind;

    CMGMarkHeap(CMG_FROM_TOP);
    helparray = (int *) CMGGetMem(n*sizeof(int), CMG_FROM_TOP);
    if (helparray == NULL) return 1;
    for(i = 0; i < n; i++) helparray[i] = start[i];
    for(i = 0; i < n; i++) start[i] = helparray[mapping[i]];
    for(i = 0; i < n; i++) helparray[i] = end[i];
    for(i = 0; i < n; i++) end[i] = helparray[mapping[i]];
 
    // reverse mapping
    for(i = 0; i < n; i++) helparray[mapping[i]] = i;
    ind = index;
    CMGIndexBitField *bf;
    for(i = 0; i < nl; i++)
    {
       bf = (CMGIndexBitField *)ind;
       bf->id = helparray[bf->id];
       ind++;
    }
    CMGReleaseHeap(CMG_FROM_TOP);
    
    return 0;
}
    




