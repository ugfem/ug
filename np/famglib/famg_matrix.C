/****************************************************************************/
/*																			*/
/* File:      famg_matrix.C													*/
/*																			*/
/* Purpose:   famg matrix classes functions									*/
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
#include "famg_misc.h"
#include "famg_matrix.h"
#include "famg_heap.h"
#include "famg_grid.h"

/* RCS_ID
$Header$
*/

static double famgzero = 0.0;


void FAMGMatrix::DevideFGDefect(double *unknown, double *defect)
{
    FAMGMatrixPtr matij;
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
 
void FAMGMatrix::VecMinusMatVec(double *defect, double *rhs, double *unknown)
{
    FAMGMatrixPtr matij;
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

void FAMGMatrix::JAC(double *vec)
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

void FAMGMatrix::FGS(double *vec)
{
    FAMGMatrixPtr matij;
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

void FAMGMatrix::BGS(double *vec)
{
    FAMGMatrixPtr matij;
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

void FAMGMatrix::SGS(double *vec)
{
    FAMGMatrixPtr matij;
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

int FAMGMatrix::GetSmallestIndex()
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


void FAMGMatrix::ModifyIndex(int f)
{
    FAMGIndexBitField bf;
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
        
void FAMGMatrix::RemodifyIndex(int f)
{
    FAMGIndexBitField *bf;
    int *ind;
    
    ind = index;
    for(int i = 0; i < nl; i++)
    {
        bf  = (FAMGIndexBitField *)ind;
        (*ind) = (bf->id+f);
        ind++;
    }

    return;
}
        
void FAMGMatrix::ModifyIndex(int *type, int f)
{
    FAMGIndexBitField bf;
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
        
void FAMGMatrix::RemodifyIndex(int *type, int f)
{
    FAMGIndexBitField *bf;
    int *ind, ii;
    
    ind = index; 
    for(int i = 0; i < nl; i++)
    {
        bf  = (FAMGIndexBitField *)ind;
        ii =  bf->id;
        type[ii] = bf->type;
        (*ind) = ii+f;
        ind++; 
    }

    return;
}
        
#ifdef FAMG_REORDERCOLUMN			
int FAMGMatrix::OrderColumns(int *map)
{
    double hd, *mat, *firstentry;
    int i, j, nc;
    FAMGIndexBitField *firstindex, *ind, hi;

    for(i = 0; i < n; i++)
    {
        nc = start[i+1]-start[i];
        firstindex = ind =  (FAMGIndexBitField *) index+start[i];
        firstentry = mat = entry+start[i];
        for(j = 0; j < nc; j++)
        {
            if(ind->id == i) break;
            ind++; mat++;
        }
        if(j == nc)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": wrong matrix structure." << endl;
            FAMGError(ostr);
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

int FAMGMatrix::OrderColumns()
{
    double hd, *mat, *firstentry;
    int i, j, nc;
    FAMGIndexBitField *firstindex, *ind, hi;

    for(i = 0; i < n; i++)
    {
        nc = start[i+1]-start[i];
        firstindex = ind =  (FAMGIndexBitField *) index+start[i];
        firstentry = mat = entry+start[i];
        for(j = 0; j < nc; j++)
        {
            if(ind->id == i) break;
            ind++; mat++;
        }
        if(j == nc)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": wrong matrix structure." << endl;
            FAMGError(ostr);
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

int FAMGMatrix::ReorderColumns(int *map)
{
    double hd, *mat, *firstentry;
    int i;
    FAMGIndexBitField *firstindex, *ind, hi;

    for(i = 0; i < n; i++)
    {
        if(map[i] == 0) continue;
        firstindex = (FAMGIndexBitField *) index+start[i];
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

int FAMGMatrix::OrderColumns2(FAMGGraph *graph)
{
    FAMGNode *node;
    double hd, *mat, *firstentry;
    int i, j, nc;
    FAMGIndexBitField *firstindex, *ind, hi;

    node = graph->GetNode();

    for(i = 0; i < n; i++)
    {
        if((node+i)->IsFGNode()) continue;
        nc = start[i+1]-start[i];
        firstindex = ind = (FAMGIndexBitField *) index+start[i];
        firstentry = mat = entry+start[i];
        for(j = 0; j < nc; j++)
        {
            if(ind->id == i) break;
            ind++; mat++;
        }
        if(j == nc)
        {
            ostrstream ostr; ostr << __FILE__ << ", line " << __LINE__ << ": wrong matrix structure." << endl;
            FAMGError(ostr);
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
#endif

int FAMGMatrix::ConstructAdjoined()
{
    FAMGMatrixPtr matij, matjk;
    int i, j, k, found;

    adjoined = (double **) FAMGGetMem(nl*sizeof(double *), FAMG_FROM_TOP);
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
                    *(matij.GetAdjoinedPtr()) = &famgzero;
                }
            }
        } while(matij.GetNext());

    }

    return 0;
}

int FAMGMatrix::ConstructAdjoinedB()
{
    FAMGMatrixPtr matij, matjk;
    int i, j, k, found;

    adjoined = (double **) FAMGGetMem(nl*sizeof(double *), FAMG_FROM_BOTTOM);
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
                    *(matij.GetAdjoinedPtr()) = &famgzero;
                }
            }
        } while(matij.GetNext());

    }

    return 0;
}


int FAMGMatrix::ConstructAdjoined2(FAMGGraph *graph)
{
    FAMGMatrixPtr matij, matjk;
    int i, j, k, found;
    FAMGNode *node = graph->GetNode();

    adjoined = (double **) FAMGGetMem(nl*sizeof(double *), FAMG_FROM_BOTTOM);
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
                    *(matij.GetAdjoinedPtr()) = &famgzero;
                    // return(1);
                }
            }
        } while(matij.GetNext());

    }

    return 0;
}


int FAMGMatrix::ConstructEnd()
{
    int i;

    end = (int *) FAMGGetMem(sizeof(int)*n,FAMG_FROM_TOP);
    if(end == NULL) return 1;
     
    for(i = 0; i < n; i++) end[i] = start[i+1];

    return 0;
}
    
int FAMGMatrix::ConstructEndB()
{
    int i;

    end = (int *) FAMGGetMem(sizeof(int)*n,FAMG_FROM_BOTTOM);
    if(end == NULL) return 1;
     
    for(i = 0; i < n; i++) end[i] = start[i+1];

    return 0;
}
    
int FAMGMatrix::ConstructEnd2()
{
    int i;

    end = (int *) FAMGGetMem(sizeof(int)*n,FAMG_FROM_BOTTOM);
    if(end == NULL) return 1;
     
    for(i = 0; i < n; i++) end[i] = start[i+1];

    return 0;
}
    
int FAMGMatrix::CGMatrix(FAMGMatrix *fgmatrix, FAMGTransfer *transfer, int *father)
{    
    FAMGMatrixPtr matik, matjk;
    FAMGTransferEntry *trans, *transij, *transks;
    double mis, tij, mjk, tks, mik, *firstentry, *mat;
    int i, j, s, k, z, nsaved, ii, offset, fnn, *fathermap, *tmpflag, *firstindex, *ind, found;
    
    fnn = fgmatrix->GetN();
    trans = transfer->GetRow();

    FAMGMarkHeap(FAMG_FROM_BOTTOM);

    fathermap = (int *) FAMGGetMem(sizeof(int)*fnn,FAMG_FROM_BOTTOM);
    if (fathermap == NULL) {FAMGReleaseHeap(FAMG_FROM_BOTTOM);  return 1;}

    tmpflag = (int *) FAMGGetMem(sizeof(int)*n,FAMG_FROM_BOTTOM);
    if (tmpflag == NULL) {FAMGReleaseHeap(FAMG_FROM_BOTTOM);  return 1;}

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
    index = (int *) FAMGGetMem(nl*sizeof(int), FAMG_FROM_TOP);
    if (index == NULL) return 1;
    
    entry = (double *) FAMGGetMem(nl*sizeof(double), FAMG_FROM_TOP);
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
            FAMGError(ostr);
            return(1);
        }
    }

    ModifyIndex(0);
    if(ConstructEnd()) return 1;
    if(OrderColumns()) return 1;
    if(ConstructAdjoined()) return 1;

    FAMGReleaseHeap(FAMG_FROM_BOTTOM);

    return 0;
}


int FAMGMatrix::TmpMatrix(FAMGMatrix *matrix, FAMGTransfer *transfer, FAMGGraph *graph)
{    
    FAMGMatrixPtr matik, matjk;
    FAMGTransferEntry *trans, *transij, *transks;
    FAMGNode *node;
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
    index = (int *) FAMGGetMem(nl*sizeof(int), FAMG_FROM_BOTTOM);
    if (index == NULL) return 1;
    
    entry = (double *) FAMGGetMem(nl*sizeof(double), FAMG_FROM_BOTTOM);
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
            FAMGError(ostr);
            return(1);
        }
    }

    ModifyIndex(0);
    if(ConstructEnd2()) return 1;
    if(OrderColumns2(graph)) return 1;
    if(ConstructAdjoined2(graph)) return 1;

    return 0;
}

void FAMGMatrix::MarkUnknowns(FAMGGraph *graph)
{
    FAMGIndexBitField *bf;
    int id;
    FAMGNode *node = graph->GetNode();
    int *ind = index;
    
    for(int i = 0; i < nl; i++)
    {
        bf = (FAMGIndexBitField *) ind;
        id = bf->id;
        if((node+id)->IsFGNode()) bf->type = 1;
        ind++;
    }

    return;
}


int FAMGMatrix::Init(int nn)
{
    n = nn;
    start = (int *) FAMGGetMem((nn+1)*sizeof(int), FAMG_FROM_TOP);
    if (start == NULL) return 1;

    nl = 0;
    entry = NULL;
    index = NULL; 
    adjoined = NULL;

    return 0;
}

int FAMGMatrix::Init2(int nn) // memory from BOTTOM
{
    n = nn;
    start = (int *) FAMGGetMem((nn+1)*sizeof(int), FAMG_FROM_BOTTOM);
    if (start == NULL) return 1;

    nl = 0;
    entry = NULL;
    index = NULL; 
    adjoined = NULL;

    return 0;
}
        
void FAMGMatrix::Mult(double *vout, double *vin)
{
    FAMGMatrixPtr matij;
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

void FAMGMatrix::MultTrans(double *vout, double *vin)
{
    FAMGMatrixPtr matij;
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


int FAMGMatrix::Order(int *mapping)
{
    int i, *helparray, *ind;

    FAMGMarkHeap(FAMG_FROM_TOP);
    helparray = (int *) FAMGGetMem(n*sizeof(int), FAMG_FROM_TOP);
    if (helparray == NULL) return 1;
    for(i = 0; i < n; i++) helparray[i] = start[i];
    for(i = 0; i < n; i++) start[mapping[i]] = helparray[i];
    for(i = 0; i < n; i++) helparray[i] = end[i];
    for(i = 0; i < n; i++) end[mapping[i]] = helparray[i];
    FAMGReleaseHeap(FAMG_FROM_TOP);

    ind = index;
    FAMGIndexBitField *bf;
    for(i = 0; i < nl; i++)
    {
        bf = (FAMGIndexBitField *)ind;
        bf->id = mapping[bf->id];
        ind++;
    }
    
    return 0;
}

int FAMGMatrix::Reorder(int *mapping)
{
    int i, *helparray, *ind;

    FAMGMarkHeap(FAMG_FROM_TOP);
    helparray = (int *) FAMGGetMem(n*sizeof(int), FAMG_FROM_TOP);
    if (helparray == NULL) return 1;
    for(i = 0; i < n; i++) helparray[i] = start[i];
    for(i = 0; i < n; i++) start[i] = helparray[mapping[i]];
    for(i = 0; i < n; i++) helparray[i] = end[i];
    for(i = 0; i < n; i++) end[i] = helparray[mapping[i]];
 
    // reverse mapping
    for(i = 0; i < n; i++) helparray[mapping[i]] = i;
    ind = index;
    FAMGIndexBitField *bf;
    for(i = 0; i < nl; i++)
    {
       bf = (FAMGIndexBitField *)ind;
       bf->id = helparray[bf->id];
       ind++;
    }
    FAMGReleaseHeap(FAMG_FROM_TOP);
    
    return 0;
}
    




