/****************************************************************************/
/*																			*/
/* File:      famg_graph.C													*/
/*																			*/
/* Purpose:   famg graph classes functions									*/
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
#include "famg_algebra.h"
#include "famg_heap.h"
#include "famg_grid.h"
#include "famg_graph.h"
#include "famg_fifo.h"

/* RCS_ID
$Header$
*/
   
    // Class List

void FAMGList::Insert(FAMGNode *node)
{   
    node->SetSucc(first);
    node->SetPred(NULL);
    if (first != NULL) first->SetPred(node);
    if (last == NULL) last = node;
    first = node;
    node->SetList(this);
}

void FAMGList::Init(FAMGList *p,FAMGList *s,int d)
{
     data = d;
     first = NULL;
     last = NULL;
     pred = p;
     succ = s;
     if(p != NULL) p->SetSucc(this);
     if(s != NULL) s->SetPred(this);
}

    // Class Node

void FAMGNode::Init(int index, const FAMGVectorEntry &i)
{
    SetData(0);
    SetId(index);
	SetLocalId(-1);
    SetVec(i);
    SetList(NULL);
    SetPaList(NULL);
    SetSucc(NULL);
    SetPred(NULL);
    SetFlag(0);
    SetFlag1(0);
    SetFlag2(0);
    SetFlag(0);
#ifdef ModelP
    SetFlagNewMarked(0);
#endif
	SetUndecidedNode();
    control.ns = 0;
    //SetLocalId is initialized with the default constructor
}
    

// Class PaList


void FAMGPaList::Init(FAMGPaList *next, int np, const int p[], double c[], double ct[], double error)
{
    int z;
    
    SetNext(next);
    SetNp(np);
    SetApprox(error);
	
    for(z = 0; z < np; z++)
    {
        SetPa(z,p[z]);
        SetCoeff(z,c[z]);
        SetCoefft(z,ct[z]);
    }

    return;
}


    // Class Graph

 
int FAMGGraph::Insert(FAMGNode *nod)
{
    FAMGList *last, *li, *pl;
    int data;

    data = nod->GetData();
    last = NULL;
    for(li = list; li != NULL; li = li->GetSucc())
    {
        if (data == li->GetData()) 
        {
            li->Insert(nod);
            return 0;
        }
        if (data < li->GetData())
        {
            if (freelist != NULL)
            {
                pl = freelist;
                freelist = pl->GetSucc();
            }
            else
            {
                pl = (FAMGList *) FAMGGetMem(sizeof(FAMGList), FAMG_FROM_BOTTOM);
            }
            if (pl == NULL) return 1;
            if(li->GetPred() == NULL) list = pl;
            pl->Init(li->GetPred(),li,data); 
            pl->Insert(nod);
            return 0;
        }
        last = li;
    }
    if (freelist != NULL)
    {
        pl = freelist;
        freelist = pl->GetSucc();
    }
    else
    {
        pl = (FAMGList *)FAMGGetMem(sizeof(class FAMGList), FAMG_FROM_BOTTOM);
    }
    if (pl == NULL) return 1;
    if(last == NULL) list = pl; 
    pl->Init(last,NULL,data);
    pl->Insert(nod);

    return 0;
}

void FAMGGraph::Remove(FAMGNode *nod)
{
    FAMGList *li;
    
    if(nod->GetList() == NULL)  return; // nod not on list
        
    if (nod->GetSucc() != NULL) 
    {
        (nod->GetSucc())->SetPred(nod->GetPred());
    }
    else
    {
        (nod->GetList())->SetLast(nod->GetPred());
    }
    if (nod->GetPred() != NULL)
    {
       (nod->GetPred())->SetSucc(nod->GetSucc()); 
    }
    else 
    {
        li = nod->GetList();
        if(nod->GetSucc() != NULL)
        {
            li->SetFirst(nod->GetSucc());
        }
        else 
        {
            if( li->GetPred() != NULL)
            {
                (li->GetPred())->SetSucc(li->GetSucc());
            }
            else
            {
                list = li->GetSucc();
            }
            if( li->GetSucc() != NULL)
            {
                (li->GetSucc())->SetPred(li->GetPred());
            }

            // add li to freelist
            li->SetSucc(freelist);
            freelist = li;
        }
    }
    nod->SetSucc(NULL);
    nod->SetPred(NULL);
    nod->SetList(NULL);
}

void FAMGGraph::Store(FAMGNode *nod)
{
    if (nod->GetList() != NULL)
    {
        Remove(nod);
        nod->SetPred(NULL);
        nod->SetSucc(helplist);
        helplist = nod;
    }

    return;
}

int FAMGGraph::InsertHelplist()
{
    FAMGNode *nextnode, *oldnode;

    nextnode = helplist;
    while(nextnode != NULL)
    {
        oldnode = nextnode;
        nextnode = oldnode->GetSucc();
        if( (!(oldnode->IsCGNode())) && (!(oldnode->IsFGNode())))
        {
            oldnode->ComputeTotalWeight();
            if(Insert(oldnode)) return 1;
        }
    }

    helplist = NULL;

    return 0;
}
       
FAMGNode *FAMGGraph::GetFirstNode()
{
    FAMGNode *nod;
    
    if(list == NULL) return NULL;
    nod = list->GetFirst();
    Remove(nod);

    return nod;
}

FAMGNode *FAMGGraph::GetLastNode()
{
    FAMGNode *nod;
    
    if(list == NULL) return NULL;
    nod = list->GetLast();
    Remove(nod);

    return nod;
}

void FAMGGraph::ClearPaList(FAMGPaList *pl)
{
    FAMGPaList *opl;

    while(pl != NULL)
    {
        opl = pl;
        pl = pl->GetNext();

        // store pl in freelist
        opl->SetNext(freepalist);
        freepalist = opl;
    }

    return;
}
 
void FAMGGraph::ClearPaListRev(FAMGPaList *&palist)
{
    FAMGPaList *opl, *pl;

    pl = palist;
    while(pl != NULL)
    {
        opl = pl;
        pl = pl->GetNext();

        // store pl in freelist
        opl->SetNext(freepalist);
        freepalist = opl;
    }
    palist = NULL;

    return;
} 

void FAMGGraph::CorrectPaList(FAMGPaList *&palist, double threshold)
{
    FAMGPaList *pl, *opl, *prev;

    pl = palist;
    prev = NULL;
    while(pl != NULL)
    {
        if(pl->GetApprox() > threshold)
        {
            opl = pl;
            pl = pl->GetNext();

            // store pl in freelist
            opl->SetNext(freepalist);
            freepalist = opl;
            
            if(prev != NULL) prev->SetNext(pl);
            else palist = pl;
        }
        else
        {
            prev = pl;
            pl = pl->GetNext();
        }
    }

    return;
} 

int FAMGGraph::SavePaList(FAMGPaList *&list, int np, const int pa[], double c[], double ct[], double error)
{
    FAMGPaList *pl;

    if (GetFreePaList() != NULL)
    {
        pl = GetFreePaList();
        SetFreePaList(pl->GetNext());
    }
    else
    {
        pl = (FAMGPaList *) FAMGGetMem(sizeof(class FAMGPaList), FAMG_FROM_BOTTOM);
        if (pl == NULL) return 1;
    }

    pl->Init(list,np,pa,c,ct,error);

    list = pl;

    return 0;
}   

void FAMGGraph::MarkFGNode(FAMGNode *fgnode)
{
    if(fgnode->IsFGNode()) return;
    fgnode->NodeMarkFG();
    // map[fgnode->GetId()] = nf;
    nf++;
	GetGridVector().SetFG(fgnode->GetVec());
	
    return;
}

void FAMGGraph::MarkCGNode(FAMGNode *cgnode)
{
    if(cgnode->IsCGNode())
		return;
    cgnode->NodeMarkCG();
	GetGridVector().SetCG(cgnode->GetVec());

    return;
}

int FAMGGraph::Init(FAMGGrid *gridptr)
{
    FAMGNode *nodei;
    int i, nrvec;
	FAMGVectorIter viter(gridptr->GetGridVector());
	FAMGVectorEntry ve;
	
    nrvec = gridptr->GetN();
    nf = 0;

    //node = (FAMGNode *) FAMGGetMem(n*sizeof(FAMGNode), FAMG_FROM_BOTTOM);
    node = new FAMGNode[nrvec];
    if (node == NULL)
		return 1;
    
	nodei = node;
	i=0;
	while(viter(ve))
		nodei++ ->Init(i++,ve);
	    
	assert(i==nrvec);
	n = nrvec;
	
#ifdef FAMG_ILU	
    map = gridptr->GetMap();
    if(map == NULL) // temporary map on level > 0
    {
        map = (int *) FAMGGetMem(n*sizeof(int),FAMG_FROM_BOTTOM);
        if (map == NULL)
			return 1;
        for(i = 0; i < nrvec; i++) map[i] = i;
    }
#endif
	
    SetList(NULL);
    SetHelpList(NULL);
    SetFreeList(NULL);
    SetFreePaList(NULL);
	SetGridVector(&(gridptr->GetGridVector()));

    // just for plotting

    //     FAMGMatrixEntry *row, *rowi, *matij;

    //     row = gridptr->GetMatrix()->GetRow();
    //     for(i = 0; i < n; i++)
    //     {
    //        rowi = row+i;
        //         rowi->SetMark(0);
    //        for(matij = rowi->GetNext(); matij != NULL; matij = matij->GetNext())
    //        {
    //            matij->SetParents(0);
    //        }
    //    }

    return 0;
}

#ifdef FAMG_ILU
int FAMGGraph::OrderILUT(FAMGMatrix *matrix)
{
    FAMGNode *nodei;
    FAMGMatrixPtr matij;
    int i, j, z, nl;
    //    double mii;

    if(n < 1) return 0;

    for(i = 0; i < n; i++)
    {
        matij = matrix->GetStart(i);
        nl = 0;
         while(matij.GetNext())
          {
             nl++;
         }
        nodei = node+i;
        nodei->Init(i);
        nodei->SetData(nl);
        if(Insert(nodei)) return 1;
    }


    FAMGMarkHeap(FAMG_FROM_BOTTOM);

    void **buffer = (void **) FAMGGetMem(sizeof(FAMGNode *)*n,FAMG_FROM_BOTTOM);
    if (buffer == NULL) {FAMGReleaseHeap(FAMG_FROM_BOTTOM);  return 1;}

    FAMGFifo fifo(buffer,sizeof(FAMGNode *)*n);

    
    fifo.In((void *) (node+0));
    (node+0)->SetFlag(1);
    while(!fifo.Empty())
    {
        nodei = (FAMGNode *) fifo.Out();
        i = nodei->GetId();
        matij = matrix->GetStart(i);
        // mii = matij.GetData();
        while(matij.GetNext())
        {
            j = Abs(matij.GetIndex());
            //   if(Abs(matij.GetData()) < 1e-15*mii) continue;
            if(!(node+j)->GetFlag())
            {
                (node+j)->SetFlag(1);
                fifo.In((void *)(node+j));
            }
        }
    }
    

    // nodei = GetFirstNode(); => nodei has smallest number of links

    for(i = 0; i < n; i++) (node+i)->SetFlag(1); // for reduzible matrices
    
    z = 0;
    while(nodei != NULL)
    {
        fifo.In((void *)nodei);
        (nodei)->SetFlag(0);
        while(!fifo.Empty())
        {        
            nodei = (FAMGNode *)fifo.Out();
            i = nodei->GetId();
            map[i] = z;
            Remove(nodei);
            z++;
            matij = matrix->GetStart(i);
            // mii = Abs(matij.GetData());
            while(matij.GetNext())
            {
                //  if(Abs(matij.GetData()) < 1e-15*mii) continue;
                j = matij.GetIndex();
                if((node+j)->GetFlag())
                {
                    (node+j)->SetFlag(0);
                    fifo.In((void *)(node+j));
                }
            }
        }
        nodei = GetFirstNode();
    }

    FAMGReleaseHeap(FAMG_FROM_BOTTOM);
 
    if(z != n) return 1;

    return 0;
}
#endif

