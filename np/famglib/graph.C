/****************************************************************************/
/*																			*/
/* File:      graph.C														*/
/*																			*/
/* Purpose:   cmg graph classes functions									*/
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
#include "graph.h"
#include "fifo.h"

/* RCS_ID
$Header$
*/
   
    // Class List

void CMGList::Insert(CMGNode *node)
{   
    node->SetSucc(first);
    node->SetPred(NULL);
    if (first != NULL) first->SetPred(node);
    if (last == NULL) last = node;
    first = node;
    node->SetList(this);
}

void CMGList::Init(CMGList *p,CMGList *s,int d)
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

void CMGNode::Init(int i)
{
    data = 0;
    id = i;
    list = NULL;
    palist = NULL;
    succ = NULL;
    pred = NULL;
    control.f0 = 0;
    control.f1 = 0;
    control.f2 = 0;
    control.nt = 0;
    control.ns = 0;
    control.lid = -1;
}
    

// Class PaList


void CMGPaList::Init(CMGPaList *nex, int n, int *p, double *c, double *ct, double error)
{
    int z;
    
    next = nex;
    np = n;
    for(z = 0; z < n; z++)
    {
        pa[z] = p[z];
        coeff[z] = c[z];
        coefft[z] = ct[z];
    }
    approx = error;

    return;
}


    // Class Graph

 
int CMGGraph::Insert(CMGNode *nod)
{
    CMGList *last, *li, *pl;
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
                pl = (CMGList *) CMGGetMem(sizeof(CMGList), CMG_FROM_BOTTOM);
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
        pl = (CMGList *)CMGGetMem(sizeof(class CMGList), CMG_FROM_BOTTOM);
    }
    if (pl == NULL) return 1;
    if(last == NULL) list = pl; 
    pl->Init(last,NULL,data);
    pl->Insert(nod);

    return 0;
}

void CMGGraph::Remove(CMGNode *nod)
{
    CMGList *li;
    
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

void CMGGraph::Store(CMGNode *nod)
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

int CMGGraph::InsertHelplist()
{
    CMGNode *nextnode, *oldnode;

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
       
CMGNode *CMGGraph::GetFirstNode()
{
    CMGNode *nod;
    
    if(list == NULL) return NULL;
    nod = list->GetFirst();
    Remove(nod);

    return nod;
}

CMGNode *CMGGraph::GetLastNode()
{
    CMGNode *nod;
    
    if(list == NULL) return NULL;
    nod = list->GetLast();
    Remove(nod);

    return nod;
}

void CMGGraph::ClearPaList(CMGPaList *pl)
{
    CMGPaList *opl;

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
 
void CMGGraph::ClearPaListRev(CMGPaList *&palist)
{
    CMGPaList *opl, *pl;

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

void CMGGraph::CorrectPaList(CMGPaList *&palist, double threshold)
{
    CMGPaList *pl, *opl, *prev;

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

int CMGGraph::SavePaList(CMGPaList *&list, int np, int *pa, double *c, double *ct, double error)
{
    CMGPaList *pl;

    if (freepalist != NULL)
    {
        pl = freepalist;
        freepalist = pl->GetNext();
    }
    else
    {
        pl = (CMGPaList *) CMGGetMem(sizeof(class CMGPaList), CMG_FROM_BOTTOM);
        if (pl == NULL) return 1;
    }

    pl->Init(list,np,pa,c,ct,error);

    list = pl;

    return 0;
}   



void CMGGraph::MarkFGNode(CMGNode *fgnode)
{
    if(fgnode->IsFGNode()) return;
    fgnode->MarkFGNode();
    // map[fgnode->GetId()] = nf;
    nf++;

    return;
}

void CMGGraph::MarkCGNode(CMGNode *cgnode)
{
    if(cgnode->IsCGNode()) return;
    cgnode->MarkCGNode();
    // map[cgnode->GetId()] = nc; 
    nc--; 

    return;
}

int CMGGraph::Init(CMGGrid *gridptr)
{
    CMGNode *nodei;
    int i;

    // noetig ?
    n = gridptr->GetN();
    nf = 0;
    nc = n-1;

    node = (CMGNode *) CMGGetMem(n*sizeof(CMGNode), CMG_FROM_BOTTOM);
    if (node == NULL) return 1;
    
    for(i = 0; i < n; i++)
    {
        nodei = node+i;
        nodei->Init(i);
    }
    
    map = gridptr->GetMap();
    if(map == NULL) // temporary map on level > 0
    {
        map = (int *) CMGGetMem(n*sizeof(int),CMG_FROM_BOTTOM);
        if (map == NULL) return 1;
        for(i = 0; i < n; i++) map[i] = i;
    }

    list = NULL;
    helplist = NULL;
    freelist = NULL;
    freepalist = NULL;


    // just for plotting

    //     CMGMatrixEntry *row, *rowi, *matij;

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

    // test
 
int CMGGraph::InsertH(CMGNode *nod)
{
    CMGList *last, *li, *pl;
    int data;

    data = nod->GetData();
    last = NULL;
    for(li = (CMGList *)helplist; li != NULL; li = li->GetSucc())
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
                pl = (CMGList *) CMGGetMem(sizeof(CMGList), CMG_FROM_BOTTOM);
            }
            if (pl == NULL) return 1;
            if(li->GetPred() == NULL) helplist = (CMGNode *)pl;
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
        pl = (CMGList *)CMGGetMem(sizeof(class CMGList), CMG_FROM_BOTTOM);
    }
    if (pl == NULL) return 1;
    if(last == NULL) helplist = (CMGNode *)pl; 
    pl->Init(last,NULL,data);
    pl->Insert(nod);

    return 0;
}

void CMGGraph::RemoveH(CMGNode *nod)
{
    CMGList *li;
    
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
                helplist = (CMGNode *)li->GetSucc();
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

CMGNode *CMGGraph::GetFirstNodeH()
{
    CMGNode *nod;
    CMGList *hlist;
    
    if(helplist == NULL) return NULL;
    hlist = (CMGList *)helplist;
    nod = hlist->GetFirst();
    RemoveH(nod);

    return nod;
}



int CMGGraph::OrderSpecial(CMGMatrix *matrix)
{
    CMGNode *nodei, *nodej;
    CMGMatrixPtr matij;
    int i, z, nl;
    double mii;

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
        nodei->SetFlag1(1);
        nodei->SetData(nl);
        if(Insert(nodei)) return 1;
    }

    z = 0;

    while(list != NULL)
    {
        nodei = GetFirstNode();
        while(nodei != NULL)
        {
            nodei->SetFlag1(0);
            i = nodei->GetId();
            map[i] = z;
            z++;
            matij = matrix->GetStart(i);
            mii = matij.GetData();
            while(matij.GetNext())
            {
                if(Abs(matij.GetData()) < 1e-5*mii) continue;
                nodej = node+matij.GetIndex();
                if((nodej->GetFlag1() == 1) && (nodej->GetFlag2() == 0))
                {
                    Remove(nodej);
                    nodej->SetFlag1(0);
                    nodej->SetFlag2(1);
                    if(InsertH(nodej)) return 1;
                }   
            }   
            nodei = GetFirstNode();
        }
        
        nodei = GetFirstNodeH();
        while(nodei != NULL)
        {
            nodei->SetFlag2(0);
            i = nodei->GetId();
            map[i] = z;
            z++;
            matij = matrix->GetStart(i);
            mii = matij.GetData();
            while(matij.GetNext())
            {
                if(Abs(matij.GetData()) < 1e-5*mii) continue;
                nodej = node+matij.GetIndex();
                if((nodej->GetFlag1() == 0) && (nodej->GetFlag2() == 1))
                {
                    RemoveH(nodej);
                    nodej->SetFlag1(1);
                    nodej->SetFlag2(0);
                    if(Insert(nodej)) return 1;
                }   
            }   
            nodei = GetFirstNodeH();
        }
    }
 
    return 0;
}
        
int CMGGraph::OrderILUT(CMGMatrix *matrix)
{
    CMGNode *nodei;
    CMGMatrixPtr matij;
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


    CMGMarkHeap(CMG_FROM_BOTTOM);

    void **buffer = (void **) CMGGetMem(sizeof(CMGNode *)*n,CMG_FROM_BOTTOM);
    if (buffer == NULL) {CMGReleaseHeap(CMG_FROM_BOTTOM);  return 1;}

    CMGFifo fifo(buffer,sizeof(CMGNode *)*n);

    
    fifo.In((void *) (node+0));
    (node+0)->SetFlag(1);
    while(!fifo.Empty())
    {
        nodei = (CMGNode *) fifo.Out();
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
            nodei = (CMGNode *)fifo.Out();
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

    CMGReleaseHeap(CMG_FROM_BOTTOM);
 
    if(z != n) return 1;

    return 0;
}


