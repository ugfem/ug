/****************************************************************************/
/*																			*/
/* File:      construct.C													*/
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
#include "system.h"

   
#ifdef UG_DRAW

extern "C"
{
#include "wpm.h"
#include "wop.h"
#include "connectuggrape.h"
#include "uginterface.h"
#include "gm.h"
} 

#endif

/* RCS_ID
$Header$
*/

int CMGGrid::AnalyseNodeSimple(int i, CMGPaList *&palist)
{
    CMGNode *nodei, *node;
    CMGPaList *pl;
    int remove, z, j;

    palist = NULL;
    node = graph->GetNode();
    nodei = node+i;
    for(pl = nodei->GetPaList(); pl != NULL; pl = pl->GetNext())
    {
        remove = 0;
        for(z = 0; z < pl->GetNp(); z++)
        {
            j = pl->GetPa(z);
            if((node+j)->IsFGNode()) 
            {
                remove = 1;
                break;
            }
        }
        if(!remove)
        {
            
            if(graph->SavePaList(palist,pl->GetNp(),pl->GetPa(),pl->GetCoeff(),pl->GetCoefft(),pl->GetApprox())) return 1;
        }
    }

    return 0;
}


double CMGPaList::TotalWeight()
{
    double weight;

    weight = 0.0*approx;
    weight += 1.0*(double) newlinks;
    weight += 10.0*(double) newcg;

    return weight;
} 


void CMGNode::ComputeTotalWeight()
{
    CMGPaList *pl;
    double weight, minweight;

    minweight = 1e+7;
    for(pl  = palist; pl != NULL; pl = pl->GetNext())
    {
        weight = pl->TotalWeight();
        if (weight < minweight)
        {
            minweight = weight;
        }
    }

    // test
    double ns = (double) GetNSons();
    if (ns < 0.5) ns = 0.5;
    //  minweight -= 1.0/ns;
   
    data = (int) ceil(100.0*minweight);
  


    return;
}


void CMGNode::CountNewCG(CMGGraph *graph)
{
    CMGNode *node;
    CMGPaList *pl;
    double nc;
    int i, np, *pa, ns;

    node = graph->GetNode();
    for(pl  = palist; pl != NULL; pl = pl->GetNext())
    {
        np = pl->GetNp();
        pa = pl->GetPa();
        nc = 0;
        for(i = 0; i < np; i++)
        {
            if((node+pa[i])->IsCGNode()) continue;
            else
            {
                ns = (node+pa[i])->GetNSons();
                if(ns > 0)
                {
                    nc += 1.0/(double)ns;
                }
                else
                {
                     nc += 1.0;
                }
            }
        }
        pl->SetNewCG(nc);
    }
    
    return;
}

int CMGNode::CountNewCG(CMGGraph *graph, int j)
{
    CMGNode *node;
    CMGPaList *pl;
    double nc;
    int i, np, *pa, found, change, ns;

    node = graph->GetNode();
    change = 0;
    for(pl  = palist; pl != NULL; pl = pl->GetNext())
    {
        np = pl->GetNp();
        pa = pl->GetPa();
        nc = 0;
        found = 0;
        for(i = 0; i < np; i++)
        {
            if(pa[i] == j) found = 1;
        }
        if(found)
        {
            for(i = 0; i < np; i++)
            {
                if((node+pa[i])->IsCGNode()) continue;
                else
                {
                    ns = (node+pa[i])->GetNSons();
                    if(ns > 0)
                    {
                        nc += 1.0/(double)ns;
                    }
                    else
                    {
                        nc += 1.0;
                    }
                }
            }
            if(pl->GetNewCG() != nc)
            {
                pl->SetNewCG(nc);
                change = 1;
            }
        }
    }
    
    return change;
}


int CMGGrid::Connected(int i, int z)
{
    CMGTransferEntry *trans, *transis, *transrj;
    CMGMatrixPtr matsr;
    CMGNode *noder, *node;
    int j,r,s;

    trans = transfer->GetRow();
    node = graph->GetNode();

    for(transis = (trans+i)->GetNext(); transis != NULL; transis = transis->GetNext())
    {
        s = transis->GetId();
        matsr = matrix->GetStart(s); 
        do
        {
            r = matsr.GetIndex();
            if(r == z) return 1;
            noder = node+r;
            if(noder->IsFGNode())
            {
                for(transrj = (trans+r)->GetNext(); transrj != NULL; transrj = transrj->GetNext())
                {
                    j = transrj->GetId();
                    if(j == z) return 1;
                }
            }
        } while(matsr.GetNext());
    }
    s = i;
    matsr = matrix->GetStart(s); 
    do
    {
        r = matsr.GetIndex();
        if(r == z) return 1;
        noder = node+r;
        if(noder->IsFGNode())
        {
            for(transrj = (trans+r)->GetNext(); transrj != NULL; transrj = transrj->GetNext())
            {
                j = transrj->GetId();
                if(j == z) return 1;
            }
        }
    } while(matsr.GetNext());
                                          
    return 0;
}

int CMGGrid::SetFlagsAndCount(int i, int f)
{
    CMGTransferEntry *trans, *transis, *transrj;
    CMGMatrixPtr matsr;
    CMGNode *nodej, *noder, *node;
    int z,j,r,s;

    trans = transfer->GetRow();
    node = graph->GetNode();

    z = 0;
    for(transis = (trans+i)->GetNext(); transis != NULL; transis = transis->GetNext())
    {
        s = transis->GetId();
        matsr = matrix->GetStart(s); 
        do
        {
            r = matsr.GetIndex();
            noder = node+r;
            if(noder->IsFGNode())
            {
                for(transrj = (trans+r)->GetNext(); transrj != NULL; transrj = transrj->GetNext())
                {
                    j = transrj->GetId();
                    nodej = node + j;
                    if(nodej->GetFlag() != f)
                    {
                        nodej->SetFlag(f);
                        z++;
                    }
                }
            }
            else
            {
                if(noder->GetFlag() != f)
                {
                    noder->SetFlag(f);
                    z++;
                }
            }
        } while(matsr.GetNext());
    }
    s = i;
    matsr = matrix->GetStart(s); 
    do
    {
        r = matsr.GetIndex();
        noder = node+r;
        if(noder->IsFGNode())
        {
            for(transrj = (trans+r)->GetNext(); transrj != NULL; transrj = transrj->GetNext())
            {
                j = transrj->GetId();
                nodej = node + j;
                if(nodej->GetFlag() != f)
                {
                    nodej->SetFlag(f);
                    z++;
                }
            }
        }
        else
        {
            if(noder->GetFlag() != f)
            {
                noder->SetFlag(f);
                z++;
            }
        }
    } while(matsr.GetNext());
                    
    return z;
}

void CMGGrid::SetFlags(int i, int f)
{
    CMGTransferEntry *trans, *transis, *transrj;
    CMGMatrixPtr matsr;
    CMGNode *nodej, *noder, *node;
    int j,r,s;

    trans = transfer->GetRow();
    node = graph->GetNode();

    for(transis = (trans+i)->GetNext(); transis != NULL; transis = transis->GetNext())
    {
        s = transis->GetId();
        matsr = matrix->GetStart(s); 
        do
        {
            r = matsr.GetIndex();
            noder = node+r;
            if(noder->IsFGNode())
            {
                for(transrj = (trans+r)->GetNext(); transrj != NULL; transrj = transrj->GetNext())
                {
                    j = transrj->GetId();
                    nodej = node + j;
                    nodej->SetFlag(f);
                }
            }
            else
            {
                noder->SetFlag(f);
            }
        } while(matsr.GetNext());
    }
    s = i;
    matsr = matrix->GetStart(s); 
    do
    {
        r = matsr.GetIndex();
        noder = node+r;
        if(noder->IsFGNode())
        {
            for(transrj = (trans+r)->GetNext(); transrj != NULL; transrj = transrj->GetNext())
            {
                j = transrj->GetId();
                nodej = node + j;
                nodej->SetFlag(f);
            }
        }
        else
        {
            noder->SetFlag(f);
        }
    } while(matsr.GetNext());
                    
                       
    return;
}

int CMGGrid::CountLinks(int i)
{
    CMGTransferEntry *trans, *transis, *transrj;
    CMGMatrixPtr matsr;
    CMGNode *nodej, *noder, *node;
    int z,j,r,s;

    trans = transfer->GetRow();
    node = graph->GetNode();

    z = 0;
    for(transis = (trans+i)->GetNext(); transis != NULL; transis = transis->GetNext())
    {
        s = transis->GetId();
        matsr = matrix->GetStart(s); 
        do
        {
            r = matsr.GetIndex();
            noder = node+r;
            if(noder->IsFGNode())
            {
                for(transrj = (trans+r)->GetNext(); transrj != NULL; transrj = transrj->GetNext())
                {
                    j = transrj->GetId();
                    nodej = node + j;
                    if(nodej->GetFlag() == 1)
                    {
                        z++;
                        nodej->SetFlag(0);
                    }
                }
            }
            else
            {
                if(noder->GetFlag() == 1)
                {
                    z++;
                    noder->SetFlag(0);
                }
            }
        } while(matsr.GetNext());
    }
    s = i;
    matsr = matrix->GetStart(s); 
    do
    {
        r = matsr.GetIndex();
        noder = node+r;
        if(noder->IsFGNode())
        {
            for(transrj = (trans+r)->GetNext(); transrj != NULL; transrj = transrj->GetNext())
            {
                j = transrj->GetId();
                nodej = node + j;
                if(nodej->GetFlag() == 1)
                {
                    z++;
                    nodej->SetFlag(0);
                }
            }
        }
        else
        {
            if(noder->GetFlag() == 1)
            {
                z++;
                noder->SetFlag(0);
            }
        }
    } while(matsr.GetNext());
                                           
    return z;
}

int CMGNode::CountNewLinks(CMGGrid *gridptr, CMGGraph *graph)
{
    CMGPaList *pl;
    CMGNode *node;
    int nnb, nl, z, np, *pa, y, newlinks;
    
    node = graph->GetNode();
    nnb = gridptr->SetFlagsAndCount(id,1);
    if(GetFlag() == 1) 
    {
        SetFlag(0); 
        nnb--;
    }
    for(pl = palist; pl != NULL; pl = pl->GetNext())
    {
        np = pl->GetNp();
        pa = pl->GetPa();
        nl = 0;
        for(z = 0; z < np; z++)
        {
            (node+pa[z])->SetFlag(0);
            nl += gridptr->CountLinks(pa[z]);
            gridptr->SetFlags(id,1);
            SetFlag(0);
            for(y = z+1; y < np; y++)
            {
                if(!gridptr->Connected(pa[z],pa[y])) nl++;
            }
        }
        newlinks = np*(nnb-1)-nl;
        if(np < 0) newlinks = 20;
        pl->SetNewLinks(newlinks);
    }
    gridptr->SetFlags(id,0);

    return 0;
}


void CMGGraph::InitNSons()
{
    CMGPaList *palist, *pl;
    int i, np, j, z, mark;

    for(i = 0; i < n; i++)
    {
        palist = (node+i)->GetPaList();
        for(pl = palist; pl != NULL; pl = pl->GetNext())
        {
            np = pl->GetNp();
            mark = 0;
            for(z = 0; z < np; z++)
            {
                j = pl->GetPa(z);
                if((node+j)->IsCGNode())
                {
                    mark = 1;
                    break;
                }
            }
            if(mark)
            {    
                for(z = 0; z < np; z++)
                {
                    j = pl->GetPa(z);
                    if(!(node+j)->IsCGNode())
                    {
                        (node+j)->SetFlag(1);
                    }
                }
            }
        }
        for(pl = palist; pl != NULL; pl = pl->GetNext())
        {
            np = pl->GetNp();
            for(z = 0; z < np; z++)
            {
                j = pl->GetPa(z);
                if((node+j)->GetFlag() == 1)
                {
                    (node+j)->SetFlag(0);
                    (node+j)->SetNSons((node+j)->GetNSons()+1);
                }
            }
        }
    }

    return;
}
            
int CMGGrid::UpdateNBNewCG(int i)
{
    CMGNode *nodej, *node;
    CMGMatrixPtr matij;
    int j;

    node = graph->GetNode();

    // check all neighbors with i as possible father, tmpmatrix is OK
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        nodej = node+j;
        if(((nodej)->IsCGNode()) || ((nodej)->IsFGNode())) continue;
        if((node+j)->CountNewCG(graph,i))
        {
            // something changed
            graph->Store(nodej);
        }
    }

    return 0;
}

void CMGGraph::UpdateNSons(CMGPaList *newlist, CMGPaList *oldlist, CMGGrid* grid)
{
    CMGPaList *pl;
    int np, z, j, mark, fmark;
    
    for(pl = oldlist; pl != NULL; pl = pl->GetNext())
    {
        np = pl->GetNp();
        mark = 0; fmark = 1;
        for(z = 0; z < np; z++)
        {
            j = pl->GetPa(z);
            if((node+j)->IsCGNode()) mark = 1;
            if((node+j)->IsFGNode()) 
            {
                fmark = 0;
                break;
            }
        }
        if(mark && fmark )
        {
            for(z = 0; z < np; z++)
            {
                j = pl->GetPa(z);
                if(!(node+j)->IsCGNode())
                {
                    (node+j)->SetFlag(1);
                }
            }
        }
    }

    for(pl = newlist; pl != NULL; pl = pl->GetNext())
    {
        np = pl->GetNp();
        mark = 0; fmark = 1;
        for(z = 0; z < np; z++)
        {
            j = pl->GetPa(z);
            if((node+j)->IsCGNode()) mark = 1;
            if((node+j)->IsFGNode()) 
            {
                fmark = 0;
                break;
            }
        }
        if(mark && fmark)
        {
            for(z = 0; z < np; z++)
            {
                j = pl->GetPa(z);
                if(!(node+j)->IsCGNode())
                {
                    (node+j)->SetFlag1(1);
                }
            }
        }
    }

    for(pl = oldlist; pl != NULL; pl = pl->GetNext())
    {
        np = pl->GetNp();
        for(z = 0; z < np; z++)
        {
            j = pl->GetPa(z);
            if(((node+j)->GetFlag())
               && (!(node+j)->GetFlag1())
               && (!(node+j)->GetFlag2()))
            {
                (node+j)->SetFlag2(1);
                (node+j)->SetNSons((node+j)->GetNSons()-1);
            }
        }
    }

    for(pl = newlist; pl != NULL; pl = pl->GetNext())
    {
        np = pl->GetNp();
        for(z = 0; z < np; z++)
        {
            j = pl->GetPa(z);
            if((!(node+j)->GetFlag()) 
               && ((node+j)->GetFlag1())
               && (!(node+j)->GetFlag2()))
            {
                
               (node+j)->SetFlag2(1);
               (node+j)->SetNSons((node+j)->GetNSons()+1);
            }
        }
    }

    for(pl = oldlist; pl != NULL; pl = pl->GetNext())
    {
        np = pl->GetNp();
        for(z = 0; z < np; z++)
        {
            j = pl->GetPa(z);
            (node+j)->SetFlag(0);
            (node+j)->SetFlag1(0);
            if((node+j)->GetFlag2())
            {
                (node+j)->SetFlag2(0);
                grid->UpdateNBNewCG(j);
            }
        }
    }

    for(pl = newlist; pl != NULL; pl = pl->GetNext())
    {
        np = pl->GetNp();
        for(z = 0; z < np; z++)
        {
            j = pl->GetPa(z);
            (node+j)->SetFlag(0);
            (node+j)->SetFlag1(0);
            if((node+j)->GetFlag2())
            {
                (node+j)->SetFlag2(0);
                grid->UpdateNBNewCG(j);
            }
        }
    }

    return;
}


int CMGGrid::UpdateNeighborsCG(int i)
{
    CMGNode *nodej, *node;
    CMGMatrixPtr matij;
    CMGPaList *pl;
    int j, z, found, k;

    node = graph->GetNode();

    // check all neighbors with i as possible father, tmpmatrix is OK
    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        nodej = node+j;
        if(((nodej)->IsCGNode()) || ((nodej)->IsFGNode())) continue;
        {
            for(pl = nodej->GetPaList(); pl != NULL; pl = pl->GetNext())
            {
                found = 0;
                for(z = 0; z < pl->GetNp(); z++)
                {
                    if(pl->GetPa(z) == i)
                    {
                        found = 1;
                        break;
                    }
                }
                if (found)
                {
                    for(z = 0; z < pl->GetNp(); z++)
                    {
                        k = pl->GetPa(z); 
                        if(k != i)
                        { 
                            (node+k)->SetFlag(1);
                            (node+k)->SetNSons((node+k)->GetNSons()+1);
                        }
                    }
                }
            }
        }
    }

    matij = tmpmatrix->GetStart(i);
    while(matij.GetNext())
    {
        j = matij.GetIndex();
        nodej = node+j;
        if(((nodej)->IsCGNode()) || ((nodej)->IsFGNode())) continue;
        {
            for(pl = nodej->GetPaList(); pl != NULL; pl = pl->GetNext())
            {
                found = 0;
                for(z = 0; z < pl->GetNp(); z++)
                {
                    if(pl->GetPa(z) == i)
                    {
                        found = 1;
                        break;
                    }
                }
                if (found)
                {
                    for(z = 0; z < pl->GetNp(); z++)
                    {
                        k = pl->GetPa(z); 
                        if((node+k)->GetFlag() == 1)
                        { 
                            (node+k)->SetFlag(0);
                            UpdateNBNewCG(k);
                        }
                    }
                }
            }
        }
    }
                    

    UpdateNBNewCG(i);

    return 0;
}

 
void CMGPaList::MarkParents(CMGGrid *grid)
{
    CMGNode *node,*cgnode;
    CMGGraph *graph;
    int i;
    
    graph = grid->GetGraph();
    node = graph->GetNode();
    for(i = 0; i < np; i++)
    {
        cgnode = node + pa[i];
        if(cgnode->IsCGNode()) continue;
        graph->Remove(cgnode);
        graph->MarkCGNode(cgnode);
        graph->UpdateNSons(NULL,cgnode->GetPaList(),grid);
        graph->ClearPaList(cgnode->GetPaList());
        cgnode->SetPaList(NULL);
        grid->UpdateNeighborsCG(pa[i]);
    }

    return;
}

        
int CMGGrid::SaveCoeffs(int i, int np, int *pa, double *coeff, double *coefft)
{
    CMGTransferEntry *trans;
    int j,z;

    trans = transfer->GetRow();
    for(z = 0; z < np; z++)
    {
        j = pa[z];
        if((trans+i)->SaveEntry(trans+j,coeff[z])) return 1;
        if((trans+j)->SaveEntry(trans+i,coefft[z])) return 1;
    }
 
    return 0;
}

int CMGNode::Eliminate(CMGGrid *grid)
{
    CMGGraph *graph;
    CMGPaList *pl, *minpl;
    double weight, minweight;

    graph = grid->GetGraph();

    minweight = 1e+10;
    for(pl = palist; pl != NULL; pl = pl->GetNext())
    {
        weight = pl->TotalWeight();
        if (weight < minweight)
        {
            minpl = pl;
            minweight = weight;
        }
    }
    
    graph->MarkFGNode(this);
    graph->UpdateNSons(NULL,palist,grid);
    graph->ClearPaList(palist);
    palist = NULL;
    minpl->MarkParents(grid);

    if(grid->SaveCoeffs(id,minpl->GetNp(),minpl->GetPa(),minpl->GetCoeff(),minpl->GetCoefft())) return 1;
   

    return 0;
}

int CMGNode::UpdateNeighborsFG(CMGGrid *grid)
{
    CMGGraph *graph;
    CMGNode *node, *nodej;
    CMGMatrixPtr matij;
    CMGPaList *palist;
    int j;
    
    graph = grid->GetGraph();
    node = graph->GetNode();

    // check all current neighbors, i.e. all neighbors and the
    // parent nodes of the FG neighbors. The parent nodes are marked,
    // thus tmpmatrix or (matrix) is fine.
    matij = grid->GetTmpMatrix()->GetStart(id);
    while(matij.GetNext())
    {
        // maybe not enough
        j = matij.GetIndex();
        nodej = node+j;
        if(((nodej)->IsCGNode()) || ((nodej)->IsFGNode())) continue;
        graph->Store(nodej);
        if(grid->AnalyseNodeSimple(j,palist)) return 1;
        graph->UpdateNSons(palist,nodej->GetPaList(),grid);
        graph->ClearPaList(nodej->GetPaList());
        nodej->SetPaList(palist);
        nodej->CountNewLinks(grid, graph);
        nodej->CountNewCG(graph);
    }

    return 0;
}
        

          
int CMGGraph::RemainingNodes(CMGGrid *gridptr)
{
    CMGNode* cgnode;
    
    cgnode = GetFirstNode();
    while(cgnode != NULL)
    {        
        // cannot be eliminated
        MarkCGNode(cgnode);
        // actually graph->UpdateNSons(NULL,NULL); (does nothing)

        // maybe not necessary
        gridptr->UpdateNeighborsCG(cgnode->GetId());

        cgnode = GetFirstNode();

        /* debug 
        PICTURE *thePic;
        char c;
        if(n < 0)
        {
            thePic = GetCurrentPicture();
            if (thePic!=NULL)
            {
                DrawUgPicture(thePic);
            } 
            cin >> c;
            } */
    }

    return 0;
}
int CMGGraph::EliminateNodes(CMGGrid *gridptr)
{
    CMGNode* fgnode;
    
    fgnode = GetFirstNode();
    while(fgnode != NULL)
    {        
        if(fgnode->GetPaList() == NULL)
        {
            // done
            if(Insert(fgnode)) return 1;
            return 0;
        }
        else
        {
            if(fgnode->Eliminate(gridptr)) return 1;
            if(fgnode->UpdateNeighborsFG(gridptr)) return 1; 
            if(InsertHelplist()) return 1;
        }
        fgnode = GetFirstNode();

        /* debug 
        PICTURE *thePic;
        char c;
        if(n < 10000)
        {
            thePic = GetCurrentPicture();
            if (thePic!=NULL)
            {
                DrawUgPicture(thePic);
            } 
            cin >> c;
            } */
    }

    return 0;
}


int CMGGraph::Construct(CMGGrid *gridptr)
{
    CMGNode *nodei;
    CMGPaList *palist;
    int i;

#ifdef UG_DRAW
    PICTURE *thePic;
    if(n > 0)
    {
        thePic = GetCurrentPicture();
        if (thePic!=NULL)
        {
            DrawUgPicture(thePic);
        } 
    } 
#endif 

    int type = CMGGetParameter()->Gettype();
    
    for(i = 0; i < n; i++)
    {
        nodei = node+i;
        if(nodei->IsFGNode()) continue;
        switch (type)
        {
        case 0: if (gridptr->AnalyseNode0(i,palist)) return 1; break;
        case 1: if (gridptr->AnalyseNode1(i,palist)) return 1; break;
        case 2: if (gridptr->AnalyseNode2(i,palist)) return 1; break;
        case 3: if (gridptr->AnalyseNode3(i,palist)) return 1; break;
        case 4: if (gridptr->AnalyseNode4(i,palist)) return 1; break;
        case 5: if (gridptr->AnalyseNode5(i,palist)) return 1; break;
        }
        nodei->SetPaList(palist);
    }
    InitNSons();
    for(i = 0; i < n; i++)
    {
        nodei = node+i;
        if(nodei->IsFGNode()) continue;
        Remove(nodei);  // might be inserted by UpdateNeighborsCG
        // if(nodei->GetPaList() == NULL)
        // {
        //     MarkCGNode(nodei);
        //     gridptr->UpdateNeighborsCG(i);
        // }
        // else
        {
            nodei->CountNewLinks(gridptr, this);
            nodei->CountNewCG(this);
            nodei->ComputeTotalWeight();
            if(Insert(nodei)) return 1;
        }
    }

    return 0;
}

int CMGGraph::Construct2(CMGGrid *gridptr)
{
    CMGNode *nodei;
    CMGPaList *palist;
    int i;

#ifdef UG_DRAW
    /* debug */
    PICTURE *thePic;
    if(n > 0)
    {
        thePic = GetCurrentPicture();
        if (thePic!=NULL)
        {
            DrawUgPicture(thePic);
        } 
    } 
#endif

    int type = CMGGetParameter()->Gettype();
    
    for(i = 0; i < n; i++)
    {
        nodei = node+i;
        if((nodei->IsFGNode()) || (nodei->IsCGNode())) continue; 
        switch (type)
        {
        case 0: if (gridptr->AnalyseNode0(i,palist)) return 1; break;
        case 1: if (gridptr->AnalyseNode1(i,palist)) return 1; break;
        case 2: if (gridptr->AnalyseNode2(i,palist)) return 1; break;
        case 3: if (gridptr->AnalyseNode3(i,palist)) return 1; break;
        case 4: if (gridptr->AnalyseNode4(i,palist)) return 1; break;
        case 5: if (gridptr->AnalyseNode5(i,palist)) return 1; break;
        }
        nodei->SetPaList(palist);
        nodei->SetNSons(0);
                                 
    }

    InitNSons();
    for(i = 0; i < n; i++)
    {
        nodei = node+i;
        if((nodei->IsFGNode()) || (nodei->IsCGNode())) continue; 
        Remove(nodei);  
        // if(nodei->GetPaList() == NULL)
        // {
        //     MarkCGNode(nodei);
        //     gridptr->UpdateNeighborsCG(i);
        // }
        // else
        {
            nodei->CountNewLinks(gridptr, this);
            nodei->CountNewCG(this);
            nodei->ComputeTotalWeight();
            if(Insert(nodei)) return 1;
        }
    }

    return 0;
}









































