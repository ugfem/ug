/****************************************************************************/
/*																			*/
/* File:      famg_construct.C												*/
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
#include "famg_system.h"

#ifdef Debug
	#ifdef USE_UG_DS
		#include "famg_uginterface.h"
	#else
		#include "famg_interface.h"
	#endif
#endif

#ifdef USE_UG_DS
extern "C"
{
#include "parallel.h"
#include "commands.h" /* for GetCurrentMultigrid for debuggung */
}
#endif

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

void prm(int level, int comp)
// for debugging
{
#ifdef USE_UG_DS
	VECTOR *v;
	MATRIX *m;
	FAMGTransferEntry *trans;
	
	GRID *g = GRID_ON_LEVEL(GetCurrentMultigrid(), level);

	printf("Matrix on level %d: component %d\n", GLEVEL(g), comp);
	for (v=PFIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
	{
#ifdef ModelP 
		printf(PFMT"(P%d)",me,PRIO(v));
#endif		
		printf("v[%4d] ", VINDEX(v));
		for (m=VSTART(v); m!=NULL; m = MNEXT(m))
		{
			printf("\t%g(->%d)",MVALUE(m,comp),VINDEX(MDEST(m)));
		}
		printf("\n");
	}
#else
	FAMGSystem &sys = *FAMG_GetSystem();
	FAMGGrid &grid = *sys.GetMultiGrid(0)->GetGrid(level);

	FAMGMatrix &M = *grid.GetMatrix();
	int n = grid.GetN(), end;
    
	FAMGMatrixPtr matij;

	printf("Matrix:\n");
    for(int i = 0; i < n; i++)
    {
		printf("v[%4d] ", i);
		end = 1;
		for( matij=M.GetStart(i); end; end=matij.GetNext() )
		{
			printf("\t%g(->%d)", matij.GetData(),matij.GetIndex());
		}
		printf("\n");
	}
#endif
}

void prim(int level)
// for debugging
{
#ifdef USE_UG_DS
	VECTOR *v;
	MATRIX *m;
	FAMGTransferEntry *trans;
	
	GRID *g = GRID_ON_LEVEL(GetCurrentMultigrid(), level);

	printf("Interpolation Matrix on level %d:\n", GLEVEL(g));
	for (v=PFIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
	{
#ifdef ModelP 
		printf(PFMT"(P%d)",me,PRIO(v));
#endif		
		printf("v[%4d] %c ", VINDEX(v), VCCOARSE(v)?'C':'F' );
		for (m=VISTART(v); m!=NULL; m = MNEXT(m))
		{
			trans = (FAMGTransferEntry*)m;
			printf("\tP=%g(->%d) R=%g", trans->GetProlongation(),VINDEX(MDEST(m)),trans->GetRestriction());
		}
		printf("\n");
	}
#else
	FAMGSystem &sys = *FAMG_GetSystem();
	FAMGGrid &grid = *sys.GetMultiGrid(0)->GetGrid(level);
	FAMGTransferEntry *transij, *resji;
	FAMGTransfer *trans;
	FAMGMatrix *matrix;
	int n = grid.GetN(), i, j;
    double resval;
	
	trans = grid.GetTransfer();
	matrix = grid.GetMatrix();
	printf("Interpolation Matrix:\n");
    for(i = 0; i < n; i++)
    {
		printf("v[%4d] ", i);
		if (matrix->GetType(i)) // FG Node
		{
			for(transij = trans->GetRow(i)->GetNext(); transij != NULL; transij = transij->GetNext())
			{
				j = transij->GetId();
				resval = -99.99;
				for(resji = trans->GetRow(j)->GetNext(); resji != NULL; resji = resji->GetNext())
				{
					if( resji->GetId() == i )
					{
						resval = resji->GetData();
							break;
					}
				}
				printf("\tP=%g(->%d) R=%g", transij->GetData(),j,resval);
			}
		}
		else
		{
			printf("\tP=%g R=%g", 1.0, 1.0);
		}
		printf("\n");
	}
#endif
}

void prv( int level, int x_nr )
/* for calling from a debugger */
{
#ifdef USE_UG_DS
	register VECTOR *v;
	DOUBLE pos[DIM];
	
	GRID *g = GRID_ON_LEVEL(GetCurrentMultigrid(), level);
	
	printf("Vector on level %d:\n",GLEVEL(g));
    for (v=PFIRSTVECTOR(g); v!= NULL; v=SUCCVC(v))
    {
		VectorPosition(v,pos);
#ifdef ModelP 
		printf(PFMT"(P%d)",me,PRIO(v));
#endif
		if(VOBJECT(v)==NULL)
		{
			printf("x= ---- y= ---- ");
#ifdef __THREEDIM__
			printf("z= ---- ");
#endif
		}
		else
		{
			printf("x=%5.2f y=%5.2f ",pos[0],pos[1]);
#ifdef __THREEDIM__
			printf("z=%5.2f ",pos[2]);
#endif
		}
		printf("  index = %d %c ", VINDEX( v ) , VCCOARSE(v)?'C':'F' );
		printf("u[%d]=%15.8f GID %08x",x_nr,VVALUE(v,x_nr), GID(v));
		/*printf("   cl %d %d sk ",VCLASS(v),VNCLASS(v));*/
		/*for (j=0; j<ncomp; j++)
			printf("%d ",((VECSKIP(v) & (1<<j))!=0));*/
		printf("\n");
	}
	return;
#else
	FAMGSystem &sys = *FAMG_GetSystem();
	FAMGGrid &grid = *sys.GetMultiGrid(0)->GetGrid(level);

	double *vec = grid.GetVector(x_nr);
	int n = grid.GetN();
    
	printf("Vector:\n");
    for(int i = 0; i < n; i++)
    {
		printf("vec[%4d] = %g\n", i, vec[i]);
	}
#endif
}

int FAMGGrid::AnalyseNodeSimple(FAMGNode* nodei, FAMGPaList *&palist)
{
    FAMGPaList *pl;
    int remove, z;

    palist = NULL;
    for(pl = nodei->GetPaList(); pl != NULL; pl = pl->GetNext())
    {
        remove = 0;
        for(z = 0; z < pl->GetNp(); z++)
        {
            if(graph->GetNode(pl->GetPa(z))->IsFGNode()) 
            {
                remove = 1;
                break;
            }
        }
        if(!remove)
        {
            if(graph->SavePaList(palist,pl->GetNp(),pl->GetPaPtr(),pl->GetCoeff(),pl->GetCoefft(),pl->GetApprox())) 
				return 1;
        }
    }

    return 0;
}


double FAMGPaList::TotalWeight()
{
    double weight;

    weight = 0.0*approx;
    weight += 1.0*(double) newlinks;
    weight += 10.0*(double) newcg;

    return weight;
} 


void FAMGNode::ComputeTotalWeight()
{
    FAMGPaList *pl;
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


void FAMGNode::CountNewCG(FAMGGraph *graph)
{
    FAMGNode *node;
    FAMGPaList *pl;
    double nc;
    int i, np, ns;

    node = graph->GetNodePtr();
    for(pl  = GetPaList(); pl != NULL; pl = pl->GetNext())
    {
        np = pl->GetNp();
        const int *pa = pl->GetPaPtr();
        nc = 0;
        for(i = 0; i < np; i++)
        {
            if(node[pa[i]].IsCGNode()) continue;
            else
            {
                ns = node[pa[i]].GetNSons();
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

int FAMGNode::CountNewCG(FAMGGraph *graph, int j)
{
    FAMGNode *node;
    FAMGPaList *pl;
    double nc;
    int i, np, found, change, ns;

    node = graph->GetNodePtr();
    change = 0;
    for(pl  = GetPaList(); pl != NULL; pl = pl->GetNext())
    {
        np = pl->GetNp();
        const int *pa = pl->GetPaPtr();
        nc = 0;
        found = 0;
        for(i = 0; i < np; i++)
        {
            if(pa[i] == j) 
				found = 1;
        }
        if(found)
        {
            for(i = 0; i < np; i++)
            {
                if(node[pa[i]].IsCGNode()) 
					continue;
                else
                {
                    ns = node[pa[i]].GetNSons();
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

// i,z are node indices: may be coarse or non-determined, but are NOT allowed to be fine
#define FAMG_TRAVERSE_RETURN_TYPE		int
#define FAMG_TRAVERSE_FCT_NAME			Connected
#define FAMG_TRAVERSE_SECOND_ARG		int z
#define FAMG_TRAVERSE_INIT(node_i)		FAMGNode *node_z = GetGraph()->GetNode(z);
#define FAMG_TRAVERSE_ACTION(node)		if((node)==node_z) return 1;
//#define FAMG_TRAVERSE_FINISH(node_i)
#define FAMG_TRAVERSE_RETURN			0
#include "famg_traverse.template"

// count the number of non-fine neighbors and set their flag=1 (with exception of i itself)
// i is node index: may be coarse or non-determined, but are NOT allowed to be fine
#define FAMG_TRAVERSE_RETURN_TYPE		int
#define FAMG_TRAVERSE_FCT_NAME			SetFlagsAndCount
#define FAMG_TRAVERSE_SECOND_ARG		int flag
#define FAMG_TRAVERSE_INIT(node_i)		int z = 0; (node_i)->SetFlag(flag);
#define FAMG_TRAVERSE_ACTION(node)		if((node)->GetFlag() != flag) {z++;(node)->SetFlag(flag);}
#define FAMG_TRAVERSE_FINISH(node_i) 	(node_i)->SetFlag(0);
#define FAMG_TRAVERSE_RETURN			z
#include "famg_traverse.template"

// for the non-fine neighbors set the flag=1 (with exception of i itself)
// i is node index: may be coarse or non-determined, but are NOT allowed to be fine
#define FAMG_TRAVERSE_RETURN_TYPE		void
#define FAMG_TRAVERSE_FCT_NAME			SetFlags
#define FAMG_TRAVERSE_SECOND_ARG		int flag
//#define FAMG_TRAVERSE_INIT(node_i)
#define FAMG_TRAVERSE_ACTION(node)		(node)->SetFlag(flag);
#define FAMG_TRAVERSE_FINISH(node_i) 	(node_i)->SetFlag(0);
//#define FAMG_TRAVERSE_RETURN
#include "famg_traverse.template"

#define FAMG_TRAVERSE_RETURN_TYPE		int
#define FAMG_TRAVERSE_FCT_NAME			CountLinks
//#define FAMG_TRAVERSE_SECOND_ARG
#define FAMG_TRAVERSE_INIT(node_i)		int z = 0;
#define FAMG_TRAVERSE_ACTION(node)		if((node)->GetFlag() == 1) {z++;(node)->SetFlag(0);}
//#define FAMG_TRAVERSE_FINISH(node_i)
#define FAMG_TRAVERSE_RETURN			z
#include "famg_traverse.template"


int FAMGNode::CountNewLinks(FAMGGrid *gridptr, FAMGGraph *graph)
{
    FAMGPaList *pl;
    int nnb, nl, z, np, y, newlinks;
    
    nnb = gridptr->SetFlagsAndCount(GetId(),1);
    for(pl = palist; pl != NULL; pl = pl->GetNext())
    {
        np = pl->GetNp();
        const int *pa = pl->GetPaPtr();
        nl = 0;
        for(z = 0; z < np; z++)
        {
            graph->GetNode(pa[z])->SetFlag(0);
            nl += gridptr->CountLinks(pa[z]);
            gridptr->SetFlags(GetId(),1);
            for(y = z+1; y < np; y++)
            {
                if(!gridptr->Connected(pa[z],pa[y])) 
					nl++;
            }
        }
        newlinks = np*(nnb-1)-nl;
        if(np < 0) 
			newlinks = 20;
        pl->SetNewLinks(newlinks);
    }
    gridptr->SetFlags(GetId(),0);

    return 0;
}


void FAMGGraph::InitNSons()
{
    FAMGPaList *palist, *pl;
	FAMGNode *node_j;
    int i, np, j, z, mark;

    for(i = 0; i < GetN(); i++)
    {
        palist = GetNode(i)->GetPaList();
        for(pl = palist; pl != NULL; pl = pl->GetNext())
        {
            np = pl->GetNp();
            mark = 0;
            for(z = 0; z < np; z++)
            {
                j = pl->GetPa(z);
                if(GetNode(j)->IsCGNode())
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
					node_j = GetNode(j);
                    if(!node_j->IsCGNode())
                    {
                        node_j->SetFlag(1);
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
				node_j = GetNode(j);
                if(node_j->GetFlag() == 1)
                {
                    node_j->SetFlag(0);
                    node_j->SetNSons(node_j->GetNSons()+1);
                }
            }
        }
    }

    return;
}
            
int FAMGGrid::UpdateNBNewCG(int i)
{
    FAMGNode *nodej;

    // check all neighbors with i as possible father, tmpmatrix is OK
	FAMGMatrixEntry matij;
	FAMGVectorEntry vec_i = GetGraph()->GetNode(i)->GetVec();
	FAMGMatrixIter matiter(*GetTmpMatrix(), vec_i);

	matiter(matij);		// skip diagonal entry
    while(matiter(matij))
    {
		nodej = GetGraph()->GetNode(matij.dest());
        if(nodej->IsCGNode() || nodej->IsFGNode()) 
			continue;
        if(nodej->CountNewCG(GetGraph(),i))
        {
            // something changed
            GetGraph()->Store(nodej);
        }
    }

    return 0;
}

void FAMGGraph::UpdateNSons(FAMGPaList *newlist, FAMGPaList *oldlist, FAMGGrid* grid)
{
    FAMGPaList *pl;
    FAMGNode *nodej;
    int np, z, j, mark, fmark;
    
    for(pl = oldlist; pl != NULL; pl = pl->GetNext())
    {
        np = pl->GetNp();
        mark = 0; fmark = 1;
        for(z = 0; z < np; z++)
        {
            j = pl->GetPa(z);
			nodej = GetNode(j);
            if(nodej->IsCGNode()) 
				mark = 1;
            if(nodej->IsFGNode()) 
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
				nodej = GetNode(j);
                if(!nodej->IsCGNode())
                {
                    nodej->SetFlag(1);
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
			nodej = GetNode(j);
            if(nodej->IsCGNode()) 
				mark = 1;
            if(nodej->IsFGNode()) 
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
				nodej = GetNode(j);
                if(!nodej->IsCGNode())
                {
                    nodej->SetFlag1(1);
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
			nodej = GetNode(j);
            if(nodej->GetFlag()
               && !nodej->GetFlag1()
               && !nodej->GetFlag2())
            {
                nodej->SetFlag2(1);
                nodej->SetNSons(nodej->GetNSons()-1);
            }
        }
    }

    for(pl = newlist; pl != NULL; pl = pl->GetNext())
    {
        np = pl->GetNp();
        for(z = 0; z < np; z++)
        {
            j = pl->GetPa(z);
			nodej = GetNode(j);
            if(!nodej->GetFlag()
               && nodej->GetFlag1()
               && !nodej->GetFlag2())
            {
               nodej->SetFlag2(1);
               nodej->SetNSons(nodej->GetNSons()+1);
            }
        }
    }

    for(pl = oldlist; pl != NULL; pl = pl->GetNext())
    {
        np = pl->GetNp();
        for(z = 0; z < np; z++)
        {
            j = pl->GetPa(z);
			nodej = GetNode(j);
            nodej->SetFlag(0);
            nodej->SetFlag1(0);
            if(nodej->GetFlag2())
            {
                nodej->SetFlag2(0);
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
			nodej = GetNode(j);
            nodej->SetFlag(0);
            nodej->SetFlag1(0);
            if(nodej->GetFlag2())
            {
                nodej->SetFlag2(0);
                grid->UpdateNBNewCG(j);
            }
        }
    }

    return;
}


int FAMGGrid::UpdateNeighborsCG(int i)
{
	FAMGNode *nodej, *nodek;
	FAMGPaList *pl;
	int z, found, k;

	// check all neighbors with i as possible father, tmpmatrix is OK
	FAMGMatrixEntry matij;
	FAMGVectorEntry vec_i = GetGraph()->GetNode(i)->GetVec();
	FAMGMatrixIter matiter(*GetTmpMatrix(), vec_i);
	
	matiter(matij);		// skip diagonal entry
	while(matiter(matij))
	{
		nodej = GetGraph()->GetNode(matij.dest());
        if(nodej->IsCGNode() || nodej->IsFGNode())
			continue;

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
						nodek = GetGraph()->GetNode(k);
						nodek->SetFlag(1);
						nodek->SetNSons(nodek->GetNSons()+1);
					}
				}
			}
		}

	}

	matiter.reset();
	matiter(matij);		// skip diagonal entry
	while(matiter(matij))
	{
		nodej = GetGraph()->GetNode(matij.dest());
        if(nodej->IsCGNode() || nodej->IsFGNode())
			continue;
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
					nodek = GetGraph()->GetNode(k); 
					if(nodek->GetFlag() == 1)
					{ 
						nodek->SetFlag(0);
						UpdateNBNewCG(k);
					}
				}
			}
		}
	}

	UpdateNBNewCG(i);

	return 0;
}

 
void FAMGPaList::MarkParents(FAMGGrid *grid)
{
	FAMGGraph *graph = grid->GetGraph();
    FAMGNode *cgnode;
    int i;
    
    for(i = 0; i < np; i++)
    {
        cgnode = graph->GetNode(pa[i]);
        if(cgnode->IsCGNode()) 
			continue;
        graph->Remove(cgnode);
        graph->MarkCGNode(cgnode);
        graph->UpdateNSons(NULL,cgnode->GetPaList(),grid);
        graph->ClearPaList(cgnode->GetPaList());
        cgnode->SetPaList(NULL);
        grid->UpdateNeighborsCG(pa[i]);
    }

    return;
}

        
int FAMGGrid::SaveCoeffs(const FAMGVectorEntry& fg_vec, int np, const int pa[], double coeff[], double coefft[])
{
    int z;

    for(z = 0; z < np; z++)
		transfer->SetEntries(fg_vec,GetGraph()->GetNode(pa[z])->GetVec(),coeff[z],coefft[z]);
 
    return 0;
}

int FAMGNode::CheckPaList(FAMGGraph *graph)
{
    FAMGPaList *pl, *ppl, *opl;
    int update, remove, z, j;
	FAMGNode *nodej;
	
    update = 0;
    pl = palist;
    ppl = NULL;
    while(pl != NULL)
    {
        opl = pl;
        pl = pl->GetNext();

        remove = 0;
        for(z = 0; z < opl->GetNp(); z++)
        {
            j = opl->GetPa(z);
			nodej = graph->GetNode(j);
            if(nodej->IsFGNode()) 
            {
                remove = 1;
                break;
            }
        }
        if(remove)
        {
            if(ppl != NULL) 
				ppl->SetNext(pl);
            else 
				palist = pl;

            // store opl in freelist
            opl->SetNext(graph->GetFreePaList());
            graph->SetFreePaList(opl);
            update = 1;
        }
        else
        {
            ppl = opl;
        }
        
    }

    return update;
}

int FAMGNode::Eliminate(FAMGGrid *grid)
{
	FAMGGraph *graph = grid->GetGraph();
    FAMGPaList *pl, *minpl = NULL;
    double weight, minweight;

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

    if(grid->SaveCoeffs(GetVec(),minpl->GetNp(),minpl->GetPaPtr(),minpl->GetCoeff(),minpl->GetCoefft())) return 1;

    return 0;
}

int FAMGNode::UpdateNeighborsFG(FAMGGrid *grid)
{
	FAMGGraph *graph = grid->GetGraph();
	FAMGNode *nodej;
	FAMGPaList *palist;
    
	// check all current neighbors, i.e. all neighbors and the
	// parent nodes of the FG neighbors. The parent nodes are marked,
	// thus tmpmatrix or (matrix) is fine.
	
	FAMGMatrixEntry matij;
	FAMGMatrixIter matiter(*(grid->GetTmpMatrix()), GetVec());

	matiter(matij);		// skip diagonal entry
    while(matiter(matij))
    {
		// maybe not enough
		nodej = graph->GetNode(matij.dest());
		if(nodej->IsCGNode() || nodej->IsFGNode()) 
			continue;
		graph->Store(nodej);
		if(grid->AnalyseNodeSimple(nodej,palist)) 
			return 1;
		graph->UpdateNSons(palist,nodej->GetPaList(),grid);
		graph->ClearPaList(nodej->GetPaList());
		nodej->SetPaList(palist);
		nodej->CountNewLinks(grid, graph);
		nodej->CountNewCG(graph);
	}

	return 0;
}
        

int FAMGGraph::RemainingNodes(void)
{
    FAMGNode* cgnode;
    
    cgnode = GetFirstNode();
    while(cgnode != NULL)
    {        
        // cannot be eliminated
        MarkCGNode(cgnode);

        cgnode = GetFirstNode();

    }

    return 0;
}
          
int FAMGGraph::EliminateNodes(FAMGGrid *gridptr)
{
    FAMGNode* fgnode;
    
    fgnode = GetFirstNode();
    while(fgnode != NULL)
    {        
        if(fgnode->GetPaList() == NULL)
        {
            // done
            if(Insert(fgnode)) return 1;
            return 0;
        }
        else if (fgnode->CheckPaList(this))
        {
            // CheckPaList is necessary for non structure
            //   symmetric matrices. It is not nice and
            //   should be avoided. 
            fgnode->ComputeTotalWeight();
            if(Insert(fgnode)) return 1;
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
            } 
		*/
    }

    return 0;
}

int FAMGGraph::InsertNode(FAMGGrid *gridptr, FAMGNode *nodei)
{
	if(nodei->IsFGNode()) 
		return 0;
	Remove(nodei);  // might be inserted by UpdateNeighborsCG
	nodei->CountNewLinks(gridptr, this);
	nodei->CountNewCG(this);
	nodei->ComputeTotalWeight();
	if(Insert(nodei))
		return 1;
	return 0;
}
		
int FAMGGraph::Construct(FAMGGrid *gridptr)
{
	FAMGGraph *graph = gridptr->GetGraph();
    FAMGNode *nodei;
    FAMGPaList *palist;
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

    int type = FAMGGetParameter()->Gettype();
    
    for(i = 0; i < n; i++)
    {
        nodei = graph->GetNode(i);
        if(nodei->IsFGNode()) 
			continue;
        switch (type)
        {
        case 0: if (gridptr->AnalyseNode0(nodei->GetVec(),palist)) return 1; break;
        case 1: if (gridptr->AnalyseNode1(nodei->GetVec(),palist)) return 1; break;
        case 2: if (gridptr->AnalyseNode2(nodei->GetVec(),palist)) return 1; break;
        case 3: if (gridptr->AnalyseNode3(nodei->GetVec(),palist)) return 1; break;
        case 4: if (gridptr->AnalyseNode4(nodei->GetVec(),palist)) return 1; break;
        case 5: if (gridptr->AnalyseNode5(nodei->GetVec(),palist)) return 1; break;
        }
        nodei->SetPaList(palist);
    }
    InitNSons();
	
#ifdef ModelP
	// in the first step eliminate only nodes in the border of the core partition
	// the remaining nodes must be processed in a second step (see FAMGGrid:ConstructTransfer)
	VECTOR *vec;
	MATRIX *mat;
	
    for(i = 0; i < n; i++)
    {
        nodei = graph->GetNode(i);
		vec = ((FAMGugVectorEntryRef*)(nodei->GetVec().GetPointer()))->myvector();
		
		if( IS_FAMG_GHOST(vec) )
			continue; // only master vectors can be in border the of the core partition
		
		// vec lies in the border of the core partition if he has a ghost or border neighbor
		for( mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat) )
			if( IS_FAMG_GHOST(MDEST(mat)) )
				break;
		
		if( mat != NULL )
			// a ghost neighbor was found
			if(InsertNode(gridptr, nodei))
				return 0;
    }
#else
#ifdef SIMULATE_HALFENING	// TODO: remove it
	VECTOR *vec;
	MATRIX *mat;
	DOUBLE pos[3];
	
    for(i = 0; i < n; i++)
    {
        nodei = graph->GetNode(i);
		vec = ((FAMGugVectorEntryRef*)(nodei->GetVec().GetPointer()))->myvector();
		
		VectorPosition(vec,pos);

		if( fabs(pos[0]-0.5)<1e-3 )	// insert the middle col of the domain into the list
			if(InsertNode(gridptr, nodei))
				return 0;
	}
#else
	// put all nodes into the list
    for(i = 0; i < n; i++)
    {
        nodei = graph->GetNode(i);
		if(InsertNode(gridptr, nodei))
			return 0;
    }
#endif
#endif
	
    return 0;
}

int FAMGGraph::Construct2(FAMGGrid *gridptr)
{
	FAMGGraph *graph = gridptr->GetGraph();
    FAMGNode *nodei;
    FAMGPaList *palist;
    int i;

#ifdef ModelP
	assert(0); // should not be called for parallel
#endif
	
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

    int type = FAMGGetParameter()->Gettype();
    
    for(i = 0; i < n; i++)
    {
        nodei = graph->GetNode(i);
        if(nodei->IsFGNode() || nodei->IsCGNode()) 
			continue; 
        switch (type)
        {
        case 0: if (gridptr->AnalyseNode0(nodei->GetVec(),palist)) return 1; break;
        case 1: if (gridptr->AnalyseNode1(nodei->GetVec(),palist)) return 1; break;
        case 2: if (gridptr->AnalyseNode2(nodei->GetVec(),palist)) return 1; break;
        case 3: if (gridptr->AnalyseNode3(nodei->GetVec(),palist)) return 1; break;
        case 4: if (gridptr->AnalyseNode4(nodei->GetVec(),palist)) return 1; break;
        case 5: if (gridptr->AnalyseNode5(nodei->GetVec(),palist)) return 1; break;
        }
        nodei->SetPaList(palist);
        nodei->SetNSons(0);
                                 
    }

    InitNSons();
    for(i = 0; i < n; i++)
    {
        nodei = graph->GetNode(i);
        if(nodei->IsFGNode() || nodei->IsCGNode())
			continue; 
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
