// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/****************************************************************************/
/*																			*/
/* File:	  handler.c														*/
/*																			*/
/* Purpose:   defines the handlers used by ddd during data management.      */
/*			  																*/
/*																			*/
/* Author:	  Stefan Lang                       				 			*/
/*			  Institut fuer Computeranwendungen III 						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: stefan@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   10.05.95 begin, ugp version 3.0								*/
/*																			*/
/* Remarks: 																*/
/*																			*/
/****************************************************************************/

/* TODO: delete this */
/* #include "conf.h" */

#ifdef ModelP

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files 									*/
/*																			*/
/****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "debug.h"
#include "compiler.h"
#include "domain.h"
#include "parallel.h"
#include "heaps.h"
#include "ugm.h"
#include "algebra.h"
#include "general.h"
#include "rm.h"
#include "refine.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/



/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/



/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)



/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


static GRID *GetGridOnDemand (MULTIGRID *mg, int level)
{
	while (level > TOPLEVEL(mg))
	{
		CreateNewLevel(mg);
		UserWriteF("CreateNewLevel %d", TOPLEVEL(mg));
		/* TODO error handling, CreateNewLevel() may return NULL */
	}

	return GRID_ON_LEVEL(mg,level);
}




/****************************************************************************/
/*																			*/
/*  	For data management during redistribution and communication 		*/
/*		DDD needs for each DDD (data) object several handlers.				*/
/*		These are:															*/
/*			HANDLER_LDATACONSTRUCTOR,  handler: init object's LDATA parts   */
/*			HANDLER_UPDATE,            handler: update objects internals    */
/*			HANDLER_OBJMKCONS,         handler: make obj consistent         */
/*			HANDLER_DESTRUCTOR,        handler: destruct object             */
/*			HANDLER_XFERCOPY,          handler: copy cmd during xfer        */
/*			HANDLER_XFERDELETE,        handler: delete cmd during xfer      */
/*			HANDLER_XFERGATHER,        handler: send additional data        */
/*			HANDLER_XFERSCATTER,       handler: recv additional data        */
/*																			*/
/*																			*/
/*	    For each of the ddd (data) objects the handlers needed for the      */
/*	    specific object are defined below in the order handlers for         */
/*																			*/
/*			DDD objects:													*/
/*				*dimension independent										*/
/*				 TypeVector, TypeIVertex, TypeBVertex, TypeNode				*/
/*				*dimension dependent										*/
/*				 2-Dim:														*/
/*				 TypeTrElem, TypeTrBElem,									*/
/*				 TypeQuElem, TypeQuBElem									*/
/*				 3-Dim:														*/
/*				 TypeTeElem, TypeTeBElem									*/
/*				 TypePyElem, TypePyBElem									*/
/*				 TypeHeElem, TypeHeBElem									*/
/*																			*/
/*			DDD data objects:												*/
/*				TypeMatrix, 												*/
/*				TypeEdge													*/
/*																			*/
/*		NOT all handlers are to be specified for each object!				*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/****************************************************************************/
/*																			*/
/*		handlers for typevector												*/
/*																			*/
/****************************************************************************/
/****************************************************************************/

void VectorUpdate (DDD_OBJ obj)
{
	VECTOR *pv = (VECTOR *)obj;
	VECTOR *after = NULL;
	GRID   *theGrid = NULL;
	int    level = DDD_InfoAttr(PARHDR(pv));
	int    prio = DDD_InfoPriority(PARHDR(pv));

	PRINTDEBUG(dddif,1,("%2d: VectorUpdate(): v=%08x/%x VEOBJ=%d\n",me,\
		DDD_InfoGlobalId(PARHDR(pv)),pv,OBJT(pv)))

	theGrid = GRID_ON_LEVEL(dddctrl.currMG,level);
	after = LASTVECTOR(theGrid);

    /* insert in vector list */
	GRID_LINK_VECTOR(theGrid,pv,prio)

/* TODO: delete this */
if (0) {
    if (after==NULL)
    {
        SUCCVC(pv) = (VECTOR*)FIRSTVECTOR(theGrid);
        PREDVC(pv) = NULL;
        if (SUCCVC(pv)!=NULL)
            PREDVC(SUCCVC(pv)) = pv;
        SFIRSTVECTOR(theGrid) = (void*)pv;
        if (LASTVECTOR(theGrid)==NULL)
            LASTVECTOR(theGrid) = (void*)pv;
    }
    else
    {
        SUCCVC(pv) = SUCCVC(after);
        PREDVC(pv) = after;
        if (SUCCVC(pv)!=NULL)
            PREDVC(SUCCVC(pv)) = pv;
        else
            LASTVECTOR(theGrid) = (void*)pv;
        SUCCVC(after) = pv;
    }
}

	VSTART(pv) = NULL;

    /* counters */
    theGrid->nVector++;
}



void VectorXferCopy (DDD_OBJ obj, int proc, int prio)
{
	int 	nmat=0;
	MATRIX	*mat;
	VECTOR  *pv = (VECTOR *)obj;
	size_t  sizeArray[30]; /* TODO: define this static global TODO: take size as
 maximum of possible connections */

    PRINTDEBUG(dddif,1,("%2d: VectorXferCopy(): v=%08x/%x proc=%d prio=%d\n",
		me,DDD_InfoGlobalId(PARHDR(pv)),pv,proc,prio))

	for(mat=VSTART(pv); mat!=NULL; mat=MNEXT(mat))
	{
		sizeArray[nmat++] = MSIZE(mat);
	}


	PRINTDEBUG(dddif,2,("%2d:  VectorXferCopy(): v=%08x/%x AddData nmat=%d\n",\
		me,DDD_InfoGlobalId(PARHDR(pv)),pv,nmat))

	DDD_XferAddDataX(nmat,TypeMatrix,sizeArray);
}



void VectorGatherMatX (DDD_OBJ obj, int cnt, DDD_TYPE type_id, void **Data)
{
	VECTOR *vec = (VECTOR *)obj;
	MATRIX *mat;
	int nmat=0;

	PRINTDEBUG(dddif,3,("%2d:  VectorGatherMatX(): v=%08x/%x ID=%d cnt=%d type=%d veobj=%d\n",\
		me,DDD_InfoGlobalId(PARHDR(vec)),vec,ID(VOBJECT(vec)),cnt,type_id,OBJT(vec)))

	if (cnt<=0) return;

	for (mat=VSTART((VECTOR *) vec); mat!=NULL; mat=MNEXT(mat))
	{
		int Size;

		IFDEBUG(dddif,0)
		if (cnt<nmat+1)
		{
			PRINTDEBUG(dddif,0,("%2d:  VectorGatherMatX(): v=%x cnt=%d nmat=%d type=%d veobj=%d\n",me,vec,cnt,nmat,type_id,OBJT(vec)))
			assert(0);
		}
		ENDDEBUG

		Size = MSIZE(mat);
		memcpy(Data[nmat],mat,Size);

		PRINTDEBUG(dddif,3,("%2d:  VectorGatherMatX(): v=%x mat=%x Size=%d nodetoID=%d\n",me,vec,mat,Size,ID(VOBJECT(MDEST(mat)))))

		nmat++;
	}
}


void VectorScatterConnX (DDD_OBJ obj, int cnt, DDD_TYPE type_id, void **Data)
{
	VECTOR *vec = (VECTOR *)obj;
	CONNECTION *first=NULL, *last=NULL;
	GRID *theGrid = NULL;
	int i, new_conns = 0;
	int  level = DDD_InfoAttr(PARHDR(vec));

	theGrid = GRID_ON_LEVEL(dddctrl.currMG,level);

	PRINTDEBUG(dddif,3,("%2d:  VectorScatterConnX(): v=%08x/%x cnt=%d type=%d veobj=%d\n",\
		me,DDD_InfoGlobalId(PARHDR(vec)),vec,cnt,type_id,OBJT(vec)))

	if (cnt<=0) return;

	for (i=0; i<cnt; i++)
	{
		MATRIX *mcopy = (MATRIX *)Data[i];

		if (MDEST(mcopy)==NULL)
		{
			/* destination vector is not on this processor */
			/* -> matrix entry is useless, throw away */
			PRINTDEBUG(dddif,4,("%2d:  VectorScatterConnX(): v=%x mat=%x Size=%d, useless\n",me,vec,mcopy,MSIZE(mcopy)))
		}
		else
		{
			MATRIX *m;
			int found=FALSE;

			/* does matrix entry already exist? */
			/* TODO not nice, linear search, change this! */
			for (m=VSTART((VECTOR *)vec); m!=NULL && (!found); m=MNEXT(m))
			{
				if (MDEST(m)==MDEST(mcopy)) found=TRUE;
			}

			if (!found)
			{
				/* matrix entry is really new */

				if (MDIAG(mcopy))
				{
					/* matrix diagonal entry, no other vector is involved */
					CONNECTION *conn = (CONNECTION *)
						GetMem(dddctrl.currMG->theHeap,
							MSIZE(mcopy), FROM_BOTTOM);
					new_conns++;

					if (conn==NULL)
					{
						UserWriteF("%2d:  VectorScatterConnX(): can't get mem for conn=%x\n",conn);
						return;
					}
	
					PRINTDEBUG(dddif,4,("%2d:  VectorScatterConnX(): v=%x conn=%x Size=%d, diag\n",me,vec,conn,MSIZE(mcopy)))

					memcpy(conn,mcopy,MSIZE(mcopy));

					if (first==NULL) first = conn;
					else MNEXT(last) = conn;
					last = conn;
				}
				else
				{
					/* matrix off-diagonal entry, another vector is involved */
					VECTOR *other = MDEST(mcopy);
					MATRIX *m, *back=NULL, *newm;
	
					/* does connection already exist for other vec? */
					/* TODO not nice, linear search, change this! */

					for (m=VSTART((VECTOR *)other); m!=NULL&&back==NULL; m=MNEXT(m))
					{
						if (MDEST(m)==vec) back=m;
					}

					if (back==NULL)
					{
						/* no backward entry, create connection */
						MATRIX *otherm;
						CONNECTION *conn = (CONNECTION *)
							GetMem(dddctrl.currMG->theHeap,
								2 * MSIZE(mcopy), FROM_BOTTOM);
						new_conns++;

						if (conn==NULL)
						{
							UserWriteF("%2d:  VectorScatterConnX(): can't get mem for mat=%x\n",mcopy);
							return;
						}
	

						if (MOFFSET(mcopy))
						{
							newm = (MATRIX *) ((char *)conn+MSIZE(mcopy));
							otherm = (MATRIX *) conn;

						PRINTDEBUG(dddif,4,("%2d:  VectorScatterConnX(): v=%x conn=%x newm=%x Size=%d vectoID=%d, getmem\n",me,vec,conn,newm, MSIZE(mcopy),ID(MDEST(mcopy))))
						}
						else
						{
							newm = (MATRIX *) conn;
							otherm = (MATRIX *) ((char *)conn+MSIZE(mcopy));

						PRINTDEBUG(dddif,4,("%2d:  VectorScatterConnX(): v=%x conn=%x newm=%x Size=%d vectoID=%d, getmem\n",me,vec,conn,newm, MSIZE(mcopy),ID(MDEST(mcopy))))
						}

						MDEST(otherm) = NULL;
					}
					else
					{
						/* backward entry found, use existing connection */
						newm = MADJ(back);

						PRINTDEBUG(dddif,4,("%2d:  VectorScatterConnX(): v=%x back=%x newm=%x Size=%d vectoID=%d, reuse\n",me,vec,back,newm,MSIZE(mcopy),ID(MDEST(mcopy))))
					}

					memcpy(newm, mcopy, MSIZE(mcopy));

					if (first==NULL) first = newm;
					else MNEXT(last) = newm;
					last = newm;
				}
			}
		}
	}

	/* enter matrix list at the beginning of existing list for this vector */
	/* ensure diagonal entry being at first position */
	if (new_conns > 0)
		if (VSTART(vec)!=NULL)
		{
			MNEXT(last) = MNEXT(VSTART(vec));
			MNEXT(VSTART(vec)) = first;
		}
		else
		{
			MNEXT(last) = VSTART(vec);
			VSTART(vec) = first;
		}

	/* count new connections */
	NC(theGrid) += new_conns;
}



void VectorObjMkCons(DDD_OBJ obj)
{
	VECTOR *vec = (VECTOR *) obj;
	MATRIX *m;

	PRINTDEBUG(dddif,2,("%2d: VectorObjMkCons(): v=%08x/%x VEOBJ=%d\n",me,DDD_InfoGlobalId(PARHDR(vec)),vec,OBJT(vec)))
	
/*
	NOTE (TODO): this might be too less. for n2n transfer, connections
	might be set up consisting of two matrix structures transfered from
	different procs. this code will NOT handle that case, the connection
	will be created with the first matrix and destructed here. when the
	second message arrives, the second matrix will lead to construction
	of a second connection, which will also be deleted here. we would
	need a mkcons after all messages to handle that case. (NIY in ddd 1.6.5)
*/

	/* find and kill useless connections */
	for (m=VSTART((VECTOR *)vec); m!=NULL; m=MNEXT(m))
	{
		if (MDEST(MADJ(m))==NULL)
		{
			PRINTDEBUG(dddif,4,("%2d:  VectorObjMkCons(): v=%x mat=%x vectoID=%d, find&kill, TODO\n",me,vec,m,ID(MDEST(m))))
			/* TODO find & kill is not done!!! */
		}
	}
}



/****************************************************************************/
/*																			*/
/* Function:  VectorDelete													*/
/*																			*/
/* Purpose:   remove vector from UG data structure.							*/
/*			  current implementation only for level 0 grids					*/
/*																			*/
/* Input:	  DDD_OBJ	obj:	the vector to handle						*/
/*																			*/
/* Output:	  void															*/
/*																			*/
/****************************************************************************/

void VectorDelete (DDD_OBJ obj)
{
	VECTOR		*pv = (VECTOR *)obj;
	GRID		*theGrid = NULL;
	int         level = DDD_InfoAttr(PARHDR(pv));

	PRINTDEBUG(dddif,2,("%2d: VectorDelete(): v=%08x/%x VOBJ=%d l=%d\n",me,\
		DDD_InfoGlobalId(PARHDR(pv)),pv,OBJT(pv),level))

	theGrid = GRID_ON_LEVEL(dddctrl.currMG,level);

	/* remove vector from its object */
	if (VTYPE(pv)==NODEVECTOR)
	{
		NVECTOR((NODE *)VOBJECT(pv))  = NULL;
	}
	else
	{
		PRINTDEBUG(dddif,0,("%2d: VectorDelete(): VTYPE=%d not implemented yet\n", me, VTYPE(pv)))
		assert(0);
	}

	/* dispose vector itself */
	DisposeVector(theGrid, pv);
}

void VectorPriorityUpdate (DDD_OBJ obj, int new)
{
	VECTOR *pv = (VECTOR *)obj;
	INT     level = DDD_InfoAttr(PARHDR(pv));
	GRID    *theGrid = GetGridOnDemand(dddctrl.currMG,level);
	INT		old = DDD_InfoPriority(PARHDR(pv));

	PRINTDEBUG(dddif,2,("%2d: VectorPriorityUpdate(): v=%08x/%x old=%d new=%d level=%d\n",me,\
		DDD_InfoGlobalId(PARHDR(pv)),pv,old,new,level))

	if (pv == NULL) return;
	if (old == new) return;

	if (old == PrioNone) {
		/* only valid for masters */
		ASSERT(new == PrioMaster);
		return;
	}
	if (new == PrioNone) {
		/* only valid when prio undefined */  
		printf("prio=%d\n",old);
		fflush(stdout);
		ASSERT(old <= 0);
		return;
	}

	GRID_UNLINK_VECTOR(theGrid,pv)

	GRID_LINK_VECTOR(theGrid,pv,new)
	
	return;
}

/****************************************************************************/
/****************************************************************************/
/*																			*/
/*		handlers for typeelement											*/
/*																			*/

/****************************************************************************/
/****************************************************************************/
/*																			*/
/*		handlers for typeivertex											*/
/*																			*/
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
/****************************************************************************/
/*																			*/
/*		handlers for typebvertex											*/
/*																			*/
/****************************************************************************/
/****************************************************************************/


void BVertexLDataConstructor (DDD_OBJ obj)
{
	V_BNDP((VERTEX *)obj) = NULL;
}


void VertexUpdate (DDD_OBJ obj)
{
	VERTEX  *pv = (VERTEX *) obj;
	VERTEX  *after;
	ELEMENT *theElement;
	GRID  *theGrid;
	int  level = DDD_InfoAttr(PARHDRV(pv));

	PRINTDEBUG(dddif,1,("%2d: VertexUpdate(): v=%x I/BVOBJ=%d\n",me,pv,OBJT(pv)))

	theGrid = GRID_ON_LEVEL(dddctrl.currMG,level);
	after = LASTVERTEX(theGrid);

        /* insert in vertex list */
        if (after==NULL)
        {
                SUCCV(pv) = FIRSTVERTEX(theGrid);
                PREDV(pv) = NULL;
                if (SUCCV(pv)!=NULL) PREDV(SUCCV(pv)) = pv;
                else LASTVERTEX(theGrid) = pv;
                FIRSTVERTEX(theGrid) = pv;
        }
        else
        {
                SUCCV(pv) = SUCCV(after);
                PREDV(pv) = after;
                if (SUCCV(pv)!=NULL) PREDV(SUCCV(pv)) = pv;
                else LASTVERTEX(theGrid) = pv;
                SUCCV(after) = pv;
        }

        /* counters */
        theGrid->nVert++;

		/* update ID of vertex */
		/* TODO: change to global id */
		ID(pv) = (theGrid->mg->vertIdCounter)++;
}


void BVertexXferCopy (DDD_OBJ obj, int proc, int prio)
{
	BVertexXferBndP (V_BNDP((VERTEX *)obj),proc,prio);
}


void BVertexGather (DDD_OBJ obj, int cnt, DDD_TYPE type_id, void *Data)
{
    BVertexGatherBndP (V_BNDP((VERTEX *)obj),cnt,Data);
}


void BVertexScatter (DDD_OBJ obj, int cnt, DDD_TYPE type_id, void *Data)
{
    BVertexScatterBndP (&(V_BNDP((VERTEX *)obj)),cnt,Data);
}


/****************************************************************************/
/****************************************************************************/
/*																			*/
/*		handlers for typenode												*/
/*																			*/
/****************************************************************************/
/****************************************************************************/


void NodeDestructor(DDD_OBJ obj)
{
	NODE *node	= (NODE *) obj;

	PRINTDEBUG(dddif,2,("%2d: NodeDestructor(): n=%x NDOBJ=%d\n",me,node,OBJT(node)))
}

void NodeObjInit(DDD_OBJ obj)
{
	NODE *node	= (NODE *) obj;

	PRINTDEBUG(dddif,2,("%2d: NodeObjInit(): n=%x NDOBJ=%d\n",me,node,OBJT(node)))
}


void NodeObjMkCons(DDD_OBJ obj)
{
	NODE *node	= (NODE *) obj;
/*
	LINK *link;
	NODE *nodeto;
*/

	PRINTDEBUG(dddif,2,("%2d: NodeObjMkCons(): n=%x NDOBJ=%d\n",me,node,OBJT(node)))

/*
	for (link=START(node); link!=NULL; link=NEXT(link))
	{
*/
		/* TODO: wird das hier noch benoetigt? */
		/*
			nodeto = NBNODE(link);
		*/

			/* restore pointer from vector to its edge */
		/*
			if (dddctrl.edgeData) 
				VOBJECT(EDVECTOR(MYEDGE(link))) = (void*)MYEDGE(link);	
		*/


		/* is nodeto really stored on this proc? */
		/* TODO: wird das hier noch benoetigt? */
		/*
		if (nodeto!=NULL)
		{
			NEXT(REVERSE(link)) = START(nodeto);
			START(nodeto) = REVERSE(link);
		}
		*/
/*
	}
*/

	/* reconstruct node pointer */
	if (dddctrl.nodeData && NVECTOR(node)) VOBJECT(NVECTOR(node)) = (void*)node;

}


/****************************************************************************/
/*																			*/
/* Function:  NodeUpdate     										        */
/*																			*/
/* Purpose:   update information related to a node.    						*/
/*			  current implementation only for level 0 grids					*/
/*																			*/
/* Input:	  DDD_OBJ	obj:	the node to handle							*/
/*																			*/
/* Output:	  void															*/
/*																			*/
/****************************************************************************/

void NodeUpdate (DDD_OBJ obj)
{
	NODE  *node = (NODE *)obj;
	NODE  *after;
	GRID  *theGrid;
	int   level = DDD_InfoAttr(PARHDR(node));
	int   prio = DDD_InfoPriority(PARHDR(node));

	PRINTDEBUG(dddif,1,("%2d: NodeUpdate(): n=%x NDOBJ=%d\n",me,node,OBJT(node)))

	theGrid = GRID_ON_LEVEL(dddctrl.currMG,level);
	after = LASTNODE(theGrid);

	GRID_LINK_NODE(theGrid,node,prio)

/* TODO: delete this */
if (0) {
        /* insert in vertex list */
        if (after==NULL)
        {
                SUCCN(node) = FIRSTNODE(theGrid);
                PREDN(node) = NULL;
                if (SUCCN(node)!=NULL) PREDN(SUCCN(node)) = node;
                /*FIRSTNODE(theGrid) = node;*/
                if (LASTNODE(theGrid)==NULL) LASTNODE(theGrid) = node;
        }
        else
        {
                SUCCN(node) = SUCCN(after);
                PREDN(node) = after;
                if (SUCCN(node)!=NULL) PREDN(SUCCN(node)) = node; else LASTNODE(theGrid) = node;
                SUCCN(after) = node;
        }
}

		START(node) = NULL;

        /* incremant counter */
        theGrid->nNode++;

		/* TODO: change to global id */
		ID(node) = (theGrid->mg->nodeIdCounter)++;
}

/****************************************************************************/
/*																			*/
/* Function:  NodeXferCopy											        */
/*																			*/
/* Purpose:   initiate dependent copy of data related to a node.	 		*/
/*																			*/
/* Input:	  DDD_OBJ	obj:	the object which is transfered to proc		*/
/*			  int		proc:	destination processor for that object		*/
/*			  int		prio:	priority of object new local object			*/
/*																			*/
/* Output:	  void															*/
/*																			*/
/****************************************************************************/

void NodeXferCopy (DDD_OBJ obj, int proc, int prio)
{
	int		nlink = 0;
	int		Size,i=0;
	LINK	*link;
	NODE	*node = (NODE *)obj;
	VECTOR	*vec = NULL;

	PRINTDEBUG(dddif,1,("%2d: NodeXferCopy(): n=%08x/%x proc=%d prio=%d\n",me,DDD_InfoGlobalId(PARHDR(node)),node,proc,prio))


	/* copy vertex */
	PRINTDEBUG(dddif,2,("%2d: NodeXferCopy(): n=%x Xfer v=%x\n",me,node,MYVERTEX(node)))

	DDD_XferCopyObj(PARHDRV(MYVERTEX(node)), proc, prio);

	/* copy vector if defined */
	if (dddctrl.nodeData)
	  {
		vec = NVECTOR(node);
		Size = sizeof(VECTOR)-sizeof(DOUBLE)+dddctrl.currMG->theFormat->VectorSizes[VTYPE(vec)];

		PRINTDEBUG(dddif,2,("%2d: NodeXferCopy(): n=%x Xfer NODEVEC=%x size=%d\n",me,node,vec,Size))

		  DDD_XferCopyObjX(PARHDR(vec), proc, prio, Size);
	  }
}


void NodeGatherEdge (DDD_OBJ n, int cnt, DDD_TYPE type_id, void *Data)
{
	char	*data;
	LINK	*link;
	EDGE	*edge;
	NODE	*node = (NODE *) n;

	data = (char *)Data;

	PRINTDEBUG(dddif,3,("%2d:NodeGatherEdge(): n=%x cnt=%d type=%d ndobj=%d\n",me,node,cnt,type_id,OBJT(node)))

	/* copy edge(s) of node */
	for (link=START(node); link!=NULL; link=NEXT(link))
	{
		PRINTDEBUG(dddif,4,("%2d:NodeGatherEdge():  n=%x link=%x XFERLINK=%d\n",me,node,link,XFERLINK(link)))

		switch (XFERLINK(link))
		{
		case COPY:
				PRINTDEBUG(dddif,4,("%2d:NodeGatherEdge():   n=%x copy link=%x\n",me,node,link))
				memcpy(data,MYEDGE(link),sizeof(EDGE));
				data += CEIL(sizeof(EDGE));

				/* clear XFERFLAG */
/*
				SETXFERLINK(link,CLEAR);
*/
				break;
		case TOUCHED:
				/* clear XFERFLAG */
				SETXFERLINK(link,CLEAR);
				break;
		default:
				break;
		}
	}
}


void NodeScatterEdge (DDD_OBJ n, int cnt, DDD_TYPE type_id, void *Data)
{
	char	*data;
	int		i;
	EDGE	*edge;
	LINK	*link,*prev;
	NODE	*node = (NODE *) n;
	GRID	*grid;
	int  level = DDD_InfoAttr(PARHDR(node));

	grid = GRID_ON_LEVEL(dddctrl.currMG,level);
	data = (char *)Data;

	/* increment counter */
	grid->nEdge+=cnt;

	PRINTDEBUG(dddif,3,("%2d:NodeScatterEdge(): n=%x cnt=%d type=%d ndobj=%d\n",me,node,cnt,type_id,OBJT(node)))

	edge = (EDGE *)GetMem(dddctrl.currMG->theHeap,sizeof(EDGE),FROM_BOTTOM);
	PRINTDEBUG(dddif,4,("%2d:NodeScatterEdge(): n=%x edge=%x size=%d\n",me,node,edge,CEIL(sizeof(EDGE))))


	/* copy data out of message */
	memcpy(edge,data,sizeof(EDGE));

	data+=CEIL(sizeof(EDGE));

	/* look which link belongs to that node 				*/
	/* TODO: change this to faster macro sequence if stable */
	if (XFERLINK(LINK0(edge))==COPY)			link = LINK0(edge);
	else if (XFERLINK(LINK1(edge))==COPY)		link = LINK1(edge);
	else PRINTDEBUG(dddif,0,("%2d NodeScatterEdge(): 	NO copy flag in edge=%x\n",me,edge))

	for (i=0,START(node)=link; i<cnt-1; i++,NEXT(prev)=link)
	{
		prev = link;

		/* CAUTION: perhaps need to change into +CEIL(SIZEOF(EDGE)) */
		edge = (EDGE *)GetMem(dddctrl.currMG->theHeap,sizeof(EDGE),FROM_BOTTOM);
		PRINTDEBUG(dddif,4,("%2d:NodeScatterEdge(): n=%x edge=%x size=%d\n",me,node,edge,CEIL(sizeof(EDGE))))

		/* copy data out of message */
		memcpy(edge,data,sizeof(EDGE));
		data+=CEIL(sizeof(EDGE));
		link++;

		/* look which link belongs to that node 				*/
		/* TODO: change this to faster macro sequence if stable */
		if (XFERLINK(LINK0(edge))==COPY)		link = LINK0(edge);
		else if (XFERLINK(LINK1(edge))==COPY)	link = LINK1(edge);
		else PRINTDEBUG(dddif,0,("%2d NodeScatterEdge(): 	NO copy flag in edge=%x\n",me,edge))
	}
	
	MNEXT(link) = NULL;
}
	
void NodePriorityUpdate (DDD_OBJ obj, int new)
{
	NODE *pn = (NODE *)obj;
	INT     level = DDD_InfoAttr(PARHDR(pn));
	GRID    *theGrid = GetGridOnDemand(dddctrl.currMG,level);
	INT		old = DDD_InfoPriority(PARHDR(pn));

	PRINTDEBUG(dddif,2,("%2d: NodePriorityUpdate(): n=%08x/%x old=%d new=%d level=%d\n",me,\
		DDD_InfoGlobalId(PARHDR(pn)),pn,old,new,level))

	if (pn == NULL) return;
	if (old == new) return;

	if (old == PrioNone) {
		/* only valid for masters */
		ASSERT(new == PrioMaster);
		return;
	}
	if (new == PrioNone) {
		/* only valid when prio undefined */  
		printf("prio=%d\n",old);
		fflush(stdout);
		ASSERT(old <= 0);
		return;
	}

	GRID_UNLINK_NODE(theGrid,pn)

	GRID_LINK_NODE(theGrid,pn,new)
	
	return;
}

/****************************************************************************/
/****************************************************************************/
/*																			*/
/*		handlers for typeelement											*/
/*																			*/
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* Function:  ElementLDataConstructor										*/
/*																			*/
/* Purpose:   update information related to an element.						*/
/*			  current implementation only for level 0 grids					*/
/*																			*/
/* Input:	  DDD_OBJ	obj:	the element to handle						*/
/*																			*/
/* Output:	  void															*/
/*																			*/
/****************************************************************************/

void ElementLDataConstructor (DDD_OBJ obj)
{
	int			i,sides;
	ELEMENT		*pe = (ELEMENT *)obj;
	ELEMENT		*after = NULL;
	ELEMENT		*before = NULL;
	VERTEX		*pv;
	EDGE		*theEdge;
	int         level = DDD_InfoAttr(PARHDRE(pe));
	int         prio = DDD_InfoPriority(PARHDRE(pe));
	GRID        *theGrid = GetGridOnDemand(dddctrl.currMG,level);

	PRINTDEBUG(dddif,2,("%2d: ElementLDataConsX(): pe=%08x/%x eID=%d EOBJ=%d l=%d\n",\
		me,DDD_InfoGlobalId(PARHDRE(pe)),pe,ID(pe),OBJT(pe),level))

	after = LASTELEMENT(theGrid);
	SETLEVEL(pe,level);
	sides = 0;

	GRID_LINK_ELEMENT(theGrid,pe,prio)

	if (OBJT(pe)==BEOBJ)
	{
		for (i=0; i<SIDES_OF_ELEM(pe); i++) 
		  SET_BNDS(pe,i,NULL);
	}

	theGrid->nElem++;

	/* TODO: in global id umrechnen */
	ID(pe) = (theGrid->mg->elemIdCounter)++;
}


/****************************************************************************/
/*																			*/
/* Function:  ElementDelete													*/
/*																			*/
/* Purpose:   remove element from UG data structure.						*/
/*			  current implementation only for level 0 grids					*/
/*																			*/
/* Input:	  DDD_OBJ	obj:	the element to handle						*/
/*																			*/
/* Output:	  void															*/
/*																			*/
/****************************************************************************/

void ElementDelete (DDD_OBJ obj)
{
	ELEMENT		*pe = (ELEMENT *)obj;
	GRID		*theGrid = NULL;
	int         level = DDD_InfoAttr(PARHDRE(pe));

	PRINTDEBUG(dddif,2,("%2d: ElementDelete(): e=%08x/%x EOBJ=%d l=%d\n",me,\
		DDD_InfoGlobalId(PARHDRE(pe)),pe,OBJT(pe),level))

	theGrid = GRID_ON_LEVEL(dddctrl.currMG,level);
	DisposeElement(theGrid, pe, FALSE);
}


/****************************************************************************/
/*																			*/
/* Function:  ElementXferCopy	  											*/
/*																			*/
/* Purpose:   initiate dependent copy of data related to an element. 		*/
/*			  this handler is implemented for an arbitrary element type.	*/
/*																			*/
/* Input:	  DDD_OBJ	obj:	the object which is transfered to proc		*/
/*			  int		proc:	destination processor for that object		*/
/*			  int		prio:	priority of object new local object			*/
/*																			*/
/* Output:	  void															*/
/*																			*/
/****************************************************************************/

void ElementXferCopy (DDD_OBJ obj, int proc, int prio)
{
	int      i,nsides;
	int		 Size;
	ELEMENT  *pe	=	(ELEMENT *)obj;
	VECTOR	 *vec;
	NODE	 *node;
	BNDS     *bnds[MAX_SIDES_OF_ELEM];

	PRINTDEBUG(dddif,1,("%d: ElementXferCopy(): "\
		"pe=%08x/%x eID=%d proc=%d prio=%d EOBJT=%d\n", me,\
		DDD_InfoGlobalId(PARHDRE(pe)),pe,ID(pe), proc, prio, OBJT(pe)))

	/* add element sides */
	/* must be done before any XferCopyObj-call! herein */
	/* or directly after XferCopyObj-call */
    
	if (OBJT(pe)==BEOBJ)
	  {
		nsides = SIDES_OF_ELEM(pe);
		for (i=0; i<nsides; i++)
		  bnds[i] = ELEM_BNDS(pe,i);
		BElementXferBndS(bnds,nsides,proc,prio);
	  }

	/* add edges of element */
	/* must be done before any XferCopyObj-call! herein */
	/* or directly after XferCopyObj-call */
	DDD_XferAddData(EDGES_OF_ELEM(pe), TypeEdge);

	/* copy corner nodes */
	for(i=0; i<CORNERS_OF_ELEM(pe); i++)
	{
		node = CORNER(pe,i);

		PRINTDEBUG(dddif,2,("%2d: ElementXferCopy():  e=%x Xfer n=%x i=%d\n",\
				me, pe, node, i))

		DDD_XferCopyObj(PARHDR(node), proc, prio);
	}



	/* send edge vectors */
	if (dddctrl.edgeData) {
		for (i=0; i<EDGES_OF_ELEM(pe); i++)
		{
			int Size;
			EDGE *edge = GetEdge(CORNER(pe,CORNER_OF_EDGE(pe,i,0)),
							 	CORNER(pe,CORNER_OF_EDGE(pe,i,1)));
			VECTOR *vec = EDVECTOR(edge);

			Size = sizeof(VECTOR)-sizeof(DOUBLE)+dddctrl.currMG->theFormat->VectorSizes[VTYPE(vec)];
			PRINTDEBUG(dddif,3,("%2d: ElementferCopy():  e=%x EDGEVEC=%x size=%d\n",me,pe,vec,Size))
			DDD_XferCopyObjX(PARHDR(vec), proc, prio, Size);
		}
	}



	/* copy element vector */
	if (dddctrl.elemData)
	  {
		vec = EVECTOR(pe);
		Size = sizeof(VECTOR)-sizeof(DOUBLE)+dddctrl.currMG->theFormat->VectorSizes[VTYPE(vec)];
		
		PRINTDEBUG(dddif,2,("%2d:ElementXferCopy(): e=%x ELEMVEC=%x size=%d\n",me,pe,vec,Size))

		  DDD_XferCopyObjX(PARHDR(vec), proc, prio, Size);
	  }

	/* copy sidevectors */
	if (dddctrl.sideData)
	  {
		for (i=0; i<SIDES_OF_ELEM(pe); i++)
		  {
			vec = SVECTOR(pe,i);
			Size = sizeof(VECTOR)-sizeof(DOUBLE)+dddctrl.currMG->theFormat->VectorSizes[VTYPE(vec)];
			if (XFERVECTOR(vec) == 0)
			  {
				PRINTDEBUG(dddif,2,("%2d:ElementXferCopy(): e=%x SIDEVEC=%x size=%d\n",me,pe,vec,Size))
				  
			    DDD_XferCopyObjX(PARHDR(vec), proc, prio, Size);
				SETXFERVECTOR(vec,1);
			  }
		  }
	  }
}


/****************************************************************************/

static void ElemGatherEdge (ELEMENT *pe, int cnt, char *data)
{
	int i;
	int size = sizeof(EDGE) - ((dddctrl.edgeData)? 0 : sizeof(VECTOR*));

	PRINTDEBUG(dddif,3,("%2d:  ElemGatherEdge(): pe=%08x/%x cnt=%d size=%d\n",me,DDD_InfoGlobalId(PARHDRE(pe)),pe,cnt,size))

	/* copy edges into message */
	for (i=0; i<EDGES_OF_ELEM(pe); i++)
	{
		EDGE *edge = GetEdge(CORNER(pe,CORNER_OF_EDGE(pe,i,0)),
							 CORNER(pe,CORNER_OF_EDGE(pe,i,1)));
		ASSERT(edge!=NULL);
		memcpy(data, (char *)edge, size);
		data += CEIL(size);

		
		PRINTDEBUG(dddif,4,("%2d:  ElemGatherEdge(): pe=%x i=%d n1=%08x n2=%08x nmid=%08x\n",me,pe,i,(edge->links[0].nbnode),edge->links[1].nbnode,edge->midnode))
/*
		PRINTDEBUG(dddif,4,("%2d:  ElemGatherEdge(): pe=%x i=%d n1=%08x n2=%08x nmid=%08x\n",me,pe,i,DDD_InfoGlobalId(PARHDR(edge->links[0].nbnode)),DDD_InfoGlobalId(PARHDR(edge->links[1].nbnode)),DDD_InfoGlobalId(PARHDR(edge->midnode))))
*/
	}
}


static void ElemScatterEdge (ELEMENT *pe, int cnt, char *data)
{
	int i;
	int size = sizeof(EDGE) - ((dddctrl.edgeData)? 0 : sizeof(VECTOR*));
	int    level = DDD_InfoAttr(PARHDRE(pe));
	GRID  *theGrid = GetGridOnDemand(dddctrl.currMG,level);


	PRINTDEBUG(dddif,3,("%2d:  ElemScatterEdge(): pe=%08x/%x cnt=%d\n",me,DDD_InfoGlobalId(PARHDRE(pe)),pe,cnt))

	/* retrieve edges from message */
	for (i=0; i<cnt; i++)
	{
		EDGE *enew, *ecopy = (EDGE *)data;
		data += CEIL(size);

		PRINTDEBUG(dddif,4,("%2d:  ElemScatterEdge(): pe=%x i=%d n1=%08x n2=%08x midnode=%08x\n",
				me,pe,i,NBNODE(LINK0(ecopy)),ecopy->links[1].nbnode,ecopy->midnode))

		enew = CreateEdge(theGrid, NBNODE(LINK0(ecopy)), NBNODE(LINK1(ecopy)), FALSE);
		PRINTDEBUG(dddif,1,("%d: ElemScatterEdge(): pe=%x create edge=%x e%d%d for n0=%x n1=%x\n",me,pe,enew,
				ID(NBNODE(LINK0(ecopy))),ID(NBNODE(LINK1(ecopy))),
				NBNODE(LINK0(ecopy)),NBNODE(LINK1(ecopy))));
		if (enew == NULL) {
			PRINTDEBUG(dddif,1,("%d:  ElemScatterEdge(): ERROR pe=%x i=%d CreateEdge returned NULL\n",me,pe,i));
			ASSERT(0);
		}
		{
		EDGE *edge0,*edge1;
			edge0 = GetEdge(NBNODE(LINK0(ecopy)),NBNODE(LINK1(ecopy)));
			edge1 = GetEdge(NBNODE(LINK1(ecopy)),NBNODE(LINK0(ecopy)));
			if (edge0 != edge1) 
				PRINTDEBUG(dddif,1,("%d: ElemScatterEdge(): n0=%x n1=%x edge0=%x BUT edge1=%x\n",me,
						NBNODE(LINK0(ecopy)),NBNODE(LINK1(ecopy)),edge0,edge1));
		}

		MIDNODE(enew) = MIDNODE(ecopy);
		if (dddctrl.edgeData)
		{
			VOBJECT(EDVECTOR(enew)) = (void *)enew;
		}

/*
		PRINTDEBUG(dddif,4,("%2d:  ElemScatterEdge(): pe=%x i=%d n1=%08x n2=%08x nmid=%08x\n",me,pe,i,DDD_InfoGlobalId(PARHDR(edge->links[0].nbnode)),DDD_InfoGlobalId(PARHDR(edge->links[1].nbnode)),DDD_InfoGlobalId(PARHDR(edge->midnode))))
*/
	}
}

/****************************************************************************/


void ElemGatherI (DDD_OBJ obj, int cnt, DDD_TYPE type_id, void *data)
{
	/* type_id is always TypeEdge */
	ElemGatherEdge((ELEMENT *)obj, cnt, (char *)data);
}


void ElemScatterI (DDD_OBJ obj, int cnt, DDD_TYPE type_id, void *data)
{
	/* type_id is always TypeEdge */
	ElemScatterEdge((ELEMENT *)obj, cnt, (char *)data);
}


void ElemGatherB (DDD_OBJ obj, int cnt, DDD_TYPE type_id, void *data)
{
	int      i,nsides;
	BNDS     *bnds[MAX_SIDES_OF_ELEM];
	ELEMENT  *pe = (ELEMENT *)obj;

	/* type_id is TypeEdge or other */
	if (type_id==TypeEdge)
	{
		ElemGatherEdge(pe, cnt, (char *)data);
	} else
	{
		nsides = SIDES_OF_ELEM(pe);
		for (i=0; i<nsides; i++)
		  bnds[i] = ELEM_BNDS(pe,i);
		BElementGatherBndS(bnds, nsides, cnt, (char *)data);
	}
}


void ElemScatterB (DDD_OBJ obj, int cnt, DDD_TYPE type_id, void *data)
{
	int      i,nsides;
	BNDS     *bnds[MAX_SIDES_OF_ELEM];
	ELEMENT  *pe = (ELEMENT *)obj;

	/* type_id is TypeEdge or other */
	if (type_id==TypeEdge)
	{
		ElemScatterEdge(pe, cnt, (char *)data);
	} else
	{
		nsides = SIDES_OF_ELEM(pe);
		for (i=0; i<nsides; i++)
		  bnds[i] = ELEM_BNDS(pe,i);
		BElementScatterBndS(bnds, nsides, cnt, (char *)data);
		for (i=0; i<nsides; i++)
		  SET_BNDS(pe,i,bnds[i]);
	}
}


/****************************************************************************/


/* two versions of ElementObjMkCons ... */


void ElementObjMkCons_Xfer (DDD_OBJ obj)
{
	int i;
	ELEMENT  *pe	=	(ELEMENT *)obj;

	/* reconstruct pointer from vectors */
	if (dddctrl.elemData) VOBJECT(EVECTOR(pe)) = (void*)pe;

	if (dddctrl.sideData)
	{
		for (i=0; i<SIDES_OF_ELEM(pe); i++) VOBJECT(SVECTOR(pe,i)) = (void*)pe;
	}
}



void ElementObjMkCons_Refine (DDD_OBJ obj)
{
	int i,j;
	ELEMENT  *pe	=	(ELEMENT *)obj;
	VERTEX *pv;

	PRINTDEBUG(dddif,0,("%2d: ElementObjMkCons_Refine(): pe=%x/%d\n",me,pe,ID(pe)))

	/* reconstruct pointer from vectors */
	if (dddctrl.elemData) VOBJECT(EVECTOR(pe)) = (void*)pe;

	if (dddctrl.sideData)
	{
		for (i=0; i<SIDES_OF_ELEM(pe); i++) VOBJECT(SVECTOR(pe,i)) = (void*)pe;
	}

	/* connect with father */
	{
		ELEMENT *father = EFATHER(pe);
		if (father != NULL) {
			assert(NSONS(father)<NSONS_OF_RULE(MARK2RULEADR(father,REFINE(father))));

			#ifdef __THREEDIM__
			/* insert only first son */
			if (SON(father,0) == NULL)
			#endif
			SET_SON(father,NSONS(father),pe);
			SETNSONS(father,NSONS(father)+1);
		}
	}

}


void ElementPriorityUpdate (DDD_OBJ obj, int new)
{
	ELEMENT *pe = (ELEMENT *)obj;
	INT     level = DDD_InfoAttr(PARHDRE(pe));
	GRID    *theGrid = GetGridOnDemand(dddctrl.currMG,level);
	INT		old = DDD_InfoPriority(PARHDRE(pe));

	PRINTDEBUG(dddif,2,("%2d: ElementPriorityUpdate(): e=%08x/%x old=%d new=%d level=%d\n",me,\
		DDD_InfoGlobalId(PARHDRE(pe)),pe,old,new,level))

	if (pe == NULL) return;
	if (old == new) return;

	if (old == PrioNone) {
		/* only valid for masters */
		ASSERT(new == PrioMaster);
		return;
	}
	if (new == PrioNone) {
		/* only valid when prio undefined */  
		printf("prio=%d\n",old);
		fflush(stdout);
		ASSERT(old <= 0);
		return;
	}

	GRID_UNLINK_ELEMENT(theGrid,pe)

	GRID_LINK_ELEMENT(theGrid,pe,new)
	
	return;
}

/****************************************************************************/
/****************************************************************************/
/*																			*/
/*		handlers for typeedge    	 										*/
/*																			*/
/****************************************************************************/
/****************************************************************************/

/* CAUTION: */
/* TODO: delete this, Update handlers are not called for DDD data objects */
void EdgeUpdate (DDD_OBJ obj)
{
	EDGE *pe = (EDGE *)obj;
	LINK *link0,*link1;
	GRID *theGrid = NULL;

	PRINTDEBUG(dddif,2,("%2d:EdgeUpdate(): edge=%x EDOBJT=%d\n",me,pe,OBJT(pe)))

	theGrid = GRID_ON_LEVEL(dddctrl.currMG,0);

	/* increment counter */
	theGrid->nEdge++;
}	




/****************************************************************************/


/* init handlers for all element */
static void ElemHandlerInit (DDD_TYPE etype, INT handlerSet)
{
	DDD_HandlerRegister(etype,
		HANDLER_LDATACONSTRUCTOR, ElementLDataConstructor,
		HANDLER_DELETE,           ElementDelete,
		HANDLER_XFERCOPY,		  ElementXferCopy,
		HANDLER_SETPRIORITY,	  ElementPriorityUpdate,
		HANDLER_END
	);


	switch (handlerSet)
	{
		case HSET_XFER:
			DDD_HandlerRegister(etype,
				HANDLER_OBJMKCONS,        ElementObjMkCons_Xfer,
				HANDLER_END
			);
			break;

		case HSET_REFINE:
			DDD_HandlerRegister(etype,
				HANDLER_OBJMKCONS,        ElementObjMkCons_Refine,
				HANDLER_END
			);
			break;
	}
}


/* init handlers for inner element */
static void IElemHandlerInit (DDD_TYPE etype, INT handlerSet)
{
	/* init standard elem handlers */
	ElemHandlerInit(etype, handlerSet);

	/* init additional handlers, necessary for inside management */
	DDD_HandlerRegister(etype,
		HANDLER_XFERGATHER,		  ElemGatherI,
		HANDLER_XFERSCATTER,	  ElemScatterI,
		HANDLER_END
	);
}


/* init handlers for boundary element */
static void BElemHandlerInit (DDD_TYPE etype, INT handlerSet)
{
	/* init standard elem handlers */
	ElemHandlerInit(etype, handlerSet);

	/* init additional handlers, necessary for boundary management */
	DDD_HandlerRegister(etype,
		HANDLER_XFERGATHER,		  ElemGatherB,
		HANDLER_XFERSCATTER,	  ElemScatterB,
		HANDLER_END
	);
}


/****************************************************************************/



/* init all handlers necessary for grid xfer */
void ddd_HandlerInit (INT handlerSet)
{
	DDD_HandlerRegister(TypeVector,
		HANDLER_UPDATE,		    VectorUpdate,
		HANDLER_DELETE,         VectorDelete,
		HANDLER_XFERCOPY,		VectorXferCopy,
		HANDLER_XFERGATHERX,	VectorGatherMatX,
		HANDLER_XFERSCATTERX,	VectorScatterConnX,
		HANDLER_OBJMKCONS,		VectorObjMkCons,
		HANDLER_SETPRIORITY,	VectorPriorityUpdate,
		HANDLER_END
	);	

	DDD_HandlerRegister(TypeIVertex,
		HANDLER_UPDATE,		    VertexUpdate,
		HANDLER_END
	);

	DDD_HandlerRegister(TypeBVertex,
		HANDLER_LDATACONSTRUCTOR,	BVertexLDataConstructor,
		HANDLER_UPDATE,		    VertexUpdate,
		HANDLER_XFERCOPY,		BVertexXferCopy,
		HANDLER_XFERGATHER,		BVertexGather,
		HANDLER_XFERSCATTER,	BVertexScatter,
		HANDLER_END
	);

	DDD_HandlerRegister(TypeNode,
		HANDLER_LDATACONSTRUCTOR,	NodeObjInit,
		HANDLER_DESTRUCTOR,			NodeDestructor,
		HANDLER_OBJMKCONS,			NodeObjMkCons,
		HANDLER_UPDATE,				NodeUpdate,
		HANDLER_XFERCOPY,			NodeXferCopy,
		HANDLER_SETPRIORITY,		NodePriorityUpdate,
/*
		HANDLER_XFERGATHER,			NodeGatherEdge,
		HANDLER_XFERSCATTER,		NodeScatterEdge,
*/
		HANDLER_END
	);


	#ifdef __TWODIM__
	IElemHandlerInit(TypeTrElem, handlerSet);
	BElemHandlerInit(TypeTrBElem, handlerSet);

	IElemHandlerInit(TypeQuElem, handlerSet);
	BElemHandlerInit(TypeQuBElem, handlerSet);
	#endif

	#ifdef __THREEDIM__
	IElemHandlerInit(TypeTeElem, handlerSet);
	BElemHandlerInit(TypeTeBElem, handlerSet);

	IElemHandlerInit(TypePyElem, handlerSet);
	BElemHandlerInit(TypePyBElem, handlerSet);

	IElemHandlerInit(TypeHeElem, handlerSet);
	BElemHandlerInit(TypeHeBElem, handlerSet);

	#endif


	DDD_HandlerRegister(TypeEdge,
		HANDLER_UPDATE,	EdgeUpdate,
		HANDLER_END
	);

    DomHandlerInit(handlerSet);
}

#endif /* ModelP */

