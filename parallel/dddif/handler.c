// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/****************************************************************************/
/*																			*/
/* File:	  handler.c														*/
/*																			*/
/* Purpose:   defines the handlers used by ddd during data management.      */
/*			  																*/
/*																			*/
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

#include "compiler.h"
#include "gm.h"
#include "parallel.h"
#include "heaps.h"
#include "debug.h"


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

MULTIGRID *DDD_currMG;

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static int nodedata;
static int edgedata;
static int sidedata;
static int elementdata;

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* Function:  InitCurrMG													*/
/*																			*/
/* Purpose:   initialize the current multigird which is handled by DDD		*/
/*																			*/
/* Input:	  MULTIGRID *MG:	the multigrid to handle						*/
/*																			*/
/* Output:	  void															*/
/*																			*/
/*																			*/
/****************************************************************************/

void InitCurrMG(MULTIGRID *MG)
{
	DDD_currMG = MG;

	nodedata = TYPE_DEF_IN_MG(theMG,NODEVECTOR);
	edgedata = TYPE_DEF_IN_MG(theMG,EDGEVECTOR);
	sidedata = TYPE_DEF_IN_MG(theMG,SIDEVECTOR);
	elementdata = TYPE_DEF_IN_MG(theMG,ELEMVECTOR);
}


/****************************************************************************/
/*																			*/
/*  	For data management during redistribution and communication 		*/
/*		DDD needs for each DDD (data) object several handlers.				*/
/*		These are:															*/
/*			HANDLER_OBJINIT,           handler: init object        		    */
/*			HANDLER_OBJUPDATE,         handler: update objects internals    */
/*			HANDLER_OBJMKCONS,         handler: make obj consistent         */
/*			HANDLER_OBJDELETE,         handler: delete object               */
/*			HANDLER_XFERCOPY,          handler: copy cmd during xfer        */
/*			HANDLER_XFERDELETE,        handler: delete cmd during xfer      */
/*			HANDLER_XFERGATHER,        handler: send additional data        */
/*			HANDLER_XFERSCATTER,       handler: recv additional data        */
/*			HANDLER_COPYMANIP,         handler: manipulate incoming copy    */
/*			HANDLER_LDATARESTORE,      handler: restore LDATA after copying */
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
/*																			*/
/*			DDD data objects:												*/
/*				TypeConnection,												*/
/*				TypeVSegment,												*/
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

void VectorUpdate (OBJECT obj)
{
	VECTOR *pv = (VECTOR *)obj;
	VECTOR *after = NULL;
	GRID *theGrid = NULL;

	PRINTDEBUG(dddif,1,("%2d: VectorUpdate(): v=%x VEOBJ=%d\n",me,pv,OBJT(pv)))

	theGrid = GRID_ON_LEVEL(DDD_currMG,0);
	after = LASTVECTOR(theGrid);

        /* insert in vector list */
        if (after==NULL)
        {
                SUCCVC(pv) = (VECTOR*)theGrid->firstVector;
                PREDVC(pv) = NULL;
                if (SUCCVC(pv)!=NULL)
                        PREDVC(SUCCVC(pv)) = pv;
                theGrid->firstVector = (void*)pv;
                if (theGrid->lastVector==NULL)
                        theGrid->lastVector = (void*)pv;
        }
        else
        {
                SUCCVC(pv) = SUCCVC(after);
                PREDVC(pv) = after;
                if (SUCCVC(pv)!=NULL)
                        PREDVC(SUCCVC(pv)) = pv;
                else
                        theGrid->lastVector = (void*)pv;
                SUCCVC(after) = pv;
        }

		VSTART(pv) = NULL;

        /* counters */
        theGrid->nVector++;
}

void VectorXferCopy (OBJECT obj, int proc, int prio)
{
	int 	nmat;
	MATRIX	*mat;
	VECTOR  *vec = (VECTOR *)obj;
	int   	sizeArray[30]; /* TODO: define this static global TODO: take size as maximum of possible connections */

	PRINTDEBUG(dddif,2,("%2d:  VectorXferCopy(): v=%x AddData nmat=%d\n",me,vec,nmat))

	nmat=0;
	for(mat=VSTART(vec); mat!=NULL; mat=MNEXT(mat))
	{
		if (XFERMATX(mat) == COPY) 
		{
			PRINTDEBUG(dddif,3,("%2d: VectorXferCopy(): v=%x COPYFLAG already set for mat=%x\n",me,vec,mat))
			continue;
		}

		if (XFERMATX(mat) == TOUCHED || MDIAG(mat))
		{
			/* set XFERMATX to copy */
			SETXFERMATX(mat,COPY);

			PRINTDEBUG(dddif,3,("%2d: VectorXferCopy():  v=%x mat=%x XFERMATX=%d\n",me,vec,mat,XFERMATX(mat)))

			sizeArray[nmat++] = (MDIAG(mat)?MSIZE(mat):2*MSIZE(mat));
		}
		else 
		{
			SETXFERMATX(MADJ(mat),TOUCHED);
		}
	}

	PRINTDEBUG(dddif,2,("%2d:  VectorXferCopy(): v=%x AddData nmat=%d\n",me,vec,nmat))

	DDD_XferAddDataX(nmat,TypeConnection,sizeArray);
}

void VectorGatherConnX (OBJECT obj, int cnt, DDD_TYPE type_id, void **Data)
{
	VECTOR *vec = (VECTOR *)obj;
	CONNECTION *conn;
	int nconn=0;
	int Size;

	PRINTDEBUG(dddif,3,("%2d:  VectorGatherConnX(): v=%x cnt=%d type=%d veobj=%d\n",me,vec,cnt,type_id,OBJT(vec)))
	if (cnt<=0) return;

	for (conn=VSTART((VECTOR *) vec); conn!=NULL; conn=MNEXT(conn))
	{
		if (XFERMATX(conn)==COPY)
		{
			IFDEBUG(dddif,0)
			if (cnt<nconn+1)
			{
				PRINTDEBUG(dddif,0,("%2d:  VectorGatherConnX(): v=%x cnt=%d nconn=%d type=%d veobj=%d\n",me,vec,cnt,nconn,type_id,OBJT(vec)))
				return;
			}
			ENDDEBUG

			PRINTDEBUG(dddif,3,("%2d:  VectorGatherConnX(): v=%x conn=%x Size=%d \n",me,vec,conn,Size))

			Size = (MDIAG(conn)?MSIZE(conn):2*MSIZE(conn));
			memcpy(Data[nconn],MMYCON(conn),Size);

			/* save pointer to destination vector */
			if (MDEST(CMATRIX0(MMYCON(conn)))==vec && !MDIAG(conn)) 
				MDEST(CMATRIX0((CONNECTION *)Data[nconn])) = MDEST(conn);

			nconn++;
		}				
	}
}

void VectorScatterConnX (OBJECT obj, int cnt, DDD_TYPE type_id, void **Data)
{
	VECTOR *vec = (VECTOR *)obj;
	CONNECTION *conn,*prev;
	GRID *theGrid = NULL;
	int nconn=0;
	INT RootType, DestType, MType, ds, Diag, Size, i;

	theGrid = GRID_ON_LEVEL(DDD_currMG,0);
	Diag = 1; /* TODO: Is diagonal element?? */
	ds;   /* TODO: How big is ds?        */

	PRINTDEBUG(dddif,3,("%2d:  VectorScatterConnX(): v=%x cnt=%d type=%d veobj=%d\n",me,vec,cnt,type_id,OBJT(vec)))
	if (cnt<=0) return;

	Diag  = MDIAG((MATRIX *)Data[nconn]);
	Size  = MSIZE((MATRIX *)Data[nconn]);
	if (!Diag) Size  *= 2;
	if (MSIZEMAX<Size)
	{
		PRINTDEBUG(dddif,0,("%2d:  VectorScatterConnX(): Size=%d but MSIZEMAX=%d\n",Size,MSIZEMAX))
		return;
	}

	conn = (CONNECTION *)GetMem(DDD_currMG->theHeap,Size,FROM_BOTTOM);
	if (conn==NULL)
	{
		UserWrite("%2d:  VectorScatterConnX(): can't get mem for conn=%x\n",conn);
		return;
	}

	PRINTDEBUG(dddif,4,("%2d:  VectorScatterConnX(): v=%x conn=%x Size=%d\n",me,vec,conn,Size))
	memcpy(conn,Data[nconn++],Size);

	/* look which matrix belongs to that vector 				*/
	/* TODO: change this to faster macro sequence if stable 	*/
	if (!MDIAG(conn))
	{
		if (XFERMATX(CMATRIX0(conn))==COPY)			conn = CMATRIX0(conn);
		else if (XFERMATX(CMATRIX1(conn))==COPY)	
		{
			conn = CMATRIX1(conn);
			/* restore pointer to destination vector */
			MDEST(conn) = MDEST(CMATRIX0(MMYCON(conn)));
			MDEST(CMATRIX0(MMYCON(conn))) = vec;
		}
		else UserWrite("%2d NodeScatterEdge(): 	NO copy flag in conn=%x\n",me,conn);
	}

	/* TODO: this loop is no nice programming */
	for (i=1,VSTART((VECTOR *)vec)=conn; i<cnt; i++,MNEXT(prev)=conn)
	{
		Diag  = MDIAG((MATRIX *)Data[nconn]);
		Size  = MSIZE((MATRIX *)Data[nconn]);
		if (!Diag) Size  *= 2;
		if (MSIZEMAX<Size)
		{
			UserWrite("%2d:  VectorScatterConnX(): Size=%d but MSIZEMAX=%d\n",Size,MSIZEMAX);
			return;
		}

		prev = conn;
		conn = (CONNECTION*)GetMem(DDD_currMG->theHeap,Size,FROM_BOTTOM);
		if (conn==NULL)
		{
			UserWrite("%2d:  VectorScatterConnX(): ERROR can't get mem for conn=%x\n",conn);
		}

		PRINTDEBUG(dddif,4,("%2d:  VectorScatterConnX(): v=%x conn=%x Size=%d\n",me,vec,conn,Size))
		memcpy(conn,Data[nconn++],Size);

		/* look which matrix belongs to that vector 				*/
		/* TODO: change this to faster macro sequence if stable 	*/
		if (!MDIAG(conn))
			{
			if (XFERMATX(CMATRIX0(conn))==COPY)			conn = CMATRIX0(conn);
			else if (XFERMATX(CMATRIX1(conn))==COPY)
			{
				conn = CMATRIX1(conn);
				/* restore pointer to destination vector */
				MDEST(conn) = MDEST(CMATRIX0(MMYCON(conn)));
				MDEST(CMATRIX0(MMYCON(conn))) = vec;
			}
			else UserWrite("%2d NodeScatterEdge(): 	NO copy flag in conn=%x\n",me,conn);
		}
	}
	
	MNEXT(conn) = NULL;

	NC(theGrid) += cnt;
}

void VectorObjMkCons(OBJECT obj)
{
	CONNECTION *conn;
	VECTOR *vector	= (VECTOR *) obj;
	VECTOR *vectorto;

	PRINTDEBUG(dddif,2,("%2d: VectorObjMkCons(): v=%x VEOBJ=%d\n",me,vector,OBJT(vector)))

	for (conn=VSTART(vector); conn!=NULL; conn=MNEXT(conn))
	{
		if (MDIAG(conn)) continue;

		/* TODO: Durch schnelle Version ersetzen */
		if (XFERMATX(conn)==COPY) vectorto = MDEST(conn);
		else
		{
			if (XFERMATX(MADJ(conn)) != COPY) 
			{
				UserWrite("%2d VectorObjMkCons():     NO copy flag in conn with matrix=%x matrix=%x\n",me,conn,MADJ(conn));
			}
			continue;
		}
		MNEXT(MADJ(conn)) = VSTART(vectorto);
		VSTART(vectorto) = MADJ(conn);
		MDEST(MADJ(conn)) = vector;
	}

}


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

void VertexUpdate (OBJECT obj)
{
	VERTEX  *pv = (VERTEX *) obj;
	VERTEX  *after;
	GRID  *theGrid;

	theGrid = GRID_ON_LEVEL(DDD_currMG,0);
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
		VERTEX_ID(pv) = (theGrid->mg->vertIdCounter)++;
}

void BVertexXferCopy (OBJECT obj, int proc, int prio)
{
	int nvseg=0;
	VERTEX   *ver = (VERTEX *)obj;
	VSEGMENT *vseg;

	for (vseg=VSEG((VERTEX *)ver); vseg!=NULL; vseg=NEXTSEG(vseg))
	{
		nvseg++;
	}

	PRINTDEBUG(dddif,2,("%2d:  BVertexXferCopy(): v=%x AddData nvseg=%d\n",me,ver,nvseg))

	if (nvseg>0)	DDD_XferAddData(nvseg,TypeVSegment); 
}

void BVertexGatherVSegment (OBJECT ver, int cnt, DDD_TYPE type_id, void *Data)
{
	char *data;
	VSEGMENT *vseg;

	data = (char *)Data;

	PRINTDEBUG(dddif,3,("%2d:  BVertexGatherVSegment(): v=%x nvseg=%d type=%d bvobj=%d\n",me,ver,cnt,type_id,OBJT(ver)))

	for (vseg=VSEG((VERTEX *)ver); vseg!=NULL; vseg=NEXTSEG(vseg))
	{
		PRINTDEBUG(dddif,4,("%2d:  BVertexGatherVSegment(): v=%x vseg=%x\n",me,ver,vseg))

		memcpy(data,vseg,sizeof(VSEGMENT));
		/* copy segment id too */
		memcpy(data+sizeof(VSEGMENT),&SEGID(BSEGDESC(vseg)),sizeof(INT));
		data += CEIL(sizeof(VSEGMENT)+sizeof(INT));
	}
}

void BVertexScatterVSegment (OBJECT ver, int cnt, DDD_TYPE type_id, void *Data)
{
	char *data;
	int i;
	INT segmentid;
	VSEGMENT *vseg,*prev;
	GRID *theGrid = NULL;

	theGrid = GRID_ON_LEVEL(DDD_currMG,0);
	data = (char *)Data;

	PRINTDEBUG(dddif,3,("%2d: BVertexScatterVSegment(): v=%x nvseg=%d\n",me,ver,cnt))

	PRINTDEBUG(dddif,4,("%2d: BVertexScatterVSegment(): v=%x vseg=%x size=%d\n",me,ver,vseg,CEIL(sizeof(VSEGMENT))))

	vseg = (VSEGMENT *)GetMem(DDD_currMG->theHeap,sizeof(VSEGMENT),FROM_BOTTOM);
	memcpy(vseg,data,sizeof(VSEGMENT));
	memcpy(&segmentid,data+sizeof(VSEGMENT),sizeof(INT));
	data += CEIL(sizeof(VSEGMENT)+sizeof(INT));
	BSEGDESC(vseg) = &DDD_currMG->segments[segmentid]; 

	for (i=1,VSEG((VERTEX *)ver)=vseg; i<cnt; i++,NEXTSEG(prev)=vseg)
	{
		PRINTDEBUG(dddif,4,("%2d: BVertexScatterVSegment(): v=%x vseg=%x size=%d\n",me,ver,vseg,CEIL(sizeof(VSEGMENT))))

		prev = vseg;
		vseg = (VSEGMENT *)GetMem(DDD_currMG->theHeap,sizeof(VSEGMENT),FROM_BOTTOM);
		memcpy(vseg,data,sizeof(VSEGMENT));
		memcpy(&segmentid,data+sizeof(VSEGMENT),sizeof(INT));
		data += CEIL(sizeof(VSEGMENT)+sizeof(INT));
		BSEGDESC(vseg) = &DDD_currMG->segments[segmentid]; 
	}
	
	MNEXT(vseg) = NULL;
}

/****************************************************************************/
/****************************************************************************/
/*																			*/
/*		handlers for typenode												*/
/*																			*/
/****************************************************************************/
/****************************************************************************/

void NodeCopyManip(OBJECT copy)
{
	NODE *node	= (NODE *) copy;

	PRINTDEBUG(dddif,2,("%2d: NodeCopyManip(): n=%x NDOBJ=%d\n",me,node,OBJT(node)))
}

void NodeObjDelete(OBJECT obj)
{
	NODE *node	= (NODE *) obj;

	PRINTDEBUG(dddif,2,("%2d: NodeObjDelete(): n=%x NDOBJ=%d\n",me,node,OBJT(node)))
}

void NodeObjInit(OBJECT obj)
{
	NODE *node	= (NODE *) obj;

	PRINTDEBUG(dddif,2,("%2d: NodeObjInit(): n=%x NDOBJ=%d\n",me,node,OBJT(node)))
}

void NodeObjMkCons(OBJECT obj)
{
	LINK *link;
	NODE *node	= (NODE *) obj;
	NODE *nodeto;

	PRINTDEBUG(dddif,2,("%2d: NodeObjMkCons(): n=%x NDOBJ=%d\n",me,node,OBJT(node)))

	for (link=START(node); link!=NULL; link=NEXT(link))
	{
		/* TODO: Durch schnelle Version ersetzen */
		if (XFERLINK(link)==COPY) nodeto = NBNODE(link);
		else
		{
			if (XFERLINK(REVERSE(link)) == COPY) 
			{
				continue;
			}
			else
			{
				PRINTDEBUG(dddif,0,("%2d NodeObjMkCons():     NO copy flag in edge with link=%x link\n",me,link,REVERSE(link)))
				continue;
			}
		}
		NEXT(REVERSE(link)) = START(nodeto);
		START(nodeto) = REVERSE(link);
	}

}


/****************************************************************************/
/*																			*/
/* Function:  NodeUpdate     										        */
/*																			*/
/* Purpose:   update information related to a node.    						*/
/*			  current implementation only for level 0 grids					*/
/*																			*/
/* Input:	  OBJECT	obj:	the node to handle							*/
/*																			*/
/* Output:	  void															*/
/*																			*/
/****************************************************************************/

void NodeUpdate (OBJECT obj)
{
	NODE  *node = (NODE *)obj;
	NODE  *after;
	GRID  *theGrid;

	PRINTDEBUG(dddif,1,("%2d: NodeUpdate(): n=%x NDOBJ=%d\n",me,node,OBJT(node)))

	theGrid = GRID_ON_LEVEL(DDD_currMG,0);
	after = LASTNODE(theGrid);

        /* insert in vertex list */
        if (after==NULL)
        {
                SUCCN(node) = theGrid->firstNode;
                PREDN(node) = NULL;
                if (SUCCN(node)!=NULL) PREDN(SUCCN(node)) = node;
                theGrid->firstNode = node;
                if (theGrid->lastNode==NULL) theGrid->lastNode = node;
        }
        else
        {
                SUCCN(node) = SUCCN(after);
                PREDN(node) = after;
                if (SUCCN(node)!=NULL) PREDN(SUCCN(node)) = node; else theGrid->lastNode = node;
                SUCCN(after) = node;
        }

		START(node) = NULL;

        /* incremant counter */
        theGrid->nNode++;

		/* TODO: change to global id */
		NODE_ID(node) = (theGrid->mg->nodeIdCounter)++;
}

/****************************************************************************/
/*																			*/
/* Function:  NodeXferCopy											        */
/*																			*/
/* Purpose:   initiate dependent copy of data related to a node.	 		*/
/*																			*/
/* Input:	  OBJECT	obj:	the object which is transfered to proc		*/
/*			  int		proc:	destination processor for that object		*/
/*			  int		prio:	priority of object new local object			*/
/*																			*/
/* Output:	  void															*/
/*																			*/
/****************************************************************************/

void NodeXferCopy (OBJECT obj, int proc, int prio)
{
	int		nlink = 0;
	int		Size;
	LINK	*link;
	NODE	*node = (NODE *)obj;
	VECTOR	*vec = NULL;

	PRINTDEBUG(dddif,1,("%2d: NodeXferCopy(): n=%x COPYFLAG already set for LINK=%x\n",me,node,link))
	/* add links of node */
	for (link=START(node); link!=NULL; link=NEXT(link))
	{

		if (XFERLINK(link) == COPY) 
		{
			PRINTDEBUG(dddif,3,("%2d: NodeXferCopy(): n=%x COPYFLAG already set for LINK=%x\n",me,node,link))
			continue;
		}


		/* check whether corresponding node of this edge is also transferred */
		if (XFERLINK(link) == TOUCHED)
		{
			/* set XFERLINK to copy */
			SETXFERLINK(link,COPY);

			PRINTDEBUG(dddif,3,("%2d: NodeXferCopy():  n=%x link=%x XFERLINK=%d\n",me,node,link,XFERLINK(link)))

			/* increment counter for links to copy for this node */
			nlink++;

			/* send vector of this edge */	 
			if (edgedata)
			  {
				vec = EDVECTOR(MYEDGE(link));

				PRINTDEBUG(dddif,3,("%2d: NodeXferCopy():  n=%x EDGEVEC=%x size=%d\n",me,node,vec,theMG->theFormat->VectorSizes[VTYPE(vec)]))

				  /* TODO: */
				  /* CAUTION: must be called after XferAddData because */
                  /* of reference                                      */
				  /* to primary element, which is here node			   */
/*
   DDD_XferCopyObjX(DDD_OBJ(vec), proc, prio, DDD_currMG->theFormat->VectorSizes[VTYPE(vec)]);
*/

			  }
		}
		else
		{
			/* set XFERFLAG of corresponding link to TOUCHED */
			SETXFERLINK(REVERSE(link),TOUCHED);
		}
	}

	/* TODO: only send link whose corresponding node is also sent */
	/* CAUTION: must be called before any XferCopy because of reference */
	/*			to primary element, which is here the node      		*/
   	if (nlink>0)
	{
		PRINTDEBUG(dddif,2,("%2d: NodeXferCopy():  n=%x AddData nlink=%d\n",me,node,nlink))

		if (nlink >0)	DDD_XferAddData(nlink,TypeEdge);    
	}

	/* copy vertex */
	PRINTDEBUG(dddif,2,("%2d: NodeXferCopy(): n=%x Xfer v=%x\n",me,node,MYVERTEX(node)))

	DDD_XferCopyObj(DDD_OBJ(MYVERTEX(node)), proc, prio);

	/* copy vector if defined */
	if (nodedata)
	  {
		vec = NVECTOR(node);
		Size = sizeof(VECTOR)-sizeof(DOUBLE)+DDD_currMG->theFormat->VectorSizes[VTYPE(vec)];

		PRINTDEBUG(dddif,2,("%2d: NodeXferCopy(): n=%x Xfer NODEVEC=%x size=%d\n",me,node,vec,Size))

		  DDD_XferCopyObjX(DDD_OBJ(vec), proc, prio, Size);
	  }
}

void NodeGatherEdge (OBJECT n, int cnt, DDD_TYPE type_id, void *Data)
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

void NodeScatterEdge (OBJECT n, int cnt, DDD_TYPE type_id, void *Data)
{
	char	*data;
	int		i;
	EDGE	*edge;
	LINK	*link,*prev;
	NODE	*node = (NODE *) n;
	GRID	*grid;

	grid = GRID_ON_LEVEL(DDD_currMG,0);
	data = (char *)Data;

	/* increment counter */
	grid->nEdge+=cnt;

	PRINTDEBUG(dddif,3,("%2d:NodeScatterEdge(): n=%x cnt=%d type=%d ndobj=%d\n",me,node,cnt,type_id,OBJT(node)))

	edge = (EDGE *)GetMem(DDD_currMG->theHeap,sizeof(EDGE),FROM_BOTTOM);
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
		edge = (EDGE *)GetMem(DDD_currMG->theHeap,sizeof(EDGE),FROM_BOTTOM);
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
	
/****************************************************************************/
/****************************************************************************/
/*																			*/
/*		handlers for typebelement											*/
/*																			*/
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* Function:  ElementUpdate  										*/
/*																			*/
/* Purpose:   update information related to an element.						*/
/*			  current implementation only for level 0 grids					*/
/*																			*/
/* Input:	  OBJECT	obj:	the element to handle						*/
/*																			*/
/* Output:	  void															*/
/*																			*/
/****************************************************************************/

void ElementUpdate (OBJECT obj)
{
	int			i,sides;
	ELEMENT		*pe = (ELEMENT *)obj;
	ELEMENT		*after = NULL;
	GRID		*theGrid = NULL;

	PRINTDEBUG(dddif,2,("%2d: ElementUpdate(): e=%x EOBJ=%d\n",me,pe,OBJT(pe)))

	theGrid = GRID_ON_LEVEL(DDD_currMG,0);
	after = LASTELEMENT(theGrid);
	sides = 0;

	/* insert in element list */
	if (after==NULL)
	{
		SUCCE(pe) = theGrid->elements;
		PREDE(pe) = NULL;
		if (SUCCE(pe)!=NULL) 
			PREDE(SUCCE(pe)) = pe;
		else 
			theGrid->lastelement = pe;
		theGrid->elements = pe;
	}
	else
	{
		SUCCE(pe) = SUCCE(after);
		PREDE(pe) = after;
		if (SUCCE(pe)!=NULL)
			PREDE(SUCCE(pe)) = pe;
		else
			theGrid->lastelement = pe;
		SUCCE(after) = pe;
	}

	if (OBJT(pe)==BEOBJ)
	{
		for (i=0; i<SIDES_OF_ELEM(pe); i++) 
			if (SIDE(pe,i)!=NULL) sides++;
		theGrid->nSide += sides;
	}

	theGrid->nElem++;

	/* TODO: in global id umrechnen */
	ELEMENT_ID(pe) = (theGrid->mg->elemIdCounter)++;
}

/****************************************************************************/
/*																			*/
/* Function:  ElementXferCopy										*/
/*																			*/
/* Purpose:   initiate dependent copy of data related to an element. 		*/
/*			  this handler is implemented for an arbitrary element type.	*/
/*																			*/
/* Input:	  OBJECT	obj:	the object which is transfered to proc		*/
/*			  int		proc:	destination processor for that object		*/
/*			  int		prio:	priority of object new local object			*/
/*																			*/
/* Output:	  void															*/
/*																			*/
/****************************************************************************/

void ElementXferCopy (OBJECT obj, int proc, int prio)
{
	int      i,nelemside;
	ELEMENT  *pe	=	(ELEMENT *)obj;
	VECTOR	 *vec;
	NODE	 *node;

	PRINTDEBUG(dddif,1,("%d: ElementXferCopy(): pe=%x proc=%d prio=%d EOBJT=%d\n", me, obj, proc, prio, OBJT(pe)))

	/* add element sides; must be done before any XferCopyObj-call! */
	nelemside = 0;
	if (OBJT(pe)==BEOBJ)
	{
		for (i=0; i<SIDES_OF_ELEM(pe); i++)
			if (SIDE(pe,i) != NULL) nelemside++;

		PRINTDEBUG(dddif,2,("%2d: ElementXferCopy():  e=%x AddData nelemside=%d\n",me,pe,nelemside))

		if (nelemside>0) DDD_XferAddData(nelemside,TypeElementSide);
	}

	/* copy corner nodes */
	for(i=0; i<CORNERS_OF_ELEM(pe); i++)
	{
		node = CORNER(pe,i);
		if (XFERNODE(node) == 0)
		{
			PRINTDEBUG(dddif,2,("%2d:ElementXferCopy():  e=%x Xfer n=%x i=%d\n",me,pe,node,i))
			DDD_XferCopyObj(DDD_OBJ(node), proc, prio);
			SETXFERNODE(node,1);
		}
	}

	/* copy element vector */
	if (elemdata)
	  {
		vec = EVECTOR(pe);
		
		PRINTDEBUG(dddif,2,("%2d:ElementXferCopy(): e=%x ELEMVEC=%x size=%d\n",me,pe,vec,DDD_currMG->theFormat->VectorSizes[VTYPE(vec)]))

		  DDD_XferCopyObjX(DDD_OBJ(vec), proc, prio, DDD_currMG->theFormat->VectorSizes[VTYPE(vec)]);
	  }

	/* copy sidevectors */
	if (sidedata)
	  {
		for (i=0; i<SIDES_OF_ELEM(pe); i++)
		  {
			vec = SVECTOR(pe,i);
			if (XFERVECTOR(vec) == 0)
			  {
				PRINTDEBUG(dddif,2,("%2d:ElementXferCopy(): e=%x SIDEVEC=%x size=%d\n",me,pe,vec,DDD_currMG->theFormat->VectorSizes[VTYPE(vec)]))
				  
				  DDD_XferCopyObjX(DDD_OBJ(vec), proc, prio, DDD_currMG->theFormat->VectorSizes[VTYPE(vec)]);
				SETXFERVECTOR(vec,1);
			  }
		  }
	  }
}

void ElemGatherElemSide (OBJECT obj, int cnt, DDD_TYPE type_id, void *Data)
{
	int i;
	char *data;
	ELEMENTSIDE *elemside;
	ELEMENT  *pe	=	(ELEMENT *)obj;

	data = (char *)Data;

	PRINTDEBUG(dddif,3,("%2d:  ElemGatherElemSide(): e=%x nelemside=%d type=%d bvobj=%d\n",me,pe,cnt,type_id,OBJT(pe)))

	for (i=0; i<SIDES_OF_ELEM(pe); i++)
	{
		if (SIDE(pe,i) != NULL) 
		{
			PRINTDEBUG(dddif,4,("%2d:  ElemGatherElemSide(): e=%x elemside=%x side=%d segid=%d\n",me,pe,SIDE(pe,i),i,SEGID(SEGDESC(SIDE(pe,i)))))
			memcpy(data,SIDE(pe,i),sizeof(ELEMENTSIDE));
			/* copy segment id too */
			memcpy(data+sizeof(ELEMENTSIDE),&SEGID(SEGDESC(SIDE(pe,i))),sizeof(INT));
			data += CEIL(sizeof(ELEMENTSIDE)+sizeof(INT));
		}
	}
}

void ElemScatterElemSide (OBJECT obj, int cnt, DDD_TYPE type_id, void *Data)
{
	char *data;
	int i;
	INT segmentid;
	ELEMENTSIDE *elemside;
	ELEMENT  *pe	=	(ELEMENT *)obj;
	GRID *theGrid = NULL;

	theGrid = GRID_ON_LEVEL(DDD_currMG,0);
	data = (char *)Data;

	PRINTDEBUG(dddif,3,("%2d: ElemScatterElemSide(): pe=%x nelemside=%x\n",me,pe,cnt))

	for (i=0; i<SIDES_OF_ELEM(pe); i++)
	{
		if (SIDE(pe,i) != NULL) 
		{
			elemside = (ELEMENTSIDE *)GetMem(DDD_currMG->theHeap,sizeof(ELEMENTSIDE),FROM_BOTTOM);
			PRINTDEBUG(dddif,4,("%2d:  ElemScatterElemSide(): e=%x elemside=%x side=%d size=%d\n",me,pe,SIDE(pe,i),i,CEIL(sizeof(ELEMENTSIDE))))

			/* copy out elementside out of message and restore SEGDESC * */ 
			memcpy(elemside,data,sizeof(ELEMENTSIDE));
			memcpy(&segmentid,data+sizeof(ELEMENTSIDE),sizeof(INT));
			data += CEIL(sizeof(ELEMENTSIDE)+sizeof(INT));
			SEGDESC(elemside) = &DDD_currMG->segments[segmentid]; 
			SET_SIDE(pe,i,elemside);

			PRINTDEBUG(dddif,4,("%2d:  ElemScatterElemSide(): e=%x elemside=%x side=%d segid=%d\n",me,pe,SIDE(pe,i),i,SEGID(SEGDESC(SIDE(pe,i)))))

			/* put into double linked list */
			PREDS(elemside) = NULL;
			SUCCS(elemside) = FIRSTELEMSIDE(theGrid);
			if (PREDS(FIRSTELEMSIDE(theGrid))) PREDS(FIRSTELEMSIDE(theGrid)) = elemside;
			FIRSTELEMSIDE(theGrid) = elemside;
		}
	}
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
void EdgeUpdate (OBJECT obj)
{
	EDGE *pe = (EDGE *)obj;
	LINK *link0,*link1;
	GRID *theGrid = NULL;

	PRINTDEBUG(dddif,2,("%2d:EdgeUpdate(): edge=%x EDOBJT=%d\n",me,pe,OBJT(pe)))

	theGrid = GRID_ON_LEVEL(DDD_currMG,0);

	/* increment counter */
	theGrid->nEdge++;
}	

void ddd_HandlerInit (void)
{

	DDD_HandlerRegister(TypeVector,
		HANDLER_OBJUPDATE,		VectorUpdate,
		HANDLER_XFERCOPY,		VectorXferCopy,
		HANDLER_XFERGATHERX,	VectorGatherConnX,
		HANDLER_XFERSCATTERX,	VectorScatterConnX,
		HANDLER_OBJMKCONS,		VectorObjMkCons,
		HANDLER_END
	);	

	DDD_HandlerRegister(TypeIVertex,
		HANDLER_OBJUPDATE,		VertexUpdate,
		HANDLER_END
	);

	DDD_HandlerRegister(TypeBVertex,
		HANDLER_OBJUPDATE,		VertexUpdate,
		HANDLER_XFERCOPY,		BVertexXferCopy,
		HANDLER_XFERGATHER,		BVertexGatherVSegment,
		HANDLER_XFERSCATTER,	BVertexScatterVSegment,
		HANDLER_END
	);

	DDD_HandlerRegister(TypeNode,
		HANDLER_COPYMANIP,		NodeCopyManip,
		HANDLER_OBJINIT,		NodeObjInit,
		HANDLER_OBJDELETE,		NodeObjDelete,
		HANDLER_OBJMKCONS,		NodeObjMkCons,
		HANDLER_OBJUPDATE,		NodeUpdate,
		HANDLER_XFERCOPY,		NodeXferCopy,
		HANDLER_XFERGATHER,		NodeGatherEdge,
		HANDLER_XFERSCATTER,	NodeScatterEdge,
		HANDLER_END
	);

/*
	DDD_HandlerRegister(TypeIElement,
		HANDLER_OBJUPDATE,	IElementUpdate,
		HANDLER_XFERCOPY,	IElementXferCopy,
		HANDLER_END
	);
*/

	#ifdef __TWODIM__
	DDD_HandlerRegister(TypeTrElem,
		HANDLER_OBJUPDATE,		ElementUpdate,
		HANDLER_XFERCOPY,		ElementXferCopy,
		HANDLER_END
	);

	DDD_HandlerRegister(TypeTrBElem,
		HANDLER_OBJUPDATE,		ElementUpdate,
		HANDLER_XFERCOPY,		ElementXferCopy,
		HANDLER_XFERGATHER,		ElemGatherElemSide,
		HANDLER_XFERSCATTER,	ElemScatterElemSide,
		HANDLER_END
	);

	DDD_HandlerRegister(TypeQuElem,
		HANDLER_OBJUPDATE,		ElementUpdate,
		HANDLER_XFERCOPY,		ElementXferCopy,
		HANDLER_END
	);

	DDD_HandlerRegister(TypeQuBElem,
		HANDLER_OBJUPDATE,		ElementUpdate,
		HANDLER_XFERCOPY,		ElementXferCopy,
		HANDLER_XFERGATHER,		ElemGatherElemSide,
		HANDLER_XFERSCATTER,	ElemScatterElemSide,
		HANDLER_END
	);
	#endif

	#ifdef __THREEDIM__
	DDD_HandlerRegister(TypeTeElem,
		HANDLER_OBJUPDATE,		ElementUpdate,
		HANDLER_XFERCOPY,		ElementXferCopy,
		HANDLER_END
	);

	DDD_HandlerRegister(TypeTeBElem,
		HANDLER_OBJUPDATE,		ElementUpdate,
		HANDLER_XFERCOPY,		ElementXferCopy,
		HANDLER_XFERGATHER,		ElemGatherElemSide,
		HANDLER_XFERSCATTER,	ElemScatterElemSide,
		HANDLER_END
	);
	#endif

	DDD_HandlerRegister(TypeEdge,
		HANDLER_OBJUPDATE,	EdgeUpdate,
		HANDLER_END
	);
}

#endif /* ModelP */
